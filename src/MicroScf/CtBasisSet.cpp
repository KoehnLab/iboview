/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#include <cmath>
#include <stdexcept>
#include <map> // for libmol export
#include <set> // for printing
#include "format.h"
#include "Ir.h"
#include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "CtBasisLibrary.h"
#include "CtMatrix.h"

// #ifdef _DEBUG
   #include "CxIo.h"
// #endif // _DEBUG

namespace ct {


FBasisSet::FBasisSet(FAtomSet const &AtomSet, FBasisContext Context_)
   : Context(Context_)
{
   LoadFromAtomSet(AtomSet, Context);
}

FBasisSet::FBasisSet(FBasisShell const *pShells_, size_t nShells_, FBasisContext Context_, std::string const &Name_)
   : Context(Context_), Name(Name_)
{
//    Shells.assign(&pShells_[0], &pShells_[0] + nShells_);
   Shells.clear();
   Shells.reserve(nShells_);
   for (size_t iSh = 0; iSh < nShells_; ++ iSh)
      Shells.push_back(pShells_[iSh]);
   Finalize();
}


FBasisSetPtr LoadNamedBasis(FRawAtomList const &Atoms, FBasisContext Context, std::string const &BasisName)
{
   FBasisSet::FBasisShellArray
      Shells;
   for (size_t iAt = 0; iAt < Atoms.size(); ++ iAt){
      FRawAtom const
         &Atom = Atoms[iAt];
      g_BasisSetLibrary.LoadBasisFunctions(Shells,
         Atom.iElement, BasisName, Atom.vPos, iAt);
   }
   return FBasisSetPtr(new FBasisSet(&Shells[0], Shells.size(), Context, BasisName));
}


class FBasisLoadError : public std::runtime_error
{
public:
   explicit FBasisLoadError(std::string const &Desc, int iAt, FAtom const &At, FBasisContext const &Context)
      : std::runtime_error(fmt::format("while instanciating basis context {} for atom {} {} at ({:.4f}, {:.4f}, {:.4f}): {}", GetBasisContextName(Context), 1+iAt, At.ElementName(), At.vPos[0], At.vPos[1], At.vPos[2], Desc)),
      m_iAt(iAt), m_Context(Context)
   {}
protected:
   int
      m_iAt;
   FBasisContext
      m_Context;
};


void FBasisSet::LoadFromAtomSet(FAtomSet const &Atoms, FBasisContext Context)
{
   std::string
      BasisName,
      UnassignedOrbBasis;
//       UnassignedOrbBasis("undefined");
   bool
      AllEqual = true;

   BasisName.reserve(64);
   for (size_t iAt = 0; iAt < Atoms.size(); ++ iAt){
      FAtom const
         &Atom = Atoms[iAt];

      // get name of basis we are supposed to be loading.
      FBasisDescs::const_iterator
         itDesc = Atom.BasisDesc.find(Context);
      if (itDesc != Atom.BasisDesc.end()) {
         // basis for context explicitly set
         BasisName = itDesc->second;
         if (BasisName.empty() && Atom.iElement > 0)
            throw FBasisLoadError("A basis name was explicitly assigned, but it was empty. If you *really* want no basis functions on the atom, set the basis to 'none'.", iAt, Atom, Context);
      } else {
         // basis for current context not set---go find one.
         std::string const
            *pOrbBasisAsRef = 0;
         // check if a main orbital basis is assigned. May need that to find a
         // suitable fitting basis set (e.g., MP2FIT sets are specific to the
         // orbital basis they go with)
         itDesc = Atom.BasisDesc.find(BASIS_Orbital);
         if (itDesc != Atom.BasisDesc.end()) {
            // found one.
            pOrbBasisAsRef = &itDesc->second;
         } else {
            if (Context == BASIS_Orbital)
               throw FBasisLoadError("No orbital basis set assigned to atom.", iAt, Atom, Context);
            // have a look if this is a context we can probably assign a reasonable
            // default basis for even without knowing the orbital basis specifically.
            bool
               NoRefBasisIsOkay = !NeedOrbitalBasisForContextDefault(Context);
            if (!NoRefBasisIsOkay) {
               throw FBasisLoadError("attempted to instanciate default basis for target context, but no orbital basis was assigned to decide on what that should be (to avoid this error, set the target context basis explicitly or set an orbital basis).", iAt, Atom, Context);
            }
            pOrbBasisAsRef = &UnassignedOrbBasis;
         }
         assert(pOrbBasisAsRef != 0);
         FindDefaultBasis(BasisName, Context, *pOrbBasisAsRef, &Atom.BasisDesc, Atom.iElement, Atom.nEcpElec);
         if (BasisName.empty())
            throw FBasisLoadError(fmt::format("No basis set for context {} assigned to atom, and no default basis assignment was programmed for orbital basis {}.", GetBasisContextName(Context), itDesc->second), iAt, Atom, Context);
      }

      if (this->Name.empty()) this->Name = BasisName;
      if (this->Name != BasisName) AllEqual = false;

      // we're importing molpro data... Molpro truncates basis names
      // at 32 characters. Do so here, too.
      BasisName.resize(std::min(size_t(BasisName.size()), size_t(32)));

      // also, convert everything to lowercase. We do that when importing sets.
      for (size_t i = 0; i < BasisName.size(); ++ i)
         BasisName[i] = ::tolower(BasisName[i]);

      // ask basis library to add the functions.
      size_t nShellsBefore = this->Shells.size();
      g_BasisSetLibrary.LoadBasisFunctions(this->Shells,
         Atom.iElement, BasisName,
         Atom.vPos, iAt);
//       xout << fmt::format("!LoadBasisFunctions for {} '{}': returned {} basis function shells.", Atom.GetElementName(), BasisName, this->Shells.size() - nShellsBefore) << std::endl;
      if (this->Shells.size() == nShellsBefore)
         throw FBasisLoadError("MicroScf internal error: BasisSetLibrary raised no exception, but also returned no basis functions.", iAt, Atom, Context);

   }

//    if (!AllEqual || this->Name.empty())
//       this->Name = GetBasisContextName(Context);
   IR_SUPPRESS_UNUSED_WARNING(AllEqual);

   // TODO:
   //   - make local copy of FGaussBfn objects and re-link this->Shells
   //     to these objects? This would reduce the memory spread of the data
   //     and prevent TLB misses during integration.

   Finalize();
}


void FBasisSet::Print(std::ostream &xout, unsigned Flags, std::string const &Ind) const
{
   bool
      PrintCentersInShortMode = false;
   char const
      *pTabRules = 0,
      *pTabRulesShort = 0;
   if (!bool(Flags & PRINT_FullDetail)) {
      pTabRules = "---------------------------------------------------------------------------------------";
      pTabRulesShort = pTabRules;
   } else {
      pTabRules = "-----------------------------------------    -----------------   ----------------------";
      pTabRulesShort = "-----------------------------------------";
   }
   if (Flags & PRINT_MetaInfo) {
//       xout << Ind << fmt::format("Basis set '{}' (for context {}):\n", Name, GetBasisContextName(Context));
      xout << Ind << fmt::format("{} basis set '{}':\n", GetBasisContextName(Context), Name);
      xout << Ind << "\n";
      xout << Ind << pTabRules << "\n";
   }
   if (!bool(Flags & PRINT_FullDetail)) {
//       xout << Ind << " iFn0   nFn   L:Num.Co.                ";
      xout << Ind << " iFn0   nFn   Num.Co./Ang.Mom.         ";
      if (PrintCentersInShortMode)
         xout << "    Center/Position               ";
      xout << "   Atom     Instanciated Library Entries\n";
      xout << Ind << pTabRules << std::endl;
      xout.flush();
      size_t
         nFnOff = 0,
         iShBeg = 0;
      while (iShBeg < Shells.size()) {
         int iCenter = Shells[iShBeg].iCenter;
         size_t iShEnd = iShBeg + 1;
         while (iShEnd < Shells.size() && iCenter == Shells[iShEnd].iCenter)
            iShEnd += 1;
         size_t nFnCen = 0;
         typedef std::set<std::string>
            FStringSet;
         FStringSet
            LibNames;
         std::stringstream
            ssFnDesc;
         for (size_t iSh = iShBeg; iSh < iShEnd; ++ iSh) {
            FBasisShell const &Sh = Shells[iSh];
            nFnCen += Sh.nFn();
            LibNames.insert(Sh.LibraryName());
            if (iSh != iShBeg)
               ssFnDesc << ":";
            ssFnDesc << Sh.nCo() << "spdfghiklmn"[Sh.l()];
         }
         {
            FBasisShell const &Sh = Shells[iShBeg];
            std::streampos p0 = xout.tellp(), p1;
            xout << fmt::format("{}{:5} {:5}   ", Ind, nFnOff, nFnCen);
            std::string sCenDesc;
            if (PrintCentersInShortMode)
               sCenDesc = fmt::format("{:8.4f} {:8.4f} {:8.4f}    ", Sh.vCenter[0], Sh.vCenter[1], Sh.vCenter[2]);
            xout << fmt::format("{:<24}   {}{:>3} {:<2}   ",
               ssFnDesc.str(), sCenDesc, Sh.iCenter+1, ElementNameFromNumber(Sh.iLibraryElement()));
            for (FStringSet::const_iterator it = LibNames.begin(); it != LibNames.end(); ++ it) {
               xout << " " << *it;
            }
            xout << std::endl;
         }
         nFnOff += nFnCen;
         iShBeg = iShEnd;
      }
   } else {
      xout << Ind << " Offs   NC/AM        Center/Position         Exponents           Contractions\n";
      xout << Ind << pTabRules << std::endl;
      xout.flush();
      size_t
         nOff = 0;
      for (size_t iSh = 0; iSh < Shells.size(); ++ iSh) {
         std::streampos p0 = xout.tellp(), p1;
         xout << fmt::format("{}{:5}  ", Ind, nOff);
         p1 = xout.tellp();
         Shells[iSh].PrintAligned(xout, Ind, uint(p1-p0));
         xout << std::endl;
         nOff += Shells[iSh].nFn();
      }
   }
   if (Flags & PRINT_MetaInfo) {
      xout << Ind << pTabRulesShort << "\n";
      xout << fmt::format("{}--> nFn = {}  nPrim = {}  nSh = {}  MaxL = {}  max(nShFn) = {}\n", Ind, nFn(), nPrimitiveGtos(), Shells.size(), nMaxL(), nFnOfLargestShell());
      xout << Ind << pTabRulesShort << std::endl;
   }
}


std::string FBasisSet::GetDesc(unsigned Flags, std::string const &Ind) const
{
   std::stringstream str;
   this->Print(str, Flags, Ind);
   return str.str();
}


size_t FBasisSet::nPrimitiveGtos() const
{
   size_t r = 0;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh)
      r += Shells[iSh].nExp() * Shells[iSh].nSh();
   return r;
}

size_t FBasisSet::nFnOfLargestShell() const
{
   size_t r = 0;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh)
      r = std::max(r, size_t(Shells[iSh].nFn()));
   return r;
}

size_t FBasisSet::nCoOfLargestShell() const
{
   size_t r = 0;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh)
      r = std::max(r, size_t(Shells[iSh].nCo()));
   return r;
}

unsigned FBasisSet::nMaxL() const
{
   unsigned r = 0;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh)
      r = std::max(unsigned(Shells[iSh].l()), unsigned(r));
   return r;
}

size_t FBasisSet::nFn() const
{
   size_t r = 0;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh)
      r += Shells[iSh].nFn();
   return r;
}



void AccGradient2ix(double *pGrad, FMatrixView const &Rdm, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor)
{
   assert(Rdm.nRows == BasisRow.nFn());
   assert(Rdm.nCols == BasisCol.nFn());

   return AccGradient2ix(pGrad, Rdm.pData, Rdm.nRowSt, Rdm.nColSt, BasisRow,
      BasisCol, Krn2i, Mem_, Prefactor);
}

void AccHessian2ix( FMatrixView Hess, FMatrixView const &Rdm, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor)
{
   assert(Rdm.nRows == BasisRow.nFn());
   assert(Rdm.nCols == BasisCol.nFn());
   assert(Hess.nRows == 3*BasisRow.nCen());
   assert(Hess.nCols == 3*BasisCol.nCen());

   return AccHessian2ix(Hess.pData, Hess.nRowSt, Hess.nColSt, Rdm.pData,
      Rdm.nRowSt, Rdm.nColSt, BasisRow, BasisCol, Krn2i, Mem_, Prefactor);
}


void MakeIntMatrix( FMatrixView &Dest, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor, bool Add )
{
   assert(Dest.nRows == BasisRow.nFn());
   assert(Dest.nCols == BasisCol.nFn());

   return MakeIntMatrix(Dest.pData, Dest.nRowSt, Dest.nColSt, BasisRow,
      BasisCol, Krn2i, Mem_, Prefactor, Add);
}

void MakeIntMatrix( FMatrixView &Dest, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor, bool Add )
{
   return MakeIntMatrix(Dest, *RowBasis.pRawBasis, *ColBasis.pRawBasis, Krn2i, Mem, Prefactor, Add);
}



size_t FBasisSet::nCenters() const
{
   size_t
      nCenters = 0;
   for ( size_t iSh = 0; iSh < Shells.size(); ++ iSh )
      if ( Shells[iSh].iCenter > 0 )
         nCenters = std::max(nCenters, 1 + static_cast<size_t>(Shells[iSh].iCenter));
   return nCenters;
}


void FBasisSet::MakeAtomOffsets( size_t *&pAtomShellOffsets, size_t *&pAtomBfnOffsets, size_t nAtoms_, FMemoryStack &Mem ) const
{
   // pAtomShellOffsets: maps atom id to first basis function shell
   //    (Shells[i]) of the atom
   // pAtomBfnOffsets: maps atom id to first basis function (AO index)
   //    of the atom
   size_t nAtoms = nCenters();
   assert( nAtoms <= nAtoms_ );
   // ^- might differ if there are some atoms without basis functions, or basis funcitons without atoms (e.g., bond functions).
   nAtoms = nAtoms_;

   Mem.Alloc(pAtomShellOffsets, nAtoms+1);
   Mem.Alloc(pAtomBfnOffsets, nAtoms+1);
   *pAtomShellOffsets = 0;
   *pAtomBfnOffsets = 0;
   // note: this code assumes that this->shells is ordered according
   // to the atom set!
   for ( size_t iAt = 0; iAt < nAtoms; ++ iAt ){
      size_t
         iSh = pAtomShellOffsets[iAt],
         iBf = pAtomBfnOffsets[iAt];
      while ( iSh < Shells.size() && Shells[iSh].iCenter == (signed)iAt ) {
         iBf += Shells[iSh].nFn();
         ++ iSh;
      }
      pAtomShellOffsets[iAt+1] = iSh;
      pAtomBfnOffsets[iAt+1] = iBf;
      assert(iBf <= nFn());
   };
}

void FBasisSet::MakeShellOffsets( size_t *&pShellOffsets, FMemoryStack &Mem ) const
{
   Mem.Alloc(pShellOffsets, Shells.size() + 1);
   *pShellOffsets = 0;
   for ( size_t iSh = 0; iSh < Shells.size(); ++ iSh )
      pShellOffsets[iSh+1] = pShellOffsets[iSh] + Shells[iSh].nFn();
}


void FBasisSet::RebuildInterfacingData()
{
   pRawBasis = MakeRawBasis();
}


bool FBasisSet::HasSameElementBases(FBasisSet const &Other) const
{
   if (this == &Other)
      return true;
   if (this->Shells.size() != Other.Shells.size())
      return false;
   for (size_t iSh = 0; iSh != Shells.size(); ++ iSh) {
      FBasisShell const
         *pShThis = &this->Shells[iSh],
         *pShOther = &Other.Shells[iSh];
      if (pShThis->iCenter != pShOther->iCenter)
         return false;
      if (pShThis->pAs.get() != pShOther->pAs.get())
         return false;
   }
   return true;
}


FRawBasisPtr FBasisSet::MakeRawBasis()
{
   if (0) {
      FRawBasisPtr
         r(new FRawBasis());
      r->Shells.reserve(Shells.size());
      r->ShellOffsets.reserve(Shells.size() + 1);
      r->ShellOffsets.push_back(0);
      r->CenterOffsets.reserve(Shells.size() + 1);
   //    r->CenterOffsets.push_back(0);
      r->nMaxFnPerShell = 0;
      for (size_t iSh = 0; iSh < Shells.size(); ++ iSh){
         while (Shells[iSh].iCenter >= (int)r->CenterOffsets.size())
            r->CenterOffsets.push_back(r->Shells.size());
         r->Shells.push_back(Shells[iSh].MakeIrShell()); // three copies...
         size_t
            nShFn = Shells[iSh].nFn();
         r->ShellOffsets.push_back(r->ShellOffsets.back() + nShFn);
         r->nMaxFnPerShell = std::max((size_t)r->nMaxFnPerShell, nShFn);
      }
      assert(r->ShellOffsets.back() == this->nFn());
      r->CenterOffsets.push_back(Shells.size());
      return r;
   } else {
      FRawBasisPtr
         r(new FRawBasis());
      r->Shells.reserve(Shells.size());

      std::vector<size_t>
         iShCen; // center ID of each shell.
      iShCen.reserve(Shells.size());
      for (size_t iSh = 0; iSh < Shells.size(); ++ iSh){
         r->Shells.push_back(Shells[iSh].MakeIrShell()); // three copies...
         iShCen.push_back(Shells[iSh].iCenter);
      }
      r->Finalize(&iShCen[0]);
      return r;
   }
}


void FBasisSet::Finalize()
{
   RebuildInterfacingData();
}

// Make InOut = R * (InOut - d).
// The funny form is chosen because this is normally required for aligning molecules
// in space. Then d is the center of mass and R gives the main axis transformation.
void Trafo3x4(FVector3 &InOut, double const *R, double const *d)
{
   FVector3
      r = InOut - FVector3(d[0], d[1], d[2]),
      Rr(0., 0., 0.);
   for (unsigned i = 0; i < 3; ++ i)
      for (unsigned j = 0; j < 3; ++ j)
         Rr[i] += R[i + 3*j] * r[j];
   InOut = Rr;
}

void FBasisSet::Transform_3x4(double const *R, double const *d, FMemoryStack &/*Mem*/)
{
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh)
      Trafo3x4(Shells[iSh].vCenter, R, d);

   // Rebuild RawShell interface.
   // Note: Currently this is not really required as IR shells store only
   // pointers to the basis function centers.
   Finalize();
}

void FAtomSet::Transform_3x4(double const *R, double const *d, FMemoryStack &/*Mem*/)
{
   for (size_t iAt = 0; iAt < size(); ++ iAt)
      Trafo3x4(m_Atoms[iAt].vPos, R, d);
}


// #ifdef INCLUDE_OPTIONALS
#ifdef PROG_IBOVIEW
void FBasisSet::TransformMoCoeffs_3x4(FMatrixView Orbs, double const *R, FMemoryStack &Mem) const
{
   assert(Orbs.nRows == nFn());
   // construct transformation amongst solid harmonic components which
   // is induced by the rotation R.
   double
      *pAllSlcT;
   size_t
      nStride;
   ir::EvalSlcXRotationTrafo(pAllSlcT, nStride, nMaxL(), R, Mem); // <- allocates pAllSlcT on Mem.

   // apply it to all the components of the MOs.
   size_t
      iShOf = 0; // offset of the current shell within the basis.
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh) {
      FBasisShell const
         &Sh = Shells[iSh];
      uint
         l = Sh.l();
      for (size_t iCo = 0; iCo < Sh.nCo(); ++ iCo) {
         size_t
            // offset of the spherical components of the current contraction
            iMoOff = iShOf + iCo * (2*l+1);
         FMatrixView
            // original MO coefficients of current spherical set
            C = FMatrixView(&Orbs(iMoOff,0), 2*l+1, Orbs.nCols, Orbs.nRowSt, Orbs.nColSt),
            // transformed MO coefficients.
            Ct = MakeStackMatrix(2*l+1, Orbs.nCols, Mem),
            // sub-matrix of pAllSlcT for transformation amongst Slc at current l.
            SlcT = FMatrixView(pAllSlcT + (l*l)*(nStride+1), (2*l+1), (2*l+1), 1, nStride);
         Mxm(Ct, SlcT, C);
         Assign(C, Ct);

         Mem.Free(Ct.pData);
      }

      iShOf += Sh.nFn();
   }
   assert(iShOf == nFn());

   Mem.Free(pAllSlcT);
}
#endif // PROG_IBOVIEW
// #endif // INCLUDE_OPTIONALS





void FBasisSet::ExportToLibmol(std::ostream &xout, std::string const &BasisName, std::string const &Comment, FAtomSet const &Atoms) const
{
   typedef std::pair<int, int>
      FElementAmKey; // first: iElement, second: AngMom.
   typedef std::map<FElementAmKey, FAtomShell const *>
      FElementShellMap;
   FElementShellMap
      ElementShells;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh) {
      FBasisShell const
         &Sh = this->Shells[iSh];
      if (Sh.iCenter < 0 || size_t(Sh.iCenter) >= Atoms.size())
         throw std::runtime_error(fmt::format("FBasisSet::ExportToLibmol: Invalid center index {} in shell #{} (there are {} atoms)", Sh.iCenter, iSh, Atoms.size()));
      int
         iElement = Atoms[Sh.iCenter].iElement;
      if (iElement <= 0)
         // dummy atom?
         continue;
      FElementAmKey
         Key(iElement, Sh.l());
      FAtomShell const
         *pAsThis = Sh.pAs.get();
      // shell data for this combination of element and angmom already present?
      FElementShellMap::const_iterator
         itElementShell = ElementShells.find(Key);
      if (itElementShell == ElementShells.end()) {
         // no, so add the current one.
         ElementShells.insert(FElementShellMap::value_type(Key, pAsThis));
      } else {
         // already present. Make sure it's the same one we already have.
         FAtomShell const
            *pAsPrev = itElementShell->second;
         if (pAsThis != pAsPrev && *pAsThis != *pAsPrev)
            throw std::runtime_error(fmt::format("FBasisSet::ExportToLibmol: Element basis for iElement={} AngMom={} is not unique. Cannot export as basis library.", Key.first, Key.second));
         // apparently it is the same. In this case we can just ignore it.
      }
      // ^- (note: there is actually no test that all elements have the same number of AM
      // shells declared, only that the ones which are declared are the same).
   }

   // by now we got all unique element/angmom basis shells collected.
   // Now just go through them one by one, format them as text, and write them to
   // the output stream
   {
      FElementShellMap::const_iterator
         itElementShell;
      for (itElementShell = ElementShells.begin(); itElementShell != ElementShells.end(); ++ itElementShell) {
         int
            iElement = itElementShell->first.first;
         itElementShell->second->ExportToLibmol(xout, ElementNameFromNumber(iElement), BasisName, Comment);
      }
   }
}




} // namespace ct
