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

#include <stdexcept>
#include <cmath>
#include <cctype> // for std::isdigit
#include <set>
#include <iostream> // FIXME: remove this.

#include "CtAtomSet.h"
#include "CxIo.h"
#include "CtMatrix.h"
#include "CtBasisSet.h"
#include "CxPhysicalUnits.h"
#include "CtBasisLibrary.h"
#include "CxParse1.h"
#include "CxAtomParamSpec.h"
#include "CxXyzFrameIo.h"


#ifdef USE_CTINT1E_H
   #include "CtInt1e.h"
#endif

#ifdef _MSC_VER
   #define strcasecmp _stricmp    // both are non-standard functions for case-insensitive string comparisons.
#endif

// #define USE_CTINT1E_H_ALWAYS

namespace ct {

double
   g_SymmetryTolerance = 1e-10;



std::ostream &operator << (std::ostream &out, FVector3 const &Pos){
   out << fmt::format("({:7.5f} {:7.5f} {:7.5f})", Pos[0], Pos[1], Pos[2]);
//    out << format("(%16.9f %16.9f %16.9f)") %  Pos[0] % Pos[1] % Pos[2];
   return out;
}


FAtomSet::FAtomSet()
{
   m_LastEnergy = 0.;
}

FAtomSet::~FAtomSet()
{}


FScalar FAtomSet::NuclearRepulsionEnergy() const
{
   FScalar
      fEnergy = 0;
   FAtomList::const_iterator
      itAtom1, itAtom2;
   for (itAtom1 = m_Atoms.begin(); itAtom1 != m_Atoms.end(); ++itAtom1)
      for (itAtom2 = itAtom1, ++itAtom2; itAtom2 != m_Atoms.end(); ++itAtom2)
         fEnergy += itAtom1->NuclearCharge() * itAtom2->NuclearCharge() *
                  1.0/std::sqrt(DistSq(itAtom1->vPos, itAtom2->vPos));
   return fEnergy;
}


FVector3 FAtomSet::NuclearDipoleMoment(FVector3 const &ExpansionPoint) const
{
   FVector3
      Result(0,0,0);
   FAtomList::const_iterator
      itAtom;
   for (itAtom = m_Atoms.begin(); itAtom != m_Atoms.end(); ++ itAtom){
      Result += static_cast<FScalar>(itAtom->NuclearCharge()) *
         (itAtom->vPos - ExpansionPoint);
   }
   return Result;
}



// ecp-reduced total nuclear charge.
int FAtomSet::NuclearCharge() const
{
   double
      Result = 0;
   FAtomList::const_iterator
      itAtom;
   for (itAtom = m_Atoms.begin(); itAtom != m_Atoms.end(); ++itAtom)
      Result += itAtom->NuclearCharge();
   if (int(Result) != Result)
      throw std::runtime_error("FAtomSet::NuclearCharge: result was not integer.");
   if (int(Result) < 0)
      throw std::runtime_error("FAtomSet::NuclearCharge: result was negative.");
   return int(Result);
}


void FAtomSet::AddAtom(FAtom const &Atom)
{
   m_Atoms.push_back(Atom);
}

// converts a map with string keys into another map of the same format but
// which's keys have been lowercased.
template<class FValue>
std::map<std::string,FValue> MakeStringKeysLowercaseInMap(
   std::map<std::string,FValue> In )
{
   std::map<std::string,FValue>
      Out;
   typename std::map<std::string,FValue>::iterator
      it;
   _for_each(it, In)
      Out[tolower(it->first)] = it->second;
   return Out;
}




void FAtomSet::AddAtomsFromIoAtomSet(xyz_io::FIoAtomSet const &IoFrame, FBasisDescs const &DefaultBases)
{
   // in case this is the only atom data in the frame (very much the default
   // case), copy in the meta information from the IoFrame data.
   if (empty()) {
      SetLastEnergy(IoFrame.GetLastEnergy());
      SetCaption(IoFrame.GetCaption());
      SetName(IoFrame.GetName());
   }
   size_t
      // remember how many atoms we had originally so that we know
      // which are the new atoms to assign the basis sets to.
      nOriginalSize = size();
   FBasisDescs
      EmptyBasisDesc; // will assign later (with support for complex str decomposition)

   for (xyz_io::FIoAtom const &IoAtom : IoFrame) {
      FAtom
         At(IoAtom.vPos, FIoElementAndTag(IoAtom.Type), EmptyBasisDesc);
      At.vGrad = IoAtom.vGrad;
      AddAtom(At);
   }
   // assign basis sets to the atoms we added.
   AssignBases(DefaultBases, true, nOriginalSize, size());
}


void FAtomSet::AddAtomsFromXyzFile(std::string const &FileName,
   FBasisDescs const &DefaultBases,
   FXyzLoadOptions const *pOtherOptions)
{
   {
      xyz_io::FIoAtomSetPtr
         pIoFrame = xyz_io::LoadXyzFrameFromFile(FileName, pOtherOptions);
      AddAtomsFromIoAtomSet(*pIoFrame, DefaultBases);
   }
}

void FAtomSet::AddAtomsFromXyzStream(std::istream &str, std::string const &FileName,
   FBasisDescs const &DefaultBases,
   FXyzLoadOptions const *pOtherOptions)
{
   {
      xyz_io::FIoAtomSetPtr
         pIoFrame(new xyz_io::FIoAtomSet);
      pIoFrame->AddAtomsFromXyzStream(str, FileName, pOtherOptions);
      AddAtomsFromIoAtomSet(*pIoFrame, DefaultBases);
   }
}


xyz_io::FIoAtomSetPtr FAtomSet::MakeIoAtomSet() const
{
   xyz_io::FIoAtomSetPtr
      pIoFrame(new xyz_io::FIoAtomSet);
   for (FAtom const &At : *this) {
      xyz_io::FIoAtom
         IoAtom(At.vPos, At.GetElementAndTag());
      IoAtom.vGrad = At.vGrad;
      pIoFrame->AddAtom(IoAtom);
   }

   pIoFrame->SetLastEnergy(this->GetLastEnergy());
   pIoFrame->SetCaption(this->GetCaption());
   pIoFrame->SetName(this->GetName());
   return pIoFrame;
}


void LoadMultiXyz(FAtomSetList &Frames, std::string const &FileName, FBasisDescs const &DefaultBases, FAtomSet::FXyzLoadOptions const *pOtherOptions)
{
   xyz_io::FIoAtomSetList
      IoFrames;
   xyz_io::AddFramesFromMultiXyz(IoFrames, FileName, pOtherOptions);
   for (xyz_io::FIoAtomSetPtr const &pIoFrame : IoFrames) {
      FAtomSetPtr
         pAtoms(new FAtomSet());
      pAtoms->AddAtomsFromIoAtomSet(*pIoFrame, DefaultBases);
      Frames.push_back(pAtoms);
   }
}


void FAtomSet::PrintAsXyz(std::ostream &out, FXyzPrintOptions const &Options) const
{
   xyz_io::FIoAtomSetPtr
      pIoFrame = MakeIoAtomSet();
   pIoFrame->PrintAsXyz(out, "", Options);
   // ^- 2nd arg: ExtraDesc
}



double FAtomSet::GetTotalGradient() const
{
   double g = 0;
   for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt)
      g += LengthSq(m_Atoms[iAt].vGrad);
   return g;
}

bool FAtomSet::HaveGradients() const
{
   for (size_t iAt = 0; iAt < size(); ++ iAt)
      if (m_Atoms[iAt].vGrad != FVector3(0., 0., 0.))
         return true;
   return false;
}



void FAtomSet::AssignBasis(FBasisContext Context, std::string const &BasisDecl, bool AutoFix, ptrdiff_t iFirstAt_, ptrdiff_t iLastAt_)
{
   if (0) {
      std::cout << fmt::format(" FAtomSet::AssignBasis({}[{}:{}], '{}', AutoFix={})", GetBasisContextName(Context), iFirstAt_, iLastAt_, BasisDecl, AutoFix) << std::endl;
   }

   assert_rt(iFirstAt_ >= 0 && (iLastAt_ == -1 || (iLastAt_ >= 0 && size_t(iLastAt_) <= size())));
   size_t
      iAtBeg = size_t(iFirstAt_),
      iAtEnd = size();
   if (iLastAt_ >= 0)
      iAtEnd = size_t(iLastAt_);

   if (!(StartsWith(BasisDecl, "{") && EndsWith(BasisDecl, "}"))) {
      // if this is a simple direct assignment, just forward the basis definitions directly.
      for (size_t iAt = iAtBeg; iAt < iAtEnd; ++ iAt) {
         m_Atoms[iAt].AssignBasis(Context, BasisDecl, AutoFix);
      }
   } else {
      // ...otherwise do CxAtomParamSpec.h thing here, in case 'BasisDecl' is
      // a more complex basis description containing detailed atom assignments.
      // FAtomParamSpec can take apart definitions such as these:
      //
      //   "{cc-pVTZ; O:aug-cc-pVTZ; H:cc-pVDZ}"
      //   "{O:aug-cc-pVTZ; H:cc-pVDZ}"
      //   "{cc-pVTZ; [O,H]:aug-cc-pVTZ}"
      //   "{cc-pVTZ; .1:aug-cc-pVTZ}"
      //   "{cc-pVTZ; .1:aug-cc-pVTZ; [#3,#7]:def2-TZVPP}"
      //
      // See CxAtomParamSpec.h for precedence rules. The ".1" denotes tagged
      // input atoms (e.g, "C1" "C2" in input geometry, to distinguish different
      // carbon atoms---not fully supported here as of now, before CtAtomSet's
      // XyzParser is merged with CxXyzFrameIo's, though...)
      //
      // Note: use scf{print{basis-orb;basis-fit}} to see what you got.
      FAtomParamSpec_String
         BasisSpec(BasisDecl, "undefined", "{;}:", FSettingParseProxy_String());
      for (size_t iAt = iAtBeg; iAt < iAtEnd; ++ iAt) {
         FAtom
            &At = m_Atoms[iAt];
         std::string const
            &BasisForAtom = BasisSpec(iAt, At.iElement, At.Tag);
         if (0) {
            xout << fmt::format(": basis spec: {} {} --> '{}'", 1+iAt, At.GetElementAndTag(), BasisForAtom) << std::endl;
         }
         At.AssignBasis(Context, BasisForAtom, AutoFix);
      }
   }
}


void FAtomSet::AssignBases(FBasisDescs const &DefaultBases, bool AutoFix, ptrdiff_t iFirstAt, ptrdiff_t iLastAt)
{
   for (FBasisDescs::const_iterator itContext = DefaultBases.begin(); itContext != DefaultBases.end(); ++ itContext) {
      AssignBasis(itContext->first, itContext->second, AutoFix, iFirstAt, iLastAt);
   }
}


// return the number of electrons we'd typically consider as core-level
// for a given element.
unsigned GetElementNumCoreElec(unsigned iElement)
{
   static size_t const
      s_RareGasAtomicNumbers[] = { 0, 2, 10, 18, 36, 54, 86, 118 },
      s_nRareGasAtoms = sizeof(s_RareGasAtomicNumbers)/sizeof(s_RareGasAtomicNumbers[0]);
   size_t
      nCoreElec1 = 0;
   for (size_t i = 0; i < s_nRareGasAtoms; ++ i) {
      if (iElement > s_RareGasAtomicNumbers[i]) {
         nCoreElec1 = s_RareGasAtomicNumbers[i];
      } else {
         break;
      }
   }
   if (iElement >= 72 && iElement <= 86)
      nCoreElec1 += 14; // treat 4f shells as core for elements beyond Lu
   if (iElement >= 104 && iElement <= 118)
      nCoreElec1 += 14; // treat 5f shells as core for elements beyond Lr
   return nCoreElec1;
}



// return the expected number of core electrons for an atom of the given element.
// If ECP is assigned, this number is reduced by the number of ECP electrons.
size_t FAtom::nCoreElec() const
{
//    static size_t const
//       s_nRareGasAtoms = 7,
//       s_RareGasAtomicNumbers[] = { 0, 2, 10, 18, 36, 54, 86 };
//    size_t
//       nCoreElectrons1 = 0;
//    for (size_t i = 0; i < s_nRareGasAtoms; ++ i) {
//       if (this->iElement > s_RareGasAtomicNumbers[i]) {
//          nCoreElectrons1 = s_RareGasAtomicNumbers[i] - this->nEcpElec;
//       } else {
//          break;
//       }
//    }
//    return nCoreElectrons1;
   return GetElementNumCoreElec(this->iElement) - this->nEcpElec;
}



// returns ecp-reduced number of orbitals in rare gas configurations.
size_t FAtomSet::nCoreElec() const
{
   size_t
      Result = 0;
   FAtomList::const_iterator
      itAtom;
   _for_each(itAtom, m_Atoms)
      Result += itAtom->nCoreElec();
   return Result;
}


size_t FAtomSet::nEcpElec() const
{
   size_t
      Result = 0;
   FAtomList::const_iterator
      itAtom;
   _for_each(itAtom, m_Atoms)
      Result += itAtom->nEcpElec;
   return Result;
}


bool FAtomSet::IsSameMolecule(FAtomSet const &Other) const
{
   if (this == &Other)
      return true;
   if (this->size() != Other.size())
      return false;
   for (size_t iAt = 0; iAt != size(); ++ iAt) {
      FAtom const
         &AtThis = (*this)[iAt],
         &AtOther = Other[iAt];
      if (AtThis.iElement != AtOther.iElement)
         return false;
      if (AtThis.nEcpElec != AtOther.nEcpElec)
         return false;
   }
   return true;
}


bool FAtomSet::IsSameGeometry(FAtomSet const &Other, double const ThrPos) const
{
   assert(ThrPos >= 0.);
   if (this == &Other)
      return true;
   if (this->size() != Other.size())
      return false;
   for (size_t iAt = 0; iAt != size(); ++ iAt) {
      FAtom const
         &AtThis = (*this)[iAt],
         &AtOther = Other[iAt];
      if (AtThis.iElement != AtOther.iElement)
         return false;
      if (DistSq(AtThis.vPos, AtOther.vPos) > sqr(ThrPos))
         return false;
   }
   return true;
}



std::string FAtom::GetBasisDesc() const
{
   FBasisDescs::const_iterator
      itDesc = BasisDesc.find(BASIS_Orbital);
   assert(itDesc != BasisDesc.end());
   if ((this->EcpDesc.empty() || this->EcpDesc == "(none)") && this->nEcpElec == 0)
      return itDesc->second;
   else
      return fmt::format("{}, {}", itDesc->second, EcpDesc, nEcpElec);
//       return fmt::format("{}, {} ({} e- in ecp)", itDesc->second, EcpDesc, nEcpElec);
}

// std::string FAtom::GetBasisDesc() const
// {
//    FBasisDescs::const_iterator
//       itBasisOrb = BasisDesc.find(BASIS_Orbital);
//    FBasisDescs::const_iterator
//       itBasisFit = BasisDesc.find(BASIS_JkFit);
//    assert(itBasisOrb != BasisDesc.end());
//    if ((this->EcpDesc.empty() || this->EcpDesc == "(none)") && this->nEcpElec == 0)
//       return fmt::format("{} // {}", itBasisOrb->second, itBasisFit->second);
// //       return itBasisOrb->second;
//    else
//       return fmt::format("{}, {} // {}", itBasisOrb->second, EcpDesc, itBasisFit->second);
// //       return fmt::format("{}, {}", itBasisOrb->second, EcpDesc, nEcpElec);
// //       return fmt::format("{}, {} ({} e- in ecp)", itBasisOrb->second, EcpDesc, nEcpElec);
// }

// FAtom::FAtom()
// {
//    Init0();
// };

FAtom::FAtom(FVector3 const &vPos_, std::string const &ElementName_, std::string
   const &BasisDesc_, std::string const &Ecp_, bool AutoFix)
{
   BasisDesc[BASIS_Orbital] = BasisDesc_;
//    throw std::runtime_error("FAtom::FAtom() called --- string only decl. You sure that is what you want? (if yes, remove this error raise)");
   Init(vPos_, ElementName_, Ecp_, AutoFix);
}

FAtom::FAtom(FVector3 const &vPos_, std::string const &ElementName_, FBasisDescs
         const &BasisDesc_, std::string const &Ecp_, bool AutoFix)
   : BasisDesc(BasisDesc_)
{
//    std::cout << fmt::format("   FAtom::FAtom({},{:.4f},{:.4f},{:.4f}) -> BasisDesc[BASIS_JkFit] == {}", ElementName_, vPos_[0], vPos_[1], vPos_[2], BasisDesc[BASIS_JkFit]) << std::endl;
   Init(vPos_, ElementName_, Ecp_, AutoFix);
//    std::cout << fmt::format("   ...after init: BasisDesc[BASIS_JkFit] == {}", BasisDesc[BASIS_JkFit]) << std::endl;
}

static std::string const
   s_DummyMarkerThatElementAndTagAreAlreadySet;

FAtom::FAtom(FVector3 const &vPos_, FIoElementAndTag const &ElementAndTag_, FBasisDescs
         const &BasisDesc_, std::string const &Ecp_, bool AutoFix)
   : iElement(ElementAndTag_.iElement), Tag(ElementAndTag_.Tag), BasisDesc(BasisDesc_)
{
   Init(vPos_, s_DummyMarkerThatElementAndTagAreAlreadySet, Ecp_, AutoFix);
   assert(iElement == ElementAndTag_.iElement);
   assert(Tag == ElementAndTag_.Tag);
   assert_rt(iElement >= 0);
}



static bool _in_range(int i, int iMin, int iMax) { return iMin <= i && i <= iMax; }


void FAtom::Init0()
{
   iElement = 0;
   Tag = ATOMTAG_None;
   nEcpElec = 0;
   vGrad = FVector3(0.,0.,0.);
   vPos = FVector3(0.,0.,0.);

   IR_SUPPRESS_UNUSED_WARNING(_in_range);
}

void FAtom::Init(FVector3 const &vPos_, std::string const &ElementName_, std::string const &EcpName_, bool AutoFix)
{
   vPos = vPos_;
   if (!(&ElementName_ == &s_DummyMarkerThatElementAndTagAreAlreadySet)) {
      iElement = ElementNumberFromName(ElementName_.c_str());
      Tag = ATOMTAG_None;
   }
   nEcpElec = 0;
//    if (EcpName_.empty())
//       EcpDesc = "(none)";
//    else
//       EcpDesc = EcpName_;
   EcpDesc = EcpName_;
   vGrad = FVector3(0.,0.,0.);
   if (AutoFix)
      this->AutoFixBasisAssociations();
}


double FAtom::NuclearCharge() const
{
   return double(iElement - nEcpElec);
}


void FAtom::AutoFixBasisAssociations()
{
   AutoAdjustBasisAssociations(this->BasisDesc, this->EcpDesc, this->nEcpElec, this->iElement);
//       if (1) {
//          std::cout << fmt::format(" *********  {:2}  nEcpElec = {}  OrbBas = {}  MINAO = {}", pElementName(), nEcpElec, BasisDesc[BASIS_Orbital], BasisDesc[BASIS_Guess]) << std::endl;
//       }
//    Charge = iElement - nEcpElec;
}


void FAtom::AssignBasis(FBasisContext Context, std::string const &BasisName, bool AutoFix)
{
   BasisDesc[Context] = BasisName;
   if (AutoFix)
      AutoFixBasisAssociations();
}



std::ostream& operator << (std::ostream &out, FAtom const &a)
{
//    double const ToAng = 1.;
//    out << boost::format( "%+2s:   (%16.10f,%16.10f,%16.10f)    %s" )
//    out << boost::format( "%+2s:   (%10.6f,%10.6f,%10.6f)    %s" )
//    out << boost::format( "%+2s:    %10.6f %10.6f %10.6f     %s" )
//       % ElementNameFromNumber( a.iElement )
//       % (ToAng * a.vPos[0])
//       % (ToAng * a.vPos[1])
//       % (ToAng * a.vPos[2])
//       % a.GetBasisDesc();
   out << fmt::format( "{:>2}:    {:10.6f} {:10.6f} {:10.6f}     {}",
      a.GetElementName(),
      (ToAng * a.vPos[0]),
      (ToAng * a.vPos[1]),
      (ToAng * a.vPos[2]),
      a.GetBasisDesc());
   return out;
}

// std::ostream& operator << ( std::ostream &out, FAtom const &a )
// {
//    out << boost::format( "%+2s: (%+.5f,%+.5f,%+.5f) basis: '%s' ecp: '%s'" )
//       % ElementNameFromNumber( a.iElement )
//       % Distance( a.vPos.x(), false )
//       % Distance( a.vPos.y(), false )
//       % Distance( a.vPos.z(), false )
//       % a.GetBasisDesc() % a.EcpDesc;
//    return out;
// }

// std::ostream& operator << ( std::ostream &out, FAtomSet const &a )
// {
//    out << "ATOMS: " << a.m_Atoms.size() << " atoms. (distance unit: "
//       << DistanceUnit() << ")\n";
//    for( uint i = 0; i < a.m_Atoms.size(); ++ i ){
//
//       out << boost::format("%+6d  ") % i << Atoms[i] << std::endl;
//    }
//    return out;
// }

// TODO: change this to a free-standing function instead of an operator, so that we can
// supply additional arguments?
std::ostream& operator << (std::ostream &xout, FAtomSet const &Atoms)
{
   xout << "   ATOM  ELM       POS/X      POS/Y      POS/Z        BASIS\n";
   for( uint i = 0; i < Atoms.size(); ++ i ){
      xout << (fmt::format("{:>6}   ", i+1)) << Atoms[i] << std::endl;
   }
   xout << "\n   (" << Atoms.size() << " atoms;"
        << " " << Atoms.NuclearCharge() - Atoms.nCoreElec() << " valence electrons;";
   {
      size_t nEcpElec = Atoms.nEcpElec();
      if (nEcpElec != 0)
         xout << fmt::format(" {} electrons absorbed in ECPs;", nEcpElec);
   }
   xout << " distance unit: " << "A" << ")\n";


   unsigned Ops[8], Gen[3]; size_t nOps, nGen;
   Atoms.FindMirrorSymmetryOps(Ops, nOps, Gen, nGen);

   xout << "\n";
   if (1) {
      // Generators    Point group
      // (null card)    $C_1$ (i.e. no point group symmetry)
      // X (or Y or Z)    $C_s$
      // XY    $C_2$
      // XYZ    $C_i$
      // X,Y    $C_{2v}$
      // XY,Z    $C_{2h}$
      // XZ,YZ    $D_2$
      // X,Y,Z    $D_{2h}$
      //
      // 1 X Y Z XY XZ YZ XYZ
      char const *pPointGroupName = "?!";
      if (nOps == 1)
         pPointGroupName = "C1";
      else if (nOps == 8)
         pPointGroupName = "D2h";
      else if (nOps == 2 && Ops[0] == 7)
         pPointGroupName = "Ci";
      else if (nOps == 2 && (Ops[0] == 1 || Ops[0] == 2 || Ops[0] == 4))
         pPointGroupName = "Cs";
      else if (nOps == 2 && (Ops[0] == 3 || Ops[0] == 6 || Ops[0] == 5))
         pPointGroupName = "C2";
      else {
         uint n1 = 0, n2 = 0;
         for (uint i = 0; i < nOps; ++i) {
            uint nBits = ((Ops[i] & 1) >> 0) + ((Ops[i] & 2) >> 1) + ((Ops[i] & 4) >> 2);
            if (nBits == 1)
               n1 += 1;
            if (nBits == 2)
               n2 += 1;
         }
         if (n1 == 2 && n2 == 1)
            pPointGroupName = "C2v"; // id X Y XY
         if (n1 == 1 && n2 == 1)
            pPointGroupName = "C2h"; // id XY Z XYZ
         if (n1 == 0 && n2 == 3)
            pPointGroupName = "D2"; // id XZ YZ XY
      }

      //       if ( Verbosity >= 0 ) {
      if (1) {
         xout << (fmt::format(" {:<31}{}", "Point group:", pPointGroupName)) << std::endl;
         xout << fmt::format(" {:<31}", "Symmetry ops:");
         for (uint iOp = 0; iOp < nOps; ++iOp) {
            if (iOp != 0)
               xout << "  ";
            for (uint i = 0; i < 3; ++i)
               if ((Ops[iOp] & (1 << i)) != 0)
                  xout << "XYZ"[i];
            for (uint i = 0; i < nGen; ++i)
               if (iOp == Gen[i])
                  xout << "*";
            if (Ops[iOp] == 0)
               xout << "id";
         }
         xout << std::endl;
      }
   }

   return xout;
}









#ifdef USE_CTINT1E_H
void FAtomSet::Add1eIntegralMatrix(FMatrixView &Dest,
      FBasisSet const &RowBasis, FBasisSet const &ColBasis,
      FDoublettIntegralFactory &IntFactory, double Factor, FMemoryStack &Mem) const
{
//    xout << "\n\nFAtomSet::Add1eIntegralMatrix WAS CALLED!!\n\n" << std::endl;
   bool
      MatrixSymmetric = (&RowBasis == &ColBasis);
//    MatrixSymmetric = false; // FIXME: remove this.

   assert(Dest.nRows == RowBasis.nFn());
   assert(Dest.nCols == ColBasis.nFn());

   FShellDoublettIntegral
      IntResult;
   size_t
      StartFnA = 0, StartFnB; // starting indices of the functions of the
                        // respective shells.
   for (size_t nShellA = 0; nShellA < RowBasis.Shells.size(); ++ nShellA) {
      FBasisShell const
         *pShellA = &RowBasis.Shells[nShellA];
      StartFnB = (!MatrixSymmetric)? 0 : StartFnA;
      for (size_t nShellB = (!MatrixSymmetric)? 0 : nShellA;
         nShellB < ColBasis.Shells.size(); ++ nShellB)
      {
         FBasisShell const
            *pShellB = &ColBasis.Shells[nShellB];
         IntFactory.EvalDoublett( IntResult, pShellA, pShellB, Mem );

         assert(IntResult.nSizeA == pShellA->nFn());
         assert(IntResult.nSizeB == pShellB->nFn());
         assert(StartFnA + IntResult.nSizeA <= RowBasis.nFn());
         assert(StartFnB + IntResult.nSizeB <= ColBasis.nFn());

         // fill data we gathered into the matrices
         for (size_t i_ = 0; i_ < IntResult.nSizeA; ++ i_) {
            for (size_t j_ = 0; j_ < IntResult.nSizeB; ++ j_) {
               size_t
                  i = StartFnA + i_,
                  j = StartFnB + j_;
               FScalar const
                  &r = IntResult.Int( i_, j_ );
               Dest(i,j) += Factor * r;
               if (MatrixSymmetric && nShellA != nShellB)
                  Dest(j,i) += Factor * r;
            }
         }
         StartFnB += pShellB->nFn();
      }
      StartFnA += pShellA->nFn();
   }
   assert(StartFnA == RowBasis.nFn());
}

void FAtomSet::Make1eIntegralMatrixMultiComp(double *pDest,
      FBasisSet const &RowBasis, FBasisSet const &ColBasis,
      FDoublettIntegralFactory &IntFactory, double Factor, FMemoryStack &Mem_) const
{
   bool
      MatrixSymmetric = (&RowBasis == &ColBasis);
//    MatrixSymmetric = false; // FIXME: remove this.

//    FShellDoublettIntegral
//       IntResult;
   size_t
      nComp = IntFactory.nComp();
   size_t
      nFnBasisA = RowBasis.nFn(), nFnBasisB = ColBasis.nFn();
   {
      FMemoryStackArray MemStacks(Mem_);
      #pragma omp parallel for schedule(dynamic)
      // ^- this doesn't help at all...
//       for ( size_t nShellA = 0; nShellA < RowBasis.Shells.size(); ++ nShellA ){
      for (int nShellA_ = 0; nShellA_ < (int)RowBasis.Shells.size(); ++ nShellA_){
         size_t nShellA = size_t(nShellA_); // <- OMP 2.0 only allows "int" as loop variable, not uint or size_t or ptrdiff_t.
         FShellDoublettIntegral
            IntResult;
         FMemoryStack &Mem = MemStacks.GetStackOfThread();
         FBasisShell const
            *pShellA = &RowBasis.Shells[nShellA];
         for (size_t nShellB = (!MatrixSymmetric)? 0 : nShellA;
            nShellB < ColBasis.Shells.size(); ++ nShellB)
         {
            FBasisShell const
               *pShellB = &ColBasis.Shells[nShellB];
            IntFactory.EvalDoublett(IntResult, pShellA, pShellB, Mem);

            size_t
               StartFnA = RowBasis.pRawBasis->iFn(nShellA),
               StartFnB = ColBasis.pRawBasis->iFn(nShellB); // starting indices of the functions of the


            // fill data we gathered into the matrices
            for (size_t iComp = 0; iComp < nComp; ++ iComp) {
               FMatrixView
                  Dest(pDest + nFnBasisA*nFnBasisB * iComp, nFnBasisA, nFnBasisB);
               for (size_t i_ = 0; i_ < IntResult.nSizeA; ++ i_) {
                  for (size_t j_ = 0; j_ < IntResult.nSizeB; ++ j_)
                  {
                     size_t
                        i = StartFnA + i_,
                        j = StartFnB + j_;
                     FScalar const
                        &r = IntResult.Int(i_, j_, iComp);
                     Dest(i,j) = Factor * r;
                     if (MatrixSymmetric && nShellA != nShellB)
                        Dest(j,i) = Factor * r;
                  }
               }
            }

         }
      }
   }
}


#endif // USE_CTINT1E_H

void FAtomSet::MakeCoreHamiltonMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const
{
   Out.Clear();
   AddKineticMatrix(Out, RowBasis, ColBasis, Mem);
   AddNuclearAttractionMatrix(Out, RowBasis, ColBasis, Mem);

   // TODO: apply other CoreH patches from Environment object etc.

//    AddDipoleMatrix( Out, RowBasis, ColBasis, FVector3(0.0,0.0,1e-5), FVector3(0,0,0), Mem );
}


void FAtomSet::AccCoreHamiltonDeriv1(FMatrixView &GradXyz, FMatrixView const &Rdm, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const
{
   if (GradXyz.nRowSt != 1 || GradXyz.nColSt != 3 || GradXyz.nRows != 3 || GradXyz.nCols != this->size())
      throw std::runtime_error("Output gradient should be matrix of format 3 x nAtoms, without strides.");
   FKrn2iPtr
      pCorePotentialKernel = MakeNuclearPotentialKernel(),
      pKineticEnergyKernel = MakeKineticEnergyKernel();

//    dynamic_cast<FKrn2i_PointMultipoles*>(pCorePotentialKernel.get())->SetDerivativeFlags(FKrn2i_PointMultipoles::DERIVFLAG_IncludePointDerivsP);
   AccGradient2ix(GradXyz.pData, Rdm, *RowBasis.pRawBasis, *ColBasis.pRawBasis, *pCorePotentialKernel, Mem);
   AccGradient2ix(GradXyz.pData, Rdm, *RowBasis.pRawBasis, *ColBasis.pRawBasis, *pKineticEnergyKernel, Mem);
}

void FAtomSet::AccNuclearRepulsionGradient(FMatrixView &GradXyz, FMemoryStack &Mem, double Factor) const
{
   // add derivative of nuclear-nuclear repulsion (i.e., of the contribution of this->NuclearRepulsionEnergy()).
   for (size_t iAtA = 0; iAtA < size(); ++ iAtA)
      for (size_t iAtB = iAtA + 1; iAtB < size(); ++ iAtB) {
         FAtom const
            *pAtA = &m_Atoms[iAtA],
            *pAtB = &m_Atoms[iAtB];
         // grad pow(rSqAB, -.5) = -.5 * pow(rSqAB, -1.5) * grad rSqAB
         // grad rSqAB = 2 * DistAB
         FVector3
            vDist = pAtA->vPos - pAtB->vPos;
         double
            rSqAB = vDist.LengthSq();
         double
            f = -.5 * 2 * pAtA->NuclearCharge() * pAtB->NuclearCharge() *
                  std::pow(rSqAB, -1.5) * Factor;
         for (uint ixyz = 0; ixyz < 3; ++ ixyz) {
            GradXyz(ixyz, iAtA) += f * vDist[ixyz];
            GradXyz(ixyz, iAtB) -= f * vDist[ixyz];
         }
      }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}





void FAtomSet::MakeOverlapMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const
{
#ifdef USE_CTINT1E_H_ALWAYS
   FDoublettIntegralFactoryOverlap
      IntFactory; // can't bind temporary to non-const reference.
   Out.Clear();
   Add1eIntegralMatrix(Out, RowBasis, ColBasis, IntFactory, 1.0, Mem);
#else
//    MakeIntMatrix(Out, RowBasis, ColBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel,0), Mem, 1.0, false);
   MakeIntMatrix(Out, RowBasis, ColBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
#endif
}

void FAtomSet::AddKineticMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const
{
#ifdef USE_CTINT1E_H_ALWAYS
   FDoublettIntegralFactoryKineticTerm
      IntFactory;
   Out.Clear();
   Add1eIntegralMatrix(Out, RowBasis, ColBasis, IntFactory, 1.0, Mem);
#else
   FKrn2iPtr
      pKineticEnergyKernel = MakeKineticEnergyKernel();
   MakeIntMatrix(Out, RowBasis, ColBasis, *pKineticEnergyKernel, Mem, 1., false);
//    MakeIntMatrix(Out, RowBasis, ColBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel,1), Mem, -.5, false);
#endif
}




FKrn2iPtr FAtomSet::MakeNuclearPotentialKernel(unsigned PotentialFlags) const
{
   FKrn2iMetaPtr
      pMetaKernel = new FKrn2i_Meta();
   if (PotentialFlags & POTENTIAL_IncludeNuclearElectrostatic) {
      // make Coulomb part of core attraction matrix.
      std::vector<FKrn2i_PointMultipoles::FPoint>
         MultiPoles(m_Atoms.size());
      std::vector<double>
         Coeffs(m_Atoms.size());
      for (uint iAt = 0; iAt < m_Atoms.size(); ++ iAt) {
         MultiPoles[iAt].vPos = m_Atoms[iAt].vPos;
         MultiPoles[iAt].l = 0;
         MultiPoles[iAt].iCenter = static_cast<ptrdiff_t>(iAt);
         Coeffs[iAt] = -m_Atoms[iAt].NuclearCharge();
         // ^- the minus is there because we effectively calculate with
         //    positive electron charge.
   //       std::cout << fmt::format(" .. wmme nucl pot krn: #{:2}  Q = {:4}  l = {}  xyz = {:10.5f} {:10.5f} {:10.5f}\n", iAt, Coeffs[iAt], MultiPoles[iAt].l = 0, MultiPoles[iAt].vPos[0], MultiPoles[iAt].vPos[1], MultiPoles[iAt].vPos[2]);
      }
      FKrn2iPtr
         pMultipoleKernel(new FKrn2i_PointMultipoles(&ir::g_IrCoulombKernel, &MultiPoles[0], &Coeffs[0], MultiPoles.size()));
      if (!NeedEcpKernel() && PotentialFlags == (POTENTIAL_IncludeNuclearElectrostatic | POTENTIAL_IncludeEcp))
         return pMultipoleKernel;
      // ^- FIXME: remove this when other core attraction patches are introduced (e.g., finite fields)

      if (0) {
         std::cout << "\n\n =======================================\n WARNING WARNING WARNING WARNING WARNING\n\n TURNING OFF COULOMB CONTRIBUTION TO VNUC KERNEL--ONLY ECP REMAIN!!\n\n WARNING WARNING WARNING WARNING WARNING\n =======================================\n\n" << std::endl;
      } else {
         pMetaKernel->AddKernel(pMultipoleKernel, 1.0);
      }
   }
   if (PotentialFlags & POTENTIAL_IncludeEcp) {
      // add kernels for ECP part of core attraction matrix, if needed. Make one for each ECP'd core
      // for a start (later on, at least the local part of ECPs could be absorbed into meta kernels which
      // are evaluated together with the Coulomb parts. For the non-local parts I am not sure atm.).
      for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt) {
         FAtom const
            &At = m_Atoms[iAt];
         if (At.nEcpElec != 0) {
#ifdef IR_ECP
            ir::FAtomEcp const
               *pIrEcp = g_BasisSetLibrary.LoadEcp(At.iElement, At.EcpDesc, At.vPos, int(iAt));

            // PrintEcpData(std::cout, pIrEcp, &At);

            if (pIrEcp->fElecAbsorbed != double(At.nEcpElec))
               throw std::runtime_error("Actual ECP electron absorbance does not match expected one.");
            if (double(At.iElement) != double(At.NuclearCharge()) + double(At.nEcpElec))
               throw std::runtime_error("Atomic charge not consistent with (iElement - nEcpElec)");
            FKrn2iPtr pEcpKrn = new FKrn2i_1CenterEcp(At.vPos, ptrdiff_t(iAt), pIrEcp);
            pMetaKernel->AddKernel(pEcpKrn, 1.0);
#else
            throw std::runtime_error(fmt::format("Atom {} has declared a ECP{}, but this version of IR was compiled without ECP support.", GetAtomLabel(iAt), At.nEcpElec));
#endif // IR_ECP
         }
      }
   }
   return pMetaKernel;
}

FKrn2iPtr FAtomSet::MakeKineticEnergyKernel() const
{
   return FKrn2iPtr(new FKrn2i_Direct(&ir::g_IrOverlapKernel,1, -.5));
}


bool FAtomSet::NeedEcpKernel() const
{
   for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt)
      if (m_Atoms[iAt].nEcpElec != 0)
         return true;
   return false;
}


void FAtomSet::AddNuclearAttractionMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem, unsigned PotentialFlags) const
{
#ifdef USE_CTINT1E_H_ALWAYS
   if (NeedEcpKernel())
      throw std::runtime_error("Sorry, CtInt1e has no ECP integrals. Cannot make nuclear attraction matrix. ");

//    Out.Clear(); // remove this!
   typedef FDoublettIntegralFactoryFieldTerms
      FNucFactory;
   FNucFactory::FPointChargeList
      PointCharges;
   FAtomSet::FAtomList::const_iterator
      itAtom;
   // make one point charge for every atom. charge: ecp-reduced core charge.
   _for_each(itAtom, m_Atoms)
      PointCharges.push_back(FNucFactory::FPointCharge(
            itAtom->vPos, static_cast<FScalar>(itAtom->Charge)));
   FNucFactory
      NuclearAttractionIntegrator(TVector3<uint>( 0, 0, 0 ), PointCharges);
   Add1eIntegralMatrix(Out, RowBasis, ColBasis, NuclearAttractionIntegrator, 1.0, Mem);
//    Out.Print(xout, "VNUC/CtInt1e");
#else
//    Out.Clear(); // remove this!
   FKrn2iPtr
      pCorePotentialKernel = FAtomSet::MakeNuclearPotentialKernel(PotentialFlags);
   MakeIntMatrix(Out, RowBasis, ColBasis, *pCorePotentialKernel, Mem, 1., true);
//    Out.Print(xout, "VNUC/MakeIntMatrix");
//    throw std::runtime_error(":/");
#endif // USE_CTINT1E_H
}


void FAtomSet::AddDipoleMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis,
   FVector3 const &Direction, FVector3 const &ExpansionPoint, FMemoryStack &Mem) const
{
#ifdef USE_CTINT1E_H
/*   FVector3
      DirNormalized = (-1.0/std::sqrt(Dot(Direction,Direction))) * Direction;*/
      // ^- -1.0: Electrons are very negatively charged. this is actually
      // the charge factor with the elementary charge (which is 1 in a.u.)
   //xout << boost::format( "Dipole Direction: %f %f %f" )
   //   % DirNormalized.x() % DirNormalized.y() % DirNormalized.z() << std::endl;

   //  make dipole matrices into the three cartesian directions.
   for (unsigned CartComp = 0; CartComp < 3; ++ CartComp){
      if (fabs(Direction[CartComp]) < 1e-15)
         continue;

      FDoublettIntegralFactoryMultipoleMoment::FMoment
         CartMoment(0, 0, 0);
      // this specifies the derivative power in cart direction Cart.
      CartMoment[CartComp] = 1;

      FDoublettIntegralFactoryMultipoleMoment
         IntFactory(CartMoment, ExpansionPoint);
      // add dot product component of current direction times dipole matrix.
      Add1eIntegralMatrix(Out, RowBasis, ColBasis, IntFactory,
         Direction[CartComp], Mem);
      //xout << "DipoleMomentOp" << CartComp << ":"
      //    << ResultDir[CartComp].View( 0, 7, 0, std::min( 40u, ResultDir[CartComp].SizeY() ) );
   }
#else
   throw std::runtime_error("dipole integrals not implemented here.");
#endif
}


void FAtomSet::MakeCartesianMomentMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis,
   TVector3<unsigned> const &CartMoment_, FVector3 const &ExpansionPoint, FMemoryStack &Mem) const
{
#ifdef USE_CTINT1E_H
   Out.Clear();

   FDoublettIntegralFactoryMultipoleMoment::FMoment
      CartMoment(CartMoment_[0], CartMoment_[1], CartMoment_[2]);
   // this specifies the derivative power in cart direction Cart.
   FDoublettIntegralFactoryMultipoleMoment
      IntFactory(CartMoment, ExpansionPoint);
   // add to output matrix.
   Add1eIntegralMatrix(Out, RowBasis, ColBasis, IntFactory, 1., Mem);
#else
   throw std::runtime_error("general cartesian moment integrals not implemented here.");
#endif
}


bool IsEquivalent(FAtom const &A, FVector3 const &vA, FAtom const &B)
{
   if (A.iElement != B.iElement)
      return false;
   if (A.Tag != B.Tag)
      return false;
   if (A.NuclearCharge() != B.NuclearCharge())
      return false;
   if (DistSq(vA, B.vPos) > sqr(g_SymmetryTolerance))
      return false;
   if (A.BasisDesc != B.BasisDesc)
      return false;
   if (A.EcpDesc != B.EcpDesc)
      return false;
   return true;
}

FVector3 ApplySymOp(FVector3 const &In, unsigned Op){
   FVector3
      MirrorPos = In;
   for (unsigned iXyz = 0; iXyz < 3; ++ iXyz)
      if (Op & (1<<iXyz))
         MirrorPos[iXyz] = -MirrorPos[iXyz];
   return MirrorPos;
}


void FAtomSet::FindEquivalentAtoms(FAtomSymGroup *&pGroups, size_t &nGroups, FMemoryStack &Mem) const
{
   Mem.Alloc(pGroups, m_Atoms.size()); // actual data will be less.
   bool
      *pDone;
   Mem.ClearAlloc(pDone, m_Atoms.size());

   unsigned Ops[8], Gen[3]; size_t nOps, nGen;
   FindMirrorSymmetryOps(Ops, nOps, Gen, nGen);


   FAtomSymGroup
      *pCurGroup = pGroups;
   for (unsigned iAt0 = 0; iAt0 < m_Atoms.size(); ++ iAt0) {
      if (pDone[iAt0])
         continue; // already seen as part of another atom symmetry group.

      pCurGroup->nEqiv = 0;

      // look up other equivalent atoms and add them to the group.
      for (unsigned iOp = 0; iOp < nOps; ++ iOp){
         FVector3
            MirrorPos = ApplySymOp(m_Atoms[iAt0].vPos, Ops[iOp]);
         for (size_t iAt1 = iAt0; iAt1 < m_Atoms.size(); ++ iAt1)
            if (!pDone[iAt1] &&
                 IsEquivalent(m_Atoms[iAt0], MirrorPos, m_Atoms[iAt1])) {
               pCurGroup->iEqiv[pCurGroup->nEqiv] = iAt1;
               pCurGroup->SymOp[pCurGroup->nEqiv] = Ops[iOp];
               pCurGroup->nEqiv += 1;
               assert(pCurGroup->nEqiv <= 8);
               pDone[iAt1] = true;
            }
      }

      // find out directions in which current atoms lie
      // on the mirror plane.
      pCurGroup->Stabilizer = 0;
      for (unsigned iXyz = 0; iXyz < 3; ++ iXyz){
         if (sqr(m_Atoms[iAt0].vPos[iXyz]) <= sqr(g_SymmetryTolerance))
            pCurGroup->Stabilizer |= (1 << iXyz);
      }

      ++ pCurGroup;
   };
   nGroups = pCurGroup - pGroups;
   assert(nGroups <= m_Atoms.size());

   Mem.Free(pDone);
   // still on stack: pGroups[1..Atoms.size()].
}



void FAtomSet::FindMirrorSymmetryOps(unsigned Ops[8], size_t &nOps, unsigned Gen[3], size_t &nGen) const
{
   nOps = 0;

   if (1) {
      // This disables symmetry everywhere.
      nOps = 1;
      nGen = 0;
      Ops[0] = 0;
      nOps = 1;
      return;
   }

   unsigned
      AllOps[8] = {0, 1, 2, 4, 3, 5, 6, 7}; // simple ones first: id X Y Z XY XZ YZ XYZ

   for (size_t iOp_ = 0; iOp_ < 8; ++ iOp_){
      unsigned
         Op = AllOps[iOp_];
      bool
         Keep = true;

      for (size_t iAt0 = 0; iAt0 < m_Atoms.size(); ++ iAt0){
         bool
            FoundEquiv = false;
         FVector3
            MirrorPos = ApplySymOp(m_Atoms[iAt0].vPos, Op);

         for (size_t iAt1 = 0; iAt1 < m_Atoms.size(); ++ iAt1)
            if (IsEquivalent(m_Atoms[iAt0], MirrorPos, m_Atoms[iAt1])) {
               FoundEquiv = true;
               break;
            }

         // no atom equivalent under symmetry xform found for iAt0.
         if (!FoundEquiv) {
            Keep = false;
            break;
         }
      }

      if (Keep) {
         Ops[nOps] = Op;
         nOps += 1;
      }
   }

   nGen = 0;
   for (size_t i = 0; i < 3; ++ i)
      Gen[i] = 0;

   for (size_t iBaseOp = 0; iBaseOp < nOps; ++ iBaseOp){
      unsigned
         BaseOp = Ops[iBaseOp];
      bool
         Redundant = false || BaseOp == 0;
      // Can the current operator can be formed as some combination
      // of the previous generators?...
      for (unsigned iPat = 0; iPat < (1u<<nGen); ++ iPat){
         unsigned
            TrialOp = 0;
         for (size_t iGen = 0; iGen < nGen; ++ iGen)
            if (iPat & (1u << iGen)) {
               // generator iGen activated in trial pattern.
               TrialOp ^= Gen[iGen];
            }
         if (TrialOp == BaseOp) {
            Redundant = true;
            break;
         }
      }
      if (Redundant == true)
         continue;

      // ...nope. Include it as new generator.
      Gen[nGen] = BaseOp;
      nGen += 1;
   }
}


std::string FAtom::GetElementAndTag() const {
   return CombineElementAndAtomTag(iElement, Tag);
}


char const *FAtom::pElementName() const {
   return ElementNameFromNumber(this->iElement);
}


std::string FAtom::GetAtomLabel(int iAtom, int iAtomWidth) const {
   char const *pFmt = "{} {:<2}";
   switch (iAtomWidth) {
      case 0: pFmt = "{} {}"; break; // smallest fit for both number and name.
      case 1: pFmt = "{} {:<2}"; break;
      case 2: pFmt = "{:>2} {:<2}"; break;
      case 3: pFmt = "{:>3} {:<2}"; break;
      case 4: pFmt = "{:>4} {:<2}"; break;
      case 5: pFmt = "{:>5} {:<2}"; break;
      default: break;
   }

   return fmt::format(pFmt, iAtom+1, ElementNameFromNumber(iElement));
}

std::string FAtomSet::GetAtomLabel(int iAtom, int iAtomWidth) const
{
   if (iAtomWidth < 0) {
      // decide on a width of the atom number field which can align all atom indices.
      iAtomWidth = 1;
      if (size() >= 10) iAtomWidth = 2;
      if (size() >= 100) iAtomWidth = 3;
      if (size() >= 1000) iAtomWidth = 4;
      if (size() >= 10000) iAtomWidth = 5;
   }
   if (iAtom < 0 || size_t(iAtom) > size())
      throw std::runtime_error(fmt::format("FAtomSet::GetAtomLabel: Requested label for non-existent atom index {} (nAt = {})", iAtom, size()));
   return m_Atoms[iAtom].GetAtomLabel(iAtom, iAtomWidth);
}




FAtomSetProperty::FAtomSetProperty(std::string Desc_, unsigned nCenters_, unsigned nStates_, unsigned nComps_, FAtomSetProperty::FDataArray const &Data_)
   : m_Desc(Desc_), m_nCenters(nCenters_), m_nStates(nStates_), m_nComps(nComps_), m_Data(Data_)
{
   SanityCheck();
}


void FAtomSetProperty::SanityCheck() const
{
   if (m_nComps * m_nStates * m_nCenters != m_Data.size())
      throw std::runtime_error(fmt::format("FAtomSetProperty::SanityCheck(): internal program error for property '{}': data size {} incompatible with stored dimensions (nComp = {}, nStates = {}, nCen = {})", m_Desc, m_Data.size(), m_nComps, m_nStates, m_nCenters));
}


double FAtomSetProperty::Get(unsigned iComp, unsigned iState, unsigned iCenter) const
{
   SanityCheck();
   if (iComp >= m_nComps || iCenter >= m_nCenters || iState >= m_nStates)
      throw std::runtime_error(fmt::format("FAtomSetProperty::Get(): property '{}': reference to non-existing component iComp = {} iCen = {} (have nComp = {}, nStates = {}, nCen = {})", m_Desc, iComp, iState, iCenter, m_nComps, m_nStates, m_nCenters));
   return m_Data[iComp + m_nComps * (iState + m_nStates * iCenter)];
}


size_t FAtomSet::nStoredPropertyEntries() const
{
   size_t nDataSize = 0;
   for (FPropertyList::const_iterator itProp = m_Properties.begin(); itProp != m_Properties.end(); ++ itProp)
      nDataSize += (*itProp)->m_Data.size();
   return nDataSize;
}


void FAtomSet::CollectPropertyData(FPropertyDataEntryList &PropertyData) const
{
   PropertyData.clear();
   PropertyData.reserve(nStoredPropertyEntries());
   for (size_t iProp = 0; iProp != m_Properties.size(); ++iProp) {
      FAtomSetProperty const
         *pProp = &*m_Properties[iProp];
      for (size_t iCenter = 0; iCenter != pProp->m_nCenters; ++ iCenter) {
         for (size_t iState = 0; iState != pProp->m_nStates; ++ iState) {
            for (size_t iComp = 0; iComp != pProp->m_nComps; ++ iComp) {
               double fValue = pProp->Get(iComp, iState, iCenter);
               FPropertyDataEntry de = {iProp, iCenter, iState, iComp, fValue};
               PropertyData.push_back(de);
            }
         }
      }
   }
}




// fix up phase a vector such it has its largest element (in absolute terms) positive.
static void FixVectorPhase(double *p, size_t n)
{
   if (n == 0)
      return;
   double
      fMax = 0.;
   size_t
      iMax = 0;
   for (size_t i = 0; i < n; ++ i) {
      double fAbs = std::abs(p[i]);
      if (fAbs > fMax) {
         iMax = i;
         fMax = fAbs;
      }
   }
   if (p[iMax] < 0)
      for (size_t i = 0; i < n; ++ i)
         p[i] *= -1;
}


void Transpose3x3(double *R) {
   std::swap(R[0 + 3*1], R[1 + 3*0]);
   std::swap(R[0 + 3*2], R[2 + 3*0]);
   std::swap(R[1 + 3*2], R[2 + 3*1]);
}


void FindAtomSetAligningTrafo(double pR[9], double pD[3], FVector3 const *pAtPos, double const *pAtMass, FVector3 const *pAtPosLast, size_t nAt, FMemoryStack &Mem)
{
#define GET_AT_MASS(iAt) (pAtMass? pAtMass[iAt] : 1.)
   double
      TotalMass = 0,
      pIm[9] = {0}; // inertial tensor.
   FVector3
      vMassCenter(0.,0.,0.);
   // compute total mass and center of mass
   for (size_t iAt = 0; iAt < nAt; ++ iAt) {
      TotalMass += GET_AT_MASS(iAt);
      vMassCenter += GET_AT_MASS(iAt) * pAtPos[iAt];
   }
   vMassCenter /= TotalMass;

   // copy center translation
   for (size_t i = 0; i < 3; ++ i)
      pD[i] = vMassCenter[i];

   if (pAtPosLast == 0) {
      // no reference frame given.

      // compute inertial tensor around the center of mass.
      for (size_t iAt = 0; iAt < nAt; ++ iAt) {
         FVector3 const
            vAtPos = pAtPos[iAt] - vMassCenter;
         double
            Rsq = Dot(vAtPos,vAtPos);
         for (size_t i = 0; i < 3; ++ i)
            for (size_t j = 0; j < 3; ++ j)
               pIm[i + 3*j] += GET_AT_MASS(iAt) * (((i==j)? 1. : 0.)*Rsq - vAtPos[i]*vAtPos[j]);
   //             I += ElementMasses[o.iElement]*(np.eye(3)*np.dot(Xyz,Xyz) - np.outer(Xyz,Xyz))
               // ^- hm, the Rsq term is diagonal; this does not do anything, does it?
      }

      // align along axes of inertia--largest moment along z (this way flat
      // molecules will be aligned in the x/y plane).
      for (size_t i = 0; i < 9; ++ i)
         pR[i] = +pIm[i];
      double
         pEw[3];
      FMatrixView
         mR(pR, 3, 3);
      Diagonalize(mR, pEw, Mem);
      // fix up phases of the eigenvectors, such that they have their largest element positive.
      for (size_t iEv = 0; iEv < 3; ++ iEv)
         FixVectorPhase(&mR(0,iEv), 3);

      // check the determinant of the transformation to see if it includes an
      // inversion. If yes, invert one of the axes. We'd get the wrong stero-isomers
      // otherwise...
      if (1) {
         FVector3
            v0(mR(0,0), mR(1,0), mR(2,0)),
            v1(mR(0,1), mR(1,1), mR(2,1)),
            v2(mR(0,2), mR(1,2), mR(2,2));
         if (Dot(v0, Cross(v1,v2)) < 0.) {
            for (size_t i = 0; i < 3; ++ i)
               mR(i,2) *= -1.;
         }
      }
   } else {
      // compute overlap matrix between current and last frame (assuming that the last one
      // has already been moved into its own center of mass)
      for (size_t iAt = 0; iAt < nAt; ++ iAt) {
         FVector3 const
            vAtPosI = pAtPos[iAt] - vMassCenter,
            vAtPosJ = pAtPosLast[iAt];
         for (size_t i = 0; i < 3; ++ i)
            for (size_t j = 0; j < 3; ++ j)
               pIm[i + 3*j] += GET_AT_MASS(iAt) * vAtPosI[i] * vAtPosJ[j];
      }
      double
         pU[9], pVt[9], pSig[3];
      FMatrixView
         mR(pR, 3, 3),
         mU(pU, 3, 3),
         mI(pIm, 3, 3),
         mVt(pVt, 3, 3);
      ComputeSvd(mU, pSig, mVt, mI, Mem);
      Mxm(mR, mU, mVt);
      // ^- I think it should be the other way around, but there is a strange Transpose3x3
      //    in the actual transformation routine...
   }

   bool
      t3x3 = true; // <- I wonder why I put this here.
   if (t3x3)
      Transpose3x3(pR);
#undef GET_AT_MASS
}






} // namespace ct


#include "CxRawAtom.h" // only for converting FAtomSet objects to these guys

namespace ct {
   FRawAtomList AsRawAtoms(FAtomSet const &Atoms)
   {
      FRawAtomList
         r;
      r.reserve(Atoms.size());
      for (size_t iAt = 0; iAt != Atoms.size(); ++ iAt) {
         r.push_back(FRawAtom(Atoms[iAt].vPos, Atoms[iAt].iElement, Atoms[iAt].Tag));
      }
      return r;
   }
}


// kate: indent-width 3;
