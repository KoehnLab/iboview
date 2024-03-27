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

#include "CtIao.h"

#include "Ir.h"
#include "CxIo.h"
#include "CxTiming.h"
#include "CxPhysicalUnits.h"


namespace ct {

// symmetrically orthogonalize a set of vectors with respect to overlap matrix S.
// void SymOrth(FMatrixView Orbs, FMatrixView const S, FMemoryStack &Mem, double fThrAbs=1e-15, double fThrRel=1e-15)
// {
//    FStackMatrix
//       SmhOrb(Orbs.nCols, Orbs.nCols, &Mem),
//       T1(Orbs.nRows, Orbs.nCols, &Mem);
//    ChainMxm(SmhOrb, Transpose(Orbs), S, Orbs, Mem);
//    CalcSmhMatrix(SmhOrb, Mem, FSmhOptions(1e-15,1e-15,0,0));
//    Mxm(T1, Orbs, SmhOrb);
//    Assign(Orbs, T1);
// }



// NOTE: CIb must be provided from outside.
//
// Description:
// - This implements the "IAO/2014" variant of the IAO construction in a
//   numerically highly efficient implementation. Both the formal definition
//   (and how it differs from the IAO/2013 defintion in the 2013 IAO paper) and
//   this algorithm were originally provided in my Python reference
//   implementation ibo_ref.py
//
// - I described exactly what this does, and how, in the SI of
//
//    Senjean, Sen, Repisky, Knizia, Visscher:
//    "Generalization of Intrinsic Orbitals to Kramers-Paired Quaternion Spinors, Molecular Fragments, and Valence Virtual Spinors"
//    J. Chem. Theory Comput. 2021, 17, pp 1337-1354
//    http://doi.org/10.1021/acs.jctc.0c00964
//
//   The article is technically on relativistic IAOs, but the SI (can be
//   accessed without subscription) describes the 2014 IAO construction
//   exactly, including the derivation of the current numerical algorithm.
void MakeIaoBasisNew(FMatrixView CIb, FBasisSet *pMinBasis,
   FAtomSet const &Atoms, FBasisSet const *pOrbBasis,
   FMatrixView const &COcc_, FMatrixView const &S1, FMatrixView const &S1cd, unsigned Flags,
   FMemoryStack &Mem)
{
   size_t
      nBf = pOrbBasis->nFn(), // number of 'real basis' functions
      nOcc = COcc_.nCols,
      nIb = pMinBasis->nFn(); // actual AOs in minimal basis
//    CIb = MakeStackMatrix(nBf, nIb, Mem);
   assert(CIb.nRows == nBf && CIb.nCols == nIb && nIb >= nOcc);

   FMatrixView
      COcc;
   if (0 == (Flags&IAO_NormalizeInput)) {
      COcc = COcc_;
   } else {
      // input has absorbed occupation numbers. Remove this.
      // Note: may also be relevant in case of input vectors with limited numerical
      // precision (e.g., stuff read from orbital files). For this reason here a
      // full orthonormalization.
      assert(COcc_.nRows == S1.nRows);
      COcc = MakeStackMatrix(COcc_.nRows, nOcc, Mem);
      Assign(COcc, COcc_);
//       SymOrth(COcc, S1, Mem);
      OrthSchmidt(COcc, S1, Mem);
   }

   {
      FStackMatrix
         S2cd(nIb, nIb, &Mem),
         P12(nBf, nIb, &Mem),
         COcc2(nIb, nOcc, &Mem),
         CTil(nIb, nOcc, &Mem),
         STilCd(nIb, nIb, &Mem),
         CTil2BarT(nOcc, nIb, &Mem),
         T4(nBf, nOcc, &Mem);
      MakeIntMatrix(P12, *pOrbBasis, *pMinBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
      MakeIntMatrix(S2cd, *pMinBasis, *pMinBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);

      // COcc2 = mdot(S21, COcc)                # O(N m^2)
      Mxm(COcc2, Transpose(P12), COcc);

      // P12 = la.solve(S1, S12)   # P12 = S1^{-1} S12
      CholeskySolve(P12, S1cd);

      // CTil = la.solve(S2, COcc2)             # O(m^3)
      CalcCholeskyFactors(S2cd);
      Assign(CTil, COcc2);
      CholeskySolve(CTil, S2cd);

      // STil = mdot(COcc2.T, CTil)             # O(m^3)
      Mxm(STilCd, Transpose(COcc2), CTil);

      // CTil2Bar = la.solve(STil, CTil.T).T    # O(m^3)
      CalcCholeskyFactors(STilCd);
      Assign(CTil2BarT, Transpose(CTil));
      CholeskySolve(CTil2BarT, STilCd);

      // T4 = COcc - mdot(P12, CTil2Bar)        # O(N m^2)
      Mxm(T4, P12, Transpose(CTil2BarT), -1.0);
      Add(T4, COcc);

      // CIb = P12 + mdot(T4, COcc2.T)          # O(N m^2)
      Mxm(CIb, T4, Transpose(COcc2));
      Add(CIb, P12);

      if ((Flags & IAO_OrthType) != 0) {
         SymOrth(CIb, S1, Mem);
         if ((Flags & IAO_OrthType) != IAO_OrthSym) {
            // ...
            throw std::runtime_error("sorry, currently only SMH orthogonalization implemented in this version of CtIao.cpp");
         }
      }
   }

   if (COcc.pData != COcc_.pData)
      Mem.Free(COcc.pData);
   IR_SUPPRESS_UNUSED_WARNING(Atoms);
}



// note: CIb is allocated on Mem.
void MakeIaoBasis(FMatrixView &CIb, size_t &nIb, FBasisSetPtr &pMinBasis,
   FAtomSet const &Atoms, FBasisSet const *pOrbBasis, FMatrixView const &COcc, FMemoryStack &Mem, unsigned Flags)
{
   pMinBasis = new FBasisSet(Atoms, BASIS_Guess);
   size_t
      nBf = pOrbBasis->nFn(), // number of 'real basis' functions
      nOcc = COcc.nCols;
   IR_SUPPRESS_UNUSED_WARNING(nOcc);
   nIb = pMinBasis->nFn(); // actual AOs in minimal basis
   CIb = MakeStackMatrix(nBf, nIb, Mem);

   FStackMatrix
      S1(nBf, nBf, &Mem),
      S1cd(nBf, nBf, &Mem);
   MakeIntMatrix(S1, *pOrbBasis, *pOrbBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
   Assign(S1cd,S1);
   CalcCholeskyFactors(S1cd);

   MakeIaoBasisNew(CIb, &*pMinBasis, Atoms, pOrbBasis, COcc, S1, S1cd, Flags, Mem);
//    MakeIaoBasisNew(CIb, nIb, &*pMinBasis, Atoms, pOrbBasis, COcc, S1, S1cd, IAO_OrthSym, Mem);
//    SymOrth(CIb, S1, Mem);
}






// TODO: merge these back into IvDocument.cpp. It has a bunch of very similar
// functions in FOrbital::MakeFullDesc and FOrbital::MakeIaoCharges.
void MakeIaoChargesRaw(double *pOut, double fOcc, double const *pIaoCoeffs, FBasisSet *pMinBasis, FAtomSet *pAtoms)
{
   size_t
      iShOf = 0;
   memset(pOut, 0, sizeof(pOut[0]) * pAtoms->size());
   for (size_t iSh = 0; iSh < pMinBasis->Shells.size(); ++ iSh) {
      ct::FBasisShell
         &Sh = pMinBasis->Shells[iSh];
      uint
         nShFn = Sh.nFn();
      double
         fChg = 0;
      for (uint iFn = 0; iFn < nShFn; ++ iFn)
         fChg += sqr(pIaoCoeffs[iFn + iShOf]);
      assert(size_t(Sh.iCenter) < pAtoms->size());
      pOut[Sh.iCenter] += fOcc * fChg;
      iShOf += nShFn;
   }
   assert(iShOf == pMinBasis->nFn());
}



} // namespace ct
