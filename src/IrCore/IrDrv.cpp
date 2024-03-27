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

#include <algorithm> // for std::swap
#include <cmath>


#include "IrDrvAux.h"
#include "IrAmrr.h"
#include "IrFactorials.h"

// #define TEST_3i2c_CODE


namespace ir {

// references:
//   [1]: PCCP 8 3072 (2006), doi: 10.1039/b605188j
//   [2]: PCCP 6 5119 (2004), doi: 10.1039/b413539c
//   [3]: for two-index two-center integrals and integral kernels:
//        Mieke Peels, Gerald Knizia - "Fast evaluation of two-center integrals over
//        Gaussian charge distributions and Gaussian orbitals with general interaction
//        kernels", J. Chem. Theory Comput. 2020
//        https://doi.org/10.1021/acs.jctc.9b01296
//   [4]: for general comments on numerical stability, and the PmA = P-A = beta/zeta * (B-A) stability matter
//        Knizia, Li, Simon, Werner, J. Chem. Theory Comput. 2011, 7, 8, 2387-2398, https://doi.org/10.1021/ct200239p



static size_t *MakeFnOffsets(FRawShell const *pCs, size_t nC, FMemoryStack &Mem)
{
   size_t *piFnC;
   Mem.Alloc(piFnC, nC+1);
   Mem.Align(16);
   piFnC[0] = 0;
   for (size_t iC = 0; iC < nC; ++ iC)
      piFnC[iC+1] = piFnC[iC] + pCs[iC].nFn();
   return piFnC;
}

static double
   dbl = 0.; // dummy.


#if 0
// this is the base version of the 2e3c integral routine. It is currently
// replaced by the generalzied version below, which separates the contraction
// work from the primitive work to support derivative integrals and inline-
// contracting drivers. It is retained here in order to illustrate what the
// generalized version does.

void EvalInt2e3c_(double *pOut, size_t *Strides,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   void
      *pBaseOfMemory = Mem.Alloc(0);
   size_t
      StrideA = Strides[0], StrideB = Strides[1], StrideC = Strides[2];
   if (pA->l < pB->l) { // <- OsrrC only implemented for la >= lb.
      std::swap(pA, pB);
      std::swap(StrideA, StrideB);
   }

   // count number of C functions and find largest lc for memory purposes.
   size_t
      *piFnC = MakeFnOffsets(pCs, nC, Mem),
      nFnC_Total = piFnC[nC];
   uint
      lc_Max = 0;
   for (size_t iC = 0; iC < nC; ++ iC)
      lc_Max = std::max(lc_Max, (uint)pCs[iC].l);
   // allocate intermediates
   size_t
      nCartX_AB = nCartX(pA->l + pB->l),
      nCartX_Am1 = nCartX(pA->l-1),
      nCartX_B = nCartX(pB->l),
      nCartX_ABmA = nCartX_AB - nCartX_Am1,
      nShA_CartXB = pA->nSh() * nCartX_B;
   // intermediates for primitive integrals
   double
      // Kernel derivatives (-d/dT)^m of (00|00) integral
      *pGm = Mem.AllocN(pA->l + pB->l + lc_Max + 1, dbl),
      // (a0|0) intermediate
      *p_A00 = Mem.AllocN(nCartX_AB, dbl),
      *p_A0C_sh_mem = Mem.AllocN(nCartX_ABmA * (2*lc_Max+1), dbl),
      *pMemOsrrB = Mem.AllocN(nCartX_AB * nCartX(lc_Max), dbl);
   // intermediates for contractions
   double
      // intermediates (a0|c) with AB primitives and C contracted, a = la..lab
      *p_A0C_ppc = Mem.AllocN(nCartX_ABmA * nFnC_Total, dbl),
      // intermediates (a0|c) with A,C contracted, a = la..lab.
      *p_A0C_cpc = Mem.AllocN(nCartX_ABmA * nFnC_Total * pA->nCo, dbl),
      // intermediates (a0|c) with A,B,C all contracted, a = la..lab.
      *p_A0C_ccc = Mem.ClearAllocN(nCartX_ABmA * nFnC_Total * pA->nCo * pB->nCo, dbl),
      // intermediates (xa|c) with A,B,C contracted and (xa| = nCartX(lb) x (2*la+1)
      *p_xAC_ccc = Mem.AllocN(nShA_CartXB * nFnC_Total * pA->nCo * pB->nCo, dbl);

//    printf("size on p_A0C_ccc = %i  expected: %i\n", p_xAC_ccc-p_A0C_ccc, nCartX_ABmA * nFnC_Total * pA->nCo * pB->nCo);
   FVec3
      vAmB = FVec3(pA->vCen) - FVec3(pB->vCen);
   double
//       fRangeKernel = sqr(pKernel->MaxRange()),
      fDistSqAB = LengthSq(vAmB);

   for (uint iExpB = 0; iExpB < pB->nExp; ++ iExpB)
   {
      memset(p_A0C_cpc, 0, nCartX_ABmA * nFnC_Total * pA->nCo * sizeof(*p_A0C_cpc));
      for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
      {
         // skip if Dist(A,B) < Range(A) + Range(B)
         if (!IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, fDistSqAB))
            continue;

         FGaussProduct
            OvAB(pA->vCen, pA->pExp[iExpA], pB->vCen, pB->pExp[iExpB]);
            // ^- P == OvAB.vCen
         double
            Sab = std::exp(-OvAB.Exp * fDistSqAB); // [1] (6)

         memset(p_A0C_ppc, 0, nCartX_ABmA * nFnC_Total * sizeof(*p_A0C_ppc));
         for (size_t iC = 0; iC < nC; ++ iC) {
            FRawShell const
               *pC = &pCs[iC];
            uint
               TotalL = pA->l + pB->l + pC->l;
            for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
            {
               FGaussProduct
                  OvPC(OvAB.vCen, OvAB.Zeta, pC->vCen, pC->pExp[iExpC]);
               double
                  *PmC = OvPC.vAmB,
                  rho = OvPC.Exp, // [1] (3)
                  T = rho * OvPC.DistSq,
                  Factor = pow15(M_PI * OvPC.InvZeta) * Sab * Prefactor; // [1] (7)

               // make I[m] = (00|0)^m, m = 0..TotalL (inclusive)
               pKernel->EvalGm(pGm, rho, T, TotalL, Factor);

               // make (a0|0)^m for a = 0..lab with lab = la+lb.
               OsrrA(p_A00, pGm + pC->l, (pA->l + pB->l),
                  OvAB.vPmA[0], OvAB.vPmA[1], OvAB.vPmA[2],
                  PmC[0], PmC[1], PmC[2], rho, OvAB.InvZeta);

               // make (a0|c) for a = la..lab, c = 0..lc.
               double
                  *p_A0C_sh;
               if (pC->l == 0) {
                  p_A0C_sh = p_A00 + nCartX_Am1;
               } else {
                  p_A0C_sh = p_A0C_sh_mem;
                  OsrrB_3c_shc(p_A0C_sh, p_A00, pMemOsrrB, pA->l, (pA->l + pB->l), pC->l,
                     PmC[0], PmC[1], PmC[2], OvPC.InvZeta, rho/pC->pExp[iExpC]);
               }

               // (a0|c) with solid harmonic c is ready now. Just need to add it to
               // its contractions.
               Contract1(&p_A0C_ppc[nCartX_ABmA * piFnC[iC]], p_A0C_sh,
                  nCartX_ABmA*pC->nSh(), pC, iExpC);
            } // c exponents
         } // c shells
         // p_A0C_ppc should be done now. Contract A and B.
         Contract1(p_A0C_cpc, p_A0C_ppc, nCartX_ABmA * nFnC_Total, pA, iExpA);
      } // a exponents
      Contract1(p_A0C_ccc, p_A0C_cpc, (nCartX_ABmA * nFnC_Total * pA->nCo), pB, iExpB);
   } // b exponents

   // transform A to solid harmonics by factoring nCartX(lab) into nCartX(lb) x Slm(A).
   ShTrA_XY(p_xAC_ccc, p_A0C_ccc, pA->l, (pA->l + pB->l), nFnC_Total * pA->nCo * pB->nCo);

   // we now have nCartX(lb) x nShA x nFnC_Total x nCoA x nCoB at p_xAC_ccc.
   // we still need to move the angular momentum from a to b and to write the
   // output integrals to their final destination.
   for (uint iCoB = 0; iCoB < pB->nCo; ++ iCoB)
      for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
         for (uint iFnC = 0; iFnC < nFnC_Total; ++ iFnC) {
            uint
               iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
            OsrrC(
               &pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC], StrideA, StrideB,
               &p_xAC_ccc[nShA_CartXB * (iFnC + (nFnC_Total * (iCoA + pA->nCo * iCoB)))],
               vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh() );
         };
      }

   Mem.Free(pBaseOfMemory);
}

#endif // 0





// this auxiliary structure covers common functionality of 3-center shell
// drivers for direct 2e3c integrals, derivative integrals, and inline-
// contracting variants of the same.
//
// TODO:
// - implement the early c-contraction as described in FTCS paper
// - consider using the FPrimitivePairData auxiliary structure from the 4c driver routine
struct FInt3cShellFactory
{
   FInt3cShellFactory(FRawShell const *pA_, FRawShell const *pB_, FRawShell const *pCs_, size_t nC_,
                      double Prefactor_, FIntegralKernel const *pKernel_, FMemoryStack &Mem_)
      : pA(pA_), pB(pB_), pCs(pCs_), nC(nC_), Prefactor(Prefactor_), pKernel(pKernel_), Mem(Mem_)
   {
      SwapAB = pA->l < pB->l;
      if (SwapAB) // <- OsrrC only implemented for la >= lb.
         std::swap(pA, pB);
   }
protected:
   IR_FORCE_INLINE void EvalDrv(uint lab_Min, uint lab_Max, bool ShTrC);
   // ^- note @IR_FORCE_INLINE: without inlining the 3c main driver function (a
   // prerequisite for inlining the virtual functions) there can be substantial
   // call overhead. Not disastrous yet, but quite noticeable. As the driver
   // function is rather large, compilers might be inclined to ignore our
   // ``inline'' directive. Try to tell them that we really mean it.

   FRawShell const
      *pA, *pB, *pCs;
   size_t
      nC;
   double
      Prefactor;
   FIntegralKernel const
      *pKernel;
   FMemoryStack
      &Mem;

   virtual void BeginContractions() = 0;
   virtual void BeginContraction1(uint iExpA) = 0;
   virtual void BeginContraction2(uint iExpA, uint iExpB) = 0;
   virtual void ContractPrimitiveA0C(double *p_A0C, uint iExpA, uint iExpB, uint iExpC, size_t iC, size_t iFnC, FRawShell const *pC) = 0;
   virtual void EndContraction2(uint iExpA, uint iExpB) = 0;
   virtual void EndContraction1(uint iExpA) = 0;
   virtual void EndContractions() = 0;
private:
   void
      *pBaseOfMemory;
   double
      *pGm, *p_A00, *p_A0C_mem, *pMemOsrrB;
   inline void InitDrv(uint lab_Min, uint lab_Max, bool ShTrC);
protected:
   bool
      SwapAB; // set if A and B were swapped due to la < lb.
   FVec3
      vAmB; // pA->vCen - pB->vCen
   size_t
      *piFnC,
      nFnC_Total;
   size_t
      nCartX_ABmx, nCartX_ABmn, nCartX_ABac;
};


void FInt3cShellFactory::InitDrv(uint lab_Min, uint lab_Max, bool ShTrC)
{
   pBaseOfMemory = Mem.Alloc(0);
   Mem.Align(16);

   vAmB = FVec3(pA->vCen) - FVec3(pB->vCen);

   // count number of C functions and find largest lc for memory purposes.
   piFnC = MakeFnOffsets(pCs, nC, Mem);
   nFnC_Total = piFnC[nC];
   uint
      lc_Max = 0;
   for (size_t iC = 0; iC < nC; ++ iC)
      lc_Max = std::max(lc_Max, (uint)pCs[iC].l);
   size_t
      nCompC_Max = ShTrC? (2*lc_Max+1) : nCartY(lc_Max);

   // allocate intermediates for primitive integrals
   nCartX_ABmx = nCartX(lab_Max); // number of (a0| components for a in 0..lab_Max
   nCartX_ABmn = nCartX(lab_Min-1); // number of (a0| components for a in 0..(lab_Min-1)
   nCartX_ABac = ALIGN_CO_SIZE_DBL(nCartX_ABmx - nCartX_ABmn); // actual number of components: difference

   // Kernel derivatives (-d/dT)^m of (00|00) integral
   pGm = Mem.AllocN(lab_Max + lc_Max + 1, dbl);

   // (a0|0) intermediate
   p_A00 = Mem.AllocN(nCartX_ABmx, dbl);
   ALIGN_CO_MEM(Mem)
   p_A0C_mem = Mem.AllocN(nCartX_ABac * nCompC_Max, dbl);
   pMemOsrrB = Mem.AllocN(nCartX_ABmx * nCartX(lc_Max), dbl);

   IR_SUPPRESS_UNUSED_WARNING(ir::IsWithinRange);
   IR_SUPPRESS_UNUSED_WARNING(ir::IsContractionWithinRange);
   IR_SUPPRESS_UNUSED_WARNING(static_cast<bool (*)(FRawShell const *, uint, FRawShell const *, uint, double)>(ir::IsPrimitiveWithinRange));
   IR_SUPPRESS_UNUSED_WARNING(static_cast<bool (*)(FRawShell const *, uint, FRawShell const *, uint)>(ir::IsPrimitiveWithinRange));
}


void FInt3cShellFactory::EvalDrv(uint lab_Min, uint lab_Max, bool ShTrC)
{
   IR_TIME_SECTION(TID_EvalInt2eNc);
   InitDrv(lab_Min, lab_Max, ShTrC);

   double
#ifdef IR_ENABLE_KERNEL_RANGE_SCREENING
      fRangeKernelSq = sqr(pKernel->MaxRange()),
#endif // IR_ENABLE_KERNEL_RANGE_SCREENING
      fDistSqAB = LengthSq(vAmB);

   // choose the function to use for the OsrrB recurrence relation. We have
   // a normal one (producing (a0|c) with nCartY(lc) cartesian cs) and an inline-
   // ShTr'ing one (producing (a0|c) with (2lc+1) solid harmonic cs)
   typedef
      void (*FOsrrBFn)(double *IR_RP pOut, double const *IR_RP pIn, double *IR_RP pMem, int la, unsigned lab, unsigned lc, double fPmQx, double fPmQy, double fPmQz, double InvEtaABC, double riz);
   FOsrrBFn
      OsrrB_3c = ShTrC? OsrrB_3c_shc : OsrrB_3c_cac;

   IR_RESUME_CLOCK(TID_PrimLoop);
   BeginContractions(); // <- in timing print it's treated as part of the prim loop, although it really isn't
   for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
   {
      BeginContraction1(iExpA);
      for (uint iExpB = 0; iExpB < pB->nExp; ++ iExpB)
      {
         // skip if Dist(A,B) < Range(A) + Range(B)
         if (!IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, fDistSqAB))
            continue;

         FGaussProduct
            OvAB(pA->vCen, pA->pExp[iExpA], pB->vCen, pB->pExp[iExpB]);
            // ^- P == OvAB.vCen
         double
            Sab = OvAB.Sab();

         BeginContraction2(iExpA, iExpB);
         for (size_t iC = 0; iC < nC; ++ iC) {
            FRawShell const
               *pC = &pCs[iC];
#ifdef IR_ENABLE_KERNEL_RANGE_SCREENING
            double
               fDistSqACK = DistSq3(pA->vCen, pC->vCen) - fRangeKernelSq,
               fDistSqBCK = DistSq3(pB->vCen, pC->vCen) - fRangeKernelSq;
#endif // IR_ENABLE_KERNEL_RANGE_SCREENING
            uint
               TotalL = lab_Max + pC->l;
            for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
            {
#ifdef IR_ENABLE_KERNEL_RANGE_SCREENING
               if (!IsPrimitiveWithinRange(pA, iExpA, pC, iExpC, fDistSqACK) ||
                   !IsPrimitiveWithinRange(pB, iExpB, pC, iExpC, fDistSqBCK))
                  continue;
#endif // IR_ENABLE_KERNEL_RANGE_SCREENING
               IR_RESUME_CLOCK(TID_PrimData);
               FGaussProduct
                  OvPC(OvAB.vCen, OvAB.Zeta, pC->vCen, pC->pExp[iExpC]);
               double
                  *PmC = OvPC.vAmB,
                  rho = OvPC.Exp, // [1] (3)
                  T = rho * OvPC.DistSq,
                  Factor = pow15(M_PI * OvPC.InvZeta) * Sab * Prefactor; // [1] (7)
               IR_PAUSE_CLOCK(TID_PrimData);

               {  IR_TIME_SECTION(TID_EvalGm);
                  // make I[m] = (00|0)^m, m = 0..TotalL (inclusive)
                  pKernel->EvalGm(pGm, rho, T, TotalL, Factor);
               }

               {  IR_TIME_SECTION(TID_OsrrA);
                  // make (a0|0)^m for a = 0..lab with lab = la+lb.
                  OsrrA(p_A00, pGm + pC->l, lab_Max,
                     OvAB.vPmA[0], OvAB.vPmA[1], OvAB.vPmA[2],
                     PmC[0], PmC[1], PmC[2], rho, OvAB.InvZeta);
               }

               // make (a0|c) for a = la..lab, c = 0..lc.
               double
                  *p_A0C;
               {  IR_TIME_SECTION(TID_OsrrB);
// #ifndef ALIGN_CO
                  if (pC->l == 0) {
                     p_A0C = p_A00 + nCartX_ABmn;
                  } else {
// #endif
                     p_A0C = p_A0C_mem;
                     OsrrB_3c(p_A0C, p_A00, pMemOsrrB, lab_Min, lab_Max, pC->l,
                        PmC[0], PmC[1], PmC[2], OvPC.InvZeta, rho/pC->pExp[iExpC]);
                  }
               }

               // (a0|c) is ready now. Just need use it in its contractions.
               ContractPrimitiveA0C(p_A0C, iExpA, iExpB, iExpC, iC, piFnC[iC], pC);
            } // c exponents
         } // c shells
         EndContraction2(iExpA, iExpB);
      } // b exponents
      EndContraction1(iExpA);
   } // a exponents
   IR_PAUSE_CLOCK(TID_PrimLoop);
   EndContractions();

   Mem.Free(pBaseOfMemory);
}


// calculates normal integrals (ab|c) and stores them at
//    pOut[ia * StrideA + ib * StrideB + ic * StrideC].
struct FInt2e3cShellFactory : public FInt3cShellFactory
{
   FInt2e3cShellFactory(
      double *pOut_, size_t StrideA_, size_t StrideB_, size_t StrideC_,
      FRawShell const *pA_, FRawShell const *pB_, FRawShell const *pCs_, size_t nC_, double Prefactor_, FIntegralKernel const *pKernel_, FMemoryStack &Mem_)
      : FInt3cShellFactory(pA_, pB_, pCs_, nC_, Prefactor_, pKernel_, Mem_),
        pOut(pOut_), StrideA(StrideA_), StrideB(StrideB_), StrideC(StrideC_)
   {}

   void Eval() {
      EvalDrv(pA->l, pA->l + pB->l, true);
   }
protected:
   double
      *pOut;
   size_t
      StrideA, StrideB, StrideC, nFnC_Out;

   // intermediates for contractions
   double
      // intermediates (a0|c) with AB primitives and C contracted, a = la..lab
      *p_A0C_ppc,
      // intermediates (a0|c) with B,C contracted, a = la..lab.
      *p_A0C_cpc,
      // intermediates (a0|c) with A,B,C all contracted, a = la..lab.
      *p_A0C_ccc;

   IR_FORCE_INLINE void SetupCoBuffers(size_t nFnC_) {
      assert(nCartX_ABac == (nCartX(pA->l + pB->l) - nCartX((signed)pA->l-1)));
      nFnC_Out = nFnC_;

      if (SwapAB)
         std::swap(StrideA, StrideB);

      // allocate intermediates for partially contracted A0C integrals.
      p_A0C_ppc = Mem.AllocN(nCartX_ABac * nFnC_Out, dbl);
      p_A0C_cpc = Mem.AllocN(nCartX_ABac * nFnC_Out * pB->nCo, dbl);
      p_A0C_ccc = Mem.AllocN(nCartX_ABac * nFnC_Out * pB->nCo * pA->nCo, dbl);
      ClearContractBuffer(p_A0C_ccc, nCartX_ABac * nFnC_Out * pB->nCo * pA->nCo);
   }

   IR_FORCE_INLINE void BeginContractions() { // override
      SetupCoBuffers(nFnC_Total);
   }

   IR_FORCE_INLINE void BeginContraction1(uint /*iExpB*/) { // override
      ClearContractBuffer(p_A0C_cpc, nCartX_ABac * nFnC_Out * pB->nCo);
   }

   IR_FORCE_INLINE void BeginContraction2(uint /*iExpA*/, uint /*iExpB*/){ // override
      ClearContractBuffer(p_A0C_ppc, nCartX_ABac * nFnC_Out);
   }

   IR_FORCE_INLINE void ContractPrimitiveA0C(double *p_A0C, uint /*iExpA*/, uint /*iExpB*/, uint iExpC, size_t /*iC*/, size_t iFnC, FRawShell const *pC) { // override
      Contract1(&p_A0C_ppc[nCartX_ABac * iFnC], p_A0C,
         nCartX_ABac*pC->nSh(), pC, iExpC);
   }
   IR_FORCE_INLINE void EndContraction2(uint /*iExpA*/, uint iExpB) { // override
      Contract1(p_A0C_cpc, p_A0C_ppc, nCartX_ABac * nFnC_Out, pB, iExpB);
   }

   IR_FORCE_INLINE void EndContraction1(uint iExpA) { // override
      Contract1(p_A0C_ccc, p_A0C_cpc, (nCartX_ABac * nFnC_Out * pB->nCo), pA, iExpA);
   }

//    IR_FORCE_INLINE void EndContractions() { // override
//       IR_TIME_SECTION(TID_Finalize);
//       // transform A to solid harmonics by factoring nCartX(lab) into nCartX(lb) x Slm(A).
//       size_t
//          nShA_CartXB = (2*pA->l+1) * nCartX(pB->l);
//       double
//          // intermediates (xa|c) with A,B,C contracted and (xa| = nCartX(lb) x (2*la+1)
//          *p_xAC_ccc = Mem.AllocN(nShA_CartXB * nFnC_Out * pB->nCo * pA->nCo, dbl);
//       {  IR_TIME_SECTION(TID_ShTr1);
//          ShTrA_XY(p_xAC_ccc, p_A0C_ccc, pA->l, (pA->l + pB->l), nFnC_Out * pB->nCo * pA->nCo);
//       }
//
//       {  IR_TIME_SECTION(TID_OsrrC1);
//          // we now have nCartX(lb) x nShA x nFnC_Out x nCoA x nCoB at p_xAC_ccc.
//          // we still need to move the angular momentum from a to b and to write the
//          // output integrals to their final destination.
//          for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB) {
//             for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
//                for (size_t iFnC = 0; iFnC < nFnC_Out; ++ iFnC) {
//                   size_t
//                      iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
//                   OsrrC(
//                      &pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC], StrideA, StrideB,
//                      &p_xAC_ccc[nShA_CartXB * (iFnC + (nFnC_Out * (iCoB + pB->nCo * iCoA)))],
//                      vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh() );
//                }
//             }
//          }
//       }
//    }
   IR_FORCE_INLINE void EndContractions() { // override
      IR_TIME_SECTION(TID_Finalize);
      // we now have nCartX(la+lb,lb) nFnC_Out x nCoA x nCoB at p_A0C_ccc. we
      // still need to transform A to spherical harmonics and move the
      // angular momentum from a to b and to write the output integrals to
      // their final destination.
      for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB) {
         for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
            for (size_t iFnC = 0; iFnC < nFnC_Out; ++ iFnC) {
               size_t
                  iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
               OsrrC_sha(
                  &pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC], StrideA, StrideB,
                  &p_A0C_ccc[nCartX_ABac * (iFnC + (nFnC_Out * (iCoB + pB->nCo * iCoA)))],
                  pA->l,
                  vAmB[0], vAmB[1], vAmB[2], pB->l);
            }
         }
      }
   }
};





void EvalInt2e3c(double *pOut, size_t *Strides,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
//    EvalInt2e3c_(pOut, Strides, pA, pB, pCs, nC, Prefactor, pKernel, Mem);
   FInt2e3cShellFactory(pOut, Strides[0], Strides[1], Strides[2],
      pA, pB, pCs, nC, Prefactor, pKernel, Mem).Eval();
}

//#ifdef INCLUDE_OPTIONALS
// transform a set of vectors given in terms of contracted solid harmonic functions
// to cartesian primitive functions. Output pCaPrimOfs[iC] gives the offset of
// pC's functions in terms of cartesian primitives.
static void TrafoVecToCaPrim(double *&pVecCaPrim_, size_t &StrideVecCaPrim, size_t *&pCaPrimOfs,
   double const *pVecShCo_, size_t StrideVec, size_t nVec,
   FRawShell const *pCs, size_t nCs, FMemoryStack &Mem)
{
   assert(nVec == 1);
   Mem.Alloc(pCaPrimOfs, nCs+1);
   Mem.Align(16);
   // make offsets of C functions in terms of cartesian primitives.
   pCaPrimOfs[0] = 0;
   for (size_t iC = 0; iC < nCs; ++ iC)
      pCaPrimOfs[iC+1] = pCaPrimOfs[iC] + nCartY(pCs[iC].l) * pCs[iC].nExp;
   StrideVecCaPrim = pCaPrimOfs[nCs];

   Mem.ClearAlloc(pVecCaPrim_, StrideVecCaPrim * nVec);
   double
      *pVecCaPrim = pVecCaPrim_;
   // transform each vector individually.
   for (size_t iVec = 0; iVec < nVec; ++ iVec){
      double const
         *pVecShCo = pVecShCo_ + iVec * StrideVec;
      for (size_t iC = 0; iC < nCs; ++ iC) {
         FRawShell const
            *pC = &pCs[iC];
         // transform to cartesians: (2lc+1) x nCo -> nCartY(lc) x nCo
         double
            *pCaCo;
         size_t
            nCa = nCartY(pC->l);
         Mem.Alloc(pCaCo, nCa * pC->nCo);
         CaTrC(pCaCo, pVecShCo, pC->nCo, pC->l);

         // uncontract
         for (size_t iCo = 0; iCo < pC->nCo; ++ iCo)
            for (size_t iExp = 0; iExp < pC->nExp; ++ iExp)
               for (size_t iCa = 0; iCa < nCa; ++ iCa)
                  pVecCaPrim[iCa + nCa*iExp] += pCaCo[iCa + nCa*iCo] * pC->fCo(unsigned(iExp),unsigned(iCo));

         Mem.Free(pCaCo);

         pVecShCo += pC->nSh() * pC->nCo;
         pVecCaPrim += nCa * pC->nExp;
      }
   }
}

// calculates the contraction j[ab,vc] = \sum_c (ab|c) gamma[c,vc] and stores the result at
//     pOut[ia * StrideA + ib * StrideB].
// gamma[c,vc] is taken from pIn[ic + vc * StrideC].
struct FInt2e3cVcaFactory : public FInt2e3cShellFactory
{
   FInt2e3cVcaFactory(
      double *pOut_, size_t StrideA_, size_t StrideB_, size_t StrideC_, double const *pIn_, size_t StrideIn_, size_t nVecC_,
      FRawShell const *pA_, FRawShell const *pB_, FRawShell const *pCs_, size_t nC_, double Prefactor_, FIntegralKernel const *pKernel_, FMemoryStack &Mem_)
      : FInt2e3cShellFactory(pOut_, StrideA_, StrideB_, StrideC_,
                             pA_, pB_, pCs_, nC_, Prefactor_, pKernel_, Mem_),
        pIn(pIn_), StrideIn(StrideIn_), nVecC(nVecC_)
   {}

   void Eval() {
      EvalDrv(pA->l, pA->l + pB->l, false);
   }
protected:
   // inputs
   double const
      *pIn;
   size_t
      StrideIn, nVecC;

   // input vector pIn transformed to cartesian primitives & assorted strides.
   double
      *pInCaPrim;
   size_t
      StrideInCaPrim, *pInCaPrimOfs;

   IR_FORCE_INLINE void BeginContractions() { // override
      // transform input vectors to cartesian primitives.
      TrafoVecToCaPrim(pInCaPrim, StrideInCaPrim, pInCaPrimOfs, // <- outputs
         pIn, StrideIn, nVecC, pCs, nC, Mem);

      SetupCoBuffers(nVecC);
   }

   IR_FORCE_INLINE void ContractPrimitiveA0C(double *p_A0C, uint iExpA, uint iExpB, uint iExpC, size_t iC, size_t iFnC, FRawShell const *pC)
   { // override
      // contract primitive (a0|c) with transformed input vectors over c.
      size_t
         nCompC = nCartY(pC->l);
      double
         *pInC = &pInCaPrim[pInCaPrimOfs[iC] + nCompC * iExpC];
      assert(nVecC == nFnC_Out);
      for (size_t iVecC = 0; iVecC < nVecC; ++ iVecC) {
         double
            *c = &p_A0C_ppc[nCartX_ABac * iVecC],
            *v = &pInC[StrideInCaPrim * iVecC];
         // TODO: check if using real mxv makes a difference here. atm I am
         // still reluctant to introduce an external dependency for this.
         for (size_t ic = 0; ic < nCompC; ++ ic)
            for (size_t iab = 0; iab < nCartX_ABac; ++ iab)
               c[iab] += v[ic] * p_A0C[iab + nCartX_ABac * ic];
      }
      IR_SUPPRESS_UNUSED_WARNING3(iExpA, iExpB, iFnC);
   }
};

void EvalInt2e3c_ContractC(double *pOut, size_t *Strides, double const *pIn, size_t StrideIn, size_t nVecIn,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
//    EvalInt2e3c_(pOut, Strides, pA, pB, pCs, nC, Prefactor, pKernel, Mem);
   FInt2e3cVcaFactory(pOut, Strides[0], Strides[1], Strides[2], pIn, StrideIn, nVecIn,
      pA, pB, pCs, nC, Prefactor, pKernel, Mem).Eval();
}


// transfer product density (ab| centered at A and B to cartesian density of ([a+b]0| centered at A (i.e., do the inverse of OsrrC.)
static void OsrrRevC(double *pOut, double const *pDenAB, size_t sa, size_t sb, unsigned la, unsigned lb, double const AmB[3], FMemoryStack &Mem)
{
   // pDenAB[ia*sa + ib*sb]: density for (2*la+1)x(2*lb+1) solid harmonic elements, centered at A and B
   // pOut: array of nCartX(la+lb) elements, containing the density expanded to cartesians centered at A.
   size_t
      nCartX_A = nCartX(la),
      iBaseCompA = nCartX(la-1),
      nCartX_B = nCartX(lb);
   size_t
      nShA = 2*la + 1;
   double
      *pTmp0,
      *pTmp1;
   // expand B from solid harmonics centerd at B into cartesians centered at A.
   Mem.Alloc(pTmp0, nCartX_B * nShA);
   for (size_t iShA = 0; iShA < nShA; ++ iShA)
      OsrrRx(pTmp0 + nCartX_B * iShA, &pDenAB[iShA * sa], sb, AmB[0], AmB[1], AmB[2], lb);

   // expand A from solid harmonics to cartesians (A is already centered at A)
   Mem.Alloc(pTmp1, nCartX_A * nCartX_B);
   for (size_t iCartB = 0; iCartB < nCartX_B; ++ iCartB) {
      double *p = &pTmp1[nCartX_A * iCartB];
      if (iBaseCompA != 0)
         memset(p, 0, sizeof(*p) * iBaseCompA);
      CaTrA(&pTmp1[iBaseCompA + nCartX_A * iCartB], &pTmp0[iCartB], nCartX_B, la);
   }

   // now have nCartX_A * nCartX_B array of density expanded
   // in terms of (r-A)^n. Still need to add up the cartesian powers.
   // note that we INCREMENT the output values!
   for (size_t iCartB = 0; iCartB != nCartX_B; ++ iCartB) {
      for (size_t iCartA = iBaseCompA; iCartA != nCartX_A; ++ iCartA)
         pOut[iv2x[ix2v[iCartB] + ix2v[iCartA]]] += pTmp1[iCartA + iCartB * nCartX_A];
   }

   Mem.Free(pTmp0);
}


// Transfer contracted product density (ab| centered at A and B to primitive
// cartesian density of Lab = la+lb centered at A.
//
// Notes:
//    - A is where the larger angular momentum is.
//    - Output is allocated on pCoFactorsAB.
//    - output data is nCartX(l_a+l_b, lab_Min) * nExpA * nExpB.
//    - by default this makes output components from lab_Min=la to lab=la+lb.
//      One can add extra lower lab_Min components by supplying lab_Min_extra > 0.
void TrafoDensityToCaPrim(double *&pCoFactorsAB, FRawShell const *pA, FRawShell const *pB, unsigned lab_Min_extra, double const *pDenAB, size_t sa, size_t sb, FMemoryStack &Mem)
{
   if (pB->l > pA->l) {
      std::swap(pA, pB);
      std::swap(sa, sb);
   }
   unsigned
      la = pA->l,
      lb = pB->l,
      lab = la+lb,
      lab_Min = pA->l - lab_Min_extra;
   if (lab_Min_extra > pA->l)
      lab_Min = 0;
   assert(lab_Min <= lab);
   size_t
      nShA = 2*la+1,
      nShB = 2*lb+1,
      nCartX_AB = nCartX(lab),
      nCartX2_AB = nCartX(lab, lab_Min);
   size_t
      nExpA = pA->nExp, nCoA = pA->nCo,
      nExpB = pB->nExp, nCoB = pB->nCo;
   double
      *pCartDenA;
   Mem.ClearAlloc(pCoFactorsAB, nCartX2_AB * nExpA * nExpB);

   // first transfer contracted (ab| density to contracted ([a+b]0| density
   Mem.ClearAlloc(pCartDenA, nCartX_AB * pA->nCo * pB->nCo);
   double
      AmB[3];
   SubVec3(&AmB[0], pA->vCen, pB->vCen);

   for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB)
      for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA)
         OsrrRevC(&pCartDenA[(iCoA + iCoB * pA->nCo) * nCartX_AB],
                  &pDenAB[iCoA*nShA*sa + iCoB*nShB*sb], sa, sb, pA->l, pB->l, AmB, Mem);

   // now uncontract from CartX(lab) x nCoA x nCoB to CartX(lab) x nExpA x nExpB
   for (size_t iCoB = 0; iCoB < nCoB; ++ iCoB)
      for (size_t iCoA = 0; iCoA < nCoA; ++ iCoA)
         for (size_t iExpB = 0; iExpB < nExpB; ++ iExpB)
            for (size_t iExpA = 0; iExpA < nExpA; ++ iExpA) {
               double
                  fCoAB = pA->pCo[iExpA + nExpA * iCoA] * pB->pCo[iExpB + nExpB * iCoB];
               if (fCoAB == 0.)
                  continue;
               double
                  *pOut0 = &pCoFactorsAB[(iExpA + nExpA * iExpB) * nCartX2_AB],
                  *pIn0 = &pCartDenA[(iCoA + nCoA * iCoB) * nCartX_AB + (nCartX_AB - nCartX2_AB)];
               Add2(pOut0, pIn0, fCoAB, nCartX2_AB);
            }
   Mem.Free(pCartDenA);

   // ^- TODO:
   // - measure this out and check if this is expensive, and if yes,
   //   possibly split it into two steps, which has lower scaling. For segmented
   //   sets I guess it would be fine. Also check loop order, maybe
   // - Additionally, we might want to consider if there is a way to do the
   //   OsrrRevC also with only nCartX2(lab,la) components, instead of all nCartX(lab)
   //   ones (but I guess in most cases it would not matter much, provided the _ContractAB
   //   are called with more than one or two C target shells.)
   //   Here we currently do the compromise of just holding the actual shell
   //   driver to the minimum la components (note difference nCartX_AB - nCartX2_AB
   //   below)
   //
   // UPDATE:
   // - had a try with preg.xyz and TZVP... took less than 1% of
   //   EvalInt2e3c_ContractAB (with 12 threads on mirage).
   // - so I guess we can just leave it as is for the time being.
}


static size_t *MakeCartesianFnOffsets(FRawShell const *pCs, size_t nC, FMemoryStack &Mem)
{
   size_t *piFnC;
   Mem.Alloc(piFnC, nC+1);
   Mem.Align(16);
   piFnC[0] = 0;
   for (size_t iC = 0; iC < nC; ++ iC)
      piFnC[iC+1] = piFnC[iC] + nCartY(pCs[iC].l) * pCs[iC].nCo;
//       piFnC[iC+1] = piFnC[iC] + (2*pCs[iC].l+1) * pCs[iC].nCo;
   return piFnC;
}


// calculates the contraction j[ab,vc] = \sum_ab (ab|c) gamma[ab,vc] and stores the result at
//     pOut[ic].
// gamma[ab] must have been made by TrafoDensityToCaPrim. Should be extended to supporting multiple
// sets of input densities gamma[ab,iDenSet] later. But not yet there.
struct FInt2e3cAcsFactory : public FInt3cShellFactory
{
   FInt2e3cAcsFactory(
      double *pOut_, size_t StrideC_DenSet_, // <- stride for output density set
      FRawShell const *pA_, FRawShell const *pB_, double const *pCoFactorsAB_, size_t sDenSet_, size_t nDenSets_,
      FRawShell const *pCs_, size_t nC_, double Prefactor_, FIntegralKernel const *pKernel_, FMemoryStack &Mem_)
      : FInt3cShellFactory(pA_, pB_, pCs_, nC_, Prefactor_, pKernel_, Mem_),
        pOut(pOut_), pCoFactorsAB(pCoFactorsAB_), StrideInDenSet(sDenSet_), nDenSets(nDenSets_), StrideC_DenSet(StrideC_DenSet_)
   {}

   void Eval() {
//       EvalDrv(pA->l, pA->l + pB->l, true);
      EvalDrv(pA->l, pA->l + pB->l, false); // false: keep C as cartesians.
//       EvalDrv(0, pA->l + pB->l, false); // false: keep C as cartesians.
//       EvalDrv(0, pA->l + pB->l, true); // true: transform C to solid harmonics
   }
protected:
   double
      *pOut;
   double const
      *pCoFactorsAB; // <- nCartX_AB x nExpA x nExpB array.
   size_t
      StrideInDenSet, nDenSets,
      StrideC_DenSet;
   size_t
      *piFnC_CartY,
      nFnC_CartY;
   double
      *p_00C_ca; // nFnC_CartY x nDenSets contracted intermediate.

   IR_FORCE_INLINE void SetupCoBuffers(size_t nFnC_) {
      assert(nCartX_ABac == (nCartX(pA->l + pB->l) - nCartX((signed)pA->l-1)));
//       assert(nCartX_ABac == (nCartX(pA->l + pB->l)));

      assert(nDenSets == 1); // need to allocate and handle more stuff in the multi-density version. Don't care enough now.
//       Mem.Alloc(p_A0C_ppp, nCartX_ABac *
      piFnC_CartY = MakeCartesianFnOffsets(pCs, nC, Mem);
      nFnC_CartY = piFnC_CartY[nC];
      Mem.ClearAlloc(p_00C_ca, nFnC_CartY * nDenSets);
      IR_SUPPRESS_UNUSED_WARNING(nFnC_);
   }

   IR_FORCE_INLINE void BeginContractions() { // override
      SetupCoBuffers(nFnC_Total);
   }

   IR_FORCE_INLINE void BeginContraction1(uint /*iExpB*/) { // override
   }

//    IR_FORCE_INLINE void BeginContraction2(uint /*iExpA*/, uint /*iExpB*/){ // override
//    }
   IR_FORCE_INLINE void BeginContraction2(uint iExpA, uint iExpB){ // override
      IR_PREFETCH_R(&pCoFactorsAB[nCartX_ABac * (iExpA + pA->nExp * iExpB)]);
      IR_PREFETCH_W(&p_00C_ca[0]);
   }

   IR_FORCE_INLINE void ContractPrimitiveA0C(double *p_A0C, uint iExpA, uint iExpB, uint iExpC, size_t iC, size_t iFnC, FRawShell const *pC) { // override
      assert(nDenSets == 1);
      size_t
         nCartY_C = nCartY(pC->l);
      TMemoryLock<double>
         pCaC(nCartY_C, &Mem);
      double const
         *IR_RP pCoFactorsAB_ = &pCoFactorsAB[nCartX_ABac * (iExpA + pA->nExp * iExpB)];
      for (size_t iCaC = 0; iCaC < nCartY_C; ++ iCaC) {
         double f = 0.;
         for (size_t iA0 = 0; iA0 < nCartX_ABac; ++ iA0)
            f += pCoFactorsAB_[iA0] * p_A0C[iA0 + nCartX_ABac * iCaC];
         pCaC[iCaC] = f;
      }
//       Mxv(pCaC, 1, p_A0C, nCartX_ABac, 1, pCoFactorsAB_, 1, nCartY_C, nCartX_ABac, false, 1.0);
      // ^- makes things worse... probably not the least because p_A0C is aligned the wrong way for an efficient outer product.

      size_t iFnC_CartY = piFnC_CartY[iC];
      Contract1(&p_00C_ca[iFnC_CartY], pCaC, nDenSets * nCartY_C, pC, iExpC);
      // ^- FIXME: this does not depend on iExpA and iExpB and probably should not be here.
      IR_SUPPRESS_UNUSED_WARNING(iFnC);
   }
   IR_FORCE_INLINE void EndContraction2(uint /*iExpA*/, uint /*iExpB*/) { // override
   }

   IR_FORCE_INLINE void EndContraction1(uint /*iExpA*/) { // override
   }

   IR_FORCE_INLINE void EndContractions() { // override
      // last thing to do: transform C from cartesians to solid harmonics.
      for (size_t iDenSet = 0; iDenSet < nDenSets; ++ iDenSet) {
         TMemoryLock<double>
            ShTmp(nFnC_Total, &Mem);
         for (size_t iC = 0; iC < nC; ++ iC)
            ShTrN_NN(&ShTmp[piFnC[iC]], &p_00C_ca[piFnC_CartY[iC]], pCs[iC].nCo, pCs[iC].l);
         Add2(&pOut[StrideC_DenSet*iDenSet], ShTmp, 1., nFnC_Total);
      }
   }
};


void EvalInt2e3c_ContractAB(double *pOut, FRawShell const *pA, FRawShell const *pB, double const *pDenAB, size_t sa, size_t sb, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   // reverse-transform density (ab| gamma[ab] into cartesian primitives centered on A.
   double
      *pCoFactorsAB;
   Mem.ClearAlloc(pCoFactorsAB, pA->nExp * pB->nExp * nCartX(pA->l + pB->l));
//    if (pA->l == 0 && pB->l == 0)
   TrafoDensityToCaPrim(pCoFactorsAB, pA, pB, 0, pDenAB, sa, sb, Mem);

   // NOTE: make a loop over the Cs here, unpacking stuff? Otherwise ppp buffers might get too large in the driver function.
   //       (well... maybe later when we actually do a proper MXV/MXM contraction. Atm there are no primitives accumulated---
   //        everything contracted in the inner loop.)
   FInt2e3cAcsFactory(pOut, 0, pA, pB, pCoFactorsAB, 0, 1, pCs, nC, Prefactor, pKernel, Mem).Eval();
}



// some auxiliary data for keeping track of (2 zeta)^n-scaled partially
// contracted integral data.
struct FZetaScaledAb {
   int lab; // actual lab of (a0| la==0..lab intemediates
   size_t nCartAB; // nCartX(lab)
   size_t Size; // total size of scaled integral class, including c and co dimensions.
   size_t Offset; // offset in full semi-contracted intermediate.

   FZetaScaledAb(){};
   FZetaScaledAb(int lab_, size_t Offset_, size_t nSets)
      : lab(lab_), nCartAB(ALIGN_CO_SIZE_DBL(nCartX(lab_))), Size(ALIGN_CO_SIZE_DBL(nSets * nCartAB)), Offset(Offset_)
   {}

   size_t SizeAligned() const { return ALIGN_CO_SIZE_DBL(Size); }
};

enum FDerivIntFlags {
   // if set, reduce (ab| components into laplace commutator form: instead
   // of storing all regular 2nd derivative components, compute only
   //
   //     ([NablaA a]b| - (a [NablaB b]|,
   //
   // were NablaA = (d/d[r1x-Ax])^2 + (d/d[r1y-Ay])^2 + (d/d[r1z-Az])^2
   // and NablaB is similarly defined.
   INTFLAG_LaplaceCommutator = 0x01
};

// nz == number of z components == derivative order
// Flags: INTFLAG_*
template<uint DerivOrder, uint Flags>
struct TInt2e3cDerivShellFactory : public FInt3cShellFactory
{

   inline bool HaveFlag(size_t iFlag) const { return 0 != (Flags & iFlag); }

   TInt2e3cDerivShellFactory(
      double *pOut_, size_t StrideA_, size_t StrideB_, size_t StrideC_, size_t StrideDeriv_,
      //uint DerivOrder_, size_t Flags_,
      FRawShell const *pA_, FRawShell const *pB_, FRawShell const *pCs_, size_t nC_, double Prefactor_,
      FIntegralKernel const *pKernel_, FMemoryStack &Mem_)
      : FInt3cShellFactory(pA_, pB_, pCs_, nC_, Prefactor_, pKernel_, Mem_),
        m_pOut(pOut_), StrideA(StrideA_), StrideB(StrideB_), StrideC(StrideC_), StrideDeriv(StrideDeriv_)//, nz(DerivOrder_), Flags(Flags_)
   {
      if (DerivOrder == 0) {
      } else if (DerivOrder == 1) {
         m_pOutA = m_pOut + 0*StrideDeriv;
         m_pOutB = m_pOut + 3*StrideDeriv;
         if (SwapAB)
            std::swap(m_pOutA, m_pOutB);
      } else if (DerivOrder == 2 && !HaveFlag(INTFLAG_LaplaceCommutator)) {
         m_pOutAA = m_pOut +  0*StrideDeriv;
         m_pOutAB = m_pOut +  6*StrideDeriv; // 9 components.
         m_pOutBB = m_pOut + 15*StrideDeriv;
         assert(!SwapAB); // <- output positions for pOutAB need to be fixed in this case.
                          // may need some modifications of the AMRR core.
         if (SwapAB)
            std::swap(m_pOutAA, m_pOutBB);
      } else if (DerivOrder == 2 && HaveFlag(INTFLAG_LaplaceCommutator)) {
         if (SwapAB)
            Prefactor *= -1.;
      } else {
         assert(0);
      }
   }

   void Eval() {
      EvalDrv(0, pA->l + pB->l + DerivOrder, true);
   }
protected:
   double
      // final output positions supplied from the outside
      *m_pOut,
      *m_pOutA, *m_pOutB, // d/dA[x,y,z] and d/dB[x,y,z], respectively
      *m_pOutAA, *m_pOutAB, *m_pOutBB; // 2nd derivatives: 6 dAA components, mixed dAdB, 6 dBB components.
   size_t
      // output addressing information supplied from the outside
      StrideA, StrideB, StrideC, StrideDeriv;

   // intermediates for contractions
//    uint
//       nz; // maximum nz in (2 ZetaX)^nz we need. That's equal to the derivative order.
   size_t
      nFnC_Out, SizeZb, SizeZaZb;
   double
      *p_A0C_ppc,
      *p_A0C_cpc,
      *p_A0C_ccc;
   FZetaScaledAb
      *pInfoZb,
      *pInfoZaZb;
//    size_t
//       Flags;

   FZetaScaledAb &Zb(int izb) { return pInfoZb[izb]; }
   FZetaScaledAb &ZaZb(int iza, int izb) { return pInfoZaZb[iza + (DerivOrder+1)*izb]; }
   double *pZb(int izb) { return p_A0C_cpc + Zb(izb).Offset; }
   double *pZaZb(int iza, int izb) { return p_A0C_ccc + ZaZb(iza,izb).Offset; }

   IR_FORCE_INLINE void SetupCoBuffers(size_t nFnC_) {
      nFnC_Out = nFnC_;

      // first derivatives: need 2 ZetaA [lab+1], 2 ZetaB [lab+1], unscaled [lab-1]
      // 2nd derivatives: need (2za)^2 [lab+2], (2za) (2zb) [lab+2], (2zb)^2 [lab+2],
      //                       (2za) [lab], (2zb) [lab]
      //                       unscaled [lab-2].
      // I didn't check, but my guess is that it just goes on like this. In
      // any case, we have to deal with scaled integral intermediates, the
      // number of which is different for each scaling factor.

      // We arrange the intermedates in matrix form of the scaling scheme:
      //
      //              |   1     2zb    (2zb)^2
      //       -----------------------------
      //       1      |  lab-2  lab    lab+2
      //       2za    |  lab    lab+2   --
      //       (2za)^2|  lab+2   --     --
      //
      // for each class, we produce a nCartX(lab+iz) x nFnC_Out x [co-components]
      // tensor in consecutive memory, which effectively looks like in the simple,
      // non-scaled case.
      //
      // now setup some auxiliary structures for accumulating these
      // contracted scaled intermediates.
      int lab = pA->l + pB->l;
      Mem.Alloc(pInfoZb, DerivOrder+1);
      Mem.Alloc(pInfoZaZb, (DerivOrder+1)*(DerivOrder+1));
      Mem.Align(16);
      SizeZb = 0;
      if (DerivOrder == 2 && HaveFlag(INTFLAG_LaplaceCommutator)) {
         // needs:
         //              |   1     2zb    (2zb)^2
         //       -----------------------------
         //       1      |  --     lab    lab+2
         //       2za    |  lab     --     --
         //       (2za)^2|  lab+2   --     --
         //
         uint iza, izb;
         izb = 0; Zb(izb) = FZetaScaledAb(lab+2, SizeZb, nFnC_Out * pB->nCo); SizeZb += Zb(izb).SizeAligned();
         izb = 1; Zb(izb) = FZetaScaledAb(lab+2, SizeZb, nFnC_Out * pB->nCo); SizeZb += Zb(izb).SizeAligned();
         izb = 2; Zb(izb) = FZetaScaledAb(lab+2, SizeZb, nFnC_Out * pB->nCo); SizeZb += Zb(izb).SizeAligned();
         // ^- hmm... it explodes if I put a lab+0 for iza=1. why?
         SizeZaZb = 0;
         iza = 1; izb = 0; ZaZb(iza,izb) = FZetaScaledAb(lab,   SizeZaZb, nFnC_Out * pB->nCo * pA->nCo); SizeZaZb += ZaZb(iza,izb).SizeAligned();
         iza = 2; izb = 0; ZaZb(iza,izb) = FZetaScaledAb(lab+2, SizeZaZb, nFnC_Out * pB->nCo * pA->nCo); SizeZaZb += ZaZb(iza,izb).SizeAligned();
         iza = 0; izb = 1; ZaZb(iza,izb) = FZetaScaledAb(lab,   SizeZaZb, nFnC_Out * pB->nCo * pA->nCo); SizeZaZb += ZaZb(iza,izb).SizeAligned();
         iza = 0; izb = 2; ZaZb(iza,izb) = FZetaScaledAb(lab+2, SizeZaZb, nFnC_Out * pB->nCo * pA->nCo); SizeZaZb += ZaZb(iza,izb).SizeAligned();
         ZaZb(0,0).Size = 0;
         ZaZb(1,1).Size = 0;
         ZaZb(1,2).Size = 0;
         ZaZb(2,1).Size = 0;
      } else {
         for (uint izb = 0; izb <= DerivOrder; ++ izb) {
   //          uint izb_max = DerivOrder - iza;
   //          int lab_ZaZb = lab - DerivOrder + 2*(iza+izb_max);
            int lab_ZaZb = lab + DerivOrder;
            assert(nCartX_ABac == nCartX(lab_ZaZb));
            // ^- otherwise the Contract1f in EndContraction2() won't work.
            Zb(izb) = FZetaScaledAb(lab_ZaZb, SizeZb, nFnC_Out * pB->nCo);
            SizeZb += Zb(izb).SizeAligned();
         }
         SizeZaZb = 0;
         for (uint iza = 0; iza <= DerivOrder; ++ iza)
            for (uint izb = 0; izb <= DerivOrder; ++ izb) {
               if (iza + izb <= DerivOrder) {
                  int lab_ZaZb = lab - DerivOrder + 2*(iza+izb);
                  ZaZb(iza,izb) = FZetaScaledAb(lab_ZaZb, SizeZaZb, nFnC_Out * pB->nCo * pA->nCo);
               } else {
                  // this one is not present: only data for scaled intermediates
                  // with iza+izb <= DerivOrder
                  ZaZb(iza,izb) = FZetaScaledAb(-1, 0, 0);
               }
               SizeZaZb += ZaZb(iza,izb).SizeAligned();
            }
      }

      if (SwapAB)
         std::swap(StrideA, StrideB);
//       assert(nCartX_ABac == nCartX(pA->l + pB->l + DerivOrder));

      // allocate intermediates for partially contracted A0C integrals.
      Mem.Alloc(p_A0C_ppc, nCartX_ABac * nFnC_Out); // c contracted
      Mem.Alloc(p_A0C_cpc, SizeZb); // c,b contracted & (2 ZetaB)^0..nz scaled.
      Mem.ClearAlloc(p_A0C_ccc, SizeZaZb);
   }

   IR_FORCE_INLINE void BeginContractions() { // override
      SetupCoBuffers(nFnC_Total);
   }

   IR_FORCE_INLINE void BeginContraction1(uint iExpA) { // override
      memset(p_A0C_cpc, 0, SizeZb * sizeof(*p_A0C_cpc));
      IR_SUPPRESS_UNUSED_WARNING(iExpA);
   }

   IR_FORCE_INLINE void BeginContraction2(uint iExpA, uint iExpB){ // override
      memset(p_A0C_ppc, 0, nCartX_ABac * nFnC_Out * sizeof(*p_A0C_cpc));
      IR_SUPPRESS_UNUSED_WARNING2(iExpA, iExpB);
   }

   IR_FORCE_INLINE void ContractPrimitiveA0C(double *p_A0C, uint iExpA, uint iExpB, uint iExpC, size_t iC, size_t iFnC, FRawShell const *pC) { // override
      // this one is just the same as for the basic integrals;
      // c need not be scaled: we can calculate c-derivatives from
      // a and b derivatives using translational invariance.
      assert(nFnC_Out == pC->nSh() * pC->nCo);
      Contract1(&p_A0C_ppc[nCartX_ABac * iFnC], p_A0C,
         nCartX_ABac*pC->nSh(), pC, iExpC);
      IR_SUPPRESS_UNUSED_WARNING3(iExpA, iExpB, iC);
   }

   IR_FORCE_INLINE void EndContraction2(uint iExpA, uint iExpB) { // override
      // contract B and scale with (2 ZetaB)^n.
      double z2b = 1.;
      for (int izb = 0; izb <= (int)DerivOrder; ++ izb) {
         Contract1f(pZb(izb), p_A0C_ppc, z2b, Zb(izb).Size/pB->nCo, pB, iExpB);
         z2b *= 2*pB->pExp[iExpB];
      }
      IR_SUPPRESS_UNUSED_WARNING(iExpA);
   }

   IR_FORCE_INLINE void EndContraction1(uint iExpA) { // override
      // contract A and scale with (2 ZetaA)^n.
      double z2a = 1.;
      for (int iza = 0; iza <= (int)DerivOrder; ++ iza) {
         for (int izb = 0; iza+izb <= (int)DerivOrder; ++ izb) {
            if (ZaZb(iza,izb).Size == 0)
               continue;
            assert(ZaZb(iza,izb).Size/(pA->nCo) <= Zb(izb).Size);
            Contract1fh(pZaZb(iza,izb), ZaZb(iza,izb).nCartAB,
               pZb(izb), Zb(izb).nCartAB, z2a, nFnC_Out * pB->nCo, pA, iExpA);
         }
         z2a *= 2*pA->pExp[iExpA];
      }
   }

   static double *FormDerivA0(double *p0Za, int la, int lab, size_t nSets, FMemoryStack &Mem)
   {
      if (lab < la)
         return 0;
      double *pOut;
      Mem.Alloc(pOut, nCartX(lab-la) * (2*la+1) * nSets);
      // transform A to solid harmonics.
      ShTrA_XfY(pOut, p0Za, la, lab, nSets);
      // ^- FIXME: p0Za step length okay? I don't think that is right
      //    in all cases.
      return pOut;
   }


   static double *FormDerivA1(double *p0Za, double *p2Za, int la, int lab, size_t nSets, FMemoryStack &Mem)
   {
      if (lab < la)
         return 0;
      // form d/dA[x,y,z] (nCartX(lab-la),0| x Slc(la) x nSets, from
      // inputs (2*ZetaA) (a0| with a=0..lab+1 and unscaled (a0| with a = 0..lab-1.
      double
         *pOut;
      Mem.Alloc(pOut, nCartX(lab-la) * (2*la+1) * 3 * nSets);
      AmrrDerivA1(pOut, p0Za, p2Za, lab, la, nSets);
      return pOut;
   }


   static double *FormDerivA2(double *p0Za, double *p2Za, double *p4Za, int la, int lab, size_t nSets, FMemoryStack &Mem)
   {
      // form (d/dA[x,y,z])^2 (nCartX(lab-la),0| x Slc(la) x nSets, from
      // inputs (2*ZetaA)^2 (a0| with a=0..lab+2,
      //        (2 ZetaA) (a0| with a=0..lab, and
      //        unscaled (a0| with a = 0..lab-2.
      if (lab < la)
         return 0;
      // form d/dA[x,y,z] (nCartX(lab-la),0| x Slc(la) x nSets, from
      // inputs (2*ZetaA) (a0| with a=0..lab+1 and unscaled (a0| with a = 0..lab-1.
      double
         *pOut;
      Mem.Alloc(pOut, nCartX(lab-la) * (2*la+1) * 6 * nSets);
      AmrrDerivA2(pOut, p0Za, p2Za, p4Za, lab, la, nSets);
      return pOut;
   }

   static double *FormDerivA0_L1(double *p2Za, double *p4Za, int la, int lab, size_t nSets, FMemoryStack &Mem)
   {
      if (lab < la)
         return 0;
      double
         *pOut;
      Mem.Alloc(pOut, nCartX(lab-la) * (2*la+1) * 1 * nSets);
      AmrrDerivA0L(pOut, p2Za, p4Za, lab, la, nSets);
      return pOut;
   }

   void xOsrrC_B0(double *pOut, double *p0Zb, uint nDerivCompA) {
      for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB)
         for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
            for (size_t iFnC = 0; iFnC < nFnC_Out; ++ iFnC) {
               size_t
                  iFnA = iCoA * pA->nSh(),
                  iFnB = iCoB * pB->nSh();
               std::size_t
                  iOut = iFnA*StrideA + iFnB*StrideB + iFnC*StrideC,
                  iSet = nDerivCompA * (iFnC + (nFnC_Out * (iCoB + pB->nCo * iCoA)));
                  // ^- nDerivCompA is more of a stride than a set index,
                  //    but fits here better.
               OsrrC(
                  &pOut[iOut], StrideA, StrideB,
                  &p0Zb[nCartX(pB->l)*(2*pA->l+1) * iSet],
                  vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh() );
            }
         }
   }

   void xOsrrC_B1(double *pOut, double *p0Zb, double *p2Zb, uint nDerivCompA) {
      for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB)
         for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
            for (size_t iFnC = 0; iFnC < nFnC_Out; ++ iFnC) {
               size_t
                  iFnA = iCoA * pA->nSh(),
                  iFnB = iCoB * pB->nSh();
               size_t
                  iOut = iFnA*StrideA + iFnB*StrideB + iFnC*StrideC,
                  iSet = nDerivCompA * (iFnC + (nFnC_Out * (iCoB + pB->nCo * iCoA)));
               OsrrC_dB1(
                  &pOut[iOut], StrideA, StrideB, StrideDeriv,
                  &p0Zb[nCartX(pB->l-1)*(2*pA->l+1) * iSet],
                  &p2Zb[nCartX(pB->l+1)*(2*pA->l+1) * iSet],
                  vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh() );
            }
         }
   }

   void xOsrrC_B2(double *pOut, double *p0Zb, double *p2Zb, double *p4Zb, uint nDerivCompA) {
      for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB)
         for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
            for (size_t iFnC = 0; iFnC < nFnC_Out; ++ iFnC) {
               size_t
                  iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
               size_t
                  iOut = iFnA*StrideA + iFnB*StrideB + iFnC*StrideC,
                  iSet = nDerivCompA * (iFnC + (nFnC_Out * (iCoB + pB->nCo * iCoA)));
               OsrrC_dB2(
                  &pOut[iOut], StrideA, StrideB, StrideDeriv,
                  &p0Zb[nCartX(pB->l-2)*(2*pA->l+1) * iSet],
                  &p2Zb[nCartX(pB->l  )*(2*pA->l+1) * iSet],
                  &p4Zb[nCartX(pB->l+2)*(2*pA->l+1) * iSet],
                  vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh() );
            }
         }
   }
   // ^- these three can be unified by baking the pNZb arrays (and
   //    StrideDeriv) into an OsrrC helper object and indirecting the calls.
   //    in any case it's not pretty, though. For performance reasons it might
   //    also be advantageous to unify them into a single set of loops,
   //    or even to defer the ShTrX and generate OsrrC-style functions
   //    taking la, lb, AmB, pZaa, pZab, pZbb, pZa, pZb, pZs as input and doing all
   //    the work. We could still provide nSets of data, we would just have to
   //    pre-calculate and store the output pointers for each iSet
   //    (just as calculated above).
   //
   //    this might even help with the regular OsrrC!

   void xOsrrC_B0L(double *pOut, double *p2Zb, double *p4Zb, uint nDerivCompA) {
      // note: OsrrC_dB0L computes the negative laplacian of B.
      for (size_t iCoB = 0; iCoB < pB->nCo; ++ iCoB)
         for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
            for (size_t iFnC = 0; iFnC < nFnC_Out; ++ iFnC) {
               size_t
                  iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
               size_t
                  iOut = iFnA*StrideA + iFnB*StrideB + iFnC*StrideC,
                  iSet = nDerivCompA * (iFnC + (nFnC_Out * (iCoB + pB->nCo * iCoA)));
               OsrrC_dB0L(
                  &pOut[iOut], StrideA, StrideB,
                  &p2Zb[nCartX(pB->l  )*(2*pA->l+1) * iSet],
                  &p4Zb[nCartX(pB->l+2)*(2*pA->l+1) * iSet],
                  vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh());
            }
         }
   }

   IR_FORCE_INLINE void EndContractions() { // override
      size_t
         // for each of those we get a nCartX(lb) x (2*la+1) x nDerivCompA matrix
         // by FormDerivA_().
         nComp = nFnC_Out * pB->nCo * pA->nCo;
      int
         lab = pA->l + pB->l;
      if (DerivOrder == 0) {
         double
            *p_d00_z00 = FormDerivA0(pZaZb(0,0), pA->l, lab, nComp, Mem);
         xOsrrC_B0(m_pOut,  p_d00_z00, 1); // base integral.
      } else if (DerivOrder == 1) {
         // 1st derivative integrals
         double
            *p_d10_z00 = FormDerivA1(pZaZb(0,0), pZaZb(1,0), pA->l, lab, nComp, Mem),
            *p_d00_z01 = FormDerivA0(pZaZb(0,1), pA->l, lab+1, nComp, Mem),
            *p_d00_z00 = FormDerivA0(pZaZb(0,0), pA->l, lab-1, nComp, Mem);
         xOsrrC_B0(m_pOutA + StrideDeriv*0,  p_d10_z00 + nCartX(pB->l)*(2*pA->l+1)*0,3); // d/d[Ax]
         xOsrrC_B0(m_pOutA + StrideDeriv*1,  p_d10_z00 + nCartX(pB->l)*(2*pA->l+1)*1,3); // d/d[Ay]
         xOsrrC_B0(m_pOutA + StrideDeriv*2,  p_d10_z00 + nCartX(pB->l)*(2*pA->l+1)*2,3); // d/d[Az]
         xOsrrC_B1(m_pOutB, p_d00_z00, p_d00_z01, 1); // d/dB[x,y,z]
      } else if (DerivOrder == 2 && !HaveFlag(INTFLAG_LaplaceCommutator)) {
         // 2nd derivative integrals (regular form)
         double
            *p_d20_z00 = FormDerivA2(pZaZb(0,0), pZaZb(1,0), pZaZb(2,0), pA->l, lab, nComp, Mem),
            *p_d10_z01 = FormDerivA1(pZaZb(0,1), pZaZb(1,1), pA->l, lab+1, nComp, Mem),
            *p_d10_z00 = FormDerivA1(pZaZb(0,0), pZaZb(1,0), pA->l, lab-1, nComp, Mem),
            *p_d00_z02 = FormDerivA0(pZaZb(0,2), pA->l, lab+2, nComp, Mem),
            *p_d00_z01 = FormDerivA0(pZaZb(0,1), pA->l, lab  , nComp, Mem),
            *p_d00_z00 = FormDerivA0(pZaZb(0,0), pA->l, lab-2, nComp, Mem);
         for (uint id2A = 0; id2A < 6; ++ id2A)
            xOsrrC_B0(m_pOutAA + StrideDeriv*id2A,   p_d20_z00 + nCartX(pB->l)*(2*pA->l+1)*id2A, 6); // d^2/dA[xx,yy,zz,xy,xz,yz]
         for (uint id1A = 0; id1A < 3; ++ id1A)
            xOsrrC_B1(m_pOutAB + StrideDeriv*3*id1A, p_d10_z00 + nCartX(pB->l-1)*(2*pA->l+1)*id1A,
                                                     p_d10_z01 + nCartX(pB->l+1)*(2*pA->l+1)*id1A,3); // d/dB[x,y,z] d/dA[x,y,z]
         // ^- loop over b derivatives? these offsets are also the same as in xOsrrC_dX. Could just pass iSet or the stride.
         xOsrrC_B2(m_pOutBB, p_d00_z00, p_d00_z01, p_d00_z02, 1); // d^2/dB[xx,yy,zz,xy,xz,yz]
      } else if (DerivOrder == 2 && HaveFlag(INTFLAG_LaplaceCommutator)) {
         // 2nd derivative integrals (Laplace commutator form)
         double
            *p_d20_z00 = FormDerivA0_L1(pZaZb(1,0), pZaZb(2,0), pA->l, lab, nComp, Mem),
            *p_d00_z02 = FormDerivA0(pZaZb(0,2), pA->l, lab+2, nComp, Mem),
            *p_d00_z01 = FormDerivA0(pZaZb(0,1), pA->l, lab  , nComp, Mem);
         // ^- notes:
         //  - OsrrC_lB1 produced shorter code if the laplace and derivative reduction were
         //    exchanged (see comments in MakeOsrrCFns_dB)... but there was something wrong
         //    then. should be checked/fixed.
         //  - find out if there really is no way of doing the laplace thing first and
         //    then running A0L1 and B0L1 through a single xOsrrC_B0.
         uint la = pA->l, lb = pB->l;
         // magic \o/.
         Add2(p_d00_z01, p_d20_z00, 1./(3 + 2*lb), nCartX(lab-la) * (2*pA->l + 1) * nComp);
         xOsrrC_B0L(m_pOut, p_d00_z01, p_d00_z02, 1); // -(d^2/dB[xx] + d^2/dB[yy] + d^2/dB[zz])
      } else {
         assert(0);
      }

      // still need to fix up the Cs from the translational invariance.
      // (or simply not return them at all)
   }
};


// nz == number of z components == derivative order
// Flags: INTFLAG_*
template<uint DerivOrder, uint Flags>
struct TInt2e3cDerivVcaFactory : public TInt2e3cDerivShellFactory<DerivOrder, Flags>
{
   typedef TInt2e3cDerivShellFactory<DerivOrder, Flags>
      FBase;

   TInt2e3cDerivVcaFactory(
      double *pOut_, size_t StrideA_, size_t StrideB_, size_t StrideC_, size_t StrideDeriv_,
      double const *pIn_, size_t StrideIn_, size_t nVecC_,
      //uint DerivOrder_, size_t Flags_,
      FRawShell const *pA_, FRawShell const *pB_, FRawShell const *pCs_, size_t nC_, double Prefactor_,
      FIntegralKernel const *pKernel_, FMemoryStack &Mem_)
      : FBase(pOut_, StrideA_, StrideB_, StrideC_, StrideDeriv_, pA_, pB_, pCs_, nC_, Prefactor_, pKernel_, Mem_),
        pIn(pIn_), StrideIn(StrideIn_), nVecC(nVecC_)
   {}

   void Eval() {
//       return FBase::Eval();
      this->EvalDrv(0, this->pA->l + this->pB->l + DerivOrder, false); // last: do not transform C to solid harmonics.
   }
protected:
   // inputs
   double const
      *pIn;
   size_t
      StrideIn, nVecC;

   // input vector pIn transformed to cartesian primitives & assorted strides.
   double
      *pInCaPrim;
   size_t
      StrideInCaPrim, *pInCaPrimOfs;

   IR_FORCE_INLINE void BeginContractions() { // override
      // transform input vectors to cartesian primitives.
      TrafoVecToCaPrim(pInCaPrim, StrideInCaPrim, pInCaPrimOfs, // <- outputs
         pIn, StrideIn, nVecC, this->pCs, this->nC, this->Mem);

      FBase::SetupCoBuffers(nVecC);
   }

   IR_FORCE_INLINE void ContractPrimitiveA0C(double *p_A0C, uint iExpA, uint iExpB, uint iExpC, size_t iC, size_t iFnC, FRawShell const *pC)
   { // override
      // contract primitive (a0|c) with transformed input vectors over c.
      size_t
         nCompC = nCartY(pC->l);
      double
         *pInC = &pInCaPrim[pInCaPrimOfs[iC] + nCompC * iExpC];
      assert(nVecC == this->nFnC_Out);
      for (size_t iVecC = 0; iVecC < nVecC; ++ iVecC) {
         double
            *c = &this->p_A0C_ppc[this->nCartX_ABac * iVecC],
            *v = &pInC[StrideInCaPrim * iVecC];
         // TODO: check if using real mxv makes a difference here. atm I am
         // still reluctant to introduce an external dependency for this.
         for (size_t ic = 0; ic < nCompC; ++ ic)
            for (size_t iab = 0; iab < this->nCartX_ABac; ++ iab)
               c[iab] += v[ic] * p_A0C[iab + this->nCartX_ABac * ic];
      }
      IR_SUPPRESS_UNUSED_WARNING3(iExpA, iExpB, iFnC);
   }
};


void EvalInt2e3c1d(double *pOut, size_t *Strides,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   TInt2e3cDerivShellFactory<1,0>(pOut, Strides[0], Strides[1], Strides[2], Strides[3], //1, 0,
      pA, pB, pCs, nC, Prefactor, pKernel, Mem).Eval();
}

void EvalInt2e3c2d(double *pOut, size_t *Strides,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   TInt2e3cDerivShellFactory<2,0>(pOut, Strides[0], Strides[1], Strides[2], Strides[3], // 2, 0,
      pA, pB, pCs, nC, Prefactor, pKernel, Mem).Eval();
}

void EvalInt2e3c_kcomm(double *pOut, size_t *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   TInt2e3cDerivShellFactory<2,INTFLAG_LaplaceCommutator>(pOut,
      Strides[0], Strides[1], Strides[2], Strides[3],
      pA, pB, pCs, nC, Prefactor, pKernel, Mem).Eval();
}

void EvalInt2e3c1d_ContractC(double *pOut, size_t *Strides, double const *pIn, size_t StrideIn, size_t nVecIn,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
//    EvalInt2e3c_(pOut, Strides, pA, pB, pCs, nC, Prefactor, pKernel, Mem);
   TInt2e3cDerivVcaFactory<1,0>(pOut, Strides[0], Strides[1], Strides[2], Strides[3], pIn, StrideIn, nVecIn,
      pA, pB, pCs, nC, Prefactor, pKernel, Mem).Eval();
}


//#endif // INCLUDE_OPTIONALS



// output: contracted kernels Fm(rho,T), format: (TotalL+1) x nCoA x nCoC
void Int2e2c_EvalCoKernels(double *pCoFmT, uint TotalL,
    FRawShell const *pA, FRawShell const *pC,
    double PrefactorExt, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   IR_RESUME_CLOCK(TID_EvalCoKernels);
   double
      t = DistSq3(pA->vCen, pC->vCen),
      *pFmT;
   Mem.Alloc(pFmT, TotalL + 1); // FmT for current primitive.

   // pre-evaluate the (1/(2 Alpha))^la and (-1/(2 Gamma))^lc prefactors
   // (yes, it does make a difference)
   IR_RESUME_CLOCK(TID_EvalCoKernels_PrimData);
   double
      *pInv2Alpha = Mem.AllocN(pA->nExp, dbl),
      *pInv2Gamma = Mem.AllocN(pC->nExp, dbl);
   for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
      pInv2Alpha[iExpA] = bool(pA->l)? std::pow(+1.0/(2*pA->pExp[iExpA]), (int)pA->l) : 1.;
      // ^- use Hermites with D Ax := [1/(2 alpha)] \partial/\partial A_i.
   for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
      pInv2Gamma[iExpC] = bool(pC->l)? std::pow(-1.0/(2*pC->pExp[iExpC]), (int)pC->l) : 1.;
      // ^- -1 because \partial_A R \propto -\partial_C R!
   IR_PAUSE_CLOCK(TID_EvalCoKernels_PrimData);

   // loop over primitives (that's all the per primitive stuff there is)
   for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
   for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
   {
      IR_RESUME_CLOCK(TID_EvalCoKernels_PrimData);
      double
         Alpha = pA->pExp[iExpA],
         Gamma = pC->pExp[iExpC],
         InvEta = 1./(Alpha + Gamma),
         Rho = (Alpha * Gamma)*InvEta, // = (Alpha * Gamma)*/(Alpha + Gamma)
         Prefactor = (M_PI*InvEta)*std::sqrt(M_PI*InvEta); // = (M_PI/(Alpha+Gamma))^{3/2}

      Prefactor *= PrefactorExt;
      Prefactor *= pInv2Alpha[iExpA] * pInv2Gamma[iExpC];
      IR_PAUSE_CLOCK(TID_EvalCoKernels_PrimData);

      // calculate derivatives (2 rho D/D[T])^m G0(rho,T) with T = rho (A-C)^2.
      IR_RESUME_CLOCK(TID_EvalCoKernels_TildeGm);
      pKernel->EvalTildeGm(pFmT, Rho, Rho*t, TotalL, Prefactor);
      IR_PAUSE_CLOCK(TID_EvalCoKernels_TildeGm);

      // contract (lamely). However, normally either nCo
      // or nExp, or TotalL (or even all of them at the same time)
      // will be small, so I guess it's okay.
      //
      // If not, one could still go and factor this, or even use
      // two real mxm to do it.
      IR_RESUME_CLOCK(TID_EvalCoKernels_Contract);
      for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
         double CoC = pC->pCo[iExpC + pC->nExp*iCoC];
         if (CoC == 0.)
            continue;
         for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
            double CoAC = CoC * pA->pCo[iExpA + pA->nExp*iCoA];
            if (CoAC == 0.)
               continue;
            Add2(&pCoFmT[(TotalL+1)*(iCoA + pA->nCo*iCoC)],
               pFmT, CoAC, (TotalL+1));
         }
      }
      IR_PAUSE_CLOCK(TID_EvalCoKernels_Contract);
   }
   IR_PAUSE_CLOCK(TID_EvalCoKernels);

   Mem.Free(pFmT);
}

// write (2*la+1) x (2*lc+1) x nCoA x nCoC matrix to final destination.
static void Scatter2e2c(double *IR_RP pOut, size_t StrideA, size_t StrideC,
   double const *IR_RP pIn, size_t la, size_t lc, size_t nComp, size_t nCoA, size_t nCoC, bool Add)
{
   size_t nShA = 2*la+1, nShC = 2*lc+1;
   if (Add) {
      for (size_t iCoC = 0; iCoC < nCoC; ++ iCoC)
         for (size_t iCoA = 0; iCoA < nCoA; ++ iCoA)
            for (size_t iShC = 0; iShC < nShC; ++ iShC)
               for (size_t iShA = 0; iShA < nShA; ++ iShA)
                  pOut[(iShA + nShA*iCoA)*StrideA + (iShC + nShC*iCoC)*StrideC]
                     += pIn[iShA + nShA * (iShC + nShC * nComp * (iCoA + nCoA * iCoC))];
   } else {
      for (size_t iCoC = 0; iCoC < nCoC; ++ iCoC)
         for (size_t iCoA = 0; iCoA < nCoA; ++ iCoA)
            for (size_t iShC = 0; iShC < nShC; ++ iShC)
               for (size_t iShA = 0; iShA < nShA; ++ iShA)
                  pOut[(iShA + nShA*iCoA)*StrideA + (iShC + nShC*iCoC)*StrideC]
                      = pIn[iShA + nShA * (iShC + nShC * nComp *(iCoA + nCoA * iCoC))];
   }
}

// Forms [CartY(TotalLab)] x nCoA x nCoC.
// Allocates memory. Free pOutR.
static void Int2e2c_EvalCoShY(double *&pOutR, unsigned &TotalCo, FRawShell const *pA, FRawShell const *pC, double Prefactor,
   unsigned TotalLab, unsigned LaplaceOrder, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   IR_RESUME_CLOCK(TID_EvalCoShY);
   IR_RESUME_CLOCK(TID_EvalCoShY_Setup);
   FVec3
      R;
   SubVec3(&R[0], pA->vCen, pC->vCen);
   TotalCo = pA->nCo * pC->nCo;

   double
      *pCoFmT, *pDataR_LapC;
   unsigned
      TotalL = TotalLab + 2 * LaplaceOrder;
   Mem.Alloc(pOutR, nCartY(TotalLab) * TotalCo);
   Mem.Alloc(pDataR_LapC, nCartY(TotalL));

   Mem.ClearAlloc(pCoFmT, (TotalL+1) * TotalCo);
   IR_PAUSE_CLOCK(TID_EvalCoShY_Setup);
   Int2e2c_EvalCoKernels(pCoFmT, TotalL, pA, pC, Prefactor, pKernel, Mem);


   IR_RESUME_CLOCK(TID_EvalCoShY_MDRR);
   for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC)
   for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
      // note: if skipping stuff here due to screening, the output must
      // be wiped unless Add == true!
      double
         *pFmT = &pCoFmT[(TotalL+1)*(iCoA + pA->nCo*iCoC)],
         *pDataR_ = &pOutR[nCartY(TotalLab) * (iCoA + pA->nCo*iCoC)];
      if (LaplaceOrder == 0)
         ShellMdrr(pDataR_, pFmT, R[0], R[1], R[2], TotalL);
      else {
         assert(LaplaceOrder == 1);
         ShellMdrr(pDataR_LapC, pFmT, R[0], R[1], R[2], TotalL);
         // note: a simple way of getting higher laplace derivatives is to
         // just apply ShellLaplace function multiple times.
         // At this moment this is not implemented since it is not required.
         ShellLaplace(pDataR_, pDataR_LapC, LaplaceOrder, TotalL - 2);
      }
   }
   IR_PAUSE_CLOCK(TID_EvalCoShY_MDRR);
   IR_PAUSE_CLOCK(TID_EvalCoShY);
}


// evaluate 2-electron 2-center integrals <a|krn * laplace^LaplaceOrder|c>
// note: to obtain the kinetic energy operator, pass an overlap kernel
//       and supply -.5 as Prefactor (ekin = -.5 laplace).
// if add is given: increment the output instead of overwriting it.
// This one is described in: https://doi.org/10.1021/acs.jctc.9b01296
void EvalInt2e2c_LaplaceC(double *pOut, size_t StrideA, size_t StrideC,
    FRawShell const *pA, FRawShell const *pC, double Prefactor, bool Add,
    unsigned LaplaceOrder, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   if (pA->l < pC->l) {
      // ^- FIXME: is that the right way around? Note also that the function works fine
      //    with either order... this is just a question of efficiency in ShTrA_YY. Need to look
      //    up which version is better (factor large la first or factor small la first).
      std::swap(pA, pC);
      std::swap(StrideA, StrideC);
   }
   IR_RESUME_CLOCK(TID_EvalInt2e2c);
   uint
      lc = pC->l, la = pA->l,
      TotalCo;
   double
      *pDataR, *pR1, *pFinal;
   Int2e2c_EvalCoShY(pDataR, TotalCo, pA, pC, Prefactor, la + lc, LaplaceOrder, pKernel, Mem);

   Mem.Alloc(pR1, nCartY(la)*(2*lc+1) * TotalCo);
   Mem.Alloc(pFinal, (2*la+1)*(2*lc+1) * TotalCo);

   IR_RESUME_CLOCK(TID_ShTrA_YY);
   ShTrA_YY(pR1, pDataR, lc, (la + lc), TotalCo);
   ShTrA_YY(pFinal, pR1, la, la, (2*lc + 1)*TotalCo);
   IR_PAUSE_CLOCK(TID_ShTrA_YY);
   // now: (2*la+1) x (2*lc+1) x nCoA x nCoC

   IR_RESUME_CLOCK(TID_Scatter2e2c);
   Scatter2e2c(pOut, StrideA, StrideC, pFinal, la, lc, 1, pA->nCo, pC->nCo, Add);
   IR_PAUSE_CLOCK(TID_Scatter2e2c);

   Mem.Free(pDataR);

   IR_PAUSE_CLOCK(TID_EvalInt2e2c);
}



// evaluate 1st derivative of 2-electron 2-center integrals <a|krn * laplace|c>
// note: to obtain the kinetic energy operator, pass an overlap kernel
//       and supply -.5 as Prefactor (ekin = -.5 laplace).
// if add is given: increment the output instead of overwriting it.
// This one is described in: https://doi.org/10.1021/acs.jctc.9b01296
void EvalInt2e2c1d_LaplaceC( double *pOutAxyz, double *pOutCxyz, size_t StrideA, size_t StrideC, size_t StrideDeriv,
    FRawShell const *pA, FRawShell const *pC, double Prefactor, bool Add,
    unsigned LaplaceOrder, FIntegralKernel const *pKernel, FMemoryStack &Mem )
{
   if (pA->l < pC->l) {
      // ^- FIXME: is that the right way around? Note also that the function works fine
      //    with either order... this is just a question of efficiency in ShTrA_YY. Need to look
      //    up which version is better (factor large la first or factor small la first).
      std::swap(pA, pC);
      std::swap(StrideA, StrideC);
      std::swap(pOutAxyz, pOutCxyz);
   }

   IR_RESUME_CLOCK(TID_EvalInt2e2c);
   uint
      lc = pC->l, la = pA->l,
      TotalCo;
   double
      *pDataR, *pDataR_xyz, *pR1, *pFinal;
   Int2e2c_EvalCoShY(pDataR, TotalCo, pA, pC, Prefactor, la + lc + 1, LaplaceOrder, pKernel, Mem);

   // factor out derivative components.
   Mem.Alloc(pDataR_xyz, nCartY(la+lc) * TotalCo * 3);
   // and now SlmY(la) and SlmY(lb).
   Mem.Alloc(pR1, nCartY(la)*(2*lc+1) * TotalCo * 3);
   Mem.Alloc(pFinal, (2*la+1)*(2*lc+1) * TotalCo * 3);


   // now: nCartY(la+lc+1) x nCoA x nCoC
   IR_RESUME_CLOCK(TID_ShTrA_YY);
   ShTrA_YY(pDataR_xyz, pDataR, 1, (la + lc + 1), TotalCo);

   // now: nCartY(la+lc) x 3 x nCoA x nCoC
   ShTrA_YY(pR1, pDataR_xyz, lc, (la + lc), TotalCo*3);
   ShTrA_YY(pFinal, pR1, la, la, (2*lc + 1)*TotalCo*3);
   IR_PAUSE_CLOCK(TID_ShTrA_YY);
   // now: (2*la+1) x (2*lc+1) x 3 x nCoA x nCoC

   IR_RESUME_CLOCK(TID_Scatter2e2c);
   size_t
      iXyzOff = (2*la+1)*(2*lc+1);
   if (pOutAxyz) {
      for (uint ixyz = 0; ixyz < 3; ++ ixyz)
         // note: this stores only one of the components (deriv wrt. A. Other would be -1 * the input derivative).
         // could this this my simply multiplying *pFinal with -1. and then going through Scatter2e2c again.
         Scatter2e2c(&pOutAxyz[StrideDeriv * ixyz], StrideA, StrideC, &pFinal[ixyz * iXyzOff], la, lc, 3, pA->nCo, pC->nCo, Add);
   }
   if (pOutCxyz) {
      for (uint i = 0; i < 3*(2*la+1)*(2*lc+1)*TotalCo; ++ i)
         pFinal[i] *= -1.;
      for (uint ixyz = 0; ixyz < 3; ++ ixyz)
         // note: this stores only one of the components (deriv wrt. A. Other would be -1 * the input derivative).
         // could this this my simply multiplying *pFinal with -1. and then going through Scatter2e2c again.
         Scatter2e2c(&pOutCxyz[StrideDeriv * ixyz], StrideA, StrideC, &pFinal[ixyz * iXyzOff], la, lc, 3, pA->nCo, pC->nCo, Add);
   }
   IR_PAUSE_CLOCK(TID_Scatter2e2c);

   Mem.Free(pDataR);
   IR_PAUSE_CLOCK(TID_EvalInt2e2c);
}




void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC,
    FRawShell const *pA, FRawShell const *pC, double Prefactor, bool Add,
    FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   if (1) {
      return EvalInt2e2c_LaplaceC(pOut, StrideA, StrideC, pA, pC,
         Prefactor, Add, 0, pKernel, Mem);
   } else {
      assert(!Add);
      double ExpB = 0., CoB = 1., CenB[3] = {0., 0., 0.};
      FRawShell
         ShB(0, &ExpB, 1, &CoB, 1, &CenB[0], 0);
      size_t Strides3[3] = {StrideA, 1, StrideC};
      return EvalInt2e3c(pOut, &Strides3[0], pA, &ShB, pC, 1, Prefactor, pKernel, Mem);
   }
}








// compute (a0|00)^m for m=0...MaxM (inclusive) and 0 <= |a| <= lab (also inclusive).
// Output array is ordered as nCartX(lab) x (MaxM+1)
void OsrrA_MultiM(double *IR_RP pOut, double const *IR_RP pIm, unsigned lab, unsigned MaxM, double PmAx, double PmAy, double PmAz, double PmQx, double PmQy, double PmQz, double rho, double InvZeta)
{
   size_t
      nCartX_AB = nCartX(lab);
   // calculate (a0|00)^m for 0 <= |a| <= lab and m == MaxM, right at its targte location.
   double
      *pA000_MaxM = pOut + MaxM * nCartX_AB;
   OsrrA(pA000_MaxM, pIm + MaxM, lab, PmAx, PmAy, PmAz, PmQx, PmQy, PmQz, rho, InvZeta);
   // for the lower 0 <= m < MaxM, step by step assemble (a0|00)^m from (a0|00)^(m+1)
   // (with |a| <= lab-1) and the scalar (00|00)^m
   for (unsigned m = MaxM - unsigned(1); m < MaxM; --m) {
      double
         *pA000_m = pOut + m * nCartX_AB,
         *pA000_m1 = pOut + (m+1) * nCartX_AB,
         Im = pIm[m];
      OsrrA_StepM(pA000_m, pA000_m1, Im, lab, PmAx, PmAy, PmAz, PmQx, PmQy, PmQz, rho, InvZeta);
   }
}



void EvalInt2e4c(double *pOut, size_t *Strides,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pC, FRawShell const *pD,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   IR_TIME_SECTION(TID_EvalInt2eNc);
   IR_RESUME_CLOCK(TID_InitialSetup);
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   size_t
      StrideA = Strides[0], StrideB = Strides[1], StrideC = Strides[2], StrideD = Strides[3];
   if (pA->l + pB->l < pC->l + pD->l) { // <- more efficient to build larger lab via OsrrA and smaller lcd via OsrrB
      std::swap(pA, pC);
      std::swap(pB, pD);
      std::swap(StrideA, StrideC);
      std::swap(StrideB, StrideD);
   }
   if (pA->l < pB->l) { // <- OsrrC only implemented for la >= lb.
      std::swap(pA, pB);
      std::swap(StrideA, StrideB);
   }
   if (pC->l < pD->l) {
      std::swap(pC, pD);
      std::swap(StrideC, StrideD);
   }
   double const
//       fThreshSab = 0.;
      fThreshSab = 1e-13;
//       fThreshSab = 1e-16;
//       fThreshSab = 1e-11;  // <- that one actually seemed fine. but don't want to risk it without controls...
   FPrimitivePairData
      PrimitivePairDataAB(pA, pB, fThreshSab, Mem),
      PrimitivePairDataCD(pC, pD, fThreshSab, Mem);
   if (PrimitivePairDataAB.nPairsToProcess() * PrimitivePairDataCD.nPairsToProcess() == 0) {
      // either on (ab| side or on |cd) side all contractions or primitive pairs
      // were screened out. write zeros to output locations and return.
      ZeroBlock4(pOut, StrideA, StrideB, StrideC, StrideD, pA->nFn(), pB->nFn(), pC->nFn(), pD->nFn());
      IR_PAUSE_CLOCK(TID_InitialSetup);
      return;
   }
   uint
      lab = pA->l + pB->l,
      lcd = pC->l + pD->l;
   uint
      TotalL = lab + lcd;
   // allocate intermediates
   size_t
      nCartX_AB = nCartX(pA->l + pB->l),
      nCartX_Am1 = nCartX(pA->l-1),
//       nCartX_B = nCartX(pB->l),
//       nShA_CartXB = pA->nSh() * nCartX_B,
      nCartX_ABmA = nCartX_AB - nCartX_Am1;
   size_t
      nCartX_CD = nCartX(pC->l + pD->l),
      nCartX_Cm1 = nCartX(pC->l-1),
//       nCartX_D = nCartX(pD->l),
//       nShC_CartXD = pC->nSh() * nCartX_D,
      nCartX_CDmC = nCartX_CD - nCartX_Cm1;
   size_t
      // size of a single set of (a0|c0), with a in CartX(la,lab) and b in CartX(lc,lcd)
      nA0C0 = nCartX_ABmA * nCartX_CDmC;
   // intermediates for primitive integrals
   double
      // (00|00)^m integrals made from kernel derivatives (-d/dT)^m of (00|00)^0 integral
      *pIm = Mem.AllocN(TotalL + 1, dbl),
      // (a0|00)^m intermediates, for m=0...lcd (inclusive)
      *p_A000m = Mem.AllocN(nCartX_AB * (lcd+1), dbl),
      // primitive (a0|c0)
      *p_A0C0 = Mem.AllocN(nA0C0, dbl),
      *pMemOsrrB = Mem.AllocN(nCartX_AB * nCartX_CD * (lcd+1), dbl);

   // Notes on contraction order:
   // - typically smaller angular momenta are contracted higher than larger ones
   // - for this reason we here start with contracting D, then B, then C, then A
   // - That is just an educated guess, though, and not necessarily the best contraction order
   //   Should probably be handled in a more clever way (such that order is decided dynamically;
   //   however, this also affects the loop order below, so it is not entirely trivial
   //   to make it dynamic...
   //   + best way to do it would probably be to make the function
   //     generate an unrolled "program" with instructions for handling all the primitives,
   //     including the order and where to contract them to with which factors.
   //   + but that's more complicated, and it does not really matter enough atm.

   // intermediates for contractions
   uint
      nCoTotal = pD->nCo * pB->nCo * pC->nCo * pA->nCo;
   IR_SUPPRESS_UNUSED_WARNING(nCoTotal);
   FCoBuffer
      // intermediates (a0|c0) with ABC primitives and D contracted, a = la..lab (size: nA0C0 * pD->nCo))
      p_A0C0_1(nA0C0, pD, Mem),
      // intermediates (a0|c0) with B,D contracted (size: nA0C0 * pD->nCo * pB->nCo)
      p_A0C0_2(p_A0C0_1.size(), pB, Mem),
      // intermediates (a0|c0) with B,C,D contracted (size: nA0C0 * pD->nCo * pB->nCo * pC->nCo)
      p_A0C0_3(p_A0C0_2.size(), pC, Mem),
      // intermediates (a0|c0) with A,B,C,D all contracted (size: nA0C0 * pD->nCo * pB->nCo * pC->nCo * pA->nCo)
      p_A0C0_4(p_A0C0_3.size(), pA, Mem);
   double
      *p_A0C0_cccc = p_A0C0_4.p;

   double const
      *vAmB = &PrimitivePairDataAB.vAmB[0],
      *vCmD = &PrimitivePairDataCD.vAmB[0];
#ifdef IR_ENABLE_KERNEL_RANGE_SCREENING
   double
      fRangeKernelSq = sqr(pKernel->MaxRange());
   double
      fDistSqACK = DistSq3(pA->vCen, pC->vCen) - fRangeKernelSq,
      fDistSqBCK = DistSq3(pB->vCen, pC->vCen) - fRangeKernelSq,
      fDistSqADK = DistSq3(pA->vCen, pD->vCen) - fRangeKernelSq,
      fDistSqBDK = DistSq3(pB->vCen, pD->vCen) - fRangeKernelSq;
#endif // IR_ENABLE_KERNEL_RANGE_SCREENING
   // ^- WARNING: not thought through very hard...

   IR_PAUSE_CLOCK(TID_InitialSetup);
   IR_RESUME_CLOCK(TID_PrimLoop);
   p_A0C0_4.Clear();
   for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
   {
      p_A0C0_3.Clear();
      for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
      {
         p_A0C0_2.Clear();
         for (uint iExpB = 0; iExpB < pB->nExp; ++ iExpB)
         {
            FPrimitivePairEntry const
               *pPairAB = PrimitivePairDataAB.pEntry(iExpA, iExpB);
            if (pPairAB == 0)
               continue; // screened out.
#ifdef IR_ENABLE_KERNEL_RANGE_SCREENING
            if (!IsPrimitiveWithinRange(pA, iExpA, pC, iExpC, fDistSqACK) ||
                !IsPrimitiveWithinRange(pB, iExpB, pC, iExpC, fDistSqBCK))
               continue;
#endif // IR_ENABLE_KERNEL_RANGE_SCREENING

            p_A0C0_1.Clear();
            for (uint iExpD = 0; iExpD < pD->nExp; ++ iExpD)
            {
               FPrimitivePairEntry const
                  *pPairCD = PrimitivePairDataCD.pEntry(iExpC, iExpD);
               if (pPairCD == 0)
                  continue; // screened out.
#ifdef IR_ENABLE_KERNEL_RANGE_SCREENING
               if (!IsPrimitiveWithinRange(pA, iExpA, pD, iExpD, fDistSqADK) ||
                   !IsPrimitiveWithinRange(pB, iExpB, pD, iExpD, fDistSqBDK))
                  continue;
#endif // IR_ENABLE_KERNEL_RANGE_SCREENING
               IR_RESUME_CLOCK(TID_PrimData);
               double
                  Sab = pPairAB->Sab,
                  Scd = pPairCD->Sab,
                  fSabScd = Sab * Scd;
               if (fSabScd < fThreshSab * 1e-2) {
                  IR_PAUSE_CLOCK(TID_PrimData);
                  continue;
               }
               FGaussProduct const
                  *pOvAB = &pPairAB->OvAB,
                  *pOvCD = &pPairCD->OvAB;

               FGaussProduct
                  OvPQ(pOvAB->vCen, pOvAB->Zeta, pOvCD->vCen, pOvCD->Zeta);
               double const
                  *PmQ = OvPQ.vAmB,
                  *PmA = pOvAB->vPmA,
                  *QmC = pOvCD->vPmA,
                  rho = OvPQ.Exp, // [1] (3)
                  T = rho * OvPQ.DistSq,
                  InvZetaEta = OvPQ.InvZeta, // 1/(zeta + eta) = 1/(alpha + beta + gamma + delta)
                  Factor = pow15(M_PI * InvZetaEta) * Sab * Scd * Prefactor; // [1] (7)
               IR_PAUSE_CLOCK(TID_PrimData);

               { IR_TIME_SECTION(TID_EvalGm);
                  // make I[m] = (00|00)^m, m = 0..TotalL (inclusive)
                  pKernel->EvalGm(pIm, rho, T, TotalL, Factor);
               }

               { IR_TIME_SECTION(TID_OsrrA);
                  // make (a0|0)^m for a = 0..lab and m = 0..lcd.
                  OsrrA_MultiM(p_A000m, pIm, lab, lcd, PmA[0], PmA[1], PmA[2],
                     PmQ[0], PmQ[1], PmQ[2], rho, pOvAB->InvZeta);
               }

               double
                  *p_A0C0_pppp;
               { IR_TIME_SECTION(TID_OsrrB);
                  // make (a0|c0) for a = la..lab, c = lc..lcd
                  if (lcd == 0) {
                     p_A0C0_pppp = p_A000m + nCartX_Am1;
                  } else {
                     p_A0C0_pppp = p_A0C0;
                     OsrrB_4c(p_A0C0_pppp, p_A000m, pMemOsrrB, PmQ[0], PmQ[1], PmQ[2],
                        QmC[0], QmC[1], QmC[2], rho, .5*pOvCD->InvZeta, .5*InvZetaEta, pA->l, lab, pC->l, lcd);
                  }
               }

               // primitive (a0|c0) with la <= |a| <= lab and lc <= |c| <= lcd is ready now.
               // add it to its target contractions
               p_A0C0_1.Contract1(p_A0C0_pppp, iExpD);
            } // d exponents
            // p_A0C0_pppc should be done now. Contract B, C, A.
            p_A0C0_2.Contract1(p_A0C0_1, iExpB);
         } // b exponents
         p_A0C0_3.Contract1(p_A0C0_2, iExpC);
      } // c exponents
      p_A0C0_4.Contract1(p_A0C0_3, iExpA);
   } // a exponents
   p_A0C0_4.Finalize1();
   IR_PAUSE_CLOCK(TID_PrimLoop);

   { IR_TIME_SECTION(TID_Finalize);
      // intermediates for OsrrC and SphTr
      double
         // intermediates (c0|ab) nCartX(lc,lcd) x (2*la+1) x x (2*lb+1) x nCo
         *p_C0AB_cccc = Mem.AllocN(nCartX_CDmC * (pA->nSh() * pB->nSh()), dbl);

      size_t
         p_C0AB_StrideA = nCartX_CDmC,
         p_C0AB_StrideB = nCartX_CDmC * pA->nSh();
// hm... one loop order is better for output alignment, the other is better for inner alignment.
// not sure. For small bases the inner-alignment one (with loops in same order as transversal
// in the co buffers above) appeared to work a bit better.
//       for (uint iCoD = 0; iCoD < pD->nCo; ++ iCoD) {
//          for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
//             for (uint iCoB = 0; iCoB < pB->nCo; ++ iCoB) {
//                for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
      for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
         for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
            for (uint iCoB = 0; iCoB < pB->nCo; ++ iCoB) {
               for (uint iCoD = 0; iCoD < pD->nCo; ++ iCoD) {
                  uint
                     iCoTotal = (iCoD + pD->nCo * (iCoB + pB->nCo * (iCoC + pC->nCo * iCoA)));
                  // we have nCartX(lab,la) * nCartX(lcd,lc) at p_A0C0_cccc_ThisCo.
                  // Transform A to solid harmonics by factoring nCartX(lab) into nCartX(lb) x Slm(A).
                  // Then move the angular momentum from a to b and transform b to solid harmonics.
                  // This finishes the (ab| part of the integral assembly, but we still need to
                  // process the |cd) side. So while we are at it, we will also use OsrrC to
                  // transpose the intermediate integrals from (xa|c0) to (c0|ab).
                  {  IR_TIME_SECTION(TID_OsrrC1);
                     for (size_t iSetCDmC = 0; iSetCDmC < nCartX_CDmC; ++ iSetCDmC) {
                        double const
                           *p_A0C0_cccc_ThisCo = &p_A0C0_cccc[nA0C0 * iCoTotal];
                        OsrrC_sha(
                           &p_C0AB_cccc[iSetCDmC], p_C0AB_StrideA, p_C0AB_StrideB,
                           &p_A0C0_cccc_ThisCo[nCartX_ABmA * iSetCDmC],
                           pA->l, vAmB[0], vAmB[1], vAmB[2], pB->l);
                     }
                  }
                  {  IR_TIME_SECTION(TID_OsrrC2);
                     // now do the same process again for the |cd) shell pair, which is now in the fast
                     // axis of the intermediate integral arrays at p_C0AB_cccc, encoded as nCartX(lc,lc+ld).
                     uint
                        iFnC = iCoC * pC->nSh(),
                        iFnD = iCoD * pD->nSh();
                     for (uint iShA = 0; iShA < pA->nSh(); ++ iShA) {
                        for (uint iShB = 0; iShB < pB->nSh(); ++ iShB) {
                           uint
                              iFnA = iShA + iCoA * pA->nSh(),
                              iFnB = iShB + iCoB * pB->nSh();
                           OsrrC_sha(
                              &pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC + iFnD*StrideD], StrideC, StrideD,
                              &p_C0AB_cccc[nCartX_CDmC * (iShA + pA->nSh() * iShB)],
                                 pC->l, vCmD[0], vCmD[1], vCmD[2], pD->l);
                        }
                     }
                  }
               }
            }
         }
      }
   } // timer lock.
}
//       for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
//          for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
//             for (uint iCoB = 0; iCoB < pB->nCo; ++ iCoB) {
//                for (uint iCoD = 0; iCoD < pD->nCo; ++ iCoD) {

//    { IR_TIME_SECTION(TID_Finalize);
//       // intermediates for OsrrC and SphTr
//       double
//          // intermediates (xa|c0) with A,B,C,D contracted and (xa| = nCartX(lb) x (2*la+1)
//          *p_xAC0_cccc = Mem.AllocN(nShA_CartXB * nCartX_CDmC * nCoTotal, dbl),
//          // intermediates (c0|ab) nCartX(lc,lcd) x (2*la+1) x x (2*lb+1) x nCo
//          *p_C0AB_cccc = Mem.AllocN(nCartX_CDmC * (pA->nSh() * pB->nSh()) * nCoTotal, dbl),
//          // intermediates (xc|ab) (nCartX(ld) x (2*lc+1)) x (2*la+1) x x (2*lb+1) x nCo
//          *p_xCAB_cccc = Mem.AllocN(nShC_CartXD * (pA->nSh() * pB->nSh()) * nCoTotal, dbl);
//
//       // transform A to solid harmonics by factoring nCartX(lab) into nCartX(lb) x Slm(A).
//       {  IR_TIME_SECTION(TID_ShTr1);
//          ShTrA_XY(p_xAC0_cccc, p_A0C0_cccc, pA->l, (pA->l + pB->l), nCartX_CDmC * nCoTotal);
//       }
//
//       // we now have nCartX(lb) x nShA x nCartX_CDmC * nCoTotal at p_xAC_ccc.
//       // Now move the angular momentum from a to b and transform b to solid harmonics.
//       // This finishes the (ab| part of the integral assembly, but we still need to
//       // process the |cd) side. So while we are at it, we will also use OsrrC to
//       // transpose the intermediate integrals from (xa|c0) to (c0|ab).
//       {  IR_TIME_SECTION(TID_OsrrC1);
//          size_t
//             p_C0AB_StrideA = nCartX_CDmC,
//             p_C0AB_StrideB = nCartX_CDmC * pA->nSh(),
//             p_C0AB_StrideCo = nCartX_CDmC * (pA->nSh() * pB->nSh());
//          for (uint iCoTotal = 0; iCoTotal < nCoTotal; ++ iCoTotal) {
//             for (size_t iSetCDmC = 0; iSetCDmC < nCartX_CDmC; ++ iSetCDmC) {
//                OsrrC(
//                   &p_C0AB_cccc[iSetCDmC + p_C0AB_StrideCo * iCoTotal], p_C0AB_StrideA, p_C0AB_StrideB,
//                   &p_xAC0_cccc[nShA_CartXB * (iSetCDmC + nCartX_CDmC * iCoTotal)],
//                   vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh());
//             }
//          }
//       }
//
//       // now do the same process again for the |cd) shell pair, which is now in the fast
//       // axis of the intermediate integral arrays at p_C0AB_cccc, encoded as nCartX(lc,lc+ld).
//
//       // transform C to solid harmonics by factoring nCartX(lcd) into nCartX(ld) x Slm(C).
//       {  IR_TIME_SECTION(TID_ShTr2);
//          ShTrA_XY(p_xCAB_cccc, p_C0AB_cccc, pC->l, (pC->l + pD->l), (pA->nSh() * pB->nSh()) * nCoTotal);
//       }
//
//       // we now have nCartX(ld) x nShC x nShA x nShB x nCoTotal at p_xCAB_cccc.
//       // final step: move angular momentum from C to D, solid-harmonic transform D,
//       // and write integrals to final output location.
//       {  IR_TIME_SECTION(TID_OsrrC2);
//          for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
//             for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
//                for (uint iCoB = 0; iCoB < pB->nCo; ++ iCoB) {
//                   for (uint iCoD = 0; iCoD < pD->nCo; ++ iCoD) {
//                      uint
//                         iCoTotal = (iCoD + pD->nCo * (iCoB + pB->nCo * (iCoC + pC->nCo * iCoA))),
//                         iFnC = iCoC * pC->nSh(),
//                         iFnD = iCoD * pD->nSh();
//                      for (uint iShA = 0; iShA < pA->nSh(); ++ iShA) {
//                         for (uint iShB = 0; iShB < pB->nSh(); ++ iShB) {
//                            uint
//                               iFnA = iShA + iCoA * pA->nSh(),
//                               iFnB = iShB + iCoB * pB->nSh();
//                            OsrrC(
//                               &pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC + iFnD*StrideD], StrideC, StrideD,
//                               &p_xCAB_cccc[nShC_CartXD * (iShA + pA->nSh() * (iShB + pB->nSh() * iCoTotal))],
//                                  vCmD[0], vCmD[1], vCmD[2], pD->l, pC->nSh());
//                         }
//                      }
//                   }
//                }
//             }
//          }
//       }
//    } // timer lock.
// }

//          if ( Sab < 1e-10 )
//             continue;
//          if (!IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, fDistSqAB) && Sab > 1e-10) {
//             printf("Problem: A[l=%i exp=%f rg=%f]  B[l=%i  exp=%f  rg=%f]  fDistAB = %10.3f  Sab = %.2e\n",
//                   pA->l, pA->pExp[iExpA], pA->ExpRange(iExpA),
//                   pB->l, pB->pExp[iExpB], pB->ExpRange(iExpB),  std::sqrt(fDistSqAB), Sab);
//          }

//          if (!IsContractionWithinRange(pA, iCoA, pB, iCoB, fDistSqAB)) {
//             uint
//                iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
//             for (uint iFnC = 0; iFnC < nFnC_Total; ++ iFnC)
//                for (uint iB = iFnB; iB != iFnB + pB->nSh(); ++ iB)
//                   for (uint iA = iFnA; iA != iFnA + pA->nSh(); ++ iA)
//                      pOut[iA*StrideA + iB*StrideB + iFnC*StrideC] = 0.;
//             continue;
//          }





unsigned GetSymmetrySignature(int l, size_t c, unsigned Flags)
{
   assert(l >= 0 && unsigned(l) <= std::max(MaxLa,MaxLc) && c <= size_t(2*l));
   size_t
      iSlc1 = nSlmX(l-1) + c;
   size_t const
      nSlcSig = sizeof(iSlcMirrorSymmetrySignatures)/sizeof(iSlcMirrorSymmetrySignatures[0]);
   assert_rt(iSlc1 < nSlcSig);
   return iSlcMirrorSymmetrySignatures[iSlc1];
   IR_SUPPRESS_UNUSED_WARNING(Flags);
}
// ^- this is here because it touches data in IrAmrr.h, the header for the
// generated part of the code.


} // namespace ir
