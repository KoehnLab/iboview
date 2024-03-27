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

#include <new> // for placement new
#include "IrDrvAux.h"

namespace ir {

void ZeroBlock2(double *pOut, size_t StrideA, size_t StrideB, unsigned nFnA, unsigned nFnB)
{
   if (StrideA == 1) {
      for (unsigned iFnB = 0; iFnB < nFnB; ++ iFnB)
         for (unsigned iFnA = 0; iFnA < nFnA; ++ iFnA)
            pOut[iFnA + iFnB*StrideB] = 0.;
   } else if (StrideB == 1) {
      for (unsigned iFnA = 0; iFnA < nFnA; ++ iFnA)
         for (unsigned iFnB = 0; iFnB < nFnB; ++ iFnB)
            pOut[iFnA*StrideA + iFnB] = 0.;
   } else {
      for (unsigned iFnB = 0; iFnB < nFnB; ++ iFnB)
         for (unsigned iFnA = 0; iFnA < nFnA; ++ iFnA)
            pOut[iFnA*StrideA + iFnB*StrideB] = 0.;
   }
}


void ZeroBlock3(double *pOut, size_t StrideA, size_t StrideB, size_t StrideC, unsigned nFnA, unsigned nFnB, unsigned nFnC)
{
   if (StrideA == 1 || StrideB == 1) {
      for (unsigned iFnC = 0; iFnC < nFnC; ++ iFnC)
         ZeroBlock2(&pOut[iFnC * StrideC], StrideA, StrideB, nFnA, nFnB);
   } else {
      for (unsigned iFnB = 0; iFnB < nFnB; ++ iFnB)
         ZeroBlock2(&pOut[iFnB * StrideB], StrideA, StrideC, nFnA, nFnC);
   }
}


// void ZeroBlock4(double *pOut, size_t StrideA, size_t StrideB, size_t StrideC, size_t StrideD, unsigned nFnA, unsigned nFnB, unsigned nFnC, unsigned nFnD)
// {
//    if (StrideA == 1 || StrideB == 1 || StrideC == 1) {
//       for (unsigned iFnD = 0; iFnD < nFnD; ++ iFnD)
//          ZeroBlock3(&pOut[iFnD * StrideD], StrideA, StrideB, StrideC, nFnA, nFnB, nFnC);
//    } else {
//       for (unsigned iFnC = 0; iFnC < nFnC; ++ iFnC)
//          ZeroBlock3(&pOut[iFnC * StrideC], StrideA, StrideB, StrideD, nFnA, nFnB, nFnD);
//    }
// }

void ZeroBlock4(double *pOut, size_t StrideA, size_t StrideB, size_t StrideC, size_t StrideD, unsigned nFnA, unsigned nFnB, unsigned nFnC, unsigned nFnD)
{
   for (unsigned iFnD = 0; iFnD < nFnD; ++ iFnD)
      for (unsigned iFnC = 0; iFnC < nFnC; ++ iFnC)
         for (unsigned iFnB = 0; iFnB < nFnB; ++ iFnB)
            for (unsigned iFnA = 0; iFnA < nFnA; ++ iFnA)
               pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC + iFnD*StrideD] = 0.;
}


FGaussProduct::FGaussProduct(double const *vA, double ExpA_, double const *vB, double ExpB_)
{
   ExpA = ExpA_;
   ExpB = ExpB_;
   Zeta = ExpA + ExpB;
   InvZeta = 1./Zeta;
   Exp = ExpA * ExpB * InvZeta;

   double vA0 = vA[0], vA1 = vA[1], vA2 = vA[2];
   double vB0 = vB[0], vB1 = vB[1], vB2 = vB[2];

   vCen[0] = InvZeta * (ExpA * vA0 + ExpB * vB0);
   vCen[1] = InvZeta * (ExpA * vA1 + ExpB * vB1);
   vCen[2] = InvZeta * (ExpA * vA2 + ExpB * vB2);
   vAmB[0] = vA0 - vB0;
   vAmB[1] = vA1 - vB1;
   vAmB[2] = vA2 - vB2;
   DistSq = sqr(vAmB[0]) + sqr(vAmB[1]) + sqr(vAmB[2]);

   // compute P - A = vCen - vA = (ExpB/(ExpA+ExpB)) * (B - A)
   // The latter form is numerically more stable if alpha/zeta is large,
   // and should always be used in all types of OSRR formulas
   // (see comment after [4] Eq 22)
   double mbi2z = -ExpB * InvZeta;
   vPmA[0] = mbi2z * vAmB[0];
   vPmA[1] = mbi2z * vAmB[1];
   vPmA[2] = mbi2z * vAmB[2];
}


FPrimitivePairData::FPrimitivePairData(FRawShell const *pA, FRawShell const *pB, double fThreshSab, FMemoryStack &Mem)
{
   IR_TIME_SECTION(TID_PrepareShellPairs);
   nExpB = pB->nExp;
   nExpA = pA->nExp;
   SubVec3(&vAmB[0], &pA->vCen[0], &pB->vCen[0]);
   assert(pA->l >= pB->l); // <- otherwise AmB will likely have to be adjusted.
   fDistSqAB = sqr(vAmB[0]) + sqr(vAmB[1]) + sqr(vAmB[2]);

   m_nPairsToProcess = 0;
   if (!IsWithinRange(pA, pB)) {
      // max contraction range of pA and pB does not overlap.
      // FIXME: not sure if this is safe. Not really tested, just added it 2020-04-10
      // because of unused warnings.
      m_ppEntries = 0;
      return;
   }
   Mem.Alloc(m_ppEntries, nExpB * nExpA);
   for (unsigned iExpB = 0; iExpB < pB->nExp; ++ iExpB)
   {
      for (unsigned iExpA = 0; iExpA < pA->nExp; ++ iExpA)
      {
         m_ppEntries[iPair(iExpA, iExpB)] = 0;

         // skip if Dist(A,B) < Range(A) + Range(B)
         if (!IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, fDistSqAB))
            continue;

         FPrimitivePairEntry
            *pEntry;
         Mem.Alloc(pEntry);
         pEntry->iExpA = iExpA;
         pEntry->iExpB = iExpB;
         new (&pEntry->OvAB) FGaussProduct(pA->vCen, pA->pExp[iExpA], pB->vCen, pB->pExp[iExpB]);
         pEntry->Sab = pEntry->OvAB.Sab();
         if (pEntry->Sab < fThreshSab) {
            Mem.Free(pEntry);
            continue;
         }

         m_ppEntries[iPair(iExpA, iExpB)] = pEntry;
         m_nPairsToProcess += 1;
      }
   }
}


#ifndef CONTRACT_WITH_BLAS

#if !defined(IR_CONTRACT_WITH_XSMM)
// accumulate matrix nSize  at pIn, corresponding to primitive iExpC,
// to nSize x nCo at pOut.
void Contract1(double *pOut, double *pIn, size_t nSize, FRawShell const *pC, uint iExpC)
{
   IR_RESUME_CLOCK(TID_ContractAcc);
   for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
      double fCoC = pC->fCo(iExpC,iCoC);
      if (fCoC != 0)
         Add2(&pOut[nSize*iCoC], pIn, fCoC, nSize);
   }
   IR_PAUSE_CLOCK(TID_ContractAcc);
}
#endif

// same as Contract1, but additionally multiply with the supplied float.
void Contract1f(double *pOut, double *pIn, double fScale, size_t nSize, FRawShell const *pC, uint iExpC)
{
   IR_RESUME_CLOCK(TID_ContractAcc);
   for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
      double fCo = fScale * pC->fCo(iExpC,iCoC);
      if (fCo != 0)
         Add2(&pOut[nSize*iCoC], pIn, fCo, nSize);
   }
   IR_PAUSE_CLOCK(TID_ContractAcc);
}

// same as Contract1f, but process multiple sets with strides and holes.
void Contract1fh(double *pOut, size_t nSize, double *pIn, size_t StrideIn,
   double fScale, size_t nSets, FRawShell const *pC, uint iExpC)
{
   IR_RESUME_CLOCK(TID_ContractAcc);
   assert(nSize <= StrideIn);
   if (StrideIn == nSize)
      return Contract1f(pOut, pIn, fScale, nSize*nSets, pC, iExpC);
   for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
      double fCo = fScale * pC->fCo(iExpC,iCoC);
      if (fCo != 0) {
         for (size_t iSet = 0; iSet < nSets; ++ iSet) {
            Add2(&pOut[nSize*(iSet + nSets*iCoC)], &pIn[StrideIn*iSet], fCo, nSize);
         }
      }
   }
   IR_PAUSE_CLOCK(TID_ContractAcc);
}

#endif // CONTRACT_WITH_BLAS

// template<unsigned n>
// static void Add2T(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n_)
// {
//    assert(n == n_); // <- ignored. These ones have static sizes.
//    for (size_t i = 0 ; i < n; ++i)
//       pOut[i] += f * pIn[i];
// }

// template<unsigned n>
// static void Assign2T(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n_)
// {
//    assert(n == n_); // <- ignored. These ones have static sizes.
//    for (size_t i = 0 ; i < n; ++i)
//       pOut[i] = f * pIn[i];
// }

// static FAddOrAssign2Fn
//    g_Add2Fns[] = {
//       Add2T< 0>, Add2T< 1>, Add2T< 2>, Add2T< 3>, Add2T< 4>, Add2T< 5>, Add2T< 6>, Add2T< 7>, Add2T< 8>, Add2T< 9>,
//       Add2T<10>, Add2T<11>, Add2T<12>, Add2T<13>, Add2T<14>, Add2T<15>, Add2T<16>, Add2T<17>, Add2T<18>, Add2T<19>,
//       Add2T<20>, Add2T<21>, Add2T<22>, Add2T<23>, Add2T<24>
//    };
// static FAddOrAssign2Fn
//    g_Assign2Fns[] = {
//       Assign2T< 0>, Assign2T< 1>, Assign2T< 2>, Assign2T< 3>, Assign2T< 4>, Assign2T< 5>, Assign2T< 6>, Assign2T< 7>, Assign2T< 8>, Assign2T< 9>,
//       Assign2T<10>, Assign2T<11>, Assign2T<12>, Assign2T<13>, Assign2T<14>, Assign2T<15>, Assign2T<16>, Assign2T<17>, Assign2T<18>, Assign2T<19>,
//       Assign2T<20>, Assign2T<21>, Assign2T<22>, Assign2T<23>, Assign2T<24>
//    };
// ^-- this made things *much* worse. interesting. I left them in here as a
//     "don't try this again" reminder.


FCoBuffer::FCoBuffer(size_t nBaseSize, FRawShell const *pSh_, FMemoryStack &Mem)
{
   m_pSh = pSh_;
   m_nCo = pSh_->nCo;
   m_nExp = pSh_->nExp;

   m_nBaseSize = nBaseSize;
   m_nTotalSize = nBaseSize * pSh_->nCo;
   Mem.Alloc(p, m_nTotalSize);

   m_CaseIsSmall = (m_nBaseSize <= 32);

   // if there are more contractions than we have bits in our bitfield to track
   // which contractions have been written to, mark this case as "small"
   // (regardless of its actual size) to disable the bitfield mechanism
   // altogether. That case is extremely unlikely in practice (would need more
   // than 64 contractions!), and this hack around it should result in correct
   // code in either case.
   m_CaseIsSmall |= (m_nCo > 8*sizeof(m_dwCoFlags));

//    Add2Fn = Add2;
// //    Assign2Fn = Assign2;
//    if (m_CaseIsSmall && m_nBaseSize < sizeof(g_Add2Fns)/sizeof(g_Add2Fns[0])) {
//       Add2Fn = g_Add2Fns[m_nBaseSize];
// //       assert(sizeof(g_Add2Fns) == sizeof(g_Assign2Fns));
// //       Assign2Fn = g_Assign2Fns[m_nBaseSize];
//    }
}


void FCoBuffer::Contract1(double const *IR_RP pIn, uint iExp)
{
   IR_TIME_SECTION(TID_ContractAcc);
   if (m_CaseIsSmall) {
      // contraction task has a small base size---don't think about it,
      // just clear out the entire buffer on Clear(), and in Contract1,
      // execute the addition operation unconditionally.
      for (uint iCo = 0; iCo < m_nCo; ++ iCo) {
         double fCo = m_pSh->fCo(iExp, iCo);
         if (fCo != 0)
            Add2(&p[m_nBaseSize*iCo], pIn, fCo, m_nBaseSize);
      }
   } else {
      // this is *not* a small contraction case; in this case it may make sense to
      // not zero out the entire contraction buffer, but rather remember which
      // parts of it we have already written to and which ones we have not.
      // (we do this with m_pCoFlags)
      // Note that in extreme cases (e.g., outer buffer in 2e4c integrals), the
      // contraction buffer at this point could have thousands of entries, and
      // there may be millions of cycles between the Clear() and the Contract1()
      // operation here.
      for (uint iCo = 0; iCo < m_nCo; ++ iCo) {
         double fCo = m_pSh->fCo(iExp, iCo);
         if (fCo != 0) {
            double *IR_RP pOut = &p[m_nBaseSize*iCo];
            if (WasCoAlreadySet(iCo)) {
               Add2(pOut, pIn, fCo, m_nBaseSize);
            } else {
               MarkCoAsSet(iCo);
               Assign2(pOut, pIn, fCo, m_nBaseSize);
            }
         }
      }
   }
}

void FCoBuffer::Clear()
{
   IR_TIME_SECTION(TID_ContractClear);
   if (m_CaseIsSmall)
      // small case. Probably Clear() close to Contract1, so just
      // memset the target unconditionally.
      ClearContractBuffer(p, m_nTotalSize);
   else {
      // a potentially large case. Do not zero out the buffers, but mark that
      // none of the target contractions in the buffer have been writen to as of
      // yet. In this case, the first contribution will be stored with an
      // Assign2 instead of an Add2.
      m_dwCoFlags = 0;
   }
}

void FCoBuffer::Finalize1()
{
   if (m_CaseIsSmall) {
      // this case does not need any special handling. Even if no contraction
      // actually ended up being emitted, since we here zero out the contraction
      // buffer at the start, this would be fine anyway.
   } else {
      IR_TIME_SECTION(TID_ContractClear);
      for (uint iCo = 0; iCo < m_nCo; ++ iCo) {
         if (!WasCoAlreadySet(iCo)) {
            // this contraction was never written to. Zero it out.
            MarkCoAsSet(iCo);
            for (size_t i = 0; i < m_nBaseSize; ++ i)
               p[m_nBaseSize*iCo + i] = 0;
         }
      }
   }
}


double *MakeTransposedCoMatrix(FRawShell const *pA, FMemoryStack &Mem)
{
   double
      *pCoT;
   Mem.Alloc(pCoT, pA->nCo * pA->nExp);
   // Matrix transposition (out-of-place form).
   // LIBXSMM_API void libxsmm_otrans(void* out, const void* in, unsigned int typesize,
   //   libxsmm_blasint m, libxsmm_blasint n, libxsmm_blasint ldi, libxsmm_blasint ldo);
// #ifdef IR_USE_XSMM
//    libxsmm_otrans(pCoT, pA->pCo, sizeof(double), pA->nExp, pA->nCo, pA->nExp, pA->nCo);
// #else
   for (unsigned iCo = 0; iCo < pA->nCo; ++ iCo)
      for (unsigned iExp = 0; iExp < pA->nExp; ++ iExp)
         pCoT[iCo + pA->nCo*iExp] = pA->pCo[iExp + pA->nExp * iCo];
// #endif // IR_USE_XSMM
   return pCoT;
}

// #ifdef IR_CONTRACT_WITH_XSMM
//    double const
//       *pCoA_T = MakeTransposedCoMatrix(pA, Mem),
//       *pCoB_T = MakeTransposedCoMatrix(pB, Mem),
//       *pCoC_T = MakeTransposedCoMatrix(pC, Mem),
//       *pCoD_T = MakeTransposedCoMatrix(pD, Mem);
//    xsmm_kernel_real_double_t
//       kernel_coD(LIBXSMM_GEMM_FLAG_NONE, nA0C0, pD->nCo, 1, nA0C0, 1, nA0C0, 1.0, 1.0),
//       kernel_coB(LIBXSMM_GEMM_FLAG_NONE, nA0C0*pD->nCo, pB->nCo, 1, nA0C0*pD->nCo, 1, nA0C0*pD->nCo, 1.0, 1.0),
//       kernel_coC(LIBXSMM_GEMM_FLAG_NONE, nA0C0*pD->nCo*pB->nCo, pC->nCo, 1, nA0C0*pD->nCo*pB->nCo, 1, nA0C0*pD->nCo*pB->nCo, 1.0, 1.0),
//       kernel_coA(LIBXSMM_GEMM_FLAG_NONE, nA0C0*pD->nCo*pB->nCo*pC->nCo, pA->nCo, 1, nA0C0*pD->nCo*pB->nCo*pC->nCo, 1, nA0C0*pD->nCo*pB->nCo*pC->nCo, 1.0, 1.0);
//       //       xsmm_kernel_real_t kernel(LIBXSMM_GEMM_FLAG_NONE, Out.rows(), Out.cols(), A.cols(), (A).outerStride(), (B).outerStride(), (Out).outerStride(), 1.0,
//    assert_rt(kernel_coD);
//    assert_rt(kernel_coB);
//    assert_rt(kernel_coC);
//    assert_rt(kernel_coA);
// #endif // IR_CONTRACT_WITH_XSMM

// #ifdef IR_CONTRACT_WITH_XSMM
//                { IR_TIME_SECTION(TID_ContractAcc);
//                  kernel_coD(p_A0C0_pppp, pCoD_T + iExpD*pD->nCo, p_A0C0_pppc);
//                }
// #else
//                Contract1(p_A0C0_pppc, p_A0C0_pppp, nA0C0, pD, iExpD);
// #endif



#ifdef IR_TIMING
   static double s_Ir34cTiming_t0;
#endif // RDTSC_TIMING

// void IrInit3c4cTiming(char const *pFuncDesc, FRawBasis const &Basis1, FRawBasis const &Basis2, FIntegralKernel const *pKernel)
void IrInit3c4cTiming(char const *pFuncDesc, FIntegralKernel const *pKernel)
{
#ifdef IR_TIMING
   ResetClocks();
//    printf("\n-- IR-TIMING: 2ix matrix <row|%s|col> shape=(%i, %i) nOmpThreads=%i\n", Krn2i.Desc().c_str(), int(BasisRow.nFn()), int(BasisCol.nFn()), omp_get_max_threads());
//    printf("\n-- IR-TIMING: %s <ab|%s|cd> bases: (%i, %i)\n", pFuncDesc, pKernel->Desc().c_str(), int(Basis1.nFn()), int(Basis2.nFn()));
   printf("\n-- IR-TIMING: %s <ab|%s|cd>\n", pFuncDesc, pKernel->Desc());
   s_Ir34cTiming_t0 = 0;
   s_Ir34cTiming_t0 -= ir::GetTicks1();
   IR_RESUME_CLOCK(TID_TimingEnvelope);
#else
   IR_SUPPRESS_UNUSED_WARNING(pFuncDesc);
   IR_SUPPRESS_UNUSED_WARNING(pKernel);
#endif // IR_TIMING
}

void IrFinish3c4cTiming(char const *pFuncDesc)
{
#ifdef IR_TIMING
   IR_PAUSE_CLOCK(TID_TimingEnvelope);
   s_Ir34cTiming_t0 += ir::GetTicks1();
   using namespace ir;
   InitTimingCode();
   SetRefClock(TID_TimingEnvelope);
   PrintClock("EvalInt2eNc", TID_EvalInt2eNc);
   SetRefClock(TID_EvalInt2eNc);
//    printf("^- I don't get it... if I divide CoShY by its fraction on the ref clock, i get the right number. But not if doing it raw.");
   // ^- UPDATE: ...that's because the first time it runs, it does the
   //    CPU frequency and overhead test. GNAAAA.
//    SetRefClock(TID_EvalInt2e2c);
   PrintClock("| InitSetup", TID_InitialSetup);
   PrintClock("| ShellPairs", TID_PrepareShellPairs);
   PrintClock("| FTCS-Xia/Wiab", TID_Ftcs_ShellPairExp);
   PrintClock("| PrimLoop", TID_PrimLoop);
   PrintClock("| | PrimData", TID_PrimData);
   PrintClock("| | KernelFn", TID_EvalGm);
//    FTscTime tPTC = TscGetClock(TID_EvalCoKernels_PrimData) + TscGetClock(TID_EvalCoKernels_TildeGm) + TscGetClock(TID_EvalCoKernels_Contract);
//    PrintClock("| | | (PTC)", tPTC, TID_EvalCoKernels);
   PrintClock("| | OsrrA", TID_OsrrA);
   PrintClock("| | OsrrB", TID_OsrrB);
   PrintClock("| | FTCS-Wabcd", TID_Ftcs_Wabcd);
   PrintClock("| | FTCS-EvalPrim", TID_Ftcs_EvalPrim);
   SetRefClock(TID_Ftcs_EvalPrim);
   PrintClock("| | | TildePVg", TID_Ftcs_TildePVg);
   PrintClock("| | | PVxPrime", TID_Ftcs_PVxPrime);
   PrintClock("| | | Ux,Uy,Uz", TID_Ftcs_Ux);
   PrintClock("| | | Uyz", TID_Ftcs_Uyz);
   PrintClock("| | | Assembly", TID_Ftcs_Assembly);
   PrintClock("| | | | DCCS deferred ImUzUy", TID_Dccs_ImUzUy);
   PrintClock("| | | | | DCCS convolution", TID_Dccs_Convolution);
   PrintClock("| | | | Inner assembly", TID_Dccs_InnerAssembly);
   SetRefClock(TID_EvalInt2eNc);

   FTscTime tCoTotal = TscGetClock(TID_ContractAcc) + TscGetClock(TID_ContractClear);
   PrintClock("| | Contract", tCoTotal, -1); // TID_Contract
   PrintClock("| | | Clear", TID_ContractClear);
   PrintClock("| | | Acc", TID_ContractAcc);
   PrintClock("| Finalize", TID_Finalize);
   SetRefClock(TID_Finalize);
   PrintClock("| | ShTr1", TID_ShTr1);
   PrintClock("| | OsrrC1", TID_OsrrC1);
   PrintClock("| | ShTr2", TID_ShTr2);
   PrintClock("| | OsrrC2", TID_OsrrC2);
   SetRefClock(TID_EvalInt2eNc);

   SetRefClock(-1);
   std::string sDesc = std::string(pFuncDesc) + " (total, real):";
   printf("   %-48s %11.3f msec\n", sDesc.c_str(), 1000.*s_Ir34cTiming_t0);
   PrintClock("Operation (total, tsc)", TID_TimingEnvelope);
   SetRefClock(-1);
#else
   IR_SUPPRESS_UNUSED_WARNING(pFuncDesc);
#endif // IR_TIMING
}


} // namespace ir
