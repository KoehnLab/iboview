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

/// Internal auxiliary/support data structures & routines for the actual shell drivers.
/// Not meant to be used outside said drivers.

#ifndef IR_DRV_AUX_H
#define IR_DRV_AUX_H

#include "CxVec3.h"
#include "Ir.h"
#include "IrInternal.h"

using ct::FMemoryStack;
using ct::TMemoryLock;

// #define CONTRACT_WITH_BLAS
#ifdef CONTRACT_WITH_BLAS
   #include "CxAlgebra.h" // for DGER (vxv)
#endif



#define ALIGN_CO_MEM(Mem)
#define ALIGN_CO_SIZE_DBL(n) n


// #define IR_ENABLE_KERNEL_RANGE_SCREENING


namespace ir {

typedef ct::TVector3<double>
   FVec3;

struct FGaussProduct
{
   double
      ExpA, ExpB, // exponents of the primitives
      Zeta, // ExpA + ExpB
      InvZeta, // 1/(ExpA+ExpB)
      Exp, // exponent of the primitive product.
      DistSq; // squared distance between A and B
   double
      vCen[4], // this is 'P' for (ab| and 'Q' for |cd) (the center of the overlap distribution)
      vAmB[4], // A-B or C-D
      vPmA[4]; // (A - vCen), in stable form
   FGaussProduct() {};
   IR_NO_INLINE FGaussProduct(double const *vA, double ExpA_, double const *vB, double ExpB_);

   double Sab() {
      return std::exp(-Exp * DistSq);
   }
};


struct FPrimitivePairEntry
{
   FGaussProduct
      OvAB;
   double
      Sab;
   unsigned
      iExpA, iExpB;
};


struct FPrimitivePairData
{
   // note: allocates all local structures on Mem; does *NOT* deallocate anything on termination
   // (assumes stack cleans up itself)
   FPrimitivePairData(FRawShell const *pA, FRawShell const *pB, double fThreshSab, FMemoryStack &Mem);

   inline FPrimitivePairEntry const *pEntry(unsigned iExpA, unsigned iExpB) const {
      return m_ppEntries[iPair(iExpA,iExpB)];
   }
   // note: this may be 0---if all pairs are screened out or the maximum contraction
   // ranges of A and B are stored and do not overlap.
   unsigned nPairsToProcess() const { return m_nPairsToProcess; }
public:
   double
      vAmB[4];
   double
      fDistSqAB;
protected:
   unsigned
      nExpB, nExpA;
   unsigned
      m_nPairsToProcess;
   FPrimitivePairEntry
      // nExpB * nExpA array of pointers to FPrimitivePairEntry structures.
      // returns 0 for screened out primitives.
      **m_ppEntries;
   inline unsigned iPair(unsigned iExpA, unsigned iExpB) const { return iExpB + nExpB * iExpA; }
};

} // namespace ir

namespace ir {



// return x^(3/2).
inline double pow15(double x) {
   return x * std::sqrt(x);
}

inline void SubVec3(double *pOut, double const *pA, double const *pB) {
   pOut[0] = pA[0] - pB[0];
   pOut[1] = pA[1] - pB[1];
   pOut[2] = pA[2] - pB[2];
}

inline bool IsVec3Equal(double const *pA, double const *pB) {
   return pA[0] == pB[0] && pA[1] == pB[1] && pA[2] == pB[2];
}


#ifdef __AVX__
   // pOut += f * pIn
   // ...this version works okay if compiled with -march=native...
   // in fact, in practice it seems to work similarly well as other versions
   // in which stuff is actually aligned.
   static void Add2(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n)
   {
      for (size_t i = 0 ; i < n; ++i)
         pOut[i] += f * pIn[i];
   }

   static void Assign2(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n)
   {
      for (size_t i = 0 ; i < n; ++i)
         pOut[i] = f * pIn[i];
   }
#else
   // pOut += f * pIn
   static void Add2(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n)
   {
      size_t i = 0;
      for ( ; i < (n & ~3); i += 4 ) {
         pOut[i]   += f * pIn[i];
         pOut[i+1] += f * pIn[i+1];
         pOut[i+2] += f * pIn[i+2];
         pOut[i+3] += f * pIn[i+3];
      }
      pOut += i;
      pIn += i;
      switch(n - i) {
         case 3: pOut[2] += f*pIn[2];
         case 2: pOut[1] += f*pIn[1];
         case 1: pOut[0] += f*pIn[0];
         default: break;
      }
   }

   static void Assign2(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n)
   {
      size_t i = 0;
      for ( ; i < (n & ~3); i += 4 ) {
         pOut[i]   = f * pIn[i];
         pOut[i+1] = f * pIn[i+1];
         pOut[i+2] = f * pIn[i+2];
         pOut[i+3] = f * pIn[i+3];
      }
      pOut += i;
      pIn += i;
      switch(n - i) {
         case 3: pOut[2] = f*pIn[2];
         case 2: pOut[1] = f*pIn[1];
         case 1: pOut[0] = f*pIn[0];
         default: break;
      }
   }
   // // pOut += f * pIn
   // static void Add2(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n)
   // {
   //    size_t i = 0;
   //    for ( ; i < (n & ~3); i += 4 ) {
   //       pOut[i]   += f * pIn[i];
   //       pOut[i+1] += f * pIn[i+1];
   //       pOut[i+2] += f * pIn[i+2];
   //       pOut[i+3] += f * pIn[i+3];
   //    }
   //    for ( ; i != n; ++ i ) {
   //       pOut[i] += f * pIn[i];
   //    }
   // }
#endif


#if false && defined(IR_CONTRACT_WITH_XSMM)
// #if defined(IR_CONTRACT_WITH_XSMM)
//    #define MXM1_NN(Out, A, B) {
//       assert(Out.rows() == A.rows() && Out.cols() == B.cols() && A.cols() == B.rows());
//       xsmm_kernel_real_t kernel(LIBXSMM_GEMM_FLAG_NONE, Out.rows(), Out.cols(), A.cols(), (A).outerStride(), (B).outerStride(), (Out).outerStride(), 1.0, 0.0);
//       assert(kernel);
//       kernel(A.data(), B.data(), Out.data());
//    }
   // accumulate matrix nSize  at pIn, corresponding to primitive iExpC,
   // to nSize x nCo at pOut.
   inline void Contract1(double *pOut, double *pIn, size_t nSize, FRawShell const *pC, uint iExpC)
   {
      IR_RESUME_CLOCK(TID_ContractAcc);
//       DGER(nSize, pC->nCo, 1.0, pIn, 1, pC->pCo + iExpC, pC->nExp, pOut, nSize);
      xsmm_kernel_real_double_t kernel(LIBXSMM_GEMM_FLAG_NONE, nSize, pC->nCo, 1, nSize, pC->nExp, nSize, 1.0, 1.0);
      assert(kernel);
      kernel(pIn, pC->pCo + iExpC, pOut);
//       for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
//          double fCoC = pC->fCo(iExpC,iCoC);
//          if (fCoC != 0)
//             Add2(&pOut[nSize*iCoC], pIn, fCoC, nSize);
//       }
//       xsmm_kernel_real_double_t add_kernel(LIBXSMM_GEMM_FLAG_NONE, nSize, 1, 1, nSize, 1, nSize, 1.0, 1.0);
//       for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC) {
//          double fCoC = pC->fCo(iExpC,iCoC);
//          if (fCoC != 0)
// //             Add2(&pOut[nSize*iCoC], pIn, fCoC, nSize);
//             add_kernel(pIn, &fCoC, &pOut[nSize*iCoC]);
//       }
      IR_PAUSE_CLOCK(TID_ContractAcc);
      // ^- this is too slow without the transposed Co matrix. Maybe make a variant which has it?
   }

   // same as Contract1, but additionally multiply with the supplied float.
   void Contract1f(double *pOut, double *pIn, double fScale, size_t nSize, FRawShell const *pC, uint iExpC);

   // same as Contract1f, but process multiple sets with strides and holes.
   void Contract1fh(double *pOut, size_t nSize, double *pIn, size_t StrideIn,
      double fScale, size_t nSets, FRawShell const *pC, uint iExpC);

//    // same as Contract1, but additionally multiply with the supplied float.
//    inline void Contract1f(double *pOut, double *pIn, double fScale, size_t nSize, FRawShell const *pC, uint iExpC)
//    {
//       IR_RESUME_CLOCK(TID_ContractAcc);
//       DGER(nSize, pC->nCo, fScale, pIn, 1, pC->pCo + iExpC, pC->nExp, pOut, nSize);
//       IR_PAUSE_CLOCK(TID_ContractAcc);
//    }
//
//    // same as Contract1f, but process multiple sets with strides and holes.
//    inline void Contract1fh(double *pOut, size_t nSize, double *pIn, size_t StrideIn,
//       double fScale, size_t nSets, FRawShell const *pC, uint iExpC)
//    {
//       IR_RESUME_CLOCK(TID_ContractAcc);
//       assert(nSize <= StrideIn);
//       if (StrideIn == nSize)
//          DGER(nSize*nSets, pC->nCo, fScale, pIn, 1, pC->pCo + iExpC, pC->nExp, pOut, nSize*nSets);
//       else {
//          for (size_t iSet = 0; iSet < nSets; ++ iSet)
//             DGER(nSize, pC->nCo, fScale, &pIn[StrideIn*iSet], 1, pC->pCo + iExpC, pC->nExp, &pOut[nSize*iSet], nSize*nSets);
//       }
//       IR_PAUSE_CLOCK(TID_ContractAcc);
//    }



#elif defined(CONTRACT_WITH_BLAS)
   // TODO: maybe try this with libxsmm? it should be pretty good at this thing.

   // accumulate matrix nSize  at pIn, corresponding to primitive iExpC,
   // to nSize x nCo at pOut.
   inline void Contract1(double *pOut, double *pIn, size_t nSize, FRawShell const *pC, uint iExpC)
   {
      IR_RESUME_CLOCK(TID_ContractAcc);
      DGER(nSize, pC->nCo, 1.0, pIn, 1, pC->pCo + iExpC, pC->nExp, pOut, nSize);
      IR_PAUSE_CLOCK(TID_ContractAcc);
   }

   // same as Contract1, but additionally multiply with the supplied float.
   inline void Contract1f(double *pOut, double *pIn, double fScale, size_t nSize, FRawShell const *pC, uint iExpC)
   {
      IR_RESUME_CLOCK(TID_ContractAcc);
      DGER(nSize, pC->nCo, fScale, pIn, 1, pC->pCo + iExpC, pC->nExp, pOut, nSize);
      IR_PAUSE_CLOCK(TID_ContractAcc);
   }

   // same as Contract1f, but process multiple sets with strides and holes.
   inline void Contract1fh(double *pOut, size_t nSize, double *pIn, size_t StrideIn,
      double fScale, size_t nSets, FRawShell const *pC, uint iExpC)
   {
      IR_RESUME_CLOCK(TID_ContractAcc);
      assert(nSize <= StrideIn);
      if (StrideIn == nSize)
         DGER(nSize*nSets, pC->nCo, fScale, pIn, 1, pC->pCo + iExpC, pC->nExp, pOut, nSize*nSets);
      else {
         for (size_t iSet = 0; iSet < nSets; ++ iSet)
            DGER(nSize, pC->nCo, fScale, &pIn[StrideIn*iSet], 1, pC->pCo + iExpC, pC->nExp, &pOut[nSize*iSet], nSize*nSets);
      }
      IR_PAUSE_CLOCK(TID_ContractAcc);
   }
#else
   // accumulate matrix nSize  at pIn, corresponding to primitive iExpC,
   // to nSize x nCo at pOut.
   void Contract1(double *pOut, double *pIn, size_t nSize, FRawShell const *pC, uint iExpC);

   // same as Contract1, but additionally multiply with the supplied float.
   void Contract1f(double *pOut, double *pIn, double fScale, size_t nSize, FRawShell const *pC, uint iExpC);

   // same as Contract1f, but process multiple sets with strides and holes.
   void Contract1fh(double *pOut, size_t nSize, double *pIn, size_t StrideIn,
      double fScale, size_t nSets, FRawShell const *pC, uint iExpC);

#endif // CONTRACT_WITH_BLAS


   inline void ClearContractBuffer(double *p, size_t n) {
      IR_RESUME_CLOCK(TID_ContractClear);
      memset(p, 0, n*sizeof(*p));
      IR_PAUSE_CLOCK(TID_ContractClear);
   }


   typedef void (*FAddOrAssign2Fn)(double *IR_RP pOut, double const *IR_RP pIn, double fCo, size_t nSize);


   struct FCoBuffer {
      double
         *p;

      // nBase: base size of C primitives (actual buffer size is nBaseSize * pSh->nCo)
      explicit FCoBuffer(size_t nBaseSize, FRawShell const *pSh_, FMemoryStack &Mem);
      void Clear();
//       void Clear() { m_HaveData = false; };
      void Contract1(double const *IR_RP pIn, uint iExp);

      inline void Contract1(FCoBuffer &PrevBuffer, uint iExp) {
         PrevBuffer.Finalize1();
         Contract1(PrevBuffer.p, iExp);
      }

      size_t size() const { return m_nTotalSize; }
      void Finalize1();
   protected:
      FRawShell const
         *m_pSh;
      size_t
         // size of one data set for a primitive of m_pSh
         // (e.g., in a 4c-HRR driver and this being the innermost contraction,
         // this might be the size of a (a0|c0) primitive block; for others it
         // would contain sizes of innermore contractions)
         m_nBaseSize,
         // total size of the contraction buffer (=m_nBaseSize * m_nCo)
         m_nTotalSize;
      unsigned
         // number of contractions and primitives.
         m_nCo,
         m_nExp;
      double
         // transpose of m_pSh->pCo
         *m_pCoT;
      bool
         // flag to mark that this contraction task is for a small buffer, so we
         // should not introduce too much logic into the execution.
         m_CaseIsSmall;
      uint64_t
         // unless m_CaseIsSmall is set, this one is used to flag contractions
         // which have already been written to. It's a bit field.
         m_dwCoFlags;
//       FAddOrAssign2Fn
//          Add2Fn;
//          , Assign2Fn;
      // ^- putting in the fixed-size contraction functions made things substantially worse. interesting.
      inline bool WasCoAlreadySet(uint iCo) const { assert(!m_CaseIsSmall); return 0 != ((m_dwCoFlags >> iCo) & 1); }
      void MarkCoAsSet(uint iCo) { assert(!m_CaseIsSmall); m_dwCoFlags |= (uint64_t(1) << iCo); }
   };

   double *MakeTransposedCoMatrix(FRawShell const *pA, FMemoryStack &Mem);


   void IrInit3c4cTiming(char const *pFuncDesc, FIntegralKernel const *pKernel);
   void IrFinish3c4cTiming(char const *pFuncDesc);



inline void dummy_to_suppress_unused_warnings_in_drv_aux() {
   IR_SUPPRESS_UNUSED_WARNING(Assign2);
   IR_SUPPRESS_UNUSED_WARNING(Add2);
}

} // namespace ir

#endif // IR_DRV_AUX_H
