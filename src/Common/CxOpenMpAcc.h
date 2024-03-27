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

/* Intended usage of TOmpAccBlock and FMemoryStackArray (both can also be used
 * independently) of each other):

   // - This object indicates that we mean to compute a block of numbers in a
   //   thread-parallel fashion. And in such a way that the result of each
   //   thread is added up at the end. Essentially, we wish to parallelize this:
   //
   //   for (size_t iJob = 0; i < nJobs; ++ iJob) {
   //      for (size_t i = 0; i <= nSize; ++ i) {
   //         Block[i] += (something computed)
   //      }
   //   }
   // - To afford this without synchronization inside the parallelized loop,
   //   TOmpAccBlock does the following:
   //
   //   + Allocates a memory block of nThreads * nSize, and zeros it out. Idea
   //     is that each thread has its independent data storage which can be
   //     handled independently of what the other threads are doing.
   //   + Provides a simple way to access the block local to each thead (via
   //     FockAcc.pTls() inside the omp parallel section)
   //   + An automated way to sum up all the intermediate values and store them
   //     to the target location (via FockAcc.Join() after the parallel
   //     section).
   size_t
      nSize = 1000;
   TMemoryLock<double>
      // target object of `nSize` doubles we wish to compute here. In this form
      // is is allocated on a FMemoryStack object, but this is not required for
      // using TOmpAccBlock---any continuous block in memory is fair game.
      pFock(nSize, &Mem);
   TOmpAccBlock<double>
      // Create one size `nSize` accumulation block per thread, and remember
      // that at the end of the parallel section, they should be accumulated to
      // target location `pFock`. Memory for this is allocated on `Mem`.
      //
      // - OMPACC_ClearTarget: indicates that the nSize*sizeof(T) length memory
      //   block at pFock should set to zero before adding the new values in
      //   Join. If not given, .Join() will add (accumulate) the results of the
      //   thread computations to whatever is in the target block at that
      //   moment (in this current case, that is just random, as pFock is not
      //   cleared).
      FockAcc(&pFock[0], nSize, OMPACC_ClearTarget, Mem);
   {
      // FMemoryStackArray splits the (entire) remaining space in `Mem` into the
      // `nThreads` slices, and attaches a new existing-memory-buffer based
      // temporary FMemoryStack2 object to each of them. All the stacks are
      // entirely independent, and therefore need no synchronization as long as
      // each thread always uses the same one.
      ct::FMemoryStackArray MemStacks(Mem);
      #pragma omp parallel for schedule(dynamic)
      for (int iGridPt_ = 0; iGridPt_ < int(Points.size()); ++ iGridPt_) {
         size_t iGridPt = size_t(iGridPt_);
         FMemoryStack
            // .GetStackOfThread() returns the sub-stack for this particular
            // thread. Warning: do not use the original `Mem` in here.
            &TheadLocalMem = MemStacks.GetStackOfThread();
         double
            // similarly, .pTls() provides a pointer to the `nSize` length
            // accumulation block for pFock which is associated with the current
            // thread (thread-local storage)
            *pThreadLocalFock = FockAcc.pTls();
         // ...
         // ...
         // add things to pThreadLocalFock in any way/shape/form you
         // need. Both `TheadLocalMem` and `pThreadLocalFock` are unique
         // objects controlling thread-local memory blocks, and therefore
         // do not interfere with the other threads.
         // ...
         // ...
      }
   }
   // Parallel section is done. So .Join now performs the following operations:
   // - If the OMPACC_ClearTarget flag was given (here it was), clear out target
   //   location `pFock` which was supplied to TOmpAccBlock at the time of its
   //   creation (that's just a `memset(pFock, 0, sizeof(pFock[0])*nSize)`-type
   //   operation.
   // - Add the content of all the thread-local memory blocks to the target
   //   location at `pFock`
   // Both of those operations are omp parallelized, too.
   FockAcc.Join();
 */
#ifndef CX_OPENMP_ACC_H
#define CX_OPENMP_ACC_H

#include <memory> // for placement new
#include "CxMemoryStack.h"
#include "CxOpenMpProxy.h"

namespace ct {

enum FOmpAccBlockFlags {
   /// if set, target buffer will be set to zero before adding the results of
   /// the thread-parallel computations in Join.
   OMPACC_ClearTarget = 0x01,
   /// for types which do not like being memset(...,0,...) ... call a placement
   /// new on each point of the target location to invoke the type's default
   /// constructor, rather than zeroing it out as normally.
   OMPACC_CallConstructor = 0x02
};

/// Allocates and manages one independent (nSize)-shape memory block per thread,
/// all of which are added to `pTarget` once Join() is called.
///
/// See longer usage example (in combination with FMemoryStackArray) at top
/// of file `CxOpenMpAcc.h`.
template<class FScalar>
struct TOmpAccBlock
{
   FScalar *pTls(); // return thread-local block of the storage. Must be called from within parallel block.
   void Join(); // copy to target by summation; must be called OUTSIDE of parallel block.
   // allocate and clear thread-local storage for all threads. Must be called OUTSIDE of parallel block.
   void Init(FScalar *pTarget, size_t nSize, unsigned Flags, FMemoryStack &Mem);

   TOmpAccBlock() : m_pTarget(0), m_pTlsData(0), m_nSize(0), m_Flags(0) {};
   TOmpAccBlock(FScalar *pTarget_, size_t nSize_, unsigned Flags_, FMemoryStack &Mem_) : m_pTarget(0), m_pTlsData(0), m_nSize(0), m_Flags(Flags_) { Init(pTarget_, nSize_, Flags_, Mem_); }
protected:
   FScalar
      *m_pTarget,
      *m_pTlsData; // m_nAlignedSize * m_nThreads
   size_t
      m_nSize,
      m_nAlignedSize,
      m_nThreads;
   unsigned
      m_Flags;
   size_t nScalarBlockSize() const { return 128; }
   size_t nTotalSize() const { return m_nAlignedSize * m_nThreads; }
//    size_t nScalarBlocksTotal() const { return nTotalSize() / nScalarBlockSize(); }
   size_t nScalarBlocks(size_t nSize) { return (nSize + nScalarBlockSize() - 1)/nScalarBlockSize(); }

   inline void ClearBlock(FScalar *p, size_t n);
};

typedef TOmpAccBlock<double>
   FOmpAccBlock;

template<class FScalar>
FScalar *TOmpAccBlock<FScalar>::pTls() {
   return &m_pTlsData[omp_get_thread_num() * m_nAlignedSize];
}

template<class FScalar>
void TOmpAccBlock<FScalar>::Init(FScalar *pTarget, size_t nSize, unsigned Flags, FMemoryStack &Mem)
{
   assert(m_nSize == 0);
   if (pTarget == 0)
       nSize = 0;
   m_nSize = nSize;
   m_nAlignedSize = AlignSizeT(nSize, CX_DEFAULT_MEM_ALIGN/sizeof(*m_pTarget));
   m_Flags = Flags;
   if (m_nSize != 0)
      m_pTarget = pTarget;
   else {
      m_pTarget = 0;
      m_nAlignedSize = 0;
      m_pTlsData = 0;
   }
   if (!m_pTarget)
      return;
   m_nThreads = omp_get_max_threads();
   Mem.Alloc(m_pTlsData, nTotalSize());

   // clear thread-local data. Or should I align it with the thread ids? may be better for cache purposes.
   #pragma omp parallel for
   for (int iBlock = 0; iBlock < int(nScalarBlocks(nTotalSize())); ++ iBlock) {
      FScalar
         *pBlock = &m_pTlsData[nScalarBlockSize() * size_t(iBlock)],
         *pBlockEnd = std::min(pBlock + nScalarBlockSize(), m_pTlsData + nTotalSize());
      ClearBlock(pBlock, pBlockEnd - pBlock);
   }
}

template<class FScalar>
static void Add2T(FScalar *IR_RP r, FScalar const *IR_RP x, FScalar f, size_t n )
{
   size_t
      i = 0;
   for ( ; i < n; ++ i ) {
      r[i] += f * x[i];
   }
}

template<class FScalar>
void TOmpAccBlock<FScalar>::ClearBlock(FScalar *p, size_t n)
{
   if (!bool(m_Flags & OMPACC_CallConstructor)) {
      // many normal scalar types are okay with this one... (in particlar, IEEE
      // doubles are). just set all bits to zero.
      memset(p, 0, sizeof(*p) * n);
   } else {
      // ...but some might not be. E.g., boost::multiprecision::cpp_bin_float
      // breaks if it is memset over. In this case, clear the elements by
      // default-constructing them.
      for (size_t i = 0; i < n; ++ i) {
         new (&p[i]) FScalar();
      }
   }
}



template<class FScalar>
void TOmpAccBlock<FScalar>::Join()
{
   if (!m_pTarget) return;

   // horizontal sum. in this loop, each block is summed up across the previous processor dimension.
   #pragma omp parallel for
   for (int iBlock = 0; iBlock < int(nScalarBlocks(m_nSize)); ++ iBlock) {
      size_t
         iBlockOffs = nScalarBlockSize() * size_t(iBlock);
      FScalar
         *pBlock = &m_pTlsData[iBlockOffs],
         *pBlockEnd = std::min(pBlock + nScalarBlockSize(), m_pTlsData + m_nSize);
      if (m_Flags & OMPACC_ClearTarget)
         ClearBlock(&m_pTarget[iBlockOffs], pBlockEnd - pBlock);

      for (size_t iAcc = 0; iAcc < m_nThreads; ++ iAcc)
         Add2T(&m_pTarget[iBlockOffs], &pBlock[iAcc * m_nAlignedSize], FScalar(1), pBlockEnd-pBlock);
   }
}


}

#endif // CX_OPENMP_ACC_H
