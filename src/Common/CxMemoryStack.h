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

#ifndef _CX_MEMORYSTACK_H
#define _CX_MEMORYSTACK_H

#include <stddef.h> // for size_t
#include <string.h> // for memset

namespace ct {

#ifndef CX_DEFAULT_MEM_ALIGN
   // FMemoryStack will ensure that all allocations returned by Alloc() are
   // aligned to multiples of this many bytes. Must be power of 2.
   #define CX_DEFAULT_MEM_ALIGN 64
#endif // CX_DEFAULT_MEM_ALIGN


// note: boundary must be size-of-2
inline size_t AlignSizeT(size_t iPos, size_t Boundary) {
   size_t
      iNew = ((iPos - 1) | (Boundary - 1)) + 1;
   return iNew;
}


/// An object managing a continuous block of memory. Used for allocation-free
/// temporary work space storage.
/// Notes:
///    - Deallocations must occur in reverse allocation order. A deallocation
///      FREES ALL MEMORY ALLOCATED AFTER THE ALLOCATION POINT.
///    - Both alloc and free are effectively non-ops in terms of computational
///      cost. Any work space, no matter how large or small, can be allocated
///      at will.
///    - There is no memory fragmentation.
struct FMemoryStack
{
    /// allocate nSizeInBytes bytes off the stack. Return pointer to start.
    virtual void* Alloc(size_t nSizeInBytes) = 0;
    /// Free pointer p on stack. Should be an address returned by Alloc.
    virtual void Free(void *p) = 0;
    /// return amount of memory left (in bytes).
    virtual size_t MemoryLeft() = 0;
    /// returns whether p is identical to what a Alloc(0) would return if executed here.
    inline bool IsOnTop(void const *p) { return p == static_cast<void const*>(GetTop()); }

    /// allocate nObjects objects of size sizeof(T). Store pointer in p.
    /// Usage:
    ///     MyType *p;
    ///     Mem.Alloc(p, 10);  // allocate array of 10 objects of type MyType.
    template<class T>
    inline void Alloc(T *&p, size_t nObjects = 1){
        p = reinterpret_cast<T*>( this->Alloc(sizeof(T) * nObjects) );
        // note: objects are not constructed. Use only for PODs!
    }

    /// as Alloc(), but set allocated memory to zero.
    template<class T>
    inline void ClearAlloc(T *&p, size_t nObjects = 1){
        Alloc(p, nObjects);
        ::memset(p, 0, sizeof(T) * nObjects);
    }

    template<class T>
    inline T* AllocN(size_t nObjects, T const &){
        return reinterpret_cast<T*>( this->Alloc(sizeof(T) * nObjects) );
    }

    template<class T>
    inline T* ClearAllocN(size_t nObjects, T const &){
        T *p;
        this->Alloc(p, nObjects);
        ::memset(p, 0, sizeof(T) * nObjects);
        return p;
    }

    /// align stack to next 'Boundary' (must be power of two)-boundary .
    virtual void Align(unsigned Boundary);

    virtual ~FMemoryStack();

   // these are diagnosis functions used to track peak memory usage.
   char *GetPeak() { return m_pPeak; }
   char *GetBase() { return m_pInitialTop; }
   char *ResetPeak() { char *p = m_pPeak; m_pPeak = GetTop(); return p; };
   void MaximizePeak( char *p ) { if ( p > m_pPeak ) m_pPeak = p; };
protected:
   FMemoryStack();
   char
      *m_pInitialTop,  // top when stack was created
      *m_pPeak; // highest top reached since ResetPeak
private:
   virtual char *GetTop();
};

/// A memory stack defined by a base pointer and a size. Can own
/// the memory managed, but can also be put on memory allocated from
/// elsewhere.
struct FMemoryStack2 : public FMemoryStack {
   /// creates a stack of 'nInitialSize' bytes on the global heap.
   /// If nInitialSize is 0, the stack is created empty and a storage
   /// area can be later assigned by AssignMemory() or Create().
   explicit FMemoryStack2(size_t nInitialSize = 0) : m_pDataStart(0), m_pDataAligned(0) { Create(nInitialSize); }
   /// if pBase != 0 && nActualSize >= nNeededSize:
   ///    creates a stack of 'nActualSize' bytes at memory pBase (i.e., at a memory location given from the outside)
   /// otherwise
   ///    construction equivalent to FMemoryStack2((nNeededSize != 0)? nNeededSize : nActualSize).
   FMemoryStack2(char *pBase, size_t nActualSize, size_t nNeededSize = 0);
   inline void PushTopMark();
   inline void PopTopMark();

   void* Alloc(size_t nSizeInBytes); // override
   using FMemoryStack::Alloc;
   void Free(void *p); // override
   size_t MemoryLeft(); // override
   bool IsOnTop(void const *p); // override
   ~FMemoryStack2(); // override

   /// create a new stack of size 'nSize' bytes on the global heap.
   void Create( size_t nSize );
   /// create a new stack in the storage area marked by 'pBase',
   /// of 'nSize' bytes. This assumes that pBase is already allocated
   /// from elsewhere.
   void AssignMemory(char *pBase_, size_t nSize);
   void Destroy();
   void Align(unsigned Boundary); // override

private:
   size_t
      m_Pos, m_Size;
   char
      *m_pDataStart,
      *m_pDataAligned; // next default-aligned address beyond m_pDataStart.
   void DefaultAlign();
   bool m_bOwnMemory; // if true, we allocated m_pData and need to free it on Destroy().

   inline void UpdatePeak();
protected:
   virtual char *GetTop(); // override
private:
    FMemoryStack2(FMemoryStack2 const &); // not implemented
    void operator = (FMemoryStack2 const &); // not implemented
};

#ifdef MOLPRO
/// This class of FMemoryStack passes *every* allocation directly to
/// Molpro's icorr/corlsr routines! It is thus not thread safe and you cannot
/// have more than one of them!
struct FMemoryStackMolproCore : public FMemoryStack {
   /// creates a stack of 'nInitialSize' bytes on the global heap.
   /// If nInitialSize is 0, the stack is created empty and a storage
   /// area can be later assigned by AssignMemory() or Create().
   explicit FMemoryStackMolproCore();

   void* Alloc( size_t nSizeInBytes ); // override
   using FMemoryStack::Alloc;
   void Free( void *p ); // override
   size_t MemoryLeft(); // override
   ~FMemoryStackMolproCore(); // override


//    // these are diagnosis functions used to track peak memory usage.
//    char *GetPeak() { return m_pPeak; }
//    char *GetBase() { return m_pInitialTop; }
//    char *ResetPeak() { char *p = m_pPeak; m_pPeak = GetTop(); return p; };
//    void MaximizePeak( char *p ) { if ( p > m_pPeak ) m_pPeak = p; };
// private:
//    char *GetTop(){
//       // allowed to be slow, typically used only in assert()s.
//       return static_cast<char*>( Alloc(0) );
//    }

//    char
//       *m_pInitialTop,  // top when stack was created
//       *m_pPeak; // highest top reached since ResetPeak
private:
   FMemoryStackMolproCore(FMemoryStackMolproCore const &); // not implemented
   void operator = (FMemoryStackMolproCore const &); // not implemented
};
#endif


enum FMemoryAllocOptions {
   ALLOC_Clear = 0x001 // clear target of memory allocation (via a memset(..,0,..)).
};

// allocate a block of memory from pMem, and free it when the object goes out of scope.
template<class FType>
struct TMemoryLock
{
   FType
      *p; // pointer to allocated data

   // allocate nCount objects from pMem, and store the pointer in this->p
   inline TMemoryLock(size_t nCount, FMemoryStack *pMem, unsigned Flags = 0);
   // allocate nCount objects from pMem, and store the pointer in this->p and pOutsidePtr.
   inline TMemoryLock(FType *&pOutsidePtr, size_t nCount, FMemoryStack *pMem, unsigned Flags = 0);
   inline ~TMemoryLock();

   operator FType* () { return p; }
   operator FType const* () const { return p; }

   inline FType* data() { return p; }
   inline FType const* data() const { return p; }

   void Clear(size_t nObjects) { ::memset(p, 0, sizeof(FType) * nObjects); }

   FMemoryStack
      *m_pMemoryStack;
//    size_t
//       m_nObjects;
};

template<class FType>
TMemoryLock<FType>::TMemoryLock(size_t nCount, FMemoryStack *pMem, unsigned Flags)
{
   m_pMemoryStack = pMem;
   m_pMemoryStack->Alloc(p, nCount);
//    m_nObjects = nCount;
   if ((Flags & ALLOC_Clear) != 0)
      Clear(nCount);
}

template<class FType>
TMemoryLock<FType>::TMemoryLock(FType *&pOutsidePtr, size_t nCount, FMemoryStack *pMem, unsigned Flags)
{
   m_pMemoryStack = pMem;
   m_pMemoryStack->Alloc(p, nCount);
   pOutsidePtr = p;
//    m_nObjects = nCount;
   if ((Flags & ALLOC_Clear) != 0)
      Clear(nCount);
}

template<class FType>
TMemoryLock<FType>::~TMemoryLock()
{
   m_pMemoryStack->Free(p);
   p = 0;
}


/// split a base memory stack into sub-stacks, one for each openmp thread.
///
/// See comments on `TOmpAccBlock` (CxOpenMpAcc.h) for a usage example and
/// explanation of how this is meant to be used.
struct FMemoryStackArray
{
    explicit FMemoryStackArray(FMemoryStack &BaseStack);
    ~FMemoryStackArray();

    FMemoryStack2
        *pSubStacks;
    FMemoryStack
        *pBaseStack;
    unsigned
        nThreads;

    // return substack for the calling openmp thread.
    FMemoryStack2 &GetStackOfThread();

    // destroy this object (explicitly) already now
    void Release();

    char
        *pStackBase;
private:
    FMemoryStackArray(FMemoryStackArray const &); // not implemented
    void operator = (FMemoryStackArray const &); // not implemented
};


// Note: Moved in templatized version to CxOpenMpAcc.
// enum FOmpAccBlockFlags {
//    OMPACC_ClearTarget = 0x01
// };
//
// struct FOmpAccBlock
// {
//    double *pTls(); // return thread-local block of the storage. Must be called from within parallel block.
//    void Join(); // copy to target by summation; must be called OUTSIDE of parallel block.
//    // allocate and clear thread-local storage for all threads. Must be called OUTSIDE of parallel block.
//    void Init(double *pTarget, size_t nSize, unsigned Flags, FMemoryStack &Mem);
//
//    FOmpAccBlock() : m_pTarget(0), m_pTlsData(0), m_nSize(0), m_Flags(0) {};
//    FOmpAccBlock(double *pTarget_, size_t nSize_, unsigned Flags_, FMemoryStack &Mem_) : m_pTarget(0), m_pTlsData(0), m_nSize(0), m_Flags(Flags_) { Init(pTarget_, nSize_, Flags_, Mem_); }
// protected:
//    double
//       *m_pTarget,
//       *m_pTlsData; // m_nAlignedSize * m_nThreads
//    size_t
//       m_nSize,
//       m_nAlignedSize,
//       m_nThreads;
//    unsigned
//       m_Flags;
//    size_t nScalarBlockSize() const { return 128; }
//    size_t nTotalSize() const { return m_nAlignedSize * m_nThreads; }
// //    size_t nScalarBlocksTotal() const { return nTotalSize() / nScalarBlockSize(); }
//    size_t nScalarBlocks(size_t nSize) { return (nSize + nScalarBlockSize() - 1)/nScalarBlockSize(); }
// };


} // namespace ct


#endif // _CX_MEMORYSTACK_H
