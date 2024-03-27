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

#ifndef CX_PODARRAY_H
#define CX_PODARRAY_H

#include <new> // for placement new
#include <utility> // for std::forward
#include <functional> // for std::less
#include <string.h> // for memcpy
#include <stdlib.h> // for malloc/free
#include <algorithm> // for swap
#include <initializer_list>
//#ifndef assert
   #include "CxDefs.h"
//#endif

namespace ct {

// swap two primitive values. Here such that we need not include <algorithm> here.
template<class FType>
void swap1(FType &A, FType &B){
   FType t = A; A = B; B = t;
}


#ifndef CX_DEFAULT_ARRAY_STATIC_STORAGE
   #define CX_DEFAULT_ARRAY_STATIC_STORAGE 0
#endif // CX_DEFAULT_ARRAY_STATIC_STORAGE


// don't want to include <algorithm> here... it drags in so much other stuff.
template<class FType>
FType _s_max(FType a, FType b) {
   if (a < b)
      return b;
   else
      return a;
}

namespace detail {
   // workaround for this:
   //
   //    CxPodArray.h:285:7: warning: ISO C++ forbids zero-size array [-Wpedantic]
   //
   // VC apparently does not like this at all, either.
   // So instead of adding the zero-length array, we wrap the static storage
   // into a separate template class (TStaticStorage), and then do a partial
   // specialization for the case StaticStorageSize == 0 does not have the zero-
   // size array member at all. I think this one should be legal C++ without any
   // caveats.

   template<class FType, size_t StaticStorageSize>
   struct TStaticStorage {
      FType m_StaticData[StaticStorageSize];
      inline FType &operator [] (size_t i) { return m_StaticData[i]; }
      inline FType const &operator [] (size_t i) const { return m_StaticData[i]; }
      inline FType *data() { return &m_StaticData[0]; }
      inline FType const *data() const { return &m_StaticData[0]; }
      size_t size() const { return StaticStorageSize; }
   };

   template<class FType>
   struct TStaticStorage<FType, size_t(0)> {
      FType &operator [] (size_t i); // not implemented.
      FType const &operator [] (size_t i) const; // not implemented
      inline FType *data() { return 0; }
      inline FType const *data() const { return 0; }
      size_t size() const { return 0; }
   };
}

// raises std::runtime_error
DECL_NORETURN void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount, size_t PrevReservedCount);
double _fSizeMb(size_t nBytes);


/// A dynamic array of POD (``plain old data'') types that can be copied via a
/// memcpy and cleared with memset(...,0,..), and require neither constructors
/// no destructors to work.
///
/// - Main point for this is that std::vector does not allow containing C-style
///   arrays (e.g., double [3]), because C-style arrays are not assignable.
/// - Additionally, std::vector can be *very* slow when allocating large amounts
///   of data, because the data is set to zero on resize. This class explicitly
///   DOES NOT DO THAT: It has RAII semantics, but non-explicitly touched data
///   is just random.
/// - Also, this class neither calls constructors nor destructors for its data
///   elements, so this should most definitely NOT be used on aggregates which
///   would be broken due to this.
///
/// Template arguments:
/// - StaticStorageSize: If 'StaticStorageSize' is zero, the class is
///   effectively a 'buffer-ptr + size' pair. Otherwise, the class contains
///   a local buffer for 'StaticStorageSize' elements which will be used
///   for allocation-less storage of the first few elements until a dynamic
///   allocation is required when the static storage size is exceeded.
///
/// TODO:
/// - The class has become complex enough to warrant splitting into a dumb
///   data implementation class (with char* ptrs) and an actual typed
///   class, I think. But it is a bit messy and dangerous, due to the
///   strict aliasing rules. Would need some actual time to think through.
///   Reconsider if I get really bored...
template<class FType, size_t StaticStorageSize = CX_DEFAULT_ARRAY_STATIC_STORAGE>
struct TArray
{
   typedef FType *iterator;
   typedef FType const *const_iterator;
   typedef ::size_t size_type;
   typedef FType value_type;

   iterator begin() { return m_pData; }
   iterator end() { return m_pDataEnd; }
   const_iterator begin() const { return m_pData; }
   const_iterator end() const { return m_pDataEnd; }

   FType &front() { return *m_pData; }
   FType &back() { return *(m_pDataEnd-1); }
   FType const &front() const { return *m_pData; }
   FType const &back() const { return *(m_pDataEnd-1); }

   FType *data() { return m_pData; }
   FType const *data() const { return m_pData; }

   FType &operator[] (size_type i) { return m_pData[i]; }
   FType const &operator[] (size_type i) const  { return m_pData[i]; }

   size_type size() const { return m_pDataEnd - m_pData; }
   size_type capacity() const { return _s_max(size_type(m_pReservedEnd - m_pData), size_type(StaticStorageSize)); }
   size_type max_size() const { return size_type(-1)/sizeof(FType); }
   bool empty() const { return m_pData == m_pDataEnd; }

   // WARNING: contrary to std::vector, this function *DOES NOT*
   // initialize (or touch, for that matter) the newly created data.
   // If you want that behavior, use the other resize function!
   // It also does not call any destructors in case the number of
   // element decreases.
   void resize(size_type n) {
      reserve(n);
      m_pDataEnd = m_pData + n;
   }


   void resize(size_type n, FType t) {
      size_type old_size = size();
      resize(n);
      if (old_size < n) {
         for (size_type i = old_size; i < n; ++ i)
            m_pData[i] = t;
      }
   }

   // memset the entire array data to 0.
   void clear_data() {
      memset(m_pData, 0, sizeof(FType)*size());
   }

   void resize_and_clear(size_type n) {
      resize(n);
      clear_data();
   }


   // change the number of controlled elements to zero, but DO NOT
   // release the array storage space they used
   void clear() {
      resize(0);
   }

   // clear array data and release memory.
   //
   // Note: This frees all controlled data, but does not invalidate the array
   // object itself. Calling release() and then refilling an array is fine.
   void release() {
      _ReleaseBuffer();
   }

   void push_back(FType const &t) {
      _ReserveSpaceForPushBack();
      m_pDataEnd[0] = t;
      ++m_pDataEnd;
   }

#ifdef CX_ARRAY_HAVE_CXX11
   template<class... Args>
   void emplace_back(Args... args) {
      _ReserveSpaceForPushBack();
      new(m_pDataEnd) FType(std::forward(args)...);
      ++m_pDataEnd;
   }
#endif

   void pop_back() {
      assert(!empty());
      m_pDataEnd -= 1;
   }


   void reserve(size_type NewReservedSize) {
      if (NewReservedSize == 0)
         return;
      if (StaticStorageSize != 0 && NewReservedSize <= StaticStorageSize) {
         // is this the initial reservation? I.e., no data in *this so far?
         if (m_pData == 0) {
            // assign reserved space to local object static buffer.
            // This will automatically deal with the reserved size. In this case,
            // the branch below will not be taken if the reservation size happens
            // to be small enough.
            // Note that the actual "reservation" in this case is always for the
            // full static storage size, independently of the requested target
            // size (because the memory is already there anyway)
            m_pData = m_StaticStorage.data();
            m_pDataEnd = m_pData;
            m_pReservedEnd = m_pData + StaticStorageSize;
            m_DataIsOnStaticBuffer = true;
         }
      }
      if (_nReserved() < NewReservedSize) {
         size_type
            nElem = size();
         FType
            *pNewData;
         if (nElem > 0 && !m_DataIsOnStaticBuffer && capacity() <= (nElem + (3+nElem)/4)) {
            // attempt a reallocation of the present memory block if...
            // ...there are already data elements stored,
            // ...current data is on dynamic memory,
            // ...and the currently reserved memory block size does not
            //    significantly exceed the memory size of the actual data.
            // The last condition is for the case in which the reallocation has to
            // resort to an alloc & copy; in that case we'd be better off memcpy'ing
            // the actual amount of amount only, rather than risking realloc copying
            // a large amount of undefined empty space betweenm_pDataEnd and m_pReservedEnd.
            pNewData = static_cast<FType*>(::realloc(m_pData, sizeof(FType) * NewReservedSize));
         } else {
            pNewData = static_cast<FType*>(::malloc(sizeof(FType) * NewReservedSize));
            if (pNewData != 0) {
               if (nElem != 0)
                  ::memcpy(pNewData, m_pData, sizeof(FType) * nElem);
               _ReleaseBuffer();
            }
         }
         if (pNewData == 0)
            return _RaiseMemoryAllocError(__func__, sizeof(value_type), NewReservedSize, _nReserved());
         m_DataIsOnStaticBuffer = false;
         m_pData = pNewData;
         m_pDataEnd = m_pData + nElem;
         m_pReservedEnd = m_pData + NewReservedSize;
      }
   }


   void shrink_to_fit(size_type NewReservedSize) {
      if (m_DataIsOnStaticBuffer)
         return; // static buffer cannot be shrunken.
      if (NewReservedSize < size())
         // don't reduce capacity below actual size needed for elements
         NewReservedSize = size();
      if (NewReservedSize == 0)
         return;
      size_t nElem = size();
      // fulfill request (via reallocation, if possible) only if either filling level is < 75% of capacity
      // or the shrink_to_fit would free up at least 16 kiB
      if (capacity() >= (nElem + (3+nElem)/4) || (capacity() - nElem) * sizeof(FType) >= (size_t(16) << 10)) {
         assert(!m_DataIsOnStaticBuffer);
         FType
            *pNewData = static_cast<FType*>(::realloc(m_pData, sizeof(FType) * NewReservedSize));
         if (pNewData == 0)
            return _RaiseMemoryAllocError(__func__, sizeof(value_type), NewReservedSize, _nReserved());
         m_pData = pNewData;
         m_pDataEnd = m_pData + nElem;
         m_pReservedEnd = m_pData + NewReservedSize;
      }
   }


   TArray()  {
      _InitDataMembers();
   }

   TArray(TArray const &other) {
      _InitDataMembers();
      *this = other;
   }

#ifdef CX_ARRAY_HAVE_CXX11
   TArray(TArray &&other) {
      _InitDataMembers();
      _MoveArrayData(other);
   }
#endif // CX_ARRAY_HAVE_CXX11

   // WARNING: array's content not initialized with this function!
   // (with intention!)
   explicit TArray(size_type n) {
      _InitDataMembers();
      resize(n);
   }

   // initializes the the array to size `n`, with all `n` entries set to `t`.
   TArray(size_type n, FType t) {
      _InitDataMembers();
      assign(n, t);
   }

#ifdef CX_ARRAY_HAVE_CXX11
   TArray(std::initializer_list<FType> init) {
      _InitDataMembers();
      assign(init.begin(), init.end());
   }
#endif // CX_ARRAY_HAVE_CXX11

   ~TArray() {
      _ReleaseBuffer();
   }

   TArray const &operator = (TArray const &other) {
      if (this == &other) return *this;
      resize(other.size());
      memcpy(m_pData, other.m_pData, sizeof(FType) * size());
      return *this;
   }

#ifdef CX_ARRAY_HAVE_CXX11
   TArray &operator = (TArray &&other) {
      _InitDataMembers();
      _MoveArrayData(other);
      return *this;
   }
#endif // CX_ARRAY_HAVE_CXX11

   void swap(TArray &other) {
      if (StaticStorageSize == 0 || (!other.m_DataIsOnStaticBuffer && !this->m_DataIsOnStaticBuffer)) {
         // no data on static buffers. Some pointer swappery is sufficient.
         assert(!other.m_DataIsOnStaticBuffer && !this->m_DataIsOnStaticBuffer);
         swap1(m_pData, other.m_pData);
         swap1(m_pDataEnd, other.m_pDataEnd);
         swap1(m_pReservedEnd, other.m_pReservedEnd);
      } else {
         // at least one of the array data sets is on static buffers. This makes
         // this complicated. We actually need to physically copy stuff around
         // in this case.
         TArray tmp(*this);
         *this = other;
         other = tmp;
      }
   }


   template<class FRandomIt>
   void assign(FRandomIt begin, FRandomIt end) {
      resize(end - begin);
      for (size_type i = 0; i < size(); ++ i)
         m_pData[i] = FType(begin[i]);
   }

   // replaces content of *this with `n` copies of value `t`.
   void assign(size_type n, FType t) {
      resize(n);
      for (size_type i = 0; i < n; ++ i)
         m_pData[i] = t;
   }


#ifdef CX_ARRAY_HAVE_CXX11
   void assign(std::initializer_list<FType> init) {
      assign(init.begin(), init.end());
   }
#endif // CX_ARRAY_HAVE_CXX11


   template<class FRandomIt>
   void insert(iterator pos, FRandomIt begin, FRandomIt end) {
      // WARNING: not checked.
      size_type num = end - begin, ipos = pos - this->begin();
      resize(num + this->size());
      for (size_type i = this->size(); i != ipos+num; -- i)
         this->m_pData[i-1] = this->m_pData[i-1-num];
      for (size_type i = 0; i < num; ++ i)
         this->m_pData[ipos+i] = FType(begin[i]);
   }

   // first compare array size; then lexicographical comparison of array elements for
   // same-size arrays. Employs member type's operator <.
   bool operator < (TArray const &other) const {
      if (this->size() < other.size()) return true;
      if (other.size() < this->size()) return false;
      assert(this->size() == other.size());
      for (size_t i = 0; i != size(); ++ i) {
         FType const
            &ai = this->m_pData[i],
            &bi = other.m_pData[i];
         if (ai < bi) return true;
         if (bi < ai) return false;
      }
      return false;
   }

   bool operator != (TArray const &other) const {
      if (other.size() != this->size())
         return true;
      for (size_t i = 0; i != size(); ++ i) {
         if (this->m_pData[i] != other.m_pData[i])
            return true;
      }
      return false;
   }

   bool operator == (TArray const &other) const {
      return !(*this != other);
   }
private:
   FType
      // start of controlled array. may be 0 if no data is contained.
      *m_pData,
      // end of actual data
      *m_pDataEnd,
      // end of reserved space (i.e., of the space *this owns)
      *m_pReservedEnd;
   bool
      // if true, m_pData points to a dynamically allocated piece of memory
      // which needs to be freed again. If false, the pointers in *this point
      // into m_StaticStorage.
      m_DataIsOnStaticBuffer;
   detail::TStaticStorage<FType, StaticStorageSize>
      // if the total number of elements is smaller than StaticStorageSize, then
      // this buffer may be used to store them.
      m_StaticStorage;
   size_t _nReserved() const {
      return static_cast<size_type>(m_pReservedEnd - m_pData);
   }

   double _fSizeMb(size_t n) {
      return _fSizeMb(sizeof(FType) * n);
   }

   void _ReleaseBuffer() {
      if ((StaticStorageSize == 0) || !m_DataIsOnStaticBuffer)
         ::free(m_pData);
      m_pData = 0;
      m_pDataEnd = 0;
      m_pReservedEnd = 0;
      // there is no data, so it's not on a static buffer. This handling simplifies
      // some operations.
      m_DataIsOnStaticBuffer = false;
   }

   void _InitDataMembers() {
      m_pData = 0;
      m_pDataEnd = 0;
      m_pReservedEnd = 0;
      m_DataIsOnStaticBuffer = false;
   }

   void _MoveArrayData(TArray &other) {
      if (StaticStorageSize == 0 || !other.m_DataIsOnStaticBuffer) {
         // other array's data is on dynamic memory, not on static buffers. just
         // move around the pointers.
         if (m_pDataEnd != 0)
            _ReleaseBuffer();

         m_pData = other.m_pData;
         m_pDataEnd = other.m_pDataEnd;
         m_pReservedEnd = other.m_pReservedEnd;
         other._InitDataMembers();
      } else {
         // other array's data is on its static memory buffer---which we obviously
         // cannot take control of in here. So there is no way to avoid physically copying
         // that other array's data.
         assign(other.begin(), other.end());
         other.clear();
      }
   }


   // exponentially increases array size if space is not enough, for average logarithmic
   // time complexity of push_back/emplace_back
   inline void _ReserveSpaceForPushBack() {
      if (this->size() + 1 > static_cast<size_type>(m_pReservedEnd - m_pData)) {
         reserve(2 * this->size() + 1);
      }
   }
};


// Note to self: @swap:
// --------------------
// The recommended usage patter for a free template swap<T> function is this:
//
// 1. When declaring swappable classes `T`, make sure that a free-standing
//    swap(T&, T&) function is declared *in the same namespace*. This way in
//    a call of `swap(a,b)` where `a` or `b` are of type `T`, this function
//    will automatically be considered as a candidate (even without explicit
//    namespace qualification) by means of argument-dependent lookup (ADL).
//
// 2. In functions using swap on unknown/template types, *do not* use qualified
//    namespaces to refer to the swap function (e.g., "std::swap"), but rather
//    bring in viable candidates into the local namespace:
//
//       template<class T>
//       void some_function(T &a, T &b) { ...
//           using std::swap;
//           swap(A, B);
//       }
//
// 3. The actual swap() resolved in there then may or may not be std::swap,
//    depending on actual #includes and, most importantly, the concrete type
//    of `T` and functions handling it.
//      In effect, the "using std::swap;" just enables *considering* all the std
//    namespace specializations as candidates in the function resolution
//    ---however, if a `swap` is defined in the same namespace as `T`, it will
//    *also* be considered as candidate (via ADL), even if `T`s namespace has
//    not been explicitly imported.
template<class FType>
void swap(TArray<FType> &A, TArray<FType> &B) {
   A.swap(B);
}


} // namespace ct



namespace ct {

/// Array class with fixed maximum size which stores its elements
/// in-place (i.e., no allocations).
///
/// For technical reasons, all MaxN array elements are default-constructed on
/// construction and only destroyed when TArrayFix is (i.e., technically,
/// all array elements are alive all the time even if size() < max_size()).
template<class FType, unsigned int MaxN, class FSizeType = std::size_t>
struct TArrayFix
{
    typedef FType value_type;
    typedef FType *iterator;
    typedef FType const *const_iterator;
    typedef FSizeType size_type;
    typedef size_type size_t;

    // compiler-generated default destructor, copy-ctor and assignment op should work.
    TArrayFix()
        : nSize(0)
    {};

    explicit TArrayFix( size_t nEntries )
        : nSize(0)
    { resize(nEntries); };

    TArrayFix( size_t nEntries, FType const &Scalar )
        : nSize(0)
    { resize(nEntries); *this = Scalar; };


    template<class FInputIt>
    TArrayFix( FInputIt first, FInputIt last ){
        nSize = 0;
        while( first != last )
            push_back(*(first++));
    }

    FType &operator[] ( size_t i ){ assert(i < nSize); return m[i]; };
    FType const &operator[] ( size_t i ) const {  assert(i < nSize); return m[i]; };

    FType &back(){ assert(nSize!=0); return m[nSize-1]; };
    FType const &back() const { assert(nSize!=0); return m[nSize-1]; };
    FType &front(){ assert(nSize!=0); return m[0]; };
    FType const &front() const { assert(nSize!=0); return m[0]; };

    bool operator == ( TArrayFix const &other ) const {
        if ( size() != other.size() )
            return false;
        for ( size_t i = 0; i != size(); ++ i )
            if ( (*this)[i] != other[i] )
                return false;
        return true;
    }

    bool operator != ( TArrayFix const &other ) const {
        return !this->operator ==(other);
    }

    size_t size() const { return nSize; }
    bool empty() const { return nSize == 0; }
    size_t capacity() const { return MaxN; }
    size_t max_size() const { return MaxN; }
    void clear() { resize(0); };

    void push_back( FType const &t ){
        assert( nSize < MaxN );
        m[nSize] = t;
        ++nSize;
    }
    void pop_back(){
        assert( nSize > 0 );
        --nSize;
    }
    void resize( size_t NewSize ){
        assert( NewSize <= MaxN );
        nSize = NewSize;
    }
    void resize( size_t NewSize, FType const &value ){
        assert( NewSize <= MaxN );
        for ( size_t i = nSize; i < NewSize; ++ i )
            m[i] = value;
        nSize = NewSize;
    }

    iterator erase( iterator itFirst, iterator itLast ){
        assert( itFirst >= begin() && itLast <= end() && itFirst <= itLast );
        nSize -= itLast - itFirst;
        for ( iterator it = itFirst; it < end(); ++ it, ++ itLast )
            *it = *itLast;
        return itFirst;
    };

    iterator erase( iterator itWhere ){
        return erase(itWhere, itWhere+1);
    };

    FType *data() { return &m[0]; };
    FType const *data() const { return &m[0]; };

    // assign Scalar to every element of *this
    void operator = (FType const &Scalar){
        for ( iterator it = begin(); it != end(); ++ it )
            *it = Scalar;
    }

    template <class FIt>
    void assign(FIt first, FIt last) {
        resize(last - first);
        for ( FSizeType i = 0; i < nSize; ++ i )
            m[i] = first[i];
    }

    iterator begin(){ return &m[0]; }
    iterator end(){ return &m[nSize]; }
    const_iterator begin() const { return &m[0]; }
    const_iterator end() const { return &m[nSize]; }

protected:
    FType
        m[MaxN];
    FSizeType
        nSize; // actual number of elements (<= MaxN).
};

// lexicographically compare two arrays. -1: A < B; 0: A == B; +1: A > B.
template<class FType, unsigned int MaxN, class FPred>
int Compare(TArrayFix<FType,MaxN> const &A, TArrayFix<FType,MaxN> const &B, FPred less = std::less<FType>())
{
    typedef TArrayFix<FType,MaxN>
        FArray;
    if (A.size() < B.size()) return -1;
    if (B.size() < A.size()) return +1;
    for (typename FArray::size_type i = 0; i < A.size(); ++ i) {
        if (less(A[i], B[i])) return -1;
        if (less(B[i], A[i])) return +1;
    }
    return 0;
}

} // namespace ct


#endif // CX_PODARRAY_H
