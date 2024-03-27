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

/// @file CxIntrusivePtr.h
///
/// A re-implementation of the boost intrusive ptr. reason is that due to
/// recent changes in gcc, some problematic template order- of-definition
/// rules are now enforced, which leads to the situation that
/// add_ref/release *must* be declared before intrusive_ptr.hpp is
/// included... something we can of course not guarantee in general!
///
/// Additionally, this allows us to drop one more dependency on boost.
///
/// Actual pointer relies on free-standing:
///      void intrusive_ptr_add_ref(T *p);
///      void intrusive_ptr_release(T *p);
/// functions, just as the boost ptr.
///
/// Objects derived from class FIntrusivePtrDest1 support those methods
/// automatically, but deriving from this class is not mandatory.
/// (as long as the corresponding add_ref/release functions are declared
/// before this header is included).
///
/// FIXME: Update this version with the updates and documentation
/// added for MqIntrusivePtr.hpp

#ifndef CX_INTRUSIVE_PTR
#define CX_INTRUSIVE_PTR

#include <algorithm> // for std::swap.
#include <functional> // for std::less

namespace ct {
   struct FIntrusivePtrDest1;
}
void intrusive_ptr_add_ref(ct::FIntrusivePtrDest1 const *p);
void intrusive_ptr_release(ct::FIntrusivePtrDest1 const *p);

namespace ct {

/// Objects derived from this class automatically support attachment
/// to both cx's ct::TIntrusivePtr and boost's intrusive_ptr.
///
/// Notes:
/// - Used as a base class, this class adds a virtual destructor and a
///   m_RefCount member variable to a class, and initializes it to zero. Apart
///   from that, the class itself does not do anything.
///
/// - There are global intrusive_ptr_add_ref and intrusive_ptr_release functions
///   declared which can access FIntrusivePtrDest1* objects.
///
/// - In case you do *not* wish to use whatever is derived from this with
///   intrusive ptrs: that is fine, too!
///   The functionality added here does not interfere with other types of uses
///   of virtual base classes. E.g., attaching FIntrusivePtrDest1-derived
///   classes to shared_ptr is totally fine (as long as the types of raw and
///   smart ptrs are not mixed).
struct FIntrusivePtrDest1
{
   FIntrusivePtrDest1() : m_RefCount(0) {}
   inline virtual ~FIntrusivePtrDest1() = 0;
   // meant as a function for debugging/info purposes only.
   ptrdiff_t GetRefCount() const { return m_RefCount; }

   // making copies of the targeted objects does not affect how many references
   // are kept to them. E.g., in this:
   //
   //   pOne = FIntrusivePtr(new ...)
   //   FObj Two = *pOne;
   //
   // there is no reference to object `Two`, regardless how many there were in
   // `pOne`. Similarly, assignment the other way around, i.e., *pOne = Two,
   // would not change the references in either `*pOne` or `Two`. So just like
   // in the regular construtor, in copy construction we need to set the ref
   // count to zero, and in assignment, we need to leave the ref count alone at
   // whatever it is before making the copy.
   FIntrusivePtrDest1(FIntrusivePtrDest1 const &/*other*/) : m_RefCount(0) {}
   void operator = (FIntrusivePtrDest1 const &/*other*/) {} // <-- must leave ref count alone on copy.
private:
   ptrdiff_t mutable
      m_RefCount;
   friend void ::intrusive_ptr_add_ref(FIntrusivePtrDest1 const *p);
   friend void ::intrusive_ptr_release(FIntrusivePtrDest1 const *p);
};

inline FIntrusivePtrDest1::~FIntrusivePtrDest1()
{
}

} // namespace ct

inline void intrusive_ptr_add_ref(ct::FIntrusivePtrDest1 const *p) {
   p->m_RefCount += 1;
}

inline void intrusive_ptr_release(ct::FIntrusivePtrDest1 const *p) {
#ifdef assert
   assert(p->m_RefCount != 0);
#endif
   p->m_RefCount -= 1;
   if (p->m_RefCount == 0)
      delete p;
}

namespace ct {

template <class T>
struct TIntrusivePtr
{
   typedef T
      *pointer;
   typedef T
      element_type;
   typedef TIntrusivePtr
      this_type;
   // we do it just like boost: this is a pointer-to-member, which can be
   // implicitly converted to bool, but not to, say, other integral types or
   // floating point numbers. returning this instead of actual bools thus prevents
   // some unwanted implicit conversions, e.g.:
   //      TIntrusivePtr<something> p;
   //      double f = p;
   // is not legal with 'pointer to member', but would be legal with 'bool'.
   typedef pointer
      this_type::*unspecified_bool_type;

//
// constructors, destructors, and assignment operators
//
   TIntrusivePtr() : m_ptr(0) {}

   // construct from raw pointer of same type
   TIntrusivePtr(pointer p, bool add_ref = true) : m_ptr(p)  {
      if (m_ptr != 0 && add_ref)
         intrusive_ptr_add_ref(m_ptr);
   }

   TIntrusivePtr(TIntrusivePtr const &p) : m_ptr(p.get())  {
      if (m_ptr != 0)
         intrusive_ptr_add_ref(m_ptr);
   }

   // construct from related intrusive pointer
   template <class U>
   TIntrusivePtr(TIntrusivePtr<U> const &p) : m_ptr(p.get())  {
      if (m_ptr != 0)
         intrusive_ptr_add_ref(m_ptr);
   }

   ~TIntrusivePtr() {
      if (m_ptr != 0)
         intrusive_ptr_release(m_ptr);
   }

   // assign via copy constructor: x = bla -> x.swap(type(bla))
   TIntrusivePtr &operator = (pointer other) {
      this_type(other).swap(*this);
      return *this;
   };

   TIntrusivePtr &operator = (TIntrusivePtr const &other) {
      this_type(other).swap(*this);
      return *this;
   }

   template <class U>
   TIntrusivePtr &operator = (TIntrusivePtr<U> const &other) {
      this_type(other).swap(*this);
      return *this;
   }

//
// accessors
//
   pointer &get() { return m_ptr; }
   pointer const &get() const { return m_ptr; }

   T &operator*() const {  return *m_ptr; }

   pointer &operator->() { return m_ptr; }
   pointer const &operator->() const { return m_ptr; }

   operator unspecified_bool_type () const {
      return m_ptr == 0? 0: &this_type::m_ptr;
   }

   bool operator !() const { return m_ptr == 0; }

   void swap(TIntrusivePtr &other) {
      std::swap(m_ptr, other.m_ptr);
   }
   void swap(TIntrusivePtr &&other) {
      std::swap(m_ptr, other.m_ptr);
   }

   // replace the currently managed object, if present, by pObj
   void reset(T *pObj) {
      this->swap(TIntrusivePtr(pObj));
   }
private:
   T *m_ptr;
};

template<class U, class V>
inline void swap(TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   a.swap(b);
}



// comparison operators: ==, !=, and < between the intrusive pointers
// and raw pointers (in all permutations).
template<class U, class V>
inline bool operator == (TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   return a.get() == b.get();
}

template<class U, class V>
inline bool operator != (TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   return a.get() != b.get();
}

template<class U, class V>
inline bool operator < (TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   return a.get() < b.get();
}

template<class U, class V>
inline bool operator == (TIntrusivePtr<U> const &a, V const *b) {
   return a.get() == b;
}

template<class U, class V>
inline bool operator != (TIntrusivePtr<U> const &a, V const *b) {
   return a.get() != b;
}
template<class U, class V>
inline bool operator < (TIntrusivePtr<U> const &a, V const *b) {
   return a.get() < b;
}

template<class U, class V>
inline bool operator == (U const *a, TIntrusivePtr<V> const &b) {
   return a == b.get();
}

template<class U, class V>
inline bool operator != (U const *a, TIntrusivePtr<V> const &b) {
   return a != b.get();
}

template<class U, class V>
inline bool operator < (U const *a, TIntrusivePtr<V> const &b) {
   return a < b.get();
}



/// An auxiliary class which allows using either smart pointer or raw pointer
/// types as map/set keys, by acting as a comparison predicate which forwards
/// comparisons to the pointed-to objects (as opposed to comparing the pointers
/// themselves, which is what the base `operator <`'s do).
///
/// Example:
///    typedef TIntrusivePtr<FObject>
///       FObjectPtr;
///    typedef std::set<FObjectPtr, TComparePtrTarget<FObject> >
///       FObjectSet;
/// ...FObjectSet may now be used to collect unique FObject instances
/// (provided FObject has a valid `operator <`)
template <class T, class FCompare = std::less<T> >
struct TComparePtrTarget
{
   TComparePtrTarget() {};
   TComparePtrTarget(FCompare Pred_) : m_Pred(Pred_) {};

//    bool operator () (TIntrusivePtr<T> const &A, TIntrusivePtr<T> const &B) const {
   template<class FPtrType>
   bool operator () (FPtrType const &pA, FPtrType const &pB) const {
      return m_Pred(*pA, *pB);
   }
protected:
   FCompare
      m_Pred;
};


} // namespace ct

#endif // CX_INTRUSIVE_PTR
