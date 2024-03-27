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

#ifndef CX_ITER_TOOLS_H
#define CX_ITER_TOOLS_H

#include <tuple>
#include <utility>
#include <functional> // for reference wrapper (more precisely, the associated unwrapping with reduced decay)
#include <type_traits>

namespace ct {
namespace iter_detail {
   // Type of the actual object's you'll get as variables in range-based for loops
   // with the enumerate() below. E.g., in loops such as
   //
   //   for (auto &&it : enumerate(my_vector)) { .. }
   //   for (auto const &it : enumerate(my_vector)) { .. }
   //
   // the base type of `it` will be one of those.
   //
   // Note: this auxiliary class is always instanciated with *reference types* for
   // FIndex and FValue! The referenced index remains in the parent iterator
   // object, and the value in the enumerated parent sequence. `i` and `v` here are
   // not copies.
   template<class FIndex, class FValue>
   struct TIndexAndValue {
      typedef FIndex index_type;
      typedef FValue value_type;
      index_type i;
      value_type v;

      index_type const &index() const { return i; }
      value_type &value() { return v; }
      value_type const &value() const { return v; }

      // Add some iterator-like dereferencing operators to get to the base value.
      // May look a bit odd at first, but I don't think it's so bad (compared with
      // the alternatives, specifically)
      typedef decltype(&v) pointer_type; // <-- can't just do value_type *, becaue value_type is instanciated as a reference type.
      value_type &operator*() { return v; }
      value_type const &operator*() const { return v; }
      pointer_type operator->() { return &v; }
      pointer_type const operator->() const { return &v; }

   #if __cplusplus >= 201703L
      // define e.get<n>() and below tuple_size<> and tuple_element<> to provide
      // support for 'structured binding'. There probably is a better way of doing
      // it (soo ugly), but I experimented a bit and did not find it. My new-age
      // C++ is not very strong.
      template<std::size_t iElemT>
      std::tuple_element_t<iElemT, TIndexAndValue> &get() {
         if constexpr (iElemT == 0) return i;
         if constexpr (iElemT == 1) return v;
      }
      template<std::size_t iElemT>
      std::tuple_element_t<iElemT, TIndexAndValue> const &get() const {
         if constexpr (iElemT == 0) return i;
         if constexpr (iElemT == 1) return v;
      }
   #endif
   };

} // namespace iter_detail
} // namespace ct


// #if __cplusplus >= 201703L
// ^-- we can actually always define that in >= C++11. It just won't work
//     for structured bindings before C++17... but there are a few other
//     uses for std::tuple_size & co.
namespace std {
   template<class I, class V>
   struct tuple_size<ct::iter_detail::TIndexAndValue<I,V> >{ enum { value = 2 }; };

   template<class I, class V>
   struct tuple_element<0, ct::iter_detail::TIndexAndValue<I,V> >{ typedef I type; };

   template<class I, class V>
   struct tuple_element<1, ct::iter_detail::TIndexAndValue<I,V> >{ typedef V type; };
}
// #endif


namespace ct {

namespace iter_detail {

   // moved out of the enumerate function because we need it for the
   // return type of enumerate() in C++11.
   template<class FSequence>
   struct TSequenceWrapper {
      struct iterator {
         typedef decltype(std::begin(std::declval<FSequence>()))
            FSequenceIt;
         typedef std::size_t
            index_type;
         typedef decltype(*std::declval<FSequenceIt>())
            value_type;
         index_type
            i; // element index relative to sequence.begin()
         FSequenceIt
            p; // actual iterator to element
         bool operator != (iterator const &other) const { return this->p != other.p; }
         void operator ++ () { ++ i; ++ p; }

         typedef ct::iter_detail::TIndexAndValue<index_type const&, value_type&>
            iav_proxy_type; // iav: index and value
         iav_proxy_type operator * () const { return iav_proxy_type{i, *p}; }
      };
      FSequence m_sequence;
      iterator begin() { return iterator{0, std::begin(m_sequence)}; }
      iterator end() { return iterator{0, std::end(m_sequence)}; }
   };

} // namespace iter_detail

// python "enumerate"-like wrapper, which adds running indices to iterable sequences.
// Base is recycled from here: https://www.reedbeta.com/blog/python-like-enumerate-in-cpp17/
// I made several adjustments, though (and I don't claim they are necessarily for the better)
//
// see ~/dev/cx/enumerate_test/enum_test_cx_include_version.cpp for tests.
// Looked like all the core things work out (no copy of container unless needed,
// and if done, it's using move, the elements do not appear to be copied in any
// variant I tested, and it works in both C++11 and C++17 now (structured
// bindings only in the latter of course)). Still might need some more in-depth
// examination before I'd swear on it (references to temporaries etc---not
// always easy to find, although a initial sanitizing run didn't find anything)
template<class FSequence>
constexpr iter_detail::TSequenceWrapper<FSequence> enumerate(FSequence &&sequence)
{
   return iter_detail::TSequenceWrapper<FSequence>{std::forward<FSequence>(sequence)};
}


#if 0
// The base version. Core is recycled from here:
// https://www.reedbeta.com/blog/python-like-enumerate-in-cpp17/
// I made several adjustments, but mostly cosmetical, and I don't claim they are necessarily for the better.
//
// This works fine as-is, in principle. The only thing irritating is that it
// returns a std::tuple (via tie), which is rather ugly because there is no
// sensible way to access its elements (the only way to get them is
// std::get<n>(..)). With structured binding that is a non-issue, but I do not
// wish to enforce C++17 support for the time being. (Some of my software should
// run on HPC clusters, which for some reason often have 10 year old compilers)
template<class FSequence, class FSequenceIt = decltype(std::begin(std::declval<FSequence>()))>
constexpr auto enumerate(FSequence &&sequence)
{
   struct iterator {
      std::size_t
         i; // element index relative to sequence.begin()
      FSequenceIt
         p; // actual iterator to element
      bool operator != (iterator const &other) const { return this->p != other.p; }
      void operator ++ () { ++ i; ++ p; }
      auto operator * () const { return std::tie(i, *p); }
   };
   struct FSequenceWrapper {
      FSequence m_sequence;
      auto begin() { return iterator{0, std::begin(m_sequence)}; }
      auto end() { return iterator{0, std::end(m_sequence)}; }
   };
   return FSequenceWrapper{std::forward<FSequence>(sequence)};
}
#endif

} // namespace ct


/// Defines a minimal iterator interface to allow use in range-based for loops
/// and generic container printing mechanism. This one is for objects which have
/// a working [] operator to access elements, but not necessarily much else.
/// Examples are raw-pointer based arrays (with non-unit stride, otherwise one
/// can just use the pointers directly).
/// Usage: Add:
///
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,const_iterator,const)
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,iterator,)
///
/// to class declaration to make both const and non-const iterators which wrap
/// calls to operator [].
#define CX_IMPLEMENT_ITERATOR_IndexInBracket(HostType,IteratorType,ConstOrNothing) \
   struct IteratorType { \
      HostType ConstOrNothing \
         *m_pArrayObj; \
      size_t \
         m_Index; \
      IteratorType &operator ++() { ++m_Index; return *this; } \
      bool operator != (IteratorType const& other) const { \
         return this->m_pArrayObj != other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      bool operator == (IteratorType const& other) const { \
         return this->m_pArrayObj == other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      auto operator *() -> decltype((*m_pArrayObj)[m_Index]) { return (*m_pArrayObj)[m_Index]; } \
   }; \
   IteratorType begin() ConstOrNothing { return IteratorType{this, 0}; } \
   IteratorType end() ConstOrNothing { return IteratorType{this, size()}; }
   // ^--- note: == operator not strictly required.
//       typename std::decay<decltype(std::declval<HostType>()[0])>::type \
//       ConstOrNothing &operator *() { return (*m_pArrayObj)[m_Index]; } \


// ...and add some template-programming stuff, too. One can argue that this is
// indeed about iteration, although the nature is a bit different, of course.
//
// Note: This is intentionally put into a separate namespace (not ct). I have
// often used completely different things called "index lists"....
namespace cx {
   // Make an implementation of something like std::index_sequence, which, unfortunately,
   // is >= C++14 only. E.g., MakeIndexList<5>() yields a return type
   // ct::TIndexList<0, 1, 2, 3, 4>, which can be used in template<class
   // TIndexList<Indices...> > in variable template argument lists or similar.
   //
   // Note to self: A common usage pattern of this works like that:
   //
   //       template<unsigned... Is, class FTupleLike>
   //       void WriteTupleImpl(std::ostream &out, FTupleLike const &A, index_list_t<Is...>) const {
   //          using std::get;
   //          int dummy[] = {
   //             ((out << ((Is == 0)? "" : get<1>(m_SequenceSeptors)) << get<Is>(A)), 0)...
   //             // ^-- the (stuff,0)... which is repeated is the array element, each of
   //             //     which evaluates to 0 because this is a comma-operator. The array
   //             //     thing is apparently the canonical way in C++11 to execute a
   //             //     statment in input order for each entry in the parameter pack.
   //          };
   //          (void)dummy; // suppress unused warning.
   //       }
   //
   //       template<class... Args>
   //       std::ostream &WriteTuple(std::ostream &out, std::tuple<Args...> const &A) const {
   //          // ...
   //          WriteTupleImpl(out, A, index_list_for_t<Args...>());
   //          // ...
   //       }
   //
   // There are two functions: The originally called one (here: WriteTuple) just
   // creates the template index list, and then defers to a second function
   // which takes it as parameter pack template argument (here: WriteTupleImpl).
   // Inside this function, one can then pack-expand the index list template
   // (here: Is) to get a comma- separate sequence for all indices.
   //
   // E.g. for Is = <0, 1, 2> (index list of three compile-time arguments), we'd
   // see those expansions:
   //
   //    fragment `m[Is]...`  --> turns into `m[0], m[1], m[2]`
   //    fragment `func(m[Is]...)` -->  turns into  `func(m[0], m[1], m[2])`
   //    fragment `func((m[Is])...)` -->  turns into  `func((m[0]), (m[1]), (m[2]))`
   //    // ...i.e., (doesn't do anything here, but sometimes the extra parenthesis are needed to disambiguate)
   //    fragment `(out << get<Is>(a))...` --> turns into `(out << get<0>(a)), (out << get<1>(a)), (out << get<2>(a))`
   //    fragment `(func(x,Is,y), 0)...` --> turns into `(func(x,0,y), 0), (func(x,1,y), 0), (func(x,2,y), 0)`
   //
   // etc. So the pack expansions are similar to text replacements: take
   // everything in the expression right before the "..." (use parenthesis to
   // define enclosed scope in detail) and call it "Expr". Then for each entry
   // in the index list expansion, one gets a new entry in which the occurrence
   // of the index list (here `Is`) in the expression is replaced by its value,
   // and these new expressions are written next to each other with a comma
   // between any two---regardless of context.
   //

   typedef unsigned index_t;
   namespace iter_detail {
      template<index_t... Integers>
      struct TIndexList {};

      template<index_t iFirst, index_t N, index_t... OtherArgs>
      struct _TIndexListHelper {
         static_assert(N <= 200000, "probably an overflow...");
         static_assert(N != 0, "should not be here if N = 0");
         typedef typename _TIndexListHelper<iFirst, (N-1), (iFirst+N-1), OtherArgs...>::type type;
      };

      template<index_t iFirst, index_t... OtherArgs>
      struct _TIndexListHelper<iFirst, 0, OtherArgs...> {
         typedef TIndexList<OtherArgs...> type;
      };
   }

   template<index_t... Args>
   constexpr index_t _IndexCountQ(iter_detail::TIndexList<Args...>) { return sizeof...(Args); }

   template<index_t... Args>
   using index_list_t = iter_detail::TIndexList<Args...>;

   template<class... Args>
   using index_list_for_t = typename iter_detail::_TIndexListHelper<0, sizeof...(Args)>::type;
   // ^-- this typedef is sometimes easier to use than the MakeIndexList
   //     alternative. A std::index_sequence_for is also in std C++, but only
   //     starting at C++17.
   template<size_t N>
   using index_list_of_size_t = typename iter_detail::_TIndexListHelper<0, N>::type;

   // if iFirst < iLast, makes integer sequence {iFirst, iFirst+1, ..., iLast-1},
   // otherwise an empty integer sequence.
   template<size_t iFirst, size_t iLast>
   using index_list_for_range_t = typename iter_detail::_TIndexListHelper<iFirst, ((iLast > iFirst) ? (iLast - iFirst) : 0)>::type;

//    template<size_t iFirst, size_t iLast>
//    using index_list_for_range_t = typename std::conditional<
//       iLast <= iFirst,
//          iter_detail::TIndexList<>,
//          typename iter_detail::_TIndexListHelper<iFirst, iLast-iFirst>::type
//       >::type;
   // ^-- hm... seems like with std::conditional both type arguments are still fully evaluated.
   //     This thing explodes with integer overflow in the template arguments if index_list_for_range_t<3, 2>
   //     or something is evaluated. I guess the other `std::conditional`'s I used could then also
   //     create some exponential expansion of object complexity and compilation time... might need to
   //     have a look at that.

   template<index_t N>
   typename iter_detail::_TIndexListHelper<0,N>::type make_index_list() {
      return typename iter_detail::_TIndexListHelper<0,N>::type();
   }

} // namespace cx


// UPCOMING: more index list stuff. List element test, and functions for various
// decompositions / modifications of the lists. Not super elegant, but all this stuff also works in C++11.
namespace cx {

   typedef unsigned index_t;

   namespace iter_detail {
      /// Provides a constexpr function
      /// ``
      ///    template _IsEqualToAny<i0, i1, i2, ...>::eval(j)
      /// ``
      /// which returns whether the function argument `j` compares equal (e.g., j == i0)
      /// to any of the template arguments i0, i1, ...
      template<index_t...> struct _IsEqualToAny;
      template<> struct _IsEqualToAny<> {
         constexpr static bool eval(index_t) { return false; }
      };
      template<index_t ref, index_t... other_ref> struct _IsEqualToAny<ref, other_ref...> {
         constexpr static bool eval(index_t trial) { return (trial == ref) || _IsEqualToAny<other_ref...>::eval(trial); }
      };
   }

   template<index_t... Integers>
   constexpr bool ContainsQ(index_t trial) { return iter_detail::_IsEqualToAny<Integers...>::eval(trial); }

   template<index_t... Integers>
   constexpr bool ContainsQ(index_t trial, index_list_t<Integers...>) { return iter_detail::_IsEqualToAny<Integers...>::eval(trial); }

   template<class FIndexList>
   constexpr bool ContainsQ(index_t trial) { return ContainsQ(trial, FIndexList()); }

   namespace iter_detail {
      /// Provides a constexpr function
      /// ``
      ///    template _GetArgN<iArg>::eval(a0, a1, a2, ...)
      /// ``
      /// which returns the iArg'th function argument. Will yield compiler error
      /// if either there are no function arguments, or less function arguments than needed for iArg.
      ///
      /// Note: no r-value refs, forwarding etc. Arguments are returned by value
      /// (meant to be used with integer lists etc), to allow constexpr
      /// evaluation (even in C++11)
      template<index_t iEntry> struct _GetArgN {
         using _GetArgN_Next = _GetArgN<iEntry-1>;
         template<class Arg0, class... OtherArgs>
         constexpr static auto eval(Arg0, OtherArgs... other_args) -> decltype(_GetArgN_Next::eval(other_args...)) {
            return _GetArgN_Next::eval(other_args...);
         }
      };
      template<> struct _GetArgN<0> {
         template<class Arg0, class... OtherArgs>
         constexpr static Arg0 eval(Arg0 arg0, OtherArgs...) { return arg0; }
      };
   }

   /// compile-time index list with basic querying functionality (C++11-compatibly) to facilitate derivative constructions
   template<index_t... Integers>
   struct TIndexListQf {
      /// Number of controlled integer template arguments
      enum { nEntries = sizeof...(Integers) };
      /// Returns whether `trial` is an element of this list (i.e., coincides
      /// with any of the indices in the template argument list)
      constexpr static bool ContainsQ(index_t trial) { return iter_detail::_IsEqualToAny<Integers...>::eval(trial); }
      /// Returns `iArg`'th element of the list
      template<index_t iArg> constexpr static index_t GetArgN() { return iter_detail::_GetArgN<iArg>::eval(Integers...); }
   };

   /// convert a basic TIndexList to a TIndexListQf with query functionality (variable version)
   template<index_t... Integers>
   TIndexListQf<Integers...> _ConvertQfil(iter_detail::TIndexList<Integers...>) { return TIndexListQf<Integers...>{}; }
   /// convert a TIndexListQf with query functionality to a basic TIndexList (with no functionality) (variable version)
   template<index_t... Integers>
   iter_detail::TIndexList<Integers...> _ConvertNfil(TIndexListQf<Integers...>) { return iter_detail::TIndexList<Integers...>{}; }
   /// convert a basic TIndexList to a TIndexListQf with query functionality (type version)
   template<class FIndexList1> using _convert_qfil_t = decltype(_ConvertQfil(std::declval<FIndexList1>()));
   /// convert a TIndexListQf with query functionality to a basic TIndexList (type version)
   template<class FIndexList1> using _convert_nfil_t = decltype(_ConvertNfil(std::declval<FIndexList1>()));

   namespace iter_detail {
      using std::get;
      template<class BaseIndices, class RetainQ, index_t N = BaseIndices::nEntries, index_t... OtherArgs>
//       template<class BaseIndices, class RetainQ, index_t N = _IndexCountQ(BaseIndices()), index_t... OtherArgs>
      struct _FilterIndices {
//          static constexpr index_t iBaseIndex = BaseIndices::template GetArgN<N-1>();
         enum { iBaseIndex = BaseIndices::template GetArgN<N-1>() };
         using type = typename std::conditional<
            RetainQ::eval(N-1, iBaseIndex),
            typename _FilterIndices <BaseIndices, RetainQ, (N-1), index_t(iBaseIndex), OtherArgs...>::type,
            typename _FilterIndices <BaseIndices, RetainQ, (N-1), OtherArgs...>::type
         >::type;
      };

      template<class BaseIndices, class RetainQ, index_t... OtherArgs>
      struct _FilterIndices <BaseIndices, RetainQ, 0, OtherArgs...> {
//          typedef TIndexListQf<OtherArgs...> type;
         typedef iter_detail::TIndexList<OtherArgs...> type;
      };

      template<class Indices> struct _IsContainedIn { constexpr static bool eval(index_t, index_t value) { return Indices::ContainsQ(value); } };
      template<class Indices> struct _IsDisjunctFrom { constexpr static bool eval(index_t, index_t value) { return !Indices::ContainsQ(value); } };
      template<index_t SlotN> struct _TakeFirstN { constexpr static bool eval(index_t pos, index_t) { return pos < SlotN; } };
      template<index_t SlotN> struct _DropFirstN { constexpr static bool eval(index_t pos, index_t) { return pos >= SlotN; } };
      // ^-- could also do things like testing for index values instead of slot indices. not sure if needed atm.
   }

   /// `index_list_reject2nd_t<B = TIndexListQf<...>, E = TIndexListQf<...> >`
   /// evaluates to another TIndexListQf index list which contains all indices
   /// of `B`, in original order, provided they are *NOT* included in `E`
   template<class BaseIndices, class Indices2>
   using index_list_reject2nd_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_IsDisjunctFrom<_convert_qfil_t<Indices2> > >::type;

   /// `index_list_retain2nd_t<B = TIndexListQf<...>, I = TIndexListQf<...> >`
   /// evaluates to another TIndexListQf index list which contains all indices
   /// of `B`, in original order, provided they are *ALSO* included in `I`
   template<class BaseIndices, class Indices2>
   using index_list_retain2nd_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_IsContainedIn<_convert_qfil_t<Indices2> > >::type;

   /// Retain the first `N` indices in `BaseIndices`.
   /// @param BaseIndices A TIndexListQf<...> object with template arguments
   ///        providing the actual index list
   /// @param N maximum number of indices to retain
   /// @returns TIndexListQf of length <= N containing up to N of the first
   ///        indices in BaseIndices
   template<class BaseIndices, index_t N>
   using index_list_take_first_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_TakeFirstN<N> >::type;

   /// Evaluates the index list remaining after dropping the first `N` indices of `BaseIndices`.
   /// `BaseIndices` and return value are TIndexListQf templates.
   template<class BaseIndices, index_t N>
   using index_list_drop_first_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_DropFirstN<N> >::type;

} // namespace cx



namespace cx {


   template<class... T> class MagicTypeRevealer;
   #ifdef INCLUDE_ABANDONED
   // template<class T> class MagicTypeRevealer;
   //
   // int bla() {
   //    auto Is = MakeIndexList<5>();
   //    MagicTypeRevealer<decltype(Is)>();
   // }
   #endif // INCLUDE_ABANDONED



   // try to get the actual element type hidden behind something returns as reference,
   // by stripping off reference-ness and const-ness. Of course, if the actual type
   // indeed *is* const, this will remove it, too.
   // (note: there is a std::remove_cvref, but only in >= C++20; I call this a bit
   // differently to reduce the chances of ambiguity with the std version if someone
   // DOES use >= C++20. I did not think about this very deeply, but atm it seems
   // quite possible that multiple variants could get into competition due to ADL
   // or using namespace ...;
   template<class FDeclType>
   using without_cvref_t = typename std::remove_cv<typename std::remove_reference<FDeclType>::type>::type;

   // no std::remove_reference_t in C++11...
   template<class FDeclType>
   using without_ref_t = typename std::remove_reference<FDeclType>::type;

   // attempts to identify the type of the actual objects referenced by iterator
   // range [sequence.begin(), sequence.end()]
   template<class FSequence>
   using sequence_element_type_t = without_cvref_t<decltype(*std::declval<FSequence>().begin())>;


   // a <= C++ shortcut for the latter std library's std::enable_if_t
   template<bool B, class T = void>
   using enable_if_t = typename std::enable_if<B,T>::type;

   // unwrap_decay_t recycled from:
   //
   //     https://en.cppreference.com/w/cpp/utility/tuple/make_tuple
   //
   // It just runs through the normal std::decay, except if doing so results in
   // a reference wrapper. In that case, the reference is unwrapped and an
   // actual ref is created.
   //
   // This is process is what make_tuple() uses to decide on output types. And
   // that is probably rather important when it comes to matters like not
   // creating dangling r-value references...
   template <class T> struct unwrap_ref{ using type = T; };
   template <class T> struct unwrap_ref<std::reference_wrapper<T> > { using type = T&; };
   template <class T> using unwrap_ref_t = typename unwrap_ref<T>::type;

   template <class T> using decay_t = typename std::decay<T>::type;
   template <class T> using decay_and_unwrap_refs_t = unwrap_ref_t<decay_t<T> >;
   // ^-- >= C++20 has std::unwrap_ref_decay_t


//    template<class FTupleLike>
//    inline constexpr size_t tuple_size_v = std::tuple_size<FTupleLike>::value;

#if __cplusplus >= 201703L
   template<class FTupleLike> inline constexpr size_t _NumElemQ = std::tuple_size<FTupleLike>::value;
   template<class FTupleLike> inline constexpr size_t _ElemNumQ = std::tuple_size<FTupleLike>::value;
   // ^-- can't decide. Also count vs num...?
#endif // __cplusplus >= 201703L
   template<class FTupleLike> using _NumElemT = typename std::tuple_size<without_cvref_t<FTupleLike> >;
   template<class FTupleLike> using _ElemNumT = typename std::tuple_size<without_cvref_t<FTupleLike> >;

   template<size_t iElem, class FTupleLike> using _ElemTypeQ = typename std::tuple_element<iElem, FTupleLike>::type;

//    template<size_t iElem, class FTupleLike>
//    auto _ElemQ(FTupleLike &&obj) -> decltype(std::get<iElem>(std::forward(obj))) { return std::get<iElem>(std::forward(obj)); };
   // ^-- _ElemValueQ?
   using std::get;

   template<size_t iElem, class FTupleLike>
   auto _ElemFwdQ(FTupleLike &&obj) -> decltype(get<iElem>(std::forward<FTupleLike>(obj))) { return get<iElem>(std::forward<FTupleLike>(obj)); }
   // ^--- I am still not good at modern C++... and I wonder if there are ever
   // cases where I should *not* use get<...>(forward)? needs some more
   // thinking, I guess...
   // maybe cross check: https://en.cppreference.com/w/cpp/utility/forward
   template<size_t iElem, class FTupleLike>
   auto _ElemValQ(FTupleLike const &obj) -> decltype(get<iElem>(obj)) { return get<iElem>(obj); }

   template<class FTupleLike, class FIndexList = index_list_of_size_t<_NumElemT<FTupleLike>::value > >
   auto _MakeElemIndexList(FTupleLike const &) -> FIndexList { return FIndexList(); }



   // FIXME: I think the full implementations of all this stuff are in ~/dev/cx/tests/list_indices, or somewhere around.

   namespace iter_detail {
      using std::get;
      template<class FTupleLike, size_t... Is>
      auto _AsTupleOfRefsImpl(FTupleLike &&obj, cx::index_list_t<Is...>)
      -> decltype(std::forward_as_tuple(get<Is>(std::forward<FTupleLike>(obj))...))
      {
      //    return {FixedSizeArray[Is]...};
         return std::forward_as_tuple(get<Is>(std::forward<FTupleLike>(obj))...);
      }

      template<class FTupleLike, size_t... Is>
      auto _AsTupleOfValsImpl(FTupleLike &&obj, cx::index_list_t<Is...>)
      -> std::tuple<std::decay<decltype(get<Is>(obj))>...>
      {
         return {get<Is>(obj)...};
      }
}
   // make an std::tuple of forwarding references for the elements of obj (e.g., a fixed size array)
   template<class FTupleLike>
   auto _AsTupleOfRefs(FTupleLike &&obj)
   -> decltype(iter_detail::_AsTupleOfRefsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj)))
   {
      return iter_detail::_AsTupleOfRefsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj));
   }

   // make an std::tuple of actual values for the input elements (i.e., not reference or
   // forwarding references)
   template<class FTupleLike>
   auto _AsTupleOfVals(FTupleLike &&obj)
   -> decltype(iter_detail::_AsTupleOfValsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj)))
   {
      return iter_detail::_AsTupleOfValsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj));
   }


   namespace iter_detail {
      template<size_t I, class T, class... Ts>
      struct _TSelectNthTypeHelper {
         using type = typename _TSelectNthTypeHelper<I-1, Ts...>::type;
      };

      template<class T, class... Ts>
      struct _TSelectNthTypeHelper<0, T, Ts...> {
         using type = T;
      };

      template<size_t I>
      struct _TSelectNthArgHelper {
         template<class T, class... Args>
         static constexpr auto eval(T&&, Args... other)
         -> decltype(_TSelectNthArgHelper<I-1>::eval(std::forward(other)...))
         {
            return _TSelectNthArgHelper<I-1>::eval(std::forward(other)...);
         }
      };

      template<>
      struct _TSelectNthArgHelper<0> {
         template<class T, class... Args>
         static constexpr T&& eval(T &&first, Args... other) { return first; }
      };
   }

   /// extracts the nth type of a parameter pack (lamely)
   template <size_t I, class... Ts>
   using nth_type = typename iter_detail::_TSelectNthTypeHelper<I, Ts...>::type;

   /// nth_arg<I>(a0, a1,...) returns the I'th of its function arguments.
   template <size_t I, class... Args>
   auto nth_arg(Args... args)
   -> decltype(iter_detail::_TSelectNthArgHelper<I>::eval(std::forward(args)...)) {
      return iter_detail::_TSelectNthArgHelper<I>::eval(std::forward(args)...);
   }

   // see here: for more clever versions (for types):
   //
   //  https://ldionne.com/2015/11/29/efficient-parameter-pack-indexing/
   // And the comment that it apparently could be implemented
   // in a few hours into an existing compiler...  there apparently are even clang and g++ compiler extensions!

//    template<size_t N, Args..., unsigned... Is>
//    auto forward_first_n(std::tuple<Args...> &&objs, index_list_t<Is...> _ = index_list_of_size_t<(N<=sizeof...(Args))? N : sizeof...(Args)>) {
//       return std::forward_as_tuple(std::get<Is>(objs));
//    }
//    template<size_t N, class FTupleLike>
//    auto forward_first_n(FTupleLike &&objs), index_list_t<Is...> _ = index_list_of_size_t<(N<=sizeof...(Args))? N : sizeof...(Args)>) {
//       return std::forward_as_tuple(std::get<Is>(objs));
//    }
//    template<size_t N, class FTupleLike, unsigned... Is>
//    auto forward_first_n(
//          FTupleLike &&objs,
//          index_list_t<Is...> _ = index_list_of_size_t<(N <= std::tuple_size<FTupleLike>::value)? N : std::tuple_size<FTupleLike>::value>) {
//       return std::forward_as_tuple(std::get<Is>(std::forward(objs))...);
//    }
// N <= std::tuple_size<FTupleLike>::value)? N : std::tuple_size<FTupleLike>::value
//    template<size_t N, class FTupleLike>
//    constexpr size_t _MinSizeOrN() { return ((N <= std::tuple_size<FTupleLike>::value)? N : std::tuple_size<FTupleLike>::value); };
//
//    template<size_t N, class FTupleLike, unsigned... Is>
//    auto forward_first_n(
//          FTupleLike &&objs,
//          index_list_t<Is...> _ = index_list_of_size_t<_MinSizeOrN<N,FTupleLike>() >()) {
//       return std::forward_as_tuple(std::get<Is>(std::forward(objs))...);
//    }

//    template<size_t N, class FTupleLike, index_t... Is>
//    auto forward_first_n(
//          FTupleLike &&objs,
//          index_list_t<Is...> = index_list_of_size_t<N>()) {
// //       MagicTypeRevealer<FTupleLike, index_list_t<Is...> >();
//       return std::forward_as_tuple(std::get<Is>(std::forward(objs))...);
//    }
//    template<size_t N, class FTupleLike, index_t... Is>
//    auto forward_first_n_impl(FTupleLike &&objs, index_list_t<Is...>) {
// //    auto forward_first_n_impl(FTupleLike &&objs, index_list_t<Is...> = index_list_of_size_t<N>()) {
// //       MagicTypeRevealer<FTupleLike, index_list_t<Is...> >();
//       using std::get;
// //       return std::forward_as_tuple(get<Is>(std::forward(objs))...);
//       return std::forward_as_tuple(get<Is>(std::forward(objs))...);
//    }
//    template<size_t N, class FTupleLike>
//    auto forward_first_n(FTupleLike &&objs) {
// //       return forward_first_n_impl<N>(std::forward(objs), index_list_of_size_t<N>());
//       return forward_first_n_impl<N>(objs, index_list_of_size_t<N>());
// //       return forward_first_n_impl<N>(objs);
//    }
   template<size_t N, class FTupleLike, index_t... Is>
   auto forward_first_n_impl(FTupleLike &&objs, index_list_t<Is...>) {
//    auto forward_first_n_impl(FTupleLike &&objs, index_list_t<Is...> = index_list_of_size_t<N>()) {
//       MagicTypeRevealer<FTupleLike, index_list_t<Is...> >();
      using std::get;
//       return std::forward_as_tuple(get<Is>(std::forward(objs))...);
      return std::forward_as_tuple(get<Is>(objs)...);
   }
   template<size_t N, class FTupleLike>
   auto forward_first_n(FTupleLike &&objs) {
//       return forward_first_n_impl<N>(std::forward(objs), index_list_of_size_t<N>());
      return forward_first_n_impl<N>(objs, index_list_of_size_t<N>());
//       return forward_first_n_impl<N>(objs);
   }

   // as in C++14...
   template<bool B, class T, class F>
   using conditional_t = typename std::conditional<B,T,F>::type;


   namespace iter_detail {
      template<class FOut, class FBinaryOp>
      FOut _FoldLeftHelper(FBinaryOp op, FOut const &neutral) {
         return neutral;
      }

      template<class FOut, class FBinaryOp, class Arg0, class Arg1, class... Args>
      FOut _FoldLeftHelper(FBinaryOp op, Arg0 a0, Arg1 a1, Args... args) {
         return _FoldLeftHelper<FOut,FBinaryOp>(op, op(std::move(a0), a1), args...);
      }
   }

   // Notes: _pp -> parameter pack. Output type T is derived from first argument.
   // (this first argument is also what you will get if there are no `real` arguments)
   template<class T, class FBinaryOp, class... Args>
   T fold_left_pp(FBinaryOp op, T neutral, Args... args) {
      return iter_detail::_FoldLeftHelper<T,FBinaryOp>(op, neutral, args...);
   }

   template<class T, class... Args>
   T sum_pp(T neutral, Args... args) {
      return fold_left_pp<T>(std::plus(), neutral, args...);
   }


}





#endif // CX_ITER_TOOLS_H
