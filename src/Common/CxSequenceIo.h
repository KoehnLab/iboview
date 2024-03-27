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

// Contains implementations of ostream-based (and maybe later fmt::BasicWriter
// based) generic implementations of formatted output operators (<<) for
// std::array / std::tuple and a simple interface to make additional ones for
// other containers (std::vector, std::set, std::map, etc.)
//
#ifndef CX_SEQUENCE_IO_H
#define CX_SEQUENCE_IO_H

#include <type_traits> // for remove_const / remove_reference
#include <array>
#include <utility>
#include <tuple>
#include <ostream>
#include <string>

#ifndef SEQUENCE_IO_PREDEFINE_LEVEL
   #define SEQUENCE_IO_PREDEFINE_LEVEL 0
#endif
// Note: SEQUENCE_IO_PREDEFINE_LEVEL may be pre-defined to control the level of
// include pollution.
//
// Level 0 (default):
//    - At level 0, this file will #include what is necessary to define ostream-based
//      formatting routines for various kinds of sequences.
//    - ...but it not actually declare & define actual generic realizations of
//      those for the standard containers, except the std::array (and
//      std::tuple, which is not really a container).
//    - Note: each container is three lines (see below), so there may not
//      actually be a terribly good reason to define them like this (apart
//      from debugging or development work)
// Level 1:
//    - At level 1, it will bring in #includes for the most widely used containers,
//      namely vector, map, set, and list, and define generic << operators for them.
//    - These still sit in namespace ct::sequence_io, so they should not by
//      themselves interfere with other code. To actually make them available,
//      you may want to add a
//
//          using ct::sequence_io::operator <<;
//
//      to a suitable local scope.
// Level 2:
//    - The file fill attempt to define working << operators (still in namespace
//      ct::sequence_io) for all standard containers, including the rather
//      rarely used ones (deque(*), unordered_multiset, forward_list, etc...)
//
// [*] (cgk's mystery: ...deque is actually a cool data structure, from a
//     algorithms & data structures point of view. But somehow, at least in my
//     past, it just happened soooo rarely that I would actualy want to use one
//     in a practical problem. And, just for the record: I *have* used plenty of
//     multisets and multimaps...)

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 1)
   // if we're trashing the include space so throughly anyway, I guess it makes
   // little difference whether one brings in all standard containers, too.
   #include <vector>
   #include <map>
   #include <set>
   #include <list>
#endif

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 2)
   #include <deque>
   #include <forward_list>
   #include <unordered_set>
   #include <unordered_map>
#endif


#include "CxIterTools.h" // for template index lists (index_list_t, index_list_for_t)


namespace ct {
namespace sequence_io {

// template<class FArray>
// void PrintArray(std::ostream &out, FArray const &A, std::array<char const *, 3> LeftMidRight = {"[",",","]"})
// {
//    using std::get;
//    out << get<0>(LeftMidRight);
//    for (size_t i = 0; i != A.size(); ++i) {
//       if (i != 0) out << get<1>(LeftMidRight);
//       out << A[i];
//    }
//    out << get<2>(LeftMidRight);
// }


// An object which makes a copy of an `std::ostream`'s flags & co (e.g., number
// formats) on construction, and restores them when it goes out of scope.
struct ostream_flag_guard_t {
   std::ostream *pout;
   int prec, width;
   std::ios::fmtflags flags;
   explicit ostream_flag_guard_t(std::ostream *pout_) : pout(pout_), prec(pout->precision()), width(pout->width()), flags(pout->flags()) {}
   ~ostream_flag_guard_t() { pout->precision(prec); pout->width(width); pout->flags(flags); }
};


// // Defines a minimal iterator interface to allow use in range-based for loops
// // and generic container printing mechanism. This one is for objects which have
// // a working [] operator to access elements, but not necessarily much else.
// // Examples are raw-pointer based arrays (with non-unit stride, otherwise one
// // can just use the pointers directly).
// //
// // Example to put in class:
// //    typedef TMinimalIterator_Bracket<FThisType const> const_iterator;
// //    const_iterator begin() const { return const_iterator{this, 0}; }
// //    const_iterator end() const { return const_iterator{this, size()}; }
// template<class FArrayLike, class FElement = typename FArrayLike::value_type>
// struct TMinimalIterator_Bracket {
//    TMinimalIterator_Bracket &operator ++() { ++m_Index; return *this; }
//    bool operator != (TMinimalIterator_Bracket const& other) const {  // <-- only != is strictly required (== optional).
//       return this->m_pArrayObj != other.m_pArrayObj || this->m_Index != other.m_Index;
//    }
//    FElement operator *() { return (*m_pArrayObj)[m_Index]; }
// protected:
//    FArrayLike
//       *m_pArrayObj;
//    size_t
//       m_Index;
// };

// ^-- UPDATE: that doesn't work because the class is not complete at that point.
//     Doing this requires a macro. I made one in CxIterTools.h. See:
//     CX_IMPLEMENT_ITERATOR_IndexInBracket

// A wrapper for giving raw-pointer based arrays an interface which is sufficient
// for WriteList (implements what is needed for the range-based for).
template<class FElement>
struct TArrayViewStrided {
   FElement const
      *m_pData;
   size_t
      // number of array elements
      m_Size;
   ptrdiff_t
      // in-memory distance between two adjacent elements we should consider
      // as adjacent. That is, actual element #i is at m_pData[i * m_Stride]
      m_Stride = 1;
   size_t size() const { return m_Size; }
   FElement const &operator[] (size_t i) { return m_pData[m_Stride * i]; }
//    typedef TMinimalIterator_Bracket<TArrayViewStrided const> const_iterator;
//    const_iterator begin() const { return const_iterator{this, 0}; }
//    const_iterator end() const { return const_iterator{this, size()}; }
#ifdef CX_IMPLEMENT_ITERATOR_IndexInBracket
   CX_IMPLEMENT_ITERATOR_IndexInBracket(TArrayViewStrided,const_iterator,const)
#endif
};


// emit a string (given in terms of iterator sequence), quoting it and (by default) performing some more
// common escape character replacements (e.g., actual newlines by the character sequence '\n')
template <class FStringIt>
void WriteQuoted(std::ostream &out, FStringIt const &first, FStringIt const &last, char QuoteLeft, char QuoteRight, bool EscapeQuotesOnly=false) {
   out << QuoteLeft;
   auto const &try_emit_escape = [&](char c, char from, char const *to, char to2) -> bool {
      if (c == from) {
         out << to;
         if (to2 != 0) out << to2;
         return true;
      }
      return false;
   };
   for (FStringIt it = first; it != last; ++ it) {
      auto c = *it;
      if (!EscapeQuotesOnly) {
         if (try_emit_escape(c, '\n', "\\n", 0)) continue;
         if (try_emit_escape(c, '\t', "\\t", 0)) continue;
         if (try_emit_escape(c, '\0', "\\0", 0)) continue;
         if (try_emit_escape(c, '\\', "\\\\", 0)) continue;
      }
      if (try_emit_escape(c, QuoteLeft, "\\", QuoteLeft)) continue;
      if (try_emit_escape(c, QuoteRight, "\\", QuoteRight)) continue;
      // still here? That means it's not one of the to-escape ones above.
      out << c;
      // btw: this is not supposed to be safe at all. Just to deal with some
      // more common escapes to help with debugging. DIY escapes are highly NOT
      // recommended in situations where you might encounter malicious data...
   }
   out << QuoteRight;
}

template<class FElement>
//static
void DefaultWriteElement_Unquoted(std::ostream &out, size_t iArg, FElement const &Arg) {
   out << Arg; (void)iArg;
}


template<class FElement>
//static
void DefaultWriteElement(std::ostream &out, size_t iArg, FElement const &Arg) {
   out << Arg; (void)iArg;
}
template<>
void DefaultWriteElement<std::string>(std::ostream &out, size_t iArg, std::string const &Arg) {
   WriteQuoted(out, Arg.begin(), Arg.end(), '\"', '\"'); (void)iArg;
}
template<>
void DefaultWriteElement<char const *>(std::ostream &out, size_t iArg, char const * const &pArg) {
   // assume 0-terminated string. count how long it is.
   size_t n = 0; for (char const *p = pArg; *p;) { n += 1; p += 1; }
   WriteQuoted(out, pArg, pArg+n, '\"', '\"'); (void)iArg;
}

// // suppress some unused-warnings. I think it is not particularly concerning if
// // certain of these functions are not, in fact, used in every single compilation unit.
// // update: took this out again and made DefaultWriteElement non-static instead. For i/o
// // stuff it doesn't much matter, I guess, although it does mean that we bring in those
// // fnuctions even if never used anywhere.
// inline void please_disregard_this() {
//    (void)DefaultWriteElement<char const *>;
//    (void)DefaultWriteElement<std::string>;
// }

template<class T>
struct remove_const_and_reference {
   // DefaultWriteElement is specilized for the base types (but takes
   // const refs as arguments to the elements). To get the right
   // function we therefore need to stripff off the reference and const-ness.
   typedef typename std::remove_reference<T>::type T1;
   typedef typename std::remove_const<T1>::type type;
};
template<class T>
using remove_const_and_reference_t = typename remove_const_and_reference<T>::type;



struct TSequenceWriter {
   typedef std::array<char const *, 3> FSeparators;
//    struct FSeparators { char const *Left, *Mid, *Right, *Key; }
   // separators between the actual direct elements of the container
   FSeparators m_SequenceSeptors = {"[",",","]"}; // LeftMidRight
   // separators for parts of an individual element; not normally used,
   // except for cases like associative containers (std::map and co),
   // in which we'd separate (key : value) pairs instead of just printing
   // them as 2-tuples.
   FSeparators m_ElementSeptors = {"",": ",""};
   typedef std::tuple<int,int,char> FFieldFormat; // width, prec, one of: f (fixed), e (scientific), g (general), x (hex)
   FFieldFormat m_FieldFormat = {-1, -1, 0};

   TSequenceWriter() {}
   explicit TSequenceWriter(FSeparators LeftMidRight_) : m_SequenceSeptors(LeftMidRight_) {}
   explicit TSequenceWriter(FSeparators LeftMidRight_, FFieldFormat FieldFormat_) : m_SequenceSeptors(LeftMidRight_), m_FieldFormat(FieldFormat_) {}
   explicit TSequenceWriter(FSeparators Sequence_LeftMidRight_, FSeparators Element_LeftMidRight_) : m_SequenceSeptors(Sequence_LeftMidRight_), m_ElementSeptors(Element_LeftMidRight_) {}

//    // arguments: target stream, index of element in sequence, element
//    template<class FElement>
//    using TWriteElementFn = std::function<void (std::ostream &, size_t, FElement const &)>;


   // change to iterator-based. Would allow sets and lists, but no longer the packed
   // index.... unless I add some sort of fake iterator class...
   template<class FArrayLike, class FWriteElementFn>
   std::ostream &WriteList(std::ostream &out, FArrayLike const &A, FWriteElementFn &&pEmitElement) const {
      ostream_flag_guard_t _flag_restore(&out);
      using std::get;
      out << get<0>(m_SequenceSeptors);
      size_t i = 0;
      for (auto const &Ai : A) {
         if (i != 0) out << get<1>(m_SequenceSeptors);
         if (get<0>(m_FieldFormat) >= 0)
            out.width(get<0>(m_FieldFormat));
         if (get<1>(m_FieldFormat) >= 0)
            out.precision(get<1>(m_FieldFormat));
         char c = get<2>(m_FieldFormat);
         if (c == 'g')
            out.unsetf(std::ios::floatfield); // switch mode to "general" (typically omits trailing zeros, unlike 'f')
         else if (c == 'f')
            out.setf(std::ios::fixed, std::ios::floatfield);
         else if (c == 'e')
            out.setf(std::ios::scientific, std::ios::floatfield);
//          out << Ai;
         pEmitElement(out, i, Ai);
         out.width(0);
         i += 1;
      }
      out << get<2>(m_SequenceSeptors);
      return out;
   }


   template<class FArrayLike>
   std::ostream &WriteList(std::ostream &out, FArrayLike const &A) const {
      return WriteList(out, A, &DefaultWriteElement<remove_const_and_reference_t<decltype(*A.begin())> >);
//       typedef typename std::remove_reference<decltype(*A.begin())>::type T1;
//       typedef typename std::remove_const<T1>::type T2;
//       // ^-- DefaultWriteElement is specilized for the base types (but takes
//       //     const refs as arguments to the elements). To get the right
//       //     function we therefore need to stripff off the reference and const-ness.
//       return WriteList(out, A, &DefaultWriteElement<T2>);
//       return WriteList(out, A, &DefaultWriteElement<typename FArrayLike::value_type>);
   }

   template<cx::index_t... Is, class FTupleLike>
   void WriteTupleImpl(std::ostream &out, FTupleLike const &A, cx::index_list_t<Is...>) const {
      using std::get;
      int dummy[] = {
         // that's apparently the recommended way before C++17 fold expressions
         // to perform statements on parameter packs...  The C-array
         // is there to ask the compiler to execute the statements in original sequence.
         // To make it an int array, we use the comma operator (..., 0).
         (
            (out << ((Is == 0)? "" : get<1>(m_SequenceSeptors)) << get<Is>(A)),
            0
         )...
      };
      (void)dummy; // suppress unused warning
   }

   template<class... Args>
   std::ostream &WriteTuple(std::ostream &out, std::tuple<Args...> const &A) const {
      ostream_flag_guard_t _flag_restore(&out);
      using std::get;
      out << get<0>(m_SequenceSeptors);
      WriteTupleImpl(out, A, cx::index_list_for_t<Args...>());
      out << get<2>(m_SequenceSeptors);
      return out;
   }

};


template<class T, size_t N> // <-- WARNING: this will not match any std::array<...> if you have another integral type than size_t there! (e.g., unsigned).
std::ostream &operator << (std::ostream &out, std::array<T,N> const &A) {
   return TSequenceWriter({"[",", ","]"}).WriteList(out, A);
}


template<class... Args>
std::ostream &operator << (std::ostream &out, std::tuple<Args...> const &A) {
   return TSequenceWriter({"(",", ",")"}).WriteTuple(out, A);
}


template<class K, class V>
void DefaultWriteMapElement(std::ostream &out, size_t i, std::pair<K,V> const &ei) {
   // map's value_type is std::pair<Key, Val>. Print them separated by `:`, similar to a python dict literal.
//    out << ei.first << ": " << ei.second;
   DefaultWriteElement<remove_const_and_reference_t<K> >(out, i, ei.first);
   out << ": ";
   DefaultWriteElement<remove_const_and_reference_t<V> >(out, i, ei.second);
}


typedef TSequenceWriter::FSeparators open_sep_close_t;


// // try some options.. wonky open/close sets and a user-defined ``I'd like to write
// // the individual elements like this: ...'' function. Here it also prints the
// // sequence index in [], in addition to the actual vector element.
// template<class T, class A> std::ostream &operator << (std::ostream &out, std::vector<T,A> const &seq) {
//    auto write1 = [&](std::ostream &out, size_t i, T const &ei) {
//       out << "[" << i << "]: " << ei;
//    };
//    return ct::sequence_io::TSequenceWriter({"|[","; ","]|"},{-1,6,'f'}).WriteList(out, seq, write1);
// }

// template<class T, class V, class A> std::ostream &operator << (std::ostream &out, std::map<T,V,A> const &seq) {
//    auto write1 = [&](std::ostream &out, size_t i, std::pair<T,V> const &ei) {
//       // map's value_type is std::pair<Key, Val>. Print them separated by `:`, similar to a python dict literal.
//       out << ei.first << ": " << ei.second;
//       (void)i;
//    };
//    return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq, write1);
// }


#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 1)
// define fully generic << operators for some more common standard containers: vector, list, set, map.

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::vector<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::list<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::set<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::multiset<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::map<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::multimap<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

#endif // if SEQUENCE_IO_PREDEFINE_LEVEL >= 1


#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 2)
   // define fully generic << operators for rarely encountered containers

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::deque<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::forward_list<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::unordered_set<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::unordered_multiset<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::unordered_map<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::unordered_multimap<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

#endif // if SEQUENCE_IO_PREDEFINE_LEVEL >= 1



} // namespace sequence_io
} // namespace ct


// template<class T, size_t N>
// std::ostream &operator << (std::ostream &out, std::array<T,N> const &A) {
//    ct::TSequenceWriter({"[",", ","]"}).WriteList(out, A);
//    return out;
// }

#endif // CX_SEQUENCE_IO_H
