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

/// @file CxArgSort.h
///
/// Provides simple support routines for ArgSort, ArgMin and ArgMax. The latter come
/// from IvAnalysis.cpp originally, from code which has since been moved to
/// CtAnalysisEos.cpp

#include <stdexcept>
#include <algorithm>

#ifndef CX_ARGSORT_H
#define CX_ARGSORT_H

namespace ct {

   template<class FSortKey, class FCompareLess>
   struct TPredArgSort {
      TPredArgSort(FSortKey const *pVals_, size_t nValSt_, bool Reverse_, FCompareLess const &CompareLess_)
         : m_pVals(pVals_), m_nValSt(nValSt_), m_Reverse(Reverse_), m_CompareLess(CompareLess_)
      {};

      bool operator () (size_t i, size_t j) const {
         if (m_Reverse)
            return m_CompareLess(Val(j), Val(i));
         else
            return m_CompareLess(Val(i), Val(j));
      }
   private:
      FSortKey const
         *m_pVals;
      size_t
         m_nValSt;
      bool
         m_Reverse;
      FCompareLess const
         &m_CompareLess;
      FSortKey const &Val(size_t i) const {
         return m_pVals[m_nValSt * i];
      }
   };

   enum FArgSortFlags {
      /// if set in Flags, ArgSortT will compute the permutation for sorting the
      /// input sequence into order of decreasing values, rather than increasing
      /// values.
      SORT_Reversed = 0x0001,
      /// if set, ArgSortT will not guarantee that input values which compare
      /// equivalent according to the sorting predicate CompareLess (i.e.,
      /// neither CompareLess(a,b) nor CompareLess(b,a) is true) remain in the
      /// relative order they had in the input sequence.
      SORT_NonStableIsOk = 0x0002,
      /// retain relative order of equivalent input sequence elements; this is
      /// the default behavior/ The flag is provided only to facilitate code
      /// verbosity in situations when either behavior matters.
      SORT_Stable = 0x0000,
   };

   /// Find permutation which sorts values in pVals[nValSt*i].
   ///
   /// Arguments:
   /// - pOrd: stores the output permutation; on exit, the sequences
   ///   pVals[pOrd[i]] for i = 0, 1, ..., (nVals - 1) will be ordered.
   /// - pVals: pointer to input sequence; input element `i` is given
   ///   by pVals[nValSt * i], and i = 0, 1, ..., (nVals - 1)
   /// - nValSt: array stride between sequence elements in pVals to
   ///   consider for sorting. In the case of dense linear packing,
   ///   use nValSt = 1.
   /// - Flags: bit field of SORT_Reversed, SORT_NonStableIsOk, SORT_Stable
   template<class FSortKey, class FCompareLess>
   void ArgSortT(size_t *pOrd, FSortKey const *pVals, size_t nValSt,
                 size_t nVals, unsigned Flags,
                 FCompareLess CompareLess = std::less<FSortKey>())
   {
      for (size_t i = 0; i < nVals; ++ i)
         pOrd[i] = i;
      TPredArgSort<FSortKey, FCompareLess>
         Pred(pVals, nValSt, bool(Flags & SORT_Reversed), CompareLess);
      if (bool(Flags & SORT_NonStableIsOk))
         std::sort(pOrd, pOrd+nVals, Pred);
      else
         std::stable_sort(pOrd, pOrd+nVals, Pred);
   }

   /// return index of first element in sequence [first,last) for which
   /// no other element is smaller (a < b). (if multiple elements compare
   /// equal, the left-most (lowest index one) is returned)
   template<class FRandomAccessIt>
   size_t ArgMin(FRandomAccessIt first, FRandomAccessIt last)
   {
      if (first == last)
         throw std::runtime_error("ArgMin failed: input range is empty.");
      FRandomAccessIt
         it = first,
         itMin = it;
      for ( ; it != last; ++ it) {
         if (*it < *itMin)
            itMin = it;
      }
      return itMin - first;
   }


   /// return index of first element in sequence [first,last) which does
   /// not compare as "less" to any other element in the sequence.
   template<class FRandomAccessIt>
   size_t ArgMax(FRandomAccessIt first, FRandomAccessIt last)
   {
      if (first == last)
         throw std::runtime_error("ArgMax failed: input range is empty.");
      FRandomAccessIt
         it = first,
         itMax = it;
      for ( ; it != last; ++ it) {
         if (*itMax < *it)
            itMax = it;
      }
      return itMax - first;
   }


   template<class FArray>
   size_t ArgMin(FArray const &a) {
      if (a.empty())
         throw std::runtime_error("ArgMin failed: input array is empty.");
      size_t
         iMin = 0;
      typename FArray::value_type
         vMin = a[0];
      for (size_t i = 1; i != a.size(); ++ i) {
         if (a[i] < vMin) {
            iMin = i;
            vMin = a[i];
         }
      }
      return iMin;
   }

}

#endif // CX_ARGSORT_H
