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

// This is a reference implementation of the "CxAssertFail" function used in
// CxDefs.h's assertion functions. It is only used to make the distribution
// self-contained: a function with signature
//
//    CxAssertFail(char const *pExpr, char const *pFile, int iLine)
//
// needs to be *somewhere* in the program, but it does not have to be this one.
// If you have a different way of handling assertion failures, you need
// not use/compile in this file, but can instead just substitute your own
// implementation.
//
// UPDATE:
// - Some other error raising function implementations also made their way here...

#include <iostream>
#include <sstream>
#include <stdexcept>
#include "CxDefs.h"
#include "format.h"
// ^-- hm... format.h was not required without the extra _RaiseMemoryAllocError
// stuff. On the other hand, I probably really don't want to make a program
// without this.

using std::size_t;


void CxAssertFail(char const *pExpr, char const *pFile, int iLine)
{
   std::stringstream
      str;
   if (pFile) {
      str << pFile << ":" << iLine << ": ";
   }
   str << "Assertion Failed: '" << pExpr << "'";
   std::cerr << str.str() << std::endl;
#ifdef _DEBUG
   DEBUG_BREAK
#endif // _DEBUG
   throw std::runtime_error(str.str());
}



namespace ct {
   double _fSizeMb(size_t nBytes) {
      return double(nBytes)/double(size_t(1) << size_t(20));
   }

   void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount) {
      throw std::runtime_error(fmt::format(
         "{}: memory allocation failed ({:.2f} MB requested; {} objects of size {}).",
         pFnName, _fSizeMb(ValueTypeSize * NewReserveCount), NewReserveCount, ValueTypeSize));
   }

   void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount, size_t PrevReservedCount) {
      throw std::runtime_error(fmt::format(
         "{}: memory allocation failed ({:.2f} MB requested; {} objects of size {}; reallocation: {}; prev size: {:.2f} MB).",
         pFnName, _fSizeMb(ValueTypeSize * NewReserveCount), NewReserveCount, ValueTypeSize,
         (PrevReservedCount != 0)? "yes" : "no", _fSizeMb(ValueTypeSize * PrevReservedCount))
      );
   }
}



#include <algorithm> // for std::stable_sort

// well... this one *most definitely* should not be here. It is declared in
// CxTypes.h, but that one does not have an implementation file, and it doesn't
// really fit anywhere else, either.
namespace ct {

struct FPredArgSort {
   FPredArgSort(double const *pVals_, size_t nValSt_, bool Reverse_) : pVals(pVals_), nValSt(nValSt_), Reverse(Reverse_) {};
   bool operator () (size_t i, size_t j) const { return Reverse? (iVal(i) > iVal(j)) : (iVal(i) < iVal(j));  }
private:
   double iVal(size_t i) const { return pVals[nValSt * i]; }
   double const *pVals;
   size_t nValSt;
   bool Reverse;
};

// find permutation which sorts values in pVals[nValSt*i].
void ArgSort1(size_t *pOrd, double const *pVals, size_t nValSt, size_t nVals, bool Reverse) {
   for (size_t i = 0; i < nVals; ++ i)
      pOrd[i] = i;
   std::stable_sort(pOrd, pOrd+nVals, FPredArgSort(pVals, nValSt, Reverse));
}


} // namespace ct
