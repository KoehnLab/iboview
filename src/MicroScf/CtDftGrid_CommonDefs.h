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

#ifndef CT_DFTGRID_COMMON_DEFS
#define CT_DFTGRID_COMMON_DEFS

#include "CxTypes.h"
#include "CxVec3.h"

namespace mig { // molecular integration grids

   typedef double
      FScalar;
   typedef ct::TVector3<FScalar>
      FVector3;

   enum {
      DFTGRID_nPeriodParams = 4 // number of per-periodic table-period parameters for DFT grids.
   };

   // returns std::isfinite if compiled with C++ version >= C++11, and indiscriminately 'true' if not.
   bool check_if_finite_if_supported(FScalar d);

   // returns the square of a number.
   template<class FNumber>
   inline FNumber sqr(FNumber const &x) {
      return x*x;
   }
} // namespace mig



#endif // CT_DFTGRID_COMMON_DEFS
