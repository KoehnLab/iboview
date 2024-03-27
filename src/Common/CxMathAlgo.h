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

/// @file CxMathAlgo.h Assortment of random support algorithms for basic math.
#ifndef CX_MATH_ALGO_H
#define CX_MATH_ALGO_H

#include <stdexcept>
#include "CxVec3.h" // for _ThreshAlmostZero only
#include "format.h"

namespace ct {

/// Allows computing values of the inverse f^{-1}(y) of a scalar function f(x) =
/// y, given the forward function f(x) and a bracketing interval
/// x ∈ [xMin, xMax].
///
/// Example:
/// ```
///    #include <cmath>
///    #include <cassert>
///    #include <CxMathAlgo.h>
///
///    // compute inverse error function erfinv(p) near p=1 by bracketing
///    // std::erfc (in cmath since C++11).
///    double ErfInv_Near1(double p) {
///       ct::TImplicitInverseFn<double>
///          InvErfc(std::erfc, -5, 10);
///       return InvErfc(1 - p);
///    }
/// ```
///
///
template<class scalar_t>
struct TImplicitInverseFn
{
   using FScalarFn = scalar_t (scalar_t);

   TImplicitInverseFn(FScalarFn *pFn, scalar_t xMin, scalar_t xMax,
      scalar_t xThr = _ThreshAlmostZero<scalar_t>(), scalar_t yThr = _ThreshAlmostZero<scalar_t>()
   )
   : m_pFn(pFn), m_BracketMin(xMin), m_BracketMax(xMax), m_xThr(xThr), m_yThr(yThr)
   {}

   scalar_t operator () (scalar_t y) const {
      using std::abs;
      scalar_t
         xMin = m_BracketMin,
         xMax = m_BracketMax,
         eyMin = m_pFn(xMin) - y,
         eyMax = m_pFn(xMax) - y;
      scalar_t
         // that is for the cases in which we look for very small ys.
         // (e.g., as in the erfc case this was meant for)
         yThrEff = m_yThr * (abs(y) + m_yThr);
      for (size_t iIter = 0; iIter < 200; iIter += 1) {
         scalar_t
            xMid = (xMin + xMax)/2,
            eyMid = m_pFn(xMid) - y;
         // remember number of iterations needed. This is not thread safe; but
         // since it's only for info purposes, I guess that does not really
         // matter.
         m_nIterLast = iIter;
         if ((abs(eyMid) <= yThrEff) || (xMax - xMin) < m_xThr * (xMax + xMin))
            return xMid;
         if ((eyMin * eyMax) > 0)
            throw std::runtime_error(fmt::format("FInverseFn::eval(): at y = {}, bracket intervals do not have opposite signs (it = {}, x0 = {}, f(x0)-y = {}, x1 = {}, f(x1)-y = {})", y, iIter, xMin, eyMin, xMax, eyMax));
         if (eyMid * eyMin < 0) {
            // bracket on left side
            eyMax = eyMid;
            xMax = xMid;
         } else {
            // bracket on right side
            eyMin = eyMid;
            xMin = xMid;
         }
      }
      assert(!"FInverseFn: bracketing failed.");
      return (xMax + xMin)/2;
   }

   size_t nIterLast() const { return m_nIterLast; }
protected:
   FScalarFn
      *m_pFn;
   scalar_t
      m_BracketMin, m_BracketMax;
   scalar_t
      m_xThr, m_yThr;
   size_t mutable
      // for info purposes only.
      m_nIterLast;
};

using FImplicitInverseFn = TImplicitInverseFn<double>;


} // namespace ct

#endif // CX_MATH_ALGO_H
