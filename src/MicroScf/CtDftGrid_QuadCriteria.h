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


#ifndef CT_DFTGRID_QUAD_CRITERIA_H
#define CT_DFTGRID_QUAD_CRITERIA_H

#include "CxMemoryStack.h"

namespace mig {

// one set of numerical quadrature criteria
struct FQuadCriteria
{
   static size_t const
      N = 8; // 6 and 7 seem to work better with multigrid than 5.
//       N = 7;
//       N = 6;
//       N = 5;
   double
      fValues[N];

   double
      m_rexps[N],
      m_wexps[N];

   FQuadCriteria();

   // evaluate with given densities
   void Eval(double *di, double *wi, double *ri, size_t n, int iElement, ct::FMemoryStack &Mem);
   // compute densities and evaluate.
   void Eval(double *ri, double *wi, size_t n, double fDensityShift, int iElement, ct::FMemoryStack &Mem);

   // combined deviation from reference values.
   double fResidual(FQuadCriteria const &Ref) const;

   void operator += (FQuadCriteria const &other);
   void operator *= (double f);
};

} // namespace mig

#endif // CT_DFTGRID_QUAD_CRITERIA_H
