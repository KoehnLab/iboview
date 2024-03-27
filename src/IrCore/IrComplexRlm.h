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

/* IrComplexRlm.h v20200511 EST [mirage, Gerald Knizia] */
#ifndef IR_RLM_H
#define IR_RLM_H
// Utilities for dealing with semi-normalized complex regular solid harmonics R^l_m(\vec r).
// Used in conjunction with ECP integration.
#include <stdexcept>
#include <complex>

namespace ir_rlm {
   typedef std::complex<double> complex_double;
   inline int nCartY(int l) { return (l+1)*(l+2)/2; }
   inline int nCartX(int l) { return (l+1)*(l+2)*(l+3)/6; }
   // inline int iSlmX(int l, int m) { return l*l + l + m; }
   inline int iSlmX(int l, int m) { assert(l >= 0 && -l <= m && m <= l); return l*l + l + m; }
   inline int nSlmX(int l) { return (l+1)*(l+1); }
   inline int nSlmY(int l) { return 2*l+1; }

   void EvalRlmX(complex_double *RlmX, double x, double y, double z, int MaxL);
   void EvalRlmX_Deriv1(complex_double *RlmX, double x, double y, double z, int MaxL);
   double GetRlmRenormFactor(int la, int ma);
   size_t _m2c(int la, int ma);
   void ConvertRlmToSlm1(complex_double *pOut, size_t so, complex_double const *pIn, size_t si, int l);
   double GetSphereIntegral3Rlm(int La, int Ma, int EcpL, int EcpM, int Lambda, int Mu);


} // namespace ir_rlm
#endif // IR_RLM_H
