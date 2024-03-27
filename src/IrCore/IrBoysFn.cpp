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

// IrBoysFn is free software. It comes without any warranty, to the extent
// permitted by applicable law. You may use it, redistribute it and/or modify
// it, in whole or in part, provided that you do so at your own risk and do not
// hold the developers or copyright holders liable for any claim, damages, or
// other liabilities arising in connection with the software.
//
// This file is part of the IR integral core, and supplied as supporting
// information of the article:
//
//  [1] Mieke Peels and Gerald Knizia - "Fast evaluation of two-center integrals
//  over Gaussian charge distributions and Gaussian orbitals with general
//  interaction kernels"
//  J. Chem. Theory Comput. 2020, https://doi.org/10.1021/acs.jctc.9b01296
//
// It evaluates the Boys function as described in the appendix of the main ESI
// of [1]. The ESI distribution of said article also includes the correponding
// data table generator RebuildBoysFnDataTables.py (for IrBoysFn.inl) which
// allows increasing the MaxM supported by this function (which may or may not
// be supplied the the file you see here).
//
// Please retain this note.
//
// - Gerald Knizia, 2019-12-28

#include <stdexcept>
#include <cmath>
#include "IrBoysFn.h"
#include "IrBoysFn.inl"

namespace ir {

static double const
   // First T which's nearest table index we will consider as lying outside the
   // range of our tablulated data. See comment at bottom regarding the -0.75.
   s_EndOfTableRangeT = (double(TableSize) - 0.75)/double(StepsPerT);

// store Boys function values Factor*Fm(T) at pOut[0...MaxM], MaxM inclusive
void IrBoysFn(double *pOut, double T, unsigned MaxM, double Factor)
{
   if (MaxM > TableMaxM)
      throw std::runtime_error("This version of IrBoysFn does not support the supplied M.");
//    if (iTab >= size_t(TableSize)) {
//    // FIXME: Needs fixing for *VERY* large T values, causing overflows in 'uint' types
   if (T >= s_EndOfTableRangeT) {
      // large T, out of tabulated range.
      // calculate F(m=0,T) directly and get the other Fms by upward recursion.
      // Note: F0(T) = sqrt(pi/4) * Erf(sqrt(T))/sqrt(T)
      //       F_{m}(T) = ((2*m-1)*F_{m-1}(T) - exp(-T))/(2*T)
      double
         ExpT = 0.,
         Erf = 1.,
         SqrtT = std::sqrt(T),
         InvSqrtT = 1./SqrtT;
      if (T < 36 + 2*MaxM)
         // may need this in the upward recursion
         ExpT = std::exp(-T);
      if (T < 36) { // <- for T >= 36 1-erf(T) < 1e-16.
         // Target: Erf = std::erf(SqrtT);
         // p is a polynomial approximation for erfc(x)*exp(x**2)*x
         // with x=T**.5 in the interval T=12..36
         double p =     -7.2445454685049530e-09;
         p = p * SqrtT + 3.8045328994203999e-07;
         p = p * SqrtT - 9.0663179364002847e-06;
         p = p * SqrtT + 1.2945938765877632e-04;
         p = p * SqrtT - 1.2313078185300554e-03;
         p = p * SqrtT + 8.1949930595442508e-03;
         p = p * SqrtT - 3.8962506318633815e-02;
         p = p * SqrtT + 1.3230181851498005e-01;
         p = p * SqrtT - 3.1344877124388754e-01;
         p = p * SqrtT + 4.8560557460742598e-01;
         p = p * SqrtT + 1.6057555044955896e-01;
         Erf = 1. - p*ExpT*InvSqrtT;
         // max error of fit for erf(sqrt(T)): 1.11e-16
      }

      double
         Inv2T = .5 * InvSqrtT * InvSqrtT; // 1/(2*T)
      ExpT *= Factor;
      pOut[0] = Factor * HalfSqrtPi * InvSqrtT * Erf;
      for (unsigned m = 1; m <= MaxM; ++ m)
         pOut[m] = ((2*m-1) * pOut[m-1] - ExpT)*Inv2T;
   } else {
//       size_t
//          // index of table entry closest to T
//          iTab = (size_t)(T * StepsPerT + .5);
      unsigned
         // index of table entry closest to T
         iTab = (unsigned)(T * StepsPerT + .5);
      // T is in tabulated range. Calculate Fm for highest order M by
      // Taylor expansion and other Ts by downward recursion.
      double const
         // difference to next tabulated point
         Delta = T - iTab*TsPerStep,
         *pTab = &BoysData[iTab*TableOrder + MaxM];
      double
         Fm;
      Fm =              pTab[8]/(double)(1*2*3*4*5*6*7*8);
      Fm = Delta * Fm - pTab[7]/(double)(1*2*3*4*5*6*7);
      Fm = Delta * Fm + pTab[6]/(double)(1*2*3*4*5*6);
      Fm = Delta * Fm - pTab[5]/(double)(1*2*3*4*5);
      Fm = Delta * Fm + pTab[4]/(double)(1*2*3*4);
      Fm = Delta * Fm - pTab[3]/(double)(1*2*3);
      Fm = Delta * Fm + pTab[2]/(double)(1*2);
      Fm = Delta * Fm - pTab[1];
      Fm = Delta * Fm + pTab[0];

      pOut[MaxM] = Factor * Fm;
      if (MaxM != 0) {
         double
            ExpT = Factor * std::exp(-T);
         for (int m = MaxM-1; m != -1; --m)
            pOut[m] = (2*T * pOut[m+1] + ExpT)*Inv2mp1[m];
      }
   }
//    printf("IrBoysFn(T=%f, M=%i, F=%f) -> [0] = %f", T, MaxM, Factor, pOut[0]);
}

} // namespace ir


// Comments on s_EndOfTableRangeT:
//
// - In exactly conforming implementations we could use -0.5 here instead of
//   -0.75 (in this case, any T-h with h > 0 would produce a valid table
//   index, and any T=s_EndOfTableRangeT or larger would be out of range).
//   But we leave a bit of buffer space for potentially 'unconventional'
//   truncation behavior on non-standard platforms, and for cases of
//   non-exactly representable table boundaries if not using StepsPerT
//   which are powers of 2.
//
// - We cannot just compute the table index and check against that. The reason
//   is that *very large* values of T could lead to an integer overflow, which
//   could potentially produce table indices inside the tabulated range despite
//   T clearly being out of said range.
