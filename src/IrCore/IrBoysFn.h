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

#ifndef IR_BOYSFN_H
#define IR_BOYSFN_H

namespace ir {

   /// Store Boys function values Factor*Fm(T) at pOut[0...MaxM], MaxM inclusive
   ///
   /// Computes values of the Boys function
   ///
   /// F_m(T) = (-d/dT)^m (1/2) sqrt(pi/T) erf(sqrt(T))
   ///
   /// at a given T for a range of m = 0,1,..., MaxM. MaxM is *inclusive*.
   ///
   /// @param[out] pOut Output buffer of at least (MaxM+1) length.
   ///   On return, pOut[m] will contain the value (Factor*F_m(T)).
   /// @param[in] T Scalar argument (typically rho*R^2) at which to evaluate F_m(T)
   /// @param[in] MaxM Highest derivative 'm' to evaluate (inclusive!).
   /// @param[in] Factor A constant scalar factor which will be multiplied into
   ///   all returned output values (e.g., if you supply (2*M_PI)/rho as Factor,
   ///   the on return you will get pOut[m] = (2*M_PI)/rho*F_m(T) in the output
   ///   array.
   void IrBoysFn(double *pOut, double T, unsigned MaxM, double Factor);

} // namespace ir

#endif // IR_BOYSFN_H
