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

#include <cmath> // for std::expm1 (well... technically a C++11 routine...)
#include <algorithm> // for std::swap
#ifndef IR_ECP_BESSELFN_TEST
   #include "CxDefs.h" // for assert only. You can turn it off if you don't have it.
#else
   #include <cassert>
#endif

using std::size_t;
using std::ptrdiff_t;

namespace ir {

static double sqr(double x) { return x*x; }

// static double const _ThrReallySmallZ = 1e-8;
// static double const _ThrReallySmallZ = 1e-3;
//  ^- either one works (1e-3 or 1e-8) for <= 5e-15 errors with the 4th order
//     series (max error does not come from here). I put in two more terms now,
//     and with these 1e-2 also is fine:
static double const _ThrReallySmallZ = 1e-2;



// compute ratio {\bar i}_{m}(z)/{\bar i}_{m-1}(z) for
// {\bar i}_m(z) = exp(-z) i_m(z) z^{-m} in the LARGE z expansion
//
// (NOTE: for 2z >= 37 -> z > 18.4, it is within 1e-16 precision for the m=0
// component... for larger m, the situation is rather more dicy, however...,
// which is why for many of those we still use rm12_z_down. The comments in
// Flores-Moreno which says you can just use this for z > 16 are definitely
// incorrect! (they cite that for the R-formula, but the R formula has exactly
// the same problem as this one).
static double rm12_iz_up(ptrdiff_t M, double iz)
{
   double
      // that's "s_{-1}" = g_{-1} + g_{0} = 0 + iz = iz (to start the recurrence)
      s_nm1 = iz,
      // that's s_0 = g_0 + g_{-1} = iz + 0 = iz
      s_n = iz;
   for (ptrdiff_t n = 1; n <= M; ++ n) {
      // upward recurrence: this one is a transform of Exp[-z]-times the
      // cosh/sinh formulas on the mathworld article mentioned below, under the
      // assumption that Exp[-2z] = 0)
      // See bessel_algo_test.py for details on how this works.
      s_nm1 += (-2*n + 1)*iz*s_n;
      // ^- that's really the new s_n. Beware of unsigned integers!
      std::swap(s_nm1, s_n);
   }
   return iz*s_n/s_nm1;
}


// compute ratio {\bar i}_{m}(z)/{\bar i}_{m-1}(z) for
// {\bar i}_m(z) = exp(-z) i_m(z) z^{-m} for SMALL and intermediate z.
//
// Notes:
// - This function is numerically stable and works fine also for large z. But
//   it can become really slow in this case.
// - So we only use it when rm12_iz_up is not trustworthy (small Z/MaxM).
//   This includes many cases with larger MaxM in which conventional wisdom
//   would (incorrectly) claim that this is not required.
static double rm12_z_down(ptrdiff_t M, double z)
{
   //double const TargetAcc = 1e-16
   double const TargetAcc = 1e-6;
   // ^- that actually works for below 1e-15 accuracy. Technically, it does
   //    not *have* to do so, but in practice is was fine in my full scan
   //    (z = Exp10Spacing(1e-8, 7e2, 2100), MaxM in range(0,30+1))

   // Estimate starting point for downward recurrence:
   // accumulate (n+m)!/m! * ((2/z)^n) = (n+m)!/m! * (2*iz)^n
   // (that's the inverse of the *direct* coefficient of i_{m+n}(z) in i_m(z).
   //  see explanations in log_m_factorial_pow_2iz_m and EstimateStartingM
   //  in the full bessel_algo_test.py regarding how this works)
   // and see at which 'n' it becomes large enough that neglecting i_{m+n}(z)
   // terms becomes unimportant to reach the target accuracy.
   ptrdiff_t
      n = 0;
   double
      AccNpM = 1., // <- technically an integer, but it could get large.
      AccZh = 1.;
   for ( ; ; ) {
      n += 1;
      AccNpM *= (n + M);
      AccZh *= .5*z;
      if (AccZh <= TargetAcc * AccNpM)
         break;
   }
   ptrdiff_t
      StartingM = M + n - 1;
   double
      // current estimates of i_m(z) and i_{m+1}(z)
      im = 1.,
      // ^- Setting this to 1 is okay because we only need a ratio in the end.
      imp1 = 0.;
   bool
      UseRatios = sqr(z) > 1e-3;
   assert(M >= 2);
   for (ptrdiff_t m = StartingM; m >= M-1; --m) {
      // ^- beware of overflow at M=1. That's why I put all signed types here.
      //    (should only ever be called with M >= 2, but who knows...)

      // Recurrence here:
      //
      //  i_{m}(z) = (2*m + 3)/z * i_{m+1}(z) + i_{m+2}(z)
      //
      // (that's https://mathworld.wolfram.com/ModifiedSphericalBesselFunctionoftheFirstKind.html's
      // recurrence relation re-indexed and solved for downward recurrence)
      imp1 = (2*m + 3)*im + sqr(z) * imp1;
      // ^- that's really the new i_m(z)
      std::swap(im, imp1);

      if (UseRatios) {
         // for very large z, re-normalize the imp1 and im during the iterations
         // --- otherwise we could run into floating point overflows because
         // *both* im and imp1 become too large.
         im /= imp1;
         imp1 = 1.;
      }
   }
   return imp1/im;
}


// compute ratio {\bar i}_{m}(z)/{\bar i}_{m-1}(z) using either the
// large-z or the small-and-intermediate-z routines to do so, depending
// on what is appropriate given the concrete set of parameters z, m.
// If != 0, iz must be the inverse of z (iz = 1/z).
static double rm12_gen(ptrdiff_t m, double z, double iz) {
   // note: "iz = 0" means "no iz, use small-z version"
   assert(iz != 0 || z <= _ThrReallySmallZ);
   if ((iz != 0.) && ((m + 2)*iz < 0.25)) {
                     // ^- m+2: otherwise there are cases with
                     //         substantial errors for MaxM=1.
      // large-z version
      return rm12_iz_up(m, iz);
   } else {
      // small and intermediate-z version
      return rm12_z_down(m, z);
   }
}


// make scaled output integrals of scaled spherical bessel functions i_m(z)
// with z^{-m} s^{m} envelopes absorbed in them:
//
//    pOut[m] = C * (i_m(z)*z^{-m}) * s^{m}
//
// for m = {0,1,...,MaxM} (MaxM inclusive). C is a undefined, but common, scale.
// This variant is for *SMALL* and moderate z.
static void MakeDownward_C_in_izn_SmallZ_v2(double *pOut, size_t MaxM, double z, double s, double ratio_m_plus_2_by_m_plus_1)
{
   double
      imp1 = ratio_m_plus_2_by_m_plus_1,
      im = 1;
   for (size_t m = MaxM; m <= MaxM; --m) {
      imp1 = ((2*m + 3)*im + sqr(z)*imp1);
      // ^- that's really the new {\bar i}_m(z)
      std::swap(im, imp1);
      pOut[m] = im;
   }
   if (s != 1.) {
      double
         acc_s = 1.;
      for (size_t m = 1; m <= MaxM; ++m) {
         acc_s *= s;
         pOut[m] *= acc_s;
      }
   }
}

// as MakeDownward_C_in_izn_SmallZ, but this is the variant for *LARGE* z.
static void MakeDownward_C_in_izn_LargeZ_v2(double *pOut, size_t MaxM, double iz, double s, double ratio_m_plus_2_by_m_plus_1)
{
   double
      imp1 = ratio_m_plus_2_by_m_plus_1,
      im = iz;
   for (size_t m = MaxM; m <= MaxM; --m) {
      imp1 += (2*m + 3)*im*iz;
      // ^- that's really the new {\bar i}_m(z)
      std::swap(im, imp1);
      pOut[m] = im;
   }
   double
      inv_a = iz*s,
      acc_iz = 1.;
   for (size_t m = 1; m <= MaxM; ++m) {
      acc_iz *= inv_a;
      pOut[m] *= acc_iz;
   }
}


/// Compute modified spherical Bessel functions of the first kind, i_m(z), with
/// envelopes. Concretely, compute for m in [0,1,...,MaxM] (MaxM inclusive!) the
/// functions
///
///       pOut[m] = i_{m}(z) * Exp[-z] * z^{-m} * s^m * Factor
///
/// where z := s * a (if you do not want a separate 's', pass a=z and s=1).
/// 'Factor' is an optional common factor multiplied into all output values.
///
/// This is a C++ port of ~/dev/ir/ecp_notes/bessel_algo_test_minimal.py.
/// Technical documentation regarding the algorithms is inside there.
void EvalModifiedSphericalBesselEnveloped(double *pOut, double s, double a, size_t MaxM, double Factor)
{
   double
      z = s*a;
   // note: z may be zero (if a is), but s is not allowed to be zero.
   assert(s != 0.);
   double
      iz, // inverse z (iz = 1/z) if available.
      OutM0;
   if (z >= _ThrReallySmallZ) {
      iz = 1./z;
      OutM0 = -.5*iz*std::expm1(-2*z);
      // ^- with expm1 this should be good for all z, including small ones
      // (as long as they are not exactly zero...).
   } else {
      // For very small input arguments 'z' we do not actually have to evaluate
      // the exponential function, but can just insert its small-z asymptotic expansion.
      // Also has the advantage of allowing evaluation of expm1(x)/x at x = 0 exactly.
      //
      // Here we use x=-2z. For this Wolfram Alpha says:
      //
      //     Series[(exp(-2*z)-1)/(-2*z),{z,0,6}]
      //     = 1 - z + (2/3) z^2 - (1/3) z^3 + (2/15) z^4 - (2/45) z^5 + (4/315) z^6 + O(z^7).
      //
      // Error introduced due to this should be O(z^7) for z < _ThrReallySmallZ
      // with the 6th order expansion
      OutM0 = 1. - z*(1. - z*((2./3.) - z*((1./3.) - z*((2./15.) - z*((2./45.) - z*(4./315.))))));
      // small z case. No iz here.
      iz = 0.;
   }

   // absorb external prefactor; all outputs will be scaled by this.
   OutM0 *= Factor;

   if (MaxM == 0) {
      pOut[0] = OutM0;
   } else {
      double ratio_m_plus_2_by_m_plus_1 = rm12_gen(ptrdiff_t(MaxM)+2, z, iz);
      if (z < 1e2) {
         MakeDownward_C_in_izn_SmallZ_v2(pOut, MaxM, z, s, ratio_m_plus_2_by_m_plus_1);
      } else {
         MakeDownward_C_in_izn_LargeZ_v2(pOut, MaxM, iz, s, ratio_m_plus_2_by_m_plus_1);
      }
      // the integrals in pOut at this point have a common, but undefined, prefactor
      // C absorbed. They currently are:
      //
      // pOut[m] = C * i_m(z) * z^{-m} * s^{m}
      //
      // Now replace this C by our target envelope function (Exp[-2z]-1)/(-2z).
      double
         f = OutM0/pOut[0];
      for (size_t m = 0; m <= MaxM; ++ m)
         pOut[m] *= f;
   }
}

} // namespace ir


#ifdef IR_ECP_BESSELFN_TEST
/*******************************************************************************
  TESTING THIS FILE:

- Debug version:

    c++ -I.. -g -Wall -Wpedantic -DIR_ECP_BESSELFN_TEST IrEcpBesselFn.cpp ../format.cpp -o /tmp/sph-bessel-test && /tmp/sph-bessel-test

- Version with full optimizations:

    c++ -I.. -O3 -march=native -DNDEBUG -Wall -Wpedantic -DIR_ECP_BESSELFN_TEST IrEcpBesselFn.cpp ../format.cpp -o /tmp/sph-bessel-test && /tmp/sph-bessel-test

*******************************************************************************/

#include <iostream>
#include <valarray>
#include <stdexcept>
#include "format.h"


namespace ir_bessel_fn_test {
   typedef std::valarray<double>
      FArrayType;


   static void _PrintTab(FArrayType const &RefData, FArrayType const &TrialData, FArrayType const &rDiff)
   {
      if (RefData.size() != TrialData.size() || TrialData.size() != rDiff.size())
         throw std::runtime_error("  IrEcpBesselFn:_PrintTab(): Input data arrays inconsistent in testing print out.");
      std::cout << "  Test values:" << std::endl;
      for (unsigned m = 0; m < RefData.size(); ++ m) {
         std::cout << fmt::format("     m ={:3}   ref = {:16.12e}   trial = {:16.12e}   rel-err = {:8.2e}", m, RefData[m], TrialData[m], rDiff[m]) << std::endl;
      }
   }

   static void _TestAgainstRef(double TrialS, double TrialA, unsigned MaxM, FArrayType const &RefData)
   {
      FArrayType
         TrialData(-1., MaxM+1); // note: MaxM inclusive.
      if (RefData.size() != MaxM+1)
         throw std::runtime_error(fmt::format("_TestAgainstRef: ref data has incorrect format. Expected {} reference entries (for MaxM={}), but found {}.", MaxM+1, MaxM, RefData.size()));

      assert(TrialData.size() == MaxM+1);
      ir::EvalModifiedSphericalBesselEnveloped(&TrialData[0], TrialS, TrialA, MaxM, 1.0);

      FArrayType
         rDiff = (TrialData - RefData)/(RefData + 1e-300); // <- +1e-300: otherwise zero division for exact 0.0 ref data.
                                                           //    (yes, I also think that is a very neat and elegant solution.
                                                           //     Thank you very much.)
      double
         fMaxErr = std::abs(rDiff).max();
      std::cout << fmt::format("  IrEcpBesselFn: Testing:  MaxM = {:<2}  |  s = {:<16.8}  |  a = {:<16.8}", MaxM, TrialS, TrialA);
      std::cout.flush();
      if (!(fMaxErr < 1e-12)) {
         std::cout << fmt::format("  ...FAILED  (fMaxErr={:.2e})", fMaxErr) << std::endl;
         _PrintTab(RefData, TrialData, rDiff);

         throw std::runtime_error(fmt::format("_TestAgainstRef: Test for s = {} a = {} MaxM = {} FAILED. fMaxErr(Relative) = {:.2e}", TrialS, TrialA, MaxM, fMaxErr));
      } else {
         std::cout << fmt::format("  ...passed  (fMaxErr={:.2e})", fMaxErr) << std::endl;
         // _PrintTab(RefData, TrialData, rDiff);
      }
   }

   // Reference data and the actual test calls and setup are in the .inl file.
   // It is generated by "python bessel_algo_test_minimal.py ref"
   #include "IrEcpBesselFnTests.inl"
} // namespace ir_bessel_fn_test


int main(int argc, char const **argv)
{
   ir_bessel_fn_test::_RunAllTests();
}

#endif // IR_ECP_BESSELFN_TEST
