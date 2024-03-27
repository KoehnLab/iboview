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

#include "CxVec3.h"
#include "CxPoly.h"
#include "CxAlgebra.h"
#include "format.h"
#include <iostream>
#include <vector>
#include <valarray>

namespace ct {

typedef std::valarray<double>
   FValueVector0;
// typedef std::valarray<TVector3<double> >
//    FValueVector1;

// squared norm of a double...  in order to make it treatable as a vector.
double LengthSq(double f) {
   return f*f;
}


double LengthSq(FValueVector0 const &v) {
   double r = 0;
   for (size_t i = 0; i != v.size(); ++ i)
      r += v[i]*v[i];
   return r;
}


// double LengthSq(FValueVector1 const &v) {
//    double r = 0;
//    for (size_t i = 0; i != v.size(); ++ i)
//       r += LengthSq(v[i]);
//    return r;
// }


template<class FScalar>
void MakePolynomialArgumentPowers(FScalar *tk, FScalar t, size_t D)
{
   FScalar
      tk0 = FScalar(1);
   for (size_t iDeg = 0; iDeg != D; ++ iDeg) {
      tk[iDeg] = tk0;
      tk0 *= t;
   }
}


template<class FScalar>
FScalar GetPolynomialDerivativeArgumentPower(FScalar const *tiPow, size_t iDeg, int iDeriv)
{
   assert(iDeriv >= 0); // integral case not implemented (although it could be)
   if (iDeg == 0) {
      assert(iDeriv >= 0);
      if (iDeriv == 0)
         return FScalar(1);  // no derivatives -- return t^0 = 1.
      else
         return FScalar(0);  // (d/dt)^d c = 0 for d != 0
   }

   if (iDeriv == 0)
      // return ti^k
      return tiPow[iDeg];
   else
      // return (d/dt)^d t^k = (d/dt)^{d-1}  (d/dt) t^k =  (d/dt)^{d-1} k * t^(k-1)
      return iDeg * GetPolynomialDerivativeArgumentPower(tiPow, iDeg-1, iDeriv-1);
}




template<class FValue, class FScalar>
TPolynomialFit<FValue, FScalar>::TPolynomialFit(FValue const *pValues, int *pDerivDegrees, FScalar const *pArguments, FScalar const *pArgumentWeights, size_t nPoints, size_t Degree, FScalar ThrSigRel)
{
   size_t
      N = nPoints,
      D = Degree + 1; // fit dimension (higher than degree because we include the 0'th order constant polynomial)
   if (N == 0 || N < D)
      throw std::runtime_error(fmt::format("TPolynomialFit: number of input points (N={}) smaller than target polynomial degree (D={}).", nPoints, Degree));

   // find and store range of argument scalars
   tFirst = pArguments[0];
   tLast = tFirst;
   FScalar
      tMinNrmSq = LengthSq(pValues[0]);
   tMinGuess = tFirst;
   for (size_t iPt = 1; iPt < N; ++ iPt) {
      FScalar
         t = pArguments[iPt];
      tFirst = std::min(tFirst, t);
      tLast = std::max(tLast, t);
      // also remember input point of minimal norm.
      FScalar
         nsq = LengthSq(pValues[iPt]);
      if (nsq < tMinNrmSq) {
         tMinGuess = t;
         tMinNrmSq = nsq;
      }
   }
//    tMinGuess = .5*(tFirst + tLast); // <- will be used as expansion point for polynomial (tshift)

   FScalar
      // expansion point of new polynomial
      tShift = tMinGuess;
   std::vector<FValue>
      // coefficients of fitted polynomial (to be determined).
      NewCoeffs(D);
   Residual = 0.;

   {
//       // compute N x D matrix of argument powers.
//       TArray<FScalar>
//          A(N * D);
//       for (size_t iPt = 0; iPt != N; ++ iPt) {
//          FScalar
//             t = (pArguments[iPt] - tShift),
//             tk = FScalar(1);
//          for (size_t iDeg = 0; iDeg != D; ++ iDeg) {
//             A[iPt + N*iDeg] = tk;
//             tk *= t;
//          }
//       }
//       assert(A.size() == N * D);

      // compute N x D matrix of argument powers (possibly including derivatives).
      TArray<FScalar>
         A(N * D),
         tiPow(D);
      std::vector<FValue>
         RhsValues1(N); // IMPORTANT: don't make TArray<>s of FValue... FValue may be non-POD! (e.g., vector-valued)
      for (size_t iPt = 0; iPt != N; ++ iPt) {
         FScalar
            t = (pArguments[iPt] - tShift);
         MakePolynomialArgumentPowers(&tiPow[0], t, D);

         for (size_t iDeg = 0; iDeg != D; ++ iDeg) {
            int iDeriv = 0;
            if (pDerivDegrees)
               iDeriv = pDerivDegrees[iPt];
            A[iPt + N*iDeg] = GetPolynomialDerivativeArgumentPower(&tiPow[0], iDeg, iDeriv);
         }
         RhsValues1[iPt] = pValues[iPt];
      }
      assert(A.size() == N * D);

      // if supplied, scale both polynomial results and target rhs by argument weight factors for fit
      if (pArgumentWeights) {
         for (size_t iPt = 0; iPt != N; ++ iPt) {
            double w = pArgumentWeights[iPt];
            for (size_t iDeg = 0; iDeg != D; ++ iDeg) {
               A[iPt + N*iDeg] *= w;
            }
            RhsValues1[iPt] *= w;
         }
      }

      // compute A's singular value decomposition:
      //    A = U * diag(Sigma) * Vt
      size_t
         nSig = std::min(N, D);
      TArray<FScalar>
         U(N * nSig),
         Vt(nSig * D),
         Sigma(nSig);

      // ComputeSvd(double *pU, size_t ldU, double *pSigma, double *pVt, size_t ldVt, double *pInAndTmp, size_t ldIn, size_t nRows, size_t nCols);
      ComputeSvd(U.data(), N, Sigma.data(), Vt.data(), nSig, A.data(), N, N, D);

      // find largest singular value (should be last; singular values should come out ordered.)
      FScalar
         MaxSig = FScalar(0);
      for (size_t iSig = 0; iSig != nSig; ++ iSig)
         MaxSig = std::max(MaxSig, Sigma[iSig]);

      // behold some super-lame matrix-vector multiplications. Done this way for
      // the sake of generality with respect to FValue vector format.
      // What this actually does is:
      //
      //   C = Vt * diag(sigma**-1) * (U.T * rhs),
      //
      // where C and rhs (the target values) are vectors, and the rest
      // are matrices of scalars.
      {
         std::vector<FValue>
            UtRhs(nSig), // = dot(U.T, rhs)
            UtRhsInvSig(nSig); // = dot(U.T, rhs) * diag(Sigma ** -1)
         for (size_t iSig = 0; iSig != nSig; ++ iSig) {
            // compute dot(U.T, rhs) and its 1/Sigma scaled version.
            FValue
//                v = FValue(0);
               v = ZeroLike<FValue,FScalar>(pValues[0]);
            for (size_t iPt = 0; iPt != N; ++ iPt)
               v += U[iPt + N*iSig] * RhsValues1[iPt];
            UtRhs[iSig] = v;

            if (Sigma[iSig] != 0. && (ThrSigRel == 0. || (Sigma[iSig] >= ThrSigRel * MaxSig))) {
               UtRhsInvSig[iSig] = (FScalar(1)/Sigma[iSig]) * UtRhs[iSig];
            } else {
               // skip contribution of current singular value in case of dependency.
//                UtRhsInvSig[iSig] = FValue(0);
               UtRhsInvSig[iSig] = ZeroLike<FValue,FScalar>(pValues[0]);
            }
         }

         // compute output coefficients from UtRhsInvSig and Vt.
         for (size_t iDeg = 0; iDeg != D; ++ iDeg) {
            FValue
//                ck = FValue(0);
               ck = ZeroLike<FValue,FScalar>(pValues[0]);
            for (size_t iSig = 0; iSig != nSig; ++ iSig) {
               ck += Vt[iSig + nSig*iDeg] * UtRhsInvSig[iSig];
            }
            NewCoeffs[iDeg] = ck;
         }

      }
   }

   p = FPolynomial(NewCoeffs.data(), NewCoeffs.size(), tShift);

   // compute rms residual
   FScalar
      fResidualSq(0);
   for (size_t iPt = 0; iPt != N; ++ iPt) {
      FValue
         // value of polynomial fit at position t[iPt]
         pi = p(pArguments[iPt]),

         // ...and its difference from the specified target value.
         delta = pi - pValues[iPt];
      if (pArgumentWeights)
         delta *= pArgumentWeights[iPt];
      fResidualSq += LengthSq(delta);
   }

   // compute rms residual --- normalized with respect to number of points N.
   Residual = std::sqrt(fResidualSq/FScalar(N));
}


// explicitly instanciate for likely target vector types.
template struct TPolynomialFit<double, double>;
template struct TPolynomialFit<TVector3<double>, double>;
template struct TPolynomialFit<FValueVector0, double>;
// template struct TPolynomialFit<FValueVector1, double>;



template<class FValue, class FScalar>
TPolynomialMinimizeParams<FValue, FScalar>::TPolynomialMinimizeParams(FScalar tFirst_, FScalar tLast_, std::ostream *pout_)
   : tFirst(tFirst_), tLast(tLast_), pFactor(+1.0), ThrGrad(1e-10), ThrArg(1e-10), MaxStep(-1), MaxIt(2048), pout(pout_)
{
   if (tFirst > tLast)
      std::swap(tFirst, tLast);
   MaxStep = 5e-2 * (tLast - tFirst);

   Print = true;
   CheckSanity();
}


template<class FValue, class FScalar>
void TPolynomialMinimizeParams<FValue, FScalar>::CheckSanity() const
{
   if (tFirst >= tLast)
      throw std::runtime_error(fmt::format("TPolynomialMinimizeParams: tFirst not smaller than tLast:  tFirst = {}   tLast = {}.", tFirst, tLast));
   if (MaxStep <= 0)
      throw std::runtime_error(fmt::format("TPolynomialMinimizeParams: MaxStep = {} is invalid. Must be > 0.", MaxStep));
}


template<class FValue, class FScalar>
FValue FindMinimum(TPolynomial<FValue, FScalar> const &p_, FValue tGuess, TPolynomialMinimizeParams<FValue, FScalar> const &Params)
{
   Params.CheckSanity();
   if (p_.Degree() < 2)
      // note: for t\in [tMin, tMax] constrained searches, technically linear
      // polynomials would be mathematically fine, too (this would result in
      // t=tMin or t=tMax as output). But the code currently does not support
      // that, because we do not need it at this moment.
      throw std::runtime_error(fmt::format("FindMinimum(): input polynomial has degree {}; currently degree >= 2 is required.", p_.Degree()));
   typedef TPolynomial<FValue, FScalar>
      FPoly1;
   std::ostream
      *pout = Params.pout;
   if (pout == 0)
      pout = &std::cout;

   FPoly1
      p = Params.pFactor * p_,
      Grad = p.GetDerivativePolynomial(),
      Hess = Grad.GetDerivativePolynomial();

   // Newton search of minimum.
   if (Params.Print && pout != 0) {
//       (*pout) << fmt::format(" {:<30}MAXSTP = {:.2e}  THRGRAD = {:.2e}\n\n", "Stationary point search:", Params.MaxStep, Params.ThrGrad);
      (*pout) << fmt::format(" {:<30}MAXSTP = {:.2e}  THRGRAD = {:.2e}  FACTOR = {:+.2f}\n\n", "Stationary point search:", Params.MaxStep, Params.ThrGrad, Params.pFactor);
//       (*pout) << fmt::format("   ITER.   COORD.VALUE/t  POLY.VALUE/p(t)       STEP     GRADIENT") << std::endl;
      (*pout) << fmt::format("   ITER.   POLY.VAL/p(t)    PARAMETER/t         STEP     GRADIENT") << std::endl;
   }

   // FIXME: cgk 2019-06-29:
   //  - This can, apparently, end up in an infinite loop, with
   //    alternating +dx/-dx, both of which are exactly equal to MaxStep:
   //      ctrace: bracketed stationary point #4: E[-2] = -26.8785095122  E[-1] = -26.8785122012  E[0] = -26.8785041296.
   //      ctrace: computing polynomial fit of E(coord) curve and its extremal geometry.
   //      Stationary point search:      MAXSTP = 2.50e-03  THRGRAD = 1.00e-10
   //
   //         ITER.   POLY.VAL/p(t)    PARAMETER/t         STEP     GRADIENT
   //           1      -26.87851220     3.76454632    -0.00250000   4.36e-05
   //           2      -26.87851231     3.76204632     0.00250000   4.44e-05
   //           3      -26.87851220     3.76454632    -0.00250000   4.36e-05
   //           4      -26.87851231     3.76204632     0.00250000   4.44e-05
   //           5      -26.87851220     3.76454632    -0.00250000   4.36e-05
   //           6      -26.87851231     3.76204632     0.00250000   4.44e-05
   //           7      -26.87851220     3.76454632    -0.00250000   4.36e-05
   //           8      -26.87851231     3.76204632     0.00250000   4.44e-05
   //           9      -26.87851220     3.76454632    -0.00250000   4.36e-05
   //        <...>
   //        2047      -26.87851220     3.76454632    -0.00250000   4.36e-05
   //        2048      -26.87851231     3.76204632     0.00250000   4.44e-05
   //
   //        /// ERROR OCCURED -- ENTERED OUTER EXCEPTION HANDLER ///
   //
   //  - Encountered in /home/cgk/calc/organic_reactions/reactions/aldol-addition-simple
   //    in a constrained ctrace reaction path optimization.
   //  - Additional info on actual numbers in gradient/hessian/step:
   //
   //            ITER.   POLY.VAL/p(t)    PARAMETER/t         STEP     GRADIENT
   //              1      -26.87492792     5.03954632    -0.00250000   g(t)=-1.1176e-04  h(t)=-1.1908e-04  (-g/h)=-9.3854e-01
   //              2      -26.87492764     5.03704632    -0.00250000   g(t)=-1.0978e-04  h(t)=-1.3275e-03  (-g/h)=-8.2697e-02
   //              3      -26.87492737     5.03454632    -0.00250000   g(t)=-1.0581e-04  h(t)=-1.7149e-03  (-g/h)=-6.1699e-02
   //              4      -26.87492711     5.03204632    -0.00250000   g(t)=-1.0190e-04  h(t)=-1.2682e-03  (-g/h)=-8.0351e-02
   //              5      -26.87492686     5.02954632     0.00250000   g(t)=-1.0017e-04  h(t)= 2.5434e-05  (-g/h)= 3.9384e+00
   //              6      -26.87492711     5.03204632    -0.00250000   g(t)=-1.0190e-04  h(t)=-1.2682e-03  (-g/h)=-8.0351e-02
   //              7      -26.87492686     5.02954632     0.00250000   g(t)=-1.0017e-04  h(t)= 2.5434e-05  (-g/h)= 3.9384e+00
   //              8      -26.87492711     5.03204632    -0.00250000   g(t)=-1.0190e-04  h(t)=-1.2682e-03  (-g/h)=-8.0351e-02
   //              9      -26.87492686     5.02954632     0.00250000   g(t)=-1.0017e-04  h(t)= 2.5434e-05  (-g/h)= 3.9384e+00
   //             10      -26.87492711     5.03204632    -0.00250000   g(t)=-1.0190e-04  h(t)=-1.2682e-03  (-g/h)=-8.0351e-02
   //
   //  - ...so we get some negative Hessian values, which invert the search direction.
   //  - Currently hacked around with the Hessian fallback variants inserted
   //    below (using a positive reference fh obtained from 2nd order polynomial
   //    coefficient in case we encounter any negative Hessian values along the way)

   {
      FScalar
         tMin = tGuess;
      FScalar
         // will be inserted if actual hessian at intermediate point turns out
         // negative
         fh_Fallback = Hess(tGuess),
         coeff_p2 = p.GetCoeff(2); // poly coeff in front of (t-ArgShift)^2
      if (fh_Fallback < std::abs(0.1*coeff_p2))
         fh_Fallback = std::abs(0.1*coeff_p2);
//       if (fh_Fallback < 0.)
//          throw std::runtime_error(fmt::format("FindMinimum(): Hessian at reference point is not positive (in polynomial line search): t = {:.6f}  g = {:.6f}  h = {:.6f}", tGuess, Grad(tGuess), Hess(tGuess));

      bool
         Converged = false;
      for (size_t iIt = 0; iIt != Params.MaxIt; ++ iIt) {
         FScalar
            fg = Grad(tMin),
            fh = Hess(tMin),
            fh_for_denom = fh;
         // ^- formally these would be vectors, but we really only need
         // the scalar variants for the current algorithm.
         if (fh_for_denom <= 0.)
            fh_for_denom = fh_Fallback;
         FScalar
            // this is the regular Newton step, unless fh has been
            // replaced by the fallback value.
            dx = -fg/fh_for_denom;

         if (std::abs(dx) > Params.MaxStep)
            dx = Params.MaxStep * ((dx > 0)? FScalar(1) : FScalar(-1));

         if (Params.Print && pout != 0) {
//             (*pout) << fmt::format("{:>6}. ->   t ={:14.8f}    g = {:9.2e}    p(t) ={:14.8f}", iIt+1, tMin, fg, p(tMin)) << std::endl;
//             (*pout) << fmt::format("{:6}{:>18.8f}{:15.8f}{:15.8f}{:11.2e}", (iIt+1), p(tMin)/Params.pFactor, tMin, dx, fg) << std::endl;
            (*pout) << fmt::format("{:6}{:>18.8f}{:15.8f}{:15.8f}   g(t)={:11.4e}  h(t)={:11.4e}  (-g/h)={:11.4e}", (iIt+1), p(tMin)/Params.pFactor, tMin, dx, fg, fh, -fg/fh) << std::endl;
         }
         if (std::abs(fg) < Params.ThrGrad) {
            Converged = true;
            break;
         }
         tMin += dx;
         if (tMin < Params.tFirst)
            tMin = Params.tFirst;
         if (tMin > Params.tLast)
            tMin = Params.tLast;
      }
      if (!Converged)
         throw std::runtime_error(fmt::format("FindMinimum(): Failed to converge after {} iterations (in polynomial line search).", Params.MaxIt));

      if (Params.Print && pout != 0) {
         (*pout) << "\n";
      }
      return tMin;
   }
}


template struct TPolynomialMinimizeParams<double, double>;
template double FindMinimum<double, double>(TPolynomial<double, double> const &, double tGuess, TPolynomialMinimizeParams<double, double> const &Params);


// template struct TPolynomialMinimizeParams<FValueVector0, double>;
// template double FindMinimum<FValueVector0, double>(TPolynomial<FValueVector0, double> const &, double tGuess, TPolynomialMinimizeParams<FValueVector0, double> const &Params);



} // namespace ct
