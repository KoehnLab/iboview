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

#ifndef CX_POLY_FIT_H
#define CX_POLY_FIT_H

#include <stdexcept>
#include <iosfwd>
#include <vector>
#include "CxTypes.h"
#include "CxPodArray.h"

namespace ct {

// represents a polynomial, mapping a scalar t of type FScalar to
//    p(t) = \sum_{k=0}^N c_k (t - tshift)^k
// where c_k is of class FValue (some sort of vector).
template<class FValue, class FScalar>
struct TPolynomial : public FIntrusivePtrDest1
{
   typedef FValue
      FValueType;
   typedef FScalar
      FScalarType;

   typedef TPolynomial<FValue, FScalar>
      FPolynomial;
   // represents a 0-polynomial.
   TPolynomial() { AssignCoeffs(0,0,FScalar(0)); }
   explicit TPolynomial(FValue const *pCoeffs_, size_t nCoeffs_, FScalar ArgShift_) { AssignCoeffs(pCoeffs_, nCoeffs_, ArgShift_); }


   // return highest monomial power in polynomial
   size_t Degree() const { return m_Coeffs.size() - 1; }

   FValue Eval(FScalar t) const;
   FValue operator () (FScalar t) const { return Eval(t); }

   // return a new polynomial object which computes the first derivative
   // of the current polynomial at point t.
   TPolynomial GetDerivativePolynomial() const;

   void operator *= (FScalar f);
   void operator += (FPolynomial const &other);

   FScalar GetArgShift() const { return m_ArgShift; }
   FValue const &GetCoeff(size_t i) const { assert(i < m_Coeffs.size()); return m_Coeffs[i]; }
protected:
   FScalar
      m_ArgShift; // will be subtracted from input 't's to make (t - tshift)^k
   typedef std::vector<FValue>
      FCoeffArray;
   FCoeffArray
      m_Coeffs;
   void AssignCoeffs(FValue const *pCoeffs_, size_t nCoeffs_, FScalar ArgShift_);
};



template<class FValue, class FScalar>
void TPolynomial<FValue,FScalar>::operator *= (FScalar f)
{
   for (size_t i = 0; i != m_Coeffs.size(); ++ i)
      m_Coeffs[i] *= f;
}


template<class FValue, class FScalar>
void TPolynomial<FValue,FScalar>::operator += (FPolynomial const &other)
{
   if (other.m_Coeffs.size() > m_Coeffs.size()) {
      m_Coeffs.resize(m_Coeffs.size(), FValueType(0));
   }
   for (size_t i = 0; i != other.m_Coeffs.size(); ++ i) {
      m_Coeffs[i] += other.m_Coeffs[i];
   }
}


template<class FValue, class FScalar>
void TPolynomial<FValue,FScalar>::AssignCoeffs(FValue const *pCoeffs_, size_t nCoeffs_, FScalar ArgShift_)
{
   m_ArgShift = ArgShift_;
   if (nCoeffs_ != 0) {
      m_Coeffs.assign(pCoeffs_, pCoeffs_ + nCoeffs_);
   } else {
//       // deal explicitly with a zero'th order polynomial returning 0.
//       // We could represent it with an empty coefficient array, but
//       // to avoid handling special cases in some instances, we just
//       // make a one-entry constant term evaluating to 0.
//       FValue Zero = FValue(0);
//       m_Coeffs.assign(&Zero, &Zero + 1);
      m_Coeffs.clear();
   }
}


template<class FValue, class FScalar>
TPolynomial<FValue,FScalar> operator * (FScalar f, TPolynomial<FValue,FScalar> const &p)
{
   TPolynomial<FValue,FScalar>
      out(p);
   out *= f;
   return out;
}


template<class FValue, class FScalar>
TPolynomial<FValue,FScalar> operator * (TPolynomial<FValue,FScalar> const &p, FScalar f) {
   return f * p;
}


template<class FValue, class FScalar>
FValue ZeroLike(FValue const &v) {
   return FScalar(0) * v;
}


template<class FValue, class FScalar>
TPolynomial<FValue,FScalar> operator + (TPolynomial<FValue,FScalar> const &a, TPolynomial<FValue,FScalar> const &b)
{
   TPolynomial<FValue,FScalar>
      out;
   if (a.size() >= b.size()) {
      out = a;
      out += b;
   } else {
      out = b;
      out += a;
   }
   return out;
}



template<class FValue, class FScalar>
FValue TPolynomial<FValue, FScalar>::Eval(FScalar t_) const
{
   if (m_Coeffs.empty())
      return FValue(0);
   FScalar
      t = (t_ - m_ArgShift),
      tk = FScalar(1.); // (t - tshift)^k
//    FValue
//       r = FValue(0); // value of the polynomial at point t
   FValue
      // value of the polynomial at point t
      r = ZeroLike<FValue,FScalar>(m_Coeffs[0]);
   for (size_t k = 0; k != m_Coeffs.size(); ++ k) {
      r += tk * m_Coeffs[k];
      tk *= t;
   }
   return r;
}


// return a new polynomial object computing the first derivative
// of the current polynomial at point t.
template<class FValue, class FScalar>
TPolynomial<FValue, FScalar> TPolynomial<FValue, FScalar>::GetDerivativePolynomial() const
{
   if (m_Coeffs.size() <= 1) {
      // constant polynomials have zero derivatives.
      return FPolynomial(0,0,FScalar(0));
   }
   // make coefficient array of new polynomial.
   // We just map \sum_{k=0}^D c_k (t-tshift)^k  to \sum_{k=1}^D c_k (k (t-tshift)^{k-1})
   // (so derivative polynomial has the same expansion point tshift)
   FCoeffArray
      NewCoeffs;
   NewCoeffs.reserve(m_Coeffs.size() - 1);
   for (size_t k = 1; k != m_Coeffs.size(); ++ k) {
      NewCoeffs.push_back(FScalar(k) * m_Coeffs[k]);
   }
   assert(NewCoeffs.size() == m_Coeffs.size() - 1);
   return FPolynomial(NewCoeffs.data(), NewCoeffs.size(), m_ArgShift);
}



template<class FVector, class FScalar>
// ct::TArray<FVector> LinSpace(FVector const &v0, FVector const &v1, size_t N)
std::vector<FVector> LinSpace(FVector const &v0, FVector const &v1, size_t N)
{
   if (N < 2)
      throw std::runtime_error("LinSpace: N must be >= 2.");
//    ct::TArray<FVector>
   std::vector<FVector>
      r;
   r.reserve(N);
   FScalar
      Nm1 = FScalar(1)/FScalar(N-1);
   for (size_t k = 0; k != N; ++ k) {
      FScalar f = k * Nm1;
      r.push_back((FScalar(1) - f)*v0 + f * v1);
   }
   return r;
}


// computes a least-squares polynomial fit of given degree to a given set of input points.
// todo: make this deal with shifts for fitting around extrema?
template<class FValue, class FScalar>
struct TPolynomialFit : public FIntrusivePtrDest1
{
   typedef TPolynomial<FValue, FScalar>
      FPolynomial;
   typedef typename FPolynomial::FValueType
      FValueType;
   FPolynomial
      // polynomial fitted to given input points
      p;
   FScalar
      // combined root-mean-square residual of least squares fit.
      Residual;
   FScalar
      // smallest and largest values of polynomial arguments 't'.
      // For tFirst <= t <= tLast the polynomial fit is an interpolation.
      // For 't' outside this range, it is an extrapolation.
      tFirst, tLast,
      // argument t = pArguments[iMin] for which the norm of the input vector norm(pValues[iMin]) is minimal.
      tMinGuess;

   // arguments:
   //    pValues[i] -> (d/dt)^d v(t_i) = pValues[i]  (vector-valued)
   //    pDerivDegrees[i] = d in (d/dt)^d. pDerivDegrees itself is allowed to be 0-pointer.
   //          In this case values v(t) themselves are being fitted, not derivatives of them (equivalent to supplying an array of 0s)
   //    pArguments[i] = t_i   (scalar arguments t in curve \vec v(t) for which reference values are given)
   //    pArgumentWeights[i]:  weight of t[i] in lstsq fit. pArgumentWeights may be 0; in this case all points are weighted equally.
   //    Degree: Degree of polynomial \vec v(t) = \sum_k \vec c_k (t-t0)^k to fit.
   explicit TPolynomialFit(FValue const *pValues, int *pDerivDegrees, FScalar const *pArguments, FScalar const *pArgumentWeights, size_t nPoints, size_t Degree, FScalar ThrSigRel = 0);
};

template<class FValue, class FScalar>
struct TPolynomialMinimizeParams
{
   FScalar
      // if set, restrict search to this range.
      tFirst, tLast;
   FScalar
      // factor to apply to p(t) before minimization.
      // Will typically be either +1 (for minimum search) or -1 (for maximum search).
      pFactor;
   FScalar
      ThrGrad, ThrArg,
      // maximum step size dt in single step of optimization. Newton steps t := t + dt
      // will be clamped to dt \in [-MaxStep,+MaxStep] interval if the original
      // dt = -f'(t)/f''(t) should exceed this interval.
      MaxStep;
   size_t
      MaxIt;
   bool
      Print;
   std::ostream
      *pout;
   explicit TPolynomialMinimizeParams(FScalar tFirst_, FScalar tLast_, std::ostream *pout_ = 0);
   void CheckSanity() const;
};

template<class FValue, class FScalar>
FValue FindMinimum(TPolynomial<FValue, FScalar> const &p, FValue tGuess, TPolynomialMinimizeParams<FValue, FScalar> const &Params);



} // namespace ct


#endif // CX_POLY_FIT_H
