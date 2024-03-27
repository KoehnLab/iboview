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

#ifndef CX_NEWTONAUX_H
#define CX_NEWTONAUX_H

// #include "EfMath.h"
#include <stddef.h> // for size_t

#include "CxDefs.h"
#include <cmath>
#include <limits> // for std::numeric_limits<T>::epsilon


namespace ct_opt {
   enum FEwType {
      // the problem in H^{-1} r is symmetric (and typically positive semi-definite),
      // and the input data in pEw are the eigenvalues (ε_i) of A.
      // In this case, we aim to compute a step length for
      //
      // δx = (H + λ I)^{-1} r
      //
      // under the constraint ‖δx‖₂ <= trust-radius in H's eigenvector basis.
      // The components are
      //
      //  [δx]_i = -1/(ε_i + λ) [r]_i
      EWTYPE_Eigh,
      // the problem in δx = J^{-1} r not symmetric. In this case, the input data
      // in pEw are singular values (σ_i) of J, and we aim to compute
      // a solution of
      //
      // (J^T r) = (J^T J + λ I) J^T δx
      //
      // (or a suitable transpose, depending on the respective dimensionalities),
      // which corresponds to solving
      //
      //   ‖J x - (-r)‖₂ --> min
      //
      // under the constraint ‖δx‖₂ <= trust-radius in J's singular basis.
      // The components are
      //
      //  [δx]_i = -σ_i/(σ_i^2 + λ^2) [r]_i
      //
      // See comments on Tikhonov regularization in info.txt.
      EWTYPE_Svd
   };

   // Compute a trust-region Newton step for quantities given in terms of the eigenvector basis
   // of the Hessian H.
   // Concretely, this computes either
   //
   //    pStep[i] = -pGrad[i]/(pEw[i] + lambda)
   //
   // (for EwType == EWTYPE_Eigh) or
   //
   //    pStep[i] = -pGrad[i]*pEw[i]/(sqr(pEw[i]) + sqr(lambda))
   //
   // (for EwType == EWTYPE_Svd) for all i for which
   // abs(pEw[i]) >= ThrSigma such that the l^2 norm of
   // of the resulting step, divided by sqrt(nDegreesOfFreedom), is either
   // smaller than TargetStep (in which case lambda = 0) or numerically equal
   // to TargetStep (in which case the coresponding damping parameter lambda
   // is computed first)
   //
   // Parameters:
   // - pStep: output step vector
   // - pEw: eigenvalues of the Hessian (for EwType==EWTYPE_Eigh) or
   //   singular values of the Jacobian (for EwType == EWTYPE_Svd)
   // - pGrad: vector of input residuals (typically the gradient of the target function), expressed in the eigenvector basis or right-singular basis.
   // - n: number of vector components; length of pStep, pEw, and pGrad arrays.
   // - ThrSigma: numerical threshold below which eigenvalues of the Hessian are ignored.
   template<class FScalar>
   void ComputeTrustRegionStep(FScalar *pStep, FScalar const *pEw, FEwType EwType, FScalar const *pGrad, size_t n, FScalar TargetStep, size_t nDegreesOfFreedom, FScalar ThrSigma);


   // hm... update formulas (BFGS, SR1, etc.) should probably also go here, shouldn't they?
   // ...but these need matrix and vector types. Maybe just keep it as is instead?


} // namespace ef


namespace ct_opt {

template <class T>
inline T sqr(T const &x) { return x*x; }

// -----------------------------------------------------------------------------
//
// UPCOMING: functionality for computing generic trust region newton steps
//
// -----------------------------------------------------------------------------

template<class FScalar>
static FScalar NtRmsd(FScalar const *p, size_t n, size_t nDegreesOfFreedom) {
   using std::sqrt;
   if (n == 0)
      return FScalar(0);
   FScalar
      f = 0;
   for (size_t i = 0; i != n; ++ i)
      f += p[i]*p[i];
   return sqrt(f/FScalar(nDegreesOfFreedom));
}


template<class FScalar>
struct TStepContext {
   TStepContext(FEwType
      EwType_, FScalar const *pEw_, FScalar const *gi_, size_t n_, FScalar fTargetStep_, size_t nDegreesOfFreedom_, FScalar fThrSig_)
      : EwType(EwType_), pEw(pEw_), gi(gi_), n(n_), fTargetStep(fTargetStep_), nDegreesOfFreedom(nDegreesOfFreedom_), fThrSig(fThrSig_)
   {}

   void MakeStep(FScalar *gi_InvEw);
protected:
   FEwType
      EwType;
   FScalar const
      *pEw, *gi;
   size_t
      n;
   FScalar
      fTargetStep;
   size_t
      nDegreesOfFreedom;
   FScalar
      fThrSig;

   // returns step rmsd
   FScalar MakeStepForDamping(FScalar *gi_InvEw, FScalar Damping);
   void BisectR(FScalar *gi_InvEw, FScalar fLeftPos, FScalar fLeftVal, FScalar fRightPos, FScalar fRightVal);
};


template<class FScalar>
void TStepContext<FScalar>::MakeStep(FScalar *gi_InvEw)
{
   FScalar
      fRightPos = 0,
      fRightVal = MakeStepForDamping(gi_InvEw, fRightPos) - fTargetStep;
   if (fRightVal < FScalar(0))
      // step with zero damping is already small enough.
      return;

   FScalar
      MaxEw = pEw[n-1];
   if (MaxEw < pEw[0])
      MaxEw = pEw[0]; // in case of singular values as input (sorted descending) instead of eigenvalues (sorted ascending)
   FScalar
      fLeftPos = FScalar(1e4)*MaxEw*FScalar(n),
      fLeftVal = MakeStepForDamping(gi_InvEw, fLeftPos) - fTargetStep;
   assert(fLeftVal < FScalar(0));

   return BisectR(gi_InvEw, fLeftPos, fLeftVal, fRightPos, fRightVal);
}


template<class FScalar>
void TStepContext<FScalar>::BisectR(FScalar *gi_InvEw, FScalar fLeftPos, FScalar fLeftVal, FScalar fRightPos, FScalar fRightVal)
{
   using std::abs;
   using std::sqrt;
   FScalar
      ThrX = sqrt(std::numeric_limits<FScalar>::epsilon()), // ~1e-8 ish for doubles
      ThrF = ThrX/100; // ~1e-10ish for doubles
   for (size_t iIt = 0; iIt < 0x10000; ++ iIt) {
      FScalar
         fCenPos = FScalar(0.5) * (fLeftPos + fRightPos),
         fCenVal = MakeStepForDamping(gi_InvEw, fCenPos) - fTargetStep;
//       if (abs(fCenVal - fTargetStep) < FScalar(1e-8) || abs(fCenVal) < FScalar(1e-10))
      if (abs(fCenVal - fTargetStep) < ThrX || abs(fCenVal) < ThrF)
         return;
      if (fCenVal < FScalar(0)) {
         fLeftPos = fCenPos;
         fLeftVal = fCenVal;
      } else {
         fRightPos = fCenPos;
         fRightVal = fCenVal;
      }
   }
   // FIXME: ...this thing has no error reporting functionality...
}


template<class FScalar>
FScalar TStepContext<FScalar>::MakeStepForDamping(FScalar *gi_InvEw, FScalar Damping)
{
   using std::abs;
   for (size_t i = 0; i < n; ++ i) {
      if (abs(pEw[i]) > fThrSig) {
         if (EwType == EWTYPE_Eigh) {
            // step for Levenberg damping in eigenvector basis
            gi_InvEw[i] = -gi[i]/(pEw[i] + Damping);
         } else if (EwType == EWTYPE_Svd) {
            // step for Tikhonov damping in singular basis
            gi_InvEw[i] = -gi[i]*pEw[i]/(sqr(pEw[i]) + sqr(Damping));
         }
      } else {
         gi_InvEw[i] = FScalar(0);
      }
   }
   return NtRmsd(gi_InvEw, n, nDegreesOfFreedom);
}




template<class FScalar>
void ComputeTrustRegionStep(FScalar *pStep, FScalar const *pEw, FEwType EwType, FScalar const *pGrad, size_t n, FScalar TargetStep, size_t nDegreesOfFreedom, FScalar ThrSigma)
{
   TStepContext<FScalar>
      StepContext(EwType, pEw, pGrad, n, TargetStep, nDegreesOfFreedom, ThrSigma);
   StepContext.MakeStep(pStep);
}



} // namespace ct_opt



#endif // CX_NEWTONAUX_H
