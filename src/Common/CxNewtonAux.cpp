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

// #include <cmath>
//
// #include "CxDefs.h"
// #include "CxNewtonAux.h"
//
//
//
// namespace ct_opt {
//
// template <class T>
// inline T sqr(T const &x) { return x*x; }
//
// // -----------------------------------------------------------------------------
// //
// // UPCOMING: functionality for computing generic trust region newton steps
// //
// // -----------------------------------------------------------------------------
//
// static FScalar NtRmsd(FScalar const *p, size_t n, size_t nDegreesOfFreedom) {
//    if (n == 0)
//       return FScalar(0);
//    FScalar
//       f = 0;
//    for (size_t i = 0; i != n; ++ i)
//       f += p[i]*p[i];
//    return std::sqrt(f/FScalar(nDegreesOfFreedom));
// }
//
//
// struct FStepContext {
//    FStepContext(FEwType
//       EwType_, FScalar const *pEw_, FScalar const *gi_, size_t n_, FScalar fTargetStep_, size_t nDegreesOfFreedom_, FScalar fThrSig_)
//       : EwType(EwType_), pEw(pEw_), gi(gi_), n(n_), fTargetStep(fTargetStep_), nDegreesOfFreedom(nDegreesOfFreedom_), fThrSig(fThrSig_)
//    {}
//
//    void MakeStep(FScalar *gi_InvEw);
// protected:
//    FEwType
//       EwType;
//    FScalar const
//       *pEw, *gi;
//    size_t
//       n;
//    FScalar
//       fTargetStep;
//    size_t
//       nDegreesOfFreedom;
//    FScalar
//       fThrSig;
//
//    // returns step rmsd
//    FScalar MakeStepForDamping(FScalar *gi_InvEw, FScalar Damping);
//    void BisectR(FScalar *gi_InvEw, FScalar fLeftPos, FScalar fLeftVal, FScalar fRightPos, FScalar fRightVal);
// };
//
//
// void FStepContext::MakeStep(FScalar *gi_InvEw)
// {
//    FScalar
//       fRightPos = 0.,
//       fRightVal = MakeStepForDamping(gi_InvEw, fRightPos) - fTargetStep;
//    if (fRightVal < 0.)
//       // step with zero damping is already small enough.
//       return;
//
//    FScalar
//       MaxEw = pEw[n-1];
//    if (MaxEw < pEw[0])
//       MaxEw = pEw[0]; // in case of singular values as input (sorted descending) instead of eigenvalues (sorted ascending)
//    FScalar
//       fLeftPos = 1e4*MaxEw*FScalar(n),
//       fLeftVal = MakeStepForDamping(gi_InvEw, fLeftPos) - fTargetStep;
//    assert(fLeftVal < 0.);
//
//    return BisectR(gi_InvEw, fLeftPos, fLeftVal, fRightPos, fRightVal);
// }
//
//
// void FStepContext::BisectR(FScalar *gi_InvEw, FScalar fLeftPos, FScalar fLeftVal, FScalar fRightPos, FScalar fRightVal)
// {
//    for (size_t iIt = 0; iIt < 0x10000; ++ iIt) {
//       FScalar
//          fCenPos = 0.5 * (fLeftPos + fRightPos),
//          fCenVal = MakeStepForDamping(gi_InvEw, fCenPos) - fTargetStep;
//       if (std::abs(fCenVal - fTargetStep) < 1e-8 || std::abs(fCenVal) < 1e-10)
//          return;
//       if (fCenVal < 0.) {
//          fLeftPos = fCenPos;
//          fLeftVal = fCenVal;
//       } else {
//          fRightPos = fCenPos;
//          fRightVal = fCenVal;
//       }
//    }
//    // FIXME: ...this thing has no error reporting functionality...
// }
//
//
// FScalar FStepContext::MakeStepForDamping(FScalar *gi_InvEw, FScalar Damping)
// {
//    for (size_t i = 0; i < n; ++ i) {
//       if (std::abs(pEw[i]) > fThrSig) {
//          if (EwType == EWTYPE_Eigh) {
//             // step for Levenberg damping in eigenvector basis
//             gi_InvEw[i] = -gi[i]/(pEw[i] + Damping);
//          } else if (EwType == EWTYPE_Svd) {
//             // step for Tikhonov damping in singular basis
//             gi_InvEw[i] = -gi[i]*pEw[i]/(sqr(pEw[i]) + sqr(Damping));
//          }
//       } else {
//          gi_InvEw[i] = 0.;
//       }
//    }
//    return NtRmsd(gi_InvEw, n, nDegreesOfFreedom);
// }
//
//
//
//
// void ComputeTrustRegionStep(FScalar *pStep, FScalar const *pEw, FEwType EwType, FScalar const *pGrad, size_t n, FScalar TargetStep, size_t nDegreesOfFreedom, FScalar ThrSigma)
// {
//    FStepContext
//       StepContext(EwType, pEw, pGrad, n, TargetStep, nDegreesOfFreedom, ThrSigma);
//    StepContext.MakeStep(pStep);
// }
//
//
//
// } // namespace ct_opt
