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

#ifndef CT_ORB_LOCALIZE_H
#define CT_ORB_LOCALIZE_H

#include "CxIo.h"
#include "CtBasisSet.h"
#include "CtMatrix.h"
#include "CtIao.h"



namespace ct {

// A scalar driving function for generalized PM-like localization methods.
// Typically represents a generic 1-argument scalar function f(n) defined on the
// real interval n\in[0,1], which is convex and at least twice continuously
// differentiable.
//
// Note: in the meantime, some scalar functions with domains and ranges larger
// than [0,1] are also used...
struct FScalarFn01 : public FIntrusivePtrDest {
   // Evaluate the represented function f(n).
   // Arguments:
   // - n: a value in real interval [0,1], to use a function argument in f(n).
   // - d1f, optional: if provided, *d1f will be set to first derivative f'(n).
   // - d2f, optional: if provided, *d1f will be set to second derivative f''(n).
   // Returns:
   //   value of f(n).
   virtual double Eval(double n, double *d1f, double *d2f) const = 0;
   // returns a description of the function h(n); e.g., the implemented formula.
   virtual std::string Desc() const = 0;
   virtual ~FScalarFn01();
   // allows the option to transform the final accumulated localization functional
   // sum into some sort of normalized form---this is only for display purposes atm.
   // It defines the value for L shown in the iteration progress.
   virtual double NormalizedL(double L) const { return L; }
};
typedef TIntrusivePtr<FScalarFn01>
   FScalarFn01Ptr;
typedef TIntrusivePtr<FScalarFn01 const>
   FScalarFn01Cptr;
typedef FScalarFn01 FLocFn1;

FScalarFn01Cptr MakeLocFn1_LogBeta(double ShiftOfN);
FScalarFn01Cptr MakeLocFn1_nPow(double Exp, double Factor=1., double Shift=0.);



struct FLocalizeOptions
{
//    uint
//       nLocExp; //  2 or 4, exponent p in localization functional: L = sum[iA] n(i,A)^p
   // define the orbital functional to maximize with respect to
   // unitary transformations of the occupied orbitals.
   // Note: *all* functionals are maximized; if you have a L[{φ_i}] which
   // *should* be minimized, just flip the sign and maximize -L[{φ_i}] instead.
   enum FLocalizeFn {
      // do not apply any changes to the input orbitals.
      LOCFN_None,
      // L = \sum_{i,A} n(i,A)^p with p=2  (--> maximize)
      //
      // Note: there is no separate p=3 entry, because the n(i,A)^2 and n(i,A)^3
      // functionals have identical maxima. The two functionals are equivalent.
      LOCFN_nA_pow_2,
      // L = \sum_{i,A} n(i,A)^p with p=4  (--> maximize)
      LOCFN_nA_pow_4,
      // L = \sum_{i,A} n(i,A)^p with p=3/2  (--> maximize)
      LOCFN_nA_pow_15,
      // L = \sum_{i,A} h(n(i,A)) with h(n) = n log2(n)
      // (=(-1)*entroy --> maximize this to minimize information entropy
      // associated with the probability distribution associated with
      // electron presence on fragments.)
      LOCFN_nA_Entropy,
//       // L = \sum_{i,A} h(n(i,A)) with h(n) = log_2 Gamma(n+2)
//       LOCFN_nA_Log2Gamma_nPlus2,
//       // L = \sum_{i,A} h(n(i,A)) with h(n) = log_2 Gamma(n+1)
//       LOCFN_nA_Log2Gamma_nPlus1,
      // L = \sum_{i,A} h(n(i,A)) with h(n) supplied by pScalarFn01.
      LOCFN_nA_gen_h,
      // L = \sum_{i} <φ_i| (f - <f>)^2 |φ_i>  (--> maximize)
      // (note that this is zero for canonical SCF orbitals)
      LOCFN_FockVariance
   };
   FLocalizeFn
      LocFn;
   uint
      nMaxIt;
   double
      ThrLoc;
//    std::ostream
//       *pxout; // if non-0, write output here.
   double
      // if != 0, apply an random unitary transformation with this many degrees
      // std deviation as perturbation to the input vectors of the localization.
      fRandomAngleDeg;
   FLog
      *pLog; // if non-0, write output here.
   uint
      Verbosity;
   std::string
      TaskDesc;
   FScalarFn01Cptr
      // ignored unless LocFn is set to LOCFN_nA_gen_h; in that case,
      // *pScalarFn01 provides the h(n) function defining a generic nA-based
      // localization functional L = \sum_{i,A} h(n(i,A)) with h(n)
      // Note:
      // - L is *maximized*. So in order to provide drive towards
      //   localization, h(n) should be a convex function.
      pScalarFn01;

//    FLocalizeOptions()
//       : nLocExp(4), nMaxIt(2048), ThrLoc(1e-8), pxout(0), Verbosity(1)
//    {}
   FLocalizeOptions();

   std::string MakeLocalizeFnDesc() const;
};



// execute a 2x2 rotation of the vectors pA and pB (both with length nSize)
// of angle phi. In-place.
void Rot2x2(double *RESTRICT pA, double *RESTRICT pB, size_t nSize, double phi);

void LocalizeVectors(FMatrixView CVec, size_t const *pGroupOffsets, size_t nGroups, FLocalizeOptions const &Opt, FMemoryStack &Mem);



// Construct positive Cayley transform: CayU(A) = 1 + A + A^2/2 + A^3/4 + A^4/8 + ...
// (that's not the self-inverse one!)
// Input matrix must be anti-symmetric. Output matrix is allocated on Mem.
FMatrixView MakeCayleyU(FMatrixView A, FMemoryStack &Mem);

// construct a random N x N anti-symmetric matrix, in which angles are approximately
// distributed according to a normal distribution with standard deviation fAngleStdDev
// (given in radians)
FMatrixView MakeRandomAntiSymmetricMatrix(size_t N, double fAngleStdDev, FMemoryStack &Mem);

// apply a random rotation (unitary transformation) to the vectors in the rows
// of CVec, with the perturbing angles being approximately normal-distributed
// around 0 with standard deviation fAngleStdDev
void RotateVectorsRandomly(FMatrixView CVec, double fAngleStdDev, FMemoryStack &Mem);



}



#endif // CT_ORB_LOCALIZE_H
