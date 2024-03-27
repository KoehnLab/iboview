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

#include <cmath>
#include <stdexcept>
#include "CxVec3.h"
#include "IrAmrr.h"
#include "Ir.h"
#include "IrEcp.h"
#include "IrFactorials.h"

#define IR_ECP_DEBUG_PRINT
#ifdef IR_ECP_DEBUG_PRINT
   #include <iostream>
   #include "format.h"
#endif

using namespace ct;
using std::size_t;
typedef ct::TVector3<double>
   FVec3;

#include "IrComplexRlm.h"

using ir_rlm::complex_double;
using ir_rlm::iSlmX;
using ir_rlm::nSlmX;
using ir_rlm::nSlmY;
using ir_rlm::EvalRlmX;
using ir_rlm::ConvertRlmToSlm1;

namespace ir {

template<class T>
static T sqr(T x) { return x*x; }


// Output: pOut[(iSt * (iM + (MaxM+1) * iCo)]
// Input: pCo[nExp * iCo + iExp]
static void EvalF(double *pOut, size_t iSt, double A, double r, size_t MaxM, FRawShell const *pSh, FMemoryStack &Mem)
{
   for (size_t iCo = 0; iCo < pSh->nCo; ++ iCo)
      for (size_t iM = 0; iM <= MaxM; ++ iM)
         pOut[iSt*(iM + (MaxM+1) * iCo)] = 0.;

   size_t
      nExp = pSh->nExp;
   TMemoryLock<double>
      K(MaxM+1, &Mem);
   for (size_t iExp = 0; iExp < nExp; ++ iExp) {
      double
         Zeta = pSh->pExp[iExp];
      double
         fExpFactor = std::exp(-Zeta*sqr(r-A));
      EvalModifiedSphericalBesselEnveloped(K, 2*Zeta*r, A, MaxM, fExpFactor);
      // ^- note: unlike the article version, this one has A^{-m} factors
      // absorbed (in such a way that the result is well defined when A=0) and
      // therefore can/must be used in conjunction with R_{lm}(A) solid
      // harmonics, NOT spherical harmonics R_{lm}(A/|A|). (note: that is good,
      // because this one is well-defined for A=0)
      for (size_t iCo = 0; iCo < pSh->nCo; ++ iCo)
         for (size_t iM = 0; iM <= MaxM; ++ iM)
            pOut[iSt*(iM + (MaxM+1) * iCo)] += K[iM] * pSh->fCo(iExp,iCo);
   }
}


// compute (-1)^m.
static double PowMinus1(int m) {
   if (m < 0)
      m = -m;
   return ((size_t(m) % 2) == 0)? 1. : -1.;
}


// Convert (2*la+1)*(2*lb+1) matrix of complex two-index integrals over Rlm(l,m) at pIn[iMa*iStInA + iMb*iStInB]
// to real two-index integrals over Slc(l,c) at pOut[iCa*iStOutA+iCb*iStOutB]. Function assumes integrals over
// output real functions to be real. Output real functions are in Slc order. Function INCREMENTS OUTPUT.
void ConvertIntMatrixRlmToSlm(double *pOut, size_t soA, size_t soB, complex_double const *pIn, size_t siA, size_t siB, int la, int lb, FMemoryStack &Mem)
{
   assert(la >= 0 && lb >= 0 && la <= 6 && lb <= 6);
   if (lb > la) {
      std::swap(la,lb);
      std::swap(soA,soB);
      std::swap(siA,siB);
   }
   TMemoryLock<complex_double>
      pHalf(size_t(2*la+1)*size_t(2*lb+1), &Mem),
      pFull(size_t(2*la+1)*size_t(2*lb+1), &Mem);
   // ^-- hmpf. g++ 9.3.1 emits scary looking warnings for the class instantiation
   // of TMemoryLock having a memset(p, 0, ..) in its Clear() routine, and appyling
   // this to a std::complex class. Apart from the fact that memset on this is very most likely
   // totally fine (for class consisting of two IEEE floats as only data members),
   // g++'s code analysis doesn't see that .Clear() is actually NEVER CALLED on this.
   // It could see it... all the stuff is inline and the control constants are compile
   // time consts.
   for (size_t imb = 0; imb < size_t(2*lb+1); ++ imb) {
      ir_rlm::ConvertRlmToSlm1(&pHalf[size_t(2*la+1)*imb], 1, &pIn[imb*siB], siA, la);
   }
   for (size_t ima = 0; ima < size_t(2*la+1); ++ ima) {
      ir_rlm::ConvertRlmToSlm1(&pFull[ima], size_t(2*la+1), &pHalf[ima], size_t(2*la+1), lb);
   }
   for (size_t ima = 0; ima < size_t(2*la+1); ++ ima) {
      for (size_t imb = 0; imb < size_t(2*lb+1); ++ imb) {
         complex_double f = pFull[ima + size_t(2*la+1)*imb];
         assert_rt(std::abs(f.imag()) < 1e-8);
         pOut[ima*soA + imb*soB] = f.real();
      }
   }
}



// Eval pOut[iRadialPt, iSlmX(lEcp,mEcp), iComplexFnA] = I^a_{lm}(r)
// = \int R^l_m(\vec r/r) \phi_a(\vec r)\,\d\Omega
//
// NOTE:
// - Not that this makes the angular parts for *ALL* ECP output components
//   (iSlmX(lEcp,mEcp)) at once.
//
// - So you still need to assemble the correct corresponding subset (multiplying
//   the appropriate U_l integrals into the appropriate (2*l+1) half-angle
//   integral components; one-sided) when building the target quantities.
static void EvalHalfAngularIntegrals(double *pOut, size_t MaxEcpL_, double const *pRadialPt, size_t nRadialPt, FRawShell const *pA, FVec3 vA, FMemoryStack &Mem)
{
   ct::TMemoryLock<double>
      pFreeMe(0, &Mem);
   double
      fNormA = vA.Length();

   int
      MaxEcpL = int(MaxEcpL_),
      la = int(pA->l),
      MaxLambda = MaxEcpL + la;
   ct::TMemoryLock<complex_double>
      pRlmXA(nSlmX(MaxLambda), &Mem);
   EvalRlmX(pRlmXA, vA[0], vA[1], vA[2], MaxLambda);

   // compute the radial factors F[iRadius,Lambda,iCoA]
   // (note: in principle we could cache them for all pShA/vC (ECP center) pairs.)
   ct::TMemoryLock<double>
      pFa(nRadialPt * (MaxLambda+1) * pA->nCo, &Mem);
   for (size_t iRadialPt = 0; iRadialPt != nRadialPt; ++ iRadialPt)
      EvalF(&pFa[iRadialPt], nRadialPt, fNormA, pRadialPt[iRadialPt], MaxLambda, pA, Mem);

   size_t
      // = (2*la+1) * nCoA, but with the m referring to complex Rlm-harmonics,
      // not Slc as usual.
      nComplexFnA = pA->nFn();
   ct::TMemoryLock<complex_double>
      pComplexInt(nRadialPt * nSlmX(MaxEcpL) * nComplexFnA, &Mem, ALLOC_Clear);

   // and assemble the radial and angular factors...
   for (int ma = -int(la); ma <= int(la); ++ ma) { // m component of output function
      for (int ltick = 0; ltick <= la; ++ ltick) {
         for (int mtick = -int(ltick); mtick <= int(ltick); ++ mtick) {
            if (ltick > la || std::abs(ma-mtick) > la-ltick)
               continue;
            if ((ltick + mtick < 0) || (ltick + mtick) > (la + ma))
               // ^- condition comes from binomial factor in Rlm addition theorem
               // for R^la_ma(r - A) into R^l'_m'(r) and R^{la-l'}_{ma-m'}(A).
               continue;
            complex_double
               FactorA = GetBinomialCoeffF(la + ma, ltick + mtick);

            FactorA *= PowMinus1(la-ltick); // (-1)^(la-l')
            assert(la - ltick <= MaxLambda);
            assert(int(la-ltick) >= 0 && std::abs(ma-mtick) <= la-ltick);
            FactorA *= pRlmXA[iSlmX(int(la-ltick), ma-mtick)];

//             FactorA *= pow(fNormA, la-ltick);
            // ^- cgk 2020-05-07: changed the bessel thing to IrEcpBesselFn.cpp version,
            // which allows absorbing the A^{-lambda} factors and can therefore be
            // used with solid harmonics instead of spherical harmonics.

            for (int EcpL = 0; EcpL <= MaxEcpL; ++ EcpL) {
               int MinLambda = std::max(ltick,EcpL) - std::min(ltick,EcpL);
               for (int lambda = MinLambda; lambda <= ltick + EcpL; ++ lambda) {
                  if ((size_t(lambda + ltick + EcpL) & 1) != 0)
                     // 3R integral zero otherwise.
                     continue;
                  for (int EcpM = -int(EcpL); EcpM <= int(EcpL); ++ EcpM) {
                     int mu = -(mtick + EcpM);
                     if (std::abs(mu) > lambda)
                        continue;

                     complex_double
                        FactorLambda = double(2*lambda+1)*PowMinus1(mu);
                     assert(lambda <= MaxLambda);
                     FactorLambda *= pRlmXA[iSlmX(int(lambda), -mu)];
                     // and the three-spherical-harmonic integral...
                     double
                        Factor3R = ir_rlm::GetSphereIntegral3Rlm(int(ltick), mtick, int(EcpL), EcpM, int(lambda), mu);
                     FactorLambda *= Factor3R;

                     if (FactorLambda == 0.)
                        continue; // there are still a bunch of zero angular factors, even with the fixed mu.

                     for (size_t iRadialPt = 0; iRadialPt != nRadialPt; ++ iRadialPt) {
                        double
                           FactorR = std::pow(pRadialPt[iRadialPt], double(ltick));
                        // ^- this one and the lambda index are a bit annoying. Without
                        //    them we could precalculate the angular stuff and put a
                        //    radial loop on the outside, each time making only a
                        //    single complex matrix and transforming to real right away.
                        //    But as it is, I did not see a good way of doing so.
                        for (size_t iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
                           complex_double
                              FinalIntegral = FactorA * FactorLambda * FactorR;
                           FinalIntegral *= pFa[iRadialPt + nRadialPt*(lambda + (MaxLambda+1) * iCoA)];
                           size_t iFnA = (la+ma) + (2*la+1)*iCoA;
                           pComplexInt[iRadialPt + nRadialPt * (iSlmX(int(EcpL), EcpM) + nSlmX(MaxEcpL) * iFnA)] += FinalIntegral;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // convert polynomial prefactors of ECP projectors and basis functions
   //
   //     phi_a(\vec r) = R^l_m(\vec r - \vec A) \left(\sum_{i} d_{ia} exp(-zeta_{ia} (\vec r - \vec a)^2) \right)
   //
   // from unnormalized complex solid harmonics R^l_m(r - A), which we currently
   // have, into real solid harmonics S^l_c(r - A) in our target component
   // order. This is done by linearly combining, renormalizing, and re-ordering
   // integrals over suitable linear combinations of R(l,+m) and R(l,-m).
   for (size_t iRadialPt = 0; iRadialPt != nRadialPt; ++ iRadialPt) {
      for (int lEcp = 0; lEcp <= MaxEcpL; ++ lEcp) {
         size_t stFn = nRadialPt*nSlmX(MaxEcpL);
         for (size_t iCoA = 0; iCoA != pA->nCo; ++ iCoA) {
            size_t
               iFnA = size_t(2*la+1)*iCoA;
            size_t
               iOff = iRadialPt + nRadialPt * size_t(nSlmX(lEcp-1)) + stFn*iFnA;
             ConvertIntMatrixRlmToSlm(&pOut[iOff], nRadialPt, stFn, &pComplexInt[iOff], nRadialPt, stFn, lEcp, la, Mem);
         }
      }
   }
   // NOTE: m >= 0 and .real/.imag vs R_{l,+m} +/- R_{l,-m}
   //
   // - We *probably* could also just make the complex m >= 0 output components,
   //   and then get the real Slc ones by just taking .real()/.imag() components
   //   of them, instead of assembling the linear combinations as
   //   ConvertRlmToSlm1 currently does.
   //
   //   UPDATE:
   //   - The m >= 0 thing can only be done for *EITHER* the iSlmX(ECP)
   //     *OR* the iSlmY(A) components of the half-angle integral, not for both.
   //
   //   - Reason is that doing so once would discard the output real/imag
   //     components, which is not acceptable if they are still needed for the
   //     second Rlm -> Slc transform.
   //
   //   - Might still be a good idea to implement it at the end. Would allow
   //     reduction of numerical work in half angle integral assembly by factor
   //     ~2. (by simply not computing the corresponding target m < 0
   //     components)
}



static double EvalScalarEcpForGivenL(int l, FAtomEcp const *pEcp, double *pPowR, size_t MaxPowR, double Prefactor, ct::FMemoryStack &Mem)
{
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   assert(MaxPowR >= 2);
   double
      Ul = 0.,
      r = pPowR[1],
      r2 = pPowR[2];
   FAtomEcp::FEcpShellFnList::const_iterator
      itSh;
   // go through ECP shells and check which ones have the requested angular momentum.
   for (itSh = pEcp->ShellFns.begin(); itSh != pEcp->ShellFns.end(); ++ itSh) {
      if (itSh->l != l)
         continue;
      for (size_t iExp = 0; iExp < itSh->nExp; ++ iExp) {
         double
            fExp = itSh->pExp[iExp],
            fCo = itSh->pCo[iExp],
            fPow = itSh->pPow[iExp];
         double
            PowR;
         size_t
            iPow = size_t(fPow);
         if (double(iPow) == fPow && iPow <= MaxPowR)
            // stored integer power of r?
            PowR = pPowR[iPow];
         else
            PowR = std::pow(r, fPow);
         // ^- note: in practice we only ever end up with fPow == 0 (and
         // therefore PowR=1) here, at least for Stuttgart-Koeln-Bonn ECPs...
         Ul += Prefactor * fCo * PowR * std::exp(-fExp * r2);
      }
   }
   return Ul;
}


static void AssembleEcpAndHalfIntegrals(double *pIntAB, double *pHalfIntA,
   FRawShell const *pA, double *pHalfIntB, FRawShell const *pB, FAtomEcp const *pEcp,
   double *pRadialPt, double *pRadialWt, size_t nRadialPt, double Prefactor, ct::FMemoryStack &Mem)
{
   size_t
      MaxEcpL = size_t(pEcp->MaxL());
   assert(MaxEcpL <= 5); // if larger there is probably a integer overflow... we should not get here for pure local ECPs.

   size_t
      nFnA = pA->nFn(), // number of non-redundant COMPLEX components
      nFnB = pB->nFn();

   for (size_t i = 0; i != nFnA*nFnB; ++ i)
      pIntAB[i] = 0;

   size_t
      nSlmX_MaxEcpL = nSlmX(int(MaxEcpL));

   // compute the contracted scalar part of the ECP for the all radii,
   // for all non-local EcpL it has (i.e., for EcpL != -1).
   TMemoryLock<double>
      pUl(nRadialPt * (MaxEcpL+1), &Mem);
   for (size_t iRadialPt = 0; iRadialPt < nRadialPt; ++ iRadialPt) {
      double
         r = pRadialPt[iRadialPt];
      double
         pPowR[3] = {1., r, r*r};

      for (int EcpL = 0; EcpL <= int(MaxEcpL); ++ EcpL) {
         double f = Prefactor * double(2*EcpL + 1) * pRadialWt[iRadialPt];
         pUl[iRadialPt + nRadialPt*size_t(EcpL)] = EvalScalarEcpForGivenL(EcpL, pEcp, &pPowR[0], 2, f, Mem);
      }
   }
   for (size_t EcpL = 0; EcpL <= MaxEcpL; ++ EcpL) {
      double const
         // get offsets of half-integral data for current EcpL subset.
         *pHalfIntA_EcpL = &pHalfIntA[nRadialPt * nSlmX(int(EcpL)-1)],
         *pHalfIntB_EcpL = &pHalfIntB[nRadialPt * nSlmX(int(EcpL)-1)];
      size_t
         stFnHlf = nRadialPt * nSlmX_MaxEcpL; // stride of iFn in pHalfIntA / pHalfIntB
      size_t
         nEcpM = 2*EcpL + 1,
         // that's the amount of values we need to contract over:
         // all radial points and the (2*EcpL+1) EcpM components for
         // the current EcpL.
         nRadialPtEcpM = nRadialPt * nEcpM;
      // absorb Ul factors into one of the half-integrals, so that we can get
      // the rest of the operation into the form of a matrix multiplication
      TMemoryLock<double>
         pHalfIntB_Ul(nRadialPt * nEcpM * nFnB, &Mem);
      for (size_t iEcpM = 0; iEcpM < nEcpM; ++ iEcpM)
         for (size_t iFnB = 0; iFnB < nFnB; ++ iFnB)
            for (size_t iRadialPt = 0; iRadialPt < nRadialPt; ++ iRadialPt)
               pHalfIntB_Ul[nRadialPt*iEcpM + iRadialPt + nRadialPtEcpM*iFnB] = pUl[iRadialPt + nRadialPt*size_t(EcpL)] * pHalfIntB_EcpL[nRadialPt*iEcpM + iRadialPt + stFnHlf*iFnB];
      // compute output matrix element. Note:
      // - This is a matrix multiplication, but I do not wish to add a BLAS/libxsmm
      //   dependency to IR at this moment---not unless required. Should be
      //   checked later on, once things work!
      // - Also, it's in TN form, which on modern machines is not very efficient.
      //   if turned into a real MxM, the pHalfIntB_Ul above should be transposed.
      // - Atm it seems to be fast enough to deal with a small number of ECP atoms
      //   in large TM complexes, even with jfit/jxfit DFT. So for I'll just leave
      //   it as is until a more pressing need to deal with it comes around.
      for (size_t iFnA = 0; iFnA < nFnA; ++ iFnA)
         for (size_t iFnB = 0; iFnB < nFnB; ++ iFnB)
            for (size_t i = 0; i < nRadialPtEcpM; ++ i)
               pIntAB[iFnA + nFnA * iFnB] += pHalfIntA_EcpL[i + stFnHlf*iFnA] * pHalfIntB_Ul[i + nRadialPtEcpM*iFnB];
   }
}




static void ScatterAcc2(double *pOut, size_t StrideA, size_t StrideB, double const *pIntAB, size_t nFnA, size_t nFnB)
{
   for (size_t iFnB = 0; iFnB != nFnB; ++ iFnB)
      for (size_t iFnA = 0; iFnA != nFnA; ++ iFnA)
         pOut[StrideA * iFnA + StrideB * iFnB] += pIntAB[iFnA + nFnA*iFnB];
}



// note: ri and wi allocated on Mem.
static void _MakeFixedRadialGrid(double *&ri, double *&wi, size_t nPt, FMemoryStack &Mem)
{
   Mem.Alloc(ri, nPt);
   Mem.Alloc(wi, nPt);
   double
      den = 1./double(nPt);

   for (size_t i = 0; i < nPt; ++ i) {
      // Mura & Knowles Log3-grid.
      double
         R = 5., // defines center position of the grid.
//          R = 1.0, // we want to integrate integrands multiplied by CORE functions! (U_l/U_loc)
         // goes from 0...1. (exclusive)
         x = ((.5 + double(i))*den);
      double
         r = -R * std::log(1.-x*x*x),
         dr = R * 3*x*x/(1.-x*x*x),  // dr/dx for dx -> dr radial volume element transform (dx = 1/nPt with the simple linear weighting scheme).
         w = den * dr * 4.*M_PI*r*r; // last: volume element (4 pi r^2)
      ri[i] = r;
      wi[i] = w;
   }
}




static void EvalEcpLocalPart(double *pOut, size_t StrideA, size_t StrideB, FRawShell const *pA, FRawShell const *pB, double Prefactor, bool Add, FEcpShellFn const *pEcpLocalSh, double const *pEcpCen, ct::FMemoryStack &Mem)
{
   // NOTE:
   // This function appears to work. It just relays the local part of the ECP to
   // a inline-contracting IR 3c integral with an overlap kernel.
   //
   // (side note: if the persons who came up with the original ECP operator
   //   forms would have stopped for a few seconds to think about what they were
   //   doing, the entire rest of this ECP mess could have been done like this,
   //   too. I see no legitimate reason for the actual ECP operator forms. This
   //   is just insanity! I do not even want to know how many year of highly
   //   qualified researchers were wasted on dealing with this *ENTIRELY*
   //   unnecessary problem...)
   //
   // ...in any case: The real-life ECPs of the Stuttgart/Koeln/Bonn form appear
   // to not actually have a local part of the ECPs. For this reason this
   // function is not actually called. However, I did try it for some cases and
   // it worked fine as long as the zero-powers are used as normally.

   ct::TMemoryLock<double>
      pFreeMe(0, &Mem);

   if (pEcpLocalSh->l != -1)
      throw std::runtime_error("attempted to issue local ECP call for non-local part of actual ECP.");
   size_t
      nExp = pEcpLocalSh->nExp;
   for (size_t i = 0; i != nExp; ++ i) {
      if (pEcpLocalSh->pPow[i] != 0.)
         throw std::runtime_error("analytic evaluation of local part of ecp supports raw Gaussian ECPs only, without non-zero monomial r-power prefactors.");
   }

   // make a Gaussian basis function shell of s-type with the exponents and coefficients
   // given by the ECP on the ECP's center. That's really all it is...
   FRawShell
      EcpRep; // = {0};
   EcpRep.l = 0;
   EcpRep.nExp = nExp;
   EcpRep.pExp = pEcpLocalSh->pExp;
   EcpRep.nCo = 1;
   EcpRep.pCo = pEcpLocalSh->pCo;
   EcpRep.pRange = 0;
   EcpRep.vCen = pEcpCen;

   double
      *pCoeffsC;
   Mem.Alloc(pCoeffsC, nExp);
   for (size_t i = 0; i != nExp; ++ i)
      pCoeffsC[i] = 1.0;

   if (0) {
      size_t
         Strides[2] = {StrideA, StrideB};
      EvalInt2e3c_ContractC(pOut, &Strides[0], pCoeffsC, 1, 1,
         pA, pB, &EcpRep, 1, Prefactor, &g_IrOverlapKernel, Mem);
      // ^- can't write to pOut with output strides directly, because this would overwrite.
      //    And this function is meant to increment, not overwrite.
      //    So we make a copy and add it over explicitly...
      // FIXME: wait... does that REALLY overwrite instead of increment?
   } else {
      size_t
         nFnA = pA->nFn(),
         nFnB = pB->nFn(),
         Strides[2] = {1, nFnA};
      double
         *pResult; // the function doesn't add.
      Mem.ClearAlloc(pResult, nFnA * nFnB);
      EvalInt2e3c_ContractC(pResult, &Strides[0], pCoeffsC, 1, 1,
         pA, pB, &EcpRep, 1, Prefactor, &g_IrOverlapKernel, Mem);

      for (size_t iFnB = 0; iFnB != nFnB; ++ iFnB)
         for (size_t iFnA = 0; iFnA != nFnA; ++ iFnA) {
            if (Add)
               pOut[iFnA * StrideA + iFnB * StrideB] += pResult[Strides[0] * iFnA + Strides[1] * iFnB];
            else
               pOut[iFnA * StrideA + iFnB * StrideB] = pResult[Strides[0] * iFnA + Strides[1] * iFnB];
         }
   }
}


// Note: should probably remove this from free standing functions and make a Meta2i-kernel.
// Such a kernel object could also get a set of "prepare" functions, which get passed both the A and B
// basis sets at the start of the MakeIntMatrix functions. In most kernels this would do nothing, but
// for ECP evaluation (which could be expensive), it might be quite useful to calculate a bunch of data
// related to the A/B shells individually, which otherwise would have to be re-computed redundantly.
void EvalEcp(double *pOut, size_t StrideA, size_t StrideB, FRawShell const *pA, FRawShell const *pB,
   double Prefactor, bool Add, FAtomEcp const *pEcp, double const *pEcpCen, ct::FMemoryStack &Mem)
{
   ct::TMemoryLock<double>
      pFreeMe(0, &Mem);
   FVec3
      vC = FVec3(pEcpCen),
      // we'll use coordinates relative to vC throughout.
      vA = FVec3(pA->vCen[0], pA->vCen[1], pA->vCen[2]) - vC,
      vB = FVec3(pB->vCen[0], pB->vCen[1], pB->vCen[2]) - vC;
   // make a radial grid for numerically evaluating the radial part of the 6D integral.
   double
      *pRadialPt, *pRadialWt;
   size_t
      // note that the radial grid really only has to cover the range of
      // functions U_l and U_loc inside the ECP. It does *not* have to cover all
      // of space! We could even estimate the outer radius using the primitives
      // (IrCore's geminal kernel has the code for that), and put a grid which
      // includes the exact nuclear position and approx and outer integration
      // limit
//       nRadialPt = 5;
//       nRadialPt = 8;
//       nRadialPt = 20;
//       nRadialPt = 40;
//       nRadialPt = 200;
      nRadialPt = 80;
      // ^- 80 with the Log3 grid seems to produce numerically effectively exact
      //    results (15 digits) for the ECPs I tested. I'll just leave it alone
      //    for the time being.

   _MakeFixedRadialGrid(pRadialPt, pRadialWt, nRadialPt, Mem);
   // ^- Note:
   // - These radial grids have the 4\pi r^2 of the radial volume element
   //   already included. For this reason we omit the corresponding r^2 and 4\pi
   //   factors later on.
   // - Concretely, our 3-R integrals have no leading 4\pi factor, and this cancels
   //   the 4\pi included here.

   int
      MaxEcpL = pEcp->MaxL();

   if (!Add) {
      // yes.. that kind of defeats the purpose.
      ZeroBlock2(pOut, StrideA, StrideB, pA->nFn(), pB->nFn());
      Add = true;
   }

   // add local part in closed form as 3c integrals.
   // Note: actual Stuttgart-type ECPs do not appear to have one...
   {
      FAtomEcp::FEcpShellFnList::const_iterator
         itSh;
      // go through ECP shells and check which ones have the requested angular momentum.
      for (itSh = pEcp->ShellFns.begin(); itSh != pEcp->ShellFns.end(); ++ itSh) {
         if (itSh->l == -1)
            EvalEcpLocalPart(pOut, StrideA, StrideB, pA, pB, Prefactor, Add, &*itSh, pEcpCen, Mem);
      }
   }

   if (MaxEcpL != -1) {
      // ECP has a non-local part.
      size_t
         nFnA = (2*pA->l + 1) * pA->nCo,
         nFnB = (2*pB->l + 1) * pB->nCo;

      TMemoryLock<double>
         pIntAB(nFnA * nFnB, &Mem);

      {
         TMemoryLock<double>
            pHalfIntA(nRadialPt * nSlmX(size_t(MaxEcpL)) * nFnA, &Mem),
            pHalfIntB(nRadialPt * nSlmX(size_t(MaxEcpL)) * nFnB, &Mem);

         EvalHalfAngularIntegrals(pHalfIntA, size_t(MaxEcpL), pRadialPt, nRadialPt, pA, vA, Mem);
         EvalHalfAngularIntegrals(pHalfIntB, size_t(MaxEcpL), pRadialPt, nRadialPt, pB, vB, Mem);

         assert_rt(Prefactor == 1.0); // doesn't appear anywere?
         AssembleEcpAndHalfIntegrals(pIntAB, pHalfIntA, pA, pHalfIntB, pB, pEcp, pRadialPt, pRadialWt, nRadialPt, Prefactor, Mem);
         // ^- warning: this CHANGES pHalfIntA and/or pHalfIntB!

         ScatterAcc2(pOut, StrideA, StrideB, pIntAB, nFnA, nFnB);
      }
   }


   // - Note that, theoretically, we *could* reduce the amount of numerical work
   //   for the half-angle integrals...
   //
   //   UPDATE:
   //   - And, technically, also of evaluating the kernel values U_l(r) on the
   //     grid points, and assembling the grids themselevs, although neither of
   //     these aspects appears to be particularly important at first glance
   //
   //   - (This actually may provide an interesting algorithm for constructing
   //     more general integrals; not only ECP ones...)
   //
   //   ...by O(N), where N is the molecule size, by for each ECP center first
   //   precomputing them *for the entire basis set* for A/B shells before
   //   assembling the ECP integrals:
   //
   //   for C in EcpCenters:
   //      // data between different EcpL is not mixed; may be best to
   //      // treat them separately (i.e., make iSlmY(ECP) output half-angle
   //      // integrals, not iSlmX(ECP) output half-angle integrals; or
   //      // at least re-order the components into:
   //      //
   //      //  pHlfX[iRadialPt,mEcp,nFnA,nDerivA,lEcp]
   //      //
   //      // with mEcp next to iRadialPt and lEcp as outer axis.
   //      // (for the MxM/Syrk at end!)
   //      for EcpL:
   //         // make a single shared radial grid pRadialPt/pRadialWt
   //         // which is suitable for accurately representing the functions
   //         // U_l(r) of the ECP center C (that's all it has to! Since
   //         // the U_l(r) will be multiplied as envelopes to whatever else
   //         // comes out on the grid points.)
   //         _MakeFixedRadialGrid(...)
   //
   //         // make all solid angle integrals of A basis functions contracted
   //         // to spherical components of Rlm, on numerical radial grid:
   //         //
   //         // pHlfA(r) = \int R_{lm} \phi_a \d\Omega
   //         for A in BasisA:
   //            for iRadialPt:
   //               pHlfA[iRadialPt,iFnA,iSlmX(ECP)] = ...
   //
   //         // if BasisB != BasisA, assemble the same integrals for BasisB
   //         for B in BasisB:
   //            for iRadialPt:
   //               pHlfB[iRadialPt,iFnB,iSlmX(ECP)] = ...
   //
   //         // compute U_l potential values on the radial grid points...
   //         for iRadialPt:
   //            Ul[iRadialPt] = ...
   //            if BasisA != BasisB:
   //               // absorb Ul into pHlfA
   //               pHlfA[iRadialPt,:,iSlmX(l,:)] *= Ul[iRadialPt]
   //            else:
   //               if ECP.IsPositiveDefinite:
   //                  assert(&pHlfA == &pHlfB)
   //                  pHlfA[iRadialPt,:,iSlmX(l,:)] *= std::sqrt(Ul[iRadialPt])
   //               else:
   //                  assert(&pHlfA != &pHlfB)
   //                  // absorb Ul into COPY(!) of pHlfA
   //                  pHlfA[iRadialPt,:,iSlmX(l,:)] *= Ul[iRadialPt]
   //
   //         // assemble output integral matrix in one step via MxM or Syrk
   //         // by contracting over the (nRadialPt * (2*EcpL+1)) real components
   //         // of the half-angle integrals...
   //         Out[iFnA,iFnB] = \sum_{i,m} pHlfA[(iRadialPt,mEcp),iFnA] * pHlfB[(iRadialPt,mEcp),iFnB]
   //
   //
}


int FAtomEcp::MaxL() const
{
   int lmax = -1;
   FEcpShellFnList::const_iterator
      itSh;
   for (itSh = ShellFns.begin(); itSh != ShellFns.end(); ++ itSh)
      lmax = std::max(int(itSh->l), lmax);
   return lmax;
}


double FAtomEcp::MaxPowR() const
{
   double MaxPowR = 0.;
   FEcpShellFnList::const_iterator
      itSh;
   for (itSh = ShellFns.begin(); itSh != ShellFns.end(); ++ itSh) {
      for (size_t iExp = 0; iExp < itSh->nExp; ++ iExp)
         MaxPowR = std::max(MaxPowR, double(itSh->pPow[iExp]));
   }
   return MaxPowR;
}



} // namespace ir

// coding: utf-8
