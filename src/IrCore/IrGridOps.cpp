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
#include "CxVec3.h"
#include "IrAmrr.h"
#include "Ir.h"
#include "IrInternal.h"

using namespace ct;
typedef ct::TVector3<double>
   FVec3;


typedef unsigned long
   FCoMask;


#ifdef _DEBUG
   #include <iostream>
   #include "format.h"
   using std::cout;
   using std::endl;


   static void PrintArray(char const *pTitle, double *pValues, size_t nValues, int nStep = -1, char const *pFmt = "{:11.5f}")
   {
      if (nStep == -1) {
         // find a reasonable step size such that we get coverage of the entire
         // array but don't print out more than ~10-15 numbers.
         nStep = (12 + nValues)/13;
      }
      cout << fmt::format("[{:>5}] {:<20}", nValues, pTitle);
      for (size_t i = 0; i < nValues; i += nStep)
         cout << fmt::format(pFmt, pValues[i]);
      cout << "\n";
   }

   void _DummyToSupressWarnings() {
      (void)PrintArray;
   }
#endif


namespace ir {

// static size_t const g_nGridDerivComp[3] = {1,4,10};
// static size_t const g_nGridDerivCompLaplace[3] = {2,5,14};
static size_t const
   g_nGridDerivComp1[2][3] = {{1,4,10}, {2,5,14}};

static bool const
   s_GridScreeningEnabled = true;

// return total number of output derivative components for the basis functions
// for the given combintaion of derivative order and 'MakeLaplaceToo'.
// e.g., for DerivOrder == 1 and MakeLaplaceToo == false, these are [1, d/dx, d/dy, d/dz] bf(r) (a total of 4),
// for DerivOrder == 1 and MakeLaplaceToo, these are [1, d/dx, d/dy, d/dz, Laplace] bf(r) (a total of 5).
static size_t GetNumDerivComponents(size_t DerivOrder_, bool MakeLaplaceToo_) {
   assert(DerivOrder_ <= 2);
   return g_nGridDerivComp1[MakeLaplaceToo_][DerivOrder_];
}


struct FEvalGridBfParams
{
   FVec3 vGridBf; // vGridPoint - vBfPos
   double fDistSq;   // norm(vGridBf)^2
   double *pPowR;    // norm(vGridBf)^l for l = 0.. pBf->l.
   size_t MaxPowL;

   size_t DerivOrder;
   size_t nScalarComp; // scalar components. 1 + DerivOrder + (possibly one more for making laplace components).
   double LogThrOrb;
   double ThrOrb;
   double *pExpV; // temporary
   bool MakeLaplaceToo;

   enum {
      GEOM_ScreenStage = 1
   };

   FEvalGridBfParams(size_t nExpMax, size_t nAmMax, size_t nCoMax, double ThrOrb, double LogThrOrb, FMemoryStack &Mem);

   void SetGeometry(double const *pGridPos, FVec3 const &vBfPos, size_t Flags);
   void SetDerivOrder(size_t DerivOrder_, bool MakeLaplaceToo_);
   IR_NO_INLINE void EvalScalar( double *IR_RP pValues, FRawShell const *pBf );
};


FEvalGridBfParams::FEvalGridBfParams(size_t nExpMax_, size_t nAmMax_, size_t /*nCoMax_*/, double ThrOrb_, double LogThrOrb_, FMemoryStack &Mem)
{
   MaxPowL = nAmMax_;
   Mem.Alloc(pExpV, nExpMax_ + (nAmMax_+1));
   pPowR = pExpV + nExpMax_;
   ThrOrb = ThrOrb_;
   LogThrOrb = LogThrOrb_;
   MakeLaplaceToo = false;
}


void FEvalGridBfParams::SetDerivOrder(size_t DerivOrder_, bool MakeLaplaceToo_)
{
   DerivOrder = DerivOrder_;
   MakeLaplaceToo = MakeLaplaceToo_;
//    nComp = GetNumDerivComponents(DerivOrder_, MakeLaplaceToo_);
   nScalarComp = DerivOrder + 1 + int(MakeLaplaceToo);
   assert(DerivOrder_ <= 2);
}


void FEvalGridBfParams::SetGeometry(double const *pGridPos, FVec3 const &vBfPos, size_t Flags)
{
   vGridBf[0] = pGridPos[0] - vBfPos[0];
   vGridBf[1] = pGridPos[1] - vBfPos[1];
   vGridBf[2] = pGridPos[2] - vBfPos[2];
   fDistSq = ct::LengthSq(vGridBf);

   // evaluate r^l. This is the leading term of the largest solid harmonic
   // centered on vBfPos. Need that for screening.
   double
      fDist = std::sqrt(fDistSq);
   pPowR[0] = fDist;
   for ( size_t l = 1; l <= MaxPowL; ++ l )
      pPowR[l] = pPowR[l-1] * fDist;

   if ( 0 != (GEOM_ScreenStage & Flags) ) {
      // these are geometry parameters for the closest grid point, as used for screening.
      // in this case we may /NOT/ remove basis functions due to small r^l values, since
      // r^l will be larger for other grid points!
      if ( pPowR[0] < 1.0 )
         for ( size_t l = 0; l <= MaxPowL; ++ l )
            pPowR[l] = 1.0;
   }
}


void FEvalGridBfParams::EvalScalar(double *IR_RP pValues, FRawShell const *pBf)
{
   //    xout << fmt::format("   !EvalScalar: fDistSq={:8.3f}  LogThr={:8.3f}", fDistSq, LogThrOrb) << std::endl;
   IR_RESUME_CLOCK(23);
   bool
      NothingThere = true;
   size_t const
      nExp = pBf->nExp,
      nScalarComp_ = this->nScalarComp;
   // evaluate the scalar part of the function. First the exponents.
   // here we assume that pPowR can only make the function larger
   // (i.e., the inner parts of high-am functions are not screened out)
   for (size_t iExp = 0; iExp < nExp; ++iExp) {
      double fExp = pBf->pExp[iExp];
      if (s_GridScreeningEnabled && (fDistSq * fExp > LogThrOrb))
         pExpV[iExp] = 0;
      else {
         pExpV[iExp] = std::exp(-fDistSq * fExp);
         NothingThere = false;
      }
   }
   IR_PAUSE_CLOCK(23);
      //    xout << fmt::format("   !EvalScalar: NothingThere? {}", (int)NothingThere) << std::endl;

   IR_RESUME_CLOCK(24);
   if (NothingThere) {
      // all primitives screened out
      for (size_t i = 0; i < nScalarComp_ * pBf->nCo; ++i)
         pValues[i] = 0;
   } else {
      double const
         fPowR = pPowR[pBf->l];
      double const
         *IR_RP pCo = pBf->pCo,
         *IR_RP pExp = pBf->pExp;
      for (size_t iCo = 0; iCo < pBf->nCo; ++iCo, pCo += nExp) {
         // now the total contracted function.
         if (nScalarComp_ == 1) {
            assert(DerivOrder == 0 && !MakeLaplaceToo);
            double
               v0 = 0;
            for (size_t iExp = 0; iExp < nExp; ++iExp)
               v0 += pCo[iExp] * pExpV[iExp];
            if (s_GridScreeningEnabled && (std::abs(v0 * fPowR) < ThrOrb))
               v0 = 0;
            pValues[iCo] = v0;
         } else if (nScalarComp_ == 2) {
            assert((DerivOrder == 1 && !MakeLaplaceToo) || (DerivOrder == 0 && MakeLaplaceToo));
            double v0, v1;
            v0 = 0; // \sum C[iZ] exp(-Z r^2)
            v1 = 0; // \sum C[iZ] (-2Z) exp(-Z r^2)
            // only the scalar parts here.
            for (size_t iExp = 0; iExp != nExp; ++iExp) {
               double
                  e = -2 * pExp[iExp],
                  g = pCo[iExp] * pExpV[iExp];
               v0 += g;
               v1 += e * g;
            }

      //             if ( std::abs(v0 * fPowR) < ThrOrb && std::abs(v1 * fPowR) < ThrOrb ){
            if (s_GridScreeningEnabled && (std::abs(v0 * fPowR) < ThrOrb)) {
               v0 = 0;
               v1 = 0;
            }
            pValues[2 * iCo] = v0; // should be "iCo*nScalarComp_"
            pValues[2 * iCo + 1] = v1;
         } else if (nScalarComp_ == 3) {
            assert((DerivOrder == 2 && !MakeLaplaceToo) || (DerivOrder ==1 && MakeLaplaceToo));
            double v0, v1, v2;
            v0 = 0; // \sum C[iZ] exp(-Z r^2)
            v1 = 0; // \sum C[iZ] (-2Z) exp(-Z r^2)
            v2 = 0; // \sum C[iZ] 4 Z^2 exp(-Z r^2)
            // only the scalar parts here. see below.
            for (size_t iExp = 0; iExp != nExp; ++iExp) {
               double
                  e = -2 * pExp[iExp],
                  g = pCo[iExp] * pExpV[iExp];
               v0 += g;
               v1 += e * g;
               v2 += e * e * g;
            }

            if (s_GridScreeningEnabled && (std::abs(v0 * fPowR) < ThrOrb)) {
               v0 = 0;
               v1 = 0;
               v2 = 0;
            }
            pValues[3 * iCo] = v0;
            pValues[3 * iCo + 1] = v1;
            pValues[3 * iCo + 2] = v2;
         } else {
            assert(nScalarComp_ == 4);
            assert(DerivOrder == 2 && MakeLaplaceToo);
            double v0, v1, v2, v3;
            v0 = 0; // \sum C[iZ] exp(-Z r^2)
            v1 = 0; // \sum C[iZ] (-2Z) exp(-Z r^2)
            v2 = 0; // \sum C[iZ] 4 Z^2 exp(-Z r^2)
            v3 = 0; // \sum C[iZ] (-8z^3) exp(-Z r^2)
            // only the scalar parts here. see below.
            for (size_t iExp = 0; iExp != nExp; ++iExp) {
               double
                  e = -2 * pExp[iExp],
                  g = pCo[iExp] * pExpV[iExp];
               v0 += g;
               v1 += e * g;
               v2 += e * e * g;
               v3 += e * e * e * g;
            }

            if (s_GridScreeningEnabled && (std::abs(v0 * fPowR) < ThrOrb)) {
               v0 = 0;
               v1 = 0;
               v2 = 0;
               v3 = 0;
            }
            pValues[4 * iCo] = v0;
            pValues[4 * iCo + 1] = v1;
            pValues[4 * iCo + 2] = v2;
            pValues[4 * iCo + 3] = v3;
         }
      } // for
   } // else of "if (NothingThere)"
   IR_PAUSE_CLOCK(24);
}


void FindClosestGridPoint(size_t &iGridPtMin, double const *pGridPos, size_t GridStride, size_t nGridPts, FVec3 const &vBfPos)
{
   double
      fDistSqMin = 1e99;
   iGridPtMin = 0xffffffff;
   for ( size_t iGridPt = 0; iGridPt < nGridPts; ++ iGridPt ) {
      FVec3
         vGridPos(pGridPos[iGridPt * GridStride],
                  pGridPos[iGridPt * GridStride + 1],
                  pGridPos[iGridPt * GridStride + 2]);
      double
         fDistSq = ct::LengthSq(vGridPos - vBfPos);
      if ( fDistSq < fDistSqMin ) {
         fDistSqMin = fDistSq;
         iGridPtMin = iGridPt;
      }
   }
}


void FindBfRangeProperties(size_t &nCoSum, size_t &nFnSum, size_t &nCoMax, size_t &nExpMax, size_t &nAmMax, FRawShell const *pFirstBf, FRawShell const *pLastBf )
{
   nCoSum = 0;
   nAmMax = 0;
   nExpMax = 0;
   nCoMax = 0;
   nFnSum = 0; // note: not actually required apart from addressing hack. should be removed!
   for ( FRawShell const *pBf = pFirstBf; pBf != pLastBf; ++ pBf ) {
      nCoSum += pBf->nCo;
      nFnSum += pBf->nCo * (2*pBf->l + 1);
      if ( pBf->nExp > nExpMax )
         nExpMax = pBf->nExp;
      if ( pBf->l > nAmMax )
         nAmMax = pBf->l;
      if ( pBf->nCo > nCoMax )
         nCoMax = pBf->nCo;
   }
}

struct FGridPtBfn
{
   FRawShell const
      *pBf;
   size_t
      iCo,     // contraction index in pBf
      iFnOff,  // offset of *function* (i.e., including spherical components) regarding iFnBase. //regarding pFirstBf
      nCoKept; // number of contractions not screened out (in which ones these are is stored in iCoMask).
   FCoMask // an size_t.
      iCoMask;
};


static void AssembleBf_Deriv0(double *IR_RP pOut0, size_t nBfStride, double const *IR_RP pSlm, double const *IR_RP pVal1, size_t nSh, size_t nCo, FCoMask iCoMask) {
   for ( size_t iCo = 0; iCo < nCo; ++ iCo ) {
      if ( 0 == (iCoMask & (1ul<<iCo)))
         continue;
      for ( size_t iSh = 0; iSh < nSh; ++ iSh )
         pOut0[nBfStride * iSh] = pSlm[iSh] * pVal1[iCo];
      pOut0 += nBfStride * nSh;
   }
}


static void AssembleBf_Deriv1(double *IR_RP pOut0, size_t nBfStride, size_t nCompStride, double const *IR_RP pSlm, double x, double y, double z, double const *IR_RP pVal1, size_t nSh, size_t nCo, FCoMask iCoMask)
{
   //       phi(r) = Slm(r - A)      * \sum C[iZ] exp(-Z (r-A)^2)
   //  d/dx phi(r) = (d/dx Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
   //                Slm(r - A) * \sum C[iZ] (-2*Z*(rx-Ax))* exp(-Z (r-A)^2)
   for ( size_t iCo = 0; iCo < nCo; ++ iCo ) {
      if ( 0 == (iCoMask & (1ul<<iCo)))
         continue;
      for ( size_t iSh = 0; iSh < nSh; ++ iSh ) {
         double
            *IR_RP pOut1 = &pOut0[iSh*nBfStride];
         double const
            *IR_RP pSlm1 = &pSlm[4*iSh];
//          assert(nComp==4);
         double
            fSlm = pSlm1[0],
            fValue = pVal1[2*iCo],
            fSlmD1 = pVal1[2*iCo+1] * fSlm;
         pOut1[0*nCompStride] = fSlm * fValue;
         pOut1[1*nCompStride] = pSlm1[1] * fValue + x * fSlmD1;
         pOut1[2*nCompStride] = pSlm1[2] * fValue + y * fSlmD1;
         pOut1[3*nCompStride] = pSlm1[3] * fValue + z * fSlmD1;
      }
      pOut0 += nBfStride * nSh;
   }
}


static void AssembleBf_Deriv2(double *IR_RP pOut0, size_t nBfStride, size_t nCompStride, double const *IR_RP pSlm, double x, double y, double z, double const *IR_RP pVal1, size_t nSh, size_t nCo, FCoMask iCoMask)
{
   //      phi(r)    = Slm(r - A)      * \sum C[iZ] exp(-Z (r-A)^2)
   // d/dx phi(r)    = (d/dx Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
   //                  Slm(r - A) * {\sum C[iZ] (-2*Z*(rx-Ax))* exp(-Z (r-A)^2)}
   // d^2/dxx phi(r) = (d^2/dxx Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
   //                  2 (d/dx Slm(r-A)) * \sum (-2*Z*(rx-Ax)) C[iZ] exp(-Z (r-A)^2) +
   //                   Slm(r - A) * {\sum C[iZ] ((-2*Z*(rx-Ax))^2 - 2 Z) * exp(-Z (r-A)^2)}
   // d^2/dxy phi(r) = (d^2/dxy Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
   //                  (d/dx Slm(r-A)) * \sum (-2*Z*(ry-Ay)) C[iZ] exp(-Z (r-A)^2) +
   //                  (d/dy Slm(r-A)) * \sum (-2*Z*(rx-Ax)) C[iZ] exp(-Z (r-A)^2) +
   //                   Slm(r - A) * {\sum C[iZ] ((-2*Z*(rx-Ax)) (-2*Z*(ry-Ay))) * exp(-Z (r-A)^2)}
   // again, only the scalar parts here.
   for ( size_t iCo = 0; iCo < nCo; ++ iCo ) {
      if ( 0 == (iCoMask & (1ul<<iCo)))
         continue;
      for ( size_t iSh = 0; iSh < nSh; ++ iSh ) {
         double
            *IR_RP pOut1 = &pOut0[iSh*nBfStride];
         double const
            *IR_RP pSlm1 = &pSlm[10*iSh];
//          assert(nComp==10);
         double
            fSlm = pSlm1[0],
            fSlmDX = pSlm1[1],
            fSlmDY = pSlm1[2],
            fSlmDZ = pSlm1[3],
            fValue = pVal1[3*iCo],
            fDer1 = pVal1[3*iCo+1],
            fDer2 = pVal1[3*iCo+2];
         pOut1[0*nCompStride] = fSlm * fValue;

         pOut1[1*nCompStride] = fSlmDX * fValue + x * fSlm * fDer1; // d/dx
         pOut1[2*nCompStride] = fSlmDY * fValue + y * fSlm * fDer1; // d/dy
         pOut1[3*nCompStride] = fSlmDZ * fValue + z * fSlm * fDer1; // d/dz

         // hmpf. apparently there are two different sets of orders for second derivatives.
         // The code was generated for xx/yy/zz/xy/xz/yz, but this function is supposed
         // to return xx/xy/xz/yy/yz/zz. this is the reason for this re-ordering business here.
         //     4  5  6  7  8  9
         //    xx/xy/xz/yy/yz/zz
         pOut1[4*nCompStride] = pSlm1[4] * fValue + 2 * fSlmDX * x * fDer1 + fSlm * (x * x * fDer2 + fDer1); // d^2/d[xx]
         pOut1[7*nCompStride] = pSlm1[5] * fValue + 2 * fSlmDY * y * fDer1 + fSlm * (y * y * fDer2 + fDer1); // d^2/d[yy]
         pOut1[9*nCompStride] = pSlm1[6] * fValue + 2 * fSlmDZ * z * fDer1 + fSlm * (z * z * fDer2 + fDer1); // d^2/d[zz]

         pOut1[5*nCompStride] = pSlm1[7] * fValue + (x * fSlmDY + y * fSlmDX) * fDer1 + fSlm * y * x * fDer2; // d^2/d[xy]
         pOut1[6*nCompStride] = pSlm1[8] * fValue + (x * fSlmDZ + z * fSlmDX) * fDer1 + fSlm * z * x * fDer2; // d^2/d[xz]
         pOut1[8*nCompStride] = pSlm1[9] * fValue + (y * fSlmDZ + z * fSlmDY) * fDer1 + fSlm * y * z * fDer2; // d^2/d[yz]
      }
      pOut0 += nBfStride * nSh;
   }
}


static void AssembleBf_Deriv1_Laplace(double *IR_RP pOut0, size_t nBfStride, size_t nCompStride, double const *IR_RP pSlm, double x, double y, double z, double const *IR_RP pVal1, size_t nSh, size_t nCo, FCoMask iCoMask)
{
   for (size_t iCo = 0; iCo < nCo; ++iCo) {
      if (0 == (iCoMask & (1ul << iCo)))
         continue;
      for (size_t iSh = 0; iSh < nSh; ++iSh) {
         double
            *IR_RP pOut1 = &pOut0[iSh*nBfStride];
         double const
            *IR_RP pSlm1 = &pSlm[4 * iSh];
         //          assert(nComp==10);
         double
            fSlm = pSlm1[0],
            fSlmDX = pSlm1[1],
            fSlmDY = pSlm1[2],
            fSlmDZ = pSlm1[3],
            fValue = pVal1[3 * iCo],
            fDer1 = pVal1[3 * iCo + 1],
            fDer2 = pVal1[3 * iCo + 2];
         pOut1[0 * nCompStride] = fSlm * fValue;

         pOut1[1 * nCompStride] = fSlmDX * fValue + x * fSlm * fDer1; // d/dx
         pOut1[2 * nCompStride] = fSlmDY * fValue + y * fSlm * fDer1; // d/dy
         pOut1[3 * nCompStride] = fSlmDZ * fValue + z * fSlm * fDer1; // d/dz

         pOut1[4 * nCompStride] = 2 * fDer1 * (fSlmDX*x + fSlmDY*y + fSlmDZ*z) + fSlm * ((x*x + y*y + z*z) * fDer2 + 3 * fDer1); // laplace
      }
      pOut0 += nBfStride * nSh;
   }
}


static void AssembleBf_Deriv2_Laplace(double *IR_RP pOut0, size_t nBfStride, size_t nCompStride, double const *IR_RP pSlm, double x, double y, double z, double const *IR_RP pVal1, size_t nSh, size_t nCo, FCoMask iCoMask)
{
   for (size_t iCo = 0; iCo < nCo; ++iCo) {
      if (0 == (iCoMask & (1ul << iCo)))
         continue;
      for (size_t iSh = 0; iSh < nSh; ++iSh) {
         double
            *IR_RP pOut1 = &pOut0[iSh*nBfStride];
         double const
            *IR_RP pSlm1 = &pSlm[10 * iSh];
         //          assert(nComp==10);
         double
            fSlm = pSlm1[0],
            fSlmDX = pSlm1[1],
            fSlmDY = pSlm1[2],
            fSlmDZ = pSlm1[3],
            fSlmDXX = pSlm1[4],
            fSlmDYY = pSlm1[5],
            fSlmDZZ = pSlm1[6],
            fSlmDXY = pSlm1[7],
            fSlmDXZ = pSlm1[8],
            fSlmDYZ = pSlm1[9],
            fValue = pVal1[4 * iCo],
            fDer1 = pVal1[4 * iCo + 1],
            fDer2 = pVal1[4 * iCo + 2],
            fDer3 = pVal1[4 * iCo + 3];

         pOut1[0 * nCompStride] = fSlm * fValue;

         pOut1[1 * nCompStride] = fSlmDX * fValue + x * fSlm * fDer1; // d/dx
         pOut1[2 * nCompStride] = fSlmDY * fValue + y * fSlm * fDer1; // d/dy
         pOut1[3 * nCompStride] = fSlmDZ * fValue + z * fSlm * fDer1; // d/dz

         pOut1[4 * nCompStride] = pSlm1[4] * fValue + 2 * fSlmDX * x * fDer1 + fSlm * (x * x * fDer2 + fDer1); // d^2/d[xx]
         pOut1[7 * nCompStride] = pSlm1[5] * fValue + 2 * fSlmDY * y * fDer1 + fSlm * (y * y * fDer2 + fDer1); // d^2/d[yy]
         pOut1[9 * nCompStride] = pSlm1[6] * fValue + 2 * fSlmDZ * z * fDer1 + fSlm * (z * z * fDer2 + fDer1); // d^2/d[zz]

         pOut1[5 * nCompStride] = pSlm1[7] * fValue + (x * fSlmDY + y * fSlmDX) * fDer1 + fSlm * y * x * fDer2; // d^2/d[xy]
         pOut1[6 * nCompStride] = pSlm1[8] * fValue + (x * fSlmDZ + z * fSlmDX) * fDer1 + fSlm * z * x * fDer2; // d^2/d[xz]
         pOut1[8 * nCompStride] = pSlm1[9] * fValue + (y * fSlmDZ + z * fSlmDY) * fDer1 + fSlm * y * z * fDer2; // d^2/d[yz]

         double
            rsq = x*x + y*y + z*z,
            dSlm_dot_xyz_2 = 2 * (x*fSlmDX + y*fSlmDY + z*fSlmDZ),
            sd0 = fDer2 * dSlm_dot_xyz_2,
            sd1 = (rsq * fDer2 + 5 * fDer1),
            sd2 = (rsq * fDer3 + 5 * fDer2) * fSlm;
         pOut1[10 * nCompStride] = fDer1 * dSlm_dot_xyz_2 + fSlm * (rsq * fDer2 + 3 * fDer1); // laplace
         pOut1[11 * nCompStride] = x*sd0 + 2*fDer1*(x*fSlmDXX + y*fSlmDXY + z*fSlmDXZ) + sd1*fSlmDX + sd2*x; // d/dx laplace
         pOut1[12 * nCompStride] = y*sd0 + 2*fDer1*(x*fSlmDXY + y*fSlmDYY + z*fSlmDYZ) + sd1*fSlmDY + sd2*y; // d/dy laplace
         pOut1[13 * nCompStride] = z*sd0 + 2*fDer1*(x*fSlmDXZ + y*fSlmDYZ + z*fSlmDZZ) + sd1*fSlmDZ + sd2*z; // d/dz laplace
      }
      pOut0 += nBfStride * nSh;
   }
}


void AssembleBfs(double *IR_RP pOut0, size_t nBfStride, size_t nCompStride, double const *IR_RP pSlmA, double x, double y, double z, double const *IR_RP pVal1, size_t nScalarComp, FGridPtBfn const *pFirstKeptBf, FGridPtBfn const *pLastKeptBf, size_t DerivOrder, bool MakeLaplaceToo)
{
   size_t
      // this never has laplace components -- only pure derivative components (unlike the output derivatives of the basis functions).
      nCompSlm = GetNumDerivComponents(DerivOrder, false);
   for ( FGridPtBfn const *pKeptBf = pFirstKeptBf; pKeptBf != pLastKeptBf; ++ pKeptBf ){
      FRawShell const
         *pBf = pKeptBf->pBf;
      int
         la = (int)pBf->l;
      size_t
         nCo = pBf->nCo,
         nSh = 2*la + 1;
      double const
         // base of current basis function components in pSlmA
         *IR_RP pSlm = &pSlmA[nCompSlm * nSlmX(la-1)];

      // TODO: here set everything to exactly zero if it lies below
      // the target threshold?
      if (!MakeLaplaceToo) {
         if (DerivOrder == 0) {
            AssembleBf_Deriv0(pOut0, nBfStride, pSlm, pVal1, nSh, nCo, pKeptBf->iCoMask);
         } else if (DerivOrder == 1) {
            AssembleBf_Deriv1(pOut0, nBfStride, nCompStride, pSlm, x, y, z, pVal1, nSh, nCo, pKeptBf->iCoMask);
         } else if (DerivOrder == 2) {
            AssembleBf_Deriv2(pOut0, nBfStride, nCompStride, pSlm, x, y, z, pVal1, nSh, nCo, pKeptBf->iCoMask);
         } else {
            assert(0);
         }
      } else {
         if (DerivOrder == 0) {
            assert(0); // AssembleBf_Deriv0(pOut0, nBfStride, pSlm, pVal1, nSh, nCo, pKeptBf->iCoMask);
         } else if (DerivOrder == 1) {
            AssembleBf_Deriv1_Laplace(pOut0, nBfStride, nCompStride, pSlm, x, y, z, pVal1, nSh, nCo, pKeptBf->iCoMask);
         } else if (DerivOrder == 2) {
            AssembleBf_Deriv2_Laplace(pOut0, nBfStride, nCompStride, pSlm, x, y, z, pVal1, nSh, nCo, pKeptBf->iCoMask);
         } else {
            assert(0);
         }
      }

      pOut0 += nBfStride * pKeptBf->nCoKept * nSh;
      pVal1 += nScalarComp * nCo;
   }
}


// Evaluate the nGridPts x #bf matrix of basis functions centered at vBfPos
// on the given grid points.
//    pOut[iPt + nCompStride * iDiff + nBfStride * iBf] = (value).
//        iDiff: derivative component.
//        iBf: basis function index starting at pFirstBf
void EvalShellGroupOnGrid( double *pOut, size_t nCompStride, size_t nBfStride,
   size_t *pMap, size_t &nMap, size_t iFnBase,
   double const *pGridPos, size_t GridStride, size_t nGridPts,
   FRawShell const *pFirstBf, FRawShell const *pLastBf, double const *pBfPos, // <- namespace names for doxygen only
   unsigned DerivOrder, bool MakeLaplaceToo, double ThrOrb, double LogThrOrb,
   FMemoryStack &Mem )
{
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   FVec3
      vBfPos(pBfPos[0], pBfPos[1], pBfPos[2]);

   // determine number of derivative components.
   assert(DerivOrder <= 2);
   size_t
//       // number of derivative components we are supposed to compute
//       nCompOut = GetNumDerivComponents(DerivOrder, MakeLaplaceToo),
      // number of derivative components in the intermediate Slm(r-Ra) solid
      // harmonic objects (these may differ from nCompOut because they have no
      // laplace components)
      nCompSlm = GetNumDerivComponents(DerivOrder, false);

   // hm... I think TECHNICALLY the functions can output stuff in transposed
   // order (i.e., nGridPt x nAo x nComp) with suitable strides. Didn't try,
   // though.
   assert(nCompStride >= nGridPts);
   assert(nBfStride >= nCompStride * GetNumDerivComponents(DerivOrder, MakeLaplaceToo));

   // find grid point closest to the current basis function center
   size_t
      iGridPtMin;
   FindClosestGridPoint(iGridPtMin, pGridPos, GridStride, nGridPts, vBfPos);


   // initialize data structure for evaluating scalar parts of basis functions
   size_t
      nCoSum, nFnSum;
   size_t
      nCoMax, nAmMax, nExpMax;
   FindBfRangeProperties(nCoSum, nFnSum, nCoMax, nExpMax, nAmMax, pFirstBf, pLastBf);


   FEvalGridBfParams
      P(nExpMax, nAmMax, nCoMax, ThrOrb, LogThrOrb, Mem);
   double
      *pValues;
   size_t
      nScalarComp = (1u + DerivOrder + size_t(MakeLaplaceToo));
   Mem.Alloc(pValues, nCoSum * nScalarComp);

   // make list of basis functions we need to keep
   FGridPtBfn
      *pFirstKeptBf, *pLastKeptBf;
   Mem.Alloc(pFirstKeptBf, nCoSum);
   pLastKeptBf = pFirstKeptBf;
   P.SetDerivOrder(0, false); // only actual function values for screening.
   P.SetGeometry(&pGridPos[iGridPtMin * GridStride], vBfPos, FEvalGridBfParams::GEOM_ScreenStage);

   size_t
      iFnOff = iFnBase;
   int
      MaxAmKept = -1;
   for (FRawShell const *pBf = pFirstBf; pBf != pLastBf; iFnOff += pBf->nCo * (2*pBf->l + 1), ++ pBf) {
      P.EvalScalar(pValues, pBf);
      FCoMask
         iCoMask = 0;
      size_t
         nCoKept = 0;
      size_t
         nSh = 2*pBf->l + 1;
      for (size_t iCo = 0; iCo < pBf->nCo; ++ iCo) {
         if (!s_GridScreeningEnabled || (pValues[iCo] != 0)) {
            iCoMask |= (size_t(1) << iCo);
            nCoKept += 1;

            for (size_t iSh = 0; iSh < nSh; ++ iSh) {
               pMap[nMap] = iFnOff + iSh + iCo * nSh;
               ++ nMap;
            }
         }
      }

      if (iCoMask != 0) { // at least one non-vanishing contraction. this function must be kept.
         pLastKeptBf->pBf = pBf;
         pLastKeptBf->iFnOff = iFnOff;
         pLastKeptBf->iCoMask = iCoMask;
         pLastKeptBf->nCoKept = nCoKept;
         if ( (int)pBf->l > MaxAmKept )
            MaxAmKept = (int)pBf->l;
         ++ pLastKeptBf;
      }
   }

   size_t
      nSlmA = nSlmX(MaxAmKept),
      nSlmA1 = nCompSlm * nSlmA;
   double
      *pSlmA;
   Mem.Alloc(pSlmA, nSlmA1);
   P.SetDerivOrder(DerivOrder, MakeLaplaceToo);
   assert(P.nScalarComp == nScalarComp);

   for (size_t iGridPt = 0; iGridPt < nGridPts; ++ iGridPt) {
      // hmhm... Maybe try a vectorization with fixed packet length along the grid point dimension?
      // all these things do ALMOST the same thing for all the grid points.
      // (apart from nMaxAmOnPt possibly differing, but for this we could just take
      // the max over the packet. Note that iCoMask is the same for all grid points.).
      // Might be worth templatizing FEvalGridBfParams on the value types (e.g., for double
      // and Vector4d) and doing the SlmX thing separately for a start.
      IR_RESUME_CLOCK(20);
      P.SetGeometry(&pGridPos[iGridPt * GridStride], vBfPos, 0);

      int
         nMaxAmOnPt = -1;
      // evaluate the scalar parts of the functions.
      double
         *pVal1 = pValues;
      for (FGridPtBfn const *pKeptBf = pFirstKeptBf; pKeptBf != pLastKeptBf; ++ pKeptBf) {
         FRawShell const *pBf = pKeptBf->pBf;
         P.EvalScalar(pVal1, pBf);
         if ((int)pBf->l > nMaxAmOnPt)
            nMaxAmOnPt = (int)pBf->l;
         pVal1 += P.nScalarComp * pBf->nCo;
      }
      IR_PAUSE_CLOCK(20);

      if (nMaxAmOnPt == -1) {
         // everything screened out.
         continue;
      }

      IR_RESUME_CLOCK(21);
      // evaluate the angular part of the function
      double const
         x = P.vGridBf[0],
         y = P.vGridBf[1],
         z = P.vGridBf[2];
      if (DerivOrder == 0)
         EvalSlcX_Deriv0(pSlmA, x, y, z, nMaxAmOnPt);
      else if (DerivOrder == 1)
         EvalSlcX_Deriv1(pSlmA, x, y, z, nMaxAmOnPt);
      else if (DerivOrder == 2)
         EvalSlcX_Deriv2(pSlmA, x, y, z, nMaxAmOnPt);
      else
         assert(0);
      IR_PAUSE_CLOCK(21);


      // multiply scalar and angular parts together to get the final grid values.
      IR_RESUME_CLOCK(22);
      AssembleBfs(&pOut[iGridPt], nBfStride, nCompStride, pSlmA, x, y, z, pValues, P.nScalarComp, pFirstKeptBf, pLastKeptBf, DerivOrder, MakeLaplaceToo);
      IR_PAUSE_CLOCK(22);
   }
}


} // namespace ir
