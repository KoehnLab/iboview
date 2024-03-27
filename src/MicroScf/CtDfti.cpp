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

#include <iostream> // FIXME: remove this.
#include "CxAlgebra.h"
#include "CtDfti.h"
#include "Ir.h"
#include "CtMatrix.h"
#include "CxOpenMpProxy.h"
#include "CxOpenMpAcc.h"
using namespace ct;


// TODO:
// - check out if Dot2, Add2, DotN, AccN etc. actually help or make things worse
//   with current compilers. Nowadays g++ can probably do a better job if the
//   partial unrolling is not done, and instead we find a way of telling it that
//   the matrices are aligned.
//
// - (additionally, we probably should try reformulating these ops from inner-
//   product based operations to outer-product based operations. Those work much
//   better on modern hardware)
//
// - FIXME: ...in particular for the Mxms and Mxvs! They are all TN at the moment! Oh my.


// #ifdef _MSC_VER
//    // OpenBLAS DGEMV crashes on me if executed in OpenMP-loop. No idea.
//    #pragma message( "WARNING: REPLACING REAL DGEMV BY MxmLame AS OPENBLAS/OPENMP WORKAROUND!")
//    #define Mxv MxvLame
// #endif
// ^- disabled since I currently mostly use eigen instead of OpenBLAS for VC++ builds.

#if 0
inline void MXMA(double const *pA, size_t iRowStA, size_t iColStA, double const *pB, size_t iRowStB, size_t iColStB, double *pOut, size_t iRowStOut, size_t iColStOut, size_t nRows, size_t nLink, size_t nCols) {
   FMatrixView const
      A(const_cast<double*>(pA), nRows, nLink, iRowStA, iColStA),
      B(const_cast<double*>(pB), nLink, nCols, iRowStB, iColStB);
   FMatrixView
      Out(pOut, nRows, nCols, iRowStOut, iColStOut);
   Mxm(Out, A, B);
}

inline void MXMB(double const *pA, size_t iRowStA, size_t iColStA, double const *pB, size_t iRowStB, size_t iColStB, double *pOut, size_t iRowStOut, size_t iColStOut, size_t nRows, size_t nLink, size_t nCols) {
   FMatrixView const
      A(const_cast<double*>(pA), nRows, nLink, iRowStA, iColStA),
      B(const_cast<double*>(pB), nLink, nCols, iRowStB, iColStB);
   FMatrixView
      Out(pOut, nRows, nCols, iRowStOut, iColStOut);
   Mxm(Out, A, B, 1.0, MXM_Add);
}
#else

IR_FORCE_INLINE void MXMA(double const *pA, size_t iRowStA, size_t iColStA, double const *pB, size_t iRowStB, size_t iColStB, double *pOut, size_t iRowStOut, size_t iColStOut, size_t nRows, size_t nLink, size_t nCols)
{
   return Mxm_InlineCalls(pOut, iRowStOut, iColStOut,
         pA, iRowStA, iColStA, pB, iRowStB, iColStB,
         nRows, nLink, nCols, false, 1.0);
}

IR_FORCE_INLINE void MXMB(double const *pA, size_t iRowStA, size_t iColStA, double const *pB, size_t iRowStB, size_t iColStB, double *pOut, size_t iRowStOut, size_t iColStOut, size_t nRows, size_t nLink, size_t nCols)
{
   return Mxm_InlineCalls(pOut, iRowStOut, iColStOut,
         pA, iRowStA, iColStA, pB, iRowStB, iColStB,
         nRows, nLink, nCols, true, 1.0);
}

#endif


namespace dfti {

using mig::FDftGrid;

enum {
   CLOCK_EvalBfn = 0,
   CLOCK_FormRho = 1,
   CLOCK_FormXcMat = 2,
   CLOCK_EvalDftfun = 3,
   CLOCK_GridWtDeriv = 4, // for gradient code.
   CLOCK_CoreGrad = 5,
   CLOCK_TimerCount = 6
};

enum {
   iNoCenter = 0xffffffff
};

struct FDftiResultSet;
struct FDensitySet;
struct FDensitySetInput;





// #define RESUME_CPU_CLOCK(iClockId)
// #define PAUSE_CPU_CLOCK(iClockId)

// note: comment with /* */ to avoid multi-line comment warnings.

#define RESUME_CPU_CLOCK(iClockId) \
   if ( Args.MeasureTime() ) \
      Timings[(iClockId)] -= Second()

#define PAUSE_CPU_CLOCK(iClockId) \
   if ( Args.MeasureTime() ) \
      Timings[(iClockId)] += Second()


struct FDensitySetInput
{
   double
      *pRdm, *pOccOrb;
   size_t
      nBf, nOcc;
   size_t
      *pMap, nMap;
   bool
      UseOrb,
      Gradient; // keep data for gradient evaluation
   size_t
      nComp,
      iLaplaceComp0;
   FEvalBfDerivType
      DerivTypeBf;
   FDensitySetInput() {};
   FDensitySetInput(double *pRdm_, double *pOccOrb_, size_t nBf_, size_t nOcc_,
         size_t *pMap_, size_t nMap_, bool UseOrb_, bool Gradient_, size_t nComp_, size_t iLaplaceComp0_, FEvalBfDerivType DerivTypeBf_)
      : pRdm(pRdm_), pOccOrb(pOccOrb_), nBf(nBf_), nOcc(nOcc_), pMap(pMap_), nMap(nMap_),
        UseOrb(UseOrb_), Gradient(Gradient_), nComp(nComp_), iLaplaceComp0(iLaplaceComp0_), DerivTypeBf(DerivTypeBf_)
   {};
};

// densities for one spin component.
struct FDensitySet
{
   double
      *pRho,      // electron density: rho(rg) = \sum_i [phi_i(rg)]^2 (len: nPt)
      *pRhoGrd,   // gradient of electron density: d/dx_i rho(rg)  (len: 3*nPt)
      *pTau,      // kinetic term: \sum_i dot(grad phi_i(rg), grad phi_i(rg)) (len: nPt)
      *pUpsilon,  // laplacian of electron density: div grad rho(rg) (len: nPt)
      *pBxRho;    // may be 0. Otherwise pOrbVal contracted with density matrix: arrays nGridPt x nBx x nMap.
   size_t
      nGridPt;
   uint
      nDiff,
      nBx; // number of components of pOrbVal stored in pBxRho (1: only density; 4: density+gradient)
   bool
      MakeTau,
      MakeUpsilon,
      AuxExpandXc;
   double
      fElec;

   void Init(size_t nGridPt_, uint nDiff_, bool MakeTau_, bool MakeUpsilon_, bool AuxExpandXc_, FMemoryStack &Mem);
   // evaluate densities (+gradients)
   void Eval(double const *pOrbVal, FDensitySetInput const &p, bool ForceNonNegative, double const *pGridWt, FMemoryStack &Mem);
   void EvalAux(double const *pOrbVal, FDensitySetInput const &p, bool ForceNonNegative, double const *pGridWt, FMemoryStack &Mem);


   // calculate and accumulate xc matrix contribution.
   void AccXcMatrix(double *pXcTriang, double *pOrbVal, FDftiResultSet const &r,
      double *pRhoGrdY, FDensitySetInput const &p, double const *pGridWt, FMemoryStack &Mem);
   void AccXcMatrixAux(double *pXcVec, double *pOrbVal, FDftiResultSet const &r,
      double *pRhoGrdY, FDensitySetInput const &p, double const *pGridWt, FMemoryStack &Mem);

   // accumulate symmetric contribution mu(r) vdrho(r) nu(r) to given xc matrix.
   // beware of mysterious flags.
   void AccXcSym1(double *pXc, double const *pVdRho, double const *pGridWt, double const *pOrbVal, size_t nMapSt, size_t nMap, uint Flags, FMemoryStack &Mem);

   double *MakePhi(double const *pIn, size_t nPts, size_t Stride, FDensitySetInput const &p, FMemoryStack &Mem) const;
   double *MakeBxRho(double const *pIn, size_t nPts, size_t Stride, FDensitySetInput const &p, FMemoryStack &Mem) const;

   double *CompressAuxVec(double const *pInput, FDensitySetInput const &p, FMemoryStack &Mem) const;
};

void FDensitySet::Init(size_t nGridPt_, uint nDiff_, bool MakeTau_, bool MakeUpsilon_, bool AuxExpandXc_, FMemoryStack &Mem)
{
   nGridPt = nGridPt_;
   nDiff = nDiff_;
   MakeTau = MakeTau_;
   MakeUpsilon = MakeUpsilon_;
   AuxExpandXc = AuxExpandXc_;

   pBxRho = 0;
   nBx = 0;

   // not required -- just here to simplify debugging.
   pRhoGrd = 0;
   pTau = 0;
   pUpsilon = 0;

   Mem.Alloc(pRho, nGridPt);
   if (nDiff >= 1)
      Mem.Alloc(pRhoGrd, 3 * nGridPt);
   if (nDiff >= 2 || MakeTau )
      Mem.Alloc(pTau, nGridPt);
   if (nDiff >= 2 || MakeUpsilon)
      Mem.Alloc(pUpsilon, nGridPt);
}

//#ifdef INCLUDE_OPTIONALS
double *FDensitySet::MakeBxRho(double const *pIn, size_t nPts, size_t Stride, FDensitySetInput const &p, FMemoryStack &Mem) const
{
   // o form contraction pOut[iPt,iMap] = pDensity[iMap,jMap] pIn[iPt,jMap].
   //   Output is nPoints x nMap array, allocated on Mem.
   // o pIn is addressed as pIn[iPt,iMap] = iPt + nStride * iMap.
   // o If UseOrb is true, the contaction is formed as a two-step process
   //   over the provided the occupied orbitals. Otherwise the density
   //   in pRdm is used.
   // o what does Bx stand for? I don't know, but the Fortran code called it like that.

   double
      *pOut;
   size_t
      nMap = p.nMap;
   size_t
      *pMap = p.pMap;
   Mem.Alloc(pOut, nPts * nMap);

   if (!p.UseOrb) {
      double
         *pRdm1;

      // unpack input density matrix from triangular into square, and compress
      // from (triangular) nBf x nBf to (square) nMap x nMap dimension.
      Mem.Alloc(pRdm1, nMap * nMap);
      for (size_t iMapCol = 0; iMapCol < nMap; ++ iMapCol) {
         size_t
            iCol = pMap[iMapCol],
            iColOff = iCol*(iCol+1)/2; // offset of i'th col in output packed triangular matrix
         for (size_t iMapRow = 0; iMapRow <= iMapCol; ++ iMapRow) {
            size_t iRow = pMap[iMapRow];
            double f = p.pRdm[iRow + iColOff];
            pRdm1[iMapRow + nMap * iMapCol] = f;
            pRdm1[iMapCol + nMap * iMapRow] = f;
         }
      }

      // multiply nMap dimension of input data by compressed density.
      MXMA(pIn,1,Stride,  pRdm1,1,nMap,   pOut,1,nPts,  nPts,nMap,nMap);

//       ir::PrintMatrixGen(ct::xout, pIn, nPts, 1, nMap, Stride, "OrbitalValues[nGridPt x nMap]", "{:12.5f}");
//       ir::PrintMatrixGen(ct::xout, pRdm1, nMap, 1, nMap, nMap, "RDM1[nMap x nMap]", "{:12.5f}");
//       ir::PrintMatrixGen(ct::xout, pOut, nPts, 1, nMap, Stride, "Out: OrbVal x Rdm", "{:12.5f}");

      Mem.Free(pRdm1);
   } else {
      double
         *pOrb1;
      // compress input orbital matrix from nOcc x nBf into nOcc x nMap.
      size_t
         nOcc = p.nOcc;
      Mem.Alloc(pOrb1, nOcc * nMap);
      for (size_t iMap = 0; iMap < nMap; ++ iMap) {
         size_t
            iBf = pMap[iMap];
         for (size_t iOcc = 0; iOcc < nOcc; ++ iOcc)
            pOrb1[iOcc + nOcc * iMap] = p.pOccOrb[iOcc + nOcc * iBf];
      }

      // transform nMap dimension of input to nOcc. Then transform again to
      // arrive at inverse nMap.
      double
         *pDummy;
      Mem.Alloc(pDummy, nPts * nOcc);
      MXMA(pIn,1,Stride,  pOrb1,nOcc,1,  pDummy,1,nPts,  nPts,nMap,nOcc);
      MXMA(pDummy,1,nPts,  pOrb1,1,nOcc,  pOut,1,nPts,  nPts,nOcc,nMap);

      Mem.Free(pOrb1);
   }

   return pOut;
}

// transforms input [nGridPt x nComp] x nMap into nOcc x [nGridPt x nComp]
double *FDensitySet::MakePhi(double const *pIn, size_t nPts, size_t nStride, FDensitySetInput const &p, FMemoryStack &Mem) const
{
   double
      *pOut;
   size_t
      nOcc = p.nOcc;
   Mem.Alloc(pOut, nOcc * nPts);

   double
      *pOrb1;
   // compress input orbital matrix from nOcc x nBf into nOcc x nMap.
   Mem.Alloc(pOrb1, nOcc * p.nMap);
   for (size_t iMap = 0; iMap < (size_t)p.nMap; ++ iMap) {
      size_t
         iBf = p.pMap[iMap];
      for (size_t iOcc = 0; iOcc < nOcc; ++ iOcc)
         pOrb1[iOcc + nOcc * iMap] = p.pOccOrb[iOcc + nOcc * iBf];
   }

   // transform nMap dimension of input to nOcc, with inline transpose of output.
   MXMA(pOrb1,1,nOcc,  pIn,nStride,1,  pOut,1,nOcc,  nOcc,p.nMap,nPts);

   Mem.Free(pOrb1);
   return pOut;
}
//#endif // INCLUDE_OPTIONALS


enum FAlgOpFlags{
   ALGOP_Add = 0x01,  // add to output instead of replacing it.
   ALGOP_Symmetrize = 0x02 // copy lower triangle to upper triangle.
};


// r = dot(x,y)
template<class FScalar>
FScalar Dot2(FScalar const *IR_RP x, FScalar const *IR_RP y, size_t n)
{
   FScalar
      r = 0;
   size_t
      i = 0;
   for ( ; i < (n & (~3)); i += 4 ) {
      r += x[i]   * y[i];
      r += x[i+1] * y[i+1];
      r += x[i+2] * y[i+2];
      r += x[i+3] * y[i+3];
   }
   for ( ; i < n; ++ i ) {
      r += x[i] * y[i];
   }
   return r;
}


template
double Dot2<double>(double const *IR_RP x, double const *IR_RP y, size_t n);

// r[i] += f * x[i]
void Add2(double *IR_RP r, double const *IR_RP x, double f, std::size_t n)
{
   std::size_t
      i = 0;
   for ( ; i < (n & (~3)); i += 4) {
      r[i]   += f * x[i];
      r[i+1] += f * x[i+1];
      r[i+2] += f * x[i+2];
      r[i+3] += f * x[i+3];
   }
   for ( ; i < n; ++ i) {
      r[i] += f * x[i];
   }
}


// r[i] += f * x[i] * y[i]
void Add2(double *IR_RP r, double const *IR_RP x, double const *IR_RP y, double f, std::size_t n)
{
   std::size_t
      i = 0;
   for ( ; i < (n & (~3)); i += 4) {
      r[i]   += f * x[i]   * y[i];
      r[i+1] += f * x[i+1] * y[i+1];
      r[i+2] += f * x[i+2] * y[i+2];
      r[i+3] += f * x[i+3] * y[i+3];
   }
   for ( ; i < n; ++ i) {
      r[i] += f * x[i] * y[i];
   }
}



// ddots multiple strided sets of vectors.
void DotN(double *pOut, double Factor, double const *pA, size_t nStrideA, double const *pB, size_t nStrideB, size_t nPoints, size_t nSets, uint Flags = 0)
{
   for (size_t iSet = 0; iSet < nSets; ++ iSet) {
      double
         f = Factor * Dot2(&pA[iSet * nStrideA], &pB[iSet * nStrideB], nPoints);

      if (Flags == 0)
         pOut[iSet] = f;
      else {
         assert(Flags == ALGOP_Add);
         pOut[iSet] += f;
      }
   }
}

// daxpys multiple strided sets of vectors.
void AccN(double *pOut, double Factor, double const *pA, size_t nStrideA, double const *pB, size_t nStrideB, size_t nPoints, size_t nSets, uint Flags = 0)
{
   if ((Flags & ALGOP_Add) == 0)
      // clear output.
      for (size_t iPt = 0; iPt < nPoints; ++ iPt)
         pOut[iPt] = 0;

   for (size_t iSet = 0; iSet < nSets; ++ iSet)
      Add2(pOut, &pA[iSet * nStrideA], &pB[iSet * nStrideB], Factor, nPoints);
}


static double Dot3(double const *IR_RP x, double const *IR_RP y, double const *IR_RP z, size_t n)
{
   double
      r = 0;
   for (size_t i = 0; i != n; ++ i)
      r += x[i] * y[i] * z[i];
   return r;
}


void FDensitySet::Eval(double const *pOrbVal, FDensitySetInput const &p, bool ForceNonNegative, double const *pGridWt, FMemoryStack &Mem)
{
   if (AuxExpandXc)
      return EvalAux(pOrbVal, p, ForceNonNegative, pGridWt, Mem);
//#ifdef INCLUDE_OPTIONALS
   void
      *pBeginOfStorage = Mem.Alloc(0); // note: NOT freed in gradient case! (for keeping bxrho!)
   size_t
      nComp = p.nComp;
   size_t
      nMap = (size_t)p.nMap,
      nOcc = p.nOcc,
      nMapSt = nGridPt * nComp; // stride between two basis function entries in OrbVal.
   if (nDiff == 0) {
      // LDA.
      if (p.UseOrb && !p.Gradient) { // for gradient case: we want to keep bxrho for later.
         // evaluate       rho(r) = \sum_i phi_i(r) phi_i(r)
         double
            *pPhi = MakePhi(pOrbVal, nGridPt, nGridPt * nComp, p, Mem);
         DotN(pRho, 1.0, pPhi, p.nOcc, pPhi, p.nOcc, p.nOcc, nGridPt);
      } else {
         pBxRho = MakeBxRho(pOrbVal, nGridPt, nGridPt * nComp, p, Mem);
         nBx = 1;
         AccN(pRho, 1.0, pOrbVal, nMapSt, pBxRho, nGridPt, nGridPt, nMap);
      }
   } else if (nDiff == 1 && !(MakeTau || MakeUpsilon)) {
      // GGA without tau. Most common case.
      nBx = 1; // contract only density
      if (p.Gradient)
         nBx = 4; // contract density & gradients.
      pBxRho = MakeBxRho(pOrbVal, nGridPt * nBx, nGridPt * nComp, p, Mem);
      AccN(pRho, 1.0, &pOrbVal[0], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
      AccN(&pRhoGrd[0 * nGridPt], 2.0, &pOrbVal[1 * nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
      AccN(&pRhoGrd[1 * nGridPt], 2.0, &pOrbVal[2 * nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
      AccN(&pRhoGrd[2 * nGridPt], 2.0, &pOrbVal[3 * nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
   } else if (nDiff == 1 && (MakeTau || MakeUpsilon)) {
      // GGA with tau.
      // evaluate       rho(r) = \sum_i phi_i(r) phi_i(r)
      //          d/dx  rho(r) = \sum_i [d/dx phi_i(r)] phi_i(r) + phi_i(r) [d/dx phi_i(r)]
      //                       = \sum_i 2 phi_i(r) d/dx phi_i(r)
      //                tau(r) = (1/2) \sum_i [d/dx phi_i(r)]^2 + [d/dy phi_i(r)]^2 + [d/dz phi_i(r)]^2
      //            upsilon(r) = [D/Dx^2 + D/Dy^2 + D/Dz^2] rho
      //                       = \sum_i Del Del (phi_i(r) phi_i(r))
      //                       = \sum_i Del ([Del phi_i(r)] phi_i(r) + phi_i(r) [Del phi_i(r)])
      //                       = \sum_i 2 Del (phi_i(r) Del phi_i(r))
      //                       = \sum_i 2 [Del phi_i(r)] [Del phi_i(r)] + 2 phi_i(r) [Del^2 phi_i(r)]
      //                       = 4 tau(r) + 2 \sum_i phi_i(r) [Del^2 phi_i(r)]
      if (p.UseOrb && !p.Gradient && !MakeUpsilon) {
         double
            // transform both bf values and bf gradients to occupied orbitals.
            // output is matrix nOcc x nGridPt x 4.
            *pPhi = MakePhi(pOrbVal, nGridPt*4, nGridPt * nComp, p, Mem);
         size_t
            n = nOcc * nGridPt;
         DotN(pRho, 1.0, &pPhi[0], nOcc, &pPhi[0], nOcc, nOcc, nGridPt);
         DotN(&pRhoGrd[0*nGridPt], 2.0, &pPhi[1*n], nOcc, &pPhi[0], nOcc, nOcc, nGridPt);
         DotN(&pRhoGrd[1*nGridPt], 2.0, &pPhi[2*n], nOcc, &pPhi[0], nOcc, nOcc, nGridPt);
         DotN(&pRhoGrd[2*nGridPt], 2.0, &pPhi[3*n], nOcc, &pPhi[0], nOcc, nOcc, nGridPt);

         // this forms tau(rg) = (1/2) \sum_{i in occ} [d/dx_j phi(rg)]^2.
         DotN(pTau, 0.5, &pPhi[1*n], nOcc, &pPhi[1*n], nOcc, nOcc, nGridPt); // X
         DotN(pTau, 0.5, &pPhi[2*n], nOcc, &pPhi[2*n], nOcc, nOcc, nGridPt, ALGOP_Add); // Y
         DotN(pTau, 0.5, &pPhi[3*n], nOcc, &pPhi[3*n], nOcc, nOcc, nGridPt, ALGOP_Add); // Z
      } else {
         nBx = 4; // contract density & gradients.
         // transform bf values and bf gradients with density matrix.
         pBxRho = MakeBxRho(pOrbVal, nGridPt * nBx, nGridPt * nComp, p, Mem);
         AccN(pRho, 1.0, &pOrbVal[0], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
         AccN(&pRhoGrd[0*nGridPt], 2.0, &pOrbVal[1*nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
         AccN(&pRhoGrd[1*nGridPt], 2.0, &pOrbVal[2*nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
         AccN(&pRhoGrd[2*nGridPt], 2.0, &pOrbVal[3*nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);

         // note: tau must be made in any case, even if pTau itself is not actually used.
         // the intermediate is used inside the computation of upsilon (and we need at least
         // ONE of them if we reached this place).
         double *pTau_ = pTau;
         if (!MakeTau) {
            assert(MakeUpsilon);
            Mem.Alloc(pTau_, nGridPt);
         }
         AccN(pTau_, 0.5, &pOrbVal[1*nGridPt], nMapSt, &pBxRho[1*nGridPt], nGridPt * nBx, nGridPt, nMap );
         AccN(pTau_, 0.5, &pOrbVal[2*nGridPt], nMapSt, &pBxRho[2*nGridPt], nGridPt * nBx, nGridPt, nMap, ALGOP_Add);
         AccN(pTau_, 0.5, &pOrbVal[3*nGridPt], nMapSt, &pBxRho[3*nGridPt], nGridPt * nBx, nGridPt, nMap, ALGOP_Add);

         if (MakeUpsilon) {
            assert(pTau_ != 0);
            AccN(pUpsilon, 2.0, &pOrbVal[p.iLaplaceComp0 * nGridPt], nMapSt, pBxRho, nGridPt * nBx, nGridPt, nMap);
            Add2(pUpsilon, pTau_, 4.0, nGridPt);
         }
      }
   } else if (nDiff == 2) {
      assert(0);
   } else {
      assert(0);
   }

   // make sure that density has not gone negative due to some
   // unfortunate numerical cancellation. Count electrons (for grid accuracy) BEFORE doing that.
   fElec = Dot(pRho, pGridWt, nGridPt);
   if (ForceNonNegative)
      for (size_t iPt = 0; iPt < nGridPt; ++iPt)
         if (pRho[iPt] < 0)
            pRho[iPt] = 0;
   // ^- UKS spin densities can be negative, this is why this enforcement is
   //    optional (ForceNonNegative is false in this case).
   //    There are other things which could go wrong, though (e.g., abs(spin
   //    density) > total density...)

   if (!p.Gradient)
      Mem.Free(pBeginOfStorage);
//#endif // INCLUDE_OPTIONALS
}


double *FDensitySet::CompressAuxVec(double const *pInput, FDensitySetInput const &p, FMemoryStack &Mem) const
{
   double *r;
   Mem.Alloc(r, p.nMap);
   for (size_t iMap = 0; iMap < p.nMap; ++ iMap)
      r[iMap] = pInput[p.pMap[iMap]];
   return r;
}

// evaluate densities via auxiliary expansion of density.
void FDensitySet::EvalAux(double const *pOrbVal, FDensitySetInput const &p, bool ForceNonNegative, double const *pGridWt, FMemoryStack &Mem)
{
   assert(!MakeTau);
   // ^- cannot make tau in auxiliary density expansion---it depends on the
   //    orbitals, which we do not have here.

   void
      *pBeginOfStorage = Mem.Alloc(0); // note: NOT freed in gradient case! (for keeping bxrho!)
   size_t
      nComp = p.nComp;
   size_t
      nMap = (size_t)p.nMap,
      nMapSt = nGridPt * nComp; // stride between two basis function entries in OrbVal.

   double
      *pAuxDen = CompressAuxVec(p.pRdm, p, Mem); // compress input density to Map dimensions.
   if (nDiff == 0) {
      // LDA.
      Mxv(pRho,1, pOrbVal,1,nMapSt, pAuxDen,1, nGridPt, nMap);
   } else if (nDiff == 1) {
      // GGA or meta-GGA (with upsilons). Most common case.
      Mxv(pRho,1, pOrbVal,1,nMapSt, pAuxDen,1, nGridPt, nMap);
      Mxv(&pRhoGrd[0*nGridPt],1, &pOrbVal[1*nGridPt],1,nMapSt, pAuxDen,1, nGridPt, nMap, false, 1.0);
      Mxv(&pRhoGrd[1*nGridPt],1, &pOrbVal[2*nGridPt],1,nMapSt, pAuxDen,1, nGridPt, nMap, false, 1.0);
      Mxv(&pRhoGrd[2*nGridPt],1, &pOrbVal[3*nGridPt],1,nMapSt, pAuxDen,1, nGridPt, nMap, false, 1.0);
      if (MakeUpsilon) {
         // make upsilon here, using directly evaluated laplacian of basis functions.
         assert(int(p.DerivTypeBf) >= int(DERIVCOMP_Value_Laplace));
         Mxv(pUpsilon,1, &pOrbVal[p.iLaplaceComp0*nGridPt],1,nMapSt, pAuxDen,1, nGridPt, nMap, false, 1.0);
      }
   } else if (nDiff == 2 || MakeUpsilon) {
      // make upsilon here, using cartesian second derivatives of the basis functions
      // (and rho and tau).
      // Ordering:
      //  0  1  2  3  4  5  6  7  8  9
      // rho x  y  z  xx xy xz yy yz zz
      // FIXME: delete this once not needed any more.
      Mxv(pRho, 1, pOrbVal, 1, nMapSt, pAuxDen, 1, nGridPt, nMap);
      Mxv(&pRhoGrd[0 * nGridPt], 1, &pOrbVal[1 * nGridPt], 1, nMapSt, pAuxDen, 1, nGridPt, nMap, false, 1.0);
      Mxv(&pRhoGrd[1 * nGridPt], 1, &pOrbVal[2 * nGridPt], 1, nMapSt, pAuxDen, 1, nGridPt, nMap, false, 1.0);
      Mxv(&pRhoGrd[2 * nGridPt], 1, &pOrbVal[3 * nGridPt], 1, nMapSt, pAuxDen, 1, nGridPt, nMap, false, 1.0);
      Mxv(pUpsilon, 1, &pOrbVal[4 * nGridPt], 1, nMapSt, pAuxDen, 1, nGridPt, nMap, false, 1.0);
      Mxv(pUpsilon, 1, &pOrbVal[7 * nGridPt], 1, nMapSt, pAuxDen, 1, nGridPt, nMap, true, 1.0);
      Mxv(pUpsilon, 1, &pOrbVal[9 * nGridPt], 1, nMapSt, pAuxDen, 1, nGridPt, nMap, true, 1.0);
   } else {
      assert(0);
   }

   // make sure that density has not gone negative due to some
   // unfortunate numerical cancellation. Count electrons (for grid accuracy) BEFORE doing that.
   fElec = Dot(pRho, pGridWt, nGridPt);
   if (ForceNonNegative)
      for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
         if (pRho[iPt] < 0)
            pRho[iPt] = 0;
   // ^- UKS spin densities can be negative, this is why this enforcement is
   //    optional (ForceNonNegative is false in this case).

   if ( !p.Gradient )
      Mem.Free(pBeginOfStorage);
}


// Note: the two TransformAbToCo functions work fine; they are disabled because
// we do not have the UKS frontend code in MicroScf at the moment, so they are not
// needed for the time being.
//
// void TransformAbToCo(double *pAC, double *pBO, size_t n)
// {
//    assert(pAC != 0 && pBO != 0);
//    for (size_t i = 0; i < n; ++ i) {
//       double c = pAC[i] + pBO[i];
//       double o = pAC[i] - pBO[i];
//       pAC[i] = c;
//       pBO[i] = o;
//    };
// }
//
// // transform a pair of density sets for alpha/beta densities to closed/open densities
// void TransformAbToCo(FDensitySet &AC, FDensitySet &BO, size_t nGridPt)
// {
//    if (AC.pRho)
//       TransformAbToCo(AC.pRho, BO.pRho, nGridPt);
//    if (AC.pRhoGrd)
//       TransformAbToCo(AC.pRhoGrd, BO.pRhoGrd, 3*nGridPt);
//    if (AC.pTau)
//       TransformAbToCo(AC.pTau, BO.pTau, nGridPt);
//    if (AC.pUpsilon)
//       TransformAbToCo(AC.pUpsilon, BO.pUpsilon, nGridPt);
// }


struct FDftiResultSet
{
   double
      *pZk,        // dft integrand Exc
      *pVdRho,      // [d/d rho] Exc
      *pVdSigma,    // [d/d sigmaxx] Exc, where xx = cc for closed case and oo for open case
      *pVdSigmaXY,  // [d/d sigmaxy] Exc, where xy = co (for both cases)
      *pVdTau,      // [d/d tau] Exc
      *pVdUpsilon;  // [d/d upsilon] Exc

   void Init(size_t nGridPt, uint nDiff, bool MakeTau, bool MakeUpsilon, double *pVdSigmaXY_, FMemoryStack &Mem);
};

void FDftiResultSet::Init(size_t nGridPt, uint nDiff, bool MakeTau, bool MakeUpsilon, double *pVdSigmaXY_, FMemoryStack &Mem)
{
   Mem.ClearAlloc(pZk, nGridPt);
   Mem.ClearAlloc(pVdRho, nGridPt);
   if (nDiff >= 1) {
      Mem.ClearAlloc(pVdSigma, nGridPt);
      if (pVdSigmaXY_)
         pVdSigmaXY = pVdSigmaXY_; // store link if provided by other case.
      else
         Mem.ClearAlloc(pVdSigmaXY, nGridPt);
   }
   if (nDiff >= 2 || MakeTau)
      Mem.ClearAlloc(pVdTau, nGridPt);
   if (nDiff >= 2 || MakeUpsilon)
      Mem.ClearAlloc(pVdUpsilon, nGridPt);
}

// void FDensitySet::MakeZm(double *&pZm, double *pOrbVal, double *pRhoGrd, size_t nGridPt, size_t nMap, size_t nMapSt, FMemoryStack &Mem)
// {
//    // form Zm(r) := (grad mu(r))*(grad rho(r))
//    Mem.Alloc(pZm, nGridPt * nMap);
//    for (size_t iMap = 0; iMap < nMap; ++ iMap)
//    {
//       double *pVal = &pOrbVal[nMapSt * iMap] + nGridPt; // start of gradient.
//       double *pOut = &pZm[nMap * iMap];
//       for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
//       {
//          pOut[iPt] =
//             pRhoGrd[iPt]             * pVal[iPt +   nGridPt] +
//             pRhoGrd[iPt +   nGridPt] * pVal[iPt +   nGridPt] +
//             pRhoGrd[iPt + 2*nGridPt] * pVal[iPt + 2*nGridPt];
//       };
// };


//#ifdef INCLUDE_OPTIONALS
void FDensitySet::AccXcSym1(double *pXc, double const *pVdRho, double const *pGridWt, double const *pOrbVal, size_t nMapSt, size_t nMap, uint Flags, FMemoryStack &Mem)
{
   // form the symmetric contribution
   //    k_{\mu\nu}(r) = \mu(r) vrho(r) \nu(r)
   // to the given exchange matrix. The integrands are given as input.
   TMemoryLock<double>
      pBeginOfStorage(0, &Mem);

   // check if we can make purely non-negative integrands.
   bool
      PositiveScalars = true;
   for (size_t iPt = 0; iPt < nGridPt && PositiveScalars; ++ iPt)
      if (pVdRho[iPt] * pGridWt[iPt] > 0)  // d/drho is supposed to be negative.
         PositiveScalars = false;
   if (PositiveScalars) {
      // absorb sqrt(w_g * zk_g) into the orbital values and use
      // symmetric matrix multiplication to form xc.
      double
         *pFac, *pRhoZkWt;
      Mem.Alloc(pFac, nGridPt);
      for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
         pFac[iPt] = std::sqrt(-pVdRho[iPt] * pGridWt[iPt]);
      Mem.Alloc(pRhoZkWt, nGridPt * nMap);
//          assert(nMapSt == nGridPt);
      // ^- only != in gradient case, and then we actually don't calcualte xc matrices.
      for (size_t iMap = 0; iMap < nMap; ++ iMap)
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
            pRhoZkWt[iPt + nGridPt * iMap] = pFac[iPt] * pOrbVal[iPt + nMapSt * iMap];
      // C := alpha * A*A^T + beta*C.
      // FD(dsyrk)( 'L', TransA, Out.nRows, A.nRows, fFactor, A.pData, lda, Beta, Out.pData, ldc );
      double
         fPrevXc = (0 != (Flags & ALGOP_Add))? 1.0 : 0.0;
      DSYRK('L', 'T', nMap, nGridPt,-1.0, pRhoZkWt, nGridPt, fPrevXc, pXc, nMap);
   } else {
      // copy w_g*zk_g-times the orbital values and multiply with
      // original orbital values, using standard matrix multiplication
      double
         *pVdRhoOrb, *pVdRhoWt;
      Mem.Alloc(pVdRhoWt, nGridPt);
      Mem.Alloc(pVdRhoOrb, nGridPt * nMap);

      for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
         pVdRhoWt[iPt] = pVdRho[iPt] * pGridWt[iPt];

      for (size_t iMap = 0; iMap < nMap; ++ iMap)
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
            pVdRhoOrb[iPt + nGridPt * iMap] = pOrbVal[iPt + nMapSt * iMap] * pVdRhoWt[iPt];

      if (0 != (Flags & ALGOP_Add))
         MXMB(pVdRhoOrb,nGridPt,1,  pOrbVal,1,nMapSt,  pXc,1,nMap,  nMap,nGridPt,nMap);
      else
         MXMA(pVdRhoOrb,nGridPt,1,  pOrbVal,1,nMapSt,  pXc,1,nMap,  nMap,nGridPt,nMap);
   }
   if (0 != (Flags & ALGOP_Symmetrize)) {
      // copy lower triangle to upper triangle.
      // TODO: accumulate to final matrix already here, and skip symmetrization later on.
      for (size_t iRow = 0; iRow < nMap; ++ iRow)
         for (size_t iCol = iRow; iCol < nMap; ++ iCol)
            pXc[iRow + nMap * iCol] = pXc[iCol + nMap * iRow];
   }
}
//#endif // INCLUDE_OPTIONALS


// accumulate xc matrix contribution to triangular full-dimension xc matrix at pXcTriang.
// X denotes the open/closed case at hand (e.g., 'c' when making closed-shell exchange), Y the other case.
// For example, when making the open-shell exchange matrix, then Sigma is SigmaOO and SigmaXY is SigmaCO.
void FDensitySet::AccXcMatrix(double *pXcTriang, double *pOrbVal, FDftiResultSet const &r,
      double *pRhoGrdY, FDensitySetInput const &p, double const *pGridWt, FMemoryStack &Mem)
{
   if (AuxExpandXc)
      return AccXcMatrixAux(pXcTriang, pOrbVal, r, pRhoGrdY, p, pGridWt, Mem);
//#ifdef INCLUDE_OPTIONALS

   void
      *pBeginOfStorage = Mem.Alloc(0);
   size_t
      nComp = p.nComp;
   size_t
      nMap = (size_t)p.nMap,
      nMapSt = nGridPt * nComp; // stride between two basis function entries in OrbVal.
   double
      *pXc;
   Mem.Alloc(pXc, nMap*nMap);
   if (nDiff == 0) {
      // LDA:
      //    k_{\mu\nu}(r) = \mu(r) vrho(r) \nu(r)
      // ... that's it already.
      AccXcSym1(pXc, r.pVdRho, pGridWt, pOrbVal, nMapSt, nMap, ALGOP_Symmetrize, Mem);
   } else if (nDiff == 1) {

      if (MakeTau || MakeUpsilon) {
         // GGA with tau: same as GGA without tau, but with additional term
         //    k_{\mu\nu}(r) += vtau(r) [\grad \mu(r)] [\grad \nu(r)]
         // Do that first (because of mystery symmetrization in AccXcSym1).
         //
         // GGA with upsilon... this one is a bit more complicated. Additional terms:
         // (c/p this to https://www.codecogs.com/latex/eqneditor.php for seeing the equations):
         //    \begin{align*}
         //    \upsilon(\vec r) &= \sum_{\mu\nu} \gamma^{\mu\nu} \Delta\big(\mu(\vec r) \nu(\vec r)\big)
         //    \\ &= \sum_{\mu\nu} \gamma^{\mu\nu} \left([\Delta \mu(\vec r)][\nu(\vec r)] + 2 [\vec\nabla \mu(\vec r)][\vec\nabla \nu(\vec r)] + [\mu(\vec r)][\Delta \nu(\vec r)]  \right)
         //    \\ v^\upsilon_{\mu\nu} &= \sum_g w_g\cdot\frac{\epsilon(\vec r_g)}{\upsilon(\vec r_g)}\cdot \frac{\upsilon(\vec r_g)}{\gamma^{\mu\nu}}
         //    \\ &= \sum_g w_g \frac{\epsilon(\vec r_g)}{\upsilon(\vec r_g)}\cdot\left([\Delta \mu(\vec r_g)][\nu(\vec r_g)] + 2 [\vec\nabla \mu(\vec r_g)][\vec\nabla \nu(\vec r_g)] + [\mu(\vec r_g)][\Delta \nu(\vec r_g)]\right)
         //    \\ &\stackrel{\mathrm{sym}^*}{=} \sum_g \bigg( 2 w_g  \frac{\epsilon(\vec r_g)}{\upsilon(\vec r_g)}\bigg)\cdot\left([\Delta \mu(\vec r_g)][\nu(\vec r_g)] + [\vec\nabla \mu(\vec r_g)][\vec\nabla \nu(\vec r_g)] \right)
         //    \end{align*}
         //
         // ...so basically, we get one term which looks just like the tau term (but with an
         // additional factor of 2), and one Laplace level term which looks different... but
         // hopefully will come out right after symmetrization (see comments on GGA without tau below).
         //
         // We deal with the  (2 vupsilon(r) + vtau(r)) [grad mu] [grad nu] terms here (together)
         // and absorb the laplace term into the non-symmetric terms below.
         TMemoryLock<double>
            pVdUpsilonTauLike(nGridPt, &Mem);
         if (MakeUpsilon && MakeTau) {
            for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt)
               pVdUpsilonTauLike[iGridPt] = 2.0 * r.pVdUpsilon[iGridPt] + r.pVdTau[iGridPt];
         } else if (MakeUpsilon && !MakeTau) {
            for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt)
               pVdUpsilonTauLike[iGridPt] = 2.0 * r.pVdUpsilon[iGridPt];
         } else if (!MakeUpsilon && MakeTau) {
            for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt)
               pVdUpsilonTauLike[iGridPt] = r.pVdTau[iGridPt];
         } else {
            assert(0);
         }
         AccXcSym1(pXc, pVdUpsilonTauLike, pGridWt, &pOrbVal[1*nGridPt], nMapSt, nMap, 0, Mem);
         AccXcSym1(pXc, pVdUpsilonTauLike, pGridWt, &pOrbVal[2*nGridPt], nMapSt, nMap, ALGOP_Add, Mem);
         AccXcSym1(pXc, pVdUpsilonTauLike, pGridWt, &pOrbVal[3*nGridPt], nMapSt, nMap, ALGOP_Add | ALGOP_Symmetrize, Mem);
      }

//       if ( MakeTau && !MakeUpsilon ) {
//          // GGA with tau: same as GGA without tau, but with additional term
//          //    k_{\mu\nu}(r) += vtau(r) [\grad \mu(r)] [\grad \nu(r)]
//          // Do that first (because of mystery symmetrization in AccXcSym1).
//          AccXcSym1(pXc, r.pVdTau, pGridWt, &pOrbVal[1*nGridPt], nMapSt, nMap, 0, Mem);
//          AccXcSym1(pXc, r.pVdTau, pGridWt, &pOrbVal[2*nGridPt], nMapSt, nMap, ALGOP_Add, Mem);
//          AccXcSym1(pXc, r.pVdTau, pGridWt, &pOrbVal[3*nGridPt], nMapSt, nMap, ALGOP_Add | ALGOP_Symmetrize, Mem);
//       }

      // GGA without tau:
      //    k_{\mu\nu}(r) = \mu(r) vrho(r) \nu(r) + 2 vsigma(r) (\grad rho(r))(\grad \mu(r) \nu(r))
      //                  = \mu(r) vrho(r) \nu(r) + 2 vsigma(r) (\grad rho(r))([\grad \mu(r)] \nu(r) + \mu(r) [\grad nu(r)])
      // where vsigma = d zk/d[sigma], vrho = dzk/d[rho].
      //    the second term results from:
      //          k_{\mu\nu}(r) = \mu(r) \nu(r) * d[Exc]/d[rho] + \grad(\mu(r) \nu(r))_i * d[Exc]/d[ grad rho_i ]
      //          k_{\mu\nu}(r) = \mu(r) \nu(r) * d[Exc]/d[rho] + \grad(\mu(r) \nu(r))_i * d[Exc]/d[sigma] * d[sigma]/d[grad rho_i]
      //    with sigma = d[rho]/d[x]^2 + d[rho]/d[y]^2 + d[rho]/d[z]^2, the last term gives a factor of 2 \grad rho_i (the input density).
      //    (note: the k_{\mu\nu} elements can be calculated as
      //      d/d[\gamma_{\mu\nu} \int_r exc(rho(r), sigma(r), ...)
      //     using chain rules.)
      // we now do the following:
      //    the last two terms, [\grad \mu(r)] \nu(r) + \mu(r) [\grad nu(r)], are
      //    symmetrizing combinations of each other. What we do is to evaluate only
      //    one of them, multiplied by two. The non-symmetric xc-matrix resulting
      //    from that is then explicitly symmetrized afterwards to account for the
      //    effect of the other term.
      // This way:
      //    k*_{\mu\nu}(r) = \mu(r) vrho(r) \nu(r) + 4 vsigma(r) (\grad rho(r))([\grad \mu(r)] \nu(r))
      //    k*_{\mu\nu}(r) = [\mu(r) vrho(r) + 4 vsigma(r) (\grad rho(r))[\grad \mu(r)]] \nu(r)
      // (*: needs explicit symmetrization)
      // In the open-shell case the above fomulas have an additional term
      //    + 2 vsigmaco (grad rhoo(r))  for the closed-shell exchange and
      //    + 2 vsigmaco (grad rhoc(r))  for the open-shell exchange.
      double
         *pLmu;
      Mem.Alloc(pLmu, nGridPt * nMap);
      for (size_t iMap = 0; iMap < nMap; ++ iMap)
      {
         double *pVal = &pOrbVal[nMapSt * iMap];
         double *pOut = &pLmu[nGridPt * iMap];
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
            double fDotGrd =
               pRhoGrd[iPt]             * pVal[iPt +   nGridPt] +
               pRhoGrd[iPt +   nGridPt] * pVal[iPt + 2*nGridPt] +
               pRhoGrd[iPt + 2*nGridPt] * pVal[iPt + 3*nGridPt];
            pOut[iPt] = pGridWt[iPt] * (r.pVdRho[iPt] * pVal[iPt] + 4 * r.pVdSigma[iPt] * fDotGrd);
         }
         if (MakeUpsilon) {
            // see comments on MakeUpsilon above. This is the left term.
            for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
               pOut[iPt] += pGridWt[iPt] * (2 * r.pVdUpsilon[iPt] * pVal[iPt + p.iLaplaceComp0*nGridPt]);
         }
         if (pRhoGrdY != 0)
            // open-shell contribution.
            for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
               double fDotGrd =
                  pRhoGrdY[iPt]             * pVal[iPt +   nGridPt] +
                  pRhoGrdY[iPt +   nGridPt] * pVal[iPt + 2*nGridPt] +
                  pRhoGrdY[iPt + 2*nGridPt] * pVal[iPt + 3*nGridPt];
               pOut[iPt] += 2 * pGridWt[iPt] * r.pVdSigmaXY[iPt] * fDotGrd;
            }
      }
      if (MakeTau || MakeUpsilon)
         // some stuff is already there in pXc. Add it it, don't overwrite.
         MXMB(pLmu,nGridPt,1,  pOrbVal,1,nMapSt,  pXc,1,nMap,  nMap,nGridPt,nMap);
      else
         MXMA(pLmu,nGridPt,1,  pOrbVal,1,nMapSt,  pXc,1,nMap,  nMap,nGridPt,nMap);
   } else {
      assert_rt(0);
      // fixme: add code here once other stuff works (MGGA).
   }

//       ir::PrintMatrixGen(ct::xout, pXc, nMap, 1, nMap, nMap, "XC (before symmetrization)", "{:12.5f}");
   // symmetrize & compress output xc matrix to triangular & add to previous.
   for (size_t iMapCol = 0; iMapCol < nMap; ++ iMapCol) {
      size_t
         iCol = p.pMap[iMapCol],
         iColOff = iCol*(iCol+1)/2; // offset of i'th col in output packed triangular matrix
      for (size_t iMapRow = 0; iMapRow <= iMapCol; ++ iMapRow) {
         size_t iRow = p.pMap[iMapRow];
         pXcTriang[iRow + iColOff] +=
            0.5 * (pXc[iMapRow + nMap * iMapCol] + pXc[iMapCol + nMap * iMapRow]);
      }
   }

   Mem.Free(pBeginOfStorage);
//#endif // INCLUDE_OPTIONALS
}


// accumulate xc matrix contribution to triangular full-dimension xc matrix at pXcTriang.
// X denotes the open/closed case at hand (e.g., 'c' when making closed-shell exchange), Y the other case.
// For example, when making the open-shell exchange matrix, then Sigma is SigmaOO and SigmaXY is SigmaCO.
void FDensitySet::AccXcMatrixAux(double *pXcVec, double *pOrbVal, FDftiResultSet const &r,
      double *pRhoGrdY, FDensitySetInput const &p, double const *pGridWt, FMemoryStack &Mem)
{
   assert(!MakeTau); // cannot be done in auxiliary density expansion.
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   size_t
      nComp = p.nComp;
   size_t
      nMap = (size_t)p.nMap,
      nMapSt = nGridPt * nComp; // stride between two basis function entries in OrbVal.
   double
      *pXc;
   Mem.Alloc(pXc, nMap);

   if (nDiff == 0) {
      // LDA:
      //    k_{\mu\nu}(r) = \mu(r) vrho(r) \nu(r)
      // ... that's it already.
      double *pVdRhoWt;
      Mem.Alloc(pVdRhoWt, nGridPt);
      for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
         pVdRhoWt[iPt] = r.pVdRho[iPt] * pGridWt[iPt];
      Mxv(pXc,1, pOrbVal,nMapSt,1, pVdRhoWt,1, nMap, nGridPt);
//       Mxv(pXc,1, pOrbVal,nMapSt,1, r.VdRho,1, nMap, nGridPt);
   } else if (nDiff == 1 || nDiff == 2) {
      // auxiliary version:
      //
      //    k[A] += \sum_g  w(g) vrho(g) A(g) + \sum_g w(g) D[Ecx(r)]/D[grad rho(r)] * grad A(r)
      //    k[A] += \sum_g  w(g) vrho(g) A(g) + \sum_g w(g) vsigma(g) 2 (\grad rho(r))*(\grad A(r))
      //
      // second term comes from derivative D[sigma]/D[grad rho(r)] = 2 grad rho(r)
      // since sigma = [grad rho(r)] x [grad rho(r)]
      double *pVdRhoWt, *pVdGrdWt;
      Mem.Alloc(pVdRhoWt, nGridPt);
      for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
         pVdRhoWt[iPt] = r.pVdRho[iPt] * pGridWt[iPt];
      Mem.Alloc(pVdGrdWt, 3*nGridPt);
      for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
         pVdGrdWt[iPt + 0*nGridPt] = 2. * pGridWt[iPt] * r.pVdSigma[iPt] * pRhoGrd[iPt + 0*nGridPt];
         pVdGrdWt[iPt + 1*nGridPt] = 2. * pGridWt[iPt] * r.pVdSigma[iPt] * pRhoGrd[iPt + 1*nGridPt];
         pVdGrdWt[iPt + 2*nGridPt] = 2. * pGridWt[iPt] * r.pVdSigma[iPt] * pRhoGrd[iPt + 2*nGridPt];
      }
      if (pRhoGrdY) {
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
            pVdGrdWt[iPt + 0*nGridPt] += pGridWt[iPt] * r.pVdSigmaXY[iPt] * pRhoGrdY[iPt + 0*nGridPt];
            pVdGrdWt[iPt + 1*nGridPt] += pGridWt[iPt] * r.pVdSigmaXY[iPt] * pRhoGrdY[iPt + 1*nGridPt];
            pVdGrdWt[iPt + 2*nGridPt] += pGridWt[iPt] * r.pVdSigmaXY[iPt] * pRhoGrdY[iPt + 2*nGridPt];
         }
      }
      Mxv(pXc,1, pOrbVal,nMapSt,1, pVdRhoWt,1, nMap, nGridPt);
      Mxv(pXc,1, &pOrbVal[1*nGridPt],nMapSt,1, &pVdGrdWt[0*nGridPt],1, nMap, nGridPt, true, 1.0);
      Mxv(pXc,1, &pOrbVal[2*nGridPt],nMapSt,1, &pVdGrdWt[1*nGridPt],1, nMap, nGridPt, true, 1.0);
      Mxv(pXc,1, &pOrbVal[3*nGridPt],nMapSt,1, &pVdGrdWt[2*nGridPt],1, nMap, nGridPt, true, 1.0);
      Mem.Free(pVdGrdWt);

      // meta-GGA with upsilon.
      if (MakeUpsilon) {
         // make upsilon here, using directly evaluated laplacian of basis functions.
         assert(int(p.DerivTypeBf) >= int(DERIVCOMP_Value_Laplace));
         TMemoryLock<double>
            pVdUpsilonWt(nGridPt, &Mem);
         for (size_t iPt = 0; iPt < nGridPt; ++iPt)
            pVdUpsilonWt[iPt] = r.pVdUpsilon[iPt] * pGridWt[iPt];

         Mxv(pXc, 1, &pOrbVal[p.iLaplaceComp0 * nGridPt], nMapSt, 1, pVdUpsilonWt, 1, nMap, nGridPt, true, 1.0);
      }

//       if (nDiff == 2 || MakeUpsilon) {
//          double *pVdUpsilonWt;
//          Mem.Alloc(pVdUpsilonWt, nGridPt);
//          for (size_t iPt = 0; iPt < nGridPt; ++iPt)
//             pVdUpsilonWt[iPt] = r.pVdUpsilon[iPt] * pGridWt[iPt];
//
//          Mxv(pXc, 1, &pOrbVal[4 * nGridPt], nMapSt, 1, pVdUpsilonWt, 1, nMap, nGridPt, true, 1.0);
//          Mxv(pXc, 1, &pOrbVal[7 * nGridPt], nMapSt, 1, pVdUpsilonWt, 1, nMap, nGridPt, true, 1.0);
//          Mxv(pXc, 1, &pOrbVal[9 * nGridPt], nMapSt, 1, pVdUpsilonWt, 1, nMap, nGridPt, true, 1.0);
//          Mem.Free(pVdUpsilonWt);
//       }
   } else {
      assert_rt(0);
      // fixme: add code here once other stuff works (MGGA).
   }

   // expand Map dimension & add to previous vxc vector.
   for (size_t iMapCol = 0; iMapCol < nMap; ++ iMapCol) {
      pXcVec[p.pMap[iMapCol]] += pXc[iMapCol];
   }
}


struct FDftiJob {
   FDftiArgs const
      &Args;
   FDftGrid::FGridBlock const
      *pGridBlock;
   // outputs.
   double
      *pFockC, // note: these are assumed to be thread-local (i.e., lock-free)
      *pFockO,
      *pGradient, // geometric gradient (in case iMakeGrad=1)
      *pDfuEnergies; // accumulated dft energies, one for each functional.
   size_t
      nGridPt,
      nDiff,    // derivative order of functional
      nDiffBf,  // derivative order of basis functions we need to evaluate (>= nDiff; required for Gradients)
      nComp,    // number of orbital components at pOrbVal (value; value,dx,dy,dz; value,dx,dy,dz,dxx,dxy,...)
      iLaplaceComp0, // index of first laplace derivative component in pOrbVal
      iGridBlockCenter;
      // ^- atom/center index on which the current block of grid points in locked.
      //    only used in gradient calculations, otherwise set to iNoCenter
   FEvalBfDerivType
      DerivTypeBf;
   FMemoryStack
      &Mem; // note: this one is assumed to be thread-local

   double
      *pGridWt;
   double const
      *pGridPt, *pGridWt_In;

   double
      Timings[CLOCK_TimerCount],
      fElec[2]; // integrated density/spin-density on the block.
   double
      // if in AddBackgroundDensities mode: a copy of rhoc *before* the
      // background densities were added (need that for Exc and vxc).
      *pEmbeddedRhoc;

   FDftiJob(double *pDfuEnergies, double *pFockC, double *pFockO, double *pGrad, FDftGrid::FGridBlock const *pGridBlock, FDftiArgs &Args, FMemoryStack &Mem);
   void Run();


   void EvalBfn(double *&pOrbVal, size_t *&pMap, size_t &nMap, ptrdiff_t *&pCenterIndices);
   void FormSigma(double *&pSigma, double const *pRhoGrdA, double const *pRhoGrdB);

   void EvalGridWtDeriv(double *&pWtGrd);

   // calculate and accumulate DFTI analytic gradient (core part)
   void AccCoreGradient(double *pGrad, double const *pOrbVal, ptrdiff_t const *pCenterIndices, FDftiResultSet const &r,
     FDensitySet &d, FDensitySetInput const &DenInp, double const *pRhoGrdY, double const *pSigmaXY, bool bGridGrad);
   void AccCoreGradientAux(double *pGrad, double const *pOrbVal, ptrdiff_t const *pCenterIndices, FDftiResultSet const &r,
     FDensitySet &d, FDensitySetInput const &DenInp, double const *pRhoGrdY, double const *pSigmaXY, bool bGridGrad);

   // calculate and accumulate DFTI analytic gradient (grid weight derivative part)
   void AccGridWtGradient(double *pGrad, double const *pZk);

   void StoreBackgroundDensities(FDensitySet const &DenC, FDensitySet const &DenO, FMemoryStack &Mem);
   void AddBackgroundDensities(FDensitySet &DenC, FDensitySet &DenO, FMemoryStack &Mem);
   void AdjustGridWtForEmbeddedCalc(FDensitySet const &DenC, FDensitySet const &DenO, FMemoryStack &Mem);
   void AdjustXcForEmbeddedCalc(FDftiResultSet &ResC, FDftiResultSet &ResO, FDensitySet const &DenC, FDensitySet const &DenO, FMemoryStack &Mem);
private:
   FDftiJob(FDftiJob const&); // not implemented
   void operator = (FDftiJob const&); // not implemented.
};


FEvalBfDerivType MakeBfDerivType(unsigned nDiff, bool MakeLaplaceToo)
{
   return FEvalBfDerivType(nDiff + (MakeLaplaceToo? unsigned(DERIVCOMP_Value_Laplace) : 0u));
}


FDftiJob::FDftiJob(double *pDfuEnergies_, double *pFockC_, double *pFockO_, double *pGrad_, FDftGrid::FGridBlock const *pGridBlock_, FDftiArgs &Args_, FMemoryStack &Mem_)
   : Args(Args_), pGridBlock(pGridBlock_), pFockC(pFockC_), pFockO(pFockO_), pGradient(pGrad_), pDfuEnergies(pDfuEnergies_), iGridBlockCenter(pGridBlock_->iAtomicCenter), Mem(Mem_)
{
   nGridPt = pGridBlock->nPt();
   pGridPt = &Args.pDftGrid->Positions[pGridBlock->iFirst][0];
   // make a copy of the grid weights---we might modify them under some conditions.
   Mem.Alloc(pGridWt, nGridPt);
   pGridWt_In = &Args.pDftGrid->Weights[pGridBlock->iFirst];
   assert(sizeof(pGridWt[0]) == sizeof(Args.pDftGrid->Weights[0]));
   memcpy(&pGridWt[0], &pGridWt_In[0], sizeof(pGridWt[0])*nGridPt);

   nDiff = 0;
   if (Args.pXcFn->NeedSigma())
      nDiff = 1;
//    if (Args.pXcFn->NeedUpsilon())
//       nDiff = 2; // FIXME: remove this once laplace level evaluation of basis functions is done.
   nDiffBf = nDiff;
   if (Args.Flags & DFTI_MakeGradient)
      nDiffBf += 1;
   DerivTypeBf = MakeBfDerivType(nDiffBf, Args.pXcFn->NeedUpsilon());

   // number of derivative components.
   nComp = GetNumBfDerivComps(DerivTypeBf, &iLaplaceComp0);

   for (uint i = 0; i < CLOCK_TimerCount; ++ i)
      Timings[i] = 0;
   fElec[0] = 0;
   fElec[1] = 0;

   pEmbeddedRhoc = 0;
}


//#ifdef INCLUDE_OPTIONALS
// calculate and accumulate DFTI analytic gradient
void FDftiJob::AccCoreGradient(double *pGrad, double const *pOrbVal, ptrdiff_t const *pCenterIndices, FDftiResultSet const &r,
     FDensitySet &d, FDensitySetInput const &DenInp, double const *pRhoGrdY, double const *pVdSigmaXY, bool bGridGrad)
{
   if (Args.UseAuxiliaryExpansion())
      return AccCoreGradientAux(pGrad, pOrbVal, pCenterIndices, r, d, DenInp, pRhoGrdY, pVdSigmaXY, bGridGrad);

   size_t
      nMap = DenInp.nMap,
      nMapSt = nGridPt * nComp;
   if (nDiff == 0) {
      // LDA. This term is dzk/drho * drho/dx_i.
      //      The latter is expanded as 2 * [\mu(r) gamma_{\mu\nu}] * d[\nu]/d[x_i]
      assert(d.pBxRho != 0 && d.nBx == 1); // make in FDenstiySet::Eval.
      double
         *pVgBxRho = d.pBxRho; // will be changed!

      assert(iGridBlockCenter != iNoCenter);
      assert(nComp == 4 && nDiffBf == 1);
      for (size_t iMap = 0; iMap < nMap; ++ iMap){
         size_t
            iCen = (size_t)pCenterIndices[iMap];
         if (iCen == iGridBlockCenter && bGridGrad)
            continue;
         double const
            *pGx = &pOrbVal[iMap * nMapSt +   nGridPt],
            *pGy = &pOrbVal[iMap * nMapSt + 2*nGridPt],
            *pGz = &pOrbVal[iMap * nMapSt + 3*nGridPt];
         double
            *pVgBx = &pVgBxRho[iMap * nGridPt];
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
            pVgBx[iPt] *= 2.0 * r.pVdRho[iPt] * pGridWt[iPt];
         double
            dx = Dot2(pGx, pVgBx, nGridPt),
            dy = Dot2(pGy, pVgBx, nGridPt),
            dz = Dot2(pGz, pVgBx, nGridPt);

         pGrad[3*iCen+0] -= dx;
         pGrad[3*iCen+1] -= dy;
         pGrad[3*iCen+2] -= dz;

         if (bGridGrad) {
            pGrad[3*iGridBlockCenter+0] += dx;
            pGrad[3*iGridBlockCenter+1] += dy;
            pGrad[3*iGridBlockCenter+2] += dz;
         }
      }

//       ir::PrintMatrixGen(ct::xout, pGrad, 3, 1, nCenters, 3, "GRD: (2 * vrho[g] * bxrho[g]) * D rho[g]", "{:14.8f}");
   } else if (nDiff == 1) {
      // order of basis function derivative components in a full set of
      // derivative components of DerivOrder == 2, including laplace components
      // and their derivatives (if made)
      enum {
         dx = 1,  dy = 2, dz = 3, dxx = 4, dxy = 5, dxz = 6, dyy = 7, dyz = 8, dzz = 9,
         laplace = 10, laplace_dx = 11, laplace_dy = 12, laplace_dz = 13
      };

      // contract bf and dbf with with density matrix
      // This forms D[\mu,\nu] bf[g,\nu] and D[\mu,\nu] [\del_i bf](g,\nu).
      // (update: now already made in FDenstiySet::Eval)
      assert(d.pBxRho != 0 && d.nBx == 4);
      double
         *pBxRho = d.pBxRho;

      // iterate over mu.
      assert(iGridBlockCenter != iNoCenter);
      for (size_t iMap = 0; iMap < nMap; ++ iMap){
         uint
            iCen = (uint)pCenterIndices[iMap];
         if (iCen == iGridBlockCenter && bGridGrad)
            continue;
         double
            dxyz[3] = {0},
            *pY,
            *pBx = &pBxRho[nGridPt*4 * iMap];
         double const
            *pOv = &pOrbVal[nMapSt * iMap];
         Mem.Alloc(pY, nGridPt * 4);

         for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
            double
               w2vrhoc  = 2.0 * pGridWt[iPt] * r.pVdRho[iPt],
               w4vsigma = 4.0 * pGridWt[iPt] * r.pVdSigma[iPt],
               *pY0 = &pY[iPt],
               *pBx0 = &pBx[iPt],
               *pRhoGrd0 = &d.pRhoGrd[iPt];
            pY0[0*nGridPt] = w2vrhoc * pBx0[0] +
               w4vsigma * (pBx0[1*nGridPt] * pRhoGrd0[0*nGridPt] +
                           pBx0[2*nGridPt] * pRhoGrd0[1*nGridPt] +
                           pBx0[3*nGridPt] * pRhoGrd0[2*nGridPt]);
            pY0[1*nGridPt] = w4vsigma * pBx0[0] * pRhoGrd0[0*nGridPt];
            pY0[2*nGridPt] = w4vsigma * pBx0[0] * pRhoGrd0[1*nGridPt];
            pY0[3*nGridPt] = w4vsigma * pBx0[0] * pRhoGrd0[2*nGridPt];
            if ( d.MakeTau ) {
               // add tau contribution
               //  \sum_j 2*(d\mu(r)]/[dx_j] \gamma_{\mu\nu}]*d/d[x_i] d[\nu(r)]/[dx_j])
               double
                  w2vtauc  = 2.0 * pGridWt[iPt] * r.pVdTau[iPt];
               assert(d.nBx >= 4);
               pY0[1*nGridPt] += w2vtauc * pBx0[1*nGridPt];
               pY0[2*nGridPt] += w2vtauc * pBx0[2*nGridPt];
               pY0[3*nGridPt] += w2vtauc * pBx0[3*nGridPt];
            }
            if (d.MakeUpsilon) {
               // add upsilon contribution:
               //  (1)  (w_g \frac{\partial \epsilon}{\partial \upsilon}) * \frac{\partial\upsilon}{\partial x_i}
               //             |--- VdUpsilon --- ..                ..---|   |--- see next ---                ---|
               //  hm... maybe shorter:
               //  (2') \upsilon(\vec r) = \Delta \rho(\vec r)
               //                        = \sum_{\mu\nu} \gamma^{\mu\nu} \Delta\left(\mu(\vec r)]\nu(\vec r)\right)
               //                       (= \sum_{\mu\nu} \gamma^{\mu\nu} \left([\Delta \mu(\vec r)]\nu(\vec r) + 2 [\vec\nabla \mu(\vec r)][\vec\nabla \nu(\vec r)] + \mu(\vec r)[\Delta \nu(\vec r)]  \right))
               //  (3') \frac{\partial \upsilon}{\partial x_i} = \sum_{\mu\nu} \gamma^{\mu\nu} \Delta\left(2 \mu(\vec r)\partial_i \nu(\vec r)\right)
               //           = 2 \sum_{\mu\nu} \gamma^{\mu\nu} \left([\Delta \mu(\vec r)][\partial_i\nu(\vec r)] + 2 [\vec\nabla \mu(\vec r)][\partial_i\vec\nabla \nu(\vec r)] + \mu(\vec r)[\partial_i\Delta \nu(\vec r)]  \right)
               //       (i.e., use symmetry of gamma here to get only i-derivatives on \nu, not on \mu, before starting to expand the laplace operator) -> less terms & simpler derivation.
               //  [ (3) \frac{\partial \upsilon}{\partial x_i} = \sum_{\mu\nu} \gamma^{\mu\nu} \left([\partial_i\Delta \mu(\vec r)]\nu(\vec r) + [\Delta \mu(\vec r)][\partial_i\nu(\vec r)] + 2 [\partial_i\vec\nabla \mu(\vec r)][\vec\nabla \nu(\vec r)] + 2 [\vec\nabla \mu(\vec r)][\partial_i\vec\nabla \nu(\vec r)] + [\partial_i\mu(\vec r)][\Delta \nu(\vec r)] + \mu(\vec r)[\partial_i\Delta \nu(\vec r)]  \right) ]

               //    \begin{align*}
               //    \frac{\partial \upsilon}{\partial x_i} &= \partial_i \sum_{\mu\nu} \gamma^{\mu\nu} \Delta\big(\mu(\vec r)\nu(\vec r)\big)
               //    \\ &= \sum_{\mu\nu} \gamma^{\mu\nu} \Delta\big(2 \mu(\vec r)\partial_i \nu(\vec r)\big)
               //    \\ &= 2 \sum_{\mu\nu} \gamma^{\mu\nu} \left([\Delta \mu(\vec r)][\partial_i\nu(\vec r)] + 2 [\vec\nabla \mu(\vec r)][\partial_i\vec\nabla \nu(\vec r)] + \mu(\vec r)[\partial_i\Delta \nu(\vec r)]  \right)
               //    \\ &= 2 \sum_{\mu\nu} \gamma^{\mu\nu} \left([\partial_i\mu(\vec r)][\Delta \nu(\vec r)] + 2 [\vec\nabla \mu(\vec r)][\partial_i\vec\nabla \nu(\vec r)] + \mu(\vec r)[\partial_i\Delta \nu(\vec r)]  \right)
               //    \end{align*}
               //  (see: https://www.codecogs.com/latex/eqneditor.php for checking online)
               double
                  w2vupsilonc = 2.0 * pGridWt[iPt] * r.pVdUpsilon[iPt],
                  w4vupsilonc = 2.0 * w2vupsilonc;

               // middle term of (3'):
               //    2 [\sum_{\mu} \gamma^{\mu\nu} \vec\nabla \mu(\vec r)][\partial_i\vec\nabla \nu(\vec r)]
               // -> has same form as tau term (but is larger by a factor of 2 (or 4 of tau is defined with
               //    the 1/2 factor))
               assert(d.nBx >= 4);
               pY0[dx*nGridPt] += w4vupsilonc * pBx0[dx*nGridPt];
               pY0[dy*nGridPt] += w4vupsilonc * pBx0[dy*nGridPt];
               pY0[dz*nGridPt] += w4vupsilonc * pBx0[dz*nGridPt];

               double const
                  *pOv0 = &pOv[iPt];
               // left term of (3'):
               // [\sum_{\mu} \gamma^{\mu\nu} \partial_i\mu(\vec r)] [\Delta \nu(\vec r)]
               double const
                  lw2vupsilonc = w2vupsilonc * pOv0[laplace*nGridPt];
               dxyz[0] += pBx0[dx*nGridPt] * lw2vupsilonc;
               dxyz[1] += pBx0[dy*nGridPt] * lw2vupsilonc;
               dxyz[2] += pBx0[dz*nGridPt] * lw2vupsilonc;

               // right term of (3'):
               // [\sum_{\mu} \gamma^{\mu\nu} \mu(\vec r)] [\partial_i\Delta \nu(\vec r)]
               double const
                  bx00w2vupsilonc = pBx0[0] * w2vupsilonc;
               dxyz[0] += bx00w2vupsilonc * pOv0[laplace_dx*nGridPt];
               dxyz[1] += bx00w2vupsilonc * pOv0[laplace_dy*nGridPt];
               dxyz[2] += bx00w2vupsilonc * pOv0[laplace_dz*nGridPt];
//                throw std::runtime_error("FDftiJob::AccCoreGradient: Upsilon terms (dfxc=0). Something seems fishy here (cf dfxc_vs_ll_vs_exact/tpss_acrylic_acid.sh &> dfxc_vs_ll_vs_exact/tpss_acrylic_acid.txt). Todo: maybe try non-self-consistent densities? Eval densities with LDA and then gradients with the various TPSS variants?");
               assert(!"FDftiJob::AccCoreGradient: Upsilon terms (dfxc=0). Something seems fishy here (cf dfxc_vs_ll_vs_exact/tpss_acrylic_acid.sh &> dfxc_vs_ll_vs_exact/tpss_acrylic_acid.txt). Todo: maybe try non-self-consistent densities? Eval densities with LDA and then gradients with the various TPSS variants?");
            }
            if (Args.OpenShell()) {
               double const
                  w2vsigmaXY = 2.0 * pGridWt[iPt] * pVdSigmaXY[iPt],
                  *pRhoGrdY0 = &pRhoGrdY[iPt];
               pY0[0*nGridPt] +=
                  w2vsigmaXY * (pBx0[1*nGridPt] * pRhoGrdY0[0*nGridPt] +
                                pBx0[2*nGridPt] * pRhoGrdY0[1*nGridPt] +
                                pBx0[3*nGridPt] * pRhoGrdY0[2*nGridPt]);
               pY0[1*nGridPt] += w2vsigmaXY * pBx0[0] * pRhoGrdY0[0*nGridPt];
               pY0[2*nGridPt] += w2vsigmaXY * pBx0[0] * pRhoGrdY0[1*nGridPt];
               pY0[3*nGridPt] += w2vsigmaXY * pBx0[0] * pRhoGrdY0[2*nGridPt];
            }
         }

         // hm... what was the point of these dot products again? Doesn't
         // this just do exactly the same thing as what I would get if I'd
         // just add to dxyz[...] directly in the iPT loop while making the
         // Y contributions here? (as in the new MakeUpsilon case)
         dxyz[0] += Dot2(&pY[ 0 * nGridPt], &pOv[dx  * nGridPt], nGridPt);
         dxyz[1] += Dot2(&pY[ 0 * nGridPt], &pOv[dy  * nGridPt], nGridPt);
         dxyz[2] += Dot2(&pY[ 0 * nGridPt], &pOv[dz  * nGridPt], nGridPt);
         dxyz[0] += Dot2(&pY[dx * nGridPt], &pOv[dxx * nGridPt], nGridPt);
         dxyz[0] += Dot2(&pY[dy * nGridPt], &pOv[dxy * nGridPt], nGridPt);
         dxyz[1] += Dot2(&pY[dx * nGridPt], &pOv[dxy * nGridPt], nGridPt);
         dxyz[0] += Dot2(&pY[dz * nGridPt], &pOv[dxz * nGridPt], nGridPt);
         dxyz[1] += Dot2(&pY[dy * nGridPt], &pOv[dyy * nGridPt], nGridPt);
         dxyz[1] += Dot2(&pY[dz * nGridPt], &pOv[dyz * nGridPt], nGridPt);
         dxyz[2] += Dot2(&pY[dy * nGridPt], &pOv[dyz * nGridPt], nGridPt);
         dxyz[2] += Dot2(&pY[dx * nGridPt], &pOv[dxz * nGridPt], nGridPt);
         dxyz[2] += Dot2(&pY[dz * nGridPt], &pOv[dzz * nGridPt], nGridPt);

         pGrad[3*iCen+0] -= dxyz[0];
         pGrad[3*iCen+1] -= dxyz[1];
         pGrad[3*iCen+2] -= dxyz[2];

         if (bGridGrad) {
            pGrad[3*iGridBlockCenter+0] += dxyz[0];
            pGrad[3*iGridBlockCenter+1] += dxyz[1];
            pGrad[3*iGridBlockCenter+2] += dxyz[2];
         }

         Mem.Free(pY);
      }
//       ir::PrintMatrixGen(ct::xout, pGrad, 3, 1, nCenters, 3, "GRD: GGA", "{:14.8f}");
   } else {
      // no code.
      assert_rt(0);
   }
}

void FDftiJob::AccCoreGradientAux(double *pGrad, double const *pOrbVal, ptrdiff_t const *pCenterIndices, FDftiResultSet const &r,
     FDensitySet &d, FDensitySetInput const &DenInp, double const *pRhoGrdY, double const *pVdSigmaXY, bool bGridGrad)
{
   assert(!Args.pXcFn->NeedTau()); // cannot be done in auxiliary density expansion.
//    if (pRhoGrdY != 0) throw std::runtime_error("open-shell not yet supported in AccCoreGradientAux.");
//    IR_SUPPRESS_UNUSED_WARNING(pVdSigmaXY);

//    void
//       *pBeginOfStorage = Mem.Alloc(0); // note: NOT freed in gradient case! (for keeping bxrho!)
   size_t
      nMap = (size_t)DenInp.nMap,
      nMapSt = nGridPt * nComp; // stride between two basis function entries in OrbVal.

   double
      *pAuxDen = d.CompressAuxVec(DenInp.pRdm, DenInp, Mem); // compress input density to Map dimensions.

   if (nDiff == 0) {
      // LDA. This term is dzk/drho * drho/dx_i.
      // w/ rho(r) = \sum[A] A(r) pAuxDen[A]  and assuming that pAuxDen is
      // stationary wrt changes in x
      assert(iGridBlockCenter != iNoCenter);
      assert(nComp == 4 && nDiffBf == 1);
      for (size_t iMap = 0; iMap < nMap; ++ iMap){
         size_t
            iCen = (size_t)pCenterIndices[iMap];
         if (iCen == iGridBlockCenter && bGridGrad)
            continue;
         double const
            *pGx = &pOrbVal[iMap * nMapSt +   nGridPt],
            *pGy = &pOrbVal[iMap * nMapSt + 2*nGridPt],
            *pGz = &pOrbVal[iMap * nMapSt + 3*nGridPt];
         double
            dx = 0., dy = 0., dz = 0.;
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
            double dzk = r.pVdRho[iPt] * pGridWt[iPt];
            dzk *= pAuxDen[iMap];
            dx += dzk * pGx[iPt];
            dy += dzk * pGy[iPt];
            dz += dzk * pGz[iPt];
         }

         pGrad[3*iCen+0] -= dx;
         pGrad[3*iCen+1] -= dy;
         pGrad[3*iCen+2] -= dz;

         if (bGridGrad) {
            pGrad[3*iGridBlockCenter+0] += dx;
            pGrad[3*iGridBlockCenter+1] += dy;
            pGrad[3*iGridBlockCenter+2] += dz;
         }
      }
   } else if (nDiff == 1) {
      // we have
      //    sigma = [grad rho] * [grad rho].
      // In the auxiliary expansion case, rho is simply
      //    rho(r) = \sum_g rho[A] * A(r)
      // so
      //    grad rho(r) = \sum_A rho[A] * [grad A](r).
      //
      // I need to compute
      //    dzk/dsigma * dsigma/dx_i.
      // This is
      //    dsigma/dx_i = d/dx_i ([grad rho] * [grad rho])
      //                = 2 [grad rho] * d[grad rho]/d x_i
      // the latter simply being a second basis function derivative. so... not very complicated.

      // Code probably needs some cleanup/optimization... need to check how long it takes.
      assert(iGridBlockCenter != iNoCenter);
      bool
         HaveUpsilon = Args.pXcFn->NeedUpsilon();
      assert((!HaveUpsilon && nComp == 10 && nDiffBf == 2 && DerivTypeBf == DERIVCOMP_Value_dFirst_dSecond) ||
              (HaveUpsilon && nComp == 14 && nDiffBf == 2 && DerivTypeBf == DERIVCOMP_Value_dFirst_dSecond_Laplace_dFirst));
      for (size_t iMap = 0; iMap < nMap; ++ iMap){
         size_t
            iCen = (size_t)pCenterIndices[iMap];
         if (iCen == iGridBlockCenter && bGridGrad)
            continue;
         double
            gx = 0., gy = 0., gz = 0.;
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt) {
            double const
               *ov = &pOrbVal[iMap * nMapSt + iPt];
            enum {
               dx = 1,  dy = 2, dz = 3, dxx = 4, dxy = 5, dxz = 6, dyy = 7, dyz = 8, dzz = 9,
               laplace = 10, laplace_dx = 11, laplace_dy = 12, laplace_dz = 13
            };
            double
               dzk = r.pVdRho[iPt] * pGridWt[iPt] * pAuxDen[iMap];
            gx += dzk * ov[dx*nGridPt];
            gy += dzk * ov[dy*nGridPt];
            gz += dzk * ov[dz*nGridPt];

            if (HaveUpsilon) {
//                double const
//                   *ov_laplace = &ov[iLaplaceComp0 * nGridPt];
//                double
//                   dupsilon = r.pVdUpsilon[iPt] * pGridWt[iPt] * pAuxDen[iMap];
//                gx += dupsilon * ov_laplace[dx*nGridPt];
//                gy += dupsilon * ov_laplace[dy*nGridPt];
//                gz += dupsilon * ov_laplace[dz*nGridPt];
               double
                  dupsilon = r.pVdUpsilon[iPt] * pGridWt[iPt] * pAuxDen[iMap];
               gx += dupsilon * ov[laplace_dx*nGridPt];
               gy += dupsilon * ov[laplace_dy*nGridPt];
               gz += dupsilon * ov[laplace_dz*nGridPt];
            }

            {
               double
                  dsigma = 2. * r.pVdSigma[iPt] * pGridWt[iPt] * pAuxDen[iMap];
               double
                  // gradient of the density at current point. They are actually independent of iMap.
                  RhoGrdX = d.pRhoGrd[            iPt],
                  RhoGrdY = d.pRhoGrd[1*nGridPt + iPt],
                  RhoGrdZ = d.pRhoGrd[2*nGridPt + iPt];
               gx += dsigma * (RhoGrdX * ov[dxx*nGridPt] + RhoGrdY * ov[dxy*nGridPt] + RhoGrdZ * ov[dxz*nGridPt]);
               gy += dsigma * (RhoGrdX * ov[dxy*nGridPt] + RhoGrdY * ov[dyy*nGridPt] + RhoGrdZ * ov[dyz*nGridPt]);
               gz += dsigma * (RhoGrdX * ov[dxz*nGridPt] + RhoGrdY * ov[dyz*nGridPt] + RhoGrdZ * ov[dzz*nGridPt]);
            }
            if (pRhoGrdY) {
               // sigma_co open-shell/closed-shell mixed contributions.
               assert(r.pVdSigmaXY == pVdSigmaXY);
               double
                  dsigmaXY = pVdSigmaXY[iPt] * pGridWt[iPt] * pAuxDen[iMap];
               double
                  RhoGrdX = pRhoGrdY[            iPt],
                  RhoGrdY = pRhoGrdY[1*nGridPt + iPt],
                  RhoGrdZ = pRhoGrdY[2*nGridPt + iPt];
               gx += dsigmaXY * (RhoGrdX * ov[dxx*nGridPt] + RhoGrdY * ov[dxy*nGridPt] + RhoGrdZ * ov[dxz*nGridPt]);
               gy += dsigmaXY * (RhoGrdX * ov[dxy*nGridPt] + RhoGrdY * ov[dyy*nGridPt] + RhoGrdZ * ov[dyz*nGridPt]);
               gz += dsigmaXY * (RhoGrdX * ov[dxz*nGridPt] + RhoGrdY * ov[dyz*nGridPt] + RhoGrdZ * ov[dzz*nGridPt]);
            }
         }

         pGrad[3*iCen+0] -= gx;
         pGrad[3*iCen+1] -= gy;
         pGrad[3*iCen+2] -= gz;

         if (bGridGrad) {
            pGrad[3*iGridBlockCenter+0] += gx;
            pGrad[3*iGridBlockCenter+1] += gy;
            pGrad[3*iGridBlockCenter+2] += gz;
         }
      }
   } else {
      // no code.
      assert_rt(0);
   }
}





// note: allocates pSigma on this->Mem.
void FDftiJob::FormSigma(double *&pSigma, double const *pRhoGrdA, double const *pRhoGrdB)
{
   // sigmaAB = [grad rhoA]*[grad rhoB]
   Mem.Alloc(pSigma, nGridPt);
   for ( size_t iPt = 0; iPt < nGridPt; ++ iPt ){
      pSigma[iPt] = pRhoGrdA[iPt            ] * pRhoGrdB[iPt            ] +
                    pRhoGrdA[iPt +   nGridPt] * pRhoGrdB[iPt +   nGridPt] +
                    pRhoGrdA[iPt + 2*nGridPt] * pRhoGrdB[iPt + 2*nGridPt];
   }
}


static void StoreDensityComponentsAtPoint(double *pOut, FDensitySet const &Den, size_t iGridPt)
{
   size_t
      nGridPt = Den.nGridPt;
   if (Den.pRho)
      pOut[0] = Den.pRho[iGridPt];
   if (Den.pRhoGrd) {
      pOut[1] = Den.pRhoGrd[0*nGridPt + iGridPt]; // d/dx
      pOut[2] = Den.pRhoGrd[1*nGridPt + iGridPt]; // d/dy
      pOut[3] = Den.pRhoGrd[2*nGridPt + iGridPt]; // d/dz
   }
   if (Den.pTau)
      pOut[4] = Den.pTau[iGridPt];
   if (Den.pUpsilon)
      pOut[5] = Den.pUpsilon[iGridPt];
}


static void AddDensityComponentsAtPoint(FDensitySet &Den, double const *pIn, size_t iGridPt)
{
   size_t
      nGridPt = Den.nGridPt;
   if (Den.pRho)
      Den.pRho[iGridPt] += pIn[0];
   if (Den.pRhoGrd) {
      Den.pRhoGrd[0*nGridPt + iGridPt] += pIn[1]; // d/dx
      Den.pRhoGrd[1*nGridPt + iGridPt] += pIn[2]; // d/dy
      Den.pRhoGrd[2*nGridPt + iGridPt] += pIn[3]; // d/dz
   }
   if (Den.pTau)
      Den.pTau[iGridPt] += pIn[4];
   if (Den.pUpsilon)
      Den.pUpsilon[iGridPt] += pIn[5];
}


void FDftiJob::StoreBackgroundDensities(FDensitySet const &DenC, FDensitySet const &DenO, FMemoryStack &Mem)
{
   // we store data interleaved (density components as fast index). Makes addressing easier.
   // No need for special handling of OpenMP, since each thread writes to a different
   // target location (since each thread covers a different subset of the grid)
   double
      *pBgDens = &Args.pBackgroundDensities[DFTI_nBackgroundDensityComponents * pGridBlock->iFirst];
   memset(pBgDens, 0, sizeof(*pBgDens) * (DFTI_nBackgroundDensityComponents * nGridPt));
   // ^- should not really be necessary.

   for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt) {
      double
         *pOut = &pBgDens[DFTI_nBackgroundDensityComponents * iGridPt];
      StoreDensityComponentsAtPoint(pOut, DenC, iGridPt);
      StoreDensityComponentsAtPoint(pOut + (DFTI_nBackgroundDensityComponents/2), DenO, iGridPt);
   }
   IR_SUPPRESS_UNUSED_WARNING(DenO);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FDftiJob::AddBackgroundDensities(FDensitySet &DenC, FDensitySet &DenO, FMemoryStack &Mem)
{
   // make a copy of the density of the current input orbitals/density matrix,
   // before adding the background densities.
   Mem.Alloc(pEmbeddedRhoc, nGridPt);
   // ^- FIXME: figure out open-shell case. For now just closed-shell.

   double const
      *pBgDens = &Args.pBackgroundDensities[DFTI_nBackgroundDensityComponents * pGridBlock->iFirst];

   for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt) {
      pEmbeddedRhoc[iGridPt] = DenC.pRho[iGridPt];
      double const
         *pIn = &pBgDens[DFTI_nBackgroundDensityComponents * iGridPt];
      AddDensityComponentsAtPoint(DenC, pIn, iGridPt);
      AddDensityComponentsAtPoint(DenO, pIn + (DFTI_nBackgroundDensityComponents/2), iGridPt);
   }
}


void FDftiJob::AdjustGridWtForEmbeddedCalc(FDensitySet const &DenC, FDensitySet const &DenO, FMemoryStack &Mem)
{
   // absorb (rho(emb)/rho(tot)) factors from derivative of embedded xc energy
   // into grid weights. This takes care of the density functional energies
   // and most of the Vd** derivative terms (all except VdRho, which needs some
   // additional hacking).
   for (size_t iGridPt = 0; iGridPt != nGridPt; ++iGridPt) {
      if (std::abs(DenC.pRho[iGridPt]) > 1e-20)
         pGridWt[iGridPt] *= pEmbeddedRhoc[iGridPt] / DenC.pRho[iGridPt];
      // (not that technically the embedded rhoc cannot be larger than DenC.pRho,
      // since the latter is the sum of the embedded rhoc and the background densities,
      // neither of which can be negative... in principle at least)
   }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   IR_SUPPRESS_UNUSED_WARNING(DenO);
}


void FDftiJob::AdjustXcForEmbeddedCalc(FDftiResultSet &ResC, FDftiResultSet &ResO,
      FDensitySet const &DenC, FDensitySet const &DenO, FMemoryStack &Mem)
{
   if (Args.OpenShell())
      throw std::runtime_error("AdjustXcForEmbeddedCalc: not implemented for open-shell case.");
   for (size_t iGridPt = 0; iGridPt != nGridPt; ++iGridPt) {
      double
         RhoTotal = DenC.pRho[iGridPt],
         RhoEmb = pEmbeddedRhoc[iGridPt];
      if (std::abs(RhoTotal) > 1e-20)
         ResC.pVdRho[iGridPt] += ResC.pZk[iGridPt] * (1./RhoEmb - 1./RhoTotal);
      // ^- ... not sure about zero embedded densities. That is not exactly unlikely to happen
      //    (stuff is computed on the full grid, even for small embedded regions).
      //    This should get cancelled by the corresponding RhoEmb factor absorbed into the
      //    weights, but is this good?
      //    Might need to be handled in a different way.
      //    (e.g., by patching up the weights only after the DFTFUN calculation, and
      //    doing something like
      //         pDfuEnergies -= Dot(zk,wt)       // counter contribution from DFTFUN
      //    (patch up weights)
      //         pDfuEnergies += Dot(zk,wt-patched-up)
      //    hm... no. I can't patch up the weights if I do this. Would probably just need
      //    to go again and manually scale all Zk and Vd contributions...
   }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   IR_SUPPRESS_UNUSED_WARNING(DenO);
   IR_SUPPRESS_UNUSED_WARNING(ResO);
}




// entry points for DFT integration for one grid block.
void FDftiJob::Run()
{
   TMemoryLock<char>
      pBaseOfStorage(0, &Mem);

   // evaluate the basis functions (+derivatives) on the grid.
   RESUME_CPU_CLOCK(CLOCK_EvalBfn);
   double
      *pOrbVal;
   size_t
      nMap, *pMap;
   ptrdiff_t
      *pCenterIndices;
   EvalBfn(pOrbVal, pMap, nMap, pCenterIndices);
   PAUSE_CPU_CLOCK(CLOCK_EvalBfn);
   if (nMap == 0)
      return; // no functions left after screening.

   bool
      UseOrbitals = false;
//#ifdef INCLUDE_OPTIONALS
   if (Args.pOccOrbC && !Args.UseAuxiliaryExpansion()) {
      // decide on which density factorization to use (orbital or density).
      size_t
         nCostDen = nMap * nMap,
         nCostOrb = nMap * (Args.nOccC + Args.nOccO);
      if (Args.OpenShell()) nCostDen *= 2; // same cost for open-shell, even if only one active orbital.
      if (nDiff == 1 && Args.NeedTau()) nCostDen *= 4; // need to transform drho additionally.
      if (nDiff >= 2 ) nCostDen *= 5; // rho, drho, and laplace rho.
      if (nDiff == 1 && !Args.NeedTau()) nCostOrb *= 2; // need to transform rho; requires two multiplications (nmap->occ, occ->nmap).
      if (nDiff == 1 && Args.NeedTau()) nCostOrb *= 4; // need to transform rho and drho (only nmap->occ).
      if (nDiff == 2) nCostOrb *= 5; // need to transform rho, drho, and laplace rho (only nmap->occ).
      UseOrbitals = (nCostOrb <= nCostDen);
   }
   // ^- note: this is mainly here because in the HaveAB case we really
   //    need to do both A and B with the same factorization, because we
   //    actually have C/O *density matrices* (but not orbitals!) in any
   //    case. These cannot be mixed with the A/B ones. In the non-HaveAB-
   //    case it would be better to decide dynamically in FDensitySet::Eval
   //    which factorization to use. Should be fixed at some point in time.
//#endif // INCLUDE_OPTIONALS

   // make electron densities (+derivatives).
   RESUME_CPU_CLOCK(CLOCK_FormRho);
   FDensitySet
      DenC = FDensitySet(), DenO = FDensitySet();
   FDensitySetInput
      DenInpC(Args.pDenC, Args.pOccOrbC, Args.nBf(), Args.nOccC,
         pMap, nMap, UseOrbitals, Args.MakeGradient(), nComp, iLaplaceComp0, DerivTypeBf),
      DenInpO;

   bool
      ForceNonNegative = true;
//    if ( Args.GridDensityFlags == 2 )
//       // just evaluate densities -- for some applications these might not be
//       // entirely positive.
//       ForceNonNegative = false;
   if (Args.UseAuxiliaryExpansion())
      ForceNonNegative = false; // see below.

   DenC.Init(nGridPt, nDiff, Args.NeedTau(), Args.NeedUpsilon(), Args.UseAuxiliaryExpansion(), Mem);
   DenC.Eval(pOrbVal, DenInpC, ForceNonNegative, pGridWt, Mem);

   if (Args.OpenShell()) {
      DenInpO = DenInpC;
      DenInpO.pRdm = Args.pDenO;
      DenInpO.pOccOrb = Args.pOccOrbO;
      DenInpO.nOcc = Args.nOccO;

      DenO.Init(nGridPt, nDiff, Args.NeedTau(), Args.NeedUpsilon(), Args.UseAuxiliaryExpansion(), Mem);
      DenO.Eval(pOrbVal, DenInpO, false, pGridWt, Mem);

//       bool
//          bHaveAB = (bool)Args.iUseAB && UseOrbitals;
//       if (bHaveAB)
//          // input was alpha/beta orbitals, thus we have alpha/beta densities now.
//          // Transform this into closed/open densities, which the rest of the
//          // code expects.
//          TransformAbToCo(DenC, DenO, nGridPt);
   }
   if (Args.UseAuxiliaryExpansion()) {
      // flip weights and densities of points with negative density. reduced gradients (sigma) and taus
      // should be unaffected by this, since they are bilinear in the density. Note that we made a copy
      // of the grid weights for this purpose.
      // WARNING: assumes closed/open input, not A/B input!
      for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt)
         if (DenC.pRho[iGridPt] < 0) {
            DenC.pRho[iGridPt] *= -1.;
            pGridWt[iGridPt] *= -1.;
            if (DenO.pRho)
               DenO.pRho[iGridPt] *= -1.;
            if (DenC.pUpsilon)
               DenC.pUpsilon[iGridPt] *= -1.;
            if (DenO.pUpsilon)
               DenO.pUpsilon[iGridPt] *= -1.;
         }
   }
   PAUSE_CPU_CLOCK(CLOCK_FormRho);

   // if requested, store densities to external array or add densities from
   // external array.
   if (Args.MakeBackgroundDensities())
      StoreBackgroundDensities(DenC, DenO, Mem);
   if (Args.AddBackgroundDensities()) {
      AddBackgroundDensities(DenC, DenO, Mem);
      AdjustGridWtForEmbeddedCalc(DenC, DenO, Mem);
   }
   // form intermediates sigma_cc/sigma_co/sigma_oo.
   double
      *pSigmaCC = 0, *pSigmaOO = 0, *pSigmaCO = 0;
   if (nDiff >= 1) {
      FormSigma(pSigmaCC, DenC.pRhoGrd, DenC.pRhoGrd);
      if (Args.OpenShell()){
         FormSigma(pSigmaCO, DenC.pRhoGrd, DenO.pRhoGrd);
         FormSigma(pSigmaOO, DenO.pRhoGrd, DenO.pRhoGrd);
      }
   }

   // allocate output data: functional kernels and derivatives
   // of kernel with respect to input data.
   // zk,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,vupsilonc,vupsilono
//    FDftiResultSet
//       ResC = {0}, ResO = {0};
// ^- why does g++ warn about this? This is perfectly well defined in both C and C++! Annoying.
   FDftiResultSet
      ResC = FDftiResultSet(), ResO = FDftiResultSet();
   ResC.Init(nGridPt, nDiff, Args.NeedTau(), Args.NeedUpsilon(), 0, Mem);
   ResO.Init(nGridPt, nDiff, Args.NeedTau(), Args.NeedUpsilon(), ResC.pVdSigmaXY, Mem);

   // branch out to fortran in order to evaluate the functional
   // and the split energy contributions (=sum zk(g) wt(g), but split according to functional)
   RESUME_CPU_CLOCK(CLOCK_EvalDftfun);
   double
      *pTmp;
   Mem.Alloc(pTmp, 30*nGridPt);
   DFTFUN_CXX(Args.pXcFn, 1, Args.OpenShell(), nDiff, nGridPt, pGridWt,
              DenC.pRho, DenO.pRho, pSigmaCC, pSigmaCO, pSigmaOO, DenC.pTau, DenO.pTau, DenC.pUpsilon, DenO.pUpsilon,
              // density functional outputs (summed over functionals)
              ResC.pZk, ResC.pVdRho, ResO.pVdRho, ResC.pVdSigma, ResC.pVdSigmaXY, ResO.pVdSigma, ResC.pVdTau, ResO.pVdTau, ResC.pVdUpsilon, ResO.pVdUpsilon,
              // energy outputs (split according to functions; one output per ndftu; accumulated)
              pDfuEnergies,
              // buffer (for accumulating zk)
              pTmp,
              // flags specifying input data
              (FORTINT)(DenC.pTau != 0), (FORTINT)(DenC.pUpsilon != 0));
   Mem.Free(pTmp);
   PAUSE_CPU_CLOCK(CLOCK_EvalDftfun);

   if (Args.Flags & DFTI_CounterDensityXcIntegral) {
      // subtract grid-computed density-vxc integral from energy contribution (experimental)
      assert(!Args.OpenShell());
      pDfuEnergies[0] -= Dot3(DenC.pRho, ResC.pVdRho, pGridWt, nGridPt);
      if (0) {
         assert(Args.pXcFn->nOrder() == 0); // rho * VdRho term only enough for LDA (apparently..)
      } else {
         assert(Args.pXcFn->nOrder() <= 1);
         if (Args.pXcFn->nOrder() >= 1) {
            assert(nDiff >= 1);
            // this works...for closed-shell PBE and TPSS. Despite that I
            // admittedly did not think about this very hard.
            //
            // It just seems like on formal reasons what we do here is subtract
            // linear terms in the 1st order Taylor approximation of exc at ρ=0
            // (see dftb.pdf writeup):
            //
            //     exc[ρ=0] ≈ exc[ρ_cur] + (0 - ρ_cur) * D[exc[ρ],ρ] |ρ=ρ_cur + O[ρ^2]
            //
            // The linear terms are already handled as part of tr(f,γ) (i.e., the band
            // structure energy = sum of orbital energies) in this formalism.
            //
            // The concrete terms then correspond to the ones in FDensitySet:AccXcMatrix.
            //
            // We could probably even get the open-shell terms going rather
            // straight-forwardly. At least the part in here. In CtRhf.cpp it
            // might be a bit more tricky to identify the right kind of Fock
            // matrix to put into tr(f,γ), and adjust the ExchO factors and γ_open
            // accordingly.
            pDfuEnergies[0] -= 2 * Dot3(pSigmaCC, ResC.pVdSigma, pGridWt, nGridPt);

            if (DenC.pTau != 0) {
               pDfuEnergies[0] -= Dot3(DenC.pTau, ResC.pVdTau, pGridWt, nGridPt);
            }
         }
      }
   }

   if (Args.NeedTau()) {
      // xc/gradient code written for tau = \sum_i <grad phi_i, grad phi_i>,
      // but the actual kinetic energy is a factor of 2 lower. Fix this.
      for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
         ResC.pVdTau[iPt] *= 0.5;
      if (Args.OpenShell()) {
         for (size_t iPt = 0; iPt < nGridPt; ++ iPt)
            ResO.pVdTau[iPt] *= 0.5;
      }
      // ^- hm.... should this not be * 2 instead of * 0.5?
      // ...does seem to work like this (can reproduce Turbomole's Ne/TPSS energy).
      // Not sure what to think of this.
   }
   if (Args.AddBackgroundDensities())
      AdjustXcForEmbeddedCalc(ResC, ResO, DenC, DenO, Mem);

   if (Args.pfElecTotal) {
      fElec[0] = DenC.fElec;
      if (Args.OpenShell())
         fElec[1] = DenO.fElec;
   }
//       fElec = Dot(DenC.pRho, pGridWt, nGridPt);

   if (Args.MakeXc()) {
      // assemble & accumulate xc matrices.
      RESUME_CPU_CLOCK(CLOCK_FormXcMat);
      DenC.AccXcMatrix(pFockC, pOrbVal, ResC, DenO.pRhoGrd, DenInpC, pGridWt, Mem);
      if (Args.OpenShell())
         DenO.AccXcMatrix(pFockO, pOrbVal, ResO, DenC.pRhoGrd, DenInpO, pGridWt, Mem);
      PAUSE_CPU_CLOCK(CLOCK_FormXcMat);
   }

//#ifdef INCLUDE_OPTIONALS
   if (Args.MakeGradient()) {
      // assemble & accumulate analytic gradient
      RESUME_CPU_CLOCK(CLOCK_CoreGrad);
      AccCoreGradient(pGradient, pOrbVal, pCenterIndices, ResC,
         DenC, DenInpC, DenO.pRhoGrd, ResC.pVdSigmaXY, Args.MakeGridGradient());
      if (Args.OpenShell())
         AccCoreGradient(pGradient, pOrbVal, pCenterIndices, ResO,
            DenO, DenInpO, DenC.pRhoGrd, ResC.pVdSigmaXY, Args.MakeGridGradient());
      PAUSE_CPU_CLOCK(CLOCK_CoreGrad);

      assert(!Args.MakeGridGradient());
//       if (Args.MakeGridGradient()) {
//          RESUME_CPU_CLOCK(CLOCK_GridWtDeriv);
//          AccGridWtGradient(pGradient, ResC.pZk);
//          PAUSE_CPU_CLOCK(CLOCK_GridWtDeriv);
//       }
   }
//#endif // INCLUDE_OPTIONALS
}


// #ifdef HAVE_GRID_KERNEL
bool Vec3Eq(double const *pA, double const *pB) {
   return pA[0] == pB[0] && pA[1] == pB[1] && pA[2] == pB[2];
}

size_t GetNumBfDerivComps(FEvalBfDerivType DerivType, size_t *iLaplaceComp0)
{
   assert(DERIVCOMP_ValueOnly == 0 && DERIVCOMP_Value_dFirst == 1 && DERIVCOMP_Value_dFirst_dSecond == 2 &&
          DERIVCOMP_Value_Laplace == 3 && DERIVCOMP_Value_dFirst_Laplace == 4 && DERIVCOMP_Value_dFirst_dSecond_Laplace_dFirst == 5);

   static const size_t
      s_nGridDerivComp1[6] = {1,4,10, 2,5,14},
      // offset for derivative components of Laplacian of basis function. The
      // derivatives of Laplace mu(r) then follow in standard order (e.g., d/dx
      // Laplace mu(r) comes as next element behind Laplace mu(r) itself)
      s_nGridDerivLaplaceComp0[6] = {size_t(-1),size_t(-1),size_t(-1), 1,4,10};

   assert(size_t(DerivType) < sizeof(s_nGridDerivComp1)/sizeof(s_nGridDerivComp1[0]));

   // note: the number of components in the full set of cartesian derivatives is
   // ((nDiffBf+1)*(nDiffBf+2)*(nDiffBf+3))/6. This is how FDftiJob used to
   // compute that quantity (nComp) before separate Laplacian routines and this
   // function were introduced.

   if (iLaplaceComp0)
      *iLaplaceComp0 = s_nGridDerivLaplaceComp0[size_t(DerivType)];

   return s_nGridDerivComp1[size_t(DerivType)];
}


// pOrbValAo:
// - input format: nGridPt * nComp * nAo (i.e., must hold room for this many scalars)
// - output format: nGridPt * nComp * nMap (screened out basis functions do not get emitted)
void IrEvalBfn(double *pOrbValAo, size_t *pMap, size_t &nMap, double const *pGridPt,
      size_t iStGridPt, size_t nGridPt, FEvalBfDerivType DerivType, FRawBasis const *pBasis,
      double ThrOrb, double LogThrOrb, FMemoryStack &Mem)
{
   TMemoryLock<char>
      pBaseOfMemory(0, &Mem);
   FMatrixView const
      Grid(const_cast<double*>(pGridPt), iStGridPt, nGridPt);
   size_t
//       nGridPt = Grid.nCols,
      nComp = GetNumBfDerivComps(DerivType),
      nAo = pBasis->nFn();
   IR_SUPPRESS_UNUSED_WARNING(nAo);

   bool
      MakeLaplaceToo = false;
   unsigned
      DerivOrder = unsigned(DerivType);
   if (DerivOrder >= unsigned(DERIVCOMP_Value_Laplace)) {
      MakeLaplaceToo = true;
      DerivOrder -= unsigned(DERIVCOMP_Value_Laplace);
   }

   // convert basis format for low level driver input...
   FRawBasis::FShellArray
      &ShellsA = const_cast<FRawBasis*>(pBasis)->Shells;

   nMap = 0;
//    assert(Grid.nRowSt == 1 && Grid.nRows == 3);
   assert(Grid.nRowSt == 1);

   size_t
      iFnBase = 0;
   ir::FRawShell
      *pLastBf = &ShellsA[0];
   for (ir::FRawShell *pFirstBf = &ShellsA[0]; pFirstBf != &ShellsA[0] + ShellsA.size(); pFirstBf = pLastBf) {
      pLastBf = pFirstBf;
      size_t nFn = 0;
      while (pLastBf != &ShellsA[0] + ShellsA.size() && Vec3Eq(pLastBf->vCen, pFirstBf->vCen)) {
         nFn += pLastBf->nFn();
         ++pLastBf;
      }
      size_t
         nBfStride = nGridPt*nComp;
      EvalShellGroupOnGrid(&pOrbValAo[nBfStride*nMap], nGridPt, nBfStride, // comp stride, bf stride
         pMap, nMap, iFnBase, &Grid(0,0), Grid.nColSt, nGridPt,
         pFirstBf, pLastBf, &pFirstBf->vCen[0],
         DerivOrder, MakeLaplaceToo, ThrOrb, LogThrOrb, Mem);
      iFnBase += nFn;
   }
   assert(iFnBase == nAo && nMap <= nAo);
   assert(pLastBf == &ShellsA[0] + ShellsA.size());
}
// #endif




void FDftiJob::EvalBfn(double *&pOrbVal, size_t *&pMap, size_t &nMap, ptrdiff_t *&pCenterIndices)
{
   size_t
      nBf = Args.nBf();
   Mem.Alloc(pCenterIndices, nBf); // used in gradient evaluation
   Mem.Alloc(pMap, nBf);
   Mem.Alloc(pOrbVal, nGridPt * nComp * nBf);
   IrEvalBfn(pOrbVal, pMap, nMap, &pGridPt[0], 3, nGridPt, DerivTypeBf, Args.pBasis, Args.ThrOrb, Args.LogThrOrb, Mem);

   // assemble indices of basis function centers.
   {
      FRawBasis
         *pBasis = Args.pBasis;
      TMemoryLock<size_t>
         AllCenters(nBf, &Mem);
      size_t iFnTotal = 0;
      for (size_t iSh = 0; iSh < pBasis->Shells.size(); ++ iSh) {
         size_t iCen = pBasis->ShellCenters[iSh];
         size_t nFn = pBasis->nFn(iSh);
         for (size_t iFn = 0; iFn < nFn; ++ iFn) {
            AllCenters[pBasis->iFn(iSh) + iFn] = iCen;
         }
         assert_rt(iFnTotal == pBasis->iFn(iSh));
         iFnTotal += nFn;
      }
      assert_rt(iFnTotal == nBf);
      for (size_t iMap = 0; iMap < (size_t)nMap; ++ iMap) {
         assert_rt(pMap[iMap] < nBf);
         pCenterIndices[iMap] = AllCenters[pMap[iMap]];
      }
   }
}



void AccXc(FDftiArgs &Args, FLog &Log, FTimerSet *pTimers, FMemoryStack &Mem_)
{
//    TIME_SECTION(pTimers, 0x37, "ACC-XC");
   TMemoryLock<double> pFreeMe(0, &Mem_);

   size_t
      nDftFunc = 1,
      nFockSize = Args.nFockSize();
   if (!Args.MakeXc())
      nFockSize = 0;

   size_t nJobs = Args.pDftGrid->GridBlocks.size();
   double DummyTimings[CLOCK_TimerCount];
   Args.pTimings = &DummyTimings[0];

   assert(Args.OpenShell() == (Args.pFockO != 0));

   // clear output & timings
   FOmpAccBlock
      FockC(Args.pFockC, nFockSize, 0, Mem_),
      FockO(Args.OpenShell()? Args.pFockO : 0, nFockSize, 0, Mem_),
      DfuEnergies(Args.pDfuEnergies, nDftFunc, OMPACC_ClearTarget, Mem_),
      Timings(Args.pTimings, CLOCK_TimerCount, OMPACC_ClearTarget, Mem_),
      fElecTotal(Args.pfElecTotal, 2, OMPACC_ClearTarget, Mem_),
      // ^- combine with other scalar quantities? Like energies?
      Gradient(Args.MakeGradient() ? Args.pGradient : 0, 3*Args.nCenters(), 0, Mem_);

   FMemoryStackArray MemStacks(Mem_);

   bool Aborted = false;
   #pragma omp parallel for schedule(dynamic,1)
   //for (size_t iJob = 0; iJob < nJobs; ++ iJob ) {
   for (int iJob__ = 0; iJob__ < int(nJobs); ++ iJob__ ) {
      size_t iJob = size_t(iJob__); // for OpenMP.
      FMemoryStack &Mem = MemStacks.GetStackOfThread();

      if (!Log.StatusOkay())
         Aborted = true;
      #pragma omp flush (Aborted)
      if (!Aborted) {
         // ^- both "break" and "continue" seem to be problemaic..
         TMemoryLock<char>
            pBaseOfMemory(0, &Mem);
         FDftGrid::FGridBlock const
            *pGridBlock = &Args.pDftGrid->GridBlocks[iJob];

         FDftiJob
            Job(DfuEnergies.pTls(), FockC.pTls(), FockO.pTls(), Gradient.pTls(), pGridBlock, Args, Mem);
         Job.Run();
         for (size_t i = 0; i < 2; ++ i)
            fElecTotal.pTls()[i] += Job.fElec[i];
      }
   }

   FockC.Join();
   FockO.Join();
   DfuEnergies.Join();
   Timings.Join();
   fElecTotal.Join();
   Gradient.Join();

   Log.CheckStatus();
   assert_rt(!Aborted); // you should not be here if this was set.

   if (Args.MeasureTime()) {
      Log.WriteTiming("DFT (basis fn.)", Args.pTimings[CLOCK_EvalBfn]);
      Log.WriteTiming("DFT (densities)", Args.pTimings[CLOCK_FormRho]);
      if (Args.MakeXc())
         Log.WriteTiming("DFT (xc matrix)", Args.pTimings[CLOCK_FormXcMat]);
      Log.WriteTiming("DFT (functional)", Args.pTimings[CLOCK_EvalDftfun]);
      if (Args.MakeGradient())
         Log.WriteTiming("DFT (core grad.)", Args.pTimings[CLOCK_CoreGrad]);
      if (Args.MakeGradient() && Args.MakeGridGradient())
         Log.WriteTiming("DFT (grid grad.)", Args.pTimings[CLOCK_GridWtDeriv]);
   }
   IR_SUPPRESS_UNUSED_WARNING(pTimers); // <- should probably move local timing code into this...
}



} // namespace dfti
