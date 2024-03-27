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

#include "CtVolumeProperty.h"

#include <stdexcept>
#include <iostream>
#include <string.h> // for memcpy

#include "format.h"
// #include "CxTiming.h" // for embdedded HF timer.
// #include "CxPhysicalUnits.h"
#include "CxMemoryStack.h"
#include "CxOpenMpProxy.h"

// #include "CxIo.h"
// #include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "Ir.h"



namespace ct {




bool Vec3Eq(double const *pA, double const *pB) {
   return pA[0] == pB[0] && pA[1] == pB[1] && pA[2] == pB[2];
}


void MakeGridValuesCore(double *pOut, size_t iGridCompSt, FMatrixView Grid, FMatrixView Orb, uint GridDxOrder, FRawBasis::FShellArray &ShellsA, size_t nGridPt, size_t iGridPt, FMemoryStack &Mem)
{
   TMemoryLock<double> pFreeMe1(0, &Mem);
   if (nGridPt == 0)
      return;
   size_t
      nComp = 1,
      nAo = Orb.nRows,
      nOrb = Orb.nCols;

   if (GridDxOrder == 1) nComp = 4;
   if (GridDxOrder == 2) nComp = 10;

   size_t
      *pMap, nMap = 0, iFnBase = 0;
   double
      *pOrbValAo, *pOrbsCompressed;
   Mem.Alloc(pOrbValAo, nGridPt * nComp * nAo);
   Mem.Alloc(pOrbsCompressed, nAo * nOrb);

   ir::FRawShell
      *pLastBf = &ShellsA[0];
   Mem.Alloc(pMap, nAo);
   for (ir::FRawShell *pFirstBf = &ShellsA[0]; pFirstBf != &ShellsA[0] + ShellsA.size(); pFirstBf = pLastBf) {
      pLastBf = pFirstBf;
      size_t nFn = 0;
      while (pLastBf != &ShellsA[0] + ShellsA.size() && Vec3Eq(pLastBf->vCen, pFirstBf->vCen)) {
         nFn += pLastBf->nFn();
         ++pLastBf;
      }
      size_t
         nBfStride = nGridPt*nComp;
      ir::EvalShellGroupOnGrid(&pOrbValAo[nBfStride*nMap], nGridPt, nBfStride, // comp stride, bf stride
         pMap, nMap, iFnBase, &Grid(0,iGridPt), Grid.nColSt, nGridPt,
         pFirstBf, pLastBf, &pFirstBf->vCen[0],
         GridDxOrder, false, 1e-10, 40, Mem);
      iFnBase += nFn;
   }
   assert(iFnBase == nAo && nMap <= nAo);
   assert(pLastBf == &ShellsA[0] + ShellsA.size());

   FMatrixView
      OrbCmpr(pOrbsCompressed, nMap, nOrb),
      OrbValAo(pOrbValAo, nGridPt*nComp, nMap);
   // compress input orbital matrix: delete basis functions not occurring in nMap.
   for (size_t iOrb = 0; iOrb != nOrb; ++ iOrb)
      for (size_t iMap_ = 0; iMap_ != nMap; ++ iMap_)
         OrbCmpr(iMap_, iOrb) = Orb(pMap[iMap_], iOrb);
   for (uint iComp = 0; iComp < nComp; ++ iComp) {
      FMatrixView
         Out(&pOut[iGridPt+iGridCompSt*iComp], nGridPt, nOrb, 1, iGridCompSt*nComp),
         OrbValAoComp = OrbValAo;
      if (nMap == 0) {
         Out.Clear();
      } else {
         OrbValAoComp.nRows = nGridPt;
         OrbValAoComp.pData += nGridPt * iComp;
         Mxm(Out, OrbValAoComp, OrbCmpr);
      }
   }
}


void MakeGridValues(double *pOut, FMatrixView Grid, FMatrixView Orb, unsigned GridDxOrder, FBasisSet const *pBasis, FMemoryStack &Mem_)
{
   TMemoryLock<double> pFreeMe0(0, &Mem_);
   size_t
      nGridPt_Total = Grid.nCols,
      nGridStep = 64;
   if (pBasis->nFn() != Orb.nRows)
      throw std::runtime_error("ct::MakeGridValues(): Input orbitals not consistent with orbital basis set.");

   // convert basis format for low level driver input...
   FRawBasis::FShellArray
      &ShellsA = const_cast<FBasisSet*>(pBasis)->pRawBasis->Shells;

   assert(Grid.nRowSt == 1 && (Grid.nRows == 3 || Grid.nRows == 4 || Grid.nRows == 5)); // grid format must be (x,y,z), (x,y,z,w), or (x,y,z,w,iat)

   if (omp_get_max_threads() == 1) {
      // non-parallel version. Used if parallel region is on the outside. Otherwise all kinds of things
      // get messed up... although technically, this should be working without this... (it doesn't).
      int nGridPt_Batches = int((nGridPt_Total + nGridStep - 1)/nGridStep);
      for ( int iGridPt_Batch = 0; iGridPt_Batch < nGridPt_Batches; ++ iGridPt_Batch ) {
         size_t iGridPt = size_t(iGridPt_Batch) * nGridStep;
         // ^- that's for OMP 2.0 which allows only 'int' variables for loop control.
         size_t
            nGridPt = nGridStep;
         if ( iGridPt + nGridPt > nGridPt_Total )
            nGridPt = nGridPt_Total - iGridPt;
         MakeGridValuesCore(pOut, nGridPt_Total, Grid, Orb, GridDxOrder, ShellsA, nGridPt, iGridPt, Mem_);
      }
   } else {
      FMemoryStackArray MemStacks(Mem_);
//       #pragma omp parallel for schedule(dynamic)
//       for (size_t iGridPt = 0; iGridPt < nGridPt_Total; iGridPt += nGridStep) {
      int nGridPt_Batches = int((nGridPt_Total + nGridStep - 1)/nGridStep);
      #pragma omp parallel for schedule(dynamic)
      for ( int iGridPt_Batch = 0; iGridPt_Batch < nGridPt_Batches; ++ iGridPt_Batch ) {
         size_t iGridPt = size_t(iGridPt_Batch) * nGridStep;
         // ^- that's for OMP 2.0 which allows only 'int' variables for loop control.
         FMemoryStack &Mem = MemStacks.GetStackOfThread();
         TMemoryLock<double> pFreeMe1(0, &Mem);
         size_t
            nGridPt = nGridStep;
         if ( iGridPt + nGridPt > nGridPt_Total )
            nGridPt = nGridPt_Total - iGridPt;
         MakeGridValuesCore(pOut, nGridPt_Total, Grid, Orb, GridDxOrder, ShellsA, nGridPt, iGridPt, Mem);
      }
   }
}





FDensityEvalContext::FDensityEvalContext(FMatrixView Orbs, double const *pOcc, FBasisSetPtr pOrbBasis)
   : m_Orbs(Orbs), m_pOrbBasis(pOrbBasis)
{
   m_Occ.assign(pOcc, pOcc + Orbs.nCols);
}


FDensityEvalContext::~FDensityEvalContext()
{
}


void FDensityEvalContext::ComputeDensity(double *pOut_, FMatrixView Grid_, FMemoryStack &Mem_)
{
   // Note: this is mostly simplified c/p code from IvIsoSurface. Main reason being that we need
   // some other things here, too...

   // disable nested parallelism... we put out parallel evaluation thing already here because otherwise
   // we might need to make too large grid blocks (nGridPtTotal, which might consist of an entire grid with
   // millions of points, times the number of orbitals times 4 components, times size of double...).
   // For this reason, we should disable the additional parallelism level in MakeGridValues.
   omp_set_nested(0);

   size_t
      nDxComp = 1,
      GridDxOrder = 0;

   size_t
      nGridPt_Total = Grid_.nCols;
   if (m_Orbs.nCols == 0) {
      // no data... make zeros and return.
      memset(pOut_, 0, sizeof(*pOut_) * nDxComp * nGridPt_Total);
      return;
   } else {
      FMemoryStackArray MemStacks(Mem_);
      size_t
         nGridStep = 64;
      int nGridPt_Batches = int((nGridPt_Total + nGridStep - 1)/nGridStep);
      #pragma omp parallel
      {
         omp_set_num_threads(1); // <- my understanding is that this is supposed to apply only to the inner region.
         #pragma omp for schedule(dynamic)
         for (int iGridPt_Batch = 0; iGridPt_Batch < nGridPt_Batches; ++ iGridPt_Batch) {

            size_t iGridPtBeg = size_t(iGridPt_Batch) * nGridStep;
            // ^- that's for OMP 2.0 which allows only 'int' variables for loop control.
            FMemoryStack &Mem = MemStacks.GetStackOfThread();
            TMemoryLock<double> pFreeMe1(0, &Mem);
            size_t
               nGridPt = nGridStep;
            if (iGridPtBeg + nGridPt > nGridPt_Total)
               nGridPt = nGridPt_Total - iGridPtBeg;
            if (nGridPt == 0)
               continue;
            double
               *pOut = &pOut_[iGridPtBeg];

            ct::FMatrixView
               Grid(&Grid_(0,iGridPtBeg), Grid_.nRows, nGridPt, Grid_.nRowSt, Grid_.nColSt);
            // evaluate the orbitals supporting the current density.
            FStackMatrix
               OrbValues(nGridPt * nDxComp, m_Orbs.nCols, &Mem);
            ct::MakeGridValues(OrbValues.pData, Grid, m_Orbs, GridDxOrder, &*m_pOrbBasis, Mem);

            // compute density
            // note: not critical since each thread handles a different block of grid points
            // on the output quantity. They do not overlap.
            for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt)
               pOut[iGridPt] = 0.;
            for (size_t iOrb = 0; iOrb < m_Orbs.nCols; ++ iOrb) {
               double w = m_Occ[iOrb];
               for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt) {
                  pOut[iGridPt] += w * sqr(OrbValues.pData[iGridPt + nGridPt * iOrb]);
               }
            }
         }
      }
   }
}



void FDensityEvalContext::AccFragmentOrbitalOverlap(FMatrixView Out, FMatrixView Grid_, FFragmentWeightFn *pWeightFn, FMemoryStack &Mem_)
{
   size_t
      nOrb = m_Orbs.nCols;
   assert(Out.nRows == nOrb && Out.nCols == nOrb);
   assert(Grid_.nRows == 4);

   if (m_Orbs.nCols == 0 || Grid_.nCols == 0) {
      // nothing to accumulate.
      return;
   }

   // disable nested parallelism.
   omp_set_nested(0);

   size_t
      nDxComp = 1,
      GridDxOrder = 0;

   size_t
      nGridPt_Total = Grid_.nCols;
   {
      assert_rt(Out.nRowSt == 1 && Out.nColSt == Out.nRows && Out.nRows == nOrb);
      FOmpAccBlock
         SijAccBlock(Out.pData, nOrb*nOrb, 0, Mem_);
      FMemoryStackArray
         MemStacks(Mem_);
      size_t
         nGridStep = 64;
      int nGridPt_Batches = int((nGridPt_Total + nGridStep - 1)/nGridStep);
      #pragma omp parallel
      {
         omp_set_num_threads(1); // <- my understanding is that this is supposed to apply only to the inner region.
         #pragma omp for schedule(dynamic)
         for (int iGridPt_Batch = 0; iGridPt_Batch < nGridPt_Batches; ++ iGridPt_Batch) {

            size_t iGridPtBeg = size_t(iGridPt_Batch) * nGridStep;
            // ^- that's for OMP 2.0 which allows only 'int' variables for loop control.
            FMemoryStack &Mem = MemStacks.GetStackOfThread();
            TMemoryLock<double> pFreeMe1(0, &Mem);
            size_t
               nGridPt = nGridStep;
            if ( iGridPtBeg + nGridPt > nGridPt_Total )
               nGridPt = nGridPt_Total - iGridPtBeg;
            if (nGridPt == 0)
               continue;

            ct::FMatrixView
               Grid(&Grid_(0,iGridPtBeg), Grid_.nRows, nGridPt, Grid_.nRowSt, Grid_.nColSt);
            // evaluate the orbitals supporting the current density.
            FStackMatrix
               OrbValues(nGridPt * nDxComp, m_Orbs.nCols, &Mem);
            ct::MakeGridValues(OrbValues.pData, Grid, m_Orbs, GridDxOrder, &*m_pOrbBasis, Mem);

            // compute values of external fragment weight contribution, if provided
            TMemoryLock<double>
               pFragmentWeights(nGridPt, &Mem);
            if (pWeightFn) {
               pWeightFn->EvalWeights(&pFragmentWeights[0], Grid, Mem);
            } else {
               // no external weight function. Just set all weights to 1.
               for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt)
                  pFragmentWeights[iGridPt] = 1.;
            }

            // absorb integration weights into orbital values.
            for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt) {
               double
                  w = Grid(3, iGridPt) * pFragmentWeights[iGridPt],
                  wsqrt = std::sqrt(w);
               assert(w >= 0);
               for (size_t iOrb = 0; iOrb != m_Orbs.nCols; ++ iOrb) {
                  OrbValues(iGridPt, iOrb) *= wsqrt;
               }
            }

            // compute contribution to S[i,j] = \sum_g w_g phi(r_g,i) phi(r_g,j)
            ct::FMatrixView
               SijAcc(SijAccBlock.pTls(), nOrb, nOrb);
            SyrkTN(SijAcc, OrbValues, 1.0, MXM_Add);
         }
      }
      SijAccBlock.Join();
   }
}



} // namespace ct
