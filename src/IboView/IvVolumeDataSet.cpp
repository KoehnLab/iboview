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

#include "Iv.h"
#include <iostream>
#include <stdexcept>

// #include <QFileInfo>
#include <cmath>
#include <algorithm>
#include <set>
#include <fstream> // for orbital basis export

#include <QTextStream>

#include "format.h"
#include "CxIo.h"
#include "CxColor.h"
#include "CxPhysicalUnits.h"
#include "CxOpenMpProxy.h"
// #include "CxXyzFrameIo.h"

#include "IvDocument.h"
#include "IvSettings.h"
// #include "IvFileConvert.h"
#include "IvVolumeDataSet.h"
#include "IvOrbital.h"

// #include "CtInt1e.h"
// #include "IrAmrr.h"
#include "CtMatrix.h"
// #include "CtRhf.h"
// #include "CtBasisLibrary.h"
#include "CtVolumeProperty.h"


double g_fSupportBox_EnlargeRms = 4.; // 2 should be good for 95% coverage.


FVolumeVisualConfig::FVolumeVisualConfig()
{
//    bColorSet = false;
//    iIsoType = ISOTYPE_Undefined;
//    fIsoValue = -1.;
   // ^- both not explicitly defined -- take whatever is currently the default when the orbtial is first traced.
}

FVolumeVisualConfig::FVolumeVisualConfig(FVolumeVisualConfigDetailsPtr pDetails_)
   : pDetails(pDetails_)
{}

FVolumeVisualConfigDetails::~FVolumeVisualConfigDetails()
{}

FVolumeVisualConfig::~FVolumeVisualConfig()
{}


FVolumeVisualConfigDetailsPtr MakeDefaultTwoPhaseIsoSurfaceConfig(FDocument *pDocument, double fIsoValue_)
{
   if (pDocument == 0) {
      IV_NOTIFY(NOTIFY_Error, "No parent document set. Very strange!");
      return FVolumeVisualConfigDetailsPtr(0);
   }
   //    uint32_t cFullAlpha = 0xff000000;
   uint32_t
      cIsoPlus, cIsoMinus;
   pDocument->GetNextOrbitalColors(cIsoPlus, cIsoMinus);
//    AssignPlusMinusColors(cIsoPlus, cIsoMinus, fIsoValue_);
   return FVolumeVisualConfigDetailsPtr(new FTwoPhaseIsoSurfaceConfig(ISOTYPE_Relative, fIsoValue_, cIsoPlus, cIsoMinus));
}


void FTwoPhaseIsoSurfaceConfig::AssignTo(FIsoSurfaceSettings &IsoSurfOpt_) const
{
   IsoSurfOpt_.IsoType = iIsoType;
   IsoSurfOpt_.IsoValues.clear();
   IsoSurfOpt_.IsoValues.push_back(FIsoThresholdEntry(fIsoValue, cIsoPlus));
   IsoSurfOpt_.IsoValues.push_back(FIsoThresholdEntry(-fIsoValue, cIsoMinus));
}


void FMultiPhaseIsoSurfaceConfig::AssignTo(FIsoSurfaceSettings &IsoSurfOpt_) const
{
   IsoSurfOpt_.IsoType = iIsoType;
   IsoSurfOpt_.IsoValues = IsoValues;
}

// void FVolumeVisualConfig::AssignDefaultColor(FDocument *pDocument, double fIsoValue_)
// {
//    if (pDocument == 0)
//       return IV_NOTIFY(NOTIFY_Error, "No parent document set. Very strange!");
//    //    uint32_t cFullAlpha = 0xff000000;
//    uint32_t
//       cIsoPlus, cIsoMinus;
//    pDocument->GetNextOrbitalColors(cIsoPlus, cIsoMinus);
//    AssignPlusMinusColors(cIsoPlus, cIsoMinus, fIsoValue_);
// }
//
// void FVolumeVisualConfig::AssignPlusMinusColors(uint32_t cIsoPlus_, uint32_t cIsoMinus_, double fIsoValue_)
// {
// //    pDetails = new FTwoPhaseIsoSurfaceConfig(FIsoType IsoType_, float fIsoValue_, uint32_t cIsoPlus_, uint32_t cIsoMinus_)
//
// //    if (fIsoValue_ != -1. && IsoValues.empty()) {
// //       IsoValues.push_back(FIsoThresholdEntry(fIsoValue_, cIsoPlus_));
// //       IsoValues.push_back(FIsoThresholdEntry(-fIsoValue_, cIsoMinus_));
// //    } else {
// //       for (size_t i = 0; i < IsoValues.size(); ++ i) {
// //          if (IsoValues[i].fIsoValue <= 0.)
// //             IsoValues[i].dwColor = cIsoPlus_;
// //          else
// //             IsoValues[i].dwColor = cIsoMinus_;
// //       }
// //    }
// //    bColorSet = true;
// }

// uint32_t FVolumeVisualConfig::cIsoPlus() const
// {
//    for (size_t i = 0; i < IsoValues.size(); ++ i) {
//       if (IsoValues[i].fIsoValue >= 0.)
//          return IsoValues[i].dwColor;
//    }
//    return 0xffff0000;
// }
//
// uint32_t FVolumeVisualConfig::cIsoMinus() const
// {
//    for (size_t i = 0; i < IsoValues.size(); ++ i) {
//       if (IsoValues[i].fIsoValue < 0.)
//          return IsoValues[i].dwColor;
//    }
//    return 0xff0000ff;
// }


void FVolumeVisualConfig::Link(FVolumeDataSet *pVolume)
{
   assert(LinkedVolumes.find(pVolume) == LinkedVolumes.end());
   LinkedVolumes.insert(pVolume);
}

void FVolumeVisualConfig::Unlink(FVolumeDataSet *pVolume)
{
   assert(LinkedVolumes.find(pVolume) != LinkedVolumes.end());
   LinkedVolumes.erase(pVolume);
}

void FVolumeVisualConfig::UpdateLinkedRepresentations(uint32_t Flags, FView3d *pView3d)
{
   FVolumeChain::iterator
      it;
   for (it = LinkedVolumes.begin(); it != LinkedVolumes.end(); ++ it) {
      if (Flags & UPDATE_InvalidateIsoSettings) {
         // completely delete the data set's mesh and rebuild it from scratch.
         (*it)->InvalidateRenderCache();
      }
      if (Flags & UPDATE_InvalidateColors) {
         (*it)->InvalidateColors();
      }
      if (Flags & UPDATE_Rebuild) {
         assert(pView3d != 0);
         (*it)->BuildRenderCache(pView3d);
      }
   }
}



FVolumeDataSet::FVolumeDataSet(QString const &Desc_, ct::FAtomSetPtr pAtoms_, FDocument *pDocument_, FVolumeVisualConfigPtr pRefVisConfig_)
   : FBase(Desc_, pAtoms_, pDocument_), pVisConfig(pRefVisConfig_)
{
   // make a new visual representation object if none was supplied from the outside.
   if (pVisConfig.get() == 0)
      pVisConfig = FVolumeVisualConfigPtr(new FVolumeVisualConfig());
   pVisConfig->Link(this);
}

FVolumeDataSet::~FVolumeDataSet()
{
   pVisConfig->Unlink(this);
}

void FVolumeDataSet::LinkVisualConfig(FVolumeVisualConfigPtr p)
{
   if (pVisConfig != p) {
      pVisConfig = p;
      pVisConfig->Link(this);
   }
}

void FVolumeDataSet::InvalidateRenderCache()
{
   pGlMesh = 0;
}


void FVolumeDataSet::InvalidateColors()
{
   FVolumeDataSet *pVolume = this;

   FTwoPhaseIsoSurfaceConfig const
      *pVisDetTwoPhase = dynamic_cast<FTwoPhaseIsoSurfaceConfig const*>(pVolume->pVisConfig->pDetails.get());
   FMultiPhaseIsoSurfaceConfig const
      *pVisDetMultiPhase = dynamic_cast<FMultiPhaseIsoSurfaceConfig const*>(pVolume->pVisConfig->pDetails.get());

   // update the color also in the mesh, if one already exists.
   if (pVolume->pGlMesh.get()) {
      for (uint i = 0; i < pVolume->pGlMesh->Vertices.size(); ++ i) {
         FBaseVertex
            &vx = pVolume->pGlMesh->Vertices[i];
         if (pVisDetTwoPhase) {
            if (vx.iRole == 0) vx.dwColor = pVisDetTwoPhase->cIsoPlus;
            if (vx.iRole == 1) vx.dwColor = pVisDetTwoPhase->cIsoMinus;
         } else if (pVisDetMultiPhase) {
            size_t iPhase = std::min(size_t(vx.iRole), size_t(pVisDetMultiPhase->IsoValues.size()));
            vx.dwColor = pVisDetMultiPhase->IsoValues[iPhase].dwColor;
         }
      }
      pVolume->pGlMesh->Invalidate();
   }
}





uint32_t FDataSet::GetBaseColor() const
{
   return 0;
}

bool FDataSet::DependsOnWaveFunction() const
{
   return false;
}

uint32_t FTwoPhaseIsoSurfaceConfig::GetBaseColor() const
{
   uint32_t dwColor = ct::irgb(ct::FColor(0.5f*(ct::FColor(cIsoPlus) + ct::FColor(cIsoMinus))).uint32());
   return dwColor;
}

uint32_t FMultiPhaseIsoSurfaceConfig::GetBaseColor() const
{
   // return first color if present, white otherwise.
   if (!IsoValues.empty())
      return IsoValues[0].dwColor;
   else
      return 0xffffff;
}



QString FVolumePropertyInfo::MakeDesc() const
{
   switch (this->Type) {
      case PROPERTY_Density: {
         return QString("Density");
         break;
      }
      case PROPERTY_SpinDensity: {
         return QString("SpinDensity");
         break;
      }
      default: {
         assert(0);
         return QString("Unknown3DProperty");
         break;
      }
   }
}

FVolumeProperty::FVolumeProperty(FFrame *pWf_, FVolumePropertyInfo const &info_, FVolumeVisualConfigPtr pRefVisConfig_)
   : FBase(info_.MakeDesc(), ct::FAtomSetPtr(pWf_->pGetAtoms()), pWf_->pGetDocument(), pRefVisConfig_), info(info_), m_pWf(pWf_)
{

   pDensityEvalData = new FDensityEvalData();

   if (info.Type == PROPERTY_Density || info.Type == PROPERTY_SpinDensity ) {
      ct::FMemoryStack2
         // for storing temporary orbtial matrix.
         Mem(size_t(m_pDocument->GetWfOptions()->GetNumThreads() * m_pDocument->GetWfOptions()->GetWorkSpaceMb())<<20);
      {
         ct::TMemoryLock<double>
            pFreeMe(0, &Mem);
         if (info.Type == PROPERTY_Density) {
            // total density --- collect all occupied orbitals.
            ct::FMatrixView
               COcc;
            TArray<FOrbital*>
               RefOrbs;
            FBasisSetPtr
               pBasis;
            m_pWf->MakeOrbitalMatrix(COcc, pBasis, &RefOrbs, FFrame::ORBMAT_OccupiedOnly, Mem);
            if (RefOrbs.size() == 0) {
               IvNotify(NOTIFY_Warning, "FVolumeProperty/Density: No occupied orbitals to trace!");
               pDensityEvalData = 0;
               return;
            }
            pDensityEvalData->pOrbs = new ct::FHeapMatrix(COcc.nRows, COcc.nCols);
            Assign(*pDensityEvalData->pOrbs, COcc);

            // remember orbital basis and occupation numbers (as weights)
            pDensityEvalData->OrbWeights.resize(COcc.nCols);
            assert(COcc.nCols == RefOrbs.size());
            for (size_t iOcc = 0; iOcc < COcc.nCols; ++ iOcc)
               pDensityEvalData->OrbWeights[iOcc] = RefOrbs[iOcc]->info.fOcc;
            pDensityEvalData->pOrbBasis = pBasis;
         } else if (info.Type == PROPERTY_SpinDensity) {
            FWfType WfType = m_pWf->GetWfType();
            if (WfType == wfi::WFTYPE_Uhf || WfType == wfi::WFTYPE_Rhf) {
               // collect (pure) alpha and (pure) beta orbitals, separately.
               ct::FMatrixView
                  COccA, COccB;
               TArray<FOrbital*>
                  RefOrbsA, RefOrbsB;
               FBasisSetPtr
                  pBasisA, pBasisB;
               m_pWf->MakeOrbitalMatrix(COccA, pBasisA, &RefOrbsA, FFrame::ORBMAT_AlphaOnly, Mem);
               m_pWf->MakeOrbitalMatrix(COccB, pBasisB, &RefOrbsB, FFrame::ORBMAT_BetaOnly, Mem);
               if (RefOrbsA.size() == 0 && RefOrbsB.size() == 0) {
                  IvNotify(NOTIFY_Warning, "FVolumeProperty/SpinDensity: No spin-carrying orbitals to trace!");
                  pDensityEvalData = 0;
                  return;
               }
               if (RefOrbsB.size() == 0)
                  pBasisB = pBasisA;
               if (RefOrbsA.size() == 0)
                  pBasisA = pBasisB;

               if (pBasisA.get() != pBasisB.get()) {
                  IvNotify(NOTIFY_Error, "Internal inconsistency in spin densities: Alpha and Beta orbitals expanded over different basis sets.");
                  pDensityEvalData = 0;
                  return;
               }
               assert(COccA.nCols == RefOrbsA.size());
               assert(COccB.nCols == RefOrbsB.size());

               pDensityEvalData->pOrbBasis = pBasisA;
               pDensityEvalData->pOrbs = new ct::FHeapMatrix(pBasisA->nFn(), COccA.nCols + COccB.nCols);
               Assign(Select(*pDensityEvalData->pOrbs, 0, 0, COccA.nRows, COccA.nCols), COccA);
               Assign(Select(*pDensityEvalData->pOrbs, 0, COccA.nCols, COccB.nRows, COccB.nCols), COccB);

               pDensityEvalData->OrbWeights.resize(COccA.nCols + COccB.nCols);
               for (size_t iOccA = 0; iOccA < COccA.nCols; ++ iOccA) {
                  assert(RefOrbsA[iOccA]->info.fOcc == 1.);
                  pDensityEvalData->OrbWeights[iOccA] = RefOrbsA[iOccA]->info.fOcc;
               }
               for (size_t iOccB = 0; iOccB < COccB.nCols; ++ iOccB) {
                  assert(RefOrbsB[iOccB]->info.fOcc == 1.);
                  pDensityEvalData->OrbWeights[COccA.nCols + iOccB] = -1. * RefOrbsB[iOccB]->info.fOcc;
               }
            } else {
               IvNotify(NOTIFY_Error, "Spin densities only supported for RHF/RKS and UHF/UKS wave functions. Input orbitals are something else.");
               pDensityEvalData = 0;
               return;
            }
         } else {
            pDensityEvalData = 0;
         }
      }

   }
}

double FVolumeProperty::GetExactValue() const
{
   if (info.Type == PROPERTY_Density || info.Type == PROPERTY_SpinDensity ) {
      if (pDensityEvalData.get() == 0)
         // normally means that there are no orbitals.
         return 0.;
      double r = 0.;
      for (size_t iOrb = 0; iOrb < pDensityEvalData->OrbWeights.size(); ++ iOrb) {
         r += pDensityEvalData->OrbWeights[iOrb];
      }
      return r;
   } else {
      return -1.;
   }
}


QString FVolumeProperty::GetType() const
{
   if (pGlMesh != 0)
      return "P*"; // to mark which ones have data?
   return "P";
}


FVolumeProperty::~FVolumeProperty()
{}

QString FVolumeProperty::GetDesc(uint Flags) const
{
   IR_SUPPRESS_UNUSED_WARNING(Flags);
   return info.MakeDesc();
}


bool FVolumeProperty::DependsOnWaveFunction() const
{
   return true;
}

FIndexedTriangleListPtr FVolumeProperty::MakeIsoSurface(FIsoSurfaceSettings const &IsoSurfOpt_)
{
   return ::MakeIsoSurface(this, IsoSurfOpt_);
}


QString FVolumeProperty::GetVolumeType() const
{
   return QString("Property");
}


// namespace ct {
//    void MakeGridValues(double *pOut, FMatrixView Grid, FMatrixView Orb, uint GridDxOrder, FBasisSet const *pBasis, FMemoryStack &Mem);
// }



void FVolumeProperty::MakeGridValues(double *pOut_, ct::FMatrixView Grid_, uint GridDxOrder, int Role, ct::FMemoryStack &Mem_) const
{
   IR_SUPPRESS_UNUSED_WARNING(Role);
   // TODO:
   // - make and store orbital matrix in constructor (via FFrame... would lose the auto-linking, but maybe still not that bad of an idea to get things working)
   // - make a 'scale/factor' vector: containing occupation numbers (for total density) and +1/-1 for alpha/beta orbitals for spin density (closed: 0).
   // - then just compute orbitals via MakeGridValues. I think it actually does support real orbital matrices.
   // - still needs assembly of gradients, though.
   // - get something to estimate the amount of memory required on the work space stack.
   // - or maybe even make a real deferred evaluation object? (like FPropertyEvalData, but virtual and more general)
   //   May be required if we start supporting more general properties, like ESP, ELF, on-top hole density, laplacian-of-the-density, fukui functions, or similar.

   if (!(info.Type == PROPERTY_Density || info.Type == PROPERTY_SpinDensity))
      throw std::runtime_error("sorry, don't know how to evaluate current property type.");
   size_t
      nDxComp = 1;
   if (GridDxOrder == 1)
      nDxComp = 4;
   if (GridDxOrder == 2)
      nDxComp = 10;
   if (GridDxOrder > 1)
      throw std::runtime_error("sorry, 2nd derivatives not supported (yet) for density trace");

   // disable nested parallelism... we put out parallel evaluation thing already here because otherwise
   // we might need to make too large grid blocks (nGridPtTotal, which might consist of an entire grid with
   // millions of points, times the number of orbitals times 4 components, times size of double...).
   // For this reason, we should disable the additional parallelism level in MakeGridValues.
   omp_set_nested(0);

   size_t
      nGridPt_Total = Grid_.nCols;
   if (pDensityEvalData.get() == 0) {
      // no data... make zeros and return.
      memset(pOut_, 0, sizeof(*pOut_) * nDxComp * nGridPt_Total);
      return;
   } else {
      ct::FMemoryStackArray MemStacks(Mem_);
      size_t
         nGridStep = 64;
      int nGridPt_Batches = int((nGridPt_Total + nGridStep - 1)/nGridStep);
      #pragma omp parallel
      {
         omp_set_num_threads(1); // <- my understanding is that this is supposed to apply only to the inner region.
         #pragma omp for schedule(dynamic)
         for ( int iGridPt_Batch = 0; iGridPt_Batch < nGridPt_Batches; ++ iGridPt_Batch ) {

            size_t iGridPtBeg = size_t(iGridPt_Batch) * nGridStep;
            // ^- that's for OMP 2.0 which allows only 'int' variables for loop control.
            ct::FMemoryStack &Mem = MemStacks.GetStackOfThread();
            ct::TMemoryLock<double> pFreeMe1(0, &Mem);
            size_t
               nGridPt = nGridStep;
            if ( iGridPtBeg + nGridPt > nGridPt_Total )
               nGridPt = nGridPt_Total - iGridPtBeg;
            if (nGridPt == 0)
               continue;
            double
               *pOut = &pOut_[iGridPtBeg];

            ct::FMatrixView
               Grid(&Grid_(0,iGridPtBeg), Grid_.nRows, nGridPt, Grid_.nRowSt, Grid_.nColSt);
            // evaluate the orbitals supporting the current density.
            size_t
               nOrb = pDensityEvalData->pOrbs->nCols;
            ct::FStackMatrix
               OrbValues(nGridPt * nDxComp, nOrb, &Mem);
            ct::MakeGridValues(OrbValues.pData, Grid, *pDensityEvalData->pOrbs, GridDxOrder, &*pDensityEvalData->pOrbBasis, Mem);

            if (GridDxOrder == 0) {
               // compute density
               for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt)
                  pOut[iGridPt] = 0.;
               for (size_t iOrb = 0; iOrb < nOrb; ++ iOrb) {
                  double w = pDensityEvalData->OrbWeights[iOrb];
                  for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt) {
                     pOut[iGridPt] += w * sqr(OrbValues.pData[iGridPt + nGridPt * iOrb]);
                  }
               }
            } else if (GridDxOrder == 1) {
               // compute density and density gradient
      //          IvEmit("Entered density gradient: nGridPt = %1  nOrb = %2  nDxComp = %3", nGridPt, nOrb, nDxComp);
               size_t
                  nOrbSt = nGridPt * nDxComp;
               for (size_t iComp = 0; iComp < 4; ++ iComp)
                  memset(&pOut[iComp*nGridPt_Total], 0, sizeof(*pOut) * nGridPt);
               for (size_t iOrb = 0; iOrb < nOrb; ++ iOrb) {
                  double
                     w = pDensityEvalData->OrbWeights[iOrb];
                  for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt) {
                     double
                        *pi = &OrbValues.pData[iGridPt + nOrbSt * iOrb],
                        *po = &pOut[iGridPt];
                     double
                        v = pi[0*nGridPt],
                        vdx = pi[1*nGridPt],
                        vdy = pi[2*nGridPt],
                        vdz = pi[3*nGridPt];
                     po[0*nGridPt_Total] += w * v * v;
                     po[1*nGridPt_Total] += 2. * w * v * vdx;
                     po[2*nGridPt_Total] += 2. * w * v * vdy;
                     po[3*nGridPt_Total] += 2. * w * v * vdz;
                  }
               }
            } else {
               assert(0);
            }
         }
      }
   }
}



ct::FAtomSetPtr FVolumeDataSet::GetSupportAtoms(ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   // default: return all atoms. That is the right thing to do for global volume
   // data sets (e.g., total density, ESP, etc.)
   return pAtoms;
}



FSupportBox::FSupportBox(double const *pDipMom, double const *pQuadMom, double fEnlargeRms)
{
   for (size_t i = 0; i < 3; ++ i)
      vCenter[i] = pDipMom[i];
   for (size_t i = 0; i < 3; ++ i)
      for (size_t j = 0; j < 3; ++ j)
         vAxes[i][j] = pQuadMom[i + 3*j];
   // diagonalize the 2nd moment matrix to get the orientations of the main axes.
   double
      Ew[3];
   ct::Diagonalize(&Ew[0], &vAxes[0][0], 3, 3);
   for (size_t i = 0; i < 3; ++ i) {
      // initial estimate of extends: from 2nd moments and weight (will be refined from actual grid data later)
      fMaxExtend[i] = fEnlargeRms * std::sqrt(Ew[i]);
      fMinExtend[i] = -fMaxExtend[i];
   }
   if (det(*(mat3d*)&vAxes[0][0]) < 0) {
      // should keep positive orientation of axes. swap y and z to restore it.
      std::swap(vAxes[1], vAxes[2]);
      std::swap(fMinExtend[1], fMinExtend[2]);
      std::swap(fMaxExtend[1], fMaxExtend[2]);
   }
}

FSupportBox FVolumeDataSet::GetSupportBox(FSupportDataPtr pSupportData, double fAbsIsoThreshold, ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(Mem);

   mig::FDftGrid
      &DftGrid = *pSupportData->pDftGrid;
   ct::TArray<double>
      &Values = pSupportData->Values;
   size_t
      nPts = pSupportData->nPts();

   FVec3d
      vDm = FVec3d(0., 0., 0.);
   double
      vQm[3][3] = {0};
   {
      // compute 3d center of quantity
      double
         fTotalWeight = 0;
      for (size_t iPt = 0; iPt < nPts; ++iPt) {
         double w = DftGrid.Weights[iPt] * GetDataWeight(Values[iPt]);
         vDm += w * DftGrid.Points[iPt].vPos;
         fTotalWeight += w;
      }
      vDm /= fTotalWeight;
      // ...and 2nd moment
      for (size_t iPt = 0; iPt < nPts; ++iPt) {
         double w = DftGrid.Weights[iPt] * GetDataWeight(Values[iPt]);
         w /= fTotalWeight;
         FVec3d vRel = DftGrid.Points[iPt].vPos - vDm;
         for (size_t i = 0; i < 3; ++ i)
            for (size_t j = 0; j < 3; ++ j)
               vQm[i][j] += w * vRel[i] * vRel[j];
      }
   }
   FSupportBox
      Box(&vDm[0], &vQm[0][0], g_fSupportBox_EnlargeRms);

//    {
//       // compute 3d center of quantity
//       Box.vCenter = FVec3d(0., 0., 0.);
//       double
//          fTotalWeight = 0;
//       for (size_t iPt = 0; iPt < nPts; ++iPt) {
//          double w = DftGrid.Weights[iPt] * GetDataWeight(Values[iPt]);
//          Box.vCenter += w * DftGrid.Positions[iPt];
//          fTotalWeight += w;
//       }
//       Box.vCenter /= fTotalWeight;
//       // ...and 2nd moment
//       double vQm[3][3] = {0};
//       for (size_t iPt = 0; iPt < nPts; ++iPt) {
//          double w = Grid.Weights[iPt] * GetDataWeight(Values[iPt]);
//          w /= fTotalWeight;
//          FVec3d vRel = DftGrid.Positions[iPt] - Box.vCenter;
//          for (size_t i = 0; i < 3; ++ i)
//             for (size_t j = 0; j < 3; ++ j)
//                vQm[i][j] += w * vRel[i] * vRel[j];
//       }
//       // diagonalize the 2nd moment matrix to get the orientations of the main axes.
//       double
//          Ew[3];
//       Diagonalize(&Ew[0], &vQm[0][0], 3, 3);
//       for (size_t i = 0; i < 3; ++ i) {
//          // initial estimate of extends: from 2nd moments and weight (will be refined from actual grid data later)
//          Box.fMaxExtend[i] = s_fEnlargeRms * std::sqrt(Ew[i]);
//          Box.fMinExtend[i] = -Box.fMaxExtend[i];
//          for (size_t j = 0; j < 3; ++ j)
//             Box.vAxes[i][j] = vQm[i][j];
//       }
//       if (det(*(mat3d*)vQm[0][0]) < 0) {
//          // should keep positive orientation of axes. swap y and z to restore it.
//          std::swap(Box.vAxes[1], Box.vAxes[2]);
//          std::swap(Box.fMinExtend[1], Box.fMinExtend[2]);
//          std::swap(Box.fMaxExtend[1], Box.fMaxExtend[2]);
//       }
//    }

   // fix extends such that all data values on the grid are smaller than the target iso value
   for (size_t iAxis = 0; iAxis < 3; ++ iAxis) {
      double fMin = 1e99, fMax = -1e99;
      double fRef = ct::Dot(&Box.vCenter[0], &Box.vAxes[iAxis][0], 3);
      for (size_t iPt = 0; iPt < nPts; ++ iPt) {
         if (std::abs(Values[iPt]) > .9 * fAbsIsoThreshold) {
            double fCoord = ct::Dot(&DftGrid.Positions[iPt][0], &Box.vAxes[iAxis][0], 3) - fRef;
            if (fCoord < fMin) fMin = fCoord;
            if (fCoord > fMax) fMax = fCoord;
         }
      }
      // add a few percent of buffer zone.
      double
         fBufferSize = 0.1 * (fMax-fMin);
      Box.fMinExtend[iAxis] = fMin - fBufferSize;
      Box.fMaxExtend[iAxis] = fMax + fBufferSize;
   }

   return Box;
}


double FVolumeDataSet::GetIntegratedValue(FSupportDataPtr pSupportData) const
{
   double r = 0.;
   for (size_t i = 0; i < pSupportData->nPts(); ++ i)
      r += pSupportData->Values[i] * pSupportData->pDftGrid->Weights[i];
   return r;
}


double FVolumeDataSet::GetIntegratedWeight(FSupportDataPtr pSupportData, bool Signed) const
{
   double r = 0.;
   for (size_t i = 0; i < pSupportData->nPts(); ++ i) {
      double f = GetDataWeight(pSupportData->Values[i]) * pSupportData->pDftGrid->Weights[i];
      if (Signed && pSupportData->Values[i] < 0) {
         f = -f;
      }
      r += f;
   }
   return r;
}


double FVolumeDataSet::GetIntegratedValueBelowThreshold(double fValueThreshold, FSupportDataPtr pSupportData) const
{
   double r = 0.;
   for (size_t i = 0; i < pSupportData->nPts(); ++ i) {
      double v = pSupportData->Values[i];
//       double v = std::abs(pSupportData->Values[i]);
      // ^- should I std::abs() this? Not sure. But whatever I do here, it needs to happen in GetAbsoluteIsoThreshold, too.
      if (v < fValueThreshold)
         r += v * pSupportData->pDftGrid->Weights[i];
   }
   return r;
}


struct FIsoValueEntry
{
   double
      fValue,
      fValueWeight,
      fValueAcc; // fValueWeight * fGridWeight

  bool operator < (FIsoValueEntry const &other) const {
     // sort by value weight---in decreasing order (i.e., large weights first).
     return other.fValueWeight < fValueWeight;
  }
};

double FVolumeDataSet::GetAbsoluteIsoThreshold(double fRelativeIsoThreshold, double fRefValue, FSupportDataPtr pSupportData, ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(Mem);

   mig::FDftGrid
      &DftGrid = *pSupportData->pDftGrid;
   TArray<double>
      &Values = pSupportData->Values;
   size_t
      nPts = pSupportData->nPts();

   if (fRefValue == 0) {
      fRefValue = GetIntegratedWeight(pSupportData);
   }

   ct::TArray<FIsoValueEntry>
      IntWt(nPts);
   for (size_t i = 0; i < nPts; ++ i) {
      IntWt[i].fValue = Values[i];
      // ^- note that data can be NEGATIVE.
      IntWt[i].fValueWeight = GetDataWeight(Values[i]);
      IntWt[i].fValueAcc = GetDataWeight(Values[i]) * DftGrid.Weights[i] / fRefValue;
   }
   // sort by weight of values
   std::stable_sort(IntWt.begin(), IntWt.end());

   // find smallest value for which
   //    sum_i w_i weight(data[i]) > fThresh.
   size_t
      iThresh;
   double
      fAcc0 = 0.;
   for (iThresh = 0; iThresh < nPts; ++ iThresh) {
      fAcc0 += IntWt[iThresh].fValueAcc;
      if (fAcc0 > std::abs(fRelativeIsoThreshold))
         break;
   }
   double
      cAbsoluteIsoThreshold;
   if (iThresh != 0)
//       cAbsoluteIsoThreshold = .5 * std::abs(IntWt[iThresh-1].fValue + IntWt[iThresh].fValue));
      cAbsoluteIsoThreshold = std::sqrt(std::abs(IntWt[iThresh-1].fValue * IntWt[iThresh].fValue));
      // ^- not sure which one is better. I guess the arithmetic mean is better if positive and negative values
      // are present, while the geometric mean is better when dealing with stuff spanning many size scales. no idea.
   else {
      IvNotify(NOTIFY_Warning, "Failed to convert relative iso-surface threshold to absolute threshold (?).");
      cAbsoluteIsoThreshold = 0.; // either no points or all values are the same!!
   }

   if (fRelativeIsoThreshold < 0)
      cAbsoluteIsoThreshold *= -1.;

   return cAbsoluteIsoThreshold;
}


size_t FVolumeProperty::GetMaxStackDataSize() const
{
   switch (info.Type) {
      case PROPERTY_Density:
      case PROPERTY_SpinDensity: {
         if (pDensityEvalData.get() == 0) {
            // nothing to trace here (happens, for example, if there are no orbitals supporting the current density type).
//             IvNotify(NOTIFY_Error, "FVolumeProperty: encountered empty density evaluation data in iso-trace. This surface will be ignored!");
            return 0;
         }
         size_t
            nOccOrb = pDensityEvalData->pOrbs->nCols,
            nAo = pDensityEvalData->pOrbBasis->nFn();
         return omp_get_max_threads() * ((sizeof(double)*4*128 + 24) * (nAo + nOccOrb) + (sizeof(double) * nAo * nOccOrb));
      }
      default: {
         IvNotify(NOTIFY_Error, "FVolumeProperty::GetMaxStackDataSize(): not implemented for current volume property type.");
         return 0;
      }
   }
}

ct::FMatrixView FSupportData::GridMatrix()
{
   // return a view of the DFT grid as 3 x nPts matrix of cartesian coordinates.
   return ct::FMatrixView(&pDftGrid->Positions[0][0], 3, nPts());
}


FSupportDataPtr FVolumeDataSet::GetSupportData(ct::FMemoryStack &Mem) const
{
   FSupportDataPtr
      pOut(new FSupportData);

   FAtomSetPtr
      pAtoms = GetSupportAtoms(Mem);
   // make a DFT grid for the atoms supporting the data set.
   pOut->pDftGrid = mig::FDftGridPtr(new mig::FDftGrid(AsRawAtoms(*pAtoms), mig::FDftGridParams(2)));

   // calculate data on the grid points.
   size_t
      nPts = pOut->pDftGrid->Points.size();
   pOut->Values.resize(nPts);
   MakeGridValues(&pOut->Values[0], pOut->GridMatrix(), 0, DATAROLE_IsoSurfaceTrace, Mem);

   assert(pOut->nPts() == nPts);
   return pOut;
}


double FVolumeDataSet::GetDataWeight(double fValue_) const {
   return std::abs(fValue_);
   // note: for some quantities we should square the data to get its weight (e.g., for orbitals).
}

QString FVolumeDataSet::GetDataUnit() const
{
   return "1/a0^3";
}

