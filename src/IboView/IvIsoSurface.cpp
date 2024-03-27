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
#include <math.h>

#include "format.h"

#include "IvMesh.h"
#include "IvIsoSurface.h"
#include "IvDocument.h" // for accessing orbitals.

#include "CtMatrix.h"
#include "CxTiming.h"
// #include "CtDftGrid_ivb.h" // FIXME: merge back with main version once done.
#include "CtDftGrid.h" // FIXME: merge back with main version once done.

#include "CxAlgebra.h" // for diagonalizing matrices.
#include "CxVec3.h"
#include "CxPodArray.h"
#include "CxPhysicalUnits.h"
#include "CxOpenMpProxy.h"

#include "IvVolumeDataSet.h" // FIXME: move this out
#include "IvOrbital.h" // FIXME: move this out


using namespace ct;

// static int s_iPrintNumberIndent = 30;
// static int s_iPrintNumberIndent = 40;
static int s_iPrintNumberIndent = 33;

static void PrintTiming(QString const &Msg, double t)
{
   IvEmit(" Time for %1%2", fmts(Msg,-s_iPrintNumberIndent), fmtf(t,8,2));
}

static void PrintNumber(QString const &Msg, double fValue, QString Info = QString())
{
   if (Info.isEmpty())
      IvEmit(" %1%2", fmts(Msg,-(s_iPrintNumberIndent+11)), fmtf(fValue,16,12));
   else
      IvEmit(" %1%2 %3", fmts(Msg,-(s_iPrintNumberIndent+11)), fmtf(fValue,16,12), Info);
}



// TODO: move to CxMesh3d?
namespace ctIsoSurf {
   // defines a not-necessarily-axis-aligned cartesian grid.
   struct FCartesianGridInfo
   {
      FVec3d
         vRef,
         dXyz[3];
      unsigned
         nPts[3];
      FVec3d vPt(unsigned iX, unsigned iY, unsigned iZ) const {
         assert(iX < nPts[0] && iY < nPts[1] && iZ < nPts[2]);
         FVec3d
            v(vRef[0] + iX * dXyz[0][0] + iY * dXyz[1][0] + iZ * dXyz[2][0],
              vRef[1] + iX * dXyz[0][1] + iY * dXyz[1][1] + iZ * dXyz[2][1],
              vRef[2] + iX * dXyz[0][2] + iY * dXyz[1][2] + iZ * dXyz[2][2]);
         return v;
      }
      FVec3d vPt(double iX, double iY, double iZ) const {
         FVec3d
            v(vRef[0] + iX * dXyz[0][0] + iY * dXyz[1][0] + iZ * dXyz[2][0],
              vRef[1] + iX * dXyz[0][1] + iY * dXyz[1][1] + iZ * dXyz[2][1],
              vRef[2] + iX * dXyz[0][2] + iY * dXyz[1][2] + iZ * dXyz[2][2]);
         return v;
      }
   };

   // Function to evaluate a (nx,ny)-shape slice of the target scalar function `f`.
   // On output, pXyLayer[ix + nx*iy] = f(pGridInfo->vPt(ix,iy,iz))  for ix \in [0,...,nx-1], iy \in [0,...,nz-1]
   typedef void (*FEvalFn)(double *pXyLayer, FCartesianGridInfo const *pGridInfo, unsigned iz, void *pEvanFnArgs_);
   FIndexedTriangleListPtr TraceIsoSurface(FCartesianGridInfo const &g, FEvalFn f, void *pEvanFnArgs_, double fIsoValue, FBaseVertex const &RefVertex);
}

FIsoSurfaceSettings::FIsoSurfaceSettings(FIsoType IsoType_, double fIsoValue_,
      double fLinearDensity_, unsigned Flags_, uint32_t dwColorPlus_, uint32_t dwColorMinus_)
   : IsoType(IsoType_), fLinearDensity(fLinearDensity_), fRelativeIsoThresholdRefValue(0.), Flags(Flags_)
{
   IsoValues.push_back(FIsoThresholdEntry(fIsoValue_, dwColorPlus_));
   if (PlusAndMinus()) {
      IsoValues.push_back(FIsoThresholdEntry(-fIsoValue_, dwColorMinus_));
   }
}



struct FDataBox
{
   FVec3d vCenter;
   FVec3d vAxes[3];
   // nPts[i] points evenly distributed from
   //      fMinExtend[i]*vAxes[i] ... fMaxExtend[i]*vAxes[i]
   double fMinExtend[3];
   double fMaxExtend[3];
   unsigned nPts[3];
   double fLinearDensity;

   void FixAxisDensity();
   void PrintInfo(std::string const &OrbitalDesc);

   size_t nPtsTotal() const { return nPts[0] * nPts[1] * nPts[2]; }

   FDataBox(FSupportBox const &SupportBox_, double fLinearDensity_);
};


FDataBox::FDataBox(FSupportBox const &SupportBox_, double fLinearDensity_)
   : fLinearDensity(fLinearDensity_)
{
   vCenter = SupportBox_.vCenter;
   for (size_t iAxis = 0; iAxis < 3; ++ iAxis) {
      vAxes[iAxis] = SupportBox_.vAxes[iAxis];
      fMinExtend[iAxis] = SupportBox_.fMinExtend[iAxis] - 2 * fLinearDensity;
      fMaxExtend[iAxis] = SupportBox_.fMaxExtend[iAxis] + 2 * fLinearDensity;
   }

   FixAxisDensity();
}


// find number of points for each axis and adjust fMinExtend/fMaxExtend
// in such away that fLinearDensity is achieved exactly.
void FDataBox::FixAxisDensity()
{
   // decide on number of points in each direction, and fix up extends
   // in such a way that they meet the linear density requirements exactly.
   for (unsigned i = 0; i < 3; ++ i) {
      // ^- why do we have matrix classes if we cannot access them as vector rows?
//       double fLengthOrig = 2 * length(Box.vAxes[i]); // 2*:   -vAxes[i] ... +vAxes[i]
      double fLengthOrig = fMaxExtend[i] - fMinExtend[i];
      nPts[i] = 1 + (unsigned)(fLengthOrig / fLinearDensity);
      if (nPts[i] < 8)
         nPts[i] = 8;
      double f = nPts[i] * fLinearDensity / fLengthOrig;
      fMaxExtend[i] *= f;
      fMinExtend[i] *= f;
   }
}

void FDataBox::PrintInfo(std::string const &/*OrbitalDesc*/)
{
   int w = 8, p = 3;
   IvEmit(" Box:  vCen = (%1,%2,%3)  nPts = (%4,%5,%6)",
      fmtf(vCenter[0],w,p), fmtf(vCenter[1],w,p), fmtf(vCenter[2],w,p),
      fmti(nPts[0],4), fmti(nPts[1],4), fmti(nPts[2],4));
   IvEmit("       vExt = (%1,%2,%3)  -> %4 points total",
      fmtf(fMaxExtend[0]-fMinExtend[0],w,p), fmtf(fMaxExtend[1]-fMinExtend[1],w,p), fmtf(fMaxExtend[2]-fMinExtend[2],w,p),
      nPtsTotal());
}


struct FGridDataRange;
struct FPropertyEvalData;
using ctIsoSurf::FCartesianGridInfo;

// Grid and orbital values for a subrange of a cartesian grid.
// idea is that we can either calculate it inline (low memory)
// or make it beforehand to share the data for multiple iso values.
struct FGridDataRange : public FIntrusivePtrDest
{
   FGridDataRange() {}
   FGridDataRange(unsigned nx, unsigned ny, unsigned iz0_, unsigned izN_, FPropertyEvalData *pData);
   ~FGridDataRange();
   unsigned
      iz0, izN;

   FMatrixView
      Values;
   double
      *pDataValues; // nx * ny * (izN-iz0) array.
private:
   FGridDataRange(FGridDataRange const &); // not implemented.
   void operator = (FGridDataRange const &); // not implemented.
};
typedef ct::TIntrusivePtr<FGridDataRange>
   FGridDataRangePtr;


struct FPropertyEvalData
{
   ct::FMemoryStack
      *pMem;
   FVolumeDataSet
      *pVolumeDataSet;
   FCartesianGridInfo
      *pCaGrid;
   FGridDataRangePtr
      pGridData;
};



FGridDataRange::FGridDataRange(unsigned nx, unsigned ny, unsigned iz0_, unsigned izN_, FPropertyEvalData *pData)
   : iz0(iz0_), izN(izN_)
{
   unsigned N = nx * ny * (izN-iz0);
   pDataValues = (double*)::malloc(N * sizeof(pDataValues[0]));
   Values = FMatrixView(pDataValues, nx * ny, (izN-iz0));

   // make the orbital values in slices. one for each z.
   for (unsigned iz = iz0; iz != izN; ++ iz) {
      // store the x/y/z coordinates of the grid points
      FStackMatrix
         Grid(3, nx * ny, pData->pMem);
      for (unsigned ix = 0; ix < nx; ++ ix)
         for (unsigned iy = 0; iy < ny; ++ iy) {
            unsigned ixy = ix + nx * iy;
            FVec3d
               v = pData->pCaGrid->vPt(ix,iy,iz);
            Grid(0, ixy) = v[0];
            Grid(1, ixy) = v[1];
            Grid(2, ixy) = v[2];
         }
      // evaluate data on the grid points for the current iz.
      pData->pVolumeDataSet->MakeGridValues(&Values(0,iz-iz0), Grid, 0, DATAROLE_IsoSurfaceTrace, *pData->pMem);
   }
}

FGridDataRange::~FGridDataRange()
{
   ::free(pDataValues);
   pDataValues = 0;
};



static void EvalVolumeData2(double *pVxy, FCartesianGridInfo const *pGridInfo, unsigned iz, void *pData_)
{
   FPropertyEvalData
      *pData = reinterpret_cast<FPropertyEvalData*>(pData_);
   unsigned
      nx = pGridInfo->nPts[0],
      ny = pGridInfo->nPts[1];

   FGridDataRangePtr
      pGridData = pData->pGridData;
   if (pGridData.get() == 0)
      // no precalculated data. make data for current slice.
      pGridData = FGridDataRangePtr(new FGridDataRange(nx, ny, iz, iz+1, pData));

   // copy to target location.
   for (unsigned ix = 0; ix < nx; ++ ix)
      for (unsigned iy = 0; iy < ny; ++ iy)
         pVxy[ix + nx*iy] = pGridData->Values(ix + nx * iy, iz - pGridData->iz0);
}




void FixVolumeDataNormals(FIndexedTriangleList *p, FPropertyEvalData *pData)
{
   size_t
      nPt = p->Vertices.size();
   FStackMatrix
      Grid(3, nPt, pData->pMem),
      Values(nPt, 4, pData->pMem);
   for (unsigned iPt = 0; iPt < nPt; ++ iPt)
      for (unsigned ixyz = 0; ixyz < 3; ++ ixyz)
         Grid(ixyz, iPt) = p->Vertices[iPt].vPos[ixyz];
   Values.Clear(); // FIXME: remove this.

   // calculate actual gradient of the volume data on the grid points
   pData->pVolumeDataSet->MakeGridValues(Values.pData, Grid, 1, DATAROLE_IsoSurfaceTrace, *pData->pMem);

   // check if we are on a positive or negative side of the iso-surface.
   // Surface construction is such that the triangles' front sides
   // point into the direction of more positive function values.
   double fSum = 0.;
   for (unsigned iPt = 0; iPt < nPt; ++ iPt)
      fSum += Values(iPt,0);
   // ^-- TODO: can't we just use the actual iso value to decide that?

   // copy back gradient, and normalize.
   for (unsigned iPt = 0; iPt < nPt; ++ iPt) {
      vec3f &vNorm = p->Vertices[iPt].vNorm;
      for (unsigned ixyz = 0; ixyz < 3; ++ ixyz)
         vNorm[ixyz] = Values(iPt, 1+ixyz);
      // -: we want them to point to the outside.
//       std::cout << fmt::format("iPt = {:5}  nx = {:10.5e} ny = {:10.5e}  nz = {:10.5e}", iPt, vNorm[0], vNorm[1], vNorm[2]) << "\n";
      if (fSum >= 0)
         vNorm *= -1;
      vNorm /= vmath::length(vNorm);
   }

   if (fSum >= 0) {
      // invert the triangles; they are currently showing backside out.
      for (size_t i = 0; i < p->Triangles.size(); ++ i)
         std::swap(p->Triangles[i][1], p->Triangles[i][2]);
   }

   // TODO: at this point it might be a good idea to also compute and store the property values
   // with DATAROLE_ValuesOnIsoSurface (if different from DATAROLE_IsoSurfaceTrace).
}


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}



FIndexedTriangleListPtr MakeIsoSurface(FVolumeDataSet const *pVolumeDataSet_, FIsoSurfaceSettings const &Options)
{
   IvEmit("\n");
   IvNotify(NOTIFY_StartWork, IvFmt("!Tracing %1 %2 [density: (%3/A)^3]", pVolumeDataSet_->GetVolumeType(),
      pVolumeDataSet_->GetDesc(), fmtf(Options.fLinearDensity,0,1)));


   ct::FMemoryStack2
//       Mem(200000000); // ~200 mb
      Mem(20000000 + pVolumeDataSet_->GetMaxStackDataSize());
      // ^- ~20 mb + whatever the volume thinks it needs for evaluating grid values.

   ct::FTimer tStart;

   FVolumeDataSet
      *pVolumeDataSet = const_cast<FVolumeDataSet*>(pVolumeDataSet_);
   // ^- yes.. I know. It is 1:39 AM... just want to get this working.

   FPropertyEvalData
      data;
   data.pMem = &Mem;
   data.pVolumeDataSet = pVolumeDataSet;


   ct::FTimer tInitialDataBox;
   FSupportDataPtr
      pSupportData = pVolumeDataSet->GetSupportData(Mem);

   double
      fPosOrNeg = pVolumeDataSet->GetIntegratedWeight(pSupportData, true);
   PrintNumber("Positive/negative balance", fPosOrNeg);

   bool
      // try to put the bigger lobe on the positive side.
      FlipIso = Options.AllowFlip() && (fPosOrNeg < 0.);

   // convert input iso tresholds to target absolute iso thresholds, if requested.
   FIsoValueList
      AbsoluteIsoValues(Options.IsoValues);
   if (Options.IsoType == ISOTYPE_Relative) {
      double fRefValue = Options.fRelativeIsoThresholdRefValue;
      if (fRefValue == 0.) {
         fRefValue = pVolumeDataSet->GetIntegratedWeight(pSupportData);
         PrintNumber(IvFmt("Integrated %1 weight", pVolumeDataSet->GetVolumeType().toLower()), fRefValue);
      }

      // convert relative iso thresholds to absolute iso thresholds.
      for (size_t iThresh = 0; iThresh < AbsoluteIsoValues.size(); ++ iThresh) {
         double
            fRelThresh = Options.IsoValues[iThresh].fIsoValue,
            fAbsThresh = pVolumeDataSet->GetAbsoluteIsoThreshold(
               (FlipIso? (-fRelThresh) : fRelThresh), fRefValue, pSupportData, Mem);
         AbsoluteIsoValues[iThresh].fIsoValue = fAbsThresh;
         if (fRelThresh > 0)
            PrintNumber(IvFmt("Iso-surface value/[%1] (%2% weight)", pVolumeDataSet->GetDataUnit(), int(100.*fRelThresh)), fAbsThresh, (FlipIso?" (flipped)":""));
      }
   } else {
      // print requested absolute iso-surface threshold values
      for (size_t iThresh = 0; iThresh < AbsoluteIsoValues.size(); ++ iThresh) {
         double
            fAbsThresh = AbsoluteIsoValues[iThresh].fIsoValue;
         PrintNumber(IvFmt("Iso-surface value/[%1] (fixed)", pVolumeDataSet->GetDataUnit()), fAbsThresh, (FlipIso?" (flipped)":""));
// For scale: an iso-value 0.001 e/a0^3 has been used to *define* atomic van der Waals-type radii (Rahm, Hoffmann, Ashcroft: "Atomic and Ionic Radii of Elements 1-96", Chem Eur J 22, 14625 (2016); 10.1002/chem.201602949)
      }
   }


   // construct cartesian bounding box of data. Should contain all data values in the given iso-value range.
   double
      fLinearDensity = 1./(ToAng*Options.fLinearDensity);
   double
      fMinAbsoluteIsoThreshold = 1e99;
   for (size_t iThresh = 0; iThresh < AbsoluteIsoValues.size(); ++ iThresh)
      fMinAbsoluteIsoThreshold = std::min(fMinAbsoluteIsoThreshold, std::abs(AbsoluteIsoValues[iThresh].fIsoValue));
   FDataBox
      ob(pVolumeDataSet->GetSupportBox(pSupportData, fMinAbsoluteIsoThreshold, Mem), fLinearDensity);
   ob.PrintInfo(pVolumeDataSet->GetDesc().toStdString());

   pSupportData = 0; // free memory. no longer required.


   FCartesianGridInfo
      CaGrid;
   CaGrid.vRef = ob.vCenter;  // should be ob.vCenter - \sum[i] ob.vAxes[i] in the end.
   for (unsigned i = 0; i < 3; ++ i) {
      double fExtend = ob.fMaxExtend[i] - ob.fMinExtend[i];
      CaGrid.dXyz[i] = (fExtend/(ob.nPts[i]-1)) * ob.vAxes[i];
//       CaGrid.vRef -= (.5*ob.nPts[i]) * CaGrid.dXyz[i];
      CaGrid.vRef += ob.fMinExtend[i] * ob.vAxes[i];
      CaGrid.nPts[i] = ob.nPts[i];
   }

   data.pCaGrid = &CaGrid;

   size_t
      nTrianglesTotal = 0;
   PrintTiming("Iso-surface setup", tStart);
   bool
      PrecalcValues = true;
   if (PrecalcValues) {
      // precalculate orbital values on the grid. This way they can be shared
      // between the plus and the minus surfaces.
      // (otherwise they'd be calculated inline. this required much less memory,
      //  but requires re-calculating the values for each surface.)
      ct::FTimer tOrbOnGrid;
      data.pGridData = FGridDataRangePtr(new FGridDataRange(CaGrid.nPts[0], CaGrid.nPts[1], 0, CaGrid.nPts[2], &data));
      PrintTiming("Data on grid (trace)", tOrbOnGrid);
   }

   FIndexedTriangleListPtr
      pTriListCombined,
      pTriList;

   for (size_t iIsoValue = 0; iIsoValue < AbsoluteIsoValues.size(); ++ iIsoValue) {
      FIsoThresholdEntry const
         &ite = AbsoluteIsoValues[iIsoValue];
      FBaseVertex RefVertex;
      RefVertex.dwColor = ite.dwColor;
      RefVertex.iRole = iIsoValue;

      ct::FTimer tSurfaceTrace;

      pTriList = ctIsoSurf::TraceIsoSurface(CaGrid, EvalVolumeData2, &data, ite.fIsoValue, RefVertex);
      IvEmit("    -> %1 vertices and %2 triangles.", pTriList->Vertices.size(), pTriList->Triangles.size());
      nTrianglesTotal += pTriList->Triangles.size();

      ct::FTimer tTriangleConv;
      FixVolumeDataNormals(&*pTriList, &data);
      PrintTiming("Iso-surface conv/normals", tTriangleConv);
      if (iIsoValue == 0)
         pTriListCombined = pTriList;
      else
         pTriListCombined->Append(*pTriList);
   }
   PrintTiming("Iso-surface total", tStart);

//    IvNotify(NOTIFY_FinishWork, IvFmt("Ready"));
   IvNotify(NOTIFY_FinishWork, IvFmt("Ready [%1k grid points and %2 triangles finished after %3 sec].", ob.nPtsTotal()/1000, nTrianglesTotal, QString("%1").arg((double)tStart,0,'f',2)));
   return pTriListCombined;
}



namespace ctIsoSurf {
   // really just a bunch of triangles...
   struct FSurface
   {
      typedef TArray<FTriangle1>
         FTriangleArray;
      FTriangleArray
         Triangles;
      void AddTriangle(FVec3d const &v0, FVec3d const &v1, FVec3d const &v2);

      // note: this does not generate normals.
      FIndexedTriangleListPtr ConvertToIndexed(FBaseVertex const &RefVertex);
   };

   void FSurface::AddTriangle(FVec3d const &v0, FVec3d const &v1, FVec3d const &v2)
   {
      Triangles.push_back(FTriangle1(v0,v1,v2));
   }

   FIndexedTriangleListPtr FSurface::ConvertToIndexed(FBaseVertex const &RefVertex)
   {
      ct::TArray<FVec3d> Pos;
      ct::TArray<unsigned> Indices;

      MakeIndexedTriangles(Pos, Indices, this->Triangles);
      return FIndexedTriangleListPtr(new FIndexedTriangleList(&Pos[0], &Pos[0], Pos.size(), &Indices[0], Indices.size(), RefVertex));
   }


   enum FVertexType {
      VERTEX_None, // "no iso-surface intersection here (vertex does not exist)"
      VERTEX_LeftSide,
      VERTEX_RightSide,
   };

   struct FIsoVertex
   {
      FVec3d
         vPos;
      FVertexType
         Type;
   };

   // one two-dimensional layer of data (i.e., nx*ny for a given iz)
   template<class T>
   struct TLayer : public TArray<T>
   {
      typedef typename ct::TArray<T>
         FBase;

      TLayer() {};

      TLayer(size_t nx_, size_t ny_)  {
         resize(nx_, ny_);
      }

      void resize(size_t nx_, size_t ny_) {
         nx = nx_;
         ny = ny_;
         FBase::resize(nx*ny);
      }

      T &operator() (size_t ix, size_t iy) {
         assert(ix < nx && iy < ny);
         return (*this)[ix + nx*iy];
      }

      size_t
         nx, ny;

      void swap(TLayer &other) {
         FBase::swap(*(FBase*)&other);
         std::swap(this->nx, other.nx);
         std::swap(this->ny, other.ny);
         // ^- not really required because we anyway use it only with equal
         //    dimensions.
      }
   };
   typedef TLayer<double>
      FDataLayer;
   typedef TLayer<FIsoVertex>
      FVertexLayer;

   // one set of iso surface/cube edge intersections
   struct FIsoSlice
   {
      FIsoSlice(size_t nx, size_t ny) {
         for (size_t ixyz = 0; ixyz < 3; ++ ixyz)
            VerticesXyz[ixyz].resize(nx, ny);
      }

      void swap(FIsoSlice &other) {
         for (size_t ixyz = 0; ixyz < 3; ++ ixyz)
            VerticesXyz[ixyz].swap(other.VerticesXyz[ixyz]);
      }

      FVertexLayer
         // vertices for iso-surface intersections along the three cartesian
         // directions.
         // Note:
         //   [0] has (nx-1) x ny entries
         //   [1] has nx * (ny-1) entries
         //   [2] has nx*ny entries.
         VerticesXyz[3];
   };


   static void FindVerticesForIsoIntersectionsOnAxis(FVertexLayer &VertexLayer,
      FCartesianGridInfo const &g,  size_t ld, size_t nx, size_t ny,
      double *pLayerA, double *pLayerB, size_t iz, size_t ixyz, double fIsoValue)
   {
      for (size_t iy = 0; iy < ny; ++ iy)
         for (size_t ix = 0; ix < nx; ++ ix) {
            double
               dA = pLayerA[ix + ld*iy] - fIsoValue,
               dB = pLayerB[ix + ld*iy] - fIsoValue;
            bool
               Left = (dA >= 0. && dB < 0),
               Right = (dA < 0 && dB >= 0.);
            FIsoVertex
               *v = &VertexLayer(ix, iy);
            if (Left || Right) {
               v->vPos = g.vPt((unsigned)ix, (unsigned)iy, (unsigned)iz);
               v->vPos += dA/(dA-dB) * g.dXyz[ixyz];
               v->Type = Left? VERTEX_LeftSide : VERTEX_RightSide;
            } else {
               v->Type = VERTEX_None;
            }
         }
   }

   static void FindVerticesForIsoIntersections(FIsoSlice &Slice, FCartesianGridInfo const &g, FDataLayer *pLayerZ0, FDataLayer *pLayerZ1, double fIsoValue, size_t iz0)
   {
      size_t ld = pLayerZ0->nx;

      // x/x+1 direction (at z = 0)
      FindVerticesForIsoIntersectionsOnAxis(Slice.VerticesXyz[0], g, ld, g.nPts[0]-1, g.nPts[1], &(*pLayerZ0)[0], &(*pLayerZ0)[1], iz0, 0, fIsoValue);
      // y/y+1 direction (at z = 0)
      FindVerticesForIsoIntersectionsOnAxis(Slice.VerticesXyz[1], g, ld, g.nPts[0], g.nPts[1]-1, &(*pLayerZ0)[0], &(*pLayerZ0)[ld], iz0, 1, fIsoValue);
      // z/z+1 direction
      if (pLayerZ1) {
         FindVerticesForIsoIntersectionsOnAxis(Slice.VerticesXyz[2], g, ld, g.nPts[0], g.nPts[1], &(*pLayerZ0)[0], &(*pLayerZ1)[0], iz0, 2, fIsoValue);
      } else {
         for (size_t iy = 0; iy < g.nPts[1]; ++ iy)
            for (size_t ix = 0; ix < g.nPts[0]; ++ ix)
               Slice.VerticesXyz[2](ix,iy).Type = VERTEX_None;
      }
   }

   // Format is as follows:
   //         {iDirection, ix,iy,iz}
   // where iDirection gives the direction along which the edge is oriented,
   // and ix/iy/iz specify the *TWO* other edge coordinates, which are *NOT* determined
   // by iDirection. That is, e.g., the edge specification
   //        {#y, ix,iy,iz}
   // translates to
   //        (ix,0,iz)--(ix,1,iz)
   // with iy being ignored.
   //
   // Note: I use the same edge order and numbering as GTS, from which the algorithm
   // idea used here is recycled (see GTS's isocube.fig).
   static unsigned char const iEdgeCoords[12][4] = {
      {0,0,0,0}, {0,0,0,1}, {0,0,1,1}, {0,0,1,0},    // <- the 4 edges along the x direction
      {1,0,0,0}, {1,0,0,1}, {1,1,0,1}, {1,1,0,0},    // <- the 4 edges along the y direction
      {2,0,0,0}, {2,1,0,0}, {2,1,1,0}, {2,0,1,0}     // <- the 4 edges along the z direction
   };

   // For each edge, this gives the indices of the other three edges which lie
   // on the same cube face. First index decides about which of the two faces
   // lying on the same edge are meant.
   static unsigned char const iEdgeLinks[2][12][3] = {
      {{9,1,8}, {6,2,5}, {10,3,11}, {7,0,4}, {3,7,0}, {11,4,8}, {2,5,1}, {10,6,9}, {5,11,4}, {1,8,0}, {6,9,7}, {2,10,3}},
      {{4,3,7}, {8,0,9}, {5,1,6}, {11,2,10}, {8,5,11}, {1,6,2}, {9,7,10}, {0,4,3}, {0,9,1}, {7,10,6}, {3,11,2}, {4,8,5}}
   };


   // make iso surface fragments (for already calculated vertices) between two layers of data.
   // I.e., make triangle indices.
   static void MakeIsoSurfaceSlice(FSurface &SurfaceOut, FIsoSlice *pSliceZ0, FIsoSlice *pSliceZ1)
   {
      size_t
         nx = pSliceZ0->VerticesXyz[0].nx,
         ny = pSliceZ0->VerticesXyz[0].ny;
      FIsoSlice
         *(pSlicesZ01[2]) = {pSliceZ0, pSliceZ1};
      FIsoVertex
         EdgeVertices[12],
         FaceVertices[12]; // vertices at edge/surface intersections

      // iterate over the (nx-1) x (ny-1) cubes in between layers at z0 and z1.
      for (size_t iy = 0; iy < ny - 1; ++ iy)
         for (size_t ix = 0; ix < nx - 1; ++ ix) {
            // copy all vertices of cube-edge/iso-surface intersections of current cube.
            for (size_t iEdge = 0; iEdge < 12; ++ iEdge) {
               unsigned char const
                  *ec = &iEdgeCoords[iEdge][0];
               EdgeVertices[iEdge] = pSlicesZ01[ec[3]]->VerticesXyz[ec[0]](ix + ec[1], iy + ec[2]);
            }
            // marker for whether the edge intersection has already been handled
            // as part of another edge trace
            bool EdgeHandled[12] = {0};

            for (size_t iEdge0 = 0; iEdge0 < 12; ++ iEdge0) {
               if (EdgeHandled[iEdge0] || EdgeVertices[iEdge0].Type == VERTEX_None)
                  continue;
               // current edge intersects with the iso surface and has not
               // already been part of a face starting at a previous edge.

               // -> make a new face and find all other edges which belong to
               // the current face.
               size_t
                  // number intersections of current cube's edges with iso surface.
                  // corresponding edges are stored in v[].
                  nEdges = 0,
                  iEdge1 = iEdge0;
               while (!EdgeHandled[iEdge1] && EdgeVertices[iEdge1].Type != VERTEX_None) {
                  FIsoVertex
                     &v1 = EdgeVertices[iEdge1];
                  FaceVertices[nEdges] = v1;
                  EdgeHandled[iEdge1] = true; // don't visit this edge again.
                  nEdges += 1;
                  size_t
                     iSide = (v1.Type == VERTEX_LeftSide)? 0 : 1;
                  // look for another edge on the current face which intersects the surface
                  unsigned char const
                     *iEdge1Links = &iEdgeLinks[iSide][iEdge1][0];
                  for (size_t iLink = 0; iLink < 3; ++ iLink) {
                     iEdge1 = iEdge1Links[iLink];
                     if (EdgeVertices[iEdge1].Type != VERTEX_None)
                        break;
                  }
               }

               if (nEdges >= 3)
                  for (size_t iTri = 0; iTri < nEdges - 2; ++ iTri)
                     SurfaceOut.AddTriangle(FaceVertices[0].vPos, FaceVertices[iTri+1].vPos, FaceVertices[iTri+2].vPos);
            }
         }
   }

   FIndexedTriangleListPtr TraceIsoSurface(FCartesianGridInfo const &g, FEvalFn f, void *pEvanFnArgs_, double fIsoValue, FBaseVertex const &RefVertex)
   {
      FSurface
         SurfaceOut;

      size_t
         nx = g.nPts[0],
         ny = g.nPts[1],
         nz = g.nPts[2];
      // we keep two layers (iz and iz+1) in memory at each time. This here
      // stores the data of them. Format: [ix + nx*iy]
      FDataLayer
         LayerA(nx, ny), // last data layer
         LayerB(nx, ny); // current data layer
      FIsoSlice
         SliceA(nx, ny),
         SliceB(nx, ny);
         // ^- if we generalize this for multiple simultaneous surface traces, then
         // these things need to be supplied once for each iso value.

      // initialize first layer.
      f(&LayerA[0], &g, 0, pEvanFnArgs_);

      for (size_t iz = 0; iz < nz; ++ iz) {
         if (iz != nz-1) {
            // evaluate function at next z
            f(&LayerB[0], &g, iz+1, pEvanFnArgs_);
            // find intersections of cube edges:
            //   - edges along x axis for current iz
            //   - edges along y axis for current iz
            //   - edges along z axis for z between (iz and iz+1)
            FindVerticesForIsoIntersections(SliceB, g, &LayerA, &LayerB, fIsoValue, iz);
         } else {
            // this is the last slice: do not evaluate f anymore,
            // but still make cube edge/iso surface intersections along x and y
            // directions.
            FindVerticesForIsoIntersections(SliceB, g, &LayerA, 0, fIsoValue, iz);
         }

         if (iz != 0)
            // Construct faces for iso surface elements between LayersA and LayerB
            // With both SliceA and SliceB we now have all edge intersections for
            // the 1x1x1 cubes with start coordinate index (ix,iy,iz-1),
            // where ix in [0..nx-1], iy in [0..ny-1].
            MakeIsoSurfaceSlice(SurfaceOut, &SliceA, &SliceB);

         LayerA.swap(LayerB);
         SliceA.swap(SliceB);
      }
      return SurfaceOut.ConvertToIndexed(RefVertex);
   }
} // namespace ctIsoSurf
