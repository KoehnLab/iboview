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

#ifndef IV_VOLUME_DATASET_H
#define IV_VOLUME_DATASET_H

#include "Iv.h"

#include <map>
#include <set>
#include <vector>
#include <list>

#include <QString>
#include <QObject>
// #include <QAbstractTableModel>
// #include <QAction>
// #include <QExplicitlySharedDataPointer>
// #include <QSharedData>


#include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "CxPodArray.h"

// #include "CxColor.h"
#include "IvMesh.h"
#include "IvGl.h"
#include "IvIsoSurface.h"

#include "IvDataOptions.h"
#include "IvDataSet.h"

// #include "IvAnalysis.h"
#include "CtMatrix.h"
// #include "IvLog.h"
// #include "IvTables.h"
// #include "IvIrc.h"
#include "CxVec3.h"
// #include "CtWfi.h"


// #include "CtDftGrid_ivb.h" // FIXME: merge back with main version once done.
#include "CtDftGrid.h"


class FFrame;
class FView3d;
class FDocument;



using ct::FAtomSetPtr;
using ct::FBasisSetPtr;
using ct::TArray;
using ct::FIntrusivePtrDest;


struct FVolumeDataSet;


// contains the actual parametric details of the visual configuration of a given volume;
// depends on concrete visualization type.
// Note: we might want to turn these things into QObjects, to allow script access.
struct FVolumeVisualConfigDetails : public FIntrusivePtrDest
{
   FVolumeVisualConfigDetails() {};
   virtual ~FVolumeVisualConfigDetails();
   virtual void AssignTo(FIsoSurfaceSettings &IsoSurfOpt_) const = 0;
   virtual uint32_t GetBaseColor() const = 0;
private:
   void operator = (FVolumeVisualConfigDetails const &); // not implemented
   FVolumeVisualConfigDetails(FVolumeVisualConfigDetails const &); // not implemented
};
typedef ct::TIntrusivePtr<FVolumeVisualConfigDetails>
   FVolumeVisualConfigDetailsPtr;

// simple iso surface: There is only one iso threshold T, and we paint surfaces
// at +/- T, with individual colors (cIsoPlus/cIsoMinus).
// This is the default choice for orbitals.
struct FTwoPhaseIsoSurfaceConfig : public FVolumeVisualConfigDetails
{
   FIsoType
      iIsoType;
   float
      fIsoValue;
   uint32_t
      cIsoPlus, cIsoMinus;
   FTwoPhaseIsoSurfaceConfig() {}
   FTwoPhaseIsoSurfaceConfig(FIsoType IsoType_, float fIsoValue_, uint32_t cIsoPlus_, uint32_t cIsoMinus_)
      : iIsoType(IsoType_), fIsoValue(fIsoValue_), cIsoPlus(cIsoPlus_), cIsoMinus(cIsoMinus_)
   {}
   void AssignTo(FIsoSurfaceSettings &IsoSurfOpt_) const; // override
   uint32_t GetBaseColor() const; // override
};

FVolumeVisualConfigDetailsPtr MakeDefaultTwoPhaseIsoSurfaceConfig(FDocument *pDocument, double fIsoValue_);

// iso-surface supporting mutiple iso-values (not necessary +/-) with arbitrary colors.
struct FMultiPhaseIsoSurfaceConfig : public FVolumeVisualConfigDetails
{
   FIsoType
      iIsoType;
   FIsoValueList
      IsoValues;
   FMultiPhaseIsoSurfaceConfig() {}
   FMultiPhaseIsoSurfaceConfig(FIsoType IsoType_, FIsoValueList IsoValues_) : iIsoType(IsoType_), IsoValues(IsoValues_) {}
   void AssignTo(FIsoSurfaceSettings &IsoSurfOpt_) const; // override
   uint32_t GetBaseColor() const; // override
};

// describes the configuration of the visual representation of an orbital or other
// volume-dependent quantity (i.e., function of 3d space).
// So far we only support iso-surfaces of the quantity itself.
// Rendering of color-coded quantities on the iso-surfaces of other quantities
// (e.g., ESP on vdW surface) will likely also be made available.
//
// Note: These visual-configuration objects are *shared* between associated
// datasets in different frames.
struct FVolumeVisualConfig : public FIntrusivePtrDest
{
   bool DetailsAssigned() const { return pDetails.get() != 0; }

   FVolumeVisualConfigDetailsPtr
      // initially 0. Becomes non-0 when the properties are first set.
      pDetails;

   FVolumeVisualConfig();
   explicit FVolumeVisualConfig(FVolumeVisualConfigDetailsPtr pDetails_);
   ~FVolumeVisualConfig();
public:
   typedef std::set<FVolumeDataSet*>
      FVolumeChain;
   FVolumeChain
      // list of all orbitals which are linked with this configuration object
      // (note: these are weak refs; the orbitals own the config object (via
      // ref), not the other way around!)
      LinkedVolumes;
   void Link(FVolumeDataSet *pVolume);
   void Unlink(FVolumeDataSet *pVolume);

   enum FUpdateFlags {
      UPDATE_InvalidateIsoSettings = 0x01,
      UPDATE_InvalidateColors = 0x02,
      UPDATE_Rebuild = 0x04 // <- for this pView3d must be supplied!
   };

   // Apply given update/change to all data sets which share this object.
   // argument: bit field of UPDATE_* flags.
   void UpdateLinkedRepresentations(uint32_t Flags, FView3d *pView3d=0);
private:
   void operator = (FVolumeVisualConfig const &); // not implemented
   FVolumeVisualConfig(FVolumeVisualConfig const &); // not implemented
   // ^- note: it is not the actual config part of which copying is the problem...
   //    the linked reference chain is. It might be useful to split this class into
   //    two parts: one "shared dataset properties" class which represents the links,
   //    and one for the actual visual configuration. This may become useful once other
   //    kinds of data sets besides orbitals and geometries are introduced.
};
typedef ct::TIntrusivePtr<FVolumeVisualConfig>
   FVolumeVisualConfigPtr;


// a cartesian box which is supposed to contain all data of a given volume, for a given minimum iso threshold.
struct FSupportBox {
   FVec3d vCenter;
   FVec3d vAxes[3];
   double fMinExtend[3];
   double fMaxExtend[3];

   FSupportBox() {}
   FSupportBox(double const *pDipMom, double const *pQuadMom, double fEnlargeRms);
};

// support structure for computing support box and iso-threshold data.
struct FSupportData : public FIntrusivePtrDest
{
   mig::FDftGridPtr
      // a numerical integration grid with supposedly reasonably complete
      // coverage of the current volume quantity.
      pDftGrid;
   TArray<double>
      // data on pDftGrid's grid points.
      Values;
   size_t nPts() const {
      assert(pDftGrid->Points.size() == Values.size());
      return pDftGrid->Points.size();
   }
   // return a view of the DFT grid as 3 x nPts matrix of cartesian coordinates.
   ct::FMatrixView GridMatrix();
};
typedef ct::TIntrusivePtr<FSupportData>
   FSupportDataPtr;

// we might want to support volume properties which are visualized on the iso-surface
// of OTHER properties (e.g., electrostatic potential (ESP) on the vdW surface).
// For this reason, we introduce two distinct roles for evaluating grid values:
// - DATAROLE_IsoSurfaceTrace: make values on the grid which determines the geometry of iso-surface
// - DATAROLE_ValuesOnIsoSurface: make final values to visualize. This will be called for points ON the
//   previously determined iso-surface.
enum FGridDataRole {
   DATAROLE_IsoSurfaceTrace, // compute data for the sake of tracing an iso-surface
   DATAROLE_ValuesOnIsoSurface // compute final data to visualize, on the previously determined iso-surface
};

struct FVolumeDataSet : public FDataSet
{
   typedef FDataSet
      FBase;
   explicit FVolumeDataSet(QString const &Desc_, ct::FAtomSetPtr pAtoms_, FDocument *pDocument_, FVolumeVisualConfigPtr pRefVisConfig_);
   ~FVolumeDataSet();

   FVolumeVisualConfigPtr
      pVisConfig;
   void LinkVisualConfig(FVolumeVisualConfigPtr p);

   // role:
   //  - 0: return values for establishing the actual iso-surface
   //  - 1: return values for coloring the iso-surface (so far not used).
   virtual void MakeGridValues(double *pOut, ct::FMatrixView Grid, uint GridDxOrder, int Role, ct::FMemoryStack &Mem) const = 0;
   virtual double GetDataWeight(double fValue) const; // default: return absolute value.
   virtual QString GetDataUnit() const; // default: 1/a0^3 (atomic unit of density)

   // make a set of atoms which supports (almost all) of the volume's data.
   virtual ct::FAtomSetPtr GetSupportAtoms(ct::FMemoryStack &Mem) const;
   virtual FSupportBox GetSupportBox(FSupportDataPtr pSupportData, double fAbsIsoThreshold, ct::FMemoryStack &Mem) const;
   virtual FSupportDataPtr GetSupportData(ct::FMemoryStack &Mem) const;
   virtual double GetIntegratedValue(FSupportDataPtr pSupportData) const;
   virtual double GetIntegratedValueBelowThreshold(double fValueThreshold, FSupportDataPtr pSupportData) const;
   virtual double GetIntegratedWeight(FSupportDataPtr pSupportData, bool Signed=false) const;
   // compute absolute iso threshold for a given relative iso threshold. Data will be normalized to fRefValue.
   // if fRefValue == 0, fRefValue is set to the numerically integrated value of data in pSupportData.
   virtual double GetAbsoluteIsoThreshold(double fRelativeIsoThreshold, double fRefValue, FSupportDataPtr pSupportData, ct::FMemoryStack &Mem) const;
   virtual size_t GetMaxStackDataSize() const = 0;

   virtual QString GetVolumeType() const = 0;
public:
   FGlMeshPtr
      pGlMesh; // should probably not be here (and wasn't before), but need it for
               // changing color data/iso surface settings.

   void InvalidateRenderCache(); // override
   void BuildRenderCache(FView3d *pView3d); // override
   // just update the colors of the mesh (i.e., rebuild the GL mesh), but leave
   // the mesh itself alone.
   void InvalidateColors();
protected:
   virtual FIndexedTriangleListPtr MakeIsoSurface(FIsoSurfaceSettings const &IsoSurfOpt_) = 0;
};



enum FVolumePropertyType {
   PROPERTY_Density,
   PROPERTY_SpinDensity
};

struct FVolumePropertyInfo {
   FVolumePropertyType
      Type;

   explicit FVolumePropertyInfo(FVolumePropertyType Type_) : Type(Type_) {};
   QString MakeDesc() const;
};


struct FDensityEvalData : public FIntrusivePtrDest
{
   ct::FHeapMatrixPtr
      // occupied orbitals (only orbitals required to evaluate the given density)
      pOrbs;
   TArray<double>
      // total density: OrbWeight = Occupation number
      // spin density: OrbWeight[i] = fOccA[i] - fOccB[i]
      OrbWeights;
   FBasisSetPtr
      // basis set over which pOrbs is expanded.
      pOrbBasis;
   // maybe make RDM, too?
};
typedef ct::TIntrusivePtr<FDensityEvalData>
   FDensityEvalDataPtr;

// This class represents general space-dependent properties of the molecule:
//  - electron density
//  - spin density
//  - later: reactivity indicators on density/vdW surfaces.
// The objects are represented implicitly, via means to compute them on certain points.
// By default they should be represented as iso-surfaces.
struct FVolumeProperty : public FVolumeDataSet
{
   typedef FVolumeDataSet
      FBase;

   explicit FVolumeProperty(FFrame *pWf_, FVolumePropertyInfo const &info_, FVolumeVisualConfigPtr pRefVisConfig_);
   ~FVolumeProperty();

   FVolumeVisualConfigPtr
      pVisConfig;
   void LinkVisualConfig(FVolumeVisualConfigPtr p);

   FVolumePropertyInfo
      info;
   FFrame
      // frame from which's wave function this property is supposed to be computed.
      // TODO: split frame and wave function?
      *m_pWf;

   QString GetDesc(uint Flags=0) const; // override
   QString GetType() const; // override

//    // just update the colors of the mesh (i.e., rebuild the GL mesh), but leave
//    // the mesh itself alone.
//    void InvalidateColors();

   virtual bool DependsOnWaveFunction() const; // override

   void MakeGridValues(double *pOut, ct::FMatrixView Grid, uint GridDxOrder, int Role, ct::FMemoryStack &Mem) const; // override
   QString GetVolumeType() const; // override
   size_t GetMaxStackDataSize() const; // override.
   double GetExactValue() const;
protected:
   FIndexedTriangleListPtr MakeIsoSurface(FIsoSurfaceSettings const &IsoSurfOpt_); // override

   FDensityEvalDataPtr
      pDensityEvalData;
};
typedef ct::TIntrusivePtr<FVolumeProperty>
   FVolumePropertyPtr;


extern double g_fSupportBox_EnlargeRms;



#endif // IV_VOLUME_DATASET_H
