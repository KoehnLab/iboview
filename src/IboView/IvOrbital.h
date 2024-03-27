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

#ifndef IV_ORBITAL_DATASET_H
#define IV_ORBITAL_DATASET_H

#include "IvVolumeDataSet.h"

// #include <map>
// #include <set>
// #include <vector>
// #include <list>
//
// #include <QString>
// #include <QObject>
// // #include <QAbstractTableModel>
// // #include <QAction>
// // #include <QExplicitlySharedDataPointer>
// // #include <QSharedData>
//
//
// #include "CtAtomSet.h"
// #include "CtBasisSet.h"
// #include "CxPodArray.h"
//
// // #include "CxColor.h"
// #include "IvMesh.h"
// #include "IvGl.h"
// #include "IvIsoSurface.h"
//
// #include "IvDataOptions.h"
// #include "IvDataSet.h"
//
// // #include "IvAnalysis.h"
// #include "CtMatrix.h"
// // #include "IvLog.h"
// // #include "IvTables.h"
// // #include "IvIrc.h"
// #include "CxVec3.h"
// // #include "CtWfi.h"
//
//
// // #include "CtDftGrid_ivb.h" // FIXME: merge back with main version once done.
// #include "CtDftGrid.h"


class FFrame;
class FView3d;
class FDocument;

using wfi::FOrbitalSpin;
using wfi::FWfType;
namespace ct {
   struct FHfOptions;
   struct FWfDecl;
}


struct FOrbitalInfo {
   double
      fEnergy,
      fOcc;
   int
      iSym;
   FOrbitalSpin
      Spin;
   int
      iOriginalIndex;
   explicit FOrbitalInfo(double fEnergy_=0, double fOcc_=0., int iSym_=0, FOrbitalSpin OrbSpin_=wfi::ORBSPIN_Unknown);
   QString MakeDesc(QString RefName) const;
   QString MakeDesc(size_t iOrb) const;

   double fChargeFactor() const;
   double fSpinFactor() const;
};


struct FOrbital : public FVolumeDataSet
{
   typedef FVolumeDataSet
      FBase;

   explicit FOrbital(QString const &Desc_, ct::FAtomSetPtr pAtoms_, FDocument *pDocument_, FBasisSetPtr pBasisSet_, double *pCoeffs_, FOrbitalInfo const &info_, FVolumeVisualConfigPtr pRefVisConfig, double *pDm=0, double *pQm=0);
   ~FOrbital();


   bool HaveMoments;
   double vDipMom[3]; // dipole moment and quadrupole moment.
   double mQuadMom[3][3];


   FOrbitalInfo
      info;
   FBasisSetPtr
      pBasisSet;
   TArray<double>
      pCoeffs;
      // ^- if this is a orbital-like-quantity, then these are its expansion
      //    coefficients. over pBasisSet.

   // minimal basis and IAO coefficients of the orbitals; used for tracking
   // orbital changes.
   FBasisSetPtr
      pMinBasis;
   TArray<double>
      pIaoCoeffs;
   TArray<double> MakeIaoCharges(bool UseOccupationNumbers, double fOccupancyOtherwise=2.) const; // make IAO charges from IAO coeffs.
   void AddChargeContributions(FChargeAnalysis &Charges) const;

   QString GetType() const; // override

   enum FFullDescOptions {
      ORBDESC_ChargesOnly = 0x01,
      ORBDESC_CompactSpace = 0x02
   };
   QString MakeFullDesc(double ThrPrint = 0.02, uint Flags = 0, int nMaxAtoms = -1) const;

   QString GetDesc(uint Flags=0) const; // override
   void UpdateDescFromInfo(size_t iOrb);

   // multiply coefficients by -1.
   void FlipPhase();

   virtual uint32_t GetBaseColor() const;
   virtual bool DependsOnWaveFunction() const; // override

   void MakeGridValues(double *pOut, ct::FMatrixView Grid, uint GridDxOrder, int Role, ct::FMemoryStack &Mem) const; // override
   QString GetVolumeType() const; // override

   ct::FAtomSetPtr GetSupportAtoms(ct::FMemoryStack &Mem) const; // override
   double GetDataWeight(double fValue) const; // override
   QString GetDataUnit() const; // override
   size_t GetMaxStackDataSize() const; // override.
protected:
   FIndexedTriangleListPtr MakeIsoSurface(FIsoSurfaceSettings const &IsoSurfOpt_); // override
};
typedef ct::TIntrusivePtr<FOrbital>
   FOrbitalPtr;


bool MakeIboChangeCurve(TArray<float> &CurveData, uint iRow, FDocument *document);

#endif // IV_ORBITAL_DATASET_H
