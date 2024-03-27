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

#include "CxTypes.h"
#include "CxVec3.h"
#include "CxMemoryStack.h"
// #include "CtAtomSet.h"
// #include "CtMatrix.h"
#include "CxPodArray.h"
#include "CxRawAtom.h"
#include "CxDensityModel.h"

namespace mig { // molecular integration grids
typedef ct::TVector3<double>
   FVec3d;


// // void MakeChiMatrixFromExplicitRadii(FMatrixView ChiAB_, double const *pAtomRadii, size_t nAtoms);
//
// // a simple atom structure used *only* to identify where atoms are, and what element they have.
// // This one is used in place of FAtomSet objects for interfacing purposes, so that
// // mostly-autarc parts of the code, such as the space partitioning stuff, can be isolated
// // from the ct base machinery
// //
// // There is a RawAtomListFromAtomSet function to convert FAtomSet objects to these guys.
// struct FRawAtom {
//    FVec3d
//       vPos;
//    ptrdiff_t
//       iElement;
//    FRawAtom() {};
//    FRawAtom(FVec3d vPos_, ptrdiff_t iElement_) : vPos(vPos_), iElement(iElement_) {}
// };
// typedef TArray<FRawAtom>
//    FRawAtomList;
//
//
// struct FXyzFrame;
// FRawAtomList RawAtomListFromAtomSet(FXyzFrame const &Atoms);


struct FVoronoiPartitionParams : public ct::FIntrusivePtrDest
{
   enum FSmoothingFnType {
      SMOOTHFN_Becke_k4,
      SMOOTHFN_PolyStep_a16,
      SMOOTHFN_PolyStep_a12,
      SMOOTHFN_PolyStep_a10,
      SMOOTHFN_Murray_k12 = SMOOTHFN_PolyStep_a12, // <- was my default for a long while
      SMOOTHFN_Murray_k10 = SMOOTHFN_PolyStep_a10,
      SMOOTHFN_PolyStep_a4sq,
      SMOOTHFN_TrigStep_k4,
      SMOOTHFN_Stratman
   };
   FSmoothingFnType
      // describes the smoothing function to use for the smooth voronoi polyhedra
      m_SmoothingFnType;

   enum FAtomSizeAdjustType {
      ATOMSIZE_None,
      ATOMSIZE_Tfvc,
      ATOMSIZE_OriginalBecke,
      ATOMSIZE_TaModOfBecke // Treutler-Ahlrichs modification of the Becke assignment. It replaces the radii by their square roots, but is otherwise identical.
   };
   FAtomSizeAdjustType
      // describes the function to use for transforming a pair of (mu,chi) into nu for making atomic size adjustments.
      m_AtomSizeAdjustType;

   double
      m_ThrPairVoronoiWeightCut;
   double
      m_fWeightCut_AtomVdwRadiusFactor; // ignored if <= 0.

   bool
      m_CutPairVoronoiAtMaxAtRange;
   ct::TArray<double>
      // length nAt
      m_MaxAtRanges;

   explicit FVoronoiPartitionParams(FSmoothingFnType SmoothingFnType_ = SMOOTHFN_PolyStep_a4sq, FAtomSizeAdjustType AtomSizeAdjustType_ = ATOMSIZE_Tfvc);

//    double const
//       // atomic size parameters Chi_{AB} are expected at pChiAB[A*iRowSt_ChiAB, B*iColSt_ChiAB]
//       *pChiAB;
//    size_t
//       iRowSt_ChiAB, iColSt_ChiAB;
//    double GetOutsideChiAB(size_t iAtA, size_t iAtB) const { return pChiAB[iAtA * iRowSt_ChiAB + iAtB * iColSt_ChiAB]; }
};
typedef ct::TIntrusivePtr<FVoronoiPartitionParams>
   FVoronoiPartitionParamsPtr;



// this is an auxiliary structure which can be used to generate the ChiAB input data
// for cases in which the atomic sizes are derived from individual atomic radii.
struct FVoronoiChiData
{
   typedef bool (*FTreatAsNeighborsFn)(double fDistAB, ct::FRawAtom const &A, double fRadA, ct::FRawAtom const &B, double fRadB);

   // note: use RSolv = -1 to use DftGridMode instead of CosmoMode
   void MakeFromAtomRadii(ct::FRawAtomList const &Atoms_, double const *pAtomRadii_, double RSolv_, FTreatAsNeighborsFn TreatAsNeighborsFn, FVoronoiPartitionParams::FAtomSizeAdjustType AtomSizeAdjustType_);
   void MakeFromDensityModel(ct::FRawAtomList const &Atoms_, ct::FDensityModel const *pDensityModel, FTreatAsNeighborsFn TreatAsNeighborsFn_, FVoronoiPartitionParams::FAtomSizeAdjustType AtomSizeAdjustType_, ct::FMemoryStack &Mem);

   ct::TArray<double>
      m_ChiAB;
   size_t
      m_iRowSt, m_iColSt;

   double &ChiAB(size_t iAt, size_t jAt) { return m_ChiAB[iAt * m_iRowSt + jAt * m_iColSt]; }
};


// note: moved from IvAnalysis.cpp. Depends on CxPoly.cpp which in turn depends
// on CxAlgebra. Could be a portability issue. But in the worst case, we could just make
// add an #ifdef around it...
// Returns line fraction Mu (-1... +1) at which the density becomes minimal.
double FindMinimumDensityOnLine(double const *pfLineFraction, double const *pfDensity, size_t nPts);

// compute Chi_AB fuzzy voronoi cell shift parameter from given A--B line density profile and line fractions.
// Line fractions go from (-1. to +1.0).
double FindTfvcChiFromDensityProfile(double const *pfDensity, double const *pfLineFraction, size_t nPts);



// Currently used for TFVC (in IboView), COSMO cavities (in MicroScf), and
// for DFT grid generator (migrid)
// TODO: change 'TreatAsNeighborsFn' to allow continuous fade out.
struct FVoronoiPartition : public ct::FIntrusivePtrDest
{
//    typedef FVoronoiPartitionParams::FSmoothingFnType
//       FSmoothingFnType;
//    typedef FVoronoiPartitionParams::FAtomSizeAdjustType
//       FAtomSizeAdjustType;
   typedef double (*FSmoothingFn)(double Nu);
   typedef double (*FAtomSizeAdjustFn)(double Mu, double Chi);

   explicit FVoronoiPartition(ct::FRawAtomList const &Atoms_, FVoronoiPartitionParams const &Params_, double const *pChiAB_, size_t iRowSt_, size_t iColSt_);

   // this function subdivides the entire space volume into atomic contributions.
   // Returns weight of atom iAtom at point vPos.
   double GetAtomWeight(FVec3d const &vGridPos, size_t iAtom, ct::FMemoryStack &Mem) const;
   // as GetAtomWeight, but consider only other atoms in pAtOrd[0..nAt] for partitioning
   double GetAtomWeight(FVec3d const &vGridPos, size_t iAtom, size_t *pAtOrd, size_t nAt, ct::FMemoryStack &Mem) const;
protected:
   ct::FRawAtomList
      m_Atoms;
   FVoronoiPartitionParams
      // note: this function makes a copy.
      m_Params;
//    double
//       m_PairVoronoiWeightCut;
//    FSmoothingFnType
//       m_SmoothingFnType;
//    FAtomSizeAdjustType
//       m_AtomSizeAdjustType;
   FSmoothingFn
      m_pSmoothingFn;
   FAtomSizeAdjustFn
      m_pAtomSizeAdjustFn;
   ct::TArray<double>
      // both nAt x nAt arrays.
      m_ChiAB,
      m_pInvDistAt;
//    void SetupVoronoiInfo(FMatrixView ChiAB_);
   void SetupVoronoiInfo(double const *pChiAB_, size_t iRowSt_, size_t iColSt_);


   double GetPairVoronoiR(double Mu, size_t iAtom, size_t iOtherAtom) const;
   void SortAtomsByDistanceToPoint(double *&pDistAg, size_t *&pAtOrd, size_t &nAtOut, FVec3d const &vGridPos, ct::FMemoryStack &Mem) const;
   double GetAtomWeight(FVec3d const &vGridPos, size_t iAtom, double *pDistAg, size_t *pAtOrd, size_t nAt, ct::FMemoryStack &Mem) const;
};
typedef ct::TIntrusivePtr<FVoronoiPartition>
   FVoronoiPartitionPtr;
typedef ct::TIntrusivePtr<FVoronoiPartition const>
   FVoronoiPartitionCptr;

} // namespace mig
