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

#ifndef _CT_DFTGRID_PARAMS_H
#define _CT_DFTGRID_PARAMS_H

#include <string>
#include <array> // for meta grid scheme ids
#include <tuple>

#include "CtDftGrid_CommonDefs.h"
#include "CxAtomParamSpec.h"


namespace ct {
   class FBasisSet; // for renormalization transform.
   struct string_slice;
}


namespace mig { // molecular integration grids


// declare a simple fixed-length "string" type which can use to identify
// different types of radial integration schemes in a common, scheme-independent
// parameter structure (FMetaGridParams).
//
// note: this class inherhits the lexicographical comparison operators
// of std::array. So its objects can be used as map keys.
static size_t const SCHEMEID_MaxLength = 8;  // enum { SCHEMEID_MaxLength = 8 };
struct FSchemeId : public std::array<char, SCHEMEID_MaxLength> {
   enum {
      // if set on init, the constructors will not complain if the input
      // string is too long---instead they will simply truncate it.
      AllowTruncation = 0x0001
   };
   FSchemeId() { this->_Init1(0,0,0); }
   explicit FSchemeId(char const *pName, unsigned Flags = 0) { this->_Init1(pName, std::string::npos, Flags); }
   explicit FSchemeId(std::string const &Name, unsigned Flags = 0) { this->_Init1(&Name[0], Name.size(), Flags); }
   operator std::string () const;
   bool IsAssigned() const;
   operator bool () const { return this->IsAssigned(); }
protected:
   void _Init1(char const *pName, size_t NameSize, unsigned Flags);
};


struct FGridDeclBase {
   typedef double
      FScalarParam;
   static FScalarParam const
      // scalar parameters will be set initialized to this. Numerical parameters
      // given by order may also be explicitly specified to ``auto`
      NotSpecified; // = std::numeric_limits<FScalarParam>::max();
   typedef unsigned
      FSizeParam;
   enum {
      DefaultSize = unsigned(-1),
      AdaptiveSize = unsigned(-2)
   };
};

// Represents a generic declaration of atomic grid defaults/overrides, which
// is used as value type for specifications which override meta grid properties
// for general defaults or atom- or element-specific defaults. The info in
// here has no intrinsic meaning, but may be converted to actual meta grids
// by the radial integration scheme, for example.
//
// Note that this object *cannot* contain any reference to the actual element.
// The reason is that there may be a default of those, which is *shared* across
// all elements, to define default overrides.
struct FAtomRadialGridDecl : public FGridDeclBase {
   FSchemeId
      // defines the radial meta-grid handler by id which should be invoked for
      // this specific grid item.
      SchemeId;
   FSizeParam
      // the following options may be used to explicitly fix the number of
      // radial grid points, the number of angular grid points, or the minimal
      // exact angular integration order (L, meaning we should use an angular
      // grid which integrates all spherical harmonics of order <= L exactly).
      // These options may also be set to DefaultSize or AdaptiveSize.
      // Defaults to DefaultSize.
      nRadialPt;

   FScalarParam
      // may be used as override for specific atoms. Defaults to NotSpecified.
      fTargetAccuracy;

   // add a number of (optional) scalar parameters for this radial grid. What,
   // if anything, these parameters actually mean is defined by the radial
   // scheme.
   // This could be used, for example, to overwrite grid center radii or grid
   // shape properties for specific atoms or elements.
   enum { RADIAL_MaxParams = 4 };
   std::array<FScalarParam, RADIAL_MaxParams>
      fParams;

   FAtomRadialGridDecl() { this->_InitDefaults(); }
   explicit FAtomRadialGridDecl(std::string::const_iterator itFirst, std::string::const_iterator itLast);
   explicit FAtomRadialGridDecl(ct::string_slice slDecl);

   // parse atomic grid specification in property list [itFirst...itLast) and update
   // values in *this with its specifications.
   void Update(std::string::const_iterator itFirst, std::string::const_iterator itLast);

   // replace settings in *this with settings from Other if in *this they are undefined
   // or set to default values.
   void Update(FAtomRadialGridDecl const &Other);

   std::string Format() const;

   inline auto SortKey() const noexcept {
      return std::tie(nRadialPt, SchemeId, fTargetAccuracy, fParams);
   }
protected:
   void _InitDefaults();
};

bool operator < (FAtomRadialGridDecl const &A, FAtomRadialGridDecl const &B);

inline void UpdateIfParamAssiged(double &fTarget, FGridDeclBase::FScalarParam const &fParam) {
   if (fParam != FGridDeclBase::NotSpecified)
      fTarget = fParam;
}

template<class FUnsignedInt>
inline void UpdateIfParamAssiged(FUnsignedInt &nTarget, FGridDeclBase::FSizeParam const &nSizeParam) {
   if (nSizeParam != FGridDeclBase::DefaultSize)
      nTarget = nSizeParam;
}



struct FDftGridParams
{
//    explicit FDftGridParams(uint iLevel_ = 3, int iMinL_=0, bool AdjustByAtomicRadii_ = true)
//       : nLevel(iLevel_), AdjustByAtomicRadii(AdjustByAtomicRadii_), iMinL(iMinL_)
   explicit FDftGridParams(int iGridLevel_ = 3);
   explicit FDftGridParams(std::string const &sGridDesc_);

   double fTargetAccuracy() const;
   // note that this does NOT override explicitly set radial or angular parameters!
   void SetTargetAccuracy(double fAcc);
   bool AdjustVoronoiCellsViaAtomicRadii() const { return m_AdjustVoronoiCellsByAtomicRadii; }
   double fWeightCut() const;

   int iMinL(int iPeriod) const;
   int iMaxL(int iPeriod) const;
   int nRadialPt(int iPeriod) const;
   bool HasFixedL(int iPeriod) const;
   // approximately translate target accuracy to Turbomole-style grid "level" (fractional variant)
   double dGridLevel() const;
   // approximately translate target accuracy to Turbomole-style grid "level" (nearest-integer variant)
   int iGridLevel() const;

   enum FVoronoiNeighborScreenType {
      VORONOI_KeepAll,
      VORONOI_ScreenByDistance,
      VORONOI_ScreenByWeightOfClosestPoint
   };
   FVoronoiNeighborScreenType VoronoiNeighborScreenType() const { return m_VoronoiNeighborScreenType; };
   double fVoronoiNeighborScreenThresholdFar() const { return m_VoronoiNeighborScreenThresholdFar; };
   double fVoronoiNeighborScreenThresholdNear() const { return m_VoronoiNeighborScreenThresholdNear; };

   enum FAngularGridAlign {
      // do not orient the angular integration grids
      ANGULARALIGN_None,
      // attempt to assemble a local descriptor of each atom's environment and
      // diagonalize it to assign the atom's grid alignment (likely the best
      // option if it works---that's the default)
      ANGULARALIGN_AtomicEnvironment,
      // generate a random 3D inv-rotation matrix and pull the angular grid
      // through it. This one is meant for grid testing purposes, obviously.
      // (it allows getting a rather good idea about the accuracy which the
      // angular integration yields).
      //
      // Note: The enum is always here, but this will not work if the
      // NO_RANDOM_GRIDS preprocessor flag is set.
      ANGULARALIGN_Randomize
   };
   FAngularGridAlign AngularGridAlign() const { return m_AngularGridAlign; }

//    // decides sets of default radial grid transforms and parameter sets.
//    // Note that parts of those may be overwritten or changed by other options.
//    enum FRadialGridScheme {
//       // radial integration as specified by Treutler & Ahlrichs (M4 T2 alpha = 0.6)
//       RADSCHEME_Ta1995,
//       // Laqua, Kussmann, Ochsenfeld revision of TA scheme
//       RADSCHEME_Lko2018,
//       // that's a cgk special. Essentially Mura-Knowles -R*log(1-x^e) grids...
//       // however, mine are calibrated and have sensible parameter sets.
//       RADSCHEME_LogE,
//       // also a cgk special. A very flexible form which encompasses Mura-Knowles
//       // grids as special cases exactly and TA grids approximately.
//       // Not fully parameterized at the moment.
//       RADSCHEME_LogBeta
//    };
//    FRadialGridScheme RadialGridScheme() const { return m_RadialGridScheme; }



   typedef ct::TAtomParamSpec<ct::TSettingParseProxy_CtorCallDirect<FAtomRadialGridDecl> >
      FRadialGridDecls;
   FAtomRadialGridDecl const &RadialGridDecl(FRadialGridDecls::FAtomIndex iAt, FRadialGridDecls::FElement iElement, FRadialGridDecls::FAtomTag iTag) const {
      return m_RadialGridDecls.get(iAt, iElement, iTag);
   }
   // returns the default value
   FAtomRadialGridDecl const &RadialGridDecl() const { return m_RadialGridDecls.GetDefault(); }

   double GetCoreRadius(int iElement) const;
   int iPrintLevel() const { return m_PrintLevel; }

   ct::FBasisSet const mutable *m_pFitBasis; // FIXME: remove this.

   enum FGridRenormType {
      GRIDRENORM_None,
      GRIDRENORM_Orbitals,
      GRIDRENORM_AuxDensity
   };
   FGridRenormType GridRenormType() const { return m_GridRenormType; };
   double fWeightCut_AtomVdwRadiusFactor() const;
   void SetWeightCut_AtomVdwRadiusFactor(double f) { m_fWeightCut_AtomVdwRadiusFactor = f; }

   bool EvalFreeAtomGridAccuracy() const { return m_EvalFreeAtomGridAccuracy; }
   void SetEvalFreeAtomGridAccuracy(bool bNew){ m_EvalFreeAtomGridAccuracy = bNew; }
protected:
//    double nLevel;
   double
      m_fTargetAccuracy;
//    FRadialGridScheme
//       m_RadialGridScheme;
   bool
      m_EvalFreeAtomGridAccuracy;
   bool
      m_AdjustVoronoiCellsByAtomicRadii;
   double
      m_fWeightCut, // -1: no weight cut. otherwise: points below this weight will be cut off.
      m_fWeightCut_AtomVdwRadiusFactor; // if > 0, smoothly phase out atomic voronoi weights close to this number f * GetVdwRadius_IsoDensity(iElement)
   double
      // if 0: no special core grids;
      // otherwise: radii smaller than this times the element's covalent radius will
      // be considered as belonging to core shells and assigned smaller integration grids.
      m_CoreRadiusThreshold;
   FVoronoiNeighborScreenType
      m_VoronoiNeighborScreenType;
   double
      // only for weight screening
      m_VoronoiNeighborScreenThresholdFar,
      m_VoronoiNeighborScreenThresholdNear; // actual threshold: 1.0 - (this). Exclude neighbor if weight higher than this for this atom and *all* previous atoms.
   int
      m_PrintLevel;
   FAngularGridAlign
      m_AngularGridAlign;
   FGridRenormType
      m_GridRenormType;

   FRadialGridDecls
      m_RadialGridDecls;

   struct FPeriodGridDef {
      // minimum L which this grid should expand exactly. -1 means 'auto'.
      int iMinL, iMaxL, nRadialPt;
      FPeriodGridDef() : iMinL(-1), iMaxL(-1), nRadialPt(-1) {}
   };
   FPeriodGridDef
      m_PeriodGridDefs[DFTGRID_nPeriodParams];
   void SetParamsViaDesc(std::string const &sGridDesc_);
   void SetParamsViaLevel(int iGridLevel_);
   void InitDefaults();
   double GetApproxLevel() const;
};

} // namespace mig

#endif // _CT_DFTGRID_PARAMS_H
