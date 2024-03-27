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

#ifndef CT_DFTGRID_RADIAL_H
#define CT_DFTGRID_RADIAL_H

#include <string> // for passing grid options
#include <map> // for option overrides of atoms and elements

#include "CtDftGrid_CommonDefs.h"
#include "CtDftGrid_Params.h" // FAtomRadialGridDecl & co
#include "CxPodArray.h"
#include "CxRawAtom.h" // for FAtomTag

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

namespace mig { // molecular integration grids

typedef ct::TArray<double>
   FScalarArray;

// full atom radial grid specification/definition/type/generator
class FRadialGridBuilder;



/// A set of radial points (i.e., radii r[i]) and associated integration weights (w[i])
/// representing the *radial* integration rule of a given atom or element.
///
/// Technical Details:
/// - The radii go span the range of 0 to infinity, with 0 on the inside (lower
///   indices), and radii monotonously increasing. Typically the integration
///   limits (0 and infinity) are not explicitly included in the rules.
///   (but for LDAs, including the 0 point might not be a bad idea, actually,
///   and specific rules may include it upon request)
///
/// - The points $r_i$ and weights $w_i$ specify a numerical approximation of
///   the 3D(!) space integral given as follows:
///
///      I[f] := \sum_{i=0}^{N-1} f(r[i]) w[i]
///           \approx \int_0^\infty f(r) 4\pi r^2\,\mathrm{d} r
///
///   That is, they include the 4*M_PI*r^2 volume element in addition to the
///   numerical representation of 'dr'.
///
/// Comments:
/// - The objects are made by one of the FRadialGridBuilder classes
///   (see below)
/// - These objects are used inside the DFT grid generator to cache the radial
///   integration rules of elements, in case an element is encountered multiple
///   times with the same parameter set.
/// - The objects are cached because there are complex advanced integration
///   rules which could, theoretically, be expensive to re-generate the points
///   and weights of (unlike for most standard radial grids, which are just defined
///   by an explicit formula)
class FRadialGrid : public ct::FIntrusivePtrDest
{
   // WARNING: this object takes ownership of the array data in Radii_ and Weights_! The inputs
   // are destroyed! (this is why this constructor is private)
   explicit FRadialGrid(FScalarArray &Radii_, FScalarArray &Weights_);
   friend class FRadialGridBuilder;
public:
   ~FRadialGrid();

   FScalar Radius(size_t iPt) const { return m_Radii[iPt]; }
   FScalar Weight(size_t iPt) const { return m_Weights[iPt]; }
   size_t size() const { assert(m_Radii.size() == m_Weights.size()); return m_Radii.size(); }
protected:
   FScalarArray
      m_Radii,
      m_Weights;
};
typedef ct::TIntrusivePtr<FRadialGrid>
   FRadialGridPtr;
typedef ct::TIntrusivePtr<FRadialGrid const>
   FRadialGridCptr;


struct FAtomSpec {
   size_t
      // number of the atom in the molecule (0-based)
      iAt;
   int
      // element number
      iElement;
   ct::FAtomTag
      // type-of-atom tag as specified in the xyz input data (see CxRawAtom.h).
      // Optional and typically neither provided nor used---but *may* be used
      // to distinguish different atoms of the same element.
      iTag;
   enum FCoreType {
      // assume all-electron atom, with core shells actually being there. Needs
      // larger grids which cover the inner core region.
      CORE_AllElectronNonRelativistic,
      // assume this atom uses a ECP association as the def2/dhf basis sets.
      // These are all-electron for elements <= 36 and are associated with various
      // ECPxxMDF (and a few ECPxxMWB, iirc) for the higher elements.
      CORE_EcpForDef2
   };
   FCoreType
      CoreType;
   FAtomSpec(size_t iAt_, int iElement_, ct::FAtomTag iTag_, FCoreType CoreType_ = CORE_EcpForDef2) : iAt(iAt_), iElement(iElement_), iTag(iTag_), CoreType(CoreType_) {}
};







struct FMetaGridParams_TA { double GridCenterR; double Alpha; };
// struct FMetaGridParams_LogE { double GridCenterR; double Exp; };
struct FMetaGridParams_LogE { double GridCenterR; double Exp; double Compress; };


struct FRadialGridScheme : public ct::FIntrusivePtrDest
{
   virtual ~FRadialGridScheme();
   /// Generate a fixed radial grid at a specific location in memory.
   virtual void MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl) = 0;
   virtual size_t EstimateSizeForTargetAccuracy(double fTargetAccu, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl);
   virtual FSchemeId Id() const = 0;
};
typedef ct::TIntrusivePtr<FRadialGridScheme>
   FRadialGridSchemePtr;



struct FRadialGridScheme_LogE : public FRadialGridScheme
{
   void MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl); // override
   size_t EstimateSizeForTargetAccuracy(double fTargetAccu, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl); // override
   FSchemeId Id() const;
protected:
   FMetaGridParams_LogE MakeDefaultAtomParams(FAtomSpec const &At, FAtomRadialGridDecl const *pDecl);
};


// As defined in Treutler, Ahlrichs - Efficient molecular numerical integration schemes (JCP 1995),
// with additions from 10.1002/jcc.23323, supp info Tab. 1 (replacement for [1], Tab 1.)
struct FRadialGridScheme_TA : public FRadialGridScheme
{
   void MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl); // override
   size_t EstimateSizeForTargetAccuracy(double fTargetAccu, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl); // override
   FSchemeId Id() const;
protected:
   FMetaGridParams_TA MakeDefaultAtomParams(FAtomSpec const &At, FAtomRadialGridDecl const *pDecl);
};

// As defined in Laqua, Kussmann, Ochsenfeld - An improved molecular partitioning scheme for numerical quadratures in density functional theory (JCP 2018)
// It makes TA M4 (alpha=0.6) grids as far as the shape and parameters are
// concerned; however, the default numbers of points is greatly revised and---
// unlike the actual TA parameters---the current setup actually works!
struct FRadialGridScheme_LKO : public FRadialGridScheme_TA
{
   size_t EstimateSizeForTargetAccuracy(double fTargetAccu, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl); // override
   FSchemeId Id() const;
};



/// Main interface for methods of building radial integration grids for concrete
/// elements and atoms.
struct FRadialGridBuilder : public ct::FIntrusivePtrDest
{
   // will register all known schemes by default.
   explicit FRadialGridBuilder();
   virtual ~FRadialGridBuilder();
   /// Generate and return a grid with *exactly* nGridPt radial grid points,
   /// specifically for atom #iAt (where iAt=0,1,...,nAt-1 indicates one of the nAt
   /// atoms of the molecule), which is of element iElement (1=H, 2=He, ...).
   ///
   /// Note:
   /// - This function is meant to be used when we have input options specifying
   ///   an explicit grid size (nradial option). The accuracy-based generator
   ///   function MakeGridForTargetAccuracy may also invoke it after determining
   ///   a suitable number of grid points to achieve a given accuracy with the
   ///   target integration rule.
   /// - iElement may be used to determine element-specific parameters of the
   ///   integration rule (e.g. related to atomic size, atomic density
   ///   distribution, or asymptotic behavior)
   /// - iAt may be used in conjunction with explicit parameter overrides for
   ///   specific atoms.
   ///   Generators need not support this option---iAt may be ignored.
   ///
   /// Default implementation calls the placement-variant of MakeFixedGrid.
   virtual FRadialGridCptr MakeFixedGrid(size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl);

   /// Generate a fixed radial grid at a specific location in memory.
   virtual void MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl);


   /// Estimate the smallest number of grid points for element iElement
   /// compatible with a *global* grid target accuracy (for small molecules with
   /// <= 10 atoms) of fTargetAccu.
   ///
   /// Notes to developers:
   /// - This should be calibrated carefully(!) for the integration rule in question,
   ///   at least for a large set of diatomics at distances between 75% and 800%
   ///   of the equilibrium distance (see diatomics trial in migrid's run_test.py).
   ///
   ///   Not getting this right will cause serious problems in geometry optimizations.
   ///   Many published default parameter sets do *VERY* poorly on this test!
   ///
   /// - This function has a default implementation, which returns default grid
   ///   sizes based on cgk's LogE grids.
   size_t EstimateSizeForTargetAccuracy(double fTargetAccu, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl);

   // ^-- note: neither function is const because we do not really need constness,
   // and the builder objects might want to cache some intermediate data in case
   // they are of a more complicated type.

   void Register(FRadialGridSchemePtr pScheme);
protected:
   typedef std::map<FSchemeId, FRadialGridSchemePtr>
      FIdToSchemeMap;
   FIdToSchemeMap
      m_IdToSchemeMap;
   // attempt to locate a scheme object implementing the scheme described by the scheme id in *pDecl
   FRadialGridScheme *FindScheme(FAtomRadialGridDecl const *pDecl, bool AssertExists = true);
};
typedef ct::TIntrusivePtr<FRadialGridBuilder>
   FRadialGridBuilderPtr;
// typedef ct::TIntrusivePtr<FRadialGridBuilder const>
//    FRadialGridBuilderCptr;



void _MakeRadialGrid_TA(double *ri, double *wi, size_t nPt, FMetaGridParams_TA const &p);
void _MakeRadialGrid_LogE(double *pGridR, double *pGridW, size_t nGridPt, FMetaGridParams_LogE const &p);

double GetSlaterBraggRadius(int iElement);
double GetGridCenter_rExpAvg1(int iElement);
double GetGridCenter_rExpAvg2(int iElement);
double GetGridCenter_TaEta(int iElement);

} // namespace mig





#endif // CT_DFTGRID_RADIAL_H
