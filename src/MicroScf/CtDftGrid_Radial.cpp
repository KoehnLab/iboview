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
#include <stdexcept>
#include "CtDftGrid_Radial.h"

#ifdef HAVE_BOOST_SPECIAL_FUNCTIONS
   #include <boost/math/special_functions/beta.hpp>
   #include <boost/math/special_functions/gamma.hpp>
   // ^- TODO: replace by homebrew? The complete beta and gamma functions actually
   // are part of C++11 and up. I'd just need the incomplete beta...
#endif // HAVE_BOOST_SPECIAL_FUNCTIONS

namespace mig {




void _MakeRadialGrid_LogE_BaseVersion(double *pGridR, double *pGridW, size_t nGridPt, FMetaGridParams_LogE const &p)
{
   // compute what we need to multiply to -log(1-x^e) in order to get the center
   // of the grid (at x = 0.5) to the externally provided GridCenterR.
   double
      fRadiusScale;
   fRadiusScale = p.GridCenterR / (-std::log1p(-std::pow(.5, p.Exp)));

   double
      den = 1./(nGridPt+1);
   for (size_t i = 0; i < nGridPt; ++ i) {
      double
         x = (1. + double(i))*den,
         dx = den;

      if (0) {
         // replace integration weight rule by Chebychev 2nd (see comments in
         // _MakeRadialGrid_TA on why this indeed is ChebT2).
         double
            t = -std::cos(x * M_PI),
            dt = dx * M_PI * std::sin(x * M_PI);
         // trafo from t in [-1., 1.] to x in [0., 1.]
         t = .5 + .5 * x; dt *= .5;
         x = t; dx = dt;
      }

      double
         powe1 = std::pow(x, p.Exp - 1.),
         r = -fRadiusScale * std::log1p(-x*powe1),
         dr = fRadiusScale * p.Exp * powe1/(1.-x*powe1),
         w = dx * dr * 4.*M_PI*r*r; // last: volume element (4 pi r^2)
      pGridR[i] = r;
      pGridW[i] = w;
   }
//    std::cout << fmt::format("RADIAL-GRID: nPt = {}  R = {}  e = {}  -->  rScale = {}\n", nGridPt, p.GridCenterR, p.Exp, fRadiusScale);
}

void _MakeRadialGrid_LogE(double *pGridR, double *pGridW, size_t nGridPt, FMetaGridParams_LogE const &p)
{
   // compute what we need to multiply to -log(1-x^e) in order to get the center
   // of the grid (at x = 0.5) to the externally provided GridCenterR.
   double
      fRadiusScale;
   fRadiusScale = p.GridCenterR / (-std::log1p(-std::pow(.5, p.Exp)));

   double
      den = 1./(nGridPt+1);
   for (size_t i = 0; i < nGridPt; ++ i) {
      double
         x = (1. + double(i))*den,
         dx = den;

      double t = x, dt = dx, s, ds;

      if (0) {
         // replace integration weight rule by Chebychev 2nd (see comments in
         // _MakeRadialGrid_TA on why this indeed is ChebT2).
         s = .5 - .5*std::cos(t * M_PI),
         ds = .5 * M_PI * std::sin(t * M_PI);
         t = s; dt *= ds;
      }

      {
         // t -> s(t) = t^e
         double e = p.Exp;
         double powe1 = std::pow(t, e - 1.);
         s = t*powe1;
         ds = e * powe1;
         t = s; dt *= ds;
      }

      {
         // t -> s(t) = -log(1-t)
         s = -std::log1p(-t);
         ds = 1/(1-t);
         t = s; dt *= ds;
      }

      if (p.Compress != 0) {
         double e = 1 - p.Compress;
         double powe1 = std::pow(t, e - 1.);
         s = t*powe1;
         ds = e * powe1;
         t = s; dt *= ds;
      }

      s = t*fRadiusScale;
      ds = fRadiusScale;
      t = s; dt *= ds;

      double
         r = t,
         w = dt * 4.*M_PI*r*r; // last: volume element (4 pi r^2)
      pGridR[i] = r;
      pGridW[i] = w;
   }
//    std::cout << fmt::format("RADIAL-GRID: nPt = {}  R = {}  e = {}  -->  rScale = {}\n", nGridPt, p.GridCenterR, p.Exp, fRadiusScale);
}


// this one is a literal implementation of the Treutler & Ahlrichs radial grids for a given fixed radial
// note: ri and wi allocated on Mem.
void _MakeRadialGrid_TA(double *ri, double *wi, size_t nPt, FMetaGridParams_TA const &p)
{
   double const
      inv_ln2 = 1./std::log(2.);
   double const
      den1 = 1/double(nPt + 1);
   for (size_t i = 0; i < nPt; ++ i) {
      FScalar
         xi, dx, r, dr_by_dx;

      // Chebychev 2nd kind
      {
         double ii = (1. + double(i))*den1;
         xi = -std::cos(ii * M_PI);  // iterates over [-1,1].
//          dx = M_PI/double(nPt+1) * sqr(std::sin(ii * M_PI)) / std::sqrt(1. - sqr(xi));
         // ^-- hm.. sqrt(1 - t^2) with t=cos(y) is just sin(y). Apparently
         //     ChebT2 is identical to a trapz rule on [0,1] with variable substitution x
         //     -> t(x) := cos(M_PI*x)... (which is just a TrigStpe scaled to
         //     [-1,1]...). And ChebT1 is the corresponding midpt rule with the same
         //     variable transform. Well... that is pointless.
         dx = M_PI * den1 * std::sin(ii * M_PI);
      }

      // note: xi now in range [-1,+1].
      {
         // Treutler/Ahlrichs M3 (alpha=0)/M4... kind of works (now that a working xi/wi is used...).
         double
            PowAlpha = std::pow(1 + xi, p.Alpha),
            Log2x = -std::log((1.-xi)/2.);
         r = (p.GridCenterR*inv_ln2) * PowAlpha * Log2x;
         // Wolfram Alpha says: D[(1+x)^α * Log[2/(1-x)],{x,1}] is R*inv_ln2*PowAlpha * ((Alpha * Log2x)/(xi+1.) + 1./(1. - xi)).
         // ...and that is exactly the same thing as I had before:
         dr_by_dx = (p.GridCenterR*inv_ln2) * PowAlpha*(1./(1.-xi) + p.Alpha*Log2x/(1.+xi));
      }

      ri[i] = r;
      wi[i] = dx * dr_by_dx* 4.*M_PI * r*r;
   }
}


// Approximately convert our grid target accuracy specifications into Turbomole
// grid levels. This one returns a float (Turbomole grid themselves are
// specified in terms of integers); see iTaGridLevel for a rounded-integer
// variant.
static double dTaGridLevel(double fTargetAccuracy) {
   return 1. + 2.*(-std::log10(fTargetAccuracy) - 4.);
   // ^-'grid 1'. That has target acc 1e-4, 'grid 3' has target acc 1e-5.
}

// return an estimate for a Turbomole-style grid level (as integer) which
// would approximately correspond to our target accuracy settings of grids
static int iTaGridLevel(double fTargetAccuracy) {
   return int(std::round(dTaGridLevel(fTargetAccuracy)));
}


// Returns the radial grid sizes as specified in the TA 1995 article
// (text and Tab III). These are the original variants, exactly as the
// article describes them. See also LKO variant below for more realistic
// estimates.
size_t _ComputeRadialGridSize_OriginalTA(double fTargetAccuracy, int iElement)
{
   // Turbomole's radial grid size offsets. We set them from Tab III, but for
   // the (unspecified) "grid 0" (so everything from grid 1 reduced by 5 pts),
   // so that we can get the first target grids by just multiplying the grid "level"
   // by 5.
   size_t
      nRadialPt_Grid0 = 45;
   if (iElement <= 54) nRadialPt_Grid0 = 40; // < Cs -- Rb-Xe
   if (iElement <= 36) nRadialPt_Grid0 = 30; // < Rb -- K-Kr
   if (iElement <= 18) nRadialPt_Grid0 = 25; // < K -- Na-Ar
   if (iElement <= 10) nRadialPt_Grid0 = 20; // < Na -- Li-Ne
   if (iElement <= 2) nRadialPt_Grid0 = 15; // < Li -- H-He

   // gridsize 1 2 3 4 5  6  7 8 9
   // radsize  1 2 3 6 8 10 14 9 3
   // H-He: 20; Li-Ne: 25; Na-Ar: 30; K-Kr: 40; Rb-Xe: 45; Cs-Lw: 50;
   double dGridLevel = dTaGridLevel(fTargetAccuracy);
   size_t nRadialPt = nRadialPt_Grid0;

   // Tab III says that for grids 1 up to grid 4 there are 5 radial points per
   // level, afterwards 10 (I guess... TA Tab III only specifies up to grid 5)
   double fThr5to10 = 4;
   double fExtraPt = 5*std::min(dGridLevel, fThr5to10) + 10*std::max(dGridLevel - fThr5to10, 0.);
   if (fExtraPt > 0) {
      nRadialPt += size_t(.5 + fExtraPt);
   }
   return nRadialPt;
};




enum FStaticTableFlags {
   // indicate that element at index 0 of the table is a dummy, so that
   // hyrdogen is at index 1.
   TABLE_HydrogenIsSecondEntry = 0x0001,
   // indicate that we should replicate the last value in the table in
   // case we run into an untablulated out of range element
   TABLE_ReplicateLastEntry = 0x0002
};

template<class FDataEntry>
static FDataEntry const &GetTabulatedValue(int iElement, unsigned Flags, FDataEntry const *pTable, size_t nTableSizeInBytes) {
   size_t
      nEntriesInTable = nTableSizeInBytes/sizeof(pTable[0]);
   assert(bool(Flags & TABLE_HydrogenIsSecondEntry) && bool(Flags & TABLE_ReplicateLastEntry));
   if (iElement <= 0)
      return pTable[1]; // replace by H
   else if (size_t(iElement) >= nEntriesInTable)
      // replace by last element we actually have
      return pTable[nEntriesInTable-1];
   else
      return pTable[size_t(iElement)]; // <- note: NOT -1! There is a dummy item at [0] such that table index and element correspond 1:1.
   (void)Flags; // suppress unused warning. Currently only applied for debug checks.
}


// I stole these from psi4/cubature.cc. psi4 says " // Bragg-Slater radii  J.C. Slater JCP 41, (1964), 3199 [bohr]", but that is most definitely not
// what these are (I have those below). Most of those elements are not in the table there.
static FScalar const s_SlaterBraggRadiiStoledFromPsi4[] = { -1,
   0.661, 0.661, // H-He
   2.740, 1.984, 1.606, 1.323, 1.228, 1.134, 0.945, 0.900, //-Ne
   3.402, 2.835, 2.362, 2.079, 1.890, 1.890, 1.890, 1.890, //-Ar
   4.157, 3.402, 3.024, 2.656, 2.551, 2.656, 2.656, 2.656, 2.551, 2.551, 2.551, 2.551, 2.457, 2.362, 2.173, 2.173, 2.173, 2.173, // -Kr
   4.441, 3.780, 3.402, 2.929, 2.740, 2.740, 2.551, 2.457, 2.551, 2.646, 3.024, 2.929, 2.929, 2.740, 2.740, 2.646, 2.646, 2.646, // -Xe
   4.913, 4.063,
            3.685, 3.496, 3.496, 3.496, 3.496, 3.496, 3.496, 3.402, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307,//*La-Lu
   2.929, 2.740, 2.551, 2.551, 2.457, 2.551, 2.551, 2.551, 2.835, 3.591, 3.024, 3.024, 3.591, 3.591, 3.591, // Rn
   4.063, 4.063,                                                                //Fr-Ra
            3.685, 3.401, 3.401, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307, 3.307,//Ac-Lr
};


// made by ~/dev/migrid/radial/grid_test.py. For H-Rn (86)
// WARNING: these guys are not just multiplied --- they need to go where x=0.5 maps to!
static double const s_GridCenters_rExpAvg1[] = {-1.00,1.5,0.927,1.67,1.53,1.36,1.19,
   1.05,0.949,0.862,0.787,0.986,1.02,1.05,1.03,0.997, 0.967,0.93,0.892,1.02,1.06,
   1.04,1,0.964,0.832,0.893,0.862,0.832,0.82,0.704,1.19,0.754,0.756,0.751,0.747,
   0.739,0.729,2.11,2.14, 2.05,1.96,1.79,1.63,1.63,1.48,1.42,1.37,1.31,1.33,1.35,
   1.34,1.32,1.3,1.29,1.25,2.49,2.53,0.883,0.893,0.902,0.883,0.863,0.843,0.822,
   0.795,0.805,0.778,0.754,0.738,0.723,0.706,0.674,1.87,1.81,1.72,1.68,1.63,
   1.57,1.45,1.44,1.44,1.48,1.47,1.46,1.45,1.43,1.4};
// WARNING: these guys are not just multiplied --- they need to go where x=0.5 maps to!
// made by ~/dev/migrid/radial/grid_test.py for rExpr=2
static double const s_GridCenters_rExpAvg2[] = {-1.00,1.73,1.09,2.49,2.08,1.78,1.51,1.31,1.18,1.06,0.962,1.57,1.57,1.6,1.51,1.42,1.35,1.27,1.2,1.64,1.68,1.64,1.55,1.48,1.21,1.35,1.3,1.25,1.22,1,1.61,1.14,1.14,1.11,1.1,1.08,1.05,2.6,2.63,2.43,2.3,2.06,1.86,1.89,1.67,1.61,1.56,1.52,1.58,1.65,1.64,1.61,1.58,1.55,1.51,2.96,3,1.01,1.05,1.06,1.04,1.02,0.996,0.97,0.944,0.988,0.94,0.895,0.878,0.862,0.835,0.786,2.16,2.07,1.96,1.9,1.85,1.77,1.61,1.61,1.64,1.74,1.74,1.72,1.7,1.67,1.64};

// "eta" radius scale parameters (=grid centers with the TA formulas), from
// 10.1002/jcc.23323, supp info Tab. 1 (replacement for [1], Tab 1.)
static double s_GridCenters_TaEta[] = { 0.,
 0.8, 0.9, 1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9, 0.9, 1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.5, 1.4, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, // H - Fe
 1.2, 1.1, 1.1, 1.1, 1.1, 1.0, 0.9, 0.9, 0.9, 0.9, 1.4, 1.4, 1.1, 1.3, 1.0, 1.2, 0.9, 0.9, 0.9, 1.0, 0.9, 1.0, 1.0, 1.3, 1.2, 1.2, // Co-Te
 0.9, 1.0, 1.7, 1.5, 1.5, 1.3, 1.3, 1.4, 1.8, 1.4, 1.2, 1.3, 1.3, 1.4, 1.1, 1.1, 1.2, 1.6, 1.4, 1.3, 1.2, 1.0, 1.0, 0.9, 1.3, 1.2, // I - Pt
 1.2, 1.0, 1.2, 1.2, 1.1, 1.2, 1.1, 2.1, 2.2, 1.8, 1.7, 1.3, 1.4, 1.2, 1.2, 1.3, 1.4, 1.4, 1.7, 1.9, 1.9, 2.0, 2.0, 1.6, 2.0 // Au-Lw
};


FScalar GetSlaterBraggRadius(int iElement) {
   return GetTabulatedValue(iElement, TABLE_HydrogenIsSecondEntry | TABLE_ReplicateLastEntry,
         &s_SlaterBraggRadiiStoledFromPsi4[0], sizeof(s_SlaterBraggRadiiStoledFromPsi4));
}

double GetGridCenter_rExpAvg1(int iElement) {
   return GetTabulatedValue(iElement, TABLE_HydrogenIsSecondEntry | TABLE_ReplicateLastEntry,
         &s_GridCenters_rExpAvg1[0], sizeof(s_GridCenters_rExpAvg1));
}

double GetGridCenter_rExpAvg2(int iElement) {
   return GetTabulatedValue(iElement, TABLE_HydrogenIsSecondEntry | TABLE_ReplicateLastEntry,
         &s_GridCenters_rExpAvg2[0], sizeof(s_GridCenters_rExpAvg2));
}


double GetGridCenter_TaEta(int iElement) {
   return GetTabulatedValue(iElement, TABLE_HydrogenIsSecondEntry | TABLE_ReplicateLastEntry,
         &s_GridCenters_TaEta[0], sizeof(s_GridCenters_TaEta));
}





FRadialGrid::FRadialGrid(FScalarArray &Radii_, FScalarArray &Weights_)
{
   if (Radii_.size() != Weights_.size())
      throw std::runtime_error("FRadialGrid::c'tor: arrays of integration points and weights are inconsistent.");
#ifdef _DEBUG
   for (size_t i = 0; i < Radii_.size(); ++ i) {
      if (Radii_[i] < 0. || !check_if_finite_if_supported(Radii_[i]))
         throw std::runtime_error("FRadialGrid::c'tor: encountered invalid negative or non-finite radius.");
   }
#endif
   m_Radii.swap(Radii_);
   m_Weights.swap(Weights_);
}


FRadialGrid::~FRadialGrid()
{
   // just to bind the destructor and vtable to this particular .cpp file.
}



FRadialGridBuilder::FRadialGridBuilder()
{
   if (1) {
      // register recognized default schemes
      Register(FRadialGridSchemePtr(new FRadialGridScheme_TA()));
      Register(FRadialGridSchemePtr(new FRadialGridScheme_LKO()));
      Register(FRadialGridSchemePtr(new FRadialGridScheme_LogE()));
   }
}


FRadialGridBuilder::~FRadialGridBuilder()
{
}


FRadialGridCptr FRadialGridBuilder::MakeFixedGrid(size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl)
{
   FRadialGridScheme
      *pScheme = FindScheme(pDecl);
   FScalarArray r(nGridPt), w(nGridPt);
   pScheme->MakeFixedGrid(&r[0], &w[0], nGridPt, At, pDecl);
   return FRadialGridCptr(new FRadialGrid(r, w));
}


void FRadialGridBuilder::MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl)
{
   FRadialGridScheme
      *pScheme = FindScheme(pDecl);
   return pScheme->MakeFixedGrid(pRadii, pWeights, nGridPt, At, pDecl);
}


size_t FRadialGridBuilder::EstimateSizeForTargetAccuracy(double fTargetAccuracy, FAtomSpec const &Atom, FAtomRadialGridDecl const *pDecl) {
   FRadialGridScheme
      *pScheme = FindScheme(pDecl);
   return pScheme->EstimateSizeForTargetAccuracy(fTargetAccuracy, Atom, pDecl);
}


void FRadialGridBuilder::Register(FRadialGridSchemePtr pScheme) {
   FSchemeId id = pScheme->Id();
   FIdToSchemeMap::iterator
      it = m_IdToSchemeMap.find(id);
   if (it != m_IdToSchemeMap.end())
      throw std::runtime_error(fmt::format("FRadialGridBuilder::Register(): tried to register scheme with id '{}', but this scheme id is already used.", static_cast<std::string>(id)));
   m_IdToSchemeMap[id] = pScheme;
}


FRadialGridScheme::~FRadialGridScheme()
{
}

FRadialGridScheme *FRadialGridBuilder::FindScheme(FAtomRadialGridDecl const *pDecl, bool AssertExists)
{
   assert(pDecl != 0);
   FSchemeId const
      &id = pDecl->SchemeId;
   // todo: do something special if the scheme id is not assigned?
   // I think we can't really mix & match here. Would be the DFT grid gen which
   // should merge defaults from different places. We do not actually have
   // the full parameter sets to do that here...
   FIdToSchemeMap::iterator
      it = m_IdToSchemeMap.find(id);
   if (it == m_IdToSchemeMap.end()) {
      if (AssertExists)
         throw std::runtime_error(fmt::format("FRadialGridBuilder::FindScheme(): scheme with id '{}' not registered.", static_cast<std::string>(id)));
      return 0;
   }
   return it->second.get();
}





size_t FRadialGridScheme::EstimateSizeForTargetAccuracy(double fTargetAccuracy, FAtomSpec const &Atom, FAtomRadialGridDecl const *pDecl) {
   // WARNING: this is the default implementation. It is only here to simplify getting
   // started by providing some very rough estimates of grid sizes one may expect to need.
   //
   // This one concretely returns a (typically untrustworthy) estimate of the
   // grid size expected for a certain accuracy based on what the TA article
   // specified. However, in practice the returned grids are often too small.
   // See the LKO revision for more realistic estimates.
   return _ComputeRadialGridSize_OriginalTA(fTargetAccuracy, Atom.iElement);
   (void)pDecl; // suppress unused warning
}


FSchemeId FRadialGridScheme_LogE::Id() const { return FSchemeId("loge"); }

size_t FRadialGridScheme_LogE::EstimateSizeForTargetAccuracy(double fTargetAccuracy, FAtomSpec const &Atom, FAtomRadialGridDecl const *pDecl)
{
   // TODO: handle ECP specials. The elements below Kr need far fewer
   // points than we make here! See fit program.
   double fRadialAccu = 1e-1 * fTargetAccuracy;
   // ^-- hm... or 1e-2? that good? maybe 1e-1 is enough?
   double ldAccu = -std::log(fRadialAccu)/std::log(10.);
   unsigned ShouldNotBeHere = 0;
   size_t
      nRadialPt = ShouldNotBeHere + int(.5 + (2 + 5*std::pow(ldAccu, (std::pow(Atom.iElement, (1./16.))))));
   return nRadialPt;
   (void)pDecl; // suppress unused warning
}

void FRadialGridScheme_LogE::MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl) {
   _MakeRadialGrid_LogE(&pRadii[0], &pWeights[0], nGridPt, MakeDefaultAtomParams(At, pDecl));
}

FMetaGridParams_LogE FRadialGridScheme_LogE::MakeDefaultAtomParams(FAtomSpec const &At, FAtomRadialGridDecl const *pDecl) {
   double
      e = 1. + std::pow(double(At.iElement), 1./4.),
      R = GetGridCenter_rExpAvg2(At.iElement) / (e-1),
      Compress = 0;
   if (pDecl != 0) {
      assert(pDecl->fParams.size() >= 2);
      UpdateIfParamAssiged(R, pDecl->fParams[0]);
      UpdateIfParamAssiged(e, pDecl->fParams[1]);
      UpdateIfParamAssiged(Compress, pDecl->fParams[2]);
   }

   return FMetaGridParams_LogE{R, e, Compress};
}



// Upcoming: Treutler/Ahlrichs-like grids.

FSchemeId FRadialGridScheme_TA::Id() const { return FSchemeId("ta"); }

void FRadialGridScheme_TA::MakeFixedGrid(double *pRadii, double *pWeights, size_t nGridPt, FAtomSpec const &At, FAtomRadialGridDecl const *pDecl) {
   _MakeRadialGrid_TA(&pRadii[0], &pWeights[0], nGridPt, MakeDefaultAtomParams(At, pDecl));
}


FMetaGridParams_TA FRadialGridScheme_TA::MakeDefaultAtomParams(FAtomSpec const &At, FAtomRadialGridDecl const *pDecl) {
   FMetaGridParams_TA
      p{GetGridCenter_TaEta(At.iElement), 0.6};
   if (pDecl != 0) {
      assert(pDecl->fParams.size() >= 2);
      UpdateIfParamAssiged(p.GridCenterR, pDecl->fParams[0]);
      UpdateIfParamAssiged(p.Alpha, pDecl->fParams[1]);
   }
   return p;
}

size_t FRadialGridScheme_TA::EstimateSizeForTargetAccuracy(double fTargetAccuracy, FAtomSpec const &Atom, FAtomRadialGridDecl const *pDecl) {
   return _ComputeRadialGridSize_OriginalTA(fTargetAccuracy, Atom.iElement);
   (void)pDecl; // suppress unused warning
}


FSchemeId FRadialGridScheme_LKO::Id() const { return FSchemeId("lko"); }

size_t FRadialGridScheme_LKO::EstimateSizeForTargetAccuracy(double fTargetAccuracy, FAtomSpec const &Atom, FAtomRadialGridDecl const *pDecl)
{
   // Ochsenfeld 2018's grid offsets (Tab III/eq A6)
   size_t
      nRadialPtOff = 30; // Cs -- open ended.
   //       if (Atom.iElement <= 86) nRadialPtOff = 30; // < Cs -- Rb-Xe
   if (Atom.iElement <= 54) nRadialPtOff = 25; // Rb-Xe
   if (Atom.iElement <= 36) nRadialPtOff = 20; // K-Kr
   if (Atom.iElement <= 18) nRadialPtOff = 10; // Na-Ar
   if (Atom.iElement <= 10) nRadialPtOff = 5; // Li-Ne
   if (Atom.iElement <= 2) nRadialPtOff = 0; // H-He
   // gridsize 1 2 3 4 5  6  7 8 9
   // radsize  1 2 3 6 8 10 14 9 3
   // H-He: 20; Li-Ne: 25; Na-Ar: 30; K-Kr: 40; Rb-Xe: 45; Cs-Lw: 50;
   size_t
      nRadialPt = 30;
   int iGridLevel = iTaGridLevel(fTargetAccuracy);
   if (iGridLevel >= 1) nRadialPt = 35;
   if (iGridLevel >= 2) nRadialPt = 40;
   if (iGridLevel >= 3) nRadialPt = 50;
   if (iGridLevel >= 4) nRadialPt = 55;
   if (iGridLevel >= 5) nRadialPt = 60;
   if (iGridLevel >= 6) nRadialPt = 70;
   if (iGridLevel >= 7) nRadialPt = 80;
   if (iGridLevel >= 8) nRadialPt = 80 + size_t(10 * (dTaGridLevel(fTargetAccuracy) - 7.));

   nRadialPt += nRadialPtOff;
   return nRadialPt;
   (void)pDecl; // suppress unused warning
}






bool check_if_finite_if_supported(FScalar d)
{
#ifdef IR_HAVE_CXX11
   return std::isfinite(d); // <- that's a >= C++11 function.
#else
   return true;
#endif
}


} // namespace mig




// TODO:
// + Implementation separation for concrete grid builder machinery:
//   - Technically, we could move everything FRadialGridBuilder-related into the .cpp
//     file, and only provide a single interface routine
//
//     FRadialGridBuilderPtr MakeRadialGridBuilder(std::string Desc);
//
//     ...which would decide on the concrete class to instanciate and build it
//     (including init options).
//   - We could still provide explicit interfaces to concrete low level grid
//     building routines (with externally supplied parameters) independently
//     of the grid builder classes.
//   - Might need some idea of checking if a concrete grid builder exists
//     and can deal with the supplied inputs already on FDftGridParams assembly
//     stage.
//
// + Interface for integration rule parameters
//   - Should these be part of the builder classes? Or part of the option classes?
//     + Putting them as part of the builders would allow more complex init
//       setups (possibly even information from external sources).
//     + But do we really want that?
//     + It would make it harder to share option sets between integration rules
//     + I guess in-class would be better for the start.
//       (do that in FParamsType constructor: pass it an iElement and a pointer to the grid generator)
//   - In case of in-buider class: we'd probably want separate functions for:
//     + assign element default params FParamsType &Opt
//     + update default parameters from parsing of option desc.
//       parse options(FParamsType &Opt) (on top of default? would mean we could not raise
//       exceptions if default options were not available...)
//   - In either case, there should be a function in TRadialGridBuilder
//     to retrieve the parameter object for a concrete atom and element.
//   - The parameter objects should be constructed only if requested by the rule,
//     or to make default objects for the specific element if overrides are specified in
//     the input.




