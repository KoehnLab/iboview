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

#include <iostream>
#include <cmath>
#include <limits> // for defining DefaultSize and AdaptiveSize
#include <stdlib.h> // for atoi

#include "CtDftGrid_Params.h"
#include "CxParse1.h"
#include "CxIo.h" // for FInputError

#include "CxPhysicalUnits.h"
#include "CxAtomData.h" // for GetCovalentRadius(int iElement)

namespace mig {

using ct::FInputError;
using ct::string_slice;
using ct::split_result;
using ct::FPropertyListStr;
using ct::FPropertyStr;


FGridDeclBase::FScalarParam const FGridDeclBase::NotSpecified = std::numeric_limits<FScalarParam>::max();


FDftGridParams::FDftGridParams(int iGridLevel_)
//    : m_RadialGridDecls(ct::string_slice(), FAtomRadialGridDecl())
{
   InitDefaults();
   SetParamsViaLevel(iGridLevel_);
}


FDftGridParams::FDftGridParams(std::string const &sGridDesc_)
{
   InitDefaults();
   SetParamsViaDesc(sGridDesc_);
}


void FDftGridParams::InitDefaults()
{
//    m_RadialGridScheme = RADSCHEME_LogE;
   std::string sRadialDefault;
//    sRadialDefault = "loge{}";
   sRadialDefault = "ta{}";
//    sRadialDefault = "loge{auto;3.;0.25}";
//    sRadialDefault = "loge{auto;auto;0.0}";
//    sRadialDefault = "lko{}";
//    m_EvalFreeAtomGridAccuracy = true;
   m_EvalFreeAtomGridAccuracy = false;
   m_RadialGridDecls.SetDefault(FAtomRadialGridDecl(string_slice(sRadialDefault)));
   m_pFitBasis = 0;
   m_AdjustVoronoiCellsByAtomicRadii = true;
   m_fWeightCut = -1; // auto.

//    m_fWeightCut = 0; std::cout << "WARNING: wiping out weight cut defaults\n";
   m_fWeightCut_AtomVdwRadiusFactor = -1; // auto
   for (size_t iPeriod = 0; iPeriod < DFTGRID_nPeriodParams; ++ iPeriod)
      m_PeriodGridDefs[iPeriod] = FPeriodGridDef();
//    SetTargetAccuracy(1e-5);
   SetTargetAccuracy(1e-4);

//    m_VoronoiNeighborScreenType = VORONOI_ScreenByDistance;
//    m_VoronoiNeighborScreenType = VORONOI_KeepAll;
   m_VoronoiNeighborScreenType = VORONOI_ScreenByWeightOfClosestPoint;
   m_VoronoiNeighborScreenThresholdFar = 0.001;
   m_VoronoiNeighborScreenThresholdNear = 0.; // disabled.
//    m_VoronoiNeighborScreenThresholdNear = m_VoronoiNeighborScreenThresholdFar;
   m_CoreRadiusThreshold = -1.0;
//    m_PrintLevel = 2;
   m_PrintLevel = 0;
   m_AngularGridAlign = ANGULARALIGN_AtomicEnvironment;
   m_GridRenormType = GRIDRENORM_None;
}


double FDftGridParams::GetCoreRadius(int iElement) const
{
   double
      CoreRadiusThreshold = m_CoreRadiusThreshold;
   if (CoreRadiusThreshold == -1.) {
      // decide by ourselves, based on target accuracy.
      CoreRadiusThreshold = 0.05;
      if (m_fTargetAccuracy >= 1e-7)
         CoreRadiusThreshold = 0.1;
      if (m_fTargetAccuracy >= 1e-6)
         CoreRadiusThreshold = 0.2;
      if (m_fTargetAccuracy >= 1e-5)
         CoreRadiusThreshold = 0.3;
//       if (iElement > 10)
//          CoreRadiusThreshold *= 0.5;
   }

   return ct::GetCovalentRadius(iElement) * CoreRadiusThreshold;
}


double FDftGridParams::fWeightCut_AtomVdwRadiusFactor() const {
   // something explicitly set?
   if (m_fWeightCut_AtomVdwRadiusFactor != -1) {
      return m_fWeightCut_AtomVdwRadiusFactor;
   } else {
      // -1 means "we decide".
      if (m_fWeightCut == 0)
         // main weight cut completely disabled. That probably means we should not cut weights, either.
         return 0;
      if (m_fTargetAccuracy >= 3e-6)
         return 2.00; // should be sufficient normally.
      else if (m_fTargetAccuracy >= 3e-7)
         return 2.50; // should be sufficient normally.
      else if (m_fTargetAccuracy >= 1e-8)
         return 3.00;
      return 0;
   }

}

double FDftGridParams::fTargetAccuracy() const {
//    // would return 1e-5 for level 3 and 1e-6 for level 5
//    return 1e-3 * std::pow(.1, double(this->nLevel+1)/2.);
   return m_fTargetAccuracy;
}


void FDftGridParams::SetTargetAccuracy(double fAcc) {
// //    fAcc = 1e-3 * std::pow(.1, double(this->nLevel+1)/2.);
// //    1e3 * fAcc = std::pow(.1, double(this->nLevel+1)/2.);
// //    std::log(1e3 * fAcc)/std::log(.1) = double(this->nLevel+1)/2.
// //    2 * std::log(1e3 * fAcc)/std::log(.1) = double(this->nLevel+1)
// //    2 * std::log(1e3 * fAcc)/std::log(.1) - 1 = this->nLevel
//    nLevel = 2. * (std::log(1e3 * fAcc)/std::log(.1)) - 1.;
//    assert_rt(std::abs(fAcc - fTargetAccuracy()) < 1e-10);
   m_fTargetAccuracy = fAcc;
}


void FDftGridParams::SetParamsViaLevel(int iGridLevel)
{
   // would return 1e-5 for level 3 and 1e-6 for level 5
   SetTargetAccuracy(1e-3 * std::pow(.1, double(iGridLevel+1)/2.));
}

double FDftGridParams::GetApproxLevel() const
{
   // a = 1e-3 * std::pow(.1, double(iGridLevel+1)/2.)
   // a/(1e-3) = std::pow(.1, double(iGridLevel+1)/2.)
   // -log10(a/1e-3) = (iGridLevel+1)/2
   // 2*(-log10(a/1e-3))-1 = iGridLevel
   return 2.*(-std::log10(m_fTargetAccuracy/1e-3)) - 1.;
}


void FDftGridParams::SetParamsViaDesc(std::string const &sGridDesc_)
{
   // check if desc is simply an integer, like '5'. In this case, we set the default parameters via the grid level.
   double
      fNumericDesc = ::atof(sGridDesc_.c_str()); // returns 0 for invalid input.
   if (fNumericDesc != 0.) {
      // is it an integer? if yes -> most likely a level spec.
      if ((fNumericDesc - int(fNumericDesc)) == 0)
         return SetParamsViaLevel(int(fNumericDesc));
      else {
         // nop. propably a grid accuracy spec.
         if (fNumericDesc < 0 || fNumericDesc > 0.1)
            throw std::runtime_error(fmt::format("failed to understand numeric grid specification. expected positive integers (grid accuracy levels) or floats in range 0 <= f <= 0.1, but got '{}'", sGridDesc_));
         return SetTargetAccuracy(fNumericDesc);
      }
   }
   // ...apparently not. Parse as a property list.
   FPropertyListStr
      GridProps = FPropertyListStr(string_slice(sGridDesc_));
   for (size_t iOpt = 0; iOpt != GridProps.size(); ++ iOpt) {
      FPropertyStr const
         &Prop = GridProps[iOpt];
      if (Prop.Name == "accu") {
         SetTargetAccuracy(Prop.Content.to_float());
      } else if (Prop.Name == "level") {
         SetParamsViaLevel(Prop.Content.to_int());
      } else if (Prop.Name == "weight-cut" || Prop.Name == "weightcut") {
         m_fWeightCut = Prop.Content.to_float();
      } else if (Prop.Name == "weight-cut-vdw-radius-factor" || Prop.Name == "weight-cut-vdw") {
         // specified as multiple of the iso-density atomic vdW radius (see GetVdwRadius_IsoDensity).
         // Note: the weight is not actually cut there---that is just where the 0.5-weight step center is placed.
         m_fWeightCut_AtomVdwRadiusFactor = Prop.Content.to_float();
      } else if (Prop.Name == "print-level" || Prop.Name=="print") {
         m_PrintLevel = Prop.Content.to_int();
      } else if (Prop.Name == "radial") {
         m_RadialGridDecls.Update(Prop.Content);
//       } else if (Prop.Name == "radial-scheme") {
//          if (Prop.Content == "ta" || Prop.Content == "ta95" || Prop.Content == "ta1995" || Prop.Content == "ahlrichs" || Prop.Content == "treutler")
//             m_RadialGridScheme = RADSCHEME_Ta1995;
//          else if (Prop.Content == "lko" || Prop.Content == "ochsenfeld" || Prop.Content == "laqua")
//             m_RadialGridScheme = RADSCHEME_Lko2018;
//          else if (Prop.Content == "loge")
//             m_RadialGridScheme = RADSCHEME_LogE;
//          else if (Prop.Content == "logb" || Prop.Content == "logbeta")
//             m_RadialGridScheme = RADSCHEME_LogBeta;
//          else
//             throw FInputError(fmt::format("invalid radial integration scheme: '{}'. 'radial' must be 'loge', 'ta', or 'lko'.", Prop.Content.to_str()));
      } else if (Prop.Name == "align") {
         if (Prop.Content == "none")
            m_AngularGridAlign = ANGULARALIGN_None;
         else if (Prop.Content == "local")
            m_AngularGridAlign = ANGULARALIGN_AtomicEnvironment;
         else if (Prop.Content == "random" || Prop.Content == "randomize")
            m_AngularGridAlign = ANGULARALIGN_Randomize;
         else
            throw FInputError(fmt::format("invalid angular grid alignment: '{}'. 'align' must be 'local', 'none', or 'random'.", Prop.Content.to_str()));
      } else if (Prop.Name == "renorm") {
         if (Prop.Content == "none")
            m_GridRenormType = GRIDRENORM_None;
         else if (Prop.Content == "orb" || Prop.Content == "orbs" || Prop.Content == "orbitals" || Prop.Content == "yes" || Prop.Content == "1" || Prop.Content == "on")
            m_GridRenormType = GRIDRENORM_Orbitals;
         else if (Prop.Content == "aux" || Prop.Content == "dens" || Prop.Content == "density" || Prop.Content == "aux-density")
            m_GridRenormType = GRIDRENORM_AuxDensity;
         else
            throw FInputError(fmt::format("invalid grid renormalization type: '{}'. 'renorm' must be 'none', 'orbitals', or 'aux'.", Prop.Content.to_str()));
      } else if (Prop.Name == "core-radius") {
         m_CoreRadiusThreshold = Prop.Content.to_float();
         if (m_CoreRadiusThreshold != -1.) {
            // note: '-1' means 'use defaults'.
            if (m_CoreRadiusThreshold < 0. || m_CoreRadiusThreshold > 1.0)
               throw FInputError(fmt::format("invalid core-radius: '{}'. 'core-radius' is given in fractions of an element's covalent radius. The value must either be '-1.0' (for using default settings) or lie in range 0.0 <= core-radius <= 1.0.", Prop.Content.to_str()));
         }
      } else if (Prop.Name == "voronoi-neighbors") {
         if (Prop.Content == "keep-all" || Prop.Content == "all") {
            m_VoronoiNeighborScreenType = FDftGridParams::VORONOI_KeepAll;
         } else if (Prop.Content == "screen-by-distance") {
            m_VoronoiNeighborScreenType = FDftGridParams::VORONOI_ScreenByDistance;
         } else if (Prop.Content.startswith("screen-by-weight")) {
            m_VoronoiNeighborScreenType = FDftGridParams::VORONOI_ScreenByWeightOfClosestPoint;
            if (Prop.Content.endswith("}")) {
               FPropertyListStr
                  ScreenProps(Prop.Content);
               for (size_t iScreenOpt = 0; iScreenOpt != ScreenProps.size(); ++ iScreenOpt) {
                  FPropertyStr const
                     &Prop = ScreenProps[iScreenOpt];
                  if (Prop.Name == "thr") {
                     m_VoronoiNeighborScreenThresholdFar = Prop.Content.to_float();
                     m_VoronoiNeighborScreenThresholdNear = m_VoronoiNeighborScreenThresholdFar;
                  } else if (Prop.Name == "thr-far") {
                     m_VoronoiNeighborScreenThresholdFar = Prop.Content.to_float();
                  } else if (Prop.Name == "thr-near") {
                     m_VoronoiNeighborScreenThresholdNear = Prop.Content.to_float();
                  } else {
                     throw FInputError(fmt::format("screen-by-weight: unrecognized option '{}'. Known properties: [thr, thr-far, thr-near].",
                        Prop.Name.to_str()));
                  }
               }
            }
         } else
            throw FInputError(fmt::format("grid-desc: unrecognized option '{}' for parameter 'voronoi-neighbors'. Known properties: ['keep-all', 'screen-by-distance', 'screen-by-weight']", Prop.Content.to_str()));
      } else if (Prop.Name == "lmin" || Prop.Name == "lmax" || Prop.Name == "lset" || Prop.Name == "nradial") {
         split_result sr;
         Prop.Content.split_list(sr);
         for (size_t iPeriod = 0; iPeriod < sr.size(); ++ iPeriod) {
            if (iPeriod >= DFTGRID_nPeriodParams)
               throw std::runtime_error(fmt::format("grid-desc: {} ... too many arguments ({}). At most {} per-period arguments supported.", Prop.Name.to_str(), sr.size(), int(DFTGRID_nPeriodParams)));
            if (Prop.Name == "lmin")
               m_PeriodGridDefs[iPeriod].iMinL = sr[iPeriod].to_int();
            else if (Prop.Name == "lmax")
               m_PeriodGridDefs[iPeriod].iMaxL = sr[iPeriod].to_int();
            else if (Prop.Name == "lset") {
               ptrdiff_t l = sr[iPeriod].to_int();
               m_PeriodGridDefs[iPeriod].iMinL = l;
               m_PeriodGridDefs[iPeriod].iMaxL = l;
            } else if (Prop.Name == "nradial")
               m_PeriodGridDefs[iPeriod].nRadialPt = sr[iPeriod].to_int();
            else {
               assert(0);
            }
         }
      } else {
         throw FInputError(fmt::format("grid-desc: unrecognized option '{}'. Known properties:"
            "\n  accu: float (target grid accuracy, e.g., 1e-5),"
            "\n  level: int (target grid accuracy level, e.g. 5),"
            "\n  lmin: [n,n,n,n] (minimum L of angular grid, for elements of period 1-4)"
            "\n  lmax: [n,n,n,n] (maximum L of angular grid, for elements of period 1-4)"
            "\n  lset: [n,n,n,n] (set L of angular grid, for elements of period 1-4)"
            "\n  nradial: [n,n,n,n] (number of radial grid points, for elements of period 1-4)",
            "\n  weight-cut: float (discard grid points below this weight; dangerous!)",
            Prop.Name.to_str()));
      }
   }
}


double FDftGridParams::dGridLevel() const
{
   return 1. + 2.*(-std::log10(this->fTargetAccuracy()) - 4.);
   // ^-'grid 1'. That has target acc 1e-4, 'grid 3' has target acc 1e-5.
}


int FDftGridParams::iGridLevel() const
{
   return int(std::round(this->dGridLevel()));
}


double FDftGridParams::fWeightCut() const
{
   // non-default params set?
   if (m_fWeightCut != -1)
      // return whatever was asked for.
      return m_fWeightCut;
   // return a semi-sensible default (note: this DOES have a substantial effect on the accuracy;
   // does not seem to affect optg performance, however, and can substantially reduce number of points).
//    return std::min(1e-12, 1e-6*m_fTargetAccuracy);
   return std::min(1e-12, 1e-3*sqr(m_fTargetAccuracy));
}


int FDftGridParams::iMinL(int iPeriod) const
{
   if (iPeriod < 1 || iPeriod > DFTGRID_nPeriodParams)
      return -1;
   return m_PeriodGridDefs[iPeriod-1].iMinL;
}


int FDftGridParams::iMaxL(int iPeriod) const
{
   if (iPeriod < 1 || iPeriod > DFTGRID_nPeriodParams)
      return -1;
   return m_PeriodGridDefs[iPeriod-1].iMaxL;
}


bool FDftGridParams::HasFixedL(int iPeriod) const
{
   if (iPeriod < 1 || iPeriod > DFTGRID_nPeriodParams)
      return false;
   return (m_PeriodGridDefs[iPeriod-1].iMaxL == m_PeriodGridDefs[iPeriod-1].iMinL) && (m_PeriodGridDefs[iPeriod-1].iMaxL >= 1);
}



int FDftGridParams::nRadialPt(int iPeriod) const
{
   if (iPeriod < 1 || iPeriod > DFTGRID_nPeriodParams)
      return -1;
   return m_PeriodGridDefs[iPeriod-1].nRadialPt;
}




void FSchemeId::_Init1(char const *pName, size_t NameSize, unsigned Flags)
{
   this->fill(0);
   // was a name given? (call from default ctor comes with pName = 0)
   if (pName) {
      auto IsNameIndexValid = [=](size_t i) -> bool { return (NameSize == std::string::npos)? (pName[i] != 0) : (i < NameSize); };
      size_t i = 0;
      for (; i < this->size() && IsNameIndexValid(i); ++i)
         (*this)[i] = pName[i];
      // still some part of the input name left which we did not have the size to copy?
      if (!bool(Flags & AllowTruncation) && IsNameIndexValid(i))
         throw std::runtime_error(fmt::format("FSchemeId()::ctor: Scheme name '{}' exceeds fixed size limit {}.", pName, this->size()));
   }
}


FSchemeId::operator std::string () const
{
   size_t
      nChars = 0;
   while (nChars < this->size() && (*this)[nChars] != 0)
      nChars += 1;
   return std::string(&(*this)[0], nChars);
}


bool FSchemeId::IsAssigned() const
{
   for (char const &v : *this) {
      if (v != 0)
         return true;
   }
   return false;
}



FAtomRadialGridDecl::FAtomRadialGridDecl(std::string::const_iterator itFirst, std::string::const_iterator itLast)
{
   _InitDefaults();
   Update(itFirst, itLast);
}


FAtomRadialGridDecl::FAtomRadialGridDecl(ct::string_slice slDecl)
{
   _InitDefaults();
   Update(slDecl.first, slDecl.last);
}


void FAtomRadialGridDecl::_InitDefaults()
{
   nRadialPt = DefaultSize;
//    nAngularOrder = DefaultSize;
   fTargetAccuracy = NotSpecified;
   fParams.fill(NotSpecified);
}


void FAtomRadialGridDecl::Update(FAtomRadialGridDecl const &other)
{
   // hm.. i could probably automate this with tie()... but i think I'll just go with
   // the explicit variant for now...
   if (other.SchemeId.IsAssigned())
      SchemeId = other.SchemeId;
   UpdateIfParamAssiged(this->nRadialPt, other.nRadialPt);
   UpdateIfParamAssiged(this->fTargetAccuracy, other.fTargetAccuracy);
   for (size_t ip = 0; ip < std::min(fParams.size(), other.fParams.size()); ++ ip) {
      UpdateIfParamAssiged(this->fParams[ip], other.fParams[ip]);
   }
}


static bool IsExplicitDefaultValue(string_slice const &sl) {
   return sl == "default" || sl == "auto" || sl == "none" || sl == "null";
}


static double fGridParam(string_slice const &sl) {
   if (sl.empty() || IsExplicitDefaultValue(sl))
      return FGridDeclBase::NotSpecified;
   return sl.to_float();
}

static unsigned nGridSizeParam(string_slice const &sl) {
   if (sl.empty() || IsExplicitDefaultValue(sl))
      return FGridDeclBase::DefaultSize;
   if (sl.startswith("adap") || sl == "auto")
      return FGridDeclBase::AdaptiveSize;
   return unsigned(sl.to_size());
}


void FAtomRadialGridDecl::Update(std::string::const_iterator itFirst, std::string::const_iterator itLast) {
   string_slice
      slFullDecl(itFirst, itLast);
   FPropertyListStr
      GridProps = FPropertyListStr(slFullDecl, 0, ct::PROPLIST_AllowEntriesWithoutName | ct::PROPLIST_IdentifierWithOptionalArgs);
   // if name is empty, leave alone whatever we currently have.
   // But if the scheme is explicitly declared as default, then clear the
   // scheme name.
   if (!GridProps.Name.empty()) {
      if (IsExplicitDefaultValue(GridProps.Name))
         this->SchemeId = FSchemeId();
      else
         this->SchemeId = FSchemeId(GridProps.Name.to_str());
   }
//    std::cout << fmt::format("!radial-grid-decl::update({}) s = '{}'", std::string(this->SchemeId), slFullDecl.to_str()) << std::endl;
   size_t
      iUnnamedArg = 0;
   for (size_t iOpt = 0; iOpt != GridProps.size(); ++ iOpt) {
      FPropertyStr const
         &Prop = GridProps[iOpt];
      if (Prop.Name.empty()) {
         if (iUnnamedArg >= fParams.size())
            throw FInputError(fmt::format("radial-grid-decl: too many numerical parameters in grid decl '{}'."
               " Only up to {} unnamed parameters are allowed.",
               slFullDecl.to_str(), fParams.size()));
         if (!Prop.Content.empty()) {
            // ^--- allow completely empty entries --- in this case leave current parameter values alone.
            if (IsExplicitDefaultValue(Prop.Content))
               fParams[iUnnamedArg] = NotSpecified;
            else
               fParams[iUnnamedArg] = fGridParam(Prop.Content);
         }
         iUnnamedArg += 1;
      } else if (Prop.Name == "accu") {
         fTargetAccuracy = fGridParam(Prop.Content);
      } else if (Prop.Name == "npt" || Prop.Name == "npts" || Prop.Name == "n" || Prop.Name == "size" || Prop.Name == "num-points") {
         nRadialPt = nGridSizeParam(Prop.Content);
//       } else if (Prop.Name == "l") {
//          nAngularOrder = nGridSizeParam(Prop.Content);
      } else {
         throw FInputError(fmt::format("radial-grid-decl: unrecognized option '{}'.",
            Prop.Name.to_str()));
      }
   }

}


std::string FAtomRadialGridDecl::Format() const
{
   auto FormatSize = [](FGridDeclBase::FSizeParam n) -> std::string {
      if (n == FGridDeclBase::DefaultSize) return "def";
      if (n == FGridDeclBase::AdaptiveSize) return "adaptive";
      return fmt::format("{}", n);
   };
   auto FormatScalar = [](FGridDeclBase::FScalarParam a, char const *pFmtStr = "{}") -> std::string {
      if (a == FGridDeclBase::NotSpecified) return "_";
      return fmt::format(pFmtStr, a);
   };
   std::stringstream
      str;
   for (size_t ip = 0; ip != fParams.size(); ++ ip) {
      if (ip != 0) str << ",";
      str << FormatScalar(fParams[ip]);
   }
   return fmt::format("{}{{n:{}; acc:{}; p:[{}]}}", static_cast<std::string>(SchemeId),
            FormatSize(nRadialPt), FormatScalar(fTargetAccuracy,"{:.2e}"), str.str());
}



inline bool operator < (FAtomRadialGridDecl const &A, FAtomRadialGridDecl const &B) {
   return A.SortKey() < B.SortKey();
}




} // namespace ct
