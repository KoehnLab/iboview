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

#include <stdexcept>
#include "CxParse1.h"
#include "CtRhfOptions.h"
#include "CxPhysicalUnits.h"
#include "CtDftFunc.h" // for checking which functionals require exact exchange.

namespace ct {



FHfOptions::FHfOptions(double ThrOrb_, double ThrDen_, unsigned MaxIt_)
   : ThrOrb(ThrOrb_), ThrDen(ThrDen_), MaxIt(MaxIt_),
      XcAlgo(XCALGO_Regular), JkAlgo(JKALGO_DfAuto), OrbType(ORBTYPE_Restricted), ComputeFlags(0),
      XcFunctionalName(""), DispCorrType(""),
      UseRefineGrid(false),
      IgnoreConvergenceError(false),
      UseBandStructureEnergyFormula(false)
{
//    LevelShifts[0] = 0.0;
//    LevelShifts[1] = 0.0;
   LevelShifts[0] = -0.4;
   LevelShifts[1] = -0.2;
   InitialDampingLevelShift = 0.;
   InitialDampingLevelShiftDecay = 0.5;
   MaxDiis = 10;
//    LevelShifts[0] = -1.0;
//    LevelShifts[1] = -0.5;
   m_RohfAlgo = ROHFALGO_OneStep;
//    m_RohfAlgo = ROHFALGO_TwoStep;
// ^- gradient broken for two step algo. And convergence is quite bad... for some reason.
   iTimingLevel = -1;
}



bool UsingDensityFitting(FJkAlgo JkAlgo) {
   return JkAlgo != JKALGO_4ixDirect;
}

// will return 'R', 'U', or 'G'.
char const *GetOrbTypePrefix(FScfOrbType OrbType)
{
   if (OrbType == ORBTYPE_Restricted) {
      return "R";
   } else if (OrbType == ORBTYPE_Unrestricted) {
      return "U";
   } else if (OrbType == ORBTYPE_Generalized) {
      return "G";
   } else {
      throw std::runtime_error(fmt::format("GetOrbTypePrefix: orbital type '{}' not recognized.", int(OrbType)));
   }
}


// returns 'true' unless JkAlgo is set to a 'pure' (Coulomb-only) density functional method.
bool FHfOptions::UseExactExchange() const
{
   if (JkAlgo == JKALGO_DfCoulombOnlyCache3ix)
      return false;
   if (JkAlgo == JKALGO_DfCoulombOnly)
      return false;
   if (JkAlgo == JKALGO_4ixDirect ||
       JkAlgo == JKALGO_DfCoulombAndExchange ||
       JkAlgo == JKALGO_DfCoulombChainOfSpheresExchange)
      return true;
   if (JkAlgo == JKALGO_DfAuto) {
      if (XcFunctionalName.empty())
         return true; // that's Hartree-Fock.
      // did we check for this particular functional already if it
      // required exact exchange?
      if (m_LastXcFuncName == XcFunctionalName)
         return m_LastXcFuncWasHybrid;
      // well...still nothing set. Have to instanciate the functional to see if
      // would require exact exchange. Maybe a bit of an overkill, but I guess
      // It's okay for the moment... and we do remember the result.
      FXcFunctional
         XcFn(XcFunctionalName);
      m_LastXcFuncName = XcFunctionalName;
      m_LastXcFuncWasHybrid = XcFn.NeedsExactExch();
      return m_LastXcFuncWasHybrid;
   }

   return true;
}


bool FHfOptions::UseXc() const
{
   return !(XcFunctionalName.empty() || XcFunctionalName == "none");
}


// make a short description of what the fitting basis is used for:
//   j: coulomb fitting (always),
//   x: exchange-correlation fitting (for DFT, i.e., DFXC=1)
//   k: exact-exchange fitting
char const *GetFitBasisType(FHfOptions const *pHfOptions, bool UpperCase, bool *pUseJkFit)
{
   bool
      UseFittedK,
      UseDfxc = (pHfOptions->XcAlgo == XCALGO_AuxiliaryExpand);
   UseFittedK = pHfOptions->UseExactExchange();
   if (pHfOptions->JkAlgo == JKALGO_4ixDirect || pHfOptions->JkAlgo == JKALGO_DfCoulombChainOfSpheresExchange)
      // COSX and 4ix make exact exchange, but not via a fitting basis
      UseFittedK = false;
   if (pUseJkFit != 0)
      *pUseJkFit = UseFittedK;
   if (!UsingDensityFitting(pHfOptions->JkAlgo)) {
      if (UseDfxc)
         throw std::runtime_error("DFXC=1 is not allowed when using non-density fitted exchange-correlation JK algorithms.");
      return 0;
   }
   char const
      *pFitType = 0;
   if (UseFittedK && UseDfxc)
      pFitType = (UpperCase? "JKXFIT" : "jkxfit");
   else if (UseFittedK && !UseDfxc)
      pFitType = (UpperCase? "JKFIT" : "jkfit");
   else if (!UseFittedK && UseDfxc)
      pFitType = (UpperCase? "JXFIT" : "jxfit");
   else {
      assert(!UseFittedK && !UseDfxc);
      pFitType = (UpperCase? "JFIT" : "jfit");
   }
   assert(pFitType != 0);
   return pFitType;
}


std::string GetMethodWfDesc(FHfOptions const *pHfOptions, FWfDecl const *pWfDecl)
{
   fmt::MemoryWriter w;
   bool first = true;
   if (pHfOptions) {
      if (UsingDensityFitting(pHfOptions->JkAlgo))
         w << "DF-";
      if (pHfOptions->XcFunctionalName == "") {
         w << GetOrbTypePrefix(pHfOptions->OrbType) << "HF";
         if (pHfOptions->DispCorrType != "" && pHfOptions->DispCorrType != "none")
            w.write("+{}", pHfOptions->DispCorrType);
      } else {
//       else if (pHfOptions->JkAlgo == JKALGO_DfCoulombOnly || pHfOptions->JkAlgo == JKALGO_DfCoulombOnlyCache3ix)
         w.write("{}KS/{}", GetOrbTypePrefix(pHfOptions->OrbType), pHfOptions->XcFunctionalName);
         if (pHfOptions->DispCorrType != "" && pHfOptions->DispCorrType != "none")
            w.write("+{}", pHfOptions->DispCorrType);
      }
//       else
//          w.write("(unknown-scf)");
//       if (pHfOptions->XcFunctionalName != "")
//          w.write("/%s", pHfOptions->XcFunctionalName);
      if (!pHfOptions->BasisName_Orbital.empty())
         w.write("/{}", pHfOptions->BasisName_Orbital);

      bool
         UseJkFit;
      char const
         *pFitType = GetFitBasisType(pHfOptions, false, &UseJkFit);
      if (pFitType != 0) { // <- 0 is returned if the calculation does not use fitting.
         w.write(" {}", pFitType);

         if (!UseJkFit && !pHfOptions->BasisName_JFit.empty())
            w.write("={}", pHfOptions->BasisName_JFit);
         else if (UseJkFit && !pHfOptions->BasisName_JkFit.empty())
            w.write("={}", pHfOptions->BasisName_JkFit);
      }

      if (pHfOptions->XcFunctionalName != "") {
         if (pHfOptions->XcAlgo != XCALGO_Gridless) {
//             if (UseDfxc)
//                w.write(" dfxc={}", 1);
            // ^- would be covered by the fitting basis desc
            if (!pHfOptions->DftGridDesc.empty())
               w.write(" grid={}",pHfOptions->DftGridDesc);
         } else {
            w << " gridless";
         }
      }
      first = false;
   }
   if (pWfDecl) {
      if (pWfDecl->iCharge != 0) {
         if (first) { first = false; w << " "; }
         w.write(" charge={}", pWfDecl->iCharge);
      }
      if (pWfDecl->Ms2 != 0) {
         if (first) { first = false; w << " "; }
         w.write(" spin={}", pWfDecl->Ms2);
      }
   }
   return w.str();
}


std::string FHfOptions::GetMethodDesc() const
{
   return GetMethodWfDesc(this, 0);
}


void FHfOptions::SetGridDesc(std::string const &NewGridDesc_1)
{
   this->DftGridDesc = NewGridDesc_1;
   std::string DftGridDesc_ = NewGridDesc_1;

   if (!DftGridDesc_.empty()) {
      if (DftGridDesc_[0] == 'm') {
         DftGridDesc_.erase(0, 1);
         UseRefineGrid = true;
      } else
         UseRefineGrid = false;
   }

   if (!DftGridDesc_.empty()) {
      try {
         DftGridParams = FDftGridParams(DftGridDesc_);
      } catch (FInputError &e) {
         throw FInputError(fmt::format("Failed to understand grid declaration '{}'.\nExamples of valid descs: '3', '5', 'm5', '{{accu: 1e-5; lmin:[15,15,15,15]}}'. Exception was:\n'{}'", DftGridDesc_, e.what()));
      }
   }

   if (UseRefineGrid) {
      DftGridParamsRefine = DftGridParams;
//       DftGridParams.SetTargetAccuracy(1e2 * DftGridParamsRefine.fTargetAccuracy());
//       DftGridParams.SetTargetAccuracy(3e1 * DftGridParamsRefine.fTargetAccuracy());
      DftGridParams.SetTargetAccuracy(1e1 * DftGridParamsRefine.fTargetAccuracy());
   }
}


void FHfOptions::SetXc(std::string const &XcName_)
{
   // most functionals are canonically written in all-upper case...
   std::string
      XcNameUpper = toupper(XcName_);
   string_slice
      slXcName(XcName_);
   slXcName.trim();
   {
      // check if there is a dispersion correction part (technically not really
      // part of a XC functional, but commonly written like that)
      cstr_it
         iDispPos = slXcName.rfind("+D");
      if (iDispPos != slXcName.end()) {
         // seems like it. extract it and set it as disp corr parameter
         // (at the time of this writing we really only have D3(BJ), actually)
         string_slice
            slDispCorr(iDispPos+1, slXcName.last);
         this->DispCorrType = slDispCorr.to_str();

         // cut out the disp corr part from the actual xc part.
         slXcName.last = iDispPos;
         slXcName.trim_right();
      }
   }


   if (slXcName == "RHF" || slXcName == "UHF" || slXcName == "HF") {
      this->XcFunctionalName = "";
//       JkAlgo = JKALGO_DfCoulombAndExchange;
   } else {
      this->XcFunctionalName = slXcName.to_str();
//       JkAlgo = JKALGO_DfCoulombOnlyCache3ix;
   }
}

void FHfOptions::SetDamping(std::string const &Args_)
{
   string_slice
      slArgs(Args_);
   try {
      bool
         IsList = slArgs.startswith('{');
      if (!IsList) {
         InitialDampingLevelShift = slArgs.to_float();
      } else {
         FPropertyListStr
            Props = FPropertyListStr(slArgs, 0, PROPLIST_AllowEntriesWithoutName);
         size_t
            iUnnamedArg = 0;
         for (size_t iOpt = 0; iOpt != Props.size(); ++ iOpt) {
            FPropertyStr
               Prop = Props[iOpt];
            std::string
               NameOverride;
            if (Prop.Name.empty()) {
               if (iUnnamedArg == 0) {
                  NameOverride = "shift";
               } else if (iUnnamedArg == 1) {
                  NameOverride = "decay";
               } else if (iUnnamedArg == 2) {
                  NameOverride = "factor";
               } else {
                  // ignore. Maybe it's from a future input.
               }
               if (!NameOverride.empty())
                  Prop.Name.assign_ref(NameOverride);
               iUnnamedArg += 1;
            }

            if (Prop.Name == "shift") {
               InitialDampingLevelShift = Prop.Content.to_float();
            } else if (Prop.Name == "decay") {
               InitialDampingLevelShiftDecay = Prop.Content.to_float();
            } else {
               // can't really report problems here, can we? It's a constructor and we have no log.
               // Make an extra function for this?
               throw std::runtime_error(fmt::format("FHfOptions::SetDamping: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
            }
         }
      }
   } catch (std::runtime_error &e) {
      // convert exceptions to input errors.
      throw FInputError(e.what());
   }
}


void FHfOptions::SetDfxc(bool NewDfxc_)
{
   if (NewDfxc_)
      XcAlgo = XCALGO_AuxiliaryExpand;
   else
      XcAlgo = XCALGO_Regular;
}


void FHfPrintOptions::SetFlag(unsigned Flag, bool Value, bool Extra)
{
   if (Value)
      m_FlagsPrintOn |= Flag;
   else
      m_FlagsPrintOn &= ~Flag;

   if (Extra)
      m_FlagsExtraDetail |= Flag;
   else
      m_FlagsExtraDetail &= ~Flag;
}


void FHfPrintOptions::SetFlag(unsigned Flag, string_slice slValue)
{
   ptrdiff_t
      iValue;
   if (slValue.try_convert_to_int(&iValue)) {
      if (iValue <= 0)
         SetFlag(Flag, false, false);
      else if (iValue == 1)
         SetFlag(Flag, true, false);
      else if (iValue >= 1)
         SetFlag(Flag, true, true);
   } else if (slValue == "extra") {
      SetFlag(Flag, true, true);
   } else {
      SetFlag(Flag, slValue.to_bool(), false);
   }
}


void FHfPrintOptions::SetArgs(std::string const &Args_)
{
   std::string
      sOn = "1";
   try {
      FPropertyListStr
         Props = FPropertyListStr(Args_, 0, PROPLIST_AllowEntriesWithoutName);
      for (size_t iOpt = 0; iOpt != Props.size(); ++ iOpt) {
         FPropertyStr
            Prop = Props[iOpt]; // note: this makes a copy, but it's just two string_slice objects (4 ptrs total).
         if (Prop.Name.empty()) {
            // by default assume that if someone puts a print key without
            // arguments, that they want the corresponding thing printed at
            // level on/1.
            Prop.Name = Prop.Content;
            Prop.Content = string_slice(sOn);
         }
         if (Prop.Name == "guess") {
            SetFlag(SCFPRINT_GuessInfo, Prop.Content);
         } else if (Prop.Name == "basis-orb" || Prop.Name == "basis") {
            SetFlag(SCFPRINT_OrbBasis, Prop.Content);
         } else if (Prop.Name == "basis-fit") {
            SetFlag(SCFPRINT_FitBasis, Prop.Content);
         } else if (Prop.Name == "orbitals") {
            SetFlag(SCFPRINT_Orbitals, Prop.Content);
         } else if (Prop.Name == "ecp" || Prop.Name == "ecps") {
            SetFlag(SCFPRINT_Ecp, Prop.Content);
         } else if (Prop.Name == "command" || Prop.Name == "job" || Prop.Name == "method") {
            SetFlag(SCFPRINT_MethodWfDesc, Prop.Content);
         } else if (Prop.Name == "eig" || Prop.Name == "ew" || Prop.Name == "orb-ew" || Prop.Name == "orb-eig") {
            SetFlag(SCFPRINT_OrbitalEnergies, Prop.Content);
         } else {
            // can't really report problems here, can we? It's a constructor and we have no log.
            // Make an extra function for this?
            throw std::runtime_error(fmt::format("FHfPrintOptions: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
         }
      }
   } catch (std::runtime_error &e) {
      // convert exceptions to input errors.
      throw FInputError(e.what());
   }
}


FMpt1e::FMpt1e()
{
   m_PropertyType = PROPTYPE_Unknown;
   m_Order = 0;
   m_Factor = 1.0;
   m_PrintLevel = 1;
   // by default do not put any coordinates. Actual evaluation routines
   // should set default center to coordinate origin if one is needed
   // but was not provided.
}

FMpt1e::~FMpt1e()
{
}


// read a single string-slice containing a list of [x,y,z] coordinates, and
// return it as a FVector3 structure. Components will be multiplied with OutputFactor.
// (e.g., for angstrom/abohr conversion)
FVector3 ParseVector3Arg(string_slice slVectorArg, double OutputFactor)
{
   FVector3
      vResult(0,0,0);
   try {
      short_split_result
         srCoords;
      slVectorArg.split_list(srCoords);
      if (srCoords.size() != 3)
         throw std::runtime_error("incorrect number of coordinates in [x,y,z]");
      for (size_t ixyz = 0; ixyz != 3; ++ ixyz)
         vResult[ixyz] = OutputFactor * srCoords[ixyz].to_float();
   } catch (std::runtime_error &e) {
      throw std::runtime_error(fmt::format("ParseVector3Arg(): expected format is '[<x>,<y>,<z>]' where <x>, <y>, <z> denote the x,y,z coordinates (floats) of the point. Actually received: '{}'; processing exception: '{}'", slVectorArg.to_str(), e.what()));
   }
   return vResult;
}


bool FMpt1e::FindPropertyType(string_slice slType, FPropertyType *pPropertyTypeOut, int *pOrderOut)
{
   FPropertyType
      PropertyType = PROPTYPE_Unknown;
   int
      Order = -1;
   if (slType == "kinetic" || slType == "ekin" || slType == "t") {
      Order = 1;
      PropertyType = PROPTYPE_KineticEnergy;
   } else if (slType == "density" || slType == "charge-density" || slType == "total-density") {
      PropertyType = PROPTYPE_ElectronDensity;
   } else if (slType == "spin-density") {
      PropertyType = PROPTYPE_SpinDensity;
   } else if (slType == "vnuc" || slType == "nuclear-potential") {
      PropertyType = PROPTYPE_NuclearAttractionPotential;
   } else if (slType == "esp" || slType == "electrostatic-potential") {
      Order = 0;
      PropertyType = PROPTYPE_Electrostatic;
   } else if (slType == "ef" || slType == "efield" || slType == "electric-field") {
      Order = 1;
      PropertyType = PROPTYPE_Electrostatic;
   } else if (slType == "efg" || slType == "electric-field-gradient") {
      Order = 2;
      PropertyType = PROPTYPE_Electrostatic;
   } else if (slType == "s" || slType == "ovl" || slType == "overlap" || slType == "zeroth-moments") {
      Order = 0;
      PropertyType = PROPTYPE_ChargeMoments_Cartesian;
   } else if (slType == "dm" || slType == "dipole-moments" || slType == "dipoles" || slType == "first-moments") {
      Order = 1;
//       PropertyType = PROPTYPE_DipoleMoments;
      // TODO: replace this back by the special dipole version once implemented.
      // (for converting stuff into debye etc.)
      PropertyType = PROPTYPE_ChargeMoments_Cartesian;
   } else if (slType == "sm" || slType == "second-moments") {
      Order = 2;
      PropertyType = PROPTYPE_ChargeMoments_Cartesian;
   } else if (slType == "tm" || slType == "third-moments") {
      Order = 3;
      PropertyType = PROPTYPE_ChargeMoments_Cartesian;
   } else {
      return false;
   }
   if (pOrderOut)
      *pOrderOut = Order;
   if (pPropertyTypeOut)
      *pPropertyTypeOut = PropertyType;
   return true;
}


void FMpt1e::SetPropertyType(string_slice slType)
{
   m_PropertyType = PROPTYPE_Unknown;
   if (!FindPropertyType(slType, &m_PropertyType, &m_Order))
      throw std::runtime_error(fmt::format("FMpt1e::FindPropertyType(): '{}' property/operator type not recognized.", slType.to_str()));
}



void FMpt1e::SetArgs(std::string const &Args_)
{
   // check if Args_ contains a { } or similar (i.e., a property list).
   // For some operators we also accept just the name, withot any properties.
   bool
      ParseAsPropertyList = false;
   for (size_t i = 0; i != Args_.size(); ++ i)
      if (Args_[i] == '{' || Args_[i] == '}')
         ParseAsPropertyList = true;
   if (!ParseAsPropertyList) {
      return SetPropertyType(string_slice(Args_.begin(), Args_.end()));
   }

   try {
      FPropertyListStr
         Props = FPropertyListStr(Args_, 0, PROPLIST_AllowEntriesWithoutName);
      if (!Props.Name.empty())
         SetPropertyType(Props.Name);
      if (m_PropertyType == PROPTYPE_Unknown)
         throw std::runtime_error(fmt::format("FMpt1e: property parameters are being set, but the operator type was never defined ('{}').", Args_));

      for (size_t iOpt = 0; iOpt != Props.size(); ++ iOpt) {
         FPropertyStr const
            &Prop = Props[iOpt];
         if (false ) {
//          if (Prop.Name.empty() && (Prop.Content == "cartesian" || Prop.Content == "cart")) {
//             AngularType = ANGTYPE_Cartesian;
//          } else if (Prop.Name.empty() && (Prop.Content == "spherical" || Prop.Content == "sph")) {
//             AngularType = ANGTYPE_RealSolidHarmonic;
//          } else if (Prop.Name.empty() && (Prop.Content == "complex")) {
//             throw FInputError("FMpt1e::SetArgs(): complex spherical harmonic multipoles are not implemented.");
         } else if (Prop.Name == "op" || Prop.Name == "type" || Prop.Name == "property" || Prop.Name == "prop") {
            SetPropertyType(Prop.Content);
         } else if (Prop.Name == "iprint" || Prop.Name == "print") {
            bool bPrint;
            if (Prop.Content.try_convert_to_bool(&bPrint)) {
               m_PrintLevel = bPrint ? 1 : 0;
            } else {
               m_PrintLevel = Prop.Content.to_int();
            }
         } else if (Prop.Name == "factor" || Prop.Name == "prefactor") {
            m_Factor = Prop.Content.to_float();
         } else if (Prop.Name == "order" || Prop.Name == "degree") {
            m_Order = Prop.Content.to_int();
         } else if (Prop.Name == "origin" || Prop.Name == "center") {
            // set a single center coordinate
            m_iCenters.clear(); m_vCenters.clear();
            m_vCenters.push_back(ParseVector3Arg(Prop.Content, 1./ToAng));
         } else if (Prop.Name == "centers" || Prop.Name == "points") {
            // set multiple center coordinates
            m_vCenters.clear();
            long_split_result srPoints;
            Prop.Content.split_list(srPoints);
            for (size_t iPt = 0; iPt != srPoints.size(); ++ iPt)
               m_vCenters.push_back(ParseVector3Arg(srPoints[iPt], 1./ToAng));
         } else if (Prop.Name == "atoms") {
            // set multiple center coordinates, by atom index
            m_iCenters.clear();
            long_split_result srPoints;
            Prop.Content.split_list(srPoints);
            for (size_t iPt = 0; iPt != srPoints.size(); ++ iPt)
               m_iCenters.push_back(srPoints[iPt].to_int() - 1);
         } else {
            // can't really report problems here, can we? It's a constructor and we have no log.
            // Make an extra function for this?
            throw std::runtime_error(fmt::format("FMpt1e: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
         }
      }
   } catch (std::runtime_error &e) {
      // convert exceptions to input errors.
      throw FInputError(e.what());
   }
}



void FMoleculePropertyTasks::SetArgs(std::string const &Args_)
{
   try {
      FPropertyListStr
         Props = FPropertyListStr(Args_, 0, PROPLIST_AllowEntriesWithoutName);
      for (size_t iOpt = 0; iOpt != Props.size(); ++ iOpt) {
         FPropertyStr
            Prop = Props[iOpt]; // note: this makes a copy, but it's just two string_slice objects (4 ptrs total).
         if (Prop.Name.empty()) {
            // assume entries without name actually denote property keys,
            // which are requested to be evaluated with default settings.
            Prop.Name = Prop.Content;
            Prop.Content = string_slice();
         }

         // start by checking if this is a known 1e property name
         if (FMpt1e::FindPropertyType(Prop.Name, 0, 0)) {
            // idea: store lists of property computation tasks, and attach their
            // results to FAtomSet objects.
            FMpt1ePtr p1eOp(new FMpt1e);
            p1eOp->SetPropertyType(Prop.Name);
            if (!Prop.Content.empty())
               p1eOp->SetArgs(Prop.Content.to_str());
            m_1eTaskList.push_back(p1eOp);
         } else {
            // can't really report problems here, can we? It's a constructor and we have no log.
            // Make an extra function for this?
            throw std::runtime_error(fmt::format("FMoleculePropertyTasks: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
         }
      }
   } catch (std::runtime_error &e) {
      // convert exceptions to input errors.
      throw FInputError(e.what());
   }
}


template<class FUnsignedIntType>
void ChangeFlagT(FUnsignedIntType &BitField, FUnsignedIntType FlagToChange, int iDir)
{
   assert(iDir == 0 || iDir == +1 || iDir == -1);
   if (iDir == +1)
      BitField |= FlagToChange;
   else if (iDir == -1)
      BitField &= ~FlagToChange;
   else if (iDir == 0) {
   } else {
      throw std::runtime_error("ChangeFlag(): iDir must be -1 (-> turn off), +1 (-> turn on), or 0 (-> leave unchanged)");
   }
}

void ChangeFlag(unsigned int &BitField, unsigned int FlagToChange, int iDir) {
   return ChangeFlagT(BitField, FlagToChange, iDir);
}

void ChangeFlag(unsigned long &BitField, unsigned long FlagToChange, int iDir) {
   return ChangeFlagT(BitField, FlagToChange, iDir);
}

// template<class FUnsignedIntType>
// void ChangeFlag(FUnsignedIntType &BitField, size_t FlagToChange, int iDir);
//
// // add explicit instanciations
// template void ChangeFlag<unsigned char>(unsigned char &, unsigned char, int);
// template void ChangeFlag<unsigned short>(unsigned short &, unsigned short, int);
// template void ChangeFlag<unsigned int>(unsigned int &, unsigned int, int);
// template void ChangeFlag<unsigned long>(unsigned long &, unsigned long, int);
//
// ^- note: un-templatized since the template would cause problems in T-parameter
//    deduction if mixing generic unsigned integer bit fields with flags from enum
//    values.



FHfGuessOptions::FHfGuessOptions()
   :  NewGuessType(INITGUESS_AtomicDensity),
      UseInputWf(USEINPUT_AlwaysIfPossible),
      PropagationType(GUESSPROP_StraightOrbitals)
{
}

void FHfGuessOptions::SetArgs(string_slice const &Args_)
{
   try {
      FPropertyListStr
         Props = FPropertyListStr(Args_, 0, PROPLIST_AllowEntriesWithoutName | PROPLIST_AllowOmitOpenClose);
      for (size_t iOpt = 0; iOpt != Props.size(); ++ iOpt) {
         FPropertyStr
            Prop = Props[iOpt];
         if (Prop.Name.empty() || Prop.Name == "type" || Prop.Name.startswith("new-guess")) {
            if (Prop.Name.empty()) {
               // if this was just specified as "scf{start: atden}" (or something to
               // this degree), assume that the user did not only mean to tell
               // use which type of guess to use for making entirely new wave
               // functions, but also that we are meant to employ this guess
               // type even if starting orbitals were provided.
               UseInputWf = USEINPUT_Never;
            }
            if (Prop.Content == "atden" || Prop.Content == "density" || Prop.Content == "atomic-density")
               NewGuessType = INITGUESS_AtomicDensity;
            else if (Prop.Content == "coreh")
               NewGuessType = INITGUESS_CoreH;
            else
               throw FInputError(fmt::format("FHfGuessOptions: failed to recognize from-scratch initial guess construction type '{}' in guess options set '{}'", Prop.Content.to_str(), Args_.to_str()));
         } else if (Prop.Name == "propagate" || Prop.Name == "prop") {
            if (Prop.Content == "orb" || Prop.Content == "orbs")
               PropagationType = GUESSPROP_StraightOrbitals;
            else if (Prop.Content == "shocc")
               PropagationType = GUESSPROP_ShOcc;
            else if (Prop.Content == "fock")
               PropagationType = GUESSPROP_Fock;
            else if (Prop.Content == "smhfock")
               PropagationType = GUESSPROP_SmhFock;
            else
               throw FInputError(fmt::format("FHfGuessOptions: initial guess propagation type '{}' not recognized.", Prop.Content.to_str()));
         } else if (Prop.Name.startswith("use-prev") || Prop.Name.startswith("use-last") ||
                    Prop.Name.startswith("use-inp") || Prop.Name.startswith("use-wf") ||
                    Prop.Name.startswith("prev") || Prop.Name.startswith("last")) {
            if (Prop.Content == "always" || Prop.Content == "auto" || Prop.Content == "yes" || Prop.Content == "on" || Prop.Content == "true") {
               UseInputWf = USEINPUT_AlwaysIfPossible;
            } else if (Prop.Content == "never" || Prop.Content == "no" || Prop.Content == "off" || Prop.Content == "false") {
               UseInputWf = USEINPUT_Never;
            } else {
               int iDir = 0; // +1 or -1 once set.
               bool bError = false;
               for (cstr_it p = Prop.Content.first; p != Prop.Content.last; ++p) {
                  if (*p == '+') {
                     if (iDir == 0) // first setting: allow nothing except what is specified.
                        UseInputWf = USEINPUT_Never;
                     iDir = +1;
                  } else if (*p == '-') {
                     if (iDir == 0) // first setting: allow everything except what is specified
                        UseInputWf = USEINPUT_AlwaysIfPossible;
                     iDir = +1;
                  } else if (iDir == 0) {
                     bError = true; break;
                  } else if (*p == 'g') {
                     ChangeFlag(UseInputWf, USEINPUT_AllowDifferentGeometry, iDir);
                  } else if (*p == 'b') {
                     ChangeFlag(UseInputWf, USEINPUT_AllowDifferentBasis, iDir);
                  } else if (*p == 'm') {
                     ChangeFlag(UseInputWf, USEINPUT_AllowDifferentMethod, iDir);
                  } else if (*p == 's') {
                     ChangeFlag(UseInputWf, USEINPUT_AllowDifferentState, iDir);
                  } else {
                     bError = true; break;
                  }
               }
               if (bError)
                  throw FInputError(fmt::format("FHfGuessOptions: option '{}: {}': input wf usage control must start with either '+' or '-': e.g., '+gbsm' allows using input data from other geometries (g), basis sets (b), electronic states (s), and scf methods (m). '-sm' forbids using input from other states or methods.", Prop.Name.to_str(), Prop.Content.to_str()));
            }
         } else {
            // can't really report problems here, can we? It's a constructor and we have no log.
            // Make an extra function for this?
            throw FInputError(fmt::format("FHfGuessOptions: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
         }
      }
   } catch (std::runtime_error &e) {
      // convert exceptions to input errors.
      throw FInputError(e.what());
   }
}



void FHfOptions::SetArgs(std::string const &CommandOptions_)
{
   // convert everything to lower case.
   std::string
      CommandOptions = tolower(CommandOptions_);
   string_slice
      slCommandOptions(CommandOptions);
   // strip whitespace if present.
   slCommandOptions.trim();
   // anything left?
   if (slCommandOptions.empty())
      // nope. nothing to do.
      return;

   // So we have some non-empty input \o/.
   // Check if there could be a property list in there ({...}).
   // If yes, split input into command part and detailed property list.
   string_slice
      slCommand,
      slPropertyList;
   {
      cstr_it
         iPropStart = slCommandOptions.find('{');
      if (iPropStart != slCommandOptions.end()) {
         // looks okay. split into command part and property list part.
         slCommand = string_slice(slCommandOptions.first, iPropStart);
         slPropertyList = string_slice(iPropStart, slCommandOptions.last);
         slCommand.trim_right();
      } else {
         // no curly braces, so no property list. Assume entire
         // input is meant to define the SCF method via moinker.
         slCommand = slCommandOptions;
         // make an empty string slice for this.
         slPropertyList = string_slice(slCommandOptions.last, slCommandOptions.last);
         slPropertyList.clear();
      }
   }

   // is there a command part?
   // (note: if we get 'scf' as "command", this is just a leftover from setting
   // a property list as scf{...}, so it should be ignored.)
   if (!slCommand.empty() && slCommand != "scf") {
      // if there *is* an explicit command, it should fully control some core
      // aspects of the method (e.g., algos for computing J/K/X and orbital/wf
      // type). For this reason, we here first reset the arguments most likely
      // to be controlled by the command (DF/JkAlgo, XcAlgo, OrbType)
      XcAlgo = XCALGO_Regular;
      JkAlgo = JKALGO_4ixDirect;
      OrbType = ORBTYPE_Restricted;

      // input could be something like "rks" or "uhf" or something more inspired
      // in which some sorts of approximations are specified first
      // (e.g., "df-", "dfj-", "dfjk-", "dfjx-", "ldf-", etc).
      // Take the command apart (at "-") and see what we get.
      split_result
         srCommandParts;
      slCommand.split(srCommandParts, '-');
      // now go through all parts and see if we can recognize everything.
      split_result::const_iterator itPart;
      for (itPart = srCommandParts.begin(); itPart != srCommandParts.end(); ++ itPart) {
         if (*itPart == "rhf" || *itPart == "hf" || *itPart == "uhf") {
            SetXc("RHF");
            if (itPart->startswith('u'))
               OrbType = ORBTYPE_Unrestricted; // we actually do not have that one at the moment.
            continue;
         } else if (*itPart == "rks" || *itPart == "ks" || *itPart == "uks") {
            if (XcFunctionalName.empty() || XcFunctionalName == "RHF")
               // if a density functional is already set, then leave it alone.
               // Otherwise, set a LDA functional... to make sure that DFT is on.
               // Will most likely be overwritten just next anyway.
               SetXc("CK16");
            if (itPart->startswith('u'))
               OrbType = ORBTYPE_Unrestricted; // we actually do not have that one at the moment.
            continue;
         } else if (itPart->startswith("df") || itPart->startswith("ri") || *itPart == "ldf" || *itPart == "ddf") {
            // setting some sort of J/K/X algorithm. Let's see what we got.
            assert(itPart->size() >= 2);
            string_slice
               // DF type without the leading "df" or "ri" (don't want to write everything twice...)
               slDfType(itPart->first+2, itPart->last);
            if (*itPart == "ldf") {
               // global density fitting Coulomb and orbital-localization based exact exchange.
               // This one affects only j and k (actually, j goes separately with this).
               // x stays as whatever it is.
               JkAlgo = JKALGO_DfCoulombDfLocalExchange;
            } else if (slDfType == "j") {
               // Coulomb fitting + regular orbital-based exchange. No K.
               JkAlgo = JKALGO_DfCoulombOnlyCache3ix;
               XcAlgo = XCALGO_Regular;
            } else if (slDfType == "jx") {
               // DF-JX: compute xc with auxiliary density from Coulomb fitting.
               // (see J. Chem. Theory Comput. 2018, 14, 1297-1303 doi.org/10.1021/acs.jctc.7b01083)
               JkAlgo = JKALGO_DfCoulombOnlyCache3ix;
               XcAlgo = XCALGO_AuxiliaryExpand;
            } else if (slDfType == "dj" || *itPart == "ddfj") {
               // integral-direct Coulomb fitting + regular orbital-based exchange. No K.
               JkAlgo = JKALGO_DfCoulombOnly;
               XcAlgo = XCALGO_Regular;
            } else if (slDfType == "cosx" || slDfType == "jcosx") {
               // regular DF-J plus Neese's real-space grid based chain-of-spheres exchange (COSX).
               JkAlgo = JKALGO_DfCoulombChainOfSpheresExchange;
               XcAlgo = XCALGO_Regular;
            } else if (slDfType == "jk") {
               // regular DF-JK, as in Hartree-Fock, and regular non-aux-expanded X/DFTI
               assert_rt(itPart->startswith("df") || itPart->startswith("ri"));
               JkAlgo = JKALGO_DfCoulombAndExchange;
               XcAlgo = XCALGO_Regular;
            } else if (slDfType == "") {
               assert_rt(*itPart == "df" || *itPart == "ri");
               // hm... this just says "use some sort of DF". kind of makes sense,
               // but it also means that we would need to decide between DF-J and DF-JK
               // based on the actually used functional (if or if not it needs exact exchange).
               // ...unfortunately, we do not actually know this at this point.
               // So, as a hack, we set DF-JK trust that FHfMethod in CtRhf.cpp will
               // sort if out if exact exchange is not actually required.
               // (might also get fixed by subsequent arguments, e.g., jxfit or jfit specs)
//                JkAlgo = JKALGO_DfCoulombAndExchange;
               JkAlgo = JKALGO_DfAuto;
               XcAlgo = XCALGO_Regular;
            } else {
               throw FInputError(fmt::format("FHfOptions: failed to recognize density fitting type '{}' in SCF command '{}'", itPart->to_str(), slCommand.to_str()));
            }
         } else {
            throw FInputError(fmt::format("FHfOptions: failed to recognize part '{}' of SCF command '{}'", itPart->to_str(), slCommand.to_str()));
         }
      }
   }

   // is there a property list part?
   if (!slPropertyList.empty()) {
      try {
         FPropertyListStr
            MethodProps = FPropertyListStr(slCommandOptions, 0, PROPLIST_AllowEntriesWithoutName);
         for (size_t iOpt = 0; iOpt != MethodProps.size(); ++ iOpt) {
            FPropertyStr const
               &Prop = MethodProps[iOpt];
            if (Prop.Name.empty()) {
               if (Prop.Content == "ignore-error") {
                  IgnoreConvergenceError = true;
               }
            } else if (Prop.Name == "thrden" || Prop.Name == "thr-den") {
               ThrDen = Prop.Content.to_float();
            } else if (Prop.Name == "throrb" || Prop.Name == "thr-orb") {
               ThrOrb = Prop.Content.to_float();
            } else if (Prop.Name == "maxit" || Prop.Name == "max-it" || Prop.Name == "max-iter") {
               MaxIt = unsigned(Prop.Content.to_int());
            } else if (Prop.Name == "maxdiis" || Prop.Name == "max-diis") {
               MaxDiis = unsigned(Prop.Content.to_int());
            } else if (Prop.Name == "dfxc" || Prop.Name == "dfjx" || Prop.Name == "df-jx") {
               SetDfxc(Prop.Content.to_bool());
            } else if (Prop.Name == "4ix") {
               if (Prop.Content.to_bool())
                  JkAlgo = JKALGO_4ixDirect;
            } else if (Prop.Name == "jkalgo" || Prop.Name == "jk") {
               if (Prop.Content == "dfjk")
                  JkAlgo = JKALGO_DfCoulombAndExchange;
               else if (Prop.Content == "3ix-direct" || Prop.Content == "df-j")
                  JkAlgo = JKALGO_DfCoulombOnly;
               else if (Prop.Content == "3ix-cached")
                  JkAlgo = JKALGO_DfCoulombOnlyCache3ix;
               else if (Prop.Content == "ldf-jk")
                  JkAlgo = JKALGO_DfCoulombDfLocalExchange;
               else if (Prop.Content == "cosx" || Prop.Content == "dfj-cosx")
                  JkAlgo = JKALGO_DfCoulombChainOfSpheresExchange;
               else if (Prop.Content == "4ix-direct")
                  JkAlgo = JKALGO_4ixDirect;
               else
                  throw std::runtime_error(fmt::format("FHfOptions: explicitly set jkalgo '{}' not recognized.", Prop.Content.to_str()));
            } else if (Prop.Name == "xc") {
               SetXc(Prop.Content.to_str());
            } else if (Prop.Name == "shift") {
               LevelShifts[0] = Prop.Content.to_float();
               LevelShifts[1] = 0.5 * LevelShifts[0];
            } else if (Prop.Name == "damping" || Prop.Name == "damp") {
               SetDamping(Prop.Content);
            } else if (Prop.Name == "guess" || Prop.Name == "start" || Prop.Name == "initial-guess") {
               InitialGuess.SetArgs(Prop.Content);
//             } else if (Prop.Name == "basis-orb") {
//                BasisName_Orbital = Prop.Content.to_str();
//             } else if (Prop.Name == "basis-fit") {
//                BasisName_Fit = Prop.Content.to_str();
// ^- not really controlled from within here.
            } else if (Prop.Name == "grid") {
               SetGridDesc(Prop.Content.to_str());
            } else if (Prop.Name == "disp") {
               DispCorrType = Prop.Content.to_str();
            } else if (Prop.Name == "timing") {
               iTimingLevel = Prop.Content.to_int();
            } else if (Prop.Name == "print") {
               Print.SetArgs(Prop.Content.to_str());
            } else if (Prop.Name == "eval" || Prop.Name == "properties") {
               PropertyTasks.SetArgs(Prop.Content.to_str());
            } else if (Prop.Name == "band-energy-formula") {
               // note: this is an experimental option; it is used for testing
               // some formal rewritings of the Kohn-Sham energy which may be
               // more suitable for tight-binding and machine-learning
               // approximations than the regular KS formula is.
               //
               // Note: not fully implemented in all variants. Closed-shell LDA
               // should work, though (at least it did at some point)
               UseBandStructureEnergyFormula = Prop.Content.to_bool();
            } else if (Prop.Name == "rohf-algo" || Prop.Name == "algo") {
               if (Prop.Content == "1" || Prop.Content == "one-step") {
                  m_RohfAlgo =  ROHFALGO_OneStep;
               } else if (Prop.Content == "2" || Prop.Content == "two-step") {
                  m_RohfAlgo =  ROHFALGO_TwoStep;
               } else {
                  throw std::runtime_error(fmt::format("FHfOptions: '{}' must be either 'one-step' or 'two-step', but got '{}'.", Prop.Name.to_str(), Prop.Content.to_str()));
               }
            } else {
               // can't really report problems here, can we? It's a constructor and we have no log.
               // Make an extra function for this?
               throw std::runtime_error(fmt::format("FHfOptions: scf option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
            }
         }
      } catch (std::runtime_error &e) {
         // convert exceptions to input errors.
         throw FInputError(e.what());
      }
   }
}




} // namespace ct
