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

// #include <iostream>
#include <cmath>
#include <algorithm> // for std::sort
#include <sstream>
#include <stdexcept> // for std::runtime_error
#include "CxDiis.h"
#include "CtRhf.h"
#include "CtFockBuild.h"
#include "CtRhfGuess.h"
// #include "CtRhfProp.h"
#include "CxIo.h"
#include "CxTiming.h"
#include "CxPodArray.h"
#if 0
   #include "RDTSC.h"
#else
   #define RESET_CLOCKS
   #define RESUME_CLOCK(x)
   #define PAUSE_CLOCK(x)
#endif

#include "CxPhysicalUnits.h"
#include "CtBasisLibrary.h" // for ECP print only

namespace ct {


FHfResult::FHfResult()
   : pAtoms(0), nIt(0), Converged(false)
{
   RohfFockIsLevelShifted = false;
}

FHfResult::~FHfResult()
{
}



FHfMethod::FHfMethod(FLog &Log_, FWfDecl const &WfDecl_, FAtomSetCptr const pAtoms_, FHfOptions const &Options_, FHfResult *pStartingGuess_, unsigned ExecFlags, FMemoryStack &Mem)
   : m_Options(Options_), m_WfDecl(WfDecl_), m_pLog(&Log_), m_ExecFlags(ExecFlags), m_pAtoms(pAtoms_.get()), m_HaveOrbitals(false)
{
   m_WfDecl.SetNuclearCharge(m_pAtoms->NuclearCharge()); // need this to compute the total number of electrons.
   if (((size_t(m_WfDecl.Ms2) % 2) != (size_t(m_WfDecl.nElec()) % 2)) || (m_WfDecl.Ms2 > m_WfDecl.nElec()))
      throw std::runtime_error(fmt::format("Ms2 {} inconsistent with total number of electrons {}.", m_WfDecl.Ms2, m_WfDecl.nElec()));
   if (!bool(m_ExecFlags & EXEC_DiscardResults))
      m_pResult = FHfResultPtr(new FHfResult);

   m_Energy = 0.;
   m_tTotal.Reset();
   m_RohfAlgo = m_Options.m_RohfAlgo;
   if (m_WfDecl.nOpen() == 0) {
      // no open shells. no point in making A/B Fock matrices.
      m_RohfAlgo = FHfOptions::ROHFALGO_OneStep;
   }


   void
      *pBaseOfMemory = Mem.Alloc(0);
   Init(pStartingGuess_, Mem);
   m_pLog->WriteLine();
   m_pLog->WriteTiming("initialization", (double)m_tTotal);
   Run(Mem);
//    FDmaDescPtr
//       pDma;
//    EvalDma(pDma, Mem);
   if (m_pResult) {
      // store main results of computation if requested.
      m_pResult->Energy = m_Energy;
      m_pResult->pAtoms = pAtoms_;
      // ^- WARNING: make sure you use the outside-provided intrusive ptr
      // for this. This is *not* the same as this->pAtoms! (raw pointer, atm).
      m_pResult->pOrbBasis = m_pOrbBasis;
      Assign(m_pResult->Orb, m_COrb);
      Assign(m_pResult->Fock, m_Fock); // <- not quite the same as S C diag(ew) C.T S.T ... since EW come from the LAST iteration.
      Assign(m_pResult->ExchO, m_ExchO);

      Assign(m_pResult->Ew, m_Ew);
      size_t nAo = _GetNumBf();
      m_pResult->Occ.Reshape(nAo,1);
      m_pResult->Occ.Clear();
      for (size_t i = 0; i < m_WfDecl.nElecB(); ++ i)
         m_pResult->Occ[i] = 2.;
      for (size_t i = m_WfDecl.nElecB(); i < m_WfDecl.nElecA(); ++ i)
         m_pResult->Occ[i] = 1.;

      if (1) {
         // delete previous wf propagation info, if present.
         // (only relevant if guess propagation type is changed between calls)
         m_pResult->ShOcc.Delete();
         m_pResult->RohfFock.Delete();
         m_pResult->SmhRohfFock.Delete();
         m_pResult->SmhExchO.Delete();
         m_pResult->RohfFockIsLevelShifted = false;
         m_pResult->ShTrafoBasisPerm.clear();
      }

      // make a copy of the settings which were used to obtain the results stored here.
      m_pResult->HfOptions = m_Options; // <-- this is actually not a small object anymore. Really a good idea?
      m_pResult->WfDecl = m_WfDecl;
   }

   if (!bool(ExecFlags & EXEC_RetainMemory))
      Mem.Free(pBaseOfMemory);
}


void SetupFockBuilderList(FFockComponentBuilderList &FockComponentBuilders, FHfOptions const &Options, FLog &Log, FTimerSet *pTimers)
{
   if (!FockComponentBuilders.empty())
      // some of the JK builders do not actually *ACCUMULATE* contributions in
      // AccFock, but simply overwrite the target matrix...  this is technically
      // a bug, but in practice I guess we should better just ignore it for the
      // time being. And simply insist that the ones which do that come in first
      // in the evaluation order. It's a bit messy.
      throw std::runtime_error("SetupFockBuilderList: JkAlgos need to come first.");

   // start with builders for the base components: J, X, and K.
   FXcOptions
      XcOptions(Options.XcFunctionalName, &Options.DftGridParams);
   XcOptions.UseBandStructureEnergyFormula = Options.UseBandStructureEnergyFormula;
   XcOptions.UseDfxc = (Options.XcAlgo == XCALGO_AuxiliaryExpand); // <- remove this?

   // check what kind of algorithm to use for J and K, in case we were
   // given a choice.
   FJkAlgo
      JkAlgo = Options.JkAlgo;
   if (JkAlgo == JKALGO_DfAuto) {
      if (!XcOptions.pXcFn || XcOptions.pXcFn->NeedsExactExch()) {
         // no functional (->Hartree Fock) or hybrid functional. Use something
         // which can make K matrices.
         JkAlgo = JKALGO_DfCoulombAndExchange;
      } else {
         // no exact exchange needed. Use the DF-JX builder with cached
         // integrals by default. Maybe some day I even get around to
         // do selective caching in low memory situations, like I implemented
         // in Molpro ten years ago...
         JkAlgo = JKALGO_DfCoulombOnlyCache3ix;
      }
   }

   FDfJkOptions
      DfJkOptions;
   DfJkOptions.UseBandStructureEnergyFormula = Options.UseBandStructureEnergyFormula;
   if (XcOptions.pXcFn) {
      DfJkOptions.fExchFactor = XcOptions.pXcFn->ExchFactor();
      if (XcOptions.pXcFn->pExchKernel() != 0)
         throw std::runtime_error("SetupFockBuilderList: setup for xc functionals with non-standard exact exchange kernels not fully implemented.");
   }

   bool
      XcAlreadyAdded = false,
      HaveExactExchange = false;

//    if (JkAlgo == JKALGO_DfCoulombAndExchange) {
//       bool brrrk = true;
// #ifdef INCLUDE_OPTIONAL_FOCK_BUILDERS
//       FockComponentBuilders.push_back(new FFockComponentBuilderDfJk(DfJkOptions, Log, pTimers));
//       brrrk = false;
// #endif // INCLUDE_OPTIONAL_FOCK_BUILDERS
//       if (brrrk) throw std::runtime_error("This version was compiled without DF-RHF support.");
//       assert_rt(Options.XcFunctionalName.empty()); // regular XC not implemented (hm... wait, why?!).
//    }


   if (JkAlgo == JKALGO_DfCoulombAndExchange) {
      FockComponentBuilders.push_back(new FFockComponentBuilderDfJk(DfJkOptions, Log, pTimers));
      HaveExactExchange = true;
   } else if (JkAlgo == JKALGO_DfCoulombOnlyCache3ix) {
      // auxiliary coulomb and exchange-correlation with expanded densities and cached integrals.
      FockComponentBuilders.push_back(new FFockComponentBuilderDfCoulXcCached(
         DfJkOptions, XcOptions, Log, pTimers));
      XcAlreadyAdded = true;
      HaveExactExchange = false;
   } else if (JkAlgo == JKALGO_4ixDirect) {
      FockComponentBuilders.push_back(new FFockComponentBuilder4ixJk(DfJkOptions, Log, pTimers));
      HaveExactExchange = true;
   } else if (JkAlgo == JKALGO_DfCoulombOnly) {
      FockComponentBuilders.push_back(new FFockComponentBuilderDfCoul(DfJkOptions, Log, pTimers));
      HaveExactExchange = false;
   } else {
      throw std::runtime_error(fmt::format("SetupFockBuilderList: JkAlgo {} not recognized.", unsigned(JkAlgo)));
//       FockComponentBuilders.push_back(new FFockComponentBuilderDfCoul(DfJkOptions));
//       FockComponentBuilders.push_back(new FFockComponentBuilderDfLx(DfJkOptions, Overlap, Scd));
   }

   if (XcOptions.pXcFn && XcOptions.pXcFn->NeedsExactExch() && !HaveExactExchange)
      throw std::runtime_error(fmt::format("SetupFockBuilderList: xc Functional '{}' requires exact exchange, but JkAlgo was set to {} which cannot provide K matrices.", XcOptions.pXcFn->Name(), unsigned(JkAlgo)));

   if (XcOptions.pXcFn && !XcAlreadyAdded) {
      // some sort of xc functional is declared. Add a regular orbital-based XC integrator
      // unless one of the Fock builders above has already done so.
      if (Options.XcAlgo == XCALGO_Regular) {
//          FockComponentBuilders.push_back(new FFockComponentBuilderXc(XcOptions, Log, pTimers));
         FockComponentBuilders.insert(FockComponentBuilders.begin(), new FFockComponentBuilderXc(XcOptions, Log, pTimers));
      } else {
         throw std::runtime_error(fmt::format("SetupFockBuilderList: XcAlgo {} not supported for current JkAlgo ({}).", unsigned(Options.XcAlgo), unsigned(JkAlgo)));
      }
   }


   if (!(Options.DispCorrType.empty() || Options.DispCorrType == "none" || Options.DispCorrType == "off")) {
      FockComponentBuilders.push_back(new FDispCorrBuilder(Options.DispCorrType, Options.XcFunctionalName, Log, pTimers));
   }
   if (!(Options.ImplicitSolvationModelDesc.empty() || Options.ImplicitSolvationModelDesc == "none" || Options.ImplicitSolvationModelDesc == "off")) {
   }
}


void PrintIndented(FLog *pLog, std::string const &sText, std::string const &Ind)
{
   if (!pLog) return;
   string_slice
      slText(sText);
   long_split_result
      srLines;
   slText.split_lines(srLines);
   for (string_slice const &slLine : srLines) {
      pLog->w << Ind << slLine.to_str() << "\n";
   }
}


void FHfMethod::Init(FHfResult *pStartingGuess, FMemoryStack &Mem)
{
   m_pLog->Write("\n" " :::: MicroScf [v20211019-RevA] -- Integrated DFT Driver                     ::::");
   m_pLog->Write(     " :::: \"...no such thing as overkill.\"      Developed by: Gerald Knizia, 2013 ::::");

   if (m_Options.UseXc())
      m_pLog->WriteProgramIntro("RESTRICTED KOHN-SHAM");
   else
      m_pLog->WriteProgramIntro("RESTRICTED HARTREE-FOCK");

//    m_pLog->Write("!{}", GetMethodWfDesc(&Options, &WfDecl));
   if (m_Options.Print[SCFPRINT_Default])
      PrintIndented(m_pLog, GetMethodWfDesc(&m_Options, &m_WfDecl), " >>> ");

   if (m_Options.OrbType != ORBTYPE_Restricted) {
      throw std::runtime_error("FHfMethod only supports restricted orbitals (RKS/RHF), not unrestricted or generalized.");
   }

   m_pOrbBasis = new FBasisSet(*m_pAtoms, BASIS_Orbital);
//    pFitBasis = new FBasisSet(*pAtoms, BASIS_JkFit);
   m_pMinBasis = new FBasisSet(*m_pAtoms, BASIS_MinAo);
   // ^-- hmm... this is currently *only* used for determining the number of
   //     valence orbitals in the valence-level orbital energy print. Even
   //     when making a density guess, it just re-instanciates the basis again.
   //     (and, technically, in the meantime the sets assigned for the MINAO and
   //     GUESS contexts could in principle be chosen differently)
   m_nOrbBasisSize = m_pOrbBasis->nFn();
   size_t nAo = _GetNumBf();
//    nFit = pFitBasis->nFn();
//    m_pLog->Write(" {:<31}{:5} (ORB) //{:5} ({})", "Number of basis functions:", nAo, nFit, GetFitBasisType(&Options, true));
   m_pLog->Write(" {:<31}{:5} (ORB) //{:5} ({})", "Number of basis functions:", nAo, m_pMinBasis->nFn(), "MIN");
   m_pLog->Write(" {:<31}{:5} (TOT) //{:5} (ALPHA) //{:5} (BETA)", "Number of electrons:", m_WfDecl.nElec(), m_WfDecl.nElecA(), m_WfDecl.nElecB());

   // print ECP data in case it was explicitly requested
   if (m_Options.Print[SCFPRINT_Ecp]) {
      std::stringstream str;
      for (size_t iAt = 0; iAt < m_pAtoms->size(); ++ iAt) {
         FAtom const
            &At = (*m_pAtoms)[iAt];
         if (At.nEcpElec != 0) {
#ifdef IR_ECP
            ir::FAtomEcp const
               *pIrEcp = g_BasisSetLibrary.LoadEcp(At.iElement, At.EcpDesc, At.vPos, int(iAt));
            str << "\n";
            PrintEcpData(str, pIrEcp, At.iElement, At.vPos, " | ");
#else
            throw std::runtime_error(fmt::format("Atom {} has declared a ECP{}, but this version of IR was compiled without ECP support.", (*m_pAtoms)[iAt].GetAtomLabel(iAt), At.nEcpElec));
#endif // IR_ECP
         }
      }
      m_pLog->w << str.str() << "\n";
   }

   // print orbital basis in case it was explicitly requested
   if (m_Options.Print[SCFPRINT_OrbBasis]) {
      m_pLog->w << "\n" << m_pOrbBasis->GetDesc(FBasisSet::PRINT_Default, " | ") << "\n";
   }

   m_pLog->WriteLine();

   FTimer Timer1e;
   m_CoreH = MakeStackMatrix(nAo, nAo, Mem);
   m_pAtoms->MakeCoreHamiltonMatrix(m_CoreH, *m_pOrbBasis, *m_pOrbBasis, Mem);

   m_Ew = MakeStackMatrix(nAo, 1, Mem);
   m_Smh = FMatrixView(0,0,0);
   m_Scd = FMatrixView(0,0,0);
   if (false)
      // not needed anymore.
      m_Smh = MakeStackMatrix(nAo, nAo, Mem);
   if (true)
      m_Scd = MakeStackMatrix(nAo, nAo, Mem);

   {
      FStackMatrix
         Overlap(nAo, nAo, &Mem);
      m_pAtoms->MakeOverlapMatrix(Overlap, *m_pOrbBasis, *m_pOrbBasis, Mem);
//       Overlap.Print(xout, "AO-Basis Overlap Matrix");

      if (m_Smh.pData) {
         Assign(m_Smh, Overlap);
         {
            std::stringstream str;
            CalcSmhMatrix(m_Smh, Mem, FSmhOptions(1e-15,1e-15," S^{-1/2} of orbital basis", &str));
            m_pLog->WriteNoNl(str.str());
         }
      }
      if (m_Scd.pData) {
         Assign(m_Scd, Overlap);
         CalcCholeskyFactors(m_Scd);
      }
   }
   m_pLog->WriteTiming("1e integrals", (double)Timer1e);


   if (m_RohfAlgo == FHfOptions::ROHFALGO_TwoStep) {
      m_RohfFock = MakeStackMatrix(nAo, nAo, Mem);
      m_FockA = MakeStackMatrix(nAo, nAo, Mem);
      m_FockB = MakeStackMatrix(nAo, nAo, Mem);
      m_FockSize = m_FockB.pData + m_FockB.GetStridedSize() - m_RohfFock.pData;
   } else {
      m_RohfFock = MakeStackMatrix(nAo, nAo, Mem);
      m_FockSize = m_RohfFock.GetStridedSize();
      m_FockA = FMatrixView(0,0,0);
      m_FockB = FMatrixView(0,0,0);
   }
   if (m_WfDecl.nOpen() != 0) {
      m_Fock = MakeStackMatrix(nAo, nAo, Mem);
   } else {
      // alias Fock to RohfFock in closed-shell case.
      m_Fock = m_RohfFock;
   }
   m_ExchO = MakeStackMatrix(nAo, nAo, Mem); // <- WARNING: this is not extrapolated! Check where used?
//    if (WfDecl.nOpen() != 0)
//       FockSpinFree = MakeStackMatrix(nAo, nAo, Mem);
   m_COrb = MakeStackMatrix(nAo, nAo, Mem);

   SetupFockBuilderList(m_FockComponentBuilders, m_Options, *m_pLog, &m_Timers);

   // note: this may allocate memory.
   // The point of separating construction/init is to allow, in principle, to setup a
   // builder chain without referencing the concrete system first.
   for (size_t i = 0; i < m_FockComponentBuilders.size(); ++ i) {
//       m_pLog->Write("...init: {}", typeid(m_FockComponentBuilders[i]).name());
      m_FockComponentBuilders[i]->Init(m_WfDecl, &*m_pOrbBasis, *m_pAtoms, m_Options, Mem);
   }

   BuildInitialFock(pStartingGuess, Mem);
}


bool FHfMethod::CanUseStartingGuess(FHfResult *pStartingGuess)
{
   if (pStartingGuess == 0)
      // nothing provided. Certainly can't use that.
      return false;
   if (m_Options.InitialGuess.UseInputWf == FHfGuessOptions::USEINPUT_Never)
      // whatever it is which was provided: we are not allowed to use it.
      return false;

   // check this actually is the same type of molecule as we have now
   // (same elements in same order). If not, we cannot use this.
   if (!m_pAtoms->IsSameMolecule(*pStartingGuess->pAtoms))
      return false;

   if (!bool(m_Options.InitialGuess.UseInputWf & FHfGuessOptions::USEINPUT_AllowDifferentState)) {
      if (!m_WfDecl.HasEqualChargeSpinSym(pStartingGuess->WfDecl))
         // we were told to not use wave functions from a different state, but the
         // previous state has a different WfDecl.
         return false;
   }
   if (!bool(m_Options.InitialGuess.UseInputWf & FHfGuessOptions::USEINPUT_AllowDifferentGeometry)) {
      if (!m_pAtoms->IsSameGeometry(*pStartingGuess->pAtoms))
         return false;
      // ^- WARNING: currently some parts of the code pass FAtomSet objects
      //    to FHfMethod, which end up in FHfResult objects, and then, afterwards,
      //    *change* the original atom sets. In this case the atom set and
      //    the rest of the FHfResult object go out of sync, and this test
      //    here will fail to perform the correct test.
      //
      //    Not sure how to best deal with this, atm.
   }
   if (!bool(m_Options.InitialGuess.UseInputWf & FHfGuessOptions::USEINPUT_AllowDifferentBasis)) {
      // note: this one checks if the same library objects for the element bases
      // are linked for all centers, and in the same order. It does not compare
      // the positions of the basis functions (i.e., it checks whether the basis
      // consists of the same basis functions, but not if they are located at
      // the same place).
      if (!m_pOrbBasis->HasSameElementBases(*pStartingGuess->pOrbBasis))
         return false;
   }

   // TODO: USEINPUT_AllowDifferentMethod is ignored.
   return true;
}




void FHfMethod::BuildInitialFock(FHfResult *pStartingGuess, FMemoryStack &Mem)
{
   // This routine prepares the RHF state for the first iteration:
   // - Our RHF loop in Run() starts with constructing a new orbital set from
   //   the current Fock matrix set. So before entering the first iteration, the
   //   Fock matrices (*not* the orbitals!) must be set up.
   //
   // - This routine constructs this initial guess. It's goal is to assemble a
   //   complete set of Fock matrices.
   //   Depending on the case, this may be:
   //
   //   + just this->Fock alone (closed-shell RHF/RKS, or ROHF/ROKS which is
   //     started from scratch, without pStartingGuess provided)
   //   + this->Fock & this->RohfFock (ROHF/ROKS with 1-step algo, starting
   //     from a previous SCF result)
   //   + Up to this->Fock, this->RohfFock, this->FockA and/or this->FockB
   //     (ROHF/ROKS with starting orbitals and two-step ROHF algo)

//    m_Timers.SetLevel(8);
   FTimer
      tGuess;
   std::string
      // short desc recording the type of initial guess constructed.
      // This being non-empty also doubles as our marker for
      // "an initial guess (already) has been set up."
      sGuessNote;
   size_t
      nAo = _GetNumBf();
   m_FockMatrixIsLevelShifted = false;

   // did we get starting information from the outside? (e.g., from previous
   // geometries in a geometry optimization). And if yes, may be use it?
   if (CanUseStartingGuess(pStartingGuess)) {
      // build guess from propagated orbitals or Sh-transformed orbitals?
      if ( (m_Options.InitialGuess.PropagationType == FHfGuessOptions::GUESSPROP_ShOcc ||
            m_Options.InitialGuess.PropagationType == FHfGuessOptions::GUESSPROP_StraightOrbitals) &&
            pStartingGuess->Orb.IsAssignedAndHasShape(nAo,nAo) ) {
         // got a compatible square matrix as input.. likely from a previous geometry step.
         // just copy and orthogonalize it, then make a Fock matrix. Hopefully it already is
         // ordered in such a way that the core orbitals come first.
         bool UseShOcc = (m_Options.InitialGuess.PropagationType == FHfGuessOptions::GUESSPROP_ShOcc);
         if (pStartingGuess->ShOcc.pData == 0 || pStartingGuess->ShOcc.nRows != nAo || pStartingGuess->ShOcc.nCols > nAo)
            // fall back to standard orbital propagation if we were told to use
            // ShOcc, but no ShOcc data was provided.
            UseShOcc = false;
         if (!UseShOcc) {
            Assign(m_COrb, pStartingGuess->Orb);
            m_Energy = pStartingGuess->Energy;
            sGuessNote = "input orbitals";
         } else {
         }

         OrthSchmidtScd(m_COrb, m_Scd, Mem);
         if (1)
            OrthSchmidtScd(m_COrb, m_Scd, Mem); // yes, this does it twice, and that is intentional.
         m_HaveOrbitals = true;
//          xout << fmt::format("!! invoke BuildRhfFock via BuildInitialFock. m_HaveOrbitals = {}  WfDecl.nOpn() == {}", m_HaveOrbitals, WfDecl.nOpen()) << std::endl;
//          _DbgCheckOrbitalOverlap(m_COrb, m_pOrbBasis.get(), "before RhfFock", Mem);
         BuildRhfFock(m_COrb, false, Mem);
         ApplyLevelShifts1(m_Options.InitialDampingLevelShift, Mem);
         // ^-- WARNING: without this some open-shell systems do not converge
         //     even with exact input orbitals, due to orbital order swapping.

//          Scale(COccC, std::sqrt(1./2.));
//          GuessDone = true;
      } else if ( m_Options.InitialGuess.PropagationType == FHfGuessOptions::GUESSPROP_Fock &&
                  pStartingGuess->RohfFock.IsAssignedAndHasShape(nAo,nAo) )
      {
         // build guess from propagated Fock matrices.
         Assign(m_RohfFock, pStartingGuess->RohfFock);
         if (m_ExchO.pData && pStartingGuess->ExchO.pData)
            Assign(m_ExchO, pStartingGuess->ExchO);
         m_FockMatrixIsLevelShifted = pStartingGuess->RohfFockIsLevelShifted;
         m_Energy = pStartingGuess->Energy;
         sGuessNote = "input Fock";
//          GuessDone = true;
         m_HaveOrbitals = false;
      } else if ( m_Options.InitialGuess.PropagationType == FHfGuessOptions::GUESSPROP_SmhFock &&
                  pStartingGuess->SmhRohfFock.IsAssignedAndHasShape(nAo,nAo) )
      {
      }
   }

   if (sGuessNote.empty()) {
      // either no starting info was provided, or the provided info could
      // not be convertd to a working starting guess (e.g., because basis
      // dimensions did not work out). Make a new guess from scratch.
      FBasisSetPtr
         pGuessBasis = new FBasisSet(*m_pAtoms, BASIS_Guess);
      if (m_Options.InitialGuess.NewGuessType == FHfGuessOptions::INITGUESS_AtomicDensity) {
         MakeDensityGuess(m_Fock, *m_pAtoms, m_Scd, m_CoreH, &*m_pOrbBasis, &*pGuessBasis, m_FockComponentBuilders, m_pLog, m_Options, Mem);
      } else if (m_Options.InitialGuess.NewGuessType == FHfGuessOptions::INITGUESS_CoreH) {
         Assign(m_Fock, m_CoreH);
      } else {
         throw std::runtime_error(fmt::format("FHfMethod::BuildInitialFock(): SCF initial guess type '{}' not recognized", int(m_Options.InitialGuess.NewGuessType)));
      }
      Assign(m_RohfFock, m_Fock);
      if (m_RohfAlgo == FHfOptions::ROHFALGO_TwoStep) {
         Assign(m_FockA, m_Fock);
         Assign(m_FockB, m_Fock);
      }
      // ^- can't use BuildRhfFock here, because we do not have a set of
      // sensible orbitals at this point---only the Fock matrices.
      m_FockMatrixIsLevelShifted = false;
      sGuessNote = "atomic density guess";
//       GuessDone = true;
      m_Energy = 0.;
      m_HaveOrbitals = false;
   }

//    if (!m_HaveOrbitals && WfDecl.nOpen() != 0) {
   if (false && !m_HaveOrbitals) {
      // make sure that some set of orbitals is available, even in the first
      // iteration. May be required if Fock matrix is diagonalized in previous
      // MO basis (e.g., for occupation locking or some open-shell orbital
      // construction algorithms), and it is assumed that those orbitals
      // are (a) complete and (b) orthonormal.
      //
      // UPDATE: current code does not require valid set of initial orbitals,
      // so this is disabled here. m_HaveOrbitals is handled in
      // UpdateOrbitals() code.
      Assign(m_COrb, m_Fock);
      TriangularSolve(m_COrb, m_Scd, 'R', 'T');
      TriangularSolve(m_COrb, m_Scd, 'L', 'N');
      Diagonalize(m_COrb, m_Ew.pData, Mem);
      TriangularSolve(m_COrb, m_Scd, 'L', 'T');
      m_HaveOrbitals = true;
   }

//    assert(GuessDone);
   assert(!sGuessNote.empty());
   m_pLog->WriteTiming("initial guess", (double)tGuess, sGuessNote);
   if (m_Options.iTimingLevel > 0) {
      std::stringstream str;
      m_Timers.PrintReport(str);
      m_pLog->WriteNoNl(str.str());
      m_Timers.Reset();
   }
}


void FHfMethod::ApplyLevelShifts1(double ExtraShift, FMemoryStack &Mem)
{
   TMemoryLock<double> pFreeMe1(0, &Mem);
   // get occupied orbitals without absorbed occupation numbers (allocated on Mem)
   FMatrixView COccC, COccO;
   // ^- I think technically we do not need all columns---only the first nClosed+nOpen ones.
   //    could possibly make a difference in initial guess generation/propagation.
   ExtractOrbitals(COccC, COccO, m_COrb, ORBEXTRACT_AllowAliasing, Mem);
//    m_pLog->Write("! Invoke: ApplyLevelShiftsToInitialGuess({})  nOccC = {}  nOccO = {}", f0, COccC.nCols, COccO.nCols);
   ApplyLevelShifts(1., COccC, 1., COccO, 1., ExtraShift, 0.5 * ExtraShift, Mem);
   // ^- Notes:
   //    - The parameters represented by ExtraShift and 0.5 * ... are
   //      the *extra* level shift in addition to the defaults given in Option.LevelShifts.
   //    - If COccC and COccO have absorbed occupation numbers, they need to be provided
   //      to ApplyLevelShifts. ApplyLevelShift1 then correctly accounts for them.
}




void BuildFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, unsigned Flags, FFockComponentBuilderList &FockBuilders, FMemoryStack &Mem, FTimerSet *pTimerSet)
{
//    TIME_SECTION(pTimerSet, 0x1ff, "Fock matrix");
//    if (pTimerSet) {
//       pTimerSet->Enter(0x1ff, "Fock matrix");
//    }
   FockC.Clear();
   if (FockO.pData)
      FockO.Clear();
   for (size_t iBuilder = 0; iBuilder < FockBuilders.size(); ++ iBuilder)
      FockBuilders[iBuilder]->AccFock(FockC, FockO, pBasis, COccC, COccO, Flags, Mem);
//    if (pTimerSet) {
//       pTimerSet->Leave(0x1ff);
//    }
   IR_SUPPRESS_UNUSED_WARNING(pTimerSet);
}


void FHfMethod::ApplyLevelShiftToSide(FMatrixView Fock0, FMatrixView Fock1, FMatrixView Fock2, FMatrixView COcc_, double fOccNumber, double LevelShift01, double LevelShift2, double fDirection, FMemoryStack &Mem)
{
   size_t
      nAo = _GetNumBf();
   if ((LevelShift01 != 0. || LevelShift2 != 0.) && COcc_.nCols != 0) {
      // add closed-shell level shift
      FStackMatrix
         T1(nAo, COcc_.nCols, &Mem);
      Assign(T1, COcc_, 1./std::sqrt(fOccNumber)); // remove absorbed occupation number
      TriangularMxm(T1, m_Scd, 'L', 'T');
      TriangularMxm(T1, m_Scd, 'L', 'N');
      if (Fock1.pData == 0) {
         SyrkNT(Fock0, T1, fDirection * LevelShift01, MXM_Add);
      } else {
         FStackMatrix
            Tmp(nAo, nAo, &Mem);
         SyrkNT(Tmp, T1, fDirection);
         Add(Fock0, Tmp, LevelShift01);
         Add(Fock1, Tmp, LevelShift01);
         if (Fock2.pData)
            Add(Fock2, Tmp, LevelShift2);
      }
      m_FockMatrixIsLevelShifted = true;
   }
}

void FHfMethod::ApplyLevelShifts(double fDirection, FMatrixView COccC, double fAbsorbedOccC, FMatrixView COccO, double fAbsorbedOccO, double fExtraShiftC, double fExtraShiftO, FMemoryStack &Mem)
{
   ApplyLevelShiftToSide(m_RohfFock, m_FockA, FMatrixView(), COccC, fAbsorbedOccC, m_Options.LevelShifts[0] + fExtraShiftC, 0., fDirection, Mem);
   ApplyLevelShiftToSide(m_RohfFock, m_FockB, m_FockA, COccO, fAbsorbedOccO, fExtraShiftO + m_Options.LevelShifts[1], fExtraShiftC + m_Options.LevelShifts[0], fDirection, Mem);
}


// /*FMatrixView COccC, FMatrixView COccO,*/
void FHfMethod::BuildRhfFock(FMatrixView const COrb_, bool InRefinementIt, FMemoryStack &Mem)
{
   TMemoryLock<double> pFreeMe0(0, &Mem);
   size_t nAo = _GetNumBf();
   {
      TMemoryLock<double> pFreeMe1(0, &Mem);
      // get occupied orbitals with absorbed occupation numbers (allocated on Mem)
      FMatrixView COccC, COccO;
      assert_rt(COrb_.IsAssignedAndHasShape(nAo, nAo));
      // ^- I think technically we do not need all columns---only the first nClosed+nOpen ones.
      //    could possibly make a difference in initial guess generation/propagation.
      ExtractOrbitals(COccC, COccO, COrb_, ORBEXTRACT_AbsorbOcc, Mem);

      assert_rt(COccC.pData != COrb_.pData);
      // ^-- below COrb_ is used for the ROHF projections. Cannot have absorbed occupation numbers in that.
      BuildFock(m_Fock, m_ExchO, &*m_pOrbBasis, COccC, COccO, 0, m_FockComponentBuilders, Mem, &m_Timers);
      Add(m_Fock, m_CoreH);
      m_FockMatrixIsLevelShifted = false;
   }

   if (m_WfDecl.nOpen() != 0 && m_HaveOrbitals) {
      FStackMatrix
         SOrb(nAo, nAo, &Mem); // = S * C
      Assign(SOrb, COrb_);
      TriangularMxm(SOrb, m_Scd, 'L', 'T');
      TriangularMxm(SOrb, m_Scd, 'L', 'N');
   //          Scale(ExchO, 1.);
      MakeRohfFockAo(m_RohfFock, m_FockA, m_FockB, m_Fock, m_ExchO, COrb_, SOrb, 0., 0., Mem);
//       xout << fmt::format("!! inside BuildRhfFock: called MakeRohfFockAo (HaveOrb={}, nOpen={})", m_HaveOrbitals, WfDecl.nOpen()) << std::endl;
   } else {
      Assign(m_RohfFock, m_Fock);
      m_FockMatrixIsLevelShifted = false;
      // ^- in the closed-shell case, RohfFock and Fock may be aliased.
      // However, assign with factor 1 does nothing for aliased in/out matrices.
//       xout << fmt::format("!! inside BuildRhfFock: bypassed MakeRohfFockAo (HaveOrb={}, nOpen={})", m_HaveOrbitals, WfDecl.nOpen()) << std::endl;
   }
   // ^-- note @ level shifts: those are now applied right before DIIS, after
   //     computing orbital gradient, energies, etc. And this only happens
   //     if there is another iteration to go; if the calculation is converged,
   //     the exit point is before that. This obviates the need to ever remove
   //     level shifts from a Fock matrix.

   IR_SUPPRESS_UNUSED_WARNING(InRefinementIt);
}

void FHfMethod::UpdateOrbitalEnergies(FMemoryStack &Mem)
{
   size_t nAo = _GetNumBf();
   FStackMatrix
      FockMo(nAo, nAo, &Mem);
   ChainMxm(FockMo, Transpose(m_COrb), m_RohfFock, m_COrb, Mem);
   for (size_t i = 0; i < nAo; ++ i)
      m_Ew[i] = FockMo(i,i);
   // ^- TODO: can't this be done without the full size mxms?
}

// void FHfMethod::MakeRohfFockAo1(FMatrixView FockRohf, FMatrixView FockA, FMatrixView FockB, FMatrixView Fock, FMatrixView ExchO, FMatrixView Orb_, FMemoryStack &Mem)
// {
//    // hm... what about the level shifts?
//    {
//       FStackMatrix
//          SOrb(nAo, nAo, &Mem);
//       Assign(SOrb, Orb_);
//       TriangularMxm(SOrb, Scd, 'L', 'T');
//       TriangularMxm(SOrb, Scd, 'L', 'N');
//    //          Scale(ExchO, 1.);
//       MakeRohfFockAo(RohfFock, FockA, FockB, Fock, ExchO, Orb_, SOrb, 0., 0., Mem);
//    }
// }


void FHfMethod::Run(FMemoryStack &Mem)
{
   size_t
      nAo = _GetNumBf();
   FStackMatrix
      OrbGrad(nAo, nAo, &Mem);
   FDiisState
      Diis(Mem, FDiisOptions(m_Options.MaxDiis));

   double
      LastEnergy = m_Energy,
      LastOrbGrad = 0.,
      CoreEnergy = 0.,
      BandEnergy = 0.;

   m_Timers.SetLevel(8);

   m_pLog->CheckStatus();
   m_pLog->WriteLine();
   if (m_Options.LevelShifts[0] != 0. || m_Options.LevelShifts[1] != 0.)
      m_pLog->Write(" {:<31}CLOSED = {:+8.2f}  ACTIVE =  {:+8.2f}", "Level shifts:", m_Options.LevelShifts[0], m_Options.LevelShifts[1]);
   if (m_Options.InitialDampingLevelShift != 0.) {
      m_pLog->Write(" {:<31}LSHIFT = {:+8.2f}  DECAY/IT ={: 8.4f}", "Initial damping:", m_Options.InitialDampingLevelShift, m_Options.InitialDampingLevelShiftDecay);
      if (m_Options.InitialDampingLevelShift > 0)
         m_pLog->Write(" Warning: initial damping level shift is postive; this accelerates initial updates instead of damping them.");
      if (1. - m_Options.InitialDampingLevelShiftDecay > 1.)
         m_pLog->Write(" Warning: damping decay is negative; this means that effective damping will INCREASE over iterations.");
   }
   m_pLog->Write(" {:<31}THRDEN = {:.2e}  THRGRAD = {:.2e}", "Convergence thresholds:", m_Options.ThrDen, m_Options.ThrOrb);
   m_pLog->Write("\n   ITER.     TOT.ENERGY    ENERGY CHANGE        DEN1     GRADIENT      TIME  DIIS");
   size_t
      iIt;
   bool
      Converged = false,
      UseRefineGrid = m_Options.UseRefineGrid && (m_Options.DftGridParamsRefine.fTargetAccuracy() < m_Options.DftGridParams.fTargetAccuracy()),
      InRefinementIt = false;
   double
      fInitialDampingLevelShift = m_Options.InitialDampingLevelShift;
   FTimer tIterations;
   for (iIt = 0; (iIt < m_Options.MaxIt) || InRefinementIt; ++ iIt) {
      FStackMatrix
//          T1(nAo, nAo, &Mem),
         OrbGrad(nAo, nAo, &Mem),
//          ExchC(nAo, nAo, &Mem),
         DenC(nAo, nAo, &Mem),
         DenO(nAo, nAo, &Mem);
      // make new orbitals
//       UpdateOrbitals(Orb, /*COccC, COccO,*/ DenC, DenO, Fock, ExchO, T1, Mem);
      UpdateOrbitals(m_COrb, m_Ew,/*COccC, COccO,*/ DenC, DenO, /*m_RohfFock,*/ /*T1,*/ Mem);
      // ^-- this computes orbitals and density matrices (i.e., DenC/DenO are outputs)


      // make new Fock matrix sets (both raw Fock/ExchO and RohfFock/FockA/FockB)
      BuildRhfFock(m_COrb, InRefinementIt, Mem);


      double
         fOrbGrad;
      ComputeOrbitalGradient(OrbGrad, &fOrbGrad, DenC, m_RohfFock, Mem);
      LastOrbGrad = fOrbGrad;
      AssembleEnergy(DenC, DenO, Mem);

      PrintTimingsAtLevel(2, TIMING_NewLineBeforeAndAfter);

      m_pLog->Write("{:6}{:18.8f}{:15.8f}{:15.8f}{:11.2e}{:10.2f}{:3}{:3}", (iIt+1), m_Energy,
         (m_Energy-LastEnergy), 0., fOrbGrad, (double)m_tTotal,
         Diis.nNextVec(), Diis.nLastDim());

      m_pLog->CheckStatus();

      Converged = (m_Options.MaxIt == 1) || InRefinementIt ||
                  (std::abs(m_Energy-LastEnergy) < m_Options.ThrDen &&
                   fOrbGrad < m_Options.ThrOrb);
      LastEnergy = m_Energy;


      if (Converged && (InRefinementIt || !UseRefineGrid)) {
         ComputeExtraEnergies(CoreEnergy, BandEnergy, DenC, DenO, Mem);
         StoreScfResultData1(DenC, DenO, Mem);
         break;
      }
      if (Converged) {
         // switch to refinement grid and do one more iteration.
         InRefinementIt = true;
         for (size_t iBuilder = 0; iBuilder < m_FockComponentBuilders.size(); ++ iBuilder)
            m_FockComponentBuilders[iBuilder]->SwitchToRefineGrid(m_Options.DftGridParamsRefine);
         // should we update orbitals again? I would say no, but it would be a bit
         // more hacky. So I just leave it like this atm.
      }

      if (1) {
         // apply level shifts: both constant ones and initial damping shifts.
         ApplyLevelShifts1(fInitialDampingLevelShift, Mem);
         // ^- Notes:
         //    - Initial dampening is currently disabled by default; but it did
         //      sometimes appear to help with geometry projections. Normally
         //      gradient goes up after first iteration. With this it does not.
         fInitialDampingLevelShift *= (1. - m_Options.InitialDampingLevelShiftDecay);
         if (std::abs(fInitialDampingLevelShift) < 0.05)
            fInitialDampingLevelShift = 0.;
      }

      { TIME_SECTION(&m_Timers, 0x104, "SCF (diis)");
         FDiisTarget1 RT(m_RohfFock.pData, m_FockSize, OrbGrad.pData, OrbGrad.GetStridedSize()); // Fock.GetStridedSize()
         Diis(RT);
      }
   }
   iIt += 1; // currently: iteration number in which the computation was stopped.
   m_pLog->WriteLine();
   m_pLog->WriteTiming("iteration", (double)tIterations, iIt);
   m_pLog->WriteLine();

   PrintTimingsAtLevel(1, TIMING_Reset | TIMING_NewLineAfter);

   if (!Converged) {
      char const *pHfOrKs = (m_Options.XcFunctionalName.empty()? "Hartree-Fock" : "Kohn-Sham");
      m_pLog->EmitWarning(fmt::format(" {} failed to converge in {} iterations.", pHfOrKs, iIt));
      if (!m_Options.IgnoreConvergenceError)
         throw ct::FConvergenceError(pHfOrKs, iIt, LastOrbGrad);
   }

   if (1) {
      m_pLog->WriteResult("Nuclear repulsion energy", m_pAtoms->NuclearRepulsionEnergy());
      m_pLog->WriteResult("1-electron energy", CoreEnergy);
      m_pLog->WriteResult("Band energy", BandEnergy);
//       m_pLog->WriteResult("band+one+rep", .5*BandEnergy + .5*CoreEnergy + m_pAtoms->NuclearRepulsionEnergy());
//    ^- interesting.. with my rohf fock this also works in open-shell case (with BandEnergy = Dot2(Fock, Denc) only).
      m_pLog->WriteLine();
   }
   for (size_t iBuilder = 0; iBuilder < m_FockComponentBuilders.size(); ++ iBuilder)
      m_FockComponentBuilders[iBuilder]->PrintEnergyContribs();

   if (m_Options.XcFunctionalName.empty())
      m_pLog->WriteResult("Hartree-Fock energy", LastEnergy);
   else
      m_pLog->WriteResult("Kohn-Sham energy", LastEnergy);
   m_pLog->Flush();

   assert(!m_FockMatrixIsLevelShifted);
   if (1) {
      // update the orbital energies---make them as diagonal values of the CURRENT Fock matrix.
      // The ones we currently have stored are from the PREVIOUS Fock matrix.
      UpdateOrbitalEnergies(Mem);
   }

   if (m_Options.Print[SCFPRINT_OrbitalEnergies])
      PrintOrbitalEnergies(m_Ew, "");
   m_pLog->Flush();

   if (m_pResult) {
      // other stuff set outside. this only here because the variables are local...
      m_pResult->nIt = iIt;
      m_pResult->Converged = Converged;
   }


   if (!m_Options.PropertyTasks.m_1eTaskList.empty()) {
      FTimer tProperties1e;
      EvalProperties1e(Mem);
      m_pLog->WriteTiming("1-electron properties", (double)tProperties1e);
   }

}


void FHfMethod::ComputeExtraEnergies(double &CoreEnergy, double &BandEnergy, FMatrixView DenC, FMatrixView DenO, FMemoryStack &Mem)
{
   // compute some extra energies for the final printing.
   // It's done in the last iteration, because afterwards the density matrices
   // are no longer around. Also, we don't always compute those because
   // we generally do not need all of these energy components.

   // BandEnergy = Dot(DenC, m_RohfFock);
   CoreEnergy = Dot(DenC, m_CoreH);
   BandEnergy = Dot(DenC, m_Fock);
   if (m_WfDecl.nOpen() != 0) {
      BandEnergy += Dot(DenO, m_Fock);
      BandEnergy -= Dot(DenO, m_ExchO);
   }
   // ^- hm.. is that right? does that evaluate to the sum of eigenvalues*occ
   //    like the thing in AssembleEnergy does?
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FHfMethod::StoreScfResultData1(FMatrixView DenC, FMatrixView DenO, FMemoryStack &Mem)
{
   if (m_pResult) {
      size_t nAo = _GetNumBf();
      m_pResult->Rdm.Reshape(nAo,nAo);
      Assign(m_pResult->Rdm, DenC);
      // ^- note: that already includes DenO; it's really the total density (A+B)
      if (false) {
         if (m_WfDecl.nOpen() != 0)
            Assign(m_pResult->RdmO, DenO);
         Assign(m_pResult->RohfFock, m_RohfFock);
         // ^-- FIXME: remove these two (and maybe the RDM, too.)
         m_pLog->Write("!CtRhf.cpp//TODO: remove storage of extra output RDM/FOCK matrices. RohfFock.Shape = {}  m_pResult->RohfFock = {}", m_RohfFock.Shape(), m_pResult->RohfFock.Shape());
      }
   }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FHfMethod::ComputeOrbitalGradient(FMatrixView OrbGrad, double *pfOrbGrad, FMatrixView DenC, FMatrixView Fock,  FMemoryStack &Mem)
{
   size_t nAo = _GetNumBf();
   assert(OrbGrad.IsAssignedAndHasShape(nAo, nAo));

   // evaluate orbital gradient:  S^{-1/2} F D S^{1/2} - h.c.
   { TIME_SECTION(&m_Timers, 0x103, "SCF (orb. grad.)");
      FStackMatrix
         T1(nAo, nAo, &Mem);
//          ChainMxm(T1, Smh, Fock, DenC, Overlap, Smh, Mem);
      Mxm(T1, Fock, DenC);
      // ^- note: this one goes with RohfFock.
//          Mxm(T1, Fock, DenC);
      TriangularMxm(T1, m_Scd, 'R', 'N');
//          TriangularMxm(T1, Scd, 'R', 'T'); // <- this one for F D S
      TriangularSolve(T1, m_Scd, 'L', 'N'); // <- this one for Smh F D Sh
      Assign(OrbGrad, T1);
      Add(OrbGrad, Transpose(T1), -1.);
      if (pfOrbGrad) {
         // compute RMSD of orb gradient from zero.
         *pfOrbGrad = Dot(OrbGrad, OrbGrad);
         *pfOrbGrad /= nAo;
         *pfOrbGrad = std::sqrt(*pfOrbGrad);
      }
   }
}


void FHfMethod::AssembleEnergy(FMatrixView DenC, FMatrixView DenO, FMemoryStack &Mem)
{
   m_Energy = m_pAtoms->NuclearRepulsionEnergy();
   size_t nAo = _GetNumBf();

   if (!m_Options.UseBandStructureEnergyFormula) {
      // compute energy in regular way... individually summing all contributions
      // (t, vnuc, j, etc)
      m_Energy += Dot(DenC, m_CoreH);
   } else {
      // compute energy from sum of orbital energies multiplied by occupations
      // (in DFT circles called ``the band structure term'') and then subtracting
      // the contributions which are doubly-counted by doing so.
      //
      // (this is also a cgk-special; it's an experimental option related to
      // tight-binding method development)
      assert(m_WfDecl.nOpen() == 0); // it's experimental. No open-shell supported.
      if (0) {
         // compute sum-of-orbital-energies from orbitals used to construct the Fock matrix
         m_Energy += Dot(DenC, m_RohfFock);
      } else {
         // Use *new* density for that, ...
         // Or, alternatively, use actual orbital eigenvalues multiplied with their occupations?
         // hack hack hack... (can't we just use the RDM for this? Maybe try...)
         FStackMatrix
            T1(nAo, nAo, &Mem);
         Assign(T1, m_RohfFock);
         TriangularSolve(T1, m_Scd, 'R', 'T');
         TriangularSolve(T1, m_Scd, 'L', 'N');
         Diagonalize(T1, m_Ew.pData, Mem);
         TriangularSolve(T1, m_Scd, 'L', 'T');

         double
            fBandEnergy = 0.;
         size_t
            nOccC = m_WfDecl.nElecB(),
            nOccO = m_WfDecl.nOpen();
         for (size_t i = 0; i != nOccC; ++ i)
            fBandEnergy += 2*m_Ew[i];
         for (size_t i = nOccC; i != nOccC+nOccO; ++ i)
            fBandEnergy += 1*m_Ew[i];
         m_Energy += fBandEnergy;
      }
   }

   if (m_WfDecl.nOpen() != 0) {
      m_Energy += 0.;
   }
   // ^- well... that seems useful.

   // collect energy contributions from the Fock builders. These include
   // two-electron energy terms, xc energy terms, and possibly other stuff
   // (e.g., the D3(BJ) dispersion correction goes in like this, and
   // continuum solvation drivers may be implementable like this, too)
   for (size_t iBuilder = 0; iBuilder < m_FockComponentBuilders.size(); ++ iBuilder)
      m_Energy += m_FockComponentBuilders[iBuilder]->Energy;

   IR_SUPPRESS_UNUSED_WARNING(DenO); //< actual open-shell contributions come from the fock builders.
}


static void PrintOrbitalEnergyTable(FLog *pLog, FMatrixView OrbEw, size_t iFirst, size_t iLast, size_t nOccA, size_t nOccB, std::string const &Desc, std::string const &Caption, unsigned Flags)
{
   pLog->w << "\n ";
   if (!Desc.empty())
      pLog->Write("{} // {}:", Desc, Caption);
   else
      pLog->Write("{}:", Caption);
//    pLog->Write("\n {}:", Caption);
   std::string
      OrbOcc_Full = "ab",
      OrbOcc_Empty = "__",
      OrbOcc = "??";
   if (!bool(Flags & FHfMethod::ORBEIG_AlphaOn))
      OrbOcc_Full[0] = OrbOcc_Empty[0];
   if (!bool(Flags & FHfMethod::ORBEIG_BetaOn))
      OrbOcc_Full[1] = OrbOcc_Empty[1];
   size_t
      nValuesPerRow = 6; // 8 looks okay in text files, but for 80 col wrap it's not so good.

   for (size_t iOrb = iFirst; iOrb < iLast; iOrb += nValuesPerRow) {
      size_t n = std::min(nValuesPerRow, iLast-iOrb);
      pLog->w << "\n ";
      for (size_t iCol = 0; iCol < n; ++ iCol) {
         size_t iOrb1 = iCol+iOrb;

         OrbOcc[0] = OrbOcc_Empty[0];
         OrbOcc[1] = OrbOcc_Empty[1];
         if (iOrb1 < nOccA)
            OrbOcc[0] = OrbOcc_Full[0];
         if (iOrb1 < nOccB)
            OrbOcc[1] = OrbOcc_Full[1];
         pLog->w.write(" {:8}:{:2} ", (iOrb1+1), OrbOcc);
      }
      pLog->w << "\n ";
      for (size_t iCol = 0; iCol < n; ++ iCol)
         pLog->w.write(" {:12.6f}", OrbEw[iCol+iOrb]);
      pLog->w << "\n";
   }
}


void FHfMethod::PrintOrbitalEnergies(FMatrixView Ew, std::string const &Desc, unsigned Flags)
{
   if (m_Options.Print[SCFPRINT_OrbitalEnergies]) {
      size_t
         nAo = _GetNumBf(),
         iOrb = m_pAtoms->nCoreElec()/2,
         nMinAo = (size_t)m_pMinBasis->nFn(),
         nPrint = std::min(nMinAo, 2 + (size_t)m_WfDecl.nElecA());
      if (m_Options.Print.ExtraDetail(SCFPRINT_OrbitalEnergies))
         // print all valence-level orbital energies if extra-level print is on
         nPrint = (size_t)m_pMinBasis->nFn();
      if (iOrb + nPrint > nAo)
         nPrint = nAo - iOrb;
      if (iOrb > 0 && m_Options.Print.ExtraDetail(SCFPRINT_OrbitalEnergies))
         PrintOrbitalEnergyTable(m_pLog, Ew, 0, iOrb, m_WfDecl.nElecA(), m_WfDecl.nElecB(), Desc, "Core orbital energies", Flags);

      PrintOrbitalEnergyTable(m_pLog, Ew, iOrb, nPrint, m_WfDecl.nElecA(), m_WfDecl.nElecB(), Desc, "Valence orbital energies", Flags);
   }
}

void FHfMethod::PrintTimingsAtLevel(int iTargetLevel, unsigned Flags)
{
   if (m_Options.iTimingLevel >= iTargetLevel) {
      std::stringstream str;
      if (0 != (Flags & TIMING_NewLineBefore))
         str << "\n";
      m_Timers.PrintReport(str);
      if (0 != (Flags & TIMING_NewLineAfter))
         m_pLog->Write(str.str());
      else
         m_pLog->WriteNoNl(str.str());
      if (0 != (Flags & TIMING_Reset)) {
         m_Timers.Reset();
      }
   }
}


void FHfMethod::GetMoSubspaceOffs(size_t &ibeg, size_t &iend, char Space)
{
   if (Space == 'c') {
      // space of occupied closed-shell-orbitals. These are the first in COrb.
      ibeg = 0;
      iend = m_WfDecl.nElecB();
   } else if (Space == 'a') {
      // space of active orbitals (in ROHF/ROKS, singly occupied wit alpha-electrons)
      ibeg = m_WfDecl.nElecB();
      iend = m_WfDecl.nElecA();
   } else if (Space == 'e') {
      // space of external orbitals (unoccupied in all reference configurations)
      ibeg = m_WfDecl.nElecA();
      iend = _GetNumBf(); // = nAo = pOrbBasis->nFn()
   } else {
      throw std::runtime_error(fmt::format("orbital subspace '{}' not recognized in FScfMethod::GetMoSubspaceOffs().", Space));
   }
}


void FHfMethod::ProjectG(FMatrixView ProjectedG, char const *pSpaces, FMatrixView InputG, FMatrixView C, FMatrixView SC, double Factor, FMemoryStack &Mem)
{
   // evaluate
   //  gp := Factor * |s1><s1| g |s2><s2|
   // where the space projectors s1 and s2 are supplied via pSpaces[0] and pSpaces[1]:
   //    'c': closed shell orbitals (nElecB),
   //    'a': active orbitals (nElecA-nElecB)
   //    'e': external orbitals (nBf - nElecA)
   // Parameters:
   // - InputG is the matrix to project. Output will be stored in ProjectedG.
   // - C is the orbital matrix
   // - SC is S*C, where S is the overlap matrix and C the orbital matrix.
   size_t ibeg0, iend0;
   size_t ibeg1, iend1;

   GetMoSubspaceOffs(ibeg0, iend0, pSpaces[0]);
   GetMoSubspaceOffs(ibeg1, iend1, pSpaces[1]);

   if ((iend0-ibeg0 == 0) || (iend1-ibeg1 == 0)) {
      // no space here. Delete g and return.
      ProjectedG.Clear();
      return;
   }
//    std::cout << fmt::format("project-g: {}  [{}:{}]  [{}:{}]", pSpaces, ibeg0, iend0, ibeg1, iend1) << std::endl;

   size_t
      nAo = _GetNumBf();
   FMatrixView
      C0 = Select(C, 0, ibeg0, nAo, iend0-ibeg0),
      C1 = Select(C, 0, ibeg1, nAo, iend1-ibeg1),
      SC0 = Select(SC, 0, ibeg0, nAo, iend0-ibeg0),
      SC1 = Select(SC, 0, ibeg1, nAo, iend1-ibeg1);
   assert(C0.nCols == SC0.nCols);
   assert(C1.nCols == SC1.nCols);


   // transform g from nBf^2 to the MO subspace asked for.
   FStackMatrix
      G_MO(C0.nCols, C1.nCols, &Mem);
   ChainMxm(G_MO, Transpose(C0), InputG, C1, Mem);
   // transform back to covariant AO basis (with lower labels _{\mu,\nu}, via via S*C).
   ChainMxm(ProjectedG, SC0, G_MO, Transpose(SC1), Mem);
   Scale(ProjectedG, Factor);
   if (pSpaces[0] != pSpaces[1])
      // this should account for the fact that we are really
      // only computing half of the matrix.
      Symmetrize(ProjectedG);
}


void FHfMethod::MakeRohfFockAo(FMatrixView FockRohf, FMatrixView FockA, FMatrixView FockB, FMatrixView Fock, FMatrixView ExchO, FMatrixView C, FMatrixView SC, double LevelShiftC, double LevelShiftO, FMemoryStack &Mem)
{
   size_t
      nAo = _GetNumBf();
   if (m_WfDecl.nOpen() == 0) {
      // closed-shell calculation. Nothing to project. Leave alone this->Fock and this->ExchO.
      assert_rt(FockRohf.pData == Fock.pData);
      return;
//       ApplyLevelShifts(1., COccC, COccO, Mem);
   } else {
//       ExchO.Clear();
//       Scale(ExchO, -1.0);
      assert_rt(FockRohf.pData != Fock.pData);
      Assign(FockRohf, Fock);
      FStackMatrix
         gcc_ao(nAo, nAo, &Mem),
         gec_ao(nAo, nAo, &Mem),
         gac_ao(nAo, nAo, &Mem);
      ProjectG(gcc_ao, "cc", ExchO, C, SC, 1., Mem);
      ProjectG(gec_ao, "ec", ExchO, C, SC, 2., Mem);
      ProjectG(gac_ao, "ca", ExchO, C, SC, 2., Mem);

      if (m_RohfAlgo == FHfOptions::ROHFALGO_TwoStep) {
         Assign(FockB, FockRohf);
         Add(FockB, ExchO, -1.);
         Add(FockB, gec_ao, +1.);
      }
      Add(FockRohf, ExchO, 1.);
      Add(FockRohf, gec_ao, -1.);
      if (m_RohfAlgo == FHfOptions::ROHFALGO_TwoStep)
         Assign(FockA, FockRohf);
      // both FockA and FockRohf now:
      //        cls   act   ext
      //      +-----+-----+-----+
      // cls  |  A  |  A  | A+B |   (note: A+B means A+B/2 etc)
      //      +-----+-----+-----+
      // act  |  A  |  A  |  A  |
      //      +-----+-----+-----+
      // ext  | A+B |  A  |  A  |
      //      +-----+-----+-----+
      Add(FockRohf, gcc_ao, -1.);
      Add(FockRohf, gac_ao, -2.);
      // FockRohf now:
      //        cls   act   ext
      //      +-----+-----+-----+
      // cls  | A+B |  B  | A+B |
      //      +-----+-----+-----+
      // act  |  B  |  A  |  A  |
      //      +-----+-----+-----+
      // ext  | A+B |  A  |  A  |
      //      +-----+-----+-----+
   }
   IR_SUPPRESS_UNUSED_WARNING(LevelShiftC);
   IR_SUPPRESS_UNUSED_WARNING(LevelShiftO);
}




void FHfMethod::UpdateOrbitals(FMatrixView Orb, FMatrixView Ew,/*FMatrixView &COccC, FMatrixView &COccO,*/
   FMatrixView &DenC, FMatrixView &DenO, /*FMatrixView const RohfFock,*/ /*FMatrixView &T1,*/ FMemoryStack &Mem)
{
   TMemoryLock<double> pFreeMe0(0, &Mem);
   m_Timers.Resume(0x101, "SCF (orbitals)");
   size_t
      nAo = _GetNumBf();
   if (0) {
      FStackMatrix
         T1(nAo, nAo, &Mem);
      // transform Fock matrix to smh basis and diagonalize to get new orbitals.
      ChainMxm(T1, Transpose(m_Smh), m_RohfFock, m_Smh, Mem);
      Diagonalize(T1, Ew.pData, Mem);
      Mxm(Orb, m_Smh, T1);
   } else if (!m_HaveOrbitals) {
//       xout << fmt::format("!! inside UpdateOrbitals: HaveOrbitals = false path --> diag(RohfFock).") << std::endl;
      Assign(Orb, m_RohfFock);
      TriangularSolve(Orb, m_Scd, 'R', 'T');
      TriangularSolve(Orb, m_Scd, 'L', 'N');
      Diagonalize(Orb, Ew.pData, Mem);
      TriangularSolve(Orb, m_Scd, 'L', 'T');
   } else {
//       xout << fmt::format("!! inside UpdateOrbitals: HaveOrbitals = true path.") << std::endl;
      // transform Fock matrix to schmidt basis and diagonalize to get new orbitals.
      if (m_RohfAlgo == FHfOptions::ROHFALGO_TwoStep)
         // first make new occupied orbitals (=alpha orbital), then diagonalize
         // beta fock matrix in basis of alpha orbitals.
         Assign(Orb, m_FockA);
      else
         // one step algorithm: make orbitals directly.
         Assign(Orb, m_RohfFock);
      TriangularSolve(Orb, m_Scd, 'R', 'T');
      TriangularSolve(Orb, m_Scd, 'L', 'N');
      Diagonalize(Orb, Ew.pData, Mem);
      TriangularSolve(Orb, m_Scd, 'L', 'T');
      if (m_RohfAlgo == FHfOptions::ROHFALGO_TwoStep) {
         size_t
            nElecA = m_WfDecl.nElecA();
         FMatrixView
            COrbAlpha = Select(Orb,0,0,nAo,nElecA);
         FStackMatrix
            Tmp(nAo, nElecA, &Mem),
            FockOccBeta(nElecA, nElecA, &Mem);
         ChainMxm(FockOccBeta, Transpose(COrbAlpha), m_FockB, COrbAlpha, Mem);
         Diagonalize(FockOccBeta, Ew.pData, Mem);
         Assign(Tmp, COrbAlpha);
         Mxm(COrbAlpha, Tmp, FockOccBeta);
      }
   }


   RebuildDensityMatrices(DenC, DenO, Mem);

   m_HaveOrbitals = true;
   m_Timers.Pause(0x101, "SCF (orbitals)");

//    IR_SUPPRESS_UNUSED_WARNING(ExchO);
}


void FHfMethod::RebuildDensityMatrices(FMatrixView &DenC, FMatrixView &DenO, FMemoryStack &Mem)
{
   {
      TMemoryLock<double> pFreeMe1(0, &Mem);
      // get occupied orbitals with absorbed occupation numbers (allocated on Mem)
      FMatrixView COccC, COccO;
      ExtractOrbitals(COccC, COccO, m_COrb, ORBEXTRACT_AbsorbOcc, Mem);

      // form density matrices.
      SyrkNT(DenC, COccC);
      SyrkNT(DenO, COccO);
      if (m_WfDecl.nOpen() != 0)
         Add(DenC, DenO);
   }
}


void FHfMethod::ExtractOrbitals(FMatrixView &COccC, FMatrixView &COccO, FMatrixView const &COrbIn, unsigned Flags, FMemoryStack &Mem)
{
   if ((Flags & ~(ORBEXTRACT_AllowAliasing | ORBEXTRACT_AbsorbOcc)) != 0)
      throw std::runtime_error("ExtractOccOrbitals: AllowAliasing and AbsorbOcc cannot be set at the same time");

   size_t
      nClosed = m_WfDecl.nElecB(),
      nOpen = m_WfDecl.nOpen();
   size_t
      nAo = _GetNumBf();

   if (COrbIn.nRows != nAo || (COrbIn.nCols < nClosed + nOpen))
      throw std::runtime_error(fmt::format("ExtractOrbitals: dimension mismatch (COrb.Shape = {} vs. nAo = {}, nClosed = {}, nOpen = {})", COrbIn.Shape(), nAo, nClosed, nOpen));

   FMatrixView
      COrbIn_OccC = Select(COrbIn, 0,0, nAo,nClosed),
      COrbIn_OccO = Select(COrbIn, 0,nClosed, nAo,nOpen);
   bool
      AbsorbOcc = bool(Flags & ORBEXTRACT_AbsorbOcc);
   if (bool(Flags & ORBEXTRACT_AllowAliasing)) {
      COccC = COrbIn_OccC;
      COccO = COrbIn_OccO;
      assert_rt(!AbsorbOcc);
   } else {
      COccC = MakeStackMatrix(nAo, nClosed, Mem);
      COccO = MakeStackMatrix(nAo, nOpen, Mem);
      Assign(COccC, COrbIn_OccC, AbsorbOcc ? std::sqrt(2.) : 1.);
      Assign(COccO, COrbIn_OccO, 1.);
   }
}





// makes strings describing cartesian moments; e.g., like 'xxy' for input (2,1,0)
std::string sCartMom(TVector3<unsigned> CartMom)
{
   if (CartMom[0] == 0 && CartMom[1] == 0 && CartMom[2] == 0)
      return "1";
   std::stringstream
      ss;
   for (unsigned ixyz = 0; ixyz < 3; ++ ixyz)
      for (size_t i = 0; i < CartMom[ixyz]; ++i)
         ss << "xyz"[ixyz];
   return ss.str();
}


// predicate for sorting cartesian moments into the same order which molpro does
struct FReverseLexiOrderPred {
   FReverseLexiOrderPred(unsigned n_) : n(ptrdiff_t(n_)) {}
   ptrdiff_t n;
   typedef TVector3<unsigned> T;
   ptrdiff_t iKey(T const &v) const { return -(n*n*ptrdiff_t(v[0]) + n*ptrdiff_t(v[1]) + v[0]); }
   bool operator () (T const &va, T const &vb) const {
      return iKey(va) < iKey(vb);
   }
};


// ...should this maybe even go into FAtomSet? Lots of other property stuff
// *is* in there...
FAtomSetPropertyPtr EvalProperty1e(FLog *pLog, FAtomSet const *pAtoms, FMatrixView const &Rdm, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, FMpt1e const *pMpt1e, FMemoryStack &Mem)
{
   // make && time microscf -t '!DF-RHF/def2-TZVPP scf{eval: {s;dm;sm;tm}; thr-orb:1e-12}' -g hessian/_hess_ch2foh-optg.xyz
   // ^- comp to molpro: seem okay in principle, but for the higher
   //    moments the agreement is numerically not so great... maybe 6 digits or so.
   //    I seem to remember that CtInt1e had some numerical stability problems
   //    in some situations... maybe that is it? Look another day. Good enough
   //    for now.
   //
   // Other TODO:
   // - Generalize and refactor this thing so that it can use the
   //   functionality also for other operators. And allow reformatting
   //   the tables. E.g., the ELEC/NUCL split could be done optional,
   //   on the other hand, some properties probably want the center coords
   //   printed.
   // - could put such things either into general options, or object-specific
   //   ones.
   // - could make a single big table, and add the extra units or other
   //   extra stuff to end, without table caption. (e.g., 0.648518 DEBYE,
   //   instead of just 0.648518; for others this could also mean info like centers)
   // - ATM the tables would also be unnecessarily hard to parse. Maybe
   //   do something about that?
   // Refactoring TODO:
   // - We actually have a file CtRhfProp.cpp... is there a reason this isn't there?
   std::string
      Desc;
   unsigned
      // note: for non-center dependent operators, set nCenters to 1.
      nCenters = 0, nStates = 1, nComps = 0;
   FAtomSetProperty::FDataArray
      // nComps x nCenters array of numerical data values for the property
      ResultData;

   // collect list of actual centers to use for current operator and current atom set.
   FMpt1e::FPointList
      Centers;
   Centers.reserve(pMpt1e->m_iCenters.size() + pMpt1e->m_vCenters.size());
   for (size_t iCen = 0; iCen != pMpt1e->m_iCenters.size(); ++ iCen) {
      int iAt = pMpt1e->m_iCenters[iCen];
      if (iAt < 0 && size_t(iAt) >= pAtoms->size())
         throw std::runtime_error(fmt::format("FHfMethod::EvalProperty1e(): while processing property of type {}: "
            "atomic center index {} is invalid (nAt={})",
            int(pMpt1e->m_PropertyType), iAt, pAtoms->size()));
      // ^- we could possibly just ignore them, but it might lead to unexpected
      // data misattribution errors, because we use these guys for output
      // property data indexing...
      Centers.push_back((*pAtoms)[iAt].vPos);
   }
   // add the centers which were explicitly given.
   Centers.insert(Centers.end(), pMpt1e->m_vCenters.begin(), pMpt1e->m_vCenters.end());
   // ...and if there is still no center given anywhere, just add the coordinate
   // origin to the list (we need at least one center, even if unused, because
   // output data is indexed as nComps x nCenters array).
   if (Centers.empty())
      Centers.push_back(FVector3(0.,0.,0.));
   nCenters = Centers.size();

   // get basis dimensions (maybe also get raw basis refs?)
   size_t
      nAoA = RowBasis.nFn(),
      nAoB = ColBasis.nFn();

   // TODO: should probably rephrase this in terms of IrMeta2i thingies...
   // (and add the other operators we actually already have...)
   // But atm I just want to get the dipole moments working.
//    if (pMpt1e->m_PropertyType == PROPTYPE_DipoleMoments) {
   if (pMpt1e->m_PropertyType == FMpt1e::PROPTYPE_ChargeMoments_Cartesian) {
      if (Centers.size() != 1) {
         throw std::runtime_error(fmt::format("FHfMethod::EvalProperty1e(): while processing property of type {} (order {}): "
            "charge distribution moments are not supposed to have exactly one expansion origin, but this one has {}",
            int(pMpt1e->m_PropertyType), pMpt1e->m_Order, pAtoms->size()));
      }
      assert(Centers.size() == 1);
      int
         iCenter = 0;
      FVector3
         vExpPoint = Centers[iCenter];
      double
         // conversion factor to use if this entry is meant to be
         // output also with an additional physical different unit
         // except A.U.
         ToExtraUnit = 0;
      char const
         *pNumFmt = "{:12.6f}",
         *pNumFmt1 = "{:16.8f}";
//          *pNumFmt1 = "{:16.8e}";
      std::string
//          pRowFmt = fmt::format("{{:6}} {{:>4}}  {{:^18}}  {0}  {0}  {0}", pNumFmt),
         pRowFmt = fmt::format("{{:6}} {{:>4}}  {{:^12}}  {{:^12}}  {1}  {0}  {0}", pNumFmt, pNumFmt1),
         pExtraUnitFmt;
      if (pLog) {
         pLog->Write(" Cartesian moments of total charge <(r-cen)^p rho(r)>:  ORDER = {}  CEN = {: .4f} {: .4f} {: .4f}",
            pMpt1e->m_Order, vExpPoint[0], vExpPoint[1], vExpPoint[2]);
         std::string
            sCaption = fmt::format("\n   CEN  COMP  {:^12}  {:^12}  {:^16}  {:^12}  {:^12}", "EXPV.STATES", "OPERATOR", "TOTAL/[A.U.]", ">FROM.ELEC.", ">FROM.NUCL."),
            sExtraUnit;
         if (pMpt1e->m_Order == 1) {
            ToExtraUnit = ToDebye;
            sExtraUnit = "TOT./[DEBYE]";
         }
         if (pMpt1e->m_Order >= 2) {
            ToExtraUnit = std::pow(ToAng, pMpt1e->m_Order);
            sExtraUnit = fmt::format("TOT./[e*A^{}]", pMpt1e->m_Order);
//             sCaption += fmt::format(" |  {:^12}", fmt::format("TOT./[e*A^{}]", pMpt1e->m_Order));
         }
         if (!sExtraUnit.empty())
            sCaption += fmt::format(" |  {:^12}", sExtraUnit);
         pLog->Write(sCaption);
      }
      if (ToExtraUnit != 0. && pExtraUnitFmt.empty())
         pExtraUnitFmt = fmt::format("  |  {0}", pNumFmt);

      Desc = fmt::format("total charge moments (cartesian, order {})", pMpt1e->m_Order);
      // define list of cartesian components to make. For simplicity, we will
      // make them one-by-one.
      // (although Make1eIntegralMatrixMultiComp would be faster, and the IR
      // meta kernel for the charge and dipole integrals also definitely...)
      std::vector<TVector3<unsigned>>
         CartMoments;
      if (pMpt1e->m_Order == 0) {
         // 1e-op: 1 (overlap) --> total nuclear charge - total electron charge
         CartMoments.push_back(TVector3<unsigned>(0,0,0));
      } else if (pMpt1e->m_Order == 1) {
         // 1e-op: x/y/z position--> nuclear charge dipole moment - electron distribution dipole moment
         CartMoments.push_back(TVector3<unsigned>(1,0,0));
         CartMoments.push_back(TVector3<unsigned>(0,1,0));
         CartMoments.push_back(TVector3<unsigned>(0,0,1));
      } else if (pMpt1e->m_Order == 2) {
         // 1e-op: xx/yy/zz/xy/xz/yz position --> nuclear charge dipole moment - electron distribution dipole moment
         CartMoments.push_back(TVector3<unsigned>(2,0,0));
         CartMoments.push_back(TVector3<unsigned>(0,2,0));
         CartMoments.push_back(TVector3<unsigned>(0,0,2));
         CartMoments.push_back(TVector3<unsigned>(1,1,0));
         CartMoments.push_back(TVector3<unsigned>(1,0,1));
         CartMoments.push_back(TVector3<unsigned>(0,1,1));
      } else if (pMpt1e->m_Order >= 2) {
         unsigned n = pMpt1e->m_Order;
         for (unsigned iz = 0; iz <= n; ++ iz) {
            for (unsigned iy = 0; iy <= n - iz; ++ iy) {
               unsigned ix = n - iz - iy;
               CartMoments.push_back(TVector3<unsigned>(ix,iy,iz));
            }
         }
         std::sort(CartMoments.begin(), CartMoments.end(), FReverseLexiOrderPred(n));
      } else {
         throw std::runtime_error(fmt::format("EvalProperty1e(): internal error. negative property order '{}'?", pMpt1e->m_Order));
      }

      nComps = CartMoments.size();
      for (unsigned iComp = 0; iComp < CartMoments.size(); ++ iComp) {
         TVector3<unsigned>
            CartPow = CartMoments[iComp];
         double
            fValueElec = 0.;
         {
            // compute electronic contribution
            FStackMatrix
               OpMat(nAoA, nAoB, &Mem);
            OpMat.Clear();
            pAtoms->MakeCartesianMomentMatrix(OpMat, RowBasis, ColBasis, CartPow, vExpPoint, Mem);
            AssertSameShape1(OpMat, Rdm);
            fValueElec = -1 * Dot(OpMat, Rdm);
            // ^- (-1) is for negative electron charge.
         }
         double
            fValueNuc = 0.;
         {
            // compute nuclear contribution
            for (size_t iAt = 0; iAt != pAtoms->size(); ++ iAt) {
               FAtom const *pAt = &(*pAtoms)[iAt];
               // get atom position relative to origin of multipole expansion
               FVector3 v = pAt->vPos - vExpPoint;
               fValueNuc += pAt->NuclearCharge() *
                  std::pow(v[0], CartPow[0]) * std::pow(v[1], CartPow[1]) * std::pow(v[2], CartPow[2]);
            }
         }
         double
            fValueTot = fValueElec + fValueNuc;
         if (pLog) {
            std::string
               sRow = fmt::format(pRowFmt, iCenter, iComp, "[0|.|0]", fmt::format("  {}*rho(r)", sCartMom(CartPow)), fValueTot, fValueElec, fValueNuc);
            if (ToExtraUnit != 0.)
               sRow += fmt::format(pExtraUnitFmt, fValueTot * ToExtraUnit);
            pLog->Write(sRow);
         }
         // store result we got for current component.
         ResultData.push_back(fValueTot);
      }
   } else {
      throw std::runtime_error(fmt::format("FHfMethod::EvalProperty1e(): encountered unrecognized property (type {})", int(pMpt1e->m_PropertyType)));
   }

   // absorb prefactor we were asked to include.
   for (size_t i = 0; i != ResultData.size(); ++ i)
      ResultData[i] *= pMpt1e->m_Factor;
   return FAtomSetPropertyPtr(new FAtomSetProperty(Desc, nCenters, nStates, nComps, ResultData));
}



void FHfMethod::EvalProperties1e(FMemoryStack &Mem)
{
   m_pLog->Write("\n Computing static property integrals.\n");
   size_t
      nAo = _GetNumBf();
   // get refs to occupied orbitals
   size_t
      nOccC = m_WfDecl.nElecB(),
      nOccO = m_WfDecl.nOpen(),
      nOcc = nOccC + nOccO;
   assert(nOcc == m_WfDecl.nElecA());
   // note: no absorbed occupation numbers here (...at least in theory...).
   FMatrixView
      COccC = Select(m_COrb, 0, 0, nAo, nOccC),
      COccO = Select(m_COrb, 0, nOccC, nAo, nOccO);
   // build 1-particle reduced density matrix.
   FStackMatrix
      Rdm(nAo, nAo, &Mem);
   SyrkNT(Rdm, COccC, 2.0, 0);
   SyrkNT(Rdm, COccO, 1.0, MXM_Add);

   // go through all the 1e properties listed and evaluate them
   FAtomSet::FPropertyList
      NewProperties;
   for (size_t iTask1e = 0; iTask1e != m_Options.PropertyTasks.m_1eTaskList.size(); ++ iTask1e) {
      FMpt1e const
         *pMpt1e = &*m_Options.PropertyTasks.m_1eTaskList[iTask1e];
      FLog
         // get ref to output log, unless we were asked to skip printing this.
         *pLog = pMpt1e->IsSilent() ? 0 : m_pLog;
      FAtomSetPropertyPtr
         pPropertyResult = EvalProperty1e(pLog, m_pAtoms, Rdm, *m_pOrbBasis, *m_pOrbBasis, pMpt1e, Mem);
      NewProperties.push_back(pPropertyResult);
      if (pLog)
         pLog->WriteLine();
   }
   // could go through response properties and other stuff next...

   // replace previous set of tabulated atom set properties (if present) by what
   // we got now.
   FAtomSet
      &AtomsOut = *const_cast<FAtomSet*>(&*m_pAtoms);
   AtomsOut.SwapPropertyList(NewProperties);
   IR_SUPPRESS_UNUSED_WARNING(nOcc);
}



} // namespace ct
