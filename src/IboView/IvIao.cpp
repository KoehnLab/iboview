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

#include "format.h"
#include <iostream>
#include <algorithm>
#include <set>

#include "CxAlgebra.h"
#include "CxTiming.h"
#include "CxPodArray.h"
#include "CtMatrix.h"
#include "IvIao.h"
#include "IvDocument.h"
// #include "CxRandom.h" // only needed when random rotation code is included.
#include "IvOrbital.h"

#include "CtIao.h"
#include "CtOrbLoc.h"


using namespace ct;



// predicate: sorts orbitals first by type, then by occupation, and finally
// by energy.
struct FSortOrbitalPred
{
   bool operator() (FDataSetPtr pA, FDataSetPtr pB) {
      FOrbital
         *pOrbA = dynamic_cast<FOrbital*>(pA.get()),
         *pOrbB = dynamic_cast<FOrbital*>(pB.get());
      if (pOrbA == 0 && pOrbB == 0)
         return false;
      if (pOrbA == 0)
         return true;
      if (pOrbB == 0)
         return false;
      // so... both are orbitals.
      if (pOrbA->info.Spin < pOrbB->info.Spin)
         return true;
      if (pOrbB->info.Spin < pOrbA->info.Spin)
         return false;
//       if (pOrbA->info.fOcc < pOrbB->info.fOcc)
//          return true;
//       if (pOrbB->info.fOcc < pOrbA->info.fOcc)
//          return false;
      return pOrbA->info.fEnergy < pOrbB->info.fEnergy;
   }
};



static FMatrixView MakeValenceVirtuals(FMatrixView CIb, FMatrixView CVir, FBasisSet *pBasis, ct::FMemoryStack &Mem)
{
   // TODO: Wait... why does this have *virtual* orbitals as input?!
   // shouldn't this just compute the SVD of mdot(CIb.T, S1, COcc.T) and then
   // take the dot(CIb,U)-orbitals with zero singular values (last CIb.nCols -
   // COcc.nCols)?! Or is this something about computing orbital energies?
   // But this routine doesn't even actually do that...
   size_t
      nBf = pBasis->nFn(),
      nVir = CVir.nCols,
      nIb = CIb.nCols;
   assert(CIb.nRows == nBf && CVir.nRows == nBf);
   FMatrixView
      CValVir = MakeStackMatrix(nBf, nIb, Mem);
   // ^- WARNING: must be before the first FStackMatrix! That's an output quantity!
   {
      FStackMatrix
         SIbVir(nIb, nVir, &Mem);

      {
         FStackMatrix
            S(nBf, nBf, &Mem);
         MakeIntMatrix(S, *pBasis, *pBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
         ChainMxm(SIbVir, Transpose(CIb), S, CVir, Mem);
      }

      FMatrixView
         U, Vt;
      double
         *pSigma;
      AllocAndComputeSvd(U, pSigma, Vt, SIbVir, Mem);

      size_t
         nSigNonzero = 0;
      while (nSigNonzero < U.nCols && std::abs(pSigma[nSigNonzero]) > 1e-8)
         nSigNonzero += 1;
      U.nCols = nSigNonzero;

   //    CValVir = MakeStackMatrix(nBf, U.nCols, Mem);
      CValVir.nCols = U.nCols;
      Mxm(CValVir, CIb, U);
   }
   // todo: maybe semi-canonicalize? Could be useful for some transformations.
   return CValVir;
}


static void MakeOrbitalEnergies_UnitaryTrafoNxN(FMatrixView &NewOrbEw, FMatrixView CIbOcc, FMatrixView PrevCIbOcc, TArray<FOrbital *> const &RefOrbitals, FMemoryStack &Mem)
{
   if ((CIbOcc.nRows != PrevCIbOcc.nRows) || (CIbOcc.nCols > PrevCIbOcc.nCols))
      throw std::runtime_error("internal error in MakeUnitaryTrafoOrbitalEnergies: supplied orbital coefficient matrices pre and post transform differ in shape.");
   assert_rt(NewOrbEw.nRows == CIbOcc.nCols && NewOrbEw.nCols == 1);
   // collect diagonal fock matrix elements of reference orbitals.
   // If the reference orbitals were canonical, this will allow us to compute
   // orbital energies for the new, transformed orbitals.
   size_t
      nNewOrb = CIbOcc.nCols,
      nPrevOrb = PrevCIbOcc.nCols,
      nRefOrb = RefOrbitals.size();
   assert_rt(nRefOrb == nPrevOrb);
   FStackMatrix
      RefOrbEw(nRefOrb, 1, &Mem);
   for (size_t iRefOrb = 0; iRefOrb < nRefOrb; ++ iRefOrb)
      RefOrbEw[iRefOrb] = RefOrbitals[iRefOrb]->info.fEnergy;

   // compute new orbital energies---assuming old orbitals were canonical.
   FStackMatrix
      U(nPrevOrb, nNewOrb, &Mem);
   Mxm(U, Transpose(PrevCIbOcc), CIbOcc);
   assert_rt(RefOrbEw.nRows == NewOrbEw.nRows);
   for (size_t iNewOrb = 0; iNewOrb < nNewOrb; ++ iNewOrb) {
      NewOrbEw[iNewOrb] = 0;
      for (size_t iPrevOrb = 0; iPrevOrb < nPrevOrb; ++ iPrevOrb)
         NewOrbEw[iNewOrb] += RefOrbEw[iPrevOrb] * sqr(U(iPrevOrb, iNewOrb));
   }
}


static void MakeOrbitalEnergies_General(FMatrixView &NewOrbEw, FMatrixView CNew, TArray<FOrbital*> const &RefOrbitals, FBasisSet const *pBasis, FMemoryStack &Mem)
{
   // we get here when making valence virtual orbitals.
   // note: theoretically only the part of the ref orbitals matters which lies
   // *inside* the IAO basis. So this *can* be represented in a nOcc*nOcc Fock matrix.
   // The issues are that
   // (a) in RefOrbEw, the energies collected there omit the
   //     transformation from CVir -> CValVir (here: in COcc). That is the main issue.
   // (b) RefOrbEw only collects the diagonal Fock matrix elements---which is okay
   //     if the reference orbitals are canonical. But latest after the CValVir
   //     transform they no longer are!
   // ...some cleanup might anyway be good. Probably I should just make
   // the stupid IAO basis Fock matrix, and handle both cases with the same code...

   assert_rt(NewOrbEw.nRows == CNew.nCols && NewOrbEw.nCols == 1);
   size_t
      nBf = pBasis->nFn();
   // make a matrix for the coefficients of the reference orbitals (CNew
   // should lie in the span of those ones, and uniquely.)
   size_t
      nRefOrb = RefOrbitals.size();
   FStackMatrix
      CRefOrb(nBf, nRefOrb, &Mem);
   for (size_t iOrb = 0; iOrb != RefOrbitals.size(); ++ iOrb) {
      FOrbital *pOrb = RefOrbitals[iOrb];
      if (pOrb->pBasisSet != pBasis)
         throw std::runtime_error("internal error in MakeValenceVirtualOrbitalEnergies: supplied main basis and reference orbital basis not identical.");
      for (size_t iFn = 0; iFn != nBf; ++ iFn)
         CRefOrb(iFn, iOrb) = pOrb->pCoeffs[iFn];
   }

   FStackMatrix
      S(nBf, nBf, &Mem);
   MakeIntMatrix(S, *pBasis, *pBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);

   FStackMatrix
      SMo(CRefOrb.nCols, CNew.nCols, &Mem);
   ChainMxm(SMo, Transpose(CRefOrb), S, CNew, Mem);

   for (size_t iNewOrb = 0; iNewOrb != CNew.nCols; ++ iNewOrb) {
      NewOrbEw[iNewOrb] = 0.;
      for (size_t iRefOrb = 0; iRefOrb != nRefOrb; ++ iRefOrb)
         NewOrbEw[iNewOrb] += RefOrbitals[iRefOrb]->info.fEnergy * sqr(SMo(iRefOrb, iNewOrb));
   }
}


void FFrame::RunIaoAnalysis(ct::FLog &RealLog, FWfOptions *pWfOptions, bool AllowLocalize, ct::FMemoryStack &Mem)
{
   // what does this do? Normally the RealLog object is coupled to the ShowText form's
   // GUI widget, and writing something to it (with flush) will cause it to re-draw.
   // Now, there are *many* individual write events inside the analysis... and this can become
   // really slow for large molecules. Not because of the chemical analysis, but because of the O(N^2)
   // scaling of the GUI widget update. So we here just buffer the entire thing into a memory stream
   // and then just write the final output to the output log once done...
   std::stringstream
      str;
   FLogStdStream
      DummyLog(str);
   RunIaoAnalysisImpl(DummyLog, pWfOptions, AllowLocalize, Mem);
// #if __GNUC__ <= 6
   RealLog.w.write("{}" + str.str(), "");
// #else
//    RealLog.w.write(str.str());
//    // ^- causes error on some very old g++? Reported by Mark Vilensky.
//    // E.g., the gcc 4.8 coming with 2019' OpenSUSE...
//    // UPDATE: the compiler version check doesn't work, either.
// #endif
   RealLog.Flush();
}


FBasisSetPtr MakeMinBasis(FAtomSet const &Atoms, FBasisSetPtr pOrbBasis)
{
   // and return make a new minimal basis by instanciating the BASIS_Guess context,
   // unless any of the atoms has the "USE-ORB-BASIS" flag set, in which case
   // a reference to the main basis is returned.
   //
   // This is a hack to enable reading and analysis of semi-empirical program
   // output, such as from GFN2-XTB. However, it is not quite right in all cases,
   // because these programs may still have polarization functions...
   for (size_t iAt = 0; iAt != Atoms.size(); ++ iAt) {
      if (Atoms[iAt].BasisDesc.find(BASIS_Guess)->second == "USE-ORB-BASIS") {
         return pOrbBasis;
      }
   }
   return FBasisSetPtr(new FBasisSet(Atoms, BASIS_Guess));
}




void FFrame::RunIaoAnalysisImpl(FLog &Log, FWfOptions *pWfOptions, bool AllowLocalize, ct::FMemoryStack &Mem)
{
   if (!HaveOrbitals())
      // nothing to make coefficients of...
      return;
   bool
      InputOrbitalsChanged = false;

   Log.WriteProgramIntro("INTRINSIC BASIS BONDING ANALYSIS (IBBA)");
   FTimer
      tIbbaTotal;

   FWfType
      WfType = GetWfType();
   if (WfType != wfi::WFTYPE_Rhf && WfType != wfi::WFTYPE_Uhf) {
      Log.Write("Sorry, wave function type not supported. Can only do RHF/RKS and UHF/UKS at this moment.");
      return;
   }

   {
      QString
         OrbDivision = pWfOptions->GetOrbDivision().toLower();
      if (OrbDivision == "closed/open" || OrbDivision == "as input wf") {
         // ^-- "Closed/Open" was the option used before 2021. For latter
         // analyses in terms of forced closed/open corresponding orbitals,
         // we may introduce a "Force closed/open" option.
         // ...in any case: this is the default behavior, so there is nothing
         // to do.
      } else if (OrbDivision == "force alpha/beta" || OrbDivision == "alpha/beta") {
         // we are asked to localize the alpha and beta orbital manifolds separately.
         // To this end we simply convert the input wave function from RHF form to UHF
         // form by splitting all spin-free orbitals into separate alpha- and beta-spin
         // components.
         if (WfType != wfi::WFTYPE_Uhf) {
            ConvertToUnrestrictedOrbitals();
            InputOrbitalsChanged = true;
         }
         // FIXME: valence-virtual for beta-spin does not work in this case!
         WfType = GetWfType(); // re-compute wf type. Should be UHF now.
         assert(WfType == wfi::WFTYPE_Uhf);
      } else {
         IvNotify(NOTIFY_Warning, IvFmt("Orbital division type '%1' not recognized. Employing 'as input wf' instead.", pWfOptions->GetOrbDivision()));
      }
   }

   using namespace ct;
   FAtomSet
      *pAtoms = pGetAtoms();

   TMemoryLock<char>
      pBaseOfMemory(0, &Mem);

   FBasisSetPtr
      pBasis,
      pMinBasis;
//       pMinBasis = new FBasisSet(*pAtoms, BASIS_Guess);
   FMatrixView
      COccAC, // alpha and closed
      COccBC; // beta and closed
   {
      TArray<FOrbital*>
         RefOrbitalsAC,
         RefOrbitalsBC;
      MakeOrbitalMatrix(COccAC, pBasis, &RefOrbitalsAC, ORBMAT_OccupiedOnly | ORBMAT_AlphaAndClosedOnly, Mem);
      assert(RefOrbitalsAC.size() == COccAC.nCols);
      if (WfType == wfi::WFTYPE_Uhf) {
         MakeOrbitalMatrix(COccBC, pBasis, &RefOrbitalsBC, ORBMAT_OccupiedOnly | ORBMAT_BetaAndClosedOnly, Mem);
         assert(RefOrbitalsBC.size() == COccBC.nCols);
      } else {
         COccBC = COccAC;
      }

   }
   if (pBasis.get() == 0 || (COccAC.nCols == 0 && COccBC.nCols == 0)) {
      Log.Write("Sorry, this frame has NO occupied orbitals. Cannot do IAO analysis.");
      return;
   }

   pMinBasis = MakeMinBasis(*pAtoms, pBasis);

   uint
      nBf = pBasis->nFn(),
      nIb = pMinBasis->nFn(),
      nOccAC = COccAC.nCols,
      nOccBC = COccBC.nCols;

   if (nOccAC > nIb || nOccBC > nIb) {
      Log.Write(" More occupied orbitals than valence orbitals. IAO analysis aborted.");
      return;
   }

   // form orthogonal IAO vectors.
   uint
      IaoFlags = IAO_NormalizeInput;
   char const
      *pOrthMethodName = "Unknown";
   if (pWfOptions->GetAoType() == "IAO (Sym Orth.)") {
      IaoFlags |= IAO_OrthSym;
      pOrthMethodName = "Symmetric (Loewdin)";
   } else if (pWfOptions->GetAoType() == "IAO (ZBD Orth.)") {
      IaoFlags |= IAO_OrthZbd;
      pOrthMethodName = "Zero-Bond-Dipole (Laikov)";
   } else {
      IaoFlags |= IAO_OrthSym;
      IvNotify(NOTIFY_Error, IvFmt("Intrinsic basis type '%1' not recognized. Using IAO (Sym Orth.) instead.", pWfOptions->GetAoType()));
      pOrthMethodName = "Symmetric (Loewdin)";
   }

   bool
      MakeAntibonds = (pWfOptions->GetOrbDisplay() != "Occupied only");


   Log.Write(" {:<34}{}", "Polarized AO type:", "IAO/2014");
   Log.Write(" {:<34}{}", "Orthogonalization method:", pOrthMethodName);
   Log.Write(" {:<34}{}", "Output orbitals:", q2s(pWfOptions->GetOrbDisplay()));

   FMatrixView
      // intrinsic basis spanning the closed & alpha orbitals. In RHF, this includes the span of the beta orbitals.
      CIbAC = MakeStackMatrix(nBf, nIb, Mem),
      // intrinsic basis spanning closed & beta orbitals
      // (used for formal almost-UHF wave function employing some closed-shell orbs;
      // e.g., constructed as corresponding orbitals between alpha- and beta manifolds)
      CIbBC;
   if (WfType == wfi::WFTYPE_Uhf) {
      CIbBC = MakeStackMatrix(nBf, nIb, Mem);
   } else {
      CIbBC = CIbAC;
   }

   {
      FStackMatrix
         S(nBf, nBf, &Mem),
         Scd(nBf, nBf, &Mem);
      MakeIntMatrix(S, *pBasis, *pBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
//       S.Print(std::cout, "MAIN BASIS OVERLAP MATRIX #24");
      if (1) {
         // apply a sanity check to the orbitals.
         FStackMatrix
            SMo(COccAC.nCols, COccAC.nCols, &Mem);
         ChainMxm(SMo, Transpose(COccAC), S, COccAC, Mem);
         double
            fRmsd = SMo.fRmsdFromIdentity();
         if (fRmsd > 1e-5) {
            fmt::MemoryWriter ss;
            ss.write("Sanity check for input orbitals failed! Rmsd from orthogonality: {:8.2e}", fRmsd);
            Log.EmitWarning(ss.str());
            Log.Write(" If the input orbitals came from a .molden file: Beware of differnt normalizations for"
                    "\n different programs. Options: (1) Check for a newer version of IboView, (2) If it"
                    "\n still happens with the current version, contact Gerald Knizia.");
            IvNotify(NOTIFY_Warning, s2q(fmt::format("Input orbitals are not orthogonal (deviation: {:.2e}). Check frame log for details.", fRmsd)));
         }
//          SMo.Print(std::cout, "MO BASIS OVERLAP MATRIX #22");
      }
      Assign(Scd, S);
      CalcCholeskyFactors(Scd);
      MakeIaoBasisNew(CIbAC, &*pMinBasis, *pAtoms, &*pBasis, COccAC, S, Scd, IaoFlags, Mem);
      if (CIbBC.pData != CIbAC.pData) {
         MakeIaoBasisNew(CIbBC, &*pMinBasis, *pAtoms, &*pBasis, COccBC, S, Scd, IaoFlags, Mem);
      }
//       TestOverlap("IB[a] vs IB[b] (after construction)", Log, CIbAC, CIbBC, pBasis.get(), Mem);
   }


   if (0) {
      FStackMatrix
         Rdm(nBf, nBf, &Mem),
         S(nBf, nBf, &Mem),
         Op(nBf, nBf, &Mem);
      Mxm(Rdm, COccAC, Transpose(COccAC));
      Op.Clear();
      pAtoms->AddKineticMatrix(Op, *pBasis, *pBasis, Mem);
      MakeIntMatrix(S, *pBasis, *pBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
//       MakeIntMatrix(Op, *pBasis, *pBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
      Log.WriteResult("Nuclear repulsion energy", pAtoms->NuclearRepulsionEnergy());
      Log.WriteResult("Kinetic energy", 2*Dot(Rdm,Op));
      Log.WriteResult("Number of electrons", Dot(Rdm,S));
   }


   // setup the localization options.
   FLocalizeOptions
      LocOpt;
   QString
      LocMethod = pWfOptions->GetLocMethod().toLower();
   if (AllowLocalize) {
      LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_pow_2;
      if (LocMethod == "ibo (exponent 2)") {
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_pow_2;
      } else if (LocMethod == "ibo (exponent 4)") {
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_pow_4;
      } else if (LocMethod == "ibo (min entropy)" || LocMethod == "ibo (entropy)") {
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_Entropy;
      } else if (LocMethod == "ibo (exponent 3/2)" || LocMethod == "ibo (exponent 1.5)") {
//          LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_pow_15;
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_gen_h;
         LocOpt.pScalarFn01 = MakeLocFn1_nPow(1.5, 1., 0.);
      } else if (LocMethod == "ibo (log gamma)" || LocMethod == "ibo (log gamma[n+2])") {
//          LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_Log2Gamma_nPlus2;
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_gen_h;
         LocOpt.pScalarFn01 = MakeLocFn1_LogBeta(2.);
      } else if (LocMethod == "ibo (log gamma)" || LocMethod == "ibo (log gamma[n+1])") {
//          LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_Log2Gamma_nPlus1;
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_gen_h;
         LocOpt.pScalarFn01 = MakeLocFn1_LogBeta(1.);
      } else if (LocMethod == "none") {
         LocOpt.LocFn = FLocalizeOptions::LOCFN_None;
      } else {
         LocOpt.LocFn = FLocalizeOptions::LOCFN_nA_pow_2;
         Log.EmitWarning(q2s(IvFmt("Localization type '%1' not recognized. Defaulting to IBO (exponent 2).", LocMethod)));
      }
      Log.Write(" {:<34}{}", "Localization functional:", LocOpt.MakeLocalizeFnDesc());
      Log.WriteLine();
   } else {
      LocOpt.LocFn = FLocalizeOptions::LOCFN_None;
      LocMethod = "none";
   }

   LocOpt.ThrLoc = std::min(double(1e-6), double(1e-2 * pWfOptions->GetThrGrad()));
   // ^- WfOptions::ThrGrad is actually the SCF gradient threshold.
   // ^- cgk 2019-10-22: why was this such a loose threshold? It is not even sufficient
   //    to find a fully reliable localization in acrylic acid if the random orbital
   //    transformation is turned on! The issue is mixing core orbitals with lone pair
   //    orbitals, and the corresponding gradient being very small for the localization
   //    functional itself.
   LocOpt.ThrLoc = 1e-12;

   typedef std::set<FOrbital*>
      FOrbitalRefSet;
   FOrbitalRefSet
      ProcessedOrbitals;

   int
      iLastCase = 1; // 0 and 1: closed/alpha (RHF) or alpha/beta (UHF)
   if (MakeAntibonds && (WfType == wfi::WFTYPE_Rhf))
      iLastCase = 2; // 0,1,2: closed/active/external
   if (MakeAntibonds && (WfType == wfi::WFTYPE_Uhf))
      iLastCase = 3; // 0,1,2,3: alpha occ // beta occ // alpha vir // beta vir
   for (int iCase = 0; iCase <= iLastCase; ++ iCase)
   {
      TMemoryLock<double>
         pFreeMeInLoop(0, &Mem);
      // need to define:
      // - RefOrbitals set and COcc
      // - Coeffs of IB vectors which span this space
      // - factors for spin and total densities?

      QString
         SpaceDesc = "";
      TArray<FOrbital*>
         RefOrbitals; // orbitals in current subspace which are localized together (i.e., freely mixable w/o changing the wave function)
      FMatrixView
         CLoc; // corresponding orbital matrix (orbitals to localize)
      double
         fChargeFactor,
         fSpinFactor;
      FMatrixView
         CIb; // IB which spanns the orbitals of the current subspace.
      bool
         // if set, assume the RefOrbitals lie in the span of CIb; in this case,
         // a simpler and faster formula can be used to compute the new orbital
         // energies (assuming RefOrbitals was canonical)
         OrbitalEw_RefOrbitalsInSpanOfIb = true;
      if (WfType == wfi::WFTYPE_Rhf) {
         CIb = CIbAC; // <- that's all the occupied orbtials in RHF.
         if (iCase == 0) {
            // closed-shell orbitals.
            SpaceDesc = "AB";
            fChargeFactor = 2.;
            fSpinFactor = 0.;
            MakeOrbitalMatrix(CLoc, pBasis, &RefOrbitals, ORBMAT_OccupiedOnly | ORBMAT_ClosedOnly, Mem);
         } else if (iCase == 1) {
            // open-shell orbitals (alpha occupied).
            SpaceDesc = "A_";
            fChargeFactor = 1.;
            fSpinFactor = 1.;
            MakeOrbitalMatrix(CLoc, pBasis, &RefOrbitals, ORBMAT_OccupiedOnly | ORBMAT_AlphaOnly, Mem);
         } else if (iCase == 2) {
            // external orbitals (unoccupied)
            SpaceDesc = "__";
            fChargeFactor = 1.;
            fSpinFactor = 1.;
            FMatrixView CVir;
            MakeOrbitalMatrix(CVir, pBasis, &RefOrbitals, ORBMAT_VirtualOnly, Mem);
            // ^- TODO: if I actually had all the virtual orbitals, would the LeastSquaresSolve be sufficient
            //    to get the span of the anti-bonds? (w/o explicitly constructing them).
            //    Would have to update the nOcc number, but if it is ONLY that...
            //    Hm... no, don't think so. Need to project. Number of orbtials will change!
            if (pBasis == 0 || CVir.nCols == 0)
               continue;
            CLoc = MakeValenceVirtuals(CIb, CVir, pBasis.get(), Mem);
            OrbitalEw_RefOrbitalsInSpanOfIb = false;
         } else {
            continue;
         }
      } else if (WfType == wfi::WFTYPE_Uhf) {
         if (iCase == 0) {
            // alpha orbitals (occupied)
            CIb = CIbAC;
            SpaceDesc = "A_";
            fChargeFactor = 1.;
            fSpinFactor = 1.;
            MakeOrbitalMatrix(CLoc, pBasis, &RefOrbitals, ORBMAT_OccupiedOnly | ORBMAT_AlphaOnly, Mem);
         } if (iCase == 1) {
            // beta orbitals (occupied)
            CIb = CIbBC;
            SpaceDesc = "_B";
            fChargeFactor = 1.;
            fSpinFactor = -1.;
            MakeOrbitalMatrix(CLoc, pBasis, &RefOrbitals, ORBMAT_OccupiedOnly | ORBMAT_BetaOnly, Mem);
         } if (iCase == 2) {
            // alpha orbitals (virtual)
            CIb = CIbAC;
            SpaceDesc = "x_";
            fChargeFactor = 1.;
            fSpinFactor = 1.;
            FMatrixView CVir;
            MakeOrbitalMatrix(CVir, pBasis, &RefOrbitals, ORBMAT_VirtualOnly | ORBMAT_AlphaOnly, Mem);
            if (pBasis == 0 || CVir.nCols == 0)
               continue;
            CLoc = MakeValenceVirtuals(CIb, CVir, pBasis.get(), Mem);
            OrbitalEw_RefOrbitalsInSpanOfIb = false;
         } if (iCase == 3) {
            // beta orbitals (virtual)
            CIb = CIbBC;
            SpaceDesc = "_x";
            fChargeFactor = 1.;
            fSpinFactor = -1.;
            FMatrixView CVir;
            MakeOrbitalMatrix(CVir, pBasis, &RefOrbitals, ORBMAT_VirtualOnly | ORBMAT_BetaOnly, Mem);
            if (pBasis == 0 || CVir.nCols == 0)
               continue;
            CLoc = MakeValenceVirtuals(CIb, CVir, pBasis.get(), Mem);
            OrbitalEw_RefOrbitalsInSpanOfIb = false;
         }
      } else {
         Log.EmitError("Unrecognized wave function type in RunIaoAnalysis (inner).");
         continue;
      }
//       Log.Write("   class: '{}'  nLoc = {}  nIb = {}", q2s(SpaceDesc), CLoc.nCols, CIb.nCols);
      IR_SUPPRESS_UNUSED_WARNING2(fSpinFactor, fChargeFactor);
      unsigned
         // number of orbitals to localize
         nLoc = CLoc.nCols;
      if (nLoc == 0)
         continue; // no orbitals in this class.

      // re-express occupied vectors in terms of the IAO vectors. This is an exact
      // transformation.
      // (note: this probably can be done in a different way than via SVD... probably via
      //      Sib^{-1} CIb.T S1 COcc, where Sib = CIb.T S1 CIb
      // I just don't want to think about it for the moment)
      // Note: BEWARE OF THE ABSORBED OCCUPATION NUMBERS!!

      // solve CIb CIbOcc = COcc for CIbOcc
      FStackMatrix
         CIbLoc(nIb, nLoc, &Mem);
      LeastSquaresSolveSafe(CIb, CIbLoc, CLoc, 1e-9, Mem);

      if (0) {
         // test if IB works.
         uint
            nBf = pBasis->nFn();
         FStackMatrix
            // COcc reconstructed from IB and CIbOcc.
            CLocRec(nBf, nIb, &Mem);
         Mxm(CLocRec, CIb, CIbLoc);
         Add(CLocRec, CLoc, -1.);
         std::cout << fmt::format("Norm(CLocRec - CLoc) = {:8.2e}   [nBf={}, nLoc={}]\n", CLocRec.fRmsdFromZero(), nBf, CLoc.nCols);
      }


      // list of new orbital energies
      FStackMatrix
         NewOrbEw(nLoc, 1, &Mem);
      NewOrbEw.Clear();


      if (LocOpt.LocFn != FLocalizeOptions::LOCFN_None) {
//          std::stringstream
//             ss;
         LocOpt.Verbosity = 1;
//          LocOpt.pxout = &ss;
         LocOpt.pLog = &Log;
         size_t
            nAt = pAtoms->size();
         TMemoryLock<size_t>
            pAtOffsets(nAt+1, &Mem);
            // ^- that is NOT the same as pRawBasis->CenterOffsets! The latter is in terms of shells!!!
         for (size_t iAt = 0; iAt <= nAt; ++ iAt)
            pAtOffsets[iAt] = pMinBasis->pRawBasis->iCenFn(iAt);
         if (nAt + 1 != pMinBasis->pRawBasis->CenterOffsets.size())
            Log.EmitError("Number of atoms in minimal basis inconsistent with frame main basis (?).");
         if (CIbLoc.nRows != pMinBasis->pRawBasis->nFn())
            Log.EmitError("Number of basis functions in minimal basis inconsistent with frame main basis (?).");

         if (0) {
            std::stringstream str;
            pMinBasis->Print(str);
            Log.Write(str.str());
            Log.Write("\ncenter offsets:");
            for (size_t i = 0; i < nAt+1; ++ i)
               Log.Write("  {:4}. {:6}", i, pAtOffsets[i]);
         }

         FStackMatrix
            PrevCIbLoc(nIb, nLoc, &Mem);
         Assign(PrevCIbLoc, CIbLoc);

         if (iCase == 0)
            Log.WriteLine();
         LocOpt.TaskDesc = "Class: " + q2s(SpaceDesc);
         LocalizeVectors(CIbLoc, pAtOffsets, nAt, LocOpt, Mem);
//          Log.Write("{}; Class: {}", ss.str(), q2s(SpaceDesc));

         // make new set of vectors by transforming CIbOcc back to the
         // main basis.
         Mxm(CLoc, CIb, CIbLoc);

         if (OrbitalEw_RefOrbitalsInSpanOfIb) {
            MakeOrbitalEnergies_UnitaryTrafoNxN(NewOrbEw, CIbLoc, PrevCIbLoc, RefOrbitals, Mem);
         } else {
            MakeOrbitalEnergies_General(NewOrbEw, CLoc, RefOrbitals, pBasis.get(), Mem);
         }
         InputOrbitalsChanged = true;
      } else {
         // I think in this case the IAO coefficients of the valence virtuals will probably
         // be garbage, if the valence-virtual case is enabled. Should probably do something
         // about that. Marker: @CaseValenceVirOnButNoLocalization
      }


      // store IAO coefficients of orbitals in the orbital data sets,
      // and also a link to the minimal basis (to distinguish what is what in the basis)
      for (uint iLoc = 0; iLoc < CLoc.nCols; ++ iLoc) {
         FOrbital *pOrb = RefOrbitals[iLoc];
         if (ProcessedOrbitals.find(pOrb) != ProcessedOrbitals.end())
            IvNotify(NOTIFY_Warning, "An orbital occured in more than one localization group!");

         ProcessedOrbitals.insert(pOrb);
         if (InputOrbitalsChanged || OrbitalEw_RefOrbitalsInSpanOfIb) {
            // FIXME: this is a hack about the problem described under @CaseValenceVirOnButNoLocalization
            pOrb->pMinBasis = pMinBasis;
            pOrb->pIaoCoeffs.clear();
            pOrb->pIaoCoeffs.insert(pOrb->pIaoCoeffs.begin(), &CIbLoc(0,iLoc), &CIbLoc(0,iLoc)+nIb);
         }
         if (InputOrbitalsChanged) {
            RefOrbitals[iLoc]->info.fEnergy = NewOrbEw[iLoc];
            pOrb->pCoeffs.clear();
            pOrb->pCoeffs.insert(pOrb->pCoeffs.begin(), &CLoc(0,iLoc), &CLoc(0,iLoc)+nBf);
         }
      }
   }
   if (InputOrbitalsChanged) {
      // remove all orbitals which were not processed. This is used to delete either
      // all or non-valence virtual orbitals after the localization procedure.
      if (1) {
         FDataSetList
            NewData;
         for (size_t i = 0; i < m_Data.size(); ++ i) {
            FOrbital
               *pOrb = dynamic_cast<FOrbital*>(m_Data[i].get());
            if (pOrb == 0 || ProcessedOrbitals.find(pOrb) != ProcessedOrbitals.end())
               // copy pointers to all data sets which either are not orbitals,
               // or are orbitals which were processed.
               NewData.push_back(m_Data[i]);
         }
         m_Data.swap(NewData);
      }

      // sort new orbitals (by type and energy) and update their descriptions
      std::stable_sort(m_Data.begin(), m_Data.end(), FSortOrbitalPred());

      size_t iOrb = 0;
      for (size_t i = 0; i < m_Data.size(); ++ i) {
         FOrbital *pOrb = dynamic_cast<FOrbital*>(m_Data[i].get());
         if (pOrb) {
            pOrb->info.iOriginalIndex = i;
            pOrb->UpdateDescFromInfo(iOrb);
//             Log.Write("   pOrb[{:4}]  E = {:12.5f}  NewDesc = '{}'", i, pOrb->info.fEnergy, q2s(pOrb->GetDesc()));
            iOrb += 1;
         }
      }

      // print orbital composition in terms of atomic partial charges (with new order)
      double ThrPrint = 0.02;
      iOrb = 0;
      Log.Write("\n Summary of localized orbital composition:   [THRPRINT={:6.3f}]", ThrPrint);
      Log.Write("\n   ORB  GRP   ORB.ENERGY   CENTERS/CHARGES");
      for (size_t i = 0; i < m_Data.size(); ++ i) {
         FOrbital *pOrb = dynamic_cast<FOrbital*>(m_Data[i].get());
         if (pOrb) {
            Log.Write("{:6}  {:>3} {:12.6f} {}", iOrb+1, pOrbDescFromType(pOrb->info.Spin, pOrb->info.fOcc),
                  pOrb->info.fEnergy, q2s(pOrb->MakeFullDesc(ThrPrint, FOrbital::ORBDESC_ChargesOnly)));
            iOrb += 1;
         }
      }
   }

   if (1) {
      Log.WriteLine();
      RunBondOrderAnalysis(Log, 0, &Mem);
   }

   if (1) {
      Log.WriteLine();
      RunChargeAnalysis(Log, 0);
   }

   Log.WriteLine();
   Log.WriteTiming("IBBA", (double)tIbbaTotal);
   Log.WriteLine();

//    GetBondOrderType
//    GetThrGrad
//    GetOrbDisplay
//    GetOrbDivision
//    GetLocMethod
}




// ct::FMatrixView FFrame::MakeIaoBasis(ct::FAtomSet *pAtoms, ct::FMemoryStack &Mem)
// {
//    IR_SUPPRESS_UNUSED_WARNING2(pAtoms, Mem);
//    return m_CIb;
// }


void FDocument::MakeHybridsForSelectedAtomGroup()
{
}
