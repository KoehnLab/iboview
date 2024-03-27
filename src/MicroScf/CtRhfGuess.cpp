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
#include <stdexcept> // for std::runtime_error
#include <set> // for keeping track of which atoms we already emitted guess info for.
#include "Ir.h"
#include "CtRhfGuess.h"


namespace ct {

extern float g_AtomicOccupations[1773];
extern short g_AtomicOccupationRanges[118][2];


// some parameters for making the correct density guesses. In cases of ECPs, we need to skip the inner shells
// for which there are no actual basis functions (or electrons) in the calculation. These are the most common ECPs.
// might need updates for some of the more funky ones (e.g., there are ECP85-92 (all of them!) for various elements).
// But I guess anyone brave enough to employ these would not be doing so using DFT with MicroScf...
static size_t const g_EcpMaxShellAm = 5;
static size_t const g_NoEcp_AbsorbedShells[g_EcpMaxShellAm] = {0,0,0,0,0}; // all electron -- nothing absorbed.
static size_t const g_Ecp2_AbsorbedShells[g_EcpMaxShellAm] = {1,0,0,0,0}; // 1s absorbed by ECP (He core)
static size_t const g_Ecp10_AbsorbedShells[g_EcpMaxShellAm] = {2,1,0,0,0}; // 1s 2s 2p absorbed by ECP (Neon core)
static size_t const g_Ecp28_AbsorbedShells[g_EcpMaxShellAm] = {3,2,1,0,0}; // 1s 2s 2p 3s 3p 3d absorbed by ECP (Argon core + next d shells)
static size_t const g_Ecp46_AbsorbedShells[g_EcpMaxShellAm] = {4,3,2,0,0}; // 1s 2s 2p 3s 3p 3d 4s 4p 4d absorbed by ECP
static size_t const g_Ecp60_AbsorbedShells[g_EcpMaxShellAm] = {4,3,2,1,0}; // 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f absorbed by ECP
static size_t const g_Ecp78_AbsorbedShells[g_EcpMaxShellAm] = {5,4,3,1,0}; // 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p absorbed by ECP

size_t const *pEcpAbsorbedShells(size_t nEcpElec) {
   if (nEcpElec == 0) return &g_NoEcp_AbsorbedShells[0];
   if (nEcpElec == 2) return &g_Ecp2_AbsorbedShells[0];
   if (nEcpElec == 10) return &g_Ecp10_AbsorbedShells[0];
   if (nEcpElec == 28) return &g_Ecp28_AbsorbedShells[0];
   if (nEcpElec == 46) return &g_Ecp46_AbsorbedShells[0];
   if (nEcpElec == 60) return &g_Ecp60_AbsorbedShells[0];
   if (nEcpElec == 78) return &g_Ecp78_AbsorbedShells[0];
   std::stringstream ss;
   ss << "Sorry, ECP shell configuration for ECP" << nEcpElec << " not tabulated. Cannot make initial guess.";
   throw std::runtime_error(ss.str());
}


// Make a closed-shell orbital matrix in the guess basis, containing
// the atomic occupations from free atoms.
void CollectFreeAtomOrbitals(FMatrixView &OrbGu, FMatrixView &fOccGu, FAtomSet const &Atoms, FBasisSet *pGuessBasis, FLog *pLog, FHfPrintOptions const &PrintOptions, FMemoryStack &Mem)
{
   size_t
      nGu = pGuessBasis->nFn();
   assert(OrbGu.IsSquare() && OrbGu.nRows == nGu);
   assert(fOccGu.nRows == nGu && fOccGu.nCols == 1);

   // in case we are asked to print initial orbital information, use this to
   // keep track of the elements for which we already emitted guess information
   typedef std::set<int>
      FElementSet;
   FElementSet
      ElementsAlreadyPrinted;
   IR_SUPPRESS_UNUSED_WARNING(ElementsAlreadyPrinted);

   OrbGu.Clear();
   fOccGu.Clear();

   std::size_t
      nOcc = 0,
      iOffGu = 0;
   for (size_t iAt = 0; iAt < Atoms.size(); ++ iAt) {
      // make an atomic orbital matrix and project it to the full basis.
      FAtom
         At = Atoms[iAt];
      FAtomSet
         OneAtom;
      OneAtom.AddAtom(At);

      FBasisSet
         BasisG(OneAtom, BASIS_Guess);
      size_t
         nBfG = BasisG.nFn();
      FMatrixView
         OrbG(&OrbGu(iOffGu, nOcc), nBfG, nBfG, 1, nGu);

      size_t const
         *pEcpShells = pEcpAbsorbedShells(At.nEcpElec);

      size_t
         nOccG = 0,
         nBfOff = 0,
         iBeg = g_AtomicOccupationRanges[At.iElement-1][0],
         iEnd = g_AtomicOccupationRanges[At.iElement-1][1];

      if (pLog != 0 && PrintOptions[SCFPRINT_GuessInfo]) {
         FElementSet::iterator
            itPrinted = ElementsAlreadyPrinted.find(At.iElement);
         if (itPrinted == ElementsAlreadyPrinted.end()) {
            ElementsAlreadyPrinted.insert(At.iElement);

            std::string Ind = "     ";
            pLog->Write("\n |-- Atomic orbital guess for atom {} {}", iAt+1,At.ElementName());
            FLogIoFlagOverride(pLog, FLog::IOFLAG_FlushAfterWrite, false);
            pLog->w << Ind << "Tabulated free-atom AO basis set:" << "\n\n";
            pLog->w << BasisG.GetDesc(FBasisSet::PRINT_Default, Ind + "| ") << "\n";
            pLog->Write("{}{:<31} [{},{},{},{},{}]", Ind, "pEcpAbsorbedShells:", pEcpShells[0], pEcpShells[1], pEcpShells[2], pEcpShells[3], pEcpShells[4]);
            pLog->WriteNoNl("{}{:<31}", Ind, "Neutral atom shell occupation:", iBeg, iEnd);
            int iAm = -1;
            double
               fElecTotal = 0;
            double
               fElecAbsorbed = 0;
            for (size_t iAtOcc = iBeg; iAtOcc != iEnd; ++ iAtOcc) {
               float f = g_AtomicOccupations[iAtOcc];
               if (f == 0 || iAtOcc == iBeg) {
                  iAm += 1;
                  if (iAtOcc != iBeg)
                     pLog->w << "]";
                  if (iAtOcc != iEnd-1)
                     pLog->WriteNoNl(" {}[", "spdfghi"[iAm]);
                  if (f == 0)
                     continue;
               }
               if (float(int(f)) == f)
                  pLog->WriteNoNl("{:x}", int(f));
               else
                  pLog->WriteNoNl("({:.2f})", f);
               fElecTotal += f;
            }
            if (0)
               pLog->Write("    (<- from g_AtOcc[{}:{}])", iBeg, iEnd);
            else
               pLog->WriteLine();
            for (size_t lEcp = 0; lEcp != 5; ++ lEcp) {
               fElecAbsorbed += 2*(2*lEcp + 1) * pEcpShells[lEcp];
            }
            pLog->Write("{}{:<31} {:12.6f}\n", Ind, "Neutral atom #elec:", fElecTotal);
            pLog->Write("{}{:<31} {:12.6f}\n", Ind, "ECP absorbed #elec:", fElecAbsorbed);
         }
      }


      // assemble orbital matrix with averaged occupation numbers
      // in basis of atom.
      OrbG.Clear();
      for (size_t l = 0; l < BasisG.Shells.size() && iBeg < iEnd; ++ l) {
         ct::FAtomShell const
            &Bf = *BasisG.Shells[l].pAs;
         assert_rt(Bf.AngMom == l);
         assert_rt(l <= g_EcpMaxShellAm);

         // skip shells of current AngMom absorbed by ECP.
         size_t
            nEcpShellsL = pEcpShells[l];
         iBeg += nEcpShellsL;

         size_t
            iEnd_ = iBeg,
            nSh = 2*l + 1;
         while (g_AtomicOccupations[iEnd_] != 0)
            ++ iEnd_;
         size_t
            nCo = iEnd_ - iBeg;
         assert_rt( Bf.nCo() >= nCo );
         for (size_t iCoOrb = 0; iCoOrb < nCo; ++ iCoOrb){
            double
               fOccSqrt = std::sqrt(g_AtomicOccupations[iBeg + iCoOrb]/static_cast<double>(nSh));
            for (size_t iSh = 0; iSh < nSh; ++ iSh) {
               OrbG(nBfOff + iSh + iCoOrb * nSh, nOccG + iSh + iCoOrb * nSh) = fOccSqrt;
               fOccGu[nOcc + nOccG + iSh + iCoOrb * nSh] = sqr(fOccSqrt);
            }
         }

         iBeg = iEnd_ + 1;
         nOccG += nSh * nCo;
         nBfOff += Bf.nFn();
      }
      nOcc += nOccG;
      iOffGu += nBfG;
   }
   OrbGu.nCols = nOcc;

   if (pLog != 0 && PrintOptions[SCFPRINT_GuessInfo] && nOcc != 0) {
      pLog->WriteLine();
   }

   IR_SUPPRESS_UNUSED_WARNING(pLog);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


// // Make a closed-shell orbital matrix in the guess basis, containing
// // the atomic occupations from free atoms.
// void CollectFreeAtomOrbitals(FMatrixView &OrbGu, FMatrixView &fOccGu, FAtomSet const &Atoms, FBasisSet *pGuessBasis, FLog *pLog, FHfPrintOptions const &PrintOptions, FMemoryStack &Mem)
// {
//    size_t
//       nGu = pGuessBasis->nFn();
//    assert(OrbGu.IsSquare() && OrbGu.nRows == nGu);
//    assert(fOccGu.nRows == nGu && fOccGu.nCols == 1);
//
//    // in case we are asked to print initial orbital information, use this to
//    // keep track of the elements for which we already emitted guess information
//    typedef std::set<int>
//       FElementSet;
//    FElementSet
//       ElementsAlreadyPrinted;
//    IR_SUPPRESS_UNUSED_WARNING(ElementsAlreadyPrinted);
//
//    OrbGu.Clear();
//    fOccGu.Clear();
//
//    std::size_t
//       nOcc = 0,
//       iOffGu = 0;
//    for (size_t iAt = 0; iAt < Atoms.size(); ++ iAt) {
//       // make an atomic orbital matrix and project it to the full basis.
//       FAtom
//          At = Atoms[iAt];
//       FAtomSet
//          OneAtom;
//       OneAtom.AddAtom(At);
//
//       FBasisSet
//          BasisG(OneAtom, BASIS_Guess);
//       size_t
//          nBfG = BasisG.nFn();
//       FMatrixView
//          OrbG(&OrbGu(iOffGu, nOcc), nBfG, nBfG, 1, nGu);
//
//       size_t const
//          *pEcpShells = pEcpAbsorbedShells(At.nEcpElec);
//
//       size_t
//          nOccG = 0,
//          nBfOff = 0,
//          iBeg = g_AtomicOccupationRanges[At.iElement-1][0],
//          iEnd = g_AtomicOccupationRanges[At.iElement-1][1];
//
//       if (PrintOptions[SCFPRINT_GuessInfo]) {
//          FElementSet::iterator
//             itPrinted = ElementsAlreadyPrinted.find(At.iElement);
//          if (itPrinted == ElementsAlreadyPrinted.end()) {
//             ElementsAlreadyPrinted.insert(At.iElement);
//
//             std::string Ind = "     ";
//             xout << fmt::format("\n |-- Atomic orbital guess for atom {} {}\n", iAt+1,At.ElementName());
//             xout << Ind << "Tabulated free-atom AO basis set:" << "\n\n";
//             xout << BasisG.GetDesc(FBasisSet::PRINT_Default, Ind + "| ") << std::endl;
//             xout << Ind << fmt::format("{:<31} [{},{},{},{},{}]\n", "pEcpAbsorbedShells:", pEcpShells[0], pEcpShells[1], pEcpShells[2], pEcpShells[3], pEcpShells[4]);
//             xout << Ind << fmt::format("{:<31}", "Neutral atom shell occupation:", iBeg, iEnd);
//             int iAm = -1;
//             double
//                fElecTotal = 0;
//             double
//                fElecAbsorbed = 0;
//             for (size_t iAtOcc = iBeg; iAtOcc != iEnd; ++ iAtOcc) {
//                float f = g_AtomicOccupations[iAtOcc];
//                if (f == 0 || iAtOcc == iBeg) {
//                   iAm += 1;
//                   if (iAtOcc != iBeg)
//                      xout << "]";
//                   if (iAtOcc != iEnd-1)
//                      xout << fmt::format(" {}[", "spdfghi"[iAm]);
//                   if (f == 0)
//                      continue;
//                }
//                if (float(int(f)) == f)
//                   xout << fmt::format("{:x}", int(f));
//                else
//                   xout << fmt::format("({:.2f})", f);
//                fElecTotal += f;
//             }
//             if (0)
//                xout << fmt::format("    (<- from g_AtOcc[{}:{}])\n", iBeg, iEnd);
//             else
//                xout << "\n";
//             for (size_t lEcp = 0; lEcp != 5; ++ lEcp) {
//                fElecAbsorbed += 2*(2*lEcp + 1) * pEcpShells[lEcp];
//             }
//             xout << Ind + fmt::format("{:<31} {:12.6f}\n", "Neutral atom #elec:", fElecTotal);
//             xout << Ind + fmt::format("{:<31} {:12.6f}\n", "ECP absorbed #elec:", fElecAbsorbed);
//          }
//       }
//
//
//       // assemble orbital matrix with averaged occupation numbers
//       // in basis of atom.
//       OrbG.Clear();
//       for (size_t l = 0; l < BasisG.Shells.size() && iBeg < iEnd; ++ l) {
//          ct::FAtomShell const
//             &Bf = *BasisG.Shells[l].pAs;
//          assert_rt(Bf.AngMom == l);
//          assert_rt(l <= g_EcpMaxShellAm);
//
//          // skip shells of current AngMom absorbed by ECP.
//          size_t
//             nEcpShellsL = pEcpShells[l];
//          iBeg += nEcpShellsL;
//
//          size_t
//             iEnd_ = iBeg,
//             nSh = 2*l + 1;
//          while (g_AtomicOccupations[iEnd_] != 0)
//             ++ iEnd_;
//          size_t
//             nCo = iEnd_ - iBeg;
//          assert_rt( Bf.nCo() >= nCo );
//          for (size_t iCoOrb = 0; iCoOrb < nCo; ++ iCoOrb){
//             double
//                fOccSqrt = std::sqrt(g_AtomicOccupations[iBeg + iCoOrb]/static_cast<double>(nSh));
//             for (size_t iSh = 0; iSh < nSh; ++ iSh) {
//                OrbG(nBfOff + iSh + iCoOrb * nSh, nOccG + iSh + iCoOrb * nSh) = fOccSqrt;
//                fOccGu[nOcc + nOccG + iSh + iCoOrb * nSh] = sqr(fOccSqrt);
//             }
//          }
//
//          iBeg = iEnd_ + 1;
//          nOccG += nSh * nCo;
//          nBfOff += Bf.nFn();
//       }
//       nOcc += nOccG;
//       iOffGu += nBfG;
//    }
//    OrbGu.nCols = nOcc;
//
//    if (PrintOptions[SCFPRINT_GuessInfo && nOcc != 0]) {
//       xout << "\n";
//    }
//
//    IR_SUPPRESS_UNUSED_WARNING(pLog);
//    IR_SUPPRESS_UNUSED_WARNING(Mem);
// }


void MakeDensityGuess(FMatrixView Fock, FAtomSet const &Atoms, FMatrixView Scd1, FMatrixView CoreH1,
   FBasisSet *pOrbBasis, FBasisSet *pGuessBasis, FFockComponentBuilderList &FockBuilders, FLog *pLog, FHfOptions const &Options, FMemoryStack &Mem)
{
   size_t
      nAo = pOrbBasis->nFn(),
      nGu = pGuessBasis->nFn();
   assert(Fock.nRows == nAo && Fock.nCols == nAo);

   // make a closed-shell orbital matrix in the guess basis, containing
   // the atomic occupations
   FStackMatrix
      OrbGu(nGu, nGu, &Mem),
      fOccGu(nGu, 1, &Mem);
   CollectFreeAtomOrbitals(OrbGu, fOccGu, Atoms, pGuessBasis, pLog, Options.Print, Mem);
   size_t
      nOcc = OrbGu.nCols;

   bool OrthogonalizeInitialOrbs = true;

   if (OrthogonalizeInitialOrbs) {
      // sort them by CoreH matrix elements (such that we get core orbitals first)
      FStackMatrix
         CoreH(nGu, nGu, &Mem),
         OrbGu2(nGu, nOcc, &Mem),
         fOccGu2(nGu, 1, &Mem),
         fVal(nOcc,1, &Mem);
      Atoms.MakeCoreHamiltonMatrix(CoreH, *pGuessBasis, *pGuessBasis, Mem);
      Assign(OrbGu2, OrbGu);
      Assign(fOccGu2, fOccGu);
      fVal.Clear();
      for (size_t iOcc = 0; iOcc < OrbGu2.nCols; ++ iOcc)
         for (size_t iGu = 0; iGu < nGu; ++ iGu)
            fVal[iOcc] += CoreH(iGu,iGu) * sqr(OrbGu2(iGu,iOcc));
      size_t
         *pOrd;
      Mem.Alloc(pOrd, nOcc);
      ArgSort1(pOrd, fVal.pData, 1, nOcc, false); // small ones (i.e., very negative ones) first.

      for (size_t iOcc = 0; iOcc < OrbGu2.nCols; ++ iOcc) {
         for (size_t iGu = 0; iGu < nGu; ++ iGu)
            OrbGu(iGu,iOcc) = OrbGu2(iGu,pOrd[iOcc]);
         fOccGu[iOcc] = fOccGu2[pOrd[iOcc]];
      }
   }


   if (0) {
      if (OrthogonalizeInitialOrbs) {
         FStackMatrix
            S2(nGu, nGu, &Mem);
         MakeIntMatrix(S2, *pGuessBasis, *pGuessBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
         CalcCholeskyFactors(S2);
         OrthSchmidtScd(OrbGu, S2, Mem);
         // put back the occupation numbers.
         for (size_t iOcc = 0; iOcc < nOcc; ++ iOcc)
            Scale(&OrbGu(0, iOcc), std::sqrt(fOccGu[iOcc]), size_t(nGu));
      }


      FStackMatrix
         COccO(nGu, 0, &Mem); // dummy open-shell orbital matrix.
      // make a Fock matrix in the guess basis.
      FStackMatrix
         CoreH(nGu, nGu, &Mem),
         ExchO(nGu, nGu, &Mem),
         FockGu(nGu, nGu, &Mem);
      Atoms.MakeCoreHamiltonMatrix(CoreH, *pGuessBasis, *pGuessBasis, Mem);
      BuildFock(FockGu, ExchO, &*pGuessBasis, OrbGu, COccO, 0, FockBuilders, Mem);
      Add(FockGu, CoreH);
//       if (0) {
//          FStackMatrix
//             DenC(nGu, nGu, &Mem);
//          Mxm(DenC, OrbGu, Transpose(OrbGu));
//          DenC.Print(xout, "Guess basis density matrix");
//          FockGu.Print(xout, "Guess basis Fock matrix");
//          xout << fmt::format(pResultFmt, "Energy of initial guess", (.5*Dot(DenC, FockGu) + .5*Dot(DenC, CoreH) + Atoms.NuclearRepulsionEnergy())) << std::endl;
//       }

      FStackMatrix
         SmhGu(nGu, nGu, &Mem);
      MakeIntMatrix(SmhGu, *pGuessBasis, *pGuessBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
      CalcSmhMatrix(SmhGu, Mem, FSmhOptions(1e-15,1e-15,0,0));

      // make projection matrix from guess basis to main basis:
      //   P12 = inv(S1) x S12  (for projecting contra-variant vectors from 2 to 1, e.g., orbitals)
      //   P12 = S12 x inv(S2)  (for projecting co-variant vectors from 2 to 1, e.g., integrals)
      FStackMatrix
         S12(nAo, nGu, &Mem),
         P12(nAo, nGu, &Mem);
      MakeIntMatrix(S12, *pOrbBasis, *pGuessBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
      ChainMxm(P12, S12, SmhGu, SmhGu, Mem);

      // project Fock matrix from the guess basis to the main basis.
      ChainMxm(Fock, P12, FockGu, Transpose(P12), Mem);
   } else {
      // project orbitals into main basis.
      FStackMatrix
         OrbAo(nAo, nOcc, &Mem),
         COccO(nAo, 0, &Mem); // dummy.
      FStackMatrix
         S12(nAo, nGu, &Mem);
      MakeIntMatrix(S12, *pOrbBasis, *pGuessBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);
//       ChainMxm(OrbAo, Smh1, Smh1, S12, OrbGu, Mem);
      Mxm(OrbAo, S12, OrbGu);
      CholeskySolve(OrbAo, Scd1);

      if (OrthogonalizeInitialOrbs) {
         OrthSchmidtScd(OrbAo, Scd1, Mem);
         // put back the occupation numbers.
         for (size_t iOcc = 0; iOcc < nOcc; ++ iOcc)
            Scale(&OrbAo(0, iOcc), std::sqrt(fOccGu[iOcc]), size_t(nAo));
      }

      FStackMatrix
         ExchO(nAo, nAo, &Mem);
      Fock.Clear();
      BuildFock(Fock, ExchO, &*pOrbBasis, OrbAo, COccO, 0, FockBuilders, Mem);
      Add(Fock, CoreH1);

//       double f0 = -1.0;
//       if (f0 != 0) {
//          ApplyLevelShift1(Fock, OrbAo, &fOccGu[0], f0, Scd1, Mem);
// //          ApplyLevelShift1(Fock, COccO, 1.0, f0/2, Scd1, Mem);
// //          ApplyLevelShift1(ExchO, COccO, 1.0, f0/2, Scd1, Mem);
//       }
   }
}


} // namespace ct
