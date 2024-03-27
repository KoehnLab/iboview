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

#include "CtBasisDesc.h"

#include "format.h"
#include "CxParse1.h"
#include "Ir.h"
#ifdef IR_ECP
   #include "IrEcp.h"
#endif // IR_ECP

#include <iostream> // FIXME: remove this.

namespace ct {

static bool s_BasisAssoocDebugPrint = false;
// static bool s_BasisAssoocDebugPrint = true;

static bool _in_range(int i, int iMin, int iMax) { return iMin <= i && i <= iMax; }

std::string GetBasisContextName(FBasisContext Context)
{
   switch(Context) {
      case BASIS_Orbital: return "ORBITAL";
      case BASIS_Guess: return "GUESS";
      case BASIS_MinAo: return "MINAO";
      case BASIS_JFit: return "JFIT";
      case BASIS_JkFit: return "JKFIT";
      case BASIS_Mp2Fit: return "MP2FIT";
      case BASIS_CcsdFit: return "EXTFIT";
      case BASIS_F12RI: return "F12RI";
      case BASIS_IaoFit: return "IAOFIT";
      case BASIS_Multipoles: return "MULTIPOLE";
default:
         assert(0);
         return "[unknown basis context]";
   }
}


bool NeedOrbitalBasisForContextDefault(FBasisContext Context)
{
   static FBasisContext const
      ContextsToIgnoreOrbs[] = {BASIS_Guess, BASIS_MinAo, BASIS_JFit, BASIS_JkFit, BASIS_Multipoles};
   for (size_t i = 0; i < sizeof(ContextsToIgnoreOrbs)/sizeof(ContextsToIgnoreOrbs[0]); ++ i)
      if (Context == ContextsToIgnoreOrbs[i])
         return false;
   return true;
}


void FindDefaultBasis(std::string &Out, FBasisContext Context, std::string const &OrbBasis, FBasisDescs const *pOtherBases, int iElement, unsigned nEcpElec)
{
   string_slice
      slOrbBasis(OrbBasis);
   if (Context == BASIS_Orbital)
      throw std::runtime_error("FindOrbitalBasis: not allowed for BASIS_Orbital context.");
   if (Context == BASIS_MinAo) {
      // no BASIS_MinAo has been declared. But maybe a BASIS_Guess? I split
      // them only in 2020, and most programs and inputs only declare Guess.
      if (pOtherBases) {
         FBasisDescs::const_iterator
            itDesc = pOtherBases->find(BASIS_Guess);
         if (itDesc != pOtherBases->end()) {
            // got a GUESS basis. Return that one for the MINAO context.
            Out = itDesc->second;
            return;
         }
      }
   }
   if (Context == BASIS_Guess || Context == BASIS_MinAo) {
      // should possibly check for ECP-aware guess basis. This one is just a guess...
      if (nEcpElec == 0)
         Out = "MINAO";
      else
         Out = "MINAO-PP";
      return;
   } else if (Context == BASIS_JkFit) {
      Out = "univ-JKFIT";
      return;
   } else if (Context == BASIS_JFit) {
      if (nEcpElec == 0)
         Out = "univ-JFIT";
      else
         Out = fmt::format("univ-ecp{}-JFIT", nEcpElec);
      return;
   } else if (Context == BASIS_Mp2Fit) {
      if (slOrbBasis.startswith("MINAO")) {
         Out = "def2-SVP-MP2FIT";
         return;
      }
   } else {
      Out = fmt::format("{}-{}", OrbBasis, GetBasisContextName(Context));
   }
   IR_SUPPRESS_UNUSED_WARNING(iElement); // might be required later on.
}



static void FixDef2BasisAndEcpDecl(std::string &OrbBasis, std::string &EcpName, unsigned &nEcpElec, int iElement)
{
   if (s_BasisAssoocDebugPrint)
      std::cout << fmt::format("    :: basis-desc: fix-def2-assoc({},'{}')", iElement, OrbBasis);

   // remove def2/dhf prefixes. Will be added back later.
   if (StartsWith(OrbBasis, "def2-")) OrbBasis = OrbBasis.substr(5);
   if (StartsWith(OrbBasis, "dhf-")) OrbBasis = OrbBasis.substr(4);

   if (!(StartsWith(OrbBasis, "SV") || StartsWith(OrbBasis, "TZV") || StartsWith(OrbBasis, "QZV")))
      throw std::runtime_error("logic error: Failed to recognize Weigend/Ahlrichs-style basis type.");
   if (iElement > 86) throw std::runtime_error("sorry, no default element associations beyond element 86 programmed.");

   // replace some def2 basis sets by the updated dhf basis sets, which are
   // associated with MDF ECPs which were not yet available when the original
   // def2 sets were made (see: Weigend, Baldes: J. Chem. Phys. 133, 174102
   // (2010)).
   //
   // Notes:
   //   - Xe (54) had already MDF ECP in original def2 sets
   //   - La was in def2, but has no MDF ECP as of 2015.
   //
   // With these replacements, all elements in range 36-86, except for Lanthan (which has a def2 set) and the lanthanides (which have no def2 sets)
   // are now associated with MDF ECPs. La is still associated with a MWB ECP.
   if ( _in_range(iElement, 37, 48) ||
       (_in_range(iElement, 53, 80) && iElement != 57 && iElement != 54)) {
      OrbBasis = "dhf-" + OrbBasis;
   } else {
      OrbBasis = "def2-" + OrbBasis;
   }
   if (_in_range(iElement, 72, 86)) { EcpName = "ECP60MDF"; nEcpElec = 60; }
   if (_in_range(iElement, 55, 71) && iElement != 57) { EcpName = "ECP46MDF"; nEcpElec = 46; }
   if (_in_range(iElement, 57, 57)) { EcpName = "ECP46MWB"; nEcpElec = 46; } // La: no MDF ECP there, and no dhf set (as of 2015).
   if (_in_range(iElement, 37, 54)) { EcpName = "ECP28MDF"; nEcpElec = 28; }
   // note: There are also ECP10MDFs for elements between K and Kr, but the def2 sets do not use them (yet?).
   // They are all-electron sets for those. If dhf sets for those come up, feel free to contact me (Gerald Knizia)
   // and I will update MicroScf.
   if (s_BasisAssoocDebugPrint)
      std::cout << fmt::format(" --> {{basis: {}, ecp: {}, nEcpElec: {} }}", OrbBasis, EcpName, nEcpElec) << std::endl;
}


// "If this element had an ECP, how many electrons would the ECP likely represent?"
size_t iLikelyEcpCharge(int iElement)
{
   // these are for ECP10MDFs used with vtz-pp. def2 sets are always all-electron for those.
   // And at the time of writing, the PP sets are not actually there for 21-28.
   if (_in_range(iElement,21,36)) return 10;

   // these are def2-nzvpp/dhf-nvzpp associations (mostly ECPxxMDF, but some MWB are left).
   if (_in_range(iElement,1,36)) return 0;
   if (_in_range(iElement,37,54)) return 28;
   if (_in_range(iElement,55,71)) return 46; // I guess? somehow I don't seem to have sets for lanthanides (except lanthan)
   if (_in_range(iElement,72,86)) return 60;

   return 0;
}


void AutoAdjustBasisAssociations(FBasisDescs &BasisDesc, std::string &EcpDesc, unsigned &nEcpElec, unsigned iElement)
{
   // fix up some basis declarations; E.g., replace generic terms like "TZVPP" by their explicit def2-/dhf- basis names,
   // and change MINAO to a suitable minimal basis set compatible with a given ECP.
   std::string
      OrbBasis = tolower(BasisDesc[BASIS_Orbital]);
   if (OrbBasis.empty()) {
      nEcpElec = 0;
      EcpDesc.clear();
      return; // no basis set. Can't actually do anything here.
   }
   if (s_BasisAssoocDebugPrint)
      std::cout << fmt::format("    :: basis-desc: auto-adjust({},'{}','{}')\n", iElement, OrbBasis, EcpDesc);

   // Mote on Stuttgart/Cologne/Bonn ECPs:
   // - Both the def2-/dhf- and the cc-pVnZ sets employ Stuttgart/Cologne/Bonn
   //   ECPs. These are also the only ones we ship with MicroScf.
   //
   // - The names comes as ECPnnXYY, where
   //    + n = number of core electrons absorbed by the ECP
   //    + X = S/M: single/multi electron fit
   //    + Y = HF/WB/DF: non/quasi/fully relativistic
   //
   // - In practice, we never want to use SHF/MHF ECPs at all. We'd like
   //   to use MDF/SDF ECPs if they are there, and SWB/MWB ECPs as fallback.
   //   Although in some cases, the SDF basis sets are apparently useful
   //   for singly ionized atoms (e.g., for alkali metals or so). But I
   //   have not fully understood the usage cases for that.
   //
   // - We are also constrainted by the basis sets: Basis sets *need*
   //   to be used with the ECPs they are made for! (or at the very least, the
   //   s and p shells needs to be uncontrated if not, and even that is only
   //   likely to work with related ECPs (e.g., a ECPnnMWB vs ECPnnSDF))
   //
   // - So we here do not give a lot of freedom, but simply assign the ECPs
   //   for the def2-/dhf- basis sets they are meant to be used with, add some
   //   similar definitions for cc-pVnZ sets, and otherwise expect that the correct
   //   ECPs are manually assigned from the outside.

   if (StartsWith(OrbBasis, "def2-") || StartsWith(OrbBasis, "dhf-") || StartsWith(OrbBasis, "sv") || StartsWith(OrbBasis, "tzv") || StartsWith(OrbBasis, "qzv"))
      FixDef2BasisAndEcpDecl(BasisDesc[BASIS_Orbital], EcpDesc, nEcpElec, iElement);
   else if (StartsWith(OrbBasis, "cc-p") || StartsWith(OrbBasis, "aug-cc-p") || StartsWith(OrbBasis, "minao")) {
      if (_in_range(iElement, 19, 36)) {
         // check for assignment of post-2000s cc-pVnZ-PP basis sets for
         // elements K -- Br, which are meant to be used with ECP10MDFs.
         //
         // Note: At the time of writing this, not all of those basis sets are
         // actually published, and we have not tabulated all of them for
         // MicroScf. I think we should still assign them, to make the program
         // generate an error, instead of silently using an incorrect basis.
         bool EndsInPP = EndsWith(OrbBasis, "-pp");
         bool UseEcp10 = false;
         if (nEcpElec == 10 || EndsInPP)
            UseEcp10 = true;
         if (UseEcp10) {
            nEcpElec = 10;
            EcpDesc = "ECP10MDF";
            if (!EndsInPP)
               // add the -pp to clarify what we chose.
               OrbBasis += "-pp";
         } else {
            nEcpElec = 0;
            EcpDesc.clear();
         }
      } else if (iElement > 36) {
         nEcpElec = iLikelyEcpCharge(iElement);
         if (!EndsWith(OrbBasis, "-pp"))
            OrbBasis += "-pp";
         // Default associations for bases which actually exist
         // at the time of wrinting:
         //
         // cc-pVDZ-PP 19 20 --> ECP : ECP10MDF
         // cc-pVDZ-PP 29 36 --> ECP : ECP10MDF
         // cc-pVDZ-PP 37 54 --> ECP : ECP28MDF
         // cc-pVDZ-PP 55 56 --> ECP : ECP46MDF
         // cc-pVDZ-PP 72 86 --> ECP : ECP60MDF
         // cc-pVDZ-PP 87 88 --> ECP : ECP78MDF
         // cc-pVDZ-PP 90 92 --> ECP : ECP60MDF
         if (_in_range(iElement, 37, 54)) {
            EcpDesc = "ECP28MDF"; nEcpElec = 28;
         } else if (_in_range(iElement, 55, 56)) { // La (Lanthan) and Ce (Cerium)
            EcpDesc = "ECP46MDF"; nEcpElec = 46;
         } else if (_in_range(iElement, 57, 71)) {
            // these are the remaining Lanthanides. At this moment we have neither
            // ECPs nor bases for those, but presumably if there were bases, that
            // would be ECP46MDF ones.
            EcpDesc = "ECP46MDF"; nEcpElec = 46;
         } else if (_in_range(iElement, 72, 86)) {
            EcpDesc = "ECP60MDF"; nEcpElec = 60;
         } else if (_in_range(iElement, 87, 88)) {
            EcpDesc = "ECP78MDF"; nEcpElec = 78;
         } else if (_in_range(iElement, 90, 92)) {
            EcpDesc = "ECP60MDF"; nEcpElec = 60;
         } else {
            EcpDesc = "CC-ECP?-FIXME";
         }

      }
      if (iElement == 1 || iElement == 2) {
         // if a core-correlation (cc-pCVnZ) or weighted core-valence
         // correlation (cc-pwCVnZ) orbital basis set was specified,
         // automatically replace it by the corresponding valence basis for H
         // and He, which for obvious reasons have no core-correlation sets.
         string_slice
            sl(OrbBasis);
         cstr_it
            it;
         it = sl.find("cc-pwcv");
         if (it != sl.end()) {
            assert_rt(OrbBasis.size() == sl.size());
            BasisDesc[BASIS_Orbital].erase(it - sl.first + 4, 2); // delete the "wC"
            OrbBasis = tolower(BasisDesc[BASIS_Orbital]);
         }
         it = sl.find("cc-pcv");
         if (it != sl.end()) {
            assert_rt(OrbBasis.size() == sl.size());
            BasisDesc[BASIS_Orbital].erase(it - sl.first + 4, 1); // delete the "C"
            OrbBasis = tolower(BasisDesc[BASIS_Orbital]);
         }
      }
   }

   if (nEcpElec != 0 || iElement > 36) {
//          xout << "ELEMENT: " << ElementName_ << "  nEcpElectrons: " << nEcpElectrons << "  iLikelyEcpCharge: " << iLikelyEcpCharge(iElement) << "  iElement: " << iElement << std::endl;
      // patch up MINAO to MINAO-PP.
      if (IsEqual_CaseInsensitive(BasisDesc[BASIS_Guess], "MINAO")) {
#ifndef PROG_IBOVIEW
         // is the number of ECP electrons we were given consistent with the VTZ-PP derived MINAO-PP basis sets?
         if (iLikelyEcpCharge(iElement) != nEcpElec && !StartsWith(OrbBasis, "!/"))  // !/: that's explicitly specified basis files in wmme.py export...
            throw std::runtime_error(fmt::format("Sorry, failed to assign ECP-dependent MINAO basis set for element {}: Declared ECP '{}' has {} electrons absorbed, but MINAO-PP expects {} absorbed electrons. Please specify MINAO set explicitly.", iElement, EcpDesc, nEcpElec, iLikelyEcpCharge(iElement)));
         // ^- FIXME: put this back in.
#endif // PROG_IBOVIEW
         // ^- FIXME: hm.. why did I turn this off again? I thought there was
         //    some reason for that...  maybe some IboView related issue with
         //    using or not using ECP10MDFs or something? I think in IboView
         //    it could cause more problems by incorrectly being enabled. So
         //    I'll switch it off there (only) for the time being...
         BasisDesc[BASIS_Guess] = "MINAO-PP";
         nEcpElec = iLikelyEcpCharge(iElement);
      }
      if (IsEqual_CaseInsensitive(BasisDesc[BASIS_JFit], "univ-JFIT")) {
         // unlike the JKFIT variants, these ones have the ECP electron number
         // they are meant to be used with as part of their name (e.g., univ-
         // ecp28-JFIT). ...at least in the way I generated them. Note that
         // doing so *does* make sense, because there could also be similar
         // fitting sets meant to be used without ECPs or smaller/larger ECPs
         // than the default ones.
         BasisDesc[BASIS_JFit] = fmt::format("univ-ecp{}-JFIT", nEcpElec);
      }
   }
//       if (1) {
//          std::cout << fmt::format(" *********  {:2}  nEcpElec = {}  OrbBas = {}  MINAO = {}", ElementName_, nEcpElectrons, BasisDesc[BASIS_Orbital], BasisDesc[BASIS_Guess]) << std::endl;
//       }
}



static void DbgWriteEcpLine(fmt::MemoryWriter &w, std::string const &Cap, double const *pData, size_t nExp, std::string const &Ind)
{
   w.write("{}       {:6}", Ind, Cap + ":");
   for (size_t iExp = 0; iExp != nExp; ++ iExp)
      w.write(" {:12.6f}", pData[iExp]);
//    w.write("\n");
   w << "\n";
}

void PrintEcpData(std::ostream &out, ir::FAtomEcp const *pIrEcp, int iElement, TVector3<double> const &vPos, std::string const &Ind)
{
   fmt::MemoryWriter
      w;
//    w.write("{}ECP for {} at [{}, {}, {}]\n", Ind, pAt->GetElementName(), pAt->vPos[0], pAt->vPos[1], pAt->vPos[2]);
   w.write("{}ECP for {} at [{}, {}, {}]\n", Ind, iElement, vPos[0], vPos[1], vPos[2]);
   w.write("{}\n", Ind);
#ifdef IR_ECP
   w.write("{}  MaxL = {}  MaxPowR = {}\n", Ind, pIrEcp->MaxL(), pIrEcp->MaxPowR());
   for (size_t iSh = 0; iSh != pIrEcp->ShellFns.size(); ++ iSh) {
      ir::FEcpShellFn const
         *pEcpSh = &pIrEcp->ShellFns[iSh];
      w.write("{}  #{}:  l = {}  nExp = {}\n", Ind, iSh, pEcpSh->l, pEcpSh->nExp);
      DbgWriteEcpLine(w, "Exp", pEcpSh->pExp, pEcpSh->nExp, Ind);
      DbgWriteEcpLine(w, "Co", pEcpSh->pCo, pEcpSh->nExp, Ind);
      DbgWriteEcpLine(w, "r^p", pEcpSh->pPow, pEcpSh->nExp, Ind);
   }
#endif
//    w.write("\n");
   // ^- ancient g++ on ubuntu doesn't like that.
//    w << "\n";
   out << w.str();
}


} // namespace ct
