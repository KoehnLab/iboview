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

#include <iostream> // fixme: remove this. Used for debugging only.
#include <QSettings>
#include "IvAnalysis.h"
#include "IvDocument.h"
#include "IvOrbital.h"

using namespace ct;

// NOTE: this used to be a big and important file... before everything was moved to ~/dev/lhf/CtWfa.cpp and co.


// void FFrame::MakeFragmentationViaAtomGroups(FFragmentationDef &Fragmentation)
// {
//    m_pDocument->MakeFragmentationViaAtomGroups(Fragmentation, pGetAtoms(), m_pMinBasis.get());
// }

void FFrame::MakeFragmentationViaAtomGroups(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis, bool GroupRemainingAtoms)
{
   m_pDocument->MakeFragmentationViaAtomGroups(Fragmentation, pAtoms, pMinBasis, GroupRemainingAtoms);
}




size_t FDocument::nAtomGroupsDefined() const
{
   return m_AtomGroups.size();
}


void FDocument::MakeFragmentationViaAtomGroups(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis, bool GroupRemainingAtoms)
{
   Fragmentation.clear();

   FAtomIdSet
      // atoms assigned by the current groups until now.
      iAssignedAtoms;
//    int
//       iLargestGroupId; // 1 + this will be used for "rest"

   FAtomGroupMap::const_iterator
      itGroup;
   for (itGroup = m_AtomGroups.begin(); itGroup != m_AtomGroups.end(); ++ itGroup) {
      FAtomIdList const
         &iGroupAtoms = itGroup->second;
      FAtomIdList
         FragmentList;
      if (iGroupAtoms.empty())
         // don't generate groups without atoms.
         continue;
      for (size_t i = 0; i < iGroupAtoms.size(); ++ i) {
         int
            iAt = iGroupAtoms[i];
         if (size_t(iAt) >= pAtoms->size())
            continue; // no such atom in the current frame.
         if (iAssignedAtoms.find(iAt) != iAssignedAtoms.end()) {
            IV_NOTIFY1(NOTIFY_Warning, "Atom %1 occurs in multiple atom groups. Fragmentation invalid!", iAt+1);
         } else {
            FragmentList.push_back(iAt);
            iAssignedAtoms.insert(iAt);
         }
      }

      Fragmentation.push_back(new FFragmentDef(FragmentList, pAtoms, pMinBasis, fmt::format("Atom Group {}", itGroup->first)));
   }

   // unless given group list is complete, create a mop up group for the remaining atoms not yet assigned anywhere.
   MakeMopUpGroupsForFragmentation(Fragmentation, pAtoms, pMinBasis, iAssignedAtoms, GroupRemainingAtoms);
//    if (iAssignedAtoms.size() != pAtoms->size()) {
//       if (GroupRemainingAtoms) {
//          FAtomIdList
//             FragmentList;
//          for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt) {
//             // atom already assigned to a group?
//             if (iAssignedAtoms.find(iAt) != iAssignedAtoms.end())
//                continue;
//             FragmentList.push_back(iAt);
//          }
//          Fragmentation.push_back(new FFragmentDef(FragmentList, pAtoms, pMinBasis, fmt::format("Other Atoms")));
//       } else {
//          // or, alternatively, make individual groups for all remaining atoms.
//          for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt) {
//             // atom already assigned to a group?
//             if (iAssignedAtoms.find(iAt) != iAssignedAtoms.end())
//                continue;
//             FAtomIdList
//                FragmentList;
//             FragmentList.push_back(iAt);
//             Fragmentation.push_back(new FFragmentDef(FragmentList, pAtoms, pMinBasis, fmt::format((*pAtoms)[iAt].GetAtomLabel(iAt))));
//          }
//       }
//    }
}


// ORBMAT_IaoBasis
FFrameWfData *FFrame::MakeWfDataForAnalysis(ct::FLog &Log, unsigned WfInfoFlags, unsigned OrbmatType, FMemoryStack *pMem)
{
   assert_rt(pMem != 0);
   FFrameWfData
      *r;
   pMem->Alloc(r);

//    TMemoryLock<double>
//       pFreeMe(0, pMem);

   r->pAtoms = pGetAtoms();
   if (r->pAtoms == 0) {
      Log.Write(" No atoms.");
      return 0;
   }
   r->WfType = GetWfType();
   if (r->WfType != wfi::WFTYPE_Rhf && r->WfType != wfi::WFTYPE_Uhf) {
      Log.Write("Sorry, wave function type not supported. Can only do RHF/RKS and UHF/UKS at this moment.");
      return 0;
   }

   FBasisSetPtr
      pBasis1;

   // assemble the closed- and open-shell 1 RDMs in the target basis.
   {
      TArray<FOrbital*>
         RefOrbitalsA,
         RefOrbitalsB,
         RefOrbitalsC;
      MakeOrbitalMatrix(r->COccA, pBasis1, &RefOrbitalsA, ORBMAT_OccupiedOnly | ORBMAT_AlphaOnly | OrbmatType, *pMem); assert(RefOrbitalsA.size() == r->COccA.nCols);
      MakeOrbitalMatrix(r->COccB, pBasis1, &RefOrbitalsB, ORBMAT_OccupiedOnly | ORBMAT_BetaOnly | OrbmatType, *pMem); assert(RefOrbitalsB.size() == r->COccB.nCols);
      MakeOrbitalMatrix(r->COccC, pBasis1, &RefOrbitalsC, ORBMAT_OccupiedOnly | ORBMAT_ClosedOnly | OrbmatType, *pMem); assert(RefOrbitalsC.size() == r->COccC.nCols);
   }
   r->pBasis = pBasis1.get();

   if (pBasis1.get() == 0) {
      IvNotify(NOTIFY_Error, "Failed to assign minimal basis for bond order analysis. No orbitals/broken orbital occupation?");
      return 0;
   }

   r->nBf = pBasis1->nFn();
   r->nElecA = 0;
   r->nElecB = 0;

   // make sure that the number of rows of all submatrices matches the expected size
   // according to the given basis set, even if no actual orbitals of the corresponding
   // type exist. This simplifies sanity checks later on.
   if (r->COccA.nRows == 0) { r->COccA.nRows = r->nBf; }
   if (r->COccB.nRows == 0) { r->COccB.nRows = r->nBf; }
   if (r->COccC.nRows == 0) { r->COccC.nRows = r->nBf; }

   if (WfInfoFlags & WFINFO_MakeRdms) {
      r->RdmC = MakeStackMatrix(r->nBf, r->nBf, *pMem);
      r->RdmO = MakeStackMatrix(r->nBf, r->nBf, *pMem);
      r->RdmC.Clear();
      r->RdmO.Clear();
   } else {
      r->RdmC = FMatrixView(0,0,0,0,0);
      r->RdmO = FMatrixView(0,0,0,0,0);
   }

   if (r->WfType == wfi::WFTYPE_Uhf) {
      if (r->COccC.nCols != 0)
         IvNotify(NOTIFY_Error, IvFmt("Found %1 closed-shell orbitals in a wave function which was assumed to be UHF/UKS type. Something went wrong!", r->COccC.nCols));
      // note: that is technically not quite 100% correct, as the alpha and beta IAO bases may be slightly different
      // due to different polarization contributions.
      // We here assume that the contributions stiff refer to the same "effective" atomic valence orbitals.
      if (WfInfoFlags & WFINFO_MakeRdms) {
         if (r->COccA.nCols != 0)
            SyrkNT(r->RdmC, r->COccA);
         if (r->COccB.nCols != 0)
            SyrkNT(r->RdmC, r->COccB, 1., MXM_Add);
         if (r->COccA.nCols != 0)
            SyrkNT(r->RdmO, r->COccA);
         if (r->COccB.nCols != 0)
            SyrkNT(r->RdmO, r->COccB, -1., MXM_Add);
      }
      r->nElecA = r->COccA.nCols;
      r->nElecB = r->COccB.nCols;
   } else {
      if (r->COccB.nCols != 0)
         IvNotify(NOTIFY_Error, IvFmt("Found %1 beta-spin orbitals in a wave function which was assumed to be RHF/RKS type. Something went wrong!", r->COccB.nCols));
      // here it is okay, as the IAO basis is made from the span of all occupied orbials.
      if (WfInfoFlags & WFINFO_MakeRdms) {
         if (r->COccC.nCols != 0)
            SyrkNT(r->RdmC, r->COccC, 2.);
         if (r->COccA.nCols != 0)
            SyrkNT(r->RdmC, r->COccA, 1., MXM_Add);
         if (r->COccA.nCols != 0)
            SyrkNT(r->RdmO, r->COccA);
      }
      r->nElecA = r->COccC.nCols + r->COccA.nCols;
      r->nElecB = r->COccC.nCols;
   }

   return r;
}


void FFrame::RunBondOrderAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms, FMemoryStack *pMem)
{
   assert_rt(pMem != 0);
   TMemoryLock<double>
      pFreeMe(0, pMem);


   // collect orbitals and make RDMs in IAO basis. All output structures go onto pMem.
   FFrameWfData
      *pIbWf = MakeWfDataForAnalysis(Log, WFINFO_MakeRdms, ORBMAT_IaoBasis, pMem);
   FFrameWfData
      *pAoWf = MakeWfDataForAnalysis(Log, 0, 0, pMem);
   if (pIbWf == 0 || pAoWf == 0)
      return;

   FBondOrderAnalysis
      BondOrderAnalysis(pIbWf->pAtoms, pIbWf->pBasis, pIbWf->RdmC, pIbWf->RdmO, *pMem, pSelectedAtoms);

   BondOrderAnalysis.MakeReport(Log);

   FGeometry
      *pGeometry = pGetGeometry();
   if (pGeometry)
      pGeometry->FindBondLines(m_pDocument, &BondOrderAnalysis);

//    if (m_pDocument->nAtomGroupsDefined() > 0) {
//       FFragmentAnalysisOptionsPtr
//          pFragmentAnalysisOptions = new FFragmentAnalysisOptions();
//       FFragmentationDef
//          Fragmentation;
//       MakeFragmentationViaAtomGroups(Fragmentation, pGetAtoms(), &*pIbWf->pBasis);
//       FFragmentAnalysis
//          FragmentAnalysis(pGetAtoms(), pIbWf, pAoWf, *pMem, Fragmentation, pFragmentAnalysisOptions);
//       FragmentAnalysis.MakeReport(Log);
//    }
}


void FFrame::RunRedoxChargeAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms, FFragmentAnalysisOptionsCptr pFragmentAnalysisOptions, FMemoryStack *pMem)
{
   assert_rt(pMem != 0);
   TMemoryLock<double>
      pFreeMe(0, pMem);


   // collect orbitals and make RDMs in IAO basis. All output structures go onto pMem.
   FFrameWfData
      *pIbWf = MakeWfDataForAnalysis(Log, WFINFO_MakeRdms, ORBMAT_IaoBasis, pMem);
   FFrameWfData
      *pAoWf = MakeWfDataForAnalysis(Log, 0, 0, pMem);
   if (pIbWf == 0 || pAoWf == 0)
      return;

   FFragmentationDef
      Fragmentation;
   if (pSelectedAtoms && pSelectedAtoms->size() != 0) {
      MakeFragmentationViaSelection(Fragmentation, pGetAtoms(), &*pIbWf->pBasis, *pSelectedAtoms, pFragmentAnalysisOptions->GroupRemainingAtoms);
   } else if (m_pDocument->nAtomGroupsDefined() > 0) {
      MakeFragmentationViaAtomGroups(Fragmentation, pGetAtoms(), &*pIbWf->pBasis, pFragmentAnalysisOptions->GroupRemainingAtoms);
   } else {
      MakeFragmentationForIndividualAtoms(Fragmentation, pGetAtoms(), &*pIbWf->pBasis);
   }
   FFragmentAnalysis
      FragmentAnalysis(pGetAtoms(), pIbWf, pAoWf, *pMem, Fragmentation, pFragmentAnalysisOptions);
   FragmentAnalysis.MakeReport(Log);
}


void FFrame::RunChargeAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms)
{
   FAtomSet
      *pAtoms = pGetAtoms();
   if (pAtoms == 0) {
      Log.Write(" No atoms.");
      return;
   }
   FChargeAnalysis
      ChargeAnalysis(pAtoms, pSelectedAtoms);
   FDataSetList::iterator
      itDataSet;
   for (itDataSet = m_Data.begin(); itDataSet != m_Data.end(); ++ itDataSet) {
      FOrbital
         *pOrb = dynamic_cast<FOrbital*>(itDataSet->get());
      if (!pOrb)
         continue;
      pOrb->AddChargeContributions(ChargeAnalysis);
   }
   ChargeAnalysis.MakeReport(Log);
}
