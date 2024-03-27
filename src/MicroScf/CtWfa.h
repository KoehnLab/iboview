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

/// @file CtWfa.h
///
/// Core code for interfaces for wave function analysis (WFA) (population
/// analysis, bond order analysis, fragmentation, etc.) and the corresponding
/// analysis reports.
///
/// Note:
/// - This used to be IboView code---this is mostly taken from IvAnalysis.h/.cpp
///   and then split up, so that we can make core reports also in MicroScf
///   directly (for batch processing, etc.).
#ifndef CT_ANALYSIS_H
#define CT_ANALYSIS_H

#include <map>
#include <set>
#include <vector>

#include "CxDefs.h"
#include "CxVec3.h"
#include "CxTiming.h"
#include "CxIo.h"
#include "CxPodArray.h"


#include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "CtMatrix.h"

#include "CtWfi.h" // wave function interface (for communicating wave function state data)


namespace ct {


struct FFrameWfData
{
   FAtomSet
      *pAtoms;
   wfi::FWfType
      WfType;
   FBasisSet
      // basis identifying the functions with respect to which the orbitals and density matrices
      // are expanded. May be either the actual orbital basis, or a minimal basis, in case an IAO
      // representation of the information was requested.
      *pBasis;
   size_t
      // number of functions in pBasis (eiher number of orbital basis functions or number of IAO basis functions)
      nBf,
      nElecA,
      nElecB;
   FMatrixView
      COccA, // alpha
      COccB, // beta
      COccC; // closed
   FMatrixView
      RdmC,
      RdmO;
   size_t nElec() const { return nElecA + nElecB; }
   ptrdiff_t iSpin() const { return ptrdiff_t(nElecA) - ptrdiff_t(nElecB); }
};


typedef std::vector<int>
   FAtomIdList;
typedef std::set<int>
   FAtomIdSet;

struct FAnalysisReport
{
   bool
      m_RestrictAtoms; // if given: exclude atoms not in m_SelectedAtoms.
   FAtomIdList
      m_SelectedAtoms;
   ct::FAtomSet
      &m_Atoms;

   virtual void MakeReport(ct::FLog &Log) = 0;
   bool IsIncluded(int iAt) const;
   std::string AtomLabel(int iAt) const;

   explicit FAnalysisReport(ct::FAtomSet *pAtoms, FAtomIdList *pSelectedAtoms=0);
   ~FAnalysisReport();
};


// encodes the result of a partial charge analysis
struct FChargeAnalysis : public FAnalysisReport
{
   struct FKey {
      int iAt; // atom number (in frame, not element)
      int AngMom; // angular momentum of the AO carrying the current entry.

      FKey() {}
      explicit FKey(int iAt_, int AngMom_) : iAt(iAt_), AngMom(AngMom_) {}
      bool operator < (FKey const &other) const {
         if (iAt < other.iAt) return true;
         if (other.iAt < iAt) return false;
         return AngMom < other.AngMom;
      }
   };
   struct FValue {
      double fCharge;
      double fSpin;
      FValue() {}
      explicit FValue(double fCharge_, double fSpin_) : fCharge(fCharge_), fSpin(fSpin_) {}
   };
   typedef std::map<FKey, FValue>
      FMap;
   FMap
      m_Data;

   void Add(int iAt, int AngMom, double fCharge, double fSpin);
   void MakeReport(ct::FLog &Log); // document only used for making labels..

   int CountAtoms() const;
   int FindMaxL() const;
   bool bNonzeroSpin() const;
   double MakeAtomTotalCharge(int iAt, double fCharge) const;

   explicit FChargeAnalysis(ct::FAtomSet *pAtoms, FAtomIdList *pSelectedAtoms=0);
   ~FChargeAnalysis();
};


struct FBondOrderAnalysis : public FAnalysisReport
{
   struct FKey {
      // atoms this bond order refers to. If jAt = -1, it encodes a valency. Could extend
      // this to multi-center bond orders, by putting in more stuff here.
      int iAt, jAt;

      FKey() {}
      explicit FKey(int iAt_, int jAt_) : iAt(iAt_), jAt(jAt_) {}
      bool operator < (FKey const &other) const {
         if (iAt < other.iAt) return true;
         if (other.iAt < iAt) return false;
         return jAt < other.jAt;
      }
      bool operator != (FKey const &other) const {
         if (iAt != other.iAt) return true;
         if (jAt != other.jAt) return true;
         return false;
      }
   };
   enum FBondOrderType {
      BONDORDER_WibergIao,
//       BONDORDER_Renorm1Iao,
//       BONDORDER_Count = 2
      BONDORDER_Count = 1
   };
   struct FValue {
      // maybe put additional stuff, like perturbative strength? the kinetic elements worked okish.
      double fBondOrders[BONDORDER_Count];
      FValue() {
         for (size_t i = 0; i < BONDORDER_Count; ++ i)
            fBondOrders[i] = 0;
      }
      explicit FValue(double fWibBo_, double fRenormBo_) {
         fBondOrders[BONDORDER_WibergIao] = fWibBo_;
//          fBondOrders[BONDORDER_Renorm1Iao] = fRenormBo_;
         (void)fRenormBo_;
      }
   };
   typedef std::map<FKey, FValue>
      FMap;
   typedef std::pair<FKey, FValue>
      FPair; // why not FMap::value_type? Because "first" is typedef'ed as const, which makes value_type's not copyable.
             // I hate C++. So very, very much. Yes, I am sure there was a good reason to do it exactly like that. In C++ there always is.
   FMap
      m_Data;

   void MakeReport(ct::FLog &Log); // document only used for making labels..
   void MakeBondOrdersForDensity(ct::FRawBasis const *pRawBasis, ct::FMatrixView Rdm, double fOrbOcc, ct::FMemoryStack &Mem);

   FValue const *FindValue(int iAt_, int jAt_) const;

   explicit FBondOrderAnalysis(ct::FAtomSet *pAtoms, ct::FBasisSetPtr pMinBasis, ct::FMatrixView IbRdmC, ct::FMatrixView IbRdmO, ct::FMemoryStack &Mem, FAtomIdList *pSelectedAtoms=0);
   ~FBondOrderAnalysis();
};


struct FFragmentAnalysisOptions : public FIntrusivePtrDest
{
   std::string
      // FDftGrid property list used to describe the target DFT grid.
      GridDesc;
   bool
      // if set, make an individual DFT grid for each fragment to integrate.
      // Otherwise one common full-molecule integration grid will be generated
      // and used to integate the density of every fragment.
      UseFragmentSpecificGrids;
   bool
      GroupRemainingAtoms;
   enum FSpacePartitionType {
      SPACEPARTITION_HilbertSpace_Iao,
      SPACEPARTITION_RealSpace_TFVC
   };
   FSpacePartitionType
      SpacePartitionType;

   explicit FFragmentAnalysisOptions();

#ifdef PROG_IBOVIEW
   void SaveState();
   void LoadState();
#endif // PROG_IBOVIEW
};
typedef ct::TIntrusivePtr<FFragmentAnalysisOptions>
   FFragmentAnalysisOptionsPtr;
typedef ct::TIntrusivePtr<FFragmentAnalysisOptions const>
   FFragmentAnalysisOptionsCptr;


// defines a fragment for a inter-fragment interaction analysis. atm derived from atom groups.
// In principle it might be better to define them by minimal basis shells, atoms might be sufficient for the moment.
struct FFragmentDef : public FIntrusivePtrDest {
   std::string
      m_FragmentName;
   FAtomIdList
      m_iAtoms;
   ct::TArray<size_t>
      // indices of the fragment's minimal basis orbitals in the parent minimal basis
      m_iMinBasisFns;
   ct::FAtomSet const
      *m_pAtoms;
   ct::FBasisSet const
      *m_pMinBasis; // <- hm.. smart ptrs might be better here.

   std::string ShortDesc() const { return fmt::format("{} ({} atoms)", m_FragmentName, m_iAtoms.size()); }

   explicit FFragmentDef(FAtomIdList const &iAtoms_, ct::FAtomSet const *pAtoms_, ct::FBasisSet const *pMinBasis_, std::string const &Name_);
};

typedef ct::TIntrusivePtr<FFragmentDef>
   FFragmentDefPtr;
typedef std::vector<FFragmentDefPtr>
   FFragmentDefPtrList;
struct FFragmentationDef : public FFragmentDefPtrList
{
};

void MakeFragmentationViaSelection(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis, FAtomIdList const &iSelectedAtoms, bool GroupRemainingAtoms);
void MakeMopUpGroupsForFragmentation(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis, FAtomIdSet const &iAssignedAtoms, bool GroupRemainingAtoms);
void MakeFragmentationForIndividualAtoms(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis);


struct FFragmentAnalysis : public FAnalysisReport
{
   FFragmentationDef
      m_FragmentDefs;
   FFragmentAnalysisOptionsCptr
      m_pOptions;

   struct FFragmentData {
      // (..hm... should I allow for fractional electrons in cases of strict degeneracy?)
      size_t nElecAb[2]; // redox charge assignment (integer charge)
      double fElecAb[2]; // actual charge (fractional partial charge)
      size_t nNuclearCharge;

      FFragmentData() { nElecAb[0] = 0; nElecAb[1] = 0; fElecAb[0] = fElecAb[1] = 0.; nNuclearCharge = 0; }
      size_t nElecA() const { return nElecAb[0]; } // number of alpha electrons
      size_t nElecB() const { return nElecAb[1]; } // number of beta electrons

      size_t nElec() const { return nElecA() + nElecB(); } // total number of electrons
      ptrdiff_t iSpin() const { return ptrdiff_t(nElecA()) - ptrdiff_t(nElecB()); } // total spin

      double fElec() const { return fElecAb[0] + fElecAb[1]; } // total number of electrons
      double fSpin() const { return fElecAb[0] - fElecAb[1]; } // total spin

      void ReplicateFromAlpha() { nElecAb[1] = nElecAb[0]; fElecAb[1] = fElecAb[0]; }
      void operator += (FFragmentData const &other) {
         nNuclearCharge += other.nNuclearCharge;
         for (size_t i = 0; i != 2; ++ i) {
            fElecAb[i] += other.fElecAb[i];
            nElecAb[i] += other.nElecAb[i];
         }
      }
   };
   typedef std::vector<FFragmentData>
      FFragmentDataList;
   FFragmentDataList
      m_FragmentData;
   std::string
      m_ElecDistComments[2],
      m_AnalysisTypeComments;
   bool
      m_HaveSpin;
   typedef std::pair<double, size_t>
      FOccPair; // maps occupation number eigenvalue to source fragment.
   typedef TArray<FOccPair>
      FOccPairList;
   ct::FTimer
      tEosaTotal;

//    explicit FFragmentAnalysis(ct::FAtomSet *pAtoms, ct::FBasisSetPtr pMinBasis, ct::FMatrixView IbRdmC, ct::FMatrixView IbRdmO, size_t nElec, ptrdiff_t iSpin, ct::FMemoryStack &Mem, FFragmentationDef const &FragmentDefs_);
   explicit FFragmentAnalysis(ct::FAtomSet *pAtoms, FFrameWfData *pIbWf, FFrameWfData *pAoWf, ct::FMemoryStack &Mem, FFragmentationDef const &FragmentDefs_, FFragmentAnalysisOptionsCptr pOptions_);
   ~FFragmentAnalysis();


   void MakeFragmentOccupationsForSpin_IB(FOccPairList &AllOccs, ct::FMatrixView IbRdmA, ct::FBasisSet const *pMinBasis, size_t iSpinCase, size_t nElecA, ct::FMemoryStack &Mem);
   void MakeFragmentOccupations_TFVC(FOccPairList &AllOccsA, FOccPairList &AllOccsB, FFrameWfData *pAoWf, ct::FMemoryStack &Mem);

//    void MakeDistributionReportForSpin(FOccPairList const &AllOccs, size_t iSpinCase, size_t nElecA, FMemoryStack &Mem);
   void AccRedoxElectronsForSpin(FOccPairList &AllOccs, size_t iSpinCase, size_t nElecA, ct::FMemoryStack &Mem);
   void MakeReport(ct::FLog &Log);
};


} // namespace ct


#endif // CT_ANALYSIS_H
