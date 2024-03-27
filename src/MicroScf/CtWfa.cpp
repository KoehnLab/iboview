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
#include <algorithm> // for std::sort and co

#ifdef PROG_IBOVIEW
   #include <QSettings>
#endif // PROG_IBOVIEW

#include "CtWfa.h"
#include "CtVoronoiPartition.h"
#include "CtVolumeProperty.h"


#include "CxArgSort.h" // for ArgMax/ArgMin/ArgSort
#include "CxPoly.h"
// #include "CxColor.h"
#include "CxOpenMpProxy.h"
#include "CxAtomData.h"
#include "format.h"

#include "CxIo.h"
#include "CxPhysicalUnits.h"
#include "CtMatrix.h"

#include "CtInt1e.h"
#include "IrAmrr.h"
#include "CtMatrix.h" // for ArgSort1
#include "CtRhf.h"


#include "CtIvStubs.h"


namespace ct {

FAnalysisReport::FAnalysisReport(ct::FAtomSet *pAtoms, FAtomIdList *pSelectedAtoms)
   : m_RestrictAtoms(pSelectedAtoms != 0),
     m_Atoms(*pAtoms)
{
   if (m_RestrictAtoms) {
      m_SelectedAtoms.insert(m_SelectedAtoms.end(), pSelectedAtoms->begin(), pSelectedAtoms->end());
      std::sort(m_SelectedAtoms.begin(), m_SelectedAtoms.end());
   }

}

FAnalysisReport::~FAnalysisReport()
{}

bool FAnalysisReport::IsIncluded(int iAt) const
{
   if (!m_RestrictAtoms)
      return true;
   FAtomIdList::const_iterator
      it = std::lower_bound(m_SelectedAtoms.begin(), m_SelectedAtoms.end(), iAt);
   return (it != m_SelectedAtoms.end() && iAt == *it);
}

std::string FAnalysisReport::AtomLabel(int iAt) const
{
   fmt::MemoryWriter
      out;
   if (iAt < 0 || size_t(iAt) >= m_Atoms.size())
      out << "?? ";
   else
      out.write("{:>2} ", m_Atoms[size_t(iAt)].GetElementName());
   out.write("{:3}", iAt+1);
   return out.str();
}



FChargeAnalysis::FChargeAnalysis(ct::FAtomSet *pAtoms, FAtomIdList *pSelectedAtoms)
   : FAnalysisReport(pAtoms, pSelectedAtoms)
{}

FChargeAnalysis::~FChargeAnalysis()
{}


int FChargeAnalysis::CountAtoms() const
{
   std::set<int>
      iAtoms;
   for (FMap::const_iterator it = m_Data.begin(); it != m_Data.end(); ++ it) {
      iAtoms.insert(it->first.iAt);
   }
   return int(iAtoms.size());
}

int FChargeAnalysis::FindMaxL() const
{
   std::set<int>
      iAngMoms;
   for (FMap::const_iterator it = m_Data.begin(); it != m_Data.end(); ++ it) {
      iAngMoms.insert(it->first.AngMom);
   }
   if (iAngMoms.empty())
      return 0;
   return *iAngMoms.rbegin();
}

bool FChargeAnalysis::bNonzeroSpin() const
{
   for (FMap::const_iterator it = m_Data.begin(); it != m_Data.end(); ++ it)
      if (it->second.fSpin != 0.)
         return true;
   return false;
}

void FChargeAnalysis::Add(int iAt, int AngMom, double fCharge, double fSpin)
{
   if (!IsIncluded(iAt))
      return;
   assert(AngMom >= 0);
   if (AngMom < 0)
      AngMom = 0;
   FKey
      key(iAt, AngMom);
   FMap::iterator
      it = m_Data.find(key);
   if (it == m_Data.end()) {
      std::pair<FMap::iterator, bool>
         itb = m_Data.insert(FMap::value_type(key, FValue(0.,0.)));
      it = itb.first;
      assert(it != m_Data.end());
   }

   it->second.fCharge += fCharge;
   it->second.fSpin += fSpin;
}

double FChargeAnalysis::MakeAtomTotalCharge(int iAt, double fCharge) const
{
   FAtom const
      *pAt = 0;
   double
      fAtomElem = 0;
   if (iAt >= 0 && size_t(iAt) < m_Atoms.size()) {
      pAt = &m_Atoms[size_t(iAt)];
//       fAtomElem = double(pAt->NuclearCharge());
      fAtomElem = double(pAt->iElement);
   } else {
      return -fCharge;
   }
   double
      fReducedCharge = fAtomElem - fCharge;
   // ^- does not take account of ECPs.
   if (pAt->nEcpElec != 0) {
      // if a certain number of ECP electrons is explicitly declared,
      // then use it.
      fReducedCharge -= pAt->nEcpElec;
   } else {
      // no ECP electrons declared, but this might just be because
      // we got the data from an orbital file which does not say anything
      // about ECPs.

      // So a apply a hack to guess the right
      //    ECP charge...  we just take the charge which is closest to 0.
      //    (of course this won't work with ECP2s etc. But the import formats
      //    often do not have any ECP information, so without it we cannot
      //    really do better.)
      double
         fReducedChargeWithPutativeEcp = fReducedCharge - (double)iLikelyEcpCharge(m_Atoms[iAt].iElement);
      if (std::abs(fReducedChargeWithPutativeEcp) < std::abs(fReducedCharge))
         fReducedCharge = fReducedChargeWithPutativeEcp;
   }
   return fReducedCharge;
}

void FChargeAnalysis::MakeReport(ct::FLog &Log)
{
   int
      nMaxL = FindMaxL(),
      nAtoms = CountAtoms();
   bool
      HaveSpin = bNonzeroSpin();
   for (int iSpinCase = 0; iSpinCase != 1 + int(HaveSpin); ++ iSpinCase)
   {
      double
         fElecTotal = 0.,
         fChargeTotal = 0.;
      if (iSpinCase == 0)
         Log.Write(" Total charge composition:\n");
      else
         Log.Write(" Spin density composition:\n");

      {  // make table caption
         fmt::MemoryWriter
            out;
         out << "   CEN ATOM  ";
         for (int i = 0; i <= nMaxL; ++ i)
            out.write("       {}   ", "spdfghiklm"[i]);
         out << "    ELECTRONS";
         if (iSpinCase == 0)
            out << " ->P.CHARGE";
         else
            out << " ->SPIN.DEN";
         Log.Write(out.str());
      }
      char const
         *pChargeFmtF = " {:10.5f}",
         *pChargeFmtS = " {:10}";
      std::vector<double>
         AmCharges(size_t(nMaxL+1), 0.);
      FMap::const_iterator
         itBeg, itEnd;
      for (itBeg = m_Data.begin(); itBeg != m_Data.end(); ) {
         // clear per-angular momentum charge data.
         for (size_t i = 0; i < AmCharges.size(); ++ i)
            AmCharges[i] = 0.;

         int
            iAt = itBeg->first.iAt,
            iAtMaxL = 0;
         double
            fAtomTotal = 0.;
         // find all angular momentum entries on the current atom.
         itEnd = itBeg;
         while (itEnd != m_Data.end() && itEnd->first.iAt == iAt) {
            double f = (iSpinCase==0)? itEnd->second.fCharge : itEnd->second.fSpin;
            iAtMaxL = std::max(iAtMaxL, int(itEnd->first.AngMom));
            AmCharges[size_t(itEnd->first.AngMom)] += f;
            fAtomTotal += f;
            ++ itEnd;
         }

         fElecTotal += fAtomTotal;
         double
            fAtomCharge;
         if (iSpinCase == 0)
            fAtomCharge = MakeAtomTotalCharge(iAt, fAtomTotal);
         else
            fAtomCharge = fAtomTotal;
         fChargeTotal += fAtomCharge;

         fmt::MemoryWriter
            out;
         out.write("  {:>4} {:>3}   ", iAt+1, m_Atoms[iAt].GetElementName());
         for (size_t i = 0; i <= size_t(iAtMaxL); ++i)
            out.write(pChargeFmtF, AmCharges[i]);
         for (size_t i = size_t(iAtMaxL + 1); i < AmCharges.size(); ++ i)
            out.write(pChargeFmtS, "");
         out << "  ";
         out.write(pChargeFmtF, fAtomTotal);
         out.write(pChargeFmtF, fAtomCharge);

         Log.Write(out.str());
         // continue with next atom.
         itBeg = itEnd;
      }
      Log.WriteLine();
      if (nAtoms != 0) {
         fmt::MemoryWriter
            out;
         out << "  -> Total" << "     ";
         for (size_t i = 0; i < AmCharges.size(); ++ i)
            out.write(pChargeFmtS, "");
         out.write(pChargeFmtF, fElecTotal);
         out.write(pChargeFmtF, fChargeTotal);
         Log.Write(out.str());
         Log.WriteLine();
      }
   }
}

// void FFrame::RunChargeAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms)
// {
//    FAtomSet
//       *pAtoms = pGetAtoms();
//    if (pAtoms == 0) {
//       Log.Write(" No atoms.");
//       return;
//    }
//    FChargeAnalysis
//       ChargeAnalysis(pAtoms, pSelectedAtoms);
//    FDataSetList::iterator
//       itDataSet;
//    for (itDataSet = m_Data.begin(); itDataSet != m_Data.end(); ++ itDataSet) {
//       FOrbital
//          *pOrb = dynamic_cast<FOrbital*>(itDataSet->get());
//       if (!pOrb)
//          continue;
//       pOrb->AddChargeContributions(ChargeAnalysis);
//    }
//    ChargeAnalysis.MakeReport(Log);
// }




FBondOrderAnalysis::FValue const *FBondOrderAnalysis::FindValue(int iAt_, int jAt_) const
{
   FKey
      Key(iAt_, jAt_);
   FMap::const_iterator
      itMap = m_Data.find(Key);
   if (itMap == m_Data.end() || itMap->first != Key)
      return 0;
   else
      return &itMap->second;
}


void FBondOrderAnalysis::MakeBondOrdersForDensity(FRawBasis const *pRawBasis, FMatrixView Rdm, double fOrbOcc, FMemoryStack &Mem)
{
   assert(Rdm.IsSquare() && Rdm.IsSymmetric());

   FStackMatrix
      // 1-particle hole density matrix.
      Hdm(Rdm.nRows, Rdm.nCols, &Mem);
   Assign(Hdm, Rdm, -1.);
   for (size_t i = 0; i < Hdm.nRows; ++ i)
      Hdm(i,i) += fOrbOcc;
   // ^- hm... is what I am doing here really right for the UHF case? Need to compute the stupd correlation matrix element again.

   for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt) {
      if (!IsIncluded(int(iAt)))
         continue;
      double fRenTotal = 0.;
      for (size_t jAt = 0; jAt < m_Atoms.size(); ++ jAt) {
         size_t
            iFnA = pRawBasis->iCenFn(iAt),
            nFnA = pRawBasis->nCenFn(iAt),
            iFnB = pRawBasis->iCenFn(jAt),
            nFnB = pRawBasis->nCenFn(jAt);
         FMatrixView
            RdmAA = Select(Rdm, iFnA, iFnA, nFnA, nFnA),
            RdmAB = Select(Rdm, iFnA, iFnB, nFnA, nFnB),
            HdmAB = Select(Hdm, iFnA, iFnB, nFnA, nFnB),
            RdmBB = Select(Rdm, iFnB, iFnB, nFnB, nFnB);
         double
            f0 = 2. / fOrbOcc, // hmmm...
            fWib = -Dot(RdmAB, HdmAB) * f0,
            fRen = 0;

         FValue
            &val = m_Data[FKey(iAt,jAt)]; // <- will default-construct (set to zero) non-present data sets.
         val.fBondOrders[BONDORDER_WibergIao] += fWib;
#ifdef INCLUDE_RENORM_BO
         { // compute renormalized bond order
            FStackMatrix
               SmhAA(nFnA, nFnA, &Mem),
               SmhBB(nFnB, nFnB, &Mem);

            FSmhOptions
               SmhOpt(1e-6,1e-6);
            Assign(SmhAA, RdmAA);
            CalcSmhMatrix(SmhAA, Mem, SmhOpt);
            Assign(SmhBB, RdmBB);
            CalcSmhMatrix(SmhBB, Mem, SmhOpt);

            FStackMatrix
               NormRdmAB(nFnA, nFnB, &Mem);
            ChainMxm(NormRdmAB, Transpose(SmhAA), RdmAB, Transpose(SmhBB), Mem);

            fRen = Dot(NormRdmAB, NormRdmAB) / f0;
            if (iAt != jAt) {
               val.fBondOrders[BONDORDER_Renorm1Iao] += fRen;
               fRenTotal += fRen;
            }
         }
#else
         IR_SUPPRESS_UNUSED_WARNING(RdmAA);
         IR_SUPPRESS_UNUSED_WARNING(RdmBB);
         IR_SUPPRESS_UNUSED_WARNING(fRen);
         IR_SUPPRESS_UNUSED_WARNING(fRenTotal);
#endif // INCLUDE_RENORM_BO
      }
#ifdef INCLUDE_RENORM_BO
      if (1) {
         FValue
            &val = m_Data[FKey(iAt,iAt)]; // <- valency entry.
         val.fBondOrders[BONDORDER_Renorm1Iao] -= fRenTotal;
      }
#endif // INCLUDE_RENORM_BO
   }
}


FBondOrderAnalysis::FBondOrderAnalysis(ct::FAtomSet *pAtoms, FBasisSetPtr pMinBasis, FMatrixView IbRdmC, FMatrixView IbRdmO, FMemoryStack &Mem, FAtomIdList *pSelectedAtoms)
   : FAnalysisReport(pAtoms, pSelectedAtoms)
{
   FRawBasis const
      *pRawBasis = pMinBasis->pRawBasis.get();
   bool HaveSpin = IbRdmO.fRmsdFromZero() > 1e-8;

   size_t
      nIb = pMinBasis->nFn();
   FStackMatrix
      IbRdmA(nIb, nIb, &Mem);
   if (HaveSpin || true) {
      // make alpha bonds.
      Assign(IbRdmA, IbRdmC, .5);
      Add(IbRdmA, IbRdmO, .5);
      MakeBondOrdersForDensity(pRawBasis, IbRdmA, 1., Mem);

      // make beta bonds.
      Assign(IbRdmA, IbRdmC, .5);
      Add(IbRdmA, IbRdmO, -.5);
      MakeBondOrdersForDensity(pRawBasis, IbRdmA, 1., Mem);
   } else {
      MakeBondOrdersForDensity(pRawBasis, IbRdmC, 2., Mem);
   }

   IR_SUPPRESS_UNUSED_WARNING(IbRdmO);
}

FBondOrderAnalysis::~FBondOrderAnalysis()
{}

struct FBondOrderSortPred {
   bool operator() (FBondOrderAnalysis::FPair const &A, FBondOrderAnalysis::FPair const &B) const {
      return A.second.fBondOrders[0] > B.second.fBondOrders[0];
   }
};

void FBondOrderAnalysis::MakeReport(ct::FLog &Log)
{
   double
      ThrPrint = 0.02;
   {
      Log.Write(" Bond order analysis:       [THRPRINT= {:.2e}]\n", ThrPrint);
//       Log.Write("   CEN ATOM   CEN ATOM     WIBERG.BO     RENORM.BO");
      Log.Write("   CEN ATOM    CEN ATOM   WIBERG.BO"
#ifdef INCLUDE_RENORM_BO
         "    RENORM.BO"
#endif // INCLUDE_RENORM_BO
      );
      char const
         *pChargeFmtF = " {:12.5f}";
//          *pChargeFmtF = " {:13.6f}";
//          *pChargeFmtS = " {:15}";
      FMap::const_iterator
         itBeg, itEnd;
      TArray<FPair>
         BondOrders;
      bool
         EmitNewLine = false;
      for (itBeg = m_Data.begin(); itBeg != m_Data.end(); itBeg = itEnd) {
         BondOrders.clear();
         int iAt = itBeg->first.iAt;

         // find all other atoms this one connects to.
         itEnd = itBeg;
         while (itEnd != m_Data.end() && itEnd->first.iAt == iAt) {
            BondOrders.push_back(FPair(itEnd->first, itEnd->second));
            ++ itEnd;
         }

         // sory entries by first bond order---large ones first.
         std::sort(BondOrders.begin(), BondOrders.end(), FBondOrderSortPred());

         // add a line between different atoms iAt.
         if (EmitNewLine) {
            Log.WriteLine();
            EmitNewLine = false;
         }

//          for (it = itBeg; it != itEnd; ++it) {
//             FValue const
//                &val = it->second;
//             FKey const
//                &key = it->first;
         for (size_t iEntry = 0; iEntry != BondOrders.size(); ++ iEntry) {
            FValue const
               &val = BondOrders[iEntry].second;
            FKey const
               &key = BondOrders[iEntry].first;

            bool Keep = false;
            assert(BONDORDER_Count == sizeof(val.fBondOrders)/sizeof(val.fBondOrders[0]));
            for (int iBo = 0; iBo < BONDORDER_Count; ++ iBo)
               if (std::abs(val.fBondOrders[iBo]) > ThrPrint)
                  Keep = true;
            if (!Keep)
               continue;

            fmt::MemoryWriter
               out;
            if (key.iAt != key.jAt) {
               out.write("  {:>4} {:>3} -  {:>4} {:>3}", key.iAt+1, m_Atoms[key.iAt].GetElementName(), key.jAt+1, m_Atoms[key.jAt].GetElementName());
               assert(BONDORDER_Count == sizeof(val.fBondOrders)/sizeof(val.fBondOrders[0]));
               for (int iBo = 0; iBo < BONDORDER_Count; ++ iBo)
                  out.write(pChargeFmtF, val.fBondOrders[iBo]);
            } else {
               out.write("  {:>4} {:>3}     valency",  key.iAt+1, m_Atoms[key.iAt].GetElementName());
//                out.write("  {:>4} {:>3}    tot.val.",  key.iAt+1, m_Atoms[key.iAt].GetElementName());
//                out.write(pChargeFmtF, -val.fBondOrders[0]);
               for (int iBo = 0; iBo < BONDORDER_Count; ++ iBo)
                  out.write(pChargeFmtF, -val.fBondOrders[iBo]);
            }

            Log.Write(out.str());
            EmitNewLine = true;
         }
      }

      Log.WriteLine();
   }
}







FFragmentAnalysisOptions::FFragmentAnalysisOptions()
{
   GridDesc = "grid{accu: 1e-5}";
//       TfvcDftGridParams("grid{accu: 1e-5; lmin: [15,15,15,15]; nradial: [45,45,45,45]}");
   UseFragmentSpecificGrids = false;
   GroupRemainingAtoms = true;
   SpacePartitionType = SPACEPARTITION_HilbertSpace_Iao;

#ifdef PROG_IBOVIEW
   LoadState();
#endif // PROG_IBOVIEW
}

#ifdef PROG_IBOVIEW
void FFragmentAnalysisOptions::SaveState()
{
   QSettings
      settings;
   settings.setValue("FragmentAnalysis/SpacePartitionType", int(SpacePartitionType));
   settings.setValue("FragmentAnalysis/GridDesc", s2q(GridDesc));
   settings.setValue("FragmentAnalysis/GroupRemainingAtoms", int(GroupRemainingAtoms));
}

void FFragmentAnalysisOptions::LoadState()
{
   QSettings
      settings;
   SpacePartitionType = FSpacePartitionType(settings.value("FragmentAnalysis/SpacePartitionType", 0).toInt());
   GridDesc = settings.value("FragmentAnalysis/GridDesc", QString("grid{accu: 1e-5}")).toString().toStdString();
   GroupRemainingAtoms = bool(settings.value("FragmentAnalysis/GroupRemainingAtoms", 1).toInt());
}
#endif // PROG_IBOVIEW



// represents a partitioning of real space into weighted fragment contributions.
// So far this is only implemented for TFVC -- topological fuzzy Voronoi cells.
struct FFragmentSpacePartitioningTfvc : public FIntrusivePtrDest
{
   explicit FFragmentSpacePartitioningTfvc(ct::FAtomSet *pAtoms, ct::FBasisSetPtr pOrbBasis, ct::FMatrixView Orbs, double const *pOrbOcc, double const *pOrbSpin, ct::FMemoryStack &Mem, FFragmentationDef const &FragmentDefs_, FFragmentAnalysisOptionsCptr pOptions_, ct::FLog *pLog_);

   void EvalFragmentWeights(double *pOut, FMatrixView Grid_, size_t iFragment, FMemoryStack &Mem) const;

   void EvalFragmentOrbitalOverlap(FMatrixView Sij, size_t iFragment, FMemoryStack &Mem);
protected:
   ct::FLog
      *m_pLog;
   FFragmentAnalysisOptionsCptr
      m_pOptions;
   ct::FAtomSet
      *m_pAtoms;
   ct::FDensityEvalContextPtr
      m_pDensityEvalContext;
   FFragmentationDef
      m_FragmentDefs;

//    struct FLinePartition {
//       double Chi; // chi_AB quantity as defined in TFVC paper
//    };
   FHeapMatrix
      m_ChiAB; // (nAt, nAt) matrix of Chi_AB quantities as defined in TFVC paper

   mig::FVoronoiPartitionParamsPtr
      m_pVoronoiPartitionParams;
   mig::FVoronoiPartitionPtr
      m_pVoronoiPartition;
   mig::FDftGridPtr
      // if using all-fragment DFT grid mode, this object will receive the
      // (shared) integration grid for the entire molecule.
      m_pAllFragmentsDftGrid;

   void MakeCellPartitions(FMemoryStack &Mem);

   FDftGridPtr MakeGridForFragment(size_t iFragment);
};
typedef ct::TIntrusivePtr<FFragmentSpacePartitioningTfvc>
   FFragmentSpacePartitioningTfvcPtr;
typedef ct::TIntrusivePtr<FFragmentSpacePartitioningTfvc const>
   FFragmentSpacePartitioningTfvcCptr;


FFragmentSpacePartitioningTfvc::FFragmentSpacePartitioningTfvc(ct::FAtomSet *pAtoms, ct::FBasisSetPtr pOrbBasis, ct::FMatrixView Orbs, double const *pOrbOcc, double const */*pOrbSpin*/, ct::FMemoryStack &Mem, FFragmentationDef const &FragmentDefs_, FFragmentAnalysisOptionsCptr pOptions_, ct::FLog *pLog_)
   : m_pLog(pLog_), m_pOptions(pOptions_), m_pAtoms(pAtoms), m_FragmentDefs(FragmentDefs_)
{
   // takes copies of orbital matrix and occupations.
   m_pDensityEvalContext = new FDensityEvalContext(Orbs, pOrbOcc, pOrbBasis);
   MakeCellPartitions(Mem);

   m_pVoronoiPartitionParams = new mig::FVoronoiPartitionParams(mig::FVoronoiPartitionParams::SMOOTHFN_Becke_k4, mig::FVoronoiPartitionParams::ATOMSIZE_Tfvc);
//    m_pVoronoiPartition = new mig::FVoronoiPartition(AsRawAtoms(*m_pAtoms), m_ChiAB, mig::FVoronoiPartition::SMOOTHFN_Becke_k4, FVoronoiPartition::ATOMSIZE_Tfvc);
   m_pVoronoiPartition = new mig::FVoronoiPartition(AsRawAtoms(*m_pAtoms), *m_pVoronoiPartitionParams, &m_ChiAB[0], 1, m_pAtoms->size());
//       explicit FVoronoiPartition(ct::FRawAtomList const &Atoms_, FVoronoiPartitionParams const &Params_, double const *pChiAB_, size_t iRowSt_, size_t iColSt_);
}


struct FFragmentWeightFnTfvc : public ct::FFragmentWeightFn
{
   explicit FFragmentWeightFnTfvc(size_t iTargetFragment_, FFragmentSpacePartitioningTfvcCptr pTfvcPartitioning_);

   virtual void EvalWeights(double *pOut, FMatrixView Grid_, FMemoryStack &Mem) const; // override
protected:
   // hm... maybe just copy over the relevant voronoi data?
   size_t
      // fragment for which the weights are supposed to be computed.
      m_iTargetFragment;
//    FAtomIdList
//       // list of atoms belonging to the fragment
//       m_iFragmentAtoms;
   FFragmentSpacePartitioningTfvcCptr
      m_pTfvcPartitioning;
};


FFragmentWeightFnTfvc::FFragmentWeightFnTfvc(size_t iTargetFragment_, FFragmentSpacePartitioningTfvcCptr pTfvcPartitioning_)
   : m_iTargetFragment(iTargetFragment_), m_pTfvcPartitioning(pTfvcPartitioning_)
{
}


void FFragmentWeightFnTfvc::EvalWeights(double *pOut, FMatrixView Grid_, FMemoryStack &Mem) const
{
   if (1) {
      return m_pTfvcPartitioning->EvalFragmentWeights(pOut, Grid_, m_iTargetFragment, Mem);
   } else {
      // debuging hack:
      // set fragment weight to 1 at all points. This should result in a numerical approximation
      // of identity matrices for the orbital overlap---individually for each fragment.
      for (size_t i = 0; i != Grid_.nCols; ++ i)
         pOut[i] = 1.;
      return;
   }
}



void FFragmentSpacePartitioningTfvc::EvalFragmentWeights(double *pOut, FMatrixView Grid_, size_t iFragment, FMemoryStack &Mem) const
{
   assert(iFragment < m_FragmentDefs.size());
   FAtomIdList const
      &iFragmentAtoms = m_FragmentDefs[iFragment]->m_iAtoms;
   // for a start, do it lamely. Just go through all grid points one by one.
   for (size_t iGridPt = 0; iGridPt != Grid_.nCols; ++ iGridPt) {
      FVec3d
         vGridPt(Grid_(0,iGridPt), Grid_(1,iGridPt), Grid_(2,iGridPt));
      // ...and add up the individual weights of all atoms in the fragment, computing
      // each of them individually. This likely could be done more efficiently, as
      // in the voronoi weigt function the total weight is also computed each time.
      // TODO: fix this.
      double w = 0;
      for (size_t iiAt = 0; iiAt != iFragmentAtoms.size(); ++iiAt) {
         w += m_pVoronoiPartition->GetAtomWeight(vGridPt, iFragmentAtoms[iiAt], Mem);
      }

      pOut[iGridPt] = w;
   }
}


FDftGridPtr FFragmentSpacePartitioningTfvc::MakeGridForFragment(size_t iFragment)
{
   FDftGridParams
      // fixme: put in minimum accuracy requirements here.
      TfvcDftGridParams(m_pOptions->GridDesc);

   if (m_pOptions->UseFragmentSpecificGrids) {
      // construct a DFT integration grid suitable for integrating
      // the electron density of (only) the target fragment's atoms.

      // make an FAtomSet object for the fragment's target atoms
      FAtomIdList
         iFragmentAtoms = m_FragmentDefs[iFragment]->m_iAtoms;
      std::sort(iFragmentAtoms.begin(), iFragmentAtoms.end());

      FAtomSetPtr
         pFragmentAtoms(new FAtomSet);
      for (size_t iiAt = 0; iiAt < iFragmentAtoms.size(); ++ iiAt) {
         size_t iAt = iFragmentAtoms[iiAt];
         assert(iAt < m_pAtoms->size());
         pFragmentAtoms->AddAtom((*m_pAtoms)[iAt]);
      }

      // construct a DFT integration grid for the target fragment's atoms
      m_pLog->Write("\n Integration grid for '{}'", m_FragmentDefs[iFragment]->m_FragmentName);
      FDftGridPtr
         pDftGrid = FDftGridPtr(new ct::FDftGrid(AsRawAtoms(*pFragmentAtoms), TfvcDftGridParams, m_pLog));
      return pDftGrid;
   } else {
      // construct a DFT integration grid suitable for integrating all fragments
      // at the same time.
      if (m_pAllFragmentsDftGrid.get() == 0) {
         m_pLog->Write("\n Sharing integration grid '{}' across fragments.", m_pOptions->GridDesc);
         m_pAllFragmentsDftGrid = FDftGridPtr(new ct::FDftGrid(AsRawAtoms(*m_pAtoms), TfvcDftGridParams, m_pLog));
      }
      return m_pAllFragmentsDftGrid;
   }
}



void FFragmentSpacePartitioningTfvc::EvalFragmentOrbitalOverlap(FMatrixView Sij, size_t iFragment, FMemoryStack &Mem)
{
   assert(iFragment < m_FragmentDefs.size());

   FFragmentWeightFnTfvc
      // make a proxy object capable of evaluating the current fragment's integration
      // weight at any point in space. This one will be passed to FDensityEvalContext's
      // ComputeFragmentOrbitalOverlap function.
      ThisFragmentsWeightFn(iFragment, FFragmentSpacePartitioningTfvcCptr(this));

   FDftGridPtr
      pDftGrid = MakeGridForFragment(iFragment);

   // reshape the DFT grid data into a linear (4,nGridPt) matrix format.
   size_t
      nGridPt = pDftGrid->Points.size();
   FHeapMatrix
      GridPts(4, nGridPt);
   for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt) {
      FDftGrid::FPoint const &p = pDftGrid->Points[iGridPt];
      GridPts(0,iGridPt) = p.vPos[0];
      GridPts(1,iGridPt) = p.vPos[1];
      GridPts(2,iGridPt) = p.vPos[2];
      GridPts(3,iGridPt) = p.fWeight;
   }

   // calculate orbital overlap data on the grid points.
   Sij.Clear();
   m_pDensityEvalContext->AccFragmentOrbitalOverlap(Sij, GridPts, &ThisFragmentsWeightFn, Mem);
}





// NOTE: TFVC chi computation moved to CtVoronoiPartition.cpp

// static double FindMinimumDensityPoint(double const *pfLineFraction, double const *pfDensity, size_t nPts)
// {
//    if (nPts == 0)
//       throw std::runtime_error("FindMinimumDensityPoint: no input data!");
//
//    // find index of lowest density input point.
//    size_t
//       iMin = ArgMin(pfDensity, pfDensity + nPts);
//    // fit a polynomial through adjacent points.
//    size_t
//       nPtsLeftRight = 4,
//       iFirst = size_t(std::max(ptrdiff_t(iMin) - ptrdiff_t(nPtsLeftRight), ptrdiff_t(0))),
//       iLast = std::min(iMin + nPtsLeftRight, nPts),
//       nPointToFit = iLast - iFirst,
//       nPolyDegree = nPointToFit - 1; // note: 0'th order polynomial (constant) has 1 degree of freedom.
//
//    ct::TPolynomialFit<double, double>
//       PolyFit(&pfDensity[iFirst], 0, &pfLineFraction[iFirst], 0, nPointToFit, nPolyDegree);
//
//    ct::TPolynomialMinimizeParams<double, double>
//       MinimizeParams(PolyFit.tFirst, PolyFit.tLast);
//    MinimizeParams.Print = false;
//    return ct::FindMinimum(PolyFit.p, PolyFit.tMinGuess, MinimizeParams);
// }
//
//
//
// // compute Chi_AB fuzzy voronoi cell shift parameter from given A--B line density profile and line fractions.
// // Line fractions go from (-1. to +1.0).
// static double FindTfvcChiFromDensityProfile(double const *pfDensity, double const *pfLineFraction, size_t nPts)
// {
//    double
//       // find mu for minimum of density
//       muMinRho = FindMinimumDensityPoint(&pfLineFraction[0], &pfDensity[0], nPts);
//    double
//       // closed-form solution for Solve[nu(muMinRho, chi) = 0, {chi}], with given muMinRho
//       Chi = (1 + muMinRho)/(1 - muMinRho);
// //    std::cout << fmt::format("Transformation factor chi_AB = {:16.8f}", Chi) << std::endl;
// //    std::cout << fmt::format("-> nu({},chi) = {}", muMinRho, NuAbSalvador(muMinRho, Chi)) << std::endl;
// //    for (size_t i = 0; i != nPts; ++ i) {
// //       double mui = pfLineFraction[i];
// //       std::cout << fmt::format("   mu = {:16.8f}      rho = {:16.8f}     ->  nu = {:16.8f}", mui, pfDensity[i], NuAbSalvador(mui, Chi)) << std::endl;
// //    }
//    return Chi;
// }



// FIXME: merge moved version in CtVoronoiPartition.h // FVoronoiChiData
void FFragmentSpacePartitioningTfvc::MakeCellPartitions(FMemoryStack &Mem)
{
   size_t
      nAt = m_pAtoms->size();
   m_ChiAB.Reshape(nAt,nAt);
   // set all Chi_AB to 1 --- this corresponds to cell boundaries at the
   // center of the line connecting both atoms.
   for (size_t jAt = 0; jAt != nAt; ++ jAt)
      for (size_t iAt = 0; iAt != nAt; ++ iAt)
         m_ChiAB(iAt,jAt) = 1.;

   // TODO: move OpenMP parallelization from FDensityEvalContext::ComputeDensity
   // to this place, into the outer loop.
   for (size_t iAt = 0; iAt != nAt; ++ iAt) {
      for (size_t jAt = iAt + 1; jAt != nAt; ++ jAt) {
         ct::FAtom const
            &AtA = (*m_pAtoms)[iAt],
            &AtB = (*m_pAtoms)[jAt];
         FVec3d
            vAtA = AtA.vPos,
            vAtB = AtB.vPos,
            vCen = 0.5 * (vAtA + vAtB);
         double
            fDistAB = Length(vAtA - vAtB),
            LinearResolution = 20.; // 20 pts per a.u. for sampling on the connecting line

         bool TreatAsNeighbors = true;

         if (fDistAB < 1e-8) {
            // too close.
            TreatAsNeighbors = false;
            throw std::runtime_error(fmt::format("illegal coincident centers ({} and {}) in FFragmentSpacePartitioningTfvc::MakeCellPartitions", iAt, jAt));
         } else if (fDistAB > 5. * (GetCovalentRadius(AtA.iElement) + GetCovalentRadius(AtB.iElement))) {
            // too far apart.
            TreatAsNeighbors = false;
         } else {
            // check if there is a third atom closer to the center on the line A--B than A or B.
            double
               fHalfDistSq = sqr(.5 * fDistAB);
            for (size_t kAt = 0; kAt != nAt; ++ kAt) {
               if (kAt == iAt || kAt == jAt)
                  continue;
               double
                  fDistSqCC = DistSq(vCen, (*m_pAtoms)[kAt].vPos);
               if (fDistSqCC < fHalfDistSq) {
                  TreatAsNeighbors = false;
                  break;
               }
            }
         }


         if (!TreatAsNeighbors) {
            // atoms are either too close, too far apart, or there is a third atom closer to their center.
            // The latter criterion was defined in the TFVC paper als excluding treating atoms as
            // neighbors.
            m_ChiAB(iAt,jAt) = 1.;
            m_ChiAB(jAt,iAt) = 1.;
         } else {
            size_t
               nGridPt = std::max(size_t(std::min(LinearResolution, LinearResolution * fDistAB)), size_t(1000));
            FStackMatrix
               // grid of points along the connecting line between both atoms
               vGridPt(3, nGridPt, &Mem),
               // total electron density along the line.
               Rho(nGridPt, 1, &Mem),
               // line fractions along the line (-1.0 -> vAtA; 1.0 -> vAtB)
               Mu(nGridPt, 1, &Mem);
            for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt) {
               double
//                      x = 2*(np.arange(nPt) + 0.5)/float(nPt) - 1.
                  mui = 2*(0.5 + double(iGridPt)) / double(nGridPt) - 1.;
                  // ^- this will exclude points lying directly on the atoms to
                  //    avoid problems with derivative singularities.
               FVec3d
//                   vi = (1. - f) * vAtA + f * vAtB;
                  vi = 0.5 * ((1. - mui) * vAtA + (1. + mui) * vAtB);

               // remember both the line fraction and the actual grid point.
               Mu[iGridPt] = mui;
               vGridPt(0,iGridPt) = vi[0];
               vGridPt(1,iGridPt) = vi[1];
               vGridPt(2,iGridPt) = vi[2];
            }
            // compute electron density on the grid points.
            m_pDensityEvalContext->ComputeDensity(&Rho[0], vGridPt, Mem);

//             double
//                ChiAB = FindTfvcChiFromDensityProfile(&Rho[0], &Mu[0], nGridPt);
            // GNAAAA... turns out ECP-ed electrons can *absolutely* have density minima
            // in the non-existing core. grmlhmpf. This is an attempted hack around that.
            // FIXME: do this properly.
            size_t i0 = 0;
            if ((*m_pAtoms)[iAt].iElement > 36) {
//                for (size_t i = 0; i != nGridPt; ++ i)
//                   std::cout << fmt::format("   mu = {:16.8f}  rho = {:16.8f}", Mu[i], Rho[i]) << std::endl;
               // find density maximum in left half of atom-connecting line.
               // The idea is that the inter-nuclear minimum should normally be expected to be
               // between the two density maxima.
               i0 = ArgMax(&Rho[0], &Rho[nGridPt/2]);
               nGridPt -= i0;
//                std::cout << fmt::format("   tfvc i0 = {}  x0 = {:16.8f}", i0, Mu[i0]) << std::endl;
//                std::cout << fmt::format("  WARNING: hacked up density profile of putatively ECP-carrying atom {} to skip minimum at inner core region.\n",AtA.GetAtomLabel(iAt));

            }
            if ((*m_pAtoms)[jAt].iElement > 36) {
               size_t i1 = nGridPt/2 + 1;
               i1 = ArgMax(&Rho[i1], &Rho[nGridPt]) + i1;
               nGridPt = i1 - i0;
//                while (Rho[nGridPt-2] > Rho[nGridPt-1])
//                while (Rho[nGridPt-1] < ThrCoreDen)
//                   nGridPt -= 1;
//                std::cout << fmt::format("  WARNING: hacked up density profile of putatively ECP-carrying atom {} to skip minimum at inner core region.\n",AtB.GetAtomLabel(jAt));
            }
            // FIXME: is there a better way to fix this? maybe search minimum from the inside out?

            double
               ChiAB = mig::FindTfvcChiFromDensityProfile(&Rho[i0], &Mu[i0], nGridPt);
            m_ChiAB(iAt,jAt) = ChiAB;
            m_ChiAB(jAt,iAt) = 1./ChiAB;
//             std::cout << fmt::format("   tfvc chi({},{}) = {:16.8f}", AtA.GetAtomLabel(iAt), AtB.GetAtomLabel(jAt), ChiAB) << std::endl;
         }
      }
   }
}




















// copy the sub-matrix defined by the given rows and column indices into a new continuous matrix, which is
// allocated on Mem.
FMatrixView ExtractMatrix(FMatrixView Parent, TArray<size_t> const &iRows, TArray<size_t> const &iCols, FMemoryStack &Mem)
{
   FMatrixView
      Out = MakeStackMatrix(iRows.size(), iCols.size(), Mem);
   for (size_t iiCol = 0; iiCol != iCols.size(); ++ iiCol) {
      for (size_t iiRow = 0; iiRow != iRows.size(); ++ iiRow) {
         size_t
            iRow = iRows[iiRow],
            iCol = iCols[iiCol];
         assert(iRow < Parent.nRows && iCol < Parent.nCols);
         Out(iiRow, iiCol) = Parent(iRow, iCol);
      }
   }
   return Out;
};



FFragmentAnalysis::~FFragmentAnalysis()
{
}

void FFragmentAnalysis::MakeReport(ct::FLog &Log)
{
   Log.WriteProgramIntro("EFFECTIVE OXIDATION STATE ANALYSIS (EOSA)");
   Log.Write(" Redox charge analysis for {} fragments.\n", m_FragmentData.size());
   Log.Write(m_AnalysisTypeComments);
   Log.Write(m_ElecDistComments[0]);
   if (!m_ElecDistComments[1].empty())
      Log.Write(m_ElecDistComments[1]);
   Log.Write(" Assignment of red(O)x-type charge/spin vs. (P)artial charge/spin:\n");
//    Log.Write(" Fragment         nAtoms      nElecA     nElecB     Charge    Spin");
//    Log.Write(" Fragment         nAtoms   nElecA  nElecB     Charge   Spin");
//    Log.Write("   FRAGMENT          ATOMS  ELEC-A  ELEC-B ->F.CHARGE  F.SPIN");
   Log.Write("   FRAGMENT          ATOMS  ELEC-A  ELEC-B ->O.CHARGE  O.SPIN  :  P.CHARGE   P.SPIN");

   assert(m_FragmentData.size() == m_FragmentDefs.size());
   size_t
      nAtomsTotal = 0;
   for (size_t iFragment = 0; iFragment != m_FragmentData.size(); ++ iFragment) {
      FFragmentDef const
         &FragmentDef = *m_FragmentDefs[iFragment];
      FFragmentData
         &FragmentData = m_FragmentData[iFragment];
//       Log.Write(" {:-32s}  {:4}  {:4}  {:4}  {:+4}", FragmentDef.ShortDesc(), FragmentData.nElecA(), FragmentData.nElecB(), FragmentData.nNuclearCharge - FragmentData.nElec(), FragmentData.iSpin());
//       Log.Write("   {:<17s} {:4}  {:6}  {:6}    {:6}  {:+6}", FragmentDef.m_FragmentName, FragmentDef.m_iAtoms.size(), FragmentData.nElecA(), FragmentData.nElecB(), ptrdiff_t(FragmentData.nNuclearCharge) - ptrdiff_t(FragmentData.nElec()), FragmentData.iSpin());
      nAtomsTotal += FragmentDef.m_iAtoms.size();
      Log.Write("   {:<17s} {:4}  {:6}  {:6}    {:6}  {:+6}     {:10.5f} {:9.5f}",
         FragmentDef.m_FragmentName, FragmentDef.m_iAtoms.size(), FragmentData.nElecA(), FragmentData.nElecB(),
         ptrdiff_t(FragmentData.nNuclearCharge) - ptrdiff_t(FragmentData.nElec()), FragmentData.iSpin(),
         double(FragmentData.nNuclearCharge) - FragmentData.fElec(), FragmentData.fSpin());
   }
   if (1) {
      FFragmentData
         FragmentData;
      for (size_t iFragment = 0; iFragment != m_FragmentData.size(); ++ iFragment) {
         FragmentData += m_FragmentData[iFragment];
      }
      Log.WriteLine();
      Log.Write("   {:<17s} {:4}  {:6}  {:6}    {:6}  {:+6}    :{:10.5f} {:9.5f}",
         "-> Total", nAtomsTotal, FragmentData.nElecA(), FragmentData.nElecB(),
         ptrdiff_t(FragmentData.nNuclearCharge) - ptrdiff_t(FragmentData.nElec()), FragmentData.iSpin(),
         double(FragmentData.nNuclearCharge) - FragmentData.fElec(), FragmentData.fSpin());
   }
   Log.WriteLine();
   Log.WriteTiming("EOSA", (double)tEosaTotal);
   Log.WriteLine();
}


// note: FragRdmA will be destroyed.
static void AddFragmentOccupationEntries(FFragmentAnalysis::FOccPairList &AllOccsA, FMatrixView FragRdmA, size_t iFragment, ct::FMemoryStack &Mem)
{
   assert(FragRdmA.IsSquare());
   size_t
      nFn = FragRdmA.nCols;

   // compute fragment natural occupation numbers
   FStackMatrix
      OccEws(nFn, 1, &Mem);
   Diagonalize_LargeEwFirst(FragRdmA, OccEws.pData, Mem);

   for (size_t iFn = 0; iFn != nFn; ++ iFn) {
      AllOccsA.push_back(FFragmentAnalysis::FOccPair(OccEws[iFn], iFragment));
   }
}


void FFragmentAnalysis::MakeFragmentOccupationsForSpin_IB(FOccPairList &AllOccs, FMatrixView IbRdmA, FBasisSet const *pMinBasis, size_t /*iSpinCase*/, size_t /*nElecA*/, FMemoryStack &Mem)
{
   AllOccs.clear();
   AllOccs.reserve(pMinBasis->nFn());

   for (size_t iFragment = 0; iFragment != m_FragmentData.size(); ++ iFragment) {
      TMemoryLock<double>
         pFreeMe(0, &Mem);
      FFragmentDef
         &Fragment = *m_FragmentDefs[iFragment];
      // contruct RDM in fragment.
      FMatrixView
         FragRdmA = ExtractMatrix(IbRdmA, Fragment.m_iMinBasisFns, Fragment.m_iMinBasisFns, Mem);
//       size_t
//          nFragIb = FragRdmA.nRows;
//       // find fragment natural occupation numbers
//       FStackMatrix
//          OccEws(nFragIb, 1, &Mem);
//       Diagonalize_LargeEwFirst(FragRdmA, OccEws.pData, Mem);
//
//       for (size_t iFragIb = 0; iFragIb != nFragIb; ++ iFragIb)
//          AllOccs.push_back(FOccPair(OccEws[iFragIb], iFragment));
      AddFragmentOccupationEntries(AllOccs, FragRdmA, iFragment, Mem);
   }
}



void FFragmentAnalysis::MakeFragmentOccupations_TFVC(FOccPairList &AllOccsA, FOccPairList &AllOccsB, FFrameWfData *pAoWf, FMemoryStack &Mem)
{
//    FMemoryLogQt
   FMemoryLog
      TempLog;
   TMemoryLock<double>
      pFreeMe1(0, &Mem);
   AllOccsA.clear();
   AllOccsB.clear();

   // next: make complete orbital matrix and orbital occupation list
   TMemoryLock<double>
      pFreeMe(0, &Mem);

   size_t
      nAo = pAoWf->nBf,
      nOrb = pAoWf->COccC.nCols + pAoWf->COccA.nCols + pAoWf->COccB.nCols; // closed, alpha, and beta orbitals
   assert(pAoWf->nBf == pAoWf->pBasis->nFn());
   assert(nAo == pAoWf->COccC.nRows);
   assert(nAo == pAoWf->COccA.nRows);
   assert(nAo == pAoWf->COccB.nRows);

   FStackMatrix
      Orbs(nAo, nOrb, &Mem),
      OrbOcc(nOrb,1, &Mem), // orbital occupation numbers
      OrbSpin(nOrb,1, &Mem);
   // copy over closed-shell orbitals, if existing
   FMatrixView
      COccC = Select(Orbs, 0, 0, nAo, pAoWf->COccC.nCols);
   Assign(COccC, pAoWf->COccC);

   FMatrixView
      COccA = Select(Orbs, 0, pAoWf->COccC.nCols, nAo, pAoWf->COccA.nCols);
   Assign(COccA, pAoWf->COccA);

   FMatrixView
      COccB = Select(Orbs, 0, pAoWf->COccC.nCols + pAoWf->COccA.nCols, nAo, pAoWf->COccB.nCols);
   Assign(COccB, pAoWf->COccB);

   size_t
      nOccC = COccC.nCols,
      nOccA = COccA.nCols,
      nOccB = COccB.nCols;

   {
      for (size_t i = 0; i != nOrb; ++i)
         OrbOcc[i] = 1.; // this works for both alpha and beta orbitals.
      for (size_t i = 0; i != nOccC; ++i)
         OrbOcc[i] = 2.; // overwrite initial ones with occupation numbers of closed-shell orbitals.
      size_t k = 0;
      for (size_t i = 0; i != nOccC; ++i)
         OrbSpin[k++] = 0.;
      for (size_t i = 0; i != nOccA; ++i)
         OrbSpin[k++] = 1.0;
      for (size_t i = 0; i != nOccB; ++i)
         OrbSpin[k++] = -1.0;
      assert(k == nOrb);
   }

   size_t
      nElecA = nOccC + nOccA,
      nElecB = nOccC + nOccB;
   assert(nElecA == pAoWf->nElecA && nElecB == pAoWf->nElecB);

   // instanciate FFragmentSpacePartitioningTfvc object
   // (this should take care of computing the TFVC cell boundaries)
   FFragmentSpacePartitioningTfvcPtr
      pTfvcPart(new FFragmentSpacePartitioningTfvc(pAoWf->pAtoms, pAoWf->pBasis, Orbs, &OrbOcc[0], &OrbSpin[0], Mem, m_FragmentDefs, m_pOptions, &TempLog));

   // go through all fragments, and for each fragment, make and diagonalize the
   // orbital overlap matrix.
   for (size_t iFragment = 0; iFragment != m_FragmentData.size(); ++ iFragment) {
      TMemoryLock<double>
         pFreeMe2(0, &Mem);
//       FFragmentDef
//          &Fragment = *m_FragmentDefs[iFragment];

      // compute orbital overlap (of ALL orbitals) within fragment.
      FStackMatrix
         Sij(nOrb, nOrb, &Mem);
      pTfvcPart->EvalFragmentOrbitalOverlap(Sij, iFragment, Mem);
//       Sij.Print(std::cout, fmt::format("FRAGMENT OVERLAP <i|P_{}|j>", iFragment));


      // from the full orbital overlap matrix, assemble a smaller orbital overlap matrix
      // for the alpha and beta orbitals separately.
      {
         FStackMatrix
            SijAlpha(nElecA, nElecA, &Mem);
         Assign(SijAlpha, Select(Sij, 0, 0, (nOccC + nOccA), (nOccC + nOccA)));

         // diagonalize to find natural occupation numbers and add to target occupancy list.
         AddFragmentOccupationEntries(AllOccsA, SijAlpha, iFragment, Mem);
      }

      {
         FStackMatrix
            SijBeta(nElecB, nElecB, &Mem);
         size_t ib = nOccC + nOccA;
         Assign(Select(SijBeta, 0, 0, nOccC, nOccC), Select(Sij, 0, 0, nOccC, nOccC));
         Assign(Select(SijBeta, nOccC, nOccC, nOccB, nOccB), Select(Sij, ib, ib, nOccB, nOccB));
         Assign(Select(SijBeta, 0, nOccC, nOccC, nOccB), Select(Sij, 0, ib, nOccC, nOccB));
         Assign(Select(SijBeta, nOccC, 0, nOccB, nOccC), Select(Sij, ib, 0, nOccB, nOccC));
         AddFragmentOccupationEntries(AllOccsB, SijBeta, iFragment, Mem);
      }
   }

   TempLog.Flush();
   {
//       std::string s = TempLog.GetText().toStdString();
      std::string s = TempLog.GetText();
      if (!s.empty() && s[s.size()-1] == '\n')
         s.resize(s.size() - 1); // otherwise too many trailing '\n'.
      m_AnalysisTypeComments += s;
   }
}





// iSpinCase: 0 -> Alpha, 1 -> Beta
// nElecA: total number of electrons of the given spin case to distribute across the fragments.
void FFragmentAnalysis::AccRedoxElectronsForSpin(FOccPairList &AllOccs, size_t iSpinCase, size_t nElecA, FMemoryStack &/*Mem*/)
{
   using namespace iv;
   assert(m_FragmentData.size() == m_FragmentDefs.size());
 // assert(iSpinCase < sizeof(FFragmentData)/sizeof(m_FragmentData[0]));
 // ^- hm??
   assert(iSpinCase < 2);

   // now distribute nElec electrons across the fragments, such that the fragments with the largest occupation
   // numbers get them first. To construct a redox charge/spin.
   std::sort(AllOccs.begin(), AllOccs.end(), std::greater<FOccPair>());

   if (nElecA > AllOccs.size()) {
      IvNotify(NOTIFY_Error, IvFmt("Cannot distribute %1 electrons across only %2 spin orbital basis functions!", nElecA, AllOccs.size()));
      return;
   }

   // distribute formal/redox electrons to highest-occupancy fractionally occupied
   // fragment-projected orbitals
   ptrdiff_t
      iLastAssignedElec = -1,
      iLastElec = std::min(ptrdiff_t(nElecA), ptrdiff_t(AllOccs.size()));

   for (ptrdiff_t iElecA = 0; iElecA < iLastElec; ++ iElecA) {
      FOccPair
         Pair = AllOccs[iElecA];
      FFragmentData
         &FragmentData = m_FragmentData[Pair.second];
      FragmentData.nElecAb[iSpinCase] += 1;
      iLastAssignedElec = ptrdiff_t(iElecA);
   }

   // also distribute fractional electrons across fragments --- to compute (non-redox) partial charge
   for (size_t iOccEntry = 0; iOccEntry != AllOccs.size(); ++ iOccEntry) {
      FOccPair
         Pair = AllOccs[iOccEntry];
      FFragmentData
         &FragmentData = m_FragmentData[Pair.second];
      FragmentData.fElecAb[iSpinCase] += Pair.first;
   }

   if (iLastAssignedElec == -1) {
      m_ElecDistComments[iSpinCase] = " No electrons to distribute.\n";
   } else if (iLastAssignedElec >= iLastElec)  {
      m_ElecDistComments[iSpinCase] = " Exactly as many electron as basis functions.\n";
   } else {
      fmt::MemoryWriter
         out;
      ptrdiff_t
         nShow = 7;
      assert(iLastAssignedElec >= 0 && iLastAssignedElec < ptrdiff_t(AllOccs.size()));
      double
         fOccGap = AllOccs[iLastAssignedElec].first - AllOccs[iLastAssignedElec+1].first;
      // ^- boundary cases (no orbitals occupied/all orbitals occupied) should have been handled
      //    by outer 'if's.
      char const
         *pSpinDesc = (iSpinCase==0)? "Alpha" : "Beta";
      out.write(" Distributed {} {} electrons into {} orbitals.\n", iLastElec, pSpinDesc, AllOccs.size());
      out.write(" Redox occupation assignments for spin {}:\n", pSpinDesc);
      out << "\n     NAT.ORB.     FRAGMENT    OCC.NUM     REDOX STATUS\n";
      ptrdiff_t
         iFirstToShow = std::max(ptrdiff_t(iLastAssignedElec-nShow), ptrdiff_t(0)),
         iLastToShow = std::min(ptrdiff_t(iLastAssignedElec+nShow+1), ptrdiff_t(AllOccs.size()));
//       iFirstToShow = 0; // FIXME: remove this.
//       out.write("!#_ OVERWRITING iFirstToShow WITH 0.\n");
      for (ptrdiff_t iElec = iFirstToShow; iElec != iLastToShow; ++ iElec) {
//       for (ptrdiff_t iElec = 0; size_t(iElec) != AllOccs.size(); ++ iElec) {
         out.write("    {:5}   ->  {:7}   {:12.5f}     {}\n", 1+iElec, AllOccs[iElec].second, AllOccs[iElec].first, (iElec <= iLastAssignedElec)? "=> occupied." : "=> empty." );
      }
      out << "\n";
//       between highest formally occupied and lowest formally unoccupied
      out.write(" Occupation gap for redox assignement:  {:12.5f} e-\n", fOccGap);
      // calculate assignment reliability indicator according to 10.1021/ct501088v, eq. 5
      // note: the range of values obtainable is from 50% to 100%
      // (i.e., "50%" denotes the worst case scenario)
      double R = 100. * std::min(1., std::max(0., 0.5 + fOccGap));
      out.write(" Redox assignment reliability (R):      {:9.2f} %\n", R);
//       out << "\n";

      m_ElecDistComments[iSpinCase] = out.str();
   }

}


// FFragmentAnalysis::FFragmentAnalysis(ct::FAtomSet *pAtoms, ct::FBasisSetPtr pMinBasis, ct::FMatrixView IbRdmC, ct::FMatrixView IbRdmO, size_t nElec, ptrdiff_t iSpin, ct::FMemoryStack &Mem, FFragmentationDef const &FragmentDefs_)
FFragmentAnalysis::FFragmentAnalysis(ct::FAtomSet *pAtoms, FFrameWfData *pIbWf, FFrameWfData *pAoWf, ct::FMemoryStack &Mem, FFragmentationDef const &FragmentDefs_, FFragmentAnalysisOptionsCptr pOptions_)
   : FAnalysisReport(pIbWf->pAtoms, 0), m_FragmentDefs(FragmentDefs_), m_pOptions(pOptions_)
{
   using namespace iv;
   assert(pIbWf == 0 || pIbWf->pAtoms == pAtoms);
   assert(pAoWf == 0 || pAoWf->pAtoms == pAtoms);
   size_t
      nFragments = m_FragmentDefs.size();
   m_FragmentData.resize(nFragments);

   // count the nuclear charges. FIXME: Handling of ECPs not correct here. Should work for def2 sets,
   // but not VnZ-PP sets if including elements between 21 and 36.
   for (size_t iFragment = 0; iFragment != m_FragmentData.size(); ++ iFragment) {
      FFragmentDef const
         &FragmentDef = *m_FragmentDefs[iFragment];
      FFragmentData
         &FragmentData = m_FragmentData[iFragment];
      for (size_t ii = 0; ii < FragmentDef.m_iAtoms.size(); ++ ii) {
         size_t
            iElement = (*pAtoms)[FragmentDef.m_iAtoms[ii]].iElement,
            iNuclearCharge = iElement;
         if (iElement >= 37)
            iNuclearCharge -= iLikelyEcpCharge(iElement);
         FragmentData.nNuclearCharge += iNuclearCharge;
      }
   }

   {
      m_HaveSpin = pIbWf->RdmO.fRmsdFromZero() > 1e-8;

      size_t
         nElecA = pIbWf->nElecA,
         nElecB = pIbWf->nElecB;
      FOccPairList
         AllOccsA, AllOccsB;
      if (m_pOptions->SpacePartitionType == FFragmentAnalysisOptions::SPACEPARTITION_HilbertSpace_Iao) {
         m_AnalysisTypeComments = " Fragment occupations computed via IAO partitioning of density matrix\n";
         size_t
            nIb = pIbWf->nBf;
         FStackMatrix
            IbRdmA(nIb, nIb, &Mem);
         if (m_HaveSpin) {
            // distribute alpha electrons.
            Assign(IbRdmA, pIbWf->RdmC, .5);
            Add(IbRdmA, pIbWf->RdmO, .5);
            MakeFragmentOccupationsForSpin_IB(AllOccsA, IbRdmA, &*pIbWf->pBasis, 0, nElecA, Mem);

            // distribute beta electrons.
            Assign(IbRdmA, pIbWf->RdmC, .5);
            Add(IbRdmA, pIbWf->RdmO, -.5);
            MakeFragmentOccupationsForSpin_IB(AllOccsB, IbRdmA, &*pIbWf->pBasis, 1, nElecB, Mem);
         } else {
            assert(nElecA == nElecB);
            Assign(IbRdmA, pIbWf->RdmC, .5);
            MakeFragmentOccupationsForSpin_IB(AllOccsA, IbRdmA, &*pIbWf->pBasis, 0, (nElecA+nElecB)/2, Mem);
         }
      } else if (m_pOptions->SpacePartitionType == FFragmentAnalysisOptions::SPACEPARTITION_RealSpace_TFVC) {
         m_AnalysisTypeComments = " Fragment occupations computed via orbital overlap &lt;i|P<sub>F</sub>|j&gt; in real space\n";
         m_AnalysisTypeComments += " Fragment density partitioning: TFVC (Topological Fuzzy Voronoi Cells)\n";

         MakeFragmentOccupations_TFVC(AllOccsA, AllOccsB, pAoWf, Mem);
      } else {
         IvNotify(NOTIFY_Error, "Encountered unrecognized type of space partitioning in fragment redox charge analysis.");
      }



      AccRedoxElectronsForSpin(AllOccsA, 0, nElecA, Mem);
      if (m_HaveSpin) {
         AccRedoxElectronsForSpin(AllOccsB, 1, nElecB, Mem);
      } else {
         AllOccsB = AllOccsA;
         for (size_t iFragment = 0; iFragment != m_FragmentData.size(); ++ iFragment)
            m_FragmentData[iFragment].ReplicateFromAlpha();
      }
   }

}


typedef std::set<int>
   FAtomIdSet;
void MakeMopUpGroupsForFragmentation(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis, FAtomIdSet const &iAssignedAtoms, bool GroupRemainingAtoms)
{
   // unless given group list is complete, create a mop up group for the remaining atoms not yet assigned anywhere.
   if (iAssignedAtoms.size() != pAtoms->size()) {
      if (GroupRemainingAtoms) {
         FAtomIdList
            FragmentList;
         for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt) {
            // atom already assigned to a group?
            if (iAssignedAtoms.find(iAt) != iAssignedAtoms.end())
               continue;
            FragmentList.push_back(iAt);
         }
         Fragmentation.push_back(new FFragmentDef(FragmentList, pAtoms, pMinBasis, fmt::format("Other Atoms")));
      } else {
         // or, alternatively, make individual groups for all remaining atoms.
         for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt) {
            // atom already assigned to a group?
            if (iAssignedAtoms.find(iAt) != iAssignedAtoms.end())
               continue;
            FAtomIdList
               FragmentList;
            FragmentList.push_back(iAt);
            Fragmentation.push_back(new FFragmentDef(FragmentList, pAtoms, pMinBasis, fmt::format((*pAtoms)[iAt].GetAtomLabel(iAt))));
         }
      }
   }
}


void MakeFragmentationForIndividualAtoms(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis)
{
   if (pAtoms == 0)
      return;
   Fragmentation.clear();
   for (size_t iAt = 0; iAt != pAtoms->size(); ++ iAt) {
      FAtomIdList
         FragmentList;
      FragmentList.push_back(iAt);
      Fragmentation.push_back(new FFragmentDef(FragmentList, pAtoms, pMinBasis, fmt::format("{}", (*pAtoms)[iAt].GetAtomLabel(iAt))));
   }
}

void MakeFragmentationViaSelection(FFragmentationDef &Fragmentation, FAtomSet const *pAtoms, FBasisSet const *pMinBasis, FAtomIdList const &iSelectedAtoms, bool GroupRemainingAtoms)
{
   if (pAtoms == 0)
      return;
   Fragmentation.clear();
   std::string
      SelectedGroupName = "Selected Atoms";
   if (iSelectedAtoms.size() == 1) {
      int iAt = iSelectedAtoms[0];
      SelectedGroupName = (*pAtoms)[iAt].GetAtomLabel(iAt,0);
   } else if (iSelectedAtoms.size() == 2) {
      int iAt = iSelectedAtoms[0];
      int jAt = iSelectedAtoms[1];
      SelectedGroupName = (*pAtoms)[iAt].GetAtomLabel(iAt,0) + ", " + (*pAtoms)[jAt].GetAtomLabel(iAt,0);
   }
   Fragmentation.push_back(new FFragmentDef(iSelectedAtoms, pAtoms, pMinBasis, SelectedGroupName));

   FAtomIdList
      iOtherAtoms;
   FAtomIdSet
      // atoms assigned by the current groups until now.
      iAssignedAtoms;
   for (size_t iiAt = 0; iiAt != iSelectedAtoms.size(); ++iiAt)
      iAssignedAtoms.insert(iSelectedAtoms[iiAt]);

   MakeMopUpGroupsForFragmentation(Fragmentation, pAtoms, pMinBasis, iAssignedAtoms, GroupRemainingAtoms);
}



FFragmentDef::FFragmentDef(FAtomIdList const &iAtoms_, FAtomSet const *pAtoms_, FBasisSet const *pMinBasis_, std::string const &Name_)
   : m_FragmentName(Name_), m_iAtoms(iAtoms_), m_pAtoms(pAtoms_), m_pMinBasis(pMinBasis_)
{
   // sort atom fragments.
   std::sort(m_iAtoms.begin(), m_iAtoms.end());
   for (size_t iiAt = 0; iiAt + 1 < m_iAtoms.size(); ++ iiAt)
      if (m_iAtoms[iiAt] == m_iAtoms[iiAt+1])
         throw std::runtime_error(fmt::format("Fragment definition contains multiple instances of atom {}.", m_iAtoms[iiAt]+1));

   FRawBasis
      *pRawBasis = &*pMinBasis_->pRawBasis;
   // make a list of all minimal-basis basis-functions associated with atoms in the fragment.
   for (size_t iiAt = 0; iiAt != m_iAtoms.size(); ++ iiAt) {
      size_t iAt = size_t(m_iAtoms[iiAt]);
      for (size_t iFn = pRawBasis->iCenFn(iAt); iFn != pRawBasis->iCenFn(iAt+1); ++ iFn)
         m_iMinBasisFns.push_back(iFn);
   }
}

} // namespace ct
