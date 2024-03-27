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
#include <iostream>
#include <cmath>
#include <algorithm>
#include <set>
#include <fstream> // for orbital basis export

#include <QTextStream>
#include <QStringList>
#include <QFileInfo>
#include <QBrush>
#include <QFont>
#include <QHeaderView>
#include <QSize>
#include <QColor>
#include <QSettings>
#include <QClipboard>
#include <QApplication>
// #include <QRegularExpression>
#include <QRegExp>

#include "IvDocument.h"
#include "IvFindOrbitalsForm.h"
#include "IvShowTextForm.h"
#include "IvScript.h"
#include "IvIrc.h"
#include "IvOrbitalFile.h"
#include "IvTables.h"
#include "IvSettings.h"
#include "IvFileConvert.h"
#include "IvVolumeDataSet.h"
#include "IvOrbital.h"

#include "format.h"
#include "CxIo.h"
#include "CxColor.h"
#include "CxPhysicalUnits.h"
#include "CxOpenMpProxy.h"
#include "CxXyzFrameIo.h"

#include "IrAmrr.h"

#include "CtMatrix.h"
// #include "CtBasisLibrary.h"
#include "CtInt1e.h"
#include "CtRhf.h"
using namespace ct;




FAtomOptions &FDocument::AtomOptions(int iAt) {
//    if (iAt >= (int)pAtoms->size()) {
//       std::cout << fmt::format("WARNING: SetAtomFlag: Atom {} does not (yet?) exist. Flags anyway stored for future reference.", iAt) << std::endl;
//    }
   return m_AtomOptions[iAt];
}

uint &FDocument::AtomFlags(int iAt) {
   return AtomOptions(iAt).Flags;
}

bool FDocument::IsAtomHidden(int iAt)
{
   return 0 != (AtomFlags(iAt) & ATOM_Hidden);
}

bool FDocument::IsAtomSelected(int iAt)
{
   return 0 != (AtomFlags(iAt) & ATOM_Selected);
}

bool FDocument::IsAtomExcludedFromAlignment(int iAt)
{
   return 0 != (AtomFlags(iAt) & ATOM_NoAlign);
}


void FDocument::ClearAtomFlags(bool EmitUpdate)
{
   m_AtomOptions.clear();
   m_SelectionSequenceId = 0;
//    m_AtomGroups.clear();
   if (EmitUpdate)
      emit SelectionChanged();
}

struct FSortSelectedAtomsPred
{
   explicit FSortSelectedAtomsPred(FDocument *pDocument_) : m_pDocument(pDocument_) {};

   bool operator() (int iAt, int jAt) const {
      return m_pDocument->AtomOptions(iAt).iSelectionSequenceId < m_pDocument->AtomOptions(jAt).iSelectionSequenceId;
   }

   FDocument mutable
      *m_pDocument;
};

FAtomIdList FDocument::GetSelectedAtoms(bool SortedBySelectionSequence)
{
   FFrame
      *pFrame = GetCurrentFrame();
   if (!pFrame) {
      IV_NOTIFY(NOTIFY_Error, "No current frame for selecting atoms");
      return FAtomIdList();
   }
   ct::FAtomSet
      *pAtoms = pFrame->pGetAtoms();
   if (!pAtoms) {
      IV_NOTIFY(NOTIFY_Error, "Frame has no associated geometry data");
      return FAtomIdList();
   }

   FAtomIdList
      iSelectedAtoms;
   FAtomOptionMap::iterator
      itAtomOptions;
   for (itAtomOptions = m_AtomOptions.begin(); itAtomOptions != m_AtomOptions.end(); ++ itAtomOptions)
      if (bool(itAtomOptions->second.Flags & ATOM_Selected) && itAtomOptions->first >= 0 && (size_t)itAtomOptions->first < pAtoms->size() )
         iSelectedAtoms.push_back(itAtomOptions->first);
   if (SortedBySelectionSequence) {
      FSortSelectedAtomsPred SortBySeqPred(this);
      std::sort(iSelectedAtoms.begin(), iSelectedAtoms.end(), SortBySeqPred);
   }
   return iSelectedAtoms;
}

QString FDocument::AtomLabel(int iAt) const
{
   if (!GetCurrentFrame())
      return "UNKN";
   FAtomSet const
      *pAtoms = GetCurrentFrame()->pGetAtoms();
   if (!pAtoms || iAt >= (int)pAtoms->size())
      return "ERR";
   return QString("%1 %2").arg(s2q((*pAtoms)[iAt].GetElementName())).arg(1+iAt);
}

// QString FDocument::FmtDistance(double fLengthInAu) const
// {
//    return QString("%1 pm").arg(fLengthInAu * ct::ToAng * 100., 0, 'f', 2);
// }



void EmitEndl(QTextStream &out) {
   int fieldWidth = out.fieldWidth();
   out.setFieldWidth(0);
   out << "\n";
   out.setFieldWidth(fieldWidth);
}

QString FDocument::GetAtomLabel(int iAt)
{
   FAtomSet
      *pAtoms = GetCurrentFrame()->pGetAtoms();
   QString
      s;
   QTextStream
      out(&s);
   if (iAt < 0 || (size_t)iAt >= pAtoms->size()) {
      IV_NOTIFY1(NOTIFY_Error, "No atom %1 exists.", iAt+1);
      return "???";
   }
   FAtom
      &At = (*pAtoms)[iAt];
   out.setFieldWidth(2);
   out << s2q(At.GetElementName()).toUpper();
   out.setFieldWidth(3);
   out << (1+iAt);
   return s;
}



void FDocument::CalcChargesOnSelectedAtoms()
{
   FFrame
      *pFrame = GetCurrentFrame();
   if (!pFrame)
      return IV_NOTIFY(NOTIFY_Error, "No current frame to search orbitals");
   FAtomIdList
      iSelectedAtoms = GetSelectedAtoms();
   FMemoryLogQt
      Log;
   FAtomSet
      *pAtoms = GetCurrentFrame()->pGetAtoms();
   if (!pAtoms)
      return IV_NOTIFY(NOTIFY_Error, "No atoms?");
//    IvEmit("\n!Search Orbitals on Atoms: [THRCHG=%1]", ThrPrint);
   QString
      sAtomList = "Charges on Atoms:";
   for (size_t iiAt = 0; iiAt < iSelectedAtoms.size(); ++ iiAt)
      sAtomList.append(QString(" %1").arg(iSelectedAtoms[iiAt]+1));

   pFrame->RunChargeAnalysis(Log, &iSelectedAtoms);

   FShowTextForm
      ShowTextForm("<pre>" + Log.GetText() + "</pre>", sAtomList, GetLogFileName("charges"), QString(), g_pMainWindow);
   ShowTextForm.exec();
}


QString FDocument::GetLogFileName(QString LogType)
{
   QString Ext;
   if (LogType.size() != 0)
      Ext = "." + LogType + "-log.txt";
   else
      Ext = ".log.txt";
   QString s = ReplaceExt(GetCurrentInputFileName(), Ext);
//    s.replace(":","-"); // for sub-frames
   for (int i = 0; i < s.size() - 1; ++ i) {
      if (s[i] == ':' && s[i+1] != '\\') // don't replace c:\ ... or similar.
         s[i] = '-';
   }
   return s;
}


// hm, this really shouldn't be here (but rather in IvMain or IvView), I think.
void FDocument::FindOrbitalsOnSelectedAtoms()
{
   FFrame
      *pFrame = GetCurrentFrame();
   if (!pFrame)
      return IV_NOTIFY(NOTIFY_Error, "No current frame to search orbitals");

   double
      ThrPrint = 0.02;
   FAtomIdList iSelectedAtoms = GetSelectedAtoms();

   FFoundOrbitalList
      FoundOrbitals;

//    IvEmit("\n!Search Orbitals on Atoms: [THRCHG=%1]", ThrPrint);
   QString
      sAtomList = "Find Orbitals on Atoms:";
   for (size_t iiAt = 0; iiAt < iSelectedAtoms.size(); ++ iiAt)
      sAtomList.append(QString(" %1").arg(iSelectedAtoms[iiAt]+1));
   sAtomList.append(QString(" [THRCHG=%1]").arg(ThrPrint));


   FDataSetList::iterator
      itDataSet;
   for (itDataSet = pFrame->m_Data.begin(); itDataSet != pFrame->m_Data.end(); ++ itDataSet) {
      FOrbital
         *pOrb = dynamic_cast<FOrbital*>(itDataSet->get());
      if (!pOrb)
         continue;
      TArray<double>
         IaoCharges = pOrb->MakeIaoCharges(false, 2.0);
      bool
         KeepThis = false;
      double
         fChgOnSelected = 0.;
      for (size_t iiAt = 0; iiAt < iSelectedAtoms.size(); ++ iiAt) {
         size_t iAt = iSelectedAtoms[iiAt];
         if (iAt >= IaoCharges.size()) // not supposed to happen.
            break;
         if (IaoCharges[iAt] > ThrPrint)
            KeepThis = true;
         fChgOnSelected += IaoCharges[iAt];
      }
      if (KeepThis) {
         FoundOrbitals.push_back(FFoundOrbital(fChgOnSelected, itDataSet - pFrame->m_Data.begin(), pOrb));
      }
   }


   FFoundOrbitalModel
      Model(FoundOrbitals, "Sel. Chg.", this);
   FFindOrbitalsForm
      Dialog(&Model, sAtomList, g_pMainWindow);
   Dialog.exec();
}

// this also shouldn't be here.
void FDocument::FindReactingOrbitals()
{
   if (!GetFrameData(0) || !GetCurrentFrameData())
      return;
   FDataSetList
      &DataList = *GetFrameData(0);

   FFoundOrbitalList
      FoundOrbitals;
   for (uint iDataSet = 0; iDataSet < DataList.size(); ++ iDataSet)
   {
      FDataSetPtr
         pActiveData = DataList[iDataSet];
      FOrbital
         *pOrbital = dynamic_cast<FOrbital*>(pActiveData.get());
      if (!pOrbital)
         continue;
      TArray<float>
         CurveData;
      if (!MakeIboChangeCurve(CurveData, iDataSet, this))
         // failed to make the delta-frame data set.
         continue;

      double fMaxChange = 0.;
      for (size_t iFrame = 0; iFrame != CurveData.size(); ++ iFrame) {
         fMaxChange = std::max(fMaxChange, double(CurveData[iFrame]));
      }
//          PrintChargeArray(std::cout, "set curve", CurveData);

      FoundOrbitals.push_back(FFoundOrbital(fMaxChange, iDataSet, pOrbital));
   }

   SortFoundOrbitalList(FoundOrbitals, true);
   FFoundOrbitalModel
      Model(FoundOrbitals, "Max. Change", this);
   FFindOrbitalsForm
      Dialog(&Model, IvFmt("Orbitals by Maximum Change along %1 Frames", GetNumFrames()), g_pMainWindow);
   Dialog.exec();
}












FFrame::FFrame(QString Desc_, FDocument *pDocument_)
   : m_pDocument(pDocument_)
{
   m_InputFileName = Desc_;
   m_pFrameLog = new FMemoryLogQt(this);

   // FIXME: should this be here?
   m_pFreeObjects.reset(new FFreeObjectSet("Free Objects", FAtomSetPtr(0), m_pDocument));
   m_pFreeObjects->Active = true;
}

IFrame::IFrame(QString Desc_, FDocument *pDocument_) : FFrame(Desc_, pDocument_) {
   setObjectName(QString("IFrame(%1)").arg(m_InputFileName));
}

int IFrame::get_num_atoms() const
{
   ct::FAtomSet const *pAtoms = pGetAtoms();
   return pAtoms? int(pAtoms->size()) : 0;
}


// void IFrame::add_objects(QList<FFreeObject*> const &fol, QString Mode)
// {
//    if (!(Mode == "link" || Mode == "clone")) {
//       IV_NOTIFY(NOTIFY_Warning, QString("IFrame::add_objects(list, mode): mode should be 'link' or 'clone', but was '%1'").arg(Mode));
//    }
//
//    for (FFreeObject const *pObj : fol) {
//       if (Mode == "link")
//          m_pFreeObjects->AddLink(const_cast<FFreeObject*>(pObj));
//       else
//          m_pFreeObjects->AddClone(pObj);
//    }
// }
void IFrame::add_objects(QVariantList fol, QString Mode)
{
   int Verbosity = 1;
   if (Verbosity >= 1)
      IvEmit("IFrame::add_objects(list, '%1') { // list.size = %2", Mode, fol.size());
   if (!(Mode == "link" || Mode == "clone")) {
      IV_NOTIFY(NOTIFY_Warning, QString("IFrame::add_objects(list, mode): mode should be 'link' or 'clone', but was '%1'").arg(Mode));
   }

   size_t iEntry = 0;
   for (auto &&var : fol) {
//       IvEmit(QString("! add_objects(): [%1] = '%2'").arg(iEntry).arg(var.toString()));
      QObject
         *pObjectBase = var.value<QObject*>();
      FFreeObject
         *pFreeObject = !bool(pObjectBase)? 0 : qobject_cast<FFreeObject*>(pObjectBase);
//       if (pObjectBase)
//          IvEmit(QString("! add_objects(): [%1] class = '%2' toString() = '%3'").arg(iEntry).arg(pObjectBase->metaObject()->className(), pFreeObject? pFreeObject->toString() : QString("null")));
      // @canConvert:
      // """A QVariant containing a pointer to a type derived from QObject will
      // also return true for this function if a qobject_cast to the type
      // described by targetTypeId would succeed. Note that this only works for
      // QObject subclasses which use the Q_OBJECT macro."""
      // ...hm... doesn't work :(
//       if (!var.canConvert<FFreeObject*>()) {
      if (!bool(pFreeObject)) {
         IV_NOTIFY(NOTIFY_Warning, QString("IFrame::add_objects(list, mode): input object of type '%1' cannot be converted to FFreeObject").arg(var.typeName()));
      } else {
//          FFreeObject *pObj = var.value<FFreeObject*>();
         if (Verbosity >= 1) {
            IvEmit("   [%1] class = '%2' -> add '%3'", iEntry, pObjectBase->metaObject()->className(), pFreeObject? pFreeObject->toString() : QString("null"));
         }
         assert(pFreeObject != 0);
         if (Mode == "link")
            m_pFreeObjects->AddLink(const_cast<FFreeObject*>(pFreeObject));
         else
            m_pFreeObjects->AddClone(pFreeObject);
      }
      iEntry += 1;
   }
   if (Verbosity >= 1)
      IvEmit("} // now: m_pFreeObjects->size() == %1;", m_pFreeObjects->size());
}

void IFrame::clear_objects()
{
   m_pFreeObjects->clear();
}



// collect basis function shells on the orbital basis of this frame, find unique ones,
// and export them into library format (libmol).
void IFrame::export_orbital_basis_set(QString const &FileName, QString const &BasisName)
{
   ct::FAtomSet const *pAtoms = pGetAtoms();
   if (pAtoms == 0) {
      IV_NOTIFY(NOTIFY_Warning, "IFrame::export_orbital_basis_set: current frame has no atoms to export basis of.");
      return;
   }
   if (!HaveOrbitals()) {
      IV_NOTIFY(NOTIFY_Warning, "IFrame::export_orbital_basis_set: current frame has no orbitals to export basis of.");
      return;
   }

   // grab first orbital we can find (note: index starts at 1), and take
   // it's orbital basis.
   FOrbital
      *pOrb = this->pGetOrbital(1);
   if (pOrb == 0) {
      IV_NOTIFY(NOTIFY_Error, "IFrame::export_orbital_basis_set: internal error: have orbitals, but cannot access them.");
      return;
   }

   std::ofstream
      out(q2s(FileName).c_str(), std::ios_base::out | std::ios_base::binary);
   // ^-- should I set it to binary mode? in case someone runs it on a mac...
   if (!out.good()) {
      IV_NOTIFY1(NOTIFY_Error, "IFrame::export_orbital_basis_set: failed to open file '%1' for writing", FileName);
   }
   FBasisSet const
      *pOrbBasis = pOrb->pBasisSet.get();
   assert_rt(pOrbBasis != 0);
   std::string
      Comment(fmt::format("Exported by IboView for {}", q2s(GetBaseInputFileName())));
   try {
      pOrbBasis->ExportToLibmol(out, q2s(BasisName), Comment, *pAtoms);
   } catch (std::runtime_error &e) {
      IV_NOTIFY1(NOTIFY_Error, "IFrame::export_orbital_basis_set: export failed. Message: '%1'", s2q(e.what()));
   }
   out.flush();
   if (!out.good()) {
      IV_NOTIFY1(NOTIFY_Error, "IFrame::export_orbital_basis_set: I/O error in writing '%1'", FileName);
   } else {
      IV_NOTIFY1(NOTIFY_Information, "IFrame::export_orbital_basis_set: Basis set data successfully written to '%1'.", FileName);
   }
}



bool FFrame::HaveOrbitals() const
{
   FDataSetList::const_iterator itDataSet;
   _for_each(itDataSet, m_Data) {
      if (dynamic_cast<FOrbital const*>(&**itDataSet) != 0)
         return true;
   }
   return false;
}

void FFrame::ClearOrbitals()
{
   // make a new list of data sets without the electronic structure things.
   FDataSetList
      New;
   FDataSetList::const_iterator itDataSet;
   _for_each(itDataSet, m_Data) {
      if (!(*itDataSet)->DependsOnWaveFunction()) {
         New.push_back(*itDataSet);
      }
   }
   m_Data.swap(New);
//    m_pOrbBasis = 0;
//    m_pMinBasis = 0;
//    m_CIaoData.clear();
//    m_CIb = FMatrixView(0,0,0);
}

void FDocument::ClearOrbitals(bool EmitSignals)
{
   if (EmitSignals)
      BeginTotalReset();
      //emit layoutAboutToBeChanged();
   for (size_t iFrame = 0; iFrame < m_Frames.size(); ++ iFrame)
      GetFrame(iFrame)->ClearOrbitals();
   if (EmitSignals)
      EndTotalReset();
   //emit layoutChanged();
}

// std::string BaseName(std::string const &FileName) {
//    // get the base file name as description, without the path.
//    size_t
//       iPathSept = FileName.find_last_of('/');
//    std::string
//       s;
//    if ( iPathSept != std::string::npos )
//       s = FileName.substr(1+iPathSept);
//    else
//       s = FileName;
//    return s;
// }

QString FFrame::GetBaseInputFileName() const
{
   return RemoveExt(m_InputFileName);
}

QString FFrame::GetFullInputFileName() const
{
   return m_InputFileName;
}


FGeometry *FFrame::pGetGeometry()
{
//    if (m_Data.empty())
//       return 0;
//    return dynamic_cast<FGeometry*>(m_Data[0].get());
   for (size_t i = 0; i < m_Data.size(); ++ i) {
      // find and return the first FGeometry-class data set.
      // Normally that should be data set [0].
      FGeometry
         *pGeometry = dynamic_cast<FGeometry*>(m_Data[i].get());
      if (pGeometry != 0)
         return pGeometry;
   }
   return 0;
}

ct::FAtomSet *FFrame::pGetAtoms()
{
   FGeometry
      *pGeometry = pGetGeometry();
   if (pGeometry == 0)
      return 0;
   return pGeometry->pAtoms.get();
}

double FFrame::GetEnergy() const
{
   ct::FAtomSet const
      *pAtoms = pGetAtoms();
   if (pAtoms)
      return pAtoms->GetLastEnergy();
   return 0.; // no atoms -> no energy.
}

double FFrame::GetGradient() const
{
   ct::FAtomSet const
      *pAtoms = pGetAtoms();
   if (pAtoms) {
      double
         fGrad = 0.;
      for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt)
         fGrad += ct::LengthSq((*pAtoms)[iAt].vGrad);
      fGrad = std::sqrt(fGrad);
      return fGrad;
   }
   return 0.;
}


FOrbital *FFrame::pGetOrbital(int iMo)
{
   if (iMo <= 0 || size_t(iMo) >= m_Data.size())
      return 0;
   return dynamic_cast<FOrbital*>(m_Data[iMo].get());
}


void IFrame::add_bond(int iAt, int jAt, QString const &Flags)
{
   FGeometry
      *pGeometry = pGetGeometry();
   if (pGeometry)
      pGeometry->AddBond(iAt, jAt, Flags);
}

void IFrame::delete_bond(int iAt, int jAt)
{
   FGeometry
      *pGeometry = pGetGeometry();
   if (pGeometry)
      pGeometry->DeleteBond(iAt, jAt);
}

void IFrame::reset_bonds()
{
   FGeometry
      *pGeometry = pGetGeometry();
   if (pGeometry) {
      pGeometry->FindBondLines(m_pDocument);
      pGeometry->UpdateVisualInfo();
   }
}

QString IFrame::get_name()
{
   return m_InputFileName;
}

void IFrame::set_name(QString const &Name)
{
   m_InputFileName = Name;
}




template<class T>
void Scale(T *p, size_t n, T Factor) {
   for (size_t i = 0; i < n; ++ i)
      p[i] *= Factor;
}


void FFrame::MakeOrbitalMatrix(FMatrixView &COrb, FBasisSetPtr &pBasis2, TArray<FOrbital*> *pRefOrbitals, uint32_t Flags, FMemoryStack &Mem)
{
   if (!HaveOrbitals()) {
      COrb = FMatrixView(0,0,0);
      pBasis2 = 0;
      if (pRefOrbitals)
         pRefOrbitals->clear();
      return;
   }

   bool bIaoBasis = (0 != (Flags & ORBMAT_IaoBasis));

   FFrame
      *pLastFrame = this;

   FBasisSetPtr
      pLastBasis2 = pBasis2;

   // count the number of orbitals and find the common basis size.
//    pBasis2 = 0;
//  ^- why removed? in some instances we call MakeOrbitalMatrix multiple times.
//     this allows cross-call asserting that basis sets are consistent (e.g., between alpha and beta orbtials).
//     Note that some of those orbitals can be 0 (e.g., all-alpha wave functions etc.).
   TArray<FOrbital*>
      RefOrbitals1,
      &RefOrbitals = (pRefOrbitals? *pRefOrbitals : RefOrbitals1);
   RefOrbitals.clear();
   RefOrbitals.reserve(pLastFrame->m_Data.size());
   for (uint i = 0; i < pLastFrame->m_Data.size(); ++ i) {
      FOrbital
         *pOrb = dynamic_cast<FOrbital*>(pLastFrame->m_Data[i].get());
      // filter out orbitals we are not supposed to return.
      if (!pOrb)
         continue;
      if (((Flags & ORBMAT_OccupiedOnly) != 0) && pOrb->info.fOcc <= 0.)
         continue;
      if (((Flags & ORBMAT_VirtualOnly) != 0) && pOrb->info.fOcc > 0.)
         continue;
      FOrbitalSpin
         OrbSpin = pOrb->info.Spin;
      if (((Flags & ORBMAT_AlphaAndClosedOnly) != 0) && (OrbSpin != wfi::ORBSPIN_Alpha && OrbSpin != wfi::ORBSPIN_SpinFree))
         continue;
      if (((Flags & ORBMAT_BetaAndClosedOnly) != 0) && (OrbSpin != wfi::ORBSPIN_Beta && OrbSpin != wfi::ORBSPIN_SpinFree))
         continue;
      if (((Flags & ORBMAT_AlphaOnly) != 0) && (OrbSpin != wfi::ORBSPIN_Alpha))
         continue;
      if (((Flags & ORBMAT_ClosedOnly) != 0) && (OrbSpin != wfi::ORBSPIN_SpinFree))
         continue;
      if (((Flags & ORBMAT_BetaOnly) != 0) && (OrbSpin != wfi::ORBSPIN_Beta))
         continue;
      FBasisSetPtr
         pThisBasis = (bIaoBasis ? pOrb->pMinBasis : pOrb->pBasisSet);
      if (pBasis2.get() == 0) {
         pBasis2 = pThisBasis;
      } else {
         if (pBasis2 != pThisBasis.get())
            throw std::runtime_error("FFrame::MakeOrbitalMatrix: Multiple orbitals in *same* frame expressed wrt. different basis sets.");
      }
      RefOrbitals.push_back(pOrb);
   }
   if (RefOrbitals.size() == 0) {
      // no orbitals of the given class exist.
//       pBasis2 = 0;
      COrb = FMatrixView(0,0,0);
      return;
   }
   uint
      nOrb2 = RefOrbitals.size(),
      nBf2 = pBasis2->nFn();

   // make a matrix for the orbitals on Mem and copy the coefficients into it
   COrb = MakeStackMatrix(nBf2, nOrb2, Mem);
   for (uint i = 0; i < RefOrbitals.size(); ++ i) {
      FOrbital
         *pOrb = RefOrbitals[i];
      TArray<double>
         *pCoeffs = (bIaoBasis? &pOrb->pIaoCoeffs : &pOrb->pCoeffs);
      assert_rt(pCoeffs->size() == nBf2);
      for (uint iBf = 0; iBf < nBf2; ++ iBf)
         COrb(iBf, i) = (*pCoeffs)[iBf];
   }

   if (pLastBasis2.get() != 0 && pBasis2.get() != 0 && pBasis2.get() != pLastBasis2.get())
      throw std::runtime_error("FFrame::MakeOrbitalMatrix: Multiple orbitals in *same* frame expressed wrt. different basis sets (across orbital types).");
}

FWfType FFrame::GetWfType()
{
   size_t
      nAlpha = 0,
      nBeta = 0,
      nClosed = 0,
      nOther = 0;
   FDataSetList::const_iterator
      itDataSet;
   _for_each(itDataSet, m_Data) {
      FOrbital const
         *pOrbital = dynamic_cast<FOrbital const*>(&**itDataSet);
      if (pOrbital == 0)
         continue; // dataset is not an orbital.
      switch (pOrbital->info.Spin) {
         case wfi::ORBSPIN_SpinFree:
            nClosed += 1;
            break;
         case wfi::ORBSPIN_Alpha:
            nAlpha += 1;
            break;
         case wfi::ORBSPIN_Beta:
            nBeta += 1;
            break;
         default:
            nOther += 1;
            break;
      }
   }
//    std::cout << fmt::format("GetWfType(): nA = {}  nB = {}  nC = {}  nO = {}\n", nAlpha, nBeta, nClosed, nOther);
   if ((nAlpha != 0 && nBeta == 0) || (nClosed != 0 && nBeta == 0))
      return wfi::WFTYPE_Rhf;
   if (nAlpha != 0 && nBeta != 0 && nClosed == 0)
      return wfi::WFTYPE_Uhf;
   if (nBeta == 0 && nAlpha == 0 && (nClosed != 0 || nOther != 0))
      // todo: should check if total number of occupied orbitals is smaller than
      // the number of valence orbitals... but I don't have a valence orbital counting
      // function coded up atm.
      return wfi::WFTYPE_Mcscf;
   return wfi::WFTYPE_Other;
}

void FFrame::ConvertToUnrestrictedOrbitals()
{
   if (GetWfType() != wfi::WFTYPE_Rhf)
      return;
   FDataSetList
      ExtraOrbs;
   FDataSetList::const_iterator
      itDataSet;
   _for_each(itDataSet, m_Data) {
      FOrbital
         *pOrbital = dynamic_cast<FOrbital*>(&**itDataSet);
      if (pOrbital == 0)
         continue; // dataset is not an orbital.
      // copy the previous orbital (...I wonder if that would copy too much? Like GlMeshes etc?)
      FOrbitalPtr
         pNewOrbital(new FOrbital(*pOrbital));
      if (pOrbital->info.Spin == wfi::ORBSPIN_SpinFree) {
         // should apply for closed-shell (fOcc==2) and external orbitals (fOcc==0) orbitals.
         // These turn into occupied (fOcc==1) and virtual (fOcc==0) alpha- and beta orbitals.
         if (pOrbital->info.fOcc != 0.0 && pOrbital->info.fOcc != 2.0)
            IvNotify(NOTIFY_Warning, IvFmt("ConvertToUnrestrictedOrbitals: Encountered unexpected orbital occupation %1 for spin-free orbital %2 in supposedly ROHF-type input orbital set", pOrbital->info.fOcc, pOrbital->GetDesc()));
         pOrbital->info.fOcc *= 0.5;
         pOrbital->info.Spin = wfi::ORBSPIN_Alpha;
         pNewOrbital->info.fOcc *= 0.5;
         pNewOrbital->info.Spin = wfi::ORBSPIN_Beta;
      } else if (pOrbital->info.Spin == wfi::ORBSPIN_Alpha) {
         // should apply for alpha-spin active orbitals (fOcc==1).
         // These turn into occupied alpha (fOcc==1) and virtual beta (fOcc==0) orbitals.
         if (pOrbital->info.fOcc != 1.0)
            IvNotify(NOTIFY_Warning, IvFmt("ConvertToUnrestrictedOrbitals: Encountered unexpected orbital occupation %1 for alpha-spin orbital %2 in supposedly ROHF-type input orbital set", pOrbital->info.fOcc, pOrbital->GetDesc()));
         pNewOrbital->info.fOcc = 1. - pOrbital->info.fOcc;
         pNewOrbital->info.Spin = wfi::ORBSPIN_Beta;
      } else {
         IvNotify(NOTIFY_Warning, IvFmt("ConvertToUnrestrictedOrbitals: Encountered unexpected orbital type %1 in supposedly ROHF-type input orbital set", pOrbital->GetDesc()));
      }
      ExtraOrbs.push_back(pNewOrbital);
   }
   m_Data.insert(m_Data.end(), ExtraOrbs.begin(), ExtraOrbs.end());
}



void FFrame::MakeOrbitalMoments(ct::FMemoryStack &Mem)
{
   if (!HaveOrbitals())
      return;
   void
      *pFreeMe = Mem.Alloc(0);
   // compute center and quadrupole moment for all orbitals in the current frame.
   // Unfortunatelly, with CtInt1e making those matrices is *very* slow for
   // larger molecules.
   ct::FAtomSet
      *pAtoms = pGetAtoms();
   ct::FMatrixView
      COrb;
   FBasisSetPtr
      pBasisSet;
   TArray<FOrbital*>
      RefOrbitals;
   MakeOrbitalMatrix(COrb, pBasisSet, &RefOrbitals, 0, Mem);

   size_t
      nOrb = COrb.nCols,
      nAo = pBasisSet->nFn();

   TArray<double>
      CartMomMatrices;
   {
      TVector3<double>
         ExpansionPoint(0., 0., 0.);
      CartMomMatrices.resize(nAo*nAo*10);
      FDoublettIntegralFactoryMultipoleMoment::FMoment
         CartMoments[10];
      for (size_t i = 0; i < 10; ++ i)
         for (size_t ixyz = 0; ixyz <= 2; ++ ixyz)
            CartMoments[i][ixyz] = ir::iCartPow[i][ixyz];
      FDoublettIntegralFactoryMultipoleMoment
         CartMomFactory(&CartMoments[0], &CartMoments[10], ExpansionPoint);
      pAtoms->Make1eIntegralMatrixMultiComp(&CartMomMatrices[0], *pBasisSet, *pBasisSet, CartMomFactory, 1.0, Mem);

      FStackMatrix
         Values(10, nOrb, &Mem);
      Values.Clear();

      {  // this has hard O(N^3) scaling but is still significantly faster than the O(n^2) variants...
         FStackMatrix
            T1(nAo, nOrb, &Mem);
         for (size_t iComp = 0; iComp < 10; ++ iComp) {
            FMatrixView
               DipN(&CartMomMatrices[nAo*nAo*(iComp)], nAo, nAo);
            Mxm(T1, DipN, COrb);
            for (size_t iOrb = 0; iOrb < nOrb; ++ iOrb)
               Values(iComp, iOrb) = Dot(&T1(0, iOrb), &COrb(0, iOrb), nAo);
         }
      }

      for (size_t iOrb = 0; iOrb < nOrb; ++ iOrb) {
         for (size_t i = 0; i < 3; ++ i)
            RefOrbitals[iOrb]->vDipMom[i] = Values(1+i,iOrb);
      }
      for (size_t i = 0; i < 3; ++ i) {
         for (size_t j = 0; j <= i; ++ j)
         {
            size_t ij;
            TVector3<unsigned> CartMom(0, 0, 0);
            CartMom[i] += 1;
            CartMom[j] += 1;
            for (ij = 4; ij <= 10; ++ ij)
               if (ir::iCartPow[ij][0] == CartMom[0] && ir::iCartPow[ij][1] == CartMom[1] && ir::iCartPow[ij][2] == CartMom[2])
                  break;
            for (size_t iOrb = 0; iOrb < nOrb; ++ iOrb) {
               TVector3<double>
                  Center(RefOrbitals[iOrb]->vDipMom[0], RefOrbitals[iOrb]->vDipMom[1], RefOrbitals[iOrb]->vDipMom[2]);
               double
//                   rsq = Dot(Center, Center),
                  mass = 1.;

               double v = Values(ij,iOrb) - Center[i] * Center[j] * mass;
//                if (i==j)
//                   v += rsq * mass;
               // ^- very strange... the "parallel axis theorem" according to wikipedia
               // gives a different transformation for diagonal elements. But this version here
               // works and their's doesn't.
               RefOrbitals[iOrb]->mQuadMom[i][j] = v;
               RefOrbitals[iOrb]->mQuadMom[j][i] = v;
            }
         }
      }

      for (size_t iOrb = 0; iOrb < nOrb; ++ iOrb)
         RefOrbitals[iOrb]->HaveMoments = true;
   }

   Mem.Free(pFreeMe);
}

int FFrame::FindRow(FDataSet* pSet)
{
   for (int i = 0; i < (int)m_Data.size(); ++ i) {
      if (m_Data[i] == pSet)
         return i;
   }
   IV_NOTIFY1(NOTIFY_Error, "Failed to find data set '%1' in the current frame.", pSet->GetDesc());
   return -1;
}


// predicate: sorts orbitals into original order, before linking them.
struct FSortOrbitalsToOriginalOrder
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
      return pOrbA->info.iOriginalIndex < pOrbB->info.iOriginalIndex;
   }
};


// this sorts & flips the orbitals into a order maximizing the overlap with the
// last frame's orbitals.
void FFrame::LinkOrbitalsToPreviousFrame(FFrame *pPrevFrame, ct::FMemoryStack &Mem)
{
   if (!HaveOrbitals())
      return;
   if (!pPrevFrame) {
      // no reference frame. Don't link the orbitals---bring them into original order
      // in which we got them (before linking visual configs)
      std::stable_sort(m_Data.begin(), m_Data.end(), FSortOrbitalsToOriginalOrder());
      return;
   }

   void
      *pTopOfStack = Mem.Alloc(0);
   // for now just do it with direct overlap. However, we first need to re-assemble
   // the last frame's orbital matrix.
   ct::FMatrixView
      CurOrb, LastOrb;
   FBasisSetPtr
      pBasis1, pBasis2;
   TArray<FOrbital*>
      RefOrbitals1, RefOrbitals2;
   // FIXME: for UHF wave functions this needs to be done for alpha and beta orbitals separately.
   //        Technically, we should probably always do it within a given subspace, but we might
   //        not be able to distinguish active orbitals with 2.0000 from closed shell-orbitals
   //        in every situation, since this information mightnot be provided in the input with
   //        sufficient accuracy...
   this->MakeOrbitalMatrix(CurOrb, pBasis1, &RefOrbitals1, 0, Mem);
   pPrevFrame->MakeOrbitalMatrix(LastOrb, pBasis2, &RefOrbitals2, 0, Mem);

   uint
      nOrb1 = CurOrb.nCols,
      nOrb2 = LastOrb.nCols, // last.
      nBf2 = pBasis2->nFn(),
      nBf1 = pBasis1->nFn();
   // for now assume that format of orbitals is equal...
   assert(nBf2 == LastOrb.nRows);
   if (nBf1 != nBf2)
      return IV_NOTIFY(NOTIFY_Error, "Current and last frame's orbital formats do not match. Currently we do not support aligning them.");
   if (nOrb1 != nOrb2)
      return IV_NOTIFY(NOTIFY_Error, "Current and last frame's number of orbitals do not match. Currently we do not support aligning them.");

   {
      // form the overlap of the last and the current orbitals with respect to the
      // current frame's overlap matrix (note that we should not directly use the
      // cross overlap: due to the geometry change, the core orbitals have moved:
      // cores of the same atoms thus have near zero overlap in different frames!)
      FStackMatrix
         SAo(nBf1, nBf1, &Mem),
         SMo(nOrb1, nOrb2, &Mem);
      pGetAtoms()->MakeOverlapMatrix(SAo, *pBasis1, *pBasis1, Mem);
      ChainMxm(SMo, Transpose(CurOrb), SAo, LastOrb, Mem);

      // delete overlap between orbitals of different spin. It's a hack to get it
      // workting for UHF and co.
      for (size_t iOrb2 = 0; iOrb2 < nOrb2; ++ iOrb2)
         for (size_t iOrb1 = 0; iOrb1 < nOrb1; ++ iOrb1)
            if (RefOrbitals1[iOrb1]->info.Spin != RefOrbitals2[iOrb2]->info.Spin)
               SMo(iOrb1, iOrb2) = 0.;

      std::vector<bool>
         bNewOrbAssigned(nOrb2, false);
      for (uint iRef = 0; iRef < nOrb2; ++ iRef) {
         // find orbital with maximum absolute overlap with respect to iRef.
         // Update: can't I just SVD the matrix to get corresponding orbitals?
         double fMaxOvl = 0.;
         int iMaxOrb = -1;
         for (uint i = 0; i < nOrb1; ++ i) {
            if (bNewOrbAssigned[i])
               continue;
            double fOvl = std::abs(SMo(i,iRef));
            if (fOvl > fMaxOvl) {
               fMaxOvl = fOvl;
               iMaxOrb = i;
            }
         }

         FOrbital
            *pLastOrb = RefOrbitals2[iRef],
            *pCurOrb = RefOrbitals1[iMaxOrb];

         // flip the orbital if this overlap is negative.
         bNewOrbAssigned[iMaxOrb] = true;
         if (fMaxOvl < 0.8) {
//             IV_NOTIFY(NOTIFY_Warning, IvFmt("poor overlap (%1) between %2/%3 and %4/%5",
//                QString::number(fMaxOvl,'g',3), pCurOrb->GetDesc(), this->GetBaseInputFileName(),
//                pLastOrb->GetDesc(), pPrevFrame->GetBaseInputFileName()));
            m_pFrameLog->Write("WARNING: poor overlap ({:.3f}) between {}/{} and {}/{}",
                      fMaxOvl,
                      q2s(pCurOrb->GetDesc()), q2s(this->GetBaseInputFileName()),
                      q2s(pLastOrb->GetDesc()), q2s(pPrevFrame->GetBaseInputFileName()));
         }
         if (SMo(iMaxOrb,iRef) < 0)
            pCurOrb->FlipPhase();
         // link up this orbital's visual config to the last one's
         pCurOrb->LinkVisualConfig(pLastOrb->pVisConfig);


         // put current orb into iRef's row index in the current frame.
         int
            iRefRow = pPrevFrame->FindRow(pLastOrb),
            iCurRow = this->FindRow(pCurOrb);
         m_Data[iCurRow].swap(m_Data[iRefRow]);
      }
   }

   Mem.Free(pTopOfStack);
}


void FDocument::Clear()
{
   BeginTotalReset();
   SetActiveCol(-1);
   SetActiveRow(-1);

   //emit layoutAboutToBeChanged();
//    m_ActiveRow = -1;
//    m_ActiveCol = -1;
   m_Frames.clear();
   m_InputFileName = "";
   m_LoadedDataFiles.clear();
   m_AtomOptions.clear();
   m_iNextOrbitalColor = 0;
   m_AtomGroups.clear();
   EndTotalReset();
//    m_iNextOrbitalColorScheme = 0;
   // leave the current setting. may be useful if someone wants to apply them
   // to multiple files.
//    delete m_pWfOptions;
//    m_pWfOptions = new FWfOptions(this);
}

void FDocument::BeginTotalReset()
{
   beginResetModel();
}
void FDocument::EndTotalReset()
{
   endResetModel();
   emit layoutAboutToBeChanged();
   emit layoutChanged();
   // ^- FIXME: fix up connections in main and then remove these two.
   //emit layoutChanged();
   emit ActiveColChanged(m_ActiveCol);
   emit ActiveRowChanged(m_ActiveRow);
}

void FDocument::RebuildWf(FLogQt &Log)
{
   if (!HaveConsistentFrames()) {
      IvNotify(NOTIFY_Error, "The loaded frames do not have compatible geometries (same number and types of atoms, in same order). RebuildWf not supported.");
      return;
   }

   ct::FHfOptions
      ScfOpt;
   m_pWfOptions->AssignScfOptions(ScfOpt);

   //emit layoutAboutToBeChanged();
   BeginTotalReset();

   if (m_pWfOptions->GetRunScf()) {
      ClearOrbitals(false); // false: no signals.
   } else {
      // delete all meshes, colors, etc.
      // UPDATE: hmpf... requires the "changing everything" signals, and we can't actually send them from here.
//       ClearOrbitalRepresentations();
   }

   // need to re-set the maximum number of threads. OMP seems to think
   // that it shouldn't spawn more threads, otherwise, if not invoked from main thread.
//    omp_set_num_threads(g_nMaxOmpThreads);
   omp_set_num_threads(m_pWfOptions->GetNumThreads());

   bool ErrorOccured = false;
   std::string ErrorDesc;

   FFrame
      *pLastFrame = 0;
   FHfResultPtr
      pLastWf; // for starting guess purposes.
   FMemoryStack2
      Mem(size_t(m_pWfOptions->GetNumThreads() * m_pWfOptions->GetWorkSpaceMb())<<20);
   for (size_t iFrame = 0; iFrame < m_Frames.size(); ++ iFrame)
   {
      FFrame
         *pFrame = GetFrame(iFrame);
      Log.Write("*** CURRENT FRAME: {} ({} of {})", q2s(RemovePath(pFrame->GetFullInputFileName())), iFrame+1, m_Frames.size());

      // forward main log to frame's local computation log.
      connect(&Log, SIGNAL(textEmitted(QString)), &pFrame->Log(), SLOT(appendText(QString)));

      try {
         ct::FAtomSet
            *pAtoms = pFrame->pGetAtoms();
         if (!pAtoms)
            continue;
         // reset the default basis sets.
         m_pWfOptions->AssignBasisSets(pAtoms);

         ct::TMemoryLock<char>
            pFreeMe(0, &Mem);
         ct::FWfDecl
            WfDecl;
         m_pWfOptions->AssignWfDecl(WfDecl, pAtoms);

         if (m_pWfOptions->GetRunScf()) {
            FHfResultPtr
               pWf = new FHfResult;

            ct::FHfMethod
               Hf(Log, WfDecl, pAtoms, ScfOpt, pLastWf.get(), 0, Mem);
            pWf = Hf.pResult();

            pLastWf = pWf;
            pAtoms->SetLastEnergy(pWf->Energy); // probably should not be here... (and neither the gradients)

            // insert the orbitals as data sets.
            // WARNING: for correct energies in anti-bonding orbtials, I
            // need to store ALL orbitals here! Including virtuals (otherwise
            // there is not enough stuff to assemble the valence-virtual basis Fock matrix!).
            size_t
               nOrbsToLoad = pWf->Orb.nCols;
            if (m_SkipVirtualOrbitals)
               nOrbsToLoad = WfDecl.nElecA();
            for ( uint i = 0; i < nOrbsToLoad; ++ i ) {
//             for ( uint i = 0; i < WfDecl.nElecA(); ++ i ) {
   //          for ( uint i = 0; i < Hf.pMinBasis->nFn(); ++ i ) {
               double
                  fOcc = 0.;
               if (i < WfDecl.nElecA())
                  fOcc = 1.;
               if (i < WfDecl.nElecB())
                  fOcc = 2.;
               double
                  fEnergy = pWf->Ew[i];
               FOrbitalInfo
                  info = FOrbitalInfo(fEnergy, fOcc, 0);
               info.iOriginalIndex = int(i);
               pFrame->m_Data.push_back( new FOrbital(info.MakeDesc(i), pAtoms, this, pWf->pOrbBasis,
                     &pWf->Orb(0,i), info, 0, 0, 0) );
            }
         }

//          pFrame->RunIaoAnalysis(pFrame->Log(), m_pWfOptions, Mem);
         pFrame->RunIaoAnalysis(Log, m_pWfOptions, m_pWfOptions->GetRunIbba(), Mem); // last: decides on whether or not to localzie orbs.
         pFrame->MakeOrbitalMoments(Mem); // <- this is not required (will be made on the fly)

         pFrame->LinkOrbitalsToPreviousFrame(pLastFrame, Mem);
         pLastFrame = pFrame;
      } catch (std::exception &e) {
         std::cout << "** ENTERED RebuildWf//ErrorHandler." << std::endl;
         ErrorDesc = e.what();
//          IvNotify(NOTIFY_Error, IvFmt("RebuildWf failed: %1", s2q(e.what())));
         // ^- can't call this outside the gui thread...
         ErrorOccured = true;
      }
      disconnect(&Log, SIGNAL(textEmitted(QString)), &pFrame->Log(), SLOT(appendText(QString)));

      if (ErrorOccured)
         break;

      Log.Write("\n");
      if (iFrame != m_Frames.size()-1)
         Log.endSection();
   }

   if (ErrorOccured) {
      ClearOrbitals(false);
      Log.EmitError("RebuildWf failed: " + ErrorDesc);
   }
//    Log.Write("\nAll Frames Done.");
   //emit layoutChanged();
   //endResetModel();
   EndTotalReset();

   // this is to update what we can select in the UI.
   if (m_ActiveCol == -1)
      SetActiveCol(0);
   SetActiveRow(0);
}



void FDocument::RunRedoxChargeAnalysis(FLogQt &Log)
{
   bool ErrorOccured = false;
   std::string ErrorDesc;

   FMemoryStack2
      Mem(size_t(m_pWfOptions->GetNumThreads() * m_pWfOptions->GetWorkSpaceMb())<<20);
   for (size_t iFrame = 0; iFrame < m_Frames.size(); ++ iFrame)
   {
      FFrame
         *pFrame = GetFrame(iFrame);
      Log.Write("*** CURRENT FRAME: {} ({} of {})", q2s(RemovePath(pFrame->GetFullInputFileName())), iFrame+1, m_Frames.size());
      Log.WriteLine();

      // forward main log to frame's local computation log.
      connect(&Log, SIGNAL(textEmitted(QString)), &pFrame->Log(), SLOT(appendText(QString)));

      try {
         ct::FAtomSet
            *pAtoms = pFrame->pGetAtoms();
         if (!pAtoms)
            continue;
         FAtomIdList
            iSelectedAtoms = GetSelectedAtoms(false);

         pFrame->RunRedoxChargeAnalysis(Log, &iSelectedAtoms, m_pFragmentAnalysisOptions, &Mem);
      } catch (std::exception &e) {
         std::cout << "** ENTERED AnalyzeWf//ErrorHandler." << std::endl;
         ErrorDesc = e.what();
//          IvNotify(NOTIFY_Error, IvFmt("RebuildWf failed: %1", s2q(e.what())));
         // ^- can't call this outside the gui thread...
         ErrorOccured = true;
      }
      disconnect(&Log, SIGNAL(textEmitted(QString)), &pFrame->Log(), SLOT(appendText(QString)));

      if (ErrorOccured)
         break;

      Log.Write("\n");
      if (iFrame != m_Frames.size()-1)
         Log.endSection();
   }

}



void FDocument::AcquireFrameObjectLock(FDataSetList &ObjectLock)
{
   ObjectLock.clear();
   FFrameList::iterator
      itFrame;
   size_t
      nObjects = 0;
   for (itFrame = m_Frames.begin(); itFrame != m_Frames.end(); ++ itFrame)
      nObjects += (*itFrame)->m_Data.size();
   ObjectLock.reserve(nObjects);
   for (itFrame = m_Frames.begin(); itFrame != m_Frames.end(); ++ itFrame) {
      FDataSetList::iterator
         itData;
      for (itData = (*itFrame)->m_Data.begin(); itData != (*itFrame)->m_Data.end(); ++ itData)
         ObjectLock.push_back(*itData); // this copies the smart pointer in order to add a reference to it.
   }
}


void FDocument::MakeOrbitalMoments()
{
   assert_rt(!"not implemented.");
}

void FDocument::MakeOrbitalCharges()
{
   assert_rt(!"not implemented.");
}

void FDocument::MakeOrbitalCachedData()
{
   // should go through all frames and fix them up.
   MakeOrbitalMoments();
   MakeOrbitalCharges();
}

void FDocument::BeginInsertFrames()
{
   // takes care of QAbstractTableModel interface requirements.
   if (m_CallEndInsertRows != -1)
      return IV_NOTIFY(NOTIFY_Error, "BeginInsertFrames already called");
   m_CallEndInsertRows = 0;

}

void FDocument::InsertFrame(FFramePtr pFrame)
{
   // check if there was a last frame; either something already loaded
   // or something inserted before during the current insert operation block.
   FFramePtr
      pLastFrame(0);
   if (!m_FramesToInsert.empty())
      pLastFrame = m_FramesToInsert.back();
   else {
      if (!m_Frames.empty())
         pLastFrame = m_Frames.back();
   }

   // do some post-processing with actual computations...
   // Temporary data for visualizations, intermediate data
   // and characterizations of the orbitals, etc.
   {
      FMemoryStack2
//          Mem(200000000); // ~200 mb
         Mem(size_t(m_pWfOptions->GetWorkSpaceMb())<<20ul);
      pFrame->RunIaoAnalysis(pFrame->Log(), m_pWfOptions, false, Mem); // false: do not make IBOs by default. Just read what's in the input.
      // (^- note: this is called for *all* loaded files.)
      pFrame->MakeOrbitalMoments(Mem);
      pFrame->LinkOrbitalsToPreviousFrame(pLastFrame.get(), Mem);
   }

   // remember the new frame, but don't insert it into the document just yet.
   // (doing that later simplifies dealing with the model interface)
   m_FramesToInsert.push_back(pFrame);
}


void FDocument::EndInsertFrames()
{
   // tell our connected clients that the data model's structure is about to change.
   bool emit1 = false;
   if (emit1) emit layoutAboutToBeChanged();

   if (m_CallEndInsertRows == -1)
      return IV_NOTIFY(NOTIFY_Error, "BeginInsertFrames not called");
   if (!m_FramesToInsert.empty()) {
      int
         nRowsOld = rowCount(),
         nRowsNew = nRowsOld;

      int
         iNewCol = (int)m_Frames.size();
      if (emit1) {
         beginInsertColumns(QModelIndex(), m_Frames.size(), m_Frames.size() + m_FramesToInsert.size() - 1); // last is exclusive.
         // count old number of rows and new number of rows.
         for (size_t iNewFrame = 0; iNewFrame < m_FramesToInsert.size(); ++iNewFrame)
            nRowsNew = std::max(nRowsOld, m_FramesToInsert[iNewFrame]->rowCount());
         if (nRowsNew != nRowsOld) {
            if (emit1) beginInsertRows(QModelIndex(), nRowsOld, nRowsNew - 1); // end exclusive.
            m_CallEndInsertRows = 1;
         }
      } else
         beginResetModel();

      // now actual insert of the new frames into the model
      m_Frames.insert(m_Frames.end(), m_FramesToInsert.begin(), m_FramesToInsert.end());
      m_FramesToInsert.clear();
      UpdateArcLengths();

      //emit dataChanged(createIndex(0, m_ActiveCol), createIndex(pFrameData->size()-1, m_ActiveCol));
      // ^- that cannot be emitted before the layout change is complete, can it?
      //    And anyway, is that required? Also, if I do the layoutChanged thing, do I
      //    still need to deal with rowsAboutToBeInserted etc?
      if (emit1) {
         if (m_CallEndInsertRows == 1)
            endInsertRows();
         endInsertColumns();
         emit layoutChanged();
      } else {
         endResetModel();
         emit layoutAboutToBeChanged();
         emit layoutChanged();
         // ^- FIXME: fix up connections in main and then remove these two.
      }
      // set the last row of the first frame as current.
//       FDataSetList
//          *pFrameData = &m_Frames.back()->m_Data;

      // switch to the last row of the first new frame as active.
//       IvEmit("!!set active col/row: %1 %2", iNewCol, GetFrame(iNewCol)->rowCount() - 1);
      SetActiveCol(iNewCol);
//       SetActiveRow(GetFrame(iNewCol)->rowCount() - 1);

//1   SetActiveCol(0);
//1   SetActiveRow(pFrameData->size() - 1);
   //    if (size_t(m_ActiveCol) >= m_Frames.size())
   //       m_ActiveCol = m_Frames.size() - 1;
   //    if (size_t(m_ActiveRow) >= GetCurrentFrame()->m_Data.size())
   //       m_ActiveRow = m_Frames.size() - 1;
   //    emit ActiveRowChanged(m_ActiveRow);
   //    emit ActiveColChanged(m_ActiveCol);
   //    int nRow1 = this->rowCount(), nCol1 = this->columnCount();
//1      emit dataChanged(createIndex(0, m_ActiveCol), createIndex(pFrameData->size()-1, m_ActiveCol));
//1      emit VisualRepresentationChanged();
   }
   m_CallEndInsertRows = -1;
}


FOrbitalSpin ConvertXmlSpinDecl(orbital_file::FOrbitalSpin s)
{
   switch (s) {
      case orbital_file::ORBSPIN_SpinFree: return wfi::ORBSPIN_SpinFree;
      case orbital_file::ORBSPIN_Alpha: return wfi::ORBSPIN_Alpha;
      case orbital_file::ORBSPIN_Beta: return wfi::ORBSPIN_Beta;
      case orbital_file::ORBSPIN_Unknown: return wfi::ORBSPIN_Unknown;
      default:
         assert(0);
         return wfi::ORBSPIN_Unknown;
   }
}

bool FDocument::LoadOrbitalFile(FFrameList &LoadedFrames, QString FileName, QString OrigFileName)
{
   // at least the basic molpro XMLs made via "{put,xml,...}" currently only
   // have one frame per file. But they clearly could support more. So we
   // should keep this reasonably flexible.
   if (OrigFileName == "")
      // might differ if this comes from a converted file which is not read directly (e.g., Orca .gbw or Turbomole 'control' files)
      OrigFileName = FileName;

   unsigned
      FileLoadFlags = 0;
   if (m_SkipVirtualOrbitals)
      FileLoadFlags |= orbital_file::LOADFILE_SkipVirtualOrbs;
   orbital_file::FMolproXmlDataPtrList
      XmlFrames;
   try {
      XmlFrames = orbital_file::LoadOrbitalFile(q2s(FileName), orbital_file::FLoadOptions(FileLoadFlags));
   } catch (orbital_file::FFileTypeUnrecognizedError &e) {
      // return that this was not an orbital file. It is called from LoadFile,
      // which may now try some other formats.
      return false;
   } catch (std::runtime_error const &e) {
      IV_NOTIFY(NOTIFY_Error, "LoadOrbitalFile failed: " + s2q(e.what()));
      return true;
   }
   if (XmlFrames.empty()) {
      IV_NOTIFY(NOTIFY_Warning, "LoadOrbitalFile confused: " + QString("File load returned no data: '%1'").arg(FileName));
      return true;
   }

   // add make IboView Frame objects from data sets in file
   for (size_t iXmlFrame = 0; iXmlFrame < XmlFrames.size(); ++ iXmlFrame) {
      orbital_file::FMolproXmlDataPtr
         pXmlData = XmlFrames[iXmlFrame];
      if (pXmlData.get() == 0) {
         IV_NOTIFY(NOTIFY_Error, "LoadOrbitalFile failed: " + QString("Unknown error during loading of '%1'").arg(FileName));
         return true;
      }

      QString
         FrameName;
      if (XmlFrames.size() > 1)
         FrameName = QString("%1_%2").arg(OrigFileName).arg((int)iXmlFrame, 4, 10, QChar('0'));
      else
         FrameName = OrigFileName;
      FFramePtr
         pFrame(new IFrame(FrameName, this));

      FDataSetList
         *pFrameData = &pFrame->m_Data;
      ct::TIntrusivePtr<FGeometry>
         pGeometry = new FGeometry(RemovePath(OrigFileName), pXmlData->pAtoms, this);
      pFrameData->push_back(pGeometry);
      pGeometry->Active = true;

      // insert the orbitals as data sets (if provided).
      if (pXmlData->pOrbSet.get() != 0) {
         for ( uint i = 0; i < pXmlData->pOrbSet->OrbInfos.size(); ++ i ) {
            orbital_file::FOrbitalInfo
               *pOrbInfo = &*pXmlData->pOrbSet->OrbInfos[i];
            FOrbitalInfo
               info = FOrbitalInfo(pOrbInfo->fEnergy, pOrbInfo->fOcc, pOrbInfo->iSym, ConvertXmlSpinDecl(pOrbInfo->Spin));
            info.iOriginalIndex = int(i);
      //       QString
      //          Desc = s2q(pOrbInfo->Desc());
            QString
               Desc = info.MakeDesc(s2q(pOrbInfo->Desc()));
            pFrameData->push_back( new FOrbital(Desc, pXmlData->pAtoms, this, pOrbInfo->pBasisSet,
                  &pOrbInfo->Orb[0], info, 0, 0, 0) );
         }
      }

      // insert the frame into the "new stuff" list.
      LoadedFrames.push_back(pFrame);
   }

   return true;
}

void FDocument::LoadXyzFile(FFrameList &LoadedFrames, QString FileName)
{
   FBasisDescs
      DefaultBases;
   DefaultBases[BASIS_Orbital] = "def2-TZVP"; // FIXME: should be taken as property from somewhere...
   DefaultBases[BASIS_Guess] = "MINAO";

   ct::FAtomSetList
      AtomSets;
   try {
      ct::LoadMultiXyz(AtomSets, q2s(FileName), DefaultBases);
   } catch (xyz_io::FXyzLoadException &e) {
      IvNotify(NOTIFY_Error, s2q(e.what()));
      return;
   }

   ct::FAtomSetList::iterator
      itAtomSet;
   size_t
      nAtomSets = AtomSets.size(),
      iAtomSet = 0;
   _for_each(itAtomSet, AtomSets) {
      // use either .xyz file's name itself as frame name, or if the .xyz
      // contains multiple frames, then its file name with the frame id appended.
      QString
         FrameName;
      if (nAtomSets > 1)
         FrameName = QString("%1_%2").arg(FileName).arg((int)iAtomSet, 4, 10, QChar('0'));
      else
         FrameName = FileName;
      // make a new frame object, consisting of only a geometry as data set.
      FFramePtr
         pFrame(new IFrame(FrameName, this));
      FDataSetList
         *pFrameData = &pFrame->m_Data;
      ct::TIntrusivePtr<FGeometry>
         pGeometry = new FGeometry(RemovePath(FrameName), *itAtomSet, this);
      pFrameData->push_back(pGeometry);
      pGeometry->Active = true;
      // insert the frame into the "new stuff" list.
      LoadedFrames.push_back(pFrame);
      ++ iAtomSet;
   }
}

char const
   // QRegExp-style regex accepting a float number. Inserted into other REs.
   // Comments:
   // - Does not contain the capturing parenthesis around the entire expression
   // - This is intentionally non-locale-aware. Otherwise scripts from one
   //   country might not work in another, due to using different decimal symbols
   //   or digit characters, etc.
   //
   // - see notes/re-test.py for tests.
   //
   // - I think this should match any standard floating point number in
   //   regular or scientific notation, including in notations such as '+.1',
   //   '-.123e+456', etc.
   //   If I am not mistaken, this will not accept any strings which would *not*
   //   be parsed as a float literal in C++.
   // Note:
   // - (?: ..) are non-capturing parenthesis
   // - might need a ^...$ to match start/end of string only
   *g_ReDecl_Float = "[-+]?(?:[0-9]+\\.[0-9]*|\\.[0-9]+|[0-9]+)(?:[eE][+-]?[0-9]+)?",
   // see notes/re-test.py for tests.
   *g_ReDecl_Int = "[-+]?[0-9]+|0x(?:[0-9]|[a-f]|[A-F])+",
   // this takes variable/function identifiers, except for allowing for '-' as part of the inside of a variable name
   *g_ReDecl_Identifier = "[a-zA-Z_](?:[\\-][a-zA-Z]|[a-zA-Z0-9_])*";
   // Note: There is also a QRegularExpression... but that's only in QT5+, not QT4.
   // And we mean to maintain QT4 compatibility here.
   //
   // Also, after hacking this up, it became clear, once again, that I should *REALLY*
   // not deal with regular expressions. This would have taken 1/10th of the time with CxParse1,
   // and be more powerful and flexible, and factor 100 faster...


struct FNameArgListRe
{
   enum FFlags {
      SkipWhitespace = 0x0001,
      CaptureArgs = 0x0002,
      CaptureName = 0x0004,
      WholeString = 0x0008,
      DefaultFlags = SkipWhitespace | CaptureName | CaptureArgs | WholeString
//       DefaultFlags = CaptureName | CaptureArgs | WholeString
   };
   QString
      reOpen, reDelim, reClose;
   unsigned
      Flags;
   explicit FNameArgListRe(QString reOpen_, QString reDelim_, QString reClose_, unsigned Flags_ = DefaultFlags)
   : reOpen(reOpen_), reDelim(reDelim_), reClose(reClose_), Flags(Flags_)
   {
      // WARNING: might need escape for literal symbols in regex!
      // e.g, use "\\(", ",",  "\\)" for (,) separators.
   }

   template<class... Args>
   QRegExp re(QString reName, Args... args)
   {
      QStringList sArgs{args...};
      QString sFn;
      if (bool(this->Flags & CaptureName))
         sFn = QString("(%1)").arg(reName);
      else
         sFn = reName;
      QString
         sBegin, sEnd,
         sArgLeft, sArgJoin, sArgRight;
      if (bool(this->Flags & WholeString)) {
         sBegin = "^";
         sEnd = "$";
      }
      if (bool(this->Flags & CaptureArgs)) {
         sArgLeft = "(";
         sArgRight = ")";
      }
      if (bool(this->Flags & SkipWhitespace)) {
         sArgLeft = "\\s*" + sArgLeft;
         sArgRight = sArgRight + "\\s*";
         sBegin = sBegin + "\\s*";
         sEnd = "\\s*" + sEnd;
      }
      sArgJoin = QString("%1%2%3").arg(sArgRight, this->reDelim, sArgLeft);
      QString sRegExp = QString("%1%2%3%4%5%6%7%8").arg(sBegin, sFn, this->reOpen, sArgLeft, sArgs.join(sArgJoin), sArgRight, this->reClose, sEnd);
      return QRegExp(sRegExp);
   }
};


void FDocument::AddAxes(QString Which, double fAxisLength_, QString Options)
{
   BeginInsertFrames();
   // ^-- hm... maybe QVariantMaps, made directly from object literals (add_axes("xyz", {length: 4., style:gray, origin: [-4.,-4.,0.]}), etc., might be a better option?
   // note: there are no named function parameters in ECMA script. Using object
   // literals for that appears to be the standard method to achieve similar
   // results.
   double
//       fAxisLength = AxisLength_,
      AxisWidth = 0.03, // axis line length & width
      AxisWeight = 1.; // dots, etc.
   double
      fAxisLabelOffs = 0.5, fAxisLabelSize = 0.6;
   FVec3d
      vAxisOrigin(0,0,0),
      // controls the fractions of the axis length attributed to the
      // negative and positive fractions of the coordinate space.
      // of 0, only the positive side will be shown. If 0.5, half the length
      // will be used for positive and half for negative parts of the axis.
      // Note: by default, setting an 'origin' argument sets this to (0,0,0).
      // To get both origin and frac, use frac after origin.
      vAxisLengthFrac(0.5, 0.5, 0.5),
      vAxisLength(fAxisLength_, fAxisLength_, fAxisLength_);
   bool
      bAxisLengthFracSetExplicitly = false;
   uint32_t
      dw_xAxis = 0xffff0000, dw_yAxis = 0xff00ff00, dw_zAxis = 0xff0000ff;
   float
      fBrightnessMod = 0.;
//       dw_xLabel = 0xffff0000, dw_yLabel = 0xff00ff00, dw_zLabel = 0xff0000ff,
//    fBrightnessMod = -.3;
   QStringList
      FlagList = Options.split("|", QString::SkipEmptyParts);

   QRegExp
      // function with one floating point argument
      reFunc1f = FNameArgListRe("\\(", "\\,", "\\)").re(g_ReDecl_Identifier, g_ReDecl_Float),
      // function with three floating point arguments
      reFunc3f = FNameArgListRe("\\(", "\\,", "\\)").re(g_ReDecl_Identifier, g_ReDecl_Float, g_ReDecl_Float, g_ReDecl_Float);

   foreach(QString Flag1, FlagList) {
      QString Flag = Flag1.trimmed();
      if (Flag == "gray" || Flag == "grey") {
         dw_xAxis = dw_yAxis = dw_zAxis = 0xffbfbfbf;
         fBrightnessMod = -.4;
         continue;
      }
      if (reFunc1f.exactMatch(Flag)) {
         QString FnName = reFunc1f.cap(1);
         double fValue = reFunc1f.cap(2).toDouble();
//          IvEmit("  '%2' matches reFunc1f ('%1')\n  (cap(1): '%3', cap(2): '%4')--> w = %5", reFunc1f.pattern(), Flag, reFunc1f.cap(1), reFunc1f.cap(2), fValue);
         if (FnName == "width") {
            AxisWidth = fValue; continue;
         } else if (FnName == "dotted" || FnName == "weight") {
            AxisWeight = fValue; continue;
         } else if (FnName == "length") {
            vAxisLength = FVec3d(fValue, fValue, fValue); continue;
         } else if (FnName == "label-size" || FnName == "label_size") {
            fAxisLabelSize = fValue; continue;
         } else if (FnName == "label-offs") {
            fAxisLabelOffs = fValue; continue;
         }
      } else if (reFunc3f.exactMatch(Flag)) {
         QString
            FnName = reFunc3f.cap(1);
         double
            fValue0 = reFunc3f.cap(2).toDouble(),
            fValue1 = reFunc3f.cap(3).toDouble(),
            fValue2 = reFunc3f.cap(4).toDouble();
         FVec3d
            vValue012 = FVec3d(fValue0, fValue1, fValue2);
//          IvEmit("  '%2' matches reFunc3f ('%1')\n  (caps: '%3')--> w = %4", reFunc3f.pattern(), Flag, QString("['%1', '%2', '%3', '%4']").arg(reFunc3f.cap(1)).arg(reFunc3f.cap(2)).arg(reFunc3f.cap(3)).arg(reFunc3f.cap(4)), QString("(%1,%2,%3,%4)").arg(FnName).arg(vValue012[0]).arg(vValue012[1]).arg(vValue012[2]));
         if (FnName == "origin" || FnName == "zero") {
            vAxisOrigin = vValue012;
            if (!bAxisLengthFracSetExplicitly)
               vAxisLengthFrac = FVec3d(0,0,0);
            continue;
         } else if (FnName == "frac" || FnName == "divide") {
            vAxisLengthFrac = vValue012;
            bAxisLengthFracSetExplicitly = true;
            continue;
         } else if (FnName == "length") {
            vAxisLength = vValue012; continue;
         }
      }
      IvNotify(NOTIFY_Warning, QString("Axis option '%1' not recognized. Ignored!").arg(Flag));
   }

   uint32_t
      dw_xLabel = ModBrightness(FColor(dw_xAxis), fBrightnessMod).uint32(),
      dw_yLabel = ModBrightness(FColor(dw_yAxis), fBrightnessMod).uint32(),
      dw_zLabel = ModBrightness(FColor(dw_zAxis), fBrightnessMod).uint32();
   {
      for (size_t iFrame = 0; iFrame < m_Frames.size(); ++iFrame) {
         FDataSetList *pFrameData = this->GetFrameData(iFrame);
         ct::TIntrusivePtr<FFreeObjectSet>
            pFreeObjects = new FFreeObjectSet("Axes", m_Frames[iFrame]->pGetAtoms(), this);
         // hm... not sure how to handle cloning/copying of data (especially if
         // from scripts), and/or sharing of data between different frames
         FVec3d
            vStart[3], vEnd[3], vLabel[3];
         for (unsigned ixyz = 0; ixyz != 3; ++ ixyz) {
            double fn = vAxisLengthFrac[ixyz], fp = 1. - fn;
            vStart[ixyz] = vAxisOrigin; vStart[ixyz][ixyz] -= fn*vAxisLength[ixyz];
            vEnd[ixyz]   = vAxisOrigin; vEnd[ixyz][ixyz] += fp*vAxisLength[ixyz];
            vLabel[ixyz] = vEnd[ixyz]; vLabel[ixyz][ixyz] += fAxisLabelOffs;
         }
         if (Which.contains("x")) {
            pFreeObjects->Take(FFreeObjectPtr(new FFreeLine(
               vStart[0], vEnd[0], dw_xAxis, AxisWidth, AxisWeight)));
            pFreeObjects->Take(FFreeObjectPtr(new FFreeLabel(
               "x", vLabel[0], fAxisLabelSize, dw_xLabel)));
         }
         if (Which.contains("y")) {
            pFreeObjects->Take(FFreeObjectPtr(new FFreeLine(
               vStart[1], vEnd[1], dw_yAxis, AxisWidth, AxisWeight)));
            pFreeObjects->Take(FFreeObjectPtr(new FFreeLabel(
               "y", vLabel[1], fAxisLabelSize, dw_yLabel)));
         }
         if (Which.contains("z")) {
            pFreeObjects->Take(FFreeObjectPtr(new FFreeLine(
               vStart[2], vEnd[2], dw_zAxis, AxisWidth, AxisWeight)));
            pFreeObjects->Take(FFreeObjectPtr(new FFreeLabel(
               "z", vLabel[2], fAxisLabelSize, dw_zLabel)));
         }
//          pFreeObjects->Take(FFreeObjectPtr(new FFreeLabel(
//             "0", FVec3d(0,0,0), fAxisLabelSize, 0xff7f00ff)));
         pFrameData->push_back(pFreeObjects);
         pFreeObjects->Active = true;
      }
   }
   EndInsertFrames();
}


void FDocument::LoadFile(FFrameList &LoadedFrames, QString FileName)
{
   if (FileName != ":/!clipboard!") {
      // guess the file type based on file name and extension
      // and forward to appropriate load routine.
      QFileInfo
         FileInfo(FileName);
      // first.. is this file actually there and readable? if we handle this
      // here we can't mess it up in the sub-routines.
      if (!FileInfo.isFile())
         return IV_NOTIFY1(NOTIFY_Error, "LoadFile failed: Input '%1' is not a file.", FileName);
      if (!FileInfo.isReadable())
         return IV_NOTIFY1(NOTIFY_Error, "LoadFile failed: Input '%1' is a file, but not readable (permissions okay?).", FileName);
      // now let's see what we got...
      // (todo: it might be helpful to have a way of overwriting
      //  the default choice based on file extensions. But for now there are
      //  more important usability tasks to take care of.)
      IvNotify(NOTIFY_StartWork, IvFmt("Loading %1...", FileName));
      m_LoadedDataFiles.append(FileName);
      QString
         FileExt = FileInfo.suffix();
//       if (FileExt == "xml" || FileExt == "molden")
//          // note: LoadOrbitalFile can also look at the first line of the file to help with identification,
//          //       maybe we should pass more stuff here and catch its exceptions?
//          return LoadOrbitalFile(LoadedFrames, FileName);
      if (FileInfo.completeBaseName() == "control") {
         // Turbomole control file (supposedly?). Convert to molden format.
         FTm2MoldenSubProcess
            SubProcess(FileName, 0);
//          if (!SubProcess.Succeeded())
         if (SubProcess.Run()) {
            LoadOrbitalFile(LoadedFrames, SubProcess.ConvertedFileName(), FileName);
         }
         return;
      } else if (FileInfo.suffix() == "gbw") {
         // Orca .gbw file
         FOrca2MklSubProcess
            SubProcess(FileName, 0);
         if (SubProcess.Run()) {
            LoadOrbitalFile(LoadedFrames, SubProcess.ConvertedFileName(), FileName);
         }
         return;
      } else if (LoadOrbitalFile(LoadedFrames, FileName)) {
         return;
      } else if (FileExt == "xyz")
         return LoadXyzFile(LoadedFrames, FileName);
      return IV_NOTIFY1(NOTIFY_Error, "LoadFile failed: File extension of '%1' not recognized.", FileName);
   } else {
      // we currently only support loading .xyz files from clipboard data.
      return LoadXyzFile(LoadedFrames, FileName);
   }
}


void FDocument::Load(QStringList FileNames)
{
   if (FileNames.empty())
      return;

//    m_InputFileName = FileName;
   m_InputFileName = "";
   BeginInsertFrames();
   {
      FFrameList
         LoadedFrames;
      foreach(QString FileName, FileNames)
         LoadFile(LoadedFrames, FileName);
      if (LoadedFrames.size() > 1)
         IvNotify(NOTIFY_StartWork, IvFmt("Processing %1 frames...", LoadedFrames.size()));
      for (size_t iFrame = 0; iFrame < LoadedFrames.size(); ++iFrame)
         InsertFrame(LoadedFrames[iFrame]);
   }

   EndInsertFrames();
   IvNotify(NOTIFY_FinishWork, IvFmt("Finished loading %1 files.", FileNames.size()));

   emit ActiveColChanged(m_ActiveCol); // that's for updating the file name in the title...
}


void FDocument::LinkOrbitals()
{
   if (!HaveOrbitals())
      return;
   ClearOrbitalRepresentations();

   FFrame
      *pLastFrame = 0;
   FMemoryStack2
      Mem(size_t(m_pWfOptions->GetNumThreads() * m_pWfOptions->GetWorkSpaceMb())<<20);
   for (size_t iFrame = 0; iFrame < m_Frames.size(); ++ iFrame) {
      FFrame
         *pFrame = GetFrame(iFrame);
      FLogQt
         &Log = pFrame->Log();
      Log.Write("*** ORBITAL LINK TO PREV. FRAME", q2s(RemovePath(pFrame->GetFullInputFileName())), iFrame+1, m_Frames.size());
      pFrame->LinkOrbitalsToPreviousFrame(pLastFrame, Mem);
      pLastFrame = pFrame;
   }
}


void FDocument::ClearOrbitalRepresentations()
{
   for (int iFrame = 0; iFrame < GetNumFrames(); ++ iFrame){
      FDataSetList
         *pData = GetFrameData(iFrame);
      if (pData) {
         for (int iRow = 0; size_t(iRow) < pData->size(); ++ iRow) {
            (*pData)[iRow]->InvalidateRenderCache();
         }
      }
   }
}


void FDocument::ReorderOrRestrictFrameSet(FFrameIndexList const &iNewIndices)
{
   // check if the currently active column is still in the new frame set.
   // if yes, find its new frame number.
   int iPrevActiveCol = GetActiveColIndex();
   int iNewActiveCol = -1;
   for (size_t i = 0; i < iNewIndices.size(); ++ i)
      if (iNewIndices[i] == iPrevActiveCol)
         iNewActiveCol = int(i);
   emit layoutAboutToBeChanged();

   FFrameList
      NewFrames;
   NewFrames.reserve(iNewIndices.size());
   for (size_t i = 0; i < iNewIndices.size(); ++ i) {
      size_t
         iFrame = size_t(iNewIndices[i]);
      if (iFrame < m_Frames.size()) {
         NewFrames.push_back(m_Frames[iFrame]);
      } else {
         IvNotify(NOTIFY_Warning, IvFmt("Attempted to include frame %1 in new frame order, but there is no such frame (have: %2). Ignored.", int(iFrame), int(m_Frames.size())));
      }
   }
   m_Frames.swap(NewFrames);
   UpdateArcLengths();

   emit layoutChanged();

   SetActiveCol(iNewActiveCol); // might still be -1. In this case deselect.
}

















// // fix up phase a vector such it has its largest element (in absolute terms) positive.
// void FixVectorPhase(double *p, size_t n)
// {
//    if (n == 0)
//       return;
//    double
//       fMax = 0.;
//    size_t
//       iMax = 0;
//    for (size_t i = 0; i < n; ++ i) {
//       double fAbs = std::abs(p[i]);
//       if (fAbs > fMax) {
//          iMax = i;
//          fMax = fAbs;
//       }
//    }
//    if (p[iMax] < 0)
//       for (size_t i = 0; i < n; ++ i)
//          p[i] *= -1;
// }
//
//
// static void FindAligningTrafo(double pR[9], double pD[3], FVector3 const *pAtPos, double const *pAtMass, FVector3 const *pAtPosLast, uint nAt, FMemoryStack &Mem)
// {
//    double
//       TotalMass = 0,
//       pIm[9] = {0}; // inertial tensor.
//    FVector3
//       vMassCenter(0.,0.,0.);
//    // compute total mass and center of mass
//    for (uint iAt = 0; iAt < nAt; ++ iAt) {
//       TotalMass += pAtMass[iAt];
//       vMassCenter += pAtMass[iAt] * pAtPos[iAt];
//    }
//    vMassCenter /= TotalMass;
//
//    // copy center translation
//    for (uint i = 0; i < 3; ++ i)
//       pD[i] = vMassCenter[i];
//
//    if (pAtPosLast == 0) {
//       // no reference frame given.
//
//       // compute inertial tensor around the center of mass.
//       for (uint iAt = 0; iAt < nAt; ++ iAt) {
//          FVector3 const
//             vAtPos = pAtPos[iAt] - vMassCenter;
//          double
//             Rsq = Dot(vAtPos,vAtPos);
//          for (uint i = 0; i < 3; ++ i)
//             for (uint j = 0; j < 3; ++ j)
//                pIm[i + 3*j] += pAtMass[iAt] * ((i==j? 1.:0.)*Rsq - vAtPos[i]*vAtPos[j]);
//    //             I += ElementMasses[o.iElement]*(np.eye(3)*np.dot(Xyz,Xyz) - np.outer(Xyz,Xyz))
//       }
//
//       // align along axes of inertia--largest moment along z (this way flat
//       // molecules will be aligned in the x/y plane).
//       for (uint i = 0; i < 9; ++ i)
//          pR[i] = +pIm[i];
//       double
//          pEw[3];
//       FMatrixView
//          mR(pR, 3, 3);
//       ct::Diagonalize(mR, pEw, Mem);
//       // fix up phases of the eigenvectors, such that they have their largest element positive.
//       for (uint iEv = 0; iEv < 3; ++ iEv)
//          FixVectorPhase(&mR(0,iEv), 3);
//
//       // check the determinant of the transformation to see if it includes an
//       // inversion. If yes, invert one of the axes. We'd get the wrong stero-isomers
//       // otherwise...
//       if (1) {
//          FVector3
//             v0(mR(0,0), mR(1,0), mR(2,0)),
//             v1(mR(0,1), mR(1,1), mR(2,1)),
//             v2(mR(0,2), mR(1,2), mR(2,2));
//          if (Dot(v0, Cross(v1,v2)) < 0.) {
//             for (size_t i = 0; i < 3; ++ i)
//                mR(i,2) *= -1.;
//          }
//       }
//    } else {
//       // compute overlap matrix between current and last frame (assuming that the last one
//       // has already been moved into its own center of mass)
//       for (uint iAt = 0; iAt < nAt; ++ iAt) {
//          FVector3 const
//             vAtPosI = pAtPos[iAt] - vMassCenter,
//             vAtPosJ = pAtPosLast[iAt];
//          for (uint i = 0; i < 3; ++ i)
//             for (uint j = 0; j < 3; ++ j)
//                pIm[i + 3*j] += pAtMass[iAt] * vAtPosI[i] * vAtPosJ[j];
//       }
//       double
//          pU[9], pVt[9], pSig[3];
//       FMatrixView
//          mR(pR, 3, 3),
//          mU(pU, 3, 3),
//          mI(pIm, 3, 3),
//          mVt(pVt, 3, 3);
//       ComputeSvd(mU, pSig, mVt, mI, Mem);
//       Mxm(mR, mU, mVt);
//       // ^- I think it should be the other way around, but there is a strange Transpose3x3
//       //    in the actual transformation routine...
//    }
//
// }

FArcLengthOptions FDocument::GetArcLengthOptions() const
{
   // assemble ArcLengthOptions object from document properties.
   QString const &WeightMode = m_AtomWeightMode;

   if (WeightMode == "charge" || WeightMode == "by_charge") {
      return FArcLengthOptions(FArcLengthOptions::ATOMWEIGHT_Charge);
   } else if (WeightMode == "mass" || WeightMode == "by_mass") {
      return FArcLengthOptions(FArcLengthOptions::ATOMWEIGHT_IrcAvgMass);
   } else if (WeightMode == "iso_mass") {
      return FArcLengthOptions(FArcLengthOptions::ATOMWEIGHT_IrcIsotopeMass);
   } else if (WeightMode == "none" || WeightMode == "coords") {
      return FArcLengthOptions(FArcLengthOptions::ATOMWEIGHT_Unity);
   } else {
     IV_NOTIFY1(NOTIFY_Warning, "Can only align frames by 'mass', 'iso-mass', 'charge', and 'coords'. Alignment mode '%1' not recognized.", WeightMode);
   }
   return FArcLengthOptions(FArcLengthOptions::ATOMWEIGHT_IrcAvgMass);
}


// void FDocument::SetArcLengthOptions(FArcLengthOptions ArcLengthOptions)
// {
//    m_ArcLengthOptions = ArcLengthOptions;
//    UpdateArcLengths();
// }
//
// FArcLengthOptions FDocument::GetArcLengthOptions() const
// {
//    return m_ArcLengthOptions;
// }


FFrameCoords::FFrameCoords(FGeometry *pGeometry)
   : pAtoms(&*pGeometry->pAtoms),
     pDocument(pGeometry->GetDocument())
{
   // copy assemble positions and weights & copy into continuous vectors.
   pAtMass.reserve(pAtoms->size());
   pAtPos.reserve(pAtoms->size());
   FArcLengthOptions
      AtomWeightOptions = pDocument->GetArcLengthOptions();
   for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt) {
      if ((pDocument->AtomFlags(iAt) & ATOM_NoAlign) != 0)
         // ignore this atom in setting up the transformation.
         continue;
      double
         AtMass = AtomWeightOptions.GetAtomWeight(int(iAt), pDocument, pAtoms);
      pAtMass.push_back(AtMass);
      pAtPos.push_back((*pAtoms)[iAt].vPos);
   }
   assert(pAtMass.size() == pAtPos.size());
}

FFrameCoords::~FFrameCoords()
{}


void FDocument::SetAtomWeightMode(QString o)
{
   if (m_AtomWeightMode != o) {
      m_AtomWeightMode = o;
      emit AtomWeightModeChanged(o);
//       UpdateArcLengths();
   }
}

void FDocument::UpdateArcLengths()
{
   TArray<FArcLength>
      ArcLengths;
   // compute arc lengths according to current prescription
   MakeIrcArcLengths(ArcLengths, this, GetArcLengthOptions());

   if (ArcLengths.empty()) {
      // that's an invalid assignment---happens, for example, if loading multiple frames
      // of different geometries.
      FArcLengthOptions
         ArcLenghtOptions(FArcLengthOptions::ATOMWEIGHT_FrameIndex);
      for (size_t iFrame = 0; iFrame < m_Frames.size(); ++ iFrame)
         m_Frames[iFrame]->SetArcLength(FArcLength(double(iFrame), ArcLenghtOptions));
   } else {
      // copy to frames.
      assert_rt(m_Frames.size() == ArcLengths.size());
      for (size_t iFrame = 0; iFrame < m_Frames.size(); ++ iFrame)
         m_Frames[iFrame]->SetArcLength(ArcLengths[iFrame]);
   }
}


void FDocument::FindAligningTrafo(double pR[9], double pD[3], FFrameCoords *pThis, FFrameCoords const *pLast, FMemoryStack &Mem)
{
   if (pThis->empty()) {
//       std::cerr << "WARNING: Current frame has no atoms. Nothing to align to!" << std::endl;
      IV_NOTIFY(NOTIFY_Warning, "Current frame has no atoms. Nothing to align to!");
      // return an identity transformation.
      memset(pR, 0, sizeof(*pR)*9);
      for (uint i = 0; i < 3; ++ i) {
         pR[4*i] = 1.;
         pD[i] = 0.;
      }
   } else {
      // make the actual transformation.
      if (pLast) {
         if (pThis->IsAlignableWith(*pLast)) {
            if (pLast->pAtMass.size() != pThis->pAtMass.size())
               IvNotify(NOTIFY_Error, "FindAligningTrafo: Size mismatch in frame alignment between current and last frame.");
         } else {
            IvNotify(NOTIFY_Warning, "FindAligningTrafo: Incompatible numbers/types of atoms between frames. Will *NOT* adjust frames to each other, but process them independently.");
            pLast = 0;
         }
      }
//       ::FindAligningTrafo(pR, pD, &pThis->pAtPos[0], &pThis->pAtMass[0], pLast? &(pLast->pAtPos[0]) : 0, uint(pThis->pAtMass.size()), Mem);
      ct::FindAtomSetAligningTrafo(pR, pD, &pThis->pAtPos[0], &pThis->pAtMass[0], pLast? &(pLast->pAtPos[0]) : 0, size_t(pThis->pAtMass.size()), Mem);
   }
}


// find index in p[i] with largest absolute value.
template <class FScalar>
static size_t ArgAbsMax(FScalar const *p, size_t n, size_t nStride = 1) {
   assert(n != 0); // not defined in this case.
   FScalar
      fMax = std::abs(p[0]);
   size_t
      iMax = 0;
   for (size_t i = 1; i != n; ++ i) {
      FScalar
         f = std::abs(p[i*nStride]);
      if (f > fMax) {
         fMax = f;
         iMax = i;
      }
   }
   return iMax;
}

static void TransformOrbBasisAndCoeffs(FBasisSet *&pCommonBasis, FBasisSetPtr &pBasis, TArray<double> &pCoeffs, double R[9], double d[3], FMemoryStack &Mem)
{
   // transform the basis set(s), unless they are shared and already have
   // been transformed.
   if (pCommonBasis && pCommonBasis != pBasis.get())
      throw std::runtime_error(fmt::format("FDocument::AlignFrames: expected basis sets to be shared, but found something strange. (pCommon = {:x}, pBasis = {:x})", intptr_t(pCommonBasis), intptr_t(pBasis.get())));
   else if (pCommonBasis == 0) {
      pCommonBasis = pBasis.get();
      pBasis->Transform_3x4(R, d, Mem);
   }

   // transform the orbital's MO coefficients.
   pBasis->TransformMoCoeffs_3x4(FMatrixView(&pCoeffs[0], pBasis->nFn(), 1), R, Mem);
}

FFrameCoordsPtr FDocument::MakeFrameCoords(int iFrame)
{
   if (iFrame < 0 || size_t(iFrame) >= m_Frames.size())
      return 0;
   FGeometry
      *pGeometry = m_Frames[iFrame]->pGetGeometry();
//    FDataSetList
//       &FrameData = m_Frames[iFrame]->m_Data;
//    if (!FrameData.empty())
//       pGeometry = dynamic_cast<FGeometry*>(FrameData[0].get());
   if (pGeometry == 0) {
      IV_NOTIFY(NOTIFY_Warning, IvFmt("WARNING: Frame %1 has no geometry to align! Transformation skipped.", iFrame));
      return 0;
   }
   assert(pGeometry != 0);

   FFrameCoordsPtr
      pThis = new FFrameCoords(pGeometry);
   return pThis;
}


void FDocument::AlignFrames(QString Mode)
{
//    IvNotify(NOTIFY_Information, IvFmt("FDocument::AlignFrames(): Enter. %1 Frames.", m_Frames.size()));
//    IvEmit("AlignFrames invoked!");
   FMemoryStack2
      Mem(20000000); // ~20 MB.

   if (!Mode.isEmpty())
      SetAtomWeightMode(Mode);

//    double
//       LastEvs[9] = {0}; // FIXME: do I still need this?
   FFrameCoordsPtr
      pLast;
   for (uint iFrame = 0; iFrame < m_Frames.size(); ++ iFrame) {
//       IvNotify(NOTIFY_Information, IvFmt("FDocument::AlignFrames(): Start process frame %1", iFrame));
      FDataSetList
         &FrameData = m_Frames[iFrame]->m_Data;
      FFrameCoordsPtr
         pThis = MakeFrameCoords(iFrame);
      if (pThis == 0) {
         IV_NOTIFY(NOTIFY_Warning, IvFmt("WARNING: Frame %1 has no geometry to align! Transformation skipped.", iFrame));
         continue;
      }

      double
         R[9], d[3];
//       ct::PrintMatrixGen(std::cout, R, 3, 1, 3, 3, "Aligning Trafo (in)");
      FindAligningTrafo(R, d, pThis.get(), pLast.get(), Mem);
//       ct::PrintMatrixGen(std::cout, R, 3, 1, 3, 3, "Aligning Trafo (out)");
#if 0
      if (iFrame != 0) {
         // in semi-degenerate cases some axes might have been swapped. align transformation with
         // previous transformation (technically we should do this in each degenerate subspace
         // individually... but atm I do not care sufficiently).
         FMatrixView
            mLastR(&LastEvs[0],3,3),
            mCurR(&R[0],3,3);
         FStackMatrix
            S(3,3, &Mem);
         Mxm(S, Transpose(mCurR), mLastR);
         bool
            iUsed[3] = {false,false,false};
         uint
            iOrd[3];
         for (uint iEv = 0; iEv < 3; ++ iEv) {
            // find other vector with largest alignment to previous ev.
            uint iLast;
            for (;;) {
               iLast = uint(ArgAbsMax(&S(0,iEv), 3));
               if (iUsed[iLast])
                  // already used for another dimension
                  S(iLast,iEv) = 0.;
               else {
                  // found one.
                  iUsed[iLast] = true;
                  break;
               }
            }
            iOrd[iEv] = iLast;
         }
         double
            R_Orderd[9];
         for (uint iEv = 0; iEv < 3; ++ iEv)
            for (uint iComp = 0; iComp != 3; ++ iComp)
               R_Orderd[iComp + 3*iEv] = R[iComp + 3*iOrd[iEv]];
         // ^- FIXME: this entire alignment stuff was written after 12hr of work.
         //    Chances are that there are lots of errors. Check when given time.
         assert(sizeof(R) == sizeof(R_Orderd));
         memcpy(R, R_Orderd, sizeof(R));
      }
#endif
      // FIXME: DO THIS ALIGNING STUFF BEFORE FINDING MATCHING ORBITALS AND MAKING THE MULTIPOLE MOMENTS!!!
      //
      // Note: We could easily store atomic flags of atoms we do not yet have loaded.
      //       The flags are stored as a map anyway!
      //
      // .. hm, I think it may be too complicated. May be better to postpone the
      //    /matching orbital visual link determination, and to allow calling it
      //    manually. Would anyway be the required way if making orbitals ourselves.

      // transform the original atom set (should be shared in all other data sets, too).
      IvEmit("* Aligned frame %1 (Weight = %3, pAtoms = %2)\n", iFrame, fmtp(pThis->pAtoms), GetAtomWeightMode());
      ct::PrintMatrixGen(std::cout, R, 3, 1, 3, 3, "Aligning trafo");
      pThis->pAtoms->Transform_3x4(R, d, Mem);

      FBasisSet
         *pCommonOrbBasis = 0,
         *pCommonMinBasis = 0;
      // transform the basis sets and the MO coefficients (basis sets should be shared, too).
      for (uint iMo = 0; iMo != FrameData.size(); ++ iMo) {
         FOrbital
            *pOrb = dynamic_cast<FOrbital*>(FrameData[iMo].get());
         if (pOrb == 0)
            continue; // that's not an orbital.
//          IvEmit("  ...processing orbital %1   (pBasis = %2, nCoeffs = %3, pMinBasis = %4, pIaoCoeffs = %5)", pOrb->GetDesc(), intptr_t(pOrb->pBasisSet.get()), pOrb->pCoeffs.size(), intptr_t(pOrb->pMinBasis.get()), pOrb->pIaoCoeffs.size());
         assert(pOrb->pCoeffs.size() != 0);
         TransformOrbBasisAndCoeffs(pCommonOrbBasis, pOrb->pBasisSet, pOrb->pCoeffs, R, d, Mem);
         if (pOrb->pIaoCoeffs.size() != 0) // <- not there for any non-valence virtual orbitals
            TransformOrbBasisAndCoeffs(pCommonMinBasis, pOrb->pMinBasis, pOrb->pIaoCoeffs, R, d, Mem);
         pOrb->HaveMoments = false;
         // ^- precomputed values are broken now.

         // fixme: call pFrame->LinkOrbitalsToPreviousFrame(pLastFrame, Mem);
         // again.
      }
//       if (t3x3)
//          ct::Transpose3x3(R);

      // update directions and positions of potential bonds
      for (size_t iDataSet = 0; iDataSet < FrameData.size(); ++ iDataSet) {
         FGeometry
            *pGeometry = dynamic_cast<FGeometry*>(FrameData[iDataSet].get());
         if (pGeometry)
            pGeometry->UpdateVisualInfo();
      }

      // TODO: add check if there are other data sets which are neither orbitals nor
      // geometries...

//       assert(sizeof(R) == sizeof(LastEvs));
//       memcpy(LastEvs, R, sizeof(R));
//       pLast = pThis;
      pLast = MakeFrameCoords(iFrame); // rebuild with updated coordinates?
      // pThis still has the ones used before the alignment!
   }
   UpdateArcLengths();
}





// Set Out := L^T In R.
void BasisChange2( FMatrixView &Out, FMatrixView const &L,
    FMatrixView const &In, FMatrixView const &R, FMemoryStack &Mem)
{
    assert( L.nRows == In.nRows && R.nRows == In.nCols );
//     Out = FMatrixView( 0, L.nCols, R.nCols );
//     Mem.Alloc(Out.pData, Out.GetStridedSize());
    assert( Out.nRows == L.nCols && Out.nCols == R.nCols );

    if ( L.nRows * L.nCols <=  R.nRows * R.nCols ) {
        // T1 := L^T In
        FStackMatrix
            T1(L.nCols, In.nCols, &Mem);
        Mxm(T1, Transpose(L), In);
        Mxm(Out, T1, R);
    } else {
        // T1 := In * R
        FStackMatrix
            T1(In.nRows, R.nCols, &Mem);
        Mxm(T1, In, R);
        Mxm(Out, Transpose(L), T1);
    }

    if ( L.pData == R.pData && In.IsSymmetric(1e-10) )
       Symmetrize(Out);
}

FDataSetPtr FDocument::GetActiveDataSet()
{
   return GetRow(m_ActiveRow, false);
};

void FDocument::ToggleActiveDataRow()
{
   ToggleDataRow(m_ActiveRow);
}

void FDocument::ToggleDataRow(int iIndex)
{
//    xout << fmt::format("FDocument::ToggleDataRow(): On entry m_ActiveRow = {}", m_ActiveRow) << std::endl;
   // toggle the row in all frames.
   bool
      ExistedSomewhere = false;
   for (unsigned iFrame = 0; iFrame < m_Frames.size(); ++ iFrame) {
      FFrame *pFrame = GetFrame(iFrame);
      if (iIndex >= 0 && size_t(iIndex) < pFrame->m_Data.size() ) {
         pFrame->m_Data[iIndex]->Active = !pFrame->m_Data[iIndex]->Active;
         ExistedSomewhere = true;
      }
   }

   if (!ExistedSomewhere) {
      xout << "document: tried to toggle an invalid data entry " << iIndex << ". Does not exist in any frame." << std::endl;
//       xout << "document: tried to toggle an invalid data entry " << iIndex << ". Have only " << m_Data.size() << "." << std::endl;
   }

   // note: emit dataChanged first to render the orbital; otherwise the visual info would not actually be there.
   // (yes, it's hacky.)
   if (ExistedSomewhere) {
      emit dataChanged(createIndex(m_ActiveRow,0),createIndex(m_ActiveRow, m_Frames.size()-1)); // note: bottomRight is INCLUSIVE, not exclusive.
   }

//    xout << fmt::format("FDocument::ToggleDataRow(): invoking SetActiveRow({}). Before: m_ActiveRow = {}", iIndex, m_ActiveRow) << std::endl;
   SetActiveRow(iIndex, true, false); // FIXME: don't update status because it would overwrite the iso trace result

//    emit ActiveDatasetChanged();
}


void FDocument::SetActiveCol(int iIndex)
{
   if (iIndex == m_ActiveCol)
      return; // do nothing---in particular do not invoke dataChanged();

   if (0 == GetFrame(iIndex,false)) {
      if (iIndex != -1)
         xout << "document: no frame #" << iIndex << " exists." << std::endl;
   } else {
//       emit layoutAboutToBeChanged();
      m_ActiveCol = iIndex;
      FDataSetList *pFrameData = GetCurrentFrameData();
      if (pFrameData) {
         emit dataChanged(createIndex(0,m_ActiveCol),createIndex(pFrameData? (pFrameData->size()-1) : 0, m_ActiveCol));
         // ^- should this be here?
   //       emit layoutChanged();
         // ^- why are these here?
         emit ActiveColChanged(iIndex);
         emit ActiveDatasetChanged();
         if (nSelectedAtoms() != 0)
            UpdateSelectedAtomsStatusText();
         else {
            FAtomSet const
               *pAtoms = m_Frames[iIndex]->pGetAtoms();
            if (pAtoms != 0 && !pAtoms->GetName().empty())
               IvNotify(NOTIFY_Information, IvFmt("Switched to Frame #%1 ['%2']", iIndex, s2q(pAtoms->GetName())));
            else
               IvNotify(NOTIFY_Information, IvFmt("Switched to Frame #%1", iIndex));
         }
      }
   }
}



void FDocument::MoveActiveCol(int iDelta)
{
   if (iDelta == 0 || GetNumFrames() <= 0)
      return;
   int
      iNewActiveCol = m_ActiveCol + iDelta;
   if (iNewActiveCol >= GetNumFrames())
      iNewActiveCol = GetNumFrames() - 1;
   if (iNewActiveCol < 0)
      iNewActiveCol = 0;
   SetActiveCol(iNewActiveCol);
}



void FDocument::SetActiveRow(int iIndex, bool ForceUpdate, bool UpdateStatus)
{
   if (iIndex == m_ActiveRow && !ForceUpdate)
      return; // do nothing---in particular do not emit signals.
   m_ActiveRow = iIndex;
   if (m_ActiveRow == -1)
      IvNotify(NOTIFY_Information, "");
   emit ActiveDatasetChanged();
   emit ActiveRowChanged(iIndex);

   if (GetActiveDataSet().get() && UpdateStatus)
      IvNotify(NOTIFY_Information, IvFmt("Active Set: %1", GetActiveDataSet()->GetDesc(FDataSet::DESC_Full)));
}



void FDocument::SelectAtomGroup(int iGroup)
{
   FAtomGroupMap::const_iterator
      itGroup = m_AtomGroups.find(iGroup);
   if (itGroup == m_AtomGroups.end()) {
      IV_NOTIFY(NOTIFY_Warning, IvFmt("No atom group %1 was defined", (int)iGroup));
      return;
   } else {
      UnselectAll(false); // don't update selection just yet.
      FAtomIdList const
         &AtGroup = itGroup->second;
      for (size_t ii = 0; ii < AtGroup.size(); ++ ii) {
         SelectAtom(AtGroup[ii], SELECT_Add, false); // false: no update just yet.
      }

      UpdateSelectedAtomsStatusText();
      emit SelectionChanged();
   }
}

void FDocument::DefineAtomGroup(int iGroup, FAtomIdList const &iAtoms)
{
   if (iAtoms.empty()) {
      // nothing selected? Delete pervious atom group.
      m_AtomGroups.erase(iGroup);
   } else {
      // store selection under given group id, either creating a new one or
      // overwriting the previous one.
      m_AtomGroups[iGroup] = iAtoms;
   }
}


void FDocument::DefineSelectionAsAtomGroup(int iGroup)
{
   FAtomIdList
      iSelectedAtoms = FDocument::GetSelectedAtoms(true); // get selected atoms---in original selection order.
   DefineAtomGroup(iGroup, iSelectedAtoms);
}

void FDocument::WriteExtendedStateScript(std::ostream &out)
{
   bool First = true;

   // check if we have atom groups defined. If yes, copy their definitions.
   FAtomGroupMap::const_iterator
      itGroup;
   for (itGroup = m_AtomGroups.begin(); itGroup != m_AtomGroups.end(); ++ itGroup) {
      if (First) { out << "\n"; First = false; }

      out << "doc.define_atom_group(" << itGroup->first << ", [";
      FAtomIdList const &AtGroup = itGroup->second;
      for (size_t ii = 0; ii < AtGroup.size(); ++ ii) {
         if (ii != 0)
            out << ",";
         out << (1 + AtGroup[ii]);
      }
      out << "]);\n";
   }
}



void FDocument::SelectAtom(int iAt, FSelectionMode SelectMode, bool EmitUpdate)
{
   if (SelectMode == SELECT_Select) {
      UnselectAll(false); // don't update selection just yet.
      SelectMode = SELECT_Add;
   }

   FAtomOptions
      &AtomOptions = this->AtomOptions(iAt);
   if (SelectMode == SELECT_Toggle)
      AtomOptions.Flags ^= ATOM_Selected;
   else if (SelectMode == SELECT_Add)
      AtomOptions.Flags |= ATOM_Selected;
   else
      IV_NOTIFY(NOTIFY_Warning, IvFmt("SelectMode %1 not recognized in FDocument::SelectAtom", (int)SelectMode));

   // remember atom's position in the sequence of stuff we selected (yes, it is hacky.).
   if (AtomOptions.Flags & ATOM_Selected) {
      AtomOptions.iSelectionSequenceId = m_SelectionSequenceId;
      ++ m_SelectionSequenceId;
   }

//    int nSelectedAtoms_ = nSelectedAtoms();
   if (EmitUpdate) {
      UpdateSelectedAtomsStatusText();
      emit SelectionChanged();
   }
}

void FDocument::UpdateSelectedAtomsStatusText()
{
   int nSelectedAtoms_ = nSelectedAtoms();
   if (nSelectedAtoms_ != 0) {
      // make a list of the currently selected atoms for status purposes.
      IFrame
         *pFrame = GetCurrentFrame();
      if (pFrame) {
         FAtomIdList
            iSelectedAtoms = FDocument::GetSelectedAtoms();
//          FAtomSet const
//             &Atoms = *pFrame->pGetAtoms();
//          if (iSelectedAtoms.size() == 1) {
//             // could get some additional info about the atom, I guess (like charges, valencies, etc?).
// //             FAtom
//          }
         QString s;
         QTextStream str(&s);
         if (nSelectedAtoms_ == 2) {
            str << FMeasureBondLength(iSelectedAtoms[0], iSelectedAtoms[1], this).MeasureFrame(pFrame);
//             if (HaveOrbitals()) {
//                str << " | " << FMeasureBondOrder(iSelectedAtoms[0], iSelectedAtoms[1], this).MeasureFrame(pFrame);
//             }

            // HMPF!! for angles/dihedrals, etc, I need the order in which atoms were selected...
         } else if (nSelectedAtoms_ == 3) {
            str << FMeasureBondAngle(iSelectedAtoms[0], iSelectedAtoms[1], iSelectedAtoms[2], this).MeasureFrame(pFrame);
         } else if (nSelectedAtoms_ == 4) {
            str << FMeasurePlaneAngle(iSelectedAtoms[0], iSelectedAtoms[1], iSelectedAtoms[2], iSelectedAtoms[3], this).MeasureFrame(pFrame);
         } else {
            str << "Selected: ";
            for (size_t ii = 0; ii != iSelectedAtoms.size(); ++ ii) {
               if (ii > 8) {
                  str << "...";
                  break;
               }
               if (ii != 0)
                  str << " | ";
               str << AtomLabel(iSelectedAtoms[ii]);
   //             int
   //                iAt = iSelectedAtoms[ii];
   //             if ((size_t)iAt >= Atoms.size())
   //                str << "?!!ERR";
   //             else
   //                str << s2q(Atoms[iAt].GetElementName());
   //             str << " " << (1+iAt);
            }
         }
         IvNotify(NOTIFY_Information, s);
      }
   } else {
      IvNotify(NOTIFY_Information, "");
   }
}


int FDocument::nSelectedAtoms()
{
   int nSelected = 0;
   FAtomOptionMap::iterator
      itAtomOptions;
   for (itAtomOptions = m_AtomOptions.begin(); itAtomOptions != m_AtomOptions.end(); ++ itAtomOptions)
      if (itAtomOptions->second.Flags & ATOM_Selected)
         nSelected += 1;
   return nSelected;
}

int FDocument::nAtomsWithFlags()
{
   int nAt = -1;
   FAtomOptionMap::iterator
      itAtomOptions;
   for (itAtomOptions = m_AtomOptions.begin(); itAtomOptions != m_AtomOptions.end(); ++ itAtomOptions)
      if (int(itAtomOptions->first) > nAt)
         nAt = int(itAtomOptions->first);
   return nAt + 1;
}


void FDocument::UnselectAll(bool EmitUpdate)
{
   FAtomOptionMap::iterator
      itAtomOptions;
   for (itAtomOptions = m_AtomOptions.begin(); itAtomOptions != m_AtomOptions.end(); ++ itAtomOptions)
      itAtomOptions->second.Flags &= ~ATOM_Selected;
   m_SelectionSequenceId = 0;

   IvNotify(NOTIFY_Information, "");
   SetActiveRow(-1);
   if (EmitUpdate)
      emit SelectionChanged();
}

void FDocument::HideSelectedAtoms()
{
   FAtomOptionMap::iterator
      itAtomOptions;
   for (itAtomOptions = m_AtomOptions.begin(); itAtomOptions != m_AtomOptions.end(); ++ itAtomOptions)
      if (itAtomOptions->second.Flags & ATOM_Selected) {
         itAtomOptions->second.Flags &= ~ATOM_Selected;
         itAtomOptions->second.Flags |= ATOM_Hidden;
      }
   emit SelectionChanged();
}

void FDocument::ChangeSelectedBond()
{
   FBondChangeAction
      *pBondAct = qobject_cast<FBondChangeAction*>(sender());
   if (pBondAct == 0) {
      IV_NOTIFY(NOTIFY_Warning, "Got a change-bond request, but sender is not a FBondChangeAction. Request ignored.");
      return;
   }
   int
      iAt = pBondAct->m_iAt + 1, // public interface has 1-based atom numbers.
      jAt = pBondAct->m_jAt + 1;
   FDocument *document = this;
   switch(pBondAct->m_Type) {
      case FBondChangeAction::ACTION_Hide: {
         for (uint iFrame = 0; iFrame < uint(document->GetNumFrames()); ++ iFrame)
            document->GetFrame(iFrame)->delete_bond(iAt, jAt);
         break;
      }
      case FBondChangeAction::ACTION_SetStyleDotted: {
         for (uint iFrame = 0; iFrame < uint(document->GetNumFrames()); ++ iFrame)
            document->GetFrame(iFrame)->add_bond(iAt, jAt, "gray|dotted");
         break;
      }
      case FBondChangeAction::ACTION_Reset: {
         for (uint iFrame = 0; iFrame < uint(document->GetNumFrames()); ++ iFrame)
            document->GetFrame(iFrame)->add_bond(iAt, jAt, "");
         break;
      }
      case FBondChangeAction::ACTION_SetBondOrder: {
         for (uint iFrame = 0; iFrame < uint(document->GetNumFrames()); ++ iFrame)
            document->GetFrame(iFrame)->add_bond(iAt, jAt, QString("bo(%1)").arg(pBondAct->m_fArg));
         break;
      }
      default: {
         IV_NOTIFY(NOTIFY_Warning, "Type of bond change not recognized. Request ignored.");
      }
   }
   emit VisualRepresentationChanged();
}

void FDocument::MakeBondLinesForSelectedAtoms()
{
   QString
      sNewStyle = "";
   QAction
      *pAct = qobject_cast<QAction*>(sender());
   if (pAct) {
      sNewStyle = pAct->data().toString();
   }

   FAtomIdList
      iSelectedAtoms = FDocument::GetSelectedAtoms();
   if (iSelectedAtoms.size() < 2u)
      return;
   int
      iAt = iSelectedAtoms[0]; // first selected is now always in place 0...
   for (size_t j = 1; j < iSelectedAtoms.size(); ++ j) {
      // make bond line from iAt to jAt.
      int jAt = iSelectedAtoms[j];
      for (uint iFrame = 0; iFrame < uint(this->GetNumFrames()); ++ iFrame)
         this->GetFrame(iFrame)->add_bond(iAt+1, jAt+1, sNewStyle);
   }
   emit VisualRepresentationChanged();
}

void FDocument::ResetBondLines()
{
   for (uint iFrame = 0; iFrame < uint(this->GetNumFrames()); ++ iFrame)
      this->GetFrame(iFrame)->reset_bonds();
   emit VisualRepresentationChanged();
}


void FDocument::CopySelectedAtomNumbers()
{
   FAtomIdList
      Atoms = GetSelectedAtoms(true); // sorted by selection sequence
   QString
      s;
   QTextStream
      str(&s);
   str << "[";
   for (FAtomIdList::const_iterator itAtom = Atoms.begin(); itAtom != Atoms.end(); ++ itAtom) {
      if (itAtom != Atoms.begin())
         str << ",";
      str << (*itAtom + 1);
   }
   str << "]";
   str.flush();
   QApplication::clipboard()->setText(*str.string());
}



void FDocument::SetNextOrbitalColorIndex(int Value)
{
   m_iNextOrbitalColor = Value;
}

void FDocument::SetNextOrbitalColorScheme(int Index)
{
   m_iNextOrbitalColorScheme = Index;
}

int FDocument::GetNextOrbitalColorScheme()
{
   return m_iNextOrbitalColorScheme;
}


void FDocument::GetNextOrbitalColors(uint32_t &cIsoPlus, uint32_t &cIsoMinus)
{
   uint32_t cAlpha = 0x99000000; // 60%
   float Value = 1.;
   switch (m_iNextOrbitalColorScheme % 3) {
      case 0: { // hue spread
         float Hue = -(m_iNextOrbitalColor+1)*120.f;
         float Spread = 50.f;
         float Saturation = .6f;
         cIsoPlus = (uint32_t) ct::Hsv(Hue + Spread/2, Saturation, Value).uint32() | cAlpha;
         cIsoMinus = (uint32_t) ct::Hsv(Hue - Spread/2, Saturation, Value).uint32() | cAlpha;
         break;
      }
      case 1: { // sat spread
         float Hue = (m_iNextOrbitalColor)*60.f;
         cIsoPlus = (uint32_t) ct::Hsv(Hue, 0.6f, Value).uint32() | cAlpha;
         cIsoMinus = (uint32_t) ct::Hsv(Hue, 0.35f, Value).uint32() | cAlpha;
         break;
      }
      case 2: { // hue spread
         uint32_t cAlpha = 0x7f000000; // 50%
         float Hue = -(m_iNextOrbitalColor+1)*60.f;
//          float Spread = 50.f;
//          float Spread = 0.f;
         float Saturation = .5f;
         float Spread = 180.f;
         cIsoPlus = (uint32_t) ct::Hsv(Hue + Spread/2, Saturation, Value).uint32() | cAlpha;
         cIsoMinus = (uint32_t) ct::Hsv(Hue - Spread/2, Saturation, Value).uint32() | cAlpha;
         break;
      }
   }
   m_iNextOrbitalColor += 1;
   emit NextOrbitalColorIndexChanged(m_iNextOrbitalColor);
}


void FDocument::AddVolumePropertyIsoSurface(FVolumePropertyInfo const &PropertyInfo, FIsoValueList const &IsoValues)
{
   BeginTotalReset();

   // FIXME:
   //  - add explicit IsoValue overrides to VisualConfig
   //    (such that surfaces will not be affected by the view's default config, but rather
   //    by what we explicitly configured)
   //  - Make it possible to edit also volume data other than orbitals via the
   //    current color interface
   //    (that is not QUITE the right thing, but might be good enough for a start.
   //    later we'd like to be able to edit color ramps/gradients for more general properties)
   FVolumeVisualConfigPtr
      // object will be shared by density sets across all frames.
      pVisConfig = new FVolumeVisualConfig();
//    pVisConfig->IsoValues = IsoValues;
//    pVisConfig->bColorSet = true;
//    pVisConfig->iIsoType = ISOTYPE_Absolute;
   pVisConfig->pDetails = FVolumeVisualConfigDetailsPtr(new FMultiPhaseIsoSurfaceConfig(ISOTYPE_Absolute, IsoValues));

   for (int iFrame = 0; size_t(iFrame) < m_Frames.size(); ++ iFrame) {
      FDataSetPtr
         pNewSet(new FVolumeProperty(&*m_Frames[iFrame], PropertyInfo, pVisConfig));
      m_Frames[iFrame]->m_Data.push_back(pNewSet);
      pNewSet->Active = true;
   }

   EndTotalReset();
}




FDocument::FDocument(QObject *parent)
   : QAbstractTableModel(parent),
     m_ActiveRow(-1),
     m_ActiveCol(-1),
     m_SelectionSequenceId(0),
     m_SkipVirtualOrbitals(!g_ShowVirtualOrbitals),
     m_CallEndInsertRows(-1)
{
   int nElements = 110;
   m_ElementOptions.reserve(nElements);
   for (int iElement = 0; iElement < nElements; ++ iElement) {
      // note: 0 is a valid index here. used for dummy atoms.
      m_ElementOptions.append(FElementOptionsPtr(new FElementOptions(iElement)));
   }

   m_iNextOrbitalColor = 0;
   m_iNextOrbitalColorScheme = 0;

   m_pWfOptions = new FWfOptions(this);
   m_pMeasures = new FDocumentMeasures(this,this);
   m_AtomWeightMode = "mass";
   m_pFragmentAnalysisOptions = new FFragmentAnalysisOptions();

   connect(this, SIGNAL(AtomWeightModeChanged(QString)), this, SLOT(UpdateArcLengths(void)));
}

IFrame *FDocument::GetCurrentFrame()
{
   return GetFrame(m_ActiveCol, false);
}

IFrame *FDocument::GetCurrentFrame() const
{
   return const_cast<FDocument*>(this)->GetFrame(m_ActiveCol, false);
}


FDataSetList *FDocument::GetFrameData(int iFrame)
{
   FFrame *pFrame = GetFrame(iFrame,false);
   if (pFrame)
      return &pFrame->m_Data;
   else
      return 0;
}

IFrame *FDocument::GetFrame(int Idx, bool AssertExists)
{
   if (Idx >= 0 && size_t(Idx) < m_Frames.size())
      return m_Frames[Idx].get();
   if (AssertExists)
      throw std::runtime_error("requested non-existent data frame.");
   return 0;
}

IFrame const *FDocument::GetFrame(int Idx, bool AssertExists) const
{
   return const_cast<FDocument*>(this)->GetFrame(Idx, AssertExists);
}

FDataSet const *FDocument::GetRow(int Idx, bool AssertExists) const
{
   return const_cast<FDocument*>(this)->GetRow(Idx, AssertExists);
}

FDataSet const *FDocument::GetRowCol(int iRow, int iCol, bool AssertExists) const
{
   return const_cast<FDocument*>(this)->GetRowCol(iRow, iCol, AssertExists);
}



FDataSetList *FDocument::GetCurrentFrameData()
{
   FFrame *pFrame = GetCurrentFrame();
   if (pFrame)
      return &pFrame->m_Data;
   else
      return 0;
}


FDataSetList FFrame::AllData() {
   FDataSetList pl = m_Data;
   if (m_pFreeObjects && !m_pFreeObjects->empty())
      pl.push_back(m_pFreeObjects);
   return pl;
}


FDataSet *FDocument::GetRowCol(int iRow, int iCol, bool AssertExists)
{
   FFrame *pFrame = GetFrame(iCol, AssertExists);
   if (pFrame) {
      if (size_t(iRow) < pFrame->m_Data.size())
         return pFrame->m_Data[iRow].get();
      if (AssertExists)
         throw std::runtime_error("requested non-existent data row.");
   }
   return 0;
}

FDataSet *FDocument::GetRow(int Idx, bool AssertExists)
{
   return GetRowCol(Idx, m_ActiveCol, AssertExists);
}


int FDocument::rowCount(const QModelIndex & /*parent*/) const
{
   int RowCount = 0;
   for (unsigned i = 0; i < m_Frames.size(); ++ i)
      RowCount = std::max(RowCount, int(m_Frames[i]->m_Data.size()));
   return RowCount;
}

int FDocument::columnCount(const QModelIndex & /*parent*/) const
{
   return m_Frames.size();
}

QString FDocument::GetCurrentInputFileName()
{
   FFrame
      *pFrame = GetCurrentFrame();
   if (pFrame)
      return pFrame->GetFullInputFileName();
   return "";
//    return "unknown_file.xml";
}

QString FDocument::GetCurrentInputBaseFileName()
{
//    std::string InputName = GetCurrentInputFileName();
//    std::size_t iext = InputName.rfind(".");
//    return InputName.substr(0, iext);
   return RemoveExt(GetCurrentInputFileName());
}

QString FDocument::GetCommonInputFileName()
{
   // if we have a script name, then return this script name.
   if (!m_InputFileName.isEmpty())
      return m_InputFileName;

   // otherwise go through all frames and check to which degree their file names
   // are the same, from the front.
   std::string Common = q2s(GetCurrentInputFileName());

   for (int iFrame = 0; size_t(iFrame) < m_Frames.size(); ++ iFrame) {
      size_t n = 0;
      std::string s = q2s(m_Frames[iFrame]->GetFullInputFileName());
      while (n < s.size() && n < Common.size() && s[n] == Common[n])
         n += 1;
      Common.erase(n);
   }
   if (RemovePath(s2q(Common)).isEmpty())
      return "unknown.xml";

   return s2q(Common);
}

size_t FDocument::iFrameId(FFrame *pFrame) const
{
   for (size_t i = 0; i < m_Frames.size(); ++ i)
      if (pFrame == &*m_Frames[i])
         return i;
   return m_Frames.size(); // invalid.
}


bool FDocument::HaveEnergies() const
{
   for (size_t i = 0; i < m_Frames.size(); ++ i)
      if (m_Frames[i]->GetEnergy() != 0.)
         return true;
   return false;
}

bool FDocument::HaveGradients() const
{
   for (size_t i = 0; i < m_Frames.size(); ++ i)
      if (m_Frames[i]->GetGradient() != 0.)
         return true;
   return false;
}

bool FDocument::HaveOrbitals() const
{
   for (size_t i = 0; i < m_Frames.size(); ++ i)
      if (m_Frames[i]->HaveOrbitals())
         return true;
   return false;
}

bool FDocument::HaveConsistentFrames() const {
   for (size_t i = 1; i < m_Frames.size(); ++ i)
      if (!m_Frames[i]->IsGeometryCompatibleWith(*m_Frames[0]))
         return false;
   return true;
}

bool FFrame::IsGeometryCompatibleWith(FFrame const &other) const
{
   FAtomSet const
      *pThis = pGetAtoms(),
      *pOther = other.pGetAtoms();
   if (pThis == 0 || pOther == 0)
      return false;
   return pThis->IsSameMolecule(*pOther);
}




QVariant FDocument::data(const QModelIndex &index, int role) const
{
   int iRow = index.row(),
       iCol = index.column();
   FDataSet const
      *pData = GetRowCol(iRow, iCol, false);
   if (!pData)
      return QVariant();

   if (role == Qt::DisplayRole) { return "[" + pData->GetType() + "] " + pData->GetDesc(); }
   if (pData->Active) {
      uint32_t dwBaseColor = pData->GetBaseColor();
      if (dwBaseColor != 0) {
         if (role == Qt::BackgroundRole) return QBrush(QColor(dwBaseColor));
         if (role == Qt::ForegroundRole) return QBrush(Qt::black);
      } else {
         if (role == Qt::BackgroundRole) return QBrush(QColor(64,96,64));
         if (role == Qt::ForegroundRole) return QBrush(Qt::white);
      }
   }

   return QVariant();
}

QVariant FDocument::headerData(int section, Qt::Orientation orientation, int role) const
{
   if (orientation == Qt::Horizontal) {
      if (role == Qt::DisplayRole) { return IvFmt("F#%1", section); }
//       if ( role == Qt::BackgroundRole ) return QBrush((section % 2 == 0)? Qt::red : QColor(128,0,0));
   }

// //    if ( orientation == Qt::Vertical ) {
// //       if ( role == Qt::DisplayRole ) return QString("[%1]").arg(section);
// //    }
//    if ( orientation == Qt::Horizontal ) {
//       if ( role == Qt::DisplayRole && section == 0 ) { return QString(""); }
//       if ( role == Qt::DisplayRole && section == 1 ) { return QString("desc"); }
//       if ( 0 ) {
//          // custom decorations...
//          if ( role == Qt::FontRole ) {
//                QFont CaptionFont;
//                CaptionFont.setBold(true);
//                return CaptionFont;
//          }
//          if ( role == Qt::BackgroundRole ) return QBrush((section % 2 == 0)? QColor(96,96,96) : QColor(128,128,128));
//          if ( role == Qt::ForegroundRole ) return QBrush(Qt::white);
//       }
// //       if ( role == Qt::BackgroundRole ) return QBrush((section % 2 == 0)? Qt::red : QColor(128,0,0));
//    }

   return QVariant();
}



FElementOptions *pElementOptions(FDocument *pDocument, int iAt, ct::FAtomSet *pAtoms)
{
   return pDocument->pElementOptions(iAt, pAtoms);
}

FElementOptions *FDocument::pElementOptions(int iAt, ct::FAtomSet *pAtoms)
{
   FAtomOptions
      &AtomOpt = AtomOptions(iAt);
   // has there been an override? If yes, return it.
   if (AtomOpt.pPropertiesOverride)
      return &*AtomOpt.pPropertiesOverride;
   // if not, return the default for the atom's element.
   // For this we need to find the element first, however.
   if (pAtoms == 0) {
      IFrame
         *pFrame = GetCurrentFrame();
      if (pFrame)
         pAtoms = pFrame->pGetAtoms();
   }
   if (pAtoms != 0) {
      if ((size_t)iAt < pAtoms->size()) {
         int iElement = (*pAtoms)[iAt].iElement;
         if (iElement >= 0 && iElement < m_ElementOptions.size())
            return &*m_ElementOptions[iElement];
      }
   }
   IvNotify(NOTIFY_Error, IvFmt("Failed to find element data for atom %1. Bad karma!", iAt));
   return 0;
}

QString RemovePath(QString const &FileName) {
   return QFileInfo(FileName).fileName();
}

QString RemoveExt(QString const &FileName) {
   return QFileInfo(FileName).baseName();
}

QString ReplaceExt(QString const &FileName, QString const &NewExt) {
   QString
      BaseName = RemoveExt(FileName);
   if (NewExt == "")
      return BaseName;
   if (NewExt.startsWith("."))
      return BaseName + NewExt;
   return QString("%1.%2").arg(BaseName).arg(NewExt);
}

