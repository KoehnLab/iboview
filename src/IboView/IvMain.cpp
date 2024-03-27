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
#include <clocale>
// ^- for setting locale to C: sometimes locale ends up at local locale.. which breaks float/string conversions.
// see: http://stackoverflow.com/questions/25661295/why-does-qcoreapplication-call-setlocalelc-all-by-default-on-unix-linux
#include <iostream>
#include <QApplication>
#include <QDialog>
#include <QTableView>
#include <QTreeView>
#include <QIcon>
#include <QFileDialog>
#include <QTimer>
#include <QProgressDialog>
#include <QSizePolicy>
#include <QStyleFactory>
#include <QClipboard>
#include <QMessageBox>
#include <QFileInfo>
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QUrl>
#include <QList>
#include <QSize>
#include <QMenu>
#include <QAction>
#include <QShortcut>
#include <QSettings>
// #include <QCommandLineParser>
// ^- it's for > QT 5.2 only :(.
#include <QThreadPool>
#include <QRunnable>
#include <QMimeData>
// #include <QMenuItem>
#include "optionparser.h" // lean mean c++ option parser (see: optionparser.sourceforge.net)

#include "QPropertyModel.h"

#include <QMetaProperty>
#include <QMetaObject>
#include <QMetaMethod>

// #include <fstream>

#include <cmath>
#include <sstream>
#include "format.h"

#include "IvMain.h"

#include "IvVolumeDataSet.h"
#include "IvOrbital.h"

#include "IvShowTextForm.h"
#include "IvSettings.h"
#include "IvLog.h"
#include "ui_MainForm2.h"
#include "ui_AboutForm.h"
#include "IvView3D.h"
#include "IvScript.h"
#include "IvIrc.h"
#include "IvComputeWfForm.h"
#include "IvTables.h"
#include "IvEditFramesForm.h"
#include "IvPreferencesForm.h"
#include "IvEditVolumeSurface.h"
#include "IvComputeEosForm.h"

#include "CtBasisLibrary.h"
#include "CxOpenMpProxy.h"
#include "CxColor.h"
// #include "CxOsInt.h"
using ct::FColor;

std::string q2s(QString const &s) {
   return s.toStdString();
}

QString s2q(std::string const &s) {
//    return QString(s.c_str());
   return QString::fromStdString(s);
}

static bool s_UseStyleFiles = false;
bool g_ShowVirtualOrbitals = true;
int g_nMaxOmpThreads = 1;

// bool ends_with(std::string const &s, std::string const &ending)
// {
//    if (s.length() >= ending.length()) {
//       return (0 == s.compare (s.length() - ending.length(), ending.length(), ending));
//    } else {
//       return false;
//    }
// }

IApplication *FMainWindow::ThisAsIApp()
{
   assert(dynamic_cast<IApplication*>(this) != 0);
   if (dynamic_cast<IApplication*>(this) == 0) {
      std::cerr << "Cast to IApplication failed. Invoked via c'tor?" << std::endl;
   }
   return dynamic_cast<IApplication*>(this);
}


int IApplication::get_orbital_color_scheme()
{
   return document->GetNextOrbitalColorScheme();
}

void IApplication::set_orbital_color_scheme(int i)
{
   return document->SetNextOrbitalColorScheme(i);
}

void IApplication::show_text(QString Caption, QString Text)
{
   FShowTextForm
      ShowTextForm(
         IvFmt("<html><pre>%1</pre></html>", Text),
         Caption, ""
      );
   ShowTextForm.exec();
}


void IApplication::notify(QString Text)
{
   IvNotify(NOTIFY_Information, Text);
}



bool isScriptFile(QString const &FileName) {
   QFileInfo
      FileInfo(FileName);
   QString
      FileExt = FileInfo.suffix();
   return FileExt == "js" || FileExt == "chai";
}

void IApplication::load_file(QString const &FileName) {
   QStringList
      L;
   L.append(FileName);
   return FMainWindow::load_files_(L);
}

void IApplication::add_axes(QString Which, double AxisLength, QString Options)
{
   // e.g.,
   //    doc.load_file("/home/cgk/molecules/acrylic_acid.xyz");
   //    doc.add_axes(4., "xy", "gray|width(0.01)|dotted(0.5)");
   //    doc.update_views();
   document->AddAxes(Which, AxisLength, Options);
}

void IApplication::close_files()
{
   return FMainWindow::close_files_();
}

void FMainWindow::load_files_(QStringList const &FileNames_)
{
   if (FileNames_.empty())
      return;

   QStringList
      DataFiles;
   foreach(QString FileName, FileNames_) {
      if (FileName == "")
         continue;
      if (isScriptFile(FileName)) {
         if (!DataFiles.empty()) {
            // first load all the data files we still have in the cache.
            document->Load(DataFiles);
            DataFiles.clear();
         }
         // now execute the script.
         ExecScript(ThisAsIApp(), this->view3d, FileName);
         document->SetInputFileName(FileName);
      } else {
         // it's a data file. We keep on piling them up until we
         // either have all of them or we encounter a script file.
         DataFiles.append(FileName);
      }
   }

   // load the remaining data files (if there should be any).
   document->Load(DataFiles);
}

void FMainWindow::close_files_()
{
   document->Clear();
}



QStringList toStringList(QScriptValue const &ScriptList)
{
   QStringList
      r;
   int Length = ScriptList.property("length").toInteger();
   for (int i = 0; i < Length; ++ i)
      r.append(ScriptList.property(i).toString());
   return r;
}

FAtomIdList toIntList(QScriptValue const &ScriptList, int iOffset=0)
{
   FAtomIdList
      r;
   int Length = ScriptList.property("length").toInteger();
   r.reserve(size_t(Length));
   for (int i = 0; i < Length; ++ i) {
      QScriptValue
         v = ScriptList.property(i);
      if (v.isNumber()) {
         r.push_back(v.toInt32() + iOffset); // <- note: in JS all numbers are FLOATSs...
      } else {
         IvNotify(NOTIFY_Warning, QString("toAtomIdList: value %1 cannot be interpreted as atom id.").arg(v.toString()));
      }
   }
   return r;
}

FAtomIdList toAtomIdList(QScriptValue const &ScriptList) {
   return toIntList(ScriptList, -1);
}


void IApplication::define_atom_group(int iAtomGroup, QScriptValue const &AtomList)
{
   document->DefineAtomGroup(iAtomGroup, toAtomIdList(AtomList));
}


void IApplication::load_files(QScriptValue const &FileList)
{
   FMainWindow::load_files_(toStringList(FileList));
//    document->Load(toStringList(FileList));
//    int Length = FileList.property("length").toInteger();
//    for (int i = 0; i < Length; ++ i)
//       load_file(FileList.property(i).toString());
}


void IApplication::orient_frames(QString const &Mode)
{
   document->AlignFrames(Mode);
   // this invalidates all already rendered data.
   onRebuildIsoSurfacesClicked();
}


QString IApplication::get_frame_name()
{
   return document->GetCurrentInputBaseFileName();
}


QStringList IApplication::get_loaded_files()
{
   return document->GetLoadedDataFileList();
}


void IApplication::set_loaded_files(QStringList const &FileList)
{
   close_files();
   load_files_(FileList);
}



FOrbital *FMainWindow::pGetOrbital(int iMo)
{
   if (document->GetNumFrames() == 0) {
      IvNotify(NOTIFY_Warning, QString("pGetOrbital: no frames loaded. Cannot get MO #%1").arg(iMo));
      return 0;
   };
   FOrbital *pOut = document->GetCurrentFrame()->pGetOrbital(iMo);
   if (pOut == 0)
      IvNotify(NOTIFY_Warning, QString("pGetOrbital: attemped to access non-existent MO #%1.").arg(iMo));
   return pOut;
}


void IApplication::set_frame(int iFrame)
{
   if (iFrame >= document->GetNumFrames()) {
      IvNotify(NOTIFY_Warning, QString("SetFrame: tried to select frame %i, but valid frame ids are only 0...%i.").arg(iFrame).arg(document->GetNumFrames()));
      return;
   }
   // ^- wtf? how could that possibly have worked in the script?!
   document->SetActiveCol(iFrame);
}


int IApplication::get_frame()
{
   return document->GetActiveColIndex();
}


void IApplication::show_mo(int iMo, /*float fColorPhase, */int cIsoPlus_, int cIsoMinus_)
{
   FOrbital
      *pOrbital = pGetOrbital(iMo);
   if (pOrbital == 0)
      return;

   FVolumeVisualConfig
      *pVisConfig = pOrbital->pVisConfig.get();
   // ^- is that really right? it might have been set beforehand if the VisConfig
   //    object is linked (on the other hand, show_mo toggles the row in all frames...).
//    pVisConfig->SetColorFromCentralHue(fColorPhase);
//    pVisConfig->SetColorFromCentralHue(-1);
   double
      fIsoValue = view3d->m_IsoThreshold/100.;
   FIsoType
      IsoType = ISOTYPE_Relative;
   pVisConfig->pDetails = new FTwoPhaseIsoSurfaceConfig(IsoType, fIsoValue, (uint32_t) cIsoPlus_, (uint32_t) cIsoMinus_);

//    pVisConfig->AssignDefaultColor(pOrbital->GetDocument(), fIsoValue);
// //    pVisConfig->AssignDefaultColor(pOrbital->GetDocument());
// // //    if (cIsoPlus_ != -1) // hm... -1 is actually a valid value (full white).
// //    pVisConfig->cIsoPlus = (uint32_t) cIsoPlus_;
// // //    if (cIsoMinus_ != -1) // hm... -1 is actually a valid value (full white).
// //    pVisConfig->cIsoMinus = (uint32_t) cIsoMinus_;
//    pVisConfig->AssignPlusMinusColors((uint32_t) cIsoPlus_, (uint32_t) cIsoMinus_, fIsoValue);

//    pVisConfig->fIsoValue = fDefaultIsoValue;
//    pVisConfig->iIsoType = DefaultIsoType;

   if (!pOrbital->Active)
      document->ToggleDataRow(iMo); // first (#0) is geometry. But we index MOs 1-based.

   // ^- hm... that's not quite right. will not do anything if orbital had been there before...
   QCoreApplication::processEvents();
}

namespace ct {
   void Rot2x2(double *RESTRICT pA, double *RESTRICT pB, size_t nSize, double phi); // IvIao.cpp
}

void IApplication::rotate_mos_2x2(int iMo, int jMo, double Angle)
{
   // hmpf... there is another one of these: void IFrame::rot_mos_2x2(int iMo, int jMo, double fAngle)

   // note: this function is intentionally hacky. No loop over frames,
   // no recalculation of orbital charges, etc. just for quickly testing some
   // things.
   FOrbital
      *pOrbitalI = pGetOrbital(iMo),
      *pOrbitalJ = pGetOrbital(jMo);
   if (pOrbitalI == 0)
      IvNotify(NOTIFY_Warning, IvFmt("iMo index %1 is not valid. Encountered in rotate_mo_2x2.", iMo));
   if (pOrbitalJ == 0)
      IvNotify(NOTIFY_Warning, IvFmt("jMo index %1 is not valid. Encountered in rotate_mo_2x2.", iMo));
   if (pOrbitalI == 0 || pOrbitalI == 0)
      return;

   double AngleInRadians = M_PI/180. * Angle;

   if (pOrbitalI->pCoeffs.size() == pOrbitalJ->pCoeffs.size()) {
      ct::Rot2x2(&pOrbitalI->pCoeffs[0], &pOrbitalJ->pCoeffs[0], pOrbitalI->pCoeffs.size(), AngleInRadians);
   } else {
      IvNotify(NOTIFY_Warning, IvFmt("Orbital coefficient arrays of iMo=%1 and jMo=%1 are inconsistent. Rot2x2 skipped. Encountered in rotate_mo_2x2.", iMo));
   }

   if (pOrbitalI->pIaoCoeffs.size() == pOrbitalJ->pIaoCoeffs.size()) {
      ct::Rot2x2(&pOrbitalI->pIaoCoeffs[0], &pOrbitalJ->pIaoCoeffs[0], pOrbitalI->pIaoCoeffs.size(), AngleInRadians);
   } else {
      IvNotify(NOTIFY_Warning, IvFmt("Orbital IAO coefficient arrays of iMo=%1 and jMo=%1 are inconsistent. Rot2x2 skipped. Encountered in rotate_mo_2x2.", iMo));
   }

   pOrbitalI->HaveMoments = false;
   pOrbitalJ->HaveMoments = false;
   pOrbitalI->InvalidateRenderCache();
   pOrbitalJ->InvalidateRenderCache();
//    update_views();

   QCoreApplication::processEvents();
}



void IApplication::hide_mo(int iMo)
{
   FOrbital
      *pOrbital = pGetOrbital(iMo);
   if (pOrbital == 0)
      return;
   if (pOrbital->Active)
      document->ToggleDataRow(iMo);
   QCoreApplication::processEvents();
}

static void Scale(double *p, std::size_t N, double f)
{
   for (std::size_t i = 0; i < N; ++ i)
      p[i] *= f;
}

void IFrame::scale_mo(int iMo, double fFactor)
{
   FOrbital
      *pOrbital = pGetOrbital(iMo);
   if (pOrbital == 0)
      return;
   Scale(&pOrbital->pCoeffs[0], pOrbital->pCoeffs.size(), fFactor);
}

template<class FScalar>
inline void rot2x2(FScalar &A, FScalar &B, FScalar cs, FScalar ss) {
   FScalar
      tA =  A * cs + B * ss,
      tB = -A * ss + B * cs;
   A = tA;
   B = tB;
}

template<class FScalar>
inline void rot2x2(FScalar *pA, FScalar *pB, std::size_t N, FScalar cs, FScalar ss) {
   for (std::size_t i = 0; i < N; ++ i)
      rot2x2(pA[i], pB[i], cs, ss);
}

void IFrame::rot_mos_2x2(int iMo, int jMo, double fAngle)
{
   FOrbital
      *pOrbitalI = pGetOrbital(iMo),
      *pOrbitalJ = pGetOrbital(jMo);

   if (pOrbitalI == 0 || pOrbitalJ == 0)
      return;
   double
      cs = std::cos(M_PI/180. * fAngle),
      ss = std::sin(M_PI/180. * fAngle);
   assert(pOrbitalI->pCoeffs.size() == pOrbitalJ->pCoeffs.size());
   assert(pOrbitalI->pIaoCoeffs.size() == pOrbitalJ->pIaoCoeffs.size());
   rot2x2(&pOrbitalI->pCoeffs[0], &pOrbitalJ->pCoeffs[0], pOrbitalI->pCoeffs.size(), cs, ss);
   rot2x2(&pOrbitalI->pIaoCoeffs[0], &pOrbitalJ->pIaoCoeffs[0], pOrbitalI->pIaoCoeffs.size(), cs, ss);
   pOrbitalI->HaveMoments = false;
   pOrbitalJ->HaveMoments = false;
   pOrbitalI->InvalidateRenderCache();
   pOrbitalJ->InvalidateRenderCache();
}


void IApplication::hide_mos()
{
   for (uint iMo = 1; iMo < document->GetCurrentFrame()->m_Data.size(); ++ iMo) {
      FOrbital
         *pOrbital = pGetOrbital(iMo);
      if (pOrbital == 0)
         return;
      if (pOrbital->Active)
         document->ToggleDataRow(iMo);
   }
   QCoreApplication::processEvents();
}


void IApplication::quit()
{
   QCoreApplication::processEvents();
//    QCoreApplication::quit();
   close(); // should work as long as this is the only window and the event loop was already entered...
   std::cout << "!Invoked quit() from script." << std::endl;
}

// // FElementOptionsList &IApplication::element_options()
// QObjectList IApplication::element_options()
// {
// //    return document->GetElementOptions();
//    QObjectList
//       q;
//    FElementOptionsList
//       ol = document->GetElementOptions();
//    for (int i = 0; i < ol.size(); ++ i)
//       q.append(&*ol[i]);
//    return q;
// }

QObject *IApplication::wf_options()
{
   return document->GetWfOptions();
}


QObject *IApplication::element_options(int iElem)
{
//    return document->GetElementOptions();
   FElementOptionsList
      &ol = document->GetElementOptions();
//    return document->pElementOptions(iAt);
   if (iElem >= 0 && iElem < ol.size())
      return ol[iElem].data();
   return 0;
}

void IApplication::reset_element_options(int iElem)
{
   FElementOptionsList
      &ol = document->GetElementOptions();
//    return document->pElementOptions(iAt);
   if (iElem >= 0 && iElem < ol.size()) {
      FElementOptions Dummy(iElem, 0);
      ol[iElem]->CopyPropertiesFrom(Dummy);
   }
}

QObject* IApplication::atom_options(int iAtom_)
{
   int
      iAtom = iAtom_ - 1; // external atom ids expected starting with 1.
   FAtomOptions
      &ao = document->AtomOptions(iAtom);
   if (ao.pPropertiesOverride)
      return &*ao.pPropertiesOverride;
   // detach atom options from element defaults.
   FElementOptions
      *pOrig = document->pElementOptions(iAtom);
   if (!pOrig)
      return 0;
   ao.pPropertiesOverride = new FElementOptions(pOrig, 0);
   return &*ao.pPropertiesOverride;
}

void IApplication::reset_atom_options(int iAtom_)
{
   int
      iAtom = iAtom_ - 1; // external atom ids expected starting with 1.
   FAtomOptions
      &ao = document->AtomOptions(iAtom);
   ao.pPropertiesOverride = 0;
}

void IApplication::update_views()
{
   view3d->update();
}



FGeometry *FMainWindow::pGetActiveGeometry()
{
   if (document->GetNumFrames() == 0)
      return 0;
   return document->GetCurrentFrame()->pGetGeometry();
}

uint IApplication::num_frames()
{
   return document->GetNumFrames();
}

QObject* IApplication::frame(int iFrame)
{
   return document->GetFrame(iFrame);
}

QObject* IApplication::frame()
{
   return document->GetCurrentFrame();
}

void IApplication::add_bond(int iAt, int jAt, QString const &Flags)
{
   for (uint iFrame = 0; iFrame < uint(document->GetNumFrames()); ++ iFrame)
      document->GetFrame(iFrame)->add_bond(iAt, jAt, Flags);
}

void IApplication::delete_bond(int iAt, int jAt)
{
   for (uint iFrame = 0; iFrame < uint(document->GetNumFrames()); ++ iFrame)
      document->GetFrame(iFrame)->delete_bond(iAt, jAt);
}

void IApplication::reset_bonds()
{
   document->ResetBondLines();
}

void IApplication::set_atom_mode(int iAt, QString const &Mode)
{
//       if (iAt == 0 || std::size_t(iAt) > pGeometry->pAtoms->size()) {
//          std::cout << fmt::format("!WARNING: Attempted to change mode of non-exitent atom #{} to '{}'.", iAt, Mode) << std::endl;
//          return;
//       }
   // note:
   //    (1) atom flags are now shared across all frames.
   //    (2) it is allowed to set them *before* the geometry is loaded.
   //        (e.g., in order to change the alignment flags, which need to be present on load).
   if (iAt == 0) {
      IvNotify(NOTIFY_Warning, "Atomic index 0 is not valid. Encountered in set_atom_mode.");
      return;
   }
   uint
      &AtomFlags = document->AtomFlags(iAt-1),
      AtomFlagsOrig = AtomFlags;
   QStringList
      FlagList = Mode.split("|", QString::SkipEmptyParts);
   foreach(QString Flag, FlagList) {
      if (Flag == "hidden")
         AtomFlags = AtomFlags | ATOM_Hidden;
      else if (Flag == "visible")
         AtomFlags = AtomFlags & (~ATOM_Hidden);
      else if (Flag == "noalign" || Flag == "no-align")
         AtomFlags = AtomFlags | ATOM_NoAlign;
      else if (Flag == "align")
         AtomFlags = AtomFlags & (~ATOM_NoAlign);
      else {
         IvNotify(NOTIFY_Warning, QString("Atom mode '%1' not recognized (attemped to set for atom %2). Ignored!").arg(Flag).arg(iAt));
      }
   }
   if (AtomFlags != AtomFlagsOrig) {
      FGeometry
         *pGeometry = pGetActiveGeometry();
      if (pGeometry && size_t(iAt-1) <= pGeometry->pAtoms->size()) {
         ct::FAtom
            &Atom = (*pGeometry->pAtoms)[iAt-1];
         IvEmit(" Changed atom mode %1 %2 to '%3'", iAt, Atom.ElementName(), Mode);
      } else {
         IvEmit(" Stored mode '%1' for use on atom #%2 once loaded.", Mode, iAt);
      }
   }
}

void IApplication::reset_atom_modes()
{
   document->ClearAtomFlags();
}



QString FMainWindow::MakeStateScript()
{
   if (!document->GetCurrentFrameData())
      return "";
   std::stringstream
      out;
//    for (int iFrame = 0; iFrame < document->GetNumFrames(); ++ iFrame)
   if (document->GetNumFrames() == 1)
      out << fmt::format("app.load_file(\"{}\");\n", q2s(document->GetCurrentInputFileName()));
   else {
      out << "app.load_files([";
      for (int iFrame = 0; iFrame < document->GetNumFrames(); ++ iFrame) {
         if (iFrame != 0)
            out << ", ";
         out << fmt::format("\"{}\"", q2s(document->GetFrame(iFrame)->GetFullInputFileName()));
      }
      out << "]);\n";
   }
   out << fmt::format("\n{}\n", q2s(view3d->GetViewDesc()));
   FDataSetList
      &DataList = *document->GetCurrentFrameData();
   bool
      First = true;
   for (uint iDataSet = 0; iDataSet < DataList.size(); ++ iDataSet)
   {
      FDataSetPtr
         pActiveData = DataList[iDataSet];
      if (!pActiveData->Active)
         continue;
      FOrbital
         *pOrbital = dynamic_cast<FOrbital*>(pActiveData.get());
      if ( pOrbital != 0 ) {
         FTwoPhaseIsoSurfaceConfig const
            *pVisDet = dynamic_cast<FTwoPhaseIsoSurfaceConfig const*>(pOrbital->pVisConfig->pDetails.get());
         if (pVisDet) {
            if (First)
               out << "\n";
            out << fmt::format("app.show_mo({},0x{:08x},0x{:08x});\n", iDataSet, pVisDet->cIsoPlus, pVisDet->cIsoMinus);
            First = false;
         }
      }
   }

   document->WriteExtendedStateScript(out);

//    char const *pCropEnabledStr = (view3d->HasFlag(VIEWFLAG_CropImages)? "true" : "false");
   if (document->GetNumFrames() == 1) {
      out << fmt::format("\n// view.save_png(\"{}.png\");\n", q2s(ReplaceExt(document->GetCurrentInputFileName(),".png"))); // % pCropEnabledStr;
   } else {
      out << "\n"
             "for (var iFrame = 0; iFrame < doc.num_frames(); ++iFrame) {\n"
             "   app.set_frame(iFrame);\n"
             "   var frame = app.frame();\n"
             "   // view.save_png(replace_ext(frame.name, \".png\"));\n"
             "}\n";
   }

   out << "\n// kate: syntax javascript;\n";
   return s2q(out.str());
}

void FMainWindow::WriteStateScript(QString FileName)
{
   if (FileName != ":/!clipboard!") {
      if (FileName.isEmpty())
         FileName = ReplaceExt(document->GetCurrentInputFileName(), ".js");
      QFile
         Out(FileName);
      Out.open(QIODevice::WriteOnly | QIODevice::Text);
      Out.write(MakeStateScript().toUtf8());
      IvEmit("* wrote state script to '%1'", FileName);
   } else {
      QApplication::clipboard()->setText(MakeStateScript());
      std::cout << "* copied state script to clipboard." << std::endl;
   }
}

void FMainWindow::dummySignalTest(double o)
{
   IvEmit("!!Received signal: New value is %1", o);
}


// QMetaMethod GetSlotMetaMethod(QObject *pReceiver, char const *pReceiverSlot)
int GetSlotIndex(QObject *pReceiver, char const *pReceiverSlot)
{
   int iMethodIndex = pReceiver->metaObject()->indexOfSlot(QMetaObject::normalizedSignature(pReceiverSlot));
//    int iMethodIndex = pReceiver->metaObject()->indexOfMethod(QMetaObject::normalizedSignature(pReceiverSlot));
//    if (iMethodIndex != -1) {
//       return pReceiver->metaObject()->method(iMethodIndex);
//    } else {
//       IvEmit("!Warning: Failed to find slot '%1' of '%2'", QString(pReceiverSlot), pReceiver->objectName());
//       return QMetaMethod();
//    }
   if (iMethodIndex == -1) {
      IvEmit("!Warning: Failed to find slot '%1' of '%2'", QString(pReceiverSlot), pReceiver->objectName());
   }
   return iMethodIndex;
}

// attempt to connect the notifty signal of a generic sender object's property change for "propertyName"
// to a given receiver object's SLOT(...).
void ConnectPropertyNotify(QObject *pSender, QString propertyName, QObject *pReceiver, char const *pReceiverSlot)
// void ConnectPropertyNotify(QObject *pSender, QString propertyName, QObject *pReceiver, QMetaMethod ReceiverSlot)
{
   int iSlotIndex = GetSlotIndex(pReceiver, pReceiverSlot);
   if (iSlotIndex == -1)
      return;
   QMetaMethod
      ReceiverSlot = pReceiver->metaObject()->method(iSlotIndex);

   const QMetaObject *pMetaObject = pSender->metaObject();
   int count = pMetaObject->propertyCount();
   for (int i = 0; i < count; ++ i) {
      QMetaProperty MetaProperty = pMetaObject->property(i);
      if (propertyName == QString(MetaProperty.name())) {
         if (MetaProperty.hasNotifySignal()) {
            QMetaMethod NotifySignal = MetaProperty.notifySignal();
//             IvEmit("!linked up: %1.notify(%2) to %3::%4", pSender->objectName(), propertyName, pReceiver->objectName(), pReceiverSlot);
//             QObject::connect(pSender, NotifySignal.signature(), pReceiver, pReceiverSlot);
            QObject::connect(pSender, NotifySignal, pReceiver, ReceiverSlot);
         }
      }
   }
}



void LinkPropertyWidgets(QObject *pTarget, QWidget *pWidgetContainer, char const *pPropertyKeyName)
{
   QPropertyDataWidgetMapper
      *viewMapper = QPropertyModel::newMapper(pTarget, pWidgetContainer);
//    IvEmit("viewMapper addr: %1",(size_t)viewMapper);
   // following are some hacks around the data model: QDataWidgetMapper has two "submit changes"
   // policies: (1) submit when the widget loses focus (they call this "auto"),
   // (2) submit manually.
   // The idea that some people might want to submit changes as soon as they happen
   // (i.e., on "valueChanges") has apparently not crossed the developer's minds.
   // So what we do is the following:
   // - change the policy to ManualSubmit
   // - connect all kinds of onChange events we can think of to the mapper's "manual"
   //   submit mechanism.
   // - This has of course the unfortunate side effect of updating *ALL* entries on
   //   *every* change of *any* widget...
   viewMapper->setSubmitPolicy(QPropertyDataWidgetMapper::ManualSubmit);
   QList<QWidget*> uiWidgets = pWidgetContainer->findChildren<QWidget*>();
   foreach(QWidget *w, uiWidgets) {
      QVariant
         vaViewOptionName = w->property(pPropertyKeyName);
      if (vaViewOptionName.isValid()) {
         viewMapper->addMapping(w, vaViewOptionName.toString());
         QAbstractButton
            *pButton = qobject_cast<QAbstractButton*>(w);
         if (pButton) {
            QObject::connect(pButton, SIGNAL(toggled(bool)), viewMapper, SLOT(submit()));
         }
         QDoubleSpinBox
            *pDoubleSpinBox = qobject_cast<QDoubleSpinBox*>(w);
         if (pDoubleSpinBox) {
//             IvEmit("...linking %1 'valueChanged' to 'submit'", pDoubleSpinBox->objectName());
            QObject::connect(pDoubleSpinBox, SIGNAL(valueChanged(double)), viewMapper, SLOT(submit()));
         }
         QSpinBox
            *pSpinBox = qobject_cast<QSpinBox*>(w);
         if (pSpinBox)
            QObject::connect(pSpinBox, SIGNAL(valueChanged(int)), viewMapper, SLOT(submit()));

         QDial
            *pDial = qobject_cast<QDial*>(w);
         if (pDial)
            QObject::connect(pDial, SIGNAL(valueChanged(int)), viewMapper, SLOT(submit()));

         QComboBox
            *pComboBox = qobject_cast<QComboBox*>(w);
         if (pComboBox)
            QObject::connect(pComboBox, SIGNAL(currentIndexChanged(int)), viewMapper, SLOT(submit()));
         // ^- hm.. the QString variant doesn't work. Not that it matters.
         QLineEdit
            *pLineEdit = qobject_cast<QLineEdit*>(w);
         if (pLineEdit)
            QObject::connect(pLineEdit, SIGNAL(editingFinished(void)), viewMapper, SLOT(submit()));
//             QObject::connect(pLineEdit, SIGNAL(textChanged(QString)), viewMapper, SLOT(submit()));
         QPlainTextEdit
            *pPlainTextEdit = qobject_cast<QPlainTextEdit*>(w);
         if (pPlainTextEdit)
            QObject::connect(pPlainTextEdit, SIGNAL(textChanged(QString)), viewMapper, SLOT(submit()));
      }
   }

//    QList<QAction*> uiActions = pWidgetContainer->findChildren<QAction*>();
//    foreach(QAction *pAction, uiActions) {
//       QVariant
//          vaViewOptionName = pAction->property(pPropertyKeyName);
//       if (vaViewOptionName.isValid()) {
//          //       QMenuItem
//          //          *pMenuItem = qobject_cast<QMenuItem*>(pAction);
//          // ^- hm... that doesn't work.
// //          viewMapper->addMapping(pAction, vaViewOptionName.toString());
//          if (pAction) {
//             IvEmit("linked up pAction: %1", pAction->objectName());
//             QObject::connect(pAction, SIGNAL(triggered(bool)), viewMapper, SLOT(submit()));
//          }
//       }
//    }
   // ^- this won't work...
   viewMapper->toFirst();
}






QMainWindow *g_pMainWindow = 0; // for usage as dialog parent.

FMainWindow::FMainWindow(QWidget *parent, Qt::WindowFlags flags)
   : FBase(parent, flags),
     ui(new Ui::MainWindow),
     document(new FDocument(this))
{
   g_pMainWindow = this;
   iUpdateLocked = 0;

   ui->setupUi(this);
   ui->tab_Frames->setEnabled(false); // shouldn't have any frames here. will be updated after OnDataChanged().
   ui->groupBox_ShadingControls->setVisible(false); // that's the secret stuff 8).
   ui->toolButton_ComputeEffectiveOxidationStates->setVisible(false);

   // we have lots of layouts inside other layouts. And, unfortunately,
   // *each layer* of this adds another set of internal margins. Depending
   // on the UI style the result looks very stupid. So we hack this up here.
   // I am not sure how bad this works with high resolution displays (i.e.,
   // if QT resizes the pixels automagically or not for high DPIs)...must
   // be tried... unfortunatelly, I don't have one available atm.
   // Note: These things cannot be controlled via style sheet.
   // docs say that the settings are inherited from the parent layout or
   // parent control... but this doesn't really seem to work for me.
   if (1) {
      QList<QWidget*> uiWidgets = this->findChildren<QWidget*>();
      foreach(QWidget *w, uiWidgets) {
         if (qobject_cast<FStatusBar*>(w) != 0)
            continue; // handles it's own layout.
         //QLayout
         //   *pLayout = qobject_cast<QLayout*>(w);
         QLayout
            *pLayout = w->layout();
         if (pLayout) {
            QGroupBox *pGroupBox = qobject_cast<QGroupBox*>(w);
            if (pGroupBox) {
               QMargins margins = pLayout->contentsMargins();
               margins.setLeft(3);
               margins.setRight(3);
               //pLayout->setContentsMargins(3, 8, 3, 6);
               pLayout->setContentsMargins(margins);
            } else
               pLayout->setContentsMargins(3, 3, 3, 3);
            //qDebug("set spacing of ", pLayout->objectName());
         }
      }
   }
   ui->tab_Orbitals->layout()->setContentsMargins(0, 0, 0, 0); // left, top, right, bottom
//    ui->frame_Toolbar->layout()->setContentsMargins(3, 0, 0, 0);
   // ^- this one leads to a strange OpenGL crash on Johannes' computer.
   // 3d view keeps on resizing itself and goes into recursive repaints, which blow up the GL stack eventually

//    ui->statusBar->setMaximumHeight(ui->toolButton_Dummy->height());

//    ui->centralwidget->layout()->setContentsMargins(3, 3, 3, 0);
   ui->centralwidget->layout()->setContentsMargins(3, 3, 3, 3);

   LoadApplicationSettings();

   dataViewFilter = new FShowOnlyCurrentColumnFilter(0, this);
   dataViewFilter->setSourceModel(document);
//    dataViewFilter->setDynamicSortFilter(true);
//    ui->datasetView->setModel(document);
   ui->datasetView->setModel(dataViewFilter); // show only data rows of currently selected frame.

   //connect(document, SIGNAL(dataChanged(const QModelIndex &,const QModelIndex &)),this, SLOT(onDataChanged()));
   //connect(document, SIGNAL(modelReset()), this, SLOT(onDataChanged()));
   connect(document, SIGNAL(ActiveColChanged(int)), this, SLOT(onDataChanged()));
   // ^- to change UI ranges (e.g., number of frames)

   show();
   QApplication::processEvents();
   // ^- try to force building and rendering the widgets... idea is to
   //    make sure that ui->frame3d is correctly set up before view3d
   //    is constructed. Otherwise there were cases where QGLWidget produced
   //    an invalid gl context.
//    std::cout << boost::format("!frame3d rect: %i %i visible? %s") % ui->frame3d->width() % ui->frame3d->height() % (ui->frame3d->isVisible()? "yes" : "no") << std::endl;
//    view3d = new FView3d(ui->frame3d, document);
   view3d = new IView3d(ui->frame3d, document);
   view3d->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
   ui->frame3d->layout()->addWidget(view3d);
//    view3d->adjustSize();
//    ui->frame3d->addWidget(view3d,0,0);
//    ui->gridLayout_3dFrame->addWidget(view3d,0,0);
//    view3d->updateGeometry(); // maybe required on MacOS X?


   QIcon
      // :/ is the resource name prefix. Note that I could also use it for shader
      // code if I adjust stuff to use QT's file loading functions.
      AppIcon(":/resources/d-lobe2.png");
   this->setWindowIcon(AppIcon);


   // link up events of data set table
   if (1) {
      ui->datasetView->verticalHeader()->hide();
      ui->datasetView->horizontalHeader()->hide();
//       ui->datasetView->resizeColumnsToContents();

      connect(ui->datasetView, SIGNAL(doubleClicked(const QModelIndex &)), this, SLOT(onDataRowDoubleClicked(const QModelIndex &)));
      connect(ui->datasetView, SIGNAL(clicked(const QModelIndex &)), this, SLOT(onDataRowClicked(const QModelIndex &)));
//       connect(document, SIGNAL(dataChanged(const QModelIndex &,const QModelIndex &)), view3d, SLOT(updateData(const QModelIndex &,const QModelIndex &)));
//       connect(document, SIGNAL(layoutChanged()), view3d, SLOT(update()));
//       connect(document, SIGNAL(SelectionChanged()), view3d, SLOT(update()));
      // ^- done in FView3d now (where it probably should have been put right in the beginning).
      connect(document, SIGNAL(ActiveColChanged(int)), this, SLOT(onTitleChanged()));
      connect(dataViewFilter, SIGNAL(layoutChanged()), ui->datasetView, SLOT(resizeColumnsToContents()));
//       connect(document, SIGNAL(layoutChanged()),
//          ui->datasetView, SLOT(resizeColumnsToContents()));
      // ^- doesn't work anymore?
      connect(document, SIGNAL(layoutChanged()), this, SLOT(onActiveIboCurveChanged()));

      ui->datasetView->setShowGrid(false);
      QHeaderView *verticalHeader = ui->datasetView->verticalHeader();
//    verticalHeader->setDefaultSectionSize(verticalHeader->fontMetrics().height()+4);
      verticalHeader->setDefaultSectionSize(int(verticalHeader->fontMetrics().lineSpacing()*1.5));
   }

   connect(ui->pushButton_ShowCurveEnergy, SIGNAL(clicked()), this, SLOT(onActiveIboCurveChanged()));
   connect(ui->pushButton_ShowCurveGradient, SIGNAL(clicked()), this, SLOT(onActiveIboCurveChanged()));

   // connect toolbar/menu actions defined in QtDesigner to local functions.
   connect(ui->actionOpen, SIGNAL(triggered()), this, SLOT(onOpenFile()));
   connect(ui->actionCloseFiles, SIGNAL(triggered()), document, SLOT(Clear()));

   connect(ui->actionSaveState, SIGNAL(triggered()), this, SLOT(onSaveState()));
   connect(ui->actionSaveStateAs, SIGNAL(triggered()), this, SLOT(onSaveStateAs()));
   connect(ui->actionCopyState, SIGNAL(triggered()), this, SLOT(onCopyState()));
   connect(ui->actionExecScriptFromClipboard, SIGNAL(triggered()), this, SLOT(onExecStateFromClipboard()));
   connect(ui->toolButton_ExportCurveData, SIGNAL(clicked()), this, SLOT(onExportCurves()));

   connect(ui->actionSavePicture, SIGNAL(triggered()), this, SLOT(onSavePicture()));
   connect(ui->actionSavePictureAs, SIGNAL(triggered()), this, SLOT(onSavePictureAs()));
   connect(ui->actionCopyPicture, SIGNAL(triggered()), this, SLOT(onCopyPicture()));
   connect(ui->actionQuit, SIGNAL(triggered()), this, SLOT(close()));
   connect(ui->actionOpenPreferences, SIGNAL(triggered()), this, SLOT(onOpenPreferences()));

   connect(ui->actionFlipSelectedOrbital, SIGNAL(triggered()), this, SLOT(onFlipOrbitalClicked()));

   connect(ui->actionShowTablesAndCurves, SIGNAL(triggered()), this, SLOT(ShowTablesAndCurves()));
   ui->actionShowTablesAndCurves->setEnabled(false);
   // ^- FIXME: not really in working condition now... need to clean up and merge first. But
   // other issues are more important for the beginning.

   connect(ui->toolButton_EditFrames, SIGNAL(clicked()), ui->actionEditFrames, SLOT(trigger()));
   connect(ui->actionEditFrames, SIGNAL(triggered()), this, SLOT(ShowEditFramesForm()));
   connect(ui->actionComputeDensityIsoSurface, SIGNAL(triggered()), this, SLOT(ShowEditVolumeSurfacesForm()));
   connect(ui->actionComputeEffectiveOxidationStates, SIGNAL(triggered()), this, SLOT(ShowComputeEosForm()));

   connect(ui->actionFindActivelyReactingOrbitals, SIGNAL(triggered()), this, SLOT(onFindReactingOrbitals()));
   connect(ui->pushButton_FindChangingOrbitals, SIGNAL(clicked()), this, SLOT(onFindReactingOrbitals()));

   connect(ui->actionAlignView_Xy, SIGNAL(triggered()), this, SLOT(onDiscreeteTrafoTriggered()));
   connect(ui->actionAlignView_Xz, SIGNAL(triggered()), this, SLOT(onDiscreeteTrafoTriggered()));
   connect(ui->actionAlignView_Yz, SIGNAL(triggered()), this, SLOT(onDiscreeteTrafoTriggered()));
   connect(ui->actionRotate_90_CW, SIGNAL(triggered()), this, SLOT(onDiscreeteTrafoTriggered()));
   connect(ui->actionRotate_90_CCW, SIGNAL(triggered()), this, SLOT(onDiscreeteTrafoTriggered()));
   connect(ui->actionHideSelectedAtoms, SIGNAL(triggered()), document, SLOT(HideSelectedAtoms()));
   connect(ui->actionResetAtomModes, SIGNAL(triggered()), document, SLOT(ClearAtomFlags()));
   connect(ui->actionResetBondLines, SIGNAL(triggered()), document, SLOT(ResetBondLines()));
   connect(ui->actionLabelElementsOther, SIGNAL(triggered()), ui->checkBox_ElementLabelsOther, SLOT(click()));
   connect(ui->actionLabelElementsC, SIGNAL(triggered()), ui->checkBox_ElementLabelsCarbon, SLOT(click()));
   connect(ui->actionLabelElementsGeneric, SIGNAL(triggered()), this, SLOT(TriggerLabelElementsGeneric()));
   connect(ui->actionLabelAtomNumbers, SIGNAL(triggered()), ui->checkBox_AtomNumbers, SLOT(click()));
   connect(ui->actionShowShadingControls, SIGNAL(triggered()), this, SLOT(toggleShadingControls()));
   QShortcut *sh = new QShortcut(QKeySequence("Ctrl+/"), this); // menu shortcuts don't work for hidden objects :(
   connect(sh, SIGNAL(activated()), this, SLOT(toggleShadingControls()));

   connect(ui->spinBox_OrbitalColorIndex, SIGNAL(valueChanged(int)), document, SLOT(SetNextOrbitalColorIndex(int)));
   connect(document, SIGNAL(NextOrbitalColorIndexChanged(int)), ui->spinBox_OrbitalColorIndex, SLOT(setValue(int)));
   connect(ui->comboBox_OrbitalColorScheme, SIGNAL(currentIndexChanged(int)), document, SLOT(SetNextOrbitalColorScheme(int)));

   connect(ui->actionTraceActiveIsoSurfaces, SIGNAL(triggered()), this, SLOT(onTraceIsoSurfacesClicked()));
   connect(ui->actionTraceAllIsoSurfaces, SIGNAL(triggered()), this, SLOT(onTraceIsoSurfacesClicked()));
   connect(ui->actionComputeWaveFunction, SIGNAL(triggered()), this, SLOT(onComputeWaveFunctionTriggered()));
   connect(ui->actionAbout, SIGNAL(triggered()), this, SLOT(onAboutClicked()));
   connect(ui->toolButton_ComputeWf, SIGNAL(clicked()), ui->actionComputeWaveFunction, SLOT(trigger()));

   // connect action buttons lying around on the UI.
   connect(ui->toolButton_ImportFiles, SIGNAL(clicked()), this, SLOT(onOpenFile()));
   connect(ui->toolButton_SaveState, SIGNAL(clicked()), this, SLOT(onSaveState()));
   connect(ui->toolButton_SaveStateAs, SIGNAL(clicked()), this, SLOT(onSaveStateAs()));
   connect(ui->toolButton_CopyState, SIGNAL(clicked()), this, SLOT(onCopyState()));
//    connect(ui->toolButton_ExecScriptFromClipboard, SIGNAL(clicked()), this, SLOT(onExecStateFromClipboard()));
   connect(ui->toolButton_SavePicture, SIGNAL(clicked()), this, SLOT(onSavePicture()));
   connect(ui->toolButton_SavePictureAs, SIGNAL(clicked()), this, SLOT(onSavePictureAs()));
   connect(ui->toolButton_CopyPicture, SIGNAL(clicked()), this, SLOT(onCopyPicture()));

   connect(ui->toolButton_AlignViewXy, SIGNAL(clicked()), ui->actionAlignView_Xy, SLOT(trigger()));
   connect(ui->toolButton_AlignViewXz, SIGNAL(clicked()), ui->actionAlignView_Xz, SLOT(trigger()));
   connect(ui->toolButton_AlignViewYz, SIGNAL(clicked()), ui->actionAlignView_Yz, SLOT(trigger()));
   connect(ui->toolButton_RotateCcw, SIGNAL(clicked()), ui->actionRotate_90_CCW, SLOT(trigger()));
   connect(ui->toolButton_RotateCw, SIGNAL(clicked()), ui->actionRotate_90_CW, SLOT(trigger()));
   connect(ui->toolButton_FlipOrbital, SIGNAL(clicked()), ui->actionFlipSelectedOrbital, SLOT(trigger()));
   connect(ui->toolButton_TraceAllIsoSurfaces, SIGNAL(clicked()), ui->actionTraceAllIsoSurfaces, SLOT(trigger()));
   connect(ui->toolButton_TraceSelectedIsoSurface, SIGNAL(clicked()), ui->actionTraceActiveIsoSurfaces, SLOT(trigger()));
   connect(ui->toolButton_ReloadShaders, SIGNAL(clicked()), view3d, SLOT(RehashShaders()));
   connect(ui->toolButton_ReloadShaders, SIGNAL(clicked()), view3d, SLOT(update()));
   connect(ui->toolButton_SetShaderPath, SIGNAL(clicked()), this, SLOT(setShaderPath()));

   connect(ui->toolButton_FrameLeft5, SIGNAL(clicked()), this, SLOT(onMoveFrameButtonClicked()));
   connect(ui->toolButton_FrameLeft1, SIGNAL(clicked()), this, SLOT(onMoveFrameButtonClicked()));
   connect(ui->toolButton_FrameRight1, SIGNAL(clicked()), this, SLOT(onMoveFrameButtonClicked()));
   connect(ui->toolButton_FrameRight5, SIGNAL(clicked()), this, SLOT(onMoveFrameButtonClicked()));

   // controls for iso surface changes
   connect(ui->toolButton_RebuildSurface, SIGNAL(clicked()), this, SLOT(onRebuildIsoSurfacesClicked()));
   // note: other controls are auto-synchronized over property model (everything linked via FView3d
   //       property model is)

   connect(ui->toolButton_FrameLog, SIGNAL(clicked()), this, SLOT(showFrameLog()));

   // link up controls for changing the bond & atom settings
   connect(ui->pushButton_ToggleTrackedOrbital, SIGNAL(clicked()), this, SLOT(onToggleTrackedOrbitalClicked()));

   // link up controls for changing orbital rendering settings.
   connect(ui->dialCentralHue, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_AlphaCenter, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_ValCenter, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_SatCenter, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_HuePlus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_HueMinus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_AlphaPlus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_AlphaMinus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_SatPlus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_SatMinus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_ValPlus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));
   connect(ui->spinBox_ValMinus, SIGNAL(valueChanged(int)), this, SLOT(onOrbitalColorChanged(int)));

   // set up a notification for *this if the active data set is changed.
   connect(document, SIGNAL(ActiveDatasetChanged()), this, SLOT(onActiveDatasetChanged()));
   connect(document, SIGNAL(ActiveDatasetChanged()), this, SLOT(onActiveIboCurveChanged()));

   // link up controls for changing the current frame/tracking orbital
   connect(document, SIGNAL(ActiveColChanged(int)), this, SLOT(onFrameIdChanged(int)));

   connect(ui->spinBox_FrameId, SIGNAL(valueChanged(int)), document, SLOT(SetActiveCol(int)));
   connect(ui->dial_FrameId, SIGNAL(valueChanged(int)), document, SLOT(SetActiveCol(int)));
   connect(ui->dial_OrbitalId, SIGNAL(valueChanged(int)), document, SLOT(SetActiveRow(int)));

   ui->checkBox_SkipVirtualOrbitals->setChecked(document->GetSkipVirtualOrbitals());
   connect(ui->checkBox_SkipVirtualOrbitals, SIGNAL(toggled(bool)), document, SLOT(SetSkipVirtualOrbitals(bool)));
   connect(document, SIGNAL(SkipVirtualOrbitalsChanged(bool)), ui->checkBox_SkipVirtualOrbitals, SLOT(setChecked(bool)));


   ui->menuView->addSeparator();
   AddPresetScript(ui->menuView, ":/resources/preset_not_very_shiny.js");
   AddPresetScript(ui->menuView, ":/resources/preset_medium_shiny.js");
   AddPresetScript(ui->menuView, ":/resources/preset_extra_shiny.js");
   AddPresetScript(ui->menuView, ":/resources/preset_sooooo_shiny.js");
   AddPresetScript(ui->menuView, ":/resources/preset_shiny_chic_21.js");
   AssignScriptFileToAction(ui->actionReloadFiles, ":/resources/script_reload_files.js");


   ui->tab_OrbitalColor->setEnabled(false);
   ui->comboBox_IsoType->addItem("density");
   ui->comboBox_IsoType->addItem("absolute");
   ui->comboBox_IsoType->setCurrentIndex(0);

   ui->toolBox_Main->adjustSize();

   LinkPropertyWidgets(view3d, this, "view_option_name");

   connect(ui->dialBondScale, SIGNAL(valueChanged(int)), this, SLOT(UpdateAtomAndBondScalesText(void)));
   connect(ui->dialAtomScale, SIGNAL(valueChanged(int)), this, SLOT(UpdateAtomAndBondScalesText(void)));

   if (g_WorkAroundAlphaCompositing) {
      // it's a work-around for Mac. See comments on g_WorkAroundAlphaCompositing.
      view3d->SetSaveAlpha(false);
      ui->pushButton_SaveAlphaChannel->setEnabled(false);
   }


   {
      // note: should probably code this similarly to the generic ConnectPropertyNotify,
      // which links up the signals/lots directly.
      QList<QAction*> uiActions = this->findChildren<QAction*>();
      foreach(QAction *pAction, uiActions) {
         QVariant
            vaViewOptionName = pAction->property("view_option_name");
         if (vaViewOptionName.isValid()) {
            QString
               viewOptionName = vaViewOptionName.toString();
//             IvEmit("linked up pAction: %1 // %s", pAction->objectName(), viewOptionName);
            QVariant
               property = view3d->property(viewOptionName.toLocal8Bit());
            if (!property.isValid()) {
               IvEmit("!WARNING: link %2 via %1 broken since target property does not exist.", pAction->objectName(), viewOptionName);
               continue;
            }
            // sync up current value of action's checked state to actual property value
            pAction->setChecked(property.toBool());

            // connect action to view property change
            QObject::connect(pAction, SIGNAL(toggled(bool)), this, SLOT(toggleViewPropertyViaAction(bool)));

//             ConnectPropertyNotify(view3d, vaViewOptionName.toString(), pAction, SLOT(setChecked(bool)));
            // connect view property change to action
            ConnectPropertyNotify(view3d, vaViewOptionName.toString(), pAction, "setChecked(bool)");
         }
      }
   }



   // something in the document changed which just requires a new 3d-rendering, but not
   // UI control updates.
   connect(document, SIGNAL(VisualRepresentationChanged()), view3d, SLOT(update()));

   MakeAtomGroupShortcuts();
   setFocusPolicy(Qt::StrongFocus); // <- we need this to be able to get key stroke input

   setAcceptDrops(true);
   onDataChanged();
}


void FMainWindow::TriggerLabelElementsGeneric()
{
   // change labels for both Cs and other atoms at the same time.
   // As indicator, we use the current information from the non-C atoms.
   bool
      OldSetting = view3d->GetLabelElements(),
      NewSetting = !OldSetting;
   view3d->SetLabelElements(NewSetting);
   view3d->SetLabelElementsC(NewSetting);
//    this->update_views();
   view3d->update();

}



void FMainWindow::toggleViewPropertyViaAction(bool newState)
{
   // triggered on QAction events which have a "view_option_name" field.
   QAction
      *pAction= qobject_cast<QAction*>(sender());
   if (pAction) {
         QVariant
            vaViewOptionName = pAction->property("view_option_name");
         if (vaViewOptionName.isValid()) {
            QByteArray
               propertyName = vaViewOptionName.toString().toLocal8Bit();
            if (!propertyName.isEmpty()) {
//                IvEmit("!set view property via QAction handler: %1 to '%2'", QString(propertyName), QString(newState? "true" : "false"));
               QVariant
                  currentState = view3d->property(propertyName.data());
               if (!currentState.isValid() || currentState.toBool() == newState) {
//                   IvEmit("!update skipped: valid = %1  old-value = %2  new-value = %3", currentState.isValid(), currentState.toBool(), newState);
               }

               view3d->setProperty(propertyName.data(), QVariant(newState));
            }
         }
   }
}

int FMainWindow::iGetActionAtomGroupId()
{
   QAction
      *pAction = qobject_cast<QAction*>(sender());
   if (pAction) {
      bool
         ok;
      int
         iGroup = pAction->data().toInt(&ok);
      if (!ok) {
         QMessageBox::warning(this, "Internal UI problem", "Attempted get atom group id via QAction, but QAction contained no group id data.");
         return -1;
      } else
         return iGroup;
   }
   return -1;
}

void FMainWindow::DefineAtomGroup()
{
   int
      iGroup = iGetActionAtomGroupId();
//    IvEmit("DefineAtomGroup %1", iGroup);
   if (iGroup != -1) {
      document->DefineSelectionAsAtomGroup(iGroup);
   }
}

void FMainWindow::SelectAtomGroup()
{
   int
      iGroup = iGetActionAtomGroupId();
   if (iGroup != -1)
      document->SelectAtomGroup(iGroup);
}

void FMainWindow::MakeAtomGroupShortcuts()
{
   // make key combinations Ctrl+0-9 / 0-9 to store and retrieve atom selections.
   for (int iAtomGroup = 0; iAtomGroup <= 9; ++ iAtomGroup) {
      QAction
         *p;
      p = new QAction(QString("Select atom group %1").arg(iAtomGroup), this);
      p->setShortcut(QKeySequence(QString("%1").arg(iAtomGroup)));
      p->setData(QVariant(iAtomGroup));
      connect(p, SIGNAL(triggered()), this, SLOT(SelectAtomGroup()));
      addAction(p); // <- add to main window's event loop.

      p = new QAction(QString("Define atom group %1").arg(iAtomGroup), this);
      p->setShortcut(QKeySequence(QString("Ctrl+%1").arg(iAtomGroup)));
      p->setData(QVariant(iAtomGroup));
      connect(p, SIGNAL(triggered()), this, SLOT(DefineAtomGroup()));
      addAction(p);
   }

   // Note: to *clear* an atom group, create an empty selection (e.g.,
   // by single clicking on an empty space next to the molecule), and
   // press Ctrl+#n (with n=0,1,...).
   // Probably also should add a menu entry for this, or something...
}



void FMainWindow::AssignScriptFileToAction(QAction *pAction, QString FileName)
{
   QString
      ScriptText = LoadTextFromFile(FileName);
   connect(pAction, SIGNAL(triggered()), this, SLOT(ExecPresetScript()));
   pAction->setData(QVariant(ScriptText));
}


void FMainWindow::AddPresetScript(QMenu *pMenu, QString FileName)
{
   QString
      ScriptText = LoadTextFromFile(FileName),
      MenuTitle = FileName;
//    IvEmit("Read script '%1': Contents:\n---%2\n---", FileName, ScriptText);
   QStringList
      Lines = ScriptText.split('\n');
   if (!Lines.isEmpty() && Lines[0].startsWith("// "))
      MenuTitle = Lines[0].midRef(3).toString();
   QAction
      *pAction = new QAction(MenuTitle, this);
   connect(pAction, SIGNAL(triggered()), this, SLOT(ExecPresetScript()));
   pAction->setData(QVariant(ScriptText));
   pMenu->addAction(pAction);
}

void FMainWindow::ExecPresetScript()
{
   QAction
      *pAction = qobject_cast<QAction*>(sender());
   if (pAction) {
      QString ScriptText = pAction->data().toString();
      if (ScriptText == "")
         QMessageBox::warning(this, "Script Execution Failed", "Attempted to execute preset script via menu action, but menu action contained no data.");
      else
         ExecScript(ThisAsIApp(), this->view3d, ScriptText, pAction->text());
   }
}

void FMainWindow::showFrameLog()
{
   FFrame
      *pFrame = document->GetCurrentFrame();
   if (pFrame) {
      FShowTextForm
         ShowTextForm(
            IvFmt("<html><pre>%1</pre></html>", pFrame->Log().GetText()),
            IvFmt("Frame Log for '%1'",RemovePath(pFrame->GetFullInputFileName())),
            document->GetLogFileName(QString("frame%1").arg(int(document->iFrameId(pFrame)), 4, 10, QLatin1Char('0')))
         );
      ShowTextForm.exec();
   }
}


void FMainWindow::onMoveFrameButtonClicked()
{
   // this event is linked to all four of the frame movement buttons on
   // the main form (<<. <, >, >>). The actual amount of frames to move
   // is stored as an integer in the dynamic property "move_frame_button_amount"
   // (-10, -1, +1, +10, currently)
   if (sender() && document->GetCurrentFrame() != 0) {
      QVariant vaAmount = sender()->property("move_frame_button_amount");
      int delta = vaAmount.toInt();
      document->MoveActiveCol(delta);
   }
}


void FMainWindow::onAboutClicked()
{
   FAboutDialog
      about(this);
   about.exec();
}

void FMainWindow::onOpenPreferences()
{
   FPreferencesForm
      PreferencesForm(document, this);
   PreferencesForm.exec();
}

FMainWindow::~FMainWindow()
{
   g_pMainWindow = 0;
      delete document; // I guess this should go down before the 3d view... GL explodes if one destroys objects after the context is gone.
   delete view3d;
   delete dataViewFilter;
   delete ui;
}



bool FShowOnlyCurrentColumnFilter::filterAcceptsColumn(int sourceColumn, const QModelIndex &/*sourceParent*/) const
{
//    std::cout << boost::format("invoked filter model for row %i.") % sourceColumn << std::endl;
   return sourceColumn == m_ShownColumn;
}

FShowOnlyCurrentColumnFilter::FShowOnlyCurrentColumnFilter(int ShownColumn, QObject *parent)
   : FBase(parent)
{
   m_ShownColumn = ShownColumn;
}

void FShowOnlyCurrentColumnFilter::setShownColumn(int NewColumn)
{
   if (NewColumn != m_ShownColumn) {
      m_ShownColumn = NewColumn;
      invalidateFilter();
   }
}


float SvBoxScale = 100.f;

float clamp(float f, float min, float max) {
   if (f < min) return min;
   if (f > max) return max;
   return f;
}

float AlphaFromCol(uint32_t c) {
   return ((c >> 24) & 0xff)/255.f;
}

uint32_t ColFromAlpha(float a) {
   int ia = clamp(::roundf(255 * a), 0., 255.);
   return ia << 24;
}

void FMainWindow::onActiveDatasetChanged()
{
//    std::cout << boost::format("ActiveDatasetChanged:  current row is %i.") % document->GetActiveRowIndex() << std::endl;
   // we'll be updating the control values. this will emit valueChanged signals,
   // which will trigger orbital updates before we are finished. stop that.
   iUpdateLocked += 1;

   FDataSetPtr
      pDataSet = document->GetActiveDataSet();
   FVolumeDataSet
      *pOrb = dynamic_cast<FVolumeDataSet*>(pDataSet.get());

   if (document->GetCurrentFrame()) {
      int nOrbs = int(document->GetCurrentFrame()->m_Data.size() - 1);
      ui->dial_OrbitalId->setEnabled(nOrbs > 0);
      if (nOrbs != 0) { // <- without this we may get an infinite loop between the slider's valueChanged and the document's ActiveColChanged...
         ui->dial_OrbitalId->setMinimum(1);
         ui->dial_OrbitalId->setMaximum(std::max(nOrbs, 1));
         if (pOrb)
            ui->dial_OrbitalId->setValue(document->GetActiveRowIndex());
      }
   }

   bool
      EnableTwoPhaseColorTab = false;
   if (pOrb && pOrb->Active) {
      // update the color values.
//       IvEmit("!Set Active Orbital: %1", pOrb->GetDesc().toStdString());
      FVolumeVisualConfig
         *pVis = pOrb->pVisConfig.get();
      if (!pVis->DetailsAssigned())
         view3d->updateGL(); // <- to force rendering the orbital such that we have a color set in the curve view...
      FTwoPhaseIsoSurfaceConfig const
         *pVisDet = dynamic_cast<FTwoPhaseIsoSurfaceConfig const*>(pOrb->pVisConfig->pDetails.get());
      if (pVisDet) {
         EnableTwoPhaseColorTab = true;
         FColor
            cPlus(pVisDet->cIsoPlus),
            cMinus(pVisDet->cIsoMinus),
            cAvg;
         cAvg = .5f * cPlus + .5f * cMinus;
         float
            hp,sp,vp,ap,
            hm,sm,vm,am,
            hc,sc,vc,ac;
         ac = 1.0;
         cPlus.ToHsv(hp,sp,vp);
         cMinus.ToHsv(hm,sm,vm);
         cAvg.ToHsv(hc,sc,vc);
         ap = AlphaFromCol(pVisDet->cIsoPlus);
         am = AlphaFromCol(pVisDet->cIsoMinus);
         ac = .5f*(ap + am);
         if (::fabs(hp-hm) > 180.) {
            if (hp > hm)
               hp -= 360.;
            else
               hm -= 360.;
         }

         hc = .5f * (hp + hm); // that is not right, is it?
         sc = .5f * (sp + sm);
         vc = .5f * (vp + vm);

   //          std::cout << boost::format("hc = %6.2f  sc = %6.2f  vc = %6.2f  ac = %6.2f") % hc % (SvBoxScale*sc) % (SvBoxScale*vc) % (SvBoxScale*ac) << std::endl;
         ui->dialCentralHue->setValue((int)hc);
         ui->label_Hue->setText(QString("Hue (%1)").arg((int)hc));
         ui->spinBox_SatCenter->setValue(SvBoxScale * sc);
         ui->spinBox_ValCenter->setValue(SvBoxScale * vc);
         ui->spinBox_AlphaCenter->setValue(SvBoxScale * ac); // not actually changed atm.

         ui->spinBox_HuePlus->setValue(hp - hc);
         ui->spinBox_SatPlus->setValue(SvBoxScale * (sp - sc));
         ui->spinBox_ValPlus->setValue(SvBoxScale * (vp - vc));
         ui->spinBox_HueMinus->setValue(hm - hc);
         ui->spinBox_SatMinus->setValue(SvBoxScale * (sm - sc));
         ui->spinBox_ValMinus->setValue(SvBoxScale * (vm - vc));
      }
   }
//    ui->tab_OrbitalColor->setEnabled(pOrb && pOrb->Active);
   ui->tab_OrbitalColor->setEnabled(EnableTwoPhaseColorTab);

   iUpdateLocked -= 1;
}

void FMainWindow::onOrbitalColorChanged(int /*iNewValue*/)
{
   // ^- can't use iNewValue... we don't actually know where it is coming from.
   if (iUpdateLocked != 0)
      // currently in the process of updating controls; that
      // is why the values are changing. don't do anything until
      // this is finished.
      return;
   FDataSetPtr
      pDataSet = document->GetActiveDataSet();
   if (pDataSet.get() && pDataSet->Active) {
      FVolumeDataSet
         *pOrb = dynamic_cast<FVolumeDataSet*>(pDataSet.get());
      FTwoPhaseIsoSurfaceConfig
         *pVisDet = dynamic_cast<FTwoPhaseIsoSurfaceConfig*>(pOrb->pVisConfig->pDetails.get());
      if (pVisDet) {
         // assemble new colors from data in controls.
         float
            hc = ui->dialCentralHue->value(),
            vc = ui->spinBox_ValCenter->value()/SvBoxScale,
            sc = ui->spinBox_SatCenter->value()/SvBoxScale,
            ac = ui->spinBox_AlphaCenter->value()/SvBoxScale;
         float
            hp = hc + ui->spinBox_HuePlus->value(),
            sp = clamp(sc + ui->spinBox_SatPlus->value()/SvBoxScale, 0.f, 1.f),
            vp = clamp(vc + ui->spinBox_ValPlus->value()/SvBoxScale, 0.f, 1.f),
            ap = clamp(ac + ui->spinBox_AlphaPlus->value()/SvBoxScale, 0.f, 1.f),
            hm = hc + ui->spinBox_HueMinus->value(),
            sm = clamp(sc + ui->spinBox_SatMinus->value()/SvBoxScale, 0.f, 1.f),
            vm = clamp(vc + ui->spinBox_ValMinus->value()/SvBoxScale, 0.f, 1.f),
            am = clamp(ac + ui->spinBox_AlphaMinus->value()/SvBoxScale, 0.f, 1.f);

         ui->label_Hue->setText(QString("Hue (%1)").arg((int)hc));

//          uint32_t
//             cPlus = ct::Hsv(hp,sp,vp).uint32() | 0xff000000,
//             cMinus = ct::Hsv(hm,sm,vm).uint32() | 0xff000000;
         uint32_t
            cPlus = ct::Hsv(hp,sp,vp).uint32() | ColFromAlpha(ap),
            cMinus = ct::Hsv(hm,sm,vm).uint32() | ColFromAlpha(am);
         FVolumeVisualConfig
            *pVis = pOrb->pVisConfig.get();
//          pVis->AssignPlusMinusColors(cPlus, cMinus, -1.);
         pVisDet->cIsoPlus = cPlus;
         pVisDet->cIsoMinus = cMinus;
//          pVisDet->bColorSet = true; // otherwise ellipsoid colors cannot be updated.
         // update the color also in the mesh, if one already exists.
         pVis->UpdateLinkedRepresentations(FVolumeVisualConfig::UPDATE_InvalidateColors, this->view3d);
         view3d->update();
      }
   }
}


void FMainWindow::toggleShadingControls()
{
   ui->groupBox_ShadingControls->setVisible(!ui->groupBox_ShadingControls->isVisible());
//    ui->groupBox_ShadingControls->setVisible(ui->actionShowShadingControls->isChecked());
}

void FMainWindow::setShaderPath()
{
   QString
      NewPath = QFileDialog::getExistingDirectory(this, "Open Path");
   if (NewPath != "") {
      view3d->SetShaderPath(NewPath);
      view3d->RehashShaders();
      view3d->update();
   }
}


void FMainWindow::onRebuildIsoSurfacesClicked()
{
   double
      fNewIsoValue = view3d->m_IsoThreshold/100.;
//       fNewIsoResolution = view3d->m_IsoResolution;
//    bool
//       bNewIsoAutoFlip = view3d->m_IsoAutoFlip;
   for (int iFrame = 0; iFrame < document->GetNumFrames(); ++ iFrame) {
      FDataSetList
         *pFrameData = document->GetFrameData(iFrame);
      if (pFrameData) {
         for (int iRow = 0; size_t(iRow) < pFrameData->size(); ++ iRow) {
            FDataSet
               *pDataSet = &*((*pFrameData)[iRow]);
            if (pDataSet == 0)
               continue;
            FOrbital
               *pOrbital = dynamic_cast<FOrbital*>(pDataSet);
            if (pOrbital != 0) {
               // this is an orbital data set. update its iso threshold value
               // and ask it to rebuild its mesh.
               FVolumeVisualConfig
                  *pVisConfig = pOrbital->pVisConfig.get();
               FTwoPhaseIsoSurfaceConfig
                  *pTwoPhaseConfig = dynamic_cast<FTwoPhaseIsoSurfaceConfig*>(pVisConfig->pDetails.get());
               if (pTwoPhaseConfig != 0) {
                  pTwoPhaseConfig->fIsoValue = fNewIsoValue;
//                   pTwoPhaseConfig->fIsoValue = fNewIsoValue;
               }
            }

            pDataSet->InvalidateRenderCache();
         }
      }
   }
   view3d->update();
}

void FMainWindow::onTraceIsoSurfacesClicked()
{
   int
      nTotalSets = 0;
   bool
      TraceAll = false;
   if (sender() == ui->actionTraceAllIsoSurfaces)
      TraceAll = true;

   if (!TraceAll) {
      // count number of active data sets, to get an estimation of the required time.
      int
         nActiveSets = 0; // note: one of them is the geometry...
      for (int iRow = 0; size_t(iRow) < document->GetCurrentFrameData()->size(); ++ iRow)
         if ((*document->GetCurrentFrameData())[iRow]->Active)
            nActiveSets += 1;
      nTotalSets = nActiveSets * document->GetNumFrames();
   } else {
      nTotalSets = document->GetCurrentFrameData()->size() * document->GetNumFrames();
   }

   QProgressDialog progress("Tracing iso-surfaces...", "Cancel", 0, nTotalSets, this);
   progress.setWindowModality(Qt::WindowModal);
   progress.setValue(0);

   int
      nDoneSets = 0;
   bool
      Cancelled = false;
   for (int iFrame = 0; iFrame < document->GetNumFrames(); ++ iFrame) {
      FDataSetList
         *pData = document->GetFrameData(iFrame);
      if (pData) {
         for (int iRow = 0; size_t(iRow) < pData->size(); ++ iRow) {
            if ((*pData)[iRow]->Active || TraceAll) {
               for (int iDummy = 0; iDummy < 100; ++ iDummy)
                  QApplication::processEvents();
               // ^- this shouldn't be here, but otherwise the "Cancel" button does not work...
               // update... well, it also doesn't work *with* it behing here..
               Cancelled = Cancelled || bool(progress.wasCanceled());
               if (Cancelled)
                  break;
               (*pData)[iRow]->BuildRenderCache(view3d);
               nDoneSets += 1;
               progress.setValue(nDoneSets);
               IvEmit("* progess: %1 of %2 iso-surfaces done.", nDoneSets, nTotalSets);
            }
         }
      }
   }
   if (Cancelled)
      std::cout << " Iso Tracing cancelled." << std::endl;
   progress.setValue(nTotalSets);
   view3d->update();
}

void FMainWindow::onFlipOrbitalClicked()
{
   if (iUpdateLocked != 0)
      // currently in the process of updating controls; that
      // is why the values are changing. don't do anything until
      // this is finished.
      return;
   FDataSetPtr
      pDataSet = document->GetActiveDataSet();
   if (pDataSet.get() && pDataSet->Active) {
      FOrbital
         *pOrb = dynamic_cast<FOrbital*>(pDataSet.get());
      if (pOrb) {
         FVolumeVisualConfig
            *pVis = pOrb->pVisConfig.get();
         FTwoPhaseIsoSurfaceConfig
            *pVisDet = dynamic_cast<FTwoPhaseIsoSurfaceConfig*>(pVis->pDetails.get());
         if (pVisDet) {
   //          pOrb->FlipPhase();
            std::swap(pVisDet->cIsoMinus, pVisDet->cIsoPlus);
//             pVis->AssignPlusMinusColors(pVis->cIsoMinus(), pVis->cIsoPlus(), -1.);
            pVis->UpdateLinkedRepresentations(FVolumeVisualConfig::UPDATE_InvalidateColors, this->view3d);

            onActiveDatasetChanged();
            // ^- to update the values in the color controls.
            view3d->update();
         }
      }
   }
}

void FMainWindow::onTitleChanged()
{
   QString
      FileName = RemovePath(document->GetCurrentInputFileName());
   if (!FileName.isEmpty())
      setWindowTitle(FileName + QString(" -- IboView "));
   else
      setWindowTitle("IboView");
}


void FMainWindow::onFrameIdChanged(int iNewFrame)
{
//    document->SetActiveCol(iNewFrame);
   dataViewFilter->setShownColumn(iNewFrame);
   ui->datasetView->resizeColumnsToContents();
   if (document->GetCurrentFrame() != 0) {
      QString const
         &FileName = document->GetCurrentFrame()->GetBaseInputFileName(),
         LabelText = IvFmt("#%1: %2", iNewFrame, FileName);
      ui->label_FrameFileName->setText(LabelText);
   }
}

void FMainWindow::onDataRowClicked(QModelIndex const &Index)
{
   document->SetActiveRow(Index.row());
   // ^- note that these are indices from 'dataViewFilter', not from 'document' itself.
}

void FMainWindow::onDataRowDoubleClicked(QModelIndex const &Index)
{
   document->ToggleDataRow(Index.row());
}

void FMainWindow::onDataChanged()
{
   int
      iFrame = document->GetActiveColIndex(),
      nFrames = document->GetNumFrames();
   //if (iFrame == -1)
   //   iFrame = 0;
   //if (nFrames == 0)
   //   nFrames = 1;
   // ^- FIXME: is that safe? can't this make an infinite loop?
   if (nFrames > 0) {
      ui->spinBox_FrameId->setMinimum(0);
      ui->dial_FrameId->setMinimum(0);

      ui->spinBox_FrameId->setMaximum(nFrames - 1);
      ui->dial_FrameId->setMaximum(nFrames - 1);

      if (iFrame != -1) {
         ui->spinBox_FrameId->setValue(iFrame);
         ui->dial_FrameId->setValue(iFrame);
      }
   }

   ui->tab_Frames->setEnabled(nFrames >= 2);
   ui->spinBox_FrameId->setEnabled(nFrames >= 2);

   // todo: update file names etc.
}


void FMainWindow::onSaveState(/*const QAction &action*/)
{
//    std::cout << "onSaveState triggered!!" << std::endl;
   WriteStateScript();
}

void FMainWindow::onCopyState(/*const QAction &action*/)
{
   WriteStateScript(":/!clipboard!");
}

// given a text, check if it kinda looks like an xyz file (admittedly, a better
// method might be to simply try to load it as an .xyz file, but this might get
// more messy).
bool DoesThisLookLikeAnXyzFile(QString Text)
{
//    IvEmit("Got into 'DoesThisLookLikeAnXyzFile(..).'");
   // check if there are at least three lines (#Atoms, Comment, First Element-X-Y-Z-line)
   QStringList
      Lines = Text.split("\n");
   if (Lines.size() < 3)
      return false;
//    IvEmit("...line number check passed.'");

   // Does the first line consists of a single integer?
   QString Line0 = Lines[0];
   bool
      ok[4] = {false,false,false,false};
   int
      nAtoms = Line0.toInt(&ok[0],10);
   // ^-  10: conversion base. sorry, no octal number of atoms for you 8).
   // (might actually lead to problems for people who prepend their ints with 0s)
   IR_SUPPRESS_UNUSED_WARNING(nAtoms);

//    if (!ok[0] || Lines.size() < 2 + nAtoms)
   // ^- actual load attempt would give more meaningfull error messages than checking consistency here.
   if (!ok[0])
      return false;
//    IvEmit("...first-line-integer check passed.'");

   // now see if there are three doubles in line 2.
   QString
      Line2 = Lines[2];
   QStringList
      ls = Line2.trimmed().split(QRegExp("\\s+")); // <- split at whitespace
   if (ls.size() >= 4) { // element x y z
//       IvEmit("...3rd line count check check passed.'");
//       IvEmit("...ls[1] = '%1'", ls[1]);
//       IvEmit("...ls[2] = '%1'", ls[2]);
//       IvEmit("...ls[3] = '%1'", ls[3]);
      ls[1].toDouble(&ok[1]);
      ls[2].toDouble(&ok[2]);
      ls[3].toDouble(&ok[3]);
   }
//    IvEmit("...stabs: %1 %2 %3 %4'", int(ok[0]), int(ok[1]), int(ok[2]), int(ok[3]));
   return ok[0] && ok[1] && ok[2] && ok[3];
}



void FMainWindow::onExecStateFromClipboard()
{
   QClipboard
      *clipboard = QApplication::clipboard();
   QString
      ScriptText = clipboard->text();
   if (ScriptText == "")
      QMessageBox::warning(this, "Script Execution Failed", "Attempted to execute script from clipboard. but clipboard does not contain text.");
   else {
      if (DoesThisLookLikeAnXyzFile(ScriptText)) {
         // try to load it as an .xyz file. Unload other stuff first.
         document->Clear();
         ThisAsIApp()->load_file(":/!clipboard!");
      } else
         // run it as a script.
         ExecScript(ThisAsIApp(), this->view3d, ScriptText, "clipboard");
   }
}

void FMainWindow::onSavePicture(/*const QAction &action*/)
{
//    view3d->HasFlag(VIEWFLAG_CropImages)
   view3d->save_png(ReplaceExt(document->GetCurrentInputFileName(), ".png"));
}

void FMainWindow::onCopyPicture(/*const QAction &action*/)
{
   view3d->save_png(":/!clipboard!");
}


// FSaveFileDirChangeProxy::FSaveFileDirChangeProxy(QFileDialog *pDialog)
//    : m_pDialog(pDialog)
// {
//    connect(pDialog, SIGNAL(fileSelected(QString const&)), this, SLOT(selectFile(QString const&)));
// }
//
// FSaveFileDirChangeProxy::~FSaveFileDirChangeProxy()
// {}
//
// void FSaveFileDirChangeProxy::selectFile(QString const &File)
// {
//    std::cout << "PROXY: Selected '" << q2s(File) << std::endl;
// }


void FMainWindow::onSavePictureAs(/*const QAction &action*/)
{
   QString
      FileNameSuggest = ReplaceExt(document->GetCurrentInputFileName(), ".png"),
      FileName;
   FileName = IvGetSaveFileName("History/SaveFiles", this, "Save Picture",
         (FileNameSuggest), "Bitmap files (*.png *.jpg);;All Files (*.*)");
//    FileName = QFileDialog::getSaveFileName(this, "Save Picture",
//          ("./"+FileNameSuggest), "Bitmap files (*.png *.jpg);;All Files (*.*)");

   if (FileName != "")
      view3d->save_png(FileName);
}

void FMainWindow::onSaveStateAs(/*const QAction &action*/)
{
   QString
      RefFileName = document->GetCommonInputFileName(),
      FileNameSuggest = ReplaceExt(RefFileName, ".js"),
      FileName;

//    FileName = QFileDialog::getSaveFileName(this, "Save State Script",
//          ("./"+FileNameSuggest), "Script files (*.js *.chai);;All Files (*.*)");
   FileName = IvGetSaveFileName("History/SaveFiles", this, "Save State Script",
         FileNameSuggest, "Script files (*.js *.chai);;All Files (*.*)");


   if (!FileName.isEmpty())
      WriteStateScript(FileName);
}

void FMainWindow::onOpenFile()
{
//    QStringList FileNames = QFileDialog::getOpenFileNames(this, "Open File",
//       "", "Molpro XMLs, Xyz Files, Molden files, IboView Scripts (*.xml *.xyz *.molden *.js *.chai);;Xyz Files (*.xyz);;All Files (*.*)");
   QStringList FileNames = IvGetOpenFileNames("History/OpenFiles", this, "Open File",
      "", "Molpro XMLs, Xyz Files, Molden files, IboView Scripts (*.xml *.xyz *.molden *.molden.input *.js *.chai);;Xyz Files (*.xyz);;All Files (*.*)");

   load_files_(FileNames);
}

void FMainWindow::dragEnterEvent(QDragEnterEvent *event)
{
//    FBase::dragEnterEvent(event);
//    if (event->mimeData()->hasFormat("text/plain"))
   // maybe I should check if these are the correct file types?
   if (event->mimeData()->hasUrls())
      event->acceptProposedAction();
}

void FMainWindow::dropEvent(QDropEvent *event)
{
   const QMimeData *mimeData = event->mimeData();

   // is it a list of files?
   if (mimeData->hasUrls()) {
      QStringList
         pathList;
      QList<QUrl>
         urlList = mimeData->urls();

      // extract the local paths of the files
      for (int i = 0; i < urlList.size(); ++i) {
         pathList.append(urlList.at(i).toLocalFile());
      }

      // call a function to open the files
//       openFiles(pathList);
//       for (int i = 0; i < pathList.size(); ++i)
//          IvEmit("!! OPEN: '%1'", pathList[i]);
      close_files_();
      load_files_(pathList);
   }
//    FBase::dropEvent(event);
}


void FMainWindow::closeEvent(QCloseEvent *event)
{
   StoreApplicationSettings();
   return QMainWindow::closeEvent(event);
}


void FMainWindow::StoreApplicationSettings()
{
   IvSaveWindowSize("MainWindow/Size", this);
   IvSaveSplitterState("MainWindow/SplitterH", ui->splitter_ViewVsTools);
}

void FMainWindow::LoadApplicationSettings()
{
   IvRestoreWindowSize("MainWindow/Size", this);
   IvRestoreSplitterState("MainWindow/SplitterH", ui->splitter_ViewVsTools);
}


void FMainWindow::ShowTablesAndCurves()
{
   // check if we already have one of those dialogs.
   FTablesForm*
      tablesForm = findChild<FTablesForm*>("TablesForm");
   if (!tablesForm) {
      // seems not. Make a new one.
      tablesForm = new FTablesForm(document, this);
   }
   tablesForm->show();
}

void FMainWindow::ShowEditFramesForm()
{
   FEditFramesForm
      *framesForm = findChild<FEditFramesForm*>("EditFramesForm");
   if (!framesForm)
      framesForm = new FEditFramesForm(document, this);
   framesForm->exec();
   framesForm->DoPostProcessingOnClose();
//    framesForm->done(0);
//    // should I call done() here to invoke the close command? Otherwise it's just hidden, I guess..
   framesForm->close();
}

void FMainWindow::ShowEditVolumeSurfacesForm()
{
//    FEditFramesForm
//       *framesForm = findChild<FEditFramesForm*>("EditFramesForm");
//    if (!framesForm)
//       framesForm = new FEditFramesForm(document, this);
//    framesForm->exec();
//    framesForm->DoPostProcessingOnClose();
// //    framesForm->done(0);
// //    // should I call done() here to invoke the close command? Otherwise it's just hidden, I guess..
   FEditVolumeSurfaceForm
      *pVolumeSurfaceForm = findChild<FEditVolumeSurfaceForm*>("EditVolumeSurfacesForm");
   if (pVolumeSurfaceForm != 0)
      IvEmit("hm... volume surfaces form already there. It shouldn't be, right?");
   FEditVolumeSurfaceForm
      EditVolumeSurfaceForm(document, this);
   EditVolumeSurfaceForm.exec();
}

void FMainWindow::ShowComputeEosForm()
{
   FComputeEosForm
      *pComputeEosForm = findChild<FComputeEosForm*>("ComputeEosForm");
   if (pComputeEosForm != 0)
      IvEmit("hm... ComputeEosCharges form already there. It shouldn't be, right?");
   FComputeEosForm
      ComputeEosForm(document, this);
   if (bool(ComputeEosForm.exec())) {
      FShowTextForm
         *pShowTextForm = new FShowTextForm(QString(), "Effective Oxidation State (EOS) Analysis", document->GetLogFileName("eos"), "Abort", this);
      pShowTextForm->setAttribute( Qt::WA_DeleteOnClose, true );
      FMemoryLogQt *pLog = new FMemoryLogQt(pShowTextForm);
      connect(pLog, SIGNAL(textEmitted(QString)), pShowTextForm, SLOT(AppendText(QString const&)), Qt::QueuedConnection);
      connect(pLog, SIGNAL(sectionEnded()), pShowTextForm, SLOT(clear()), Qt::QueuedConnection);
      connect(pLog, SIGNAL(reportEnded(QString)), pShowTextForm, SLOT(setText(QString)), Qt::QueuedConnection);
      connect(pShowTextForm, SIGNAL(redButtonPressed()), pLog, SLOT(emitAbortSignal()), Qt::QueuedConnection);

      document->RunRedoxChargeAnalysis(*pLog);

      pShowTextForm->exec();
   }
}

void FMainWindow::onDiscreeteTrafoTriggered()
{
   QAction
      *pAction = qobject_cast<QAction*>(sender());
   // note: could also use pAction->setData() to discern what is what.
   if (pAction == ui->actionAlignView_Xy) {
      view3d->modify("align xy", 0.);
   } else if (pAction == ui->actionAlignView_Xz) {
      view3d->modify("align xz", 0.);
   } else if (pAction == ui->actionAlignView_Yz) {
      view3d->modify("align yz", 0.);
   } else if (pAction == ui->actionRotate_90_CW) {
      view3d->modify("roll", +90.);
   } else if (pAction == ui->actionRotate_90_CCW) {
      view3d->modify("roll", -90.);
   } else {
      IvNotify(NOTIFY_Warning, IvFmt("WARNTING: Action '%1' not recognized in dtraf. Ignored!", pAction->objectName()));
   }
}



void FMainWindow::onFindReactingOrbitals()
{
   document->FindReactingOrbitals();
}


// FIXME: this stuff should all be merged into IvTables.
void FMainWindow::onExportCurves()
{
   QString
      RefFileName = document->GetCommonInputFileName(),
      FileNameSuggest = ReplaceExt(RefFileName, ".csv"),
      FileName;
//    FileName = QFileDialog::getSaveFileName(this, "Export Curve Data",
//          ("./"+FileNameSuggest), "Comma Separated Values (*.csv);;All Files (*.*)");
   FileName = IvGetSaveFileName("History/SaveCurves", this, "Export Curve Data",
         FileNameSuggest, "Comma Separated Values (*.csv);;All Files (*.*)");
   if (FileName.isEmpty())
      return;
   size_t
      nFrames = document->GetNumFrames();

   TArray<float>
      ArcLengths;
   {
      document->UpdateArcLengths();
      ArcLengths.resize(nFrames);
      for (size_t iFrame = 0; iFrame < nFrames; ++ iFrame)
         ArcLengths[iFrame] = document->GetFrame(iFrame)->GetArcLength().Value();
   }

   QList<TArray<float> >
      CurveData;
   QString
      CaptionLine = QString("\"%1\"").arg(document->GetArcLengthOptions().GetArcCaption());
//    if ((ArcLengthFlags & ARCLENGTH_NoMassWeighting) != 0)
//       CaptionLine = "\"Arc Length (bohr)\"";
//    else
//       CaptionLine = "\"Reaction Coordinate (bohr$\\cdot$amu$^{1/2}$)\"";

   for (int iCurve = 0; iCurve < ui->curveView->getNumCurves(); ++ iCurve) {
      TArray<float>
         Data;
      QString Title;
      QColor Color;
      uint32_t Flags;
      ui->curveView->getCurve(Data, Title, Flags, Color, iCurve);
      if (Title == "" || Title == "(current)")
         continue;
      CurveData.append(Data);
      uint32_t dwColor = ((Color.red() << 16) + (Color.green() << 8) + Color.blue());
      QString Caption = QString::fromStdString(fmt::format("{} [#{:06x}]", Title.toStdString(), dwColor));
      if (!CaptionLine.isEmpty())
         CaptionLine += ",";
      CaptionLine.append(IvFmt("\"%1\"", Caption));
   }
   QFile
      OutFile(FileName);
   if (OutFile.open(QFile::WriteOnly | QFile::Truncate)) {
      QTextStream Out(&OutFile);

      Out << CaptionLine << "\n";
      Out.setRealNumberNotation(QTextStream::FixedNotation);
      Out.setRealNumberPrecision(6);
      for (uint iFrame = 0; iFrame < nFrames; ++ iFrame) {
         Out << ArcLengths[iFrame];
         for (int iCurve = 0; iCurve < CurveData.size(); ++ iCurve) {
//             if (iCurve != 0)
//                Out << ",";
            Out << ",";
            float f = CurveData.at(iCurve)[iFrame];
            Out << f;
         }
         Out << "\n";
      }
      IvEmit("* wrote curve data to '%1'", FileName);
   } else {
      IvNotify(NOTIFY_Error, IvFmt("Failed to open file '%1' for writing.", FileName));
   }
}


template<class FFloat>
static void PrintChargeArray(std::ostream &xout, std::string const &Desc, TArray<FFloat> &CurveData)
{
   xout << fmt::format("!chgs {}  npts = {}  [", Desc, CurveData.size());
   for (uint i = 0; i < CurveData.size(); ++ i) {
      if (i != 0)
         xout << ", ";
      xout << fmt::format("{:5.2f}", CurveData[i]);
   }
   xout << "]\n";
}

static void NormalizeCurve(TArray<float> &CurveData, float fMin, float fMax, bool LowerTo0=false)
{
   if (CurveData.empty())
      return;
   float
      fDataMin = CurveData[0],
      fDataMax = CurveData[0];
   for (size_t i = 1; i < CurveData.size(); ++ i) {
      fDataMax = std::max(float(fDataMax), float(CurveData[i]));
      fDataMin = std::min(float(fDataMin), float(CurveData[i]));
   }
   if (LowerTo0)
      fDataMin = 0.;
   if (fDataMax - fDataMin <= 0) {
      for (size_t i = 0; i < CurveData.size(); ++ i)
         CurveData[i] = (fMax + fMin) / 2;
   } else {
      for (size_t i = 0; i < CurveData.size(); ++ i)
         CurveData[i] =  ((CurveData[i] - fDataMin)/(fDataMax - fDataMin)) * (fMax - fMin) + fMin;
   }
}


// FIXME: this stuff should all be merged into IvTables.
void FMainWindow::onActiveIboCurveChanged()
{
   if (document->GetActiveRowIndex() == 0)
      return;
//    std::cout << "onActiveIboCurveChanged [iUpd = " << iUpdateLocked << "]" << std::endl;
   if (iUpdateLocked != 0)
      return;
   if (!document->GetFrameData(0) || !document->GetCurrentFrameData())
      return;
   iUpdateLocked += 1;
   FDataSetList
      &DataList = *document->GetFrameData(0);
//    ui->dial_OrbitalId->setMinimum(1);
//    ui->dial_OrbitalId->setMaximum(int(DataList.size())-1);

//    uint32_t
//       dwEnabledColor = 0xffa0a0a0, // color of data sets which are enabled but not current (or use active color?)
//       dwActiveColor = 0xffffffff; // color of current data set
   TArray<float>
      CurveData;
   uint
      nCurves = 0,
      iFrame = document->GetActiveColIndex();
   ui->curveView->LockUpdates();
   ui->curveView->setBackgroundBrush(QBrush(QColor(0xff404040)));
   ui->curveView->clearCurves();
//    Qt::DotLine
   QPen AxisPen = QPen(QBrush(QColor(0xff171717)), 3.41, Qt::SolidLine);
   ui->curveView->addHline(float(0.f), AxisPen, HVCURVE_EndArrow);
   ui->curveView->addVline(float(0.f), AxisPen, HVCURVE_EndArrow);
   uint32_t
      dwFrameIndicatorColor = 0xff70e070;
//       dwFrameIndicatorColor = 0xffff6060;
   ui->curveView->addVline(float(iFrame), QPen(QBrush(QColor(dwFrameIndicatorColor)), 2.41, Qt::DotLine));
   if (document->HaveEnergies() && ui->pushButton_ShowCurveEnergy->isChecked()) {
      CurveData.resize(document->GetNumFrames());
      for (int iFrame = 0; iFrame < document->GetNumFrames(); ++ iFrame)
         CurveData[iFrame] = float(document->GetFrame(iFrame)->GetEnergy());
      NormalizeCurve(CurveData, 0., 1.);
      uint iCurveId = nCurves + 1;
      uint32_t
         dwColor = 0xffc0c0,
         dwFlags = CURVE_Outline | CURVE_Ellipses;
      QPen pen = QPen(QBrush(QColor(dwColor)), 1.4, Qt::SolidLine); // Qt::DashLine
      ui->curveView->setCurve(iCurveId, CurveData, "Energy", dwFlags, -1, pen, QBrush(dwColor));
      nCurves += 1;
   }
   if (document->HaveGradients() && ui->pushButton_ShowCurveGradient->isChecked()) {
      CurveData.resize(document->GetNumFrames());
      for (int iFrame = 0; iFrame < document->GetNumFrames(); ++ iFrame)
         CurveData[iFrame] = float(document->GetFrame(iFrame)->GetGradient());
      for (size_t i = 0; i < CurveData.size(); ++ i) {
         if (CurveData[i] < 1e-10)
            CurveData[i] = 1e-10;
         CurveData[i] = std::log(CurveData[i]);
      }
      NormalizeCurve(CurveData, 0., 1.);
      uint iCurveId = nCurves + 1;
      uint32_t
         dwColor = 0xc0c0ff,
         dwFlags = CURVE_Outline | CURVE_Ellipses;
      QPen pen = QPen(QBrush(QColor(dwColor)), 1.4, Qt::SolidLine); // Qt::DashLine
      ui->curveView->setCurve(iCurveId, CurveData, "Gradient", dwFlags, -1, pen, QBrush(dwColor));
      nCurves += 1;
   }

   for (uint iDataSet = 0; iDataSet < DataList.size(); ++ iDataSet)
   {
      FDataSetPtr
         pActiveData = DataList[iDataSet];
      bool
         IsCurrent_ = int(iDataSet) == document->GetActiveRowIndex();
//          IsCurrent_ = int(iDataSet) == ui->dial_OrbitalId->value();
      if (!pActiveData->Active && !IsCurrent_)
         continue;
      FOrbital
         *pOrbital = dynamic_cast<FOrbital*>(pActiveData.get());
      if (!pOrbital)
         continue;
      if (IsCurrent_) {
         ui->pushButton_ToggleTrackedOrbital->setChecked(pActiveData->Active);
         ui->pushButton_ToggleTrackedOrbital->setText(IvFmt("Show #%1", iDataSet));
      }
      for (uint IsCurrent = 1 - uint(pActiveData->Active); IsCurrent != 1 + uint(IsCurrent_); ++ IsCurrent) {
         QPen
            pen;
         uint32_t
            dwColor,
            dwFlags = CURVE_Outline | CURVE_Ellipses;
//          dwFlags |= CURVE_Bold;
         int
            zOrder = 0;
         uint
            iCurveId = nCurves + 1;
         if (IsCurrent == 1) {
            dwColor = 0xffffffff;
//             pen = QPen(QBrush(QColor(dwColor)), 3.41, Qt::DotLine); // Qt::DashLine
            pen = QPen(QBrush(QColor(dwColor)), 2.41, Qt::DotLine); // Qt::DashLine
            zOrder = 1; // paint on top.
//             dwFlags = CURVE_Fill;
//             dwFlags |= CURVE_Bold;
            dwFlags &= ~CURVE_Ellipses;
            iCurveId = 0;
         } else {
//             if (pOrbital->pVisConfig->bColorSet)
// //                dwColor = (FColor(0.5f*(FColor(pOrbital->pVisConfig->cIsoPlus) + FColor(pOrbital->pVisConfig->cIsoMinus))).uint32());
//                dwColor = ct::irgb(FColor(0.5f*(FColor(pOrbital->pVisConfig->cIsoPlus) + FColor(pOrbital->pVisConfig->cIsoMinus))).uint32());
//             else
//                dwColor = 0xffffff;
//                // ^-- not set before rendering is finished...
//             dwColor |= 0xff000000;
            dwColor = 0xff000000 | pOrbital->GetBaseColor();
            pen = QPen(QBrush(QColor(dwColor)), 2.01, Qt::SolidLine);
         }

         if (!MakeIboChangeCurve(CurveData, iDataSet, document))
            // failed to make the delta-frame data set.
            continue;

//          PrintChargeArray(std::cout, "set curve", CurveData);
         QString
            Title = IvFmt("Orb. %1", iDataSet);
         if (IsCurrent)
            Title = "(current)";
         ui->curveView->setCurve(iCurveId, CurveData, Title, dwFlags, zOrder, pen, QBrush(dwColor));
         if (IsCurrent != 1) {
            nCurves += 1;
         }
      }
   }
   ui->curveView->UnlockUpdates();
   ui->curveView->update();
   iUpdateLocked -= 1;
}

void FMainWindow::onToggleTrackedOrbitalClicked()
{
   if (iUpdateLocked != 0)
      return;
   if (!document->GetFrameData(0) || !document->GetCurrentFrameData())
      return;
   iUpdateLocked += 1;
   FDataSetList
      *pDataList = document->GetCurrentFrameData();
//    ui->pushButton_ToggleTrackedOrbital->setChecked(pActiveData->Active);
   // ui->pushButton_ToggleTrackedOrbital->isChecked();
   int
      iRow = ui->dial_OrbitalId->value();
   if (pDataList && size_t(iRow) < pDataList->size()) {
      FDataSetPtr
         pActiveData = (*pDataList)[iRow];
      document->ToggleDataRow(iRow);
      view3d->updateGL(); // <- to force rendering the orbital such that we have a color set in the curve view...
      ui->pushButton_ToggleTrackedOrbital->setChecked(pActiveData->Active);
   }
//    document->GetCurrentFrame()->document->GetFrame(0)
   iUpdateLocked -= 1;
   onActiveIboCurveChanged();
   onActiveDatasetChanged(); // to enable/disable the orbital color controls.
}

void FMainWindow::UpdateAtomAndBondScalesText()
{
   ui->groupBox_AtomAndBondScales->setTitle(QString("Atom && Bond Size (%1, %2)").arg(ui->dialAtomScale->value()).arg(ui->dialBondScale->value()));
}


class QMakeWfThread : public QRunnable
{
public:
   QMakeWfThread(FLogQt &Log_, FDocument *pDocument_, QAbstractButton *pAbortButton_)
      : m_Log(Log_), m_pDocument(pDocument_), m_pAbortButton(pAbortButton_) {}
   void run();
protected:
//    ct::FLog &m_Log;
   FLogQt &m_Log;
   FDocument *m_pDocument;
   QAbstractButton *m_pAbortButton;
};

QString HtmlHightlight(QString s) { return "<pre style=\"color:white;\">" + s + "</pre>";  }

void QMakeWfThread::run()
{
   m_Log.Write(q2s(HtmlHightlight("*** WAVE FUNCTION COMPUTATION STARTING")));
   m_pDocument->RebuildWf(m_Log);
   m_Log.Write(q2s(HtmlHightlight("*** WAVE FUNCTION COMPUTATION FINISHED.")));
   m_Log.endReport();
   //if (m_pAbortButton)
   //   m_pAbortButton->setEnabled(false);
   // ^- hmpf... it's in a different thread. Can't do that.
}


void FMainWindow::onComputeWaveFunctionTriggered()
{
   if (!document->HaveConsistentFrames()) {
      IvNotify(NOTIFY_Error, "The loaded frames do not have compatible geometries (same number and types of atoms, in same order). RebuildWf not supported.");
      return;
   }
   FComputeWfForm
      WfDialog(document);
   if (bool(WfDialog.exec()) && (WfDialog.GetRunScf() || WfDialog.GetRunIbba())) {
      if (!WfDialog.IsMemoryOkay()) {
         QMessageBox::StandardButton
            btn = QMessageBox::question(this, "Potential Memory Overload",
               "This computation may need more memory "
               "than installed on this computer. We strongly recommend to perform it "
               "with a standard quantum chemistry package instead of IboView itself, and to just load "
               "the exported result for analysis."
               "\n\n"
               "Do you really wish to run this calculation now?",
               QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
         if (btn != QMessageBox::Yes)
            return;
      }

      FDataSetList
         ObjectLock;
      document->AcquireFrameObjectLock(ObjectLock);
      // ^- what is this? it is a list of *all* data set objects within
      // frames... because some of them have GL bindings, and GL cannot be
      // accessed from anything but the main thread. This is not restricted to
      // rendering. Trying to delete the GL resource bindings from another
      // thread crashes the application. I am not kidding. So we here keep smart
      // pointer refs to them to defer their deletion until the lock object goes
      // out of scope (which it does in the main thread).

//       IvNotify(NOTIFY_Warning, IvFmt("Compute WF triggered! Basis = %1", document->GetWfOptions()->GetOrbBasis()));
//       document->RebuildWf();
      FShowTextForm
         *pShowTextForm = new FShowTextForm(QString(), "SCF Progress", document->GetLogFileName("scf"), "Abort", this);
      pShowTextForm->setAttribute(Qt::WA_DeleteOnClose, true);
//       ct::FLogStdStream xxLog(std::cout);
//       ct::FLog *pLog = &xxLog;
//       FLogQt *pLog = new FLogQt(pShowTextForm);
      FMemoryLogQt *pLog = new FMemoryLogQt(pShowTextForm);
      connect(pLog, SIGNAL(textEmitted(QString)), pShowTextForm, SLOT(AppendText(QString const&)), Qt::QueuedConnection);
      connect(pLog, SIGNAL(sectionEnded()), pShowTextForm, SLOT(clear()), Qt::QueuedConnection);
      connect(pLog, SIGNAL(reportEnded(QString)), pShowTextForm, SLOT(setText(QString)), Qt::QueuedConnection);
      connect(pShowTextForm, SIGNAL(redButtonPressed()), pLog, SLOT(emitAbortSignal()), Qt::QueuedConnection);
//       connect(document, SIGNAL(layoutChanged()), pShowTextForm, SLOT(close()), Qt::QueuedConnection);
      QThreadPool::globalInstance()->start(new QMakeWfThread(*pLog, document, pShowTextForm->GetRedButton()));
//       QMakeWfThread test(*pLog, document, 0);
//       test.run();
      // ^- should this take ownership?
      pShowTextForm->exec();
      // FIXME: this will soooo blow up if the form is closed before the thread is done...

      // invalidate cached rendering data
      onRebuildIsoSurfacesClicked();
   }
}




struct FAboutDialogImpl
{
};

FAboutDialog::FAboutDialog(QWidget *parent)
   : QDialog(parent),
     ui(new Ui::AboutDialog),
     p(new FAboutDialogImpl)
{
   ui->setupUi(this);
   ui->svgWidget->sizePolicy().setHeightForWidth(true);
   ui->svgWidget->load(QString(":/resources/iboview_logo.svg"));
   ui->textBrowser->setHtml(LoadTextFromFile(":/resources/credits.htm"));
   if (!IvRestoreWindowSize("AboutDialog/Size", this))
      IvGuessSubDialogSize(this);
}

FAboutDialog::~FAboutDialog()
{
   IvSaveWindowSize("AboutDialog/Size", this);
   delete ui;
   delete p;
}

// void FAboutDialog::closeEvent(QCloseEvent *)
// {
//    return QDialog::closeEvent(event);
// }


// #include <QScriptEngine>

class FMyApplication;
static FMyApplication *g_papp = 0;


class FMyApplication : public QApplication
{
public:
   FMyApplication(int &argc, char **argv) : QApplication(argc, argv)
   {
      // reset numeric locale to "C", because otherwise string/float conversions in import/export
      // may get messed up. (TODO: still needs testing to see if this properly does the trick,
      // By default QCoreApplication sets locale to (LC_ALL, ""), which means "take local variants".)
      std::setlocale(LC_NUMERIC, "C");

      g_papp = this;
//       MainWindow.show();
#if QT_VERSION >= 0x050000
      // must be done after QApplication object is constructed. Otherwise undeployable
      // on windows ('could not load or find the Qt platform plugin "windows"')
      if (s_UseStyleFiles) {
         QApplication::setDesktopSettingsAware(false);
         QApplication::setStyle(QStyleFactory::create("Fusion"));
      }
#endif
      if (s_UseStyleFiles) {
         // load and apply our super cool style sheet 8).
         // (I could not resist)
         QFile StyleFile(":/resources/style.qss");
         StyleFile.open(QFile::ReadOnly);
         QString Text(StyleFile.readAll());
         this->setStyleSheet(Text);
      }

      pMainWindow = new IApplication();

      // check if we have a startup script to execute. If yes, execut it.
      if (1) {
         QSettings
            settings;
         QString
            StartupScriptText,
            StartupScriptFile = settings.value("IboView/StartupScriptFile").toString();
         if (!StartupScriptFile.isEmpty()) {
            StartupScriptText = LoadTextFileViaQt(StartupScriptFile);

            if (!StartupScriptText.isEmpty()) {
               ExecScript(pMainWindow, pMainWindow->view3d, StartupScriptText, QString("Startup Script: %1").arg(StartupScriptFile));
            }
         }
      }
   }
//    IApplication MainWindow;
   IApplication
      *pMainWindow;

   ~FMyApplication() {
      g_papp = 0;
      delete pMainWindow; // not in this->children.
      pMainWindow = 0;
   }
};

static QStringList
   s_CommandLineArgs;

void FMainWindow::processInputs()
{
   // I guess at this point only file names are left?
   load_files_(s_CommandLineArgs);

   IvNotify(NOTIFY_FinishWork, "");
}



void IvWarn(std::string const &Text)
{
   IvNotify(NOTIFY_Warning, s2q(Text));
}
void IvWarn(QString const &Text)
{
   IvNotify(NOTIFY_Warning, Text);
}

void IvEmit(QString const &Text)
{
//    if (Text.startsWith("!") || Text.startsWith("*")) {
//       FMainWindow
//          *pWin = qobject_cast<FMainWindow*>(g_pMainWindow);
//       std::cout << "^^^ STR: " << pWin << std::endl;
//       if (pWin) {
//          pWin->ui->statusBar->SetStatus(STATUS_Working, Text);
//       }
//    }

   std::cout << Text.toStdString() << std::endl;
}

void IvEmit(QString const &Text);

void IvEmit(QString Text, FEmitArg const &a0) { return IvEmit(Text.arg(a0.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1) { return IvEmit( Text.arg(a0.m, a1.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2) { return IvEmit( Text.arg(a0.m, a1.m, a2.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3) { return IvEmit( Text.arg(a0.m, a1.m, a2.m, a3.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4) { return IvEmit( Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a5) { return IvEmit( Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m, a5.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a5, FEmitArg const &a6) { return IvEmit( Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m, a5.m, a6.m)); }
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a5, FEmitArg const &a6, FEmitArg const &a7) { return IvEmit( Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m, a5.m, a6.m, a7.m)); }

QString IvFmt(QString Text, FEmitArg const &a0) { return Text.arg(a0.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1) { return Text.arg(a0.m, a1.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2) { return Text.arg(a0.m, a1.m, a2.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3) { return Text.arg(a0.m, a1.m, a2.m, a3.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4) { return Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a5) { return Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m, a5.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a5, FEmitArg const &a6) { return Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m, a5.m, a6.m); }
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a5, FEmitArg const &a6, FEmitArg const &a7) { return Text.arg(a0.m, a1.m, a2.m, a3.m, a4.m, a5.m, a6.m, a7.m); }

namespace fmt {
/*
   #define MAKE_FMT_X_FN(FmtName, FType) \
      QString FmtName(QString const &fmt, FType arg) { \
         return QString::fromStdString(boost::str(boost::format(fmt.toStdString()) % arg)); \
      }
      // ^- that must be the most expensive string formatting routine ever.
   MAKE_FMT_X_FN(fmti, long)
   MAKE_FMT_X_FN(fmti, unsigned long)
   MAKE_FMT_X_FN(fmti, int)
   MAKE_FMT_X_FN(fmti, unsigned int)
   MAKE_FMT_X_FN(fmtf, float)
   MAKE_FMT_X_FN(fmtf, double)
   #undef MAKE_FMT_X_FN

   QString fmts(QString const &fmt, QString s) {
      return QString::fromStdString(boost::str(boost::format(fmt.toStdString()) % s.toStdString()));
   }
*/
   inline QString fmt_width(QString const &s, int width) {
      if (width == 0)
         return s;
      else
         return QString("%1").arg(s, width);
   }
   // these ones use explicit width/prec/base specifications (note: still
   // options missing... e.g., field alignment, leading zeros, 'always use sign',
   // etc..).
   QString fmti(int i, int width, int base)                { return fmt_width(QString::number(i, (base==0? 10:base)), width); }
   QString fmti(unsigned int i, int width, int base)       { return fmt_width(QString::number(i, (base==0? 10:base)), width); }
   QString fmti(long i, int width, int base)                 { return fmt_width(QString::number(i, (base==0? 10:base)), width); }
   QString fmti(unsigned long i, int width, int base)        { return fmt_width(QString::number(i, (base==0? 10:base)), width); }
#ifdef NEED_EXTRA_SIZE_T_TYPES
   QString fmti(ptrdiff_t i, int width, int base)                 { return fmt_width(QString::number(i, (base==0? 10:base)), width); }
   QString fmti(size_t i, int width, int base)        { return fmt_width(QString::number(i, (base==0? 10:base)), width); }
#endif

   QString fmtf(double f, int width, int prec, char format) { return fmt_width(QString::number(f, format, prec), width); }
   QString fmtf(float f, int width, int prec, char format)  { return fmt_width(QString::number(f, format, prec), width); }

   QString fmts(QString const &s, int width) {
      return fmt_width(s, width);
      // ^- note: QString uses negative field widths to indicate left alignment.
      //          For strings that is usually what we want.
   }
   QString fmts(char const *p, int width) { return fmts(QString(p), width); }
   QString fmtp(void *p) { return "0x"+fmti(reinterpret_cast<size_t>(p),8,16); }
   // ^- should really be uintptr_t, but that's not standard C++.
}


FStatusClass NotifyClassToStatusClass(FNotificationClass NotifyClass) {
   switch (NotifyClass) {
      case NOTIFY_Information: return STATUS_Idle;
      case NOTIFY_StartWork: return STATUS_Working;
      case NOTIFY_FinishWork: return STATUS_Idle;
      case NOTIFY_Warning: return STATUS_Warning;
      case NOTIFY_Error: return STATUS_Error;
      default: return STATUS_Unknown;
   }
}


// FIXME: make a log class etc. This is just a hack for now.
static void UpdateStatus(FNotificationClass NotifyClass, QString const &Text)
{
   FMainWindow
      *pWin = qobject_cast<FMainWindow*>(g_pMainWindow);
   if (pWin) {
      pWin->ui->statusBar->SetStatus(NotifyClassToStatusClass(NotifyClass), Text);
   }
}

//    if (document->columnCount() == 0)
//       ui->statusBar->SetStatus(STATUS_Unknown, "Ready");
//    else
//       ui->statusBar->SetStatus(STATUS_Idle, "Ready");


void FMainWindow::customEvent(QEvent *pEvent)
{
   if (pEvent->type() == g_QtNotifyEventType) {
      FNotifyEvent
         *pne = dynamic_cast<FNotifyEvent*>(pEvent);
      assert_rt(pne != 0);
//       processNotifyEvent(pne);
      if (pne != 0)
         notifyUser(pne->notifyClass(), pne->message(), pne->title());
      pne->accept();
   } else
      QMainWindow::customEvent(pEvent);
}


void FMainWindow::notifyUser(FNotificationClass NotifyClass, QString const &Text, QString const &ExplicitTitle)
{
   if (NotifyClass == NOTIFY_StartWork || NotifyClass == NOTIFY_Information || NotifyClass == NOTIFY_FinishWork) {
      if (NotifyClass == NOTIFY_FinishWork || NotifyClass == NOTIFY_Information) {
         if (Text.isEmpty())
            return UpdateStatus(NotifyClass, "Ready");
      }
      return UpdateStatus(NotifyClass, Text);
   }
   {
      QString
         Title,
         Message = Text;
      if (!ExplicitTitle.isEmpty()) {
         Title = ExplicitTitle;
      } else {
         Title = "IboView Message";
         int iContextStart = Text.indexOf(":\t");
         if (iContextStart != -1) {
            Title = Text.mid(0, iContextStart);
            Message = Text.mid(iContextStart+2);
         }

         if (NotifyClass == NOTIFY_Error) {
            Title = "Error in " + Title;
         }
         else if (NotifyClass == NOTIFY_Warning) {
            Title = "Warning from " + Title;
         }
      }

      UpdateStatus(NotifyClass, Text);
      if (NotifyClass == NOTIFY_Error) {
         QMessageBox::critical(this, Title, Message);
      }
   //    else if (NotifyClass == NOTIFY_Warning)
   //       QMessageBox::warning(this, Title, Message);
   //    else
   //       QMessageBox::information(this, Title, Message);
   // ^- don't make message boxes for those by default.
   }
}


void IvNotify(FNotificationClass NotifyClass, QString const &Text)
{
   bool
      IsError = NotifyClass == NOTIFY_Error;
   std::ostream
      *pout = IsError? &std::cerr : &std::cout;
//    if (!(NotifyClass == NOTIFY_FinishWork || NotifyClass == NOTIFY_Information))
   if (NotifyClass == NOTIFY_Error) {
      *pout << "\n!ERROR: ";
   } else if (NotifyClass == NOTIFY_Warning) {
      *pout << "!WARNING: ";
   }
   if (!Text.isEmpty() || (NotifyClass == NOTIFY_Information || NotifyClass == NOTIFY_Warning || NotifyClass == NOTIFY_Error))
      *pout << q2s(Text) << std::endl;

   FMainWindow
      *parent = g_papp? g_papp->pMainWindow : 0;
   if (parent) {
      bool
         // if set, we can't access GUI events from here. Need to defer notification
         // over event queue.
         DeferMessagePost = (QThread::currentThread() != parent->thread());
      if (DeferMessagePost) {
         QCoreApplication::postEvent(parent, new FNotifyEvent(NotifyClass, Text));
      } else {
         parent->notifyUser(NotifyClass, Text);
      }
   }
}

int
   g_QtNotifyEventType = 0; // goes to event::type (comes from QEvent::registerEventType, and will make something between QEvent::User and QEvent::MaxUser)

FNotifyEvent::FNotifyEvent(FNotificationClass NotifyClass, QString const &Message, QString const &Title)
   : QEvent(static_cast<QEvent::Type>(g_QtNotifyEventType)), m_NotifyClass(NotifyClass), m_Message(Message), m_Title(Title)
{
}

FNotifyEvent::~FNotifyEvent()
{
}


// void IvNotify(FNotificationClass NotifyClass, QString const &Text)
// {
//    bool
//       IsError = NotifyClass == NOTIFY_Error;
//    std::ostream
//       *pout = IsError? &std::cerr : &std::cout;
//
//    if (NotifyClass == NOTIFY_StartWork || NotifyClass == NOTIFY_Information || NotifyClass == NOTIFY_FinishWork) {
//       if (!(NotifyClass == NOTIFY_FinishWork || NotifyClass == NOTIFY_Information))
//          *pout << q2s(Text) << std::endl;
//       else {
//          if (Text.isEmpty())
//             return UpdateStatus(NotifyClass, "Ready");
//       }
//       return UpdateStatus(NotifyClass, Text);
//    }
//
//    QString
//       Title = "IboView Message",
//       Message = Text;
//    int iContextStart = Text.indexOf(":\t");
//    if (iContextStart != -1) {
//       Title = Text.mid(0, iContextStart);
//       Message = Text.mid(iContextStart+2);
//    }
//
//    if (NotifyClass == NOTIFY_Error) {
//       *pout << "\n!ERROR: ";
//       Title = "Error in " + Title;
//    }
//    else if (NotifyClass == NOTIFY_Warning) {
//       *pout << "!WARNING: ";
//        Title = "Warning from " + Title;
//    }
//    *pout << q2s(Text) << std::endl;
// //    if (LongText == "")
// //       *pout << q2s(Text) << std::endl;
// //    else
// //       *pout << q2s(Text) << "\n" << q2s(LongText) << std::endl;
//    UpdateStatus(NotifyClass, Text);
//
//    QWidget
//       *parent = g_papp? g_papp->pMainWindow : 0;
//
//    if (NotifyClass == NOTIFY_Error) {
//       //       https://stackoverflow.com/questions/35018713/detect-that-im-running-in-qt-gui-event-thread
//       if (bool(parent)) {
//          if (QThread::currentThread() == parent->thread())
//             QMessageBox::critical(parent, Title, Message);
//          else
//             QCoreApplication::postEvent(parent, new FNotifyEvent(NOTIFY_Error, Title, Message));
//       }
//    }
// //    else if (NotifyClass == NOTIFY_Warning)
// //       QMessageBox::warning(parent, Title, Message);
// //    else
// //       QMessageBox::information(parent, Title, Message);
// // ^- don't make message boxes for those by default.
// }


// void IvNotify(FNotificationClass NotifyClass, QString const &Text, QString const &LongText)
// {
//    QString
//       Title = Text,
//       Message = LongText;
//    if (LongText == "") {
//       Title = "IboView";
//       Message = Text;
//    }
//
//    bool
//       IsError = NotifyClass == NOTIFY_Error;
//    std::ostream
//       *pout = IsError? &std::cerr : &std::cout;
//    if (NotifyClass == NOTIFY_Error)
//       *pout << "\n!ERROR: ";
//    else if (NotifyClass == NOTIFY_Warning)
//       *pout << "!WARNING: ";
//    if (LongText == "")
//       *pout << q2s(Text) << std::endl;
//    else
//       *pout << q2s(Text) << "\n" << q2s(LongText) << std::endl;
//
//    QWidget
//       *parent = g_papp? g_papp->pMainWindow : 0;
//
//    if (NotifyClass == NOTIFY_Error)
//       QMessageBox::critical(parent, Title, Message);
// //    else if (NotifyClass == NOTIFY_Warning)
// //       QMessageBox::warning(parent, Title, Message);
// //    else
// //       QMessageBox::information(parent, Title, Message);
// // ^- don't make message boxes for those by default.
// }

// QString GetUserName() {
//    // recycled from http://qt-project.org/forums/viewthread/11951.
//    QString
//       name = qgetenv("USER"); // get the user name in Linux
//    if (name.isEmpty())
//       name = qgetenv("USERNAME"); // get the name in Windows
//    return name;
// }


struct FCmdLineArg: public option::Arg
{
   static void printError(const char* msg1, const option::Option& opt, const char* msg2)
   {
      std::fprintf(stderr, "ERROR: %s", msg1);
      std::fwrite(opt.name, opt.namelen, 1, stderr);
      std::fprintf(stderr, "%s", msg2);
   }

   static option::ArgStatus Required(const option::Option& option, bool msg)
   {
      if (option.arg != 0)
         return option::ARG_OK;

      if (msg) printError("Option '", option, "' requires an argument\n");
      return option::ARG_ILLEGAL;
   }

   static option::ArgStatus Numeric(const option::Option& option, bool msg)
   {
      char* endptr = 0;
      if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
      if (endptr != option.arg && *endptr == 0)
         return option::ARG_OK;

      if (msg) printError("Option '", option, "' requires a numeric argument\n");
      return option::ARG_ILLEGAL;
   }
};


enum FOptionIndex { OPTION_Unknown, OPTION_Help, OPTION_ShowVirtuals, OPTION_SkipVirtuals, OPTION_NumThreads, OPTION_ViewAlpha, OPTION_AllowBackfaces, OPTION_UseStyleFile };
enum FOptionType { OPTTYPE_Other = 0, OPTTYPE_Enable = 1, OPTTYPE_Disable = 2};
static const option::Descriptor s_OptionDesc[] =
{
   {OPTION_Unknown, 0, "", "", FCmdLineArg::None, "USAGE: example [options]\n\n"
                                          "Options:" },
   {OPTION_Help, 0,"", "help", FCmdLineArg::None, "  --help  \tPrint usage and exit." },
//    {OPTION_ShowVirtuals, 0,"", "virt-orb", FCmdLineArg::None, "  --virt-orb  \tIf set, also load virtual orbitals from files, not only occupied orbitals." },
   {OPTION_SkipVirtuals, 0,"", "skip-virt", FCmdLineArg::None, "  --skip-virt  \tIf set, virtual orbitals will not be loaded from files." },
   {OPTION_NumThreads, 0, "t", "num-threads", FCmdLineArg::Numeric, "  -t <n>, --num-threads=<n>  \tSet the number of threads to use for computations (requires OpenMP)." },
   {OPTION_UseStyleFile, OPTTYPE_Disable, "", "disable-style", FCmdLineArg::None, "  --disable-style  \tIf set, disable custom style sheet for widgets. Some colors may look off." },
   {OPTION_ViewAlpha, OPTTYPE_Disable, "", "disable-alpha", FCmdLineArg::None, "  --disable-alpha  \tIf set, disable destination alpha in 3d window (workaround for MacOS window compositing)." },
   {OPTION_AllowBackfaces, OPTTYPE_Disable, "", "disable-backfaces", FCmdLineArg::None, "  --disable-backfaces  \tIf set, disable option to turn on orbital back faces (workaround for MacOS OpenGL driver bugs)." },
   {0,0,0,0,0,0}
};



int wrapped_main(int &argc, char *argv[])
{
   // make a "notification event" type id.
   assert(g_QtNotifyEventType == 0);
   g_QtNotifyEventType = QEvent::registerEventType();


   // argc/argv without program name.
   int argc1 = (argc > 0)? (argc - 1) : 0;
   char **argv1 = (argc > 0)? (argv + 1) : 0;

   // parse command line
   option::Stats
      stats(s_OptionDesc, argc1 , argv1);
   std::vector<option::Option>
      options(stats.options_max);
   std::vector<option::Option>
      buffer(stats.buffer_max);
   option::Parser
      parse(s_OptionDesc, argc1 , argv1, &options[0], &buffer[0]);

   if (parse.error())
      return 1;

   if (options[OPTION_Help]) {
      option::printUsage(std::cout, s_OptionDesc);
      return 0;
   }

   if (options[OPTION_ViewAlpha].last()->type() == OPTTYPE_Disable)
      g_WorkAroundAlphaCompositing = true;
   if (options[OPTION_AllowBackfaces].last()->type() == OPTTYPE_Disable)
      g_WorkAroundGlFrontFacingBug = true;

   for (int i = 0; i < parse.nonOptionsCount(); ++ i)
      s_CommandLineArgs.append(QString(parse.nonOption(i)));

   if (options[OPTION_NumThreads]) {
      g_nMaxOmpThreads = QString(options[OPTION_NumThreads].arg).toInt();
   } else {
      // take everything there is.
      g_nMaxOmpThreads = omp_get_num_procs();
      // ^- note that this ignores OMP_NUM_THREADS (unlike omp_get_max_threads)
   }
   omp_set_num_threads(g_nMaxOmpThreads);

//    g_ShowVirtualOrbitals = options[OPTION_ShowVirtuals];
   g_ShowVirtualOrbitals = !options[OPTION_SkipVirtuals];

//#if QT_VERSION >= 0x050000
//   QApplication::setDesktopSettingsAware(false);
//   QApplication::setStyle(QStyleFactory::create("Fusion"));
//#endif

   if (0) {
      // set a standard scheme in order to check if the style sheet is complete.
      QApplication::setDesktopSettingsAware(false);
      QApplication::setStyle(QStyleFactory::create("plastique"));
      //QApplication::setStyle(QStyleFactory::create("Fusion"));
      QStringList keys = QStyleFactory::keys();
      //QMessageBox::information(0, "known styles", keys.join("\n"));
      //qDebug("known styles:\n" + keys.join("\n"));
      for (int i = 0; i < keys.size(); ++i)
          std::cout << keys.at(i).toLocal8Bit().constData() << std::endl;
   }
   if (!options[OPTION_UseStyleFile] || options[OPTION_UseStyleFile].last()->type() != OPTTYPE_Disable) {
      s_UseStyleFiles = true;
   }

   ct::g_BasisSetLibrary.ImportMolproLib(":/bases/minao.libmol");
   ct::g_BasisSetLibrary.ImportMolproLib(":/bases/def2-nzvpp-orb.libmol");
   ct::g_BasisSetLibrary.ImportMolproLib(":/bases/def2-nzvpp-jfit.libmol");
   ct::g_BasisSetLibrary.ImportMolproLib(":/bases/def2-nzvpp-jkfit.libmol");
   ct::g_BasisSetLibrary.ImportMolproLib(":/bases/ecp-mdf.libmol");

   // these ones are for instanciating QSettings objects for persistent storage of application settings
   // (in particular, window size, splitter size, and position).
   QCoreApplication::setOrganizationName("Knizia-Research-Group");
   QCoreApplication::setOrganizationDomain("iboview.org");
   QCoreApplication::setApplicationName("IboView");
   QCoreApplication::setApplicationVersion("v20211019-RevA");
#ifdef Q_WS_X11
   // meant to allow message boxes (in particular, as created by IvNotify(NOTIFY_Error, ...) from other threads.
   // ...but it just crashes somewhere else.
   QCoreApplication::setAttribute(Qt::AA_X11InitThreads, true);
#endif
   FMyApplication
      app(argc, argv); // <- that's the main QApplication object.
   app.setWheelScrollLines(1); // mainly for IBO tracking. others work more or less fine with default value of 3.

   // process command line arguments, but not before the event loop
   // is started (otherwise we cannot quit the application via script!)
   QTimer::singleShot(0, app.pMainWindow, SLOT(processInputs()));
//    QTimer::singleShot(1000, app.pMainWindow, SLOT(processInputs()));
   return app.exec();
}



void _RegisterMetaTypes();

int main(int argc, char *argv[])
{
   _RegisterMetaTypes();
   if (0) {
      try {
         return wrapped_main(argc, argv);
      } catch (std::exception &e) {
//          std::cerr << boost::format("\n# Application crashed. Received Unhandled Exception:\n%s") % e.what() << std::endl;
         IvNotify(NOTIFY_Error, "Application crashed:\t" + QString("Received Unhandled Exception:\n%1").arg(QString(e.what())));
      }
   } else {
      return wrapped_main(argc, argv);
   }
   return -1;
}
