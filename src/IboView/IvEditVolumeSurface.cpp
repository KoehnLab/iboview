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

#include <QProgressDialog>
#include <QMessageBox>
#include <QtAlgorithms>
#include <QSettings>
#include <QVariant>
#include <fstream>
#include <QColorDialog>
#include <QRegExp>
#include "CxColor.h"

#include "Iv.h"
#include "IvEditVolumeSurface.h"
#include "ui_EditVolumeSurfaceForm.h"
#include "IvDocument.h"
#include "IvSettings.h"
#include "IvVolumeDataSet.h"


enum FVolumeIsoSurfaceType {
   SURFACETYPE_ElectronDensity,
   SURFACETYPE_SpinDensity,
   SURFACETYPE_Count
};


enum FSurfaceNameType {
   NAMETYPE_AsKey, // use as identifier key---no special characters
   NAMETYPE_AsDesc // use as description for user.
};


struct FEditVolumeSurfaceFormImpl
{
   FVolumePropertyPtr
      m_pVolumeProperty;
   FSupportDataPtr
      m_pSupportData;
   double
      m_VolumePropertyExactValue;
   uint32_t
      m_dwNextColor;
   bool
      m_SelectionChanging;
   FEditVolumeSurfaceFormImpl() {
//       m_dwNextColor = 0x70ff007f;
      m_dwNextColor = 0x7fff5e5b;
      m_SelectionChanging = false;
   }
};



static QString GetSurfaceTypeName(int SurfaceType, FSurfaceNameType NameType) {
   QString s;
   switch (SurfaceType) {
      case SURFACETYPE_ElectronDensity: s = QString("Electron Density"); break;
      case SURFACETYPE_SpinDensity: s = QString("Spin Density"); break;
      default:
         IvWarn(IvFmt("in GetSurfaceTypeName(): surface type index '%1' not recognized", SurfaceType));
         s = QString("Unknown Surface Type");
   }
   if (NameType == NAMETYPE_AsKey)
      // delete whitespace if used as keys.
      s.replace(QChar(' '), QString());
   return s;
}


// static int GetSurfaceTypeIndex(QString Key) {
//    QString
//       // normalize the input key
//       s = Key.toLower().trimmed().replace(QChar(' '),QString());
//    // for though all known surface types and check if their normalized entry agrees with the current one.
//    for (int i = 0; i < SURFACETYPE_Count; ++ i)
//       if (GetSurfaceTypeName(i, NAMETYPE_AsKey).toLower() == s)
//          return i;
//    return -1;
// }


static QString GetIsoThresholdsStorageKeyName(int SurfaceType) {
   return IvFmt("EditVolumeSurfaceWindow/IsoThresholds%1", GetSurfaceTypeName(SurfaceType, NAMETYPE_AsKey));
}


FEditVolumeSurfaceForm::FEditVolumeSurfaceForm(FDocument *document, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::EditVolumeSurfaceForm),
     p(new FEditVolumeSurfaceFormImpl),
     m_pDocument(document),
     m_SurfaceType(-1)
{
   ui->setupUi(this);

   restoreSavedIsoThresholds();

   connect(ui->toolButton_AddIsoThreshold, SIGNAL(clicked()), this, SLOT(addIsoThreshold()));
   connect(ui->toolButton_DeleteIsoThreshold, SIGNAL(clicked()), this, SLOT(deleteIsoThreshold()));
   connect(ui->toolButton_SetIsoColor, SIGNAL(clicked()), this, SLOT(changeSurfaceColor()));

   connect(ui->doubleSpinBox_IsoThreshold, SIGNAL(valueChanged(double)), this, SLOT(changeIsoThreshold(double)));
   connect(ui->tabWidget_SurfaceType, SIGNAL(currentChanged(int)), this, SLOT(changeSurfaceType(int)));

   changeSurfaceType(0); // select density by default

   ui->buttonBox->setFocus();

   if (!IvRestoreWindowSize("EditVolumeSurfaceWindow/Size", this))
      IvGuessSubDialogSize(this);
   IvRestoreSplitterState("EditVolumeSurfaceWindow/Splitter", ui->splitter);


   // hmpf. This stuff doesn't work. Can use the variants themselves with custom types,
   // but if I put them into QList<QVariant>s, they cannot be serialized. No idea. Giving
   // up now and storing them as strings...
//    {
//       IvEmit("trying QVariant::store...");
//       FIsoSliceDecl
//          d;
//       QVariant stored;
//       stored.setValue(d);
//
//       IvEmit("trying QVariant::load...");
//       FIsoSliceDecl
//          e = qvariant_cast<FIsoSliceDecl>(stored);
//       IvEmit("...out.");
//    }
//    static bool MetaTypeRegistered = false;
//    if (!MetaTypeRegistered) {
//       MetaTypeRegistered = true;
//       qRegisterMetaType<FIsoSliceDecl>();
//    }
   connect(ui->listWidget_IsoThresholds->selectionModel(), SIGNAL(currentChanged(const QModelIndex &, const QModelIndex &)), this, SLOT(currentThresholdChanged(const QModelIndex &, const QModelIndex &)));
}


FEditVolumeSurfaceForm::~FEditVolumeSurfaceForm()
{
   IvSaveSplitterState("EditVolumeSurfaceWindow/Splitter", ui->splitter);
   IvSaveWindowSize("EditVolumeSurfaceWindow/Size", this);
   delete p;
   delete ui;
}


void FEditVolumeSurfaceForm::restoreSavedIsoThresholds()
{
   QSettings
      settings;
   m_IsoThresholdsForSurfaceTypes.clear();
   // make one table of iso-thresholds for each iso-surface type
   for (int i = 0; i < SURFACETYPE_Count; ++ i) {
      FIsoThresholdList
         L; // <- empty list

      // if present, restore iso-threshold list last computed for this surface type
      QVariant
         vaTresholdList = settings.value(GetIsoThresholdsStorageKeyName(i));
      if (vaTresholdList.isValid()) {
//          QList<QVariant>
//             vaList = vaTresholdList.toList();
         QStringList
            vaList = vaTresholdList.toStringList();
         for (int k = 0; k < vaList.size(); ++ k) {
//             bool ok;
//             double d = vaList[k].toDouble(&ok);
//             if (vaList[k].canConvert<FIsoSliceDecl>()) {
//                FIsoSliceDecl
//                   decl = vaList[k].value<FIsoSliceDecl>();
//                L.append(decl);
//             } else
//                IvWarn(IvFmt("Failed to load iso-threshold list '%1' from configuration: '%2' is not a iso-slice declaration.", GetIsoThresholdsStorageKeyName(i), vaList[k].toString()));
            bool ok;
            FIsoSliceDecl decl = FIsoSliceDecl::fromString(vaList[k], &ok);
//             double d = vaList[k].toDouble(&ok);
            if (ok)
               L.append(decl);
            else
               IvWarn(IvFmt("Failed to load iso-threshold list '%1' from configuration: '%2' is not a iso-slice declaration.", GetIsoThresholdsStorageKeyName(i), vaList[k]));
         }
      }

      if (L.empty()) {
         // no previous data stored. init to something semi-sensible.
         if (i == SURFACETYPE_ElectronDensity) {
            L.append(FIsoSliceDecl::fromString(QString("(0.05,7fff5e5b)")));
            L.append(FIsoSliceDecl::fromString(QString("(0.005,50ffc65a)")));
         }
         if (i == SURFACETYPE_SpinDensity) {
            L.append(FIsoSliceDecl::fromString(QString("(-0.05,7f5b5eff)")));
            L.append(FIsoSliceDecl::fromString(QString("(0.05,7fff5e5b)")));
            L.append(FIsoSliceDecl::fromString(QString("(-0.005,505ac6ff)")));
            L.append(FIsoSliceDecl::fromString(QString("(0.005,50ffc65a)")));
         }
      }

      m_IsoThresholdsForSurfaceTypes.append(L);
   }
}

FIsoSliceDecl FIsoSliceDecl::fromString(QString s, bool *ok)
{
   FIsoSliceDecl
      r;
   IvEmit("read: '%1'", s);
   QRegExp
//       re("^\\((@[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?@),([0123456789abcdef]{0,8})\\)$");
//       re("\\((@[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?@)\,([0123456789abcdef]*)\\)");
      re("\\(([0123456789+-.eE]*),([0123456789abcdef]*)\\)");
      // ^- I hate regex.
   if (re.exactMatch(s)) {
      QStringList
         L = re.capturedTexts();
      if (L.size() == 3) {
         if (ok)
            *ok = true;

         r.fIsoValue = L[1].toDouble(ok);
         if (ok && *ok) {
            r.dwColor = L[2].toUInt(ok, 16);
         }
      }
   } else {
      if (ok)
         *ok = false;
   }
   return r;
}

QString FIsoSliceDecl::toString() const
{
   return QString("(%1,%2)").arg(fIsoValue).arg(QString::number(dwColor,16));
}


// Initial color: r = 127  g = 0  b = 255  a = 112
// New color: 1358933596 (50ffae5c)
// Initial color: r = 92  g = 174  b = 255  a = 80
// New color: 2197777497 (82ff6c59)


void FEditVolumeSurfaceForm::saveIsoThresholds(int SurfaceType)
{
//    IvEmit("trying saveIsoThresholds...");
   QSettings
      settings;
   assert(m_IsoThresholdsForSurfaceTypes.size() == SURFACETYPE_Count);
   for (int i = 0; i < SURFACETYPE_Count; ++ i) {
      if (SurfaceType != -1 && i != SurfaceType)
         continue;

      FIsoThresholdList const
         &L = m_IsoThresholdsForSurfaceTypes[i]; // <- current list

//       // convert to list of QVariants.
//       QList<QVariant>
//          vaList;
//       for (int k = 0; k < L.size(); ++ k) {
//          vaList.append(QVariant(L[k]));
//          vaList.append(QVariant::fromValue(L[k]));
//          QVariant va;
//          va.setValue(L[k]);
//          vaList.append(va);
//       }
//       settings.setValue(GetIsoThresholdsStorageKeyName(i), QVariant(vaList));
//       IvEmit("setting LIST of variants now.");
//       settings.setValue(GetIsoThresholdsStorageKeyName(i), QVariant::fromValue(vaList));
//       settings.setValue(GetIsoThresholdsStorageKeyName(i), vaList);

      // convert to list of QStrings
//       QList<QString>
      QStringList
         vaList;
      for (int k = 0; k < L.size(); ++ k) {
         vaList.append(L[k].toString());
      }
      settings.setValue(GetIsoThresholdsStorageKeyName(i), vaList);
   }
//    IvEmit("...out.");
}


void FEditVolumeSurfaceForm::updateThresholdListInUi()
{
   assert(m_SurfaceType < m_IsoThresholdsForSurfaceTypes[m_SurfaceType].size());

   // get current selection
   QModelIndexList
      selected = ui->listWidget_IsoThresholds->selectionModel()->selectedIndexes();

   FIsoThresholdList
      &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];
   ui->listWidget_IsoThresholds->clear();
   for (int i = 0; i < L.size(); ++ i) {
//       ui->listWidget_IsoThresholds->addItem(QString("rho(x) == %1").arg(L[i],16,'f',8));
      QString
         ItemDesc = QString("rho(x): %1  col: 0x%2").arg(L[i].fIsoValue).arg(QString::number(ct::irgb(L[i].dwColor),16).rightJustified(8,'0'));
      QListWidgetItem
         *pNewItem = new QListWidgetItem(ItemDesc, 0);
      pNewItem->setBackground(QBrush(QColor(ct::irgb(L[i].dwColor))));
      if ((L[i].dwColor & 0xff) < 0x7f)
         pNewItem->setForeground(Qt::white);
      else
         pNewItem->setForeground(Qt::black);
      ui->listWidget_IsoThresholds->addItem(pNewItem);
   }

   // re-apply the previous selection as far as possible
   QItemSelection
      newSelection;
   foreach(QModelIndex const &index, selected) {
      if (index.row() < L.size())
         newSelection.select(index, index);
   }
   ui->listWidget_IsoThresholds->selectionModel()->select(newSelection, QItemSelectionModel::Select);
}


void FEditVolumeSurfaceForm::addIsoThreshold()
{
   double
      d = ui->doubleSpinBox_IsoThreshold->value();
   if (d == 0) {
      QMessageBox::warning(this, "Add Iso Threshold", "Zero is not a valid iso-surface threshold. Command ignored.");
      return;
   }
   if (d < 0 && m_SurfaceType == SURFACETYPE_ElectronDensity) {
      // only allow non-negative tresholds for total densities
      QMessageBox::warning(this, "Add Iso Threshold", "Only positive iso-thresholds allowed for total densities. Command ignored.");
      return;
   }
//    if (m_SurfaceType == SURFACETYPE_SpinDensity || m_SurfaceType == SURFACETYPE_ElectronDensity) {
//    if (m_SurfaceType == SURFACETYPE_ElectronDensity) {
//       // only allow non-negative tresholds for total densities
//       if (d < 0) d = -d;
//    }
   assert(m_SurfaceType < m_IsoThresholdsForSurfaceTypes.size());
   FIsoThresholdList
      &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];

//    // value already in there?
//    if (qFind(L, d) == L.end()) {
//       // nope. Add it.
//       L.append(d);
//       qSort(L);
//       updateThresholdListInUi();
//    }
   L.append(FIsoSliceDecl(d, p->m_dwNextColor));
   qSort(L);
   updateThresholdListInUi();
}

bool FIsoSliceDecl::operator < (FIsoSliceDecl const &other) const
{
   // large absolute values (inner surfaces) first.
   if (std::abs(fIsoValue) > std::abs(other.fIsoValue)) return true;
   if (std::abs(other.fIsoValue) > std::abs(fIsoValue)) return false;
   // then minus first, then plus
   if (fIsoValue < other.fIsoValue) return true;
   if (other.fIsoValue < fIsoValue) return false;
   // then sort by color.
   if (dwColor < other.dwColor) return true;
   if (other.dwColor < dwColor) return false;
   return false;
}

void FEditVolumeSurfaceForm::deleteIsoThreshold()
{
   assert(m_SurfaceType < m_IsoThresholdsForSurfaceTypes[m_SurfaceType].size());

   QModelIndexList
      selected = ui->listWidget_IsoThresholds->selectionModel()->selectedIndexes();

   // get indices of selected thresholds
   QList<int>
      indices;
   foreach(QModelIndex const &index, selected) {
      indices << index.row();
   }
   // sort in reverse order (for deletion -- this keeps the remaining indices constant).
   qSort(indices.begin(), indices.end(), qGreater<int>());

   // deselect all selected rows
   foreach(QModelIndex const &index, selected) {
      ui->listWidget_IsoThresholds->selectionModel()->select(index, QItemSelectionModel::Deselect);
   }


//    int iThresholdToDelete = ui->listWidget_IsoThresholds->currentRow();
   // remove backend version of the selected thresholds.
   foreach(int iThresholdToDelete, indices) {
      if (iThresholdToDelete == -1) {
         QMessageBox::warning(this, "Delete Iso Threshold", "No iso-threshold selected in list.");
         return;
      }
      FIsoThresholdList
         &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];
      if (iThresholdToDelete > L.size()) {
         QMessageBox::warning(this, "Delete Iso Threshold", "internal inconsistency in UI: selected row does not exist in backend! Canceling deletion.");
         return;
      }

      // delete the item from the list.
      L.takeAt(iThresholdToDelete);
   }
   updateThresholdListInUi();
}


static FVolumePropertyType SurfaceTypeToVolumePropertyType(int SurfaceType) {
   switch (SurfaceType) {
      case SURFACETYPE_ElectronDensity: {
         return PROPERTY_Density;
      }
      case SURFACETYPE_SpinDensity: {
         return PROPERTY_SpinDensity;
      }
      default: {
         IvNotify(NOTIFY_Error, IvFmt("failed to convert surface type '%1' to FVolumePropertyType.", SurfaceType));
         return FVolumePropertyType(0);
      }
   }
}


void FEditVolumeSurfaceForm::changeIsoThreshold(double newValue)
{
   if (m_SurfaceType == SURFACETYPE_ElectronDensity || m_SurfaceType == SURFACETYPE_SpinDensity) {
      // compute number of electrons excluded with the current threshold.
//       double fElecExcluded = 0.1235467890123456 + newValue;

      if (p->m_pVolumeProperty.get() != 0 && p->m_pSupportData.get() != 0) {

         double fElecExcluded = p->m_pVolumeProperty->GetIntegratedValueBelowThreshold(newValue, p->m_pSupportData);

         ui->label_ThresholdInfo->setText(QString("(^~ %1 electrons in density below this threshold)").arg(fElecExcluded,0,'f',6));
      } else {
         ui->label_ThresholdInfo->setText("(no data)");
      }
   }

   if (!p->m_SelectionChanging) {
      // apply new threshold to selected rows -- if any.
      // (don't do it if currently switching between rows... because at this point the old rows
      // are still selected, and without the guard we'd assign the new row's data to the old rows)
      QModelIndexList
         selected = ui->listWidget_IsoThresholds->selectionModel()->selectedIndexes();
      if (!selected.empty()) {
         assert(m_SurfaceType < m_IsoThresholdsForSurfaceTypes.size());
         FIsoThresholdList
            &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];

         foreach(QModelIndex const &index, selected) {
            int iThreshold = index.row();
            if (iThreshold < L.size())
               L[iThreshold].fIsoValue = newValue;
         }
         updateThresholdListInUi();
      }
   }
}


void FEditVolumeSurfaceForm::changeSurfaceColor()
{
   uint32_t
      dwInitialColor = p->m_dwNextColor;
   QColor
      InitialColor = QColor::fromRgba(ct::irgb(dwInitialColor));
   IvEmit("Initial color: r = %1  g = %2  b = %3  a = %4", InitialColor.red(), InitialColor.green(), InitialColor.blue(), InitialColor.alpha());

   QColorDialog
      ColorDialog(InitialColor, 0);
   ColorDialog.setOption(QColorDialog::ShowAlphaChannel, true);
   ColorDialog.setCurrentColor(InitialColor);
   // ^- without this initial alpha gets discarded (likely because we set color with alpha
   //    before turning on the alpha channel).
   int
      Result = ColorDialog.exec();
   if (Result == QDialog::Accepted) {
      // update the color value.
      uint32_t
         NewColor = ct::irgb(ColorDialog.selectedColor().rgba());

      IvEmit("New color: %1", NewColor);
      p->m_dwNextColor = NewColor;

      {
         // apply new color to selected rows -- if any.
         QModelIndexList
            selected = ui->listWidget_IsoThresholds->selectionModel()->selectedIndexes();
         if (!selected.empty()) {
            assert(m_SurfaceType < m_IsoThresholdsForSurfaceTypes.size());
            FIsoThresholdList
               &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];

            foreach(QModelIndex const &index, selected) {
               int iThreshold = index.row();
               if (iThreshold < L.size())
                  L[iThreshold].dwColor = NewColor;
            }
            updateThresholdListInUi();
         }
      }
   }
}


void FEditVolumeSurfaceForm::currentThresholdChanged(const QModelIndex &current, const QModelIndex &previous)
{
   IR_SUPPRESS_UNUSED_WARNING(previous);
   p->m_SelectionChanging = true;

   assert(m_SurfaceType < m_IsoThresholdsForSurfaceTypes.size());
   FIsoThresholdList
      &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];
   if (current.row() < L.size() && current.row() >= 0) {
      FIsoSliceDecl
         &Decl = L[current.row()];
      ui->doubleSpinBox_IsoThreshold->setValue(Decl.fIsoValue);
      p->m_dwNextColor = Decl.dwColor;
   }

   p->m_SelectionChanging = false; // should i mutex here? not sure if this can multithread...
}


void FEditVolumeSurfaceForm::resetVolumePropertyData()
{
   FVolumePropertyInfo
      info(SurfaceTypeToVolumePropertyType(m_SurfaceType));
   FFrame
      *pFrame = m_pDocument->GetCurrentFrame();
   if (pFrame == 0) {
      ui->label_ThresholdInfo->setText("(no current frame)");
      p->m_pSupportData = 0;
      p->m_pVolumeProperty = 0;
      p->m_VolumePropertyExactValue = 0;
      return;
   } else {
      p->m_pVolumeProperty = new FVolumeProperty(pFrame, info, 0);
      p->m_VolumePropertyExactValue = p->m_pVolumeProperty->GetExactValue();
   }

   {
      size_t
         nMemNeeded = p->m_pVolumeProperty->GetMaxStackDataSize();
      ct::FMemoryStack2
         Mem((size_t(20) << size_t(20)) + nMemNeeded);
         // ^- 20 mb + whatever the volume property would like to have.
      p->m_pSupportData = p->m_pVolumeProperty->GetSupportData(Mem);
   }

   // just to update the text...
   changeIsoThreshold(ui->doubleSpinBox_IsoThreshold->value());
}


void FEditVolumeSurfaceForm::updateWfInformation(int SurfaceType)
{
   {
      QProgressDialog
         progress(IvFmt("Computing %1...", GetSurfaceTypeName(SurfaceType, NAMETYPE_AsDesc)), QString(), 0, 10, this);
      progress.setMinimumDuration(0);
      progress.setWindowModality(Qt::WindowModal);
      progress.setValue(0);

      // do stuff here.
      resetVolumePropertyData();
      double fGridValue = p->m_pVolumeProperty->GetIntegratedValue(p->m_pSupportData);
      ui->label_ElecFromGridIntegration->setText(QString::number(fGridValue,'f',6));
      ui->label_NumElec->setText(QString::number(p->m_VolumePropertyExactValue));

      if (m_pDocument->GetCurrentFrame()) {
         FWfType
            WfType = m_pDocument->GetCurrentFrame()->GetWfType();
         QString s;
         switch (WfType) {
            case wfi::WFTYPE_Rhf: {
               s = "RHF/RKS";
               break;
            }
            case wfi::WFTYPE_Uhf: {
               s = "UHF/UKS";
               break;
            }
            case wfi::WFTYPE_Mcscf: {
               s = "MCSCF";
               break;
            }
//             case WFTYPE_Other: {
//                s = "Other";
//                break;
//             }
            default: {
               s = "(unknown)";
            }
         }
         ui->label_WfType->setText(s);
      } else {
         ui->label_WfType->setText("(no frame)");
      }

      progress.setValue(10);
   }
}


void FEditVolumeSurfaceForm::changeSurfaceType(int newIndex)
{
   if (newIndex < 0)
      newIndex = 0;
   if (newIndex > int(SURFACETYPE_Count))
      newIndex = int(SURFACETYPE_Count)-1;
   if (m_SurfaceType == newIndex)
      return;
   m_SurfaceType = newIndex;
   updateWfInformation(newIndex);
   updateThresholdListInUi();
}


void FEditVolumeSurfaceForm::accept()
{
   // collect list of iso-surface thresholds.
   // ...
   if (m_IsoThresholdsForSurfaceTypes[m_SurfaceType].size() == 0) {
      QMessageBox::warning(this, "Add Iso Surfaces", "No iso-thresholds in list. Cannot add surface.\n(use 'Add' button to add surfaces to the list, or use 'Cancel' to close the dialog without adding surfaces.).");
      return; // <- this cancels closing the dialog.
   }

   // add a corresponding iso-surface data set to the document (for all frames)
   // ...

//    if (m_pDocument->GetCurrentFrame()) {
//       // note: this is NOT the right way to do it! should make a new script-accessible
//       // function in FDocument which adds the surfaces to all frames.
//       FFrame
//          *pFrame = m_pDocument->GetCurrentFrame();
//       FVolumePropertyInfo
//          info(SurfaceTypeToVolumePropertyType(m_SurfaceType));
//       pFrame->m_Data.push_back(new FVolumeProperty(pFrame, info, 0));
//    }
   {
      FVolumePropertyInfo
         info(SurfaceTypeToVolumePropertyType(m_SurfaceType));
      FIsoValueList
         IsoValues;
      FIsoThresholdList const
         &L = m_IsoThresholdsForSurfaceTypes[m_SurfaceType];
      for (int i = 0; i < L.size(); ++ i) {
//          IsoValues.push_back(FIsoThresholdEntry(L[i], (L[i] < 0)? 0x4000ff00 : 0x40ff00ff));
         IsoValues.push_back(FIsoThresholdEntry(L[i].fIsoValue, L[i].dwColor));
      }

      m_pDocument->AddVolumePropertyIsoSurface(info, IsoValues);
   }


   // remember current settings
   saveIsoThresholds(m_SurfaceType);

   // close the dialog.
   FBase::accept();
}
