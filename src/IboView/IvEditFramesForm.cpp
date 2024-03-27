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

#include <algorithm>
#include <QTextStream>
#include <QItemSelectionModel>
#include <QProgressDialog>
#include <fstream>

#include "Iv.h"
#include "IvEditFramesForm.h"
#include "ui_EditFramesForm.h"
#include "IvDocument.h"
#include "IvSettings.h"



FFrameListModel::FFrameListModel(FDocument *pDocument_, QObject *Parent_)
   : QAbstractTableModel(Parent_), m_pDocument(pDocument_)
{
   m_pDocument->UpdateArcLengths();
   ConnectForwardSignals();
}

FFrameListModel::~FFrameListModel()
{
}


int FFrameListModel::rowCount(const QModelIndex &parent) const
{
   return m_pDocument->GetNumFrames();
   IR_SUPPRESS_UNUSED_WARNING(parent);
}

int FFrameListModel::columnCount(const QModelIndex &parent) const
{
   return 2; // ID, measure, name
   IR_SUPPRESS_UNUSED_WARNING(parent);
}

QVariant FFrameListModel::headerData(int section, Qt::Orientation orientation, int role) const
{
   if (role == Qt::DisplayRole && orientation == Qt::Vertical) {
      FFrame
         *pFrame = m_pDocument->GetFrame(section, false);
      if (pFrame == 0)
         return QVariant("#UNK");
      return IvFmt("%1", section); // frame number.
   }
   if (role == Qt::DisplayRole && orientation == Qt::Horizontal) {
      if (section == 0)
         return QVariant("s (Arc)");
      if (section == 1)
         return QVariant("Frame");
      return QVariant("#UNK");
   }
   if ( role == Qt::FontRole ) {
         QFont CaptionFont;
         CaptionFont.setBold(true);
         return CaptionFont;
   }
   return QVariant();
}

QVariant FFrameListModel::data(const QModelIndex &index, int role) const
{
   if (role == Qt::DisplayRole) {
//       FMeasure const
//          *pMeasure = GetMeasure(index);
//       FFrame
//          *pFrame = m_pDocument->GetFrame(index.row(),false);
//       if (pFrame == 0 || pMeasure == 0)
//          return QVariant("#UNK");
//       return QVariant(pMeasure->MeasureFrame(pFrame,FMeasure::FORMAT_ValueOnly));
      FFrame
         *pFrame = m_pDocument->GetFrame(index.row(), false);
      if (index.column() == 1 && pFrame != 0)
         return QVariant(pFrame->GetFullInputFileName());
      if (index.column() == 0 && pFrame != 0)
         return QVariant(QString("%1").arg(pFrame->GetArcLength().Value(), 12, 'f', 5));
   }
   if (role == Qt::TextAlignmentRole) {
      if (index.column() == 0)
         return QVariant(Qt::AlignRight);
      else
         return QVariant(Qt::AlignLeft);
   }
   return QVariant();
}

bool FFrameListModel::removeRows(int row, int count, const QModelIndex &parent)
{
   if (count < 0 || row + count >= rowCount()) {
      IV_NOTIFY(NOTIFY_Warning, "Attempted to remove non-existent frame.");
      return false;
   }
//    beginRemoveRows(parent, row, row + count - 1);
//    FMeasureList::iterator
//       itFirst = m_List.begin() + row,
//       itLast = itFirst + count;
//    m_List.erase(itFirst, itLast);
//    endRemoveRows();
//    return true;
   return false;
   IR_SUPPRESS_UNUSED_WARNING(parent);
}

QModelIndex FFrameListModel::makeRowIndex(int iRow)
{
   return createIndex(iRow, 0, (void*)0);
}


void FFrameListModel::ConnectForwardSignals()
{
//    // forward those things... to deal with newly inserted frames etc.
//    connect(m_pDocument, SIGNAL(layoutAboutToBeChanged()), this, SLOT(parentLayoutAboutToBeChanged()));
//    connect(m_pDocument, SIGNAL(layoutChanged()), this, SLOT(parentLayoutChanged()));
   // (handled via model reset now.)
}

void FFrameListModel::parentLayoutAboutToBeChanged()
{
   emit layoutAboutToBeChanged();
}

void FFrameListModel::parentLayoutChanged()
{
   emit layoutChanged();
}


void FFrameListModel::beginReset1()
{
   emit beginResetModel();
}

void FFrameListModel::endReset1()
{
   emit endResetModel();
}

void FFrameListModel::updateData1()
{
   emit dataChanged(QModelIndex(), QModelIndex()); // <- claim that all data changed.
//    emit beginResetModel();
//    emit endResetModel();
}



static QString GetWeightModeForIndex(int iMode) {
   QString
      AlignMode = "??";
   switch(iMode) {
      case 0: AlignMode = "mass"; break;
      case 1: AlignMode = "iso_mass"; break;
      case 2: AlignMode = "charge"; break;
      case 3: AlignMode = "coords"; break; // every atom gets the same weight (except excluded ones)
   }
   return AlignMode;
}


static int GetWeightModeIndexForDesc(QString Desc) {
   if (Desc == "mass") return 0;
   if (Desc == "iso_mass") return 1;
   if (Desc == "charge") return 2;
   if (Desc == "coords") return 3;
   IV_NOTIFY(NOTIFY_Warning, "Encountered unexpected atomic weight mode in FEditFramesForm.");
   return 0;
}


FEditFramesForm::FEditFramesForm(FDocument *document, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::EditFramesForm),
     m_pDocument(document),
     m_pFrameList(new FFrameListModel(document, this)),
     m_OrbitalRelinkNeeded(false)
{
   ui->setupUi(this);
//    ui->label_WorkSpace->setVisible(false);
//    ui->spinBox_WorkSpacePerThread->setVisible(false);
//    ui->spinBox_WorkSpacePerThread->setEnabled(false);

//    ui->tabWidget->layout()->setContentsMargins(0, 0, 0, 0); // left, top, right, bottom
//    connect(ui->checkBox_RunIbba, SIGNAL(toggled(bool)), this, SLOT(ToggleIbbaPage(bool)));

   ui->tableView_FrameList->setModel(m_pFrameList);

   UpdateExcludedAtoms();
   ui->comboBox_AlignWeight->setCurrentIndex(GetWeightModeIndexForDesc(m_pDocument->GetAtomWeightMode()));

   connect(ui->toolButton_MoveUp, SIGNAL(clicked()), this, SLOT(moveFramesUp()));
   connect(ui->toolButton_MoveDown, SIGNAL(clicked()), this, SLOT(moveFramesDown()));
   connect(ui->toolButton_MoveToTop, SIGNAL(clicked()), this, SLOT(moveFramesToTop()));
   connect(ui->toolButton_MoveToBottom, SIGNAL(clicked()), this, SLOT(moveFramesToBottom()));
   connect(ui->toolButton_ReverseFrames, SIGNAL(clicked()), this, SLOT(reverseFrameOrder()));
   connect(ui->toolButton_DeleteSelected, SIGNAL(clicked()), this, SLOT(deleteFrames()));
   connect(ui->toolButton_AlignFrames, SIGNAL(clicked()), this, SLOT(alignFrames()));
   connect(ui->comboBox_AlignWeight, SIGNAL(currentIndexChanged(int)), this, SLOT(changeWeightMode(int)));
   connect(ui->toolButton_DeleteCloseFrames, SIGNAL(clicked()), this, SLOT(deleteCloseFrames()));
   connect(ui->toolButton_ExportXyz, SIGNAL(clicked()), this, SLOT(exportFramesAsXyz()));

   ui->buttonBox->setFocus();
//    setAttribute(Qt::WA_DeleteOnClose); // should I do this?

   ui->tableView_FrameList->resizeColumnsToContents();

   if (!IvRestoreWindowSize("EditFramesWindow/Size", this))
      IvGuessSubDialogSize(this);
}

FEditFramesForm::~FEditFramesForm()
{
   IvSaveWindowSize("EditFramesWindow/Size", this);
   delete ui;
}


void FEditFramesForm::UpdateExcludedAtoms()
{
   // todo: go though atom flags of document and check which ones have the flag set.
   int
      nAt = m_pDocument->nAtomsWithFlags();
   QString
      s;
   QTextStream
      str(&s);
   str << "Atoms Excluded: ";
   int nAtomsExcluded = 0;
   for (int iAt = 0; iAt < nAt; ++ iAt) {
      if (m_pDocument->IsAtomExcludedFromAlignment(iAt)) {
         if (nAtomsExcluded != 0)
            str << " | ";
         str << m_pDocument->AtomLabel(iAt);
         nAtomsExcluded += 1;
      }
   }
//    if (nAtomsExcluded <= 0) {
//       ui->label_ExcludedAtoms->setText("(no atoms excluded.)");
//    } else {
//       ui->label_ExcludedAtoms->setText(s);
//    }
   if (nAtomsExcluded <= 0)
      str << "None";
   ui->label_ExcludedAtoms->setText(s);
}


static void InvertOrder(FFrameIndexList &iAll_, FFrameIndexList &iSelected_)
{
   FFrameIndexList
      iAll,
      iSelected;
   iAll.reserve(iAll_.size());
   size_t
      n = iAll_.size(),
      s = iSelected_.size();
   for (size_t i = n-1; i < n; -- i)
      iAll.push_back(iAll_[i]);
   iSelected.reserve(iSelected_.size());
   for (size_t i = s-1; i < s; -- i)
      iSelected.push_back(n - iSelected_[i] - 1);
   iAll_.swap(iAll);
   iSelected_.swap(iSelected);
}

static void MoveSelectedUp(FFrameIndexList &iAll, FFrameIndexList const &iSelected)
{
   if (iAll.empty() || iSelected.empty())
      return; // no frames there or nothing to move.
//    assert(iAll[0] == 0 || iAll.empty());
   assert(iAll[0] == 0 || iAll[0] == int(iAll.size()-1));

   // note: this logic is for handling cases like that:
   //   #Frame   Selected?
   //     0        yes
   //     1        yes
   //     2        no
   //     3        no
   //     4        yes
   //     5        yes
   // In this case everything which can be moved up, should be. The new order would
   // become:  0 1 2 4 5 3. But Frames 0 and 1 cannot be moved up, since they already
   // are at the top positions.
   int i = 0;
//    while (iSelected[i] == iAll[i] && i < int(iSelected.size()))
   while (iSelected[i] == i && i < int(iSelected.size()))
      ++ i;
   for (; i < int(iSelected.size()); ++ i)
   {
      int
         iOld = iSelected[i],
         iNew = iOld - 1;
      assert(iNew > 0);
      std::swap(iAll[iOld], iAll[iNew]);
   }
}

static bool isSelected(FFrameIndexList::value_type iVal, FFrameIndexList const &iSelected)
{
   FFrameIndexList::const_iterator
      itSel = std::lower_bound(iSelected.begin(), iSelected.end(), iVal);
   return itSel != iSelected.end() && *itSel == iVal;
}

static void MoveSelectedToTop(FFrameIndexList &iAll_, FFrameIndexList const &iSelected)
{
   // algorithm requires selected items to be sorted and mutually distinct.
// #ifdef _DEBUG
   for (size_t i = 1; i < iSelected.size(); ++ i)
      assert(iSelected[i-1] < iSelected[i]);
// #endif
   // start out with the selected items.
   FFrameIndexList
      iAll;
   iAll.reserve(iAll_.size());
   for (size_t i = 0; i < iSelected.size(); ++ i)
      iAll.push_back(iAll_[size_t(iSelected[i])]);

   // add all other items in original order, unless they were already included in the
   // sorted subset.
   for (size_t i = 0; i < iAll_.size(); ++ i) {
      if (!isSelected(int(i), iSelected))
         iAll.push_back(iAll_[i]);
   }
   iAll_.swap(iAll);
}


void FEditFramesForm::moveFramesUp()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iSelected = GetSelectedFrameIndices();
   MoveSelectedUp(iAll, iSelected);
   SetNewFrameOrder(iAll);
}

void FEditFramesForm::moveFramesDown()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iSelected = GetSelectedFrameIndices();
   InvertOrder(iAll, iSelected);
   MoveSelectedUp(iAll, iSelected);
   InvertOrder(iAll, iSelected);
   SetNewFrameOrder(iAll);
}

void FEditFramesForm::moveFramesToTop()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iSelected = GetSelectedFrameIndices();
   MoveSelectedToTop(iAll, iSelected);
   SetNewFrameOrder(iAll);
}

void FEditFramesForm::moveFramesToBottom()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iSelected = GetSelectedFrameIndices();
   InvertOrder(iAll, iSelected);
   MoveSelectedToTop(iAll, iSelected);
   InvertOrder(iAll, iSelected);
   SetNewFrameOrder(iAll);
}

void FEditFramesForm::reverseFrameOrder()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iSelected = GetSelectedFrameIndices();
   if (iAll.empty())
      return;
   if (iSelected.empty() || iSelected.size() == 1u)
      // nothing selected? invert entire list.
      iSelected = iAll;
   assert(iAll[0] == 0 && iAll.back() == int(iAll.size() - 1));

   // go through the first half of the selected elements, and exchange
   // them with the second half. In case of odd number, center element stays
   // where it is.
   for (size_t i = 0; i != iSelected.size()/2; ++ i) {
      size_t
         j = iSelected.size() - i - 1;
      assert(i < j);
      size_t
         iFrame = iSelected[i],
         jFrame = iSelected[j];
      assert(iFrame << iAll.size() && jFrame < iAll.size());
      std::swap(iAll[iFrame], iAll[jFrame]);
   }
//    for (size_t i = 0; i != iAll.size(); ++ i)
//       IvEmit("  NewFrame %1:   %2", i, iAll[i]  );
   SetNewFrameOrder(iAll);
}

void FEditFramesForm::deleteFrames()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iNew,
      iSelected = GetSelectedFrameIndices();
   if (iAll.empty() || iSelected.empty())
      return;

   // go through the first half of the selected elements, and exchange
   // them with the second half. In case of odd number, center element stays
   // where it is.
   for (size_t i = 0; i != iAll.size(); ++ i) {
      if (!isSelected(i, iSelected))
         iNew.push_back(iAll[i]);
   }
   SetNewFrameOrder(iNew, false);
   // ^- false: do not update selection. In practice
   // this means that the cursor stays in the row it was in before, which
   // now points to the next frame. (good for deleting stuff sequentially)
}

void FEditFramesForm::deleteCloseFrames()
{
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iNew,
      iSelected = GetSelectedFrameIndices();
   if (iSelected.empty())
      iSelected = iAll;
   if (iSelected.empty())
      return;

   double
      fTargetArcLength = std::abs(ui->doubleSpinBox_FrameDistance->value());

   for (size_t i = 0; i != iAll.size(); ++ i) {
      // keep all frames which are not selected, and both the first and the last selected frame.
      if (iAll[i] == iSelected.front() || iAll[i] == iSelected.back() || !isSelected(i, iSelected)) {
         iNew.push_back(iAll[i]);
         continue;
      }
      assert(iNew.size() != 0);
      double
         fLastArcLength = m_pDocument->GetFrame(iNew.back())->GetArcLength().Value(),
         fThisArcLength = m_pDocument->GetFrame(iAll[i])->GetArcLength().Value();
      if (std::abs(fLastArcLength - fThisArcLength) < fTargetArcLength)
         continue;
      else
         iNew.push_back(iAll[i]);
   }

   SetNewFrameOrder(iNew, true);
   // ^- true: update selection.
}

void FEditFramesForm::alignFrames()
{
   if (GetWeightModeForIndex(ui->comboBox_AlignWeight->currentIndex()) != m_pDocument->GetAtomWeightMode()) {
      IV_NOTIFY(NOTIFY_Warning, "Atomic weight unexpectedly inconsistent between document and edit frames form.");
      m_pDocument->SetAtomWeightMode(GetWeightModeForIndex(ui->comboBox_AlignWeight->currentIndex()));
   }

   m_pDocument->AlignFrames(QString());
   m_OrbitalRelinkNeeded = true;

   m_pFrameList->updateData1(); // <- to update the arc lengths
}

void FEditFramesForm::exportFramesAsXyz()
{
//    IvNotify(NOTIFY_Error, "REINSTANTIATE exportFramesAsXyz -- bad git commit.");
   FFrameIndexList
      iAll = GetBaseFrameOrder(),
      iSelected = GetSelectedFrameIndices();
   if (iSelected.empty())
      iSelected = iAll;
   if (iSelected.empty())
      return IvNotify(NOTIFY_Error, "No frames---nothing to export.");

   QString
      FileNameSuggest = ReplaceExt(m_pDocument->GetCurrentInputFileName(), ".xyz"),
      FileName;
   FileName = IvGetSaveFileName("History/SaveFiles", this, "Export Frames to .xyz",
         (FileNameSuggest), "xyz with gradients (*.xyz);;All Files (*.*)");

   if (FileName == "")
      return;
   std::ofstream
      out(q2s(FileName).c_str());
   if (!out.good())
      return IvNotify(NOTIFY_Error, IvFmt("Failed to open file '%1' for writing.", FileName));
   for (size_t ii = 0; ii != iSelected.size(); ++ii) {
      m_pDocument->GetFrame(iSelected[ii])->pGetAtoms()->PrintAsXyz(out);
      if (out.fail())
         return IvNotify(NOTIFY_Error, IvFmt("Failed write data to file '%1'.", FileName));
   }
}


void FEditFramesForm::changeWeightMode(int iNewMode)
{
   m_pDocument->SetAtomWeightMode(GetWeightModeForIndex(iNewMode));
   m_pFrameList->updateData1();
//    IvEmit("Atom weight mode changed to '%1'.", GetWeightModeForIndex(iNewMode));
}


void FEditFramesForm::SetNewFrameOrder(FFrameIndexList const &NewIndices, bool UpdateSelection)
{
   FFrameIndexList
      iSelectedOld = GetSelectedFrameIndices();

   bool
      FrameNumberChanges = (NewIndices.size() != size_t(m_pDocument->GetNumFrames()));
   if (FrameNumberChanges)
      m_pFrameList->beginReset1();

   m_pDocument->ReorderOrRestrictFrameSet(NewIndices);

   if (UpdateSelection) {
      // keep track of selection: Select new subset of frames which was selected before.
      QItemSelection
         newSelection;
      for (size_t i = 0; i < NewIndices.size(); ++ i) {
         if (isSelected(NewIndices[i], iSelectedOld))
            newSelection.push_back(QItemSelectionRange(m_pFrameList->makeRowIndex(int(i))));
      }
      QItemSelectionModel
         *selectionModel = ui->tableView_FrameList->selectionModel();
      selectionModel->select(newSelection, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
   }
   m_OrbitalRelinkNeeded = true;

   if (FrameNumberChanges)
      m_pFrameList->endReset1();
   else
      m_pFrameList->updateData1();
}

FFrameIndexList FEditFramesForm::GetBaseFrameOrder()
{
   FFrameIndexList r;
   r.reserve(m_pDocument->GetNumFrames());
   for (int iFrame = 0; iFrame < m_pDocument->GetNumFrames(); ++ iFrame)
      r.push_back(iFrame);
   return r;
}

FFrameIndexList FEditFramesForm::GetSelectedFrameIndices()
{
   QItemSelectionModel
      *selectionModel = ui->tableView_FrameList->selectionModel();
   QModelIndexList
      selectedRows = selectionModel->selectedRows();
   FFrameIndexList
      r;
   r.reserve(selectedRows.size());
   for (int i = 0; i < selectedRows.size(); ++ i)
      r.push_back(selectedRows[i].row());
   std::sort(r.begin(), r.end());
   return r;
}


// this doesn't work reliably...
// only called if the window is closed by pressing X, not via "ok"/"cancel",
// or other means.

// void FEditFramesForm::closeEvent(QCloseEvent *event)
// {
//    IvEmit("ENTERED CLOSE-EVENT.");
//    if (m_OrbitalRelinkNeeded) {
//       IvEmit("ENTERED RE-LINK ORBITALS.");
//       QProgressDialog
//          progress("Re-Linking corresponding orbitals...", QString(), 0, m_pDocument->GetNumFrames(), this);
//       progress.setMinimumDuration(0);
//       progress.setWindowModality(Qt::WindowModal);
//       progress.setValue(0);
//       m_pDocument->LinkOrbitals();
//       progress.setValue(m_pDocument->GetNumFrames());
//    }
//
//    QDialog::closeEvent(event);
// }

void FEditFramesForm::DoPostProcessingOnClose()
{
   IvEmit("ENTERED POST-PROCESS-EVENT.");
   if (m_OrbitalRelinkNeeded) {
      IvEmit("ENTERED RE-LINK ORBITALS.");
      QProgressDialog
         progress("Re-Linking corresponding orbitals...", QString(), 0, m_pDocument->GetNumFrames(), this);
      progress.setMinimumDuration(0);
      progress.setWindowModality(Qt::WindowModal);
      progress.setValue(0);
      m_pDocument->LinkOrbitals();
      progress.setValue(m_pDocument->GetNumFrames());
   }
}
