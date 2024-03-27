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

#include <algorithm> // for std::sort
#include <QItemSelectionModel>
#include "IvFindOrbitalsForm.h"
#include "IvSettings.h"
#include "IvOrbital.h"
#include "ui_FindOrbitalsForm.h"

// represent a selection of orbitals found via "Find Orbitals"
// on selected atoms.
FFoundOrbitalModel::FFoundOrbitalModel(FFoundOrbitalList const &L, QString CharacterName_, FDocument *pDocument, QObject *parent_)
   : QAbstractTableModel(parent_), m_pDocument(pDocument), m_FoundOrbitals(L), m_CharacterName(CharacterName_)
{
   connect(m_pDocument, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex &)), this, SLOT(documentDataChanged(const QModelIndex &, const QModelIndex &)));
}

int FFoundOrbitalModel::rowCount(const QModelIndex &/*parent*/) const
{
   return (int)m_FoundOrbitals.size();
}

int FFoundOrbitalModel::columnCount(const QModelIndex &/*parent*/) const
{
   return 3;
}

QVariant FFoundOrbitalModel::headerData(int section, Qt::Orientation orientation, int role) const
{
//    char const *pCaptions[] = {"Sel. Chg.", "Desc.", "Centers/Charges"};
   QString Captions[] = {m_CharacterName, QString("Desc."), QString("Centers/Charges")};
   if ( orientation == Qt::Horizontal ) {
      if ( role == Qt::DisplayRole ) {
         if ((size_t)section < (sizeof(Captions)/sizeof(Captions[0])))
//             return QString(pCaptions[section]);
            return Captions[section];
      }
      if ( role == Qt::FontRole ) {
            QFont CaptionFont;
            CaptionFont.setBold(true);
            return CaptionFont;
      }
   }

   return QVariant();
}

QVariant FFoundOrbitalModel::data(const QModelIndex &index, int role) const
{
   int
      iRow = index.row(),
      iCol = index.column();
   if (iRow < 0 && (size_t)iRow < m_FoundOrbitals.size())
      return QVariant();
   FFoundOrbital const
      *pFoundOrbital = &m_FoundOrbitals[iRow];
   FOrbital const
      *pOrbital = pFoundOrbital->pOrbital;

   if ( role == Qt::DisplayRole ) {
      if (iCol == 0) {
         QString
            sChgOnSelected(QString("%1").arg(pFoundOrbital->fCharacter,7,'f',3));
         if (pFoundOrbital->fCharacter > 1.5)
            sChgOnSelected.append("*");
         else
            sChgOnSelected.append(" ");
         return QVariant(sChgOnSelected);
      } else if (iCol == 1) {
         return QVariant(pOrbital->GetDesc());
      }else if (iCol == 2) {
         return QVariant(pOrbital->MakeFullDesc(0.02, FOrbital::ORBDESC_ChargesOnly));
      }
      return QVariant();
   }

//    if ( role == Qt::FontRole && iCol == 2 ) {
//       QFont MonoFont("");
//       MonoFont.setStyleHint(QFont::TypeWriter);
//       // ^- there is also QFont::Monospace, but that does not seem to work.
//       return MonoFont;
//    }

   if ( pOrbital->Active && iCol == 1 ) {
      uint32_t dwBaseColor = pOrbital->GetBaseColor();
      if (dwBaseColor != 0) {
         if ( role == Qt::BackgroundRole ) return QBrush(QColor(dwBaseColor));
         if ( role == Qt::ForegroundRole ) return QBrush(Qt::black);
      } else {
         if ( role == Qt::BackgroundRole ) return QBrush(QColor(64,96,64));
         if ( role == Qt::ForegroundRole ) return QBrush(Qt::white);
      }
   }
   return QVariant();
}

void FFoundOrbitalModel::documentDataChanged(const QModelIndex &/*topLeft*/, const QModelIndex &/*bottomRight*/)
{
//    IvEmit("forwarding FFoundOrbitalModel::documentDataChanged.");
   // well.. we COULD sort out the right range... or just emit "update all"...
   emit dataChanged(createIndex(0,0), createIndex(3, m_FoundOrbitals.size()));
}






FFindOrbitalsForm::FFindOrbitalsForm(FFoundOrbitalModel *pModel_, QString Title_, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::FindOrbitalsForm),
     m_pModel(pModel_)
{
   ui->setupUi(this);
   ui->orbitalTable->setWordWrap(false);
   ui->orbitalTable->setSelectionBehavior(QAbstractItemView::SelectRows);
#if QT_VERSION >= 0x050000
   ui->orbitalTable->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
   ui->orbitalTable->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
#else
   // qt5 doc's say there is a compatibiltiy layer for this.. but it doesn't seem to exist
   // in my windows qt 5.3.2 version.
   ui->orbitalTable->horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
   ui->orbitalTable->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);
#endif
   ui->orbitalTable->setModel(m_pModel);
   setWindowTitle(Title_);
   // todo: connect actions etc (toggle data row, sort).

   connect(ui->orbitalTable, SIGNAL(doubleClicked(const QModelIndex&)), this, SLOT(toggleRow(const QModelIndex&)));
   connect(ui->pushButton_ToggleRows, SIGNAL(clicked()), this, SLOT(toggleSelectedRows()));

   if (!IvRestoreWindowSize("FindOrbitalsWindow/Size", this))
      IvGuessSubDialogSize(this);
}

FFindOrbitalsForm::~FFindOrbitalsForm()
{
   IvSaveWindowSize("FindOrbitalsWindow/Size", this);
   delete ui;
}

void FFindOrbitalsForm::toggleRow(const QModelIndex &index)
{
   int
      iRow = index.row();
   if ((size_t)iRow > m_pModel->m_FoundOrbitals.size())
      return;
   FFoundOrbital const
      *pFoundOrbital = &m_pModel->m_FoundOrbitals[(size_t)iRow];
//    IvEmit("trigger row %1", pFoundOrbital->iDataRow);

   m_pModel->m_pDocument->ToggleDataRow(pFoundOrbital->iDataRow);
}

void FFindOrbitalsForm::toggleSelectedRows()
{
   QItemSelectionModel
      *pSelectionModel = ui->orbitalTable->selectionModel();
   QModelIndexList
      SelectedRows = pSelectionModel->selectedRows();
   QModelIndexList::const_iterator
      itSel;
   for (itSel = SelectedRows.begin(); itSel != SelectedRows.end(); ++ itSel) {
      toggleRow(*itSel);
   }
}


struct FFoundOrbitalSortPred
{
   bool
      m_Descending;
   explicit FFoundOrbitalSortPred(bool Descending_) : m_Descending(Descending_) {}

   bool operator() (FFoundOrbital const &a, FFoundOrbital const &b) {
      if (m_Descending)
         return a.fCharacter > b.fCharacter;
      else
         return a.fCharacter < b.fCharacter;
   }
};


void SortFoundOrbitalList(FFoundOrbitalList &L, bool Descending)
{
   std::stable_sort(L.begin(), L.end(), FFoundOrbitalSortPred(Descending));
}
