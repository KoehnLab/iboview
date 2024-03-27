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

#ifndef IV_EDIT_FRAMES_FORM
#define IV_EDIT_FRAMES_FORM

#include <QDialog>
#include <QCloseEvent>
#include <QAbstractTableModel>
// #include "IvStatusBar.h"

namespace Ui{
    class EditFramesForm;
}

class FDocument;
class FFrameListModel;

typedef std::vector<int>
   FFrameIndexList;

class FFrameListModel : public QAbstractTableModel
{
   Q_OBJECT
public:
   typedef QAbstractTableModel
      FBase;
   int rowCount(const QModelIndex &parent = QModelIndex()) const;
   int columnCount(const QModelIndex &parent = QModelIndex()) const;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
   bool removeRows(int row, int count, const QModelIndex &parent = QModelIndex());

   QModelIndex makeRowIndex(int iRow);

   void beginReset1();
   void endReset1();
   void updateData1();
public:
   FFrameListModel(FDocument *pDocument_, QObject *Parent_);
   ~FFrameListModel();
public slots:
   // forwarded from FDocument
   void parentLayoutAboutToBeChanged();
   void parentLayoutChanged();
protected:
   FDocument
      *m_pDocument;
   void ConnectForwardSignals();
};


class FEditFramesForm : public QDialog
{
   Q_OBJECT

   Ui::EditFramesForm *ui;
   FDocument
      *m_pDocument;
public:
   explicit FEditFramesForm(FDocument *document, QWidget *parent=0);
   ~FEditFramesForm();

   void DoPostProcessingOnClose();
public slots:
   void moveFramesUp();
   void moveFramesDown();
   void moveFramesToTop();
   void moveFramesToBottom();
   void reverseFrameOrder();
   void deleteFrames();

   void alignFrames();
   void deleteCloseFrames();
   void exportFramesAsXyz();

   void UpdateExcludedAtoms();

   void changeWeightMode(int iNewMode);
protected:
   FFrameListModel
      *m_pFrameList;
   bool
      // set if we changed the order or orientation of frames such that it
      // may be required to re-link the orbital visual configurations
      m_OrbitalRelinkNeeded;

   void SetNewFrameOrder(FFrameIndexList const &NewIndices, bool UpdateSelection=true);
   FFrameIndexList GetBaseFrameOrder();
   FFrameIndexList GetSelectedFrameIndices();

//    void closeEvent(QCloseEvent *event); // override
};


#endif // IV_EDIT_FRAMES_FORM
