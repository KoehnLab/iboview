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

#ifndef IV_TABLES_H
#define IV_TABLES_H

#include <QAbstractTableModel>
#include <QString>
#include <vector>
#include <QDialog>
#include "CxTypes.h"
#include "CxPodArray.h"

class FDocument;
class FFrame;
namespace Ui{
    class TablesForm;
}

class FDocument;



// support for tables/curves of user-definable measurements on the
// loaded frames (e.g., IRC, bond lengths, charges, IBO displacements, etc.)

typedef ct::TArray<double>
   FMeasureSet; // one entry for each frame.
typedef ct::TArray<FFrame*>
   FFrameMeasurePtrList;

// should I make them QObjects instead? Not sure.
// also not sure if "make everything given a frame list"
// is really that good of an idea. It might be expensive
// and/or re-calculate stuff often. Might need cache
// or invalidate settings or similar.
struct FMeasure : public ct::FIntrusivePtrDest
{
   enum FFmtFlags {
      FORMAT_ValueOnly = 0x01
   };

   explicit FMeasure(FDocument *pDocument_, QString const &UnitName_, int nDigits_, double fAuToOutFactor_);
   virtual void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const = 0;
   virtual QString Name() const = 0;
   virtual ~FMeasure();
   virtual QString UnitName() const; // default implementation returns m_UnitName.
   virtual QString FmtValue(double f, uint Flags=0) const;

   // a convenience function for computing values for single frames for info purposes.
   // Flags: FORMAT_*
   QString MeasureFrame(FFrame *pFrame, uint Flags=0) const;
protected:
   FDocument *m_pDocument;
   QString m_UnitName;
   int m_nDigits; // maybe put in a format.h-format instead? would also be better for curves and changing units etc.
   double m_AuToOut; // factor for converting from atomic units to target unit.
};
typedef ct::TIntrusivePtr<FMeasure>
   FMeasurePtr;

struct FMeasureBondLength : public FMeasure
{
   explicit FMeasureBondLength(int iAt_, int jAt_, FDocument *pDocument_);

   void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const; // override.
   QString Name() const; // override.
protected:
   int
      m_iAt, m_jAt;
};

struct FMeasureBondAngle : public FMeasure
{
   explicit FMeasureBondAngle(int iAt_, int jAt_, int kAt_, FDocument *pDocument_);

   void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const; // override.
   QString Name() const; // override.
protected:
   int
      m_iAt, m_jAt, m_kAt;
};

struct FMeasurePlaneAngle : public FMeasure
{
   explicit FMeasurePlaneAngle(int iAt_, int jAt_, int kAt_, int lAt_, FDocument *pDocument_);

   void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const; // override.
   QString Name() const; // override.
protected:
   int
      m_iAt, m_jAt, m_kAt, m_lAt;
};


struct FMeasureFrameId : public FMeasure
{
   explicit FMeasureFrameId(FDocument *pDocument_);

   void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const; // override.
   QString Name() const; // override.
};

struct FMeasureFrameGradient : public FMeasure
{
   explicit FMeasureFrameGradient(FDocument *pDocument_);

   void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const; // override.
   QString Name() const; // override.
};

struct FMeasureFrameEnergy : public FMeasure
{
   explicit FMeasureFrameEnergy(FDocument *pDocument_);

   void MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const; // override.
   QString Name() const; // override.
};


typedef std::vector<FMeasurePtr>
   FMeasureList;
// columns == measurements, rows == frames.
class FDocumentMeasures : public QAbstractTableModel
{
   Q_OBJECT
public:
   typedef QAbstractTableModel
      FBase;
   int rowCount(const QModelIndex &parent = QModelIndex()) const;
   int columnCount(const QModelIndex &parent = QModelIndex()) const;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
   bool removeColumns(int row, int count, const QModelIndex &parent = QModelIndex());

   void AddMeasure(FMeasurePtr pMeasure);

   void ConnectForwardSignals();
public:
   FDocumentMeasures(FDocument *pDocument_, QObject *Parent_);
   ~FDocumentMeasures();
public slots:
   // forwarded from FDocument
   void parentLayoutAboutToBeChanged();
   void parentLayoutChanged();
protected:
   FDocument
      *m_pDocument;
   FMeasureList
      m_List;
   FMeasure *GetMeasure(QModelIndex const &index);
   FMeasure *GetMeasure(int iRow);
   FMeasure const *GetMeasure(QModelIndex const &index) const {return const_cast<FDocumentMeasures*>(this)->GetMeasure(index);}
   FMeasure const *GetMeasure(int iRow) const {return const_cast<FDocumentMeasures*>(this)->GetMeasure(iRow);}
};


struct FTablesForm : public QDialog
{
   Q_OBJECT

   Ui::TablesForm
      *ui;
   FDocument
      *m_pDocument;
public:
   typedef QDialog FBase;
   explicit FTablesForm(FDocument *pDocument_, QWidget *pParent_ = 0);
};


#endif // IV_TABLES_H
