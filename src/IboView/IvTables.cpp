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

#include <QAbstractProxyModel>
#include <algorithm>
#include <math.h>
#include "CxVec3.h"
#include "Iv.h"
#include "IvDocument.h"
#include "IvTables.h"
#include "CxPhysicalUnits.h"
#include "ui_TablesForm.h"
using std::sqrt;

// // c/p'd from http://www.qtcentre.org/wiki/index.php?title=Transpose_proxy_model
// class FTransposeProxyModel : public QAbstractProxyModel
// {
// public:
//    FTransposeProxyModel(QObject *p = 0) : QAbstractProxyModel(p){}
//
//    QModelIndex mapFromSource( const QModelIndex & sourceIndex ) const{
//       return index(sourceIndex.column(), sourceIndex.row());
//    }
//    QModelIndex mapToSource( const QModelIndex & proxyIndex ) const{
//       return sourceModel()->index(proxyIndex.column(), proxyIndex.row());
//    }
//    QModelIndex index(int r, int c, const QModelIndex &index=QModelIndex()) const{
//       return createIndex(r,c);
//       IR_SUPPRESS_UNUSED_WARNING(index);
//    }
//    QModelIndex parent(const QModelIndex&) const {
//       return QModelIndex();
//    }
//    int rowCount(const QModelIndex &) const{
//       return sourceModel()->columnCount();
//    }
//    int columnCount(const QModelIndex &) const{
//       return sourceModel()->rowCount();
//    }
//    QVariant data(const QModelIndex &ind, int role) const {
//       return sourceModel()->data(mapToSource(ind), role);
//    }
// };

// ^- hmmm.. this doesn't work. headerData arguments are wrong (section always supplied as 0),
//    and orientation is not patched up. maybe I just switch everything manually...









FMeasure::FMeasure(FDocument *pDocument_, QString const &UnitName_, int nDigits_, double fAuToOutFactor_)
   : m_pDocument(pDocument_), m_UnitName(UnitName_), m_nDigits(nDigits_), m_AuToOut(fAuToOutFactor_)
{
}

FMeasure::~FMeasure()
{
}

QString FMeasure::UnitName() const
{
   return m_UnitName;
}

QString FMeasure::FmtValue(double f, uint Flags) const
{
   if (0 != (Flags & FORMAT_ValueOnly))
      return QString("%2").arg(f,0,'f',m_nDigits);
   else
      return QString("%1 = %2 %3").arg(Name()).arg(f,0,'f',m_nDigits).arg(UnitName());
}

QString FMeasure::MeasureFrame(FFrame *pFrame, uint Flags) const
{
   FFrameMeasurePtrList
      DummyList;
//    DummyList.push_back(FFramePtr(pFrame));
   DummyList.push_back(pFrame);
   FMeasureSet
      Value;
   MakeValues(Value, DummyList);
   if (Value.size() != 1) {
      IvNotify(NOTIFY_Error, IvFmt("Something went wrong while computing measure %1.", this->Name()));
      return "!MEASURE-ERROR";
   }
   return FmtValue(Value[0], Flags);
}



FMeasureBondLength::FMeasureBondLength(int iAt_, int jAt_, FDocument *pDocument_)
   : FMeasure(pDocument_,"pm",2,ct::ToAng*100.), m_iAt(iAt_), m_jAt(jAt_)
{}

QString FMeasureBondLength::Name() const
{
   return QString("r(%1--%2)").arg(m_pDocument->AtomLabel(m_iAt), m_pDocument->AtomLabel(m_jAt));
}

void FMeasureBondLength::MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const
{
   FrameValues.clear();
   FrameValues.reserve(Frames.size());
//    foreach(FFrameMeasurePtrList::value_type const *itFrame, Frames) {
   FFrameMeasurePtrList::const_iterator itFrame;
   _for_each(itFrame, Frames) {
      ct::FAtomSet *pAtoms = (*itFrame)->pGetAtoms();
      if (pAtoms && size_t(m_iAt) < pAtoms->size() && size_t(m_jAt) < pAtoms->size()) {
         ct::FVector3
            vI = (*pAtoms)[size_t(m_iAt)].vPos,
            vJ = (*pAtoms)[size_t(m_jAt)].vPos;
         FrameValues.push_back(m_AuToOut * Dist(vI, vJ));
      } else {
         FrameValues.push_back(0);
      }
   }
}

// return angle between Left--Center and Center--Right (in radians)
static double Angle(ct::FVector3 const &Left, ct::FVector3 const &Center, ct::FVector3 const &Right)
{
   ct::FVector3
      vA = Normalized(Left - Center),
      vB = Normalized(Right - Center);
   return std::acos(Dot(vA, vB));
}


FMeasureBondAngle::FMeasureBondAngle(int iAt_, int jAt_, int kAt_, FDocument *pDocument_)
   : FMeasure(pDocument_,QString::fromWCharArray(L"\u00b0"),2,180./M_PI), m_iAt(iAt_), m_jAt(jAt_), m_kAt(kAt_)
{}

QString FMeasureBondAngle::Name() const
{
   return QString(QString::fromWCharArray(L"\u03b1(%1 | %2 | %3)")).arg(m_pDocument->AtomLabel(m_iAt), m_pDocument->AtomLabel(m_jAt), m_pDocument->AtomLabel(m_kAt));
}

void FMeasureBondAngle::MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const
{
   FrameValues.clear();
   FrameValues.reserve(Frames.size());
//    foreach(FFrameMeasurePtrList::value_type const *itFrame, Frames) {
   FFrameMeasurePtrList::const_iterator itFrame;
   _for_each(itFrame, Frames) {
      ct::FAtomSet *pAtoms = (*itFrame)->pGetAtoms();
      if (pAtoms && size_t(m_iAt) < pAtoms->size() && size_t(m_jAt) < pAtoms->size()) {
         ct::FVector3
            vI = (*pAtoms)[size_t(m_iAt)].vPos,
            vJ = (*pAtoms)[size_t(m_jAt)].vPos,
            vK = (*pAtoms)[size_t(m_kAt)].vPos;
         FrameValues.push_back(m_AuToOut * Angle(vI, vJ, vK));
      } else {
         FrameValues.push_back(0);
      }
   }
}


// return angle between Left--Center and Center--Right (in radians)
static double Dihedral(ct::FVector3 const &vA, ct::FVector3 const &vB, ct::FVector3 const &vC, ct::FVector3 const &vD)
{
   ct::FVector3
      v0 = vB - vA,
      v1 = vC - vB,
      v2 = vD - vC,
      vCross01 = Cross(v0,v1),
      vCross12 = Cross(v1,v2),
      vCross0112 = Cross(vCross01, vCross12);
   return std::atan2(Dot(vCross0112, v1.Normalized()), Dot(vCross01, vCross12));
   // ^- confused yet? 8)
}


FMeasurePlaneAngle::FMeasurePlaneAngle(int iAt_, int jAt_, int kAt_, int lAt_, FDocument *pDocument_)
   : FMeasure(pDocument_,QString::fromWCharArray(L"\u00b0"),2,180./M_PI), m_iAt(iAt_), m_jAt(jAt_), m_kAt(kAt_), m_lAt(lAt_)
{}

QString FMeasurePlaneAngle::Name() const
{
   return QString(QString::fromWCharArray(L"\u03d5(%1 | %2 | %3 | %4)")).arg(m_pDocument->AtomLabel(m_iAt), m_pDocument->AtomLabel(m_jAt), m_pDocument->AtomLabel(m_kAt), m_pDocument->AtomLabel(m_lAt));
}


void FMeasurePlaneAngle::MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const
{
   FrameValues.clear();
   FrameValues.reserve(Frames.size());
//    foreach(FFrameMeasurePtrList::value_type const *itFrame, Frames) {
   FFrameMeasurePtrList::const_iterator itFrame;
   _for_each(itFrame, Frames) {
      ct::FAtomSet *pAtoms = (*itFrame)->pGetAtoms();
      if (pAtoms && size_t(m_iAt) < pAtoms->size() && size_t(m_jAt) < pAtoms->size()) {
         ct::FVector3
            vI = (*pAtoms)[size_t(m_iAt)].vPos,
            vJ = (*pAtoms)[size_t(m_jAt)].vPos,
            vK = (*pAtoms)[size_t(m_kAt)].vPos,
            vL = (*pAtoms)[size_t(m_lAt)].vPos;
         FrameValues.push_back(m_AuToOut * Dihedral(vI, vJ, vK, vL));
      } else {
         FrameValues.push_back(0);
      }
   }
}




FMeasureFrameGradient::FMeasureFrameGradient(FDocument *pDocument_)
   : FMeasure(pDocument_,"mgu",6,1.)
{
   // mgu: mystery gradient units (different programs write different units into
   // the gradient fields of .xyz files, and particularly what Turbomole was
   // doing... no idea.)
}

void FMeasureFrameGradient::MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const
{
   FrameValues.clear();
   FrameValues.reserve(Frames.size());
//    foreach(FFrameMeasurePtrList::value_type const *itFrame, Frames) {
   FFrameMeasurePtrList::const_iterator itFrame;
   _for_each(itFrame, Frames) {
//       ct::FAtomSet
//          *pAtoms = (*itFrame)->pGetAtoms();
//       double
//          fGrad = 0.;
//       if (pAtoms) {
//          for (size_t iAt = 0; iAt < pAtoms->size(); ++ iAt)
//             fGrad += ct::LengthSq((*pAtoms)[iAt].vGrad);
//       }
//       FrameValues.push_back(sqrt(fGrad));
      FrameValues.push_back((*itFrame)->GetGradient());
   }
}

QString FMeasureFrameGradient::Name() const
{
   return QString::fromUtf8("|grad E|");
}


FMeasureFrameEnergy::FMeasureFrameEnergy(FDocument *pDocument_)
   : FMeasure(pDocument_,"Eh",8,1.)
{
}

void FMeasureFrameEnergy::MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const
{
   FrameValues.clear();
   FrameValues.reserve(Frames.size());
//    foreach(FFrameMeasurePtrList::value_type const *itFrame, Frames) {
   FFrameMeasurePtrList::const_iterator itFrame;
   _for_each(itFrame, Frames) {
      ct::FAtomSet
         *pAtoms = (*itFrame)->pGetAtoms();
      double
         fEnergy = 0.;
      if (pAtoms)
         fEnergy = pAtoms->GetLastEnergy();
      FrameValues.push_back(fEnergy);
   }
}

QString FMeasureFrameEnergy::Name() const
{
   return QString::fromUtf8("Energy");
}


FMeasureFrameId::FMeasureFrameId(FDocument *pDocument_)
   : FMeasure(pDocument_,"",0,1.)
{
}

void FMeasureFrameId::MakeValues(FMeasureSet &FrameValues, FFrameMeasurePtrList &Frames) const
{
   FrameValues.resize(Frames.size());
   for (size_t i = 0; i < Frames.size(); ++ i)
//       FrameValues[i] = double(i);
      FrameValues[i] = m_pDocument->iFrameId(&*Frames[i]);
}

QString FMeasureFrameId::Name() const
{
   return QString("#Frame");
}



FDocumentMeasures::FDocumentMeasures(FDocument *pDocument_, QObject *Parent_)
   : FBase(Parent_), m_pDocument(pDocument_)
{
   AddMeasure(new FMeasureFrameId(m_pDocument));
   AddMeasure(new FMeasureFrameEnergy(m_pDocument));
   AddMeasure(new FMeasureFrameGradient(m_pDocument));

   // note: if this doesn't work, the Q_OBJECT is missing. Would complain about non-existent slots
   // in the base class.
   ConnectForwardSignals();
}

void FDocumentMeasures::ConnectForwardSignals()
{
   // forward those things... to deal with newly inserted frames etc.
   connect(m_pDocument, SIGNAL(layoutAboutToBeChanged()), this, SLOT(parentLayoutAboutToBeChanged()));
   connect(m_pDocument, SIGNAL(layoutChanged()), this, SLOT(parentLayoutChanged()));
}

FDocumentMeasures::~FDocumentMeasures()
{
}


int FDocumentMeasures::rowCount(const QModelIndex &parent) const
{
   return m_pDocument->columnCount();
   IR_SUPPRESS_UNUSED_WARNING(parent);
}

int FDocumentMeasures::columnCount(const QModelIndex &parent) const
{
   return static_cast<int>(m_List.size());
   IR_SUPPRESS_UNUSED_WARNING(parent);
}

QVariant FDocumentMeasures::headerData(int section, Qt::Orientation orientation, int role) const
{
//    IvEmit("HeaderData: Section = %1  orientation = %2  role = %3", section, int(orientation), role);
   if (role == Qt::DisplayRole && orientation == Qt::Vertical) {
      FFrame
         *pFrame = m_pDocument->GetFrame(section, false);
      if (pFrame == 0)
         return QVariant("#UNK");
      return IvFmt("Frame %1", section);
   }
   if (role == Qt::DisplayRole && orientation == Qt::Horizontal) {
      FMeasure const
         *pMeasure = GetMeasure(int(section));
      if (pMeasure == 0)
         return QVariant("#UNK");
      QString
         UnitName = pMeasure->UnitName();
      if (UnitName.isEmpty())
         return pMeasure->Name();
      else
         return IvFmt("%1 [%2]", pMeasure->Name(), UnitName);
   }
   if ( role == Qt::FontRole ) {
         QFont CaptionFont;
         CaptionFont.setBold(true);
         return CaptionFont;
   }
   return QVariant();
}

QVariant FDocumentMeasures::data(const QModelIndex &index, int role) const
{
   if (role == Qt::DisplayRole) {
      FMeasure const
         *pMeasure = GetMeasure(index);
      FFrame
         *pFrame = m_pDocument->GetFrame(index.row(),false);
      if (pFrame == 0 || pMeasure == 0)
         return QVariant("#UNK");
      return QVariant(pMeasure->MeasureFrame(pFrame,FMeasure::FORMAT_ValueOnly));
   }
   if (role == Qt::TextAlignmentRole)
      return QVariant(Qt::AlignRight);
   return QVariant();
}

bool FDocumentMeasures::removeColumns(int column, int count, const QModelIndex &parent)
{
   if (count < 0 || column + count >= columnCount()) {
      IV_NOTIFY(NOTIFY_Warning, "Attempted to remove non-existent measurement rows.");
      return false;
   }
   beginRemoveColumns(parent, column, column + count - 1);
   FMeasureList::iterator
//       itFirst = std::next(m_List.begin(), column),
//       itLast = std::next(itFirst, count);
      itFirst = m_List.begin() + column,
      itLast = itFirst + count;
   m_List.erase(itFirst, itLast);
   endRemoveColumns();
   return true;
}

void FDocumentMeasures::AddMeasure(FMeasurePtr pMeasure)
{
   beginInsertColumns(QModelIndex(), rowCount(), rowCount()); // 'last' is inclusive!.
   m_List.push_back(pMeasure);
   endInsertColumns();
}

FMeasure *FDocumentMeasures::GetMeasure(QModelIndex const &index)
{
   FFrame
      *pFrame = m_pDocument->GetFrame(index.row(),false);
   if (pFrame == 0)
      return 0;
   return GetMeasure(int(index.column()));
}

FMeasure *FDocumentMeasures::GetMeasure(int iColumn)
{
   if (iColumn < 0 || iColumn >= columnCount())
      return 0;
   return &*m_List[size_t(iColumn)];
}

void FDocumentMeasures::parentLayoutAboutToBeChanged()
{
   emit layoutAboutToBeChanged();
}

void FDocumentMeasures::parentLayoutChanged()
{
   emit layoutChanged();
}


FTablesForm::FTablesForm(FDocument *pDocument_, QWidget *pParent_)
   : FBase(pParent_),
     ui(new Ui::TablesForm),
     m_pDocument(pDocument_)
{
   ui->setupUi(this);
   ui->tableView_Values->setModel(m_pDocument->GetMeasures());

#if QT_VERSION >= 0x050000
//    ui->tableView_Values->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
   ui->tableView_Values->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
#else
   // qt5 doc's say there is a compatibiltiy layer for this.. but it doesn't seem to exist
   // in my windows qt 5.3.2 version.
//    ui->tableView_Values->horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
   ui->tableView_Values->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);
#endif
   ui->tableView_Values->resizeColumnsToContents();

//    FTransposeProxyModel
//       *pProxyModel = new FTransposeProxyModel(this);
//    pProxyModel->setSourceModel(m_pDocument->GetMeasures());
//    ui->tableView_Values->setModel(pProxyModel);
}
