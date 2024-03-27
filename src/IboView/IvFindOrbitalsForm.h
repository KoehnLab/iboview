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

#ifndef IV_FIND_ORBITALS_FORM
#define IV_FIND_ORBITALS_FORM

#include <QDialog>
#include <QAbstractTableModel>
#include "IvDocument.h"

namespace Ui{
    class FindOrbitalsForm;
}


struct FFoundOrbital {
   double fCharacter; // scalar criterion to search orbitals by, e.g., charge on selected atoms (if finding orbitals on atoms), or maximum orbital change (if finding orbitals reacting along a reaction path).
   int iDataRow; // dataset number.
   FOrbital *pOrbital;
   FFoundOrbital(double fCharacter_, int iDataRow_, FOrbital *pOrbital_): fCharacter(fCharacter_), iDataRow(iDataRow_), pOrbital(pOrbital_) {}
};
typedef std::vector<FFoundOrbital>
   FFoundOrbitalList;

// sort items in L by value in FFoundOrbital.fCharacter
void SortFoundOrbitalList(FFoundOrbitalList &L, bool Descending);

// represent a selection of orbitals found via "Find Orbitals"
// on selected atoms/ordered by change along reaction path.
class FFoundOrbitalModel : public QAbstractTableModel
{
   Q_OBJECT
public:
   FFoundOrbitalModel(FFoundOrbitalList const &L, QString CharacterName_, FDocument *pDocument, QObject *parent=0);

   int rowCount(const QModelIndex &parent = QModelIndex()) const;
   int columnCount(const QModelIndex &parent = QModelIndex()) const;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;

public slots:
   void documentDataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight);
public:
   FDocument
      *m_pDocument;
   FFoundOrbitalList
      m_FoundOrbitals;
   QString
      // name to use for criterion in table caption (e.g., "Sel.Chg. or Max.Change")
      m_CharacterName;
};



class FFindOrbitalsForm : public QDialog
{
    Q_OBJECT
public:
    Ui::FindOrbitalsForm *ui;
    FFoundOrbitalModel *m_pModel;

    FFindOrbitalsForm(FFoundOrbitalModel *pModel_, QString Title_, QWidget *parent = 0);
   ~FFindOrbitalsForm();
public slots:
   void toggleRow(const QModelIndex &index);
   void toggleSelectedRows();
};



#endif // IV_FIND_ORBITALS_FORM
