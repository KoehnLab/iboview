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

#ifndef IV_EDIT_VOLUME_SURFACE_FORM
#define IV_EDIT_VOLUME_SURFACE_FORM

#include "CxTypes.h"
#include <QDialog>
#include <QList>
#include <QMetaType>
#include <QItemSelectionModel>
// #include "IvStatusBar.h"

class FIsoSliceDecl {
public:
   FIsoSliceDecl() : fIsoValue(0.), dwColor(0xffffffff) {}
   FIsoSliceDecl(double fIsoValue_, uint32_t dwColor_) : fIsoValue(fIsoValue_), dwColor(dwColor_) {}
   FIsoSliceDecl(FIsoSliceDecl const &other) { fIsoValue = other.fIsoValue; dwColor = other.dwColor; }
   ~FIsoSliceDecl() {}
public:
   double
      fIsoValue;
   uint32_t
      dwColor;
   bool operator < (FIsoSliceDecl const &other) const;

   static FIsoSliceDecl fromString(QString s, bool *ok=0);
   QString toString() const;
};

// for storing this in QVariants (note: must be in header, for moc).
Q_DECLARE_METATYPE(FIsoSliceDecl);


namespace Ui{
    class EditVolumeSurfaceForm;
}

class FDocument;

typedef QList<FIsoSliceDecl>
   FIsoThresholdList;
   // ^- QList because this is what script interface would make for js objects of type [0.001, 0.01, ...].
   // Also, we can put it into a Variant directly.

struct FEditVolumeSurfaceFormImpl;

class FEditVolumeSurfaceForm : public QDialog
{
   Q_OBJECT

   typedef QDialog
      FBase;
   Ui::EditVolumeSurfaceForm
      *ui;
public:
   explicit FEditVolumeSurfaceForm(FDocument *document, QWidget *parent=0);
   ~FEditVolumeSurfaceForm();

   void DoPostProcessingOnClose();
public slots:
   void deleteIsoThreshold();
   void addIsoThreshold();
   void changeIsoThreshold(double newValue);
   void changeSurfaceType(int newIndex);
   void accept();
   void updateWfInformation(int SurfaceType);
   void changeSurfaceColor();
   void currentThresholdChanged(const QModelIndex & current, const QModelIndex & previous);
protected:
   FEditVolumeSurfaceFormImpl
      *p;
   FDocument
      *m_pDocument;
   int
      m_SurfaceType;
   QList<FIsoThresholdList>
      // list of iso-thresholds to use (one separate list for each surface type)
      m_IsoThresholdsForSurfaceTypes;

//    // get threshold list for current surface type
//    FIsoThresholdList GetThresholds();
//    // set threshold list for current surface type
//    void SetThresholds(FIsoThresholdList const &IsoThresholds);
   void updateThresholdListInUi();

   void restoreSavedIsoThresholds();
   // if not -1: store settings only for given surfaceType
   void saveIsoThresholds(int surfaceType=-1);

   void resetVolumePropertyData();
};


#endif // IV_EDIT_VOLUME_SURFACE_FORM
