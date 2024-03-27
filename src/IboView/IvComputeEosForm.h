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

#ifndef IV_COMPUTE_EOS_FORM
#define IV_COMPUTE_EOS_FORM

#include "CxTypes.h"
#include <QDialog>
#include <QList>
#include <QMetaType>
#include <QItemSelectionModel>
// #include "IvStatusBar.h"


namespace Ui{
    class ComputeEosForm;
}

class FDocument;

class FComputeEosForm : public QDialog
{
   Q_OBJECT

   typedef QDialog
      FBase;
   Ui::ComputeEosForm
      *ui;
public:
   explicit FComputeEosForm(FDocument *document, QWidget *parent=0);
   ~FComputeEosForm();

   void DoPostProcessingOnClose();
public slots:
   void accept();
protected:
   FDocument
      *m_pDocument;
};


#endif // IV_COMPUTE_EOS_FORM
