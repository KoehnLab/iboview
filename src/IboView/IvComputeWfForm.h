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

#ifndef IV_COMPUTE_WF_FORM
#define IV_COMPUTE_WF_FORM

#include <QDialog>
#include <QCloseEvent>
#include "IvStatusBar.h"

namespace Ui{
    class ComputeWfForm;
}

class FDocument;

class FComputeWfForm : public QDialog
{
   Q_OBJECT

   Ui::ComputeWfForm *ui;
   FDocument
      *m_pDocument;
public:
   explicit FComputeWfForm(FDocument *document, QWidget *parent=0);
   ~FComputeWfForm();

   bool GetRunScf() const;
   bool GetRunIbba() const;
   bool IsMemoryOkay() const;
public slots:
   void ToggleScfPage(bool Checked);
   void ToggleIbbaPage(bool Checked);
   void RecomputeMemory();
//    void RecomputeWfInfo();
// protected:
//    void closeEvent(QCloseEvent *event); // override
protected:
   void CheckEcps();
   void SetMemoryText(QString s, FStatusClass Class);
   void SetWfInfoText(QString s);
   bool
      m_bMemoryOkay,
      m_bOtherError;
};


#endif // IV_COMPUTE_WF_FORM
