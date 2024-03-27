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

#ifndef IV_PREFERENCES_FORM
#define IV_PREFERENCES_FORM

#include <QDialog>
#include <QString>

namespace Ui{
    class PreferencesForm;
}

class FDocument;
class QAbstractButton;
class QLineEdit;

class FPreferencesForm : public QDialog
{
   Q_OBJECT

   Ui::PreferencesForm *ui;
   FDocument
      *m_pDocument;
protected:
   template <class TWidget>
   TWidget *FindToolWidget(QString Which);
public:
   explicit FPreferencesForm(FDocument *document, QWidget *parent=0);
   ~FPreferencesForm();
public slots:
   void ClearStartupScriptFile();
   void SearchStartupScriptFile();
   void accept();
   void SearchExternalTool();
};


#endif // IV_PREFERENCES_FORM
