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

#ifndef IV_SHOW_TEXT_FORM_H
#define IV_SHOW_TEXT_FORM_H

#include <QDialog>
#include <QString>
#include <QTextStream>
#include <QAbstractButton>
#include <QCloseEvent>
#include "CxIo.h"

namespace Ui{
    class ShowTextForm;
}


class FShowTextForm: public QDialog
{
    Q_OBJECT
public:
    Ui::ShowTextForm *ui;

    FShowTextForm(QString Text_, QString Title_, QString RefFileName_, QString RedButtonText_ = QString(), QWidget *parent = 0);
   ~FShowTextForm();

   QAbstractButton *GetRedButton();
signals:
   void redButtonPressed();
public slots:
   void copyAll();
   void saveText();
   void mysteryAction();
   void clear();
   void setText(QString s);

   void AppendText(QString s);
protected:
   bool m_RedButtonAssigned;
   void closeEvent(QCloseEvent *event); // override.
   QString
      m_RefFileName;
};


#endif // IV_SHOW_TEXT_FORM_H
