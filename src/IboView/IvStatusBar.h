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

#ifndef IV_STATUS_H
#define IV_STATUS_H

#include <QStatusBar>
#include <QFrame>
#include <QLabel>
#include <QString>

// should log class should also go here? or into a separate file to avoid
// recompiling *everything* whenever the status bar implementation changes?

// this defines the colors used.
enum FStatusClass {
   STATUS_Idle,
   STATUS_Confused,
   STATUS_Working,
   STATUS_WorkingBg,
   STATUS_Warning,
   STATUS_Error,
   STATUS_Unknown
};

// class FStatusBar : public QStatusBar
class FStatusBar : public QFrame
{
   Q_OBJECT
public:
   typedef QFrame FBase;
   explicit FStatusBar(QWidget *parent);

   void SetStatus(FStatusClass Class, QString Text);
//    void Finished();
protected:
   QLabel *m_pStatusText;
   FStatusClass m_LastStatus;
};

QString GetStatusStyle(FStatusClass Class);



#endif // IV_STATUS_H
