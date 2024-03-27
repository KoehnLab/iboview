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

#ifndef IV_SETTINGS_H
#define IV_SETTINGS_H

#include <QWidget>
#include <QSplitter>
#include <QFileDialog>

// all of those routines deal with saving some kinds of user settings in a persistent location.

void IvSaveWindowSize(QString Key, QWidget *pWindow);
bool IvRestoreWindowSize(QString Key, QWidget *pWindow);
void IvSaveSplitterState(QString Key, QSplitter *pSplitter);
bool IvRestoreSplitterState(QString Key, QSplitter *pSplitter);

void IvGuessSubDialogSize(QWidget *pWindow, double fDefaultVerticalScale = -1.);

QString IvGetSaveFileName(QString HistoryName, QWidget *parent = 0, const QString &caption = QString(), const QString &dir = QString(), const QString &filter = QString(), QString *selectedFilter = 0, QFileDialog::Options options = 0);
QStringList IvGetOpenFileNames(QString HistoryName, QWidget *parent = 0, const QString &caption = QString(), const QString &dir = QString(), const QString &filter = QString(), QString *selectedFilter = 0, QFileDialog::Options options = 0);

enum FExecOptions {
   GETEXEC_ShellReplace = 0x01,  // if set, automatically replace environment variables.
   GETEXEC_AssertExists = 0x02   // if set, make a error box if the executable file does not exist.
};
QString GetExecutableName(QString Which, unsigned Flags = GETEXEC_ShellReplace | GETEXEC_AssertExists);
void SetExecutableName(QString Which, QString ExeName);
double GetExternalToolTimeout();
void SetExternalToolTimeout(double fTimeOutInSeconds);


#endif // IV_SETTINGS_H
