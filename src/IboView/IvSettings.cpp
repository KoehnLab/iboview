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

#include <QSize>
#include <QSettings>
#include <QMainWindow>
#include <QMessageBox>
#include <QProcessEnvironment>
#include <QFileInfo>
#include <QDir>
#include <cmath>
#include "Iv.h"
#include "IvSettings.h"

void IvSaveWindowSize(QString Key, QWidget *pWindow)
{
   QSettings
      settings;
   settings.setValue(Key, pWindow->size());
}


bool IvRestoreWindowSize(QString Key, QWidget *pWindow)
{
   // restore size of the dialog if it has been saved before.
   QSettings
      settings;
   QSize
      StoredSize = settings.value(Key).toSize();
   if (StoredSize.isValid()) {
      pWindow->resize(StoredSize);
      // remember that this was not just made up, but actually set up by the user like this.
      // otherwise there is some logic in some places to GUESS a window size
      // (see SetDefaultDialogSize)
      pWindow->setProperty("UserSizeRestored", true);
      return true;
   } else
      return false;
}

void IvSaveSplitterState(QString Key, QSplitter *pSplitter)
{
   QSettings
      settings;
   settings.setValue(Key, pSplitter->saveState());
}

bool IvRestoreSplitterState(QString Key, QSplitter *pSplitter)
{
   QSettings
      settings;
   QByteArray
      StoredSplitterConfig = settings.value(Key).toByteArray();
   if (!StoredSplitterConfig.isEmpty()) {
      pSplitter->restoreState(StoredSplitterConfig);
      return true;
   } else
      return false;
}

extern QMainWindow *g_pMainWindow; // for usage as dialog parent and for getting reference sizes.

void IvGuessSubDialogSize(QWidget *pWindow, double fDefaultVerticalScale)
{
   // make some adjustments to size -- make it relative to parent window size,
   // to better deal with high DPI displays.

   // look for a parent object to take the size of. Either the dialog's real parent, or
   // the main window if it does not have one.
   QWidget
      *pParent = pWindow->parentWidget();
   if (pParent == 0)
      pParent = g_pMainWindow;
   if (pParent) {
      // If the widget has an explicitly set size, do not replace it by our guesswork.
      QVariant
         UserSizeRestored = pWindow->property("UserSizeRestored");
      if (!UserSizeRestored.isValid() || UserSizeRestored.toBool() == false) {
         if (fDefaultVerticalScale < 0.)
            fDefaultVerticalScale = 0.8;

         QSize
            sz = pWindow->size();
         if (sz.height() > 0) {
            double
               fv = fDefaultVerticalScale*double(pParent->height())/double(sz.height()),
               fh = std::pow(fv,0.75);

            // do some further adjustments to prevent making windows wider than their parents.
            double
               fDefaultHorizontalScale = 0.7;
            if (fh*double(sz.width()) > fDefaultHorizontalScale * double(pParent->width())) {
               fh = fDefaultHorizontalScale*double(pParent->width())/double(sz.width());
               fv = fh;
            }

            if (fh > 1. && fv > 1.) {
               sz.setHeight(int(sz.height() * fv));
               sz.setWidth(int(sz.width() * fh));
               pWindow->resize(sz);
            }
         }
      }
   }
}


QStringList IvMakeFileDialogWithHistory(QFileDialog::FileMode FileMode, QFileDialog::AcceptMode AcceptMode, QFileDialog::ViewMode ViewMode,  QString HistoryName, QWidget *parent, const QString &caption, const QString &dir, const QString &filter, QString *selectedFilter, QFileDialog::Options options)
{
   QStringList
      FileNames;
   QFileDialog
      FileDialog(parent, caption, QString(), filter);
//       FileDialog(parent, caption, QString()); // <- no filter.
//       ("./"+dir)
//    FileDialog.setFileMode(QFileDialog::AnyFile);
//    FileDialog.setAcceptMode(QFileDialog::AcceptSave);
//    FileDialog.setViewMode(QFileDialog::Detail);
   FileDialog.setFileMode(FileMode);
   FileDialog.setAcceptMode(AcceptMode);
   FileDialog.setViewMode(ViewMode);
   FileDialog.setOptions(options);
   FileDialog.selectFile(dir);
   QSettings
      settings;
   QStringList
      History;
   if (!HistoryName.isEmpty()) {
      History = settings.value(HistoryName).toStringList();
      if (!History.empty()) {
         FileDialog.setHistory(History);
//          QMessageBox::information(0, "history-in", History.join("\n"));
      }
      // ^- this doesn't work... history does not appear anywhere.
      // my QT version seems to have some other strange behaviors. Like not changing to directories
      // when directories are selected, but instead asking me to overwrite them. Even though the code in
      // qfiledialog.cpp clearly checks for this case and looks like it does the right thing
      // (update: that's because it uses the native dialog... and just emits QDialog::accept() in this case. Gna.).
      // Maybe this could be hacked via the filesSelected signal, but that could just as well make everything
      // worse...
      // UPDATE: doesn't even work. only called after the dialog closes. too late to change at that point.
   }
//    FSaveFileDirChangeProxy
//       Dummy(&FileDialog);
   QString
      LastPath = settings.value(HistoryName+"Path").toString();
   if (!LastPath.isEmpty()) {
      FileDialog.setDirectory(LastPath);
//       QMessageBox::information(0, "last path / " + HistoryName+"Path", LastPath);
   }
   if (FileDialog.exec()) {
      FileNames = FileDialog.selectedFiles();

      if (!HistoryName.isEmpty()) {
         History = FileDialog.history(); // <- doesn't seem to add files automatically.
         QString
            FirstPath;
         for (int i = 0; i < FileNames.size(); ++ i) {
            QFileInfo
               info(FileNames[i]);
            if (i == 0) FirstPath = info.path();
            History.insert(0, info.path());
         }
         History.removeDuplicates();
         while (History.size() > 5)
            History.pop_back(); // QList has no resize method.

         settings.setValue(HistoryName, History);
         if (!FirstPath.isEmpty())
            settings.setValue(HistoryName+"Path", FirstPath);
      }
//          QMessageBox::information(0, "history-out", History.join("\n"));
      if (selectedFilter)
         *selectedFilter = FileDialog.selectedNameFilter();
   }
   return FileNames;
}

QString IvGetSaveFileName(QString HistoryName, QWidget *parent, const QString &caption, const QString &dir, const QString &filter, QString *selectedFilter, QFileDialog::Options options)
{
   QStringList
      FileNames = IvMakeFileDialogWithHistory(QFileDialog::AnyFile, QFileDialog::AcceptSave, QFileDialog::Detail, HistoryName, parent, caption, dir, filter, selectedFilter, options);
   if (FileNames.size() == 1)
      return FileNames[0];
   // cancelled. Note: In mode QFileDialog::AnyFile theoretically the dialog should never return more than one selected file.
   return QString();
}

QStringList IvGetOpenFileNames(QString HistoryName, QWidget *parent, const QString &caption, const QString &dir, const QString &filter, QString *selectedFilter, QFileDialog::Options options)
{
   return IvMakeFileDialogWithHistory(QFileDialog::ExistingFiles, QFileDialog::AcceptOpen, QFileDialog::Detail, HistoryName, parent, caption, dir, filter, selectedFilter, options);
}

static QString GetDefaultExecutableName(QString Which)
{
   QString
      s;
   if (Which == "tm2molden")
      s = QString("$(TURBODIR)/bin/$(TURBOMOLE_SYSNAME)/tm2molden");
   else
      s = Which;
#ifdef Q_WS_WIN
   s = s + ".exe";
#endif
   return s;
}

QString ReplaceEnvironmentVariables(QString in) {
   QString
      out = in;
   QProcessEnvironment const
      &env = QProcessEnvironment::systemEnvironment();
   QStringList
      keys = env.keys();
   foreach(QString const &key, keys) {
//       IvEmit("REV: %1 -> %2", QString("$(%1)").arg(key), env.value(key));
      out = out.replace(QString("$(%1)").arg(key), env.value(key));
   }
//    IvEmit("REV: -> %1", out);
   out = QDir::cleanPath(out);
   return out;
}


void SetExecutableName(QString Which, QString ExeName)
{
   if (ExeName != GetDefaultExecutableName(Which)) {
      QSettings
         settings;
       settings.setValue(QString("ExternalTools/%1").arg(Which), ExeName);
   }
}


QString GetExecutableName(QString Which, unsigned Flags)
{
   QString
      ExeName;
   {
      QSettings
         settings;
      QVariant
         v = settings.value(QString("ExternalTools/%1").arg(Which));
      if (v.isNull() || !v.isValid() || !v.canConvert<QString>() ) {
         ExeName = GetDefaultExecutableName(Which);
      } else {
         ExeName = v.toString();
      }
   }

   if (0 != (Flags & GETEXEC_ShellReplace)) {
      ExeName = ReplaceEnvironmentVariables(ExeName);
   }

   if (0 != (Flags & GETEXEC_AssertExists)) {
      QFileInfo
         info(ExeName);
      if (!info.exists() || !info.isFile()) {
         IvNotify(NOTIFY_Error, QString("While invoking '%2': Executable '%1' does not exist (or is not a file). Please use Edit/Preferences/Tools to setup path for %2.").arg(ExeName).arg(Which));
         ExeName = "";
      } else if (!info.isExecutable()) {
         IvNotify(NOTIFY_Error, QString("While invoking '%2': File '%1' exists, but is not executable. Permission problems?").arg(ExeName).arg(Which));
         ExeName = "";
      }
   }
//    IvEmit("Searching %1: Return '%2'", Which, ExeName);
   return ExeName;
}


double GetExternalToolTimeout()
{
   QSettings
      settings;
   QVariant
      v = settings.value(QString("ExternalTools/TimeoutInSeconds"));
   if (v.isNull() || !v.isValid() || !v.canConvert<double>() ) {
      return -1.;
   } else {
      bool
         ok;
      double
         f = v.toDouble(&ok);
      if (!ok)
         return -1.;
      else
         return f;
   }
}

void SetExternalToolTimeout(double fTimeOutInSeconds)
{
   QSettings
      settings;
   settings.setValue(QString("ExternalTools/TimeoutInSeconds"), fTimeOutInSeconds);
}



