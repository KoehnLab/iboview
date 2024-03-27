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

#include "Iv.h"
#include "IvPreferencesForm.h"
#include "IvSettings.h"
// #include "QPropertyModel.h"
#include "ui_PreferencesForm.h"
#include "IvDocument.h" // sigh... that one is only for resetting work space size on accept.
#include <QSettings>
#include <QFileDialog>
#include <QAbstractButton>
#include <QLineEdit>

static bool IsPropertyEqual(QWidget *w, char const *PropertyName, QString ToWhat) {
   QVariant
      vaPropName = w->property(PropertyName);
   if (vaPropName.isValid()) {
      if (vaPropName.toString() == ToWhat)
         return true;
   }
   return false;
};

// find a widget of the templated type which has the dynamic 'ext_tool_property' set to the corresponding value.
template<class TWidget>
TWidget *FPreferencesForm::FindToolWidget(QString Which)
{
   QList<TWidget*> uiWidgets = this->findChildren<TWidget*>();
   foreach(TWidget *w, uiWidgets) {
      if (IsPropertyEqual(w, "ext_tool_filename", Which))
         return w;
   }
   return 0;
};


FPreferencesForm::FPreferencesForm(FDocument *document, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::PreferencesForm),
     m_pDocument(document)
{
   ui->setupUi(this);
   connect(ui->toolButton_ClearScriptFile, SIGNAL(clicked()), this, SLOT(ClearStartupScriptFile()));
   connect(ui->toolButton_FindScriptFile, SIGNAL(clicked()), this, SLOT(SearchStartupScriptFile()));

   QSettings
      settings;
   ui->lineEdit_StartupScriptFile->setText(settings.value("IboView/StartupScriptFile").toString());

//    ui->lineEdit_Executable_tm2molden->setText(GetExecutableName("tm2molden", 0));
//    ui->lineEdit_Executable_orca2_mkl->setText(GetExecutableName("orca2_mkl", 0));
   QList<QWidget*> uiWidgets = this->findChildren<QWidget*>();
   foreach(QWidget *w, uiWidgets) {
      QVariant
         vaToolName = w->property("ext_tool_filename");
      if (vaToolName.isValid()) {
         QString
            sToolName = vaToolName.toString();
         QAbstractButton
            *pButton = qobject_cast<QAbstractButton*>(w);
         QLineEdit
            *pLineEdit = qobject_cast<QLineEdit*>(w);
         if (pButton != 0) {
            // it's the button to 'find' the executable for corresponding tool.
            connect(w, SIGNAL(clicked()), this, SLOT(SearchExternalTool()));
         } else if (pLineEdit != 0) {
            pLineEdit->setText(GetExecutableName(sToolName, 0));
         } else {
            IvNotify(NOTIFY_Error, QString("Coding error: Encountered unexpected widget in preferences dialog (for '%1')").arg(sToolName));
         }
      }
   }
   ui->doubleSpinBox_ExternalToolTimeout->setValue(GetExternalToolTimeout());

   ui->spinBox_WorkspaceMemoryMb->setValue(settings.value("IboView/MemoryStackSize").toInt());
}


FPreferencesForm::~FPreferencesForm()
{
   delete ui;
}


void FPreferencesForm::SearchExternalTool()
{
   QString
      Which;
   if (sender()) {
      QVariant
         vaPropName = sender()->property("ext_tool_filename");
      if (vaPropName.isValid())
         Which = vaPropName.toString();
   }
   if (Which == "") {
      IvNotify(NOTIFY_Error, QString("Coding error: FindExternalTool signal from broken sender (no sender/no ext_tool_name property)"));
   } else {
      // get a pointer to the respective line edit control.
      QLineEdit *lineEdit = FindToolWidget<QLineEdit>(Which);
      if (lineEdit == 0) {
         return IvNotify(NOTIFY_Error, QString("Coding error: Failed to find path control edit for tool '%1'").arg(Which));
      }
      // get a new file name
      QString
#ifdef Q_WS_WIN
         ExeFilter = "Executables (*.exe);;All Files (*.*)";
#else
         ExeFilter = "Executables (*);;All Files (*.*)";
#endif
      QString
         FileName = QFileDialog::getOpenFileName(this, QString("Select Executable For %1").arg(Which),
            lineEdit->text(), ExeFilter);
      if (!FileName.isEmpty())
         lineEdit->setText(FileName);
   }
}



void FPreferencesForm::accept()
{
   QSettings
      settings;
   settings.setValue("IboView/StartupScriptFile", ui->lineEdit_StartupScriptFile->text());
   settings.setValue("IboView/MemoryStackSize", ui->spinBox_WorkspaceMemoryMb->value());
   m_pDocument->GetWfOptions()->SetWorkSpaceMb(ui->spinBox_WorkspaceMemoryMb->value());

//    SetExecutableName("tm2molden", ui->lineEdit_Executable_tm2molden->text());
//    SetExecutableName("orca2_mkl", ui->lineEdit_Executable_orca2_mkl->text());
   QList<QLineEdit*> uiWidgets = this->findChildren<QLineEdit*>();
   foreach(QLineEdit *w, uiWidgets) {
      QVariant
         vaToolName = w->property("ext_tool_filename");
      if (vaToolName.isValid()) {
         QString
            sToolName = vaToolName.toString();
         SetExecutableName(sToolName, w->text());
      }
   }
   SetExternalToolTimeout(ui->doubleSpinBox_ExternalToolTimeout->value());

   return QDialog::accept();
}

void FPreferencesForm::ClearStartupScriptFile()
{
   ui->lineEdit_StartupScriptFile->setText("");
}

void FPreferencesForm::SearchStartupScriptFile()
{
   QString FileName = QFileDialog::getOpenFileName(this, "Select Script File",
      "", "IboView Scripts (*.js);;All Files (*.*)");
   if (!FileName.isEmpty())
      ui->lineEdit_StartupScriptFile->setText(FileName);
}
