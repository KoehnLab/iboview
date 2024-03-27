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

#include <QProgressDialog>
#include <QMessageBox>
#include <QtAlgorithms>
#include <QSettings>
#include <QVariant>
#include <fstream>

#include "Iv.h"
#include "IvComputeEosForm.h"
#include "ui_ComputeEosForm.h"
#include "IvDocument.h"
#include "IvSettings.h"



FComputeEosForm::FComputeEosForm(FDocument *document, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::ComputeEosForm),
     m_pDocument(document)
{
   ui->setupUi(this);

//    restoreSavedIsoThresholds();
//
//    connect(ui->toolButton_AddIsoThreshold, SIGNAL(clicked()), this, SLOT(addIsoThreshold()));
//    connect(ui->toolButton_DeleteIsoThreshold, SIGNAL(clicked()), this, SLOT(deleteIsoThreshold()));
//    connect(ui->toolButton_SetIsoColor, SIGNAL(clicked()), this, SLOT(changeSurfaceColor()));
//
//    connect(ui->doubleSpinBox_IsoThreshold, SIGNAL(valueChanged(double)), this, SLOT(changeIsoThreshold(double)));
//    connect(ui->tabWidget_SurfaceType, SIGNAL(currentChanged(int)), this, SLOT(changeSurfaceType(int)));
//
//    changeSurfaceType(0); // select density by default
//
//    ui->buttonBox->setFocus();
   connect(ui->buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
   connect(ui->buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

   if (!IvRestoreWindowSize("ComputeEosForm/Size", this))
      IvGuessSubDialogSize(this);

   FFragmentAnalysisOptionsPtr
      pFragmentOptions = m_pDocument->GetFragmentAnalysisOptions();
   switch(pFragmentOptions->SpacePartitionType) {
      case FFragmentAnalysisOptions::SPACEPARTITION_HilbertSpace_Iao: {
         ui->radioButton_IaoPartition->setChecked(true);
         break;
      }
      case FFragmentAnalysisOptions::SPACEPARTITION_RealSpace_TFVC: {
         ui->radioButton_TfvcPartition->setChecked(true);
         break;
      }
      default:
         IvNotify(NOTIFY_Error, "FComputeEosForm: unrecognized space partitioning state encountered in pFragmentOptions->SpacePartitionType");
   }

   ui->checkBox_GroupRemainingAtoms->setChecked(pFragmentOptions->GroupRemainingAtoms);
   ui->lineEdit_IntegrationGrid->setText(s2q(pFragmentOptions->GridDesc));
}


FComputeEosForm::~FComputeEosForm()
{
   IvSaveWindowSize("ComputeEosForm/Size", this);
   delete ui;
}



void FComputeEosForm::accept()
{
   // remember current settings

   FFragmentAnalysisOptionsPtr
      pFragmentOptions = m_pDocument->GetFragmentAnalysisOptions();
   if (ui->radioButton_IaoPartition->isChecked())
      pFragmentOptions->SpacePartitionType = FFragmentAnalysisOptions::SPACEPARTITION_HilbertSpace_Iao;
   if (ui->radioButton_TfvcPartition->isChecked())
      pFragmentOptions->SpacePartitionType = FFragmentAnalysisOptions::SPACEPARTITION_RealSpace_TFVC;

   pFragmentOptions->GroupRemainingAtoms = ui->checkBox_GroupRemainingAtoms->isChecked();
   pFragmentOptions->GridDesc = q2s(ui->lineEdit_IntegrationGrid->text());

   pFragmentOptions->SaveState();

   // close the dialog.
   FBase::accept();
}
