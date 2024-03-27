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

#include <QTextEdit>
#include <QMessageBox>
#include <QApplication>
#include <QClipboard>
#include "IvShowTextForm.h"
#include "ui_ShowTextForm.h"
#include <QSettings>
#include "IvSettings.h"

FShowTextForm::FShowTextForm(QString Text_, QString Title_, QString RefFileName_, QString RedButtonText_, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::ShowTextForm),
     m_RedButtonAssigned(false),
     m_RefFileName(RefFileName_)
{
   ui->setupUi(this);
//    ui->textControl->setWordWrap(false);
   setWindowTitle(Title_);
   ui->textEdit->setLineWrapMode(QTextEdit::NoWrap);
   ui->textEdit->setText(Text_);
//    ui->textEdit->setPlainText(Text_);
//    AppendText(Text_);

   if (!RedButtonText_.isEmpty()) {
      m_RedButtonAssigned = true;
      ui->pushButton_Mystery->setText(RedButtonText_);
   }
   ui->pushButton_Mystery->setStyleSheet("QPushButton {background-color: #6f0000; color: white;} QPushButton:disabled {background-color: #5f4040; color: #999;}");
   connect(ui->pushButton_CopyToClipboard, SIGNAL(clicked(bool)), this, SLOT(copyAll()));
   connect(ui->pushButton_Mystery, SIGNAL(clicked(bool)), this, SLOT(mysteryAction()));
   connect(ui->pushButton_SaveText, SIGNAL(clicked(bool)), this, SLOT(saveText()));
   setStyleSheet("QTextEdit {background-color: #000; color: #ccc; font-family: monospace;}");

   if (!IvRestoreWindowSize("TextWindow/Size", this))
      IvGuessSubDialogSize(this);
}

void FShowTextForm::closeEvent(QCloseEvent *event)
{
   // doesn't work. closeEvent is never called. Apparently that is intended behavior. No idea.
//    QMessageBox::information(this, "If you can see this, closeEvent() was called.");
//    IvSaveWindowSize("TextWindow/Size", this);
   return QDialog::closeEvent(event);
}


FShowTextForm::~FShowTextForm()
{
//    QMessageBox::information(this, "yay.", "If you can see this, destructor was called.");
   IvSaveWindowSize("TextWindow/Size", this);
}

QAbstractButton *FShowTextForm::GetRedButton()
{
   return ui->pushButton_Mystery;
}


void FShowTextForm::copyAll()
{
   QClipboard *clipboard = QApplication::clipboard();
   clipboard->setText(ui->textEdit->toPlainText());
}

void FShowTextForm::mysteryAction()
{
   if (m_RedButtonAssigned) {
      emit redButtonPressed();
   } else
      QMessageBox::information(this, "Unsolicited Life Advice", "Do not press mysterious red buttons.\nYou never know what they do.");
}

void FShowTextForm::saveText()
{
   QString
      FileNameSuggest = m_RefFileName,
      FileName;
   FileName = IvGetSaveFileName("History/SaveFiles", this, "Save Text",
         (FileNameSuggest), "Text files (*.txt *.log);;All Files (*.*)");
   if (FileName.size() != 0) {
      QFile
         Out(FileName);
      Out.open(QIODevice::WriteOnly | QIODevice::Text);
      Out.write(ui->textEdit->toPlainText().toUtf8());
   }
}

void FShowTextForm::AppendText(QString s)
{
//    if (s.endsWith('\n'))
//       s.chop(1); // append adds one newline automatically.
   s.replace("\n","<br>");
//    s.replace("\a", "<pre style=\"color:white;\">");
//    s.replace("\b", "</pre>");
//    ui->textEdit->append("<pre>"+s+"</pre>");
//    ui->textEdit->append("<pre>"+s+"</pre>");
//    ui->textEdit->append(s);
   // ^- that is probably not very fast...
   ui->textEdit->moveCursor(QTextCursor::End);
   ui->textEdit->textCursor().insertHtml("<pre>" + s + "</pre>");
   //ui->textEdit->textCursor().insertHtml("<pre style=\"color:#ccc;\">" + s + "</pre>");
   ui->textEdit->moveCursor(QTextCursor::End);
}


void FShowTextForm::clear()
{
   ui->textEdit->clear();
}

void FShowTextForm::setText(QString s)
{
   QString x(s);
   x.replace("\n","<br>");
   ui->textEdit->setText(x);
}
