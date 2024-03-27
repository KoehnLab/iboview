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

#include <QDir>
#include <QFile>
#include <QTemporaryFile>
#include <QMainWindow>

#include "Iv.h"
#include "IvFileConvert.h"
#include "IvSettings.h"
#include "IvShowTextForm.h"


int GetTimeoutMs()
{
   double
      fTimeOut = GetExternalToolTimeout();
   int
      nTimeoutMs;
   if (fTimeOut <= 0)
      nTimeoutMs = 20000;
   else
      nTimeoutMs = int(fTimeOut * 1e3 + 0.5);
   return nTimeoutMs;
}


FTm2MoldenSubProcess::FTm2MoldenSubProcess(QString InputFileName, QObject *parent)
   : QProcess(parent), m_pTempFile(0)
{
   QFileInfo
      FileInfo(InputFileName);
   setProcessChannelMode(QProcess::MergedChannels);
   setWorkingDirectory(FileInfo.absolutePath());
}

FTm2MoldenSubProcess::~FTm2MoldenSubProcess()
{
   delete m_pTempFile;
}

bool FTm2MoldenSubProcess::Run()
{
   QString
      ConversionExe = GetExecutableName("tm2molden");
   if (ConversionExe == "")
      return false;

   m_pTempFile = new QTemporaryFile(QDir::cleanPath(QDir::tempPath() + QDir::separator() + "tm2_XXXXXX.molden"));
   if (!m_pTempFile->open()) {
      IvNotify(NOTIFY_Error, IvFmt("Failed to open temporary file '%1'.", m_pTempFile->fileName()));
      return false;
   }

//       // invoked when tm2molden says something (it does not appear to have a batch mode).
//       connect(this, SIGNAL(readyRead()), this, SLOT(readText()));
   connect(this, SIGNAL(started()), this, SLOT(writeInputs(void)));
   start(ConversionExe);
   if (!waitForFinished(GetTimeoutMs())) {
      IV_NOTIFY(NOTIFY_Error, "tm2molden conversion process timed out.");
      return false;
   } else {
      QByteArray sOutput = readAllStandardOutput();
      if (sOutput.contains("tm2molden ended normally")) {
         return true;
      } else {
         QString s = "tm2molden failed. Output was:\n" + QString(sOutput);
         FShowTextForm
            ShowTextForm("<pre>" + s + "</pre>", "tm2molden failed.", QString(), QString(), g_pMainWindow);
         ShowTextForm.exec();
//             IV_NOTIFY(NOTIFY_Error, "tm2molden failed. Output was:\n");
         return false;
      }
   }
   return false;
}


QString FTm2MoldenSubProcess::ConvertedFileName() {
      return (m_pTempFile != 0)? m_pTempFile->fileName() : "";
}


void FTm2MoldenSubProcess::writeInputs() {
   write(q2s(ConvertedFileName()+"\n\n\n\n\n\n\n\n\n\n\n\n").c_str());
}


FOrca2MklSubProcess::FOrca2MklSubProcess(QString InputFileName, QObject *parent)
   : QProcess(parent), m_InputFileName(InputFileName), m_pTempFileGbw(0), m_pTempFileMolden(0)
{
}

FOrca2MklSubProcess::~FOrca2MklSubProcess()
{
//    delete m_pTempFileMolden;
//    delete m_pTempFileGbw;
   if (m_FileNameGbw != "" && QFileInfo(m_FileNameGbw).exists())
      QFile(m_FileNameGbw).remove();
   if (m_FileNameMolden != "" && QFileInfo(m_FileNameMolden).exists())
      QFile(m_FileNameMolden).remove();
}

bool FOrca2MklSubProcess::Run()
{
   m_pTempFileGbw = new QTemporaryFile(QDir::cleanPath(QDir::tempPath() + QDir::separator() + "o2m_XXXXXX.gbw"), this);
   // copy input to temp dir (because there seems to be no way to tell orca_2mkl
   // where to copy its stuff -- it makes names automatically)
   m_pTempFileGbw->open();
   m_FileNameGbw = m_pTempFileGbw->fileName();
   m_pTempFileGbw->remove();
   {
      QFile
         InputFile(m_InputFileName);
      InputFile.copy(m_FileNameGbw);
   }
   QString
      BaseName;
   {
      QFileInfo
         info(m_FileNameGbw);
      BaseName = QDir::cleanPath(info.absolutePath() + QDir::separator() + info.completeBaseName());
      if (!info.exists() || !info.isReadable()) {
         IvNotify(NOTIFY_Error, IvFmt("Failed to copy '%1' to '%2'. Dest not readable.", m_InputFileName, info.absoluteFilePath()));
         return false;
      }
   }
   // that's the name orca2_mkl generates, from the basename. (without the .gbw!).
   // if the .gbw is provided, it will not do anything.
   m_FileNameMolden = BaseName + ".molden.input";
   IvEmit("attempted conversion: '%1' -> '%2'", m_FileNameGbw, m_FileNameMolden);

   QStringList
      args;
   args << BaseName;
   args << "-molden";
   start(GetExecutableName("orca_2mkl"), args);
   if (!waitForFinished(GetTimeoutMs())) {
      IV_NOTIFY(NOTIFY_Error, "orca_2mkl conversion process timed out.");
      return false;
   } else {
      QByteArray sOutput = readAllStandardOutput();
      if (exitCode() != 0) {
         QString s = "orca_2mkl failed. Output was:\n" + QString(sOutput);
         FShowTextForm
            ShowTextForm("<pre>" + s + "</pre>", "orca_2mkl failed.", QString(), QString(), g_pMainWindow);
         ShowTextForm.exec();
         return false;
      }
   }
   return true;
}

QString FOrca2MklSubProcess::ConvertedFileName()
{
   return m_FileNameMolden;
//    if (m_pTempFileGbw) {
//       QFileInfo
//          info(m_pTempFileGbw->fileName());
//    } else
//       return "";
}
