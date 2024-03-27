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

#ifndef IV_FILECONVERT_H
#define IV_FILECONVERT_H

#include <QProcess>
#include <QString>

class QTemporaryFile;

class FTm2MoldenSubProcess : public QProcess
{
   Q_OBJECT
public:
   explicit FTm2MoldenSubProcess(QString InputFileName, QObject *parent=0);
   ~FTm2MoldenSubProcess();

   bool Run();

   QString ConvertedFileName();
protected:
   QTemporaryFile
      *m_pTempFile;
public slots:
   void writeInputs();
};


class FOrca2MklSubProcess : public QProcess
{
   Q_OBJECT
public:
   explicit FOrca2MklSubProcess(QString InputFileName, QObject *parent=0);
   ~FOrca2MklSubProcess();

   bool Run();
   QString ConvertedFileName();
protected:
   QString
      m_InputFileName,
      m_FileNameGbw,
      m_FileNameMolden;
   QTemporaryFile
      *m_pTempFileGbw,
      *m_pTempFileMolden;
};


#endif // IV_FILECONVERT_H
