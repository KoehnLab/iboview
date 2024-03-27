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

#include <QMessageBox>
#include <QApplication>
// #include <iostream> // FIXME: remove this.
#include "IvLog.h"


FLogQt::FLogQt(QObject *parent)
   : QObject(parent), m_AbortSignalled(false)
{
}

FLogQt::~FLogQt()
{
}

void FLogQt::Flush()
{
   if (w.size() != 0) {
      emit textEmitted(QString(w.c_str()));
      w.clear();
   }
}

// tell the running process to emit an abort signal the next time
// CheckStatus() is called.
void FLogQt::emitAbortSignal()
{
   m_AbortSignalled = true;
}

ct::FLog::FRunStatus FLogQt::GetStatus()
{
   if (m_AbortSignalled)
      return STATUS_AbortSignalled;
   else
      return STATUS_Okay;
}

void FLogQt::endSection()
{
   emit sectionEnded();
}

void FLogQt::endReport()
{
   emit reportEnded(QString());
}


FMemoryLogQt::FMemoryLogQt(QObject *parent)
   : FLogQt(parent), m_LogText(""), m_LogStream(&m_LogText)
{
}

QString const &FMemoryLogQt::GetText()
{
   return m_LogText;
}

void FMemoryLogQt::Clear()
{
   m_LogText.clear();
}

void FMemoryLogQt::Flush()
{
   if (w.size() != 0) {
      QString
         s(w.c_str());
      emit textEmitted(s);
      m_LogStream << s;
      m_LogStream.flush();
//       std::cerr << w.str();

      w.clear();
   }
}

void FMemoryLogQt::appendText(QString s)
{
//    w << q2s(s);
//    Flush();
   Flush();
   m_LogStream << s;
}


FMemoryLogQt::~FMemoryLogQt()
{
}

void FMemoryLogQt::endReport()
{
   Flush();
   emit reportEnded("<html><pre>" + m_LogText + "</pre></html>");
}
