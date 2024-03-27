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

#ifndef IV_LOG_H
#define IV_LOG_H

#include <QString>
#include <QTextStream>
#include "CxIo.h"

class FLogQt : public QObject, public ct::FLog
{
   Q_OBJECT
public:
   void Flush(); // override

   explicit FLogQt(QObject *parent=0);
   ~FLogQt();

   FRunStatus GetStatus(); // override
signals:
   void textEmitted(QString s);
   void sectionEnded();
   // ^- that is really just a work-around for QTextBrowser getting really slow when
   //    appending stuff to long texts...
   void reportEnded(QString Log);
public slots:
   // tell the running process to emit an abort signal the next time
   // CheckStatus() is called. Works like this:
   //
   // Gui Thread: send (queued) abort signal to FLogQt. This flags
   //     the AbortSignalled flag in the FLogQt object.
   // Worker thread: call CheckStatus() at convenient point during
   //     calculations. If CheckStatus sees that an abort was
   //     signalled, it raises the corresponding FLogError
   //     exception.
   void emitAbortSignal();
   void endSection();
   virtual void endReport();
protected:
   bool m_AbortSignalled;
};


class FMemoryLogQt : public FLogQt
{
   Q_OBJECT
public:
   explicit FMemoryLogQt(QObject *parent=0);
   ~FMemoryLogQt();

   QString const &GetText();
   void Clear();
   void Flush(); // override
public slots:
   void appendText(QString s);
   void endReport(); // overwrite
protected:
   QString
      m_LogText;
   QTextStream
      m_LogStream;
};


#endif // IV_LOG_H
