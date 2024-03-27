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

#include <QApplication>
#include <QStyle>
#include <QSizeGrip>
#include "IvStatusBar.h"
#include <QLayout>
#include <QHBoxLayout>
#include <QSizePolicy>
#include <QEventLoop>

// QColor GetStatusColor(FStatusClass Class)
// {
//    // these are some pure colors from http://www.w3schools.com/tags/ref_colorpicker.asp
//    switch (Class) {
//       case STATUS_Idle: return QRgb(0x0066CC);
//       case STATUS_Working: return QRgb(0x009900);
//       case STATUS_WorkingBg: return QRgb(0x003300);
//       case STATUS_Warning: return QRgb(0xFF9900);
//       case STATUS_Error: return QRgb(0xCC0000);
//       default: return QRgb(0x000000);
//    }
// }

QString GetStatusStyle(FStatusClass Class)
{
   // these are some pure colors from http://www.w3schools.com/tags/ref_colorpicker.asp
   // hm... maybe I take the android default colors? http://developer.android.com/design/style/color.html
   //   seems to be really close to what I chose. Wow, I have style 8).
//    font-weight: normal;
#define BGSTYLE(a) "color: #ffffff; background: " a
   switch (Class) {
      case STATUS_Idle: return BGSTYLE("#0066CC");
//       case STATUS_Working: return BGSTYLE("#009900");
//       case STATUS_WorkingBg: return BGSTYLE("#003300");
//       case STATUS_Warning: return BGSTYLE("#ff9900");
//       case STATUS_Error: return BGSTYLE("#cc0000");
//       default: return BGSTYLE("#000000");
//       case STATUS_Idle: return BGSTYLE("#0099CC"); // bright: 33b5e5
      case STATUS_Working: return BGSTYLE("#669900");// bright: 99cc00
      case STATUS_WorkingBg: return BGSTYLE("#003300");
      case STATUS_Warning: return BGSTYLE("#ff8800"); // bright: ffbb33
      case STATUS_Error: return BGSTYLE("#cc0000");
      case STATUS_Confused:
      default: return BGSTYLE("#9933cc"); // bright: aa66cc
   }
#undef BGSTYLE
}

QString GetStatusStyle1(FStatusClass Class)
{
   return QString("FStatusBar{margin: -10px; %1}").arg(GetStatusStyle(Class));
//    // these are some pure colors from http://www.w3schools.com/tags/ref_colorpicker.asp
//    // hm... maybe I take the android default colors? http://developer.android.com/design/style/color.html
//    //   seems to be really close to what I chose. Wow, I have style 8).
// #define BGSTYLE(a) "FStatusBar{margin: -10px; color: #ffffff; font-weight: normal; background: " a "}"
//    switch (Class) {
//       case STATUS_Idle: return BGSTYLE("#0066CC");
// //       case STATUS_Working: return BGSTYLE("#009900");
// //       case STATUS_WorkingBg: return BGSTYLE("#003300");
// //       case STATUS_Warning: return BGSTYLE("#ff9900");
// //       case STATUS_Error: return BGSTYLE("#cc0000");
// //       default: return BGSTYLE("#000000");
// //       case STATUS_Idle: return BGSTYLE("#0099CC"); // bright: 33b5e5
//       case STATUS_Working: return BGSTYLE("#669900");// bright: 99cc00
//       case STATUS_WorkingBg: return BGSTYLE("#003300");
//       case STATUS_Warning: return BGSTYLE("#ff8800"); // bright: ffbb33
//       case STATUS_Error: return BGSTYLE("#cc0000");
//       case STATUS_Confused:
//       default: return BGSTYLE("#9933cc"); // bright: aa66cc
//    }
// //    switch (Class) {
// //       case STATUS_Idle: return "QStatusBar{color: #ffffff; background: #0066CC}";
// //       case STATUS_Working: return "QStatusBar{color: #ffffff; background: #009900}";
// //       case STATUS_WorkingBg: return "QStatusBar{color: #ffffff; background: #003300}";
// //       case STATUS_Warning: return "QStatusBar{color: #ffffff; background: #ff9900}";
// //       case STATUS_Error: return "QStatusBar{color: #ffffff; background: #cc0000}";
// //       default: return "QStatusBar{color: #ffffff; background: #000000}";
// //    }
}


FStatusBar::FStatusBar(QWidget *parent)
   : FBase(parent), m_LastStatus(STATUS_Unknown)
{
   // hm... problem: since we use style sheets, all the palette stuff
   // is ignored, since the inherited cascading style sheet gets preference.
   // Setting an empty style sheet does not help.
   // It seems that the only way to change the appearance of a styled widget is
   // to set a new style sheet. Maybe that is the wrong approach altogether, and
   // I should use a graphicsView object? or make an own paint function? (might
   // need that for the progress indicators anyway...)
//    setStyleSheet(GetStatusStyle(STATUS_Working));
//    layout()->setContentsMargins(0,-4,0,-6);

   QHBoxLayout *pLayout = new QHBoxLayout;
   // hm.. something is off. Is there something which sets negative margins?
   pLayout->setContentsMargins(15,0,15,0);

   m_pStatusText = new QLabel(this);
   pLayout->addWidget(m_pStatusText);
   pLayout->addWidget(new QSizeGrip(this));
   setLayout(pLayout);
//    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

   SetStatus(STATUS_Working, "Starting");
}

void FStatusBar::SetStatus(FStatusClass Class, QString Text)
{
//    QPalette palette = this->palette();
//    palette.setColor( this->foregroundRole(), Qt::white );
//    palette.setColor( this->backgroundRole(), QColor(0xf33333) );
//    palette.setColor( QPalette::Background, QColor(0xf3ff33) );
//    palette.setColor( QPalette::Window, QColor(0xf3fff3) );
//    this->setPalette(palette);
//    this->setAutoFillBackground(true);

   if (m_LastStatus != Class) {
      setStyleSheet(GetStatusStyle1(Class));
      m_LastStatus = Class;
   }
//    showMessage(Text);
   m_pStatusText->setText(Text);
//    if (Class == STATUS_Working) {
//       for (int iDummy = 0; iDummy < 100; ++ iDummy)
//          QCoreApplication::processEvents();
//          QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
//       // FIXME: ^- this is a REAALLLY bad idea. Could lead to re-entrant calls
//       // of all kinds of functions which are not expecting this.
//    }
}

// void FStatusBar::Finished()
// {
//    SetStatus(STATUS_Idle, "Ready");
// }
