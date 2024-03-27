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

#ifndef IV_CURVEVIEW_H
#define IV_CURVEVIEW_H

#include "Iv.h"
#include <QGraphicsView>
#include "CxPodArray.h"
#include <vector>


enum FCurveFlags {
   CURVE_Visible = 0x01,
   CURVE_Fill = 0x02,
   CURVE_Outline = 0x04,
   CURVE_Bold = 0x08,
   CURVE_Ellipses = 0x10
};

struct FCurveData {
   FCurveData() {}
   explicit FCurveData(TArray<float> const &Data_, uint32_t Flags_, int zOrder, QPen Stroke_, QBrush Fill_);
   TArray<float>
      Data;
   float GetMin();
   float GetMax();

   QString Title;

   uint32_t Flags;
   QPen Stroke;
   QBrush Fill;
   int zOrder;

   bool isVisible() const { return 0 != (Flags & CURVE_Visible); }
   bool isFilled() const { return 0 != (Flags & CURVE_Fill); }
   bool isOutlined() const { return 0 != (Flags & CURVE_Outline); }
   bool isBold() const { return 0 != (Flags & CURVE_Bold); }
};

enum FHvLineDataType {
   HVCURVE_Hline = 0x00,
   HVCURVE_Vline = 0x01,
   HVCURVE_TypeMask = 0x01,
   HVCURVE_StartArrow = 0x02,
   HVCURVE_EndArrow = 0x04
};

struct FHvLineData {
   float
      fPos;
   uint32_t
      Flags;
   QPen
      Pen;
   FHvLineData(float fPos_, uint32_t Flags_, QPen const &Pen_) : fPos(fPos_), Flags(Flags_), Pen(Pen_) {}
};

class FCurveView : public QGraphicsView
{
    Q_OBJECT

public:
   FCurveView(QWidget *parent = 0);
   ~FCurveView();
   typedef QGraphicsView FBase;

   // lock updates: prevent widget from automatically updating its UI representation
   // on change of data (setCurve/clearCurves). Will only update UI once "UnlockUpdates()"
   // is called.
   void LockUpdates();
   void UnlockUpdates();

   void clearCurves();
   void setCurve(unsigned iCurveId, TArray<float> const &Data, QString Title, uint32_t Flags = CURVE_Outline, int zOrder = 0, QPen const &Stroke = QPen(QColor(0xffffffff)), QBrush const &Fill = QBrush());
   int getNumCurves() { return (int)m_Curves.size(); }
   void getCurve(TArray<float> &Data, QString &Title, uint32_t &Flags, QColor &Color, unsigned iCurveId);
   void addHline(float fPos, QPen const &Pen, uint32_t Flags = 0);
   void addVline(float fPos, QPen const &Pen, uint32_t Flags = 0);
protected:
   void keyPressEvent(QKeyEvent *event); // override
   void wheelEvent(QWheelEvent *event); // override
   void drawBackground(QPainter *painter, const QRectF &rect); // override?
   void resizeEvent(QResizeEvent *event);

   void mousePressEvent(QMouseEvent *event);
   void mouseReleaseEvent(QMouseEvent *event);
   void mouseMoveEvent(QMouseEvent *event);

   void scaleView(qreal scaleFactor);

   std::vector<FCurveData>
      m_Curves;
   std::vector<FHvLineData>
      m_HvLines;
   QTransform
      // transform from data space to window space
      m_TrafoData;
   int
      m_iUpdateLocks;
   bool
      m_UpdateNeeded;

   QRectF getCurveBounds();
   void RebuildScene();
private:
   unsigned LastX, LastY;
};


#endif
