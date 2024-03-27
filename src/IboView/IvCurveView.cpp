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

#include "IvCurveView.h"
#include <QMouseEvent>
#include <QColor>
#include <iostream>
#include <QGraphicsPathItem>

FCurveView::FCurveView(QWidget *parent)
   : QGraphicsView(parent)
{
   m_iUpdateLocks = 0;
   m_UpdateNeeded = true;

   setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
   setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

   QGraphicsScene *scene = new QGraphicsScene(this);
   // ^- FIXME: should I delete this or does this belong to the view? I think it's mine...
   scene->setItemIndexMethod(QGraphicsScene::NoIndex);
//     scene->setSceneRect(0, 0, 1000, 1000);
//     scene->setSceneRect(0, 0, width(), height());
   setScene(scene);

//    setCacheMode(CacheBackground);
   // ^- doesn't work with widget-fixed backgrounds :(.
//    setBackgroundBrush(QBrush(QRgb(0xff404040))); // dark gray
   setViewportUpdateMode(BoundingRectViewportUpdate);
   setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
                  QPainter::HighQualityAntialiasing | QPainter::SmoothPixmapTransform);

//    setTransformationAnchor(AnchorUnderMouse);
// //     scale(qreal(0.8), qreal(0.8));
// //     setMinimumSize(400, 400);
//    setAlignment(0);
// //    setAlignment(0);
//    setDragMode(ScrollHandDrag);
//
//
//    setResizeAnchor(AnchorViewCenter);
//    setInteractive(true);
//    setTransformationAnchor(AnchorUnderMouse);
// //    setDragMode(QGraphicsView::RubberBandDrag);
   setDragMode(QGraphicsView::ScrollHandDrag);


//    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
//    setMinimumHeight(100);
//    setMaximumHeight(100);

//    {
//       TArray<float> hmpf;
//       for (unsigned i = 0; i < 50; ++ i)
//          hmpf.push_back(rand() % 40);
//       setCurve(0, hmpf, CURVE_Fill | CURVE_Outline, 0xffffffff, 0xff703030);
//    }
}

FCurveView::~FCurveView()
{
   delete scene();
   // ^- FIXME hmhmhm...
}


void FCurveView::LockUpdates()
{
   assert(m_iUpdateLocks >= 0);
   m_iUpdateLocks += 1;
}


void FCurveView::UnlockUpdates()
{
   m_iUpdateLocks -= 1;
   assert(m_iUpdateLocks >= 0);
   if (m_iUpdateLocks == 0)
      // done. Do the scene rebuild now.
      RebuildScene();
}


QRectF FCurveView::getCurveBounds()
{
   if (m_Curves.empty())
      return QRectF(0., 0., 50, 1.);
   unsigned
      nPts = 0;
   float
      Min=1e9, Max=-1e9;
   for (unsigned i = 0; i < m_Curves.size(); ++ i) {
      nPts = std::max(nPts, unsigned(m_Curves[i].Data.size()));
      Min = std::min(Min, m_Curves[i].GetMin());
      Max = std::max(Max, m_Curves[i].GetMax());
   }
   return QRectF(0, Min, nPts, Max-Min);
}

// this converts a QRectF to a QPolygonF... I *think* QPolygonF::QPolygonF(QRectF)
// should do exactly the same, but for some reason the thus constructed QPolygonF
// objects do not work in quatToQuad on my machine.
static QPolygonF rectToPoly(QRectF const &rect) {
   QPolygonF
      poly;
   poly << rect.topLeft()
        << rect.topRight()
        << rect.bottomRight()
        << rect.bottomLeft();
   // ^- note: this does NOT close the polygon. Apparently this is the problem
   //    with the QRectF constructor... it makes a polygon with 5 vertices.
   return poly;
}

static const double Pi = 3.14159265358979323846264338327950288419717;
// static double TwoPi = 2.0 * Pi;

static void Orth2(QPointF &u, QPointF &v, QPointF const &dxy)
{
   float f = 1.f/std::sqrt(sqr(dxy.x()) + sqr(dxy.y()));
   u = f * dxy;
   v = QPointF(u.y(), -u.x());
}

static void AddArrowToPath(QPainterPath *path, QPointF p0, QPointF p1, float ArrowAngle, float ArrowSize)
{
   QPointF
      u, v;
   Orth2(u, v, p0-p1);

   float
      du = ArrowSize,
      dv = -std::tan(ArrowAngle/2) * du;

   QPolygonF
      Arrow;
   Arrow << p0 + du * u
         << p0 + dv * v
         << p0 - dv * v;
   path->addPolygon(Arrow);
   path->closeSubpath();
}

static void DrawLineWithArrows(QGraphicsScene *scene, QLineF const &line, QPen const &pen, bool bStartArrow, bool bEndArrow)
{
   float
//       LineAngle = std::atan2(line.dy(), line.dx()),
      ArrowAngle = Pi/3,
      ArrowSize = 3. * pen.widthF();
   {
      QPainterPath
         path;
      path.moveTo(line.p1());
      path.lineTo(line.p2());
      scene->addPath(path, pen, QBrush());
   }

   if (bStartArrow || bEndArrow) {
      QPainterPath
         path;
      if (bStartArrow)
         AddArrowToPath(&path, line.p1(), line.p2(), ArrowAngle, ArrowSize);
      if (bEndArrow)
         AddArrowToPath(&path, line.p2(), line.p1(), ArrowAngle, ArrowSize);
      scene->addPath(path, QPen(Qt::NoPen), pen.brush());
   }
}


void FCurveView::RebuildScene()
{
   m_UpdateNeeded = true;
   if (m_iUpdateLocks > 0)
      // do not do actual updates of the scene at this point;
      // wait until data is complete. Otherwise O(N^2) scaling in some cases...
      return;

   scene()->clear();
   QRectF
      cb = getCurveBounds();
   float
//       fMarginH = 0.03 * cb.width(),
//       fMarginV = 0.08 * cb.height();
      fMarginH = 10.f, // pixels.
      fMarginV = 10.f;
   bool
//       bTransWorked = QTransform::quadToQuad(QPolygonF(getCurveBounds()), QPolygonF(QRectF(this->rect())), m_TrafoData);
//       bTransWorked = QTransform::quadToQuad(QPolygonF(QRectF(cb.left(),0.f,cb.width(),1.f)), QPolygonF(QRectF(this->rect())), m_TrafoData);
//       bTransWorked = QTransform::quadToQuad(rectToPoly(QRectF(cb.left()-fMarginH,1.f+fMarginV,cb.width()+2*fMarginH,-1.f-2*fMarginV)), rectToPoly(QRectF(0,0,width(),height())), m_TrafoData);
      bTransWorked = QTransform::quadToQuad(rectToPoly(QRectF(cb.left(),1.f,cb.width(),-1.f)), rectToPoly(QRectF(fMarginH,fMarginV,width()-2*fMarginH,height()-2*fMarginV)), m_TrafoData);
   if (!bTransWorked) {
      std::cerr << "\nWARNING: Failed to set up data transformation in FCurveView.\n" << std::endl;
   }

#define DMAP(x,y) m_TrafoData.map(QPointF((x),(y)))
   for (uint iLine = 0; iLine < m_HvLines.size(); ++ iLine) {
      FHvLineData
         &HvLine = m_HvLines[iLine];
      QPainterPath path;
      QLineF line;
      if ((HvLine.Flags & HVCURVE_TypeMask) == HVCURVE_Hline) {
         line = QLineF(QPointF(DMAP(cb.left(), HvLine.fPos)),
                       QPointF(DMAP(cb.right(), HvLine.fPos)));
      } else {
         line = QLineF(QPointF(DMAP(HvLine.fPos, 0.f)),
                       QPointF(DMAP(HvLine.fPos, 1.f)));
      }
//       scene()->addPath(path, HvLine.Pen, QBrush());
      DrawLineWithArrows(scene(), line, HvLine.Pen, bool(HvLine.Flags & HVCURVE_StartArrow), bool(HvLine.Flags & HVCURVE_EndArrow));
   }

   for (uint iCurveId = 0; iCurveId < m_Curves.size(); ++ iCurveId)
   {
      FCurveData
         &Curve = m_Curves[iCurveId];

      QPen &pen = Curve.Stroke;
      QBrush &brush = Curve.Fill;

      if (Curve.isBold())
         pen.setWidthF(3.00);

      // see also: /opt/prog/qt-everywhere-opensource-src-4.8.5/examples/graphicsview/elasticnodes
      {
         QPainterPath path;

         if (Curve.Flags & CURVE_Ellipses) {
            size_t nPt = Curve.Data.size();
            for (unsigned i = 0; i < Curve.Data.size(); ++ i) {
               QPointF vCen = DMAP(i, Curve.Data[i]);
               if (nPt < 100) {
                  path.addEllipse(DMAP(i, Curve.Data[i]), 2., 2.);
                  // ^- this is pretty, but makes things *REAAAAALLLLLLYYYY* slow for long IRCs...
                  //    So if there are many points, we will first draw rects instead of ellipses
                  //    (these render much faster), and at some point nothing at all (won't be visually
                  //    distinguishable anyway).
               } else if (nPt < 500) {
                  qreal r = 1.8;
                  path.addRect(vCen.x(), vCen.y(), r, r);
               } else {
                  // ignore request for dots if too many points.
               }
            }
         }

         for (unsigned i = 0; i < Curve.Data.size(); ++ i)
            if (i == 0)
               path.moveTo(DMAP(i, Curve.Data[i]));
            else
               path.lineTo(DMAP(i, Curve.Data[i]));
         if (Curve.isOutlined())
            scene()->addPath(path, pen, QBrush())->setZValue(Curve.zOrder);
         if (Curve.isFilled()) {
            path.lineTo(DMAP(Curve.Data.size()-1, 0.));
            path.lineTo(DMAP(0, 0.));
            scene()->addPath(path, QPen(), brush)->setZValue(Curve.zOrder);
         }
      }
   }
#undef DMAP

   m_UpdateNeeded = false;
}

void FCurveView::clearCurves()
{
   m_HvLines.clear();
   m_Curves.clear();
   RebuildScene();
}

void FCurveView::addHline(float fPos, QPen const &Pen, uint32_t Flags)
{
   m_HvLines.push_back(FHvLineData(fPos, HVCURVE_Hline | Flags, Pen));
}

void FCurveView::addVline(float fPos, QPen const &Pen, uint32_t Flags)
{
//    std::cout << "\n!vline at " << fPos << "\n" << std::endl;
   m_HvLines.push_back(FHvLineData(fPos, HVCURVE_Vline | Flags, Pen));
}

void FCurveView::getCurve(TArray<float> &Data, QString &Title, uint32_t &Flags, QColor &Color, unsigned iCurveId)
{
   if (iCurveId >= m_Curves.size())
      throw std::runtime_error("requested data of non-existent curve.");
   FCurveData
      &Curve = m_Curves[iCurveId];
   Data = Curve.Data;
   Flags = Curve.Flags;
   Color = Curve.Stroke.color();
   Title = Curve.Title;
}



void FCurveView::setCurve(unsigned iCurveId, TArray<float> const &Data, QString Title, uint32_t Flags, int zOrder, QPen const &Stroke, QBrush const &Fill)
{
   if (m_Curves.size() < iCurveId+1)
      m_Curves.resize(iCurveId+1);
   FCurveData
      &Curve = m_Curves[iCurveId];
   Curve = FCurveData(Data, Flags, zOrder, Stroke, Fill);
   Curve.Title = Title;
   RebuildScene();



//    {
//       QPainterPath path;
//       float dy = -height()/(Max-Min);
//       float y0 = height() - dy * Min;
//       float dx = float(width())/nPts;
// #define F(x) (y0+(x)*dy)
//       FCurveData
//          &Curve = m_Curves[iCurveId];
//
//       path.moveTo(0, F(Min));
//       for (unsigned i = 0; i < Curve.Data.size(); ++ i)
//          path.lineTo(dx*i, F(Curve.Data[i]));
//       path.lineTo(dx*(Curve.Data.size()-1), F(Min));
//       scene()->addPath(path, pen, brush);
// #undef F
//    }
}


void FCurveView::resizeEvent(QResizeEvent */*event*/)
{
//    scene()->setSceneRect(0, 0, width(), height());
//
//    QRectF cb = getCurveBounds();
//    fitInView(QRectF(qreal(cb.left()), 0., qreal(cb.width()),1.), Qt::IgnoreAspectRatio);
//    fitInView(QRectF(qreal(cb.left()), cb.top(), qreal(cb.width()), cb.height()), Qt::IgnoreAspectRatio);
//    std::cout << fmt::format("l: {}  t: {}  w: {}  h: {}", cb.left(), cb.top(), cb.width(), cb.height()) << std::endl;
//    fitInView(QRectF(cb.left(),0.,50.,1.), Qt::IgnoreAspectRatio);
   RebuildScene();
}

// FCurveData::FCurveData(TArray<float> const &Data_, uint32_t Flags_, uint32_t Stroke_, uint32_t Fill_)
//    : Data(Data_), Flags(Flags_), Stroke(Stroke_), Fill(Fill_)
// {}

FCurveData::FCurveData(TArray<float> const &Data_, uint32_t Flags_, int zOrder_, QPen Stroke_, QBrush Fill_)
   : Data(Data_), Flags(Flags_), Stroke(Stroke_), Fill(Fill_), zOrder(zOrder_)
{}


float FCurveData::GetMin() {
   if (Data.empty())
      return 0.;
   float Min = Data[0];
   for (size_t i = 0; i < Data.size(); ++ i)
      Min = std::min(Data[i], Min);
   return Min;
}

float FCurveData::GetMax() {
   if (Data.empty())
      return 0.;
   float Max = Data[0];
   for (size_t i = 0; i < Data.size(); ++ i)
      Max = std::max(Data[i], Max);
   return Max;
}


void FCurveView::keyPressEvent(QKeyEvent */*event*/)
{
}

void FCurveView::wheelEvent(QWheelEvent *event)
{
   scaleView(pow((double)2, event->delta() / 240.0));
}


void FCurveView::mousePressEvent(QMouseEvent *event)
{
   LastX = event->x();
   LastY = event->y();
   FBase::mousePressEvent(event);
}

void FCurveView::mouseReleaseEvent(QMouseEvent *event)
{
   FBase::mouseReleaseEvent(event);
//    translate(-1000, 0);
//    std::cout << "mouse release event!" << std::endl;
}

void FCurveView::mouseMoveEvent(QMouseEvent *event)
{
   FBase::mouseMoveEvent(event);
//    if (event->buttons() != 0) {
//       float fScale = 1.;
//       float fDeltaX = fScale * (event->x() - LastX);
//       float fDeltaY = fScale * (event->y() - LastY);
//
// //       if ( event->buttons() == (Qt::LeftButton | Qt::RightButton) ) {
//          translate(fDeltaX, fDeltaY);
//          update();
// //       }
//    }
   LastX = event->x();
   LastY = event->y();
}


void FCurveView::drawBackground(QPainter *painter, const QRectF &rect)
{
//    FBase::drawBackground(painter, rect);
//    Q_UNUSED(rect);

//    QRectF
//       sceneRect = this->sceneRect();

//    QLinearGradient
//       gradient(sceneRect.topLeft(), sceneRect.bottomRight());
   QLinearGradient
      gradient(rect.topLeft(), rect.bottomRight());
   gradient.setColorAt(0, QColor(0xff606060));
//    gradient.setColorAt(1, Qt::lightGray);
   gradient.setColorAt(1, QColor(0xff202020));
   painter->fillRect(rect, gradient);
//    painter->fillRect(rect.intersect(sceneRect), gradient);
//    painter->setBrush(Qt::NoBrush);
//    painter->drawRect(sceneRect);
   // ^- yes... it is stupid. I could not resist.
}

void FCurveView::scaleView(qreal scaleFactor)
{
   qreal
      factor = transform().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
   if (factor < 0.07 || factor > 100)
      return;

   scale(scaleFactor, scaleFactor);
   // ^- transformation anchor is under mouse... pretty cool, imo.
}
