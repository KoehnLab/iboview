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

#ifndef _3D_VIEW_WIDGET_H
#define _3D_VIEW_WIDGET_H

#include "Iv.h"
#include <QAbstractTableModel>
#include <QGLWidget>
#include <QTimer>
#include <IvDataOptions.h>
// #include "IvScript.h"

class FDocument;
class FViewImpl;


// class FView3d : public QGLWidget, public IView3d {
// class FView3d : public IView3d {
class FView3d : public QGLWidget {

   Q_OBJECT // must include this if you use Qt signals/slots

public:
#include "prop_FView3d.h.inl"

public:
   typedef QGLWidget
      FBase;
   FView3d(QWidget *parent, FDocument *document);
   ~FView3d();

public:
   // return a script setting up the current view.
   QString GetViewDesc();

   bool RenderAnyLabels() const;
public slots:
   void updateData(const QModelIndex &topleft,const QModelIndex &bottomright);

   void UpdateBondMesh();
//    void executeUpdate();
   void ResetFrameBuffers();

   void PrepareDataLayoutChange();
   void FinishDataLayoutChange();
   void RehashShaders();
protected:
   void initializeGL();
   void resizeGL(int w, int h);
   void paintGL();
   void mousePressEvent(QMouseEvent *event);
   void mouseReleaseEvent(QMouseEvent *event);
   void mouseMoveEvent(QMouseEvent *event);
   void wheelEvent(QWheelEvent *event);
   void keyPressEvent(QKeyEvent *event);

//    QTimer m_UpdateTimer;
public:
   FDocument *d;
   FViewImpl *v;
   friend class FViewImpl;

   QImage grabFrameBuffer1(bool withAlpha);
   friend class IView3d;
};



// script interface
class IView3d : public FView3d {
   Q_OBJECT
public slots:
   virtual void set_camera_pos(float x, float y, float z); // = 0;
   virtual void set_camera_dir(float x, float y, float z); // = 0;
   virtual void set_camera_vup(float x, float y, float z); // = 0;
   virtual void set_camera_zoom(float z); // = 0;
   virtual void save_png(QString const &FileName); // = 0;
   virtual void set_option(QString const &OptionName, QVariant f); // = 0;
//    virtual void set_option(QString const &OptionName, double f); // = 0;
//    virtual double get_option(QString const &OptionName); // = 0;
   virtual void modify(QString const &How, float f); // = 0;

//    virtual void set_option(QString const &OptionName, bool f); // default impl calls set_view_optioni with +1 or -1.
//    virtual void set_option(QString const &OptionName, int i); // default impl calls set_view_optionf(double(i)).

   virtual void set_size(int width, int height);
public: // here for technical reasons. Not part of script interface.
   IView3d(QWidget *parent_, FDocument *document_);
   ~IView3d(); // does nothing---just to fix the vtable.
};

// script interface
class ICamera : public QObject
{
   Q_OBJECT
public slots:
   virtual void set_pos(float x, float y, float z); // = 0;
   virtual void set_dir(float x, float y, float z); // = 0;
   virtual void set_vup(float x, float y, float z); // = 0;
   virtual void set_zoom(float z); // = 0;
public:
   explicit ICamera(QObject *pView3d);
private:
   IView3d *pView3d();
};



#endif  /* _3D_VIEW_WIDGET_H */
