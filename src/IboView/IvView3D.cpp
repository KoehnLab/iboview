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

// /usr/bin/qmake-qt4 CONFIG+=debug -o Makefile main.pro && make && ./wheee
// #include <boost/math/special_functions/erf.hpp>
#include "Iv.h"
#include "IvGl.h"

#include <map>
#include <set>
#include <iostream>
#include <cmath> // for std::trunc and erfc
#include <QtGlobal> // for QT_VERSION
#include <QRect>
#include <QImage>
#include <QApplication>
#include <QClipboard>
#include <QFile>
#include <QTextStream>
#include <QMenu>
#include <QApplication>
// #include <QEvent> // for QUpdateLaterEvent (hmpf.. not exported anymore)

#include <QMimeData>
#include <QByteArray>
#include <QBuffer>

#include <stdio.h>
#include <stdexcept>
#include <sstream>
#include <cstdlib> // for rand
#include <QMouseEvent>
#include <QMutex>
#include <QMutexLocker>

#include "format.h"
#include "IvView3D.h"
#include "IvScript.h"
#include "IvDocument.h"
#include "IvMesh.h"
#include "IvIsoSurface.h"
#include "IvVolumeDataSet.h"
#include "IvOrbital.h"

#include "CtAtomSet.h"
#include "CxIo.h"
#include "CxAtomData.h"

#include "CxColor.h"
// #include "CxOsInt.h"
#include "CxVec3.h"
#include "CxMathAlgo.h" // for implicit ErfInv

#include "CxTiming.h"
using namespace ct;


// static const float fCameraDist = 80.f;
// static const float fCameraDist = 80.f;
// static const float fCameraDist = 80.f;
static const float fCameraDist = 20.f;
static const float fCameraDistActual = 100.f;


bool gluInvertMatrix(const float m[16], float invOut[16]);

enum FTransparencyStyle {
   TRANSSTYLE_DepthPeeling
};
FTransparencyStyle
//    TransparencyStyle = TRANSSTYLE_AverageColor;
   TransparencyStyle = TRANSSTYLE_DepthPeeling;


enum FPickingBaseIds {
   PICKID_Atom = 0x80000000,
   PICKID_DataRow = 0x90000000,
   PICKID_BondLine = 0xa0000000,
   PICKID_KeyMask = 0xff000000,
   PICKID_ValueMask = 0x00ffffff
};

static const unsigned PICKID_BondLineStride = 0x1000;


enum FRenderFlags {
   RENDER_Opaque = 0x01,
   RENDER_Transparent = 0x02
};


class FViewImpl
{
public:
//    uint32_t
//       ViewFlags; // combination of RENDERFLAG_*
   int LastX, LastY, FirstX, FirstY;
   bool SimpleClick; // used to distinguish left-click from view left+right view transform.
   FVec3f
      vCameraPos, vCameraDir, vCameraUp;
   float
      fZoomFactor;
   FDocument *d;
   FView3d *v;

   FViewImpl(FDocument *pDocument_, FView3d *pView_);

   QMutex
      m_MutexPainting;
   int
      m_UpdateLocked;

   FGlMatrixStack
      TrafoView,
      TrafoProj;
   FShaderSet
      AtomAndBondShader,
      OrbitalShader,
      WriteZShader, // trivial shader just producing color 0 and correct z.
      CombShader, // combination shader. combine accumulated image into output buffer.
      CombShaderG, // combination shader with actual vertex transformations.
      PickShader; // flat shader for rendering object IDs into the color buffer.
   FShaderSet
      Filter5H, // horizontal part of separable 5x5 filter
      Filter5V, // vertical part of separable 5x5 filter
      FakeAA;

   FFrameBufferPtr
      pMainFbo, // the actual window backbuffer
      pPickFbo, // back buffer we render object IDs into, if needed, because I don't feel like coding the ray tracing.
      pDpLayer0, // dummy buffers for depth peeling. Both share a single color surface
      pDpLayer1; // but have individual depth buffers.
   FFrameBufferPtr
      pAccFbo2; // accumulation buffers for final image in transparent surface assembly
   FFullScreenQuad
      FullScreenQuad;

   FImageFilterPtr
      pBlurH, pBlurV, pFakeAA;
   uint
      nMainFboSamples;

   FGlMeshPtr
      pLineMesh,
      pAtomMesh,
      pBondMesh,
      pPartialBondMesh,
      pSelectedAtomMesh;
   FGlTextBufferPtr
      pTextBuffer;
   ct::FTimerSet
      Timers;



   void SetSceneShaderUniforms(uint RenderFlags, FShaderSet &ShaderSet);
   void RenderObjects(uint RenderFlags, FShaderSet &ShaderSet);
   void RenderScene();
   void RenderPickBuffer(bool KeepPickBufferBound);

   uint width() const { return v->width(); }
   uint height() const { return v->height(); }
   QSize size() const { return v->size(); }

   void RenderFreeLabel(FFreeLabel *pFreeLabel, uint RenderFlags, FShaderSet &ShaderSet);
   void RenderFreeLine(FFreeLine *pFreeLine, uint RenderFlags, FShaderSet &ShaderSet);
   void RenderFreeObject(FFreeObject *pFreeObject, uint RenderFlags, FShaderSet &ShaderSet);
   void RenderFreeObjects(FFreeObjectSet &FreeObjects, uint RenderFlags, FShaderSet &ShaderSet);
   void RenderAtomSet(FGeometry &Geometry, uint RenderFlags, FShaderSet &ShaderSet);
   void RenderAtom(ct::FAtom const &At, uint AtomFlags, size_t iAt, uint RenderFlags, FGeometry &Geometry, FShaderSet &ShaderSet);
   // render half bond from iAt to jAt (fixme: make parameters less redundant...).
   void RenderHalfBond(ct::FAtom const &At0, size_t iAt, ct::FAtom const &At1, size_t jAt, FBondLine const &bl, uint BondFlags, FGeometry &Geometry, FShaderSet &ShaderSet);
   void RenderVolumeDataSet(FVolumeDataSet &Orb, FShaderSet &ShaderSet);
   void SetupSampleWeightBuffer(GLenum &IdWeightTexture, GLenum &IdWeightBuffer, uint nSamples);

   void SetupTrafos(int w = -1, int h = -1);
   void NormalizeCameraPos(int iCase=0);
   void ResetProjectionAndZoom(int w = -1, int h = -1);
   void AlignView(FVec3f Right, FVec3f Up);
   void RollCamera(float fAngle, int x=-1, int y=-1);

   void ClickPosition(int x, int y, Qt::MouseButton MouseButton, Qt::KeyboardModifiers KeyboardModifiers, QPoint globalPos);
   void SelectRect(int x0, int y0, int x1, int y1, Qt::MouseButton MouseButton, Qt::KeyboardModifiers KeyboardModifiers, QPoint globalPos);

   void UnselectAll();
   void SelectAtom(int iAt);

   // convert screen coordinates in terms of pixels (x,y,z) to world coordinates.
   FVec3f ScreenToWorld(FVec3f ScreenCoords);

   FIsoSurfaceSettings MakeIsoSurfaceSettings() const;

   float fAtomScale() const { return 0.4 * 0.01 * v->m_AtomScale; }
   float fBondScaleOuter() const { return .35*(.01 * v->m_BondScale * fAtomScale()/.5); }
   float fBondThinning() const { return v->m_BondThinning; }
//    float fPartialBondThresh() const { return 0.3333; }
//    float fPartialBondThresh() const { return 0.25; }
   float fPartialBondThresh() const { return v->m_FullBondThresh; }

   QWidget *FindTopmostParent();

   QPoint iPickBufferPos(int x, int y);

   int
      m_DevicePixelRatio;
};



FViewImpl::FViewImpl(FDocument *pDocument_, FView3d *pView_)
   : TrafoView(TRANSFORM_ModelView),
     TrafoProj(TRANSFORM_Projection)
{
   d = pDocument_;
   v = pView_;
   vCameraPos = FVec3f(0.0f, -20.0f, 0.0f);
   vCameraDir = FVec3f(0.0f, 1.0f, 0.0f);
   vCameraUp = FVec3f(0.0f, 0.0f, 1.0f);
   fZoomFactor = 1.0f;
   SimpleClick = false;

   NormalizeCameraPos();
   m_UpdateLocked = 0;
   m_DevicePixelRatio = 1;
#if QT_VERSION >= 0x050000
   // for Mac retina scaling. Doesn't appear to be needed on either windows or linux.
   m_DevicePixelRatio = int(v->devicePixelRatio());
#endif
}

FIsoSurfaceSettings FViewImpl::MakeIsoSurfaceSettings() const
{
   FIsoSurfaceSettings
      IsoSurfOpt;
//    IsoSurfOpt.fIsoValue = v->m_IsoThreshold/100.;
   double
      fIsoValue = v->m_IsoThreshold/100.;
   IsoSurfOpt.IsoValues.clear();
   IsoSurfOpt.IsoValues.push_back(FIsoThresholdEntry(fIsoValue, 0xffff0000));
   IsoSurfOpt.IsoValues.push_back(FIsoThresholdEntry(-fIsoValue, 0xff0000ff));
   IsoSurfOpt.fLinearDensity = v->m_IsoResolution;
   return IsoSurfOpt;
}

void FViewImpl::NormalizeCameraPos(int iCase)
{
   // always view object from same distance from center. Not sure if this is the right thing to do.
   float
      fDist;
   if (iCase == 1)
      fDist = fCameraDistActual;
   else
      fDist = fCameraDist;
   vCameraPos = vCameraPos * (-fDist/Dot(vCameraDir, vCameraPos));
}


static void InitGlew1()
{
   CheckGlError("Before GLEW Init");
   glewExperimental = GL_TRUE;
   glewInit();
   // eat up invalid enum errors.
   while (glGetError() != GL_NO_ERROR) {
   }
   CheckGlError("GLEW Init");
//    std::cout << "GLEW Init passed." << std::endl;
}

static QGLFormat MakeGlFormat()
{
   QGLFormat
      GlFormat(QGL::AlphaChannel | QGL::DepthBuffer | QGL::DoubleBuffer);
   if (TransparencyStyle != TRANSSTYLE_DepthPeeling)
      GlFormat.setSampleBuffers(true);
   else {
      // In the default depth peeling mode we technically need neither a depth buffer
      // nor double buffering on the main FBO, as all we ever do with it is use it as
      // target for blitting from the higher resolution accumulation FBs.
      GlFormat.setDepth(false);
      GlFormat.setDoubleBuffer(false);
   }

   QGLFormat::OpenGLVersionFlags
      GlVersionFlags = QGLFormat::openGLVersionFlags();
   if (GlVersionFlags & QGLFormat::OpenGL_Version_4_0) {
      g_GlVersion = 40;
      GlFormat.setVersion(4,0);
      GlFormat.setProfile(QGLFormat::CoreProfile);
   } else if (GlVersionFlags & QGLFormat::OpenGL_Version_3_2) {
      g_GlVersion = 32;
      GlFormat.setVersion(3,2);
      GlFormat.setProfile(QGLFormat::CoreProfile);
   } else {
//       throw std::runtime_error("OpenGL Version 3.0 not supported.");
      IvNotify(NOTIFY_Warning, "This program requires OpenGL version >=3.2. This is not supported by the current hardware/software setup."
        " The program will switch into partial compatibility mode and try to go on, but there is a chance that it will not work and visual quality will be impacted.");
      g_GlVersion = 30;
      GlFormat.setVersion(3,0);
      GlFormat.setProfile(QGLFormat::CoreProfile);
   }
//    GlFormat.setProfile(QGLFormat::CompatibilityProfile);
   // ^- core profile works on my card, but apparently not everywhere...
   //    update: on other machines, the compatibility profile does not work...
   return GlFormat;
}

FView3d::FView3d(QWidget *parent, FDocument *document)
   : FBase(MakeGlFormat(), parent)
//    QGL::StencilBuffer
//    : QGLWidget(QGLFormat(QGL::AlphaChannel), parent)
{
//    setMouseTracking(true);
#ifdef OPENGL_DEBUG
   glEnable(GL_DEBUG_OUTPUT);
   if (glGetError() != GL_NO_ERROR)
      throw std::runtime_error("Failed to enable debug output");
#endif
   // otherwise our alpha channel get's blended with the background of the form on mac...
   // does this work?
   setAttribute(Qt::WA_TranslucentBackground, false);
   setAutoFillBackground(false);
   d = document;
   InitProperties();

//    if (g_GlVersion < 40) {
//       SetFakeAntiAliasing(false); // doesn't seem to work on intel chips in linux.
//    }
// (fixed.)
   v = new FViewImpl(d, this);

   connect(this, SIGNAL(BondThinningChanged(double)), this, SLOT(UpdateBondMesh()));
   connect(this, SIGNAL(SuperSampleChanged(bool)), this, SLOT(ResetFrameBuffers()));
//    connect(&m_UpdateTimer, SIGNAL(timeout()), this, SLOT(update()));
//    m_UpdateTimer.start(20);


   connect(document, SIGNAL(dataChanged(const QModelIndex &,const QModelIndex &)), this, SLOT(updateData(const QModelIndex &,const QModelIndex &)));
   connect(document, SIGNAL(layoutAboutToBeChanged()), this, SLOT(PrepareDataLayoutChange()));
   connect(document, SIGNAL(layoutChanged()), this, SLOT(FinishDataLayoutChange()));
   connect(document, SIGNAL(SelectionChanged()), this, SLOT(update()));

   // this is what allows the view3d to accept key press events.
   // Default is Qt::NoFocus, which makes it impossible to receive keyboard events.
   setFocusPolicy(Qt::ClickFocus);
}

FView3d::~FView3d()
{
   CALL_GL( glUseProgram(0) );
   delete v;
//    delete d; // not safe if init of d fails... should use RAII smart pointer.
}



void FView3d::initializeGL() {
   InitGlew1();
#ifdef OPENGL_DEBUG
   glDebugMessageCallback(OutputGlDebugMessage, 0);
   glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, 0, GL_TRUE);
   glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);

   std::cout << "!OpenGL debug output enabled?" << std::endl;
//    glDisable(GL_TEXTURE_2D); // not allowed in core profile. See if this works.
#endif


   IvEmit(QString("Context initialized to GL %1.%2.").arg(format().majorVersion()).arg(format().minorVersion()));
   IvEmit(QString("GL says its context is version: %1").arg((char*)glGetString(GL_VERSION)));

   CALL_GL( glEnable(GL_DEPTH_TEST) );
   CALL_GL( glEnable(GL_CULL_FACE) );
   CALL_GL( glCullFace(GL_BACK) );

   v->FullScreenQuad.Init(); // must be done after glewInit

//    v->nMainFboSamples = 4; // set this to 0 to disable MSAA in accumulation stuff.
//       nMainFboSamples = 0; // set this to 0 to disable MSAA in accumulation stuff.
//    GLenum
//       GL_TEX_2D_MS = (nMainFboSamples==0)?  GL_TEXTURE_2D : GL_TEXTURE_2D_MULTISAMPLE;


//    GLint bla;
//    glGetFramebufferParameteriv(GL_FRAMEBUFFER, GL_FRAMEBUFFER_DEFAULT_SAMPLES, &bla);
//    nMainFboSamples = (uint) bla;
//    std::cout << fmt::format("OpenGL says my main FB has {} multi-sample samples.", nMainFboSamples) << std::endl;
   // ^- says 48. surely not!!
   // this->format() returns a QGLFormat object.
//    if (nMainFboSamples != 0)
   if (this->format().sampleBuffers())
      v->nMainFboSamples = (uint)this->format().samples();
   else
      v->nMainFboSamples = 0;
//    std::cout << fmt::format("QGLWidget says my main FB has {} multi-sample samples.", v->nMainFboSamples) << std::endl;

   RehashShaders();

   // make meshes for atoms and bonds. These will be displaced and colored (instanced)
   // to draw the actual main geometry.
   {
      // ...hm... the mesh looks okay and the tesselation degree of the level 3
      // sphere is quite similar to the 32/32 gluSphere (which has 2048 triangles
      // compared with the 1280 triangles of the subdiv sphere). But somehow it
      // does not look quite as spherical as I would like. Some triangles are
      // rather distorted, and the edges can be seen in some cases (despite the
      // exact normals...). Maybe I should subdivide something else than a
      // isocahedron? Or simply use the lateral/longitudinal grid subdivision?
      // And why does it have 642 vertices? Shouldn't it have 640?
      FBaseVertex
         vRefVertex;
      vRefVertex.dwColor = 0xffffffff;
      vRefVertex.iRole = 0;
      v->pAtomMesh = new FGlMesh(*MakeSubdivSphere(vec3f(0.,0.,0.), 1.0, 3, vRefVertex));
      v->pLineMesh = new FGlMesh(*MakeCylinder(1., 1., 1.0, 16, vRefVertex)); // <-- TODO: maybe add option to add cylinder caps?
      v->pBondMesh = new FGlMesh(*MakeCylinder(1., v->fBondThinning(), 1.0, 16, vRefVertex));
      v->pPartialBondMesh = new FGlMesh(*MakeCylinder(1., 1., 1.0, 16, vRefVertex));
      v->pSelectedAtomMesh = new FGlMesh(*MakeIcosahedron(1.8, vRefVertex));
   }

   {
      size_t nMaxCharacters = 65536;
      v->pTextBuffer = new FGlTextBuffer(nMaxCharacters);
   }
}

void FView3d::RehashShaders()
{
//    printf("! INVOKING ARB SHADER STUFF...\n");
   v->AtomAndBondShader.Create("vertex5.glsl", "pixel5.glsl", m_ShaderPath);
   if (TransparencyStyle == TRANSSTYLE_DepthPeeling) {
      v->OrbitalShader.Create("vertex5.glsl", "pixel5_orb_dp.glsl", m_ShaderPath);
      v->CombShader.Create("vertex5_pass.glsl", "pixel5_combine_dp.glsl", m_ShaderPath);
//       v->CombShaderG.Create("vertex5.glsl", "pixel5_combine_dp2.glsl", m_ShaderPath);
   }
   v->Filter5H.Create("vertex5_pass.glsl", "pixel5_sfilter_h.glsl", m_ShaderPath);
   v->Filter5V.Create("vertex5_pass.glsl", "pixel5_sfilter_v.glsl", m_ShaderPath);
//    v->pBlurH = new FImageFilter(v->Filter5H, FImageFilter::FILTER_Blur, m_ShaderPath);
//    v->pBlurV = new FImageFilter(v->Filter5V, FImageFilter::FILTER_Blur, m_ShaderPath);
//    v->FakeAA.Create("vertex5_pass.glsl", "pixel5_fakeAA.glsl", m_ShaderPath);
   v->FakeAA.Create("vertex5_pass.glsl", "pixel5_fakeAA2.glsl", m_ShaderPath);
   v->pFakeAA = new FImageFilter(v->FakeAA, FImageFilter::FILTER_Other);

   v->PickShader.Create("vertex5.glsl", "pixel5_pick.glsl", m_ShaderPath);

//    if (1) {
//       GLint
//          IdFilterStep = v->Filter5H.GetUniformLocation("FilterStep"),
//          IdFilterWeight = v->Filter5H.GetUniformLocation("FilterWeight");
// //       std::cout << fmt::format("Filter uniforms: {} {}", IdFilterStep, IdFilterWeight) << std::endl;
//    }
}


void FView3d::UpdateBondMesh()
{
   if (v->pBondMesh.get()) {
      // ^- if this isn't there most likely GL is not yet initialized.
      FBaseVertex
         vRefVertex;
      vRefVertex.dwColor = 0xffffffff;
      vRefVertex.iRole = 0;
      v->pBondMesh = new FGlMesh(*MakeCylinder(1., v->fBondThinning(), 1.0, 16, vRefVertex));
   }
}


static FMat4f MakeOrthoProjectionMatrix(float left_, float right_, float bottom_, float top_, float near_, float far_)
{
   // does the same as glOrtho.
   // see also: http://en.wikipedia.org/wiki/Orthographic_projection (it is translation into center of viewport followed by scaling)
   //           I am a bit confused by some of the 2s, but since this is what glOrtho does and it works, I guess it is correct.
   FMat4f r = vmath::identity4<float>();
   r(0,0) = 2./(right_ - left_);
   r(1,1) = 2./(top_ - bottom_);
   r(2,2) = -2./(far_ - near_);
   r(0,3) = -(right_ + left_)/(right_ - left_);
   r(1,3) = -(top_ + bottom_)/(top_ - bottom_);
   r(2,3) = -(far_ + near_)/(far_ - near_);
   r(3,3) = 1.;
   return r;
   // ^- note: beware of funky '#define's for 'far' and 'near' on win32.
}

static FMat4f MakePerspectiveProjectionMatrix(float fovy_, float aspect_, float near_, float far_)
{
   // does the same as gluPerspective.
   FMat4f r = vmath::identity4<float>();
   float f = 1./std::tan((M_PI/180.)*fovy_/2);
   r(0,0) = f/aspect_;
   r(1,1) = f;
   r(2,2) = (near_ + far_)/(near_ - far_);
   r(2,3) = (2*near_ * far_)/(near_ - far_);
   r(3,2) = -1;
   r(3,3) = 0.;
   return r;
}


void FViewImpl::ResetProjectionAndZoom(int w, int h)
{
   CheckGlError("ResetProjectionAndZoom/enter");
   if (w == -1) w = v->width()*m_DevicePixelRatio;
   if (h == -1) h = v->height()*m_DevicePixelRatio;

   CALL_GL( glViewport(0, 0, w, h) );
   CheckGlError("ResetProjectionAndZoom/glViewport");

   if (1) {
      // use ortho projection
      float halfh = 8. / fZoomFactor;
      TrafoProj.Set(MakeOrthoProjectionMatrix(-(halfh*w)/h, (halfh*w)/h, -halfh, halfh, 0.02f * fCameraDistActual, 2*fCameraDistActual));
   } else {
      // use perspective projection
      TrafoProj.Set(MakePerspectiveProjectionMatrix(70.f/(fZoomFactor*fCameraDist/30.f), (float)w/(float)h, 0.02f * fCameraDistActual, 2*fCameraDistActual));
   }
}


void FView3d::ResetFrameBuffers()
{

   CheckGlError("FView3d::ResetFrameBuffers (enter)");
   GLenum
//    DepthFormat = GL_DEPTH_COMPONENT;
      DepthFormat = GL_DEPTH_COMPONENT32F;
   // ^- doesn't work without full format... needs to be accurate enough
   //    to fully represent gl_FragCoord.z.

   // make new frame buffers. in theory this should delete the old ones automagically.
   // (unless I did something wrong)
   if (1) {
      v->pMainFbo = new FFrameBuffer("MainFbo", size()*v->m_DevicePixelRatio, true);
#if 1
      QSize
         FboSize = size()*v->m_DevicePixelRatio;
      if (TransparencyStyle == TRANSSTYLE_DepthPeeling) {
         // it cannot be done with proper MSAA---at least not without processing
         // every sample individually (which would save nothing). Do old-school
         // AA by rendering everything with double resolution explicitly, and
         // then blit it down.
         if (m_SuperSample)
            FboSize = QSize(2*width(), 2*height());
      }

      if (TransparencyStyle == TRANSSTYLE_DepthPeeling) {
         // we make two frame buffer objects. Both refer to the same color
         // surface, but they have individual depth buffers.
         uint
//             nPeelingSamples = nMainFboSamples;
            nPeelingSamples = 0;

         v->pDpLayer0 = new FFrameBuffer("DpFbo0", FboSize);
         v->pDpLayer0->Attach(new FFrameBufferAttachment("Dp0/Color", FboSize, GL_COLOR_ATTACHMENT0, nPeelingSamples, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE));
         v->pDpLayer0->Attach(new FFrameBufferAttachment("Dp0/Depth", FboSize, GL_DEPTH_ATTACHMENT, nPeelingSamples, DepthFormat, GL_DEPTH_COMPONENT, GL_FLOAT));
         v->pDpLayer0->Finalize();

         v->pDpLayer1 = new FFrameBuffer("DpFbo1", FboSize);
         v->pDpLayer1->Attach(v->pDpLayer0->Attachments[0]);
//          v->pDpLayer1->Attach(new FFrameBufferAttachment("Dp1/Color", FboSize, GL_COLOR_ATTACHMENT0, nPeelingSamples, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE));
         v->pDpLayer1->Attach(new FFrameBufferAttachment("Dp1/Depth", FboSize, GL_DEPTH_ATTACHMENT, nPeelingSamples, DepthFormat, GL_DEPTH_COMPONENT, GL_FLOAT));
         v->pDpLayer1->Finalize();

         // color only.
         v->pAccFbo2 = new FFrameBuffer("DpFbo1", FboSize);
         v->pAccFbo2->Attach(new FFrameBufferAttachment("Dp1/Color", FboSize, GL_COLOR_ATTACHMENT0, nPeelingSamples, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE));
         v->pAccFbo2->Finalize();
      }
      v->pMainFbo->Bind();
      CheckGlError("frame buffer setup");
#endif
   }
   // make a smaller resolution frame buffer which we may use to render object IDs for picking.
   if (1) {
      QSize
         PickBufferSize = size()/2;
      v->pPickFbo = new FFrameBuffer("PickFbo", PickBufferSize);
      v->pPickFbo->Attach(new FFrameBufferAttachment("Pick/Color", PickBufferSize, GL_COLOR_ATTACHMENT0, 0, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE));
      v->pPickFbo->Attach(new FFrameBufferAttachment("Pick/Depth", PickBufferSize, GL_DEPTH_ATTACHMENT, 0, DepthFormat, GL_DEPTH_COMPONENT, GL_FLOAT));
      v->pPickFbo->Finalize();
      CheckGlError("pick buffer setup");
   }
   v->pMainFbo->Bind();
   CheckGlError("FView3d::ResetFrameBuffers (leave)");

   // FIXME: technically we would need to setup the projection and zoom again now,
   // because the projection matrix depends on w/h. However, this is actually done
   // in RenderScene automatically, so rendering is not affected. It does, however,
   // affect the inverse picking transformation (ScreenToWorld) which is used in some
   // places.
}


void FView3d::resizeGL(int w, int h) {
   CheckGlError("FView3d::resizeGL (enter)");
   v->ResetProjectionAndZoom(w,h);
   ResetFrameBuffers();

   CheckGlError("FView3d::resizeGL (leave)");
   // may cause a crash due to current re-entrant behavior...
   IvNotify(NOTIFY_Information, IvFmt("Viewport: %1 x %2", w, h));
}

// that's the rasmol CPKnew colors, according to jmol homepage.
// See MakeColorCodes.py
// uint32_t ElementColors[109] = {0xffffff,0xffc0cb,0xb22121,0xff1493,0x00ff00,0xd3d3d3,0x87cee6,0xff0000,0xdaa520,0xff1493,0x0000ff,0x228b22,0x696969,0xdaa520,0xffaa00,0xffff00,0x00ff00,0xff1493,0xff1493,0x696969,0xff1493,0x696969,0xff1493,0x696969,0x696969,0xffaa00,0xff1493,0x802828,0x802828,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0x696969,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xa020f0,0xff1493,0xff1493,0xffaa00,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xdaa520,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493};

// as above, but replaced carbon by 0x999999 (0.6/0.6/0.6)
// static uint32_t ElementColors[109] = {0xffffff,0xffc0cb,0xb22121,0xff1493,0x00ff00,0x999999,0x87cee6,0xff0000,0xdaa520,0xff1493,0x0000ff,0x228b22,0x696969,0xdaa520,0xffaa00,0xffff00,0x00ff00,0xff1493,0xff1493,0x696969,0xff1493,0x696969,0xff1493,0x696969,0x696969,0xffaa00,0xff1493,0x802828,0x802828,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0x696969,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xa020f0,0xff1493,0xff1493,0xffaa00,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xdaa520,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493};
//
// FColor GetAtomColor(FAtom const &At) {
//    FVec3f
//       v;
//    float &r = v[0], &g = v[1], &b = v[2];
//    switch ( At.iElement ) {
//       case 1: r = 1.0; g = 1.0; b = 1.0; break; // H
//       case 6: r = 0.6; g = 0.6; b = 0.6; break; // C
//       default:
//          uint32_t dw;
//          dw = ElementColors[At.iElement-1];
//          r = (1.f/255.f)*((dw>>16)&0xff);
//          g = (1.f/255.f)*((dw>>8)&0xff);
//          b = (1.f/255.f)*((dw>>0)&0xff);
//    }
//    return v;
// }
//
// void SetGlAtomColor1(GLint IdColor, FAtom const &At)
// {
//    if (IdColor == -1)
//       return;
//    FColor rgb = GetAtomColor(At);
//    CALL_GL( glUniform4f(IdColor, rgb[0], rgb[1], rgb[2], 1.0f) ); // rgba
// }


// static uint32_t ElementColors[109] = {0xffffff,0xffc0cb,0xb22121,0xff1493,0x00ff00,0x999999,0x87cee6,0xff0000,0xdaa520,0xff1493,0x0000ff,0x228b22,0x696969,0xdaa520,0xffaa00,0xffff00,0x00ff00,0xff1493,0xff1493,0x696969,0xff1493,0x696969,0xff1493,0x696969,0x696969,0xffaa00,0xff1493,0x802828,0x802828,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0x696969,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xa020f0,0xff1493,0xff1493,0xffaa00,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xdaa520,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493};

FColor MakeGlColor(uint32_t dw) {
   FVec3f
      v;
   float &r = v[0], &g = v[1], &b = v[2];
   r = (1.f/255.f)*((dw>>16)&0xff);
   g = (1.f/255.f)*((dw>>8)&0xff);
   b = (1.f/255.f)*((dw>>0)&0xff);
   return v;
}

void SetGlAtomColor1(GLint IdColor, uint32_t dwColor)
{
   if (IdColor == -1)
      return;
   FColor rgb = MakeGlColor(dwColor);
   CALL_GL( glUniform4f(IdColor, rgb[0], rgb[1], rgb[2], 1.0f) ); // rgba
}


void SetGlPickId(GLint IdObjectId, uint32_t ObjectId)
{
   if (IdObjectId == -1)
      return;
   float
      a = (1.f/255.f)*((ObjectId>>24)&0xff),
      b = (1.f/255.f)*((ObjectId>>16)&0xff),
      g = (1.f/255.f)*((ObjectId>>8)&0xff),
      r = (1.f/255.f)*((ObjectId>>0)&0xff);
   CALL_GL( glUniform4f(IdObjectId, r, g, b, a) );
   // ^- This works for me, but I wonder if it is save...  there is a
   //    glUniform4ui, but my understanding of the docs is that I cannot use it
   //    to set floating point buffer data. And apparently there is lots of
   //    incompatibility with uint vectors on various GL versions.
   //
   // Note also that the actual format of the GL buffer is RGBA, apparently. Not sure about
   // endian-ness etc. R and B are definitely exchanged with QT's color defs.
}



FUniformProxy::~FUniformProxy()
{}

struct FSetSceneUniformProxy : public FUniformProxy
{
   void SetUniforms(FShaderSet &ShaderSet); // override.

   FSetSceneUniformProxy(FView3d *v_, uint RenderFlags_) : v(v_), RenderFlags(RenderFlags_) {}
   ~FSetSceneUniformProxy() {};
protected:
   FView3d *v;
   uint RenderFlags;
};

// void FViewImpl::SetSceneShaderUniforms(uint RenderFlags, FShaderSet &ShaderSet)
void FSetSceneUniformProxy::SetUniforms(FShaderSet &ShaderSet)
{
   // feed some parameters controlling the appearance into the pixel shader.
   char const
      *pShaderRegNames[4] = {"ShaderReg0", "ShaderReg1", "ShaderReg2", "ShaderReg3"};
   double
      fShaderReg[4] = {0., 0., 0., 0.};
   if ((RenderFlags & RENDER_Opaque) != 0) {
      // atoms and bonds.
      fShaderReg[0] = v->m_ShaderRegA0; fShaderReg[1] = v->m_ShaderRegA1;
      fShaderReg[2] = v->m_ShaderRegA2; fShaderReg[3] = v->m_ShaderRegA3;
   } else {
      // orbitals and selection markers.
      fShaderReg[0] = v->m_ShaderRegO0; fShaderReg[1] = v->m_ShaderRegO1;
      fShaderReg[2] = v->m_ShaderRegO2; fShaderReg[3] = v->m_ShaderRegO3;
   }
   for (size_t i = 0; i != 4; ++ i)
      CALL_GL( glUniform1f(ShaderSet.GetUniformLocation(pShaderRegNames[i]), fShaderReg[i]) );
   if (v->m_FadeType == 0) {
      CALL_GL( glUniform1f(ShaderSet.GetUniformLocation("FadeWidth"), 0.) );
      CALL_GL( glUniform1f(ShaderSet.GetUniformLocation("FadeBias"), 0.) );
   } else {
      CALL_GL( glUniform1f(ShaderSet.GetUniformLocation("FadeWidth"), v->m_FadeWidth) );
      CALL_GL( glUniform1f(ShaderSet.GetUniformLocation("FadeBias"), v->m_FadeBias) );
   }
}

void FViewImpl::SetSceneShaderUniforms(uint RenderFlags, FShaderSet &ShaderSet)
{
   FSetSceneUniformProxy
      proxy(v, RenderFlags);
   proxy.SetUniforms(ShaderSet);
}



void FViewImpl::RenderObjects(uint RenderFlags, FShaderSet &ShaderSet)
{
   CALL_GL( glUseProgram(ShaderSet.hProgram) );
   SetSceneShaderUniforms(RenderFlags, ShaderSet);

   bool
      Wireframe = false;
   if (Wireframe) {
      CALL_GL( glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) );
   }

   // export view and projection transformations to their respective places
   // in the shader's uniform records.
   TrafoProj.Actualize(ShaderSet);
   TrafoView.Actualize(ShaderSet);

//    FDataSetList
//       *pDataList = v->d->GetCurrentFrameData();
   FDataSetList
      DataList;
   FDataSetList
      *pDataList = &DataList;
   if (bool(v->d->GetCurrentFrame())) {
      // may differ from GetCurrentFrameData() in that it may contain additional
      // `internal` data sets, such as free-object sets if not handled as part
      // of the common data.
      DataList = v->d->GetCurrentFrame()->AllData();
   }
   if (pDataList) {
      for (size_t iDataSet = 0; iDataSet < pDataList->size(); ++ iDataSet)
      {
         FDataSetPtr
            pActiveData = (*pDataList)[iDataSet];
         if (!pActiveData->Active)
            continue; // current set is hidden.

         FGeometry
            *pGeometryData = dynamic_cast<FGeometry*>(pActiveData.get());
//          if (pGeometryData != 0 && (RenderFlags & RENDER_Opaque) != 0) {
         if (pGeometryData != 0) {
            RenderAtomSet(*pGeometryData, RenderFlags, ShaderSet);
         }

         if (!pGeometryData)
            // we don't want to make IDs for "the entirety of the geometry" here.
            // atoms get their own Ids.
            SetGlPickId(ShaderSet.hObjectId, PICKID_DataRow + iDataSet);

         // maybe redesign this for double dispatch? would need #including
         // IvView3d in the data set stuff then, however. also quite ugly.
         FFreeObjectSet
            *pFreeObjects = dynamic_cast<FFreeObjectSet*>(pActiveData.get());
//          if (pFreeObjects != 0 && (RenderFlags & RENDER_Transparent) != 0)
//          if (pFreeObjects != 0 && (RenderFlags & RENDER_Opaque) != 0)
//             RenderFreeObjects(*pFreeObjects, RenderFlags, ShaderSet);
         if (pFreeObjects != 0 && (RenderFlags & RENDER_Opaque) != 0)
            RenderFreeObjects(*pFreeObjects, RenderFlags, ShaderSet);

         FVolumeDataSet
            *pVolumeDataSet = dynamic_cast<FVolumeDataSet*>(pActiveData.get());
         if (pVolumeDataSet != 0 && (RenderFlags & RENDER_Transparent) != 0) {
            RenderVolumeDataSet(*pVolumeDataSet, ShaderSet);
         }
      }
   }

   if (RenderFlags == RENDER_Opaque) {
      // render text (atom labels).
      if (v->RenderAnyLabels()) {
         CALL_GL( glEnable(GL_BLEND) );
         CALL_GL( glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) );
         FSetSceneUniformProxy
            proxy(v, RenderFlags);
         pTextBuffer->DrawText1(&TrafoView, &TrafoProj, &proxy);
         CALL_GL( glDisable(GL_BLEND) );
      }
   }


//    if (0) {
//       // for checking picking transformation (inverse view trafo).
//       FVec3f vx = ScreenToWorld(FVec3f(LastX, LastY, 0.));
//       vx -= Dot(vCameraDir, vx)*vCameraDir;
//       FAtom dummy(ct::FVector3(vx[0],vx[1],vx[2]), "He");
//       TrafoView.Push();
//       RenderAtom(dummy, 0, 0, 0, ShaderSet);
//       TrafoView.Pop();
//    }

   if (Wireframe) {
      CALL_GL( glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) );
   }

   CALL_GL( glUseProgram(0) );
}


void FViewImpl::RenderPickBuffer(bool KeepPickBufferBound)
{
   pPickFbo->Bind();
   CALL_GL( glClearDepth(1.) );
   CALL_GL( glClearColor(0., 0., 0., 0.) );
   CALL_GL( glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) );

   // disable all sampling and fancy stuff? we just want to render flat colors.

   SetupTrafos(pPickFbo->Size.width(), pPickFbo->Size.height());

   RenderObjects(RENDER_Opaque | RENDER_Transparent, PickShader);

   if (!KeepPickBufferBound)
      pMainFbo->Bind();
}


void FViewImpl::SetupTrafos(int w, int h)
{
   // compute projection matrix and setup viewport
   ResetProjectionAndZoom(w, h);

   NormalizeCameraPos(); // FIXME: probably should not be here...

   {
      // make a view matrix representing the current camera settings.
      FVec3f
         vUp = vCameraUp,
         vDir = vCameraDir,
         vRight = Cross(vDir, vUp),
         vPos = vCameraPos;
      vPos = vPos - Dot(vPos,vDir)*vDir - fCameraDistActual*vDir;
      mat4f
         mViewT = hStack(vRight, vUp, -vDir, FVec3f(0.,0.,0.)),
         mView = vmath::transpose(mViewT);
         // ^- we need the inverse (=transpose) of mViewT as the view matrix.
         //    stuff is moved *from* the directions we give *to* x/y/z, after
         //    all.

      TrafoView.Clear();
      TrafoView.Set(mView);
      TrafoView.Translate(-vPos);
   }
}


void FViewImpl::RenderScene()
{
   CheckGlError("FViewImpl::RenderScene (enter)");
   bool
      BindSucceeded = pMainFbo->Bind(GL_DRAW_FRAMEBUFFER); // should not be required here!
   CALL_GL( glBindFramebuffer(GL_READ_FRAMEBUFFER, 0) );  // should not be required here!
   if (!BindSucceeded) {
      // this one is for MacOS... where apparently one cannot actually
      // initialize any frame buffers before the window is first physically shown
      const_cast<QGLContext*>(this->v->context())->makeCurrent();
      v->ResetFrameBuffers();
      if (!pMainFbo->Bind(GL_DRAW_FRAMEBUFFER))
         // still doesn't work? Probably can't fix it now.
         return;
   }

   GLfloat fClearAlpha = 0.0; // that's what it SHOULD be.
   if (g_WorkAroundAlphaCompositing || !v->m_SaveAlpha)
      fClearAlpha = 1.0; // see comments in IvGl.h

//    CALL_GL( glClearColor(1., 1., 1., 0.) );
   CALL_GL( glClearColor(1., 1., 1.000001, fClearAlpha) );
   // ^- glClearColor(1.,1.,1.,.0.) does not work on intel with mesa 9.2.3.
   // Any other color does, including glClearColor(1., 1., 1.000001, 0.) (which
   // resolves to the very same 32bit RGBA value).
   // I wish I was making this up.

//    CALL_GL( glClearColor(0.015, 0.027, 0.09, fClearAlpha) );

//    glClearColor(0., 0., 0., 0.);
//    glClearColor(0.8, 0.8, 0.8, 0);

   CheckGlError("RenderScene()","Clear Buffers");
   pTextBuffer->Clear();

   SetupTrafos();

   // render non-transparent things
   if (TransparencyStyle != TRANSSTYLE_DepthPeeling) {
      CALL_GL( glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) );
      RenderObjects(RENDER_Opaque, AtomAndBondShader);
      CheckGlError("RenderOpaque");
   }


   if (TransparencyStyle == TRANSSTYLE_DepthPeeling) {
      if (1) {
         pDpLayer0->Bind(GL_DRAW_FRAMEBUFFER);
         CALL_GL( glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT) );
         RenderObjects(RENDER_Opaque, AtomAndBondShader);
         pAccFbo2->BlitFrom(&*pDpLayer0, GL_COLOR_BUFFER_BIT, GL_NEAREST);
      }

      // note: at this moment the main depth layer is stored in pDpLayer0's depth buffer.
      CALL_GL( glClearDepth(0.) );  // note: this clears to 'everything very close.'. default is 1.0 (far away)
//       glClearColor(0., 0.2, 0.4, 0.);
      CALL_GL( glClearColor(0., 0., 0., 0.) );
      CALL_GL( glDepthFunc(GL_GEQUAL) );
      CheckGlError("DP", "setup 0");

      // render transparency layers back to front.
      for (uint iLayer = 0; iLayer < (uint)v->m_DepthPeelingLayers; ++ iLayer) {
         CALL_GL( glUseProgram(OrbitalShader.hProgram) );

         // put reference depth from pDpLayer0 to pDpLayer1
         std::swap(pDpLayer0, pDpLayer1);

         // render to (new) pAccFbo.
         pDpLayer0->Bind();

         // disable blending.
         CALL_GL( glDisable(GL_BLEND) );
         CALL_GL( glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) );
         CheckGlError("DP", "clear buffers");

         // render stuff which lies before depth in pDpLayer1,
         // but otherwise is as far back as possible.
//          std::cout << "! DP TEX 0: " <<  pDpLayer1->Attachments[1]->pDesc << std::endl;
         pDpLayer1->Attachments[1]->BindAsTexture(GL_TEXTURE0); // attachment 0 -> Depth1
         CALL_GL( glUniform1i(OrbitalShader.GetUniformLocation("Depth1"), 0) ); // set texture indices.
         CheckGlError("DP", "link depth1");

         if (v->m_RenderBacksides)
            CALL_GL( glDisable(GL_CULL_FACE) );

         RenderObjects(RENDER_Transparent, OrbitalShader);

         if (v->m_RenderBacksides)
            CALL_GL( glEnable(GL_CULL_FACE) );

         pDpLayer1->Attachments[1]->UnbindAsTexture(GL_TEXTURE0); // should not be required.

         // render color layer in pDpLayer0 into main buffer.
         // (...).
         if (0) {
            if (iLayer == 1) {
               glBindFramebuffer(GL_READ_FRAMEBUFFER, pDpLayer0->hFbo);
               glBindFramebuffer(GL_DRAW_FRAMEBUFFER, pMainFbo->hFbo);
//                glBlitFramebuffer(0, 0, width(), height(), 0, 0, width(), height(), GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);
               glBlitFramebuffer(0, 0, width(), height(), 0, 0, width(), height(), GL_COLOR_BUFFER_BIT, GL_NEAREST);
               CheckGlError("copy acc->main");
            }
         } else {
//             pMainFbo->Bind(GL_DRAW_FRAMEBUFFER);
//             pAccFbo1->Bind(GL_DRAW_FRAMEBUFFER);
            pAccFbo2->Bind(GL_DRAW_FRAMEBUFFER);

            CALL_GL( glUseProgram(CombShader.hProgram) );
            CALL_GL( glDepthMask(false) ); // disable z-write
            CALL_GL( glDisable(GL_DEPTH_TEST) );
            CALL_GL( glEnable(GL_BLEND) );
            CALL_GL( glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE) );
//             glEnable(GL_TEXTURE_2D);
            pDpLayer0->Attachments[0]->BindAsTexture(GL_TEXTURE0); // attachment 0 -> LayerColor
            CALL_GL( glUniform1i(CombShader.GetUniformLocation("LayerColor"), 0) ); // set texture indices.

            FullScreenQuad.Draw(CombShader.hProgram);
            CheckGlError("re-assemble image");

            CALL_GL( glEnable(GL_DEPTH_TEST) );
//             glDisable(GL_TEXTURE_2D);
            CALL_GL( glDisable(GL_BLEND) );
            CALL_GL( glDepthMask(true) ); // re-enable z-write
         }
//          break; // just first layer for a start.
      }

      // final image is now in pAccFbo2. Anti-alias/downsample (if asked for)
      // and move to front buffer
      if (v->m_SuperSample) {
         if (v->m_FakeAntiAliasing) {
            // fake-anti-alias the (normally supersampled)
            // accumulation frame buffer, which should now contain
            // our final image.
            pFakeAA->Apply(pDpLayer1, pAccFbo2, pAccFbo2->width(), pAccFbo2->height(), FullScreenQuad);
            // ..and blit the final image to our (single buffered, depth-less) front buffer.
            pMainFbo->BlitFrom(&*pDpLayer1, GL_COLOR_BUFFER_BIT, GL_LINEAR);
         } else {
            // blit directly.
            pMainFbo->BlitFrom(&*pAccFbo2, GL_COLOR_BUFFER_BIT, GL_LINEAR);
         }
      } else {
         if (v->m_FakeAntiAliasing) {
            pFakeAA->Apply(pMainFbo, pAccFbo2, pAccFbo2->width(), pAccFbo2->height(), FullScreenQuad);
         } else {
            pMainFbo->BlitFrom(&*pAccFbo2, GL_COLOR_BUFFER_BIT, GL_NEAREST);
            // ^- this does not work on intel. pAccFbo2 content is definitely
            // correct (drawing with full screen quad works fine), data formats
            // are consistent, all FBOs and attachments are complete, and there
            // is no GL error of any kind. There is simply nothing. The call
            // seems to be completely ignored. Since the other paths do seem
            // to work I'll just ignore this problem for now.
         }
      }

      // reset depth function and stuff.
      CALL_GL( glDepthFunc(GL_LEQUAL) );
      CALL_GL( glClearDepth(1.) );  // reset clearing depth.
   }
   pMainFbo->Bind(GL_FRAMEBUFFER);
   CALL_GL( glBindFramebuffer(GL_READ_FRAMEBUFFER, 0) );  // should not be required here! it's an attempt to get Macs to cooperate.

   CheckGlError("FViewImpl::RenderScene (leave)");
}

void FView3d::paintGL() {
   if (v->m_UpdateLocked)
      return;

//    std::cout << "!inside paintGL." << std::endl;
//    QMutexLocker
//       PaintMutexLock(&v->m_MutexPainting);
   if (v->m_MutexPainting.tryLock()) {
      // ^- this is not about threading... it is an ugly hack to prevent
      //    re-entrant calling of paintGL (paintGL may cause iso surfaces to be
      //    rendered, which might in turn produce some status updates causing
      //    QCoreApplication to process events.). This is, of course, all a horribly
      //    ugly hack and it is quite surprising that it works at all, to some degree.
      //    Must definitely be fixed later on.

      // hm... I think I need to remove additional update events from the queue...
      // otherwise we may simply get too many of them. atm the render time is typically
      // less than 2ms, but still the output can seem laggy, because the frequency of mouse
      // events is even higher.

   //    static double fLastUpdate = -1e30;
   //    double fCurrentTime = ct::Second();
   //    if (fCurrentTime < fLastUpdate + 0.03) {
   //       return;
   //    }
   //    fLastUpdate = fCurrentTime;
   //
   //    v->Timers.SetLevel(99);
   //    v->Timers.Reset();
   //    v->Timers.Resume(1, "Render");



      v->RenderScene();
      CALL_GL( glFinish() );
      // ^- interesting. That makes everything much *less* sluggish in complex scenes.
      //    maybe this allows for better event compression in QT?

   //    v->Timers.Pause(1, "Render");
   //    v->Timers.PrintReport(std::cerr, 0.);
   //    std::cout.flush();
   //    m_UpdateTimer.stop();
      v->m_MutexPainting.unlock();
   }
}


void FView3d::PrepareDataLayoutChange()
{
   if (v->m_UpdateLocked)
      return IvNotify(NOTIFY_Error, "re-entrant call of PrepareDataLayoutChange()?");
   v->m_UpdateLocked = 1;
//    IvEmit("!Setting update lock = 1");
}

void FView3d::FinishDataLayoutChange()
{
   v->m_UpdateLocked = 0;
//    IvEmit("!Setting update lock = 0");
   update();
}




void FViewImpl::RenderAtom(ct::FAtom const &At, uint AtomFlags, size_t iAt, uint RenderFlags, FGeometry &Geometry, FShaderSet &ShaderSet)
{
   if (0 != (AtomFlags & ATOM_Hidden))
      return;
   if (!v->GetShowHydrogens() && (At.iElement == 1))
      return;

   if (((RenderFlags & RENDER_Opaque) == 0) && ((AtomFlags & ATOM_Selected) == 0))
      // pure opaque pass. Atoms only have transparent contributions if selected.
      return;

   SetGlPickId(ShaderSet.hObjectId, PICKID_Atom + iAt);

   FElementOptions
      *pElementOptions = d->pElementOptions(int(iAt), &*Geometry.pAtoms);

   // hack up the view transformation in order to move the uniform
   // sphere at the origin to the current atom's location. Then render
   // the uniform sphere. (the same one for all atoms).
   // Note: if there are lots of atoms, geometry instancing and possibly
   // geometry shader subdivision might help considerably.
   TrafoView.Push();
   float
//       fRadius = fAtomScale()*GetAtomDrawRadius(At.iElement);
      fRadius = fAtomScale() * pElementOptions->GetDrawRadius();
   TrafoView.Translate(FVec3f(At.vPos));
   TrafoView.Scale(FVec3f(fRadius,fRadius,fRadius));
   TrafoView.Actualize(ShaderSet);
   if (RenderFlags & RENDER_Opaque) {
//       SetGlAtomColor1(ShaderSet.hDiffuseColor, At);
      SetGlAtomColor1(ShaderSet.hDiffuseColor, pElementOptions->GetColor());
      pAtomMesh->Draw();

      if (v->RenderAnyLabels()) {
         bool IsC = (At.iElement == 6);
         bool IsH = (At.iElement == 1);
         std::string c;
         if ((v->m_LabelElements && !(IsH || IsC)) ||
            (v->m_LabelElementsC && IsC) ||
            (v->m_LabelElementsH && IsH))
         {
//             c = At.GetElementName();
            c = q2s(pElementOptions->GetDrawName());
         }
         if (v->m_LabelAtomNumbers) {
            if (!c.empty())
               c += " ";
            c += fmt::format("{}", (iAt+1));
         }
         if (!c.empty()) {
            // FIXME: why not just leave it as wchar_t?! or at least use std::wstring...
            wchar_t wc[64];
            mbstowcs(wc, c.c_str(), 64);

            uint32_t dw;
   //          dw = 0xff000000 | irgb(pElementOptions->GetColor());
            dw = irgb(ModBrightness(FColor(pElementOptions->GetColor()), v->m_LabelBrightness).uint32()) | 0xff000000;

      //       *sqrt(fZoomFactor)
            FVec3f x = FVec3f(At.vPos) - fRadius * vCameraDir;
            float fSz = v->m_LabelSize;
            if (At.iElement == 1)
               fSz *= v->m_LabelSizeH;
            pTextBuffer->Print(vec3f(x[0],x[1],x[2]), fSz/100.f, dw, wc, FGlTextBuffer::HALIGN_Center | FGlTextBuffer::VALIGN_Center);
//             IvEmit("▒  emit-text((%1,%2,%3), '%4', %5, %6)", x[0], x[1], x[2], "??", fmt::format("0x{:08x}",dw), fSz/100.f);
         }
      }
   }
   if (bool(RenderFlags & RENDER_Transparent) && bool(AtomFlags & ATOM_Selected)) {
//       glUniform4f(ShaderSet.hDiffuseColor, 0.7, 1.0, 0.7, 0.5);
//       FVec3f rgb = .4f * FVec3f(GetAtomColor(At)) + .6f * FVec3f(.7, 1., .7);
//       FVec3f rgb = .4f * FVec3f(GetAtomColor(At)) + .6f * FVec3f(1., 1., 1.);
      FVec3f rgb = .4f * FVec3f(MakeGlColor(pElementOptions->GetColor())) + .6f * FVec3f(1., 1., 1.);
      CALL_GL( glUniform4f(ShaderSet.hDiffuseColor, rgb[0], rgb[1], rgb[2], 0.5) );
      pSelectedAtomMesh->Draw();
   }

   TrafoView.Pop();
   SetGlPickId(ShaderSet.hObjectId, 0);
}

// Make two directions orthogonal to vDir,
// such that the coordinate system (vOrth0, vOrth1, vIn) is right-handed.
void MakeOrthDirections(FVec3f &vOrth0, FVec3f &vOrth1, FVec3f const &vIn)
{
   FVec3f
      vDir = vIn,
      vRef0(1.,0.,0.),
      vRef1(0.,1.,0.);
   vDir.Normalize();
   Cross(vOrth0, vDir, vRef0);
   if ( vOrth0.LengthSq() < 1e-1 )
//       Cross(vOrth0, vRef1, vDir);
      Cross(vOrth0, vDir, vRef1);
   vOrth0.Normalize();
   Cross(vOrth1, vDir, vOrth0);
}


// void FViewImpl::RenderHalfBond(ct::FAtom const &At0, size_t iAt, ct::FAtom const &At1, size_t jAt, FBondLine const &bl, uint BondFlags, FGeometry &Geometry, FShaderSet &ShaderSet)
// {
//    (void)bl;
//    bool
//       DottedBond = (BondFlags & BOND_Partial) != 0;
//
//    FElementOptions
//       *pElementOptions0 = d->pElementOptions(int(iAt), &*Geometry.pAtoms);
//    if (0 == (BOND_Grey & BondFlags)) {
//       // set to source atom color.
// //       SetGlAtomColor1(ShaderSet.hDiffuseColor, At0);
//       SetGlAtomColor1(ShaderSet.hDiffuseColor, pElementOptions0->GetBondColor1());
//
//    } else {
//       // set bond line as gray.
//       CALL_GL( glUniform4f(ShaderSet.hDiffuseColor, 0.2, 0.2, 0.2, 1.0f) );
//    }
//
//    float fBondScaleOuter = this->fBondScaleOuter();
//
//    if (DottedBond) {
//       fBondScaleOuter *= 0.2;
//    }
//
//    float
//       r0 = pElementOptions0->GetDrawRadius();
// //       r0 = GetAtomDrawRadius(At0.iElement);
//
//    // create a random orthogonal direction to the bond line.
//    FVec3f
//       v0(At0.vPos[0], At0.vPos[1], At0.vPos[2]),
//       v1(At1.vPos[0], At1.vPos[1], At1.vPos[2]),
//       vOrth0, vOrth1,
//       vDiff = v1 - v0,
//       vDir = vDiff;
//    vDir.Normalize();
//    MakeOrthDirections(vOrth0, vOrth1, vDir);
//
//    // relative offset of the base of the cylinder from the atomic position.
//    double fOffs = 0.3 * fAtomScale();
//
//    FMat4f mView = hStack(vOrth0, vOrth1, vDir, v0 + float(fOffs * r0) * vDir);
//
//    float fHeight = 0.5*std::sqrt(vDiff.LengthSq()) - r0*fOffs;
//
//    TrafoView.Push();
//    TrafoView.Multiply(mView);
//    // ^- z axis is now along bond direction, with z=0 on base position.
//
//    if (!DottedBond) {
//       // regular bond
//       TrafoView.Scale(FVec3f(fBondScaleOuter,fBondScaleOuter,fHeight));
//       TrafoView.Actualize(ShaderSet);
//       pBondMesh->Draw();
//       // ^- note: this kind of transformation is highly non-unitary. For this reason
//       //    we do some special care in the tranformation of normals, and actually
//       //    calculate and set the normal transformations properly whenever uploading
//       //    the view transformation (admittedly, it is done in a somewhat hacky way).
//    } else {
//       // bond forming/breaking-bond.
//       int
//          nSegments = 1 + fHeight / 0.7;
//       float
//          fDotScale = 0.4, // smaller: smaller holes.
//          fSegmentHeight = fHeight / (nSegments - fDotScale/2); // with -.5 they join exactly (no space), with -0. the space is twice too large.
//       TrafoView.Scale(FVec3f(fBondScaleOuter,fBondScaleOuter,(1-fDotScale)*fSegmentHeight));
//       TrafoView.Actualize(ShaderSet);
//       for (int i = 0; i < nSegments; ++ i) {
//          pPartialBondMesh->Draw();
//          TrafoView.Translate(FVec3f(0., 0., 1./(1.-fDotScale)));
//          TrafoView.Actualize(ShaderSet);
//       }
//    }
//
//    TrafoView.Pop();
//
//    IR_SUPPRESS_UNUSED_WARNING(jAt);
// }

#if 0
static double pow3(double x) {
   return x*x*x;
}

// recycled from: http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// returns |x|^p * sign(x), with sign equal to the sign of x.
static double signed_pow(double x, double p) {
   return sign(x) * std::pow(std::abs(x), p);
}

static double fBondDotWeightTrafo(double x) {
   // this is a p'th order polynomial centered at 0.5, such that
   //   f(0.) = 0.,
   //   f(.5) = .5
   //   f(1.) = 1.
   // and the in-between values get mapped either closer to 0.5 (p > 1)
   // or closer to the boundaries 0.0/1.0 (p < 1.0).
   // See bond_scale.py.
   double p = 0.75;
   return .5*(signed_pow(2*(x-.5), p) + 1);
}

// this one puts in a 3rd attractor at 0.5.
static double fBondDotWeightTrafo3(double fBondWeight) {
   if (fBondWeight > .5)
      return .5 + .5*fBondDotWeightTrafo(2.*(fBondWeight - .5));
   else
      return .5*fBondDotWeightTrafo(2.*fBondWeight);
}
#endif

static double fBondDotWeightTrafo4(double fBondWeight) {
   return 1. - sqr(1. - fBondWeight);
}



static void DrawBond1(double fBondScaleOuter, double fHeight, double fBondWeight, FGlMatrixStack &TrafoView, FGlMeshPtr pBondMesh, FGlMeshPtr pPartialBondMesh, FShaderSet &ShaderSet) {
//    if (!DottedBond) {
   assert(fBondWeight >= 0. && fBondWeight <= 1.);
   if (fBondWeight >= 1.) {
      // regular bond
      TrafoView.Scale(FVec3f(fBondScaleOuter,fBondScaleOuter,fHeight));
      TrafoView.Actualize(ShaderSet);
      pBondMesh->Draw();
      // ^- note: this kind of transformation is highly non-unitary. For this reason
      //    we do some special care in the tranformation of normals, and actually
      //    calculate and set the normal transformations properly whenever uploading
      //    the view transformation (admittedly, it is done in a somewhat hacky way).
   } else {
      // bond forming/breaking-bond.
      double
         fDotWeight;
      if (0) {
         fDotWeight = fBondWeight;
      } else {
//          fDotWeight = fBondDotWeightTrafo(fBondWeight);
//          fDotWeight = fBondDotWeightTrafo3(fBondWeight);
         fDotWeight = fBondDotWeightTrafo4(fBondWeight);
//          fBondScaleOuter *= std::sqrt(fBondWeight/fDotWeight);
         fBondScaleOuter *= std::sqrt(fBondWeight);
         // note: ^- sqrt() here because visual weight is determined by the volume, which is
         //       proportional to the radius squared.
      }

      int
         nSegments = std::floor(1 + fHeight / 0.5);
//          nSegments = std::floor(.5 + 1 + fHeight / 0.7);
      float
         fDotScale = 1. - fDotWeight, // <- note: used to default to 0.4
         fSegmentHeight = fHeight / (nSegments - fDotScale/2); // with -.5 they join exactly (no space), with -0. the space is twice too large.
      TrafoView.Scale(FVec3f(fBondScaleOuter,fBondScaleOuter,(1-fDotScale)*fSegmentHeight));
      TrafoView.Actualize(ShaderSet);
      for (int i = 0; i < nSegments; ++ i) {
         pPartialBondMesh->Draw();
         TrafoView.Translate(FVec3f(0., 0., 1./(1.-fDotScale)));
         TrafoView.Actualize(ShaderSet);
      }
   }
}



void FViewImpl::RenderHalfBond(ct::FAtom const &At0, size_t iAt, ct::FAtom const &At1, size_t jAt, FBondLine const &bl, uint BondFlags, FGeometry &Geometry, FShaderSet &ShaderSet)
{
   (void)At0;
   (void)At1;
   (void)bl;
//    double
//       fMultiBond = 1. + ((std::min(iAt,jAt) % 7)/2.);

   double
      fMultiBond = bl.fBondOrder;
   // if bond order is close to an integer, use nearest integer.
   double
      ThrFullBond = fPartialBondThresh();
   if (std::abs(fMultiBond - std::floor(fMultiBond+.5)) < ThrFullBond)
      fMultiBond = std::floor(fMultiBond+.5);

   int
      nMultiBond = std::ceil(fMultiBond);
   if (nMultiBond == 0)
      return;


   bool
      DottedBond = (BondFlags & BOND_Partial) != 0;

   FElementOptions
      *pElementOptions0 = d->pElementOptions(int(iAt), &*Geometry.pAtoms);
   if (0 == (BOND_Grey & BondFlags)) {
      // set to source atom color.
//       SetGlAtomColor1(ShaderSet.hDiffuseColor, At0);
      SetGlAtomColor1(ShaderSet.hDiffuseColor, pElementOptions0->GetBondColor1());

   } else {
      // set bond line as gray.
      CALL_GL( glUniform4f(ShaderSet.hDiffuseColor, 0.2, 0.2, 0.2, 1.0f) );
   }

   float fBondScaleOuter = this->fBondScaleOuter();
   double fMultiBondScaleFactor = this->v->GetMultiBondScale()/100.; // used to be fixed to 1.2 (120)
   double fMultiBondPosFactor = this->v->GetMultiBondPos()/100.; // used to be fixed to 0.7
   int iMultiBondType = this->v->GetMultiBondType(); // used to be fixed to 0 (linear)

   if (DottedBond) {
      fBondScaleOuter *= 0.2;
   }

   float
      r0 = pElementOptions0->GetDrawRadius();
//       r0 = GetAtomDrawRadius(At0.iElement);

   FBondVisualInfo const
      *pVisInfo = Geometry.FindBondVisualInfo(int(iAt), int(jAt));
   if (pVisInfo == 0)
      return; // not supposed to happen.


   // relative offset of the base of the cylinder from the atomic position.
   double
      fOffs = 0.3 * fAtomScale();

   // get precomputed base orthogonal directions to form the bond line.
   bool
      Flip = pVisInfo->iAt != int(iAt);
   FVec3f
      v0 = Flip ? pVisInfo->vAtPosJ : pVisInfo->vAtPosI,
      vTanU = pVisInfo->vTanU,
      vTanV = pVisInfo->vTanV,
      vDir = pVisInfo->vNorm;
   if (Flip) {
      vDir *= -1.;
      vTanV *= -1.;
   }
   FMat4f
      mView = hStack(vTanU, vTanV, vDir, v0 + float(fOffs * r0) * vDir);

   float
      fHeight = 0.5*pVisInfo->fDist - r0*fOffs;
   TrafoView.Push();
   TrafoView.Multiply(mView);
   // ^- z axis is now along bond direction, with z=0 on base position.


   if (nMultiBond != 0) {
      for (int iMultiBond = 0; iMultiBond < nMultiBond; ++ iMultiBond) {
         TrafoView.Push();
//          double
//             cx, sx,
// //             fMultiBondScale = 0.4,
//             fMultiBondScale = 1.2/double(nMultiBond),
// //             fMultiBondPos = 0.7 * fBondScaleOuter;
//             fMultiBondPos = 0.7 * fBondScaleOuter,
//             fBondWeight = 1.;
         double
            cx, sx,
//             fMultiBondScale = 1.2/double(nMultiBond),
//             fMultiBondScale = 1.2/double(nMultiBond),
// //             fMultiBondPos = 0.7 * fBondScaleOuter;
//             fMultiBondScale = 1.4/double(nMultiBond),
//             fMultiBondPos = 0.7 * fBondScaleOuter,
            fMultiBondScale = fMultiBondScaleFactor/double(nMultiBond),
            fMultiBondPos = fMultiBondPosFactor * fBondScaleOuter,
//             fMultiBondPos = 1.1 * fBondScaleOuter,
            fBondWeight = 1.;
         if (nMultiBond == 1) fMultiBondScale = 1.;
         if (iMultiBondType == 1) {
            // linear alignment
            if (nMultiBond == 1)
               cx = 0.;
            else
               cx = 2. * double(iMultiBond)/double(nMultiBond-1) - 1.;
            sx = 0.;
         } else {
            // radial alignment
            if (nMultiBond == 1) {
               cx = 0.;
               sx = 0.;
            } else {
               double fPhase = 2*M_PI*iMultiBond/double(nMultiBond) * (Flip? -1. : 1.);
               cx = std::cos(fPhase);
               sx = std::sin(fPhase);
            }
         }
         double
            fPartialBond = 1.;
         if (iMultiBond == nMultiBond - 1) {
            fPartialBond = fMultiBond - double(nMultiBond-1);
         }
         if (fPartialBond < 1. - ThrFullBond) {
//             DottedBond = true;
            fBondWeight = fPartialBond;
//             fMultiBondScale *= (.5 + .5 *fBondThinning());
            fMultiBondScale *= (.25 + .75 *fBondThinning());
            // ^- well... it is half correct.
         }

         TrafoView.Translate(FVec3f(fMultiBondPos*cx,fMultiBondPos*sx,0.f));
         DrawBond1(fMultiBondScale*fBondScaleOuter, fHeight, fBondWeight, TrafoView, pBondMesh, pPartialBondMesh, ShaderSet);
         TrafoView.Pop();
      }
   } else {
      DrawBond1(fBondScaleOuter, fHeight, (DottedBond? 0.4 : 1.), TrafoView, pBondMesh, pPartialBondMesh, ShaderSet);
//       if (!DottedBond) {
//          // regular bond
//          TrafoView.Scale(FVec3f(fBondScaleOuter,fBondScaleOuter,fHeight));
//          TrafoView.Actualize(ShaderSet);
//          pBondMesh->Draw();
//          // ^- note: this kind of transformation is highly non-unitary. For this reason
//          //    we do some special care in the tranformation of normals, and actually
//          //    calculate and set the normal transformations properly whenever uploading
//          //    the view transformation (admittedly, it is done in a somewhat hacky way).
//       } else {
//          // bond forming/breaking-bond.
//          int
//             nSegments = 1 + fHeight / 0.7;
//          float
//             fDotScale = 0.4, // smaller: smaller holes.
//             fSegmentHeight = fHeight / (nSegments - fDotScale/2); // with -.5 they join exactly (no space), with -0. the space is twice too large.
//          TrafoView.Scale(FVec3f(fBondScaleOuter,fBondScaleOuter,(1-fDotScale)*fSegmentHeight));
//          TrafoView.Actualize(ShaderSet);
//          for (int i = 0; i < nSegments; ++ i) {
//             pPartialBondMesh->Draw();
//             TrafoView.Translate(FVec3f(0., 0., 1./(1.-fDotScale)));
//             TrafoView.Actualize(ShaderSet);
//          }
//       }
   }
   TrafoView.Pop();

   IR_SUPPRESS_UNUSED_WARNING(jAt);
}


void FViewImpl::RenderFreeLine(FFreeLine *pFreeLine, uint RenderFlags, FShaderSet &ShaderSet)
{
   // this might well be the slowest possible way to draw a line.
   SetGlAtomColor1(ShaderSet.hDiffuseColor, pFreeLine->GetColor());
   bool
      Flip = false;
   FVec3f
      v0(pFreeLine->GetFrom()),
      v1(pFreeLine->GetTo()),
      vNorm = v1 - v0;
   float
      fDist = Length(vNorm);
   if (fDist != 0.)
      vNorm /= fDist;

   // make orthonormal directions.
   FVec3f
      vDir = vNorm,
      vTanU, vTanV;
   MakeOrthDirections(vTanU, vTanV, vDir);
   std::swap(vTanU, vTanV);
   vTanU *= -1.;

   if (Flip) {
      vDir *= -1.;
      vTanV *= -1.;
   }
   FMat4f
      mView = hStack(vTanU, vTanV, vDir, v0);

   float
      fHeight = fDist,
      fWeight = pFreeLine->GetWeight();
   TrafoView.Push();
   TrafoView.Multiply(mView);
   // ^- z axis is now along bond direction, with z=0 on base position.

   DrawBond1(pFreeLine->GetWidth(), fHeight, fWeight, TrafoView, pLineMesh, pPartialBondMesh, ShaderSet);

   TrafoView.Pop();

   (void)ShaderSet; // suppress unused warnings
   (void)RenderFlags;
}


void FViewImpl::RenderFreeLabel(FFreeLabel *pFreeLabel, uint RenderFlags, FShaderSet &ShaderSet)
{
   std::wstring wc = pFreeLabel->GetText().toStdWString();
   if (!wc.empty()) {
      uint32_t dw;
      dw = irgb(pFreeLabel->GetColor()) | 0xff000000;
      FVec3f x(pFreeLabel->GetPos());
      pTextBuffer->Print(vec3f(x[0],x[1],x[2]), float(pFreeLabel->GetSize()), dw, wc.c_str(),
         FGlTextBuffer::HALIGN_Center | FGlTextBuffer::VALIGN_Center);
   }
   (void)ShaderSet; // suppress unused warnings
   (void)RenderFlags;
}


void FViewImpl::RenderFreeObject(FFreeObject *pFreeObject, uint RenderFlags, FShaderSet &ShaderSet)
{
   {
      FFreeLine *pFreeLine = dynamic_cast<FFreeLine*>(pFreeObject);
      if (pFreeLine != 0)
         return RenderFreeLine(pFreeLine, RenderFlags, ShaderSet);
   }

   {
      FFreeLabel *pFreeLabel = dynamic_cast<FFreeLabel*>(pFreeObject);
      if (pFreeLabel != 0)
         return RenderFreeLabel(pFreeLabel, RenderFlags, ShaderSet);
   }
}


void FViewImpl::RenderFreeObjects(FFreeObjectSet &FreeObjects, uint RenderFlags, FShaderSet &ShaderSet)
{
   for (FFreeObject &obj : FreeObjects)
      RenderFreeObject(&obj, RenderFlags, ShaderSet);
}


void FViewImpl::RenderAtomSet(FGeometry &Geometry, uint RenderFlags, FShaderSet &ShaderSet)
{
   FAtomSet
      &Atoms = *Geometry.pAtoms;
   TrafoView.Push();

   if (1) {
      for ( uint iAt = 0; iAt < Atoms.size(); ++ iAt )
         RenderAtom(Atoms[iAt], d->AtomFlags(iAt), iAt, RenderFlags, Geometry, ShaderSet);

      if (RenderFlags & RENDER_Opaque) {
         for ( uint iBond = 0; iBond < Geometry.m_BondLines.size(); ++ iBond ) {
            FBondLine
               bl = Geometry.m_BondLines[iBond]; // make a copy.
            if (bl.fBondOrder < fPartialBondThresh())
               continue;
            bool
               iHidden = d->IsAtomHidden(bl.iAt),
               jHidden = d->IsAtomHidden(bl.jAt),
               BothHidden = iHidden && jHidden,
               OneHidden = (iHidden || jHidden) && !BothHidden;
            if (BothHidden)
               // skip bond lines between two hidden atoms
               continue;
            if (OneHidden && (!v->m_IndicateHiddenAtoms || (bl.Flags != 0)))
               // skip bond lines between hidden and non-hidden atoms, unless
               // either hidden atom connections are supposed to be indicated ("IndicateHiddenAtoms)",
               // or some non-standard flags have been manually set for the current bond line (bl.Flags != 0).
               continue;
            if (!v->GetShowHydrogens() && (Atoms[bl.iAt].iElement == 1 || Atoms[bl.jAt].iElement == 1))
               // (at least) one of the atoms is a hydrogen, and we have the
               // "hide hydrogens completely" flag set (independently of regular hides).
               // (not sure if it would be better to simply put regular hydrogen hide flags,
               // but this is certainly simpler, and does not mess up static settings, if one
               // were just to try to get a better overview over how a molecule looks.)
               continue;

            if (OneHidden && bl.Flags == 0) {
               bl.Flags |= BOND_Partial;
//                bl.Flags |= BOND_Grey;
               // ^- hm.. does this look better or worse without "grey"?
               bl.fBondOrder = 0.5; // <- required to draw the bond as dotted.
            }

            SetGlPickId(ShaderSet.hObjectId, PICKID_BondLine + bl.iAt + PICKID_BondLineStride * bl.jAt);
            RenderHalfBond(Atoms[bl.iAt], bl.iAt, Atoms[bl.jAt], bl.jAt, bl, bl.Flags, Geometry, ShaderSet);
            RenderHalfBond(Atoms[bl.jAt], bl.jAt, Atoms[bl.iAt], bl.iAt, bl, bl.Flags, Geometry, ShaderSet);
         }
         SetGlPickId(ShaderSet.hObjectId, 0);
      }
   }

   // restore input transformations.
   TrafoView.Pop();
   TrafoView.Actualize(ShaderSet);
}

void FVolumeDataSet::BuildRenderCache(FView3d *pView3d)
{
   FIsoSurfaceSettings
      IsoSurfOpt = pView3d->v->MakeIsoSurfaceSettings();
      // ^- make a copy of the view's default settings and patch in
      //    explicitly set information for the current orbital.
   if (pGlMesh.get() == 0) {
      // no mesh yet. Make a mesh for the orbital's iso-surface.
      FVolumeVisualConfig
         *pVis = this->pVisConfig.get();
//       std::cout << fmt::format("*make volume mesh: {} (pVisCfg = {})\n", GetDesc().toStdString(), (void*)(pVis));

//       if (pVis->IsoValues.empty())
//          pVis->IsoValues = IsoSurfOpt.IsoValues;
//       else
//          IsoSurfOpt.IsoValues = pVis->IsoValues;
// //       double fIsoValue = pView3d->m_IsoThreshold/100.;
// //       if (pVis->fIsoValue != -1.)
// //          fIsoValue = pVis->fIsoValue;
//       if (!pVis->bColorSet) {
// //          pVis->SetColorFromCentralHue();
//          pVis->AssignDefaultColor(m_pDocument);
//          IsoSurfOpt.IsoValues = pVis->IsoValues;
//       }
// //       IsoSurfOpt.IsoValues.clear();
// //       IsoSurfOpt.IsoValues.push_back(FIsoThresholdEntry(fIsoValue, pVis->cIsoPlus));
// //       IsoSurfOpt.IsoValues.push_back(FIsoThresholdEntry(-fIsoValue, pVis->cIsoMinus));
      if (!pVis->DetailsAssigned()) {
         pVis->pDetails = MakeDefaultTwoPhaseIsoSurfaceConfig(m_pDocument, pView3d->m_IsoThreshold/100.);
      }
      pVis->pDetails->AssignTo(IsoSurfOpt);

      QApplication::setOverrideCursor(Qt::WaitCursor);
      FIndexedTriangleListPtr
         pSurfCombined = this->MakeIsoSurface(IsoSurfOpt);
      QApplication::restoreOverrideCursor();

      pGlMesh = new FGlMesh(*pSurfCombined);
   }
}


// compute inverse error function erfinv(p) near p=1 by bracketing
// std::erfc (in cmath since C++11).
double ErfInv_Near1(double p) {
   ct::TImplicitInverseFn<double>
      InvErfc(std::erfc, -5, 10);
   return InvErfc(1 - p);
}


static double GetOrbitalSquareMomentFactor(double fIsoThreshold)
{
   using std::sqrt;
   assert(fIsoThreshold >= 0 && fIsoThreshold < 1.);
   if (fIsoThreshold >= 1.)
      fIsoThreshold = 0.999;
   // P(mu - n sigma < x < mu + n sigma) = erf(n/sqrt(2))
   // -> n/sqrt(2) = erfinv(p)
   //    n = sqrt(2) * erfinv(p)
//    return boost::math::erf_inv(fIsoThreshold) * sqrt(2.);
//    // ^-- I guess I could just bracket erfc(1-fIsoThreshold) to get the inverse of this.
//    //     The compute time for this does not really matter, and it would get rid of the
//    //     boost function.
   return ErfInv_Near1(fIsoThreshold) * sqrt(2.);
}



void FViewImpl::RenderVolumeDataSet(FVolumeDataSet &VolumeDataSet, FShaderSet &ShaderSet)
{
   FOrbital
      *pOrb = dynamic_cast<FOrbital*>(&VolumeDataSet);
   if (pOrb == 0 || !v->GetOrbitalEllipsoidsOnly()) {
      if (VolumeDataSet.pGlMesh.get() == 0)
         VolumeDataSet.BuildRenderCache(v);

      // just render the mesh. State changes should, technically,
      // be supplied from the outside.
      if (VolumeDataSet.pGlMesh != 0) {
         CALL_GL( glUniform4f(ShaderSet.hDiffuseColor, 1.0f, 1.0f, 1.0f, 1.0f) );
         VolumeDataSet.pGlMesh->Draw();
      }
   } else {
      // don't render the iso surface, just an ellipsoid at the orbital's center
      if (!pOrb->HaveMoments) {
         IvEmit("!view: can't render orbital ellipsoid---no moments made.");
      } else {
         FVolumeVisualConfig
            *pVisConfig = pOrb->pVisConfig.get();
         // hmhmhm...
         // FIXME: this really should not be here, right?
         if (!pVisConfig->DetailsAssigned()) {
//             pVisConfig->AssignDefaultColor(pOrb->GetDocument());
            pVisConfig->pDetails = MakeDefaultTwoPhaseIsoSurfaceConfig(pOrb->GetDocument(), v->m_IsoThreshold/100.);
            IvEmit("!unexpected absence of visual configuration details in FViewImpl::RenderVolumeDataSet.");
         }

         TrafoView.Push();
         TrafoView.Translate(FVec3f(float(pOrb->vDipMom[0]), float(pOrb->vDipMom[1]), float(pOrb->vDipMom[2])));
         FMat4f
            Scale3x3 = vmath::identity4<float>();
         double
            fBaseScale = GetOrbitalSquareMomentFactor(v->GetIsoThreshold()/100.);
//          IvEmit("  ellipse: iso-thr %1  base-scale %2", v->GetIsoThreshold(), fBaseScale);
         if (v->GetOrientOrbitalEllipsoids()) {
            double mQuadMomEv[9], QuadMomEw[3];
            for (size_t i = 0; i < 3; ++ i)
               for (size_t j = 0; j < 3; ++ j)
                  mQuadMomEv[i+3*j] = pOrb->mQuadMom[i][j];
            Diagonalize(&QuadMomEw[0], &mQuadMomEv[0], 3, 3);

            for (size_t i = 0; i < 3; ++ i)
               for (size_t j = 0; j < 3; ++ j)
                  Scale3x3(i,j) = float(fBaseScale * std::sqrt(QuadMomEw[j]) * mQuadMomEv[i+3*j]);
            if (det(Scale3x3) < 0.) {
               for (size_t i = 0; i < 3; ++ i) {
                  std::swap(Scale3x3(i,1), Scale3x3(i,2));
               }
            }
            TrafoView.Multiply(Scale3x3);
         } else {
            TrafoView.Scale(FVec3f(fBaseScale,fBaseScale,fBaseScale));
         }

         TrafoView.Actualize(ShaderSet);
//          CALL_GL( glUniform4f(ShaderSet.hDiffuseColor, 1.0f, 1.0f, 1.0f, 1.0f) );
         uint32_t dwColor0 = irgb(pOrb->GetBaseColor());
         SetGlAtomColor1(ShaderSet.hDiffuseColor, dwColor0);
         // note: it's drawn with the orbital shader!
         pAtomMesh->Draw();
         TrafoView.Pop();
      }
   }
}

static void Orth1(FVec3f &InOut, FVec3f &vDir)
{
   FVec3f vNorm = vDir;
   vNorm.Normalize();
   InOut -= Dot(InOut, vNorm) * vNorm;
   InOut.Normalize();
}


// convert screen coordinates in terms of pixels (x,y,z) to world coordinates.
FVec3f FViewImpl::ScreenToWorld(FVec3f ScreenCoords)
{
   ScreenCoords[0] = ( 2*(ScreenCoords[0]/width()) - 1.);
   ScreenCoords[1] = (-2*(ScreenCoords[1]/height()) + 1.);

   if (std::abs(TrafoProj.Top()(3,2) + 1) < 1e-5) {
      // undo perspective transformation.
      // for perspective trafo, w = -z, and for ortho trafo w = 1.
      // (this is approximate -- the actual screen space position depends
      //  on the point on the picking ray in this case. we here just assume that
      //  the clicked part is somewhere close to the camera distance)
      float w = fCameraDistActual;
      ScreenCoords[0] *= w;
      ScreenCoords[1] *= w;
   }

   ScreenCoords[2] = 0.;
   vec4f v(ScreenCoords[0], ScreenCoords[1], -ScreenCoords[2], 1.);

   FMat4f
      InvTrafo = vmath::transpose(vmath::inverse(TrafoView.Top() * TrafoProj.Top()));
      // ^- that cannot be right.. is vmath multiplying the matrices in a different
      //    order than I think it should? (row major vs col major?) Or is inverse() broken?
   v = InvTrafo * v;

//    v = vmath::inverse(vmath::transpose(TrafoProj.Top())) * v;
//    v[2] = -1.; // forward
//    v[3] =  0.; // ray, not point.
//    v = vmath::inverse(vmath::transpose(TrafoView.Top())) * v;

   return FVec3f(v[0], v[1], v[2]);
}


long _lame_rint(double f) {
   using namespace std; // without this it doesn't work on g++ 4.8, becaue neither trunc nor round are in namespace std:: in cmath.
   // fun fact: did you know that there is no "round to nearest integer" function in the
   // C++ standard library before C++11?
//    double t = std::trunc(f);
   double t = ::trunc(f);
   double d = f - t;
   if (d > 0.5)
      t += 1.;
   else if (d < -0.5)
      t -= 1.;
   return int(t);
}

int GetFrameMoveDeltaAmount(double fDir, Qt::KeyboardModifiers KeyboardModifiers)
{
   double
      delta = 0.; // nothing if neither 'control' nor 'shift' is pressed.
   if (KeyboardModifiers & Qt::ShiftModifier) {
      delta = fDir; // +1 or -1
   }
   if (KeyboardModifiers & Qt::ControlModifier) {
      delta *= 10.;
   }
   return int(_lame_rint(delta));
}



void FView3d::wheelEvent(QWheelEvent *event) {
//    std::cout << fmt::format("zoom: {:8.5f}   wheel: {}   delta = {}", v->fZoomFactor, (int)event->buttons(), event->delta()) << std::endl;
//    if ( (event->buttons() == (Qt::RightButton)) || (event->buttons()==0 && event->modifiers() == ControlModifier) ) {
   if ( event->buttons() == (Qt::RightButton) ) {
      v->fZoomFactor *= std::pow(1.05f, float(event->delta())/120.f);
      update();
   } else if ( event->buttons() == (Qt::RightButton | Qt::LeftButton) ) {
      float fRate = 5.0;
      if (event->modifiers() & Qt::ControlModifier)
         fRate *= .1f;
      v->RollCamera(fRate * event->delta()/120.f, event->x(), event->y());
      update();
   } else if (event->modifiers() & Qt::ShiftModifier) {
      d->MoveActiveCol(GetFrameMoveDeltaAmount(double(event->delta())/120., event->modifiers()));
   } else {
      event->ignore();
   }
//    std::cout << "\n" + q2s(GetViewDesc()) << std::endl;
}

void FView3d::mousePressEvent(QMouseEvent *event) {
   // remember if this was a simple click---if there is more than one mouse button pressed,
   // or the mouse is moved, it is not.
   v->SimpleClick = (event->buttons() == event->button());
//    std::cout << fmt::format("press: btn = {}  btns = {}   simple? {}", (int)event->button(), (int)event->buttons(), (int)v->SimpleClick) << std::endl;

   v->LastX = event->x();
   v->LastY = event->y();

   if (v->SimpleClick) {
      v->FirstX = v->LastX;
      v->FirstY = v->LastY;
   }
//    grabMouse();
//    QApplication::setOverrideCursor(Qt::BlankCursor);
}

void FView3d::mouseReleaseEvent(QMouseEvent *event) {
//    releaseMouse();
//    QApplication::restoreOverrideCursor();

   int delta = std::max(std::abs(v->FirstX - event->x()), std::abs(v->FirstY - event->y()));
//    if (delta > 4)
//       v->SimpleClick = false;

//    if (event->button() == Qt::LeftButton && v->SimpleClick) {
   if (v->SimpleClick) {
      if (delta <= 4)
         v->ClickPosition(event->x(), event->y(), event->button(), event->modifiers(), event->globalPos());
      else if (event->button() == Qt::LeftButton)
         // is that right? will this not get triggered if left+right moving the mouse?
         v->SelectRect(v->FirstX, v->FirstY, event->x(), event->y(), event->button(), event->modifiers(), event->globalPos());
   }
}



FVec3f Extract3(vec4f const &v) { return FVec3f(v[0], v[1], v[2]); }
vec4f v4(FVec3f const &xyz, float w) { return vec4f(xyz[0], xyz[1], xyz[2], w); }


bool IsOrthoNorm(FVec3f const &vA, FVec3f const &vB) {
   if (!IsEq(vA.LengthSq(), 1.f) || !IsEq(vB.LengthSq(), 1.f))
      return false;
   if (!IsEq(Dot(vA, vB), 0.f))
      return false;
   return true;
}

static FMat4f MakeRotMatrix(uint iAxisX, uint iAxisY, float fAngleXY)
{
   FMat4f
      mRot = vmath::identity4<float>();
   float
      cs = std::cos(fAngleXY/180. * M_PI),
      ss = std::sin(fAngleXY/180. * M_PI);
   mRot(iAxisX,iAxisX) = cs;
   mRot(iAxisY,iAxisX) = -ss;
   mRot(iAxisX,iAxisY) = ss;
   mRot(iAxisY,iAxisY) = cs;
   return mRot;
}

// make a rotation matrix for rotating in the plane spanned by vAxisX and vAxisY around angle fAngleXY,
// around the center of the coordinate system.
static FMat4f MakeRotMatrix(FVec3f const &vAxisX, FVec3f const &vAxisY, float fAngleXY)
{
   assert(IsOrthoNorm(vAxisX, vAxisY));
   FVec3f
      vAxisZ = Cross(vAxisX, vAxisY);
   FMat4f
      mAxis = hStack(vAxisX, vAxisY, vAxisZ, FVec3f(0.f,0.f,0.f));
   // ^- FIXME: this is supposed to be transposed. get order of stuff in vector_math
   //    fixed (or consistent with what I think it should be).
   return vmath::transpose(mAxis) * MakeRotMatrix(0, 1, fAngleXY) * mAxis;
}


void FViewImpl::RollCamera(float fAngle, int x, int y)
{
   // find the axis around which to rotate
   FVec3f
      v0,
      vAxis = vCameraDir, // <- that is not quite right in the perspective case. should use
                          //    difference of ScreenToWorld with different z then.
      vOrth0, vOrth1;

   vAxis.Normalize();
   MakeOrthDirections(vOrth0, vOrth1, vAxis);

   if (x != -1 && true){
      // rotate around mouse position
      v0 = ScreenToWorld(FVec3f(float(x),float(y),0.0));
   } else {
      // rotate around center of view
      v0 = vCameraPos;
   }

   v0 = v0 - Dot(v0, vAxis)*vAxis;
//    std::cout << fmt::format("v0: ({:.4f}, {:.4f}, {:.4f})\n", v0[0], v0[1], v0[2]);

   FMat4f
      mRot = MakeRotMatrix(vOrth0, vOrth1, -fAngle);

   vCameraDir = Extract3(mRot * v4(vCameraDir, 0.));
   vCameraUp = Extract3(mRot * v4(vCameraUp, 0.));
   vCameraPos = Extract3(mRot * v4(vCameraPos - v0, 0.)) + v0;

   // reorthonormalize view system; This *should* not be required, but somehow is.
   vCameraDir.Normalize();
   Orth1(vCameraUp, vCameraDir);

   NormalizeCameraPos();
}


void FView3d::mouseMoveEvent(QMouseEvent *event) {
   if ( event->buttons() != 0 ) {
      float fScale = 20./this->width();
      float fPosScale = v->fZoomFactor * fCameraDist/30.;
      float fDeltaX = fScale * (event->x() - v->LastX);
      float fDeltaY = fScale * (event->y() - v->LastY);
      FVec3f
         vRight = Cross(v->vCameraDir, v->vCameraUp);

      if ( event->buttons() == (Qt::LeftButton | Qt::RightButton) ) {
         if (0) {
            v->vCameraPos -= (fPosScale * fDeltaY) * v->vCameraUp;
            v->vCameraPos += (fPosScale * fDeltaX) * vRight;
         } else {
            FVec3f
               v0 = v->ScreenToWorld(FVec3f(event->x(),event->y(),1.0)),
               v1 = v->ScreenToWorld(FVec3f(v->LastX,v->LastY,1.0)),
               dv = v1 - v0;
//             dv -= Dot(dv,v->vCameraDir) * v->vCameraDir;
            v->vCameraPos += dv;
         }
         // put position up to a certain distance from the cam.
         v->NormalizeCameraPos();
         update();
      }

      if ( event->buttons() == (Qt::MiddleButton) ) {
         v->vCameraUp -= (.1f * fDeltaX) * v->vCameraUp;
         v->vCameraUp += (.1f * fDeltaX) * vRight;
         v->vCameraUp.Normalize();
         Orth1(v->vCameraDir, v->vCameraUp);
         update();
      }

      if ( event->buttons() == (Qt::RightButton) ) {
         float
            fRotScale = 0.1f * 180.f/M_PI;
         FVec3f
            &vDir = v->vCameraDir,
            &vUp = v->vCameraUp,
            &vPos = v->vCameraPos,
            vRight = Cross(vDir, vUp);
         FVec3f
            vRotCenter;
         if (0 == (event->modifiers() & Qt::ShiftModifier)) {
            // rotate view around object
            vRotCenter = FVec3f(0.f, 0.f, 0.f);
         } else {
            // rotate camera (rotate view around camera position)
            vRotCenter = vPos;
         }
         FMat4f
            mRotX = MakeRotMatrix(vRight, vDir, -fRotScale * fDeltaX),
            mRotY = MakeRotMatrix(vUp, vDir, fRotScale * fDeltaY),
            mTot = mRotY * mRotX;
         vDir = Extract3(mTot * v4(vDir, 0.));
         vUp = Extract3(mTot * v4(vUp, 0.));
         vPos = Extract3(mTot * v4(vPos - vRotCenter, 0.)) + vRotCenter;

         // reorthonormalize view system; This *should* not be required, but somehow is.
         v->vCameraDir.Normalize();
         Orth1(v->vCameraUp, v->vCameraDir);
         // ^- hmmm... something strange is going on. can the combination of orthogonal
         //    transforms really cause numerical  problems? it seems to do.
         v->NormalizeCameraPos();
         update();
      }
   }
   v->LastX = event->x();
   v->LastY = event->y();
}


void FView3d::keyPressEvent(QKeyEvent* event) {
//   IvEmit("// FView3d::keyPressEvent  key = %1  nativeModifiers = %2  text = %3", event->key(), event->nativeModifiers(), event->text());
//    std::cout << fmt::format("zoom: {:8.5f}", v->fZoomFactor) << std::endl;
//    switch(event->key()) {
//    case Qt::Key_Escape:
//       close();
//       break;
//    case Qt::Key_PageUp: {
//       v->fZoomFactor *= std::pow(1.05f, +1.f);
//       update();
//       break;
//    }
//    case Qt::Key_PageDown: {
//       v->fZoomFactor *= std::pow(1.05f, -1.f);
//       update();
//       break;
//    }
//    default:
//       event->ignore();
//       break;
//    }
// ^- I don't think this is ever called.

   switch(event->key()) {
//       case Qt::Key_Escape:
//          close();
//          break;
      case Qt::Key_PageUp: {
         v->fZoomFactor *= std::pow(1.05f, +1.f);
         update();
         break;
      }
      case Qt::Key_PageDown: {
         v->fZoomFactor *= std::pow(1.05f, -1.f);
         update();
         break;
      }
      case Qt::Key_Left:
      case Qt::Key_Right: {
         double dir = (event->key() == Qt::Key_Left? -1. : +1.);
         int delta = GetFrameMoveDeltaAmount(dir, event->modifiers());
//          IvNotify(NOTIFY_Information, IvFmt("IvView3D:: key press dir = %1 delta = %2", dir, delta));
         if (delta != 0) {
            d->MoveActiveCol(delta);
         } else {
            FBase::keyPressEvent(event);
         }
         break;
      }
      default:
         FBase::keyPressEvent(event);
         // ^- this propagates to main form event handler if this widget does not
         // do anything with the key. 'event->ignore()' is not good here, unless
         // we wish to eat the event.
         // event->ignore();
         break;
   }
}


void FView3d::updateData(const QModelIndex &/*topleft*/, const QModelIndex &/*bottomright*/)
{
//    std::cout << "called FView3d::updateData." << std::endl;
   update();
}


static FSelectionMode GetModifierSelectMode(Qt::KeyboardModifiers KeyboardModifiers)
{
   FSelectionMode
      SelectMode = SELECT_Select;
   if (KeyboardModifiers & Qt::ShiftModifier)
      SelectMode = SELECT_Toggle;
//       if (KeyboardModifiers & Qt::AltModifier)
//          SelectMode = SELECT_Toggle;
   return SelectMode;
}


QPoint FViewImpl::iPickBufferPos(int x, int y)
{
   if (x < 0)
      x = 0;
   if (y < 0)
      y = 0;
   // translate window position into picking buffer position.
   int
      pw = pPickFbo->Size.width(),
      ph = pPickFbo->Size.height(),
      px = (x * pw) / v->width(),
      py = (y * ph) / v->height();
   if (px >= pw) px = pw - 1;
   if (py >= ph) py = ph - 1; // maybe possible due to rounding? didn't check.

   py = ph - (py + 1); // OpenGL's (0,0) is the LOWER left corner of the buffer.
   return QPoint(px,py);
}


void FViewImpl::SelectRect(int x0, int y0, int x1, int y1, Qt::MouseButton MouseButton, Qt::KeyboardModifiers KeyboardModifiers, QPoint globalPos)
{
   (void)globalPos;
   bool Verbose = false;
   if (MouseButton != Qt::LeftButton)
      // otherwise also called for view transformation stuff, where we do not want it.
      return;

   if (x0 > x1) std::swap(x0,x1);
   if (y0 < y1) std::swap(y0,y1);
   // ^- note: x0 must be the LOWER x coordinate, and y0 the HIGHER y coordinate,
   //    because OpenGL begins counting in the lower-left, not the top-left.

   RenderPickBuffer(true); // keep picking buffer bound.

   // convert mouse coords to downsampled coordinate system of pick buffer
   QPoint
      p0 = iPickBufferPos(x0,y0),
      p1 = iPickBufferPos(x1,y1);
   int
      pw = p1.x() - p0.x(),
      ph = p1.y() - p0.y();
   if (pw <= 0 || ph <= 0)
      // empty selection rectangle selected---do nothing.
      return;
   assert_rt(pw > 0 && ph > 0);

   TArray<uint32_t>
      PickData;
   PickData.resize(pw*ph, 0xffffffff);
   if (Verbose)
      IvEmit("SelectRect xy = (%1 / %2)  wh = (%3 / %4)", p0.x(), p0.y(), pw, ph);
   CALL_GL( glReadPixels(p0.x(), p0.y(), pw, ph, GL_RGBA, GL_UNSIGNED_BYTE, &PickData[0]) );

   if (Verbose)
      IvEmit("SelectRect(%1,%2,%3,%4) -> got %5 pixels.", x0, y0, x1, y1, PickData.size());

   // find unique object IDs in the given range.
   typedef std::set<uint32_t>
      FIdSet;
   FIdSet
      PickedIds;
   for (size_t iDw = 0; iDw != PickData.size(); ++ iDw)
      PickedIds.insert(PickData[iDw]);

   FIdSet::const_iterator
      it;
   FSelectionMode
      SelectMode = GetModifierSelectMode(KeyboardModifiers);
   for (it = PickedIds.begin(); it != PickedIds.end(); ++ it) {
      uint32_t
         dwPickData = *it,
         dwPickKey = dwPickData & PICKID_KeyMask,
         dwPickValue = dwPickData & PICKID_ValueMask;
      if (Verbose)
         IvEmit("   object id:   %1 / %2", dwPickKey, dwPickValue);

      if (dwPickKey == PICKID_Atom) {
         d->SelectAtom(dwPickValue, SelectMode);
         if (SelectMode == SELECT_Select)
            SelectMode = SELECT_Add;
      }
   }
}

void FViewImpl::ClickPosition(int x, int y, Qt::MouseButton MouseButton, Qt::KeyboardModifiers KeyboardModifiers, QPoint globalPos)
{
   if (x < 0 || y < 0 || x >= v->width() || y >= v->height())
      return;
   RenderPickBuffer(true); // keep picking buffer bound.

   uint32_t
      dwPickData = 0xffffffff;
   {
      // convert mouse coords to downsampled coordinate system of pick buffer
      QPoint
         p0 = iPickBufferPos(x, y);
      CALL_GL( glReadPixels(p0.x(), p0.y(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &dwPickData) );
   }

//    std::cout << fmt::format("! ObjectId at ({:3},{:3}): 0x{:08x}"), x, y, dwPickData) << std::endl;
   uint32_t
      dwPickKey = dwPickData & PICKID_KeyMask,
      dwPickValue = dwPickData & PICKID_ValueMask;

   if (MouseButton == Qt::LeftButton) {
      // left button -- select stuff.
      FSelectionMode
         SelectMode = GetModifierSelectMode(KeyboardModifiers);

      if (dwPickKey == PICKID_DataRow && SelectMode == SELECT_Select) {
         d->SetActiveRow(dwPickValue);
      }

      if (dwPickKey == PICKID_Atom) {
         d->SelectAtom(dwPickValue, SelectMode);
      }

      if (dwPickData == 0 && SelectMode == SELECT_Select)
         d->UnselectAll(true);
   } else if (MouseButton == Qt::RightButton) {
      if (dwPickKey == PICKID_Atom) {
         bool
            SelectedForMenu = false;
         if (!d->IsAtomSelected(dwPickValue)) {
            d->SelectAtom(dwPickValue, SELECT_Select);
            SelectedForMenu = true;
         }

         QAction
            ActCalcCharges("Make Charges", v),
            ActFindOrbitals("Find Orbitals", v),
            ActHideAtoms("Hide Atoms", v),
            ActMakeHybrids("Make Hybrids", v),
            ActCopyAtomNumbers("Copy Atom Numbers", v),
            ActAddBonds0("Add Bond Lines (regular)", v),
            ActAddBonds1("Add Bond Lines (dotted)", v);
         d->connect(&ActHideAtoms, SIGNAL(triggered()), d, SLOT(HideSelectedAtoms()));
         d->connect(&ActCalcCharges, SIGNAL(triggered()), d, SLOT(CalcChargesOnSelectedAtoms()));
         d->connect(&ActFindOrbitals, SIGNAL(triggered()), d, SLOT(FindOrbitalsOnSelectedAtoms()));
         d->connect(&ActMakeHybrids, SIGNAL(triggered()), d, SLOT(MakeHybridsForSelectedAtomGroup()));
         d->connect(&ActCopyAtomNumbers, SIGNAL(triggered()), d, SLOT(CopySelectedAtomNumbers()));
         d->connect(&ActAddBonds0, SIGNAL(triggered()), d, SLOT(MakeBondLinesForSelectedAtoms()));
         d->connect(&ActAddBonds1, SIGNAL(triggered()), d, SLOT(MakeBondLinesForSelectedAtoms()));
         ActAddBonds1.setData(QVariant(QString("gray|dotted")));
         QMenu
            menu(v);
         menu.addAction(&ActFindOrbitals);
         menu.addAction(&ActCalcCharges);
         menu.addAction(&ActCopyAtomNumbers);
         menu.addSeparator();
         menu.addAction(&ActHideAtoms);
         if (d->nSelectedAtoms() > 1) {
            menu.addAction(&ActAddBonds0);
            menu.addAction(&ActAddBonds1);
         }
//          if (0) {
//             // add measures etc.
//             menu.addSeparator();
//          }
         QAction *pAct = menu.exec(globalPos);
         if (pAct == 0 && SelectedForMenu)
            // we selected an atom just for the menu click, but the user did not actually
            // do anything with it. Unselect it again.
            d->UnselectAll();
      } else if (dwPickKey == PICKID_DataRow) {
         d->SetActiveRow(dwPickValue);
         QAction
            ActToggleRow("Hide Orbital", v);
         QMenu
            menu(v);
         d->connect(&ActToggleRow, SIGNAL(triggered()), d, SLOT(ToggleActiveDataRow()));
         menu.addAction(&ActToggleRow);
         menu.exec(globalPos);
      } else if (dwPickKey == PICKID_BondLine) {
         int
            iAt = dwPickValue % PICKID_BondLineStride,
            jAt = dwPickValue / PICKID_BondLineStride;
         double
            fCurrentBo = -1.;

         // check if we currently have a bond order set---if yes, retrieve it.
         IFrame
            *pFrame = d->GetCurrentFrame();
         if (pFrame) {
            FGeometry
               *pGeometry = pFrame->pGetGeometry();
            if (pGeometry) {
               for (size_t iBondLine = 0; iBondLine != pGeometry->m_BondLines.size(); ++ iBondLine) {
                  FBondLine const &bl = pGeometry->m_BondLines[iBondLine];
                  if ((bl.iAt == iAt && bl.jAt == jAt) || (bl.iAt == jAt && bl.jAt == iAt)) {
                     fCurrentBo = bl.fBondOrder;
                  }
               }
            }
         }

         FBondChangeAction
            ActHideBond(iAt, jAt, FBondChangeAction::ACTION_Hide, "Hide Bond", v),
            ActPartializeBond(iAt, jAt, FBondChangeAction::ACTION_SetStyleDotted, "Set Style: Dotted", v),
            ActResetBond(iAt, jAt, FBondChangeAction::ACTION_Reset, "Reset Bond", v),
            ActSetBondOrder05(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 0.5", v, 0.5),
            ActSetBondOrder10(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 1.0", v, 1.0),
            ActSetBondOrder15(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 1.5", v, 1.5),
            ActSetBondOrder20(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 2.0", v, 2.0),
            ActSetBondOrder25(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 2.5", v, 2.5),
            ActSetBondOrder30(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 3.0", v, 3.0),
            ActSetBondOrder40(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 4.0", v, 4.0),
            ActSetBondOrder50(iAt, jAt, FBondChangeAction::ACTION_SetBondOrder, "Set Bond Order: 5.0", v, 5.0);
         QMenu
            CurrentBoInfo(s2q(fmt::format("[BO = {:.4f}]", fCurrentBo)), v);
         QMenu
            menu(v);
         d->connect(&ActHideBond, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActPartializeBond, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActResetBond, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder05, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder10, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder15, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder20, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder25, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder30, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder40, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         d->connect(&ActSetBondOrder50, SIGNAL(triggered()), d, SLOT(ChangeSelectedBond()));
         menu.addAction(&ActHideBond);
         menu.addAction(&ActPartializeBond);
         menu.addAction(&ActResetBond);
         menu.addSeparator();
         if (fCurrentBo != -1.)
            menu.addMenu(&CurrentBoInfo);
         menu.addAction(&ActSetBondOrder05);
         menu.addAction(&ActSetBondOrder10);
         menu.addAction(&ActSetBondOrder15);
         menu.addAction(&ActSetBondOrder20);
         menu.addAction(&ActSetBondOrder25);
         menu.addAction(&ActSetBondOrder30);
         menu.addAction(&ActSetBondOrder40);
         menu.addAction(&ActSetBondOrder50);
         menu.exec(globalPos);
      }
   }


   pMainFbo->Bind();
}







// script interface
ICamera::ICamera(QObject *pView3d_)
   : QObject(pView3d_)
{
   this->setObjectName("camera");
}

IView3d *ICamera::pView3d()
{
   return dynamic_cast<IView3d*>(parent());
}

void ICamera::set_pos(float x, float y, float z) {
   return pView3d()->set_camera_pos(x,y,z);
}
void ICamera::set_dir(float x, float y, float z) {
   return pView3d()->set_camera_dir(x,y,z);
}
void ICamera::set_vup(float x, float y, float z) {
   return pView3d()->set_camera_vup(x,y,z);
}
void ICamera::set_zoom(float z) {
   return pView3d()->set_camera_zoom(z);
}

IView3d::IView3d(QWidget *parent_, FDocument *document_)
   : FView3d(parent_, document_)
{
   new ICamera(this);
   // ^- my understanding is that this automatically adds the new object
   // to *this's children, and that is then automatically destroyed when
   // the current object is.
//    std::cout << fmt::format("!! IView3d/camera: '{}'", this->findChild<QObject*>("camera")->objectName().toStdString()) << std::endl;
};


void IView3d::set_camera_pos(float x, float y, float z) {
   v->vCameraPos = FVec3f(x,y,z);
   v->NormalizeCameraPos();
   update();
}

void IView3d::set_camera_dir(float x, float y, float z){
   v->vCameraDir = FVec3f(x,y,z);
   v->vCameraDir.Normalize();
   Orth1(v->vCameraUp, v->vCameraDir);
   update();
}

void IView3d::set_camera_vup(float x, float y, float z){
   v->vCameraUp = FVec3f(x,y,z);
   v->vCameraUp.Normalize();
   Orth1(v->vCameraDir, v->vCameraUp);
   update();
}

void IView3d::set_camera_zoom(float z){
   v->fZoomFactor = z;
   update();
}

void IView3d::set_option(QString const &OptionName, QVariant f)
{
   // that's a compatibility function for scripts setting properties in
   // the previously used way.
   QByteArray
      baOptionName = OptionName.toLocal8Bit();
   baOptionName.replace("-","_");
   if (baOptionName == "allow_iso_flip")
      baOptionName = "iso_auto_flip";
   if (!setProperty(baOptionName, f))
      IvNotify(NOTIFY_Warning, QString("Script tried to set view option '%1', but this option is not recognized.").arg(OptionName));
}



void FViewImpl::AlignView(FVec3f Right, FVec3f Up)
{
   Right.Normalize();
   Orth1(Up, Right);
   Up.Normalize();
   FVec3f
      Dir = Cross(Up, Right);
   vCameraDir = Dir;
   vCameraUp = Up;
   vCameraPos = -1.0f * Dir;
   NormalizeCameraPos();
   v->update();
}

void IView3d::modify(QString const &How, float f)
{
   IvEmit("Modify View: '%1' [%2]", How, f);
//    << fmt::format("modv: {} = {}", q2s(How), f) << std::endl;
   if (How == "align xy")
      v->AlignView(FVec3f(1.,0.,0.), FVec3f(0.,1.,0.));
   else if (How == "align xz")
      v->AlignView(FVec3f(1.,0.,0.), FVec3f(0.,0.,-1.));
   else if (How == "align yz")
      v->AlignView(FVec3f(0.,1.,0.), FVec3f(0.,.0,1.));
   else if (How == "roll") {
      v->RollCamera(-f,-1,-1);
      update();
   }
}


QRect GetBoundingBox(QImage &img)
{
   // c/p'd from http://stackoverflow.com/questions/3720947/does-qt-have-a-way-to-find-bounding-box-of-an-image
   int
      l = img.width(),
      r = 0,
      t = img.height(),
      b = 0;
   for (int y = 0; y < img.height()-1; ++ y )
      // ^-   -1: for some reason there sometimes are errors in the last row,
      //    which can be marked as filled wrongly. might be a bug in the post-processing
      //    step. for the moment we simply ignore the last row and hope that is sufficient.
   {
      QRgb
         *row = (QRgb*)img.scanLine(y);
      bool
         rowFilled = false;
      for (int x = 0; x < img.width(); ++ x) {
//          if (qAlpha(row[x])) {
//          if (x == 200)
//             std::cout << fmt::format("x = {:4} y = {:4}  col = {:8x}\n", x, y, row[x]);
         if ((row[x]&0xffffff) != 0xffffff) { // <- background color white?
            rowFilled = true;
            r = std::max(r, x);
            if (l > x) {
               l = x;
               x = r; // shortcut to only search for new right bound from here
            }
//             if (y==523)
//                std::cout << fmt::format("* col/row filled: y = {},  x = {}\n", y, x);
         }
      }
      if (rowFilled) {
         t = std::min(t, y);
         b = y;
//          std::cout << fmt::format("* row filled: y = {}\n", y);
      }
//       std::cout << fmt::format("* y = {:4}  -> lt=({:4} {:4}) rb=({:4} {:4})\n", y, l, t, r, b);
   }
//    std::cout << fmt::format("* bounding box: lt=({:4} {:4}) rb=({:4} {:4})\n", l, t, r, b);
   return QRect(QPoint(l,t), QPoint(r,b));
}

// convert from what QT thinks is pre-multiplied alpha with a black background
// to what it actually is---premultiplied alpha with white background. Convert
// image format to ARGB32 exlicit instead of Format_ARGB32_Premultiplied.
QImage ConvertToExplicitAlpha(QImage &img)
{
   QImage
      out(img.size(), QImage::Format_ARGB32);
   int
      w = img.width(),
      h = img.height();
   for (int y = 0; y < h; ++ y) {
      QRgb
         *pRowIn = (QRgb*)img.scanLine(y),
         *pRowOut = (QRgb*)out.scanLine(y);
      for (int x = 0; x < w; ++ x) {
         // current color is blended with our background (white)
         // "unblend" the color component to restore original
         // without the background. Idea is that this will be combined
         // with some other background again, if the image is used elsewhere.
         // (there is some loss of precision here, but I really don't want
         //  to render everything with a black background in IboView...it's
         //  not pretty).
         uint32_t p = pRowIn[x];
         float a = ((p >> 24)&0xff)/255.f;
#define TX(comp)  (a==0? 0xff : int(255.*(1.f - std::min(((1.f - (comp)/255.f)/a),1.f))))
         uint32_t q =  (  ((p >> 24)&0xff) << 24) |  // a
                       (TX((p >> 16)&0xff) << 16) |  // r
                       (TX((p >>  8)&0xff) <<  8) |  // g
                       (TX((p >>  0)&0xff) <<  0);   // b
         pRowOut[x] = q;
#undef TX
      }
   }
   return out;
}


void IView3d::save_png(QString const &FileName){
   QCoreApplication::processEvents(); // yes.. I know. UI coding elegance at its best.
   update();
   QCoreApplication::processEvents();

   IvEmit("* write png '%1'", FileName);
   bool
      WriteAlpha = m_SaveAlpha,
      Crop = m_CropImages;
   QImage
      img = this->grabFrameBuffer(WriteAlpha);
//       img = this->renderPixmap(200,200,false).toImage();
      // ^- doesn't work... calls initializeGL again, but does NOT make a new
      //    view3d object! That means that all the resident objects (off-screen FBOs)
      //    are overwritten!
   if (WriteAlpha)
      // unblend the background color and convert from supposed premultiplied alpha
      // to explicit alpha (we have what is effectively premultiplied alpha in the frame buffer,
      // but it is premultiplied with a different background color than QT thinks).
      img = ConvertToExplicitAlpha(img);
   QImage
      img_out;
   if (!Crop) {
      img_out = img;
   } else {
      img_out = img.copy(GetBoundingBox(img));
   }
   if (FileName == ":/!clipboard!") {
      QClipboard *clipboard = QApplication::clipboard();
#ifndef _MSC_VER
      clipboard->setImage(img_out);
#else
      // on Windows, just using setImage loses the alpha channel.
      // I'll hack around it using this option:
      // http://stackoverflow.com/questions/1260253/how-do-i-put-an-qimage-with-transparency-onto-the-clipboard-for-another-applicati
      QMimeData* mimeData = new QMimeData();

      QByteArray data;
      QBuffer buffer(&data);
      buffer.open(QIODevice::WriteOnly);
      img_out.save(&buffer, "PNG");
      buffer.close();
      mimeData->setData("PNG", data);
      clipboard->setMimeData(mimeData);
      // wheee, it works \o/.
#endif
   } else {
//       IvEmit("!!calling img_out.save('%1')", FileName);
      img_out.save(FileName); //, const char *format = 0, int quality = -1 ) c
   }
}



QString FView3d::GetViewDesc()
{
   std::stringstream str;
   str << q2s(GetOptionsDesc());

   str << fmt::format("view.set_size({},{});\n", width(), height());
   str << fmt::format("view.camera.set_pos({:.5f},{:.5f},{:.5f});", v->vCameraPos[0], v->vCameraPos[1], v->vCameraPos[2]);
   str << fmt::format("\nview.camera.set_dir({:.5f},{:.5f},{:.5f});", v->vCameraDir[0], v->vCameraDir[1], v->vCameraDir[2]);
   str << fmt::format("\nview.camera.set_vup({:.5f},{:.5f},{:.5f});", v->vCameraUp[0], v->vCameraUp[1], v->vCameraUp[2]);
   str << fmt::format("\nview.camera.set_zoom({:.5f});", v->fZoomFactor);
   return s2q(str.str());
}

QWidget *FViewImpl::FindTopmostParent()
{
   QWidget *p = v;
   while (qobject_cast<QWidget*>(p->parent()) != 0)
      p = qobject_cast<QWidget*>(p->parent());
   return p;
}

void IView3d::set_size(int width_, int height_)
{
   QWidget
      *pMainWindow = v->FindTopmostParent();
   QSize
      CurrentSize = size(),
      ParentSize = pMainWindow->size(),
      NewParentSize = QSize(width_ + (ParentSize.width() - CurrentSize.width()),
                            height_ + (ParentSize.height() - CurrentSize.height()));
   pMainWindow->resize(NewParentSize);
}

bool FView3d::RenderAnyLabels() const
{
   return m_LabelAtomNumbers || m_LabelElements || m_LabelElementsC || m_LabelElementsH;
}



#include "prop_FView3d.cpp.inl"
