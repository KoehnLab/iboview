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

#include "IvGl.h"

#include <map>
#include <iostream>
#include <QRect>
#include <QImage>
#include <QApplication>
#include <QClipboard>
#include <QFile>
#include <QTextStream>
#include <QDir>

#include <stdio.h>
#include <stdexcept>
#include <sstream>
#include <cstdlib> // for rand

#include "CxColor.h"
// #include "CxOsInt.h"
#include "CxVec3.h"
using namespace ct;

int g_GlVersion = 30;
#ifdef __APPLE__
   // well, that may be overkill... it's probably not *all* apples, just some, which have this particular
   // driver bug. should be controllable via command line argument, once we introduce sensible handling of those.
   // (hm.. there is a QCommandLineParser. That might be useful).
   bool g_WorkAroundGlFrontFacingBug = true;
#else
   bool g_WorkAroundGlFrontFacingBug = false;
#endif
bool g_WorkAroundAlphaCompositing = false;


FMat4f hStack(FVec3f const &Row0, FVec3f const &Row1, FVec3f const &Row2, FVec3f const &Row3)
{
   FMat4f r;
   r(0,0) = Row0[0]; r(1,0) = Row0[1]; r(2,0) = Row0[2]; r(3,0) = 0.;
   r(0,1) = Row1[0]; r(1,1) = Row1[1]; r(2,1) = Row1[2]; r(3,1) = 0.;
   r(0,2) = Row2[0]; r(1,2) = Row2[1]; r(2,2) = Row2[2]; r(3,2) = 0.;
   r(0,3) = Row3[0]; r(1,3) = Row3[1]; r(2,3) = Row3[2]; r(3,3) = 1.;
   return r;
}

// vec4f vertices of a full screen quad given as two triangles.
static const GLfloat s_QuadPoints[] = {
   -1.0f, -1.0f, 0.0f, 1.0f,
    1.0f, -1.0f, 0.0f, 1.0f,
   -1.0f,  1.0f, 0.0f, 1.0f,
    1.0f, -1.0f, 0.0f, 1.0f,
    1.0f,  1.0f, 0.0f, 1.0f,
   -1.0f,  1.0f, 0.0f, 1.0f
};

void FFullScreenQuad::Init()
{
   // make a vertex descriptor-- we have just positions.
   CALL_GL( glGenVertexArrays(1, &hVertexDesc) );
   CALL_GL( glBindVertexArray(hVertexDesc) );

   CALL_GL( glGenBuffers(1, &IdVertexBuffer) );
   CALL_GL( glBindBuffer(GL_ARRAY_BUFFER, IdVertexBuffer) )
   CALL_GL( glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*4, 0) ); // must be done after GL_ARRAY_BUFFER is bound.
   CALL_GL( glBufferData(GL_ARRAY_BUFFER, sizeof(s_QuadPoints), s_QuadPoints, GL_STATIC_DRAW) )
   CheckGlError("FullScreenQuadInit");
}

FFullScreenQuad::~FFullScreenQuad()
{
   CALL_GL( glDeleteVertexArrays(1, &hVertexDesc) );
   CALL_GL( glDeleteBuffers(1, &IdVertexBuffer) );
   hVertexDesc = 0;
   IdVertexBuffer = 0;
}

void FFullScreenQuad::Draw(GLenum Program)
{
   CALL_GL( glUseProgram(Program) );
   CALL_GL( glBindVertexArray(hVertexDesc) );
   GLuint
      IdPos = glGetAttribLocation(Program, "in_Position");
   CALL_GL( glEnableVertexAttribArray(IdPos) );
   CALL_GL( glBindBuffer(GL_ARRAY_BUFFER, IdVertexBuffer) );
   CALL_GL( glVertexAttribPointer(IdPos, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*4, 0) );
   CALL_GL( glDrawArrays(GL_TRIANGLES, 0, sizeof(s_QuadPoints)/sizeof(s_QuadPoints[0])) );
   CheckGlError("FFullScreenQuad::Draw");
}

FFrameBufferAttachment::FFrameBufferAttachment(char const *pDesc_, QSize Size_, GLenum Attachment_, uint nSamples_, GLuint InternalFormat, GLuint Format, GLenum Type)
{
   CheckGlError(pDesc, "FFrameBufferAttachment/c'tor (enter)");
   nSamples = nSamples_;
   GlAttachment = Attachment_;
   pDesc = pDesc_;
   Size = Size_;

   if (nSamples == 0)
      GlTarget = GL_TEXTURE_2D;
   else
      GlTarget = GL_TEXTURE_2D_MULTISAMPLE;

   CALL_GL( glGenTextures(1, &this->hTexture) );
   CALL_GL( glBindTexture(GlTarget, hTexture) );
   if (nSamples == 0) {
//       glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_FALSE);
//       CheckGlError(pDesc, "glTexParameteri/GL_GENERATE_MIPMAP");
      // ^- not present in core profile, apparently, and anyway disabled by default (use glGenerateMipmap​.)
      CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST) );
      CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST) );
      CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE) );
      CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE) );
      CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0) ); // <- otherwise lots of "Waste of memory: Texture 0 has mipmaps, while its min filter is inconsistent with mipmaps."-warnings.
      CALL_GL( glTexImage2D(GlTarget, 0, InternalFormat, width(), height(), 0, Format, Type, 0) );
   } else {
      CALL_GL( glTexImage2DMultisample(GlTarget, nSamples, InternalFormat, width(), height(), GL_FALSE) );
   }
   CALL_GL( glBindTexture(GlTarget, 0) ); // unbind.
   CheckGlError(pDesc, "FFrameBufferAttachment/c'tor (leave)");
   // note: this makes the attachment. It is, however, not yet attached to anything.
   // that only happens in FFrameBuffer::Attach.
}

FFrameBufferAttachment::~FFrameBufferAttachment()
{
   CALL_GL( glBindTexture(GL_TEXTURE_2D, 0) );
   // ^- on Macs the glDeleteTextures gives "INVALID OPERATION". I don't think it is supposed to.
   //    but I read online that some other people had luck with unbinding the attachments first.
   CALL_GL( glDeleteTextures(1, &hTexture) );
}

void FFrameBufferAttachment::BindAsTexture(GLenum GlTexture)
{
   CheckGlError(pDesc, "FFrameBufferAttachment::BindAsTexture (enter)");
   CALL_GL( glActiveTexture(GlTexture) );
   CALL_GL( glBindTexture(GlTarget, hTexture) );
}

void FFrameBufferAttachment::UnbindAsTexture(GLenum GlTexture)
{
   CheckGlError(pDesc, "FFrameBufferAttachment::Unbind (enter)");
   CALL_GL( glActiveTexture(GlTexture) );
   CALL_GL( glBindTexture(GlTarget, 0) );
}

FFrameBuffer::FFrameBuffer(char const *pDesc_, QSize Size_, bool ReferCurrent)
   : pDesc(pDesc_), Size(Size_)
{
   CheckGlError("FFrameBuffer::c'tor (enter)");
   if (!ReferCurrent) {
      // make a new frame buffer
      CALL_GL( glGenFramebuffers(1, &hFbo) );
      Bind(GL_FRAMEBUFFER, false);
      NeedCleanup = true;
   } else {
      // make an empty object referring to the current frame buffer.
      CALL_GL( glGetIntegerv(GL_FRAMEBUFFER_BINDING, (GLint*)&hFbo) );
      if (hFbo != 0)
         IvNotify(NOTIFY_Warning, "Default frame buffer does not have id 0. I think it should have!");
      NeedCleanup = false;
   }
   CheckGlError("FFrameBuffer::c'tor (leave)");
}

FFrameBuffer::~FFrameBuffer()
{
   if (NeedCleanup && hFbo != 0) {
      CALL_GL( glDeleteFramebuffers(1, &hFbo) );
      hFbo = 0;
   }
}

bool FFrameBuffer::Bind(GLenum Target, bool AssertCompleteness)
{
   CheckGlError("bind frame buffer (enter)");
   CALL_GL( glBindFramebuffer(Target, hFbo) );
   CALL_GL( glViewport(0, 0, width(), height()) );

   GLenum
      dwFrameBufferStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
   bool
      IsComplete = (dwFrameBufferStatus == GL_FRAMEBUFFER_COMPLETE);

   if (AssertCompleteness && IsComplete)
         IvNotify(NOTIFY_Error, QString("FFrameBuffer::Bind(%1): Bind succeeded, but frame buffer is incomplete (fb status: %2).").arg(pDesc).arg(dwFrameBufferStatus));
   CheckGlError("bind frame buffer (leave)");
   return IsComplete;
}

void FFrameBuffer::Attach(FFrameBufferAttachmentPtr p)
{
   // attach the surface to the current FBO.
   CheckGlError("attach frame buffer (enter)");
   CALL_GL( glBindFramebuffer(GL_FRAMEBUFFER, hFbo) );
   CALL_GL( glFramebufferTexture2D(GL_FRAMEBUFFER, p->GlAttachment, p->GlTarget, p->hTexture, 0) );
   Attachments.push_back(p);
}


void FFrameBuffer::BlitFrom(FFrameBuffer *pSrc, GLbitfield Layers, GLenum Filter)
{
   CheckGlError("blit frame buffer (enter)");
   CALL_GL( glBindFramebuffer(GL_READ_FRAMEBUFFER, pSrc->hFbo) );
   CALL_GL( glBindFramebuffer(GL_DRAW_FRAMEBUFFER, this->hFbo) );
   // Layers: e.g., GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT
   // Filter: e.g., GL_NEAREST
   CALL_GL( glBlitFramebuffer(0, 0, pSrc->width(), pSrc->height(), 0, 0, this->width(), this->height(), Layers, Filter) );
   CheckGlError("blit frame buffer (leave)");
}


void FFrameBuffer::Finalize()
{
   CheckGlError("FFrameBuffer::Finalize (enter)");
//    GLenum DrawBuffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
   GLenum
      DrawBuffers[16] = {0};
   uint
      nDrawBuffers = 0;
   for (uint i = 0; i < Attachments.size(); ++ i) {
      GLuint g = Attachments[i]->GlAttachment;
      // apparently I'm only supposed to register color attachments,
      // not depth or stencil.
      if (g >= GL_COLOR_ATTACHMENT0 && g <= GL_COLOR_ATTACHMENT15) {
         assert(nDrawBuffers < sizeof(DrawBuffers)/sizeof(DrawBuffers[0]));
         DrawBuffers[nDrawBuffers] = g;
         nDrawBuffers += 1;
      }
   }
   CALL_GL(glDrawBuffers(nDrawBuffers, &DrawBuffers[0]));

   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      IvNotify(NOTIFY_Error, QString("Failed to finalize deferred frame buffer '%1'.").arg(pDesc));
   CheckGlError("FFrameBuffer::Finalize (leave)");
}


FImageFilter::FImageFilter(FShaderSet &ShaderSet, FFilterType Type)
   : Filter5HV(ShaderSet)
{
   CheckGlError("FImageFilter::c'tor (enter)");
   CALL_GL( glUseProgram(ShaderSet.hProgram) );
   IdFilterStepH = ShaderSet.GetUniformLocation("FilterStepH");
   IdFilterStepV = ShaderSet.GetUniformLocation("FilterStepV");
   IdFilterWeight = ShaderSet.GetUniformLocation("FilterWeight");
   IdSourceImage = ShaderSet.GetUniformLocation("SourceImage");
   InitKernel(Type);
   CheckGlError("FImageFilter::c'tor (leave)");
}

void FImageFilter::InitKernel(FImageFilter::FFilterType Type)
{
   if (Type == FILTER_Blur) {
      GLfloat
//          BlurWeights[] = {1., -3., 9., -3., 1.},
//          BlurWeights[] = {-1., 3., 9., 3., -1.},
//          BlurWeights[] = {-1., 2, 5., 2, -1.},
         BlurWeights[5],
         w = 0;
      for (uint i = 0; i < nFilterSize; ++ i) {
         float x = (double)i - 2.;
//          float sigma = 1.3;
         float sigma = 0.7;
         BlurWeights[i] = std::exp(-sqr(x)/(2*sqr(sigma)));
      }
      for (uint i = 0; i < nFilterSize; ++ i)
         w += BlurWeights[i];
      for (uint i = 0; i < nFilterSize; ++ i)
         Kernel[i] = BlurWeights[i] / w;
   }
}


static char const *GetGlErrorString(GLenum ErrorCode) {
   switch (ErrorCode) {
      case GL_NO_ERROR: return "No error";
      case GL_INVALID_ENUM: return "Invalid enum value.";
      case GL_INVALID_VALUE: return "Invalid value (numeric argument out of range)";
      case GL_INVALID_OPERATION: return "Invalid operation (op not allowed in current state)";
      case GL_STACK_OVERFLOW: return "Command would cause stack overflow";
      case GL_STACK_UNDERFLOW: return "Command would cause stack underflow";
      case GL_OUT_OF_MEMORY: return "Out of memory (not enough memory to execute the command)";
      case GL_TABLE_TOO_LARGE: return "Table size exceed implementation's maximum allowed table size";
      default: return "(GL Error code not recognized)";
   }
}

void CheckGlError(char const *pDesc, char const *pDesc2) {
   bool
      Brrrk = false;
   size_t nErrors = 0;
   for ( ; ; ) { // hm... it's an error stack.
      GLenum
         ErrorCode = glGetError();
      if (ErrorCode == GL_NO_ERROR)
         break;
      QString
         Desc(pDesc);
      if (pDesc2)
         Desc = Desc + "/" + QString(pDesc2);
      // __debugbreak;
      // *(int*)0 = 0;
      // throw std::runtime_error("OpenGL broke :(.");

      IvNotify(NOTIFY_Error, "Problem with OpenGL: " + IvFmt("Error in %1: %2", Desc, GetGlErrorString(ErrorCode)));

      Brrrk = true;
#ifdef _MSC_VER
      break; // on windows the errors don't (always?) seem go away if you keep on eating them with glGetError().
#else
      nErrors += 1;
	  if (nErrors > 4)
	     break;
#endif
   }
//    if (Brrrk) {
// //       *(int*)0 = 0;
// //       __builtin_trap();
//       throw std::runtime_error("OpenGL broke :(.");
//    }
   (void)Brrrk;
}


// extern "C" {
//    void glClearTexImage(GLuint texture,  GLint level,  GLenum format,  GLenum type,  const void * data);
// }
// ^- not there :(.

QString LoadTextFromFile(QString FileName)
{
   QFile
      File(FileName);
   if (!File.open(QFile::ReadOnly | QFile::Text))
      return "";
   QTextStream TextStream(&File);
   return TextStream.readAll();
}


enum FGlObjectType {
   GLOBJECT_Shader,
   GLOBJECT_Program
};

GLuint GlCheckShaderStatus(GLuint hObject, FGlObjectType ObjectType, QString FileName)
{
   // check shader compilation.
   uint const nBuf = 0xffff;
   char infobuffer[nBuf];
   int infobufferlen = 0;
   GLint CompileResult = 0;
   if (ObjectType == GLOBJECT_Shader) {
      CALL_GL( glGetShaderiv(hObject, GL_COMPILE_STATUS, &CompileResult) );
      // ^- on the mac implementation glGetShaderInfoLog may return non-empty messages
      //    even if there was no error (including for things like performance warnings).
      //    So we check this here and only ask for the full result if it returns false.
      //    might otherwise confuse users
      CALL_GL( glGetShaderInfoLog(hObject, nBuf-1, &infobufferlen, infobuffer) );
   } else if (ObjectType == GLOBJECT_Program) {
      CALL_GL( glGetProgramiv(hObject, GL_LINK_STATUS, &CompileResult) );
      CALL_GL( glGetProgramInfoLog(hObject, nBuf-1, &infobufferlen, infobuffer) );
   }
   infobuffer[infobufferlen] = 0;
   QString s(&infobuffer[0]);
   if (CompileResult != GL_TRUE) {
      IvNotify(NOTIFY_Error, "Failed to compile GLSL shader: " + QString("OpenGL says for file '%1': %2").arg(FileName, s));
      return 0;
   }
   if (!s.isEmpty()) {
      IvEmit(QString("INFO: OpenGL successfully compiled '%1', but raised complaints:\n%2").arg(FileName, s));
   }

   return hObject;
}

GLuint GlLoadShader(QString FileName, QString Prefix, int Type)
{
   GLenum
      handle;
   QString
      SourceCode = LoadTextFromFile(FileName);
   if ("" == SourceCode) {
      IvNotify(NOTIFY_Error, "Failed to compile GLSL shader: " + QString("File '%1' could not be loaded").arg(FileName));
      return 0;
   }
   if (Prefix != "")
      SourceCode = Prefix + SourceCode;

   if (g_WorkAroundGlFrontFacingBug && Type == GL_FRAGMENT_SHADER) {
      SourceCode.replace("gl_FrontFacing", "true");
   }

   handle = glCreateShader(Type);
   CheckGlError("GlLoadShader","glCreateShader");
   QByteArray
      SourceCodeAsAscii = SourceCode.toLocal8Bit();
   char const *pSrc = SourceCodeAsAscii.data();

   CALL_GL( glShaderSource(handle, 1, &pSrc, 0) );
   CALL_GL( glCompileShader(handle) );

   // check shader compilation.
   return GlCheckShaderStatus(handle, GLOBJECT_Shader, FileName);
}


// see also: http://www.altdev.co/2011/06/23/improving-opengl-error-messages/
//           http://virtrev.blogspot.de/2012/12/custom-opengl-context-with-qt.html
// however, I could not get the debug context running on QT4.
// maybe I can fake it with hijacking WGL_create_context via LD_PRELOAD and putting
// in a CONTEXT_DEBUG_BIT?
//
// Or I could try it with QT5: http://qt-project.org/doc/qt-5/qopengldebuglogger.html
// It is also ugly, but maybe does the trick.

#ifdef OPENGL_DEBUG
// it's a callback for OpenGL. Needs debug context support---which QT does not
// have by default (but which can be hacked in by hijacking glXCreateContextAttribsARB
// via LD_PRELOAD and sneaking in the GLX_CONTEXT_DEBUG_BIT_ARB flag.
// See ~/dev/gl_debug_context_preload/)
void OutputGlDebugMessage(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *message, void */*userParam*/)
{
   std::string DebugMsg(message, length);

   char const *pSource = "?", *pSeverity = "?", *pType = "?";
   switch(source) {
      case GL_DEBUG_SOURCE_API_ARB:             pSource = "API"; break;
      case GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB:   pSource = "WINDOW_SYSTEM"; break;
      case GL_DEBUG_SOURCE_SHADER_COMPILER_ARB: pSource = "SHADER_COMPILER"; break;
      case GL_DEBUG_SOURCE_THIRD_PARTY_ARB:     pSource = "THIRD_PARTY"; break;
      case GL_DEBUG_SOURCE_APPLICATION_ARB:     pSource = "APPLICATION"; break;
      case GL_DEBUG_SOURCE_OTHER_ARB:           pSource = "OTHER"; break;
   }
   switch(type) {
      case GL_DEBUG_TYPE_ERROR_ARB:               pType = "ERROR"; break;
      case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB: pType = "DEPRECATED_BEHAVIOR"; break;
      case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB:  pType = "UNDEFINED_BEHAVIOR"; break;
      case GL_DEBUG_TYPE_PORTABILITY_ARB:         pType = "PORTABILITY"; break;
      case GL_DEBUG_TYPE_PERFORMANCE_ARB:         pType = "PERFORMANCE"; break;
      case GL_DEBUG_TYPE_OTHER_ARB:               pType = "OTHER"; break;
   }
   switch(severity) {
      case GL_DEBUG_SEVERITY_HIGH_ARB:   pSeverity = "HIGH";   break;
      case GL_DEBUG_SEVERITY_MEDIUM_ARB: pSeverity = "MEDIUM"; break;
      case GL_DEBUG_SEVERITY_LOW_ARB:    pSeverity = "LOW"; break;
   }
   if (severity != GL_DEBUG_SEVERITY_LOW_ARB)
   IvEmit("GL DEBUG [%2/%3]: %1 (source: %4, id: %5)", DebugMsg, pType, pSeverity, pSource, id);
//    if (severity != GL_DEBUG_SEVERITY_LOW_ARB)
//       asm("int $3");

//       __builtin_trap();
   // ^- hmpf... there is one illegal op in QGLWidget itself. can't break like this.
}
#endif
// #endif // INCLUDE_OPTIONALS



void FImageFilter::Apply(FFrameBufferPtr pDest, FFrameBufferPtr pSrc, uint Width, uint Height, FFullScreenQuad &FullScreenQuad)
{
   CheckGlError("FImageFilter::Apply (enter)");

   // this thing is only here because GetUniformLocation does not work
   // in draw mode.
   CALL_GL( glUseProgram(Filter5HV.hProgram) );

   assert(nFilterSize == sizeof(Kernel)/sizeof(Kernel[0]));
   CALL_GL( glUniform1fv(IdFilterWeight, nFilterSize, Kernel) );
   CALL_GL( glUniform1f(IdFilterStepH, 1./Width) );
   CALL_GL( glUniform1f(IdFilterStepV, 1./Height) );

   pSrc->Attachments[0]->BindAsTexture(GL_TEXTURE0);

   // set input filters to bi-linear (defaults to nearest).
   GLint MinFilter, MagFilter;
   CALL_GL( glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, &MinFilter) );
   CALL_GL( glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, &MagFilter) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) );

   CALL_GL( glUniform1i(IdSourceImage, 0) );
   pDest->Bind(GL_DRAW_FRAMEBUFFER);

   CALL_GL( glDisable(GL_DEPTH_TEST) ); // FIXME: move to FFullScreenQuad::Draw()
   FullScreenQuad.Draw(Filter5HV.hProgram);
   CALL_GL( glEnable(GL_DEPTH_TEST) ); // FIXME: move to FFullScreenQuad::Draw()

   // reset filters.
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, MinFilter) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, MagFilter) );

   // unbind the frame buffer texture.
   pSrc->Attachments[0]->UnbindAsTexture(GL_TEXTURE0);
   CheckGlError("FImageFilter::Apply (leave)");
}


FShaderSet::FShaderSet()
{
   ZeroOutHandles();
}

FShaderSet::FShaderSet(QString const &sVertexFile, QString const &sPixelFile, QString const &BasePath)
{
   ZeroOutHandles();
   Create(sVertexFile, sPixelFile, BasePath);
}

void FShaderSet::ZeroOutHandles()
{
   hProgram = 0;
   hVertexShader = 0;
   hPixelShader = 0;
   hPixelCommon = 0;
}

void FShaderSet::Create(QString const &sVertexFile, QString const &sPixelFile, QString const &sBasePath)
{
   QDir
      BasePath(sBasePath);
   CheckGlError("FShaderSet::Create", "(enter)");
   if (hProgram != 0)
      Destroy();
   ZeroOutHandles();
//    static std::string BasePath = GetExePath() + "/";
//    QString BasePath = ":/"; // <- embedded resource root.

   hProgram = glCreateProgram();
   CheckGlError("FShaderSet::Create", "CreateProgram");
#ifndef __APPLE__
   QString
      ShaderPrefix = "#version 400\n\n";
//       ShaderPrefix = "#version 300 es\nprecision mediump float;\n\n";
      // ^- we add this in software now... to make it work with Intel linux
      //    drivers, which can do a number of OpenGL 4.0 things, but only really
      //    accepts "3.0 es" as #version.
   if (g_GlVersion < 40) {
      ShaderPrefix = "#version 300 es\nprecision mediump float;\n\n";
   }
#else
   // Apple. This and #version 150 were the ONLY ones it accepted.
   // Nothing else (no 140, no 130, no 300, no 330, no 300 es, ...).
   // ...In a OpenGL 3.2 core context.
   QString
      ShaderPrefix = "#version 400\n\n";
#endif
   hVertexShader = GlLoadShader(BasePath.filePath(sVertexFile), ShaderPrefix, GL_VERTEX_SHADER);
   hPixelShader = GlLoadShader(BasePath.filePath(sPixelFile), ShaderPrefix, GL_FRAGMENT_SHADER);
   hPixelCommon = GlLoadShader(BasePath.filePath(QString("pixel_common.glsl")), ShaderPrefix, GL_FRAGMENT_SHADER);

   CALL_GL( glAttachShader(hProgram, hVertexShader) );
   CALL_GL( glAttachShader(hProgram, hPixelShader) );
   CALL_GL( glAttachShader(hProgram, hPixelCommon) );
   CALL_GL( glLinkProgram(hProgram) );

   GlCheckShaderStatus(hProgram, GLOBJECT_Program, QString("LINK [%1 // %2]").arg(sVertexFile, sPixelFile));

   // query and store some uniform locations which most shaders have
   for (uint i = 0; i < sizeof(hMatrix)/sizeof(hMatrix[0]); ++ i)
      hMatrix[i] = -1;
   hMatrix[TRANSFORM_ModelView] = glGetUniformLocation(hProgram, "mView");
   hMatrix[TRANSFORM_Projection] = glGetUniformLocation(hProgram, "mProj");
   hMatrix[TRANSFORM_Normals] = glGetUniformLocation(hProgram, "mNorm");
   hDiffuseColor = glGetUniformLocation(hProgram, "DiffuseColor");
   hObjectId = glGetUniformLocation(hProgram, "ObjectId");

   CheckGlError("FShaderSet::Create", "(exit)");
}


void FShaderSet::Destroy()
{
   CALL_GL( glUseProgram(0) );
   CALL_GL( glDeleteProgram(hProgram) );
   CALL_GL( glDeleteShader(hVertexShader) );
   CALL_GL( glDeleteShader(hPixelShader) );
   CALL_GL( glDeleteShader(hPixelCommon) );
   // note: deleting 0 objects is fine. Result is ignored.

   ZeroOutHandles();
}

FShaderSet::~FShaderSet()
{
   Destroy();
}

void SetGlNormalTransformationMatrix(GLint hTrafoNorm, FMat4f const &mView)
{
   // make and set the transformation for normal vectors. may be non-trivial if the
   // main transformation is not unitary (and we use non-unitary trafos to
   // place atoms and bonds).
   FMat3f
      mNorm1;
   for (uint i = 0; i != 3; ++ i)
      for (uint j = 0; j != 3; ++ j)
         mNorm1(i,j) = mView(i,j);
   mNorm1 = vmath::transpose(vmath::inverse(mNorm1));
   CALL_GL( glUniformMatrix3fv(hTrafoNorm, 1, GL_FALSE, &mNorm1(0,0)) );
}

void FGlMatrixStack::Actualize(FShaderSet &ShaderSet)
{
   GLint
      hMatrix = ShaderSet.hMatrix[this->DefaultRole];
   if (hMatrix == -1)
      throw std::runtime_error("FGlMatrixStack::Actualize: Current shader set does not have a matrix of the current role.");

   CALL_GL( glUniformMatrix4fv(hMatrix, 1, GL_FALSE, &Top()(0,0)) );
   // ^- false: do not transpose (openGL expects column major order, i.e., row stride = 1, column stride = nrows).
   //    that is what I use, too.
   CheckGlError("FGlMatrixStack::Actualize", "set matrices.");
   if (DefaultRole == TRANSFORM_ModelView)
      SetGlNormalTransformationMatrix(ShaderSet.hMatrix[TRANSFORM_Normals], Top());
}

void FGlMatrixStack::Translate(FVec3f const &v)
{
//    FMat4f
//       Dummy = vmath::identity4<float>();
//    for (uint i = 0; i < 3; ++ i)
//       Dummy(i,3) = v[i];
//    Multiply(Dummy);
   float
      fRowDisp[3] = {0, 0, 0};
   for (uint iRow = 0; iRow != 3; ++ iRow)
      for (uint iComp = 0; iComp != 3; ++ iComp)
         fRowDisp[iRow] += v[iComp] * Top()(iRow, iComp);
   for (uint iRow = 0; iRow != 3; ++ iRow)
      Top()(iRow,3) += fRowDisp[iRow];
}

void FGlMatrixStack::Scale(FVec3f const &s)
{
   // Top = Scale * Top
   // where
   //           ( s[0]              )
   //   Scale = (      s[1]         )
   //           (           s[2]    )
   //           (                1. )
   for (uint iDir = 0; iDir < 3; ++ iDir)
      for (uint iComp = 0; iComp < 4; ++ iComp)
         Top()(iComp,iDir) *= s[iDir];
}

void FGlMatrixStack::SetIdentity()
{
   memset(&Top()(0,0), 0, sizeof(FMat4f));
   for (uint i = 0; i < 4; ++ i)
      Top()(i,i) = 1.f;
}

void FGlMatrixStack::Clear()
{
   FBase::resize(1);
   SetIdentity();
}

void FGlMatrixStack::Set(FMat4f const &m)
{
   Top() = m;
}

void FGlMatrixStack::Multiply(FMat4f const &m)
{
   Top() = m * Top();
}




#include "fn_LiberationSans.h"

template<class FVertex>
TGlMesh<FVertex>::TGlMesh(TIndexedTriangleList<FVertex> const &List, GLenum UsageType_)
{
   m_UsageType = UsageType_;
   Vertices = List.Vertices;
   Triangles = List.Triangles;
   CreateGlBindings();
}

template<class FVertex>
TGlMesh<FVertex>::TGlMesh(GLenum UsageType_)
{
   m_UsageType = UsageType_;
}


template<class FVertex>
TGlMesh<FVertex>::~TGlMesh()
{
   DestroyGlBindings();
}

template<class FVertex>
void TGlMesh<FVertex>::Invalidate()
{
   DestroyGlBindings();
   CreateGlBindings(true); // rebuild interfaces and upload data
}


void FBaseVertex::AssignVertexAttributes()
{
   FBaseVertex v;

   if (1) {
      // GL4? required?
      for (uint i = 0; i < 4; ++ i ) {
         CALL_GL( glEnableVertexAttribArray(i) );
      }
      CALL_GL( glDisableVertexAttribArray(1) );

      // position
      CALL_GL( glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(v), (void*)((char*)&v.vPos - (char*)&v)) );
      // normal
      CALL_GL( glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(v), (void*)((char*)&v.vNorm - (char*)&v)) );
      // color (normalize: convert 0..255 -> 0...1)
      CALL_GL( glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(v), (void*)((char*)&v.dwColor - (char*)&v)) );
   } else {
      // this worked in GL3
      glVertexPointer(3, GL_FLOAT, sizeof(v), (void*)((char*)&v.vPos - (char*)&v));
      glNormalPointer(GL_FLOAT, sizeof(v), (void*)((char*)&v.vNorm - (char*)&v));
   //    glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(v), (void*)((char*)&v.vNorm - (char*)&v));
      glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(v), (void*)((char*)&v.dwColor - (char*)&v));
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glEnableClientState(GL_COLOR_ARRAY);
   }
}

void FTextVertex::AssignVertexAttributes()
{
   FTextVertex v;

   // GL4? required?
   for (uint i = 0; i < 4; ++ i ) {
      CALL_GL( glEnableVertexAttribArray(i) );
   }
//    glDisableVertexAttribArray(2);

   // position (3d)
   CALL_GL( glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(v), (void*)((char*)&v.vPos - (char*)&v)) );
   // position (2d displacement)... not sure if stupid to do it like this.
   CALL_GL( glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(v), (void*)((char*)&v.vCoord2d - (char*)&v)) );
   // texture coordinate
   CALL_GL( glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(v), (void*)((char*)&v.vTex - (char*)&v)) );
   // color (normalize: convert 0..255 -> 0...1)
   CALL_GL( glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(v), (void*)((char*)&v.dwColor - (char*)&v)) );
}


template<class FVertex>
void TGlMesh<FVertex>::CreateGlBindings(bool UploadData)
{
   CALL_GL( glGenVertexArrays(1, &m_hVertexArrayObject) ); // <- that's a descriptor of vertex data.
   CALL_GL( glBindVertexArray(m_hVertexArrayObject) );
   // ^- defines layout of vertex buffers, effectively. also stores
   //    which buffers were bound to what.
   CheckGlError("FGlMesh::CreateGlBindings", "create vertex descriptor");

   FVertex
      *pVertices = 0;
   vec3ui
      *pTriangles = 0;
   if (UploadData) {
      m_MaxTriangles = Triangles.size();
      m_MaxVertices = Vertices.size();
      pTriangles = &Triangles[0];
      pVertices = &Vertices[0];
   }

   CALL_GL( glGenBuffers(1, &m_hVertexBuffer) );
   CALL_GL( glBindBuffer(GL_ARRAY_BUFFER, m_hVertexBuffer) );
   CALL_GL( glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices[0])*m_MaxVertices, pVertices, m_UsageType) );
   CheckGlError("FGlMesh::CreateGlBindings", "create vertex buffer");

   // setup the vertex descriptor (didn't work on some machines if done before the GL_ARRAY_BUFFER stuff)
   // note: only allowed if a "vertex array object" is bound (glBindVertexArray)
   FVertex::AssignVertexAttributes();
   CheckGlError("FGlMesh::CreateGlBindings", "setup vertex descriptor");

   CALL_GL( glBindVertexArray(0) );
// ^- hm.. why was this here? doesn't this un-associate the vertex descriptor with the vertex array object?
   CALL_GL( glGenBuffers(1, &m_hIndexBuffer) );
   CALL_GL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_hIndexBuffer) );
   CALL_GL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Triangles[0]) * m_MaxTriangles, pTriangles, m_UsageType) );

   CheckGlError("FGlMesh::CreateGlBindings", "create index buffer");
}

template<class FVertex>
void TGlMesh<FVertex>::UploadData()
{
   if (Triangles.size() > m_MaxTriangles || Vertices.size() > m_MaxVertices) {
      IvNotify(NOTIFY_Error, "Ran out of buffer space in TGlMesh<FVertex>::UploadData().");
      return;
      // alternatively I could invalidate and re-upload. But I think this cannot
      // be done in all render states, can it?
   }
   CALL_GL( glBindVertexArray(m_hVertexArrayObject) );
   // ^- I think this binds the other arrays automagically (hm.. or maybe it doesn't. crashes if I remove the
   //    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,...) in the Draw() routine.
   CALL_GL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_hIndexBuffer) );
   CALL_GL( glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(Triangles[0]) * Triangles.size(), &Triangles[0]) );
   CALL_GL( glBindBuffer(GL_ARRAY_BUFFER, m_hVertexBuffer) );
   CALL_GL( glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertices[0])*Vertices.size(), &Vertices[0]) );
}

template<class FVertex>
void TGlMesh<FVertex>::DestroyGlBindings()
{
   CheckGlError("TGlMesh<FVertex>::DestroyGlBindings() (enter)");
   CALL_GL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
   CALL_GL( glDeleteBuffers(1, &m_hIndexBuffer) );
   m_hIndexBuffer = 0;

   CALL_GL( glBindBuffer(GL_ARRAY_BUFFER, 0) );
   CALL_GL( glDeleteBuffers(1, &m_hVertexBuffer) );
   m_hVertexBuffer = 0;

//    for (uint i = 0; i < 4; ++ i )
//       glDisableVertexAttribArray(i);
// ^- may not be allowed unless the vertex array object is bound. also, doesn't actually do anything here.

   CALL_GL( glBindVertexArray(0) );
   CALL_GL( glDeleteVertexArrays(1, &m_hVertexArrayObject) );
   m_hVertexArrayObject = 0;
   CheckGlError("TGlMesh<FVertex>::DestroyGlBindings() (leave)");
}

template<class FVertex>
void TGlMesh<FVertex>::Draw()
{
   CheckGlError("FGlMesh::Draw", "(enter)");
//    glBindVertexArray(hVertexArrayObject);     // defines layout of vertex buffers, effectively.
//    glBindBuffer(GL_ARRAY_BUFFER, hVertexBuffer);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, hIndexBuffer);
//    for (uint i = 0; i < 3; ++ i )
//       glEnableVertexAttribArray(i);

   CALL_GL( glBindVertexArray(m_hVertexArrayObject) );
   CALL_GL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_hIndexBuffer) );

//    std::cout << "Triangles.size: " << Triangles.size() << std::endl;
   CALL_GL( glDrawElements(GL_TRIANGLES, 3*Triangles.size(), GL_UNSIGNED_INT, (GLvoid*)0) );
   // ^- why * 3? Do I need to specify the number of indices instead
   //    of the number of primitives? The way the docs are phrased suggests 'no'.

   CALL_GL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
   CALL_GL( glBindVertexArray(0) );
//    glBindBuffer(GL_ARRAY_BUFFER, 0);
   CheckGlError("FGlMesh::Draw", "(leave)");
}

template struct TGlMesh<FBaseVertex>;


FGlTextBuffer::FGlTextBuffer(size_t MaxCharacters)
//    : FBase(GL_STREAM_DRAW)
   : FBase(GL_DYNAMIC_DRAW)
{
   m_pFont = &embedded_fonts::efont_sans;
   m_Dirty = true;
   m_MaxVertices = 4 * MaxCharacters;
   m_MaxTriangles = 2 * MaxCharacters;
   Triangles.reserve(m_MaxTriangles);
   Vertices.reserve(m_MaxVertices);
   CreateGlBindings(false); // make buffers, but don't upload any data.
}

FGlTextBuffer::~FGlTextBuffer()
{
}

void FGlTextBuffer::Clear()
{
   Triangles.clear();
   Vertices.clear();
   m_Dirty = true;
}


void FGlTextBuffer::Print(vec3f vPos, float Scale, uint32_t dwColor, wchar_t const *pText, uint Flags)
{
   size_t
      iFirstVertex = Vertices.size();
   wchar_t
      // last emitted character (for building kerning pairs)
      LastChar = 0;
   vec3f
      // position of text. x,y in units of pixels. Will be fixed up and moved
      // to vPos later.
      vCur(0.f,0.f,0.f);
//    float
//       fMaxHeight = 0.;
   size_t
      nChars = wcslen(pText);
   if (2*nChars + Triangles.size() > m_MaxTriangles) {
      IvNotify(NOTIFY_Warning, "Text buffer too small. Some text will not be rendered (FGlTextBuffer::Print).");
      return;
   }
   for (size_t i = 0; i < nChars; ++ i)
   {
      embedded_fonts::texture_glyph_t
         *pGlyph = 0;
      for (size_t j = 0; j < m_pFont->glyphs_count; ++ j)
         if (m_pFont->glyphs[j].charcode == pText[i]) {
            pGlyph = &m_pFont->glyphs[j];
            break;
         }
      if (!pGlyph)
         continue;

      float
         kerning = 0.f;
      for (size_t ikern = 0; ikern < pGlyph->kerning_count; ++ ikern)
         if (pGlyph->kerning[ikern].charcode == LastChar)
            kerning = pGlyph->kerning[ikern].kerning;
      vCur[0] += kerning;

      float
         x = vCur[0] + (float)pGlyph->offset_x,
         y = vCur[1] + (float)pGlyph->offset_y,
         w = (float)pGlyph->width,
         h = (float)pGlyph->height;
      size_t iVert = Vertices.size();
//       Vertices.push_back(FTextVertex(vec3f(x,  y,  vPos.z), vec2f(pGlyph->s0, pGlyph->t0), dwColor));
//       Vertices.push_back(FTextVertex(vec3f(x,  y-h,vPos.z), vec2f(pGlyph->s0, pGlyph->t1), dwColor));
//       Vertices.push_back(FTextVertex(vec3f(x+w,y-h,vPos.z), vec2f(pGlyph->s1, pGlyph->t1), dwColor));
//       Vertices.push_back(FTextVertex(vec3f(x+w,y,  vPos.z), vec2f(pGlyph->s1, pGlyph->t0), dwColor));

      Vertices.push_back(FTextVertex(vPos, vec2f(x,   y),   vec2f(pGlyph->s0, pGlyph->t0), dwColor));
      Vertices.push_back(FTextVertex(vPos, vec2f(x,   y-h), vec2f(pGlyph->s0, pGlyph->t1), dwColor));
      Vertices.push_back(FTextVertex(vPos, vec2f(x+w, y-h), vec2f(pGlyph->s1, pGlyph->t1), dwColor));
      Vertices.push_back(FTextVertex(vPos, vec2f(x+w, y),   vec2f(pGlyph->s1, pGlyph->t0), dwColor));

//       Vertices.push_back(FTextVertex(vPos, vec2f(x,   y),   vec2f(0., 0.), dwColor));
//       Vertices.push_back(FTextVertex(vPos, vec2f(x,   y-h), vec2f(0., 1.), dwColor));
//       Vertices.push_back(FTextVertex(vPos, vec2f(x+w, y-h), vec2f(1., 1.), dwColor));
//       Vertices.push_back(FTextVertex(vPos, vec2f(x+w, y),   vec2f(1., 0.), dwColor));
      Triangles.push_back(vec3ui(iVert+0, iVert+1, iVert+2));
      Triangles.push_back(vec3ui(iVert+0, iVert+2, iVert+3));
//       fMaxHeight = std::max(fMaxHeight, h);

      vCur[0] += pGlyph->advance_x;
      vCur[1] += pGlyph->advance_y;

      LastChar = pText[i];
   }

   // apply scaling and transformation to x and y positions of glyphs.
   vec2f
      vOff(0., 0.);
   if ((Flags & HALIGN_Mask) == HALIGN_Center)
      vOff[0] -= vCur[0]/2;
   if ((Flags & HALIGN_Mask) == HALIGN_Right)
      vOff[0] -= vCur[0];
   if ((Flags & VALIGN_Mask) == VALIGN_Center)
      vOff[1] -= (m_pFont->ascender + m_pFont->descender)/2;
   if ((Flags & VALIGN_Mask) == VALIGN_Top)
      vOff[1] -= (m_pFont->ascender + m_pFont->descender);

   for (size_t iVert = iFirstVertex; iVert != Vertices.size(); ++ iVert) {
//       vec3f &v = Vertices[iVert].vPos;
//       v = (v + vOff)*(Scale/m_pFont->size);
      vec2f &v = Vertices[iVert].vCoord2d;
      v = (v + vOff)*(Scale/m_pFont->size);
//       std::cout << format("text v2f: (%8.3f, %8.3f)\n") % v[0] % v[1];
   }

   m_Dirty = true;
}

void FGlTextBuffer::DrawText1(FGlMatrixStack *pView, FGlMatrixStack *pProj, FUniformProxy *pUniformProxy)
{
   if (m_Dirty) {
      UploadData();
      m_Dirty = false;
   }
   CALL_GL( glUseProgram(pShader->hProgram) );
   CALL_GL( glUniform1i(pShader->GetUniformLocation("FontTexture" ), 0) );
   CALL_GL( glActiveTexture(GL_TEXTURE0) );
   CALL_GL( glBindTexture(GL_TEXTURE_2D, m_hFontTexture) );

   if (pView)
      pView->Actualize(*pShader);
   if (pProj)
      pProj->Actualize(*pShader);

   if (pUniformProxy)
      pUniformProxy->SetUniforms(*pShader);

//    {
//       FGlMatrixStack mIdProj(TRANSFORM_Projection);
//       FGlMatrixStack mIdView(TRANSFORM_ModelView);
//       mIdProj.SetIdentity();
//       mIdView.SetIdentity();
//       CALL_GL( glUniformMatrix4fv(pShader->hMatrix[TRANSFORM_ModelView], 1, GL_FALSE, &mIdView.Top()(0,0)) );
//       CALL_GL( glUniformMatrix4fv(pShader->hMatrix[TRANSFORM_Projection], 1, GL_FALSE, &mIdProj.Top()(0,0)) );
//    }
   CALL_GL( glHint(GL_FRAGMENT_SHADER_DERIVATIVE_HINT, GL_NICEST) );

//    UploadData(); // FIXME: should this be here?
//    std::cout << format("drawing %i triangles and %i vertices via text buffer.\n") % Triangles.size() % Vertices.size();
   FBase::Draw();

   CALL_GL( glHint(GL_FRAGMENT_SHADER_DERIVATIVE_HINT, GL_DONT_CARE) );
}

void FGlTextBuffer::Invalidate()
{
   FBase::Invalidate();
   m_Dirty = false;
}

void FGlTextBuffer::CreateGlBindings(bool UploadData)
{
   FBase::CreateGlBindings(UploadData);
   // create texture object.

   CALL_GL( glGenTextures(1, &m_hFontTexture) );
   CALL_GL( glActiveTexture(GL_TEXTURE0) );
   CALL_GL( glBindTexture(GL_TEXTURE_2D, m_hFontTexture) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR) );
//    CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) );
   CALL_GL( glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 9) );

   CALL_GL( glTexImage2D( GL_TEXTURE_2D, 0, GL_R8, m_pFont->tex_width, m_pFont->tex_height,
                  0,  GL_RED, GL_UNSIGNED_BYTE, &m_pFont->tex_data[0] ) );
//    CALL_GL( glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST) );
// ^- apparently not allowed in the core profile. There doesn't seem to be another way
//    of controlling this, either.
   CALL_GL( glGenerateMipmap(GL_TEXTURE_2D) );


   pShader = new FShaderSet("vertex_text.glsl", "pixel_text.glsl");
}

void FGlTextBuffer::DestroyGlBindings()
{
   pShader = 0; // should destroy objects automatically.
   FBase::DestroyGlBindings();
}

void FGlTextBuffer::UploadData()
{
   FBase::UploadData();
}
