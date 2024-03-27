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

#ifndef IV_GL_H
#define IV_GL_H

#include "GL/glew.h"
#define QT_NO_OPENGL_ES_2
#define GL_GLEXT_PROTOTYPES

#include "Iv.h"

#include <QRect>
#include <vector>
#include <cmath>

#include "IvMesh.h" // for FBaseVertex and the indexed triangle lists.

#ifdef __GNUC__
   #define IV_FUNCTION_NAME __PRETTY_FUNCTION__
#else
   #define IV_FUNCTION_NAME __FUNCTION__
   // ^- VC doesn't have it.
#endif

// note: OPENGL_DEBUG needs export LD_PRELOAD=/home/cgk/dev/gl_debug_context_preload/replace_context.so
//       in order to hack in the debug context creation flags.
// #define OPENGL_DEBUG
void CheckGlError(char const *pDesc, char const *pDesc2=0);
// well... it is less ugly than explicit call/check pairs for *every single* GL function call.
#ifdef _MSC_VER
   #define IV_LINE_STR2(x) #x
   #define IV_LINE_STR(x) IV_LINE_STR2(x)
   #define FUNCTION_LINE_NAME IV_FUNCTION_NAME " (" __FILE__ ":" IV_LINE_STR(__LINE__) ")"
   #define CALL_GL(X) { X; CheckGlError(FUNCTION_LINE_NAME, #X); }
#else
   #define CALL_GL(X) { X; CheckGlError(IV_FUNCTION_NAME, #X); }
#endif

#ifdef OPENGL_DEBUG
void OutputGlDebugMessage(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *message, void */*userParam*/);
#endif

namespace ct {
void PrintMatrixGen( std::ostream &out, float const *pData,
        uint nRows, uint nRowStride, uint nCols, uint nColStride,
        std::string const &Name );
}

extern int g_GlVersion;

typedef TVector3<float>
   FVec3f;
typedef mat4f
   FMat4f;
typedef mat3f
   FMat3f;

static const float fEpsilon = 1e-6f;
static const float dEpsilon = 1e-12;

inline bool IsEq(float a, float b, float Epsilon = fEpsilon) {
   return (std::abs(a-b) < Epsilon);
}

inline bool IsEq(double a, double b, double Epsilon = dEpsilon) {
   return (std::abs(a-b) < Epsilon);
}

// if set, gl_FrontFacing will be text-replaced by "true" and viewing back faces will be disabled (this one is for Mac)
extern bool g_WorkAroundGlFrontFacingBug;
// if set, dest alpha will be cleared to 1.0 instead of 0.0, and copying of alpha channels in pictures will be disabled
// (on Mac the main FBO alpha buffer seems to be used for blending the QGlWidget with the background. I put in a
// setAttribute(Qt::WA_TranslucentBackground, false);, but I am not sure if it will work)
extern bool g_WorkAroundAlphaCompositing;


enum FTransformType {
   TRANSFORM_ModelView = 0,
   TRANSFORM_Projection = 1,
   TRANSFORM_Normals = 2,
   TRANSFORM_Count = 3
};

struct FShaderSet;

struct FGlMatrixStack : protected ct::TArray<FMat4f>
{
   explicit FGlMatrixStack(FTransformType DefaultRole_) : DefaultRole(DefaultRole_) { Clear(); }

   FTransformType
      DefaultRole;
   typedef ct::TArray<FMat4f>
      FBase;
   FMat4f &Top() { return FBase::back(); };
   FMat4f const &Top() const { return FBase::back(); };

   void Push() { FBase::push_back(Top()); }
   void Pop() { FBase::pop_back(); }

   void Translate(FVec3f const &v);
   void Scale(FVec3f const &s);
   // set current level to Identity transform
   void SetIdentity();
   // clear all levels (except first) and set first level to Identity transform.
   void Clear();
   void Set(FMat4f const &m);
   void Multiply(FMat4f const &m);

   // export the current matrix state to the given shader
   void Actualize(FShaderSet &ShaderSet);
};

// assemble a 4x4 matrix by stacking the given row vectors horizontally:
//        ( r0x r1x r2x r3z )
//  Out = ( r0y r1y r2y r3y )
//        ( r0z r1z r2z r3z )
//        (  0   0   0   1  )
FMat4f hStack(FVec3f const &Row0, FVec3f const &Row1, FVec3f const &Row2, FVec3f const &Row3);


struct FShaderSet : public ct::FIntrusivePtrDest
{
   GLenum
      hProgram,
      hVertexShader,
      hPixelShader,
      hPixelCommon;
   GLint
      // uniform locations of model/view, projection matrix, and diffuse color.
      // All of those might not be present and will be set to -1 in this case.
      // These are here because most shaders have them (but not all).
      hMatrix[TRANSFORM_Count],
      hDiffuseColor,
      hObjectId;
   void Create(QString const &sVertexFile, QString const &sPixelFile, QString const &BasePath=":/shader");
   void Destroy();
   void ZeroOutHandles();
   // todo: move stuff for getting uniform locations here.

   FShaderSet();
   FShaderSet(QString const &sVertexFile, QString const &sPixelFile, QString const &BasePath=":/shader");
   ~FShaderSet();


   GLint GetUniformLocation(char const *pName) {
      CheckGlError("GetUniformLocation (enter)", pName);
      GLint res = glGetUniformLocation(this->hProgram, pName);
      CheckGlError("GetUniformLocation", pName);
      return res;
   }
private:
   void operator = (FShaderSet const &); // not implemented.
   FShaderSet(FShaderSet const &); // not implemented.
};

typedef ct::TIntrusivePtr<FShaderSet>
   FShaderSetPtr;


// callback for setting shader uniforms.
struct FUniformProxy
{
   virtual void SetUniforms(FShaderSet &ShaderSet) = 0;
   ~FUniformProxy();
};


struct FFullScreenQuad
{
   GLenum
      IdVertexBuffer,
      hVertexDesc;
   void Init();
   ~FFullScreenQuad();

   void Draw(GLenum Program);
};


struct FFrameBufferAttachment : public ct::FIntrusivePtrDest
{
   GLuint
      hTexture;
   GLenum
      GlTarget,
      GlAttachment;
   uint
      nSamples;
   QSize
      Size;
   char const
      *pDesc; // for error reporting and debugging.
   // construct a texture and add it as attachment to the current frame buffer.
   // 'target' is either GL_TEXTURE_2D (in case of nSamples==0) or GL_TEXTURE_2D_MULTISAMPLE
   // (otherwise).
   // Attachment: e.g., GL_COLOR_ATTACHMENT1
   //
   // What is the difference between InternalFormat, Format, and Type? No idea.
   // Here are some examples which work in this combination:
   //
   // InternalFormat: e.g., GL_R32F.
   // Format: e.g., GL_RED.
   // Type: e.g., GL_FLOAT.
   //
   // note: this attaches to the currently bound frame buffer, whatever that might be!
   FFrameBufferAttachment(char const *pDesc, QSize Size, GLenum Attachment, uint nSamples_, GLuint InternalFormat, GLuint Format, GLenum Type);
   ~FFrameBufferAttachment();

   uint width() const { return Size.width(); };
   uint height() const { return Size.height(); };

   // GlTexture: e.g., GL_TEXTURE0
   void BindAsTexture(GLenum GlTexture);
   void UnbindAsTexture(GLenum GlTexture);
};
typedef ct::TIntrusivePtr<FFrameBufferAttachment>
   FFrameBufferAttachmentPtr;


struct FFrameBuffer : public ct::FIntrusivePtrDest
{
   typedef std::vector<FFrameBufferAttachmentPtr>
      FAttachmentList;
   FAttachmentList
      Attachments;
   GLuint
      hFbo;
   char const
      *pDesc;
   bool
      NeedCleanup;
   QSize
      Size;

   // this constructs and binds a new frame buffer.
   FFrameBuffer(char const *pDesc_, QSize Size_, bool ReferCurrent = false);
   ~FFrameBuffer();

   bool Bind(GLenum Target = GL_FRAMEBUFFER, bool AssertCompleteness = false);
   void Attach(FFrameBufferAttachmentPtr Attachment);

   void BlitFrom(FFrameBuffer *pSrc, GLbitfield Layers, GLenum Filter);

   // collect all attachments and activate them as draw buffers.
   void Finalize();

//    uint width() const { return Attachments[0]->Size.width(); }
//    uint height() const { return Attachments[0]->Size.height(); }
   uint width() const { return Size.width(); }
   uint height() const { return Size.height(); }
};
typedef ct::TIntrusivePtr<FFrameBuffer>
   FFrameBufferPtr;


struct FImageFilter : public ct::FIntrusivePtrDest
{
   enum FFilterType {
      FILTER_Blur,
      FILTER_Other
   };

   static const uint
      nFilterSize = 5;
   GLint
      IdFilterStepH,
      IdFilterStepV,
      IdFilterWeight,
      IdSourceImage;
   FShaderSet
      &Filter5HV;
   GLfloat
      Kernel[5],
      Step;

   FImageFilter(FShaderSet &ShaderSet, FFilterType Type);
   void InitKernel(FFilterType Type);
   void Apply(FFrameBufferPtr pDest, FFrameBufferPtr pSrc, uint Width, uint Height, FFullScreenQuad &FullScreenQuad);
};
typedef ct::TIntrusivePtr<FImageFilter>
   FImageFilterPtr;

void CheckGlError(char const *pDesc, char const *pDesc2);


QString LoadTextFromFile(QString FileName);
GLenum GlLoadShader(QString FileName, QString Prefix, int Type);


namespace embedded_fonts {
   struct texture_font_t;
}


// vertex format used for drawing text (via distance-mapped textures)
struct FTextVertex {
   vec3f vPos; // center position
   vec2f vCoord2d; // xy displacement in view space.
   vec2f vTex;
   uint32_t dwColor;
   static void AssignVertexAttributes();
   FTextVertex() {}
   FTextVertex(vec3f const &vPos_, vec2f vCoord2d_, vec2f const &vTex_, uint32_t dwColor_) : vPos(vPos_), vCoord2d(vCoord2d_), vTex(vTex_), dwColor(dwColor_) {}
};

template<class FVertex>
struct TGlMesh : public ct::FIntrusivePtrDest
{
   ct::TArray<FVertex>
      Vertices;
   ct::TArray<vec3ui>
      Triangles;

   ~TGlMesh();
//    FGlMesh(FBaseVertex *pVertices_, unsigned nVertices, FIndexType *pIndices_, unsigned nIndices);
   explicit TGlMesh(TIndexedTriangleList<FVertex> const &List, GLenum UsageType_=GL_STATIC_DRAW);
   virtual void Draw();
   virtual void Invalidate();
protected:
   explicit TGlMesh(GLenum UsageType_=GL_STATIC_DRAW);

   GLuint
      m_hVertexArrayObject,
      m_hIndexBuffer,
      m_hVertexBuffer;
   GLenum
      m_UsageType;

   virtual void CreateGlBindings(bool UploadData=true);
   virtual void DestroyGlBindings();
   virtual void UploadData(); // copy data from this->Triangles and this->Vertices into the GL buffer objects.
   size_t
      // maximum number of vertices and triangles for which we have space
      // in the GL objects.
      m_MaxVertices,
      m_MaxTriangles;
private:
   void operator = (TGlMesh const &); // not implemented.
   TGlMesh(TGlMesh const &); // not implemented.
};
typedef TGlMesh<FBaseVertex>
   FGlMesh;

typedef ct::TIntrusivePtr<FGlMesh>
   FGlMeshPtr;


struct FGlTextBuffer : public TGlMesh<FTextVertex>
{
   typedef TGlMesh<FTextVertex>
      FBase;
   enum FPrintFlags{
      HALIGN_Left = 0x00,
      HALIGN_Center = 0x01, // center text horizontally
      HALIGN_Right = 0x02,
      HALIGN_Mask = 0x03,
      VALIGN_Top = 0x00,    // given position to top of text
      VALIGN_Center = 0x04, // center text vertically
      VALIGN_Baseline = 0x08, // align baseline to given coordinate
      VALIGN_Mask = 0x0c
   };

   explicit FGlTextBuffer(size_t MaxCharacters);
   ~FGlTextBuffer();
   void Clear();
   void Print(vec3f Pos, float Scale, uint32_t dwColor, wchar_t const *pText, uint Flags);
   // note:
   // - it might be a good idea to keep vectors as member variables which denote the positive x and y directions.
   // - flag affecting the text appearance (e.g. color, boldness) could also be stored as state.
   //   note that with my distance fields I can make letters bolder or less bold comparatively easily,
   //   and I can apply outline colors and foll colors separately.

   void CreateGlBindings(bool UploadData=true); // override
   void DestroyGlBindings(); // override
   void UploadData(); // override

   void DrawText1(FGlMatrixStack *pView, FGlMatrixStack *pProj, FUniformProxy *pUniformProxy=0);
   // ^- beware of funky #define of DrawText to DrawTextW on win32.
   void Invalidate(); // override
   FShaderSetPtr
      pShader;
protected:
   bool
      // set if this->Triangles or this->Vertices are out of sync with the
      // GL objects.
      m_Dirty;

   GLuint
      m_hFontTexture;
   embedded_fonts::texture_font_t
      *m_pFont;
private:
   void operator = (FGlTextBuffer const &); // not implemented.
   FGlTextBuffer(FGlTextBuffer const &); // not implemented.
};

typedef ct::TIntrusivePtr<FGlTextBuffer>
   FGlTextBufferPtr;




#endif // IV_GL_H
