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

// main header of IboView... define some stuff used in many places.
#ifndef CT_IV_H
#define CT_IV_H

// #ifdef _WIN32_
//   // otherwise lots of strange errors in qopengles2ext
//   #define WIN32_LEAN_AND_MEAN
//   #include <windows.h>
// #endif
// (update: that seemed to be a problem with my QT version, which emulated
//  GL2ES on top of DirectX, as it seems.)

#if defined(_MSC_VER) && defined (_M_X64)
	#define NEED_EXTRA_SIZE_T_TYPES
	// ^- in windows size_t/ptrdiff_t are neither defined with base types
	// int nor long (because long is 32bit there!). In this case we need to
	// define some functions taking extra versions of integer arguments
	// in order to avoid risking ambiguous conversion errors. Fortunately,
    // at least 'int' and 'long' are always treated as different types (even
    // if they are identical)
#endif

#include "GL/glew.h"
// ^- only needed in a few places, but MUST be included before gl.h,
//    if gl.h is used anywhere (e.g., by QGLWidget...)

#include <cmath>
#include <stdint.h>
#include <cstddef>
#include <stdexcept>
#include <stdarg.h>

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

#include "vector_math.h"
#include "CxTypes.h"
#include "CxVec3.h"
#include "CxPodArray.h"

// extern bool
//    // global flag set when starting the program to tell it that we are analyzing files
//    // exported from the semi-empirical GFN-XTB program of Grimme. This requires some
//    // global adjustments (e.g., extracting the minimal basis from the main basis,
//    // which in this mode is near-minimal, rather than using MINAO; or for calculating
//    // charges taking into account the non-existing electrons in the core).
//    g_SemiEmpiricalFileMode;

using ct::TArray;
using ct::TVector3;
using std::size_t;
using std::ptrdiff_t;

typedef vmath::vec3<double> vec3d;
typedef vmath::vec3<float> vec3f;
typedef vmath::mat3<double> mat3d;
typedef vmath::mat3<float> mat3f;

// for homogeneous trafos.
typedef vmath::vec4<double> vec4d;
typedef vmath::vec4<float> vec4f;
typedef vmath::mat4<double> mat4d;
typedef vmath::mat4<float> mat4f;

// for vertex formats
typedef vmath::vec2<float> vec2f;
typedef vmath::vec3<unsigned int> vec3ui;

// hm... probably should get rid of either vector_math.h or CxVec3.h...
typedef ct::TVector3<double> FVec3d;
typedef ct::TVector3<float> FVec3f;
typedef ct::TVector3<unsigned int> FVec3ui;
using ct::sqr;

#include <QString>
#include <string>
std::string q2s(QString const &s);
QString s2q(std::string const &s);
// std::string RemovePath(std::string const &FileName);
// std::string RemoveExt(std::string const &FileName);
// std::string ReplaceExt(std::string const &FileName, std::string const &NewExt);
QString RemovePath(QString const &FileName);
QString RemoveExt(QString const &FileName);
QString ReplaceExt(QString const &FileName, QString const &NewExt);

class QMainWindow;
extern QMainWindow *g_pMainWindow; // for usage as dialog parent.

extern bool g_ShowVirtualOrbitals;
extern int g_nMaxOmpThreads;

enum FNotificationClass {
   NOTIFY_Information,
   NOTIFY_StartWork,
   NOTIFY_FinishWork,
   NOTIFY_Warning,
   NOTIFY_Error,
   NOTIFY_Unknown
};
// void IvNotify(FNotificationClass NotifyClass, std::string const &Text, std::string const &LongText = "");
void IvNotify(FNotificationClass NotifyClass, QString const &Text);
void IvWarn(std::string const &Text);
void IvWarn(QString const &Text);
void IvEmit(QString const &Text);

// an actual log class might be more useful...  could even inherit it from QStringList.
#define IV_NOTIFY(Class,Text) IvNotify((Class), QString(__FUNCTION__) + ":\t" + QString(Text))
#define IV_NOTIFY1(Class,Text,A1) IV_NOTIFY(Class, QString((Text)).arg((A1)))
#define IV_NOTIFY2(Class,Text,A1,A2) IV_NOTIFY(Class, QString((Text)).arg((A1)).arg((A2)))
#define IV_NOTIFY3(Class,Text,A1,A2,A3) IV_NOTIFY(Class, QString((Text)).arg((A1)).arg((A2)).arg((A3)))


namespace fmt {
   // these ones use boost formats. Maybe I should switch to cppformat instead, globally?
   // that still would require all the iostream stuff and lead to problems with utf8 vs char*,
   // but it would solve many other problems.
   //
   // note: for a start I could just use the trivial "replace {0} by arg0"-thing I put into
   //       script.format...

//    QString fmti(QString const &fmt, long i);
//    QString fmti(QString const &fmt, unsigned long i);
//    QString fmti(QString const &fmt, int i);
//    QString fmti(QString const &fmt, unsigned int i);
//    QString fmtf(QString const &fmt, float f);
//    QString fmtf(QString const &fmt, double f);
//    QString fmts(QString const &fmt, QString const &s);

   // these ones use explicit width/prec/base specifications (note: still
   // options missing... e.g., field alignment, leading zeros, 'always use sign',
   // etc..).
   QString fmti(int i, int width=0, int base=0);
   QString fmti(unsigned int i, int width=0, int base=0);
   QString fmti(long i, int width=0, int base=0);
   QString fmti(unsigned long i, int width=0, int base=0);
#ifdef NEED_EXTRA_SIZE_T_TYPES
   QString fmti(size_t i, int width=0, int base=0);
   QString fmti(ptrdiff_t i, int width=0, int base=0);
#endif
   QString fmtf(float f, int width=0, int prec=-1, char format='f');
   QString fmtf(double f, int width=0, int prec=-1, char format='f');
   QString fmts(QString const &s, int width=0);
   QString fmts(char const *p, int width=0);

   // format pointer
   QString fmtp(void *p);
   // format float as exponential (with explicit width/base)
   inline QString fmte(double f, int width=0, int prec=-1) { return fmtf(f,width,prec,'e'); }
   // format int as hexadecimal
   inline QString fmtx(long i, int width=0, int base=0) { return fmti(i,width,base); };
   inline QString fmtx(unsigned long i, int width=0, int base=0) { return fmti(i,width,base); };
}
using namespace fmt;


// dummy class which auto-converts stuff to QString via implicit constructors.
// It is type safe, but not very flexible (e.g., no field width, no precision, etc).
struct FEmitArg {
   FEmitArg(int n) : m(QString::number(n)) {}
   FEmitArg(uint n) : m(QString::number(n)) {}
//    FEmitArg(std::ptrdiff_t n) : m(QString::number(ulong(n))) {}
//    FEmitArg(std::size_t n) : m(QString::number(long(n))) {}
   FEmitArg(long n) : m(QString::number(n)) {}
   FEmitArg(unsigned long n) : m(QString::number(n)) {}
#ifdef NEED_EXTRA_SIZE_T_TYPES
   FEmitArg(size_t n) : m(QString::number(n)) {}
   FEmitArg(ptrdiff_t n) : m(QString::number(n)) {}
#endif
   FEmitArg(double n) : m(QString::number(n)) {}
   FEmitArg(float n) : m(QString::number(n)) {}
   FEmitArg(QString s) : m(s) {}
   FEmitArg(std::string s) : m(QString::fromStdString(s)) {}
   FEmitArg(char const *p) : m(QString(p)) {}
   QString
      m;
};
QString IvFmt(QString Text, FEmitArg const &a0);
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1);
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2);
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3);
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4);
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a6);
QString IvFmt(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a6, FEmitArg const &a7);
void IvEmit(QString Text, FEmitArg const &a0);
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1);
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2);
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3);
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4);
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a6);
void IvEmit(QString Text, FEmitArg const &a0, FEmitArg const &a1, FEmitArg const &a2, FEmitArg const &a3, FEmitArg const &a4, FEmitArg const &a6, FEmitArg const &a7);
// void IvEmit(QString const &Text, ...);
// void IvEmit(QString const &Text, va_list va);

QString LoadTextFromFile(QString FileName);


/// NOTE: This is c/p code from CxIterTools.h.
//
/// Defines a minimal iterator interface to allow use in range-based for loops
/// and generic container printing mechanism. This one is for objects which have
/// a working [] operator to access elements, but not necessarily much else.
/// Examples are raw-pointer based arrays (with non-unit stride, otherwise one
/// can just use the pointers directly).
/// Usage: Add:
///
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,const_iterator,const)
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,iterator,)
///
/// to class declaration to make both const and non-const iterators which wrap
/// calls to operator [].
#define IV_IMPLEMENT_ITERATOR_IndexInBracket(HostType,IteratorType,ConstOrNothing) \
   struct IteratorType { \
      HostType ConstOrNothing \
         *m_pArrayObj; \
      size_t \
         m_Index; \
      IteratorType &operator ++() { ++m_Index; return *this; } \
      bool operator != (IteratorType const& other) const { \
         return this->m_pArrayObj != other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      bool operator == (IteratorType const& other) const { \
         return this->m_pArrayObj == other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      auto operator *() -> decltype((*m_pArrayObj)[m_Index]) { return (*m_pArrayObj)[m_Index]; } \
   }; \
   IteratorType begin() ConstOrNothing { return IteratorType{this, 0}; } \
   IteratorType end() ConstOrNothing { return IteratorType{this, size()}; }

namespace iv {
   template<class... T> class MagicTypeRevealer;
}


#endif // CT_IV_H
