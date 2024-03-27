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

/// @file CtIvStubs.h
///
/// Code for replacing some IboView functions by stubs/dummies/replacements
/// in order to simplify moving IboView code to MicroScf and/or using
/// MicroScf code in IboView. And yes, this code will not win any prices.

#ifdef PROG_IBOVIEW
   #include "Iv.h"
   // make sure the namespace is there (actual IboView currently does not put
   // anything into a namespace), so that we can get the Iv* stuff with a
   // function-local `using namespace iv;`, regardless of whether or not the
   // functions are actually inside it.
   namespace iv { }
#else

#ifndef CT_IBOVIEW_STUBS
#define CT_IBOVIEW_STUBS

#include <sstream>

namespace iv {
typedef std::string QString; // <- super hack!

static std::string const &q2s(QString const &s) { return s; }
static QString const &s2q(std::string const &s) { return s; }

enum FNotificationClass {
   NOTIFY_Information,
   NOTIFY_StartWork,
   NOTIFY_FinishWork,
   NOTIFY_Warning,
   NOTIFY_Error,
   NOTIFY_Unknown
};
void IvNotify(FNotificationClass NotifyClass, QString const &Text);
void IvWarn(std::string const &Text);

#define IV_NOTIFY(Class,Text) IvNotify((Class), QString(__FUNCTION__) + ":\t" + QString(Text))
#define IV_NOTIFY1(Class,Text,A1) IV_NOTIFY(Class, IvFmt((Text), (A1)))
#define IV_NOTIFY2(Class,Text,A1,A2) IV_NOTIFY(Class, IvFmt((Text), (A1), (A2)))
#define IV_NOTIFY3(Class,Text,A1,A2,A3) IV_NOTIFY(Class, IvFmt((Text), (A1), (A2), (A3)))


namespace iv_fmt {
   template<class FIntType>
   QString fmti(FIntType const &i, int width=0, int base=0) {
      std::stringstream str;
      if (base == 16)
         str << std::hex;
         // ^- we just ignore all others...
      if (width != 0)
         str.width(width);
      str << i;
      return str.str();
   }

   template<class FFloatType>
   QString fmtf(FFloatType const &f, int width=0, int prec=-1, char format='f') {
      std::stringstream str;
      if (format == 'f')
         str.setf(std::ios_base::fixed, std::ios_base::floatfield);
      else if (format == 'e')
         str.setf(std::ios_base::scientific, std::ios_base::floatfield);
      if (width != 0)
         str.width(width);
      if (prec != -1)
         str.precision(prec);
      str << f;
      return str.str();
   }
   template<class FStrType>
   QString fmts(FStrType const &s, int width=0) {
      if (width == 0)
         return s;
      else {
         std::stringstream str;
         str.width(width);
         str << s;
         return str.str();
      }
   }

   // format pointer
   QString fmtp(void *p);
   // format float as exponential (with explicit width/base)
   inline QString fmte(double f, int width=0, int prec=-1) { return fmtf(f,width,prec,'e'); }
   // format int as hexadecimal
   inline QString fmtx(long i, int width=0) { return fmti(i,width,16); };
   inline QString fmtx(unsigned long i, int width=0) { return fmti(i,width,16); };
}
using namespace iv_fmt;


struct FEmitArg {
   FEmitArg(int val) { ConvertToString(val); }
   FEmitArg(unsigned val) { ConvertToString(val); }
   FEmitArg(long val) { ConvertToString(val); }
   FEmitArg(unsigned long val) { ConvertToString(val); }
#ifdef NEED_EXTRA_SIZE_T_TYPES
   FEmitArg(size_t val) { ConvertToString(val); }
   FEmitArg(ptrdiff_t val) { ConvertToString(val); }
#endif
   FEmitArg(double val) { ConvertToString(val); }
   FEmitArg(float val) { ConvertToString(val); }
   FEmitArg(QString const &val) { ConvertToString(q2s(val)); }
   FEmitArg(char const *val) { ConvertToString(val); }
   QString
      m;
protected:
   template<class T>
   void ConvertToString(T const &value) {
      std::stringstream str;
      str << value;
      m = str.str();
   }
};

void ReplaceArgInPlace(QString &Text, FEmitArg const &arg);

template<typename... Ts> QString IvFmt(QString const &Text, Ts... Args) {
   QString TextWithArgs(Text);
   int dummy[sizeof...(Ts)] = { (ReplaceArgInPlace(TextWithArgs, Args), 0)... };
   // ^-- construction used to ensure ReplaceArgInPlace is called in sequence for the parameter pack
   //     expansion.
   (void)dummy; // <- to suppress "unused variable" warning
   return TextWithArgs;
}

void _IvEmit(QString const &Text);

template<typename... Ts> void IvEmit(QString const &Text, Ts... Args) {
   _IvEmit(IvFmt(Text, Args...));
}



} // namespace iv

#endif // CT_IBOVIEW_STUBS

#endif // PROG_IBOVIEW
