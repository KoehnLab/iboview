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

// https://doc.qt.io/archives/qt-4.8/ecmascript.html
// https://lost-contact.mit.edu/afs/cs.wisc.edu/i386_rh62/site/unsup/qt-4.8.0/i386_rhel5/doc/html/scripting.html#ecmascript-compatibility
#include <list>
#include "Iv.h"
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <sstream>

#include <QScriptEngine>
#include <QScriptValueIterator>
#include <QTextStream>
#include <QFile>

#include "IvScript.h"
#include "CxColor.h"
#include "IvDataOptions.h" // <-- we'll make c'tors from some of those.


// some other info, for ideas on processing vectors, lines, etc:
// https://stackoverflow.com/questions/2652390/qtscript-passing-an-array-of-objects-to-c


IApplication::~IApplication()
{}

IView3d::~IView3d()
{}

// void IView3d::set_option(QString const &OptionName, bool f) {
//    if (f)
//       set_option(OptionName, int(1));
//    else
//       set_option(OptionName, int(-1));
// }
// void IView3d::set_option(QString const &OptionName, int i) {
//    set_option(OptionName, double(i));
// }


// // will load an entire file into the stringstream str. Returns false if failed.
// bool LoadFileIntoMemory( std::string &sFileContent,
//         std::string const &FileName, unsigned *pFileLength )
// {
//     // read the entire file into an stringstream object.
//     std::ifstream
//         File( FileName.c_str() );
//     std::size_t
//         FileLength;
//     if ( false == File.good() )
//         return false;
//     File.seekg( 0, std::ios::end );
//     FileLength = File.tellg();
//     if ( 0 != pFileLength )
//         *pFileLength = FileLength;
//     sFileContent.clear();
//     sFileContent.resize(FileLength, 0);
//     File.seekg( 0, std::ios::beg );
//     File.read( &sFileContent[0], FileLength );
// //     pFileContent.resize(2 + FileLength);
// //     memset( &pFileContent[0], 0, 2 + FileLength );
// //     File.seekg( 0, std::ios::beg );
// //     File.read( &pFileContent[0], FileLength );
//     return true;
// };

uint32_t hsv_uint(float h, float s, float v) {
   return 0xff000000 | (uint32_t)ct::Hsv(h,s,v).uint32();
//    return 0xff000000 + (uint32_t)ct::FColor(h,s,v);
//    return 0xfffffffe;
//    return (uint32_t)ct::Hsv(h,s,v);
}

uint32_t ColFromAlpha(float a);

uint32_t hsva_uint(float h, float s, float v, float a) {
   return ColFromAlpha(a) | (uint32_t)ct::Hsv(h,s,v).uint32();
//    return 0xff000000 + (uint32_t)ct::FColor(h,s,v);
//    return 0xfffffffe;
//    return (uint32_t)ct::Hsv(h,s,v);
}

// static std::string s_format_i(std::string const &Format, int val) {
//    std::stringstream str;
//    str << boost::format(Format) % val;
//    return str.str();
// }
//
// static std::string s_format_f(std::string const &Format, double val) {
//    std::stringstream str;
//    str << boost::format(Format) % val;
//    return str.str();
// }
//
// static std::string s_format_s(std::string const &Format, std::string const &val) {
//    std::stringstream str;
//    str << boost::format(Format) % val;
//    return str.str();
// }

static QString s_format_i(QString const &Format, int val) {
   QString r;
   r.sprintf(Format.toLatin1().constData(), val);
   return r;
}

static QString s_format_f(QString const &Format, double val) {
   QString r;
   r.sprintf(Format.toLatin1().constData(), val);
   return r;
}

static QString s_format_s(QString const &Format, QString const &val) {
   QString r;
   r.sprintf(Format.toLatin1().constData(), val.toUtf8().constData());
   return r;
}


QString FormatArray(QScriptValue const &ScriptList, bool AddQuotesToStrs, QString Left, QString Separator, QString Right)
{
//    assert(ScriptList.isArray());
   assert(ScriptList.property("length").isValid());
   QString
      sout;
   QTextStream
      str(&sout);
   str << Left;
   int Length = ScriptList.property("length").toInteger();
   for (int i = 0; i < Length; ++ i) {
      if (i != 0) {
         str << Separator;
      }
      QScriptValue
         v = ScriptList.property(i);
      if (v.isString() && AddQuotesToStrs) {
         str << "\"" << v.toString() << "\"";
      } else if (v.isArray()) {
         str << FormatArray(v, AddQuotesToStrs, Left, Separator, Right);
      } else {
         str << v.toString();
      }
   }
   str << Right;
   str.flush();
//    IvNotify(NOTIFY_Error, sout);
   return sout;
}


static QScriptValue s_Invoke_format(QScriptContext *context, QScriptEngine * /*engine*/)
{
   int
      iArg0,
      nArgs = context->argumentCount();
   QString
      s;
//    if (context->thisObject().isUndefined()) {
   if (!context->thisObject().isString()) {
      // called as stand-alone function. Take first argument as base string.
      // (the 'thisObject' is still defined: it is returned as "[object global]".)
      s = qscriptvalue_cast<QString>(context->argument(0));
      // IvEmit("!! s_Invoke_format: Taking free-standing version. s = '%1'", s);
      iArg0 = 1;
      nArgs -= 1;
   } else {
      // take "this" as base string and rest as arguments.
      s = qscriptvalue_cast<QString>(context->thisObject());
      // IvEmit("!! s_Invoke_format: Taking object version. s = '%1'", s);
      iArg0 = 0;
   }
//    for (int i = iArg0; i < nArgs; ++ i)
//       s = s.arg(qscriptvalue_cast<QString>(context->argument(i)));
#if 0
   // can only do basic replacements.
   for (int i = 0; i < nArgs; ++ i) {
      s.replace(QString("{%1}").arg(i), qscriptvalue_cast<QString>(context->argument(iArg0+i)));
   }
#else
   fmt::MemoryWriter w;
   uint64_t
      Types = 0,
      Shift = 0;
   std::vector<fmt::internal::Value>
      Args;
#define ADD_ARG(v) {\
      fmt::internal::Value value = fmt::internal::MakeValue<char>(v); \
      Args.push_back(value); \
      Types |= fmt::internal::make_type(v) << Shift; \
      Shift += 4; \
   }
   std::list<std::string>
      StrBuf;
   if (nArgs > fmt::ArgList::MAX_ARGS)
      IvNotify(NOTIFY_Warning, "Script invoked format() function with more than supported number of arguments. Some conversions will not be done. Sorry :(");
   for (int i = 0; i < nArgs; ++ i) {
      QScriptValue arg = context->argument(iArg0+i);
      if (arg.isNumber()) {
         qsreal q = arg.toNumber();
         // ECMA script does not distinguish between floats and integers...
         if (q == long(q)) {
            ADD_ARG(long(q));
         } else {
            ADD_ARG(double(q));
         }
      } else if (arg.isBool()) {
         ADD_ARG(bool(arg.toBool()));
//       } else if (arg.isArray()) {
      } else if (arg.isString()) {
         StrBuf.push_back(arg.toString().toStdString());
         ADD_ARG(StrBuf.back());
      } else if (arg.property("length").isValid()) {
         StrBuf.push_back(FormatArray(arg,true,"[",", ","]").toStdString());
         ADD_ARG(StrBuf.back());
      } else {
         StrBuf.push_back(arg.toString().toStdString());
         ADD_ARG(StrBuf.back());
      }
   }
#undef ADD_ARG
   fmt::ArgList
      args(Types, &Args[0]);
   w.write(s.toStdString(), args);
   s = QString::fromStdString(w.str());


//    fmt::BasicFormatter<char> const &cstr = w.Format(s.toStdString());
//    fmt::BasicFormatter<char> &str = const_cast<fmt::BasicFormatter<char>&>(cstr);
//
//    for (int i = 0; i < nArgs; ++ i) {
//       QScriptValue arg = context->argument(iArg0+i);
//       if (arg.isNumber()) {
//          qsreal q = arg.toNumber();
//          // ECMA script does not distinguish between floats and integers...
//          if (q == long(q))
//             str << long(q);
//          else
//             str << q;
//       } else {
//          str << arg.toString().toStdString();
//       }
//       s = QString::fromStdString(w.str());
//    }
#endif
   return QScriptValue(s);
}

static QScriptValue s_Invoke_join(QScriptContext *context, QScriptEngine * /*engine*/)
{
   int
      iArg0,
      nArgs = context->argumentCount();
   QString
      s;
//    if (context->thisObject().isUndefined()) {
   if (!context->thisObject().isString()) {
      // called as stand-alone function. Take first argument as base string.
      // (the 'thisObject' is still defined: it is returned as "[object global]".)
      s = qscriptvalue_cast<QString>(context->argument(0));
      iArg0 = 1;
      nArgs -= 1;
   } else {
      // take "this" as base string and rest as arguments.
      s = qscriptvalue_cast<QString>(context->thisObject());
      iArg0 = 0;
   }
   QScriptValue arg = context->argument(iArg0);
   if (nArgs != 1 || !arg.isValid() || !arg.property("length").isValid()) {
      IvNotify(NOTIFY_Error, "s_Invoke_join called with unsupported arguments");
   }
   return QScriptValue(FormatArray(arg, false, "", s, ""));
}


#define FARG(FType,iArg) qscriptvalue_cast<FType>(context->argument((iArg)))

//     QScriptValue callee = context->callee();
//      if (context->argumentCount() == 1) // writing?
//          callee.setProperty("value", context->argument(0));
//      return callee.property("value");
//    return QScriptValue(hsva_uint(context->argument(0).toUInt32()));
static QScriptValue s_Invoke_hsva(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(hsva_uint(FARG(float, 0), FARG(float, 1), FARG(float,2), FARG(float,3)));
}

// static QScriptValue s_Invoke_fmti(QScriptContext *context, QScriptEngine * /*engine*/) {
//    return QScriptValue(QString::fromStdString(s_format_i(FARG(QString,0).toStdString(), FARG(int,1))));
// }
// static QScriptValue s_Invoke_fmtf(QScriptContext *context, QScriptEngine * /*engine*/) {
//    return QScriptValue(QString::fromStdString(s_format_f(FARG(QString,0).toStdString(), FARG(double,1))));
// }
//
// static QScriptValue s_Invoke_fmts(QScriptContext *context, QScriptEngine * /*engine*/) {
//    return QScriptValue(QString::fromStdString(s_format_s(FARG(QString,0).toStdString(), FARG(QString,1).toStdString())));
// }
//
// ^- note: we can probably do this much better now with QtScript...

static QScriptValue s_Invoke_fmti(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(s_format_i(FARG(QString,0), FARG(int,1)));
}

static QScriptValue s_Invoke_fmtf(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(s_format_f(FARG(QString,0), FARG(double,1)));
}

static QScriptValue s_Invoke_fmts(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(s_format_s(FARG(QString,0), FARG(QString,1)));
}

static QScriptValue s_Invoke_irgb(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(ct::irgb(FARG(quint32,0)));
}

static QScriptValue s_Invoke_replace_ext(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(ReplaceExt(FARG(QString,0), FARG(QString,1)));
}

static QScriptValue s_Invoke_remove_path(QScriptContext *context, QScriptEngine * /*engine*/) {
   return QScriptValue(RemovePath(FARG(QString,0)));
}



QString ArgsToString(QScriptContext *context) {
   QString out;
   QTextStream str(&out);
//    str << "(";
   for (int i = 0; i < int(context->argumentCount()); ++ i) {
      if (i != 0) str << ", ";
      str << context->argument(i).toString();
   }
//    str << ")";
   return out;
}

FVec3d FVec3d_FromScriptValue(QScriptValue array, QScriptContext *context) {
   typedef FVec3d
      FVecN;
   FVecN v;
   int length = array.property("length").toInteger();
   if (length != (int)v.size()) {
      context->throwError(QScriptContext::RangeError, "FVec3d.ctor(array): expected array of length three.");
      return FVec3d(0,0,0);
   }
   for (int i = 0; i < length; i++)
      v[i] = array.property(i).toNumber();
   return v;
}

QScriptValue FVec3d_ctor(QScriptContext *context, QScriptEngine *engine)
{
   typedef FVec3d
      FVecN;
   FVecN v;
   if (0) {
      IvEmit("| #invk: FVec3d_ctor(%1)", ArgsToString(context));
   }
   if (size_t(context->argumentCount()) == v.size()) {
      for (int i = 0; i < (int)v.size(); ++ i)
         v[i] = context->argument(i).toNumber();
      return engine->toScriptValue(v);
   } else if (context->argumentCount() == 1) {
// //       QVariantList vl = qscriptvalue_cast<QVariantList>(context->argument(0));
// //       if (vl.size() != (int)v.size())
// //          return context->throwError(QScriptContext::RangeError, "FVec3d.ctor(array): expected array of length three.");
// //       for (int i = 0; i < (int)v.size(); ++ i)
// //          v[i] = vl[i].toNumber();
// //       return engine->toScriptValue(v);
//       QScriptValue array = context->argument(0);
//       int length = array.property("length").toInteger();
//       if (length != (int)v.size())
//          return context->throwError(QScriptContext::RangeError, "FVec3d.ctor(array): expected array of length three.");
//       for (int i = 0; i < length; i++)
//          v[i] = array.property(i).toNumber();
      v = FVec3d_FromScriptValue(context->argument(0), context);
//          qDebug() << array.property(i).toString();
      return engine->toScriptValue(v);
   } else {
      return context->throwError(QScriptContext::TypeError, "FVec3d.ctor: expected components or list of components as arguments.");
   }
}

// function to get/set a vector component. Index of component is stored as data object on the prototype function.
QScriptValue FVec3d_prototype_gsetX(QScriptContext *context, QScriptEngine *engine)
{
   size_t iEntry = size_t(qscriptvalue_cast<int>(context->callee().data()));

   // Cast to a pointer to be able to modify the underlying C++ value
   FVec3d *vptr = qscriptvalue_cast<FVec3d*>(context->thisObject());
   if (!vptr)
      return context->throwError(QScriptContext::TypeError, "FVec3d.prototype.gsetX: this object is not a FVec3d");
   if (iEntry >= vptr->size())
      return context->throwError(QScriptContext::RangeError, QString("FVec3d.prototype.getSetX: attempted access to non-existent vector component %1").arg(iEntry));
   if (context->argumentCount() == 1) {
      // called as "setter" (1 argument)
      (*vptr)[iEntry] = context->argument(0).toNumber();
      return engine->undefinedValue();
   } else {
      return (*vptr)[iEntry];
   }
}

QScriptValue FVec3d_prototype_getLength(QScriptContext *context, QScriptEngine *engine)
{
   typedef FVec3d FVecN;
   FVecN v;
   return engine->toScriptValue(v.size());
   (void)context; // suppress unused warning
}



QScriptValue FVec3d_prototype_toString(QScriptContext *context, QScriptEngine *engine) {
   FVec3d *vptr = qscriptvalue_cast<FVec3d*>(context->thisObject());
   if (!vptr)
      return context->throwError(QScriptContext::TypeError, "FVec3d.prototype.toString: this object is not a FVec3d");
   QString
      res;
   QTextStream
      str(&res);
//    bool AsArray = true;
   bool AsArray = false;
   str << (AsArray ? "[" : "Vec3(");
   for (size_t i = 0; i < vptr->size(); ++ i) {
      if (i != 0) str << ",";
      str << (*vptr)[i];
   }
   str << (AsArray ? "]" : ")");
   return engine->toScriptValue(res);
}


// """A single Qt Script function can act as both getter and setter for a property. When it is called as a getter, the argument count is 0. When it is called as a setter, the argument count is 1; the argument is the new value of the property. In the following example, we define a native combined getter/setter that transforms the value slightly:"""
// ^-- @https://doc.qt.io/qt-5/qtscript-index.html
//
// QScriptValue FVec3d_prototype_setX(QScriptContext *context, QScriptEngine *engine)
// {
//    // Cast to a pointer to be able to modify the underlying C++ value
//    FVec3d *vptr = qscriptvalue_cast<FVec3d*>(context->thisObject());
//    if (!vptr)
//       return context->throwError(QScriptContext::TypeError, "FVec3d.prototype.setX: this object is not a FVec3d");
//    (*vptr)[0] = context->argument(0).toNumber();
//    return engine->undefinedValue();
// }

// https://stackoverflow.com/questions/9265288/casting-a-list-as-qvariant-or-qvariant-list
//
// // Note QVariant itself is already able to convert from QList<T> to QVariantList. This means you could also use QVariant to perform the conversion.
// https://www.qtcentre.org/threads/21161-QtScript-and-arrays



// Q_DECLARE_METATYPE(QSharedPointer<FFreeLine>)
//
// QScriptValue FFreeLine_ctor(QScriptContext *context, QScriptEngine *engine)
// {
//    int nArgs = context->argumentCount();
//    if (nArgs < 2)
//       return context->throwError(QScriptContext::RangeError, "FreeLine.ctor(): expected at least two arguments (from,to)");
// //    QSharedPointer<FFreeLine> pobj(new FFreeLine(dummy.GetFrom(), dummy.GetTo(), dummy.GetColor(), dummy.GetWidth()));
//    QSharedPointer<FFreeLine> pobj(new FFreeLine());
//    int iArg = 0;
//    if (nArgs > iArg) { pobj->SetFrom(qscriptvalue_cast<FVec3d>(context->argument(iArg))); iArg += 1; }
//    if (nArgs > iArg) { pobj->SetTo(qscriptvalue_cast<FVec3d>(context->argument(iArg))); iArg += 1; }
//    if (nArgs > iArg) { pobj->SetColor(qscriptvalue_cast<uint32_t>(context->argument(iArg))); iArg += 1; }
//    if (nArgs > iArg) { pobj->SetWidth(qscriptvalue_cast<double>(context->argument(iArg))); iArg += 1; }
//
//    // store the shared pointer in the script object that we are constructing
//    return engine->newVariant(context->thisObject(), QVariant::fromValue(pobj));
//    // ^-- note: didn't think about this very deeply. Just c/p'd from the QSharedPointer example
//    // in https://doc.qt.io/qt-5/qtscript-index.html
// }

Q_DECLARE_METATYPE(QSharedPointer<FFreeLine>)

// https://dreamswork.github.io/qt4/classQScriptEngine.html#aca0d5d8173b864830147da04467ebab9
// template<typename T >
// T qscriptvalue_cast	(	const QScriptValue & 	value	)
// {
//      T t;
//      const int id = qMetaTypeId<T>();
//
//      if (qscriptvalue_cast_helper(value, id, &t))
//          return t;
//      else if (value.isVariant())
//          return qvariant_cast<T>(value.toVariant());
//
//      return T();
//  }


// void _DebugPrintScriptValueProperties(QString desc, QScriptValue obj, bool DataMembersOnly=true) {
// //    QScriptValueIterator it(object);
// //    while (it.hasNext()) {
// //       it.next();
// //       qDebug() << it.name() << ": " << it.value().toString();
// //    }
//    QScriptValue proto = obj;
//    while (proto.isObject()) {
//       QScriptValueIterator it(proto);
//       while (it.hasNext()) {
//          it.next();
//          QScriptValue val = obj.property(it.name());
//          if (DataMembersOnly && (!val.isFunction() || bool(it.flags() & QScriptValue::PropertyGetter)))
//             IvEmit("  [%1] %2::%3 (fn? {})= %4", desc, it.scriptName().toString(), it.name(), val.toString());
// //             IvEmit("  [%1] %2::%3 (fn? {})= %4", desc, it.scriptName().toString(), it.name(), it.value().toString());
//       }
//       proto = 0;
//       // ^-- do not iterate over properties from prototypes.
// //       proto = proto.prototype();
//    }
//
// }
void _DebugPrintScriptValueProperties(QString desc, QScriptValue obj, bool DataMembersOnly=true) {
   QScriptValueIterator it(obj);
   while (it.hasNext()) {
      it.next();
//       QScriptValue val = obj.property(it.name());
      QScriptValue val = it.value();
      if (DataMembersOnly && (!val.isFunction() || bool(it.flags() & QScriptValue::PropertyGetter)))
//          IvEmit("  [%1] %2::%3 (fn? %4) = %5", desc, it.scriptName().toString(), it.name(), val.isFunction(), val.toString());
         IvEmit("  [%1] %2 (fn? %4) = %5", desc, it.name(), val.isFunction(), val.toString());
   }
}


struct FVariantFmt {
//    enum {
//       delta_ind = 2
//    };
   void _DebugFmt(QTextStream &ss, QString desc, QVariantMap vm, int level) {
      ss << "{";
      for (QVariantMap::const_iterator it = vm.begin(); it != vm.end(); ++ it) {
         if (it != vm.begin())
            ss << ", ";
         ss << it.key() << ": ";
         _DebugFmt(ss, desc, it.value(), level + 1);
      }
      ss << "}";
   }
   void _DebugFmt(QTextStream &ss, QString desc, QVariantList vm, int level) {
      ss << "[";
      for (QVariantList::const_iterator it = vm.begin(); it != vm.end(); ++ it) {
         if (it != vm.begin())
            ss << ", ";
         _DebugFmt(ss, desc, *it, level + 1);
      }
      ss << "]";
   }

   void _DebugFmt(QTextStream &ss, QString desc, QVariant variant, int level = 0) {
   // see: https://doc.qt.io/qt-5/qvariant.html
   // hmpf... QAssociativeIterable is QT 5.2+ only.
//    if (variant.canConvert<QVariantMap>()) {
//       QAssociativeIterable iterable = variant.value<QAssociativeIterable>();
//       // Can use foreach over the values:
//       foreach (const QVariant &v, iterable) {
//          qDebug() << v;
//       }
//       // Can use C++11 range-for over the values:
//       for (const QVariant &v : iterable) {
//          qDebug() << v;
//       }
//       // Can use iterators:
//       QAssociativeIterable::const_iterator it = iterable.begin();
//       const QAssociativeIterable::const_iterator end = iterable.end();
//       for ( ; it != end; ++it) {
//          qDebug() << *it; // The current value
//          qDebug() << it.key();
//          qDebug() << it.value();
//       }
//    }
      if (variant.canConvert<QVariantMap>()) {
         _DebugFmt(ss, desc, variant.value<QVariantMap>(), level+1);
      } else if (variant.canConvert<QVariantList>()) {
         _DebugFmt(ss, desc, variant.value<QVariantList>(), level+1);
      } else if (variant.canConvert<FVec3d>()) {
         FVec3d const &v = variant.value<FVec3d>();
         ss << "Vec3(" << v[0] << "," << v[1] << "," << v[2] << ")";
//       } else if (variant.canConvert<uint32_t>()) {
//          uint32_t const &v = variant.value<uint32_t>();
//          if (v > 20) {
//             ss << s2q(fmt::format("0x{:08x}", v));
//          }
//          ss << s2q(fmt::format("0x{:08x}", v));
//       } else if (variant.canConvert<QString>()) {
//          ss << "\"" << variant.toString() << "\"";

      } else if (variant.type() == QVariant::Type(QMetaType::QString)) {
         // stored as actual string. Print it with quotes (should probably add escapes...)
//          ss << "\"" << variant.toString().replace("\\\"","") << "\"";
         ss << "\"" << variant.toString().replace("\"","\\\"") << "\"";
      } else if (variant.canConvert<double>()) {
         // ^-- interestingly, variant.canConvert<double>() returns 'True' for a
         // variant containing a QString with value "whooo?". Not exactly what I
         // would have expected...
         double const &v = variant.value<double>();
         if (v >= 0 && double(uint32_t(v)) == v) {
            uint32_t const &v_ui = uint32_t(v);
            if (v_ui >= 0xff) {
               // large unsigned integer? maybe it's a color?
               if (v_ui > 0xffffff)
                  ss << s2q(fmt::format("0x{:08x}", v_ui));
               else
                  ss << s2q(fmt::format("0x{:06x}", v_ui));
            } else {
               ss << v_ui;
            }
         } else {
            ss << v;
//             ss << variant.toString();
         }
      } else if (variant.canConvert<QString>()) {
         // not a string, but something which can be converted to string. stringify it, but do not add quotes.
         ss << variant.toString();
      } else {
         ss << "unrecognized_variant_type";
      }
   }

   void _DebugPrint(QString desc, QVariant variant) {
      QString out;
      QTextStream ss(&out);
      _DebugFmt(ss, desc, variant);
      IvEmit("_DebugPrint(v) = '%1'", out);
   }
};


QScriptValue FFreeLine_ctor(QScriptContext *context, QScriptEngine *engine)
{
   int nArgs = context->argumentCount();
   if (nArgs < 2)
      return context->throwError(QScriptContext::RangeError, "FreeLine.ctor(): expected at least two arguments (from,to)");
   FFreeLine *pobj(new FFreeLine(engine));
   int iArg = 0;
//    FVec3d dummy_from = qscriptvalue_cast<FVec3d>(context->argument(0));
//    IvEmit("dummy_from = %1", dummy_from[0]);
//    if (nArgs > iArg) { pobj->SetFrom(qscriptvalue_cast<FVec3d>(context->argument(iArg))); iArg += 1; }
//    if (nArgs > iArg) { pobj->SetTo(qscriptvalue_cast<FVec3d>(context->argument(iArg))); iArg += 1; }
   if (nArgs == 1) {
      pobj->setp(qscriptvalue_cast<QVariantMap>(context->argument(iArg)));
   } else if (nArgs >= 3) {
      if (nArgs > iArg) { pobj->SetFrom(FVec3d_FromScriptValue(context->argument(iArg), context)); iArg += 1; }
      if (nArgs > iArg) { pobj->SetTo(FVec3d_FromScriptValue(context->argument(iArg), context)); iArg += 1; }
   //    if (nArgs > iArg) { pobj->SetColor(qscriptvalue_cast<uint32_t>(context->argument(iArg))); iArg += 1; }
   //    if (nArgs > iArg) { pobj->SetWidth(qscriptvalue_cast<double>(context->argument(iArg))); iArg += 1; }
      if (nArgs > iArg) { pobj->setp(qscriptvalue_cast<QVariantMap>(context->argument(iArg))); iArg += 1; }
   }

//    _DebugPrintScriptValueProperties("prop/last-arg", context->argument(nArgs-1));
//    IvEmit("context->argument(0).toString() = '%1'", context->argument(0).toString());
//    _DebugPrintScriptValueProperties("prop/arg0", context->argument(0));
//    IvEmit("context->argument(0).toString() = '%1'", context->argument(0).toString());
   if (0) {
      FVariantFmt fmt;
      fmt._DebugPrint("prop/last-arg", context->argument(nArgs-1).toVariant());
   }


   // let the engine manage the new object's lifetime.
   return engine->newQObject(pobj, QScriptEngine::ScriptOwnership);
}



// QScriptValue Person_ctor(QScriptContext *ctx, QScriptEngine *eng)
// {
//     QScriptValue object;
//     if (ctx->isCalledAsConstructor()) {
//         object = ctx->thisObject();
//     } else {
//         object = eng->newObject();
//         object.setPrototype(ctx->callee().property("prototype"));
//     }
//     object.setProperty("name", ctx->argument(0));
//     return object;
// }


#undef FARG


static void AddScriptObject(QScriptEngine &ScriptEngine, QObject *pObject, QString const &Name)
{
   QScriptValue
//       ObjectValue = ScriptEngine.newQObject(pObject, QScriptEngine::QtOwnership, QScriptEngine::ExcludeChildObjects | QScriptEngine::ExcludeSuperClassContents);
//       ObjectValue = ScriptEngine.newQObject(pObject, QScriptEngine::QtOwnership, QScriptEngine::ExcludeSuperClassContents);
      ObjectValue = ScriptEngine.newQObject(pObject, QScriptEngine::QtOwnership, QScriptEngine::ExcludeSuperClassMethods);
   ScriptEngine.globalObject().setProperty(Name, ObjectValue);
}

static void AddScriptFunction(QScriptEngine &ScriptEngine, QScriptEngine::FunctionSignature pFn, QString const &Name)
{
   QScriptValue
      FunctionValue = ScriptEngine.newFunction(pFn);
   ScriptEngine.globalObject().setProperty(Name, FunctionValue);
}

void ExecScript(IApplication *app, IView3d *view, QString const &ScriptText, QString const &ScriptName) {
//    std::cout << QString("--- exec script '%1'").arg(ScriptName).toStdString() << std::endl;
   IvEmit("--- exec script '%1'", ScriptName);
   QScriptEngine
      ScriptEngine;

   AddScriptObject(ScriptEngine, app, "app");
   // FIXME: added IApplication as both "app" and "doc" until things get sorted out.
   AddScriptObject(ScriptEngine, app, "doc");
   AddScriptObject(ScriptEngine, view, "view");
   AddScriptFunction(ScriptEngine, s_Invoke_hsva, "hsva");
   AddScriptFunction(ScriptEngine, s_Invoke_irgb, "irgb");
   AddScriptFunction(ScriptEngine, s_Invoke_fmti, "fmti");
   AddScriptFunction(ScriptEngine, s_Invoke_fmtf, "fmtf");
   AddScriptFunction(ScriptEngine, s_Invoke_fmts, "fmts");
   AddScriptFunction(ScriptEngine, s_Invoke_format, "format");
   AddScriptFunction(ScriptEngine, s_Invoke_join, "join");
   AddScriptFunction(ScriptEngine, s_Invoke_replace_ext, "replace_ext"); // replace a file extension by another file extension
   AddScriptFunction(ScriptEngine, s_Invoke_remove_path, "remove_path"); // strip off path from a file name

   // add 'format' also to prototype of standard class 'String'.
   // This should then allow for "blabla {0} blabla {1}".format("wheee","whooo")-type calls.
   QScriptValue
      FormatFunction = ScriptEngine.newFunction(s_Invoke_format),
      JoinFunction = ScriptEngine.newFunction(s_Invoke_join);
   QScriptValue
      StringProto = ScriptEngine.globalObject().property("String").property("prototype");
   StringProto.setProperty("format", FormatFunction);
   StringProto.setProperty("join", JoinFunction);


   ScriptEngine.globalObject().setProperty("Vec3", ScriptEngine.newFunction(FVec3d_ctor));
   ScriptEngine.globalObject().setProperty("FreeLine", ScriptEngine.newFunction(FFreeLine_ctor));

   {
      QScriptValue Vec3Proto = ScriptEngine.newObject();

      char const *CompNames[] = {"x","y","z","w"};
      size_t nCompNames = sizeof(CompNames)/sizeof(CompNames[0]);

      // make & register functions to get/set the vector components in array syntax.
      for (size_t i = 0; i < FVec3d::static_size; ++ i) {
         QScriptValue gsetComp_i = ScriptEngine.newFunction(FVec3d_prototype_gsetX);
         gsetComp_i.setData(ScriptEngine.toScriptValue(int(i)));
//          ScriptEngine.globalObject().setProperty(QString("gsetVec_i%1").arg(i), gsetComp_i); // <-- not really needed, is it?

         // also assign a array-index based property (v[0], v[1], ...).
         Vec3Proto.setProperty(quint32(i), gsetComp_i, QScriptValue::PropertyGetter | QScriptValue::PropertySetter);

         // also assign a name based propert (.x, .y, .z, .w... maybe add .rgba, too?)
         if (i < nCompNames)
            Vec3Proto.setProperty(CompNames[i], gsetComp_i, QScriptValue::PropertyGetter | QScriptValue::PropertySetter);
      }

      Vec3Proto.setProperty("length", ScriptEngine.newFunction(FVec3d_prototype_getLength), QScriptValue::PropertyGetter);
      Vec3Proto.setProperty("toString", ScriptEngine.newFunction(FVec3d_prototype_toString));

      ScriptEngine.setDefaultPrototype(qMetaTypeId<FVec3d*>(), Vec3Proto);
      ScriptEngine.setDefaultPrototype(qMetaTypeId<FVec3d>(), Vec3Proto);
//    https://doc.qt.io/qt-5/qtscript-script-defaultprototypes-example.html

      // hmpf. Do I need this? -> https://doc.qt.io/qt-5/qscriptengine.html#qScriptRegisterSequenceMetaType
      // won't quite work. Needs a push_back... not for fixed-size arrays :(
   }


//    std::string sInputScript;
//    LoadFileIntoMemory(sInputScript, FileName.toStdString(), 0);

   QScriptValue
//       res = ScriptEngine.evaluate(QString::fromStdString(sInputScript));
      res = ScriptEngine.evaluate(ScriptText);
   if (ScriptEngine.hasUncaughtException()) {
      int iErrorLine = ScriptEngine.uncaughtExceptionLineNumber();
//       std::cerr << QString("\n!ERROR during execution of script:\n%1:%2: %3\n").arg(FileName, QString::number(iErrorLine), res.toString()).toStdString() << std::endl;
      IvNotify(NOTIFY_Error, "Error during script execution:\n" + QString("%1:%2: %3\n").arg(ScriptName, QString::number(iErrorLine), res.toString()));
   }
   IvEmit("--- script terminated.");
}

QString LoadTextFileViaQt(QString const &FileName)
{
   QFile
      File(FileName);
   if (!File.open(QIODevice::ReadOnly | QIODevice::Text)) {
//       std::cerr << "!ERROR: Failed to open script file '%s'" % FileName << std::endl;
      IvNotify(NOTIFY_Error, QString("Failed to open script file '%1'").arg(FileName));
      return "";
   }
   QTextStream
      Stream(&File);
   Stream.setAutoDetectUnicode(true);
   // ^- deals with UTF-16 and UTF-32. I guess it will default to UTF8?
   //    I didn't quite get this in the docs.
   QString
      Result = Stream.readAll();
   File.close();
   return Result;
}


// void ExecScript(FMainWindow *ui, std::string FileName) {
void ExecScript(IApplication *app, IView3d *view, QString const &FileName)
{
//    IvEmit("--- exec script '%1'", FileName);
   QString
      ScriptText = LoadTextFileViaQt(FileName);
   ExecScript(app, view, ScriptText, FileName);
}
