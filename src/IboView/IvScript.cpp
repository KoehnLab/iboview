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

// https://doc.qt.io/qt-6/qtqml-javascript-qtjavascript.html
// https://doc.qt.io/qt-6/qjsengine.html
#include <list>
#include "Iv.h"
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <sstream>

#include <QJSEngine>
#include <QJSValueIterator>
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
//       set_option(OptionName, double(i));
// }
// void IView3d::set_option(QString const &OptionName, int i) {
//    set_option(OptionName, double(i));
// }


uint32_t hsv_uint(float h, float s, float v) {
   return 0xff000000 | (uint32_t)ct::Hsv(h,s,v).uint32();
}

uint32_t ColFromAlpha(float a);

uint32_t hsva_uint(float h, float s, float v, float a) {
   return ColFromAlpha(a) | (uint32_t)ct::Hsv(h,s,v).uint32();
}

static QString s_format_i(QString const &Format, int val) {
   return QString::asprintf(Format.toLatin1().constData(), val);
}

static QString s_format_f(QString const &Format, double val) {
   return QString::asprintf(Format.toLatin1().constData(), val);
}

static QString s_format_s(QString const &Format, QString const &val) {
   return QString::asprintf(Format.toLatin1().constData(), val.toUtf8().constData());
}


// stringify a QVariant which is known to be a list (originates from a JS array
// argument marshalled through a QVariantList-typed Q_INVOKABLE parameter).
static QString FormatArray(QVariant const &ArrayValue, bool AddQuotesToStrs, QString Left, QString Separator, QString Right)
{
   QString
      sout;
   QTextStream
      str(&sout);
   str << Left;
   QVariantList
      L = ArrayValue.toList();
   for (int i = 0; i < L.size(); ++ i) {
      if (i != 0) {
         str << Separator;
      }
      QVariant const
         &v = L[i];
      if (v.typeId() == QMetaType::QString && AddQuotesToStrs) {
         str << "\"" << v.toString() << "\"";
      } else if (v.typeId() == QMetaType::QVariantList) {
         str << FormatArray(v, AddQuotesToStrs, Left, Separator, Right);
      } else {
         str << v.toString();
      }
   }
   str << Right;
   str.flush();
   return sout;
}


static bool IsNumericVariant(QVariant const &v) {
   switch (v.typeId()) {
      case QMetaType::Int: case QMetaType::UInt: case QMetaType::Double:
      case QMetaType::LongLong: case QMetaType::ULongLong: case QMetaType::Float:
         return true;
      default:
         return false;
   }
}


// core of the script-visible format(baseStr, args...)/"baseStr".format(args...) function.
static QString Invoke_format(QString const &s0, QVariantList const &args)
{
   QString
      s = s0;
   int
      nArgs = args.size();
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
      QVariant const &arg = args[i];
      if (IsNumericVariant(arg)) {
         double q = arg.toDouble();
         // ECMA script does not distinguish between floats and integers...
         if (q == long(q)) {
            ADD_ARG(long(q));
         } else {
            ADD_ARG(double(q));
         }
      } else if (arg.typeId() == QMetaType::Bool) {
         ADD_ARG(bool(arg.toBool()));
      } else if (arg.typeId() == QMetaType::QString) {
         StrBuf.push_back(arg.toString().toStdString());
         ADD_ARG(StrBuf.back());
      } else if (arg.typeId() == QMetaType::QVariantList) {
         StrBuf.push_back(FormatArray(arg,true,"[",", ","]").toStdString());
         ADD_ARG(StrBuf.back());
      } else {
         StrBuf.push_back(arg.toString().toStdString());
         ADD_ARG(StrBuf.back());
      }
   }
#undef ADD_ARG
   fmt::ArgList
      fargs(Types, &Args[0]);
   w.write(s.toStdString(), fargs);
   s = QString::fromStdString(w.str());
   return s;
}


// core of the script-visible join(sepStr, array)/"sepStr".join(array) function.
static QString Invoke_join(QString const &s, QVariantList const &args)
{
   if (args.size() != 1 || args[0].typeId() != QMetaType::QVariantList) {
      IvNotify(NOTIFY_Error, "join() called with unsupported arguments");
      return QString();
   }
   return FormatArray(args[0], false, "", s, "");
}


// convert a JS array (marshalled as QVariantList) or a Vec3-like JS object
// (marshalled as QVariantMap, see ToScriptValue() below) into a FVec3d.
static FVec3d FromVariant(QVariant const &v)
{
   FVec3d
      r(0.,0.,0.);
   if (v.typeId() == QMetaType::QVariantList) {
      QVariantList L = v.toList();
      for (int i = 0; i < (int)r.size() && i < L.size(); ++ i)
         r[i] = L[i].toDouble();
   } else if (v.typeId() == QMetaType::QVariantMap) {
      QVariantMap m = v.toMap();
      if (m.contains("x")) r[0] = m.value("x").toDouble();
      if (m.contains("y")) r[1] = m.value("y").toDouble();
      if (m.contains("z")) r[2] = m.value("z").toDouble();
   }
   return r;
}


Q_DECLARE_METATYPE(QSharedPointer<FFreeLine>)

void _DebugPrintScriptValueProperties(QString desc, QJSValue obj, bool DataMembersOnly=true) {
   QJSValueIterator it(obj);
   while (it.hasNext()) {
      it.next();
      QJSValue val = it.value();
      if (DataMembersOnly && !val.isCallable())
         IvEmit("  [%1] %2 (fn? %4) = %5", desc, it.name(), val.isCallable(), val.toString());
   }
}


struct FVariantFmt {
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
      if (variant.canConvert<QVariantMap>()) {
         _DebugFmt(ss, desc, variant.value<QVariantMap>(), level+1);
      } else if (variant.canConvert<QVariantList>()) {
         _DebugFmt(ss, desc, variant.value<QVariantList>(), level+1);
      } else if (variant.canConvert<FVec3d>()) {
         FVec3d const &v = variant.value<FVec3d>();
         ss << "Vec3(" << v[0] << "," << v[1] << "," << v[2] << ")";
      } else if (variant.typeId() == QMetaType::QString) {
         // stored as actual string. Print it with quotes (should probably add escapes...)
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


// QJSEngine has no live-prototype mechanism for a custom C++ value type, so build
// a fresh plain JS object with x/y/z (and numeric-index) properties on every crossing.
static QJSValue ToScriptValue(QJSEngine &Engine, FVec3d const &v)
{
   QJSValue
      obj = Engine.newObject();
   for (size_t i = 0; i < v.size(); ++ i)
      obj.setProperty(quint32(i), v[i]);
   obj.setProperty("x", v[0]);
   obj.setProperty("y", v[1]);
   obj.setProperty("z", v[2]);
   obj.setProperty("length", int(v.size()));
   obj.setProperty("toString", Engine.globalObject().property("__vec3_toString"));
   return obj;
}


// Exposes free functions/constructors callable from script as Q_INVOKABLE methods.
// Instantiated once per ExecScript() call and installed under a hidden global
// name; the actual script-visible names (format, join, Vec3, FreeLine, ...) are
// either aliased directly onto the global object (fixed-arity functions) or
// reached via small JS wrapper functions set up in ExecScript() (variable-arity
// functions, which need to gather the JS 'arguments' object into a list first).
class FScriptGlobals : public QObject
{
   Q_OBJECT
   QJSEngine
      *m_pEngine;
public:
   explicit FScriptGlobals(QJSEngine *pEngine, QObject *parent = 0)
      : QObject(parent), m_pEngine(pEngine)
   {}

   Q_INVOKABLE quint32 hsva(float h, float s, float v, float a) const { return hsva_uint(h,s,v,a); }
   Q_INVOKABLE quint32 irgb(quint32 c) const { return ct::irgb(c); }
   Q_INVOKABLE QString fmti(QString const &Format, int val) const { return s_format_i(Format, val); }
   Q_INVOKABLE QString fmtf(QString const &Format, double val) const { return s_format_f(Format, val); }
   Q_INVOKABLE QString fmts(QString const &Format, QString const &val) const { return s_format_s(Format, val); }
   Q_INVOKABLE QString replace_ext(QString const &FileName, QString const &NewExt) const { return ReplaceExt(FileName, NewExt); }
   Q_INVOKABLE QString remove_path(QString const &FileName) const { return RemovePath(FileName); }

   Q_INVOKABLE QString format_impl(QString const &s, QVariantList const &args) const { return Invoke_format(s, args); }
   Q_INVOKABLE QString join_impl(QString const &s, QVariantList const &args) const { return Invoke_join(s, args); }

   Q_INVOKABLE QJSValue vec3_impl(QVariantList const &args) const
   {
      FVec3d v(0.,0.,0.);
      if (args.size() == 3) {
         for (int i = 0; i < 3; ++ i)
            v[i] = args[i].toDouble();
      } else if (args.size() == 1) {
         v = FromVariant(args[0]);
      } else {
         return m_pEngine->newErrorObject(QJSValue::TypeError, "Vec3(...): expected components or list of components as arguments.");
      }
      return ToScriptValue(*m_pEngine, v);
   }

   Q_INVOKABLE QJSValue freeline_impl(QVariantList const &args) const
   {
      int nArgs = args.size();
      if (nArgs < 2)
         return m_pEngine->newErrorObject(QJSValue::RangeError, "FreeLine.ctor(): expected at least two arguments (from,to)");
      FFreeLine *pobj(new FFreeLine(m_pEngine));
      int iArg = 0;
      if (nArgs == 1) {
         pobj->setp(args[iArg].toMap());
      } else if (nArgs >= 3) {
         if (nArgs > iArg) { pobj->SetFrom(FromVariant(args[iArg])); iArg += 1; }
         if (nArgs > iArg) { pobj->SetTo(FromVariant(args[iArg])); iArg += 1; }
         if (nArgs > iArg) { pobj->setp(args[iArg].toMap()); iArg += 1; }
      }
      // let the engine manage the new object's lifetime.
      QJSEngine::setObjectOwnership(pobj, QJSEngine::JavaScriptOwnership);
      return m_pEngine->newQObject(pobj);
   }
};


static void AddScriptFunction(QJSEngine &ScriptEngine, QJSValue &GlobalsObject, QString const &Name)
{
   ScriptEngine.globalObject().setProperty(Name, GlobalsObject.property(Name));
}

void ExecScript(IApplication *app, IView3d *view, QString const &ScriptText, QString const &ScriptName) {
   IvEmit("--- exec script '%1'", ScriptName);
   QJSEngine
      ScriptEngine;

   ScriptEngine.globalObject().setProperty("app", ScriptEngine.newQObject(app));
   // FIXME: added IApplication as both "app" and "doc" until things get sorted out.
   ScriptEngine.globalObject().setProperty("doc", ScriptEngine.newQObject(app));
   ScriptEngine.globalObject().setProperty("view", ScriptEngine.newQObject(view));
   // 'app'/'view' are externally owned, long-lived objects (they outlive this
   // function and this per-call ScriptEngine). QJSEngine::newQObject() defaults
   // to taking JavaScript ownership of any wrapped QObject that has no parent,
   // which for 'app' (the top-level main window, parent-less by design) means
   // its garbage collector would delete the main window itself once this
   // ScriptEngine is destroyed below. Force C++ ownership to prevent that.
   QJSEngine::setObjectOwnership(app, QJSEngine::CppOwnership);
   QJSEngine::setObjectOwnership(view, QJSEngine::CppOwnership);

   // note: the FScriptGlobals instance is owned by the engine (see setObjectOwnership
   // below), and must be constructed after the engine, so it is fine to reference it here.
   FScriptGlobals
      *pGlobals = new FScriptGlobals(&ScriptEngine);
   QJSEngine::setObjectOwnership(pGlobals, QJSEngine::JavaScriptOwnership);
   QJSValue
      GlobalsObject = ScriptEngine.newQObject(pGlobals);
   ScriptEngine.globalObject().setProperty("__script_globals", GlobalsObject);

   // fixed-arity functions: can be aliased onto the global object directly.
   AddScriptFunction(ScriptEngine, GlobalsObject, "hsva");
   AddScriptFunction(ScriptEngine, GlobalsObject, "irgb");
   AddScriptFunction(ScriptEngine, GlobalsObject, "fmti");
   AddScriptFunction(ScriptEngine, GlobalsObject, "fmtf");
   AddScriptFunction(ScriptEngine, GlobalsObject, "fmts");
   AddScriptFunction(ScriptEngine, GlobalsObject, "replace_ext"); // replace a file extension by another file extension
   AddScriptFunction(ScriptEngine, GlobalsObject, "remove_path"); // strip off path from a file name

   // variable-arity functions/constructors: need small JS-side wrappers to gather
   // the call's 'arguments' into an array before forwarding to the Q_INVOKABLE
   // implementation (QJSEngine has no way to bind a C++ callback directly as a
   // variadic script function). This keeps the call syntax ("format(s, a, b)",
   // "s.format(a, b)", "Vec3(x,y,z)", ...) unchanged for existing scripts.
   QJSValue BootstrapResult = ScriptEngine.evaluate(
      "function format(s) { return __script_globals.format_impl(s, Array.prototype.slice.call(arguments, 1)); }\n"
      "String.prototype.format = function() { return __script_globals.format_impl(this, Array.prototype.slice.call(arguments)); };\n"
      "function join(s) { return __script_globals.join_impl(s, Array.prototype.slice.call(arguments, 1)); }\n"
      "String.prototype.join = function() { return __script_globals.join_impl(this, Array.prototype.slice.call(arguments)); };\n"
      "function __vec3_toString() { return 'Vec3(' + this.x + ',' + this.y + ',' + this.z + ')'; }\n"
      "function Vec3() { return __script_globals.vec3_impl(Array.prototype.slice.call(arguments)); }\n"
      "function FreeLine() { return __script_globals.freeline_impl(Array.prototype.slice.call(arguments)); }\n",
      "<script-bootstrap>");
   if (BootstrapResult.isError()) {
      IvNotify(NOTIFY_Error, "Internal error setting up script engine bootstrap functions:\n" + BootstrapResult.toString());
      return;
   }

   QJSValue
      res = ScriptEngine.evaluate(ScriptText, ScriptName);
   if (res.isError()) {
      IvNotify(NOTIFY_Error, "Error during script execution:\n" + QString("%1:%2: %3\n").arg(ScriptName, QString::number(res.property("lineNumber").toInt()), res.toString()));
   }
   IvEmit("--- script terminated.");
}

QString LoadTextFileViaQt(QString const &FileName)
{
   QFile
      File(FileName);
   if (!File.open(QIODevice::ReadOnly | QIODevice::Text)) {
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


void ExecScript(IApplication *app, IView3d *view, QString const &FileName)
{
   QString
      ScriptText = LoadTextFileViaQt(FileName);
   ExecScript(app, view, ScriptText, FileName);
}

#include "IvScript.moc"
