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

#ifndef IV_DOCUMENT_OPTIONS_H
#define IV_DOCUMENT_OPTIONS_H

#include <string>
#include <sstream>

#include <QString>
#include <QObject>

#include <QExplicitlySharedDataPointer>
#include <QSharedData>
#include <QSharedPointer>
#include <QVariant>
#include <QMetaProperty>

#include "CxVec3.h"
#include "format.h"
typedef ct::TVector3<double>
   FVec3d;


namespace ct {
   struct FHfOptions;
   struct FWfDecl;
   struct FAtomSet;
}
// using ct::FAtomSetPtr;
// using ct::FBasisSetPtr;
// using ct::TArray;
// using ct::FIntrusivePtrDest;

class FView3d;
class FDocument;
class QTextStream;

void IvWarn(std::string const &Text);
void IvWarn(QString const &Text);
void IvEmit(QString const &Text);
std::string q2s(QString const &s);
QString s2q(std::string const &s);
QTextStream &operator << (QTextStream &str, FVec3d const &v);


char const *BoolToCstr(bool o);
Q_DECLARE_METATYPE(FVec3d)
Q_DECLARE_METATYPE(FVec3d*)


// constructs a FVec3d object from either a variant already holding a FVec3d object,
// or from a variant list type (from conversion of script arrays)
FVec3d FVec3d_FromVariant(QVariant const &v);


enum FOptionsDescType {
   // Generate list of property assignments for a concretely known named object. E.g.,
   //
   //    view.atom_scale = 75.00000;
   //    view.depth_peeling_layers = 8;
   //
   OPTIONSDESC_PropertyAssign,
   // Generate ECMA-script object literal (cf Python dict) for the non-default properties
   //
   //    {atom_scale: 75.00000, depth_peeling_layers: 8}
   //
   OPTIONSDESC_ObjectLiteral,
   // Code for constructing an object. By default, specifying a corresponding script class name,
   // this will generate a constructor call with the object literal as argument:
   //
   //    FreeLine({from: [0.,0.,0.], to: [1.,0.,0.], color: 0xff3f3f3f})
   OPTIONSDESC_ObjectCtor
};

// helper class for generating script state representations for property objects
struct FOptionsDesc {
   // TODO: it would probably be a way better idea to just create the text in
   // QStrings to begin with, rather than converting it to std::string, and,
   // most likely, back at the end.
   explicit FOptionsDesc();
   explicit FOptionsDesc(FOptionsDescType Type, char const *pScriptClassName, char const *pScriptVarName);

   void start(FOptionsDescType Type, char const *pScriptClassName, char const *pScriptVarName);
   void finish();

#ifdef INCLUDE_ABENDONED
//    template<class T>
//    void setp(char const *PropName_, char const *ValueFmt_, T const &PropValue) {
//       if (m_Finished)
//          throw std::runtime_error("FOptionsDesc: attempted to modify this object after it was already finalized.");
//       if (m_nAssigned != 0)
//          m_str << m_Delim;
//       m_str << fmt::format(m_PropFmt, fmt::format(ValueFmt_, PropValue));
//       m_nAssigned += 1;
//    };
#endif // INCLUDE_ABENDONED
   // PropIndex: index of the property in case of fixed-order properties (e.g., as c'tor arguments).
   // Use -1 if the property does not have an index. If indices are used, they must be consecutive in
   // order 0,1,..,(max),-1, and after the first -1, no other indices may follow.
   void setp1(int PropIndex, char const *PropName_, std::string const &PropValue) {
      if (m_Finished)
         throw std::runtime_error(fmt::format("FOptionsDesc: attempted to modify this object (set .{} := {}) after it was already finalized.", PropName_, PropValue));
      if (PropIndex != -1 && PropIndex != m_NextPositionalArg)
         throw std::runtime_error(fmt::format("FOptionsDesc: attempted to modify this object (set .{} := {}) after it was already finalized.", PropName_, PropValue));
      if (m_nAssigned != 0)
         m_str << m_Delim;
      if (PropIndex < 0) {
         if (m_NextPositionalArg != -1) {
            // first named argument. no more positional arguments allowed.
            m_NextPositionalArg = -1;
            m_str << m_OpenNamed;
            m_nAssigned = 0;
         }
         PropIndex = -1;
      } else {
         assert(PropIndex == m_NextPositionalArg);
         m_NextPositionalArg = PropIndex + 1;
      }
      if (PropIndex >= 0 && !m_PropFmtPositional.empty())
         m_str << fmt::format(m_PropFmtPositional, PropName_, PropValue);
      else
         m_str << fmt::format(m_PropFmtNamed, PropName_, PropValue);
      m_nAssigned += 1;
      (void)PropIndex; // suppress unusued warning
   };
   // ^-- TODO: implement CtorArgIndex?

   // set property using default formatting (Format1), applying it recursively on complex objects.
   template<class T>
   void setp(int PropIndex, char const *PropName_, T const &PropValue) {
      setp1(PropIndex, PropName_, this->Format1(PropValue));
   };

   // set simple property using explicit fmt::format() format string.
   template<class T>
   void setp(int PropIndex, char const *PropName_, char const *ValueFmt_, T const &PropValue) {
      setp1(PropIndex, PropName_, fmt::format(ValueFmt_, PropValue));
   };

//    std::string str() const;
   std::string str();
protected:
   std::string
      m_PropFmtNamed, m_PropFmtPositional,
      m_Open, m_Delim, m_Close,
      m_OpenNamed, m_CloseNamed;
   std::stringstream
      m_str;
   size_t
      m_nAssigned;
   int
      m_NextPositionalArg;
   bool
      m_Finished;

   void _InitData();
public:
   // default formatting routines for creating js-code representations of objects.
   // Note: one can overwrite those by directly formatting the property as string.

   std::string Format1(unsigned long dw) const { return fmt::format("0x{:08x}", dw); } // <-- mostly used for colors in our context
   std::string Format1(unsigned int dw) const { return fmt::format("0x{:08x}", dw); } // <-- mostly used for colors in our context
//    std::string Format1(uint32_t dw) const { return fmt::format("0x{:08x}", dw); } // <-- mostly used for colors in our context
   std::string Format1(int d) const { return fmt::format("{}", d); }
   std::string Format1(double d) const { return fmt::format("{:.8g}", d); }
   std::string Format1(bool d) const { return BoolToCstr(d); }
   std::string Format1(QString const &s) const { return s.toStdString(); }
   std::string Format1(std::string const &s) const { return s; }

   template<class T>
   std::string Format1(ct::TVector3<T> const &v) const {
      return fmt::format("Vec3({},{},{})", Format1(v[0]), Format1(v[1]), Format1(v[2]));
   }

   std::string Format1(QObject *obj) const {
      // check if this has a toString slot we could invoke to format it. Script objects should have those.
      QString sRetVal;
//       Q_ARG(int, param)
      if (!QMetaObject::invokeMethod(obj, "toString", Qt::DirectConnection, Q_RETURN_ARG(QString, sRetVal))) {
         IvWarn(QString("FOptionsDesc: failed to transcribe object '%1' of type '%2'.").arg(obj->objectName()).arg(obj->metaObject()->className()));
         return "null";
      }
      return sRetVal.toStdString();
   }
   std::string Format1(QObject &obj) const { return Format1(&obj); }

   template<class K, class V>
   std::string Format1(QMap<K,V> &qm) const {
      std::stringstream ss;
      ss << "{";
      for (typename QMap<K,V>::iterator it = qm.begin(); it != qm.end(); ++it) {
         if (it != qm.begin()) ss << ", ";
         ss << Format1(it.key()) << ": " << Format1(it.value());
      }
      ss << "}";
      return ss.str();
   }
};


class FWfOptions : public QObject
{
   Q_OBJECT
public:
#include "prop_FWfOptions.h.inl"

public:
   FWfOptions(QObject *parent = 0);

   // assign basis sets from HfOptions to pAtoms
   void AssignBasisSets(ct::FAtomSet *pAtoms);

   void AssignScfOptions(ct::FHfOptions &HfOptions);
   void AssignWfDecl(ct::FWfDecl &WfDecl, ct::FAtomSet *pAtoms);
};



int const ColorNotSet = 0x04030201; // if you are thinking of using this: Use black instead 8).

// this class doubles as element properties and per-atom settings container.
// It is used for exposing these settings to user modification and to scripts.
// It is a bit more complicated than the other property classes, as all the properties
// just default to "not set"...
class FElementOptions : public QObject, public QSharedData {
   Q_OBJECT
public:
#include "prop_FElementOptions.h.inl"

   explicit FElementOptions(int iElement_, QObject *parent=0);
   explicit FElementOptions(FElementOptions const *other, QObject *parent=0);

   int iElement() const { return m_iElement; }
   // return the bond color if explicitly set, and otherwise the atom color.
   uint32_t GetBondColor1() const;

   // -1: no override (these things are used in IRC and alignment weighting. Currently
   // they cannot be set from script.).
   double GetMassOverride() const { return -1.; }
   double GetNuclearChargeOverride() const { return -1.; }
//    double GetVdwRadius() const;
protected:
   // these ones get stuff from tables.
   uint32_t GetDefaultColor() const;
   uint32_t GetDefaultBondColor() const;
   double GetDefaultCovalentRadius() const;
   double GetDefaultDrawRadius() const;
   double GetDefaultVdwRadius() const;
   char const *ElementName() const;

   int
      m_iElement;
   void InitForElement(int iElement_);
};
typedef QExplicitlySharedDataPointer<FElementOptions>
   FElementOptionsPtr;
typedef QList<FElementOptionsPtr>
   FElementOptionsList;

struct FAtomOptions {
   uint
      // bitfield of FAtomFlag type. (ATOM_*)
      Flags;
   int64_t
      // sequence number in order to keep track of order in which things were selected
      // (important for bond angle and dihedral measures)
      iSelectionSequenceId;
   FElementOptionsPtr
      // used to provide values to the atom which differ from the default element properties.
      pPropertiesOverride;

   FAtomOptions() : Flags(0), iSelectionSequenceId(-1) {}
};


// Some sort of additional geometry object which can be added to the scene
class FFreeObject : public QObject
{
   Q_OBJECT
public:
   FFreeObject(QObject *parent);
   typedef QSharedPointer<FFreeObject> FNewFreeObjectPtr;

   virtual FNewFreeObjectPtr newLinkedObject(QObject *parent) = 0;
   virtual FNewFreeObjectPtr newClonedObject(QObject *parent) const = 0;
public slots:
   virtual QString toString() const;
};
// typedef FFreeObject
//    *FFreeObjectPtr;
typedef QSharedPointer<FFreeObject>
   FFreeObjectPtr;
// ^--
// The actual objects are split into:
// - this QObject derived class, with the actual interfaces. However, it really only
//   has a QExplicitlySharedDataPointer to the data part as class member.
//   That is the object we access in scripts, and/or link to UI interfaces.
// - a Data part (derived from QSharedData),
//   which is copyable/cloneable and used with either implicit copy-on-write (QSharedPointer)
//   or explicit copy-on-write (QExplicitlySharedDataPointer; using detach() manually if needed);
//
// So while the objects are all individual, their actual data is/may be shared.
// I think this should allow cleanly taking care of the ownership issues of script objects
// vs UI objects vs serialized objects:
//
// - e.g., script objects belong to the script engine, and will be destroyed with it.
//   however, their actual data members can be shared with the FDocument / UI versions
//   of those objects; in particular, they can still be modified via script
//
// See info here on how QSharedData is meant to be used:
//
// - https://doc.qt.io/qt-5/implicit-sharing.html
// - https://doc.qt.io/qt-5/qexplicitlyshareddatapointer.html
//
// TODO:
// - If we're updating make_properties.py anyway, we might also want to add
//   automatic constructors from a QScriptObject or QVariantMap parameter. And
//   code to generate something equivalent to the QVariantMap to represent the
//   changed parameters.

Q_DECLARE_METATYPE(FFreeObject*)

// A straight line or arrow between two points in 3D space
class FFreeLine : public FFreeObject
{
   Q_OBJECT
public:
#include "prop_FFreeLine.h.inl"

public slots:
   QString toString() const;

public:
//    FFreeLine(QObject *parent = 0);
   FFreeLine(FVec3d vFrom, FVec3d vTo, uint32_t dwColor = 0xff3f3f3f, double fWidth = 0.02, double fWeight = 1., QObject *parent = 0);
};

// A free, compact text label which can be attached to objects in 3D space
class FFreeLabel : public FFreeObject
{
   Q_OBJECT
public:
#include "prop_FFreeLabel.h.inl"

public slots:
   QString toString() const;

public:
   FFreeLabel(QObject *parent = 0);
   FFreeLabel(QString Text, FVec3d vPos, double fSize, uint32_t dwColor, QObject *parent = 0);
};

// don't ask.
class FInvalidFreeObject : public FFreeObject
{
   Q_OBJECT
public:
   FInvalidFreeObject(QObject *parent);
   FNewFreeObjectPtr newLinkedObject(QObject *parent); // override
   FNewFreeObjectPtr newClonedObject(QObject *parent) const; // override
};

Q_DECLARE_METATYPE(uint32_t);

void _RegisterMetaTypes();

#endif // IV_DOCUMENT_OPTIONS_H
