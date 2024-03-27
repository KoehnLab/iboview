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

// TODO:
// - in make_properties.py, add generator for c'tor from QScriptEngine /
//   QScriptContext and setp(QScriptValue) (optimally both in the positional
//   arguments and named arguments fashion). Currently there is quite a bit of
//   code duplication between here and IvScript.cpp. And if we mean to extend the
//   scripting capabilities, there will be many more similar applications.

#include <QSettings>
#include <QVariant>
#include <QTextStream>

#include "IvDataOptions.h"
#include "Iv.h"

#include "CxOpenMpProxy.h"
#include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "CtRhf.h"


FVec3d FVec3d_FromVariant(QVariant const &variant) {
   for ( ; ; ) {
      if (variant.canConvert<FVec3d>()) {
         // that already is a FVec3d. just return it as-is.
         return variant.value<FVec3d>();
      }

      if (variant.canConvert<QVariantList>()) {
         // if it's a length-3 array of some numeric types, that is
         // also fine for making a vector out of it.
         QVariantList ql = variant.value<QVariantList>();
         FVec3d v;
         if (size_t(ql.size()) != size_t(v.size()))
            break;
         for (size_t i = 0; i < v.size(); ++ i) {
            bool ok = false;
            v[i] = ql[i].toDouble(&ok);
            if (!ok) break;
         }
         return v;
      }

      // not an actual vector or a number list. Don't know what to do with this at this moment.
      break;
   }
   IvNotify(NOTIFY_Error, IvFmt("failed to convert input variant '%1' to FVec3d", variant.toString()));
   return FVec3d(0,0,0);
}


void FOptionsDesc::_InitData() {
   m_nAssigned = 0;
   m_Finished = false;
   m_NextPositionalArg = 0;
}

FOptionsDesc::FOptionsDesc() {
   _InitData();
}

FOptionsDesc::FOptionsDesc(FOptionsDescType Type, char const *pScriptClassName, char const *pScriptVarName) {
   _InitData();
   start(Type, pScriptClassName, pScriptVarName);
}

void FOptionsDesc::start(FOptionsDescType Type, char const *pScriptClassName, char const *pScriptVarName) {
   m_OpenNamed = std::string();
   m_CloseNamed = std::string();
   // unless explicitly set: positional arguments not allowed.
   m_PropFmtPositional = std::string();
   if (Type == OPTIONSDESC_ObjectLiteral || Type == OPTIONSDESC_ObjectCtor) {
      if (bool(pScriptClassName) && Type == OPTIONSDESC_ObjectCtor) {
         m_Open = fmt::format("{}(", pScriptClassName);
         m_Close = ")";
         m_OpenNamed = "{";
         m_CloseNamed = "}";
      } else {
         m_Open = "{";
         m_Close = "}";
      }
      m_Delim = ", ";
      m_PropFmtNamed = "{}: {}";
      if (Type == OPTIONSDESC_ObjectCtor)
         // ignore name, just format value
         m_PropFmtPositional = "{1}";
   } else if (Type == OPTIONSDESC_PropertyAssign) {
      if (!bool(pScriptVarName))
         pScriptVarName = "";
      m_Open = "";
      m_Close = "";
      m_Delim = "\n";
      m_PropFmtNamed = fmt::format("{}.{{}} = {{}};", pScriptVarName);
   } else {
      throw std::runtime_error(fmt::format("FOptionsDesc: Desc type {} not recognized.", int(Type)));
   }
   m_str.str(std::string());
   m_str << m_Open;
};


void FOptionsDesc::finish() {
   if (!m_Finished) {
      m_Finished = true;
      if (m_NextPositionalArg == -1 && !m_CloseNamed.empty())
         m_str << m_CloseNamed;
      m_str << m_Close;
   }
}


std::string FOptionsDesc::str() {
   this->finish();
   return m_str.str();
}



extern int g_nMaxOmpThreads;

FWfOptions::FWfOptions(QObject *parent)
   : QObject(parent)
{
   setObjectName("WfOptions");
   InitProperties();
   m_NumThreads = g_nMaxOmpThreads;
   if (1) {
      // check if we have a global override for the memory stack size.
      QSettings settings;
      QVariant v = settings.value("IboView/MemoryStackSize");
      bool ok;
      int MemOverride = v.toInt(&ok);
      if (ok && MemOverride > 0)
         m_WorkSpaceMb = MemOverride;

   }
}

void FWfOptions::AssignBasisSets(ct::FAtomSet *pAtoms)
{
   std::string
      sOrbBasis = q2s(GetOrbBasis()),
      sFitBasis = q2s(GetFitBasis());
//    IvEmit(s2q("// SET BASIS: " + sOrbBasis + " // " + sFitBasis));
//    for (size_t iAt = 0; iAt != pAtoms->size(); ++ iAt) {
//       (*pAtoms)[iAt].BasisDesc[ct::BASIS_Orbital] = sOrbBasis;
//       (*pAtoms)[iAt].BasisDesc[ct::BASIS_JFit] = sFitBasis;
//       (*pAtoms)[iAt].BasisDesc[ct::BASIS_JkFit] = sFitBasis;
//    }
   pAtoms->AssignBasis(ct::BASIS_Orbital, sOrbBasis, true); // true: AutoFix=on.
   // hm... should I really keep these? I mean the explicit assignments of the fitting bases.
   if (!(sFitBasis.empty() || sFitBasis == "auto")) {
      pAtoms->AssignBasis(ct::BASIS_JFit, sFitBasis);
      pAtoms->AssignBasis(ct::BASIS_JkFit, sFitBasis);
   }
}

void FWfOptions::AssignScfOptions(ct::FHfOptions &ScfOpt)
{
   ScfOpt.SetGridDesc("1e-4");
   ScfOpt.ThrOrb = GetThrGrad();
   ScfOpt.ThrDen = GetThrEnergy();
   ScfOpt.MaxIt = (uint)GetMaxIter();
   ScfOpt.XcAlgo = ct::XCALGO_AuxiliaryExpand;
   ScfOpt.JkAlgo = ct::JKALGO_DfCoulombOnlyCache3ix;
   ScfOpt.ComputeFlags = 0;

   // set xc functional first... otherwise, if SCF command is RHF,
   // this will run a DFJK-RKS with exchange factor 0.0 and actually
   // compute the functional...
   {
      // remove the explanation parts of the functional so that
      // FXcFunctional has a chance to recognize it.
      QString
         sFunctional1 = GetFunctional();
      int iParens = sFunctional1.indexOf("(");
      if (iParens != -1) {
         sFunctional1.chop(sFunctional1.size() - iParens);
         sFunctional1 = sFunctional1.trimmed();
      }
   //    IvEmit("Set xc: '%1'",sFunctional1);
//       ScfOpt.XcFunctionalName = q2s(sFunctional1);
      ScfOpt.SetXc(q2s(sFunctional1));
   }

   {
      // pass *only* the explanation part to SetArgs...
      QString
         sScfMethod = GetScfMethod();
      int iParens0 = sScfMethod.indexOf("(");
      if (iParens0 != -1) {
         int iParens1 = sScfMethod.indexOf(")", iParens0);
//          QStringRef sScfCommand(&sScfMethod, iParens0+1, iParens1-iParens0-1);
         QString sScfCommand = sScfMethod.mid(iParens0+1, iParens1-iParens0-1);
         ScfOpt.SetArgs(q2s(sScfCommand));
//          IV_NOTIFY2(NOTIFY_Warning, "Set SCF Command '%1' xc: '%2'", sScfCommand, s2q(ScfOpt.XcFunctionalName));
      }
   }
//    if (sScfMethod == "Kohn-Sham (DFJX-RKS)" || sScfMethod == "Kohn-Sham (DF-RKS)") {
//       ScfOpt.XcAlgo = ct::XCALGO_AuxiliaryExpand;
//       ScfOpt.JkAlgo = ct::JKALGO_DfCoulombOnlyCache3ix;
//    } else if (sScfMethod == "Kohn-Sham (DFJ-RKS, cached)") {
//       ScfOpt.XcAlgo = ct::XCALGO_Regular;
//       ScfOpt.JkAlgo = ct::JKALGO_DfCoulombOnlyCache3ix;
//    } else if (sScfMethod == "Kohn-Sham (DFJ-RKS, direct)") {
//       ScfOpt.XcAlgo = ct::XCALGO_Regular;
//       ScfOpt.JkAlgo = ct::JKALGO_DfCoulombOnly;
//    } else if (sScfMethod == "Hartree-Fock (DF-RHF)") {
//       ScfOpt.XcAlgo = ct::XCALGO_Regular;
//       ScfOpt.JkAlgo = ct::JKALGO_DfCoulombAndExchange;
//    } else if (sScfMethod == "Kohn-Sham/Hybrid (DFJK-RKS)") {
//       ScfOpt.XcAlgo = ct::XCALGO_Regular;
//       ScfOpt.JkAlgo = ct::JKALGO_DfCoulombAndExchange;
//    }

   // grid params?
   {
      // we provide the ability to set free options. But since there is a chance this won't
      // work, we wrap it
      QString
         sScfExtraOptions = GetScfOptions();
//       try {
         ScfOpt.SetArgs(q2s(sScfExtraOptions));
//       } catch (ct::FInputError const &e) {
//          IvNotify(NOTIFY_Error, IvFmt("while processing free scf option input:\n'%1':\n%2", sScfExtraOptions, s2q(e.what())));
//          return false;
//       }
   }
//    return true;
}

void FWfOptions::AssignWfDecl(ct::FWfDecl &WfDecl, ct::FAtomSet *pAtoms)
{
   int
      iCharge = GetCharge(),
      iExtraSpin = GetExtraSpin();
   int
      nElec = pAtoms->NuclearCharge() - iCharge,
      Ms2;
   if (nElec % 2 == 0)
      Ms2 = 0 + 2 * iExtraSpin;
   else
      Ms2 = 1 + 2 * iExtraSpin;
   WfDecl = ct::FWfDecl(iCharge, Ms2);
   WfDecl.SetNuclearCharge(pAtoms->NuclearCharge());
   if (nElec != int(WfDecl.nElec()))
      IV_NOTIFY2(NOTIFY_Error, "Something went wrong in WfDecl assignment nElec=%1 vs WfDecl.nElec=%2", nElec, int(WfDecl.nElec()));
}




void FElementOptions::InitForElement(int iElement_)
{
   m_iElement = iElement_;
   InitProperties(); // must come AFTER m_iElement assignment.
}

FElementOptions::FElementOptions(int iElement_, QObject *parent)
   : QObject(parent)
{
   InitForElement(iElement_); // must come AFTER m_iElement assignment.
   assert(parent == 0);
   // ^- intended to be used with smart pointers, due to
   //    explicit sharing. Not with ownership hierarchies.
}

FElementOptions::FElementOptions(FElementOptions const *other, QObject *parent)
   : QObject(parent)
{
   InitForElement(other->iElement()); // must come AFTER m_iElement assignment.
   CopyPropertiesFrom(*other);
   assert(parent == 0);
   // ^- intended to be used with smart pointers, due to
   //    explicit sharing. Not with ownership hierarchies.
}


// that's the rasmol CPKnew colors, according to jmol homepage.
// See MakeColorCodes.py
// uint32_t ElementColors[109] = {0xffffff,0xffc0cb,0xb22121,0xff1493,0x00ff00,0xd3d3d3,0x87cee6,0xff0000,0xdaa520,0xff1493,0x0000ff,0x228b22,0x696969,0xdaa520,0xffaa00,0xffff00,0x00ff00,0xff1493,0xff1493,0x696969,0xff1493,0x696969,0xff1493,0x696969,0x696969,0xffaa00,0xff1493,0x802828,0x802828,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0x696969,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xa020f0,0xff1493,0xff1493,0xffaa00,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xdaa520,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493};

// as above, but replaced carbon by 0x999999 (0.6/0.6/0.6)
// also inserted a number for dummies (element 0).
static uint32_t ElementColors[110] = {0x404040, 0xffffff,0xffc0cb,0xb22121,0xff1493,0x00ff00,0x999999,0x87cee6,0xff0000,0xdaa520,0xff1493,0x0000ff,0x228b22,0x696969,0xdaa520,0xffaa00,0xffff00,0x00ff00,0xff1493,0xff1493,0x696969,0xff1493,0x696969,0xff1493,0x696969,0x696969,0xffaa00,0xff1493,0x802828,0x802828,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0x802828,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0x696969,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xa020f0,0xff1493,0xff1493,0xffaa00,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xdaa520,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493};


uint32_t FElementOptions::GetDefaultColor() const {
   if (m_iElement < 0 || size_t(m_iElement) > sizeof(ElementColors)/sizeof(ElementColors[0])) {
      assert(0); // no entry for this number...
      return 0xff00ff;
   }
   return ElementColors[size_t(m_iElement)];
}

uint32_t FElementOptions::GetDefaultBondColor() const
{
   return ColorNotSet;
//    return GetDefaultColor();
}

uint32_t FElementOptions::GetBondColor1() const
{
   if (m_BondColor != ColorNotSet)
      return m_BondColor;
   return m_Color;
}

double FElementOptions::GetDefaultCovalentRadius() const
{
   return ct::GetCovalentRadius(m_iElement);
}

double FElementOptions::GetDefaultVdwRadius() const
{
   return ct::GetVdwRadius(m_iElement, 0); // 0: flag that we do not want to crash for non-existent radii (will return -1 for these).
}


// float UffBondRadii[103] = {0.87,1.60,2.52,2.03,1.58,1.43,1.32,1.29,1.26,1.74,2.91,2.69,2.35,2.11,2.08,2.04,1.97,1.95,3.69,3.33,2.86,2.67,2.65,2.54,2.61,2.52,2.35,2.20,2.46,2.25,2.38,2.26,2.29,2.25,2.25,2.17,4.27,3.88,3.21,2.96,2.78,2.80,2.50,2.79,2.52,2.53,2.62,2.65,2.76,2.64,2.66,2.62,2.61,2.39,4.86,4.30,3.67,3.48,3.44,3.43,3.40,3.36,3.35,3.28,3.27,3.23,3.20,3.16,3.14,3.09,3.16,3.04,2.86,2.88,2.59,2.59,2.59,2.58,2.38,2.53,2.87,2.76,2.86,2.83,2.92,2.68,5.44,4.75,3.75,3.25,3.23,3.18,3.15,3.13,3.14,3.40,3.33,3.31,3.26,3.24,3.19,3.17,3.21};
// // static float AtomicRadii[103] = {0.76,0.64,2.68,1.80,1.64,1.54,1.50,1.46,1.42,1.38,3.08,2.60,2.36,2.22,2.12,2.04,1.98,1.94,3.92,3.48,2.88,2.72,2.50,2.54,2.78,2.50,2.52,2.42,2.76,2.62,2.52,2.44,2.38,2.32,2.28,2.20,4.22,3.84,3.24,2.96,2.74,2.90,3.12,2.52,2.70,2.62,3.06,2.96,2.88,2.82,2.76,2.70,2.66,2.60,4.50,3.96,3.38,3.48,3.44,3.43,3.40,3.36,3.35,3.28,3.27,3.23,3.20,3.16,3.14,3.09,3.20,3.00,2.76,2.92,3.18,2.56,2.74,2.56,2.88,2.98,2.96,2.94,2.92,2.83,2.92,2.90,5.44,4.75,3.75,3.25,3.23,3.18,3.15,3.13,3.14,3.40,3.33,3.31,3.26,3.24,3.19,3.17,3.21}; // see make_atomic_radii.py
static float AtomicRadii[104] = {0, 0.87,1.60,2.52,2.03,1.58,1.43,1.32,1.29,1.26,1.74,2.91,2.69,2.35,2.11,2.08,2.04,1.97,1.95,3.69,3.33,2.86,2.67,2.65,2.54,2.61,2.52,2.35,2.20,2.46,2.25,2.38,2.26,2.29,2.25,2.25,2.17,4.27,3.88,3.21,2.96,2.78,2.80,2.50,2.79,2.52,2.53,2.62,2.65,2.76,2.64,2.66,2.62,2.61,2.39,4.86,4.30,3.67,3.48,3.44,3.43,3.40,3.36,3.35,3.28,3.27,3.23,3.20,3.16,3.14,3.09,3.16,3.04,2.86,2.88,2.59,2.59,2.59,2.58,2.38,2.53,2.87,2.76,2.86,2.83,2.92,2.68,5.44,4.75,3.75,3.25,3.23,3.18,3.15,3.13,3.14,3.40,3.33,3.31,3.26,3.24,3.19,3.17,3.21};

inline float GetAtomDrawRadius(int iElement) { return AtomicRadii[iElement]; } // hmmm... looks rather non-similar to the covalent radius.
// float GetAtomDrawRadius(int iElement) { return GetCovalentRadius[iElement]; }

double FElementOptions::GetDefaultDrawRadius() const
{
   return GetAtomDrawRadius(m_iElement);
}

char const *FElementOptions::ElementName() const
{
   return ct::ElementNameFromNumber(m_iElement);
}

char const *BoolToCstr(bool o) { return o? "true" : "false"; }


FFreeObject::FFreeObject(QObject *parent)
   : QObject(parent)
{
}

QString FFreeObject::toString() const
{
   return "FFreeObject()";
}


// FFreeLine::FFreeLine(QObject *parent)
//    : FFreeObject(parent)
// {
//    InitProperties();
// }

FFreeLine::FFreeLine(FVec3d vFrom, FVec3d vTo, uint32_t dwColor, double fWidth, double fWeight, QObject *parent)
   : FFreeObject(parent), d(new FFreeLineData)
{
   InitProperties();
   SetFrom(vFrom);
   SetTo(vTo);
   SetWidth(fWidth);
   SetWeight(fWeight);
   SetColor(dwColor);
}


QTextStream &operator << (QTextStream &str, FVec3d const &v) {
   bool AsArray = false;
   str << (AsArray ? "[" : "Vec3(");
   for (size_t i = 0; i < v.size(); ++ i) {
      if (i != 0) str << ",";
      str << v[i];
   }
   str << (AsArray ? "]" : ")");
   return str;
}


QString FFreeLine::toString() const
{
//    QString res;
//    QTextStream str(&res);
//    str << "FreeLine(" << GetFrom() << ", " << GetTo() << ", " << GetColor() << ", " << GetWidth() << ")";
//    return res;
   return this->GetOptionsDesc(OPTIONSDESC_ObjectCtor);
}



FFreeLabel::FFreeLabel(QObject *parent)
   : FFreeObject(parent)
{
   InitProperties();
}


FFreeLabel::FFreeLabel(QString Text, FVec3d vPos, double fSize, uint32_t dwColor, QObject *parent)
   : FFreeObject(parent), d(new FFreeLabelData)
{
   InitProperties();
   SetText(Text);
   SetPos(vPos);
   SetSize(fSize);
   SetColor(dwColor);
}


QString FFreeLabel::toString() const
{
   return this->GetOptionsDesc(OPTIONSDESC_ObjectCtor);
//    QString res;
//    QTextStream str(&res);
//    str << "FreeLabel(\"" << GetText() << "\", " << GetPos() << ", " << GetSize() << ", " << GetColor() << ")";
//    return res;
}

FInvalidFreeObject::FInvalidFreeObject(QObject *parent)
   : FFreeObject(parent)
{}

FFreeObject::FNewFreeObjectPtr FInvalidFreeObject::newLinkedObject(QObject *parent) {
   return FNewFreeObjectPtr(new FInvalidFreeObject(parent));
}
FFreeObject::FNewFreeObjectPtr FInvalidFreeObject::newClonedObject(QObject *parent) const {
   return FNewFreeObjectPtr(new FInvalidFreeObject(parent));
}


void _RegisterMetaTypes()
{
   // these ones are needed to make this work with the property interface.
   // see: https://doc.qt.io/qt-5/qmetatype.html#qRegisterMetaType-1
   // without this, e.g., uint32_t types cannot be initialized from doubles (or anything else, for that matter).
   // I guess I might need to register conversion routines, too...
   // UPDATE: after adding
   //
   //    qRegisterMetaType<uint32_t>("uint32_t")
   //
   // the uint32_t colors work. At least on QT 4.8... didn't try QT 5 yet.
   // Also not sure about the pointer types.
   //
   // In QT 5.2+, one can register implicit conversion routines:
   //
   //    https://doc.qt.io/qt-5/qmetatype.html#registerConverter
   //
   // However, this does not work in QT 4.8
   qRegisterMetaType<uint32_t>("uint32_t");
   qRegisterMetaType<FVec3d>("FVec3d");
//    qRegisterMetaType<FFreeObject*>("FFreeObject*");
//    qRegisterMetaType<FVec3d*>("FVec3d*");
}


#include "prop_FWfOptions.cpp.inl"
#include "prop_FElementOptions.cpp.inl"
#include "prop_FFreeLine.cpp.inl"
#include "prop_FFreeLabel.cpp.inl"
