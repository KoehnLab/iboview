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

// This file provides support for parsing and representing parameters/settings
// which apply to atoms (e.g., basis set specifications, but really any type of
// settings).
//
// The core idea here is that we provide a unified interface for handling both
// defaults of the corresponding parameters, but also for setting diverging
// parameters for:
//
// - specific elements: e.g.,
//
//     basis='cc-pVTZ,O:aug-cc-pVTZ'
//
//   denotes that a default basis of 'cc-pVTZ' is set, but all oxygen atoms get
//   assigned an aug-cc-pVTZ instead. In
//
//     basis='cc-pVTZ,[C,O,N]:aug-cc-pVTZ'
//
//   the alternative basis setting applies to all carbon, oxygen, and nitrogen
//   atoms.
//
// - atoms marked with a specific tag in the input file: e.g.,
//
//     basis='cc-pVDZ,.2:cc-pVTZ,O2:aug-cc-pVTZ'
//
//   to denote that the default basis is cc-pVDZ, but all atoms with a tag of
//   '2' (".2") get cc-pVTZ instead, except Os with tag '2' ("O2"), which get
//   aug-cc-pVTZ.
//
// - atoms with certain indices (1-based) in the total frame: e.g.,
//
//     basis='cc-pVDZ,[#1,#2,#3,#4]:cc-pVTZ'
//
//   to denote that the default basis is cc-pVDZ, but the first four atoms in
//   the frame (#1, ..., #4) get a cc-pVTZ instead.
//
//   These can also be specified in numerical ranges, e.g.,
//
//     basis='cc-pVDZ,#1--4:cc-pVTZ'
//
// Precedence of settings is as follows (highest to lowest):
//
// - settings for explicit atom indices
// - settings for combinations of elements and atom tags
// - settings for atom tags
// - settings for elements
// - global default
//
// Provisions for specifying settings for ranges of atom indices, elements, or
// tags, may be added later.
#ifndef CX_ATOM_PARAM_SPEC_H
#define CX_ATOM_PARAM_SPEC_H
// #define CX_PARSE_ATOM_SETTINGS_H

#include <cctype> // for isalpha + co
#include <map>
#include <tuple>
#include <stdexcept>
#include <string>

#include "format.h"
#include "CxTypes.h"
#include "CxParse1.h"
#include "CxAtomData.h" // for ElementNumberFromName

namespace ct {


/// Manages per-atom settings, which can be provided via default values, per-
/// element overrides, per-atom-tag override, per-combined-(element, atom tag)
/// overrides, or per atom-index overrides. This class is mainly meant for user
/// interface purposes (for allowing end users to communicate their intended
/// settings to use for all or specific atoms to the host program, e.g., basis
/// set choices).
template<class FSettingProxy_>
class TAtomParamSpec : public FIntrusivePtrDest
{
public:
   typedef FSettingProxy_
      FSettingProxy;
   typedef typename FSettingProxy::FStoredType
      FStoredType;
   typedef typename FSettingProxy::FSettingType
      FSettingType;

   typedef int
      FElement; // number of chemical element (e.g., 1->H, 2->He, ...)
   typedef double
      FAtomTag; // explicit subtype tag of an atom (if given in input)
   typedef size_t
      FAtomIndex; // index of an atom inside the geometry frame (first atom of molecule is 0, second atom 1, ..., last atom is nAt-1)

   TAtomParamSpec(string_slice slFullSpec_ = string_slice(), FStoredType const &Default_ = FStoredType(), char const *pOpenDelimCloseDef_=0, FSettingProxy const &Proxy_ = FSettingProxy());
   virtual ~TAtomParamSpec() {}; // not very pretty, but probably can't be helped.

   /// Add specifications from slFullSpec_ to current set of settings (possibly
   /// overwriting previous ones).
   ///
   /// If provided, pOpenDelimCloseDef_ must be a c-string of length 4 defining
   /// the input string format; e.g., "{,}:" means:
   /// - '{' is the (optional) opening parenthesis around the entire specification
   /// - ',' is used to delimit different sub-specifications
   /// - ':' is used to separate key-value pair definitions
   /// - '}' is the (optional) closing parenthesis around the entire specification
   ///  With this, e.g.
   ///
   ///    "{cc-pVDZ,.2:cc-pVTZ,O2:aug-cc-pVTZ}"
   ///
   /// is a valid specification (for string types).
   /// If instead pOpenDelimCloseDef_ = "{;}=" were used, the same specification
   /// would read as
   ///
   ///    "{cc-pVDZ;.2=cc-pVTZ;O2=aug-cc-pVTZ}"
   ///
   /// The default is '{;}:'.
   void Update(string_slice slFullSpec_, char const *pOpenDelimCloseDef_=0);

//    void Update1(string_slice slKey_, string_slice slValue_);
   void Update1(string_slice slKey_, FStoredType const &Value_);

   /// Find and return the value of the target setting for an atom with index iAt
   /// (0-based!), element iElement, and atom tag iTag; this is done by tracing
   /// the setting overrides from the highest priority to the lowest (eventually
   /// returning the value given as default if another one was not specified).
   FSettingType const &get(FAtomIndex iAt, FElement iElement, FAtomTag iTag) const {
      return m_Proxy.SettingValue(this->FindStoredValue(iAt, iElement, iTag));
   }
   /// same as this->get(iAt, iElement, iTag).
   FSettingType const &operator()(FAtomIndex iAt, FElement iElement, FAtomTag iTag) const {
      return this->get(iAt, iElement, iTag);
   }

   FSettingType const &GetDefault() const { return m_Default; }
   void SetDefault(FSettingType const &NewDefault) { m_Default = NewDefault; }
protected:
   char const
      *m_pDefaultOpenDelimCloseDef;
   FSettingProxy
      m_Proxy;
   FStoredType
      m_Default;

   typedef std::map<FElement, FStoredType>
      FElementSpecs;
   FElementSpecs
      m_ElementSpecs;

   typedef std::map<FAtomTag, FStoredType>
      FAtomTagSpecs;
   FAtomTagSpecs
      m_AtomTagSpecs;

   typedef std::pair<FElement, FAtomTag>
      FElementAndTag;
   typedef std::map<FElementAndTag, FStoredType>
      FElementAndTagSpecs;
   FElementAndTagSpecs
      m_ElementAndTagSpecs;

   typedef std::map<FAtomIndex, FStoredType>
      FAtomIndexSpecs;
   FAtomIndexSpecs
      m_AtomIndexSpecs;

   FStoredType const &FindStoredValue(FAtomIndex iAt, FElement iElement, FAtomTag iTag) const;
};


template<class FSettingProxy>
TAtomParamSpec<FSettingProxy>::TAtomParamSpec(string_slice slFullSpec_, typename FSettingProxy::FStoredType const &Default_, char const *pOpenDelimCloseDef_, FSettingProxy const &Proxy_)
   : m_Proxy(Proxy_), m_Default(Default_)
{
   if (pOpenDelimCloseDef_)
      m_pDefaultOpenDelimCloseDef = pOpenDelimCloseDef_;
   else
      m_pDefaultOpenDelimCloseDef = "{;}:";

   Update(slFullSpec_);
}


template<class FSettingProxy>
typename FSettingProxy::FStoredType const &TAtomParamSpec<FSettingProxy>::FindStoredValue(FAtomIndex iAt, FElement iElement, FAtomTag iTag) const
{
   // check overrides by priority, from highest to lowest.
   typename FAtomIndexSpecs::const_iterator
      itAtomIndexSpec = m_AtomIndexSpecs.find(iAt);
   if (itAtomIndexSpec != m_AtomIndexSpecs.end())
      return itAtomIndexSpec->second;

   typename FElementAndTagSpecs::const_iterator
      itElementAndTagSpec = m_ElementAndTagSpecs.find(FElementAndTag(iElement, iTag));
   if (itElementAndTagSpec != m_ElementAndTagSpecs.end())
      return itElementAndTagSpec->second;

   typename FAtomTagSpecs::const_iterator
      itAtomTagSpec = m_AtomTagSpecs.find(iTag);
   if (itAtomTagSpec != m_AtomTagSpecs.end())
      return itAtomTagSpec->second;

   typename FElementSpecs::const_iterator
      itElementSpec = m_ElementSpecs.find(iElement);
   if (itElementSpec != m_ElementSpecs.end())
      return itElementSpec->second;

   return m_Default;
}


// a predicate for string_slice::find_c
struct is_nonalpha_c {
   inline bool operator() (char c) const {
      return !std::isalpha(static_cast<unsigned char>(c));
   }
};


template<class FSettingProxy>
// void TAtomParamSpec<FSettingProxy>::Update1(string_slice slKey_, string_slice slValue_)
void TAtomParamSpec<FSettingProxy>::Update1(string_slice slKey_, FStoredType const &Value_)
{
   slKey_.trim();
   if (slKey_.empty()) {
      // no key given---set new default.
      m_Default = Value_;
   } else if (slKey_.startswith('[') && slKey_.endswith(']')) {
      // given a list of keys to assign the same value to.
      // take apart the key list and set the given entries by calling Update1
      // again for each key entry.
      long_split_result
         sr;
      slKey_.split_list(sr, '[', ']', ',');
      for (long_split_result::const_iterator itSubKey = sr.begin(); itSubKey != sr.end(); ++ itSubKey) {
         Update1(*itSubKey, Value_);
      }
   } else if (slKey_.startswith('.')) {
      // specifies an atom tag (.123).
      string_slice sAtomTag(slKey_);
      sAtomTag.first += 1;
      FAtomTag iAtomTag;
      sAtomTag.convert(&iAtomTag);
      m_AtomTagSpecs[iAtomTag] = Value_;
   } else if (slKey_.startswith('#')) {
      // specifies an 1-based atom index in the molecule (#1, #2, ...)
      string_slice sAtomIndex(slKey_);
      sAtomIndex.first += 1;
      // check if this is meant to specify a numerical range of atom indices
      cstr_it itRs = sAtomIndex.find("--");
      if (itRs == sAtomIndex.last) {
         // no. just a single index. try to convert it.
         FAtomIndex iAtomIndex;
         sAtomIndex.convert(&iAtomIndex);
         iAtomIndex -= 1; // convert from 1-based iAt (input) to 0-based iAt (internal)
         m_AtomIndexSpecs[iAtomIndex] = Value_;
      } else {
         // yes, index range. restrict first index slice to part
         // before the range separator, and make second index slice
         // for the rest of the string
         sAtomIndex.last = itRs;
         string_slice
            sAtomIndexEnd(itRs+2, slKey_.last);
         // convert to 1-based (external) numerical indices
         FAtomIndex
            iAtBeg, iAtEnd;
         sAtomIndex.convert(&iAtBeg);
         sAtomIndexEnd.convert(&iAtEnd); // note: inclusive for once!
         // assign values to all individual atoms (see comment below on
         // storing non-overlapping ranges if interested...)
         for (FAtomIndex iAt = iAtBeg; iAt <= iAtEnd; ++ iAt)
            m_AtomIndexSpecs[iAt-1] = Value_;
      }
   } else {
      // should be either just element or element + atom tag. Find out what it is.
      cstr_it
         iFirstNonAlpha = slKey_.find_c(is_nonalpha_c());
      // assume first alphabetic (NOT alphanumeric!) characters make up the
      // element symbol.
      string_slice
         slElement(slKey_.first, iFirstNonAlpha),
         slRest(iFirstNonAlpha, slKey_.last);
      FElement
         iElement;
      try {
         iElement = ElementNumberFromName(slElement.to_str().c_str());
      } catch(EDataRequestError const &e) {
         throw std::runtime_error(fmt::format("expected '{}' to define an element symbol or element+numerical tag pair, but failed to look up element index of '{}' ('{}')", slKey_.to_str(), slElement.to_str(), e.what()));
      }
      if (iFirstNonAlpha == slKey_.last) {
         // *all* characters are alphabetic. So we should have *only* an element
         // symbol in this case.
         m_ElementSpecs[iElement] = Value_;
      } else if (slRest.startswith("--")) {
         // range of elements, from iElement to next element after it.
         slRest.first += 2; // skip the "--"
         FElement
            iElementEnd; // note: inclusive!
         try {
            iElementEnd = ElementNumberFromName(slRest.to_str().c_str());
         } catch(EDataRequestError const &e) {
            throw std::runtime_error(fmt::format("expected '{}' to define a range Ab--Cd of element symbols, but failed to look up element index of '{}' ('{}')", slKey_.to_str(), slRest.to_str(), e.what()));
         }
         // set target value to all the elements in range individually
         // (note that if we really wanted to, we totally could map *ranges* of
         // numerical assignments instead of individual assignments (as long as
         // the ranges do not overlap; i.e., adding a new range should
         // ltrim/rtrim/split all previous ranges it overlapped with... maybe
         // later)
         for (FElement i = iElement; i <= iElementEnd; ++ i)
            m_ElementSpecs[i] = Value_;
      } else {
         // assume starting at the first non-alpha character, the rest
         // make up a numeric atom tag.
         string_slice
            sAtomTag = slRest;
         FAtomTag iAtomTag;
         try {
            sAtomTag.convert(&iAtomTag);
         } catch (std::runtime_error const &e) {
            throw std::runtime_error(fmt::format("expected '{}' to define an element symbol or element+numerical tag pair, but failed to interpret '{}' as a numerical atom tag ('{}')", slKey_.to_str(), sAtomTag.to_str(), e.what()));
         }
         m_ElementAndTagSpecs[FElementAndTag(iElement, iAtomTag)] = Value_;
      }
   }
}


template<class FSettingProxy>
void TAtomParamSpec<FSettingProxy>::Update(string_slice slFullSpec_, char const *pOpenDelimCloseDef_)
{
   slFullSpec_.trim();
   if (!slFullSpec_.empty()) {
      if (pOpenDelimCloseDef_ == 0)
         pOpenDelimCloseDef_ = m_pDefaultOpenDelimCloseDef;
      unsigned
         PropFlags = PROPLIST_AllowOneEntryWithoutName | PROPLIST_AllowOmitOpenClose | PROPLIST_AllowComplexKeys;

      FPropertyListStr
         Props(slFullSpec_, 0, PropFlags, pOpenDelimCloseDef_);
      for (FPropertyListStr::const_iterator itProp = Props.begin(); itProp != Props.end(); ++ itProp) {
         Update1(itProp->Name, m_Proxy.Parse(itProp->Content));
   }
//       if (itProp->Name.empty()) {
//          // no key given---set new default.
//          m_Default = m_Proxy.Parse(itProp->Content);
//       } else if (itProp->Name.startswith('[') && itProp->Name.endsswith(']')) {
//          // given a list of keys to assign the same value to.
//          // take apart the key list and set the given entries.
//          FStoredType
//             Value = m_Proxy.Parse(itProp->Content);
//
//       } else {
//          // assume single ket set---forward to Update1 to find out what key
//          // this might be and assign it
//          Update1(itProp->Name, m_Proxy.Parse(itProp->Content));
//       }
   }
}



// base class for settings which are stored by-value
// (note: we might want to add Format() functions to those (calling either
// Desc()/MakeDesc() or format() directly) in case we need functions for
// assembling summary descriptions of the atom parameter objects.)
template<class FType>
struct TSettingProxy_ByValue
{
   typedef FType FSettingType;
   typedef FType FStoredType;
   FSettingType const &SettingValue(FStoredType const &StoredValue) const { return StoredValue; }
   FSettingType &SettingValue(FStoredType &StoredValue) const { return StoredValue; }
};


// base class for settings which are stored indirectly via TIntrusivePtr<FType>
// (so FSettingType must be derived from FIntrusivePtrDest)
template<class FType>
struct TSettingProxy_ByPtr
{
   typedef TIntrusivePtr<FType>
      FTypePtr;
   typedef FType FSettingType;
   typedef FTypePtr FStoredType;
   FSettingType const &SettingValue(FStoredType const &StoredValue) const { return *StoredValue; }
   FSettingType &SettingValue(FStoredType &StoredValue) const { return *StoredValue; }
};


// proxy object which generates FSettingType objects by calling their
// FSettingType(string_slice s)-compatible constructor (only works on
// types which have such a constructor)
template<class FType>
struct TSettingParseProxy_CtorCallDirect : public TSettingProxy_ByValue<FType>
{
   inline FType Parse(string_slice sl) const { return FType(sl); }
};

// similar to TSettingParseProxy_CtorCallDirect, but the settings are not stored
// by-value but rather indirectly via an intrusive ptr.
// (so FSettingType must be derived from FIntrusivePtrDest, and must
// additionally have a constructor which can be called with a string_slice
// argument)
template<class FType>
struct TSettingParseProxy_CtorCallPtr : public TSettingProxy_ByPtr<FType>
{
   typedef TSettingProxy_ByPtr<FType> FBase;
   typedef typename FBase::FTypePtr FTypePtr;
   inline FTypePtr Parse(string_slice sl) const { return FTypePtr(new FType(sl)); }
};


// proxy objects which convert the input string slice to basic types (strings,
// signed/unsigned integers, floats, or bools), via invoking
// string_slice::convert directly.
template<class FType>
struct TSettingParseProxy_BaseType : public TSettingProxy_ByValue<FType> {
   typedef typename TSettingProxy_ByValue<FType>::FStoredType FStoredType;
   FStoredType Parse(string_slice sl) const { FType val; sl.convert(&val); return val; }
};

typedef TSettingParseProxy_BaseType<std::string> FSettingParseProxy_String;
typedef TSettingParseProxy_BaseType<ptrdiff_t> FSettingParseProxy_Int;
typedef TSettingParseProxy_BaseType<size_t> FSettingParseProxy_Size;
typedef TSettingParseProxy_BaseType<double> FSettingParseProxy_Float;
typedef TSettingParseProxy_BaseType<bool> FSettingParseProxy_Bool;


// represents atom parameters given by a raw character string; the object just
// splits the input string into individual sub-strings, but does not attempt to
// process/parse those further.
typedef TAtomParamSpec<FSettingParseProxy_String> FAtomParamSpec_String;
typedef TIntrusivePtr<FAtomParamSpec_String> FAtomParamSpecPtr_String;
typedef TIntrusivePtr<FAtomParamSpec_String const> FAtomParamSpecCptr_String;

// represents atom parameters in the form of signed integers
typedef TAtomParamSpec<FSettingParseProxy_Int> FAtomParamSpec_Int;
typedef TIntrusivePtr<FAtomParamSpec_Int> FAtomParamSpecPtr_Int;
typedef TIntrusivePtr<FAtomParamSpec_Int const> FAtomParamSpecCptr_Int;

// represents atom parameters in the form of unsigned integers
typedef TAtomParamSpec<FSettingParseProxy_Size> FAtomParamSpec_Size;
typedef TIntrusivePtr<FAtomParamSpec_Size> FAtomParamSpecPtr_Size;
typedef TIntrusivePtr<FAtomParamSpec_Size const> FAtomParamSpecCptr_Size;

// represents atom parameters in the form of floating point numbers
typedef TAtomParamSpec<FSettingParseProxy_Float> FAtomParamSpec_Float;
typedef TIntrusivePtr<FAtomParamSpec_Float> FAtomParamSpecPtr_Float;
typedef TIntrusivePtr<FAtomParamSpec_Float const> FAtomParamSpecCptr_Float;

// represents atom parameters in the form of on/off switches
typedef TAtomParamSpec<FSettingParseProxy_Bool> FAtomParamSpec_Bool;
typedef TIntrusivePtr<FAtomParamSpec_Bool> FAtomParamSpecPtr_Bool;
typedef TIntrusivePtr<FAtomParamSpec_Bool const> FAtomParamSpecCptr_Bool;


// typedef TAtomParamSpec<std::string, typename TSettingParseProxy_CtorCall<string_slice> >
//    FAtomParamSpec;





}

#endif // CX_ATOM_PARAM_SPEC_H
