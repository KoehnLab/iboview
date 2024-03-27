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

#ifndef CX_PRINT_LEVEL_H
#define CX_PRINT_LEVEL_H

#include <string>
#include <array>
// #include <functional>

namespace ct {
   struct string_slice; // CxParse1.h/.cpp
}

namespace ct {

// A wrapped integer we use to control print levels of different
// parts of the program
struct FPrintLevel {
   typedef int FBaseType;
   // this is a pointer-to-member, which can be implicitly converted
   // to bool, but not many other things.
   typedef FBaseType FPrintLevel::*unspecified_bool_type;

//    static constexpr FBaseType
   enum FDefaultLevel {
      Suppress = -1, Off = 0, Basic = 1, More = 2, Most = 3, All = 4
   };

   // convert from / to integer level
   FPrintLevel(FBaseType i = Off) : m_iLevel(i) {} // <-- note: *not* explicit!
   explicit operator FBaseType () const { return m_iLevel; }
   // construct/init from string or string_slice
   explicit FPrintLevel(std::string const &sValue);
   explicit FPrintLevel(string_slice slValue);
   void SetArgs(string_slice slValue);

   // assignment/increment/decrement with integer levels
   FPrintLevel &operator = (FBaseType i) { m_iLevel = i; return *this; }
   void operator += (FBaseType i) { m_iLevel += i; }
   void operator -= (FBaseType i) { m_iLevel -= i; }
   // distinct object for higher or lower detail, for easy control of sub-processes.
   FPrintLevel operator + (FBaseType i) { return FPrintLevel(m_iLevel + i); }
   FPrintLevel operator - (FBaseType i) { return FPrintLevel(m_iLevel - i); }

   // on/off level tests (`if (PrintLevel) { ... }` and `if (!PrintLevel) { ... }`)
   operator unspecified_bool_type () const { return (m_iLevel >= 1)? &FPrintLevel::m_iLevel : 0; }
   bool operator !() const { return m_iLevel <= 0; }

   // queries for comparison to standard levels
   bool SupressQ() const { return m_iLevel <= Suppress; }
   bool OffQ() const { return m_iLevel <= Off; }
   bool BasicQ() const { return m_iLevel >= Basic; }
   bool MoreQ() const { return m_iLevel >= More; }
   bool MostQ() const { return m_iLevel >= Most; }
   bool AllQ() const { return m_iLevel >= All; }
   bool BasicOnlyQ() const { return m_iLevel == Basic; }
   bool MoreOnlyQ() const { return m_iLevel == More; }
   bool MostOnlyQ() const { return m_iLevel == Most; }

   // comparison to integer levels
   bool operator == (FBaseType i) const { return m_iLevel == i; }
   bool operator != (FBaseType i) const { return m_iLevel != i; }
   bool operator <= (FBaseType i) const { return m_iLevel <= i; }
   bool operator >= (FBaseType i) const { return m_iLevel >= i; }
   bool operator <  (FBaseType i) const { return m_iLevel < i; }
   bool operator >  (FBaseType i) const { return m_iLevel > i; }

   bool AtLeastQ(FBaseType i) const { return m_iLevel >= i; }
   bool ExactlyQ(FBaseType i) const { return m_iLevel == i; }
protected:
   FBaseType m_iLevel;
};


struct FPrintOptionsBase
{
   enum { DefaultLevel = FPrintLevel::Off };
   typedef signed char FStoredType;

   // a function which should map the string argument to a flag id in range 0
   // ... NumOptionsT-1. E.g. if the program has options like
   // "print:{orbitals:1, timing:2}", this function should map strings
   // "orbitals" or "timing" to their id. If unrecognized, should return
   // InvalidOption
   typedef unsigned (*FGetOptionIdFn)(string_slice const &sl);
   static unsigned const InvalidOption;
protected:
   static void _ParsePrintOptions(FStoredType *pOptionValues, size_t NumOptions, string_slice const &Args, FGetOptionIdFn pOptionId);
   static void _Set1(unsigned iOption, string_slice const &slValue, FStoredType *pOptionValues, size_t NumOptions);
   static void _Set1(unsigned iOption, FPrintLevel const &pl, FStoredType *pOptionValues, size_t NumOptions);
};


template<size_t NumOptionsT>
struct TPrintOptions : public FPrintOptionsBase {
   enum { nOptions = NumOptionsT };

   explicit TPrintOptions(FPrintLevel Default = DefaultLevel) { this->SetAll(Default); }
   explicit TPrintOptions(std::string const &Args, FGetOptionIdFn pOptionId);

   FPrintLevel operator [](unsigned iOption) const { return FPrintLevel(m_OptionValues.at(iOption)); }

   void SetArgs(std::string const &Args, FGetOptionIdFn pOptionId);
   void SetArgs(string_slice const &Args, FGetOptionIdFn pOptionId);
   void Set(unsigned iOption, FPrintLevel const &pl) { _Set1(iOption, pl, &m_OptionValues[0], m_OptionValues.size()); };
   void Set(unsigned iOption, string_slice const &slValue) { _Set1(iOption, slValue, &m_OptionValues[0], m_OptionValues.size()); };
   void SetAll(FPrintLevel pl) { m_OptionValues.fill(FStoredType(int(pl))); }

protected:
   typedef std::array<FStoredType, NumOptionsT>
      FSettingsArray;
   FSettingsArray
      m_OptionValues;
};


template<size_t NumOptionsT>
TPrintOptions<NumOptionsT>::TPrintOptions(std::string const &Args, FGetOptionIdFn pOptionId) {
   m_OptionValues.fill(DefaultLevel);
   SetArgs(Args, pOptionId);
}


template<size_t NumOptionsT>
void TPrintOptions<NumOptionsT>::SetArgs(std::string const &Args, FGetOptionIdFn pOptionId) {
   _ParsePrintOptions(&m_OptionValues[0], NumOptionsT, string_slice(Args), pOptionId);
}


template<size_t NumOptionsT>
void TPrintOptions<NumOptionsT>::SetArgs(string_slice const &Args, FGetOptionIdFn pOptionId) {
   _ParsePrintOptions(&m_OptionValues[0], NumOptionsT, Args, pOptionId);
}


} // namespace ct

#endif // CX_PRINT_LEVEL_H
