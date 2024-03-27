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

#include <limits>
#include "CxPrintLevel.h"
#include "CxParse1.h"
#include "format.h"

namespace ct {

unsigned const FPrintOptionsBase::InvalidOption = unsigned(-1); // <-- what GetOptionId functions should return for unrecognized options.


FPrintLevel::FPrintLevel(std::string const &sValue)
   : m_iLevel(0)
{
   this->SetArgs(ct::string_slice(sValue));
}

FPrintLevel::FPrintLevel(ct::string_slice slValue)
   : m_iLevel(0)
{
   this->SetArgs(slValue);
}


void FPrintLevel::SetArgs(string_slice slValue)
{
   ptrdiff_t
      iValue;
   bool
      bValue;
   if (slValue.try_convert_to_int(&iValue)) {
      m_iLevel = FBaseType(iValue);
      if (m_iLevel != iValue)
         throw FInputError(fmt::format("FPrintLevel::SetArgs(): argument '{}' out of range", slValue.to_str()));
   } else if (slValue.try_convert_to_bool(&bValue)) {
      m_iLevel = bValue ? FPrintLevel::Basic : FPrintLevel::Off;
   } else if (slValue == "basic") {
      m_iLevel = FPrintLevel::Basic;
   } else if (slValue == "more" || slValue == "extra") {
      m_iLevel = FPrintLevel::More;
   } else if (slValue == "most") {
      m_iLevel = FPrintLevel::Most;
   } else if (slValue == "all" || slValue == "full") {
      m_iLevel = FPrintLevel::All;
   } else if (slValue == "suppress" || slValue == "quiet") {
      m_iLevel = FPrintLevel::Suppress;
   } else {
      throw FInputError(fmt::format("FPrintLevel::SetArgs(): argument '{}' not recognized. Should be integer or one of: 'suppress', 'off', 'basic', 'more', 'most', 'all'.", slValue.to_str()));
   }
}


void FPrintOptionsBase::_ParsePrintOptions(FStoredType *pOptionValues, size_t NumOptions, string_slice const &Args, FGetOptionIdFn pOptionId)
{
   std::string
      sOn = "1";
   try {
      FPropertyListStr
         Props = FPropertyListStr(Args, 0, PROPLIST_AllowEntriesWithoutName);
      for (size_t iProp = 0; iProp != Props.size(); ++ iProp) {
         FPropertyStr
            Prop = Props[iProp]; // note: this makes a copy, but it's just two string_slice objects (4 ptrs total).
         if (Prop.Name.empty()) {
            // by default assume that if someone puts a print key without
            // arguments, that they want the corresponding thing printed at
            // level on/1.
            Prop.Name = Prop.Content;
            Prop.Content = string_slice(sOn);
         }
         // process the option name
         unsigned iOption = pOptionId(Prop.Name);
         if (iOption == InvalidOption)
            throw FInputError(fmt::format("TPrintOptions: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
         if (iOption >= NumOptions)
            throw FInputError(fmt::format("TPrintOptions: option '{}: {}' recognized, but the option's flag id was mapped to {}, which exceeds the number of supported flags ({}). This is a programming error.", Prop.Name.to_str(), Prop.Content.to_str(), iOption, NumOptions));

         // process the value and store
         _Set1(iOption, Prop.Content, pOptionValues, NumOptions);
      }
   } catch (FInputError const &e) {
      throw;
   } catch (std::runtime_error const &e) {
      // convert other exceptions to input errors.
      throw FInputError(e.what());
   }
}


void FPrintOptionsBase::_Set1(unsigned iOption, string_slice const &slValue, FStoredType *pOptionValues, size_t NumOptions)
{
   FPrintLevel e;
   e.SetArgs(slValue);
   _Set1(iOption, e, pOptionValues, NumOptions);
}

void FPrintOptionsBase::_Set1(unsigned iOption, FPrintLevel const &pl, FStoredType *pOptionValues, size_t NumOptions)
{
   FPrintLevel::FBaseType
      ei(pl); // <- convert to integer
   if (iOption >= NumOptions)
      throw FInputError(fmt::format("TPrintOptions::Set(): option id '{}' does not exist (attempted to set to '{}'; have only option ids 0 ... {})", iOption, ei, NumOptions));
   FStoredType
      iValue(ei);
   if (iValue != ei)
      throw FInputError(fmt::format("TPrintOptions::Set(): option {} numeric value ({}) exceeds supported range [{}, {}]", iOption, iValue, std::numeric_limits<FStoredType>::min(), std::numeric_limits<FStoredType>::max()));

   // well, we arrived here. All good, I guess?
   pOptionValues[iOption] = iValue;
}


} // namespace ct
