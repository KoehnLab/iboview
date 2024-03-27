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

#include <cctype> // for isalpha, toupper, to lower
#include <cstdlib> // for strtod
#include "CxRawAtom.h"
#include "CxAtomData.h" // for ElementNameFromNumber
#include "format.h"

namespace ct {


void FIoElementAndTag::Parse(char const *pFirst, char const *pLast)
{
   this->iElement = 0;
   this->Tag = ATOMTAG_None;

   assert(pFirst <= pLast);
   std::string
      sElement;
   sElement.reserve(pLast - pFirst);
   // copy over characters from [pFirst,pLast) from left until we find a
   // non-alphabetic one.
   char const *it = pFirst;
   for (; it != pLast; ++ it) {
      int
         c = int(*it);
      if (!std::isalpha(c))
         // first non-alphabetic character.
         break;
      if (it == pFirst)
         // upper-case-ify first letter of element symbol
         c = std::toupper(c);
      else
         // lower-case-ify all other letters of element symbol
         c = std::tolower(c);
      sElement.push_back(char(c));
   }

   // we should now be ready to interpret the element symbol.
   this->iElement = ct::ElementNumberFromName(sElement.c_str());

   // `it` should now point to the end of the initial alpha-character sequence
   // making up the element symbol. Check for a following numeric atom tag next.
   // E.g., "C1", "Na30", etc. We actually allow doubles for those.

   if (it == pLast) {
      // that's it---all characters of input string are processed, and there is
      // nothing left to possible represent a tag. We're done.
      return;
   } else {
      // some characters are left. Check if we can interpret those as a double
      // float. If yes, take this as atom tag.
      char const
         *pTagBeg = it,
         *pTagEnd = 0;
      this->Tag = std::strtod(pTagBeg, const_cast<char**>(&pTagEnd));
      if (pTagEnd != pLast)
         throw std::runtime_error(fmt::format("FIoElementAndTag::Parse(): failed to convert atom subtype tag '{}' to a number while parsing atom type '{}'.", std::string(pTagBeg,pLast), std::string(pFirst, pLast)));
   }
}


std::string CombineElementAndAtomTag(int iElement, FAtomTag const &Tag)
{
   if (Tag != ATOMTAG_None) {
      return fmt::format("{}{}", ElementNameFromNumber(iElement), Tag);
   } else {
      return ElementNameFromNumber(iElement);
   }
}


std::string FIoElementAndTag::Format() const
{
   return CombineElementAndAtomTag(this->iElement, this->Tag);
}


char const *FIoElementAndTag::Element() const
{
   return ElementNameFromNumber(iElement);
}


} // namespace ct
