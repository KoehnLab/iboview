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

/* CxAtomData.h v20200519 EST [mirage, Gerald Knizia] */
#ifndef ATOMIC_DATA_H
#define ATOMIC_DATA_H

#include <stdexcept>

namespace ct {
   class EDataRequestError : public std::runtime_error {
   public:
      EDataRequestError(std::string const &s);
   private:
      typedef std::runtime_error
         FBase;
   };

   struct FIsotopeInfo {
      int
         iElement, nNucleons;
      double
         fMass, fAbundance;
      char const
         *pName;
      // returns true if this elements does not occur in nature
      // (and thus has no data on natural abundance)
      inline bool Synthetic() const { return fAbundance == -1.; };
   };

   enum FAtomicMassType {
      ATMASS_MostCommonIsotope,
      ATMASS_StandardAtomicWeight,
      ATMASS_ChargeInsteadOfMass
   };

   enum FDataRequestOptions {
      DATAREQUEST_AssertDataExists = 0x01 // if set, raise EDataRequestError if requested data is not tabulated.
   };

   void CheckElement(int iElement);
   double GetAtomicMass(int iElement, FAtomicMassType Type);
   char const *ElementNameFromNumber(int iElement);
   int ElementNumberFromName(char const *pName);
   double GetCovalentRadius(int iElement);
    // WARNING: not set for all elements!
   double GetVdwRadius(int iElement, unsigned Flags = DATAREQUEST_AssertDataExists);
   double GetVdwRadius_IsoDensity(int iElement);

   // note: [0] is dummy. I.e., g_ElementNames[1] is "H".
   int const g_nMaxElement = 118; // highest recognized iElement
   extern char const *g_ElementNames[119];
   extern double const g_ElementAverageMasses[119];
   extern double const g_ElementMostCommonIsotopeMasses[119];

   extern FIsotopeInfo const g_IsotopeInfo[355];
   extern unsigned const g_IsotopeInfoElementOffsets[120];

   extern double const g_CovalentRadii[119];
   extern double const g_VdwRadii[119];

} // namespace ct

#endif // ATOMIC_DATA_H
