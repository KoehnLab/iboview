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

#ifndef IV_IRC_H
#define IV_IRC_H

// #include "IvDocument.h"
class FDocument;
namespace ct {
   class FAtomSet;
}

bool MakeIboChangeCurve(TArray<float> &CurveData, uint iRow, FDocument *document);


// maybe other options? e.g., on curve approximations, use of gradients, etc?
struct FArcLengthOptions
{
   enum FWeightType {
      // if set, do not do the mass weighting. This will not produce approximate IRC coordinates, but
      // actual arc lengths as X-axis.
      ATOMWEIGHT_Unity,
      // if set, make approximate IRC arc length based on average atomic masses.
      ATOMWEIGHT_IrcAvgMass,
      // if set, make approximate IRC arc length based on most-common-isotope atomic masses.
      ATOMWEIGHT_IrcIsotopeMass,
      // if set, weight atoms by nuclear charge.
      ATOMWEIGHT_Charge,
      // if set, this is not actually an arc length, but just a frame index (used when
      // no consistent arc length could be assigned)
      ATOMWEIGHT_FrameIndex
   };
   FWeightType
      WeightType;

   explicit FArcLengthOptions(FWeightType WeightType_);
   double GetAtomWeight(int iElement, double fMassOverride = -1., double fNuclearChargeOverride = -1.) const;
   // pAtoms: if given: take from there; otherwise: take from current frame.
   double GetAtomWeight(int iAt, FDocument *document, ct::FAtomSet *pAtoms = 0) const;

   QString GetArcUnit() const;
   QString GetArcCaption() const;
};

// enum FIrcArcLengthFlags {
//    ARCLENGTH_NoMassWeighting = 0x1,
//    // if set, use most-common-isotope masses instead of average masses.
//    ARCLENGTH_UseIsotopeMasses = 0x2
// };

struct FArcLength
{
   FArcLength();
   FArcLength(double fValue_, FArcLengthOptions Type_);

   double Value() const { return m_Value; }
   FArcLengthOptions ArcType() const { return m_ArcType; }
protected:
   double
      m_Value;
   FArcLengthOptions
      m_ArcType;
};


void MakeIrcArcLengths(TArray<FArcLength> &ArcLengths, FDocument *document, FArcLengthOptions ArcOptions);
TArray<float> ExtractFloats(TArray<FArcLength> const &ArcLengths);


#endif // IV_IRC_H
