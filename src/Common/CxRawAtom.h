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

#ifndef CX_RAW_ATOM_H
#define CX_RAW_ATOM_H
// #include "CxPodArray.h" // <- moved to std::vector for the time being. doesn't really matter for this.
#include <string>
#include <vector>
#include <limits>
#include <utility> // for std::pair
#include "CxTypes.h"
#include "CxVec3.h"

namespace ct {
typedef TVector3<double>
   FVec3d;

/// Type of optional value denoting the numeric sub-type of atoms of a certain
/// element, or ATOMTAG_None if not present.
typedef double
   FAtomTag;
/// used to mark the absence of explicit atom tags.
FAtomTag const ATOMTAG_None = std::numeric_limits<FAtomTag>::min();

/// A simple atom structure used *only* to identify where atoms are, which
/// chemical element they are of, and possibly which sub-type each atom has been
/// labelled with (e.g. to assign different parameter sets to different sets of
/// atoms of the same element).
///
/// This structure is used in place of FAtomSet/FXyzFrame objects
/// for interfacing purposes, so that mostly-autarc parts of the code, such as
/// the space partitioning stuff, can be isolated from the base machinery
///
/// There are MakeRawAtomList() function to convert FAtomSet or FXyzFrame
/// objects to these guys.
struct FRawAtom {
   FVec3d
      vPos;
   int
      iElement;
   FAtomTag
      iTag;
   FRawAtom() {};
   FRawAtom(FVec3d vPos_, int iElement_, FAtomTag iTag_) : vPos(vPos_), iElement(iElement_), iTag(iTag_) {}
};
/// A minimalistic specification of a molecular geometry for interfacing purposes,
/// as an array of atoms (the "list" is meant in everyday-speak).
typedef std::vector<FRawAtom>
   FRawAtomList;

// typedef TArray<FRawAtom>
//    FRawAtomList;


// interface functions for building FRawAtomList objects from more complex
// geometry representations (note: used to be called "RawAtomListFromAtomSet")

struct FXyzFrame;
/// Build a FRawAtomList object from a CxBaseLib FXyzFrame object
FRawAtomList AsRawAtoms(FXyzFrame const &Atoms); // in CxXyzFrame.cpp

struct FAtomSet;
/// Build a FRawAtomList object from a MicroScf FAtomSet object
FRawAtomList AsRawAtoms(FAtomSet const &Atoms); // in CtAtomSet.cpp


typedef std::pair<char const *, char const *>
   pchar_range;
// ^-- string_slice objects can be converted to those (via type case operator
//     pchar_range and string_slice::to_pchar_range() explicitly if needed)

/// Auxiliary data structure for parsing and formatting element names (with an
/// optional number tag) into an element number combined with a numeric tag.
/// For completeness there is a `operator <`, too, in case these are needed
/// as map/set keys.
struct FIoElementAndTag {
   int
      /// element number (e.g., 1->H, 6->C, ...)
      iElement;
   FAtomTag
      /// extra numerical tag (normally FAtomTag = int or double) which can be
      /// provided to indicate a specific subtype of an atom/element (e.g., in an
      /// .xyz file, having both C1 and C2 types of carbons, and assigning
      /// different parameter sets/basis sets to them). Defaults to ATOMTAG_None
      /// if not present.
      Tag;

   FIoElementAndTag()
      : iElement(0), Tag(ATOMTAG_None)
   {}

   explicit FIoElementAndTag(int iElement_, FAtomTag const &Tag_ = ATOMTAG_None)
      : iElement(iElement_), Tag(Tag_)
   {}

   explicit FIoElementAndTag(char const *pFirst, char const *pLast) { Parse(pFirst, pLast); };
   // ^-- can pass &*string_slice.first and &*string_slice.last for that

   explicit FIoElementAndTag(std::string const &s) { Parse(s); };
   explicit FIoElementAndTag(pchar_range const &pcr) { Parse(pcr); };

   /// generates a std::string with element symbol and appended numerical tag (if present)
   std::string Format() const;
   /// returns (0-terminated) c-str of element symbol (e.g., "Ne" for iElement == 10)
   char const *Element() const;

   bool operator < (FIoElementAndTag const &other) const {
      if (iElement < other.iElement) return true;
      if (other.iElement < iElement) return false;
      return Tag < other.Tag;
   }
public:
   // parse the string slice in [first,last) as an element+tag combination.
   void Parse(char const *pFirst, char const *pLast);
   void Parse(std::string const &s) { return this->Parse(&s[0], &s[0] + s.size()); }
   void Parse(pchar_range const &pcr) { return this->Parse(pcr.first, pcr.second); }
};

std::string CombineElementAndAtomTag(int iElement, FAtomTag const &Tag);





} // namespace ct

#endif // CX_RAW_ATOM_H
