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

/// This file deals with constructions for loading and storing geometry data
/// from/to files, primarily .xyz files.
///
/// Comments:
/// - The file is meant to carry a minimum of dependencies and offer a maximum
///   of flexibility. In particular, its usage is not restricted to
///   MicroScf/wmme/IboView (or other ab initio programs).
///
/// - For this reason, it does not represent atom symbols as the (iElement,
///   iTag) combinations in the rest of wmme-framework code---indeed, it has no
///   means of parsing atom tags like this, because it does not reference
///   CxAtomData.cpp, and cannot actually map element symbols to/from names!
#ifndef CT_XYZ_FRAME_IO_H
#define CT_XYZ_FRAME_IO_H

#include <stdexcept>
#include <vector>
#include <string>
#include <istream>

#include "CxTypes.h"
#include "CxPodArray.h"
#include "CxVec3.h"


namespace xyz_io {

   typedef ct::TVector3<double>
      FVector3;
   typedef std::string
      FAtomType;

   extern double const g_XyzDistanceScale; // factor between a.u. atomic positions (bohr) and what goes into .xyz files. Normally Angstroms, but the constant changes every few years...
   extern double const g_XyzGradientScale; // factor between a.u. gradients (Hartree/bohr) and what goes into .xyz files. It is not standardized.


   struct FIoAtom
   {
      FAtomType
         Type; // normally a string containing the element name (possibly followed by a substring starting with a non-alphabetic character indicating the subtype of the element, e.g., "Na12")
      FVector3
         vPos; // center of atom, in a_bohr (always in a_bohr)

      FVector3
         // if set: gradient dE/d[x,y,z] of the atom in a.u. (TODO: should this be here?
         // it is an /input. in certain kinds of xyz files, so in this sense it would be
         // an actual property...)
         vGrad;

      FIoAtom();
      FIoAtom(FVector3 const &vPos_, std::string const &AtomType_);

      // return Element part of this->Type. Works only if this->Type actually starts with the element symbol
      std::string Element() const;
//       int iElement() const; // works only if Element() works.
      // returns substring in this->Type from first non-alphabetic character to end of string (e.g., in "Na12" it would return "12")
      std::string Tag(bool ConvertToLowerCase=true) const;

      // if present, returns this->Tag() converted to a number, if not, returns iDefaultTag.
      double NumericTag(double DefaultTag) const;
   private:
      void Init(FVector3 const &vPos_, std::string const &AtomType_);
      // returns the index of first character in this->Type belonging to the
      // tag part of the atom type (assuming the atom type starts with the element symbol),
      // or std::string::npos if no tag part is there.
      std::string::size_type FindTagStartIndex() const;
   };


   struct FXyzIoOptions
   {
      double
         // input coordinates will be multiplied by this number before being
         // entered into .xyz file.
         //
         // Default is xyz_io::g_XyzDistanceScale; this correponds to assuming
         // FIoAtom::vPos is given in abohrs and .xyz coords are given in Angstrom.
         m_OutputCoordsFactor,
         m_OutputGradientFactor;

      FXyzIoOptions();

      // Sets scalar conversion factors between internal position/gradient
      // representation (as in FIoAtom::vPos and FIoAtom::vGrad) and data as written in output .xyz file.
      //
      // Concretely:
      // - the internal repr should *multiplied* by these factors to get the .xyz repr
      // - the .xyz repr should *divided* by these factors to get the internal repr
      //
      // Setting both factors to 1.0 disabled unit conversion.
      void SetOutputFactors(double fCoordsFactor, double fGradientFactor);
   public:
      // convert AtomType as in .xyz file into local representation used in FIoAtom
      std::string GetNormalizedAtomType(std::string const &AtomType_) const;
      // convert xyz position as in .xyz file into local representation used in FIoAtom (i.e., to abohr)
      FVector3 ConvertToInternalCoords(FVector3 const &vPos) const;
      // convert gradient as in .xyz file into local representation used in FIoAtom (i.e., to Hartree/abohr)
      FVector3 ConvertToInternalGradient(FVector3 const &vGrad) const;

      // convert xyz position as in FIoAtom (i.e., in abohr) into .xyz file representation
      FVector3 ConvertToOutputCoords(FVector3 const &vPos) const;
      // convert xyz position as in FIoAtom (i.e., in Hartree/abohr) into .xyz file representation
      FVector3 ConvertToOutputGradient(FVector3 const &vGrad) const;
   };

   // data for modifying behavior of .xyz loading routines.
   //
   // TODO: should I merge this and FXyzPrintOptions into FXyzIoOptions? There is considerable overlap.
   struct FXyzLoadOptions : public FXyzIoOptions
   {
      enum FAtomTypeNormalization {
         ATOMTYPE_ToLower = 0x0001, // if set, convert atom type string to lower case
         ATOMTYPE_Normalize = ATOMTYPE_ToLower | 0x0002,  // if set, convert atom type string to first letter upper case and others lower case
         ATOMTYPE_RemoveTrailingNumbers = 0x0100 // if set, remove trailing numbers from atom types (e.g., "C34" becomes just "C")
      };
      unsigned
         m_LoadFlags;

      // other things to consider adding:
      // - how to interpret caption lines (energies etc)
      // - how to read/convert gradient data
      // - how to read/ignore charges or other info.
      explicit FXyzLoadOptions(unsigned LoadFlags_ = ATOMTYPE_Normalize);
      explicit FXyzLoadOptions(FXyzIoOptions const &Base, unsigned LoadFlags_ = ATOMTYPE_Normalize);

   public:
      // convert AtomType as in .xyz file into local representation used in FIoAtom
      std::string GetNormalizedAtomType(std::string const &AtomType_) const;

      bool HasLoadFlag(unsigned Flag) const { return 0 != (m_LoadFlags & Flag); }
   };


   struct FXyzPrintOptions : public FXyzIoOptions
   {
      enum FXyzPrintFlags {
         XYZPRINT_IncludeGradients = 0x01,
         XYZPRINT_IncludeGradientNorm = 0x02,
         XYZPRINT_Default = XYZPRINT_IncludeGradients | XYZPRINT_IncludeGradientNorm
      };
      unsigned
         m_PrintFlags;
      std::string
         // will be added to front of all output lines (typically a set of spaces)
         m_OutputLineLead;

      explicit FXyzPrintOptions(unsigned PrintFlags_ = XYZPRINT_Default);
      explicit FXyzPrintOptions(FXyzIoOptions const &Base, unsigned PrintFlags_ = XYZPRINT_Default);
   public:
      bool HasPrintFlag(unsigned Flag) const { return 0 != (m_PrintFlags & Flag); }

      std::string MakeDescLine(std::string const &FrameDesc, std::string const &ExtraDesc, double Energy, double Gradient) const;
   private:
      void InitPrintOptions(unsigned PrintFlags_);
   };



   struct FIoAtomSet : public ct::FIntrusivePtrDest
   {
   public:
      FIoAtomSet();
      ~FIoAtomSet();

      void AddAtom(FIoAtom const &IoAtom);

      void AddAtomsFromXyzFile(std::string const &FileName, FXyzLoadOptions const *pOtherOptions);
      void AddAtomsFromXyzStream(std::istream &str, std::string const &FileName, FXyzLoadOptions const *pOtherOptions);

      void PrintAsXyz(std::ostream &out, std::string const &ExtraDesc, FXyzPrintOptions const &Options = FXyzPrintOptions()) const;

   public:
      // these are mainly here for simplifying dealing with xyz files describing
      // reaction paths, geometry optimizations and similar things.
      double GetLastEnergy() const { return m_LastEnergy; }
      void SetLastEnergy(double Energy) { m_LastEnergy = Energy; }

      std::string const &GetCaption() const { return m_CaptionLine; }
      void SetCaption(std::string const &Caption_) { m_CaptionLine = Caption_; }

      std::string const &GetName() const { return m_InputName; };
      void SetName(std::string const &InputName_) { m_InputName = InputName_; };

      // returns sum of squared gradient norms of all atoms
      double GetTotalGradientSq() const;
      // returns whether at least one gradient entry in the set is non-zero
      bool HaveGradients() const;

   protected:
      typedef std::vector<FIoAtom>
         FAtomList;
      FAtomList
         m_Atoms;

      double
         m_LastEnergy; // last energy computed for the set
      std::string
         m_InputName; // if read from an .xyz file: remember name of the file here.
      std::string
         m_CaptionLine; // if read from an .xyz file: the xyz caption line.
   public:
      size_t size() const { return m_Atoms.size(); }
      FIoAtom const &operator [] (size_t iAt) const { return m_Atoms[iAt]; }
      FIoAtom &operator [] (size_t iAt) { return m_Atoms[iAt]; }
      void clear() { m_Atoms.clear(); }
      bool empty() const { return m_Atoms.empty(); }

      typedef FAtomList::iterator iterator;
      typedef FAtomList::const_iterator const_iterator;
      iterator begin() { return m_Atoms.begin(); }
      iterator end() { return m_Atoms.end(); }
      const_iterator begin() const { return m_Atoms.begin(); }
      const_iterator end() const { return m_Atoms.end(); }
   };
   typedef ct::TIntrusivePtr<FIoAtomSet>
      FIoAtomSetPtr;
   typedef std::vector<FIoAtomSetPtr>
      FIoAtomSetList;




   /// exception type raised when .xyz loading fails.
   class FXyzLoadException : public std::runtime_error
   {
   public:
      FXyzLoadException(std::string const &Reason, std::string const &FileName);
      FXyzLoadException(std::string const &Reason, std::string const &CurLine, std::string const &FileName);
   };


   /// Load frame data from .xyz file FileName_ and add it to 'Frames'.
   /// This is the main and most flexible .xyz input routine of this module.
   ///
   /// Notes:
   ///
   /// - New geometries are *ADDED* to whatever already is in 'Frames'!
   ///   You may want to clear the array before calling this in case this is not desired.
   ///
   /// - FileName may contain specification of which specific frames to load.
   ///   For example, to select either the single frame 123 from statpt.xyz, one could use:
   ///
   ///       statpt.xyz:123
   ///
   ///   or the following expr to select frames [0,1,2].
   ///
   ///       statpt.xyz:[0,1,2]
   ///
   ///   The part before ':' denotes the xyz file name and
   ///   the part after ':' the frame index/indices (starting at 0).
   void AddFramesFromMultiXyz(FIoAtomSetList &Frames, std::string const &FileName_, FXyzLoadOptions const *pOtherOptions);

   /// simpler-interface loading routine: This one can only load one frame from a .xyz file (it does
   /// not have to be the first one, though; e.g., "statpt.xyz:4" should work!).
   /// The routine raises an exception if less than or more than one frame are contained in the given file specification.
   FIoAtomSetPtr LoadXyzFrameFromFile(std::string const &FileName_, xyz_io::FXyzLoadOptions const *pOtherOptions = 0);
}



#endif // CT_XYZ_FRAME_IO_H
