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

#ifndef IV_ORBITAL_FILE_H
#define IV_ORBITAL_FILE_H

// Code for loading export formats of quantum chemistry programs:
//
// - Molpro xml (preferred format)
// - Molden (worst format)


#include <map>
#include <string>
#include <stdexcept>
#include <vector>
#include "CxTypes.h"
#include "CxPodArray.h"
#include "CtBasisSet.h"
#include "CtAtomSet.h"

namespace orbital_file {

class FFileLoadError : public std::runtime_error
{
public:
   typedef std::runtime_error
      FBase;
   FFileLoadError(std::string const &What, std::string const &Where = "");
};

class FFileTypeUnrecognizedError : public FFileLoadError
{
public:
   FFileTypeUnrecognizedError(std::string const &What, std::string const &Where = "");
};

enum FLoadFlags {
   LOADFILE_SkipVirtualOrbs = 0x0001
};

struct FLoadOptions {
   FLoadOptions(unsigned Flags) : m_Flags(Flags) {};

   bool SkipVirtuals() const { return bool(m_Flags & LOADFILE_SkipVirtualOrbs); }
protected:
   unsigned m_Flags;
};

// keeps a list of the stuff stored under `[Molpro variables]`
struct FVariableSet
{
   typedef std::map<std::string, std::string>
      FVariableMap;
   FVariableMap
      // keeps assignments of all variables -- in string format.
      Map;
   // return if a variable is assigned.
   bool Has(std::string const &Name);
   // try to convert the content of a variable to a float and to
   // return it. raises var_not_found_error if not possible.
   double GetFloat(std::string const &Name);
};



typedef ct::TIntrusivePtr<ct::FAtomSet>
   FAtomSetPtr;
typedef ct::TArray<double>
   FScalarData;
using ct::FBasisSetPtr;

enum FOrbitalSpin {
   ORBSPIN_SpinFree,
   ORBSPIN_Alpha,
   ORBSPIN_Beta,
   ORBSPIN_Unknown
};

char const *sOrbitalSpin(FOrbitalSpin Spin);

struct FOrbitalInfo : public ct::FIntrusivePtrDest
{
   double
      fEnergy;
   double
      fOcc;
   int
      iSym;
   FScalarData
      Orb;
   FBasisSetPtr
      pBasisSet;
   std::string
      sDesc;
   FOrbitalSpin
      Spin;
   std::string Desc() const { return sDesc; };
   std::string FullDesc() const;
   FOrbitalInfo();
};

typedef ct::TIntrusivePtr<FOrbitalInfo>
   FOrbitalInfoPtr;

struct FOrbitalSet : public ct::FIntrusivePtrDest
{
   typedef std::vector<FOrbitalInfoPtr>
      FOrbitalInfoList;
   // ^- FIXME: make this a FOrbitalInfoPtr.
   // note also: in principle there could be more than one orbital set in a file
   // (e.g., for multiple states in MCSCF)
   FOrbitalInfoList
      OrbInfos;
   FOrbitalSpin
      // set to ORBSPIN_Unknown by default; in particular, whenever *NOT* all OrbInfos[i]->Spin are identical.
      CommonSpin;
   FOrbitalSet() : CommonSpin(ORBSPIN_Unknown), MethodDesc(""), OrbType("") {}

   std::string
      // e.g., RKS/UHF/MCSCF
      MethodDesc,
      // e.g., 'canonical', 'natural', 'local(ibo)'
      OrbType;
};

typedef ct::TIntrusivePtr<FOrbitalSet>
   FOrbitalSetPtr;


struct FMolproXmlData : public ct::FIntrusivePtrDest
{
   FOrbitalSetPtr
      pOrbSet;
   FAtomSetPtr
      pAtoms;
   FBasisSetPtr
      pBasisOrb;
   bool
      SemiEmpiricalMode; // no full basis sets.

   FMolproXmlData();
};

typedef ct::TIntrusivePtr<FMolproXmlData>
   FMolproXmlDataPtr;
typedef std::vector<FMolproXmlDataPtr>
   FMolproXmlDataPtrList;

// FMolproXmlDataPtr LoadMolproXmlFile(std::string const &FileName);
// FMolproXmlDataPtr LoadOrbitalFile(std::string const &FileName, FLoadOptions const &LoadOptions);

// returns data set(s) from FileName, one for each Molecule (could be more than one, depending on file type)
FMolproXmlDataPtrList LoadOrbitalFile(std::string const &FileName, FLoadOptions const &LoadOptions);


} // namespace orbital_file



#endif // IV_ORBITAL_FILE_H
