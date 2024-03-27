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

#include "pugixml.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <vector>
#include "format.h"
#include <sstream>
#include <cctype> // for 1-argument std::tolower.
#include <boost/algorithm/string.hpp> // for boost::split. The only one left...
#include "CxParse1.h"
#include "CtBasisShell.h"
#include "CtMatrix.h"
#include "CxPhysicalUnits.h"
#include "IvOrbitalFile.h"

namespace ctimp {
   void Vec_Ca2Sh(double *pSh, double const *pCa, size_t *ls, size_t N);
   void Vec_CaMolden2CaMolpro(double *pOrb, size_t *ls, size_t N);
   void Vec_ShMolden2ShMolpro(double *pOrb, size_t *ls, size_t N);
   void Vec_Ca2Ca_XSqrtNorm(double *pOrb, size_t *ls, size_t N);
}


namespace orbital_file {


char const *sOrbitalSpin(FOrbitalSpin Spin)
{
   switch (Spin) {
      case ORBSPIN_SpinFree: return "X";
      case ORBSPIN_Alpha: return "A";
      case ORBSPIN_Beta: return "B";
      case ORBSPIN_Unknown:
      default:
         return "?";
   }
}


std::string FOrbitalInfo::FullDesc() const
{
   return fmt::format("{} {{spin:{}, nocc:{:.6f}, eps:{:+.6f}}}", this->Desc(), sOrbitalSpin(this->Spin), this->fOcc, this->fEnergy);
}



// class FXmlParseError : public FFileLoadError
// {
// public:
//    typedef FFileLoadError
//       FBase;
// //    FXmlParseError(std::string const &What, std::string const &Where = "");
// };

FFileTypeUnrecognizedError::FFileTypeUnrecognizedError(std::string const &What, std::string const &Where)
   : FFileLoadError(What, Where)
{}

class FUnexpectedFormatError : public FFileLoadError
{
public:
   typedef FFileLoadError
      FBase;
   FUnexpectedFormatError(std::string const &Message)
      : FBase("Unexpected format: '" + Message)
   {}
};

class FVariableNotFoundError : public FFileLoadError
{
public:
   typedef FFileLoadError
      FBase;
   FVariableNotFoundError( std::string const &VarName )
      : FBase("Variable not found: '" + VarName + "'")
   {}
};

static std::string FmtParseError(std::string const &What, std::string const &Where)
{
   if (Where.empty())
      return What;
   else {
      std::stringstream str;
      str << fmt::format("* ERROR while loading file '{}':\n{}", Where, What) << std::endl;
      return str.str();
   }
}

FFileLoadError::FFileLoadError(std::string const &What, std::string const &Where)
   : FBase(FmtParseError(What,Where))
{
//    raise(SIGTRAP);
//    __builtin_trap();
}


// FXmlParseError::FXmlParseError(std::string const &What, std::string const &Where)
//    : FBase(FmtParseError(What,Where))
// {}


static void str_tolower(std::string::iterator itFirst, std::string::iterator itLast) {
   for (std::string::iterator it = itFirst; it != itLast; ++ it)
      *it = std::tolower(*it);
//    std::transform(str.begin(), str.end(), str.begin(), std::tolower);
//    for (size_t i = 0; i != str.size(); ++ i)
//       *const_cast<char *>(&str[i]) = std::tolower(str[i]);
   // ...
}

static void str_tolower(std::string &str) {
   str_tolower(str.begin(), str.end());
}

// static void str_tolower(ct::string_slice &str) {
//    str_tolower(*const_cast<ct::str_it*>(&str.first), *const_cast<ct::str_it*>(&str.last));
// //    str_tolower(str.first, str.last);
// }

// static void str_lower_and_trim(ct::string_slice &str)
// {
//    str_tolower(str.begin(), str.end());
//    str.trim();
// }

static void str_trim(std::string &sinout) {
   ct::string_slice sl(sinout);
   sl.trim();
   sinout = sl.to_str();
}

static void str_trim_right(std::string &sinout) {
   ct::string_slice sl(sinout);
   sl.trim_right();
   sinout = sl.to_str();
}

// static void str_lower_and_trim(ct::string_slice &str)
// {
//    str_tolower(str.begin(), str.end());
//    str.trim();
// }

static void str_lower_and_trim(std::string &str)
{
   str_tolower(str.begin(), str.end());
   str_trim(str);
}


FOrbitalInfo::FOrbitalInfo()
{
   fEnergy = 0.;
   fOcc = 0.;
   iSym = 0;
   Spin = ORBSPIN_Unknown;
}

enum FLineReadFlags {
   LINE_Trim = 0x01,           // if set, output lines are trimmed on both left and right sides (otherwise only right)
   LINE_ToLower = 0x02,        // if set, output lines are converted to lower case
   LINE_KeepEmpty = 0x04,      // if set, empty lines are not skipped
   LINE_AcceptEof = 0x08,      // if set, EOF will not raise an exception, but just return 'false'.
   LINE_SkipTurboCommentLines = 0x10  // if set, lines starting with '#' will be ignored.
};

static bool get_next_line(std::string &line, std::istream &in, unsigned flags=0) {
   // read first line. skip empty lines.
   line = "";
   while (in.good() && line.empty()) {
      if (!std::getline(in, line)) {
         if (flags & LINE_AcceptEof)
            return false;
         else
            throw FUnexpectedFormatError("Unexpected end of file.");
      }
//       ct::string_slice
//          slLine(line.begin(), line.end());
      if (flags & LINE_ToLower)
         str_tolower(line);
      if (flags & LINE_Trim)
         str_trim(line);
//          slLine.trim();
      else
         str_trim_right(line);
//          slLine.trim_right();
//       line = slLine.to_str();
      if (flags & LINE_KeepEmpty)
         break;
      if (flags & LINE_SkipTurboCommentLines) {
         if (!line.empty() && line[0] == '#') {
            line = "";
         }
      }
   }
   return true;
}

bool ends_with(std::string const &str, std::string const &suffix)
{
   if (str.size() < suffix.size())
      return false;
   return str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}


bool starts_with(std::string const &str, std::string const &suffix)
{
   if (str.size() < suffix.size())
      return false;
   return str.compare(0, suffix.size(), suffix) == 0;
}


static void replace_all_1(std::string &target, std::string const &from, std::string const &to)
{
//    boost::algorithm::replace_all(target, from, to);
//    return;

   if (!(from.size() == to.size() || to.empty()))
      throw std::runtime_error("replace_all_1 is kind of stupid and can only replace strs by either same-length strings or empty strings.");
   if (from.size() == 0)
      throw std::runtime_error("replace_all_1: 'from' str cannot be empty.");
   if (to.size() == from.size()) {
      // replace in-place.
      size_t
         n = from.size();
      for (size_t i0 = 0; i0 + n < target.size(); ++ i0) {
         bool eq = true;
         for (size_t k = 0; k != n; ++ k)
            if (target[i0+k] != from[k])
               eq = false;
         if (eq)
            for (size_t k = 0; k != n; ++ k)
               target[i0+k] = to[k];
      }
   } else if (to.size() == 0) {
//       std::string sincpy = target;
      size_t
         i0 = 0,
         i1 = 0;
      size_t
         n = from.size();
      assert(n != 0); // should have been filtered out above.
      for ( ; i0 < target.size(); ) {
         bool eq = false;
         // is target[i0:i0+n] equal to the 'from' string?
         if (i0 + n < target.size()) {
            eq = true;
            for (size_t k = 0; k != n; ++ k)
               if (target[i0+k] != from[k])
                  eq = false;
         }
         if (eq) {
            // yes, skip this part.
            i0 += n;
         } else {
            // no, copy this part.
            target[i1] = target[i0];
            i0 += 1;
            i1 += 1;
         }
      }
      target.resize(i1);
//       std::cout << fmt::format("replace1: |{}|\n .. --> ..|{}| (delete |{}|)", sincpy, target, from) << std::endl;
//    boost::algorithm::replace_all(target, from, to);
//    return;
   } else {
      throw std::runtime_error("replace_all_1: should not have arrived here.");
   }
}


void replace_fortran_exponents(std::string &line, bool AssumeFloat=false)
{
   // replace fortran exponent format by rest-of-the-world exponent format.
   if (AssumeFloat) {
      replace_all_1(line, "d", "e");
      replace_all_1(line, "D", "E");
   } else {
      replace_all_1(line, "D+", "e+");
      replace_all_1(line, "D-", "e-");
      replace_all_1(line, "d+", "e+");
      replace_all_1(line, "d-", "e-");
   }
}

bool try_read_fortran_float(double &f, std::string s)
{
   replace_fortran_exponents(s, true);
   std::stringstream ss2(s);
   std::string MaybeFloat;

   // this is indeed the most compact C++ way I could come up with
   // for checking if something "cleanly" converts to 'double'...
   // (and even this is not really clean, because strtod uses "errno",
   // which is not threadsafe).
   ss2 >> MaybeFloat;
   char const *pBeg = MaybeFloat.c_str();
   char *pEnd;
   f = std::strtod(pBeg, &pEnd);
   if (pEnd - pBeg == ptrdiff_t(MaybeFloat.size()))
      return true;
   return false;
}

double read_fortran_float(std::string const &s)
{
   double f;
   if (!try_read_fortran_float(f,s))
      throw FUnexpectedFormatError("Expected float, but found '" + s + "'.");
   return f;
}


enum FOrbitalSource {
   ORBSOURCE_Xml,
   ORBSOURCE_Molden
};

static void FixOrbitalTypes1(FOrbitalSet::FOrbitalInfoList &OrbList, FOrbitalSource Source)
{
   size_t
      nAlpha = 0,
      nBeta = 0,
      nClosed = 0,
      nOther = 0,
      nNonIntegerOcc = 0;
   FOrbitalSet::FOrbitalInfoList::iterator
      it;
   for (it = OrbList.begin(); it != OrbList.end(); ++ it) {
      switch((*it)->Spin) {
         case ORBSPIN_Alpha: nAlpha += 1; continue;
         case ORBSPIN_Beta: nBeta += 1; continue;
         case ORBSPIN_SpinFree: nClosed += 1; continue;
         default: nOther += 1; continue;
      }
      if (std::abs(std::modf((*it)->fOcc, 0) > 1e-8))
         nNonIntegerOcc += 1;
   }

   bool RebuildWithOccupations = false;
   if (Source == ORBSOURCE_Molden) {
      // there are *only* Alpha orbitals. Most likely we got something R(o)HFy.
      // Adjust spins based on occupation numbers.
      if (nBeta == 0 && nOther == 0 && nClosed == 0)
         RebuildWithOccupations = true;
   } else if (Source == ORBSOURCE_Xml) {
      // if there are no non-integer occupations, we probably also have a ROHFy wave function,
      // except if we actually got a UKS/UHF one (in which case we need to leave the virtual
      // alpha/beta orbitals alone)
      if (nNonIntegerOcc == 0 && !((nAlpha != 0 || nBeta != 0) && nClosed == 0))
         RebuildWithOccupations = true;
   }

   if (RebuildWithOccupations) {
      for (it = OrbList.begin(); it != OrbList.end(); ++ it) {
         double fOcc = (*it)->fOcc;
         if (fOcc == 2. || fOcc == 0.)
            (*it)->Spin = ORBSPIN_SpinFree;
         if (fOcc == 1. && (*it)->Spin != ORBSPIN_Beta)
            (*it)->Spin = ORBSPIN_Alpha;
      }
   }
}


namespace molpro_xml {

typedef ::orbital_file::FFileLoadError
   FXmlParseError;

using ct::TArray;

void EnumerateChildren(std::ostream &out, pugi::xml_node Node)
{
   out << "Node '" << Node.name() << "':\n";
   for (pugi::xml_node_iterator it = Node.begin(); it != Node.end(); ++it)
   {
      out << "-> '" << it->name() << "' Attr: [";

      for (pugi::xml_attribute_iterator ait = it->attributes_begin(); ait != it->attributes_end(); ++ait)
      {
         if (ait != it->attributes_begin())
            out << " ";
         out << ait->name() << "=" << ait->value();
      }
      out << "]";

      out << std::endl;
   };
};

// append value list contained in a xml node to output list.
template<class FType>
void AppendToList(TArray<FType> &Out, std::string const &s)
{
//    std::cout << "rd list: " << s << std::endl;
   std::stringstream
      ss(s);
   while (1) {
      FType fValue;
      ss >> fValue;
      if (ss.fail())
         break;
      Out.push_back(fValue);
//       std::cout << format("  Val: %10.5f\n") % fValue;
   }
}

static void AppendAssociationList(TArray<int> &Out, std::string sText, std::string const StartHere)
{
   // example strings:
   //   #xpointer(//molecule[@id='2AA7']//basisGroup[@id='1' or @id='2' or @id='3'])]
   //   #xpointer(//molecule[@id='2AA7']//atom[@id='1' or @id='3' or @id='5' or @id='7' or @id='9' or @id='11'])]
   // skip everything until '//basisGroup['. Then erase "@id" and similar things.
   sText.erase(0, sText.find(StartHere.c_str())+StartHere.size());
   sText.erase(sText.find("]"));
   replace_all_1(sText, "@id='", "");
   replace_all_1(sText, " or", "");
   replace_all_1(sText, "'", "");
   replace_all_1(sText, "a", "");
//    std::cout << format("sText: '%s'\n") % sText;
   // should now just be a linear list of integers.
   AppendToList<int>(Out, sText);
}

template<class T>
static std::ostream &PrintArray(std::ostream &out, T const &Array, std::string const &Caption, std::string const &Fmt)
{
   out << Caption << ": ";
   for (typename T::const_iterator it = Array.begin(); it != Array.end(); ++ it) {
      if (it != Array.begin())
         out << " ";
      out << fmt::format(Fmt, *it);
   }
   return out;
}

ct::FBasisSetPtr LoadBasisSet(pugi::xml_node BasisNode, ct::FAtomSet const &Atoms)
{
   std::string
      Name = BasisNode.attribute("id").value(),
      AngularType = BasisNode.attribute("angular").value();
//    uint
//       nShells = BasisNode.attribute("groups").as_int();
//    std::cout << format("Loading basis '%s' with %i groups.\n") % Name % nShells;
   // gna... need to convert to spherical anyway.
//    if (AngularType != "spherical")
//       throw FXmlParseError("Sorry, can only use spherical harmonic basis sets, not cartesians.");
//    ct::FBasisSetPtr pBasis(new F
//    pugi::xml_node Node

   // load the shell information: only unique ones are stored.
   typedef std::map<int, ct::FAtomShellPtr>
      FIdToShellMap;
   FIdToShellMap
      // id -> object.
      AtomShells;
   for (pugi::xml_node GroupNode = BasisNode.child("basisGroup"); GroupNode; GroupNode = GroupNode.next_sibling("basisGroup"))
   {
//       EnumerateChildren(std::cout, GroupNode);
      TArray<double>
         Exp, Co;
      int
         AngMom = GroupNode.attribute("minL").as_int(),
         nExp = GroupNode.attribute("primitives").as_int(),
         nCo = GroupNode.attribute("contractions").as_int(),
         iGroup = GroupNode.attribute("id").as_int();
      if (AngMom != GroupNode.attribute("maxL").as_int())
         throw FXmlParseError("basis contains minL != maxL shells. Cannot deal with those.");
      Exp.reserve(nExp);
      Co.reserve(nExp * nCo);
      // read exponents.
      AppendToList<double>(Exp, GroupNode.child_value("basisExponents"));
      // read contraction coefficients.
      int iCo = 1;
      for (pugi::xml_node CoNode = GroupNode.child("basisContraction"); CoNode; CoNode = CoNode.next_sibling("basisContraction"))
      {
//          std::cout << "contraction " << iCo << "\n";
         int iCoId = CoNode.attribute("id").as_int();
         if (iCoId != iCo)
            throw FXmlParseError(fmt::format("error while reading basis group {}:"
               "expected contraction with id {}, but found {}.", iGroup, iCo, iCoId));
         AppendToList<double>(Co, CoNode.child_value());
         iCo += 1;
      }
      if (nExp != (signed)Exp.size() || nExp*nCo != (signed)Co.size())
         throw FXmlParseError(fmt::format("contraction {}: exponents or co-coefficient number does not match up.", iGroup));
      if (AtomShells.find(iGroup) != AtomShells.end())
         throw FXmlParseError(fmt::format("multiple occurrences of group id {}", iGroup));
      // that's all we need. make a atom basis object
      // and remember it.
      ct::FAtomShellPtr
         pFn(new ct::FAtomShell(AngMom, &Exp[0], nExp, &Co[0], nCo));
      AtomShells[iGroup] = pFn;
   }
//    std::cout << fmt::format("read {} unique atom basis shells.\n", AtomShells.size());


   typedef std::multimap<int, int>
      FAtomToShellIdMap;
   typedef std::pair<FAtomToShellIdMap::iterator, FAtomToShellIdMap::iterator>
      FAtomToShellRangeIt;
   FAtomToShellIdMap
      // atom index -> list of shell ids for the atom (stored in AtomShells)
      AtomShellIds;

   // now this part is... more complicated. we get a list of basis/atom
   // associations, in the form of 'list of unique shells' to 'list of atoms'.
   // However, this information is stored the form of cross-referenced xlinks...
//    EnumerateChildren(std::cout, BasisNode.child("associ  ation"));
   for (pugi::xml_node AssocNode = BasisNode.child("association"); AssocNode; AssocNode = AssocNode.next_sibling("association"))
   {
      // I here hope that I can just skip some of this information, and assume
      // that there is only one molecule and basis sets between different
      // molecules are not mixed (technically they could).
      std::string const
         sBases = AssocNode.child("bases").attribute("xlink:href").value(),
         sAtoms = AssocNode.child("atoms").attribute("xlink:href").value();
      // example strings:
      //   #xpointer(//molecule[@id='2AA7']//basisGroup[@id='1' or @id='2' or @id='3'])]
      //   #xpointer(//molecule[@id='2AA7']//atom[@id='1' or @id='3' or @id='5' or @id='7' or @id='9' or @id='11'])]
      // skip everything until '//basisGroup['. Then erase "@id" and similar things.
      TArray<int>
         iShells, iAtoms;
      AppendAssociationList(iShells, sBases, "//basisGroup[");
      AppendAssociationList(iAtoms, sAtoms, "//atom[");
//       PrintArray(std::cout, iShells, "Shell ids", "{}") << " <-> ";
//       PrintArray(std::cout, iAtoms, "Atom ids", "{}") << "\n";

//     <bases xlink:type="locator" xlink:label="bases"
//       xlink:href="#xpointer(//molecule[@id='436']//basisGroup[@id='1' or @id='2' or @id='3'])"/>
//     <atoms xlink:label="atoms"
//       xlink:href="#xpointer(//molecule[@id='436']//atom[@id='a1' or @id='a3' or @id='a5' or @id='a7' or @id='a9' or @id='a11'])"/>


//       std::cout << ""
      for (uint i = 0; i < iAtoms.size(); ++ i) {
         for (uint j = 0; j < iShells.size(); ++ j)
            AtomShellIds.insert(FAtomToShellIdMap::value_type(iAtoms[i], iShells[j]));
         if (iAtoms[i] > (int)Atoms.size() || iAtoms[i] == 0)
            throw FXmlParseError("encountered basis shell association for non-exitant atom.");
      }
   }

   // in principle we now have all the required information to build the basis set (as long
   // as all functions come in the order I imagine they do...)
   // Go through the atoms and make a basis shell for each of its associated atom shells.
   ct::FBasisSet::FBasisShellArray
      Shells;
   for (uint iAt = 0; iAt < Atoms.size(); ++ iAt) {
      ct::FAtom const
         &At = Atoms[iAt];
      FAtomToShellRangeIt itEq = AtomShellIds.equal_range(iAt+1);
      for (FAtomToShellIdMap::iterator it = itEq.first; it != itEq.second; ++ it) {
         FIdToShellMap::iterator
            itSh = AtomShells.find(it->second);
         if (itSh == AtomShells.end())
            throw FXmlParseError("encountered a non-extant atom shell assigned to an atom.");
         ct::FAtomShellPtr
            pFn = itSh->second;
         Shells.push_back(ct::FBasisShell(At.vPos, iAt, pFn));
//             FBasisShell(FVector3 vCenter_, int iCenter_, FAtomShellCptr pAs_)

      }
   }

//       FBasisSet(FBasisShell const *pShells_, uint nShells_, FBasisContext Context_, std::string const &Name);
   return ct::FBasisSetPtr(new ct::FBasisSet(&Shells[0], Shells.size(), ct::BASIS_Orbital, "(external)"));
}


FAtomSetPtr LoadGeometry(pugi::xml_node GeometryNode, double Factor = 1./ct::ToAng)
{
   FAtomSetPtr
      pAtoms = new ct::FAtomSet();
   int
      iAtom = 1; // starts at 1.
   for (pugi::xml_node_iterator it = GeometryNode.begin(); it != GeometryNode.end(); ++it, ++iAtom)
   {
      std::string
         id = it->attribute("id").value(),
         element = it->attribute("elementType").value();
      std::string
         expected_id = "a" + fmt::format("{}", iAtom);
      if (id != expected_id)
         throw FXmlParseError("problem reading atoms: expected atom id '" + expected_id + "', but found '" + id + "'.");
      try {
         double
            x = ct::string_slice(it->attribute("x3").value()).to_float(),
            y = ct::string_slice(it->attribute("y3").value()).to_float(),
            z = ct::string_slice(it->attribute("z3").value()).to_float();
//          pAtoms->Atoms.push_back(ct::FAtom(Factor * ct::FVector3(x,y,z), element, "(external)"));
         pAtoms->AddAtom(ct::FAtom(Factor * ct::FVector3(x,y,z), element, "(external)"));
      } catch(std::runtime_error &e) {
         throw FXmlParseError(fmt::format("Failed to read geometry at xml node: {}\nError: '{}'", it->name(), e.what()));
      }

//       std::cout << fmt::format("!AT  id = {:4}  elem = {:4}  x = {:10.5f} y = {:10.5f} z = {:10.5f}", id, element, x, y, z);
//       std::cout << std::endl;
   }
//    std::cout << *pAtoms;
   return pAtoms;
}


template<class T>
void ReadArrayT(typename ct::TArray<T> &Out, std::istream &str)
{
//    for (;;) {
//       T x;
//       str >> x;
//       if (str.good())
//          Out.push_back(x);
//       else
//          break;
//    }
   for (;;) {
      T x;
      if (str >> x)
         Out.push_back(x);
      else
         break;
   }
   // ^-- UPDATE: 2020-12-08 cgk.. i think .good returns false if the op was successful, but reached EOF...
}


double ReadNodeTextAsDouble(std::string const &s, double DefaultValueIfNotThere)
{
   // FIXME: get rid of lexical casts. They use locale, which is super-wrong! (and they have no sensible error reporting mechanism)
   if (s != "")
      return ct::string_slice(s).to_float();
   else
      return DefaultValueIfNotThere;
}

FOrbitalSetPtr LoadOrbitals(pugi::xml_node OrbitalsNode, ct::FBasisSetPtr pBasisOrb, FLoadOptions const &LoadOptions)
{
   FOrbitalSetPtr
      pOrbSet = new FOrbitalSet();
   uint
      iOrb = 0;
   uint
      OrbsForSym[8] = {0};
   TArray<size_t>
      ls;
   unsigned
      nBfSph = pBasisOrb->nFn(),
      nBfCa = 0; // expected number of cartesian basis functions.
   bool
      ConvertToSpherical = std::string(OrbitalsNode.attribute("angular").value()) != std::string("spherical");
      // ^- that's the default setting...
//    if (ConvertToSpherical)
//       std::cout << "Expecting Cartesian orbital input.\n";
//    else {
//       std::cout << "Expecting spherical orbital input.\n";
//    }
   for (uint iSh = 0; iSh < pBasisOrb->Shells.size(); ++ iSh)
      for (uint iCo = 0; iCo < pBasisOrb->Shells[iSh].pAs->nCo(); ++ iCo) {
         unsigned l = pBasisOrb->Shells[iSh].pAs->AngMom;
         ls.push_back(l);
//          def nCartY(l): return ((l+1)*(l+2))/2
         nBfCa += ((l+1)*(l+2))/2;
//          std::cout << fmt::format("iSh = {}  l = {} nBfCa = {}\n", iSh, l, nBfCa);
      }

   FOrbitalSpin
      OrbitalSetSpin = ORBSPIN_SpinFree;

   std::string
      OrbitalSet_Method(OrbitalsNode.attribute("method").value()),
      OrbitalSet_Type(OrbitalsNode.attribute("type").value()),
      OrbitalSet_Spin(OrbitalsNode.attribute("spin").value());
   str_tolower(OrbitalSet_Method);
   str_tolower(OrbitalSet_Type);
   str_tolower(OrbitalSet_Spin);
   // so... how does UKS orbital export work? In this case molpro makes multiple sets
   // of output orbitals (so different <orbitals> nodes), and in the current case there
   // will be one alpha set and beta set. Interestingly, at least the 2020 development
   // version I have here DOES NOT export the orbitals as type=CANONICAL and spin=ALPHA,
   // but rather as type=ALPHA and spin=closed (??!!!). So we need to hack around this,
   // and optimally also support the way this was probably originally meant to work,
   // in case this is ever fixed.
   // Also, there appears to be no way of telling Molpro to *NOT* store UHF/UKS natural
   // orbitals. Annoying.
   if ((OrbitalSet_Type == "alpha" && OrbitalSet_Spin == "closed") ||
       (OrbitalSet_Spin == "alpha")) {
      OrbitalSetSpin = ORBSPIN_Alpha;
   }
   if ((OrbitalSet_Type == "beta" && OrbitalSet_Spin == "closed") ||
       (OrbitalSet_Spin == "beta")) {
      OrbitalSetSpin = ORBSPIN_Beta;
   }

   if (OrbitalSetSpin == ORBSPIN_Alpha || OrbitalSetSpin == ORBSPIN_Beta)
      pOrbSet->CommonSpin = OrbitalSetSpin;
   else
      pOrbSet->CommonSpin = ORBSPIN_Unknown;

   pOrbSet->MethodDesc = OrbitalSet_Method;
   if ((OrbitalSet_Type == "alpha" || OrbitalSet_Type == "beta"))
      pOrbSet->OrbType = "canonical";
      // ^- note: "alpha" and "beta" are *NOT* orbital types (natural, canonical, local, ...),
      //    despite being exported as such.
      //    That UKS orbitals marked like that is most likely an error in Molpro export.
      //    So if we see that, we just stay with our default guess of "canonical"... which of
      //    course can be terribly wrong, too. But since the data is missing there is no way
      //    to tell.
   else
      pOrbSet->OrbType = OrbitalSet_Type;


   for (pugi::xml_node OrbNode = OrbitalsNode.child("orbital"); OrbNode; ++iOrb, OrbNode = OrbNode.next_sibling("orbital"))
   {
//       EnumerateChildren(std::cout, OrbNode);
      FOrbitalInfoPtr
         pInfo = new FOrbitalInfo;
      try {
         pInfo->fEnergy = ReadNodeTextAsDouble(OrbNode.attribute("energy").value(), 0.);
         pInfo->fOcc = ReadNodeTextAsDouble(OrbNode.attribute("occupation").value(), 2.);
         pInfo->iSym = int(ReadNodeTextAsDouble(OrbNode.attribute("symmetryID").value(), 1.));
      } catch(std::runtime_error &e) {
         throw FXmlParseError(fmt::format("Failed to orbital info at xml node: {}\nError: '{}'", OrbNode.child_value(), e.what()));
      }
      pInfo->pBasisSet = pBasisOrb;
      pInfo->Spin = OrbitalSetSpin;
      if (pInfo->fOcc < 0.00001 && LoadOptions.SkipVirtuals())
         continue;
      // ^- why not break? in symmetry calculations there may be additional
      //    occupied orbitals later.
//       std::cout << fmt::format("ORB: E[{:3}] = {:12.6f} H   Occ = {:12.6f}  Sym = {}\n", iOrb, fEnergy, fOcc, iSym);
//       std::cout << OrbNode.child_value() << std::endl;
      if (pInfo->iSym <= 8) {
         OrbsForSym[pInfo->iSym-1] += 1;
//          pInfo->sDesc = fmt::format("{:4}.{} [E={:8.4f}  O= {:6.4f}]", OrbsForSym[pInfo->iSym-1], pInfo->iSym, pInfo->fEnergy, pInfo->fOcc);
         pInfo->sDesc = fmt::format("{:4}.{}", OrbsForSym[pInfo->iSym-1], pInfo->iSym);
      } else {
//          pInfo->sDesc = fmt::format("{:4}.? [E={:8.4f}  O= {:6.4f}]", iOrb, pInfo->fEnergy, pInfo->fOcc);
         pInfo->sDesc = fmt::format("{:4}.?", iOrb);
      }

      std::stringstream
         str(OrbNode.child_value());
      TArray<double>
         OrbCa;
      ReadArrayT(OrbCa, str);
      pInfo->Orb.resize(nBfSph, 777.777);
      if (ConvertToSpherical) {
         // transform from cartesian back to spherical. (...hopefully...)
         if (OrbCa.size() != nBfCa)
            throw FXmlParseError(fmt::format("import orbital {}: expected {} cartesian components, but got {}.", (iOrb+1), nBfCa, OrbCa.size()));
         ctimp::Vec_Ca2Sh(&pInfo->Orb[0], &OrbCa[0], &ls[0], ls.size());
//          throw FXmlParseError("Cartesian basis input currently not supported. Use {put,xml,...; keepspherical}");
      } else {
         // hopefully already is spherical...
         if (OrbCa.size() != nBfSph)
            throw FXmlParseError(fmt::format("import orbital {}: expected {} spherical components, but got {}.", (iOrb+1), nBfSph, OrbCa.size()));
         pInfo->Orb = OrbCa;
      }
      pOrbSet->OrbInfos.push_back(pInfo);
//       std::cout << pInfo->Desc() << std::endl;
//       std::cout << fmt::format("   nCoeff: {}  (of {})", pInfo->Orb.size(), pBasisOrb->nFn()) << std::endl;
   }
   FixOrbitalTypes1(pOrbSet->OrbInfos, ORBSOURCE_Xml);
   return pOrbSet;
}

int TestOrbitals(ct::FBasisSetPtr pBasis, FAtomSetPtr pAtoms, FOrbitalSetPtr pOrbSet)
{
   using namespace ct;

   ct::FMemoryStack2 Mem(200000000);
   uint
      nAo = pBasis->nFn(),
      nOrb = pOrbSet->OrbInfos.size();
   FStackMatrix
      S(nAo, nAo, &Mem),
      COrb(nAo, nOrb, &Mem);
   pAtoms->MakeOverlapMatrix(S, *pBasis, *pBasis, Mem);
//    S.Print(std::cout, "AO overlap");
   for (uint iOrb = 0; iOrb < pOrbSet->OrbInfos.size(); ++ iOrb) {
      FOrbitalInfo
         &oi = *pOrbSet->OrbInfos[iOrb];
      assert(nAo == oi.Orb.size());
      for (uint iBf = 0; iBf < nAo; ++ iBf)
         COrb(iBf, iOrb) = oi.Orb[iBf];
   }
//    COrb.Print(std::cout, "MO coeffs");
   {
      double rmsd = 0.;
      FStackMatrix
         OrbOvl(nOrb, nOrb, &Mem);
      ChainMxm(OrbOvl, Transpose(COrb), S, COrb, Mem);
      for (size_t i = 0; i < OrbOvl.nRows; ++ i)
         for (size_t j = 0; j < OrbOvl.nCols; ++ j) {
            double f = OrbOvl(i,j) - ((i==j)? 1. : 0.);
            rmsd += f*f;
         }
      rmsd = std::sqrt(rmsd);
      OrbOvl.Print(std::cout, "Orbital MO overlap");
      std::cout << fmt::format("  RMSD from orthogonality: {:8.2e}", rmsd) << std::endl;
   }
   return 0;
}


static void DbgPrintOrbitalSet(FOrbitalSetPtr pThisSet)
{
//    std::cout << fmt::format("* ORBITAL SET: {}/{}/{}  (Note: HaveAlpha={}, HaveBeta={}, SkipThis?={})\n", pThisSet->MethodDesc, pThisSet->OrbType, sOrbitalSpin(pThisSet->CommonSpin), HaveAlphaOrbs, HaveBetaOrbs, SkipThisSet);
   std::cout << fmt::format("* ORBITAL SET: {}/{}/{}\n", pThisSet->MethodDesc, pThisSet->OrbType, sOrbitalSpin(pThisSet->CommonSpin));
   size_t i = 0;
   for (FOrbitalInfoPtr &pOrbInfo : pThisSet->OrbInfos) {
      std::cout << fmt::format("   [{:3}] {}\n", i, pOrbInfo->FullDesc());
      i += 1;
   }
}


FMolproXmlDataPtr LoadMolproXmlFile(std::string const &FileName, FLoadOptions const &LoadOptions)
{
   try {
      FMolproXmlDataPtr
         pXmlData = new FMolproXmlData();
      pugi::xml_document doc;
      pugi::xml_parse_result result = doc.load_file(FileName.c_str());
      if (result) {
//          std::cout << fmt::format("* read '{}'.", FileName) << std::endl;
      } else {
         throw FXmlParseError(fmt::format("XML load result: '{}'", result.description()));
//          return FMolproXmlDataPtr(0);
      }

      pugi::xml_node MoleculeNode = doc.child("molpro").child("molecule");
      pugi::xml_node GeometryNode = MoleculeNode.child("cml:atomArray");
      if (!GeometryNode) {
         // see if we have the geometry as a sub node in a <cml:molecule> node instead
         // of directly under <molecule> as in previous versions of Molpro master
         pugi::xml_node CmlMoleculeNode = MoleculeNode.child("cml:molecule");
         if (CmlMoleculeNode)
            GeometryNode = CmlMoleculeNode.child("cml:atomArray");
      }
      if (!GeometryNode)
         throw FXmlParseError("Failed to find <cml:atomArray> node.");
      pXmlData->pAtoms = LoadGeometry(GeometryNode);

      pXmlData->pAtoms->SetLastEnergy(ReadNodeTextAsDouble(MoleculeNode.attribute("energy").value(), 0.));

      pugi::xml_node OrbBasisNode = MoleculeNode.find_child_by_attribute("basisSet", "id", "ORBITAL");
      if (!OrbBasisNode)
         throw FXmlParseError("Failed to find orbital basis node (<basisSet id='ORBITAL' ...>).");

      pXmlData->pBasisOrb = LoadBasisSet(OrbBasisNode, *pXmlData->pAtoms);
   //    std::cout << *pBasisOrb;

//       pugi::xml_node OrbitalsNode = MoleculeNode.child("orbitals");
//       if (!OrbitalsNode)
//          throw FXmlParseError("Failed to find orbital node (<orbitals ...>).");
//       pXmlData->pOrbSet = LoadOrbitals(OrbitalsNode, pXmlData->pBasisOrb, LoadOptions);

      {
         std::vector<FOrbitalSetPtr>
            OrbitalSets;
         bool
            HaveAlphaOrbs = false,
            HaveBetaOrbs = false;
         for (pugi::xml_node OrbitalsNode = MoleculeNode.child("orbitals"); OrbitalsNode; OrbitalsNode = OrbitalsNode.next_sibling("orbitals"))
         {
            FOrbitalSetPtr
               pOrbSet = LoadOrbitals(OrbitalsNode, pXmlData->pBasisOrb, LoadOptions);
            OrbitalSets.push_back(pOrbSet);
            if (pOrbSet->CommonSpin == ORBSPIN_Alpha)
               HaveAlphaOrbs = true;
            if (pOrbSet->CommonSpin == ORBSPIN_Beta)
               HaveBetaOrbs = true;
         }
         if (OrbitalSets.empty())
            throw FXmlParseError("Failed to find any orbital set node (<orbitals ...>).");

         // next, merge all the relevant orbitals we got into a new orbital set.
         // However... there are some caveats... e.g., for some reason beyond my
         // comprehension Molpro insists on exporting Alpha, Beta, *AND* Natural(!!)
         // orbitals when using put,xml on a UHF/UKS record. And if there is any
         // way to turn it off... well, I did not manage to find it. So, if we *do*
         // have a UKS/UHF state, we need to IGNORE the natural orbitals. They represent
         // the same electrons as the alpha/beta orbitals, after all, and both together
         // would create havoc.
         FOrbitalSetPtr
            pMergedOrbSet(new FOrbitalSet);
         for (FOrbitalSetPtr &pThisSet : OrbitalSets) {
            bool
               SkipThisSet = false;
            if ((HaveAlphaOrbs || HaveBetaOrbs) && pThisSet->OrbType == "natural")
               // skip superfluous natural orbital sets if we already have canonical
               // orbital sets for the individual spins.
               SkipThisSet = true;

            if (0) {
               DbgPrintOrbitalSet(pThisSet);
               std::cout << fmt::format("  (^-- Note: HaveAlpha={}, HaveBeta={}, SkipThis?={})\n", HaveAlphaOrbs, HaveBetaOrbs, SkipThisSet);
            }
            if (SkipThisSet)
               continue;
            pMergedOrbSet->MethodDesc = pThisSet->MethodDesc;
            if (pMergedOrbSet->MethodDesc.empty())
               pMergedOrbSet->MethodDesc = pThisSet->MethodDesc;
            else if (pMergedOrbSet->MethodDesc != pThisSet->MethodDesc)
               std::cerr << fmt::format("\nWARNING: while processing xml file '{}': found multiple orbital sets of methods. These are now merged. That cannot be good...\n", FileName) << std::endl;
            if (pMergedOrbSet->OrbType.empty())
               pMergedOrbSet->OrbType = pThisSet->OrbType;
            else if (pMergedOrbSet->OrbType != pThisSet->OrbType)
               std::cerr << fmt::format("\nWARNING: while processing xml file '{}': found multiple orbital sets of orbital types. These are now merged. That cannot be good...\n", FileName) << std::endl;
            pMergedOrbSet->CommonSpin = ORBSPIN_Unknown;
            pMergedOrbSet->OrbInfos.insert(pMergedOrbSet->OrbInfos.end(), pThisSet->OrbInfos.begin(), pThisSet->OrbInfos.end());
         }

         if (0) {
            DbgPrintOrbitalSet(pMergedOrbSet);
         }

         pXmlData->pOrbSet = pMergedOrbSet;
      }
      return pXmlData;
   } catch(FXmlParseError &e) {
      throw FXmlParseError(e.what(), FileName);
   }
   return 0;
}


} // namespace molpro_xml


namespace molden_file {
   typedef ::orbital_file::FFileLoadError
      FMoldenParseError;

   using ct::TArray;

// the molden files made by different programs look almost the same, but they
// are all slightly incompatible due to various different bugs.
enum FMoldenSource {
   MOLDENSOURCE_Molpro,
   // Molcas thinks the .molden format should have different headings than the Molden spec. E.g. '(AU)' tags instead of 'AU' tags.
   MOLDENSOURCE_Molcas,
   // Turbomole uses a different basis function normalization than Molpro and Molcas
   MOLDENSOURCE_Turbomole,
   // Orca again uses a different normalization than either Molpro/Molcas or Turbomole
   MOLDENSOURCE_Orca,
   // This is for the Grimme GFNx-XTB program's '--molden' option.
   // As of version xtb-190318, special adjustments are required because:
   // - The program prints an empty title section, but does it incorrectly
   //   such that the [atoms]-section header is regarded as title, and therefore
   //   the entire atoms section is skipped.
   // - For the basis functions, it combines the Turbomole bug regarding the solid
   //   harmonic normalization (the sqrt(2l-1) factors) with the Orca bug regarding
   //   the primitive Gaussian normalization
   // - As of GFN2-XTB with version 190318, there are only s,p,d functions in the
   //   XTB basis sets. So the present code has only been checked for up to these.
   // - Additionally, XTB uses floating occupation number SCF, and therefore will
   //   occasionally generate orbital files which do not actually correspond to
   //   any N-electron wave function (this is a feature, though, not a bug).
   //   Therefore some other parts of the code need adjusted regarding wf handling.
   //   Also minimal basis sets and nuclear electron numbers need adjustment,
   //   but this is in other parts of the code.
   MOLDENSOURCE_GrimmeXtb
};


struct FMoldenAtom
{
   ct::FVector3
      vPos; // stored in atomic units
   std::string
      Label; // normally corresponds to element name. But might contain suffices
   int
      // I don't quite get it, but center numbers are explicitly stored.
      // I assume therefore they can be out of order, or with holes, although
      // in the examples I've seen they never were.
      iCenter,
      iElement; // 1: H, 2: He, ...
   int
      // index this atom will get in the loaded file. We require in various places
      // that atoms are numbered continously, starting at 0, and that the basis sets
      // are given in the same order as the atoms. We thus might have to shovel things
      // around...
      iCenterOut;
   std::string MakeDesc() const {
      return fmt::format("FMoldenAtom(Label='{}', iCenter={}, iElement={}, vPos=({:10.5f}, {:10.5f}, {:10.5f}))",
         Label, iCenter, iElement, vPos[0], vPos[1], vPos[2]);
   }
};

// maps center indices to atom objects. Note comments on iCenter
// in FMoldenAtom.
struct FMoldenAtomSet : public std::map<int, FMoldenAtom>
{
   FMoldenAtom &operator [] (int iCenter) {
      return (*this)[iCenter];
   }
   FMoldenAtom const &operator [] (int iCenter) const;

   size_t size() const { return (*this).size(); }
   void AddAtom(FMoldenAtom const &A);

   FMoldenAtom &Get(int iCenter);

   ct::FAtomSetPtr Convert();
};


struct FMoldenGaussShell
{
   ct::TArray<double>
      Co,  // nExp x nCo contraction matrix
      Exp; // nExp primitive exponents
   size_t
      nExp,
      nCo;
   int
      iAngMom,
      iCenter;
//    FMoldenGaussShell() {};
   size_t
      // position of the function in the original basis order. Output basis
      // will be ordered by center index and angular momentum.
      iOffsetIn;
   bool
      Spherical;
   inline size_t nFn() const;
   inline size_t nFnSph() const;

   bool operator < (FMoldenGaussShell const &other) const;
};

struct FMoldenBasis
{
   typedef std::vector<FMoldenGaussShell>
      FGaussFnList;
   FGaussFnList
      Shells;
   size_t nFn() const; // as declared in actual shell objects
   size_t nFnSph() const; // as if all shells were declared as spherical
   size_t nFnCart() const; // as if all shells were declared as cartesian

   // convert Molden's basis format into our own by rebuilding the general contractions and
   // joining equivalent shell functions into common FBasisFn objects and deferred FBasisShell objects.
   FBasisSetPtr Convert(FMoldenAtomSet *pMoldenAtomSet);

   TArray<size_t>
      // iOrderedFn = m_InputToOutputOrder[iInputFn].
      // Make by this->Sort(). This takes care only of the center/AM re-alignment.
      // individual function orders are dealt with separately.
      m_InputToOutputOrder;
   // sort basis functions by center index and angular momentum.
   void Sort();

   bool Sorted() const { return m_IsSorted; };
   FMoldenBasis() : m_IsSorted(false) {}
private:
   bool m_IsSorted;
};

// // that's just the same thing... as FVariableSet. Apparently a random list of
// // stuff we can assign to various properties associated with an orbital. Not
// // sure if they are supposed to be predefined. Molpro stores stuff like 'Ene'
// // (orbital energy), 'Occup' (orbital occupancy), or 'Spin' (for unknown reason
// // ``Alpha'' for closed orbitals).
// struct FMoldenOrbInfo : public FVariableSet
// {
//    FScalarData
//       Coeffs;
// };


struct FMoldenFile
{
   FVariableSet
      Vars;
   FMoldenAtomSet
      Atoms;
   FMoldenBasis
      Basis;

   typedef std::vector<FOrbitalInfoPtr>
      FOrbInfoList;
   FOrbInfoList
      // information provided in molden file on the individual molecular orbitals.
      // There is one object for each MO in Orbs, even if no data is associated with
      // the orbital.
      OrbInfo;
//    size_t
//       // number of *spherical* basis functions and of molecular orbitals
//       nAo, nOrb;
//    ct::TArray<double>
//       // nAo x nOrb matrix.
//       Orbs;

   bool
      m_InputOrbitalsAreSpherical,
      m_InputOrbitalsAreSphericalForG,
      m_AllowLargerThanG;
   FMoldenSource
      m_MoldenSource;
   FLoadOptions const
      *m_pLoadOptions;

   explicit FMoldenFile(std::istream &in, FLoadOptions const *pLoadOptions_);
private:
   void SkipCurrentSection(std::istream &in);
   void AddAtom(FMoldenAtom const &A);
   void ReadAtomsSection(std::istream &in, double fCoordFactor);
   void ReadBasisSection(std::istream &in);
   void ReadMoSection(std::istream &in);
   void FixOrbitalTypes();
   void FixBasis(); // fix spherical declarations and basis function order.
};

// number of cartesians components with angular momentum == l
static size_t nCartY(int l) { return static_cast<size_t>((l+1)*(l+2)/2); }
// number of solid harmonic components with angular momentum == l
static size_t nSlmY(int l) { return static_cast<size_t>(2*l+1); }

static size_t const npos = std::string::npos;

size_t FMoldenGaussShell::nFn() const
{
   if (Spherical)
      return nSlmY(iAngMom);
   else
      return nCartY(iAngMom);
}

size_t FMoldenGaussShell::nFnSph() const {
   return nSlmY(iAngMom);
}

size_t FMoldenBasis::nFnSph() const
{
   size_t
      r = 0;
   FGaussFnList::const_iterator
      it;
   _for_each(it, Shells)
      r += nSlmY(it->iAngMom);
   return r;
}

size_t FMoldenBasis::nFnCart() const
{
   size_t
      r = 0;
   FGaussFnList::const_iterator
      it;
   _for_each(it, Shells)
      r += nCartY(it->iAngMom);
   return r;
}

size_t FMoldenBasis::nFn() const
{
   size_t
      r = 0;
   FGaussFnList::const_iterator
      it;
   _for_each(it, Shells)
      r += it->nFn();
   return r;
}

FMoldenAtom &FMoldenAtomSet::Get(int iCenter) {
   iterator
      itAtom = this->find(iCenter);
   if ( itAtom == this->end() )
      throw std::runtime_error(fmt::format("atomic center not found: {}", iCenter));
   return itAtom->second;
}

FMoldenAtom const &FMoldenAtomSet::operator [] (int iCenter) const {
   const_iterator
      itAtom = this->find(iCenter);
   if ( itAtom == this->end() )
      throw std::runtime_error(fmt::format("atomic center not found: {}", iCenter));
   return itAtom->second;
}

FMoldenFile::FMoldenFile(std::istream &in, FLoadOptions const *pLoadOptions_)
   : m_InputOrbitalsAreSpherical(false), // that's the default
     m_InputOrbitalsAreSphericalForG(false), // that's also a default---independent of the d and f functions.
     m_AllowLargerThanG(false),
     m_MoldenSource(MOLDENSOURCE_Molpro),
     m_pLoadOptions(pLoadOptions_)
{
   std::string
      line;
   get_next_line(line, in, LINE_Trim | LINE_ToLower);
   if (!in.good() || line != "[molden format]")
      throw FUnexpectedFormatError("expected '[Molden Format]'");
   for (; in.good(); ) {
      // read section header
      if (!get_next_line(line, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof))
         break; // the end.
      // std::cout << fmt::format("\n NEXT MOLDEN SECTION: '{}'\n\n", line);
      if (line == "[molpro variables]") {
         SkipCurrentSection(in);
         continue;
      } else if (line == "[title]") {
         // FIXME: Turbomole molden files have a "[Title]" section and Molpro
         // files don't. Molcas files also don't. Until recently Molpro files
         // also  "[Molpro variables]", but these were removed at some point in
         // late 2014 (and, apparently, this was also backported to 2012.1...).
         // Not good! Need to find a better way of telling the programs apart.
         std::string
            title;
         get_next_line(title, in, LINE_Trim | LINE_ToLower | LINE_KeepEmpty);

         // check if what we just read looks like a section header, rather
         // than a title. Grimme's xtb program emits a '[Title]' section,
         // which is followed *directly* by '[Atoms] AU' in the next line, rather
         // than an empty title. This causes the entire section to be skipped
         // if not handled explicitly, because it would be read as a file title
         // being '[Atoms] AU' followed by a load of junk.
         // (it also means we can't have file titles starting with square brackets,
         // but I guess that is acceptable..)
         bool ThisIsASectionHeader = false;
         if (!title.empty() && title[0] == '[') {
            // don't process this line as title. But rather fall through
            // and proces the title line as next section header.
            ThisIsASectionHeader = true;
            line = title;
            // we're most likely here from Grimme's xtb program. But we just chose a
            // set of basis normalization options compatible with what it appears to
            // export.
            m_MoldenSource = MOLDENSOURCE_GrimmeXtb;
         }

         if (!ThisIsASectionHeader) {
            // files generated by orca_2mkl xxx -molden have a tag saying
            // "created by orca_2mkl" in the title. Very helpful! So we now apply some heuristics...
            // if there is a title and it says "Hi, I'm from Orca" we apply orca hacks, otherwise
            // we apply Turbomole hacks.
            if (std::string::npos != title.find("created by orca_2mkl")) {
               m_MoldenSource = MOLDENSOURCE_Orca;
   //             std::cout << "Molden file format: Orca" << std::endl;
            } else if (std::string::npos != title.find("molden2aim,")) {
               // dunno... let's go with Molpro?  ... yes, that seems to work.
               m_MoldenSource = MOLDENSOURCE_Molpro;
            } else {
               // note: at least in my case the title is empty for files made by turbomole.
               // but I wouldn't bet on it staying that way.
               m_MoldenSource = MOLDENSOURCE_Turbomole;
   //             std::cout << "Molden file format: Turbomole (title='" << title << "')" << std::endl;
            }
            SkipCurrentSection(in);
            continue;
         }
      };

      if (line == "[atoms] angs") {
//          std::cout << "Read atoms: ANGS" << std::endl;
         ReadAtomsSection(in, 1./ct::ToAng);
      } else if (line == "[atoms] au" || line == "[atoms] (au)") { // second variant used off-spec by Molcas.
//          std::cout << "Read atoms: AU" << std::endl;
         ReadAtomsSection(in, 1.);
      } else if (line == "[gto] (au)") {
         // Molcas files can be identified by the parenthesis around "au".
         // There are not supposed to be any according to spec. Especially not after "[gto]".
         m_MoldenSource = MOLDENSOURCE_Molcas;
         ReadBasisSection(in);
      } else if (line == "[gto]") {
         ReadBasisSection(in);
      } else if (line == "[5d]" || line == "[5d7f]" || line == "[7f]") {
         m_InputOrbitalsAreSpherical = true;
         // ^- WARNING: Orca emits these AFTER emitting the basis.
         // That means that we can only fix the basis right before reading the MOs,
         // not after reading the basis as I'd do otherwise.
      } else if (line == "[5d10f]]") {
         throw FMoldenParseError("Mixed cartesian and spherical input functions are not supported. Please use *either* cartesian *or* spherical functions, but not both!");
      } else if (line == "[9g]") {
         m_InputOrbitalsAreSphericalForG = true;
      } else if (line == "[allow_larger_g]") {
         // Note: This is not an "officially" sactioned tag. But there is little reason
         // to prevent this if both IboView and the exporting program can do it, only
         // because Molden can't.
         // With this we will allow larger-than-g functions to be
         // specified. For the spherical versions the order is not very
         // involved and most likely right. For the cartesians we assume
         // Molpro order (which for d and g functions miraculously seems
         // to be equal to Molden order).
         m_AllowLargerThanG = true;
      } else if (line == "[mo]") {
         if (m_InputOrbitalsAreSphericalForG && !m_InputOrbitalsAreSpherical)
            throw FMoldenParseError("Input specifies spherical g functions, but does not specify spherical lower functions. Mixed cartesian and spherical functions are not supported by this loader.");
         // (^- note: can't check the other error variant here (non-spherical g
         //  and spherical d/f), since the input might not have any g functions,
         //  thus there would not be any need to specify them as spherical)
         FixBasis();
         ReadMoSection(in);
      } else if (line.find('[') == npos) {
         throw FMoldenParseError("Expected section header, but found '" + line + "'");
      } else
         SkipCurrentSection(in);
   }
}


static bool IsSectionStart(std::string const &line)
{
   size_t iopen = line.find('[');
   if (iopen != npos && line.find(']',iopen+1) != npos)
      return true;
   return false;
}


void FMoldenFile::SkipCurrentSection(std::istream &in)
{
   // skipp content of current section. That is, everything until the next section
   // header (which should be indicated by a '[..]' in the line)
   std::string
      line;
   while (in.good()) {
      std::streampos pos = in.tellg();
      if (!get_next_line(line, in, LINE_AcceptEof))
         return;
      if (IsSectionStart(line)) {
         // rewind to start of line.
         in.seekg(pos);
         return;
      }
   }
}

void FMoldenFile::AddAtom(FMoldenAtom const &A)
{
   FMoldenAtomSet::iterator
      itAtom = Atoms.find(A.iCenter);
   if ( itAtom != Atoms.end() )
      throw FMoldenParseError(fmt::format("center {} already assigned to another atom", A.iCenter));
   Atoms.insert(itAtom, FMoldenAtomSet::value_type(A.iCenter,A));
}

void FMoldenFile::ReadAtomsSection(std::istream &in, double fCoordFactor)
{
   std::string
      line;
   while (in.good()) {
      std::streampos pos = in.tellg();
      if (!get_next_line(line, in, LINE_Trim | LINE_AcceptEof))
         break;
      if (IsSectionStart(line)) {
         // rewind to start of line.
         in.seekg(pos);
         break;
      }
      // expected format:
      //  C     1    6         0.0000000000        0.0000000000        0.0000000000
      //             ^- charge?
      //        ^- center index (1 based)
      //  ^- label?
      FMoldenAtom
         At;
      std::stringstream ss(line);
      ss >> At.Label >> At.iCenter >> At.iElement >> At.vPos[0] >> At.vPos[1] >> At.vPos[2];
      At.vPos *= fCoordFactor;
//       bool is_bad = ss.bad();
//       bool is_fail = ss.fail();
//       bool is_eol = ss.eof();
//       bool is_good = ss.good();
      if (ss.bad() || ss.fail())
         throw FMoldenParseError("Failed to understand atom declaration '" + line + "'");
      AddAtom(At);
//       std::cout << fmt::format("   read atom: {}\n", At.MakeDesc());
   }

   // make a new set of indices which is ordered continously and comes without holes.
   size_t
      iCenterOut = 0;
   for (FMoldenAtomSet::iterator itAtom = Atoms.begin(); itAtom != Atoms.end(); ++ itAtom) {
      itAtom->second.iCenterOut = iCenterOut;
      ++ iCenterOut;
   }
}

int iAngMomFromStr(std::string const &s_) {
   std::string s(s_);
   str_tolower(s);
   if (s == "sp")
      throw FMoldenParseError("This parser does not supported mixed 'sp' shells. Please use a different basis set (e.g., def2 sets).");
   if (s.size() != 1)
      throw FMoldenParseError("Expected angular momentum, but found '" + s + "'.");
   static char const pAngMomNames[] = "spdfghikl";
   for (int i = 0; size_t(i) < sizeof(pAngMomNames)/sizeof(pAngMomNames[0]); ++ i)
      if (pAngMomNames[i] == s[0])
         return i;
   throw FMoldenParseError("'" + s + "' is not a recognized angular momentum.");
}

static double OrcaRenormFactor(double fExp, unsigned l)
{
   // note: this differs from RawGaussNorm(fExp,l) by a missing DoubleFactR(2*l-1) inside the square root.
   // I guess they just use a different solid harmonics definition.
//    return pow(M_PI/(2*fExp),.75) * sqrt(1./pow(4.*fExp,l));

   double f = pow(M_PI/(2*fExp),.75) * sqrt(1./pow(4.*fExp,l));
   if (l == 4) return f * std::sqrt(3); // .oO(what.)
   return f;
}


void FMoldenFile::ReadBasisSection(std::istream &in)
{
   Basis.Shells.clear();

   size_t
      iInputBasisOffset = 0;
   std::string
      line;
   while (in.good()) {
      std::streampos pos = in.tellg();
      if (!get_next_line(line, in, LINE_Trim | LINE_AcceptEof))
         break;
      if (IsSectionStart(line)) {
         // rewind to start of line.
         in.seekg(pos);
         break;
      }
      // expected format:
      //   1 0                                    // center index, ??? (always 0?)
      //  s    5  1.00                            // number of primitives, ??? (always 1? number of contractions?)
      //   0.1238401694D+04  0.5517365048D-02     // exponents and contraction coefficients.
      //   0.1862900499D+03  0.4108882856D-01
      //   0.4225117635D+02  0.1822538213D+00     // beware of fortran exponent format!! (D+... instead of E+...)
      //   0.1167655793D+02  0.4682845944D+00
      //   0.3593050648D+01  0.4457581734D+00
      //  s    1  1.00
      //   0.4024514736D+00  0.1000000000D+01
      //  s    1  1.00
      //   0.1309018267D+00  0.1000000000D+01
      //  p    3  1.00
      //   0.9468097062D+01  0.5688833201D-01
      //   0.2010354514D+01  0.3129405934D+00
      //   0.5477100471D+00  0.7606501651D+00
      //  p    1  1.00
      //   0.1526861379D+00  0.1000000000D+01
      //  d    1  1.00
      //   0.8000000000D+00  0.1000000000D+01
      //
      int
         iCenter, iDunno;
      {
         std::stringstream
            str1(line);
         str1 >> iCenter >> iDunno;
      }
      while (in.good()) {
         if (!get_next_line(line, in, LINE_Trim | LINE_KeepEmpty | LINE_AcceptEof))
            break;
         if (line.empty())
            break; // end of basis functions for current center.
         FMoldenGaussShell
            Fn = FMoldenGaussShell(),
            *pFn = &Fn;
         {
            std::stringstream
               str2(line);
            std::string
               sAngMom, Dummy;
            str2 >> sAngMom >> pFn->nExp >> Dummy;
            pFn->iAngMom = iAngMomFromStr(sAngMom);
//             if (pFn->iAngMom >= 5 && !m_AllowLargerThanG)
//                throw FMoldenParseError("Molden itself does not support larger-than-g functions. You can ask this program to load them anyway by adding [allow_larger_g] as section to the molden input file.");
//             if (pFn->iAngMom == 4 && m_InputOrbitalsAreSpherical && !m_InputOrbitalsAreSphericalForG)
//                throw FMoldenParseError("Mixed spherical and cartesian functions are not supported. (g was specified as default cartesian, while df were defined as spherical).");
         }
         pFn->nCo = 1; // <- I don't think this can handle general contractions, or can it?
         pFn->Exp.resize(pFn->nExp);
         pFn->Co.resize(pFn->nCo * pFn->nExp);
         for (unsigned iExp = 0; iExp < pFn->nExp; ++ iExp) {
            get_next_line(line, in, LINE_Trim | LINE_KeepEmpty);
            replace_fortran_exponents(line);
            std::stringstream str3(line);
            str3 >> pFn->Exp[iExp] >> pFn->Co[iExp];
         }

         // orca exports contraction coefficients with respect to raw primitive Gaussians,
         // not normalized primitive Gaussians. So we need to fix up the coefficients before converting them back.
         // WARNING: at this moment we do not yet know if or if not the functions are spherical!
         if (m_MoldenSource == MOLDENSOURCE_Orca || m_MoldenSource == MOLDENSOURCE_GrimmeXtb) {
            for (unsigned iExp = 0; iExp < pFn->nExp; ++ iExp)
               pFn->Co[iExp] *= OrcaRenormFactor(pFn->Exp[iExp], int(pFn->iAngMom));
         }

         FMoldenAtomSet::iterator
            itAt = Atoms.find(iCenter);
         if (itAt == Atoms.end()) {
            fmt::MemoryWriter w;
            for (FMoldenAtomSet::iterator itAt_ = Atoms.begin(); itAt_ != Atoms.end(); ++ itAt_) {
               if (itAt_ != Atoms.begin())
                  w << " ";
               w << itAt_->first;
            }
//             throw FMoldenParseError(fmt::format("Encountered basis function referring to non-declared center index '{}'. Encountered {} atom indices:\n{}", iCenter, int(Atoms.size()), w.str()));
//             throw FMoldenParseError(fmt::format("Encountered basis function referring to non-declared center index '{}'. Encountered {} atom indices.", iCenter, int(Atoms.size())));
            throw FMoldenParseError(fmt::format("Encountered basis function referring to non-declared center index '{}'. Molden Source: {}, Encountered atom indices:\n{}", iCenter, int(m_MoldenSource), w.str()));
//             FIXME: with pecimens/xtb-molden.input... it hangs if I count the number of atoms encountered. Memory corruption?
         }
         pFn->iCenter = iCenter;
//          pFn->Spherical = m_InputOrbitalsAreSpherical;
//          pFn->iOffsetIn = iInputBasisOffset;
         iInputBasisOffset += pFn->nFn();
         Basis.Shells.push_back(*pFn);
//          Basis.Shells.push_back(FGaussShell(pFn, iCenter, itAt->second.vPos));
      }
   }

 //  Basis.Sort();
 // ^- moved to extra function since data on whether basis is spherical or not
 //    was may or may not have already been given at this point.
}

void FMoldenFile::FixBasis()
{
   if (Basis.Sorted())
      throw FMoldenParseError("called FixBasis more than once. Something went wrong.");

   // make input order and check for spherical vs cartesian functions...
   size_t
      iInputBasisOffset = 0;

   for (size_t iSh = 0; iSh < Basis.Shells.size(); ++ iSh) {
      FMoldenGaussShell
         *pFn = &Basis.Shells[iSh];
      {
         if (pFn->iAngMom >= 5 && !m_AllowLargerThanG)
            throw FMoldenParseError("Molden itself does not support larger-than-g functions. You can ask this program to load them anyway by adding [allow_larger_g] as section to the molden input file.");
         if (pFn->iAngMom == 4 && m_InputOrbitalsAreSpherical && !m_InputOrbitalsAreSphericalForG)
            throw FMoldenParseError("Mixed spherical and cartesian functions are not supported. (g was specified as default cartesian, while df were defined as spherical).");
      }

      pFn->Spherical = m_InputOrbitalsAreSpherical;
      pFn->iOffsetIn = iInputBasisOffset;
      iInputBasisOffset += pFn->nFn();
   }

   Basis.Sort();
}


bool FMoldenGaussShell::operator < (FMoldenGaussShell const &other) const
{
   if (iCenter < other.iCenter) return true;
   if (other.iCenter < iCenter) return false;
   if (iAngMom < other.iAngMom) return true;
   if (other.iAngMom < iAngMom) return false;
   return false;
}


void FMoldenBasis::Sort()
{
   if (Sorted())
      throw FMoldenParseError("called FMoldenBasis::Sort more than once. Something went wrong.");

   size_t
      nAo = nFn(); // may be spherical or cartesian or mixed at this point.
   m_InputToOutputOrder.clear();
   m_InputToOutputOrder.resize(nAo);

   // sort basis shells such that they are ordered as the atoms are.
   std::stable_sort(Shells.begin(), Shells.end());

   size_t
      iOffsetOut = 0;
   // assign new positions to the basis functions.
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh) {
      FMoldenGaussShell
         &Sh = Shells[iSh];
      size_t
         nShFn = Sh.nFn();
      for (size_t iShFn = 0; iShFn < nShFn; ++ iShFn)
         m_InputToOutputOrder[Sh.iOffsetIn + iShFn] = iOffsetOut + iShFn;
      iOffsetOut += nShFn;
   }
   if (0) {
      std::cout << "Basis function order:" << std::endl;
      for (size_t i = 0; i < nAo; ++ i) {
         std::cout << fmt::format("  In[{:4}] -> internal [{:4}]", i, m_InputToOutputOrder[i]) << std::endl;
      }
   }
   m_IsSorted = true;
}

static size_t DoubleFactR(ptrdiff_t l) {
   if (l <= 1)
      return 1;
   else
      return l * DoubleFactR(l - 2);
}


void FMoldenFile::ReadMoSection(std::istream &in)
{
   if (!Basis.Sorted())
      throw FMoldenParseError("Basis not sorted on entry of ReadMoSection.");

   TArray<size_t>
      ls;
   ls.reserve(Basis.Shells.size());
   // make a list of the angular momenta in the basis (need that for converting MO coeffs
   // between cartesian orders and cartesians/sphericals).
   for (size_t iSh = 0; iSh < Basis.Shells.size(); ++ iSh)
      ls.push_back(size_t(Basis.Shells[iSh].iAngMom));
   if (0) {
      std::cout << "Output l order:" << std::endl;
      for (size_t iSh = 0; iSh < Basis.Shells.size(); ++ iSh) {
         std::cout << fmt::format("  Sh[{:3}] -> l = {}", iSh, Basis.Shells[iSh].iAngMom) << std::endl;
      }
   }

   // these ones are required for program-specific hacks (note: both this
   // and ls are in output order! but functions are NOT yet converted to
   // spherical (unless they already were spherical in the input basis)).
   TArray<FMoldenGaussShell*>
      pShellForBf;
   pShellForBf.reserve(Basis.nFn());
   for (size_t iSh = 0; iSh < Basis.Shells.size(); ++ iSh)
      for (size_t iFn = 0; iFn < Basis.Shells[iSh].nFn(); ++ iFn)
         pShellForBf.push_back(&Basis.Shells[iSh]);
   assert(pShellForBf.size() == Basis.nFn());

   // On the input format:
   // the data format is rather... confusing. There is no discernable separator
   // between the data for the different MOs, and orbital data and variables
   // occur intermixed. Additionally, data looks sparse: there is no
   // strict reasons why there actually *IS* data for a MO.
   // Molpro puts a variable section before each individual
   // MO, and the variables can be detected by containing an assignment sign or
   // having a space as first character. I do not know if this is general,
   // however. Here I assume that:
   //    - there is one variable section before each MO
   //    - each orbital contains at least one data line (otherwise it can't
   //      be normalized)
   FScalarData
      // currently read orbital -- still in cartesian coordinates!
      Row;
   size_t
      nFnCart = Basis.nFnCart(),
      nFnSph = Basis.nFnSph();
   if (m_InputOrbitalsAreSpherical)
      Row.resize(nFnSph);
   else
      Row.resize(nFnCart);
//    FMoldenOrbInfo
//       CurrentInfo;
   std::string
      line;
   for ( ; ; ) {
      std::istream &in_file = in;
      std::streampos pos = in_file.tellg();
      if (!get_next_line(line, in_file, LINE_Trim | LINE_AcceptEof))
         break;
      if (IsSectionStart(line)) {
         // rewind to start of line.
         in_file.seekg(pos);
         break;
      }
      in_file.seekg(pos);
//       std::stringstream in(line);
//       std::cout << fmt::format(" ...orb-line: '{}'", line) << std::endl;

      FOrbitalInfoPtr
         pInfo = new FOrbitalInfo;

      // clear data row for the current MO.
      Row.clear_data();

      if ( !in.good() )
         throw FMoldenParseError("Expected variables, but found '?'.");

      // read variables
      std::string
         Key, Val;
      for ( ; ; ) {
         // read variable lines
         std::streampos
            iLast = in.tellg();
//          in >> Key >> Val;
         if (!get_next_line(line, in_file, LINE_Trim | LINE_AcceptEof))
            break;
         size_t
            iEquals = line.find('=');
//          std::cout << fmt::format("     ...k: '{}'  v: '{}'  ok? {}", Key, Val, iEquals) << std::endl;
         if ( !in.good() || iEquals == npos ) {
            // that's not a variable line.
            in.seekg(iLast, std::ios::beg);
            break;
         }
         Key = line.substr(0, iEquals);
         str_trim(Key);
         Val = line.substr(iEquals+1);
         str_trim(Val);
//          std::cout << fmt::format("     ...k: '{}'  v: '{}'", Key, Val) << std::endl;
//          if ( CurrentInfo.Map.find(Key) != CurrentInfo.Map.end() )
//             throw FMoldenParseError(fmt::format("Variable '{}' already assigned while parsing MO {}.", Key, OrbRows.size()));
//          CurrentInfo.Map[Key] = Val;
         str_lower_and_trim(Key);
         if (Key == "sym")
            pInfo->sDesc = Val;
         else if (Key == "ene")
            pInfo->fEnergy = read_fortran_float(Val);
         else if (Key == "spin") {
            pInfo->Spin = ORBSPIN_Unknown;
            str_lower_and_trim(Val);
            if (Val == "alpha")
               // this is also used for spin free orbitals!! Molden does not
               // distinguish them (which means that we can't read files with
               // all orbitals open-shell alpha, since the file format does not
               // allow for separating them from all-orbitals closed-shell)
               pInfo->Spin = ORBSPIN_Alpha;
//                pInfo->Spin = ORBSPIN_Alpha; // guess from occupation.
            if (Val == "beta")
               pInfo->Spin = ORBSPIN_Beta;
         } else if (Key == "occup") {
            pInfo->fOcc = read_fortran_float(Val);
         } else {
            // not recognized. Should we warn?
//             std::cout << fmt::format("     ...unrecognized: '{}' = '{}' ", Key, Val) << std::endl;
         }
      }

      if ( !in.good() )
         throw FMoldenParseError("Expected data, but found '?' while parsing molden file.");

      for ( ; ; ) {
         // read data lines
         std::streampos
            iLast = in.tellg();
         size_t
            iBf;
         double
            fData;
         std::string
            sData;
         in >> iBf >> sData;

         iBf -= 1; // translate 1-based indexing to 0-based
         fData = read_fortran_float(sData);
         if ( !in.good() ) {
            // that's not a data line.
//             in.clear(std::ios::failbit);
            // ^- interesting fact: stream::clear(bit) does not *clear* the given bit,
            //    but rather *sets* the stream state to the given bit. I guess, given the
            //    rest of the C++ stream API, this was to be expected. I guess the real mistake
            //    was trying to use it.
            in.clear();
            in.seekg(iLast, std::ios::beg);
//             bool is_bad = in.bad();
//             bool is_fail = in.fail();
//             bool is_eol = in.eof();
//             bool is_good = in.good();
//             std::string bla;
//             in >> bla;
            break;
         }
//          std::cout << fmt::format("     ...iAo: '{}'  sData: '{}'  fData: {:12.6f}", iBf, sData, fData) << std::endl;
         if ( iBf >= Row.size() )
            throw FMoldenParseError(fmt::format("AO index '{}' too large while parsing MO {}.", (1+iBf), OrbInfo.size()));
         iBf = Basis.m_InputToOutputOrder[iBf];
         if ( Row[iBf] != 0 )
            throw FMoldenParseError(fmt::format("AO component Orb[{},{}] already assigned before.", (1+iBf), OrbInfo.size()));

         if (m_MoldenSource == MOLDENSOURCE_Turbomole || m_MoldenSource == MOLDENSOURCE_GrimmeXtb) {
            size_t l = pShellForBf[iBf]->iAngMom;
            fData *= std::sqrt(DoubleFactR(2*l-1));
         }

         Row[iBf] = fData;
      }
//             std::cout << "Reading file. Input spherical? " << m_InputOrbitalsAreSpherical << " Input from Orca?" << (m_MoldenSource == MOLDENSOURCE_Orca) << std::endl;
      // re-order/convert and store the row
      if (m_InputOrbitalsAreSpherical) {
         assert(Row.size() == nFnSph);

         if (m_MoldenSource == MOLDENSOURCE_Orca) {
            // orca seems to use a definition of solid harmonics which differs by some -1 factors
            // for higher m. I am not sure if this is quite right...
            double *pOrb = &Row[0];
            for (size_t i = 0; i < ls.size(); ++ i) {
               int l = ls[i];
               assert_rt(l == pShellForBf[pOrb-&Row[0]]->iAngMom);
               if (l == 3) { // f
                  for (size_t ic = 5; ic < 7; ++ ic)
                     pOrb[ic] *= -1.;
               }
               if (l == 4) { // g        a b b c c d d e e
                  bool switch_this[9] = {0,0,0,0,0,1,1,1,1};
                  for (size_t ic = 0; ic < 9; ++ ic)
                     if (switch_this[ic])
                        pOrb[ic] *= -1.;
//                   std::swap(pOrb[8],pOrb[7]);
//                   std::swap(pOrb[6],pOrb[5]);
//                   std::swap(pOrb[4],pOrb[3]);
//                   std::swap(pOrb[2],pOrb[1]);
                  // doesn't work :(.
                  //   ...yes, apart from the phase thing, there is also a mystery sqrt(3) factor in the g functions (see OrcaRenormFactor())
                  //   It is running now. I don't even want to know.
               }
               pOrb += size_t(2*l + 1);

               if (l >= 5)
                  throw FMoldenParseError("Sorry, currently only up to 'g' functions are supported for import from ORCA. This input has a higher harmonic.");
            }
         }

         ctimp::Vec_ShMolden2ShMolpro(&Row[0], &ls[0], ls.size());
         pInfo->Orb = Row;
      } else {
         assert(Row.size() == nFnCart);
         ctimp::Vec_CaMolden2CaMolpro(&Row[0], &ls[0], ls.size());
         if (m_MoldenSource == MOLDENSOURCE_Orca || m_MoldenSource == MOLDENSOURCE_Molcas)
            // unlike Molpro, which always converts to Cartesian for export, these programs export spherical functions.
            // So if we encounter a Cartesian .molden file, it most likely means that the actual computational basis
            // was Cartesian.
            throw FMoldenParseError("Sorry, this program does not support Cartesian basis functions. Please use a spherical harmonic basis (e.g., def2-* or cc-pVnZ bases)");

//          if (m_MoldenSource == MOLDENSOURCE_Orca) {
//             std::cout << "RENORMALIZING STUFF!!" << std::endl;
//             ctimp::Vec_Ca2Ca_XSqrtNorm(&Row[0], &ls[0], ls.size());
//          }
         pInfo->Orb.resize(nFnSph);
//          ctimp::Vec_Ca2Ca_XSqrtNorm(&Row[0], &ls[0], ls.size());
         ctimp::Vec_Ca2Sh(&pInfo->Orb[0], &Row[0], &ls[0], ls.size());
      }

      if (m_pLoadOptions && m_pLoadOptions->SkipVirtuals() && pInfo->fOcc < 0.00001)
         continue;
      OrbInfo.push_back(pInfo);
   }
//    std::cout << fmt::format(" ...read {} MOs from file.", OrbInfo.size()) << std::endl;

   FixOrbitalTypes();
}

void FMoldenFile::FixOrbitalTypes()
{
   FixOrbitalTypes1(OrbInfo, ORBSOURCE_Molden);
}

// void FMoldenFile::FixOrbitalTypes()
// {
//    size_t
//       nAlpha = 0,
//       nBeta = 0,
//       nClosed = 0,
//       nOther = 0;
//    FOrbInfoList::iterator
//       it;
//    for (it = OrbInfo.begin(); it != OrbInfo.end(); ++ it) {
//       switch((*it)->Spin) {
//          case ORBSPIN_Alpha: nAlpha += 1; continue;
//          case ORBSPIN_Beta: nBeta += 1; continue;
//          case ORBSPIN_SpinFree: nClosed += 1; continue;
//          default: nOther += 1; continue;
//       }
//    }
//
//    // there are *only* Alpha orbitals. Most likely we got something R(o)HFy.
//    // Adjust spins based on occupation numbers.
//    if (nBeta == 0 && nOther == 0 && nClosed == 0) {
//       for (it = OrbInfo.begin(); it != OrbInfo.end(); ++ it) {
//          double fOcc = (*it)->fOcc;
//          if (fOcc == 2. || fOcc == 0.)
//             (*it)->Spin = ORBSPIN_SpinFree;
//       }
//    }
// }



// struct FMoldenGaussShell
// {
//    ct::TArray<double>
//       Co,  // nExp x nCo contraction matrix
//       Exp; // nExp primitive exponents
//    size_t
//       nExp,
//       nCo;
//    int
//       iAngMom,
//       iCenter;
// //    FMoldenGaussShell() {};
// };


template <class FSequence>
int CmpArrays(FSequence const &A, FSequence const &B)
{
   size_t
      nA = A.size(),
      nB = B.size();
   if (nA < nB) return -1;
   if (nB < nA) return +1;
   typename FSequence::const_iterator
      itA = A.begin(),
      itB = B.begin();
   for ( ; itA != A.end(); ) {
      if (*itA < *itB) return -1;
      if (*itB < *itA) return +1;
      ++ itA;
      ++ itB;
   }
   return 0;
}

// compare two IR gauss shell objects, in order to find equivalent ones on different atoms.
struct FGaussFnCmp
{
   bool operator () (ct::FAtomShellPtr const &pA, ct::FAtomShellPtr const &pB) const {
      if (pA->AngMom < pB->AngMom) return true;
      if (pB->AngMom < pA->AngMom) return false;
      int iCmp;
      iCmp = CmpArrays(pA->Exponents, pB->Exponents);
      if (iCmp < 0) return true;
      if (iCmp > 0) return false;
      iCmp = CmpArrays(pA->CoMatrix, pB->CoMatrix);
      if (iCmp < 0) return true;
      if (iCmp > 0) return false;
      return false;
   }
};


ct::FAtomSetPtr FMoldenAtomSet::Convert()
{
   ct::FAtomSetPtr
      pAtoms = new ct::FAtomSet();
   // it's ordered in iCenterOut order.
   for (FMoldenAtomSet::iterator itAtom = this->begin(); itAtom != this->end(); ++ itAtom) {
      FMoldenAtom
         &At = itAtom->second;
      // note: technically we also have element names from input... (both labels and numbers).
      // I guess the numbers might be safer? Because programs might use strange pre/suffices or
      // spellings of elements...
//       pAtoms->Atoms.push_back(ct::FAtom(At.vPos, ct::ElementNameFromNumber(At.iElement), "(external)"));
      pAtoms->AddAtom(ct::FAtom(At.vPos, ct::ElementNameFromNumber(At.iElement), "(external)"));
   }
   return pAtoms;
}

ct::FBasisSetPtr FMoldenBasis::Convert(FMoldenAtomSet *pMoldenAtomSet)
{
   typedef std::set<ct::FAtomShellPtr, FGaussFnCmp>
      FAtomShellSet;
   FAtomShellSet
      AtomShells;
   std::vector<ct::FBasisShell>
      NewShells;
   NewShells.reserve(Shells.size());

   // Go through the molden shells. At this point they should already have been ordered
   // by center and angular momentum, so we just need to join them and fix the output
   // atom indices.
   size_t
      iShBeg = 0, iShEnd;
   for ( ; iShBeg != Shells.size(); iShBeg = iShEnd) {
      // find all shells on the same center with the same AngMom
      iShEnd = iShBeg + 1;
      while (iShEnd != Shells.size()
          && Shells[iShEnd].iCenter == Shells[iShBeg].iCenter
          && Shells[iShEnd].iAngMom == Shells[iShBeg].iAngMom)
         iShEnd += 1;
      // make a sorted list of all unique exponents. We here assume that exponents come
      // out at the same binary version of IEEE doubles if they are equivalent (should be fine)
      typedef std::map<double, size_t>
         FExpMap;
      FExpMap
         // maps negative of exponent to exponent index.
         // (negative such that we get the large ones first)
         ExpMap;
      for (size_t iSh = iShBeg; iSh != iShEnd; ++ iSh)
         for (size_t iExp = 0; iExp < Shells[iSh].Exp.size(); ++ iExp)
            ExpMap[-Shells[iSh].Exp[iExp]] = 0;
      // assign exponent indices and make a single ordered list of exponents
      TArray<double>
         Exps;
      size_t
         nExp = ExpMap.size();
      Exps.reserve(nExp);
      for (FExpMap::iterator it = ExpMap.begin(); it != ExpMap.end(); ++ it) {
         it->second = Exps.size();
         Exps.push_back(-it->first);
      }
      // make a constraction matrix.
      size_t
         nCo = iShEnd - iShBeg;
      TArray<double>
         Cos;
      Cos.resize(nExp * nCo);
      Cos.clear_data();
      for (size_t iSh = iShBeg; iSh != iShEnd; ++ iSh) {
         size_t
            iCo = iSh - iShBeg;
         FMoldenGaussShell
            &Sh = Shells[iSh];
         for (size_t iExp = 0; iExp < Sh.Exp.size(); ++ iExp) {
            // find the index of the current exponent in the output exponent list.
            FExpMap::const_iterator
               itExp = ExpMap.find(-Sh.Exp[iExp]);
            if (itExp == ExpMap.end())
               throw FMoldenParseError("Internal indexing error in re-assembly of imported basis shells.");
            // and add the current coefficient to it (technically an exponent might occur multiple
            // times, this is why we add instead of replacing)
            Cos[nExp * iCo + itExp->second] += Sh.Co[iExp];
         }
      }

      // make a CT/IR shell object
      ct::FAtomShellPtr
         pIrSh = new ct::FAtomShell(unsigned(Shells[iShBeg].iAngMom), &Exps[0], unsigned(Exps.size()), &Cos[0], unsigned(nCo), 0);
      // check if we had one of those things already before.
      std::pair<FAtomShellSet::iterator,bool>
         itInsResult = AtomShells.insert(pIrSh);
      if (itInsResult.second)
         // take the previous one.
         pIrSh = *itInsResult.first;

      // now make a single BasisShell object referring to the generally contracted atom shell.
      FMoldenAtom
         &At = pMoldenAtomSet->Get(Shells[iShBeg].iCenter);
      NewShells.push_back(ct::FBasisShell(At.vPos, At.iCenterOut, pIrSh));

   }
   return ct::FBasisSetPtr(new ct::FBasisSet(&NewShells[0], NewShells.size(), ct::BASIS_Orbital, "(external)"));
}


// return number of electrons treated as core (and therefore not occurring explicitly
// in the wave function) in the GFN2-XTB method of Grimme and coworkers.
// These are needed to compute partial charges etc. correctly.
unsigned GetNumCoreElec_Gfn2Xtb(unsigned iElement)
{
   // WARNING: this isn't quite right for all the elements.
   // It's just to get things started.
   return ct::GetElementNumCoreElec(iElement);
}


FMolproXmlDataPtr LoadMoldenFile(std::string const &FileName, FLoadOptions const &LoadOptions)
{
   FMolproXmlDataPtr
      pXmlData = new FMolproXmlData();
   std::ifstream
      in(FileName.c_str(), std::ifstream::in | std::ifstream::binary);
   // ^- @binary: the files are TEXT files, but otherwise we get line ending hazard,
   //    such that in windows we cannot read lines which are terminated by lf only, not cr lf.
   //    There may be better hacks around this (see
   //    http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf),
   //    but so far I have not been overly impressed with the suggestions. So hack for now: read as
   //    binary and discard any \r's occuring anywhere.
   //    I guess one of the savest things to do would be to load the entire file, check if there
   //    are any \n's at all, if yes: replace all \r's by nothing, if not: replace all \r's by \n's.
   //    But that is also not intensely pretty.
   if (!in.good())
      throw FMoldenParseError("File load failed.");
   FMoldenFile
      MoldenData(in, &LoadOptions);

   // convert to output format.
   pXmlData->pAtoms = MoldenData.Atoms.Convert();
   pXmlData->pBasisOrb = MoldenData.Basis.Convert(&MoldenData.Atoms);
   if (MoldenData.m_MoldenSource == MOLDENSOURCE_GrimmeXtb) {
      // this one won't be quite correct if we'd use the ab-initio MINAO minimal basis
      // for IAO construction or analysis. To analyze such files, we should use
      // the minimal basis subset of the AO basis it employs itself.
      // We'll mark this by replacing the assigned basis name for the BASIS_Guess
      // context for these atoms.
      // Additionally, we assign GFN2-XTB numbers as ECP charge to account
      // for the omitted core electrons in the valence-only treatment.
      for (size_t iAt = 0; iAt != pXmlData->pAtoms->size(); ++ iAt) {
         ct::FAtom
            &At = (*pXmlData->pAtoms)[iAt];
//          At.BasisDesc[ct::BASIS_Guess] = "GFN2-XTB-MIN";
         At.BasisDesc[ct::BASIS_Guess] = "USE-ORB-BASIS";
         At.EcpDesc = "GFN2-XTB-CORE";
         // std::cout << fmt::format("!set atom {} guess basis to {}.\n", iAt+1, At.BasisDesc[ct::BASIS_Guess]);
         At.nEcpElec = GetNumCoreElec_Gfn2Xtb(At.iElement);
      }
   }

   pXmlData->pOrbSet = new FOrbitalSet();
   FOrbitalSet::FOrbitalInfoList
      &OrbInfos = pXmlData->pOrbSet->OrbInfos;
   OrbInfos.reserve(MoldenData.OrbInfo.size());
   bool AddOrbitalNumbers = false;
   for (size_t iOrb = 0; iOrb < MoldenData.OrbInfo.size(); ++ iOrb) {
      OrbInfos.push_back(MoldenData.OrbInfo[iOrb]);
      FOrbitalInfoPtr &oi = OrbInfos.back();
      oi->pBasisSet = pXmlData->pBasisOrb;

      // check if the orbital description contains some sort of string which
      // can be interpreted as the target orbital number in IboView counting format.
      std::string
         sOrbNumber = fmt::format("{}", iOrb+1);
      if (oi->sDesc.find(sOrbNumber) == std::string::npos)
         AddOrbitalNumbers = true;
//       oi->sDesc = fmt::format("{:4}.{}",(1+iOrb), oi.iSym);

   }
   if (AddOrbitalNumbers) {
      // at least one of the orbital description strings does not contain the orbital ID.
      // so add orbital IDs to all of them (in some cases SOME can contain valid IDs, but
      // not all; e.g., when orbitals are numbered per symmetry. Then 1.1, 2.1, 3.1 etc.
      // would correspond to IboView orbital IDs 1,2,3, but 1.2, 2.2, 3.2, etc would not)
      for (size_t iOrb = 0; iOrb < OrbInfos.size(); ++ iOrb) {
         FOrbitalInfoPtr &oi = OrbInfos[iOrb];
         str_trim(oi->sDesc);
         oi->sDesc = fmt::format("{:3}. {}", iOrb+1, oi->sDesc);
      }
   }
   return pXmlData;
}


} // namespace molden_file




namespace turbo_file {

typedef std::vector<std::string>
   FStringList;

enum {
   TURBOREAD_TryReadHeaderLine = 0x01,
   TURBOREAD_TryReadGradient = 0x02
};


// this reads a SINGLE geometry from the current $coord or $grad section.
FAtomSetPtr ReadGeometryFromCoordsSection(std::istream &in, uint Flags, FLoadOptions const &LoadOptions)
{
   // This involves some guesswork, as we here do not know how many atoms there are supposed to be.
   // This works as follows:
   //   - if TryReadHeaderLine is set, attempt to read energy and cycle number (from 'gradient' files),
   //     and then continue with next line.
   //   - read xyz coords of atoms as long as the format is 'x y z elementname [optinal f/k/etc.]'
   //     stop reading if we either get a 'x y z' (without element name---this would then be
   //     interpreted as a gradient) or a section id (something starting with '$') or
   //     an EOF.
   //   - rewind to first line we have not consumed.
   FAtomSetPtr
      pAtoms( new ct::FAtomSet() );
   std::string
      FrameDesc;
   double
      FrameEnergy = 0.,
      FrameGrad = 0.;
   long
      FrameCycle = 0;
   FStringList
      ls;
   if (Flags & TURBOREAD_TryReadHeaderLine) {
      std::streampos
         pos = in.tellg();
      std::string
         Header;
      if (!get_next_line(Header, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof | LINE_SkipTurboCommentLines))
         return 0;
//        cycle =      1    SCF energy =    -1018.6937684340   |dE/dxyz| =  0.019278
//       boost::split(ls, Header, boost::is_any_of(" \t"), boost::token_compress_on);
      if (Header == "$end")
         return 0;
      boost::split(ls, Header, boost::is_any_of("= \t"), boost::token_compress_on);
      if (!ls.empty())
         str_trim_right(ls.front());
      if (!ls.empty() && ls.front() == "cycle") {
         FStringList::iterator
            it = ls.end();
         if (it != ls.begin()) --it; // 0.019278
         if (it != ls.end())
            try_read_fortran_float(FrameGrad, it->c_str());
         if (it != ls.begin()) --it; // |dE/dxyz|
         if (it != ls.begin()) --it; // -1018.6937684340
         if (it != ls.begin() && it != ls.end())
         if (it != ls.end())
            try_read_fortran_float(FrameEnergy, it->c_str());
         it = ls.begin();
         if (it != ls.end())
            ++ it;
         if (it != ls.end())
            FrameCycle = std::strtol(it->c_str(), 0, 10);
//          std::cout << fmt::format("reading $grad section: cycle: {}  energy: {:16.8f}", FrameCycle, FrameEnergy) << std::endl;

         pAtoms->SetLastEnergy(FrameEnergy);
         pAtoms->SetName(fmt::format("cycle: {}; energy: {:.8f}; grad: {:.2e}", FrameCycle, FrameEnergy, FrameGrad));
      } else {
         // not a header line?
         in.seekg(pos);
      }
   }

   // add atoms positions.
   for ( ; ; ) {
      std::streampos
         pos = in.tellg();
      std::string
         Line;
      if (!get_next_line(Line, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof | LINE_SkipTurboCommentLines))
         break;
      boost::split(ls, Line, boost::is_any_of(" \t"), boost::token_compress_on);
      if (!(ls.size() == 4 || ls.size() == 5) || (!ls.empty() && starts_with(ls[0], "$"))) {
         // not an atom line anymore. Either next section, or something else.
         in.seekg(pos);
         break;
      }

      ct::FVector3 vPos;
      std::string ElementName;
      vPos[0] = read_fortran_float(ls[0]);
      vPos[1] = read_fortran_float(ls[1]);
      vPos[2] = read_fortran_float(ls[2]);
      // ^- I think these come in Bohr units. No conversion needed then.
      ElementName = ls[3];
      // if ls.size() == 5, ls[4] would contain additional turbomole coord specs
      // like 'f' (fixed). I think we do not actually need those here.
      pAtoms->AddAtom(ct::FAtom(vPos, ElementName, ""));
   }

   // add atom gradients, if provided.
   if (Flags & TURBOREAD_TryReadGradient) {
      size_t
         iCurrentAtom = 0;
      for ( ; ; ) {
         std::streampos
            pos = in.tellg();
         std::string
            Line;
         if (!get_next_line(Line, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof | LINE_SkipTurboCommentLines))
            break;
         boost::split(ls, Line, boost::is_any_of(" \t"), boost::token_compress_on);
         if (ls.size() != 3 || (!ls.empty() && starts_with(ls[0], "$"))) {
            // not an gradient line anymore. Either next section, or something else.
            in.seekg(pos);
            break;
         }

         ct::FVector3 vGrad;
         vGrad[0] = read_fortran_float(ls[0]);
         vGrad[1] = read_fortran_float(ls[1]);
         vGrad[2] = read_fortran_float(ls[2]);
         // FIXME: convert Turbomole gradients to something else? They came in a funky unit iirc.
         if (iCurrentAtom >= pAtoms->size())
            throw FFileLoadError(fmt::format("While reading gradient information: Too much data -- more gradient ({}) than atoms ({}) while reading '{}'.'", iCurrentAtom+1, pAtoms->size(), Line));
         (*pAtoms)[iCurrentAtom].vGrad = vGrad;
         iCurrentAtom += 1;
      }
   }

   return pAtoms;
   IR_SUPPRESS_UNUSED_WARNING(LoadOptions);
}


FMolproXmlDataPtrList LoadTurboFile(std::string const &FileName, FLoadOptions const &LoadOptions)
{
   FMolproXmlDataPtrList
      XmlFrames;
   std::ifstream
      in(FileName.c_str());
   std::string
      Line;
   for ( ; ; ) {
      if (!get_next_line(Line, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof))
         break;
      FStringList
         ls;
      boost::split(ls, Line, boost::is_any_of(" \t"), boost::token_compress_on);

//       std::cout << fmt::format("ls.size = {:4}  ls.front = '{}'", ls.size(), ls.front()) << std::endl;
      if (!ls.empty() && ls.front() == "$coord") {
//          std::cout << "reading $coord section" << std::endl;
         FMolproXmlDataPtr
            pXmlData(new FMolproXmlData);
         pXmlData->pAtoms = ReadGeometryFromCoordsSection(in, 0, LoadOptions);
         XmlFrames.push_back(pXmlData);
         return XmlFrames;
      } else if (!ls.empty() && ls.front() == "$grad") {
//          std::cout << "reading $grad section" << std::endl;
         for ( ; ; ) {
            FMolproXmlDataPtr
               pXmlData(new FMolproXmlData);
            pXmlData->pAtoms = ReadGeometryFromCoordsSection(in, TURBOREAD_TryReadGradient | TURBOREAD_TryReadHeaderLine, LoadOptions);
            if (pXmlData->pAtoms.get() == 0 || pXmlData->pAtoms->size() == 0)
               break;
            XmlFrames.push_back(pXmlData);
         }
         return XmlFrames;
      }
   }
   return XmlFrames;

   IR_SUPPRESS_UNUSED_WARNING(LoadOptions);
}





} // namespace turbo_file








// int main(int argc, char *argv[])
// {
//    using namespace molpro_xml;
//    FMolproXmlDataPtr pXmlData = LoadMolproXmlFile("test1.xml");
//    TestOrbitals(pXmlData->pBasisOrb, pXmlData->pAtoms, pXmlData->pOrbSet);
// };

// int main(int argc, char *argv[])
// {
//    using namespace molpro_xml;
//    pugi::xml_document doc;
//    pugi::xml_parse_result result = doc.load_file("test1.xml");
//
//    std::cout << "XML Load result: " << result.description() << std::endl;
//
//
//    pugi::xml_node MoleculeNode = doc.child("molpro").child("molecule");
// //    EnumerateChildren(std::cout, MoleculeNode);
//    pugi::xml_node GeometryNode = MoleculeNode.child("cml:atomArray");
// //    EnumerateChildren(std::cout, GeometryNode);
//    FAtomSetPtr
//       pAtoms = LoadGeometry(GeometryNode);
// //    std::cout << *pAtoms;
//
//    pugi::xml_node OrbBasisNode = MoleculeNode.find_child_by_attribute("basisSet", "id", "ORBITAL");
// //    EnumerateChildren(std::cout, OrbBasisNode);
//
//    ct::FBasisSetPtr
//       pBasisOrb = LoadBasisSet(OrbBasisNode, *pAtoms);
// //    std::cout << *pBasisOrb;
//
//    pugi::xml_node
//       OrbitalsNode = MoleculeNode.child("orbitals");
// //    EnumerateChildren(std::cout, OrbitalsNode);
//    FOrbitalSetPtr
//       pOrbSet = LoadOrbitals(OrbitalsNode, pBasisOrb);
//    // do I still need to normalize the functions?
//
// //    TestOrbitals(pBasisOrb, pAtoms, pOrbSet);
//
//
//
//
// //    std::cout << doc.child("molpro").child("molecule").child("cml:symmetry").value() << std::endl;
// //    for (pugi::xml_node tool = doc.child("Tool"); tool; tool = tool.next_sibling("Tool"))
// //    {
// //       std::cout << "Tool " << tool.attribute("Filename").value();
// //       std::cout << ": AllowRemote " << tool.attribute("AllowRemote").as_bool();
// //       std::cout << ", Timeout " << tool.attribute("Timeout").as_int();
// //       std::cout << ", Description '" << tool.child_value("Description") << "'\n";
// //    }
//
//
// //    std::cout << doc.child("mesh").child("molecule").child("cml:symmetry").value() << std::endl;
// //    for (pugi::xml_node tool = tools.child("Tool"); tool; tool = tool.next_sibling("Tool"))
// //    {
// //       std::cout << "Tool " << tool.attribute("Filename").value();
// //       std::cout << ": AllowRemote " << tool.attribute("AllowRemote").as_bool();
// //       std::cout << ", Timeout " << tool.attribute("Timeout").as_int();
// //       std::cout << ", Description '" << tool.child_value("Description") << "'\n";
// //    }
// //    for (pugi::xml_node tool = tools.child("Tool"); tool; tool = tool.next_sibling("Tool"))
// //    {
// //       std::cout << "Tool " << tool.attribute("Filename").value();
// //       std::cout << ": AllowRemote " << tool.attribute("AllowRemote").as_bool();
// //       std::cout << ", Timeout " << tool.attribute("Timeout").as_int();
// //       std::cout << ", Description '" << tool.child_value("Description") << "'\n";
// //    }
//
// //    << ", mesh name: " << doc.child("mesh").attribute("name").value() << std::endl;
// }


//    FMolproXmlDataPtr LoadOrbitalFile(std::string const &FileName, FLoadOptions const &LoadOptions)
//    {
//       try {
//          if (ends_with(FileName, ".xml"))
//             return molpro_xml::LoadMolproXmlFile(FileName, LoadOptions);
//          std::string
//             FirstLine;
//          {
//             std::ifstream
//                in(FileName.c_str());
//             if (!get_next_line(FirstLine, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof))
//                FirstLine = "!2empty44";
//          }
//          if (ends_with(FileName, ".molden") || ends_with(FileName, ".molden.input") || FirstLine == "[molden format]")
//             return molden_file::LoadMoldenFile(FileName, LoadOptions);
//
//          throw FFileTypeUnrecognizedError("Failed to recognize file extension.");
//       } catch (FFileTypeUnrecognizedError &e) {
//          throw FFileTypeUnrecognizedError(e.what(), FileName);
//       } catch (FFileLoadError &e) {
//          throw FFileLoadError(e.what(), FileName);
//       }
//    }

   FMolproXmlDataPtrList LoadOrbitalFile(std::string const &FileName, FLoadOptions const &LoadOptions)
   {
      FMolproXmlDataPtrList
         Frames;
      try {
         if (ends_with(FileName, ".xml")) {
            // FIXME: make it possible to load multiple frames here.
            Frames.push_back(molpro_xml::LoadMolproXmlFile(FileName, LoadOptions));
            return Frames;
         }
         std::string
            FirstLine;
         {
            std::ifstream
               in(FileName.c_str());
            if (!get_next_line(FirstLine, in, LINE_Trim | LINE_ToLower | LINE_AcceptEof))
               FirstLine = "!2empty44";
         }
         if (ends_with(FileName, ".molden") || ends_with(FileName, ".molden.input") || FirstLine == "[molden format]") {
            // FIXME: Molden format may contain multiple sets of geometries (but only one set of orbtials).
            // Allow loading them here?
            Frames.push_back(molden_file::LoadMoldenFile(FileName, LoadOptions));
            return Frames;
         }
//          if (FirstLine == "$coord" || starts_with(FirstLine, "$grad ") || FirstLine == "$control") {
         if (FirstLine == "$coord" || starts_with(FirstLine, "$grad ")) {
            return turbo_file::LoadTurboFile(FileName, LoadOptions);
         }
         throw FFileTypeUnrecognizedError("Failed to recognize file extension.");
      } catch (FFileTypeUnrecognizedError &e) {
         throw FFileTypeUnrecognizedError(e.what(), FileName);
      } catch (FFileLoadError &e) {
         throw FFileLoadError(e.what(), FileName);
      }
      // should never get here.
      assert(!"in LoadOrbitalFile(): should not be here.");
      return FMolproXmlDataPtrList();
   }


   FMolproXmlData::FMolproXmlData()
   {
      SemiEmpiricalMode = false;
   }

} // namespace orbital_file
