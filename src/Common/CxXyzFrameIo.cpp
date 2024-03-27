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

#include <stdexcept>
#include <cmath>
#include <cctype> // for std::isdigit, std::isalpha, and std::tolower/std::toupper and so on.
#include <set>
#include <sstream>
#include <cstdlib> // for strtol/strtod

#include "CxXyzFrameIo.h"
#include "CxParse1.h" // for string_slice
#include "CxPhysicalUnits.h"

// #define USE_FSTREAM_IN_XYZIO
#ifdef USE_FSTREAM_IN_XYZIO
   #include <fstream>
#else
   #include "CxIo.h" // for LoadFileIntoMemory, tolower (of strings), Trim, and StripLineComment.
   // ^- Not sure if we really want that dependency for LoadFileIntoMemory only...
   //    (not the least since in both instances it is used, it is anyway
   //    attached to a std::stringstream object. We could just attach it to an
   //    ifstream in text mode directly, and apart from error reporting, not
   //    being able to defer to QT if loaded in a QT context, and possibly
   //    efficiency, it would not change much).
   //
   //    The string functions which it uses should probably anyway be replaced
   //    by CxParse1's slice versions... (slow!), but changing the xyz parse
   //    code would itself likely require several hours to get right, due to all
   //    the possible complications and variants of .xyz files. Not ready
   //    to deal with this at the moment.
#endif // USE_FSTREAM_IN_XYZIO


using ct::string_slice;

namespace xyz_io {


double const
   g_XyzDistanceScale = ct::ToAng,
   g_XyzGradientScale = -1. * (ct::ToEv/ct::ToAng); // we store negative gradients, in eV/A.


//
// =============================================================================
//
// UPCOMING: .xyz exception reporting
//
// -----------------------------------------------------------------------------

FXyzLoadException::FXyzLoadException(std::string const &Reason, std::string const &FileName)
   : std::runtime_error(fmt::format("Failed to load '{}' as .xyz: {}",FileName, Reason))
{}

FXyzLoadException::FXyzLoadException(std::string const &Reason, std::string const &CurLine, std::string const &FileName)
   : std::runtime_error(fmt::format("Failed to load '{}' as .xyz: {} (offending line: '{}')",FileName, Reason, CurLine))
{}



//
// =============================================================================
//
// UPCOMING: .xyz file loading options / conversions
//
// -----------------------------------------------------------------------------


FXyzIoOptions::FXyzIoOptions()
   : m_OutputCoordsFactor(g_XyzDistanceScale),
     m_OutputGradientFactor(g_XyzGradientScale)
{
}


void FXyzIoOptions::SetOutputFactors(double fCoordsFactor, double fGradientFactor)
{
   m_OutputCoordsFactor = fCoordsFactor;
   m_OutputGradientFactor = fGradientFactor;
}



FVector3 FXyzIoOptions::ConvertToInternalCoords(FVector3 const &vPos) const
{
   // our internal data structures are actually input in a_bohrs.
   return (1./m_OutputCoordsFactor) * FVector3(vPos[0], vPos[1], vPos[2]);
}


FVector3 FXyzIoOptions::ConvertToInternalGradient(FVector3 const &vGrad) const
{
   return (1./m_OutputGradientFactor) * FVector3(vGrad[0], vGrad[1], vGrad[2]);
}


FVector3 FXyzIoOptions::ConvertToOutputCoords(FVector3 const &vPos) const
{
   return m_OutputCoordsFactor * FVector3(vPos[0], vPos[1], vPos[2]);
}


FVector3 FXyzIoOptions::ConvertToOutputGradient(FVector3 const &vGrad) const
{
   return m_OutputGradientFactor * FVector3(vGrad[0], vGrad[1], vGrad[2]);
}



FXyzLoadOptions::FXyzLoadOptions(unsigned LoadFlags_)
   : m_LoadFlags(LoadFlags_)
{
}


FXyzLoadOptions::FXyzLoadOptions(FXyzIoOptions const &Base, unsigned LoadFlags_)
   : FXyzIoOptions(Base),
     m_LoadFlags(LoadFlags_)
{
}


std::string FXyzLoadOptions::GetNormalizedAtomType(std::string const &AtomType_) const
{
   std::string
      AtomType;
   if (0 != (ATOMTYPE_ToLower & m_LoadFlags)) {
//       AtomType = ct::tolower(AtomType_);
      AtomType = AtomType_;
      for (size_t i = 0; i != AtomType.size(); ++ i)
         AtomType[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(AtomType[i])));
   } else {
      AtomType = AtomType_;
   }

   if (0 != (ATOMTYPE_Normalize & m_LoadFlags)) {
      if (!AtomType.empty())
         AtomType[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(AtomType[0])));
   }


   if (0 != (ATOMTYPE_RemoveTrailingNumbers & m_LoadFlags)) {
      // in some formats elements are marked with numbers in order to distinguish
      // different versions of them (e.g., there might be symbols C1 and C2 for
      // carbons using differnt basis sets). We cannot deal with this here, so we
      // simply strip off any digits which might be in there.
      while (!AtomType.empty() && std::isdigit(AtomType[AtomType.size()-1]))
         AtomType.resize(AtomType.size() - 1);
   }

   return AtomType;
}



FXyzPrintOptions::FXyzPrintOptions(unsigned PrintFlags_)
{
   InitPrintOptions(PrintFlags_);
}


FXyzPrintOptions::FXyzPrintOptions(FXyzIoOptions const &Base, unsigned PrintFlags_)
   : FXyzIoOptions(Base)
{
   InitPrintOptions(PrintFlags_);
}


void FXyzPrintOptions::InitPrintOptions(unsigned PrintFlags_)
{
   m_PrintFlags = PrintFlags_;
   m_OutputLineLead = "";
}




// FXyzPrintOptions::MakeDesc()
// {
//    fmt::MemoryWriter w;
//    w.write("energy = {:.12f} ", Energy);
//    if ((0 != (XyzFlags & XYZFLAGS_IncludeGradientNorm)) && !Grad.empty()) {
//       w.write("grad = {:.2e} ", fGradNorm());
//    }
//    if (!FrameDesc.empty())
//       w << FrameDesc << " ";
//    if (!Desc1.empty())
//       w << Desc1 << " ";
//    std::string Desc = w.str();
//    out << fmt::format("{}{}\n{}{}\n", Lead, size(), Lead, Desc);
//    for (unsigned i = 0; i < size(); ++i) {
//       out << fmt::format("{}{:<4} {:16.10f} {:16.10f} {:16.10f}",
//          Lead,
//          ElementNameFromNumber(Elements[i]),
//          (xyz_io::g_XyzDistanceScale * Coords[i][0]),
//          (xyz_io::g_XyzDistanceScale * Coords[i][1]),
//          (xyz_io::g_XyzDistanceScale * Coords[i][2]));
//       if (!Grad.empty() && (0 != (XyzFlags & XYZFLAGS_IncludeGradients))) {
//          out << fmt::format(" {:8.2f} {:14.8f} {:14.8f} {:14.8f}",
//             0.0, // partial charge
//             (Grad[i][0]*xyz_io::g_XyzGradientScale),
//             (Grad[i][1]*xyz_io::g_XyzGradientScale),
//             (Grad[i][2]*xyz_io::g_XyzGradientScale));
//       }
//       out << "\n";
//    }
//    out.flush();
// }

std::string FXyzPrintOptions::MakeDescLine(std::string const &FrameDesc, std::string const &ExtraDesc, double Energy, double Gradient) const
{
   fmt::MemoryWriter w;
   w.write("energy = {:.12f} ", Energy);
   if (HasPrintFlag(XYZPRINT_IncludeGradientNorm) && Gradient >= 0.) {
      w.write("grad = {:.2e} ", Gradient);
   }
   if (!FrameDesc.empty())
      w << FrameDesc << " ";
   if (!ExtraDesc.empty())
      w << ExtraDesc << " ";
   return w.str();
}


//
// =============================================================================
//
// UPCOMING: auxiliary classes for representing the loaded .xyz data
//
// -----------------------------------------------------------------------------

FIoAtom::FIoAtom()
{
}


FIoAtom::FIoAtom(FVector3 const &vPos_, std::string const &AtomType_)
{
   Init(vPos_, AtomType_);
}


void FIoAtom::Init(FVector3 const &vPos_, std::string const &AtomType_)
{
   this->Type = AtomType_;
   this->vPos = vPos_;
   this->vGrad = FVector3(0., 0., 0.);
}


std::string FIoAtom::Element() const
{
   std::string
      s;
   s.reserve(Type.size());
   // copy over characters of 'Type' from left until we find a non-alphabetic one.
   //
   // This should work for descs like "C1" or "N1e" or "N2e" (nitrogens with one or two
   // pi electrons) for example, and just return the alphabetic part of this.
   // If the atomic subtype also starts alphabetically (as in some force fields) it may
   // be required to put some extra letter in-between (e.g., "C-terminal").
   for (size_t i = 0; i != Type.size(); ++ i) {
      int
         c = int(Type[i]);
      if (!std::isalpha(c))
         // first non-alphabetic character.
         break;
      if (i == 0)
         // upper-case-ify first letter of element symbol
         c = std::toupper(c);
      else
         // lower-case-ify all other letters of element symbol
         c = std::tolower(c);
      s.push_back(char(c));
   }
   return s;
}

// int FIoAtom::iElement() const
// {
//    return ct::ElementNumberFromName(Element());
// }


std::string FIoAtom::Tag(bool ConvertToLowerCase) const
{
   std::string::size_type
      iTagStart = this->FindTagStartIndex();
   if (iTagStart == std::string::npos) {
      // no tag in atom type string
      return std::string();
   } else {
      // extract substr from here to the end of the type.
      std::string
         sTag = Type.substr(iTagStart);
      if (ConvertToLowerCase) {
         for (size_t k = 0; k != sTag.size(); ++ k)
            sTag[k] = char(std::tolower(int(sTag[k])));
      }
      return sTag;
   }
}


double FIoAtom::NumericTag(double DefaultTag) const
{
   std::string::size_type
      iTagStart = this->FindTagStartIndex();
   if (iTagStart == std::string::npos) {
      // no tag in atom type string
      return DefaultTag;
   } else {
      assert(iTagStart < this->Type.size());
      char const
         *pTypeBeg = this->Type.c_str(),
         *pTagBeg = pTypeBeg + iTagStart,
         *pTagEnd = 0;
      double
         r_ = std::strtod(pTagBeg, const_cast<char**>(&pTagEnd));
      if (pTagEnd != pTypeBeg + this->Type.size())
//          throw FXyzLoadException(fmt::format("Atom {} {}: failed to convert atom subtype tag '{}' to an integer while parsing atom type '{}'.", iAt+1, sElement, sTag, Atoms[iAt].Type), m_InputName);
         throw std::runtime_error(fmt::format("FIoAtom::NumericTag(): failed to convert atom subtype tag '{}' to a number while parsing atom type '{}'.", pTagBeg, Type));
      return r_;
   }
}


std::string::size_type FIoAtom::FindTagStartIndex() const
{
   for (size_t i = 0; i != Type.size(); ++ i) {
      int
         c = int(Type[i]);
      if (!std::isalpha(c))
         // found first non-alphabetic character. That's where we assume
         // the tag part (typically a single integer) starts.
         return i;
   }
   // loop terminated without encountering a non-alphabetic character.
   // Assume that means there are no tags present (in case of alphabetic
   // tags they should be indicated with a delimiter like '-' or ":")
   return std::string::npos;
}



//
// =============================================================================
//
// UPCOMING: actual .xyz data loading routines
//
// -----------------------------------------------------------------------------


FIoAtomSet::FIoAtomSet()
{
}


FIoAtomSet::~FIoAtomSet()
{
}


void FIoAtomSet::AddAtom(FIoAtom const &IoAtom)
{
   m_Atoms.push_back(IoAtom);
}


// make a local stream object (either an ifstream or a stringstream) which is
// attached to an input file FileName in text mode. Used for reading .xyz
// files from.
#ifdef USE_FSTREAM_IN_XYZIO
   // make a std::ifstream and return it directly (avoids dependency on
   // ct::LoadFileIntoMemory)
   #define MAKE_ISTREAM_FOR_FILE(str, FileName) \
      std::ifstream \
         str((FileName).c_str()); \
      if (!str.good()) \
         throw FXyzLoadException("Could not read file from disk", (FileName));
#else
   // use ct::LoadFileIntoMemory to load the entire file into a character buffer
   // first and then attach a std::stringstream to it (has dependency, but
   // should be faster and safer regarding IO, and allows substituting the
   // file loading routine by QT's (to access resource and clipboard content)
   // if used in a QT context (done in CxIo.cpp)
   #define MAKE_ISTREAM_FOR_FILE(str, FileName) \
      ct::TArray<char> \
         pFileContent; \
      if (false == ct::LoadFileIntoMemory(pFileContent, (FileName))) \
         throw FXyzLoadException("Could not read file from disk", (FileName)); \
      std::stringstream \
         str(&pFileContent[0], std::stringstream::in);
#endif // USE_FSTREAM_IN_XYZIO


// Adds the molecules out of the Rasmol XYZ file to the current atom set.
void FIoAtomSet::AddAtomsFromXyzFile(std::string const &FileName, FXyzLoadOptions const *pOtherOptions)
{
// #ifdef USE_FSTREAM_IN_XYZIO
//    std::ifstream
//       str(FileName.c_str());
//    if (!str.good())
//       throw FXyzLoadException("Could not read file from disk", FileName);
// #else
//    ct::TArray<char>
//       pFileContent;
//    if (false == ct::LoadFileIntoMemory(pFileContent, FileName))
//       throw FXyzLoadException("Could not read file from disk", FileName);
//    std::stringstream
//       str(&pFileContent[0], std::stringstream::in);
// #endif // USE_FSTREAM_IN_XYZIO

   MAKE_ISTREAM_FOR_FILE(str, FileName);
   AddAtomsFromXyzStream(str, FileName, pOtherOptions);
   SetName(FileName);
}


// Adds the atoms at the current position of the .xyz input stream to current atom set
// (typically called when a new frame starts)
void FIoAtomSet::AddAtomsFromXyzStream(std::istream &str, std::string const &FileName, FXyzLoadOptions const *pOtherOptions)
{
   // the .xyz file format is quite popular because it is simple... at least in principle.
   //
   // Format is as follows:
   // - first line: number of atoms =: N
   // - second line: some title or information
   // - N x [Element Name Xcoord Ycoord Zcoord]. Elements may be separated
   //   by whitespace. Following the lines some other data may or may not stand.
   //
   // In practice it turned out to be a bit more messy. It can contain additional gradient information,
   // and/or multiple different geometry frames.
   FXyzLoadOptions
      LoadOptions;
   if (0 != pOtherOptions) {
      LoadOptions = *pOtherOptions;
   }

   bool
      bNumAtomsSet = false,
      bProcessCurrentLine = false;
   size_t
      nAtoms = 0;
   std::string
      CurLine;
//    xout << "PROCESS FRAME: Cl0 = '" << CurLine << "'"<< std::endl;

   // read first line. Should be number of atoms, OR, in some cases, the first
   // line of the xyz file itself (in this case there is no pre-determined number
   // of atoms and no comment line).
   while (str.good()) {
      if (std::getline(str, CurLine)) {
         ct::Trim(CurLine);
         // ignore empty lines.
         if (!CurLine.empty())
            break;
      } else {
         if (str.eof())
            // encountered only empty lines while scanning for the header...
            // that means the file is empty.
            return;
         throw FXyzLoadException("Encountered stream error while looking for header line", FileName);
      }
   }

//    xout << "PROCESS FRAME: Cl1 = '" << CurLine << "'"<< std::endl;
   {
      std::stringstream ss(CurLine);
      ss >> nAtoms;
      if (!ss.fail()) {
         // okay, it starts with an integer, apparently.
         bNumAtomsSet = true;

         // read comment line.
         std::getline(str, CurLine);
         // check if there is anything in there which looks like
         // an energy. We take the first thing which cleanly converts
         // to double (or would it be better to take the last?).
         std::stringstream ss2(CurLine);
         SetCaption(CurLine);
         std::string MaybeEnergy;
         while (ss2.good()) {
            // this is indeed the most compact C++ way I could come up with
            // for checking if something "cleanly" converts to 'double'...
            // (and even this is not really clean, because strtod uses "errno",
            // which is not threadsafe).
            ss2 >> MaybeEnergy;
            char const *pBeg = MaybeEnergy.c_str();
            char *pEnd;
            double E = std::strtod(pBeg, &pEnd);
            if (pEnd - pBeg == ptrdiff_t(MaybeEnergy.size())) {
               SetLastEnergy(E);
               break;
            }
         }
      } else {
         // doesn't start with an integer. So guess this line is part of
         // the data already. Remember to not discard the line which we
         // already read.
         bProcessCurrentLine = true;
      }
   }
   bool
      GradSet = false;

//    xout << "PROCESS FRAME: nAtSet? " << bNumAtomsSet << "  nAt = " << nAtoms << std::endl;
   for (size_t iAt = 0; !bNumAtomsSet || iAt < nAtoms; ){
      bool GetLineWorked = true;
      if (bProcessCurrentLine) {
         // first line --- process the line in CurLine which we already loaded
         // for putatively dealing with the file header (which turned out to be
         // not actually there)
         bProcessCurrentLine = false;
      } else {
         // read next data line
         GetLineWorked = bool(std::getline(str, CurLine));
         ct::Trim(CurLine);
      }
//       if (str.bad()) {
      if (!GetLineWorked) {
         if (bNumAtomsSet)
            // if a given number of atoms was set than there'd
            // better be enough of them.
            throw FXyzLoadException("There are fewer atom data lines than were declared.", FileName);
         break;
      }

//       xout << "L: '" << CurLine << "'" << std::endl;

      // we put in some lentinence: we will ignore empty lines and lines
      // starting with typical comment characters (all of which would be invalid
      // in the data body of a regular .xyz file).
      ct::StripLineComment(CurLine, "#");
      ct::StripLineComment(CurLine, "//");
      if (CurLine.empty())
         continue;
//       xout << "  >'" << CurLine << "'" << std::endl;

      // So... let's interpret it as actual data.
      std::stringstream
         ss(CurLine);
      std::string
         AtomType;
      double
         x, y, z;
      ss >> AtomType >> x >> y >> z;
      if ( ss.bad() ) {
         // if number of atoms was set in the file, then there'd better
         // be the right number of atoms in there.
         if (bNumAtomsSet)
            throw FXyzLoadException("Data line not understood (should be: <Element> <x> <y> <z> [...])", CurLine, FileName);
         // if the number was not set, we just reached EOL or some other
         // stuff and just quit here.
         break;
      }

      double
         q(0.), gx(0.), gy(0.), gz(0.);
      ss >> q >> gx >> gy >> gz; // check if there is some other stuff (mysterious double, gradient)
      if (iAt == 0 && !ss.bad())
         GradSet = true; // so we have some gradient information.
      if (iAt != 0 && (!ss.bad() != GradSet))
         throw FXyzLoadException("Data line inconsistent: Some data lines contain gradient info, and some do not.", CurLine, FileName);
      if (!GradSet) {
         // if there is none on the first atom, we don't want it on any others, either.
         q = 0.; gx = 0.; gy = 0.; gz = 0.;
      }

//       xout << "LINE: " << CurLine << "  GradSet: " << GradSet << "  Grad:" << gx << " " << gy << " " << gz << std::endl;

      try {
         FVector3
            vPos = LoadOptions.ConvertToInternalCoords(FVector3(x,y,z));
         FIoAtom
            At(vPos, LoadOptions.GetNormalizedAtomType(AtomType));
         if (GradSet)
            At.vGrad = LoadOptions.ConvertToInternalGradient(FVector3(gx,gy,gz));
         AddAtom(At);
      } catch (std::runtime_error &e) {
         throw FXyzLoadException(e.what(), CurLine, FileName);
      }
      ++ iAt;
//       str.ignore( 0xffff, '\n');
   }
//    xout << "PROCESS FRAME: -> Out" << std::endl;
}


void AddFramesFromMultiXyz(FIoAtomSetList &Frames, std::string const &FileName_, FXyzLoadOptions const *pOtherOptions)
{
   bool
      LoadAllFrames = true;
   std::string
      FileName = FileName_;
   string_slice
      slFileName(FileName.begin(), FileName.end());

   // check if we were asked to only load a subset of the input frames.
   typedef std::set<ptrdiff_t>
      FFrameIdSet;
   FFrameIdSet
      iFramesToAdd;
   if (1) {
      // if given, extract the specified frame ID.
      // This would be specified as, for example,
      //
      //    statpt.xyz:123
      //
      // or
      //
      //    statpt.xyz:[0,1,2]
      //
      // where the part before ':' denotes the xyz file and
      // the part after ':' the frame index/indices (starting at 0).

      string_slice::FPair
         slp1 = slFileName.rsplit1(':');
//       std::cout << fmt::format(":slp1 '{}' '{}'", slp1.first.to_str(), slp1.second.to_str()) << std::endl;
      string_slice
         slFrameIds = slp1.second;
      slFrameIds.trim();
      if (!slFrameIds.empty() && slFrameIds.find_c(ct::is_path_separator_c()) == slFrameIds.last) {
         // is this a list of frames?
         if (slFrameIds.startswith('[')) {
            ct::long_split_result
               srFrameIdList;
            slFrameIds.split_list(srFrameIdList);
            for (size_t iSrFrameId = 0; iSrFrameId != srFrameIdList.size(); ++ iSrFrameId) {
               iFramesToAdd.insert(srFrameIdList[iSrFrameId].to_int());
            }
         } else {
            // nope. Assume just a single integer representing the frame ID..
            iFramesToAdd.insert(slFrameIds.to_int());
         }
         LoadAllFrames = false;
         FileName = slp1.first.to_str();
      } else {
         // do nothing. ':' is not at the end of the string.
         // Might be a directory separator or something else.
      }
   }


//    ct::TArray<char>
//       pFileContent;
//    if (false == LoadFileIntoMemory(pFileContent, FileName))
//       throw FXyzLoadException("Could not read file from disk", FileName);
//    std::stringstream
//       str(&pFileContent[0], std::stringstream::in);
   MAKE_ISTREAM_FOR_FILE(str, FileName);
   size_t
      iFrame = 0;
   while (str.good()) {
      FIoAtomSetPtr
         pAtoms(new FIoAtomSet());
      pAtoms->AddAtomsFromXyzStream(str, FileName, pOtherOptions);
      pAtoms->SetName(fmt::format("{}:{}", FileName, iFrame));

      if (LoadAllFrames || iFramesToAdd.find(iFrame) != iFramesToAdd.end()) {
         iFramesToAdd.erase(iFrame);
         if (!pAtoms->empty())
            Frames.push_back(pAtoms);
      }
      iFrame += 1;
   }

   if (!iFramesToAdd.empty()) {
      fmt::MemoryWriter w;
//       w << fmt::format("LoadMultiXyz: While processing input '{}': ", FileName_);
      w << "Failed to load frame(s) ";
      for (FFrameIdSet::const_iterator itLostFrame = iFramesToAdd.begin(); itLostFrame != iFramesToAdd.end(); ++ itLostFrame) {
         if (itLostFrame != iFramesToAdd.begin())
            w << ", ";
         w << *itLostFrame;
      }
      w << " because input .xyz file contains only " << iFrame << " frames total.";
      w << " (note: frame indices are 0-based!)";
      throw FXyzLoadException(w.str(), FileName_);
   }
}

double FIoAtomSet::GetTotalGradientSq() const
{
   double g = 0;
   for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt)
      g += ct::LengthSq(m_Atoms[iAt].vGrad);
   return g;
}

bool FIoAtomSet::HaveGradients() const
{
   for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt)
      if (m_Atoms[iAt].vGrad != FVector3(0., 0., 0.))
         return true;
   return false;
}

void FIoAtomSet::PrintAsXyz(std::ostream &out, std::string const &ExtraDesc, FXyzPrintOptions const &Options) const
//                             std::string const &Lead, double DistanceFactor, double GradientFactor, bool IncludeGradients)
//       void PrintAsXyz(std::ostream &out, FXyzPrintOptions const &Options = FXyzPrintOptions);
{
   bool
      IncludeGradients = Options.HasPrintFlag(FXyzPrintOptions::XYZPRINT_IncludeGradients);
   if (!HaveGradients())
      IncludeGradients = false;

//    fmt::MemoryWriter w;
//    w.write("E = {:18.12f}", GetLastEnergy());
//    std::string Desc = w.str();
   std::string
      Desc = Options.MakeDescLine(GetCaption(), ExtraDesc, GetLastEnergy(), std::sqrt(GetTotalGradientSq()));

   out << fmt::format("{}{}\n{}{}\n", Options.m_OutputLineLead, size(), Options.m_OutputLineLead, Desc);
   for (size_t i = 0; i < size(); ++i) {
      FIoAtom const
         &At = m_Atoms[i];
      FVector3
         vPosXyz = Options.ConvertToOutputCoords(At.vPos);
      out << fmt::format("{}{:<4s} {:16.10f} {:16.10f} {:16.10f}",
            Options.m_OutputLineLead, At.Type, vPosXyz[0], vPosXyz[1], vPosXyz[2]);
      if (IncludeGradients) {
         FVector3
            vGradXyz = Options.ConvertToOutputGradient(At.vGrad);
         out << fmt::format("  0.00 {:12.7f} {:12.7f} {:12.7f}",
            vGradXyz[0], vGradXyz[1], vGradXyz[2]);
      }
      out << "\n";
   }
   out.flush();
}


FIoAtomSetPtr LoadXyzFrameFromFile(std::string const &FileName_, xyz_io::FXyzLoadOptions const *pOtherOptions)
{
   FIoAtomSetList
      Frames;
   AddFramesFromMultiXyz(Frames, FileName_, pOtherOptions);
   if (Frames.size() != 1)
      throw FXyzLoadException(fmt::format("expected a single-frame .xyz file, but this one has {} frames.", Frames.size()), FileName_);
   return Frames.front();
}

} // namespace xyz_io
