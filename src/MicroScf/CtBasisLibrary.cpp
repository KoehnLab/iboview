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

// This file is released under the GNU General Public License ("GPL", version 2)
// as part of the CT8K program. Copying, modification, creating derivative works
// and redistribution of the file is allowed, but _only_ subject to the terms
// of that GPL. You should have received a version of this license along
// with this source code. The program comes "as is", without any kind of
// warranty.
//
// Authors/Copyright holders:
//  - Gerald Knizia, 2006 (tag: cgk, contact: cgk.d@gmx.net)

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm> // for std::stable_sort / std::sort
#include <utility> // for std::pair
#include <memory>
#include <set>
#include <map>
#include <list>
#include <set>
#include <sstream>
#include <iostream> // for cerr.
#include <cctype> // for isdigit

// #include "CtCommon.h"
#include "CtBasisLibrary.h"
#include "CxIo.h"
#include "CxParse1.h" // <- atm only for the legacy functions (TrimRight, tolower, etc)... but would be better if we'd convert this to string slice completely.
#include "CtAtomSet.h" // for ElementNumberFromName
#ifdef IR_ECP
   #include "IrEcp.h"
   namespace ct {
      using ir::FAtomEcp;
      using ir::FAtomEcpPtr;
   }
#endif // IR_ECP

// #define BASIS_LIB_DBG


namespace ct {
typedef std::vector<double>
   FScalarArray;
typedef std::vector<int>
   FIntArray;
typedef std::list<std::string>
   FStringList;

typedef unsigned int
   uint;

struct FBasisSetLibraryImpl
{
protected:
   FBasisSetLibraryImpl
      *p;

   struct FBasisKey {
      std::string const
         *pName; // pointer into this->BasisNames.
      uint
         iElement; // element number
      bool operator < (FBasisKey const &other) const {
         if (pName < other.pName) return true;
         if (pName > other.pName) return false;
         return iElement < other.iElement;
      }
   };
   FBasisKey MakeKey(std::string const &Name, uint iElement, bool AssertExists = true) const;
   void TryMakeBasis(std::string const &BasisDesc, int iElement);

   typedef std::multimap<FBasisKey, FAtomShellPtr>
      FBasisFnMap;
   FBasisFnMap
      m_BasisFns;
//    typedef FBasisFnMap::iterator
//       FBfnIt;
   typedef FBasisFnMap::const_iterator
      FBfnCit;

   FAtomShellPtr ParseBasisFn(std::stringstream &line, std::string &Type, std::istream &str, int iElement);

#ifdef IR_ECP
   typedef std::multimap<FBasisKey, FAtomEcpPtr>
      FEcpFnMap;
   FEcpFnMap
      m_EcpFns;
   typedef FEcpFnMap::const_iterator
      FEcpCit;
   FAtomEcpPtr ParseEcpFn(std::stringstream &line, std::string &Type, std::istream &str, int iElement);
#endif // IR_ECP

   typedef std::set<std::string>
      FBasisNameSet;
   FBasisNameSet
      // list of all names of imported entries (AO sets, Aux sets, Guesses, ECPs...)
      m_BasisNames;

   typedef std::set<FElementShellInfoCptr, TComparePtrTarget<FElementShellInfo> >
      FElementShellInfoSet;
   FElementShellInfoSet
      // contains all unique FElementShellInfo objects we collected. These may be
      // shared across different elements.
      m_ElementShellInfos;
   // either inserts p's target into m_ElementShellInfos and returns p, or returns
   // the previous object equivialent to *p in case one exists.
   FElementShellInfoCptr Inserted(FElementShellInfoCptr p);
public:
   // imports a .libmol file into memory (does not yet write any files)
   void ImportMolproLib(std::string const &FileName, std::ostream *pxout);

   // return false if failed. otherwise: *ADD*s basis functions to the
   // array, does not replace them.
   //   - nAtomIdx: index of atom in atom set allowed to be broken when not used
   //     anyway.
   //   - nElement: nuclear charge of atom.
   void LoadBasisFunctions(std::vector<FBasisShell> &Shells,
            int iElement, std::string const &BasisDesc,
            FVector3 const &vAtomPos, int iAtomIdx) const;
   ir::FAtomEcp const *LoadEcp(int iElement, std::string const &EcpDesc, FVector3 const &vAtomPos, int iAtomIdx) const;
   bool HaveBasis(int iElement, std::string const &Name) const;
   bool HaveEcp(int iElement, std::string const &Name) const;
};


FBasisSetLibrary
   g_BasisSetLibrary;

FBasisSetLibraryImpl::FBasisKey FBasisSetLibraryImpl::MakeKey(std::string const &Name, uint iElement, bool AssertExists) const
{
   FBasisNameSet::const_iterator
      itName = m_BasisNames.find(Name);
   FBasisKey
      r;
   if (itName == m_BasisNames.end()) {
      if (Name.size() > 16)
         // don't ask. Molpro truncates them.
         return MakeKey(Name.substr(0,16), iElement, AssertExists);

      if (AssertExists) {
         std::cerr << "Recognized names:\n";
         _for_each( itName, m_BasisNames )
            std::cerr << "   " << *itName << "\n";
         std::cerr << "\n";
         throw std::runtime_error("FBasisSetLibraryImpl: Basis/ECP entry '" + Name + "' not found.");
      }
      r.pName = 0;
      r.iElement = 0;
      return r;
   }
   r.pName = &*itName;
   r.iElement = iElement;
   return r;
}


FElementShellInfoCptr FBasisSetLibraryImpl::Inserted(FElementShellInfoCptr p)
{
   typedef FElementShellInfoSet::iterator
      FShellInfoSetIt;
   std::pair<FShellInfoSetIt, bool>
      it = m_ElementShellInfos.insert(p);
   return *it.first;
}


void FBasisSetLibraryImpl::ImportMolproLib(std::string const &FileName, std::ostream *pxout)
{
   // how this files look:
   //   Note: Lines with * are comments.
      // He s STO-3G STO3G : 3 1 1.3
      // STO-3G
      // 6.3624214 1.158923 0.31364979 0.15432897 0.53532814 0.44463454
      // Li s STO-3G STO3G : 6 2 1.3 4.6
      // STO-3G
      // 16.119575 2.9362007 0.7946505 0.6362897 0.1478601 0.0480887 0.15432897
      // 0.53532814 0.44463454 -0.09996723 0.39951283 0.70011547
      // Li p STO-3G STO3G : 3 1 1.3
      // STO-3G
      // 0.6362897 0.1478601 0.0480887 0.15591627 0.60768372 0.39195739
      // B s cc-pVDZ VDZ : 9 3 1.9 1.9 9.9
      // cc-pVDZ
      // 4570 685.9 156.5 44.47 14.48 5.131 1.898 0.3329 0.1043 0.000696
      // 0.005353 0.027134 0.10138 0.272055 0.448403 0.290123 0.014322 -0.003486
      // -0.000139 -0.001097 -0.005444 -0.021916 -0.059751 -0.138732 -0.131482
      // 0.539526 0.580774 1
   // interpretation: First Line:
   //      ElementName OrbitalType N*[AlternativeNameI] :
   //            #PrimOrbitals #OrbitalsTheyAreContractedTo N*[ContractFrom.ContractTo]
   //      Any string, usually reference to be cited for the basis set
   //      M*[PrimitiveOrbitalExponents] L*[PrimitiveOrbitalCoefficients]
   //   for each contraction contraction length coefficients are supplied.
   // and data sets like this are also seen:
      // H  P 6-311G2P :   2  0
      // G92 6-311G polarization
      //    .15000000D+01   .37500000D+00
      // H  P 6-311G3P :   3  0
      // G92 6-311G polarization
      //    .30000000D+01   .75000000D+00   .18750000D+00
      // H  S 6-31G** 6-31G* 6-31G :   4  1      1.3
      // G92 6-31G**
       // .18731137D+02   .28253944D+01   .64012169D+00   .16127776D+00   .33494604D-01
      // .23472695D+00   .81375733D+00
   // all/some primitive orbitals are probably left uncontracted. Therefore no
   // contraction coefficients are supplied for these (only exponents).

   // Also I originally tought that the names of the basis set after the first
   // one were alternative Names of the set (like cc-pVDZ = VDZ). This is
   // however not how it works, at least not in all cases: It seems that names
   // which are supplied actually mean that the currently described elements
   // have to be inserted into the basis sets of ALL names provided, which
   // may or may not be identical basis sets (for example, for the s and p
   // orbitals of the 931G basis set, also the starred names are listed, which
   // means that these basis sets share these orbitals, altough in general they
   // are different).

   // so.. let's begin the mess.
   // read the entire file into a stringstream object (we need to modify it).
   TArray<char>
      pFileContent;
   if (!LoadFileIntoMemory(pFileContent, FileName))
      throw std::runtime_error("FBasisSetLibraryImpl: Failed to open file '"+FileName+"' for basis library import.");

   // okay, now this is very bad. This FORTRAN stuff (sometimes) denotes
   // exponents in scientific notation not as "23423e-3" but as "23423D-03". We
   // do a lame attempt to convert that. This might break in some situations.
   // FIXME: correct this somehow.
   // FIXME: seriously... this is super-dangerous! At the very least it should
   // be done only for the numerical blocks inside the parsing functions.
   for (size_t i = 1; i < pFileContent.size() - 2; ++i) {
      if ((pFileContent[i-1]=='.' || std::isdigit(pFileContent[i-1])) &&
          pFileContent[i]=='D' &&
          (std::isdigit(pFileContent[i+1]) ||
           ((pFileContent[i+1]=='+' || pFileContent[i+1]=='-') && std::isdigit(pFileContent[i+2]))))
         pFileContent[i] = 'E';
   }
//    for (uint i = 0; i < pFileContent.size() - 2; ++i) {
//       if (pFileContent[i]=='D' &&
//           (pFileContent[i+1]=='+' || pFileContent[i+1]=='-') &&
//           pFileContent[i+2]=='0')
//          pFileContent[i] = 'E';
//    }

   FBasisNameSet
      AllBasisNames;
   std::stringstream
      str(&pFileContent[0], std::stringstream::in);
   // ^- NOTE: "in" means from stream, into local data, not
   // into stream (file naming convention)

   std::string
      ProcessedEntry = "(none)";
//    try {
      while(str.good())
      {
         // clear exception mask, now stream will not throw() when something
         // unexpected happens.
         str.exceptions(std::ios::goodbit);
         std::string
            s;
         std::getline(str, s);
         ProcessedEntry = s;
         TrimRight(s);
         if (s.size() == 0 || s[0] == '*' || s[0] == '!')
            // empty or comment line, throw it away and go on with the next
            continue;
         if (!str.good()) // eof, bad etc.
            break;
         str.exceptions(std::ios::failbit);
         // ^- when something fails, throw an exception. this will happen
         // if the actual file format does not match the one I had in mind
         // when coding this.

         // expected format: ElementName[w]OrbitalType[w]AlternativeNames[w] :
         // using namespace std;
         // cout << "** read line: '" << s << "'" << endl;
         std::stringstream
            line(s, std::stringstream::in);
         line.exceptions(std::ios::badbit | std::ios::failbit);
         std::string
            Element, Type;
         line >> Element >> Type;
         Type = tolower(Type);
         int
            iElement = ElementNumberFromName(Element.c_str());

         FStringList
            BasisNames; // all names of basis sets in which the
                        // current entry is to be inserted.
         for (line >> s; s != ":"; line >> s) {
      //             std::cout << "Alternative Name:" << s << std::endl;
            BasisNames.push_back(tolower(stripwhitespace(s)));
            AllBasisNames.insert(stripwhitespace(s));
         }

         // import all names of the basis function/ECP
         FStringList::const_iterator
            itName;
         _for_each(itName, BasisNames)
            m_BasisNames.insert(*itName);

         // read & make the actual basis function or ECP object, and link it to all
         // all the names it occurs under.
         try {
            if (Type != "ecp") {
               FAtomShellPtr
                  pBfn = ParseBasisFn(line, Type, str, iElement);
               _for_each(itName, BasisNames)
                  m_BasisFns.insert(FBasisFnMap::value_type(MakeKey(*itName, iElement), pBfn));
            } else {
#ifdef IR_ECP
               FAtomEcpPtr
                  pEcp = ParseEcpFn(line, Type, str, iElement);
               _for_each(itName, BasisNames)
                  m_EcpFns.insert(FEcpFnMap::value_type(MakeKey(*itName, iElement), pEcp));
#else
               // compiled without ECP support. Ignore entry.
            throw std::runtime_error("Encountered ECP in libmol, but this version of CtBasisLibrary was compiled without ECP support.");
#endif // IR_ECP
            }
         } catch (std::runtime_error &e) {
            std::stringstream ss;
            ss << "During parse of '" << ProcessedEntry << "': " << e.what();
            throw std::runtime_error(ss.str());
         }

         // chew the EOL marker if present, leave loop otherwise.
         str.exceptions(std::ios::goodbit);
         str.ignore(0xbad, '\n');
      };
//    } catch (std::ios_base::failure &e){
//       // this is not exactly something i would usually
//       // call "error handling" but i really hate this string-
//       // fiddling stuff and we can't really do anything better
//       // about it anyway.
//       std::cerr << "PARSER EXCEPTION:" << e.what() << std::endl;
//       throw std::runtime_error("Parsing of LIBMOL file FAILED because the actual syntax did not match the expected one. Failed entry: '" + ProcessedEntry + "'");
//    } catch (std::exception &e){
//       std::cerr << "Exception during LibmolFile parsing: " << e.what() << std::endl;
//       throw;
//    };


   // if provided, write some imporant looking comments about what
   // we loaded to the standard output. Makes things look so much more sciency!
   if (pxout) {
      std::ostream
         &xout = *pxout;
      std::size_t
         iDirSep = FileName.rfind('/');
      if (iDirSep == std::string::npos)
         iDirSep = 0;
      else
         iDirSep += 1;
      xout << fmt::format(" Loaded {:<25}", FileName.substr(iDirSep));
      if (true) {
         xout << "[";
         FBasisNameSet::const_iterator
            itSet;
         size_t
            nLen = 0;
         _for_each(itSet, AllBasisNames) {
            if (nLen >= 40) {
               xout << ",...";
               break;
            }

            if (itSet != AllBasisNames.begin())
               xout << ", ";
            xout << *itSet;
            nLen += itSet->size();
         }
         xout << "]";
      }
      xout << std::endl;
   }
}

FAtomShellPtr FBasisSetLibraryImpl::ParseBasisFn(std::stringstream &line, std::string &Type, std::istream &str, int iElement)
{
#ifdef BASIS_LIB_DBG
   std::cout << "Reading basis fn: '" << line.str() << "... ";
   std::cout.flush();
#endif
   if (Type.size() != 1)
      throw std::runtime_error(fmt::format("Parsing error, cannot interpret orbital type '{}'", Type));

   int
      AngMom;
   std::vector<std::pair<int,int> >
      Cos;
   std::vector<double>
      Exps;
   char
      cAngMom = ::tolower(Type[0]);
   for (AngMom = 0; AngMom < 9; ++ AngMom)
      if (cAngMom == "spdfghikl"[AngMom])
         break;
   if (AngMom == 9)
      throw std::runtime_error(fmt::format("Failed to understand angular momentum '{:c}'.", cAngMom));


   // expected format: #prim orbitals #contractions (#contr.)*[a.b]
   // denoting indices of begin and end of a contraction with the
   // following exponents/contraction coefficients.
   int
      nExp,
      nCo,
      nCoeff(0), // total number of contraction coefficients to read (in all contractions).
      nHighestExpInCo(0); // 1-based index.
   line >> nExp;
   // some files have no number of contractions declared---we are
   // supposed to interpret them as pure primitive sets.
   if (line.eof())
      nCo = 0;
   else
      line >> nCo;

//    try {
//       line >> nCo;
//    } catch (std::ios_base::failure &e){
//       nCo = 0;
//    }
//          std::cout << "#Prim " << nExp << " #Co " << nCo << std::endl;
   Cos.reserve(nCo);
   for (int i = 0; i < nCo; ++ i){
      std::pair<int,int>
         iCo;
      char Dot;
      line >> iCo.first >> Dot >> iCo.second;
      iCo.first -= 1; // convert to 0-based [begin,end).
      if (Dot != '.')
         throw std::runtime_error("GTO-Contraction read format error.");
//             std::cout << "  Co: #" << iCo.first << "-#" << iCo.second << std::endl;
      if (iCo.second <= iCo.first || iCo.second > nExp)
         throw std::runtime_error("GTO-Contraction logical error.");
      nCoeff += iCo.second - iCo.first;
      nHighestExpInCo = std::max(nHighestExpInCo, iCo.second);
      Cos.push_back(iCo);
   }

   std::string
      EntryComment;
   do { // read name, maybe skip comments.
      getline(str, EntryComment);
   } while (!EntryComment.empty() && EntryComment[0] == '*');
   // cout << "Entry Comment: " << EntryComment << endl;

   // now read exponents and contraction coefficients;
   // (this will break if comments are present in between)
   Exps.resize(nExp);
   std::vector<double>
      Coeffs(nCoeff, 0);
   for (int i = 0; i < nExp; ++ i){
      str >> Exps[i];
      // cout << "Exp: " << Exps.back() << endl;
   }

   // read in contraction coefficients.
   if (nCo != 0)
      for (int i = 0; i < nCoeff; ++ i){
         double Coeff;
         str >> Coeff;
         // cout << "Coeff: " << Coeff << endl;
         Coeffs[i] = Coeff;
      }
   // copy over the contraction coefficients to the contractions.
   std::vector<double>
      CoMatrix(Exps.size() * nCo, 0.);
   int
      iCoeff = 0;
   for (int i = 0; i < nCo; ++ i){
      for (int j = Cos[i].first; j != Cos[i].second; ++ j)
         CoMatrix[j + i*Exps.size()] = Coeffs[iCoeff + j - Cos[i].first];
      iCoeff += (Cos[i].second - Cos[i].first);
   }

   // in some files some primitive orbitals are left uncontracted.
   // but these are not stored as 1-GTO contractions but the
   // coefficients for these are just not present in the file.
   // Make 1-GTO contractions for them.
   int
      nAdditionalCo = nExp - nHighestExpInCo;
   if (0 != nAdditionalCo)
   {   // generate 1-GTO-each contractions manually.
      int nCoExplicit = nCo;
      nCo += nAdditionalCo;
      CoMatrix.resize(Exps.size() * nCo, 0.);
      int
         iStart = nHighestExpInCo;
         // ^- 0 based index, the rhs one is 1-based.
      for (int i = 0; i < nAdditionalCo; ++ i) {
         int iCo = i + nCoExplicit;
         CoMatrix[(i + iStart) + Exps.size() * iCo] = 1.;
      }
   }

#ifdef BASIS_LIB_DBG
   std::cout << "  ... ok." << std::endl;
#endif

   // make and return the actual basis function object.
   FAtomShellPtr
      pBf(new FAtomShell(AngMom, &Exps[0], Exps.size(), &CoMatrix[0], CoMatrix.size()/Exps.size()));
   if (!EntryComment.empty()) {
      FElementShellInfoCptr
         pInfo(new FElementShellInfo("", "", EntryComment));
      pBf->AttachInfo(iElement, Inserted(pInfo));
   }

   return pBf;
}

#ifdef IR_ECP
static size_t _AsSize(double f) {
   size_t i = size_t(f);
   if (double(i) != f) {
      std::stringstream ss;
      ss << "Expected integer, but found " << f << ".";
      throw std::runtime_error(ss.str());
   }
   return i;
}

FAtomEcpPtr FBasisSetLibraryImpl::ParseEcpFn(std::stringstream &line, std::string &Type, std::istream &str, int iElement)
{
   if (Type != "ecp")
      throw std::runtime_error( "Parsing error, cannot interpret ECP type '" + Type + "'" );
#ifdef BASIS_LIB_DBG
   std::cout << "Reading ecp: '" << line.str() << "... ";
   std::cout.flush();
#endif

   // expected format:
   //   Au ECP ECP60MDF : 60 5 4 91
   // ! -> 60 core electrons,
   //      lmax=5 for semi-local part (from l=0 to 5, inclusive; the last term l=5 is the LOCAL term (not s=0)),
   //      lmax=4 for SO part (from l=1 to 4, inclusive, no local term),
   //      total number of numerical parameters following (in this case 91).

   int
      nEcpElec, AngMomMaxScalar, AngMomMaxSpinOrbit, nParamsTotal;

   line >> nEcpElec >> AngMomMaxScalar >> AngMomMaxSpinOrbit >> nParamsTotal;
   if (line.fail())
      throw std::runtime_error("Failed to interpret ECP declaration. Expected: <nEcpElec> <lMax> <lMaxSo> <nParams>");
   if (nParamsTotal < 0 || AngMomMaxScalar < 0 || AngMomMaxSpinOrbit < 0)
      throw std::runtime_error("Failed to interpret ECP declaration. lMax, lMaxSo, and nParams cannot be negative.");

   std::string
      EntryComment;
   do { // read comment (at least one line), maybe skip comments (then additional lines).
      getline(str, EntryComment);
   } while (EntryComment.size() != 0 && EntryComment[0] == '*');
   // cout << "Entry Comment: " << EntryComment << endl;

   // read numerical parameters.
   std::vector<double>
      AllParams;
   AllParams.reserve(size_t(nParamsTotal));
   for (int i = 0; i < nParamsTotal; ++ i) {
      double f;
      str >> f;
      if (str.fail()) {
         std::stringstream ss;
         ss << "Expected " << nParamsTotal << " scalar parameters, but failed to read parameter #" << (AllParams.size() + 1) << ".";
         throw std::runtime_error(ss.str());
      }
      AllParams.push_back(f);
   }

   // setup the output object.
   FAtomEcpPtr
      pEcp = new FAtomEcp();
   pEcp->fElecAbsorbed = double(nEcpElec);
   pEcp->fEffectiveCoreCharge = double(iElement) - pEcp->fElecAbsorbed;

   FAtomEcp::FEcpShellFnList
      &ShellFnsScalar = pEcp->ShellFns,
      &ShellFnsSo = pEcp->ShellFnsSpinOrbit;
   FAtomEcp::FParamList
      &Params = pEcp->Params;
   Params.reserve(AllParams.size()); // we actually need a few less. But that's okay.
   double
      *pOrigParamsStart = &Params[0];
   ShellFnsScalar.reserve(1 + AngMomMaxScalar); // 0 to lmax (but lmax is the LOCAL part, which is added to everything, including l=0...infty)
   ShellFnsSo.reserve(AngMomMaxSpinOrbit); // 1 to lmax (has no local part)

   size_t
      iParamOff = 0;
   for (int il = 0; il < 1 + AngMomMaxScalar + AngMomMaxSpinOrbit; ++ il) {
      // note: l for scalar part goes from 0 to AngMomMaxScalar (inclusive),
      //       l for SO part goes from 1 to AngMomMaxSpinOrbit (inclusive);
      FAtomEcp::FEcpShellFnList
         *pShellFns = (il <= AngMomMaxScalar)? &ShellFnsScalar : &ShellFnsSo;
      ir::FEcpShellFn
         EcpFn;
//       EcpFn.l = (il <= AngMomMaxScalar)? il : (il - AngMomMaxScalar + 1);
//       if (il == AngMomMaxScalar)
//          EcpFn.l = -1; // local part.

      if (iParamOff >= AllParams.size())
         throw std::runtime_error("Expected and required number of scalar parameters does not match.");

      // read number of exponents (= number of contractions).
      EcpFn.nExp = _AsSize(AllParams[iParamOff++]);
      if (iParamOff + 3*EcpFn.nExp > AllParams.size())
         throw std::runtime_error("Expected and required number of scalar parameters does not match.");
      EcpFn.pExp = &Params[Params.size() + 0*EcpFn.nExp];
      EcpFn.pCo  = &Params[Params.size() + 1*EcpFn.nExp];
      EcpFn.pPow = &Params[Params.size() + 2*EcpFn.nExp];
      Params.resize(Params.size() + 3*EcpFn.nExp);
      for (size_t iExp = 0; iExp < EcpFn.nExp; ++ iExp) {
         const_cast<double*>(EcpFn.pPow)[iExp] = AllParams[iParamOff++] - 2.;
         // ^- @ -2.: They store "exponent + 2" in the libmol file:
         //   See: http://www.tc.uni-koeln.de/cgi-bin/pp.pl?language=en,job=getreadme
         //   That may be related to the r^2 in the 4 \pi r^2 \d r of the volume
         //   element, which we absorb in grid weights and should not do again.
         const_cast<double*>(EcpFn.pExp)[iExp] = AllParams[iParamOff++];
         const_cast<double*>(EcpFn.pCo)[iExp] =  AllParams[iParamOff++];
      }

      EcpFn.l = (il <= AngMomMaxScalar)? (il - 1) : (il - AngMomMaxScalar + 1);
      if (0) {
         // apparently the normal ECPs in Molpro's .libmol files do not actually have a local part.
         // This construction here applies the highest-AM part to the local part.
         // But my undestanding is now that this is not actually what the ECP descs mean.
         // So this part is disabled and the local parts are now used as verbatim.
         //
         // I leave this in because the actual reference treatment appears very fishy to me.
         // doing this one directly is also not quite right, however, because the higher l
         // ecp entries are supposed to specify the *DIFFERENCE* between the l-level ECP and
         // the local ECP. And using the highest-l ECP only (without adding these subtractions)
         // would result in some ECP contributions occuring twice.
         if (il == 0) {
            // explicit local part(?), which is apparently ignored in Molpro ECP descs.
            assert(EcpFn.l == -1);
            if (EcpFn.nExp != 1 || EcpFn.pCo[0] != 0.0)
               throw std::runtime_error(fmt::format("encountered unexpected local part of ECP type '{}' for element '{}'", Type, iElement));
            // don't store this.
            continue;
         }
         if (il == AngMomMaxScalar) {
            EcpFn.l = -1; // max am part to be used as local part.
         }
      } else {
         // leave as is. First one (local/lmax) part is zero now.
         // We just skip it. TODO: ask Peterson if that possibly can be correct...  seems shady.
         if (EcpFn.nExp == 1 || EcpFn.pCo[0] == 0.0)
            continue;
      }

      if (0) {
         if (!(EcpFn.l == 0 || EcpFn.l == 1 || EcpFn.l == 2)) {
//          if (!(EcpFn.l == 0 || EcpFn.l == 1)) {
//          if (!(EcpFn.l == 0 || EcpFn.l == 0)) {
//          if (EcpFn.l >= 3) {
            std::cout << fmt::format("* WARNING: skipping component l = {} of ECP type '{}' for element '{}'\n", EcpFn.l, Type, iElement);
            continue;
         }
      }

      pShellFns->push_back(EcpFn);
   }
   std::sort(ShellFnsScalar.begin(), ShellFnsScalar.end());
   std::sort(ShellFnsSo.begin(), ShellFnsSo.end());
   if (pOrigParamsStart != &pEcp->Params[0])
      // I think this should not happen in practice, but I think this is not guaranteed by the standard.
      // And if it does, then all the internal pointers in the EcpFns are now wrong.
      // At least produce a diagnostic message in this case.
      throw std::runtime_error("unexpected reallocation. Fix hacky ECP parse.");

#ifdef BASIS_LIB_DBG
   std::cout << "  ... ok." << std::endl;
#endif
   return pEcp;
}
#endif // IR_ECP




// converts something like "sto3g>>spd" into "sto3g" and ['s','p','d'].
// First two parameters are Out parameters.
void ParseBasisDesc(std::string &BasisName, std::set<char> &OrbitalTypes,
   FBasisDesc const &BasisDesc)
{
   OrbitalTypes.clear();
   std::size_t
      nSeparatorIdx = BasisDesc.find(">>");
   // separator present in string?
   if ( nSeparatorIdx == std::string::npos ){
      // no->primitive declaration, "import all orbitals from [...]".
      BasisName = tolower(stripwhitespace(BasisDesc));
   } else { // yes. find out which orbitals to import.
      BasisName = tolower(stripwhitespace(BasisDesc.substr(0,nSeparatorIdx)));
      for ( std::size_t i = nSeparatorIdx+2; i < BasisDesc.size(); ++i ){
         char c = ::tolower( BasisDesc[i] );
         if ( ( c >= 'a' ) && ( c <= 'z' ) )
            // ^- we should consider using more b- and q-orbitals.
            OrbitalTypes.insert(c);
      }
   }
}

bool starts_with(std::string const &s, char const *p)
{
   std::size_t
      N = s.size();
   for ( std::size_t i = 0; i < N && *p && s[i] == *p; ++ i, ++p ){
   }
   return *p == 0;
}

int AngMomFromChar(char c) {
   switch (c) {
      case 's': return 0;
      case 'p': return 1;
      case 'd': return 2;
      case 'f': return 3;
      case 'g': return 4;
      case 'h': return 5;
      case 'i': return 6;
      case 'k': return 7;
      case 'l': return 8;
      default: {
         std::stringstream str;
         str << "Angular momentum not recognized: '" << c << "'. Should be one of spdfghik.";
         throw std::runtime_error(str.str());
      }
   };
}

void FBasisSetLibraryImpl::TryMakeBasis(std::string const &BasisDesc, int iElement)
{
   if (starts_with(BasisDesc, ".et[")) {
      // syntax: .ET[spdf,<center>, <ratio>, <powmin>, <powmax>]
      char
         AngMoms[20];
      double
         fCenter, fRatio;
      int
         iPowMin, iPowMax, iScan;
      iScan = std::sscanf(BasisDesc.c_str(), ".et[%10[spdfghik],%lf,%lf,%i,%i]",
         &AngMoms[0], &fCenter, &fRatio, &iPowMax, &iPowMin);
      if (!((iScan == 4 && iPowMin > 0) || iScan == 5)) {
         std::stringstream str;
         str << "Failed to understand even tempered basis declaration '" << BasisDesc << "'."
            << " Expected: .et[<spd...>, <center>, <ratio>, <powmax>(, <powmin>)].";
         throw std::runtime_error(str.str());
      }
      if (iScan == 4)
         iPowMin = iPowMax;
      iPowMin *= -1;

      // make the basis function and store it.
      std::vector<double>
         Exps;
      Exps.reserve(iPowMax - iPowMin + 1);
      for (int i = iPowMax; i >= iPowMin; -- i)
         Exps.push_back(fCenter * std::pow(fRatio, (double)i));

      m_BasisNames.insert(BasisDesc);

      for (char *pAm = &AngMoms[0]; *pAm != 0; ++ pAm) {
         FAtomShellPtr
            pBfn(new FAtomShell(AngMomFromChar(*pAm), &Exps[0], Exps.size()));
         m_BasisFns.insert(FBasisFnMap::value_type(MakeKey(BasisDesc, iElement), pBfn));
      }
   } else if (starts_with(BasisDesc, "aug-" ) ||
            starts_with(BasisDesc, "daug-" ) ||
            starts_with(BasisDesc, "taug-" )) {
//       xout << "  in .aug generation code." << std::endl;
      uint
         iAugLevel = 1;
      if (starts_with(BasisDesc, "daug-"))
         iAugLevel = 2;
      if (starts_with(BasisDesc, "taug-"))
         iAugLevel = 3;
      // look up the parent basis.
      std::string
         ParentDesc(BasisDesc.substr((iAugLevel==1)? 4 : 5)); // includes the '-'.
      std::pair<FBfnCit,FBfnCit>
         itBfs = m_BasisFns.equal_range(MakeKey(ParentDesc, iElement));
      if (itBfs.first == m_BasisFns.end())
         // not there.
         return;

      m_BasisNames.insert(BasisDesc);

      uint
         AmMask = 0; // bitmask of angular momenta already handled.
      for (FBfnCit it = itBfs.first; it != itBfs.second; ++ it){
         // note: this code assumes that all functions within a angular
         // momentum are stored inside one single contracted shell!
         FAtomShellCptr
            pFn = it->second;
         if ((AmMask & (1 << pFn->AngMom)) != 0) {
            std::stringstream str;
            str << "sorry, diffuse augmentation code got confused. Encountered multiple shells with angular momentum l=" << pFn->AngMom << ".";
            throw std::runtime_error(str.str());
         }
         AmMask |= 1 << pFn->AngMom;

         // make exponents.
         TArray<double>
            Exps = pFn->Exponents;
         if (Exps.empty())
            throw std::runtime_error("Encountered a basis function without exponents during augmentation.");
         double
            fCen = Exps.back(),
            fRatio = 1./2.5;
         if (Exps.size() > 1) {
            fRatio = Exps[Exps.size()-1]/Exps[Exps.size()-2];
         }
         for (uint iAug = 0; iAug < iAugLevel; ++ iAug)
            Exps.push_back(fCen * std::pow(fRatio, (double)(1+iAug)));
         // make new contraction matrix.
         TArray<double>
            CoMatrix(Exps.size() * (pFn->nCo() + iAugLevel), 0.);
         // copy original data. Note that the number of exponents has changed,
         // so we need a strided copy.
         for (uint iCo = 0; iCo < pFn->nCo(); ++ iCo)
            for (uint iExp = 0; iExp < pFn->nExp(); ++ iExp)
               CoMatrix[iExp + iCo * Exps.size()] = pFn->CoMatrix[iExp + pFn->nExp() * iCo];
         // make coefficients for new primitives.
         for (uint iAug = 0; iAug < iAugLevel; ++ iAug)
            CoMatrix[Exps.size() * (pFn->nCo() + iAug) + iAug] =
               1./ir::RawGaussNorm(Exps[pFn->nExp() + iAug + 1], pFn->AngMom);
         FAtomShellPtr
            pBfn(new FAtomShell(pFn->AngMom, &Exps[0], Exps.size(),
               &CoMatrix[0], CoMatrix.size()/Exps.size(), FAtomShell::TYPE_Unnormalized));
         m_BasisFns.insert(FBasisFnMap::value_type(MakeKey(BasisDesc, iElement), pBfn));
      }
   }
}


struct FSortPredAngMom
{
   bool operator() (FBasisShell const &a, FBasisShell const &b) const {
      return a.l() < b.l();
   }
};


void FBasisSetLibraryImpl::LoadBasisFunctions(std::vector<FBasisShell> &Shells,
   int iElement, std::string const &BasisDesc_,
   ct::FVector3 const &vAtomPos, int iAtomIdx) const
{
   if (BasisDesc_ == "none")
      // user explicitly requested no basis functions for this atom
      return;

   // todo: tokenize BasisDesc at ';'s here (and strip whitespace).
   // do the same as we he have below for the individual parts.

   std::string BasisDesc = BasisDesc_;
   if (BasisDesc == "cc-pvtz-ig")
      BasisDesc = "cc-pvtz";
   if (BasisDesc == "def2-tzvp-ig")
      BasisDesc = "def2-tzvp";
   if (BasisDesc == "def2-sv(p)-jfit-ig")
      BasisDesc = "def2-sv(p)-jfit";

   FBasisKey
      Key = MakeKey(BasisDesc, iElement, false);
   std::pair<FBfnCit,FBfnCit>
      itBfs = m_BasisFns.equal_range(Key);
   if (Key.pName == 0 || itBfs.first == itBfs.second) {
      // basis not yet officially registered, but maybe it is one we
      // can make by modifying another one (e.g., aug-cc-pVTZ from cc-pVTZ)
      const_cast<FBasisSetLibraryImpl*>(this)->TryMakeBasis(BasisDesc, iElement);
      Key = MakeKey(BasisDesc, iElement);
      itBfs = m_BasisFns.equal_range(Key);
   }

   if (itBfs.first == itBfs.second) {
      std::stringstream str;
      str << "Requested basis '" << BasisDesc << "' for element '"
          << ElementNameFromNumber(iElement) << "' was not found.";
      throw std::runtime_error(str.str());
   }

   // note that we do NOT clear Shells but rather append only.
   size_t
      iShellBeg = Shells.size();
   for (FBfnCit it = itBfs.first; it != itBfs.second; ++ it){
      FAtomShellCptr
         pFn = it->second;
      it->second->Finalize();
      Shells.push_back(FBasisShell(vAtomPos, iAtomIdx, pFn, Key.pName));
   }
   // re-order basis shells by angular momentum, in case they have
   // not been specified in order of angular momentum in the basis set library
   // (e.g., because some shells were shared with other basis sets, and some
   // were not; happens in minao.libmol)
   std::stable_sort(Shells.begin() + iShellBeg, Shells.end(), FSortPredAngMom());
}

ir::FAtomEcp const *FBasisSetLibraryImpl::LoadEcp(int iElement, std::string const &EcpName, FVector3 const &vAtomPos, int iAtomIdx) const
{
#ifdef IR_ECP
   std::string
      NormalizedName = tolower(stripwhitespace(EcpName));
   std::pair<FEcpCit, FEcpCit>
      itEcps = m_EcpFns.equal_range(MakeKey(NormalizedName, iElement, false));
   if (itEcps.first == itEcps.second) {
      std::stringstream str;
      str << "Requested ECP '" << EcpName << "' for element '"
          << ElementNameFromNumber(iElement) << "' was not found.";
      throw std::runtime_error(str.str());
   }
   -- itEcps.second;
   if (itEcps.first != itEcps.second) {
      std::stringstream str;
      str << "More than one ECP of name '" << EcpName << "' registered for element '"
          << ElementNameFromNumber(iElement) << "'.";
      throw std::runtime_error(str.str());
   }

   return itEcps.first->second.get();

   IR_SUPPRESS_UNUSED_WARNING(vAtomPos);
   IR_SUPPRESS_UNUSED_WARNING(iAtomIdx);
#else
   return 0;
   IR_SUPPRESS_UNUSED_WARNING(iElement);
   IR_SUPPRESS_UNUSED_WARNING(EcpName);
   IR_SUPPRESS_UNUSED_WARNING(vAtomPos);
   IR_SUPPRESS_UNUSED_WARNING(iAtomIdx);
#endif
}

bool FBasisSetLibraryImpl::HaveBasis(int iElement, std::string const &Name) const
{
   std::string
      NormalizedName = tolower(stripwhitespace(Name));
   std::pair<FBfnCit, FBfnCit>
      itBfs = m_BasisFns.equal_range(MakeKey(NormalizedName, iElement, false));
   return itBfs.first != itBfs.second;
}

bool FBasisSetLibraryImpl::HaveEcp(int iElement, std::string const &Name) const
{
#ifdef IR_ECP
   std::string
      NormalizedName = tolower(stripwhitespace(Name));
   std::pair<FEcpCit, FEcpCit>
      itEcps = m_EcpFns.equal_range(MakeKey(NormalizedName, iElement, false));
   return itEcps.first != itEcps.second;
#else
   return false;
#endif
}


FBasisSetLibrary::FBasisSetLibrary()
   : p(new FBasisSetLibraryImpl())
{}

FBasisSetLibrary::~FBasisSetLibrary() {
   delete p;
}

void FBasisSetLibrary::ImportMolproLib(std::string const &FileName, std::ostream *pxout) {
   return p->ImportMolproLib(FileName, pxout);
}

void FBasisSetLibrary::LoadBasisFunctions(std::vector<FBasisShell> &Shells,
   int iElement, std::string const &BasisDesc,
   FVector3 const &vAtomPos, int iAtomIdx) const
{
   return p->LoadBasisFunctions(Shells, iElement, BasisDesc, vAtomPos, iAtomIdx);
}



bool FBasisSetLibrary::HaveBasis(int iElement, std::string const &BasisDesc) const
{
   return p->HaveBasis(iElement, BasisDesc);
}

bool FBasisSetLibrary::HaveEcp(int iElement, std::string const &EcpDesc) const
{
   return p->HaveEcp(iElement, EcpDesc);
}

ir::FAtomEcp const *FBasisSetLibrary::LoadEcp(int iElement, std::string const &EcpDesc, FVector3 const &vAtomPos, int iAtomIdx) const
{
   return p->LoadEcp(iElement, EcpDesc, vAtomPos, iAtomIdx);
}




} // namespace ct







// kate: indent-width 3; indent-mode normal;
