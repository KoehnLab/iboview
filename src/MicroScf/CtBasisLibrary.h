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

#ifndef CT_BASIS_LIBRARY_H
#define CT_BASIS_LIBRARY_H

#include <vector>
#include <string>
#include <iosfwd>

#include "CtBasisShell.h"
#ifdef IR_ECP
   #include "IrEcp.h"
#else
   namespace ir {
      struct FAtomEcp {
         double fElecAbsorbed;
      }; // dummy.
   } // namespace ir
#endif

namespace ct {

struct FBasisSet;

struct FBasisSetLibraryImpl;

struct FBasisSetLibrary
{
public:
   // imports a .libmol file into memory. If given, a ``xxx loaded'' message
   // is prited to pxout.
   void ImportMolproLib(std::string const &FileName, std::ostream *pxout=0);

   /// return false if failed. Otherwise: *ADDs* basis functions to the
   /// array, does not replace them.
   ///   - iAtomIdx: index of atom in atom set. Will be forwarded (unchecked)
   ///     to the created FBasisShell objects.
   ///   - iElement: nuclear charge of atom.
   void LoadBasisFunctions(std::vector<FBasisShell> &Shells,
            int iElement, std::string const &BasisDesc,
            FVector3 const &vAtomPos, int iAtomIdx) const;

   ir::FAtomEcp const *LoadEcp(int iElement, std::string const &EcpDesc, FVector3 const &vAtomPos, int iAtomIdx) const;

   bool HaveBasis(int iElement, std::string const &BasisDesc) const;
   bool HaveEcp(int iElement, std::string const &EcpDesc) const;

   FBasisSetLibrary();
   ~FBasisSetLibrary();
protected:
   FBasisSetLibraryImpl
      *p;
private:
   FBasisSetLibrary(FBasisSetLibrary const&); // not implemented
   void operator = (FBasisSetLibrary const&); // not implemented
};

extern FBasisSetLibrary
   g_BasisSetLibrary; // <- no point in having more than one of these things around.

} // namespace ct

#endif // CT_BASIS_LIBRARY_H
