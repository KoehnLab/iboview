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

#ifndef CT_WFDECL_H
#define CT_WFDECL_H

#include "CxTypes.h" // for size_t and ptrdiff_t
#include <string>

namespace ct {

// represents the declaration of a type of wave function --- mostly charge and spin.
// Atm silently assumes basic RHF/RKS type of calculation. Might be extended by more concrete
// information (e.g., active orbitals, spatial symmetry, other kinds of restrictions)
struct FWfDecl
{
   int
      iCharge;
   unsigned
      Ms2; // total number of open shells (=2 * total-M quantum number of electrons)
   std::string
      FileName_WfRead,
      FileName_WfWrite;

   // total number of electrons (can only be called after SetNuclearCharge has been called!)
   unsigned nElec() const;
   // number of up-spin electrons (alpha-electrons)
   unsigned nElecA() const { return (nElec() + Ms2)/2; }
   // number of up-spin electrons (alpha-electrons)
   unsigned nElecB() const { return (nElec() - Ms2)/2; }
   // number of open-shell orbitals
   unsigned nOpen() const { return Ms2; } // nOccA - nOccB
   // number of closed-shell orbitals
   unsigned nClosed() const { return nElecB(); } // nOccA - nOccB

//    int iCharge() const { return nElec - m_nNuclearCharge; }
   explicit FWfDecl(int iCharge_=0, unsigned Ms2_=0);
   explicit FWfDecl(std::string const &WfDesc);

   // set the total nuclear charge (sum of atomic numbers of atoms - ECP charge).
   // Required to compute the total number of electrons from the relative charge iCharge.
   void SetNuclearCharge(unsigned NuclearCharge_) { m_nNuclearCharge = int(NuclearCharge_); };

   void SetArgs(std::string const &CommandOptions);
   void SetCharge(ptrdiff_t iCharge_);
   void SetMs2(ptrdiff_t iCharge_);

   bool HasEqualChargeSpinSym(FWfDecl const &other) const;
protected:
   int m_nNuclearCharge;

   void AssertNotInstanciated();
};

} // namespace ct


#endif // CT_WFDECL_H
