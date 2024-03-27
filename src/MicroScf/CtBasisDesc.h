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

#ifndef CT_BASISDESC_H
#define CT_BASISDESC_H

#include <map>
#include <string>
#include <stdexcept>

#include "CxTypes.h"
#include "CxVec3.h"

namespace ir {
   struct FAtomEcp;
}

namespace ct {

enum FBasisContext {
   // main orbital basis: the primary basis set used to expand molecular
   // orbitals in (e.g., def2-TZVPP or cc-pVTZ)
   BASIS_Orbital,
   // fitting basis set for 3-center density fitting for computing SCF Coulomb
   // integral operators (j) and/or coulomb and exchange-correlation potential
   // operators (j + vxc). For use in DF-RKS with pure DFT functionals (e.g.,
   // CK16, PBE, TPSS) with or without auxiliary xc expansion (DF-J or DF-JX).
   //
   // Example basis sets:
   // - univ-JFIT (F. Weigend: PCCP 8, 1057 (2006); 10.1039/B515623H)
   //
   // For DF-JX, see:
   //    Bienvenu, Knizia: J. Chem. Theory Comput. 2018, 14, pp.1297-1303
   //    10.1021/acs.jctc.7b01083
   BASIS_JFit,
   // fitting basis set for 3-center density fitting for computing SCF Coulomb
   // (j) and Hartree-Fock exchange (k) integral operators.
   // For use in density-fitted Kohn-Sham with hybrid functionals (e.g., PBE0,
   // B3LYP) and in density fitted actual Hartree-Fock (DF-RKS, DF-RHF).
   //
   // Example basis sets:
   // - univ-JKFIT (F. Weigend: J. Comput. Chem. 29, 167 (2008); 10.1002/jcc.20702)
   BASIS_JkFit,
   // fitting basis set for 3-center density fitting for computing K^ij_ab
   // integral operator matrices, as used in most correlated wave function
   // methods. These are not only used in MP2 --- it's called MP2FIT because the
   // K^ij_ab operators the sets are typically designed for are meant for
   // computing MP2-like exchange integral operators.
   BASIS_Mp2Fit,
   BASIS_CcsdFit,
   // resolution of the identitity (RI) basis sets for constructing auxiliary
   // complementary basis sets (CABS) in explicitly correlated wave function
   // methods (such as DF-MP2-F12 or RHF-UCCSD(T)-F12a).
   //
   // Note:
   // - Sometimes the density fitting approximation used to compute integral
   //   operators is also is also called "resolution of the identity" approximation
   //   (e.g., RI-RHF instead of DF-RHF). However, that is something *different*.
   // - The current RI basis sets are unrelated to integral fitting.
   // - The F12RI context may be split into a CARI (complementary auxiliary RI)
   //   and and DRI (direct RI) context at a later point.
   BASIS_F12RI,
   // A minimal orbital basis used for constructing SCF initial guesses.
   BASIS_Guess,
   // A minimal orbital basis used as reference orbital set for interpretative
   // purposes (in particular: for IAO construction). Used only for explicit
   // overrides; if not set explicitly, defaults to the same basis as
   // BASIS_Guess.
   BASIS_MinAo,
   // note: IAOFIT is a fitting basis for local exchange drivers; it is *not*
   // the MINAO reference for making IAOs (we actually use either BASIS_MinAo
   // or BASIS_Guess for those).
   BASIS_IaoFit,
   BASIS_Multipoles
};


typedef std::string
   FBasisDesc;
typedef std::string
   FEcpDesc;
typedef std::map<FBasisContext, FBasisDesc>
   FBasisDescs;


std::string GetBasisContextName(FBasisContext Context);

// returns whether an explicit orbital basis needs to be assigned and supplied for
// FindDefaultBasis to do its job, or whether this is a context for which it can
// likely find a reasonable default without knowing the orbital basis it goes with.
bool NeedOrbitalBasisForContextDefault(FBasisContext Context);

// attempt to find a default basis set for the given context which goes with OrbBasis
// for the given element. Theoretically some extra information could be required to
// make the best choice (e.g., the exact type of ECP), but for the moment that is enough.
void FindDefaultBasis(std::string &Out, FBasisContext Context, std::string const &OrbBasis, FBasisDescs const *pOtherBases, int iElement, unsigned nEcpElec);

// void PrintEcpData(std::ostream &out, ir::FAtomEcp const *pIrEcp, FAtom const *pAt, std::string const &Ind = "   ");
void PrintEcpData(std::ostream &out, ir::FAtomEcp const *pIrEcp, int iElement, TVector3<double> const &vPos, std::string const &Ind);

void AutoAdjustBasisAssociations(FBasisDescs &BasisDesc, std::string &EcpDesc, unsigned &nEcpElectrons, unsigned iElement);
size_t iLikelyEcpCharge(int iElement);


} // namespace ct

#endif // CT_BASISDESC_H
