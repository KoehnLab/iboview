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

#ifndef IR_ECP_H
#define IR_ECP_H

#include "Ir.h"
#include <vector> // for storing ECP parameters.
#include "CxTypes.h"

namespace ir {
   /// represents one ECP-type scalar core shell:
   ///
   ///       U_l(r) = \sum_k Co[k] r^pPow[k] Exp[-Exp[k] * r^2]
   ///
   /// which will get combined into one of the integral forms:
   ///
   ///       U_L(\vec r) = U_L(|r|)            (harmless: here stored as l=-1)
   ///
   /// and the super-awkward
   ///
   ///       \sum_{m=-l}^{m=l} |Y_lm> U_l(r) <Y_lm|
   ///
   /// where the projection onto the S_lm represents ANGULAR COORDINATES ONLY (i.e.,
   /// the integal CANNOT be factorized into overlap-type integrals over standard
   /// functions).
   ///
   /// Note: unlike for FRawShell objects, here the center is supplied externally
   /// (i.e., not part of the object).
   struct FEcpShellFn  {
      int
         l;
      double const
         *pExp, *pCo, *pPow;
      size_t
         nExp;
      FEcpShellFn() {}
//       FEcpShell(int l_, double const *pExp_, double const *pCo_, size_t nExp_, double const *vCen_) : l(l_), pExp(pExp_), pCo(pCo_), nExp(nExp_), vCen(vCen_) {}
      FEcpShellFn(int l_, double const *pExp_, double const *pCo_, double const *pPow_, size_t nExp_) : l(l_), pExp(pExp_), pCo(pCo_), pPow(pPow_), nExp(nExp_) {}

      bool IsLocal() const { return l == -1; }
      bool IsGaussian() const { for (size_t i = 0; i < nExp; ++i) if (pPow[i] != 0.) return false; return true; }
      bool operator < (FEcpShellFn const &other) const { return l < other.l; }
   };

   /// represents the parameters of a Köln-Stuttgart-type ECP for a given element,
   /// but *not* its center (center supplied externally). Unlike raw shell objects, this object
   /// stores its parameters itself, so they cannot be trivially provided from the outside. This is
   /// done for practical reasons.
   struct FAtomEcp : public ct::FIntrusivePtrDest {
      typedef std::vector<FEcpShellFn>
         FEcpShellFnList;
      FEcpShellFnList
         // Scalar-relativistic part. sorted by ls. Local part (l=-1) comes first.
         ShellFns,
         // Spin-orbit part. sorted by ls. Local part (l=-1) comes first. Has no l = 0 term.
         ShellFnsSpinOrbit;
      double
         // number of electrons this ECP represents (e.g., 10 for a ECP10MDF)
         fElecAbsorbed,
         // Zeff: Prefactor of the REMAINING coulomb potential (typically iElement - fElecAbsorbed)
         // note: beware of prefactor. It's a core charge, not an electron charge.
         fEffectiveCoreCharge;
      typedef std::vector<double>
         FParamList;
      FParamList
         Params; // the pointer entries in FEcpShellFn point into this array.
      int MaxL() const;
      double MaxPowR() const;
   };
   typedef ct::TIntrusivePtr<FAtomEcp>
      FAtomEcpPtr;

   void EvalEcp(double *pOut, size_t StrideA, size_t StrideB, FRawShell const *pA, FRawShell const *pB,
      double Prefactor, bool Add, FAtomEcp const *pEcp, double const *pEcpCen, ct::FMemoryStack &Mem);


} // namespace ir

#endif // IR_ECP_H
