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

#ifndef _CT8K_DFT_FUNCTIONAL_H
#define _CT8K_DFT_FUNCTIONAL_H

// #include "Ir.h"
#include "CxTypes.h"
#include <string>

namespace ir {
   struct FIntegralKernel;
} // namespace ir


namespace ct {

struct FXcFunctionalImpl;



struct FXcFunctional : public FIntrusivePtrDest
{
   explicit FXcFunctional(std::string const &Name);
   ~FXcFunctional();

   bool NeedSigma() const; // GGA?
   bool NeedTau() const; // tau-style m-GGA?
   bool NeedUpsilon() const; // laplace-style m-GGA?
   bool TranslateUpsilonToTau() const;

   // input: rhoc, rhoo, sigmacc, sigmaco, sigmaoo49
   // output: zk, vrhoc, vrhoo, vsigmacc, vsigmaco, vsigmaoo (case nDiff == 1)
   // output: zk only (case nDiff == 0)
   void Eval(double *pOut, size_t nOutSt, double const *pIn, size_t nInSt, size_t nPts, size_t nDiff) const;
   // return derivative order:
   // 0 == LDA, 1 == GGA or m-GGA (tau-style), 2 == m-GGA (upsilon-style)
   uint nOrder() const;

   // factor for Hartree-Fock exchange (commonly called "exact exchange") to use
   // with this functional, relative to Hartree-Fock exchange (so a "1.0" here
   // means "same factor as in Hartree-Fock"... which is -0.5 total).
   double ExchFactor() const;
   bool NeedsExactExch() const { return ExchFactor() != 0.; }
   // if != 0, returns the integral kernel to use for computing the Hartree-
   // Fock-exchange-like contributions this functional is associated with.
   // if = 0, either there is no exact exchange, or the exact exchange is
   // computed with the Coulomb kernel.
   ir::FIntegralKernel const *pExchKernel() const;

   std::string const &Name() const;
   std::string Desc() const;
private:
   FXcFunctionalImpl
      *p;
private:
   FXcFunctional(FXcFunctional const &other); // not implemented.
   void operator = (FXcFunctional const &other); // not implemented.
};

typedef TIntrusivePtr<FXcFunctional>
   FXcFunctionalPtr;

// typedef std::pair<FXcFunctionalPtr, FScalar>
//    // maps a dft functional to its prefactor.
//    FXcFunctionalSet;



} // namespace ct

#endif // DFT_FUNCTIONAL_H
