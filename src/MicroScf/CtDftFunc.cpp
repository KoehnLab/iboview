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

#include <string>
#include <iostream>
#include <vector>
#include <sstream>

#include "format.h"
#include "Ir.h"
#include "CtDftFunc.h"
#include "CxIo.h"
#include "CxParse1.h" // for toupper

#include "xc_meta.h"
typedef void (*FXcFn1)(double &E, double *dE, double const *D, double const Factor);


namespace ct {

struct FXcFunctionalEntry
{
   char const
      *pName;
   FXcFn1
      pFn1; // function to evaluate value + 1st derivative of xc functional.
   double
      Factor;
   FXcFunctionalEntry(char const *pName_, FXcFn1 pFn1_, double Factor_) : pName(pName_), pFn1(pFn1_), Factor(Factor_) {}
};

struct FXcFunctionalImpl
{
   explicit FXcFunctionalImpl(std::string const &Name);
   ~FXcFunctionalImpl();

   void Eval(double *pOut, size_t nOutSt, double const *pIn, size_t nInSt, size_t nPts, size_t nDiff) const;

   bool
      NeedTau,
      NeedSigma,
      NeedUpsilon,
      TranslateUpsilonToTau;
   bool
      m_SetupFailed;
   double
      // factor for Hartree-Fock exchange (commonly called "exact exchange") to
      // use, relative to Hartree-Fock exchange (so a "1.0" here means "same
      // factor as in Hartree-Fock"... which is -0.5 total).
      m_ExchFactor;
   ir::FIntegralKernelPtr
      // integral kernel to use for "exact exchange"-type contributions
      // to Fock matrix, if different from Coulomb kernel
      // (e.g. for functionals coupled to special short-range/long-range
      // exchange)
      m_pExchKernel;
   std::string
      // formal reference name of the functional (in correct capitalization)
      m_Name;

   inline uint nInputLength() const;

   std::vector<FXcFunctionalEntry>
      Components;
   std::string Desc() const;
};


uint FXcFunctionalImpl::nInputLength() const
{
   assert(!NeedTau || NeedSigma); // tau only supported if sigma is also given.
   if (NeedUpsilon) return 2 + 3 + 2 + 2; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo, tauc, tauo, upsilonc, upsilono
   // ^- (note: places for tau also reserved if upsilon is used, but tau is not)
   if (NeedTau) return 2 + 3 + 2; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo, tauc, tauo
   if (NeedSigma) return 2 + 3; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo
   return 2; // rhoc, rhoo
}


FXcFunctionalImpl::FXcFunctionalImpl(std::string const &Name_)
   : m_SetupFailed(false), m_ExchFactor(0), m_pExchKernel(0)
{
   NeedSigma = false;
   NeedTau = false;
   NeedUpsilon = false;
   TranslateUpsilonToTau = false;
   m_Name = toupper(Name_);

   // note: may need to change m_Name below to correct captialization,
   // unless the official name is all-caps. Or other changes may be
   // needed.
   if (m_Name == "LDA" || m_Name == "CK16") {
      Components.push_back(FXcFunctionalEntry("UEG-X", rx::xc_diracx, 1.0));
      Components.push_back(FXcFunctionalEntry("CK16-C", rx::xc_ldac_ck16, 1.0));
      m_Name = "CK16";
   } else if (m_Name == "C16") {
      Components.push_back(FXcFunctionalEntry("UEG-X", rx::xc_diracx, 1.0));
      Components.push_back(FXcFunctionalEntry("C16-C", rx::xc_ldac_c16, 1.0));
      m_Name = "C16";
   } else if (m_Name == "LDA" || m_Name == "PW92") {
      Components.push_back(FXcFunctionalEntry("UEG-X", rx::xc_diracx, 1.0));
      Components.push_back(FXcFunctionalEntry("PW92-C", rx::xc_pw92c, 1.0));
   } else if (m_Name == "UEG-X" || m_Name == "UEG-X" || m_Name == "DIRACX" || m_Name == "DIRAC-X") {
      // uniform electron gas, exchange only.
      Components.push_back(FXcFunctionalEntry("UEG-X", rx::xc_diracx, 1.0));
   } else if (m_Name == "PBE") {
      Components.push_back(FXcFunctionalEntry("PBE-X", rx::xc_pbex, 1.0));
      Components.push_back(FXcFunctionalEntry("PBE-C", rx::xc_pbec, 1.0));
      NeedSigma = true;
   } else if (m_Name == "PBE0") {
      m_ExchFactor = 0.25; // 10.1063/1.478522 <- funky that there should be a reference for that...
      Components.push_back(FXcFunctionalEntry("PBE-X", rx::xc_pbex, 0.75));
      Components.push_back(FXcFunctionalEntry("PBE-C", rx::xc_pbec, 1.0));
      NeedSigma = true;
   } else if (m_Name == "PBEX" || m_Name == "PBE-X") {
      Components.push_back(FXcFunctionalEntry("PBE-X", rx::xc_pbex, 1.0));
      NeedSigma = true;
   } else if (m_Name == "MN15-L") {
      Components.push_back(FXcFunctionalEntry("MN15-L", rx::xc_mn15l, 1.0));
      NeedSigma = true;
      NeedTau = true;
   } else if (m_Name == "LL-MN15-L") {
      Components.push_back(FXcFunctionalEntry("MN15-L", rx::xc_mn15l, 1.0));
      NeedSigma = true;
      NeedUpsilon = true;
      TranslateUpsilonToTau = true;
   } else if (m_Name == "TPSS") {
      Components.push_back(FXcFunctionalEntry("TPSS-X", rx::xc_tpssx, 1.0));
      Components.push_back(FXcFunctionalEntry("TPSS-C", rx::xc_tpssc, 1.0));
      NeedSigma = true;
      NeedTau = true;
   } else if (m_Name == "LL-TPSS") {
      if (1) {
         Components.push_back(FXcFunctionalEntry("TPSS-X", rx::xc_tpssx, 1.0));
         //          Components.push_back(FXcFunctionalEntry("PBEX", rx::xc_pbex, 1.0));
         Components.push_back(FXcFunctionalEntry("TPSS-C", rx::xc_tpssc, 1.0));
         NeedSigma = true;
         NeedUpsilon = true;
         TranslateUpsilonToTau = true;
      } else {
         // note: this does NOT need tau! (but it requires upsilon)
         Components.push_back(FXcFunctionalEntry("LL-TPSS", rx::xc_ll_tpss, 1.0));
         NeedSigma = true;
         NeedUpsilon = true;
         TranslateUpsilonToTau = false;
      }
   } else if (m_Name == "LL-TPSSX" || m_Name == "LL-TPSS-X") {
      Components.push_back(FXcFunctionalEntry("TPSS-X", rx::xc_tpssx, 1.0));
      NeedSigma = true;
      NeedUpsilon = true;
      TranslateUpsilonToTau = true;
   } else if (m_Name == "TPSSX" || m_Name == "TPSS-X") {
      Components.push_back(FXcFunctionalEntry("TPSS-X", rx::xc_tpssx, 1.0));
      NeedSigma = true;
      NeedTau = true;
      //    } else if ( m_Name == "PC07" ) {
      //       // note: that's a kinetic energy functional, not a xc functional (here for testing purposes).
      //       Components.push_back(FXcFunctionalEntry("PC07", rx::xc_ek_pc07, 1.0));
      //       NeedSigma = true;
      //       NeedUpsilon = true;
   } else if (m_Name == "NONE" || m_Name.empty() || m_Name == "NULL") {
      // no DFT functional. Zero DFXC.
      m_Name = "null";
   } else {
      m_SetupFailed = true;
      throw std::runtime_error(fmt::format("FXcFunctional: failed to set up for xc functional '{}'", Name_));
   }
}

std::string FXcFunctionalImpl::Desc() const
{
   std::stringstream str;
   str << " Density functional:         ";
   if (m_SetupFailed) {
      str << "ERROR/XC-SETUP-FAILED";
   } else {
      char const
         *pCompFmt =  "  {} ={:6.3f}",
         *pOtherFmt = "  {} = {}";
      for (uint iComp = 0; iComp < Components.size(); ++iComp) {
         FXcFunctionalEntry const &e = Components[iComp];
         str << fmt::format(pCompFmt, e.pName, e.Factor);
      }
      if (m_ExchFactor != 0.) {
//          str << fmt::format(pCompFmt, "HF-EXCH", m_ExchFactor);
         str << fmt::format(pCompFmt, "EXCH", m_ExchFactor);
//          str << fmt::format(pCompFmt, "EXACT", m_ExchFactor);
         if (m_pExchKernel)
            str << fmt::format(pOtherFmt, "EXCH-KER", m_pExchKernel->Desc());
      }
      if (TranslateUpsilonToTau) {
//          str << "  TAU = PC07";
         str << fmt::format(pOtherFmt, "TAU", "PC07");
      }
   }
   return str.str();
}



void FXcFunctionalImpl::Eval(double *pOut, size_t nOutSt, double const *pIn, size_t nInSt, size_t nPts, size_t nDiff) const
{
   size_t
      nInp = nInputLength(),
      nOut = 1 + nInp; // Energy + dEnergy/dInputs.
   assert(nOutSt >= nOut && nInSt >= nInp);
   assert_rt(nDiff == 1);
   for (size_t iPt = 0; iPt < nPts; ++ iPt) {
      double const
         *pInPt = &pIn[nInSt * iPt];
      double
         *pOutPt = &pOut[nOutSt * iPt];
      for (size_t iOut = 0; iOut < nOut; ++ iOut)
         pOutPt[iOut] = 0;
      for (size_t iComp = 0; iComp < Components.size(); ++ iComp) {
         FXcFunctionalEntry const &e = Components[iComp];
         e.pFn1(pOutPt[0], &pOutPt[1], pInPt, e.Factor);
         // ^- this adds to the previous results.
      }
   }
}


FXcFunctionalImpl::~FXcFunctionalImpl()
{
}

FXcFunctional::FXcFunctional(std::string const &Name)
   : p(new FXcFunctionalImpl(Name))
{
}

FXcFunctional::~FXcFunctional()
{
   delete p;
}

void FXcFunctional::Eval(double *pOut, size_t nOutSt, double const *pIn, size_t nInSt, size_t nPts, size_t nDiff) const
{
   return p->Eval(pOut, nOutSt, pIn, nInSt, nPts, nDiff);
}

bool FXcFunctional::NeedSigma() const
{
   return p->NeedSigma;
}

bool FXcFunctional::NeedTau() const
{
   return p->NeedTau;
}

bool FXcFunctional::NeedUpsilon() const
{
   return p->NeedUpsilon;
}

bool FXcFunctional::TranslateUpsilonToTau() const
{
   return p->TranslateUpsilonToTau;
}


uint FXcFunctional::nOrder() const
{
   if (NeedUpsilon())
      return 2;
   return NeedSigma()? 1 : 0;
}

std::string FXcFunctional::Desc() const
{
   return p->Desc();
}

double FXcFunctional::ExchFactor() const {
   return p->m_ExchFactor;
}

ir::FIntegralKernel const *FXcFunctional::pExchKernel() const {
   return p->m_pExchKernel.get();
}

std::string const &FXcFunctional::Name() const {
   return p->m_Name;
}


} // namespace ct.



// #define XC_FUNC_TEST
#ifdef XC_FUNC_TEST

static void PrintArray1(char const *pTitle, double *pValues, uint nValues, int nStep = -1, char const *pFmt = "{:11.5f}")
{
   using ct::xout;
   if (nStep == -1) {
      // find a reasonable step size such that we get coverage of the entire
      // array but don't print out more than ~10-15 numbers.
      nStep = (12+nValues)/13;
   }
   xout << fmt::format("[{:5}] {:>20}", nValues, pTitle);
   for ( uint i = 0; i < nValues; i += nStep )
      xout << format(pFmt, pValues[i]);
   xout << "\n";
}


int main_xc_func_test(int argc, char *argv[])
{
   using namespace ct;
   FXcFunctionalPtr
//       pRef = new FXcFunctional("DIRAC");
//       pRef = new FXcFunctional("PW92C");
//       pRef = new FXcFunctional("PBEX");
      pRef = new FXcFunctional("PBEC");

   double In[5] = {0.123, -0.05, 0.07, 0.01, 0.02}; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo
//    double In[5] = {0.123, -0.05, 0., 0., 0.}; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo
//    double In[5] = {0.123, -0.05, 1e-8, 0., 1e-8}; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo
//    double In[5] = {0.123, 0., 0.05, 0., 0.}; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo
//    double In[5] = {0.01, 0., 0.01, 0., 0.}; // rhoc, rhoo, sigmacc, sigmaco, sigmaoo
   double OutRef[6] = {0}; // energy and potentials.
   double OutTst[6] = {0}; // energy and potentials.

   pRef->Eval(&OutRef[0], 6, &In[0], 5, 1, 1);

   typedef void (*FXcFn1)(double &E, double *dE, double const *D, double const Factor);
   FXcFn1 xcfn1 = rx::xc_pbec;
   xcfn1(OutTst[0], &OutTst[1], &In[0], 1.0);

   uint nVar = pRef->NeedSigma()? 5 : 2;
   PrintArray1("Inputs", &In[0], nVar, 1, "{:16.8f}");
   PrintArray1("Output (ref)", &OutRef[0], nVar+1, 1, "{:16.8f}");
   PrintArray1("Output (tst)", &OutTst[0], nVar+1, 1, "{:16.8f}");
   return 0;
}

#endif // XC_FUNC_TEST
