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

#include <cmath>
#include <algorithm> // for std::max
#include <stdexcept>
#include "Ir.h"
#include "IrInternal.h"

#ifdef IR_INCLUDE_PRINT_FUNCTIONS
   // these are for the printing functions IrPrintRow, IrPrintMatrixGen, and operator << (std::ostream, FIrRawShell const&).
   // If this external dependency is not wanted, the corresponding functions can simply be deleted.
   #include "format.h"
   #include <iostream>
#endif // IR_INCLUDE_PRINT_FUNCTIONS

// #if defined(INCLUDE_OPTIONALS) || defined(IR_ECP)
   #include "IrFactorials.h"
   #include "IrFactorials.inl"
// #endif // INCLUDE_OPTIONALS

namespace ir {


double InvRawRawGaussNorm(double fExp, unsigned l)
{
   // purple book (6.6.14). (factor 4*pi/(2*l+1) converting between Ylm and Slm omitted)
   return std::pow(2*fExp/M_PI,.75) * std::sqrt(std::pow(4.*fExp,l) / GetDoubleFactorialI(2*signed(l)-1));
}

double RawGaussNorm(double fExp, unsigned l)
{
//    return 1./InvRawGaussNorm(fExp, l);
   return std::pow(M_PI/(2*fExp),.75) * std::sqrt(GetDoubleFactorialI(2*signed(l)-1)/std::pow(4.*fExp,l));
}



// those things do not actually carry data, and they are very likely to be needed.
FOverlapKernel
   g_IrOverlapKernel;
FCoulombKernel
   g_IrCoulombKernel;


void FIntegralKernel::EvalTildeGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // first compute standard kernel integrals Gm(rho,T) = (- D/D[T]) G0(rho,T)
   EvalGm(pOut, rho, T, MaxM, Factor);

   // now convert from Gm(rho,T) to \tilde Gm(rho,T) = (2 rho D/D[T]) G0(rho,T)
   // by absorbing powers of (-2 rho).
   //
   // These rho powers would normally be present in the R of a textbook MDRR.
   // But our two-center Shell-MDRR has all of the exponent quantities absorbed
   // in kernel integrals, which is what calls for the modified base kernels
   // (\tilde Gm(rho,T) instead of regular Gm(rho,T))
   double
      RhoPow = 1.;
   for (unsigned i = 0; i < MaxM + 1; ++ i){
      pOut[i] *= RhoPow;
      RhoPow *= -2*rho;
   }
}


void FCoulombKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   IrBoysFn(pOut, T, MaxM, Factor*(2*M_PI)/rho);
}


double FCoulombKernel::MaxRange() const
{
   return 1e30;
}

char const *FCoulombKernel::Desc() const
{
   return "1/r";
}


void FOverlapKernel::EvalGm(double *pOut, double /*rho*/, double T, unsigned MaxM, double Factor) const
{
   pOut[0] = Factor * std::exp(-T);
   for (unsigned i = 1; i <= MaxM; ++i)
      pOut[i] = pOut[i-1];
}

void FOverlapKernel::EvalTildeGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // note: this one is not required (default implementation works fine).
   // but wanted to try if it makes a difference for computing overlap
   // matrices in semi-empirical stuff.
   double
      f = Factor * std::exp(-T);
   for (unsigned i = 0; i <= MaxM; ++i) {
      pOut[i] = f;
      f *= -2*rho;
   }
}


double FOverlapKernel::MaxRange() const
{
   return 0.;
} // it's a delta function.

char const *FOverlapKernel::Desc() const
{
   return "delta(r)";
}


void FRsqKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // This one: Eq. (F20) in the supporting info of
   // J. Chem. Theory Comput. 2020, 16, 4, 2570-2583 https://doi.org/10.1021/acs.jctc.9b01296
   double f = Factor * std::sqrt(pow3(M_PI)/std::pow(rho, 5));
   pOut[0] = f * (T + 1.5);
   pOut[1] = -f; // (-d/d T) G0
   // ...for this one, all the other derivatives just vanish.
   for (unsigned i = 2; i <= MaxM; ++i)
      pOut[i] = 0;
}

double FRsqKernel::MaxRange() const
{
   return 1e30;
} // ^- this one actually *increases* in magnitude with distance...

char const *FRsqKernel::Desc() const
{
   return "r^2";
}





FIntegralKernel::~FIntegralKernel() {}
FCoulombKernel::~FCoulombKernel() {}
FOverlapKernel::~FOverlapKernel() {}
FRsqKernel::~FRsqKernel() {}


// // #ifdef INCLUDE_OPTIONALS
unsigned const MaxJ = 40;

static void EvalErfCoulombGm(double *pOut, double rho, double T, unsigned MaxM, double PrefactorFm, double Omega)
{
   // dx.doi.org/10.1039/b605188j eq.52, 2nd term.
   double
      f = Omega/std::sqrt(Omega*Omega + rho);
   IrBoysFn(pOut, f*f*T, MaxM, f*PrefactorFm);
   double
      Acc = 1.,
//       ifsq = 1./(f*f);
      ifsq = (f*f);
   // hack up the derivatives. Note:
   //    Gm(rho,T) := [-d/dT]^m G0(rho,T)
   // so we just get some factors for the T we put in the erf'd Boys function.
   // The normal one calculates derivatives with respect to (f*f*T),
   // and we need to convert it to derivatives with respect to T.
   // FIXME: shold this be f*f or 1/f*f?
   for (unsigned i = 1; i <= MaxM; ++ i){
      Acc *= ifsq;
      pOut[i] *= Acc;
   }
}

FErfCoulombKernel::FErfCoulombKernel(double Rc_)
   : m_Omega(1./(Rc_*Rc_)), m_Rc(Rc_)
{}

void FErfCoulombKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   EvalErfCoulombGm(pOut, rho, T, MaxM, (2*M_PI)*Factor/rho, m_Omega);
}

double FErfCoulombKernel::MaxRange() const // that's the long range part.
{
   return 1e30;
}

char const *FErfCoulombKernel::Desc() const
{
   return "erf(r/rc)/r";
}


FErfcCoulombKernel::FErfcCoulombKernel(double Rc_)
   : m_Omega(1./(Rc_*Rc_)), m_Rc(Rc_)
{}

void FErfcCoulombKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // form regular coulomb kernel and subtract the long-range kernel.
   double
      Fm[MaxJ],
      FactorFm = (2 * M_PI) * Factor/rho;
   IrBoysFn(&Fm[0], T, MaxM, FactorFm);
   EvalErfCoulombGm(pOut, rho, T, MaxM, FactorFm, m_Omega);
   for (unsigned i = 0; i <= MaxM; ++ i)
      pOut[i] -= Fm[i];
}

double FErfcCoulombKernel::MaxRange() const
{
   return 6.*m_Rc;
} // erfc(6) \approx 2e-17

char const *FErfcCoulombKernel::Desc() const
{
   return "erfc(r/rc)/r";
}


FErfCoulombKernel::~FErfCoulombKernel() {}
FErfcCoulombKernel::~FErfcCoulombKernel() {}



void FGaussGeminal::Finalize()
{
   if (nExp == 0) {
      // I guess we have a '1' kernel then? Actually, not sure... it is
      // not multiplicative.
      assert(0); // <- we probably do not want those.
      m_MaxRange = 1e30;
      return;
   } else {
      if (m_MaxRange <= 0.) {
         // estimate maximum range of the kernel, if not provided from outside.
         //
         // atm: just take the smallest exponent and look where its exp(-z
         // r12^2) becomes zero. Could do this better by considering the
         // contracted kernel, like done for the basis function range, but in
         // this case we should probably provide some holder objects, as doing
         // this has non-negligible cost which would prevent instantiating the
         // kernels on the fly.

         double fMinExp = pExp[0];
         for (size_t i = 1; i < nExp; ++ i)
            fMinExp = std::min(pExp[i], fMinExp);
         //    exp(-z r^2) = 1e-18
         //    -z r^2 = ln(1e-18)
         //    r = sqrt(-ln(1e-18)/z)
         m_MaxRange = sqrt(41.5/fMinExp);
      }
   }
}


FGaussKernel::~FGaussKernel()
{
}

double FGaussKernel::MaxRange() const
{
   return m_pGaussGeminal->MaxRange();
}

char const *FGaussKernel::Desc() const
{
   return "F(r)";
}

void FGaussKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // see PCCP 10 3390 (2008).
   double const
      *IR_RP Omega = &m_pGaussGeminal->pExp[0],
      *IR_RP pCo = &m_pGaussGeminal->pCo[0];
   size_t
      nExp = m_pGaussGeminal->nExp;

   for (size_t m = 0; m < MaxM+1; ++ m)
      pOut[m] = 0;

   for (size_t iExp = 0; iExp < nExp; ++ iExp) {
      double
         f = 1.0/(rho + Omega[iExp]),
         RhoTilde = f * Omega[iExp];

      // [1], eq. (24):
      //    G_m(rho,T) = \sum_i c[i] RhoTilde[i]^m (pi/(rho+omega_i))^(3/2) exp(-RhoTilde[i] * T).
      double
         u = M_PI * f,
         Base = Factor * pCo[iExp] * u * std::sqrt(u) * std::exp(-RhoTilde * T),
         RhoTildeM = 1.0;
      for (size_t m = 0; m <= MaxM; ++ m) {
         pOut[m] += Base * RhoTildeM;
         RhoTildeM *= RhoTilde;
      }
   }
}



FGaussCoulombKernel::~FGaussCoulombKernel()
{
}

double FGaussCoulombKernel::MaxRange() const {
   return m_pGaussGeminal->MaxRange();
}

char const *FGaussCoulombKernel::Desc() const
{
   return "F(r)/r";
}


void FGaussCoulombKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // see PCCP 10 3390 (2008).
   double const
      *IR_RP Omega = &m_pGaussGeminal->pExp[0],
      *IR_RP pCo = &m_pGaussGeminal->pCo[0];
   size_t
      nExp = m_pGaussGeminal->nExp;

   for (size_t m = 0; m < MaxM+1; ++ m)
      pOut[m] = 0;

   size_t const MaxJ = 24;
   double
      Fm[MaxJ+1];
   assert(MaxM <= MaxJ);

   for (size_t iExp = 0; iExp < nExp; ++ iExp ) {
      double
         f = 1.0/(rho + Omega[iExp]),
         RhoTilde = f * Omega[iExp],
         RhoHat = f * rho,
         InvRhoHat = 1.0/RhoHat;

      IrBoysFn(&Fm[0], RhoHat * T, MaxM, Factor*(2*M_PI)*f*std::exp(-RhoTilde * T) * pCo[iExp]);

      double
         RhoHatM = 1.0; // RhoHat^m
      for (size_t m = 0; m <= MaxM; ++ m) {
         double
            r = 0,
            t = RhoHatM,
            u = 1.0;
         // this is supposed to form
         //   r = \sum_k BinomialCoeff(m,k) * RhoTilde**k * RhoHat**(m-k) * Fm[m-k];
         //   (second line of eq. (30))
         for (size_t k = 0; k <= m; ++ k) {
            r += GetBinomialCoeffF(m,k) * u * t * Fm[m-k];
            u *= RhoTilde;
            t *= InvRhoHat;
         }

         RhoHatM *= RhoHat;
         pOut[m] += r;
      }
   }
}


FGaussKineticKernel::~FGaussKineticKernel()
{
}

double FGaussKineticKernel::MaxRange() const {
   return m_pGaussGeminal->MaxRange();
}

char const *FGaussKineticKernel::Desc() const
{
   return "[grad F(r)]^2";
}


void FGaussKineticKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   // see PCCP 10 3390 (2008).
   double const
      *IR_RP Omega = &m_pGaussGeminal->pExp[0],
      *IR_RP pCo = &m_pGaussGeminal->pCo[0];
   size_t
      nExp = m_pGaussGeminal->nExp;

   for (size_t m = 0; m < MaxM+1; ++ m)
      pOut[m] = 0;

   // eq. (32).
   for (size_t j = 0; j < nExp; ++ j) {
      for (size_t i = j; i < nExp; ++ i) {
         double
            wi = Omega[i],
            wj = Omega[j],
            wij = wi + wj,
            f = 1.0/(rho + wij),
            RhoTilde = f * wij,
            RhoHat = f * rho,
            Pref = Factor * 4.0 * wi * wj *
                   f*f*sqrt((M_PI*M_PI*M_PI)*f) *
                   exp(-RhoTilde * T);
         if (i != j)
            Pref *= 2.0;
         double
            t = RhoTilde * (1.5 + RhoHat * T);
         double
            PrefRhoTildeM1 = pCo[i] * pCo[j] * Pref/RhoTilde;
         for (size_t m = 0; m <= MaxM; ++ m) {
            pOut[m] += PrefRhoTildeM1 * (t - double(m) * RhoHat);
            PrefRhoTildeM1 *= RhoTilde;
         }
      }
   }
}


FMetaKernel::FMetaKernel(FIntegralKernel const **pKernels, double const *pFactors, unsigned nKernels)
   : m_pKernels(pKernels), m_pFactors(pFactors), m_nKernels(nKernels)
{
   m_Range = 0.;
   for (unsigned i = 0; i < m_nKernels; ++ i)
      m_Range = std::max(m_Range, m_pKernels[i]->MaxRange());
}

double FMetaKernel::MaxRange() const {
   return m_Range;
}

char const *FMetaKernel::Desc() const {
   return "meta";
}

void FMetaKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   double
      Gm[MaxJ];
   for (unsigned i = 0; i <= MaxM; ++ i)
      pOut[i] = 0;
   for (unsigned iKernel = 0; iKernel < m_nKernels; ++ iKernel) {
      m_pKernels[iKernel]->EvalGm(&Gm[0], rho, T, MaxM, Factor);
      for (unsigned i = 0; i <= MaxM; ++ i)
         pOut[i] += Gm[i];
   }
}

// #endif // INCLUDE_OPTIONALS



#ifdef IR_INCLUDE_PRINT_FUNCTIONS
// some printing stuff which may be helpful in debug mode.
// Not included in default mode because it induces unnecessary dependencies, and
// the host programs probably have own means of printing this kind of stuff anyway.
// (you can enabled them without _DEBUG by setting -DIR_INCLUDE_PRINT_FUNCTIONS)

template<class T>
void IrPrintRow(std::ostream &xout, T const *pData, ptrdiff_t iColSt, size_t nCols, std::string const &FmtF, char const *pEndLn="\n")
{
   for (size_t iCol = 0; iCol < nCols; ++ iCol)
      if (FmtF.empty())
         xout << " " << fmt::format("{:14.6f}", pData[iColSt * iCol]);
      else
         xout << " " << fmt::format(FmtF, pData[iColSt * iCol]);
   xout << pEndLn;
}


void IrPrintMatrixGen(std::ostream &xout, double *pData, size_t nRows, size_t iRowSt, size_t nCols, size_t iColSt, std::string const &Caption, char const *pNumFmt)
{
   xout << fmt::format("  Matrix {}, {} x {}.\n", Caption, nRows, nCols);
//    std::string
//       FmtS = "{:14}",
//       FmtI = "{:>11}   ",
//       FmtF = "{:14.6f}";
   std::string
      FmtF = pNumFmt,
      FmtS = FmtF;
   // make a string format by copying FmtF to FmtS, and then deleting the decimal separator and the
   // floating point format specification (g,f,e) if present.
   size_t
      iDecimal = FmtS.find('.');
   if (iDecimal != std::string::npos) {
      size_t iClosingBrace = FmtS.find('}', iDecimal);
      size_t n = FmtS.size();
      if (iClosingBrace != std::string::npos)
         n = iClosingBrace - iDecimal;
      FmtS.erase(iDecimal, n);
   }

   std::string
      FmtI = FmtS;
   std::string
      IndFmt = "{:10}";
   xout << fmt::format(IndFmt, "");
   for (size_t iCol = 0; iCol < nCols; ++ iCol)
      xout << " " << fmt::format(FmtI, iCol);
   xout << std::endl;
   for (size_t iRow = 0; iRow < nRows; ++ iRow) {
//       xout << fmt::format("{}   ", fmt::format(FmtI, iRow));
      xout << fmt::format("{}   ", fmt::format(IndFmt, iRow));
      IrPrintRow(xout, &pData[iRow*iRowSt], iColSt, nCols, FmtF);
   }
   xout << std::endl;
}


void operator << (std::ostream &xout, FRawShell const &rs)
{
   xout << fmt::format("IrRawShell: l = {}  cen = [{:+8.3f},{:+8.3f},{:+8.3f}]  nCo = {}  nExp = {}\n",
      rs.l, rs.vCen[0], rs.vCen[1], rs.vCen[2], rs.nCo, rs.nExp);
   xout << ("      exps:"); IrPrintRow(xout, rs.pExp, 1, rs.nExp, "{:12.4e}");
   if (rs.pRange) {
      xout << ("    rg/exp:"); IrPrintRow(xout, rs.pRange+1, 1, rs.nExp, "{:12.4f}");
   }
   for (size_t iCo = 0; iCo < rs.nCo; ++ iCo) {
      xout << fmt::format("     {:>2}-co:", iCo);
      IrPrintRow(xout, &rs.pCo[rs.nExp*iCo], rs.nExp, rs.nExp, "{:12.6f}", "");
      if (rs.pRange)
         xout << fmt::format(" | rg/co: {:12.4f}", rs.CoRange(iCo));
      else
         xout << "\n";
   }
   xout.flush();
}

#endif // IR_INCLUDE_PRINT_FUNCTIONS


} // namespace ir
