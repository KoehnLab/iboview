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

#include <iostream>
#include <cmath>
#include <stdlib.h> // for atoi

#include "CxDefs.h"
#include "CxTypes.h"
#include "CxOpenMpProxy.h"
#include "CxParse1.h"

#include "CxPhysicalUnits.h"
#include "CtAtomDensity.h"
#include "CxTiming.h"
#include "CxAtomData.h" // for GetCovalentRadius(int iElement)

#include "CtDftGrid_QuadCriteria.h"

using ct_free_atom_density_data::EvalFreeAtomDensity;
using ct_free_atom_density_data::EvalFreeAtomNumElec;

// #define NO_ADAPTIVE_INTEGRATION

namespace mig {

template<class T> inline T sqr(T const &x) { return x*x; }

FQuadCriteria::FQuadCriteria()
{
   // initialize arrays of density exponents and weights for the contributions
   // of to the quadrature residuals
   for (size_t i = 1; i < N; ++ i) {
      double fCurExp = 4./6. + double(i-1)/6.;
      // this is 0 for exp = 7/6 (halfway between rho and ex) and 1 for exp = 1 (rho) and exp = 8/6 (ex).
      double dexp = 6.*(fCurExp - 7./6.);
      // this is 2 for rho^(7/6), 1 for rho^1 (density) and rho^(4/3) (exchange energy), and then quickly
      // falls off for exponents around that.
//       double wexp = 2./(1. + sqr(sqr(sqr(dexp))));
      double wexp;
      if (1) {
         wexp = 2./(1. + sqr(sqr(dexp)));
   //       double wexp = 2./(1. + sqr(dexp));

         // increase the weight of heigher density powers. This increases the representation of the energy.
         wexp *= std::max(1., sqr(fCurExp));
         wexp *= std::max(1., sqr(fCurExp));
         wexp *= std::max(1., sqr(fCurExp));
//       wexp *= std::max(1., sqr(sqr(fCurExp)));
      } else {
         // I think that's what I had until 92204443876fcf6e26b57be42fa650f8fa0b9c87
         // (With N=7 and starting at /*5/6*/)
         wexp = 0.4 * sqr(sqr(fCurExp));
      }

      m_rexps[i] = fCurExp;
      m_wexps[i] = wexp;
//       m_wexps[i] = wexp * std::max(1., sqr(sqr(fCurExp)));
   }
   m_rexps[0] = 0;
   m_wexps[0] = 1.; // this is for the correlation energy.

//          fValues[i] += w * rhop * sqr(sqr(fCurExp));
//          fValues[i] += w * rhop * sqr(sqr(fCurExp)) * .5;
//          fValues[i] += w * rhop * sqr(fCurExp);
//          fValues[i] += w * rhop * 1.*std::exp(-1.*sqr(8./6. - fCurExp));
//          fValues[i] += w * rhop * .5;

}



double FQuadCriteria::fResidual(FQuadCriteria const &Ref) const
{
   double
      r = 0.;
   for (size_t i = 0; i < N; ++ i) {
      r += sqr(fValues[i] - Ref.fValues[i]);
   }
   return std::sqrt(r/double(N));
//    return std::sqrt(r)/double(N);
}

void FQuadCriteria::operator += (FQuadCriteria const &other)
{
   for (size_t i = 0; i < N; ++ i)
      fValues[i] += other.fValues[i];
}

void FQuadCriteria::operator *= (double f)
{
   for (size_t i = 0; i < N; ++ i)
      fValues[i] *= f;
}





static void pw91_lda_xc(double &Ec, double &Ex, double &rho16, double rho) {
   // see pw91_lsdac.max. This here is Eu (the unpolarized version).
   double const
      A = 0.0310907,
      Alpha1 = 0.21370,
      Beta1 = 7.5957,
      Beta2 = 3.5876,
      Beta3 = 1.6382,
      Beta4 = 0.49294,
//       cc = std::pow(3./(4.*M_PI), (1./3.)),    // prefactor for Wigner-Seitz radius rs
//       cf = std::pow(-.75*(3./M_PI), (1./3.));  // prefactor of Dirac exchange energy (should have minus sign)
      cc = 0.62035049089940009,
      cf = -0.73855876638202234;

   // Direac exchange
   rho16 = std::pow(rho, (1./6.));
   double
      rho13 = rho16*rho16;
   Ex = cf * sqr(sqr(rho13));

   Ec = 0.;
   if (rho13 > 1e-16 ) {
      // PW91 LDA corr (unpolarized). Note: Ecu is the correlation energy *per
      // electron*, not per volume. We thus still need to multiply it by the
      // density.
      double
         Rs = cc/rho13,
         rs12 = std::sqrt(Rs),
         Ecu = -2.*A*(1. + Alpha1*Rs) * std::log(1. + 1./(2.*A*rs12*(Beta1 + rs12*(Beta2 + rs12*(Beta3 + rs12*Beta4)))));
         // ^- 10.1103/PhysRevB.45.13244, eq. 10 and Table I (row for Ec(rs,0), 4th column)
      Ec = Ecu * rho;
   }
}


void FQuadCriteria::Eval(double *ri, double *wi, size_t n, double fDensityShift, int iElement, ct::FMemoryStack &Mem)
{

   ct::TMemoryLock<double>
      LockD(n, &Mem);
   double
      *di = LockD.p;
   EvalFreeAtomDensity(di, ri, n, iElement, 0);
   return Eval(di, wi, ri, n, iElement, Mem);
   IR_SUPPRESS_UNUSED_WARNING(fDensityShift); // <-- don't remember what that was for...
}


void FQuadCriteria::Eval(double *di, double *wi, double *ri, size_t n, int iElement, ct::FMemoryStack &Mem)
{
   for (size_t i = 0; i < N; ++ i)
      fValues[i] = 0;

   for (size_t i = 0; i < n; ++ i) {
      double
         w = wi[i],
         rho = di[i],
         Ec, Ex, rho16;
      pw91_lda_xc(Ec, Ex, rho16, rho);

      fValues[0] += 100. * w * Ec;

#if 0
      assert(m_rexps[1] == 5./6.);
      double
         rho56 = rho16*rho16*rho16*rho16*rho16, // this one also works for rho16 == 0 ...unlike rho/rhop.
         rhop = rho56;
#else
      assert(m_rexps[1] == 4./6.);
      double
         rho46 = rho16*rho16*rho16*rho16,
         rhop = rho46;
#endif

      // absorb radial weight into accumulator of density power.
      rhop *= w;

      // large powers are good for representing the energy accurately,
      // small powers help representing the density.
      for (size_t i = 1; i < N; ++ i) {
         fValues[i] += m_wexps[i] * rhop;
         rhop *= rho16;
      }
   }

   IR_SUPPRESS_UNUSED_WARNING(Mem);
   IR_SUPPRESS_UNUSED_WARNING(ri);
   IR_SUPPRESS_UNUSED_WARNING(iElement); // already used to compute input density.
}



} // namespace ct
