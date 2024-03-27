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
#include <math.h>
#include <stdexcept>
#include "Ir.h"
#include "CtDft.h"
#include "CtDftFunc.h"
#include "CxAlgebra.h"
#include "CxFortranInt.h"
#include "xc_meta.h"

#include <iostream> // FIXME: remove this
#include "format.h" // FIXME: remove this

namespace ct {


#define GRID_OBTAIN_GRADWT_THRD FORT_Extern(grid_obtain_gradwt_thrd,GRID_OBTAIN_GRADWT_THRD)
void GRID_OBTAIN_GRADWT_THRD(FINTARG icatom, FINTARG npt, double const (*pXyz)[3], double const *pGridWt, FINTARG jwt, double *pWtGrd, double *buf)
{
   assert_rt(0); // not implemented here.
   IR_SUPPRESS_UNUSED_WARNING4(icatom,npt,jwt,buf);
   IR_SUPPRESS_UNUSED_WARNING3(pXyz,pGridWt,pWtGrd);
}

// calculate buffer size for grid_obtain_gradwt_thrd; defined in dftgrid.f
#define GRID_OBTAIN_GRADWT_BUFSIZE FORT_Extern(grid_obtain_gradwt_bufsize,GRID_OBTAIN_GRADWT_BUFSIZE)
void GRID_OBTAIN_GRADWT_BUFSIZE(FORTINT *nbuf, FINTARG icatom, FINTARG npt, FINTARG jwt)
{
   assert_rt(0); // not implemented here.
   IR_SUPPRESS_UNUSED_WARNING4(nbuf,icatom,npt,jwt);
}

static double sqr(double d) { return d*d; }


void DFTFUN_CXX(FXcFunctional const *pXcFn,
         FINTARG fderiv, FINTARG open, FINTARG igrad, FINTARG nGridPt_, double const *wt,
         // density inputs.
         double const *rhoc, double const *rhoo, double const *sigmacc, double const *sigmaco, double const *sigmaoo, double const *tauc, double const *tauo, double const *upsilonc, double const *upsilono,
         // density functional outputs (summed over functionals)
         double *zk, double *vrhoc, double *vrhoo, double *vsigmacc, double *vsigmaco, double *vsigmaoo, double *vtauc, double *vtauo, double *vupsilonc, double *vupsilono,
         // energy outputs (split according to functions; one output per ndftu)
         double *dfu_energies,
         // temporary buffer (currently: len npt).
         double *tmp,
         // flags indicating whether tau/upsilon arrays were passed (to control our rather hacky way of scaling dftfun inputs..)
         FINTARG iHaveTau, FINTARG iHaveUpsilon)
{
   size_t
      nGridPt = nGridPt_,
      nGridPtsKept = 0; // number of points we actually forwarded to the functional
   double
      *pZeros = 0;
   if (open == 0) {
      // pure closed shell case -- code is written in such a way that these
      // quantities are never allocated in this case.
      pZeros = tmp;
      memset(pZeros, 0, sizeof(pZeros[0]) * nGridPt);
      rhoo = pZeros;
      sigmaco = pZeros;
      sigmaoo = pZeros;
      tauo = pZeros;
      upsilono = pZeros;
   }
//    std::cout << fmt::format("     DFTFUN: fderiv = {}  igrad = {}  open = {}  rhoo = {}  vrhoo = {}", fderiv, igrad, open, (void*)rhoo, (void*)vrhoo) << std::endl;

//    FXcFunctional
//       *pXcFn = g_pXcFn.get();

   // copy input into linearized format expected by FXcFunctional.
   double
      *pIn,
      *pOut,
      // indices of points we keep (double instead of int since that's all we have)
      // indices 0 .. nGridPtsKept
      *pGridPtWeKept;
//    if ( iHaveUpsilon )
//       throw std::runtime_error("Upsilon-style Meta-GGAs not supported."  " DFTFUN_CXX_PROXY");
   uint
      nComp = 2;
   bool
      NeedSigma = pXcFn->NeedSigma(),
      NeedTau = pXcFn->NeedTau(),
      NeedUpsilon = pXcFn->NeedUpsilon(),
      TranslateUpsilonToTau = pXcFn->TranslateUpsilonToTau();
   assert(!NeedTau || (bool)iHaveTau);
   assert(!NeedUpsilon || (bool)iHaveUpsilon);
   assert(NeedSigma || !NeedTau); // tau implies sigma.
   assert(!TranslateUpsilonToTau || !NeedTau); // tau cannot be a legitimate input IF we make the tau ourselves here.

//    memset(zk, 0, sizeof(*zk)*nGridPt);
//    if (fderiv) {
//       memset(vrhoc, 0, sizeof(*vrhoc)*nGridPt);
//       memset(vrhoo, 0, sizeof(*vrhoo)*nGridPt);
//       if (NeedSigma) {
//          memset(vsigmacc, 0, sizeof(*vsigmacc)*nGridPt);
//          memset(vsigmaco, 0, sizeof(*vsigmaco)*nGridPt);
//          memset(vsigmaoo, 0, sizeof(*vsigmaoo)*nGridPt);
//       }
//       if (NeedTau) {
//          memset(vtauc, 0, sizeof(*vtauc)*nGridPt);
//          memset(vtauo, 0, sizeof(*vtauo)*nGridPt);
//       }
//       if (NeedUpsilon) {
//          memset(vupsilonc, 0, sizeof(*vupsilonc)*nGridPt);
//          memset(vupsilono, 0, sizeof(*vupsilono)*nGridPt);
//       }
//    }

   if (TranslateUpsilonToTau)
      NeedTau = true;
   if (NeedSigma)
      nComp = 5;
   if (NeedTau)
      nComp = 7;
   if (NeedUpsilon)
      nComp = 9;

   if ((igrad == 0 && pXcFn->NeedSigma()) ||
       (igrad == 1 && !pXcFn->NeedSigma()) ||
       (igrad == 2 && !pXcFn->NeedUpsilon()) ||
       (NeedUpsilon && !((bool)iHaveUpsilon || TranslateUpsilonToTau)) ||
       (igrad >= 3))
      throw std::runtime_error("internal error in hacky DFT functional management.");

   size_t
      nInSt = nComp,
      nOutSt = nInSt + 1; // one for the energies, rest for derivatives wrt input
//    Mem.Alloc(pIn, nInSt*nGridPt);
//    Mem.Alloc(pOut, nOutSt*nGridPt);
   if ( fderiv == 0 )
      nOutSt = 1; // evaluate energy only, not potentials.
   pIn = tmp + nGridPt; // first one is for zeros in closed-shell case.
   pOut = pIn + nInSt * nGridPt;
   pGridPtWeKept = pOut + nOutSt * nGridPt;
   // FIXME: maybe just skip points with zero densities?
   // and map rhoo to [-rhoc,+rhoc]
   memset(pIn, 0, sizeof(pIn[0])*nGridPt*nInSt);
   enum {
      COMP_RhoC = 0,
      COMP_RhoO = 1,
      COMP_SigmaCC = 2,
      COMP_SigmaCO = 3,
      COMP_SigmaOO = 4,
      COMP_TauC = 5,
      COMP_TauO = 6,
      COMP_UpsilonC = 7,
      COMP_UpsilonO = 8
   };
   for (size_t iPtIn = 0; iPtIn < nGridPt; ++ iPtIn) {
      size_t iPtOut = size_t(nGridPtsKept);
      // patch up the densities to fulfull some conditions: |rhoo| must be <= rhoc,
      // and we can't have exactly zero densities or exactly zero or one spin polarization,
      // because otherwise some DFT functionals might blow up.
      double
         *pPtDen = &pIn[iPtOut*nInSt];
      double min_rho = 1e-15;
//       double min_rho = 1e-12;
//       double min_tau = min_rho;
      double min_sigma = 1e-10*sqr(min_rho);
//       double min_sigma = 0.;
      double min_tau = 0.;
      double msp = 1. - 1e-8; // max spin polarization.
      {
         double rhoc_ = rhoc[iPtIn];
         double rhoo_ = rhoo[iPtIn];
         if (rhoc_ < min_rho) // totally skip very low density values.
            continue;
         if (rhoc_ < min_rho) rhoc_ = min_rho;
         if (rhoo_ > rhoc_*msp) rhoo_ = rhoc_*msp;
         if (-rhoo_ > rhoc_*msp) rhoo_ = -rhoc_*msp;

         pPtDen[COMP_RhoC] = rhoc_;
         pPtDen[COMP_RhoO] = rhoo_;
      }
      if (NeedSigma) {
         double sigmacc_ = sigmacc[iPtIn];
         double sigmaco_ = sigmaco[iPtIn];
         double sigmaoo_ = sigmaoo[iPtIn];
         // enforce norm relation ships (all sigma >= 0) and triangle inequalities.
         //    simgaaa = (sigmacc + 2*sigmaco + sigmaoo) >= 0
         //    simgabb = (sigmacc - 2*sigmaco + sigmaoo) >= 0
         if (sigmacc_ < min_sigma) sigmacc_ = min_sigma;
         if (sigmaoo_ < min_sigma) sigmaoo_ = min_sigma;
         if (sigmaco_ > msp*.5*(sigmacc_ + sigmaoo_)) sigmaco_ = msp*.5*(sigmacc_ + sigmaoo_);
         if (-sigmaco_ > msp*.5*(sigmacc_ + sigmaoo_)) sigmaco_ = -(msp*.5*(sigmacc_ + sigmaoo_));
         // ^- should msp be squared here?
         pPtDen[COMP_SigmaCC] = sigmacc_;
         pPtDen[COMP_SigmaCO] = sigmaco_;
         pPtDen[COMP_SigmaOO] = sigmaoo_;
      }
      if (NeedUpsilon) {
         // I don't think there are any special constraints on these, are there?
         pPtDen[COMP_UpsilonC] = upsilonc[iPtIn];
         pPtDen[COMP_UpsilonO] = upsilono[iPtIn];
      }

      if (NeedTau) {
         double tauc_, tauo_;
         if (TranslateUpsilonToTau) {
            // atm: just compute them two times (here and in the derivative translation)...
            // (we could store it easily)
            double TauA_ = 0., dTauA_[12] = {0}, TauB_ = 0, dTauB_[12] = {0};
            rx::xc_ek_pc07_a(TauA_, &dTauA_[0], pPtDen, 1.);
            rx::xc_ek_pc07_b(TauB_, &dTauB_[0], pPtDen, 1.);
            tauc_ = TauA_ + TauB_;
            tauo_ = TauA_ - TauB_;
            assert_rt(TauA_ >= 0. && TauB_ >= 0.);
         } else {
            tauc_ = tauc[iPtIn];
            tauo_ = tauo[iPtIn];
         }

         if (tauc_ < min_tau) tauc_ = min_tau;
         if (tauo_ > msp*tauc_) tauo_ = msp*tauc_;
         if (-tauo_ > msp*tauc_) tauo_ = -(msp*tauc_);
         // taua = tauc + tauo, taub = tauc - tauo ... both also have to be non-negative.
         pPtDen[COMP_TauC] = tauc_; // must be non-negative.
         pPtDen[COMP_TauO] = tauo_;
      }
      pGridPtWeKept[iPtOut] = iPtIn;
      nGridPtsKept += 1;
   }

   memset(pOut, 0, sizeof(pOut[0])*nGridPtsKept*nOutSt);
   pXcFn->Eval(pOut, nOutSt, pIn, nInSt, nGridPtsKept, fderiv);

   // copy back outputs into target format.
   size_t iNextPtOut = 0;
   for (size_t iPtIn = 0; iPtIn < nGridPtsKept; ++ iPtIn) {
      size_t
         iPtOut = size_t(pGridPtWeKept[iPtIn]);
      for (size_t iPtSkipped = iNextPtOut; iPtSkipped < iPtOut; ++ iPtSkipped) {
         // remove energy and derivative contributions of points we left out completely.
         zk[iPtSkipped] = 0.;
         vrhoc[iPtSkipped] = 0.;
         vrhoo[iPtSkipped] = 0.;
         if (NeedSigma) {
            vsigmacc[iPtSkipped] = 0.;
            vsigmaco[iPtSkipped] = 0.;
            vsigmaoo[iPtSkipped] = 0.;
         }
         if (NeedTau && !TranslateUpsilonToTau) {
            vtauc[iPtSkipped] = 0.;
            vtauo[iPtSkipped] = 0.;
         }
         if (NeedUpsilon) {
            vupsilonc[iPtSkipped] = 0.;
            vupsilono[iPtSkipped] = 0.;
         }
      }

      double
         *pPtE = &pOut[iPtIn*nOutSt], // energy values
         *pPtV = pPtE + 1; // derivatives
      double
         *pPtDen = &pIn[iPtIn*nInSt];
      zk[iPtOut] = pPtE[0];
      if (fderiv != 0) {
         if (TranslateUpsilonToTau) {
            // get dtau_alpha/d[rhoc,rhoo,...,upsilono] & dtau_beta/d[rhoc,rhoo,...,upsilono]
            double TauA_ = 0., dTauA_[12] = {0}, TauB_ = 0, dTauB_[12] = {0};
            rx::xc_ek_pc07_a(TauA_, &dTauA_[0], pPtDen, 1.);
            rx::xc_ek_pc07_b(TauB_, &dTauB_[0], pPtDen, 1.);
            // clear out original initial upsilon derivatives (should be 0 already)
            // will be re-made from tau's and other stuff.
            pPtV[COMP_UpsilonC] = 0.;
            pPtV[COMP_UpsilonO] = 0.;
            double
               dEdTauC_ = pPtV[COMP_TauC],
               dEdTauO_ = pPtV[COMP_TauO];
            for (size_t iComp = 0; iComp != nComp; ++ iComp)
               pPtV[iComp] += dEdTauC_ * (dTauA_[iComp] + dTauB_[iComp]) +
                              dEdTauO_ * (dTauA_[iComp] - dTauB_[iComp]);
            // clear out original tau derivatives.
            pPtV[COMP_TauC] = 0.;
            pPtV[COMP_TauO] = 0.;
         }

         vrhoc[iPtOut] = pPtV[COMP_RhoC];
         vrhoo[iPtOut] = pPtV[COMP_RhoO];
         if (NeedSigma) {
            vsigmacc[iPtOut] = pPtV[COMP_SigmaCC];
            vsigmaco[iPtOut] = pPtV[COMP_SigmaCO];
            vsigmaoo[iPtOut] = pPtV[COMP_SigmaOO];
         }
         if (NeedTau && !TranslateUpsilonToTau) {
            vtauc[iPtOut] = pPtV[COMP_TauC];
            vtauo[iPtOut] = pPtV[COMP_TauO];
         }
         if (NeedUpsilon) {
            vupsilonc[iPtOut] = pPtV[COMP_UpsilonC];
            vupsilono[iPtOut] = pPtV[COMP_UpsilonO];
         }
      }
      iNextPtOut = iPtOut + 1;
   }
//    for (size_t i = 0; i < nGridPt; ++ i) {
//       double rho = pIn[i*nInSt];
//       double h = rho * 0.0001;
//       double dout[3] = {0};
//       double di[2] = {0};
//       di[0] = rho + h;
//       pXcFn->Eval(&dout[1], 1, &di[0], 2, 1, 0);
//       di[0] = rho - h;
//       pXcFn->Eval(&dout[0], 1, &di[0], 2, 1, 0);
//       vrhoc[i] = (dout[1] - dout[0])/(2*h);
//       vrhoo[i] = 0;
//       if (isnan(vrhoc[i]))
//          __asm__("int $0x03");
//    }

   dfu_energies[0] += Dot(wt, zk, nGridPt);
//    Mem.Free(pIn);
   IR_SUPPRESS_UNUSED_WARNING(iHaveTau);
//    IR_SUPPRESS_UNUSED_WARNING4(upsilonc,upsilono,vupsilonc,vupsilono);
}




} // namespace ct
