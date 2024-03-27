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

#ifndef CT_DFTI_H
#define CT_DFTI_H

#include "CxTypes.h"
#include "CxDefs.h"
#include "CtDftGrid.h"
#include "CtDftFunc.h"
#include "CtDft.h"
#include "CxIo.h"
#include "CxTiming.h"
#include "CtBasisSet.h"

namespace dfti {

enum FDftiFlags {
   DFTI_AuxiliaryExpand = 0x01,
   DFTI_MakeXc = 0x02,
   DFTI_MakeGradient = 0x04,
   // if set, include the gradient of the integration grid in the
   // case of a gradient computation
   DFTI_MakeGridGradient = 0x08,
   DFTI_MeasureTime = 0x10,
   // if set, compute the densities (rhoc,grad rhoc,tauc,rhoo,grad rhoo,tauo)
   // for the given input density matrix/orbitals and store them at
   // pBackgroundDensities.
   DFTI_MakeBackgroundDensities = 0x20,
   // if set, add densities from pBackgroundDensities to the total
   // densities before inputting them into the functional. Also, compute
   // the Exc and vxc contributions with weighting (rho/(rho + rho_background))*epsilon.
   DFTI_AddBackgroundDensities = 0x40,
   // subtract \int \rho(\vec r) v_\mathrm{xc}(\vec r)\,\d^3r from density functional energy.
   // (experimental option for allowing inclusion of full band structure term in total energy
   // computation of DFT)
   DFTI_CounterDensityXcIntegral = 0x80
};


enum {
   // [rhoc, d/dx rhoc, d/dy rhoc, d/dz rhoc, (tauc or upsilonc), (same for rhoo)]
   DFTI_nBackgroundDensityComponents = 12
};

struct FDftiArgs
{
   uint32_t
      Flags; // combination of DFTI_*.

   // energies, energy gradients, and triangular fock matrices (output) or Fock potential expansion coefficients.
   double *pDfuEnergies; double *pGradient; double *pFockC; double *pFockO;
   // triangular density matrices (input) in case of regular XC, expansion coefficients
   // in terms of auxiliary basis in other case.
   // if pDenO == 0, open-shell contributions are not evaluated.
   double *pDenC; double *pDenO;
   double *pOccOrbC; double *pOccOrbO;
   size_t nOccC; size_t nOccO;

   // basis over which the orbitals are expanded (regular) or auxiliary basis.
   ct::FRawBasis *pBasis;

   // the functional.
   ct::FXcFunctional *pXcFn;

   // molecular integration grid.
   mig::FDftGrid const *pDftGrid;

   double ThrOrb; double LogThrOrb; // threshold for orbital-on-grid screening

   double *pTimings;
   double *pfElecTotal;

   // a set of however-many-density-components-we-have which may be used in
   // embedded DFT calculations. Must hold room for
   // DFTI_nBackgroundDensityComponents*nGridPt scalars.
   double mutable *pBackgroundDensities;
public:
   bool OpenShell() const { return pDenO != 0; };
   bool MakeGradient() const { return Flags & DFTI_MakeGradient; };
   bool MakeXc() const { return Flags & DFTI_MakeXc; };
   bool MakeGridGradient() const { return Flags & DFTI_MakeGridGradient; }
   bool UseAuxiliaryExpansion() const { return Flags & DFTI_AuxiliaryExpand; };
   bool MeasureTime() const { return Flags & DFTI_MeasureTime; }
   bool NeedTau() const { return pXcFn->NeedTau(); }
   bool NeedSigma() const { return pXcFn->NeedSigma(); }
   bool NeedUpsilon() const { return pXcFn->NeedUpsilon(); }
   bool MakeBackgroundDensities() const { return Flags & DFTI_MakeBackgroundDensities; }
   bool AddBackgroundDensities() const { return Flags & DFTI_AddBackgroundDensities; }
   size_t nBf() const { return pBasis->nFn(); }; // number of basis functions in pBasis
   size_t nCenters() const { return pBasis->nCen(); }; // number of centers (required to determine gradient dimension)
   size_t nFockSize() const {
      size_t nBf_ = nBf();
      if (UseAuxiliaryExpansion())
         return nBf_;
      else
         return (nBf_*(nBf_+1))/2;
   }
};

void AccXc(FDftiArgs &Args, ct::FLog &Log, ct::FTimerSet *pTimers, ct::FMemoryStack &Mem_);


enum FEvalBfDerivType {
   // evaluate values of the actual basis functions mu(r)
   // (i.e., no derivative components)
   DERIVCOMP_ValueOnly = 0,
   // evaluate [mu(r), d/dx mu(r), d/dy mu(r), d/dz mu(r)]
   DERIVCOMP_Value_dFirst = 1,
   // evaluate [mu(r), d/dx mu(r), d/dy mu(r), d/dz mu(r), d^2/d[xx] mu(r), ..., d^2/d[yz] mu(r)]
   // (i.e., value of bf and all first and 2nd cartesian derivatives)
   DERIVCOMP_Value_dFirst_dSecond = 2,
   // evaluate [mu(r), Laplace mu(r)]
   // (i.e., the values of the basis functions, and their scalar Laplace derivatives)
   DERIVCOMP_Value_Laplace = 3,
   // evaluate [mu(r), d/dx mu(r), d/dy mu(r), d/dz mu(r), Laplace mu(r)]
   // (i.e., value of bf and all its first cartesian derivatives, and the the bf's scalar Laplace derivatives)
   DERIVCOMP_Value_dFirst_Laplace = 4,
   // evaluate [mu(r), d/dx mu(r), ..., d^2/d[yz] mu(r), Laplace mu(r), d/dx Laplace mu(r), d/dy Laplace mu(r), d/dz Laplace mu(r)]
   // (i.e., evaluate mu(r), all its first and 2nd cartesian derivatives, Laplace mu(r), and all of Laplace mu(r)'s first Cartesian derivatives.)
   DERIVCOMP_Value_dFirst_dSecond_Laplace_dFirst = 5
};
// ^- WARNING: if changing this, adjust GetNumBfDerivComps()!

size_t GetNumBfDerivComps(FEvalBfDerivType DerivType, size_t *iLaplaceComp0 = 0);
FEvalBfDerivType MakeBfDerivType(unsigned nDiff, bool MakeLaplaceToo);

void IrEvalBfn(double *pOrbValAo, size_t *pMap, size_t &nMap, double const *pGridPt, size_t iStGridPt, size_t nGridPt, FEvalBfDerivType DerivType, ct::FRawBasis const *pBasis, double ThrOrb, double LogThrOrb, ct::FMemoryStack &Mem);


} // namespace dfti


#endif // CT_DFTI_H
