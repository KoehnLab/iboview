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

// UPDATE: version with renormalization in /home/cgk/dev/sphere_grids/grid-accuracy-test/src/CtAtomDensity.cpp (2021-02)
#include <stddef.h> // for size_t
#include <cmath>
#include <assert.h>
#include <stdexcept>

namespace ct_free_atom_density_data {

static size_t const
   nMaxExp = 32,
   nMaxL = 6;

struct FFreeAtomDensityShell
{
   size_t
      iType, l, nExp, nCo;
   double const
      // note: contractions are normalized for *raw* primitives (not normalized primitives)
      *pExps, *pCos;
   double const
      // occupation numbers for each contracted function (could technically be absorbed into pCos,
      // but I prefer them separately at the moment). This has the (2*l+1) from the number of
      // functions in the shell absorbed.
      *pOcc;

   double EvalDensity(double r) const;
   double fElec() const;
};

struct FFreeAtomDensityInfo
{
   int
      iElement;
   size_t
      nShells;
   FFreeAtomDensityShell const
      *pGaussShells;

   double EvalDensity(double r) const;
   double fElec() const;
};

#include "CtAtomDensity.inl"

double FFreeAtomDensityShell::fElec() const
{
   double fElecOut = 0.;
   for (size_t iCo = 0; iCo < nCo; ++ iCo) {
      fElecOut += pOcc[iCo];
   }
   return fElecOut;
}

double FFreeAtomDensityShell::EvalDensity(double r) const
{
   assert(nExp < nMaxExp);
   assert(l < nMaxL);

   if (iType == 0) {
      double
         fExpVal[nMaxExp];
      for (size_t iExp = 0; iExp < nExp; ++ iExp) {
         double fArg = -pExps[iExp]*r*r;
         if (fArg > -37.) // exp(-37) ~= 9e-17
            fExpVal[iExp] = std::exp(fArg);
         else
            fExpVal[iExp] = 0.;
      }

      double
         Density = 0.,
         rl = 1.; // r^l.

      for (size_t i = 1; i <= l; ++ i)
         rl *= r;

      for (size_t iCo = 0; iCo < nCo; ++ iCo) {
         // evaluate basis function minus the r^l (added at end since equal for all functions)
         double const
            *pCo = &pCos[iCo * nExp];
         double
            fFn = 0.;
         for (size_t iExp = 0; iExp < nExp; ++ iExp) {
            fFn += pCo[iExp] * fExpVal[iExp];
         }
         // and add to density (note that the basis functions are orthonormal).
         Density += pOcc[iCo] * fFn * fFn;
      }
      Density *= rl * rl;  // r^l is squared because it belongs to the basis function (from Slm(r)), not the density.
      return Density;
   } else {
      assert(iType == 1);
      assert(nCo == 1);
      double
         Density = 0.,
         rSq = r*r,
         rSql = 1.; // r^(2l)
      for (size_t i = 1; i <= l; ++ i)
         rSql *= rSq;
      for (size_t iExp = 0; iExp < nExp; ++ iExp) {
         double
            fArg = -pExps[iExp]*r;
         if (fArg < -37.)
            continue;
         Density += std::exp(fArg) * pCos[iExp];
      }
      Density *= rSql;
      return Density;
   }
}


double FFreeAtomDensityInfo::fElec() const
{
   double
      fElec = 0.;
   for (size_t iShell = 0; iShell < nShells; ++ iShell)
      fElec += pGaussShells[iShell].fElec();
   return fElec;
}


double FFreeAtomDensityInfo::EvalDensity(double r) const
{
   double
      Density = 0.;
   for (size_t iShell = 0; iShell < nShells; ++ iShell)
      Density += pGaussShells[iShell].EvalDensity(r);
   return Density;
}


static FFreeAtomDensityInfo const *GetInfo(int iElement, int iFitType)
{
   if (iElement < 1)
      return 0; // dummy. has no density.
   if (size_t(iElement) > sizeof(s_FreeElementDensitiesG)/sizeof(s_FreeElementDensitiesG[0]))
      throw std::runtime_error("ran out of elements in CalcFreeAtomDensity.");
   FFreeAtomDensityInfo const
      *pInfo;
   if (iFitType == 0)
      pInfo = &s_FreeElementDensitiesG[size_t(iElement-1)];
   else
      pInfo = &s_FreeElementDensitiesS[size_t(iElement-1)];
   assert(pInfo->iElement == iElement);
   return pInfo;
}


void EvalFreeAtomDensity(double *pDensity, double const *pDist, size_t n, int iElement, int iFitType)
{
   FFreeAtomDensityInfo const
      *pInfo = GetInfo(iElement, iFitType);
   if (!pInfo) {
      // dummy atom. Has no density. So clear output.
      for (size_t i = 0; i < n; ++ i)
         pDensity[i] = 0;
   } else {
      for (size_t i = 0; i < n; ++ i)
         pDensity[i] = pInfo->EvalDensity(pDist[i]);
   }
}


double EvalFreeAtomNumElec(int iElement, int iFitType)
{
   FFreeAtomDensityInfo const
      *pInfo = GetInfo(iElement, iFitType);
   return pInfo->fElec();
}


} // namespace ct_free_atom_density_data






#ifdef MOLPRO
#define MACHINES_H_DEFINES_ONLY
#include "util/machines.h"

#define CT_EVAL_FREE_ATOM_DENSITY FORT_Extern(ct_eval_free_atom_density,CT_EVAL_FREE_ATOM_DENSITY)
extern "C" {
   void CT_EVAL_FREE_ATOM_DENSITY(double *pDensity, double *pDist, FORTINT const &nPoints, FORTINT const &iElement, FORTINT const &iFitType) {
      return ct_free_atom_density_data::EvalFreeAtomDensity(pDensity, pDist, size_t(nPoints), int(iElement), int(iFitType));
   }
}

#endif

