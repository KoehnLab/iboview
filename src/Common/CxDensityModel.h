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

#ifndef CX_DENSITY_MODEL_H
#define CX_DENSITY_MODEL_H

#include "CxIntrusivePtr.h"
#include "CxMemoryStack.h"
#include "CxRawAtom.h"

namespace ct {

// TODO: integrate with CtVolumeProperty // FDensityEvalContext... (which evaluates densities
// for *real* SCF solutions). This here is just supposed to be the corresponding abstract
// base class
struct FDensityModel : public FIntrusivePtrDest1
{
   virtual ~FDensityModel();

   // compute total electron density at points given by pGridPt[ixyz + stGridPt * iGridPt].
   // Results stored in output array pOut. The optional parameter pAtomsOfInterest may be
   // used in a model-specific way to restrict/modify the density computation (e.g., if
   // only density from the given specific atoms are meant to be considered)
   virtual void ComputeDensity(double *pOut, double const *pGridPt, size_t stGridPt, size_t nGridPt, FMemoryStack &Mem, int const *pAtomsOfInterest = 0, size_t nAtomsOfInterest = 0) const = 0;
};
typedef ct::TIntrusivePtr<FDensityModel>
   FDensityModelPtr;
typedef ct::TIntrusivePtr<FDensityModel const>
   FDensityModelCptr;


struct FDensityModel_FreeAtomSuperposition : public FDensityModel
{
   explicit FDensityModel_FreeAtomSuperposition(FRawAtomList const &Atoms, int iFreeAtomFitType = 0, bool OnlyAtomsOfInterest = false);
   ~FDensityModel_FreeAtomSuperposition();

   void ComputeDensity(double *pOut, double const *pGridPts, size_t stGridPt, size_t nGridPt, FMemoryStack &Mem, int const *pAtomsOfInterest = 0, size_t nAtomsOfInterest = 0) const; // override
protected:
   FRawAtomList
      m_Atoms;
   int
      m_iFreeAtomFitType;
   bool
      m_AddOnlySelectedAtoms;
};


} // namespace ct

#endif // #define CX_DENSITY_MODEL_H
