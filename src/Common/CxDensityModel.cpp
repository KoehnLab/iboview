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

#include "CxDensityModel.h"
#include "CtAtomDensity.h"
#include "CxAtomData.h" // for ElementNameFromNumber
#include "CxVec3.h"


using ct_free_atom_density_data::EvalFreeAtomDensity;
using ct_free_atom_density_data::EvalFreeAtomNumElec;

namespace ct {

FDensityModel::~FDensityModel()
{
}


FDensityModel_FreeAtomSuperposition::FDensityModel_FreeAtomSuperposition(FRawAtomList const &Atoms, int iFreeAtomFitType, bool OnlyAtomsOfInterest)
   : m_Atoms(Atoms), m_iFreeAtomFitType(iFreeAtomFitType), m_AddOnlySelectedAtoms(OnlyAtomsOfInterest)
{
}

FDensityModel_FreeAtomSuperposition::~FDensityModel_FreeAtomSuperposition()
{
}


typedef ct::TVector3<double>
   FVec3d;
static double EvalFreeAtomDensityAtPoint(FVec3d const &vGridPt, FVec3d const &vAtPos, int iElement, int iFreeAtomFitType)
{
   double fDensityAt(0);
   double rAt = double(ct::Dist(vAtPos, vGridPt));
   EvalFreeAtomDensity(&fDensityAt, &rAt, 1,iElement, iFreeAtomFitType);
   return fDensityAt;
}


void FDensityModel_FreeAtomSuperposition::ComputeDensity(double *pOut, double const *pGridPts, size_t stGridPt, size_t nGridPt, FMemoryStack &Mem, int const *pAtomsOfInterest, size_t nAtomsOfInterest) const
{
   for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt) {
      double
         fDensity(0);
      double const
         *pGridPt = &pGridPts[iGridPt * stGridPt];
      FVec3d
         vGridPt(pGridPt[0], pGridPt[1], pGridPt[2]);
      if (!m_AddOnlySelectedAtoms) {
//       if (false) {
         for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt)
            fDensity += EvalFreeAtomDensityAtPoint(vGridPt, m_Atoms[iAt].vPos, m_Atoms[iAt].iElement, m_iFreeAtomFitType);
      } else {
         for (size_t iiAt = 0; iiAt < nAtomsOfInterest; ++ iiAt) {
            int iAt = pAtomsOfInterest[iiAt];
            if (size_t(iAt) >= m_Atoms.size())
               throw std::runtime_error("invalid atom index in FDensityModel_FreeAtomSuperposition::ComputeDensity (0-based!)");
            fDensity += EvalFreeAtomDensityAtPoint(vGridPt, m_Atoms[iAt].vPos, m_Atoms[iAt].iElement, m_iFreeAtomFitType);
         }
      }
//       for (size_t iAt = 0; iAt < m_Atoms.size(); ++ iAt) {
//          // TODO: add options for skipping all contributions except for the ones from the "atoms of interest"?
//          double fDensityAt(0);
//          double rAt = double(ct::Dist(m_Atoms[iAt].vPos, vGridPt));
//          EvalFreeAtomDensity(&fDensityAt, &rAt, 1, m_Atoms[iAt].iElement, m_iFreeAtomFitType);
//          fDensity += fDensityAt;
//       }
      pOut[iGridPt] = fDensity;
   }
//    IR_SUPPRESS_UNUSED_WARNING(Mem);
   (void)Mem;
}


} // namespace ct
