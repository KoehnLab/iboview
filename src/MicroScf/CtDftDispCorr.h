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
#include <algorithm>

#include "CxTypes.h"
#include "CtAtomSet.h"
#include "CtMatrix.h" // only for Gradient.


namespace ct {

struct FD3DispersionContextImpl;

struct FD3DispersionContext
{
   // note: this allocates memory on Mem
   explicit FD3DispersionContext(FLog *pLog_, FAtomSet const *pAtoms, std::string const &DispCorrType, std::string const &XcName, FMemoryStack &Mem);
   ~FD3DispersionContext();

   // compute and return D3(BJ) dispersion energy for initialized atom set.
   // if Gradient.pData != 0, also compute analytic gradient of D3(BJ)
   // dispersion energy and add it to 'Gradient'. Gradient should be a 3 x
   // nAtoms matrix.
   double CalcEnergyAndGradient(FMatrixView Gradient, FMemoryStack &Mem) const;
   // estimate D3(BJ) contributions to the Hessian. Calculation need not be exact (e.g., C6/C8 frozen at ref geometry).
   void EstimateHessian(FMatrixView Hessian, FMemoryStack &Mem, bool Add) const;
protected:
   FD3DispersionContextImpl
      *p;
};

} // namespace ct
