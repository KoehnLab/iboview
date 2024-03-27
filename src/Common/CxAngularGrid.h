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

#ifndef CX_ANGULAR_GRID_H
#define CX_ANGULAR_GRID_H

#include <stddef.h>

namespace aig { // angular integration grids
   enum FAngularGridSymmetry{
      SYMMETRY_Icosahedral,
      SYMMETRY_Octahedral,
      SYMMETRY_Tetrahedral
   };

   typedef void (*FGridGeneratorFn)(double *pOut);

   struct FAngularGridEntry {
      size_t
         // maximum L of spherical harmonics still integrated exactly.
         MaxL,
         // number of points in the angular grid
         nPoints;
      FAngularGridSymmetry
         Symmetry;
      double
         // ratio between largest grid weight and smallest grid weight in the grid
         // (ratios close to 1.0 are good.)
         WeightSpread;
      FGridGeneratorFn
         // function pointer to the function which makes the output points of the grid.
         // Functions of this type store an array of [x0, y0, z0, w0, x1, y1, z1, w1, ..., w_{N-1}]
         // at *pOut, where N is this->nPoints.
         pGeneratorFn;
      void MakeGrid(double (*pOut)[4]) const {
         return pGeneratorFn(&pOut[0][0]);
      }
   };

   // return number of tabulated angular grids
   size_t GetNumAngularGrids();
   // return pointer to iGrid'th tabulated angular grid information
   FAngularGridEntry const *GetAngularGridInfo(size_t iGrid);


} // namespace aig

#endif // CX_ANGULAR_GRID_H
