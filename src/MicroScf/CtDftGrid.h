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

#ifndef CT_DFTGRID_H
#define CT_DFTGRID_H

#include <vector>
#include <string>

#include "CtDftGrid_CommonDefs.h"
#include "CtDftGrid_Params.h"

#include "CxPodArray.h"
#include "CxIo.h" // for output log (FLog)
#include "CxRawAtom.h" // for specifying the input geometry


namespace mig { // molecular integration grids

/// A non-uniform three-dimensional molecular integration grid, based on atom
/// positions, suitable for density functional and related real-space
/// integration purposes.
struct FDftGrid : public ct::FIntrusivePtrDest
{
   struct FPoint{
      FVector3
         vPos;
      FScalar
         fWeight;
      ptrdiff_t
         iAtomicCenter;
   };

   typedef std::vector<FPoint>
      FPointList;
   FPointList
      // \sum_p p.fWeight integrates to total volume.
      Points;

   struct FGridBlock {
      size_t
         iFirst, iLast; // grid block spans points [iFirst,iLast).
      FVector3
         // average of spanned points
         vCenter;
      FScalar
         // radius around vCenter such that all points spanned by
         // the block lie inside the sphere of fRadius around vCenter.
         fRadius,
         // largest weight of any point within the block.
         fLargestWeight;
      ptrdiff_t
         // atom with which the grid point moves. -1 if not moving with a center.
         iAtomicCenter;
      size_t nPt() const { return iLast - iFirst; }
   };
   typedef std::vector<FGridBlock>
      FGridBlockList;
   FGridBlockList
      GridBlocks;

   FDftGrid(ct::FRawAtomList const &Atoms, FDftGridParams const &Params, ct::FLog *pLog = 0);
   ~FDftGrid();

   // data in compatibility format -- DFTI-CXX expects positions and
   // weights as separate entities.
   ct::TArray<double[3]>
      Positions;
   ct::TArray<double>
      Weights;

   // returns total number of grid points
   size_t size() const { return Points.size(); }
   // Maybe we should compute and return the total covered 3D volume, too?
private:
   // copy stuff from this->Points to this->Positions and this->Weights.
   void MakeAdditionalRepresentations();
};

typedef ct::TIntrusivePtr<FDftGrid>
   FDftGridPtr;
typedef ct::TIntrusivePtr<FDftGrid const>
   FDftGridCptr;


} // namespace mig

#endif // CT_DFTGRID_H
