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

/// @file CtVolumeProperty.h
///
/// Support code for evaluation and analysis of properties which are defined as
/// a function of space (e.g., the electron density).
///
/// Note: Original version of code here comes from various placed:
/// - FDensityEvalContext & co are moved here from IvAnalysis.cpp
///   (from the TFVC fragmentation analysis)
///
/// - MakeGridValues somehow was in ct but in IvIntInit.cpp (???).
///   It seems the file is now empty---i guess this was one of the first
///   IboView files I made.
///
/// TODO:
/// - merge with code in IvIsoSurface for evaluating orbitals and
///   densities on a grid (split off actual evaluation and
///   core storage routines of FVolumeProperty in
///   IvDocument.h/.cpp and IvIsoSurface.cpp)

#ifndef CT_VOLUME_PROPERTIES_H
#define CT_VOLUME_PROPERTIES_H

#include "CtMatrix.h"
#include "CxPodArray.h"
#include "CtBasisSet.h"
// #include "CxIo.h"



namespace ct {


// an object mapping each point in space to the numerical weight (0...1) its represented
// fragment has at this point.
struct FFragmentWeightFn
{
   virtual void EvalWeights(double *pOut, FMatrixView Grid_, FMemoryStack &Mem) const = 0;
};


// TODO: split off into abstract sub-type: one for orbital-based actual densities,
// and one for superposition-of-free-atoms reference densities.
// FIXME: There is a new version of this in CtVoronoiPartition.h/.cpp
// FIXME: There is also very similar code in FVolumeProperty of IvDocument.cpp
struct FDensityEvalContext : public FIntrusivePtrDest
{
   // note: makes a copy of input orbital matrix and occupation numbers.
   explicit FDensityEvalContext(FMatrixView Orbs, double const *pOcc, FBasisSetPtr pOrbBasis);
   ~FDensityEvalContext();

   // compute total electron density at points given by pvPos[0...nPos-1]. Results
   // stored in output array pOut.
   void ComputeDensity(double *pOut, FMatrixView Grid_, FMemoryStack &Mem);

   // note: pWeightFn may be 0. Grid_ comes as a (4,nGridPt) array of grid points, the 4th denoting
   // the regular integration weight (without pWeightFn's factor, which will be multiplied in)
   // Note: Output is INCREMENTED, not set to zero!
   void AccFragmentOrbitalOverlap(FMatrixView Out, FMatrixView Grid_, FFragmentWeightFn *pWeightFn, FMemoryStack &Mem);

protected:
   ct::FHeapMatrix
      m_Orbs;
   TArray<double>
      m_Occ; // total orbital occupations of the orbitals in m_Orbs
   FBasisSetPtr
      m_pOrbBasis;
};
typedef ct::TIntrusivePtr<FDensityEvalContext>
   FDensityEvalContextPtr;
typedef ct::TIntrusivePtr<FDensityEvalContext const>
   FDensityEvalContextCptr;


void MakeGridValues(double *pOut, FMatrixView Grid, FMatrixView Orb, unsigned GridDxOrder, FBasisSet const *pBasis, FMemoryStack &Mem_);

}

#endif // CT_VOLUME_PROPERTIES_H
