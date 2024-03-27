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

#ifndef CT_IAO_H
#define CT_IAO_H

#include "CtBasisSet.h"
#include "CtMatrix.h"
#include "CtAtomSet.h"

namespace ct {

   enum FIaoBasisFlags {
      IAO_OrthType = 0x03,
      IAO_OrthNone = 0x00,
      IAO_OrthSym = 0x01,
      IAO_OrthZbd = 0x02,
      IAO_NormalizeInput = 0x04
   };

   // void SymOrth(FMatrixView Orbs, FMatrixView const S, FMemoryStack &Mem);
   void MakeIaoBasisNew(FMatrixView CIb, FBasisSet *pMinBasis,
      FAtomSet const &Atoms, FBasisSet const *pOrbBasis,
      FMatrixView const &COcc_, FMatrixView const &S1, FMatrixView const &S1cd, unsigned Flags,
      FMemoryStack &Mem);

   void MakeIaoBasis(FMatrixView &CIb, size_t &nIb, FBasisSetPtr &pMinBasis,
      FAtomSet const &Atoms, FBasisSet const *pOrbBasis, FMatrixView const &COcc, FMemoryStack &Mem, unsigned Flags = IAO_OrthSym|IAO_NormalizeInput);

   void MakeIaoChargesRaw(double *pOut, double fOcc, double const *pIaoCoeffs, FBasisSet *pMinBasis, FAtomSet *pAtoms);


} // namespace ct



#endif // CT_IAO_H
