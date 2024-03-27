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

#ifndef CT_RHF_GUESS_H
#define CT_RHF_GUESS_H

#include "CtRhf.h"

namespace ct {

void MakeDensityGuess(FMatrixView Fock, FAtomSet const &Atoms, FMatrixView Scd1, FMatrixView CoreH1,
   FBasisSet *pOrbBasis, FBasisSet *pGuessBasis, FFockComponentBuilderList &FockBuilders, FLog *pLog, FHfOptions const &Options, FMemoryStack &Mem);

} // namespace ct

#endif // CT_RHF_GUESS_H