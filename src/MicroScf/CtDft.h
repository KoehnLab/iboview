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

#ifndef CT_DFT_H
#define CT_DFT_H

// #include "CtDftFunc.h"
// #include "CtDftGrid.h"

// namespace ct {
//    extern FXcFunctionalPtr
//       g_pXcFunctional;
// }
namespace ct {
   struct FXcFunctional;
}


#include "CxFortranInt.h"

namespace ct {

extern "C" {
   // proxy for dftfun (evaluation of density functional (+gradients))
   #define DFTFUN_CXX FORT_Extern(dftfun_cxx,DFTFUN_CXX)
   void DFTFUN_CXX(FXcFunctional const *pXcFn,
            FINTARG fderiv, FINTARG open, FINTARG igrad, FINTARG npt, double const *wt,
            // density inputs.
            double const *rhoc, double const *rhoo, double const *sigmacc, double const *sigmaco, double const *sigmaoo, double const *tauc, double const *tauo, double const *upsilonc, double const *upsilono,
            // density functional outputs (summed over functionals)
            double *zk, double *vrhoc, double *vrhoo, double *vsigmacc, double *vsigmaco, double *vsigmaoo, double *vtauc, double *vtauo, double *vupsilonc, double *vupsilono,
            // energy outputs (split according to functions; one output per ndftu)
            double *dfu_energies,
            // temporary buffer (currently: len 20*npt).
            double *tmp,
            // flags indicating whether tau/upsilon arrays were passed (to control our rather hacky way of scaling dftfun inputs..)
            FINTARG iHaveTau, FINTARG iHaveUpsilon);

   // calculate derivatives of grid weights with respect to atomic coordinates of atom icatom
   #define GRID_OBTAIN_GRADWT_THRD FORT_Extern(grid_obtain_gradwt_thrd,GRID_OBTAIN_GRADWT_THRD)
   void GRID_OBTAIN_GRADWT_THRD(FINTARG icatom, FINTARG npt, double const (*pXyz)[3], double const *pGridWt, FINTARG jwt, double *pWtGrd, double *buf);

   // calculate buffer size for grid_obtain_gradwt_thrd
   #define GRID_OBTAIN_GRADWT_BUFSIZE FORT_Extern(grid_obtain_gradwt_bufsize,GRID_OBTAIN_GRADWT_BUFSIZE)
   void GRID_OBTAIN_GRADWT_BUFSIZE(FORTINT *nbuf, FINTARG icatom, FINTARG npt, FINTARG jwt);
}

}


#endif // CT_DFT_H
