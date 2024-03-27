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

#include <cstddef> // for size_t
#include <math.h>
#include "coords_meta.h"

namespace rx {
   inline double sqr(double x) { return x*x; }
   #define FORT_Extern(lower,UPPER) lower##_

   #include "pc_pos_x.inl"
   #include "pc_pos_x_hess.inl"
   #include "pc_pos_y.inl"
   #include "pc_pos_y_hess.inl"
   #include "pc_pos_z.inl"
   #include "pc_pos_z_hess.inl"
   #include "pc_stretch.inl"
   #include "pc_stretch_hess.inl"
   #include "pc_bend.inl"
   #include "pc_bend_hess.inl"
   #include "pc_torsion.inl"
   #include "pc_torsion_hess.inl"
   #include "pc_dist_cop21.inl"
   #include "pc_dist_cop21_hess.inl"
   #include "pc_dist_cop22.inl"
   #include "pc_dist_cop22_hess.inl"
   #include "eh_stretch.inl"
   #include "eh_stretch_hess.inl"
   #include "eh_bend.inl"
   #include "eh_bend_hess.inl"
   #include "eh_torsion.inl"
   #include "eh_torsion_hess.inl"
   #include "ef_d3bj.inl"
   #include "ef_d3bj_hess.inl"
   #include "ef_coul_gauss.inl"
   #include "ef_coul_gauss_hess.inl"
} // namespace rx
