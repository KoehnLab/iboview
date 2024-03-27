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

#ifndef COORDS_META_H
#define COORDS_META_H

namespace rx {
   typedef std::size_t
      index_t;

   void pc_pos_x(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_pos_x_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_pos_y(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_pos_y_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_pos_z(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_pos_z_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_stretch(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_stretch_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_bend(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_bend_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_torsion(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_torsion_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_dist_cop21(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_dist_cop21_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_dist_cop22(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom);
   void pc_dist_cop22_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom);
   void eh_stretch(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double re, double k);
   void eh_stretch_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom, double re, double k);
   void eh_bend(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double theta0, double kijk);
   void eh_bend_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom, double theta0, double kijk);
   void eh_torsion(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double theta0, double k);
   void eh_torsion_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom, double theta0, double k);
   void ef_d3bj(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double c6, double r2r4_term, double a1, double a2, double s6, double s18);
   void ef_d3bj_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom, double c6, double r2r4_term, double a1, double a2, double s6, double s18);
   void ef_coul_gauss(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double q1, double q2, double sigma1, double sigma2);
   void ef_coul_gauss_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom, double q1, double q2, double sigma1, double sigma2);
} // namespace rx

#endif // COORDS_META_H
