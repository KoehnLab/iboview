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

#ifndef XC_META_H
#define XC_META_H

namespace rx {
   typedef std::size_t
      index_t;

   void xc_diracx(double &E, double *dE, double const *D, double const Factor);
   void xc_ldac_c16(double &E, double *dE, double const *D, double const Factor);
   void xc_ldac_ck16(double &E, double *dE, double const *D, double const Factor);
   void xc_pw92c(double &E, double *dE, double const *D, double const Factor);
   void xc_pbex(double &E, double *dE, double const *D, double const Factor);
   void xc_pbec(double &E, double *dE, double const *D, double const Factor);
   void xc_tpssx(double &E, double *dE, double const *D, double const Factor);
   void xc_tpssc(double &E, double *dE, double const *D, double const Factor);
   void xc_ek_pc07_a(double &E, double *dE, double const *D, double const Factor);
   void xc_ek_pc07_b(double &E, double *dE, double const *D, double const Factor);
   void xc_ll_tpss(double &E, double *dE, double const *D, double const Factor);
   void xc_mn15l(double &E, double *dE, double const *D, double const Factor);
} // namespace rx

#endif // XC_META_H
