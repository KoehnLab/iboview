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
#include <cmath>
#include <math.h>
#include "xc_meta.h"

namespace rx {
   inline double sqr(double x) { return x*x; }
   #define FORT_Extern(lower,UPPER) lower##_
}


namespace rx {
   // Fab(z) Function for LL_TPSS
   static double PC07_Fab(double z) {
      double const a = 0.5389;
      double const b = 3.;
      if (z <= 0) {
         return 0.;
      } else if (z >= a) {
         return 1.;
      } else {
         double az0 = a / (a - z);
         double az1 = a / z;
         if (az0 > 45) return 1.; // required to avoid NaNs and floating point overflows.
         if (az1 > 45) return 0.; // required to avoid NaNs and floating point overflows.
         double e0 = std::exp(az0);
         double e1 = std::exp(az1);
//          double e0 = std::exp(a / (a - z));   // very large for z ~ a (where F becomes 1)
//          double e1 = std::exp(a / z);         // very large for z ~ 0 (where F becomes 0)
         return std::pow((1. + e0) / (e1 + e0), b);
      }
   }

   // ...and PC07_Fab's first derivative with respect to z.
   static double dPC07_Fab(double z) {
      double const a = 0.5389;
      double const b = 3.;
      if (z <= 0) {
         return 0.;
      } else if (z >= a) {
         return 0.;
      } else {
         double az0 = a / (a - z);
         double az1 = a / z;
         if (az0 > 45) return 0.; // required to avoid NaNs and floating point overflows.
         if (az1 > 45) return 0.; // required to avoid NaNs and floating point overflows.
         double e0 = std::exp(az0);
         double e1 = std::exp(az1);
//          double e0 = std::exp(a / (a - z));   // very large for z ~ a (where F becomes 1)
//          double e1 = std::exp(a / z);         // very large for z ~ 0 (where F becomes 0)
//          return std::pow((1. + e0) / (e1 + e0), b - 1) * b * (((a * e0) / ((e1 + e0) * (std::pow(a - z, 2)))) - (((1 + e0) * ((a * e0 / (sqr(a - z))) - (a * e1 / (sqr(z))))) / (sqr((e1) + (e0)))));
//          assert(b == 3.);
         return sqr((1. + e0) / (e1 + e0)) * b * (((a * e0) / ((e1 + e0) * (sqr(a - z)))) - (((1 + e0) * ((a * e0 / (sqr(a - z))) - (a * e1 / (sqr(z))))) / (sqr((e1) + (e0)))));
      }
   }
}

// kate: syntax c++;


namespace rx {
   #include "xc_diracx.inl"
   #include "xc_ldac_c16.inl"
   #include "xc_ldac_ck16.inl"
   #include "xc_pw92c.inl"
   #include "xc_pbex.inl"
   #include "xc_pbec.inl"
   #include "xc_tpssx.inl"
   #include "xc_tpssc.inl"
   #include "xc_ek_pc07_a.inl"
   #include "xc_ek_pc07_b.inl"
   #include "xc_ll_tpss.inl"
   #include "xc_mn15l.inl"
} // namespace rx
