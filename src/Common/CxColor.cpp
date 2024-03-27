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

#include "CxColor.h"
#include <assert.h>
#include <algorithm> // for max/min

#include <iostream>
#include <boost/format.hpp>

namespace ct {

// def hex2col(h_): # hex color
//    h = h_
//    if type(h_) is int:
//       h = '%06x' % h_
//    if h.startswith('0x'):
//       h = h[2:]
//    assert(type(h) is str and len(h)==6)
//    c = lambda s: int(s,16)/255.
//    return (c(h[0:2]), c(h[2:4]), c(h[4:6]))
//
// def col2hex(col):
//    def c2i(c):
//       return max(min(int(.5 + 255.*c),255),0)
//    return "%02x%02x%02x" % (c2i(col[0]),c2i(col[1]),c2i(col[2]))

inline unsigned f2c(float f) {
   if (f < 0.) return   0;
   if (f > 1.) return 255;
   return (unsigned)(f * 255. + .5);
}

// FColor::operator uint32_t () const {
// //    return (f2c(m[0]) << 16) | (f2c(m[1]) << 8) | f2c(m[2]);
//    return (f2c(m[0]) << 0) | (f2c(m[1]) << 8) | (f2c(m[2]) << 16);
// }

uint32_t FColor::uint32 () const {
//    return (f2c(m[0]) << 16) | (f2c(m[1]) << 8) | f2c(m[2]);
   return (f2c(m[0]) << 0) | (f2c(m[1]) << 8) | (f2c(m[2]) << 16);
}


inline float c2f(unsigned t) {
   return ((float)t)/255.;
}

FColor::FColor(uint32_t dwColor)
//    : base_type(c2f((dwColor >> 16)&0xff), c2f((dwColor >> 8)&0xff), c2f((dwColor >> 0)&0xff))
   : base_type(c2f((dwColor >> 0)&0xff), c2f((dwColor >> 8)&0xff), c2f((dwColor >> 16)&0xff))
{}


void FColor::ModBrightness(float f)
{
   if (f == 0)
      return;
   FColor d;
   if (f > 0)
      d = FColor(1.f, 1.f, 1.f);
   else {
      d = FColor(0.f, 0.f, 0.f);
      f *= -1;
   }
   base_type::operator = ((1.f-f)*(*this) + f * (base_type)d);
}

FColor Hsv(float h, float s, float v)
{
   while (h >= 360.f) h -= 360.f;
   while (h < 0.f) h += 360.f;
   assert(0. <= s && s <= 1.);
   assert(0. <= v && v <= 1.);

   float H = h / 60.f;
   int i = (int)H;
   assert(0 <= i && i <= 5);
   float f = H - (float)i;

   int a[6] = {0,2,1,1,2,0};
   int b[6] = {2,0,0,2,1,1};
   int c[6] = {1,1,2,0,0,2};

   if (i % 2 == 0) f = 1.f-f;
   float V[3] = {v, v*(1.f-s), v*(1.f-s*f)};

   return FColor(V[a[i]], V[b[i]], V[c[i]]);
};

void FColor::ToHsv(float &h, float &s, float &v)
{
   float
      r = m[0], g = m[1], b = m[2];
   assert(0. <= r && r <= 1.);
   assert(0. <= g && g <= 1.);
   assert(0. <= b && b <= 1.);
   float
      M = std::max(std::max(r,g), b),
      m = std::min(std::min(r,g), b);
   if (M == m) {
      h = 0.f;
   } else {
      float denom = 1.f/(M-m);
      if (M == r) {
         h = 60.f*(g-b)*denom;
         if (g < b)
            h += 360.f;
      } else if (M == g) {
         h = 60.f*(b-r)*denom + 120.f;
      } else {
         h = 60.f*(r-g)*denom + 240.f;
      }
   }
   if ( M != 0.f ) {
      s = 1.f - m/M;
   } else {
      s = 0.f;
   }
   v = M;
}


// invert the components. Should make color compatible with, e.g., inkscape format.
uint32_t irgb(uint32_t rgb)
{
   return (((rgb >> 24)&0xff) << 24) |
          (((rgb >> 16)&0xff) <<  0) |
          (((rgb >>  8)&0xff) <<  8) |
          (((rgb >>  0)&0xff) << 16);
//    return (((rgb >> 24)&0xff) <<  0) |
//           (((rgb >> 16)&0xff) <<  8) |
//           (((rgb >>  8)&0xff) << 16) |
//           (((rgb >>  0)&0xff) << 24);
//    return rgb;
}

} // namespace ct
