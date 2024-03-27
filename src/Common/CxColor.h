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

#ifndef CX_COLOR_H
#define CX_COLOR_H

#include <stdint.h>
#include "CxVec3.h"

namespace ct {

//
struct FColor : public TVector3<float>
{
   typedef TVector3<float> base_type;

   FColor() {};
   FColor(float r, float g, float b) : base_type(r,g,b) {}
   FColor(base_type const &a) : base_type(a) {};

   explicit FColor(uint32_t dwColor);

   // +: mix with white, -: mix with black
   void ModBrightness(float fAmount);

//    operator uint32_t () const;
   uint32_t uint32() const;

   void ToHsv(float &h, float &s, float &v);
};

FColor Hsv(float h, float s, float v);

inline FColor ModBrightness(FColor c, float f) {
   FColor r = c;
   r.ModBrightness(f);
   return r;
}

// invert the components. Should make color compatible with, e.g., inkscape format.
uint32_t irgb(uint32_t rgb);




} // namespace ct


#endif // CX_COLOR_H
