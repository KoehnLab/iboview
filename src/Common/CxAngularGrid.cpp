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

// g++ -Wall -pedantic -c -DNEDEBUG -O3 CxAngularGrid.cpp -o CxAngularGrid.o && nm --demangle --size-sort --print-size --format=b --radix=d CxAngularGrid.o

#include <stdexcept>
#include "CxAngularGrid.h"

#define IR_RP

#ifdef USE_LEBEDEV_GRIDS
   #include <math.h> // need sqrt()
#endif

namespace aig { // angular integration grids

#include "CxAngularGrid_Orbits.inl"

#ifdef USE_LEBEDEV_GRIDS
   #include "CxLebedevGrid.inl"
   #include "CxAngularGrid_Grids_Lebedev.inl"
#else
   #include "CxAngularGrid_Grids.inl"
#endif

size_t GetNumAngularGrids()
{
   return sizeof(s_AngularGridInfos)/sizeof(s_AngularGridInfos[0]);
}


FAngularGridEntry const *GetAngularGridInfo(size_t iGrid)
{
   if (iGrid >= GetNumAngularGrids())
      throw std::runtime_error("attempted to instanciate non-existent angular integration grid.");
   return &s_AngularGridInfos[iGrid];
}


} // namespace aig
