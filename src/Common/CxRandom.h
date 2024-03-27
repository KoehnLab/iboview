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

#ifndef CX_RANDOM_H
#define CX_RANDOM_H

#include <stdint.h>
namespace ct {
   struct FRandomNumberGeneratorImpl;
   typedef unsigned int uint;
   double Second();

   extern int64_t g_LastSeedOffset;

   // holds the random number generator.
   struct FRandomNumberGenerator {
      // returns a random number in the range [LowerBound,UpperBound]
      // (including the lower and the upper bounds in the range.)
      size_t GetUniformU(size_t LowerBound, size_t UpperBound);
      double GetUniformF(double LowerBound, double UpperBound);
      // sample standard normal distribution with mean 0 and variance 1.
      double GetStdNormal();

      // seed the RNG on construction.
      explicit FRandomNumberGenerator(int64_t SeedValue_ = (g_LastSeedOffset++) + static_cast<int64_t>((1000000.0)*Second()));
      ~FRandomNumberGenerator();
   private:
      FRandomNumberGeneratorImpl
         *p;
      FRandomNumberGenerator(FRandomNumberGenerator const &other); // not implemented
      void operator = (FRandomNumberGenerator const &other); // not implemented
   };


#ifdef ENABLE_GLOBAL_RNG
   extern FRandomNumberGenerator
      // WARNING: not thread safe, and not automatically initialized to anything sensible!
      g_SharedRng;
#endif // ENABLE_GLOBAL_RNG
} // namespace ct

#endif // CX_RANDOM_H
