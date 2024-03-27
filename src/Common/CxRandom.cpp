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

#define _USE_MATH_DEFINES
#include <cmath>
#include "CxTypes.h"
#include "CxRandom.h"

#ifndef M_PI
   #define M_PI 3.14159265358979323846264338327950288419
#endif


namespace ct {

int64_t g_LastSeedOffset = 0;

// holds the random number generator.
struct FRandomNumberGeneratorImpl
{
    uint32_t x,y,z,c; // Seed variables

    // public domain JKISS PRNG recycled from: David Jones - ``Good practice in
    //    (pseudo) random number generation for bioinformatics applications''.
    // This PRNG is much less bad than it looks.
    uint32_t JKISS()
    {
        unsigned long long t;
        x = 314527869 * x + 1234567; y ^= y << 5; y ^= y >> 7; y ^= y << 22;
        t = 4294584393ULL * z + c;
        c = static_cast<uint32_t>(t >> 32); z = static_cast<uint32_t>(t);
        return x + y + z;
   }

    uint32_t Get();
    inline void SetStartValue(int64_t SeedValue_);

    bool HaveNextNormal;
    double NextNormal;
};

FRandomNumberGenerator::FRandomNumberGenerator(int64_t SeedValue_)
   : p(new FRandomNumberGeneratorImpl)
{
    p->SetStartValue(SeedValue_);
    p->HaveNextNormal = false;
}

FRandomNumberGenerator::~FRandomNumberGenerator()
{
    delete p;
    p = 0;
}


uint FRandomNumberGeneratorImpl::Get()
{
    return JKISS();
}

void FRandomNumberGeneratorImpl::SetStartValue(int64_t SeedValue_)
{
    x = 123456789;
    y = 987654321;
    z = 43219876;
    c = 6543217;
    x += (uint32_t)SeedValue_;
    y += (uint32_t)(SeedValue_>>32);
}

size_t FRandomNumberGenerator::GetUniformU(size_t LowerBound, size_t UpperBound) {
    assert(UpperBound >= LowerBound);
    return static_cast<size_t>(p->Get() % (UpperBound - LowerBound + 1) + LowerBound);
}


double FRandomNumberGenerator::GetUniformF(double LowerBound, double UpperBound) {
    assert(UpperBound >= LowerBound);
    // hm.. not enough bits for a double :/.
    uint32_t Mask = (uint32_t(1)<<31u);
    double f = double(p->Get()&(Mask-1)) / (Mask-1);
    return static_cast<double>((1.-f)*LowerBound  + f*UpperBound);
}


double FRandomNumberGenerator::GetStdNormal()
{
   if (p->HaveNextNormal) {
      p->HaveNextNormal = false;
      return p->NextNormal;
   } else {
      // Box-Muller transform.
      double
         u1 = 0., u2 = 0.;
      while (u1 == 0.)
         u1 = GetUniformF(0., 1.);
      while (u2 == 0.)
         u2 = GetUniformF(0., 1.);

      double
         pref_u1 = std::sqrt(-2.*std::log(u1)),
         two_pi_u2 = (2.*M_PI)*u2,
         cos_u2 = std::cos(two_pi_u2),
         sin_u2 = std::sin(two_pi_u2);
      p->NextNormal = pref_u1 * sin_u2;
      p->HaveNextNormal = true;
      return pref_u1 * cos_u2;
   }
}

#ifdef ENABLE_GLOBAL_RNG
FRandomNumberGenerator
   g_SharedRng;
#endif // ENABLE_GLOBAL_RNG

} // namespace ct
