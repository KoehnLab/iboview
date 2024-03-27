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

/* IrImportTrafo.cpp v20150312 EST [storm, Gerald Knizia] */
#include <stddef.h> // for size_t
namespace ctimp {
static const double sd0 = 5.e-01;
static const double sd1 = 1.7320508075688772;
static const double sd2 = 8.660254037844386e-01;
static const double sd3 = 6.1237243569579458e-01;
static const double sd4 = 2.4494897427831779;
static const double sd5 = 1.5;
static const double sd6 = 7.9056941504209488e-01;
static const double sd7 = 2.3717082451262845;
static const double sd8 = 3.872983346207417;
static const double sd9 = 1.9364916731037085;
static const double sda = 3.75e-01;
static const double sdb = 7.5e-01;
static const double sdc = 3.;
static const double sdd = 1.1180339887498947;
static const double sde = 6.7082039324993676;
static const double sdf = 3.1622776601683791;
static const double sd10 = 7.3950997288745202e-01;
static const double sd11 = 4.4370598373247123;
static const double sd12 = 5.5901699437494734e-01;
static const double sd13 = 3.3541019662496838;
static const double sd14 = 2.9580398915498081;
static const double sd15 = 2.0916500663351889;
static const double sd16 = 6.2749501990055672;
static const double sd17 = 4.8412291827592718e-01;
static const double sd18 = 9.6824583655185437e-01;
static const double sd19 = 5.809475019311126;
static const double sd1a = 2.5617376914898995;
static const double sd1b = 5.1234753829797981;
static const double sd1c = 5.2291251658379723e-01;
static const double sd1d = 1.0458250331675947;
static const double sd1e = 4.1833001326703778;
static const double sd1f = 1.5687375497513918;
static const double sd20 = 1.2549900398011134e+01;
static const double sd21 = 8.8741196746494246;
static const double sd22 = 2.2185299186623562;
static const double sd23 = 1.3311179511974137e+01;
static const double sd24 = 3.5078038001005702;
static const double sd25 = 7.0156076002011396;
static const double sd26 = 7.0156076002011403e-01;
static const double sd27 = 1.8750000000000002;
static const double sd28 = 3.7500000000000004;
static const double sd29 = 5.;
static const double sd2a = 1.0246950765959596e+01;
static const double sd2b = 6.7169328938139616e-01;
static const double sd2c = 1.0075399340720942e+01;
static const double sd2d = 9.0571104663683977e-01;
static const double sd2e = 1.8114220932736795;
static const double sd2f = 1.4491376746189438e+01;
static const double sd30 = 2.3268138086232857;
static const double sd31 = 2.3268138086232856e+01;
static const double sd32 = 1.1634069043116428e+01;
static const double sd33 = 4.9607837082461076e-01;
static const double sd34 = 2.4803918541230536;
static const double sd35 = 4.9607837082461073;
static const double sd36 = 2.9764702249476645e+01;
static const double sd37 = 4.5285552331841988e-01;
static const double sd38 = 7.245688373094719;
static const double sd39 = 4.0301597362883772;
static const double sd3a = 1.3433865787627923e+01;
static const double sd3b = 2.7171331399105201;
static const double sd3c = 5.434266279821041;
static const double sd3d = 8.1513994197315611;
static const double sd3e = 2.1737065119284161e+01;
static const double sd3f = 1.984313483298443;
static const double sd40 = 1.9843134832984429e+01;
static const double sd41 = 3.125e-01;
static const double sd42 = 9.375e-01;
static const double sd43 = 5.625;
static const double sd44 = 1.125e+01;
static const double sd45 = 7.4999999999999991;
static const double sd46 = 2.8641098093474002;
static const double sd47 = 5.7282196186948005;
static const double sd48 = 1.1456439237389599e+01;
static const double sd49 = 4.5825756949558407;
// transformation from Molpro's obscure cartesian format to Molpro's equally obscure spherical format.
// (...I hope...)
void Vec_Ca2Sh(double *pSh, double const *pCa, size_t *ls, size_t N)
{
   for (size_t i = 0; i < N; ++ i) {
      switch (ls[i]) {
         case 0: {
            pSh[0] = (1.0*pCa[0]);
            pCa += 1;
            pSh += 1;
            continue;
         }
         case 1: {
            pSh[0] = (1.0*pCa[0]);
            pSh[1] = (1.0*pCa[1]);
            pSh[2] = (1.0*pCa[2]);
            pCa += 3;
            pSh += 3;
            continue;
         }
         case 2: {
            pSh[0] = -(0.6666666666666666*pCa[0])*sd0 - (0.6666666666666666*pCa[1])*sd0 + (0.6666666666666666*pCa[2]);
            pSh[1] = (0.5773502691896258*pCa[3])*sd1;
            pSh[2] = (0.5773502691896258*pCa[4])*sd1;
            pSh[3] = (0.6666666666666666*pCa[0])*sd2 - (0.6666666666666666*pCa[1])*sd2;
            pSh[4] = (0.5773502691896258*pCa[5])*sd1;
            pCa += 6;
            pSh += 5;
            continue;
         }
         case 3: {
            pSh[0] = -(0.29814239699997197*pCa[5])*sd3 + (0.29814239699997197*pCa[7])*sd4 - (0.4*pCa[0])*sd3;
            pSh[1] = -(0.29814239699997197*pCa[3])*sd3 + (0.29814239699997197*pCa[8])*sd4 - (0.4*pCa[1])*sd3;
            pSh[2] = -(0.29814239699997197*pCa[4])*sd5 - (0.29814239699997197*pCa[6])*sd5 + (0.4*pCa[2]);
            pSh[3] = -(0.29814239699997197*pCa[5])*sd7 + (0.4*pCa[0])*sd6;
            pSh[4] = (0.2581988897471611*pCa[9])*sd8;
            pSh[5] = (0.29814239699997197*pCa[3])*sd7 - (0.4*pCa[1])*sd6;
            pSh[6] = (0.29814239699997197*pCa[4])*sd9 - (0.29814239699997197*pCa[6])*sd9;
            pCa += 10;
            pSh += 7;
            continue;
         }
         case 4: {
            pSh[0] = -(0.1301200097264711*pCa[10])*sdc - (0.1301200097264711*pCa[11])*sdc + (0.1301200097264711*pCa[9])*sdb + (0.22857142857142856*pCa[0])*sda + (0.22857142857142856*pCa[1])*sda + (0.22857142857142856*pCa[2]);
            pSh[1] = (0.1126872339638022*pCa[14])*sde - (0.15118578920369088*pCa[3])*sdd - (0.15118578920369088*pCa[5])*sdd;
            pSh[2] = -(0.1126872339638022*pCa[13])*sd7 - (0.15118578920369088*pCa[4])*sd7 + (0.15118578920369088*pCa[7])*sdf;
            pSh[3] = -(0.1301200097264711*pCa[9])*sd11 + (0.22857142857142856*pCa[0])*sd10 + (0.22857142857142856*pCa[1])*sd10;
            pSh[4] = -(0.1126872339638022*pCa[12])*sd7 - (0.15118578920369088*pCa[6])*sd7 + (0.15118578920369088*pCa[8])*sdf;
            pSh[5] = (0.1301200097264711*pCa[10])*sd13 - (0.1301200097264711*pCa[11])*sd13 - (0.22857142857142856*pCa[0])*sd12 + (0.22857142857142856*pCa[1])*sd12;
            pSh[6] = (0.15118578920369088*pCa[3])*sd14 - (0.15118578920369088*pCa[5])*sd14;
            pSh[7] = -(0.1126872339638022*pCa[13])*sd16 + (0.15118578920369088*pCa[4])*sd15;
            pSh[8] = (0.1126872339638022*pCa[12])*sd16 - (0.15118578920369088*pCa[6])*sd15;
            pCa += 15;
            pSh += 9;
            continue;
         }
         case 5: {
            pSh[0] = -(0.04337333657549037*pCa[12])*sd19 + (0.05819143739626463*pCa[3])*sd18 - (0.05819143739626463*pCa[5])*sd19 + (0.0761904761904762*pCa[10])*sd17 + (0.0761904761904762*pCa[14])*sd8 + (0.12698412698412698*pCa[0])*sd17;
            pSh[1] = -(0.04337333657549037*pCa[8])*sd19 - (0.05819143739626463*pCa[17])*sd19 + (0.05819143739626463*pCa[6])*sd18 + (0.0761904761904762*pCa[19])*sd8 + (0.0761904761904762*pCa[1])*sd17 + (0.12698412698412698*pCa[15])*sd17;
            pSh[2] = -(0.05819143739626463*pCa[18])*sd1b + (0.05819143739626463*pCa[9])*sd1b + (0.0761904761904762*pCa[16])*sd1a - (0.0761904761904762*pCa[2])*sd1a;
            pSh[3] = -(0.04337333657549037*pCa[12])*sd20 + (0.05819143739626463*pCa[3])*sd1d + (0.05819143739626463*pCa[5])*sd1e + (0.0761904761904762*pCa[10])*sd1f - (0.12698412698412698*pCa[0])*sd1c;
            pSh[4] = -(0.05039526306789696*pCa[11])*sd21 + (0.05039526306789696*pCa[4])*sd21;
            pSh[5] = (0.04337333657549037*pCa[8])*sd20 - (0.05819143739626463*pCa[17])*sd1e - (0.05819143739626463*pCa[6])*sd1d - (0.0761904761904762*pCa[1])*sd1f + (0.12698412698412698*pCa[15])*sd1c;
            pSh[6] = -(0.04337333657549037*pCa[7])*sd23 + (0.0761904761904762*pCa[16])*sd22 + (0.0761904761904762*pCa[2])*sd22;
            pSh[7] = -(0.05819143739626463*pCa[6])*sd25 + (0.0761904761904762*pCa[1])*sd24 + (0.12698412698412698*pCa[15])*sd26;
            pSh[8] = (0.04337333657549037*pCa[7])*sd28 - (0.05819143739626463*pCa[18])*sd29 - (0.05819143739626463*pCa[9])*sd29 + (0.0761904761904762*pCa[16])*sd27 + (0.0761904761904762*pCa[2])*sd27 + (0.12698412698412698*pCa[20]);
            pSh[9] = -(0.05819143739626463*pCa[3])*sd25 + (0.0761904761904762*pCa[10])*sd24 + (0.12698412698412698*pCa[0])*sd26;
            pSh[10] = -(0.05039526306789696*pCa[11])*sd1b + (0.05039526306789696*pCa[13])*sd2a - (0.05039526306789696*pCa[4])*sd1b;
            pCa += 21;
            pSh += 11;
            continue;
         }
         case 6: {
            pSh[0] = (0.026526119002773005*pCa[10])*sd2c - (0.026526119002773005*pCa[3])*sd2c + (0.06926406926406926*pCa[0])*sd2b - (0.06926406926406926*pCa[21])*sd2b;
            pSh[1] = -(0.017545378532260507*pCa[17])*sd2f - (0.017545378532260507*pCa[8])*sd2f + (0.022972292920210562*pCa[19])*sd2f + (0.02353959545345999*pCa[6])*sd2e + (0.03828715486701761*pCa[15])*sd2d + (0.03828715486701761*pCa[1])*sd2d;
            pSh[2] = -(0.017545378532260507*pCa[7])*sd31 + (0.022972292920210562*pCa[16])*sd32 + (0.03828715486701761*pCa[2])*sd30;
            pSh[3] = -(0.015100657524077793*pCa[12])*sd36 + (0.026526119002773005*pCa[10])*sd34 + (0.026526119002773005*pCa[23])*sd35 + (0.026526119002773005*pCa[3])*sd34 + (0.026526119002773005*pCa[5])*sd35 - (0.06926406926406926*pCa[0])*sd33 - (0.06926406926406926*pCa[21])*sd33;
            pSh[4] = -(0.017545378532260507*pCa[11])*sd31 + (0.022972292920210562*pCa[4])*sd32 + (0.03828715486701761*pCa[22])*sd30;
            pSh[5] = -(0.026526119002773005*pCa[10])*sd37 + (0.026526119002773005*pCa[14])*sd38 + (0.026526119002773005*pCa[23])*sd38 - (0.026526119002773005*pCa[25])*sd38 + (0.026526119002773005*pCa[3])*sd37 - (0.026526119002773005*pCa[5])*sd38 + (0.06926406926406926*pCa[0])*sd37 - (0.06926406926406926*pCa[21])*sd37;
            pSh[6] = -(0.02353959545345999*pCa[6])*sd3a + (0.03828715486701761*pCa[15])*sd39 + (0.03828715486701761*pCa[1])*sd39;
            pSh[7] = -(0.017545378532260507*pCa[18])*sd3e + (0.017545378532260507*pCa[7])*sd3c + (0.022972292920210562*pCa[16])*sd3d + (0.02353959545345999*pCa[9])*sd38 - (0.03828715486701761*pCa[2])*sd3b;
            pSh[8] = -(0.017545378532260507*pCa[17])*sd40 + (0.017545378532260507*pCa[8])*sd40 + (0.03828715486701761*pCa[15])*sd3f - (0.03828715486701761*pCa[1])*sd3f;
            pSh[9] = (0.015100657524077793*pCa[12])*sd44 - (0.026526119002773005*pCa[10])*sd42 - (0.026526119002773005*pCa[14])*sd45 + (0.026526119002773005*pCa[23])*sd43 - (0.026526119002773005*pCa[25])*sd45 - (0.026526119002773005*pCa[3])*sd42 + (0.026526119002773005*pCa[5])*sd43 - (0.06926406926406926*pCa[0])*sd41 - (0.06926406926406926*pCa[21])*sd41 + (0.06926406926406926*pCa[27]);
            pSh[10] = -(0.017545378532260507*pCa[11])*sd3c + (0.017545378532260507*pCa[13])*sd3e - (0.022972292920210562*pCa[4])*sd3d - (0.02353959545345999*pCa[24])*sd38 + (0.03828715486701761*pCa[22])*sd3b;
            pSh[11] = (0.017545378532260507*pCa[11])*sd47 - (0.017545378532260507*pCa[13])*sd48 + (0.022972292920210562*pCa[4])*sd46 - (0.02353959545345999*pCa[24])*sd48 + (0.03828715486701761*pCa[22])*sd46 + (0.03828715486701761*pCa[26])*sd49;
            pSh[12] = -(0.017545378532260507*pCa[18])*sd48 + (0.017545378532260507*pCa[7])*sd47 + (0.022972292920210562*pCa[16])*sd46 - (0.02353959545345999*pCa[9])*sd48 + (0.03828715486701761*pCa[20])*sd49 + (0.03828715486701761*pCa[2])*sd46;
            pCa += 28;
            pSh += 13;
            continue;
         }
      }
   }
   // assert(0);
}

// transformation from Molden's cartesian function order to Molpro's cartesian function order.
// (note: the g-order is specified as such. It really seems to be equal in both programs)
void Vec_CaMolden2CaMolpro(double *pOrb, size_t *ls, size_t N)
{
   for (size_t i = 0; i < N; ++ i) {
      switch (ls[i]) {
         case 0: {
            pOrb += 1;
            continue;
         }
         case 1: {
            pOrb += 3;
            continue;
         }
         case 2: {
            pOrb += 6;
            continue;
         }
         case 3: {
            double c3 = pOrb[4]; double c4 = pOrb[5]; double c5 = pOrb[3]; double c6 = pOrb[8]; double c7 = pOrb[6]; double c8 = pOrb[7];
            pOrb[3] = c3; pOrb[4] = c4; pOrb[5] = c5; pOrb[6] = c6; pOrb[7] = c7; pOrb[8] = c8;
            pOrb += 10;
            continue;
         }
         case 4: {
            pOrb += 15;
            continue;
         }
         case 5: {
            // WARNING: not actually specified in Molden format (only defined up to g)
            pOrb += 21;
            continue;
         }
         case 6: {
            // WARNING: not actually specified in Molden format (only defined up to g)
            pOrb += 28;
            continue;
         }
      }
   }
   // assert(0);
}

// transformation from Molden's spherical function order to Molpro's spherical function order.
void Vec_ShMolden2ShMolpro(double *pOrb, size_t *ls, size_t N)
{
   for (size_t i = 0; i < N; ++ i) {
      switch (ls[i]) {
         case 0: {
            pOrb += 1;
            continue;
         }
         case 1: {
            pOrb += 3;
            continue;
         }
         case 2: {
            double c1 = pOrb[4]; double c2 = pOrb[1]; double c4 = pOrb[2];
            pOrb[1] = c1; pOrb[2] = c2; pOrb[4] = c4;
            pOrb += 5;
            continue;
         }
         case 3: {
            double c0 = pOrb[1]; double c1 = pOrb[2]; double c2 = pOrb[0]; double c3 = pOrb[5]; double c5 = pOrb[6]; double c6 = pOrb[3];
            pOrb[0] = c0; pOrb[1] = c1; pOrb[2] = c2; pOrb[3] = c3; pOrb[5] = c5; pOrb[6] = c6;
            pOrb += 7;
            continue;
         }
         case 4: {
            double c1 = pOrb[4]; double c2 = pOrb[1]; double c3 = pOrb[7]; double c4 = pOrb[2]; double c5 = pOrb[3]; double c6 = pOrb[8]; double c7 = pOrb[5]; double c8 = pOrb[6];
            pOrb[1] = c1; pOrb[2] = c2; pOrb[3] = c3; pOrb[4] = c4; pOrb[5] = c5; pOrb[6] = c6; pOrb[7] = c7; pOrb[8] = c8;
            pOrb += 9;
            continue;
         }
         case 5: {
            // WARNING: not actually specified in Molden format (only defined up to g)
            double c0 = pOrb[1]; double c1 = pOrb[2]; double c2 = pOrb[3]; double c3 = pOrb[5]; double c4 = pOrb[8]; double c5 = pOrb[6]; double c6 = pOrb[7]; double c7 = pOrb[10]; double c8 = pOrb[0]; double c10 = pOrb[4];
            pOrb[0] = c0; pOrb[1] = c1; pOrb[2] = c2; pOrb[3] = c3; pOrb[4] = c4; pOrb[5] = c5; pOrb[6] = c6; pOrb[7] = c7; pOrb[8] = c8; pOrb[10] = c10;
            pOrb += 11;
            continue;
         }
         case 6: {
            // WARNING: not actually specified in Molden format (only defined up to g)
            double c0 = pOrb[11]; double c1 = pOrb[4]; double c2 = pOrb[9]; double c3 = pOrb[7]; double c4 = pOrb[10]; double c5 = pOrb[3]; double c6 = pOrb[12]; double c7 = pOrb[5]; double c9 = pOrb[0]; double c10 = pOrb[6]; double c11 = pOrb[2]; double c12 = pOrb[1];
            pOrb[0] = c0; pOrb[1] = c1; pOrb[2] = c2; pOrb[3] = c3; pOrb[4] = c4; pOrb[5] = c5; pOrb[6] = c6; pOrb[7] = c7; pOrb[9] = c9; pOrb[10] = c10; pOrb[11] = c11; pOrb[12] = c12;
            pOrb += 13;
            continue;
         }
      }
   }
   // assert(0);
}

// transformation from cartesian functions in Molpro order to carteian functions in same order...
// but multiplied with obscure monomial dependent renormalization functions. It is as fun as it looks.
// (i.e., don't ask).
void Vec_Ca2Ca_XSqrtNorm(double *pOrb, size_t *ls, size_t N)
{
   for (size_t i = 0; i < N; ++ i) {
      switch (ls[i]) {
         case 0: {
            pOrb += 1;
            continue;
         }
         case 1: {
            pOrb += 3;
            continue;
         }
         case 2: {
            pOrb[0] = 1.224744871391589 * pOrb[0];
            pOrb[1] = 1.224744871391589 * pOrb[1];
            pOrb[2] = 1.224744871391589 * pOrb[2];
            pOrb[3] = 1.3160740129524924 * pOrb[3];
            pOrb[4] = 1.3160740129524924 * pOrb[4];
            pOrb[5] = 1.3160740129524924 * pOrb[5];
            pOrb += 6;
            continue;
         }
         case 3: {
            pOrb[0] = 1.5811388300841898 * pOrb[0];
            pOrb[1] = 1.5811388300841898 * pOrb[1];
            pOrb[2] = 1.5811388300841898 * pOrb[2];
            pOrb[3] = 1.8314207507423532 * pOrb[3];
            pOrb[4] = 1.8314207507423532 * pOrb[4];
            pOrb[5] = 1.8314207507423532 * pOrb[5];
            pOrb[6] = 1.8314207507423532 * pOrb[6];
            pOrb[7] = 1.8314207507423532 * pOrb[7];
            pOrb[8] = 1.8314207507423532 * pOrb[8];
            pOrb[9] = 1.9679896712654306 * pOrb[9];
            pOrb += 10;
            continue;
         }
         case 4: {
            pOrb[0] = 2.091650066335189 * pOrb[0];
            pOrb[1] = 2.091650066335189 * pOrb[1];
            pOrb[2] = 2.091650066335189 * pOrb[2];
            pOrb[3] = 2.571843361805201 * pOrb[3];
            pOrb[4] = 2.571843361805201 * pOrb[4];
            pOrb[5] = 2.571843361805201 * pOrb[5];
            pOrb[6] = 2.571843361805201 * pOrb[6];
            pOrb[7] = 2.571843361805201 * pOrb[7];
            pOrb[8] = 2.571843361805201 * pOrb[8];
            pOrb[9] = 2.772221685664712 * pOrb[9];
            pOrb[10] = 2.772221685664712 * pOrb[10];
            pOrb[11] = 2.772221685664712 * pOrb[11];
            pOrb[12] = 2.9789460677644746 * pOrb[12];
            pOrb[13] = 2.9789460677644746 * pOrb[13];
            pOrb[14] = 2.9789460677644746 * pOrb[14];
            pOrb += 15;
            continue;
         }
         case 5: {
            pOrb[0] = 2.806243040080456 * pOrb[0];
            pOrb[1] = 3.6228441865473595 * pOrb[1];
            pOrb[2] = 3.6228441865473595 * pOrb[2];
            pOrb[3] = 4.145438318933765 * pOrb[3];
            pOrb[4] = 4.454563371755354 * pOrb[4];
            pOrb[5] = 4.145438318933765 * pOrb[5];
            pOrb[6] = 4.145438318933765 * pOrb[6];
            pOrb[7] = 4.801628809415519 * pOrb[7];
            pOrb[8] = 4.801628809415519 * pOrb[8];
            pOrb[9] = 4.145438318933765 * pOrb[9];
            pOrb[10] = 3.6228441865473595 * pOrb[10];
            pOrb[11] = 4.454563371755354 * pOrb[11];
            pOrb[12] = 4.801628809415519 * pOrb[12];
            pOrb[13] = 4.454563371755354 * pOrb[13];
            pOrb[14] = 3.6228441865473595 * pOrb[14];
            pOrb[15] = 2.806243040080456 * pOrb[15];
            pOrb[16] = 3.6228441865473595 * pOrb[16];
            pOrb[17] = 4.145438318933765 * pOrb[17];
            pOrb[18] = 4.145438318933765 * pOrb[18];
            pOrb[19] = 3.6228441865473595 * pOrb[19];
            pOrb[20] = 2.806243040080456 * pOrb[20];
            pOrb += 21;
            continue;
         }
         case 6: {
            pOrb[0] = 3.799671038392666 * pOrb[0];
            pOrb[1] = 5.110618379809705 * pOrb[1];
            pOrb[2] = 5.110618379809705 * pOrb[2];
            pOrb[3] = 6.139926088146812 * pOrb[3];
            pOrb[4] = 6.597779957941507 * pOrb[4];
            pOrb[5] = 6.139926088146812 * pOrb[5];
            pOrb[6] = 6.51779208550841 * pOrb[6];
            pOrb[7] = 7.549507637978121 * pOrb[7];
            pOrb[8] = 7.549507637978121 * pOrb[8];
            pOrb[9] = 6.51779208550841 * pOrb[9];
            pOrb[10] = 6.139926088146812 * pOrb[10];
            pOrb[11] = 7.549507637978121 * pOrb[11];
            pOrb[12] = 8.137707412866751 * pOrb[12];
            pOrb[13] = 7.549507637978121 * pOrb[13];
            pOrb[14] = 6.139926088146812 * pOrb[14];
            pOrb[15] = 5.110618379809705 * pOrb[15];
            pOrb[16] = 6.597779957941507 * pOrb[16];
            pOrb[17] = 7.549507637978121 * pOrb[17];
            pOrb[18] = 7.549507637978121 * pOrb[18];
            pOrb[19] = 6.597779957941507 * pOrb[19];
            pOrb[20] = 5.110618379809705 * pOrb[20];
            pOrb[21] = 3.799671038392666 * pOrb[21];
            pOrb[22] = 5.110618379809705 * pOrb[22];
            pOrb[23] = 6.139926088146812 * pOrb[23];
            pOrb[24] = 6.51779208550841 * pOrb[24];
            pOrb[25] = 6.139926088146812 * pOrb[25];
            pOrb[26] = 5.110618379809705 * pOrb[26];
            pOrb[27] = 3.799671038392666 * pOrb[27];
            pOrb += 28;
            continue;
         }
      }
   }
   // assert(0);
}

} // namespace ctimp
