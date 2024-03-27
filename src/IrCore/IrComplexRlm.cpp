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

/* IrComplexRlm.cpp v20200511 EST [mirage, Gerald Knizia] */
#include "CxDefs.h" // for assert
#include <algorithm> // for std::swap
#include "IrComplexRlm.h"

namespace ir_rlm {

template<class T> inline T sqr(T const &x) { return x*x; }

// Evaluate Rlm for l == 0 ... MaxL (inclusive)
// (not very pretty or efficient, especially the Deriv1 function... could be
// made more efficient if computed recursively, like the EvalSlmX functions. But
// don't want to deal with it now. Getting it right is not as easy as it
// looks...)
void EvalRlmX(complex_double *RlmX, double x, double y, double z, int MaxL)
{
   assert(MaxL >= 0);
   RlmX[0] = 1.0;
   if (0 == MaxL) return;
   complex_double xpiy = complex_double(x,y);   // (x + i y)^1
   complex_double xmiy = complex_double(x,-y);  // (x - i y)^1
   RlmX[1] = -1.0/2.0*xmiy;
   RlmX[2] = z;
   RlmX[3] = xpiy;
   if (1 == MaxL) return;
   complex_double xpiy2 = xpiy*xpiy;
   complex_double xmiy2 = xmiy*xmiy;
   double z2 = z*z;
   RlmX[4] = (1.0/8.0)*xmiy2;
   RlmX[5] = -1.0/2.0*xmiy*z;
   RlmX[6] = -1.0/2.0*xmiy*xpiy + z2;
   RlmX[7] = 3.0*xpiy*z;
   RlmX[8] = 3.0*xpiy2;
   if (2 == MaxL) return;
   complex_double xpiy3 = xpiy2*xpiy;
   complex_double xmiy3 = xmiy2*xmiy;
   double z3 = z2*z;
   RlmX[9] = -1.0/48.0*xmiy3;
   RlmX[10] = (1.0/8.0)*xmiy2*z;
   RlmX[11] = (1.0/8.0)*xmiy*(xmiy*xpiy - 4.0*z2);
   RlmX[12] = z*(-3.0/2.0*xmiy*xpiy + z2);
   RlmX[13] = (3.0/2.0)*xpiy*(-xmiy*xpiy + 4.0*z2);
   RlmX[14] = 15.0*xpiy2*z;
   RlmX[15] = 15.0*xpiy3;
   if (3 == MaxL) return;
   complex_double xpiy4 = xpiy3*xpiy;
   complex_double xmiy4 = xmiy3*xmiy;
   double z4 = z3*z;
   RlmX[16] = (1.0/384.0)*xmiy4;
   RlmX[17] = -1.0/48.0*xmiy3*z;
   RlmX[18] = (1.0/48.0)*xmiy2*(-xmiy*xpiy + 6.0*z2);
   RlmX[19] = (1.0/8.0)*xmiy*z*(3.0*xmiy*xpiy - 4.0*z2);
   RlmX[20] = (3.0/8.0)*xmiy2*xpiy2 - 3.0*xmiy*xpiy*z2 + z4;
   RlmX[21] = (5.0/2.0)*xpiy*z*(-3.0*xmiy*xpiy + 4.0*z2);
   RlmX[22] = (15.0/2.0)*xpiy2*(-xmiy*xpiy + 6.0*z2);
   RlmX[23] = 105.0*xpiy3*z;
   RlmX[24] = 105.0*xpiy4;
   if (4 == MaxL) return;
   complex_double xpiy5 = xpiy4*xpiy;
   complex_double xmiy5 = xmiy4*xmiy;
   double z5 = z4*z;
   RlmX[25] = -1.0/3840.0*xmiy5;
   RlmX[26] = (1.0/384.0)*xmiy4*z;
   RlmX[27] = (1.0/384.0)*xmiy3*(xmiy*xpiy - 8.0*z2);
   RlmX[28] = (1.0/16.0)*xmiy2*z*(-xmiy*xpiy + 2.0*z2);
   RlmX[29] = (1.0/16.0)*xmiy*(-xmiy2*xpiy2 + 12.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[30] = z*((15.0/8.0)*xmiy2*xpiy2 - 5.0*xmiy*xpiy*z2 + z4);
   RlmX[31] = (15.0/8.0)*xpiy*(xmiy2*xpiy2 - 12.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[32] = (105.0/2.0)*xpiy2*z*(-xmiy*xpiy + 2.0*z2);
   RlmX[33] = (105.0/2.0)*xpiy3*(-xmiy*xpiy + 8.0*z2);
   RlmX[34] = 945.0*xpiy4*z;
   RlmX[35] = 945.0*xpiy5;
   if (5 == MaxL) return;
   complex_double xpiy6 = xpiy5*xpiy;
   complex_double xmiy6 = xmiy5*xmiy;
   double z6 = z5*z;
   RlmX[36] = (1.0/46080.0)*xmiy6;
   RlmX[37] = -1.0/3840.0*xmiy5*z;
   RlmX[38] = (1.0/3840.0)*xmiy4*(-xmiy*xpiy + 10.0*z2);
   RlmX[39] = (1.0/384.0)*xmiy3*z*(3.0*xmiy*xpiy - 8.0*z2);
   RlmX[40] = (1.0/128.0)*xmiy2*(xmiy2*xpiy2 - 16.0*xmiy*xpiy*z2 + 16.0*z4);
   RlmX[41] = (1.0/16.0)*xmiy*z*(-5.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[42] = -5.0/16.0*xmiy3*xpiy3 + (45.0/8.0)*xmiy2*xpiy2*z2 - 15.0/2.0*xmiy*xpiy*z4 + z6;
   RlmX[43] = (21.0/8.0)*xpiy*z*(5.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[44] = (105.0/8.0)*xpiy2*(xmiy2*xpiy2 - 16.0*xmiy*xpiy*z2 + 16.0*z4);
   RlmX[45] = (315.0/2.0)*xpiy3*z*(-3.0*xmiy*xpiy + 8.0*z2);
   RlmX[46] = (945.0/2.0)*xpiy4*(-xmiy*xpiy + 10.0*z2);
   RlmX[47] = 10395.0*xpiy5*z;
   RlmX[48] = 10395.0*xpiy6;
   if (6 == MaxL) return;
   complex_double xpiy7 = xpiy6*xpiy;
   complex_double xmiy7 = xmiy6*xmiy;
   double z7 = z6*z;
   RlmX[49] = -1.0/645120.0*xmiy7;
   RlmX[50] = (1.0/46080.0)*xmiy6*z;
   RlmX[51] = (1.0/46080.0)*xmiy5*(xmiy*xpiy - 12.0*z2);
   RlmX[52] = (1.0/3840.0)*xmiy4*z*(-3.0*xmiy*xpiy + 10.0*z2);
   RlmX[53] = xmiy3*(-1.0/1280.0*xmiy2*xpiy2 + (1.0/64.0)*xmiy*xpiy*z2 - 1.0/48.0*z4);
   RlmX[54] = (1.0/384.0)*xmiy2*z*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   RlmX[55] = (1.0/128.0)*xmiy*(5.0*xmiy3*xpiy3 - 120.0*xmiy2*xpiy2*z2 + 240.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[56] = (1.0/16.0)*z*(-35.0*xmiy3*xpiy3 + 210.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 16.0*z6);
   RlmX[57] = (7.0/16.0)*xpiy*(-5.0*xmiy3*xpiy3 + 120.0*xmiy2*xpiy2*z2 - 240.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[58] = (63.0/8.0)*xpiy2*z*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   RlmX[59] = (315.0/8.0)*xpiy3*(3.0*xmiy2*xpiy2 - 60.0*xmiy*xpiy*z2 + 80.0*z4);
   RlmX[60] = (3465.0/2.0)*xpiy4*z*(-3.0*xmiy*xpiy + 10.0*z2);
   RlmX[61] = (10395.0/2.0)*xpiy5*(-xmiy*xpiy + 12.0*z2);
   RlmX[62] = 135135.0*xpiy6*z;
   RlmX[63] = 135135.0*xpiy7;
   if (7 == MaxL) return;
   complex_double xpiy8 = xpiy7*xpiy;
   complex_double xmiy8 = xmiy7*xmiy;
   double z8 = z7*z;
   RlmX[64] = (1.0/10321920.0)*xmiy8;
   RlmX[65] = -1.0/645120.0*xmiy7*z;
   RlmX[66] = (1.0/645120.0)*xmiy6*(-xmiy*xpiy + 14.0*z2);
   RlmX[67] = (1.0/15360.0)*xmiy5*z*(xmiy*xpiy - 4.0*z2);
   RlmX[68] = (1.0/15360.0)*xmiy4*(xmiy2*xpiy2 - 24.0*xmiy*xpiy*z2 + 40.0*z4);
   RlmX[69] = (1.0/768.0)*xmiy3*z*(-3.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 16.0*z4);
   RlmX[70] = (1.0/256.0)*xmiy2*(-xmiy3*xpiy3 + 30.0*xmiy2*xpiy2*z2 - 80.0*xmiy*xpiy*z4 + 32.0*z6);
   RlmX[71] = (1.0/128.0)*xmiy*z*(35.0*xmiy3*xpiy3 - 280.0*xmiy2*xpiy2*z2 + 336.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[72] = (35.0/128.0)*xmiy4*xpiy4 - 35.0/4.0*xmiy3*xpiy3*z2 + (105.0/4.0)*xmiy2*xpiy2*z4 - 14.0*xmiy*xpiy*z6 + z8;
   RlmX[73] = (9.0/16.0)*xpiy*z*(-35.0*xmiy3*xpiy3 + 280.0*xmiy2*xpiy2*z2 - 336.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[74] = (315.0/16.0)*xpiy2*(-xmiy3*xpiy3 + 30.0*xmiy2*xpiy2*z2 - 80.0*xmiy*xpiy*z4 + 32.0*z6);
   RlmX[75] = (3465.0/8.0)*xpiy3*z*(3.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 16.0*z4);
   RlmX[76] = (10395.0/8.0)*xpiy4*(xmiy2*xpiy2 - 24.0*xmiy*xpiy*z2 + 40.0*z4);
   RlmX[77] = (135135.0/2.0)*xpiy5*z*(-xmiy*xpiy + 4.0*z2);
   RlmX[78] = (135135.0/2.0)*xpiy6*(-xmiy*xpiy + 14.0*z2);
   RlmX[79] = 2027025.0*xpiy7*z;
   RlmX[80] = 2027025.0*xpiy8;
   if (8 == MaxL) return;
   complex_double xpiy9 = xpiy8*xpiy;
   complex_double xmiy9 = xmiy8*xmiy;
   double z9 = z8*z;
   RlmX[81] = -1.0/185794560.0*xmiy9;
   RlmX[82] = (1.0/10321920.0)*xmiy8*z;
   RlmX[83] = (1.0/10321920.0)*xmiy7*(xmiy*xpiy - 16.0*z2);
   RlmX[84] = (1.0/645120.0)*xmiy6*z*(-3.0*xmiy*xpiy + 14.0*z2);
   RlmX[85] = (1.0/215040.0)*xmiy5*(-xmiy2*xpiy2 + 28.0*xmiy*xpiy*z2 - 56.0*z4);
   RlmX[86] = (1.0/3072.0)*xmiy4*z*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[87] = (1.0/3072.0)*xmiy3*(xmiy3*xpiy3 - 36.0*xmiy2*xpiy2*z2 + 120.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[88] = (1.0/256.0)*xmiy2*z*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   RlmX[89] = (1.0/256.0)*xmiy*(-7.0*xmiy4*xpiy4 + 280.0*xmiy3*xpiy3*z2 - 1120.0*xmiy2*xpiy2*z4 + 896.0*xmiy*xpiy*z6 - 128.0*z8);
   RlmX[90] = (1.0/128.0)*z*(315.0*xmiy4*xpiy4 - 3360.0*xmiy3*xpiy3*z2 + 6048.0*xmiy2*xpiy2*z4 - 2304.0*xmiy*xpiy*z6 + 128.0*z8);
   RlmX[91] = (45.0/128.0)*xpiy*(7.0*xmiy4*xpiy4 - 280.0*xmiy3*xpiy3*z2 + 1120.0*xmiy2*xpiy2*z4 - 896.0*xmiy*xpiy*z6 + 128.0*z8);
   RlmX[92] = (495.0/16.0)*xpiy2*z*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   RlmX[93] = (3465.0/16.0)*xpiy3*(-xmiy3*xpiy3 + 36.0*xmiy2*xpiy2*z2 - 120.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[94] = (135135.0/8.0)*xpiy4*z*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[95] = (135135.0/8.0)*xpiy5*(xmiy2*xpiy2 - 28.0*xmiy*xpiy*z2 + 56.0*z4);
   RlmX[96] = (675675.0/2.0)*xpiy6*z*(-3.0*xmiy*xpiy + 14.0*z2);
   RlmX[97] = (2027025.0/2.0)*xpiy7*(-xmiy*xpiy + 16.0*z2);
   RlmX[98] = 34459425.0*xpiy8*z;
   RlmX[99] = 34459425.0*xpiy9;
   if (9 == MaxL) return;
   complex_double xpiy10 = xpiy9*xpiy;
   complex_double xmiy10 = xmiy9*xmiy;
   double z10 = z9*z;
   RlmX[100] = (1.0/3715891200.0)*xmiy10;
   RlmX[101] = -1.0/185794560.0*xmiy9*z;
   RlmX[102] = (1.0/185794560.0)*xmiy8*(-xmiy*xpiy + 18.0*z2);
   RlmX[103] = (1.0/10321920.0)*xmiy7*z*(3.0*xmiy*xpiy - 16.0*z2);
   RlmX[104] = xmiy6*((1.0/3440640.0)*xmiy2*xpiy2 - 1.0/107520.0*xmiy*xpiy*z2 + (1.0/46080.0)*z4);
   RlmX[105] = (1.0/645120.0)*xmiy5*z*(-15.0*xmiy2*xpiy2 + 140.0*xmiy*xpiy*z2 - 168.0*z4);
   RlmX[106] = (1.0/43008.0)*xmiy4*(-xmiy3*xpiy3 + 42.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 112.0*z6);
   RlmX[107] = (1.0/3072.0)*xmiy3*z*(7.0*xmiy3*xpiy3 - 84.0*xmiy2*xpiy2*z2 + 168.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[108] = (1.0/3072.0)*xmiy2*(7.0*xmiy4*xpiy4 - 336.0*xmiy3*xpiy3*z2 + 1680.0*xmiy2*xpiy2*z4 - 1792.0*xmiy*xpiy*z6 + 384.0*z8);
   RlmX[109] = (1.0/256.0)*xmiy*z*(-63.0*xmiy4*xpiy4 + 840.0*xmiy3*xpiy3*z2 - 2016.0*xmiy2*xpiy2*z4 + 1152.0*xmiy*xpiy*z6 - 128.0*z8);
   RlmX[110] = -63.0/256.0*xmiy5*xpiy5 + (1575.0/128.0)*xmiy4*xpiy4*z2 - 525.0/8.0*xmiy3*xpiy3*z4 + (315.0/4.0)*xmiy2*xpiy2*z6 - 45.0/2.0*xmiy*xpiy*z8 + z10;
   RlmX[111] = (55.0/128.0)*xpiy*z*(63.0*xmiy4*xpiy4 - 840.0*xmiy3*xpiy3*z2 + 2016.0*xmiy2*xpiy2*z4 - 1152.0*xmiy*xpiy*z6 + 128.0*z8);
   RlmX[112] = (495.0/128.0)*xpiy2*(7.0*xmiy4*xpiy4 - 336.0*xmiy3*xpiy3*z2 + 1680.0*xmiy2*xpiy2*z4 - 1792.0*xmiy*xpiy*z6 + 384.0*z8);
   RlmX[113] = (6435.0/16.0)*xpiy3*z*(-7.0*xmiy3*xpiy3 + 84.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[114] = (45045.0/16.0)*xpiy4*(-xmiy3*xpiy3 + 42.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 112.0*z6);
   RlmX[115] = (135135.0/8.0)*xpiy5*z*(15.0*xmiy2*xpiy2 - 140.0*xmiy*xpiy*z2 + 168.0*z4);
   RlmX[116] = xpiy6*((2027025.0/8.0)*xmiy2*xpiy2 - 8108100.0*xmiy*xpiy*z2 + 18918900.0*z4);
   RlmX[117] = (11486475.0/2.0)*xpiy7*z*(-3.0*xmiy*xpiy + 16.0*z2);
   RlmX[118] = (34459425.0/2.0)*xpiy8*(-xmiy*xpiy + 18.0*z2);
   RlmX[119] = 654729075.0*xpiy9*z;
   RlmX[120] = 654729075.0*xpiy10;
   if (10 == MaxL) return;
   // suppress potential unused warning(s)
   (void)z10;
}

// Evaluate Rlm for l == 0 ... MaxL (inclusive), and derivative components [1,d/dx,d/dy,d/dz]
void EvalRlmX_Deriv1(complex_double *RlmX, double x, double y, double z, int MaxL)
{
   complex_double const I(0,1);
   assert(MaxL >= 0);
   RlmX[4*0+0] = 1.0;
   RlmX[4*0+1] = 0;
   RlmX[4*0+2] = 0;
   RlmX[4*0+3] = 0;
   if (0 == MaxL) return;
   complex_double xpiy = complex_double(x,y);   // (x + i y)^1
   complex_double xmiy = complex_double(x,-y);  // (x - i y)^1
   RlmX[4*1+0] = -1.0/2.0*xmiy;
   complex_double dR1m1_xmiy = -1.0/2.0;
   RlmX[4*1+1] = dR1m1_xmiy;
   RlmX[4*1+2] = -I*dR1m1_xmiy;
   RlmX[4*1+3] = 0;
   RlmX[4*2+0] = z;
   RlmX[4*2+1] = 0;
   RlmX[4*2+2] = 0;
   RlmX[4*2+3] = 1.0;
   RlmX[4*3+0] = xpiy;
   RlmX[4*3+1] = 1.0;
   RlmX[4*3+2] = I;
   RlmX[4*3+3] = 0;
   if (1 == MaxL) return;
   complex_double xpiy2 = xpiy*xpiy;
   complex_double xmiy2 = xmiy*xmiy;
   double z2 = z*z;
   RlmX[4*4+0] = (1.0/8.0)*xmiy2;
   complex_double dR2m2_xmiy = (1.0/4.0)*xmiy;
   RlmX[4*4+1] = dR2m2_xmiy;
   RlmX[4*4+2] = -I*dR2m2_xmiy;
   RlmX[4*4+3] = 0;
   RlmX[4*5+0] = -1.0/2.0*xmiy*z;
   complex_double dR2m1_xmiy = -1.0/2.0*z;
   RlmX[4*5+1] = dR2m1_xmiy;
   RlmX[4*5+2] = -I*dR2m1_xmiy;
   RlmX[4*5+3] = -1.0/2.0*xmiy;
   RlmX[4*6+0] = -1.0/2.0*xmiy*xpiy + z2;
   complex_double dR2n0_xpiy = -1.0/2.0*xmiy;
   complex_double dR2n0_xmiy = -1.0/2.0*xpiy;
   RlmX[4*6+1] = dR2n0_xmiy + dR2n0_xpiy;
   RlmX[4*6+2] = -I*dR2n0_xmiy + I*dR2n0_xpiy;
   RlmX[4*6+3] = 2.0*z;
   RlmX[4*7+0] = 3.0*xpiy*z;
   complex_double dR2p1_xpiy = 3.0*z;
   RlmX[4*7+1] = dR2p1_xpiy;
   RlmX[4*7+2] = I*dR2p1_xpiy;
   RlmX[4*7+3] = 3.0*xpiy;
   RlmX[4*8+0] = 3.0*xpiy2;
   complex_double dR2p2_xpiy = 6.0*xpiy;
   RlmX[4*8+1] = dR2p2_xpiy;
   RlmX[4*8+2] = I*dR2p2_xpiy;
   RlmX[4*8+3] = 0;
   if (2 == MaxL) return;
   complex_double xpiy3 = xpiy2*xpiy;
   complex_double xmiy3 = xmiy2*xmiy;
   double z3 = z2*z;
   RlmX[4*9+0] = -1.0/48.0*xmiy3;
   complex_double dR3m3_xmiy = -1.0/16.0*xmiy2;
   RlmX[4*9+1] = dR3m3_xmiy;
   RlmX[4*9+2] = -I*dR3m3_xmiy;
   RlmX[4*9+3] = 0;
   RlmX[4*10+0] = (1.0/8.0)*xmiy2*z;
   complex_double dR3m2_xmiy = (1.0/4.0)*xmiy*z;
   RlmX[4*10+1] = dR3m2_xmiy;
   RlmX[4*10+2] = -I*dR3m2_xmiy;
   RlmX[4*10+3] = (1.0/8.0)*xmiy2;
   RlmX[4*11+0] = (1.0/8.0)*xmiy*(xmiy*xpiy - 4.0*z2);
   complex_double dR3m1_xpiy = (1.0/8.0)*xmiy2;
   complex_double dR3m1_xmiy = (1.0/4.0)*xmiy*xpiy - 1.0/2.0*z2;
   RlmX[4*11+1] = dR3m1_xmiy + dR3m1_xpiy;
   RlmX[4*11+2] = -I*dR3m1_xmiy + I*dR3m1_xpiy;
   RlmX[4*11+3] = -xmiy*z;
   RlmX[4*12+0] = z*(-3.0/2.0*xmiy*xpiy + z2);
   complex_double dR3n0_xpiy = -3.0/2.0*xmiy*z;
   complex_double dR3n0_xmiy = -3.0/2.0*xpiy*z;
   RlmX[4*12+1] = dR3n0_xmiy + dR3n0_xpiy;
   RlmX[4*12+2] = -I*dR3n0_xmiy + I*dR3n0_xpiy;
   RlmX[4*12+3] = -3.0/2.0*xmiy*xpiy + 3.0*z2;
   RlmX[4*13+0] = (3.0/2.0)*xpiy*(-xmiy*xpiy + 4.0*z2);
   complex_double dR3p1_xpiy = -3.0*xmiy*xpiy + 6.0*z2;
   complex_double dR3p1_xmiy = -3.0/2.0*xpiy2;
   RlmX[4*13+1] = dR3p1_xmiy + dR3p1_xpiy;
   RlmX[4*13+2] = -I*dR3p1_xmiy + I*dR3p1_xpiy;
   RlmX[4*13+3] = 12.0*xpiy*z;
   RlmX[4*14+0] = 15.0*xpiy2*z;
   complex_double dR3p2_xpiy = 30.0*xpiy*z;
   RlmX[4*14+1] = dR3p2_xpiy;
   RlmX[4*14+2] = I*dR3p2_xpiy;
   RlmX[4*14+3] = 15.0*xpiy2;
   RlmX[4*15+0] = 15.0*xpiy3;
   complex_double dR3p3_xpiy = 45.0*xpiy2;
   RlmX[4*15+1] = dR3p3_xpiy;
   RlmX[4*15+2] = I*dR3p3_xpiy;
   RlmX[4*15+3] = 0;
   if (3 == MaxL) return;
   complex_double xpiy4 = xpiy3*xpiy;
   complex_double xmiy4 = xmiy3*xmiy;
   double z4 = z3*z;
   RlmX[4*16+0] = (1.0/384.0)*xmiy4;
   complex_double dR4m4_xmiy = (1.0/96.0)*xmiy3;
   RlmX[4*16+1] = dR4m4_xmiy;
   RlmX[4*16+2] = -I*dR4m4_xmiy;
   RlmX[4*16+3] = 0;
   RlmX[4*17+0] = -1.0/48.0*xmiy3*z;
   complex_double dR4m3_xmiy = -1.0/16.0*xmiy2*z;
   RlmX[4*17+1] = dR4m3_xmiy;
   RlmX[4*17+2] = -I*dR4m3_xmiy;
   RlmX[4*17+3] = -1.0/48.0*xmiy3;
   RlmX[4*18+0] = (1.0/48.0)*xmiy2*(-xmiy*xpiy + 6.0*z2);
   complex_double dR4m2_xpiy = -1.0/48.0*xmiy3;
   complex_double dR4m2_xmiy = (1.0/16.0)*xmiy*(-xmiy*xpiy + 4.0*z2);
   RlmX[4*18+1] = dR4m2_xmiy + dR4m2_xpiy;
   RlmX[4*18+2] = -I*dR4m2_xmiy + I*dR4m2_xpiy;
   RlmX[4*18+3] = (1.0/4.0)*xmiy2*z;
   RlmX[4*19+0] = (1.0/8.0)*xmiy*z*(3.0*xmiy*xpiy - 4.0*z2);
   complex_double dR4m1_xpiy = (3.0/8.0)*xmiy2*z;
   complex_double dR4m1_xmiy = (1.0/4.0)*z*(3.0*xmiy*xpiy - 2.0*z2);
   RlmX[4*19+1] = dR4m1_xmiy + dR4m1_xpiy;
   RlmX[4*19+2] = -I*dR4m1_xmiy + I*dR4m1_xpiy;
   RlmX[4*19+3] = -xmiy*z2 + (1.0/8.0)*xmiy*(3.0*xmiy*xpiy - 4.0*z2);
   RlmX[4*20+0] = (3.0/8.0)*xmiy2*xpiy2 - 3.0*xmiy*xpiy*z2 + z4;
   complex_double dR4n0_xpiy = (3.0/4.0)*xmiy*(xmiy*xpiy - 4.0*z2);
   complex_double dR4n0_xmiy = (3.0/4.0)*xpiy*(xmiy*xpiy - 4.0*z2);
   RlmX[4*20+1] = dR4n0_xmiy + dR4n0_xpiy;
   RlmX[4*20+2] = -I*dR4n0_xmiy + I*dR4n0_xpiy;
   RlmX[4*20+3] = -6.0*xmiy*xpiy*z + 4.0*z3;
   RlmX[4*21+0] = (5.0/2.0)*xpiy*z*(-3.0*xmiy*xpiy + 4.0*z2);
   complex_double dR4p1_xpiy = -15.0*xmiy*xpiy*z + 10.0*z3;
   complex_double dR4p1_xmiy = -15.0/2.0*xpiy2*z;
   RlmX[4*21+1] = dR4p1_xmiy + dR4p1_xpiy;
   RlmX[4*21+2] = -I*dR4p1_xmiy + I*dR4p1_xpiy;
   RlmX[4*21+3] = 20.0*xpiy*z2 + (5.0/2.0)*xpiy*(-3.0*xmiy*xpiy + 4.0*z2);
   RlmX[4*22+0] = (15.0/2.0)*xpiy2*(-xmiy*xpiy + 6.0*z2);
   complex_double dR4p2_xpiy = (45.0/2.0)*xpiy*(-xmiy*xpiy + 4.0*z2);
   complex_double dR4p2_xmiy = -15.0/2.0*xpiy3;
   RlmX[4*22+1] = dR4p2_xmiy + dR4p2_xpiy;
   RlmX[4*22+2] = -I*dR4p2_xmiy + I*dR4p2_xpiy;
   RlmX[4*22+3] = 90.0*xpiy2*z;
   RlmX[4*23+0] = 105.0*xpiy3*z;
   complex_double dR4p3_xpiy = 315.0*xpiy2*z;
   RlmX[4*23+1] = dR4p3_xpiy;
   RlmX[4*23+2] = I*dR4p3_xpiy;
   RlmX[4*23+3] = 105.0*xpiy3;
   RlmX[4*24+0] = 105.0*xpiy4;
   complex_double dR4p4_xpiy = 420.0*xpiy3;
   RlmX[4*24+1] = dR4p4_xpiy;
   RlmX[4*24+2] = I*dR4p4_xpiy;
   RlmX[4*24+3] = 0;
   if (4 == MaxL) return;
   complex_double xpiy5 = xpiy4*xpiy;
   complex_double xmiy5 = xmiy4*xmiy;
   double z5 = z4*z;
   RlmX[4*25+0] = -1.0/3840.0*xmiy5;
   complex_double dR5m5_xmiy = -1.0/768.0*xmiy4;
   RlmX[4*25+1] = dR5m5_xmiy;
   RlmX[4*25+2] = -I*dR5m5_xmiy;
   RlmX[4*25+3] = 0;
   RlmX[4*26+0] = (1.0/384.0)*xmiy4*z;
   complex_double dR5m4_xmiy = (1.0/96.0)*xmiy3*z;
   RlmX[4*26+1] = dR5m4_xmiy;
   RlmX[4*26+2] = -I*dR5m4_xmiy;
   RlmX[4*26+3] = (1.0/384.0)*xmiy4;
   RlmX[4*27+0] = (1.0/384.0)*xmiy3*(xmiy*xpiy - 8.0*z2);
   complex_double dR5m3_xpiy = (1.0/384.0)*xmiy4;
   complex_double dR5m3_xmiy = (1.0/96.0)*xmiy2*(xmiy*xpiy - 6.0*z2);
   RlmX[4*27+1] = dR5m3_xmiy + dR5m3_xpiy;
   RlmX[4*27+2] = -I*dR5m3_xmiy + I*dR5m3_xpiy;
   RlmX[4*27+3] = -1.0/24.0*xmiy3*z;
   RlmX[4*28+0] = (1.0/16.0)*xmiy2*z*(-xmiy*xpiy + 2.0*z2);
   complex_double dR5m2_xpiy = -1.0/16.0*xmiy3*z;
   complex_double dR5m2_xmiy = (1.0/16.0)*xmiy*z*(-3.0*xmiy*xpiy + 4.0*z2);
   RlmX[4*28+1] = dR5m2_xmiy + dR5m2_xpiy;
   RlmX[4*28+2] = -I*dR5m2_xmiy + I*dR5m2_xpiy;
   RlmX[4*28+3] = (1.0/4.0)*xmiy2*z2 + (1.0/16.0)*xmiy2*(-xmiy*xpiy + 2.0*z2);
   RlmX[4*29+0] = (1.0/16.0)*xmiy*(-xmiy2*xpiy2 + 12.0*xmiy*xpiy*z2 - 8.0*z4);
   complex_double dR5m1_xpiy = (1.0/8.0)*xmiy2*(-xmiy*xpiy + 6.0*z2);
   complex_double dR5m1_xmiy = -3.0/16.0*xmiy2*xpiy2 + (3.0/2.0)*xmiy*xpiy*z2 - 1.0/2.0*z4;
   RlmX[4*29+1] = dR5m1_xmiy + dR5m1_xpiy;
   RlmX[4*29+2] = -I*dR5m1_xmiy + I*dR5m1_xpiy;
   RlmX[4*29+3] = (1.0/16.0)*xmiy*(24.0*xmiy*xpiy*z - 32.0*z3);
   RlmX[4*30+0] = z*((15.0/8.0)*xmiy2*xpiy2 - 5.0*xmiy*xpiy*z2 + z4);
   complex_double dR5n0_xpiy = (5.0/4.0)*xmiy*z*(3.0*xmiy*xpiy - 4.0*z2);
   complex_double dR5n0_xmiy = (5.0/4.0)*xpiy*z*(3.0*xmiy*xpiy - 4.0*z2);
   RlmX[4*30+1] = dR5n0_xmiy + dR5n0_xpiy;
   RlmX[4*30+2] = -I*dR5n0_xmiy + I*dR5n0_xpiy;
   RlmX[4*30+3] = (15.0/8.0)*xmiy2*xpiy2 - 5.0*xmiy*xpiy*z2 + z4 + z*(-10.0*xmiy*xpiy*z + 4.0*z3);
   RlmX[4*31+0] = (15.0/8.0)*xpiy*(xmiy2*xpiy2 - 12.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dR5p1_xpiy = (45.0/8.0)*xmiy2*xpiy2 - 45.0*xmiy*xpiy*z2 + 15.0*z4;
   complex_double dR5p1_xmiy = (15.0/4.0)*xpiy2*(xmiy*xpiy - 6.0*z2);
   RlmX[4*31+1] = dR5p1_xmiy + dR5p1_xpiy;
   RlmX[4*31+2] = -I*dR5p1_xmiy + I*dR5p1_xpiy;
   RlmX[4*31+3] = (15.0/8.0)*xpiy*(-24.0*xmiy*xpiy*z + 32.0*z3);
   RlmX[4*32+0] = (105.0/2.0)*xpiy2*z*(-xmiy*xpiy + 2.0*z2);
   complex_double dR5p2_xpiy = (105.0/2.0)*xpiy*z*(-3.0*xmiy*xpiy + 4.0*z2);
   complex_double dR5p2_xmiy = -105.0/2.0*xpiy3*z;
   RlmX[4*32+1] = dR5p2_xmiy + dR5p2_xpiy;
   RlmX[4*32+2] = -I*dR5p2_xmiy + I*dR5p2_xpiy;
   RlmX[4*32+3] = 210.0*xpiy2*z2 + (105.0/2.0)*xpiy2*(-xmiy*xpiy + 2.0*z2);
   RlmX[4*33+0] = (105.0/2.0)*xpiy3*(-xmiy*xpiy + 8.0*z2);
   complex_double dR5p3_xpiy = 210.0*xpiy2*(-xmiy*xpiy + 6.0*z2);
   complex_double dR5p3_xmiy = -105.0/2.0*xpiy4;
   RlmX[4*33+1] = dR5p3_xmiy + dR5p3_xpiy;
   RlmX[4*33+2] = -I*dR5p3_xmiy + I*dR5p3_xpiy;
   RlmX[4*33+3] = 840.0*xpiy3*z;
   RlmX[4*34+0] = 945.0*xpiy4*z;
   complex_double dR5p4_xpiy = 3780.0*xpiy3*z;
   RlmX[4*34+1] = dR5p4_xpiy;
   RlmX[4*34+2] = I*dR5p4_xpiy;
   RlmX[4*34+3] = 945.0*xpiy4;
   RlmX[4*35+0] = 945.0*xpiy5;
   complex_double dR5p5_xpiy = 4725.0*xpiy4;
   RlmX[4*35+1] = dR5p5_xpiy;
   RlmX[4*35+2] = I*dR5p5_xpiy;
   RlmX[4*35+3] = 0;
   if (5 == MaxL) return;
   complex_double xpiy6 = xpiy5*xpiy;
   complex_double xmiy6 = xmiy5*xmiy;
   double z6 = z5*z;
   RlmX[4*36+0] = (1.0/46080.0)*xmiy6;
   complex_double dR6m6_xmiy = (1.0/7680.0)*xmiy5;
   RlmX[4*36+1] = dR6m6_xmiy;
   RlmX[4*36+2] = -I*dR6m6_xmiy;
   RlmX[4*36+3] = 0;
   RlmX[4*37+0] = -1.0/3840.0*xmiy5*z;
   complex_double dR6m5_xmiy = -1.0/768.0*xmiy4*z;
   RlmX[4*37+1] = dR6m5_xmiy;
   RlmX[4*37+2] = -I*dR6m5_xmiy;
   RlmX[4*37+3] = -1.0/3840.0*xmiy5;
   RlmX[4*38+0] = (1.0/3840.0)*xmiy4*(-xmiy*xpiy + 10.0*z2);
   complex_double dR6m4_xpiy = -1.0/3840.0*xmiy5;
   complex_double dR6m4_xmiy = (1.0/768.0)*xmiy3*(-xmiy*xpiy + 8.0*z2);
   RlmX[4*38+1] = dR6m4_xmiy + dR6m4_xpiy;
   RlmX[4*38+2] = -I*dR6m4_xmiy + I*dR6m4_xpiy;
   RlmX[4*38+3] = (1.0/192.0)*xmiy4*z;
   RlmX[4*39+0] = (1.0/384.0)*xmiy3*z*(3.0*xmiy*xpiy - 8.0*z2);
   complex_double dR6m3_xpiy = (1.0/128.0)*xmiy4*z;
   complex_double dR6m3_xmiy = (1.0/32.0)*xmiy2*z*(xmiy*xpiy - 2.0*z2);
   RlmX[4*39+1] = dR6m3_xmiy + dR6m3_xpiy;
   RlmX[4*39+2] = -I*dR6m3_xmiy + I*dR6m3_xpiy;
   RlmX[4*39+3] = -1.0/24.0*xmiy3*z2 + (1.0/384.0)*xmiy3*(3.0*xmiy*xpiy - 8.0*z2);
   RlmX[4*40+0] = (1.0/128.0)*xmiy2*(xmiy2*xpiy2 - 16.0*xmiy*xpiy*z2 + 16.0*z4);
   complex_double dR6m2_xpiy = (1.0/64.0)*xmiy3*(xmiy*xpiy - 8.0*z2);
   complex_double dR6m2_xmiy = (1.0/32.0)*xmiy*(xmiy2*xpiy2 - 12.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[4*40+1] = dR6m2_xmiy + dR6m2_xpiy;
   RlmX[4*40+2] = -I*dR6m2_xmiy + I*dR6m2_xpiy;
   RlmX[4*40+3] = (1.0/128.0)*xmiy2*(-32.0*xmiy*xpiy*z + 64.0*z3);
   RlmX[4*41+0] = (1.0/16.0)*xmiy*z*(-5.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 8.0*z4);
   complex_double dR6m1_xpiy = (5.0/8.0)*xmiy2*z*(-xmiy*xpiy + 2.0*z2);
   complex_double dR6m1_xmiy = (1.0/16.0)*z*(-15.0*xmiy2*xpiy2 + 40.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[4*41+1] = dR6m1_xmiy + dR6m1_xpiy;
   RlmX[4*41+2] = -I*dR6m1_xmiy + I*dR6m1_xpiy;
   RlmX[4*41+3] = (1.0/16.0)*xmiy*z*(40.0*xmiy*xpiy*z - 32.0*z3) + (1.0/16.0)*xmiy*(-5.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[4*42+0] = -5.0/16.0*xmiy3*xpiy3 + (45.0/8.0)*xmiy2*xpiy2*z2 - 15.0/2.0*xmiy*xpiy*z4 + z6;
   complex_double dR6n0_xpiy = (15.0/16.0)*xmiy*(-xmiy2*xpiy2 + 12.0*xmiy*xpiy*z2 - 8.0*z4);
   complex_double dR6n0_xmiy = (15.0/16.0)*xpiy*(-xmiy2*xpiy2 + 12.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[4*42+1] = dR6n0_xmiy + dR6n0_xpiy;
   RlmX[4*42+2] = -I*dR6n0_xmiy + I*dR6n0_xpiy;
   RlmX[4*42+3] = (45.0/4.0)*xmiy2*xpiy2*z - 30.0*xmiy*xpiy*z3 + 6.0*z5;
   RlmX[4*43+0] = (21.0/8.0)*xpiy*z*(5.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dR6p1_xpiy = 21.0*z*((15.0/8.0)*xmiy2*xpiy2 - 5.0*xmiy*xpiy*z2 + z4);
   complex_double dR6p1_xmiy = (105.0/4.0)*xpiy2*z*(xmiy*xpiy - 2.0*z2);
   RlmX[4*43+1] = dR6p1_xmiy + dR6p1_xpiy;
   RlmX[4*43+2] = -I*dR6p1_xmiy + I*dR6p1_xpiy;
   RlmX[4*43+3] = (21.0/8.0)*xpiy*z*(-40.0*xmiy*xpiy*z + 32.0*z3) + (21.0/8.0)*xpiy*(5.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[4*44+0] = (105.0/8.0)*xpiy2*(xmiy2*xpiy2 - 16.0*xmiy*xpiy*z2 + 16.0*z4);
   complex_double dR6p2_xpiy = (105.0/2.0)*xpiy*(xmiy2*xpiy2 - 12.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dR6p2_xmiy = (105.0/4.0)*xpiy3*(xmiy*xpiy - 8.0*z2);
   RlmX[4*44+1] = dR6p2_xmiy + dR6p2_xpiy;
   RlmX[4*44+2] = -I*dR6p2_xmiy + I*dR6p2_xpiy;
   RlmX[4*44+3] = (105.0/8.0)*xpiy2*(-32.0*xmiy*xpiy*z + 64.0*z3);
   RlmX[4*45+0] = (315.0/2.0)*xpiy3*z*(-3.0*xmiy*xpiy + 8.0*z2);
   complex_double dR6p3_xpiy = 1890.0*xpiy2*z*(-xmiy*xpiy + 2.0*z2);
   complex_double dR6p3_xmiy = -945.0/2.0*xpiy4*z;
   RlmX[4*45+1] = dR6p3_xmiy + dR6p3_xpiy;
   RlmX[4*45+2] = -I*dR6p3_xmiy + I*dR6p3_xpiy;
   RlmX[4*45+3] = 2520.0*xpiy3*z2 + (315.0/2.0)*xpiy3*(-3.0*xmiy*xpiy + 8.0*z2);
   RlmX[4*46+0] = (945.0/2.0)*xpiy4*(-xmiy*xpiy + 10.0*z2);
   complex_double dR6p4_xpiy = (4725.0/2.0)*xpiy3*(-xmiy*xpiy + 8.0*z2);
   complex_double dR6p4_xmiy = -945.0/2.0*xpiy5;
   RlmX[4*46+1] = dR6p4_xmiy + dR6p4_xpiy;
   RlmX[4*46+2] = -I*dR6p4_xmiy + I*dR6p4_xpiy;
   RlmX[4*46+3] = 9450.0*xpiy4*z;
   RlmX[4*47+0] = 10395.0*xpiy5*z;
   complex_double dR6p5_xpiy = 51975.0*xpiy4*z;
   RlmX[4*47+1] = dR6p5_xpiy;
   RlmX[4*47+2] = I*dR6p5_xpiy;
   RlmX[4*47+3] = 10395.0*xpiy5;
   RlmX[4*48+0] = 10395.0*xpiy6;
   complex_double dR6p6_xpiy = 62370.0*xpiy5;
   RlmX[4*48+1] = dR6p6_xpiy;
   RlmX[4*48+2] = I*dR6p6_xpiy;
   RlmX[4*48+3] = 0;
   if (6 == MaxL) return;
   complex_double xpiy7 = xpiy6*xpiy;
   complex_double xmiy7 = xmiy6*xmiy;
   double z7 = z6*z;
   RlmX[4*49+0] = -1.0/645120.0*xmiy7;
   complex_double dR7m7_xmiy = -1.0/92160.0*xmiy6;
   RlmX[4*49+1] = dR7m7_xmiy;
   RlmX[4*49+2] = -I*dR7m7_xmiy;
   RlmX[4*49+3] = 0;
   RlmX[4*50+0] = (1.0/46080.0)*xmiy6*z;
   complex_double dR7m6_xmiy = (1.0/7680.0)*xmiy5*z;
   RlmX[4*50+1] = dR7m6_xmiy;
   RlmX[4*50+2] = -I*dR7m6_xmiy;
   RlmX[4*50+3] = (1.0/46080.0)*xmiy6;
   RlmX[4*51+0] = (1.0/46080.0)*xmiy5*(xmiy*xpiy - 12.0*z2);
   complex_double dR7m5_xpiy = (1.0/46080.0)*xmiy6;
   complex_double dR7m5_xmiy = (1.0/7680.0)*xmiy4*(xmiy*xpiy - 10.0*z2);
   RlmX[4*51+1] = dR7m5_xmiy + dR7m5_xpiy;
   RlmX[4*51+2] = -I*dR7m5_xmiy + I*dR7m5_xpiy;
   RlmX[4*51+3] = -1.0/1920.0*xmiy5*z;
   RlmX[4*52+0] = (1.0/3840.0)*xmiy4*z*(-3.0*xmiy*xpiy + 10.0*z2);
   complex_double dR7m4_xpiy = -1.0/1280.0*xmiy5*z;
   complex_double dR7m4_xmiy = (1.0/768.0)*xmiy3*z*(-3.0*xmiy*xpiy + 8.0*z2);
   RlmX[4*52+1] = dR7m4_xmiy + dR7m4_xpiy;
   RlmX[4*52+2] = -I*dR7m4_xmiy + I*dR7m4_xpiy;
   RlmX[4*52+3] = (1.0/192.0)*xmiy4*z2 + (1.0/3840.0)*xmiy4*(-3.0*xmiy*xpiy + 10.0*z2);
   RlmX[4*53+0] = xmiy3*(-1.0/1280.0*xmiy2*xpiy2 + (1.0/64.0)*xmiy*xpiy*z2 - 1.0/48.0*z4);
   complex_double dR7m3_xpiy = (1.0/640.0)*xmiy4*(-xmiy*xpiy + 10.0*z2);
   complex_double dR7m3_xmiy = (1.0/256.0)*xmiy2*(-xmiy2*xpiy2 + 16.0*xmiy*xpiy*z2 - 16.0*z4);
   RlmX[4*53+1] = dR7m3_xmiy + dR7m3_xpiy;
   RlmX[4*53+2] = -I*dR7m3_xmiy + I*dR7m3_xpiy;
   RlmX[4*53+3] = xmiy3*((1.0/32.0)*xmiy*xpiy*z - 1.0/12.0*z3);
   RlmX[4*54+0] = (1.0/384.0)*xmiy2*z*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   complex_double dR7m2_xpiy = (5.0/192.0)*xmiy3*z*(3.0*xmiy*xpiy - 8.0*z2);
   complex_double dR7m2_xmiy = (1.0/32.0)*xmiy*z*(5.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[4*54+1] = dR7m2_xmiy + dR7m2_xpiy;
   RlmX[4*54+2] = -I*dR7m2_xmiy + I*dR7m2_xpiy;
   RlmX[4*54+3] = (1.0/384.0)*xmiy2*z*(-160.0*xmiy*xpiy*z + 192.0*z3) + (1.0/384.0)*xmiy2*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   RlmX[4*55+0] = (1.0/128.0)*xmiy*(5.0*xmiy3*xpiy3 - 120.0*xmiy2*xpiy2*z2 + 240.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dR7m1_xpiy = (15.0/128.0)*xmiy2*(xmiy2*xpiy2 - 16.0*xmiy*xpiy*z2 + 16.0*z4);
   complex_double dR7m1_xmiy = (5.0/32.0)*xmiy3*xpiy3 - 45.0/16.0*xmiy2*xpiy2*z2 + (15.0/4.0)*xmiy*xpiy*z4 - 1.0/2.0*z6;
   RlmX[4*55+1] = dR7m1_xmiy + dR7m1_xpiy;
   RlmX[4*55+2] = -I*dR7m1_xmiy + I*dR7m1_xpiy;
   RlmX[4*55+3] = (1.0/128.0)*xmiy*(-240.0*xmiy2*xpiy2*z + 960.0*xmiy*xpiy*z3 - 384.0*z5);
   RlmX[4*56+0] = (1.0/16.0)*z*(-35.0*xmiy3*xpiy3 + 210.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 16.0*z6);
   complex_double dR7n0_xpiy = (21.0/16.0)*xmiy*z*(-5.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 8.0*z4);
   complex_double dR7n0_xmiy = (21.0/16.0)*xpiy*z*(-5.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[4*56+1] = dR7n0_xmiy + dR7n0_xpiy;
   RlmX[4*56+2] = -I*dR7n0_xmiy + I*dR7n0_xpiy;
   RlmX[4*56+3] = -35.0/16.0*xmiy3*xpiy3 + (105.0/8.0)*xmiy2*xpiy2*z2 - 21.0/2.0*xmiy*xpiy*z4 + z6 + (1.0/16.0)*z*(420.0*xmiy2*xpiy2*z - 672.0*xmiy*xpiy*z3 + 96.0*z5);
   RlmX[4*57+0] = (7.0/16.0)*xpiy*(-5.0*xmiy3*xpiy3 + 120.0*xmiy2*xpiy2*z2 - 240.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dR7p1_xpiy = -35.0/4.0*xmiy3*xpiy3 + (315.0/2.0)*xmiy2*xpiy2*z2 - 210.0*xmiy*xpiy*z4 + 28.0*z6;
   complex_double dR7p1_xmiy = (105.0/16.0)*xpiy2*(-xmiy2*xpiy2 + 16.0*xmiy*xpiy*z2 - 16.0*z4);
   RlmX[4*57+1] = dR7p1_xmiy + dR7p1_xpiy;
   RlmX[4*57+2] = -I*dR7p1_xmiy + I*dR7p1_xpiy;
   RlmX[4*57+3] = (7.0/16.0)*xpiy*(240.0*xmiy2*xpiy2*z - 960.0*xmiy*xpiy*z3 + 384.0*z5);
   RlmX[4*58+0] = (63.0/8.0)*xpiy2*z*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   complex_double dR7p2_xpiy = (189.0/2.0)*xpiy*z*(5.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dR7p2_xmiy = (315.0/4.0)*xpiy3*z*(3.0*xmiy*xpiy - 8.0*z2);
   RlmX[4*58+1] = dR7p2_xmiy + dR7p2_xpiy;
   RlmX[4*58+2] = -I*dR7p2_xmiy + I*dR7p2_xpiy;
   RlmX[4*58+3] = (63.0/8.0)*xpiy2*z*(-160.0*xmiy*xpiy*z + 192.0*z3) + (63.0/8.0)*xpiy2*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   RlmX[4*59+0] = (315.0/8.0)*xpiy3*(3.0*xmiy2*xpiy2 - 60.0*xmiy*xpiy*z2 + 80.0*z4);
   complex_double dR7p3_xpiy = (4725.0/8.0)*xpiy2*(xmiy2*xpiy2 - 16.0*xmiy*xpiy*z2 + 16.0*z4);
   complex_double dR7p3_xmiy = (945.0/4.0)*xpiy4*(xmiy*xpiy - 10.0*z2);
   RlmX[4*59+1] = dR7p3_xmiy + dR7p3_xpiy;
   RlmX[4*59+2] = -I*dR7p3_xmiy + I*dR7p3_xpiy;
   RlmX[4*59+3] = (315.0/8.0)*xpiy3*(-120.0*xmiy*xpiy*z + 320.0*z3);
   RlmX[4*60+0] = (3465.0/2.0)*xpiy4*z*(-3.0*xmiy*xpiy + 10.0*z2);
   complex_double dR7p4_xpiy = (17325.0/2.0)*xpiy3*z*(-3.0*xmiy*xpiy + 8.0*z2);
   complex_double dR7p4_xmiy = -10395.0/2.0*xpiy5*z;
   RlmX[4*60+1] = dR7p4_xmiy + dR7p4_xpiy;
   RlmX[4*60+2] = -I*dR7p4_xmiy + I*dR7p4_xpiy;
   RlmX[4*60+3] = 34650.0*xpiy4*z2 + (3465.0/2.0)*xpiy4*(-3.0*xmiy*xpiy + 10.0*z2);
   RlmX[4*61+0] = (10395.0/2.0)*xpiy5*(-xmiy*xpiy + 12.0*z2);
   complex_double dR7p5_xpiy = 31185.0*xpiy4*(-xmiy*xpiy + 10.0*z2);
   complex_double dR7p5_xmiy = -10395.0/2.0*xpiy6;
   RlmX[4*61+1] = dR7p5_xmiy + dR7p5_xpiy;
   RlmX[4*61+2] = -I*dR7p5_xmiy + I*dR7p5_xpiy;
   RlmX[4*61+3] = 124740.0*xpiy5*z;
   RlmX[4*62+0] = 135135.0*xpiy6*z;
   complex_double dR7p6_xpiy = 810810.0*xpiy5*z;
   RlmX[4*62+1] = dR7p6_xpiy;
   RlmX[4*62+2] = I*dR7p6_xpiy;
   RlmX[4*62+3] = 135135.0*xpiy6;
   RlmX[4*63+0] = 135135.0*xpiy7;
   complex_double dR7p7_xpiy = 945945.0*xpiy6;
   RlmX[4*63+1] = dR7p7_xpiy;
   RlmX[4*63+2] = I*dR7p7_xpiy;
   RlmX[4*63+3] = 0;
   if (7 == MaxL) return;
   complex_double xpiy8 = xpiy7*xpiy;
   complex_double xmiy8 = xmiy7*xmiy;
   double z8 = z7*z;
   RlmX[4*64+0] = (1.0/10321920.0)*xmiy8;
   complex_double dR8m8_xmiy = (1.0/1290240.0)*xmiy7;
   RlmX[4*64+1] = dR8m8_xmiy;
   RlmX[4*64+2] = -I*dR8m8_xmiy;
   RlmX[4*64+3] = 0;
   RlmX[4*65+0] = -1.0/645120.0*xmiy7*z;
   complex_double dR8m7_xmiy = -1.0/92160.0*xmiy6*z;
   RlmX[4*65+1] = dR8m7_xmiy;
   RlmX[4*65+2] = -I*dR8m7_xmiy;
   RlmX[4*65+3] = -1.0/645120.0*xmiy7;
   RlmX[4*66+0] = (1.0/645120.0)*xmiy6*(-xmiy*xpiy + 14.0*z2);
   complex_double dR8m6_xpiy = -1.0/645120.0*xmiy7;
   complex_double dR8m6_xmiy = (1.0/92160.0)*xmiy5*(-xmiy*xpiy + 12.0*z2);
   RlmX[4*66+1] = dR8m6_xmiy + dR8m6_xpiy;
   RlmX[4*66+2] = -I*dR8m6_xmiy + I*dR8m6_xpiy;
   RlmX[4*66+3] = (1.0/23040.0)*xmiy6*z;
   RlmX[4*67+0] = (1.0/15360.0)*xmiy5*z*(xmiy*xpiy - 4.0*z2);
   complex_double dR8m5_xpiy = (1.0/15360.0)*xmiy6*z;
   complex_double dR8m5_xmiy = (1.0/7680.0)*xmiy4*z*(3.0*xmiy*xpiy - 10.0*z2);
   RlmX[4*67+1] = dR8m5_xmiy + dR8m5_xpiy;
   RlmX[4*67+2] = -I*dR8m5_xmiy + I*dR8m5_xpiy;
   RlmX[4*67+3] = -1.0/1920.0*xmiy5*z2 + (1.0/15360.0)*xmiy5*(xmiy*xpiy - 4.0*z2);
   RlmX[4*68+0] = (1.0/15360.0)*xmiy4*(xmiy2*xpiy2 - 24.0*xmiy*xpiy*z2 + 40.0*z4);
   complex_double dR8m4_xpiy = (1.0/7680.0)*xmiy5*(xmiy*xpiy - 12.0*z2);
   complex_double dR8m4_xmiy = xmiy3*((1.0/2560.0)*xmiy2*xpiy2 - 1.0/128.0*xmiy*xpiy*z2 + (1.0/96.0)*z4);
   RlmX[4*68+1] = dR8m4_xmiy + dR8m4_xpiy;
   RlmX[4*68+2] = -I*dR8m4_xmiy + I*dR8m4_xpiy;
   RlmX[4*68+3] = (1.0/15360.0)*xmiy4*(-48.0*xmiy*xpiy*z + 160.0*z3);
   RlmX[4*69+0] = (1.0/768.0)*xmiy3*z*(-3.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 16.0*z4);
   complex_double dR8m3_xpiy = (1.0/384.0)*xmiy4*z*(-3.0*xmiy*xpiy + 10.0*z2);
   complex_double dR8m3_xmiy = (1.0/768.0)*xmiy2*z*(-15.0*xmiy2*xpiy2 + 80.0*xmiy*xpiy*z2 - 48.0*z4);
   RlmX[4*69+1] = dR8m3_xmiy + dR8m3_xpiy;
   RlmX[4*69+2] = -I*dR8m3_xmiy + I*dR8m3_xpiy;
   RlmX[4*69+3] = (1.0/768.0)*xmiy3*z*(40.0*xmiy*xpiy*z - 64.0*z3) + (1.0/768.0)*xmiy3*(-3.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 16.0*z4);
   RlmX[4*70+0] = (1.0/256.0)*xmiy2*(-xmiy3*xpiy3 + 30.0*xmiy2*xpiy2*z2 - 80.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dR8m2_xpiy = (1.0/256.0)*xmiy3*(-3.0*xmiy2*xpiy2 + 60.0*xmiy*xpiy*z2 - 80.0*z4);
   complex_double dR8m2_xmiy = (1.0/256.0)*xmiy*(-5.0*xmiy3*xpiy3 + 120.0*xmiy2*xpiy2*z2 - 240.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[4*70+1] = dR8m2_xmiy + dR8m2_xpiy;
   RlmX[4*70+2] = -I*dR8m2_xmiy + I*dR8m2_xpiy;
   RlmX[4*70+3] = (1.0/256.0)*xmiy2*(60.0*xmiy2*xpiy2*z - 320.0*xmiy*xpiy*z3 + 192.0*z5);
   RlmX[4*71+0] = (1.0/128.0)*xmiy*z*(35.0*xmiy3*xpiy3 - 280.0*xmiy2*xpiy2*z2 + 336.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dR8m1_xpiy = (7.0/128.0)*xmiy2*z*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   complex_double dR8m1_xmiy = (1.0/32.0)*z*(35.0*xmiy3*xpiy3 - 210.0*xmiy2*xpiy2*z2 + 168.0*xmiy*xpiy*z4 - 16.0*z6);
   RlmX[4*71+1] = dR8m1_xmiy + dR8m1_xpiy;
   RlmX[4*71+2] = -I*dR8m1_xmiy + I*dR8m1_xpiy;
   RlmX[4*71+3] = (1.0/128.0)*xmiy*z*(-560.0*xmiy2*xpiy2*z + 1344.0*xmiy*xpiy*z3 - 384.0*z5) + (1.0/128.0)*xmiy*(35.0*xmiy3*xpiy3 - 280.0*xmiy2*xpiy2*z2 + 336.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[4*72+0] = (35.0/128.0)*xmiy4*xpiy4 - 35.0/4.0*xmiy3*xpiy3*z2 + (105.0/4.0)*xmiy2*xpiy2*z4 - 14.0*xmiy*xpiy*z6 + z8;
   complex_double dR8n0_xpiy = (7.0/32.0)*xmiy*(5.0*xmiy3*xpiy3 - 120.0*xmiy2*xpiy2*z2 + 240.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dR8n0_xmiy = (7.0/32.0)*xpiy*(5.0*xmiy3*xpiy3 - 120.0*xmiy2*xpiy2*z2 + 240.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[4*72+1] = dR8n0_xmiy + dR8n0_xpiy;
   RlmX[4*72+2] = -I*dR8n0_xmiy + I*dR8n0_xpiy;
   RlmX[4*72+3] = -35.0/2.0*xmiy3*xpiy3*z + 105.0*xmiy2*xpiy2*z3 - 84.0*xmiy*xpiy*z5 + 8.0*z7;
   RlmX[4*73+0] = (9.0/16.0)*xpiy*z*(-35.0*xmiy3*xpiy3 + 280.0*xmiy2*xpiy2*z2 - 336.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dR8p1_xpiy = (9.0/4.0)*z*(-35.0*xmiy3*xpiy3 + 210.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 16.0*z6);
   complex_double dR8p1_xmiy = (63.0/16.0)*xpiy2*z*(-15.0*xmiy2*xpiy2 + 80.0*xmiy*xpiy*z2 - 48.0*z4);
   RlmX[4*73+1] = dR8p1_xmiy + dR8p1_xpiy;
   RlmX[4*73+2] = -I*dR8p1_xmiy + I*dR8p1_xpiy;
   RlmX[4*73+3] = (9.0/16.0)*xpiy*z*(560.0*xmiy2*xpiy2*z - 1344.0*xmiy*xpiy*z3 + 384.0*z5) + (9.0/16.0)*xpiy*(-35.0*xmiy3*xpiy3 + 280.0*xmiy2*xpiy2*z2 - 336.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[4*74+0] = (315.0/16.0)*xpiy2*(-xmiy3*xpiy3 + 30.0*xmiy2*xpiy2*z2 - 80.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dR8p2_xpiy = (315.0/16.0)*xpiy*(-5.0*xmiy3*xpiy3 + 120.0*xmiy2*xpiy2*z2 - 240.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dR8p2_xmiy = (315.0/16.0)*xpiy3*(-3.0*xmiy2*xpiy2 + 60.0*xmiy*xpiy*z2 - 80.0*z4);
   RlmX[4*74+1] = dR8p2_xmiy + dR8p2_xpiy;
   RlmX[4*74+2] = -I*dR8p2_xmiy + I*dR8p2_xpiy;
   RlmX[4*74+3] = (315.0/16.0)*xpiy2*(60.0*xmiy2*xpiy2*z - 320.0*xmiy*xpiy*z3 + 192.0*z5);
   RlmX[4*75+0] = (3465.0/8.0)*xpiy3*z*(3.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 16.0*z4);
   complex_double dR8p3_xpiy = (3465.0/8.0)*xpiy2*z*(15.0*xmiy2*xpiy2 - 80.0*xmiy*xpiy*z2 + 48.0*z4);
   complex_double dR8p3_xmiy = (3465.0/4.0)*xpiy4*z*(3.0*xmiy*xpiy - 10.0*z2);
   RlmX[4*75+1] = dR8p3_xmiy + dR8p3_xpiy;
   RlmX[4*75+2] = -I*dR8p3_xmiy + I*dR8p3_xpiy;
   RlmX[4*75+3] = (3465.0/8.0)*xpiy3*z*(-40.0*xmiy*xpiy*z + 64.0*z3) + (3465.0/8.0)*xpiy3*(3.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 16.0*z4);
   RlmX[4*76+0] = (10395.0/8.0)*xpiy4*(xmiy2*xpiy2 - 24.0*xmiy*xpiy*z2 + 40.0*z4);
   complex_double dR8p4_xpiy = xpiy3*((31185.0/4.0)*xmiy2*xpiy2 - 155925.0*xmiy*xpiy*z2 + 207900.0*z4);
   complex_double dR8p4_xmiy = (10395.0/4.0)*xpiy5*(xmiy*xpiy - 12.0*z2);
   RlmX[4*76+1] = dR8p4_xmiy + dR8p4_xpiy;
   RlmX[4*76+2] = -I*dR8p4_xmiy + I*dR8p4_xpiy;
   RlmX[4*76+3] = (10395.0/8.0)*xpiy4*(-48.0*xmiy*xpiy*z + 160.0*z3);
   RlmX[4*77+0] = (135135.0/2.0)*xpiy5*z*(-xmiy*xpiy + 4.0*z2);
   complex_double dR8p5_xpiy = 135135.0*xpiy4*z*(-3.0*xmiy*xpiy + 10.0*z2);
   complex_double dR8p5_xmiy = -135135.0/2.0*xpiy6*z;
   RlmX[4*77+1] = dR8p5_xmiy + dR8p5_xpiy;
   RlmX[4*77+2] = -I*dR8p5_xmiy + I*dR8p5_xpiy;
   RlmX[4*77+3] = 540540.0*xpiy5*z2 + (135135.0/2.0)*xpiy5*(-xmiy*xpiy + 4.0*z2);
   RlmX[4*78+0] = (135135.0/2.0)*xpiy6*(-xmiy*xpiy + 14.0*z2);
   complex_double dR8p6_xpiy = (945945.0/2.0)*xpiy5*(-xmiy*xpiy + 12.0*z2);
   complex_double dR8p6_xmiy = -135135.0/2.0*xpiy7;
   RlmX[4*78+1] = dR8p6_xmiy + dR8p6_xpiy;
   RlmX[4*78+2] = -I*dR8p6_xmiy + I*dR8p6_xpiy;
   RlmX[4*78+3] = 1891890.0*xpiy6*z;
   RlmX[4*79+0] = 2027025.0*xpiy7*z;
   complex_double dR8p7_xpiy = 14189175.0*xpiy6*z;
   RlmX[4*79+1] = dR8p7_xpiy;
   RlmX[4*79+2] = I*dR8p7_xpiy;
   RlmX[4*79+3] = 2027025.0*xpiy7;
   RlmX[4*80+0] = 2027025.0*xpiy8;
   complex_double dR8p8_xpiy = 16216200.0*xpiy7;
   RlmX[4*80+1] = dR8p8_xpiy;
   RlmX[4*80+2] = I*dR8p8_xpiy;
   RlmX[4*80+3] = 0;
   if (8 == MaxL) return;
   complex_double xpiy9 = xpiy8*xpiy;
   complex_double xmiy9 = xmiy8*xmiy;
   double z9 = z8*z;
   RlmX[4*81+0] = -1.0/185794560.0*xmiy9;
   complex_double dR9m9_xmiy = -1.0/20643840.0*xmiy8;
   RlmX[4*81+1] = dR9m9_xmiy;
   RlmX[4*81+2] = -I*dR9m9_xmiy;
   RlmX[4*81+3] = 0;
   RlmX[4*82+0] = (1.0/10321920.0)*xmiy8*z;
   complex_double dR9m8_xmiy = (1.0/1290240.0)*xmiy7*z;
   RlmX[4*82+1] = dR9m8_xmiy;
   RlmX[4*82+2] = -I*dR9m8_xmiy;
   RlmX[4*82+3] = (1.0/10321920.0)*xmiy8;
   RlmX[4*83+0] = (1.0/10321920.0)*xmiy7*(xmiy*xpiy - 16.0*z2);
   complex_double dR9m7_xpiy = (1.0/10321920.0)*xmiy8;
   complex_double dR9m7_xmiy = (1.0/1290240.0)*xmiy6*(xmiy*xpiy - 14.0*z2);
   RlmX[4*83+1] = dR9m7_xmiy + dR9m7_xpiy;
   RlmX[4*83+2] = -I*dR9m7_xmiy + I*dR9m7_xpiy;
   RlmX[4*83+3] = -1.0/322560.0*xmiy7*z;
   RlmX[4*84+0] = (1.0/645120.0)*xmiy6*z*(-3.0*xmiy*xpiy + 14.0*z2);
   complex_double dR9m6_xpiy = -1.0/215040.0*xmiy7*z;
   complex_double dR9m6_xmiy = (1.0/30720.0)*xmiy5*z*(-xmiy*xpiy + 4.0*z2);
   RlmX[4*84+1] = dR9m6_xmiy + dR9m6_xpiy;
   RlmX[4*84+2] = -I*dR9m6_xmiy + I*dR9m6_xpiy;
   RlmX[4*84+3] = (1.0/23040.0)*xmiy6*z2 + (1.0/645120.0)*xmiy6*(-3.0*xmiy*xpiy + 14.0*z2);
   RlmX[4*85+0] = (1.0/215040.0)*xmiy5*(-xmiy2*xpiy2 + 28.0*xmiy*xpiy*z2 - 56.0*z4);
   complex_double dR9m5_xpiy = (1.0/107520.0)*xmiy6*(-xmiy*xpiy + 14.0*z2);
   complex_double dR9m5_xmiy = (1.0/30720.0)*xmiy4*(-xmiy2*xpiy2 + 24.0*xmiy*xpiy*z2 - 40.0*z4);
   RlmX[4*85+1] = dR9m5_xmiy + dR9m5_xpiy;
   RlmX[4*85+2] = -I*dR9m5_xmiy + I*dR9m5_xpiy;
   RlmX[4*85+3] = (1.0/215040.0)*xmiy5*(56.0*xmiy*xpiy*z - 224.0*z3);
   RlmX[4*86+0] = (1.0/3072.0)*xmiy4*z*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dR9m4_xpiy = (1.0/1536.0)*xmiy5*z*(xmiy*xpiy - 4.0*z2);
   complex_double dR9m4_xmiy = (1.0/1536.0)*xmiy3*z*(3.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 16.0*z4);
   RlmX[4*86+1] = dR9m4_xmiy + dR9m4_xpiy;
   RlmX[4*86+2] = -I*dR9m4_xmiy + I*dR9m4_xpiy;
   RlmX[4*86+3] = (1.0/3072.0)*xmiy4*z*(-16.0*xmiy*xpiy*z + 32.0*z3) + (1.0/3072.0)*xmiy4*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[4*87+0] = (1.0/3072.0)*xmiy3*(xmiy3*xpiy3 - 36.0*xmiy2*xpiy2*z2 + 120.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dR9m3_xpiy = (1.0/1024.0)*xmiy4*(xmiy2*xpiy2 - 24.0*xmiy*xpiy*z2 + 40.0*z4);
   complex_double dR9m3_xmiy = (1.0/512.0)*xmiy2*(xmiy3*xpiy3 - 30.0*xmiy2*xpiy2*z2 + 80.0*xmiy*xpiy*z4 - 32.0*z6);
   RlmX[4*87+1] = dR9m3_xmiy + dR9m3_xpiy;
   RlmX[4*87+2] = -I*dR9m3_xmiy + I*dR9m3_xpiy;
   RlmX[4*87+3] = (1.0/3072.0)*xmiy3*(-72.0*xmiy2*xpiy2*z + 480.0*xmiy*xpiy*z3 - 384.0*z5);
   RlmX[4*88+0] = (1.0/256.0)*xmiy2*z*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dR9m2_xpiy = (7.0/256.0)*xmiy3*z*(-3.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 16.0*z4);
   complex_double dR9m2_xmiy = (1.0/256.0)*xmiy*z*(-35.0*xmiy3*xpiy3 + 280.0*xmiy2*xpiy2*z2 - 336.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[4*88+1] = dR9m2_xmiy + dR9m2_xpiy;
   RlmX[4*88+2] = -I*dR9m2_xmiy + I*dR9m2_xpiy;
   RlmX[4*88+3] = (1.0/256.0)*xmiy2*z*(140.0*xmiy2*xpiy2*z - 448.0*xmiy*xpiy*z3 + 192.0*z5) + (1.0/256.0)*xmiy2*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   RlmX[4*89+0] = (1.0/256.0)*xmiy*(-7.0*xmiy4*xpiy4 + 280.0*xmiy3*xpiy3*z2 - 1120.0*xmiy2*xpiy2*z4 + 896.0*xmiy*xpiy*z6 - 128.0*z8);
   complex_double dR9m1_xpiy = (7.0/64.0)*xmiy2*(-xmiy3*xpiy3 + 30.0*xmiy2*xpiy2*z2 - 80.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dR9m1_xmiy = -35.0/256.0*xmiy4*xpiy4 + (35.0/8.0)*xmiy3*xpiy3*z2 - 105.0/8.0*xmiy2*xpiy2*z4 + 7.0*xmiy*xpiy*z6 - 1.0/2.0*z8;
   RlmX[4*89+1] = dR9m1_xmiy + dR9m1_xpiy;
   RlmX[4*89+2] = -I*dR9m1_xmiy + I*dR9m1_xpiy;
   RlmX[4*89+3] = (1.0/256.0)*xmiy*(560.0*xmiy3*xpiy3*z - 4480.0*xmiy2*xpiy2*z3 + 5376.0*xmiy*xpiy*z5 - 1024.0*z7);
   RlmX[4*90+0] = (1.0/128.0)*z*(315.0*xmiy4*xpiy4 - 3360.0*xmiy3*xpiy3*z2 + 6048.0*xmiy2*xpiy2*z4 - 2304.0*xmiy*xpiy*z6 + 128.0*z8);
   complex_double dR9n0_xpiy = (9.0/32.0)*xmiy*z*(35.0*xmiy3*xpiy3 - 280.0*xmiy2*xpiy2*z2 + 336.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dR9n0_xmiy = (9.0/32.0)*xpiy*z*(35.0*xmiy3*xpiy3 - 280.0*xmiy2*xpiy2*z2 + 336.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[4*90+1] = dR9n0_xmiy + dR9n0_xpiy;
   RlmX[4*90+2] = -I*dR9n0_xmiy + I*dR9n0_xpiy;
   RlmX[4*90+3] = (315.0/128.0)*xmiy4*xpiy4 - 105.0/4.0*xmiy3*xpiy3*z2 + (189.0/4.0)*xmiy2*xpiy2*z4 - 18.0*xmiy*xpiy*z6 + z8 + (1.0/128.0)*z*(-6720.0*xmiy3*xpiy3*z + 24192.0*xmiy2*xpiy2*z3 - 13824.0*xmiy*xpiy*z5 + 1024.0*z7);
   RlmX[4*91+0] = (45.0/128.0)*xpiy*(7.0*xmiy4*xpiy4 - 280.0*xmiy3*xpiy3*z2 + 1120.0*xmiy2*xpiy2*z4 - 896.0*xmiy*xpiy*z6 + 128.0*z8);
   complex_double dR9p1_xpiy = (1575.0/128.0)*xmiy4*xpiy4 - 1575.0/4.0*xmiy3*xpiy3*z2 + (4725.0/4.0)*xmiy2*xpiy2*z4 - 630.0*xmiy*xpiy*z6 + 45.0*z8;
   complex_double dR9p1_xmiy = (315.0/32.0)*xpiy2*(xmiy3*xpiy3 - 30.0*xmiy2*xpiy2*z2 + 80.0*xmiy*xpiy*z4 - 32.0*z6);
   RlmX[4*91+1] = dR9p1_xmiy + dR9p1_xpiy;
   RlmX[4*91+2] = -I*dR9p1_xmiy + I*dR9p1_xpiy;
   RlmX[4*91+3] = (45.0/128.0)*xpiy*(-560.0*xmiy3*xpiy3*z + 4480.0*xmiy2*xpiy2*z3 - 5376.0*xmiy*xpiy*z5 + 1024.0*z7);
   RlmX[4*92+0] = (495.0/16.0)*xpiy2*z*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dR9p2_xpiy = (495.0/16.0)*xpiy*z*(-35.0*xmiy3*xpiy3 + 280.0*xmiy2*xpiy2*z2 - 336.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dR9p2_xmiy = (3465.0/16.0)*xpiy3*z*(-3.0*xmiy2*xpiy2 + 20.0*xmiy*xpiy*z2 - 16.0*z4);
   RlmX[4*92+1] = dR9p2_xmiy + dR9p2_xpiy;
   RlmX[4*92+2] = -I*dR9p2_xmiy + I*dR9p2_xpiy;
   RlmX[4*92+3] = (495.0/16.0)*xpiy2*z*(140.0*xmiy2*xpiy2*z - 448.0*xmiy*xpiy*z3 + 192.0*z5) + (495.0/16.0)*xpiy2*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   RlmX[4*93+0] = (3465.0/16.0)*xpiy3*(-xmiy3*xpiy3 + 36.0*xmiy2*xpiy2*z2 - 120.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dR9p3_xpiy = (10395.0/8.0)*xpiy2*(-xmiy3*xpiy3 + 30.0*xmiy2*xpiy2*z2 - 80.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dR9p3_xmiy = (10395.0/16.0)*xpiy4*(-xmiy2*xpiy2 + 24.0*xmiy*xpiy*z2 - 40.0*z4);
   RlmX[4*93+1] = dR9p3_xmiy + dR9p3_xpiy;
   RlmX[4*93+2] = -I*dR9p3_xmiy + I*dR9p3_xpiy;
   RlmX[4*93+3] = (3465.0/16.0)*xpiy3*(72.0*xmiy2*xpiy2*z - 480.0*xmiy*xpiy*z3 + 384.0*z5);
   RlmX[4*94+0] = (135135.0/8.0)*xpiy4*z*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dR9p4_xpiy = (135135.0/4.0)*xpiy3*z*(3.0*xmiy2*xpiy2 - 20.0*xmiy*xpiy*z2 + 16.0*z4);
   complex_double dR9p4_xmiy = (135135.0/4.0)*xpiy5*z*(xmiy*xpiy - 4.0*z2);
   RlmX[4*94+1] = dR9p4_xmiy + dR9p4_xpiy;
   RlmX[4*94+2] = -I*dR9p4_xmiy + I*dR9p4_xpiy;
   RlmX[4*94+3] = (135135.0/8.0)*xpiy4*z*(-16.0*xmiy*xpiy*z + 32.0*z3) + (135135.0/8.0)*xpiy4*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   RlmX[4*95+0] = (135135.0/8.0)*xpiy5*(xmiy2*xpiy2 - 28.0*xmiy*xpiy*z2 + 56.0*z4);
   complex_double dR9p5_xpiy = (945945.0/8.0)*xpiy4*(xmiy2*xpiy2 - 24.0*xmiy*xpiy*z2 + 40.0*z4);
   complex_double dR9p5_xmiy = (135135.0/4.0)*xpiy6*(xmiy*xpiy - 14.0*z2);
   RlmX[4*95+1] = dR9p5_xmiy + dR9p5_xpiy;
   RlmX[4*95+2] = -I*dR9p5_xmiy + I*dR9p5_xpiy;
   RlmX[4*95+3] = (135135.0/8.0)*xpiy5*(-56.0*xmiy*xpiy*z + 224.0*z3);
   RlmX[4*96+0] = (675675.0/2.0)*xpiy6*z*(-3.0*xmiy*xpiy + 14.0*z2);
   complex_double dR9p6_xpiy = (14189175.0/2.0)*xpiy5*z*(-xmiy*xpiy + 4.0*z2);
   complex_double dR9p6_xmiy = -2027025.0/2.0*xpiy7*z;
   RlmX[4*96+1] = dR9p6_xmiy + dR9p6_xpiy;
   RlmX[4*96+2] = -I*dR9p6_xmiy + I*dR9p6_xpiy;
   RlmX[4*96+3] = 9459450.0*xpiy6*z2 + (675675.0/2.0)*xpiy6*(-3.0*xmiy*xpiy + 14.0*z2);
   RlmX[4*97+0] = (2027025.0/2.0)*xpiy7*(-xmiy*xpiy + 16.0*z2);
   complex_double dR9p7_xpiy = 8108100.0*xpiy6*(-xmiy*xpiy + 14.0*z2);
   complex_double dR9p7_xmiy = -2027025.0/2.0*xpiy8;
   RlmX[4*97+1] = dR9p7_xmiy + dR9p7_xpiy;
   RlmX[4*97+2] = -I*dR9p7_xmiy + I*dR9p7_xpiy;
   RlmX[4*97+3] = 32432400.0*xpiy7*z;
   RlmX[4*98+0] = 34459425.0*xpiy8*z;
   complex_double dR9p8_xpiy = 275675400.0*xpiy7*z;
   RlmX[4*98+1] = dR9p8_xpiy;
   RlmX[4*98+2] = I*dR9p8_xpiy;
   RlmX[4*98+3] = 34459425.0*xpiy8;
   RlmX[4*99+0] = 34459425.0*xpiy9;
   complex_double dR9p9_xpiy = 310134825.0*xpiy8;
   RlmX[4*99+1] = dR9p9_xpiy;
   RlmX[4*99+2] = I*dR9p9_xpiy;
   RlmX[4*99+3] = 0;
   if (9 == MaxL) return;
   complex_double xpiy10 = xpiy9*xpiy;
   complex_double xmiy10 = xmiy9*xmiy;
   double z10 = z9*z;
   RlmX[4*100+0] = (1.0/3715891200.0)*xmiy10;
   complex_double dRama_xmiy = (1.0/371589120.0)*xmiy9;
   RlmX[4*100+1] = dRama_xmiy;
   RlmX[4*100+2] = -I*dRama_xmiy;
   RlmX[4*100+3] = 0;
   RlmX[4*101+0] = -1.0/185794560.0*xmiy9*z;
   complex_double dRam9_xmiy = -1.0/20643840.0*xmiy8*z;
   RlmX[4*101+1] = dRam9_xmiy;
   RlmX[4*101+2] = -I*dRam9_xmiy;
   RlmX[4*101+3] = -1.0/185794560.0*xmiy9;
   RlmX[4*102+0] = (1.0/185794560.0)*xmiy8*(-xmiy*xpiy + 18.0*z2);
   complex_double dRam8_xpiy = -1.0/185794560.0*xmiy9;
   complex_double dRam8_xmiy = (1.0/20643840.0)*xmiy7*(-xmiy*xpiy + 16.0*z2);
   RlmX[4*102+1] = dRam8_xmiy + dRam8_xpiy;
   RlmX[4*102+2] = -I*dRam8_xmiy + I*dRam8_xpiy;
   RlmX[4*102+3] = (1.0/5160960.0)*xmiy8*z;
   RlmX[4*103+0] = (1.0/10321920.0)*xmiy7*z*(3.0*xmiy*xpiy - 16.0*z2);
   complex_double dRam7_xpiy = (1.0/3440640.0)*xmiy8*z;
   complex_double dRam7_xmiy = (1.0/1290240.0)*xmiy6*z*(3.0*xmiy*xpiy - 14.0*z2);
   RlmX[4*103+1] = dRam7_xmiy + dRam7_xpiy;
   RlmX[4*103+2] = -I*dRam7_xmiy + I*dRam7_xpiy;
   RlmX[4*103+3] = -1.0/322560.0*xmiy7*z2 + (1.0/10321920.0)*xmiy7*(3.0*xmiy*xpiy - 16.0*z2);
   RlmX[4*104+0] = xmiy6*((1.0/3440640.0)*xmiy2*xpiy2 - 1.0/107520.0*xmiy*xpiy*z2 + (1.0/46080.0)*z4);
   complex_double dRam6_xpiy = (1.0/1720320.0)*xmiy7*(xmiy*xpiy - 16.0*z2);
   complex_double dRam6_xmiy = (1.0/430080.0)*xmiy5*(xmiy2*xpiy2 - 28.0*xmiy*xpiy*z2 + 56.0*z4);
   RlmX[4*104+1] = dRam6_xmiy + dRam6_xpiy;
   RlmX[4*104+2] = -I*dRam6_xmiy + I*dRam6_xpiy;
   RlmX[4*104+3] = xmiy6*(-1.0/53760.0*xmiy*xpiy*z + (1.0/11520.0)*z3);
   RlmX[4*105+0] = (1.0/645120.0)*xmiy5*z*(-15.0*xmiy2*xpiy2 + 140.0*xmiy*xpiy*z2 - 168.0*z4);
   complex_double dRam5_xpiy = (1.0/64512.0)*xmiy6*z*(-3.0*xmiy*xpiy + 14.0*z2);
   complex_double dRam5_xmiy = (1.0/6144.0)*xmiy4*z*(-xmiy2*xpiy2 + 8.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[4*105+1] = dRam5_xmiy + dRam5_xpiy;
   RlmX[4*105+2] = -I*dRam5_xmiy + I*dRam5_xpiy;
   RlmX[4*105+3] = (1.0/645120.0)*xmiy5*z*(280.0*xmiy*xpiy*z - 672.0*z3) + (1.0/645120.0)*xmiy5*(-15.0*xmiy2*xpiy2 + 140.0*xmiy*xpiy*z2 - 168.0*z4);
   RlmX[4*106+0] = (1.0/43008.0)*xmiy4*(-xmiy3*xpiy3 + 42.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 112.0*z6);
   complex_double dRam4_xpiy = (1.0/14336.0)*xmiy5*(-xmiy2*xpiy2 + 28.0*xmiy*xpiy*z2 - 56.0*z4);
   complex_double dRam4_xmiy = (1.0/6144.0)*xmiy3*(-xmiy3*xpiy3 + 36.0*xmiy2*xpiy2*z2 - 120.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[4*106+1] = dRam4_xmiy + dRam4_xpiy;
   RlmX[4*106+2] = -I*dRam4_xmiy + I*dRam4_xpiy;
   RlmX[4*106+3] = (1.0/43008.0)*xmiy4*(84.0*xmiy2*xpiy2*z - 672.0*xmiy*xpiy*z3 + 672.0*z5);
   RlmX[4*107+0] = (1.0/3072.0)*xmiy3*z*(7.0*xmiy3*xpiy3 - 84.0*xmiy2*xpiy2*z2 + 168.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dRam3_xpiy = (7.0/1024.0)*xmiy4*z*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dRam3_xmiy = (1.0/512.0)*xmiy2*z*(7.0*xmiy3*xpiy3 - 70.0*xmiy2*xpiy2*z2 + 112.0*xmiy*xpiy*z4 - 32.0*z6);
   RlmX[4*107+1] = dRam3_xmiy + dRam3_xpiy;
   RlmX[4*107+2] = -I*dRam3_xmiy + I*dRam3_xpiy;
   RlmX[4*107+3] = (1.0/3072.0)*xmiy3*z*(-168.0*xmiy2*xpiy2*z + 672.0*xmiy*xpiy*z3 - 384.0*z5) + (1.0/3072.0)*xmiy3*(7.0*xmiy3*xpiy3 - 84.0*xmiy2*xpiy2*z2 + 168.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[4*108+0] = (1.0/3072.0)*xmiy2*(7.0*xmiy4*xpiy4 - 336.0*xmiy3*xpiy3*z2 + 1680.0*xmiy2*xpiy2*z4 - 1792.0*xmiy*xpiy*z6 + 384.0*z8);
   complex_double dRam2_xpiy = (7.0/768.0)*xmiy3*(xmiy3*xpiy3 - 36.0*xmiy2*xpiy2*z2 + 120.0*xmiy*xpiy*z4 - 64.0*z6);
   complex_double dRam2_xmiy = (1.0/512.0)*xmiy*(7.0*xmiy4*xpiy4 - 280.0*xmiy3*xpiy3*z2 + 1120.0*xmiy2*xpiy2*z4 - 896.0*xmiy*xpiy*z6 + 128.0*z8);
   RlmX[4*108+1] = dRam2_xmiy + dRam2_xpiy;
   RlmX[4*108+2] = -I*dRam2_xmiy + I*dRam2_xpiy;
   RlmX[4*108+3] = (1.0/3072.0)*xmiy2*(-672.0*xmiy3*xpiy3*z + 6720.0*xmiy2*xpiy2*z3 - 10752.0*xmiy*xpiy*z5 + 3072.0*z7);
   RlmX[4*109+0] = (1.0/256.0)*xmiy*z*(-63.0*xmiy4*xpiy4 + 840.0*xmiy3*xpiy3*z2 - 2016.0*xmiy2*xpiy2*z4 + 1152.0*xmiy*xpiy*z6 - 128.0*z8);
   complex_double dRam1_xpiy = (9.0/64.0)*xmiy2*z*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dRam1_xmiy = (1.0/256.0)*z*(-315.0*xmiy4*xpiy4 + 3360.0*xmiy3*xpiy3*z2 - 6048.0*xmiy2*xpiy2*z4 + 2304.0*xmiy*xpiy*z6 - 128.0*z8);
   RlmX[4*109+1] = dRam1_xmiy + dRam1_xpiy;
   RlmX[4*109+2] = -I*dRam1_xmiy + I*dRam1_xpiy;
   RlmX[4*109+3] = (1.0/256.0)*xmiy*z*(1680.0*xmiy3*xpiy3*z - 8064.0*xmiy2*xpiy2*z3 + 6912.0*xmiy*xpiy*z5 - 1024.0*z7) + (1.0/256.0)*xmiy*(-63.0*xmiy4*xpiy4 + 840.0*xmiy3*xpiy3*z2 - 2016.0*xmiy2*xpiy2*z4 + 1152.0*xmiy*xpiy*z6 - 128.0*z8);
   RlmX[4*110+0] = -63.0/256.0*xmiy5*xpiy5 + (1575.0/128.0)*xmiy4*xpiy4*z2 - 525.0/8.0*xmiy3*xpiy3*z4 + (315.0/4.0)*xmiy2*xpiy2*z6 - 45.0/2.0*xmiy*xpiy*z8 + z10;
   complex_double dRan0_xpiy = (45.0/256.0)*xmiy*(-7.0*xmiy4*xpiy4 + 280.0*xmiy3*xpiy3*z2 - 1120.0*xmiy2*xpiy2*z4 + 896.0*xmiy*xpiy*z6 - 128.0*z8);
   complex_double dRan0_xmiy = (45.0/256.0)*xpiy*(-7.0*xmiy4*xpiy4 + 280.0*xmiy3*xpiy3*z2 - 1120.0*xmiy2*xpiy2*z4 + 896.0*xmiy*xpiy*z6 - 128.0*z8);
   RlmX[4*110+1] = dRan0_xmiy + dRan0_xpiy;
   RlmX[4*110+2] = -I*dRan0_xmiy + I*dRan0_xpiy;
   RlmX[4*110+3] = (1575.0/64.0)*xmiy4*xpiy4*z - 525.0/2.0*xmiy3*xpiy3*z3 + (945.0/2.0)*xmiy2*xpiy2*z5 - 180.0*xmiy*xpiy*z7 + 10.0*z9;
   RlmX[4*111+0] = (55.0/128.0)*xpiy*z*(63.0*xmiy4*xpiy4 - 840.0*xmiy3*xpiy3*z2 + 2016.0*xmiy2*xpiy2*z4 - 1152.0*xmiy*xpiy*z6 + 128.0*z8);
   complex_double dRap1_xpiy = (55.0/128.0)*z*(315.0*xmiy4*xpiy4 - 3360.0*xmiy3*xpiy3*z2 + 6048.0*xmiy2*xpiy2*z4 - 2304.0*xmiy*xpiy*z6 + 128.0*z8);
   complex_double dRap1_xmiy = (495.0/32.0)*xpiy2*z*(7.0*xmiy3*xpiy3 - 70.0*xmiy2*xpiy2*z2 + 112.0*xmiy*xpiy*z4 - 32.0*z6);
   RlmX[4*111+1] = dRap1_xmiy + dRap1_xpiy;
   RlmX[4*111+2] = -I*dRap1_xmiy + I*dRap1_xpiy;
   RlmX[4*111+3] = (55.0/128.0)*xpiy*z*(-1680.0*xmiy3*xpiy3*z + 8064.0*xmiy2*xpiy2*z3 - 6912.0*xmiy*xpiy*z5 + 1024.0*z7) + (55.0/128.0)*xpiy*(63.0*xmiy4*xpiy4 - 840.0*xmiy3*xpiy3*z2 + 2016.0*xmiy2*xpiy2*z4 - 1152.0*xmiy*xpiy*z6 + 128.0*z8);
   RlmX[4*112+0] = (495.0/128.0)*xpiy2*(7.0*xmiy4*xpiy4 - 336.0*xmiy3*xpiy3*z2 + 1680.0*xmiy2*xpiy2*z4 - 1792.0*xmiy*xpiy*z6 + 384.0*z8);
   complex_double dRap2_xpiy = (1485.0/64.0)*xpiy*(7.0*xmiy4*xpiy4 - 280.0*xmiy3*xpiy3*z2 + 1120.0*xmiy2*xpiy2*z4 - 896.0*xmiy*xpiy*z6 + 128.0*z8);
   complex_double dRap2_xmiy = (3465.0/32.0)*xpiy3*(xmiy3*xpiy3 - 36.0*xmiy2*xpiy2*z2 + 120.0*xmiy*xpiy*z4 - 64.0*z6);
   RlmX[4*112+1] = dRap2_xmiy + dRap2_xpiy;
   RlmX[4*112+2] = -I*dRap2_xmiy + I*dRap2_xpiy;
   RlmX[4*112+3] = (495.0/128.0)*xpiy2*(-672.0*xmiy3*xpiy3*z + 6720.0*xmiy2*xpiy2*z3 - 10752.0*xmiy*xpiy*z5 + 3072.0*z7);
   RlmX[4*113+0] = (6435.0/16.0)*xpiy3*z*(-7.0*xmiy3*xpiy3 + 84.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dRap3_xpiy = (19305.0/8.0)*xpiy2*z*(-7.0*xmiy3*xpiy3 + 70.0*xmiy2*xpiy2*z2 - 112.0*xmiy*xpiy*z4 + 32.0*z6);
   complex_double dRap3_xmiy = (135135.0/16.0)*xpiy4*z*(-xmiy2*xpiy2 + 8.0*xmiy*xpiy*z2 - 8.0*z4);
   RlmX[4*113+1] = dRap3_xmiy + dRap3_xpiy;
   RlmX[4*113+2] = -I*dRap3_xmiy + I*dRap3_xpiy;
   RlmX[4*113+3] = (6435.0/16.0)*xpiy3*z*(168.0*xmiy2*xpiy2*z - 672.0*xmiy*xpiy*z3 + 384.0*z5) + (6435.0/16.0)*xpiy3*(-7.0*xmiy3*xpiy3 + 84.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 64.0*z6);
   RlmX[4*114+0] = (45045.0/16.0)*xpiy4*(-xmiy3*xpiy3 + 42.0*xmiy2*xpiy2*z2 - 168.0*xmiy*xpiy*z4 + 112.0*z6);
   complex_double dRap4_xpiy = (315315.0/16.0)*xpiy3*(-xmiy3*xpiy3 + 36.0*xmiy2*xpiy2*z2 - 120.0*xmiy*xpiy*z4 + 64.0*z6);
   complex_double dRap4_xmiy = (135135.0/16.0)*xpiy5*(-xmiy2*xpiy2 + 28.0*xmiy*xpiy*z2 - 56.0*z4);
   RlmX[4*114+1] = dRap4_xmiy + dRap4_xpiy;
   RlmX[4*114+2] = -I*dRap4_xmiy + I*dRap4_xpiy;
   RlmX[4*114+3] = (45045.0/16.0)*xpiy4*(84.0*xmiy2*xpiy2*z - 672.0*xmiy*xpiy*z3 + 672.0*z5);
   RlmX[4*115+0] = (135135.0/8.0)*xpiy5*z*(15.0*xmiy2*xpiy2 - 140.0*xmiy*xpiy*z2 + 168.0*z4);
   complex_double dRap5_xpiy = (14189175.0/8.0)*xpiy4*z*(xmiy2*xpiy2 - 8.0*xmiy*xpiy*z2 + 8.0*z4);
   complex_double dRap5_xmiy = (675675.0/4.0)*xpiy6*z*(3.0*xmiy*xpiy - 14.0*z2);
   RlmX[4*115+1] = dRap5_xmiy + dRap5_xpiy;
   RlmX[4*115+2] = -I*dRap5_xmiy + I*dRap5_xpiy;
   RlmX[4*115+3] = (135135.0/8.0)*xpiy5*z*(-280.0*xmiy*xpiy*z + 672.0*z3) + (135135.0/8.0)*xpiy5*(15.0*xmiy2*xpiy2 - 140.0*xmiy*xpiy*z2 + 168.0*z4);
   RlmX[4*116+0] = xpiy6*((2027025.0/8.0)*xmiy2*xpiy2 - 8108100.0*xmiy*xpiy*z2 + 18918900.0*z4);
   complex_double dRap6_xpiy = 2027025.0*xpiy5*(xmiy2*xpiy2 - 28.0*xmiy*xpiy*z2 + 56.0*z4);
   complex_double dRap6_xmiy = (2027025.0/4.0)*xpiy7*(xmiy*xpiy - 16.0*z2);
   RlmX[4*116+1] = dRap6_xmiy + dRap6_xpiy;
   RlmX[4*116+2] = -I*dRap6_xmiy + I*dRap6_xpiy;
   RlmX[4*116+3] = xpiy6*(-16216200.0*xmiy*xpiy*z + 75675600.0*z3);
   RlmX[4*117+0] = (11486475.0/2.0)*xpiy7*z*(-3.0*xmiy*xpiy + 16.0*z2);
   complex_double dRap7_xpiy = 45945900.0*xpiy6*z*(-3.0*xmiy*xpiy + 14.0*z2);
   complex_double dRap7_xmiy = -34459425.0/2.0*xpiy8*z;
   RlmX[4*117+1] = dRap7_xmiy + dRap7_xpiy;
   RlmX[4*117+2] = -I*dRap7_xmiy + I*dRap7_xpiy;
   RlmX[4*117+3] = 183783600.0*xpiy7*z2 + (11486475.0/2.0)*xpiy7*(-3.0*xmiy*xpiy + 16.0*z2);
   RlmX[4*118+0] = (34459425.0/2.0)*xpiy8*(-xmiy*xpiy + 18.0*z2);
   complex_double dRap8_xpiy = (310134825.0/2.0)*xpiy7*(-xmiy*xpiy + 16.0*z2);
   complex_double dRap8_xmiy = -34459425.0/2.0*xpiy9;
   RlmX[4*118+1] = dRap8_xmiy + dRap8_xpiy;
   RlmX[4*118+2] = -I*dRap8_xmiy + I*dRap8_xpiy;
   RlmX[4*118+3] = 620269650.0*xpiy8*z;
   RlmX[4*119+0] = 654729075.0*xpiy9*z;
   complex_double dRap9_xpiy = 5892561675.0*xpiy8*z;
   RlmX[4*119+1] = dRap9_xpiy;
   RlmX[4*119+2] = I*dRap9_xpiy;
   RlmX[4*119+3] = 654729075.0*xpiy9;
   RlmX[4*120+0] = 654729075.0*xpiy10;
   complex_double dRapa_xpiy = 6547290750.0*xpiy9;
   RlmX[4*120+1] = dRapa_xpiy;
   RlmX[4*120+2] = I*dRapa_xpiy;
   RlmX[4*120+3] = 0;
   if (10 == MaxL) return;
   // suppress potential unused warning(s)
   (void)z10;
}

// Renormalization factors for complex Rlm:
// Inverse square roots (Nsq^(-1/2)) of solid angle integrals Nsq = (1/((2*l+1) 4pi))*Integrate[conj(Rlm)*Rlm, dOmega]
static double const fNR0z0 = 1.;
static double const fNR1m1 = 1.4142135623730951;
static double const fNR1z0 = 1.;
static double const fNR1p1 = 7.0710678118654757e-01;
static double const fNR2m2 = 4.8989794855663558;
static double const fNR2m1 = 2.4494897427831779;
static double const fNR2z0 = 1.;
static double const fNR2p1 = 4.0824829046386302e-01;
static double const fNR2p2 = 2.0412414523193151e-01;
static double const fNR3m3 = 2.6832815729997478e+01;
static double const fNR3m2 = 1.0954451150103322e+01;
static double const fNR3m1 = 3.4641016151377544;
static double const fNR3z0 = 1.;
static double const fNR3p1 = 2.8867513459481287e-01;
static double const fNR3p2 = 9.1287092917527679e-02;
static double const fNR3p3 = 3.7267799624996496e-02;
static double const fNR4m4 = 2.0079840636817812e+02;
static double const fNR4m3 = 7.0992957397195397e+01;
static double const fNR4m2 = 1.8973665961010276e+01;
static double const fNR4m1 = 4.4721359549995796;
static double const fNR4z0 = 1.;
static double const fNR4p1 = 2.2360679774997896e-01;
static double const fNR4p2 = 5.2704627669472988e-02;
static double const fNR4p3 = 1.4085904245475277e-02;
static double const fNR4p4 = 4.9801192055599734e-03;
static double const fNR5m5 = 1.9049409439665053e+03;
static double const fNR5m4 = 6.0239521910453436e+02;
static double const fNR5m3 = 1.4198591479439079e+02;
static double const fNR5m2 = 2.8982753492378876e+01;
static double const fNR5m1 = 5.4772255750516612;
static double const fNR5z0 = 1.;
static double const fNR5p1 = 1.8257418583505536e-01;
static double const fNR5p2 = 3.4503277967117711e-02;
static double const fNR5p3 = 7.0429521227376385e-03;
static double const fNR5p4 = 1.660039735186658e-03;
static double const fNR5p5 = 5.2495065695726002e-04;
static double const fNR6m6 = 2.1886105181141756e+04;
static double const fNR6m5 = 6.3179743589223281e+03;
static double const fNR6m4 = 1.3469966592386188e+03;
static double const fNR6m3 = 2.4592681838303037e+02;
static double const fNR6m2 = 4.0987803063838392e+01;
static double const fNR6m1 = 6.4807406984078604;
static double const fNR6z0 = 1.;
static double const fNR6p1 = 1.5430334996209191e-01;
static double const fNR6p2 = 2.4397501823713332e-02;
static double const fNR6p3 = 4.0662503039522215e-03;
static double const fNR6p4 = 7.4239233864562329e-04;
static double const fNR6p5 = 1.5827857841616382e-04;
static double const fNR6p6 = 4.5691089927761737e-05;

// Return renormalization factors for complex Rlm:
// Inverse square roots (Nsq^(-1/2)) of solid angle integrals Nsq = (1/((2*l+1) 4pi))*Integrate[conj(Rlm)*Rlm, dOmega]
double GetRlmRenormFactor(int la, int ma)
{
   assert(la >= 0 && ma <= la && -ma <= la);
   assert(la <= 6);
   switch(iSlmX(la,ma)) {
      case 0: return fNR0z0; // la = 0, ma = 0
      case 1: return fNR1m1; // la = 1, ma = -1
      case 2: return fNR1z0; // la = 1, ma = 0
      case 3: return fNR1p1; // la = 1, ma = 1
      case 4: return fNR2m2; // la = 2, ma = -2
      case 5: return fNR2m1; // la = 2, ma = -1
      case 6: return fNR2z0; // la = 2, ma = 0
      case 7: return fNR2p1; // la = 2, ma = 1
      case 8: return fNR2p2; // la = 2, ma = 2
      case 9: return fNR3m3; // la = 3, ma = -3
      case 10: return fNR3m2; // la = 3, ma = -2
      case 11: return fNR3m1; // la = 3, ma = -1
      case 12: return fNR3z0; // la = 3, ma = 0
      case 13: return fNR3p1; // la = 3, ma = 1
      case 14: return fNR3p2; // la = 3, ma = 2
      case 15: return fNR3p3; // la = 3, ma = 3
      case 16: return fNR4m4; // la = 4, ma = -4
      case 17: return fNR4m3; // la = 4, ma = -3
      case 18: return fNR4m2; // la = 4, ma = -2
      case 19: return fNR4m1; // la = 4, ma = -1
      case 20: return fNR4z0; // la = 4, ma = 0
      case 21: return fNR4p1; // la = 4, ma = 1
      case 22: return fNR4p2; // la = 4, ma = 2
      case 23: return fNR4p3; // la = 4, ma = 3
      case 24: return fNR4p4; // la = 4, ma = 4
      case 25: return fNR5m5; // la = 5, ma = -5
      case 26: return fNR5m4; // la = 5, ma = -4
      case 27: return fNR5m3; // la = 5, ma = -3
      case 28: return fNR5m2; // la = 5, ma = -2
      case 29: return fNR5m1; // la = 5, ma = -1
      case 30: return fNR5z0; // la = 5, ma = 0
      case 31: return fNR5p1; // la = 5, ma = 1
      case 32: return fNR5p2; // la = 5, ma = 2
      case 33: return fNR5p3; // la = 5, ma = 3
      case 34: return fNR5p4; // la = 5, ma = 4
      case 35: return fNR5p5; // la = 5, ma = 5
      case 36: return fNR6m6; // la = 6, ma = -6
      case 37: return fNR6m5; // la = 6, ma = -5
      case 38: return fNR6m4; // la = 6, ma = -4
      case 39: return fNR6m3; // la = 6, ma = -3
      case 40: return fNR6m2; // la = 6, ma = -2
      case 41: return fNR6m1; // la = 6, ma = -1
      case 42: return fNR6z0; // la = 6, ma = 0
      case 43: return fNR6p1; // la = 6, ma = 1
      case 44: return fNR6p2; // la = 6, ma = 2
      case 45: return fNR6p3; // la = 6, ma = 3
      case 46: return fNR6p4; // la = 6, ma = 4
      case 47: return fNR6p5; // la = 6, ma = 5
      case 48: return fNR6p6; // la = 6, ma = 6
   }
   assert(0);
   return 0.;
}

size_t _m2c(int la, int ma)
{
   assert(la >= 0 && ma <= la && -ma <= la);
   assert(la <= 6);
   switch(iSlmX(la,ma)) {
      case 0: return 0; // la = 0, ma = 0
      case 1: return 1; // la = 1, ma = -1
      case 2: return 2; // la = 1, ma = 0
      case 3: return 0; // la = 1, ma = 1
      case 4: return 1; // la = 2, ma = -2
      case 5: return 4; // la = 2, ma = -1
      case 6: return 0; // la = 2, ma = 0
      case 7: return 2; // la = 2, ma = 1
      case 8: return 3; // la = 2, ma = 2
      case 9: return 5; // la = 3, ma = -3
      case 10: return 4; // la = 3, ma = -2
      case 11: return 1; // la = 3, ma = -1
      case 12: return 2; // la = 3, ma = 0
      case 13: return 0; // la = 3, ma = 1
      case 14: return 6; // la = 3, ma = 2
      case 15: return 3; // la = 3, ma = 3
      case 16: return 6; // la = 4, ma = -4
      case 17: return 8; // la = 4, ma = -3
      case 18: return 1; // la = 4, ma = -2
      case 19: return 4; // la = 4, ma = -1
      case 20: return 0; // la = 4, ma = 0
      case 21: return 2; // la = 4, ma = 1
      case 22: return 5; // la = 4, ma = 2
      case 23: return 7; // la = 4, ma = 3
      case 24: return 3; // la = 4, ma = 4
      case 25: return 7; // la = 5, ma = -5
      case 26: return 4; // la = 5, ma = -4
      case 27: return 5; // la = 5, ma = -3
      case 28: return 10; // la = 5, ma = -2
      case 29: return 1; // la = 5, ma = -1
      case 30: return 8; // la = 5, ma = 0
      case 31: return 0; // la = 5, ma = 1
      case 32: return 2; // la = 5, ma = 2
      case 33: return 3; // la = 5, ma = 3
      case 34: return 6; // la = 5, ma = 4
      case 35: return 9; // la = 5, ma = 5
      case 36: return 6; // la = 6, ma = -6
      case 37: return 4; // la = 6, ma = -5
      case 38: return 8; // la = 6, ma = -4
      case 39: return 10; // la = 6, ma = -3
      case 40: return 1; // la = 6, ma = -2
      case 41: return 11; // la = 6, ma = -1
      case 42: return 9; // la = 6, ma = 0
      case 43: return 12; // la = 6, ma = 1
      case 44: return 5; // la = 6, ma = 2
      case 45: return 7; // la = 6, ma = 3
      case 46: return 3; // la = 6, ma = 4
      case 47: return 2; // la = 6, ma = 5
      case 48: return 0; // la = 6, ma = 6
   }
   assert(0);
   return 0;
}

// Convert one set of (2l+1) data entries [d(l,m) for m in (-l,..,l)]
// denoting Rlm integrals to corresponding Slc integrals.
// Input data is taken from pIn[si * (l+m)], output is written to pOut[so * c].
// Input and output arrays may overlap if and only if strides match (so == si)
void ConvertRlmToSlm1(complex_double *pOut, size_t so, complex_double const *pIn, size_t si, int l)
{
   assert(pOut != pIn || si == so);
   assert(l >= 0 && l <= 6);
   complex_double const I(0,1); // imaginary unit
   switch(l) {
      case 0: {
         complex_double r0 = pIn[0*si];
         pOut[0*so] = r0;
         break;
      }
      case 1: {
         complex_double rn1 = pIn[0*si], r0 = pIn[1*si], rp1 = pIn[2*si];
         pOut[2*so] = r0;
         pOut[0*so] = -rn1 + 0.5*rp1;
         pOut[1*so] = -I*rn1 - 0.5*I*rp1;
         break;
      }
      case 2: {
         complex_double rn2 = pIn[0*si], rn1 = pIn[1*si], r0 = pIn[2*si], rp1 = pIn[3*si], rp2 = pIn[4*si];
         pOut[0*so] = r0;
         pOut[2*so] = -1.7320508075688772*rn1 + 0.28867513459481287*rp1;
         pOut[4*so] = -1.7320508075688772*I*rn1 - 0.28867513459481287*I*rp1;
         pOut[3*so] = 3.4641016151377544*rn2 + 0.14433756729740643*rp2;
         pOut[1*so] = 3.4641016151377544*I*rn2 - 0.14433756729740643*I*rp2;
         break;
      }
      case 3: {
         complex_double rn3 = pIn[0*si], rn2 = pIn[1*si], rn1 = pIn[2*si], r0 = pIn[3*si], rp1 = pIn[4*si], rp2 = pIn[5*si], rp3 = pIn[6*si];
         pOut[2*so] = r0;
         pOut[0*so] = -2.4494897427831779*rn1 + 0.20412414523193148*rp1;
         pOut[1*so] = -2.4494897427831779*I*rn1 - 0.20412414523193148*I*rp1;
         pOut[6*so] = 7.745966692414834*rn2 + 0.064549722436790288*rp2;
         pOut[4*so] = 7.745966692414834*I*rn2 - 0.064549722436790288*I*rp2;
         pOut[3*so] = -18.973665961010276*rn3 + 0.026352313834736494*rp3;
         pOut[5*so] = -18.973665961010276*I*rn3 - 0.026352313834736494*I*rp3;
         break;
      }
      case 4: {
         complex_double rn4 = pIn[0*si], rn3 = pIn[1*si], rn2 = pIn[2*si], rn1 = pIn[3*si], r0 = pIn[4*si], rp1 = pIn[5*si], rp2 = pIn[6*si], rp3 = pIn[7*si], rp4 = pIn[8*si];
         pOut[0*so] = r0;
         pOut[2*so] = -3.1622776601683795*rn1 + 0.158113883008419*rp1;
         pOut[4*so] = -3.1622776601683795*I*rn1 - 0.158113883008419*I*rp1;
         pOut[5*so] = 13.416407864998739*rn2 + 0.037267799624996496*rp2;
         pOut[1*so] = 13.416407864998739*I*rn2 - 0.037267799624996496*I*rp2;
         pOut[7*so] = -50.19960159204453*rn3 + 0.0099602384111199486*rp3;
         pOut[8*so] = -50.19960159204453*I*rn3 - 0.0099602384111199486*I*rp3;
         pOut[3*so] = 141.98591479439079*rn4 + 0.0035214760613688193*rp4;
         pOut[6*so] = 141.98591479439079*I*rn4 - 0.0035214760613688193*I*rp4;
         break;
      }
      case 5: {
         complex_double rn5 = pIn[0*si], rn4 = pIn[1*si], rn3 = pIn[2*si], rn2 = pIn[3*si], rn1 = pIn[4*si], r0 = pIn[5*si], rp1 = pIn[6*si], rp2 = pIn[7*si], rp3 = pIn[8*si], rp4 = pIn[9*si], rp5 = pIn[10*si];
         pOut[8*so] = r0;
         pOut[0*so] = -3.872983346207417*rn1 + 0.12909944487358058*rp1;
         pOut[1*so] = -3.872983346207417*I*rn1 - 0.12909944487358058*I*rp1;
         pOut[2*so] = 20.493901531919196*rn2 + 0.024397501823713332*rp2;
         pOut[10*so] = 20.493901531919196*I*rn2 - 0.024397501823713332*I*rp2;
         pOut[3*so] = -100.39920318408906*rn3 + 0.0049801192055599743*rp3;
         pOut[5*so] = -100.39920318408906*I*rn3 - 0.0049801192055599743*I*rp3;
         pOut[6*so] = 425.95774438317238*rn4 + 0.0011738253537896064*rp4;
         pOut[4*so] = 425.95774438317238*I*rn4 - 0.0011738253537896064*I*rp4;
         pOut[9*so] = -1346.9966592386188*rn5 + 0.00037119616932281165*rp5;
         pOut[7*so] = -1346.9966592386188*I*rn5 - 0.00037119616932281165*I*rp5;
         break;
      }
      case 6: {
         complex_double rn6 = pIn[0*si], rn5 = pIn[1*si], rn4 = pIn[2*si], rn3 = pIn[3*si], rn2 = pIn[4*si], rn1 = pIn[5*si], r0 = pIn[6*si], rp1 = pIn[7*si], rp2 = pIn[8*si], rp3 = pIn[9*si], rp4 = pIn[10*si], rp5 = pIn[11*si], rp6 = pIn[12*si];
         pOut[9*so] = r0;
         pOut[12*so] = -4.5825756949558398*rn1 + 0.10910894511799618*rp1;
         pOut[11*so] = -4.5825756949558398*I*rn1 - 0.10910894511799618*I*rp1;
         pOut[5*so] = 28.982753492378876*rn2 + 0.017251638983558856*rp2;
         pOut[1*so] = 28.982753492378876*I*rn2 - 0.017251638983558856*I*rp2;
         pOut[7*so] = -173.89652095427326*rn3 + 0.0028752731639264759*rp3;
         pOut[10*so] = -173.89652095427326*I*rn3 - 0.0028752731639264759*I*rp3;
         pOut[3*so] = 952.47047198325265*rn4 + 0.00052495065695726002*rp4;
         pOut[8*so] = 952.47047198325265*I*rn4 - 0.00052495065695726002*I*rp4;
         pOut[2*so] = -4467.4825125567086*rn5 + 0.00011191985611463616*rp5;
         pOut[4*so] = -4467.4825125567086*I*rn5 - 0.00011191985611463616*I*rp5;
         pOut[0*so] = 15475.813387347367*rn6 + 3.2308479527724681e-5*rp6;
         pOut[6*so] = 15475.813387347367*I*rn6 - 3.2308479527724681e-5*I*rp6;
         break;
      }
   }
}

static double const s_3RlmSphereInts[3617] = {
   1.0, -1.0/3.0, 1.0/3.0, -1.0/3.0, 1.0/5.0, -1.0/5.0, 1.0/5.0, -1.0/5.0,
   1.0/5.0, -1.0/7.0, 1.0/7.0, -1.0/7.0, 1.0/7.0, -1.0/7.0, 1.0/7.0, -1.0/7.0,
   1.0/9.0, -1.0/9.0, 1.0/9.0, -1.0/9.0, 1.0/9.0, -1.0/9.0, 1.0/9.0, -1.0/9.0,
   1.0/9.0, -1.0/11.0, 1.0/11.0, -1.0/11.0, 1.0/11.0, -1.0/11.0, 1.0/11.0, -1.0/11.0,
   1.0/11.0, -1.0/11.0, 1.0/11.0, -1.0/11.0, 1.0/13.0, -1.0/13.0, 1.0/13.0, -1.0/13.0,
   1.0/13.0, -1.0/13.0, 1.0/13.0, -1.0/13.0, 1.0/13.0, -1.0/13.0, 1.0/13.0, -1.0/13.0,
   1.0/13.0, 0, 0, -1.0/3.0, 0, 0, 1.0/15.0, -1.0/5.0,
   2.0/5.0, 0, 0, -1.0/35.0, 3.0/35.0, -6.0/35.0, 2.0/7.0, -3.0/7.0,
   0, 0, 1.0/63.0, -1.0/21.0, 2.0/21.0, -10.0/63.0, 5.0/21.0, -1.0/3.0,
   4.0/9.0, 0, 0, -1.0/99.0, 1.0/33.0, -2.0/33.0, 10.0/99.0, -5.0/33.0,
   7.0/33.0, -28.0/99.0, 4.0/11.0, -5.0/11.0, 0, 0, 1.0/143.0, -3.0/143.0,
   6.0/143.0, -10.0/143.0, 15.0/143.0, -21.0/143.0, 28.0/143.0, -36.0/143.0, 45.0/143.0, -5.0/13.0,
   6.0/13.0, 2.0/5.0, -1.0/5.0, 1.0/15.0, -3.0/7.0, 2.0/7.0, -6.0/35.0, 3.0/35.0,
   -1.0/35.0, 4.0/9.0, -1.0/3.0, 5.0/21.0, -10.0/63.0, 2.0/21.0, -1.0/21.0, 1.0/63.0,
   -5.0/11.0, 4.0/11.0, -28.0/99.0, 7.0/33.0, -5.0/33.0, 10.0/99.0, -2.0/33.0, 1.0/33.0,
   -1.0/99.0, 6.0/13.0, -5.0/13.0, 45.0/143.0, -36.0/143.0, 28.0/143.0, -21.0/143.0, 15.0/143.0,
   -10.0/143.0, 6.0/143.0, -3.0/143.0, 1.0/143.0, -7.0/15.0, 2.0/5.0, -22.0/65.0, 11.0/39.0,
   -3.0/13.0, 12.0/65.0, -28.0/195.0, 7.0/65.0, -1.0/13.0, 2.0/39.0, -2.0/65.0, 1.0/65.0,
   -1.0/195.0, 0, 1.0/3.0, 0, 0, -1.0/15.0, 2.0/15.0, -1.0/5.0,
   0, 0, 1.0/35.0, -2.0/35.0, 3.0/35.0, -4.0/35.0, 1.0/7.0, 0,
   0, -1.0/63.0, 2.0/63.0, -1.0/21.0, 4.0/63.0, -5.0/63.0, 2.0/21.0, -1.0/9.0,
   0, 0, 1.0/99.0, -2.0/99.0, 1.0/33.0, -4.0/99.0, 5.0/99.0, -2.0/33.0,
   7.0/99.0, -8.0/99.0, 1.0/11.0, 0, 0, -1.0/143.0, 2.0/143.0, -3.0/143.0,
   4.0/143.0, -5.0/143.0, 6.0/143.0, -7.0/143.0, 8.0/143.0, -9.0/143.0, 10.0/143.0, -1.0/13.0,
   0, -1.0/5.0, 2.0/15.0, -1.0/15.0, 1.0/7.0, -4.0/35.0, 3.0/35.0, -2.0/35.0,
   1.0/35.0, -1.0/9.0, 2.0/21.0, -5.0/63.0, 4.0/63.0, -1.0/21.0, 2.0/63.0, -1.0/63.0,
   1.0/11.0, -8.0/99.0, 7.0/99.0, -2.0/33.0, 5.0/99.0, -4.0/99.0, 1.0/33.0, -2.0/99.0,
   1.0/99.0, -1.0/13.0, 10.0/143.0, -9.0/143.0, 8.0/143.0, -7.0/143.0, 6.0/143.0, -5.0/143.0,
   4.0/143.0, -3.0/143.0, 2.0/143.0, -1.0/143.0, 1.0/15.0, -4.0/65.0, 11.0/195.0, -2.0/39.0,
   3.0/65.0, -8.0/195.0, 7.0/195.0, -2.0/65.0, 1.0/39.0, -4.0/195.0, 1.0/65.0, -2.0/195.0,
   1.0/195.0, -1.0/3.0, 0, 0, 1.0/15.0, -1.0/15.0, 1.0/15.0, 0,
   0, -1.0/35.0, 1.0/35.0, -1.0/35.0, 1.0/35.0, -1.0/35.0, 0, 0,
   1.0/63.0, -1.0/63.0, 1.0/63.0, -1.0/63.0, 1.0/63.0, -1.0/63.0, 1.0/63.0, 0,
   0, -1.0/99.0, 1.0/99.0, -1.0/99.0, 1.0/99.0, -1.0/99.0, 1.0/99.0, -1.0/99.0,
   1.0/99.0, -1.0/99.0, 0, 0, 1.0/143.0, -1.0/143.0, 1.0/143.0, -1.0/143.0,
   1.0/143.0, -1.0/143.0, 1.0/143.0, -1.0/143.0, 1.0/143.0, -1.0/143.0, 1.0/143.0, 0,
   0, 1.0/15.0, -1.0/15.0, 1.0/15.0, -1.0/35.0, 1.0/35.0, -1.0/35.0, 1.0/35.0,
   -1.0/35.0, 1.0/63.0, -1.0/63.0, 1.0/63.0, -1.0/63.0, 1.0/63.0, -1.0/63.0, 1.0/63.0,
   -1.0/99.0, 1.0/99.0, -1.0/99.0, 1.0/99.0, -1.0/99.0, 1.0/99.0, -1.0/99.0, 1.0/99.0,
   -1.0/99.0, 1.0/143.0, -1.0/143.0, 1.0/143.0, -1.0/143.0, 1.0/143.0, -1.0/143.0, 1.0/143.0,
   -1.0/143.0, 1.0/143.0, -1.0/143.0, 1.0/143.0, -1.0/195.0, 1.0/195.0, -1.0/195.0, 1.0/195.0,
   -1.0/195.0, 1.0/195.0, -1.0/195.0, 1.0/195.0, -1.0/195.0, 1.0/195.0, -1.0/195.0, 1.0/195.0,
   -1.0/195.0, 0, 0, 0, 0, 1.0/5.0, 0, 0,
   0, 0, -1.0/35.0, 1.0/7.0, -3.0/7.0, 0, 0, 0,
   0, 1.0/105.0, -1.0/21.0, 1.0/7.0, -1.0/3.0, 2.0/3.0, 0, 0,
   0, 0, -1.0/231.0, 5.0/231.0, -5.0/77.0, 5.0/33.0, -10.0/33.0, 6.0/11.0,
   -10.0/11.0, 0, 0, 0, 0, 1.0/429.0, -5.0/429.0, 5.0/143.0,
   -35.0/429.0, 70.0/429.0, -42.0/143.0, 70.0/143.0, -10.0/13.0, 15.0/13.0, 0, 0,
   -2.0/35.0, 3.0/35.0, -2.0/35.0, 0, 0, 1.0/21.0, -2.0/21.0, 4.0/35.0,
   -2.0/21.0, 1.0/21.0, 0, 0, -4.0/99.0, 1.0/11.0, -10.0/77.0, 100.0/693.0,
   -10.0/77.0, 1.0/11.0, -4.0/99.0, 0, 0, 5.0/143.0, -12.0/143.0, 56.0/429.0,
   -70.0/429.0, 25.0/143.0, -70.0/429.0, 56.0/429.0, -12.0/143.0, 5.0/143.0, 0, 0,
   -2.0/65.0, 1.0/13.0, -18.0/143.0, 24.0/143.0, -28.0/143.0, 147.0/715.0, -28.0/143.0, 24.0/143.0,
   -18.0/143.0, 1.0/13.0, -2.0/65.0, 2.0/3.0, -1.0/3.0, 1.0/7.0, -1.0/21.0, 1.0/105.0,
   -10.0/11.0, 6.0/11.0, -10.0/33.0, 5.0/33.0, -5.0/77.0, 5.0/231.0, -1.0/231.0, 15.0/13.0,
   -10.0/13.0, 70.0/143.0, -42.0/143.0, 70.0/429.0, -35.0/429.0, 5.0/143.0, -5.0/429.0, 1.0/429.0,
   -7.0/5.0, 1.0, -9.0/13.0, 6.0/13.0, -42.0/143.0, 126.0/715.0, -14.0/143.0, 7.0/143.0,
   -3.0/143.0, 1.0/143.0, -1.0/715.0, 28.0/17.0, -21.0/17.0, 77.0/85.0, -11.0/17.0, 99.0/221.0,
   -66.0/221.0, 42.0/221.0, -126.0/1105.0, 14.0/221.0, -7.0/221.0, 3.0/221.0, -1.0/221.0, 1.0/1105.0,
   0, 0, 0, -1.0/5.0, 0, 0, 0, 0,
   1.0/35.0, -4.0/35.0, 2.0/7.0, 0, 0, 0, 0, -1.0/105.0,
   4.0/105.0, -2.0/21.0, 4.0/21.0, -1.0/3.0, 0, 0, 0, 0,
   1.0/231.0, -4.0/231.0, 10.0/231.0, -20.0/231.0, 5.0/33.0, -8.0/33.0, 4.0/11.0, 0,
   0, 0, 0, -1.0/429.0, 4.0/429.0, -10.0/429.0, 20.0/429.0, -35.0/429.0,
   56.0/429.0, -28.0/143.0, 40.0/143.0, -5.0/13.0, 0, 0, 2.0/35.0, -1.0/35.0,
   -1.0/35.0, 2.0/35.0, 0, -1.0/21.0, 1.0/21.0, -2.0/105.0, -2.0/105.0, 1.0/21.0,
   -1.0/21.0, 0, 4.0/99.0, -5.0/99.0, 3.0/77.0, -10.0/693.0, -10.0/693.0, 3.0/77.0,
   -5.0/99.0, 4.0/99.0, 0, -5.0/143.0, 7.0/143.0, -20.0/429.0, 14.0/429.0, -5.0/429.0,
   -5.0/429.0, 14.0/429.0, -20.0/429.0, 7.0/143.0, -5.0/143.0, 0, 2.0/65.0, -3.0/65.0,
   7.0/143.0, -6.0/143.0, 4.0/143.0, -7.0/715.0, -7.0/715.0, 4.0/143.0, -6.0/143.0, 7.0/143.0,
   -3.0/65.0, 2.0/65.0, -1.0/3.0, 4.0/21.0, -2.0/21.0, 4.0/105.0, -1.0/105.0, 4.0/11.0,
   -8.0/33.0, 5.0/33.0, -20.0/231.0, 10.0/231.0, -4.0/231.0, 1.0/231.0, -5.0/13.0, 40.0/143.0,
   -28.0/143.0, 56.0/429.0, -35.0/429.0, 20.0/429.0, -10.0/429.0, 4.0/429.0, -1.0/429.0, 2.0/5.0,
   -4.0/13.0, 3.0/13.0, -24.0/143.0, 84.0/715.0, -56.0/715.0, 7.0/143.0, -4.0/143.0, 2.0/143.0,
   -4.0/715.0, 1.0/715.0, -7.0/17.0, 28.0/85.0, -22.0/85.0, 44.0/221.0, -33.0/221.0, 24.0/221.0,
   -84.0/1105.0, 56.0/1105.0, -7.0/221.0, 4.0/221.0, -2.0/221.0, 4.0/1105.0, -1.0/1105.0, 0,
   0, 1.0/5.0, 0, 0, 0, 0, -1.0/35.0, 3.0/35.0,
   -6.0/35.0, 0, 0, 0, 0, 1.0/105.0, -1.0/35.0, 2.0/35.0,
   -2.0/21.0, 1.0/7.0, 0, 0, 0, 0, -1.0/231.0, 1.0/77.0,
   -2.0/77.0, 10.0/231.0, -5.0/77.0, 1.0/11.0, -4.0/33.0, 0, 0, 0,
   0, 1.0/429.0, -1.0/143.0, 2.0/143.0, -10.0/429.0, 5.0/143.0, -7.0/143.0, 28.0/429.0,
   -12.0/143.0, 15.0/143.0, 0, 0, -2.0/35.0, -1.0/35.0, 2.0/35.0, -1.0/35.0,
   -2.0/35.0, 1.0/21.0, 0, -1.0/35.0, 4.0/105.0, -1.0/35.0, 0, 1.0/21.0,
   -4.0/99.0, 1.0/99.0, 8.0/693.0, -17.0/693.0, 20.0/693.0, -17.0/693.0, 8.0/693.0, 1.0/99.0,
   -4.0/99.0, 5.0/143.0, -2.0/143.0, -1.0/429.0, 2.0/143.0, -3.0/143.0, 10.0/429.0, -3.0/143.0,
   2.0/143.0, -1.0/429.0, -2.0/143.0, 5.0/143.0, -2.0/65.0, 1.0/65.0, -2.0/715.0, -1.0/143.0,
   2.0/143.0, -1.0/55.0, 14.0/715.0, -1.0/55.0, 2.0/143.0, -1.0/143.0, -2.0/715.0, 1.0/65.0,
   -2.0/65.0, 1.0/7.0, -2.0/21.0, 2.0/35.0, -1.0/35.0, 1.0/105.0, -4.0/33.0, 1.0/11.0,
   -5.0/77.0, 10.0/231.0, -2.0/77.0, 1.0/77.0, -1.0/231.0, 15.0/143.0, -12.0/143.0, 28.0/429.0,
   -7.0/143.0, 5.0/143.0, -10.0/429.0, 2.0/143.0, -1.0/143.0, 1.0/429.0, -6.0/65.0, 1.0/13.0,
   -9.0/143.0, 36.0/715.0, -28.0/715.0, 21.0/715.0, -3.0/143.0, 2.0/143.0, -6.0/715.0, 3.0/715.0,
   -1.0/715.0, 7.0/85.0, -6.0/85.0, 66.0/1105.0, -11.0/221.0, 9.0/221.0, -36.0/1105.0, 28.0/1105.0,
   -21.0/1105.0, 3.0/221.0, -2.0/221.0, 6.0/1105.0, -3.0/1105.0, 1.0/1105.0, 0, -1.0/5.0,
   0, 0, 0, 0, 1.0/35.0, -2.0/35.0, 3.0/35.0, 0,
   0, 0, 0, -1.0/105.0, 2.0/105.0, -1.0/35.0, 4.0/105.0, -1.0/21.0,
   0, 0, 0, 0, 1.0/231.0, -2.0/231.0, 1.0/77.0, -4.0/231.0,
   5.0/231.0, -2.0/77.0, 1.0/33.0, 0, 0, 0, 0, -1.0/429.0,
   2.0/429.0, -1.0/143.0, 4.0/429.0, -5.0/429.0, 2.0/143.0, -7.0/429.0, 8.0/429.0, -3.0/143.0,
   0, 0, 0, 3.0/35.0, -1.0/35.0, -1.0/35.0, 3.0/35.0, 0,
   -1.0/21.0, 1.0/35.0, -1.0/105.0, -1.0/105.0, 1.0/35.0, -1.0/21.0, 0, 1.0/33.0,
   -5.0/231.0, 1.0/77.0, -1.0/231.0, -1.0/231.0, 1.0/77.0, -5.0/231.0, 1.0/33.0, 0,
   -3.0/143.0, 7.0/429.0, -5.0/429.0, 1.0/143.0, -1.0/429.0, -1.0/429.0, 1.0/143.0, -5.0/429.0,
   7.0/429.0, -3.0/143.0, 0, 1.0/65.0, -9.0/715.0, 7.0/715.0, -1.0/143.0, 3.0/715.0,
   -1.0/715.0, -1.0/715.0, 3.0/715.0, -1.0/143.0, 7.0/715.0, -9.0/715.0, 1.0/65.0, 0,
   -1.0/21.0, 4.0/105.0, -1.0/35.0, 2.0/105.0, -1.0/105.0, 1.0/33.0, -2.0/77.0, 5.0/231.0,
   -4.0/231.0, 1.0/77.0, -2.0/231.0, 1.0/231.0, -3.0/143.0, 8.0/429.0, -7.0/429.0, 2.0/143.0,
   -5.0/429.0, 4.0/429.0, -1.0/143.0, 2.0/429.0, -1.0/429.0, 1.0/65.0, -2.0/143.0, 9.0/715.0,
   -8.0/715.0, 7.0/715.0, -6.0/715.0, 1.0/143.0, -4.0/715.0, 3.0/715.0, -2.0/715.0, 1.0/715.0,
   -1.0/85.0, 12.0/1105.0, -11.0/1105.0, 2.0/221.0, -9.0/1105.0, 8.0/1105.0, -7.0/1105.0, 6.0/1105.0,
   -1.0/221.0, 4.0/1105.0, -3.0/1105.0, 2.0/1105.0, -1.0/1105.0, 1.0/5.0, 0, 0,
   0, 0, -1.0/35.0, 1.0/35.0, -1.0/35.0, 0, 0, 0,
   0, 1.0/105.0, -1.0/105.0, 1.0/105.0, -1.0/105.0, 1.0/105.0, 0, 0,
   0, 0, -1.0/231.0, 1.0/231.0, -1.0/231.0, 1.0/231.0, -1.0/231.0, 1.0/231.0,
   -1.0/231.0, 0, 0, 0, 0, 1.0/429.0, -1.0/429.0, 1.0/429.0,
   -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, 0, 0,
   0, 0, -2.0/35.0, 2.0/35.0, -2.0/35.0, 0, 0, 2.0/105.0,
   -2.0/105.0, 2.0/105.0, -2.0/105.0, 2.0/105.0, 0, 0, -2.0/231.0, 2.0/231.0,
   -2.0/231.0, 2.0/231.0, -2.0/231.0, 2.0/231.0, -2.0/231.0, 0, 0, 2.0/429.0,
   -2.0/429.0, 2.0/429.0, -2.0/429.0, 2.0/429.0, -2.0/429.0, 2.0/429.0, -2.0/429.0, 2.0/429.0,
   0, 0, -2.0/715.0, 2.0/715.0, -2.0/715.0, 2.0/715.0, -2.0/715.0, 2.0/715.0,
   -2.0/715.0, 2.0/715.0, -2.0/715.0, 2.0/715.0, -2.0/715.0, 0, 0, 1.0/105.0,
   -1.0/105.0, 1.0/105.0, -1.0/105.0, 1.0/105.0, -1.0/231.0, 1.0/231.0, -1.0/231.0, 1.0/231.0,
   -1.0/231.0, 1.0/231.0, -1.0/231.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0,
   -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/715.0, 1.0/715.0, -1.0/715.0, 1.0/715.0,
   -1.0/715.0, 1.0/715.0, -1.0/715.0, 1.0/715.0, -1.0/715.0, 1.0/715.0, -1.0/715.0, 1.0/1105.0,
   -1.0/1105.0, 1.0/1105.0, -1.0/1105.0, 1.0/1105.0, -1.0/1105.0, 1.0/1105.0, -1.0/1105.0, 1.0/1105.0,
   -1.0/1105.0, 1.0/1105.0, -1.0/1105.0, 1.0/1105.0, 0, 0, 0, 0,
   0, 0, -1.0/7.0, 0, 0, 0, 0, 0,
   0, 1.0/63.0, -1.0/9.0, 4.0/9.0, 0, 0, 0, 0,
   0, 0, -1.0/231.0, 1.0/33.0, -4.0/33.0, 4.0/11.0, -10.0/11.0, 0,
   0, 0, 0, 0, 0, 5.0/3003.0, -5.0/429.0, 20.0/429.0,
   -20.0/143.0, 50.0/143.0, -10.0/13.0, 20.0/13.0, 0, 0, 0, 0,
   2.0/105.0, -1.0/21.0, 1.0/21.0, 0, 0, 0, 0, -1.0/77.0,
   10.0/231.0, -6.0/77.0, 1.0/11.0, -2.0/33.0, 0, 0, 0, 0,
   4.0/429.0, -5.0/143.0, 75.0/1001.0, -50.0/429.0, 20.0/143.0, -18.0/143.0, 10.0/143.0, 0,
   0, 0, 0, -1.0/143.0, 4.0/143.0, -28.0/429.0, 49.0/429.0, -70.0/429.0,
   28.0/143.0, -28.0/143.0, 2.0/13.0, -1.0/13.0, 0, 0, -2.0/33.0, 1.0/11.0,
   -6.0/77.0, 10.0/231.0, -1.0/77.0, 0, 0, 10.0/143.0, -18.0/143.0, 20.0/143.0,
   -50.0/429.0, 75.0/1001.0, -5.0/143.0, 4.0/429.0, 0, 0, -1.0/13.0, 2.0/13.0,
   -28.0/143.0, 28.0/143.0, -70.0/429.0, 49.0/429.0, -28.0/429.0, 4.0/143.0, -1.0/143.0, 0,
   0, 7.0/85.0, -3.0/17.0, 54.0/221.0, -60.0/221.0, 630.0/2431.0, -2646.0/12155.0, 392.0/2431.0,
   -252.0/2431.0, 135.0/2431.0, -5.0/221.0, 6.0/1105.0, 20.0/13.0, -10.0/13.0, 50.0/143.0, -20.0/143.0,
   20.0/429.0, -5.0/429.0, 5.0/3003.0, -7.0/3.0, 4.0/3.0, -28.0/39.0, 14.0/39.0, -70.0/429.0,
   28.0/429.0, -28.0/1287.0, 7.0/1287.0, -1.0/1287.0, 56.0/17.0, -35.0/17.0, 21.0/17.0, -12.0/17.0,
   84.0/221.0, -42.0/221.0, 210.0/2431.0, -84.0/2431.0, 28.0/2431.0, -7.0/2431.0, 1.0/2431.0, -84.0/19.0,
   56.0/19.0, -616.0/323.0, 385.0/323.0, -231.0/323.0, 132.0/323.0, -924.0/4199.0, 462.0/4199.0, -210.0/4199.0,
   84.0/4199.0, -28.0/4199.0, 7.0/4199.0, -1.0/4199.0, 0, 0, 0, 0,
   0, 1.0/7.0, 0, 0, 0, 0, 0, 0,
   -1.0/63.0, 2.0/21.0, -1.0/3.0, 0, 0, 0, 0, 0,
   0, 1.0/231.0, -2.0/77.0, 1.0/11.0, -8.0/33.0, 6.0/11.0, 0, 0,
   0, 0, 0, 0, -5.0/3003.0, 10.0/1001.0, -5.0/143.0, 40.0/429.0,
   -30.0/143.0, 60.0/143.0, -10.0/13.0, 0, 0, 0, 0, -2.0/105.0,
   1.0/35.0, 0, -1.0/21.0, 0, 0, 0, 1.0/77.0, -1.0/33.0,
   8.0/231.0, -1.0/77.0, -1.0/33.0, 2.0/33.0, 0, 0, 0, -4.0/429.0,
   1.0/39.0, -40.0/1001.0, 125.0/3003.0, -10.0/429.0, -2.0/143.0, 8.0/143.0, -10.0/143.0, 0,
   0, 0, 1.0/143.0, -3.0/143.0, 16.0/429.0, -7.0/143.0, 7.0/143.0, -14.0/429.0,
   0, 6.0/143.0, -1.0/13.0, 1.0/13.0, 0, 2.0/33.0, -1.0/33.0, -1.0/77.0,
   8.0/231.0, -1.0/33.0, 1.0/77.0, 0, -10.0/143.0, 8.0/143.0, -2.0/143.0, -10.0/429.0,
   125.0/3003.0, -40.0/1001.0, 1.0/39.0, -4.0/429.0, 0, 1.0/13.0, -1.0/13.0, 6.0/143.0,
   0, -14.0/429.0, 7.0/143.0, -7.0/143.0, 16.0/429.0, -3.0/143.0, 1.0/143.0, 0,
   -7.0/85.0, 8.0/85.0, -15.0/221.0, 6.0/221.0, 30.0/2431.0, -504.0/12155.0, 686.0/12155.0, -140.0/2431.0,
   9.0/187.0, -80.0/2431.0, 19.0/1105.0, -6.0/1105.0, -10.0/13.0, 60.0/143.0, -30.0/143.0, 40.0/429.0,
   -5.0/143.0, 10.0/1001.0, -5.0/3003.0, 1.0, -8.0/13.0, 14.0/39.0, -28.0/143.0, 14.0/143.0,
   -56.0/1287.0, 7.0/429.0, -2.0/429.0, 1.0/1287.0, -21.0/17.0, 14.0/17.0, -9.0/17.0, 72.0/221.0,
   -42.0/221.0, 252.0/2431.0, -126.0/2431.0, 56.0/2431.0, -21.0/2431.0, 6.0/2431.0, -1.0/2431.0, 28.0/19.0,
   -336.0/323.0, 231.0/323.0, -154.0/323.0, 99.0/323.0, -792.0/4199.0, 462.0/4199.0, -252.0/4199.0, 126.0/4199.0,
   -56.0/4199.0, 21.0/4199.0, -6.0/4199.0, 1.0/4199.0, 0, 0, 0, 0,
   -1.0/7.0, 0, 0, 0, 0, 0, 0, 1.0/63.0,
   -5.0/63.0, 5.0/21.0, 0, 0, 0, 0, 0, 0,
   -1.0/231.0, 5.0/231.0, -5.0/77.0, 5.0/33.0, -10.0/33.0, 0, 0, 0,
   0, 0, 0, 5.0/3003.0, -25.0/3003.0, 25.0/1001.0, -25.0/429.0, 50.0/429.0,
   -30.0/143.0, 50.0/143.0, 0, 0, 0, 0, 2.0/105.0, -1.0/105.0,
   -1.0/35.0, 1.0/21.0, 1.0/21.0, 0, 0, -1.0/77.0, 4.0/231.0, -1.0/231.0,
   -5.0/231.0, 10.0/231.0, -1.0/33.0, -2.0/33.0, 0, 0, 4.0/429.0, -7.0/429.0,
   43.0/3003.0, -5.0/3003.0, -5.0/273.0, 16.0/429.0, -6.0/143.0, 2.0/143.0, 10.0/143.0, 0,
   0, -1.0/143.0, 2.0/143.0, -7.0/429.0, 5.0/429.0, 0, -7.0/429.0, 14.0/429.0,
   -6.0/143.0, 5.0/143.0, 0, -1.0/13.0, -2.0/33.0, -1.0/33.0, 10.0/231.0, -5.0/231.0,
   -1.0/231.0, 4.0/231.0, -1.0/77.0, 10.0/143.0, 2.0/143.0, -6.0/143.0, 16.0/429.0, -5.0/273.0,
   -5.0/3003.0, 43.0/3003.0, -7.0/429.0, 4.0/429.0, -1.0/13.0, 0, 5.0/143.0, -6.0/143.0,
   14.0/429.0, -7.0/429.0, 0, 5.0/429.0, -7.0/429.0, 2.0/143.0, -1.0/143.0, 7.0/85.0,
   -1.0/85.0, -29.0/1105.0, 9.0/221.0, -96.0/2431.0, 354.0/12155.0, -14.0/935.0, 14.0/12155.0, 23.0/2431.0,
   -37.0/2431.0, 191.0/12155.0, -1.0/85.0, 6.0/1105.0, 50.0/143.0, -30.0/143.0, 50.0/429.0, -25.0/429.0,
   25.0/1001.0, -25.0/3003.0, 5.0/3003.0, -5.0/13.0, 10.0/39.0, -70.0/429.0, 14.0/143.0, -70.0/1287.0,
   35.0/1287.0, -5.0/429.0, 5.0/1287.0, -1.0/1287.0, 7.0/17.0, -5.0/17.0, 45.0/221.0, -30.0/221.0,
   210.0/2431.0, -126.0/2431.0, 70.0/2431.0, -35.0/2431.0, 15.0/2431.0, -5.0/2431.0, 1.0/2431.0, -140.0/323.0,
   105.0/323.0, -77.0/323.0, 55.0/323.0, -495.0/4199.0, 330.0/4199.0, -210.0/4199.0, 126.0/4199.0, -70.0/4199.0,
   35.0/4199.0, -15.0/4199.0, 5.0/4199.0, -1.0/4199.0, 0, 0, 0, 1.0/7.0,
   0, 0, 0, 0, 0, 0, -1.0/63.0, 4.0/63.0,
   -10.0/63.0, 0, 0, 0, 0, 0, 0, 1.0/231.0,
   -4.0/231.0, 10.0/231.0, -20.0/231.0, 5.0/33.0, 0, 0, 0, 0,
   0, 0, -5.0/3003.0, 20.0/3003.0, -50.0/3003.0, 100.0/3003.0, -25.0/429.0, 40.0/429.0,
   -20.0/143.0, 0, 0, 0, 0, -2.0/105.0, -1.0/105.0, 4.0/105.0,
   -2.0/105.0, -2.0/21.0, 0, 0, 1.0/77.0, -1.0/231.0, -1.0/77.0, 2.0/77.0,
   -5.0/231.0, -1.0/77.0, 1.0/11.0, 0, 0, -4.0/429.0, 1.0/143.0, 2.0/1001.0,
   -38.0/3003.0, 20.0/1001.0, -19.0/1001.0, 2.0/429.0, 4.0/143.0, -12.0/143.0, 0, 0,
   1.0/143.0, -1.0/143.0, 1.0/429.0, 2.0/429.0, -5.0/429.0, 7.0/429.0, -7.0/429.0, 4.0/429.0,
   1.0/143.0, -5.0/143.0, 1.0/13.0, 0, 1.0/11.0, -1.0/77.0, -5.0/231.0, 2.0/77.0,
   -1.0/77.0, -1.0/231.0, 1.0/77.0, -12.0/143.0, 4.0/143.0, 2.0/429.0, -19.0/1001.0, 20.0/1001.0,
   -38.0/3003.0, 2.0/1001.0, 1.0/143.0, -4.0/429.0, 1.0/13.0, -5.0/143.0, 1.0/143.0, 4.0/429.0,
   -7.0/429.0, 7.0/429.0, -5.0/429.0, 2.0/429.0, 1.0/429.0, -1.0/143.0, 1.0/143.0, -6.0/85.0,
   42.0/1105.0, -16.0/1105.0, -3.0/2431.0, 126.0/12155.0, -172.0/12155.0, 168.0/12155.0, -129.0/12155.0, 14.0/2431.0,
   -6.0/12155.0, -48.0/12155.0, 7.0/1105.0, -6.0/1105.0, -20.0/143.0, 40.0/429.0, -25.0/429.0, 100.0/3003.0,
   -50.0/3003.0, 20.0/3003.0, -5.0/3003.0, 5.0/39.0, -40.0/429.0, 28.0/429.0, -56.0/1287.0, 35.0/1287.0,
   -20.0/1287.0, 10.0/1287.0, -4.0/1287.0, 1.0/1287.0, -2.0/17.0, 20.0/221.0, -15.0/221.0, 120.0/2431.0,
   -84.0/2431.0, 56.0/2431.0, -35.0/2431.0, 20.0/2431.0, -10.0/2431.0, 4.0/2431.0, -1.0/2431.0, 35.0/323.0,
   -28.0/323.0, 22.0/323.0, -220.0/4199.0, 165.0/4199.0, -120.0/4199.0, 84.0/4199.0, -56.0/4199.0, 35.0/4199.0,
   -20.0/4199.0, 10.0/4199.0, -4.0/4199.0, 1.0/4199.0, 0, 0, -1.0/7.0, 0,
   0, 0, 0, 0, 0, 1.0/63.0, -1.0/21.0, 2.0/21.0,
   0, 0, 0, 0, 0, 0, -1.0/231.0, 1.0/77.0,
   -2.0/77.0, 10.0/231.0, -5.0/77.0, 0, 0, 0, 0, 0,
   0, 5.0/3003.0, -5.0/1001.0, 10.0/1001.0, -50.0/3003.0, 25.0/1001.0, -5.0/143.0, 20.0/429.0,
   0, 0, 0, 0, 2.0/105.0, 1.0/35.0, -1.0/35.0, -2.0/105.0,
   4.0/35.0, 0, 0, -1.0/77.0, -2.0/231.0, 4.0/231.0, -1.0/77.0, -1.0/231.0,
   8.0/231.0, -6.0/77.0, 0, 0, 4.0/429.0, 1.0/429.0, -9.0/1001.0, 32.0/3003.0,
   -2.0/273.0, -1.0/1001.0, 43.0/3003.0, -14.0/429.0, 8.0/143.0, 0, 0, -1.0/143.0,
   0, 2.0/429.0, -1.0/143.0, 1.0/143.0, -2.0/429.0, 0, 1.0/143.0, -7.0/429.0,
   4.0/143.0, -6.0/143.0, 0, 0, -6.0/77.0, 8.0/231.0, -1.0/231.0, -1.0/77.0,
   4.0/231.0, -2.0/231.0, -1.0/77.0, 8.0/143.0, -14.0/429.0, 43.0/3003.0, -1.0/1001.0, -2.0/273.0,
   32.0/3003.0, -9.0/1001.0, 1.0/429.0, 4.0/429.0, -6.0/143.0, 4.0/143.0, -7.0/429.0, 1.0/143.0,
   0, -2.0/429.0, 1.0/143.0, -1.0/143.0, 2.0/429.0, 0, -1.0/143.0, 36.0/1105.0,
   -2.0/85.0, 191.0/12155.0, -111.0/12155.0, 46.0/12155.0, 4.0/12155.0, -3.0/935.0, 59.0/12155.0, -64.0/12155.0,
   54.0/12155.0, -29.0/12155.0, -1.0/1105.0, 6.0/1105.0, 20.0/429.0, -5.0/143.0, 25.0/1001.0, -50.0/3003.0,
   10.0/1001.0, -5.0/1001.0, 5.0/3003.0, -5.0/143.0, 4.0/143.0, -28.0/1287.0, 7.0/429.0, -5.0/429.0,
   10.0/1287.0, -2.0/429.0, 1.0/429.0, -1.0/1287.0, 6.0/221.0, -5.0/221.0, 45.0/2431.0, -36.0/2431.0,
   28.0/2431.0, -21.0/2431.0, 15.0/2431.0, -10.0/2431.0, 6.0/2431.0, -3.0/2431.0, 1.0/2431.0, -7.0/323.0,
   6.0/323.0, -66.0/4199.0, 55.0/4199.0, -45.0/4199.0, 36.0/4199.0, -28.0/4199.0, 21.0/4199.0, -15.0/4199.0,
   10.0/4199.0, -6.0/4199.0, 3.0/4199.0, -1.0/4199.0, 0, 1.0/7.0, 0, 0,
   0, 0, 0, 0, -1.0/63.0, 2.0/63.0, -1.0/21.0, 0,
   0, 0, 0, 0, 0, 1.0/231.0, -2.0/231.0, 1.0/77.0,
   -4.0/231.0, 5.0/231.0, 0, 0, 0, 0, 0, 0,
   -5.0/3003.0, 10.0/3003.0, -5.0/1001.0, 20.0/3003.0, -25.0/3003.0, 10.0/1001.0, -5.0/429.0, 0,
   0, 0, 0, 0, -1.0/21.0, 0, 1.0/21.0, -2.0/21.0,
   0, 0, 0, 5.0/231.0, -2.0/231.0, -1.0/231.0, 4.0/231.0, -1.0/33.0,
   10.0/231.0, 0, 0, 0, -5.0/429.0, 20.0/3003.0, -5.0/3003.0, -10.0/3003.0,
   25.0/3003.0, -40.0/3003.0, 5.0/273.0, -10.0/429.0, 0, 0, 0, 1.0/143.0,
   -2.0/429.0, 1.0/429.0, 0, -1.0/429.0, 2.0/429.0, -1.0/143.0, 4.0/429.0, -5.0/429.0,
   2.0/143.0, 0, 0, 0, 10.0/231.0, -1.0/33.0, 4.0/231.0, -1.0/231.0,
   -2.0/231.0, 5.0/231.0, 0, -10.0/429.0, 5.0/273.0, -40.0/3003.0, 25.0/3003.0, -10.0/3003.0,
   -5.0/3003.0, 20.0/3003.0, -5.0/429.0, 0, 2.0/143.0, -5.0/429.0, 4.0/429.0, -1.0/143.0,
   2.0/429.0, -1.0/429.0, 0, 1.0/429.0, -2.0/429.0, 1.0/143.0, 0, -2.0/221.0,
   19.0/2431.0, -16.0/2431.0, 1.0/187.0, -10.0/2431.0, 7.0/2431.0, -4.0/2431.0, 1.0/2431.0, 2.0/2431.0,
   -5.0/2431.0, 8.0/2431.0, -1.0/221.0, 0, -5.0/429.0, 10.0/1001.0, -25.0/3003.0, 20.0/3003.0,
   -5.0/1001.0, 10.0/3003.0, -5.0/3003.0, 1.0/143.0, -8.0/1287.0, 7.0/1287.0, -2.0/429.0, 5.0/1287.0,
   -4.0/1287.0, 1.0/429.0, -2.0/1287.0, 1.0/1287.0, -1.0/221.0, 10.0/2431.0, -9.0/2431.0, 8.0/2431.0,
   -7.0/2431.0, 6.0/2431.0, -5.0/2431.0, 4.0/2431.0, -3.0/2431.0, 2.0/2431.0, -1.0/2431.0, 1.0/323.0,
   -12.0/4199.0, 11.0/4199.0, -10.0/4199.0, 9.0/4199.0, -8.0/4199.0, 7.0/4199.0, -6.0/4199.0, 5.0/4199.0,
   -4.0/4199.0, 3.0/4199.0, -2.0/4199.0, 1.0/4199.0, -1.0/7.0, 0, 0, 0,
   0, 0, 0, 1.0/63.0, -1.0/63.0, 1.0/63.0, 0, 0,
   0, 0, 0, 0, -1.0/231.0, 1.0/231.0, -1.0/231.0, 1.0/231.0,
   -1.0/231.0, 0, 0, 0, 0, 0, 0, 5.0/3003.0,
   -5.0/3003.0, 5.0/3003.0, -5.0/3003.0, 5.0/3003.0, -5.0/3003.0, 5.0/3003.0, 0, 0,
   0, 0, 0, 0, 1.0/21.0, -1.0/21.0, 1.0/21.0, 0,
   0, 0, 0, -1.0/77.0, 1.0/77.0, -1.0/77.0, 1.0/77.0, -1.0/77.0,
   0, 0, 0, 0, 5.0/1001.0, -5.0/1001.0, 5.0/1001.0, -5.0/1001.0,
   5.0/1001.0, -5.0/1001.0, 5.0/1001.0, 0, 0, 0, 0, -1.0/429.0,
   1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0,
   0, 0, 0, 0, -1.0/77.0, 1.0/77.0, -1.0/77.0, 1.0/77.0,
   -1.0/77.0, 0, 0, 5.0/1001.0, -5.0/1001.0, 5.0/1001.0, -5.0/1001.0, 5.0/1001.0,
   -5.0/1001.0, 5.0/1001.0, 0, 0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0,
   -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 0, 0, 3.0/2431.0,
   -3.0/2431.0, 3.0/2431.0, -3.0/2431.0, 3.0/2431.0, -3.0/2431.0, 3.0/2431.0, -3.0/2431.0, 3.0/2431.0,
   -3.0/2431.0, 3.0/2431.0, 0, 0, 5.0/3003.0, -5.0/3003.0, 5.0/3003.0, -5.0/3003.0,
   5.0/3003.0, -5.0/3003.0, 5.0/3003.0, -1.0/1287.0, 1.0/1287.0, -1.0/1287.0, 1.0/1287.0, -1.0/1287.0,
   1.0/1287.0, -1.0/1287.0, 1.0/1287.0, -1.0/1287.0, 1.0/2431.0, -1.0/2431.0, 1.0/2431.0, -1.0/2431.0,
   1.0/2431.0, -1.0/2431.0, 1.0/2431.0, -1.0/2431.0, 1.0/2431.0, -1.0/2431.0, 1.0/2431.0, -1.0/4199.0,
   1.0/4199.0, -1.0/4199.0, 1.0/4199.0, -1.0/4199.0, 1.0/4199.0, -1.0/4199.0, 1.0/4199.0, -1.0/4199.0,
   1.0/4199.0, -1.0/4199.0, 1.0/4199.0, -1.0/4199.0, 0, 0, 0, 0,
   0, 0, 0, 0, 1.0/9.0, 0, 0, 0,
   0, 0, 0, 0, 0, -1.0/99.0, 1.0/11.0, -5.0/11.0,
   0, 0, 0, 0, 0, 0, 0, 0,
   1.0/429.0, -3.0/143.0, 15.0/143.0, -5.0/13.0, 15.0/13.0, 0, 0, 0,
   0, 0, 0, -2.0/231.0, 1.0/33.0, -4.0/99.0, 0, 0,
   0, 0, 0, 0, 5.0/1001.0, -10.0/429.0, 8.0/143.0, -12.0/143.0,
   10.0/143.0, 0, 0, 0, 0, 0, 0, -4.0/1287.0,
   7.0/429.0, -20.0/429.0, 40.0/429.0, -20.0/143.0, 2.0/13.0, -4.0/39.0, 0, 0,
   0, 0, 2.0/143.0, -5.0/143.0, 45.0/1001.0, -5.0/143.0, 2.0/143.0, 0,
   0, 0, 0, -2.0/143.0, 6.0/143.0, -10.0/143.0, 35.0/429.0, -10.0/143.0,
   6.0/143.0, -2.0/143.0, 0, 0, 0, 0, 3.0/221.0, -10.0/221.0,
   210.0/2431.0, -294.0/2431.0, 980.0/7293.0, -294.0/2431.0, 210.0/2431.0, -10.0/221.0, 3.0/221.0, 0,
   0, -4.0/39.0, 2.0/13.0, -20.0/143.0, 40.0/429.0, -20.0/429.0, 7.0/429.0, -4.0/1287.0,
   0, 0, 7.0/51.0, -4.0/17.0, 56.0/221.0, -140.0/663.0, 350.0/2431.0, -196.0/2431.0,
   784.0/21879.0, -28.0/2431.0, 5.0/2431.0, 0, 0, -56.0/323.0, 105.0/323.0, -126.0/323.0,
   120.0/323.0, -1260.0/4199.0, 882.0/4199.0, -5880.0/46189.0, 3024.0/46189.0, -1260.0/46189.0, 35.0/4199.0, -6.0/4199.0,
   70.0/17.0, -35.0/17.0, 49.0/51.0, -7.0/17.0, 35.0/221.0, -35.0/663.0, 35.0/2431.0, -7.0/2431.0,
   7.0/21879.0, -126.0/19.0, 70.0/19.0, -630.0/323.0, 315.0/323.0, -147.0/323.0, 63.0/323.0, -315.0/4199.0,
   105.0/4199.0, -315.0/46189.0, 63.0/46189.0, -7.0/46189.0, 10.0, -6.0, 66.0/19.0, -110.0/57.0,
   330.0/323.0, -165.0/323.0, 77.0/323.0, -33.0/323.0, 165.0/4199.0, -55.0/4199.0, 15.0/4199.0, -3.0/4199.0,
   1.0/12597.0, 0, 0, 0, 0, 0, 0, 0,
   -1.0/9.0, 0, 0, 0, 0, 0, 0, 0,
   0, 1.0/99.0, -8.0/99.0, 4.0/11.0, 0, 0, 0, 0,
   0, 0, 0, 0, -1.0/429.0, 8.0/429.0, -12.0/143.0, 40.0/143.0,
   -10.0/13.0, 0, 0, 0, 0, 0, 0, 2.0/231.0,
   -5.0/231.0, 1.0/99.0, 4.0/99.0, 0, 0, 0, 0, 0,
   -5.0/1001.0, 5.0/273.0, -14.0/429.0, 4.0/143.0, 2.0/143.0, -10.0/143.0, 0, 0,
   0, 0, 0, 4.0/1287.0, -17.0/1287.0, 1.0/33.0, -20.0/429.0, 20.0/429.0,
   -2.0/143.0, -2.0/39.0, 4.0/39.0, 0, 0, 0, -2.0/143.0, 3.0/143.0,
   -10.0/1001.0, -10.0/1001.0, 3.0/143.0, -2.0/143.0, 0, 0, 0, 2.0/143.0,
   -4.0/143.0, 4.0/143.0, -5.0/429.0, -5.0/429.0, 4.0/143.0, -4.0/143.0, 2.0/143.0, 0,
   0, 0, -3.0/221.0, 7.0/221.0, -100.0/2431.0, 84.0/2431.0, -98.0/7293.0, -98.0/7293.0,
   84.0/2431.0, -100.0/2431.0, 7.0/221.0, -3.0/221.0, 0, 4.0/39.0, -2.0/39.0, -2.0/143.0,
   20.0/429.0, -20.0/429.0, 1.0/33.0, -17.0/1287.0, 4.0/1287.0, 0, -7.0/51.0, 5.0/51.0,
   -4.0/221.0, -28.0/663.0, 490.0/7293.0, -14.0/221.0, 980.0/21879.0, -532.0/21879.0, 23.0/2431.0, -5.0/2431.0,
   0, 56.0/323.0, -49.0/323.0, 21.0/323.0, 6.0/323.0, -300.0/4199.0, 378.0/4199.0, -294.0/3553.0,
   168.0/2717.0, -1764.0/46189.0, 875.0/46189.0, -29.0/4199.0, 6.0/4199.0, -35.0/17.0, 56.0/51.0, -28.0/51.0,
   56.0/221.0, -70.0/663.0, 280.0/7293.0, -28.0/2431.0, 56.0/21879.0, -7.0/21879.0, 56.0/19.0, -560.0/323.0,
   315.0/323.0, -168.0/323.0, 84.0/323.0, -504.0/4199.0, 210.0/4199.0, -840.0/46189.0, 252.0/46189.0, -56.0/46189.0,
   7.0/46189.0, -4.0, 48.0/19.0, -88.0/57.0, 880.0/969.0, -165.0/323.0, 88.0/323.0, -44.0/323.0,
   264.0/4199.0, -110.0/4199.0, 40.0/4199.0, -12.0/4199.0, 8.0/12597.0, -1.0/12597.0, 0, 0,
   0, 0, 0, 0, 1.0/9.0, 0, 0, 0,
   0, 0, 0, 0, 0, -1.0/99.0, 7.0/99.0, -28.0/99.0,
   0, 0, 0, 0, 0, 0, 0, 0,
   1.0/429.0, -7.0/429.0, 28.0/429.0, -28.0/143.0, 70.0/143.0, 0, 0, 0,
   0, 0, 0, -2.0/231.0, 1.0/77.0, 8.0/693.0, -5.0/99.0, -4.0/99.0,
   0, 0, 0, 0, 5.0/1001.0, -40.0/3003.0, 43.0/3003.0, 2.0/429.0,
   -6.0/143.0, 8.0/143.0, 10.0/143.0, 0, 0, 0, 0, -4.0/1287.0,
   1.0/99.0, -2.0/117.0, 7.0/429.0, 0, -14.0/429.0, 28.0/429.0, -2.0/39.0, -4.0/39.0,
   0, 0, 2.0/143.0, -1.0/143.0, -1.0/91.0, 20.0/1001.0, -1.0/91.0, -1.0/143.0,
   2.0/143.0, 0, 0, -2.0/143.0, 2.0/143.0, 0, -7.0/429.0, 10.0/429.0,
   -7.0/429.0, 0, 2.0/143.0, -2.0/143.0, 0, 0, 3.0/221.0, -4.0/221.0,
   23.0/2431.0, 16.0/2431.0, -14.0/663.0, 196.0/7293.0, -14.0/663.0, 16.0/2431.0, 23.0/2431.0, -4.0/221.0,
   3.0/221.0, -4.0/39.0, -2.0/39.0, 28.0/429.0, -14.0/429.0, 0, 7.0/429.0, -2.0/117.0,
   1.0/99.0, -4.0/1287.0, 7.0/51.0, 2.0/51.0, -53.0/663.0, 40.0/663.0, -14.0/561.0, -28.0/7293.0,
   406.0/21879.0, -448.0/21879.0, 25.0/1683.0, -18.0/2431.0, 5.0/2431.0, -56.0/323.0, -7.0/323.0, 28.0/323.0,
   -27.0/323.0, 222.0/4199.0, -6.0/323.0, -336.0/46189.0, 966.0/46189.0, -84.0/3553.0, 889.0/46189.0, -556.0/46189.0,
   23.0/4199.0, -6.0/4199.0, 49.0/51.0, -28.0/51.0, 196.0/663.0, -98.0/663.0, 490.0/7293.0, -196.0/7293.0,
   196.0/21879.0, -49.0/21879.0, 7.0/21879.0, -392.0/323.0, 245.0/323.0, -147.0/323.0, 84.0/323.0, -588.0/4199.0,
   294.0/4199.0, -1470.0/46189.0, 588.0/46189.0, -196.0/46189.0, 49.0/46189.0, -7.0/46189.0, 28.0/19.0, -56.0/57.0,
   616.0/969.0, -385.0/969.0, 77.0/323.0, -44.0/323.0, 308.0/4199.0, -154.0/4199.0, 70.0/4199.0, -28.0/4199.0,
   28.0/12597.0, -7.0/12597.0, 1.0/12597.0, 0, 0, 0, 0, 0,
   -1.0/9.0, 0, 0, 0, 0, 0, 0, 0,
   0, 1.0/99.0, -2.0/33.0, 7.0/33.0, 0, 0, 0, 0,
   0, 0, 0, 0, -1.0/429.0, 2.0/143.0, -7.0/143.0, 56.0/429.0,
   -42.0/143.0, 0, 0, 0, 0, 0, 0, 2.0/231.0,
   -1.0/231.0, -17.0/693.0, 3.0/77.0, 1.0/11.0, 0, 0, 0, 0,
   -5.0/1001.0, 25.0/3003.0, -1.0/1001.0, -19.0/1001.0, 16.0/429.0, -2.0/143.0, -18.0/143.0, 0,
   0, 0, 0, 4.0/1287.0, -1.0/143.0, 1.0/143.0, 1.0/1287.0, -7.0/429.0,
   14.0/429.0, -14.0/429.0, -2.0/143.0, 2.0/13.0, 0, 0, -2.0/143.0, -1.0/143.0,
   18.0/1001.0, -9.0/1001.0, -9.0/1001.0, 18.0/1001.0, -1.0/143.0, -2.0/143.0, 0, 2.0/143.0,
   0, -2.0/143.0, 7.0/429.0, -1.0/143.0, -1.0/143.0, 7.0/429.0, -2.0/143.0, 0,
   2.0/143.0, 0, -3.0/221.0, 1.0/221.0, 21.0/2431.0, -3.0/187.0, 106.0/7293.0, -14.0/2431.0,
   -14.0/2431.0, 106.0/7293.0, -3.0/187.0, 21.0/2431.0, 1.0/221.0, -3.0/221.0, 2.0/13.0, -2.0/143.0,
   -14.0/429.0, 14.0/429.0, -7.0/429.0, 1.0/1287.0, 1.0/143.0, -1.0/143.0, 4.0/1287.0, -3.0/17.0,
   9.0/221.0, 1.0/51.0, -86.0/2431.0, 70.0/2431.0, -322.0/21879.0, 14.0/7293.0, 41.0/7293.0, -163.0/21879.0,
   1.0/187.0, -5.0/2431.0, 63.0/323.0, -21.0/323.0, -1.0/323.0, 129.0/4199.0, -144.0/4199.0, 1194.0/46189.0,
   -630.0/46189.0, 126.0/46189.0, 203.0/46189.0, -333.0/46189.0, 303.0/46189.0, -1.0/247.0, 6.0/4199.0, -7.0/17.0,
   56.0/221.0, -98.0/663.0, 196.0/2431.0, -98.0/2431.0, 392.0/21879.0, -49.0/7293.0, 14.0/7293.0, -7.0/21879.0,
   147.0/323.0, -98.0/323.0, 63.0/323.0, -504.0/4199.0, 294.0/4199.0, -1764.0/46189.0, 882.0/46189.0, -392.0/46189.0,
   147.0/46189.0, -42.0/46189.0, 7.0/46189.0, -28.0/57.0, 112.0/323.0, -77.0/323.0, 154.0/969.0, -33.0/323.0,
   264.0/4199.0, -154.0/4199.0, 84.0/4199.0, -42.0/4199.0, 56.0/12597.0, -7.0/4199.0, 2.0/4199.0, -1.0/12597.0,
   0, 0, 0, 0, 1.0/9.0, 0, 0, 0,
   0, 0, 0, 0, 0, -1.0/99.0, 5.0/99.0, -5.0/33.0,
   0, 0, 0, 0, 0, 0, 0, 0,
   1.0/429.0, -5.0/429.0, 5.0/143.0, -35.0/429.0, 70.0/429.0, 0, 0, 0,
   0, 0, 0, -2.0/231.0, -1.0/231.0, 20.0/693.0, -10.0/693.0, -10.0/77.0,
   0, 0, 0, 0, 5.0/1001.0, -10.0/3003.0, -2.0/273.0, 20.0/1001.0,
   -5.0/273.0, -10.0/429.0, 20.0/143.0, 0, 0, 0, 0, -4.0/1287.0,
   5.0/1287.0, 0, -10.0/1287.0, 20.0/1287.0, -7.0/429.0, 0, 20.0/429.0, -20.0/143.0,
   0, 0, 2.0/143.0, 3.0/143.0, -1.0/91.0, -9.0/1001.0, 18.0/1001.0, -9.0/1001.0,
   -1.0/91.0, 3.0/143.0, 2.0/143.0, -2.0/143.0, -2.0/143.0, 2.0/143.0, -1.0/429.0, -4.0/429.0,
   2.0/143.0, -4.0/429.0, -1.0/429.0, 2.0/143.0, -2.0/143.0, -2.0/143.0, 3.0/221.0, 2.0/221.0,
   -32.0/2431.0, 18.0/2431.0, 1.0/663.0, -64.0/7293.0, 28.0/2431.0, -64.0/7293.0, 1.0/663.0, 18.0/2431.0,
   -32.0/2431.0, 2.0/221.0, 3.0/221.0, -20.0/143.0, 20.0/429.0, 0, -7.0/429.0, 20.0/1287.0,
   -10.0/1287.0, 0, 5.0/1287.0, -4.0/1287.0, 30.0/221.0, -40.0/663.0, 115.0/7293.0, 16.0/2431.0,
   -28.0/1989.0, 280.0/21879.0, -5.0/663.0, 40.0/21879.0, 46.0/21879.0, -8.0/2431.0, 5.0/2431.0, -42.0/323.0,
   22.0/323.0, -116.0/4199.0, 15.0/4199.0, 30.0/3553.0, -564.0/46189.0, 504.0/46189.0, -329.0/46189.0, 10.0/3553.0,
   30.0/46189.0, -116.0/46189.0, 11.0/4199.0, -6.0/4199.0, 35.0/221.0, -70.0/663.0, 490.0/7293.0, -98.0/2431.0,
   490.0/21879.0, -245.0/21879.0, 35.0/7293.0, -35.0/21879.0, 7.0/21879.0, -49.0/323.0, 35.0/323.0, -315.0/4199.0,
   210.0/4199.0, -1470.0/46189.0, 882.0/46189.0, -490.0/46189.0, 245.0/46189.0, -105.0/46189.0, 35.0/46189.0, -7.0/46189.0,
   140.0/969.0, -35.0/323.0, 77.0/969.0, -55.0/969.0, 165.0/4199.0, -110.0/4199.0, 70.0/4199.0, -42.0/4199.0,
   70.0/12597.0, -35.0/12597.0, 5.0/4199.0, -5.0/12597.0, 1.0/12597.0, 0, 0, 0,
   -1.0/9.0, 0, 0, 0, 0, 0, 0, 0,
   0, 1.0/99.0, -4.0/99.0, 10.0/99.0, 0, 0, 0, 0,
   0, 0, 0, 0, -1.0/429.0, 4.0/429.0, -10.0/429.0, 20.0/429.0,
   -35.0/429.0, 0, 0, 0, 0, 0, 0, 2.0/231.0,
   1.0/77.0, -17.0/693.0, -10.0/693.0, 100.0/693.0, 0, 0, 0, 0,
   -5.0/1001.0, -5.0/3003.0, 32.0/3003.0, -38.0/3003.0, -5.0/3003.0, 125.0/3003.0, -50.0/429.0, 0,
   0, 0, 0, 4.0/1287.0, -1.0/1287.0, -5.0/1287.0, 10.0/1287.0, -10.0/1287.0,
   1.0/1287.0, 7.0/429.0, -20.0/429.0, 40.0/429.0, 0, 0, 0, -5.0/143.0,
   -10.0/1001.0, 20.0/1001.0, -9.0/1001.0, -9.0/1001.0, 20.0/1001.0, -10.0/1001.0, -5.0/143.0, 0,
   4.0/143.0, 0, -5.0/429.0, 5.0/429.0, -2.0/429.0, -2.0/429.0, 5.0/429.0, -5.0/429.0,
   0, 4.0/143.0, 0, -5.0/221.0, 10.0/2431.0, 14.0/2431.0, -5.0/561.0, 53.0/7293.0,
   -20.0/7293.0, -20.0/7293.0, 53.0/7293.0, -5.0/561.0, 14.0/2431.0, 10.0/2431.0, -5.0/221.0, 0,
   40.0/429.0, -20.0/429.0, 7.0/429.0, 1.0/1287.0, -10.0/1287.0, 10.0/1287.0, -5.0/1287.0, -1.0/1287.0,
   4.0/1287.0, -50.0/663.0, 25.0/561.0, -163.0/7293.0, 164.0/21879.0, 28.0/21879.0, -115.0/21879.0, 125.0/21879.0,
   -86.0/21879.0, 2.0/1683.0, 3.0/2431.0, -5.0/2431.0, 20.0/323.0, -10.0/247.0, 101.0/4199.0, -555.0/46189.0,
   174.0/46189.0, 60.0/46189.0, -175.0/46189.0, 199.0/46189.0, -160.0/46189.0, 86.0/46189.0, -5.0/46189.0, -5.0/4199.0,
   6.0/4199.0, -35.0/663.0, 280.0/7293.0, -196.0/7293.0, 392.0/21879.0, -245.0/21879.0, 140.0/21879.0, -70.0/21879.0,
   28.0/21879.0, -7.0/21879.0, 14.0/323.0, -140.0/4199.0, 105.0/4199.0, -840.0/46189.0, 588.0/46189.0, -392.0/46189.0,
   245.0/46189.0, -140.0/46189.0, 70.0/46189.0, -28.0/46189.0, 7.0/46189.0, -35.0/969.0, 28.0/969.0, -22.0/969.0,
   220.0/12597.0, -55.0/4199.0, 40.0/4199.0, -28.0/4199.0, 56.0/12597.0, -35.0/12597.0, 20.0/12597.0, -10.0/12597.0,
   4.0/12597.0, -1.0/12597.0, 0, 0, 1.0/9.0, 0, 0, 0,
   0, 0, 0, 0, 0, -1.0/99.0, 1.0/33.0, -2.0/33.0,
   0, 0, 0, 0, 0, 0, 0, 0,
   1.0/429.0, -1.0/143.0, 2.0/143.0, -10.0/429.0, 5.0/143.0, 0, 0, 0,
   0, 0, 0, -2.0/231.0, -5.0/231.0, 8.0/693.0, 3.0/77.0, -10.0/77.0,
   0, 0, 0, 0, 5.0/1001.0, 20.0/3003.0, -9.0/1001.0, 2.0/1001.0,
   43.0/3003.0, -40.0/1001.0, 75.0/1001.0, 0, 0, 0, 0, -4.0/1287.0,
   -1.0/429.0, 2.0/429.0, -5.0/1287.0, 0, 1.0/143.0, -2.0/117.0, 1.0/33.0, -20.0/429.0,
   0, 0, 0, 0, 45.0/1001.0, -10.0/1001.0, -1.0/91.0, 18.0/1001.0,
   -1.0/91.0, -10.0/1001.0, 45.0/1001.0, 0, 0, -4.0/143.0, 5.0/429.0, 0,
   -1.0/143.0, 4.0/429.0, -1.0/143.0, 0, 5.0/429.0, -4.0/143.0, 0, 0,
   45.0/2431.0, -24.0/2431.0, 23.0/7293.0, 4.0/2431.0, -1.0/221.0, 40.0/7293.0, -1.0/221.0, 4.0/2431.0,
   23.0/7293.0, -24.0/2431.0, 45.0/2431.0, 0, 0, -20.0/429.0, 1.0/33.0, -2.0/117.0,
   1.0/143.0, 0, -5.0/1287.0, 2.0/429.0, -1.0/429.0, -4.0/1287.0, 75.0/2431.0, -54.0/2431.0,
   25.0/1683.0, -64.0/7293.0, 29.0/7293.0, -10.0/21879.0, -1.0/561.0, 20.0/7293.0, -53.0/21879.0, 2.0/2431.0,
   5.0/2431.0, -90.0/4199.0, 69.0/4199.0, -556.0/46189.0, 381.0/46189.0, -18.0/3553.0, 115.0/46189.0, -24.0/46189.0,
   -3.0/3553.0, 74.0/46189.0, -81.0/46189.0, 60.0/46189.0, -1.0/4199.0, -6.0/4199.0, 35.0/2431.0, -28.0/2431.0,
   196.0/21879.0, -49.0/7293.0, 35.0/7293.0, -70.0/21879.0, 14.0/7293.0, -7.0/7293.0, 7.0/21879.0, -42.0/4199.0,
   35.0/4199.0, -315.0/46189.0, 252.0/46189.0, -196.0/46189.0, 147.0/46189.0, -105.0/46189.0, 70.0/46189.0, -42.0/46189.0,
   21.0/46189.0, -7.0/46189.0, 7.0/969.0, -2.0/323.0, 22.0/4199.0, -55.0/12597.0, 15.0/4199.0, -12.0/4199.0,
   28.0/12597.0, -7.0/4199.0, 5.0/4199.0, -10.0/12597.0, 2.0/4199.0, -1.0/4199.0, 1.0/12597.0, 0,
   -1.0/9.0, 0, 0, 0, 0, 0, 0, 0,
   0, 1.0/99.0, -2.0/99.0, 1.0/33.0, 0, 0, 0, 0,
   0, 0, 0, 0, -1.0/429.0, 2.0/429.0, -1.0/143.0, 4.0/429.0,
   -5.0/429.0, 0, 0, 0, 0, 0, 0, 0,
   1.0/33.0, 1.0/99.0, -5.0/99.0, 1.0/11.0, 0, 0, 0, 0,
   0, -5.0/429.0, 1.0/429.0, 1.0/143.0, -7.0/429.0, 1.0/39.0, -5.0/143.0, 0,
   0, 0, 0, 0, 7.0/1287.0, -1.0/429.0, -1.0/1287.0, 5.0/1287.0,
   -1.0/143.0, 1.0/99.0, -17.0/1287.0, 7.0/429.0, 0, 0, 0, 0,
   0, -5.0/143.0, 3.0/143.0, -1.0/143.0, -1.0/143.0, 3.0/143.0, -5.0/143.0, 0,
   0, 0, 7.0/429.0, -5.0/429.0, 1.0/143.0, -1.0/429.0, -1.0/429.0, 1.0/143.0,
   -5.0/429.0, 7.0/429.0, 0, 0, 0, -21.0/2431.0, 49.0/7293.0, -35.0/7293.0,
   7.0/2431.0, -7.0/7293.0, -7.0/7293.0, 7.0/2431.0, -35.0/7293.0, 49.0/7293.0, -21.0/2431.0, 0,
   0, 0, 7.0/429.0, -17.0/1287.0, 1.0/99.0, -1.0/143.0, 5.0/1287.0, -1.0/1287.0,
   -1.0/429.0, 7.0/1287.0, 0, -21.0/2431.0, 161.0/21879.0, -133.0/21879.0, 35.0/7293.0, -7.0/1989.0,
   49.0/21879.0, -7.0/7293.0, -7.0/21879.0, 35.0/21879.0, -7.0/2431.0, 0, 21.0/4199.0, -203.0/46189.0,
   175.0/46189.0, -147.0/46189.0, 7.0/2717.0, -7.0/3553.0, 63.0/46189.0, -35.0/46189.0, 7.0/46189.0, 21.0/46189.0,
   -49.0/46189.0, 7.0/4199.0, 0, -7.0/2431.0, 56.0/21879.0, -49.0/21879.0, 14.0/7293.0, -35.0/21879.0,
   28.0/21879.0, -7.0/7293.0, 14.0/21879.0, -7.0/21879.0, 7.0/4199.0, -70.0/46189.0, 63.0/46189.0, -56.0/46189.0,
   49.0/46189.0, -42.0/46189.0, 35.0/46189.0, -28.0/46189.0, 21.0/46189.0, -14.0/46189.0, 7.0/46189.0, -1.0/969.0,
   4.0/4199.0, -11.0/12597.0, 10.0/12597.0, -3.0/4199.0, 8.0/12597.0, -7.0/12597.0, 2.0/4199.0, -5.0/12597.0,
   4.0/12597.0, -1.0/4199.0, 2.0/12597.0, -1.0/12597.0, 1.0/9.0, 0, 0, 0,
   0, 0, 0, 0, 0, -1.0/99.0, 1.0/99.0, -1.0/99.0,
   0, 0, 0, 0, 0, 0, 0, 0,
   1.0/429.0, -1.0/429.0, 1.0/429.0, -1.0/429.0, 1.0/429.0, 0, 0, 0,
   0, 0, 0, 0, 0, -4.0/99.0, 4.0/99.0, -4.0/99.0,
   0, 0, 0, 0, 0, 0, 4.0/429.0, -4.0/429.0,
   4.0/429.0, -4.0/429.0, 4.0/429.0, 0, 0, 0, 0, 0,
   0, -4.0/1287.0, 4.0/1287.0, -4.0/1287.0, 4.0/1287.0, -4.0/1287.0, 4.0/1287.0, -4.0/1287.0,
   0, 0, 0, 0, 0, 0, 2.0/143.0, -2.0/143.0,
   2.0/143.0, -2.0/143.0, 2.0/143.0, 0, 0, 0, 0, -2.0/429.0,
   2.0/429.0, -2.0/429.0, 2.0/429.0, -2.0/429.0, 2.0/429.0, -2.0/429.0, 0, 0,
   0, 0, 14.0/7293.0, -14.0/7293.0, 14.0/7293.0, -14.0/7293.0, 14.0/7293.0, -14.0/7293.0,
   14.0/7293.0, -14.0/7293.0, 14.0/7293.0, 0, 0, 0, 0, -4.0/1287.0,
   4.0/1287.0, -4.0/1287.0, 4.0/1287.0, -4.0/1287.0, 4.0/1287.0, -4.0/1287.0, 0, 0,
   28.0/21879.0, -28.0/21879.0, 28.0/21879.0, -28.0/21879.0, 28.0/21879.0, -28.0/21879.0, 28.0/21879.0, -28.0/21879.0,
   28.0/21879.0, 0, 0, -28.0/46189.0, 28.0/46189.0, -28.0/46189.0, 28.0/46189.0, -28.0/46189.0,
   28.0/46189.0, -28.0/46189.0, 28.0/46189.0, -28.0/46189.0, 28.0/46189.0, -28.0/46189.0, 0, 0,
   7.0/21879.0, -7.0/21879.0, 7.0/21879.0, -7.0/21879.0, 7.0/21879.0, -7.0/21879.0, 7.0/21879.0, -7.0/21879.0,
   7.0/21879.0, -7.0/46189.0, 7.0/46189.0, -7.0/46189.0, 7.0/46189.0, -7.0/46189.0, 7.0/46189.0, -7.0/46189.0,
   7.0/46189.0, -7.0/46189.0, 7.0/46189.0, -7.0/46189.0, 1.0/12597.0, -1.0/12597.0, 1.0/12597.0, -1.0/12597.0,
   1.0/12597.0, -1.0/12597.0, 1.0/12597.0, -1.0/12597.0, 1.0/12597.0, -1.0/12597.0, 1.0/12597.0, -1.0/12597.0,
   1.0/12597.0
};
static unsigned const s_i3RlmSphereInts[5] = {
   0, 49, 337, 1012, 2132
};

static unsigned _3RlmSphereIntIndex(int la, int ma, int lb, int mb, int lamb)
{
   int const MaxLa = 6;
   // if this is not a pretty indexing function, I do not know what is...
   assert(0 <= la && la <= MaxLa && 0 <= lb && 0 <= lamb);
   assert(lb <= la);
   assert(la - lb <= lamb && lamb <= la + lb);
   assert(-la <= ma && ma <= la);
   assert(-lb <= mb && mb <= lb);
   int
      lb2 = sqr(lb),
      nSlmX_MaxA = sqr(MaxLa+1);
   int
      idx;
   // find starting index of data under index lb (all of it!)
   // idx = -lb*(lb + 1)*(nSlmX_MaxA*(10 - 40*lb) + 6 + lb*(9 - lb*(39 - 24*lb)))/60;
   // ...same quantity
   assert(unsigned(lb) < sizeof(s_i3RlmSphereInts)/sizeof(s_i3RlmSphereInts[0]));
   idx = s_i3RlmSphereInts[lb];

   // add index for (la,ma) components, for la >= lb.
   idx += (sqr(la) + la + ma - lb2);  // ... = (iSlmX(la,ma) - iSlmX(lb));

   // get total number of (la,ma) components for lb <= la <= MaxLa
   int nCompA = (nSlmX_MaxA - lb2);    // ... = nSlmX(MaxLa,lb) = iSlmX(MaxLa+1) - iSlmX(lb) = nSlmx(MaxLa) - nSlmX(lb-1)

   // add final index for mb and lambda components
   idx += nCompA * ((lb + lamb - la)/2 + (lb + 1)*(lb + mb));
   assert(idx >= 0 && size_t(idx) < sizeof(s_3RlmSphereInts)/sizeof(s_3RlmSphereInts[0]));
   return unsigned(idx);

}


// Return sphere surface integrals over complex harmonics:
//    int R[La,Ma] R[EcpL,EcpM] R[Lambda,Mu] dOmega
// Notes:
//   - Leading 4pi factor is omitted
//   - The integral itself is invariant to permutations of the Rlms, but this function is not.
//     The range of tabulated values differs. We have:
//        0 <= La <= 6"
//        0 <= EcpL <= 4"
//        La-EcpL <= Lambda <= La+EcpL
//   - This integral is manifestly zero unless Ma + EcpM + Mu == 0
//   - As a result, all non-vanishing integrals are real, as the complex phases are only in the phi integral
//   - There are still a bunch of other zeros around. So it may help to check for that.
double GetSphereIntegral3Rlm(int La, int Ma, int EcpL, int EcpM, int Lambda, int Mu)
{
   assert(La >= 0 && -La <= Ma && Ma <= La);
   assert(EcpL >= 0 && -EcpL <= EcpM && EcpM <= EcpL);
   assert(Lambda >= 0 && -Lambda <= Mu && Mu <= Lambda);
   if (Ma + EcpM + Mu != 0) return 0.;
   if (Lambda > La + EcpL) return 0.;
   if (EcpL > La) {
      std::swap(EcpL,La);
      std::swap(EcpM,Ma);
   }
   if (Lambda < La - EcpL) return 0.;
   if (EcpM < -EcpL || EcpM > EcpL) return 0.;
   if ((unsigned(Lambda + EcpL + La) & 1) != 0) return 0.;
   if (La > 6 || EcpL > 4)
      throw std::runtime_error("IR/ECP calculation: ran out of angular momenta.");
   return s_3RlmSphereInts[_3RlmSphereIntIndex(La, Ma, EcpL, EcpM, Lambda)];
}

} // namespace ir_rlm
