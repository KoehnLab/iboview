/* CxAngularGrid_Orbits.inl v20170313 EST [storm, Gerald Knizia] */
static void make_orbit_oct_g(double *IR_RP p, double x, double y, double z, double w)
{
   // octahedral general orbit, length 24
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = x;    p[ 5] = -y;   p[ 6] = -z;   p[ 7] = w;
   p[ 8] = -x;   p[ 9] = y;    p[10] = -z;   p[11] = w;
   p[12] = -x;   p[13] = -y;   p[14] = z;    p[15] = w;
   p[16] = z;    p[17] = x;    p[18] = y;    p[19] = w;
   p[20] = -z;   p[21] = -x;   p[22] = y;    p[23] = w;
   p[24] = -z;   p[25] = x;    p[26] = -y;   p[27] = w;
   p[28] = z;    p[29] = -x;   p[30] = -y;   p[31] = w;
   p[32] = y;    p[33] = z;    p[34] = x;    p[35] = w;
   p[36] = -y;   p[37] = z;    p[38] = -x;   p[39] = w;
   p[40] = -y;   p[41] = -z;   p[42] = x;    p[43] = w;
   p[44] = y;    p[45] = -z;   p[46] = -x;   p[47] = w;
   p[48] = x;    p[49] = -z;   p[50] = y;    p[51] = w;
   p[52] = x;    p[53] = z;    p[54] = -y;   p[55] = w;
   p[56] = z;    p[57] = y;    p[58] = -x;   p[59] = w;
   p[60] = -z;   p[61] = y;    p[62] = x;    p[63] = w;
   p[64] = -y;   p[65] = x;    p[66] = z;    p[67] = w;
   p[68] = y;    p[69] = -x;   p[70] = z;    p[71] = w;
   p[72] = y;    p[73] = x;    p[74] = -z;   p[75] = w;
   p[76] = -y;   p[77] = -x;   p[78] = -z;   p[79] = w;
   p[80] = z;    p[81] = -y;   p[82] = x;    p[83] = w;
   p[84] = -z;   p[85] = -y;   p[86] = -x;   p[87] = w;
   p[88] = -x;   p[89] = z;    p[90] = y;    p[91] = w;
   p[92] = -x;   p[93] = -z;   p[94] = -y;   p[95] = w;
}

static void make_orbit_oct_coe(double *IR_RP p, double x, double y, double z, double w)
{
   // octahedral center of edge orbit, length 12
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = x;    p[ 5] = -y;   p[ 6] = -z;   p[ 7] = w;
   p[ 8] = -x;   p[ 9] = y;    p[10] = -z;   p[11] = w;
   p[12] = -x;   p[13] = -y;   p[14] = z;    p[15] = w;
   p[16] = z;    p[17] = x;    p[18] = y;    p[19] = w;
   p[20] = -z;   p[21] = -x;   p[22] = y;    p[23] = w;
   p[24] = -z;   p[25] = x;    p[26] = -y;   p[27] = w;
   p[28] = z;    p[29] = -x;   p[30] = -y;   p[31] = w;
   p[32] = y;    p[33] = z;    p[34] = x;    p[35] = w;
   p[36] = -y;   p[37] = z;    p[38] = -x;   p[39] = w;
   p[40] = -y;   p[41] = -z;   p[42] = x;    p[43] = w;
   p[44] = y;    p[45] = -z;   p[46] = -x;   p[47] = w;
}

static void make_orbit_oct_cof(double *IR_RP p, double x, double y, double z, double w)
{
   // octahedral center of face orbit, length 6
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = -x;   p[ 5] = y;    p[ 6] = -z;   p[ 7] = w;
   p[ 8] = z;    p[ 9] = x;    p[10] = y;    p[11] = w;
   p[12] = -z;   p[13] = -x;   p[14] = y;    p[15] = w;
   p[16] = y;    p[17] = z;    p[18] = x;    p[19] = w;
   p[20] = -y;   p[21] = z;    p[22] = -x;   p[23] = w;
}

static void make_orbit_oct_v(double *IR_RP p, double x, double y, double z, double w)
{
   // octahedral vertex orbit, length 8
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = x;    p[ 5] = -y;   p[ 6] = -z;   p[ 7] = w;
   p[ 8] = -x;   p[ 9] = y;    p[10] = -z;   p[11] = w;
   p[12] = -x;   p[13] = -y;   p[14] = z;    p[15] = w;
   p[16] = x;    p[17] = -z;   p[18] = y;    p[19] = w;
   p[20] = x;    p[21] = z;    p[22] = -y;   p[23] = w;
   p[24] = -z;   p[25] = y;    p[26] = x;    p[27] = w;
   p[28] = -y;   p[29] = -x;   p[30] = -z;   p[31] = w;
}

static void make_orbit_tet_g(double *IR_RP p, double x, double y, double z, double w)
{
   // tetrahedral general orbit, length 12
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = z;    p[ 5] = x;    p[ 6] = y;    p[ 7] = w;
   p[ 8] = -z;   p[ 9] = -x;   p[10] = y;    p[11] = w;
   p[12] = -z;   p[13] = x;    p[14] = -y;   p[15] = w;
   p[16] = z;    p[17] = -x;   p[18] = -y;   p[19] = w;
   p[20] = y;    p[21] = z;    p[22] = x;    p[23] = w;
   p[24] = -y;   p[25] = z;    p[26] = -x;   p[27] = w;
   p[28] = -y;   p[29] = -z;   p[30] = x;    p[31] = w;
   p[32] = y;    p[33] = -z;   p[34] = -x;   p[35] = w;
   p[36] = x;    p[37] = -y;   p[38] = -z;   p[39] = w;
   p[40] = -x;   p[41] = y;    p[42] = -z;   p[43] = w;
   p[44] = -x;   p[45] = -y;   p[46] = z;    p[47] = w;
}

static void make_orbit_tet_coe(double *IR_RP p, double x, double y, double z, double w)
{
   // tetrahedral center of edge orbit, length 6
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = z;    p[ 5] = x;    p[ 6] = y;    p[ 7] = w;
   p[ 8] = -z;   p[ 9] = -x;   p[10] = y;    p[11] = w;
   p[12] = y;    p[13] = z;    p[14] = x;    p[15] = w;
   p[16] = -y;   p[17] = -z;   p[18] = x;    p[19] = w;
   p[20] = x;    p[21] = -y;   p[22] = -z;   p[23] = w;
}

static void make_orbit_tet_cof(double *IR_RP p, double x, double y, double z, double w)
{
   // tetrahedral center of face orbit, length 4
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = z;    p[ 5] = x;    p[ 6] = y;    p[ 7] = w;
   p[ 8] = -z;   p[ 9] = x;    p[10] = -y;   p[11] = w;
   p[12] = z;    p[13] = -x;   p[14] = -y;   p[15] = w;
}

static void make_orbit_tet_v(double *IR_RP p, double x, double y, double z, double w)
{
   // tetrahedral vertex orbit, length 4
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = -z;   p[ 5] = -x;   p[ 6] = y;    p[ 7] = w;
   p[ 8] = -z;   p[ 9] = x;    p[10] = -y;   p[11] = w;
   p[12] = z;    p[13] = -x;   p[14] = -y;   p[15] = w;
}

static void make_orbit_ico_g(double *IR_RP p, double x, double y, double z, double w)
{
   // icosahedral general orbit, length 60
   double const b = 0.809016994374947424102293; // beta = 1/4 + sqrt(5)/4
   double const g = -0.309016994374947424102293; // gamma = -sqrt(5)/4 + 1/4
   double const bx = b*x;
   double const by = b*y;
   double const bz = b*z;
   double const gx = g*x;
   double const gy = g*y;
   double const gz = g*z;
   double const hx = x/2;
   double const hy = y/2;
   double const hz = z/2;
   p[  0] = x;               p[  1] = y;               p[  2] = z;               p[  3] = w;
   p[  4] = -x;              p[  5] = y;               p[  6] = -z;              p[  7] = w;
   p[  8] = -x;              p[  9] = -y;              p[ 10] = z;               p[ 11] = w;
   p[ 12] = x;               p[ 13] = -y;              p[ 14] = -z;              p[ 15] = w;
   p[ 16] = y;               p[ 17] = z;               p[ 18] = x;               p[ 19] = w;
   p[ 20] = -z;              p[ 21] = x;               p[ 22] = -y;              p[ 23] = w;
   p[ 24] = y;               p[ 25] = -z;              p[ 26] = -x;              p[ 27] = w;
   p[ 28] = z;               p[ 29] = x;               p[ 30] = y;               p[ 31] = w;
   p[ 32] = -y;              p[ 33] = -z;              p[ 34] = x;               p[ 35] = w;
   p[ 36] = -z;              p[ 37] = -x;              p[ 38] = y;               p[ 39] = w;
   p[ 40] = -y;              p[ 41] = z;               p[ 42] = -x;              p[ 43] = w;
   p[ 44] = z;               p[ 45] = -x;              p[ 46] = -y;              p[ 47] = w;
   p[ 48] = bz - gy - hx;    p[ 49] = -by - gx + hz;   p[ 50] = bx - gz + hy;    p[ 51] = w;
   p[ 52] = -bz - gy - hx;   p[ 53] = -by - gx - hz;   p[ 54] = -bx - gz - hy;   p[ 55] = w;
   p[ 56] = -bz + gy - hx;   p[ 57] = -by + gx + hz;   p[ 58] = -bx - gz + hy;   p[ 59] = w;
   p[ 60] = bz + gy - hx;    p[ 61] = -by + gx - hz;   p[ 62] = bx - gz - hy;    p[ 63] = w;
   p[ 64] = by - gx - hz;    p[ 65] = bx + gz - hy;    p[ 66] = -bz + gy - hx;   p[ 67] = w;
   p[ 68] = -by - gx + hz;   p[ 69] = -bx + gz - hy;   p[ 70] = -bz + gy + hx;   p[ 71] = w;
   p[ 72] = -by - gx - hz;   p[ 73] = -bx - gz - hy;   p[ 74] = -bz - gy - hx;   p[ 75] = w;
   p[ 76] = -bx - gz + hy;   p[ 77] = bz - gy + hx;    p[ 78] = by - gx - hz;    p[ 79] = w;
   p[ 80] = -bx + gz + hy;   p[ 81] = -bz - gy + hx;   p[ 82] = -by + gx - hz;   p[ 83] = w;
   p[ 84] = -bx + gz - hy;   p[ 85] = bz - gy - hx;    p[ 86] = by + gx - hz;    p[ 87] = w;
   p[ 88] = -bx - gz - hy;   p[ 89] = -bz - gy - hx;   p[ 90] = -by - gx - hz;   p[ 91] = w;
   p[ 92] = by - gx + hz;    p[ 93] = bx - gz - hy;    p[ 94] = -bz - gy + hx;   p[ 95] = w;
   p[ 96] = bz + gy - hx;    p[ 97] = by - gx + hz;    p[ 98] = -bx + gz + hy;   p[ 99] = w;
   p[100] = -bz + gy - hx;   p[101] = by - gx - hz;    p[102] = bx + gz - hy;    p[103] = w;
   p[104] = bx - gz + hy;    p[105] = -bz + gy + hx;   p[106] = by + gx - hz;    p[107] = w;
   p[108] = bx + gz - hy;    p[109] = -bz + gy - hx;   p[110] = by - gx - hz;    p[111] = w;
   p[112] = -by + gx + hz;   p[113] = bx + gz - hy;    p[114] = bz - gy + hx;    p[115] = w;
   p[116] = by + gx - hz;    p[117] = -bx + gz - hy;   p[118] = bz - gy - hx;    p[119] = w;
   p[120] = -by + gx - hz;   p[121] = bx - gz - hy;    p[122] = bz + gy - hx;    p[123] = w;
   p[124] = by + gx + hz;    p[125] = -bx - gz - hy;   p[126] = bz + gy + hx;    p[127] = w;
   p[128] = bx - gz - hy;    p[129] = bz + gy - hx;    p[130] = -by + gx - hz;   p[131] = w;
   p[132] = bx + gz + hy;    p[133] = bz + gy + hx;    p[134] = -by - gx - hz;   p[135] = w;
   p[136] = -bz - gy - hx;   p[137] = by + gx + hz;    p[138] = bx + gz + hy;    p[139] = w;
   p[140] = bz - gy - hx;    p[141] = by + gx - hz;    p[142] = -bx + gz - hy;   p[143] = w;
   p[144] = -by - gx + hz;   p[145] = bx - gz + hy;    p[146] = bz - gy - hx;    p[147] = w;
   p[148] = by - gx + hz;    p[149] = -bx + gz + hy;   p[150] = bz + gy - hx;    p[151] = w;
   p[152] = -by - gx - hz;   p[153] = bx + gz + hy;    p[154] = bz + gy + hx;    p[155] = w;
   p[156] = by - gx - hz;    p[157] = -bx - gz + hy;   p[158] = bz - gy + hx;    p[159] = w;
   p[160] = bx + gz - hy;    p[161] = bz - gy + hx;    p[162] = -by + gx + hz;   p[163] = w;
   p[164] = bx - gz + hy;    p[165] = bz - gy - hx;    p[166] = -by - gx + hz;   p[167] = w;
   p[168] = bz - gy + hx;    p[169] = by - gx - hz;    p[170] = -bx - gz + hy;   p[171] = w;
   p[172] = -bz + gy + hx;   p[173] = by + gx - hz;    p[174] = bx - gz + hy;    p[175] = w;
   p[176] = bz + gy + hx;    p[177] = by + gx + hz;    p[178] = -bx - gz - hy;   p[179] = w;
   p[180] = -bz - gy + hx;   p[181] = by - gx + hz;    p[182] = bx - gz - hy;    p[183] = w;
   p[184] = bx + gz + hy;    p[185] = -bz - gy - hx;   p[186] = by + gx + hz;    p[187] = w;
   p[188] = bx - gz - hy;    p[189] = -bz - gy + hx;   p[190] = by - gx + hz;    p[191] = w;
   p[192] = by + gx + hz;    p[193] = bx + gz + hy;    p[194] = -bz - gy - hx;   p[195] = w;
   p[196] = -by + gx - hz;   p[197] = -bx + gz + hy;   p[198] = -bz - gy + hx;   p[199] = w;
   p[200] = -bx - gz - hy;   p[201] = bz + gy + hx;    p[202] = by + gx + hz;    p[203] = w;
   p[204] = -bx - gz + hy;   p[205] = -bz + gy - hx;   p[206] = -by + gx + hz;   p[207] = w;
   p[208] = -bz + gy + hx;   p[209] = -by - gx + hz;   p[210] = -bx + gz - hy;   p[211] = w;
   p[212] = bz - gy + hx;    p[213] = -by + gx + hz;   p[214] = bx + gz - hy;    p[215] = w;
   p[216] = by + gx - hz;    p[217] = bx - gz + hy;    p[218] = -bz + gy + hx;   p[219] = w;
   p[220] = -by + gx + hz;   p[221] = -bx - gz + hy;   p[222] = -bz + gy - hx;   p[223] = w;
   p[224] = -bx + gz + hy;   p[225] = bz + gy - hx;    p[226] = by - gx + hz;    p[227] = w;
   p[228] = -bx + gz - hy;   p[229] = -bz + gy + hx;   p[230] = -by - gx + hz;   p[231] = w;
   p[232] = -bz - gy + hx;   p[233] = -by + gx - hz;   p[234] = -bx + gz + hy;   p[235] = w;
   p[236] = bz + gy + hx;    p[237] = -by - gx - hz;   p[238] = bx + gz + hy;    p[239] = w;
}

static void make_orbit_ico_coe(double *IR_RP p, double x, double y, double z, double w)
{
   // icosahedral center of edge orbit, length 30
   double const b = 0.809016994374947424102293; // beta = 1/4 + sqrt(5)/4
   double const g = -0.309016994374947424102293; // gamma = -sqrt(5)/4 + 1/4
   double const bx = b*x;
   double const by = b*y;
   double const bz = b*z;
   double const gx = g*x;
   double const gy = g*y;
   double const gz = g*z;
   double const hx = x/2;
   double const hy = y/2;
   double const hz = z/2;
   p[  0] = x;               p[  1] = y;               p[  2] = z;               p[  3] = w;
   p[  4] = -x;              p[  5] = y;               p[  6] = -z;              p[  7] = w;
   p[  8] = y;               p[  9] = z;               p[ 10] = x;               p[ 11] = w;
   p[ 12] = -z;              p[ 13] = x;               p[ 14] = -y;              p[ 15] = w;
   p[ 16] = y;               p[ 17] = -z;              p[ 18] = -x;              p[ 19] = w;
   p[ 20] = -z;              p[ 21] = -x;              p[ 22] = y;               p[ 23] = w;
   p[ 24] = bz - gy - hx;    p[ 25] = -by - gx + hz;   p[ 26] = bx - gz + hy;    p[ 27] = w;
   p[ 28] = -bz - gy - hx;   p[ 29] = -by - gx - hz;   p[ 30] = -bx - gz - hy;   p[ 31] = w;
   p[ 32] = -bz + gy - hx;   p[ 33] = -by + gx + hz;   p[ 34] = -bx - gz + hy;   p[ 35] = w;
   p[ 36] = bz + gy - hx;    p[ 37] = -by + gx - hz;   p[ 38] = bx - gz - hy;    p[ 39] = w;
   p[ 40] = by - gx - hz;    p[ 41] = bx + gz - hy;    p[ 42] = -bz + gy - hx;   p[ 43] = w;
   p[ 44] = -by - gx + hz;   p[ 45] = -bx + gz - hy;   p[ 46] = -bz + gy + hx;   p[ 47] = w;
   p[ 48] = -by - gx - hz;   p[ 49] = -bx - gz - hy;   p[ 50] = -bz - gy - hx;   p[ 51] = w;
   p[ 52] = -bx - gz + hy;   p[ 53] = bz - gy + hx;    p[ 54] = by - gx - hz;    p[ 55] = w;
   p[ 56] = -bx + gz + hy;   p[ 57] = -bz - gy + hx;   p[ 58] = -by + gx - hz;   p[ 59] = w;
   p[ 60] = -bx + gz - hy;   p[ 61] = bz - gy - hx;    p[ 62] = by + gx - hz;    p[ 63] = w;
   p[ 64] = -bx - gz - hy;   p[ 65] = -bz - gy - hx;   p[ 66] = -by - gx - hz;   p[ 67] = w;
   p[ 68] = by - gx + hz;    p[ 69] = bx - gz - hy;    p[ 70] = -bz - gy + hx;   p[ 71] = w;
   p[ 72] = bx - gz + hy;    p[ 73] = -bz + gy + hx;   p[ 74] = by + gx - hz;    p[ 75] = w;
   p[ 76] = bx + gz - hy;    p[ 77] = -bz + gy - hx;   p[ 78] = by - gx - hz;    p[ 79] = w;
   p[ 80] = -by + gx + hz;   p[ 81] = bx + gz - hy;    p[ 82] = bz - gy + hx;    p[ 83] = w;
   p[ 84] = by + gx - hz;    p[ 85] = -bx + gz - hy;   p[ 86] = bz - gy - hx;    p[ 87] = w;
   p[ 88] = -by + gx - hz;   p[ 89] = bx - gz - hy;    p[ 90] = bz + gy - hx;    p[ 91] = w;
   p[ 92] = by + gx + hz;    p[ 93] = -bx - gz - hy;   p[ 94] = bz + gy + hx;    p[ 95] = w;
   p[ 96] = bx - gz - hy;    p[ 97] = bz + gy - hx;    p[ 98] = -by + gx - hz;   p[ 99] = w;
   p[100] = bx + gz + hy;    p[101] = bz + gy + hx;    p[102] = -by - gx - hz;   p[103] = w;
   p[104] = bz - gy + hx;    p[105] = by - gx - hz;    p[106] = -bx - gz + hy;   p[107] = w;
   p[108] = -bz + gy + hx;   p[109] = by + gx - hz;    p[110] = bx - gz + hy;    p[111] = w;
   p[112] = bz + gy + hx;    p[113] = by + gx + hz;    p[114] = -bx - gz - hy;   p[115] = w;
   p[116] = -bz - gy + hx;   p[117] = by - gx + hz;    p[118] = bx - gz - hy;    p[119] = w;
}

static void make_orbit_ico_cof(double *IR_RP p, double x, double y, double z, double w)
{
   // icosahedral center of face orbit, length 20
   double const b = 0.809016994374947424102293; // beta = 1/4 + sqrt(5)/4
   double const g = -0.309016994374947424102293; // gamma = -sqrt(5)/4 + 1/4
   double const bx = b*x;
   double const by = b*y;
   double const bz = b*z;
   double const gx = g*x;
   double const gy = g*y;
   double const gz = g*z;
   double const hx = x/2;
   double const hy = y/2;
   double const hz = z/2;
   p[ 0] = x;                p[ 1] = y;                p[ 2] = z;                p[ 3] = w;
   p[ 4] = -x;               p[ 5] = y;                p[ 6] = -z;               p[ 7] = w;
   p[ 8] = -x;               p[ 9] = -y;               p[10] = z;                p[11] = w;
   p[12] = x;                p[13] = -y;               p[14] = -z;               p[15] = w;
   p[16] = y;                p[17] = z;                p[18] = x;                p[19] = w;
   p[20] = -z;               p[21] = x;                p[22] = -y;               p[23] = w;
   p[24] = y;                p[25] = -z;               p[26] = -x;               p[27] = w;
   p[28] = z;                p[29] = x;                p[30] = y;                p[31] = w;
   p[32] = -y;               p[33] = -z;               p[34] = x;                p[35] = w;
   p[36] = -z;               p[37] = -x;               p[38] = y;                p[39] = w;
   p[40] = -y;               p[41] = z;                p[42] = -x;               p[43] = w;
   p[44] = z;                p[45] = -x;               p[46] = -y;               p[47] = w;
   p[48] = -bz + gy - hx;    p[49] = -by + gx + hz;    p[50] = -bx - gz + hy;    p[51] = w;
   p[52] = bz + gy - hx;     p[53] = -by + gx - hz;    p[54] = bx - gz - hy;     p[55] = w;
   p[56] = by - gx - hz;     p[57] = bx + gz - hy;     p[58] = -bz + gy - hx;    p[59] = w;
   p[60] = -bx - gz + hy;    p[61] = bz - gy + hx;     p[62] = by - gx - hz;     p[63] = w;
   p[64] = -bx + gz + hy;    p[65] = -bz - gy + hx;    p[66] = -by + gx - hz;    p[67] = w;
   p[68] = by - gx + hz;     p[69] = bx - gz - hy;     p[70] = -bz - gy + hx;    p[71] = w;
   p[72] = bx + gz - hy;     p[73] = -bz + gy - hx;    p[74] = by - gx - hz;     p[75] = w;
   p[76] = bx - gz - hy;     p[77] = bz + gy - hx;     p[78] = -by + gx - hz;    p[79] = w;
}

static void make_orbit_ico_v(double *IR_RP p, double x, double y, double z, double w)
{
   // icosahedral vertex orbit, length 12
   p[ 0] = x;    p[ 1] = y;    p[ 2] = z;    p[ 3] = w;
   p[ 4] = -x;   p[ 5] = y;    p[ 6] = -z;   p[ 7] = w;
   p[ 8] = -x;   p[ 9] = -y;   p[10] = z;    p[11] = w;
   p[12] = x;    p[13] = -y;   p[14] = -z;   p[15] = w;
   p[16] = y;    p[17] = z;    p[18] = x;    p[19] = w;
   p[20] = -z;   p[21] = x;    p[22] = -y;   p[23] = w;
   p[24] = y;    p[25] = -z;   p[26] = -x;   p[27] = w;
   p[28] = z;    p[29] = x;    p[30] = y;    p[31] = w;
   p[32] = -y;   p[33] = -z;   p[34] = x;    p[35] = w;
   p[36] = -z;   p[37] = -x;   p[38] = y;    p[39] = w;
   p[40] = -y;   p[41] = z;    p[42] = -x;   p[43] = w;
   p[44] = z;    p[45] = -x;   p[46] = -y;   p[47] = w;
}
