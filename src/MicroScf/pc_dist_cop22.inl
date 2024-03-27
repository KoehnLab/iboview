// distance between center of positions (AB) and (CD)
void pc_dist_cop22(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom)
{
   index_t iAx = 3.0*iAtom[0] + 0;
   index_t iAy = 3.0*iAtom[0] + 1;
   index_t iAz = 3.0*iAtom[0] + 2;
   index_t iBx = 3.0*iAtom[1] + 0;
   index_t iBy = 3.0*iAtom[1] + 1;
   index_t iBz = 3.0*iAtom[1] + 2;
   index_t iCx = 3.0*iAtom[2] + 0;
   index_t iCy = 3.0*iAtom[2] + 1;
   index_t iCz = 3.0*iAtom[2] + 2;
   index_t iDx = 3.0*iAtom[3] + 0;
   index_t iDy = 3.0*iAtom[3] + 1;
   index_t iDz = 3.0*iAtom[3] + 2;
   double Ax = pXyz[iAx];
   double Ay = pXyz[iAy];
   double Az = pXyz[iAz];
   double Bx = pXyz[iBx];
   double By = pXyz[iBy];
   double Bz = pXyz[iBz];
   double Cx = pXyz[iCx];
   double Cy = pXyz[iCy];
   double Cz = pXyz[iCz];
   double Dx = pXyz[iDx];
   double Dy = pXyz[iDy];
   double Dz = pXyz[iDz];
   double r1 = Bx + Ax;
   double r2 = Dx + Cx;
   double r3 = 0.5*r1 - 0.5*r2;
   double r4 = By + Ay;
   double r5 = Dy + Cy;
   double r6 = 0.5*r4 - 0.5*r5;
   double r7 = Bz + Az;
   double r8 = Dz + Cz;
   double r9 = 0.5*r7 - 0.5*r8;
   double r10 = sqrt(sqr(r9) + sqr(r6) + sqr(r3));
   double r11 = 1.0/r10;
   double r12 = 0.5*r11*r3;
   double r13 = 0.5*r11*r6;
   double r14 = 0.5*r11*r9;
   double r15 = 0.5*r11*(0.5*r2 - 0.5*r1);
   double r16 = 0.5*r11*(0.5*r5 - 0.5*r4);
   double r17 = 0.5*r11*(0.5*r8 - 0.5*r7);
   V += Factor*(r10);
   dVdX[iAx] += Factor*(r12);
   dVdX[iAy] += Factor*(r13);
   dVdX[iAz] += Factor*(r14);
   dVdX[iBx] += Factor*(r12);
   dVdX[iBy] += Factor*(r13);
   dVdX[iBz] += Factor*(r14);
   dVdX[iCx] += Factor*(r15);
   dVdX[iCy] += Factor*(r16);
   dVdX[iCz] += Factor*(r17);
   dVdX[iDx] += Factor*(r15);
   dVdX[iDy] += Factor*(r16);
   dVdX[iDz] += Factor*(r17);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz; (void)Cx; (void)Cy;
   (void)Cz; (void)Dx; (void)Dy; (void)Dz;
}

/* kate: syntax c++; */

