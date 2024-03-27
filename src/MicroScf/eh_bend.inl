// bond bend harmonic energy
void eh_bend(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double theta0, double kijk)
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
   double Ax = pXyz[iAx];
   double Ay = pXyz[iAy];
   double Az = pXyz[iAz];
   double Bx = pXyz[iBx];
   double By = pXyz[iBy];
   double Bz = pXyz[iBz];
   double Cx = pXyz[iCx];
   double Cy = pXyz[iCy];
   double Cz = pXyz[iCz];
   double r1 = -1.0*Bx;
   double r2 = r1 + Ax;
   double r3 = -1.0*By;
   double r4 = r3 + Ay;
   double r5 = -1.0*Bz;
   double r6 = r5 + Az;
   double r7 = sqr(r6) + sqr(r4) + sqr(r2);
   double r8 = sqrt(r7);
   double r9 = 1.0/r8;
   double r10 = Cx + r1;
   double r11 = Cy + r3;
   double r12 = Cz + r5;
   double r13 = r6*r12 + r4*r11 + r2*r10;
   double r14 = sqr(r12) + sqr(r11) + sqr(r10);
   double r15 = sqrt(r14);
   double r16 = 1.0/r15;
   double r17 = acos(r9*r13*r16) - 1.0*theta0;
   double r18 = 1.0/sqrt(1.0 - (1.0*sqr(r13))/(r14*r7));
   double r19 = 1.0/pow(r8,3);
   double r20 = 1.0/pow(r15,3);
   V += Factor*(0.5*sqr(r17)*kijk);
   dVdX[iAx] += Factor*(-1.0*r17*r18*(r9*r10*r16 - 1.0*r13*r16*r19*r2)*kijk);
   dVdX[iAy] += Factor*(-1.0*r17*r18*(r9*r11*r16 - 1.0*r13*r16*r19*r4)*kijk);
   dVdX[iAz] += Factor*(-1.0*r17*r18*(r9*r12*r16 - 1.0*r13*r16*r19*r6)*kijk);
   dVdX[iBx] += Factor*(-1.0*r17*r18*(r16*r9*((-1.0*Cx) + 2.0*Bx - 1.0*Ax) + r9*r10*r13*r20 + r2*r19*r13*r16)*kijk);
   dVdX[iBy] += Factor*(-1.0*r17*r18*(r16*r9*((-1.0*Cy) + 2.0*By - 1.0*Ay) + r9*r11*r13*r20 + r4*r19*r13*r16)*kijk);
   dVdX[iBz] += Factor*(-1.0*r17*r18*(r16*r9*((-1.0*Cz) + 2.0*Bz - 1.0*Az) + r9*r12*r13*r20 + r19*r6*r13*r16)*kijk);
   dVdX[iCx] += Factor*(-1.0*r17*r18*(r2*r9*r16 - 1.0*r10*r13*r20*r9)*kijk);
   dVdX[iCy] += Factor*(-1.0*r17*r18*(r4*r9*r16 - 1.0*r11*r13*r20*r9)*kijk);
   dVdX[iCz] += Factor*(-1.0*r17*r18*(r9*r6*r16 - 1.0*r12*r13*r20*r9)*kijk);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz; (void)Cx; (void)Cy;
   (void)Cz; (void)theta0; (void)kijk;
}

/* kate: syntax c++; */

