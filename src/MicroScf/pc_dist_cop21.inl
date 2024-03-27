// distance between center of positions (AB) and C
void pc_dist_cop21(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom)
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
   double r1 = 0.5*(Bx + Ax) - 1.0*Cx;
   double r2 = 0.5*(By + Ay) - 1.0*Cy;
   double r3 = 0.5*(Bz + Az) - 1.0*Cz;
   double r4 = sqrt(sqr(r3) + sqr(r2) + sqr(r1));
   double r5 = 1.0/r4;
   double r6 = 0.5*r1*r5;
   double r7 = 0.5*r2*r5;
   double r8 = 0.5*r3*r5;
   V += Factor*(r4);
   dVdX[iAx] += Factor*(r6);
   dVdX[iAy] += Factor*(r7);
   dVdX[iAz] += Factor*(r8);
   dVdX[iBx] += Factor*(r6);
   dVdX[iBy] += Factor*(r7);
   dVdX[iBz] += Factor*(r8);
   dVdX[iCx] += Factor*(-1.0*r1*r5);
   dVdX[iCy] += Factor*(-1.0*r2*r5);
   dVdX[iCz] += Factor*(-1.0*r3*r5);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz; (void)Cx; (void)Cy;
   (void)Cz;
}

/* kate: syntax c++; */

