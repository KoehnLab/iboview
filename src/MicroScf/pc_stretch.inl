// bond stretch coordinate: Distance between A and B
void pc_stretch(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom)
{
   index_t iAx = 3.0*iAtom[0] + 0;
   index_t iAy = 3.0*iAtom[0] + 1;
   index_t iAz = 3.0*iAtom[0] + 2;
   index_t iBx = 3.0*iAtom[1] + 0;
   index_t iBy = 3.0*iAtom[1] + 1;
   index_t iBz = 3.0*iAtom[1] + 2;
   double Ax = pXyz[iAx];
   double Ay = pXyz[iAy];
   double Az = pXyz[iAz];
   double Bx = pXyz[iBx];
   double By = pXyz[iBy];
   double Bz = pXyz[iBz];
   double r1 = Ax - 1.0*Bx;
   double r2 = Ay - 1.0*By;
   double r3 = Az - 1.0*Bz;
   double r4 = sqrt(sqr(r3) + sqr(r2) + sqr(r1));
   double r5 = 1.0/r4;
   V += Factor*(r4);
   dVdX[iAx] += Factor*(r1*r5);
   dVdX[iAy] += Factor*(r2*r5);
   dVdX[iAz] += Factor*(r5*r3);
   dVdX[iBx] += Factor*(-1.0*r1*r5);
   dVdX[iBy] += Factor*(-1.0*r2*r5);
   dVdX[iBz] += Factor*(-1.0*r3*r5);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz;
}

/* kate: syntax c++; */

