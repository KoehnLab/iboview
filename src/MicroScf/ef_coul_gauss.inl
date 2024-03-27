// Coulomb interaction energy between two s-type Gaussian charges
void ef_coul_gauss(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double q1, double q2, double sigma1, double sigma2)
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
   double r4 = sqr(r3) + sqr(r2) + sqr(r1);
   double r5 = sqrt(r4);
   double r6 = 1.414213562373095;
   double r7 = sqr(sigma2) + sqr(sigma1);
   double r8 = 1.0/sqrt(r7);
   double r9 = erf((1.0*r5*r8)/r6);
   double r10 = 0.5641895835477563;
   double r11 = 1.0/r4;
   double r12 = 1.0/pow(2.718281828459045,(r4/(2.0*r7)));
   double r13 = 1.0/pow(r5,3);
   V += Factor*((1.0*r9*q1*q2)/r5);
   dVdX[iAx] += Factor*(r6*r10*r1*r11*q1*q2*r8*r12 - 1.0*r1*r13*r9*q1*q2);
   dVdX[iAy] += Factor*(r6*r10*r2*r11*q1*q2*r8*r12 - 1.0*r13*r2*r9*q1*q2);
   dVdX[iAz] += Factor*(r6*r10*r11*r3*q1*q2*r8*r12 - 1.0*r13*r3*r9*q1*q2);
   dVdX[iBx] += Factor*(r1*r13*q1*q2*r9 - 1.0*r1*r10*r11*r12*r6*r8*q1*q2);
   dVdX[iBy] += Factor*(r2*r13*q1*q2*r9 - 1.0*r10*r11*r12*r2*r6*r8*q1*q2);
   dVdX[iBz] += Factor*(r13*r3*q1*q2*r9 - 1.0*r10*r11*r12*r3*r6*r8*q1*q2);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz; (void)q1; (void)q2;
   (void)sigma1; (void)sigma2;
}

/* kate: syntax c++; */

