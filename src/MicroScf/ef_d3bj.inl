// approximate hessian for D3(BJ) dispersion correction (approx: c6 held fixed)
void ef_d3bj(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom, double c6, double r2r4_term, double a1, double a2, double s6, double s18)
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
   double r5 = 1.732050807568877*a1*sqrt(r2r4_term) + a2;
   double r6 = pow(r5,8) + pow(r4,4);
   double r7 = pow(r4,3);
   double r8 = pow(r5,6) + r7;
   double r9 = 1.0/sqr(r6);
   double r10 = sqr(r4);
   double r11 = 1.0/sqr(r8);
   V += Factor*((1.0*c6*s6)/r8 + (3.0*c6*r2r4_term*s18)/r6);
   dVdX[iAx] += Factor*((-6.0*r1*r10*r11*c6*s6) - 24.0*r1*r7*r9*c6*r2r4_term*s18);
   dVdX[iAy] += Factor*((-6.0*r10*r11*r2*c6*s6) - 24.0*r2*r7*r9*c6*r2r4_term*s18);
   dVdX[iAz] += Factor*((-6.0*r10*r11*r3*c6*s6) - 24.0*r3*r7*r9*c6*r2r4_term*s18);
   dVdX[iBx] += Factor*(6.0*r1*r10*r11*c6*s6 + 24.0*r1*r7*r9*c6*r2r4_term*s18);
   dVdX[iBy] += Factor*(6.0*r10*r11*r2*c6*s6 + 24.0*r2*r7*r9*c6*r2r4_term*s18);
   dVdX[iBz] += Factor*(6.0*r10*r11*r3*c6*s6 + 24.0*r3*r7*r9*c6*r2r4_term*s18);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz; (void)c6;
   (void)r2r4_term; (void)a1; (void)a2; (void)s6; (void)s18;
}

/* kate: syntax c++; */

