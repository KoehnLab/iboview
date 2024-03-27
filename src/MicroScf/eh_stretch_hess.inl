// bond stretch harmonic energy
void eh_stretch_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom, double re, double k)
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
   double r2 = sqr(r1);
   double r3 = Ay - 1.0*By;
   double r4 = sqr(r3);
   double r5 = Az - 1.0*Bz;
   double r6 = sqr(r5);
   double r7 = 1.0/(r6 + r4 + r2);
   double r8 = r2*r7*k;
   double r9 = r1*r3*r7*k;
   double r10 = r1*r7*r5*k;
   double r11 = -1.0*r2*r7*k;
   double r12 = -1.0*r1*r3*r7*k;
   double r13 = -1.0*r1*r5*r7*k;
   double r14 = r4*r7*k;
   double r15 = r3*r7*r5*k;
   double r16 = -1.0*r4*r7*k;
   double r17 = -1.0*r3*r5*r7*k;
   double r18 = r7*r6*k;
   double r19 = -1.0*r6*r7*k;
   dVdXX[iAx + nb*(iAx)] += Factor*(r8);
   dVdXX[iAx + nb*(iAy)] += Factor*(r9);
   dVdXX[iAx + nb*(iAz)] += Factor*(r10);
   dVdXX[iAx + nb*(iBx)] += Factor*(r11);
   dVdXX[iAx + nb*(iBy)] += Factor*(r12);
   dVdXX[iAx + nb*(iBz)] += Factor*(r13);
   dVdXX[iAy + nb*(iAx)] += Factor*(r9);
   dVdXX[iAy + nb*(iAy)] += Factor*(r14);
   dVdXX[iAy + nb*(iAz)] += Factor*(r15);
   dVdXX[iAy + nb*(iBx)] += Factor*(r12);
   dVdXX[iAy + nb*(iBy)] += Factor*(r16);
   dVdXX[iAy + nb*(iBz)] += Factor*(r17);
   dVdXX[iAz + nb*(iAx)] += Factor*(r10);
   dVdXX[iAz + nb*(iAy)] += Factor*(r15);
   dVdXX[iAz + nb*(iAz)] += Factor*(r18);
   dVdXX[iAz + nb*(iBx)] += Factor*(r13);
   dVdXX[iAz + nb*(iBy)] += Factor*(r17);
   dVdXX[iAz + nb*(iBz)] += Factor*(r19);
   dVdXX[iBx + nb*(iAx)] += Factor*(r11);
   dVdXX[iBx + nb*(iAy)] += Factor*(r12);
   dVdXX[iBx + nb*(iAz)] += Factor*(r13);
   dVdXX[iBx + nb*(iBx)] += Factor*(r8);
   dVdXX[iBx + nb*(iBy)] += Factor*(r9);
   dVdXX[iBx + nb*(iBz)] += Factor*(r10);
   dVdXX[iBy + nb*(iAx)] += Factor*(r12);
   dVdXX[iBy + nb*(iAy)] += Factor*(r16);
   dVdXX[iBy + nb*(iAz)] += Factor*(r17);
   dVdXX[iBy + nb*(iBx)] += Factor*(r9);
   dVdXX[iBy + nb*(iBy)] += Factor*(r14);
   dVdXX[iBy + nb*(iBz)] += Factor*(r15);
   dVdXX[iBz + nb*(iAx)] += Factor*(r13);
   dVdXX[iBz + nb*(iAy)] += Factor*(r17);
   dVdXX[iBz + nb*(iAz)] += Factor*(r19);
   dVdXX[iBz + nb*(iBx)] += Factor*(r10);
   dVdXX[iBz + nb*(iBy)] += Factor*(r15);
   dVdXX[iBz + nb*(iBz)] += Factor*(r18);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz; (void)re; (void)k;
}

/* kate: syntax c++; */

