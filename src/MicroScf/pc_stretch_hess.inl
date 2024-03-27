// bond stretch coordinate: Distance between A and B
void pc_stretch_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom)
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
   double r7 = sqrt(r6 + r4 + r2);
   double r8 = 1.0/pow(r7,3);
   double r9 = 1.0/r7;
   double r10 = r9 - 1.0*r2*r8;
   double r11 = -1.0*r1*r3*r8;
   double r12 = -1.0*r1*r5*r8;
   double r13 = -1.0*r9;
   double r14 = r13 + r2*r8;
   double r15 = r1*r3*r8;
   double r16 = r1*r8*r5;
   double r17 = r9 - 1.0*r4*r8;
   double r18 = -1.0*r3*r5*r8;
   double r19 = r13 + r4*r8;
   double r20 = r3*r8*r5;
   double r21 = r9 - 1.0*r6*r8;
   double r22 = r8*r6 + r13;
   dVdXX[iAx + nb*(iAx)] += Factor*(r10);
   dVdXX[iAx + nb*(iAy)] += Factor*(r11);
   dVdXX[iAx + nb*(iAz)] += Factor*(r12);
   dVdXX[iAx + nb*(iBx)] += Factor*(r14);
   dVdXX[iAx + nb*(iBy)] += Factor*(r15);
   dVdXX[iAx + nb*(iBz)] += Factor*(r16);
   dVdXX[iAy + nb*(iAx)] += Factor*(r11);
   dVdXX[iAy + nb*(iAy)] += Factor*(r17);
   dVdXX[iAy + nb*(iAz)] += Factor*(r18);
   dVdXX[iAy + nb*(iBx)] += Factor*(r15);
   dVdXX[iAy + nb*(iBy)] += Factor*(r19);
   dVdXX[iAy + nb*(iBz)] += Factor*(r20);
   dVdXX[iAz + nb*(iAx)] += Factor*(r12);
   dVdXX[iAz + nb*(iAy)] += Factor*(r18);
   dVdXX[iAz + nb*(iAz)] += Factor*(r21);
   dVdXX[iAz + nb*(iBx)] += Factor*(r16);
   dVdXX[iAz + nb*(iBy)] += Factor*(r20);
   dVdXX[iAz + nb*(iBz)] += Factor*(r22);
   dVdXX[iBx + nb*(iAx)] += Factor*(r14);
   dVdXX[iBx + nb*(iAy)] += Factor*(r15);
   dVdXX[iBx + nb*(iAz)] += Factor*(r16);
   dVdXX[iBx + nb*(iBx)] += Factor*(r10);
   dVdXX[iBx + nb*(iBy)] += Factor*(r11);
   dVdXX[iBx + nb*(iBz)] += Factor*(r12);
   dVdXX[iBy + nb*(iAx)] += Factor*(r15);
   dVdXX[iBy + nb*(iAy)] += Factor*(r19);
   dVdXX[iBy + nb*(iAz)] += Factor*(r20);
   dVdXX[iBy + nb*(iBx)] += Factor*(r11);
   dVdXX[iBy + nb*(iBy)] += Factor*(r17);
   dVdXX[iBy + nb*(iBz)] += Factor*(r18);
   dVdXX[iBz + nb*(iAx)] += Factor*(r16);
   dVdXX[iBz + nb*(iAy)] += Factor*(r20);
   dVdXX[iBz + nb*(iAz)] += Factor*(r22);
   dVdXX[iBz + nb*(iBx)] += Factor*(r12);
   dVdXX[iBz + nb*(iBy)] += Factor*(r18);
   dVdXX[iBz + nb*(iBz)] += Factor*(r21);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az; (void)Bx; (void)By; (void)Bz;
}

/* kate: syntax c++; */

