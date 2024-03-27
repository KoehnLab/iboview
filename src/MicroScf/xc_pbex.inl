// xc functional: PBE GGA exchange functional. (http://dx.doi.org/10.1103/PhysRevLett.77.3865)
void xc_pbex(double &E, double *dE, double const *D, double const Factor)
{
   double rhoc = D[0];
   double rhoo = D[1];
   double sigmacc = D[2];
   double sigmaco = D[3];
   double sigmaoo = D[4];
   double r1 = 1.442249570307408;
   double r2 = 1.464591887561523;
   double r3 = 1.0/r2;
   double r4 = rhoc - 1.0*rhoo;
   double r5 = pow(r4,(4.0/3.0));
   double r6 = 0.39685026299205;
   double r7 = 0.30285343213869;
   double r8 = 2.145029397111025;
   double r9 = sigmaoo - 2.0*sigmaco + sigmacc;
   double r10 = (0.02766376451077943*r6*r7*r8*r9)/pow(r4,(8.0/3.0)) + 1;
   double r11 = 1.804 - 0.804/r10;
   double r12 = rhoo + rhoc;
   double r13 = pow(r12,(4.0/3.0));
   double r14 = sigmaoo + 2.0*sigmaco + sigmacc;
   double r15 = (0.02766376451077943*r14*r6*r7*r8)/pow(r12,(8.0/3.0)) + 1;
   double r16 = 1.804 - 0.804/r15;
   double r17 = 3.174802103936399;
   double r18 = 1.0/pow(r4,(7.0/3.0));
   double r19 = 1.0/sqr(r10);
   double r20 = 0.4807498567691361;
   double r21 = pow(r4,(1.0/3.0));
   double r22 = 1.0/sqr(r15);
   double r23 = (0.002780208333333333*r1*r14*r17*r2*r22*r7)/pow(r12,(7.0/3.0));
   double r24 = -1.5*r20*r3*pow(r12,(1.0/3.0))*r16;
   double r25 = 1.0/r5;
   double r26 = 1.0/r13;
   double r27 = (-0.008340625*r6*r1*r7*r2*r26*r22) - 0.008340625*r6*r1*r7*r2*r25*r19;
   double r28 = 0.7937005259840997;
   E += Factor*((-0.375*r1*r3*r13*r16) - 0.375*r1*r3*r5*r11);
   dE[0] += Factor*(r24 + r23 - 1.5*r20*r3*r21*r11 + 0.002780208333333333*r17*r1*r7*r2*r18*r9*r19);
   dE[1] += Factor*(r24 + r23 + 1.5*r20*r3*r21*r11 - 0.002780208333333333*r17*r1*r7*r2*r18*r9*r19);
   dE[2] += Factor*(r27);
   dE[3] += Factor*(0.008340625*r28*r1*r7*r2*r25*r19 - 0.008340625*r28*r1*r7*r2*r26*r22);
   dE[4] += Factor*(r27);
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo; (void)sigmacc; (void)sigmaco; (void)sigmaoo;
}

/* kate: syntax c++; */

