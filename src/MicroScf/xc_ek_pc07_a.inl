// ek functional: PC07 Laplace-Level kinetic energy functional (alpha-spin component tau_alpha)
// Functional Form: Pure functional; Semi-local meta-GGA (upsilon-form)
// Parameterization: Electron gas + empirical functional form compatible with some exact scaling limits.
// Original article: Perdew, Constantin; Phys. Rev. B 75 155109 (2006), http://dx.doi.org/10.1103/PhysRevB.75.155109
// 
// Implemented by Alyssa Bienvenu, June 2016
void xc_ek_pc07_a(double &E, double *dE, double const *D, double const Factor)
{
   double rhoc = D[0];
   double rhoo = D[1];
   double sigmacc = D[2];
   double sigmaco = D[3];
   double sigmaoo = D[4];
   double tauc = D[5];
   double tauo = D[6];
   double upsilonc = D[7];
   double upsilono = D[8];
   double r1 = 6.240251469155714;
   double r2 = 4.601151114470489;
   double r3 = rhoo + rhoc;
   double r4 = pow(r3,(5.0/3.0));
   double r5 = 1.0/r1;
   double r6 = 1.0/r2;
   double r7 = 1.0/pow(r3,(8.0/3.0));
   double r8 = sigmaoo + 2.0*sigmaco + sigmacc;
   double r9 = 1.25*r5*r6*r7*r8;
   double r10 = 0.01780555025070874;
   double r11 = 0.0009511128591915428;
   double r12 = 0.04723533569227512;
   double r13 = 1.0/pow(r3,(16.0/3.0));
   double r14 = sqr(r8);
   double r15 = 0.5*r11*r12*r13*r14;
   double r16 = 0.05341665075212624;
   double r17 = 1.0/r4;
   double r18 = upsilono + upsilonc;
   double r19 = 0.02568004719817164;
   double r20 = 1.0/pow(r3,(13.0/3.0));
   double r21 = -0.0625*r12*r18*r19*r20*r8;
   double r22 = 0.002853338577574629;
   double r23 = 1.0/pow(r3,(10.0/3.0));
   double r24 = sqr(r18);
   double r25 = 0.5*r12*r22*r23*r24;
   double r26 = 1.25*r10*r6*r7*r8 + 5.0*r16*r17*r18*r6 + r25 + r21 + r15 + 1;
   double r27 = r9 + 1;
   double r28 = 1.0/sqr(r27);
   double r29 = r25 + r21 + r15;
   double r30 = sqr(r29);
   double r31 = sqrt(r28*r30 + 1);
   double r32 = 1.0/r31;
   double r33 = r26*r32 - 1.25*r5*r6*r7*r8;
   double r34 = PC07_Fab(r33);
   double r35 = r33*r34 + r9;
   double r36 = 1.0/pow(r3,(11.0/3.0));
   double r37 = -(0.002536300957844113*r12*r14)/pow(r3,(19.0/3.0));
   double r38 = 0.006955012782838153*r12*r13*r18*r8;
   double r39 = -5.0*r11*r12*r20*r24;
   double r40 = 1.0/pow(r27,3);
   double r41 = 1.0/pow(r31,3);
   double r42 = (-0.5*r26*r41*(20.0*r16*r30*r36*r40*r6*r8 + 2.0*r28*r29*(r39 + r38 + r37))) + r32*((-0.05935183416902912*r36*r6*r8) - 25.0*r10*r18*r6*r7 + r39 + r38 + r37) + 10.0*r16*r36*r6*r8;
   double r43 = dPC07_Fab(r33);
   double r44 = 0.05*r1*r2*r4*((-10.0*r16*r36*r6*r8) + r42*r33*r43 + r42*r34) + 0.520020955762976*r2*pow(r3,(2.0/3.0))*r35;
   double r45 = r11*r12*r13*r8;
   double r46 = -0.0625*r12*r18*r19*r20;
   double r47 = (-0.5*r26*r41*(2.0*r28*r29*(r46 + r45) - 2.5*r30*r40*r5*r6*r7)) + r32*(1.25*r10*r6*r7 + r46 + r45) - 1.25*r5*r6*r7;
   double r48 = 0.05*r1*r2*r4*(1.25*r5*r6*r7 + r47*r33*r43 + r47*r34);
   double r49 = 2.0*r11*r12*r13*r8;
   double r50 = -0.125*r12*r18*r19*r20;
   double r51 = (-0.5*r26*r41*(2.0*r28*r29*(r50 + r49) - 5.0*r30*r40*r5*r6*r7)) + r32*(2.5*r10*r6*r7 + r50 + r49) - 2.5*r5*r6*r7;
   double r52 = -0.0625*r12*r19*r20*r8;
   double r53 = r22*r12*r23*r18;
   double r54 = r32*(5.0*r16*r17*r6 + r53 + r52) - 1.0*r26*r28*r29*r41*(r53 + r52);
   double r55 = 0.05*r1*r2*r4*(r54*r33*r43 + r54*r34);
   E += Factor*(0.05*r1*r2*r35*r4);
   dE[0] += Factor*(r44);
   dE[1] += Factor*(r44);
   dE[2] += Factor*(r48);
   dE[3] += Factor*(0.05*r1*r2*r4*(2.5*r5*r6*r7 + r51*r33*r43 + r51*r34));
   dE[4] += Factor*(r48);
   dE[5] += Factor*(0);
   dE[6] += Factor*(0);
   dE[7] += Factor*(r55);
   dE[8] += Factor*(r55);
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo; (void)sigmacc; (void)sigmaco; (void)sigmaoo; (void)tauc;
   (void)tauo; (void)upsilonc; (void)upsilono;
}

/* kate: syntax c++; */

