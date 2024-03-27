// ek functional: PC07 Laplace-Level kinetic energy functional (beta-spin component tau_beta)
// Functional Form: Pure functional; Semi-local meta-GGA (upsilon-form)
// Parameterization: Electron gas + empirical functional form compatible with some exact scaling limits.
// Original article: Perdew, Constantin; Phys. Rev. B 75 155109 (2006), http://dx.doi.org/10.1103/PhysRevB.75.155109
// 
// Implemented by Alyssa Bienvenu, June 2016
void xc_ek_pc07_b(double &E, double *dE, double const *D, double const Factor)
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
   double r3 = rhoc - 1.0*rhoo;
   double r4 = pow(r3,(5.0/3.0));
   double r5 = 1.0/r1;
   double r6 = 1.0/r2;
   double r7 = 1.0/pow(r3,(8.0/3.0));
   double r8 = sigmaoo - 2.0*sigmaco + sigmacc;
   double r9 = 1.25*r5*r6*r7*r8;
   double r10 = r9 + 1;
   double r11 = 1.0/sqr(r10);
   double r12 = 0.0009511128591915428;
   double r13 = 0.04723533569227512;
   double r14 = 1.0/pow(r3,(16.0/3.0));
   double r15 = sqr(r8);
   double r16 = 0.5*r12*r13*r14*r15;
   double r17 = 0.002853338577574629;
   double r18 = 1.0/pow(r3,(10.0/3.0));
   double r19 = upsilonc - 1.0*upsilono;
   double r20 = sqr(r19);
   double r21 = 0.5*r13*r17*r18*r20;
   double r22 = 0.02568004719817164;
   double r23 = 1.0/pow(r3,(13.0/3.0));
   double r24 = -0.0625*r13*r19*r22*r23*r8;
   double r25 = r24 + r21 + r16;
   double r26 = sqr(r25);
   double r27 = sqrt(r11*r26 + 1);
   double r28 = 1.0/r27;
   double r29 = 0.01780555025070874;
   double r30 = 0.05341665075212624;
   double r31 = 1.0/r4;
   double r32 = 1.25*r29*r6*r7*r8 + 5.0*r19*r30*r31*r6 + r24 + r21 + r16 + 1;
   double r33 = r28*r32 - 1.25*r5*r6*r7*r8;
   double r34 = PC07_Fab(r33);
   double r35 = r34*r33 + r9;
   double r36 = 2.080083823051904;
   double r37 = pow(r3,(2.0/3.0));
   double r38 = 1.0/pow(r3,(11.0/3.0));
   double r39 = -10.0*r30*r38*r6*r8;
   double r40 = 10.0*r30*r38*r6*r8;
   double r41 = 0.005935183416902912;
   double r42 = 0.0003170376197305141;
   double r43 = 1.0/pow(r3,(19.0/3.0));
   double r44 = -8.0*r13*r15*r42*r43;
   double r45 = -5.0*r12*r13*r20*r23;
   double r46 = 0.00856001573272388;
   double r47 = 0.8125*r13*r14*r19*r46*r8;
   double r48 = 1.0/pow(r27,3);
   double r49 = 1.0/pow(r10,3);
   double r50 = (-0.5*r32*r48*(20.0*r26*r30*r38*r49*r6*r8 + 2.0*r11*r25*(r47 + r45 + r44))) + r28*((-10.0*r38*r41*r6*r8) - 25.0*r19*r29*r6*r7 + r47 + r45 + r44) + r40;
   double r51 = dPC07_Fab(r33);
   double r52 = 8.0*r13*r15*r42*r43;
   double r53 = 5.0*r12*r13*r20*r23;
   double r54 = -0.8125*r13*r14*r19*r46*r8;
   double r55 = (-0.5*r32*r48*(2.0*r11*r25*(r54 + r53 + r52) - 20.0*r26*r30*r38*r49*r6*r8)) + r28*(10.0*r38*r41*r6*r8 + 25.0*r19*r29*r6*r7 + r54 + r53 + r52) + r39;
   double r56 = r12*r13*r14*r8;
   double r57 = -0.0625*r13*r19*r22*r23;
   double r58 = (-0.5*r32*r48*(2.0*r11*r25*(r57 + r56) - 2.5*r26*r49*r5*r6*r7)) + r28*(1.25*r29*r6*r7 + r57 + r56) - 1.25*r5*r6*r7;
   double r59 = 0.05*r1*r2*r4*(1.25*r5*r6*r7 + r34*r58 + r51*r33*r58);
   double r60 = -2.0*r12*r13*r14*r8;
   double r61 = 0.125*r13*r19*r22*r23;
   double r62 = (-0.5*r32*r48*(5.0*r26*r49*r5*r6*r7 + 2.0*r11*r25*(r61 + r60))) + r28*((-2.5*r29*r6*r7) + r61 + r60) + 2.5*r5*r6*r7;
   double r63 = -0.0625*r13*r22*r23*r8;
   double r64 = r17*r13*r18*r19;
   double r65 = r28*(r64 + r63 + 5.0*r30*r31*r6) - 1.0*r11*r25*r32*r48*(r64 + r63);
   double r66 = 0.0625*r13*r22*r23*r8;
   double r67 = -1.0*r13*r17*r18*r19;
   double r68 = r28*(r67 + r66 - 5.0*r30*r31*r6) - 1.0*r11*r25*r32*r48*(r67 + r66);
   E += Factor*(0.05*r1*r2*r35*r4);
   dE[0] += Factor*(0.05*r1*r2*r4*(r51*r33*r50 + r34*r50 + r39) + 0.25*r2*r35*r36*r37);
   dE[1] += Factor*(0.05*r1*r2*r4*(r51*r33*r55 + r34*r55 + r40) - 0.25*r2*r35*r36*r37);
   dE[2] += Factor*(r59);
   dE[3] += Factor*(0.05*r1*r2*r4*((-2.5*r5*r6*r7) + r34*r62 + r51*r33*r62));
   dE[4] += Factor*(r59);
   dE[5] += Factor*(0);
   dE[6] += Factor*(0);
   dE[7] += Factor*(0.05*r1*r2*r4*(r51*r33*r65 + r34*r65));
   dE[8] += Factor*(0.05*r1*r2*r4*(r51*r33*r68 + r34*r68));
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo; (void)sigmacc; (void)sigmaco; (void)sigmaoo; (void)tauc;
   (void)tauo; (void)upsilonc; (void)upsilono;
}

/* kate: syntax c++; */

