// xc functional: PW92c LSDA correlation functional. (http://dx.doi.org/10.1103/PhysRevB.45.13244)
void xc_pw92c(double &E, double *dE, double const *D, double const Factor)
{
   double const Au = 0.0310907;
   double const Alpha1u = 0.2137;
   double const Beta1u = 7.5957;
   double const Beta2u = 3.5876;
   double const Beta3u = 1.6382;
   double const Beta4u = 0.49294;
   double const Ap = 0.01554535;
   double const Alpha1p = 0.20548;
   double const Beta1p = 14.1189;
   double const Beta2p = 6.1977;
   double const Beta3p = 3.3662;
   double const Beta4p = 0.62517;
   double const Aa = 0.0168869;
   double const Alpha1a = 0.11125;
   double const Beta1a = 10.357;
   double const Beta2a = 3.6231;
   double const Beta3a = 0.88026;
   double const Beta4a = 0.49671;
   double rhoc = D[0];
   double rhoo = D[1];
   double r1 = 2.080083823051904;
   double r2 = 0.39685026299205;
   double r3 = 0.4661940770354117;
   double r4 = 1.0/pow(rhoc,(2.0/3.0));
   double r5 = 1.732050807568877;
   double r6 = 0.5641895835477563;
   double r7 = sqrt(rhoc);
   double r8 = 1.0/r7;
   double r9 = 1.442249570307408;
   double r10 = 0.6299605249474366;
   double r11 = 0.6827840632552957;
   double r12 = 1.0/pow(rhoc,(1.0/3.0));
   double r13 = 0.7937005259840997;
   double r14 = 1.200936955176003;
   double r15 = 0.8263074871107581;
   double r16 = 1.0/pow(rhoc,(1.0/6.0));
   double r17 = 0.5*r5*r6*r8*Beta3u + r1*r2*r3*Beta4u*r4 + r13*r14*r15*Beta1u*r16 + r9*r10*r11*Beta2u*r12;
   double r18 = 0.5/(r17*Au) + 1;
   double r19 = log(r18);
   double r20 = r9*r10*r11*Alpha1u*r12 + 1;
   double r21 = -2.0*r19*r20*Au;
   double r22 = 2.519842099789746;
   double r23 = 1.0/(r22 - 2.0);
   double r24 = 0.5*r5*r6*r8*Beta3p + r1*r2*r3*Beta4p*r4 + r13*r14*r15*Beta1p*r16 + r9*r10*r11*Beta2p*r12;
   double r25 = 0.5/(r24*Ap) + 1;
   double r26 = log(r25);
   double r27 = r9*r10*r11*Alpha1p*r12 + 1;
   double r28 = 2.0*r19*r20*Au - 2.0*r26*r27*Ap;
   double r29 = 1.0/pow(rhoc,4);
   double r30 = pow(rhoo,4);
   double r31 = 1.0/rhoc;
   double r32 = 1 - 1.0*r31*rhoo;
   double r33 = r31*rhoo + 1;
   double r34 = pow(r33,(4.0/3.0)) + pow(r32,(4.0/3.0)) - 2;
   double r35 = r23*r28*r29*r30*r34;
   double r36 = 0.5*r5*r6*r8*Beta3a + r1*r2*r3*Beta4a*r4 + r13*r14*r15*Beta1a*r16 + r9*r10*r11*Beta2a*r12;
   double r37 = 0.5/(r36*Aa) + 1;
   double r38 = log(r37);
   double r39 = r9*r10*r11*Alpha1a*r12 + 1;
   double r40 = 1 - 1.0*r29*r30;
   double r41 = 2.25*r34*r38*r39*r40*Aa;
   double r42 = 1.0/r18;
   double r43 = 1.0/r9;
   double r44 = 1.0/pow(rhoc,(5.0/3.0));
   double r45 = 1.0/pow(r7,3);
   double r46 = 1.0/r1;
   double r47 = 1.0/pow(rhoc,(4.0/3.0));
   double r48 = 1.0/r22;
   double r49 = 0.4003123183920009;
   double r50 = 1.0/pow(rhoc,(7.0/6.0));
   double r51 = (-2.0*r2*r3*r43*r44*Beta4u) - 0.25*r45*r5*r6*Beta3u - 1.0*r10*r11*r46*r47*Beta2u - 1.0*r15*r48*r49*r50*Beta1u;
   double r52 = 1.0/sqr(r17);
   double r53 = 1.0/sqr(rhoc);
   double r54 = pow(r32,(1.0/3.0));
   double r55 = pow(r33,(1.0/3.0));
   double r56 = 1.333333333333333*r53*r54*rhoo - 1.333333333333333*r53*r55*rhoo;
   double r57 = 1.0/pow(rhoc,5);
   double r58 = 1.333333333333333*r31*r55 - 1.333333333333333*r31*r54;
   double r59 = pow(rhoo,3);
   E += Factor*(rhoc*(r41 + r35 + r21));
   dE[0] += Factor*((r23*r29*r30*r34*((1.0*r27*((-2.0*r2*r3*r43*r44*Beta4p) - 0.25*r45*r5*r6*Beta3p - 1.0*r10*r11*r46*r47*Beta2p - 1.0*r15*r48*r49*r50*Beta1p))/(sqr(r24)*r25) - 2.0*r10*r11*r19*r46*r47*Alpha1u*Au + 2.0*r10*r11*r26*r46*r47*Alpha1p*Ap - 1.0*r20*r42*r51*r52) - (1.125*r34*r39*r40*((-2.0*r2*r3*r43*r44*Beta4a) - 0.25*r45*r5*r6*Beta3a - 1.0*r10*r11*r46*r47*Beta2a - 1.0*r15*r48*r49*r50*Beta1a))/(sqr(r36)*r37) + 2.0*r10*r11*r19*r46*r47*Alpha1u*Au - 0.6814202223120525*r11*r34*r38*r40*r47*Aa*Alpha1a + 9.0*r30*r34*r38*r39*r57*Aa + 2.25*r38*r39*r40*r56*Aa - 4.0*r23*r28*r30*r34*r57 + r23*r28*r29*r30*r56 + r42*r51*r20*r52)*rhoc + r41 + r35 + r21);
   dE[1] += Factor*(((-9.0*r29*r34*r38*r39*r59*Aa) + 2.25*r38*r39*r40*r58*Aa + 4.0*r23*r28*r29*r34*r59 + r23*r28*r29*r30*r58)*rhoc);
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo; (void)Au; (void)Alpha1u; (void)Beta1u; (void)Beta2u;
   (void)Beta3u; (void)Beta4u; (void)Ap; (void)Alpha1p; (void)Beta1p; (void)Beta2p;
   (void)Beta3p; (void)Beta4p; (void)Aa; (void)Alpha1a; (void)Beta1a; (void)Beta2a;
   (void)Beta3a; (void)Beta4a;
}

/* kate: syntax c++; */

