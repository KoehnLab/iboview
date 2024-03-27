// xc functional: TPSS correlation functional
// Functional Form: Pure functional; Semi-local meta-GGA (tau-form)
// Parameterization: Electron gas + empirical functional form compatible with some exact scaling limits and avoiding 1e self-interaction errors in the correlation part.
// Original article: http://dx.doi.org/10.1103/PhysRevLett.91.146401
// 
// Implemented by Alyssa Bienvenu, June 2016
void xc_tpssc(double &E, double *dE, double const *D, double const Factor)
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
   double sigmacc = D[2];
   double sigmaco = D[3];
   double sigmaoo = D[4];
   double tauc = D[5];
   double tauo = D[6];
   double r1 = 1.0/pow(rhoc,3);
   double r2 = pow(sigmacc,3);
   double r3 = 1.0/sqr(rhoc);
   double r4 = sqr(sigmacc);
   double r5 = sqr(rhoo);
   double r6 = 1.0/pow(rhoc,4);
   double r7 = pow(rhoo,4);
   double r8 = 1.0/pow(rhoc,6);
   double r9 = pow(rhoo,6);
   double r10 = 2.26*r8*r9 + 0.5*r6*r7 + 0.87*r3*r5 + 0.53;
   double r11 = 2.080083823051904;
   double r12 = 1.0/r11;
   double r13 = 0.2173369174628993;
   double r14 = 1.0/pow(rhoc,(8.0/3.0));
   double r15 = 1.0/rhoc;
   double r16 = 1 - 1.0*r15*rhoo;
   double r17 = pow(r16,(4.0/3.0));
   double r18 = r15*rhoo + 1;
   double r19 = pow(r18,(4.0/3.0));
   double r20 = 1.0/r19 + 1.0/r17;
   double r21 = sigmaoo - 2.0*r15*rhoo*sigmaco + r3*r5*sigmacc;
   double r22 = 0.125*r12*r13*r14*r20*r21 + 1;
   double r23 = 1.0/pow(r22,4);
   double r24 = r10*r23 + 1;
   double r25 = 0.5*rhoc;
   double r26 = r25 - 0.5*rhoo;
   double r27 = 1.0/Ap;
   double r28 = 1.587401051968199;
   double r29 = 1.0/r28;
   double r30 = 2.145029397111025;
   double r31 = 1.0/r30;
   double r32 = rhoc - 1.0*rhoo;
   double r33 = 1.414213562373095;
   double r34 = 1.0/r33;
   double r35 = 1.732050807568877;
   double r36 = 0.5641895835477563;
   double r37 = sqrt(r32);
   double r38 = 1.259921049894873;
   double r39 = 1.0/r38;
   double r40 = 1.442249570307408;
   double r41 = 1.464591887561523;
   double r42 = 1.0/r41;
   double r43 = 1.0/pow(r32,(1.0/3.0));
   double r44 = 0.8908987181403393;
   double r45 = 1.200936955176003;
   double r46 = 0.8263074871107581;
   double r47 = (1.0*r11*r29*r31*Beta4p)/pow(r32,(2.0/3.0)) + (1.0*r34*r35*r36*Beta3p)/r37 + (1.0*r44*r45*r46*Beta1p)/pow(r32,(1.0/6.0)) + r39*r40*r42*Beta2p*r43;
   double r48 = (0.5*r27)/r47 + 1;
   double r49 = log(r48);
   double r50 = r39*r40*r42*Alpha1p*r43 + 1;
   double r51 = 1.0/r40;
   double r52 = 1.0/pow(r32,(7.0/3.0));
   double r53 = sigmaoo - 2.0*sigmaco + sigmacc;
   double r54 = pow(2.718281828459045,(128.6558737716593*Ap*r49*r50));
   double r55 = r54 - 1;
   double r56 = 1.0/r55;
   double r57 = 0.2682657924959205*r51*r41*r52*r56*r53;
   double r58 = r57 + 1;
   double r59 = 1.0/pow(r32,(14.0/3.0));
   double r60 = 1.0/sqr(r55);
   double r61 = sqr(r53);
   double r62 = 0.07196653542346429*r12*r30*r59*r60*r61 + r57 + 1;
   double r63 = 1.0/r62;
   double r64 = 0.2682657924959205*r51*r41*r52*r53*r58*r63 + 1;
   double r65 = 0.01554534543482745*log(r64) - 2.0*r49*r50*Ap;
   double r66 = 0.5*rhoo + r25;
   double r67 = rhoo + rhoc;
   double r68 = 1.0/pow(r67,(1.0/3.0));
   double r69 = r39*r40*r42*Alpha1p*r68 + 1;
   double r70 = sqrt(r67);
   double r71 = (1.0*r11*r29*r31*Beta4p)/pow(r67,(2.0/3.0)) + (1.0*r34*r35*r36*Beta3p)/r70 + (1.0*r44*r45*r46*Beta1p)/pow(r67,(1.0/6.0)) + r39*r40*r42*Beta2p*r68;
   double r72 = (0.5*r27)/r71 + 1;
   double r73 = log(r72);
   double r74 = 1.0/pow(r67,(7.0/3.0));
   double r75 = sigmaoo + 2.0*sigmaco + sigmacc;
   double r76 = pow(2.718281828459045,(128.6558737716593*Ap*r69*r73));
   double r77 = r76 - 1;
   double r78 = 1.0/r77;
   double r79 = 0.2682657924959205*r51*r41*r74*r78*r75;
   double r80 = r79 + 1;
   double r81 = 1.0/pow(r67,(14.0/3.0));
   double r82 = 1.0/sqr(r77);
   double r83 = sqr(r75);
   double r84 = 0.07196653542346429*r12*r30*r81*r82*r83 + r79 + 1;
   double r85 = 1.0/r84;
   double r86 = 0.2682657924959205*r51*r41*r74*r75*r80*r85 + 1;
   double r87 = 0.01554534543482745*log(r86) - 2.0*r69*r73*Ap;
   double r88 = r15*r66*r87 + r15*r26*r65;
   double r89 = 1.0/sqr(tauc);
   double r90 = 0.39685026299205;
   double r91 = 1.0/pow(rhoc,(2.0/3.0));
   double r92 = sqrt(rhoc);
   double r93 = 1.0/r92;
   double r94 = 0.6299605249474366;
   double r95 = 1.0/pow(rhoc,(1.0/3.0));
   double r96 = 1.0/pow(rhoc,(1.0/6.0));
   double r97 = 0.5*r35*r36*r93*Beta3u + r39*r45*r46*Beta1u*r96 + r40*r94*r42*Beta2u*r95 + r11*r90*r31*Beta4u*r91;
   double r98 = 0.5/(r97*Au) + 1;
   double r99 = log(r98);
   double r100 = r40*r94*r42*Alpha1u*r95 + 1;
   double r101 = 2.519842099789746;
   double r102 = 1.0/(r101 - 2.0);
   double r103 = 0.5*r35*r36*r93*Beta3p + r39*r45*r46*Beta1p*r96 + r40*r94*r42*Beta2p*r95 + r11*r90*r31*Beta4p*r91;
   double r104 = (0.5*r27)/r103 + 1;
   double r105 = log(r104);
   double r106 = r40*r94*r42*Alpha1p*r95 + 1;
   double r107 = 2.0*r100*r99*Au;
   double r108 = r107 - 2.0*r105*r106*Ap;
   double r109 = r19 + r17 - 2;
   double r110 = 0.5*r35*r36*r93*Beta3a + r39*r45*r46*Beta1a*r96 + r40*r94*r42*Beta2a*r95 + r11*r90*r31*Beta4a*r91;
   double r111 = 0.5/(r110*Aa) + 1;
   double r112 = log(r111);
   double r113 = r40*r94*r42*Alpha1a*r95 + 1;
   double r114 = 1 - 1.0*r6*r7;
   double r115 = pow(r18,(2.0/3.0)) + pow(r16,(2.0/3.0));
   double r116 = pow(r115,3);
   double r117 = 1.0/pow(rhoc,(7.0/3.0));
   double r118 = sqr(r115);
   double r119 = 1.0/r118;
   double r120 = 1.0/r116;
   double r121 = (-2.25*r109*r112*r113*r114*Aa) - 1.0*r102*r108*r109*r6*r7 + r107;
   double r122 = pow(2.718281828459045,(257.3117475433186*r120*r121));
   double r123 = r122 - 1;
   double r124 = 1.0/r123;
   double r125 = 0.5365315849918411*r51*r41*r117*r119*r124*sigmacc;
   double r126 = r125 + 1;
   double r127 = 1.0/pow(rhoc,(14.0/3.0));
   double r128 = 1.0/pow(r115,4);
   double r129 = 1.0/sqr(r123);
   double r130 = 0.2878661416938572*r12*r30*r127*r128*r129*r4 + r125 + 1;
   double r131 = 1.0/r130;
   double r132 = 0.5365315849918411*r51*r41*r117*r119*sigmacc*r126*r131 + 1;
   double r133 = log(r132);
   double r134 = (-2.0*r100*r99*Au) + 2.25*r109*r112*r113*r114*Aa + 0.003886336358706862*r116*r133 + r102*r108*r6*r7*r109;
   double r135 = 0.015625*r10*r23*r3*r4*r89 + 1;
   double r136 = r134*r135 - 0.015625*r24*r3*r4*r88*r89;
   double r137 = 1.0/pow(tauc,3);
   double r138 = 0.00546875*r1*r2*r136*r137 + 1;
   double r139 = 1.0/r48;
   double r140 = 1.0/pow(r32,(5.0/3.0));
   double r141 = 1.0/pow(r33,3);
   double r142 = 1.0/pow(r37,3);
   double r143 = 1.0/pow(r32,(4.0/3.0));
   double r144 = 0.4454493590701696;
   double r145 = 0.4003123183920009;
   double r146 = 1.0/pow(r32,(7.0/6.0));
   double r147 = (-1.0*r140*r31*r38*r51*Beta4p) - 1.0*r141*r142*r35*r36*Beta3p - 1.0*r12*r143*r39*r42*Beta2p - 1.0*r144*r145*r146*r46*Beta1p;
   double r148 = 1.0/sqr(r47);
   double r149 = 1.0/r64;
   double r150 = 1.0/sqr(r62);
   double r151 = 4.326748710922225;
   double r152 = 1.0/r151;
   double r153 = 1.0/pow(r32,(10.0/3.0));
   double r154 = -1.877860547471444*r152*r41*r153*r56*r53;
   double r155 = (-128.6558737716593*r39*r12*r42*Alpha1p*Ap*r49*r143) - 64.32793688582964*r139*r147*r50*r148;
   double r156 = -0.2682657924959205*r51*r41*r155*r52*r60*r54*r53;
   double r157 = 0.1602499522563787;
   double r158 = 1.0/pow(r32,(17.0/3.0));
   double r159 = 1.0/pow(r55,3);
   double r160 = 1.0/pow(r67,(4.0/3.0));
   double r161 = (-(1.0*r31*r38*r51*Beta4p)/pow(r67,(5.0/3.0))) - (1.0*r141*r35*r36*Beta3p)/pow(r70,3) - 1.0*r12*r160*r39*r42*Beta2p - (1.0*r144*r145*r46*Beta1p)/pow(r67,(7.0/6.0));
   double r162 = 1.0/sqr(r71);
   double r163 = 1.0/r72;
   double r164 = 1.0/r86;
   double r165 = 1.0/sqr(r84);
   double r166 = 1.0/pow(r67,(10.0/3.0));
   double r167 = -1.877860547471444*r152*r41*r166*r78*r75;
   double r168 = (-128.6558737716593*r39*r12*r42*Alpha1p*Ap*r160*r73) - 64.32793688582964*r161*r69*r162*r163;
   double r169 = -0.2682657924959205*r51*r41*r74*r82*r76*r168*r75;
   double r170 = r15*r66*(0.01554534543482745*r164*((-1.877860547471444*r152*r41*r166*r75*r80*r85) + 0.2682657924959205*r51*r41*r74*r75*(r169 + r167)*r85 - 0.2682657924959205*r165*r41*r51*r74*r75*r80*((-(1.0075314959285*r157*r30*r82*r83)/pow(r67,(17.0/3.0))) - (0.1439330708469286*r12*r168*r30*r76*r81*r83)/pow(r77,3) + r169 + r167)) + r28*r12*r42*Alpha1p*Ap*r160*r73 + r161*r69*r162*r163);
   double r171 = 0.5*r15*r87;
   double r172 = 1.0/pow(r16,(7.0/3.0));
   double r173 = 1.0/pow(r18,(7.0/3.0));
   double r174 = 0.125*r12*r13*r14*r20*(2.0*r3*rhoo*sigmaco - 2.0*r1*r5*sigmacc) + 0.125*r12*r13*r14*r21*(1.333333333333333*r173*r3*rhoo - 1.333333333333333*r172*r3*rhoo) - (1.0*r13*r157*r20*r21)/pow(rhoc,(11.0/3.0));
   double r175 = 1.0/pow(r22,5);
   double r176 = 1.0/pow(rhoc,5);
   double r177 = (-(13.56*r9)/pow(rhoc,7)) - 2.0*r176*r7 - 1.74*r1*r5;
   double r178 = 1.0/r98;
   double r179 = 1.0/pow(rhoc,(5.0/3.0));
   double r180 = 1.0/pow(r92,3);
   double r181 = 1.0/pow(rhoc,(4.0/3.0));
   double r182 = 1.0/r101;
   double r183 = 1.0/pow(rhoc,(7.0/6.0));
   double r184 = (-2.0*r179*r31*r51*r90*Beta4u) - 0.25*r180*r35*r36*Beta3u - 1.0*r12*r181*r42*r94*Beta2u - 1.0*r145*r182*r183*r46*Beta1u;
   double r185 = 1.0/sqr(r97);
   double r186 = pow(r16,(1.0/3.0));
   double r187 = pow(r18,(1.0/3.0));
   double r188 = 1.333333333333333*r186*r3*rhoo - 1.333333333333333*r187*r3*rhoo;
   double r189 = -1.0*r100*r178*r184*r185;
   double r190 = -2.0*r12*r181*r42*r94*r99*Alpha1u*Au;
   double r191 = (1.0*r106*((-2.0*r179*r31*r51*r90*Beta4p) - 0.25*r180*r35*r36*Beta3p - 1.0*r12*r181*r42*r94*Beta2p - 1.0*r145*r182*r183*r46*Beta1p))/(sqr(r103)*r104) + 2.0*r105*r12*r181*r42*r94*Alpha1p*Ap + r190 + r189;
   double r192 = 1.0/r111;
   double r193 = (-2.0*r179*r31*r51*r90*Beta4a) - 0.25*r180*r35*r36*Beta3a - 1.0*r12*r181*r42*r94*Beta2a - 1.0*r145*r182*r183*r46*Beta1a;
   double r194 = 1.0/sqr(r110);
   double r195 = 0.1574901312368591;
   double r196 = 1.0/r132;
   double r197 = 1.0/sqr(r130);
   double r198 = 1.0/r186;
   double r199 = 1.0/r187;
   double r200 = 0.6666666666666666*r198*r3*rhoo - 0.6666666666666666*r199*r3*rhoo;
   double r201 = -1.073063169983682*r51*r41*r117*r200*r120*r124*sigmacc;
   double r202 = 1.0/pow(rhoc,(10.0/3.0));
   double r203 = -3.755721094942888*r152*r41*r202*r119*r124*sigmacc;
   double r204 = 257.3117475433186*r120*((-9.0*r109*r112*r113*r176*r7*Aa) - 2.25*r112*r113*r114*r188*Aa - 1.0*r102*r109*r191*r6*r7 - 1.0*r102*r108*r188*r6*r7 + 4.0*r102*r108*r109*r176*r7 + 1.125*r109*r113*r114*r192*r193*r194 + r190 + r189 + r151*r195*r42*Aa*Alpha1a*r112*r181*r114*r109) - 771.9352426299557*r200*r128*r121;
   double r205 = -0.5365315849918411*r51*r41*r117*r119*r204*r129*r122*sigmacc;
   double r206 = 1.0/pow(r115,5);
   double r207 = 1.0/pow(rhoc,(17.0/3.0));
   double r208 = 1.0/pow(r123,3);
   double r209 = r135*(2.0*r12*r181*r42*r94*r99*Alpha1u*Au - 1.0*r109*r112*r114*r151*r181*r195*r42*Aa*Alpha1a + 9.0*r109*r112*r113*r176*r7*Aa + 2.25*r112*r113*r114*r188*Aa - 4.0*r102*r108*r109*r176*r7 - 1.125*r109*r113*r114*r192*r193*r194 + r102*r108*r6*r7*r188 + r178*r184*r100*r185 + 0.01165900907612059*r200*r118*r133 + 0.003886336358706862*r116*r196*(0.5365315849918411*r51*r41*r117*r119*sigmacc*(r205 + r203 + r201)*r131 - 3.755721094942888*r152*r41*r202*r119*sigmacc*r126*r131 - 1.073063169983682*r51*r41*r117*r200*r120*sigmacc*r126*r131 - 0.5365315849918411*r51*r41*r117*r119*sigmacc*r126*r197*((-0.5757322833877143*r12*r30*r127*r128*r204*r208*r122*r4) - 4.030125983714*r157*r30*r207*r128*r129*r4 - 1.151464566775429*r12*r30*r127*r200*r206*r129*r4 + r205 + r203 + r201)) + r102*r191*r6*r7*r109) + r134*(0.015625*r177*r23*r3*r4*r89 - 0.0625*r10*r174*r175*r3*r4*r89 - 0.03125*r1*r10*r23*r4*r89) - 0.015625*(r177*r23 - 4.0*r10*r174*r175)*r3*r4*r88*r89 + 0.03125*r1*r24*r4*r88*r89 - 0.015625*r24*r3*r4*((-1.0*r3*r66*r87) - 1.0*r26*r3*r65 + 0.5*r15*r65 + r15*r26*(0.01554534543482745*r149*(0.2682657924959205*r51*r41*r52*r53*(r156 + r154)*r63 - 1.877860547471444*r152*r41*r153*r53*r58*r63 - 0.2682657924959205*r51*r41*r52*r53*r58*r150*((-0.1439330708469286*r12*r30*r155*r59*r159*r54*r61) - 1.0075314959285*r157*r30*r158*r60*r61 + r156 + r154)) + r28*r12*r42*Alpha1p*Ap*r49*r143 + r139*r147*r50*r148) + r171 + r170)*r89;
   double r210 = r144*r145*r46*Beta1p*r146 + r39*r12*r42*Beta2p*r143 + r141*r35*r36*Beta3p*r142 + r38*r51*r31*Beta4p*r140;
   double r211 = 1.877860547471444*r152*r41*r153*r56*r53;
   double r212 = 128.6558737716593*r39*r12*r42*Alpha1p*Ap*r49*r143 - 64.32793688582964*r139*r210*r50*r148;
   double r213 = -0.2682657924959205*r51*r41*r212*r52*r60*r54*r53;
   double r214 = 0.125*r12*r13*r14*r20*(2.0*r3*rhoo*sigmacc - 2.0*r15*sigmaco) + 0.125*r12*r13*r14*(1.333333333333333*r15*r172 - 1.333333333333333*r15*r173)*r21;
   double r215 = pow(rhoo,3);
   double r216 = 13.56*r8*pow(rhoo,5) + 2.0*r6*r215 + 1.74*r3*rhoo;
   double r217 = 1.333333333333333*r15*r187 - 1.333333333333333*r15*r186;
   double r218 = 0.6666666666666666*r15*r199 - 0.6666666666666666*r15*r198;
   double r219 = -1.073063169983682*r51*r41*r117*r218*r120*r124*sigmacc;
   double r220 = 257.3117475433186*r120*(9.0*r109*r112*r113*r215*r6*Aa - 2.25*r112*r113*r114*r217*Aa - 1.0*r102*r108*r217*r6*r7 - 4.0*r102*r108*r109*r215*r6) - 771.9352426299557*r218*r128*r121;
   double r221 = -0.5365315849918411*r51*r41*r117*r119*r220*r129*r122*sigmacc;
   double r222 = (-0.015625*r24*r3*r4*r89*(r15*r26*((-1.0*r12*r143*r28*r42*r49*Alpha1p*Ap) + 0.01554534543482745*r149*(0.2682657924959205*r51*r41*r52*r53*(r213 + r211)*r63 + 1.877860547471444*r152*r41*r153*r53*r58*r63 - 0.2682657924959205*r51*r41*r52*r53*r58*r150*((-0.1439330708469286*r12*r30*r212*r59*r159*r54*r61) + 1.0075314959285*r157*r30*r158*r60*r61 + r213 + r211)) + r139*r210*r50*r148) - 0.5*r15*r65 + r171 + r170)) + r135*((-9.0*r109*r112*r113*r215*r6*Aa) + 2.25*r112*r113*r114*r217*Aa + 4.0*r102*r108*r109*r215*r6 + r102*r108*r6*r7*r217 + 0.01165900907612059*r218*r118*r133 + 0.003886336358706862*r116*r196*(0.5365315849918411*r51*r41*r117*r119*sigmacc*(r221 + r219)*r131 - 1.073063169983682*r51*r41*r117*r218*r120*sigmacc*r126*r131 - 0.5365315849918411*r51*r41*r117*r119*sigmacc*r126*r197*((-0.5757322833877143*r12*r30*r127*r128*r220*r208*r122*r4) - 1.151464566775429*r12*r30*r127*r218*r206*r129*r4 + r221 + r219))) + r134*(0.015625*r216*r23*r3*r4*r89 - 0.0625*r10*r175*r214*r3*r4*r89) - 0.015625*(r216*r23 - 4.0*r10*r175*r214)*r3*r4*r88*r89;
   double r223 = -0.015625*(0.01554534543482745*r15*r66*(0.2682657924959205*r51*r41*r74*r80*r85 + 0.07196653542346429*r12*r30*r81*r78*r75*r85 - 0.2682657924959205*r51*r41*r74*r75*(0.1439330708469286*r12*r30*r81*r82*r75 + 0.2682657924959205*r51*r41*r74*r78)*r80*r165)*r164 + 0.01554534543482745*r15*r26*(0.2682657924959205*r51*r41*r52*r58*r63 + 0.07196653542346429*r12*r30*r59*r56*r53*r63 - 0.2682657924959205*r51*r41*r52*r53*(0.1439330708469286*r12*r30*r59*r60*r53 + 0.2682657924959205*r51*r41*r52*r56)*r58*r150)*r149)*r24*r3*r4*r89;
   double r224 = 1.0/pow(rhoc,(20.0/3.0));
   double r225 = r134*(0.03125*r10*r23*r3*r89*sigmacc - 0.0078125*r10*r12*r13*r175*r20*r224*r4*r5*r89) - 0.03125*r24*r3*r88*r89*sigmacc + 0.0078125*r10*r12*r13*r175*r20*r224*r4*r5*r88*r89 + r223 + 0.003886336358706862*r116*(0.5365315849918411*r51*r41*r117*r119*r126*r131 + 0.2878661416938572*r12*r30*r127*r128*r124*sigmacc*r131 - 0.5365315849918411*r51*r41*r117*r119*sigmacc*(0.5757322833877143*r12*r30*r127*r128*r129*sigmacc + 0.5365315849918411*r51*r41*r117*r119*r124)*r126*r197)*r196*r135;
   double r226 = (-0.015625*r10*r12*r13*r175*r20*r207*r4*r88*r89*rhoo) + 0.015625*r10*r12*r13*r134*r175*r20*r207*r4*r89*rhoo - 0.015625*(0.01554534543482745*r15*r66*(0.5365315849918411*r51*r41*r74*r80*r85 + 0.1439330708469286*r12*r30*r81*r78*r75*r85 - 0.2682657924959205*r51*r41*r74*r75*(0.2878661416938572*r12*r30*r81*r82*r75 + 0.5365315849918411*r51*r41*r74*r78)*r80*r165)*r164 + 0.01554534543482745*r15*r26*((-0.5365315849918411*r51*r41*r52*r58*r63) - 0.1439330708469286*r12*r30*r59*r56*r53*r63 - 0.2682657924959205*r51*r41*r52*r53*((-0.2878661416938572*r12*r30*r59*r60*r53) - 0.5365315849918411*r51*r41*r52*r56)*r58*r150)*r149)*r24*r3*r4*r89;
   double r227 = 0.0078125*r10*r12*r127*r13*r175*r20*r4*r88*r89 - 0.0078125*r10*r12*r127*r13*r134*r175*r20*r4*r89 + r223;
   double r228 = 0.03125*r137*r24*r3*r4*r88 - 0.03125*r10*r134*r137*r23*r3*r4;
   E += Factor*(rhoc*r138*r136);
   dE[0] += Factor*(r138*r136 + rhoc*(0.00546875*r1*r2*r209*r137 - 0.01640625*r6*r2*r136*r137)*r136 + rhoc*r138*r209);
   dE[1] += Factor*(rhoc*r138*r222 + 0.00546875*r3*r2*r222*r136*r137);
   dE[2] += Factor*(rhoc*(0.01640625*r1*r4*r136*r137 + 0.00546875*r1*r2*r225*r137)*r136 + rhoc*r138*r225);
   dE[3] += Factor*(rhoc*r138*r226 + 0.00546875*r3*r2*r226*r136*r137);
   dE[4] += Factor*(rhoc*r138*r227 + 0.00546875*r3*r2*r227*r136*r137);
   dE[5] += Factor*(r136*rhoc*(0.00546875*r1*r2*r228*r137 - (0.01640625*r1*r136*r2)/pow(tauc,4)) + rhoc*r228*r138);
   dE[6] += Factor*(0);
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo; (void)sigmacc; (void)sigmaco; (void)sigmaoo; (void)tauc;
   (void)tauo; (void)Au; (void)Alpha1u; (void)Beta1u; (void)Beta2u; (void)Beta3u;
   (void)Beta4u; (void)Ap; (void)Alpha1p; (void)Beta1p; (void)Beta2p; (void)Beta3p;
   (void)Beta4p; (void)Aa; (void)Alpha1a; (void)Beta1a; (void)Beta2a; (void)Beta3a;
   (void)Beta4a;
}

/* kate: syntax c++; */
