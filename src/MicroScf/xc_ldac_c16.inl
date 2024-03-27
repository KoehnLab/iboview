// xc functional: C16: original Chachiyo C16 LSDA correlation functional (J. Chem. Phys. 145, 021101 (2016); https://doi.org/10.1063/1.4958669) without the updated parameters from Karasiev (J. Chem. Phys. 145, 157101 (2016); https://doi.org/10.1063/1.4964758)**/
// /* #C2F-OutName: xc_ldac_c16  */
// /* #C2F-IncrementOutput: True */
// /* #C2F-Inputs: rhoc,rhoo */
// 
// /* note also its improvement described in:
// Comment on “Communication: Simple and accurate uniform electron gas correlation energy for the full range of densities” [J. Chem. Phys. 145, 021101 (2016)]
// J. Chem. Phys. 145, 157101 (2016); https://doi.org/10.1063/1.4964758
// Valentin V. Karasieva
// */
// 
// 
// 
// /* seitz radius */
// SeitzRs(rhoc) := (3/(4*%pi*rhoc))^(1/3);
// Rs(rhoc) := SeitzRs(rhoc);
// 
// /* relative spin polarization (-1...+1), Zeta = (nAlpha-nBeta)/(nAlpha+nBeta) */
// SpinPolarizationZeta(rhoo,rhoc) := rhoo/rhoc;
// Zeta(rhoo,rhoc) := SpinPolarizationZeta(rhoo,rhoc);
// 
// 
// 
// /* Cf: prefactor of Dirac exchange energy (should have minus sign) */
// /* parameters: ['Cf = -.75*(3/%pi)^(1/3)]; */
// Cf_DiracX : -.75*(3/%pi)^(1/3);
// 
// /* unpolarized free electron gas exchange energy (per volume, not per particle!)*/
// E_DiracX_u(rhoc) := Cf_DiracX * rhoc^(4/3);
// 
// /* polarized energy from scaling:  (http://dx.doi.org/10.1103/PhysRevA.20.397, Eq. 2.13, 2.14)
// * Ex(rhoa,rhob) = 1/2 (Ex(2 rhoa) + Ex(2 rhob))
// * Ex(rhoc,Zeta) = (1/2) ((1 + Zeta)**(4/3) + (1 - Zeta)^(4/3)) * Ex(rhoc)
// */
// 
// /* note: wikipedia says this is the exact exchange energy scaling for general exchange functionals, not only LDAX */
// FSpinX(Zeta) := ((1 + Zeta)^(4/3) + (1 - Zeta)^(4/3))/2;
// E_DiracX_p(rhoc,Zeta) := FSpinX(Zeta) * E_DiracX_u(rhoc);
// 
// 
// 
// 
// 
// /* kate: syntax maxima; */
// 
// /* C16, eq. (3) (values for paramagnetic case/zero polarization),
// * and eq. (12) (ferromagnetic/full polarization)*/
// /*parameters: ['C16_a_zerop = log(2)/(2*%pi^2), 'C16_b_zerop = 20.4562557,
// 'C16_a_fullp = (log(2)-1)/(4*%pi^2), 'C16_b_fullp = 27.4203609 ];*/
// C16_a_zerop: (log(2)-1)/(2*%pi^2);
// C16_a_fullp: (log(2)-1)/(4*%pi^2);
// C16_b_zerop: 20.4562557;
// C16_b_fullp: 27.4203609;
// /* ^- these are the original values */
// /* C16_b_zerop: 21.7392245; */
// /* C16_b_fullp: 28.3559732; */
// /* ^- these are the updated values from Karasiev */
// 
// 
// /* C16, eq. (1) (note: technically we could use more constants in the logarithm, and fit them explicitly.
// There is no theoretical need for the rs^-1 and rs^-2 coefficients to be identical (see discussion in paper
// around eq. (9)). It is just chosen for simplicity by assuming that rs=1 is the dividing mark between
// low-density behavior and high-density behavior (based on some discussion of the hydrogen atom)).
// 
// See also discussion around eq. (10): The unpolarized form is taken from analytically computed quantities
// using perturbation theory around the paramagnetic limit. The ferromagnetic limit has different a/b.
// */
// C16_LDAC_unpolarized_rs(rs) := C16_a_zerop * log(1 + C16_b_zerop/rs + C16_b_zerop/(rs^2));
// C16_LDAC_unpolarized(rhoc) := C16_LDAC_unpolarized_rs(SeitzRs(rhoc));
// 
// C16_LDAC_fully_polarized_rs(rs) := C16_a_fullp * log(1 + C16_b_fullp/rs + C16_b_fullp/(rs^2));
// 
// /* paragraph between eq. (11) and (12) */
// C16_Fzeta(Zeta) := ((1 + Zeta)^(4/3) + (1 - Zeta)^(4/3) - 2)/(2*(2^(1/3) - 1));
// 
// 
// C16_LDAC_polarized_rs(rs,Zeta) := C16_LDAC_unpolarized_rs(rs) + C16_Fzeta(Zeta) * (C16_LDAC_fully_polarized_rs(rs) - C16_LDAC_unpolarized_rs(rs));
// C16_LDAC_polarized(rhoc,rhoo) := C16_LDAC_polarized_rs(SeitzRs(rhoc),rhoo/rhoc);
// 
// 
// /* note: The function described in the paper is the energy
// *per electron*, not energy per volume. We thus still multiply it by the density.*/
// E : rhoc * C16_LDAC_polarized(rhoc,rhoo);
// /* E : rhoc * C16_LDAC_unpolarized(rhoc); */
// 
// 
// /* kate: syntax maxima; */
// 
// dE_drhoc : diff(E,rhoc);
// dE_drhoo : diff(E,rhoo);
// fname: 'xc_ldac_c16;
// outputs: ['E = E,
// 'dE_drhoc = dE_drhoc,
// 'dE_drhoo = dE_drhoo];
// inputs: [rhoc,
// rhoo]
void xc_ldac_c16(double &E, double *dE, double const *D, double const Factor)
{
   double rhoc = D[0];
   double rhoo = D[1];
   double r1 = 0.1013211836423378;
   double r2 = -0.3068528194400547;
   double r3 = 0.6933612743506348;
   double r4 = 1.587401051968199;
   double r5 = 1.464591887561523;
   double r6 = pow(rhoc,(1.0/3.0));
   double r7 = 0.4807498567691361;
   double r8 = 2.519842099789746;
   double r9 = 2.145029397111025;
   double r10 = pow(rhoc,(2.0/3.0));
   double r11 = 20.4562557*r3*r4*r5*r6 + 20.4562557*r7*r8*r9*r10 + 1;
   double r12 = log(r11);
   double r13 = 0.5*r1*r12*r2;
   double r14 = 3.847322101863072;
   double r15 = 27.4203609*r3*r4*r5*r6 + 27.4203609*r7*r8*r9*r10 + 1;
   double r16 = 0.25*r1*log(r15)*r2 - 0.5*r1*r12*r2;
   double r17 = 1.0/rhoc;
   double r18 = 1 - 1.0*r17*rhoo;
   double r19 = r17*rhoo + 1;
   double r20 = pow(r19,(4.0/3.0)) + pow(r18,(4.0/3.0)) - 2;
   double r21 = 0.5*r14*r16*r20;
   double r22 = 0.2311204247835449;
   double r23 = 1.0/r10;
   double r24 = 0.1602499522563787;
   double r25 = 1.0/r6;
   double r26 = 40.9125114*r24*r8*r9*r25 + 20.4562557*r22*r4*r5*r23;
   double r27 = 1.0/r11;
   double r28 = 1.0/sqr(rhoc);
   double r29 = pow(r18,(1.0/3.0));
   double r30 = pow(r19,(1.0/3.0));
   E += Factor*(rhoc*(r21 + r13));
   dE[0] += Factor*(rhoc*(0.5*r14*r16*(1.333333333333333*r28*r29*rhoo - 1.333333333333333*r28*r30*rhoo) + 0.5*r14*r20*((0.25*r1*r2*(54.8407218*r24*r8*r9*r25 + 27.4203609*r22*r4*r5*r23))/r15 - 0.5*r1*r2*r26*r27) + 0.5*r1*r2*r26*r27) + r21 + r13);
   dE[1] += Factor*(0.5*r14*r16*(1.333333333333333*r17*r30 - 1.333333333333333*r17*r29)*rhoc);
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo;
}

/* kate: syntax c++; */

