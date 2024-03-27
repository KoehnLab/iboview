// xc functional: Dirac-Slater LSDA exchange functional (http://dx.doi.org/10.1103/PhysRev.81.385).**/
// /* #C2F-OutName: xc_diracx  */
// /* #C2F-IncrementOutput: True */
// /* #C2F-Inputs: rhoc,rhoo */
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
// 
// /* E : (1/2)*E_DiracX_u(2*rhoc); */
// /* E : E_DiracX_u(rhoc); */
// 
// E : E_DiracX_p(rhoc, Zeta(rhoo,rhoc));
// 
// /* E : (1/2)*(E_DiracX_u(rhoc+rhoo)+E_DiracX_u(rhoc-rhoo)); */
// 
// /* kate: syntax maxima; */
// dE_drhoc : diff(E,rhoc);
// dE_drhoo : diff(E,rhoo);
// fname: 'xc_diracx;
// outputs: ['E = E,
// 'dE_drhoc = dE_drhoc,
// 'dE_drhoo = dE_drhoo];
// inputs: [rhoc,
// rhoo]
void xc_diracx(double &E, double *dE, double const *D, double const Factor)
{
   double rhoc = D[0];
   double rhoo = D[1];
   double r1 = 1.442249570307408;
   double r2 = 0.6827840632552957;
   double r3 = pow(rhoc,(4.0/3.0));
   double r4 = 1.0/rhoc;
   double r5 = 1 - 1.0*r4*rhoo;
   double r6 = r4*rhoo + 1;
   double r7 = pow(r6,(4.0/3.0)) + pow(r5,(4.0/3.0));
   double r8 = 1.0/sqr(rhoc);
   double r9 = pow(r5,(1.0/3.0));
   double r10 = pow(r6,(1.0/3.0));
   E += Factor*(-0.375*r1*r2*r3*r7);
   dE[0] += Factor*((-0.375*r1*r2*r3*(1.333333333333333*r8*r9*rhoo - 1.333333333333333*r10*r8*rhoo)) - 0.7211247851537042*r2*r7*pow(rhoc,(1.0/3.0)));
   dE[1] += Factor*(-0.375*r1*r2*r3*(1.333333333333333*r10*r4 - 1.333333333333333*r4*r9));
   // suppress potential unused variable warnings.
   (void)rhoc; (void)rhoo;
}

/* kate: syntax c++; */

