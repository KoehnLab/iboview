// atom cartesian coordinate of atom A
void pc_pos_x_hess(void *, double *dVdXX, index_t nb, double const &Factor, double const *pXyz, index_t const *iAtom)
{
   index_t iAx = 3.0*iAtom[0] + 0;
   index_t iAy = 3.0*iAtom[0] + 1;
   index_t iAz = 3.0*iAtom[0] + 2;
   double Ax = pXyz[iAx];
   double Ay = pXyz[iAy];
   double Az = pXyz[iAz];
   dVdXX[iAx + nb*(iAx)] += Factor*(0);
   dVdXX[iAx + nb*(iAy)] += Factor*(0);
   dVdXX[iAx + nb*(iAz)] += Factor*(0);
   dVdXX[iAy + nb*(iAx)] += Factor*(0);
   dVdXX[iAy + nb*(iAy)] += Factor*(0);
   dVdXX[iAy + nb*(iAz)] += Factor*(0);
   dVdXX[iAz + nb*(iAx)] += Factor*(0);
   dVdXX[iAz + nb*(iAy)] += Factor*(0);
   dVdXX[iAz + nb*(iAz)] += Factor*(0);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az;
}

/* kate: syntax c++; */

