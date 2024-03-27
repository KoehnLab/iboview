// atom cartesian coordinate of atom A
void pc_pos_z(void *, double &V, double *dVdX, double const &Factor, double const *pXyz, index_t const *iAtom)
{
   index_t iAx = 3.0*iAtom[0] + 0;
   index_t iAy = 3.0*iAtom[0] + 1;
   index_t iAz = 3.0*iAtom[0] + 2;
   double Ax = pXyz[iAx];
   double Ay = pXyz[iAy];
   double Az = pXyz[iAz];
   V += Factor*(Az);
   dVdX[iAx] += Factor*(0);
   dVdX[iAy] += Factor*(0);
   dVdX[iAz] += Factor*(1);
   // suppress potential unused variable warnings.
   (void)Ax; (void)Ay; (void)Az;
}

/* kate: syntax c++; */

