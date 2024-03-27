/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#include <cmath>
#include <stdint.h>
#include <algorithm>  // for std::max.

#include "Ir.h"
#include "IrAmrr.h"
#include "CtInt1e.h"
#include "CtBasisShell.h"

namespace ct {


// return x^n for low n.
template<class FScalar>
inline FScalar intpow( FScalar x, int n ){
   FScalar r = 1;
   for ( ; n > 0; -- n )
      r *= x;
   return r;
}

// template<class FScalar>
// inline FScalar sqr(FScalar x) {
//    return x*x;
// }


// ####################################################################
// ********************************************************************
//
//      Now starting: integral code for 1e^- integrals (i.e. overlap,
//      multipoles, kinetic terms, nuclear attraction and similiar stuff)
//
// --------------------------------------------------------------------

FDoublettIntegralFactory::~FDoublettIntegralFactory()
{}


// most integrals have just one component.
uint FDoublettIntegralFactory::nComp() const
{
   return 1;
}


// returns: RealTarget[ j*(AngMomA+1) + i ] = overlap for cartesian power i
// in shell a and cartesian power j in shell b. RealTarget must provide size
// for (AngMomA+1) * (AngMomB+1) scalars.
void EvalOverlapIntegrals1D( FScalar *RealTarget, uint AngMomA, uint AngMomB,
   FScalar ZetaSum, FScalar ZetaGeomMean, FScalar fAtoB,
   FScalar fAtoWC, FScalar fBtoWC, FMemoryStack &Mem )
{
   FScalar
      *pTempStorage;
      // this version works but seems to be numerically unstable. I
      // therefore try to implement the realtions again in exactly the way
      // they are written down in here ((9.3.8) + (9.3.9)).
   /*
   Mem.ClearAlloc(pTempStorage, (AngMomA + AngMomB + 1) * (AngMomB + 1));
   FScalar
      *Target = &pTempStorage[0];

   // S_00 overlap integral (first index: a-cartesian power, second: b's),
   // i.e. overlap integral over s-functions.
   *Target = std::sqrt( M_PI/ZetaSum ) * std::exp( -ZetaGeomMean * fAtoB * fAtoB );

   // Obara-Saika recurrence relations ((9.3.8) in purple book):
   //   S[i+1,j] = X_PA S[i,j] + 1/(2p) ( i S[i-1,j] + j S[i,j-1] )
   //   S[i,j+1] = X_PB S[i,j] + 1/(2p) ( i S[i-1,j] + j S[i,j-1] )
   // horitontal recurrence relation (9.3.7):
   //   S[i,j+1] - S[1+i,j] = X_AB S_ij
   // a note: in the purple book X_UV means U_x - V_x, i.e. VtoW_x. therefore
   // the minus signs here and there.

   // first generate all overlap integrals of type S_i0. i = 1..AngMomA.
   // note: we need one additional i-terms for the horizontal recurrence
   // relations for each j we calculate -> we need S_i0 with
   // i = 1.. (AngMomA+AngMomB).
   for ( uint i = 0; i < AngMomA + AngMomB; ++ i ){
      // note that j=0 -> the j*S[i,j-1] term does not apply (not even exist).
      Target[1+i] = (fAtoWC) * Target[i];
      if ( i != 0 )
         Target[1+i] += 1.0/(2.0*ZetaSum)*( i * Target[i-1] );
   }

   uint
      StrideJ = 1 + AngMomA + AngMomB;

   #define S(i,j) (Target[((j)*StrideJ) + (i)])
   // now use the horizontal recurrence relations to build the rest of the
   // overlap terms from left to right (S_ij with j = 1 .. AngMomB)
   for ( uint j = 0; j < AngMomB; ++ j ){
      // S[i,j+1] = X_AB S_ij + S[1+i,j]
      for ( uint i = 0; i <= AngMomA + AngMomB - j; ++ i ){
         S(i,j+1) = (-fAtoB) * S(i,j) + S(1+i,j);
      }
   }

   // copy stuff around to get rid of intermediate data of recurrence
   // relations.
   for ( uint j = 0; j <= AngMomB; ++ j )
      for ( uint i = 0; i <= AngMomA; ++ i ){
         RealTarget[ (1+AngMomA) * j + i ] = S(i,j);
      }

   #undef S
   */

   // hm... nop, this did not make any difference at all. the upper version
   // produces exactly the same results.
   uint
      AngMomMax = std::max( AngMomA, AngMomB );
   Mem.ClearAlloc(pTempStorage, (AngMomMax + 1) * (AngMomMax + 1));
   FScalar
      *Target = &pTempStorage[0];
   // S_00 overlap integral (first index: a-cartesian power, second: b's),
   // i.e. overlap integral over s-functions.
   *Target = std::sqrt( M_PI/ZetaSum ) * std::exp( -ZetaGeomMean * sqr(fAtoB) );

   uint
      StrideJ = 1 + AngMomMax;
   #define S(i,j) (Target[((j)*StrideJ) + (i)])

   //   S[i+1,j] = X_PA S[i,j] + 1/(2p) ( i S[i-1,j] + j S[i,j-1] )
   //   S[i,j+1] = X_PB S[i,j] + 1/(2p) ( i S[i-1,j] + j S[i,j-1] )
   FScalar
      f = 1/(2.0 * ZetaSum);

   for ( uint n = 0; n < AngMomMax; ++ n ){
      // All S[i,j] with i,j <= n present.
      // Make n+1 th "rectangular shell" out of this.
      for ( uint j = 0; j <= n; ++ j ){
         S(1+n,j) = fAtoWC * S(n,j);
         if ( n != 0 )
            S(1+n,j) += n * f * S(n-1,j);
         if ( j != 0 )
            S(1+n,j) += j * f * S(n,j-1);
      }

      for ( uint i = 0; i <= n+1; ++ i ){
         S(i,1+n) = fBtoWC * S(i,n);
         if ( i != 0 )
            S(i,1+n) += i * f * S(i-1,n);
         if ( n != 0 )
            S(i,1+n) += n * f * S(i,n-1);
      }
   }



   // copy stuff around to get rid of intermediate data of recurrence
   // relations.
   for ( uint j = 0; j <= AngMomB; ++ j )
      for ( uint i = 0; i <= AngMomA; ++ i ){
         RealTarget[ (1+AngMomA) * j + i ] = S(i,j);
      }

   #undef S

   Mem.Free(pTempStorage);
}

// Target must have space for (1+AngMomA)(1+AngMomB)(1+Moment) scalars.
// Target[i + j*(1+AngMomA) + m*(1+AngMomA)(1+AngMomB)] will receive the
// 1d-moment integral S[m;i,j]
void EvalPrimitiveMultipoleMoment1D( FScalar *Target,
   uint Moment, uint AngMomA, uint AngMomB, FScalar ZetaSum, FScalar ZetaGeomMean,
   FScalar fAtoB, FScalar fAtoWC, FScalar fBtoWC, FScalar fMtoWC, FMemoryStack &Mem )
{
   uint
      StrideJ = (1 + AngMomA),
      StrideM = (1 + AngMomA) * (1 + AngMomB);
   #define S(Mom,PowA,PowB) ( Target[(PowB)*StrideJ + (PowA) + (Mom)*StrideM ]  )

   // first generate all shell-moments for Moment = 0 (i.e. 1d-overlap
   // integrals). the higher moment terms are then constructed out of this
   // using the diagonal recurrence relations as in (9.3.19).
   EvalOverlapIntegrals1D( Target, AngMomA, AngMomB, ZetaSum, ZetaGeomMean,
      fAtoB, fAtoWC, fBtoWC, Mem );

   // (9.3.19): S[e+1; i,j] = X_PC S[e;i,j] +
   //   1/(2p)*(i S[e;i-1,j] + j S[e,i,j-1] + e S[e-1;i,j])
   for ( uint e = 0; e < Moment; ++ e ){
      for ( uint j = 0; j <= AngMomB; ++ j )
         for ( uint i = 0; i <= AngMomA; ++ i ){
            FScalar
               f = 1.0/(2.0*ZetaSum);
            FScalar
               &s = S(1+e,i,j);
            s = (fMtoWC) * S(e,i,j);
            if ( e != 0 )
               s += f * e * S(e-1,i,j);
            if ( i != 0 )
               s += f * i * S(e,i-1,j);
            if ( j != 0 )
               s += f * j * S(e,i,j-1);
         }
   }
   #undef S
}


// Target must have space for (1+AngMomA)(1+AngMomB)(1+Derivative) scalars.
// Target[i + j*(1+AngMomA) + e*(1+AngMomA)(1+AngMomB)] will receive the
// 1d-derivative integral S[e;i,j]
void EvalPrimitiveDerivativeIntegrals1D( FScalar *Target,
   uint Derivative, uint AngMomA, uint AngMomB, FScalar ZetaA, FScalar ZetaB,
   FScalar ZetaSum, FScalar ZetaGeomMean,
   FScalar fAtoB, FScalar fAtoWC, FScalar fBtoWC,
   FMemoryStack &Mem )
{
   /*FScalarArray
      TempX;
   TempX.resize( (1 + Derivative)*(AngMomB+1)*(AngMomA+Derivative+1), 32.0 );
   #define D(Deriv,PowA,PowB) ( TempX[ (PowA) + ( 1 + AngMomA + Derivative )*( (PowB) + (1 + AngMomB)*(Deriv) ) ] )

   assert_rt( Derivative == 2 );
   EvalOverlapIntegrals1D( &TempX[0], AngMomA + Derivative, AngMomB,
      ZetaSum, ZetaGeomMean, fAtoB, fAtoWC, fBtoWC, TempStorage );

   uint
      StrideJ = (1 + AngMomA),
      StrideD = (1 + AngMomA) * (1 + AngMomB);
   #define Dx(Deriv,PowA,PowB) ( Target[(PowB)*StrideJ + (PowA) + (Deriv)*StrideD ] )

   for ( uint j = 0; j <= AngMomB; ++ j )
      for ( uint i = 0; i <= AngMomA; ++ i ){
         Dx(0,i,j) = D(0,i,j);
         Dx(1,i,j) = 2*ZetaA*D(0,i+1,j);
         if ( i != 0 )
            Dx(1,i,j) = i * D(0,i-1,j);
         Dx(2,i,j) = 4*sqr(ZetaA)*D(0,2+i,j) - 2*ZetaA*(2*i+1)*D(0,i,j);
         if ( i > 1 )
            Dx(2,i,j) += i*(i-1)*D(0,i-2,j);
      }
   #undef Dx
   #undef D*/

   // this one in fact did make a difference in numerical stability. This
   // has to be investigated further by analytically evaluating some integrals.
   // Note: the actual energy results (RHF, Corr) are pretty stable regarding
   // slight changes here, but properties might be not.

   uint
      StrideJ = (1 + AngMomA),
      StrideD = (1 + AngMomA) * (1 + AngMomB);
   #define D(Deriv,PowA,PowB) ( Target[(PowB)*StrideJ + (PowA) + (Deriv)*StrideD ] )

   // first generate all shell-moments for Derivative = 0 (i.e. 1d-overlap
   // integrals). the higher Derivative terms are then constructed out of this
   // using the diagonal recurrence relations as in (9.3.28) where D[e;i+1,j]
   // has been substituted by (9.3.26).
   EvalOverlapIntegrals1D( Target, AngMomA, AngMomB, ZetaSum, ZetaGeomMean,
      fAtoB, fAtoWC, fBtoWC, Mem );

   // (9.3.28) with (9.3.26): D[e+1;i,j] = 2 ZetaA X_PA D[e;i,j] +
   //      2 ZetaA/(2p)*( i D[e;i-1,j] + j D[e;i,j-1] - 2 ZetaB e D[e-1;i,j] )
   //      - i D[e;i-1,j]
   for ( uint e = 0; e < Derivative; ++ e ){
      for ( uint j = 0; j <= AngMomB; ++ j )
         for ( uint i = 0; i <= AngMomA; ++ i ){
            FScalar
               f = ZetaA/ZetaSum;
            FScalar
               &s = D(1+e,i,j), x;
            s = 2.0 * ZetaA * (fAtoWC) * D(e,i,j);
//            if ( i != 0 )
//               s += f * i * D(e,i-1,j);
            x = 0;
            if ( i != 0 )
               x += (f - 1.0) * i * D(e,i-1,j);
            if ( j != 0 )
               x += f * j * D(e,i,j-1);
            if ( e != 0 )
               x -= f * 2.0 * ZetaB * e * D(e-1,i,j);
            s += x;
//            if ( i != 0 )
//               s -= i * D(e,i-1,j);

         }
   }
   #undef D
}




// sets contribution of the given primitive cartesian shell to integrals at
// Result[AngCompA + AngCompB*NumAngCompA].
void FDoublettIntegralFactoryMultipoleMoment::EvalPrimitiveIntegral(
      FScalar *Result, FPrimitiveIntegralInfo const &info )
{
//   xout << "-------- begin" << std::endl;
//   xout << "AngCompA[0]: " << (int)AngCompsA[0].m[0] << " " << (int)AngCompsA[0].m[1] << " " << (int)AngCompsA[0].m[2] << ";" << std::endl;
//   xout << "AngCompB[0]: " << (int)AngCompsB[0].m[0] << " " << (int)AngCompsB[0].m[1] << " " << (int)AngCompsB[0].m[2] << ";\n" << std::endl;

   // This is a naive implementation of the Obara-Saika scheme for multipole
   // moment evaluation as described in the purple book on p.347 and prior.
   // (note especially the def. of the quantities in the recurrence
   // relations on p.341)


   // we'll need the results for each cartesian direction and for all cartesian
   // powers of the shell functions <= AngMomA/AngMomB to build the complete
   // shell integral. As all the lower cartesian powers are produced as
   // by-products in the evaluation of higher powers anyway, it is very
   // beneficial to do this shell-wise.
   FAngularMomentum const
      &AngMomA = info.AngMomA,
      &AngMomB = info.AngMomB;
   FVector3
      vMtoWC = info.vWeightedCenter - this->vMomentCenter;


   size_t
      StrideM = (1 + AngMomA) * (1 + AngMomB),
      BaseM[3];

   FMoment
      vMaxMoment(0,0,0);
   for (size_t iComp = 0; iComp < vMoments.size(); ++ iComp)
      for (size_t iCart = 0; iCart <= 2; ++ iCart)
         vMaxMoment[iCart] = std::max(vMaxMoment[iCart], vMoments[iComp][iCart]);

   FScalar
      *pMultipoles1D;
   info.Mem.ClearAlloc(pMultipoles1D, StrideM * (3 + vMaxMoment[0] + vMaxMoment[1] + vMaxMoment[2]));

   // get the 1-d integrals for all directions.
   for ( size_t nCart = 0; nCart <= 2; ++ nCart ){
      BaseM[nCart] = 0;
      for ( size_t m = 0; m < nCart; ++ m )
         BaseM[nCart] += StrideM * (1 + vMaxMoment[m]);
      FScalar
         *Target = &pMultipoles1D[BaseM[nCart]];

      EvalPrimitiveMultipoleMoment1D(Target, vMaxMoment[nCart],
         info.AngMomA, info.AngMomB, info.ZetaSum, info.ZetaGeomMean,
         info.vAtoB[nCart], info.vAtoWC[nCart], info.vBtoWC[nCart],
         vMtoWC[nCart], info.Mem);
   }


   size_t
      StrideJ = (1 + AngMomA);
   #define S(Dir,Mom,PowA,PowB) ( pMultipoles1D[(PowB)*StrideJ + (PowA) + (Mom)*StrideM + BaseM[Dir]] )

   // total multipole moment: just product of multipole moments in each single
   // direction.
   size_t
      nAcA = info.AngCompsA.size(),
      nAcB = info.AngCompsB.size();
   for ( size_t iComp = 0; iComp < vMoments.size(); ++ iComp)
   {
      FMoment
         vMoment = vMoments[iComp];
      for ( size_t j = 0; j < info.AngCompsB.size(); ++ j )
      {
         for ( size_t i = 0; i < info.AngCompsA.size(); ++ i )
         {
            FAngularComp const
               &AngCompA = info.AngCompsA[i],
               &AngCompB = info.AngCompsB[j];
            FScalar
               r = 1.0;
            for ( size_t nCart = 0; nCart <= 2; ++ nCart )
               r *= S( nCart, vMoment[nCart], AngCompA[nCart], AngCompB[nCart] );

//         Result[(i)+(j)*info.AngCompsA.size()] += r;
            Result[(i)+nAcA*((j) + nAcB*iComp)] = r;
         }
      }
   }
   #undef S
   info.Mem.Free(pMultipoles1D);
}

void FDoublettIntegralFactoryKineticTerm::EvalPrimitiveIntegral( FScalar *Result,
      FPrimitiveIntegralInfo const &info )
{
   // for this we need overlap integrals and 2-snd derivative integrals.
   // what we want to calculate is:
   //      T = <a| -0.5 (\del_x^2 + \del_y^2 + \del_z^2) |b>
   //        = -0.5 ( <a|\del_x^2|b> + .. + <a|\del_z^2|b> ).
   // since the derivatives are only applied in one cartesian direction each,
   // a term like <a|\del_x^2|b> factors into the 1-D integrals
   // <ax|\del_x^2|bx>*<ay|by>*<az|bz>.

   FAngularMomentum const
      &AngMomA = info.AngMomA,
      &AngMomB = info.AngMomB;
   TVector3<uint>
      vDerivative( 2, 2, 2 ); // we currently need only 2nd derivatives, but
      // we try to keep this general so that scalar relativistic terms can
      // easily be integrated later.
   uint
      StrideJ = (1 + AngMomA),
      StrideD = (1 + AngMomA) * (1 + AngMomB),
      BaseD[3];
   #define D(Dir,Deriv,PowA,PowB) ( pTempStorage[(PowB)*StrideJ + (PowA) + (Deriv)*StrideD + BaseD[Dir]] )

   FScalar
      *pTempStorage;
   info.Mem.ClearAlloc(pTempStorage, StrideD *
         (3 + vDerivative[0] + vDerivative[1] + vDerivative[2]) );

   // get the 1-d integrals for all directions and all derivatives \leq
   // vDerivative[nCart] for that direction.
   for ( uint nCart = 0; nCart <= 2; ++ nCart ){
      BaseD[nCart] = 0;
      for ( uint m = 0; m < nCart; ++ m )
         BaseD[nCart] += StrideD * (1 + vDerivative[m]);
      FScalar
         *Target = &pTempStorage[BaseD[nCart]];

      EvalPrimitiveDerivativeIntegrals1D( Target, vDerivative[nCart],
         info.AngMomA, info.AngMomB, info.ZetaA, info.ZetaB,
         info.ZetaSum, info.ZetaGeomMean,
         info.vAtoB[nCart], info.vAtoWC[nCart], info.vBtoWC[nCart],
         info.Mem );
   }

   // build kinetic term by summing up the 3-products (2 overlap, 1 deriv)
   // for all three cartesian directions
   for ( uint j = 0; j < info.AngCompsB.size(); ++ j )
      for ( uint i = 0; i < info.AngCompsA.size(); ++ i )
      {
         FAngularComp const
            &AngCompA = info.AngCompsA[i],
            &AngCompB = info.AngCompsB[j];
         FScalar
            Derivs2[3], // <ar|\del_r^2|ar>   r = x,y,z
            Overlap[3]; // <ar|1|ar>
         for ( uint nCart = 0; nCart <= 2; ++ nCart ){
            Derivs2[nCart] = D( nCart, 2, AngCompA[nCart], AngCompB[nCart] );
            Overlap[nCart] = D( nCart, 0, AngCompA[nCart], AngCompB[nCart] );
         }
         Result[(i)+(j)*info.AngCompsA.size()] =
            -0.5 * ( Derivs2[0]*Overlap[1]*Overlap[2] +
                   Overlap[0]*Derivs2[1]*Overlap[2] +
                   Overlap[0]*Overlap[1]*Derivs2[2] );
      }

   #undef D
   info.Mem.Free(pTempStorage);
}


// used to evaluate electron densities
//         <a| \delta(r - R) |b>
// with R provided.
void FDoublettIntegralFactoryDensity::EvalPrimitiveIntegral( FScalar *Result,
   FPrimitiveIntegralInfo const &info )
{
   // this is a particularily primitive primitive integral driver...
   struct loc{
      // vDiff: EvalPoint - function center
      // ExponentialTerm: std::exp( -Zeta * vDiff.LengthSq() )
      static inline FScalar EvalGaussFn( FVector3 const &vDiff,
         FAngularComp const &Type, FScalar const &ExponentialTerm )
      {
         return intpow(vDiff[0u], Type[0u]) * intpow(vDiff[1u], Type[1u]) *
            intpow(vDiff[2u], Type[2u]) * ExponentialTerm;
      }
   };

   FVector3
      vDiffA = m_Pos - info.vCenterA,
      vDiffB = m_Pos - info.vCenterB;
   FScalar
      fExponentialA = std::exp( -info.ZetaA * LengthSq(vDiffA) ),
      fExponentialB = std::exp( -info.ZetaB * LengthSq(vDiffB) );
   // somehow i feel an urge to tabulate the vDiff[x] powers. It is hard to
   // resist, but ... this time i'll manage to do it! (cgk)

   for ( uint j = 0; j < info.AngCompsB.size(); ++ j )
      for ( uint i = 0; i < info.AngCompsA.size(); ++ i )
      {
         FScalar
            // result is just conj(a(R))*b(R), and as we're doing real
            // arithmetic...
            ThisInt = loc::EvalGaussFn( vDiffA, info.AngCompsA[i],
               fExponentialA ) * loc::EvalGaussFn( vDiffB,
               info.AngCompsB[j], fExponentialB );

         Result[(i)+(j)*info.AngCompsA.size()] = ThisInt;
      }
}



// will calculate expansion coefficients of gaussian overlap distribution
// into hermite gaussians in one dimension, i.e. the factors E^ij_t in the
// purple book. Target must hold room for (max_i + max_j + 1)*(max_i+1)*(max_j+1)
// scalars and E^ij_t is stored in
//      Target[ t + (i)*(1+max_t) + (1+max_i)*(1+max_t)*(j) ]
// with max_t = max_i + max_j. This is therefore the core for Davidson-
// McMurchie schemes.
void DecomposeCartesianOverlapToHermiteGaussians1D( FScalar *Target,
      uint max_i, uint max_j,
      FScalar ZetaSum, FScalar ZetaGeomMean,
      FScalar AtoB, FScalar AtoWC, FScalar BtoWC,
      FMemoryStack &Mem )
{
   uint
      max_idx = std::max(max_i,max_j),
      max_t = max_i + max_j;

   FScalar
      *pTempStorage;
   Mem.ClearAlloc(pTempStorage, (1+max_idx) * (1+max_idx));
   #define E0(i,j) pTempStorage[ (i) + (1+max_idx)*(j) ]

   // now proceed with step 1, i.e. calculate the cartesian overlap
   // distribution in terms of   hermite gaussians.
   // For this we have the following recurrence relations:
   //    (9.5.8)  E^00_0 = K[dir]
   // as starting point with
   //    (9.2.15) K[dir] = exp(-ZetaGeomMean AtoB[dir]^2 )
   // and for the recursion:
   //   E[i+1,j;0] = X_PA E^ij_0 + E^ij_1
   //   E[i,j+1;0] = X_PB E^ij_0 + E^ij_1
   //   E[i,j;t] = 1/(2pt)( i E[i-1,j;t-1] + j E[i,j-1;t-1] ) for t > 0.

   // inserting the 3rd eq into the 1st and 2nd we arrive at expressions
   // for E[x,x;0] which we can actually calculate:
   //   E[i+1,j;0] = X_PA E^ij_0 + 1/(2p)( i E[i-1,j;0] + j E[i,j-1;0] )
   //   E[i,j+1;0] = X_PB E^ij_0 + 1/(2p)( i E[i-1,j;0] + j E[i,j-1;0] ).
   // do that now.
   E0(0,0) = std::exp( -ZetaGeomMean * AtoB*AtoB );
   for ( uint Radius = 0; Radius < max_idx; ++ Radius ){
      // due to the dependence on prior values of the level 0 recursion
      // relations we need to proceed evaluating them in a kind of unusual
      // manner. after each iteration of the outer loop, all values E(i,j,0)
      // with i, j <= Radius should be in place. these are precisely the
      // ones required for the next iteration.
      for ( uint i = 0; i <= Radius; ++ i )
      {
         uint const &j = Radius;
         FScalar &e = E0(i,j+1);
         e = BtoWC * E0(i,j);
         if ( i != 0 )
            e += (0.5/ZetaSum) * i * E0(i-1,j);
         if ( j != 0 )
            e += (0.5/ZetaSum) * j * E0(i,j-1);
      }
      for ( uint j = 0; j <= Radius + 1; ++ j ) // note the +1 for the corner.
      {
         uint const &i = Radius;
         FScalar &e = E0(i+1,j);
         e = AtoWC * E0(i,j);
         if ( i != 0 )
            e += (0.5/ZetaSum) * i * E0(i-1,j);
         if ( j != 0 )
            e += (0.5/ZetaSum) * j * E0(i,j-1);
      }
   }

   #define E(i,j,t) Target[ t + (i)*(1+max_t) + (1+max_i)*(1+max_t)*(j) ]
   // copy intermediate E0 data into final array with correct format.
   for ( uint i = 0; i <= max_i; ++ i )
      for ( uint j = 0; j <= max_j; ++ j )
         E(i,j,0) = E0(i,j);
   #undef E0

   for ( uint t = 1; t <= max_i + max_j; ++ t )
   { // the rest of the coefficients can be calculated in straight-forward
     // fashion using the 3rd equation alone:
     //   E[i,j;t] = 1/(2pt)( i E[i-1,j;t-1] + j E[i,j-1;t-1] ) for t > 0.
      for ( uint i = 0; i <= max_i; ++ i )
         for ( uint j = 0; j <= max_j; ++ j ){
            FScalar &e = E(i,j,t);
            e = 0.0;
            if ( i != 0 )
               e += 0.5/(t*ZetaSum) * ( i * E(i-1,j,t-1) );
            if ( j != 0 )
               e += 0.5/(t*ZetaSum) * ( j * E(i,j-1,t-1) );
         }
   }

   #undef E
   Mem.Free(pTempStorage);
}


void FDoublettIntegralFactoryFieldTerms::EvalPrimitiveIntegral( FScalar *Result,
   FPrimitiveIntegralInfo const &info )
{
   // this one applies the McMurchie-Davidson scheme for Coulomb integrals,
   // again as described in the purple book (this time not the Obara-Saika
   // scheme because the McMurchie-Davidson thing looks a bit simpler and
   // mainly because it makes it easy to calculate some more integral classes
   // with the same code basically for free).
   // We will need to do the following steps:
   //   - decompose the input cartesian GTOs overlap distribution (i.e. the
   //    product of the two GTOs, which is a sum of GTOs at
   //     info.vWeightedCenter) into Hermite GTOs. This is
   //     accomplished by formulas (9.5.15-.17) for the E^ij_t terms.
   //     the result of this process will be a set of expansion factors E^ij_t
   //     with t = 0..i+j for each cartesian direction (i,j: cartesian powers
   //     of the input GTOs in that direction).
   //  - calculae hermite coulomb integrals using Boys functions of various
   //     kind
   //  - multi-linearily transform the hermite coulomb integrals to the
   //    cartesian expansions.

   uint
      max_t = info.AngMomA + info.AngMomB + m_FieldType[0],
      max_u = info.AngMomA + info.AngMomB + m_FieldType[1],
      max_v = info.AngMomA + info.AngMomB + m_FieldType[2],
      max_idx = std::max( max_t, std::max( max_u, max_v ) ),
      num_idx = 1 + max_idx,
      N = max_t + max_u + max_v;


   // allocate space for values of Boys function and temporary hermite
   // integrals R_tuv in the same array. also we'll need space for the
   // overlap hermite expansion coefficients (more on that later).
   FScalar
      *pTempStorage;
   info.Mem.ClearAlloc(pTempStorage, (1+N) + (N+2)*num_idx*num_idx*num_idx +
      3*(1 + info.AngMomA + info.AngMomB)*(1 + info.AngMomA)*(1 + info.AngMomB));
   //   ^- this might look a bit intimidating but it actually worked the first
   //      time i tried.

   FScalar
      *const HermiteIntegralBase = &pTempStorage[1+N], // just behind the Boys-fn values.
      *const ExpansionBase = &pTempStorage[(1+N) + (N+2)*num_idx*num_idx*num_idx];
#define R(n, t,u,v) HermiteIntegralBase[(v) + num_idx * ( (u) + num_idx * ( (t) + num_idx * (n) ) )]
#define Rf(t,u,v) R(N+1,t,u,v)   // final contracted hermite integrals
   // ^- we waste some space and computation time. not all integrals are really
   // necessary and present. But we don't aim for a super-fast implementation
   // anyway because these are only one-electron integrals. For really fast
   // stuff one would need to write a code generator similar to LIBINT.

   // zero out our contracted hermite integral targets. the point-charge
   // contraction will be done instantly during the calculation of the
   // R_tuv hermite integrals for each single C.
   for ( uint t = 0; t <= max_t; ++ t )
      for ( uint u = 0; u <= max_u; ++ u )
         for ( uint v = 0; v <= max_v; ++ v )
            Rf(t,u,v) = 0.0;


   // we begin with the second step. Calculate hermite integrals
   //   R_tuv(p, R_PC) = R_tuv(info.ZetaSum, vCtoP) for all required
   // t,u,v in the final transformation. This has to happen for each single
   // point charge. The back-transformations can be combined for all at once
   // however.
   FPointChargeList::iterator
      itCharge;
   _for_each( itCharge, m_PointCharges )
   {
      FVector3
         vCtoP = info.vWeightedCenter - itCharge->Pos;

      // eq. (9.9.14); n = 0 .. N; Writes values R[n;0,0,0] for n=0..N
      // (incl N.).
      ir::IrBoysFn(&pTempStorage[0], info.ZetaSum * LengthSq(vCtoP), N, 1.0);
      for ( uint n = 0; n <= N; ++ n )
         R(n, 0,0,0) = itCharge->Charge * std::pow(-2.0 * info.ZetaSum,
               static_cast<double>(n)) * pTempStorage[n];
         // ^- note the charge prefactor. I think it's okay to put it here,
         // it will just propagate linearly to all R-terms of this C.

      // recurrence relations:
      // R[n;t+1,u,v] = t R[n+1;t-1,u,v] + X_PC R[n+1;t,u,v]
      // R[n;t,u+1,v] = t R[n+1;t,u-1,v] + Y_PC R[n+1;t,u,v]
      // R[n;t,u,v+1] = t R[n+1;t,u,v-1] + Z_PC R[n+1;t,u,v]
      //
      // note the subtle index shifts. we go down in a 4d-fashion, increasing
      // the t,u,v ranges at each step and decreasing the N, ranges.
      //   starting at R[N;0,0,0] we're trying to down to
      // R[0;0..max_t,0..max_u,0..max_v], which are the integrals we're really
      // looking for.

      for ( uint i = 0; i < N; ++ i ){
         uint n = static_cast<uint>( N - i - 1 );
         for ( uint t = 0; (t <= i) && (t <= max_t); ++ t )
            for ( uint u = 0; (u <= i) && (u <= max_u); ++ u )
               for ( uint v = 0; (v <= i) && (v <= max_v); ++ v )
               {
                  FScalar
                     rn1 = R(n+1, t,u,v);
                  if ( t != max_t ){
                     FScalar
                        &rt = R(n, t+1,u,v);
                     rt = vCtoP[0] * rn1;
                     if ( t != 0 )
                        rt += t * R(n+1, t-1,u,v);
                  }
                  if ( u != max_u ){
                     FScalar
                        &ru = R(n, t,u+1,v);
                     ru = vCtoP[1] * rn1;
                     if ( u != 0 )
                        ru += u * R(n+1, t,u-1,v);
                  }
                  if ( v != max_v ){
                     FScalar
                        &rv = R(n, t,u,v+1);
                     rv = vCtoP[2] * rn1;
                     if ( v != 0 )
                        rv += v * R(n+1, t,u,v-1);
                  }
               }
      }

      // required R-integrals for this C should be in place now. Add them
      // the the others. (note that the charge prefactor has already
      // been accounted for as a prefactor to the Boys function)
      for ( uint t = 0; t <= max_t; ++ t )
         for ( uint u = 0; u <= max_u; ++ u )
            for ( uint v = 0; v <= max_v; ++ v )
               Rf(t,u,v) += R(0,t,u,v);
   }

   // now proceed with step 1, i.e. calculate the cartesian overlap
   // distribution in terms of   hermite gaussians.
   //      Target[ t + (i)*(1+max_t) + (1+max_i)*(1+max_t)*(j) ]
   // with max_t = max_i + max_j.
   uint
      StrideI = info.AngMomA + info.AngMomB + 1,
      StrideJ = StrideI * ( info.AngMomA + 1 ),
      StrideDir = StrideJ * ( info.AngMomB + 1 );
   FScalar
      *BaseE[3] = { ExpansionBase, ExpansionBase + StrideDir,
                 ExpansionBase + 2*StrideDir };

#define E(Dir,i,j,t) (BaseE[Dir][(t) + (i)*StrideI + (j)*StrideJ ])
   for ( uint nCart = 0; nCart <= 2; ++ nCart )
      DecomposeCartesianOverlapToHermiteGaussians1D( BaseE[nCart],
         info.AngMomA, info.AngMomB,
         info.ZetaSum, info.ZetaGeomMean,
         info.vAtoB[nCart], info.vAtoWC[nCart], info.vBtoWC[nCart],
         info.Mem );

   // now for the final step, contract the hermitian integrals with the
   // hermitian expansion coefficients of the primitive.
   // For this we have:
   //  (9.9.32) V^efg_ab = (-1)^(e+f+g) 2 Pi/p \sum_tuv E^ab_tuv R[t+e,u+f,v+g]
   // with E^ab_tuv = E^ij_t E^kl_u E^mn_v with i,j cartesian powers of the
   // GTOs in x-direction, k,l cartesian powers in y-direction and so on.
   // (e,f,g) are the components of the field type we're evaluating, i.e.
   // (0,0,0) for nuclear attraction integrals.

   // for now we try a brute force approach. this procedure could be speeded
   // up by using temporaries and partial transformations, similar to what
   // one does in the MP2 integral transformations, but it most probably
   // is not worth the effort.
   FScalar
      Prefactor = intpow(-1.0, m_FieldType[0] + m_FieldType[1] +
            m_FieldType[2]) * 2.0*M_PI/info.ZetaSum *
            (-1.0);
   // ^- one additional -1 for the negativity of the electron charge (i hope).

   for ( uint j = 0; j < info.AngCompsB.size(); ++ j )
      for ( uint i = 0; i < info.AngCompsA.size(); ++ i )
      {
         FAngularComp const
            &AngCompA = info.AngCompsA[i],
            &AngCompB = info.AngCompsB[j];
         FScalar
            r = 0.0;
         for ( uint t = 0; t <= max_t; ++ t )
            for ( uint u = 0; u <= max_u; ++ u )
               for ( uint v = 0; v <= max_v; ++ v )
                  r += E(0, AngCompA[0u], AngCompB[0u], t) *
                      E(1, AngCompA[1u], AngCompB[1u], u) *
                      E(2, AngCompA[2u], AngCompB[2u], v) *
                      Rf(t, u, v);

         Result[ i + j * info.AngCompsA.size() ] = Prefactor * r;
      }

#undef E
#undef Rf
#undef R
   info.Mem.Free(pTempStorage);
}




// void PrintShellAngularConfig( FGaussShell *a )
// {
//    xout << "AngMom " << (int)a->AngularMomentum << " Components: ";
//    for ( uint i = 0; i < a->AngularComps.size(); ++ i )
//       xout << (int)a->AngularComps[i][0] << " "
//           << (int)a->AngularComps[i][1] << " "
//           << (int)a->AngularComps[i][2] << " " << " | ";
//    xout << std::endl;
// }
//
// #include <iostream>
//
// calculates zeta-dependent data in *info and calls EvalPrimitiveIntegral
bool FDoublettIntegralFactory::CalcDataAndEvalPrimitiveIntegral( FScalar *Result,
   FScalar ZetaA, FScalar ZetaB, FPrimitiveIntegralInfo &info )
{
   info.ZetaA = ZetaA;
   info.ZetaB = ZetaB;
   info.ZetaSum = ZetaA + ZetaB;
   info.ZetaGeomMean = (ZetaA * ZetaB) / (ZetaA + ZetaB);

   if (1) {
      // Screen on primitive overlap between s-functions on A and B.
      // FIXME: this is not quite right when using long ranged kernels...
      // However, in practice this function is currently used only for
      // Multipole moments, where it should work.
      double
         fSqAtoB = info.vAtoB.LengthSq();
      if (-info.ZetaGeomMean * fSqAtoB < -40.) {
//          memset(Result, 0, sizeof(Result[0]) * info.AngCompsA.size() * info.AngCompsB.size());
         return false;
      }
   }

   info.vWeightedCenter = (1.0/info.ZetaSum) *
         (ZetaA * info.vCenterA + ZetaB * info.vCenterB);
   info.vAtoWC = info.vWeightedCenter - info.vCenterA;
   info.vBtoWC = info.vWeightedCenter - info.vCenterB;

//    std::cout << "PrimIntegral: |W|^2 = " << LengthSq(info.vWeightedCenter) << std::endl;

   EvalPrimitiveIntegral(Result, info);
   return true;
}


void TransformCartesianDoublettToSpherical( FScalar *Target, FScalar const *Source,
   size_t l0, size_t l1, FMemoryStack &Mem ); // forward.

void ContractPrimitive(FShellDoublettIntegral &Result, FBasisShell const *a, FBasisShell const *b, size_t PrimA, size_t PrimB, size_t iComp, FScalar *pSphericalPrimitive, FMemoryStack &/*Mem*/)
{
   size_t
      nAngularCompsA = (2 * a->l() + 1),
      nAngularCompsB = (2 * b->l() + 1);
#define ContrInt(ConA,ConB,AngCompA,AngCompB) Result.Int(\
   (ConA) * nAngularCompsA + (AngCompA), \
   (ConB) * nAngularCompsB + (AngCompB), \
   iComp)

#define LOOP_CONTR(x,X)   \
      for (size_t iCon##X = 0; iCon##X < x->nCo(); ++ iCon##X){\
         FScalar \
            Coeff##X = x->fCo(Prim##X, iCon##X); \
         if (Coeff##X == 0.) \
            continue;

   LOOP_CONTR(a,A)
      LOOP_CONTR(b,B)
         FScalar
            // weight of current primitive in the
            // contraction integral
            CoeffAB = CoeffA * CoeffB;
         for ( size_t i = 0; i < nAngularCompsA; ++ i )
            for ( size_t j = 0; j < nAngularCompsB; ++ j )
               ContrInt(iConA, iConB, i, j) += CoeffAB *
                  pSphericalPrimitive[i + j * nAngularCompsA];
      }
   }
#undef LOOP_CONTR
#undef ContrInt
}


void FDoublettIntegralFactory::EvalDoublett(
      FShellDoublettIntegral &Result,
      FBasisShell const *a, FBasisShell const *b,
      FMemoryStack &Mem )
{
   TMemoryLock<char> pFreeMe(0, &Mem);
   Mem.Align(32);
// #ifdef _DEBUG
//    a->SanityCheck();
//    b->SanityCheck();
// #endif // _DEBUG

   FAngularCompList
      AngularCompsCartA = AngularComponentsCart(a->l()),
      AngularCompsCartB = AngularComponentsCart(b->l());
   // TODO: get rid of this -^ (with new 2e integral core..)

   // we use this structure to pass information to the primitive calculations
   // in order to reduce code duplication and parameter passing overhead.
   FPrimitiveIntegralInfo
      PrimInfo(AngularCompsCartA, AngularCompsCartB, Mem);
   PrimInfo.AngMomA = a->l();
   PrimInfo.AngMomB = b->l();
   PrimInfo.vAtoB = b->vCenter - a->vCenter;
   PrimInfo.vCenterA = a->vCenter;
   PrimInfo.vCenterB = b->vCenter;

   size_t
      nIntComps = this->nComp();


   // allocate space for result integral doublet. primitive integral data
   // we calculate will be contracted inline into this result.
   Result.nSizeA = a->nFn();
   Result.nSizeB = b->nFn();
//    Result.Data.resize(0); // <- to make sure that the next line really zeros
//                       // out all result values.
//    Result.Data.resize(Result.nSizeA * Result.nSizeB, 0.0);
   Result.Data.resize_and_clear(Result.nSizeA * Result.nSizeB * nIntComps);

   size_t
      nAngularCompsA = (2 * a->l() + 1),
      nAngularCompsB = (2 * b->l() + 1),
      // for each contraction, we have one integral for each angular
      // momentum component in each shell.
      nIntegralsPerPrimitive = nAngularCompsA * nAngularCompsB,
      nIntegralsPerPrimitiveCart = AngularCompsCartA.size()
            * AngularCompsCartB.size();

   // generate the primitive integrals one by one and instantly contract them.
   // the integrals are always generated in cartesian form and stored in
   // TempStorage0 and onward. If solid harmonics GTOs are used, additional
   // space for the transformation is required.
//    FScalar
//       *pTempStorage0;
//    Mem.ClearAlloc(pTempStorage0, nIntegralsPerPrimitiveCart + nIntegralsPerPrimitive);
//    FScalar
//       // target for the cart primitive -> spherical primitive transformation.
//       *pSphericalPrimitive = &pTempStorage0[nIntegralsPerPrimitiveCart];
   FScalar
      *pTempStorage0,
      // target for the cart primitive -> spherical primitive transformation.
      *pSphericalPrimitive;
   Mem.ClearAlloc(pTempStorage0, nIntegralsPerPrimitiveCart * nIntComps);
   Mem.ClearAlloc(pSphericalPrimitive, nIntegralsPerPrimitive);

   for (size_t PrimB = 0; PrimB < b->nExp(); ++PrimB) {
      for (size_t PrimA = 0; PrimA < a->nExp(); ++PrimA) {
         if (!CalcDataAndEvalPrimitiveIntegral(&pTempStorage0[0],
            a->fExp(PrimA), b->fExp(PrimB), PrimInfo))
            continue; // screened out.
         // transform primitive to spherical basis if that should be wished.
         // the same comment applies as in EvalQuartet()
         for (size_t iComp = 0; iComp < nIntComps; ++ iComp) {
            TransformCartesianDoublettToSpherical(pSphericalPrimitive,
               &pTempStorage0[nIntegralsPerPrimitiveCart * iComp], a->l(), b->l(), Mem);
            ContractPrimitive(Result, a, b, PrimA, PrimB, iComp, pSphericalPrimitive, Mem);
         }
      }
   }

   Mem.Free(pTempStorage0);
}


void TransformCartesianDoublettToSpherical(FScalar *pTarget, FScalar const *pSource,
   size_t la, size_t lb, FMemoryStack &Mem)
{
   FScalar
      *pTmp;
   size_t
      nShA = (2*la+1),
      nShB = (2*lb+1),
      nCartA = (la+1)*(la+2)/2;

   Mem.Alloc(pTmp, nCartA * nShB);
   ir::ShTrN(pTmp, pSource, nCartA, lb);
   for ( size_t i = 0; i < nShB; ++ i )
      ir::ShTrN(pTarget + i * nShA, pTmp + i * nCartA, 1, la);
   Mem.Free(pTmp);
}


static FAngularCompList MakeAngularComponentsCart(FAngularMomentum const &l)
{
   FAngularCompList
      r;
   r.reserve(ir::nCartY(l));
   uint
      n = ir::nCartX((int)l-1);
   for ( uint i = 0; i < ir::nCartY(l); ++ i )
      r.push_back(FAngularComp(ir::iCartPow[n+i][0],
                               ir::iCartPow[n+i][1],
                               ir::iCartPow[n+i][2]));
   return r;
}

struct FAngularCompListHolder
{
   std::vector<FAngularCompList>
      AngCompsForL;
   FAngularCompListHolder(uint MaxL) {
      AngCompsForL.reserve(MaxL + 1);
      for (uint i = 0; i <= MaxL; ++ i)
         AngCompsForL.push_back(MakeAngularComponentsCart(i));
   }
};

static FAngularCompListHolder const
   s_AngCompHolder(9);

FAngularCompList const &AngularComponentsCart(FAngularMomentum const &l)
{
   //static FAngularCompListHolder const
   //   s_AngCompHolder(9);
   // ^- needs to be global due to thread safety issues.
   return s_AngCompHolder.AngCompsForL.at(l);
}





} // namespace ct

// kate: space-indent on; tab-indent off; indent-width 3; mixedindent off; indent-mode normal;
