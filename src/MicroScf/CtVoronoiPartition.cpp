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

// #define CT_VORONOI_DEBUG 3
#define CT_VORONOI_DEBUG 0
#ifdef CT_VORONOI_DEBUG
   #include <iostream>
#endif // CT_VORONOI_DEBUG

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

#include "format.h"
#include "CxAtomData.h" // for ElementNameFromNumber
#include "CtVoronoiPartition.h"
#include "CxArgSort.h"
#include "CxPoly.h"

namespace mig {

double GetGridCenter_rExpAvg1(int iElement);
double GetGridCenter_rExpAvg2(int iElement);

template<class FScalar>
static FScalar ClampMin(FScalar const &x, FScalar const &min) {
   return (x < min)? min : x;
}

template<class FScalar>
static FScalar ClampMax(FScalar const &x, FScalar const &max) {
   return (x > max)? max : x;
}

template<class FScalar>
static FScalar Clamp(FScalar const &x, FScalar const &min, FScalar const &max) {
   return ClampMax(ClampMin(x, min), max);
}


// returns the square of a number.
template<class FNumber>
inline FNumber sqr(FNumber const &x) {
   return x*x;
}


// return n!! = n * (n-2) * ... (terminates at (-1)!! := 1 or 0!! := 1 (also 1!! := 1))
constexpr size_t Factorial2(ptrdiff_t n) {
    return n <= 0 ? 1 : n * Factorial2(n-2);
}


namespace step_fn {



double InvPolyStep01(double x) { // "f[x -> 1/2 - sin(asin(1-2x)/3)]"
   return .5 - std::sin(std::asin(1-2*x)/3);
}



// Up next: implementation of PolyStepT_m11 in various variants. This is an explicit form of
// the symmetric beta step function
// Beta[x, 1+nAlpha, 1+nAlpha] for integer nAlpha >= 0,
// and translated from [0,1]-->[0,1] to [-1,1]-->[-1,1]
// (i.e., this evaluates -1 + 2*Beta[(1/2 + x/2), 1+nAlpha, 1+nAlpha])
namespace detail_PolyStepT {
   template<typename FScalar>
   constexpr FScalar PolyStepT_m11_HornerCoeff(int n) {
      return Factorial2(2*ptrdiff_t(n)-1)/FScalar(Factorial2(2*ptrdiff_t(n)));
   }

   template<int n, typename FScalar>
   struct FPolyStepT_m11_HornerCoeff { static FScalar constexpr value = PolyStepT_m11_HornerCoeff<FScalar>(n); };

   template<int n, typename FScalar>
   struct FPolyStepT {
      static FScalar AccTerm(FScalar Accu, FScalar f1mx2) {
         return FPolyStepT<n-1,FScalar>::AccTerm(FPolyStepT_m11_HornerCoeff<n, FScalar>::value + f1mx2 * Accu, f1mx2);
      }
   };
   template<typename FScalar>
   struct FPolyStepT<-1,FScalar> {
      static FScalar AccTerm(FScalar Accu, FScalar f1mx2) {
         return Accu;
         (void)f1mx2; // suppress unused warning.
      }
   };
}

template<unsigned nAlpha, typename FScalar>
double PolyStepT_m11_v1(FScalar x) {
   // note: this ones doesn't need any template stuff. Depends on compiler for figuring
   // out things, however.
   FScalar
      f1mx2 = 1 - x*x,
      fProdTerm = x,
      fAccu = fProdTerm;
   for (unsigned i = 0; i < nAlpha; ++i) {
      FScalar const fCoeff = ((2*i+1)/FScalar(2*i+2));
      fProdTerm *= f1mx2 * fCoeff;
      fAccu += fProdTerm;
   }
   return fAccu;
}


template<unsigned nAlpha, typename FScalar>
double PolyStepT_m11_v2(FScalar x) {
   using detail_PolyStepT::PolyStepT_m11_HornerCoeff;
   FScalar
      f1mx2 = 1 - x*x,
      fAccu = PolyStepT_m11_HornerCoeff<FScalar>(nAlpha);
   for (unsigned i = nAlpha-1; i < nAlpha; --i) {
      fAccu = f1mx2 * fAccu + PolyStepT_m11_HornerCoeff<FScalar>(i);
   }
   return x*fAccu;
}


template<unsigned nAlpha, typename FScalar>
double PolyStepT_m11(FScalar x) {
   using namespace detail_PolyStepT;
   FScalar
      f1mx2 = 1 - x*x;
   return x * FPolyStepT<int(nAlpha)-1,FScalar>::AccTerm(FPolyStepT_m11_HornerCoeff<nAlpha, FScalar>::value, f1mx2);
}

double TrigStep_m11_kgen(double znu, size_t k) {
   // maps [-1,1] to [-1,1] (such that f(-1) = -1 and f(+1)=+1), UNLIKE the Voronoi step functions!).
   // TrigStep1[x_]:=Sin[Pi/2*x]
   // TrigStep[x_,k_]:=If[k==0,x,TrigStep[TrigStep1[x],k-1]]
   // Plot[Table[TrigStep[x,k],{k,0,4}],{x,-1,1}]
   using std::sin;
   double const coeff = M_PI/2;
   double g = znu;
   for (size_t it = 0; it < k; ++ it)
      g = sin(coeff * g);
   return g;
}

double TrigStep_01_kgen(double x, size_t k) {
   // maps [0,1] to [0,1].
   assert_rt(x >= 0 && x <= 1);
   double f = .5 + .5 * TrigStep_m11_kgen(2*x-1, k);
   assert_rt(f >= 0 && f <= 1);
   return f;
}


double PolyStep_01_kgen(double x, size_t k) {
   // maps [0,1] to [0,1].
   assert_rt(x >= 0 && x <= 1);
   for (size_t it = 0; it < k; ++ it) {
      x = (3 - 2*x)*x*x;
   }
   return x;
}


double PolyStep_m11_kgen(double x, size_t k) {
   // maps [-1,1] to [-1,1] (such that f(-1) = -1 and f(+1)=+1), UNLIKE the Voronoi step functions!).
   assert_rt(x >= -1 && x <= 1);
   for (size_t it = 0; it < k; ++ it) {
      x = 1.5*x - .5*x*x*x;
   }
   return x;
}




double PolyStep1(double x) {
   double t = 3*x*x - 2*x*x*x;
//    dt_dx = (3 - 2*x)*(x*x);
   return t;
}

double TrigStep1(double x) {
   using std::sin;
   double const pi = M_PI;
   double const f05_pi = pi/2;
   double t = sqr(sin(f05_pi*x)); // I think the sin form should be more stable (no subtraction)
//    dt_dx = f05_pi*sin(pi*x);
   // ^- FullSimplify[2*Sin[Pi*x/2]*Cos[Pi*x/2]*Pi/2] evaluates to 1/2 Pi Sin[Pi x]
   return t;
}

double InvTrigStep1(double x) {
   using std::acos; using std::sqrt;
   double const inv_pi = 1/M_PI;
   double t = inv_pi * acos(1 - 2*x);
//    dt_dx = inv_pi/sqrt(x - x*x);
   // ^- that is equal to inv_pi/(sqrt(x) * sqrt(1-x)).
   return t;
}

double InvPolyStep1(double x) {
   using std::acos; using std::asin; using std::sin; using std::sqrt;
   double const inv_3 = 1/double(3);
   double inv3_asin_1m2x = inv_3 * asin(1 - 2*x);
   double t = .5 - sin(inv3_asin_1m2x);
//    dt_dx = cos(inv3_asin_1m2x)*inv_3/sqrt(x - x*x);
   return t;
}

} // namespace step_fn




static double NuAbSalvador(double mu, double chi) {
   // JCP 139 071103 (2013); doi: 10.1063/1.4818751 eq. (9)
   double t = chi * (1. - mu);
   return (1. + mu - t) / (1. + mu + t);
}


static double NuAbOriginalBecke(double Mu, double Chi) {
//    assert_rt(0); // until the multiplicative vs additive trafo has been figured out.
   return Mu + Chi * (1.0 - Mu*Mu); // [3], eq. (A2)
   // Note: ^- For this to fit to the a_ij // u_ij // wtf transforms in eqs 9 -- 12 of TA,
   // the "Chi" here must be replaced by "a_ij" in the paper. That is,
   // We'd get:  Chi(here) = uij/(sqr(uij) - 1)  with uij = (xij - 1)/(xij + 1)  and xij = R_i/R_j
   // or xij = sqrt(R_i/R_j)  (in the TA modification).
   // Note: Neither says anything about "do this only if atoms are close together".
}


static double NuAbNoAdjustment(double Mu, double Chi) {
   return Mu;
   (void)Chi; // suppress unused warning.
}



static double VoronoiWeight_Becke_k4(double znu) {
   // Compute the smoothed & truncated voronoi polynomial
   // according to Becke's original prescription (JCP 88 2547 (1988)).
   // Here done with k=4 (stiffness parameter) as recommended in JCP 139 071103 (p.2, bottom right).
   // cgk note: This is similar, but not identical, to a scaled PolyStep[x,a,a] with a=2^n
   double g = znu;
   g = ((3./2.)*g - .5*g*g*g);
   g = ((3./2.)*g - .5*g*g*g);
   g = ((3./2.)*g - .5*g*g*g);
   g = ((3./2.)*g - .5*g*g*g);
   return .5 * (1. - g);
}






static double VoronoiWeight_PolyStep_a16(double znu) {
   return .5 - .5 * step_fn::PolyStepT_m11<16, double>(znu);
}


static double VoronoiWeight_PolyStep_a4sq(double znu) {
   using step_fn::PolyStepT_m11;
   return .5 - .5 * PolyStepT_m11<4, double>(PolyStepT_m11<4, double>(znu));
}

static double VoronoiWeight_PolyStep_a12(double znu) {
   using step_fn::PolyStepT_m11;
   return .5 - .5 * PolyStepT_m11<12, double>(znu);
}

static double VoronoiWeight_PolyStep_a10(double znu) {
   using step_fn::PolyStepT_m11;
   return .5 - .5 * PolyStepT_m11<10, double>(znu);
//    return .5 - .5 * PolyStepT_m11<7, double>(znu);
}


static double VoronoiWeight_TrigStep_k4(double znu) {
   // iterated trigstep... this one is a cgk special. It initially
   // maps [-1,1] to [-1,1]. Then it gets changes to map
   // -1 to +1 and +1 to 0, for compatibility with the other voronoi
   // thingies.
   using std::sin;
   double const coeff = M_PI/2;
   double g = znu;
   g = sin(coeff * g);
   g = sin(coeff * g);
   g = sin(coeff * g);
   g = sin(coeff * g);
   return .5 * (1. - g);
}





static double VoronoiWeight_Stratman(double znu) {
   // Stratmann's modified Becke scheme [6]
   double
      g = znu,
      a = 0.64; // comment after eq. 14
   // [6], eq.11
   if (g <= -a)
      g = -1.;
   else if (g >= a)
      g = +1;
   else {
      // [6], eq.14
      double ma = g/a; // mu/a
      g = (1/16.)*(ma*(35 + ma*ma*(-35 + ma*ma*(21 - 5 *ma*ma))));
   }
   return .5 * (1. - g);
}


static double GetNeutralChi(FVoronoiPartitionParams::FAtomSizeAdjustType AtomSizeAdjustType_) {
   switch (AtomSizeAdjustType_) {
      case FVoronoiPartitionParams::ATOMSIZE_Tfvc:
         return 1.;
      case FVoronoiPartitionParams::ATOMSIZE_OriginalBecke:
      case FVoronoiPartitionParams::ATOMSIZE_TaModOfBecke:
      case FVoronoiPartitionParams::ATOMSIZE_None:
         return 0;
      default:
         throw std::runtime_error(fmt::format("GetNeutralChi(): size adjustment type {} not recognized.", unsigned(AtomSizeAdjustType_)));
   }
}


static double GetReverseChi(FVoronoiPartitionParams::FAtomSizeAdjustType AtomSizeAdjustType_, double ChiAB) {
   // computes ChiBA from ChiAB
   switch (AtomSizeAdjustType_) {
      case FVoronoiPartitionParams::ATOMSIZE_Tfvc:
         return 1/ChiAB;
      case FVoronoiPartitionParams::ATOMSIZE_OriginalBecke:
      case FVoronoiPartitionParams::ATOMSIZE_TaModOfBecke:
      case FVoronoiPartitionParams::ATOMSIZE_None:
         return -ChiAB;
      default:
         throw std::runtime_error(fmt::format("GetReverseChi(): size adjustment type {} not recognized.", unsigned(AtomSizeAdjustType_)));
   }
}



FVoronoiPartitionParams::FVoronoiPartitionParams(FSmoothingFnType SmoothingFnType_, FAtomSizeAdjustType AtomSizeAdjustType_)
   : m_SmoothingFnType(SmoothingFnType_), m_AtomSizeAdjustType(AtomSizeAdjustType_)
{
   m_ThrPairVoronoiWeightCut = 1e-12;
   m_fWeightCut_AtomVdwRadiusFactor = 0.; // <-- must be set explicitly if wanted. Would break TFVC if auto-enabled.
   m_CutPairVoronoiAtMaxAtRange = false;
}


FVoronoiPartition::FVoronoiPartition(ct::FRawAtomList const &Atoms_, FVoronoiPartitionParams const &Params_, double const *pChiAB_, size_t iRowSt_, size_t iColSt_)
   : m_Atoms(Atoms_), m_Params(Params_)
//    m_SmoothingFnType(Params_.m_SmoothingFnType), m_AtomSizeAdjustType(Params_.m_AtomSizeAdjustType)
{
//    if (pAtomRadii_)
//       m_ExplicitAtomRadii.assign(pAtomRadii_, pAtomRadii_ + m_Atoms.size());

//    SetupVoronoiInfo(ChiAB_);
   switch (Params_.m_SmoothingFnType) {
      case FVoronoiPartitionParams::SMOOTHFN_Becke_k4: { m_pSmoothingFn = VoronoiWeight_Becke_k4; break; }
      case FVoronoiPartitionParams::SMOOTHFN_PolyStep_a16: { m_pSmoothingFn = VoronoiWeight_PolyStep_a16; break; }
      case FVoronoiPartitionParams::SMOOTHFN_PolyStep_a12: { m_pSmoothingFn = VoronoiWeight_PolyStep_a12; break; }
      case FVoronoiPartitionParams::SMOOTHFN_PolyStep_a10: { m_pSmoothingFn = VoronoiWeight_PolyStep_a10; break; }
      case FVoronoiPartitionParams::SMOOTHFN_Stratman: { m_pSmoothingFn = VoronoiWeight_Stratman; break; }
      case FVoronoiPartitionParams::SMOOTHFN_PolyStep_a4sq: { m_pSmoothingFn = VoronoiWeight_PolyStep_a4sq; break; }
      case FVoronoiPartitionParams::SMOOTHFN_TrigStep_k4: { m_pSmoothingFn = VoronoiWeight_TrigStep_k4; break; }

      default:
         throw std::runtime_error(fmt::format("FVoronoiPartition::c'tor: encountered unknown smoothing fn type '{}'.", int(Params_.m_SmoothingFnType)));
   }

   switch (Params_.m_AtomSizeAdjustType) {
      case FVoronoiPartitionParams::ATOMSIZE_Tfvc: { m_pAtomSizeAdjustFn = NuAbSalvador; break; }
      case FVoronoiPartitionParams::ATOMSIZE_OriginalBecke: { m_pAtomSizeAdjustFn = NuAbOriginalBecke; break; }
      case FVoronoiPartitionParams::ATOMSIZE_TaModOfBecke: { m_pAtomSizeAdjustFn = NuAbOriginalBecke; break; }
      case FVoronoiPartitionParams::ATOMSIZE_None: { m_pAtomSizeAdjustFn = NuAbNoAdjustment; break; }

      default:
         throw std::runtime_error(fmt::format("FVoronoiPartition::c'tor: encountered unknown atomic size adjustment fn type '{}'.", int(Params_.m_AtomSizeAdjustType)));
   }


   SetupVoronoiInfo(pChiAB_, iRowSt_, iColSt_);
   if (m_Params.m_CutPairVoronoiAtMaxAtRange && Atoms_.size() != m_Params.m_MaxAtRanges.size())
      throw std::runtime_error("FVoronoiPartition: if thresholded pair voronoi weighting is to be used, the atomic ranges must be provided in FVoronoiPartitionParams::m_MaxAtRanges.");
}


// void FVoronoiPartition::SetupVoronoiInfo(FMatrixView ChiAB_)
void FVoronoiPartition::SetupVoronoiInfo(double const *pChiAB_, size_t iRowSt_, size_t iColSt_)
{
   size_t
      nAt = m_Atoms.size();

   // compute inverse distance between all atoms.
   m_pInvDistAt.resize(nAt*nAt);
   for (size_t iAt = 0; iAt < nAt; ++ iAt)
      for (size_t jAt = 0; jAt < nAt; ++ jAt) {
         double rij = Dist(m_Atoms[iAt].vPos, m_Atoms[jAt].vPos);
         if (m_Params.m_CutPairVoronoiAtMaxAtRange) {
            // A variant of the core idea in Laqua, Kussmann, Ochsenfeld,
            // J. Chem. Phys. 149, 204111 (2018); https://doi.org/10.1063/1.5049435
            // ...very clever. We just truncate the maximum smoothing distance
            // between atomic pair voronoi partitions, and assign 100% or 0%
            // weight to atoms not within the buffer zone. This is still an
            // exact partition, and has no effect if atoms are close.
            // The Ochsenfeld paper above suggests an universal cutoff of 5 bohr.
            // We give a bit more freedom, and choose the maximum buffer size based
            // on (typically) scaled Rahm-style iso-density vdW radii.
            rij = std::min(rij, m_Params.m_MaxAtRanges[iAt] + m_Params.m_MaxAtRanges[jAt]);
         }
         if (rij != 0)
            rij = 1./rij;
         else
            rij = 0.;
         m_pInvDistAt[nAt*iAt + jAt] = rij;
         m_pInvDistAt[nAt*jAt + iAt] = rij;
      }

   m_ChiAB.resize(nAt * nAt);
   // flatten and transpose input matrix for better access pattern.
   double const NeutralChi = GetNeutralChi(m_Params.m_AtomSizeAdjustType);
   for (size_t iAt = 0; iAt != nAt; ++ iAt) {
      m_ChiAB[nAt*iAt + iAt] = NeutralChi;
      for (size_t jAt = 0; jAt != nAt; ++ jAt) {
         double chi = pChiAB_? pChiAB_[iRowSt_ * iAt + iColSt_ * jAt] : NeutralChi;
         m_ChiAB[iAt*nAt + jAt] = chi;
//          m_ChiAB[jAt*nAt + iAt] = -chi;
      }
   }
}


void FVoronoiPartition::SortAtomsByDistanceToPoint(double *&pDistAg, size_t *&pAtOrd, size_t &nAtOut, FVec3d const &vGridPos, ct::FMemoryStack &Mem) const
{
   size_t
      nAt = m_Atoms.size();
   Mem.Alloc(pAtOrd, nAt);
   Mem.Alloc(pDistAg, nAt);

   ct::TMemoryLock<double>
      pDistAg1(nAt, &Mem);
   // compute distance from grid point to all other atoms.
   for (size_t iAt = 0; iAt < nAt; ++ iAt)
      pDistAg1[iAt] = Dist(vGridPos, m_Atoms[iAt].vPos);

   // sort atoms by distance from grid point. close ones first.
   ct::ArgSort1(pAtOrd, pDistAg1, 1, nAt, false);
   for (size_t iAt = 0; iAt < nAt; ++ iAt)
      pDistAg[iAt] = pDistAg1[pAtOrd[iAt]];

   nAtOut = nAt;
}


double FVoronoiPartition::GetAtomWeight(FVec3d const &vGridPos, size_t iAtom, ct::FMemoryStack &Mem) const
{
   ct::TMemoryLock<char> pFreeMe(0, &Mem);
   double
      *pDistAg;
   size_t
      *pAtOrd, nAt;
   // compute distance from grid point to all other atoms.
   SortAtomsByDistanceToPoint(pDistAg, pAtOrd, nAt, vGridPos, Mem);
   return GetAtomWeight(vGridPos, iAtom, pDistAg, pAtOrd, nAt, Mem);
}


double FVoronoiPartition::GetAtomWeight(FVec3d const &vGridPos, size_t iAtom, size_t *pAtOrd, size_t nAt, ct::FMemoryStack &Mem) const
{
//    s_vDummy = vGridPos;
   ct::TMemoryLock<double>
      pDistAgDirect(nAt, &Mem);
   for (size_t jAt = 0; jAt < nAt; ++ jAt)
      pDistAgDirect[jAt] = Dist(vGridPos, m_Atoms[pAtOrd[jAt]].vPos);
   return GetAtomWeight(vGridPos, iAtom, pDistAgDirect.p, pAtOrd, nAt, Mem);
}

double FVoronoiPartition::GetAtomWeight(FVec3d const &vGridPos, size_t iAtom, double *pDistAg, size_t *pAtOrd, size_t nAt, ct::FMemoryStack &/*Mem*/) const
{
   using step_fn::PolyStep_01_kgen;
   using std::pow;
   using std::exp;
   double
      ThrPairVoronoiWeightCut = m_Params.m_ThrPairVoronoiWeightCut;
//    ThrPairVoronoiWeightCut = 0;
   double
      wOut = 0.,
      wTotal = 0.;

   for (size_t iCenterAtom_ = 0; iCenterAtom_ != nAt; ++ iCenterAtom_) {
      double const
         r1 = pDistAg[iCenterAtom_];
      double
         wCen = 1.;
      size_t
         iCenterAtom = pAtOrd[iCenterAtom_];
      if (m_Params.m_fWeightCut_AtomVdwRadiusFactor > 0.) {
         double
            R1 = ct::GetVdwRadius_IsoDensity(m_Atoms[iCenterAtom].iElement);
         double
            fCut = m_Params.m_fWeightCut_AtomVdwRadiusFactor,
            RCut = 0;
         R1 *= fCut; // cut at fCut*atomic 0.001 e/a0^3 iso-density radius
         R1 *= 2;    // note @2*R1: the step is at x=0.5
         if (r1 >= RCut + R1)
            wCen = 0;
         else if (r1 > RCut)
            wCen *= 1 - PolyStep_01_kgen((r1-RCut)/R1, 4);
      }

      // TODO: could I cut half of those?
      // UPDATE: yes, I think so. I could first compute the Mu as a nAt x nAt matrix (note that Mu[AB] = 1. - Mu[BA])
      //         and then form the weight loops. It might get a bit messy, though...
      for (size_t iOtherAtom_ = 0; iOtherAtom_ != nAt; ++ iOtherAtom_)  {
         if (wCen < ThrPairVoronoiWeightCut)
            break;
         double const
            r2 = pDistAg[iOtherAtom_]; // <-- note addressing via _-quantity (distance-ordered atoms), not total atom index.
         size_t
            iOtherAtom = pAtOrd[iOtherAtom_];

         if (iCenterAtom == iOtherAtom)
            continue;
         {
            if (1) {
               double
                  // note: in the case of MaxAtRange cutoffs, InvR12 may contain the truncated
                  // inverse distance, rather than the actual inverse atomic distance
                  InvR12 = m_pInvDistAt[iCenterAtom * m_Atoms.size() + iOtherAtom],
                  // ``confocal elliptical coordinates''. A great word.
                  Mu = (r1-r2)*InvR12;
               if (m_Params.m_CutPairVoronoiAtMaxAtRange)
                  Mu = Clamp(Mu, -1.0, 1.0);
               wCen *= GetPairVoronoiR(Mu, iCenterAtom, iOtherAtom);
            } else {
               double const e = 8;
               double z1 = std::pow(r1,e), z2 = std::pow(r2,e);
               double Mu = -1 + 2 * z1/(z1 + z2);
               wCen *= .5 - .5*Mu;
            }
         }
      }
      if (wCen < ThrPairVoronoiWeightCut)
         wCen = 0;
      else
         wCen = std::pow(wCen, .75);
      if (iCenterAtom == iAtom)
         wOut = wCen;
      wTotal += wCen;
   }
   if (wTotal == 0.)
      return 0.;
   return wOut/wTotal;
   (void)vGridPos; // suppress unused warning (already used to compute pDistAg)
}

#ifdef INCLUDE_ABANDONED
//                if (1) {
//                   FVec3d
//                      v1 = m_Atoms[iCenterAtom].vPos,
//                      v2 = m_Atoms[iOtherAtom].vPos;
//                   double
//                      d12 = Dist(v1, v2),
//                      InvD12 = 1/d12;
//                   // compute plane equation for regular voronoi cell.
//                   FVec3d
//                      v12 = v2 - v1,
//                      n12 = v12 * InvD12,
//                      vCen12 = .5*(v1 + v2);
//                   double fPlaneCen = Dot(vCen12, n12);
//                   double fDistGridToPlane = Dot(vGridPos, n12) - fPlaneCen;
//                   double fPlaneStepWidth = 1.*std::min(d12, double(5.));
//                   double // -1: left boundary, +1: right boundary.
//                      MuPlane = fDistGridToPlane/(.5*fPlaneStepWidth);
//                   double xPlane = MuPlane;
//                   xPlane = Clamp(xPlane,-1.,1.);
//                   xPlane = PolyStep_m11_kgen(xPlane,3);
//
// //                   FVec3d vgc = vGridPos - vCen12;
// //                   xPlane = Dot(vgc, n12) / (.25 * std::min(d12, 5.));
// //                   xPlane = Clamp(xPlane,-1.,1.);
// //                   xPlane = PolyStep_m11_kgen(xPlane,3);
//
// //                   r1 = Dot(vGridPos - v1, n12);
// //                   r2 = Dot(vGridPos - v2, n12);
//
//                   double xRadial = 0;
//
//                   if (1) {
//                      double e; int nPolyStepIt = 0;
//                      e = 8; AccWeightExp = 1.;
// //                      e = 16; AccWeightExp = .5;
// //                      auto GetElementRadiusForPairWt = GetGridCenter_rExpAvg1;
// //                      auto GetElementRadiusForPairWt = GetGridCenter_rExpAvg2;
//       //                auto GetElementRadiusForPairWt = ct::GetVdwRadius_IsoDensity;
//       //                auto GetElementRadiusForPairWt = ct::GetCovalentRadius;
//                      auto GetElementRadiusForPairWt = [](double iElement) { return 1; };
// //                      auto gTrafoDist = [&](double r, size_t iAt) { return r/GetElementRadiusForPairWt(m_Atoms[iAt].iElement); };
//                      auto gTrafoDist = [&](double r, size_t iAt) { return r*GetElementRadiusForPairWt(m_Atoms[iAt].iElement); };
//                      // ^-- hm... I think I *should* be dividing those!
//                      { r1 = gTrafoDist(r1,iCenterAtom); r2 = gTrafoDist(r2,iOtherAtom); }
//                      auto gDistFn = [&](double r) { return std::pow(r,e); };
//                      double MuRadial = -1 + 2 * gDistFn(r1)/(gDistFn(r1)+gDistFn(r2));
//
//                      xRadial = MuRadial;
//                   }
//
//                   {
//                      double xTotal;
//                      auto lerp = [](double a, double b, double f){ return (1-f)*a + f*b; };
// //                      xTotal = lerp(xPlane, xRadial, std::sqrt(1-sqr(xPlane)));
// //                      xTotal = lerp(xRadial, xPlane, std::pow(1-sqr(xRadial), .25));
// //                      xTotal = xPlane;
//                      xTotal = xRadial;
//                      xTotal = Clamp(xTotal,-1.,1.);
//                      wCen *= (.5 - .5*xTotal);
//                      AccWeightExp *= 1.;
//                   }
#endif // INCLUDE_ABANDONED



double FVoronoiPartition::GetPairVoronoiR(double Mu, size_t iAtom, size_t iOtherAtom) const
{
   double
      Chi = m_ChiAB[m_Atoms.size() * iAtom + iOtherAtom];
//       Chi = m_ChiAB[m_Atoms.size() * iOtherAtom + iAtom];
      // ^- FIXME: is this in the right order? UPDATE: well... actual TFVC is
      //    working. I don't think it would if this was wrong. unless the ChiAB
      //    computation from atomic radii gets the order wrong... which is quite
      //    possible, too, of course. Argh. Maybe shoud check this out sometime...
      //    (atm size corrections are anyway disabled---they tend to break complicated
      //    molecules the way they are usually done. So not top-priority.)
   double
      Nu = m_pAtomSizeAdjustFn(Mu, Chi);

   return m_pSmoothingFn(Nu);
}



double FindTfvcChiFromLineFraction(double MuOfMinDensity)
{
   double
      // closed-form solution for Solve[nu(muMinRho, chi) = 0, {chi}], with given MuOfMinDensity
      // (note: this requires the Salvador atom size adjustment in the Voronoi partition!)
      ChiAB = (1. + MuOfMinDensity)/(1. - MuOfMinDensity);
   return ChiAB;
}



// Returns line fraction Mu (-1... +1) at which the density becomes minimal.
// note: moved from IvAnalysis.cpp (was: FindMinimumDensityPoint)
double FindMinimumDensityOnLine(double const *pfLineFraction, double const *pfDensity, size_t nPts)
{
   if (nPts == 0)
      throw std::runtime_error("FindMinimumDensityOnLine: no input data!");

   // find index of lowest density input point.
   size_t
      iMin = ct::ArgMin(pfDensity, pfDensity + nPts);
   // fit a polynomial through adjacent points.
   size_t
      nPtsLeftRight = 4,
      iFirst = size_t(std::max(ptrdiff_t(iMin) - ptrdiff_t(nPtsLeftRight), ptrdiff_t(0))),
      iLast = std::min(iMin + nPtsLeftRight, nPts),
      nPointToFit = iLast - iFirst,
      nPolyDegree = nPointToFit - 1; // note: 0'th order polynomial (constant) has 1 degree of freedom.

   ct::TPolynomialFit<double, double>
      PolyFit(&pfDensity[iFirst], 0, &pfLineFraction[iFirst], 0, nPointToFit, nPolyDegree);

   ct::TPolynomialMinimizeParams<double, double>
      MinimizeParams(PolyFit.tFirst, PolyFit.tLast);
   MinimizeParams.Print = false;
   return ct::FindMinimum(PolyFit.p, PolyFit.tMinGuess, MinimizeParams);
}


// compute Chi_AB fuzzy voronoi cell shift parameter from given A--B line density profile and line fractions.
// Line fractions go from (-1. to +1.0).
double FindTfvcChiFromDensityProfile(double const *pfDensity, double const *pfLineFraction, size_t nPts)
{
   double
      // find mu for minimum of density
      muMinRho = FindMinimumDensityOnLine(&pfLineFraction[0], &pfDensity[0], nPts);
   double
      Chi = FindTfvcChiFromLineFraction(muMinRho);
   return Chi;
}



double FindMinimumDensityLineFraction(ct::FRawAtom const &AtA, int iAtA, ct::FRawAtom const &AtB, int iAtB, ct::FDensityModel const *pDensityModel, ct::FMemoryStack &Mem)
{
   FVec3d
      vAtA = AtA.vPos,
      vAtB = AtB.vPos;
//       vCen = 0.5 * (vAtA + vAtB);
   double
      fDistAB = Length(vAtA - vAtB),
      LinearResolution = 20.; // 20 pts per a.u. for sampling on the connecting line

   if (fDistAB < 1e-8) {
      // too close.
      throw std::runtime_error(fmt::format("FindMinimumDensityLineFraction(): illegal coincident centers ({} and {})", iAtA, iAtB));
   }
   size_t
      nGridPt = std::max(size_t(std::min(LinearResolution, LinearResolution * fDistAB)), size_t(1000));
   ct::TMemoryLock<double>
      // grid of points along the connecting line between both atoms
      vGridPt(3 * nGridPt, &Mem),
      // total electron density along the line.
      Rho(nGridPt, &Mem),
      // line fractions along the line (-1.0 -> vAtA; 1.0 -> vAtB)
      Mu(nGridPt, &Mem);
   for (size_t iGridPt = 0; iGridPt != nGridPt; ++ iGridPt) {
      double
         // x = 2*(np.arange(nPt) + 0.5)/float(nPt) - 1.
         mui = 2*(0.5 + double(iGridPt)) / double(nGridPt) - 1.;
         // ^- this will exclude points lying directly on the atoms to
         //    avoid problems with derivative singularities.
      FVec3d
         // vi = (1. - f) * vAtA + f * vAtB;
         vi = 0.5 * ((1. - mui) * vAtA + (1. + mui) * vAtB);

      // remember both the line fraction and the actual grid point.
      Mu[iGridPt] = mui;
      vGridPt[0 + 3*iGridPt] = vi[0];
      vGridPt[1 + 3*iGridPt] = vi[1];
      vGridPt[2 + 3*iGridPt] = vi[2];
   }
   // compute electron density on the grid points.
   int AtomsOfInterest[2] = {iAtA, iAtB};
   pDensityModel->ComputeDensity(&Rho[0], &vGridPt[0], 3, nGridPt, Mem, &AtomsOfInterest[0], 2);

   // GNAAAA... turns out ECP-ed electrons can *absolutely* have density minima
   // in the non-existing core. grmlhmpf. This is an attempted hack around that.
   size_t i0 = 0;
   if (AtA.iElement > 36) {
      // find density maximum in left half of atom-connecting line.
      // The idea is that the inter-nuclear minimum should normally be expected to be
      // between the two density maxima.
      i0 = ct::ArgMax(&Rho[0], &Rho[nGridPt/2]);
      nGridPt -= i0;
   }
   if (AtB.iElement > 36) {
      size_t i1 = nGridPt/2 + 1;
      i1 = ct::ArgMax(&Rho[i1], &Rho[nGridPt]) + i1;
      nGridPt = i1 - i0;
   }
   // Hm, not sure. This seems to work for all cases I have tested. But there still may be a
   // beter way to go about it. Maybe search minimum from the inside out?
   return FindMinimumDensityOnLine(&Mu[i0], &Rho[i0], nGridPt);
}



void FVoronoiChiData::MakeFromDensityModel(ct::FRawAtomList const &Atoms_, ct::FDensityModel const *pDensityModel, FTreatAsNeighborsFn TreatAsNeighborsFn_, FVoronoiPartitionParams::FAtomSizeAdjustType AtomSizeAdjustType_, ct::FMemoryStack &Mem)
{
   size_t
      nAt = Atoms_.size();
   // allocate output array
   m_ChiAB.resize(nAt * nAt);
   m_iRowSt = 1;
   m_iColSt = nAt;

   for (size_t jAt = 0; jAt != nAt; ++ jAt)
      for (size_t iAt = 0; iAt != nAt; ++ iAt)
         ChiAB(iAt, jAt) = 1.;

   assert_rt(AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_Tfvc);

   for (size_t iAt = 0; iAt != nAt; ++ iAt) {
      for (size_t jAt = iAt + 1; jAt != nAt; ++ jAt) {
         ct::FRawAtom const
            &AtA = Atoms_[iAt],
            &AtB = Atoms_[jAt];
         FVec3d
            vAtA = AtA.vPos,
            vAtB = AtB.vPos;
         double
            fDistAB = Length(vAtA - vAtB);
         bool
            TreatAsNeighbors = true;
         if (TreatAsNeighborsFn_ != 0)
            TreatAsNeighbors = TreatAsNeighborsFn_(fDistAB, AtA, -1, AtB, -1);
         if (!TreatAsNeighbors) {
            // atoms too far apart. Put plane at center line.
            ChiAB(iAt, jAt) = 1.;
            ChiAB(jAt, iAt) = 1.;
         } else {
            double
               // line fractions along the line (-1.0 -> vAtA; 1.0 -> vAtB)
               MuDivider = 0.;

            MuDivider = FindMinimumDensityLineFraction(AtA, iAt, AtB, jAt, pDensityModel, Mem);

            for (int iFlatten = 0; iFlatten < 0; ++ iFlatten)
               // with the real density it's just too sharp... e.g., O-H bond is almost completely on O side.
               MuDivider = -1+2*step_fn::InvTrigStep1(.5+.5*MuDivider);

            double
               // convert Mu (line fraction) to Chi (multiplicative parameter of the
               // size adjustment function). Note: MUST fit to the actual size adjustment type!
               Chi = FindTfvcChiFromLineFraction(MuDivider);
            assert(AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_Tfvc);

#if defined(CT_VORONOI_DEBUG) && (CT_VORONOI_DEBUG >= 1)
//             double chi3 = InvTrigStep1(Chi);
//             double mu01 = .5 + .5*MuDivider;
            std::cout << fmt::format("     [{:3} {:2}  {:3} {:2}]  min-rho/chi   dist = {:12.6f}  // mu = {:16.8f}  chi = {:16.8f}  [tr: {:16.8f}  {:16.8f}]", iAt, ct::ElementNameFromNumber(AtA.iElement), jAt,  ct::ElementNameFromNumber(AtB.iElement), fDistAB, MuDivider, Chi, -1+2*InvTrigStep1(.5+.5*MuDivider), -1+2*InvTrigStep1(.5+.5*MuDivider)) << std::endl;
#endif
//             Chi = std::sqrt(Chi); // FIXME: remove this?
//             Chi = std::pow(Chi,3); // FIXME: remove this?
//             Chi = 1/Chi;

            ChiAB(iAt, jAt) = Chi;
            ChiAB(jAt, iAt) = 1./Chi;
         }
      }
   }
}


void FVoronoiChiData::MakeFromAtomRadii(ct::FRawAtomList const &Atoms_, double const *pAtomRadii_, double RSolv_, FTreatAsNeighborsFn TreatAsNeighborsFn_, FVoronoiPartitionParams::FAtomSizeAdjustType AtomSizeAdjustType_)
{
   size_t
      nAt = Atoms_.size();
   // allocate output array
   m_ChiAB.resize(nAt * nAt);
   m_iRowSt = 1;
   m_iColSt = nAt;

   // note: this code based on the code from CtCosmo.cpp, which itself is based
   // on IvAnalysis.cpp's FFragmentSpacePartitioningTfvc.

   double const
      fNeutralChi = GetNeutralChi(AtomSizeAdjustType_);

   // set all Chi_AB to 1 --- this corresponds to cell boundaries at the
   // center of the line connecting both atoms.
   // Apart from the diagonal elements, all other elements will later be re-set
   // in the loop below.
   for (size_t jAt = 0; jAt != nAt; ++ jAt)
      for (size_t iAt = 0; iAt != nAt; ++ iAt)
         ChiAB(iAt, jAt) = fNeutralChi;


   for (size_t iAt = 0; iAt != nAt; ++ iAt) {
      for (size_t jAt = iAt + 1; jAt != nAt; ++ jAt) {
         ct::FRawAtom const
            &AtA = Atoms_[iAt],
            &AtB = Atoms_[jAt];
         FVec3d
            vAtA = AtA.vPos,
            vAtB = AtB.vPos;
//             vCen = 0.5 * (vAtA + vAtB);
         double
            fDistAB = Length(vAtA - vAtB);
         bool
            TreatAsNeighbors = true;
         double
            fRadA = pAtomRadii_[iAt],
            fRadB = pAtomRadii_[jAt];
//          if (fDistAB > fRadA + fRadB)
//             TreatAsNeighbors = false;
         if (TreatAsNeighborsFn_ != 0)
            TreatAsNeighbors = TreatAsNeighborsFn_(fDistAB, AtA, fRadA, AtB, fRadB);
         if (!TreatAsNeighbors) {
            // atoms too far apart. Put plane at center line.
            ChiAB(iAt, jAt) = fNeutralChi;
            ChiAB(jAt, iAt) = fNeutralChi;
         } else {
            double
               Chi = fNeutralChi;
            if (AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_None) {
               Chi = fNeutralChi;
            } else if (AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_OriginalBecke || AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_TaModOfBecke) {
               using std::sqrt;
               double xij = fRadA/fRadB; // // [TA] Eq. 12 (called chi there)
               if (AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_TaModOfBecke)
                  xij = sqrt(xij); // [TA] Eq. 13
               double uij = (xij - 1)/(xij + 1); // // [TA] Eq. 11
               Chi = uij/(sqr(uij) - 1); // [TA] Eq. 10 (called aij there)
               if (Chi < -.5)
                     Chi = -0.5;
               if (Chi > .5)
                     Chi = 0.5;
            } else if (AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_Tfvc) {
               double
                  muRing = 0.;
   //                RSolv_ = 0; // FIXME: REMOVE THIS.
   //                RSolv_ = -2.0;
               if (RSolv_ != -1.) {
   //                assert_rt(RSolv_ >= 0.);
                  // put the plane through the outer ring, i.e., the intersection of the atomic spheres
                  // with radius b = fRadA + RSolv around atom A, and a = fRadB + RSolv around atom B.
                  // (see cosmo_cavity.svgz)
                  //
                  // This gives us a triangle defined by three lengths to get the plane:
                  double
                     a = fRadA + RSolv_, // note: That's the line connecting B and C! (side opposite to A)
                     b = fRadB + RSolv_, // note: That's the line connecting A and C! (side opposite to B)
                     c = fDistAB;
                  // (^- note: changed from drawing, to have pts and sides opposite to each other; trig convention)
                  // FIXME: The a = fRadA + .. and b = rRadB + ... are not consistent with what the
                  //        text says those quantities SHOULD be! Check drawing and formulas again.

                  // and the ought quantity is:
                  //   - First, the length 'x', the projection of the sphere intersection point C
                  //     onto the line 'c' between A and B.
                  //   - And, based on that, the scalar line fraction Mu = 2 * x/c - 1
                  //     (with Mu=-1 landing on A, and Mu = +1 landing on B). That's the
                  //     convention used in the Voronoi Partition code.
                  //
                  // To get those, we need the actual position of C (the sphere intersection point
                  // defining the outer ring position) in the plane. So, let's compute...
                  // Def:
                  //     A = (0,0)
                  //     B = (c,0)
                  //     C = (x,y)
                  // Triangle side lengths:
                  //    |C-A|^2 = x^2 + y^2 = b
                  //    |C-B|^2 = (x - c)^2 + y^2 = a
                  // Given that, Solve[{x^2 + y^2 == b, (x - c)^2 + y^2 == a}, {x,y}] says:
                  //
                  //    x = (c^2 + b - a)/(2c)
                  //    y = +/- sqrt(b - x^2)
                  //
                  // Hooray for triangles \o/.
                  double
                     x = (c*c + b - a)/(2*c);
                  muRing = 2. * x/c - 1.;
               } else {
                  // FIXME: CosmoMode and DftGrid mode give different results even with RSolv_ = 0.
                  // Maybe something wrong? Check trafos!
                  muRing = 2.*(fRadA / (fRadA + fRadB)) - 1.;

                  // ^-- note: mu is supposed to be a line fraction (-1 to 1). I *think* what this
                  // computes is right, at least for two spheres with fixed radii lying right next
                  // to each other and touching.
               }
               Chi = FindTfvcChiFromLineFraction(muRing);
               assert(AtomSizeAdjustType_ == FVoronoiPartitionParams::ATOMSIZE_Tfvc);
               // ^-- that's what FindTfvcChiFromLineFraction is the mu-inverse from
//             Chi = std::sqrt(fRadA/fRadB); // FIXME: remove this.
            } else {
               throw std::runtime_error(fmt::format("FVoronoiChiData::MakeFromAtomRadii(): atomic size adjustment type '{}' not recognized.", AtomSizeAdjustType_));
            }
#if defined(CT_VORONOI_DEBUG) && (CT_VORONOI_DEBUG >= 1)
            std::cout << fmt::format("     [{:3} {:2}  {:3} {:2}]  radii  radA = {:12.6f}  radB = {:12.6f}  dist = {:12.6f}  // mu = {:16.8f}  chi = {:16.8f}", iAt, ct::ElementNameFromNumber(AtA.iElement), jAt,  ct::ElementNameFromNumber(AtB.iElement), fRadA, fRadB, fDistAB, muRing, Chi) << std::endl;
#endif

            ChiAB(iAt, jAt) = Chi;
            ChiAB(jAt, iAt) = GetReverseChi(AtomSizeAdjustType_, Chi);
         }
      }
   }

//    (void)TrigStep1; // suppress unused warning
//    (void)InvPolyStep1; // suppress unused warning
//    (void)PolyStep1; // suppress unused warning
}


} // namespace mig
