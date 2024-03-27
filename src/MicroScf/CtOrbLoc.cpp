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
#include "format.h"
#include <ostream>

#include <iostream>
#include <algorithm>
#include <set>
#include <stdexcept>

#include "CtOrbLoc.h"

#include "CxAlgebra.h"
#include "CxTiming.h"
#include "CxPodArray.h"
#include "CtMatrix.h"
#include "CxRandom.h" // only needed when random rotation code is included.

#ifdef HAVE_BOOST_SPECIAL_FUNCTIONS
   #include <boost/math/special_functions/gamma.hpp>
   #include <boost/math/special_functions/digamma.hpp>
   #include <boost/math/special_functions/trigamma.hpp>
#endif // INCLUDE_OPTIONALS


namespace ct {



FLocalizeOptions::FLocalizeOptions()
   : LocFn(LOCFN_nA_pow_4), nMaxIt(2048), ThrLoc(1e-8), fRandomAngleDeg(18.0), pLog(0), Verbosity(1)
{}


std::string FLocalizeOptions::MakeLocalizeFnDesc() const
{
   switch (this->LocFn) {
      case LOCFN_None: return "L = 0 --> localization disabled";
//       case LOCFN_nA_pow_2: return "L = sum[A,i] (i|n_A|i)^2 --> max";
//       case LOCFN_nA_pow_4: return "L = sum[A,i] (i|n_A|i)^4 --> max";
//       case LOCFN_nA_pow_2: return "L = \\sum_{A,i} (i|n_A|i)^2 --> max";
//       case LOCFN_nA_pow_4: return "L = \\sum_{A,i} (i|n_A|i)^4 --> max";
      case LOCFN_nA_pow_15: return "L = \\sum_{{A,i}} (i|P_A|i)^(3/2) --> max";
      case LOCFN_nA_pow_2: return "L = \\sum_{{A,i}} (i|P_A|i)^2 --> max";
      case LOCFN_nA_pow_4: return "L = \\sum_{{A,i}} (i|P_A|i)^4 --> max";
      // ^- no idea where those {{}} come from... I do not see how they would
      //    possibly end up in a format key. FIXME: find out what this is.
      case LOCFN_nA_gen_h: {
         if (this->pScalarFn01 == 0)
            throw std::runtime_error("FLocalizeOptions::MakeLocalizeFnDesc: general h(n) localization, but pScalarFn01 was not set.");
         return fmt::format("L = \\sum_{{{{A,i}}}} h[(i|P_A|i)], h(n)={} --> max", pScalarFn01->Desc());
      }
//       case LOCFN_nA_Entropy: return "L = sum[A,i] (-n(A,i)*log2(n(A,i))) -> min";
      case LOCFN_nA_Entropy: return "L = \\sum_{{A,i}} -p(A,i)*\\log_2(p(A,i)) --> min";
//       case LOCFN_FockVariance: return "L = \\sum_{i} [(i|f^2|i) - (i|f|i)^2] --> max";
      case LOCFN_FockVariance: return "L = \\sum_{{i}} (<i|f^2|i> - <i|f|i>^2) --> max";
   }
   return fmt::format("FLocalizeOptions::MakeLocalizeFnDesc: LocFn = {} not recognized.", (int)this->LocFn);
}

typedef double (*FRealFn)(double); // a real function taking one double argument and returning a double.
static double pow15(double x) { return x*std::sqrt(x); }
static double pow2(double x) { return x*x; }
static double pow3(double x) { return x*x*x; }
static double pow4(double x) { double x2 = x*x; return x2*x2; }

static double const
   // N[{Log[2],1/Log[2]},73]... 73 should be enough for up to 256bit floats (with 71.3 decimal digits)
   X_LOG2 = 1.442695040888963407359924681001892137426645954152985934135449406931109219,
   X_INVLOG2 = 0.6931471805599453094172321214581765680755001343602552541206800094933936220;
// information entropy increment in S({p_i}) = \sum_i (-p_i) \log_2(p_i)
static double minus_p_log2_p(double p) {
   if (p <= 0.)
      return 0.;
   return -p*std::log(p) * X_INVLOG2; // <- there is a std::log2, but only since C++11.
}
static double identity_fn(double x) { return x; } // returns input argument unchanged
static double pow05(double x) { return std::sqrt(x); }
static double pow025(double x) { return std::sqrt(std::sqrt(x)); }
static double pow066(double x) { return std::pow(x, 2./3.); }

/*
Table[Subscript[g,d]==FullSimplify[D[h[n],{n,d}]],{h,{Function[n,4*n*(n-1)],Function[n,(4*n*(n-1))^2],Function[n,-n*Log2[n]-(1-n)*Log2[1-n]]}},{d,0,2}]

\text{TexForm}\left(\left(
\begin{array}{ccc}
 g_0=4 (n-1) n & g_1=8 n-4 & g_2=8 \\
 g_0=16 (n-1)^2 n^2 & g_1=32 (n-1) n (2 n-1) & g_2=32 (6 (n-1) n + 1) \\
 g_0=\frac{(n-1) \log (1-n)-n \log (n)}{\log (2)} & g_1=\frac{2 \tanh ^{-1}(1-2
   n)}{\log (2)} & g_2=\frac{1}{(n-1) n \log (2)} \\
\end{array}
\right)\right)
*/
FScalarFn01::~FScalarFn01()
{ // just to bind the destructor and vtable symbols.
}

// FIXME: remove this?
struct FLocFn1_h1 : public FLocFn1 {
   double Eval(double s, double *d1f, double *d2f) const { // override
      double f = 4*s*(s-1);
      if (d1f) *d1f = 8*s - 4;
      if (d2f) *d2f = 8;
      return f;
   }
   std::string Desc() const {
      return "4n(n-1)";
   }
};

// FIXME: remove this?
struct FLocFn1_h2 : public FLocFn1 {
   double Eval(double s, double *d1f, double *d2f) const { // override
      double pref = -1;
      double f = pref*sqr(4*s*(s-1));
      if (d1f) *d1f = pref*32*(s-1)*s*(2*s-1);
      if (d2f) *d2f = pref*32*(1 + 6*(s-1)*s);
      return f;
   }
   std::string Desc() const {
      return "(4n(n-1))^2";
   }
};

struct FLocFn1_Entropy : public FLocFn1 {
   double Eval(double s, double *d1f, double *d2f) const { // override
      double pref = -1;
      double f = pref*(-X_INVLOG2*(s*std::log(s) + (1-s)*std::log(1-s)));
      if (d1f) *d1f = pref*X_INVLOG2*2*std::atanh(1 - 2*s);
      if (d2f) *d2f = pref*X_INVLOG2/(s*(s-1));
      return f;
   }
   std::string Desc() const {
      return "-(n log2(n) + (1-n) log2(1-n))";
   }
};

struct FLocFn1_LogBeta : public FLocFn1 {
   explicit FLocFn1_LogBeta(double shift_of_n=2.)
      : m_shift_of_n(shift_of_n)
   {}

   double Eval(double n, double *d1f, double *d2f) const { // override
#if HAVE_BOOST_SPECIAL_FUNCTIONS
      double pref = X_INVLOG2;
      double f = pref*(boost::math::lgamma(m_shift_of_n + n));
      if (d1f) *d1f = pref*(boost::math::digamma(m_shift_of_n + n));
      if (d2f) *d2f = pref*(boost::math::trigamma(m_shift_of_n + n));
      return f;
#else
      // FIXME: make a polynomial fit of it, if it works? on [0,1] for log2 Gamma(n+2)
      // or log2 Gamma (n+1) that should not be a problem at all...
      throw std::runtime_error("this version of IboView was compiled without support for special functions from boost::math. Please recompile with -DHAVE_BOOST_SPECIAL_FUNCTIONS");
#endif // HAVE_BOOST_SPECIAL_FUNCTIONS
   }
   std::string Desc() const {
      return fmt::format("log2 Gamma({}+n)", m_shift_of_n);
   }
protected:
   double
      m_shift_of_n;
};

struct FLocFn1_nPow : public FLocFn1 {
   explicit FLocFn1_nPow(double Exp=2., double Factor=1., double Shift=0.)
      : m_Exp(Exp), m_Factor(Factor), m_Shift(Shift)
   {}

   double Eval(double n, double *d1f, double *d2f) const { // override
      double pref = m_Factor;
      double f = pref*(std::pow(n + m_Shift, m_Exp));
      if (d1f) *d1f = pref*(m_Exp*std::pow(n + m_Shift, m_Exp-1));
      if (d2f) *d2f = pref*(m_Exp*(m_Exp-1.)*std::pow(n + m_Shift, m_Exp-2.));
      return f;
   }
   double NormalizedL(double L) const { // override
      if (L > 0.)
         return std::pow(L, 1./m_Exp)/m_Factor;
      else
         return L/m_Factor;
   }
   std::string Desc() const {
      std::string sFnArg = "n";
      if (m_Shift < 0)
         sFnArg = fmt::format("(n-{})", -m_Shift);
      else if (m_Shift > 0)
         sFnArg = fmt::format("(n+{})", m_Shift);
      return fmt::format("{}^{}", sFnArg, m_Exp);
   }
protected:
   double
      m_Exp, m_Factor, m_Shift;
};


FScalarFn01Cptr MakeLocFn1_LogBeta(double ShiftOfN) {
   return FScalarFn01Cptr(new FLocFn1_LogBeta(ShiftOfN));
}

FScalarFn01Cptr MakeLocFn1_nPow(double Exp, double Factor, double Shift) {
   return FScalarFn01Cptr(new FLocFn1_nPow(Exp, Factor, Shift));
}






// // symmetrically orthogonalize a set of vectors with respect to overlap matrix S.
// void SymOrth(FMatrixView Orbs, FMatrixView const S, FMemoryStack &Mem)
// {
//    FStackMatrix
//       SmhOrb(Orbs.nCols, Orbs.nCols, &Mem),
//       T1(Orbs.nRows, Orbs.nCols, &Mem);
//    ChainMxm(SmhOrb, Transpose(Orbs), S, Orbs, Mem);
//    CalcSmhMatrix(SmhOrb, Mem, FSmhOptions(1e-15,1e-15,0,0));
//    Mxm(T1, Orbs, SmhOrb);
//    Assign(Orbs, T1);
// }


// def CayleyU(A):
//    # positive Cayley transform: CayU(A) = 1 + A + A^2/2 + A^3/4 + A^4/8 + ...
//    # (that's not the self-inverse one!)
//    assert(A.shape[0] == A.shape[1] and np.allclose(A + A.T, 0.))
//    I = np.eye(A.shape[0])
//    return la.solve(I - .5*A, I + .5*A)
//
// def CayleyA(U):
//    """inverse Caylet transform: Calculate an A which fulfills U == CayleyU(A)."""
//    I = np.eye(U.shape[0])
//    assert(U.shape[0] == U.shape[1] and np.allclose(np.dot(U, U.T) - I, 0.))
//    return 2.*la.solve(I + U, U - I)

// Construct positive Cayley transform: CayU(A) = 1 + A + A^2/2 + A^3/4 + A^4/8 + ...
// (that's not the self-inverse one!)
// Input matrix must be anti-symmetric. Output matrix is allocated on Mem.
FMatrixView MakeCayleyU(FMatrixView A, FMemoryStack &Mem)
{
   assert(A.nRows == A.nCols);
   FMatrixView
      U = MakeStackMatrix(A.nRows, A.nCols, Mem);
   FStackMatrix
      Ip(A.nRows, A.nCols, &Mem),
      Im(A.nRows, A.nCols, &Mem);
   // Ip := I + .5 * A
   Ip.SetIdentity();
   Add(Ip, A, 0.5);
   // Im := I - .5 * A
   Im.SetIdentity();
   Add(Im, A, -0.5);

   // U := la.solve(I - .5*A, I + .5*A)
   LeastSquaresSolveSafe(Im, U, Ip, 1e-12, Mem);
   // ^- LU solve would be enough... but apparently I don't have a primitive for this atm.
   return U;
}


// construct a random N x N anti-symmetric matrix, in which angles are approximately
// distributed according to a normal distribution with standard deviation fAngleStdDev
// (given in radians)
FMatrixView MakeRandomAntiSymmetricMatrix(size_t N, double fAngleStdDev, FMemoryStack &Mem)
{
   FMatrixView
      A = MakeStackMatrix(N, N, Mem);
   FRandomNumberGenerator
      rng;
   for (size_t i = 0; i != N; ++ i) {
      A(i,i) = 0.;
      for (size_t j = 0; j < i; ++ j) {
         double f = fAngleStdDev * rng.GetStdNormal();
         A(i,j) = f;
         A(j,i) = -f;
      }
   }
   return A;
}


// apply a random rotation (unitary transformation) to the vectors in the rows
// of CVec, with the perturbing angles being approximately normal-distributed
// around 0 with standard deviation fAngleStdDev
void RotateVectorsRandomly(FMatrixView CVec, double fAngleStdDev, FMemoryStack &Mem)
{
   TMemoryLock<double>
      pFreeMe0(0, &Mem);
   size_t
      nCols = CVec.nCols;
   FMatrixView
//       A = MakeRandomAntiSymmetricMatrix(nCols, 0.1, Mem),
      A = MakeRandomAntiSymmetricMatrix(nCols, fAngleStdDev, Mem),
      U = MakeCayleyU(A, Mem);
   FStackMatrix
      OriginalCVec(CVec.nRows, CVec.nCols, &Mem);
   Assign(OriginalCVec, CVec);
   Mxm(CVec, OriginalCVec, U);
}







// execute a 2x2 rotation of the vectors pA and pB (both with length nSize)
// of angle phi. In-place.
void Rot2x2(double *RESTRICT pA, double *RESTRICT pB, size_t nSize, double phi)
{
   double
      cs = std::cos(phi),
      ss = std::sin(phi);
   for (std::size_t i = 0; i < nSize; ++ i) {
      double
         tempA =  cs * pA[i] + ss * pB[i],
         tempB = -ss * pA[i] + cs * pB[i];
      pA[i] = tempA;
      pB[i] = tempB;
   }
}


// it's a uglyness-hiding function for interfacing with the localization functionals
// and FLocFn1. Some of them (e.g., n^2, n^4) currently have no generic implementation,
// for no particularly good reason. Others (e.g., entropy) might need some special
// handling in the construction of the Aij/Bij accumulators (singularities etc).
//
// This thing should probably eventually be fixed by a better generic FLocFn1
// implementation. But for the moment this is good enough.
//
// ...So:
// TODO: this entire construction is quite ugly. We should implement
// also n2 and n4 with generic methods, and instead of having the Eval
// function as-is, split it into EvalD0 and EvalD12 (because for the functional
// we need only the value itself, and for the gradient/hessian only the first
// and second derivatives)
struct FLocFn1Proxy {
protected:
   FLocalizeOptions::FLocalizeFn
      m_LocFn;
   FLocFn1 const
      *m_pLocFn1;
   FRealFn
      // function for accumulating contributions to the localization functional
      // and for normalizing the final sum of FnAcc(n(A,i)) contributions.
      m_pLocFnAcc, m_pLocFnNormalize;
public:
   explicit FLocFn1Proxy(FLocalizeOptions const &Opt) {
      m_LocFn = Opt.LocFn;
      m_pLocFn1 = 0;
      m_pLocFnAcc = 0;
      m_pLocFnNormalize = 0;

      switch(m_LocFn) {
         case FLocalizeOptions::LOCFN_nA_pow_2:
            m_pLocFnAcc = pow2; m_pLocFnNormalize = pow05; break;
         case FLocalizeOptions::LOCFN_nA_pow_4:
            m_pLocFnAcc = pow4; m_pLocFnNormalize = pow025; break;
         case FLocalizeOptions::LOCFN_nA_pow_15:
            m_pLocFnAcc = pow15; m_pLocFnNormalize = pow066; break;
         case FLocalizeOptions::LOCFN_nA_Entropy:
            m_pLocFnAcc = minus_p_log2_p; m_pLocFnNormalize = identity_fn; break;
            // ^- warning: this accumulates the *actual* entropy to what is shown
            //    as L (which is MINIMIZED), rather than (-1) times the entropy
            //    which is what we actually optimize (because this routine always
            //    MAXIMIZES). I.e., this relies on pLocFnAcc being purey cosmetic.
         case FLocalizeOptions::LOCFN_nA_gen_h:
            m_pLocFn1 = Opt.pScalarFn01.get();
            if (m_pLocFn1 == 0)
               throw std::runtime_error("LocalizeVectors: localization functional was set to LOCFN_nA_gen_h; but the entry FLocalizeOptions::pScalarFn01, which identifies the general functional h(n), is null.");
            break;
         default:
            throw std::runtime_error(fmt::format("LocalizeVectors: localization functional '{}' not recognized.", int(Opt.LocFn)));
      }
   }

   double NormalizedL(double L) {
      // if specified, apply a normalization function to the final L after all
      // contributions are accumulated. This is meant to simplify comparisons
      // between localization results from different LocFns, and is purely
      // cosmetic.
      if (m_pLocFn1)
         return m_pLocFn1->NormalizedL(L);
      else
         return m_pLocFnNormalize(L); // e.g., L = std::pow(L, 1./double(Opt.nLocExp));
   }

   void AccL(double &L, double nA) {
      if (m_pLocFn1)
         L += m_pLocFn1->Eval(nA,0,0);
      else
         L += m_pLocFnAcc(nA);   // e.g., L += std::pow(nA, (int)Opt.nLocExp);
   }

   void AccAijBij(double &Aij, double &Bij, double QiiA, double QijA, double QjjA) {
      if (m_pLocFn1 != 0) {
         // generic version
         double
            h_Qii, d1h_Qii, d2h_Qii,
            h_Qjj, d1h_Qjj, d2h_Qjj;
         h_Qii = m_pLocFn1->Eval(QiiA, &d1h_Qii, &d2h_Qii);
         h_Qjj = m_pLocFn1->Eval(QjjA, &d1h_Qjj, &d2h_Qjj);
         IR_SUPPRESS_UNUSED_WARNING(h_Qii);
         IR_SUPPRESS_UNUSED_WARNING(h_Qjj);
         Bij += 0.5 * QijA * (d1h_Qii - d1h_Qjj);
         Aij += 0.25 * sqr(QijA)*(d2h_Qii + d2h_Qjj) - 0.125*(QiiA - QjjA)*(d1h_Qii - d1h_Qjj);
      } else if (m_LocFn == FLocalizeOptions::LOCFN_nA_pow_2) {
         Aij += sqr(QijA) - .25*sqr(QiiA - QjjA);
         Bij += QijA*(QiiA - QjjA);
      } else if (m_LocFn == FLocalizeOptions::LOCFN_nA_pow_4) {
         Aij += -pow4(QiiA) - pow4(QjjA) + 6.*(pow2(QiiA) + pow2(QjjA))*pow2(QijA) + pow3(QiiA)*QjjA + QiiA*pow3(QjjA);
         Bij += 4.*QijA*(pow3(QiiA) - pow3(QjjA));
      } else if (m_LocFn == FLocalizeOptions::LOCFN_nA_pow_15) {
         // ...this also doesn't do the benzene as I wanted.
         // I guess I should just implement the general h and try with transfer functions
         // like in the BO paper. I think I possibly *should* treat n~0 and n~1 equally,
         // but none of the functions here do. But I really wonder. Shouldn't the symmetry
         // in entropy, for example, come from the other sum terms?!
         Aij += (-3/16.)*(QiiA - QjjA)*(pow05(QiiA) - pow05(QjjA)) + (3/16.)*pow2(QijA)*(1./pow05(QiiA) + 1./pow05(QjjA));
         Bij += (3./4.)*QijA*(pow05(QiiA) - pow05(QjjA));
      } else if (m_LocFn == FLocalizeOptions::LOCFN_nA_gen_h) {
         throw std::runtime_error("LocalizeVectors(): internal programming error; LOCFN_nA_gen_h set, but pLocFn1 == 0.");
      } else if (m_LocFn == FLocalizeOptions::LOCFN_nA_Entropy) {
//          double ThrQd = 1e-2*Opt.ThrLoc;
         if (0) {
         // if (sqr(QiiA) < sqr(ThrQd) || sqr(QjjA) < sqr(ThrQd)) {
            // if either Qii == 0 or Qjj == 0, the expressions diverge.
            // They diverge in such a way that for a single A,
            // Bij/Aij -> 0. However, there is still the F sum...
         } else {
            // well... it works for acrylic acid and benzene. But not terribly
            // well, and it just gets the qualitatively same localization as output
            // as h(n)=n^4 does...
            double Log2_Qii_minus_Qjj = X_INVLOG2 * (std::log(QiiA) - std::log(QjjA)); // phrase as log(QiiA/QjjA)?
            Aij += .25*sqr(QijA)/(1./QiiA + 1./QjjA) - 0.125*(QiiA - QjjA)*Log2_Qii_minus_Qjj;
            Bij += .5*Log2_Qii_minus_Qjj*QijA;
         }
      } else {
         throw std::runtime_error("LocalizeVectors(): internal programming error; unrecognized update branch in localization functional.");
      }
   }
};


// Find an unitary transformation of the input vectors such that on output,
// the vectors maximize a localization functional
//
//    L = \sum_{i}\sum_{A} LocFn(n(A,i)) -> max
//
// where the fraction of orbital $i$ assigned to group $A$ is given by
//
//    n(A,i) := \sum_{\mu \in A} CVec[\mu,i]^2.
//
// Arguments:
//    - CVec: (nFn,nVec)-shape matrix containing the input vectors on entry
//      and the unitarily transformed output vectors on output.
//
//    - pGroupOffsets: (nGroups+1)-length array of integers identifying the
//      boundaries between different fragments/groups (in orbital localization,
//      these will typically be the boundaries of orthonormal basis functions
//      associated with either individual atoms or groups of atoms).
//
//      Concretely, the vector components CVec[\mu,i] for \mu fulfilling
//      pGroupOffsets[A] <= \mu < pGroupOffsets[A+1]
//      are assiciated with group A.
//
//    - nGroups: number of groups/fragments in pGroupOffsets.
//
// Notes:
// - This function could be extended to localize on other criteria
//   than charge, provided a replacement of the calculation of the Cii/Cjj/Cij
//   matrix elements is provided.
//
//   Most useful in the local exchange case is likely the optimization of <i|(r-A)^4|i>.
//   It can also be used as a replacement for ZBD orthogonalization then (which is
//   equivalent to optimizing on (r-A)^2).
//
// References:
// - [1] Gerald Knizia; "Intrinsic atomic orbitals: An unbiased bridge between quantum
//       theory and chemical concepts;"
//       JCTC 2013, 9, 4834-4843, https://doi.org/10.1021/ct400687b
//
//       + This is the original reference for the IAO and IBO methods.
//
//       + The localization algorithm which this function executes is described
//         in the appendix for the special cases LocFn(n) = n^p p\in{2,3,4}
//
//       + It is somewhat terse regarding technical details. For full details,
//         see [2] below.
//
// - [2] Bruno Senjean, Souloke Sen, Michal Repisky, Gerald Knizia, Lucas Visscher;
//       "Generalization of intrinsic orbitals to Kramers-paired quaternion spinors,
//       molecular fragments and valence virtual spinors"
//       (submitted to JCTC; preprint: https://arxiv.org/abs/2009.08671)
//
//       + Apart from introducing relativistic generalizations of the IAO and
//         IBO methods, this paper also contains full derivations and technical
//         explanations of both the IAO construction itself and the localization
//         procedure executed here. These can be found in its Supporting Information,
//         and they also cover the non-relativistic case.
//
//       + This is the paper which derived the formulas for optimizing on general
//         1-argument localization functions LocFn and introduced the minimum entropy
//         localization. Additionally, it introduced some tricks regarding convergence
//         stabilization (all in the SI).
void LocalizeVectors(FMatrixView CVec, size_t const *pGroupOffsets, size_t nGroups, FLocalizeOptions const &Opt, FMemoryStack &Mem)
{
   IR_SUPPRESS_UNUSED_WARNING(Mem);

   if (1) {
//       double fAngleDeg = 180.;
      double fAngleDeg = 18.;
//       double fAngleDeg = 45.;
      if (Opt.pLog) {
         Opt.pLog->Write(" Apply random perturbation to input vectors: {:.2f} deg.", fAngleDeg);
      }
      RotateVectorsRandomly(CVec, (M_PI/180.) * fAngleDeg, Mem);
   }


// #ifdef INCLUDE_OPTIONALS
//    FRandomNumberGenerator
//       rng;
// #endif // INCLUDE_OPTIONALS
   // compute inverse vector norms, in case there are absorbed occupation
   // numbers or other types of non-standard normalization. The pow2 and pow4
   // functionals are not terribly interested in that, but the entropy
   // functional will not work correctly unless vectors are normalized to 1.
   //
   // FIXME: no idea what I was thinking. This cannot possibly work unless
   // all vector norms are equal! The vectors are mixed with each other below...
   TMemoryLock<double>
      pVecReNorm(CVec.nCols, &Mem);
   for (size_t iVec = 0; iVec < CVec.nCols; ++ iVec) {
      double r = 0;
      for (size_t iFn = 0; iFn != CVec.nRows; ++ iFn)
         r += pow2(CVec(iFn,iVec));
      pVecReNorm[iVec] = 1./std::sqrt(r);
   }

   FLocFn1Proxy
      LocFn1Proxy(Opt);
//       LocFn1Proxy(Opt.LocFn, Opt.pScalarFn01);

   size_t
      iIt;
   double
      L = -1., fVar2 = -1.;
   if (Opt.pLog && Opt.Verbosity >= 2) {
      Opt.pLog->Write("   ITER.       LOC(Orb)      GRADIENT");
   }
   for (iIt = 0; iIt < Opt.nMaxIt; ++ iIt) {
      // calculate the value of the functional.
      L = 0.;
      // loop over groups (normally a ``group'' is all functions on an atom)
      for (size_t iGrp = 0; iGrp < nGroups; ++ iGrp)
         for (size_t iVec = 0; iVec < CVec.nCols; ++ iVec) {
            double nA = 0.; // population on the atom.
            for (size_t iFn = pGroupOffsets[iGrp]; iFn != pGroupOffsets[iGrp+1]; ++ iFn)
               nA += pow2(CVec(iFn,iVec));
            nA *= pow2(pVecReNorm[iVec]);
            LocFn1Proxy.AccL(L, nA); // e.g., L += std::pow(nA, (int)Opt.nLocExp);
         }
      // apply a (purely cosmetic) normalization of the value of the functional
      // to simplify comparision between different localizations. Note: not
      // possible/sensible for all driving functions.
      L = LocFn1Proxy.NormalizedL(L);  // e.g., L = std::pow(L, 1./double(Opt.nLocExp));

      fVar2 = 0.;
      // loop over vector pairs.
      for (size_t iVec = 0; iVec < CVec.nCols; ++ iVec)
//          for (unsigned nMic = 0; nMic < 5; ++ nMic)
         for (size_t jVec = 0; jVec < iVec; ++ jVec)
//          for (unsigned nMic = 0; nMic < 3; ++ nMic)
         {
//             if (iIt > 50 && rng.GetUniformF(0.,1.) < 0.5)
//                continue;
            double
               Aij = 0.0, // hessian for 2x2 rotation (with variable atan(phi/4) substituted)
               Bij = 0.0; // gradient for 2x2 rotation (with variable atan(phi/4) substituted)
            // loop over groups (normally atoms)
            for (size_t iGrp = 0; iGrp < nGroups; ++ iGrp) {
               // calculate the charge matrix elements
               //     QiiA = <i|P_A|i>,
               //     QjjA = <j|P_A|j>,
               // and QijA = <i|P_A|j>.
               double
                  QiiA = 0., QijA = 0., QjjA = 0.;
               for (size_t iFn = pGroupOffsets[iGrp]; iFn != pGroupOffsets[iGrp+1]; ++ iFn) {
                  QiiA += CVec(iFn,iVec) * CVec(iFn,iVec);
                  QjjA += CVec(iFn,jVec) * CVec(iFn,jVec);
                  QijA += CVec(iFn,iVec) * CVec(iFn,jVec);
               }
               QiiA *= pVecReNorm[iVec] * pVecReNorm[iVec];
               QijA *= pVecReNorm[iVec] * pVecReNorm[jVec];
               QjjA *= pVecReNorm[jVec] * pVecReNorm[jVec];

               LocFn1Proxy.AccAijBij(Aij, Bij, QiiA, QijA, QjjA);
            }
            // filter out rotations which have a very small Hessian and gradient
            // at the same time. In PM-like optimizations this happens, for example,
            // if there are multiple lone pairs or on the same center.
//             if (std::abs(Aij) > Opt.ThrLoc && pow2(Bij)/std::abs(Aij) > 1e-2*Opt.ThrLoc) {
               // ^- note: Bij * (Bij/Aij) is an estimate for the change in L:
               //    delta_x = Bij/Aij, and this is multiplied by dL/dx = Bij.
               //    this worked for all variants in FeNO(CO)3, but it makes
               //    the entropy functional fail completely in benzene.
            if (std::abs(Aij) > Opt.ThrLoc) {
               double
                  phi = .25 * std::atan2(Bij,-Aij);
               // FIXME: this does not update the vector norms!!
               // 2x2 rotate vectors iVec and jVec.
               Rot2x2(&CVec(0,iVec), &CVec(0,jVec), CVec.nRows, phi);
//                fVar2 += pow2(phi);
               fVar2 += pow2(4*phi*Bij); // angle times angle D[L,x]; note near x=0, 4 phi\approx x.
//                fVar2 += pow2(4*phi*Bij) + pow4(phi); // 4*phi*Bij: angle times angle D[L,x]; note near x=0, 4 phi\approx x.
            }


#ifdef INCLUDE_OPTIONALS
//             // non-degenerate rotation? This condition is not quite right,
//             // and in Python it was not required. However, the Fortran atan2
//             // sometimes did mysterious things without it.
//             //
//             // note: for the L_2 functional, the second derivatives at the
//             // maxima come out as -16 \sqrt{A_{ij}^2 + B_{ij}^2}
//             // UPDATE: trying to re-localize examples/feco3no_ibos_exp2.xml
//             // required ungodly amounts of iterations. Why? Very small
//             // angles combined with very small Hessians, maybe?
// //             if (pow2(Aij) > pow2(Opt.ThrLoc * 1e-2)) {
// //             if (pow2(Aij) + pow2(Bij) > pow2(Opt.ThrLoc * 1e-2)) {
// //             if (pow2(Aij) + pow2(Bij) > Opt.ThrLoc * 1e-2) {
// //             if (1) {
// //             if (std::abs(Aij) > Opt.ThrLoc) {
// //             if (pow2(Aij) > Opt.ThrLoc && std::abs(Bij) > Opt.ThrLoc) {
//             if (pow2(Aij) > Opt.ThrLoc && std::abs(Bij) > Opt.ThrLoc) {
//                double
//                   phi = .25 * std::atan2(Bij,-Aij);
//                // 2x2 rotate vectors iVec and jVec.
//                Rot2x2(&CVec(0,iVec), &CVec(0,jVec), CVec.nRows, phi);
//                fVar2 += pow2(phi);
// //                if (iIt >= 200 && iIt <= 210 && pow2(phi) > 1e-12)
//                if (iIt >= 200 && iIt <= 210 && std::abs(phi) > 1e-10)
// //                if (iIt >= 200 && iIt <= 210)
//                   Opt.pLog->Write("       IT:{:6}    iOrb+1={:<4}  jOrb+1={:<4}  // Aij={:9.2e}   Bij={:9.2e}  --> phi={:8.2e}   THRLOC: {:8.2e}", (1+iIt), (1+iVec), (1+jVec), Aij, Bij, phi, Opt.ThrLoc);
//             }
#endif // INCLUDE_OPTIONALS
         }
      fVar2 = std::sqrt(fVar2 / pow2(CVec.nCols));

      if (Opt.pLog && Opt.Verbosity >= 2) {
         Opt.pLog->Write("{:6}     {:12.6f}      {:8.2e}", (1+iIt), L, fVar2);
      }

      if (fVar2 < Opt.ThrLoc)
         break;
   }

   if (Opt.pLog && Opt.Verbosity >= 1) {
      if (Opt.Verbosity >= 2)
         Opt.pLog->WriteLine();
//       Opt.pLog->Write(" Iterative localization: IB/PM, {} iter; Final gradient {:8.2e}{}", iIt, fVar2, Opt.TaskDesc.empty()? "" : ("; " + Opt.TaskDesc));
      Opt.pLog->Write(" Iterative localization: IB/PM, {} iter -> L = {:.6f}; Final gradient {:8.2e}{}", iIt, L, fVar2, Opt.TaskDesc.empty()? "" : ("; " + Opt.TaskDesc));
      Opt.pLog->Flush();
   }
}




































//
//
// // NOTE: There is a working version of this in IvIao.cpp in IboView. Should be merged.
//
// namespace ct {
//
// struct FLocalizeOptions
// {
//    uint
//       nLocExp; //  2 or 4, exponent p in localization functional: L = sum[iA] n(i,A)^p
//    uint
//       nMaxIt;
//    double
//       ThrLoc;
//    std::ostream
//       *pxout; // if non-0, write output here.
//    uint
//       Verbosity;
//
//    FLocalizeOptions()
//       : nLocExp(4), nMaxIt(2048), ThrLoc(1e-8), pxout(0), Verbosity(1)
//    {}
// };
//
// inline double pow2(double x) { return x*x; }
// inline double pow3(double x) { return x*x*x; }
// inline double pow4(double x) { double x2 = x*x; return x2*x2; }
//
// // execute a 2x2 rotation of the vectors pA and pB (both with length nSize)
// // of angle phi. In-place.
// void Rot2x2(double *RESTRICT pA, double *RESTRICT pB, size_t nSize, double phi)
// {
//    double
//       cs = std::cos(phi),
//       ss = std::sin(phi);
//    for (std::size_t i = 0; i < nSize; ++ i) {
//       double
//          tempA =  cs * pA[i] + ss * pB[i],
//          tempB = -ss * pA[i] + cs * pB[i];
//       pA[i] = tempA;
//       pB[i] = tempB;
//    }
// }
//
// // localize a set of vectors onto
// //  L = \sum[i,A] \sum[mu \in A] CVec(mu,i)^p -> min.
// // in a PM-like fashion.
// // Parameters:
// //    - pGroupOffsets: Array of length nGroups+1.
// //      [pGroupOffsets[A], pGroupOffsets[A+1]) gives the association of
// //      vector components to groups.
// // Note: This function can be extended to localize on other criteria
// // than charge, provided a replacement of the calculation of the Cii/Cjj/Cij
// // matrix elements is provided. Most useful in the local exchange case is likely
// // the optimization of <i|(r-A)^4|i>. It can also be used as a replacement for
// // ZBD orthogonalization then (which is equivalent to optimizing on (r-A)^2).
// void LocalizeVectors(FMatrixView CVec, size_t *pGroupOffsets, size_t nGroups, FLocalizeOptions const &Opt, FMemoryStack2 const &Mem)
// {
//    IR_SUPPRESS_UNUSED_WARNING(Mem);
//
//    size_t
//       iIt;
//    double
//       L = -1., fVar2 = -1.;
//    if (Opt.pxout && Opt.Verbosity >= 2) {
//       *Opt.pxout << "  ITER.   LOC(Orb)      GRADIENT" << std::endl;
//    }
//    for (iIt = 0; iIt < Opt.nMaxIt; ++ iIt) {
//       // calculate the value of the functional.
//       L = 0.;
//       // loop over groups (normally a ``group'' is all functions on an atom)
//       for (size_t iGrp = 0; iGrp < nGroups; ++ iGrp)
//          for (size_t iVec = 0; iVec < CVec.nCols; ++ iVec) {
//             double nA = 0.; // population on the atom.
//             for (size_t iFn = pGroupOffsets[iGrp]; iFn != pGroupOffsets[iGrp+1]; ++ iFn)
//                nA += pow2(CVec(iFn,iVec));
//             L += std::pow(nA, Opt.nLocExp);
//          }
//       L = std::pow(L, 1./Opt.nLocExp);  // <- easier to compare different powers that way.
//
//       fVar2 = 0.;
//       // loop over vector pairs.
//       for (size_t iVec = 0; iVec < CVec.nCols; ++ iVec)
//          for (size_t jVec = 0; jVec < iVec; ++ jVec) {
//             double
//                Aij = 0.0, // hessian for 2x2 rotation (with variable atan(phi/4) substituted)
//                Bij = 0.0; // gradient for 2x2 rotation (with variable atan(phi/4) substituted)
//             // loop over groups (normally atoms)
//             for (size_t iGrp = 0; iGrp < nGroups; ++ iGrp) {
//                // calculate the charge matrix elements
//                //     Cii = <i|A|i>,
//                //     Cjj = <j|A|j>,
//                // and Cij = <i|A|j>.
//                double
//                   Cii = 0., Cij = 0., Cjj = 0.;
//                for (size_t iFn = pGroupOffsets[iGrp]; iFn != pGroupOffsets[iGrp+1]; ++ iFn) {
//                   Cii += CVec(iFn,iVec) * CVec(iFn,iVec);
//                   Cjj += CVec(iFn,jVec) * CVec(iFn,jVec);
//                   Cij += CVec(iFn,iVec) * CVec(iFn,jVec);
//                }
//
//                // see Supporting Information for `` Intrinsic atomic orbitals:
//                // An unbiased bridge between quantum theory and chemical concepts''
//                // for a description of what this does.
//                if (Opt.nLocExp == 2 || Opt.nLocExp == 3) {
//                   Aij += pow2(Cij) - .25*pow2(Cii - Cjj);
//                   Bij += Cij*(Cii - Cjj);
//                } else if (Opt.nLocExp == 4) {
//                   Aij += -pow4(Cii) - pow4(Cjj) + 6.*(pow2(Cii) + pow2(Cjj))*pow2(Cij) + pow3(Cii)*Cjj + Cii*pow3(Cjj);
//                   Bij += 4.*Cij*(pow3(Cii) - pow3(Cjj));
//                }
//             }
//
//             // non-degenerate rotation? This condition is not quite right,
//             // and in Python it was not required. However, the Fortran atan2
//             // sometimes did mysterious things without it.
//             if (pow2(Aij) + pow2(Bij) > Opt.ThrLoc * 1e-2) {
//                double
//                   phi = .25 * std::atan2(Bij,-Aij);
//                // 2x2 rotate vectors iVec and jVec.
//                Rot2x2(&CVec(0,iVec), &CVec(0,jVec), CVec.nRows, phi);
//                fVar2 += pow2(phi);
//             }
//          }
//       fVar2 = std::sqrt(fVar2 / pow2(CVec.nCols));
//
//       if (Opt.pxout && Opt.Verbosity >= 2) {
//          *Opt.pxout << fmt::format("{:>5}     {:12.6f}      {:8.2e}\n", (1+iIt), L, fVar2);
//       }
//
//       if (fVar2 < Opt.ThrLoc)
//          break;
//    }
//
//    if (Opt.pxout && Opt.Verbosity >= 1) {
//       if (Opt.Verbosity >= 2)
//          *Opt.pxout << "\n";
//       *Opt.pxout << fmt::format("Iterative localization: IB/PM, {} iter; Final gradient {:8.2e}", iIt, fVar2) << std::endl;
//    }
// }
//
//
//
// void LocalizeOrbitals(FBasisSetPtr pBasis, FMatrixView COrb, FLocalizeOptions const &Opt);
//
//
// } // namespace ct





} // namespace ct
