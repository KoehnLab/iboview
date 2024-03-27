// clear && c++ -fmax-errors=3 -std=c++11 -Wpedantic -Wall -Wpedantic -g ibeta.cpp && a.out
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip> // for std::setw & co
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>


/// Auxiliary class to evaluate continued fractions
///
///   CfFwd[L] := a[0]/(b[0] + a[1]/(b[1] + a[2]/(b[2] + a[3]/(b[3] + a[4]/( .... + a[N-1]/b[N-1])))))
///
/// with successive coefficient pairs (a_i,b_i) provided incrementally by a
/// generator object.
///
/// Details:
/// - Incrementally build the full result of ContFracFwd[L] as
/// 
///      (CfFwd(L[0:2])/CfFwd(L[0:1])) * (CfFwd(L[0:3])/CfFwd(L[0:2])) * (CfFwd(L[0:4])/CfFwd(L[0:3])) * ...
/// 
///    by computing the fractions in this expressions via the relation
///
///      cfFwd0/cfFwd1 == cfRev1/cfRev0
///
///   from the reverse continued fraction transformation. (Unlike for the direct
///   forward computation, the reverse variant allows incremental evaluation
///   starting at earlier coefficient list elements and proceeding to later
///   ones. So we do not need to know a-priorily how many terms we require)
///
/// - CfFwd[L] stands for ContFracFwd[L], and CfRev[L] for ContFracRevN[L].
///
///   Concrete definitions of continued fraction functions (Mathematica syntax;
///   for details see continued_fractions_algo.nb):
///  
///      ContFracFwd[L_] := 
///         If [Length[L] == 0, 0,
///            Block[{am,bm},
///               {am, bm} = First[L];                (* get pair in list element 1*)
///               am / (bm + ContFracFwd[Rest[L]])    (* apply trafo and recurse on elements {2,3,...,N} *)
///            ]
///         ]
///  
///      ContFracRevN[L_] := 
///         If [Length[L] == 0, Infinity,
///            Block[{am, bm},
///               {am, bm} = Last[L];                 (* get pair in list element N *)
///               bm + am / ContFracRevN[Drop[L,-1]]  (* apply trafo and recurse on elements {1,2,...,N-1} *)
///            ]
///         ]
///
/// - In comments the list slices are described in python syntax (L[ibeg:iend],
///   with iend exclusive and ibeg starting at 0).
template<class FScalar>
struct TContFrac {
   
   /// compute continued fraction
   ///
   /// CfFwd[L] := a[0]/(b[0] + a[1]/(b[1] + a[2]/(b[2] + a[3]/(b[3] + a[4]/( .... + a[N-1]/b[N-1])))))
   ///
   template<class FCoeffGenerator>
   FScalar Eval(FCoeffGenerator *pGen)
   {
      bool
         Converged = false;
      FScalar
         // intermediate multiplicative accumulator for target output
         // ContFracFwd[L] value; concretely, it is
         // `OutputConst * CfFwd[L[0:k]]`
         // (global factors in this will be propagated)
         cfFwdK,
         // intermediate accumulator for CfRev[L[0:k]], but in the 1/x version.
         cfRev0InverseK,
         // intermediate accumulator for CfRev[L[1:k]], but in the 1/x version.
         cfRev1InverseK;
      {
         FScalar a1, b1, b1Inv;
         // get first set of continued fraction coefficients a_1, and b_1.
         pGen->MakeNextCoeffs(&a1, &b1);
         b1Inv = 1/b1;
         // set up initial values. (note: the incremental ratio formula only
         // works for k >= 1)
         cfFwdK = a1*b1Inv;
         cfRev0InverseK = b1Inv;
         cfRev1InverseK = 0;
      }
      // now iteratively compute CfFwd[L] from the inputs...
      size_t k = 1;
      for (; k < 4096; ++k) {
         FScalar
            // continued fraction coefficients a_k and b_k for k = current cf depth
            ak, bk,
            // the current depth's intermediate results CfRev(L[0:k]) and
            // CfRev(L[1:k]) for the intermediate negated reverse continued
            // fractions used to compute incremental ratios of cfFwdK/cfFwdKm1.
            cfRev0K, cfRev1K;

         // collect continued fraction coefficients (a_k, b_k) for the current level.
         pGen->MakeNextCoeffs(&ak, &bk);

         // At this moment, cfRevXInverseK are still 1/cfRevX_Km1; that is, for the last
         // iteration's `k`, which is now k-1.
         // So the next two expressions are requivalent to:
         //
         //      cfRev0K = b_k + a_k / cfRev0Km1  // <- compute CfRev(L[0:k]) from CfRev(L[0:(k-1)])
         //      cfRev1K = b_k + a_k / cfRev1Km1  // <- compute CfRev(L[1:k]) from CfRev(L[1:(k-1)])
         //
         // however, we formulate these two relations partially in terms of
         // their 1/x quantities to simplify dealing with the starting value
         // cfRev1 == Infinity and to reduce the number of divisions.
         cfRev0K = bk + ak * cfRev0InverseK;
         cfRev1K = bk + ak * cfRev1InverseK;
         // update the (1/x) inverse quantities (rev0 for ratio, and both for next iteration)
         cfRev0InverseK = 1/cfRev0K;
         cfRev1InverseK = 1/cfRev1K;
         
         // for k < len(L) and len(L) >= 2, there is a continued fraction identity
         //
         //    cfFwdK/cfFwdKm1 == cfRev1K/cfRev0K
         //
         // relating the ratio of successive forward continued fractions:
         //
         //    cfFwdK/cfFwdKm1 = CfFwd(L[0:k]) / CfFwd(L[0:(k-1)])
         //
         // to ratio of the negated reverse continued fractions:
         //
         //    cfRev1K/cfRev0K = CfRev(L[1:k]) / CfRev(L[0:k]).
         //
         // Now use this to evaluate ratio := cfFwdK/cfFwdKm1 via cfRev1K/cfRev0K.
         // (for details, see continued_fractions_algo.nb)
         FScalar
            cfFwdK_div_cfFwdKm1 = cfRev1K * cfRev0InverseK;

         // compute CfFwd(L[0:k]) from cfFwdKm1 = CfFwd(L[0:(k-1)]) and the ratio
         // cfFwd0/cfFwd1, which we just evaluated via cfFwd0/cfFwd1 == cfRev1/cfRev0.
         // (at this moment cfFwdK is cfFwdKm1)
         cfFwdK *= cfFwdK_div_cfFwdKm1;
         
         {
            using std::abs;
            if (abs(FScalar(1) - cfFwdK_div_cfFwdKm1) <= 4*std::numeric_limits<FScalar>::epsilon()) {
               Converged = true;
               std::cout << "      ...Cf::Eval: stopped at k = " << k << " and ratio " << cfFwdK_div_cfFwdKm1 << " (last ak = " << ak << ", last bk = " << bk << ", " << "last cfRev0K = " << cfRev0K << ", " << "last cfRev1K = " << cfRev1K << ")" << std::endl;
               break;
            }
         }
         // note: Current k will be (k-1) in next iteration; so in the next iteration,
         // cfFwdK, cfRev0InverseK, and cfRev1InverseK  will really start out as
         // cfFwdKm1, cfRev0InverseKm1, and cfRev1InverseKm1
      }
      if (!Converged) {
         std::stringstream str;
         str << "stopped at k = " << k << std::setprecision(16) << std::scientific << " (" << "last 1/cfRev0K = " << cfRev0InverseK << ", " << "last 1/cfRev1K = " << cfRev1InverseK << ", " << "(ratio-1) = " << (cfRev0InverseK/cfRev1InverseK - 1) << ")";
         throw std::runtime_error("TContFrac::Eval() failed to converge: " + str.str());
      }
      return cfFwdK;
   }
};



template<class FScalar, class FInteger = std::ptrdiff_t>
struct TBetaContFracGenerator
{
   TBetaContFracGenerator(FScalar alpha_, FScalar beta_, FScalar x_)
   : m_NextK(0), alpha(alpha_), beta(beta_), x(x_)
   {}

   void MakeNextCoeffs(FScalar *pak, FScalar *pbk) {
      FInteger
         k = m_NextK;
      ++ m_NextK;
   
      // see: https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/math_toolkit/sf_beta/ibeta_function.html
      // note: we start at k=0, not at m=1.
      // ...but the formulas are given for a_{m+1} and b_{m+1}, so we
      // anyway can replace m by k without further offsets.
      FScalar
         f1d_aP2mM1 = 1/(alpha + FScalar(2*k - 1)),
         f1d_aP2mP1 = 1/(alpha + FScalar(2*k + 1)),
         k_betaMk_x_d = k*(beta - k)*x * f1d_aP2mM1;
      if (k == 0) {
         *pak = 1;
         *pbk = f1d_aP2mP1 * (alpha*(alpha - (alpha+beta)*x + 1));
      } else {
         *pak = (f1d_aP2mM1 * (alpha + k - 1) * (alpha + beta + k - 1)) * (x*k_betaMk_x_d);
         *pbk = k + k_betaMk_x_d + f1d_aP2mP1 * ((alpha + k)*(alpha - (alpha+beta)*x + 1 + k*(2-x)));
      }
   }
private:
   FInteger
      m_NextK = 0;
   FScalar
      alpha, beta, x;
};


/// Implements an non-regularized incomplete beta function via continued fraction
/// expansion. Very straight-forward---no shortcuts, no asymptotic expansions, etc.
///
/// Note:
/// - Integer variants of regularized ibeta (=generalized PolyStep)
///   can be made by
///
///     HornerForm[FunctionExpand[Beta[x,alpha+1,beta+1]],x]
///
///   for concrete alpha/beta. Although the direct variant with the y terms
///   might be a better call. Should check (well, i probably want either beta
///   or betac, but not mix them...).
///   I did implement those in /quadlib/unpublished_experimental/PolyStep.inl
/// - I also had a binomial expansion for those which could be run at runtime.
///   It's in the quadlib writeup, but I did not implement it in code (could
///   probably do it with reasonable amount of work using IrFactorials).
template<class FScalar>
FScalar MyLameBeta(FScalar const &alpha, FScalar const &beta, FScalar const &x)
{
   using std::pow;
   TBetaContFracGenerator<FScalar>
      cfGen(alpha, beta, x);
   TContFrac<FScalar>
      cf;
   FScalar
      y = 1 - x,
      cfVal = cf.Eval(&cfGen),
      xyPow = pow(x, alpha) * pow(y, beta); // todo: check if alpha == beta!
   return xyPow * cfVal;
}


template<class FScalar>
void WriteResult(std::string const &sTitle, FScalar const &val, FScalar const *pRef = 0)
{
   std::ostream &out = std::cout;
   out << std::left << std::setw(40) << sTitle << " ";
   out << std::fixed << std::setw(24) << std::setprecision(15) << val;
   if (pRef) {
      FScalar fErr = (val - *pRef) / *pRef;
      out << "   // d(ref) = " << std::scientific << std::setprecision(2) << fErr;
   }
   out << std::endl;
}


int main()
{
   std::ostream &out = std::cout;
   using std::exp;
   using std::tgamma;
   using std::lgamma;
   typedef double FScalar;
   FScalar
      alpha = 2.5, beta = 1.37, x = 0.37;
      
   FScalar
      ibeta_ref = boost::math::ibeta(alpha, beta, x), // <- that (ibeta) is the regularized one.
      beta_ref = boost::math::beta(alpha, beta, x), // <- that (beta) is the non-regularized one.
      cbeta_ref = boost::math::beta(alpha, beta); // that's the complete beta fn.
   FScalar
      cbeta_via_lgamma = exp(lgamma(alpha) + lgamma(beta) - lgamma(alpha+beta)),
      cbeta_via_tgamma = tgamma(alpha)/tgamma(alpha + beta) * tgamma(beta);
   
//       dt_dx = std::pow(x, alpha-1) * std::pow(1.-x, beta-1) / boost::math::beta(alpha, beta);
   
   WriteResult("α", alpha);
   WriteResult("β", beta);
   WriteResult("x", x);
   out << "\n";
   WriteResult("B(α,β) // boost", cbeta_ref);
   WriteResult("B(α,β) // lgamma", cbeta_via_lgamma, &cbeta_ref);
   WriteResult("B(α,β) // tgamma", cbeta_via_tgamma, &cbeta_ref);
   out << "\n";
   WriteResult("B(α,β,x) // boost", beta_ref);
   FScalar
      beta_test = MyLameBeta(alpha, beta, x);
   WriteResult("B(α,β,x) // mine", beta_test, &beta_ref);
   
   out << "\n";
   FScalar
      // make the regularized variant (note: in our application we can precompute 1/cbeta(alpha,beta),
      // so there is no point in computing it efficiently. If the tgamma/lgamma thing works and is
      // resonably accurate, then that is good enough).
      ibeta_test = MyLameBeta(alpha, beta, x) / cbeta_via_lgamma;
   WriteResult("RegB(α,β,x) // boost", ibeta_ref);
   WriteResult("RegB(α,β,x) // mine", ibeta_test, &ibeta_ref);
}
