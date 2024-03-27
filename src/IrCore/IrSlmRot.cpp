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

#include "Ir.h"
#include "IrAmrr.h"
#include "CxMemoryStack.h"
#include "CxAlgebra.h"
#include <algorithm> // for std::min.
#include <stdexcept> // for std::runtime_error

namespace ct {
   size_t LeastSquaresSolve(double *pSol, size_t ldSol, double const *pA, size_t ldA, double const *pRhs, size_t ldRhs, size_t nRows, size_t nCols, size_t nVec, double ThrSig);
}

namespace ir {

void EvalSlcX_Deriv0(double *IR_RP Out, double x, double y, double z, unsigned L);


/// For a given 3x3 unitary rotation matrix R[3][3], construct the matrix T[lc,l'c'] which makes
///    S[l,c](R r) = \sum_{c',l'} T[lc,l'c'] S[l',c'](r)
/// for all 3-vectors 'r' (i.e., compute the transformation between solid harmonic components
/// induced by rotations in real space).
///
/// Notes:
///  - The matrix elements T[lc,l'c'] for l != l' vanish.
///  - T[lc,l'c'] is stored at pOut[iSlcX(l,c) + nSlmX(MaxL) * iSlcX(l',c')].
///    nStride is set to nSlmX(MaxL); pOut is allocated on Mem by this routine.
///  - Order of trafo not checked. If it doesn't work, use transpose(T) instead of T.
void EvalSlcXRotationTrafo(double *&pOut, size_t &nStride, unsigned MaxL, double const *R, ct::FMemoryStack &Mem)
{
   nStride = nSlmX(MaxL);
   Mem.Alloc(pOut, nStride*nStride);
   ct::TMemoryLock<double>
      pFreeMe(0, &Mem); // <-- auto-free memory up to here
   unsigned
      nSl = nSlmX(MaxL),
      nCa = nCartX(MaxL);
   double
      *pA, *pB;
   Mem.Alloc(pA, nSl * nCa);
   Mem.Alloc(pB, nSl * nCa);

   // while this could probably be done by employing some funky solid harmonic
   // transformation theorems, we here just set up an equation system such that
   //
   //     a_i = T * b_i,
   //
   // where a_i = Slc(R r_i) and b_i = Slc(r_i) for a sufficient basis
   // of vectors r_i (with Slc components in the rows).
   //
   // The equation system is then linearly solved via LAPACK.
   //
   // Such a set of vectors is given, for example, by the cartesian vectors
   //     r_i = (x_i,y_i,z_i)
   // with x_i + y_i + z_i = MaxL (for individual l), or
   // with all x_i + y_i + z_i <= MaxL (for all l at once).
   //
   // UPDATE:
   // - actually, x_i + y_i + z_i in {MaxL, MaxL-1} should be sufficient
   //   for the second case (all l <= MaxL). This would even be a non-redundant
   //   equation system and could be LU'd instead of SVD'd.)
   //
   // UPDATE:
   // - See longer comment below function

   // setup the matrices B = (b_i) and A = (a_i).
   unsigned
      iComp = 0;
   for (unsigned i = 0; i <= MaxL; ++i)
      for (unsigned j = 0; j <= i; ++j)
         for (unsigned k = 0; k <= j; ++k) {
            // make trial vector x and it's rotated counterpart R x.
            double
               r[3] = {double(i), double(j), double(k)},
               Rr[3] = {0., 0., 0.};
            for (unsigned a = 0; a < 3; ++ a)
               for (unsigned b = 0; b < 3; ++ b)
                  Rr[a] += R[a + 3*b] * r[b];
            // evaluate SlcX(x) and SlcX(R x).
            EvalSlcX_Deriv0(&pB[iComp * nSl], r[0], r[1], r[2], MaxL);
            EvalSlcX_Deriv0(&pA[iComp * nSl], Rr[0], Rr[1], Rr[2], MaxL);
            iComp += 1;
         }
   assert(iComp == nCa);

   // solve for T such that A = T B.
   // By transposing both sides we get this into standard form:
   //   A^T = B^T T,
   // where we solve for T. Due to LAPACK restrictions we do the
   // transposition explicitly here.
   double
      *pAt, *pBt;
   Mem.Alloc(pAt, nCa * nSl);
   Mem.Alloc(pBt, nCa * nSl);
   for (unsigned iSl = 0; iSl < nSl; ++ iSl)
      for (unsigned iCa = 0; iCa < nCa; ++ iCa) {
         pAt[iCa + nCa*iSl] = pA[iSl + nSl*iCa];
         pBt[iCa + nCa*iSl] = pB[iSl + nSl*iCa];
      }

   // FIXME: do at least the SVDs individually per l. Results should be the same,
   // but this would make it faster.
//    double
//       *pSig, *pWork;
//    Mem.Alloc(pSig, nCa);
//    FORTINT
//       nRank = -1,
//       info = 0,
//       // minimal work space required according to docs.
//       nWork = 3*std::min(nCa,nSl) + std::max(std::max(2*std::min(nCa, nSl), std::max(nCa, nSl)), nSl);
//    nWork += 2*nCa*nSl*nSl; // <- we don't really care...
//    Mem.Alloc(pWork, nWork);
//    DGELSS(nCa, nSl, nSl, pAt, nCa, pBt, nCa, pSig, 1e-10, &nRank, pWork, nWork, &info);
//    if (info != 0)
//       throw std::runtime_error("MakeRotatedSlmXTransform: dgelss failed.");
//    if (nRank != nSl)
//       throw std::runtime_error("MakeRotatedSlmXTransform: equation system solution went wrong.");
//
//    // copy solution to output matrix.
//    for (unsigned iSl = 0; iSl < nSl; ++ iSl)
//       for (unsigned jSl = 0; jSl < nSl; ++ jSl)
//          pOut[iSl + nStride * jSl] = pBt[iSl + nCa * jSl];
   size_t
      nRank = ct::LeastSquaresSolve(pOut, nStride, pAt, nCa, pBt, nCa, nCa, nSl, nSl, 1e-10);
   if (nRank != nSl)
      throw std::runtime_error("MakeRotatedSlmXTransform: equation system solution went wrong.");


//    for (unsigned l = 0; l <= MaxL; ++ l) {
//       double
//          *IR_RP pT = &pOut[what?].
//       (use addressing as above.. just make svds separately)
//    }
}

   // TODO @ Update (2021-08-10) of EvalSlcXRotationTrafo:
   // ----------------------------------------------------
   //
   // - I did a cursory scan of the literature. There are formulas for Ylm
   //   rotations around, ...but
   //   (a) published ones are usually for complex harmonics. In principle
   //   straight forward to convert, but still nasty.
   //   (b) I found many variants. They are all horrible (Euler angles, Wigner-D
   //   matrices and similar insanities),
   //   (c) I found several published approaches, some with many citations,
   //   which I would consider even WORSE than mine here.
   //
   // - So the equation system variant may not be that bad, after all. In
   //   principle, for the smaller system some symbolic solution should be quite
   //   feasible. But I played around with it for a bit and did not manage to
   //   get Mathematica to produce it in a reasonable form on the first try.
   //   Certainly not in a form scalable to L=8 at least (we use those for
   //   fitting sets)
   //
   // - For the numerical variant, I think we need to replace the SVD with a
   //   basic SVD or LU code (IR really cannot have a LAPACK dependency); I
   //   think I had made one, somewhere, not so long ago.
   //
   // - Regarding the representation of the result:
   //
   //   + Result should probably be represented with a more suitable data
   //     structure than the nSlcX(MaxL)^2 block-diagonal matrix I have here at
   //     the moment.
   //   + ...but I presently do not have a clean overview over where and how
   //     exactly the trafos are used (definitely in IboView, and I think
   //     MicroScf had something like this for initial guess propagation,
   //     somewhere). Need to investigate and update.
   //   + I think a higher meta data structure would be in order, even if
   //     we get sets of (2l+1)^2 block matrices of fixed structure. A little
   //     meta structure to simplify access to them.
   //
   // - Regarding the concrete setup of the equation systems:
   //
   //   + The rotation trafos SHOULD be solved for in the different L
   //     subspaces individually.
   //
   //   + Technically, we just need a set of a sufficient number (2l+1) of
   //     distinct real-space points, and their R-transforms, to derive the
   //     equation system uniquely, following the ideas above.
   //
   //     One such set of 2l+1 distinct points would be:
   //
   //           p_0     = (  l,   0,  0)
   //           p_1     = (l-1,   1,  0)
   //           p_2     = (l-2,   2,  0)
   //           ...       ...
   //           p_l     = (  1, l-1,  0)
   //           p_l     = (  0,   l,  0)
   //           p_{l+1} = (  0, l-1,  1)
   //           p_{l+2} = (  0, l-2,  2)
   //           ...       ...
   //           p_{2l}  = (  0,   0,  l)
   //
   //     - Note that these are points, not vectors. And linear dependence is not
   //       actually a measure, because the S^l_m are not linear functions.
   //     - I guess things like conditioning would still exist, so not all kinds
   //       of point sets would be numerically quite equivalent; but it is hard
   //       to estimate what effect that would have, and what better choices
   //       should look like (apart form normalized points on the unit sphere
   //       probaby being a good idea)
   //
   //   + The transformation T_l  we mean to find is a (2l+1, 2l+1)-shape
   //     matrix, which is supposed to have the following outcome:
   //
   //           S[l,c](R r) = \sum_{c'=0}^{2l} T_l[c,c'] S[l',c'](r)    (A1)
   //
   //     So once we have those trial points, we just use them to set up
   //     multiple rows of (A1), just enough to uniquely determine T:
   //
   //           A_{iRow,iCol} := S[l,iRow](p_{iCol})
   //           B_{iRow,iCol} := S[l,iRow](R⋅p_{iCol})
   //
   //           B = T * A
   //           ⇔ B^⊹ = A^⊹ * T^⊹
   //
   //     ...and that's what we'd solve for. T = solve(A^⊹, B^⊹)^⊹
   //
   //   + A bit of an issue is that we'd actually get *all* low Ls, too. In the Slc
   //     components. Might need to do something about that. E.g., a different ordering
   //     of the points which can be shared across all Ls. Can that be done?
   //
   //           p_{  0} = (0, 0, 1)
   //           p_{  1} = (1, 0, 0)
   //           p_{  2} = (0, 1, 0)        // <-- L = 1 boundary
   //
   //           p_{  3} = (1, 0, 1)/√2
   //           p_{  4} = (0, 1, 1)/√2     // <-- L = 2 boundary
   //
   //           p_{  5} = (2, 0, 1)/√5
   //           p_{  6} = (0, 2, 1)/√5     // <-- L = 3 boundary
   //
   //           p_{  5} = (3, 0, 1)/√10
   //           p_{  6} = (0, 3, 1)/√10    // <-- L = 4 boundary
   //
   //     ...apparently it can. And it does become clearer what is going on.
   //     Instead of fixing all the z elements, we could also employ the following
   //     variant:
   //
   //           p_{  0} = (0, 0, 1)
   //
   //           p_{  1} = (1, 0, 0)
   //           p_{  2} = (0, 1, 0)        // <-- L = 1 boundary
   //
   //           p_{  3} = (1, 0, 1)/√2
   //           p_{  4} = (0, 1, 1)/√2     // <-- L = 2 boundary
   //
   //           p_{  5} = (1, 0, 2)/√5
   //           p_{  6} = (0, 1, 2)/√5     // <-- L = 3 boundary
   //
   //           p_{  5} = (1, 0, 3)/√10
   //           p_{  6} = (0, 1, 3)/√10    // <-- L = 4 boundary
   //
   //     (note: the sqrts are Sqrt[(L-1)^2 + 1]. I guess just dividing by (L-1)
   //     would also do the trick if I use solid harmonics. But maybe better go
   //     with unit sphere and do full normalization...)
   //
   //   + Well... all in all, this seems quite doable. Need to have a closer
   //     look into how to do the transpositions, though---I get the c/c' components
   //     as rows, but would need them as cols (logical cols, that is, whether
   //     that agrees with the storage order is another matter).
   //     Oops... the code above also does the matrix copy. In my case I'd
   //     need to preserve the original matrices, rather than overriding rhs,
   //     so a copy may anyway be needed.
   //



} // namespace ir
