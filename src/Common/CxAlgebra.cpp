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

#include <algorithm> // for std::min
#include <stdexcept>
#include <stdlib.h>
#include "format.h"
#include "CxPodArray.h"

// #define USE_EIGEN_LA
#ifdef USE_EIGEN_LA
   #include <Eigen/Cholesky>
   #include <Eigen/SVD>
   #include <Eigen/Eigenvalues>

   #define EIGEN_SVD_ALGO JacobiSVD // straight up jacobi SVD (more accurate)
   // #define EIGEN_SVD_ALGO BDCSVD // bidirectional divide&conquer SVD (faster)

   using FMatrix1 = Eigen::MatrixXd;
   enum {
      EIGEN_DEFAULT_ALIGN = Eigen::Unaligned
//       EIGEN_DEFAULT_ALIGN = 8
//       EIGEN_DEFAULT_ALIGN = 64
   };
   using FMatrix1Map = Eigen::Map<FMatrix1, EIGEN_DEFAULT_ALIGN, Eigen::OuterStride<> >; // <-- should probably put some assert around that....
   using FVector1Map = Eigen::Map<Eigen::VectorXd>;
#endif


#include "CxAlgebra.h"
#include "CxDefs.h" // for assert.

namespace ct {


void Mxm(double *pOut, ptrdiff_t iRowStO, ptrdiff_t iColStO,
         double const *pA, ptrdiff_t iRowStA, ptrdiff_t iColStA,
         double const *pB, ptrdiff_t iRowStB, ptrdiff_t iColStB,
         size_t nRows, size_t nLink, size_t nCols, bool AddToDest, double fFactor)
{
   return Mxm_InlineCalls(pOut, iRowStO, iColStO, pA, iRowStA, iColStA, pB, iRowStB, iColStB, nRows, nLink, nCols, AddToDest, fFactor);
}


// this one is used if neither the column nor the row stride in pMat is unity.
void MxvLameG(double *RESTRICT pOut, ptrdiff_t iStO, double const *RESTRICT pMat, ptrdiff_t iRowStM, ptrdiff_t iColStM,
   double const *RESTRICT pIn, ptrdiff_t iStI, size_t nRows, size_t nLink, bool AddToDest, double fFactor)
{
   for (size_t iRow = 0; iRow < nRows; ++ iRow) {
      double
         d = 0;
      double const
         *RESTRICT pM = &pMat[iRowStM * iRow];
      for (size_t iLink = 0; iLink < nLink; ++ iLink) {
         d += pIn[iStI * iLink] * pM[iColStM * iLink];
      }
      d *= fFactor;

      double *RESTRICT r = &pOut[iStO * iRow];
      if (AddToDest)
         *r += d;
      else
         *r = d;
   }
}

// this one is used if neither the column nor the row stride in pMat is unity.
void MxvLame(double *RESTRICT pOut, ptrdiff_t iStO, double const *RESTRICT pMat, ptrdiff_t iRowStM, ptrdiff_t iColStM,
   double const *RESTRICT pIn, ptrdiff_t iStI, size_t nRows, size_t nLink, bool AddToDest, double fFactor)
{
   if (iStO != 1 || iStI != 1)
      return MxvLameG(pOut, iStO, pMat, iRowStM, iColStM, pIn, iStI, nRows, nLink, AddToDest, fFactor);
   for (size_t iRow = 0; iRow < nRows; ++iRow) {
      double
         d = 0;
      double const
         *RESTRICT pM = &pMat[iRowStM * iRow];
      if (iColStM == 1) {
         for (size_t iLink = 0; iLink < nLink; ++iLink) {
            d += pIn[iLink] * pM[iLink];
         }
      } else {
         for (size_t iLink = 0; iLink < nLink; ++iLink) {
            d += pIn[iLink] * pM[iColStM * iLink];
         }
      }
      d *= fFactor;

      double *RESTRICT r = &pOut[iRow];
      if (AddToDest)
         *r += d;
      else
         *r = d;
   }
}




#if defined(USE_EIGEN_LA) && true
   void Diagonalize(double *pEw, double *pH, unsigned ldH, unsigned N)
   {
      FMatrix1Map
         H(pH, N, N, Eigen::OuterStride<>(ldH));
      Eigen::SelfAdjointEigenSolver<FMatrix1>
         sd(H, Eigen::ComputeEigenvectors); // compute spectral decomposition
      if (sd.info() != Eigen::Success)
         throw std::runtime_error(fmt::format("Eigen::SelfAdjointEigenSolver: failed to compute factorization A = V diag(ew) V^T for input matrix of shape ({},{})", N, N));
      FVector1Map
         Ew(pEw, N);
      Ew = sd.eigenvalues();
      H = sd.eigenvectors();
   }
#else
   // using LAPACK
   void Diagonalize(double *pEw, double *pH, unsigned ldH, unsigned N)
   {
      if (N == 0)
         return;
      FORTINT info = 0;
      // workspace query.
      double fWork = 0;
      FORTINT nWork, niWork[2] = {0};
      DSYEVD('V', 'L', N, pH, ldH, pEw, &fWork, -1, &niWork[0], -1, info);
      if (info != 0) throw std::runtime_error("dsyevd workspace query failed.");
      nWork = FORTINT(fWork);

      double *pWork = (double*)::malloc(sizeof(double)*nWork+1000);
      FORTINT *piWork = (FORTINT*)::malloc(sizeof(FORTINT)*niWork[0]+1000);
      DSYEVD('V', 'L', N, pH, ldH, pEw, pWork, nWork, piWork, niWork[0], info);
      ::free(piWork);
      ::free(pWork);
      if (info != 0) throw std::runtime_error("dsyevd failed.");
   }
#endif // USE_EIGEN_LA


#if defined(USE_EIGEN_LA) && true
   void CalcCholeskyFactors(double *pM, size_t ldM, size_t N, bool UtU_insteadOf_LLt) {
      FMatrix1Map
         M(pM, N, N, Eigen::OuterStride<>(ldM));
      if (!UtU_insteadOf_LLt) {
         Eigen::LLT<FMatrix1, Eigen::Lower> cd(M);
         if (cd.info() != Eigen::Success)
            throw std::runtime_error(fmt::format("Eigen::LLT: failed to compute factorization A = L*L^T for input matrix of shape ({},{})", N, N));
         M = cd.matrixL();
      } else {
         Eigen::LLT<FMatrix1, Eigen::Upper> cd(M);
         if (cd.info() != Eigen::Success)
            throw std::runtime_error(fmt::format("Eigen::LLT: failed to compute factorization A = U^T*U for input matrix of shape ({},{})", N, N));
         M = cd.matrixU();
      }
   }

   void InvertTriangularMatrix(double *pL, size_t ldL, size_t N) {
      FMatrix1
         I = FMatrix1::Identity(N,N);
      FMatrix1Map
         L(pL, N, N, Eigen::OuterStride<>(ldL));
      L.triangularView<Eigen::Lower>().solveInPlace(I);
      L = I; // this should also be a lower triangular matrix.
   }
#else
   void CalcCholeskyFactors(double *pM, size_t ldM, size_t N, bool UtU_insteadOf_LLt) {
      FORTINT
         info = 0;
      // factorize M into M = L * L^T or M = U^T * U
      DPOTRF(UtU_insteadOf_LLt? 'U' : 'L', N, pM, ldM, &info);
      if (info != 0) {
         throw std::runtime_error(fmt::format("CalcCholeskyFactors: dpotrf failed on input matrix of shape ({},{}) with strides [{}, {}]. LAPACK returned info code {}.", N, N, 1, ldM, info));
      }
   }

   void InvertTriangularMatrix(double *pL, size_t ldL, size_t N) {
      FORTINT
         info = 0;
      DTRTRI('L', 'N', N, pL, ldL, &info);
      if (info != 0)
         throw std::runtime_error("InvertTriangularMatrix: dtrtri failed.");
   }
#endif // USE_EIGEN_LA



#if defined(USE_EIGEN_LA) && true
   // In: (nRows, nCols), U: (nRows, nSig), Vt: (nSig, nCols)
   // where nSig = min(nRows, nCols)
   void ComputeSvd(double *pU, size_t ldU, double *pSigma, double *pVt, size_t ldVt, double *pInAndTmp, size_t ldIn, size_t nRows, size_t nCols)
   {
      size_t
         nSig = std::min(nRows, nCols);
      assert(ldU >= nRows && ldVt >= nSig && ldIn >= nRows);
      FMatrix1Map
         U(pU, nRows, nSig, Eigen::OuterStride<>(ldU)),
         Vt(pVt, nSig, nCols, Eigen::OuterStride<>(ldVt)),
         In(pInAndTmp, nRows, nCols, Eigen::OuterStride<>(ldIn));
      FVector1Map
         sigma(pSigma, nSig);
      Eigen::EIGEN_SVD_ALGO<FMatrix1>
         svd(In, Eigen::ComputeThinU | Eigen::ComputeThinV);
      if (svd.info() != Eigen::Success) {
         throw std::runtime_error(fmt::format("EIGEN_SVD_ALGO failed for input matrix of shape ({},{})", nRows, nCols));
      }
      assert(size_t(svd.singularValues().size()) == nSig);
      assert(Vt.rows() == svd.matrixV().cols() && Vt.cols() == svd.matrixV().rows());
      assert(size_t(svd.matrixU().rows()) == nRows && size_t(svd.matrixU().cols()) == nSig);
      assert(size_t(svd.matrixV().rows()) == nCols && size_t(svd.matrixV().cols()) == nSig);
      sigma = svd.singularValues();
      U = svd.matrixU();
      Vt = svd.matrixV().transpose();
   }

#else
   // In: (nRows, nCols), U: (nRows, nSig), Vt: (nSig, nCols)
   // where nSig = min(nRows, nCols)
   void ComputeSvd(double *pU, size_t ldU, double *pSigma, double *pVt, size_t ldVt, double *pInAndTmp, size_t ldIn, size_t nRows, size_t nCols)
   {
      size_t
         nSig = std::min(nRows, nCols);
      assert(ldU >= nRows && ldVt >= nSig && ldIn >= nRows);
      FORTINT
         lWork = 0,
         info = 0;
      FORTINT
         *piWork = (FORTINT*)::malloc(sizeof(FORTINT)*8*nSig);
      // workspace query.
      double
         flWork = 0;
      DGESDD('S', nRows, nCols, pInAndTmp, ldIn, pSigma, pU, ldU, pVt, ldVt, &flWork, -1, piWork, &info);
      if (info != 0) {
         ::free(piWork); // my understanding of the docs is that piWork needs to be valid for the workspace query...
         throw std::runtime_error("dgesdd workspace query failed.");
      }
      lWork = FORTINT(flWork);
      double
         *pWork = (double*)::malloc(sizeof(double)*lWork);
      DGESDD('S', nRows, nCols, pInAndTmp, ldIn, pSigma, pU, ldU, pVt, ldVt, pWork, lWork, piWork, &info);
      ::free(pWork);
      ::free(piWork);
      if (info != 0)
         throw std::runtime_error("dgesdd failed.");
   }
#endif // USE_EIGEN_LA


#if defined(USE_EIGEN_LA) && true
   // solve ||A * sol - rhs||_2 -> min
   // A: (nRows, nCols)
   // sol: (nCols, nVec)
   // rhs: (nRows, nVec)
   // returns rank of A.
   size_t LeastSquaresSolve(double *pSol, size_t ldSol, double const *pA, size_t ldA, double const *pRhs, size_t ldRhs, size_t nRows, size_t nCols, size_t nVec, double ThrSig) {
      // should probably also compute the rank...
      assert_rt(ldSol >= nCols && ldRhs >= nRows && ldA >= nRows);
      FMatrix1Map
         A(const_cast<double*>(pA), nRows, nCols, Eigen::OuterStride<>(ldA)),
         Sol(pSol, nCols, nVec, Eigen::OuterStride<>(ldSol)),
         Rhs(const_cast<double*>(pRhs), nRows, nVec, Eigen::OuterStride<>(ldRhs));
      Eigen::EIGEN_SVD_ALGO<FMatrix1>
         svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
      if (svd.info() != Eigen::Success) {
         throw std::runtime_error(fmt::format("EIGEN_SVD_ALGO failed for input matrix of shape ({},{})", nRows, nCols));
      }
      svd.setThreshold(ThrSig);
      Sol = svd.solve(Rhs);
      return size_t(svd.rank());
   }

#else
   // solve linear least squares ||M * x - b||_2 -> min subject to ||x||_2 ->min
   // using singular value decomposition of A. Singular values below fThrRel are
   // treated as zero. A negative value indicates the usage of machine precision.
   //
   // solve ||A * sol - rhs||_2 -> min
   // A: (nRows, nCols)
   // sol: (nCols, nVec)
   // rhs: (nRows, nVec)
   // returns rank of A.
   size_t LeastSquaresSolve(double *pSol, size_t ldSol, double const *pA, size_t ldA, double const *pRhs, size_t ldRhs, size_t nRows, size_t nCols, size_t nVec, double ThrSig) {
      // should probably also compute the rank...
      using std::max;
      using std::min;
      assert(ldSol >= nCols && ldRhs >= nRows && ldA >= nRows);
      size_t
         ldRhsSol = std::max(nRows, nCols);
      TArray<double>
         ACopy(nRows * nCols),
         RhsSol(ldRhsSol * nVec);
      for (size_t c = 0; c < nCols; ++ c)
         for (size_t r = 0; r < nRows; ++ r)
            ACopy[r + nRows*c] = pA[r + ldA*c];
      for (size_t v = 0; v < nVec; ++ v)
         for (size_t r = 0; r < nRows; ++ r)
            RhsSol[r + ldRhsSol*v] = pRhs[r + ldRhs*v];

      // TODO:
      // - get rid of DGELSS and rephrase it in terms of ComputeSvd?
      // - The netlib DGELSS source actually does *not* appear to use the SVD routines
      //   directly in many cases, but rather does some funky QR factorizations
      //   (didn't look into it in detail).
      //   I just wonder if it is worth the extra LAPACK dependencies...
      //   the SVD version should probably be good enough in the vast majority of cases.
      FORTINT
         nRank = -1,
         info = 0,
         // minimal work space required according to docs.
         nWork = 3*min(nRows,nCols) + max(max(2*min(nRows, nCols), max(nRows, nCols)), nVec);
      double
         fWork = 0;
      TArray<double>
         pSig(min(nRows, nCols));  // singular values (decreasing order).
      DGELSS(nRows, nCols, nVec, &ACopy[0], nRows,
         &RhsSol[0], ldRhsSol, &pSig[0], ThrSig, &nRank,
         &fWork, -1, &info);
      if (info != 0) throw std::runtime_error("dgelss failed in workspace query.");
      nWork = static_cast<FORTINT>(fWork + 0.5);
      TArray<double>
         pWork(nWork);
      DGELSS(nRows, nCols, nVec, &ACopy[0], nRows,
         &RhsSol[0], ldRhsSol, &pSig[0], ThrSig, &nRank,
         &pWork[0], nWork, &info);
      if (info != 0) throw std::runtime_error("dgelss failed.");
      // output @ RhsSol: solution. M.nCols x nRhs

      for (size_t v = 0; v < nVec; ++ v)
         for (size_t c = 0; c < nCols; ++ c)
            pSol[c + ldSol*v] = RhsSol[c + ldRhsSol*v];
      return size_t(nRank);
      // wow... that is sooooo ugly.
   };

#endif // USE_EIGEN_LA


void LinearSolveSymSvd(double *pOut, double *pMat, uint nStr, double *pIn, uint nDim, double Thr, FMemoryStack &Mem)
{
    TMemoryLock<double>
        pEvs(nDim*nDim, &Mem),
        pEws(nDim, &Mem),
        pXv(nDim, &Mem); // input vectors in EV basis.
    for (size_t j = 0; j < nDim; ++ j)
        for (size_t i = 0; i < nDim; ++ i)
            pEvs[i + nDim*j] = -pMat[i + nStr*j]; // <- beware of minus.
    Diagonalize(pEws, pEvs, nDim, nDim);

    Mxv(pXv,1, pEvs,nDim,1,  pIn,1, nDim,nDim);
    for (size_t iEw = 0; iEw != nDim; ++ iEw)
        if (std::abs(pEws[iEw]) >= Thr)
            pXv[iEw] /= -pEws[iEw];
            // ^- note that this only screens by absolute value.
            // no positive semi-definiteness is assumed!
        else
            pXv[iEw] = 0.;
    Mxv(pOut,1, pEvs,1,nDim, pXv,1, nDim,nDim);
}



}

// kate: indent-width 3
