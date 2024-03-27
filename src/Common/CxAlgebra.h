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

#ifndef CX_ALGEBRA_H
#define CX_ALGEBRA_H

// #include <iostream> // for debug print. TODO: remove this.
// #include "format.h" // for debug print. TODO: remove this.

#include <cstddef>
#include "CxDefs.h"
#include "CxFortranInt.h"
#include "CxMemoryStack.h" // only used in LinearSolveSymSvd atm, but if we move some other general stuff here, it could be more general.

using std::ptrdiff_t;
using std::size_t;
typedef unsigned int uint;

// behold the elegance of FORTRAN77/C interfaces!
extern "C"{
      // declaration of FORTRAN77 blas/lapack routines... (to exist in ACML/MKL/etc.)

      //  C := alpha*op( A )*op( B ) + beta*C,
      // Trans*: 'N' (notranspose), 'T' (transpose) or 'C' (conjugate-transpose)(=T)
      #define DGEMM FORT_Extern(dgemm,DGEMM)
      void DGEMM( char const &TransA, char const &TransB, FINTARG M, FINTARG N, FINTARG K,
          FDBLARG alpha, double const *A, FINTARG lda, double const *B, FINTARG ldb,
          FDBLARG beta, double *C, FINTARG ldc );

      //  C := alpha*op(A)*op(A^T) + beta*C,
      #define DSYRK FORT_Extern(dsyrk,DSYRK)
      void DSYRK( char const &UpLo, char const &TransA, FINTARG N, FINTARG K,
          FDBLARG alpha, double const *A, FINTARG lda,
          FDBLARG beta, double *C, FINTARG ldc );

      // y += alpha * M * x
      #define DGEMV FORT_Extern(dgemv,DGEMV)
      void DGEMV(char const &Trans, FINTARG M, FINTARG N, FDBLARG Alpha, double const *RESTRICT A, FINTARG lda, double const *RESTRICT X, FINTARG incx, FDBLARG Beta, double *RESTRICT y, FINTARG incy);

      // A += alpha * x y^T
      #define DGER FORT_Extern(dger,DGER)
      void DGER(FINTARG M, FINTARG N, double const &Alpha, double const *RESTRICT X, FINTARG incx, double const *RESTRICT Y, FINTARG incy, double *RESTRICT A, FINTARG ldA);

      // computes eigenvalues and eigenvectors of a symmetric matrix:
      //  jobz: 'N'/'V': compute eigenvalues only/compute also eigenvectors
      //  uplo: 'U'/'L': upper/lower triangle of A is stored.
      //  N: matrix size. A is N x N matrix.
      //  A: input matrix, output vectors. LDA: row stride A.
      //  W: output, eigenvectors in ascending order. (vector of length N).
      //  Work: work space
      //  lWork: Work space size. "For optimal efficiency, LWORK >= (NB+2)*N,
      // where NB is the blocksize for DSYTRD returned by ILAENV."
      #define DSYEV FORT_Extern(dsyev,DSYEV)
      void DSYEV(char const &jobz, char const &uplo, FINTARG N, double *A, FINTARG LDA, double *W, double *WORK, FINTARG LWORK, FORTINT &INFO);
      // divide & conquer spectral decomposition
      // WARNING:
      // - in my f2c-LAPACK hack with Eigen this one fails sometimes on win32!
      // - Don't know if it's a porting problem or something else.
      // - Or possibly OpenGL resetting the floating point control word to
      //   'single precision' or something. It used to do that in earlier times...
      // - Should be checked and repaired!
      #define DSYEVD FORT_Extern(dsyevd,DSYEVD)
      void DSYEVD(char const &jobz, char const &uplo, FINTARG N, double *A, FINTARG LDA, double *W, double *WORK, FINTARG LWORK, FORTINT *IWORK, FINTARG LIWORK, FORTINT &INFO);

//       // solves real generalized eigensystem
//       //    S C = H C
//       // with S, H symmetric and S positive definite
//       #define DSYGV FORT_Extern(dsygv,DSYGV)
//       void DSYGV(FINTARG ITYPE, char const &JOBZ, char const &UPLO, FINTARG N,
//          double *A, FINTARG LDA, double const *B, FINTARG LDB, double *EW,
//          double *WORK, FORTINT &LWORK, FORTINT &INFO );

      // calculate Cholesky factorization of A
      #define DPOTRF FORT_Extern(dpotrf,DPOTRF)
      void DPOTRF(char const &UpLo, FINTARG n, double *A, FINTARG lda, FORTINT *info);
//       #define DPOTRS FORT_Extern(dpotrs,DPOTRS)
//       void DPOTRS(char const &UpLo, FINTARG n, FINTARG nRhs, double *A, FINTARG lda, double *B, FINTARG ldb, FORTINT *info);
//       #define DTRTRS FORT_Extern(dtrtrs,DTRTRS)
//       void DTRTRS(char const &UpLo, char const &Trans, char const &Diag, FINTARG N, FINTARG NRHS, double *A, FINTARG lda, double *B, FINTARG ldb, FORTINT *info);
//       // ^- gna.. dtrtrs is rather useless. It just does some argument checks and
//       // then calls dtrsm with side == 'L' (which probably does the same checks again).
      // solve triangular equation system op(L) X = B
      #define DTRSM FORT_Extern(dtrsm,DTRSM)
      void DTRSM(char const &Side, char const &UpLo, char const &Trans, char const &Diag, FINTARG nRowsB, FINTARG nColsB, double const &Alpha, double *A, FINTARG lda, double *B, FINTARG ldb, FORTINT *info);
      // multiply with triangular matrix: B := op(L) * B or B := B * op(L)
      #define DTRMM FORT_Extern(dtrmm,DTRMM)
      void DTRMM(char const &Side, char const &UpLo, char const &Trans, char const &Diag, FINTARG nRowsB, FINTARG nColsB, double const &Alpha, double *A, FINTARG lda, double *B, FINTARG ldb);
      // compute inverse of a real upper or lower triangular matrix
      #define DTRTRI FORT_Extern(dtrtri,DTRTRI)
      void DTRTRI(char const &Uplo, char const &Diag, FINTARG N, double *pA, FINTARG LDA, FORTINT *info);

      // compute m x n matrix LU factorization.
      // info: =0 success. > 0: matrix is singular, factorization cannot be used
      // to solve linear systems.
      #define DGETRF FORT_Extern(dgetrf,DGETRF)
      void DGETRF(FINTARG M, FINTARG N, double const *pA, FINTARG LDA, FORTINT *ipiv, FORTINT *INFO );
      // solves A * X = B for X. n: number of equations (order of A).
      // needs LU decomposition as input (LU+partial pivoting; e.g., DGETRF).
      #define DGESV FORT_Extern(dgesv,DGESV)
      void DGESV( FINTARG n, FINTARG nrhs, double *A, FINTARG lda, FINTARG ipivot, double *B,
          FINTARG ldb, FORTINT &info );

//       // compute SVD (simple driver for small cases)
//       #define DGESVD FORT_Extern(dgesvd,DGESVD)
//       void DGESVD(char const &JobU, char const &JobVT, FINTARG M, FINTARG N, double *A, FINTARG ldA, double *S, double *U, FINTARG ldU, double *Vt, FINTARG ldVt, double *Work, FINTARG lWork, FORTINT *info);
//           double *WORK, FINTARG LWORK, FORTINT *INFO );
      // divide & conquer SVD.
      #define DGESDD FORT_Extern(dgesdd,DGESDD)
      void DGESDD(char const &JobZ, FINTARG M, FINTARG N, double *A, FINTARG ldA, double *S, double *U, FINTARG ldU, double *Vt, FINTARG ldVt, double *Work, FINTARG lWork, FORTINT *piWork, FORTINT *info);
      // linear least squares using divide & conquer SVD.
      #define DGELSS FORT_Extern(dgelss,DGELSS)
      void DGELSS(FINTARG M, FINTARG N, FINTARG NRHS, double *A, FINTARG LDA, double *B, FINTARG LDB, double *S, double const &RCOND, FORTINT *RANK,
          double *WORK, FINTARG LWORK, FORTINT *INFO);

// FIXME: This one is entirely unnecessary, but was used in ItfMcaUtils.cpp. Should be removed.
      #define DAXPY FORT_Extern(daxpy,DAXPY)
      void DAXPY(FINTARG M, FDBLARG Alpha, double const *A, FINTARG lda, double *B, FINTARG ldb);

}





namespace ct {

// if TransX says 'no transpose' -> return 'transpose',
// if TransX says 'transpose' -> return 'no transpose'
inline char InvertBlasTransposeArg(char TransX) {
   assert(TransX == 'N' || TransX == 'T');
   return TransX == 'N'? 'T' : 'N';
}

void Mxm(double *pOut, ptrdiff_t iRowStO, ptrdiff_t iColStO,
         double const *pA, ptrdiff_t iRowStA, ptrdiff_t iColStA,
         double const *pB, ptrdiff_t iRowStB, ptrdiff_t iColStB,
         size_t nRows, size_t nLink, size_t nCols, bool AddToDest = false, double fFactor = 1.0);

// note: both H and S are overwritten. Eigenvectors go into H.
void DiagonalizeGen(double *pEw, double *pH, unsigned ldH, double *pS, unsigned ldS, unsigned N);
void Diagonalize(double *pEw, double *pH, unsigned ldH, unsigned N);
// note: pInAndTmp will be overwritten. Output is U: (nRows, nSig), Vt : (nSig, nCols), nSig = min(nRows,nCols).
void ComputeSvd(double *pU, size_t ldU, double *pSigma, double *pVt, size_t ldVt, double *pInAndTmp, size_t ldIn, size_t nRows, size_t nCols);
void LinearSolveSymSvd(double *pOut, double *pMat, unsigned nStr, double *pIn, unsigned nDim, double Thr, FMemoryStack &Mem);
void CalcCholeskyFactors(double *pM, size_t ldM, size_t N, bool UtU_insteadOf_LLt);
void InvertTriangularMatrix(double *pL, size_t ldL, size_t N);
size_t LeastSquaresSolve(double *pSol, size_t ldSol, double const *pA, size_t ldA, double const *pRhs, size_t ldRhs, size_t nRows, size_t nCols, size_t nVec, double ThrSig);



void MxvLame(double *RESTRICT pOut, ptrdiff_t iStO, double const *RESTRICT pMat, ptrdiff_t iRowStM, ptrdiff_t iColStM,
    double const *RESTRICT pIn, ptrdiff_t iStI, size_t nRows, size_t nLink, bool AddToDest = false, double fFactor = 1.0);


inline void Mxv(double *RESTRICT pOut, ptrdiff_t iStO, double const *RESTRICT pMat, ptrdiff_t iRowStM, ptrdiff_t iColStM,
    double const *RESTRICT pIn, ptrdiff_t iStI, size_t nRows, size_t nLink, bool AddToDest = false, double fFactor = 1.0)
{
//    if (iRowStM != 1 && iColStM != 1) {
//       return MxvLame(pOut, iStO, pMat, iRowStM, iColStM, pIn, iStI, nRows, nLink, AddToDest, fFactor);
//    };

   char OpM;
   if (iRowStM == 1 && iColStM == 1) {
      // pMat either has only one column or only one row. In this case
      // both strides can be 1 at the same time, so we cannot use them
      // to distinguish the 'N' vs 'T' case. See comments on Mxm_InlineCalls
      // how the 'N' vs 'T' business works in this case.
      OpM = (nRows <= 1) ? 'N' : 'T';
   } else if (iRowStM == 1) {
      OpM = 'N';
   } else if (iColStM == 1) {
      OpM = 'T';
   } else {
      assert_rt(!"Mxv(...): encountered iRowStM and iColStM both being non-unity. Currently not allowed.");
      return;
   }

   if (OpM == 'N') {
      assert(iRowStM == 1);
      DGEMV('N', nRows, nLink, fFactor, pMat, iColStM,  pIn,iStI,  AddToDest? 1.0 : 0.0, pOut,iStO);
   } else {
      assert(iColStM == 1);
      DGEMV('T', nLink, nRows, fFactor, pMat, iRowStM,  pIn,iStI,  AddToDest? 1.0 : 0.0, pOut,iStO);
   }
}


// Out = f * A * B
// Note: inline version is primarily for putting it into CtMatrix.h; doing this
// might not actually be a good idea---might want to put more complex drivers
// (e.g., libxsmm calls if present) and case distinctions into this.
IR_FORCE_INLINE void Mxm_InlineCalls(double *pOut, ptrdiff_t iRowStO, ptrdiff_t iColStO,
         double const *pA, ptrdiff_t iRowStA, ptrdiff_t iColStA,
         double const *pB, ptrdiff_t iRowStB, ptrdiff_t iColStB,
         size_t nRows, size_t nLink, size_t nCols, bool AddToDest, double fFactor)
{
   if (nRows == 0 || nCols == 0)
      // no output matrix. And could easily confuse MKL's DGEMM. Better get out
      // of here asap...
      return;

   assert(iRowStO == 1 || iColStO == 1);
   assert(iRowStA == 1 || iColStA == 1);
   assert(iRowStB == 1 || iColStB == 1);
   // ^- otherwise dgemm directly not applicable. Would need local copy
   // of matrix/matrices with compressed strides.

   double
      Beta = AddToDest? 1.0 : 0.0;
   char
      TransA, TransB;
   FORTINT
      lda, ldb;

   // now for the fun part: deciding how to communicate the memory layout of our
   // matrices to dgemm, by claiming they are/are not to be transposed.

   // First 'A':
   if (iRowStA == 1 && iColStA == 1) {
      // 'A' either has only one row or only one col... in which case both strides
      // can be 1 at the same time. In either case, we get lda=1:
      lda = 1;
      // However, this way we cannot use the strides for deciding whether 'A'
      // should be "transposed" or not (and MKL complains bitterly if you get it
      // wrong, not that it matters in this case...). So here we have to use
      // nRows ("M"), nLink ("K"), nCols ("N"), to decide the 'N' vs 'T'
      // situation: op(A) is in either case a (nRows,nLink)-shape matrix. This
      // means:
      if (nRows <= 1) {
         // ...if nRows==1, MKL will be happy with A being a (1,nLink) matrix
         // with lda==1, so we need to set TransA='N'.
         TransA = 'N';
      } else {
         // if nRows is not 1, however, we need A to be a (nLink,nRows) matrix
         // with nLink=1 (so that op(A) is a (nRows,nLink) matrix). So if we set
         // TransA='T' and nLink is indeed <= 1, everything is fine.
         assert(nLink <= 1);
         TransA = 'T';
      }
   } else if (iRowStA == 1) {
      TransA = 'N'; lda = iColStA;
   } else {
      TransA = 'T'; lda = iRowStA;
   }

   // Next in line is 'B':
   if (iRowStB == 1 && iColStB == 1) {
      ldb = 1;
      // ...with the same arguments to be considered if both strides are 1:
      // - Op(B) is a (nLink,nCols)-shape matrix.
      // - So if Op(B)=B (i.e., TransB='N'), B is itself a (nLink,nCols)-shape
      //   matrix. So if nLink<=1 while ldb==1, everything should be fine.
      // - Similarly, with Op(B)=B^T (i.e., TransB='T'), B is a
      //   (nCols,nLink)-shape matrix, which for which havind ldb=1 should be
      //   acceptable to BLAS if nCols <= 1.
      if (nCols <= 1) {
         TransB = 'T';
      } else {
         assert(nLink <= 1);
         TransB = 'N';
      }
   } else if (iRowStB == 1) {
      TransB = 'N'; ldb = iColStB;
   } else {
      TransB = 'T'; ldb = iRowStB;
   }

   assert(lda != 0 && ldb != 0);

   char TransC;
   if (iRowStO == 1 && iColStO == 1) {
      // Out either has only one row or only one column. This time the
      // N vs T game says:
      if (nRows == 1) {
         // ...if TransC=='N', then C is a (nRows,nCols)-shape matrix directly.
         // for ldc=1 to be acceptable to BLAS, this needs nRows==1.
         TransC = 'N';
      } else {
         // ...if TransC=='T', then C is a (nCols,nRows)-shape matrix (because
         // op(C) has (nRows,nCols)-shape). for ldc=1 to be acceptable to BLAS,
         // this needs nCols==1.
         assert(nCols <= 1);
         TransC = 'T';
      }
   } else if (iRowStO == 1) {
      TransC = 'N';
   } else {
      assert(iColStO == 1);
      TransC = 'T';
   }
   if (TransC == 'N') {
      assert(iColStO != 0);
      assert(nRows == 0 || nCols == 0 || (pOut != pA && pOut != pB)); // input and output matrices cannot be aliased.
      // Out is in col-major alignment with unit stride between rows.
      // That's the standard alignment---pass data as we made it.
      DGEMM(TransA, TransB, nRows, nCols, nLink, fFactor,
         pA, lda, pB, ldb, Beta, pOut, iColStO);
   } else {
      assert_rt(!"out.T path not tested.");
      // Out is in row-major alignment, with unit stride between columns,
      // not rows. To pass this, we need some transformations:
      //
      // - we have Out = A * B
      // -> Out^T = B^T * A^T
      //
      // So instead of Out with row major alignment, we pass
      // Out.T with standard col-major alignment and we compute
      // B^T * A^T.
      assert(iColStO == 1);
      assert(iRowStO != 0);
      assert(nRows == 0 || nCols == 0 || (pOut != pA && pOut != pB)); // input and output matrices cannot be aliased.
      char
         TransAt = InvertBlasTransposeArg(TransA),
         TransBt = InvertBlasTransposeArg(TransB);
      DGEMM(TransBt, TransAt, nCols, nRows, nLink, fFactor,
         pB, ldb, pA, lda, Beta, pOut, iRowStO);
   }
}



template<class FScalar>
inline FScalar Dot(FScalar const *a, FScalar const *b, size_t n)
{
   FScalar
      r = 0;
   for (size_t i = 0; i < n; ++i)
      r += a[i] * b[i];
   return r;
}

// r += f * x
template<class FScalar>
inline void Add(FScalar * RESTRICT r, FScalar const * RESTRICT x, FScalar f, size_t n)
{
   if (f != 1.0) {
      for (size_t i = 0; i < n; ++i)
         r[i] += f * x[i];
   } else {
      for (size_t i = 0; i < n; ++i)
         r[i] += x[i];
   }
}

// r *= f
template<class FScalar>
inline void Scale(FScalar *r, FScalar f, size_t n)
{
   for (size_t i = 0; i < n; ++i)
      r[i] *= f;
}


// r = f * x
template<class FScalar>
inline void Move(FScalar * RESTRICT r, FScalar const * RESTRICT x, FScalar f, size_t n)
{
   for (size_t i = 0; i < n; ++i)
      r[i] = f * x[i];
}


} // namespace ct

#endif // CX_ALGEBRA_H

// kate: indent-width 3
