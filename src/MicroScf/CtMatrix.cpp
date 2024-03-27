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

#include <stdexcept>
#include <cmath>
#include <fstream> // only for matrix import/export
#include <algorithm> // for std::sort

// #include "fmt.h"
#include "CxTypes.h"
#include "CtMatrix.h"
#include "CxIo.h"

#include "CxAlgebra.h" // fortran BLAS/LAPACK externals.


namespace ct {

inline void AssertSameShape(FMatrixView const &A, FMatrixView const &B)
{
   assert(A.nRows == B.nRows);
   assert(A.nCols == B.nCols);
   IR_SUPPRESS_UNUSED_WARNING(A);
   IR_SUPPRESS_UNUSED_WARNING(B);
}

void AssertSameShape1(FMatrixView const &A, FMatrixView const &B)
{
   AssertSameShape(A, B);
}

FMatrixView Transpose(FMatrixView const &In)
{
   return FMatrixView(In.pData, In.nCols, In.nRows, In.nColSt, In.nRowSt);
}

FMatrixView Select(FMatrixView const &In, size_t iRow, size_t iCol, size_t nRows, size_t nCols)
{
   assert(iRow + nRows <= In.nRows && iCol + nCols <= In.nCols);
   return FMatrixView(In.pData + iRow * In.nRowSt + iCol * In.nColSt,
      nRows, nCols, In.nRowSt, In.nColSt);
}

bool FMatrixView::IsSymmetric(FScalar Thresh) const
{
   if (nRows != nCols)
      return false;
   for (size_t iRow = 0; iRow < nRows; ++iRow)
      for (size_t iCol = 0; iCol < iRow; ++iCol)
         if (std::abs((*this)(iRow, iCol) - (*this)(iCol, iRow)) > Thresh)
            return false;
   return true;
}

double FMatrixView::fRmsdFromIdentity() const
{
   assert(nRows == nCols);
   if (nRows == 0)
      return 0.;
   double r = 0.;
   for (size_t iCol = 0; iCol < nCols; ++iCol) {
      double c = 0.;
      for (size_t iRow = 0; iRow < nRows; ++iRow) {
         double f = (*this)(iRow, iCol);
         if (iRow == iCol)
            f -= 1.;
         c += f * f;
      }
      r += c;
   }
   return std::sqrt(r) / double(nRows);
}

double FMatrixView::fRmsdFromSymmetry() const
{
   assert(nRows == nCols);
   double r = 0.;
   for (size_t iCol = 0; iCol < nCols; ++iCol) {
      double c = 0.;
      for (size_t iRow = 0; iRow < iCol; ++iRow) {
         double f = (*this)(iRow, iCol) - (*this)(iCol, iRow);
         c += f * f;
      }
      r += c;
   }
   return std::sqrt(r) / double(nRows);
}


// note: this was once in IvIao.cpp and once in CtFockBuild.cpp
// It can be replaced by M.fRmsdFromZero().
// double CalcMatrixRmsNorm(FMatrixView M)
// {
//    double d = 0;
//    for (size_t iCol = 0; iCol < M.nCols; ++ iCol)
//       for (size_t iRow = 0; iRow < M.nRows; ++ iRow)
//          d += sqr(M(iRow,iCol));
//    return std::sqrt(d);
// }



bool FMatrixView::IsDataContinuous() const
{
   return (nRowSt == 1 && (nColSt == nRows || nCols == 1)) ||
          (nColSt == 1 && (nRowSt == nCols || nRows == 1));
}

// element-wise operators on two matrices (FElementWiseOp2*).

// matrix operators (unary or binary) which do not accumulate anything and do not return anything.
struct FElementWiseOp_VoidResult {
   // these operators do not accumulate anything and do not return anything.
   typedef void FAccResult;
   FAccResult AccResult() const {};
};

// Dest[i,j] := Src[i,j]
struct FElementWiseOp2_Assign : public FElementWiseOp_VoidResult {
   inline void exec1(FScalar &Dest, FScalar Src) { Dest = Src; }
   inline void execN(FScalar *IR_RP pDest, FScalar const *IR_RP pSrc, size_t n) { memcpy(pDest, pSrc, n*sizeof(*pDest)); }
};

// Dest[i,j] := f*Src[i,j]
struct FElementWiseOp2_AssignScaled : public FElementWiseOp_VoidResult {
   FScalar Factor;
   FElementWiseOp2_AssignScaled(FScalar Factor_) : Factor(Factor_) {}
   inline void exec1(FScalar &Dest, FScalar Src) { Dest = Factor * Src; }
   inline void execN(FScalar *IR_RP pDest, FScalar const *IR_RP pSrc, size_t n) { for (size_t i = 0; i < n; ++ i) pDest[i] = Factor * pSrc[i]; }
};

// Dest[i,j] += Src[i,j]
struct FElementWiseOp2_Add : public FElementWiseOp_VoidResult  {
   inline void exec1(FScalar &Dest, FScalar Src) { Dest += Src; }
   inline void execN(FScalar *IR_RP pDest, FScalar const *IR_RP pSrc, size_t n) { for (size_t i = 0; i < n; ++ i) pDest[i] += pSrc[i]; }
};

// Dest[i,j] += f * Src[i,j]
struct FElementWiseOp2_AddScaled : public FElementWiseOp_VoidResult  {
   FScalar Factor;
   FElementWiseOp2_AddScaled(FScalar Factor_) : Factor(Factor_) {}
   inline void exec1(FScalar &Dest, FScalar Src) { Dest += Factor * Src; }
   inline void execN(FScalar *IR_RP pDest, FScalar const *IR_RP pSrc, size_t n) { for (size_t i = 0; i < n; ++ i) pDest[i] += Factor * pSrc[i]; }
};

// Dest[i,j] -= Src[i,j]
struct FElementWiseOp2_Subtract : public FElementWiseOp_VoidResult  {
   inline void exec1(FScalar &Dest, FScalar Src) { Dest -= Src; }
   inline void execN(FScalar *IR_RP pDest, FScalar const *IR_RP pSrc, size_t n) { for (size_t i = 0; i < n; ++ i) pDest[i] -= pSrc[i]; }
};

// matrix operators (unary or binary) which accumulate and return a single scalar value
struct FElementWiseOp_ScalarResult {
   typedef FScalar FAccResult;
   FAccResult m_AccResult;
   FElementWiseOp_ScalarResult() : m_AccResult(0.) {}
   FAccResult AccResult() const { return m_AccResult; }
};

// fAcc += A[i,j] * B[i,j]
struct FElementWiseOp2_Dot2 : public FElementWiseOp_ScalarResult {
   inline void exec1(FScalar Aij, FScalar Bij) {
      m_AccResult += Aij * Bij;
   }
   inline void execN(FScalar *IR_RP pA, FScalar const *IR_RP pB, size_t n) {
      for (size_t i = 0; i < n; ++ i)
         m_AccResult += pA[i] * pB[i];
   }
};


// execute element-wise operator Op(A,B), for two matrices Out, In (or A, B).
template<class FElementWiseOp2>
static typename FElementWiseOp2::FAccResult ExecElementWiseMatrixOp(FMatrixView Out, FMatrixView const &In, FElementWiseOp2 Op)
{
   AssertSameShape(Out, In);
   if (In.nRowSt == Out.nRowSt && In.nColSt == Out.nColSt && In.IsDataContinuous() && Out.IsDataContinuous()) {
      // both In and Out continuous memory blocks, and do both have the same
      // memory alignment? (either both col major or both row major)
      // -> process all elements in one loop
      Op.execN(Out.pData, In.pData, In.nRows * In.nCols);
   } else if (In.nRowSt == 1 && Out.nRowSt == 1) {
      // are columns of both In and Out continuous in memory?
      // -> execute one continuous array operation per column
      for (size_t iCol = 0; iCol != In.nCols; ++ iCol)
         Op.execN(Out.pData + Out.nColSt * iCol, In.pData + In.nColSt * iCol, In.nRows);
   } else if (In.nColSt == 1 && Out.nColSt == 1) {
      // are rows of both In and Out continuous in memory?
      // -> execute one continuous array operation per row
      for (size_t iRow = 0; iRow != In.nRows; ++ iRow)
         Op.execN(Out.pData + Out.nRowSt * iRow, In.pData + In.nRowSt * iRow, In.nCols);
   } else {
      // general case (e.g., inline transpose operation, or non-unit stride for
      // both axes)
      // -> execute operation one by one on all elements.
      for (size_t iCol = 0; iCol < Out.nCols; ++ iCol)
         for (size_t iRow = 0; iRow < Out.nRows; ++ iRow)
            Op.exec1(Out(iRow, iCol), In(iRow, iCol));
   }
   return Op.AccResult();
}


// Dest[i,j] := 0
struct FElementWiseOp1_Clear : public FElementWiseOp_VoidResult {
   inline void exec1(FScalar &Dest) { Dest = 0; }
//    inline void execN(FScalar *IR_RP pDest, size_t n) { memset(pDest, 0, n*sizeof(*pDest)); }
   inline void execN(FScalar *IR_RP pDest, size_t n) { for (size_t i = 0; i < n; ++ i) pDest[i] = 0; }
   // ^- I guess nowadays these types of things get translated to memset anyway?
};

// Dest[i,j] *= f
struct FElementWiseOp1_Scale : public FElementWiseOp_VoidResult {
   FScalar Factor;
   FElementWiseOp1_Scale(FScalar Factor_) : Factor(Factor_) {}
   inline void exec1(FScalar &Dest) { Dest *= Factor; }
   inline void execN(FScalar *IR_RP pDest, size_t n) { for (size_t i = 0; i < n; ++ i) pDest[i] *= Factor; }
};

// Acc += sqr(A[i,j])
struct FElementWiseOp1_AccNormSq : public FElementWiseOp_ScalarResult {
   inline void exec1(FScalar Aij) { m_AccResult += Aij*Aij; }
   inline void execN(FScalar *IR_RP pA, size_t n) { for (size_t i = 0; i < n; ++ i) exec1(pA[i]); }
};


// execute element-wise operator Op(A) on a single matrix A
template<class FElementWiseOp1>
static typename FElementWiseOp1::FAccResult ExecElementWiseMatrixOp(FMatrixView Out, FElementWiseOp1 Op)
{
   if (Out.IsDataContinuous()) {
      // -> process all elements in one loop
      Op.execN(Out.pData, Out.nRows * Out.nCols);
   } else if (Out.nRowSt == 1) {
      // are columns continuous in memory?
      // -> execute one continuous array operation per column
      for (size_t iCol = 0; iCol != Out.nCols; ++ iCol)
         Op.execN(Out.pData + Out.nColSt * iCol, Out.nRows);
   } else if (Out.nColSt == 1) {
      // are rows continuous in memory?
      // -> execute one continuous array operation per row
      for (size_t iRow = 0; iRow != Out.nRows; ++ iRow)
         Op.execN(Out.pData + Out.nRowSt * iRow, Out.nCols);
   } else {
      // general case (e.g., inline transpose operation, or non-unit stride for
      // both axes)
      // -> execute operation one by one on all elements.
      for (size_t iCol = 0; iCol < Out.nCols; ++ iCol)
         for (size_t iRow = 0; iRow < Out.nRows; ++ iRow)
            Op.exec1(Out(iRow, iCol));
   }
   return Op.AccResult();
}



// Out += f * In
void Add(FMatrixView Out, FMatrixView const &In, FScalar fFactor)
{
   if (fFactor == 1)
      ExecElementWiseMatrixOp(Out, In, FElementWiseOp2_Add());
   else if (fFactor == -1)
      ExecElementWiseMatrixOp(Out, In, FElementWiseOp2_Subtract());
   else
      ExecElementWiseMatrixOp(Out, In, FElementWiseOp2_AddScaled(fFactor));
}

// Out = f * In
void Assign(FMatrixView Out, FMatrixView const &In, FScalar fFactor)
{
   if (fFactor == 1) {
      if (Out.pData == In.pData) {
         // Out and In are aliased. The operation would not do anyhing.
         AssertSameShape(Out, In);
         assert_rt(Out.nRowSt == In.nRowSt && Out.nColSt == In.nColSt);
         return;
      }
      ExecElementWiseMatrixOp(Out, In, FElementWiseOp2_Assign());
   } else {
      ExecElementWiseMatrixOp(Out, In, FElementWiseOp2_AssignScaled(fFactor));
   }
}

// dot-product of two matrices. dot(A^T B) = Tr(A * B).
FScalar Dot(FMatrixView const &A, FMatrixView const &B)
{
   return ExecElementWiseMatrixOp(A, B, FElementWiseOp2_Dot2());
}


static double sqr(double f) { return f*f; }

FScalar Rmsd(FMatrixView const &A, FMatrixView const &B)
{
   AssertSameShape(A, B);
   double fSqDev = 0.;
   for (size_t iCol = 0; iCol < A.nCols; ++ iCol)
      for (size_t iRow = 0; iRow < A.nRows; ++ iRow)
         fSqDev += sqr(A(iRow,iCol) - B(iRow,iCol));
   size_t
      n = A.nRows * A.nCols;
   if (n == 0) {
      assert(fSqDev == 0.);
      return 0.;
   } else {
      return std::sqrt(fSqDev / double(n));
   }
}


FScalar Trace(FMatrixView const &In)
{
   assert(In.nRows == In.nCols);
   FScalar res = 0.0;
   for (size_t n = 0; n < In.nRows; ++n)
      res += In(n, n);
   return res;
}

// scale matrix by a factor in place
void Scale(FMatrixView &Out, FScalar fFactor)
{
   if (fFactor == 1.0)
      return;
   ExecElementWiseMatrixOp(Out, FElementWiseOp1_Scale(fFactor));
}

void FMatrixView::Clear()
{
   ExecElementWiseMatrixOp(*this, FElementWiseOp1_Clear());
}

void FMatrixView::SetIdentity()
{
   assert(nRows == nCols);
   Clear();
   for (size_t i = 0; i < nRows; ++i)
      (*this)(i, i) = 1.0f;
}

double FMatrixView::fRmsdFromZero(size_t nParamsForDivisor) const
{
   if (nRows * nCols == 0)
      return 0.;
   if (nParamsForDivisor == 0)
      nParamsForDivisor = nRows * nCols;
   double
      r = ExecElementWiseMatrixOp(*this, FElementWiseOp1_AccNormSq());
   return std::sqrt(r / (double(nParamsForDivisor)));
}







// Out = f * A * B
void Mxm(FMatrixView &Out, FMatrixView const &A, FMatrixView const &B,
    FScalar fFactor, uint Flags)
{
   // hm... if this is including CxAlgebra.h anyway, why don't we just defer the
   // call to there?
   assert(Out.nRows == A.nRows);
   assert(Out.nCols == B.nCols);
   assert(A.nCols == B.nRows); // == nLink

   if (Out.nRows == 0 || Out.nCols == 0)
      // output matrix is empty. Nothing to do.
      return;
   if (A.nCols == 0) {
      // no link dimension -- nothing to multiply.
      // That normally means that the update is a zero matrix.
      if (0 == (Flags & MXM_AddToDest))
         // output is a zero matrix and we are supposed to overwrite it.
         Out.Clear();
      return;
   }

   return Mxm_InlineCalls(Out.pData, Out.nRowSt, Out.nColSt,
         A.pData, A.nRowSt, A.nColSt,
         B.pData, B.nRowSt, B.nColSt,
         Out.nRows, A.nCols, Out.nCols,
         (0 != (Flags & MXM_AddToDest)), fFactor);
}

// Out += A * B, Out,B being vectors
void Mxva(FScalar *RESTRICT pOut, FMatrixView const &A, FScalar const *RESTRICT pIn)
{
   return Mxv(pOut, 1, A.pData, A.nRowSt, A.nColSt, pIn, 1, A.nRows, A.nCols, false, 1.0);
}

// Out += A * B, Out,B being vectors
void Mxvb(FScalar *RESTRICT pOut, FMatrixView const &A, FScalar const *RESTRICT pIn)
{
   return Mxv(pOut, 1, A.pData, A.nRowSt, A.nColSt, pIn, 1, A.nRows, A.nCols, true, 1.0);
}





// Out = A * B
void Mxm(FMatrixView &Out, FMatrixView const &A, FMatrixView const &B, uint Flags)
{
   return Mxm(Out, A, B, 1.0, Flags);
}

// Out = f * A^T * A for symmetric matrix Out.
void SyrkTN(FMatrixView &Out, FMatrixView const &A, FScalar fFactor, uint Flags)
{
//    return Mxm(Out, Transpose(A), A, fFactor, AddToDest);
   assert(Out.nRows == Out.nCols);
   assert(A.nCols == Out.nCols);

   assert(A.nRowSt == 1 || A.nColSt == 1);
   assert(Out.nRowSt == 1 || Out.nColSt == 1);
   // ^- otherwise dsyrk directly not applicable. Would need local copy
   // of matrix/matrices with compressed strides.

   if (A.nRows == 0) {
      if (0 == (Flags & MXM_AddToDest))
         Out.Clear();
      return;
   }

   double
      Beta = (0 != (Flags & MXM_AddToDest)) ? 1.0 : 0.0;
   char
      TransA = ((A.nRowSt == 1) ? 'T' : 'N');
   FINTARG
      lda = (A.nRowSt == 1) ? A.nColSt : A.nRowSt,
      ldc = (Out.nRowSt == 1) ? Out.nColSt : Out.nRowSt;
   // note: since Out is symmetric, just exchanging ldc here in case
   // Out.nRowSt != 1 is sufficient (unlike for general Mxm)
   DSYRK('L', TransA, Out.nRows, A.nRows, fFactor, A.pData, lda, Beta, Out.pData, ldc);
   // fix up the upper triangle of the matrix.
   for (size_t iCol = 0; iCol < Out.nCols; ++iCol)
      for (size_t iRow = 0; iRow < iCol; ++iRow)
         Out(iRow, iCol) = Out(iCol, iRow);
}

// Out = f * A * A^T for symmetric matrix Out.
void SyrkNT(FMatrixView &Out, FMatrixView const &A, FScalar fFactor, uint Flags)
{
   return SyrkTN(Out, Transpose(A), fFactor, Flags);
}




void Diagonalize(FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack & /*Mem*/)
{
   assert_rt(InOut.nRowSt == 1);
   assert(InOut.nRows == InOut.nCols);
   Diagonalize(pEigenValues, InOut.pData, InOut.nColSt, InOut.nRows);
   // ^- in CxAlgebra.cpp
}

// As Diagonalize(), but eigenvalues are returned in descending order.
void Diagonalize_LargeEwFirst(FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem)
{
   Scale(InOut, -1.0);
   Diagonalize(InOut, pEigenValues, Mem);
   Scale(pEigenValues, -1.0, InOut.nCols);
}


void ComputeSvd(FMatrixView &U, FScalar *pSigma, FMatrixView &Vt, FMatrixView In, FMemoryStack &Mem)
{
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   size_t
      nRows = In.nRows,
      nCols = In.nCols,
      nSig = std::min(nRows, nCols);
   assert(U.nRows >= nRows && U.nCols >= nSig);
   assert(Vt.nRows >= nSig && Vt.nCols >= nCols);
   U.nRows = nRows;
   U.nCols = nSig;
   Vt.nRows = nSig;
   Vt.nCols = nCols;
   assert_rt(U.nRowSt == 1 && Vt.nRowSt == 1);
   FMatrixView
      Tmp = MakeStackMatrix(nRows, nCols, Mem);
   Assign(Tmp, In); // to get rid of strides. Also, DGESDD destroys the input matrix.
   ComputeSvd(U.pData, U.nColSt, pSigma, Vt.pData, Vt.nColSt, Tmp.pData, Tmp.nColSt, nRows, nCols);
}


// As ComputeSvd, but allocate output quantities (U, pSigma, and Vt) on Mem.
void AllocAndComputeSvd(FMatrixView &U, FScalar *&pSigma, FMatrixView &Vt, FMatrixView In, FMemoryStack &Mem)
{
   size_t nRows = In.nRows, nCols = In.nCols, nSig = std::min(nRows, nCols);
   U = MakeStackMatrix(nRows, nSig, Mem);
   Vt = MakeStackMatrix(nSig, nCols, Mem);
   Mem.Alloc(pSigma, nSig);
   ComputeSvd(U, pSigma, Vt, In, Mem);
}

void AllocAndComputeSvd(FMatrixView &U, FMatrixView &Sigma, FMatrixView &Vt, FMatrixView In, FMemoryStack &Mem)
{
   double *pSigma;
   AllocAndComputeSvd(U, pSigma, Vt, In, Mem);
   Sigma = FMatrixView(pSigma, U.nCols, 1);
}


void FMatrixView::Print(std::ostream &out, std::string const &Name) const
{
   PrintMatrixGen(out, pData, nRows, nRowSt, nCols, nColSt, Name);
}

void Assign(FStackMatrix &Out, FMatrixView const &In, FScalar fFactor)
{
   Out.Reshape(In.nRows, In.nCols);
   Assign(*static_cast<FMatrixView *>(&Out), In, fFactor);
}

void Assign(FHeapMatrix &Out, FMatrixView const &In, FScalar fFactor)
{
   Out.Reshape(In.nRows, In.nCols);
   Assign(*static_cast<FMatrixView *>(&Out), In, fFactor);
}

void Mxm(FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, uint Flags)
{
   assert((0 == (Flags & MXM_AddToDest)) || (Out.nRows == A.nRows && Out.nCols == B.nCols));
   Out.Reshape(A.nRows, B.nCols);
   Mxm(*static_cast<FMatrixView *>(&Out), A, B, Flags);
}

void Mxm(FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags)
{
   assert((0 == (Flags & MXM_AddToDest)) || (Out.nRows == A.nRows && Out.nCols == B.nCols));
   Out.Reshape(A.nRows, B.nCols);
   Mxm(*static_cast<FMatrixView *>(&Out), A, B, fFactor, Flags);
}

std::size_t FMatrixView::TriangularStorageSize(int Sign) const
{
   assert(nRows == nCols);
   assert(Sign == +1 || Sign == -1);
//    if (Sign == 1)
//       return (nRows*(nRows+1))/2;
//    else
//       return (nRows*(nRows-1))/2;
   // ^- Molpro does not actually omit storing the diagonal of antisymmetric
   //    matrices. By doing it also this way we can keep binary compatibility.
   //    The current TriangularExpand/Reduce functions also assume the diagonal
   //    elements to be present even in the Sign=-1 case.
   return (nRows*(nRows+1))/2;
   IR_SUPPRESS_UNUSED_WARNING(Sign);
}


void FMatrixView::TriangularExpand(int Sign, FScalar const *pTriangData)
{
   assert(nRowSt == 1 && nColSt >= nRows); // <- required?
   assert(nRows == nCols);
   if (pTriangData == 0)
      pTriangData = pData; // in-place.
   size_t
      iOff = (nRows * (nRows + 1)) / 2;
   size_t
      nRowSt_ = nRowSt,
      nColSt_ = nColSt;
   assert(nRowSt == 1 || nColSt == 1);
   if (nRowSt_ > nColSt_)
      std::swap(nRowSt_, nColSt_);
   // unpack explicitly stored triangle
   for (size_t iCol = nCols - 1; iCol < nCols; --iCol)
      for (size_t iRow = iCol; iRow <= iCol; --iRow) {
         iOff -= 1;
         pData[iRow * nRowSt_ + iCol * nColSt_] = pTriangData[iOff];
      }
   assert(iOff == 0);
   // copy restored elements to other triangle
   for (size_t iCol = 0; iCol != nCols; ++iCol)
      for (size_t iRow = 0; iRow < iCol; ++iRow)
         // ^- note: assumes C++ standard-conformant underflow behavior of unsigned types
         pData[iCol * nRowSt_ + iRow * nColSt_] = Sign * pData[iCol * nColSt_ + iRow * nRowSt_];
   ClearEmptySpace();
}

void FMatrixView::TriangularReduce(int Sign, FScalar *pTriangData, FScalar fOffDiagFactor, bool TransposeSource)
{
   assert(nRowSt == 1 && nColSt >= nRows); // <- required?
   assert(nRows == nCols);
   if (pTriangData == 0)
      pTriangData = pData; // in-place.
   size_t
      iOff = 0;
   size_t
      nRowSt_ = nRowSt,
      nColSt_ = nColSt;
   assert(nRowSt == 1 || nColSt == 1);
   if (nRowSt_ > nColSt_)
      std::swap(nRowSt_, nColSt_);
#ifdef _DEBUG
   for (size_t iCol = 0; iCol < nCols; ++iCol)
      for (size_t iRow = 0; iRow <= iCol; ++iRow) {
         assert(std::abs(pData[iRow * nRowSt_ + iCol * nColSt_] -
            (double)Sign * pData[iCol * nRowSt_ + iRow * nColSt_]) < 1e-10);
      }
   iOff = 0;
#endif
   if (!TransposeSource) {
      for (size_t iCol = 0; iCol < nCols; ++iCol)
         for (size_t iRow = 0; iRow <= iCol; ++iRow) {
            // xout << fmt::format("pack: T[{:3}] <- Q[{:3}]  iRow={:2}  iCol={:2}  F={:f}", iOff, (iRow *
            //             nRowSt_ + iCol * nColSt_), iRow, iCol, pData[iRow * nRowSt_ + iCol * nColSt_]) << std::endl;
            double dij = pData[iRow * nRowSt_ + iCol * nColSt_];
            if (fOffDiagFactor != 1. && iRow != iCol)
               dij *= fOffDiagFactor;
            pTriangData[iOff] = dij;
            iOff += 1;
         }
   } else {
      for (size_t iCol = 0; iCol < nCols; ++iCol)
         for (size_t iRow = 0; iRow <= iCol; ++iRow) {
            double dij = pData[iRow * nColSt_ + iCol * nRowSt_];
            if (fOffDiagFactor != 1. && iRow != iCol)
               dij *= fOffDiagFactor;
            pTriangData[iOff] = dij;
            iOff += 1;
         }
   }
   // PrintMatrixGen(xout, pTriangData, 1, 1, iOff, 1, "TRIANGULAR DATA AFTER PACKING:");
   IR_SUPPRESS_UNUSED_WARNING(Sign);
}

void FMatrixView::TriangularExpandOrReduce(int Sign, FScalar *pTriangData, FTriangularTransformFlags Direction)
{
   assert(Direction == TRIANG_Expand || Direction == TRIANG_Reduce);
   if (Direction == TRIANG_Expand)
      return TriangularExpand(Sign, pTriangData);
   else
      return TriangularReduce(Sign, pTriangData);
}

// executes Out += Factor * TriangularExpand(Sign,pDataTr,fOffDiagFactor), where
// TriangularExpand(Mat,Sign,fOffDiagFactor) denotes the square-expanded version
// of packed triangular matrix with Sign +1 (symmetric) or Sign -1 (anti-symmetric).
//
// WARNING: unlike the FMatrixView.TriangularExpand function, this one
// *DOES NOT* allow pMatTr and Out.pData to overlap!
void Add_TriangularExpanded(FMatrixView Out, int Sign, double const *pDataTr, double fFactor, double fOffDiagFactor)
{
   assert(Out.nRows == Out.nCols);
   assert(pDataTr != Out.pData);
   assert(Sign == +1 || Sign == -1);
   size_t
      nRows = Out.nRows;
   fOffDiagFactor *= fFactor;
   for (size_t iCol = 0; iCol < nRows; ++ iCol) {
      size_t
         iColOff = (iCol*(iCol+1))/2;
      for (size_t iRow = 0; iRow < iCol; ++ iRow) {
         // off-diagonal elements.
         assert(iRow < iCol);
         double
            fData = fOffDiagFactor * pDataTr[iRow + iColOff];
         Out(iRow,iCol) += fData;
         Out(iCol,iRow) += Sign*fData;
      }
      // process diagonal element -- in symmetric case only
      // (in the anti-symmetric case it is zero, so there is nothing to increment the target with)
      if (Sign == +1)
         Out(iCol,iCol) += fFactor * pDataTr[iColOff + iCol];
   }
}


void FMatrixView::ClearEmptySpace()
{
   if (nRowSt == 1) {
      if (nColSt == nRows)
         return;
      for (size_t iCol = 0; iCol < nCols; ++iCol)
         for (size_t iEmptyRow = nRows; iEmptyRow != nColSt; ++iEmptyRow)
            pData[iEmptyRow + iCol * nColSt] = 0;

   } else {
      assert(nColSt == 1);
      if (nRowSt == nCols)
         return;
      for (size_t iRow = 0; iRow < nRows; ++iRow)
         for (size_t iEmptyCol = nCols; iEmptyCol != nRowSt; ++iEmptyCol)
            pData[iEmptyCol + iRow * nRowSt] = 0;
   }
}



void CalcSmhMatrix(FMatrixView M, FMemoryStack &Mem, FSmhOptions const &Opt) {
   return CalcSymMatrixPower(M, -0.5, Mem, Opt);
}

void CalcSymMatrixPower(FMatrixView M, double fExponent, FMemoryStack &Mem, FSmhOptions const &Opt)
{
   assert_rt(M.nRows == M.nCols);
   if (M.nRows == 0)
       return;
   size_t
       N = M.nRows;
   FStackMatrix
       T1(N, N, &Mem),
       T2(N, N, &Mem);
   Assign(T1, M); // should get rid of possible row strides of M (diag2 doesn't like them).
   TMemoryLock<double>
      // eigenvalues
      pEw(N, &Mem);

   Diagonalize(T1, pEw, Mem);

   if (false && Opt.pXout && Opt.pDelMsg) {
      std::ostream &xout = *Opt.pXout;
      xout << "Spectrum for " << Opt.pDelMsg << "\n    ";
      for (size_t i = 0; i < N; ++ i)
         xout << fmt::format("  {:10.3e}", pEw[i]);
      xout << std::endl;
   }
   if (Opt.pXout && Opt.pDelMsg) {
      std::ostream &xout = *Opt.pXout;
      xout << "" << Opt.pDelMsg << ": ";
      size_t
         nLen = 0;
      for (char const *p = Opt.pDelMsg; *p; ++p)
         nLen += 1;
      for (size_t i = 27; i > nLen; -- i)
         xout << " ";
      if (N == 0) {
         xout << "EwNone Size Zero Matrix";
      } else {
         xout << fmt::format("  EwMin ={:10.3e}  EwMax={:10.3e}e  Ratio={:10.3e}",
                  pEw[0], pEw[N-1], (pEw[0]/pEw[N-1]))
              << std::endl;
      }
   }
   size_t
      nDel = 0;
   bool
      // Actually negative eigenvalues will be ignored indiscriminately.
      // However, zero and very-small eigenvalues are only a problem
      // if we mean to compute a *negative* power of the matrix. So
      // for positive exponents we do not apply the thresholds.
      SkipSmallEw = (fExponent <= 0.);
   for (size_t iCol = 0; iCol < N; ++ iCol) {
      double
         Ew = pEw[iCol],
         f = 0;
      if (Ew < 0. || (SkipSmallEw && (Ew < Opt.ThrAbs || Ew < Opt.ThrRel * pEw[N-1]))) {
         f = 0;
         nDel += 1;
      } else {
         if (fExponent == -0.5)
            // make a Smh = S^{-1/2} matrix.
            f = 1.0/std::sqrt(std::sqrt(Ew));
         else if (fExponent == +0.5)
            // make a Sh = S^{+1/2} matrix.
            f = std::sqrt(std::sqrt(Ew));
         else
            // make a different matrix power
            f = std::pow(Ew, .5*fExponent);
      }
//       f = std::pow(Ew, .5*fExponent);
      for (size_t iRow = 0; iRow < N; ++ iRow) {
         T2(iRow,N-iCol-1) = f * T1(iRow,iCol);
         // ^- invert order: smallest coeffs first.
         // (that's what setting T1 = -M; ... ew *= -1 was for originally)
      }
   }

   if (nDel != 0 && Opt.pDelMsg) {
      std::ostream
         &out = *(Opt.pXout? Opt.pXout : &xout);
      out << "Warning: Deleted " << nDel
          << " vectors while constructing " << Opt.pDelMsg << "." << std::endl;
   }

   // = Mxm(M, T2, Transpose(T2));
   SyrkNT(M, T2);
}



#if 0
// solve linear least squares ||M * x - b||_2 -> min subject to ||x||_2 ->min
// using singular value decomposition of A. Singular values below fThrRel are
// treated as zero. A negative value indicates the usage of machine precision.
// Both M and RhsSol are overwritten.
void LeastSquaresSolveUnsafe(FMatrixView M, FMatrixView &RhsSol, double fThrRel, FMemoryStack &Mem)
{
   using std::min;
   using std::max;
   assert(M.nRowSt == 1);
   assert(RhsSol.nRowSt == 1);
   assert(RhsSol.nRows == M.nRows); // input: rhs. M.nRows x nRhs
   assert(RhsSol.nColSt >= std::max(M.nRows, M.nCols));
   // ^- output also goes here. must have enough space for this.

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
      nWork = 3*min(M.nRows,M.nCols) + max(max(2*min(M.nRows, M.nCols), max(M.nRows, M.nCols)), RhsSol.nCols);
   double
      fWork = 0,
      *pWork,
      *pSig; // singular values (decreasing order).
   Mem.Alloc(pSig, min(M.nRows, M.nCols));
   DGELSS(M.nRows, M.nCols, RhsSol.nCols, M.pData, M.nColSt,
      RhsSol.pData, RhsSol.nColSt, pSig, fThrRel, &nRank,
      &fWork, -1, &info);
   if (info != 0) throw std::runtime_error("dgelss failed in workspace query.");
   nWork = static_cast<FORTINT>(fWork + 0.5);
   Mem.Alloc(pWork, nWork);
   DGELSS(M.nRows, M.nCols, RhsSol.nCols, M.pData, M.nColSt,
      RhsSol.pData, RhsSol.nColSt, pSig, fThrRel, &nRank,
      pWork, nWork, &info);
   if (info != 0) throw std::runtime_error("dgelss failed.");
   RhsSol.nRows = M.nCols; // output: solution. M.nCols x nRhs
   Mem.Free(pWork);
   Mem.Free(pSig);
}

void LeastSquaresSolveSafe(FMatrixView const M, FMatrixView Sol, FMatrixView const Rhs, double fThrRel, FMemoryStack &Mem)
{
   assert(M.nRows == Rhs.nRows && M.nCols == Sol.nRows && Sol.nCols == Rhs.nCols);
   size_t
      ldRhsSol = std::max(Rhs.nRows, Sol.nRows);
   FStackMatrix
      MCopy(M.nRows, M.nCols, &Mem),
      RhsSol(ldRhsSol, Rhs.nCols, &Mem);
   RhsSol.nRows = Rhs.nRows;
   Assign(*(FMatrixView *)&RhsSol, Rhs);
   Assign(MCopy, M);
   //     MCopy.Print(xout, "MCopy/IN");
   //     RhsSol.Print(xout, "RhsSol/IN");
   LeastSquaresSolveUnsafe(MCopy, RhsSol, fThrRel, Mem);
   //     RhsSol.Print(xout, "RhsSol/OUT");
   Assign(Sol, RhsSol);
};
#else
void LeastSquaresSolveSafe(FMatrixView const M, FMatrixView Sol, FMatrixView const Rhs, double fThrRel, FMemoryStack &Mem)
{
   assert_rt(M.nRows == Rhs.nRows && M.nCols == Sol.nRows && Sol.nCols == Rhs.nCols);
   assert_rt(M.nRowSt == 1 && Sol.nRowSt == 1 && Rhs.nRowSt == 1);
   size_t
      nRank = LeastSquaresSolve(Sol.pData, Sol.nColSt, M.pData, M.nColSt, Rhs.pData, Rhs.nColSt, M.nRows, M.nCols, Rhs.nCols, fThrRel);
   (void)nRank;
   (void)Mem; // suppress unused warning.
};
#endif // 0






// schmidt-orthogonalize the nAo x nOcc orbitals in C with respect to the nAo x nAo overlap matrix S1.
void OrthSchmidt(FMatrixView Orbs, FMatrixView S1, FMemoryStack &Mem)
{
   size_t
      nOcc = Orbs.nCols;
   FStackMatrix
      SOrb(nOcc, nOcc, &Mem);
   ChainMxm(SOrb, Transpose(Orbs), S1, Orbs, Mem);
   CalcCholeskyFactors(SOrb);
   TriangularSolve(Transpose(Orbs), SOrb);
}


// schmidt-orthogonalize the nAo x nOcc orbitals in C with respect to the nAo x nAo overlap matrix S1,
// which is implicitly given by its lower Cholesky decomposition: S1 = Scd1 * Scd1.T
void OrthSchmidtScd(FMatrixView C, FMatrixView Scd1, FMemoryStack &Mem)
{
   size_t
      nAo = C.nRows,
      nOcc = C.nCols;
   FStackMatrix
      SOrb(nOcc, nOcc, &Mem),
      T1(nAo, nOcc, &Mem);
   assert(nAo == Scd1.nRows && nAo == Scd1.nCols);
   Assign(T1, C);
   TriangularMxm(T1, Scd1, 'L', 'T');
   SyrkTN(SOrb, T1);
   CalcCholeskyFactors(SOrb);
   TriangularSolve(Transpose(C), SOrb, 'L', 'N');
}



// symmetrically orthogonalize a set of vectors with respect to overlap matrix S.
void SymOrth(FMatrixView Orbs, FMatrixView const S, FMemoryStack &Mem, double fThrAbs, double fThrRel)
{
   FStackMatrix
      SmhOrb(Orbs.nCols, Orbs.nCols, &Mem),
      // temp to hold result of Orbs * Smh
      T1(Orbs.nRows, Orbs.nCols, &Mem);
   ChainMxm(SmhOrb, Transpose(Orbs), S, Orbs, Mem);
   CalcSmhMatrix(SmhOrb, Mem, FSmhOptions(fThrAbs,fThrRel,0,0));
   Mxm(T1, Orbs, SmhOrb);
   Assign(Orbs, T1);
}


// calculate Cholesky factorization M = L * L^T
void CalcCholeskyFactors(FMatrixView M)
{
//    assert_rt(M.nRowSt == 1 && M.nRows == M.nCols && M.nRowSt == 1 && M.nColSt >= M.nRows && bool("void CalcCholeskyFactors(FMatrixView M)"));
   assert_rt(M.nRowSt == 1 && M.nRows == M.nCols && (M.nRowSt == 1 || M.nColSt == 1) && bool("void CalcCholeskyFactors(FMatrixView M)"));
   if (M.nRows == 0)
      return; // nothing to do...
   bool
      TransposedL = false;
   if (M.nColSt == 1 && M.nRowSt > 1) {
      // if the input matrix does not have default memory alignment (nRowSt ==
      // 1, and nColSt >= nRows), we can still Cholesky-decompose it.
      //
      // First, we transpose M by flipping the row and column strides. As M is
      // supposed to be symmetric on input, this has no further consequence.
      //
      // However, we do need to formally flip the CD type from L L.T to U.T U
      // because 'M' is also used for referencing the output matrix, and we just
      // transposed it: so we're making the U.T U decomposition with U
      // implicitly standing for L.T at the same storage location..
      assert(M.nRows == M.nCols);
      std::swap(M.nColSt, M.nRowSt);
      TransposedL = true;
      assert(!"transposed case of CalcCholeskyFactors not yet tested");
   };
   assert(M.IsSymmetric());
//    FORTINT
//       info = 0;
//    // factorize M into M = L * L^T.
//    DPOTRF(TransposedL? 'U' : 'L', M.nRows, M.pData, M.nColSt, &info);
//    if (info != 0) {
//       throw std::runtime_error(fmt::format("CalcCholeskyMatrix: dpotrf failed on input matrix of shape {} with strides [{}, {}]. LAPACK returned info code {}.", M.Shape(), M.nRowSt, M.nColSt, info));
//    }
   CalcCholeskyFactors(M.pData, M.nColSt, M.nRows, TransposedL);

   // restore the original matrix order.
   if (TransposedL) {
      std::swap(M.nColSt, M.nRowSt);
      TransposedL = false;
   }

   // clear out strict upper triangle. Neither DPOTRF nor
   // TriangularMxm/TriangularSolve reference it, but clearing it will allow us
   // to use L like a reguar matrix if needed (e.g., to re-assemble the original
   // M with a SyrkNT)
   for (size_t i = 0; i < M.nRows; ++i)
      for (size_t j = i + 1; j < M.nRows; ++j)
         M(i, j) = 0;
}

// invert lower triangular matrix L. Upper part ignored.
void InvertTriangularMatrix(FMatrixView L)
{
   assert(L.nRowSt == 1 && L.nRows == L.nCols);
//    FORTINT
//       info = 0;
//    DTRTRI('L', 'N', L.nRows, L.pData, L.nColSt, &info);
//    if (info != 0)
//       throw std::runtime_error("InvertTriangularMatrix: dtrtri failed.");
   InvertTriangularMatrix(L.pData, L.nColSt, L.nRows);
}


void TriangularSolve(FMatrixView RhsSol, FMatrixView const &L, char Side, char TransL, double fFactor)
{
   assert(L.nRowSt == 1 && L.nRows == L.nCols);
   assert(Side == 'L' || Side == 'R');
   assert(TransL == 'N' || TransL == 'T');
   if (L.nColSt == 0) {
      assert(L.nRows == 0);
      return;
   }
   FORTINT
   info = 0;
   char RevTransL = (TransL == 'T')? 'N' : 'T'; // <- reverse of of TransL
   if (Side == 'L') {
      assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
      assert(L.nRows == RhsSol.nRows);
      if (RhsSol.nRowSt == 1)
         DTRSM('L', 'L', TransL, 'N', RhsSol.nRows, RhsSol.nCols, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt, &info);
      //         DTRSM('L', 'L', 'N', 'N', RhsSol.nRows, RhsSol.nCols, 1.0, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt, &info);
      else
         // A X^T == B^T
         // <-> X A^T == B
         DTRSM('R', 'L', RevTransL, 'N', RhsSol.nCols, RhsSol.nRows, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt, &info);
      //         DTRSM('R', 'L', 'T', 'N', RhsSol.nCols, RhsSol.nRows, 1.0, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt, &info);
   } else {
      assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
      assert(L.nRows == RhsSol.nCols);

      if (RhsSol.nRowSt == 1)
         DTRSM('R', 'L', TransL, 'N', RhsSol.nRows, RhsSol.nCols, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt, &info);
      else
         // A X^T == B^T
         // <-> X A^T == B
         DTRSM('L', 'L', RevTransL, 'N', RhsSol.nCols, RhsSol.nRows, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt, &info);
   }
   if (info != 0)
      throw std::runtime_error("TriangularSolve: DTRSM failed.");
}


// set A := L * A where L is a lower triangular matrix (LR = 'L').
// or  A := A * L where L is a lower triangular matrix (LR = 'R')
void TriangularMxm(FMatrixView RhsSol, FMatrixView const &L, char Side, char TransL, double fFactor)
{
   assert(L.nRowSt == 1 && L.nRows == L.nCols);
   assert(Side == 'L' || Side == 'R');
   assert(TransL == 'N' || TransL == 'T');
   if (L.nColSt == 0) {
      assert(L.nRows == 0);
      return;
   }
   char RevTransL = (TransL == 'T')? 'N' : 'T'; // <- reverse of of TransL
   if (Side == 'L') {
      assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
      assert(L.nRows == RhsSol.nRows);
      if (RhsSol.nRowSt == 1)
         DTRMM('L', 'L', TransL, 'N', RhsSol.nRows, RhsSol.nCols, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt);
      else
         // A X^T == B^T
         // <-> X A^T == B
         DTRMM('R', 'L', RevTransL, 'N', RhsSol.nCols, RhsSol.nRows, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt);
   } else {
      assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
      assert(L.nRows == RhsSol.nCols);

      if (RhsSol.nRowSt == 1)
         DTRMM('R', 'L', TransL, 'N', RhsSol.nRows, RhsSol.nCols, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt);
      else
         // A X^T == B^T
         // <-> X A^T == B
         DTRMM('L', 'L', RevTransL, 'N', RhsSol.nCols, RhsSol.nRows, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt);
   }
}



// Solve A X = B, where A = L*L^T and L is a lower triangular matrix.
void CholeskySolve(FMatrixView RhsSol, FMatrixView const &L)
{
   // A X = B
   // -> L L^T X = B
   // -> L (L^T X) = B
   // so get X1 = solve(L, B), where X1 = L^T X
   TriangularSolve(RhsSol, L, 'L', 'N', 1.0);
   // then get X = solve(L^T, X1)
   TriangularSolve(RhsSol, L, 'L', 'T', 1.0);
}

// Multiply X = A * B, where A = L*L^T and L is a lower triangular matrix.
void CholeskyMxm(FMatrixView RhsSol, FMatrixView const &L)
{
   // X := A * B = (L * L^T) * B
   // so first left ('L')-multiply L^T to get X1 := L^T B
   TriangularMxm(RhsSol, L, 'L', 'T', 1.0);
   // then left-multyiply with L itself to get X = L X1 = L (L^T B) = A B
   TriangularMxm(RhsSol, L, 'L', 'N', 1.0);
}


void CholeskySolve(double *pRhsAndSolution, FMatrixView const &L)
{
   CholeskySolve(FMatrixView(pRhsAndSolution, L.nRows, 1), L);
}

void CholeskyMxm(double *pRhsAndSolution, FMatrixView const &L)
{
   return CholeskyMxm(FMatrixView(pRhsAndSolution, L.nRows, 1), L);
}


void Symmetrize(FMatrixView M, double Phase)
{
   assert(M.IsSquare());
   for (size_t i = 0; i < M.nRows; ++i) {
      for (size_t j = 0; j <= i; ++j) {
         FScalar f = .5*(M(i,j) + Phase * M(j,i));
         M(i,j) = f;
         M(j,i) = Phase * f;
      }
   }
}


FBasicPerm Inverted(FBasicPerm const &in)
{
   size_t
      N = in.size();
   size_t const
      iInvalid = N;
   FBasicPerm
      out(N, iInvalid); // initialize all entries to the (invalid) index `N`
   bool
      brrked = false;
   for (size_t i = 0; i < N; ++ i) {
      size_t j = in[i];
      if (size_t(j) >= size_t(N) || out[j] != iInvalid) {
         // invalid permutation: either in[i] is not in the valid range {0,1,..,N-1},
         // or the value of in[i] already occurred before (we initialized all entries
         // in out[:] to iInvalid == N; so if out[in[i]] is NOT `iInvalid` anymore,
         // it means we already overwrote the entry for a lower `i` iteration.
         brrked = true;
         break;
      }
      out[j] = i;
   }
   if (brrked)
      // note: that above is the only check we need to perform. If the loop goes
      // through then the output permutation is automatically valid.
      throw std::runtime_error(fmt::format("Inverted(FBasicPerm const&): input permutation of length {} is invalid (contains repeated indices or indices outside the allowed set {0, 1,..., in.size()-1}).", in.size()));
   return out;
}


FBasicPerm ArgSort(double const *pVals, size_t iValSt, size_t nVals, bool Reverse)
{
   FBasicPerm
      out(nVals, nVals); // initialize all entries to the (invalid) index `nVals`
   ArgSort1(&out[0], pVals, iValSt, nVals, Reverse);
   return out;
}


void ArgSort(FBasicPerm &Indices, double const *pVals, size_t iValSt, size_t nVals, bool Reverse)
{
   Indices.assign(nVals, nVals);
   ArgSort1(&Indices[0], pVals, iValSt, nVals, Reverse);
}



// permute rows of matrix such that
//    NewM[iRow,iCol] = OldM[P[iRow], iCol]
void PermuteRows(FMatrixView M, size_t const *P, FMemoryStack &Mem)
{
   for (size_t iCol = 0; iCol < M.nCols; ++iCol) {
      double *p;
      Mem.Alloc(p, M.nRows);
      for (size_t iRow = 0; iRow < M.nRows; ++iRow)
         p[iRow] = M(P[iRow], iCol);
      for (size_t iRow = 0; iRow < M.nRows; ++iRow)
         M(iRow, iCol) = p[iRow];
      Mem.Free(p);
   }
}

void Permute(FMatrixView &M, size_t const *pRowPerm, size_t const *pColPerm, FMemoryStack &Mem)
{
   PermuteRows(M, pRowPerm, Mem);
   PermuteRows(Transpose(M), pColPerm, Mem);
}



// chained matrix product: Out = pM[0] x pM[1] x ...
// WARNING: this function modifies pM in place!
void ChainMxmR(FMatrixView &Out, FMatrixView **pM, size_t nM, FMemoryStack &Mem, double fFactor, uint Flags)
{
   assert(nM > 0);
   if (nM == 1) {
      if (Flags & MXM_AddToDest)
         return Add(Out, *pM[0], fFactor);
      else
         return Assign(Out, *pM[0], fFactor);
   }
   if (nM == 2)
      return Mxm(Out, *pM[0], *pM[1], fFactor, Flags);

   // more than two matrices: check which matrices to multiply first. we seek to
   // find the smallest intermediate (that's a heuristic).
   //
   // (note: there are exact solutions for this problem, in small time
   // complexity; but I do not need this often, and the heuristic tends to work
   // well enough. Should still be revisited some time)
   assert(nM < 16);
   std::size_t nSzMin = pM[0]->nRows * pM[1]->nCols, iMin = 0;
   for (size_t i = 1; i < nM-1; ++ i) {
      std::size_t
        nSz = pM[i]->nRows * pM[i+1]->nCols;
      if (nSz < nSzMin) {
         nSzMin = nSz;
         iMin = i;
      }
   }

   FStackMatrix
      MI(pM[iMin]->nRows, pM[iMin+1]->nCols, &Mem);
   Mxm(MI, *pM[iMin], *pM[iMin+1]);

   // remove pM[iMin] and pM[iMin+1] from the chain and replace by their product.
   pM[iMin] = &MI;
   for (size_t i = iMin+1; i < nM-1; ++ i)
      pM[i] = pM[i+1];

   return ChainMxmR(Out, pM, nM-1, Mem, fFactor, Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMemoryStack &Mem, uint Flags) {
   FMatrixView *M[1] = {&A0};
   return ChainMxmR(Out, &M[0], 1, Mem, 1., Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMemoryStack &Mem, uint Flags) {
   FMatrixView *M[2] = {&A0, &A1};
   return ChainMxmR(Out, &M[0], 2, Mem, 1., Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMemoryStack &Mem, uint Flags) {
   FMatrixView *M[3] = {&A0, &A1, &A2};
   return ChainMxmR(Out, &M[0], 3, Mem, 1., Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMemoryStack &Mem, uint Flags) {
   FMatrixView *M[4] = {&A0, &A1, &A2, &A3};
   return ChainMxmR(Out, &M[0], 4, Mem, 1., Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMemoryStack &Mem, uint Flags) {
   FMatrixView *M[5] = {&A0, &A1, &A2, &A3, &A4};
   return ChainMxmR(Out, &M[0], 5, Mem, 1., Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMatrixView A5, FMemoryStack &Mem, uint Flags) {
   FMatrixView *M[6] = {&A0, &A1, &A2, &A3, &A4, &A5};
   return ChainMxmR(Out, &M[0], 6, Mem, 1., Flags);
}


void ChainMxm(FMatrixView Out, FMatrixView A0, FMemoryStack &Mem, double fFactor, uint Flags) {
   FMatrixView *M[1] = {&A0};
   return ChainMxmR(Out, &M[0], 1, Mem, fFactor, Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMemoryStack &Mem, double fFactor, uint Flags) {
   FMatrixView *M[2] = {&A0, &A1};
   return ChainMxmR(Out, &M[0], 2, Mem, fFactor, Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMemoryStack &Mem, double fFactor, uint Flags) {
   FMatrixView *M[3] = {&A0, &A1, &A2};
   return ChainMxmR(Out, &M[0], 3, Mem, fFactor, Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMemoryStack &Mem, double fFactor, uint Flags) {
   FMatrixView *M[4] = {&A0, &A1, &A2, &A3};
   return ChainMxmR(Out, &M[0], 4, Mem, fFactor, Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMemoryStack &Mem, double fFactor, uint Flags) {
   FMatrixView *M[5] = {&A0, &A1, &A2, &A3, &A4};
   return ChainMxmR(Out, &M[0], 5, Mem, fFactor, Flags);
}

void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMatrixView A5, FMemoryStack &Mem, double fFactor, uint Flags) {
   FMatrixView *M[6] = {&A0, &A1, &A2, &A3, &A4, &A5};
   return ChainMxmR(Out, &M[0], 6, Mem, fFactor, Flags);
}


FHeapMatrix::FHeapMatrix()
   : FMatrixView(0, 0, 0)
{
}

FHeapMatrix::FHeapMatrix(FMatrixView v)
{
   Reshape(v.nRows, v.nCols);
   Assign(*this, v);
}

FHeapMatrix::FHeapMatrix(size_t nRows, size_t nCols)
{
   Reshape(nRows, nCols);
}

FHeapMatrix::~FHeapMatrix()
{
}

void FHeapMatrix::Reshape(size_t nRows_, size_t nCols_)
{
   this->nRows = nRows_;
   this->nCols = nCols_;
   this->nRowSt = 1;
   this->nColSt = nRows_;
   m_Buffer.resize(nRows_ * nCols_);
   this->pData = &m_Buffer[0];
}

FHeapMatrix::FHeapMatrix(FHeapMatrix const &other)
   : FMatrixView(0, 0, 0), FIntrusivePtrDest(other) // <- note: FIntrusivePtrDest c'tor just sets ref count to 0. It does *NOT* copy the ref count (that would be wrong)
{
   Assign(*this, other);
}

void FHeapMatrix::operator=(FHeapMatrix const &other)
{
   Assign(*this, other);
}


std::ostream &operator << (std::ostream &out, FMatrixShape const &Shape) {
   out << "(" << Shape.nRows << ", " << Shape.nCols << ")";
   return out;
}


// prints a general rectangular matrix specified by memory layout.
// Element m(r,c) is given by pData[r * nRowStride + c * nColStride].
template<class FFloat>
void PrintMatrixGenT(std::ostream &out, FFloat const *pData, size_t nRows, size_t nRowStride, size_t nCols,
    size_t nColStride, std::string const &Name, uint nFloatWidth, uint nFloatPrec)
{
   std::string
      ff = fmt::format("{{:{}.{}f}}", nFloatWidth, nFloatPrec),
      fi = fmt::format("{{:>{}}}", nFloatWidth - 3);
   if (Name != "") {
      out << fmt::format("  Matrix {}, {}x{}.", Name, nRows, nCols) << std::endl;
      if (nRows * nCols == 0)
         return;
   }
   out << "           ";
   for (size_t iCol = 0; iCol < nCols; ++iCol)
      out << fmt::format(fi, iCol);
   out << "\n";
   for (size_t iRow = 0; iRow < nRows; ++iRow) {
      out << fmt::format("    {:>4}   ", iRow);
      for (size_t iCol = 0; iCol < nCols; ++iCol) {
         FFloat const &f = pData[iRow * nRowStride + iCol * nColStride];
         out << fmt::format(ff, f);
      }
      out << "\n";
   }
   out.flush();
}


void PrintMatrixGen(std::ostream &out, double const *pData, size_t nRows, size_t nRowStride, size_t nCols, size_t nColStride, std::string const &Name)
{
   uint
      nFloatWidth = 14, // 25, //20, //14,
      nFloatPrec = 8;    // 16; //15; //8;
   return PrintMatrixGenT(out, pData, nRows, nRowStride, nCols, nColStride, Name, nFloatWidth, nFloatPrec);
}

void PrintMatrixGen(std::ostream &out, float const *pData, size_t nRows, size_t nRowStride, size_t nCols, size_t nColStride, std::string const &Name)
{
   uint
      nFloatWidth = 11, // 25, //20, //14,
      nFloatPrec = 5;    // 16; //15; //8;
   return PrintMatrixGenT(out, pData, nRows, nRowStride, nCols, nColStride, Name, nFloatWidth, nFloatPrec);
}


void ReadMatrixFromFile(FMatrixView &M, std::string const &FileName, FMemoryStack &Mem)
{
   // format:
   //    <matrix name>
   //    <nRows> x <nCols> <SYMMETRIC/RECTANGULAR>
   //    <irow> <icol> <m[irow,icol]>
   std::string
      MatrixName, Dummy, Type;
   std::ifstream
      File(FileName.c_str());
   if (!File.good())
      throw std::runtime_error("failed to open: " + FileName + ".");
   File >> MatrixName;
   size_t
      nRows, nCols;
   std::getline(File, Dummy); // matrix name
   File >> nRows >> Dummy >> nCols >> Type;
   if (Dummy != "x" || (Type != "RECTANGULAR" && Type != "SYMMETRIC"))
      throw std::runtime_error("unexpected format: input matrix " + FileName + ".");
   bool
      Symmetric = false;
   Symmetric = Type == "SYMMETRIC";
   if (Symmetric && nRows != nCols)
      throw std::runtime_error("symmetric matrices must have equal number of columns and rows: input matrix " + FileName + ".");

   M = FMatrixView(0, nRows, nCols);
   Mem.ClearAlloc(M.pData, M.GetStridedSize());
   while (true) {
      int
         iRow, iCol;
      double
         f;
      File >> iRow >> iCol >> f;
      if (!File.good())
         break;
      if (iRow < 1 || (unsigned)iRow > nRows || iCol < 1 || (unsigned)iCol > nCols) {
         std::stringstream str;
         str << fmt::format("index out of bounds: iRow = {}  iCol = {} while reading {} x {} matrix {}.",
            iRow, iCol, nRows, nCols, FileName);
         throw std::runtime_error(str.str());
      }
      M(iRow-1, iCol-1) = f;
      if (Symmetric)
         M(iCol-1, iRow-1) = f;
   }
}

void WriteMatrixToFile(std::string const &FileName, FMatrixView const &M, std::string const &MatrixName, size_t *pRowPerm, size_t *pColPerm)
{
   // format:
   //    <matrix name>
   //    <nRows> x <nCols> <SYMMETRIC/RECTANGULAR>
   //    <irow> <icol> <m[irow,icol]>
   std::ofstream
      File(FileName.c_str());
   if (!File.good())
      throw std::runtime_error("failed to open: " + FileName + ".");
   bool
      Symmetric = M.IsSymmetric(1e-10);
   File << fmt::format("{}\n{} x {} {}\n", MatrixName, M.nRows, M.nCols, (Symmetric? "SYMMETRIC" : "RECTANGULAR"));

   for (size_t iCol = 0; iCol < M.nCols; ++ iCol)
      for (size_t iRow = 0; iRow < M.nRows; ++ iRow) {
         if (Symmetric && iCol < iRow)
            continue;
         size_t
            iRow_ = pRowPerm? pRowPerm[iRow] : iRow,
            iCol_ = pColPerm? pColPerm[iCol] : iCol;
         FScalar
            f = M(iRow_,iCol_);
         if (Symmetric)
            f = .5 * (f + M(iCol_,iRow_));
         if (std::abs(f) > 1e-25)
            File << fmt::format("{:>5} {:>5} {:24.16e}\n", (iRow+1), (iCol+1), f);
      }
}


} // namespace ct

// kate: space-indent on; tab-indent off; backspace-indent off; indent-width 3; mixedindent off; indent-mode normal;
