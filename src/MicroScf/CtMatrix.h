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

#ifndef CT8k_MATRIX_H
#define CT8k_MATRIX_H

#include "CxTypes.h" // for size_t.
#include "CxAlgebra.h"
#include "CxMemoryStack.h"
#include "CxPodArray.h" // for TArray<>; FHeapMatrix data is stored in those.

#include <string>
#include <sstream>

namespace ct {

typedef double
    FScalar;

using std::size_t;

// size_t const
//    SIZE_Invalid = size_t(-1);
struct FMatrixShape {
   size_t nRows, nCols;
//    // says whether matrix from which this comes has a memory location assigned;
//    // if not, one can argue that it does not have any valid shape at all.
//    // Not sure if this a great way of dealing with this...
//    bool IsAssigned;
   FMatrixShape(size_t nRows_, size_t nCols_) : nRows(nRows_), nCols(nCols_) {}
   bool operator == (FMatrixShape const &other) const { return nRows == other.nRows && nCols == other.nCols; }
   bool operator != (FMatrixShape const &other) const { return nRows != other.nRows || nCols != other.nCols; }
   bool IsSquare() const { return nRows == nCols; }
};

std::ostream &operator << (std::ostream &out, FMatrixShape const &Shape);

/// Describes the layout of a matrix somewhere in memory, which this object
/// DOES NOT own.
struct FMatrixView
{
   typedef FScalar
       FScalarType;
   typedef FScalar
       value_type;
   FScalar
       *pData; ///< start of data of this matrix.
   size_t
       nRows, ///< number of rows (first index, fast)
       nCols, ///< number of columns (second index, slow)
       nRowSt, ///< number of words between adjacent rows. In default-alignment, nRowSt == 1
       nColSt;  ///< number of words between adjacent cols. In default-alignment, nColSt == nRows
   // Note: One of both strides should be 1.

   /// element retrieval operator: m(nRow,nCol) -> element stored at that
   /// part of the matrix. Note: indices are 0-based!.
   inline FScalar& operator () (size_t iRow, size_t iCol){
       assert(iRow <= nRows);
       assert(iCol <= nCols);
       return pData[nRowSt * iRow + nColSt * iCol];
   };
   inline FScalar const& operator () (size_t iRow, size_t iCol) const {
       return const_cast<FMatrixView*>(this)->operator()(iRow, iCol);
   };

   FScalar &operator [] (size_t iEntry) { return pData[iEntry]; }
   FScalar const &operator []  (size_t iEntry) const { return pData[iEntry]; }

   FMatrixView()
       : pData(0)
   {}

   FMatrixView(double *pData_, size_t nRows_, size_t nCols_, size_t nRowSt_ = 1, size_t nColSt_ = 0)
       : pData(pData_), nRows(nRows_), nCols(nCols_),
         nRowSt(nRowSt_), nColSt((nColSt_ != 0)? nColSt_ : nRows)
   {};

   inline size_t GetStridedSize() const;

   void Print(std::ostream &out, std::string const &Name = "") const;
   /// sets matrix to zero
   void Clear();
   /// sets matrix to identity
   void SetIdentity();
//     /// returns copy of diagonal elements of the matrix. Only valid for square matrices.
//     std::vector<FScalar> GetDiagonal() const;

   bool IsSquare() const { return nRows == nCols; };
   bool IsSymmetric(FScalar Thresh=1e-10) const;
   /// returns whether the data referenced by this view is a single continuous block of memory,
   /// without any gaps between data in rows or columns.
   bool IsDataContinuous() const;
   FMatrixShape Shape() const { return FMatrixShape(nRows, nCols); }
   bool IsAssignedAndHasShape(size_t nRows1, size_t nCols1) const { return pData != 0 && nRows == nRows1 && nCols == nCols1; }
   bool IsViewIdentical(FMatrixView const &other) const { return pData == other.pData && nRows == other.nRows && nCols == other.nCols && nRowSt == other.nRowSt && nColSt == other.nColSt; }
   bool HasVecs() const { return pData != 0 && nCols != 0; }
//    bool NoVecs() const { return pData == 0 || nCols == 0; }


   double fRmsdFromIdentity() const;
   double fRmsdFromZero(size_t nParamsForDivisor=0) const;
   double fRmsdFromSymmetry() const;

   /// if nRowSt==1, sets space between nRow and nColSt to zero (in order to
   /// not have possible NaNs/Denorms in spaces in the matrices which are only
   /// there for alignment purposes; e.g., ntb vs nt strides.
   void ClearEmptySpace();

   /// returns size of this matrix would require in triangular storage (e.g.,
   /// storing only the lower triangular part). If Sign == -1, size is
   /// calculated for anti-symmetric matrix.
   size_t TriangularStorageSize(int Sign = 1) const;
   size_t TriangularStorageSize1() const { return TriangularStorageSize(1); };
   void TriangularExpand(int Sign = 1, FScalar const *pTriangData = 0);
   // fOffDiagFactor: if given, off diagonal elements will be multiplied by this factor before reduction.
   void TriangularReduce(int Sign = 1, FScalar *pTriangData = 0, FScalar fOffDiagFactor = 1., bool TransposeSource = false);
   enum FTriangularTransformFlags {
       TRIANG_Expand,
       TRIANG_Reduce
   };
   void TriangularExpandOrReduce(int Sign, FScalar *pTriangData, FTriangularTransformFlags Direction);
};


size_t FMatrixView::GetStridedSize() const
{
   if (nRowSt == 1) {
      return nCols * nColSt;
   } else {
      if (nColSt == 1) {
         return nRows * nRowSt;
      } else {
         assert(0); // not supported.
         return 0;
      }
   }
}


/// asserts that A and B have the same number of rows and columns. (used to be called AssertSameShape1)
void AssertSameShape1(FMatrixView const &A, FMatrixView const &B);

/// allocate a nRows x nCols matrix on Mem and return the given object.
inline FMatrixView MakeStackMatrix(size_t nRows, size_t nCols, FMemoryStack &Mem)
{
   FMatrixView r(0, nRows, nCols);
   Mem.Alloc(r.pData, nRows * nCols);
   return r;
}

uint const MXM_AddToDest = 0x1u;
uint const MXM_Add = MXM_AddToDest;

/// Creates a view representing the transpose of the input matrix.
/// The returned view still references the original data (i.e., In and Out are
/// aliased)
FMatrixView Transpose(FMatrixView const &In);
/// Out += f * In
void Add(FMatrixView Out, FMatrixView const &In, FScalar fFactor = 1.0);
/// Out = f * In   [Note: used to be called "Move", but always copied data, rather than moving it; maybe call it 'Assign' or 'Set?'?]
void Assign(FMatrixView Out, FMatrixView const &In, FScalar fFactor = 1.0);
/// Out = A * B
void Mxm(FMatrixView &Out, FMatrixView const &A, FMatrixView const &B, uint Flags = 0);
/// Out = f * A * B
void Mxm(FMatrixView &Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags = 0);
/// Out = f * A * A^T for symmetric matrix Out [Note: Out will be symmetrized!].
void SyrkNT(FMatrixView &Out, FMatrixView const &A, FScalar fFactor = 1.0, uint Flags = 0);
/// Out = f * A^T * A for symmetric matrix Out [Note: Out will be symmetrized!].
void SyrkTN(FMatrixView &Out, FMatrixView const &A, FScalar fFactor = 1.0, uint Flags = 0);


/// Out += Factor * TriangularExpand(Sign,pDataTr,fOffDiagFactor)
/// Note: Out and pDataTr may not overlap.
void Add_TriangularExpanded(FMatrixView Out, int Sign, double const *pDataTr, double fFactor, double fOffDiagFactor=1.);


/// Chained matrix product: Out = A0 * A1 * A2....
void ChainMxm(FMatrixView Out, FMatrixView A0, FMemoryStack &Mem, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMemoryStack &Mem, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMemoryStack &Mem, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMemoryStack &Mem, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMemoryStack &Mem, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMatrixView A5, FMemoryStack &Mem, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMemoryStack &Mem, double fFactor, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMemoryStack &Mem, double fFactor, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMemoryStack &Mem, double fFactor, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMemoryStack &Mem, double fFactor, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMemoryStack &Mem, double fFactor, uint Flags = 0);
void ChainMxm(FMatrixView Out, FMatrixView A0, FMatrixView A1, FMatrixView A2, FMatrixView A3, FMatrixView A4, FMatrixView A5, FMemoryStack &Mem, double fFactor, uint Flags = 0);


/// Out = A * B, Out,B being vectors
void Mxva(FScalar *RESTRICT pOut, FMatrixView const &A, FScalar const *RESTRICT pIn);
/// Out += A * B, Out,B being vectors
void Mxvb(FScalar *RESTRICT pOut, FMatrixView const &A, FScalar const *RESTRICT pIn);
/// trace of a matrix
FScalar Trace(FMatrixView const &In);
/// dot-product of two matrices. dot(A^T B) = Tr(A * B).
FScalar Dot(FMatrixView const &A, FMatrixView const &B);
/// root mean square deviation of the matrix elements
FScalar Rmsd(FMatrixView const &A, FMatrixView const &B);
/// scale matrix by a factor in place
void Scale(FMatrixView &Out, FScalar fFactor);
/// diagonalize a real symmetric matrix in place and store the eigenvalues at
/// pEigenValues. Input matrix must be real symmetric and pEigenValues must hold
/// room for nRows==nCols entries.
void Diagonalize(FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem);
/// As Diagonalize(), but eigenvalues are returned in descending order.
void Diagonalize_LargeEwFirst(FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem);

/// Factorize `In` into `U * diag(sigma) * Vt`, where `U` and `V = transpose(Vt)` are
/// have orthonormal column vectors (i.e., `U.T * U = id` and `V.T * V = id`)
/// Inputs:
/// - In: (nRows, nCols) matrix.
///   Content of In is left alone (input is copied).
/// Outputs:
/// - pSigma: pointer to nSig-length array of singular values, where nSig = min(nRows,nCols)
/// - U: (nRows, nSig)-shape matrix, where nSig = min(nRows,nCols)
/// - Vt: (nSig, nCols)-shape matrix, where nSig = min(nRows,nCols)
///   (so V is a (nCols, nSig)-shape matrix, but this routine returns Vt = transpose(V),
///   not V itself! That's for LAPACK compatibility, which follows the same convention)
/// Notes:
/// - This routine *does not* allocate the output matrices!
///   Output matrices U and Vt are assumed to have sufficient size to support the SVD.
///   If unsure, use AllocAndComputeSvd
void ComputeSvd(FMatrixView &U, FScalar *pSigma, FMatrixView &Vt, FMatrixView In, FMemoryStack &Mem);

/// As ComputeSvd, but allocate output quantities (U, pSigma, and Vt) on Mem.
void AllocAndComputeSvd(FMatrixView &U, FScalar *&pSigma, FMatrixView &Vt, FMatrixView In, FMemoryStack &Mem);
void AllocAndComputeSvd(FMatrixView &U, FMatrixView &Sigma, FMatrixView &Vt, FMatrixView In, FMemoryStack &Mem);

/// calculate 'M's lower holesky factorization: M = L * L^T; this overwrites M with L.
void CalcCholeskyFactors(FMatrixView M);
/// invert lower triangular matrix L. Upper part ignored.
void InvertTriangularMatrix(FMatrixView L);
/// Solve L X = B where L is a lower triangular matrix.
void TriangularSolve(FMatrixView RhsSol, FMatrixView const &L, char Side = 'L', char TransL = 'N', double fFactor = 1.0);
/// set A := op(L) * A where L is a lower triangular matrix (Side == 'L').
/// or  A := A * op(L) where L is a lower triangular matrix (Side == 'R')
/// op(L) is identity if TransL == 'N' and transpose if TransL=='T'.
void TriangularMxm(FMatrixView A, FMatrixView const &L, char Side, char TransL = 'N', double fFactor = 1.0);
/// solve linear least squares ||M * x - b||_2 -> min subject to ||x||_2 ->min
/// using singular value decomposition of A. Singular values below fThrRel are
/// treated as zero. A negative value indicates the usage of machine precision.
/// Both M and RhsSol are overwritten.
///
/// Notes:
///    - On input, RhsSol is the RHS and must have format M.nRows x nRhs
///      On output, RhsSol is the solution and is set to format M.nCols x nRhs
///    - This happens IN PLACE. That means that RhsSol must have a column stride
///      of at least std::max(M.nCols, M.nRows).
///    - Number of rows on RhsSol is adjusted as required on output.
void LeastSquaresSolveUnsafe(FMatrixView M, FMatrixView &RhsSol, double fThrRel, FMemoryStack &Mem);
/// Solve M * Sol == Rhs.
/// Similar to LeastSquaresSolve, but Rhs and Sol are provided independently (and are copied
/// if required for reasons of strides). Neither M nor Rhs is overwritten.
void LeastSquaresSolveSafe(FMatrixView const M, FMatrixView Sol, FMatrixView const Rhs, double fThrRel, FMemoryStack &Mem);


// // Solve X L = B where L is a lower triangular matrix.
// void TriangularSolveRight(FMatrixView RhsSol, FMatrixView const &L);

/// Solve A X = B, where A = L*L^T is implicitly given by the input lower triangular matrix L.
void CholeskySolve(FMatrixView RhsAndSolution, FMatrixView const &L);
/// Solve A X = B, where A = L*L^T is implicitly given by the input lower triangular matrix L.
void CholeskySolve(double *pRhsAndSolution, FMatrixView const &L);
/// Multiply X = A * B, where A = L*L^T is implicitly given by the input lower triangular matrix L.
void CholeskyMxm(FMatrixView RhsSol, FMatrixView const &L);
/// Multiply B = A * X, where A = L*L^T is implicitly given by the input lower triangular matrix L.
void CholeskyMxm(double *pRhsAndSolution, FMatrixView const &L);

void Symmetrize(FMatrixView M, double Phase = 1);

/// creates a view representing the sub-matrix of size nRow x nCol
/// starting at (iRow,iCol).
/// The returned view still references the original data (i.e., In and Out are aliased)
FMatrixView Select(FMatrixView const &In, size_t iRow, size_t iCol, size_t nRows, size_t nCols);

// exactly same as above, but allowing for temporaries as first argument (C++ ftw...).
inline void Add0(FMatrixView Out, FMatrixView const &In, FScalar fFactor = 1.0) {
   return Add(Out, In, fFactor);
}
inline void Assign0(FMatrixView Out, FMatrixView const &In, FScalar fFactor = 1.0) {
   return Assign(Out, In, fFactor);
}
inline void Mxm0(FMatrixView Out, FMatrixView const &A, FMatrixView const &B, uint Flags = 0) {
   return Mxm(Out, A, B, Flags);
}
inline void Mxm0(FMatrixView Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags = 0) {
   return Mxm(Out, A, B, fFactor, Flags);
}
inline void Scale0(FMatrixView Out, FScalar fFactor) {
   return Scale(Out, fFactor);
}



/// A matrix object with FMatrixView interface, but which actually owns its
/// data; the data is placed on a FMemoryStack object upon construction and
/// release upon destruction.
///
/// Notes:
/// - Max data size must be specified on construction; will adjust to input dimensions
///   of operations as long as the actual matrix size does not exceed nMaxSize.
///
///   Typically the actual number of rows and columns is given on construction and
///   not changed afterwards, e.g.:
///
///       FStackMatrix M(123, 456, &Mem);
///
/// WARNING:
/// - This will free the memory stack up to its own data allocation upon
///   destruction. This implies that all data allocated on the stack after the
///   matrix will be freed, too!
class FStackMatrix : public FMatrixView
{
public:
   FStackMatrix(size_t nMaxSize_, FMemoryStack *pMemoryStack_)
      : nMaxSize(nMaxSize_), pMemory(pMemoryStack_)
   {
      pMemory->Align(64);
      pMemory->Alloc(pData, nMaxSize_);
      Reshape(0, 0);
   }

   FStackMatrix(size_t nRows_, size_t nCols_, FMemoryStack *pMemoryStack_)
      : nMaxSize(nRows_ * nCols_), pMemory(pMemoryStack_)
   {
      pMemory->Align(64);
      pMemory->Alloc(pData, nMaxSize);
      Reshape(nRows_, nCols_);
   }

   ~FStackMatrix() { Free(); }
   void Free()
   {
      if (pData)
         pMemory->Free(pData);
      pData = 0;
   }

   void Reshape(size_t nRows_, size_t nCols_)
   {
      assert(nRows_ * nCols_ <= nMaxSize);
      nRows = nRows_;
      nCols = nCols_;
      nRowSt = 1;
      nColSt = nRows_;
   }

private:
   size_t
      nMaxSize;
   FMemoryStack
      *pMemory;
};


/// A matrix object with FMatrixView interface, but which actually owns its
/// data; the data is placed on the main heap.
///
/// WARNING:
/// - This one has *NO* virtual destructor! and FMatrixView is not going to get one.
/// FIXME:
/// ../lhf/CtMatrix.cpp: In copy constructor ‘ct::FHeapMatrix::FHeapMatrix(const ct::FHeapMatrix&)’:
/// ../lhf/CtMatrix.cpp:1734:1: warning: base class ‘struct ct::FIntrusivePtrDest1’ should be explicitly initialized in the copy constructor [-Wextra]
/// 1734 | FHeapMatrix::FHeapMatrix(FHeapMatrix const &other)
class FHeapMatrix : public FMatrixView, public FIntrusivePtrDest
{
public:
   FHeapMatrix();
   FHeapMatrix(FMatrixView v);
   FHeapMatrix(size_t nRows, size_t nCols);
   virtual ~FHeapMatrix();
   void Reshape(size_t nRows_, size_t nCols_);

   FHeapMatrix(FHeapMatrix const &other);
   void operator=(FHeapMatrix const &other);

   void Delete() { Reshape(0,0); pData = 0; }
protected:
   TArray<FScalar>
      m_Buffer;
};
typedef TIntrusivePtr<FHeapMatrix>
   FHeapMatrixPtr;


/// Reshape Out to A*B's dimensions, and assign Out = f * In
void Assign(FStackMatrix &Out, FMatrixView const &In, FScalar fFactor = 1.0);
/// Reshape Out to A*B's dimensions, and assign Out = A * B
void Mxm(FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, uint Flags = 0);
/// Reshape Out to A*B's dimensions, and assign Out = f * A * B
void Mxm(FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags = 0);
/// Reshape Out to In's dimensions, and assign Out := f * In
void Assign(FHeapMatrix &Out, FMatrixView const &In, FScalar fFactor = 1.0);
// ^- note: the other operations from FMatrixView should work unchanged
// on FStackMatrix, because they do not require a change in the shape.

/// Print a nRows x nCols matrix to stream out, where M[iRow,iCol] = pData[nRowStride * iRow + nColStride * iCol]
void PrintMatrixGen(std::ostream &out, double const *pData, size_t nRows, size_t nRowStride, size_t nCols,
    size_t nColStride, std::string const &Name);
/// Print a nRows x nCols matrix to stream out, where M[iRow,iCol] = pData[nRowStride * iRow + nColStride * iCol]
void PrintMatrixGen(std::ostream &out, float const *pData, size_t nRows, size_t nRowStride, size_t nCols,
    size_t nColStride, std::string const &Name);


struct FSmhOptions {
    double
        // ignore eigenvectors with eigenvalue < std::max(ThrAbs, MaxEw * ThrRel)
        ThrAbs, ThrRel;
    char const
        *pDelMsg; // may be 0.
    std::ostream
        *pXout; // may be 0 if pDelMsg is 0.
    FSmhOptions(double ThrAbs_ = 0, double ThrRel_ = 0, char const *pDelMsg_ = 0,
            std::ostream *pXout_ = 0)
        : ThrAbs(ThrAbs_), ThrRel(ThrRel_), pDelMsg(pDelMsg_), pXout(pXout_)
    {}
};

/// Replace M by M^{-1/2} (matrix inverse square root).
void CalcSmhMatrix(FMatrixView M, FMemoryStack &Mem, FSmhOptions const &Opt);
/// Replace M by M^p (matrix raised to power `fExp`). M must be symmetric and
/// positive semi-definite (if `fExp` is negative, then strictly positive definite).
/// Opt.Thr*are only used if `fExp` < 0.
void CalcSymMatrixPower(FMatrixView M, double fExp, FMemoryStack &Mem, FSmhOptions const &Opt);
// // orthogonalizes C in place, where S is C's overlap matrix.
// // both C and S are overwritten.
// void OrthSchmidt1(FMatrixView S, FMatrixView C, FMemoryStack &Mem);

/// Schmidt-orthogonalize the nAo x nOcc orbitals in C, with respect to the nAo x nAo overlap matrix S1.
void OrthSchmidt(FMatrixView Orbs, FMatrixView S1, FMemoryStack &Mem);
/// Symmetrically orthogonalize the nAo x nOcc orbitals in C, with respect to the nAo x nAo overlap matrix S1.
void SymOrth(FMatrixView Orbs, FMatrixView const S1, FMemoryStack &Mem, double fThrAbs = 1e-15, double fThrRel = 1e-15);
/// schmidt-orthogonalize the nAo x nOcc orbitals in C with respect to the nAo x nAo overlap matrix S1,
/// which is implicitly given by its lower Choleksy decomposition: S1 = Scd1 * Scd1.T
void OrthSchmidtScd(FMatrixView C, FMatrixView Scd1, FMemoryStack &Mem);


/// represents a stored permutation of N indices: An array of values from the
/// set {0, 1, ..., N-1} in which every value occurs exactly once.
///
/// Convention is such that:
///
///   New[iFn] = Old[P[iFn]]
typedef TArray<size_t>
   FBasicPerm;
FBasicPerm Inverted(FBasicPerm const &in);
/// it's just a wrapper for ArgSort1. Finds and returns the set permutation of
/// indices I such that [Val[i] for i in I] is a sorted array. Val[i] is
/// accessed as pVals[i*iValSt]. The sort in ArgSort is stable.
FBasicPerm ArgSort(double const *pVals, size_t iValSt, size_t nVals, bool Reverse = false);
/// it's just a wrapper for ArgSort1. On output, `Indices` is set to the index
/// permutation for which [Val[i] for i in Indices] is a sorted array. Val[i] is
/// accessed as pVals[i*iValSt]. The sort in ArgSort is stable.
void ArgSort(FBasicPerm &Indices, double const *pVals, size_t iValSt, size_t nVals, bool Reverse = false);



/// Permute rows of matrix such that
///    NewM[iRow,iCol] = OldM[P[iRow], iCol]
void PermuteRows(FMatrixView M, size_t const *P, FMemoryStack &Mem);
/// Permute both rows and columns in convention given by PermuteRows.
void Permute(FMatrixView &M, size_t const *pRowPerm, size_t const *pColPerm, FMemoryStack &Mem);

/// Write a human-readably formatted matrix to a file
void WriteMatrixToFile(std::string const &FileName, FMatrixView const &M, std::string const &MatrixName,
    size_t *pRowPerm = 0, size_t *pColPerm = 0);

/// Read back a matrix written by WriteMatrixToFile; M is allocated on 'Mem'.
void ReadMatrixFromFile(FMatrixView &M, std::string const &FileName, FMemoryStack &Mem);

// void ArgSort1(size_t *pOrd, double const *pVals, size_t nValSt, size_t nVals, bool Reverse = false);




} // namespace ct

#endif // CT8k_MATRIX_H

// kate: space-indent on; tab-indent off; backspace-indent off; indent-width 3; mixedindent off; indent-mode normal;
