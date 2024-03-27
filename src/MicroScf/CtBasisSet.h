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

#ifndef CT8K_BASISSET_H
#define CT8K_BASISSET_H

#include <vector>
#include <list>
#include "Ir.h"
#ifdef IR_ECP
   #include "IrEcp.h"
#endif // IR_ECP
#include "IrMeta2i.h"
#include "CxPodArray.h"
#include "CtBasisShell.h"
#include "CxRawAtom.h"
#include "CtBasisDesc.h"

namespace ct {


struct FMatrixView;
struct FAtomSet;


// describes a single symmetry adapted basis function (SO) by giving its
// linear combination in terms of AOs.
struct FSymOrb {
   static const unsigned
      iNegWeight = 0x80000000,
      iAoMask = iNegWeight - 1;
   unsigned
      // nEqiv[iSo]:  nti = number of equivalent AOs forming SO iSo.
      nEqivAo,
      nEqivAoCase, // 0, 1, 2 or 3. nEqiv = 2^nEqivAoCase
      // iEqiv[0..nti-1, iSo]:  indices of the equivalent AOs. If >= iNegWeight,
      //   actual index is iEqiv & iAoMask but weight in stabilizer is negative.
      iEqivAo[8];
   double
      // Weight[iSo]:       1/sqrt(nti) -> weight in stabilizer (equal for all AOs).
      Weight;
};


struct FUniqueShell { // a non-redundant basis function shell
   unsigned iSh;
   double fWeight; // number of times it occurs symmetry-equivalently.
};


// describes the transformation between a raw basis set and a symmetry
// adapted basis set
struct FSymTrans {
   unsigned
      nIrreps,
      // offset of first SO transforming according to that irrep.
      // Length: nIrrep+1. Last entry contains total number.
      IrrepOffsets[9],
      // number of symmetry adapted basis functions (== number of basis functions)
      nSo,
      nUniqueShells;
   FSymOrb
      // for each SO, gives the AO linear combination forming it.
      *pSo;
   FUniqueShell
      *pUniqueShells;
};



// program native format of a basis set.
struct FBasisSet : FIntrusivePtrDest
{
   typedef std::vector<FBasisShell>
      FBasisShellArray;
   FBasisShellArray
      Shells;  // Gaussian shells centered on different atoms.
   FBasisContext
      Context; // for informative purposes only.
   std::string
      Name; // for informative purposes only.

   enum FPrintFlags {
      PRINT_MetaInfo = 0x0001,
      PRINT_FullDetail = 0x0002,
      PRINT_Default = PRINT_MetaInfo
   };

   void Print(std::ostream &out, unsigned Flags = PRINT_Default, std::string const &Ind = "  ") const;
   std::string GetDesc(unsigned Flags = PRINT_Default, std::string const &Ind = "  ") const;
   size_t nFn() const;
   size_t nPrimitiveGtos() const; // for decoration purposes only.
   size_t nFnOfLargestShell() const;
   size_t nCoOfLargestShell() const;
   unsigned nMaxL() const; // find largest angular momentum in the basis.
   size_t nCenters() const; // highest center ID of a basis function occuring in the set
   // ^- FIXME: some parts of the code assume that the number of centers is equal to the
   // number of atoms. This is not necessarily the case, as
   // (a) there might be dummy atoms or other centers carrying basis functions, but no
   //     atoms (e.g., bond-centered functions) and
   // (b) there might be atoms without basis functions (e.g., implicit QM Hs)
   size_t nCen() const { return nCenters(); };

   // rotate and/or translate the basis (i.e., the positions of the atoms) in real space. R is a
   // 3x3 matrix and 'd' is a 3-vector such that r' = R (r - d) gives the transformed atomic positions.
   void Transform_3x4(double const *R, double const *d, FMemoryStack &Mem);

   // transform a set of MO coefficients (in C) in such a way that if the molecule and its
   // basis are rotated in real space by the 3x3 matrix R (via FBasisSet Transform_3x4), then
   // this function computes the MO coefficients corresponding to the new alignment.
   void TransformMoCoeffs_3x4(FMatrixView C, double const *R, FMemoryStack &Mem) const;

   // this identifies the unique element bases in the atoms of this basis set,
   // writes a corresponding basis specification in Molpro basis library format
   // to output stream xout. Atoms is only needed to identify the element for
   // each center.
   //
   // Note: This assumes that basis sets for different atoms of the same element
   // are consistent with each other. This is not checked!
   void ExportToLibmol(std::ostream &xout, std::string const &BasisName, std::string const &Comment, FAtomSet const &Atoms) const;

   // Checks if the same FAtomShell objects are linked for all shells in `*this`
   // and `Other`. It does not compare the positions of the basis functions.
   bool HasSameElementBases(FBasisSet const &Other) const;

   // construct from atom set and context
   FBasisSet(FAtomSet const &AtomSet, FBasisContext Context_);
   // construct from explicit list of shells. Name and Context for informative purposes only.
   FBasisSet(FBasisShell const *pShells_, size_t nShells_, FBasisContext Context_, std::string const &Name);
public:
   FRawBasisPtr
      // Only valid if RebuildInterfacingData has been called after the basis
      // was last modified. If constructed regularly, this is done
      // automatically.
      pRawBasis;

   // rebuilds the data for low-level interfaces (i.e., the data in this section)
   // from this->Shells
   void RebuildInterfacingData();

   FRawBasisPtr MakeRawBasis();
public:
   void MakeSymmetryTransformation(FSymTrans &st, FAtomSet const &AtomSet, FMemoryStack &Mem);
   void PrintSoBasis(std::ostream &xout, FAtomSet const &AtomSet, FMemoryStack &Mem);
private:
   // makes *this a basis set corresponding to the basis descriptions in
   // `Atoms`. Attempts to load the neccessary external data. Throws
   // std::runtime_error if failed.
   void LoadFromAtomSet(FAtomSet const &Atoms, FBasisContext Context);

   void MakeAtomOffsets( size_t *&pAtomShellOffsets, size_t *&pAtomBfnOffsets, size_t nAtoms_, FMemoryStack &Mem ) const;
   void MakeShellOffsets( size_t *&pShellOffsets, FMemoryStack &Mem ) const;
   void Finalize();
};

typedef TIntrusivePtr<FBasisSet>
   FBasisSetPtr;

// construct a basis set from the library which is *exactly* named as "BasisName" for all elements.
// There is no way to override the basis to load for specific atoms or elements.
// Use with AsRawAtoms() in case you need it for a FXyzFrame or a FAtomSet.
//
// WARNING: This does not apply any of the "you said basis xxx, but you
// probably meant..." logic in CtBasisDesc.cpp. Use at your own risk.
FBasisSetPtr LoadNamedBasis(FRawAtomList const &Atoms, FBasisContext Context, std::string const &BasisName);


void Trafo3x4(FVector3 &InOut, double const *R, double const *d);


inline std::ostream &operator << (std::ostream &out, FBasisSet const &BasisSet) {
   BasisSet.Print(out);
   return out;
}


void MakeIntMatrix( FMatrixView &Out, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor=1., bool Add=false);
void MakeIntMatrix( FMatrixView &Out, FRawBasis const &RowBasis,
   FRawBasis const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor=1., bool Add=false);
// increment gradient in pGrad[ixyz + 3 *iCenter] by \sum{a,b} Rdm[a,b] d/d[Eu] Krn2i[a,b].
// If &RowBasis == &ColBasis, Rdm is assumed to be symmetric. If Rdm is not symmetric, pls
// symmetrize first (not checked!).
void AccGradient2ix( double *pGrad, FMatrixView const &Rdm, FRawBasis const &RowBasis,
   FRawBasis const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor=1.);

// increment Hessian in Hess(ixyzA + 3 *iCenterA, ixyzB + 3 *iCenterB) by
// \sum{a,b} Rdm[a,b] d^2/)d[Au]d[Bv]) Krn2i[a,b].
// If &RowBasis == &ColBasis, Rdm is assumed to be symmetric. If Rdm is not symmetric, pls
// symmetrize first (not checked!).
void AccHessian2ix( FMatrixView Hess, FMatrixView const &Rdm, FRawBasis const &RowBasis,
   FRawBasis const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor=1.);



enum FSymTransDir {
   ST_AoToSo = 0,
   ST_SoToAo = 1
};

/// transform, in place, the row dimension of a matrix between raw basis
/// functions (AOs) and symmetry adapted basis functions (SOs)
/// \p pData   Matrix; Input & Output
/// \p Dir     Determines if AO->SO or SO->AO transform is done
/// \p nRowSt  Stride between two consecutive rows of the matrix
/// \p nColSt  Stride between two consecutive columns of the matrix
/// \p nCols   Number of columns
/// \p st      Determines the symmetry space for the rows (by this also nRows!).
/// \p Mem     Space for temporaries
void SymTrans1(double *pData, uint nRowSt, uint nColSt, uint nCols, FSymTrans const &st, FSymTransDir Dir, FMemoryStack &Mem);
void SymTrans2(double *pData, uint nRowSt, uint nColSt, FSymTrans const &RowSy, FSymTrans const &ColSy, FSymTransDir Dir, FMemoryStack &Mem);
void SymTrans2(FMatrixView &M, FSymTrans const &RowSy, FSymTrans const &ColSy, FSymTransDir Dir, FMemoryStack &Mem);



} // namespace ct

#endif // CT8K_BASISSET_H
