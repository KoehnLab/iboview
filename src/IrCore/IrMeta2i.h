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

#ifndef IR_META_2I
#define IR_META_2I

#include "CxTypes.h"
#include <list>
#include <string> // for Desc() output
#include "Ir.h"
#ifdef IR_ECP
   #include "IrEcp.h"
#endif
#include "CxPodArray.h"
#include "CxVec3.h"
#include "CxOpenMpAcc.h"

namespace ct {
typedef
   TVector3<double> FVector3;


// IR driver-level interface to a basis set. Semi-lightweight.
// supposed to be generated from the program's native basis set
// format (FBasisSet here) on demand for interfacing between IR
// and the program.
struct FRawBasis : public FIntrusivePtrDest
{
   typedef TArray<ir::FRawShell>
      FShellArray;
   FShellArray
      Shells;
   TArray<size_t>
      // offset of Shell[i] with respect to underlying basis in terms of
      // individual functions. Contains one additional element giving the
      // total number of functions.
      ShellOffsets,
      // shells of atom #i start at CenterOffsets[i].
      CenterOffsets,
      // [i]: center index of shell [i].
      ShellCenters;
   size_t
      // largest number of functions in a shell occuring in the basis.
      // May be useful for allocating intermediates.
      nMaxFnPerShell;
   size_t nFn() const { return ShellOffsets.back(); }
   size_t nSh() const { return Shells.size(); }
   size_t nFn(size_t iSh) const { return ShellOffsets[iSh+1]-ShellOffsets[iSh]; }
   size_t iFn(size_t iSh) const { return ShellOffsets[iSh]; }
   // index of first shell on center iCen
   size_t iCenSh(size_t iCen) const { return CenterOffsets[iCen]; }
   size_t nCenSh(size_t iCen) const { return CenterOffsets[iCen+1] - CenterOffsets[iCen]; }
   // index of first function center iCen
   size_t iCenFn(size_t iCen) const { return iFn(iCenSh(iCen)); }
   size_t nCenFn(size_t iCen) const { return iCenFn(iCen+1) - iCenFn(iCen); }
   // total number of centers
   size_t nCen() const { return CenterOffsets.size() - 1; }

   // center index of a given shell.
   size_t iShCen(size_t iSh) const { return ShellCenters[iSh]; }

   // construct derived data (ShellOffsets, CenterOffsets, nMaxFnPerShell) from
   // this->Shells and the shell center indices given as argument. pShCen[iSh]
   // should give the center id for Shells[iSh], and those must be consecutive.
   void Finalize(size_t const *pShCen);
};
typedef TIntrusivePtr<FRawBasis>
   FRawBasisPtr;

// generalized 2-index integral kernel: effectively a functor which
// can calculate shell doublets for two shell objects. Used for making
// one-electron matrices.
struct FKrn2i : public FIntrusivePtrDest
{
   FKrn2i() {}
   virtual ~FKrn2i();
   virtual void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const = 0;
   // Accumulate derivative contraction
   //     Grad[Nuc,xyz] += Prefactor * \sum_{ac} rdm[a,c] d/d[Nuc,xyz] (a|krn|c).
   // Parameters:
   //   - pGrad[3*iCenter + ixyz]: Location in which the output gradient for center iCenter and ixyz=0...2 is stored.
   //   - pRdmAC[StrideA * ia + StrideC * ic]: gradient coefficients for *pA, *pC,
   //     where ia/ic go over the contracted functions of *pA, *pC.
   //   - iCenA, iCenC: center indices in pGrad itself which belong
   // default implementation just crashes and says "can't do this".
   virtual void AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;
   virtual void AccCoreDeriv2d(double *pHess, size_t iStHessA, size_t iStHessC, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;

   virtual std::string Desc() const = 0;
private:
   FKrn2i(FKrn2i const&); // not implemented
   void operator = (FKrn2i const&); // not implemented
};
typedef TIntrusivePtr<FKrn2i>
   FKrn2iPtr;


// evaluate <a|krn Laplace^LaplaceOrder|b> by directly calling
// a IR integral kernel.
struct FKrn2i_Direct : public FKrn2i
{
   FKrn2i_Direct(ir::FIntegralKernel const *pIrKernel_, uint LaplaceOrder=0, double Prefactor=1.);
   ~FKrn2i_Direct();
   void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
//    void EvalInt2e2c1d(double *pOutDa, double *pOutDb, size_t StrideA, size_t StrideC, size_t StrideDeriv, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
   void AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;
   std::string Desc() const; // override
protected:
   ir::FIntegralKernel const
      *m_pIrKernel;
   uint
      m_LaplaceOrder;
   double
      m_Prefactor;
};

// this kernel = linear combination of other kernel objects.
struct FKrn2i_Meta : public FKrn2i
{
   FKrn2i_Meta();
   ~FKrn2i_Meta();
   void AddKernel(FKrn2iPtr pKrn2i, double fPrefactor);
   void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
   void AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;
   std::string Desc() const; // override
protected:
   struct FKernelEntry {
      FKrn2iPtr
         p;
      double
         Factor;
      FKernelEntry() {}
      FKernelEntry(FKrn2iPtr p_, double Factor_) : p(p_), Factor(Factor_) {}
   };
   typedef std::list<FKernelEntry>
      FKernelList;
   FKernelList
      m_Kernels;
};
typedef TIntrusivePtr<FKrn2i_Meta>
   FKrn2iMetaPtr;

// evaluate <a| v(r) |b> where v(r) is sum of the fields emitted via pIrKernel
// for the sum of the supplied spherical point multipoles.
// Can be used, for example, to make nuclear attraction integrals.
struct FKrn2i_PointMultipoles : public FKrn2i
{
   // each point should provide 2*l+1 coefficients; one for each of its spherical
   // components. Slm order is equal to the order used for basis functions.
   // (in particular, dipole components are x,y,z, not x,z,y.).
   struct FPoint {
      FVector3 vPos;
      uint l;
      ptrdiff_t iCenter;
      // ^- for gradients: identify the atomic center index on which this multipole sits. Set to -1 if None.
      //    note: even if having external lattices, one could just append their gradient array to the QM system gradient array
      //    and then do both in one go.
   };
   FKrn2i_PointMultipoles(ir::FIntegralKernel const *pIrKernel_, FPoint const *pPoints, double const *pCoeffs, size_t nPoints);
   ~FKrn2i_PointMultipoles();
   void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
// #ifdef INCLUDE_OPTIONALS
   void AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;
   void AccCoreDeriv2d(double *pHess, size_t iStHessA, size_t iStHessC, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;
// #endif // INCLUDE_OPTIONALS
   enum FDerivativeFlags {
      DERIVFLAG_IncludeShellDerivsA = 0x01, // if set, include derivatives with respect to the row basis function centers (A)
      DERIVFLAG_IncludeShellDerivsB = 0x02, // if set, include derivatives with respect to the col basis function centers (B)
      DERIVFLAG_IncludePointDerivsP = 0x04, // if set, include derivatives with respect to the point multipole centers (C)
      DERIVFLAG_IncludeShellDerivsAB = DERIVFLAG_IncludeShellDerivsA | DERIVFLAG_IncludeShellDerivsB, // if set, include derivatives with respect to the basis function centers (AB)
      DERIVFLAG_AbcMask = DERIVFLAG_IncludeShellDerivsA | DERIVFLAG_IncludeShellDerivsB | DERIVFLAG_IncludePointDerivsP
   };
   void SetDerivativeFlags(unsigned Flags);
   std::string Desc() const; // override
protected:
   ir::FIntegralKernel const
      *m_pIrKernel;
   TArray<ir::FRawShell>
      m_PointShells;
   TArray<double>
      m_Data;
   TArray<ptrdiff_t>
      m_PointCenters;
   double
      m_Exp;
   uint
      m_MaxL;
   unsigned
      // bit field of DERIVFLAG_*
      m_DerivativeFlags;
   bool bDerivA() const { return 0 != ((m_DerivativeFlags & DERIVFLAG_AbcMask) & DERIVFLAG_IncludeShellDerivsA); }
   bool bDerivB() const { return 0 != ((m_DerivativeFlags & DERIVFLAG_AbcMask) & DERIVFLAG_IncludeShellDerivsB); }
   bool bDerivP() const { return 0 != ((m_DerivativeFlags & DERIVFLAG_AbcMask) & DERIVFLAG_IncludePointDerivsP); }
};

#ifdef IR_ECP
struct FKrn2i_1CenterEcp : public FKrn2i
{
   FKrn2i_1CenterEcp(FVector3 vEcpCenter, ptrdiff_t iEcpCenter, ir::FAtomEcp const *pIrEcp);
   ~FKrn2i_1CenterEcp();
   void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
// #ifdef INCLUDE_OPTIONALS
   void AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const;
// #endif // INCLUDE_OPTIONALS
   std::string Desc() const; // override
protected:
   ir::FAtomEcp const
      *m_pIrEcp;
   ptrdiff_t
      m_iEcpCenter;
   FVector3
      m_vEcpCenter;
};
#endif // IR_ECP


void MakeIntMatrix( double *pDest, size_t iRowSt, size_t iColSt, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor, bool Add );
void AccGradient2ix(double *pGrad, double *pRdm, size_t iRowSt, size_t iColSt, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor);
void AccHessian2ix(double *pHess, size_t iRowStH, size_t iColStH,  double *pRdm, size_t iRowStD, size_t iColStD,
   FRawBasis const &BasisRow, FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor);


} // namespace ct

#endif // IR_META_2I
