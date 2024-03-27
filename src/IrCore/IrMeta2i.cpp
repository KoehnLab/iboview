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
#include "IrMeta2i.h"
#include "IrInternal.h"
// #include "CxOpenMpProxy.h"

#ifdef IR_TIMING
   #include <stdio.h>
   #include "format.h"
#endif // IR_TIMING

namespace ct {

FKrn2i::~FKrn2i()
{}

void FKrn2i::AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(pGrad);
   IR_SUPPRESS_UNUSED_WARNING(pRdmAC);
   IR_SUPPRESS_UNUSED_WARNING(StrideA);
   IR_SUPPRESS_UNUSED_WARNING(StrideC);
   IR_SUPPRESS_UNUSED_WARNING(pA);
   IR_SUPPRESS_UNUSED_WARNING(iCenterA);
   IR_SUPPRESS_UNUSED_WARNING(pC);
   IR_SUPPRESS_UNUSED_WARNING(iCenterC);
   IR_SUPPRESS_UNUSED_WARNING(Prefactor);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   throw std::runtime_error("FKrn2i: 1st derivative integrals apparently not implemented for current integral type. Overwrite AccCoreDeriv1d.");
}


void FKrn2i::AccCoreDeriv2d(double *pHess, size_t iStHessA, size_t iStHessC, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(pHess);
   IR_SUPPRESS_UNUSED_WARNING(iStHessA);
   IR_SUPPRESS_UNUSED_WARNING(iStHessC);
   IR_SUPPRESS_UNUSED_WARNING(pRdmAC);
   IR_SUPPRESS_UNUSED_WARNING(StrideA);
   IR_SUPPRESS_UNUSED_WARNING(StrideC);
   IR_SUPPRESS_UNUSED_WARNING(pA);
   IR_SUPPRESS_UNUSED_WARNING(iCenterA);
   IR_SUPPRESS_UNUSED_WARNING(pC);
   IR_SUPPRESS_UNUSED_WARNING(iCenterC);
   IR_SUPPRESS_UNUSED_WARNING(Prefactor);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   throw std::runtime_error("FKrn2i: 2nd derivative integrals apparently not implemented for current integral type. Overwrite AccCoreDeriv2d.");
}



std::string FKrn2i_PointMultipoles::Desc() const
{
   return "vext";
}

FKrn2i_PointMultipoles::FKrn2i_PointMultipoles(ir::FIntegralKernel const *pIrKernel_, FPoint const *pPoints, double const *pCoeffs, size_t nPoints)
   : m_pIrKernel(pIrKernel_), m_Exp(1e20)
{
   m_PointShells.reserve(nPoints);
   m_PointCenters.reserve(nPoints);
   // count number of coefficients and find max angular momentum.
   size_t
      nCoeff = 0;
   m_MaxL = 0;
   for (size_t iPt = 0; iPt < nPoints; ++ iPt) {
      nCoeff += 2*pPoints[iPt].l + 1;
      m_MaxL = std::max(m_MaxL, pPoints[iPt].l);
   }
   // store the point centers and their strenghts here.
   m_Data.resize(3*nPoints + nCoeff);
   // copy centers
   for (size_t iPt = 0; iPt < nPoints; ++ iPt)
      for (uint i = 0; i < 3; ++ i)
         m_Data[3*iPt + i] = pPoints[iPt].vPos[i];
   // copy coefficients and make RawShell objects.
   std::size_t
      iCoeff = 0;
   for (size_t iPt = 0; iPt < nPoints; ++ iPt) {
      FPoint const &Pt = pPoints[iPt];
      m_PointShells.push_back(ir::FRawShell(Pt.l, &m_Exp, 1,
         &m_Data[3*nPoints+iCoeff], 1, &m_Data[3*iPt]));
      m_PointCenters.push_back(Pt.iCenter);

      double
         // This works for point charges. I think it should be the
         // right factor for higher multipoles, too.
         fInvNorm = std::pow(m_Exp/M_PI, 3./2.);
      if (Pt.l != 0)
         throw std::runtime_error("case Pt.l != 0 not tested. Not 100% sure if norm factor is correct!");
      for (size_t iSh = 0; iSh < 2*Pt.l + 1; ++ iSh)
         m_Data[3*nPoints + iCoeff + iSh] = fInvNorm * pCoeffs[iCoeff + iSh];
      iCoeff += 2*Pt.l+1;
   }
   assert(iCoeff == nCoeff);

   // by default construct all derivatives.
   SetDerivativeFlags(DERIVFLAG_IncludeShellDerivsAB | DERIVFLAG_IncludePointDerivsP);
}

FKrn2i_PointMultipoles::~FKrn2i_PointMultipoles()
{}


void FKrn2i_PointMultipoles::SetDerivativeFlags(unsigned Flags_)
{
   m_DerivativeFlags = Flags_;
}


void FKrn2i_PointMultipoles::EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideB, ir::FRawShell const *pA, ir::FRawShell const *pB, double Prefactor, bool Add, ct::FMemoryStack &Mem) const
{
   if (!Add)
      for (size_t iB = 0; iB != pB->nFn(); ++ iB)
         for (size_t iA = 0; iA != pA->nFn(); ++ iA)
            pOut[StrideA*iA + StrideB*iB] = 0.;
   double
      *pIntData = Mem.AllocN(pA->nFn() * pB->nFn() * (2*m_MaxL+1), (double)0.);
   for (size_t iPt = 0; iPt < m_PointShells.size(); ++ iPt) {
      size_t
         Strides3[3] = {1, pA->nFn(), pA->nFn()*pB->nFn()};
      // TODO:
      //   - do this with a specialized routine which takes ZetaC == +inf
      //     into account. Can be done much cheaper then.
      //   - use the C-inline-contracting version of EvalInt2e3c: EvalInt2e3c_ContractC.
      //     we actually have it implemented!!
      ir::FRawShell const
         *pShPt = &m_PointShells[iPt];
      EvalInt2e3c(pIntData, Strides3,
         pA, pB, pShPt, 1, Prefactor, m_pIrKernel, Mem);
      for (size_t iB = 0; iB != pB->nFn(); ++ iB)
         for (size_t iA = 0; iA != pA->nFn(); ++ iA)
            for (size_t iCo = 0; iCo < pShPt->nSh(); ++ iCo)
               pOut[iA*StrideA + iB*StrideB] += pIntData[iA*Strides3[0] + iB*Strides3[1] + iCo*Strides3[2]];
   }
   Mem.Free(pIntData);
}

// #ifdef INCLUDE_OPTIONALS
void FKrn2i_PointMultipoles::AccCoreDeriv1d(double *pGrad, double const *pRdmAB, size_t StrideA, size_t StrideB, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pB, size_t iCenterB, double Prefactor, ct::FMemoryStack &Mem) const
{
   double
      *pIntData = Mem.AllocN(pA->nFn() * pB->nFn() * (2*m_MaxL+1) * 6, (double)0.);
   size_t
      nFnA = pA->nFn(),
      nFnB = pB->nFn();
   for (size_t iPt = 0; iPt < m_PointShells.size(); ++ iPt) {
      ir::FRawShell const
         *pShPt = &m_PointShells[iPt];
      size_t
         nFnPt = pShPt->nFn(),
         Strides4[4] = {1, nFnA, nFnA*nFnB, nFnA*nFnB*nFnPt };
//       memset(pIntData, 0, sizeof(double) * pA->nFn() * pB->nFn() * 6);
      // FIXME: use inline contracting version (see comments in FKrn2i_PointMultipoles::EvalInt2e2c)
      EvalInt2e3c1d(pIntData, Strides4,
         pA, pB, pShPt, 1, Prefactor, m_pIrKernel, Mem);

      ptrdiff_t iCenterPt = m_PointCenters[iPt];
//       xout << format("iCenterPt = %i\n") % iCenterPt;
      for (size_t iFnB = 0; iFnB != nFnB; ++ iFnB)
         for (size_t iFnA = 0; iFnA != nFnA; ++ iFnA) {
            double RdmAB = pRdmAB[StrideA * iFnA + StrideB * iFnB];
            for (size_t iFnPt = 0; iFnPt < pShPt->nSh(); ++ iFnPt) {
               for (size_t ixyz = 0; ixyz != 3; ++ ixyz) {
                  double dA = RdmAB * pIntData[iFnA + nFnA*(iFnB + nFnB*(iFnPt + nFnPt*(0 + ixyz)))];
                  double dB = RdmAB * pIntData[iFnA + nFnA*(iFnB + nFnB*(iFnPt + nFnPt*(3 + ixyz)))];
                  double dPt = -dA - dB;
                  if (bDerivA())
                     pGrad[3*iCenterA + ixyz] += dA;
                  if (bDerivB())
                     pGrad[3*iCenterB + ixyz] += dB;
                  if (iCenterPt != -1 && bDerivP())
                     pGrad[3*iCenterPt + ixyz] += dPt;
               }
            }
         }
   }
   Mem.Free(pIntData);
}


void FKrn2i_PointMultipoles::AccCoreDeriv2d(double *pHess, size_t iStHessA, size_t iStHessC, double const *pRdmAB, size_t StrideA, size_t StrideB, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pB, size_t iCenterB, double Prefactor, ct::FMemoryStack &Mem) const
{
   if (bDerivA() || bDerivB() || !bDerivP())
      throw std::runtime_error("sorry, AccCoreDeriv2d not quite correctly implemented for the non-C case.");
   enum {
      iAxx,iAyy,iAzz,iAxy,iAxz,iAyz,
      iAxBx,iAyBx,iAzBx,  iAxBy,iAyBy,iAzBy,  iAxBz,iAyBz,iAzBz,
      iBxx,iByy,iBzz,iBxy,iBxz,iByz,
   };
   if (pA->l < pB->l) {
      // FIXME: fix swap case in 2d kernel and then remove this.
      std::swap(pA, pB);
      std::swap(iStHessA, iStHessC);
      std::swap(StrideA, StrideB);
      std::swap(iCenterA, iCenterB);
   }
   // translational invariance gives us: (u,v in {x,y,z})
   //    (d[Au] + d[Bu] + d[Cu]) = 0    ->   d[Cu] = -(D[Au] + D[Bu])
   //    (d[Av] + d[Bv] + d[Cv]) = 0    ->   d[Cv] = -(D[Av] + D[Bv])
   //
   // so we can get the 2nd C derivatives as:
   //
   // -> d[Cu] d[Cv] int = +1. * [(d[Au] + d[Bu]) (d[Av] + d[Bv])] int
   //                    = +1. * [d[Au] d[Av] + d[Av] d[Bu] + d[Au] d[Bv] + d[Bu] d[Bv]] int
   double
      *pIntData = Mem.AllocN(pA->nFn() * pB->nFn() * (2*m_MaxL+1) * (6+9+6), (double)0.);
   size_t
      nFnA = pA->nFn(),
      nFnB = pB->nFn();
   for (size_t iPt = 0; iPt < m_PointShells.size(); ++ iPt) {
      ir::FRawShell const
         *pShPt = &m_PointShells[iPt];
      size_t
         nFnPt = pShPt->nFn(),
         Strides4[4] = {1, nFnA, nFnA*nFnB, nFnA*nFnB*nFnPt };
//       memset(pIntData, 0, sizeof(double) * pA->nFn() * pB->nFn() * 6);
      // FIXME: use inline contracting version (see comments in FKrn2i_PointMultipoles::EvalInt2e2c)
      EvalInt2e3c2d(pIntData, Strides4,
         pA, pB, pShPt, 1, Prefactor, m_pIrKernel, Mem);

      ptrdiff_t iCenterPt = m_PointCenters[iPt];
//       xout << format("iCenterPt = %i\n") % iCenterPt;
      for (size_t iFnB = 0; iFnB != nFnB; ++ iFnB)
         for (size_t iFnA = 0; iFnA != nFnA; ++ iFnA) {
            double RdmAB = pRdmAB[StrideA * iFnA + StrideB * iFnB];
            for (size_t iFnPt = 0; iFnPt < pShPt->nSh(); ++ iFnPt) {
               if (iCenterPt == -1)
                  continue;
               size_t
                  iStD = Strides4[3];
               double
                  *pD = &pIntData[iFnA + nFnA*(iFnB + nFnB*(iFnPt + nFnPt))],
                  *pH = &pHess[iStHessA * 3*iCenterPt + iStHessC * 3*iCenterPt];
                  // ^- yes, this has block-diagonal terms only in C-only mode.
#define D(i) pD[iStD*(i)]
               pH[iStHessA*0 + 0*iStHessC] += RdmAB*( D(iAxx) + D(iAxBx) + D(iAxBx) + D(iBxx) );
               pH[iStHessA*1 + 0*iStHessC] += RdmAB*( D(iAxy) + D(iAxBy) + D(iAyBx) + D(iBxy) );
               pH[iStHessA*2 + 0*iStHessC] += RdmAB*( D(iAxz) + D(iAzBx) + D(iAxBz) + D(iBxz) );

               pH[iStHessA*0 + 1*iStHessC] += RdmAB*( D(iAxy) + D(iAxBy) + D(iAyBx) + D(iBxy) );
               pH[iStHessA*1 + 1*iStHessC] += RdmAB*( D(iAyy) + D(iAyBy) + D(iAyBy) + D(iByy) );
               pH[iStHessA*2 + 1*iStHessC] += RdmAB*( D(iAyz) + D(iAyBz) + D(iAzBy) + D(iByz) );

               pH[iStHessA*0 + 2*iStHessC] += RdmAB*( D(iAxz) + D(iAzBx) + D(iAxBz) + D(iBxz) );
               pH[iStHessA*1 + 2*iStHessC] += RdmAB*( D(iAyz) + D(iAyBz) + D(iAzBy) + D(iByz) );
               pH[iStHessA*2 + 2*iStHessC] += RdmAB*( D(iAzz) + D(iAzBz) + D(iAzBz) + D(iBzz) );
#undef D
            }
         }
   }
   Mem.Free(pIntData);
}

// #endif // INCLUDE_OPTIONALS



FKrn2i_Direct::FKrn2i_Direct(ir::FIntegralKernel const *pIrKernel_, uint LaplaceOrder, double Prefactor)
   : m_pIrKernel(pIrKernel_), m_LaplaceOrder(LaplaceOrder), m_Prefactor(Prefactor)
{}

FKrn2i_Direct::~FKrn2i_Direct()
{}

std::string FKrn2i_Direct::Desc() const
{
   if (m_LaplaceOrder == 0) {
      return m_pIrKernel->Desc();
   } else if (m_LaplaceOrder == 1) {
      return "Laplace " + std::string(m_pIrKernel->Desc());
   } else {
      return "Laplace^n " + std::string(m_pIrKernel->Desc());
   }
}

void FKrn2i_Direct::EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const
{
   ir::EvalInt2e2c_LaplaceC(pOut, StrideA, StrideC, pA, pC, m_Prefactor * Prefactor, Add,
      m_LaplaceOrder, m_pIrKernel, Mem);
}

void FKrn2i_Direct::AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const
{
   void
      *pBaseOfMemory = Mem.Alloc(0);
   size_t
      nFnA = pA->nFn(),
      nFnC = pC->nFn();
   double
      *pDerivA, *pDerivC;
   Mem.Alloc(pDerivA, nFnA * nFnC * 3);
   Mem.Alloc(pDerivC, nFnA * nFnC * 3);

   // make the actual gradient integrals (yes, we actually calculate them explcitly.).
   ir::EvalInt2e2c1d_LaplaceC(pDerivA, pDerivC, 1, nFnA, nFnA*nFnC, pA, pC, m_Prefactor * Prefactor, false, m_LaplaceOrder, m_pIrKernel, Mem);

   // contract them with their section of the density matrix.
   for (size_t iFnC = 0; iFnC < nFnC; ++ iFnC)
      for (size_t iFnA = 0; iFnA < nFnA; ++ iFnA) {
         double RdmAC = pRdmAC[StrideA * iFnA + StrideC * iFnC];
         for (size_t ixyz = 0; ixyz != 3; ++ ixyz) {
            pGrad[3*iCenterA + ixyz] += RdmAC * pDerivA[iFnA + nFnA*(iFnC + nFnC * ixyz)];
            pGrad[3*iCenterC + ixyz] += RdmAC * pDerivC[iFnA + nFnA*(iFnC + nFnC * ixyz)];
//             pGrad[3*iCenterC + ixyz] -= RdmAC * pDerivA[iFnA + nFnA*(iFnC + nFnC * ixyz)];
         }
      }

   Mem.Free(pBaseOfMemory);
}

FKrn2i_Meta::FKrn2i_Meta()
{}

FKrn2i_Meta::~FKrn2i_Meta()
{}

std::string FKrn2i_Meta::Desc() const
{
   return "2iMeta";
}


void FKrn2i_Meta::AddKernel(FKrn2iPtr pKrn2i, double fPrefactor)
{
   if (fPrefactor != 0.)
      m_Kernels.push_back(FKernelEntry(pKrn2i, fPrefactor));
}


void FKrn2i_Meta::EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const
{
   FKernelList::const_iterator
      itKrn2i;
   for (itKrn2i = m_Kernels.begin(); itKrn2i != m_Kernels.end(); ++ itKrn2i) {
      itKrn2i->p->EvalInt2e2c(pOut, StrideA, StrideC, pA, pC, Prefactor * itKrn2i->Factor, Add, Mem);
      Add = true;
   }
   if (!Add) {
      // if we reached this point and 'Add' is still false, then we have
      // an empty kernel list and no EvalInt2e2c was actually called
      // (and therefore nothing overwrite the target location).
      // That means we still have to zero out the target memory now
      // before we return.
      ir::ZeroBlock2(pOut, StrideA, StrideC, pA->nFn(), pC->nFn());
   }
}


void FKrn2i_Meta::AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const
{
   FKernelList::const_iterator
      itKrn2i;
   for (itKrn2i = m_Kernels.begin(); itKrn2i != m_Kernels.end(); ++ itKrn2i) {
      itKrn2i->p->AccCoreDeriv1d(pGrad, pRdmAC, StrideA, StrideC, pA, iCenterA, pC, iCenterC, Prefactor * itKrn2i->Factor, Mem);
   }
}



#ifdef IR_ECP
FKrn2i_1CenterEcp::FKrn2i_1CenterEcp(FVector3 vEcpCenter, ptrdiff_t iEcpCenter, ir::FAtomEcp const *pIrEcp)
   : m_pIrEcp(pIrEcp), m_iEcpCenter(iEcpCenter), m_vEcpCenter(vEcpCenter)
{}

FKrn2i_1CenterEcp::~FKrn2i_1CenterEcp()
{}

std::string FKrn2i_1CenterEcp::Desc() const
{
   return "ecp1c";
}

void FKrn2i_1CenterEcp::EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const
{
// #ifdef INCLUDE_OPTIONALS
   ir::EvalEcp(pOut, StrideA, StrideC, pA, pC, Prefactor, Add, m_pIrEcp, &m_vEcpCenter[0], Mem);
   return;
// #endif // INCLUDE_OPTIONALS
   throw std::runtime_error("This version of MicroScf was compiled without ECP support.");
}

// #ifdef INCLUDE_OPTIONALS
void FKrn2i_1CenterEcp::AccCoreDeriv1d(double *pGrad, double const *pRdmAC, size_t StrideA, size_t StrideC, ir::FRawShell const *pA, size_t iCenterA, ir::FRawShell const *pC, size_t iCenterC, double Prefactor, ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(pGrad);
   IR_SUPPRESS_UNUSED_WARNING(pRdmAC);
   IR_SUPPRESS_UNUSED_WARNING(StrideA);
   IR_SUPPRESS_UNUSED_WARNING(StrideC);
   IR_SUPPRESS_UNUSED_WARNING(pA);
   IR_SUPPRESS_UNUSED_WARNING(iCenterA);
   IR_SUPPRESS_UNUSED_WARNING(pC);
   IR_SUPPRESS_UNUSED_WARNING(iCenterC);
   IR_SUPPRESS_UNUSED_WARNING(Prefactor);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   throw std::runtime_error("ECP derivatives not implemented.");
}
// #endif // INCLUDE_OPTIONALS
#endif // IR_ECP

void FRawBasis::Finalize(size_t const *pShCen)
{
   ShellOffsets.clear();
   CenterOffsets.clear();
   ShellCenters.clear();

   ShellOffsets.reserve(Shells.size() + 1);
   ShellOffsets.push_back(0);
   CenterOffsets.reserve(Shells.size() + 1);
   ShellCenters.reserve(Shells.size());
   nMaxFnPerShell = 0;
   size_t
      iShellOff = 0;
   for (size_t iSh = 0; iSh < Shells.size(); ++ iSh){
      while (pShCen[iSh] >= (size_t)CenterOffsets.size())
         CenterOffsets.push_back(iShellOff);
      iShellOff += 1;
      size_t
         nShFn = Shells[iSh].nFn();
      ShellOffsets.push_back(ShellOffsets.back() + nShFn);
      nMaxFnPerShell = std::max((size_t)nMaxFnPerShell, nShFn);
      ShellCenters.push_back(pShCen[iSh]);
   }
   assert(ShellOffsets.back() == this->nFn());
   CenterOffsets.push_back(Shells.size());
}


#ifdef IR_TIMING
   static double s_Ir2ixTiming_t0;
#endif // RDTSC_TIMING

static void IrInit2ixTiming(char const *pFuncDesc, FRawBasis const &BasisRow, FRawBasis const &BasisCol, FKrn2i const &Krn2i)
{
#ifdef IR_TIMING
   ResetClocks();
//    printf("\n-- IR-TIMING: 2ix matrix <row|%s|col> shape=(%i, %i) nOmpThreads=%i\n", Krn2i.Desc().c_str(), int(BasisRow.nFn()), int(BasisCol.nFn()), omp_get_max_threads());
   printf("\n-- IR-TIMING: %s <row|%s|col> shape=(%i, %i)\n", pFuncDesc, Krn2i.Desc().c_str(), int(BasisRow.nFn()), int(BasisCol.nFn()));
   s_Ir2ixTiming_t0 = 0;
   s_Ir2ixTiming_t0 -= ir::GetTicks1();
#else
   IR_SUPPRESS_UNUSED_WARNING(pFuncDesc);
   IR_SUPPRESS_UNUSED_WARNING(BasisRow);
   IR_SUPPRESS_UNUSED_WARNING(BasisCol);
   IR_SUPPRESS_UNUSED_WARNING(Krn2i);
#endif // IR_TIMING
}


static void IrFinish2ixTiming(char const *pFuncDesc, FRawBasis const &BasisRow, FRawBasis const &BasisCol, FKrn2i const &Krn2i)
{
#ifdef IR_TIMING
   s_Ir2ixTiming_t0 += ir::GetTicks1();
   using namespace ir;
   InitTimingCode();
//    PrintClock("EvalInt2e2c(get)", TscGetClock(TID_EvalInt2e2c), TID_EvalInt2e2c);
   SetRefClock(4);
   PrintClock("EvalInt2e2c", TID_EvalInt2e2c);
//    printf("^- I don't get it... if I divide CoShY by its fraction on the ref clock, i get the right number. But not if doing it raw.");
   // ^- UPDATE: ...that's because the first time it runs, it does the
   //    CPU frequency and overhead test. GNAAAA.
//    SetRefClock(TID_EvalInt2e2c);
   PrintClock("| CoShY", TID_EvalCoShY);
   PrintClock("| | Setup", TID_EvalCoShY_Setup);
   PrintClock("| | CoKernels", TID_EvalCoKernels);
   FTscTime tPTC = TscGetClock(TID_EvalCoKernels_PrimData) + TscGetClock(TID_EvalCoKernels_TildeGm) + TscGetClock(TID_EvalCoKernels_Contract);
   PrintClock("| | | (PTC)", tPTC, TID_EvalCoKernels);
   PrintClock("| | | PrimData", TID_EvalCoKernels_PrimData);
   PrintClock("| | | TildeGm", TID_EvalCoKernels_TildeGm);
   PrintClock("| | | Contract", TID_EvalCoKernels_Contract);
   PrintClock("| | MDRR", TID_EvalCoShY_MDRR);
   PrintClock("| ShTrA_YY", TID_ShTrA_YY);
   PrintClock("| Scatter2e2c", TID_Scatter2e2c);

   SetRefClock(-1);
   printf("   %s (total, real): %31.3f msec\n", pFuncDesc, 1000.*s_Ir2ixTiming_t0);
   PrintClock("Operation (total, tsc)", 4);
   SetRefClock(4);
   PrintClock("Krn2i.EvalInt2e2c (tot)", 1);
   PrintClock("write output (in-loop)", 2);
   PrintClock("write output (sym fix/grad join)", 3);
   SetRefClock(-1);
#else
   IR_SUPPRESS_UNUSED_WARNING(pFuncDesc);
   IR_SUPPRESS_UNUSED_WARNING(BasisRow);
   IR_SUPPRESS_UNUSED_WARNING(BasisCol);
   IR_SUPPRESS_UNUSED_WARNING(Krn2i);
#endif // IR_TIMING
}


void MakeIntMatrix( double *pDest, size_t iRowSt, size_t iColSt, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor, bool Add )
{
   bool
      MatrixSymmetric = (&BasisRow == &BasisCol);

   FMemoryStackArray MemStacks(Mem_);

//    IrInit2ixTiming("2ix matrix", BasisRow, BasisCol, Krn2i);
   IrInit2ixTiming("MakeIntMatrix", BasisRow, BasisCol, Krn2i);
   IR_RESUME_CLOCK(4);

//          Mem.Alloc(pIntResult, BasisRow.nMaxFnPerShell * BasisCol.nMaxFnPerShell);
   #pragma omp parallel for schedule(dynamic)
   for (int iShellB_ = 0; iShellB_ < int(BasisCol.Shells.size()); ++ iShellB_) {
      FMemoryStack &Mem = MemStacks.GetStackOfThread();
      size_t iShellB = size_t(iShellB_);
      ir::FRawShell const
         *pShellB = &BasisCol.Shells[iShellB];
      for (size_t iShellA = (!MatrixSymmetric)? 0 : iShellB;
           iShellA < BasisRow.Shells.size(); ++ iShellA)
      {
         ir::FRawShell const
            *pShellA = &BasisRow.Shells[iShellA];
         size_t
            nSizeA = pShellA->nFn(),
            nSizeB = pShellB->nFn();
         assert(nSizeA == BasisRow.nFn(iShellA) &&
                nSizeB == BasisCol.nFn(iShellB));
         TMemoryLock<double>
            pIntResult(nSizeA * nSizeB, &Mem);

         IR_RESUME_CLOCK(1);
         Krn2i.EvalInt2e2c(pIntResult, 1, nSizeA, pShellA, pShellB, Prefactor, false, Mem);
         IR_PAUSE_CLOCK(1);

         // fill data we gathered into the matrices
         IR_RESUME_CLOCK(2);
         for ( size_t j_ = 0; j_ < nSizeB; ++ j_ ) {
            for ( size_t i_ = 0; i_ < nSizeA; ++ i_ ) {
               size_t
                  i = BasisRow.iFn(iShellA) + i_,
                  j = BasisCol.iFn(iShellB) + j_;
               double const
                  &r = pIntResult[i_ + nSizeA * j_];
               if (!Add) {
                  pDest[i*iRowSt + j*iColSt] = r;
               } else {
                  pDest[i*iRowSt + j*iColSt] += r;
               }
            }
         }
         IR_PAUSE_CLOCK(2);
      }
   }
   IR_RESUME_CLOCK(3);
   if (MatrixSymmetric) {
      // copy newly computed triangle to other triangle (this version should have continuous
      // writes for standard layout matrices).
//       size_t nFn = BasisRow.nFn();
//       for (size_t i = 0; i < nFn; ++ i)
//          for (size_t j = 0; j < i; ++ j)
//             pDest[j*iRowSt + i*iColSt] = pDest[i*iRowSt + j*iColSt];
      #pragma omp parallel for schedule(dynamic)
      for (int iShellB_ = 0; iShellB_ < int(BasisCol.Shells.size()); ++ iShellB_) {
         size_t
            iShellB = size_t(iShellB_);
         size_t
            iFnB_Beg = BasisCol.ShellOffsets[iShellB],
            iFnB_End = BasisCol.ShellOffsets[iShellB+1];
         for (size_t i = iFnB_Beg; i < iFnB_End; ++ i)
            for (size_t j = 0; j < i; ++ j)
               pDest[j*iRowSt + i*iColSt] = pDest[i*iRowSt + j*iColSt];
      }
   }
   IR_PAUSE_CLOCK(3);

   IR_PAUSE_CLOCK(4);
   IrFinish2ixTiming("MakeIntMatrix", BasisRow, BasisCol, Krn2i);

}


void AccGradient2ix(double *pGrad, double *pRdm, size_t iRowSt, size_t iColSt, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor)
{
   bool
      MatrixSymmetric = (&BasisRow == &BasisCol);

   IrInit2ixTiming("AccGradient2ix", BasisRow, BasisCol, Krn2i);
   IR_RESUME_CLOCK(4);

   assert(BasisRow.nCen() == BasisCol.nCen());
   FOmpAccBlock
      Gradient(pGrad, 3*BasisRow.nCen(), 0, Mem_);
   FMemoryStackArray MemStacks(Mem_);

   #pragma omp parallel for schedule(dynamic)
   for (int iShellA_ = 0; iShellA_ < int(BasisRow.Shells.size()); ++ iShellA_) {
      size_t iShellA = size_t(iShellA_);
      FMemoryStack &Mem = MemStacks.GetStackOfThread();
      ir::FRawShell const
         *pShellA = &BasisRow.Shells[iShellA];
      for (size_t iShellB = (!MatrixSymmetric)? 0 : iShellA;
           iShellB < BasisCol.Shells.size(); ++ iShellB)
      {
         ir::FRawShell const
            *pShellB = &BasisCol.Shells[iShellB];
         double
            Factor1 = Prefactor;
         if (MatrixSymmetric && iShellA != iShellB)
            Factor1 *= 2.0;
         double const
            *pRdmAB = &pRdm[iRowSt * BasisRow.iFn(iShellA) + iColSt * BasisCol.iFn(iShellB)];
         IR_RESUME_CLOCK(1);
         Krn2i.AccCoreDeriv1d(Gradient.pTls(), pRdmAB, iRowSt, iColSt,
            pShellA, BasisRow.iShCen(iShellA), pShellB, BasisCol.iShCen(iShellB), Factor1, Mem);
         IR_PAUSE_CLOCK(1);
      }
   }
   IR_RESUME_CLOCK(3);
   Gradient.Join();
   IR_PAUSE_CLOCK(3);

   IR_PAUSE_CLOCK(4);
   IrFinish2ixTiming("AccGradient2ix", BasisRow, BasisCol, Krn2i);
}


void AccHessian2ix(double *pHess, size_t iRowStH, size_t iColStH,  double *pRdm, size_t iRowStD, size_t iColStD,
   FRawBasis const &BasisRow, FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem_, double Prefactor)
{
   bool
      MatrixSymmetric = (&BasisRow == &BasisCol);

   assert(BasisRow.nCen() == BasisCol.nCen());
   assert_rt(iColStH == BasisCol.nCen()*3);
//    size_t nHessSize = iColStH * (BasisCol.nCen()*3);

//    FOmpAccBlock
//       Hessian(pHess, nHessSize, 0, Mem_);
//    FMemoryStackArray MemStacks(Mem_);

//    #pragma omp parallel for schedule(dynamic)
   for (int iShellA_ = 0; iShellA_ < int(BasisRow.Shells.size()); ++ iShellA_) {
      size_t iShellA = size_t(iShellA_);
//       FMemoryStack &Mem = MemStacks.GetStackOfThread();
      FMemoryStack &Mem = Mem_;
      ir::FRawShell const
         *pShellA = &BasisRow.Shells[iShellA];
      for (size_t iShellB = (!MatrixSymmetric)? 0 : iShellA;
           iShellB < BasisCol.Shells.size(); ++ iShellB)
      {
         ir::FRawShell const
            *pShellB = &BasisCol.Shells[iShellB];
         double
            Factor1 = Prefactor;
         if (MatrixSymmetric && iShellA != iShellB)
            Factor1 *= 2.0;
         double const
            *pRdmAB = &pRdm[iRowStD * BasisRow.iFn(iShellA) + iColStD * BasisCol.iFn(iShellB)];
//          double
//             *pHessAB = &pHess[iRowStH * 3*BasisRow.iShCen(iShellA) + iColStH * 3*BasisCol.iShCen(iShellB)];
//             *pHessAB = &Hessian.pTls()[iRowStH * BasisRow.iShCen(iShellA) + iColStH * BasisCol.iShCen(iShellB)];
         Krn2i.AccCoreDeriv2d(pHess, iRowStH, iColStH, pRdmAB, iRowStD, iColStD,
            pShellA, BasisRow.iShCen(iShellA), pShellB, BasisCol.iShCen(iShellB), Factor1, Mem);
      }
   }
//    Hessian.Join();
}





} // namespace ct
