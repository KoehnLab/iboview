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

#include <iostream>
#include <cmath>
#include "Ir.h"
#include "CxDiis.h"
#include "CxIo.h"
#include "CxTiming.h"
#include "CxPodArray.h"
#include "CtRhf.h"
#if 0
   #include "RDTSC.h"
#else
   #define RESET_CLOCKS
   #define RESUME_CLOCK(x)
   #define PAUSE_CLOCK(x)
#endif

#include "CxPhysicalUnits.h"
#include "CtDft.h"


#include "CtDfti.h"
#include "CtDftDispCorr.h"

namespace ct {

// double g_fThrOrb = 1e-14;
// double g_fThrOrb = 1e-12;
//    double g_fThrOrb = 1e-10;

double MakeLogThrOrb(double ThrOrb) {
   // hm... this just used a constant 40 before. Seems a bit much.
   // set it to ThrOrb * 1e-3 or something?
   return -std::log(.1 * ThrOrb);
//    return 40.;
}


FXcOptions::FXcOptions(std::string const &XcName, FDftGridParams const *pGridParams_)
{
   if (XcName.empty() || XcName == "none") {
      // note: there is also a "null" functional... but that does a somewhat different
      // thing: it *does* evaluate all the input/output contributions... it just
      // returns 0 for them.
      pXcFn = 0;
   } else {
      pXcFn = new FXcFunctional(XcName);
   }
   pGridParams = pGridParams_;

   fThrOrb = 1e-14;
   UseDfxc = false;
   UseBandStructureEnergyFormula = false;
   // ^-- note: DFXC and band structure formula are set in ScfOptions, not GridOptions.
   //     Therefore not available here.
   GridRenormType = pGridParams->GridRenormType();
}



void InitFittingBasisAndJcd(FBasisSetPtr &pFitBasis, FMatrixView &Jcd, FAtomSet const *pAtoms, FBasisContext FitContext, FHfPrintOptions const *pPrintOptions, FLog *pLog_, FMemoryStack &Mem, char const *pComment=0)
{
   pFitBasis = new FBasisSet(*pAtoms, FitContext);
//    xout << *pFitBasis;
   char const
      *pContextName;
   if (FitContext == BASIS_JFit)
      pContextName = "JFIT";
   else if (FitContext == BASIS_JkFit)
      pContextName = "JKFIT";
   else if (FitContext == BASIS_Mp2Fit)
      pContextName = "MP2FIT";
   else if (FitContext == BASIS_CcsdFit)
      pContextName = "CCFIT";
   else
      pContextName = "???";

   size_t
      nFit = pFitBasis->nFn();
   if (pLog_) {
      pLog_->Write(" Generated {} fitting basis {}{}with {} functions.", pContextName, pFitBasis->Name, pFitBasis->Name.empty()? "" : " ", nFit);
      if (pPrintOptions && (*pPrintOptions)[SCFPRINT_FitBasis]) {
         pLog_->w << "\n" << pFitBasis->GetDesc(FBasisSet::PRINT_Default, " | ") << "\n";
      }
   }
   FTimer
      TimerJcd;
   Jcd = MakeStackMatrix(nFit, nFit, Mem);
   MakeIntMatrix(Jcd, *pFitBasis, *pFitBasis, FKrn2i_Direct(&ir::g_IrCoulombKernel), Mem);
   CalcCholeskyFactors(Jcd);
   if (pLog_) {
      if (pComment)
         pLog_->WriteTiming(fmt::format("fitting metric ({})", pComment), (double)TimerJcd);
      else
         pLog_->WriteTiming("fitting metric", (double)TimerJcd);
   }
}


FFockComponentBuilder::~FFockComponentBuilder()
{}


void FFockComponentBuilder::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem)
{
   size_t
      nBf = pBasis->nFn();
   if (!((FockC.nRows == nBf && FockC.nCols == nBf) &&
         (FockO.pData == 0 || (FockO.nRows == nBf && FockO.nCols == nBf)) &&
         (COccC.nRows == nBf) &&
         (COccO.pData == 0 || COccO.nRows == nBf)))
      throw std::runtime_error("FFockComponentBuilder: Input orbitals not consistent with orbital basis sets.");
   IR_SUPPRESS_UNUSED_WARNING(Flags);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FFockComponentBuilder::PrintEnergyContribs()
{
   m_Log.Write("FFockComponentBuilder::PrintEnergyContribs: not implemented for current Fock builder. Fix this.");
}


void FFockComponentBuilder::AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   m_Log.Write("FFockComponentBuilder::AccGradient: not implemented for current Fock builder. Fix this.");
   IR_SUPPRESS_UNUSED_WARNING(Gradient);
   IR_SUPPRESS_UNUSED_WARNING(COccC);
   IR_SUPPRESS_UNUSED_WARNING(COccO);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


bool FFockComponentBuilder::SwitchToRefineGrid(FDftGridParams const &NewGridParams)
{
   return false;
   IR_SUPPRESS_UNUSED_WARNING(NewGridParams);
}


void FFockComponentBuilder::ComputeBackgroundDensities(FMatrixView CEnvOccC, FMatrixView CEnvOccO, FMemoryStack &Mem)
{
   IR_SUPPRESS_UNUSED_WARNING(CEnvOccC);
   IR_SUPPRESS_UNUSED_WARNING(CEnvOccO);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}



// if i,j are two symmetric quantities to be stored in an array, SymIdx will
// return an unique index for the pair (unique modulo symmetry of course).
// (note: these guys are c/p'd from the 2006 ct8k...)
inline size_t SymIdx(size_t i, size_t j)
{
   if (i > j)
      std::swap(i,j);
   return i + (j*(j+1))/2;
}


inline size_t SymIdx(size_t i, size_t j, size_t k, size_t l)
{
   return SymIdx(SymIdx(i, j), SymIdx(k, l));
}


// if ik = SymIdx(i,k) then (i,k) = UnpackIndex(ik) with i<k
inline void UnpackIndex(size_t &i, size_t &k, size_t ik)
{

   k = static_cast<size_t>(-0.5 + 0.5 * std::sqrt(static_cast<double>(1 + 8*ik)));
   i = ik - (k+1)*k/2;
   assert(SymIdx(i,k) == ik);
   assert(i <= k);
}




// accumulate a coulomb matrix: j[ab] += (ab|cd) d[cd].
void AccCoul(double *pCoulTr, double const *pDensityTr, size_t nAo, double const *IR_RP pIntData, size_t StrideA, size_t StrideB, size_t StrideC, size_t StrideD, size_t iShA, size_t iShB, size_t iShC, size_t iShD, FRawBasis const *pOrbBasis, double fFactor, FMemoryStack &Mem)
{
   size_t
      iFnA = pOrbBasis->iFn(iShA), nFnA = pOrbBasis->nFn(iShA),
      iFnB = pOrbBasis->iFn(iShB), nFnB = pOrbBasis->nFn(iShB),
      iFnC = pOrbBasis->iFn(iShC), nFnC = pOrbBasis->nFn(iShC),
      iFnD = pOrbBasis->iFn(iShD), nFnD = pOrbBasis->nFn(iShD);
   for (size_t iA = 0; iA < nFnA; ++ iA) {
      for (size_t iB = 0; iB < nFnB; ++ iB) {
         for (size_t iC = 0; iC < nFnC; ++ iC) {
            for (size_t iD = 0; iD < nFnD; ++ iD) {
               pCoulTr[SymIdx(iFnA + iA, iFnB + iB)] += fFactor * pIntData[StrideA * iA + StrideB * iB + StrideC * iC + StrideD * iD] * pDensityTr[SymIdx(iFnC + iC, iFnD + iD)];
            }
         }
      }
   }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   IR_SUPPRESS_UNUSED_WARNING(nAo);
}


// accumulate a set of exchange matrices: k_i[ac] += (ab|cd) d_i[bd].
// The matrices k_i and d_i are provided in *packed triangular form*: you supply
// a (nAoTr,nExchSets) matrix for each, where nAoTr is the packed triangular size
// (should be nAo*(nAo+1)/2) and nExchSets is the number of exchange operators to
// make (and associated number of input densities)
//
// NOTE:
// - this and AccCoul effectively do the very same thing just with the middle labels
//   exchanged. And a different interface.
// - the logical operation performed by this is similar to what we'd need for kext.
//   However, there (a) we'd certainly want to use real linear algebra routines to
//   execute it (optimally TBLIS or similar) and (b) we need to handle anti-symmetry,
//   too.
//   While neither aspect is there, I think the function is interesting from a
//   conceptual perspective (combined with its call pattern in AccJk4i, I mean)
void AccExch(FMatrixView ExchTr, FMatrixView const DensityTr, size_t nAo, double const *IR_RP pIntData,
   size_t StrideA, size_t StrideB, size_t StrideC, size_t StrideD, size_t iShA, size_t iShB, size_t iShC, size_t iShD,
   FRawBasis const *pOrbBasis, double fFactor, FMemoryStack &Mem)
{
   AssertSameShape1(ExchTr, DensityTr);
   assert(ExchTr.nRows == nAo*(nAo+1)/2);
   size_t
      iFnA = pOrbBasis->iFn(iShA), nFnA = pOrbBasis->nFn(iShA),
      iFnB = pOrbBasis->iFn(iShB), nFnB = pOrbBasis->nFn(iShB),
      iFnC = pOrbBasis->iFn(iShC), nFnC = pOrbBasis->nFn(iShC),
      iFnD = pOrbBasis->iFn(iShD), nFnD = pOrbBasis->nFn(iShD);
   size_t
      nSets = ExchTr.nCols;
   for (size_t iD = 0; iD < nFnD; ++ iD) {
      for (size_t iB = 0; iB < nFnB; ++ iB) {
         size_t
            iTrBD = SymIdx(iFnB + iB, iFnD + iD);
         TMemoryLock<double>
            pAccData(nSets, &Mem);
         for (size_t iSet = 0; iSet != nSets; ++ iSet)
            pAccData[iSet] = 0.;
         for (size_t iC = 0; iC < nFnC; ++ iC) {
            for (size_t iA = 0; iA < nFnA; ++ iA) {
               size_t
                  iTrAC = SymIdx(iFnA + iA, iFnC + iC);
               double
                  fThisInt = fFactor * pIntData[StrideA * iA + StrideB * iB + StrideC * iC + StrideD * iD];
               for (size_t iSet = 0; iSet != nSets; ++ iSet)
                  pAccData[iSet] += fThisInt * DensityTr(iTrAC, iSet);
            }
         }
         for (size_t iSet = 0; iSet != nSets; ++ iSet)
            ExchTr(iTrBD, iSet) += pAccData[iSet];
      }
   }
   IR_SUPPRESS_UNUSED_WARNING(nAo);
}


// this INCREMENTS Coul, ExchC, and ExchO
void AccJk4i(FLog &Log, double &fEnergyCoul, double &fEnergyExch, FMatrixView Coul, FMatrixView ExchC, FMatrixView ExchO,
   ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis,
   FMatrixView OrbC, FMatrixView OrbO, double fExchFactor, FMemoryStack &Mem_)
{
   RESET_CLOCKS

   RESUME_CLOCK(200)
   size_t
      nAo = pOrbBasis->nFn();
   if (nAo != OrbC.nRows || nAo != OrbO.nRows)
      throw std::runtime_error("AccJk4i: Input orbitals not consistent with orbital basis set.");

   FStackMatrix
      Density(nAo, nAo, &Mem_),  // total density: OrbC x OrbC.T + OrbO x OrbO.T
      DensityO(nAo, nAo, &Mem_); // ^- (note: orbs have absorbed occupation numbers)
   SyrkNT(Density, OrbC);
   SyrkNT(DensityO, OrbO);
   Add(Density, DensityO);

   // make packed triangular versions of input density matrices and output
   // coulomb/exchange matrices
   size_t
      nAoTr = nAo*(nAo+1)/2,
      nExchSets = 1;
   if (OrbO.nCols != 0)
      nExchSets = 2;
   FStackMatrix
      CoulTr(nAoTr, 1, &Mem_),
      ExchTr(nAoTr, nExchSets, &Mem_), // closed and open
      DensityTr(nAoTr, nExchSets, &Mem_); // closed and open
   Density.TriangularReduce(1, &DensityTr(0,0), 1.0);
   if (nExchSets >= 2)
      DensityO.TriangularReduce(1, &DensityTr(0,1), 1.0);

   size_t
      nAoShells = pOrbBasis->Shells.size();
   // allocate thread-local accumulation buffers and clear output
   FOmpAccBlock
      OmpAcc_CoulTr(CoulTr.pData, CoulTr.GetStridedSize(), OMPACC_ClearTarget, Mem_),
      OmpAcc_ExchTr(ExchTr.pData, ExchTr.GetStridedSize(), OMPACC_ClearTarget, Mem_);

   // should probably make something to enumerate shell pairs...
   // (or, better, actual integrals coming from an abstract source, like in ct8k...)
   FMemoryStackArray MemStacks(Mem_);
   #pragma omp parallel for schedule(dynamic)
   for (int iShAB_omp = 0; iShAB_omp < int(nAoShells*(nAoShells+1)/2); ++ iShAB_omp)
   {
      FMemoryStack
         &Mem = MemStacks.GetStackOfThread();
      FMatrixView
         CoulTr_Local(OmpAcc_CoulTr.pTls(), CoulTr.nRows, CoulTr.nCols, CoulTr.nRowSt, CoulTr.nColSt),
         ExchTr_Local(OmpAcc_ExchTr.pTls(), ExchTr.nRows, ExchTr.nCols, ExchTr.nRowSt, ExchTr.nColSt);
      size_t
         iShAB = size_t(iShAB_omp),
         iShA, iShB;
      UnpackIndex(iShB, iShA, iShAB);
      {
         ir::FRawShell const
            &ShA = pOrbBasis->Shells[iShA],
            &ShB = pOrbBasis->Shells[iShB];
         size_t
            nFnA = ShA.nFn(),
            nFnB = ShB.nFn();
         double
            fDistSqAB = DistSq(FVector3(ShA.vCen), FVector3(ShB.vCen));
         if (ShA.pRange && ShB.pRange && sqr(ShA.MaxCoRange() + ShB.MaxCoRange()) < fDistSqAB)
            continue;

         for (size_t iShC = 0; iShC < nAoShells; ++ iShC) {
            for (size_t iShD = 0; iShD <= iShC; ++ iShD) {
               size_t
                  iShCD = SymIdx(iShC, iShD);
               if (iShAB < iShCD)
                  continue;
               ir::FRawShell const
                  &ShC = pOrbBasis->Shells[iShC],
                  &ShD = pOrbBasis->Shells[iShD];
               size_t
                  nFnC = ShC.nFn(),
                  nFnD = ShD.nFn();
               double
                  fDistSqCD = DistSq(FVector3(ShC.vCen), FVector3(ShD.vCen));
               if (ShC.pRange && ShD.pRange && sqr(ShC.MaxCoRange() + ShD.MaxCoRange()) < fDistSqCD)
                  continue;

               // compute the unique integral block.
               size_t
                  Strides[4] = {1, nFnA, nFnA * nFnB, nFnA * nFnB * nFnC};
               TMemoryLock<double>
                  pIntData(nFnA * nFnB * nFnC * nFnD, &Mem);
               RESUME_CLOCK(100)
               ir::EvalInt2e4c(&pIntData[0], Strides, &ShA, &ShB, &ShC, &ShD, 1.0, pIntKernel, Mem);
               PAUSE_CLOCK(100)

               // accumulate Coulomb matrix
               RESUME_CLOCK(201)
               double fFactorJ = 1.;
               if (iShA != iShB)
                  fFactorJ *= 2.0;
               if (iShC != iShD)
                  fFactorJ *= 2.0;
               AccCoul(CoulTr_Local.pData, &DensityTr(0,0), nAo, pIntData, Strides[0], Strides[1], Strides[2], Strides[3], iShA, iShB, iShC, iShD, pOrbBasis, fFactorJ, Mem);
               if (iShAB != iShCD)
                  AccCoul(CoulTr_Local.pData, &DensityTr(0,0), nAo, pIntData, Strides[2], Strides[3], Strides[0], Strides[1], iShC, iShD, iShA, iShB, pOrbBasis, fFactorJ, Mem);
               PAUSE_CLOCK(201)

               // accumulate Exchange matrices
               if (fExchFactor != 0.) {
                  RESUME_CLOCK(202)
                  double fFactorK = 1.0;
                  if (iShAB != iShCD)
                     fFactorK *= 2.0;
                  AccExch(ExchTr_Local, DensityTr, nAo, pIntData, Strides[0], Strides[1], Strides[2], Strides[3], iShA, iShB, iShC, iShD, pOrbBasis, fFactorK, Mem);
                  if (iShA != iShB)
                     AccExch(ExchTr_Local, DensityTr, nAo, pIntData, Strides[1], Strides[0], Strides[2], Strides[3], iShB, iShA, iShC, iShD, pOrbBasis, fFactorK, Mem);
                  if (iShC != iShD)
                     AccExch(ExchTr_Local, DensityTr, nAo, pIntData, Strides[0], Strides[1], Strides[3], Strides[2], iShA, iShB, iShD, iShC, pOrbBasis, fFactorK, Mem);
                  if (iShA != iShB && iShC != iShD)
                     AccExch(ExchTr_Local, DensityTr, nAo, pIntData, Strides[1], Strides[0], Strides[3], Strides[2], iShB, iShA, iShD, iShC, pOrbBasis, fFactorK, Mem);
                  // ^- note: the routines AccCoul and AccExch essentially do the same thing
                  // apart from having the middle two shells and strides exchanged (and
                  // the exchange variant supporting multiple sets in parallel. Might
                  // want to use this fact if this is ever implemented properly...
                  PAUSE_CLOCK(202)
               }
            }
         }
      }
   }

   // collect intermediate data from threads and add results to output objects
   OmpAcc_ExchTr.Join();
   OmpAcc_CoulTr.Join();

   Add_TriangularExpanded(Coul, +1, &CoulTr(0,0), 1.0, 0.5); // sign, pDataTr, factor, off-diag-factor
   Add_TriangularExpanded(ExchC, +1, &ExchTr(0,0), fExchFactor, 0.5);
   if (nExchSets >= 2)
      Add_TriangularExpanded(ExchO, +1, &ExchTr(0,1), fExchFactor, 0.5);
//    // ^-- actually... if I do this based on density matrices, like here, do I need
//    // to get rid of some part of ExchC?

//    xout << fmt::format("fExchFactor : {:12.6f}", fExchFactor) << std::endl;
   fEnergyCoul = 0.5 * Dot(Density, Coul);
//    Add(ExchC, ExchO);
   // ^-- actually... if I do this based on density matrices, like here, that is not needed. Is it?
   fEnergyExch = -0.25 * Dot(Density, ExchC);
   if (OrbO.nCols != 0) {
      fEnergyExch += -0.25 * Dot(DensityO, ExchO);
   }
   PAUSE_CLOCK(200)

   IR_SUPPRESS_UNUSED_WARNING(Log);
}


FFockComponentBuilder4ixJk::FFockComponentBuilder4ixJk(FJkOptions const &Options_, FLog &Log_, FTimerSet *pTimers_)
   : FFockComponentBuilder(Log_, pTimers_), Options(Options_), EnergyCoulomb(0), EnergyExch(0), EnergyXc(0)
{
}


void FFockComponentBuilder4ixJk::Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem)
{
   pAtoms = &Atoms_;
   pOrbBasis = pOrbBasis_;
   m_WfDecl = WfDecl_;
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   IR_SUPPRESS_UNUSED_WARNING(Options_);
}


void FFockComponentBuilder4ixJk::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem)
{
   FFockComponentBuilder::AccFock(FockC, FockO, pBasis, COccC, COccO, Flags, Mem); // consistency checks.

   size_t
      nBf = pBasis->nFn();
   FStackMatrix
      Coul(nBf, nBf, &Mem),
      ExchC(nBf, nBf, &Mem),
      ExchO(nBf, nBf, &Mem);
   Coul.Clear();
   ExchC.Clear();
   ExchO.Clear(); // <- otherwise error in exch energy accumulation...
   AccJk4i(m_Log, EnergyCoulomb, EnergyExch, Coul, ExchC, ExchO, &ir::g_IrCoulombKernel, &*pBasis->pRawBasis,
      COccC, COccO, Options.fExchFactor, Mem);
//    Scale(FockO, -.5);
   Add(FockC, Coul, 1.0);
   Add(FockO, ExchO, -.5);
   // note: the default factor is -.5
   Add(FockC, ExchC, -.5);
//    Scale(FockO, -.5);
//    EnergyExch *= Options.fExchFactor;
   if (Options.UseBandStructureEnergyFormula) {
      // flip sign of regular coulomb and exchange energy to counter the double-
      // counted 2e terms in the band structure term.
      EnergyExch *= -1.;
      EnergyCoulomb *= -1.;
      // ^- not sure about open-shell stuff.
   }
   Energy = EnergyExch + EnergyCoulomb + EnergyXc;
}


void FFockComponentBuilder4ixJk::PrintEnergyContribs()
{
   m_Log.WriteResult("Coulomb energy", EnergyCoulomb);
   m_Log.WriteResult("Exchange energy", EnergyExch);
   if (EnergyXc != 0)
      m_Log.WriteResult("Density functional energy", EnergyXc);
}


FFockComponentBuilder4ixJk::~FFockComponentBuilder4ixJk()
{
}





// forms the integral matrix set (\mu\nu|F) for a single set of F shells.
// \mu comes from OrbBasis1, \nu comes from OrbBasis2. If both pointers are identical,
// symmetric integrals will be assumed [(\mu\nu|F) = (\nu\mu|F)]  and this symmetry will be used.
double *FormIntMNF(ir::FRawShell const &ShF,
   ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasisA, FRawBasis const *pOrbBasisB, FRawBasis const *pFitBasis, FMemoryStack &Mem, FMatrixView ScrDen, double fThr, unsigned Flags)
{
   size_t
      nAoA = pOrbBasisA->nFn(),
      nAoB = pOrbBasisB->nFn(),
      nFnF = ShF.nFn();
   double
      *pMNF;
   bool
      Symmetric = (pOrbBasisA == pOrbBasisB);
   Mem.Alloc(pMNF, nAoA * nAoB * nFnF);

   size_t nIntTotal = 0, nIntRetained = 0;

   size_t const
      StrideA = 1, StrideB = nAoA, StrideF = nAoA*nAoB;
//       StrideA = nAoB, StrideB = 1, StrideF = nAoA*nAoB;
   size_t
      Strides[3] = {StrideA, StrideB, StrideF};

   for ( size_t iShB = 0; iShB < pOrbBasisB->Shells.size(); ++ iShB ) {
      ir::FRawShell const &ShB = pOrbBasisB->Shells[iShB];
      size_t nFnB = ShB.nFn();
//       size_t iShA_First = Symmetric ? iShB : 0;
//       for ( size_t iShA = iShA_First; iShA < pOrbBasisA->Shells.size(); ++ iShA ) {
      // ^- this variant is faster, but i found no sensible way of combining it with the format used in TriangularReduce.
      size_t iShA_Last = Symmetric ? iShB + 1 : pOrbBasisA->Shells.size();
      for ( size_t iShA = 0; iShA < iShA_Last; ++ iShA ) {
         ir::FRawShell const &ShA = pOrbBasisA->Shells[iShA];
         size_t nFnA = ShA.nFn();
         nIntTotal += nFnA * nFnB;

//          if (ScrDen(iShA, iShB) < fThr)
//             continue;
         double
            fDistSqAB = DistSq(FVector3(ShA.vCen), FVector3(ShB.vCen));
         double
            *pIntData = &pMNF[StrideA*pOrbBasisA->iFn(iShA) + StrideB*pOrbBasisB->iFn(iShB)];
         if (ShA.pRange && ShB.pRange && sqr(ShA.MaxCoRange() + ShB.MaxCoRange()) < fDistSqAB) {
            // screened out.
            for ( size_t iF = 0; iF < nFnF; ++ iF )
               for ( size_t iB = 0; iB < nFnB; ++ iB )
                  for ( size_t iA = 0; iA < nFnA; ++ iA )
                     pIntData[iA*StrideA + iB*StrideB + iF*StrideF] = 0.;
         } else {
            nIntRetained += nFnA * nFnB;
            ir::EvalInt2e3c(pIntData, Strides, &ShA, &ShB, &ShF,1, 1.0, pIntKernel, Mem);
         }
      }
   }

   // fix up upper triangle.
   if (Symmetric && ((Flags * INTFLAGS_SkipTriangleFix)==0)) {
      if (1) {
         for (size_t iF = 0; iF < nFnF; ++ iF)
            for (size_t iB = 0; iB < nAoB; ++ iB)
               for (size_t iA = iB; iA < nAoA; ++ iA)
                  pMNF[iA * StrideA + iB * StrideB + iF * StrideF] = pMNF[iB * StrideA + iA * StrideB + iF * StrideF];
      } else {
         for (size_t iF = 0; iF < nFnF; ++ iF)
            for (size_t iB = 0; iB < nAoB; ++ iB)
               for (size_t iA = 0; iA <= iB; ++ iA)
                  pMNF[iA * StrideA + iB * StrideB + iF * StrideF] = pMNF[iB * StrideA + iA * StrideB + iF * StrideF];
      }
   }

   return pMNF;
   IR_SUPPRESS_UNUSED_WARNING(pFitBasis);
   IR_SUPPRESS_UNUSED_WARNING(ScrDen);
   IR_SUPPRESS_UNUSED_WARNING(fThr);
}




double *FormIntMNF(ir::FRawShell const &ShF,
   ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis, FRawBasis const *pFitBasis, FMemoryStack &Mem, FMatrixView ScrDen, double fThr)
{
   return FormIntMNF(ShF, pIntKernel, pOrbBasis, pOrbBasis, pFitBasis, Mem, ScrDen, fThr);
}





// contract 3-index gradient with factoriable 3-index density matrix.
//
//    grad[q] += (µν|A)^q D[μν] v[A]
//
// assumes D[μν] == D[νμ] if pBasisM == pBasisN.
void AccGradient3ix_Factored(double *pGrad,
   double const *pD_MN, size_t StrideM, size_t StrideN,
   double const *pV_F,
   FRawBasis const *pBasisM, size_t iShM_, size_t nShM,
   FRawBasis const *pBasisN, size_t iShN_, size_t nShN,
   FRawBasis const *pBasisF, size_t iShF_, size_t nShF,
   ir::FIntegralKernel *pIntKernel,
   FMemoryStack &Mem_, double Prefactor_ = 1.)
{
   bool Symmetric = pBasisM == pBasisN;

   FMemoryStackArray MemStacks(Mem_);
   #pragma omp parallel for schedule(dynamic)
   for (int iShM__ = iShN_; iShM__ < int(iShN_ + nShM); ++ iShM__) {
      size_t iShM = size_t(iShM__);
      FMemoryStack &Mem = MemStacks.GetStackOfThread();

      ir::FRawShell const &ShM = pBasisM->Shells[iShM];
      size_t nFnM = ShM.nFn();

      for (size_t iShN = iShN_; iShN < iShN_ + nShN; ++ iShN) {
         ir::FRawShell const &ShN = pBasisN->Shells[iShN];
         size_t nFnN = ShN.nFn();

         if (Symmetric && iShN < iShM)
            continue;
         double
            Prefactor = Prefactor_;
         if (Symmetric && iShM != iShN)
            Prefactor *= 2.;

         // iterate through F shells. We can do all shells sitting on the same
         // center in one go ...but not more. Need C's center id for the
         // gradient.
         size_t iShF_Next;
         for (size_t iShF = iShF_; iShF != iShF_ + nShF; iShF = iShF_Next) {
            size_t
               iCenM = pBasisM->iShCen(iShM),
               iCenN = pBasisN->iShCen(iShN),
               iCenF = pBasisF->iShCen(iShF);

            // find first shell on next center.
            iShF_Next = iShF + 1;
            while (iShF_Next != iShF_ + nShF && iCenF == pBasisF->iShCen(iShF_Next))
               iShF_Next += 1;

            size_t
               IntStrides[4] = {1, nFnM, 1, nFnM*nFnN}; // note: only one C vector. Therefore c-stride 1.
            FStackMatrix
               IntDeriv(nFnM*nFnN, 6, &Mem);
//             IntDeriv.Clear(); // FIXME: remove this.
            EvalInt2e3c1d_ContractC(IntDeriv.pData, IntStrides,
               &pV_F[pBasisF->iFn(iShF) - pBasisF->iFn(iShF_)], 1, 1, // 1,1: stride (between vectors), nVecC.
               &ShM, &ShN, &pBasisF->Shells[iShF], iShF_Next - iShF,
               Prefactor, pIntKernel, Mem);

            double const
               *pD = &pD_MN[StrideM * (pBasisM->iFn(iShM) - pBasisM->iFn(iShM_)) +
                            StrideN * (pBasisN->iFn(iShN) - pBasisN->iFn(iShN_))];

            // fixme: accumulate stuff in per-thread gradient. can this even be done with OpenMP?
            #pragma omp critical
            {
               for (size_t iM = 0; iM < nFnM; ++ iM)
                  for (size_t iN = 0; iN < nFnN; ++ iN) {
                     double
                        Rdm = pD[iM * StrideM + iN * StrideN];
                     for (size_t ixyz = 0; ixyz < 3; ++ ixyz) {
                        double
                           dM = Rdm * IntDeriv(iM + nFnM * iN, (ixyz  )),
                           dN = Rdm * IntDeriv(iM + nFnM * iN, (ixyz+3)),
                           dF = -dM - dN;
                        pGrad[ixyz + 3*iCenM] += dM;
                        pGrad[ixyz + 3*iCenN] += dN;
                        pGrad[ixyz + 3*iCenF] += dF;
                     }
                  }
            }
         }
      }
   }
}


// contract 3-index gradient with non-factoriable 3-index density matrix.
//
//    grad[q] += (µν|A)^q D[μνA]
//
// assumes D[μνA] == D[νμA] if pBasisM == pBasisN.
void AccGradient3ix_NonFactored(double *pGrad,
   double const *pD_MNF, size_t StrideM, size_t StrideN, size_t StrideF,
   FRawBasis const *pBasisM, size_t iShM_, size_t nShM,
   FRawBasis const *pBasisN, size_t iShN_, size_t nShN,
   FRawBasis const *pBasisF, size_t iShF_, size_t nShF,
   ir::FIntegralKernel *pIntKernel,
   FMemoryStack &Mem, double Prefactor_ = 1.)
{
   bool Symmetric = pBasisM == pBasisN;

   for (size_t iShF = iShF_; iShF < iShF_ + nShF; ++ iShF) {
      ir::FRawShell const &ShF = pBasisF->Shells[iShF];
      size_t nFnF = ShF.nFn();

      for (size_t iShM = iShN_; iShM < iShN_ + nShM; ++ iShM) {
         ir::FRawShell const &ShM = pBasisM->Shells[iShM];
         size_t nFnM = ShM.nFn();

         for (size_t iShN = iShN_; iShN < iShN_ + nShN; ++ iShN) {
            ir::FRawShell const &ShN = pBasisN->Shells[iShN];
            size_t nFnN = ShN.nFn();

            if (Symmetric && iShN < iShM)
               continue;
            double
               Prefactor = Prefactor_;
            if (Symmetric && iShM != iShN)
               Prefactor *= 2.;

            size_t
               IntStrides[4] = {1, nFnM, nFnM*nFnN, nFnM*nFnN*nFnF};
            FStackMatrix
               IntDeriv(nFnM*nFnN, nFnF*6, &Mem);
            EvalInt2e3c1d(IntDeriv.pData, IntStrides, &ShM, &ShN, &ShF, 1,
               Prefactor, pIntKernel, Mem);

            size_t
               iCenM = pBasisM->iShCen(iShM),
               iCenN = pBasisN->iShCen(iShN),
               iCenF = pBasisF->iShCen(iShF);
            double const
               *pD = &pD_MNF[StrideM * (pBasisM->iFn(iShM) - pBasisM->iFn(iShM_)) +
                             StrideN * (pBasisN->iFn(iShN) - pBasisN->iFn(iShN_)) +
                             StrideF * (pBasisF->iFn(iShF) - pBasisF->iFn(iShF_))];

            for (size_t iF = 0; iF < nFnF; ++ iF)
               for (size_t iM = 0; iM < nFnM; ++ iM)
                  for (size_t iN = 0; iN < nFnN; ++ iN) {
                     double
                        Rdm = pD[iF * StrideF + iM * StrideM + iN * StrideN];
                     for (size_t ixyz = 0; ixyz < 3; ++ ixyz) {
                        double
                           dM = Rdm * IntDeriv(iM + nFnM * iN, iF + (ixyz  )*nFnF),
                           dN = Rdm * IntDeriv(iM + nFnM * iN, iF + (ixyz+3)*nFnF),
                           dF = -dM - dN;
                        pGrad[ixyz + 3*iCenM] += dM;
                        pGrad[ixyz + 3*iCenN] += dN;
                        pGrad[ixyz + 3*iCenF] += dF;
                     }
                  }
         }
      }
   }
}



// forms the contraction Out[\mu,\nu] := \sum_F (\mu\nu|F) Gamma[F]
// Note: this INCREMENTS the output matrix!
// (FIXME: but it appears as if it is not doing it quite correctly;
// this generates incorrect results unless I clear the target and add...
// this is why at this moment the CoulC matrices exist in FFockComponentBuilderDfJk and FFockComponentBuilderDfCoul...)
void FormIntMNF_ContractF(FMatrixView Out, ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis,
   FRawBasis const *pFitBasis, double const *pGamma, FMemoryStack &Mem_)
{
   size_t
      nAo = pOrbBasis->nFn(),
      nFit = pFitBasis->nFn();
   assert(Out.nRows == nAo && Out.nCols == nAo);

   FMemoryStackArray MemStacks(Mem_);
   #pragma omp parallel for schedule(dynamic)
   for (int iShB__ = 0; size_t(iShB__) < pOrbBasis->Shells.size(); ++ iShB__ ){
      size_t iShB = size_t(iShB__);
      FMemoryStack &Mem = MemStacks.GetStackOfThread();

      ir::FRawShell const &ShB = pOrbBasis->Shells[iShB];
      size_t nFnB = ShB.nFn();
      for (size_t iShA = iShB; iShA < pOrbBasis->Shells.size(); ++ iShA) {
         ir::FRawShell const &ShA = pOrbBasis->Shells[iShA];
         double
            fDistSqAB = DistSq(FVector3(ShA.vCen), FVector3(ShB.vCen));
         if (ShA.pRange && ShB.pRange && sqr(ShA.MaxCoRange() + ShB.MaxCoRange()) < fDistSqAB)
            continue;
         size_t
            iFnA = pOrbBasis->iFn(iShA),
            iFnB = pOrbBasis->iFn(iShB),
            nFnA = ShA.nFn(),
            Strides[3] = {1, nAo, nAo * nAo};

         ir::EvalInt2e3c_ContractC(&Out(iFnA,iFnB),Strides, pGamma,nFit,1, &ShA, &ShB,
            &pFitBasis->Shells[0],pFitBasis->Shells.size(), 1.0, pIntKernel, Mem);

         for (size_t iB = iFnB; iB < iFnB + nFnB; ++ iB)
            for (size_t iA = iFnA; iA < iFnA + nFnA; ++ iA)
               // assign to (\mu\nu| and (\nu\mu|. (int has perm symmetry).
               Out(iB,iA) = Out(iA,iB);
      }
   }
}


void Solve3ixFittingEq(double *pDF_NFi, size_t nAo, size_t nFit, size_t nOcc, FMatrixView Jcd)
{
   for (size_t i = 0; i < nOcc; ++ i) {
      FMatrixView
         D_NF(pDF_NFi + nAo * nFit * i, nAo, nFit);
      TriangularSolve(Transpose(D_NF), Jcd);
   }
}


// NOTE: OrbC and OrbO are supposed to have their respective occupation numbers absorbed,
// such that the closed/open densities are obtained by DenC := OrbC x OrbC.T and DenO := OrbO x OrbO.T
// FIXME: rephrase that from TN matrix multiplication to NT matrix multiplication(!!)
void MakeJk(double &fEnergyCoul, double &fEnergyExch, FMatrixView Coul, FMatrixView ExchC, FMatrixView ExchO,
   ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis, FRawBasis const *pFitBasis,
   FMatrixView OrbC, FMatrixView OrbO, FMatrixView Jcd, FMemoryStack &Mem_)
{
   RESET_CLOCKS

   assert_rt(ExchO.pData != 0 || OrbO.nCols == 0);
   size_t
      nAo = pOrbBasis->nFn(),
      nFit = pFitBasis->nFn();
   if (nAo != OrbC.nRows || nAo != OrbO.nRows || nFit != Jcd.nRows)
      throw std::runtime_error("MakeJk: Input orbitals not consistent with orbital basis set.");

   // first, make the exchange matrices and j[A] = (A|rs) gamma[r,s].
//    FStackMatrix
//       k1c_NFi(nAo, nFit * OrbC.nCols, &Mem_), // (mu|A i)
//       k1o_NFi(nAo, nFit * OrbO.nCols, &Mem_);
   FHeapMatrix
//       // ^-- UPDATE: 2021-05-24 considered changing to heap matrices. Note that for the
//       // cached 3ix integrals, the fitting intermediates are also placed on the
//       // heap, not the stack. We should probably keep some consistency in this
//       // regard. On the other hand, this here is really no permanent data, and if
//       // we keep it on stack, we can just assign everything to the stack elements.
//       // so... not sure how best to proceed.
      k1c_NFi(nAo, nFit * OrbC.nCols), // (mu|A i)
      k1o_NFi(nAo, nFit * OrbO.nCols);
   FStackMatrix
      Density(nAo, nAo, &Mem_), // OrbC x OrbC.T + OrbO x OrbO.T
      jgamma(nFit, 1, &Mem_); // (A|rs) gamma[r,s].
   SyrkNT(Density, OrbC);
   SyrkNT(Density, OrbO, 1., MXM_Add);
   jgamma.Clear();

   FMatrixView ScrDen;

   {
      FMemoryStackArray MemStacks(Mem_);
      #pragma omp parallel for schedule(dynamic)
      for ( int iShF__ = 0; iShF__ < int(pFitBasis->Shells.size()); ++ iShF__ ){
         size_t iShF = size_t(iShF__);
         FMemoryStack &Mem = MemStacks.GetStackOfThread();

         ir::FRawShell const &ShF = pFitBasis->Shells[iShF];
         size_t nFnF = ShF.nFn();

         double
            // (\mu\nu|F): nAo x nAo x nFnF
            *pMNF = FormIntMNF(ShF, pIntKernel, pOrbBasis, pFitBasis, Mem, ScrDen, 0.);

         // D[\nu A i] = (i\nu|A) = C(\mu,i) (\mu [\nu|A)]
         RESUME_CLOCK(200)
         FMatrixView
            DF_NFIc(k1c_NFi.pData + nAo*pFitBasis->iFn(iShF), nAo * nFnF, OrbC.nCols, 1, nFit * nAo),
            DF_NFIo(k1o_NFi.pData + nAo*pFitBasis->iFn(iShF), nAo * nFnF, OrbO.nCols, 1, nFit * nAo);
   //    FTimer TimerMxm;
         Mxm(DF_NFIc,
            FMatrixView(pMNF, nAo * nFnF, nAo, nAo, 1),
            OrbC);
   //       xout << fmt::format("mxm: {} x {} x {}:  {:12.3f} GFLOP/s\n", (nAo*nFnF), nAo, OrbC.nCols, (2e-6*(nAo*nFnF) * nAo * OrbC.nCols / (double)TimerMxm));

         Mxm(DF_NFIo,
            FMatrixView(pMNF, nAo * nFnF, nAo, nAo, 1),
            OrbO);
         PAUSE_CLOCK(200)
         // jgamma[A] = (\mu\nu|A) Density[\mu,\nu]
         Mxv(&jgamma(pFitBasis->iFn(iShF),0), 1,
            pMNF, nAo*nAo, 1,
            &Density(0,0), 1, nFnF, nAo*nAo);

         Mem.Free(pMNF);
      }
   }

   //   Solve[Jcd[AB] D[\nu B i]] -> D[\nu A I]
   RESUME_CLOCK(201)
   CholeskySolve(jgamma, Jcd);
   Solve3ixFittingEq(&k1c_NFi(0,0), nAo, nFit, OrbC.nCols, Jcd);
   Solve3ixFittingEq(&k1o_NFi(0,0), nAo, nFit, OrbO.nCols, Jcd);
   PAUSE_CLOCK(201)

   // assemble exchange matrices.
   RESUME_CLOCK(202)
   SyrkNT(ExchC, k1c_NFi);
   if (ExchO.pData != 0)
      SyrkNT(ExchO, k1o_NFi);
   PAUSE_CLOCK(202)

   // recalculate integrals to form the coulomb matrix.
   RESUME_CLOCK(300)

//    fEnergyCoul = 0.;
   if (false) {
      Coul.Clear();
      for (size_t iShF = 0; iShF != pFitBasis->Shells.size(); ++ iShF){
         ir::FRawShell const &ShF = pFitBasis->Shells[iShF];
         size_t nFnF = ShF.nFn();
         double
            // (\mu\nu|F): nAo x nAo x nFnF
            *pMNF = FormIntMNF(ShF, pIntKernel, pOrbBasis, pFitBasis, Mem_, ScrDen, 0.);

         // j[\mu,\nu] = (\mu\nu,A) jgamma[A]
         Mxv(&Coul(0,0), 1,
            pMNF, 1, nAo*nAo,
            &jgamma(pFitBasis->iFn(iShF),0), 1,
            nAo*nAo, nFnF, true, 1.0);

         Mem_.Free(pMNF);
      }
   } else {
      Coul.Clear();
//       fEnergyCoul -= 0.5 * Dot(Density, Coul);
      FormIntMNF_ContractF(Coul, pIntKernel, pOrbBasis, pFitBasis, jgamma.pData, Mem_);
   }
   PAUSE_CLOCK(300)

//    fEnergyCoul += 0.5 * Dot(Density, Coul);
   fEnergyCoul = 0.5 * Dot(Density, Coul);
   if (OrbO.nCols != 0)
      Add(ExchC, ExchO);
   fEnergyExch = -0.25 * Dot(Density, ExchC);
   if (OrbO.nCols != 0) {
      FStackMatrix
         DensityO(nAo, nAo, &Mem_);
      SyrkNT(DensityO, OrbO, 1.);
      fEnergyExch += -0.25 * Dot(DensityO, ExchO);
   }
}

// form the transformed integrals (ij|A).
// iSolveFitEq == 0: Do nothing with A.
// iSolveFitEq == 1: Solve half-sided fitting equations (L side)
// iSolveFitEq == 2: Solve full L L^T fitting equations.
FMatrixView MakeF_ijF(ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis, FRawBasis const *pFitBasis,
   FMatrixView OrbI, FMatrixView OrbJ, FMatrixView Jcd, int iSolveFitEq, FMemoryStack &Mem)
{
   size_t
      nOrbI = OrbI.nCols,
      nOrbJ = OrbJ.nCols,
      nAo = pOrbBasis->nFn(),
      nFit = pFitBasis->nFn();
   if (nAo != OrbI.nRows || nAo != OrbJ.nRows || nFit != Jcd.nRows)
      throw std::runtime_error("MakeJk: Input orbitals not consistent with orbital basis set.");

   FMatrixView
      F_Aij = MakeStackMatrix(nOrbI * nOrbJ, nFit, Mem);

   FMatrixView ScrDen;
   for (size_t iShF = 0; iShF != pFitBasis->Shells.size(); ++ iShF) {
      ir::FRawShell const &ShF = pFitBasis->Shells[iShF];
      size_t nFnF = ShF.nFn();

      double
         // (\mu\nu|F): nAo x nAo x nFnF
         *pMNF = FormIntMNF(ShF, pIntKernel, pOrbBasis, pFitBasis, Mem, ScrDen, 0.);

      // D[ijA] = (ij|A) = C(\mu,i) C(\nu,j) (\mu \nu|A)
      for (uint iF = 0; iF < nFnF; ++ iF) {
         FMatrixView
            MN(&pMNF[nAo * nAo * iF], nAo, nAo),
            ij(&F_Aij(0, iF + pFitBasis->iFn(iShF)), nOrbI, nOrbJ);
         ChainMxm(ij, Transpose(OrbI), MN, OrbJ, Mem);
         // ^- could be refactored to include nFnF dimension in one of the MxMs.
         //    don't care about it at this moment, though.
      }

      Mem.Free(pMNF);
   }


   if (iSolveFitEq == 1)
      TriangularSolve(Transpose(F_Aij), Jcd);
   else if (iSolveFitEq == 2)
      CholeskySolve(Transpose(F_Aij), Jcd);
   return F_Aij;
}



FFockComponentBuilderDfJk::FFockComponentBuilderDfJk(FDfJkOptions const &Options_, FLog &Log_, FTimerSet *pTimers_)
   : FFockComponentBuilder(Log_, pTimers_), Options(Options_), EnergyCoulomb(0), EnergyExch(0), EnergyXc(0)
{
   if (Options.fExchFactor == 0.) {
      Log_.EmitWarning("FFockComponentBuilderDfJk instanciated with ExchFactor == 0. This is not efficient! Use one of the Coulomb algorithms instead!");
   }
}


FFockComponentBuilderDfJk::~FFockComponentBuilderDfJk()
{
}


void FFockComponentBuilderDfJk::Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem)
{
   pAtoms = &Atoms_;
   pOrbBasis = pOrbBasis_;
   m_WfDecl = WfDecl_;

   InitFittingBasisAndJcd(pFitBasis, Jcd, &Atoms_, BASIS_JkFit, &Options_.Print, &m_Log, Mem);
   nFit = pFitBasis->nFn();
}


void FFockComponentBuilderDfJk::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem)
{
   FFockComponentBuilder::AccFock(FockC, FockO, pBasis, COccC, COccO, Flags, Mem); // consistency checks.

   size_t
      nBf = pBasis->nFn();
   // TODO: add options to MakeJk to *update* the fock matrix instead of overwriting it?
   // TODO: maybe do these with some less clearing and copying and adding of matrices?
   // Note:
   // - ...it is not so trivial, however, because there may be previous results from
   //   other Fock component builders, in particular, from FFockComponentBuilderXc.
   // - At the very least we'd have to accept some additional numerical rounding errors
   //   due to collecting first EnergyCoulomb (to get the baseline), and then re-computing it.
   {
      TMemoryLock<double>
         pFreeMe(0, &Mem);
      FMatrixView
         CoulC = MakeStackMatrix(nBf, nBf, Mem),
         ExchC = MakeStackMatrix(nBf, nBf, Mem),
         ExchO;
      if (COccO.nCols != 0)
         ExchO = MakeStackMatrix(nBf, nBf, Mem);
      MakeJk(EnergyCoulomb, EnergyExch, CoulC, ExchC, ExchO,
         &ir::g_IrCoulombKernel, &*pBasis->pRawBasis, &*pFitBasis->pRawBasis,
         COccC, COccO, Jcd, Mem);
      Add(FockC, CoulC);
      // note: the default factor is -.5
      Add(FockC, ExchC, -.5*Options.fExchFactor);
      if (COccO.nCols != 0)
         Add(FockO, ExchO, -.5*Options.fExchFactor);
      EnergyExch *= Options.fExchFactor;
      if (Options.UseBandStructureEnergyFormula) {
         // flip sign of regular coulomb and exchange energy to counter the double-
         // counted 2e terms in the band structure term.
         EnergyExch *= -1.;
         EnergyCoulomb *= -1.;
         // ^- not sure about open-shell stuff.
      }
      Energy = EnergyExch + EnergyCoulomb + EnergyXc;
   }
}


void FFockComponentBuilderDfJk::PrintEnergyContribs()
{
   m_Log.WriteResult("Coulomb energy", EnergyCoulomb);
   m_Log.WriteResult("Exchange energy", EnergyExch);
   if (EnergyXc != 0)
      m_Log.WriteResult("Density functional energy", EnergyXc);
}


static void ScaleDijF_ActiveActive(FMatrixView DijF, double f, size_t nOccC, size_t nOccO)
{
   size_t
      nOcc = nOccC + nOccO;
   assert(DijF.nRows == nOcc*nOcc);
   for (size_t iFit = 0; iFit != DijF.nCols; ++ iFit)
      for (size_t j = nOccC; j != nOccC + nOccO; ++ j)
         for (size_t i = nOccC; i != nOccC + nOccO; ++ i)
            DijF(i + nOcc*j, iFit) *= f;
}


void FFockComponentBuilderDfJk::AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   // combine the two occupied orbital subsets into a single matrix
   // (note: they are already occupation number scaled).
   size_t
      nOcc = COccC.nCols + COccO.nCols,
      nAo = pOrbBasis->nFn();
   assert(COccO.nRows == COccC.nRows);
   assert(nAo == COccO.nRows);
   FStackMatrix
      COcc(nAo, nOcc, &Mem);
   Assign(Select(COcc,0,0,nAo,COccC.nCols), COccC);
   Assign(Select(COcc,0,COccC.nCols,nAo,COccO.nCols), COccO);

   // make the density matrix.
   // (note: this one could be saved by passing Rdm from the outside.)
   FStackMatrix
      Rdm(nAo, nAo, &Mem);
   SyrkNT(Rdm, COcc); // <- includes both closed and open shells.

   TMemoryLock<char>
      pFreeMe(0, &Mem);
   {
      // okay.. now the absorbed occupation numbers are in the closed-shell
      // orbitals. The correspondig energy formula is:
      //
      //   Ec + Ex = \sum_{ij}  .5 * (ii|jj) - .25 * (ij|ij)
      //
      // One of the .5's stems from the self-interaction of orbitals (technically
      // the sum should be over i, j >= i.
      // This gives us a total factor of 4 compared to my implementation in Molpro.

      FRawBasis const
         *pOrbBasisRaw = &*pOrbBasis->pRawBasis,
         *pFitBasisRaw = &*pFitBasis->pRawBasis;

      // 2e terms (see exchg.pdf).
      //
      // Coulomb part:
      //
      //   v[A] = J_AB^{-1} (μν|B) D^{μν}
      //   Ecoul^q = D^{μν} (µν|A)^q v[A] - (1/2) v_A (A|B)^q v_B
      //
      // (note: should use two-sided inline contracting driver. Or at least 1-sided
      //  inline contracting driver (on v[A] side))
      //
      // TODO: also, we could just store the density from the last Fock build;
      // would save one evaluation of integrals.
      // FFockComponentBuilderDfCoulXcCached does it like this.
      FMatrixView
         // (note: needed for exchange)
         D_ijA = MakeF_ijF(&ir::g_IrCoulombKernel, pOrbBasisRaw, pFitBasisRaw,
            COcc, COcc, Jcd, 2, Mem);

      // make v -- expansion coefficients of fitted coulomb potential.
      FStackMatrix
         v(nFit, 1, &Mem);
      for (size_t iFit = 0; iFit < nFit; ++ iFit) {
         v[iFit] = 0.;
         for (size_t iOcc = 0; iOcc < nOcc; ++ iOcc)
            v[iFit] += D_ijA(iOcc * (nOcc+1), iFit);
      }


      ScaleDijF_ActiveActive(D_ijA, 2., COccC.nCols, COccO.nCols);

      double
         ExchFac = Options.fExchFactor,
         CoulFac = 1.;
      // ! 1st contribution: 2 D[\mu,\nu,A] (\mu\nu|A)^q
      for (size_t iShF = 0; iShF < pFitBasisRaw->nSh(); ++ iShF) {
         size_t
            iFnF = pFitBasisRaw->iFn(iShF),
            nFnF = pFitBasisRaw->nFn(iShF);
         double
            DFac = -ExchFac * .5;
         // transform d[ijA] back to MO basis.
         FStackMatrix
            D_MNF(nAo*nAo, nFnF, &Mem);
         if (ExchFac != 0.) {
            // add exchange contribution:
            //   -0.5 C[\mu,i] C[\nu,i]
            for (size_t iF = 0; iF < nFnF; ++ iF) {
               FMatrixView
                  MN(&D_MNF(0,iF), nAo, nAo),
                  ij(&D_ijA(0,iF + iFnF), nOcc, nOcc);
               ChainMxm(MN, COcc, ij, Transpose(COcc), Mem);
            }
         } else {
            // you probably should not be here in case ExchFac = 0.0,
            // but let's make it get the right result anyway...
            D_MNF.Clear();
            DFac = 1.0;
            // (^- note: why the /DFac business? these guys below there have
            // inline argument scaling anyway, but ChainMxm doesn't, and the
            // coefficient object D[\mu\nu F] is large).
         }
         // add coulomb contribution:
         //   +1.0 D[\mu\nu] d[A] (\mu\nu|A)^q
         DGER(nAo*nAo, nFnF, CoulFac/DFac, Rdm.pData,1, &v[iFnF],1, D_MNF.pData,D_MNF.nColSt);

         // contract with derivative integrals.
         AccGradient3ix_NonFactored(Gradient.pData,
            D_MNF.pData, 1, nAo, nAo*nAo,
            pOrbBasisRaw, 0, pOrbBasisRaw->nSh(),
            pOrbBasisRaw, 0, pOrbBasisRaw->nSh(),
            pFitBasisRaw, iShF, 1,
            &ir::g_IrCoulombKernel, Mem, DFac);
      }

      // ! 2nd contribution: +D[ijA] J^q[AB] D[ijB]    (exchange)
      // !                   -.5 * D[A] J^q[AB] D[B]   (coulomb)
      {
         FStackMatrix
            DAB(nFit, nFit, &Mem);
         ScaleDijF_ActiveActive(D_ijA, std::sqrt(2.)/2., COccC.nCols, COccO.nCols);
         SyrkTN(DAB, D_ijA);
         if (ExchFac != 1.)
            Scale(DAB, ExchFac);
         // add coulomb contribution: -.5 * D[A] * D[B]  (note the .25 on the integrals)
         if (CoulFac != 0.)
            DGER(nFit, nFit, -2.*CoulFac, v.pData,1, v.pData,1, DAB.pData, nFit);
         AccGradient2ix(Gradient.pData, DAB, *pFitBasisRaw,
            *pFitBasisRaw, FKrn2i_Direct(&ir::g_IrCoulombKernel), Mem, .25);
      }

   }
}






FFockComponentBuilderDfCoul::FFockComponentBuilderDfCoul(FDfJkOptions const &JkOptions_, FLog &Log, FTimerSet *pTimers_)
   : FFockComponentBuilder(Log, pTimers_), JkOptions(JkOptions_)
{
}

FFockComponentBuilderDfCoul::~FFockComponentBuilderDfCoul()
{
}


void FFockComponentBuilderDfCoul::Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem)
{
   IR_SUPPRESS_UNUSED_WARNING(WfDecl_);
   pAtoms = &Atoms_;
   pOrbBasis = pOrbBasis_;

   // Compute cholesky decomposition of fitting metric J_AB
   InitFittingBasisAndJcd(pFitBasis, Jcd, &Atoms_, BASIS_JFit, &Options_.Print, &m_Log, Mem, "coul");
   nFit = pFitBasis->nFn();
}


void MakeCoul(FMatrixView Coul, double &EnergyCoul, ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis, FRawBasis const *pFitBasis,
   FMatrixView OrbC, FMatrixView OrbO, FMatrixView Jcd, FMemoryStack &Mem_, FTimerSet *m_pTimers)
{
//    Coul.Clear(); // FIXME: required?
   size_t
      nAo = pOrbBasis->nFn(),
      nFit = pFitBasis->nFn();
   if (nAo != OrbC.nRows || nAo != OrbO.nRows || nFit != Jcd.nRows)
      throw std::runtime_error("MakeCoul: Input orbitals not consistent with orbital basis set.");

   FStackMatrix
      Density(nAo, nAo, &Mem_), // OrbC x OrbC.T + OrbO x OrbO.T
      jgamma(nFit, 1, &Mem_);   // (A|rs) gamma[r,s].

   { TIME_SECTION(m_pTimers, 0x201, "DF-J (RDM)");
      SyrkNT(Density, OrbC);
      SyrkNT(Density, OrbO, 1.0, MXM_Add);
   }

   FMatrixView ScrDen;

   RESUME_CLOCK(300) // <- note: those are IR timers.
   if (0) {
      jgamma.Clear();
      for (size_t iShF = 0; iShF != pFitBasis->Shells.size(); ++ iShF){
         ir::FRawShell const &ShF = pFitBasis->Shells[iShF];
         size_t nFnF = ShF.nFn();

         double
            // (\mu\nu|F): nAo x nAo x nFnF
            *pMNF = FormIntMNF(ShF, pIntKernel, pOrbBasis, pFitBasis, Mem_, ScrDen, 0.);

         // jgamma[A] = (\mu\nu|A) Density[\mu,\nu]
         Mxv(&jgamma(pFitBasis->iFn(iShF),0), 1,
            pMNF, nAo*nAo, 1,
            &Density(0,0), 1, nFnF, nAo*nAo);

         Mem_.Free(pMNF);
      }
   } else { TIME_SECTION(m_pTimers, 0x202, "DF-J (density fit)");
      FMemoryStackArray MemStacks(Mem_);
      FOmpAccBlock
         jgamma_acc(jgamma.pData, nFit, OMPACC_ClearTarget, Mem_);
      #pragma omp parallel for schedule(dynamic,1)
      for (int iShA_ = 0; iShA_ < int(pOrbBasis->Shells.size()); ++ iShA_) {
         FMemoryStack
            &Mem = MemStacks.GetStackOfThread();
         size_t iShA = size_t(iShA_);
//          for (size_t iShB = 0; iShB <= iShA; ++ iShB) {
         for (size_t iShB = iShA; iShB < pOrbBasis->Shells.size(); ++ iShB) {
            ir::FRawShell
               *pShA = const_cast<ir::FRawShell*>(&pOrbBasis->Shells[iShA]),
               *pShB = const_cast<ir::FRawShell*>(&pOrbBasis->Shells[iShB]);
            double const
               *pDenAB = &Density(pOrbBasis->iFn(iShA), pOrbBasis->iFn(iShB));
            // hm... this one is *really* show... takes almost twice as much time
            // as the ContractF below. Wasn't this one supposed to be the fast one?
            // Should probably have a look at that...
            EvalInt2e3c_ContractAB(jgamma_acc.pTls(), pShA, pShB, pDenAB, 1, nAo,
               const_cast<ir::FRawShell*>(&pFitBasis->Shells[0]), pFitBasis->Shells.size(), (iShA==iShB)? 1. : 2., &ir::g_IrCoulombKernel, Mem);
         }
      }
      jgamma_acc.Join();
   }


   //   Solve[Jcd[AB] D[\nu B i]] -> D[\nu A I]
   CholeskySolve(jgamma, Jcd);

   // recalculate integrals to form the coulomb matrix.
   { TIME_SECTION(m_pTimers, 0x203, "DF-J (coul. contrib.)");
      Coul.Clear();
      FormIntMNF_ContractF(Coul, pIntKernel, pOrbBasis, pFitBasis, jgamma.pData, Mem_);
      EnergyCoul = 0.5 * Dot(Coul, Density);
   }
   PAUSE_CLOCK(300)
}


void FFockComponentBuilderDfCoul::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem)
{
   FFockComponentBuilder::AccFock(FockC, FockO, pBasis, COccC, COccO, Flags, Mem); // consistency checks.

   // todo: add options to MakeJk to *update* the fock matrix instead of overwriting it.
   size_t
      nBf = pOrbBasis->nFn();
   FStackMatrix
      CoulC(nBf, nBf, &Mem);
   MakeCoul(CoulC, this->Energy, &ir::g_IrCoulombKernel, &*pBasis->pRawBasis, &*pFitBasis->pRawBasis, COccC, COccO, Jcd, Mem, m_pTimers);
   if (JkOptions.UseBandStructureEnergyFormula) {
      this->Energy *= -1.;
   }
   Add(FockC, CoulC);
}

void FFockComponentBuilderDfCoul::PrintEnergyContribs()
{
   m_Log.WriteResult("Coulomb energy", Energy);
}





// #ifdef INCLUDE_OPTIONALS

bool FGridRenorm::IsGoodFor(FDftGrid const *pGrid_, FRawBasis const *pBasis_) const
{
   return pBasis_ == m_pBasis && pGrid_ == m_pGrid;
}

FGridRenorm::FGridRenorm(FRawBasis *pBasis, FDftGrid *pGrid, FXcOptions const *pXcOpt, FLog *pLog, FMemoryStack &Mem)
   : m_pGrid(pGrid), m_pBasis(pBasis), nBf(pBasis->nFn()), m_pLog(pLog)
{
   m_pSaCd = new FHeapMatrix(nBf, nBf);

   // first make the analytic overlap matrix.
   MakeIntMatrix(*m_pSaCd, *m_pBasis, *m_pBasis, FKrn2i_Direct(&ir::g_IrOverlapKernel), Mem);

   // now the numerical approximation.
   m_pSgCd = new FHeapMatrix(nBf, nBf);
   m_pSgCd->Clear();
   for (int iGridBlock = 0; iGridBlock < (int)m_pGrid->GridBlocks.size(); ++ iGridBlock) {
      FDftGrid::FGridBlock
         *pGridBlock = &m_pGrid->GridBlocks[iGridBlock];
      size_t
         nGridPt = pGridBlock->nPt();
      double const
         *pGridPt = &m_pGrid->Positions[pGridBlock->iFirst][0],
         *pGridWt = &m_pGrid->Weights[pGridBlock->iFirst];

      FStackMatrix
         // actual size: nGridPt x nMap (with nMap <= nBf);
         OrbVal(nGridPt, nBf, &Mem);
      TMemoryLock<size_t>
         pMap(nBf, &Mem);
      size_t
         nMap;
      dfti::IrEvalBfn(OrbVal.pData, pMap, nMap, pGridPt, 3, nGridPt, dfti::DERIVCOMP_ValueOnly,
         m_pBasis, pXcOpt->fThrOrb, MakeLogThrOrb(pXcOpt->fThrOrb), Mem);
      if (nMap != 0) {
         OrbVal.nCols = nMap;
         // absorb weights into grid dimension.
         for (size_t iGridPt = 0; iGridPt < nGridPt; ++ iGridPt) {
            double fWtNrm = std::sqrt(pGridWt[iGridPt]);
            for (size_t iMap = 0; iMap < nMap; ++ iMap)
               OrbVal(iGridPt, iMap) *= fWtNrm;
         }
         // and add to overlap matrix.
         FStackMatrix
            gS(nMap, nMap, &Mem);
         SyrkTN(gS, OrbVal);
   //       SyrkTN(gS, OrbVal, 1.0, MXM_Add);
         for (size_t jj = 0; jj < nMap; ++ jj) {
            for (size_t ii = 0; ii <= jj; ++ ii) {
               size_t
                  i = pMap[ii],
                  j = pMap[jj];
               double v = gS(ii,jj);
               (*m_pSgCd)(i,j) += v;
               if (i != j)
                  (*m_pSgCd)(j,i) += v;
            }
         }
      }
   }
//    m_pSgCd->Print(xout, "Numerical overlap matrix.");
//    m_pSaCd->Print(xout, "Analytical overlap matrix.");

   if (1) {
      double
         fRmsd = 0.;
      for (size_t j = 0; j < nBf; ++ j)
         for (size_t i = 0; i < nBf; ++ i)
            fRmsd += sqr((*m_pSaCd)(i,j) - (*m_pSgCd)(i,j));
      fRmsd = std::sqrt(fRmsd)/double(nBf);
      m_fOverlapRmsd = fRmsd;
   }


   // make cholesky decompositions.
   CalcCholeskyFactors(*m_pSaCd);
   CalcCholeskyFactors(*m_pSgCd);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FGridRenorm::RenormalizeRows(FMatrixView InAndOut, FRenormDirection Dir, FCoeffType CoeffType, FRowVectorType RowVectorType, FMemoryStack &Mem)
{
   assert(Dir == RENORM_GridToAnalytic || Dir == RENORM_AnalyticToGrid);
   assert(CoeffType == COEFF_Covariant || CoeffType == COEFF_Contravariant);
   assert(RowVectorType == VECTOR_Orbital|| RowVectorType == VECTOR_Density);
   assert_rt(InAndOut.nRows == nBf);
   FMatrixView
      *pS1cd, *pS2cd;
   pS1cd = m_pSaCd.get(),
   pS2cd = m_pSgCd.get();
   assert(pS1cd->nRows == nBf && pS1cd->IsSquare() && pS2cd->nRows == nBf && pS2cd->IsSquare());
   if (Dir == RENORM_GridToAnalytic)
      std::swap(pS1cd, pS2cd);
//    pS1cd = pS2cd; // FIXME: remove this (testing: should leave vector invariant)

   if (RowVectorType == VECTOR_Density) {
      if (CoeffType == COEFF_Contravariant) {
         // transform expansion vector coefficients: C' = S'^{-1} S C
         CholeskyMxm(InAndOut, *pS1cd);
         CholeskySolve(InAndOut, *pS2cd);
      } else {
         // transform potential vector coefficients: V' = S' S^{-1} V
         assert(CoeffType == COEFF_Covariant);
         CholeskySolve(InAndOut, *pS1cd);
         CholeskyMxm(InAndOut, *pS2cd);
      }
   } else {
//       return; // FIXME: remove this.
      assert(RowVectorType == VECTOR_Orbital);
      if (CoeffType == COEFF_Contravariant) {
         // transform expansion vector coefficients: C' = solve(L'^T, L^T C)
         TriangularMxm(InAndOut, *pS1cd, 'L', 'T');
         TriangularSolve(InAndOut, *pS2cd, 'L', 'T');
      } else {
         // transform potential vector coefficients: V' = L' solve(L, V)
         // note:
         // - (L L^T)^{-1} = (L^T)^{-1} L^{-1}
         // - goal: C^T V = C'^T V'
         //   Let's calculate:
         //   C'^T = (solve(L'^T, L^T C))^T
         //        = ((L'^T)^{-1} L^T C)^T
         //        = C^T L L'^{-1}
         //   So if V' = L' solve(L, V), we get for the contraction:
         //   C'^T V' = (C^T L L'^{-1}) (L' L^{-1} V)
         //           = C^T (L (L'^{-1} L') L^{-1}) V
         //           = C^T V
         // --> all okay. Seems correct. But it doesn't work particularly well.
         // In fact, this appears to do very little at all.
         // I wonder if I could be missing some derivative-related terms? I don't *think* so...
         //
         // test:
         // make && microscf -t '!DFDJ-RKS/PBE/def2-TZVPP scf{thr-orb:1e-10} grid{renorm:orb;level:2}' -g ~/molecules/acrylic_acid.xyz

         assert(CoeffType == COEFF_Covariant);
         TriangularSolve(InAndOut, *pS1cd, 'L', 'N');
         TriangularMxm(InAndOut, *pS2cd, 'L', 'N');
      }
   }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


// #endif // INCLUDE_OPTIONALS




FFockComponentBuilderXc::FFockComponentBuilderXc(FXcOptions const &XcOptions_, FLog &Log_, FTimerSet *pTimers_)
   : FFockComponentBuilder(Log_, pTimers_), XcOptions(XcOptions_)
{
   if (XcOptions.pXcFn.get()) {
      m_Log.Write("\n"+XcOptions.pXcFn->Desc());
      if (XcOptions.pXcFn->NeedTau() && XcOptions.UseDfxc)
         throw std::runtime_error(fmt::format("XC functional '{}' needs tau inputs (kinetic energy density). Sorry, but tau cannot be computed with DFXC=1 (use a laplace-level (LL-*) functional instead).", XcOptions.pXcFn->Name()));
   }
   EnergyXc = 0.;
   fElecTotal[0] = 0.;
   fElecTotal[1] = 0.;
   AuxExpandXc = XcOptions.UseDfxc;
   m_BaseAccXcFlags = 0;
   if (AuxExpandXc)
      m_BaseAccXcFlags |= dfti::DFTI_AuxiliaryExpand;
//    AuxExpandXc = false;
//    m_Log.Write(" !!xc/j builder initialized with grid level: {}/{}", GridParams.nLevel, GridParams_.nLevel);
   if (XcOptions.UseBandStructureEnergyFormula)
      m_BaseAccXcFlags |= dfti::DFTI_CounterDensityXcIntegral;
}


FFockComponentBuilderDfCoulXcCached::FFockComponentBuilderDfCoulXcCached(FDfJkOptions const &JkOptions_, FXcOptions const &XcOptions_, FLog &Log_, FTimerSet *pTimers_)
   : FFockComponentBuilderXc(XcOptions_, Log_, pTimers_), JkOptions(JkOptions_)
{
   EnergyCoulomb = 0.;
}


FFockComponentBuilderXc::~FFockComponentBuilderXc()
{
}

FFockComponentBuilderDfCoulXcCached::~FFockComponentBuilderDfCoulXcCached()
{
}



void FFockComponentBuilderXc::Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem)
{
   m_WfDecl = WfDecl_;
   fElecTotalAnalytic[0] = double(WfDecl_.nElec());
   fElecTotalAnalytic[1] = double(WfDecl_.nOpen());

   pAtoms = &Atoms_;
   pOrbBasis = pOrbBasis_;
   nAo = pOrbBasis->nFn();
   // compute size of packed triangular matrices as used in DFTI interface
   nAoTr = (nAo*(nAo+1))/2;

   // instanciate the molecular integration grid.
   SwitchToRefineGrid(*XcOptions.pGridParams);

   if (XcOptions.GridRenormType == FDftGridParams::GRIDRENORM_Orbitals) {
      FTimer TimerGridRenorm;
      m_pGridRenormOrb = new FGridRenorm(pOrbBasis->pRawBasis.get(), pDftGrid.get(), &XcOptions, &m_Log, Mem);
      m_Log.WriteTiming("grid renorm trafo", (double)TimerGridRenorm, fmt::format("rmsd(S): {:.2e}", m_pGridRenormOrb->fOverlapRmsd()));
   }

   m_Log.CheckStatus(); // may raise exception.
   IR_SUPPRESS_UNUSED_WARNING(Options_);
   IR_SUPPRESS_UNUSED_WARNING(WfDecl_);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FFockComponentBuilderDfCoulXcCached::Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem)
{
   FFockComponentBuilderXc::Init(WfDecl_, pOrbBasis_, Atoms_, Options_, Mem);

   InitFittingBasisAndJcd(pFitBasis, Jcd, &Atoms_, BASIS_JFit, &Options_.Print, &m_Log, Mem);
   nFit = pFitBasis->nFn();

   FRawBasis
      *pOrbBasisRaw = pOrbBasis->pRawBasis.get(),
      *pFitBasisRaw = pFitBasis->pRawBasis.get();



   // make the cached integrals themselves.
   assert(nAoTr == (nAo*(nAo+1))/2);
   nAoTrSt = AlignSizeT(nAoTr, CX_DEFAULT_MEM_ALIGN/sizeof(double));
   Int3ixStorage.resize(nAoTrSt * nFit);
   Int3ix = FMatrixView(&Int3ixStorage[0], nAoTr, nFit, 1, nAoTrSt);

   FTimer Timer3ix;
   FMatrixView ScrDen(0,0,0);
   {
      FMemoryStackArray MemStacks(Mem);
      #pragma omp parallel for schedule(dynamic)
      for (int iShF__ = 0; iShF__ < int(pFitBasisRaw->Shells.size()); ++ iShF__) {
         size_t iShF = size_t(iShF__); // that one's for OpenMP.
         if (m_Log.StatusOkay()) {
            FMemoryStack &Mem1 = MemStacks.GetStackOfThread();
            ir::FRawShell const &ShF = pFitBasisRaw->Shells[iShF];
            size_t nFnF = ShF.nFn();

            double
               // (\mu\nu|F): nAo x nAo x nFnF
               *pMNF = FormIntMNF(ShF, &ir::g_IrCoulombKernel, pOrbBasisRaw, pOrbBasisRaw, pFitBasisRaw, Mem1, ScrDen, 0., INTFLAGS_SkipTriangleFix);
//                *pMNF = FormIntMNF(ShF, &ir::g_IrCoulombKernel, pOrbBasisRaw, pOrbBasisRaw, pFitBasisRaw, Mem1, ScrDen, 0., 0);

            // hmpf. copying around all this stuff takes longer than evaluating the integrals themselves.
            // And the profiler did not even show it. Note that this is the slower variant, with inverted loops in FormIntMNF.
            // However, this part should later anyway be replaced by cutting the integrals into (ab|-block x F pieces x nAbBlocks
            // pieces. Current format is inherently problematic for CPU scaling (also in Mvx used for forming the matrices).
            assert(nAoTrSt == Int3ix.nColSt);
            size_t iTriOffs = 0;
            for (size_t iCol = 0; iCol < nAo; ++ iCol) {
               double
                  *IR_RP s = &pMNF[iCol * nAo],
                  *IR_RP d = &Int3ix(iTriOffs, pFitBasisRaw->iFn(iShF));
               for (size_t iF = 0; iF < nFnF; ++ iF)
                  for ( size_t iRow = 0; iRow <= iCol; ++ iRow )
                     d[iRow + nAoTrSt * iF] = s[iRow + nAo*nAo * iF];
               iTriOffs += iCol + 1;
            }
            Mem1.Free(pMNF);
         }
      }
   }
   m_Log.CheckStatus(); // may raise exception.

   double
      fIntSizeMb = double(Int3ix.GetStridedSize()) * double(sizeof(Int3ix[0])) / double(1<<20);
   m_Log.WriteTiming(fmt::format("3-index integrals ({} MB)", (size_t)fIntSizeMb), (double)Timer3ix);
   IR_SUPPRESS_UNUSED_WARNING(Options_);
   IR_SUPPRESS_UNUSED_WARNING(WfDecl_);
}





void FFockComponentBuilderXc::MakeDensityMatrices(FMatrixView &Rdm, FMatrixView &RdmO, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   size_t
      nAo = pOrbBasis->nFn();

   assert((nAo == COccC.nRows) && (COccO.pData == 0 || COccO.nRows == nAo));
   // make total and open-shell density matrix, and re-pack orbitals for DFTI (if required)
   Rdm = MakeStackMatrix(nAo, nAo, Mem);
   RdmO = FMatrixView(0,0,0);

   SyrkNT(Rdm, COccC);
   if (COccO.nCols != 0) {
      RdmO = MakeStackMatrix(nAo, nAo, Mem);
      SyrkNT(RdmO, COccO);
      Add(Rdm, RdmO);
   }
//    xout << fmt::format("FFockComponentBuilderXc::MakeDensityMatrices: nOccC = {}  nOccO = {}", COccC.nCols, COccO.nCols) << std::endl;
}


// fixme: this should probably make separate open- and closed-shell densities. open-shell needed for xc (but not for coulomb).
void FFockComponentBuilderDfCoulXcCached::Make1ixDensity1(FMatrixView jgamma, FMatrixView AuxDen, FMatrixView RdmC, FMatrixView RdmO, size_t nDenVec, FMemoryStack &Mem)
{
   assert(AuxDen.nCols == nDenVec && jgamma.nCols == nDenVec);
   FStackMatrix
      DensityTr(nAoTr, nDenVec, &Mem);
   m_pTimers->Resume(0x201, "DF-J/XC (RDM)");

   RdmC.TriangularReduce(1, &DensityTr(0,0), 2.); // reduce to triangular form; multiply off diagonal elements by 2.0
   if (nDenVec >= 2)
      RdmO.TriangularReduce(1, &DensityTr(0,1), 2.);
   m_pTimers->Pause(0x201, "DF-J/XC (RDM)");

   m_pTimers->Resume(0x202, "DF-J/XC (1x density fit)");
   if (nDenVec == 1)
      Mxva(jgamma.pData, Transpose(Int3ix), DensityTr.pData);
   else
      Mxm(jgamma, Transpose(Int3ix), DensityTr);

   Assign(AuxDen, jgamma);

   //   Solve[Jcd[AB] D[\nu B i]] -> D[\nu A I]  to get density coefficients
   CholeskySolve(AuxDen, Jcd);
   m_pTimers->Pause(0x202, "DF-J/XC (Density fit)");
}


void FFockComponentBuilderXc::ComputeBackgroundDensities(FMatrixView CEnvOccC, FMatrixView CEnvOccO, FMemoryStack &Mem)
{
   using namespace dfti;
   TMemoryLock<double>
      pFreeMe(0, &Mem);

   if (!m_BackgroundDensities.empty())
      throw std::runtime_error("ComputeBackgroundDensities: background densities already assigned.");
   m_BackgroundDensities.resize(pDftGrid->Points.size() * DFTI_nBackgroundDensityComponents);

   if (AuxExpandXc)
      throw std::runtime_error("ComputeBackgroundDensities: sorry, not implemented with DFXC=1.");

   if (nAo != CEnvOccC.nRows)
      throw std::runtime_error("ComputeBackgroundDensities: Input orbitals not consistent with orbital basis set.");

   m_Log.Write(" Computing environment background density for {} closed-shell electrons and {} open-shell electrons.", 2*CEnvOccC.nCols, CEnvOccO.nCols);

   // compute total densities of the environment orbitals on all grid points
   // and store them in m_BackgroundDensities.
   FMatrixView
      EnvRdm, EnvRdmO;
   MakeDensityMatrices(EnvRdm, EnvRdmO, CEnvOccC, CEnvOccO, Mem);

   if (XcOptions.pXcFn && !AuxExpandXc) {
      using namespace dfti;
      TMemoryLock<double> pFreeMe4(0, &Mem);
      FDftiArgs DftiArgs = {
         DFTI_MakeBackgroundDensities | m_BaseAccXcFlags,
         &EnergyXc, 0, 0, 0,
         0, 0, 0, 0, 0, 0, // density/orbs (assigned next)
         pOrbBasis->pRawBasis.get(),
         XcOptions.pXcFn.get(),
         pDftGrid.get(),
         XcOptions.fThrOrb, MakeLogThrOrb(XcOptions.fThrOrb), 0, &fElecTotal[0], 0 // 1e-1 * ThrDen?
      };
      PackAndAssignDensities(DftiArgs, EnvRdm, EnvRdmO, CEnvOccC, CEnvOccO, Mem);
      AccXc(DftiArgs, m_Log, m_pTimers, Mem);
      m_Log.Write(" Note: Background electron number integrates to {:.6f} numerically", fElecTotal[0]);
   }

   // tell AccXc to add the background densities in subsequent DFT integrations
   m_BaseAccXcFlags |= DFTI_AddBackgroundDensities;
}




void FFockComponentBuilderXc::PackAndAssignDensities(dfti::FDftiArgs &DftiArgs, FMatrixView RdmC, FMatrixView RdmO,
      FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   double
      *pRdmTrC = 0,
      *pRdmTrO = 0;
   Mem.Alloc(pRdmTrC, nAoTr);
   RdmC.TriangularReduce(1, pRdmTrC, 1.);
   if (COccO.nCols != 0) {
      Mem.Alloc(pRdmTrO, nAoTr);
//       Scale(RdmO, std::sqrt(0.5));
      RdmO.TriangularReduce(1, pRdmTrO, 1.);
//       memset(pRdmTrO, 0, sizeof(double) * nAoTr);
   }
   FMatrixView
      OccOrbC(0,0,0),
      OccOrbO(0,0,0); // note: nOcc x nAo with absorbed occupation numbers.
   OccOrbC = MakeStackMatrix(COccC.nCols+COccO.nCols, COccC.nRows, Mem);
   Assign(Select(OccOrbC,0,0,COccC.nCols,COccC.nRows), Transpose(COccC));
   Assign(Select(OccOrbC,COccC.nCols,0,COccO.nCols,COccC.nRows), Transpose(COccO));
   if (COccO.nCols != 0) {
      OccOrbO = MakeStackMatrix(COccO.nCols, COccO.nRows, Mem);
      Assign(OccOrbO, Transpose(COccO));
   }

   DftiArgs.pDenC = pRdmTrC;
   DftiArgs.pDenO = pRdmTrO;
   if (1) {
      // FIXME: atm open-shell is broken...
      DftiArgs.pOccOrbC = OccOrbC.pData;
      DftiArgs.pOccOrbO = OccOrbO.pData;
      DftiArgs.nOccC = COccC.nCols + COccO.nCols;
      DftiArgs.nOccO = COccO.nCols;
   }
   DftiArgs.pBackgroundDensities = m_BackgroundDensities.data(); // may be 0.
   if (DftiArgs.Flags & dfti::DFTI_AddBackgroundDensities)
      m_Log.Write(" Including background densities in xc calculation for {} closed-shell electrons and {} open-shell electrons.", 2*COccC.nCols, COccO.nCols);
}


double FFockComponentBuilderXc::ComputeAuxiliaryElectronNumbers(FMatrixView AuxDen, FRawBasis const *pFitBasisRaw)
{
   TIME_SECTION(m_pTimers, 0x203, "DF-J/XC (nelec)");
   // compute number of electrons in the auxiliary basis. Since we fit densities,
   // this number may not be exact. Note that only s functions carry electrons.
   // Mathematica says: Assuming[\[Alpha] > 0, Integrate[Exp[-\[Alpha]*r^2]*r^2*4*Pi, {r, 0, Infinity}]] = (pi/alpha)^(3/2)
   fElecTotalAnalytic[0] = 0.;
   fElecTotalAnalytic[1] = 0.;
   for (size_t iSh = 0; iSh < pFitBasisRaw->Shells.size(); ++ iSh) {
      ir::FRawShell const &Sh = pFitBasisRaw->Shells[iSh];
      if (Sh.l != 0) continue;
      for (size_t iExp = 0; iExp < Sh.nExp; ++ iExp) {
         double fDen = std::pow(M_PI/Sh.pExp[iExp], 1.5);
         for (size_t iDenVec = 0; iDenVec < AuxDen.nCols; ++ iDenVec)
         {
            double const *pCoeff = &AuxDen(pFitBasisRaw->iFn(iSh), iDenVec);
            for (size_t iCo = 0; iCo < Sh.nCo; ++ iCo)
               fElecTotalAnalytic[iDenVec] += pCoeff[iCo] * Sh.pCo[iExp + Sh.nExp*iCo] * fDen;
         }
      }
   }
   return fElecTotalAnalytic[0];
}


void ExpandRenormalizeAndAdd_Fock(FMatrixView FockOut, double const *pFockInTr, double Factor, FGridRenorm *pGridRenorm, FMemoryStack &Mem)
{
   int iSign = +1; // <-- +1: symmetric, -1: anti-symmetric.
   if (pGridRenorm == 0) {
      Add_TriangularExpanded(FockOut, iSign, pFockInTr, Factor);
   } else {
      assert(FockOut.nRows == FockOut.nCols);
      FStackMatrix
         FockInSq(FockOut.nRows, FockOut.nCols, &Mem);
      FockInSq.TriangularExpand(iSign, pFockInTr);
      assert_rt(FockInSq.nRows == FockInSq.nCols);
      pGridRenorm->RenormalizeRows(FockInSq, FGridRenorm::RENORM_GridToAnalytic, FGridRenorm::COEFF_Covariant, FGridRenorm::VECTOR_Orbital, Mem);
      pGridRenorm->RenormalizeRows(Transpose(FockInSq), FGridRenorm::RENORM_GridToAnalytic, FGridRenorm::COEFF_Covariant, FGridRenorm::VECTOR_Orbital, Mem);
      Add(FockOut, FockInSq, Factor);
   }
}


void FFockComponentBuilderXc::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC_, FMatrixView const &COccO_, uint Flags, FMemoryStack &Mem)
{
   TMemoryLock<double>
      pBaseOfMemory(0, &Mem);
   TIME_SECTION(m_pTimers, 0x200, "XC");
   FFockComponentBuilder::AccFock(FockC, FockO, pBasis, COccC_, COccO_, Flags, Mem); // consistency checks.

   if (pBasis != pOrbBasis)
      throw std::runtime_error("FFockComponentBuilderXc must now be used with orb-projected initial guess (not Fock-projected guess!)");

   if (nAo != COccC_.nRows || nAo != COccO_.nRows)
      throw std::runtime_error("AccFock: Input orbitals not consistent with orbital basis set.");

   FMatrixView
      COccC = COccC_, COccO = COccO_;
// #ifdef INCLUDE_OPTIONALS
   if (m_pGridRenormOrb.get()) {
      // copy input orbital matrices and renormalize them
      assert(COccO.pData == 0 || COccO.nRows == COccC.nRows);
      FMatrixView COcc = MakeStackMatrix(nAo, COccC.nCols + COccO.nCols, Mem);
      COccC = Select(COcc, 0, 0, nAo, COccC.nCols);
      COccO = Select(COcc, 0, COccC.nCols, nAo, COccO.nCols);
      Assign(COccC, COccC_);
      Assign(COccO, COccO_);
      m_pGridRenormOrb->RenormalizeRows(COcc, FGridRenorm::RENORM_AnalyticToGrid, FGridRenorm::COEFF_Contravariant, FGridRenorm::VECTOR_Orbital, Mem);
   }
// #endif // INCLUDE_OPTIONALS

   FMatrixView
      Rdm, RdmO;
   { TIME_SECTION(m_pTimers, 0x201, "XC (RDM)");
      MakeDensityMatrices(Rdm, RdmO, COccC, COccO, Mem);
   }

   bool
      OpenShell = (COccO.nCols != 0);
   size_t
      nFock = (OpenShell? 2 : 1);
   if (XcOptions.pXcFn && AuxExpandXc)
      throw std::runtime_error("FFockComponentBuilderXc: AccFock of this should not be called in AuxExpandXc==false case.");

   if (!XcOptions.pXcFn) {
      // nothing to do in this case.
      EnergyXc = 0;
   } else {
      TIME_SECTION(m_pTimers, 0x204, "XC (xc contrib.)");
      TMemoryLock<double> pFreeMe4(0, &Mem);

      FStackMatrix
         // resulting triangular packed Fock matrices from DFTI
         FockTr(nAoTr, nFock, &Mem);
      FockTr.Clear();

      using namespace dfti;

      FDftiArgs DftiArgs = {
         DFTI_MakeXc | m_BaseAccXcFlags,
         &EnergyXc, 0, &FockTr(0,0), (OpenShell? (&FockTr(0,1)) : 0),
         0, 0, 0, 0, 0, 0, // density/orbs (assigned next)
         pOrbBasis->pRawBasis.get(),
         XcOptions.pXcFn.get(),
         pDftGrid.get(),
         XcOptions.fThrOrb, MakeLogThrOrb(XcOptions.fThrOrb), 0, &fElecTotal[0], 0 // 1e-1 * ThrDen?
      };
      PackAndAssignDensities(DftiArgs, Rdm, RdmO, COccC, COccO, Mem);

      AccXc(DftiArgs, m_Log, m_pTimers, Mem);

      // unpack the fock matrices and add them to their target locations.
      ExpandRenormalizeAndAdd_Fock(FockC, &FockTr(0,0), 1.0, m_pGridRenormOrb.get(), Mem);
      if (OpenShell)
         ExpandRenormalizeAndAdd_Fock(FockO, &FockTr(0,1), 1.0, m_pGridRenormOrb.get(), Mem);
   }
   Energy = EnergyXc;
}


void FFockComponentBuilderDfCoulXcCached::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem)
{
   TIME_SECTION(m_pTimers, 0x200, "DF-J/XC");
   // NOTE:
   // - this *DOES NOT* call FFockComponentBuilder::AccFock!
   // - we could do it if AuxExpandXc is false, but we'd need to find some
   //   way of passing intermediate data in this case, with some sort
   //   of indirected interface (on input, the RDMs, on output, the triangular
   //   Fock matrices).
   // - ...probably not worth it, and easier to read with the code duplications.
   //   not that big anyway.
   FFockComponentBuilder::AccFock(FockC, FockO, pBasis, COccC, COccO, Flags, Mem);

   if (pBasis != pOrbBasis)
      throw std::runtime_error("FFockComponentBuilderDfCoulXcCached must now be used with orb-projected initial guess (not Fock-projected guess!)");

   if (nAo != COccC.nRows || nAo != COccO.nRows || nFit != Jcd.nRows)
      throw std::runtime_error("AccFock: Input orbitals not consistent with orbital basis set.");
   if (m_pGridRenormOrb)
      throw std::runtime_error("AccFock: Grid renormalization is an experimental option. It is not implemented in the DF-J/DF-JX (cached 3ix) Fock builder.");

   FMatrixView
      Rdm, RdmO;
   MakeDensityMatrices(Rdm, RdmO, COccC, COccO, Mem);

   bool
      OpenShell = (COccO.nCols != 0);
   size_t
      nDenVec = (OpenShell? (AuxExpandXc? 2 : 1) : 1),
      nFock = (OpenShell? 2 : 1);
   FStackMatrix
      jgamma(nFit, nDenVec, &Mem),   // (A|rs) gamma[r,s].
      AuxDen(nFit, nDenVec, &Mem);   // (A|B)^{-1} jgamma[B]
   Make1ixDensity1(jgamma, AuxDen, Rdm, RdmO, nDenVec, Mem);

   if (1 && AuxExpandXc)
      ComputeAuxiliaryElectronNumbers(AuxDen, pFitBasis->pRawBasis.get());

   if (1) {
      // make copies of the auxiliary densities --- we might still need
      // them in the gradient evaluation.
      m_LastDensity.resize(nFit * nDenVec);
      assert(sizeof(AuxDen[0]) == sizeof(m_LastDensity[0]));
      memcpy(&m_LastDensity[0], &AuxDen[0], sizeof(AuxDen[0]) * nFit * nDenVec);
   }

   if (1) {
      TIME_SECTION(m_pTimers, 0x203, "DF-J/XC (coul. energy)");
      // make coulomb energy from 2ix integrals:
      // Ecoulomb = .5 j[A] (A|B) j[B]
      // (since we use a robust fit this is equal to the density matrix contraction.
      // the advantage here is that we can do it also in the xc auxiliary expansion case,
      // where we will not get a pure j matrix).
      EnergyCoulomb = .5 * Dot(AuxDen.pData, jgamma.pData, nFit);
   }

   EnergyXc = 0;
   if (XcOptions.pXcFn && AuxExpandXc) {
      TIME_SECTION(m_pTimers, 0x204, "DF-J/XC (xc contrib.)");
      using namespace dfti;
      FStackMatrix
         vxc(nFit, nDenVec, &Mem);;

      vxc.Clear();

      FDftiArgs DftiArgs = {
         DFTI_MakeXc | m_BaseAccXcFlags,
         &EnergyXc, 0, &vxc(0,0), (OpenShell? (&vxc(0,1)) : 0),
         &AuxDen(0,0), (OpenShell? (&AuxDen(0,1)) : 0), 0, 0, 0, 0, // no orbitals provided (not needed in dfxc=1 case).
         pFitBasis->pRawBasis.get(),
         XcOptions.pXcFn.get(),
         pDftGrid.get(),
         XcOptions.fThrOrb, MakeLogThrOrb(XcOptions.fThrOrb), 0, &fElecTotal[0], 0 // 1e-1 * ThrDen?
      };
      AccXc(DftiArgs, m_Log, m_pTimers, Mem);

      CholeskySolve(vxc, Jcd);
      // add auxiliary density to coulomb potential, for total-density part.
      // for other part, assign vxc to AuxDen.
      Add(AuxDen, vxc);
      if (nDenVec >= 2) {
         for (size_t i = 0; i < nFit; ++ i)
            AuxDen(i,1) = vxc(i,1);
      }
   }

   // recalculate integrals to form the coulomb matrix.
   FStackMatrix
      FockTr(nAoTr, nFock, &Mem);

   m_pTimers->Resume(0x205, "DF-J/XC (j/vxc matrix)");
   if (nDenVec == 1) {
      Mxva(FockTr.pData, Int3ix, AuxDen.pData);
      if (OpenShell)
         memset(&FockTr(0,1), 0, sizeof(double)*FockTr.nRows);
   } else
      Mxm(FockTr, Int3ix, AuxDen);
   m_pTimers->Pause(0x205, "DF-J/XC (j/vxc matrix)");


   if (XcOptions.pXcFn && !AuxExpandXc) {
      TIME_SECTION(m_pTimers, 0x204, "DF-J/XC (xc contrib.)");
      using namespace dfti;
      TMemoryLock<double> pFreeMe4(0, &Mem);

      FDftiArgs DftiArgs = {
         DFTI_MakeXc | m_BaseAccXcFlags,
         &EnergyXc, 0, &FockTr(0,0), (OpenShell? (&FockTr(0,1)) : 0),
         0, 0, 0, 0, 0, 0, // density/orbs (assigned next)
         pOrbBasis->pRawBasis.get(),
         XcOptions.pXcFn.get(),
         pDftGrid.get(),
         XcOptions.fThrOrb, MakeLogThrOrb(XcOptions.fThrOrb), 0, &fElecTotal[0], 0 // 1e-1 * ThrDen?
      };
      PackAndAssignDensities(DftiArgs, Rdm, RdmO, COccC, COccO, Mem);

      AccXc(DftiArgs, m_Log, m_pTimers, Mem);
   }

   if (XcOptions.UseBandStructureEnergyFormula)
      // in this case coulomb energy is used to counter twice-counted e-e
      // energy from band structure term of total energy. So now subtract
      // it instead of adding it.
      EnergyCoulomb *= -1.0;

   Energy = EnergyCoulomb + EnergyXc;

   // unpack the fock matrices and add them to their target locations.
   Add_TriangularExpanded(FockC, +1, &FockTr(0,0), 1.0);
   if (OpenShell)
      Add_TriangularExpanded(FockO, +1, &FockTr(0,1), 1.0);
}

void FFockComponentBuilderXc::PrintEnergyContribs()
{
   if (fElecTotalAnalytic[1] == 0) {
      m_Log.WriteInfoExpf("Lost electrons", fElecTotal[0] - fElecTotalAnalytic[0]);
      if (pAtoms->size() != 0)
         m_Log.WriteInfoExpf(".../sqrt(nAtoms)", (fElecTotal[0] - fElecTotalAnalytic[0])/std::sqrt(double(pAtoms->size())));
      m_Log.WriteLine();
   } else {
      m_Log.WriteInfoExpf("Lost electrons", fElecTotal[0] - fElecTotalAnalytic[0], "total density");
      m_Log.WriteInfoExpf("Lost electrons", fElecTotal[1] - fElecTotalAnalytic[1], "spin density");
      m_Log.WriteLine();
   }

   m_Log.WriteResult("Density functional energy", EnergyXc);
}


void FFockComponentBuilderDfCoulXcCached::PrintEnergyContribs()
{
   FFockComponentBuilderXc::PrintEnergyContribs();
   m_Log.WriteResult("Coulomb energy", EnergyCoulomb);
}


void FFockComponentBuilderXc::AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   if (AuxExpandXc)
      throw std::runtime_error("FFockComponentBuilderXc: AccGradient of this should not be called in AuxExpandXc==false case.");

   if (XcOptions.pXcFn) {
      FMatrixView
         Rdm, RdmO;

      { TIME_SECTION(m_pTimers, 0x201, "XC (RDM)");
         MakeDensityMatrices(Rdm, RdmO, COccC, COccO, Mem);
      }

      { TIME_SECTION(m_pTimers, 0x204, "XC (xc contrib.)");
         TMemoryLock<double> pFreeMe4(0, &Mem);
         using namespace dfti;

         FDftiArgs DftiArgs = {
            uint(DFTI_MakeGradient | m_BaseAccXcFlags),
            &EnergyXc, Gradient.pData, 0, 0,
            0, 0,
            0, 0, 0, 0,
            pOrbBasis->pRawBasis.get(),
            XcOptions.pXcFn.get(),
            pDftGrid.get(),
            XcOptions.fThrOrb, MakeLogThrOrb(XcOptions.fThrOrb), 0, &fElecTotal[0], 0 // 1e-1 * ThrDen?
         };
         // pack input densities/orbitals for DFTI
         PackAndAssignDensities(DftiArgs, Rdm, RdmO, COccC, COccO, Mem);

         AccXc(DftiArgs, m_Log, m_pTimers, Mem);
      }
   }
}



void FFockComponentBuilderDfCoulXcCached::AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   TMemoryLock<double>
      pFreeMe(0, &Mem);
   bool
      OpenShell = (COccO.nCols != 0);
   size_t
      nDenVec = (OpenShell? (AuxExpandXc? 2 : 1) : 1);
   size_t
      nAo = pOrbBasis->nFn();
   FMatrixView
      Rdm, RdmO;
   MakeDensityMatrices(Rdm, RdmO, COccC, COccO, Mem);

   {
      FRawBasis const
         *pOrbBasisRaw = &*pOrbBasis->pRawBasis,
         *pFitBasisRaw = &*pFitBasis->pRawBasis;

      // 2e terms (see exchg.pdf).
      //
      // Coulomb part:
      //
      //   v[A] = J_AB^{-1} (μν|B) D^{μν}
      //   Ecoul^q = D^{μν} (µν|A)^q v[A] - (1/2) v_A (A|B)^q v_B
      //
      // (note: should use two-sided inline contracting driver. Or at least 1-sided
      //  inline contracting driver (on v[A] side))

      // make v.
      FStackMatrix
         jgamma(nFit, nDenVec, &Mem), // (A|rs) gamma[r,s].
         v3ix(nFit, nDenVec, &Mem), // (A|B)^{-1} jgamma[B]
         v2ixLeft(nFit, nDenVec, &Mem),
         v2ixRight(nFit, nDenVec, &Mem);
//       Make1ixDensity(jgamma, v3ix, COccC, COccO, Mem);
      memcpy(&v3ix[0], &m_LastDensity[0], sizeof(m_LastDensity[0]) * nFit * nDenVec);
      assert_rt(m_LastDensity.size() == nFit * nDenVec);
      // ^- hm... now, is this more or less consistent than making a new density?
      //    I guess one of them could have a first-order error and the other a second
      //    order error (in the accuracy of the orbitals), but I don't know.
      //    update: the way this is coded now... doesn't matter. ends anyway with
      //    a fock build. so I get identical results.
      //    In the other case it seems like making a new density is beneficial for
      //    bad densities, due to the better densities this produces in the xc core gradient.

      Assign(v2ixLeft,v3ix);
      Assign(v2ixRight,v3ix);
//       Transpose(v).Print(xout, "aux density");


      if (XcOptions.pXcFn) {
         TMemoryLock<double> pFreeMe4(0, &Mem);
         m_pTimers->Resume(0x204, "DF-J/XC (xc contrib.)");
         using namespace dfti;
         FStackMatrix
            vxc(nFit, nDenVec, &Mem);

         vxc.Clear();
//          AuxExpandXc = true;
         FDftiArgs DftiArgs = {
            uint((AuxExpandXc? (DFTI_AuxiliaryExpand | DFTI_MakeXc) : 0u) | DFTI_MakeGradient | m_BaseAccXcFlags),
            &EnergyXc, Gradient.pData, &vxc(0,0), (OpenShell? (&vxc(0,1)) : 0),
            AuxExpandXc? &v3ix(0,0) : 0,
            (OpenShell? (AuxExpandXc? &v3ix(0,1) : 0) : 0),
            0, 0, 0, 0,
            const_cast<FRawBasis*>(AuxExpandXc? pFitBasisRaw : pOrbBasisRaw),
            XcOptions.pXcFn.get(),
            pDftGrid.get(),
            XcOptions.fThrOrb, MakeLogThrOrb(XcOptions.fThrOrb), 0, &fElecTotal[0], 0 // 1e-1 * ThrDen?
         };
         if (!AuxExpandXc)
            // pack input densities/orbitals for DFTI, if required.
            PackAndAssignDensities(DftiArgs, Rdm, RdmO, COccC, COccO, Mem);

         AccXc(DftiArgs, m_Log, m_pTimers, Mem);
         if (AuxExpandXc) {
            CholeskySolve(vxc, Jcd);
            Add(v3ix, vxc);
            Add(v2ixLeft, vxc, 2.);
//             Add(v2ixRight, vxc, 0.);
            // ^- this needs to be added to v... but only on one side for the 2ix integrals.
         }
         m_pTimers->Pause(0x204, "DF-J/XC (xc contrib.)");
      }

      m_pTimers->Resume(0x224, "DF-J/XC (3ix gradient)");
      double
//          ExchFac = 1.,
//          DFac = -ExchFac * .5,
         CoulFac = 1.;

      AccGradient3ix_Factored(Gradient.pData,
         Rdm.pData, 1, nAo, &v3ix[0],
         pOrbBasisRaw, 0, pOrbBasisRaw->nSh(),
         pOrbBasisRaw, 0, pOrbBasisRaw->nSh(),
         pFitBasisRaw, 0, pFitBasisRaw->nSh(),
         &ir::g_IrCoulombKernel, Mem, CoulFac);

//       Gradient.Print(xout, "NUCLEAR GRADIENT (after 3ix part of e-e repulsion)");
      m_pTimers->Pause(0x224, "DF-J/XC (3ix gradient)");

      m_pTimers->Resume(0x225, "DF-J/XC (2ix gradient)");
      // ! 2nd contribution: +D[ijA] J^q[AB] D[ijB]    (exchange)
      // !                   -.5 * D[A] J^q[AB] D[B]   (coulomb)
      {
         FStackMatrix
            DAB(nFit, nFit, &Mem);
         DAB.Clear();
         // add coulomb contribution: -.5 * D[A] * D[B]  (note the .25 on the integrals)
         DGER(nFit, nFit, -2.*CoulFac, v2ixLeft.pData,1, v2ixRight.pData,1, DAB.pData, nFit);
         Symmetrize(DAB);
         AccGradient2ix(Gradient.pData, DAB, *pFitBasisRaw,
            *pFitBasisRaw, FKrn2i_Direct(&ir::g_IrCoulombKernel), Mem, .25);
         // hm... I could inline-contract this one too. Should I care?
      }
      m_pTimers->Pause(0x225, "DF-J/XC (2ix gradient)");
//       Gradient.Print(xout, "NUCLEAR GRADIENT (after 2ix part of e-e repulsion)");
   }
}


bool FFockComponentBuilderXc::SwitchToRefineGrid(FDftGridParams const &NewGridParams)
{
   // build a new DFT grid with the refinement grid parameters.
   // (but only if there actually is a functional to use it with)
   if (XcOptions.pXcFn) {
      pDftGrid = new FDftGrid(AsRawAtoms(*pAtoms), NewGridParams, &m_Log);
   }
   return true;
}






FDispCorrBuilder::FDispCorrBuilder(std::string const &DispCorrType_, std::string const &XcFunctionalName_, FLog &Log_, FTimerSet *pTimers_)
   : FFockComponentBuilder(Log_, pTimers_), m_pD3DispersionContext(0), m_DispCorrType(DispCorrType_), m_XcFunctionalName(XcFunctionalName_)
{
}


FDispCorrBuilder::~FDispCorrBuilder()
{
   delete m_pD3DispersionContext;
}


void FDispCorrBuilder::Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem)
{
   pAtoms = &Atoms_;
   assert(m_pD3DispersionContext == 0);
   m_pD3DispersionContext = new FD3DispersionContext(&m_Log, pAtoms, m_DispCorrType, m_XcFunctionalName, Mem);

   Energy = m_pD3DispersionContext->CalcEnergyAndGradient(FMatrixView(0,0,0), Mem);
   // ^- doesn't change between iterations.

   IR_SUPPRESS_UNUSED_WARNING(WfDecl_);
   IR_SUPPRESS_UNUSED_WARNING(pOrbBasis_);
   IR_SUPPRESS_UNUSED_WARNING(Options_);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FDispCorrBuilder::AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem)
{
   // Energy already set in Init() --- no need to do anything here.
   IR_SUPPRESS_UNUSED_WARNING3(FockC,FockO,pBasis);
   IR_SUPPRESS_UNUSED_WARNING3(COccC,COccO,Flags);
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}


void FDispCorrBuilder::AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem)
{
   m_pD3DispersionContext->CalcEnergyAndGradient(Gradient, Mem);
   IR_SUPPRESS_UNUSED_WARNING(COccC);
   IR_SUPPRESS_UNUSED_WARNING(COccO);
}


void FDispCorrBuilder::PrintEnergyContribs()
{
   m_Log.WriteResult("Dispersion correction", Energy, fmt::format("D3(BJ) for {}", m_XcFunctionalName));
}




} // namespace ct
