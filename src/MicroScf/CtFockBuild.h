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

#ifndef CT_RHF_FOCKBUILD_H
#define CT_RHF_FOCKBUILD_H

#include <string>
#include "CtBasisSet.h"
#include "CtAtomSet.h"
#include "CtMatrix.h"
// #include "CtDma.h"
#include "CtDftGrid.h" // for grid params.
#include "CtDftFunc.h"
#include "CtRhfOptions.h"
#include "CxTiming.h"
#include "CxIo.h"

namespace dfti {
   class FDftiArgs;
}

namespace ct {

enum FFockBuilderFlags {
   FOCKBUILD_ClosedOnly //
};

using mig::FDftGridParams;
using mig::FDftGrid;
using mig::FDftGridPtr;
using mig::FDftGridCptr;


struct FFockComponentBuilder : public FIntrusivePtrDest
{
   explicit FFockComponentBuilder(FLog &Log_, FTimerSet *pTimers_=0) : Energy(0.), m_Log(Log_), m_pTimers(pTimers_) {};

   // TODO: should the HfOptions really be here? I guess things like thresholds,
   // ExchFactor, the band energy formula thing etc are relevant here, but the
   // structure contains many other things most builders probably will not need.
   // And it impedes sharing with the semi-empirical version. Maybe split?
   virtual void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem) = 0;
   // Accumulate fock matrix contribution to FockC/FockO.
   // Notes:
   //   - OrbC and OrbO are supposed to have their respective occupation numbers absorbed,
   //     such that the closed/open densities are obtained by DenC := OrbC x OrbC.T and DenO := OrbO x OrbO.T
   //   - Both orbitals and Fock matrices are expressed in terms of pBasis. pBasis may differ from pOrbBasis given in Init().
   virtual void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   virtual ~FFockComponentBuilder() = 0;

   virtual void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);

   double Energy;
   virtual void PrintEnergyContribs();

   // returns true if this fock builder did anything here. Default implementation returns false (command ignored---no grids)
   virtual bool SwitchToRefineGrid(FDftGridParams const &NewGridParams);

   // Compute and store internally the densities of the environment orbitals.
   // Used in DFT-based embeddding calculations. Default implementation does nothing.
   virtual void ComputeBackgroundDensities(FMatrixView CEnvOccC, FMatrixView CEnvOccO, FMemoryStack &Mem);
protected:
   FLog
      &m_Log;
   FTimerSet
      *m_pTimers;
   FWfDecl
      m_WfDecl;
private:
   FFockComponentBuilder(FFockComponentBuilder const &); // not implemented
   void operator = (FFockComponentBuilder const &); // not implemented
};
typedef TIntrusivePtr<FFockComponentBuilder>
   FFockComponentBuilderPtr;

typedef std::vector<FFockComponentBuilderPtr>
   FFockComponentBuilderList;

void BuildFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FFockComponentBuilderList &BuilderList, FMemoryStack &Mem, FTimerSet *pTimers=0);
void SetupFockBuilderList(FFockComponentBuilderList &FockComponentBuilders, FHfOptions const &Options, FLog &Log, FTimerSet *pTimers);


struct FJkOptions
{
   FJkOptions() : fExchFactor(1.0), UseBandStructureEnergyFormula(false) {}
   double
      fExchFactor;
   bool
      UseBandStructureEnergyFormula;
};


// build both Coulomb and Exchange using 4-index integrals (i.e., no density fitting)
struct FFockComponentBuilder4ixJk : public FFockComponentBuilder
{
   FFockComponentBuilder4ixJk(FJkOptions const &Options, FLog &Log_, FTimerSet *pTimers_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
//    void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   void PrintEnergyContribs();
   ~FFockComponentBuilder4ixJk();
public:
   FJkOptions
      Options;
   FAtomSet const
      *pAtoms;
   FBasisSet
      *pOrbBasis;
   double
      EnergyCoulomb, EnergyExch, EnergyXc;
};



struct FDfJkOptions : public FJkOptions
{
};


struct FXcOptions
{
//    std::string
//       XcFunctionalName;
   FXcFunctionalPtr
      pXcFn;
   FDftGridParams const
      // parameter set for the initial DFT grid (may be changed in refine grid).
      // note: target object must stay alive for as long as the corresponding
      // Fock builders are alive.
      *pGridParams;
   bool
      UseBandStructureEnergyFormula;
   bool
      UseDfxc;
   double
      fThrOrb;
   FDftGridParams::FGridRenormType
      GridRenormType;

   explicit FXcOptions(std::string const &XcName, FDftGridParams const *pGridParams_);
};


// build both DF Coulomb and Exchange
struct FFockComponentBuilderDfJk : public FFockComponentBuilder
{
   FFockComponentBuilderDfJk(FDfJkOptions const &Options, FLog &Log_, FTimerSet *pTimers_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   void PrintEnergyContribs();
   ~FFockComponentBuilderDfJk();
public:
   FDfJkOptions
      Options;
   FAtomSet const
      *pAtoms;
   FBasisSetPtr
      pFitBasis;
   FBasisSet
      *pOrbBasis;
   FMatrixView
      Jcd;
   size_t
      nFit;
   double
      EnergyCoulomb, EnergyExch, EnergyXc;
};


// build DF Coulomb only
struct FFockComponentBuilderDfCoul : public FFockComponentBuilder
{
   FFockComponentBuilderDfCoul(FDfJkOptions const &JkOptions, FLog &Log_, FTimerSet *pTimers_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   void PrintEnergyContribs();
   ~FFockComponentBuilderDfCoul();
public:
   FDfJkOptions
      JkOptions;
   FAtomSet const
      *pAtoms;
   FBasisSetPtr
      pFitBasis;
   FBasisSet
      *pOrbBasis;
   FMatrixView
      Jcd;
   size_t
      nFit;
};



struct FGridRenorm : public FIntrusivePtrDest
{
   enum FRenormDirection {
      RENORM_AnalyticToGrid,
      RENORM_GridToAnalytic
   };
   enum FCoeffType {
      COEFF_Covariant,
      COEFF_Contravariant
   };
   enum FRowVectorType {
      // vectors squared give a density
      VECTOR_Orbital,
      // vectors themselves give a density
      VECTOR_Density
   };

   FGridRenorm(FRawBasis *pBasis, FDftGrid *pGrid, FXcOptions const *pXcOpt, FLog *pLog, FMemoryStack &Mem);

   bool IsGoodFor(FDftGrid const *pGrid_, FRawBasis const *pBasis_) const;
   // iDir +1: trafo from analytic rep to grid rep. -1: reverse.
   // iCov > 0: trafo contravariant coordinates, iCov < 0: covariant coordinates (functional)
   void RenormalizeRows(FMatrixView InAndOut, FRenormDirection Dir, FCoeffType CoeffType, FRowVectorType RowVectorType, FMemoryStack &Mem);

   double fOverlapRmsd() const { return m_fOverlapRmsd; }
protected:
   FDftGrid
      *m_pGrid;
   FRawBasis
      *m_pBasis;
   FHeapMatrixPtr
      m_pSgCd, // Cholesky decomposition of <A|B> computed on the grid
      m_pSaCd; // Cholesky decomposition of <A|B> computed analytically.
      // ^- note: wastes some space.. not stored in packed form.
   size_t
      nBf;
   double
      // rmsd over matrix elements of (AnalyticS - NumericalS).
      m_fOverlapRmsd;
   FLog
      *m_pLog;

};
typedef TIntrusivePtr<FGridRenorm>
   FGridRenormPtr;



// build standard or aux-expanded XC contributions to Fock matrix
struct FFockComponentBuilderXc : public FFockComponentBuilder
{
   FFockComponentBuilderXc(FXcOptions const &XcOptions_, FLog &Log_, FTimerSet *pTimers_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   ~FFockComponentBuilderXc();
   void PrintEnergyContribs();
   bool SwitchToRefineGrid(FDftGridParams const &NewGridParams);
   void ComputeBackgroundDensities(FMatrixView CEnvOccC, FMatrixView CEnvOccO, FMemoryStack &Mem);
protected:
   FXcOptions
      XcOptions;
   FAtomSet const
      *pAtoms;
   FBasisSet
      *pOrbBasis;
//    FRawBasis
//       *pOrbBasisRaw;
   size_t
      nAo,
      nAoTr; // triangular packed storage size of a symmetric (nAo,nAo)-shape matrix (as used in DFTI)
   FDftGridPtr
      pDftGrid;
   FGridRenormPtr
      // grid renormalization transform for orbital-basis vectors (orbitals and corresponding
      // potentials). If 0 -> no orbital-type grid renormalization.
      m_pGridRenormOrb;
   TArray<double>
      // electronic densities of frozen environment orbitals. Used in embedding calculations only.
      m_BackgroundDensities;
   bool
      AuxExpandXc;
   unsigned
      m_BaseAccXcFlags;

   double
      EnergyXc,
      // closed/open.
      fElecTotal[2], fElecTotalAnalytic[2];
   void PackAndAssignDensities(dfti::FDftiArgs &DftiArgs, FMatrixView RdmC, FMatrixView RdmO, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   void MakeDensityMatrices(FMatrixView &Rdm, FMatrixView &RdmO, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   double ComputeAuxiliaryElectronNumbers(FMatrixView AuxDen, FRawBasis const *pFitBasisRaw);
private:
   FFockComponentBuilderXc(FFockComponentBuilderXc const&); // not implemented
   void operator = (FFockComponentBuilderXc const &); // not implemented
};

// build df-coulomb and auxiliary expanded xc. Cache all integrals in memory.
struct FFockComponentBuilderDfCoulXcCached : public FFockComponentBuilderXc
{
   typedef FFockComponentBuilderXc
      FBase;
   FFockComponentBuilderDfCoulXcCached(FDfJkOptions const &JkOptions_, FXcOptions const &XcOptions_, FLog &Log_, FTimerSet *pTimers_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   ~FFockComponentBuilderDfCoulXcCached();
   void PrintEnergyContribs();
//    bool SwitchToRefineGrid(FDftGridParams const &NewGridParams);
//    void ComputeBackgroundDensities(FMatrixView CEnvOccC, FMatrixView CEnvOccO, FMemoryStack &Mem);
protected:
   FDfJkOptions
      JkOptions;
   FBasisSetPtr
      pFitBasis;
//    FRawBasis
//       *pFitBasisRaw;
   FMatrixView
      Jcd;
   size_t
      nFit,
      nAoTrSt;
   TArray<double>
      // format: Align([(nAo * (nAo+1))/2]) x nFit.
      // This version is for *small* cases and stores *all* integrals.
      Int3ixStorage;
   TArray<double>
      // last fitted auxiliary density coefficients; for use in gradient calculations.
      m_LastDensity;
   FMatrixView
      Int3ix;
   double
      EnergyCoulomb;
//    void Make1ixDensity(FMatrixView Jgamma, FMatrixView AuxDen, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   void Make1ixDensity1(FMatrixView Jgamma, FMatrixView AuxDen, FMatrixView RdmC, FMatrixView RdmO, size_t nDenVec, FMemoryStack &Mem);
//    double ComputeAuxiliaryElectronNumbers(FMatrixView AuxDen);
private:
   FFockComponentBuilderDfCoulXcCached(FFockComponentBuilderDfCoulXcCached const&); // not implemented
   void operator = (FFockComponentBuilderDfCoulXcCached const &); // not implemented
};




// // build df-coulomb and auxiliary expanded xc. Cache all integrals in memory.
// // TODO: rebuild this as base class of FFockComponentBuilderXc
// struct FFockComponentBuilderDfCoulXcCached : public FFockComponentBuilder
// {
//    FFockComponentBuilderDfCoulXcCached(FDfJkOptions const &JkOptions_, FXcOptions const &XcOptions_, FLog &Log_, FTimerSet *pTimers_);
//    void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
//    void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
//    void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
//    ~FFockComponentBuilderDfCoulXcCached();
//    void PrintEnergyContribs();
//    bool SwitchToRefineGrid(FDftGridParams const &NewGridParams);
//    void ComputeBackgroundDensities(FMatrixView CEnvOccC, FMatrixView CEnvOccO, FMemoryStack &Mem);
// protected:
// //    std::string
// //       XcFunctionalName;
// //    bool
// //       ApplyBandStructureEnergyFormula; // experimental option. see comments in CtRhfOptions.h
//    FDfJkOptions
//       JkOptions;
//    FXcOptions
//       XcOptions;
//    FAtomSet const
//       *pAtoms;
//    FBasisSetPtr
//       pFitBasis;
//    FBasisSet
//       *pOrbBasis;
//    FRawBasis
//       *pOrbBasisRaw,
//       *pFitBasisRaw;
//    FMatrixView
//       Jcd;
//    size_t
//       nAo,
//       nFit,
//       nAoTr,
//       nAoTrSt;
// //    FDftGridParams
// //       GridParams;
//    FDftGridPtr
//       pDftGrid;
// //    FXcFunctionalPtr
// //       pXcFn;
//    TArray<double>
//       // format: Align([(nAo * (nAo+1))/2]) x nFit.
//       // This version is for *small* cases and stores *all* integrals.
//       Int3ixStorage;
//    FMatrixView
//       Int3ix;
// //    FFockComponentBuilderDfCoulXcCachedImpl
// //       *p;
//    TArray<double>
//       m_LastDensity;
//    TArray<double>
//       // electronic densities of frozen environment orbitals. Used in embedding calculations only.
//       m_BackgroundDensities;
//    bool
//       AuxExpandXc;
//    unsigned
//       m_BaseAccXcFlags;
//
//    double
//       EnergyCoulomb, EnergyXc,
//       // closed/open.
//       fElecTotal[2], fElecTotalAnalytic[2];
// //    void Make1ixDensity(FMatrixView Jgamma, FMatrixView AuxDen, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
//    void Make1ixDensity1(FMatrixView Jgamma, FMatrixView AuxDen, FMatrixView RdmC, FMatrixView RdmO, size_t nDenVec, FMemoryStack &Mem);
//
//    void PackAndAssignDensities(dfti::FDftiArgs &DftiArgs, FMatrixView RdmC, FMatrixView RdmO, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
//    void MakeDensityMatrices(FMatrixView &Rdm, FMatrixView &RdmO, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
//    double ComputeAuxiliaryElectronNumbers(FMatrixView AuxDen);
// private:
//    FFockComponentBuilderDfCoulXcCached(FFockComponentBuilderDfCoulXcCached const&); // not implemented
//    void operator = (FFockComponentBuilderDfCoulXcCached const &); // not implemented
// };


class FD3DispersionContext;

struct FDispCorrBuilder : public FFockComponentBuilder
{
   FDispCorrBuilder(std::string const &DispCorrType, std::string const &XcFunctionalName, FLog &Log_, FTimerSet *pTimers_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   void AccGradient(FMatrixView Gradient, FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   ~FDispCorrBuilder();
   void PrintEnergyContribs();
protected:
   FD3DispersionContext
      *m_pD3DispersionContext;
   std::string
      m_DispCorrType,
      m_XcFunctionalName;
   FAtomSet const
      *pAtoms;
private:
   FDispCorrBuilder(FDispCorrBuilder const&); // not implemented
   void operator = (FDispCorrBuilder const &); // not implemented
};





// build IAO local exchange
struct FFockComponentBuilderDfLx : public FFockComponentBuilder
{
   FFockComponentBuilderDfLx(FDfJkOptions const &Options_, FMatrixView S1_, FMatrixView S1cd_, FLog &Log_);
   void Init(FWfDecl const &WfDecl_, FBasisSet *pOrbBasis_, FAtomSet const &Atoms_, FHfOptions const &Options_, FMemoryStack &Mem);
   void AccFock(FMatrixView &FockC, FMatrixView &FockO, FBasisSet *pBasis, FMatrixView const &COccC, FMatrixView const &COccO, uint Flags, FMemoryStack &Mem);
   ~FFockComponentBuilderDfLx();
public:
   FDfJkOptions
      Options;
   FAtomSet const
      *pAtoms;
   FBasisSetPtr
      pFitBasis,
      pMinBasis,
      pIaoBasis;
   FBasisSet
      *pOrbBasis;
   FMatrixView
      Jcd,  // cholesky decomposition of full JKFIT basis (may not be required, depending on the mode of use)
      S1,   // overlap matrix of orbital basis
      S1cd; // choleksy decomposition of orbital basis.
};

enum FIntMnfFlags {
   INTFLAGS_SkipTriangleFix = 0x01 // if set, do NOT fix up the upper triangle a>b of (ab|c) integrals.
};

double *FormIntMNF(ir::FRawShell const &ShF,
   ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasisA, FRawBasis const *pOrbBasisB, FRawBasis const *pFitBasis, FMemoryStack &Mem, FMatrixView ScrDen, double fThr, unsigned Flags=0);
double *FormIntMNF(ir::FRawShell const &ShF,
   ir::FIntegralKernel *pIntKernel, FRawBasis const *pOrbBasis, FRawBasis const *pFitBasis, FMemoryStack &Mem, FMatrixView ScrDen, double fThr);


} // namespace ct

#endif // CT_RHF_FOCKBUILD_H
