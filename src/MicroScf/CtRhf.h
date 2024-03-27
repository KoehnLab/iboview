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

#ifndef CT_RHF_H
#define CT_RHF_H

#include <string>

#include "CxTiming.h"
#include "CxIo.h"
#include "CxPodArray.h"

#include "CtBasisSet.h"
#include "CtAtomSet.h"
#include "CtMatrix.h"
// #include "CtDma.h"
#include "CtDftGrid.h" // for grid params.
#include "CtFockBuild.h" // for FFockComponentBuilderList only. Maybe maybe turn into something forward-declarable?
#include "CtRhfOptions.h"

namespace ct {


// Stores stuff which came out of a calculation. Can also be used for initial
// guesses of later calculations.
struct FHfResult : public FIntrusivePtrDest
{
   FHeapMatrix
      Orb, // nAo x nAo array.
      Fock, // nAo x nAo array.
      ExchO, // nAo x nAo array.
      Occ, // nAo x 1 array; occupation numbers (not occupied orbitals)
      Ew; // nAo x 1 array: orbital eigenvalues
   FHeapMatrix
      // theoretically this can be compued from Orb and Occ.
      // Keep anyway? Needs lots of memory.
      Rdm, RdmO;
   FBasisSetPtr
      pOrbBasis;
   double
      Energy;
   FAtomSetCptr
      // hm... maybe we should make an actual copy of this? This wouldn't be
      // good if atoms are actually moved from the outside (e.g., during optg or
      // similar). In that case this pAtoms would not be equivalent to the one
      // we made the wf for anymore.
      pAtoms;
   size_t
      nIt; // actual number of iterations used.
   bool
      Converged;
   FHfOptions
      // HF options and WF decl on which these results are based.
      HfOptions;
   FWfDecl
      WfDecl;

   FHfResult();
   ~FHfResult();

   // The following sets of matrices may be used for starting guess propagation
   // if the corresponding options for InitialGuessPropagationType are set.
   // If not, these entries cannot be assumed to be present.
   FHeapMatrix
      // nAo x nOccOld array (or nAo x nAo);
      ShOcc,
      // a RohfFock matrix (including ROHF projections and possibly level
      // shifts) which may diagonalized to yield orbitals.
      RohfFock,
      // an orthonormal variant of RohfFock (something similar to S^{-1/2} RohfFock S^{-1/2})
      SmhRohfFock,
      // an orthonormal variant of ExchO.
      SmhExchO;
   bool
      RohfFockIsLevelShifted;
   FBasicPerm
      // may be used in conjuction with initial guess propagation via ShOcc
      // and/or SmhRohfFock. Otherwise not used.
      ShTrafoBasisPerm;
};
typedef TIntrusivePtr<FHfResult>
   FHfResultPtr;


struct FHfMethod
{
   FHfOptions const
      &m_Options;
   FWfDecl
      m_WfDecl;

   enum FHfRunOptions {
      // if set, no FHfResult object will be assembled to retain a copy of
      // the final orbitals, Fock matrix, etc.
      EXEC_DiscardResults = 0x0001,
      // if set, the final memory stack will not be reset. May be used to
      // temporarily retain some information inside the Fock builders.
      EXEC_RetainMemory = 0x0002,
      EXEC_Default = 0x0000
   };

   // note: pAtoms_ will be linked to output pResult_, if provided. I.e., the pResult_ entry will keep a reference.
   FHfMethod(FLog &Log_, FWfDecl const &WfDecl_, FAtomSetCptr const pAtoms_, FHfOptions const &Options_, FHfResult *pStartingGuess_, unsigned ExecFlags, FMemoryStack &Mem);

   FHfResultPtr pResult() { return m_pResult; } // intentionally not const (one may modify the result)
protected:
   FLog
      *m_pLog;
   FHfResultPtr
      m_pResult;
   unsigned
      m_ExecFlags;

   FAtomSet const
      // Geometry & environment we are calculating:
      // required to instantiate basis sets and 1e integrals
      *m_pAtoms;
//    FAtomSetCptr const
//       pAtoms;
   FBasisSetPtr
      m_pOrbBasis,
//       pFitBasis,
      m_pMinBasis;
   size_t
      m_FockSize; // size of fock matrix block to extrapolate.
   FMatrixView
      m_Ew,      // fock matrix eigenvalues.
//       Overlap, // S := <mu|nu> of orbital basis
      m_Smh,     // S^-{1/2}. WARNING: not built anymore.
      m_Scd,     // Cholesky decomposition of S
      m_CoreH,   // 1e Hamiltonian
      m_Fock,    // Fock matrix (regular total Fock matrix)
      m_ExchO;   // open-shell exchange.
   FHfOptions::FRohfAlgorithm
      m_RohfAlgo;
   FMatrixView
      // ROHF total Fock matrix which was last used to assemble an orbital set
      // by diagonalization (includes funky ROHF mix & match block partitioning,
      // and possibly level shifts). Aliased to 'Fock' in the closed-shell case.
      m_RohfFock,
      // Alpha- and Beta-Fock matrices. Only used if m_RohfAlgo == ROHFALGO_TwoStep, and in
      // this case, RohfFock, FockA, and FockB are allocated in a single continuous block (for DIIS).
      m_FockA,
      m_FockB;
   FMatrixView
      m_COrb;     // full (nAo,nAo)-shape orbital coefficient matrix
   FTimerSet
      m_Timers;
   FTimer
      m_tTotal;
   double
      m_Energy; // total energy of result.
   bool
      // set once we have a set of orbitals which is reasonably consistent with
      // the current Fock matrix, as computed with a RHF coupling pattern. This is
      // not the case if (a) we are in the 1st iteration, and (b) we used anything
      // except for a previous orbital set as initial guess.
      m_HaveOrbitals;
   bool
      m_FockMatrixIsLevelShifted;
   size_t m_nOrbBasisSize; // was: nao
   size_t _GetNumBf() { return m_nOrbBasisSize; }

   FFockComponentBuilderList
      m_FockComponentBuilders;

   void Init(FHfResult *pStartingGuess, FMemoryStack &Mem);
   bool CanUseStartingGuess(FHfResult *pStartingGuess);
   void BuildInitialFock(FHfResult *pStartingGuess, FMemoryStack &Mem);
   void Run(FMemoryStack &Mem);
   void UpdateOrbitals(FMatrixView Orb, FMatrixView Ew,/*FMatrixView &COccC, FMatrixView &COccO,*/ FMatrixView &DenC, FMatrixView &DenO,
      /*FMatrixView const Fock,*/ /*FMatrixView const &ExchO,*/ /*FMatrixView &T1,*/ FMemoryStack &Mem);
   void RebuildDensityMatrices(FMatrixView &DenC, FMatrixView &DenO, FMemoryStack &Mem);
   enum FExtractOrbitalsFlags {
      ORBEXTRACT_AllowAliasing = 0x0001,
      ORBEXTRACT_AbsorbOcc = 0x0002
   };
   // Extracts occupied subset of orbitals from COrbIn (e.g., this->Orb) according to this->WfDecl.
   // - if ORBEXTRACT_AllowAliasing is set, COccC and COccO may be aliased into COrbIn.
   //   Otherwise, COccC and COccO get allocated on Mem.
   // - if ORBEXTRACT_AbsorbOcc is set, COccC and COccO will be returned with
   //   absorbed occupation numbers; that is, they get scaled by the square
   //   roots of their occupation numbers (or in more general cases later, may
   //   get transformed by the matrix square root of the density matrix).
   void ExtractOrbitals(FMatrixView &COccC, FMatrixView &COccO, FMatrixView const &COrbIn, unsigned Flags, FMemoryStack &Mem);
   void MakeInitialGuess(FMemoryStack &Mem);
   // output: 3 x pAtoms->size() array.
   void MakeGradient(double *pOut, FMemoryStack &Mem);
   void PrintGradient(double *pGrad);
   // note: fAbsorbedOccC/O denotes the occupation numbers which have been absorbed into the respective
   // orbital matrices. Set to 1.0 if no occupation numbers have been absorbed.
   void ApplyLevelShifts(double fDirection, FMatrixView COccC, double fAbsorbedOccC, FMatrixView COccO, double fAbsorbedOccO, double fExtraShiftC, double fExtraShiftO, FMemoryStack &Mem);
   // this one extracts COccC/COccO automatically from this->Orb; otherwise does the same as ApplyLevelShifts.
   void ApplyLevelShifts1(double ExtraShift, FMemoryStack &Mem);

   void BuildRhfFock(FMatrixView const COrb_, bool InRefinementIt, FMemoryStack &Mem);
//    void RemoveLevelShifts(FMatrixView COccC, FMatrixView COccO, FMemoryStack &Mem);
   void UpdateOrbitalEnergies(FMemoryStack &Mem);
//    void EvalDma(FDmaDescPtr &pDmaResult, FMemoryStack &Mem);

   void PrintOrbitalEnergies(FMatrixView Ew, std::string const &Desc, unsigned Flags = ORBEIG_AlphaOn | ORBEIG_BetaOn);

   void EvalProperties1e(FMemoryStack &Mem); // this one is the driver routine
//    FAtomSetPropertyPtr EvalProperty1e(FMpt1e const *pMpt1e, FMemoryStack &Mem);

   void AssembleEnergy(FMatrixView DenC, FMatrixView DenO, FMemoryStack &Mem);
   void ComputeOrbitalGradient(FMatrixView OrbGrad, double *pfOrbGrad, FMatrixView DenC, FMatrixView Fock, FMemoryStack &Mem);
   void ComputeExtraEnergies(double &CoreEnergy, double &BandEnergy, FMatrixView DenC, FMatrixView DenO, FMemoryStack &Mem);
   void StoreScfResultData1(FMatrixView DenC, FMatrixView DenO, FMemoryStack &Mem);

public:
   enum FOrbEwPrintFlags {
      ORBEIG_AlphaOn = 0x01,
      ORBEIG_BetaOn = 0x02
   };
protected:
   void GetMoSubspaceOffs(size_t &ibeg, size_t &iend, char Space);
   void ProjectG(FMatrixView ProjectedG, char const *pSpaces, FMatrixView InputG, FMatrixView C, FMatrixView SC, double Factor, FMemoryStack &Mem);
   // C is the orbital matrix, SC is the S*C matrix.
   // If WfDecl.nOpen() == 0, it is assumed that FockRohf and FockC are aliased (closed-shell Fock matrix)
   void MakeRohfFockAo(FMatrixView FockRohf, FMatrixView FockA, FMatrixView FockB, FMatrixView Fock, FMatrixView ExchO, FMatrixView C, FMatrixView SC, double LevelShiftC, double LevelShiftO, FMemoryStack &Mem);
   void MakeRohfFockAo1(FMatrixView FockRohf, FMatrixView FockA, FMatrixView FockB, FMatrixView Fock, FMatrixView ExchO, FMatrixView C, FMemoryStack &Mem);
   void ApplyLevelShiftToSide(FMatrixView Fock0, FMatrixView Fock1, FMatrixView Fock2, FMatrixView COcc, double fOccNumber, double LevelShift01, double LevelShift2, double fDirection, FMemoryStack &Mem);
//    void ApplyLevelShiftsToInitialGuess(double f0, FMemoryStack &Mem);
   enum FTimingFlags {
      TIMING_Reset = 0x01,
      TIMING_NewLineAfter = 0x02,
      TIMING_NewLineBefore = 0x04,
      TIMING_NewLineBeforeAndAfter = TIMING_NewLineBefore | TIMING_NewLineAfter
   };
   void PrintTimingsAtLevel(int iTargetLevel, unsigned Flags);
};

} // namespace ct

#endif // CT_RHF_H
