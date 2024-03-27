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

// this file exists for the sole purpose of avoiding circular header dependencies.
#ifndef CT_RHF_OPTIONS_H
#define CT_RHF_OPTIONS_H

#include <string>
#include "CtDftGrid.h" // for grid params.
#include "CtWfDecl.h"
#include "CxParse1.h"
#include "CxVec3.h"

namespace ct {
typedef TVector3<double>
   FVector3;

using mig::FDftGridParams;
using mig::FDftGrid;

// Changes a specific flag (or flag combination) in BitField, in the direction
// specified by `iDir`. Concretely:
// - if iDir == +1, FlagToChange is turned on in BitField,
// - if iDir == -1, FlagToChange is turned off in BitField,
// - if iDir == 0, FlagToChange is unchanged in BitField.
void ChangeFlag(unsigned int &BitField, unsigned int FlagToChange, int iDir);
void ChangeFlag(unsigned long &BitField, unsigned long FlagToChange, int iDir);


enum FXcAlgo {
   XCALGO_Regular, // standard computation of xc contributions, with full density.
   XCALGO_AuxiliaryExpand, // computation of xc contributions via auxiliary expansion of density computed in Coulomb.
   XCALGO_Gridless  // DFT without meshes and grids.
};

enum FJkAlgo {
   JKALGO_DfAuto, // choose one of DF-JK, DF-J (cached) or DF-J (direct), depending on functional and environment
   JKALGO_DfCoulombAndExchange, // DF-JK, as in regular density fitted RHF
   JKALGO_DfCoulombOnly, // DF-J, as in regular density fitted RKS with pure functionals
   JKALGO_DfCoulombOnlyCache3ix, // for small molecules: also DF-J, but attempt to cache *all* 3ix integrals in memory.
   JKALGO_DfCoulombDfLocalExchange, // full DF-J, and orbital localization based LDF-K.
   JKALGO_DfCoulombChainOfSpheresExchange, // full DF-J, and Neese's COSX---chain of spheres (grid based) local exchange.
   JKALGO_4ixDirect // regular non-density fitted Coulomb and Exchange (with 4-index ERIs etc); this one is integral-direct.
};

enum FHfComputeFlags {
   COMPUTE_Gradient = 0x01, // compute analytic gradient?
   COMPUTE_IboAnalysis = 0x02
//    COMPUTE_Multipoles = 0x04 // this one is for automatically computing charge/dipole/quadpole in SCF.
//    // ^- this one gets special treatment. Eventually I should make a separate
//    // real property interface to compute stuff like field gradients, ESP,
//    // distributed multipoles, interactions with charges, etc, and parse and
//    // attach the properties with actual legends to the atom set objects (like we
//    // treat energies and gradient at the moment)
};

enum FScfOrbType {
   ORBTYPE_Restricted, // spin-restricted (shared spatial orbitals for alpha- and beta spin),
   ORBTYPE_Unrestricted, // spin-unrestricted (alpha and beta orbitals with separate spatial parts)
   ORBTYPE_Generalized // generalized (allowing spin-orbitals with linear combinations between alpha
                       // and beta spin components; this corresponds to non-colinear spin of allowing complex amplitudes)
};

enum FPrintFlags {
   // in extra detail mode, probably should print the full basis data (the code
   // is already there, just switched off atm).
   SCFPRINT_GuessInfo = 0x00000001,
   SCFPRINT_OrbBasis = 0x00000002,
   SCFPRINT_FitBasis = 0x00000004,
   // @orbitals: could separate this into levels (e.g., valence-only, core+valence, all), using separate int to specify details
   // ...or projections onto the IAO basis or similar.
   SCFPRINT_Orbitals = 0x00000008,
   SCFPRINT_Ecp = 0x00000010, // print effective core potentials
   SCFPRINT_MethodWfDesc = 0x00000020,
   SCFPRINT_OrbitalEnergies = 0x00000040, // print valence orbital energies by default; in extra mode, also core level energies.
   // maybe data on timings? on densities? fock matrix? orbital occupation info inside iterations?
   SCFPRINT_Default = SCFPRINT_MethodWfDesc | SCFPRINT_OrbitalEnergies
};


struct FHfPrintOptions {
   explicit FHfPrintOptions(unsigned Flags_=SCFPRINT_Default) : m_FlagsPrintOn(Flags_), m_FlagsExtraDetail(0) {}
   bool operator [](unsigned Flag) const { return 0 != (m_FlagsPrintOn & Flag); }
   bool ExtraDetail(unsigned Flag) const { return 0 != (m_FlagsPrintOn & Flag) && 0 != (m_FlagsExtraDetail & Flag); }

   void SetArgs(std::string const &Args);
   void SetFlag(unsigned Flag, bool Value, bool Extra=false);
   void SetFlag(unsigned Flag, string_slice slValue);
protected:
   unsigned
      m_FlagsPrintOn,
      m_FlagsExtraDetail;
   // ^- maybe it would be better to just keep an array of signed chars around
   // instead of dealing with multi-level flags. Sounded like a good idea
   // originally, but even moving over the timing flags is not that easy as-is.
   // can't currently move the timing stuff because of that.
};




// MPT = Molecular property tasks.
// These ones define properties based on generic one-electron operators (and
// nuclear contributions) to evaluate.
struct FMpt1e : public FIntrusivePtrDest {
   enum FPropertyType {
      // this entry is invalid.
      PROPTYPE_Unknown,
      // evaluate total dipole moment of charge distribution (including nuclear
      // point charge contributions). This one gets special handling from the
      // other moment integrals because we need it for IR intensities, and it
      // is often more interesting.
      PROPTYPE_DipoleMoments,
      PROPTYPE_KineticEnergy,
      // For a given set of points \vec R_j in space, this operator evaluates either
      //
      // - (if Order==0): the total electron density (alpha+beta)
      // - (if Order==1): the total electron density (alpha+beta) + its x/y/z
      //                  gradients (four components total)
      // - (if Order==2): the total electron density (alpha+beta) +
      //                  x/y/z density gradients +
      //                  six unique hessian components
      //                  (order defined by FSlcX in IrAmrr.h; at the time of writing,
      //                  for IR's default configuration this is:
      //                  xx, yy, zz, xy, xz, yz
      // as PROPTYPE_ChargeDensity, but for the difference between alpha- and
      // beta-spin density
      PROPTYPE_ElectronDensity,
      // as PROPTYPE_ElectronDensity, but for the difference between alpha- and
      // beta-spin density
      PROPTYPE_SpinDensity,
      PROPTYPE_NuclearAttractionPotential,
      // For a given set of points \vec R_j in space, this operator evaluates either
      //
      // - (if Order==0): the electrostatic potential (ESP),
      // - (if Order==1): the components E_x, E_y, E_z of the electric field,
      // - (if Order==2): the electric field gradients: The derivatives
      //                  d(E_i) / d(R_j) for i,j \in \{x,y,z\}
      //                  of the electric field at the given points.
      //
      // The reported contributions include both electronic and nuclear contributions.
      // For the electronic part, this means evaluating integrals of the form
      //
      //   <\mu|(\partial {R_x})^e (\partial_{R_y})^f (\partial_{R_z})^g) 1/(r - R)|\nu>
      //
      // (Not sure what to do with COSMO or MD charges/lattice charges, etc. if
      // such things are implemented)
      //
      // NOTE: I don't think we need CtInt1e for that... point multipole thing
      // should be just fine.
      PROPTYPE_Electrostatic,
      // General Cartesian moments x^i y^j z^k of the charge distribution
      // around the origin. Includes contributions from nuclear point charges.
      //
      // Google Molpro manual "One-electron operators and expectation values (GEXPEC)"
      // for difference between Second Moments (SM) and Quadrupole Moments (QM).
      PROPTYPE_ChargeMoments_Cartesian
//       OPTYPE_MultipoleMoments,
      // note: we probably could do basic perturbative scalar relativistic corrections,
      // too, if needed... (darwin and mass-velocity at least)
   };
   FPropertyType
      m_PropertyType;

   int
      // if the operator can come in multiple different "orders" (e.g.,
      // Cartesian moments, derivative orders of field gradients), this value
      // specifies what this order is (or what order it is "up to", depending on
      // the operator in question).
      m_Order;
   double
      // results will be multiplied with this factor.
      m_Factor;

   // Note @center arguments:
   // - If the operator does not actually depend on a center
   //   coordinate (e.g., kinetic energy), these should be left empty.
   //   Behavior is undefined in this case (may be ignored or raise an error).
   // - If both a list of atom indices and of 3D center points is given, both
   //   will be processed, with the atom centers going first.
   // - We might want to add functions (or meta objects) for generating
   //   points on larger sets of geometries (e.g., regular grids or cubes)
   typedef std::vector<int>
      FIndexList;
   typedef std::vector<FVector3>
      FPointList;
   FIndexList
      // list of atom indices (0-based) on which's position to evaluate the operator
      // (e.g., for electric field gradient this would be of interest)
      m_iCenters;
   FPointList
      // list of points in 3D space on which to evaluate the operator.
      m_vCenters;

   int
      m_PrintLevel; // 0 = silent. (still goes to .npy property files, but not main output)

   FMpt1e();
   ~FMpt1e(); // only there to bind it. copying these guys should work just fine.

//    bool TryParse(string_slice slArgs, bool AssertThisWorks = true);
   void SetArgs(std::string const &Args_);
   void SetPropertyType(string_slice slType);

   // query if we are requested to *not* print this one to the main
   // output log (e.g., because it would be too long or annoying).
   bool IsSilent() const { return m_PrintLevel <= 0; }

   static bool FindPropertyType(string_slice slType, FPropertyType *pPropertyTypeOut, int *pOrderOut);
};
typedef ct::TIntrusivePtr<FMpt1e>
   FMpt1ePtr;
typedef ct::TIntrusivePtr<FMpt1e const>
   FMpt1eCptr;


// Defines properties of the system/molecule we would like to have evaluated by
// the program.
//
// - This currently encompasses only a very small number of very simple
//   properties, which can be evaluated directly from the wave function's
//   density matrix (and even then not all which might be needed; for example,
//   things like density at points/cubes, ESP on points or vdW surface,
//   interactions with finite electrostatic charges, etc. could all be quite
//   useful)
//
// - However, we might want to add qualitatively different things, too, at a
//   later point:
//
//   + Trivial things like rotational constants, molecule mass, etc.
//
//   + Qualitative analysis tasks (IAO analysis, IAO/IBO analysis, etc.; possibly
//     also in variants requiring more complex calculations (e.g., bond
//     polarization tensors)),
//
//   + Properties requiring calculations beyond single points, but still with
//     fixed schedules (e.g., finite field calculations for dipole moment
//     derivatives with respect to atomic positions, or static multipole
//     polarizabilities),
//
//   + Specialty applications for making simplified representations of the
//     molecule for MD and the like (e.g., distributed multipoles, point-like or
//     Gaussian, fits of the electron density or exchange potential), or for
//     theoretical purposes (e.g., evaluation of electronic pair distribution
//     functions or intracules or similar)
//
//   + ...or even calculations requiring complex schedule control (e.g.,
//     harmonic vibrational spectra and thermo corrections)
//
// - But none of those things are there yet. Next to the actual coding of the
//   properties themselves, this also requires more complex setup and design of
//   interfaces, calculation schedules, and means of storing the results.
//
//   ...but if we ever do these things, we need a place to store the commands
//   and options in, and that might be right here! Yay. So the structure carries
//   a more grand sounding name than it lives up to at the moment...
struct FMoleculePropertyTasks {
//    FMptMultipoles
//       Multipoles;
//    FMpt1ePtr
//       m_pMultipoles;
   typedef std::vector<FMpt1ePtr>
      F1eTaskList;
   F1eTaskList
      // list of (physical) properties which can be directly evaluated from a
      // given combination of atom set and molecular wave function.
      m_1eTaskList;

   void SetArgs(std::string const &Args_);
};


struct FHfGuessOptions
{
   enum FInitialGuessType {
      INITGUESS_AtomicDensity,
      INITGUESS_CoreH
      // ^-- probably should add also atomic potentials, and maybe some
      // Extended-Hückel like thing to combine with AtDen guess.
   };
   FInitialGuessType
      // type of initial guess to use when building an initial Fock matrix from
      // scratch (i.e., without externally supplied starting data).
      NewGuessType;

   enum FUsePreviousScfResultFlags {
      // if set, ignore FHfResult starting information provided to
      // FHfMethod::BuildInitialFock, and build the initial start from scratch
      USEINPUT_Never = 0x0000,
      // if set, allow using information calculated for the same molecule,
      // but at a different geometry (e.g., during geometry optimization)
      USEINPUT_AllowDifferentGeometry = 0x0001,
      // if set, allow using information calculated for the same molecule,
      // but with a different basis (e.g., during basis set extrapolation)
      USEINPUT_AllowDifferentBasis = 0x0002,
      // if set, allow using information calculated for the same molecule,
      // but for a different electronic structure method (e.g., starting RHF
      // from a converged RKS result).
      // WARNING: currently ignored.
      USEINPUT_AllowDifferentMethod = 0x0004,
      // if set, allow using information calculated for the same molecule,
      // but for a different electronic state (e.g., different charge or
      // spin)
      USEINPUT_AllowDifferentState = 0x0008,
      // if set, attempt to use the provided starting orbitals whenever this can
      // be done.
      USEINPUT_AlwaysIfPossible = USEINPUT_AllowDifferentGeometry | USEINPUT_AllowDifferentBasis | USEINPUT_AllowDifferentMethod | USEINPUT_AllowDifferentState
   };
   unsigned
      // controls when provided starting orbitals may be used to construct a
      // starting point for SCF (as opposed to assembling an SCF guess from
      // scratch via the method in InitialGuessType)
      UseInputWf;

   enum FGuessPropagationType {
      // if set, start a follow-up calculation using the last orbital set
      GUESSPROP_StraightOrbitals,
      // if set, employ a cgk-special for propagating initial guesses
      // via orbitals from one frame to another.
      GUESSPROP_ShOcc,
      // if set, start new calculations based on the last calculation's ROHF
      // Fock matrix diagonalized to make orbitals
      // (i.e., including ROHF block partitioning and possibly level shifts)
      GUESSPROP_Fock,
      // as GUESSPROP_Fock, but use a cgk-special for its representation
      // (store something similar to S^{-1/2} F S^{-1/2})
      GUESSPROP_SmhFock
   };
   FGuessPropagationType
      // controls which and how SCF restart data is stored.
      PropagationType;

   FHfGuessOptions();
   void SetArgs(string_slice const &Args_);
};



struct FHfOptions : public FIntrusivePtrDest {
   FHfPrintOptions
      Print;
   FHfGuessOptions
      InitialGuess;
   double
      ThrOrb, ThrDen; // threshold for orbital gradient (residual) and energy change, respectively
   unsigned
      MaxIt, // maximum number of iterations.
      MaxDiis; // maximum number of subspace vectors to keep around

//    double
//       ExactExchFactor; // prefactor for exact exchange. If 0, no exact exchange is computed.
   FXcAlgo
      // how to compute xc functional contributions, provided that an xc functional is set.
      XcAlgo;
   FJkAlgo
      // how to compute coulomb and exact exchange contributions.
      JkAlgo;
   FScfOrbType
      // restricted, unrestricted, generalized
      OrbType;
   uint32_t
      // decides on what we compute (COMPFLAG_*).
      ComputeFlags;
   double
      LevelShifts[2]; // [0]: for closed shells, [1]: for occupied open shells.
   double
      // additional closed-shell level shift to use in initial iterations
      // (half of that goes to active orbitals)
      // note: in order to actually dampen the initial iterations, this
      // value should be negative!
      InitialDampingLevelShift,
      // '1 - decay' is used as factor by which InitialDampingLevelShift is
      // multiplied each iteration 1 - decay should be <= 1.
      InitialDampingLevelShiftDecay;

   std::string
      XcFunctionalName,
      DispCorrType,
      // note: these for informative purposes only. actual basis sets taken
      // from basis assignment in pAtoms! By default they are not assigned.
      BasisName_Orbital,
      BasisName_JFit,
      BasisName_JkFit,
      DftGridDesc,
      ImplicitSolvationModelDesc;
   FDftGridParams
      DftGridParams,
      DftGridParamsRefine;
   bool
      // if set, switch to FinalGrid after convergence and make one more iteration.
      UseRefineGrid,
      // if set, there will be no exception if SCF fails to converge.
      // (a warning will still be emitted)
      IgnoreConvergenceError;
   bool
      // experimental option to reformulate KS total energy equation in such a way that
      // the sum of the orbital energies (multiplied by occupation numbers) is included.
      // *might* have different grid convergence, basis convergence, and SCF convergence
      // properties than the standard KS formulation (not yet sure if in good way or bad
      // way). Formally the reformulation is exact, but only supposed to be fully equivalent
      // to standard KS in the limit of full SCF convergence, infinite basis and possibly
      // infinite integration grid. Not sure yet.
      // Has no influence on Hartree-Fock; in the HF case, both total energy formulations
      // are mathematically equivalent.
      UseBandStructureEnergyFormula;
   enum FRohfAlgorithm {
      // make a single ROHF Fock matrix and diagonalize it to get the
      // next set of closed/active/external orbitals
      ROHFALGO_OneStep,
      // make a FockA/FockB. Diagonalize FockA to find the next space of (all!) occupied orbitals.
      // Diagonalize FockB in the subspace of FockA occupied orbitals to split the occupied
      // orbitals into closed/active.
      ROHFALGO_TwoStep
   };
   FRohfAlgorithm
      m_RohfAlgo;
   int
      iTimingLevel;
   FMoleculePropertyTasks
      PropertyTasks;

   // make a short string summarizing the options regarding what we are doing (method, parameters, etc.).
   std::string GetMethodDesc() const;

   explicit FHfOptions(double ThrOrb_=1e-5, double ThrDen_=1e-6, unsigned MaxIt_=120);
   void SetArgs(std::string const &CommandOptions);
   // note: this handles the refine grid stuff.
   void SetGridDesc(std::string const &NewGridDesc_);
   void SetXc(std::string const &XcName_);
   void SetDfxc(bool NewDfxc_);
   void SetDamping(std::string const &Args_);

   // returns 'true' unless JkAlgo is set to a 'pure' (Coulomb-only) density
   // functional method, or JkAlgo is set to JKALGO_DfAuto and XcFunctionalName
   // refers to a pure functional.
   //
   // Note: "exact" exchange means that some sort of k-matrix is computed (in
   // addition to j and possibly x), not that this is done without further
   // approximations---DFJK-RHF, LDF-RHF, RIJCOSX etc. all have "exact" exchange
   // featuring various levels of accuracy)
   bool UseExactExchange() const;
   // returns whether a xc functional calculation is necessary
   bool UseXc() const;
private:
   // next two are used for caching the results of UseExactExchange(). It may
   // have to instanciate the functional to make a call.
   std::string mutable
      m_LastXcFuncName;
   bool mutable
      m_LastXcFuncWasHybrid;
};

std::string GetMethodWfDesc(FHfOptions const *pHfOptions, FWfDecl const *pWfDecl);
char const *GetFitBasisType(FHfOptions const *pHfOptions, bool UpperCase);
// will return 'R', 'U', or 'G'.
char const *GetOrbTypePrefix(FScfOrbType OrbType);



} // namespace ct


#endif // CT_RHF_OPTIONS_H
