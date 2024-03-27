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

#ifndef IR_H
#define IR_H

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

#if defined(_DEBUG) && !defined(IR_INCLUDE_PRINT_FUNCTIONS)
   #define IR_INCLUDE_PRINT_FUNCTIONS
#endif

#ifdef IR_INCLUDE_PRINT_FUNCTIONS
   #include <ostream> // for debug printing of shell information.
   #include <stdio.h> // for printf
#endif

#include <stddef.h> // for size_t

#ifdef MOLPRO
   #define INCLUDE_OPTIONALS
#endif

#ifdef IR_KERNEL_PTRS
   #include "CxTypes.h" // for intrusive ptr base class of integral kernels
#endif // IR_KERNEL_PTRS

#include "CxDefs.h"
#include "CxMemoryStack.h"
#include "IrBoysFn.h"

namespace ir {
   /// A low level shell data structure the integral drivers work with
   /// (for simplifying talking to existing Fortran and C++ programs!)
   struct FRawShell
   {
      unsigned
         /// angular momentum, number of primitive exponents, number of contractions
         l, nExp, nCo;
      double const
         /// pExp[iExp], iExp = 0..nExp-1: Exponent of primitive #i.
         *pExp,
         /// pCo[iExp + nExp*iCo]: coefficient of primitive #iExp in contraction #iCo
         ///
         /// Note: Contraction coefficients are stored with respect to raw,
         /// unnormalized Gaussians. Use the RawGaussNorm() function to convert
         /// between contraction coefficients in terms of normalized Gaussians
         /// (as used in libraries, for example), and raw coefficients (as used
         /// in most programs)
         *pCo,
         /// pCenter[i]: components 0..2: x,y,z position in space.
         *vCen;
      double const
         /// Can be 0. If provided: Screening information.
         /// [0]: largest range of any contracted function in the shell,
         /// [1+iExp]: range of primitive #iExp
         /// [1+nExp+iCo]: range of contraction #iCo
         *pRange;

      FRawShell() {}
      FRawShell(unsigned l_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_,
                double const *pvCen_, double const *pRange_=0)
         : l(l_), nExp(nExp_), nCo(nCo_), pExp(pExp_), pCo(pCo_),
           vCen(pvCen_), pRange(pRange_)
      {}

      /// return number of spherical components per contraction
      inline unsigned nSh() const { return (2*l+1); }
      /// return number of functions represented by the shell
      inline unsigned nFn() const { return (2*l+1) * nCo; }
      /// return contraction coefficient of primitive iExp in contraction iCo.
      inline double fCo(unsigned iExp, unsigned iCo) const { assert(iCo < nCo && iExp < nExp); return pCo[iExp + nExp*iCo]; }

      inline double MaxCoRange() const { assert(pRange); return pRange[0]; }
      inline double CoRange(unsigned iCo) const { assert(pRange && iCo < nCo); return pRange[1+nExp+iCo]; }
      inline double ExpRange(unsigned iExp) const { assert(pRange && iExp < nExp); return pRange[1+iExp]; }

      /// return number of primitive Gaussian functions of which the shell is composed
      /// (note: this is mainly for display/statistics purposes and for IR
      /// internal use; nFn and strides of input/output quantities to IR drivers
      /// exclusively refers to contracted functions (as nFn computes), not to
      /// nPrimFn!)
      inline unsigned nPrimFn() const { return (2*l+1) * nExp; }
   };

   /// Norm of the radial functions (incl. r^l) of a primitive Gaussian of
   /// angular momentum 'l'. Divide contraction coefficients CoMatrixNormalized
   /// from integral library format (e.g., libmol files, or data from basis set
   /// exchange) by this as follows:
   ///
   ///   for (unsigned iExp = 0; iExp < nExp; ++ iExp) {
   ///      double
   ///         fInvNorm = 1./RawGaussNorm(pExp[iExp], AngMom);
   ///      for (unsigned iCo = 0; iCo < nCo; ++ iCo)
   ///         CoMatrixRaw[iExp + iCo*nExp] = fInvNorm * CoMatrixNormalized[iExp + iCo*nExp];
   ///   }
   ///
   /// to convert them contraction specifications in terms of raw Gaussian
   /// primitives (CoMatrixRaw), which are used as the internal format by IR
   /// (e.g., in the FRawShell data structure).
   double RawGaussNorm(double fExp, unsigned l);
   /// returns 1/RawGaussNorm(fExp,l).
   double InvRawRawGaussNorm(double fExp, unsigned l);


   /// integral kernels represent scalar functions K(r_{12}), to be used in matrix
   /// elements like \int A(r_1) K(r_1 - r_2) C(r_2) d^3r_1 d^3r_2.
   ///
   /// Note: full details on the description of what these kernels do and how
   /// the concrete variants are derived, and how to develop new ones for different
   /// integral kernels K(r_{12}), if needed, are provided in the following article,
   /// and Appendices E and F of its supporting information:
   ///
   /// [3] Mieke Peels, Gerald Knizia - "Fast evaluation of two-center integrals over
   /// Gaussian charge distributions and Gaussian orbitals with general interaction
   /// kernels", J. Chem. Theory Comput. 2020
   /// https://doi.org/10.1021/acs.jctc.9b01296
   struct FIntegralKernel
#ifdef IR_KERNEL_PTRS
      : public ct::FIntrusivePtrDest
#endif // IR_KERNEL_PTRS
   {
      virtual double MaxRange() const = 0;
      /// evaluate the Ahlrichs-style scalar kernel integrals
      ///
      ///     Gm(rho,T) = (-D/D[T])^m G0(rho,T)
      ///
      /// for all m = 0...MaxM (MaxM included), where T = rho |P-Q|^2.
      /// All output data is multiplied by Factor.
      virtual void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const = 0;

      /// compute the modified kernel integrals
      ///
      ///     \tilde Gm(rho,T) = (+2 rho D/D[T])^m G0(rho,T).
      ///
      /// These are used in two-center integration routines [3].
      ///
      /// The default implemenation just computes standard Gm(rho,T) kernel
      /// integrals and then transforms them, and there is no need to modify it.
      /// But in some integral types slight performance improvements could be
      /// achieved by a direct re-implementation (here done for overlap and Coulomb).
      virtual void EvalTildeGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const;
      virtual ~FIntegralKernel() = 0;

      /// return a short description of what this kernel is/does
      virtual char const *Desc() const = 0;
   };
#ifdef IR_KERNEL_PTRS
   typedef ct::TIntrusivePtr<FIntegralKernel>
      FIntegralKernelPtr;
   typedef ct::TIntrusivePtr<FIntegralKernel const>
      FIntegralKernelCptr;
#endif // IR_KERNEL_PTRS

   /// 2e interaction kernel for \f[ K(r_{12}) = \delta(r_{12}) \f]
   struct FOverlapKernel : public FIntegralKernel
   {
      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      void EvalTildeGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FOverlapKernel();
   };

   /// 2e interaction kernel for \f[ K(r_{12}) = 1/r_{12} \f]
   struct FCoulombKernel : public FIntegralKernel
   {
      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FCoulombKernel();
   };


   /// 2e interaction kernel for \f[ K(r_{12}) = r_{12}^2 \f]
   ///
   /// Note:
   /// - I have derived formulas for general \f[ K(r_{12}) = r_{12}^{2k} exp(-\omega r_{12}^2), \f]
   ///   of which this particular kernel here would be a special case.
   /// - They are not implemented here, but if you need them: feel free to write a note to cgk
   struct FRsqKernel : public FIntegralKernel
   {
      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FRsqKernel();
   };


// #ifdef INCLUDE_OPTIONALS
   /// 2e interaction kernel for \f[ K(r_{12}) = \frac{\operatorname{erf}(r/r_c)}{r} \f]
   /// WARNING: not tested.
   struct FErfCoulombKernel : public FIntegralKernel
   {
      explicit FErfCoulombKernel(double Rc);
      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FErfCoulombKernel();
   protected:
      double
         m_Omega,
         m_Rc; // Omega == (1/Rc^2)
   };

   /// 2e interaction kernel for \f[ K(r_{12}) = \frac{\operatorname{erfc}(r/r_c)}{r} \f]
   /// WARNING: not tested.
   struct FErfcCoulombKernel : public FIntegralKernel
   {
      explicit FErfcCoulombKernel(double Rc);
      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FErfcCoulombKernel();
   protected:
      double
         m_Omega,
         m_Rc; // Omega == (1/Rc^2)
   };

   /// 2e interaction kernel given as sum of some other kernels, with
   /// prefactors.
   /// WARNING: This object keeps referencing the original kernels
   /// and factors!!
   struct FMetaKernel : public FIntegralKernel
   {
      FMetaKernel(FIntegralKernel const **pKernels, double const *pFactors, unsigned nKernels);
      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
   protected:
      FIntegralKernel const
         **m_pKernels;
      double const
         *m_pFactors;
      unsigned
         m_nKernels;
      double
         m_Range;
   };
// #endif // INCLUDE_OPTIONALS

// #ifdef INCLUDE_OPTIONALS
   /// Describes a single general geminal, expressed as a contraction over
   /// primitive Gaussian geminals:
   ///
   /// \f[  K(r_{12}) := F(r_{12}) = \sum{i=0}^{n-1} c_i \exp\bigl(-\zeta_i * r_{12}^2\bigr) \f]
   ///
   /// Evaluating integrals over these is reasonably efficient. This can be used to
   /// fit arbitrary scalar integral kernels.
   ///
   /// Note:
   /// - I have derived formulas for general \f[ K(r_{12}) = r_{12}^{2k} exp(-\omega r_{12}^2), \f]
   ///   of which this particular kernel here would be a special case (if Gaussian contractions are added).
   /// - They are not implemented here, but if you need them: feel free to write a note to cgk
   struct FGaussGeminal
   {
      double const
         /// pointer to externally provided array of exponents \zeta_i
         *IR_RP pExp,
         /// pointer to externally provided array of coefficients c_i
         *IR_RP pCo;
      size_t
         /// number of contractions (i.e., number of (\zeta_i, c_i) summed over)
         nExp;
      /// Create a contracted Gaussian geminal kernel, representing the two-
      /// electron interaction given by the function
      ///
      /// \f[     K(r_{12}) = \sum_{i}^{\mathrm{nExp}} c_i * e^{-\zeta_i r_{12}^2}  \f]
      ///
      /// Arguments:
      /// - nExp: number of exponents in the contraction (use nExp=1 for a primitive geminal)
      /// - pCo_: pCo_[0], ..., pCo_[nExp-1]: contraction coefficients c_i
      /// - pExp_: pExp_[0], ..., pExp_[nExp-1]: exponents \zeta_i
      /// - MaxRange: if given, integral routines may approximate K(r_{12}) as
      ///   zero for inter-electron ranges r_{12} >= MaxRange.
      ///   (If not given, MaxRange will be estimated based on the lowest exponent
      ///   in the contraction)
      ///
      /// WARNING:
      /// - This object DOES NOT COPY the input data, but keeps on referencing
      ///   the array pointers you provide in here!
      ///
      /// - pCo_ and pExp_ therefore must stay valid as long as this object is
      ///   used. Do not pass pointers to temporaries with shorter lifetime
      ///   than this geminal object!
      FGaussGeminal(double const *IR_RP pExp_, double const *IR_RP pCo_, size_t nExp_, double MaxRange_ = -1.) : pExp(pExp_), pCo(pCo_), nExp(nExp_), m_MaxRange(MaxRange_) { Finalize(); }

      double MaxRange() const { return m_MaxRange; };
   private:
      double m_MaxRange;
      void Finalize();
   };

   /// 2e interaction with generic scalar kernel provided as a contracted Gaussian geminal: \f[ F(r_{12}) = \sum{i=0}^{n-1} c_i \exp\bigl(-\zeta_i * r_{12}^2\bigr)\f]
   struct FGaussKernel : public FIntegralKernel
   {
      explicit FGaussKernel(FGaussGeminal const *pGaussGeminal_) : m_pGaussGeminal(pGaussGeminal_) {};

      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FGaussKernel();
   private:
      FGaussGeminal const
         *m_pGaussGeminal;
   };

   /// 2e interaction kernel \f[ F(r_{12})/r_{12} \f]
   struct FGaussCoulombKernel : public FIntegralKernel
   {
      explicit FGaussCoulombKernel(FGaussGeminal const *pGaussGeminal_) : m_pGaussGeminal(pGaussGeminal_) {};

      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FGaussCoulombKernel();
   private:
      FGaussGeminal const
         *m_pGaussGeminal;
   };

   /// 2e interaction kernel \f[ \bigl[\nabla F(r_{12})\bigr]^2 \f]
   struct FGaussKineticKernel : public FIntegralKernel
   {
      explicit FGaussKineticKernel(FGaussGeminal const *pGaussGeminal_) : m_pGaussGeminal(pGaussGeminal_) {};

      double MaxRange() const; // override
      char const *Desc() const; // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FGaussKineticKernel();
   private:
      FGaussGeminal const
         *m_pGaussGeminal;
   };
// #endif // INCLUDE_OPTIONALS

   /// the following refer to global overlap and Coulomb kernel objects.
   /// We keep around global versions because these kernels have no parameters,
   /// and there is no point in having more than one (but you can just make
   /// a new one if you want, anyway. Wouldn't hurt or help)
   extern FOverlapKernel g_IrOverlapKernel;
   extern FCoulombKernel g_IrCoulombKernel;

   /// Evaluate 2-electron, 3-center integrals:
   /// \f[  \pOut[a,b,c] := \prefac (ab|K_{12}|c) \f]
   void EvalInt2e3c(double *pOut, size_t *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   /// Evaluate |c)-contraction of 2-electron, 3-center integrals (ab|kernel|c):
   /// \f[  \pOut[a,b,V] := \prefac \sum_c (ab|K_{12}|c) \den[c,V]  \f]
   /// with \f$\den[c,V]\f$ `= pIn[c + StrideIn*V]`
   void EvalInt2e3c_ContractC(double *pOut, size_t *Strides, double const *pIn, size_t StrideIn, size_t nVecIn,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   /// Evaluate (ab|-contraction of 2-electron, 3-center integrals (ab|kernel|c):
   /// \f[  \pOut[c] := \prefac \sum_{ab} (ab|K_{12}|c) \den[a,b]  \f]
   /// with \f$\den[a,b]\f$ `= pDenAB[a*sa + b*sb]`
   void EvalInt2e3c_ContractAB(double *pOut, FRawShell const *pA, FRawShell const *pB,
      double const *pDenAB, size_t sa, size_t sb, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);

   /// Evaluate 2-electron, 4-center integrals:
   /// \f[  \pOut[a,b,c,d] := \prefac (ab|K_{12}|cd) \f]
   void EvalInt2e4c(double *pOut, size_t *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pC, FRawShell const *pD,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);

   /// Evaluate or accumulate 2-electron, 2-center integrals:
   /// \f[   \pOut[a,c] := \prevdata + \prefac (a|K_{12}|c) \f]
   /// If `Add==true`, data at output location is incremented by new integrals, otherwise replaced.
   void EvalInt2e2c(double *pOut, size_t StrideA, size_t StrideC, FRawShell const *pA, FRawShell const *pC,
      double Prefactor, bool Add, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   /// Evaluate or accumulate 2-electron, 2-center integrals:
   /// \f[   \pOut[a,c] := \prevdata + \prefac (a|K_{12} \Delta^n|c)   \f]
   /// If `Add==true`, data at output location is incremented by new integrals, otherwise replaced.
   void EvalInt2e2c_LaplaceC(double *pOut, size_t StrideA, size_t StrideC, FRawShell const *pA, FRawShell const *pC,
      double Prefactor, bool Add, unsigned LaplaceOrder, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   /// Evaluate or accumulate 1st derivative (with respect to atomic coordinates) of 2-electron, 2-center integrals:
   /// \f[   \var{pOutAxyz}[a,c,\kappa] := \prevdata + \frac{\partial}{\partial A_\kappa} (a|K_{12} \Delta^n|c)\f]
   /// \f[   \var{pOutCxyz}[a,c,\kappa] := \prevdata + \frac{\partial}{\partial C_\kappa} (a|K_{12} \Delta^n|c)\f]
   /// with $\kappa\in\{0,1,2\}$ indexing the Cartesian 1st derivative components.
   /// If `Add==true`, data at output location is incremented by new integrals, otherwise replaced.
   void EvalInt2e2c1d_LaplaceC(double *pOutAxyz, double *pOutCxyz, size_t StrideA, size_t StrideC, size_t StrideDeriv, FRawShell const *pA, FRawShell const *pC,
      double Prefactor, bool Add, unsigned LaplaceOrder, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);

   /// Compute 1st derivative integrals
   ///
   ///       d/d[Ru]  (a b|kernel|c)
   ///
   /// Where Ru (R=A,B, u=x,y,z) denotes the positions of the atomic centers on
   /// which the basis functions pA and pB sit.
   ///
   /// Notes:
   /// - The output of the derivative components is organized as
   ///
   ///       pOut[... + iComp*Strides[3]]
   ///
   ///   with iComp = 0...5 in order: d/d[Ax], d/d[Ay], d/d[Az], d/d[Bx], d/d[By], d/d[Bz].
   ///
   /// - This function *DOES NOT* explicitly compute the derivatives with respect
   ///   to C. These can be obtained from the A and B derivatives via
   ///   translational invariance:
   ///
   ///       (d/d[Ax] + d/d[Bx] + d/d[Cx]) I = 0.
   ///
   /// - TODO: consider extending this function to indirect addressing of output
   ///   components by passing eiter iCenterA/iCenterB as input or a double
   ///   **pOut denoting the output addresses directly. Currently I always
   ///   seem to be storing stuff to buffers and munching on it afterwards.
   void EvalInt2e3c1d(double *pOut, size_t *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   void EvalInt2e3c2d(double *pOut, size_t *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   /// evaluate laplace op commutator integral:  ([Nabla a]b|c) - (a [Nabla b]|c)
   void EvalInt2e3c_kcomm(double *pOut, size_t *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   void EvalInt2e3c1d_ContractC(double *pOut, size_t *Strides, double const *pIn, size_t StrideIn, size_t nVecIn,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, size_t nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);

   /// Evaluate the nGridPts x #bf matrix of basis functions (all must be placed
   /// at the same center vBfPos) on the given set of grid points. Output:
   ///
   ///    pOut[iPt + nCompStride * iDiff + nBfStride * iMap] = (value).
   ///        iDiff: derivative component.
   ///        iMap: non-screened basis function index (i.e., grid values
   ///              for iBf = pMap[iMap], where iBf is a basis function
   ///              index relative to the first bf in pFirstBf)
   ///
   /// Input grid point positions in 3D space are read from
   ///    pGridPos[iGridPt * GridStride + ixyz]
   ///        iGridPt: index of grid point
   ///        ixyz: indexes Cartesian components x (ixyz=0), y (ixyz=1), z (ixyz=2).
   ///
   /// On output, pMap/nMap will contain the list of indices of basis functions
   /// between pFirstBf and pLastBf, which have *NOT* been completely screened
   /// away on all grid points.
   ///
   /// The algorithm implemented here is described in:
   ///    [4] Alyssa Bienvenu, Gerald Knizia - "Efficient Treatment of Local
   ///       Meta-generalized Gradient Density Functionals via Auxiliary Density
   ///       Expansion: The Density Fitting J + X Approximation"
   ///       J. Chem. Theory Comput. 2018, 14, 3, 1297-1303
   ///       https://doi.org/10.1021/acs.jctc.7b01083
   ///
   /// TODO: make a function with similar interface which will work if functions
   /// are not all on the same center. pBfPos then goes away. Screening interface
   /// could also be adjusted to IR-style.
   void EvalShellGroupOnGrid(double *pOut, size_t nCompStride, size_t nBfStride,
      size_t *pMap, size_t &nMap, size_t iFnBase,
      double const *pGridPos, size_t GridStride, size_t nGridPts,
      FRawShell const *pFirstBf, FRawShell const *pLastBf, double const *pBfPos,
      unsigned DerivOrder, bool MakeLaplaceToo, double ThrOrb, double LogThrOrb, ct::FMemoryStack &Mem);

   /// For a given 3x3 unitary rotation matrix R[3][3], construct the matrix T[lc,l'c'] which makes
   ///    S[l,c](R r) = \sum_{c',l'} T[lc,l'c'] S[l',c'](r)
   /// for all 3-vectors 'r' (i.e., compute the transformation between solid harmonic components
   /// induced by rotations in real space).
   ///
   /// Notes:
   ///  - The matrix elements T[lc,l'c'] for l != l' vanish.
   ///  - T[lc,l'c'] is stored at pOut[iSlcX(l,c) + nSlmX(MaxL) * iSlcX(l',c')].
   ///    nStride is set to nSlmX(MaxL); pOut is allocated on Mem by this routine.
   ///  - Order of trafo not checked. If it doesn't work, use transpose(T) instead of T.
   void EvalSlcXRotationTrafo(double *&pOut, size_t &nStride, unsigned MaxL, double const *R, ct::FMemoryStack &Mem);

   /// Returns an integer bitfield describing the behavior of the 'c'th angular
   /// component basis function (Slc) in a shell of angular momentum 'l', upon
   /// applying mirror symmetries of space.
   ///
   /// Background:
   /// - Upon application of mirror symmetries of space (i.e., mirror at x=0
   ///   plane, y=0 plane, ...), all basis functions we use have the properties
   ///
   ///     Slc(x,y,z) = fx * Slc(-x,y,z)
   ///     Slc(x,y,z) = fy * Slc(x,-y,z)
   ///     Slc(x,y,z) = fz * Slc(x,y,-z)
   ///
   ///   where fx,fy,fz \in {+1,-1} (i.e., all functions either are invariant
   ///   to x/y/z mirror symmetry operations, or they flip their sign).
   ///
   /// - This function describes the factors fx,fy,fz. Concretely, if
   ///
   ///   iSig = GetSymmetrySignature(l,c,0)
   ///
   ///   then:
   ///
   ///   - If bit 0 is set ((iSig & (1u<<0)) != 0), fx = -1, otherwise fx = +1
   ///   - If bit 1 is set ((iSig & (1u<<1)) != 0), fy = -1, otherwise fy = +1
   ///   - If bit 2 is set ((iSig & (1u<<2)) != 0), fz = -1, otherwise fz = +1
   ///
   /// - 'Flags' is currently unused, but may later be used to distinguish
   ///   different basis function types. For Flags == 0, the returned signatures
   ///   are for the standard solid harmonics used in IR (in the order defined
   ///   at compile-time).
   unsigned GetSymmetrySignature(int l, size_t c, unsigned Flags);

   /// Compute modified spherical Bessel functions of the first kind, i_m(z),
   /// with envelopes. Concretely, compute for m in [0,1,...,MaxM] (MaxM
   /// inclusive!) the functions
   ///
   ///       pOut[m] = i_{m}(z) * Exp[-z] * z^{-m} * s^m * Factor
   ///
   /// where z := s * a (if you do not want a separate 's', pass a=z and s=1).
   /// 'Factor' is an optional common factor multiplied into all output values.
   void EvalModifiedSphericalBesselEnveloped(double *pOut, double s, double a, size_t MaxM, double Factor);

   void ZeroBlock2(double *pOut, size_t StrideA, size_t StrideB, unsigned nFnA, unsigned nFnB);
   void ZeroBlock3(double *pOut, size_t StrideA, size_t StrideB, size_t StrideC, unsigned nFnA, unsigned nFnB, unsigned nFnC);
   void ZeroBlock4(double *pOut, size_t StrideA, size_t StrideB, size_t StrideC, size_t StrideD, unsigned nFnA, unsigned nFnB, unsigned nFnC, unsigned nFnD);


#ifdef IR_INCLUDE_PRINT_FUNCTIONS
   void IrPrintMatrixGen(std::ostream &xout, double *pData, size_t nRows, size_t iRowSt, size_t nCols, size_t iColSt, std::string const &Caption, char const *pNumFmt="{:14.6f}");
   void operator << (std::ostream &xout, FRawShell const &rs);
#endif // IR_INCLUDE_PRINT_FUNCTIONS

} // namespace ir


namespace ir {
   // template<class T>
   // inline T sqr(T x) {
   //    return x*x;
   // }
   // ^- this causes gcc 4.9.2 to emit a call to __ieee754_pow_sse2 (libm-2.20.so), which takes up >10% of the total integration time. wtf?
   inline double sqr(double x) { return x*x; }
   inline unsigned sqr(unsigned x) { return x*x; }
   inline int sqr(int x) { return x*x; }
   inline double DistSq3(double const *pA, double const *pB) {
      return sqr(pA[0] - pB[0]) + sqr(pA[1] - pB[1]) + sqr(pA[2] - pB[2]);
   }

   template<class T>
   inline T pow3(T const &x) {
      return x*x*x;
   }

   static bool IsWithinRange(FRawShell const *pA, FRawShell const *pB) {
      if (!pA->pRange || !pB->pRange)
         return true; // shells have no screening data.
      return ir::sqr(pA->MaxCoRange() + pB->MaxCoRange()) >= DistSq3(pA->vCen, pB->vCen);
   }

   inline bool IsPrimitiveWithinRange(FRawShell const *pA, uint iExpA, FRawShell const *pB, uint iExpB, double fDistSqAB)
   {
      if (!pA->pRange || !pB->pRange)
         return true; // shells have no screening data.
      return ir::sqr(pA->ExpRange(iExpA) + pB->ExpRange(iExpB)) >= fDistSqAB;
   }

   static bool IsContractionWithinRange(FRawShell const *pA, uint iCoA, FRawShell const *pB, uint iCoB, double fDistSqAB)
   {
      if (!pA->pRange || !pB->pRange)
         return true; // shells have no screening data.
      return ir::sqr(pA->CoRange(iCoA) + pB->CoRange(iCoB)) >= fDistSqAB;
   }

   inline bool IsPrimitiveWithinRange(FRawShell const *pA, uint iExpA, FRawShell const *pB, uint iExpB)
   {
      return IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, DistSq3(pA->vCen, pB->vCen));
   }

   inline void dummy_to_suppress_unused_warnings() {
      (void)IsWithinRange;
      (void)IsContractionWithinRange;
   }
}

#endif // IR_H
