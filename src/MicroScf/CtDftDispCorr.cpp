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

// c++ -c -O3 -DNDEBUG -I../wmme CtDftDispCorr.cpp && nm --demangle --size-sort --print-size --format=b --radix=d d3bj.o
#include <string>
#include <algorithm>
#include "format.h"

#include "CxTypes.h"
#include "CtAtomSet.h"
#include "CtMatrix.h" // only for Gradient.
#include "CxIo.h"


// stores lots of parameters as big static tables...
#include "CtDftDispCorr.h"
#include "CtDftDispCorr.inl"
#include "coords_meta.h"

namespace ct {

// note: this is an implementation of Grimme's D3(BJ) correction for MicroScf, based on parameter files from dftd3.f:
//    [(V3.2 Rev 0), obtained from http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=getd3&lang=english, 2016-06-28]

static double const
    // global ad hoc parameters (from common/param)
    d3bj_k1 = 16.0,
    d3bj_k2 = 4./3., // shouldn't this go into the D3 coordination numbers? (d3cn). It's there in D3 itself. Or was it removed in D3(BJ)?
    d3bj_k3 = -4.;

static double const
    // threshold for neglect of contributions to coordination numbers (in rij^2 form).
    d3bj_cn_thr_rjisq = 1600.0, // <- called 'cn_thr' in dftd3.f
    d3bj_rthr2 = 1600.0,
    d3bj_rthr = 9000.0;


struct FD3DispersionContextImpl
{
   FAtomSet const
      &Atoms;

   explicit FD3DispersionContextImpl(FLog *pLog_, FAtomSet const *pAtoms, std::string const &DispCorrType, std::string const &XcName, FMemoryStack &Mem);
   ~FD3DispersionContextImpl();

   // compute and return D3(BJ) dispersion energy for initialized atom set.
   // if Gradient.pData != 0, also compute analytic gradient of D3(BJ)
   // dispersion energy and add it to 'Gradient'. Gradient should be a 3 x
   // nAtoms matrix.
   double CalcEnergyAndGradient(FMatrixView Gradient, FMemoryStack &Mem) const;
   void EstimateHessian(FMatrixView Hessian, FMemoryStack &Mem, bool Add) const;
protected:
   FLog
      *m_pLog;
   FMemoryStack
      *m_pOrigMem;
   double
      // allocated on input memory stack.
      *m_pEffectiveCn,
      *m_pC6Params, // big array of reference data... (as in copyc6 of dftd3.f)
      // functional-specific constants.
      m_s6, m_s18, // prefactors for C6/C8 energy terms
      m_rs6, m_rs18; // parameters for damping (here: a1/a2 parameters for D3(BJ))

   static size_t const
      MaxElem = 94, // highest tabulated element number
      MaxRefCn = 5; // highest tabulated number of reference variants of coordination numbers/element (mxc/maxc in dftd3.f).
   std::string
      m_XcName;
   int
      m_Verbosity;

   void MakeCoordinationNumberIncrement(double &Cn_ij, double &dCndRijsq, int iAt, int jAt) const;
   void MakeCoordinationNumbers(double *&pEffectiveCn, FMemoryStack &Mem);
   double GetSqrtR2R4(int iElemA, int iElemB) const;
   void GetFunctionalParameters(std::string const &XcName);
   void GetC6(double &C6, double &dC6dCn_i, double &dC6dCn_j, bool MakeGradient,
      int iZ, int jZ, double cni, double cnj, FMemoryStack &Mem) const;
};


FD3DispersionContextImpl::FD3DispersionContextImpl(FLog *pLog_, FAtomSet const *pAtoms_, std::string const &DispCorrType, std::string const &XcName, FMemoryStack &Mem)
   : Atoms(*pAtoms_), m_pLog(pLog_), m_pOrigMem(&Mem), m_pEffectiveCn(0), m_rs6(-1.), m_rs18(-1.)
{
   m_Verbosity = 2;
   if (pLog_ == 0)
      m_Verbosity = -1;
   if (!(DispCorrType == "D3(BJ)" || DispCorrType == "d3(bj)" || DispCorrType == "on"))
      throw std::runtime_error(fmt::format("DFT-D3 correction: Unrecognized dispersion correction type '{}'. Know 'on' and 'D3(BJ)' (use 'off' or 'none' to disable).", DispCorrType));
   GetFunctionalParameters(XcName);
   // this allocates m_pEffectiveCn.
   MakeCoordinationNumbers(m_pEffectiveCn, Mem);
}


FD3DispersionContextImpl::~FD3DispersionContextImpl()
{
//    m_pOrigMem->Free(m_pEffectiveCn);
}


void FD3DispersionContextImpl::GetFunctionalParameters(std::string const &XcName)
{
   m_XcName = XcName;
   m_s6 = 1.;
   if (XcName == "PBE") {
      m_rs6 = 0.4289;
      m_s18 = 0.7875;
      m_rs18 = 4.4407;
   } else if (XcName == "TPSS" || XcName == "LL-TPSS") {
      m_rs6 = 0.4535;
      m_s18 = 1.9435;
      m_rs18 = 4.4752;
      m_XcName = "TPSS";
   } else if (XcName == "PBEh-3c") {
      m_rs6 = 0.4860;
      m_s18 = 0.0000;
      m_rs18 = 4.5000;
   } else if (XcName == "B-LYP" || XcName == "BLYP") {
      m_rs6 = 0.4298;
      m_s18 = 2.6996;
      m_rs18 = 4.2359;
      m_XcName = "BLYP";
   } else if (XcName == "B-P" || XcName == "B-P86" || XcName == "BP86") {
      m_rs6 =0.3946;
      m_s18 =3.2822;
      m_rs18 = 4.8516;
      m_XcName = "B-P86";
   } else if (XcName == "HF" || XcName == "RHF" || XcName == "") {
      m_rs6 = 0.3385;
      m_s18 = 0.9171;
      m_rs18 = 2.8830;
      m_XcName = "HF";
   } else if (XcName == "TPSS0") {
      m_rs6 = 0.3768;
      m_s18 = 1.2576;
      m_rs18 = 4.5865;
   } else if (XcName == "TPSSh") {
      m_rs6 = 0.4529;
      m_s18 = 2.2382;
      m_rs18= 4.6550;
   } else if (XcName == "PBE0") {
      m_rs6 = 0.4145;
      m_s18 = 1.2177;
      m_rs18 = 4.8593;
   } else {
      m_rs6 = 0.;
      m_s18 = 0.;
      m_rs18 = 0.;
      throw std::runtime_error(fmt::format("DFT-D3 correction: no parameters for functional '{}'.", XcName));
   }

   // absorb -1 factor of total dispersion energy (it is supposed to be negative!) into s6/s18 terms.
   m_s6 *= -1.;
   m_s18 *= -1.;
}



// D3(BJ) two-body dispersion energy contribution (coordinates2/edisp_d3bj_e.max and lhf/d3-bj/edisp_d3bj1.max)
static void edisp_d3bj(double &E, double &dEdRijsq, double &dEdC6, double rijsq, double c6, double r2r4_term, double a1, double a2, double s6, double s18)
{
   double r1 = 1.732050807568877*a1*sqrt(r2r4_term) + a2;
   double r1_pow2 = r1*r1;
   double r1_pow6 = r1_pow2*r1_pow2*r1_pow2;
   double rijsq_pow3 = rijsq*rijsq*rijsq;
   double r2 = rijsq_pow3*rijsq + r1_pow6*r1_pow2;
   double r3 = 1.0/r2;
   double r4 = rijsq_pow3;
   double r5 = r4 + r1_pow6;
   double r6 = 1.0/r5;
   E = c6*r6*s6 + 3.0*r3*c6*r2r4_term*s18;
   dEdRijsq = -3.0*c6*sqr(rijsq)*s6/sqr(r5) - 12.0*r4*c6*r2r4_term*s18/sqr(r2);
   dEdC6 = r6*s6 + 3.0*r3*r2r4_term*s18;
}


// D3(BJ) effective coordination number increment (lhf/d3-bj/edisp_d3cn.max)
static void edisp_d3cn(double &Cn, double &dCndRijsq, double rijsq, double rcovij, double d3bj_k1)
{
   double r1 = sqrt(rijsq);
   double r2 = 1.0/std::exp((rcovij/r1 - 1.0)*d3bj_k1);
   double r3 = r2 + 1;
   Cn = 1.0/r3;
   dCndRijsq = -0.5*r2*d3bj_k1*rcovij/(r1*rijsq*sqr(r3));
}


void FD3DispersionContextImpl::MakeCoordinationNumberIncrement(double &Cn_ij, double &dCndRijsq, int iAt, int jAt) const
{
   assert(iAt != jAt);
   // get link to tabulated covalent radii.
   double const
      *pCovR = &dftd3::g_C6RcovTable[0];
   size_t
      nCovR = sizeof(dftd3::g_C6RcovTable)/sizeof(dftd3::g_C6RcovTable[0]);

   // compute actual increment in dependence of atoms I and J.
   FVector3
      dXyz = Atoms[iAt].vPos - Atoms[jAt].vPos;
   double
      rijsq = dXyz.LengthSq();
   if (rijsq > d3bj_cn_thr_rjisq) {
      Cn_ij = 0.;
      dCndRijsq = 0.;
   } else {
      int
         iElemI = Atoms[iAt].iElement,
         iElemJ = Atoms[jAt].iElement;
      if (!(size_t(iElemI) < nCovR && size_t(iElemJ) < nCovR))
         throw std::runtime_error(fmt::format("DFT-D3 correction: ran out of covalent radii. Zi = {}, Zj = {}. Have: Z = 0...{}.", iElemI, iElemJ, nCovR-1));
      double
         // covalent distance in Bohr radii
         rcovij = pCovR[iElemI] + pCovR[iElemJ];
      edisp_d3cn(Cn_ij, dCndRijsq, rijsq, rcovij, d3bj_k1);
      // ^- note: formula in D3 paper has also the d3bj_k2 constant in there...
      //    but my program doesn't, and neither does Grimme's own dftd3.F. And both
      //    mine and theirs generated identical results...
//       m_pLog->Write("   {:<} -- {:<}   r = {:16.8f}  rcov = {:16.8f}  -> CN_ij = {:16.8f}", Atoms[iAt].GetAtomLabel(iAt), Atoms[jAt].GetAtomLabel(jAt), std::sqrt(rijsq), rcovij, Cn_ij);
   }
}


void FD3DispersionContextImpl::MakeCoordinationNumbers(double *&pEffectiveCn, FMemoryStack &Mem)
{
   // port of dftd3.f:ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
   // note on variable renames:
   //    cn(i) -> pEffectiveCn(i)
   //    i -> iAt
   //    iAt -> jAt
   //    cn_thr -> cn_thr_rjisq
   size_t
      nAt = Atoms.size();
   Mem.Alloc(pEffectiveCn, nAt);
   for (size_t iAt = 0; iAt != nAt; ++ iAt) {
      double cn_i = 0.; // called 'xn' in dftd3.f
      for (size_t jAt = 0; jAt != nAt; ++ jAt) {
         if (iAt == jAt)
            continue;
         double
            Cn_ij, dCndRijsq;
         MakeCoordinationNumberIncrement(Cn_ij, dCndRijsq, iAt, jAt);
         cn_i += Cn_ij;
      }
      pEffectiveCn[iAt] = cn_i;
//       if (m_Verbosity >= 2) {
//          m_pLog->Write("   {:<}   CN: {:16.8f}", Atoms[iAt].GetAtomLabel(iAt), cn_i);
//          m_pLog->WriteLine();
//       }
   }
}


static void GetTabulatedC6(FLog *m_pLog, double &c6, double &cni, double &cnj, int iZ_, int jZ_, int icn_, int jcn_)
{
   bool swapped;
   int iZ, jZ, icn, jcn;
   if ((iZ_ > jZ_) || (iZ_ == jZ_ && (icn_ >= jcn_))) {
      iZ = iZ_; jZ = jZ_; icn = icn_; jcn = jcn_;
      swapped = false;
   } else {
      iZ = jZ_; jZ = iZ_; icn = jcn_; jcn = icn_;
      swapped = true;
   }
   size_t
      iPair = size_t((iZ*(iZ-1))/2 + (jZ-1));
   if (iPair >= sizeof(dftd3::g_C6ParameterIndexTable)/sizeof(dftd3::g_C6ParameterIndexTable[0]))
      throw std::runtime_error(fmt::format("DFT-D3 correction: no data entry for element pair Za = {} // Zb = {}", iZ, jZ));
   dftd3::FC6EntryIndex const
      &PairIndex = dftd3::g_C6ParameterIndexTable[iPair];
   assert(PairIndex.iz == iZ && PairIndex.jz == jZ);

   for (size_t iDataLine = PairIndex.ibeg; iDataLine != PairIndex.iend; ++ iDataLine) {
      dftd3::FC6Entry const &DataEntry = dftd3::g_C6ParameterTable[iDataLine];
      assert(DataEntry.iz == iZ && DataEntry.jz == jZ);
      if (DataEntry.icn == icn && DataEntry.jcn == jcn) {
         c6 = DataEntry.c6;
         if (swapped) {
            cni = DataEntry.cnj;
            cnj = DataEntry.cni;
         } else {
            cni = DataEntry.cni;
            cnj = DataEntry.cnj;
         }
//          m_pLog->Write("              getc6(iZ={:2},jZ={:2}; {}/{})   ->  c6 ={:14.6}  cni={:14.6}  cnj={:14.6}",iZ_, jZ_, icn_, jcn_, c6, cni, cnj);
         return;
      }
   }
   // indicate non-existence of data (not all data pairs asked for are there).
   c6 = 0.;
   cni = 0.;
   cnj = 0.;
//    m_pLog->Write("              getc6(iZ={:2},jZ={:2}; {}/{})   ->  no entry.", iZ_, jZ_, icn_, jcn_);
   IR_SUPPRESS_UNUSED_WARNING(m_pLog);
}


void FD3DispersionContextImpl::GetC6(double &C6, double &dC6dCn_i, double &dC6dCn_j, bool MakeGradient,
      int iZ, int jZ, double cni, double cnj, FMemoryStack &Mem) const
{
   C6 = 0.;
   dC6dCn_i = 0.;
   dC6dCn_j = 0.;

   if (iZ < 1 || jZ < 1)
      return; // no dispersion interaction with dummy atoms.
   size_t
      nCiTabEntries = sizeof(dftd3::g_C6MaxCiTable)/sizeof(dftd3::g_C6MaxCiTable[0]);
   if (size_t(iZ) >= nCiTabEntries || size_t(jZ) >= nCiTabEntries)
      throw std::runtime_error(fmt::format("DFT-D3 correction: ran out of c6 coefficient co-indices. Zi = {}, Zj = {}. Have: Z = 0...{}.", iZ, jZ, nCiTabEntries-1));

   // compute c6 from interpolation over reference c6 values for various coordination numbers.
   // D3 paper, eq. 16
   double
      fWeightSum = 0.,
      dWeightSumdCn_i = 0.,
      dWeightSumdCn_j = 0.;
   for (size_t iCnRef = 0; iCnRef != dftd3::g_C6MaxCiTable[size_t(iZ)]; ++ iCnRef) {
      for (size_t jCnRef = 0; jCnRef != dftd3::g_C6MaxCiTable[size_t(jZ)]; ++ jCnRef) {
         double
            fRefC6,
            fRefCnI,
            fRefCnJ;
         GetTabulatedC6(m_pLog, fRefC6, fRefCnI, fRefCnJ, iZ, jZ, iCnRef, jCnRef);
         if (fRefC6 == 0.) // no data.
            continue;
         double
            fWeight1 = d3bj_k3*(sqr(cni - fRefCnI) + sqr(cnj - fRefCnJ)),
            fWeight = 0.;
         if (fWeight1 < -36.)
            continue;
         fWeight = std::exp(fWeight1);
         fWeightSum += fWeight;
         C6 += fWeight * fRefC6;

         if (MakeGradient) {
            double dWeightdCn_i = fWeight * (d3bj_k3*2*(cni - fRefCnI));
            double dWeightdCn_j = fWeight * (d3bj_k3*2*(cnj - fRefCnJ));
            dC6dCn_i += fRefC6 * dWeightdCn_i;
            dC6dCn_j += fRefC6 * dWeightdCn_j;
            dWeightSumdCn_i += dWeightdCn_i;
            dWeightSumdCn_j += dWeightdCn_j;
         }
      }
   }

   if (fWeightSum != 0.) {
      double fInvWeight = 1./fWeightSum;
      C6 *= fInvWeight;
      if (MakeGradient) {
         dC6dCn_i = fInvWeight * (dC6dCn_i - C6 * dWeightSumdCn_i);
         dC6dCn_j = fInvWeight * (dC6dCn_j - C6 * dWeightSumdCn_j);
      }
   } else {
      C6 = 0.;
      dC6dCn_i = 0.;
      dC6dCn_j = 0.;
   }
   IR_SUPPRESS_UNUSED_WARNING(Mem);
}



double FD3DispersionContextImpl::CalcEnergyAndGradient(FMatrixView Gradient, FMemoryStack &Mem) const
{
   TMemoryLock<double> pFreeMe(0, &Mem);
   assert(m_pEffectiveCn != 0);
   bool
      MakeGradient = (Gradient.pData != 0);
   assert(!MakeGradient || (Gradient.nRows == 3 && Gradient.nCols == Atoms.size()));
   size_t
      nAt = Atoms.size();

   // for the gradient we are going to accumulate:
   // - \frac{\partial E_\mathrm{disp}{\partial R_i} directly (direct rijsq term, with c6 fixed)
   // - and \frac{\partial E_\mathrm{disp}{\partial cn_i}, via dE/dc6 and separate accs for dc6/dcni // dc6/dcnj.
   // In a second step, we then compute the derivatives of the coordination numbers and mix it with the intermediate dEdcn derivatives.
   double
      *dEdCn = 0;  // d[total E_disp]/d[cn_i]
   if (MakeGradient)
      Mem.ClearAlloc(dEdCn, nAt);

   double const
      a1 = m_rs6, // functional specific constants (should have been loaded already)
      a2 = m_rs18;

   double
      fEnergy = 0.;
   for (size_t iAt = 0; iAt < nAt; ++ iAt) {
      for (size_t jAt = 0; jAt < iAt; ++ jAt) {
         FVector3
            dXyz = Atoms[iAt].vPos - Atoms[jAt].vPos;
         double
            rijsq = dXyz.LengthSq();
         int
            iElemI = Atoms[iAt].iElement,
            iElemJ = Atoms[jAt].iElement;

         if (rijsq > d3bj_cn_thr_rjisq)
            continue;
         double
            c6 = 0., dC6dCn_i = 0., dC6dCn_j = 0.;
         GetC6(c6, dC6dCn_i, dC6dCn_j, MakeGradient,
            iElemI, iElemJ, m_pEffectiveCn[iAt], m_pEffectiveCn[jAt], Mem);
         if (c6 == 0.)
            continue;
         double
            r2r4_term = GetSqrtR2R4(iElemI, iElemJ);
         {
            double Eij, dEijdRijsq, dEijdC6;
            edisp_d3bj(Eij, dEijdRijsq, dEijdC6, rijsq, c6, r2r4_term, a1, a2, m_s6, m_s18);
            fEnergy += Eij;
//             m_pLog->Write(" edisp(ij):   {:4}   {:4}   c6 ={:14.8f}   c8 ={:14.8f}   eij ={:14.8f}", iAt+1, jAt+1, c6, 3*c6*r2r4_term, Eij);
            if (MakeGradient) {
               // compute direct contribution to gradient (via rijsq).
               for (size_t ixyz = 0; ixyz != 3; ++ ixyz) {
                  Gradient(ixyz,iAt) += 2*dXyz[ixyz] * dEijdRijsq;
                  Gradient(ixyz,jAt) -= 2*dXyz[ixyz] * dEijdRijsq;
               }
               // remember contributions to dEtotal/dcn[i]
               dEdCn[iAt] += dEijdC6 * dC6dCn_i;
               dEdCn[jAt] += dEijdC6 * dC6dCn_j;
            }
         }
      }
   }

   if (MakeGradient) {
      // add contributions to energy via differentials of effective coordination numbers.
      for (size_t iAt = 0; iAt < nAt; ++ iAt) {
         for (size_t jAt = 0; jAt < iAt; ++ jAt) {
            FVector3
               dXyz = Atoms[iAt].vPos - Atoms[jAt].vPos;
            double
               Cn_ij, dCndRijsq;
            MakeCoordinationNumberIncrement(Cn_ij, dCndRijsq, iAt, jAt);
            double
               rg = (dEdCn[iAt] + dEdCn[jAt]) * 2 * dCndRijsq;
            for (size_t ixyz = 0; ixyz != 3; ++ ixyz) {
               Gradient(ixyz,iAt) += rg*dXyz[ixyz];
               Gradient(ixyz,jAt) -= rg*dXyz[ixyz];
            }
         }
      }
   }
   return fEnergy;
   // note: by default D3 and D3(BJ) have no three-body terms. One can turn them on,
   // but they are off by default. Analytic gradients are not implemented for those
   // in dftd3.f as of version 3.2 (June 2016). This suggests that they are not widely used.
}



void FD3DispersionContextImpl::EstimateHessian(FMatrixView Hessian, FMemoryStack &Mem, bool Add) const
{
   TMemoryLock<double> pFreeMe(0, &Mem);
   double const
      a1 = m_rs6, // functional specific constants (should have been loaded already)
      a2 = m_rs18;
   if (!Add)
      Hessian.Clear();
   size_t
      nAt = Atoms.size();
   assert(Hessian.nRows == 3*nAt && Hessian.nCols == 3*nAt);

   // re-pack atomic coordinates into format used in coords_meta.h
   // (as a single flat position array)
   TMemoryLock<double>
      pAtXyz(3*nAt, &Mem);
   for (size_t iAt = 0; iAt < nAt; ++ iAt)
      for (size_t ixyz = 0; ixyz < 3; ++ ixyz)
         pAtXyz[ixyz + 3*iAt] = Atoms[iAt].vPos[ixyz];

   for (size_t iAt = 0; iAt < nAt; ++ iAt) {
      for (size_t jAt = 0; jAt < iAt; ++ jAt) {
         FVector3
            dXyz = Atoms[iAt].vPos - Atoms[jAt].vPos;
         double
            rijsq = dXyz.LengthSq();
         int
            iElemI = Atoms[iAt].iElement,
            iElemJ = Atoms[jAt].iElement;
         if (rijsq > d3bj_cn_thr_rjisq)
            continue;
         double
            c6 = 0., dC6dCn_i = 0., dC6dCn_j = 0.;
         GetC6(c6, dC6dCn_i, dC6dCn_j, false,
            iElemI, iElemJ, m_pEffectiveCn[iAt], m_pEffectiveCn[jAt], Mem);
         if (c6 == 0.)
            continue;
         double
            r2r4_term = GetSqrtR2R4(iElemI, iElemJ);
         rx::index_t
            iAtoms[2] = {iAt, jAt};
         assert(Hessian.nRowSt == 1);
         // this evaluates the derivative of the energy under the side-condition that the c6 coefficient is NOT
         // geometry dependent. I.e., the use the c6 from the reference geometry for this.
         // Since this is only meant to improve geometry convergence a bit, I don't think this is a super-bad
         // approximation.
         rx::ef_d3bj_hess(0, Hessian.pData, Hessian.nColSt, 1., &pAtXyz[0], &iAtoms[0], c6, r2r4_term, a1, a2, m_s6, m_s18);
      }
   }
}




double FD3DispersionContextImpl::GetSqrtR2R4(int iElemA, int iElemB) const
{
   size_t nEntriesR2R4 = sizeof(dftd3::g_C6R2r4TermTable) / sizeof(dftd3::g_C6R2r4TermTable[0]);
   if (size_t(iElemA) > nEntriesR2R4 || size_t(iElemB) > nEntriesR2R4)
      throw std::runtime_error(fmt::format("DFT-D3 correction: ran out of r2r4 terms. Zi = {}, Zj = {}. Have: Z = 0...{}.", iElemA, iElemB, nEntriesR2R4-1));
   return dftd3::g_C6R2r4TermTable[size_t(iElemA)] * dftd3::g_C6R2r4TermTable[size_t(iElemB)];
}



// // c cut-off radii for all element pairs
// //       real*8 r0ab(max_elem,max_elem)
// // c C6 for all element pairs
// //       real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
//
//
//
// static double const GetR0ab(int iZ, int jZ) {
//    static const double autoang = 0.52917726;
//
//    if (iZ < 1 || iZ > d3_MaxElem)
//       return 0.;
//    if (jZ < 1 || jZ > d3_MaxElem)
//       return 0.;
//    size_t
//       i_ = std::min(size_t(iZ-1), size_t(jZ-1)),
//       j_ = std::max(size_t(iZ-1), size_t(jZ-1));
//    return d3_r0ab[((j_*(j_+1))/2 + i_)]/autoang;
// }

// ^- not needed for D3(BJ).




FD3DispersionContext::FD3DispersionContext(FLog *pLog_, FAtomSet const *pAtoms, std::string const &DispCorrType, std::string const &XcName, FMemoryStack &Mem)
   : p(new FD3DispersionContextImpl(pLog_, pAtoms, DispCorrType, XcName, Mem))
{
}


FD3DispersionContext::~FD3DispersionContext()
{
   delete p;
}


double FD3DispersionContext::CalcEnergyAndGradient(FMatrixView Gradient, FMemoryStack &Mem) const
{
   return p->CalcEnergyAndGradient(Gradient, Mem);
}

void FD3DispersionContext::EstimateHessian(FMatrixView Hessian, FMemoryStack &Mem, bool Add) const
{
   return p->EstimateHessian(Hessian, Mem, Add);
}


} // namespace ct
