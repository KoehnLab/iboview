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

#include <sstream> // for PrintAligned
#include <ostream> // for PrintAligned
#include <vector> // for libmol export
#include "format.h" // for PrintAligned

#include <cmath>
#include "Ir.h"
#include "CtBasisShell.h"


namespace ct {

// Note: moved to IrCore.cpp
//
// static unsigned DoubleFactR(int l) {
//    unsigned r = 1;
//    while (l > 1) {
//       r *= l;
//       l -= 2;
//    }
//    return r;
// }
//
// double InvRawRawGaussNorm(double fExp, unsigned l)
// {
//    // purple book (6.6.14). (factor 4*pi/(2*l+1) converting between Ylm and Slm omitted)
//    return pow(2*fExp/M_PI,.75) * sqrt(pow(4.*fExp,l) / DoubleFactR(2*signed(l)-1));
// }
//
// double RawGaussNorm(double fExp, unsigned l)
// {
// //    return 1./InvRawGaussNorm(fExp, l);
//    return pow(M_PI/(2*fExp),.75) * sqrt(DoubleFactR(2*signed(l)-1)/pow(4.*fExp,l));
// }


FElementShellInfo::FElementShellInfo(std::string const &Definition_, std::string const &/*LibraryName_*/, std::string const &LibraryComment_)
   : Definition(Definition_), /*LibraryName(LibraryName_),*/ LibraryComment(LibraryComment_)
{
}


bool FElementShellInfo::operator < (FElementShellInfo const &other) const
{
   FElementShellInfo const
      &A = *this,
      &B = other;
#define COMPARE_STRING_FIELD(FieldName) \
   { \
      int icmp = A.FieldName.compare(B.FieldName); \
      if (icmp < 0) return true; \
      if (icmp > 0) return false; \
   }
   COMPARE_STRING_FIELD(Definition)
//    COMPARE_STRING_FIELD(LibraryName)
   COMPARE_STRING_FIELD(LibraryComment)
#undef COMPARE_STRING_FIELD
   return false;
}

bool FElementShellInfo::operator == (FElementShellInfo const &other) const
{
   FElementShellInfo const
      &A = *this,
      &B = other;
   return A.Definition == B.Definition /*&& A.LibraryName == B.LibraryName*/ && A.LibraryComment == B.LibraryComment;
}

bool FElementShellInfo::operator != (FElementShellInfo const &other) const
{
   return !this->operator == (other);
}



FAtomShell::FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_)
{
   Init(AngMom_, pExp_, nExp_, pCo_, nCo_, InitFlags_);
}

FAtomShell::FAtomShell(unsigned AngMom_, double const fExp_,  unsigned InitFlags_)
{
   double fCo = 1.;
   Init(AngMom_, &fExp_, 1, &fCo, 1, InitFlags_);
}

FAtomShell::FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, unsigned InitFlags_)
{
   TArray<double>
      CoMatrix_(nExp_ * nExp_, 0.);
   for (unsigned i = 0; i != nExp_; ++ i)
      CoMatrix_[(nExp_+1)*i] = 1.;
   Init(AngMom_, pExp_, nExp_, &CoMatrix[0], nExp_, InitFlags_);
}

void FAtomShell::Init(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_=0)
{
   AngMom = AngMom_;
   Exponents.resize(nExp_);
   CoMatrix.resize(nExp_ * nCo_);
   for (unsigned iExp = 0; iExp < nExp_; ++ iExp) {
      double
         fNorm = 1.;
      if (0 == (InitFlags_ & TYPE_Unnormalized))
         fNorm = 1./ir::RawGaussNorm(pExp_[iExp], AngMom);
      Exponents[iExp] = pExp_[iExp];
      for (unsigned iCo = 0; iCo < nCo_; ++ iCo)
         CoMatrix[iExp + iCo*nExp_] = fNorm * pCo_[iExp + iCo*nExp_];
   }
   iLibraryElement = 0; // not set.
}


// format a float in decimal exponential notation, possibly cutting out
// unnecessary digits for simple numbers
std::string FmtFloatForExport(double f)
{
//    if (f <= 1000 && f >= 0.0001) {
//       std::string
//          r = fmt::format("{:.14g}", f);
//       while (r.back() == '0')  // <- note: that can't be done in exponential notation!
//          r.pop_back();
//       return r;
//    }
   return fmt::format("{:.12e}", f);
}


void FAtomShell::ExportToLibmol(std::ostream &xout, std::string const &ElementName, std::string const &BasisName, std::string const &Comment) const
{
   // An example of the target format:
   // --
   // LI  s cc-pVDZ VDZ         :    9   3  1.09  1.09  9.09
   // basis set from gbasis database
   //   0.14690000E+04  0.22050000E+03  0.50260000E+02  0.14240000E+02  0.45810000E+01
   //   0.15800000E+01  0.56400000E+00  0.73450000E-01  0.28050000E-01  0.76600000E-03
   //   0.58920000E-02  0.29671000E-01  0.10918000E+00  0.28278900E+00  0.45312300E+00
   //   0.27477400E+00  0.97510000E-02 -0.31800000E-02 -0.12000000E-03 -0.92300000E-03
   //  -0.46890000E-02 -0.17682000E-01 -0.48902000E-01 -0.96009000E-01 -0.13638000E+00
   //   0.57510200E+00  0.51766100E+00  0.10000000E+01
   // --

   // first find primitive ranges spanned by the coefficients of the individual contractions.
   typedef std::pair<size_t, size_t>
      FExpRange;
   typedef std::vector<FExpRange>
      FExpRangeList;
   FExpRangeList
      ExpRanges;
   ExpRanges.reserve(nCo());
   for (size_t iCo = 0; iCo < nCo(); ++ iCo) {
      FExpRange
         ExpRange(0, nExp());
      // skip zero coefficients for leading primitives
      while (ExpRange.first < nExp() && fCo(ExpRange.first, iCo) == 0.)
         ExpRange.first += 1;
      // skip zero coefficients for trailing primitives
      while (ExpRange.second > ExpRange.first && fCo(ExpRange.second-1, iCo) == 0.)
         ExpRange.second -= 1;
      ExpRanges.push_back(ExpRange);
   }

   // make constraction declarations for caption line
   fmt::MemoryWriter
      wCo;
   for (size_t iCo = 0; iCo < nCo(); ++ iCo) {
      FExpRange const
         &ExpRange = ExpRanges[iCo];
      if (iCo != 0)
         wCo << " ";
      else
         wCo.write("{}.{}", (1 + ExpRange.first), ExpRange.second);
         // ^- note at +1 on .first: libmol starts indexing primitive exponents
         //    at 1, not 0, and the last index is inclusive, not exclusive, so
         //    .second doesn't need the +1
   }

   fmt::MemoryWriter w;
   // make caption line
   w.write("{:<4} {} {} : {:4} {:4} {}\n{}\n", ElementName, "spdfghiklm"[l()], BasisName, nExp(), nCo(), wCo.str(), Comment.c_str());

   // make a vector of floating point data entries to collect exponents and
   // coefficients into
   std::vector<double>
      Data;
   Data.reserve(nExp() + nExp() * nCo());

   // collect primitive exponents
   for (size_t iExp = 0; iExp != nExp(); ++ iExp) {
      Data.push_back(fExp(iExp));
   }

   // collect contraction coefficients
   for (size_t iCo = 0; iCo < nCo(); ++ iCo) {
      FExpRange const
         &ExpRange = ExpRanges[iCo];
      for (size_t iExp = ExpRange.first; iExp != ExpRange.second; ++ iExp) {
         double
            fThisCo = fCo(iExp, iCo),
            fThisExp = fExp(iExp);
         // renormalize contraction coefficients from IR format (contractions
         // are given in terms of raw Gaussian primitives) into library format
         // (contractions are given in terms of normalized Gaussian primitives)
         fThisCo *= ir::RawGaussNorm(fThisExp, l());
         Data.push_back(fThisCo);
      }
   }

   // export exponent & contraction coefficient data.
   {
      size_t
         nCharsOnLine = 0;
      char const
         *pNextSep = "";
      for (size_t iDataEntry = 0; iDataEntry != Data.size(); ++ iDataEntry) {
         std::string
            s = FmtFloatForExport(Data[iDataEntry]);
         if (s.size() + 1 + nCharsOnLine > 80) {
            pNextSep = "\n";
            nCharsOnLine = 0;
         }
         w.write("{}{}", pNextSep, s);
         nCharsOnLine += 1 + s.size();
         pNextSep = " ";
      }
   }

   // write to output stream and add final newline
   xout << w.str() << std::endl;
}


bool FAtomShell::operator < (FAtomShell const &other) const
{
   FAtomShell const
      &A = *this,
      &B = other;
   if (A.AngMom < B.AngMom) return true;
   if (B.AngMom < A.AngMom) return false;

   if (A.Exponents < B.Exponents) return true;
   if (B.Exponents < A.Exponents) return false;

   if (A.CoMatrix < B.CoMatrix) return true;
   if (B.CoMatrix < A.CoMatrix) return false;

   return false;

   // note: not compared: RangeInfo. Technically can be generated
   // if needed, and does not actually reflect the "logical" meaning
   // of the basis function shell.
}

bool FAtomShell::operator != (FAtomShell const &other) const
{
   return (*this < other) || (other < *this);
}

bool FAtomShell::operator == (FAtomShell const &other) const
{
   return !(*this != other);
}







ir::FRawShell FBasisShell::MakeIrShell()
{
   return ir::FRawShell(pAs->AngMom, &pAs->Exponents[0], pAs->nExp(),
      &pAs->CoMatrix[0], pAs->nCo(), &vCenter[0], pAs->RangeInfo.empty()? 0 : &pAs->RangeInfo[0]);
}


void FBasisShell::PrintAligned(std::ostream &xout, std::string const &Ind, unsigned TotalIndent) const
{
   std::streampos
      p0 = xout.tellp(),
      p1;
   p0 = xout.tellp();
   xout << fmt::format("{:3}: {:c}   {:8.4f} {:8.4f} {:8.4f}    ",
            nCo(), ("spdfghiklm"[l()]), vCenter[0], vCenter[1], vCenter[2]);
   p1 = xout.tellp();

   // ^- hm... this doesn't work. tellp() I mean.
//    p0 = 0;
//    p1 = 47;

   for (size_t iExp = 0; iExp < pAs->Exponents.size(); ++iExp){
      if (iExp != 0) {
         xout << "\n" << Ind;
         for (ptrdiff_t i = 0; i < ptrdiff_t(TotalIndent) + ptrdiff_t(p1) - ptrdiff_t(p0) - ptrdiff_t(Ind.size()); ++ i)
            xout << " ";
      }
      xout << fmt::format("{:16.7f}  ", pAs->Exponents[iExp]);

      double
         fRenorm = ir::RawGaussNorm(pAs->Exponents[iExp], l());
      std::stringstream
         str;
      for (size_t iCo = 0; iCo < nCo(); ++ iCo){
         double
            fCo = pAs->CoMatrix[nExp() * iCo + iExp];
         if (fCo != 0.)
            str << fmt::format(" {:9.5f}", (fCo*fRenorm));
         else
            str << fmt::format(" {:>9}", "  - - - -");
      }
      std::string
         s = str.str();
      if (0) {
         while(!s.empty() && (s[s.size()-1] == ' ' || s[s.size()-1] == '-' ))
            s.resize(s.size() - 1);
      }
      xout << s;
   }
}


// #ifdef INCLUDE_OPTIONALS

// Evaluate the radial part of a contracted Gauss function:
//    Out = r^l * \sum_i c[i] Exp[-zeta_i r^2]
// for a number of (input) radii. Output is a nRadii x nCo matrix,
// pOut[iPt + nRadii*iCo] giving the value of the radial basis function at point
// iPt for function iCo.
void EvalBfnRadial(double *RESTRICT pOut,
   double const *RESTRICT pExp, uint nExp, double const *RESTRICT pCo, uint nCo, uint l,
   double const *RESTRICT pRadii, uint nRadii)
{
   // exp(-50) ~= 2e-22.  36 would be good for 2e-16, but let's put some leeway into it.
   double const
      LogThrOrb = 50;
   uint const
      nMaxExp = 64; // not nice, but this way we need to deal with memory from the outside.
   if ( nExp > nMaxExp )
      throw std::runtime_error("ir_eval_bfn_radial ran out of memory. Too many exponents in contracted basis function!");

   double
//       fNrm = 1.;
      fNrm = std::sqrt(1./(2*l+1));
   for ( uint iPt = 0; iPt < nRadii; ++ iPt ) {
      double
         pExpV[nMaxExp],
         fDist = pRadii[iPt], // <- need that unsquared for r^l with odd l.
         fDistSq = fDist*fDist;
      for ( uint iExp = 0; iExp < nExp; ++ iExp ){
         double fExp = pExp[iExp];
         if ( fDistSq * fExp > LogThrOrb )
            pExpV[iExp] = 0;
         else
            pExpV[iExp] = std::exp(-fDistSq * fExp);
      };
      double const
         fPowR = std::pow(fDist, (int)l);
      for ( uint iCo = 0; iCo < nCo; ++ iCo) {
         // now the total contracted function.
         double
            v0 = 0;
         for ( uint iExp = 0; iExp < nExp; ++ iExp )
            v0 += pCo[iExp + nExp*iCo] * pExpV[iExp];
         pOut[iPt + nRadii*iCo] = fNrm * fPowR * v0;
      }
   }
};


// c/p'd from CtDftGrid.cpp
void MakeAtomRadialGrid(double *r, double *w, uint n, double AtomicScale)
{
   // main references:
   //   [1] JCP 102 346 (1995)   (Treutler & Ahlrichs)
   //   [2] JCP 108 3226 (1998)  (Krack & Koester)
   //   [3] JCP 88 2547 (1988)   (Becke)
   //   [4] JCP 104 9848 (1996)  (Mura & Knowles)
   double
      den = 1./(n+1.),
      ln05 = 1./std::log(2.),
      R = AtomicScale,
      Alpha = 0.6;
   // [1]; (T2) and (M4).
   // formula for T2: [2], eq.(9) and (10)
   for ( uint i = 1; i <= n; ++ i ){
      double
         // xi: -1 .. +1, wi: weights in that range.
         // will later be mapped to actual range.
         xi,wi,
         // derivative d[ri]/d[xi] (for weight)
         dri;
      if ( 0 ) {
         // chebychev pos&weights: beware of i starting at 1!
         double
            SinPhase = std::sin(i*M_PI*den),
            CosPhase = std::cos(i*M_PI*den),
            SinPhase2 = SinPhase*SinPhase;
            // T2 Chebychev abscissa x[i] (goes from (-1..+1), exclusive)
         xi = ((int)n + 1 - 2*(int)i)*den + (2./M_PI)*(1. + (2./3.)*SinPhase2)*
            CosPhase * SinPhase,
         // weight w[i].
         wi = (16./3.) * den * (SinPhase2 * SinPhase2);
      } else {
         // uniform weights. Should work fine for logarithmic grids.
         xi = ((int)n + 1 - 2*(int)i)*den;
         wi = 2.*den;
      }

      // note on a: the strange comments in [1] (eq. 20) refer to the
      // integration scheme: if integrating with a x = [0..1] quadrature scheme,
      // we should use a = 0, if using a x = [-1..1] quadrature scheme (like the T2
      // we employ here), then a = 1. A does thus not occur in the formulas

      // make mapping of x to actual radius, and init the grid weight to actual radii.
      if ( 1 ) {
         double r2 = 5.;
         double x = .5*xi + 0.5;
         r[i-1] = -r2 * std::log(1.0-x*x*x); // Mura & Knowles Log3 grid.
         dri = .5 * r2*3*x*x/(1.-x*x*x);
      } else if ( 0 ) {
         // simple logarithmic grid. works well with T2.
         r[i-1] = R*ln05*std::log(2./(1.-xi));
         dri = R*ln05/(1. - xi); // that is : dri/dxi = R/(1. - xi)
      } else {
         // Ahlrichs M4... works with simple linear weighting.
         double
            PowAlpha = std::pow(1 + xi, Alpha),
            Log2x = std::log(2/(1.-xi));
         r[i-1] = R*ln05* PowAlpha * Log2x;
         dri = R*ln05*PowAlpha*(1/(1.-xi) + Alpha * Log2x / (1 + xi));
      }
      w[i-1] = dri * wi * 4.*M_PI*r[i-1]*r[i-1]; // and that is the radial volume element.
   }
}

static double sqr(double x) { return x*x; }

// Find the effective range of a contracted basis function shell by explicitly
// evaluating the basis functions on a radial grid.
//
//  pOutMax[iCo] = (range r for which the integral Int[r,Infty] mu^2(r) d^3r <= ThrDen).
//
// i.e., ThrDen is a threshold on the neglected *electron number*. If we put one
// electron onto the basis function, r is the radius at which we lose less than
// ThrDen electrons if we neglect the function beyond that r. For this to work
// the input functions should best be of split-valence type.
void FindBfnRange(double *RESTRICT pOutMax, double *RESTRICT pOutMin,
   double const *RESTRICT pExp, uint nExp, double const *RESTRICT pCo, uint nCo, uint l,
   double const &ThrDen, uint nResolution, FMemoryStack &Mem)
{
   uint
      nPt = nResolution;
   double
      *pWork = Mem.AllocN(nPt*(2+nCo), 0.);
   double
      *pRadii = pWork,
      *pWeights = pWork + 1 * nPt,
      *pValues = pWork + 2 * nPt;
   // set up a dft-style logarithmic grid for the radial points
   MakeAtomRadialGrid(pRadii, pWeights, nPt, 1.0);
   // evaluate the basis functions on the grid...
   EvalBfnRadial(pValues, pExp, nExp, pCo, nCo, l, pRadii, nPt);

#if 0
   xout << fmt::format("Radial grids, weights, and values:\n");
   xout << fmt::format("rad[0] = {:.2e}  rad[nPt=1] = {:.2e}:\n", pRadii[0], pRadii[nPt-1]);
   PrintArray("Radii", pRadii, nPt);
   PrintArray("Weights", pWeights, nPt);
   for ( uint iCo = 0; iCo < nCo; ++ iCo )
      PrintArray("Values", &pValues[nPt*iCo], nPt);
#endif

   // ...and now go through the grid to find the effective ranges.
   double
      *fDen = Mem.AllocN(nCo, 0.);
   if (1) {
      if (pOutMax) {
         for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
            pOutMax[iCo] = 1e30;
            fDen[iCo] = 0.;
         }
         for ( uint iPt = 0; iPt < nPt; ++ iPt ) {
            for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
               fDen[iCo] += pWeights[iPt] * sqr(pValues[iPt + iCo*nPt]);
               if ( fDen[iCo] < ThrDen )
                  pOutMax[iCo] = pRadii[iPt];
            }
         }
      }
      if (pOutMin) {
         // now do the same for the inner radius.
         for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
            pOutMin[iCo] = 0.;
            fDen[iCo] = 0.;
         }
         for ( uint iPt_ = 0; iPt_ < nPt; ++ iPt_ ) {
            uint iPt = nPt - iPt_ - 1; // <- iterate in reverse order
            for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
               fDen[iCo] += pWeights[iPt] * sqr(pValues[iPt + iCo*nPt]);
               if ( fDen[iCo] < ThrDen )
                  pOutMin[iCo] = pRadii[iPt];
            }
         }
      }
   } else {
      // calculate sqrt(<r^2>)
      for ( uint iCo = 0; iCo < nCo; ++ iCo )
         pOutMax[iCo] = 0.;
      for ( uint iPt = 0; iPt < nPt; ++ iPt )
         for ( uint iCo = 0; iCo < nCo; ++ iCo )
            pOutMax[iCo] += pWeights[iPt] * sqr(pRadii[iPt]) * sqr(pValues[iPt + iCo*nPt]);
      for ( uint iCo = 0; iCo < nCo; ++ iCo )
         pOutMax[iCo] = -std::log(ThrDen) * std::sqrt(pOutMax[iCo]);
   }
#if 0
   // check if this works: if you put in normalized functions, fDen should now be about 1
   xout << fmt::format("FindBfnRange(nExp={}, nCo={}, l={}):\n", nExp, nCo, l);
   PrintArray("Integrated densities", &fDen[0], nCo);
   if (pOutMin)
      PrintArray("Effective ranges (Min)", pOutMin, nCo);
   if (pOutMax)
      PrintArray("Effective ranges (Max)", pOutMax, nCo);
#endif
   Mem.Free(pWork);
}

template<class T> T FindMax(T *pStart, size_t n, T StartValue = T(0)) {
   T fMax = StartValue;
   for (size_t i = 0; i < n; ++ i)
      fMax = std::max(pStart[i], fMax);
   return fMax;
}

template<class T> T FindMin(T *pStart, size_t n, T StartValue) {
   T fMin = StartValue;
   for (size_t i = 0; i < n; ++ i)
      fMin = std::min(pStart[i], fMin);
   return fMin;
}

// ThrEl: maximum fraction of an electron per basis function we are willing to lose.
// nRes: resolution for radial integration.
void FAtomShell::FindMinMaxRange(double *pMaxRange, double *pMinRange, uint nRes, double ThrEl) const
{
   // calculate screening ranges.
   FMemoryStack2
      Mem(1000 + ((nExp() + nRes) * (nExp() + nCo() + 2) * sizeof(double)));
   // make a contaction matrix containing first the primitives themselves,
   // and then the actual contractions. This is for calculating screening
   // information for the primitives, too.
   double
      *pCoX;
   Mem.ClearAlloc(pCoX, nExp() * (nExp() + nCo()));
   for (uint iExp = 0; iExp < nExp(); ++ iExp)
      pCoX[(nExp()+1)*iExp] = 1./ir::RawGaussNorm(fExp(iExp), AngMom);
   for (uint iCo = 0; iCo < nCo(); ++ iCo)
      for (uint iExp = 0; iExp < nExp(); ++ iExp)
         pCoX[nExp() * (nExp() + iCo) + iExp] = fCo(iExp,iCo);
//    RangeInfo.resize(1 + nExp() + nCo());
   FindBfnRange(&pMaxRange[1], pMinRange? &pMinRange[1] : 0, &Exponents[0], nExp(),
      pCoX, nExp() + nCo(), AngMom, ThrEl, nRes, Mem);

   pMaxRange[0] = FindMax(&pMaxRange[1+nExp()], nCo(), 0.);
   if (pMinRange)
      pMinRange[0] = FindMin(&pMinRange[1+nExp()], nCo(), 1e30);

      // hm.. something is wrong.
//       for (uint iExp = 0; iExp < nExp(); ++ iExp)
//          pMaxRange[1+iExp] = std::sqrt(-std::log(ThrEl)/Exponents[iExp]);
      // ^-   exp(-z*r^2) = thr
      //      -z*r^2 = log(thr)
      //      r^2 = -log(thr)/z
      //      r = sqrt(-log(thr)/z)
//       PrintArray("Output ranges", &pMaxRange[0], nCo()+nExp()+1);
#if 0
      FBasisShell dummy(FVector3(0.,0.,0.),0, this);
      xout << "\n";
      xout << "\n";
      xout << "\n";
      dummy.PrintAligned(xout,12);
      xout << "\n";
//       PrintMatrixGen(xout, &pMaxRange[0], 1, 1, nCo()+nExp()+1, 1, "SCREEN INFO" );
      PrintArray("Output ranges", &pMaxRange[0], nCo()+nExp()+1);
      PrintMatrixGen(xout, pCoX, nExp(), 1, nCo()+nExp(), nExp(), "COX" );
#endif
}


// #endif // INCLUDE_OPTIONALS


void FAtomShell::Finalize()
{
// #ifdef INCLUDE_OPTIONALS
   if (RangeInfo.empty()) {
      // calculate screening ranges.
      RangeInfo.resize(1 + nExp() + nCo());
      FindMinMaxRange(&RangeInfo[0], 0, 2000, 1e-15);
//       FindMinMaxRange(&RangeInfo[0], 0, 2000, 1e-10);
   }
// #endif // INCLUDE_OPTIONALS
}


static std::string const
   s_EmptyString;

std::string const &FAtomShell::FAtomShell::Definition() const
{
   if (pInfo)
      return pInfo->Definition;
   else
      return s_EmptyString;
}


// std::string const &FAtomShell::LibraryName() const
// {
//    if (pInfo)
//       return pInfo->Definition;
//    else
//       return s_EmptyString;
// }


std::string const &FAtomShell::LibraryComment() const
{
   if (pInfo)
      return pInfo->Definition;
   else
      return s_EmptyString;
}


void FAtomShell::AttachInfo(int iElement_, FElementShellInfoCptr pInfo_)
{
   assert(iElement_ >= 0);
   this->iLibraryElement = iElement_;
   this->pInfo = pInfo_;
}


std::string const &FBasisShell::LibraryName() const
{
   if (this->pLibraryName)
      return *this->pLibraryName;
   else
      return s_EmptyString;
}


std::string const &FBasisShell::LibraryComment() const {
   return pAs->LibraryComment();
}


int FBasisShell::iLibraryElement() const {
   return pAs->iLibraryElement;
}


} // namespace ct
