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

#include "Iv.h"
#include <iostream>
#include <stdexcept>

// #include <QFileInfo>
#include <cmath>
#include <algorithm>
#include <set>
#include <fstream> // for orbital basis export

#include <QTextStream>

#include "format.h"
#include "CxIo.h"
#include "CxColor.h"
#include "CxPhysicalUnits.h"
#include "CxOpenMpProxy.h"
// #include "CxXyzFrameIo.h"

#include "IvDocument.h"
#include "IvSettings.h"
// #include "IvFileConvert.h"
#include "IvOrbital.h"

// #include "CtInt1e.h"
// #include "IrAmrr.h"
#include "CtMatrix.h"
#include "CtRhf.h"
#include "CtVolumeProperty.h"
// #include "CtBasisLibrary.h"


// TODO: more-refactoring. Related code in:
// - IvIao.h/.cpp
// - IvAnalysis.h/.cpp
// - IvDocument.h/.cpp


static const bool g_UseOccupationNumbersForIboChangeCurves = false;



FOrbitalInfo::FOrbitalInfo(double fEnergy_, double fOcc_, int iSym_, wfi::FOrbitalSpin OrbSpin_)
   : fEnergy(fEnergy_), fOcc(fOcc_), iSym(iSym_), iOriginalIndex(-1)
{
   using namespace wfi;
   // fix up occupation numbers close to 0.0, 1.0, and 2.0.
   // These might have been messed up due to input reading.
   double const fEps = 1e-6;
   if (std::abs(fOcc - 2.) < fEps)
      fOcc = 2.;
   if (std::abs(fOcc - 1.) < fEps)
      fOcc = 1.;
   if (std::abs(fOcc - 0.) < fEps)
      fOcc = 0.;
   // should we do something about negative occupations? We could see some in transition
   // densities, or perturbative RDMs.

   if (OrbSpin_ != wfi::ORBSPIN_Unknown)
      Spin = OrbSpin_;
   else {
      // Nothing explicit set, so we guess. Assume high-spin open-shell WF by
      // default (with closed/alpha/empty orbitals).
      if (std::abs(fOcc - 2.) < fEps) {
         Spin = wfi::ORBSPIN_SpinFree;
         fOcc = 2.;
      } else if (std::abs(fOcc - 0.) < fEps) {
         Spin = wfi::ORBSPIN_SpinFree;
         fOcc = 0.;
      } else if (std::abs(fOcc - 1.) < fEps) {
         Spin = wfi::ORBSPIN_Alpha;
         fOcc = 1.;
      } else
         Spin = wfi::ORBSPIN_SpinFree;
   }
}

FOrbital::FOrbital(QString const &Desc_, ct::FAtomSetPtr pAtoms_, FDocument *pDocument_, FBasisSetPtr pBasisSet_, double *pCoeffs_, FOrbitalInfo const &info_, FVolumeVisualConfigPtr pRefVisConfig, double *pDm, double *pQm)
   : FBase(Desc_, pAtoms_, pDocument_, pRefVisConfig), info(info_), pBasisSet(pBasisSet_)
{
   // copy over coefficients.
   pCoeffs.insert(pCoeffs.end(), pCoeffs_, pCoeffs_ + pBasisSet->nFn());

   // if dipole and quadrupole moments were supplied, store those too. These
   // are used for setting up the viewing box and orientation of the orbital.
   HaveMoments = false;
   assert((pDm != 0) == (pQm != 0));
   if (pDm) {
      for (uint i = 0; i < 3; ++ i)
         vDipMom[i] = pDm[i];
      for (uint i = 0; i < 3; ++ i)
         for (uint j = 0; j < 3; ++ j)
            mQuadMom[i][j] = pQm[3*i + j];
      HaveMoments = true;
   }
}

FOrbital::~FOrbital()
{
}


QString FOrbital::GetType() const {
   if (pGlMesh != 0)
      return "O*"; // to mark which ones have data?
   return "O"; // orbital
};


uint32_t FOrbital::GetBaseColor() const
{
   FOrbital const *pOrbital = this;
   uint32_t dwColor;

   if (pOrbital->pVisConfig->pDetails.get()) {
//    dwColor = (FColor(0.5f*(FColor(pOrbital->pVisConfig->cIsoPlus) + FColor(pOrbital->pVisConfig->cIsoMinus))).uint32());
      dwColor = pVisConfig->pDetails->GetBaseColor();
   } else
      dwColor = 0xffffff;
      // ^-- not set before rendering is finished...
   dwColor |= 0xff000000;
   return dwColor;
}

bool FOrbital::DependsOnWaveFunction() const
{
   return true;
}

void FOrbital::FlipPhase()
{
   for (size_t i = 0; i < pCoeffs.size(); ++ i)
      pCoeffs[i] *= -1;
   for (size_t i = 0; i < pIaoCoeffs.size(); ++ i)
      pIaoCoeffs[i] *= -1;
   InvalidateRenderCache();
}

// void FOrbital::InvalidateColors()
// {
//    FOrbital *pOrb = this;
//    // update the color also in the mesh, if one already exists.
//    if (pOrb->pGlMesh.get()) {
//       for (uint i = 0; i < pOrb->pGlMesh->Vertices.size(); ++ i) {
//          FBaseVertex
//             &vx = pOrb->pGlMesh->Vertices[i];
//          if (vx.iRole == 0) vx.dwColor = pOrb->pVisConfig->cIsoPlus;
//          if (vx.iRole == 1) vx.dwColor = pOrb->pVisConfig->cIsoMinus;
//       }
//       pOrb->pGlMesh->Invalidate();
//    }
// }
//
//
// void FOrbital::InvalidateRenderCache() {
//    pGlMesh = 0;
// }


// QString FOrbitalInfo::MakeDesc(size_t iOrb) const
QString FOrbitalInfo::MakeDesc(QString RefName) const
{
//    return QString("%1 [E =%2  O =%3]").arg(RefName.trimmed(),3).arg(fEnergy,8,'f',4).arg(fOcc,7,'f',4);
   QString
      s;
   QTextStream
      out(&s);
   out.setFieldAlignment(QTextStream::AlignRight);
   out << QString("%1 [").arg(RefName, 5);
   out << QString("E =%1 ").arg(fEnergy,8,'f',4);
   out << "O =" << wfi::pOrbDescFromType(Spin, fOcc);
   out << "]";
   return s;
//    return QString("%1.1 [E =%2  O =%3]").arg(iOrb+1,3).arg(fEnergy,8,'f',4).arg(fOcc,7,'f',4);
//     [E=%8.4f  O= %6.4f]
}

QString FOrbitalInfo::MakeDesc(size_t iOrb) const
{
   return MakeDesc(QString("%1.1").arg(iOrb+1));
//    return QString("%1.1 [E =%2  O =%3]").arg(iOrb+1,3).arg(fEnergy,8,'f',4).arg(fOcc,7,'f',4);
}

void FOrbital::UpdateDescFromInfo(size_t iOrb)
{
   m_Desc = info.MakeDesc(iOrb);
}

QString FOrbital::GetDesc(uint Flags) const
{
   if (Flags == DESC_Full) {
      return MakeFullDesc(0.02, ORBDESC_CompactSpace, 4);
   }
   return FBase::GetDesc(Flags);
}


QString FOrbital::MakeFullDesc(double ThrPrint, uint Flags, int nMaxAtoms) const
{
   bool
      Compact = bool(Flags & ORBDESC_CompactSpace);
   QString
      s;
   if (0 == (Flags & ORBDESC_ChargesOnly))
      s = this->GetDesc() + (Compact? "   " : "   ");
   if (Compact)
      s.replace("  ", " ");
   TArray<double>
      Charges = MakeIaoCharges(false, 2.0);
   TArray<size_t>
      Indices;
   Indices.resize(Charges.size());
   ct::ArgSort1(&Indices[0], &Charges[0], 1, Charges.size(), true); // last: reverse. Largest first.

   {
      QTextStream
         out(&s);
//       out << "  CENTERS/CHARGES:";
//       out << "   ";
      out.setRealNumberNotation(QTextStream::FixedNotation);
      double
         fChgRest = 0.;
      int
         nAtomsEmitted = 0;
      for (size_t ii = 0; ii < Indices.size(); ++ ii) {
         size_t
            iAt = Indices[ii];
         if (iAt >= pAtoms->size())
            continue; // not supposed to happen.
         double
            fAtChg = Charges[iAt];
         if (fAtChg < ThrPrint || (nMaxAtoms >= 0 && (nAtomsEmitted >= nMaxAtoms))) {
            fChgRest += fAtChg;
            continue;
         }
         ct::FAtom
            &At = (*pAtoms)[iAt];
         out.setRealNumberPrecision(0);
         if (Compact) {
            if (ii != 0)
               out << " | ";
            out << s2q(At.GetElementName()) << " " << (1+iAt) << ": ";
            out.setRealNumberPrecision(2);
            out << Charges[iAt];
         } else {
            out << "  ";
            out.setFieldWidth(2);
            out << s2q(At.GetElementName()).toUpper();
            out.setFieldWidth(3);
            out << (1+iAt);
            out.setFieldWidth(0);
            out << " ";
            out.setFieldWidth(6);
            out.setRealNumberPrecision(3);
            out << Charges[iAt];
         }

         out.setFieldWidth(0);
         nAtomsEmitted += 1;
      }
      if (fChgRest > ThrPrint) {
         out.setRealNumberPrecision(3);
         out << "   (other: " << fChgRest << ")";
      }
   }
   return s;
}


FIndexedTriangleListPtr FOrbital::MakeIsoSurface(FIsoSurfaceSettings const &IsoSurfOpt_)
{
   return ::MakeIsoSurface(this, IsoSurfOpt_);
}


double FOrbitalInfo::fChargeFactor() const
{
   return fOcc;
}

double FOrbitalInfo::fSpinFactor() const
{
   if (Spin == wfi::ORBSPIN_Alpha)
      return +fOcc;
   if (Spin == wfi::ORBSPIN_Beta)
      return -fOcc;
   return 0.;
   // ^- note: this is not quite right for active orbitals with MCSCF.. but we
   //    can't do much about it with the information we get from imported files.
}


void FOrbital::AddChargeContributions(FChargeAnalysis &Charges) const
{
   if (pIaoCoeffs.empty() || pMinBasis.get() == 0)
      return;
   double
      fCharge = info.fChargeFactor(),
      fSpin = info.fSpinFactor();
   if (fCharge == 0. && fSpin == 0.)
      return; // nothing to add (e.g., virtual orbital.)
   uint
      iShOf = 0;
   for (uint iSh = 0; iSh < pMinBasis->Shells.size(); ++ iSh) {
      ct::FBasisShell
         &Sh = pMinBasis->Shells[iSh];
      uint
         nShFn = Sh.nFn();
      double
         fChg = 0;
      for (uint iFn = 0; iFn < nShFn; ++ iFn)
         fChg += sqr(pIaoCoeffs[iFn + iShOf]);
      assert(size_t(Sh.iCenter) < pAtoms->size());
      Charges.Add(Sh.iCenter, Sh.l(), fChg*fCharge, fChg*fSpin);
      iShOf += nShFn;
   }
   assert(iShOf == pMinBasis->nFn());
}


TArray<double> FOrbital::MakeIaoCharges(bool UseOccupationNumbers, double fOccupancyOtherwise) const
{
   if (pIaoCoeffs.empty() || pMinBasis.get() == 0)
      return TArray<double>();
   TArray<double>
      Out;
   double
//       fOcc1 = 1.;
      fOcc1 = fOccupancyOtherwise;
   if (UseOccupationNumbers)
      fOcc1 = info.fOcc;
   Out.resize(pAtoms->size(), 0.);
   uint
      iShOf = 0;
   for (uint iSh = 0; iSh < pMinBasis->Shells.size(); ++ iSh) {
      ct::FBasisShell
         &Sh = pMinBasis->Shells[iSh];
      uint
         nShFn = Sh.nFn();
      double
         fChg = 0;
      for (uint iFn = 0; iFn < nShFn; ++ iFn)
         fChg += sqr(pIaoCoeffs[iFn + iShOf]);
      assert(size_t(Sh.iCenter) < pAtoms->size());
//       Out[Sh.iCenter] += info.fOcc * fChg;
      Out[Sh.iCenter] += fOcc1 * fChg;
      iShOf += nShFn;
   }
   assert(iShOf == pMinBasis->nFn());
   return Out;
}

QString FOrbital::GetVolumeType() const
{
   return QString("Orbital");
}


// make the difference in charge for IBOs along all frames.
bool MakeIboChangeCurve(TArray<float> &CurveData, uint iRow, FDocument *document)
{
   // WARNING: due to the way the document loading is done, this function is
   // invoking at least O(N^2) cost (and maybe O(N^3)) when loading lots of frames,
   // since it is called after every frame.
   // Should be fixed... might get expensive?
   FOrbital
      *pRefOrb = dynamic_cast<FOrbital*>(document->GetRowCol(iRow, 0, false));
   if (!pRefOrb)
      // data set is not there or is not an orbital.
      return false;
   uint
      nCols = document->GetNumFrames();
   CurveData.clear();
   CurveData.resize(nCols, 0.);
   TArray<double>
      RefChg = pRefOrb->MakeIaoCharges(g_UseOccupationNumbersForIboChangeCurves, 1.0);
//    PrintChargeArray(std::cout, "ref charge", RefChg);
   if (RefChg.empty())
      return false; // IAO coeffs not made.
   for (uint iCol = 0; iCol < nCols; ++ iCol) {
      FOrbital
         *pCurOrb = dynamic_cast<FOrbital*>(document->GetRowCol(iRow, iCol, false));
      TArray<double>
         CurChg = pCurOrb->MakeIaoCharges(g_UseOccupationNumbersForIboChangeCurves, 1.0);
//       PrintChargeArray(std::cout, "cur charge", RefChg);
      double f = 0;
      assert(RefChg.size() == CurChg.size());
      for (uint iAt = 0; iAt < RefChg.size(); ++ iAt)
         f += sqr(RefChg[iAt] - CurChg[iAt]);
      CurveData[iCol] = 1.0 * std::sqrt(f);
   }
   return true;
}


size_t FOrbital::GetMaxStackDataSize() const
{
   return omp_get_max_threads() * (sizeof(double)*4*128 + 24) * pBasisSet->nFn();
   // ^- some per-thread temporaries in MakeGridValues (4 deriv components, 128 grid block, max number of basis fn).
}


void FOrbital::MakeGridValues(double *pOut, ct::FMatrixView Grid, uint GridDxOrder, int Role, ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(Role);
   ct::FMatrixView
      CoeffMatrix = ct::FMatrixView(const_cast<double*>(&this->pCoeffs[0]), this->pCoeffs.size(), 1);
   return ct::MakeGridValues(pOut, Grid, CoeffMatrix, GridDxOrder, &*this->pBasisSet, Mem);
}


ct::FAtomSetPtr FOrbital::GetSupportAtoms(ct::FMemoryStack &Mem) const
{
   IR_SUPPRESS_UNUSED_WARNING(Mem);
   using ct::FVector3;
   if (!HaveMoments)
      // no dipole/quadrupole moments stored. return list of all atoms.
      // (note: in principle we could also estimate the support atoms based on
      //  partial charges, if we have them).
      return pAtoms;
   else {
      // have first and 2nd moments of orbital. Return only atoms which
      // lie within a certain radius of the orbital centroid, depending
      // on the orbital 2nd moments (around its centroid).
      ct::FAtomSetPtr
         pSupportAtoms(new ct::FAtomSet);
      FSupportBox
         Box(&vDipMom[0], &mQuadMom[0][0], g_fSupportBox_EnlargeRms);
      for (size_t iAtom = 0; iAtom < pAtoms->size(); ++ iAtom) {
         bool Keep = true;
         for (size_t iAxis = 0; iAxis < 3; ++ iAxis) {
            FVector3 c(Box.vCenter[0], Box.vCenter[1], Box.vCenter[2]);
            FVector3 vAxisPos = (*pAtoms)[iAtom].vPos - c;
            double f = ct::Dot(&vAxisPos[0], &Box.vAxes[iAxis][0], 3);
            if (f < Box.fMinExtend[iAxis] || f > Box.fMaxExtend[iAxis])
               Keep = false;
         }
         if (Keep) {
//             pSupportAtoms->Atoms.push_back((*pAtoms)[iAtom]);
            pSupportAtoms->AddAtom((*pAtoms)[iAtom]);
         }
      }
      if (pSupportAtoms->empty())
         // in case there are no atoms at all, most likely something went wrong.
         // return the entire set in this case.
         return pAtoms;
      else
         return pSupportAtoms;
   }

}


double FOrbital::GetDataWeight(double fValue_) const {
   return sqr(fValue_);
}

QString FOrbital::GetDataUnit() const
{
   return "a0^{-3/2}";
}
