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
#include <math.h> // for sqrt.

#include "CxIo.h"
#include "CtMatrix.h"
#include "CxColor.h"

#include "IvDocument.h"
#include "IvScript.h"
#include "IvIrc.h"
// #include "IvOrbitalFile.h"

#include "CxAtomData.h"
#include "CtBasisLibrary.h"
#include "CxPodArray.h"
using namespace ct;


FArcLength::FArcLength()
   : m_Value(-1.), m_ArcType(FArcLengthOptions::ATOMWEIGHT_Unity)
{
}


FArcLength::FArcLength(double fValue_, FArcLengthOptions Type_)
   : m_Value(fValue_), m_ArcType(Type_)
{
}

TArray<float> ExtractFloats(TArray<FArcLength> const &ArcLengths)
{
   TArray<float> r;
   r.resize(ArcLengths.size());
   for (size_t i = 0; i < ArcLengths.size(); ++ i)
      r[i] = float(ArcLengths[i].Value());
   return r;
}

FArcLengthOptions::FArcLengthOptions(FWeightType WeightType_)
   : WeightType(WeightType_)
{
}


double FArcLengthOptions::GetAtomWeight(int iElement, double fMassOverride, double fNuclearChargeOverride) const
{
   switch (WeightType) {
      case ATOMWEIGHT_Unity:
         return 1.;
      case ATOMWEIGHT_IrcAvgMass:
         if (fMassOverride == -1.)
            return GetAtomicMass(iElement, ATMASS_StandardAtomicWeight);
         else
            return fMassOverride;
      case ATOMWEIGHT_IrcIsotopeMass:
         if (fMassOverride == -1.)
            return GetAtomicMass(iElement, ATMASS_MostCommonIsotope);
         else
            return fMassOverride;
      case ATOMWEIGHT_Charge:
         if (fNuclearChargeOverride == -1.)
            return double(iElement);
         else
            return fNuclearChargeOverride;
      default:
         throw std::runtime_error("FArcLengthOptions::GetAtomWeight: encountered unknown atomic weight type in arc length computation.");
   }
}

QString FArcLengthOptions::GetArcUnit() const
{
   return QString();
}

QString FArcLengthOptions::GetArcCaption() const
{
   switch (WeightType) {
      case ATOMWEIGHT_Unity:
         return "Arc Length (bohr)";
      case ATOMWEIGHT_IrcAvgMass:
      case ATOMWEIGHT_IrcIsotopeMass:
         return "Reaction Coordinate (bohr$\\cdot$amu$^{1/2}$)";
      case ATOMWEIGHT_Charge:
         return "Reaction Coordinate (bohr$\\cdot$e$^{1/2}$)";
      case ATOMWEIGHT_FrameIndex:
         return "Frame index (0-based)";
      default:
         return "Arc Length (?)";
   }
}


// pAtoms: if given: take from there; otherwise: take from current frame.
double FArcLengthOptions::GetAtomWeight(int iAt, FDocument *document, FAtomSet *pAtoms) const
{
   if ((document->AtomFlags(iAt) & ATOM_NoAlign) != 0)
      // ignore this atom in setting up the transformation & coordinates.
      return 0.;
   FElementOptions const
      *pAtOptions = document->pElementOptions(iAt, pAtoms);
   double
      fWeight = GetAtomWeight(pAtOptions->iElement(), pAtOptions->GetMassOverride(), pAtOptions->GetNuclearChargeOverride());
//    IvEmit("  !GetAtomWeight(%1) = %2 ", iAt, fWeight);
   return fWeight;
}



// root mean square deviation between two vectors (i.e., 2-norm of their difference)
static double Rmsd(double const *pA, double const *pB, size_t N) {
   double RmsdSq = 0;
   for (size_t i = 0; i < N; ++ i)
      RmsdSq += sqr(pA[i] - pB[i]);
   return sqrt(RmsdSq);
}


// estimate the IRC arc positions for all loaded frames. Optionally do not do the mass weighting.
void MakeIrcArcLengths(TArray<FArcLength> &ArcLengths, FDocument *document, FArcLengthOptions ArcOptions)
{
   ArcLengths.clear();

   // get a pointer to the initial geometry.
   IFrame
      *pFrame0 = document->GetFrame(0, false); // false: don't raise error if frame #0 does not exist
   if (pFrame0 == 0)
      // no frames.
      return;
   FAtomSet
      *pAtoms0 = pFrame0->pGetAtoms();
   if (pAtoms0 == 0)
      return;

   // make a list of the atomic masses
   uint
      nAt = pAtoms0->size();
   TArray<double>
      Masses(nAt, 1.);
   for (int iAt = 0; size_t(iAt) < nAt; ++ iAt)
      Masses[size_t(iAt)] = ArcOptions.GetAtomWeight(iAt, document, pAtoms0);

   // assemble a matrix of the mass weighted coordinates for all frames.
   // The numbers and types of atoms in all frames must be equal.
   uint
      nFrames = (unsigned)document->GetNumFrames();
   FMemoryStack2
      Mem(3 * nAt * nFrames * sizeof(double) + 2000000);
   FStackMatrix
      Mwcs(3 * nAt, nFrames, &Mem);
   bool
      AssignmentInvalid = false;
   for (uint iFrame = 0; iFrame < nFrames; ++ iFrame) {
      IFrame
         *pFrame = document->GetFrame(iFrame);
      assert(pFrame != 0);
      FAtomSet
         *pAtoms = pFrame->pGetAtoms();
      // check if all other geometries have the same number and types of atoms
      if (pAtoms == 0 || pAtoms->size() != nAt) {
         IvNotify(NOTIFY_Warning, "IrcArcLength: Number of atoms not consistent across different frames.");
         AssignmentInvalid = true;
         break;
      }
      for (uint iAt = 0; iAt < nAt; ++ iAt) {
         if ((*pAtoms)[iAt].iElement != (*pAtoms0)[iAt].iElement) {
            IvNotify(NOTIFY_Warning, "IrcArcLength: Type of atoms not consistent across different frames.");
            AssignmentInvalid = true;
            break;
         }
         for (uint ixyz = 0; ixyz < 3; ++ ixyz)
            Mwcs(3*iAt + ixyz, iFrame) = sqrt(Masses[iAt]) * (*pAtoms)[iAt].vPos[ixyz];
      }
   }
   if (AssignmentInvalid) {
      ArcLengths.clear();
   } else {
      // approximate the IRC arc lengths via the differences of the 2-norm of the mass-weighted
      // coordinates (note: if we have gradients we could calculate them more accurately. should be
      // done at some point. Even without gradients we could at least put a spline through the points etc.).
      ArcLengths.clear();
      ArcLengths.reserve(nFrames);
      double ArcLength = 0.;
      ArcLengths.push_back(FArcLength(ArcLength, ArcOptions));
      for (uint iFrame = 1; iFrame < nFrames; ++ iFrame) {
         double Delta = Rmsd(&Mwcs(0, iFrame), &Mwcs(0, iFrame-1), Mwcs.nRows);
         ArcLength += Delta;
         ArcLengths.push_back(FArcLength(ArcLength, ArcOptions));
      }
      assert(ArcLengths.size() == nFrames);
      // that's it already.
   }
}
