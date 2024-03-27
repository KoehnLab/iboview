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
// #include "CxColor.h"

#include "IvDataSet.h"
#include "IvAnalysis.h" // for FBondOrderAnalysis. Used to derive bond lines in FGeometry
// #include "CxPhysicalUnits.h"
// #include "CxOpenMpProxy.h"
// // #include "CxXyzFrameIo.h"
//
// #include "IvDocument.h"
// #include "IvSettings.h"
// #include "IvFileConvert.h"
//
// // #include "CtInt1e.h"
// // #include "IrAmrr.h"
// #include "CtMatrix.h"
// #include "CtRhf.h"

#include "CtBasisLibrary.h"

// #include "IvDocument.h" // only for document->pDocument


using namespace ct;


FDataSet::FDataSet(QString const &Desc_, FAtomSetPtr pAtoms_, FDocument *pParent)
   : pAtoms(pAtoms_), m_Desc(Desc_), m_pDocument(pParent)
{
   Active = false;
}

QString FDataSet::GetDesc(uint) const { return m_Desc; };


void FDataSet::InvalidateRenderCache()
{
}

void FDataSet::BuildRenderCache(FView3d * /*pView3d*/)
{
}


FGeometry::FGeometry(QString const &Desc_, FAtomSetPtr pAtoms_, FDocument *pDocument_)
   : FDataSet(Desc_, pAtoms_, pDocument_)
{
   FindBondLines(pDocument_);
   UpdateVisualInfo();
   FixMinBasis();
}

void FGeometry::FixMinBasis()
{
   for (unsigned iAt = 0; iAt < pAtoms->size(); ++ iAt) {
      ct::FAtom
         &At = (*pAtoms)[iAt];
      // is some sort of non-trivial minimal basis already set explictly?
      // (e.g., during orbital file loading for semi-empirical methods)
      {
         std::string const
            &sGuessBasis = At.BasisDesc[BASIS_Guess];
         if (!(sGuessBasis == "MINAO" || sGuessBasis == "")) {
            // if yes, leave it alone.
            if (sGuessBasis != "USE-ORB-BASIS" && !ct::g_BasisSetLibrary.HaveBasis(At.iElement, sGuessBasis)) {
               IV_NOTIFY3(NOTIFY_Warning, "Atom %1 (%2) has minimal basis '%3' explicitly declared, but it could not be instanciated.", iAt, s2q(At.ElementName()), s2q(sGuessBasis));
            }
            continue;
         }
      }

      if (ct::g_BasisSetLibrary.HaveBasis(At.iElement, "MINAO"))
         At.BasisDesc[BASIS_Guess] = "MINAO";
      else if (ct::g_BasisSetLibrary.HaveBasis(At.iElement, "MINAO-PP"))
         At.BasisDesc[BASIS_Guess] = "MINAO-PP";
      else {
         IV_NOTIFY2(NOTIFY_Warning, "Failed to find a minimal basis for atom %1 (%2)", iAt, s2q(At.ElementName()));
         At.BasisDesc[BASIS_Guess] = "MINAO";
      }
   }
}


void FGeometry::AddBond(FBondLine const &bl)
{
   if (bl.iAt >= (int)pAtoms->size() || bl.jAt >= (int)pAtoms->size()) {
//       std::cout << fmt::format("!Error: Ignoring bond {:3}--{:3}. Have only {} atoms!", (1+bl.iAt), (1+bl.jAt), pAtoms->size()) << std::endl;
      IV_NOTIFY3(NOTIFY_Warning, "Ignoring bond %1--%2. Have only %3 atoms!", 1+bl.iAt, 1+bl.jAt, (int)pAtoms->size());
      return;
   }
   DeleteBond(bl.iAt+1, bl.jAt+1, false);
   m_BondLines.push_back(bl);
}

void FGeometry::DeleteBond(int iAt, int jAt, bool ComplainIfNotThere)
{
   iAt -= 1; // convert 1-based (external) to 0-based (internal)
   jAt -= 1;
   bool
      FoundSomething = false;
   for (uint iBond = 0; iBond < m_BondLines.size(); ++ iBond) {
      FBondLine &bl = m_BondLines[iBond];
      if ( (bl.iAt == iAt && bl.jAt == jAt) ||
           (bl.iAt == jAt && bl.jAt == iAt) ) {
         m_BondLines.erase(m_BondLines.begin() + iBond);
         -- iBond;
         FoundSomething = true;
      }
   };
   if (ComplainIfNotThere && !FoundSomething)
      IV_NOTIFY(NOTIFY_Warning, QString("No bond %1--%2 present which could be deleted.").arg(iAt).arg(jAt));
}

void FGeometry::AddBond(int iAt, int jAt, QString const &sFlags)
{
   uint Flags = 0;
   QStringList FlagList = sFlags.split("|", QString::SkipEmptyParts);
   double
      BondOrder = 1.;
   foreach(QString const &Flag, FlagList) {
      if (Flag == "gray" || Flag == "grey")
         Flags |= BOND_Grey;
      else if (Flag == "dotted")
         Flags |= BOND_Partial;
      else if (Flag.startsWith("bo(") && Flag.endsWith(")")) {
         bool ok;
         BondOrder = Flag.mid(3, Flag.size()-4).toDouble(&ok);
         if (!ok) {
            BondOrder = 1.;
            IV_NOTIFY1(NOTIFY_Warning, "Bond order '%1' not recognized.", Flag);
         }
      } else {
   //       std::cerr << "!WARNING: bond flag '" << q2s(sFlags) << "' not recognized." << std::endl;
         IV_NOTIFY1(NOTIFY_Warning, "Bond flag '%1' not recognized.", Flag);
      }
   }
   if ((Flags & BOND_Partial) && BondOrder == 1.)
      BondOrder = 0.5;
   AddBond(FBondLine(iAt-1, jAt-1, Flags, BondOrder));
}


// forward (IvDocument.cpp). proxy to pDocument->pElementOptions(iAt, pAtoms)
FElementOptions *pElementOptions(FDocument *pDocument, int iAt, ct::FAtomSet *pAtoms);


void FGeometry::FindBondLines(FDocument *pDocument, FBondOrderAnalysis *pBondOrderAnalysis)
{
   m_BondLines.clear();
   for (unsigned  iAt0 = 0; iAt0 < pAtoms->size(); ++ iAt0) {
      FAtom const
         &At0 = (*pAtoms)[iAt0];
      for (unsigned iAt1 = iAt0+1; iAt1 < pAtoms->size(); ++ iAt1) {
         FAtom const
            &At1 = (*pAtoms)[iAt1];
         if (pBondOrderAnalysis) {
            // use computed bond orders for deciding where the bonds are, and how many.
            FBondOrderAnalysis::FValue const
               *pBoValue = pBondOrderAnalysis->FindValue(iAt0, iAt1);
            if (!pBoValue)
               continue;
            double
//                fBondOrder = pBoValue->fBondOrders[FBondOrderAnalysis::BONDORDER_Renorm1Iao];
               fBondOrder = pBoValue->fBondOrders[FBondOrderAnalysis::BONDORDER_WibergIao];
               // ^-- hm... wasn't this configurable at some point?
            AddBond(FBondLine(iAt0, iAt1, 0, fBondOrder));
         } else {
            // find bond lines based on geometric distances.
            FElementOptions
               *pEo0 = pElementOptions(pDocument, iAt0, &*pAtoms),
               *pEo1 = pElementOptions(pDocument, iAt1, &*pAtoms);
            double
               r0 = pEo0->GetCovalentRadius(),
               r1 = pEo1->GetCovalentRadius(),
               rij = Dist(At0.vPos, At1.vPos);
            if (rij > 0.5 * (pEo0->GetBondRadiusFactor() + pEo1->GetBondRadiusFactor()) * (r0 + r1))
               continue;
            AddBond(FBondLine(iAt0, iAt1, 0));
         }
      }
   }
}
#ifdef INCLUDE_ABANDONED
//             FElementOptions
//                *pEo0 = pDocument->pElementOptions(iAt0, &*pAtoms),
//                *pEo1 = pDocument->pElementOptions(iAt1, &*pAtoms);
#endif // INCLUDE_ABANDONED



FBondLine::FBondLine(int iAt_, int jAt_, uint Flags_, double fBondOrder_)
   : iAt(iAt_), jAt(jAt_), Flags(Flags_), fBondOrder(fBondOrder_)
{
}

FBondVisualInfo::FBondVisualInfo(int iAt_, int jAt_)
   : iAt(iAt_), jAt(jAt_)
{
   if (jAt < iAt)
      std::swap(iAt,jAt);
}

typedef TVector3<float>
   FVec3f;
void MakeOrthDirections(FVec3f &vOrth0, FVec3f &vOrth1, FVec3f const &vIn);

static double fWeightedDistSq(int iAt, int jAt, FAtomSet const &Atoms) {
   assert(size_t(iAt) < Atoms.size() && size_t(jAt) < Atoms.size());
   double
      fCovI = ct::GetCovalentRadius(Atoms[iAt].iElement),
      fCovJ = ct::GetCovalentRadius(Atoms[jAt].iElement);

   // give precedence to some elements typically occuring in hetero cycles.
   int iElemJ = Atoms[jAt].iElement;
   double fPiEnh = 1.;
   if ((iElemJ >= 5 && iElemJ <= 8) || (iElemJ >= 15 && iElemJ <= 16))
      fPiEnh = 1.5;

   return DistSq(Atoms[iAt].vPos, Atoms[jAt].vPos)/sqr(.5*(fCovI + fCovJ)*fPiEnh);
}

FBondVisualInfo::FBondVisualInfo(int iAt_, int jAt_, FAtomSet const &Atoms)
  : iAt(iAt_), jAt(jAt_)
{
   if (jAt < iAt)
      std::swap(iAt,jAt);
   if (iAt == jAt)
      throw std::runtime_error(fmt::format("FBondVisualInfo: cannot make a bond between atom {} and itself.", iAt));
   if (iAt < 0 || size_t(iAt) >= Atoms.size())
      throw std::runtime_error(fmt::format("FBondVisualInfo: explicit reference to non-existing atom {} in frame {}", iAt, Atoms.GetName()));
   if (jAt < 0 || size_t(jAt) >= Atoms.size())
      throw std::runtime_error(fmt::format("FBondVisualInfo: explicit reference to non-existing atom {} in frame {}", jAt, Atoms.GetName()));
   vAtPosI = FVec3f(Atoms[iAt].vPos);
   vAtPosJ = FVec3f(Atoms[jAt].vPos);
   vNorm = vAtPosJ - vAtPosI;
   fDist = Length(vNorm);
   if (fDist != 0.)
      vNorm /= fDist;

   fDistScaled = fDist / (ct::GetCovalentRadius(Atoms[iAt].iElement) * ct::GetCovalentRadius(Atoms[jAt].iElement));

   // make an initial set of orthonormal directions.
   // Better: align them with other atoms, such that we get predictable placement of
   // multi-bonds which are displaced along the vTanU/vTanV directions.
   MakeOrthDirections(vTanU, vTanV, vNorm);
   std::swap(vTanU, vTanV);
   vTanU *= -1.;

   if (1) {
      // find a direction which more or less aligns double bonds with their plane.
      // Since we do not actually have any electronic structure info at this point, we just
      // search for the closest non-colinear atom(s) to form the sum.
      // A second problem is the U alignment: In benzene rings, for example, we want the
      // partial bonds shown on the INSIDE, so some additional hackery is needed.
      // Should at some point be updated with electronic structure info, once available.
      size_t
         nAt = Atoms.size();
      TArray<std::pair<double, size_t> >
         AtomDist;
      AtomDist.resize(nAt);
      for (size_t kAt = 0; kAt < nAt; ++ kAt) {
//          AtomDist[kAt].first = std::min(DistSq(Atoms[kAt].vPos, Atoms[iAt].vPos), DistSq(Atoms[kAt].vPos, Atoms[jAt].vPos));
         AtomDist[kAt].first = std::min(fWeightedDistSq(iAt,kAt,Atoms), fWeightedDistSq(jAt,kAt,Atoms));
         AtomDist[kAt].second = kAt;
      }
      std::sort(AtomDist.begin(), AtomDist.end());

      double const
         ThrColinearAngle = 10. * (M_PI/180.),
         CosThrColinearAngle = std::cos(ThrColinearAngle);
      for (size_t kkAt = 0; kkAt < nAt; ++ kkAt) {
         size_t
            kAt = AtomDist[kkAt].second;

         // loop still includes iAt and jAt itself. Skip them.
         FVec3f
            vAtPosK = FVec3f(Atoms[kAt].vPos);
         if (DistSq(vAtPosK, vAtPosI) < sqr(1e-3) || DistSq(vAtPosK, vAtPosJ) < sqr(1e-3))
            continue;

         FVec3f
            vNormJK = Normalized(vAtPosK - vAtPosI);
         // check if AtK lies on a line between I and J. If yes, it will not be useful for
         // determining the plane alignment.
         double
            fCosPhi = Dot(vNorm, vNormJK);
         if (fCosPhi < -CosThrColinearAngle || fCosPhi > CosThrColinearAngle)
            // well.. that would be near colinear.
            continue;
         vTanV = Normalized(Cross(vNorm, vNormJK));
         vTanU = Cross(vTanV, vNorm);
         break;
      }
   }

}


// predicate for finding bond visualization data belonging to a given pair of atoms (iAt,jAt).
struct FBondVisualInfoAtomSortPred {
   bool operator () (FBondVisualInfo const &A, FBondVisualInfo const &B) const {
      if (A.iAt < B.iAt) return true;
      if (B.iAt < A.iAt) return false;
      return A.jAt < B.jAt;
   }
};

void FGeometry::UpdateVisualInfo()
{
   m_BondVisualInfo.clear();
   size_t
      nAt = pAtoms->size();
   m_BondVisualInfo.reserve((nAt*(nAt-1))/2);
   for (int iAt = 0; size_t(iAt) < nAt; ++ iAt)
      for (int jAt = iAt+1; size_t(jAt) < nAt; ++ jAt)
         m_BondVisualInfo.push_back(FBondVisualInfo(iAt, jAt, *pAtoms));

   // note: should be in correct order already, but let us make sure anyway.
   std::sort(m_BondVisualInfo.begin(), m_BondVisualInfo.end(), FBondVisualInfoAtomSortPred());
}

FBondVisualInfo const *FGeometry::FindBondVisualInfo(int iAt, int jAt) const
{
   FBondVisualInfo
      DummyRef(iAt, jAt);
   FBondVisualInfoList::const_iterator
      itVisualInfo = std::lower_bound(m_BondVisualInfo.begin(), m_BondVisualInfo.end(), DummyRef, FBondVisualInfoAtomSortPred());
   if (itVisualInfo == m_BondVisualInfo.end())
      return 0;
   if (itVisualInfo->iAt == DummyRef.iAt && itVisualInfo->jAt == DummyRef.jAt)
      return &*itVisualInfo;
   return 0;
}

QString FGeometry::GetType() const { return "G"; }; // geometry






FFreeObjectSet::FFreeObjectSet(QString const &Desc_, FAtomSetPtr pAtoms_, FDocument *pDocument_)
   : FDataSet(Desc_, pAtoms_, pDocument_)
{
}


FFreeObjectSet::~FFreeObjectSet()
{
}


QString FFreeObjectSet::GetType() const { return "g"; }; // geometry supplement: additional free objects


void FFreeObjectSet::Take(FFreeObjectPtr p) {
   m_Objects.push_back(p);
}


void FFreeObjectSet::AddLink(FFreeObject *p) {
   QObject *parent = 0;
   m_Objects.push_back(p->newLinkedObject(parent));
}


void FFreeObjectSet::AddClone(FFreeObject const *p) {
   QObject *parent = 0;
   m_Objects.push_back(p->newClonedObject(parent));
}


void FFreeObjectSet::UpdateVisualInfo()
{
}

bool FFreeObjectSet::_IndexValidQ(size_t i) const
{
   return i < size();
}

void FFreeObjectSet::_AssertIndexIsValid(size_t i) const
{
   if (!_IndexValidQ(i)) {
      throw std::out_of_range(fmt::format("FFreeObjectSet: attempt to access element at invalid index #{} (not 0 <= i < {})", i, this->size()));
   }
}


FFreeObject *FFreeObjectSet::_MakeInvalidIndexObject(size_t i) const {
   // this is a somewhat more ugly and elaborate construction to return
   // something not entirely likely to crash if scripts attempt to access a non-
   // existent free object. Normally one would just raise a
   // runtime_error/out_of_range, but I guess in this particular case it might
   // be a better choice to not crash the application, but rather just complain
   // and ignore to invalid command.
   if (!_IndexValidQ(i)) {
      IvWarn(QString("FFreeObjectSet: attempt to access element at invalid index #%1 (not 0 <= i < %2)").arg(i).arg(this->size()));
   }
   static FInvalidFreeObject invalidObject(0);
   return &invalidObject;
};

FFreeObject const &FFreeObjectSet::operator[] (size_t i) const {
   if (!_IndexValidQ(i))
      return *_MakeInvalidIndexObject(i);
   _AssertIndexIsValid(i);
   return *m_Objects.at(i);
}
FFreeObject &FFreeObjectSet::operator[] (size_t i) {
   if (!_IndexValidQ(i))
      return *_MakeInvalidIndexObject(i);
   _AssertIndexIsValid(i);
   return *m_Objects[i];
}







// FRenderCache::~FRenderCache()
// {
// };


// void FDataSet::InvalidateRenderCache()
// {
//    this->pRenderCache = 0;
// };

