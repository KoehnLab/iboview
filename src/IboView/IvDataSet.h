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

#ifndef IV_DATASET_H
#define IV_DATASET_H

#include "Iv.h"

#include <map>
#include <set>
#include <vector>
#include <list>

#include <QString>
#include <QObject>
#include <QStringList>
// #include <QAbstractTableModel>
// #include <QAction>
// #include <QExplicitlySharedDataPointer>
// #include <QSharedData>


#include "CtAtomSet.h"
// #include "CtBasisSet.h"
#include "CxPodArray.h"

// #include "CxColor.h"
#include "IvMesh.h"
#include "IvGl.h"
// #include "IvIsoSurface.h"
#include "IvDataOptions.h"
// #include "IvAnalysis.h"
// #include "CtMatrix.h"
// #include "IvLog.h"
// #include "IvTables.h"
// #include "IvIrc.h"
// #include "CxVec3.h"
// #include "CtWfi.h"


// #include "CtDftGrid_ivb.h" // FIXME: merge back with main version once done.
// #include "CtDftGrid.h"


class FFrame;
class FView3d;
class FDocument;


// using wfi::FOrbitalSpin;
// using wfi::FWfType;

using ct::FAtomSetPtr;
using ct::FBasisSetPtr;
using ct::TArray;
using ct::FIntrusivePtrDest;
namespace ct {
   struct FHfOptions;
   struct FWfDecl;
}


// TODO: Decide QObject vs FIntrusivePtrDest
// - Would it be a good or bad idea to make all the FDataSet classes QObjects,
//   instead of FIntrusivePtrDest, as they are now?
//
// - As is, in principle my data sets can be copied (even if associated resource
//   management, say, GL stuff, might make it less trivial in practice).
//   Not sure if this freedom is used anywhere at the moment.
//   As QObjects, copying will become impossible.
//
// - Making them QObjects would, at least in principle, allow exposing them
//   to the scripting interface and UI signals directly.
//
//   + A script interface could also be achieved by creating proxy object
//     accessors only (which are QObjects), rather than converting the actual
//     data set class. I guess these would have to be dynamically created and
//     put as childs to the script engine, if that is done. How to do UI
//     interfaces with clean ownership would be more complicated.
//
// - Making them QObjects would simplify adding QObjects as children
//   (e.g., as the FFreeObject lists I used for FFreeObjectSet)
//
// - Leaving them explicitly managed, as they are, yields a more transparent
//   memory/resource model. One does not need to really understand QObject
//   hierarchies to see how and when the objects will be deleted.

struct FDataSet : public FIntrusivePtrDest
{
   explicit FDataSet(QString const &Desc_, FAtomSetPtr pAtoms_ = 0, FDocument *pDocument_ = 0);

   bool
      Active;
   FAtomSetPtr
      pAtoms;
   enum FDescFlags {
      DESC_Full = 0x01,
      DESC_Compact = 0x02
   };

//    FRenderCachePtr
//       pRenderCache; // may be 0. may be set to 0 at any point in time.
   virtual QString GetType() const = 0;
   virtual QString GetDesc(uint Flags=0) const; // default: return m_Desc

   // rebuild the visual representation of the object on next rendering.
   virtual void InvalidateRenderCache();
   virtual void BuildRenderCache(FView3d *pView3d);

   virtual uint32_t GetBaseColor() const;

   virtual bool DependsOnWaveFunction() const;

   // return pointer to its parent object.
   FDocument *GetDocument() { return m_pDocument; };
protected:
   QString m_Desc;
   FDocument *m_pDocument;
// private:
};
typedef ct::TIntrusivePtr<FDataSet>
   FDataSetPtr;




enum FBondFlags {
   BOND_Partial = 0x01,
   BOND_Grey = 0x02 // make the bond grey instead of the usual half-bond with their attached atom colors.
};

struct FBondLine {
   int
      iAt, jAt;
   uint
      Flags;
   // additional visualization data.
   double
      fBondOrder;
   FBondLine() {};
   FBondLine(int iAt_, int jAt_, uint Flags_, double fBondOrder_ = 1.);
};

enum FAtomFlag {
   // do not show the atom in renderings
   ATOM_Hidden = 0x01,
   // if aligning the molecule in space, ignore this atom when determining
   // the aligning transformation.
   ATOM_NoAlign = 0x02,
   // atom is part of an currently active selection
   ATOM_Selected = 0x04
};

// some additional data used in the visualization of potential bonds between atoms.
struct FBondVisualInfo {
   int
      iAt, jAt;
   TVector3<float>
      vAtPosI, vAtPosJ, vNorm, vTanU, vTanV;
   float
      fDist,
      fDistScaled;
   explicit FBondVisualInfo(int iAt_, int jAt_, ct::FAtomSet const &Atoms);
   explicit FBondVisualInfo(int iAt_, int jAt_); // this one used as dummy only.
};

namespace ct {
   struct FBondOrderAnalysis;
}

// Describes a molecular geometry (i.e., positions of atoms, bonds, etc.)
struct FGeometry : public FDataSet
{
   explicit FGeometry(QString const &Desc_, FAtomSetPtr pAtoms_, FDocument *pDocument_);
   QString GetType() const; // override
   // add a dotted bond to indicate the formation/breaking of a bond between the indexed atoms.
   void AddBond(FBondLine const &bl);
   void AddBond(int iAt, int jAt, QString const &Flags);
   void DeleteBond(int iAt, int jAt, bool ComplainIfNotThere=true);
   void FindBondLines(FDocument *pDocument, ct::FBondOrderAnalysis *pBondOrderAnalysis = 0);
   void UpdateVisualInfo();
public:
   std::vector<FBondLine>
      m_BondLines;
   typedef std::vector<FBondVisualInfo>
      FBondVisualInfoList;
   FBondVisualInfoList
      m_BondVisualInfo;
   FBondVisualInfo const *FindBondVisualInfo(int iAt, int jAt) const;
protected:
   void FixMinBasis();
};




// https://doc.qt.io/qt-5/implicit-sharing.html
// https://doc.qt.io/qt-5/qexplicitlyshareddatapointer.html
// typedef QSharedPointer<FFreeObject>
//    FFreeObjectPtr;
// typedef QList<FFreeObjectPtr>
//    FFreeObjectPtrList;
typedef QList<FFreeObjectPtr> // <-- those are `QObject`s now, linked to shared data blocks.
   FFreeObjectPtrList;

// Sets of additional free objects (not attached to a specific geometry) in 3D space.
// Used to represent user-defined labels, lines, planes, meshes, etc.
struct FFreeObjectSet : public FDataSet
{
   explicit FFreeObjectSet(QString const &Desc_, FAtomSetPtr pAtoms_, FDocument *pDocument_);
   ~FFreeObjectSet();
   QString GetType() const; // override
   void UpdateVisualInfo();
// protected:

   // adds *p to the objects controlled by *this. Assumed to take ownership.
   void Take(FFreeObjectPtr p);
   // adds a new FFreeObject to *this, however, its actual indirect data records
   // will be linked to the shared data section of *p
   void AddLink(FFreeObject *p);
   // adds a new FFreeObject to *this, with data section cloned from *p.
   void AddClone(FFreeObject const *p);
protected:
   FFreeObjectPtrList
      m_Objects;
public:
//    typedef FFreeObjectPtrList::iterator iterator;
//    typedef FFreeObjectPtrList::const_iterator const_iterator;
//    iterator begin() { return m_Objects.begin(); }
//    iterator end() { return m_Objects.end(); }
   bool _IndexValidQ(size_t i) const;
   void _AssertIndexIsValid(size_t i) const;
   size_t size() const { return m_Objects.size(); }
   bool empty() const { return m_Objects.empty(); }
   void clear() { m_Objects.clear(); }
   FFreeObject const &operator[] (size_t i) const;
   FFreeObject &operator[] (size_t i);

   IV_IMPLEMENT_ITERATOR_IndexInBracket(FFreeObjectSet,const_iterator,const)
   IV_IMPLEMENT_ITERATOR_IndexInBracket(FFreeObjectSet,iterator,)
protected:
   FFreeObject *_MakeInvalidIndexObject(size_t i) const;
};

typedef ct::TIntrusivePtr<FFreeObjectSet>
   FFreeObjectSetPtr;





#endif // IV_DATASET_H
