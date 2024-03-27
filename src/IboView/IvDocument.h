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

#ifndef IV_DOCUMENT_H
#define IV_DOCUMENT_H

#include "Iv.h"

#include <map>
#include <set>
#include <vector>
#include <list>

#include <QString>
#include <QAbstractTableModel>
#include <QAction>
// #include <QExplicitlySharedDataPointer>
// #include <QSharedData>


#include "CxPodArray.h"
#include "CtMatrix.h"
#include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "CtWfi.h"
// #include "CtDftGrid_ivb.h" // FIXME: merge back with main version once done.
// #include "CtDftGrid.h"

// #include "CxColor.h"
// #include "IvMesh.h"
// #include "IvGl.h"
// #include "IvIsoSurface.h"
#include "IvAnalysis.h"
#include "IvLog.h"
#include "IvTables.h"
#include "IvIrc.h"

#include "IvDataSet.h"
#include "IvDataOptions.h"
#include "IvIsoSurface.h"
class FOrbital;
class FVolumePropertyInfo;



using wfi::FOrbitalSpin;
using wfi::FWfType;

using ct::FAtomSetPtr;
using ct::FBasisSetPtr;
using ct::TArray;
using ct::FIntrusivePtrDest;
// namespace ct {
//    struct FHfOptions;
//    struct FWfDecl;
// }




// // abstract base class for storing intermediate data for rendering,
// // which is associated with data sets
// struct FRenderCache : public ct::FIntrusivePtrDest
// {
//    virtual ~FRenderCache();
// private:
//    FRenderCache(FRenderCache const &); // not implemented
//    void operator = (FRenderCache const &); // not implemented
// };
//
// typedef ct::TIntrusivePtr<FRenderCache>
//    FRenderCachePtr;

class FView3d;
class FDocument;

// struct FFragmentationDef;
// struct FFrameWfData;

typedef std::vector<FDataSetPtr>
   FDataSetList;
struct FFrame : public QObject, FIntrusivePtrDest
{
   explicit FFrame(QString Desc_, FDocument *pDocument_);

   enum FOrbitalMatrixFlags {
      ORBMAT_OccupiedOnly = 0x01,
      ORBMAT_VirtualOnly = 0x02,
      ORBMAT_AlphaAndClosedOnly = 0x04,
      ORBMAT_BetaAndClosedOnly = 0x08,
      ORBMAT_AlphaOnly = 0x10,
      ORBMAT_ClosedOnly = 0x20,
      ORBMAT_BetaOnly = 0x40,
      ORBMAT_IaoBasis = 0x1000 // return pMinBasis as basis and make matrix of IAO basis coefficients, instead of full basis coefficients.
   };
   enum FWfInfoFlags {
      WFINFO_MakeRdms = 0x01
   };


   FGeometry *pGetGeometry();
   FGeometry const *pGetGeometry() const{ return const_cast<FFrame*>(this)->pGetGeometry(); };
   ct::FAtomSet *pGetAtoms();
   ct::FAtomSet const *pGetAtoms() const { return const_cast<FFrame*>(this)->pGetAtoms(); };
   FOrbital *pGetOrbital(int iMo); // note: this is an "orbital number" (i.e., starts a 1, goes to nMo).
   double GetEnergy() const;
   double GetGradient() const;
   // returns whether both frames have the same number and types of atoms
   bool IsGeometryCompatibleWith(FFrame const &other) const;

   QString GetBaseInputFileName() const;
   QString GetFullInputFileName() const;

   // collect all coefficients (and the current basis) of orbitals in *this.
   // notes:
   //   - COrb (output) is allocated on Mem.
   //   - If pRefOrbitals != 0, then a pointer to COrb(:,iOrb) will be stored in (*pRefOrbitals)[iOrb].
   void MakeOrbitalMatrix(ct::FMatrixView &COrb, FBasisSetPtr &pBasis2, TArray<FOrbital*> *pRefOrbitals, uint32_t Flags, ct::FMemoryStack &Mem);
   // splits all spin-free orbitals into separate alpha- and beta spin-orbitals
   void ConvertToUnrestrictedOrbitals();
//    void MakeIaoCoeffs(FWfOptions *pWfOptions, ct::FMemoryStack &Mem);
   void RunIaoAnalysis(ct::FLog &Log, FWfOptions *pWfOptions, bool AllowLocalize, ct::FMemoryStack &Mem);
   void RunIaoAnalysisImpl(ct::FLog &Log, FWfOptions *pWfOptions, bool AllowLocalize, ct::FMemoryStack &Mem);
   void MakeOrbitalMoments(ct::FMemoryStack &Mem);
   void LinkOrbitalsToPreviousFrame(FFrame *pPrevFrame, ct::FMemoryStack &Mem);
//    bool CanDoIaoAnalysis();

//    ct::FMatrixView MakeIaoBasis(ct::FAtomSet *pAtoms, ct::FMemoryStack &Mem);

   FFrameWfData *MakeWfDataForAnalysis(ct::FLog &Log, unsigned WfInfoFlags, unsigned OrbmatType, ct::FMemoryStack *pMem);
   void RunChargeAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms = 0);
   void RunBondOrderAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms = 0, ct::FMemoryStack *pMem = 0);
   void RunRedoxChargeAnalysis(ct::FLog &Log, FAtomIdList *pSelectedAtoms, FFragmentAnalysisOptionsCptr pFragmentAnalysisOptions, ct::FMemoryStack *pMem);

   void MakeFragmentationViaAtomGroups(FFragmentationDef &Fragmentation, ct::FAtomSet const *pAtoms, ct::FBasisSet const *pMinBasis, bool GroupRemainingAtoms);

//    void MakeFragmentationViaAtomGroups(FFragmentationDef &Fragmentation);

   // returns whether or not any orbital data sets are present.
   bool HaveOrbitals() const;
   // delete all electronic structure stuff we might be keeping.
   void ClearOrbitals();

   FWfType GetWfType();

   // find the data set row in which a given object is currently stored.
   int FindRow(FDataSet* pSet);
   int rowCount() const { return (int)m_Data.size(); }

   FDataSetList
      m_Data;
   FFreeObjectSetPtr
      // should this be here, rather than being in m_Data? Honestly, no idea.
      m_pFreeObjects;

   FDataSetList AllData();

   FMemoryLogQt &Log() { return *m_pFrameLog; }

   FArcLength GetArcLength() const { return m_ArcLength; }
   void SetArcLength(FArcLength const &ArcLength_) { m_ArcLength = ArcLength_; }

   FDocument *pGetDocument() { return m_pDocument; };
protected:
   FMemoryLogQt
      *m_pFrameLog;
   QString
      m_InputFileName;
//    FBasisSetPtr
//       m_pOrbBasis,
//       m_pMinBasis;
//    TArray<double>
//       m_CIaoData; // coefficients of IAOs in terms of main basis
//    ct::FMatrixView
//       m_CIb;
// ^- Note: actual data stored in orbitals now (use MakeOrbitalMatrix())
   FDocument
      *m_pDocument; // parent of the frame. used to take settings from.
   FArcLength
      m_ArcLength; // (approximate) position along the IRC.

//    // creates a new FBasisSet object representing a minimal basis for the atoms in
//    // this frame. By default this just makes a MINAO basis, but there may be
//    // adjustments if we are running in semi-empirical mode.
//    FBasisSetPtr MakeMinimalBasis();
};


// script interface
class IFrame : public FFrame
{
   Q_OBJECT

   Q_PROPERTY(QString name READ get_name WRITE set_name);
   Q_PROPERTY(int natoms READ get_num_atoms);
   Q_PROPERTY(int num_atoms READ get_num_atoms);
public slots:
   // add/remove bond lines. Atom indices are 1-based here.
   virtual void add_bond(int iAt, int jAt, QString const &Flags); // = 0;
   virtual void delete_bond(int iAt, int jAt); // = 0;
   virtual void reset_bonds(); // = 0; // reset all bonds to normal.
   virtual QString get_name();
   virtual void set_name(QString const &Name);

   virtual void scale_mo(int iMo, double fFactor); // = 0;
   virtual void rot_mos_2x2(int iMo, int jMo, double fAngle); // = 0;
   virtual int get_num_atoms() const;

   // collect basis function shells on the orbital basis of this frame, find unique ones,
   // and export them into library format (libmol).
   virtual void export_orbital_basis_set(QString const &FileName, QString const &BasisName);

//    virtual void add_objects(QList<FFreeObject*> const &fol, QString Mode = "clone");
   virtual void add_objects(QVariantList fol, QString Mode = "clone");
   virtual void clear_objects();
public:
   explicit IFrame(QString Desc_, FDocument *pDocument_);
};

typedef ct::TIntrusivePtr<IFrame>
   FFramePtr;


enum FSelectionMode {
   SELECT_Select, // delete previous selection and make a new selection with given objects
   SELECT_Toggle, // toggle selection state of given objects
   SELECT_Add // add to previous selection
};


struct FBondChangeAction : public QAction
{
   Q_OBJECT
public:
   enum FBondChangeActionType {
      ACTION_Hide,
      ACTION_SetStyleDotted,
      ACTION_Reset,
      ACTION_SetBondOrder
   };

   FBondChangeAction(int iAt_, int jAt_, FBondChangeActionType Type_, QString Text_, QObject *Parent_, double fArg = 0)
      : QAction(Text_, Parent_), m_Type(Type_), m_iAt(iAt_), m_jAt(jAt_), m_fArg(fArg)
   {}

   FBondChangeActionType
      m_Type;
   int
      m_iAt, m_jAt;
   double
      m_fArg;
};

typedef std::vector<int>
   FFrameIndexList;



QString ReplaceExt(QString const &FileName, QString const &NewExt);
QString RemovePath(QString const &FileName);
QString RemoveExt(QString const &FileName);


// used for alignment and IRC/other arc-length purposes.
struct FFrameCoords : public FIntrusivePtrDest
{
   ct::FAtomSet// const
      *pAtoms;
   std::vector<double>
      pAtMass;
   std::vector<ct::FVector3>
      pAtPos;
   FDocument
      *pDocument;
   explicit FFrameCoords(FGeometry *pGeometry);
   virtual ~FFrameCoords();

   bool empty() const { return pAtMass.empty(); }
   bool IsAlignableWith(FFrameCoords const &other) const { return pAtMass == other.pAtMass; }
};
typedef ct::TIntrusivePtr<FFrameCoords>
   FFrameCoordsPtr;

// idea: rows == associated data sets per frame,
//       cols == frames.
// (note that they may come in non-natural orders for later frames and that
//  some cells might be missing)
class FDocument : public QAbstractTableModel
{
   Q_OBJECT

public:
   Q_PROPERTY(QString atom_weight_mode READ GetAtomWeightMode WRITE SetAtomWeightMode NOTIFY AtomWeightModeChanged)
   Q_SLOT void SetAtomWeightMode(QString o);
//    Q_SLOT void SetAtomWeightMode(QString o) { m_AtomWeightMode = o; }
   Q_SIGNAL void AtomWeightModeChanged(QString o);
   QString GetAtomWeightMode() const { return m_AtomWeightMode; };

   Q_PROPERTY(bool skip_virtual_orbitals READ GetSkipVirtualOrbitals WRITE SetSkipVirtualOrbitals NOTIFY SkipVirtualOrbitalsChanged)
   Q_SLOT void SetSkipVirtualOrbitals(bool o) { m_SkipVirtualOrbitals = o; }
   Q_SIGNAL void SkipVirtualOrbitalsChanged(bool o);
   bool GetSkipVirtualOrbitals() const { return m_SkipVirtualOrbitals; };

   // hmpf... doc is not actually exposed to script currently. Script's 'doc' variable actually links to IApplication.
   // needs bigger re-design to work..
   Q_PROPERTY(FElementOptionsList elements READ GetElementOptions) // note: this COPIES the list!
   FElementOptionsList const &GetElementOptions() const { return m_ElementOptions; };
   FElementOptionsList &GetElementOptions() { return m_ElementOptions; };

   FFragmentAnalysisOptionsPtr GetFragmentAnalysisOptions() { return m_pFragmentAnalysisOptions; };

//    Q_PROPERTY(FWfOptions wf READ GetWfOptions) // note: this COPIES the list!
   FWfOptions *GetWfOptions() { return m_pWfOptions; }; // should this be here?
 public:
   FDocument(QObject *parent);
   int rowCount(const QModelIndex &parent = QModelIndex()) const;
   int columnCount(const QModelIndex &parent = QModelIndex()) const;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;

   uint &AtomFlags(int iAt);
   FAtomOptions &AtomOptions(int iAt);
   bool IsAtomHidden(int iAt);
   bool IsAtomSelected(int iAt);
   bool IsAtomExcludedFromAlignment(int iAt);
   FElementOptions *pElementOptions(int iAt, ct::FAtomSet *pAtoms = 0);
   // ^- pAtoms: if given, take data for atoms from this atom set. Otherwise: from current frame.

   void UnselectAll(bool EmitUpdate=true);
   void SelectAtom(int iAt, FSelectionMode SelectMode, bool EmitUpdate=true);
   void UpdateSelectedAtomsStatusText();
   void SelectAtomGroup(int iGroup);
   void DefineAtomGroup(int iGroup, FAtomIdList const &iAtoms);
   void DefineSelectionAsAtomGroup(int iGroup);
   // count the number of atoms which presently are selected.
   int nSelectedAtoms();
   int nAtomsWithFlags(); // count largest number of defined centers.

   void WriteExtendedStateScript(std::ostream &out);

   QString AtomLabel(int iAt) const;

   void MakeFragmentationViaAtomGroups(FFragmentationDef &Fragmentation, ct::FAtomSet const *pAtoms, ct::FBasisSet const *pMinBasis, bool GroupRemainingAtoms);
   size_t nAtomGroupsDefined() const;

   void AddVolumePropertyIsoSurface(FVolumePropertyInfo const &PropertyInfo, FIsoValueList const &IsoValues);

   QStringList GetLoadedDataFileList() { return m_LoadedDataFiles; };
public:
   typedef std::vector<FFramePtr>
      FFrameList;
   FDataSetPtr GetActiveDataSet();
//    ct::FAtomSet
//       Atoms;
   void Load(QStringList FileNames);

   // this will change the document by either re-ordering frames or deleting a
   // subset of the previous frames. iNewIndices gives the indices of the new
   // frames in the new order they are supposd to appear in. Any frame not
   // contained herein will be removed.
   void ReorderOrRestrictFrameSet(FFrameIndexList const &iNewIndices);

   IFrame *GetCurrentFrame();
   IFrame *GetCurrentFrame() const;
   FDataSetList *GetCurrentFrameData();
   FDataSetList *GetFrameData(int iFrame);
   QString GetLogFileName(QString LogType);

   IFrame *GetFrame(int Idx, bool AssertExists=true);
   FDataSet *GetRow(int Idx, bool AssertExists=true);
   FDataSet *GetRowCol(int iRow, int iCol, bool AssertExists=true);
   IFrame const *GetFrame(int Idx, bool AssertExists=true) const;
   FDataSet const *GetRow(int Idx, bool AssertExists=true) const;
   FDataSet const *GetRowCol(int iRow, int iCol, bool AssertExists=true) const;
   int GetActiveColIndex() const { return m_ActiveCol; };
   int GetActiveRowIndex() const { return m_ActiveRow; };
   int GetNumFrames() const { return int(m_Frames.size()); }

   void GetNextOrbitalColors(uint32_t &cIsoPlus, uint32_t &cIsoMinus);

   QString GetCurrentInputFileName();
   QString GetCurrentInputBaseFileName();
   QString GetCommonInputFileName();
// public slots:
//    void onToggleDataset(const QModelIndex &index);
//    void onSelectDataset(const QModelIndex &index);

   void AlignFrames(QString Mode);
   void SetInputFileName(QString FileName_) { m_InputFileName = FileName_; }
   QString GetAtomLabel(int iAt);

   void ClearOrbitals(bool EmitSignals=true);

   FDocumentMeasures *GetMeasures() { return m_pMeasures; };

   size_t iFrameId(FFrame *pFrame) const;
   bool HaveEnergies() const;
   bool HaveGradients() const;
   bool HaveOrbitals() const;

   // returns if either there is only a single frame, or there are multiple
   // frames, but all of them share the same number and types of atoms
   bool HaveConsistentFrames() const;

//    void SetArcLengthOptions(FArcLengthOptions ArcLengthOptions);
   FArcLengthOptions GetArcLengthOptions() const;
public slots:
   void ToggleDataRow(int iIndex);
   void SetActiveCol(int iIndex);
   // ForceUpdate: force emit signal even if new row is equal to old row
   // (e.g., if properties of the row have changed, like its activity status)
   void SetActiveRow(int iIndex, bool ForceUpdate = false, bool UpdateStatus=true);
   // move active column---up (+) or down (-) by amount iDelta. Function takes care
   // of limiting movements to existing frames.
   void MoveActiveCol(int iDelta);

   void HideSelectedAtoms();
   void ClearAtomFlags(bool EmitUpdate=true);
   void FindOrbitalsOnSelectedAtoms();
   void FindReactingOrbitals();
   void CalcChargesOnSelectedAtoms();
   void MakeHybridsForSelectedAtomGroup();
   void CopySelectedAtomNumbers();
   void ToggleActiveDataRow();
   void MakeBondLinesForSelectedAtoms();
   void ResetBondLines();
   void RunRedoxChargeAnalysis(FLogQt &Log);

   // link most-similar orbitals between all frames
   void LinkOrbitals();

   // data comes in over the QAction (sender), which should be a FBondChangeAction
   void ChangeSelectedBond();

   // these control the way new orbital colors are assigned.
   void SetNextOrbitalColorIndex(int Value);
   void SetNextOrbitalColorScheme(int Index);
   int GetNextOrbitalColorScheme();

   void RebuildWf(FLogQt &Log);
   void Clear(); // delete everything.

   void AddAxes(QString Which, double AxisLength, QString Options);

   // make a list of the smart pointers (in order to keep references to them)
   // to defer their point of deletion to a later point. Required due to GL-weiredness.
   // See rebuid wf code in IvMain.cpp.
   void AcquireFrameObjectLock(FDataSetList &ObjectLock);

   void UpdateArcLengths();
signals:
   void ActiveDatasetChanged();
   void ActiveColChanged(int iIndex);
   void ActiveRowChanged(int iIndex);
   void SelectionChanged();

   // i.e., "please update the 3d view"
   void VisualRepresentationChanged();
   void NextOrbitalColorIndexChanged(int Value);
protected:
   FElementOptionsList
      m_ElementOptions;

   FFrameList
      m_Frames;
   int
      m_ActiveRow,
      m_ActiveCol;
   QString
      m_InputFileName; // only set if this is a script?
   QStringList
      // set & appended by Load().
      m_LoadedDataFiles;
   int64_t
      m_SelectionSequenceId;

   typedef std::map<int, FAtomOptions>
      FAtomOptionMap;
   FAtomOptionMap
      m_AtomOptions;
   QString
      // default mode for aligning frame geometries to each other
      // and/or weighting atoms for arc lengths.
      m_AtomWeightMode;
   bool
      // if set, do not load or compute virtual orbitals
      m_SkipVirtualOrbitals;

   FFragmentAnalysisOptionsPtr
      m_pFragmentAnalysisOptions;

   FDocumentMeasures
      *m_pMeasures;

   // these two control the way new orbital colors are chosen.
   // (used when orbitals are first rendered iff they have no color sets
   // assigned otherwise)
   ptrdiff_t
      m_iNextOrbitalColor;
   size_t
      m_iNextOrbitalColorScheme;

   typedef std::map<int, FAtomIdList>
      FAtomGroupMap;
   FAtomGroupMap
      m_AtomGroups;

//    FArcLengthOptions
//       // controls how arc length coordinates (IRC etc.) are generated & exported.
//       m_ArcLengthOptions;

//    void FindAligningTrafo(double pR[9], double pD[3], ct::FAtomSet *pAtoms, std::string const &Mode, ct::FMemoryStack &Mem);
   FFrameCoordsPtr MakeFrameCoords(int iFrame);
   void FindAligningTrafo(double pR[9], double pD[3], FFrameCoords *pThis, FFrameCoords const *pLast, ct::FMemoryStack &Mem);

   // delete all cached geometric data of orbitals
   void ClearOrbitalRepresentations();

   void LoadFile(FFrameList &LoadedFrames, QString FileName);
   // returns true if the given file was identified as an orbital file (independent of whether loading it worked or not).
   bool LoadOrbitalFile(FFrameList &LoadedFrames, QString FileName, QString OrigFileName = "");
   void LoadXyzFile(FFrameList &LoadedFrames, QString FileName);

   void MakeOrbitalMoments();
   void MakeOrbitalCharges();
   void MakeOrbitalCachedData();

   // not: may be called only during document load.
   void BeginInsertFrames();
   void InsertFrame(FFramePtr pFrame);
   void EndInsertFrames();

   void BeginTotalReset();
   void EndTotalReset();

   int
      m_CallEndInsertRows; // -1: BeginInsertFrames not called.
   FFrameList
      m_FramesToInsert; // used only during InsertFrame operatations.
   FWfOptions
      *m_pWfOptions;

//    typedef std::vector<int>
//       FAtomIdList;
   FAtomIdList GetSelectedAtoms(bool SortedBySelectionSequence = true);
};

// just returns pDocument->pElementOptions(iAt, pAtoms); There as proxy to avoid #include of IvDocument.h for this only.
FElementOptions *pElementOptions(FDocument *pDocument, int iAt, ct::FAtomSet *pAtoms = 0);


// ^- note: automoc requires Q_OBJECT derived classes to be defined in
//    HEADER FILES. Otherwise it won't find them. Will result in "undefined
//    reference to 'vtable for ...'" errors.




#endif // IV_DOCUMENT_H
