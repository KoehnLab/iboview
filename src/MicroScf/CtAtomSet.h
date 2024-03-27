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

#ifndef CT8K_ATOMSET_H
#define CT8K_ATOMSET_H

#include <map>
#include <string>
#include <vector>
#include <istream>
#include <stdexcept>

#include "CxTypes.h"
#include "CxVec3.h"
#include "CxMemoryStack.h"
#include "CtBasisSet.h"
#include "CtBasisDesc.h"
#include "CxAtomData.h"
#include "CxPodArray.h"
#include "CxRawAtom.h" // for conversion functions and FAtomTag
#include "CxXyzFrameIo.h"


namespace ct {


typedef TVector3<double>
   FVector3;

struct FDoublettIntegralFactory;
struct FMatrixView;


struct FAtom
{
   FVector3
      vPos; // center of atom's nuclear point charge, in a_bohr (always in a_bohr)
   int
      iElement;
   FAtomTag
      // optional explicit atom sub-type label defined in the input (e.g., to
      // distinguish two C atoms, tagging them as C1 and C2 in the input, to
      // give them different basis sets)
      Tag;
   unsigned
      // number of electrons implicit due to effective
      // core potential ("pseudopotential")
      nEcpElec;
   FBasisDescs
      BasisDesc;
   FEcpDesc
      EcpDesc;

   FVector3
      // if set: gradient dE/d[x,y,z] of the atom in a.u. (TODO: should this be here?
      // it is an /input. in certain kinds of xyz files, so in this sense it would be
      // an actual property...)
      vGrad;

   std::string GetElementName() const { return std::string(this->pElementName()); }
   std::string ElementName() const { return std::string(this->pElementName()); }
   char const *pElementName() const;
   // ^-- hm... I somehow ended up with three of them, all doing the same operation.

   // combines element with the optinal atom-specific numeric tag in this->Tag.
   std::string GetElementAndTag() const;
   // makes labels like '1 Se' (output atom labels start at 1, i.e., iAtom is increased by 1 for making the label)
   std::string GetAtomLabel(int iAtom, int iAtomWidth=3) const;
   std::string GetBasisDesc() const;

   FAtom() { Init0(); };
   // AutoFix: change/derive some standard parameters based on the input basis name; e.g., choose
   // ECPs or change def2-TZVPP to dhf-TZVPP for suitable elements.
   FAtom(FVector3 const &vPos_, std::string const &ElementName_, std::string
            const &BasisDesc_, std::string const &Ecp_ = "", bool AutoFix = true);
   FAtom(FVector3 const &vPos_, std::string const &ElementName_, FBasisDescs
            const &BasisDesc_, std::string const &Ecp_ = "", bool AutoFix = true);
   FAtom(FVector3 const &vPos_, FIoElementAndTag const &ElementAndTag_, FBasisDescs
            const &BasisDesc_, std::string const &Ecp_ = "", bool AutoFix = true);
   FAtom(FVector3 const &vPos_, FIoElementAndTag const &ElementAndTag_)
      { Init0(); vPos = vPos_; iElement = ElementAndTag_.iElement; Tag = ElementAndTag_.Tag; }


   // assigns a particular basis to the atom, via name, and if AutoFix=true, adjusts relevant
   // variant basis sets (e.g., dhf vs def2 or xx vs xx-PP).
   void AssignBasis(FBasisContext Context, std::string const &BasisDecl, bool AutoFix = true);

   // return the expected number of core electrons for an atom of the given element.
   // If ECP is assigned, this number is reduced by the number of ECP electrons.
   size_t nCoreElec() const;

   double NuclearCharge() const; // might differ from iElement if ECPs are used.
private:
   void Init0(); // zeros out all data members to semi-valid arguments.
   void Init(FVector3 const &vPos_, std::string const &ElementName_, std::string const &EcpName_, bool AutoFix);
   void AutoFixBasisAssociations();
};


// return the number of electrons we'd typically consider as core-level
// for a given element.
unsigned GetElementNumCoreElec(unsigned iElement);



// describes global effects acting on the current calculation. For example,
// external fields. If added at some point, should probably be done in a way
// compatible with the property evaluation interface (which actually is there
// already now)
struct FEnvironment
{

};

struct FBasisSet;

// represents a group of symmetry-equivalent atoms (i.e., "group" as a in
// 'collection of things', not in the mathematical sense)
struct FAtomSymGroup {
   // number of equivalnt atoms
   size_t nEqiv;
   // index into Atoms[]
   size_t iEqiv[8];
   // symmetry op moving iEqiv[0] to iEqiv[i]
   unsigned SymOp[8];
   // bit pattern of XYZ symmetry op not affecting the atom positions,
   // since the relevant coordinates are zero.
   unsigned Stabilizer;
};
FVector3 ApplySymOp(FVector3 const &In, unsigned Op);


// represents an evaluated numerical property (e.g., the expectation value
// of some operator) associated with the frame.
struct FAtomSetProperty : public FIntrusivePtrDest
{
   std::string
      m_Desc;
   unsigned
      // number of centers for which the property has been evaluated
      // (if not center-depedent, set to 1)
      m_nCenters,
      // number of electronic state combinations (in case of transition moments)
      // which have been stored.
      m_nStates,
      // number of components (e.g., dipole moment x,y,z -> 3) per center entry
      m_nComps;
   typedef std::vector<double>
      FDataArray;
   FDataArray
      // stored as nComps x m_nStates x nCenters dense array (comps in fast dimension)
      m_Data;
   FAtomSetProperty(std::string Desc_, unsigned nCenters_, unsigned nStates_, unsigned nComps_, FDataArray const &Data_);
   double Get(unsigned iComp, unsigned iState, unsigned iCenter) const;
   void SanityCheck() const;
};
typedef TIntrusivePtr<FAtomSetProperty>
   FAtomSetPropertyPtr;


// An atom set specifies which atoms are where and which properties these
// atoms do have in calculations (i.e. which basis/ecp to use and so on).
// Additionally,
// So an instance of this class basically describes the problem which
// we currently handle.
struct FAtomSet : public FIntrusivePtrDest
{
public:
   typedef std::vector<FAtom>
      FAtomList;

   FAtomSet();
   ~FAtomSet();

protected:
   FAtomList
      m_Atoms;
   FEnvironment
      m_Environment;
public:
   // --------------------------------------------------------------------------
   // Container interface for treating FAtomSet as a sequence of atoms.
   // --------------------------------------------------------------------------
   size_t size() const { return m_Atoms.size(); }
   FAtom const &operator [] (size_t iAt) const { return m_Atoms[iAt]; }
   FAtom &operator [] (size_t iAt) { return m_Atoms[iAt]; }
   void clear() { m_Atoms.clear(); }
   bool empty() const { return m_Atoms.empty(); }
   void reserve(size_t n) { m_Atoms.reserve(n); }
   // ^-- note: no "push_back" but AddAtom (see below), because AddAtom may
   // perform additional actions besides just adding an atom to the controlled
   // set (e.g., consistency checks, internal structure adjustments)

   typedef FAtomList::iterator iterator;
   typedef FAtomList::const_iterator const_iterator;
   iterator begin() { return m_Atoms.begin(); }
   iterator end() { return m_Atoms.end(); }
   const_iterator begin() const { return m_Atoms.begin(); }
   const_iterator end() const { return m_Atoms.end(); }

public:
   // --------------------------------------------------------------------------
   // Interfaces related to probing and modifying the molecular geometry
   // --------------------------------------------------------------------------

   // adds a single atom to the controlled set.
   void AddAtom(FAtom const &Atom);

   // returns whether this is the same "molecule for electronic structure purposes"
   // as `Other` (i.e., both have same elements, in same order, with same number of
   // ECP electrons per for each atom).
   // TODO: check dummy tags if introduced? Dummy atoms should be ignored.
   bool IsSameMolecule(FAtomSet const &Other) const;
   // returns whether *this and `Other` have the same elements, in the same order,
   // at the same physical location. If set, ThrPos > 0 may be used to give some
   // leeway for the position comparison.
   // Note: does *not* check for ECPs, basis sets, etc.
   bool IsSameGeometry(FAtomSet const &Other, double const ThrPos=0) const;

   // rotate and/or translate the basis (i.e., the positions of the atoms) in real space. R is a
   // 3x3 matrix and 'd' is a 3-vector such that r' = R (r - d) gives the transformed atomic positions.
   void Transform_3x4(double const *R, double const *d, FMemoryStack &Mem);


public:
   // --------------------------------------------------------------------------
   // Interfaces related to reading/writing .xyz files for the molecular geometry,
   // and other geometry load/store operations
   // --------------------------------------------------------------------------

   typedef xyz_io::FXyzLoadOptions FXyzLoadOptions;
   typedef xyz_io::FXyzPrintOptions FXyzPrintOptions;
   typedef xyz_io::FXyzLoadException FXyzLoadException;

   // Adds the molecules out of the Rasmol XYZ file to the current atom set.
   // Atoms will be added with basis string StandardBasis unless some other
   // basis is specified in *pOtherOptions.
   void AddAtomsFromXyzFile(std::string const &FileName, FBasisDescs const &DefaultBases, FXyzLoadOptions const *pOtherOptions = 0);
   void AddAtomsFromXyzStream(std::istream &str, std::string const &FileName, FBasisDescs const &DefaultBases, FXyzLoadOptions const *pOtherOptions);

   void PrintAsXyz(std::ostream &out, FXyzPrintOptions const &Options = FXyzPrintOptions()) const;

   // make a FIoAtomSet structure for the current geometry; this is a minimalistic
   // representation of the molecular geometry for input/output purposes (-> CxXyzFrameIo.h).
   xyz_io::FIoAtomSetPtr MakeIoAtomSet() const;
   // convert FIoAtomSet geometry specification into FAtom objects, and assign them
   // to *this. After adding the atoms, AssignBases(DefaultBases, ...) is applied
   // to the newly added atoms.
   void AddAtomsFromIoAtomSet(xyz_io::FIoAtomSet const &IoAtoms, FBasisDescs const &DefaultBases);

   // assigns a particular basis to all atom, via name, and if AutoFix=true, adjusts relevant
   // variant basis sets (e.g., dhf vs def2 or xx vs xx-PP).
   void AssignBasis(FBasisContext Context, std::string const &BasisName, bool AutoFix = true, ptrdiff_t iFirstAt = 0, ptrdiff_t iLastAt = -1);
   void AssignBases(FBasisDescs const &DefaultBases, bool AutoFix = true, ptrdiff_t iFirstAt = 0, ptrdiff_t iLastAt = -1);

public:
   // interfaces related to point-group symmetry

   // Get up to 8 symmetry operators of the subgroups of D2h.
   // Each op: bit 0 -> x, bit 1->y, bit 2->z; E.g., Op 0b011 means invariance
   // to *simultaneous* mirroring around x and y.
   // Gen[i] indexes into Ops. It designates which of the operators form the minimal
   // generator set on which the symmetry adapted basis functions are supposed to
   // be based. See also `FBasisSet::MakeSymmetryTransformation`.
   void FindMirrorSymmetryOps(unsigned Ops[8], size_t &nOps, unsigned Gen[3], size_t &nGen) const;
   void FindEquivalentAtoms(FAtomSymGroup *&pGroups, size_t &nGroups, FMemoryStack &Mem) const;
public:
   // --------------------------------------------------------------------------
   // Interfaces related to nucler interactions and interactions between nuclei
   // and electrons
   // --------------------------------------------------------------------------

   enum FNuclearPotentialFlags {
      POTENTIAL_IncludeNuclearElectrostatic = 0x0001,
      POTENTIAL_IncludeEcp = 0x0002,
      POTENTIAL_Default = POTENTIAL_IncludeNuclearElectrostatic | POTENTIAL_IncludeEcp
   };
   // ecp-reduced total nuclear charge (in a neutral molecule, that's also the
   // number of explicit electrons which have to be there; in wmme's FAtomSet
   // this is called `nElecNeutral`).
   int    NuclearCharge() const;

   // returns ecp-reduced number of electrons in rare gas configurations.
   size_t nCoreElec() const;
   // returns total number of electrons absorbed in ECPs
   size_t nEcpElec() const;


   // ecp-aware nuclear repulsion energy.
   double NuclearRepulsionEnergy() const;
   // ecp-aware nuclear contribution to the molecualar dipole moment around
   // ExpansionPoint.
   FVector3 NuclearDipoleMoment(FVector3 const &ExpansionPoint = FVector3(0,0,0)) const;


   // makes a matrix representing a given one electron operator, and adds it
   // to Dest (which must either be 0x0 or have compatible dimensions).
   // The operator here can be specified generally using a fitting integral
   // driver object for it. For some operators of interest, direct evaluation
   // functions are provided (see next functions).
   // FDoublettIntegralFactory and derived classes are described in integrals.h
   void Add1eIntegralMatrix( FMatrixView &Dest,
      FBasisSet const &RowBasis, FBasisSet const &ColBasis,
      FDoublettIntegralFactory &IntFactory, double Factor, FMemoryStack &Mem ) const;
   void Make1eIntegralMatrixMultiComp( double *pDest,
      FBasisSet const &RowBasis, FBasisSet const &ColBasis,
      FDoublettIntegralFactory &IntFactory, double Factor, FMemoryStack &Mem ) const;
   // the following functions are more or less dummy functions that just call
   // Make1eIntegralMatrix.
   void MakeCoreHamiltonMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const;
   void MakeOverlapMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const;
   void AddKineticMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem) const;
   void AddNuclearAttractionMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem, unsigned PotentialFlags=POTENTIAL_Default) const;
   void AddDipoleMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis,
      FVector3 const &Direction, FVector3 const &ExpansionPoint, FMemoryStack &Mem) const;
   // for quadrupoles etc. CartMoment = [i,j,k] gives the powers of <\mu| x^i y^j z^k |\nu>.
   void MakeCartesianMomentMatrix(FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis,
      TVector3<unsigned> const &CartMoment, FVector3 const &ExpansionPoint, FMemoryStack &Mem) const;

   // form 1st derivative integral contraction GradXyz[Ax] = \sum{\mu,\nu} Rdm[\mu,\nu] d/d[Ax] (\mu|CoreH|\nu)
   // where Ax denotes the combination of atomic center iNuc and xyz cartesian derivative.
   void AccCoreHamiltonDeriv1( FMatrixView &GradXyz, FMatrixView const &Rdm, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const;
   // add gradient of Nuclear-nuclear repulsion energy to GradXyz
   void AccNuclearRepulsionGradient(FMatrixView &GradXyz, FMemoryStack &Mem, double Factor=1.) const;
   // add Hessian of Nuclear-nuclear repulsion energy to HessXyz ((3N,3N)-matrix)
   void AccNuclearRepulsionHessian(FMatrixView &HessXyz, FMemoryStack &Mem, double Factor=1.) const;

   // idea of those things is to adjust them when more complex external fields
   // (e.g., point charges or ECPs) or other kinetic energy kernels (e.g.,
   // scalar-relativistic ones) are used. But maybe other things need to be
   // adjuted anyway.
   FKrn2iPtr MakeNuclearPotentialKernel(unsigned PotentialFlags=POTENTIAL_Default) const;
   FKrn2iPtr MakeKineticEnergyKernel() const;
   bool NeedEcpKernel() const; // returns true if there are ECPs defined on any atom.
public:
   // --------------------------------------------------------------------------
   // Interfaces related to summary properties and input/output
   // --------------------------------------------------------------------------

   // these are mainly here for simplifying dealing with xyz files describing
   // reaction paths, geometry optimizations and similar things.
   double GetLastEnergy() const { return m_LastEnergy; }
   void SetLastEnergy(double Energy) { m_LastEnergy = Energy; }
   bool HaveGradients() const;
   double GetTotalGradient() const;

   std::string const &GetCaption() const { return m_CaptionLine; }
   void SetCaption(std::string const &Caption_) { m_CaptionLine = Caption_; }
   std::string const &GetName() const { return m_InputName; };
   void SetName(std::string const &InputName_) { m_InputName = InputName_; };

   // makes labels like '1 Se' (output atom labels start at 1, i.e., iAtom is increased by 1 for making the label)
   // iAtomWidth: width of field for atom number. If -1, decide automatically (based on number of atoms in set).
   // If 0, use smallest fitting. If 1,2,3,4: use specified witdth.
   std::string GetAtomLabel(int iAtom, int iAtomWidth=-1) const;
protected:
   double
      m_LastEnergy; // last energy computed for the set
   std::string
      m_InputName; // if read from an .xyz file: remember name of the file here. should probably not be here...
   std::string
      m_CaptionLine; // if read from an .xyz file: the xyz caption line.
public:
   // --------------------------------------------------------------------------
   // Interface for evaluated physical properties (e.g., dipole moments, energies,
   // etc) attached to the atom set.
   // --------------------------------------------------------------------------
   //
   // TODO: needs cleanup. Probably shouldn't be all public, and come in two
   // different versions...
   typedef std::vector<FAtomSetPropertyPtr>
      FPropertyList;
   FPropertyList
      m_Properties;
   void ClearProperties() { m_Properties.clear(); }
   void AddProperty(FAtomSetPropertyPtr pProp) { m_Properties.push_back(pProp); }
   void SwapPropertyList(FPropertyList &NewProperties) { m_Properties.swap(NewProperties); }
   size_t nStoredPropertyEntries() const;
   struct FPropertyDataEntry {
      size_t iProperty, iCenter, iState, iComp;
      double fValue;
   };
   typedef std::vector<FPropertyDataEntry>
      FPropertyDataEntryList;
   void CollectPropertyData(FPropertyDataEntryList &PropertyData) const;
   FPropertyDataEntryList CollectPropertyData() const { FPropertyDataEntryList r; this->CollectPropertyData(r); return r; }
};
typedef TIntrusivePtr<FAtomSet>
   FAtomSetPtr;
typedef TIntrusivePtr<FAtomSet const>
   FAtomSetCptr;

typedef std::vector<FAtomSetPtr>
   FAtomSetList;

// load the xyz frame(s) from the given file and add them to Frames.
void LoadMultiXyz(FAtomSetList &Frames, std::string const &FileName, FBasisDescs const &DefaultBases, FAtomSet::FXyzLoadOptions const *pOtherOptions = 0);

std::ostream &operator << (std::ostream &out, FAtomSet const&);
std::ostream &operator << (std::ostream &out, FVector3 const &Pos);

// make a alignment transformation to be used with FAtomSet::Transform_3x4 and/or FBasisSet::Transform_3x4.
// That probably should not be here; at least not in this form (note: this used to be in IboView).
void FindAtomSetAligningTrafo(double pR[9], double pD[3], FVector3 const *pAtPos, double const *pAtMass, FVector3 const *pAtPosLast, size_t nAt, FMemoryStack &Mem);
void Transpose3x3(double *R);
void Trafo3x4(FVector3 &InOut, double const *R, double const *d);




} // namespace ct


#endif // CT8K_ATOMSET_H

// kate: indent-width 3; indent-mode normal;
