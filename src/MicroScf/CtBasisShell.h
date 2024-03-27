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

#ifndef CT_BASIS_SHELL_H
#define CT_BASIS_SHELL_H

#include "CxTypes.h"

#include <iosfwd>
#include <string>

#include "CxPodArray.h"
#include "CxVec3.h"
#include "Ir.h" // only for making ir:FRawShell instances.


namespace ct {

typedef TVector3<double>
   FVector3;

/// A structure containing auxiliary information about a concrete basis set,
/// particularly information about where the data originates from.
///
/// Warning: These objects may be shared across FAtomShell entries, including
/// FAtomShell entries of different elements!
struct FElementShellInfo : public FIntrusivePtrDest {
   std::string
      /// optional: actual instanciation command used to generate this object
      /// in case this is *not* just a standard library shell.
      Definition,
      /// optional: actual basis set name in the source record of a basis set
      /// library this object was instancated as (e.g., univ-ecp28-JFIT for
      /// a univ-JFIT loaded for Rb).
      ///
      /// UPDATE: no... this does not work. Almost all of the basis sets
      /// have partially identical data sets under many different names...
      /// which are presently shared. I'd have to unlink those objects for
      /// the name thing to work in here.
      /// I replaced this with the pLibraryName mechanism in FBasisShell
      /// objects. That is possibly a bit too fragile, but as long as all
      /// the shells come from FBasisLibrary objects, and the FBasisLibrary
      /// stays around (it kind of has to), that may be aceptable.
//       LibraryName,
      /// optional: comment line for this basis in the basis set library
      LibraryComment;

   FElementShellInfo() {};
   FElementShellInfo(std::string const &Definition_, std::string const &LibraryName_, std::string const &LibraryComment_);

   /// these ones compare actual definition data. They are meant to allow
   /// using FElementShellInfo objects as map/set keys to identify unique ones
   /// if needed.
   bool operator < (FElementShellInfo const &other) const;
   bool operator == (FElementShellInfo const &other) const;
   bool operator != (FElementShellInfo const &other) const;
};
typedef TIntrusivePtr<FElementShellInfo>
   FElementShellInfoPtr;
typedef TIntrusivePtr<FElementShellInfo const>
   FElementShellInfoCptr;



/// A generally contracted shell of basis functions of an atom:
///
/// - A shell contains all contracted radial functions which are defined by the
///   fixed linear combination in CoMatrix
///
/// - A shell contains all spherical components of the solid harmonic S(l,m).
///   E.g., for a d-shell, all five components d2- d1- d0 d1+ and d2+ are
///   included, for each radial function. (e.g., a d shell with two contractions
///   has ten functions).
///
/// - This object DOES NOT contain a basis function center. These objects are
///   intended to be managed by a basis set library and shared across FBasisShell
///   objects (see below) of the same atom type. For example, two C atoms sitting on
///   different centers will typically share the same FAtomShell object,
///   and this object will typically be controlled by a basis set library object.
///
/// - FIXME: rename to FElementShell (and possibly the one down (FBasisShell)
///   into FAtomShell?)
struct FAtomShell : public FIntrusivePtrDest {
   TArray<double>
      /// one exponent for each primitive.
      Exponents,
      /// nExp x nCo contraction matrix. Contractions are stored
      /// with respect to unnormalized Gaussians.
      CoMatrix,
      /// range of contracted and primitive functions in IR format, or 0.
      RangeInfo;
   unsigned
      /// Angular momentum l (0 -> s, 1 -> p, 2 -> d, ...)
      AngMom;
   FElementShellInfoCptr
      /// optional data entry to allow retaining data on where the object
      /// originated.
      pInfo;
   int
      /// optional data entry specifying for which element this basis set
      /// was instanciated for. Defaults to 0 (element not set).
      iLibraryElement;

   unsigned l() const { return AngMom; }
   size_t nExp() const { return Exponents.size(); }
   size_t nCo() const { return CoMatrix.size() / Exponents.size(); }
   size_t nSh() const { return 2 * l() + 1; }
   size_t nFn() const { return nCo() * nSh(); }
   double fCo(unsigned iExp, unsigned iCo) const { assert(iExp < nExp() && iCo < nCo()); return CoMatrix[iExp + nExp() * iCo]; }
   double fExp(unsigned iExp) const { assert(iExp < nExp()); return Exponents[iExp]; }

   /// returns pInfo->Definition if set, otherwise an empty string.
   std::string const &Definition() const;
//    /// returns pInfo->LibraryName if set, otherwise an empty string.
//    std::string const &LibraryName() const;
   /// returns pInfo->LibraryComment if set, otherwise an empty string.
   std::string const &LibraryComment() const;


   enum FInitType {
      // if set, input co matrix is considered as referring to raw Gaussians
      // instead of normalized Gaussians.
      TYPE_Unnormalized = 0x01
   };

   /// create an empty shell. Fill components yourself.
   FAtomShell() {};
   /// create a generally contracted shell. By defauly pCo is a nExp x nCo matrix
   /// given in library format (i.e., referring to normalized primitive Gaussians).
   FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_=0);
   /// create a single normalized primitive.
   FAtomShell(unsigned AngMom_, double const fExp_, unsigned InitFlags_=0);
   /// create a number of primitives.
   FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, unsigned InitFlags_=0);


   /// attach source information about the basis shell to the object.
   void AttachInfo(int iElement_, FElementShellInfoCptr pInfo_);
   /// Perform internal consistency checks and build auxiliary data structures
   /// required for using the object. By default this also computes the RangeInfo
   /// entry, unless already there.
   void Finalize();
   /// calculate screening ranges (max/min) in IR format. pMinRange may be 0. pMaxRange and
   /// pMinRange must hold room for (1 + nCo + nExp) scalars.
   /// ThrEl specifies the maximum fraction of an electron per basis function we are willing to lose.
   /// nRes denotes the resolution for radial integration.
   void FindMinMaxRange(double *pMaxRange, double *pMinRange, uint nRes, double ThrEl) const;

   /// creates a text representation of the atomic basis function shell data this
   /// object represents, in Molpro basis library format (libmol), and appends it to
   /// the end of output stream xout. Notes:
   /// - ElementName (e.g., 'Au') and BasisName (e.g., 'cc-pVTZ') denote which element
   ///   and basis export name the exported set will be filed under.
   /// - BasisName should not contain non-alphanumeric characters except for '-'. Case
   ///   is ignored in libmol format.
   /// - "Comment" is what goes into the comment/reference line entry
   /// - It is possible to specify multiple export names, by separating them with a
   ///   space (e.g., "cc-pVTZ VTZ" would make the basis recognized under both those
   ///   names)
   /// - This function can be used to export basis set data obtained from other sources
   ///   than libmol (e.g., from xml files or molden files) into libmol format for later reuse.
   void ExportToLibmol(std::ostream &xout, std::string const &ElementName, std::string const &BasisName, std::string const &Comment) const;

   /// these ones compare actual basis shell data. They are meant to allow
   /// use FAtomShell objects as map/set keys to identify unique ones if needed.
   bool operator < (FAtomShell const &other) const;
   bool operator == (FAtomShell const &other) const;
   bool operator != (FAtomShell const &other) const;
private:
   void Init(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_);
};

typedef TIntrusivePtr<FAtomShell>
   FAtomShellPtr;
typedef TIntrusivePtr<FAtomShell const>
   FAtomShellCptr;

// A shell of generally contracted molecular basis functions. Contrary to
// FAtomShell, this object is centered in space and represents an actual basis
// function used in the current calculation.
struct FBasisShell
{
   FAtomShellCptr
      pAs;
   FVector3
      vCenter;
   int
      // center index; not used by the low-level integral driver, but may be
      // used by the program or by mid-level functions to identify atoms.
      iCenter;

   FBasisShell() {};
   FBasisShell(FVector3 vCenter_, int iCenter_, FAtomShellCptr pAs_, std::string const *pLibraryName_ = 0)
      : pAs(pAs_), vCenter(vCenter_), iCenter(iCenter_), pLibraryName(pLibraryName_)
   {}

   ir::FRawShell MakeIrShell();

   unsigned l() const { return pAs->l(); }
   unsigned nExp() const { return pAs->nExp(); }
   unsigned nCo() const { return pAs->nCo(); }
   unsigned nSh() const { return pAs->nSh(); }
   unsigned nFn() const { return pAs->nFn(); }
   double fCo(unsigned iExp, unsigned iCo) const { return pAs->fCo(iExp, iCo); }
   double fExp(unsigned iExp) const { return pAs->fExp(iExp); }

   // returns actual basis set library name this object was instanciated as, if
   // known, and otherwise an empty string.
   std::string const &LibraryName() const;
   // returns basis set library comment line associated with this object if known,
   // otherwise an empty string.
   std::string const &LibraryComment() const;
   // returns chemical element this basis set was listed under in the basis set
   // library, if known, and otherwise zero.
   //
   // WARNING: the library element is not necessarily identical with the actual
   // chemical element of the atom the basis is placed on (one may place the
   // basis designed for one atom on another).
   int iLibraryElement() const;

   // prints parameters of shell. For the intend of making tables, for all lines after the first,
   // this routine first emits "Ind" and then spaces up to a total line length of 'TotalIndent',
   // before emitting the actual shell parameters.
   void PrintAligned(std::ostream &xout, std::string const &Ind, unsigned TotalIndent) const;

protected:
   std::string const
      // optional: if set, a pointer into the FBasisLibrary object which
      // indicates the actual basis set name under which this object was
      // instanciated.
      *pLibraryName;
};

// // calculate and return the integral Sqrt[<mu|mu>] for mu being a raw primitive Gaussian.
// // Used for conversion of contraction coefficients between raw Gaussians (integral
// // driver format) and normalized Gaussians (library format).
// //
// // This function accounts for the radial part. The angular normalization is
// // included in the coefficients of Slm(r).
// double RawGaussNorm(double fExp, unsigned l);
//
// ^- note: moved to IrCore.cpp

} // namespace ct

#endif // CT_BASIS_SHELL_H
