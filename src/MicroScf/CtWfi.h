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

/// @file CtWfi.h
///
/// MicroScf Wave Function Interface (wfi): Include support for
/// import/export/transfer of wave function data (in particular, orbital sets
/// and related data) for the sake of import/export/transfer (between programs),
/// conversion, wave function analysis (wfa) etc.
///
/// Well... that was the idea. In practice the code is still widely scattered
/// and not actually moved here yet. See comments in notes/notes-on-localization-interface.txt
///
/// TODO:
/// - merge FWfType FFrame::GetWfType()
/// - merge wfa's FFrameWfData
#ifndef CT_WFDATA_H
#define CT_WFDATA_H

namespace wfi {

enum FOrbitalSpin {
   ORBSPIN_SpinFree,
   ORBSPIN_Alpha,
   ORBSPIN_Beta,
   ORBSPIN_Unknown
};
char const *pOrbDescFromType(FOrbitalSpin Type, double fOcc);
// ^- maybe I should even use bitfields for Alpha and Beta? Then closed would be ALPHA | BETA.


enum FWfType {
   WFTYPE_Rhf, // only orbitals with occupancy closed or alpha
   WFTYPE_Uhf, // only orbitals with occupancy alpha or beta
   WFTYPE_Mcscf, // spin-free orbitals, but with occupation numbers other than 1.0 or 2.0 (can still make IAO charges, but not bond orders or IBOs for the active space)
   WFTYPE_Other
};
// ^-- TODO: the stuff in here needs some cleaning and orthogonalizing. See notes-on-localization-interface.txt
// (e.g., restricted/unrestricted/general and single-determinant-scf/mcscf/non-scf are really orthogonal concepts...)


} // namespace wfi



#endif // CT_WFDATA_H
