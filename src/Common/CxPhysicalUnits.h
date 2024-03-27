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

#ifndef CT_CONSTANTS_H
#define CT_CONSTANTS_H

namespace ct {
   static double const
      ToAng2006 = 0.529177249, // 2006 value.
//       ToAng = 0.5291772108,
//       ToAng = 0.529177209, // molpro default (up to development version 2014-08-27).
      ToAng = 0.52917721092, // molpro default (from 2014-08-27; CODATA 2010 value)
      ToAngCp2k = 0.529177211,
      ToEv2006 = 27.2113961,
      ToKcal = 627.5096,
      ToEv = 27.21138505,
      ToDebye = 2.54158;
} // namespace ct



#endif // CT_CONSTANTS_H
