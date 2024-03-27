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

#ifndef ISO_SURFACE_H
#define ISO_SURFACE_H

#include "CxTypes.h"
#include "CxPodArray.h"
#include "IvMesh.h"

enum FIsoSurfaceFlags {
   ISOFLAGS_PlusAndMinus = 0x01,
   ISOFLAGS_AllowFlip = 0x02 // if true, put the larger lobe on the positive side.
};

enum FIsoType {
   ISOTYPE_Undefined,
   ISOTYPE_Absolute,
   ISOTYPE_Relative
};

struct FIsoThresholdEntry {
   // either the absolute iso value c, or the value relative to the integrated data weight.
   // In this case, 0.80 would be 80%.
   // note: this is SIGNED. If you want -x and +x, you need to make two entries!
   double
      fIsoValue;
   uint32_t
      dwColor;
   FIsoThresholdEntry() {};
   FIsoThresholdEntry(double fIsoValue_, uint32_t dwColor_) : fIsoValue(fIsoValue_), dwColor(dwColor_) {}
};
typedef TArray<FIsoThresholdEntry>
   FIsoValueList;

struct FIsoSurfaceSettings : public ct::FIntrusivePtrDest
{
   // absolute or relative based iso values?
   FIsoType
      IsoType;

   FIsoValueList
      IsoValues;

   double
      // number of points per Angstrom in each direction. surface grid is derived from that.
      fLinearDensity,
      // value to use for integrated weight of property if relative absolute thresholds are used.
      // if 0: automatic (take from numerically integrated data weight)
      fRelativeIsoThresholdRefValue;

   // bitfield of ISOFLAGS_*
   unsigned
      Flags;

   bool PlusAndMinus() const { return 0 != (Flags & ISOFLAGS_PlusAndMinus); }
   bool AllowFlip() const { return PlusAndMinus() && (0 != (Flags & ISOFLAGS_AllowFlip)); }

//    FIsoSurfaceSettings(FIsoType IsoType_ = ISOTYPE_Relative, double fLinearDensity_ = 20., uint Flags_ = ISOFLAGS_PlusAndMinus);

   FIsoSurfaceSettings(FIsoType IsoType_ = ISOTYPE_Relative, double fIsoValue_ = 0.8,
      double fLinearDensity_ = 20., unsigned Flags_ = ISOFLAGS_PlusAndMinus,
      uint32_t dwColorPlus_ = 0xffff0000, uint32_t dwColorMinus_ = 0xff0000ff);
};
typedef ct::TIntrusivePtr<FIsoSurfaceSettings>
   FIsoSurfaceSettingsPtr;

struct FVolumeDataSet;

// FIndexedTriangleListPtr MakeIsoSurface(FOrbital const *pOrbital, FIsoSurfaceSettings const &Options);
FIndexedTriangleListPtr MakeIsoSurface(FVolumeDataSet const *pVolume, FIsoSurfaceSettings const &Options);


#endif // ISO_SURFACE_H
