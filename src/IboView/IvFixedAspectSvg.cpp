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

#include <QSizePolicy>
#include "IvFixedAspectSvg.h"
// #include "Iv.h"

FFixedAspectSvg::FFixedAspectSvg(QWidget *parent_)
   : QSvgWidget(parent_)
{
   QSizePolicy p = sizePolicy();
   p.setHeightForWidth(true);
   setSizePolicy(p);
//    IvEmit("Entered FFixedAspectSvg::FFixedAspectSvg.");
//    IvEmit("this->sizePolicy().hasHeightForWidth(): %1", this->sizePolicy().hasHeightForWidth());
}

int FFixedAspectSvg::heightForWidth(int w) const
{
//    IvEmit("Entered FFixedAspectSvg::heightForWidth. w = %1", w);
   QSize svgSize = sizeHint();
//    QSize svgSize(930,360);
   return int(double(w)/double(svgSize.width())*double(svgSize.height()));
   // ^- wheee...works. That was more annoying than expected...
}
