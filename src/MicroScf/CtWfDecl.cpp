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

#include <stdexcept>
#include "CxParse1.h"
#include "CtWfDecl.h"
#include "CxIo.h" // for FInputError

namespace ct {

FWfDecl::FWfDecl(int iCharge_, unsigned Ms2_)
   : iCharge(iCharge_), Ms2(Ms2_), m_nNuclearCharge(-1)
{
}


unsigned FWfDecl::nElec() const {
   if (m_nNuclearCharge == -1)
      throw std::runtime_error("FWfDecl: cannot use nElec() before setting nuclear charge!");
   return unsigned(int(m_nNuclearCharge) - iCharge);
}


void FWfDecl::AssertNotInstanciated()
{
   if (m_nNuclearCharge != -1)
      throw std::runtime_error("FWfDecl: this object has already been assigned a charge. Further argument modifications may break things! Target usage: (a) make generic FWfDecl object (without assigned atom set), (b) in concrete applications, copy it and assign atom set to the copy, leaving the original intact.");
}

void FWfDecl::SetCharge(ptrdiff_t iCharge_)
{
   AssertNotInstanciated();
   this->iCharge = int(iCharge_);
   assert_rt(ptrdiff_t(this->iCharge) == iCharge_);
}


void FWfDecl::SetMs2(ptrdiff_t iMs2_)
{
   AssertNotInstanciated();
   if (iMs2_ < 0 || ptrdiff_t(unsigned(iMs2_)) != iMs2_)
      throw std::runtime_error(fmt::format("FWfDecl: Ms2/Spin argument must be a non-negative integer, but got '{}'.", iMs2_));
   this->Ms2 = unsigned(iMs2_);
}

void FWfDecl::SetArgs(std::string const &CommandOptions)
{
   AssertNotInstanciated();
   try {
      if (!CommandOptions.empty()) {
         FPropertyListStr
            OptgProps = FPropertyListStr(string_slice(CommandOptions), 0);
         for (size_t iOpt = 0; iOpt != OptgProps.size(); ++ iOpt) {
            FPropertyStr const
               &Prop = OptgProps[iOpt];
            if (Prop.Name == "charge") {
               SetCharge(Prop.Content.to_int());
            } else if (Prop.Name == "spin" || Prop.Name == "ms2") {
               ptrdiff_t ims2 = Prop.Content.to_int();
               SetMs2(ims2);
            } else if (Prop.Name == "read") {
               FileName_WfRead = Prop.Content.to_str();
            } else if (Prop.Name == "write") {
               FileName_WfWrite = Prop.Content.to_str();
            } else if (Prop.Name == "dump") {
               FileName_WfRead = Prop.Content.to_str();
               FileName_WfWrite = FileName_WfRead;
            } else {
               throw std::runtime_error(fmt::format("FWfDecl: option '{}: {}' not recognized.", Prop.Name.to_str(), Prop.Content.to_str()));
            }
         }
      }
   } catch (std::runtime_error &e) {
      // convert exceptions to input errors.
      throw FInputError(e.what());
   }
}

bool FWfDecl::HasEqualChargeSpinSym(FWfDecl const &other) const
{
   return iCharge == other.iCharge && Ms2 == other.Ms2 && nElec() == other.nElec();
   // ^-- note: having both iCharge and the nElec checks is not redundant. They
   // may differ if one wf uses ECPs and the other does not. We generally want
   // to treat such cases as different molecules, even if they describe the same
   // physical system (basically, for electronic structure methods purposes, we
   // mostly need to treat "Fe" and "Fe with ECP10MDF" as two different
   // elements)
}


} // namespace ct
