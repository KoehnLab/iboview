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

#ifndef FORTRAN_INT_H
#define FORTRAN_INT_H

#include <stddef.h>
#include "CxTypes.h"
// #include <boost/cstdint.hpp>
//#include <stdint.h>
// ^- get those from C99 stdint.h (technically not C++ standard)
//    or from boost cstdint.hpp from outside.

// basic fortran integer type.
#ifdef _I4_
   typedef int32_t FORTINT;
#else
//    typedef boost::int64_t FORTINT;
   typedef ptrdiff_t FORTINT;
   // ^- (should default to 32bit for 32bit systems. Note that long does not have
   //     this property... it's 32bit on 64bit windows...)
#endif

typedef FORTINT const
   &FINTARG;
typedef double const
   &FDBLARG;

// import declarations for extern "C" functions without (FC_FUNC)
// and with (FC_FUNC_) underscores in their symbol name.
#define FC_FUNC(name,NAME) name ## _
#define FC_FUNC_(name,NAME) name ## _


// number of trailing underscores functions with FORTRAN signature get.
#ifndef FORT_UNDERSCORES
   #define FORT_UNDERSCORES 1 /* most linux fortran compilers do it like that. */
#endif

// macro for defining C functions which are callable from FORTRAN side
// (and reverse direction):
// A function
//     void FORT_Extern(bla,BLA)(FORTINT &a, double &b);
// can be called as a function
//     subroutine bla(a, b)
//       integer :: a
//       double precision :: b
// from Fortran side.
#ifdef FORT_UPPERCASE
  #if FORT_UNDERSCORES == 0
    #define FORT_Extern(lowercasef,UPPERCASEF) UPPERCASEF
  #elif FORT_UNDERSCORES == 1
    #define FORT_Extern(lowercasef,UPPERCASEF) UPPERCASEF##_
  #elif FORT_UNDERSCORES == 2
    #define FORT_Extern(lowercasef,UPPERCASEF) UPPERCASEF##__
  #else
    #error "should define FORT_UNDERSCORES for fortran function signatures."
  #endif
#else
  #if FORT_UNDERSCORES == 0
    #define FORT_Extern(lowercasef,UPPERCASEF) lowercasef
  #elif FORT_UNDERSCORES == 1
    #define FORT_Extern(lowercasef,UPPERCASEF) lowercasef##_
  #elif FORT_UNDERSCORES == 2
    #define FORT_Extern(lowercasef,UPPERCASEF) lowercasef##__
  #else
    #error "should define FORT_UNDERSCORES for fortran function signatures."
  #endif
#endif

#endif /* FORTRAN_INT_H */
