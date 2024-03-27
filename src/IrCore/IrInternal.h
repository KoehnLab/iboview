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

#ifndef IR_INTERNAL_H
#define IR_INTERNAL_H

#ifdef IR_TIMING
   #include "RDTSC.h"
   #define IR_RESUME_CLOCK(x) RESUME_CLOCK(x)
   #define IR_PAUSE_CLOCK(x) PAUSE_CLOCK(x)
#else
   #define IR_RESUME_CLOCK(x) (void)(0)
   #define IR_PAUSE_CLOCK(x) (void)(0)
#endif

namespace ir {
   typedef int
      timer_id_t;
#define FTimerId static const timer_id_t
   // timers for 2e2c integral driver routines
   FTimerId TID_MemsetTargetMatrix = 299; // a single full-size memset on the target matrix.
   FTimerId TID_EvalInt2e2c = 300;
   FTimerId TID_EvalCoShY = 301;
   FTimerId TID_EvalCoShY_Setup = 302;
   FTimerId TID_EvalCoShY_MDRR = 303;
   FTimerId TID_EvalCoKernels = 310;
   FTimerId TID_EvalCoKernels_PrimData = 311;
   FTimerId TID_EvalCoKernels_TildeGm = 312;
   FTimerId TID_EvalCoKernels_Contract = 313;
   FTimerId TID_ShTrA_YY = 315;
   FTimerId TID_Scatter2e2c = 316;

   // timers for 3c and 4c integral driver routines
   FTimerId TID_TimingEnvelope = 4; // outer loop enclosing all the driver calls.
   FTimerId TID_EvalInt2eNc = 200;
   FTimerId TID_InitialSetup = 204;
   FTimerId TID_PrepareShellPairs = 205;
   FTimerId TID_PrimLoop = 210; // everything inside the primitive loop combined
   // calc of rho/T/R/Sab/Scd etc (as done inside the main primitive loop, does not cover code under TID_PrepareShellPairs).
   FTimerId TID_PrimData = 212;
   FTimerId TID_EvalGm = 213; // evaluation of scalar integral kernel derivatives (EvalGm/EvalTildeGm)
   FTimerId TID_OsrrA = 214;
   FTimerId TID_OsrrB = 215;
   FTimerId TID_Ftcs_ShellPairExp = 220; // per-shell (ab| or |cd) FTCS data (CartXia or HermWiab; note: technically could be made free by precalc)
   FTimerId TID_Ftcs_Wabcd = 221; // real expansion vectors Wa0c0 or Wabcd (before PVTil contraction)
   FTimerId TID_Ftcs_EvalPrim = 222; // everything under (ab|cd) construction, given the Ws are known.
   FTimerId TID_Ftcs_TildePVg = 223; // assemble & Fourier-transform PV_kappa and g
   FTimerId TID_Ftcs_PVxPrime = 224; // PVxPrime & co
   FTimerId TID_Ftcs_Ux = 225; // construction final \Tilde U[k,a,b,c,d,kappa] matrices going into the FTCS formula
   FTimerId TID_Ftcs_Uyz = 226; // combining intermediates Uy and Uz into Uyz.
   FTimerId TID_Ftcs_Assembly = 229; // evaluation of the target primitive integral via the FTCS formula
   FTimerId TID_Dccs_ImUzUy = 230; // deferred evaluation of convoluted intermediates in standard order.
   FTimerId TID_Dccs_Convolution = 231; // evaluation of modified direct convolutions itself in the DCCS variant of the method
   FTimerId TID_Dccs_InnerAssembly = 232; // the dot-product part of the inner DCCS assembly (everything after the direct convolutions---*theoretically* this should not be more expensive than Rys assembly...)
//    FTimerId TID_Ftcs_ = 222; // everything under (ab|cd) construction, given the Ws are known.
   FTimerId TID_Finalize = 240;
   FTimerId TID_ShTr1 = 241;
   FTimerId TID_OsrrC1 = 242;
   FTimerId TID_ShTr2 = 243;
   FTimerId TID_OsrrC2 = 244;
   FTimerId TID_Contract = 250;
   // memsets on accumulation buffers for primitives to be contracted (could be considered a part of contract, too)
   FTimerId TID_ContractClear = 251;
//    FTimerId TID_ContractClear = TID_Contract;
   // accumulation operations in contraction (i.e., the actual contract)
//    FTimerId TID_ContractAcc = TID_Contract;
   FTimerId TID_ContractAcc = 252;
#undef FTimerId
}
// ^- todo: add these guys into the actual target routines instead of the driver routines?
//    would make code more readable, but would also mean that we have to recompile IrAmrr.cpp
//    everytime any timing stuff changes.

namespace ir {
   struct FTimerLock {
#ifdef IR_TIMING
      inline explicit FTimerLock(timer_id_t id) : m_id(id) { IR_RESUME_CLOCK(m_id); };
      inline ~FTimerLock() { IR_PAUSE_CLOCK(m_id); }
      timer_id_t m_id;
#else
      inline explicit FTimerLock(timer_id_t id) { IR_SUPPRESS_UNUSED_WARNING(id); };
      inline ~FTimerLock() {}
#endif // IR_TIMING
   private:
      FTimerLock(const FTimerLock&); // not implemented
      void operator = (const FTimerLock&); // not implemented
   };

// make a IR_TIME_SECTION(id) macro which will attribute everything between its
// instanciation and leaving the enclosing scope to a given timer id.
// (i.e., this one stops the timer automatically if the scope if left)
#ifdef IR_TIMING
#define IR_TIME_SECTION(id) ir::FTimerLock IR_UNIQUE_NAME(id)
#else
#define IR_TIME_SECTION(id) (void)(0)
#endif // IR_TIMING
}


#ifdef IR_TIMING
// FIXME: remove this
#include <sys/time.h>
namespace ir {
   inline double GetTicks1()
   {
      timeval
         tv;
      gettimeofday( &tv, 0 );
      return tv.tv_sec + 1e-6 * tv.tv_usec;
   }
}
#endif // IR_TIMING


#endif // IR_INTERNAL_H
