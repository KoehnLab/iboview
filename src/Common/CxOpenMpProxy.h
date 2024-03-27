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

// include OpenMP header if available or define inline dummy functions
// for OMP primitives if not.
#ifndef OPENMP_PROXY_H
#define OPENMP_PROXY_H

#ifdef _OPENMP
   #include <omp.h>
#else
   inline int omp_get_thread_num() { return 0; } // current thread id
   inline void omp_set_num_threads(int) {}
   inline int omp_get_max_threads() { return 1; } // total number of threads supposed to be running.
   inline int omp_get_num_procs() { return 1; } // total number of "virtual" CPU cores (includes Hyperthreading cores)

   struct omp_lock_t {};
   inline void omp_destroy_lock(omp_lock_t *){}
   inline void omp_init_lock(omp_lock_t *){}
   inline void omp_set_lock(omp_lock_t *){}
   inline void omp_set_nested(int){}
#endif

#endif // OPENMP_PROXY_H
