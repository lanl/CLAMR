/*
 *  Copyright (c) 2011-2012, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2011-2012. Los Alamos National Security, LLC. This software was produced 
 *  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
 *  ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
 *  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 *  to produce derivative works, such modified software should be clearly marked,
 *  so as not to confuse it with the version available from LANL.
 *
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE LOS ALAMOS NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
 *  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *  
 *  CLAMR -- LA-CC-11-094
 *  This research code is being developed as part of the 
 *  2011 X Division Summer Workshop for the express purpose
 *  of a collaborative code for development of ideas in
 *  the implementation of AMR codes for Exascale platforms
 *  
 *  AMR implementation of the Wave code previously developed
 *  as a demonstration code for regular grids on Exascale platforms
 *  as part of the Supercomputing Challenge and Los Alamos 
 *  National Laboratory
 *  
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Neal Davis              davis68@lanl.gov, davis68@illinois.edu
 *           David Nicholaeff        dnic@lanl.gov, mtrxknight@aol.com
 *           Dennis Trujillo         dptrujillo@lanl.gov, dptru10@gmail.com
 * 
 */
#ifndef _ZORDER_H
#define _ZORDER_H

#if !defined(REG_INTEGER) && !defined(SHORT_INTEGER) && !defined(MIN_INTEGER)
#define REG_INTEGER
#endif

#if defined(MIN_INTEGER)
   // define all to needed ranges and then typedef or define to actual
   typedef unsigned short ushort_t; // 0 to 65,535
   typedef short          short_t;  // -32,768 to 32,767
   typedef unsigned char  uchar_t;  // 0 to 255
   typedef char           char_t;   // -128 to 127 
#elif defined(SHORT_INTEGER)
   typedef unsigned short ushort_t;
   typedef short          short_t;
   typedef unsigned short uchar_t;
   typedef short          char_t;
#elif defined(REG_INTEGER)
   typedef unsigned int ushort_t;
   typedef int          short_t;
   typedef unsigned int uchar_t;
   typedef int          char_t;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

void calc_zorder(int size, int *i, int *j, uchar_t *level, int levmx, int ibase, int *z_index, int *z_order);
unsigned long long index_to_bit(unsigned long long index, int lev, int levmx, int ibase);
unsigned long long twobit_to_index(unsigned long long ibit, unsigned long long jbit);
void printbits(int n);

#ifdef __cplusplus
}
#endif

#endif  /* _ZORDER_H */

