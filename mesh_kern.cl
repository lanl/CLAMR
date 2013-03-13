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

#define HASH_SETUP_OPT_LEVEL 4

#ifndef GPU_DOUBLE_SUPPORT
#define GPU_DOUBLE_SUPPORT
#ifdef HAVE_CL_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double  real;
typedef double2 real2;
typedef double4 real4;
#define ZERO 0.0
#define HALF 0.5
#define ONE  1.0
#define GRAVITATIONAL_CONSTANT 9.80
#else
typedef float   real;
typedef float2  real2;
typedef float4  real4;
#define ZERO 0.0f
#define HALF 0.5f
#define ONE  1.0f
#define GRAVITATIONAL_CONSTANT 9.80f
#endif
#endif

#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef fabs
#define fabs(a) ( (a) < 0 ? -(a) : a)
#endif

#ifndef INT_MIN
#define INT_MAX 2147483647
#define INT_MIN (–2147483647–1)
#endif

#define TWO 2
//#define __OLD_STENCIL__
#define __NEW_STENCIL__

enum boundary
{  REAL_CELL      =  1,         //  Denotes cell type of real cell.
   LEFT_BOUNDARY  = -1,         //  Denotes left boundary ghost cell.
   RIGHT_BOUNDARY = -2,         //  Denotes right boundary ghost cell.
   BOTTOM_BOUNDARY= -3,         //  Denotes bottom boundary ghost cell.
   TOP_BOUNDARY   = -4,         //  Denotes top boundary ghost cell.
   FRONT_BOUNDARY = -5,         //  Denotes front boundary ghost cell.
   BACK_BOUNDARY  = -6 };       //  Denotes back boundary ghost cell.

enum orientation
{  SW,                          //  SW quadrant.
   NW,                          //  NW quadrant.
   NE,                          //  NE quadrant.
   SE };                        //  SE quadrant.

int is_lower_left(int i, int j)  { return(i % 2 == 0 && j % 2 == 0); }
int is_lower_right(int i, int j) { return(i % 2 == 1 && j % 2 == 0); }
int is_upper_left(int i, int j)  { return(i % 2 == 0 && j % 2 == 1); }
int is_upper_right(int i, int j) { return(i % 2 == 1 && j % 2 == 1); }

#define hashval(j,i) hash[(j)*imaxsize+(i)]
#define hashval_local(j,i) hash[(j)*(imaxsize-iminsize)+(i)]

__constant ulong prime = 4294967291;
__constant uint hash_jump_prime = 41;

void write_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const uint  giX,
            const ulong hashkey,
   __global       int   *hash)
{
   uint hashloc;
   int icount = 0;
   uint jump;
   int old_key;
   int MaxTries = 1000;

#ifndef __APPLE_CC__
   switch (hash_method) {
   case -1:
#endif
      hash[hashkey] = giX;
#ifndef __APPLE_CC__
      break;
   case 0:
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (int icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc++;
         hashloc %= hash_table_size;
      
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   case 1:
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (int icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc+=(icount*icount);
         hashloc %= hash_table_size;
      
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   case 2:
      jump = 1+hashkey%hash_jump_prime;
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (int icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc += (icount*jump);
         hashloc %= hash_table_size;
      
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   }
#endif
}

int read_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const ulong hashkey,
   __global const int   *hash)
{
   int hashval = -1;
   uint hashloc;
   int icount = 0;
   uint jump;

#ifndef __APPLE_CC__
   switch (hash_method) {
   case -1:
#endif
      return(hash[hashkey]);
#ifndef __APPLE_CC__
      break;
   case 0:
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc++,hashloc %= hash_table_size);
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   case 1:
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc %= hash_table_size){
         icount++;
      }
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   case 2:
      jump = 1+hashkey%hash_jump_prime;
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc %= hash_table_size){
         icount++;
      }
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   }
#endif
}

__kernel void hash_init_cl(
                          const int isize,              // 0
                 __global       int *hash)              // 1
{
   const uint giX  = get_global_id(0);

   if (giX >= isize) return;

   hash[giX] = -1;
}

__kernel void hash_init_corners_cl(
                          const int  isize,      // 0
                          const int  levmx,      // 1
                          const int  imax,       // 2
                          const int  jmax,       // 3
                 __global const int  *levtable,  // 4
                 __global const int  *corners_i, // 5
                 __global const int  *corners_j, // 6
                 __global const int4 *sizes,     // 7
                 __global       int  *hash)      // 8
{
   const uint giX  = get_global_id(0);

   if (giX >= isize) return;

   // Expand size by 2*coarse_cells
   sizes[0].s0 = max(sizes[0].s0-2*levtable[levmx],0);
   sizes[0].s1 = min(sizes[0].s1+2*levtable[levmx],(imax+1)*levtable[levmx]);
   sizes[0].s2 = max(sizes[0].s2-2*levtable[levmx],0);
   sizes[0].s3 = min(sizes[0].s3+2*levtable[levmx],(jmax+1)*levtable[levmx]);

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int ii = corners_i[giX]-iminsize;
   int jj = corners_j[giX]-jminsize;

   if (ii >= 0 && ii < imaxsize-iminsize &&
       jj >= 0 && jj < jmaxsize-jminsize) {
      hash[jj*(imaxsize-iminsize)+ii] = INT_MIN;
   }

}

__kernel void hash_setup_cl(
                          const int  isize,            // 0
                          const int  levmx,            // 1
                          const int  imax,             // 2
                          const int  jmax,             // 3
                          const int  imaxsize,         // 4
                          const int  hash_method,      // 5
                          const ulong hash_table_size, // 6
                          const ulong AA,              // 7
                          const ulong BB,              // 8
                 __global const int  *levtable,        // 9
                 __global const int  *lev_ibeg,        // 10
                 __global const int  *lev_iend,        // 11
                 __global const int  *lev_jbeg,        // 12
                 __global const int  *lev_jend,        // 13
                 __global const int  *level,           // 14
                 __global const int  *i,               // 15
                 __global const int  *j,               // 16
                 __global       int  *hash)            // 17
      //do_compact_hash = do_compact_hash_in;
      //hash_method     = hash_method_in;
      //hash_table_size = hash_table_size_in;
      //AA              = AA_in;
      //BB              = BB_in;
{

   const uint giX  = get_global_id(0);

   if (giX >= isize) return;

   int lev = level[giX];
   int ii = i[giX];
   int jj = j[giX];
   int levmult = levtable[levmx-lev];

#if HASH_SETUP_OPT_LEVEL == 0

   int iimin =  ii   *levmult;
   int iimax = (ii+1)*levmult;
   int jjmin =  jj   *levmult;
   int jjmax = (jj+1)*levmult;

   if      (ii < lev_ibeg[lev]) iimin = 0;                        // left boundary
   else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]; // right boundary
   else if (jj < lev_jbeg[lev]) jjmin = 0;                        // bottom boundary
   else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]; // top boundary

   for (   int jjj = jjmin; jjj < jjmax; jjj++) {
      for (int iii = iimin; iii < iimax; iii++) {
         hashval(jjj, iii) = giX;
      }
   }

#elif HASH_SETUP_OPT_LEVEL == 1

   int iimin =  ii   *levmult;
   int iimax = (ii+1)*levmult;
   int jjmin =  jj   *levmult;
   int jjmax = (jj+1)*levmult;

   if      (ii < lev_ibeg[lev]) iimin = 0;                        // left boundary
   else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]; // right boundary
   else if (jj < lev_jbeg[lev]) jjmin = 0;                        // bottom boundary
   else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]; // top boundary

   for (int iii = iimin; iii < iimax; iii++) {
      hashval(jjmin,   iii) = giX;
      hashval(jjmax-1, iii) = giX;
   }
   for (int jjj = jjmin+1; jjj < jjmax-1; jjj++) {
      hashval(jjj, iimin) = giX;
      hashval(jjj, iimax-1) = giX;
   }

#elif HASH_SETUP_OPT_LEVEL == 2
   if (lev > 0 && (ii < lev_ibeg[lev] || ii > lev_iend[lev] || jj < lev_jbeg[lev] || jj > lev_jend[lev]) ) {

      int iimin =  ii   *levmult;
      int iimax = (ii+1)*levmult;
      int jjmin =  jj   *levmult;
      int jjmax = (jj+1)*levmult;

      if      (ii < lev_ibeg[lev]) iimin = 0;                        // left boundary
      else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]; // right boundary
      else if (jj < lev_jbeg[lev]) jjmin = 0;                        // bottom boundary
      else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]; // top boundary

      for (   int jjj = jjmin; jjj < jjmax; jjj++) {
         for (int iii = iimin; iii < iimax; iii++) {
            hashval(jjj, iii) = giX;
         }
      }
      return;
   }

   if(lev == levmx) {
      hashval(jj,ii) = giX;
      return;
   }

   jj *= levmult;
   ii *= levmult;
   hashval(jj,ii) = giX;

   ii += levmult/2;
   hashval(jj,ii) = giX;
   if(levmult > 2) {
      ii += levmult/2 - 1;
      hashval(jj,ii) = giX;
      ii -= levmult/2 - 1;
   }
   ii -= levmult/2;
   jj += levmult/2;
   hashval(jj,ii) = giX;
   ii += levmult - 1;
   hashval(jj,ii) = giX;

   if(levmult > 2) {
      ii -= levmult - 1;
      jj += levmult/2 - 1;
      hashval(jj,ii) = giX;
      ii += levmult/2;
      hashval(jj,ii) = giX;
   }
#elif HASH_SETUP_OPT_LEVEL == 3
   jj *= levmult;
   ii *= levmult;
   hashval(jj,ii) = giX;
#elif HASH_SETUP_OPT_LEVEL == 4
   jj *= levmult;
   ii *= levmult;
   write_hash(hash_method, hash_table_size, AA, BB, giX, jj*imaxsize+ii, hash);
#endif
}

__kernel void hash_setup_local_cl(
                          const int  isize,     // 0
                          const int  levmx,     // 1
                          const int  imax,      // 2
                          const int  jmax,      // 3
                          const int  noffset,   // 4
                 __global const int4 *sizes,    // 5
                 __global const int  *levtable, // 6
                 __global const int  *lev_ibeg, // 7
                 __global const int  *lev_iend, // 8
                 __global const int  *lev_jbeg, // 9
                 __global const int  *lev_jend, // 10
                 __global const int  *level,    // 11
                 __global const int  *i,        // 12
                 __global const int  *j,        // 13
                 __global       int  *hash)     // 14
{

   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   //int jmaxsize = sizes[0].s3;

   int ii = i[giX];
   int jj = j[giX];
   int lev = level[giX];
   int levmult = levtable[levmx-lev];
   int cell_number = giX+noffset;

#if HASH_SETUP_OPT_LEVEL == 0
/* Original Hash Setup */

   int iimin =  ii   *levmult-iminsize;
   int iimax = (ii+1)*levmult-iminsize;
   int jjmin =  jj   *levmult-jminsize;
   int jjmax = (jj+1)*levmult-jminsize;

   if      (ii < lev_ibeg[lev]) iimin = 0;                                 // left boundary
   else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]-iminsize; // right boundary
   else if (jj < lev_jbeg[lev]) jjmin = 0;                                 // bottom boundary
   else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]-jminsize; // top boundary

   for (   int jjj = jjmin; jjj < jjmax; jjj++) {
      for (int iii = iimin; iii < iimax; iii++) {
         hashval_local(jjj, iii) = cell_number;
      }
   }

#elif HASH_SETUP_OPT_LEVEL == 1
/* Just the outer cells */

   int iimin =  ii   *levmult-iminsize;
   int iimax = (ii+1)*levmult-iminsize;
   int jjmin =  jj   *levmult-jminsize;
   int jjmax = (jj+1)*levmult-jminsize;

   if      (ii < lev_ibeg[lev]) iimin = 0;                                 // left boundary
   else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]-iminsize; // right boundary
   else if (jj < lev_jbeg[lev]) jjmin = 0;                                 // bottom boundary
   else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]-jminsize; // top boundary

   for (int iii = iimin; iii < iimax; iii++) {
      hashval_local(jjmin,   iii) = cell_number;
      hashval_local(jjmax-1, iii) = cell_number;
   }
   for (int jjj = jjmin+1; jjj < jjmax-1; jjj++) {
      hashval_local(jjj, iimin) = cell_number;
      hashval_local(jjj, iimax-1) = cell_number;
   }

#elif HASH_SETUP_OPT_LEVEL >= 2
/* Optimized Hash Setup */
   if (lev > 0 && (ii < lev_ibeg[lev] || ii > lev_iend[lev] || jj < lev_jbeg[lev] || jj > lev_jend[lev]) ) {

      int iimin =  ii   *levmult-iminsize;
      int iimax = (ii+1)*levmult-iminsize;
      int jjmin =  jj   *levmult-jminsize;
      int jjmax = (jj+1)*levmult-jminsize;

      if      (ii < lev_ibeg[lev]) iimin = 0;                                 // left boundary
      else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]-iminsize; // right boundary
      else if (jj < lev_jbeg[lev]) jjmin = 0;                                 // bottom boundary
      else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]-jminsize; // top boundary

      for (   int jjj = jjmin; jjj < jjmax; jjj++) {
         for (int iii = iimin; iii < iimax; iii++) {
            hashval_local(jjj, iii) = cell_number;
         }
      }
      return;
   }

   if(lev == levmx) {
      hashval_local(jj-jminsize,ii-iminsize) = cell_number;
      return;
   }

   jj = jj*levmult - jminsize;
   ii = ii*levmult - iminsize;

   hashval_local(jj,ii) = cell_number; // lower left corner

   ii += levmult/2;
   hashval_local(jj,ii) = cell_number; // lower boundary mid-point
   if(levmult > 2) {
      ii += levmult/2 - 1;
      hashval_local(jj,ii) = cell_number; // lower right corner
      ii -= levmult/2 - 1;
   }
   ii -= levmult/2;
   jj += levmult/2;
   hashval_local(jj,ii) = cell_number; // left boundary mid-point
   ii += levmult - 1;
   hashval_local(jj,ii) = cell_number; // right boundary mid-point

   if(levmult > 2) {
      ii -= levmult - 1;
      jj += levmult/2 - 1;
      hashval_local(jj,ii) = cell_number; // upper left boundary
      ii += levmult/2;
      hashval_local(jj,ii) = cell_number; // upper boundary mid-point
   }
#endif

}

__kernel void calc_neighbors_cl(
                          const int  isize,            // 0 
                          const int  levmx,            // 1 
                          const int  imax,             // 2 
                          const int  jmax,             // 3 
                          const int  imaxsize,         // 4 
                          const int  jmaxsize,         // 5 
                          const int  hash_method,      // 6
                          const ulong hash_table_size, // 7
                          const ulong AA,              // 8
                          const ulong BB,              // 9
                 __global const int  *levtable,        // 10
                 __global const int  *level,           // 11
                 __global const int  *i,               // 12
                 __global const int  *j,               // 13
                 __global       int  *nlft,            // 14
                 __global       int  *nrht,            // 15
                 __global       int  *nbot,            // 16
                 __global       int  *ntop,            // 17
                 __global const int  *hash)            // 18
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int ii, jj, iii, jjj, lev, levmult;

   ii = i[giX];
   jj = j[giX];
   lev = level[giX];
   levmult = levtable[levmx-lev];

#if HASH_SETUP_OPT_LEVEL <= 2
   int nlftval = hashval(      jj   *levmult               , max(  ii   *levmult-1, 0         ));
   int nrhtval = hashval(      jj   *levmult               , min( (ii+1)*levmult,   imaxsize-1));
   int nbotval = hashval(max(  jj   *levmult-1, 0)         ,       ii   *levmult               );
   int ntopval = hashval(min( (jj+1)*levmult,   jmaxsize-1),       ii   *levmult               );

   // Handles the four boundary corners
   if (nlftval < 0){
      iii = max( ii*levmult-1, 0);
      jjj = jj*levmult;
      if ( (jjj < 1*levtable[levmx] || jjj > (jmax-1)*levtable[levmx] ) &&
            iii < 1*levtable[levmx] ) iii = ii*levmult;
      nlftval = hashval(jjj, iii);
   }
   if (nrhtval < 0){
      iii = min( (ii+1)*levmult, imaxsize-1);
      jjj = jj*levmult;
      if ( (jjj < 1*levtable[levmx] || jjj > jmax*levtable[levmx]-1 ) &&
            iii > imax*levtable[levmx]-1 ) iii = ii*levmult;
      nrhtval = hashval(jjj, iii);
   }
   if (nbotval < 0) {
      iii = ii*levmult;
      jjj = max( jj*levmult-1, 0);
      if ( (iii < 1*levtable[levmx] || iii > (imax-1)*levtable[levmx] ) &&
            jjj < 1*levtable[levmx] ) jjj = jj*levmult;
      nbotval = hashval(jjj, iii);
   }
   if (ntopval < 0) {
      iii = ii*levmult;
      jjj = min( (jj+1)*levmult, jmaxsize-1);
      if ( (iii < 1*levtable[levmx] || iii > imax*levtable[levmx]-1 ) &&
            jjj > jmax*levtable[levmx]-1 ) jjj = jj*levmult;
      ntopval = hashval(jjj, iii);
   }
#elif HASH_SETUP_OPT_LEVEL == 3

   int iicur = ii*levmult;
   int iilft = max( (ii-1)*levmult, 0         );   
   int iirht = min( (ii+1)*levmult, imaxsize-1);
   int jjcur = jj*levmult;
   int jjbot = max( (jj-1)*levmult, 0         );   
   int jjtop = min( (jj+1)*levmult, jmaxsize-1);

   int nlftval = -1;
   int nrhtval = -1;
   int nbotval = -1;
   int ntopval = -1;

   // Taking care of boundary cells
   // Force each boundary cell to point to itself on its boundary direction
   if (iicur <    1*levtable[levmx]  ) nlftval = giX;
   if (jjcur <    1*levtable[levmx]  ) nbotval = giX;
   if (iicur > imax*levtable[levmx]-1) nrhtval = giX;
   if (jjcur > jmax*levtable[levmx]-1) ntopval = giX;
   // Boundary cells next to corner boundary need special checks
   if (iicur ==    1*levtable[levmx] &&  (jjcur < 1*levtable[levmx] || jjcur > (jmax-1)*levtable[levmx] ) ) nlftval = giX;
   if (jjcur ==    1*levtable[levmx] &&  (iicur < 1*levtable[levmx] || iicur > (imax-1)*levtable[levmx] ) ) nbotval = giX;
   if (iirht == imax*levtable[levmx] &&  (jjcur < 1*levtable[levmx] || jjcur > (jmax-1)*levtable[levmx] ) ) nrhtval = giX;
   if (jjtop == jmax*levtable[levmx] &&  (iicur < 1*levtable[levmx] || iicur > (imax-1)*levtable[levmx] ) ) ntopval = giX;

   // need to check for finer neighbor first
   if (lev != levmx) {
      int iilftfiner = iicur-(iicur-iilft)/2;
      int jjbotfiner = jjcur-(jjcur-jjbot)/2;
      if (nlftval < 0) nlftval = hashval(jjcur,iilftfiner);
      if (nbotval < 0) nbotval = hashval(jjbotfiner,iicur);
   }    

   // same size neighbor
   if (nlftval < 0) nlftval = hashval(jjcur,iilft);
   if (nrhtval < 0) nrhtval = hashval(jjcur,iirht);
   if (nbotval < 0) nbotval = hashval(jjbot,iicur);
   if (ntopval < 0) ntopval = hashval(jjtop,iicur);

   // Now we need to take care of special case where bottom and left boundary need adjustment since
   // expected cell doesn't exist on these boundaries if it is finer than current cell
   if (lev != levmx) {
      if (jjcur < 1*levtable[levmx]) {
         if (nrhtval < 0) {
            int jjtopfiner = (jjcur+jjtop)/2;
            nrhtval = hashval(jjtopfiner,iirht);
         }
         if (nlftval < 0) {
            int iilftfiner = iicur-(iicur-iilft)/2;
            int jjtopfiner = (jjcur+jjtop)/2;
            nlftval = hashval(jjtopfiner,iilftfiner);
         }
      }

      if (iicur < 1*levtable[levmx]) {
         if (ntopval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            ntopval = hashval(jjtop,iirhtfiner);
         }
         if (nbotval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            int jjbotfiner = jjcur-(jjcur-jjbot)/2;
            nbotval = hashval(jjbotfiner,iirhtfiner);
         }
      }
   }

   // coarser neighbor
   if (lev != 0){
      if (nlftval < 0) { 
         iilft -= iicur-iilft;
         int jjlft = (jj/2)*2*levmult;
         nlftval = hashval(jjlft,iilft);
      }    
      if (nrhtval < 0) {
         int jjrht = (jj/2)*2*levmult;
         nrhtval = hashval(jjrht,iirht);
      }    
      if (nbotval < 0) { 
         jjbot -= jjcur-jjbot;
         int iibot = (ii/2)*2*levmult;
         nbotval = hashval(jjbot,iibot);
      }                
      if (ntopval < 0) {
         int iitop = (ii/2)*2*levmult;
         ntopval = hashval(jjtop,iitop);
      }    
   }
#elif HASH_SETUP_OPT_LEVEL == 4

   int iicur = ii*levmult;
   int iilft = max( (ii-1)*levmult, 0         );   
   int iirht = min( (ii+1)*levmult, imaxsize-1);
   int jjcur = jj*levmult;
   int jjbot = max( (jj-1)*levmult, 0         );   
   int jjtop = min( (jj+1)*levmult, jmaxsize-1);

   int nlftval = -1;
   int nrhtval = -1;
   int nbotval = -1;
   int ntopval = -1;

   // Taking care of boundary cells
   // Force each boundary cell to point to itself on its boundary direction
   if (iicur <    1*levtable[levmx]  ) nlftval = giX;
   if (jjcur <    1*levtable[levmx]  ) nbotval = giX;
   if (iicur > imax*levtable[levmx]-1) nrhtval = giX;
   if (jjcur > jmax*levtable[levmx]-1) ntopval = giX;
   // Boundary cells next to corner boundary need special checks
   if (iicur ==    1*levtable[levmx] &&  (jjcur < 1*levtable[levmx] || jjcur > (jmax-1)*levtable[levmx] ) ) nlftval = giX;
   if (jjcur ==    1*levtable[levmx] &&  (iicur < 1*levtable[levmx] || iicur > (imax-1)*levtable[levmx] ) ) nbotval = giX;
   if (iirht == imax*levtable[levmx] &&  (jjcur < 1*levtable[levmx] || jjcur > (jmax-1)*levtable[levmx] ) ) nrhtval = giX;
   if (jjtop == jmax*levtable[levmx] &&  (iicur < 1*levtable[levmx] || iicur > (imax-1)*levtable[levmx] ) ) ntopval = giX;

   // need to check for finer neighbor first
   if (lev != levmx) {
      int iilftfiner = iicur-(iicur-iilft)/2;
      int jjbotfiner = jjcur-(jjcur-jjbot)/2;
      if (nlftval < 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur     *imaxsize+iilftfiner, hash);
      if (nbotval < 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*imaxsize+iicur,      hash);
   }    

   // same size neighbor
   if (nlftval < 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*imaxsize+iilft, hash);
   if (nrhtval < 0) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*imaxsize+iirht, hash);
   if (nbotval < 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbot*imaxsize+iicur, hash);
   if (ntopval < 0) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*imaxsize+iicur, hash);

   // Now we need to take care of special case where bottom and left boundary need adjustment since
   // expected cell doesn't exist on these boundaries if it is finer than current cell
   if (lev != levmx) {
      if (jjcur < 1*levtable[levmx]) {
         if (nrhtval < 0) {
            int jjtopfiner = (jjcur+jjtop)/2;
            nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*imaxsize+iirht, hash);
         }
         if (nlftval < 0) {
            int iilftfiner = iicur-(iicur-iilft)/2;
            int jjtopfiner = (jjcur+jjtop)/2;
            nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*imaxsize+iilftfiner, hash);
         }
      }

      if (iicur < 1*levtable[levmx]) {
         if (ntopval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*imaxsize+iirhtfiner, hash);
         }
         if (nbotval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            int jjbotfiner = jjcur-(jjcur-jjbot)/2;
            nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*imaxsize+iirhtfiner, hash);
         }
      }
   }

   // coarser neighbor
   if (lev != 0){
      if (nlftval < 0) { 
         iilft -= iicur-iilft;
         int jjlft = (jj/2)*2*levmult;
         nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjlft*imaxsize+iilft, hash);
      }    
      if (nrhtval < 0) {
         int jjrht = (jj/2)*2*levmult;
         nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjrht*imaxsize+iirht, hash);
      }    
      if (nbotval < 0) { 
         jjbot -= jjcur-jjbot;
         int iibot = (ii/2)*2*levmult;
         nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbot*imaxsize+iibot, hash);
      }                
      if (ntopval < 0) {
         int iitop = (ii/2)*2*levmult;
         ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*imaxsize+iitop, hash);
      }    
   }
#endif

   nlft[giX] = nlftval;
   nrht[giX] = nrhtval;
   nbot[giX] = nbotval;
   ntop[giX] = ntopval;
}

__kernel void calc_neighbors_local_cl(
                          const int  isize,     // 0 
                          const int  levmx,     // 1 
                          const int  imax,      // 2 
                          const int  jmax,      // 3 
                 __global const int4 *sizes,    // 4 
                 __global const int  *levtable, // 5 
                 __global const int  *level,    // 6 
                 __global const int  *i,        // 7 
                 __global const int  *j,        // 8 
                 __global       int  *nlft,     // 9 
                 __global       int  *nrht,     // 10
                 __global       int  *nbot,     // 11
                 __global       int  *ntop,     // 12
                 __global const int  *hash)     // 13
{
                
   const uint giX  = get_global_id(0);

   if (giX >= isize) return;

   int iii, jjj;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   //int jmaxsize = sizes[0].s3;

   int imaxcalc = (imax+1)*levtable[levmx];
   int jmaxcalc = (jmax+1)*levtable[levmx];

   int ii = i[giX];
   int jj = j[giX];
   int lev = level[giX];
   int levmult = levtable[levmx-lev];
   int nlftval = hash[ ( (      jj   *levmult               )-jminsize) *(imaxsize-iminsize) + ( (max(  ii   *levmult-1, 0         ))-iminsize)];
   int nrhtval = hash[ ( (      jj   *levmult               )-jminsize) *(imaxsize-iminsize) + ( (min( (ii+1)*levmult,   imaxcalc-1))-iminsize)];
   int nbotval = hash[ ( (max(  jj   *levmult-1, 0)         )-jminsize) *(imaxsize-iminsize) + ( (      ii   *levmult               )-iminsize)];
   int ntopval = hash[ ( (min( (jj+1)*levmult,   jmaxcalc-1))-jminsize) *(imaxsize-iminsize) + ( (      ii   *levmult               )-iminsize)];

   // Handles the four boundary corners
   if (nlftval == INT_MIN){
      iii = ii*levmult;
      jjj = jj*levmult;
      nlftval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (nrhtval == INT_MIN){
      iii = ii*levmult;
      jjj = jj*levmult;
      nrhtval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (nbotval == INT_MIN) {
      iii = ii*levmult;
      jjj = jj*levmult;
      nbotval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (ntopval == INT_MIN) {
      iii = ii*levmult;
      jjj = jj*levmult;
      ntopval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }

   nlft[giX] = nlftval;
   nrht[giX] = nrhtval;
   nbot[giX] = nbotval;
   ntop[giX] = ntopval;
}

__kernel void calc_border_cells_cl(
                          const int isize,      // 0
                          const int noffset,    // 1
                 __global const int *nlft,      // 2
                 __global const int *nrht,      // 3
                 __global const int *nbot,      // 4
                 __global const int *ntop,      // 5
                 __global const int *level,     // 6
                 __global       uint *border_cell_out //7
                                  )
{
   const unsigned int giX = get_global_id(0);

   if (giX >= isize) return;

   int border_cell = 0;

   if (nlft[giX] == -1 || (level[nlft[giX]-noffset] > level[giX] && ntop[nlft[giX]-noffset] == -1) ){
      border_cell |= 0x0001;
   }
   if (nrht[giX] == -1 || (level[nrht[giX]-noffset] > level[giX] && ntop[nrht[giX]-noffset] == -1) ){
      border_cell |= 0x0002;
   }
   if (nbot[giX] == -1 || (level[nbot[giX]-noffset] > level[giX] && nrht[nbot[giX]-noffset] == -1) ) {
      border_cell |= 0x0004;
   }
   // top neighbor is undefined -- or -- if top is at finer level check top right for undefined
   if (ntop[giX] == -1 || ( level[ntop[giX]-noffset] > level[giX] && nrht[ntop[giX]-noffset] == -1) ) {
      border_cell |= 0x0008;
   }

   border_cell_out[giX] = border_cell;
}

__kernel void calc_border_cells2_cl(
                          const int   isize,            // 0
                          const int   noffset,          // 1
                 __global const int   *nlft,            // 2
                 __global const int   *nrht,            // 3
                 __global const int   *nbot,            // 4
                 __global const int   *ntop,            // 5
                 __global const int   *level,           // 6
                 __global const uint  *border_cell_in,  // 7
                 __global       uint  *border_cell_out, // 8
                 __global       uint  *ioffset,         // 9
                 __global       int   *nbsize,          // 10
                 __local        int   *itile)           // 11
{
   const uint giX = get_global_id(0);
   const uint tiX = get_local_id(0);
   const uint group_id = get_group_id(0);
   const uint ntX = get_local_size(0);

   itile[tiX] = 0;

   if (giX >= isize) return;

   int border_cell = border_cell_in[giX];

   if (border_cell == 0) {

      int nl = nlft[giX]-noffset;
      if (nl >= 0 && nl < isize) {
         if ((border_cell_in[nl] & 0x0001) == 0x0001) {
            border_cell |= 0x0016;
         } else if (level[nl] > level[giX]){
            int ntl = ntop[nl]-noffset;
            if (ntl >= 0 && ntl < isize && (border_cell_in[ntl] & 0x0001) == 0x0001) {
               border_cell |= 0x0016;
            }
         }
      }
      int nr = nrht[giX]-noffset;
      if (nr >= 0 && nr < isize) {
         if ((border_cell_in[nrht[giX]-noffset] & 0x0002) == 0x0002) {
            border_cell |= 0x0032;
         } else if (level[nr] > level[giX]){
            int ntr = ntop[nr]-noffset;
            if (ntr >= 0 && ntr < isize && (border_cell_in[ntr] & 0x0002) == 0x0002) {
               border_cell |= 0x0032;
            }
         }
      }
      int nb = nbot[giX]-noffset;
      if (nb >= 0 && nb < isize) {
         if ((border_cell_in[nb] & 0x0004) == 0x0004) {
            border_cell |= 0x0064;
         } else if (level[nb] > level[giX]){
            int nrb = nrht[nb]-noffset;
            if (nrb >= 0 && nrb < isize && (border_cell_in[nrb] & 0x0004) == 0x0004) {
               border_cell |= 0x0064;
            }
         }
      }
      int nt = ntop[giX]-noffset;
      if (nt >= 0 && nt < isize) {
         if ((border_cell_in[nt] & 0x0008) == 0x0008) {
            border_cell |= 0x0128;
         } else if (level[nt] > level[giX]){
            int nrt = nrht[nt]-noffset;
            if (nrt >= 0 && nrt < isize && (border_cell_in[nrt] & 0x0008) == 0x0008) {
               border_cell |= 0x0128;
            }
         }
      }
   }

// Setting up for scan of border cells
   if (border_cell) itile[tiX] = 1;

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset = ntX >> 1; offset > 32; offset >>= 1) { 
      if (tiX < offset) {
         itile[tiX] += itile[tiX+offset]; 
      }    
      barrier(CLK_LOCAL_MEM_FENCE);
   }    

   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  itile[tiX] += itile[tiX+32];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+1]; }

   if (tiX == 0) { 
     ioffset[group_id] = itile[0];
     (*nbsize) = itile[0];
   }    

   border_cell_out[giX] = border_cell;
}

inline uint scan_warp_exclusive(__local volatile uint *input, const uint idx, const uint lane) {
    if (lane > 0 ) input[idx] += input[idx - 1];
    if (lane > 1 ) input[idx] += input[idx - 2];
    if (lane > 3 ) input[idx] += input[idx - 4];
    if (lane > 7 ) input[idx] += input[idx - 8];
    if (lane > 15) input[idx] += input[idx - 16];

    return (lane > 0) ? input[idx-1] : 0;
}

inline uint scan_warp_inclusive(__local volatile uint *input, const uint idx, const uint lane) {
    if (lane > 0 ) input[idx] += input[idx - 1];
    if (lane > 1 ) input[idx] += input[idx - 2];
    if (lane > 3 ) input[idx] += input[idx - 4];
    if (lane > 7 ) input[idx] += input[idx - 8];
    if (lane > 15) input[idx] += input[idx - 16];

    return input[idx];
}

inline uint scan_workgroup_exclusive(
    __local uint* itile,
    const uint tiX,
    const uint lane,
    const uint warpID) {

    // Step 1: scan each warp
    uint val = scan_warp_exclusive(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 2: Collect per-warp sums
    if (lane == 31) itile[warpID] = itile[tiX];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 3: Use 1st warp to scan per-warp sums
    if (warpID == 0) scan_warp_inclusive(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 4: Accumulate results from Steps 1 and 3
    if (warpID > 0) val += itile[warpID-1];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 6: Write and return the final result
    itile[tiX] = val;
    barrier(CLK_LOCAL_MEM_FENCE);

    return val;
}

__kernel void finish_scan_cl(
                          const int   size,       // 0
                 __global       uint  *ioffset,   // 1
                 __global       int   *nbsize,    // 2
                 __local        uint  *itile)     // 3
{
   const uint tiX = get_local_id(0);
   const uint group_id = get_group_id(0);
   const uint ntX = get_local_size(0);

   const uint lane = tiX & 31; 
   const uint warpID = tiX >> 5;
   const uint EPT = (size+ntX-1)/ntX; //elements_per_thread;

   uint reduceValue = 0;

//  #pragma unroll 4
   for(uint i = 0; i < EPT; ++i)
   {
      uint offsetIdx = i * ntX + tiX;

#ifdef IS_NVIDIA
//    if (offsetIdx >= size) return;
#endif

      // Step 1: Read ntX elements from global (off-chip) memory to local memory (on-chip)
      uint input = 0;
      if (offsetIdx < size) input = ioffset[offsetIdx];    
      itile[tiX] = input;    
      barrier(CLK_LOCAL_MEM_FENCE);

      // Step 2: Perform scan on ntX elements
      uint val = scan_workgroup_exclusive(itile, tiX, lane, warpID);

      // Step 3: Propagate reduced result from previous block of ntX elements
      val += reduceValue;

      // Step 4: Write out data to global memory
      if (offsetIdx < size) ioffset[offsetIdx] = val;

      // Step 5: Choose reduced value for next iteration
      if (tiX == (ntX-1)) itile[tiX] = input + val;
      barrier(CLK_LOCAL_MEM_FENCE);

      reduceValue = itile[ntX-1];
      barrier(CLK_LOCAL_MEM_FENCE);
   }
   (*nbsize) = itile[ntX-1];
}

__kernel void get_border_data_cl(
                          const int  isize,              // 0
                          const int  noffset,            // 1
                 __global       uint *ioffset,           // 2
                 __global const int  *border_cell,       // 3
                 __global const int  *i,                 // 4
                 __global const int  *j,                 // 5
                 __global const int  *level,             // 6
                 __global       int  *border_cell_i,     // 8
                 __global       int  *border_cell_j,     // 9
                 __global       int  *border_cell_level, // 10
                 __global       int  *border_cell_num,   // 11
                 __local        uint *itile)             // 12
{
   const uint giX = get_global_id(0);
   const uint tiX = get_local_id(0);
   const uint group_id = get_group_id(0);

   const uint lane   = tiX & 31; 
   const uint warpid = tiX >> 5;

   int cell_num = giX-noffset;

   // Step 1: load global data into tile
   int temp_val = 0;
   if (giX < isize) temp_val = border_cell[giX];
   itile[tiX] = 0;
   if (temp_val > 0) itile[tiX] = 1;
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 2: scan each warp
   uint val = scan_warp_exclusive(itile, tiX, lane);
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 3: Collect per-warp sums
   if (lane == 31) itile[warpid] = itile[tiX];
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 4: Use 1st warp to scan per-warp sums
   if (warpid == 0) scan_warp_inclusive(itile, tiX, lane);
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 5: Accumulate results from Steps 2 and 4
   if (warpid > 0) val += itile[warpid-1];
   barrier(CLK_LOCAL_MEM_FENCE);

   if (giX >= isize || temp_val <= 0) return;

   // Step 6: Write and return the final result
   //itile[tiX] = val;
   //barrier(CLK_LOCAL_MEM_FENCE);

   val += ioffset[group_id];   //index to write to for each thread

   border_cell_i[val]     = i[giX];
   border_cell_j[val]     = j[giX];
   border_cell_level[val] = level[giX];
   border_cell_num[val]   = giX+noffset;
}


__kernel void calc_layer1_cl (
                          const int  isize,               // 0 
                          const int  levmx,               // 1 
                          const int  imax,                // 2 
                          const int  jmax,                // 3 
                 __global const int4 *sizes,              // 4 
                 __global const int  *levtable,           // 5 
                 __global const int  *border_cell_i,      // 6
                 __global const int  *border_cell_j,      // 7
                 __global const int  *border_cell_level,  // 8
                 __global       int  *border_cell_needed, // 9
                 __global const int  *hash)               // 10
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int imaxcalc = (imax+1)*levtable[levmx];
   int jmaxcalc = (jmax+1)*levtable[levmx];

   int ii = border_cell_i[giX];
   int jj = border_cell_j[giX];
   int lev = border_cell_level[giX];
   int levmult = levtable[levmx-lev];

   int iborder = 0;

   if (max(ii*levmult-1, 0)-iminsize >= 0 && min((jj+1)*levmult -1,jmaxcalc-1) < jmaxsize) {  // Test for cell to left
      if (hash[(    (jj   *levmult)              -jminsize)*(imaxsize-iminsize)+(max(ii*levmult-1, 0)-iminsize)] >= 0 ||
          hash[(min(jj*levmult + levmult/2,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(max(ii*levmult-1, 0)-iminsize)] >= 0 ) {
         iborder |= 0x0001;
      }
   }
   if (min( (ii+1)*levmult,   imaxcalc-1) < imaxsize && min((jj+1)*levmult -1,jmaxcalc-1) < jmaxsize){ // Test for cell to right
      if (hash[(    (jj   *levmult)              -jminsize)*(imaxsize-iminsize)+(min( (ii+1)*levmult,imaxcalc-1)-iminsize)] >= 0 ||
          hash[(min((jj+1)*levmult -1,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min( (ii+1)*levmult,imaxcalc-1)-iminsize)] >= 0 ) {
         iborder |= 0x0002;
      }
   }
   if (max(jj*levmult-1, 0)-jminsize >= 0 && min((ii+1)*levmult -1,imaxcalc-1) < imaxsize){ // Test for cell to bottom
      if (hash[((max(jj*levmult-1, 0) )-jminsize)*(imaxsize-iminsize)+(    (ii   *levmult)              -iminsize)] >= 0 ||
          hash[((max(jj*levmult-1, 0) )-jminsize)*(imaxsize-iminsize)+(min(ii*levmult+levmult/2,imaxcalc-1)-iminsize)] >= 0 ) {
         iborder |= 0x0004;
      }
   }
   if ((min( (jj+1)*levmult,   jmaxcalc-1)) < jmaxsize && min((ii+1)*levmult -1,imaxcalc-1) < imaxsize){ // Test for cell to top
      if (hash[((min( (jj+1)*levmult, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+(    (ii   *levmult)              -iminsize)] >= 0 ||
          hash[((min( (jj+1)*levmult, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+(min((ii+1)*levmult -1,imaxcalc-1)-iminsize)] >= 0 ) {
         iborder |= 0x0008;
      }
   }
   if (iborder) border_cell_needed[giX] = iborder;
}

__kernel void calc_layer1_sethash_cl (
                          const int  isize,               // 0 
                          const int  ncells,              // 1 
                          const int  noffset,             // 2 
                          const int  levmx,               // 3 
                 __global const int4 *sizes,              // 4 
                 __global const int  *levtable,           // 5 
                 __global const int  *border_cell_i,      // 6
                 __global const int  *border_cell_j,      // 7
                 __global const int  *border_cell_level,  // 8
                 __global       int  *border_cell_needed, // 9
                 __global       int  *hash)               // 10
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   int iborder = border_cell_needed[giX];

   if (iborder) {
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;

      int ii = border_cell_i[giX];
      int jj = border_cell_j[giX];
      int lev = border_cell_level[giX];
      int levmult = levtable[levmx-lev];

      if (lev == levmx) {
         hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)] = ncells+noffset+giX;
      } else {
         for (int    j = max(jj*levmult-jminsize,0); j < min((jj+1)*levmult,jmaxsize)-jminsize; j++) {
            for (int i = max(ii*levmult-iminsize,0); i < min((ii+1)*levmult,imaxsize)-iminsize; i++) {
               hash[j*(imaxsize-iminsize)+i] = ncells+noffset+giX;
            }
         }
      }
   }
}

__kernel void calc_layer2_cl (
                          const int isize,                   // 0 
                          const int   ncells,                  // 1 
                          const int   noffset,                 // 2 
                          const int   levmx,                   // 3 
                          const int   imax,                    // 4 
                          const int   jmax,                    // 5 
                 __global const int4  *sizes,                  // 6 
                 __global const int   *levtable,               // 7 
                 __global const int   *border_cell_i,          // 8
                 __global const int   *border_cell_j,          // 9
                 __global const int   *border_cell_level,      // 10
                 __global const int   *border_cell_needed,     // 11
                 __global       int   *border_cell_needed_out, // 12
                 __global const int   *hash,                   // 13
                 __global       int   *ioffset,                // 14
                 __global       int   *nbsize_new,             // 15
                 __local        int   *itile)                  // 16
{
   const uint giX = get_global_id(0);
   const uint tiX = get_local_id(0);
   const uint group_id = get_group_id(0);
   const uint ntX = get_local_size(0);

   itile[tiX] = 0;

   if (giX >= isize) return;

   int iborder = border_cell_needed[giX];

   if (iborder <= 0) {
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;

      int imaxcalc = (imax+1)*levtable[levmx];
      int jmaxcalc = (jmax+1)*levtable[levmx];

      int ii = border_cell_i[giX];
      int jj = border_cell_j[giX];
      int lev = border_cell_level[giX];
      int levmult = levtable[levmx-lev];

      if (max(ii*levmult-1, 0)-iminsize >= 0) {  // Test for cell to left
         int nl  = hash[(    (jj   *levmult)            -jminsize)*(imaxsize-iminsize)+(max(ii*levmult-1, 0)-iminsize)];
         if ( nl  >= (int)(ncells+noffset) && (border_cell_needed[nl -ncells-noffset] & 0x0001) == 0x0001) {
            iborder = 0x0001;
         } else if (min((jj+1)*levmult -1,jmaxcalc-1) < jmaxsize) {
            int nlt = hash[(min((jj+1)*levmult -1,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(max(ii*levmult-1, 0)-iminsize)];
            if ( nlt >= (int)(ncells+noffset) && (border_cell_needed[nlt-ncells-noffset] & 0x0001) == 0x0001) iborder = 0x0001;
         }
      }
      if (min( (ii+1)*levmult, imaxcalc-1) < imaxsize){ // Test for cell to right
         int nr  = hash[((jj *levmult)-jminsize)*(imaxsize-iminsize)+((min( (ii+1)*levmult,   imaxcalc-1))-iminsize)];
         if ( nr  >= (int)(ncells+noffset) && (border_cell_needed[nr -ncells-noffset] & 0x0002) == 0x0002) {
            iborder = 0x0002;
         } else if (min((jj+1)*levmult -1,jmaxcalc-1) < jmaxsize) {
            int nrt = hash[(min((jj+1)*levmult -1,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min( (ii+1)*levmult,imaxcalc-1)-iminsize)];
            if ( nrt >= (int)(ncells+noffset) && (border_cell_needed[nrt-ncells-noffset] & 0x0002) == 0x0002) iborder = 0x0002;
         }
      }
      if (max(jj*levmult-1, 0)-jminsize >= 0){ // Test for cell to bottom
         int nb  = hash[((max(jj*levmult-1, 0) )-jminsize)*(imaxsize-iminsize)+((ii*levmult)-iminsize)];
         if ( nb  >= (int)(ncells+noffset) && (border_cell_needed[nb -ncells-noffset] & 0x0004) == 0x0004) {
            iborder = 0x0004;
         } else if (min((ii+1)*levmult -1,imaxcalc-1) < imaxsize) {
            int nbr = hash[((max(jj*levmult-1, 0) )-jminsize)*(imaxsize-iminsize)+(min((ii+1)*levmult -1,imaxcalc-1)-iminsize)];
            if ( nbr >= (int)(ncells+noffset) && (border_cell_needed[nbr-ncells-noffset] & 0x0004) == 0x0004) iborder = 0x0004;
         }
      }
      if (min( (jj+1)*levmult, jmaxcalc-1) < jmaxsize){ // Test for cell to top
         int nt  = hash[((min( (jj+1)*levmult, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((ii*levmult)-iminsize)];
         if ( nt  >= (int)(ncells+noffset) && (border_cell_needed[nt -ncells-noffset] & 0x0008) == 0x0008) {
            iborder = 0x0008;
         } else if (min((ii+1)*levmult -1,imaxcalc-1) < imaxsize) {
            int ntr = hash[((min( (jj+1)*levmult, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+(min((ii+1)*levmult -1,imaxcalc)-iminsize)];
            if ( ntr >= (int)(ncells+noffset) && (border_cell_needed[ntr-ncells-noffset] & 0x0008) == 0x0008) iborder = 0x0008;
         }
      }
      if (iborder) iborder |= 0x0016;
   }

// Setting up for scan of border cells
   if (iborder) itile[tiX] = 1; 

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset = ntX >> 1; offset > 32; offset >>= 1) { 
      if (tiX < offset) {
         itile[tiX] += itile[tiX+offset]; 
      }    
      barrier(CLK_LOCAL_MEM_FENCE);
   }    

   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  itile[tiX] += itile[tiX+32];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+1]; }

   if (tiX == 0) { 
     ioffset[group_id] = itile[0];
     (*nbsize_new) = itile[0];
   }    

   border_cell_needed_out[giX] = iborder;
}

__kernel void get_border_data2_cl(
                          const int   isize,                  // 0
                 __global       uint  *ioffset,               // 1
                 __global const int   *border_cell,           // 2
                 __global const int   *border_cell_i,         // 3
                 __global const int   *border_cell_j,         // 4
                 __global const int   *border_cell_level,     // 5
                 __global const int   *border_cell_num,       // 6
                 __global       int   *border_cell_i_new,     // 7
                 __global       int   *border_cell_j_new,     // 8
                 __global       int   *border_cell_level_new, // 9
                 __global       int   *indices_needed,        // 10
                 __local        uint  *itile)                 // 11
{
   const uint giX = get_global_id(0);
   const uint tiX = get_local_id(0);
   const uint group_id = get_group_id(0);

   const uint lane   = tiX & 31; 
   const uint warpid = tiX >> 5;

   // Step 1: load global data into tile
   int temp_val = 0;
   if (giX < isize) temp_val = border_cell[giX];
   itile[tiX] = 0;
   if (temp_val > 0) itile[tiX] = 1;
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 2: scan each warp
   uint val = scan_warp_exclusive(itile, tiX, lane);
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 3: Collect per-warp sums
   if (lane == 31) itile[warpid] = itile[tiX];
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 4: Use 1st warp to scan per-warp sums
   if (warpid == 0) scan_warp_inclusive(itile, tiX, lane);
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 5: Accumulate results from Steps 2 and 4
   if (warpid > 0) val += itile[warpid-1];
   barrier(CLK_LOCAL_MEM_FENCE);

   if (giX >= isize || temp_val <= 0) return;

   // Step 6: Write and return the final result
   //itile[tiX] = val;
   //barrier(CLK_LOCAL_MEM_FENCE);

   val += ioffset[group_id];   //index to write to for each thread

   border_cell_i_new[val]     = border_cell_i[giX];
   border_cell_j_new[val]     = border_cell_j[giX];
   border_cell_level_new[val] = border_cell_level[giX];
   indices_needed[val]        = border_cell_num[giX];
}


__kernel void calc_layer2_sethash_cl (
                          const int  isize,               // 0 
                          const int  ncells,              // 1 
                          const int  noffset,             // 2 
                          const int  levmx,               // 3 
                          const int  imax,                // 4 
                          const int  jmax,                // 5 
                 __global const int4 *sizes,              // 6 
                 __global const int  *levtable,           // 7 
                 __global const int  *lev_ibeg,           // 8
                 __global const int  *lev_iend,           // 9
                 __global const int  *lev_jbeg,           // 10
                 __global const int  *lev_jend,           // 11
                 __global const int  *border_cell_i,      // 12
                 __global const int  *border_cell_j,      // 13
                 __global const int  *border_cell_level,  // 14
                 __global const int  *border_cell_num,    // 15
                 __global       int  *border_cell_needed, // 16
                 __global       int  *hash)               // 17
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int ii = border_cell_i[giX];
   int jj = border_cell_j[giX];
   int lev = border_cell_level[giX];
   int cell_number = border_cell_num[giX];
   int levmult = levtable[levmx-lev];

#if HASH_SETUP_OPT_LEVEL == 0
   /* Original Hash Setup */

   int iimin =  ii   *levmult-iminsize;
   int iimax = (ii+1)*levmult-iminsize;
   int jjmin =  jj   *levmult-jminsize;
   int jjmax = (jj+1)*levmult-jminsize;

   if      (ii < lev_ibeg[lev]) iimin = 0;                                 // left boundary
   else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]-iminsize; // right boundary
   else if (jj < lev_jbeg[lev]) jjmin = 0;                                 // bottom boundary
   else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]-jminsize; // top boundary

   for (   int jjj = jjmin; jjj < jjmax; jjj++) {
      for (int iii = iimin; iii < iimax; iii++) {
         hashval_local(jjj, iii) = -(ncells+giX);
      }
   }

#elif HASH_SETUP_OPT_LEVEL >= 1
   /* Just the outer cells */

   int iimin =  ii   *levmult-iminsize;
   int iimax = (ii+1)*levmult-iminsize;
   int jjmin =  jj   *levmult-jminsize;
   int jjmax = (jj+1)*levmult-jminsize;

   if      (ii < lev_ibeg[lev]) iimin = 0;                                 // left boundary
   else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]-iminsize; // right boundary
   else if (jj < lev_jbeg[lev]) jjmin = 0;                                 // bottom boundary
   else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]-jminsize; // top boundary

   for (int iii = iimin; iii < iimax; iii++) {
      hashval_local(jjmin,   iii) = -(ncells+giX);
      hashval_local(jjmax-1, iii) = -(ncells+giX);
   }
   for (int jjj = jjmin+1; jjj < jjmax-1; jjj++) {
      hashval_local(jjj, iimin) = -(ncells+giX);
      hashval_local(jjj, iimax-1) = -(ncells+giX);
   }
#endif


#ifdef XXX
   if (ii < lev_ibeg[lev]) { // left boundary
      for (int    jjj = jj*levmult-jminsize; jjj < (jj+1)*levmult-jminsize; jjj++) {
         for (int iii = 0;                   iii < (ii+1)*levmult-iminsize; iii++) {
            hash[jjj*(imaxsize-iminsize)+iii] = cell_number;
         }
      }
   } else if (ii > lev_iend[lev]) { // right boundary
      for (int    jjj = jj*levmult-jminsize; jjj < (jj+1)*levmult-jminsize;           jjj++) {
         for (int iii = ii*levmult-iminsize; iii < (imax+1)*levtable[levmx]-iminsize; iii++) {
            hash[jjj*(imaxsize-iminsize)+iii] = cell_number;
         }
      }
   } else if (jj < lev_jbeg[lev]) { // bottom boundary
      for (int    jjj = 0;                   jjj < (jj+1)*levmult-jminsize; jjj++) {
         for (int iii = ii*levmult-iminsize; iii < (ii+1)*levmult-iminsize; iii++) {
            hash[jjj*(imaxsize-iminsize)+iii] = cell_number;
         }
      }
   } else if (jj > lev_jend[lev]) { // top boundary
      for (int    jjj = jj*levmult-jminsize; jjj < (jmax+1)*levtable[levmx]-jminsize; jjj++) {
         for (int iii = ii*levmult-iminsize; iii < (ii+1)*levmult-iminsize;           iii++) {
            hash[jjj*(imaxsize-iminsize)+iii] = cell_number;
         }
      }
   } else if (lev == levmx) {
      hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)] = cell_number;
   } else {
      for (int    j = max(jj*levmult-jminsize,0); j < min((jj+1)*levmult,jmaxsize)-jminsize; j++) {
         for (int i = max(ii*levmult-iminsize,0); i < min((ii+1)*levmult,imaxsize)-iminsize; i++) {
            hash[j*(imaxsize-iminsize)+i] = cell_number;
         }
      }
   }
#endif

}

__kernel void copy_mesh_data_cl(
                          const int  isize,         // 0
                 __global const int  *celltype_old, // 1
                 __global       int  *celltype,     // 2
                 __global const int  *i_old,        // 3
                 __global       int  *i,            // 4
                 __global const int  *j_old,        // 5
                 __global       int  *j,            // 6
                 __global const int  *level_old,    // 7
                 __global       int  *level,        // 8
                 __global const int  *nlft_old,     // 9
                 __global       int  *nlft,         // 10
                 __global const int  *nrht_old,     // 11
                 __global       int  *nrht,         // 12
                 __global const int  *nbot_old,     // 13
                 __global       int  *nbot,         // 14
                 __global const int  *ntop_old,     // 15
                 __global       int  *ntop)         // 16
{

   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   celltype[giX] = celltype_old[giX];
   i[giX]        = i_old[giX];
   j[giX]        = j_old[giX];
   level[giX]    = level_old[giX];
   nlft[giX]     = nlft_old[giX];
   nrht[giX]     = nrht_old[giX];
   nbot[giX]     = nbot_old[giX];
   ntop[giX]     = ntop_old[giX];
}

__kernel void fill_mesh_ghost_cl (
                          const int  isize,              // 0 
                          const int  ncells,             // 1 
                 __global const int  *lev_ibegin,        // 2
                 __global const int  *lev_iend,          // 3
                 __global const int  *lev_jbegin,        // 4
                 __global const int  *lev_jend,          // 5
                 __global const int  *border_cell_i,     // 6
                 __global const int  *border_cell_j,     // 7
                 __global const int  *border_cell_level, // 8
                 __global       int  *i,                 // 9 
                 __global       int  *j,                 // 10
                 __global       int  *level,             // 11
                 __global       int  *celltype,          // 12
                 __global       int  *nlft,              // 13
                 __global       int  *nrht,              // 14
                 __global       int  *nbot,              // 15
                 __global       int  *ntop)              // 16
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   int ncout = ncells+giX;

   int ii  = border_cell_i[giX];
   int jj  = border_cell_j[giX];
   int lev = border_cell_level[giX];

   if (ii < lev_ibegin[lev]) celltype[ncout] = LEFT_BOUNDARY;
   if (ii > lev_iend[lev])   celltype[ncout] = RIGHT_BOUNDARY;
   if (jj < lev_jbegin[lev]) celltype[ncout] = BOTTOM_BOUNDARY;
   if (jj > lev_jend[lev])   celltype[ncout] = TOP_BOUNDARY;

   i[ncout]     = ii;
   j[ncout]     = jj;
   level[ncout] = lev;
   nlft[ncout]  = -1;
   nrht[ncout]  = -1;
   nbot[ncout]  = -1;
   ntop[ncout]  = -1;
}

__kernel void fill_neighbor_ghost_cl (
                          const int  ncells_ghost, // 0 
                          const int  levmx,        // 1 
                          const int  imax,         // 2 
                          const int  jmax,         // 3 
                 __global const int4 *sizes,       // 4 
                 __global const int  *levtable,    // 5 
                 __global const int  *i,           // 6 
                 __global const int  *j,           // 7 
                 __global const int  *level,       // 8 
                 __global const int  *hash,        // 9 
                 __global       int  *nlft,        // 10
                 __global       int  *nrht,        // 11
                 __global       int  *nbot,        // 12
                 __global       int  *ntop)        // 13
{
   const uint giX  = get_global_id(0);

   if (giX >= ncells_ghost) return;

   int nlftval = nlft[giX];
   int nrhtval = nrht[giX];
   int nbotval = nbot[giX];
   int ntopval = ntop[giX];

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int imaxcalc = (imax+1)*levtable[levmx];
   int jmaxcalc = (jmax+1)*levtable[levmx];

   int ii = i[giX];
   int jj = j[giX];
   int lev = level[giX];
   int levmult = levtable[levmx-lev];

   if (nlftval == -1){
      int iii = max(ii*levmult-1, 0)-iminsize;
      int jjj = (jj*levmult)-jminsize;
      if (iii >= 0 && iii < imaxsize-iminsize && jjj >= 0 && jjj < jmaxsize-jminsize){  // Test for cell to left
         nlft[giX] = hash[jjj*(imaxsize-iminsize)+iii];
      }
   }
   if (nrhtval == -1){
      int iii = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
      int jjj = (jj*levmult)-jminsize;
      if (iii >= 0 && iii < imaxsize-iminsize && jjj >= 0 && jjj < jmaxsize-jminsize){ // Test for cell to right
         nrht[giX] = hash[jjj*(imaxsize-iminsize)+iii];
      }
   }
   if (nbotval == -1){
      int iii = (ii*levmult)-iminsize;
      int jjj = max(jj*levmult-1, 0)-jminsize;
      if (iii >= 0 && iii < imaxsize-iminsize && jjj >= 0 && jjj < jmaxsize-jminsize){ // Test for cell to bottom
         nbot[giX] = hash[jjj*(imaxsize-iminsize)+iii];
      }
   }
   if (ntopval == -1) {
      int iii = (ii*levmult)-iminsize;
      int jjj = min((jj+1)*levmult,jmaxcalc-1)-jminsize;
      if (iii >= 0 && iii < imaxsize-iminsize && jjj >= 0 && jjj < jmaxsize-jminsize){ // Test for cell to top
         ntop[giX] = hash[jjj*(imaxsize-iminsize)+iii];
      }
   }
}

__kernel void set_corner_neighbor_cl(
                          const int  nghost,    // 0 
                          const int  ncells,    // 1 
                          const int  levmx,     // 2 
                 __global const int4 *sizes,    // 3 
                 __global const int  *levtable, // 4 
                 __global const int  *i,        // 5 
                 __global const int  *j,        // 6 
                 __global const int  *level,    // 7 
                 __global const int  *hash,     // 8 
                 __global       int  *nlft,     // 9 
                 __global       int  *nrht,     // 10
                 __global       int  *nbot,     // 11
                 __global       int  *ntop)     // 12
{
   uint giX  = get_global_id(0);

   if (giX >= nghost) return;

   giX += ncells;

   int nlftval = nlft[giX];
   int nrhtval = nrht[giX];
   int nbotval = nbot[giX];
   int ntopval = ntop[giX];
  
   if (nlftval == INT_MIN){
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;
      int lev = level[giX];
      int levmult = levtable[levmx-lev];
      int jj = j[giX]*levmult;
      int ii = i[giX]*levmult;
      nlft[giX] = hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)];
   }
   if (nrhtval == INT_MIN){
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;
      int lev = level[giX];
      int levmult = levtable[levmx-lev];
      int jj = j[giX]*levmult;
      int ii = i[giX]*levmult;
      nrht[giX] = hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)];
   }
   if (nbotval == INT_MIN) {
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;
      int lev = level[giX];
      int levmult = levtable[levmx-lev];
      int jj = j[giX]*levmult;
      int ii = i[giX]*levmult;
      nbot[giX] = hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)];
   }
   if (ntopval == INT_MIN) {
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;
      int lev = level[giX];
      int levmult = levtable[levmx-lev];
      int jj = j[giX]*levmult;
      int ii = i[giX]*levmult;
      ntop[giX] = hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)];
   }
}

__kernel void adjust_neighbors_local_cl(
                          const int  ncells_ghost,    // 0 
                          const int  ncells,          // 1 
                          const int  noffset,         // 2 
                 __global const int  *indices_needed, // 3 
                 __global       int  *nlft,           // 4 
                 __global       int  *nrht,           // 5
                 __global       int  *nbot,           // 6
                 __global       int  *ntop)           // 7
{
   const uint giX  = get_global_id(0);

   if (giX >= ncells_ghost) return;

   int nghost = ncells_ghost - ncells;
   int nlftval = nlft[giX];
   int nrhtval = nrht[giX];
   int nbotval = nbot[giX];
   int ntopval = ntop[giX];

   if (nlftval <= -ncells && nlftval > -ncells_ghost){
      nlftval = abs(nlftval);
   } else if (nlftval >= noffset && nlftval < noffset+ncells) {
      nlftval -= noffset;
   }
   if (nrhtval <= -ncells && nrhtval > -ncells_ghost){
      nrhtval = abs(nrhtval);
   } else if (nrhtval >= noffset && nrhtval < noffset+ncells) {
      nrhtval -= noffset;
   }
   if (nbotval <= -ncells && nbotval > -ncells_ghost){
      nbotval = abs(nbotval);
   } else if (nbotval >= noffset && nbotval < noffset+ncells) {
      nbotval -= noffset;
   }
   if (ntopval <= -ncells && ntopval > -ncells_ghost){
      ntopval = abs(ntopval);
   } else if (ntopval >= noffset && ntopval < noffset+ncells) {
      ntopval -= noffset;
   }

   nlft[giX] = nlftval;
   nrht[giX] = nrhtval;
   nbot[giX] = nbotval;
   ntop[giX] = ntopval;
}

__kernel void finish_reduction_scan_cl(
                 const int   isize,
        __global       int  *ioffset,
        __global       int  *result,
        __local        int  *itile_scratch,
        __local        uint *itile)
{
   const uint tiX = get_local_id(0);
   const uint gID = get_group_id(0);
   const uint ntX = get_local_size(0);

   const uint lane = tiX & 31;
   const uint warpID = tiX >> 5;
   const uint EPT = (isize+ntX-1)/ntX; //elements_per_thread;
   const uint BLOCK_SIZE = EPT * ntX;

   uint reduceValue = 0;

// #pragma unroll 4
   for(uint i = 0; i < EPT; ++i)
   {
      uint offsetIdx = i * ntX + (gID * BLOCK_SIZE) + tiX;

#ifdef IS_NVIDIA
//    if (offsetIdx >= isize) return;
#endif

      // Step 1: Read ntX elements from global (off-chip) memory to local memory (on-chip)
      uint input = 0;
      if (offsetIdx < isize) input = ioffset[offsetIdx];           
      itile[tiX] = input;           
      barrier(CLK_LOCAL_MEM_FENCE);

      // Step 2: Perform scan on ntX elements
      uint val = scan_workgroup_exclusive(itile, tiX, lane, warpID);

      // Step 3: Propagate reduced result from previous block of ntX elements
      val += reduceValue;

      // Step 4: Write out data to global memory
      if (offsetIdx < isize) ioffset[offsetIdx] = val;

      // Step 5: Choose reduced value for next iteration
      if (tiX == (ntX-1)) {
         itile[tiX] = input + val;
         if (i == EPT-1) result[0] = input + val;
      }
      barrier(CLK_LOCAL_MEM_FENCE);

      reduceValue = itile[ntX-1];
      barrier(CLK_LOCAL_MEM_FENCE);
   }

}

void reduction_minmax_within_tile4(__local  int4  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

    for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        if (tile[tiX+offset].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+offset].s0;
        if (tile[tiX+offset].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+offset].s1;
        if (tile[tiX+offset].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+offset].s2;
        if (tile[tiX+offset].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+offset].s3;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tiX < 32){
      if (tile[tiX+32].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+32].s0;
      if (tile[tiX+32].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+32].s1;
      if (tile[tiX+32].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+32].s2;
      if (tile[tiX+32].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+32].s3;
      barrier(CLK_LOCAL_MEM_FENCE);

      if (tile[tiX+16].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+16].s0;
      if (tile[tiX+16].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+16].s1;
      if (tile[tiX+16].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+16].s2;
      if (tile[tiX+16].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+16].s3;
      barrier(CLK_LOCAL_MEM_FENCE);

      if (tile[tiX+8].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+8].s0;
      if (tile[tiX+8].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+8].s1;
      if (tile[tiX+8].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+8].s2;
      if (tile[tiX+8].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+8].s3;
      barrier(CLK_LOCAL_MEM_FENCE);

      if (tile[tiX+4].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+4].s0;
      if (tile[tiX+4].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+4].s1;
      if (tile[tiX+4].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+4].s2;
      if (tile[tiX+4].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+4].s3;
      barrier(CLK_LOCAL_MEM_FENCE);

      if (tile[tiX+2].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+2].s0;
      if (tile[tiX+2].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+2].s1;
      if (tile[tiX+2].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+2].s2;
      if (tile[tiX+2].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+2].s3;
      barrier(CLK_LOCAL_MEM_FENCE);

      if (tile[tiX+1].s0 < tile[tiX].s0) tile[tiX].s0 = tile[tiX+1].s0;
      if (tile[tiX+1].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+1].s1;
      if (tile[tiX+1].s2 < tile[tiX].s2) tile[tiX].s2 = tile[tiX+1].s2;
      if (tile[tiX+1].s3 > tile[tiX].s3) tile[tiX].s3 = tile[tiX+1].s3;
    }
}

__kernel void calc_hash_size_cl(
                          const int   ncells,      // 0 
                          const int   levmx,       // 1 
                 __global const int   *levtable,   // 2 
                 __global const int   *level,      // 3 
                 __global const int   *i,          // 4 
                 __global const int   *j,          // 5 
                 __global       int4  *redscratch, // 6
                 __global       int4  *sizes,      // 7
                 __local        int4  *tile)       // 8
{
   const unsigned int giX  = get_global_id(0);
   const unsigned int tiX  = get_local_id(0);

   tile[tiX].s0 = 2147483647;
   tile[tiX].s1 = 0;
   tile[tiX].s2 = 2147483647;
   tile[tiX].s3 = 0;

   if (giX >= ncells) return;

   int group_id = get_group_id(0);

   int lev = level[giX];

   tile[tiX].s0 = i[giX]   *levtable[levmx-lev]; // imincalc
   tile[tiX].s1 =(i[giX]+1)*levtable[levmx-lev]; // imaxcalc

   tile[tiX].s2 = j[giX]   *levtable[levmx-lev]; // jmincalc
   tile[tiX].s3 =(j[giX]+1)*levtable[levmx-lev]; // jmaxcalc

   barrier(CLK_LOCAL_MEM_FENCE);
   reduction_minmax_within_tile4(tile);

   //  Write the local value back to an array size of the number of groups
   if (tiX == 0){
      redscratch[group_id].s0123 = tile[0].s0123;
      if (group_id == 0) sizes[0].s0123 = tile[0].s0123;
   }
}

__kernel void finish_reduction_minmax4_cl(
                 const int    isize,        // 0
        __global       int4   *redscratch,  // 1
        __global       int4   *sizes,       // 2
        __local        int4   *tile)        // 3
{  
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX].s0 = 2147483647;
   tile[tiX].s1 = 0;
   tile[tiX].s2 = 2147483647;
   tile[tiX].s3 = 0;

   if (tiX < isize) tile[tiX].s0123 = redscratch[giX].s0123;

   for (giX += ntX; giX < isize; giX += ntX) {
     if (redscratch[giX].s0 < tile[tiX].s0) tile[tiX].s0 = redscratch[giX].s0;
     if (redscratch[giX].s1 > tile[tiX].s1) tile[tiX].s1 = redscratch[giX].s1;
     if (redscratch[giX].s2 < tile[tiX].s2) tile[tiX].s2 = redscratch[giX].s2;
     if (redscratch[giX].s3 > tile[tiX].s3) tile[tiX].s3 = redscratch[giX].s3;
     barrier(CLK_LOCAL_MEM_FENCE);
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_minmax_within_tile4(tile);

   if (tiX == 0) {
     redscratch[0].s0123 = tile[0].s0123;
     sizes[0].s0123 = tile[0].s0123;
   }
}

__kernel void calc_spatial_coordinates_cl(
          const    int   isize,
          const    real   xmin,
          const    real   ymin,
          __global real   *lev_deltax,
          __global real   *lev_deltay,
          __global real   *x,
          __global real   *dx,
          __global real   *y,
          __global real   *dy,
          __global int    *level,
          __global int    *i,
          __global int    *j)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int lev = level[giX];

   x[giX] = isize;
   x[giX]  = xmin + lev_deltax[lev] * (real)i[giX];
   dx[giX] =        lev_deltax[lev];
   y[giX]  = ymin + lev_deltay[lev] * (real)j[giX];
   dy[giX] =        lev_deltay[lev];
}

__kernel void do_load_balance_cl(
                  const int   ncells,
                  const int   lower_block_size,
                  const int   middle_block_size,
                  const int   middle_block_start,
         __global const real *state_var_old,
         __global const real *state_var_lower,
         __global const real *state_var_upper,
         __global       real *state_var_new)
{
   const uint giX = get_global_id(0);

   //if(giX >= ncells) return;

   //state_var_new[giX] = 0.0;

   if(giX < lower_block_size) {
      state_var_new[giX] = state_var_lower[giX];
   } else if(giX < lower_block_size + middle_block_size) {
      //const uint destgiX = giX;
      const uint srcgiX = middle_block_start + giX - lower_block_size;

      state_var_new[giX] = state_var_old[srcgiX];
   } else if(giX < ncells) {
      const uint srcgiX = giX - (lower_block_size + middle_block_size);

      state_var_new[giX] = state_var_upper[srcgiX];
   }
}

__kernel void do_load_balance_lower_cl(
         __global       int  *dev_i_new,
         __global       int  *dev_j_new,
         __global       int  *dev_level_new,
         __global       int  *dev_celltype_new,
         __global const int  *dev_i_lower,
         __global const int  *dev_j_lower,
         __global const int  *dev_level_lower,
         __global const int  *dev_celltype_lower,
         const int            lower_block_size)
{

   const unsigned int giX = get_global_id(0);

   if(giX >= lower_block_size) return;

   dev_i_new[giX]        = dev_i_lower[giX];
   dev_j_new[giX]        = dev_j_lower[giX];
   dev_level_new[giX]    = dev_level_lower[giX];
   dev_celltype_new[giX] = dev_celltype_lower[giX];

}


__kernel void do_load_balance_middle_cl(
         __global       int  *dev_i_new,
         __global       int  *dev_j_new,
         __global       int  *dev_level_new,
         __global       int  *dev_celltype_new,
         __global const int  *dev_i,
         __global const int  *dev_j,
         __global const int  *dev_level,
         __global const int  *dev_celltype,
         const int            lower_block_size,
         const int            middle_block_size,
         const int            middle_block_start)
{

   const unsigned int giX = get_global_id(0);

   if(giX >= middle_block_size) return;

   const unsigned int rgiX = lower_block_size + giX;
   const unsigned int dgiX = middle_block_start + giX;

   dev_i_new[rgiX]        = dev_i[dgiX];
   dev_j_new[rgiX]        = dev_j[dgiX];
   dev_level_new[rgiX]    = dev_level[dgiX];
   dev_celltype_new[rgiX] = dev_celltype[dgiX];

}


__kernel void do_load_balance_upper_cl(
         __global       int  *dev_i_new,
         __global       int  *dev_j_new,
         __global       int  *dev_level_new,
         __global       int  *dev_celltype_new,
         __global const int  *dev_i_upper,
         __global const int  *dev_j_upper,
         __global const int  *dev_level_upper,
         __global const int  *dev_celltype_upper,
         const int            lower_block_size,
         const int            middle_block_size,
         const int            upper_block_size)
{

   const unsigned int giX = get_global_id(0);

   if(giX >= upper_block_size) return;

   const unsigned int rgiX = lower_block_size + middle_block_size + giX;

   dev_i_new[rgiX]        = dev_i_upper[giX];
   dev_j_new[rgiX]        = dev_j_upper[giX];
   dev_level_new[rgiX]    = dev_level_upper[giX];
   dev_celltype_new[rgiX] = dev_celltype_upper[giX];

}

void reduction_sum_int_within_tile(__local  int  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

    for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        tile[tiX] += tile[tiX+offset];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tiX < 32){
      tile[tiX] += tile[tiX+32];
      tile[tiX] += tile[tiX+16];
      tile[tiX] += tile[tiX+8];
      tile[tiX] += tile[tiX+4];
      tile[tiX] += tile[tiX+2];
      tile[tiX] += tile[tiX+1];
    }

}

__kernel void count_BCs_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const int   *i,         // 1
        __global const int   *j,         // 2
        __global const int   *level,     // 3
        __global const int   *lev_ibeg,  // 4
        __global const int   *lev_iend,  // 5
        __global const int   *lev_jbeg,  // 6
        __global const int   *lev_jend,  // 7
        __global       int   *scratch,   // 8
        __local        int   *tile)      // 9
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);

    const unsigned int group_id = get_group_id(0);

    tile[tiX] = 0;
    if (giX < isize) {
      if (i[giX] == lev_ibeg[level[giX]]) tile[tiX]++;
      if (i[giX] == lev_iend[level[giX]]) tile[tiX]++;
      if (j[giX] == lev_jbeg[level[giX]]) tile[tiX]++;
      if (j[giX] == lev_jend[level[giX]]) tile[tiX]++;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_sum_int_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      scratch[group_id] = tile[0];
    }
}

