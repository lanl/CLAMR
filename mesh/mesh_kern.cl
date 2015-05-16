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

#if !defined(FULL_PRECISION) && !defined(MIXED_PRECISION) && !defined(MINIMUM_PRECISION)
#define FULL_PRECISION
#endif
#ifdef NO_CL_DOUBLE
#undef  FULL_PRECISION
#undef  MIXED_PRECISION
#define MINIMUM_PRECISION
#endif

#if defined(MINIMUM_PRECISION)
typedef float  real_t;
typedef float  spatial_t;
#define HALF 0.5f
#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double  real_t;
typedef float  spatial_t;
#define HALF 0.5
#elif defined(FULL_PRECISION)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double real_t;
typedef double spatial_t;
#define HALF 0.5
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

__kernel void hash_adjust_sizes_cl(
                          const int  isize,      // 0
                          const int  levmx,      // 1
                          const int  imax,       // 2
                          const int  jmax,       // 3
                 __global const int  *levtable,  // 4
                 __global const int4 *sizes)     // 5
{
   const uint giX  = get_global_id(0);

   if (giX >= isize) return;

   // Expand size by 2*coarse_cells
   sizes[0].s0 = max(sizes[0].s0-2*levtable[levmx],0);
   sizes[0].s1 = min(sizes[0].s1+2*levtable[levmx],(imax+1)*levtable[levmx]);
   sizes[0].s2 = max(sizes[0].s2-2*levtable[levmx],0);
   sizes[0].s3 = min(sizes[0].s3+2*levtable[levmx],(jmax+1)*levtable[levmx]);
}

__kernel void neighbor_init_cl(
                          const int   isize,            // 0
                 __global       int   *nlft,            // 1
                 __global       int   *nrht,            // 2
                 __global       int   *nbot,            // 3
                 __global       int   *ntop)            // 4
{
   const uint ic  = get_global_id(0);

   if (ic >= isize) return;

   nlft[ic] = -1;
   nrht[ic] = -1;
   nbot[ic] = -1;
   ntop[ic] = -1;
}
                          
__kernel void hash_setup_cl(
                          const int   isize,            // 0
                          const int   levmx,            // 1
                          const int   imaxsize,         // 2
                 __global const int   *levtable,        // 3
                 __global const int   *level,           // 4
                 __global const int   *i,               // 5
                 __global const int   *j,               // 6
                 __global const int   *nlft,            // 7
                 __global const int   *nrht,            // 8
                 __global const int   *nbot,            // 9
                 __global const int   *ntop,            // 10
                 __global const ulong *hash_header,     // 11
                 __global       int   *hash)            // 12
{

   const uint ic  = get_global_id(0);

   if (ic >= isize) return;

   int lev = level[ic];

   bool need_hash = (nlft[ic] == -1 || nrht[ic] == -1 || nbot[ic] == -1 || ntop[ic] == -1) ? true : false;

   if (! need_hash) { 
      if ( (level[nlft[ic]] > lev && ntop[nlft[ic]] == -1) ||
           (level[nrht[ic]] > lev && ntop[nrht[ic]] == -1) ||
           (level[nbot[ic]] > lev && nrht[nbot[ic]] == -1) ||
           (level[ntop[ic]] > lev && nrht[ntop[ic]] == -1) ) need_hash = true;
   }

   if (need_hash) {
      const int hash_method       = (int)hash_header[0];
      const ulong hash_table_size =      hash_header[1];
      const ulong AA              =      hash_header[2];
      const ulong BB              =      hash_header[3];

      int levmult = levtable[levmx-lev];
      int ii = i[ic]*levmult;
      int jj = j[ic]*levmult;

      write_hash(hash_method, hash_table_size, AA, BB, ic, jj*imaxsize+ii, hash);
   }
}

__kernel void hash_setup_local_cl(
                          const int   isize,           // 0
                          const int   levmx,           // 1
                          const int   imax,            // 2
                          const int   jmax,            // 3
                          const int   noffset,         // 4
                 __global const int4  *sizes,          // 5
                 __global const int   *levtable,       // 6
                 __global const int   *level,          // 7
                 __global const int   *i,              // 8
                 __global const int   *j,              // 9
                 __global const ulong *hash_header,    // 10
                 __global       int   *hash)           // 11
{

   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

/*
   // Expand size by 2*coarse_cells
   sizes[0].s0 = max(sizes[0].s0-2*levtable[levmx],0);
   sizes[0].s1 = min(sizes[0].s1+2*levtable[levmx],(imax+1)*levtable[levmx]);
   sizes[0].s2 = max(sizes[0].s2-2*levtable[levmx],0);
   sizes[0].s3 = min(sizes[0].s3+2*levtable[levmx],(jmax+1)*levtable[levmx]);
*/

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int lev = level[giX];
   int levmult = levtable[levmx-lev];
   int cell_number = giX+noffset;
   int ii = i[giX]*levmult-iminsize;
   int jj = j[giX]*levmult-jminsize;

   write_hash(hash_method, hash_table_size, AA, BB, cell_number, jj*(imaxsize-iminsize)+ii, hash);
}

__kernel void calc_neighbors_cl(
                          const int   isize,            // 0 
                          const int   levmx,            // 1 
                          const int   imax,             // 2 
                          const int   jmax,             // 3 
                          const int   imaxsize,         // 4 
                          const int   jmaxsize,         // 5 
                 __global const int   *levtable,        // 6
                 __global const int   *level,           // 7
                 __global const int   *i,               // 8
                 __global const int   *j,               // 9
                 __global       int   *nlft,            // 10
                 __global       int   *nrht,            // 11
                 __global       int   *nbot,            // 12
                 __global       int   *ntop,            // 13
                 __global const ulong *hash_header,     // 14
                 __global const int   *hash)            // 15
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

   int ii, jj, iii, jjj, lev, levmult;

   ii = i[giX];
   jj = j[giX];
   lev = level[giX];
   levmult = levtable[levmx-lev];

   int iicur = ii*levmult;
   int iilft = max( (ii-1)*levmult, 0         );   
   int iirht = min( (ii+1)*levmult, imaxsize-1);
   int jjcur = jj*levmult;
   int jjbot = max( (jj-1)*levmult, 0         );   
   int jjtop = min( (jj+1)*levmult, jmaxsize-1);

   int nlftval = nlft[giX];
   int nrhtval = nrht[giX];
   int nbotval = nbot[giX];
   int ntopval = ntop[giX];

   // Taking care of boundary cells
   // Force each boundary cell to point to itself on its boundary direction
   if (nlftval < 0 && iicur <    1*levtable[levmx]  ) nlftval = giX;
   if (nbotval < 0 && jjcur <    1*levtable[levmx]  ) nbotval = giX;
   if (nrhtval < 0 && iicur > imax*levtable[levmx]-1) nrhtval = giX;
   if (ntopval < 0 && jjcur > jmax*levtable[levmx]-1) ntopval = giX;
   // Boundary cells next to corner boundary need special checks
   if (nlftval < 0 && iicur ==    1*levtable[levmx] &&  (jjcur < 1*levtable[levmx] || jjcur >= jmax*levtable[levmx] ) ) nlftval = giX;
   if (nbotval < 0 && jjcur ==    1*levtable[levmx] &&  (iicur < 1*levtable[levmx] || iicur >= imax*levtable[levmx] ) ) nbotval = giX;
   if (nrhtval < 0 && iirht == imax*levtable[levmx] &&  (jjcur < 1*levtable[levmx] || jjcur >= jmax*levtable[levmx] ) ) nrhtval = giX;
   if (ntopval < 0 && jjtop == jmax*levtable[levmx] &&  (iicur < 1*levtable[levmx] || iicur >= imax*levtable[levmx] ) ) ntopval = giX;

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

   nlft[giX] = nlftval;
   nrht[giX] = nrhtval;
   nbot[giX] = nbotval;
   ntop[giX] = ntopval;
}

__kernel void calc_neighbors_local_cl(
                          const int   ncells,          // 0 
                          const int   levmx,           // 1 
                          const int   imax,            // 2 
                          const int   jmax,            // 3 
                          const int   noffset,         // 4
                 __global const int4  *sizes,          // 5
                 __global const int   *levtable,       // 6
                 __global const int   *level,          // 7
                 __global const int   *i,              // 8
                 __global const int   *j,              // 9
                 __global       int   *nlft,           // 10
                 __global       int   *nrht,           // 11
                 __global       int   *nbot,           // 12
                 __global       int   *ntop,           // 13
                 __global const ulong *hash_header,    // 14
                 __global const int   *hash)           // 15
{
                
   const uint giX  = get_global_id(0);

   if (giX >= ncells) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

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

   int iicur = ii*levmult-iminsize;
   int iilft = max( (ii-1)*levmult, 0         )-iminsize;
   int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
   int jjcur = jj*levmult-jminsize;
   int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
   int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

   int nlftval = -1;
   int nrhtval = -1;
   int nbotval = -1;
   int ntopval = -1;

   // Taking care of boundary cells
   // Force each boundary cell to point to itself on its boundary direction
   if (iicur <    1*levtable[levmx]  -iminsize) nlftval = giX+noffset;
   if (jjcur <    1*levtable[levmx]  -jminsize) nbotval = giX+noffset;
   if (iicur > imax*levtable[levmx]-1-iminsize) nrhtval = giX+noffset;
   if (jjcur > jmax*levtable[levmx]-1-jminsize) ntopval = giX+noffset;
   // Boundary cells next to corner boundary need special checks
   if (iicur ==    1*levtable[levmx]-iminsize &&  (jjcur < 1*levtable[levmx]-jminsize || jjcur >= jmax*levtable[levmx]-jminsize ) ) nlftval = giX+noffset;
   if (jjcur ==    1*levtable[levmx]-jminsize &&  (iicur < 1*levtable[levmx]-iminsize || iicur >= imax*levtable[levmx]-iminsize ) ) nbotval = giX+noffset;
   if (iirht == imax*levtable[levmx]-iminsize &&  (jjcur < 1*levtable[levmx]-jminsize || jjcur >= jmax*levtable[levmx]-jminsize ) ) nrhtval = giX+noffset;
   if (jjtop == jmax*levtable[levmx]-jminsize &&  (iicur < 1*levtable[levmx]-iminsize || iicur >= imax*levtable[levmx]-iminsize ) ) ntopval = giX+noffset;

   // need to check for finer neighbor first
   // Right and top neighbor don't change for finer, so drop through to same size
   // Left and bottom need to be half of same size index for finer test
   if (lev != levmx) {
      int iilftfiner = iicur-(iicur-iilft)/2;
      int jjbotfiner = jjcur-(jjcur-jjbot)/2;
      if (nlftval < 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur     *(imaxsize-iminsize)+iilftfiner, hash);
      if (nbotval < 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iicur,      hash);
   }

   // same size neighbor
   if (nlftval < 0) {
      int nlfttry = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilft, hash);
      if (nlfttry >= 0 && nlfttry < ncells && level[nlfttry] == lev) nlftval = nlfttry;
   }
   if (nrhtval < 0) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iirht, hash);
   if (nbotval < 0) { 
      int nbottry = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iicur, hash);
      if (nbottry >= 0 && nbottry < ncells && level[nbottry] == lev) nbotval = nbottry;
   }
   if (ntopval < 0) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iicur, hash);

   // Now we need to take care of special case where bottom and left boundary need adjustment since
   // expected cell doesn't exist on these boundaries if it is finer than current cell
   if (lev != levmx) {
      if (jjcur < 1*levtable[levmx]) {
         if (nrhtval < 0) {
            int jjtopfiner = (jjcur+jjtop)/2;
            nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iirht, hash);
         }
         if (nlftval < 0) {
            int iilftfiner = iicur-(iicur-iilft)/2;
            int jjtopfiner = (jjcur+jjtop)/2;
            nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
         }
      }

      if (iicur < 1*levtable[levmx]) {
         if (ntopval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
         }
         if (nbotval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            int jjbotfiner = jjcur-(jjcur-jjbot)/2;
            nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
         }
      }
   }

   // coarser neighbor
   if (lev != 0){
      if (nlftval < 0) {
         iilft -= iicur-iilft;
         int jjlft = (jj/2)*2*levmult-jminsize;
         int nlfttry = read_hash(hash_method, hash_table_size, AA, BB, jjlft*(imaxsize-iminsize)+iilft, hash);
         if (nlfttry >= 0 && nlfttry < ncells && level[nlfttry] == lev-1) nlftval = nlfttry;
      }
      if (nrhtval < 0) {
         int jjrht = (jj/2)*2*levmult-jminsize;
         int nrhttry = read_hash(hash_method, hash_table_size, AA, BB, jjrht*(imaxsize-iminsize)+iirht, hash);
         if (nrhttry >= 0 && nrhttry < (int)ncells && level[nrhttry] == lev-1) nrhtval = nrhttry;
      }
      if (nbotval < 0) {
         jjbot -= jjcur-jjbot;
         int iibot = (ii/2)*2*levmult-iminsize;
         int nbottry = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iibot, hash);
         if (nbottry >= 0 && nbottry < ncells && level[nbottry] == lev-1) nbotval = nbottry;
      }
      if (ntopval < 0) {
         int iitop = (ii/2)*2*levmult-iminsize;
         int ntoptry = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iitop, hash);
         if (ntoptry >= 0 && ntoptry < (int)ncells && level[ntoptry] == lev-1) ntopval = ntoptry;
      }
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

inline int2 scan_warp_exclusive2(__local volatile int2 *input, const uint idx, const uint lane) {
    int2 zero = 0;
    if (lane > 0 ) input[idx].s01 += input[idx - 1].s01;
    if (lane > 1 ) input[idx].s01 += input[idx - 2].s01;
    if (lane > 3 ) input[idx].s01 += input[idx - 4].s01;
    if (lane > 7 ) input[idx].s01 += input[idx - 8].s01;
    if (lane > 15) input[idx].s01 += input[idx - 16].s01;

    return (lane > 0) ? input[idx-1].s01 : zero;
}

inline int2 scan_warp_inclusive2(__local volatile int2 *input, const uint idx, const uint lane) {
    if (lane > 0 ) input[idx].s01 += input[idx - 1].s01;
    if (lane > 1 ) input[idx].s01 += input[idx - 2].s01;
    if (lane > 3 ) input[idx].s01 += input[idx - 4].s01;
    if (lane > 7 ) input[idx].s01 += input[idx - 8].s01;
    if (lane > 15) input[idx].s01 += input[idx - 16].s01;

    return input[idx].s01;
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

inline int2 scan_workgroup_exclusive2(
    __local int2* itile,
    const uint tiX,
    const uint lane,
    const uint warpID) {

    // Step 1: scan each warp
    int2 val = scan_warp_exclusive2(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 2: Collect per-warp sums
    if (lane == 31) itile[warpID].s01 = itile[tiX].s01;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 3: Use 1st warp to scan per-warp sums
    if (warpID == 0) scan_warp_inclusive2(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 4: Accumulate results from Steps 1 and 3
    if (warpID > 0) val.s01 += itile[warpID-1].s01;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 6: Write and return the final result
    itile[tiX].s01 = val.s01;
    barrier(CLK_LOCAL_MEM_FENCE);

    return val.s01;
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
                          const int   isize,               // 0 
                          const int   ncells,              // 1 
                          const int   levmx,               // 2 
                          const int   imax,                // 3 
                          const int   jmax,                // 4 
                          const int   noffset,             // 5 
                 __global const int4  *sizes,              // 6
                 __global const int   *levtable,           // 7
                 __global const int   *level,              // 8
                 __global const int   *border_cell_i,      // 9
                 __global const int   *border_cell_j,      // 10
                 __global const int   *border_cell_level,  // 11
                 __global       int   *border_cell_needed, // 12
                 __global const ulong *hash_header,        // 13
                 __global const int   *hash)               // 14
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

   border_cell_needed[giX] = 0;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int imaxcalc = (imax+1)*levtable[levmx];
   int jmaxcalc = (jmax+1)*levtable[levmx];

   int iborder = 0;

   // Layer 1
   int jj = border_cell_j[giX];
   int ii = border_cell_i[giX];
   int lev = border_cell_level[giX];
   int levmult = levtable[levmx-lev];

   int iicur = ii*levmult-iminsize;
   int iilft = max( (ii-1)*levmult, 0         )-iminsize;
   int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
   int jjcur = jj*levmult-jminsize;
   int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
   int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

   // Test for cell to left
   if (iicur-(iicur-iilft)/2 >= 0 && iicur-(iicur-iilft)/2 < imaxsize-iminsize && jjcur >= 0 && (jjcur+jjtop)/2 < jmaxsize-jminsize){
      int nlftval = -1;
      // Check for finer cell left and bottom side
      if (lev != levmx){                                // finer neighbor
         int iilftfiner = iicur-(iicur-iilft)/2;
         nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilftfiner, hash);
         // Also check for finer cell left and top side
         if (nlftval < 0) {
            int jjtopfiner = (jjcur+jjtop)/2;
            nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
         }
      }

      if (nlftval < 0 && iilft >= 0) {  // same size
         int nlfttry = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilft, hash);
         // we have to test for same level or it could be a finer cell one cell away that it is matching
         if (nlfttry-noffset >= 0 && nlfttry-noffset < ncells && level[nlfttry-noffset] == lev) {
            nlftval = nlfttry;
         }
      }

      if (lev != 0 && nlftval < 0 && iilft-(iicur-iilft) >= 0){      // coarser neighbor
         iilft -= iicur-iilft;
         int jjlft = (jj/2)*2*levmult-jminsize;
         int nlfttry = read_hash(hash_method, hash_table_size, AA, BB, jjlft*(imaxsize-iminsize)+iilft, hash);
         // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
         if (nlfttry-noffset >= 0 && nlfttry-noffset < ncells && level[nlfttry-noffset] == lev-1) {
           nlftval = nlfttry;
         }
      }
      if (nlftval >= 0) iborder |= 0x0001;
   }

   // Test for cell to right
   if (iirht < imaxsize-iminsize && iirht >= 0 && jjcur >= 0 && jjtop < jmaxsize-jminsize) {
      int nrhtval = -1;
      // right neighbor -- finer, same size and coarser
      nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iirht, hash);
      // right neighbor -- finer right top test
      if (nrhtval < 0 && lev != levmx){
         int jjtopfiner = (jjcur+jjtop)/2;
         nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iirht, hash);
      }
      if (nrhtval < 0 && lev != 0) { // test for coarser, but not directly above
         int jjrhtcoarser = (jj/2)*2*levmult-jminsize;
         if (jjrhtcoarser != jjcur) {
            int nrhttry = read_hash(hash_method, hash_table_size, AA, BB, jjrhtcoarser*(imaxsize-iminsize)+iirht, hash);
            if (nrhttry-noffset >= 0 && nrhttry-noffset < ncells && level[nrhttry-noffset] == lev-1) {
               nrhtval = nrhttry;
            }
         }
      }
      if (nrhtval > 0)  iborder |= 0x0002;
   }

   // Test for cell to bottom
   if (iicur >= 0 && (iicur+iirht)/2 < imaxsize-iminsize && jjcur-(jjcur-jjbot)/2 >= 0 && jjcur-(jjcur-jjbot)/2 < jmaxsize-jminsize){
      int nbotval = -1;
      // Check for finer cell below and left side
      if (lev != levmx){                                // finer neighbor
         int jjbotfiner = jjcur-(jjcur-jjbot)/2;
         nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iicur, hash);
         // Also check for finer cell below and right side
         if (nbotval < 0) {
            int iirhtfiner = (iicur+iirht)/2;
            nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
         }
      }

      if (nbotval < 0 && jjbot >= 0) {  // same size
         int nbottry = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iicur, hash);
         // we have to test for same level or it could be a finer cell one cell away that it is matching
         if (nbottry-noffset >= 0 && nbottry-noffset < ncells && level[nbottry-noffset] == lev) {
            nbotval = nbottry;
         }
      }

      if (lev != 0 && nbotval < 0 && jjbot-(jjcur-jjbot) >= 0){      // coarser neighbor
         jjbot -= jjcur-jjbot;
         int iibot = (ii/2)*2*levmult-iminsize;
         int nbottry = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iibot, hash);
         // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
         if (nbottry-noffset >= 0 && nbottry-noffset < ncells && level[nbottry-noffset] == lev-1) {
           nbotval = nbottry;
         }
      }
      if (nbotval >= 0) iborder |= 0x0004;
   }

   // Test for cell to top
   if (iirht < imaxsize-iminsize && iicur >= 0 && jjtop >= 0 && jjtop < jmaxsize-jminsize) {
      int ntopval = -1;
      // top neighbor -- finer, same size and coarser
      ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iicur, hash);
      // top neighbor -- finer top right test
      if (ntopval < 0 && lev != levmx){
         int iirhtfiner = (iicur+iirht)/2;
         ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
      }
      if (ntopval < 0 && lev != 0) { // test for coarser, but not directly above
         int iitopcoarser = (ii/2)*2*levmult-iminsize;
         if (iitopcoarser != iicur) {
            int ntoptry = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iitopcoarser, hash);
            if (ntoptry-noffset >= 0 && ntoptry-noffset < ncells && level[ntoptry-noffset] == lev-1) {
               ntopval = ntoptry;
            }
         }
      }
      if (ntopval > 0)  iborder |= 0x0008;
   }

   if (iborder) border_cell_needed[giX] = iborder;
}

__kernel void calc_layer1_sethash_cl (
                          const int   isize,               // 0 
                          const int   ncells,              // 1 
                          const int   noffset,             // 2 
                          const int   levmx,               // 3 
                 __global const int4  *sizes,              // 4 
                 __global const int   *levtable,           // 5 
                 __global const int   *border_cell_i,      // 6
                 __global const int   *border_cell_j,      // 7
                 __global const int   *border_cell_level,  // 8
                 __global       int   *border_cell_needed, // 9
                 __global const ulong *hash_header,        // 10
                 __global       int   *hash)               // 11
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

   int iborder = border_cell_needed[giX];

   if (iborder) {
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;

      int lev = border_cell_level[giX];
      int levmult = levtable[levmx-lev];
      int ii = border_cell_i[giX]*levmult-iminsize;
      int jj = border_cell_j[giX]*levmult-jminsize;

      write_hash(hash_method, hash_table_size, AA, BB, ncells+noffset+giX, jj*(imaxsize-iminsize)+ii, hash);
   }
}

__kernel void calc_layer2_cl (
                          const int   nbsize_local,            // 0 
                          const int   ncells,                  // 1 
                          const int   noffset,                 // 2 
                          const int   levmx,                   // 3 
                          const int   imax,                    // 4 
                          const int   jmax,                    // 5 
                 __global const int4  *sizes,                  // 6
                 __global const int   *levtable,               // 7
                 __global const int   *level,                  // 8
                 __global const int   *border_cell_i,          // 9
                 __global const int   *border_cell_j,          // 10
                 __global const int   *border_cell_level,      // 11
                 __global const int   *border_cell_needed,     // 12
                 __global       int   *border_cell_needed_out, // 13
                 __global const ulong *hash_header,            // 14
                 __global const int   *hash,                   // 15
                 __global       int   *ioffset,                // 16
                 __global       int   *nbsize_new,             // 17
                 __local        int   *itile)                  // 18
{
   const uint giX = get_global_id(0);
   const uint tiX = get_local_id(0);
   const uint group_id = get_group_id(0);
   const uint ntX = get_local_size(0);

   itile[tiX] = 0;

   if (giX >= nbsize_local) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

   int iborder = border_cell_needed[giX];

   if (iborder <= 0) {
      int iminsize = sizes[0].s0;
      int imaxsize = sizes[0].s1;
      int jminsize = sizes[0].s2;
      int jmaxsize = sizes[0].s3;

      int imaxcalc = (imax+1)*levtable[levmx];
      int jmaxcalc = (jmax+1)*levtable[levmx];

      int jj = border_cell_j[giX];
      int ii = border_cell_i[giX];
      int lev = border_cell_level[giX];
      int levmult = levtable[levmx-lev];

      int iicur = ii*levmult-iminsize;
      int iilft = max( (ii-1)*levmult, 0         )-iminsize;
      int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
      int jjcur = jj*levmult-jminsize;
      int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
      int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

      // Test for cell to left
      if (iicur-(iicur-iilft)/2 >= 0 && iicur-(iicur-iilft)/2 < imaxsize-iminsize && jjcur >= 0 &&      (jjcur+jjtop)/2 < jmaxsize-jminsize){
         // Check for finer cell left and bottom side
         if (lev != levmx){                                // finer neighbor
            int iilftfiner = iicur-(iicur-iilft)/2;
            int nl = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilftfiner, hash);
            if (nl >= (int)(ncells+noffset) && (border_cell_needed[nl-ncells-noffset] & 0x0001) == 0x0001) {
               iborder = 0x0001;
            } else {
               // Also check for finer cell left and top side
               int jjtopfiner = (jjcur+jjtop)/2;
               int nlt = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
               if ( nlt >= (int)(ncells+noffset) && (border_cell_needed[nlt-ncells-noffset] & 0x0001) == 0x0001) {
                  iborder = 0x0001;
               }
            }
         }
         if ( (iborder & 0x0001) == 0 && iilft >= 0) { //same size
            int nl = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilft, hash);
            int levcheck = -1;
            if (nl-noffset >= 0 && nl-noffset < ncells) {
               levcheck = level[nl-noffset];
            } else if (nl >= 0 && nl-ncells-noffset >= 0 && nl-ncells-noffset < nbsize_local) {
               levcheck = border_cell_level[nl-ncells-noffset];
            }
            if (nl >= (int)(ncells+noffset) && levcheck == lev && (border_cell_needed[nl-ncells-noffset] & 0x0001) == 0x0001) {
               iborder = 0x0001;
            } else if (lev != 0 && iilft-(iicur-iilft) >= 0){      // coarser neighbor
               iilft -= iicur-iilft;
               int jjlft = (jj/2)*2*levmult-jminsize;
               nl = read_hash(hash_method, hash_table_size, AA, BB, jjlft*(imaxsize-iminsize)+iilft, hash);
               levcheck = -1;
               if (nl-noffset >= 0 && nl-noffset < ncells) {
                  levcheck = level[nl-noffset];
               } else if (nl >= 0 && nl-ncells-noffset >= 0 && nl-ncells-noffset < nbsize_local) {
                  levcheck = border_cell_level[nl-ncells-noffset];
               }
               // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
               if (nl  >= (int)(ncells+noffset) && levcheck == lev-1 && (border_cell_needed[nl-ncells-noffset] & 0x0001) == 0x0001) {
                  iborder = 0x0001;
               }
            }
         }
      }

      // Test for cell to right
      if (iirht < imaxsize-iminsize && iirht >= 0 && jjcur >= 0 && jjtop < jmaxsize-jminsize) {
         // right neighbor -- finer, same size and coarser
         int nr = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iirht, hash);
         if (nr >= (int)(ncells+noffset) && (border_cell_needed[nr-ncells-noffset] & 0x0002) == 0x0002) {
            iborder = 0x0002;
         } else if (lev != levmx){
            // right neighbor -- finer right top test
            int jjtopfiner = (jjcur+jjtop)/2;
            int nrt = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iirht, hash);
            if (nrt >= (int)(ncells+noffset) && (border_cell_needed[nrt-ncells-noffset] & 0x0002) == 0x0002) {
               iborder = 0x0002;
            }
         }
         if ( (iborder & 0x0002) == 0  && lev != 0) { // test for coarser, but not directly right
            int jjrhtcoarser = (jj/2)*2*levmult-jminsize;
            if (jjrhtcoarser != jjcur) {
               int nr = read_hash(hash_method, hash_table_size, AA, BB, jjrhtcoarser*(imaxsize-iminsize)+iirht, hash);
               int levcheck = -1;
               if (nr-noffset >= 0 && nr-noffset < ncells) {
                  levcheck = level[nr-noffset];
               } else if (nr >= 0  && nr-ncells-noffset >= 0 && nr-ncells-noffset < nbsize_local) {
                  levcheck = border_cell_level[nr-ncells-noffset];
               }
               if (nr >= (int)(ncells+noffset) && levcheck == lev-1 && (border_cell_needed[nr-ncells-noffset] & 0x0002) == 0x0002) {
                  iborder = 0x0002;
               }
            }
         }
      }

      // Test for cell to bottom
      if (iicur >= 0 && (iicur+iirht)/2 < imaxsize-iminsize && jjcur-(jjcur-jjbot)/2 >= 0 && jjcur-(jjcur-jjbot)/2 < jmaxsize-jminsize){
         // Check for finer cell below and left side
         if (lev != levmx){                                // finer neighbor
            int jjbotfiner = jjcur-(jjcur-jjbot)/2;
            int nb = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iicur, hash);
            if (nb >= (int)(ncells+noffset) && (border_cell_needed[nb-ncells-noffset] & 0x0004) == 0x0004) {
               iborder = 0x0004;
            } else {
               // Also check for finer cell below and right side
               int iirhtfiner = (iicur+iirht)/2;
               int nbr = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
               if (nbr >= (int)(ncells+noffset) && (border_cell_needed[nbr-ncells-noffset] & 0x0004) == 0x0004) {
                  iborder = 0x0004;
               }
            }
         }
         if ( (iborder & 0x0004) == 0 && jjbot >= 0) { //same size
            int nb = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iicur, hash);
            int levcheck = -1;
            if (nb-noffset >= 0 && nb-noffset < ncells) {
               levcheck = level[nb-noffset];
            } else if (nb >= 0  && nb-ncells-noffset >= 0 && nb-ncells-noffset < nbsize_local) {
               levcheck = border_cell_level[nb-ncells-noffset];
            }
            if (nb >= (int)(ncells+noffset) && levcheck == lev && (border_cell_needed[nb-ncells-noffset] & 0x0004) == 0x0004) {
               iborder = 0x0004;
            } else if (lev != 0 && jjbot-(jjcur-jjbot) >= 0){      // coarser neighbor
               jjbot -= jjcur-jjbot;
               int iibot = (ii/2)*2*levmult-iminsize;
               nb = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iibot, hash);
               levcheck = -1;
               if (nb-noffset >= 0 && nb-noffset < ncells) {
                  levcheck = level[nb-noffset];
               } else if (nb >= 0  && nb-ncells-noffset >= 0 && nb-ncells-noffset < nbsize_local) {
                  levcheck = border_cell_level[nb-ncells-noffset];
               }
               // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
               if (nb >= (int)(ncells+noffset) && levcheck == lev-1 && (border_cell_needed[nb-ncells-noffset] & 0x0004) == 0x0004) {
                  iborder = 0x0004;
               }
            }
         }
      }

      // Test for cell to top
      if (iirht < imaxsize-iminsize && iicur >= 0 && jjtop >= 0 && jjtop < jmaxsize-jminsize) {
         // top neighbor -- finer, same size and coarser
         int nt = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iicur, hash);
         if (nt  >= (int)(ncells+noffset) && (border_cell_needed[nt-ncells-noffset] & 0x0008) == 0x0008) {
            iborder = 0x0008;
         } else if (lev != levmx){
            int iirhtfiner = (iicur+iirht)/2;
            int ntr = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
            if ( ntr >= (int)(ncells+noffset) && (border_cell_needed[ntr-ncells-noffset] & 0x0008) == 0x0008) {
               iborder = 0x0008;
            }
         }
         if ( (iborder & 0x0008) == 0  && lev != 0) { // test for coarser, but not directly above
            int iitopcoarser = (ii/2)*2*levmult-iminsize;
            if (iitopcoarser != iicur) {
               int nb = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iitopcoarser, hash);
               int levcheck = -1;
               if (nb-noffset >= 0 && nb-noffset < ncells) {
                  levcheck = level[nb-noffset];
               } else if (nb >= 0 && nb-ncells-noffset >= 0 && nb-ncells-noffset < nbsize_local) {
                  levcheck = border_cell_level[nb-ncells-noffset];
               }
               if (nb-noffset >= (int)(ncells-noffset) && levcheck == lev-1 && (border_cell_needed[nb-ncells-noffset] & 0x0008) == 0x0008) {
                  iborder = 0x0008;
               }
            }
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
                          const int   isize,               // 0 
                          const int   ncells,              // 1 
                          const int   noffset,             // 2 
                          const int   levmx,               // 3 
                          const int   imax,                // 4 
                          const int   jmax,                // 5 
                 __global const int4  *sizes,              // 6
                 __global const int   *levtable,           // 7
                 __global const int   *lev_ibeg,           // 8
                 __global const int   *lev_iend,           // 9
                 __global const int   *lev_jbeg,           // 10
                 __global const int   *lev_jend,           // 11
                 __global const int   *border_cell_i,      // 12
                 __global const int   *border_cell_j,      // 13
                 __global const int   *border_cell_level,  // 14
                 __global const int   *border_cell_num,    // 15
                 __global       int   *border_cell_needed, // 16
                 __global const ulong *hash_header,        // 17
                 __global       int   *hash)               // 18
{
   const uint giX = get_global_id(0);

   if (giX >= isize) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int lev = border_cell_level[giX];
   int levmult = levtable[levmx-lev];
   int ii = border_cell_i[giX]*levmult-iminsize;
   int jj = border_cell_j[giX]*levmult-jminsize;

   write_hash(hash_method, hash_table_size, AA, BB, -(ncells+giX), jj*(imaxsize-iminsize)+ii, hash);
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
                          const int   ncells_ghost,    // 0 
                          const int   levmx,           // 1 
                          const int   imax,            // 2 
                          const int   jmax,            // 3 
                 __global const int4  *sizes,          // 4 
                 __global const int   *levtable,       // 5 
                 __global const int   *i,              // 6
                 __global const int   *j,              // 7
                 __global const int   *level,          // 8
                 __global const ulong *hash_header,    // 9
                 __global const int   *hash,           // 10
                 __global       int   *nlft,           // 11
                 __global       int   *nrht,           // 12
                 __global       int   *nbot,           // 13
                 __global       int   *ntop)           // 14
{
   const uint giX  = get_global_id(0);

   if (giX >= ncells_ghost) return;

   const int hash_method       = (int)hash_header[0];
   const ulong hash_table_size =      hash_header[1];
   const ulong AA              =      hash_header[2];
   const ulong BB              =      hash_header[3];

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

   int nlftval = nlft[giX];
   int nrhtval = nrht[giX];
   int nbotval = nbot[giX];
   int ntopval = ntop[giX];

   int iicur = ii*levmult-iminsize;
   int iilft = max( (ii-1)*levmult, 0         )-iminsize;
   int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
   int jjcur = jj*levmult-jminsize;
   int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
   int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

   if (nlftval == -1){
      // Taking care of boundary cells
      // Force each boundary cell to point to itself on its boundary direction
      if (iicur <    1*levtable[levmx]  -iminsize) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);

      // Boundary cells next to corner boundary need special checks
      if (iicur ==    1*levtable[levmx]-iminsize &&  (jjcur < 1*levtable[levmx]-jminsize || jjcur >= jmax*levtable[levmx]-jminsize ) ) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);

      // need to check for finer neighbor first
      // Right and top neighbor don't change for finer, so drop through to same size
      // Left and bottom need to be half of same size index for finer test
      if (lev != levmx) {
         int iilftfiner = iicur-(iicur-iilft)/2;
         if (nlftval == -1 && iilftfiner >= 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilftfiner, hash);
      }

      // same size neighbor
      if (nlftval == -1 && iilft >= 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iilft, hash);

      // Now we need to take care of special case where bottom and left boundary need adjustment since
      // expected cell doesn't exist on these boundaries if it is finer than current cell
      if (lev != levmx && jjcur < 1*levtable[levmx]) {
         if (nlftval == -1) {
            int iilftfiner = iicur-(iicur-iilft)/2;
            int jjtopfiner = (jjcur+jjtop)/2;
            if (jjtopfiner < jmaxsize-jminsize && iilftfiner >= 0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
         }
      }

      // coarser neighbor
      if (lev != 0){
         if (nlftval == -1) {
            int iilftcoarser = iilft - (iicur-iilft);
            int jjlft = (jj/2)*2*levmult-jminsize;
            if (iilftcoarser >=0) nlftval = read_hash(hash_method, hash_table_size, AA, BB, jjlft*(imaxsize-iminsize)+iilftcoarser, hash);
         }
      }

      if (nlftval != -1) nlft[giX] = nlftval;
   }

   if (nrhtval == -1) {
      // Taking care of boundary cells
      // Force each boundary cell to point to itself on its boundary direction
      if (iicur > imax*levtable[levmx]-1-iminsize) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);

      // Boundary cells next to corner boundary need special checks
      if (iirht == imax*levtable[levmx]-iminsize &&  (jjcur < 1*levtable[levmx]-jminsize || jjcur >= jmax*levtable[levmx]-jminsize ) ) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);

      // same size neighbor
      if (nrhtval == -1 && iirht < imaxsize-iminsize) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iirht, hash);

      // Now we need to take care of special case where bottom and left boundary need adjustment since
      // expected cell doesn't exist on these boundaries if it is finer than current cell
      if (lev != levmx && jjcur < 1*levtable[levmx]) {
         if (nrhtval == -1) {
            int jjtopfiner = (jjcur+jjtop)/2;
            if (jjtopfiner < jmaxsize-jminsize && iirht < imaxsize-iminsize) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjtopfiner*(imaxsize-iminsize)+iirht, hash);
         }
      }

      // coarser neighbor
      if (lev != 0){
         if (nrhtval == -1) {
            int jjrht = (jj/2)*2*levmult-jminsize;
            if (iirht < imaxsize-iminsize) nrhtval = read_hash(hash_method, hash_table_size, AA, BB, jjrht*(imaxsize-iminsize)+iirht, hash);
         }
      }
      if (nrhtval != -1) nrht[giX] = nrhtval;
   }

   if (nbotval == -1) {
      // Taking care of boundary cells
      // Force each boundary cell to point to itself on its boundary direction
      if (jjcur <    1*levtable[levmx]  -jminsize) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);
      // Boundary cells next to corner boundary need special checks
      if (jjcur ==    1*levtable[levmx]-jminsize &&  (iicur < 1*levtable[levmx]-iminsize || iicur >= imax*levtable[levmx]-iminsize ) ) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);

      // need to check for finer neighbor first
      // Right and top neighbor don't change for finer, so drop through to same size
      // Left and bottom need to be half of same size index for finer test
      if (lev != levmx) {
         int jjbotfiner = jjcur-(jjcur-jjbot)/2;
         if (nbotval == -1 && jjbotfiner >= 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iicur, hash);
      }

      // same size neighbor
      if (nbotval == -1 && jjbot >=0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbot*(imaxsize-iminsize)+iicur, hash);

      // Now we need to take care of special case where bottom and left boundary need adjustment since
      // expected cell doesn't exist on these boundaries if it is finer than current cell
      if (lev != levmx && iicur < 1*levtable[levmx]) {
         if (nbotval == -1) {
            int iirhtfiner = (iicur+iirht)/2;
            int jjbotfiner = jjcur-(jjcur-jjbot)/2;
            if (jjbotfiner >= 0 && iirhtfiner < imaxsize-iminsize) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
         }
      }

      // coarser neighbor
      if (lev != 0){
         if (nbotval == -1) {
            int jjbotcoarser = jjbot - (jjcur-jjbot);
            int iibot = (ii/2)*2*levmult-iminsize;
            if (jjbotcoarser >= 0 && iibot >= 0) nbotval = read_hash(hash_method, hash_table_size, AA, BB, jjbotcoarser*(imaxsize-iminsize)+iibot, hash);
         }
      }
      if (nbotval != -1) nbot[giX] = nbotval;
   }
   
   if (ntopval == -1) {
      // Taking care of boundary cells
      // Force each boundary cell to point to itself on its boundary direction
      if (jjcur > jmax*levtable[levmx]-1-jminsize) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);
      // Boundary cells next to corner boundary need special checks
      if (jjtop == jmax*levtable[levmx]-jminsize &&  (iicur < 1*levtable[levmx]-iminsize || iicur >= imax*levtable[levmx]-iminsize ) ) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjcur*(imaxsize-iminsize)+iicur, hash);

      // same size neighbor
      if (ntopval == -1 && jjtop < jmaxsize-jminsize) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iicur, hash);
  
      if (iicur < 1*levtable[levmx]) {
         if (ntopval == -1) {
            int iirhtfiner = (iicur+iirht)/2;
            if (jjtop < jmaxsize-jminsize && iirhtfiner < imaxsize-iminsize) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
         }
      }
  
      // coarser neighbor
      if (lev != 0){
         if (ntopval == -1) {
            int iitop = (ii/2)*2*levmult-iminsize;
            if (jjtop < jmaxsize-jminsize && iitop < imaxsize-iminsize) ntopval = read_hash(hash_method, hash_table_size, AA, BB, jjtop*(imaxsize-iminsize)+iitop, hash);
         }
      }
      if (ntopval != -1) ntop[giX] = ntopval;
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

__kernel void finish_reduction_count2_cl(
                 const int   isize,
        __global       int2 *redscratch,
        __global       int2 *result,
        __local        int2 *itile)
{
   const uint tiX = get_local_id(0);
   const uint ntX = get_local_size(0);

   int giX = tiX;

   itile[tiX].s01 = 0;

   if (tiX < isize) itile[tiX].s01 = redscratch[giX].s01;

   for (giX += ntX; giX < isize; giX += ntX) {
     itile[tiX].s01 += redscratch[giX].s01;
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        itile[tiX].s01 += itile[tiX+offset].s01;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }

   if (tiX < 32){
      itile[tiX].s01 += itile[tiX+32].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+16].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+8].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+4].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+2].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+1].s01;
   }

   if (tiX == 0) {
     redscratch[0].s01 = itile[0].s01;
     result[0].s01 = itile[0].s01;
   }
}

__kernel void finish_reduction_count_cl(
                 const int   isize,
        __global       int *redscratch,
        __global       int *result,
        __local        int *itile)
{
   const uint tiX = get_local_id(0);
   const uint ntX = get_local_size(0);

   int giX = tiX;

   itile[tiX] = 0;

   if (tiX < isize) itile[tiX] = redscratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
     itile[tiX] += redscratch[giX];
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        itile[tiX] += itile[tiX+offset];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }

   if (tiX < 32){
      itile[tiX] += itile[tiX+32];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+1];
   }

   if (tiX == 0) {
     redscratch[0] = itile[0];
     result[0] = itile[0];
   }
}

__kernel void finish_reduction_scan2_cl(
                 const int   isize,
        __global       uint  *ioffset,
        __global       uint  *result,
        __local        uint  *itile)
{
   const uint tiX = get_local_id(0);
   const uint gID = get_group_id(0);
   const uint ntX = get_local_size(0);

   const uint lane = tiX & 31;
   const uint warpID = tiX >> 5;
   const uint EPT = (isize+ntX-1)/ntX; //elements_per_thread;
   const uint BLOCK_SIZE = EPT * ntX;

   int reduceValue = 0;

// #pragma unroll 4
   for(uint i = 0; i < EPT; ++i)
   {
      uint offsetIdx = i * ntX + (gID * BLOCK_SIZE) + tiX;

#ifdef IS_NVIDIA
//    if (offsetIdx >= isize) return;
#endif

      // Step 1: Read ntX elements from global (off-chip) memory to local memory (on-chip)
      int input = 0;
      if (offsetIdx < isize) input = ioffset[offsetIdx];           
      itile[tiX] = input;           
      barrier(CLK_LOCAL_MEM_FENCE);

      // Step 2: Perform scan on ntX elements
      int val = scan_workgroup_exclusive(itile, tiX, lane, warpID);

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

   tile[tiX].s0 = i[giX]   *levtable[levmx-lev]  ; // imincalc
   tile[tiX].s1 =(i[giX]+1)*levtable[levmx-lev]-1; // imaxcalc

   tile[tiX].s2 = j[giX]   *levtable[levmx-lev]  ; // jmincalc
   tile[tiX].s3 =(j[giX]+1)*levtable[levmx-lev]-1; // jmaxcalc

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
          const    int       isize,
          const    real_t    xmin,
          const    real_t    ymin,
          __global real_t    *lev_deltax,
          __global real_t    *lev_deltay,
          __global spatial_t *x,
          __global spatial_t *dx,
          __global spatial_t *y,
          __global spatial_t *dy,
          __global int       *level,
          __global int       *i,
          __global int       *j)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int lev = level[giX];

   x[giX] = isize;
   x[giX]  = xmin + lev_deltax[lev] * (real_t)i[giX];
   dx[giX] =        lev_deltax[lev];
   y[giX]  = ymin + lev_deltay[lev] * (real_t)j[giX];
   dy[giX] =        lev_deltay[lev];
}

#ifndef MINIMUM_PRECISION
__kernel void do_load_balance_double_cl(
                  const int    ncells,
                  const int    lower_block_size,
                  const int    middle_block_size,
                  const int    middle_block_start,
         __global const double *state_var_old,
         __global const double *state_var_lower,
         __global const double *state_var_upper,
         __global       double *state_var_new)
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
#endif

__kernel void do_load_balance_float_cl(
                  const int   ncells,
                  const int   lower_block_size,
                  const int   middle_block_size,
                  const int   middle_block_start,
         __global const float *state_var_old,
         __global const float *state_var_lower,
         __global const float *state_var_upper,
         __global       float *state_var_new)
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

// XXX Refinement Smoothing XXX
__kernel void refine_smooth_cl(
                 const int    ncells,       // 0  Total number of cells.
                 const int    ncells_ghost, // 1  Total number of cells.
                 const int    levmx,        // 2  Maximum level
        __global const int   *nlft,         // 3  Array of left neighbors.
        __global const int   *nrht,         // 4  Array of right neighbors.
        __global const int   *nbot,         // 5  Array of bottom neighbors.
        __global const int   *ntop,         // 6  Array of top neighbors.
        __global const int   *level,        // 7  Array of level information.
        __global const int   *celltype,     // 8  Array of celltype information.
        __global const int   *mpot_old,     // 9  Array of mesh potential information.
        __global       int   *mpot,         // 10 Array of mesh potential information.
        __global       int   *redscratch,   // 11  Array of number of cells smoothed per tile
        __global       int   *result,       // 12
        __local        int   *itile)        // 13  Tile size in int2.
{

   const uint giX  = get_global_id(0);
   const uint tiX  = get_local_id(0);

   const uint group_id = get_group_id(0);

   const uint ntX  = get_local_size(0);

   itile[tiX] = 0;

   if(giX >= ncells)
      return;

//   setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);

   barrier (CLK_LOCAL_MEM_FENCE);

   int nl, nr, nt, nb;
   int nlt, nrt, ntr, nbr;
   int lev, ll, lr, lt, lb;
   int llt, lrt, ltr, lbr;

   int mpotval = mpot_old[giX];
   int ctype   = celltype[giX];

   int new_count = 0;
   if (mpotval > 0){
      new_count = (ctype == REAL_CELL) ? 3 : 1;
   }

   int ic = giX;
   mpot[giX] = mpot_old[giX];

   lev = level[ic];
   if(mpot_old[ic] <= 0) {

      nl = nlft[ic];
      if (nl >= 0 && nl < ncells_ghost) {
         ll = level[nl];
         if(mpot_old[nl] > 0) ll++;
   
         if(ll - lev > 1) {
            mpot[ic]=1;
            new_count++;
            lev = levmx;
         }

         ll = level[nl];
         if (ll > lev) {
            nlt = ntop[nl];
            if (nlt >= 0 && nlt < ncells_ghost) {
               llt = level[nlt];
               if(mpot_old[nlt] > 0) llt++;

               if(llt - lev > 1) {
                  mpot[ic]=1;
                  new_count++;
                  lev = levmx;
               }
            }
         }
      }

      nr = nrht[ic];
      if (nr >= 0 && nr < ncells_ghost) {
         lr = level[nr];
         if(mpot_old[nr] > 0) lr++;
   
         if(lr - lev > 1) {
            mpot[ic]=1;
            new_count++;
            lev = levmx;
         }

         lr = level[nr];
         if (lr > lev) {
            nrt = ntop[nr];
            if (nrt >= 0 && nrt < ncells_ghost) {
               lrt = level[nrt];
               if(mpot_old[nrt] > 0) lrt++;

               if(lrt - lev > 1) {
                  mpot[ic]=1;
                  new_count++;
                  lev = levmx;
               }
            }
         }
      }

      nt = ntop[ic];
      if (nt >= 0 && nt < ncells_ghost) {
         lt = level[nt];
         if(mpot_old[nt] > 0) lt++;
   
         if(lt - lev > 1) {
            mpot[ic]=1;
            new_count++;
            lev = levmx;
         }

         lt = level[nt];
         if (lt > lev) {
            ntr = nrht[nt];
            if (ntr >= 0 && ntr < ncells_ghost) {
               ltr = level[ntr];
               if(mpot_old[ntr] > 0) ltr++;

               if(ltr - lev > 1) {
                  mpot[ic]=1;
                  new_count++;
                  lev = levmx;
               }
            }
         }
      }

      nb = nbot[ic];
      if (nb >= 0 && nb < ncells_ghost) {
         lb = level[nb];
         if(mpot_old[nb] > 0) lb++;
   
         if(lb - lev > 1) {
            mpot[ic]=1;
            new_count++;
            lev = levmx;
         }

         lb = level[nb];
         if (lb > lev) {
            nbr = nrht[nb];
            if (nbr >= 0 && nbr < ncells_ghost) {
               lbr = level[nbr];
               if(mpot_old[nbr] > 0) lbr++;

               if(lbr - lev > 1) {
                  mpot[ic]=1;
                  new_count++;
                  lev = levmx;
               }
            }
         }
      }
   }

   itile[tiX] = new_count;

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
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */
      itile[tiX] += itile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */
      itile[tiX] += itile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */
      itile[tiX] += itile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */
      itile[tiX] += itile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */
      itile[tiX] += itile[tiX+1]; }

   if (tiX == 0) {
     redscratch[group_id] = itile[0];
     result[0] = itile[0];
   }

}

__kernel void coarsen_smooth_cl(
                 const int    ncells,       // 0  Total number of cells.
        __global const int   *nlft,         // 1  Array of left neighbors.
        __global const int   *nrht,         // 2  Array of right neighbors.
        __global const int   *nbot,         // 3  Array of bottom neighbors.
        __global const int   *ntop,         // 4  Array of top neighbors.
        __global const int   *i,            // 5  Array of i values.
        __global const int   *j,            // 6  Array of j values.
        __global const int   *level,        // 7  Array of level information.
        __global const int   *mpot_old,     // 8  Array of mpot_old information.
        __global       int   *mpot)         // 9  Array of mpot information.
{
   const uint giX  = get_global_id(0);

   if(giX >= ncells)
      return;

   int ival    = i[giX];
   int jval    = j[giX];
   int lev     = level[giX];
   int mpotval = mpot_old[giX];

   if (mpotval < 0) {
      if (        is_upper_right(ival,jval) ) {
         int nr = nrht[giX];
         int lr = level[nr];
         if (mpot_old[nr] > 0) lr++;
         int nt = ntop[giX];
         int lt = level[nt];
         if (mpot_old[nt] > 0) lt++;
         if (lr > lev || lt > lev) mpotval = 0;
      } else if ( is_upper_left(ival,jval) ) {
         int nl = nlft[giX];
         int ll = level[nl];
         if (mpot_old[nl] > 0) ll++;
         int nt = ntop[giX];
         int lt = level[nt];
         if (mpot_old[nt] > 0) lt++;
         if (ll > lev || lt > lev) mpotval = 0;
      } else if ( is_lower_right(ival,jval) ) {
         int nr = nrht[giX];
         int lr = level[nr];
         if (mpot_old[nr] > 0) lr++;
         int nb = nbot[giX];
         int lb = level[nb];
         if (mpot_old[nb] > 0) lb++;
         if (lr > lev || lb > lev) mpotval = 0;
      } else if ( is_lower_left(ival,jval) ) {
         int nl = nlft[giX];
         int ll = level[nl];
         if (mpot_old[nl] > 0) ll++;
         int nb = nbot[giX];
         int lb = level[nb];
         if (mpot_old[nb] > 0) lb++;
         if (ll > lev || lb > lev) mpotval = 0;
      }
   }

   mpot[giX] = mpotval;
}

__kernel void coarsen_check_block_cl(
                 const int    ncells,       // 0  Total number of cells.
        __global const int   *nlft,         // 1  Array of left neighbors.
        __global const int   *nrht,         // 2  Array of right neighbors.
        __global const int   *nbot,         // 3  Array of bottom neighbors.
        __global const int   *ntop,         // 4  Array of top neighbors.
        __global const int   *i,            // 5  Array of i values.
        __global const int   *j,            // 6  Array of j values.
        __global const int   *level,        // 7  Array of level information.
        __global const int   *celltype,     // 8  Array of celltype information.
        __global const int   *mpot_old,     // 9  Array of mpot_old information.
        __global       int   *mpot,         // 10 Array of mpot information.
        __global       int   *redscratch,   // 10 Reduction scratch
        __global       int   *result,       // 11 Reduction result
        __local        int   *itile)        // 12 itile scratch.
{
   const uint giX      = get_global_id(0);
   const uint tiX      = get_local_id(0);
   const uint group_id = get_group_id(0);
   const uint ntX      = get_local_size(0);

   itile[tiX] = 0;

   if(giX >= ncells)
      return;

   int mpotval = mpot_old[giX];

   if (mpotval < 0) {
      int n1, n2, n3;
      int ival = i[giX];
      int jval = j[giX];
      int lev  = level[giX];
      if ( is_upper_right(ival,jval) ) {
         n1 = nbot[giX];
         n2 = nlft[giX];
         n3 = nlft[n1];
      } else if ( is_upper_left(ival,jval) ) {
         n1 = nbot[giX];
         n2 = nrht[giX];
         n3 = nrht[n1];
      } else if ( is_lower_right(ival,jval) ) {
         n1 = ntop[giX];
         n2 = nlft[giX];
         n3 = nlft[n1];
      } else if ( is_lower_left(ival,jval) ) {
         n1 = ntop[giX];
         n2 = nrht[giX];
         n3 = nrht[n1];
      }
      if (n3 < 0) {
         mpotval = 0;
      } else {
         int lev1 = level[n1];
         int lev2 = level[n2];
         int lev3 = level[n3];
         if (mpot_old[n1] > 0) lev1++;
         if (mpot_old[n2] > 0) lev2++;
         if (mpot_old[n3] > 0) lev3++;

         if (mpot_old[n1] != -1 || lev1 != lev ||
             mpot_old[n2] != -1 || lev2 != lev ||
             mpot_old[n3] != -1 || lev3 != lev) {
            mpotval = 0;
         }
      }

      if (mpotval < 0) {
         int ctype = celltype[giX];
         if (ctype == REAL_CELL) {
            if (! is_lower_left(ival,jval) ) itile[tiX] = 1;
         } else {
            if (! is_upper_right(ival,jval) || is_lower_left(ival,jval) ) itile[tiX] = 1;
         }
      }
   }

   mpot[giX] = mpotval;

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset = ntX >> 1; offset > 32; offset >>= 1) {
      if (tiX < offset) {
         itile[tiX] += itile[tiX+offset];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }

   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32) {
      itile[tiX] += itile[tiX+32];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX] += itile[tiX+1];
   }

   if (tiX == 0) {
     redscratch[group_id] = itile[0];
     result[0] = itile[0];
   }
}

__kernel void rezone_all_cl(
                 const int    isize,         // 0 
                 const int    stencil,       // 1 
                 const int    levmx,         // 2
        __global const int    *mpot,         // 3   Array of mesh potential information.
        __global const int    *level,        // 4
        __global const int    *i,            // 5
        __global const int    *j,            // 6
        __global const int    *celltype,     // 7
        __global       int    *level_new,    // 8
        __global       int    *i_new,        // 9
        __global       int    *j_new,        // 10
        __global       int    *celltype_new, // 11
        __global const uint   *ioffset,      // 12   
        __global       uint   *indexoffset,  // 13   
        __global const real_t *lev_dx,       // 14
        __global const real_t *lev_dy,       // 15
        __global const int    *levtable,     // 16
        __global const int    *ijadd,        // 17
        __local        uint   *itile)        // 18
        //__local        real2  *tile)         // 19
{ 
   uint tiX = get_local_id(0);
   uint giX = get_global_id(0);

   uint ntX = get_local_size(0);
   uint group_id = get_group_id(0);

   const uint lane   = tiX & 31;
   const uint warpid = tiX >> 5;

   // Step 1: load global data into tile
   uint indexval=0;
   itile[tiX] = 0;

   if (giX >= isize) return;

   int mpotval = mpot[giX];
   int ctype = celltype[giX];
   int ival = i[giX];
   int jval = j[giX];

   int do_coarsening = 0;
   if (mpotval < 0){
      indexval = 0;
      if (is_lower_left(ival, jval) ) do_coarsening = 1;
      if (ctype != REAL_CELL && is_upper_right(ival, jval) ) do_coarsening = 1;
      if (do_coarsening) indexval = 1;
   } else {
      if (celltype[giX] == REAL_CELL){
         indexval = mpot[giX] ? 4 : 1;
      } else {
         indexval = mpot[giX] ? 2 : 1;
      }
   }
   itile[tiX] = indexval;
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 2: scan each warp
   indexval = scan_warp_exclusive(itile, tiX, lane);
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 3: Collect per-warp sums
   if (lane == 31) itile[warpid] = itile[tiX];
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 4: Use 1st warp to scan per-warp sums
   if (warpid == 0) scan_warp_inclusive(itile, tiX, lane);
   barrier(CLK_LOCAL_MEM_FENCE);

   // Step 5: Accumulate results from steps 2 and 4
   if (warpid > 0) indexval += itile[warpid-1];
   barrier(CLK_LOCAL_MEM_FENCE);

   indexval += ioffset[group_id];

   int lev = level[giX];
   if (mpotval > 0) lev++;

   if (mpotval == 0) {
      i_new[indexval] = ival;
      j_new[indexval] = jval;
      level_new[indexval] = lev;
      celltype_new[indexval] = ctype;
   } else if (mpotval < 0) {
      if (do_coarsening) {
         i_new[indexval] = ival/2;
         j_new[indexval] = jval/2;
         level_new[indexval] = lev - 1;
         celltype_new[indexval] = ctype;
      }
   } else if (mpotval > 0) {
      level_new[indexval] = lev;
      level_new[indexval+1] = lev;
      celltype_new[indexval] = ctype;
      celltype_new[indexval+1] = ctype;
      if (ctype == REAL_CELL) {
         level_new[indexval+2] = lev;
         level_new[indexval+3] = lev;
         celltype_new[indexval+2] = ctype;
         celltype_new[indexval+3] = ctype;
          
         int order[4] = {SW, SE, NW, NE};
         if (stencil) {
#ifdef __OLD_STENCIL__
            spatial_t nx[3], ny[3];
            int ifirst = ijadd[0];
            int iend   = ijadd[1];
            int jfirst = ijadd[2];
            int jend   = ijadd[3];
            int level_first  = ijadd[4];
            int level_end    = ijadd[5];

//             x[ic]  = xmin + lev_deltax[level[ic]] * (spatial_t)(i[ic] - ibase);

            if (giX != 0) {
               nx[0] = lev_dx[level[giX-1]] * (spatial_t)i[giX-1];
               ny[0] = lev_dy[level[giX-1]] * (spatial_t)j[giX-1];
            } else {
               nx[0] = lev_dx[level_first] * (spatial_t)ifirst;
               ny[0] = lev_dy[level_first] * (spatial_t)jfirst;
            }

            nx[1] = lev_dx[level[giX  ]] * (spatial_t)i[giX  ];
            ny[1] = lev_dy[level[giX  ]] * (spatial_t)j[giX  ];
            if (giX != isize-1) {
               nx[2] = lev_dx[level[giX+1]] * (spatial_t)i[giX+1];
               ny[2] = lev_dy[level[giX+1]] * (spatial_t)j[giX+1];
            } else {
               nx[2] = lev_dx[level_end] * (spatial_t)iend;
               ny[2] = lev_dy[level_end] * (spatial_t)jend;
            }

            //  Figure out relative orientation of the neighboring cells.  We are
            //  are aided in this because the Hilbert curve only has six possible
            //  ways across the cell:  four Ls and two straight lines.  Then
            //  refine the cell according to the relative orientation and order
            //  according to the four-point Hilbert stencil.
            if      (nx[0] < nx[1] and ny[2] < ny[1])   //  southwest L, forward order
            {  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE; }
            else if (nx[2] < nx[1] and ny[0] < ny[1])   //  southwest L, reverse order
            {  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW; }
            else if (nx[0] > nx[1] and ny[2] < ny[1])   //  southeast L, forward order
            {  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW; }
            else if (nx[2] > nx[1] and ny[0] < ny[1])   //  southeast L, reverse order
            {  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE; }
            else if (nx[0] > nx[1] and ny[2] > ny[1])   //  northeast L, forward order
            {  order[0] = SE; order[1] = SW; order[2] = NW; order[3] = NE; }
            else if (nx[2] > nx[1] and ny[0] > ny[1])   //  northeast L, reverse order
            {  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE; }
            else if (nx[0] < nx[1] and ny[2] > ny[1])   //  northwest L, forward order
            {  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW; }
            else if (nx[2] < nx[1] and ny[0] > ny[1])   //  northwest L, reverse order
            {  order[0] = NW; order[1] = NE; order[2] = SE; order[3] = SW; }
            else if (nx[0] > nx[1] and nx[1] > nx[2])   //  straight horizontal, forward order
            {  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW; }
            else if (nx[0] < nx[1] and nx[1] < nx[2])   //  straight horizontal, reverse order
            {  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE; }
            else if (ny[0] > ny[1] and ny[1] > ny[2])   //  straight vertical, forward order
            {  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE; }
            else if (ny[0] < ny[1] and ny[1] < ny[2])   //  straight vertical, reverse order
            {  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW; }
            else                                        //  other, default to z-order
            {  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE; }
#endif

#ifdef __NEW_STENCIL__
// XXX Need conditional check on giX != 0, ncells XXX
               int ir[3],   // First i index at finest level of the mesh
                   jr[3];   // First j index at finest level of the mesh
               // Cell's Radius at the Finest level of the mesh
               int crf = levtable[levmx-level[giX]];
               int ifirst = ijadd[0];
               int iend   = ijadd[1];
               int jfirst = ijadd[2];
               int jend   = ijadd[3];
               int level_first  = ijadd[4];
               int level_end    = ijadd[5];
               if (giX != 0) {
                  ir[0] = i[giX - 1] * levtable[levmx-level[giX - 1]];
                  jr[0] = j[giX - 1] * levtable[levmx-level[giX - 1]];
               } else {
                  ir[0] = ifirst * levtable[levmx-level_first];
                  jr[0] = jfirst * levtable[levmx-level_first];
               }
               ir[1] = i[giX    ] * crf;
               jr[1] = j[giX    ] * crf;
               if (giX != isize-1) {
                  ir[2] = i[giX + 1] * levtable[levmx-level[giX + 1]];
                  jr[2] = j[giX + 1] * levtable[levmx-level[giX + 1]];
               } else {
                  ir[2] = iend * levtable[levmx-level_end];
                  jr[2] = jend * levtable[levmx-level_end];
               }

               int dir_in  = ir[1] - ir[0];
               int dir_out = ir[1] - ir[2];
               int djr_in  = jr[1] - jr[0];
               int djr_out = jr[1] - jr[2];

               char  in_direction = 'X';
               char out_direction = 'X';

               // Left In
               if( (djr_in == 0 && (dir_in == crf*HALF || dir_in == crf || dir_in == crf*TWO)) || (djr_in == -crf*HALF && dir_in == crf*HALF) || (djr_in == crf && dir_in == crf*TWO) ) {
                  in_direction = 'L';
               }
               // Bottom In
               else if( (dir_in == 0 && (djr_in == crf*HALF || djr_in == crf || djr_in == crf*TWO)) || (dir_in == -crf*HALF && djr_in == crf*HALF) || (dir_in == crf && djr_in == crf*TWO) ) {
                  in_direction = 'B';
               }
               // Right In
               else if( (dir_in == -crf && (djr_in == -crf*HALF || djr_in == 0 || (djr_in == crf && level[giX-1] < level[giX]))) ) {
                  in_direction = 'R';
               }
               // Top In
               else if( (djr_in == -crf && (dir_in == -crf*HALF || dir_in == 0 || (dir_in == crf && level[giX-1] < level[giX]))) ) {
                  in_direction = 'T';
               }
               // Further from the left
               else if( dir_in > 0 && djr_in == 0 ) {
                  in_direction = 'L';
               }
               // Further from the right
               else if( dir_in < 0 && djr_in == 0 ) {
                  in_direction = 'R';
               }
               // Further from the bottom
               else if( djr_in > 0 && dir_in == 0 ) {
                  in_direction = 'B';
               }
               // Further from the top
               else if( djr_in < 0 && dir_in == 0 ) {
                  in_direction = 'T';
               }
               // SW in; 'M'
               else if( dir_in > 0 && djr_in > 0) {
                  in_direction = 'M';
               }
               // NW in; 'W'
               else if( dir_in > 0 && djr_in < 0) {
                  in_direction = 'W';
               }
               // SE in; 'F'
               else if( dir_in < 0 && djr_in > 0) {
                  in_direction = 'F';
               }
               // NE in; 'E'
               else if( dir_in < 0 && djr_in < 0) {
                  in_direction = 'E';
               }

               // Left Out
               if( (djr_out == 0 && (dir_out == crf*HALF || dir_out == crf || dir_out == crf*TWO)) || (djr_out == -crf*HALF && dir_out == crf*HALF) || (djr_out == crf && dir_out == crf*TWO) ) {
                  out_direction = 'L';
               }
               // Bottom Out
               else if( (dir_out == 0 && (djr_out == crf*HALF || djr_out == crf || djr_out == crf*TWO)) || (dir_out == -crf*HALF && djr_out == crf*HALF) || (dir_out == crf && djr_out == crf*TWO) ) {
                  out_direction = 'B';
               }
               // Right Out
               else if( (dir_out == -crf && (djr_out == -crf*HALF || djr_out == 0 || (djr_out == crf && level[giX+1] < level[giX]))) ) {
                  out_direction = 'R';
               }
               // Top Out
               else if( (djr_out == -crf && (dir_out == -crf*HALF || dir_out == 0 || (dir_out == crf && level[giX+1] < level[giX]))) ) {
                  out_direction = 'T';
               }
               // Further from the left
               else if( dir_out > 0 && djr_out == 0 ) {
                  out_direction = 'L';
               }
               // Further from the right
               else if( dir_out < 0 && djr_out == 0 ) {
                  out_direction = 'R';
               }
               // Further from the bottom
               else if( djr_out > 0 && dir_out == 0 ) {
                  out_direction = 'B';
               }
               // Further from the top
               else if( djr_out < 0 && dir_out == 0 ) {
                  out_direction = 'T';
               }
               // SW out; 'M'
               else if( dir_out > 0 && djr_out > 0) {
                  out_direction = 'M';
               }
               // NW out; 'W'
               else if( dir_out > 0 && djr_out < 0) {
                  out_direction = 'W';
               }
               // SE out; 'F'
               else if( dir_out < 0 && djr_out > 0) {
                  out_direction = 'F';
               }
               // NE out; 'E'
               else if( dir_out < 0 && djr_out < 0) {
                  out_direction = 'E';
               }

               // Set the Stencil
               if(in_direction == 'L' && (out_direction == 'B' || out_direction == 'R' || out_direction == 'F')) {
                  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE;
               }
               else if(in_direction == 'L' && (out_direction == 'T' || out_direction == 'W' )) {
                  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW;
               }
               else if(in_direction == 'L' && out_direction == 'M') {
                  order[0] = NW; order[1] = NE; order[2] = SE; order[3] = SW;
               }
               else if(in_direction == 'L' && out_direction == 'E') {
                  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE;
               }

               else if(in_direction == 'B' && (out_direction == 'R' || out_direction == 'F' )) {
                  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE;
               }
               else if(in_direction == 'B' && (out_direction == 'L' || out_direction == 'T' || out_direction == 'W' )) {
                  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW;
               }
               else if(in_direction == 'B' && out_direction == 'M') {
                  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW;
               }
               else if(in_direction == 'B' && out_direction == 'E') {
                  order[0] = SW; order[1] = NW; order[2] = SE; order[3] = NE;
               }
               
               else if(in_direction == 'R' && (out_direction == 'T' || out_direction == 'L' || out_direction == 'W' )) {
                  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW;
               }
               else if(in_direction == 'R' && (out_direction == 'B' || out_direction == 'F' )) {
                  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE;
               }
               else if(in_direction == 'R' && out_direction == 'M') {
                  order[0] = NE; order[1] = NW; order[2] = SE; order[3] = SW;
               }
               else if(in_direction == 'R' && out_direction == 'E') {
                  order[0] = SE; order[1] = SW; order[2] = NW; order[3] = NE;
               }

               else if(in_direction == 'T' && (out_direction == 'L' || out_direction == 'W' )) {
                  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW;
               }
               else if(in_direction == 'T' && (out_direction == 'R' || out_direction == 'B' || out_direction == 'F' )) {
                  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE;
               }
               else if(in_direction == 'T' && out_direction == 'M') {
                  order[0] = NE; order[1] = SE; order[2] = NW; order[3] = SW;
               }
               else if(in_direction == 'T' && out_direction == 'E') {
                  order[0] = NW; order[1] = SW; order[2] = SE; order[3] = NE;
               }

               else if(in_direction == 'M' && (out_direction == 'L' || out_direction == 'W' || out_direction == 'T') ) {
                  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW;
               }
               else if(in_direction == 'M' && (out_direction == 'R' || out_direction == 'F' || out_direction == 'B') ) {
                  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE;
               }
               else if(in_direction == 'M' && out_direction == 'E') {
                  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE;
               }
 
               else if(in_direction == 'W' && (out_direction == 'L' || out_direction == 'M' || out_direction == 'B') ) {
                  order[0] = NW; order[1] = NE; order[2] = SE; order[3] = SW;
               }
               else if(in_direction == 'W' && (out_direction == 'R' || out_direction == 'E' || out_direction == 'T') ) {
                  order[0] = NW; order[1] = SW; order[2] = SE; order[3] = NE;
               }
               else if(in_direction == 'W' && out_direction == 'F') {
                  order[0] = NW; order[1] = NE; order[2] = SW; order[3] = SE;
               }

               else if(in_direction == 'F' && (out_direction == 'L' || out_direction == 'M' || out_direction == 'B') ) {
                  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW;
               }
               else if(in_direction == 'F' && (out_direction == 'R' || out_direction == 'E' || out_direction == 'T') ) {
                  order[0] = SE; order[1] = SW; order[2] = NW; order[3] = NE;
               }
               else if(in_direction == 'F' && out_direction == 'W') {
                  order[0] = SE; order[1] = NE; order[2] = SW; order[3] = NW;
               }

               else if(in_direction == 'E' && (out_direction == 'L' || out_direction == 'W' || out_direction == 'T') ) {
                  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW;
               }
               else if(in_direction == 'E' && (out_direction == 'R' || out_direction == 'F' || out_direction == 'B') ) {
                  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE;
               }
               else if(in_direction == 'E' && out_direction == 'M') {
                  order[0] = NE; order[1] = SE; order[2] = NW; order[3] = SW;
               }

               else { // Default to a knot 
                  order[0] = NW; order[1] = SE; order[2] = SW; order[3] = NE;
               }
               //  Determine the relative orientation of the neighboring cells.
               //  There are 12 possible ways across the cell: 4 Ls and 2 straight
               //  lines, each with 2 directions of traversal.
               //  Then the cell is refined and ordered according to the relative
               //  orientation and four-point Hilbert stencil.

               // XXX NOTE that the four-point stencil varies depending upon
               // the starting and ending point of the global Hilbert curve.
               // The stencil applied here assumes the start at (0,0) and the end
               // at (0,y_max). XXX WRONG
#endif

         } // End local stencil version

         for (int ii = 0; ii < 4; ii++){
#ifdef IS_NVIDIA
            switch (order[ii]) {
            case SW:
                // lower left
                i_new[indexval+ii] = ival*2;
                j_new[indexval+ii] = jval*2;
                break;
            case SE:
                // lower right
                i_new[indexval+ii] = ival*2+1;
                j_new[indexval+ii] = jval*2;
                break;
            case NW:
                // upper left
                i_new[indexval+ii] = ival*2;
                j_new[indexval+ii] = jval*2+1;
                break;
            case NE:
                // upper right
                i_new[indexval+ii] = ival*2+1;
                j_new[indexval+ii] = jval*2+1;
                break;
            }
#endif
         }
      } else if (ctype == LEFT_BOUNDARY) {
         i_new[indexval]   = ival*2+1;
         j_new[indexval]   = jval*2;
         i_new[indexval+1] = ival*2+1;
         j_new[indexval+1] = jval*2+1;
      } else if (ctype == RIGHT_BOUNDARY) {
         i_new[indexval]   = ival*2;
         j_new[indexval]   = jval*2;
         i_new[indexval+1] = ival*2;
         j_new[indexval+1] = jval*2+1;
      } else if (ctype == BOTTOM_BOUNDARY) {
         i_new[indexval]   = ival*2;
         j_new[indexval]   = jval*2+1;
         i_new[indexval+1] = ival*2+1;
         j_new[indexval+1] = jval*2+1;
      } else if (ctype == TOP_BOUNDARY) {
         i_new[indexval] = ival*2+1;
         j_new[indexval] = jval*2;
         i_new[indexval+1]   = ival*2;
         j_new[indexval+1]   = jval*2;
      }
   }

   indexoffset[giX] = indexval;
}


__kernel void rezone_neighbors_cl(
                 const int    isize,        // 0 
        __global const int   *mpot,         // 1
        __global const int   *index_offset, // 2
        __global const int   *nlft,         // 3
        __global const int   *nrht,         // 4
        __global const int   *nbot,         // 5
        __global const int   *ntop,         // 6
        __global const int   *celltype_new, // 7
        __global       int   *nlft_new,     // 8
        __global       int   *nrht_new,     // 9
        __global       int   *nbot_new,     // 10
        __global       int   *ntop_new)     // 11
{
   uint ic = get_global_id(0);

   if (ic >= isize) return;

   int nc = index_offset[ic];

   if (mpot[ic] == 0){
      nlft_new[nc] = (mpot[nlft[ic]] == 0) ? index_offset[nlft[ic]] : -1;
      nrht_new[nc] = (mpot[nrht[ic]] == 0) ? index_offset[nrht[ic]] : -1;
      nbot_new[nc] = (mpot[nbot[ic]] == 0) ? index_offset[nbot[ic]] : -1;
      ntop_new[nc] = (mpot[ntop[ic]] == 0) ? index_offset[ntop[ic]] : -1;
   } else if (mpot[ic] <= -2) {
      nlft_new[nc]  = -1;
      nrht_new[nc]  = -1;
      nbot_new[nc]  = -1;
      ntop_new[nc]  = -1;
   } else if (mpot[ic] > 0){
      nlft_new[nc]    = -1;
      nlft_new[nc+1]  = -1;
      nrht_new[nc]    = -1;
      nrht_new[nc+1]  = -1;
      nbot_new[nc]    = -1;
      nbot_new[nc+1]  = -1;
      ntop_new[nc]    = -1;
      ntop_new[nc+1]  = -1;
      if (celltype_new[nc] == REAL_CELL){
         nlft_new[nc+2]  = -1;
         nlft_new[nc+3]  = -1;
         nrht_new[nc+2]  = -1;
         nrht_new[nc+3]  = -1;
         nbot_new[nc+2]  = -1;
         nbot_new[nc+3]  = -1;
         ntop_new[nc+2]  = -1;
         ntop_new[nc+3]  = -1;
      }
   }
   switch(celltype_new[nc]){
   case LEFT_BOUNDARY:
      nlft_new[nc] = nc;
      break;
   case RIGHT_BOUNDARY:
      nrht_new[nc] = nc;
      break;
   case BOTTOM_BOUNDARY:
      nbot_new[nc] = nc;
      break;
   case TOP_BOUNDARY:
      ntop_new[nc] = nc;
      break;
   }
   if (mpot[ic] > 0){
      nc++;
      switch(celltype_new[nc]){
      case LEFT_BOUNDARY:
         nlft_new[nc] = nc;
         break;
      case RIGHT_BOUNDARY:
         nrht_new[nc] = nc;
         break;
      case BOTTOM_BOUNDARY:
         nbot_new[nc] = nc;
         break;
      case TOP_BOUNDARY:
         ntop_new[nc] = nc;
         break;
      }
   }
}

__kernel void rezone_one_float_cl(
                 const int    isize,        // 0
        __global const int   *i,            // 1
        __global const int   *j,            // 2
        __global const int   *nlft,         // 3
        __global const int   *nrht,         // 4
        __global const int   *nbot,         // 5
        __global const int   *ntop,         // 6
        __global const int   *celltype,     // 7
        __global const int   *mpot,         // 8   Array of mesh potential information.
        __global const uint  *indexoffset,  // 9
        __global const float *Var,          // 10
        __global       float *Var_new)      // 11
{
   uint giX = get_global_id(0);

   if (giX >= isize) return;

   float Varold = Var[giX];
   int mpotval = mpot[giX];
   int ctype = celltype[giX];
   int indexval = indexoffset[giX];

   if (mpotval == 0) {
      Var_new[indexval] = Varold;
   } else if (mpotval < 0) {
      if (is_lower_left(i[giX],j[giX]) ) {
         int nr = nrht[giX];
         int nt = ntop[giX];
         int nrt = nrht[nt];
         Var_new[indexval] = (Varold+Var[nr]+Var[nt]+Var[nrt])*0.25;
      }
      if (ctype != REAL_CELL && is_upper_right(i[giX],j[giX]) ) {
         int nl = nlft[giX];
         int nb = nbot[giX];
         int nlb = nlft[nb];
         Var_new[indexval] = (Varold+Var[nl]+Var[nb]+Var[nlb])*0.25;
      }
   } else if (mpotval > 0) {
      Var_new[indexval]   = Varold;
      Var_new[indexval+1] = Varold;
      if (ctype == REAL_CELL) {
         Var_new[indexval+2] = Varold;
         Var_new[indexval+3] = Varold;
      }
   }

}

#ifndef MINIMUM_PRECISION
__kernel void rezone_one_double_cl(
                 const int     isize,        // 0
        __global const int    *i,            // 1
        __global const int    *j,            // 2
        __global const int    *nlft,         // 3
        __global const int    *nrht,         // 4
        __global const int    *nbot,         // 5
        __global const int    *ntop,         // 6
        __global const int    *celltype,     // 7
        __global const int    *mpot,         // 8   Array of mesh potential information.
        __global const uint   *indexoffset,  // 9
        __global const double *Var,          // 10
        __global       double *Var_new)      // 11
{
   uint giX = get_global_id(0);

   if (giX >= isize) return;

   double Varold = Var[giX];
   int mpotval = mpot[giX];
   int ctype = celltype[giX];
   int indexval = indexoffset[giX];

   if (mpotval == 0) {
      Var_new[indexval] = Varold;
   } else if (mpotval < 0) {
      if (is_lower_left(i[giX],j[giX]) ) {
         int nr = nrht[giX];
         int nt = ntop[giX];
         int nrt = nrht[nt];
         Var_new[indexval] = (Varold+Var[nr]+Var[nt]+Var[nrt])*0.25;
      }
      if (ctype != REAL_CELL && is_upper_right(i[giX],j[giX]) ) {
         int nl = nlft[giX];
         int nb = nbot[giX];
         int nlb = nlft[nb];
         Var_new[indexval] = (Varold+Var[nl]+Var[nb]+Var[nlb])*0.25;
      }
   } else if (mpotval > 0) {
      Var_new[indexval]   = Varold;
      Var_new[indexval+1] = Varold;
      if (ctype == REAL_CELL) {
         Var_new[indexval+2] = Varold;
         Var_new[indexval+3] = Varold;
      }
   }

}
#endif

__kernel void copy_mpot_ghost_data_cl(
                          const int  ncells,        // 0
                          const int  nghost,        // 1
                 __global       int  *mpot,         // 2
                 __global       int  *mpot_add)     // 3
{
   const unsigned int giX  = get_global_id(0);

   if (giX >= nghost) return;

   mpot[ncells+giX] = mpot_add[giX];
}

__kernel void set_boundary_refinement(
               const int ncells,       // 0  Total number of cells
      __global const int  *nlft,       // 1
      __global const int  *nrht,       // 2
      __global const int  *nbot,       // 3
      __global const int  *ntop,       // 4
      __global const int  *i,          // 5
      __global const int  *j,          // 6
      __global const int  *celltype,   // 7
      __global       int  *mpot,       // 8  refinement flag
      __global       int2 *redscratch, // 9  reduction array
      __global       int  *ioffset,    // 10 prefix scan offset array
      __global       int2 *result,     // 11 new cell count
      __local        int2 *itile       // 12 Scratch to do reduction
                                 )
{
   const uint tiX = get_local_id(0);
   const uint giX = get_global_id(0);
   const uint gID = get_group_id(0);
   const uint ntX = get_local_size(0);
   const uint num_groups = get_num_groups(0);

   itile[tiX].s01 = 0;

   if (giX >= ncells) return;

   int ctype = celltype[giX];
   int mpotval = mpot[giX];
   int mpotval_save = mpotval;

   if (ctype < 0) {
      switch (ctype) {

      case LEFT_BOUNDARY:
         mpotval = mpot[nrht[giX]];
         break;
      case RIGHT_BOUNDARY:
         mpotval = mpot[nlft[giX]];
         break;
      case BOTTOM_BOUNDARY:
         mpotval = mpot[ntop[giX]];
         break;
      case TOP_BOUNDARY:
         mpotval = mpot[nbot[giX]];
         break;
      }
   }
   if (mpotval != mpotval_save) mpot[giX] = mpotval;

   if (mpotval < 0){
      int ival = i[giX];
      int jval = j[giX];
      if (ctype == REAL_CELL) {
         if (! is_lower_left(ival,jval) ) itile[tiX].s1 = 1;
      } else {
         if (! (is_upper_right(ival,jval) || is_lower_left(ival,jval) ) ) itile[tiX].s1 = 1;
      }
   } else if (mpotval > 0){
      if (ctype == REAL_CELL) {
         itile[tiX].s0 = 3;
      } else {
         itile[tiX].s0 = 1;
      }
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset = ntX >> 1; offset > 32; offset >>= 1) {
      if (tiX < offset) {
         itile[tiX].s01 += itile[tiX+offset].s01;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }

   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32) {
      itile[tiX].s01 += itile[tiX+32].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+16].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+8].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+4].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+2].s01;
      barrier(CLK_LOCAL_MEM_FENCE);
      itile[tiX].s01 += itile[tiX+1].s01;
   }

   if (tiX == 0) {
     redscratch[gID].s01 = itile[0].s01;
     if (gID == num_groups-1){
        ioffset[gID] = (ncells - (num_groups-1)*ntX)  + itile[0].s0 - itile[0].s1;
     } else {
        ioffset[gID] = ntX + itile[0].s0 - itile[0].s1;
     }
     result[0] = itile[0].s01;
   }
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

