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

#define HASH_SETUP_OPT_LEVEL 2

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

void reduction_max_within_tile1(__local  real  *tile);
void reduction_min_within_tile1(__local  real  *tile);
void reduction_max_within_tile2(__local  real2  *tile);
void reduction_minmax_within_tile4(__local  int4  *tile);

__kernel void set_timestep_cl(
                 const int    ncells,    // 0  Total number of cells.
                 const real   sigma,     // 1
        __global const real  *H_in,      // 2  
        __global const real  *U_in,      // 3  
        __global const real  *V_in,      // 4  
        __global const int   *level,     // 5  Array of level information.
        __global const int   *celltype,  // 6 
        __global const real  *lev_dx,    // 7  
        __global const real  *lev_dy,    // 8  
        __global       real  *redscratch,// 9 
        __global       real  *deltaT,    // 10
        __local        real  *tile)      // 11
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);
    
    tile[tiX] = 1000.0;

    if (giX >= ncells) return;
    
    const unsigned int group_id = get_group_id(0);
    const unsigned int ntX  = get_local_size(0);

    //  Set physical constants.
    const real   g     = GRAVITATIONAL_CONSTANT;   // gravitational constant
    
    //--MEMORY MANAGEMENT-------------------------------------------------------
    //  Set values for the main cell.
    real H       = H_in[giX];
    real U       = U_in[giX];
    real V       = V_in[giX];
    int lev      = level[giX];
    int type     = celltype[giX];
    
    //--CALCULATIONS------------------------------------------------------------
    if (type == REAL_CELL){
      real wavespeed = sqrt(g * H);
      real xspeed  = (fabs(U) + wavespeed)/lev_dx[lev];
      real yspeed  = (fabs(V) + wavespeed)/lev_dy[lev];
      tile[tiX] = sigma/(xspeed+yspeed);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    reduction_min_within_tile1(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      redscratch[group_id] = tile[0];
      (*deltaT)            = tile[0];
    }
}

/* finish_reduction */

__kernel void finish_reduction_min_cl(
                 const int   isize,
        __global       real  *redscratch,
        __global       real  *deltaT,
        __local        real  *tile)
{       
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX] = 1.0e20;

   if (tiX < isize) tile[tiX] = redscratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
     if (redscratch[giX] < tile[tiX]) tile[tiX] = redscratch[giX];
   }
 
   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_min_within_tile1(tile);

   if (tiX == 0) {
     (*deltaT) = tile[0];
   }
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

__kernel void rezone_all_cl(
                 const int    isize,        // 0 
                 const int    stencil,      // 1 
                 const int    levmx,        // 2
        __global const int   *mpot,         // 3   Array of mesh potential information.
        __global const int   *level,        // 4
        __global const int   *i,            // 5
        __global const int   *j,            // 6
        __global const int   *celltype,     // 7
        __global const real  *H,            // 8
        __global const real  *U,            // 9
        __global const real  *V,            // 10
        __global       int   *level_new,    // 11
        __global       int   *i_new,        // 12
        __global       int   *j_new,        // 13
        __global       int   *celltype_new, // 14
        __global       real  *H_new,        // 15
        __global       real  *U_new,        // 16
        __global       real  *V_new,        // 17
        __global const int   *ioffset,      // 18   
        __global const real  *lev_dx,       // 19
        __global const real  *lev_dy,       // 20
        __global const int   *levtable,     // 21
        __global const int   *ijadd,        // 22
        __local        uint  *itile,        // 23
        __local        real2 *tile)         // 24
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

   if (celltype[giX] == REAL_CELL){
      indexval = mpot[giX] ? 4 : 1;
   } else {
      indexval = mpot[giX] ? 2 : 1;
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
   int mpotval = mpot[giX];
   if (mpotval > 0) lev++;
   real Hval = H[giX];
   real Uval = U[giX];
   real Vval = V[giX];
   int ctype = celltype[giX];
   int ival = i[giX];
   int jval = j[giX];

   if (mpotval == 0) {
      i_new[indexval] = ival;
      j_new[indexval] = jval;
      level_new[indexval] = lev;
      celltype_new[indexval] = ctype;
      H_new[indexval] = Hval;
      U_new[indexval] = Uval;
      V_new[indexval] = Vval;
   } else if (mpotval < 0) {
   } else if (mpotval > 0) {
      level_new[indexval] = lev;
      level_new[indexval+1] = lev;
      celltype_new[indexval] = ctype;
      celltype_new[indexval+1] = ctype;
      H_new[indexval] = Hval;
      H_new[indexval+1] = Hval;
      U_new[indexval] = Uval;
      U_new[indexval+1] = Uval;
      V_new[indexval] = Vval;
      V_new[indexval+1] = Vval;
      if (ctype == REAL_CELL) {
         level_new[indexval+2] = lev;
         level_new[indexval+3] = lev;
         celltype_new[indexval+2] = ctype;
         celltype_new[indexval+3] = ctype;
         H_new[indexval+2] = Hval;
         H_new[indexval+3] = Hval;
         U_new[indexval+2] = Uval;
         U_new[indexval+3] = Uval;
         V_new[indexval+2] = Vval;
         V_new[indexval+3] = Vval;
          
         int order[4] = {SW, SE, NW, NE};
         if (stencil) {
#ifdef __OLD_STENCIL__
            real nx[3], ny[3];
            int ifirst = ijadd[0];
            int iend   = ijadd[1];
            int jfirst = ijadd[2];
            int jend   = ijadd[3];
            int level_first  = ijadd[4];
            int level_end    = ijadd[5];

//             x[ic]  = xmin + lev_deltax[level[ic]] * (real)(i[ic] - ibase);

            if (giX != 0) {
               nx[0] = lev_dx[level[giX-1]] * (real)i[giX-1];
               ny[0] = lev_dy[level[giX-1]] * (real)j[giX-1];
            } else {
               nx[0] = lev_dx[level_first] * ifirst;
               ny[0] = lev_dy[level_first] * jfirst;
            }

            nx[1] = lev_dx[level[giX  ]] * (real)i[giX  ];
            ny[1] = lev_dy[level[giX  ]] * (real)j[giX  ];
            if (giX != isize-1) {
               nx[2] = lev_dx[level[giX+1]] * (real)i[giX+1];
               ny[2] = lev_dy[level[giX+1]] * (real)j[giX+1];
            } else {
               nx[2] = lev_dx[level_end] * iend;
               ny[2] = lev_dy[level_end] * jend;
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
                  //printf("Nonlocal case for the stencil.\n");
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

}


#define hashval(j,i) hash[(j)*imaxsize+(i)]
#define hashval_local(j,i) hash[(j)*(imaxsize-iminsize)+(i)]

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

void reduction_max_within_tile1(__local  real  *tile) 
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

    for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        if (tile[tiX+offset] > tile[tiX]) tile[tiX] = tile[tiX+offset];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tiX < 32){
      if (tile[tiX+32] > tile[tiX]) tile[tiX] = tile[tiX+32];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+16] > tile[tiX]) tile[tiX] = tile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+8] > tile[tiX]) tile[tiX] = tile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+4] > tile[tiX]) tile[tiX] = tile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+2] > tile[tiX]) tile[tiX] = tile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+1] > tile[tiX]) tile[tiX] = tile[tiX+1];
    }

}

void reduction_min_within_tile1(__local  real  *tile) 
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

    for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        if (tile[tiX+offset] < tile[tiX]) tile[tiX] = tile[tiX+offset];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tiX < 32){
      if (tile[tiX+32] < tile[tiX]) tile[tiX] = tile[tiX+32];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+16] < tile[tiX]) tile[tiX] = tile[tiX+16];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+8] < tile[tiX]) tile[tiX] = tile[tiX+8];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+4] < tile[tiX]) tile[tiX] = tile[tiX+4];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+2] < tile[tiX]) tile[tiX] = tile[tiX+2];
      barrier(CLK_LOCAL_MEM_FENCE);
      if (tile[tiX+1] < tile[tiX]) tile[tiX] = tile[tiX+1];
    }

}

void reduction_max_within_tile2(__local  real2  *tile) 
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

    for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        if (tile[tiX+offset].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+offset].s0;
        if (tile[tiX+offset].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+offset].s1;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tiX < 32){
      if (tile[tiX+32].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+32].s0;
      if (tile[tiX+16].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+16].s0;
      if (tile[tiX+8].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+8].s0;
      if (tile[tiX+4].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+4].s0;
      if (tile[tiX+2].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+2].s0;
      if (tile[tiX+1].s0 > tile[tiX].s0) tile[tiX].s0 = tile[tiX+1].s0;

      if (tile[tiX+32].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+32].s1;
      if (tile[tiX+16].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+16].s1;
      if (tile[tiX+8].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+8].s1;
      if (tile[tiX+4].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+4].s1;
      if (tile[tiX+2].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+2].s1;
      if (tile[tiX+1].s1 > tile[tiX].s1) tile[tiX].s1 = tile[tiX+1].s1;
    }

}

void reduction_sum_within_tile(__local  real *tile);

__kernel void reduce_sum_mass_stage1of2_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const real  *array,     // 1
        __global const int   *level,     // 2
        __global const real  *levdx,     // 3
        __global const real  *levdy,     // 4
        __global const int   *celltype,  // 5
        __global       real  *mass_sum,  // 6
        __global       real  *scratch,   // 7
        __local        real  *tile)      // 8
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);

    const unsigned int group_id = get_group_id(0);

    tile[tiX] = ZERO;
    if (giX < isize && celltype[giX] == REAL_CELL) {
      tile[tiX] = array[giX]*levdx[level[giX]]*levdy[level[giX]];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_sum_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      scratch[group_id] = tile[0];
      (*mass_sum) = tile[0];
    }
}

void reduction_epsum_within_tile(__local  real2 *tile);

__kernel void reduce_epsum_mass_stage1of2_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const real  *array,     // 1
        __global const int   *level,     // 2
        __global const real  *levdx,     // 3
        __global const real  *levdy,     // 4
        __global const int   *celltype,  // 5
        __global       real2  *mass_sum, // 6
        __global       real2 *scratch,   // 7
        __local        real2 *tile)      // 8
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);

    const unsigned int group_id = get_group_id(0);

    // Going to do a Kahan sum -- an enhanced precision sum, so load
    // tile s0 with data and s1 with ZERO. For sizes greater than 
    // giX, need to initialize to ZERO for tile reduction below
    tile[tiX].s01 = ZERO;
    if (giX < isize && celltype[giX] == REAL_CELL) {
      tile[tiX].s0 = array[giX]*levdx[level[giX]]*levdy[level[giX]];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_epsum_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    // s0 contains the sum and s1 has the correction term
    if (tiX == 0){
      scratch[group_id].s01 = tile[0].s01;
      (*mass_sum).s01 = tile[0].s01;
    }
}

__kernel void reduce_sum_mass_stage2of2_cl(
                 const int    isize,
        __global       real   *mass_sum,
        __global       real   *scratch,
        __local        real   *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX] = 0.0;

   // load the sum from scratch
   if (tiX < isize) tile[tiX] = scratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
      tile[tiX] += scratch[giX];
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_sum_within_tile(tile);

   if (tiX == 0) {
     (*mass_sum) = tile[0];
   }
}
__kernel void reduce_epsum_mass_stage2of2_cl(
                 const int    isize,
        __global       real2  *mass_sum,
        __global       real2  *scratch,
        __local        real2  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   real corrected_next_term, new_sum;

   int giX = tiX;

   tile[tiX].s01 = 0.0;

   // load both the sum s0 and the correction term s1 from scratch
   if (tiX < isize) tile[tiX].s01 = scratch[giX].s01;

   for (giX += ntX; giX < isize; giX += ntX) {
      // Kahan sum
      corrected_next_term = scratch[giX].s0 + (scratch[giX].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_epsum_within_tile(tile);

   if (tiX == 0) {
     // s0 contains the sum and s1 has the correction term
     (*mass_sum).s01 = tile[0].s01;
   }
}

void reduction_sum_within_tile(__local  real  *tile)
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

void reduction_epsum_within_tile(__local  real2  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   real corrected_next_term, new_sum;

    for (int offset=ntX>>1; offset > 32; offset >>= 1){
      if (tiX < offset){
        // Kahan sum
        corrected_next_term = tile[tiX+offset].s0 + (tile[tiX+offset].s1 +tile[tiX].s1);
        new_sum = tile[tiX].s0 + corrected_next_term;
        tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
        tile[tiX].s0 = new_sum;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tiX < 32){
      // Kahan sum -- unrolled 
      corrected_next_term = tile[tiX+32].s0 + (tile[tiX+32].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */

      corrected_next_term = tile[tiX+16].s0 + (tile[tiX+16].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */

      corrected_next_term = tile[tiX+8].s0 + (tile[tiX+8].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */

      corrected_next_term = tile[tiX+4].s0 + (tile[tiX+4].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */

      corrected_next_term = tile[tiX+2].s0 + (tile[tiX+2].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
      barrier(CLK_LOCAL_MEM_FENCE);         /* Fix for Cuda 4.1 */

      corrected_next_term = tile[tiX+1].s0 + (tile[tiX+1].s1 +tile[tiX].s1);
      new_sum = tile[tiX].s0 + corrected_next_term;
      tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
      tile[tiX].s0 = new_sum;
    }

}

