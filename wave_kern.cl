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


__kernel void do_load_balance_cl_lower(
         __global       real *dev_H_new,
         __global       real *dev_U_new,
         __global       real *dev_V_new,
         __global       int  *dev_i_new,
         __global       int  *dev_j_new,
         __global       int  *dev_level_new,
         __global       int  *dev_celltype_new,
         __global const real *dev_H_lower,
         __global const real *dev_U_lower,
         __global const real *dev_V_lower,
         __global const int  *dev_i_lower,
         __global const int  *dev_j_lower,
         __global const int  *dev_level_lower,
         __global const int  *dev_celltype_lower,
         const int            lower_block_size)
{

   const unsigned int giX = get_global_id(0);

   if(giX >= lower_block_size) return;

   dev_H_new[giX]        = dev_H_lower[giX];
   dev_U_new[giX]        = dev_U_lower[giX];
   dev_V_new[giX]        = dev_V_lower[giX];
   dev_i_new[giX]        = dev_i_lower[giX];
   dev_j_new[giX]        = dev_j_lower[giX];
   dev_level_new[giX]    = dev_level_lower[giX];
   dev_celltype_new[giX] = dev_celltype_lower[giX];

}


__kernel void do_load_balance_cl_middle(
         __global       real *dev_H_new,
         __global       real *dev_U_new,
         __global       real *dev_V_new,
         __global       int  *dev_i_new,
         __global       int  *dev_j_new,
         __global       int  *dev_level_new,
         __global       int  *dev_celltype_new,
         __global const real *dev_H,
         __global const real *dev_U,
         __global const real *dev_V,
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

   dev_H_new[rgiX]        = dev_H[dgiX];
   dev_U_new[rgiX]        = dev_U[dgiX];
   dev_V_new[rgiX]        = dev_V[dgiX];
   dev_i_new[rgiX]        = dev_i[dgiX];
   dev_j_new[rgiX]        = dev_j[dgiX];
   dev_level_new[rgiX]    = dev_level[dgiX];
   dev_celltype_new[rgiX] = dev_celltype[dgiX];

}


__kernel void do_load_balance_cl_upper(
         __global       real *dev_H_new,
         __global       real *dev_U_new,
         __global       real *dev_V_new,
         __global       int  *dev_i_new,
         __global       int  *dev_j_new,
         __global       int  *dev_level_new,
         __global       int  *dev_celltype_new,
         __global const real *dev_H_upper,
         __global const real *dev_U_upper,
         __global const real *dev_V_upper,
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

   dev_H_new[rgiX]        = dev_H_upper[giX];
   dev_U_new[rgiX]        = dev_U_upper[giX];
   dev_V_new[rgiX]        = dev_V_upper[giX];
   dev_i_new[rgiX]        = dev_i_upper[giX];
   dev_j_new[rgiX]        = dev_j_upper[giX];
   dev_level_new[rgiX]    = dev_level_upper[giX];
   dev_celltype_new[rgiX] = dev_celltype_upper[giX];

}

#define hashval(j,i) hash[(j)*imaxsize+(i)]

__kernel void hash_init_cl(
                          const ulong isize,     // 0 
                 __global       int   *hash)     // 1 
{
   const ulong giX  = get_global_id(0);

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
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   sizes[0].s0 = max(sizes[0].s0-2*levtable[levmx],0);
   sizes[0].s1 = min(sizes[0].s1+2*levtable[levmx],(imax+1)*levtable[levmx]);
   sizes[0].s2 = max(sizes[0].s2-2*levtable[levmx],0);
   sizes[0].s3 = min(sizes[0].s3+2*levtable[levmx],(jmax+1)*levtable[levmx]);

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int ii = corners_i[giX];
   int jj = corners_j[giX];

   if (ii >= iminsize && ii < imaxsize &&
       jj >= jminsize && jj < jmaxsize) {
      hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)] = -99;
   }

}

__kernel void hash_setup_cl(
                          const int  isize,     // 0 
                          const int  levmx,     // 1 
                          const int  imax,      // 2 
                          const int  jmax,      // 3 
                          const int  imaxsize,  // 4 
                 __global const int  *levtable, // 5 
                 __global const int  *lev_ibeg, // 6
                 __global const int  *lev_iend, // 7
                 __global const int  *lev_jbeg, // 8
                 __global const int  *lev_jend, // 9
                 __global const int  *level,    // 10
                 __global const int  *i,        // 11
                 __global const int  *j,        // 12
                 __global       int  *hash)     // 13
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int lev = level[giX];
   int ii = i[giX];
   int jj = j[giX];

#define HASH_SETUP_NEW

#ifdef HASH_SETUP_NEW
   if (lev > 0 && (ii < lev_ibeg[lev] || ii > lev_iend[lev] || jj < lev_jbeg[lev] || jj > lev_jend[lev]) ) {
#endif
      int levdiff = levmx-lev;
   
      int iimin =  ii   *levtable[levdiff];
      int iimax = (ii+1)*levtable[levdiff];
      int jjmin =  jj   *levtable[levdiff];
      int jjmax = (jj+1)*levtable[levdiff];

      if      (ii < lev_ibeg[lev]) iimin = 0;                        // left boundary
      else if (ii > lev_iend[lev]) iimax = (imax+1)*levtable[levmx]; // right boundary
      else if (jj < lev_jbeg[lev]) jjmin = 0;                        // bottom boundary
      else if (jj > lev_jend[lev]) jjmax = (jmax+1)*levtable[levmx]; // top boundary 
   
      for (   int jjj = jjmin; jjj < jjmax; jjj++) {
         for (int iii = iimin; iii < iimax; iii++) {
            hashval(jjj, iii) = giX;
         }
      }
#ifdef HASH_SETUP_NEW
      return;
   }

   if(lev == levmx) {
      hashval(jj,ii) = giX;
      return;
   }

   int wid = levtable[levmx-lev];
   jj *= wid;
   ii *= wid;
   hashval(jj,ii) = giX;

   ii += wid/2;
   hashval(jj,ii) = giX;
   if(wid > 2) {
      ii = ii + wid/2 - 1;
      hashval(jj,ii) = giX;
      ii = ii - wid/2 + 1;
   }
   ii -= wid/2;
   jj += wid/2;
   hashval(jj,ii) = giX;
   ii = ii + wid - 1;
   hashval(jj,ii) = giX;

   if(wid > 2) {
      ii = ii - wid + 1;
      jj = jj + wid/2 - 1;
      hashval(jj,ii) = giX;
      ii += wid/2;
      hashval(jj,ii) = giX;
   }
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

   int lev = level[giX];

   if (i[giX] < lev_ibeg[lev]) { // left boundary
      for (int    jj = j[giX]*levtable[levmx-lev]-jminsize; jj < (j[giX]+1)*levtable[levmx-lev]-jminsize; jj++) {
         for (int ii = 0;                                   ii < (i[giX]+1)*levtable[levmx-lev]-iminsize; ii++) {
            hash[jj*(imaxsize-iminsize)+ii] = giX+noffset;
         }
      }
   } else if (i[giX] > lev_iend[lev]) { // right boundary
      for (int    jj = j[giX]*levtable[levmx-lev]-jminsize; jj < (j[giX]+1)*levtable[levmx-lev]-jminsize; jj++) {
         for (int ii = i[giX]*levtable[levmx-lev]-iminsize; ii < (imax+1)*levtable[levmx]-iminsize;       ii++) {
            hash[jj*(imaxsize-iminsize)+ii] = giX+noffset;
         }
      }
   } else if (j[giX] < lev_jbeg[lev]) { // bottom boundary
      for (int    jj = 0;                                   jj < (j[giX]+1)*levtable[levmx-lev]-jminsize; jj++) {
         for (int ii = i[giX]*levtable[levmx-lev]-iminsize; ii < (i[giX]+1)*levtable[levmx-lev]-iminsize; ii++) {
            hash[jj*(imaxsize-iminsize)+ii] = giX+noffset;
         }
      }
   } else if (j[giX] > lev_jend[lev]) { // top boundary
      for (int    jj = j[giX]*levtable[levmx-lev]-jminsize; jj < (jmax+1)*levtable[levmx]-jminsize;       jj++) {
         for (int ii = i[giX]*levtable[levmx-lev]-iminsize; ii < (i[giX]+1)*levtable[levmx-lev]-iminsize; ii++) {
            hash[jj*(imaxsize-iminsize)+ii] = giX+noffset;
         }
      }
   } else if (lev == levmx) {
      hash[(j[giX]-jminsize)*(imaxsize-iminsize)+(i[giX]-iminsize)] = giX+noffset;
   } else {
/* Original Hash Setup
      for (   int jj = j[giX]*levtable[levmx-lev]-jminsize; jj < (j[giX]+1)*levtable[levmx-lev]-jminsize; jj++) {
         for (int ii = i[giX]*levtable[levmx-lev]-iminsize; ii < (i[giX]+1)*levtable[levmx-lev]-iminsize; ii++) {
            hash[jj*(imaxsize-iminsize)+ii] = giX+noffset;
         }
      }
*/
/* Optimized Hash Setup */
      int wid = levtable[levmx-lev];
      int jj = j[giX]*wid - jminsize;
      int ii = i[giX]*wid - iminsize;
      #define hashv(jj,ii) ( hash[(jj)*(imaxsize-iminsize)+(ii)] )

      hashv(jj,ii) = giX + noffset;

      ii += wid/2;
      hashv(jj,ii) = giX + noffset;
      if(wid > 2) {
         ii = ii + wid/2 - 1;
         hashv(jj,ii) = giX + noffset;
         ii = ii - wid/2 + 1;
      }
      ii -= wid/2;
      jj += wid/2;
      hashv(jj,ii) = giX + noffset;
      ii = ii + wid - 1;
      hashv(jj,ii) = giX + noffset;

      if(wid > 2) {
         ii = ii - wid + 1;
         jj = jj + wid/2 - 1;
         hashv(jj,ii) = giX + noffset;
         ii += wid/2;
         hashv(jj,ii) = giX + noffset;
      }

   }

}

__kernel void hash_setup_border_cl(
                          const int  isize,         // 0 
                          const int  levmx,         // 1 
                          const int  noffset,       // 2 
                 __global const int4 *sizes,        // 3 
                 __global const int  *levtable,     // 4 
                 __global const int  *border_level, // 5 
                 __global const int  *border_i,     // 6 
                 __global const int  *border_j,     // 7 
                 __global const int  *border_num,   // 8 
                 __global       int  *hash)         // 9 
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   int jmaxsize = sizes[0].s3;

   int lev = border_level[giX];
   int ii  = border_i[giX];
   int jj  = border_j[giX];
   int num = border_num[giX];

   if (lev == levmx) {
      hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)] = num;
   } else {
      for (   int jjj = max(jj*levtable[levmx-lev]-jminsize,0); jjj < min((jj+1)*levtable[levmx-lev],jmaxsize)-jminsize; jjj++) {
         for (int iii = max(ii*levtable[levmx-lev]-iminsize,0); iii < min((ii+1)*levtable[levmx-lev],imaxsize)-iminsize; iii++) {
            hash[jjj*(imaxsize-iminsize)+iii] = num;
         }
      }
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

   tile[tiX].s0 =  9999;
   tile[tiX].s1 = -9999;
   tile[tiX].s2 =  9999;
   tile[tiX].s3 = -9999;

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

__kernel void finish_reduction_minmax4_cl(
                 const int    isize,        // 0
        __global       int4   *redscratch,  // 1
        __global       int4   *sizes,       // 2
        __local        int4   *tile)        // 3
{  
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX].s0 =  9999;
   tile[tiX].s1 = -9999;
   tile[tiX].s2 =  9999;
   tile[tiX].s3 = -9999;

   if (tiX < isize) tile[tiX].s0123 = redscratch[giX].s0123;

   for (giX += ntX; giX < isize; giX += ntX) {
     if (redscratch[giX].s0 < tile[tiX].s0) tile[tiX].s0 = redscratch[giX].s0;
     if (redscratch[giX].s1 > tile[tiX].s1) tile[tiX].s1 = redscratch[giX].s1;
     if (redscratch[giX].s2 < tile[tiX].s2) tile[tiX].s2 = redscratch[giX].s2;
     if (redscratch[giX].s3 > tile[tiX].s3) tile[tiX].s3 = redscratch[giX].s3;
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_minmax_within_tile4(tile);

   if (tiX == 0) {
     redscratch[0].s0123 = tile[0].s0123;
     sizes[0].s0123 = tile[0].s0123;
   }
}

__kernel void calc_neighbors_cl(
                          const int  isize,     // 0 
                          const int  levmx,     // 1 
                          const int  imax,      // 2 
                          const int  jmax,      // 3 
                          const int  imaxsize,  // 2 
                          const int  jmaxsize,  // 3 
                 __global const int  *levtable, // 4 
                 __global const int  *level,    // 5 
                 __global const int  *i,        // 6 
                 __global const int  *j,        // 7 
                 __global       int  *nlft,     // 8 
                 __global       int  *nrht,     // 9 
                 __global       int  *nbot,     // 10
                 __global       int  *ntop,     // 11
                 __global const int  *hash)     // 12
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int ii, jj, iii, jjj, lev, levmult;
   int nlftval, nrhtval, nbotval, ntopval;

   ii = i[giX];
   jj = j[giX];
   lev = level[giX];
   levmult = levtable[levmx-lev];
   nlftval = hashval(      jj   *levmult               , max(  ii   *levmult-1, 0         ));
   nrhtval = hashval(      jj   *levmult               , min( (ii+1)*levmult,   imaxsize-1));
   nbotval = hashval(max(  jj   *levmult-1, 0)         ,       ii   *levmult               );
   ntopval = hashval(min( (jj+1)*levmult,   jmaxsize-1),       ii   *levmult               );

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
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int ii, jj, iii, jjj, lev, levmult;
   int nlftval;
   int nrhtval;
   int nbotval;
   int ntopval;

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   //int jmaxsize = sizes[0].s3;

   int imaxcalc = (imax+1)*levtable[levmx];
   int jmaxcalc = (jmax+1)*levtable[levmx];

   ii = i[giX];
   jj = j[giX];
   lev = level[giX];
   levmult = levtable[levmx-lev];
   nlftval = hash[ ( (      jj   *levmult               )-jminsize) *(imaxsize-iminsize) + ( (max(  ii   *levmult-1, 0         ))-iminsize)];
   nrhtval = hash[ ( (      jj   *levmult               )-jminsize) *(imaxsize-iminsize) + ( (min( (ii+1)*levmult,   imaxcalc-1))-iminsize)];
   nbotval = hash[ ( (max(  jj   *levmult-1, 0)         )-jminsize) *(imaxsize-iminsize) + ( (      ii   *levmult               )-iminsize)];
   ntopval = hash[ ( (min( (jj+1)*levmult,   jmaxcalc-1))-jminsize) *(imaxsize-iminsize) + ( (      ii   *levmult               )-iminsize)];

   // Handles the four boundary corners
   if (nlftval == -99){
      iii = ii*levmult;
      jjj = jj*levmult;
      nlftval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (nrhtval == -99){
      iii = ii*levmult;
      jjj = jj*levmult;
      nrhtval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (nbotval == -99) {
      iii = ii*levmult;
      jjj = jj*levmult;
      nbotval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (ntopval == -99) {
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
                          const int isize,      // 0
                          const int noffset,    // 1
                 __global const int *nlft,      // 2
                 __global const int *nrht,      // 3
                 __global const int *nbot,      // 4
                 __global const int *ntop,      // 5
                 __global const int *level,     // 6
                 __global const uint *border_cell_in, //7
                 __global       uint *border_cell_out //8
                                  )
{
   const unsigned int giX = get_global_id(0);

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

   border_cell_out[giX] = border_cell;
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
          hash[(min((jj+1)*levmult -1,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(max(ii*levmult-1, 0)-iminsize)] >= 0 ) {
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
          hash[((max(jj*levmult-1, 0) )-jminsize)*(imaxsize-iminsize)+(min((ii+1)*levmult -1,imaxcalc-1)-iminsize)] >= 0 ) {
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
                          const int  isize,                   // 0 
                          const int  ncells,                  // 1 
                          const int  noffset,                 // 2 
                          const int  levmx,                   // 3 
                          const int  imax,                    // 4 
                          const int  jmax,                    // 5 
                 __global const int4 *sizes,                  // 6 
                 __global const int  *levtable,               // 7 
                 __global const int  *border_cell_i,          // 8
                 __global const int  *border_cell_j,          // 9
                 __global const int  *border_cell_level,      // 10
                 __global const int  *border_cell_needed,     // 11
                 __global       int  *border_cell_needed_out, // 12
                 __global const int  *hash)                   // 13
{
   const uint giX = get_global_id(0);

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
   border_cell_needed_out[giX] = iborder;
}

__kernel void calc_layer2_sethash_cl (
                          const int  isize,               // 0 
                          const int  ncells,              // 1 
                          const int  noffset,             // 2 
                          const int  levmx,               // 3 
                 __global const int4 *sizes,              // 4 
                 __global const int  *levtable,           // 5 
                 __global const int  *border_cell_i,      // 6
                 __global const int  *border_cell_j,      // 7
                 __global const int  *border_cell_level,  // 8
                 __global const int  *border_cell_num,    // 9
                 __global       int  *border_cell_needed, // 10
                 __global       int  *hash)               // 11
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
      int cell_number = border_cell_num[giX];
      int levmult = levtable[levmx-lev];

      if (lev == levmx) {
         hash[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)] = cell_number;
      } else {
         for (int    j = max(jj*levmult-jminsize,0); j < min((jj+1)*levmult,jmaxsize)-jminsize; j++) {
            for (int i = max(ii*levmult-iminsize,0); i < min((ii+1)*levmult,imaxsize)-iminsize; i++) {
               hash[j*(imaxsize-iminsize)+i] = cell_number;
            }
         }
      }
   }
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
   nlft[ncout]  = -98;
   nrht[ncout]  = -98;
   nbot[ncout]  = -98;
   ntop[ncout]  = -98;
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

   if (nlftval < 0){
      if (max(ii*levmult-1, 0)-iminsize >= 0) {  // Test for cell to left
         nlft[giX] = hash[((      jj   *levmult               )-jminsize)*(imaxsize-iminsize)+((max(  ii   *levmult-1, 0         ))-iminsize)];
      }
   }
   if (nrhtval < 0){
      if (min( (ii+1)*levmult, imaxcalc-1) < imaxsize){ // Test for cell to right
         nrht[giX] = hash[((      jj   *levmult               )-jminsize)*(imaxsize-iminsize)+((min( (ii+1)*levmult,   imaxcalc-1))-iminsize)];
      }
   }
   if (nbotval < 0){
      if (max(jj*levmult-1, 0)-jminsize >= 0){ // Test for cell to bottom
         nbot[giX] = hash[((max(  jj   *levmult-1, 0)         )-jminsize)*(imaxsize-iminsize)+((      ii   *levmult               )-iminsize)];
      }
   }
   if (ntopval < 0) {
      if (min( (jj+1)*levmult, jmaxcalc-1) < jmaxsize){ // Test for cell to top
         ntop[giX] = hash[((min( (jj+1)*levmult,   jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((      ii   *levmult               )-iminsize)];
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
  
   if (nlftval == -99){
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
   if (nrhtval == -99){
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
   if (nbotval == -99) {
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
   if (ntopval == -99) {
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

   if (nlftval >= noffset && nlftval < noffset+ncells) {
      nlftval -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (nlftval ==indices_needed[ig]) {nlftval = ig+ncells; break;}
      }
   }
   if (nrhtval >= noffset && nrhtval < noffset+ncells) {
      nrhtval -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (nrhtval==indices_needed[ig]) {nrhtval = ig+ncells; break;}
      }
   }
   if (nbotval >= noffset && nbotval < noffset+ncells) {
      nbotval -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (nbotval==indices_needed[ig]) {nbotval = ig+ncells; break;}
      }
   }
   if (ntopval >= noffset && ntopval < noffset+ncells) {
      ntopval -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (ntopval==indices_needed[ig]) {ntopval = ig+ncells; break;}
      }
   }

   nlft[giX] = nlftval;
   nrht[giX] = nrhtval;
   nbot[giX] = nbotval;
   ntop[giX] = ntopval;
}

__kernel void calc_neighbors_local2_cl(
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
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   int ii, jj, iii, jjj, lev, levmult;
   int nlftval = nlft[giX];
   int nrhtval = nrht[giX];
   int nbotval = nbot[giX];
   int ntopval = ntop[giX];

   int iminsize = sizes[0].s0;
   int imaxsize = sizes[0].s1;
   int jminsize = sizes[0].s2;
   //int jmaxsize = sizes[0].s3;

   int imaxcalc = (imax+1)*levtable[levmx];
   int jmaxcalc = (jmax+1)*levtable[levmx];

   ii = i[giX];
   jj = j[giX];
   lev = level[giX];
   levmult = levtable[levmx-lev];
   if (nlftval == -1) nlftval = hash[ ( (      jj   *levmult               )-jminsize) *(imaxsize-iminsize) + ( (max(  ii   *levmult-1, 0         ))-iminsize)];
   if (nrhtval == -1) nrhtval = hash[ ( (      jj   *levmult               )-jminsize) *(imaxsize-iminsize) + ( (min( (ii+1)*levmult,   imaxcalc-1))-iminsize)];
   if (nbotval == -1) nbotval = hash[ ( (max(  jj   *levmult-1, 0)         )-jminsize) *(imaxsize-iminsize) + ( (      ii   *levmult               )-iminsize)];
   if (ntopval == -1) ntopval = hash[ ( (min( (jj+1)*levmult,   jmaxcalc-1))-jminsize) *(imaxsize-iminsize) + ( (      ii   *levmult               )-iminsize)];

   // Handles the four boundary corners
   if (nlftval == -99){
      iii = ii*levmult;
      jjj = jj*levmult;
      nlftval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (nrhtval == -99){
      iii = ii*levmult;
      jjj = jj*levmult;
      nrhtval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (nbotval == -99) {
      iii = ii*levmult;
      jjj = jj*levmult;
      nbotval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }
   if (ntopval == -99) {
      iii = ii*levmult;
      jjj = jj*levmult;
      ntopval = hash[(jjj-jminsize)*(imaxsize-iminsize)+(iii-iminsize)];
   }

   nlft[giX] = nlftval;
   nrht[giX] = nrhtval;
   nbot[giX] = nbotval;
   ntop[giX] = ntopval;
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

__kernel void copy_ghost_data_cl(
                          const int  ncells,        // 0 
                          const int  nghost,        // 1 
                 __global       int  *celltype,     // 2
                 __global const int  *celltype_add) // 3
{
   const unsigned int giX  = get_global_id(0);

   if (giX >= nghost) return;

   celltype[ncells+giX] = celltype_add[giX];
}

__kernel void adjust_neighbors_cl(
                          const int  ncells_ghost,    // 0 
                          const int  nghost,          // 1 
                          const int  noffset,         // 2 
                 __global const int  *indices_needed, // 3 
                 __global       int  *nlft,           // 4
                 __global       int  *nrht,           // 5
                 __global       int  *nbot,           // 6
                 __global       int  *ntop)           // 7
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= ncells_ghost) return;

   int ncells = ncells_ghost - nghost;

   if (nlft[giX] >= noffset && nlft[giX] <noffset+ncells){
      nlft[giX] -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (nlft[giX]==indices_needed[ig]) {nlft[giX] = ig+ncells; break;}
      }
   }

   if (nrht[giX] >= noffset && nrht[giX] <noffset+ncells){
      nrht[giX] -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (nrht[giX]==indices_needed[ig]) {nrht[giX] = ig+ncells; break;}
      }
   }

   if (nbot[giX] >= noffset && nbot[giX] <noffset+ncells){
      nbot[giX] -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (nbot[giX]==indices_needed[ig]) {nbot[giX] = ig+ncells; break;}
      }
   }

   if (ntop[giX] >= noffset && ntop[giX] <noffset+ncells){
      ntop[giX] -= noffset;
   } else {
      for (int ig=0; ig<nghost; ig++){
         if (ntop[giX]==indices_needed[ig]) {ntop[giX] = ig+ncells; break;}
      }
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

void reduction_sum_int_within_tile(__local  int  *tile);

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

