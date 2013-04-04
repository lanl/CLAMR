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
#define ZERO 0.0
#define MINUS_INFINITY -1.0e30
#define PLUS_INFINITY 1.0e30
#else
typedef float   real;
#define ZERO 0.0f
#define MINUS_INFINITY -1.0e30
#define PLUS_INFINITY 1.0e30
#endif
#endif

/*
		   MPI_BAND                MPI_Op
		   MPI_BOR                 MPI_Op
		   MPI_BXOR                MPI_Op
		   MPI_LAND                MPI_Op
		   MPI_LOR                 MPI_Op
		   MPI_LXOR                MPI_Op
		   MPI_MAX                 MPI_Op
		   MPI_MAXLOC              MPI_Op
		   MPI_MIN                 MPI_Op
		   MPI_MINLOC              MPI_Op
		   MPI_OP_NULL             MPI_Op
		   MPI_PROD                MPI_Op
		   MPI_SUM                 MPI_Op
*/

__kernel void reduce_sum_cl(
                 const int    isize,    // 0   Array length.
        __global       real  *array,    // 1   Array to be reduced.
        __global       real  *result,   // 2   Final result of operation.
        __local        real  *tile)     // 3   
{  const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   
   int giX = tiX;
   
   tile[tiX] = array[giX];
   
   for (giX += ntX; giX < isize; giX += ntX)
   {  tile[tiX] += array[giX]; }
   barrier(CLK_LOCAL_MEM_FENCE);
   
   for (int offset = ntX >> 1; offset > 32; offset >>= 1)
   {  if (tiX < offset)
      {  tile[tiX] += tile[tiX+offset]; }
      barrier(CLK_LOCAL_MEM_FENCE); }
   
   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  tile[tiX] += tile[tiX+32];
      tile[tiX] += tile[tiX+16];
      tile[tiX] += tile[tiX+8];
      tile[tiX] += tile[tiX+4];
      tile[tiX] += tile[tiX+2];
      tile[tiX] += tile[tiX+1]; }
   
   if (tiX == 0)
   {  result[0] = tile[0]; }
}

__kernel void reduce_product_cl(
                 const int    isize,    // 0   Array length.
        __global       real  *array,    // 1   Array to be reduced.
        __global       real  *result,   // 2   Final result of operation.
        __local        real  *tile)     // 3   
{  const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   
   int giX = tiX;
   
   tile[tiX] = array[giX];
   
   for (giX += ntX; giX < isize; giX += ntX)
   {  tile[tiX] *= array[giX]; }
   barrier(CLK_LOCAL_MEM_FENCE);
   
   for (int offset = ntX >> 1; offset > 32; offset >>= 1)
   {  if (tiX < offset)
      {  tile[tiX] *= tile[tiX+offset]; }
      barrier(CLK_LOCAL_MEM_FENCE); }
   
   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  tile[tiX] *= tile[tiX+32];
      tile[tiX] *= tile[tiX+16];
      tile[tiX] *= tile[tiX+8];
      tile[tiX] *= tile[tiX+4];
      tile[tiX] *= tile[tiX+2];
      tile[tiX] *= tile[tiX+1]; }
   
   if (tiX == 0)
   {  result[0] = tile[0]; }
}

__kernel void reduce_max_cl(
                 const int    isize,    // 0   Array length.
        __global       real  *array,    // 1   Array to be reduced.
        __global       real  *result,   // 2   Final result of operation.
        __local        real  *tile)     // 3   
{  const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   
   int giX = tiX;
   
   tile[tiX] = array[giX];
   
   for (giX += ntX; giX < isize; giX += ntX)
   {  if (array[giX] > tile[tiX]) tile[tiX] = array[giX]; }
   barrier(CLK_LOCAL_MEM_FENCE);
   
   for (int offset = ntX >> 1; offset > 32; offset >>= 1)
   {  if (tiX < offset)
      {  if (array[giX+offset] > tile[tiX]) tile[tiX] = tile[tiX+offset]; }
      barrier(CLK_LOCAL_MEM_FENCE); }
   
   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  if (array[giX+32] > tile[tiX]) tile[tiX] = tile[tiX+32];
      if (array[giX+16] > tile[tiX]) tile[tiX] = tile[tiX+16];
      if (array[giX+8] > tile[tiX]) tile[tiX] = tile[tiX+8];
      if (array[giX+4] > tile[tiX]) tile[tiX] = tile[tiX+4];
      if (array[giX+2] > tile[tiX]) tile[tiX] = tile[tiX+2];
      if (array[giX+1] > tile[tiX]) tile[tiX] = tile[tiX+1]; }
   
   if (tiX == 0)
   {  result[0] = tile[0]; }
}

//XXX:  doesn't take abs() into account
__kernel void reduce_min_cl(
                 const int    isize,    // 0   Array length.
        __global       real  *array,    // 1   Array to be reduced.
        __global       real  *result,   // 2   Final result of operation.
        __local        real  *tile)     // 3   
{  const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   
   int giX = tiX;
   
   tile[tiX] = array[giX];
   
   for (giX += ntX; giX < isize; giX += ntX)
   {  if (array[giX] < tile[tiX]) tile[tiX] = array[giX]; }
   barrier(CLK_LOCAL_MEM_FENCE);
   
   for (int offset = ntX >> 1; offset > 32; offset >>= 1)
   {  if (tiX < offset)
      {  if (array[giX+offset] < tile[tiX]) tile[tiX] = tile[tiX+offset]; }
      barrier(CLK_LOCAL_MEM_FENCE); }
   
   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  if (array[giX+32] < tile[tiX]) tile[tiX] = tile[tiX+32];
      if (array[giX+16] < tile[tiX]) tile[tiX] = tile[tiX+16];
      if (array[giX+8] < tile[tiX]) tile[tiX] = tile[tiX+8];
      if (array[giX+4] < tile[tiX]) tile[tiX] = tile[tiX+4];
      if (array[giX+2] < tile[tiX]) tile[tiX] = tile[tiX+2];
      if (array[giX+1] < tile[tiX]) tile[tiX] = tile[tiX+1]; }
   
   if (tiX == 0)
   {  result[0] = tile[0]; }
}
/* //XXX:will need to be put in reduce.h, reduce.c as well.
__kernel void reduce_maxloc_cl(
                 const int    isize,    // 0   Array length.
        __global       real  *array,    // 1   Array to be reduced.
        __global       real  *result,   // 2   Final result of operation.
        __local        real  *tile)     // 3   
{  const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   
   int giX = tiX;
   
   tile[tiX] = array[giX]; //XXX:IC?
   
   for (giX += ntX; giX < isize; giX += ntX)
   {  if (array[giX] > tile[tiX]) tile[tiX] = giX; }
   barrier(CLK_LOCAL_MEM_FENCE);
   
   for (int offset = ntX >> 1; offset > 32; offset >>= 1)
   {  if (tiX < offset)
      {  if (array[giX+offset] > tile[tiX]) tile[tiX] = tiX+offset; }
      barrier(CLK_LOCAL_MEM_FENCE); }
   
   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  if (array[giX+32] > tile[tiX]) tile[tiX] = tiX+32;
      if (array[giX+16] > tile[tiX]) tile[tiX] = tiX+16;
      if (array[giX+8] > tile[tiX]) tile[tiX] = tiX+8;
      if (array[giX+4] > tile[tiX]) tile[tiX] = tiX+4;
      if (array[giX+2] > tile[tiX]) tile[tiX] = tiX+2;
      if (array[giX+1] > tile[tiX]) tile[tiX] = tiX+1; }
   
   if (tiX == 0)
   {  result[0] = tile[0]; }
}

//XXX:  doesn't take abs() into account
__kernel void reduce_minloc_cl(
                 const int    isize,    // 0   Array length.
        __global       real  *array,    // 1   Array to be reduced.
        __global       real  *result,   // 2   Final result of operation.
        __local        real  *tile)     // 3   
{  const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   
   int giX = tiX;
   
   tile[tiX] = array[giX]; //XXX:IC?
   
   for (giX += ntX; giX < isize; giX += ntX)
   {  if (array[giX] < tile[tiX]) tile[tiX] = giX; }
   barrier(CLK_LOCAL_MEM_FENCE);
   
   for (int offset = ntX >> 1; offset > 32; offset >>= 1)
   {  if (tiX < offset)
      {  if (array[giX+offset] < tile[tiX]) tile[tiX] = tiX+offset; }
      barrier(CLK_LOCAL_MEM_FENCE); }
   
   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  if (array[giX+32] < tile[tiX]) tiX+32;
      if (array[giX+16] < tile[tiX]) tiX+16;
      if (array[giX+8] < tile[tiX]) tiX+8;
      if (array[giX+4] < tile[tiX]) tiX+4;
      if (array[giX+2] < tile[tiX]) tiX+2;
      if (array[giX+1] < tile[tiX]) tiX+1; }
   
   if (tiX == 0)
   {  result[0] = tile[0]; }
}
*/

void reduction_sum_within_tile(__local  real  *tile);
void reduction_sum_int_within_tile(__local  int  *tile);
void reduction_max_within_tile(__local  real  *tile);
void reduction_min_within_tile(__local  real  *tile);

__kernel void reduce_sum_stage1of2_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const real  *array,     // 1  
        __global       real  *scratch,   // 2 
        __local        real  *tile)      // 3
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);
    
    const unsigned int group_id = get_group_id(0);

    tile[tiX] = ZERO;
    if (giX < isize) tile[tiX] = array[giX];

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_sum_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      scratch[group_id] = tile[0];
    }
}

__kernel void reduce_sum_stage2of2_cl(
        const    int    isize,
        __global real  *scratch,
        __local  real  *tile)
{       
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX] = 0.0;

   if (tiX < isize) tile[tiX] = scratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
     tile[tiX] += scratch[giX];
   }
 
   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_sum_within_tile(tile);

   if (tiX == 0) {
     scratch[0] = tile[0];
   }
}

__kernel void reduce_sum_int_stage1of2_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const int   *array,     // 1  
        __global       int   *scratch,   // 2 
        __local        int   *tile)      // 3
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);
    
    const unsigned int group_id = get_group_id(0);

    tile[tiX] = ZERO;
    if (giX < isize) tile[tiX] = array[giX];

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_sum_int_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      scratch[group_id] = tile[0];
    }
}

__kernel void reduce_sum_int_stage2of2_cl(
        const    int    isize,
        __global int   *scratch,
        __local  int   *tile)
{       
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX] = 0;

   if (tiX < isize) tile[tiX] = scratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
     tile[tiX] += scratch[giX];
   }
 
   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_sum_int_within_tile(tile);

   if (tiX == 0) {
     scratch[0] = tile[0];
   }
}

__kernel void reduce_max_stage1of2_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const real  *array,     // 1  
        __global       real  *scratch,   // 2 
        __local        real  *tile)      // 3
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);
    
    const unsigned int group_id = get_group_id(0);

    tile[tiX] = MINUS_INFINITY;
    if (giX < isize) tile[tiX] = array[giX];

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_max_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      scratch[group_id] = tile[0];
    }
}

__kernel void reduce_max_stage2of2_cl(
        const    int    isize,
        __global real  *scratch,
        __local  real  *tile)
{       
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX] = -1.0e20;

   if (tiX < isize) tile[tiX] = scratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
     if (scratch[giX] > tile[tiX]) tile[tiX] = scratch[giX];
   }
 
   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_max_within_tile(tile);

   if (tiX == 0) {
     scratch[0] = tile[0];
   }
}

__kernel void reduce_min_stage1of2_cl(
                 const int    isize,     // 0  Total number of cells.
        __global const real  *array,     // 1  
        __global       real  *scratch,   // 2 
        __local        real  *tile)      // 3
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);
    
    const unsigned int group_id = get_group_id(0);

    tile[tiX] = PLUS_INFINITY;
    if (giX < isize) tile[tiX] = array[giX];

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_min_within_tile(tile);

    //  Write the local value back to an array size of the number of groups
    if (tiX == 0){
      scratch[group_id] = tile[0];
    }
}

__kernel void reduce_min_stage2of2_cl(
        const    int    isize,
        __global real  *scratch,
        __local  real  *tile)
{       
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   int giX = tiX;

   tile[tiX] = 1.0e20;

   if (tiX < isize) tile[tiX] = scratch[giX];

   for (giX += ntX; giX < isize; giX += ntX) {
     if (scratch[giX] < tile[tiX]) tile[tiX] = scratch[giX];
   }
 
   barrier(CLK_LOCAL_MEM_FENCE);

   reduction_min_within_tile(tile);

   if (tiX == 0) {
     scratch[0] = tile[0];
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

void reduction_max_within_tile(__local  real  *tile) 
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
      if (tile[tiX+16] > tile[tiX]) tile[tiX] = tile[tiX+16];
      if (tile[tiX+8] > tile[tiX]) tile[tiX] = tile[tiX+8];
      if (tile[tiX+4] > tile[tiX]) tile[tiX] = tile[tiX+4];
      if (tile[tiX+2] > tile[tiX]) tile[tiX] = tile[tiX+2];
      if (tile[tiX+1] > tile[tiX]) tile[tiX] = tile[tiX+1];
    }

}

void reduction_min_within_tile(__local  real  *tile) 
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
      if (tile[tiX+16] < tile[tiX]) tile[tiX] = tile[tiX+16];
      if (tile[tiX+8] < tile[tiX]) tile[tiX] = tile[tiX+8];
      if (tile[tiX+4] < tile[tiX]) tile[tiX] = tile[tiX+4];
      if (tile[tiX+2] < tile[tiX]) tile[tiX] = tile[tiX+2];
      if (tile[tiX+1] < tile[tiX]) tile[tiX] = tile[tiX+1];
    }

}
