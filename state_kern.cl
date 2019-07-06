/*
 *  Copyright (c) 2011-2019, Triad National Security, LLC.
 *  All rights Reserved.
 *
 *  CLAMR -- LA-CC-11-094
 *
 *  Copyright 2011-2019. Triad National Security, LLC. This software was produced 
 *  under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Triad National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 *  TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
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
 *     * Neither the name of the Triad National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE TRIAD NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL
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

#if !defined(REG_INTEGER) && !defined(SHORT_INTEGER) && !defined(MIN_INTEGER)
#define REG_INTEGER
#endif

#if defined(MIN_INTEGER)
   // define all to needed ranges and then typedef or define to actual
   typedef unsigned short ushort_t; // 0 to 65,535
   typedef short          short_t;  // -32,768 to 32,767
   typedef unsigned char  uchar_t;  // 0 to 255
   typedef char           char_t;   // -128 to 127 
#ifdef HAVE_OPENCL
   typedef cl_ushort cl_ushort_t;
   typedef cl_short  cl_short_t;
   typedef cl_uchar  cl_uchar_t;
   typedef cl_char   cl_char_t;
#endif

#elif defined(SHORT_INTEGER)
   typedef unsigned short ushort_t;
   typedef short          short_t;
   typedef unsigned short uchar_t;
   typedef short          char_t;
#ifdef HAVE_OPENCL
   typedef cl_ushort cl_ushort_t;
   typedef cl_short  cl_short_t;
   typedef cl_short  cl_uchar_t;
   typedef cl_short  cl_char_t;
#endif

#elif defined(REG_INTEGER)
   typedef unsigned int ushortt_t;
   typedef int          short_t;
   typedef unsigned int uchar_t;
   typedef int          char_t;
#ifdef HAVE_OPENCL
   typedef cl_uint cl_ushort_t;
   typedef cl_int  cl_short_t;
   typedef cl_uint cl_uchar_t;
   typedef cl_int  cl_char_t;
#endif
#endif

#if !defined(FULL_PRECISION) && !defined(MIXED_PRECISION) && !defined(MINIMUM_PRECISION)
#define FULL_PRECISION
#endif
#ifdef NO_CL_DOUBLE
#undef  FULL_PRECISION
#undef  MIXED_PRECISION
#define MINIMUM_PRECISION
#endif

#if defined(MINIMUM_PRECISION)
   typedef float  state_t; // this is for physics state variables ncell in size
   typedef float4 state4_t;
   typedef float  real_t; // this is used for intermediate calculations
   typedef float2 real2_t;
   typedef float4 real4_t;
#define ZERO 0.0f
#define HALF 0.5f
#define QUARTER 0.25f
#define ONE  1.0f
#define GRAVITATIONAL_CONSTANT 9.80f
#define THOUSAND 1000.0f
#define EPSILON 1.0f-30
#define STATE_EPS        15.0f
#define CONSERVATION_EPS    15.0f
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_ONE 1.0f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
   typedef float   state_t;
   typedef float4  state4_t;
   typedef double  real_t;
   typedef double2 real2_t;
   typedef double4 real4_t;
#define ZERO 0.0
#define HALF 0.5
#define QUARTER 0.25
#define ONE  1.0
#define GRAVITATIONAL_CONSTANT 9.80
#define THOUSAND 1000.0
#define EPSILON 1.0e-30
#define STATE_EPS        .02
#define CONSERVATION_EPS    .02
// calc refine is done in single precision
#define REFINE_ONE 1.0f
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(FULL_PRECISION)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
   typedef double  state_t;
   typedef double4 state4_t;
   typedef double  real_t;
   typedef double2 real2_t;
   typedef double4 real4_t;
#define ZERO 0.0
#define HALF 0.5
#define QUARTER 0.25
#define ONE  1.0
#define GRAVITATIONAL_CONSTANT 9.80
#define THOUSAND 1000.0
#define EPSILON 1.0e-30
#define STATE_EPS        .02
#define CONSERVATION_EPS    .02
// calc refine is done in single precision
#define REFINE_ONE 1.0f
#define REFINE_GRADIENT  0.10
#define COARSEN_GRADIENT 0.05
#define REFINE_HALF 0.5
#define REFINE_NEG_THOUSAND -1000.0

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

int SUM_INT(int a, int b)
{
    return a + b;
}

real_t SUM(real_t a, real_t b)
{
    return a + b;
}

real_t MIN(real_t a, real_t b)
{
    return min(a, b);
}

#define REDUCE_IN_TILE(operation, _tile_arr)                                    \
    for (int offset = ntX >> 1; offset > MIN_REDUCE_SYNC_SIZE; offset >>= 1)    \
    {                                                                           \
        if (tiX < offset)                                                       \
        {                                                                       \
            _tile_arr[tiX] = operation(_tile_arr[tiX], _tile_arr[tiX+offset]);  \
        }                                                                       \
        barrier(CLK_LOCAL_MEM_FENCE);                                           \
    }                                                                           \
    if (tiX < MIN_REDUCE_SYNC_SIZE)                                             \
    {                                                                           \
        for (int offset = MIN_REDUCE_SYNC_SIZE; offset > 1; offset >>= 1)       \
        {                                                                       \
            _tile_arr[tiX] = operation(_tile_arr[tiX], _tile_arr[tiX+offset]);  \
            barrier(CLK_LOCAL_MEM_FENCE);                                       \
        }                                                                       \
        _tile_arr[tiX] = operation(_tile_arr[tiX], _tile_arr[tiX+1]);           \
    }

void reduction_min_within_tile1(__local  real_t  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   REDUCE_IN_TILE(MIN, tile);
}

void reduction_sum_within_tile(__local  real_t  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

   REDUCE_IN_TILE(SUM, tile);
}

void reduction_sum_int2_within_tile(__local  int8  *itile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);

    for (int offset = ntX >> 1; offset > MIN_REDUCE_SYNC_SIZE; offset >>= 1)
    {
        if (tiX < offset)
        {
            itile[tiX].s01 += itile[tiX+offset].s01;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (tiX < MIN_REDUCE_SYNC_SIZE)
    {
        for (int offset = MIN_REDUCE_SYNC_SIZE; offset > 1; offset >>= 1)
        {
            itile[tiX].s01 += itile[tiX+offset].s01;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        itile[tiX].s01 += itile[tiX+1].s01;
    }
}

__kernel void set_timestep_cl(
                 const int     ncells,     // 0  Total number of cells.
                 const real_t  sigma,      // 1
        __global const state_t *H_in,      // 2  
        __global const state_t *U_in,      // 3  
        __global const state_t *V_in,      // 4  
        __global const uchar_t *level,     // 5  Array of level information.
        __global const char_t  *celltype,  // 6 
        __global const real_t  *lev_dx,    // 7  
        __global const real_t  *lev_dy,    // 8  
        __global       real_t  *redscratch,// 9 
        __global       real_t  *deltaT,    // 10
        __local        real_t  *tile)      // 11
{
    const unsigned int giX  = get_global_id(0);
    const unsigned int tiX  = get_local_id(0);
    
    tile[tiX] = 1000.0;

    if (giX >= ncells) return;
    
    const unsigned int group_id = get_group_id(0);
    const unsigned int ntX  = get_local_size(0);

    //  Set physical constants.
    const real_t g     = GRAVITATIONAL_CONSTANT;   // gravitational constant
    
    //--MEMORY MANAGEMENT-------------------------------------------------------
    //  Set values for the main cell.
    real_t H       = H_in[giX];
    real_t U       = U_in[giX];
    real_t V       = V_in[giX];
    uchar_t lev    = level[giX];
    char_t type    = celltype[giX];
    
    //--CALCULATIONS------------------------------------------------------------
    if (type == REAL_CELL){
      real_t wavespeed = sqrt(g * H);
      real_t xspeed  = (fabs(U) + wavespeed)/lev_dx[lev];
      real_t yspeed  = (fabs(V) + wavespeed)/lev_dy[lev];
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
                 const int    isize,
        __global       real_t *redscratch,
        __global       real_t *deltaT,
        __local        real_t *tile)
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

//#ifdef __APPLE_CC__
//#define max(a,b) ((a) > (b) ? (a) : (b))
//#define fabs(a) ( (a) < 0 ? -(a) : a)
//#endif
void setup_tile(__local        state4_t *tile,
                __local        int8     *itile,
                         const int      isize,
                __global const state_t  *H,
                __global const state_t  *U,
                __global const state_t  *V,
                __global const int      *nlft,
                __global const int      *nrht,
                __global const int      *ntop,
                __global const int      *nbot,
                __global const uchar_t  *level
                );

void setup_refine_tile(
                __local        state_t *tile,
                __local        int8    *itile,
                         const int     isize,
                __global const state_t *H,
                __global const int     *nlft,
                __global const int     *nrht,
                __global const int     *ntop,
                __global const int     *nbot,
                __global const uchar_t *level
                );

void setup_xface(
                __local           int8        *xface,
                __global    const int         *map_xface2cell_lower,   
                __global    const int         *map_xface2cell_upper
                );

void setup_yface(
                __local           int8        *yface,
                __global    const int         *map_yface2cell_lower,       
                __global    const int         *map_yface2cell_upper
                );

__kernel void copy_state_data_cl(
                          const int    isize,         // 0 
                 __global      state_t *H,            // 1 
                 __global      state_t *U,            // 2 
                 __global      state_t *V,            // 3 
                 __global      state_t *H_new,        // 4 
                 __global      state_t *U_new,        // 5 
                 __global      state_t *V_new)        // 6 
{
                
   const uint giX  = get_global_id(0);

   if (giX >= isize) return;

   H_new[giX]    = H[giX];
   U_new[giX]    = U[giX];
   V_new[giX]    = V[giX];
}

__kernel void copy_state_ghost_data_cl(
                          const int   ncells,        // 0 
                          const int   nghost,        // 1 
                 __global      state_t *H,            // 2 
                 __global      state_t *H_add,        // 3 
                 __global      state_t *U,            // 4 
                 __global      state_t *U_add,        // 5 
                 __global      state_t *V,            // 6 
                 __global      state_t *V_add)        // 7 
{
   const uint giX  = get_global_id(0);

   if (giX >= nghost) return;

   H[ncells+giX] = H_add[giX];
   U[ncells+giX] = U_add[giX];
   V[ncells+giX] = V_add[giX];
}


#ifndef SET_TILE_VARIABLES
#define SET_TILE_VARIABLES
//  Define macros for local tile access.

#define Hval(i)     ( tile[i].s0 )
#define Uval(i)     ( tile[i].s1 )
#define Vval(i)     ( tile[i].s2 )

#define Hrefval(i)     ( tile[i] )

#define nlftval(i)  ( itile[i].s0 )
#define nrhtval(i)  ( itile[i].s1 )
#define ntopval(i)  ( itile[i].s2 )
#define nbotval(i)  ( itile[i].s3 )
#define levelval(i) ( itile[i].s4 )
#define mpotval(i)  ( itile[i].s5 )

#endif

#ifndef SET_FACE_VARIABLES
#define SET_FACE_VARIABLES
// Define macros for local face access

//#define xfacelower(i)   ( xface[i].s0 )
//#define xfaceupper(i)   ( xface[i].s1 )
//#define yfacelower(i)   ( yface[i].s0 )
//#define yfaceupper(i)   ( yface[i].s1 )

#endif

#define SQ(x)      ( (x)*(x) )
#define MIN3(a,b,c) ( min(min((a),(b)),(c)) )

#define HXFLUX(ic)  ( U[ic] )
#define UXFLUX(ic)  ( SQ(U[ic])/H[ic] + ghalf*SQ(H[ic]) )
#define UVFLUX(ic)  ( U[ic]*V[ic]/H[ic] )

#define HXFLUXIC ( Uic )
#define HXFLUXNL ( Ul )
#define HXFLUXNR ( Ur )
#define HXFLUXNB ( Ub )
#define HXFLUXNT ( Ut )

#define UXFLUXIC ( SQ(Uic)/Hic + ghalf*SQ(Hic) )
#define UXFLUXNL ( SQ(Ul)/Hl + ghalf*SQ(Hl) )
#define UXFLUXNR ( SQ(Ur)/Hr + ghalf*SQ(Hr) )
#define UXFLUXNB ( SQ(Ub)/Hb + ghalf*SQ(Hb) )
#define UXFLUXNT ( SQ(Ut)/Ht + ghalf*SQ(Ht) )

#define UVFLUXIC ( Uic*Vic/Hic )
#define UVFLUXNL ( Ul*Vl/Hl )
#define UVFLUXNR ( Ur*Vr/Hr )
#define UVFLUXNB ( Ub*Vb/Hb )
#define UVFLUXNT ( Ut*Vt/Ht )

#define HYFLUX(ic)  ( V[ic] )
#define VUFLUX(ic)  ( V[ic]*U[ic]/H[ic] )
#define VYFLUX(ic)  ( SQ(V[ic])/H[ic] + ghalf*SQ(H[ic]) )

#define HYFLUXIC ( Vic )
#define HYFLUXNL ( Vl )
#define HYFLUXNR ( Vr )
#define HYFLUXNB ( Vb )
#define HYFLUXNT ( Vt )

#define VUFLUXIC  ( Vic*Uic/Hic )
#define VUFLUXNL  ( Vl*Ul/Hl )
#define VUFLUXNR  ( Vr*Ur/Hr )
#define VUFLUXNB  ( Vb*Ub/Hb )
#define VUFLUXNT  ( Vt*Ut/Ht )

#define VYFLUXIC  ( SQ(Vic)/Hic + ghalf*SQ(Hic) )
#define VYFLUXNL  ( SQ(Vl)/Hl + ghalf*SQ(Hl) )
#define VYFLUXNR  ( SQ(Vr)/Hr + ghalf*SQ(Hr) )
#define VYFLUXNB  ( SQ(Vb)/Hb + ghalf*SQ(Hb) )
#define VYFLUXNT  ( SQ(Vt)/Ht + ghalf*SQ(Ht) )

#define HNEWXFLUXMINUS  ( Uxminus )
#define HNEWXFLUXPLUS   ( Uxplus )
#define UNEWXFLUXMINUS  ( SQ(Uxminus)/Hxminus + ghalf*SQ(Hxminus) )
#define UNEWXFLUXPLUS   ( SQ(Uxplus) /Hxplus +  ghalf*SQ(Hxplus)  )
#define UVNEWFLUXMINUS  ( Uxminus*Vxminus/Hxminus )
#define UVNEWFLUXPLUS   ( Uxplus *Vxplus /Hxplus  )

#define HNEWYFLUXMINUS  ( Vyminus )
#define HNEWYFLUXPLUS   ( Vyplus  )
#define VNEWYFLUXMINUS  ( SQ(Vyminus)/Hyminus + ghalf*SQ(Hyminus) )
#define VNEWYFLUXPLUS   ( SQ(Vyplus) /Hyplus  + ghalf*SQ(Hyplus)  )
#define VUNEWFLUXMINUS  ( Vyminus*Uyminus/Hyminus )
#define VUNEWFLUXPLUS   ( Vyplus *Uyplus /Hyplus )

#define HXFLUXFACE (Ux)
#define UXFLUXFACE (SQ(Ux)/Hx + ghalf*SQ(Hx))
#define VXFLUXFACE (Ux*Vx/Hx)

#define HYFLUXFACE (Vy)
#define UYFLUXFACE (Vy*Uy/Hy)
#define VYFLUXFACE (SQ(Vy)/Hy + ghalf*SQ(Hy))

// XXX Added XXX
#define HXFLUXNLT ( Ult )
#define HXFLUXNRT ( Urt )
#define UXFLUXNLT ( SQ(Ult)/Hlt + ghalf*SQ(Hlt) )
#define UXFLUXNRT ( SQ(Urt)/Hrt + ghalf*SQ(Hrt) )
#define UVFLUXNLT ( Ult*Vlt/Hlt )
#define UVFLUXNRT ( Urt*Vrt/Hrt )
#define HYFLUXNBR ( Vbr )
#define HYFLUXNTR ( Vtr )
#define VUFLUXNBR  ( Vbr*Ubr/Hbr )
#define VUFLUXNTR  ( Vtr*Utr/Htr )
#define VYFLUXNBR  ( SQ(Vbr)/Hbr + ghalf*SQ(Hbr) )
#define VYFLUXNTR  ( SQ(Vtr)/Htr + ghalf*SQ(Htr) )
#define HNEWXFLUXMINUS2  ( Uxminus2 )
#define HNEWXFLUXPLUS2   ( Uxplus2 )
#define UNEWXFLUXMINUS2  ( SQ(Uxminus2)/Hxminus2 + ghalf*SQ(Hxminus2) )
#define UNEWXFLUXPLUS2   ( SQ(Uxplus2) /Hxplus2 +  ghalf*SQ(Hxplus2)  )
#define UVNEWFLUXMINUS2  ( Uxminus2*Vxminus2/Hxminus2 )
#define UVNEWFLUXPLUS2   ( Uxplus2 *Vxplus2 /Hxplus2  )
#define HNEWYFLUXMINUS2  ( Vyminus2 )
#define HNEWYFLUXPLUS2   ( Vyplus2  )
#define VNEWYFLUXMINUS2  ( SQ(Vyminus2)/Hyminus2 + ghalf*SQ(Hyminus2) )
#define VNEWYFLUXPLUS2   ( SQ(Vyplus2) /Hyplus2  + ghalf*SQ(Hyplus2)  )
#define VUNEWFLUXMINUS2  ( Vyminus2*Uyminus2/Hyminus2 )
#define VUNEWFLUXPLUS2   ( Vyplus2 *Uyplus2 /Hyplus2 )

#define U_halfstep(deltaT, U_i, U_n, F_i, F_n, r_i, r_n, A_i, A_n, V_i, V_n) (( (( r_i*U_n + r_n*U_i ) / ( r_i + r_n )) - HALF*deltaT*(( F_n*A_n*min(ONE, A_i/A_n) - F_i*A_i*min(ONE, A_n/A_i) ) / ( V_n*min(HALF, V_i/V_n) + V_i*min(HALF, V_n/V_i) )) ))

real_t U_halfstep_ORIG(// XXX Fix the subindices to be more intuitive XXX
        real_t  deltaT,     // Timestep
        real_t  U_i,        // Initial cell's (downwind's) state variable
        real_t  U_n,        // Next cell's    (upwind's)   state variable
        real_t  F_i,        // Initial cell's (downwind's) state variable flux
        real_t  F_n,        // Next cell's    (upwind's)   state variable flux
        real_t  r_i,        // Initial cell's (downwind's) center to face dist
        real_t  r_n,        // Next cell's    (upwind's)   center to face dist
        real_t  A_i,        // Cell's            face surface area
        real_t  A_n,        // Cell's neighbor's face surface area
        real_t  V_i,        // Cell's            volume
        real_t  V_n) {      // Cell's neighbor's volume

   return (( r_i*U_n + r_n*U_i ) / ( r_i + r_n ))
          - HALF*deltaT*(( F_n*A_n*min(ONE, A_i/A_n) - F_i*A_i*min(ONE, A_n/A_i) )
                    / ( V_n*min(HALF, V_i/V_n) + V_i*min(HALF, V_n/V_i) ));

}

real_t U_halfstep_BD(// XXX Fix the subindices to be more intuitive XXX
        real_t  deltaT,     // Timestep
        real_t  U_i,        // Initial cell's (downwind's) state variable
        real_t  U_n,        // Next cell's    (upwind's)   state variable
        real_t  F_i,        // Initial cell's (downwind's) state variable flux
        real_t  F_n,        // Next cell's    (upwind's)   state variable flux
        real_t  r_i,        // Initial cell's (downwind's) center to face distance
        real_t  r_n,        // Next cell's    (upwind's)   center to face distance
        real_t  A_i,        // Cell's            face surface area
        real_t  A_n,        // Cell's neighbor's face surface area
        real_t  V_i,        // Cell's            volume
        real_t  V_n) {      // Cell's neighbor's volume

   return ((U_i*r_n+U_n*r_i)/(r_n+r_i) + (deltaT/(r_n+r_i))*(F_n-F_i));

}

#define U_fullstep(deltaT, dr, U, F_plus, F_minus, G_plus, G_minus) (( (U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus)) ))

real_t U_fullstep_ORIG(
        real_t  deltaT,
        real_t  dr,
        real_t  U,
        real_t  F_plus,
        real_t  F_minus,
        real_t  G_plus,
        real_t  G_minus) {

   return (U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus));

}

//#define w_corrector(deltaT, dr, U_eigen, grad_half, grad_minus, grad_plus) (( HALF*(HALF*U_eigen*deltaT/dr)*(ONE-(HALF*U_eigen*deltaT/dr))*(ONE- max(MIN3(ONE, (grad_plus*grad_half/max(SQ(grad_half),EPSILON)), (grad_minus*grad_half/max(SQ(grad_half),EPSILON))), ZERO)) ))

real_t w_corrector(//_ORIG(
        real_t  deltaT,       // Timestep
        real_t  dr,           // Cell's center to face distance
        real_t  U_eigen,      // State variable's eigenvalue (speed)
        real_t  grad_half,    // Centered gradient
        real_t  grad_minus,   // Downwind gradient
        real_t  grad_plus) {  // Upwind gradient

   real_t nu     = HALF * U_eigen * deltaT / dr;
   nu          = nu * (ONE - nu);

   real_t rdenom = ONE / max(SQ(grad_half), EPSILON);
   real_t rplus  = (grad_plus  * grad_half) * rdenom;
   real_t rminus = (grad_minus * grad_half) * rdenom;

   return HALF*nu*(ONE- max(MIN3(ONE, rplus, rminus), ZERO));

}

__kernel void apply_boundary_conditions_local_cl(
                 const int      ncells,   // 0  Total number of cells
        __global const char_t   *celltype, // 1  Array of left neighbors
        __global const int     *nlft,     // 2  Array of left neighbors
        __global const int     *nrht,     // 3  Array of right neighbors
        __global const int     *ntop,     // 4  Array of top neighbors
        __global const int     *nbot,     // 5  Array of bottom neighbors
        __global       state_t *H,        // 6  H array
        __global       state_t *U,        // 7  U array
        __global       state_t *V)        // 8  V array
{
   const uint giX  = get_global_id(0);
   const uint tiX  = get_local_id(0);

   // Ensure the executing thread is not extraneous
   if(giX >= ncells)
      return;

   char_t ctype = celltype[giX];

   if (ctype == LEFT_BOUNDARY){
      int nr = nrht[giX];
      if (nr < (int)ncells) {
         H[giX] =  H[nr];
         U[giX] = -U[nr];
         V[giX] =  V[nr];
      }
   }
   if (ctype == RIGHT_BOUNDARY){
      int nl = nlft[giX];
      if (nl < (int)ncells) {
         H[giX] =  H[nl];
         U[giX] = -U[nl];
         V[giX] =  V[nl];
      }
   }
   if (ctype == BOTTOM_BOUNDARY){
      int nt = ntop[giX];
      if (nt < (int)ncells) {
         H[giX] =  H[nt];
         U[giX] =  U[nt];
         V[giX] = -V[nt];
      }
   }
   if (ctype == TOP_BOUNDARY){
      int nb = nbot[giX];
      if (nb < (int)ncells) {
         H[giX] =  H[nb];
         U[giX] =  U[nb];
         V[giX] = -V[nb];
      }
   }
}

__kernel void apply_boundary_conditions_ghost_cl(
                 const int     ncells,   // 0  Total number of cells
        __global const char_t  *celltype, // 1  Array celltypes
        __global const int    *nlft,     // 2  Array of left neighbors
        __global const int    *nrht,     // 3  Array of right neighbors
        __global const int    *ntop,     // 4  Array of top neighbors
        __global const int    *nbot,     // 5  Array of bottom neighbors
        __global       real_t *H,        // 6  H array
        __global       real_t *U,        // 7  U array
        __global       real_t *V)        // 8  V array
{
   const unsigned int giX  = get_global_id(0);
   const unsigned int tiX  = get_local_id(0);

   // Ensure the executing thread is not extraneous
   if(giX >= ncells)
      return;

   char_t ctype = celltype[giX];

   if (ctype == LEFT_BOUNDARY){
      int nr = nrht[giX];
      if (nr >= (int)ncells) {
         H[giX] =  H[nr];
         U[giX] = -U[nr];
         V[giX] =  V[nr];
      }
   }
   if (ctype == RIGHT_BOUNDARY){
      int nl = nlft[giX];
      if (nl >= (int)ncells) {
         H[giX] =  H[nl];
         U[giX] = -U[nl];
         V[giX] =  V[nl];
      }
   }
   if (ctype == BOTTOM_BOUNDARY){
      int nt = ntop[giX];
      if (nt >= (int)ncells) {
         H[giX] =  H[nt];
         U[giX] =  U[nt];
         V[giX] = -V[nt];
      }
   }
   if (ctype == TOP_BOUNDARY){
      int nb = nbot[giX];
      if (nb >= (int)ncells) {
         H[giX] =  H[nb];
         U[giX] =  U[nb];
         V[giX] = -V[nb];
      }
   }
}

__kernel void apply_boundary_conditions_cl(
                 const int      ncells,   // 0  Total number of cells
        __global const char_t   *celltype, // 1  Array of left neighbors
        __global const int     *nlft,     // 2  Array of left neighbors
        __global const int     *nrht,     // 3  Array of right neighbors
        __global const int     *ntop,     // 4  Array of top neighbors
        __global const int     *nbot,     // 5  Array of bottom neighbors
        __global       state_t *H,        // 6  H array
        __global       state_t *U,        // 7  U array
        __global       state_t *V)        // 8  V array
{
   const uint giX  = get_global_id(0);
   const uint tiX  = get_local_id(0);

   // Ensure the executing thread is not extraneous
   if(giX >= ncells)
      return;

   char_t ctype = celltype[giX];

   if (ctype == LEFT_BOUNDARY){
      int nr = nrht[giX];
      H[giX] =  H[nr];
      U[giX] = -U[nr];
      V[giX] =  V[nr];
   }
   if (ctype == RIGHT_BOUNDARY){
      int nl = nlft[giX];
      H[giX] =  H[nl];
      U[giX] = -U[nl];
      V[giX] =  V[nl];
   }
   if (ctype == BOTTOM_BOUNDARY){
      int nt = ntop[giX];
      H[giX] =  H[nt];
      U[giX] =  U[nt];
      V[giX] = -V[nt];
   }
   if (ctype == TOP_BOUNDARY){
      int nb = nbot[giX];
      H[giX] =  H[nb];
      U[giX] =  U[nb];
      V[giX] = -V[nb];
   }
}



__kernel void calc_finite_difference_cl(
                 const int       ncells,   // 0  Total number of cells
                 const int       levmx,    // 1  Maximum level
        __global const state_t  *H,        // 2  
        __global const state_t  *U,        // 3  
        __global const state_t  *V,        // 4  
        __global       state_t  *H_new,    // 5  
        __global       state_t  *U_new,    // 6  
        __global       state_t  *V_new,    // 7  
        __global const int      *nlft,     // 8   Array of left neighbors
        __global const int      *nrht,     // 9   Array of right neighbors
        __global const int      *ntop,     // 10  Array of top neighbors
        __global const int      *nbot,     // 11  Array of bottom neighbors
        __global const uchar_t  *level,    // 12  Array of level information
                 const real_t    deltaT,   // 13  Size of time step.
        __global const real_t   *lev_dx,   // 14
        __global const real_t   *lev_dy,   // 15
        __local        state4_t *tile,     // 16  Tile size in state4_t
        __local        int8     *itile){   // 17  Tile size in int8

   /////////////////////////////////////////////
   /// Get thread identification information ///
   /////////////////////////////////////////////

   const uint giX  = get_global_id(0);
   const uint tiX  = get_local_id(0);
   
   const uint ngX  = get_global_size(0);
   const uint ntX  = get_local_size(0);
   
   const uint group_id = get_group_id(0);
    
   // Ensure the executing thread is not extraneous
   if(giX >= ncells)
      return;

   /////////////////////////////////////////////
   /// Set local tile & apply boundary conds ///
   /////////////////////////////////////////////

   setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);

   barrier(CLK_LOCAL_MEM_FENCE);

   /////////////////////////////////////////////////
   /// Declare all constants and local variables ///
   /////////////////////////////////////////////////

   const real_t g     = GRAVITATIONAL_CONSTANT;   // gravitational constant
   const real_t ghalf = HALF*g;

   // Left, right, ... left-left, right-right, ... left-top, right-top neighbor
   int nl, nr, nt, nb;
   int nll, nrr, ntt, nbb;

   // Level
   int lvl, lvl_nl, lvl_nr, lvl_nt, lvl_nb;
   int lvl_nll, lvl_nrr, lvl_ntt, lvl_nbb;

   // Left-top, right-top, top-right, bottom-right neighbor 
   int nlt, nrt, ntr, nbr;

   // State variables at x-axis control volume face
   real_t Hxminus, Hxplus;
   real_t Uxminus, Uxplus;
   real_t Vxminus, Vxplus;

   // State variables at y-axis control volume face
   real_t Hyminus, Hyplus;
   real_t Uyminus, Uyplus;
   real_t Vyminus, Vyplus;

   // Variables for artificial viscosity/flux limiting
   real_t wminusx_H, wminusx_U;
   real_t wplusx_H, wplusx_U;
   real_t wminusy_H, wminusy_V;
   real_t wplusy_H, wplusy_V;

   int nltl;
   real_t Hll2;

   int nrtr;
   real_t Hrr2;

   real_t Ull2;
   real_t Urr2;

   int ntrt;
   real_t Htt2;

   int nbrb;
   real_t Hbb2;

   real_t Vtt2;
   real_t Vbb2;

   real_t Hxminus2, Hxplus2;
   real_t Uxminus2, Uxplus2;
   real_t Vxminus2, Vxplus2;

   real_t Hyminus2, Hyplus2;
   real_t Uyminus2, Uyplus2;
   real_t Vyminus2, Vyplus2;

   real_t Hxfluxminus;
   real_t Uxfluxminus;
   real_t Vxfluxminus;

   real_t Hxfluxplus;
   real_t Uxfluxplus;
   real_t Vxfluxplus;

   real_t Hyfluxminus;
   real_t Uyfluxminus;
   real_t Vyfluxminus;

   real_t Hyfluxplus;
   real_t Uyfluxplus;
   real_t Vyfluxplus;


   // XXX Assuming square cells! XXX
   // State variables and cell widths and lengths
   real_t dric, drl, drr, drt, drb;
//   real_t drlt, drrt, drtr, drbr;

   real_t Hic, Hl, Hr, Ht, Hb;
   real_t Hll, Hrr, Htt, Hbb;

   real_t Uic, Ul, Ur, Ut, Ub;
   real_t Ull, Urr;

   real_t Vic, Vl, Vr, Vt, Vb;
   real_t Vtt, Vbb;

   real_t Hlt, Hrt, Htr, Hbr;
   real_t Ult, Urt, Utr, Ubr;
   real_t Vlt, Vrt, Vtr, Vbr;


   // Local values for the state variables and cell widths and heights for the local cell as well
   // as its neighboring cells
   real_t dxic, dxl, dxr, dyic, dyt, dyb;

   //////////////////////////
   /// Set the local tile ///
   //////////////////////////

   int start_idx = group_id * ntX;
//   int end_idx = (group_id + 1) * ntX;

   lvl  = levelval(tiX);

   dxic = lev_dx[lvl];
   dyic = lev_dy[lvl];

   nl = nlftval(tiX);
   nr = nrhtval(tiX);
   nt = ntopval(tiX);
   nb = nbotval(tiX);

//   nl = nlft[ic];
//   nr = nrht[ic];
//   nt = ntop[ic];
//   nb = nbot[ic];


   dric = dxic;

   Hic  = Hval(tiX);
   Uic  = Uval(tiX);
   Vic  = Vval(tiX);


   int glob_flag = 0;
   // Storing all values associated with the left neighbors in local variables
   if(nl < 0) {
      nl      = abs(nl+1);
      lvl_nl  = level[nl];
      nll     = nlft[nl];
      Hl      = H[nl];
      Ul      = U[nl];
      Vl      = V[nl];
      dxl     = lev_dx[level[nl]]; 
      nlt     = ntop[nl];
      glob_flag = 1;
   }
   else {
      lvl_nl  = levelval(nl);
      nll     = nlftval(nl);
      Hl      = Hval(nl);
      Ul      = Uval(nl);
      Vl      = Vval(nl);
      dxl     = lev_dx[levelval(nl)];
      nlt     = ntopval(nl);
   }
   drl = dxl; // lev_dx[level[nl]]; 
   if(nll < 0 || glob_flag == 1) {
      if (nll < 0) nll     = abs(nll+1);
      lvl_nll = level[nll];
      Hll     = H[nll];
      Ull     = U[nll];
   }
   else {
      lvl_nll = levelval(nll);
      Hll     = Hval(nll);
      Ull     = Uval(nll);
      nll    += start_idx;
   }


   if(nlt < 0 || glob_flag == 1) {
      if (nlt < 0) nlt     = abs(nlt+1);
      glob_flag = 1;
   }

   if(lvl < lvl_nl) {
      if(glob_flag == 1) {
         Hlt     = H[nlt];
         Ult     = U[nlt];
         Vlt     = V[nlt];
         nltl    = nlft[nlt];
      }
      else {
         Hlt     = Hval(nlt);
         Ult     = Uval(nlt);
         Vlt     = Vval(nlt);
         nltl    = nlftval(nlt);
      }

      if(nltl < 0 || glob_flag == 1) {
         if (nltl < 0) nltl    = abs(nltl+1);
         Hll2    = H[nltl];
         Ull2    = U[nltl];
      }
      else {
         Hll2    = Hval(nltl);
         Ull2    = Uval(nltl);
         nltl   += start_idx;
      }
   }


   glob_flag = 0;
   // Storing all values associated with the right neighbors in local variables
   if(nr < 0) {
      nr      = abs(nr+1);
      lvl_nr  = level[nr];
      nrr     = nrht[nr];
      Hr      = H[nr];
      Ur      = U[nr];
      Vr      = V[nr];
      dxr     = lev_dx[level[nr]]; 
      nrt     = ntop[nr];
      glob_flag = 1;
   }
   else {
      lvl_nr  = levelval(nr);
      nrr     = nrhtval(nr);
      Hr      = Hval(nr);
      Ur      = Uval(nr);
      Vr      = Vval(nr);
      dxr     = lev_dx[levelval(nr)];
      nrt     = ntopval(nr);
   }
   drr = dxr; // lev_dx[level[nr]]; 
   if(nrr < 0 || glob_flag == 1) {
      if (nrr < 0) nrr     = abs(nrr+1);
      lvl_nrr = level[nrr];
      Hrr     = H[nrr];
      Urr     = U[nrr];
   }
   else {
      lvl_nrr = levelval(nrr);
      Hrr     = Hval(nrr);
      Urr     = Uval(nrr);
      nrr    += start_idx;
   }

   if(nrt < 0 || glob_flag == 1) {
      if (nrt < 0) nrt     = abs(nrt+1);
      glob_flag = 1;
   }


   if(lvl < lvl_nr) {
      if(glob_flag == 1) {
         Hrt     = H[nrt];
         Urt     = U[nrt];
         Vrt     = V[nrt];
         nrtr    = nrht[nrt];
      }
      else {
         Hrt     = Hval(nrt);
         Urt     = Uval(nrt);
         Vrt     = Vval(nrt);
         nrtr    = nrhtval(nrt);
      }

      if(nrtr < 0 || glob_flag == 1) {
         if (nrtr < 0) nrtr    = abs(nrtr+1);
         Hrr2    = H[nrtr];
         Urr2    = U[nrtr];
      }
      else {
         Hrr2    = Hval(nrtr);
         Urr2    = Uval(nrtr);
         nrtr   += start_idx;
      }
   }


   glob_flag = 0;
   // Storing all values associated with the top neighbors in local variables
   if(nt < 0) {
      nt      = abs(nt+1);
      lvl_nt  = level[nt];
      ntt     = ntop[nt];
      Ht      = H[nt];
      Ut      = U[nt];
      Vt      = V[nt];
      dyt     = lev_dy[level[nt]]; 
      ntr     = nrht[nt];
      glob_flag = 1;
   }
   else {
      lvl_nt  = levelval(nt);
      ntt     = ntopval(nt);
      Ht      = Hval(nt);
      Ut      = Uval(nt);
      Vt      = Vval(nt);
      dyt     = lev_dy[levelval(nt)];
      ntr     = nrhtval(nt);
   }
   drt = dyt; // lev_dy[level[nt]]; 
   if(ntt < 0 || glob_flag == 1) {
      if (ntt < 0) ntt     = abs(ntt+1);
      lvl_ntt = level[ntt];
      Htt     = H[ntt];
      Vtt     = V[ntt];
   }
   else {
      lvl_ntt = levelval(ntt);
      Htt     = Hval(ntt);
      Vtt     = Vval(ntt);
      ntt    += start_idx;
   }

   if(ntr < 0 || glob_flag == 1) {
      if (ntr < 0) ntr     = abs(ntr+1);
      glob_flag = 1;
   }

   if(lvl < lvl_nt) {
      if(glob_flag == 1) {
         Htr     = H[ntr];
         Utr     = U[ntr];
         Vtr     = V[ntr];
         ntrt    = ntop[ntr];
      }
      else {
         Htr     = Hval(ntr);
         Utr     = Uval(ntr);
         Vtr     = Vval(ntr);
         ntrt    = ntopval(ntr);
      }

      if(ntrt < 0 || glob_flag == 1) {
         if (ntrt < 0) ntrt    = abs(ntrt+1);
         Htt2    = H[ntrt];
         Vtt2    = V[ntrt];
      }
      else {
         Htt2    = Hval(ntrt);
         Vtt2    = Vval(ntrt);
         ntrt   += start_idx;
      }
   }


   glob_flag = 0;
   // Storing all values associated with the bottom neighbors in local variables
   if(nb < 0) {
      nb      = abs(nb+1);
      lvl_nb  = level[nb];
      nbb     = nbot[nb];
      Hb      = H[nb];
      Ub      = U[nb];
      Vb      = V[nb];
      dyb     = lev_dy[level[nb]]; 
      nbr     = nrht[nb];
      glob_flag = 1;
   }
   else {
      lvl_nb  = levelval(nb);
      nbb     = nbotval(nb);
      Hb      = Hval(nb);
      Ub      = Uval(nb);
      Vb      = Vval(nb);
      dyb     = lev_dy[levelval(nb)];
      nbr     = nrhtval(nb);
   }
   drb = dyb; // lev_dy[level[nb]]; 
   if(nbb < 0 || glob_flag == 1) {
      if (nbb < 0) nbb     = abs(nbb+1);
      lvl_nbb = level[nbb];
      Hbb     = H[nbb];
      Vbb     = V[nbb];
   }
   else {
      lvl_nbb = levelval(nbb);
      Hbb     = Hval(nbb);
      Vbb     = Vval(nbb);
      nbb    += start_idx;
   }

   if(nbr < 0 || glob_flag == 1) {
      if (nbr < 0) nbr     = abs(nbr+1);
      glob_flag = 1;
   }

   if(lvl < lvl_nb) {
      if(glob_flag == 1) {
         Hbr     = H[nbr];
         Ubr     = U[nbr];
         Vbr     = V[nbr];
         nbrb    = nbot[nbr];
      }
      else {
         Hbr     = Hval(nbr);
         Ubr     = Uval(nbr);
         Vbr     = Vval(nbr);
         nbrb    = nbotval(nbr);
      }

      if(nbrb < 0 || glob_flag == 1) {
         if (nbrb < 0) nbrb    = abs(nbrb+1);
         Hbb2    = H[nbrb];
         Vbb2    = V[nbrb];
      }
      else {
         Hbb2    = Hval(nbrb);
         Vbb2    = Vval(nbrb);
         nbrb   += start_idx;
      }
   }


   /////////////////////////////////////////////
   ///    Half time-step for Lax-Wendroff    ///
   /////////////////////////////////////////////

//#define V_delta(V_i, V_n) (( HALF*deltaT / ((V_n)*min(HALF, (V_i)/(V_n))+(V_i)*min(HALF, (V_n)/(V_i))) ))

//   real_t V_stag = V_delta( SQ(dxl), SQ(dxic) );

   Hxminus = U_halfstep(deltaT, Hl, Hic, HXFLUXNL, HXFLUXIC, 
                        dxl, dxic, dxl, dxic, SQ(dxl), SQ(dxic));
   Uxminus = U_halfstep(deltaT, Ul, Uic, UXFLUXNL, UXFLUXIC,
                        dxl, dxic, dxl, dxic, SQ(dxl), SQ(dxic));
   Vxminus = U_halfstep(deltaT, Vl, Vic, UVFLUXNL, UVFLUXIC,
                        dxl, dxic, dxl, dxic, SQ(dxl), SQ(dxic));

//   V_stag = V_delta( SQ(dxic), SQ(dxr) );

   Hxplus  = U_halfstep(deltaT, Hic, Hr, HXFLUXIC, HXFLUXNR,
                        dxic, dxr, dxic, dxr, SQ(dxic), SQ(dxr));
   Uxplus  = U_halfstep(deltaT, Uic, Ur, UXFLUXIC, UXFLUXNR,
                        dxic, dxr, dxic, dxr, SQ(dxic), SQ(dxr));
   Vxplus  = U_halfstep(deltaT, Vic, Vr, UVFLUXIC, UVFLUXNR,
                        dxic, dxr, dxic, dxr, SQ(dxic), SQ(dxr));

//   V_stag = V_delta( SQ(dyb), SQ(dyic) );

   Hyminus = U_halfstep(deltaT, Hb, Hic, HYFLUXNB, HYFLUXIC,
                        dyb, dyic, dyb, dyic, SQ(dyb), SQ(dyic));
   Uyminus = U_halfstep(deltaT, Ub, Uic, VUFLUXNB, VUFLUXIC,
                        dyb, dyic, dyb, dyic, SQ(dyb), SQ(dyic));
   Vyminus = U_halfstep(deltaT, Vb, Vic, VYFLUXNB, VYFLUXIC,
                        dyb, dyic, dyb, dyic, SQ(dyb), SQ(dyic));

//   V_stag = V_delta( SQ(dyic), SQ(dyt) );

   Hyplus  = U_halfstep(deltaT, Hic, Ht, HYFLUXIC, HYFLUXNT,
                        dyic, dyt, dyic, dyt, SQ(dyic), SQ(dyt));
   Uyplus  = U_halfstep(deltaT, Uic, Ut, VUFLUXIC, VUFLUXNT,
                        dyic, dyt, dyic, dyt, SQ(dyic), SQ(dyt));
   Vyplus  = U_halfstep(deltaT, Vic, Vt, VYFLUXIC, VYFLUXNT,
                        dyic, dyt, dyic, dyt, SQ(dyic), SQ(dyt));


   Hxfluxminus = HNEWXFLUXMINUS;
   Uxfluxminus = UNEWXFLUXMINUS;
   Vxfluxminus = UVNEWFLUXMINUS;

   Hxfluxplus  = HNEWXFLUXPLUS;
   Uxfluxplus  = UNEWXFLUXPLUS;
   Vxfluxplus  = UVNEWFLUXPLUS;

   Hyfluxminus = HNEWYFLUXMINUS;
   Uyfluxminus = VUNEWFLUXMINUS;
   Vyfluxminus = VNEWYFLUXMINUS;

   Hyfluxplus  = HNEWYFLUXPLUS;
   Uyfluxplus  = VUNEWFLUXPLUS;
   Vyfluxplus  = VNEWYFLUXPLUS;


   if(lvl < lvl_nl) {

//      V_stag = V_delta( SQ(drl), SQ(dric) );

      Hxminus2 = U_halfstep(deltaT, Hlt, Hic, HXFLUXNLT, HXFLUXIC,
                            drl, dric, drl, dric, SQ(drl), SQ(dric));
      Uxminus2 = U_halfstep(deltaT, Ult, Uic, UXFLUXNLT, UXFLUXIC,
                            drl, dric, drl, dric, SQ(drl), SQ(dric));
      Vxminus2 = U_halfstep(deltaT, Vlt, Vic, UVFLUXNLT, UVFLUXIC,
                            drl, dric, drl, dric, SQ(drl), SQ(dric));

      Hxfluxminus = (Hxfluxminus + HNEWXFLUXMINUS2) * HALF;
      Uxfluxminus = (Uxfluxminus + UNEWXFLUXMINUS2) * HALF;
      Vxfluxminus = (Vxfluxminus + UVNEWFLUXMINUS2) * HALF;

   }

   if(lvl < lvl_nr) {

//      V_stag = V_delta( SQ(dric), SQ(drr) );

      Hxplus2  = U_halfstep(deltaT, Hic, Hrt, HXFLUXIC, HXFLUXNRT,
                            dric, drr, dric, drr, SQ(dric), SQ(drr));
      Uxplus2  = U_halfstep(deltaT, Uic, Urt, UXFLUXIC, UXFLUXNRT,
                            dric, drr, dric, drr, SQ(dric), SQ(drr));
      Vxplus2  = U_halfstep(deltaT, Vic, Vrt, UVFLUXIC, UVFLUXNRT,
                            dric, drr, dric, drr, SQ(dric), SQ(drr));

      Hxfluxplus  = (Hxfluxplus + HNEWXFLUXPLUS2) * HALF;
      Uxfluxplus  = (Uxfluxplus + UNEWXFLUXPLUS2) * HALF;
      Vxfluxplus  = (Vxfluxplus + UVNEWFLUXPLUS2) * HALF;

   }

   if(lvl < lvl_nb) {

//      V_stag = V_delta( SQ(drb), SQ(dric) );

      Hyminus2 = U_halfstep(deltaT, Hbr, Hic, HYFLUXNBR, HYFLUXIC,
                            drb, dric, drb, dric, SQ(drb), SQ(dric));
      Uyminus2 = U_halfstep(deltaT, Ubr, Uic, VUFLUXNBR, VUFLUXIC,
                            drb, dric, drb, dric, SQ(drb), SQ(dric));
      Vyminus2 = U_halfstep(deltaT, Vbr, Vic, VYFLUXNBR, VYFLUXIC,
                            drb, dric, drb, dric, SQ(drb), SQ(dric));

      Hyfluxminus = (Hyfluxminus + HNEWYFLUXMINUS2) * HALF;
      Uyfluxminus = (Uyfluxminus + VUNEWFLUXMINUS2) * HALF;
      Vyfluxminus = (Vyfluxminus + VNEWYFLUXMINUS2) * HALF;

   }

   if(lvl < lvl_nt) {

//      V_stag = V_delta( SQ(dric), SQ(drt) );

      Hyplus2  = U_halfstep(deltaT, Hic, Htr, HYFLUXIC, HYFLUXNTR,
                            dric, drt, dric, drt, SQ(dric), SQ(drt));
      Uyplus2  = U_halfstep(deltaT, Uic, Utr, VUFLUXIC, VUFLUXNTR,
                            dric, drt, dric, drt, SQ(dric), SQ(drt));
      Vyplus2  = U_halfstep(deltaT, Vic, Vtr, VYFLUXIC, VYFLUXNTR,
                            dric, drt, dric, drt, SQ(dric), SQ(drt));

      Hyfluxplus  = (Hyfluxplus + HNEWYFLUXPLUS2) * HALF;
      Uyfluxplus  = (Uyfluxplus + VUNEWFLUXPLUS2) * HALF;
      Vyfluxplus  = (Vyfluxplus + VNEWYFLUXPLUS2) * HALF;

   }



   ///////////////////////////////////////////////////////////////////////
   ///                    XXX Flux Correction Terms XXX                ///
   ///////////////////////////////////////////////////////////////////////

   if(lvl_nl < lvl_nll) {
      Hll = (Hll + H[ ntop[nll] ]) * HALF;
      Ull = (Ull + U[ ntop[nll] ]) * HALF;
   }

   real_t Hr2 = Hr;
   real_t Ur2 = Ur;
   if(lvl < lvl_nr) {
      Hr2 = (Hr2 + Hrt) * HALF;
      Ur2 = (Ur2 + Urt) * HALF;
   }

   real_t numinusx = fabs(Uxminus/Hxminus) + sqrt(g*Hxminus);

   wminusx_H = w_corrector(deltaT, (dric+drl)*HALF, numinusx, Hic-Hl, Hl-Hll, Hr2-Hic);
   wminusx_H *= Hic - Hl;

   if(lvl < lvl_nl) {
      if(lvl_nl < level[nltl]) // Note that lvl_nl is the same as lvl_nlt
         Hll2 = (Hll2 + H[ ntop[nltl] ]) * HALF;
      wminusx_H = ((w_corrector(deltaT, (dric+drl)*HALF, fabs(Uxminus2/Hxminus2) +
                               sqrt(g*Hxminus2), Hic-Hlt, Hlt-Hll2, Hr2-Hic) *
                   (Hic - Hlt)) + wminusx_H)*QUARTER;
   }


   if(lvl_nr < lvl_nrr) {
      Hrr = (Hrr + H[ ntop[nrr] ]) * HALF;
      Urr = (Urr + U[ ntop[nrr] ]) * HALF;
   }

   real_t Hl2 = Hl;
   real_t Ul2 = Ul;
   if(lvl < lvl_nl) {
      Hl2 = (Hl2 + Hlt) * HALF;
      Ul2 = (Ul2 + Ult) * HALF;
   }

   real_t nuplusx = fabs(Uxplus/Hxplus) + sqrt(g*Hxplus);

   wplusx_H = w_corrector(deltaT, (dric+dxr)*HALF, nuplusx, Hr-Hic, Hic-Hl2, Hrr-Hr);
   wplusx_H *= Hr - Hic;

   if(lvl < lvl_nr) {
      if(lvl_nr < level[nrtr]) // Note that lvl_nr is the same as lvl_nrt
         Hrr2 = (Hrr2 + H[ ntop[nrtr] ]) * HALF;
      wplusx_H = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                               sqrt(g*Hxplus2), Hrt-Hic, Hic-Hl2, Hrr2-Hrt) *
                   (Hrt - Hic))+wplusx_H)*QUARTER;
   }


   wminusx_U = w_corrector(deltaT, (dric+dxl)*HALF, numinusx, Uic-Ul, Ul-Ull, Ur2-Uic);
   wminusx_U *= Uic - Ul;

   if(lvl < lvl_nl) {
      if(lvl_nl < level[nltl]) // Note that lvl_nl is the same as lvl_nlt
         Ull2 = (Ull2 + U[ ntop[nltl] ]) * HALF;
      wminusx_U = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                               sqrt(g*Hxminus2), Uic-Ult, Ult-Ull2, Ur2-Uic) *
                   (Uic - Ult))+wminusx_U)*QUARTER;
   }

   wplusx_U = w_corrector(deltaT, (dric+dxr)*HALF, nuplusx, Ur-Uic, Uic-Ul2, Urr-Ur);
   wplusx_U *= Ur - Uic;

   if(lvl < lvl_nr) {
      if(lvl_nr < level[nrtr]) // Note that lvl_nr is the same as lvl_nrt
         Urr2 = (Urr2 + U[ ntop[nrtr] ]) * HALF;
      wplusx_U = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                               sqrt(g*Hxplus2), Urt-Uic, Uic-Ul2, Urr2-Urt) *
                   (Urt - Uic))+wplusx_U)*QUARTER;
   }


   if(lvl_nb < lvl_nbb) {
      Hbb = (Hbb + H[ nrht[nbb] ]) * HALF;
      Vbb = (Vbb + V[ nrht[nbb] ]) * HALF;
   }

   real_t Ht2 = Ht;
   real_t Vt2 = Vt;
   if(lvl < lvl_nt) {
      Ht2 = (Ht2 + Htr) * HALF;
      Vt2 = (Vt2 + Vtr) * HALF;
   }

   real_t numinusy = fabs(Vyminus/Hyminus) + sqrt(g*Hyminus);

   wminusy_H = w_corrector(deltaT, (dric+dyb)*HALF, numinusy, Hic-Hb, Hb-Hbb, Ht2-Hic);
   wminusy_H *= Hic - Hb;

   if(lvl < lvl_nb) {
      if(lvl_nb < level[nbrb]) // Note that lvl_nb is the same as lvl_nbr
         Hbb2 = (Hbb2 + H[ nrht[nbrb] ]) * HALF;
      wminusy_H = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                               sqrt(g*Hyminus2), Hic-Hbr, Hbr-Hbb2, Ht2-Hic) *
                   (Hic - Hbr))+wminusy_H)*QUARTER;
   }


   if(lvl_nt < lvl_ntt) {
      Htt = (Htt + H[ nrht[ntt] ]) * HALF;
      Vtt = (Vtt + V[ nrht[ntt] ]) * HALF;
   }

   real_t Hb2 = Hb;
   real_t Vb2 = Vb;
   if(lvl < lvl_nb) {
      Hb2 = (Hb2 + Hbr) * HALF;
      Vb2 = (Vb2 + Vbr) * HALF;
   }

   real_t nuplusy = fabs(Vyplus/Hyplus) + sqrt(g*Hyplus);

   wplusy_H = w_corrector(deltaT, (dric+dyt)*HALF, nuplusy, Ht-Hic, Hic-Hb2, Htt-Ht);
   wplusy_H *= Ht - Hic;

   if(lvl < lvl_nt) {
      if(lvl_nt < level[ntrt]) // Note that lvl_nt is the same as lvl_ntr
         Htt2 = (Htt2 + H[ nrht[ntrt] ]) * HALF;
      wplusy_H = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                               sqrt(g*Hyplus2), Htr-Hic, Hic-Hb2, Htt2-Htr) *
                   (Htr - Hic))+wplusy_H)*QUARTER;
   }

   wminusy_V = w_corrector(deltaT, (dric+dyb)*HALF, numinusy, Vic-Vb, Vb-Vbb, Vt2-Vic);
   wminusy_V *= Vic - Vb;

   if(lvl < lvl_nb) {
      if(lvl_nb < level[nbrb]) // Note that lvl_nb is the same as lvl_nbr
         Vbb2 = (Vbb2 + V[ nrht[nbrb] ]) * HALF;
      wminusy_V = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                               sqrt(g*Hyminus2), Vic-Vbr, Vbr-Vbb2, Vt2-Vic) *
                   (Vic - Vbr))+wminusy_V)*QUARTER;
   }

   wplusy_V = w_corrector(deltaT, (dric+dyt)*HALF, nuplusy, Vt-Vic, Vic-Vb2, Vtt-Vt);
   wplusy_V *= Vt - Vic;

   if(lvl < lvl_nt) {
      if(lvl_nt < level[ntrt]) // Note that lvl_nt is the same as lvl_ntr
         Vtt2 = (Vtt2 + V[ nrht[ntrt] ]) * HALF;
      wplusy_V = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                               sqrt(g*Hyplus2), Vtr-Vic, Vic-Vb2, Vtt2-Vtr) *
                   (Vtr - Vic))+wplusy_V)*QUARTER;
   }

   ////////////////////////////////////////////////////////////////////////
   ///        ACTUAL CALCULATIONS COMBINING LAX_WENDROFF AND TVD        ///
   ////////////////////////////////////////////////////////////////////////

   Hic = U_fullstep(deltaT, dxic, Hic,
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                  - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
   Uic = U_fullstep(deltaT, dxic, Uic,
                       Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                  - wminusx_U + wplusx_U;
   Vic = U_fullstep(deltaT, dxic, Vic,
                       Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                  - wminusy_V + wplusy_V;

#if 0
   // XXX AMD seems to be fine on this bit
   // ATI fails on correction terms
   Hic += -(deltaT / dxic)*(Hxfluxplus - Hxfluxminus + Hyfluxplus - Hyfluxminus);
   //Hic += -wminusx_H + wplusx_H - wminusy_H + wplusy_H;

   Uic += -(deltaT / dxic)*(Uxfluxplus - Uxfluxminus + Uyfluxplus - Uyfluxminus);
   //Uic += - wminusx_U + wplusx_U;

   Vic += -(deltaT / dxic)*(Vxfluxplus - Vxfluxminus + Vyfluxplus - Vyfluxminus);
   //Vic += - wminusy_V + wplusy_V;
#endif

// XXX How does this barrier effect the finite difference calculation? XXX
   barrier(CLK_LOCAL_MEM_FENCE);
   

    //  Write values for the main cell back to global memory.
    H_new[giX] = Hic;
    U_new[giX] = Uic;
    V_new[giX] = Vic;

//////////////////////////////////////////////////////////////////
////////////////////          END           //////////////////////
////////////////////   calc_one_cycle_cl    //////////////////////
//////////////////////////////////////////////////////////////////

}

__kernel void calc_finite_difference_via_faces_face_comps_cl(
            __global          int       *nface,                     // 0 Number array of faces
                        const int       levmx,                      // 1 Maximum level
            __global    const state_t   *H,                         // 2
            __global    const state_t   *U,                         // 3
            __global    const state_t   *V,                         // 4
            __global    const uchar_t   *level,                     // 5 Array of level information
                        const real_t    deltaT,                     // 6 Size of time step
            __global    const real_t    *lev_deltax,                // 7
            __global    const real_t    *lev_deltay,                // 8
            __local           state4_t  *tile,                      // 9Tile size in state4_t
            __local           int8      *itile,                     // 10 Tile size in int8
            __local           int8      *xface,                     // 11 xFace size in int8
            __local           int8      *yface,                     // 12 yFace size in int8 
            __global    const int       *map_xface2cell_lower,      // 13 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 14 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 15 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 16 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 17 
            __global    const int       *map_xcell2face_left2,      // 18
            __global    const int       *map_xcell2face_right1,     // 19 
            __global    const int       *map_xcell2face_right2,     // 20 
            __global    const int       *map_ycell2face_bot1,       // 21 
            __global    const int       *map_ycell2face_bot2,       // 22 
            __global    const int       *map_ycell2face_top1,       // 23 
            __global    const int       *map_ycell2face_top2,       // 24 
            __global          state_t   *HxFlux,                    // 25
            __global          state_t   *UxFlux,                    // 26
            __global          state_t   *VxFlux,                    // 27
            __global          state_t   *HyFlux,                    // 28
            __global          state_t   *UyFlux,                    // 29
            __global          state_t   *VyFlux,                    // 30
            __global          state_t   *Wx_H,                      // 31
            __global          state_t   *Wx_U,                      // 32
            __global          state_t   *Wy_H,                      // 33
            __global          state_t   *Wy_V,                      // 34
            __global          int       *nlft,                      // 35
            __global          int       *nrht,                      // 36
            __global          int       *nbot,                      // 37
            __global          int       *ntop) {                    // 38

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    // Ensure the executing thread is not extraneous
    if (giX >= max(nface[0], nface[1]))
        return;


  //if (giX < ncells) { // only the workers equal to the number of cells
  //  setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);
  //}
  //if (giX < nface[0]) {
  //  setup_xface(xface, nface[0], map_xface2cell_lower, map_xface2cell_upper);
  //}
  //if (giX < nface[1]) {
  //  setup_yface(yface, nface[1], map_yface2cell_lower, map_yface2cell_upper);
  //}

//if (0) {
  real_t g = 9.80;
  real_t ghalf = 0.5*g;
  if (giX < nface[0]) {
      int iface = giX;
      int cell_lower = map_xface2cell_lower[iface];
      int cell_upper = map_xface2cell_upper[iface];
      real_t Hx, Ux, Vx;
      if (level[cell_lower] == level[cell_upper]) {

         // set the two faces
         int fl = map_xcell2face_left1[cell_lower];
         int fr = map_xcell2face_right1[cell_upper];
         // set the two cells away
         int nll = map_xface2cell_lower[fl];
         int nrr = map_xface2cell_upper[fr];
 
         int lev = level[cell_lower];
         real_t dxic = lev_deltax[lev];
         real_t Cxhalf = 0.5*deltaT/dxic;
 
         real_t Hic = H[cell_lower];
         real_t Hr  = H[cell_upper];
         real_t Hl  = H[nll];
         real_t Hrr = H[nrr];
         real_t Uic = U[cell_lower];
         real_t Ur  = U[cell_upper];
         real_t Ul  = U[nll];
         real_t Urr = U[nrr];
 
         Hx=HALF*(H[cell_upper]+H[cell_lower]) - Cxhalf*( HXFLUX(cell_upper)-HXFLUX(cell_lower) );
         Ux=HALF*(U[cell_upper]+U[cell_lower]) - Cxhalf*( UXFLUX(cell_upper)-UXFLUX(cell_lower) );
         Vx=HALF*(V[cell_upper]+V[cell_lower]) - Cxhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );
 
         real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);
 
         Wx_H[iface] = w_corrector(deltaT, dxic, U_eigen, Hr-Hic, Hic-Hl, Hrr-Hr) * (Hr - Hic);
 
         Wx_U[iface] = w_corrector(deltaT, dxic, U_eigen, Ur-Uic, Uic-Ul, Urr-Ur) * (Ur - Uic);

         HxFlux[iface] = HXFLUXFACE;
         UxFlux[iface] = UXFLUXFACE;
         VxFlux[iface] = VXFLUXFACE;
      } else {

         real_t dx_lower = lev_deltax[level[cell_lower]];
         real_t dx_upper = lev_deltax[level[cell_upper]];

         real_t FA_lower = dx_lower;
         real_t FA_upper = dx_upper;
         real_t FA_lolim = FA_lower*min(ONE, FA_upper/FA_lower);
         real_t FA_uplim = FA_upper*min(ONE, FA_lower/FA_upper);

         real_t CV_lower = SQ(dx_lower);
         real_t CV_upper = SQ(dx_upper);
         real_t CV_lolim = CV_lower*min(HALF, CV_upper/CV_lower);
         real_t CV_uplim = CV_upper*min(HALF, CV_lower/CV_upper);

         // Weighted half-step calculation
         //
         // (dx_lower*H[cell_upper]+dx_upper*H[cell_lower])
         // -----------------------------------------------   -
         //             (dx_lower+dx_upper)
         //
         //                ( (FA_uplim*HXFLUX(cell_upper))-(FA_lolim*HXFLUX(cell_lower)) )
         // 0.5*deltaT  *  ----------------------------------------------------------------
         //                                    (CV_uplim+CV_lolim)
         //

         Hx=(dx_lower*H[cell_upper]+dx_upper*H[cell_lower])/(dx_lower+dx_upper) -
                   HALF*deltaT*( (FA_uplim*HXFLUX(cell_upper))-(FA_lolim*HXFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Ux=(dx_lower*U[cell_upper]+dx_upper*U[cell_lower])/(dx_lower+dx_upper) -
                   HALF*deltaT*( (FA_uplim*UXFLUX(cell_upper))-(FA_lolim*UXFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Vx=(dx_lower*V[cell_upper]+dx_upper*V[cell_lower])/(dx_lower+dx_upper) -
                   HALF*deltaT*( (FA_uplim*UVFLUX(cell_upper))-(FA_lolim*UVFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim); 

         // set the two faces
         int fl = map_xcell2face_left1[cell_lower];
         int fr = map_xcell2face_right1[cell_upper];
         // set the two cells away
         int nll = map_xface2cell_lower[fl];
         int nrr = map_xface2cell_upper[fr];

         uchar_t lev = level[cell_lower];
         uchar_t levr = level[cell_upper];
         real_t dxic = lev_deltax[lev]; 
         real_t dxr = lev_deltax[levr];

         real_t Hic = H[cell_lower];
         real_t Hr  = H[cell_upper];
         real_t Hl  = H[nll];
         real_t Hrr = H[nrr];
         real_t Uic = U[cell_lower];
         real_t Ur  = U[cell_upper];
         real_t Ul  = U[nll];
         real_t Urr = U[nrr];

         real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);
         real_t dx_avg = (dxic+dxr)*HALF;

         if(level[cell_upper] < level[nrr]) {
            //Hrr = (Hrr + H[map_yface2cell_upper[map_ycell2face_top1[nrr]]]) * HALF;
            //Urr = (Urr + U[map_yface2cell_upper[map_ycell2face_top1[nrr]]]) * HALF;
            Hrr = (Hrr + H[ntop[nrr]]) * HALF;
            Urr = (Urr + U[ntop[nrr]]) * HALF;
         }

         real_t Hl2 = Hl;
         real_t Ul2 = Ul;
         if(lev < level[nll]) {
            //Hl2 = (Hl2 + H[map_yface2cell_upper[map_ycell2face_top1[nll]]]) * HALF;
            //Ul2 = (Ul2 + U[map_yface2cell_upper[map_ycell2face_top1[nll]]]) * HALF;
            Hl2 = (Hl2 + H[ntop[nll]]) * HALF;
            Ul2 = (Ul2 + U[ntop[nll]]) * HALF;
         }

         Wx_H[iface] = w_corrector(deltaT, dx_avg, U_eigen, Hr-Hic, Hic-Hl2, Hrr-Hr) * (Hr - Hic);
         Wx_U[iface] = w_corrector(deltaT, dx_avg, U_eigen, Ur-Uic, Uic-Ul2, Urr-Ur) * (Ur - Uic);
 
         HxFlux[iface] = HXFLUXFACE;
         UxFlux[iface] = UXFLUXFACE;
         VxFlux[iface] = VXFLUXFACE;
      }
   }

   if (giX < nface[1]) {
      int iface = giX;
      int cell_lower = map_yface2cell_lower[iface];
      int cell_upper = map_yface2cell_upper[iface];
      real_t Hy, Uy, Vy;
      if (level[cell_lower] == level[cell_upper]) {

         // set the two faces
         int fb = map_ycell2face_bot1[cell_lower];
         int ft = map_ycell2face_top1[cell_upper];
         // set the two cells away
         int nbb = map_yface2cell_lower[fb];
         int ntt = map_yface2cell_upper[ft];
	
         int lev = level[cell_lower];
         real_t dyic    = lev_deltay[lev];
         real_t Cyhalf = 0.5*deltaT/lev_deltay[lev];

         real_t Hic = H[cell_lower];
         real_t Ht  = H[cell_upper];
         real_t Hb  = H[nbb];
         real_t Htt = H[ntt];
         real_t Vic = V[cell_lower];
         real_t Vt  = V[cell_upper];
         real_t Vb  = V[nbb];
         real_t Vtt = V[ntt];

         Hy=HALF*(H[cell_upper]+H[cell_lower]) - Cyhalf*( HYFLUX(cell_upper)-HYFLUX(cell_lower) );
         Uy=HALF*(U[cell_upper]+U[cell_lower]) - Cyhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );
         Vy=HALF*(V[cell_upper]+V[cell_lower]) - Cyhalf*( VYFLUX(cell_upper)-VYFLUX(cell_lower) );

         real_t U_eigen = fabs(Vy/Hy) + sqrt(g*Hy);
 
         Wy_H[iface] = w_corrector(deltaT, dyic, U_eigen, Ht-Hic, Hic-Hb, Htt-Ht) * (Ht - Hic);
         Wy_V[iface] = w_corrector(deltaT, dyic, U_eigen, Vt-Vic, Vic-Vb, Vtt-Vt) * (Vt - Vic);

         HyFlux[iface] = HYFLUXFACE;
         UyFlux[iface] = UYFLUXFACE;
         VyFlux[iface] = VYFLUXFACE;
      } else {
         real_t dy_lower = lev_deltay[level[cell_lower]];
         real_t dy_upper = lev_deltay[level[cell_upper]];

         real_t FA_lower = dy_lower;
         real_t FA_upper = dy_upper;
         real_t FA_lolim = FA_lower*min(ONE, FA_upper/FA_lower);
         real_t FA_uplim = FA_upper*min(ONE, FA_lower/FA_upper);

         real_t CV_lower = SQ(dy_lower);
         real_t CV_upper = SQ(dy_upper);
         real_t CV_lolim = CV_lower*min(HALF, CV_upper/CV_lower);
         real_t CV_uplim = CV_upper*min(HALF, CV_lower/CV_upper);

         // Weighted half-step calculation
         //
         // (dy_lower*H[cell_upper]+dy_upper*H[cell_lower])
         // -----------------------------------------------   -
         //             (dy_lower+dy_upper)
         //
         //                ( (FA_uplim*HYFLUX(cell_upper))-(FA_lolim*HYFLUX(cell_lower)) )
         // 0.5*deltaT  *  ----------------------------------------------------------------
         //                                    (CV_uplim+CV_lolim)
         //

         Hy=(dy_lower*H[cell_upper]+dy_upper*H[cell_lower])/(dy_lower+dy_upper) -
                   HALF*deltaT*( (FA_uplim*HYFLUX(cell_upper))-(FA_lolim*HYFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Uy=(dy_lower*U[cell_upper]+dy_upper*U[cell_lower])/(dy_lower+dy_upper) -
                   HALF*deltaT*( (FA_uplim*UVFLUX(cell_upper))-(FA_lolim*UVFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Vy=(dy_lower*V[cell_upper]+dy_upper*V[cell_lower])/(dy_lower+dy_upper) -
                   HALF*deltaT*( (FA_uplim*VYFLUX(cell_upper))-(FA_lolim*VYFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);

         // set the two faces
         int fb = map_ycell2face_bot1[cell_lower];
         int ft = map_ycell2face_top1[cell_upper];
         // set the two cells away
         int nbb = map_yface2cell_lower[fb];
         int ntt = map_yface2cell_upper[ft];
 
         uchar_t lev = level[cell_lower];
         uchar_t levt = level[cell_upper];
         real_t dyic = lev_deltay[lev];
         real_t dyt = lev_deltay[levt];

         real_t Hic = H[cell_lower];
         real_t Ht  = H[cell_upper];
         real_t Hb  = H[nbb];
         real_t Htt = H[ntt];
         real_t Vic = V[cell_lower];
         real_t Vt  = V[cell_upper];
         real_t Vb  = V[nbb];
         real_t Vtt = V[ntt];

         real_t V_eigen = fabs(Vy/Hy) + sqrt(g*Hy);
         real_t dy_avg = (dyic+dyt)*HALF; 

         if(level[cell_upper] < level[ntt]) {
            //Htt = (Htt + H[map_xface2cell_upper[map_xcell2face_right1[ntt]]]) * HALF;
            //Vtt = (Vtt + V[map_xface2cell_upper[map_xcell2face_right1[ntt]]]) * HALF;
            Htt = (Htt + H[nrht[ntt]]) * HALF;
            Vtt = (Vtt + V[nrht[ntt]]) * HALF;
         } 

         real_t Hb2 = Hb;
         real_t Vb2 = Vb;
         if(lev < level[nbb]) {
            //Hb2 = (Hb2 + H[map_xface2cell_upper[map_xcell2face_right1[nbb]]]) * HALF;
            //Vb2 = (Vb2 + V[map_xface2cell_upper[map_xcell2face_right1[nbb]]]) * HALF;
            Hb2 = (Hb2 + H[nrht[nbb]]) * HALF;
            Vb2 = (Vb2 + V[nrht[nbb]]) * HALF;
         }

         Wy_H[iface] = w_corrector(deltaT, dy_avg, V_eigen, Ht-Hic, Hic-Hb2, Htt-Ht) * (Ht - Hic);
         Wy_V[iface] = w_corrector(deltaT, dy_avg, V_eigen, Vt-Vic, Vic-Vb2, Vtt-Vt) * (Vt - Vic);
 
         HyFlux[iface] = HYFLUXFACE;
         UyFlux[iface] = UYFLUXFACE;
         VyFlux[iface] = VYFLUXFACE;
      }
   }
//}
}


__kernel void calc_finite_difference_via_faces_cell_comps_cl (
            __global    const state_t   *H,                         // 0
            __global    const state_t   *U,                         // 1
            __global    const state_t   *V,                         // 2
            __global    const uchar_t   *level,                     // 3 Array of level information
            __global    const int       *map_xface2cell_lower,      // 4 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 5 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 6 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 7 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 8
            __global    const int       *map_xcell2face_left2,      // 9 
            __global    const int       *map_xcell2face_right1,     // 10 
            __global    const int       *map_xcell2face_right2,     // 11 
            __global    const int       *map_ycell2face_bot1,       // 12 
            __global    const int       *map_ycell2face_bot2,       // 13 
            __global    const int       *map_ycell2face_top1,       // 14 
            __global    const int       *map_ycell2face_top2,       // 15 
            __global          state_t   *HxFlux,                    // 16
            __global          state_t   *UxFlux,                    // 17
            __global          state_t   *VxFlux,                    // 18
            __global          state_t   *HyFlux,                    // 19
            __global          state_t   *UyFlux,                    // 20
            __global          state_t   *VyFlux,                    // 21
            __global          state_t   *Wx_H,                      // 22
            __global          state_t   *Wx_U,                      // 23
            __global          state_t   *Wy_H,                      // 24
            __global          state_t   *Wy_V,                      // 25
            __global          state_t   *H_new,                     // 26
            __global          state_t   *U_new,                     // 27
            __global          state_t   *V_new,                     // 28
                        const int       ncells,                     // 29  Total number of cells
                        const real_t    deltaT,                     // 30 Size of time step
            __global    const real_t    *lev_deltax,                // 31
            __global    const real_t    *lev_deltay) {              // 32

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

   /////////////////////////////////////////////
   /// Set local tile & apply boundary conds ///
   /////////////////////////////////////////////

   //setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);

   //barrier(CLK_LOCAL_MEM_FENCE);

    // Ensure the executing thread is not extraneous
    if (giX >= ncells)
        return;

      int ic = giX;
      real_t dxic    = lev_deltax[level[ic]];
      // set the four faces
      int fl = map_xcell2face_left1[ic];
      int fr = map_xcell2face_right1[ic];
      int fb = map_ycell2face_bot1[ic];
      int ft = map_ycell2face_top1[ic];
      int fl2 = map_xcell2face_left2[ic];
      int fr2 = map_xcell2face_right2[ic];
      int fb2 = map_ycell2face_bot2[ic];
      int ft2 = map_ycell2face_top2[ic];

      // set the four neighboring cells
      int nl = map_xface2cell_lower[fl];
      int nr = map_xface2cell_upper[fr];
      int nb = map_yface2cell_lower[fb];
      int nt = map_yface2cell_upper[ft];

      if (nb == ic  || nt == ic || nl == ic || nr == ic) 
          return;

      real_t Hic     = H[ic];
      real_t Uic     = U[ic];
      real_t Vic     = V[ic];

      real_t Hxfluxminus = HxFlux[fl];
      real_t Uxfluxminus = UxFlux[fl];
      real_t Vxfluxminus = VxFlux[fl];

      real_t Hxfluxplus  = HxFlux[fr];
      real_t Uxfluxplus  = UxFlux[fr];
      real_t Vxfluxplus  = VxFlux[fr];

      real_t Hyfluxminus = HyFlux[fb];
      real_t Uyfluxminus = UyFlux[fb];
      real_t Vyfluxminus = VyFlux[fb];

      real_t Hyfluxplus  = HyFlux[ft];
      real_t Uyfluxplus  = UyFlux[ft];
      real_t Vyfluxplus  = VyFlux[ft];

      real_t wminusx_H = Wx_H[fl];
      real_t wminusx_U = Wx_U[fl];

      real_t wplusx_H = Wx_H[fr];
      real_t wplusx_U = Wx_U[fr];

      real_t wminusy_H = Wy_H[fb];
      real_t wminusy_V = Wy_V[fb];

      real_t wplusy_H = Wy_H[ft];
      real_t wplusy_V = Wy_V[ft];

      if (level[ic] < level[nl]) {
         Hxfluxminus = (Hxfluxminus + HxFlux[fl2]) * HALF;
         Uxfluxminus = (Uxfluxminus + UxFlux[fl2]) * HALF;
         Vxfluxminus = (Vxfluxminus + VxFlux[fl2]) * HALF;
         wminusx_H = (wminusx_H + Wx_H[fl2]) * HALF * HALF;
         wminusx_U = (wminusx_U + Wx_U[fl2]) * HALF * HALF;
      }
   
      if (level[ic] < level[nr]) {
         Hxfluxplus = (Hxfluxplus + HxFlux[fr2]) * HALF;
         Uxfluxplus = (Uxfluxplus + UxFlux[fr2]) * HALF;
         Vxfluxplus = (Vxfluxplus + VxFlux[fr2]) * HALF;
         wplusx_H = (wplusx_H + Wx_H[fr2]) * HALF * HALF;
         wplusx_U = (wplusx_U + Wx_U[fr2]) * HALF * HALF;
      }
   
      if (level[ic] < level[nb]) {
         Hyfluxminus = (Hyfluxminus + HyFlux[fb2]) * HALF;
         Uyfluxminus = (Uyfluxminus + UyFlux[fb2]) * HALF;
         Vyfluxminus = (Vyfluxminus + VyFlux[fb2]) * HALF;
         wminusy_H = (wminusy_H + Wy_H[fb2]) * HALF * HALF;
         wminusy_V = (wminusy_V + Wy_V[fb2]) * HALF * HALF;
      }
   
      if (level[ic] < level[nt]) {
         Hyfluxplus = (Hyfluxplus + HyFlux[ft2]) * HALF;
         Uyfluxplus = (Uyfluxplus + UyFlux[ft2]) * HALF;
         Vyfluxplus = (Vyfluxplus + VyFlux[ft2]) * HALF;
         wplusy_H = (wplusy_H + Wy_H[ft2]) * HALF * HALF;
         wplusy_V = (wplusy_V + Wy_V[ft2]) * HALF * HALF;
      }

      barrier(CLK_LOCAL_MEM_FENCE);

      H_new[ic] = U_fullstep(deltaT, dxic, Hic,
                      Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                 - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
      U_new[ic] = U_fullstep(deltaT, dxic, Uic,
                      Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                 - wminusx_U + wplusx_U;
      V_new[ic] = U_fullstep(deltaT, dxic, Vic,
                      Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                 - wminusy_V + wplusy_V;
}

__kernel void calc_finite_difference_in_place_cell_comps_cl (
                        const int       ncells,                     // 0 Number of cells (not including phantom)
            __global    const int       *nfaces,                    // 1 Number of x faces
                        const int       levmx,                      // 2 Maximum level
            __global    const state_t   *H,                         // 3
            __global    const state_t   *U,                         // 4
            __global    const state_t   *V,                         // 5
            __global    const uchar_t   *level,                     // 6 Array of level information
                        const real_t    deltaT,                     // 7 Size of time step
            __global    const real_t    *lev_dx,                    // 8
            __global    const real_t    *lev_dy,                    // 9
            __local           state4_t  *tile,                      // 10 Tile size in state4_t
            __local           int8      *itile,                     // 11 Tile size in int8
            __local           int8      *xface,                     // 12 xFace size in int8
            __local           int8      *yface,                     // 13 yFace size in int8 
            __global    const int       *map_xface2cell_lower,      // 14 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 15 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 16 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 17 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 18 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 19 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 20 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 21 A cell's top primary face 
            __global          state_t   *Hxfluxplus,                // 22
            __global          state_t   *Hxfluxminus,               // 23
            __global          state_t   *Uxfluxplus,                // 24
            __global          state_t   *Uxfluxminus,               // 25
            __global          state_t   *Vxfluxplus,                // 26
            __global          state_t   *Vxfluxminus,               // 27
            __global          state_t   *Hyfluxplus,                // 28
            __global          state_t   *Hyfluxminus,               // 29
            __global          state_t   *Uyfluxplus,                // 30
            __global          state_t   *Uyfluxminus,               // 31
            __global          state_t   *Vyfluxplus,                // 32
            __global          state_t   *Vyfluxminus,               // 33
            __global          state_t   *wplusx_H,                  // 34
            __global          state_t   *wminusx_H,                 // 35
            __global          state_t   *wplusx_U,                  // 36
            __global          state_t   *wminusx_U,                 // 37
            __global          state_t   *wplusy_H,                  // 38
            __global          state_t   *wminusy_H,                 // 39
            __global          state_t   *wplusy_V,                  // 40
            __global          state_t   *wminusy_V) {               // 41

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    
    if (giX >= ncells) 
        return;

    real_t   g     = 9.80; 
    real_t   ghalf = HALF*g;

    int ic = giX;
    int lev = level[ic];
    real_t dxic    = lev_dx[lev];
    real_t dyic    = lev_dy[lev];
    real_t Cxhalf = 0.5*deltaT/dxic;
    real_t Cyhalf = 0.5*deltaT/dyic;

    int fl = map_xcell2face_left1[ic];
    int fr = map_xcell2face_right1[ic];
    int fb = map_ycell2face_bot1[ic];
    int ft = map_ycell2face_top1[ic];

    int nl = map_xface2cell_lower[fl];
    int nr = map_xface2cell_upper[fr];
    int nb = map_yface2cell_lower[fb];
    int nt = map_yface2cell_upper[ft];

    //if (ic == nl || ic == nr || ic == nb || ic == nt) return;

    real_t Hic     = H[ic];
    real_t Uic     = U[ic];
    real_t Vic     = V[ic];

    int nll     = map_xface2cell_lower[map_xcell2face_left1[nl]];
    real_t Hl      = H[nl];
    real_t Ul      = U[nl];
    real_t Vl      = V[nl];

    int nrr     = map_xface2cell_upper[map_xcell2face_right1[nr]];
    real_t Hr      = H[nr];
    real_t Ur      = U[nr];
    real_t Vr      = V[nr];

    int ntt     = map_yface2cell_upper[map_ycell2face_top1[nt]];
    real_t Ht      = H[nt];
    real_t Ut      = U[nt];
    real_t Vt      = V[nt];

    int nbb     = map_yface2cell_lower[map_ycell2face_bot1[nb]];
    real_t Hb      = H[nb];
    real_t Ub      = U[nb];
    real_t Vb      = V[nb];

    real_t Hll     = H[nll];
    real_t Ull     = U[nll];

    real_t Hrr     = H[nrr];
    real_t Urr     = U[nrr];

    real_t Htt     = H[ntt];
    real_t Vtt     = V[ntt];

    real_t Hbb     = H[nbb];
    real_t Vbb     = V[nbb];

    real_t Hxminus = HALF*(Hic+Hl)-Cxhalf*(HXFLUXIC-HXFLUXNL);
    real_t Uxminus = HALF*(Uic+Ul)-Cxhalf*(UXFLUXIC-UXFLUXNL);
    real_t Vxminus = HALF*(Vic+Vl)-Cxhalf*(UVFLUXIC-UVFLUXNL);

    real_t Hxplus = HALF*(Hr+Hic)-Cxhalf*(HXFLUXNR-HXFLUXIC);
    real_t Uxplus = HALF*(Ur+Uic)-Cxhalf*(UXFLUXNR-UXFLUXIC);
    real_t Vxplus = HALF*(Vr+Vic)-Cxhalf*(UVFLUXNR-UVFLUXIC);

    real_t Hyminus = HALF*(Hic+Hb)-Cyhalf*(HYFLUXIC-HYFLUXNB);
    real_t Uyminus = HALF*(Uic+Ub)-Cyhalf*(VUFLUXIC-VUFLUXNB);
    real_t Vyminus = HALF*(Vic+Vb)-Cyhalf*(VYFLUXIC-VYFLUXNB);

    real_t Hyplus = HALF*(Ht+Hic)-Cyhalf*(HYFLUXNT-HYFLUXIC);
    real_t Uyplus = HALF*(Ut+Uic)-Cyhalf*(VUFLUXNT-VUFLUXIC);
    real_t Vyplus = HALF*(Vt+Vic)-Cyhalf*(VYFLUXNT-VYFLUXIC);

    Hxfluxminus[ic] = HNEWXFLUXMINUS;
    Uxfluxminus[ic] = UNEWXFLUXMINUS;
    Vxfluxminus[ic] = UVNEWFLUXMINUS;

    Hxfluxplus[ic]  = HNEWXFLUXPLUS;
    Uxfluxplus[ic]  = UNEWXFLUXPLUS;
    Vxfluxplus[ic]  = UVNEWFLUXPLUS;

    Hyfluxminus[ic] = HNEWYFLUXMINUS;
    Uyfluxminus[ic] = VUNEWFLUXMINUS;
    Vyfluxminus[ic] = VNEWYFLUXMINUS;

    Hyfluxplus[ic]  = HNEWYFLUXPLUS;
    Uyfluxplus[ic]  = VUNEWFLUXPLUS;
    Vyfluxplus[ic]  = VNEWYFLUXPLUS;

    real_t U_eigen = fabs(Uxminus/Hxminus) + sqrt(g*Hxminus);
    wminusx_H[ic] = w_corrector(deltaT, dxic, U_eigen, Hic-Hl, Hl-Hll, Hr-Hic) * (Hic - Hl);
    wminusx_U[ic] = w_corrector(deltaT, dxic, U_eigen, Uic-Ul, Ul-Ull, Ur-Uic) * (Uic - Ul);

    U_eigen = fabs(Uxplus/Hxplus) + sqrt(g*Hxplus);
    wplusx_H[ic] = w_corrector(deltaT, dxic, U_eigen, Hr-Hic, Hic-Hl, Hrr-Hr) * (Hr - Hic);
    wplusx_U[ic] = w_corrector(deltaT, dxic, U_eigen, Ur-Uic, Uic-Ul, Urr-Ur) * (Ur - Uic);

    U_eigen = fabs(Vyminus/Hyminus) + sqrt(g*Hyminus);
    wminusy_H[ic] = w_corrector(deltaT, dyic, U_eigen, Hic-Hb, Hb-Hbb, Ht-Hic) * (Hic - Hb);
    wminusy_V[ic] = w_corrector(deltaT, dyic, U_eigen, Vic-Vb, Vb-Vbb, Vt-Vic) * (Vic - Vb);

    U_eigen = fabs(Vyplus/Hyplus) + sqrt(g*Hyplus);
    wplusy_H[ic] = w_corrector(deltaT, dyic, U_eigen, Ht-Hic, Hic-Hb, Htt-Ht) * (Ht - Hic);
    wplusy_V[ic] = w_corrector(deltaT, dyic, U_eigen, Vt-Vic, Vic-Vb, Vtt-Vt) * (Vt - Vic);
}

__kernel void calc_finite_difference_in_place_fixup_cl(
                        const int        nxfixup,                   // 0
                        const int        nyfixup,                   // 1
            __global    const int       *xrecvCIdx,                 // 2
            __global    const int       *xplusCell2Idx,             // 3
            __global    const int       *xminusCell2Idx,            // 4
            __global    const int       *xsendIdx1,                 // 5
            __global    const int       *xsendIdx2,                 // 6
            __global    const int       *yrecvCIdx,                 // 7
            __global    const int       *yplusCell2Idx,             // 8
            __global    const int       *yminusCell2Idx,            // 9
            __global    const int       *ysendIdx1,                 // 10
            __global    const int       *ysendIdx2,                 // 11
            __global    const int       *map_xface2cell_lower,      // 12
            __global    const int       *map_xface2cell_upper,      // 13
            __global    const int       *map_yface2cell_lower,      // 14 
            __global    const int       *map_yface2cell_upper,      // 15 
            __global          state_t   *Hxfluxplus,                // 16
            __global          state_t   *Hxfluxminus,               // 17
            __global          state_t   *Uxfluxplus,                // 18
            __global          state_t   *Uxfluxminus,               // 19
            __global          state_t   *Vxfluxplus,                // 20
            __global          state_t   *Vxfluxminus,               // 21
            __global          state_t   *Hyfluxplus,                // 22
            __global          state_t   *Hyfluxminus,               // 23
            __global          state_t   *Uyfluxplus,                // 24
            __global          state_t   *Uyfluxminus,               // 25
            __global          state_t   *Vyfluxplus,                // 26
            __global          state_t   *Vyfluxminus,               // 27
            __global          state_t   *wplusx_H,                  // 28
            __global          state_t   *wminusx_H,                 // 29
            __global          state_t   *wplusx_U,                  // 30
            __global          state_t   *wminusx_U,                 // 31
            __global          state_t   *wplusy_H,                  // 32
            __global          state_t   *wminusy_H,                 // 33
            __global          state_t   *wplusy_V,                  // 34
            __global          state_t   *wminusy_V) {               // 35

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    
    if (giX >= max(nxfixup, nyfixup)) 
        return;

    int ifix = giX;
    if (giX < nxfixup) {
      int ic = xrecvCIdx[ifix];

      if (xplusCell2Idx[ic] > -1) {
         int ifixup = xplusCell2Idx[ic];

         int ns1 = map_xface2cell_upper[xsendIdx1[ifixup]];
         int ns2 = map_xface2cell_upper[xsendIdx2[ifixup]];

         Hxfluxplus[ic] = (Hxfluxminus[ns1] + Hxfluxminus[ns2]) * HALF;
         Uxfluxplus[ic] = (Uxfluxminus[ns1] + Uxfluxminus[ns2]) * HALF;
         Vxfluxplus[ic] = (Vxfluxminus[ns1] + Vxfluxminus[ns2]) * HALF;
         wplusx_H[ic] = (wminusx_H[ns1] + wminusx_H[ns2]) * 0.25;
         wplusx_U[ic] = (wminusx_U[ns1] + wminusx_U[ns2]) * 0.25;
      }

      if (xminusCell2Idx[ic] > -1) {
         int ifixup = xminusCell2Idx[ic];

         int ns1 = map_xface2cell_lower[xsendIdx1[ifixup]];
         int ns2 = map_xface2cell_lower[xsendIdx2[ifixup]];

         Hxfluxminus[ic] = (Hxfluxplus[ns1] + Hxfluxplus[ns2]) * HALF;
         Uxfluxminus[ic] = (Uxfluxplus[ns1] + Uxfluxplus[ns2]) * HALF;
         Vxfluxminus[ic] = (Vxfluxplus[ns1] + Vxfluxplus[ns2]) * HALF;
         wminusx_H[ic] = (wplusx_H[ns1] + wplusx_H[ns2]) * 0.25;
         wminusx_U[ic] = (wplusx_U[ns1] + wplusx_U[ns2]) * 0.25;
      }
    }

    if (giX < nyfixup) {
      int ic = yrecvCIdx[ifix];

      if (yplusCell2Idx[ic] > -1) {
         int ifixup = yplusCell2Idx[ic];

         int ns1 = map_yface2cell_upper[ysendIdx1[ifixup]];
         int ns2 = map_yface2cell_upper[ysendIdx2[ifixup]];

         Hyfluxplus[ic] = (Hyfluxminus[ns1] + Hyfluxminus[ns2]) * HALF;
         Uyfluxplus[ic] = (Uyfluxminus[ns1] + Uyfluxminus[ns2]) * HALF;
         Vyfluxplus[ic] = (Vyfluxminus[ns1] + Vyfluxminus[ns2]) * HALF;
         wplusy_H[ic] = (wminusy_H[ns1] + wminusy_H[ns2]) * 0.25;
         wplusy_V[ic] = (wminusy_V[ns1] + wminusy_V[ns2]) * 0.25;
      }

      if (yminusCell2Idx[ic] > -1) {
         int ifixup = yminusCell2Idx[ic];

         int ns1 = map_yface2cell_lower[ysendIdx1[ifixup]];
         int ns2 = map_yface2cell_lower[ysendIdx2[ifixup]];

         Hyfluxminus[ic] = (Hyfluxplus[ns1] + Hyfluxplus[ns2]) * HALF;
         Uyfluxminus[ic] = (Uyfluxplus[ns1] + Uyfluxplus[ns2]) * HALF;
         Vyfluxminus[ic] = (Vyfluxplus[ns1] + Vyfluxplus[ns2]) * HALF;
         wminusy_H[ic] = (wplusy_H[ns1] + wplusy_H[ns2]) * 0.25;
         wminusy_V[ic] = (wplusy_V[ns1] + wplusy_V[ns2]) * 0.25;
      }
    
    }
}

__kernel void calc_finite_difference_in_place_fill_new_cl(
                        const int       ncells,                     // 0 Number of cells (not including phantom)
                        const real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global    const state_t   *Hxfluxplus,                // 4
            __global    const state_t   *Hxfluxminus,               // 5
            __global    const state_t   *Uxfluxplus,                // 6
            __global    const state_t   *Uxfluxminus,               // 7
            __global    const state_t   *Vxfluxplus,                // 8
            __global    const state_t   *Vxfluxminus,               // 9
            __global    const state_t   *Hyfluxplus,                // 10
            __global    const state_t   *Hyfluxminus,               // 11
            __global    const state_t   *Uyfluxplus,                // 12
            __global    const state_t   *Uyfluxminus,               // 13
            __global    const state_t   *Vyfluxplus,                // 14
            __global    const state_t   *Vyfluxminus,               // 15
            __global    const state_t   *wplusx_H,                  // 16
            __global    const state_t   *wminusx_H,                 // 17
            __global    const state_t   *wplusx_U,                  // 18
            __global    const state_t   *wminusx_U,                 // 19
            __global    const state_t   *wplusy_H,                  // 20
            __global    const state_t   *wminusy_H,                 // 21
            __global    const state_t   *wplusy_V,                  // 22
            __global    const state_t   *wminusy_V,                 // 23
            __global    const int       *level,                     // 24
            __global    const int       *map_xface2cell_lower,      // 25 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 26 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 27 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 28 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 29 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 30 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 31 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 32 A cell's top primary face 
            __global    const state_t   *H,                         // 33
            __global    const state_t   *U,                         // 34
            __global    const state_t   *V,                         // 35
            __global          state_t   *H_new,                     // 36
            __global          state_t   *U_new,                     // 37
            __global          state_t   *V_new) {                   // 38

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    
    if (giX >= ncells) 
        return;

      int ic = giX;
      int lev = level[ic];
      real_t dxic    = lev_dx[lev];
      real_t dyic    = lev_dy[lev];

      int fl = map_xcell2face_left1[ic];
      int fr = map_xcell2face_right1[ic];
      int fb = map_ycell2face_bot1[ic];
      int ft = map_ycell2face_top1[ic];

      int nl = map_xface2cell_lower[fl];
      int nr = map_xface2cell_upper[fr];
      int nb = map_yface2cell_lower[fb];
      int nt = map_yface2cell_upper[ft];

      //if (ic == nl || ic == nr || ic == nb || ic == nt) return;

      H_new[ic] = U_fullstep(deltaT, dxic, H[ic],
                       Hxfluxplus[ic], Hxfluxminus[ic], Hyfluxplus[ic], Hyfluxminus[ic])
                  - wminusx_H[ic] + wplusx_H[ic] - wminusy_H[ic] + wplusy_H[ic];
      U_new[ic] = U_fullstep(deltaT, dxic, U[ic],
                       Uxfluxplus[ic], Uxfluxminus[ic], Uyfluxplus[ic], Uyfluxminus[ic])
                  - wminusx_U[ic] + wplusx_U[ic];
      V_new[ic] = U_fullstep(deltaT, dxic, V[ic],
                       Vxfluxplus[ic], Vxfluxminus[ic], Vyfluxplus[ic], Vyfluxminus[ic])
                  - wminusy_V[ic] + wplusy_V[ic];

}

__kernel void calc_finite_difference_via_face_in_place_face_comps_cl(
            __global    const int       *nfaces,                    // 0 Number of faces
                        const int       levmx,                      // 1 Maximum level
            __global    const state_t   *H,                         // 2
            __global    const state_t   *U,                         // 3
            __global    const state_t   *V,                         // 4
            __global    const uchar_t   *level,                     // 5 Array of level information
                        const real_t    deltaT,                     // 6 Size of time step
            __global    const real_t    *lev_dx,                    // 7
            __global    const real_t    *lev_dy,                    // 8
            __global    const int       *map_xface2cell_lower,      // 9 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 10 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 11 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 12 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 13 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 14 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 15 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 16 A cell's top primary face 
            __global          state_t   *HxFlux,                    // 17
            __global          state_t   *UxFlux,                    // 18
            __global          state_t   *VxFlux,                    // 19
            __global          state_t   *HyFlux,                    // 20
            __global          state_t   *UyFlux,                    // 21
            __global          state_t   *VyFlux,                    // 22
            __global          state_t   *Wx_H,                      // 23
            __global          state_t   *Wx_U,                      // 24
            __global          state_t   *Wy_H,                      // 25
            __global          state_t   *Wy_V) {                    // 26

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    
    if (giX >= max(nfaces[0], nfaces[1])) 
        return;

    real_t   g     = 9.80; 
    real_t   ghalf = HALF*g;

    int iface = giX;

    if (giX < nfaces[0]) {
        int cell_lower = map_xface2cell_lower[iface];
        int cell_upper = map_xface2cell_upper[iface];
  
        int fl = map_xcell2face_left1[cell_lower];
        int fr = map_xcell2face_right1[cell_upper];
        if (fl != -1 && fr != -1) {
        real_t Hx, Ux, Vx;
  
        int nll = map_xface2cell_lower[fl];
        int nrr = map_xface2cell_upper[fr];
  
        uchar_t lev = level[cell_lower];
        real_t dxic    = lev_dx[lev];
        real_t Cxhalf = 0.5*deltaT/dxic;
  
        real_t Hic = H[cell_lower];
        real_t Hr  = H[cell_upper];
        real_t Hl  = H[nll];
        real_t Hrr = H[nrr];
        real_t Uic = U[cell_lower];
        real_t Ur  = U[cell_upper];
        real_t Ul  = U[nll];
        real_t Urr = U[nrr];
  
        Hx=HALF*(H[cell_upper]+H[cell_lower]) - Cxhalf*( HXFLUX(cell_upper)-HXFLUX(cell_lower) );
        Ux=HALF*(U[cell_upper]+U[cell_lower]) - Cxhalf*( UXFLUX(cell_upper)-UXFLUX(cell_lower) );
        Vx=HALF*(V[cell_upper]+V[cell_lower]) - Cxhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );

        real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);

        Wx_H[iface] = w_corrector(deltaT, dxic, U_eigen, Hr-Hic, Hic-Hl, Hrr-Hr) * (Hr - Hic);
        Wx_U[iface] = w_corrector(deltaT, dxic, U_eigen, Ur-Uic, Uic-Ul, Urr-Ur) * (Ur - Uic);

        HxFlux[iface] = HXFLUXFACE;
        UxFlux[iface] = UXFLUXFACE;
        VxFlux[iface] = VXFLUXFACE;
        }
    }

    if (giX < nfaces[1]) {
        int cell_lower = map_yface2cell_lower[iface];
        int cell_upper = map_yface2cell_upper[iface];
  
        int fb = map_ycell2face_bot1[cell_lower];
        int ft = map_ycell2face_top1[cell_upper];
        if (fb != -1 && ft != -1) {
        real_t Hy, Uy, Vy;
  
        int nbb = map_yface2cell_lower[fb];
        int ntt = map_yface2cell_upper[ft];
  
        uchar_t lev = level[cell_lower];
        real_t dyic    = lev_dy[lev];
        real_t Cyhalf = 0.5*deltaT/dyic;
  
        real_t Hic = H[cell_lower];
        real_t Ht  = H[cell_upper];
        real_t Hb  = H[nbb];
        real_t Htt = H[ntt];
        real_t Vic = V[cell_lower];
        real_t Vt  = V[cell_upper];
        real_t Vb  = V[nbb];
        real_t Vtt = V[ntt];
  
        Hy=HALF*(H[cell_upper]+H[cell_lower]) - Cyhalf*( HYFLUX(cell_upper)-HYFLUX(cell_lower) );
        Uy=HALF*(U[cell_upper]+U[cell_lower]) - Cyhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );
        Vy=HALF*(V[cell_upper]+V[cell_lower]) - Cyhalf*( VYFLUX(cell_upper)-VYFLUX(cell_lower) );

        real_t U_eigen = fabs(Vy/Hy) + sqrt(g*Hy);
  
        Wy_H[iface] = w_corrector(deltaT, dyic, U_eigen, Ht-Hic, Hic-Hb, Htt-Ht) * (Ht - Hic);
        Wy_V[iface] = w_corrector(deltaT, dyic, U_eigen, Vt-Vic, Vic-Vb, Vtt-Vt) * (Vt - Vic);
  
        HyFlux[iface] = HYFLUXFACE;
        UyFlux[iface] = UYFLUXFACE;
        VyFlux[iface] = VYFLUXFACE;
        }
    }

}

__kernel void calc_finite_difference_via_face_in_place_fixup_cl(
                        const int        nxfixup,                   // 0
                        const int        nyfixup,                   // 1
            __global    const int       *xrecvIdx,                  // 2
            __global    const int       *xsendIdx1,                 // 3
            __global    const int       *xsendIdx2,                 // 4
            __global    const int       *yrecvIdx,                  // 5
            __global    const int       *ysendIdx1,                 // 6
            __global    const int       *ysendIdx2,                 // 7
            __global    const int       *map_xface2cell_lower,      // 8
            __global    const int       *map_xface2cell_upper,      // 9
            __global    const int       *map_yface2cell_lower,      // 10 
            __global    const int       *map_yface2cell_upper,      // 11 
            __global          state_t   *HxFlux,                    // 12
            __global          state_t   *UxFlux,                    // 13
            __global          state_t   *VxFlux,                    // 14
            __global          state_t   *HyFlux,                    // 15
            __global          state_t   *UyFlux,                    // 16
            __global          state_t   *VyFlux,                    // 17
            __global          state_t   *Wx_H,                      // 18
            __global          state_t   *Wx_U,                      // 19
            __global          state_t   *Wy_H,                      // 20
            __global          state_t   *Wy_V) {                    // 21

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    
    if (giX >= max(nxfixup, nyfixup)) 
        return;

    int ifixup = nxfixup;

    if (giX < nxfixup) {
        int ir  = xrecvIdx[ifixup];
        int is1 = xsendIdx1[ifixup];
        int is2 = xsendIdx2[ifixup];
        HxFlux[ir] = (HxFlux[is1] + HxFlux[is2]) * HALF;
        UxFlux[ir] = (UxFlux[is1] + UxFlux[is2]) * HALF;
        VxFlux[ir] = (VxFlux[is1] + VxFlux[is2]) * HALF;
        Wx_H[ir] = (Wx_H[is1] + Wx_H[is2]) * 0.25;
        Wx_U[ir] = (Wx_U[is1] + Wx_U[is2]) * 0.25;
    }

    ifixup = nyfixup;

    if (giX < nyfixup) {
        int ir  = yrecvIdx[ifixup];
        int is1 = ysendIdx1[ifixup];
        int is2 = ysendIdx2[ifixup];
        HyFlux[ir] = (HyFlux[is1] + HyFlux[is2]) * HALF;
        UyFlux[ir] = (UyFlux[is1] + UyFlux[is2]) * HALF;
        VyFlux[ir] = (VyFlux[is1] + VyFlux[is2]) * HALF;
        Wy_H[ir] = (Wy_H[is1] + Wy_H[is2]) * 0.25;
        Wy_V[ir] = (Wy_V[is1] + Wy_V[is2]) * 0.25;
    }

}

__kernel void calc_finite_difference_via_face_in_place_fill_new_cl(
                        const int       ncells,                     // 0 Number of cells (not including phantom)
                        const real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global    const int       *level,                     // 4
            __global    const int       *map_xcell2face_left1,      // 5 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 6 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 7 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 8 A cell's top primary face 
            __global    const state_t   *HxFlux,                    // 9
            __global    const state_t   *UxFlux,                    // 10
            __global    const state_t   *VxFlux,                    // 11
            __global    const state_t   *HyFlux,                    // 12
            __global    const state_t   *UyFlux,                    // 13
            __global    const state_t   *VyFlux,                    // 14
            __global    const state_t   *Wx_H,                      // 15
            __global    const state_t   *Wx_U,                      // 16
            __global    const state_t   *Wy_H,                      // 17
            __global    const state_t   *Wy_V,                      // 18
            __global    const state_t   *H,                         // 19
            __global    const state_t   *U,                         // 20
            __global    const state_t   *V,                         // 21
            __global          state_t   *H_new,                     // 22
            __global          state_t   *U_new,                     // 23
            __global          state_t   *V_new) {                   // 24

    /////////////////////////////////////////////
    /// Get thread identification information ///
    /////////////////////////////////////////////

    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);

    const uint ngX = get_global_size(0);
    const uint ntX = get_local_size(0);

    const uint group_id = get_group_id(0);

    
    if (giX >= ncells) 
        return;

      int ic = giX;

      real_t dxic = lev_dx[level[ic]];
      int fl = map_xcell2face_left1[ic];
      int fr = map_xcell2face_right1[ic];
      int fb = map_ycell2face_bot1[ic];
      int ft = map_ycell2face_top1[ic];

      //if (fl == fr || fb == ft) return;

      H_new[ic] = U_fullstep(deltaT,dxic,H[ic],
                  HxFlux[fr], HxFlux[fl], HyFlux[ft], HyFlux[fb])
                - Wx_H[fl] + Wx_H[fr] - Wy_H[fb] + Wy_H[ft];
      U_new[ic] = U_fullstep(deltaT,dxic,U[ic],
                  UxFlux[fr], UxFlux[fl], UyFlux[ft], UyFlux[fb])
                - Wx_U[fl] + Wx_U[fr];
      V_new[ic] = U_fullstep(deltaT,dxic,V[ic],
                  VxFlux[fr], VxFlux[fl], VyFlux[ft], VyFlux[fb])
                - Wy_V[fb] + Wy_V[ft];

}

__kernel void refine_potential_cl(
                 const int      ncells,     // 0  Total number of cells.
                 const int      levmx,      // 1  Maximum level
        __global const state_t *H,          // 2
        __global const state_t *U,          // 3
        __global const state_t *V,          // 4
        __global const int     *nlft,       // 5  Array of left neighbors.
        __global const int     *nrht,       // 6  Array of right neighbors.
        __global const int     *ntop,       // 7  Array of bottom neighbors.
        __global const int     *nbot,       // 8  Array of top neighbors.
        __global const int     *i,          // 9  Array of i values.
        __global const int     *j,          // 10 Array of j values.
        __global const uchar_t *level,      // 11 Array of level information.
        __global const char_t  *celltype,   // 12  Array of celltype information.
        __global const real_t  *lev_dx,     // 13
        __global const real_t  *lev_dy,     // 14
        __global       char_t  *mpot,       // 15  Array of mesh potential information.
        __global       int2    *redscratch, // 16  Array of new giX offsets.
        __global       int2    *result,     // 17
        __local        state_t *tile,       // 18  Tile size in state_t.
        __local        int8    *itile)      // 19  Tile size in int8.
{

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
////////////////////                       ///////////////////////
////////////////////   calc_gradients_cl   ///////////////////////
////////////////////                       ///////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

   state_t qpot = ZERO;
    
   /////////////////////////////////////////////
   /// Get thread identification information ///
   /////////////////////////////////////////////

   const uint giX  = get_global_id(0);
   const uint tiX  = get_local_id(0);

   const uint group_id = get_group_id(0);

   const uint ntX  = get_local_size(0);

   itile[tiX].s01 = 0;

   if(giX >= ncells)
      return;

   setup_refine_tile(tile, itile, ncells, H, nlft, nrht, ntop, nbot, level);

   char_t ctype = celltype[giX];

   barrier (CLK_LOCAL_MEM_FENCE);

   int nlt, nrt, nbr, ntr;
   int nl, nr, nb, nt;
   state_t Hic;
   real_t Hl; //, Ul, Vl;
   state_t Hr; //, Ur, Vr;
   state_t Hb; //, Ub, Vb;
   state_t Ht; //, Ut, Vt;

   state_t duminus1;
   state_t duplus1;
   state_t duhalf1;

   nl = nlftval(tiX);
   nr = nrhtval(tiX);
   nb = nbotval(tiX);
   nt = ntopval(tiX);

   Hic  = Hrefval(tiX);

   uchar_t lvl = levelval(tiX);

   //////////////////////////
   //////////////////////////
   //////////////////////////

   //barrier(CLK_LOCAL_MEM_FENCE);

   //////////////////////////////
   //////////////////////////////
   //////////////////////////////

   int mpotval = 0;

   if (ctype == REAL_CELL){
      // Setting the left and left-left neighbor state values for the control volume based
      // on the state variables of the actual left, left-left, and left-top neighbor cells

      // Using global access for the left neighbor values
      if(nl < 0) {
         nl = abs(nl+1);
         Hl = H[nl];
         //Ul = U[nl];
         //Vl = V[nl];

         if (level[nl] > lvl) {
            nlt = ntop[nl];
            Hl = REFINE_HALF * (Hl + H[nlt]);
            //Ul = REFINE_HALF * (Ul + U[nlt]);
            //Vl = REFINE_HALF * (Vl + V[nlt]);
         }
      }
      // Using local access for the left neighbor
      else {
         Hl = Hrefval(nl);
         //Ul = Uval(nl);
         //Vl = Vval(nl);

         // The left neighbor is more refined than the current cell
         if (levelval(nl) > lvl) {
            nlt = ntopval(nl);
            if(nlt >= 0) {
               Hl = REFINE_HALF * (Hl + Hrefval(nlt));
               //Ul = REFINE_HALF * (Ul + Uval(nlt));
               //Vl = REFINE_HALF * (Vl + Vval(nlt));
            }
            else {
               nlt = abs(nlt+1);
               Hl = REFINE_HALF * (Hl + H[nlt]);
               //Ul = REFINE_HALF * (Ul + U[nlt]);
               //Vl = REFINE_HALF * (Vl + V[nlt]);
            }
         }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////

      // Setting the right and right-right neighbor state values for the control volume based
      // on the state variables of the actual right, right-right, and right-top neighbor cells

      // Using global access for the right neighbor values
      if(nr < 0) {
         nr = abs(nr+1);
         Hr = H[nr];
         //Ur = U[nr];
         //Vr = V[nr];
   
         if (level[nr] > lvl) {
            nrt = ntop[nr];
            Hr = REFINE_HALF * (Hr + H[nrt]);
            //Ur = REFINE_HALF * (Ur + U[nrt]);
            //Vr = REFINE_HALF * (Vr + V[nrt]);
         }
      }
      // Using local access for the right neighbor
      else {
         Hr = Hrefval(nr);
         //Ur = Uval(nr);
         //Vr = Vval(nr);

         if (levelval(nr) > lvl) {
            nrt = ntopval(nr);
            if(nrt >= 0) {
               Hr = REFINE_HALF * (Hr + Hrefval(nrt));
               //Ur = REFINE_HALF * (Ur + Uval(nrt));
               //Vr = REFINE_HALF * (Vr + Vval(nrt));
            }
            else {
               nrt = abs(nrt+1);
               Hr = REFINE_HALF * (Hr + H[nrt]);
               //Ur = REFINE_HALF * (Ur + U[nrt]);
               //Vr = REFINE_HALF * (Vr + V[nrt]);
            }
         }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////



      // Setting the bottom and bottom-bottom neighbor state values for the control volume based
      // on the state variables of the actual bottom, bottom-bottom, and bottom-right neighbor cells

      // Using global access for the bottom neighbor values
      if (nb < 0) {
         nb = abs(nb+1);
         Hb = H[nb];
         //Ub = U[nb];
         //Vb = V[nb];

         if (level[nb] > lvl) {
            nbr = nrht[nb];
            Hb = REFINE_HALF * (Hb + H[nbr]);
            //Ub = REFINE_HALF * (Ub + U[nbr]);
            //Vb = REFINE_HALF * (Vb + V[nbr]);
         }
      }
      // Using local access for the bottom neighbor
      else {
         Hb = Hrefval(nb);
         //Ub = Uval(nb);
         //Vb = Vval(nb);

         if (levelval(nb) > lvl) {
            nbr = nrhtval(nb);
            if(nbr >= 0) {
               Hb = REFINE_HALF * (Hb + Hrefval(nbr));
               //Ub = REFINE_HALF * (Ub + Uval(nbr));
               //Vb = REFINE_HALF * (Vb + Vval(nbr));
            }
            else {
               nbr = abs(nbr+1);
               Hb = REFINE_HALF * (Hb + H[nbr]);
               //Ub = REFINE_HALF * (Ub + U[nbr]);
               //Vb = REFINE_HALF * (Vb + V[nbr]);
            }
         }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////


      // Setting the top and top-top neighbor state values for the control volume based
      // on the state variables of the actual top, top-top, and top-right neighbor cells
  
      // Using global access for the top neighbor values
      if (nt < 0) {
         nt = abs(nt+1);
         Ht = H[nt];
         //Ut = U[nt];
         //Vt = V[nt];

         if (level[nt] > lvl) {
            ntr = nrht[nt];
            Ht = REFINE_HALF * (Ht + H[ntr]);
            //Ut = REFINE_HALF * (Ut + U[ntr]);
            //Vt = REFINE_HALF * (Vt + V[ntr]);
         }
      }
      // Using local access for the top neighbor
      else {
         Ht = Hrefval(nt);
         //Ut = Uval(nt);
         //Vt = Vval(nt);

         if (levelval(nt) > lvl) {
            ntr = nrhtval(nt);
            if(ntr >= 0) {
               Ht = REFINE_HALF * (Ht + Hrefval(ntr));
               //Ut = REFINE_HALF * (Ut + Uval(ntr));
               //Vt = REFINE_HALF * (Vt + Vval(ntr));
            }
            else {
               ntr = abs(ntr+1);
               Ht = REFINE_HALF * (Ht + H[ntr]);
               //Ut = REFINE_HALF * (Ut + U[ntr]);
               //Vt = REFINE_HALF * (Vt + V[ntr]);
            }
         }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////

       //--CALCULATIONS------------------------------------------------------------
       //  Calculate the gradient between the right and left neighbors and the
       //  main cell.
       state_t invHic = REFINE_ONE / Hic; //  For faster math.
       state_t qmax = REFINE_NEG_THOUSAND;          //  Set the default maximum low to catch the real one.
    
       duplus1 = Hr - Hic;
       duhalf1 = Hic - Hl;
       qpot = max(fabs(duplus1 * invHic), fabs(duhalf1 * invHic));
       if (qpot > qmax) qmax = qpot;
    
       duminus1= Hic - Hl;
       duhalf1 = Hr - Hic;
       qpot = max(fabs(duminus1 * invHic), fabs(duhalf1 * invHic));
       if (qpot > qmax) qmax = qpot;
    
       //  Calculate the gradient between the top and bottom neighbors and the
       //  main cell.
       duplus1 = Ht - Hic;
       duhalf1 = Hic - Hb;
       qpot = max(fabs(duplus1 * invHic), fabs(duhalf1 * invHic));
       if (qpot > qmax) qmax = qpot;
    
       duminus1= Hic - Hb;
       duhalf1 = Ht - Hic;
       qpot = max(fabs(duminus1 * invHic), fabs(duhalf1 * invHic));
       if (qpot > qmax) qmax = qpot;


       //--CALCULATIONS------------------------------------------------------------
       //  Refine the mesh if the gradient is large
       if (qmax > REFINE_GRADIENT && lvl < levmx) {
          mpotval = 1;
       } else if (qmax < COARSEN_GRADIENT && lvl > 0) {
          mpotval = -1;
       }
    } 

    barrier(CLK_LOCAL_MEM_FENCE);

    // Reusing itile, so we have to reset it
    itile[tiX].s01 = 0;

    if (mpotval > 0){
       itile[tiX].s0 = (ctype == REAL_CELL) ? 3 : 1;
    } else if (mpotval < 0) {
       int ival = i[giX];
       int jval = j[giX];
       if (ctype == REAL_CELL) {
          if (! is_lower_left(ival,jval) ) itile[tiX].s1 = 1;
       } else {
          if (! is_upper_right(ival,jval) || is_lower_left(ival,jval) ) itile[tiX].s1 = 1;
       }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    reduction_sum_int2_within_tile(itile);

    if (tiX == 0) {
      redscratch[group_id].s01 = itile[0].s01;
      result[0].s01 = itile[0].s01;
    }

    /////////////////////
    /// GLOBAL WRITES ///
    /////////////////////

    //  Put the mesh potential on the global array.
    mpot[giX] = mpotval; /**/
}

void setup_tile(__local        state4_t *tile, 
                __local        int8     *itile, 
                         const int      isize,
                __global const state_t  *H,
                __global const state_t  *U,
                __global const state_t  *V,
                __global const int      *nlft,
                __global const int      *nrht,
                __global const int      *ntop,
                __global const int      *nbot,
                __global const uchar_t  *level
                )
{
   const uint giX = get_global_id (0);
   const uint tiX = get_local_id (0);

   const uint ntX = get_local_size (0);

   const uint group_id = get_group_id (0);

   int start_idx = group_id * ntX;
   int end_idx = (group_id + 1) * ntX;
   end_idx = min(end_idx, isize);

   Hval(tiX) = H[giX];
   Uval(tiX) = U[giX];
   Vval(tiX) = V[giX];

   if (nlft[giX] >= start_idx && nlft[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      nlftval(tiX) =  nlft[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      nlftval(tiX) = -nlft[giX]-1;
   }

   if (nrht[giX] >= start_idx && nrht[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      nrhtval(tiX) =  nrht[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      nrhtval(tiX) = -nrht[giX]-1;
   }

   if (ntop[giX] >= start_idx && ntop[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      ntopval(tiX) =  ntop[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      ntopval(tiX) = -ntop[giX]-1;
   }

   if (nbot[giX] >= start_idx && nbot[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      nbotval(tiX) =  nbot[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      nbotval(tiX) = -nbot[giX]-1;
   }

   levelval(tiX) = level[giX];
}

void setup_refine_tile(
                __local        state_t *tile, 
                __local        int8    *itile, 
                         const int     isize,
                __global const state_t *H,
                __global const int     *nlft,
                __global const int     *nrht,
                __global const int     *ntop,
                __global const int     *nbot,
                __global const uchar_t *level
                )
{
   const unsigned int giX = get_global_id (0);
   const unsigned int tiX = get_local_id (0);

   const unsigned int ntX = get_local_size (0);

   const unsigned int group_id = get_group_id (0);

   int start_idx = group_id * ntX;
   int end_idx = (group_id + 1) * ntX;
   end_idx = min(end_idx, isize);

   Hrefval(tiX) = H[giX];

   if (nlft[giX] >= start_idx && nlft[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      nlftval(tiX) =  nlft[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      nlftval(tiX) = -nlft[giX]-1;
   }

   if (nrht[giX] >= start_idx && nrht[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      nrhtval(tiX) =  nrht[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      nrhtval(tiX) = -nrht[giX]-1;
   }

   if (ntop[giX] >= start_idx && ntop[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      ntopval(tiX) =  ntop[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      ntopval(tiX) = -ntop[giX]-1;
   }

   if (nbot[giX] >= start_idx && nbot[giX] < end_idx) {
      // If on block, offset to local index by subtracting start index
      nbotval(tiX) =  nbot[giX] - start_idx;
   } else {
      // If off block, set to negative to indicate global data
      nbotval(tiX) = -nbot[giX]-1;
   }

   levelval(tiX) = level[giX];
}

/*void setup_xface(
                __local           int8        *xface,
                __local           int         isize;
                __global    const int         *map_xface2cell_lower,   
                __global    const int         *map_xface2cell_upper,      
                )
{
    const unsigned int giX = get_global_id (0);
    const unsigned int tiX = get_local_id (0);

    const unsigned int ntX = get_local_size (0);

    const unsigned int group_id = get_group_id (0);

    int start_idx = group_id * ntX;
    int end_idx = (group_id + 1) * ntX;
    end_idx = min(end_idx, isize);

    if (map_xface2cell_lower[giX] >= start_idx && map_xface2cell_lower[giX] < end_idx) {
        xfacelower(tiX) = map_xface2cell_lower[giX];
    }
    else {
        xfacelower(tiX) = map_xface2cell_lower[giX];
    }

    if (map_xface2cell_upper[giX] >= start_idx && map_xface2cell_upper[giX] < end_idx) {
        xfaceupper(tiX) = map_xface2cell_upper[giX];
    }
    else {
        xfaceupper(tiX) = map_xface2cell_upper[giX];
    }
}

void setup_yface(
                __local           int8        *yface,
                __local           int         isize;
                __global    const int         *map_yface2cell_lower,   
                __global    const int         *map_yface2cell_upper,      
                )
{
    const unsigned int giX = get_global_id (0);
    const unsigned int tiX = get_local_id (0);

    const unsigned int ntX = get_local_size (0);

    const unsigned int group_id = get_group_id (0);

    int start_idx = group_id * ntX;
    int end_idx = (group_id + 1) * ntX;
    end_idx = min(end_idx, isize);

    if (map_yface2cell_lower[giX] >= start_idx && map_yface2cell_lower[giX] < end_idx) {
        yfacelower(tiX) = map_yface2cell_lower[giX];
    }
    else {
        yfacelower(tiX) = map_yface2cell_lower[giX];
    }

    if (map_xface2cell_upper[giX] >= start_idx && map_xface2cell_upper[giX] < end_idx) {
        yfaceupper(tiX) = map_yface2cell_upper[giX];
    }
    else {
        yfaceupper(tiX) = map_yface2cell_upper[giX];
    }
}*/

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

__kernel void reduce_sum_mass_stage1of2_cl(
                 const int      isize,     // 0  Total number of cells.
        __global const state_t *array,     // 1
        __global const uchar_t *level,     // 2
        __global const real_t  *levdx,     // 3
        __global const real_t  *levdy,     // 4
        __global const char_t   *celltype,  // 5
        __global       real_t  *mass_sum,  // 6
        __global       real_t  *scratch,   // 7
        __local        real_t  *tile)      // 8
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

void reduction_epsum_within_tile(__local  real2_t  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   real_t corrected_next_term, new_sum;

    for (int offset = ntX >> 1; offset > MIN_REDUCE_SYNC_SIZE; offset >>= 1)
    {
        if (tiX < offset)
        {
            corrected_next_term = tile[tiX+offset].s0 + (tile[tiX+offset].s1 +tile[tiX].s1);
            new_sum = tile[tiX].s0 + corrected_next_term;
            tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
            tile[tiX].s0 = new_sum;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (tiX < MIN_REDUCE_SYNC_SIZE)
    {
        for (int offset = MIN_REDUCE_SYNC_SIZE; offset > 1; offset >>= 1)
        {
            corrected_next_term = tile[tiX+offset].s0 + (tile[tiX+offset].s1 +tile[tiX].s1);
            new_sum = tile[tiX].s0 + corrected_next_term;
            tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
            tile[tiX].s0 = new_sum;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        corrected_next_term = tile[tiX+1].s0 + (tile[tiX+1].s1 +tile[tiX].s1);
        new_sum = tile[tiX].s0 + corrected_next_term;
        tile[tiX].s1 = corrected_next_term - (new_sum - tile[tiX].s0);
        tile[tiX].s0 = new_sum;
    }
}

__kernel void reduce_epsum_mass_stage1of2_cl(
                 const int      isize,     // 0  Total number of cells.
        __global const state_t *array,     // 1
        __global const uchar_t *level,     // 2
        __global const real_t  *levdx,     // 3
        __global const real_t  *levdy,     // 4
        __global const char_t   *celltype,  // 5
        __global       real2_t *mass_sum, // 6
        __global       real2_t *scratch,   // 7
        __local        real2_t *tile)      // 8
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
        __global       real_t *mass_sum,
        __global       real_t *scratch,
        __local        real_t *tile)
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
                 const int      isize,
        __global       real2_t  *mass_sum,
        __global       real2_t  *scratch,
        __local        real2_t  *tile)
{
   const unsigned int tiX  = get_local_id(0);
   const unsigned int ntX  = get_local_size(0);
   real_t corrected_next_term, new_sum;

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

