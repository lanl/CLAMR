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

void reduction_min_within_tile1(__local  real_t  *tile)
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

__kernel void set_timestep_cl(
                 const int     ncells,     // 0  Total number of cells.
                 const real_t  sigma,      // 1
        __global const state_t *H_in,      // 2  
        __global const state_t *U_in,      // 3  
        __global const state_t *V_in,      // 4  
        __global const int     *level,     // 5  Array of level information.
        __global const int     *celltype,  // 6 
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
    int lev      = level[giX];
    int type     = celltype[giX];
    
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
                __global const int      *level
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
                __global const int     *level
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

#define SQ(x)      ( (x)*(x) )
#define MIN3(a,b,c) ( min(min((a),(b)),(c)) )

#define HXFLUX(ic)  ( U_old[ic] )
#define UXFLUX(ic)  ( SQ(U_old[ic])/H_old[ic] + ghalf*SQ(H_old[ic]) )
#define UVFLUX(ic)  ( U_old[ic]*V_old[ic]/H_old[ic] )

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

#define HYFLUX(ic)  ( V_old[ic] )
#define VUFLUX(ic)  ( V_old[ic]*U_old[ic]/H_old[ic] )
#define VYFLUX(ic)  ( SQ(V_old[ic])/H_old[ic] + ghalf*SQ(H_old[ic]) )

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

// XXX Added XXX
#define HXFLUXNLT ( Ult )
#define HXFLUXNRT ( Urt )

// XXX Added XXX
#define UXFLUXNLT ( SQ(Ult)/Hlt + ghalf*SQ(Hlt) )
#define UXFLUXNRT ( SQ(Urt)/Hrt + ghalf*SQ(Hrt) )

// XXX Added XXX
#define UVFLUXNLT ( Ult*Vlt/Hlt )
#define UVFLUXNRT ( Urt*Vrt/Hrt )

// XXX Added XXX
#define HYFLUXNBR ( Vbr )
#define HYFLUXNTR ( Vtr )

// XXX Added XXX
#define VUFLUXNBR  ( Vbr*Ubr/Hbr )
#define VUFLUXNTR  ( Vtr*Utr/Htr )

// XXX Added XXX
#define VYFLUXNBR  ( SQ(Vbr)/Hbr + ghalf*SQ(Hbr) )
#define VYFLUXNTR  ( SQ(Vtr)/Htr + ghalf*SQ(Htr) )

// XXX Added XXX
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
        __global const int     *celltype, // 1  Array of left neighbors
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

   int ctype = celltype[giX];

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
        __global const int    *celltype, // 1  Array of left neighbors
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

   int ctype = celltype[giX];

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
        __global const int     *celltype, // 1  Array of left neighbors
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

   int ctype = celltype[giX];

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
        __global const int      *level,    // 12  Array of level information
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
      if(lvl_nt < level[ntrt]) // Note that lvl_nt is the same as lvl_ntr)
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

#ifdef ATI
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
#ifdef IS_NVIDIA
    H_new[giX] = Hic;
    U_new[giX] = Uic;
    V_new[giX] = Vic;
#endif

//////////////////////////////////////////////////////////////////
////////////////////          END           //////////////////////
////////////////////   calc_one_cycle_cl    //////////////////////
//////////////////////////////////////////////////////////////////

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
        __global const int     *level,      // 11 Array of level information.
        __global const int     *celltype,   // 12  Array of celltype information.
        __global const real_t  *lev_dx,     // 13
        __global const real_t  *lev_dy,     // 14
        __global       int     *mpot,       // 15  Array of mesh potential information.
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

   int ctype = celltype[giX];

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

   int lvl = levelval(tiX);

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

    for (int offset = ntX >> 1; offset > 32; offset >>= 1) {
       if (tiX < offset) {
          itile[tiX].s01 += itile[tiX+offset].s01; 
       }
       barrier(CLK_LOCAL_MEM_FENCE);
    }

    //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
    if (tiX < 32)
    {  itile[tiX].s01 += itile[tiX+32].s01;
       barrier(CLK_LOCAL_MEM_FENCE);
       itile[tiX].s01 += itile[tiX+16].s01;
       barrier(CLK_LOCAL_MEM_FENCE);
       itile[tiX].s01 += itile[tiX+8].s01;
       barrier(CLK_LOCAL_MEM_FENCE);
       itile[tiX].s01 += itile[tiX+4].s01;
       barrier(CLK_LOCAL_MEM_FENCE);
       itile[tiX].s01 += itile[tiX+2].s01;
       barrier(CLK_LOCAL_MEM_FENCE);
       itile[tiX].s01 += itile[tiX+1].s01; }

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
                __global const int      *level
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
                __global const int     *level
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

void reduction_sum_within_tile(__local  real_t  *tile)
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


__kernel void reduce_sum_mass_stage1of2_cl(
                 const int      isize,     // 0  Total number of cells.
        __global const state_t *array,     // 1
        __global const int     *level,     // 2
        __global const real_t  *levdx,     // 3
        __global const real_t  *levdy,     // 4
        __global const int     *celltype,  // 5
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

__kernel void reduce_epsum_mass_stage1of2_cl(
                 const int      isize,     // 0  Total number of cells.
        __global const state_t *array,     // 1
        __global const int     *level,     // 2
        __global const real_t  *levdx,     // 3
        __global const real_t  *levdy,     // 4
        __global const int     *celltype,  // 5
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

