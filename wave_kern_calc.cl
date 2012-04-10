/*
 *  Copyright (c) 2011, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2011. Los Alamos National Security, LLC. This software was produced 
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
typedef double4 real4;
typedef double8 real8;
#define ZERO 0.0
#define HALF 0.5
#define QUARTER 0.25
#define ONE  1.0
#define GRAVITATIONAL_CONSTANT 9.80
#define THOUSAND 1000.0
#define EPSILON 1.0e-30
#else
typedef float   real;
typedef float4  real4;
typedef float8  real8;
#define ZERO 0.0f
#define HALF 0.5f
#define QUARTER 0.25f
#define ONE  1.0f
#define GRAVITATIONAL_CONSTANT 9.80f
#define THOUSAND 1000.0f
#define EPSILON 1.0f-30
#endif
#endif

//#ifdef __APPLE_CC__
//#define max(a,b) ((a) > (b) ? (a) : (b))
//#define fabs(a) ( (a) < 0 ? -(a) : a)
//#endif
void setup_tile(__local        real4  *tile,
                __local        int8   *itile,
                         const int    isize,
                __global const real   *H,
                __global const real   *U,
                __global const real   *V,
                __global const int    *nlft,
                __global const int    *nrht,
                __global const int    *ntop,
                __global const int    *nbot,
                __global const int    *level
                );

void apply_BCs(__local  real4        *tile,
               __local  int8         *itile,
               __global const real   *H,
               __global const real   *U,
               __global const real   *V,
               __global const int    *nlft,
               __global const int    *nrht,
               __global const int    *ntop,
               __global const int    *nbot);

__kernel void copy_state_data_cl(
                          const int  isize,         // 0 
                 __global      real  *H,            // 1 
                 __global      real  *U,            // 2 
                 __global      real  *V,            // 3 
                 __global      real  *H_new,        // 4 
                 __global      real  *U_new,        // 5 
                 __global      real  *V_new)        // 6 
{
                
   const unsigned int giX  = get_global_id(0);

   if (giX >= isize) return;

   H_new[giX]    = H[giX];
   U_new[giX]    = U[giX];
   V_new[giX]    = V[giX];
}

__kernel void copy_state_ghost_data_cl(
                          const int  ncells,        // 0 
                          const int  nghost,        // 1 
                 __global      real  *H,            // 2 
                 __global      real  *H_add,        // 3 
                 __global      real  *U,            // 4 
                 __global      real  *U_add,        // 5 
                 __global      real  *V,            // 6 
                 __global      real  *V_add)        // 7 
{
   const unsigned int giX  = get_global_id(0);

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



#define Hptrval(i)     ( (*tile)[i].s0 )

#define nlftval(i)  ( itile[i].s0 )
#define nrhtval(i)  ( itile[i].s1 )
#define ntopval(i)  ( itile[i].s2 )
#define nbotval(i)  ( itile[i].s3 )
#define levelval(i) ( itile[i].s4 )
#define mpotval(i)  ( itile[i].s5 )
#endif

#define SQR(x)      ( (x)*(x) )
#define MIN3(a,b,c) ( min(min((a),(b)),(c)) )

#define HXFLUX(ic)  ( U_old[ic] )
#define UXFLUX(ic)  ( SQR(U_old[ic])/H_old[ic] + ghalf*SQR(H_old[ic]) )
#define UVFLUX(ic)  ( U_old[ic]*V_old[ic]/H_old[ic] )

#define HXFLUXIC ( Uic )
#define HXFLUXNL ( Ul )
#define HXFLUXNR ( Ur )
#define HXFLUXNB ( Ub )
#define HXFLUXNT ( Ut )

#define UXFLUXIC ( SQR(Uic)/Hic + ghalf*SQR(Hic) )
#define UXFLUXNL ( SQR(Ul)/Hl + ghalf*SQR(Hl) )
#define UXFLUXNR ( SQR(Ur)/Hr + ghalf*SQR(Hr) )
#define UXFLUXNB ( SQR(Ub)/Hb + ghalf*SQR(Hb) )
#define UXFLUXNT ( SQR(Ut)/Ht + ghalf*SQR(Ht) )

#define UVFLUXIC ( Uic*Vic/Hic )
#define UVFLUXNL ( Ul*Vl/Hl )
#define UVFLUXNR ( Ur*Vr/Hr )
#define UVFLUXNB ( Ub*Vb/Hb )
#define UVFLUXNT ( Ut*Vt/Ht )

#define HYFLUX(ic)  ( V_old[ic] )
#define VUFLUX(ic)  ( V_old[ic]*U_old[ic]/H_old[ic] )
#define VYFLUX(ic)  ( SQR(V_old[ic])/H_old[ic] + ghalf*SQR(H_old[ic]) )

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

#define VYFLUXIC  ( SQR(Vic)/Hic + ghalf*SQR(Hic) )
#define VYFLUXNL  ( SQR(Vl)/Hl + ghalf*SQR(Hl) )
#define VYFLUXNR  ( SQR(Vr)/Hr + ghalf*SQR(Hr) )
#define VYFLUXNB  ( SQR(Vb)/Hb + ghalf*SQR(Hb) )
#define VYFLUXNT  ( SQR(Vt)/Ht + ghalf*SQR(Ht) )

#define HNEWXFLUXMINUS  ( Uxminus )
#define HNEWXFLUXPLUS   ( Uxplus )
#define UNEWXFLUXMINUS  ( SQR(Uxminus)/Hxminus + ghalf*SQR(Hxminus) )
#define UNEWXFLUXPLUS   ( SQR(Uxplus) /Hxplus +  ghalf*SQR(Hxplus)  )
#define UVNEWFLUXMINUS  ( Uxminus*Vxminus/Hxminus )
#define UVNEWFLUXPLUS   ( Uxplus *Vxplus /Hxplus  )

#define HNEWYFLUXMINUS  ( Vyminus )
#define HNEWYFLUXPLUS   ( Vyplus  )
#define VNEWYFLUXMINUS  ( SQR(Vyminus)/Hyminus + ghalf*SQR(Hyminus) )
#define VNEWYFLUXPLUS   ( SQR(Vyplus) /Hyplus  + ghalf*SQR(Hyplus)  )
#define VUNEWFLUXMINUS  ( Vyminus*Uyminus/Hyminus )
#define VUNEWFLUXPLUS   ( Vyplus *Uyplus /Hyplus )

// XXX Added XXX
#define HXFLUXNLT ( Ult )
#define HXFLUXNRT ( Urt )

// XXX Added XXX
#define UXFLUXNLT ( SQR(Ult)/Hlt + ghalf*SQR(Hlt) )
#define UXFLUXNRT ( SQR(Urt)/Hrt + ghalf*SQR(Hrt) )

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
#define VYFLUXNBR  ( SQR(Vbr)/Hbr + ghalf*SQR(Hbr) )
#define VYFLUXNTR  ( SQR(Vtr)/Htr + ghalf*SQR(Htr) )

// XXX Added XXX
#define HNEWXFLUXMINUS2  ( Uxminus2 )
#define HNEWXFLUXPLUS2   ( Uxplus2 )
#define UNEWXFLUXMINUS2  ( SQR(Uxminus2)/Hxminus2 + ghalf*SQR(Hxminus2) )
#define UNEWXFLUXPLUS2   ( SQR(Uxplus2) /Hxplus2 +  ghalf*SQR(Hxplus2)  )
#define UVNEWFLUXMINUS2  ( Uxminus2*Vxminus2/Hxminus2 )
#define UVNEWFLUXPLUS2   ( Uxplus2 *Vxplus2 /Hxplus2  )

#define HNEWYFLUXMINUS2  ( Vyminus2 )
#define HNEWYFLUXPLUS2   ( Vyplus2  )
#define VNEWYFLUXMINUS2  ( SQR(Vyminus2)/Hyminus2 + ghalf*SQR(Hyminus2) )
#define VNEWYFLUXPLUS2   ( SQR(Vyplus2) /Hyplus2  + ghalf*SQR(Hyplus2)  )
#define VUNEWFLUXMINUS2  ( Vyminus2*Uyminus2/Hyminus2 )
#define VUNEWFLUXPLUS2   ( Vyplus2 *Uyplus2 /Hyplus2 )

#define U_halfstep(deltaT, U_i, U_n, F_i, F_n, r_i, r_n, A_i, A_n, V_i, V_n) (( (( r_i*U_n + r_n*U_i ) / ( r_i + r_n )) - HALF*deltaT*(( F_n*A_n*min(ONE, A_i/A_n) - F_i*A_i*min(ONE, A_n/A_i) ) / ( V_n*min(HALF, V_i/V_n) + V_i*min(HALF, V_n/V_i) )) ))

real U_halfstep_ORIG(// XXX Fix the subindices to be more intuitive XXX
        real    deltaT,     // Timestep
        real    U_i,        // Initial cell's (downwind's) state variable
        real    U_n,        // Next cell's    (upwind's)   state variable
        real    F_i,        // Initial cell's (downwind's) state variable flux
        real    F_n,        // Next cell's    (upwind's)   state variable flux
        real    r_i,        // Initial cell's (downwind's) center to face dist
        real    r_n,        // Next cell's    (upwind's)   center to face dist
        real    A_i,        // Cell's            face surface area
        real    A_n,        // Cell's neighbor's face surface area
        real    V_i,        // Cell's            volume
        real    V_n) {      // Cell's neighbor's volume

   return (( r_i*U_n + r_n*U_i ) / ( r_i + r_n ))
          - HALF*deltaT*(( F_n*A_n*min(ONE, A_i/A_n) - F_i*A_i*min(ONE, A_n/A_i) )
                    / ( V_n*min(HALF, V_i/V_n) + V_i*min(HALF, V_n/V_i) ));

}

real U_halfstep_BD(// XXX Fix the subindices to be more intuitive XXX
        real    deltaT,     // Timestep
        real    U_i,        // Initial cell's (downwind's) state variable
        real    U_n,        // Next cell's    (upwind's)   state variable
        real    F_i,        // Initial cell's (downwind's) state variable flux
        real    F_n,        // Next cell's    (upwind's)   state variable flux
        real    r_i,        // Initial cell's (downwind's) center to face distance
        real    r_n,        // Next cell's    (upwind's)   center to face distance
        real    A_i,        // Cell's            face surface area
        real    A_n,        // Cell's neighbor's face surface area
        real    V_i,        // Cell's            volume
        real    V_n) {      // Cell's neighbor's volume

   return ((U_i*r_n+U_n*r_i)/(r_n+r_i) + (deltaT/(r_n+r_i))*(F_n-F_i));

}

#define U_fullstep(deltaT, dr, U, F_plus, F_minus, G_plus, G_minus) (( (U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus)) ))

real U_fullstep_ORIG(
        real    deltaT,
        real    dr,
        real    U,
        real    F_plus,
        real    F_minus,
        real    G_plus,
        real    G_minus) {

   return (U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus));

}

//#define w_corrector(deltaT, dr, U_eigen, grad_half, grad_minus, grad_plus) (( HALF*(HALF*U_eigen*deltaT/dr)*(ONE-(HALF*U_eigen*deltaT/dr))*(ONE- max(MIN3(ONE, (grad_plus*grad_half/max(SQR(grad_half),EPSILON)), (grad_minus*grad_half/max(SQR(grad_half),EPSILON))), ZERO)) ))

real w_corrector(//_ORIG(
        real    deltaT,       // Timestep
        real    dr,           // Cell's center to face distance
        real    U_eigen,      // State variable's eigenvalue (speed)
        real    grad_half,    // Centered gradient
        real    grad_minus,   // Downwind gradient
        real    grad_plus) {  // Upwind gradient

   real nu     = HALF * U_eigen * deltaT / dr;
   nu          = nu * (ONE - nu);

   real rdenom = ONE / max(SQR(grad_half), EPSILON);
   real rplus  = (grad_plus  * grad_half) * rdenom;
   real rminus = (grad_minus * grad_half) * rdenom;

   return HALF*nu*(ONE- max(MIN3(ONE, rplus, rminus), ZERO));

}



__kernel void calc_finite_difference_cl(
                 const int    ncells,   // 0  Total number of cells
                 const int    levmx,    // 1  Maximum level
        __global const real  *H,//_old,    // 2  
        __global const real  *U,//_old,    // 3  
        __global const real  *V,//_old,    // 4  
        __global       real  *H_new,    // 5  
        __global       real  *U_new,    // 6  
        __global       real  *V_new,    // 7  
        __global const int   *nlft,     // 8   Array of left neighbors
        __global const int   *nrht,     // 9   Array of right neighbors
        __global const int   *ntop,     // 10  Array of top neighbors
        __global const int   *nbot,     // 11  Array of bottom neighbors
        __global const int   *level,    // 12  Array of level information
                 const real   deltaT,   // 13  Size of time step.
        __global const real  *lev_dx,   // 14
        __global const real  *lev_dy,   // 15
        __local        real4 *tile,     // 16  Tile size in real4
        __local        int8  *itile){   // 17  Tile size in int8

   /////////////////////////////////////////////
   /// Get thread identification information ///
   /////////////////////////////////////////////

   const unsigned int giX  = get_global_id(0);
   const unsigned int tiX  = get_local_id(0);
   
   const unsigned int ngX  = get_global_size(0);
   const unsigned int ntX  = get_local_size(0);
   
   const unsigned int group_id = get_group_id(0);
    
   // Ensure the executing thread is not extraneous
   if(giX >= ncells)
      return;

   /////////////////////////////////////////////
   /// Set local tile & apply boundary conds ///
   /////////////////////////////////////////////

   setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);

   apply_BCs(tile, itile, H, U, V, nlft, nrht, ntop, nbot);

   barrier(CLK_LOCAL_MEM_FENCE);

   /////////////////////////////////////////////////
   /// Declare all constants and local variables ///
   /////////////////////////////////////////////////

   const real   g     = GRAVITATIONAL_CONSTANT;   // gravitational constant
   const real   ghalf = HALF*g;

   // Left, right, ... left-left, right-right, ... left-top, right-top neighbor
   int nl, nr, nt, nb;
   int nll, nrr, ntt, nbb;

   // Level
   int lvl, lvl_nl, lvl_nr, lvl_nt, lvl_nb;
   int lvl_nll, lvl_nrr, lvl_ntt, lvl_nbb;

   // Left-top, right-top, top-right, bottom-right neighbor 
   int nlt, nrt, ntr, nbr;

   // State variables at x-axis control volume face
   real Hxminus, Hxplus;
   real Uxminus, Uxplus;
   real Vxminus, Vxplus;

   // State variables at y-axis control volume face
   real Hyminus, Hyplus;
   real Uyminus, Uyplus;
   real Vyminus, Vyplus;

   // Variables for artificial viscosity/flux limiting
   real wminusx_H, wminusx_U;
   real wplusx_H, wplusx_U;
   real wminusy_H, wminusy_V;
   real wplusy_H, wplusy_V;

   int nltl;
   real Hll2;

   int nrtr;
   real Hrr2;

   real Ull2;
   real Urr2;

   int ntrt;
   real Htt2;

   int nbrb;
   real Hbb2;

   real Vtt2;
   real Vbb2;

   real Hxminus2, Hxplus2;
   real Uxminus2, Uxplus2;
   real Vxminus2, Vxplus2;

   real Hyminus2, Hyplus2;
   real Uyminus2, Uyplus2;
   real Vyminus2, Vyplus2;

   real Hxfluxminus;
   real Uxfluxminus;
   real Vxfluxminus;

   real Hxfluxplus;
   real Uxfluxplus;
   real Vxfluxplus;

   real Hyfluxminus;
   real Uyfluxminus;
   real Vyfluxminus;

   real Hyfluxplus;
   real Uyfluxplus;
   real Vyfluxplus;


   // XXX Assuming square cells! XXX
   // State variables and cell widths and lengths
   real dric, drl, drr, drt, drb;
//   real drlt, drrt, drtr, drbr;

   real Hic, Hl, Hr, Ht, Hb;
   real Hll, Hrr, Htt, Hbb;

   real Uic, Ul, Ur, Ut, Ub;
   real Ull, Urr, Utt, Ubb;

   real Vic, Vl, Vr, Vt, Vb;
   real Vll, Vrr, Vtt, Vbb;

   real Hlt, Hrt, Htr, Hbr;
   real Ult, Urt, Utr, Ubr;
   real Vlt, Vrt, Vtr, Vbr;


   // Local values for the state variables and cell widths and heights for the local cell as well
   // as its neighboring cells
   real dxic, dxl, dxr, dyic, dyt, dyb;

   //////////////////////////
   /// Set the local tile ///
   //////////////////////////

   int start_idx = group_id * ntX;
//   int end_idx = (group_id + 1) * ntX;

   int ic = giX;

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
      Vll     = V[nll];
   }
   else {
      lvl_nll = levelval(nll);
      Hll     = Hval(nll);
      Ull     = Uval(nll);
      Vll     = Vval(nll);
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
      Vrr     = V[nrr];
   }
   else {
      lvl_nrr = levelval(nrr);
      Hrr     = Hval(nrr);
      Urr     = Uval(nrr);
      Vrr     = Vval(nrr);
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
      Utt     = U[ntt];
      Vtt     = V[ntt];
   }
   else {
      lvl_ntt = levelval(ntt);
      Htt     = Hval(ntt);
      Utt     = Uval(ntt);
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
      Ubb     = U[nbb];
      Vbb     = V[nbb];
   }
   else {
      lvl_nbb = levelval(nbb);
      Hbb     = Hval(nbb);
      Ubb     = Uval(nbb);
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

//   real V_stag = V_delta( SQR(dxl), SQR(dxic) );

   Hxminus = U_halfstep(deltaT, Hl, Hic, HXFLUXNL, HXFLUXIC, 
                        dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
   Uxminus = U_halfstep(deltaT, Ul, Uic, UXFLUXNL, UXFLUXIC,
                        dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
   Vxminus = U_halfstep(deltaT, Vl, Vic, UVFLUXNL, UVFLUXIC,
                        dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));

//   V_stag = V_delta( SQR(dxic), SQR(dxr) );

   Hxplus  = U_halfstep(deltaT, Hic, Hr, HXFLUXIC, HXFLUXNR,
                        dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
   Uxplus  = U_halfstep(deltaT, Uic, Ur, UXFLUXIC, UXFLUXNR,
                        dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
   Vxplus  = U_halfstep(deltaT, Vic, Vr, UVFLUXIC, UVFLUXNR,
                        dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));

//   V_stag = V_delta( SQR(dyb), SQR(dyic) );

   Hyminus = U_halfstep(deltaT, Hb, Hic, HYFLUXNB, HYFLUXIC,
                        dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
   Uyminus = U_halfstep(deltaT, Ub, Uic, VUFLUXNB, VUFLUXIC,
                        dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
   Vyminus = U_halfstep(deltaT, Vb, Vic, VYFLUXNB, VYFLUXIC,
                        dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));

//   V_stag = V_delta( SQR(dyic), SQR(dyt) );

   Hyplus  = U_halfstep(deltaT, Hic, Ht, HYFLUXIC, HYFLUXNT,
                        dyic, dyt, dyic, dyt, SQR(dyic), SQR(dyt));
   Uyplus  = U_halfstep(deltaT, Uic, Ut, VUFLUXIC, VUFLUXNT,
                        dyic, dyt, dyic, dyt, SQR(dyic), SQR(dyt));
   Vyplus  = U_halfstep(deltaT, Vic, Vt, VYFLUXIC, VYFLUXNT,
                        dyic, dyt, dyic, dyt, SQR(dyic), SQR(dyt));


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

//      V_stag = V_delta( SQR(drl), SQR(dric) );

      Hxminus2 = U_halfstep(deltaT, Hlt, Hic, HXFLUXNLT, HXFLUXIC,
                            drl, dric, drl, dric, SQR(drl), SQR(dric));
      Uxminus2 = U_halfstep(deltaT, Ult, Uic, UXFLUXNLT, UXFLUXIC,
                            drl, dric, drl, dric, SQR(drl), SQR(dric));
      Vxminus2 = U_halfstep(deltaT, Vlt, Vic, UVFLUXNLT, UVFLUXIC,
                            drl, dric, drl, dric, SQR(drl), SQR(dric));

      Hxfluxminus = (Hxfluxminus + HNEWXFLUXMINUS2) * HALF;
      Uxfluxminus = (Uxfluxminus + UNEWXFLUXMINUS2) * HALF;
      Vxfluxminus = (Vxfluxminus + UVNEWFLUXMINUS2) * HALF;

   }

   if(lvl < lvl_nr) {

//      V_stag = V_delta( SQR(dric), SQR(drr) );

      Hxplus2  = U_halfstep(deltaT, Hic, Hrt, HXFLUXIC, HXFLUXNRT,
                            dric, drr, dric, drr, SQR(dric), SQR(drr));
      Uxplus2  = U_halfstep(deltaT, Uic, Urt, UXFLUXIC, UXFLUXNRT,
                            dric, drr, dric, drr, SQR(dric), SQR(drr));
      Vxplus2  = U_halfstep(deltaT, Vic, Vrt, UVFLUXIC, UVFLUXNRT,
                            dric, drr, dric, drr, SQR(dric), SQR(drr));

      Hxfluxplus  = (Hxfluxplus + HNEWXFLUXPLUS2) * HALF;
      Uxfluxplus  = (Uxfluxplus + UNEWXFLUXPLUS2) * HALF;
      Vxfluxplus  = (Vxfluxplus + UVNEWFLUXPLUS2) * HALF;

   }

   if(lvl < lvl_nb) {

//      V_stag = V_delta( SQR(drb), SQR(dric) );

      Hyminus2 = U_halfstep(deltaT, Hbr, Hic, HYFLUXNBR, HYFLUXIC,
                            drb, dric, drb, dric, SQR(drb), SQR(dric));
      Uyminus2 = U_halfstep(deltaT, Ubr, Uic, VUFLUXNBR, VUFLUXIC,
                            drb, dric, drb, dric, SQR(drb), SQR(dric));
      Vyminus2 = U_halfstep(deltaT, Vbr, Vic, VYFLUXNBR, VYFLUXIC,
                            drb, dric, drb, dric, SQR(drb), SQR(dric));

      Hyfluxminus = (Hyfluxminus + HNEWYFLUXMINUS2) * HALF;
      Uyfluxminus = (Uyfluxminus + VUNEWFLUXMINUS2) * HALF;
      Vyfluxminus = (Vyfluxminus + VNEWYFLUXMINUS2) * HALF;

   }

   if(lvl < lvl_nt) {

//      V_stag = V_delta( SQR(dric), SQR(drt) );

      Hyplus2  = U_halfstep(deltaT, Hic, Htr, HYFLUXIC, HYFLUXNTR,
                            dric, drt, dric, drt, SQR(dric), SQR(drt));
      Uyplus2  = U_halfstep(deltaT, Uic, Utr, VUFLUXIC, VUFLUXNTR,
                            dric, drt, dric, drt, SQR(dric), SQR(drt));
      Vyplus2  = U_halfstep(deltaT, Vic, Vtr, VYFLUXIC, VYFLUXNTR,
                            dric, drt, dric, drt, SQR(dric), SQR(drt));

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

   real Hr2 = Hr;
   real Ur2 = Ur;
   if(lvl < lvl_nr) {
      Hr2 = (Hr2 + Hrt) * HALF;
      Ur2 = (Ur2 + Urt) * HALF;
   }

   real numinusx = fabs(Uxminus/Hxminus) + sqrt(g*Hxminus);

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

   real Hl2 = Hl;
   real Ul2 = Ul;
   if(lvl < lvl_nl) {
      Hl2 = (Hl2 + Hlt) * HALF;
      Ul2 = (Ul2 + Ult) * HALF;
   }

   real nuplusx = fabs(Uxplus/Hxplus) + sqrt(g*Hxplus);

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

   real Ht2 = Ht;
   real Vt2 = Vt;
   if(lvl < lvl_nt) {
      Ht2 = (Ht2 + Htr) * HALF;
      Vt2 = (Vt2 + Vtr) * HALF;
   }

   real numinusy = fabs(Vyminus/Hyminus) + sqrt(g*Hyminus);

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

   real Hb2 = Hb;
   real Vb2 = Vb;
   if(lvl < lvl_nb) {
      Hb2 = (Hb2 + Hbr) * HALF;
      Vb2 = (Vb2 + Vbr) * HALF;
   }

   real nuplusy = fabs(Vyplus/Hyplus) + sqrt(g*Hyplus);

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

#define REAL_CELL 1

__kernel void refine_potential_cl(
                 const int    ncells,   // 0  Total number of cells.
                 const int    levmx,    // 1  Maximum level
        __global const real  *H,        // 2
        __global const real  *U,        // 3
        __global const real  *V,        // 4
        __global const int   *nlft,     // 5  Array of left neighbors.
        __global const int   *nrht,     // 6  Array of right neighbors.
        __global const int   *ntop,     // 7  Array of bottom neighbors.
        __global const int   *nbot,     // 8  Array of top neighbors.
        __global const int   *level,    // 9  Array of level information.
        __global const int   *celltype, // 10  Array of celltype information.
        __global       int   *mpot,     // 11  Array of mesh potential information.
        __global       int   *ioffset,  // 12  Array of new giX offsets.
        __global const real  *lev_dx,   // 13
        __global const real  *lev_dy,   // 14
        __global       int   *result,   // 15
        __local        real4 *tile,     // 16  Tile size in real4.
        __local        int8  *itile)    // 17  Tile size in int8.
{

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
////////////////////                       ///////////////////////
////////////////////   calc_gradients_cl   ///////////////////////
////////////////////                       ///////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

   real qpot = ZERO;
    
   /////////////////////////////////////////////
   /// Get thread identification information ///
   /////////////////////////////////////////////

   const unsigned int giX  = get_global_id(0);
   const unsigned int tiX  = get_local_id(0);

   const unsigned int group_id = get_group_id(0);

   const unsigned int ntX  = get_local_size(0);

   itile[tiX].s0 = 0;

   if(giX >= ncells)
      return;

   setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);

   int cell_add = (celltype[giX] == REAL_CELL) ? 4 : 2;

   barrier (CLK_LOCAL_MEM_FENCE);

   int nlt, nrt, nbr, ntr;
   int nl, nr, nb, nt;
   real Hic, Uic, Vic;
   real Hl, Ul, Vl;
   real Hr, Ur, Vr;
   real Hb, Ub, Vb;
   real Ht, Ut, Vt;
   real dxl, dxr, dyb, dyt;

   real duminus1, duminus2;
   real duplus1, duplus2;
   real duhalf1, duhalf2;

   nl = nlftval(tiX);
   nr = nrhtval(tiX);
   nb = nbotval(tiX);
   nt = ntopval(tiX);

   Hic  = Hval(tiX);
   Uic  = Uval(tiX);
   Vic  = Vval(tiX);

   int lvl = levelval(tiX);

   //////////////////////////
   //////////////////////////
   //////////////////////////

   barrier(CLK_LOCAL_MEM_FENCE);

   //////////////////////////////
   //////////////////////////////
   //////////////////////////////


   // Setting the left and left-left neighbor state values for the control volume based
   // on the state variables of the actual left, left-left, and left-top neighbor cells

   // Using global access for the left neighbor values
   if(nl < 0) {
      nl = abs(nl+1);
      dxl = lev_dx[level[nl]];
      Hl = H[nl];
      Ul = U[nl];
      Vl = V[nl];

      if (level[nl] > lvl) {
         nlt = ntop[nl];
         Hl = HALF * (Hl + H[nlt]);
         Ul = HALF * (Ul + U[nlt]);
         Vl = HALF * (Vl + V[nlt]);
      }
   }
   // Using local access for the left neighbor
   else {
      dxl = lev_dx[level[nl]];
      Hl = Hval(nl);
      Ul = Uval(nl);
      Vl = Vval(nl);

      // The left neighbor is more refined than the current cell
      if (levelval(nl) > lvl) {
         nlt = ntopval(nl);
         if(nlt >= 0) {
            Hl = HALF * (Hl + Hval(nlt));
            Ul = HALF * (Ul + Uval(nlt));
            Vl = HALF * (Vl + Vval(nlt));
         }
         else {
            nlt = abs(nlt+1);
            Hl = HALF * (Hl + H[nlt]);
            Ul = HALF * (Ul + U[nlt]);
            Vl = HALF * (Vl + V[nlt]);
         }
      }
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////

   // Setting the right and right-right neighbor state values for the control volume based
   // on the state variables of the actual right, right-right, and right-top neighbor cells

   // Using global access for the right neighbor values
   if(nr < 0) {
      nr = abs(nr+1);
      dxr = lev_dx[level[nr]] ;
      Hr = H[nr];
      Ur = U[nr];
      Vr = V[nr];
   
      if (level[nr] > lvl) {
         nrt = ntop[nr];
         Hr = HALF * (Hr + H[nrt]);
         Ur = HALF * (Ur + U[nrt]);
         Vr = HALF * (Vr + V[nrt]);
      }
   }
   // Using local access for the right neighbor
   else {
      dxr = lev_dx[level[nr]] ;
      Hr = Hval(nr);
      Ur = Uval(nr);
      Vr = Vval(nr);

      if (levelval(nr) > lvl) {
         nrt = ntopval(nr);
         if(nrt >= 0) {
            Hr = HALF * (Hr + Hval(nrt));
            Ur = HALF * (Ur + Uval(nrt));
            Vr = HALF * (Vr + Vval(nrt));
         }
         else {
            nrt = abs(nrt+1);
            Hr = HALF * (Hr + H[nrt]);
            Ur = HALF * (Ur + U[nrt]);
            Vr = HALF * (Vr + V[nrt]);
         }
      }
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////



   // Setting the bottom and bottom-bottom neighbor state values for the control volume based
   // on the state variables of the actual bottom, bottom-bottom, and bottom-right neighbor cells

   // Using global access for the bottom neighbor values
   if (nb < 0) {
      nb = abs(nb+1);
      dyb = ONE / lev_dy[level[nb]];
      Hb = H[nb];
      Ub = U[nb];
      Vb = V[nb];

      if (level[nb] > lvl) {
         nbr = nrht[nb];
         Hb = HALF * (Hb + H[nbr]);
         Ub = HALF * (Ub + U[nbr]);
         Vb = HALF * (Vb + V[nbr]);
      }
   }
   // Using local access for the bottom neighbor
   else {
      dyb = ONE / lev_dy[levelval(nb)];
      Hb = Hval(nb);
      Ub = Uval(nb);
      Vb = Vval(nb);

      if (levelval(nb) > lvl) {
         nbr = nrhtval(nb);
         if(nbr >= 0) {
            Hb = HALF * (Hb + Hval(nbr));
            Ub = HALF * (Ub + Uval(nbr));
            Vb = HALF * (Vb + Vval(nbr));
         }
         else {
            nbr = abs(nbr+1);
            Hb = HALF * (Hb + H[nbr]);
            Ub = HALF * (Ub + U[nbr]);
            Vb = HALF * (Vb + V[nbr]);
         }
      }
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////


   // Setting the top and top-top neighbor state values for the control volume based
   // on the state variables of the actual top, top-top, and top-right neighbor cells
  
   // Using global access for the top neighbor values
   if (nt < 0) {
      nt = abs(nt+1);
      dyt = ONE / lev_dy[level[nt]];
      Ht = H[nt];
      Ut = U[nt];
      Vt = V[nt];

      if (level[nt] > lvl) {
         ntr = nrht[nt];
         Ht = HALF * (Ht + H[ntr]);
         Ut = HALF * (Ut + U[ntr]);
         Vt = HALF * (Vt + V[ntr]);
      }
   }
   // Using local access for the top neighbor
   else {
      dyt = ONE / lev_dy[levelval(nt)];
      Ht = Hval(nt);
      Ut = Uval(nt);
      Vt = Vval(nt);

      if (levelval(nt) > lvl) {
         ntr = nrhtval(nt);
         if(ntr >= 0) {
            Ht = HALF * (Ht + Hval(ntr));
            Ut = HALF * (Ut + Uval(ntr));
            Vt = HALF * (Vt + Vval(ntr));
         }
         else {
            ntr = abs(ntr+1);
            Ht = HALF * (Ht + H[ntr]);
            Ut = HALF * (Ut + U[ntr]);
            Vt = HALF * (Vt + V[ntr]);
         }
      }
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////

    //--CALCULATIONS------------------------------------------------------------
    //  Calculate the gradient between the right and left neighbors and the
    //  main cell.
    real invHic = ONE / Hic; //  For faster math.
    real qmax = -THOUSAND;          //  Set the default maximum low to catch the real one.
    
    duplus1 = Hr - Hic;     duplus2 = Ur - Uic;
    duhalf1 = Hic - Hl;     duhalf2 = Uic - Ul;
    qpot = max(fabs(duplus1 * invHic), fabs(duhalf1 * invHic));
    if (qpot > qmax) qmax = qpot;
    
    duminus1= Hic - Hl;     duminus2= Uic - Ul;
    duhalf1 = Hr - Hic;     duhalf2 = Ur - Uic;
    qpot = max(fabs(duminus1 * invHic), fabs(duhalf1 * invHic));
    if (qpot > qmax) qmax = qpot;
    
    //  Calculate the gradient between the top and bottom neighbors and the
    //  main cell.
    duplus1 = Ht - Hic;     duplus2 = Vt - Vic;
    duhalf1 = Hic - Hb;     duhalf2 = Vic - Vb;
    qpot = max(fabs(duplus1 * invHic), fabs(duhalf1 * invHic));
    if (qpot > qmax) qmax = qpot;
    
    duminus1= Hic - Hb;     duminus2= Vic - Vb;
    duhalf1 = Ht - Hic;     duhalf2 = Vt - Vic;
    qpot = max(fabs(duminus1 * invHic), fabs(duhalf1 * invHic));
    if (qpot > qmax) qmax = qpot;


    //--CALCULATIONS------------------------------------------------------------
    //  Refine the mesh if the gradient is large
    int mpotval = 0;
    if (qmax > 0.10 && lvl < levmx)  //  XXX:  eliminate hard-coded vars
    {   mpotval = 1; }
    
    itile[tiX].s0 = mpotval ? cell_add : 1;

    barrier(CLK_LOCAL_MEM_FENCE);

    for (int offset = ntX >> 1; offset > 32; offset >>= 1) {
       if (tiX < offset) {
          itile[tiX].s0 += itile[tiX+offset].s0; 
       }
       barrier(CLK_LOCAL_MEM_FENCE);
    }

    //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
    if (tiX < 32)
    {  itile[tiX].s0 += itile[tiX+32].s0;
       itile[tiX].s0 += itile[tiX+16].s0;
       itile[tiX].s0 += itile[tiX+8].s0;
       itile[tiX].s0 += itile[tiX+4].s0;
       itile[tiX].s0 += itile[tiX+2].s0;
       itile[tiX].s0 += itile[tiX+1].s0; }

    if (tiX == 0) {
      ioffset[group_id] = itile[0].s0;
      (*result) = itile[0].s0;
    }

    /////////////////////
    /// GLOBAL WRITES ///
    /////////////////////

    //  Put the mesh potential on the global array.
    mpot[giX] = mpotval; /**/

}

// XXX Refinement Smoothing XXX
__kernel void refine_smooth_cl(
                 const int    ncells,   // 0  Total number of cells.
                 const int    levmx,    // 1  Maximum level
        __global const int   *nlft,     // 2  Array of left neighbors.
        __global const int   *nrht,     // 3  Array of right neighbors.
        __global const int   *ntop,     // 4  Array of bottom neighbors.
        __global const int   *nbot,     // 5  Array of top neighbors.
        __global const int   *level,    // 6  Array of level information.
        __global const int   *celltype, // 7  Array of celltype information.
        __global       int   *mpot,     // 8  Array of mesh potential information.
        __global       int   *ioffset,  // 9  Array of new giX offsets.
        __global       int   *newcount, // 10  Array of number of cells smoothed per tile
        __global       int   *result,   // 11
        __local        int2  *itile)    // 12  Tile size in int2.
{

   const unsigned int giX  = get_global_id(0);
   const unsigned int tiX  = get_local_id(0);

   const unsigned int group_id = get_group_id(0);

   const unsigned int ntX  = get_local_size(0);

   itile[tiX].s0 = 0;
   itile[tiX].s1 = 0;

   if(giX >= ncells)
      return;

//   setup_tile(tile, itile, ncells, H, U, V, nlft, nrht, ntop, nbot, level);

   barrier (CLK_LOCAL_MEM_FENCE);

   int nl, nr, nt, nb;
   int nlt, nrt, ntr, nbr;
   int lev, ll, lr, lt, lb;
   int llt, lrt, ltr, lbr;
   int new_count = 0;

   int ic = giX;

   lev = level[ic];
   if(mpot[ic] > 0) lev++;

   nl = nlft[ic];
   ll = level[nl];
   if(mpot[nl] > 0) ll++;
   
   if(ll - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   nlt = ntop[nl];
   llt = level[nlt];
   if(mpot[nlt] > 0) llt++;

   if(llt - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   nr = nrht[ic];
   lr = level[nr];
   if(mpot[nr] > 0) lr++;
   
   if(lr - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   nrt = ntop[nr];
   lrt = level[nrt];
   if(mpot[nrt] > 0) lrt++;

   if(lrt - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   nt = ntop[ic];
   lt = level[nt];
   if(mpot[nt] > 0) lt++;
   
   if(lt - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   ntr = nrht[nt];
   ltr = level[ntr];
   if(mpot[ntr] > 0) ltr++;

   if(ltr - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   nb = nbot[ic];
   lb = level[nb];
   if(mpot[nb] > 0) lb++;
   
   if(lb - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   nbr = nrht[nb];
   lbr = level[nbr];
   if(mpot[nbr] > 0) lbr++;

   if(lbr - lev > 1) {
      mpot[ic]++;
      new_count++;
      lev = levmx;
   }

   itile[tiX].s1 = new_count;

   // XXX NOTE the global writes to mpot[ic] above... change this eventually XXX
   int cell_add = (celltype[giX] == REAL_CELL) ? 4 : 2;
   itile[tiX].s0 = mpot[ic] ? cell_add : 1;

   barrier(CLK_LOCAL_MEM_FENCE);

   for (int offset = ntX >> 1; offset > 32; offset >>= 1) {
      if (tiX < offset) {
         itile[tiX].s1 += itile[tiX+offset].s1; 
         itile[tiX].s0 += itile[tiX+offset].s0; 
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }

   //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
   if (tiX < 32)
   {  itile[tiX].s1 += itile[tiX+32].s1;
      itile[tiX].s0 += itile[tiX+32].s0;
      itile[tiX].s1 += itile[tiX+16].s1;
      itile[tiX].s0 += itile[tiX+16].s0;
      itile[tiX].s1 += itile[tiX+8].s1;
      itile[tiX].s0 += itile[tiX+8].s0;
      itile[tiX].s1 += itile[tiX+4].s1;
      itile[tiX].s0 += itile[tiX+4].s0;
      itile[tiX].s1 += itile[tiX+2].s1;
      itile[tiX].s0 += itile[tiX+2].s0;
      itile[tiX].s1 += itile[tiX+1].s1;
      itile[tiX].s0 += itile[tiX+1].s0; }

   if (tiX == 0) {
     newcount[group_id] = itile[0].s1;
     ioffset[group_id] = itile[0].s0;
     (*result) = itile[0].s0;
   }

}



/* finish_reduction_scan */

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

__kernel void finish_reduction_scan_cl(
                 const int   isize,
        __global       int  *ioffset,
        __global       int  *result,
        __local        int  *itile_scratch,
        __local        uint  *itile)
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



void setup_tile(__local        real4  *tile, 
                __local        int8   *itile, 
                         const int    isize,
                __global const real   *H,
                __global const real   *U,
                __global const real   *V,
                __global const int    *nlft,
                __global const int    *nrht,
                __global const int    *ntop,
                __global const int    *nbot,
                __global const int    *level
                )
{
   const unsigned int giX = get_global_id (0);
   const unsigned int tiX = get_local_id (0);

   const unsigned int ntX = get_local_size (0);

   const unsigned int group_id = get_group_id (0);

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

void apply_BCs(__local  real4        *tile, 
               __local  int8         *itile, 
               __global const real   *H,
               __global const real   *U,
               __global const real   *V,
               __global const int    *nlft,
               __global const int    *nrht,
               __global const int    *ntop,
               __global const int    *nbot)
{

   int nr, nl, nt, nb;

   const unsigned int giX = get_global_id (0);
   const unsigned int tiX = get_local_id (0);

   const unsigned int ntX = get_local_size (0);

   const unsigned int group_id = get_group_id (0);

   int start_idx = group_id * ntX;

   ////////////////////////////////////////////////////////////////////////////////////////////
   ///                  Setting the local values for the neighbors                          ///
   ///        If a global read is necessary, the value of neighbor is negated.              ///
   /// If it's a local read, the value is mapped by the relative offset to the start index. ///
   ////////////////////////////////////////////////////////////////////////////////////////////

   //sets left boundary conditions

   int lft_bound = 0;

   // Test for global index -- it will be negative if so
   if (nlftval(tiX) < 0) {
      // set left boundary value if equal to self
      if (abs(nlftval(tiX)+1) == giX) lft_bound = 1;
   } else {
      // Check local index, but add offset of start_idx
      if (nlftval(tiX) + start_idx == giX) lft_bound = 1;
   }

   if (lft_bound) {
      // if left boundary get right neigbor to update boundary value from
      nr = nrhtval(tiX);

      // Checking for global index or local -- negative is global
      if (nr < 0) {
        // Copy from global data
        Hval(tiX) =  H[abs(nr+1)];
        Uval(tiX) = -U[abs(nr+1)];
        Vval(tiX) =  V[abs(nr+1)];
      } else {
        // Copy from local data
        Hval(tiX) =  Hval(nr);
        Uval(tiX) = -Uval(nr);
        Vval(tiX) =  Vval(nr);
      }
   }

   //sets right boundary conditions

   int rht_bound = 0;

   // Test for global index -- it will be negative if so
   if (nrhtval(tiX) < 0) {
      // set left boundary value if equal to self
      if (abs(nrhtval(tiX)+1) == giX) rht_bound = 1;
   } else {
      // Check local index, but add offset of start_idx
      if (nrhtval(tiX) + start_idx == giX) rht_bound = 1;
   }

   if (rht_bound) {
      nl = nlftval(tiX);

      if (nl < 0) {
        // Copy from global data
         Hval(tiX) =  H[abs(nl+1)];
         Uval(tiX) = -U[abs(nl+1)];
         Vval(tiX) =  V[abs(nl+1)];
      } else {
        // Copy from local data
         Hval(tiX) =  Hval(nl);
         Uval(tiX) = -Uval(nl);
         Vval(tiX) =  Vval(nl);
      }
   }

   //sets bottom boundary conditions

   int bot_bound = 0;

   // Test for global index -- it will be negative if so
   if (nbotval(tiX) < 0) {
      // set left boundary value if equal to self
      if (abs(nbotval(tiX)+1) == giX) bot_bound = 1;
   } else {
      // Check local index, but add offset of start_idx
      if (nbotval(tiX) + start_idx == giX) bot_bound = 1;
   }

   if (bot_bound) {
      nt = ntopval(tiX);

      if (nt < 0) {
        // Copy from global data
         Hval(tiX) =  H[abs(nt+1)];
         Uval(tiX) =  U[abs(nt+1)];
         Vval(tiX) = -V[abs(nt+1)];
      } else {
        // Copy from local data
         Hval(tiX) =  Hval(nt);
         Uval(tiX) =  Uval(nt);
         Vval(tiX) = -Vval(nt);
      }
   }

   //sets top boundary conditions

   int top_bound = 0;

   if (ntopval(tiX) < 0) {
      if (abs(ntopval(tiX)+1) == giX) top_bound = 1;
   } else {
      if (ntopval(tiX) + start_idx == giX) top_bound = 1;
   }

   if (top_bound) {
      nb = nbotval(tiX);

      if (nb < 0) {
         // Copy from global data
         Hval(tiX) =  H[abs(nb+1)];
         Uval(tiX) =  U[abs(nb+1)];
         Vval(tiX) = -V[abs(nb+1)];
      } else {
         // Copy from local data
         Hval(tiX) =  Hval(nb);
         Uval(tiX) =  Uval(nb);
         Vval(tiX) = -Vval(nb);
      }
   }
}

