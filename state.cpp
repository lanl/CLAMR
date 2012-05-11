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
#include <unistd.h>
#include <stdio.h>
#include <mpi.h>
#include "mesh.h"
#include "state.h"
#include "kdtree/KDTree.h"
#include "reorder.h"
#include "timer/timer.h"

#ifdef HAVE_CL_DOUBLE
typedef double      real;
typedef struct
{
   double s0;
   double s1;
}  real2;
typedef cl_double2  cl_real2;
#define L7_REAL L7_DOUBLE
#else
typedef struct
{
   float s0;
   float s1;
}  real2;
typedef cl_float2   cl_real2;
#define L7_REAL L7_FLOAT
#endif

typedef unsigned int uint;

struct esum_type{
   double sum;
   double correction;
};
MPI_Datatype MPI_TWO_DOUBLES;
MPI_Op KAHAN_SUM;
int commutative = 1;
void kahan_sum(struct esum_type *in, struct esum_type *inout, int *len, MPI_Datatype *MPI_TWO_DOUBLES);

int save_ncells;

#define CONSERVED_EQNS
#define REFINE_GRADIENT 0.10

#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define SQR(x) ( x*x )
#define EPSILON 1.0e-30
#define MIN3(x,y,z) ( min( min(x,y), z) )

cl_kernel kernel_set_timestep;
cl_kernel kernel_reduction_min;
cl_kernel kernel_copy_state_data;
cl_kernel kernel_copy_state_ghost_data;
cl_kernel kernel_calc_finite_difference;
cl_kernel kernel_refine_potential;
cl_kernel kernel_refine_smooth;
cl_kernel kernel_rezone_all;
cl_kernel kernel_reduce_sum_mass_stage1of2;
cl_kernel kernel_reduce_sum_mass_stage2of2;
cl_kernel kernel_reduce_epsum_mass_stage1of2;
cl_kernel kernel_reduce_epsum_mass_stage2of2;

real U_halfstep(// XXX Fix the subindices to be more intuitive XXX
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

   return (( r_i*U_n + r_n*U_i ) / ( r_i + r_n )) 
          - HALF*deltaT*(( F_n*A_n*min((real)ONE, A_i/A_n) - F_i*A_i*min((real)ONE, A_n/A_i) )
                    / ( V_n*min((real)HALF, V_i/V_n) + V_i*min((real)HALF, V_n/V_i) ));

}

real U_fullstep(
        real    deltaT,
        real    dr,
        real    U,
        real    F_plus,
        real    F_minus,
        real    G_plus,
        real    G_minus) {

   return (U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus));

}


real w_corrector(
        real    deltaT,       // Timestep
        real    dr,           // Cell's center to face distance
        real    U_eigen,      // State variable's eigenvalue (speed)
        real    grad_half,    // Centered gradient
        real    grad_minus,   // Downwind gradient
        real    grad_plus) {  // Upwind gradient

   real nu     = HALF * U_eigen * deltaT / dr;
   nu          = nu * (ONE - nu);

   real rdenom = ONE / max(SQR(grad_half), (real)EPSILON);
   real rplus  = (grad_plus  * grad_half) * rdenom;
   real rminus = (grad_minus * grad_half) * rdenom;

   return HALF*nu*(ONE- max(MIN3((real)ONE, rplus, rminus), (real)ZERO));

}


State::State(int ncells, cl_context context)
{
   cpu_time_apply_BCs          = 0;
   cpu_time_set_timestep       = 0;
   cpu_time_finite_difference  = 0;
   cpu_time_refine_potential   = 0;
   cpu_time_rezone_all         = 0;
   cpu_time_mass_sum           = 0;

   gpu_time_apply_BCs          = 0;
   gpu_time_set_timestep       = 0;
   gpu_time_finite_difference  = 0;
   gpu_time_refine_potential   = 0;
   gpu_time_rezone_all         = 0;
   gpu_time_mass_sum           = 0;
   gpu_time_read               = 0;
   gpu_time_write              = 0;
}

void State::init(int ncells, cl_context context, int do_gpu_calc)
{
   if (do_gpu_calc) {
      kernel_set_timestep      = ezcl_create_kernel(context, "wave_kern.cl",      "set_timestep_cl",    0);
      kernel_reduction_min     = ezcl_create_kernel(context, "wave_kern.cl",      "finish_reduction_min_cl",  0);
      kernel_copy_state_data   = ezcl_create_kernel(context, "wave_kern_calc.cl", "copy_state_data_cl",  0);
      kernel_copy_state_ghost_data   = ezcl_create_kernel(context, "wave_kern_calc.cl", "copy_state_ghost_data_cl",  0);
      kernel_calc_finite_difference = ezcl_create_kernel(context, "wave_kern_calc.cl", "calc_finite_difference_cl",  0);
      kernel_refine_potential  = ezcl_create_kernel(context, "wave_kern_calc.cl", "refine_potential_cl", 0);
      kernel_refine_smooth     = ezcl_create_kernel(context, "wave_kern_calc.cl", "refine_smooth_cl", 0);
      kernel_rezone_all        = ezcl_create_kernel(context, "wave_kern.cl",      "rezone_all_cl",      0);
      kernel_reduce_sum_mass_stage1of2 = ezcl_create_kernel(context, "wave_kern.cl", "reduce_sum_mass_stage1of2_cl",  0);
      kernel_reduce_sum_mass_stage2of2 = ezcl_create_kernel(context, "wave_kern.cl", "reduce_sum_mass_stage2of2_cl",  0);
      kernel_reduce_epsum_mass_stage1of2 = ezcl_create_kernel(context, "wave_kern.cl", "reduce_epsum_mass_stage1of2_cl",  0);
      kernel_reduce_epsum_mass_stage2of2 = ezcl_create_kernel(context, "wave_kern.cl", "reduce_epsum_mass_stage2of2_cl",  0);
   }

   H.resize(ncells);
   U.resize(ncells);
   V.resize(ncells);

   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init){
      MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_TWO_DOUBLES);
      MPI_Type_commit(&MPI_TWO_DOUBLES);
      MPI_Op_create((MPI_User_function *)kahan_sum, commutative, &KAHAN_SUM);
   }
}

void kahan_sum(struct esum_type *in, struct esum_type *inout, int *len, MPI_Datatype *MPI_TWO_DOUBLES)
{
   double corrected_next_term, new_sum;

   corrected_next_term = in->sum +(in->correction+inout->correction);
   new_sum = inout->sum + corrected_next_term;
   inout->correction = corrected_next_term - (new_sum - inout->sum);
   inout->sum = new_sum;
}

void State::add_boundary_cells(Mesh *mesh)
{
   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   // This is for a mesh with no boundary cells -- they are added and
   // the mesh sizes increased
   size_t &ncells        = mesh->ncells;
   vector<int>  &i        = mesh->i;
   vector<int>  &j        = mesh->j;
   vector<int>  &level    = mesh->level;
   vector<int>  &index    = mesh->index;
   vector<int>  &celltype = mesh->celltype;
   vector<int>  &nlft     = mesh->nlft;
   vector<int>  &nrht     = mesh->nrht;
   vector<int>  &nbot     = mesh->nbot;
   vector<int>  &ntop     = mesh->ntop;
   vector<real> &x        = mesh->x;
   vector<real> &dx       = mesh->dx;
   vector<real> &y        = mesh->y;
   vector<real> &dy       = mesh->dy;

   vector<int> &lev_ibegin = mesh->lev_ibegin;
   vector<int> &lev_iend   = mesh->lev_iend;
   vector<int> &lev_jbegin = mesh->lev_jbegin;
   vector<int> &lev_jend   = mesh->lev_jend;

   // Pre-count number of cells to add
   int icount = 0;
   for (uint ic=0; ic<ncells; ic++) {
      if (i[ic] == lev_ibegin[level[ic]]) icount++; // Left boundary
      if (i[ic] == lev_iend[level[ic]])   icount++; // Right boundary
      if (j[ic] == lev_jbegin[level[ic]]) icount++; // Bottom boundary
      if (j[ic] == lev_jend[level[ic]])   icount++; // Top boundary
   }
      
   int new_ncells = ncells + icount;
   // Increase the arrays for the new boundary cells
   H.resize(new_ncells);
   U.resize(new_ncells);
   V.resize(new_ncells);
   i.resize(new_ncells);
   j.resize(new_ncells);
   nlft.resize(new_ncells);
   nrht.resize(new_ncells);
   nbot.resize(new_ncells);
   ntop.resize(new_ncells);
   level.resize(new_ncells);
   celltype.resize(new_ncells);
   index.resize(new_ncells);
   x.resize(new_ncells);
   dx.resize(new_ncells);
   y.resize(new_ncells);
   dy.resize(new_ncells);
   for (int nc=ncells; nc<new_ncells; nc++) {
      nlft[nc] = -1;
      nrht[nc] = -1;
      nbot[nc] = -1;
      ntop[nc] = -1;
   }
      
   // In the first pass, set two of the neighbor indices and all
   // the other data to be brought across. Set the inverse of the
   // the velocity to enforce the reflective boundary condition
   uint nc=ncells;
   for (uint ic=0; ic<ncells; ic++) {
      if (i[ic] == lev_ibegin[level[ic]]) {
         nlft[ic] = nc;
         nlft[nc] = nc;
         nrht[nc] = ic;
         i[nc] = lev_ibegin[level[ic]]-1;
         j[nc] = j[ic];
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic]-dx[ic];
         y[nc] = y[ic];
         H[nc] =  H[ic];
         U[nc] = -U[ic];
         V[nc] =  V[ic];
         nc++;
      }
      if (i[ic] == lev_iend[level[ic]]) {
         nrht[ic] = nc;
         nrht[nc] = nc;
         nlft[nc] = ic;
         i[nc] = lev_iend[level[ic]]+1;
         j[nc] = j[ic];
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic]+dx[ic];
         y[nc] = y[ic];
         H[nc] =  H[ic];
         U[nc] = -U[ic];
         V[nc] =  V[ic];
         nc++;
      }
      if (j[ic] == lev_jbegin[level[ic]]) {
         nbot[ic] = nc;
         nbot[nc] = nc;
         ntop[nc] = ic;
         i[nc] = i[ic];
         j[nc] = lev_jbegin[level[ic]]-1;
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic];
         y[nc] = y[ic]-dy[ic];
         H[nc] =  H[ic];
         U[nc] =  U[ic];
         V[nc] = -V[ic];
         nc++;
      }
      if (j[ic] == lev_jend[level[ic]]) {
         ntop[ic] = nc;
         ntop[nc] = nc;
         nbot[nc] = ic;
         i[nc] = i[ic];
         j[nc] = lev_jend[level[ic]]+1;
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic];
         y[nc] = y[ic]+dy[ic];
         H[nc] =  H[ic];
         U[nc] =  U[ic];
         V[nc] = -V[ic];
         nc++;
      }
   }

   // Now set the other two neighbor indices
   for (int nc=ncells; nc<new_ncells; nc++) {
      if (i[nc] == lev_ibegin[level[nc]]-1) {
         // Need to check if also a bottom boundary cell
         if (j[nc] == lev_jbegin[level[nc]]){
           nbot[nc] = nc;
         } else {
           nbot[nc] = nlft[nbot[nrht[nc]]];
         }
         if (j[nc] == lev_jend[level[nc]]){
           ntop[nc] = nc;
         } else {
           ntop[nc] = nlft[ntop[nrht[nc]]];
         }
      }
      if (i[nc] == lev_iend[level[nc]]+1)   {
         if (level[nc] <= level[nbot[nlft[nc]]]){
            if (j[nc] == lev_jbegin[level[nc]]){
               nbot[nc] = nc;
            } else {
               nbot[nc] = nrht[nbot[nlft[nc]]];
            }
            if (j[nc] == lev_jend[level[nc]]){
               ntop[nc] = nc;
            } else {
               ntop[nc] = nrht[ntop[nlft[nc]]];
            }
         // calculation is a little different if going through a
         // finer zoned region
         } else {
            nbot[nc] = nrht[nrht[nbot[nlft[nc]]]];
            ntop[nc] = nrht[nrht[ntop[nlft[nc]]]];
         }
      }
      if (j[nc] == lev_jbegin[level[nc]]-1) {
         if (i[nc] == lev_ibegin[level[nc]]){
            nlft[nc] = nc;
         } else {
            nlft[nc] = nbot[nlft[ntop[nc]]];
         }
         if (i[nc] == lev_iend[level[nc]]){
            nrht[nc] = nc;
         } else {
            nrht[nc] = nbot[nrht[ntop[nc]]];
         }
      }
      if (j[nc] == lev_jend[level[nc]]+1)   {
         if (level[nc] <= level[nlft[nbot[nc]]]){
            if (i[nc] == lev_ibegin[level[nc]]){
               nlft[nc] = nc;
            } else {
               nlft[nc] = ntop[nlft[nbot[nc]]];
            }
            if (i[nc] == lev_iend[level[nc]]){
               nrht[nc] = nc;
            } else {
               nrht[nc] = ntop[nrht[nbot[nc]]];
            }
         } else {
            nlft[nc] = ntop[ntop[nlft[nbot[nc]]]];
            nrht[nc] = ntop[ntop[nrht[nbot[nc]]]];
         }
      }
   }
   save_ncells = ncells;
   ncells = new_ncells;

   cpu_time_apply_BCs += cpu_timer_stop(tstart_cpu);
}

void State::apply_boundary_conditions(Mesh *mesh)
{
   struct timeval tstart_cpu;
   int nl, nr, nb, nt;

   size_t &ncells = mesh->ncells;
   vector<int> &nlft = mesh->nlft;
   vector<int> &nrht = mesh->nrht;
   vector<int> &nbot = mesh->nbot;
   vector<int> &ntop = mesh->ntop;

   cpu_timer_start(&tstart_cpu);

   // This is for a mesh with boundary cells
   for (uint ic=0; ic<ncells; ic++) {
      if (mesh->is_left_boundary(ic)) {
         nr = nrht[ic];
         H[ic] =  H[nr];
         U[ic] = -U[nr];
         V[ic] =  V[nr];
      }
      if (mesh->is_right_boundary(ic))  {
         nl = nlft[ic];
         H[ic] =  H[nl];
         U[ic] = -U[nl];
         V[ic] =  V[nl];
      }
      if (mesh->is_bottom_boundary(ic)) {
         nt = ntop[ic];
         H[ic] =  H[nt];
         U[ic] =  U[nt];
         V[ic] = -V[nt];
      }
      if (mesh->is_top_boundary(ic)) {
         nb = nbot[ic];
         H[ic] =  H[nb];
         U[ic] =  U[nb];
         V[ic] = -V[nb];
      }
   }

   cpu_time_apply_BCs += cpu_timer_stop(tstart_cpu);
}

void State::remove_boundary_cells(Mesh *mesh)
{
   size_t &ncells = mesh->ncells;
   vector<int> &nlft     = mesh->nlft;
   vector<int> &nrht     = mesh->nrht;
   vector<int> &nbot     = mesh->nbot;
   vector<int> &ntop     = mesh->ntop;
   vector<int> &level    = mesh->level;
   vector<int> &i        = mesh->i;
   vector<int> &j        = mesh->j;
   vector<int> &celltype = mesh->celltype;
   vector<int> &index    = mesh->index;
   vector<real> &x       = mesh->x;
   vector<real> &dx      = mesh->dx;
   vector<real> &y       = mesh->y;
   vector<real> &dy      = mesh->dy;

   if(mesh->have_boundary) return;

   // Resize to drop all the boundary cells
   ncells = save_ncells;
   H.resize(save_ncells);
   U.resize(save_ncells);
   V.resize(save_ncells);
   nlft.resize(save_ncells);
   nrht.resize(save_ncells);
   nbot.resize(save_ncells);
   ntop.resize(save_ncells);
   level.resize(save_ncells);
   i.resize(save_ncells);
   j.resize(save_ncells);
   celltype.resize(save_ncells);
   index.resize(save_ncells);
   x.resize(save_ncells);
   dx.resize(save_ncells);
   y.resize(save_ncells);
   dy.resize(save_ncells);

   // Reset the neighbors due to the dropped boundary cells
   for (uint ic=0; ic<ncells; ic++) {
      if (i[ic] == mesh->lev_ibegin[level[ic]]) nlft[ic] = ic;
      if (i[ic] == mesh->lev_iend[level[ic]])   nrht[ic] = ic;
      if (j[ic] == mesh->lev_jbegin[level[ic]]) nbot[ic] = ic;
      if (j[ic] == mesh->lev_jend[level[ic]])   ntop[ic] = ic;
   }

}

double State::set_timestep(Mesh *mesh, double g, double sigma)
{
   int lev;
   double wavespeed, xspeed, yspeed;
   double deltaT, globalmindeltaT;
   double mindeltaT = 1000.0;
   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   size_t &ncells        = mesh->ncells;
   int &parallel         = mesh->parallel;
   vector<int> &celltype = mesh->celltype;
   vector<int> &level    = mesh->level;

   for (uint ic=0; ic<ncells; ic++) {
      if (celltype[ic] == REAL_CELL) {
         lev = level[ic];
         wavespeed = sqrt(g*H[ic]);
         xspeed = (fabs(U[ic])+wavespeed)/mesh->lev_deltax[lev];
         yspeed = (fabs(V[ic])+wavespeed)/mesh->lev_deltay[lev];
         deltaT=sigma/(xspeed+yspeed);
         if (deltaT < mindeltaT) mindeltaT = deltaT;
      }
   }

   if (parallel) MPI_Allreduce(&mindeltaT, &globalmindeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

   cpu_time_set_timestep += cpu_timer_stop(tstart_cpu);

   if (parallel) {
      return(globalmindeltaT);
   } else {
      return(mindeltaT);
   }
}

double State::gpu_set_timestep(cl_command_queue command_queue, Mesh *mesh, double sigma)
{
   cl_event set_timestep_event;
   cl_event reduction_min_event;
   cl_event start_read_event;

   double deltaT, globalmindeltaT;

   size_t &ncells       = mesh->ncells;
   int &parallel        = mesh->parallel;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_level);
   assert(dev_celltype);
   assert(dev_levdx);
   assert(dev_levdy);

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   cl_mem dev_redscratch = ezcl_malloc(NULL, &block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

      /*
      __kernel void set_timestep_cl(
                       const int    ncells,     // 0  Total number of cells.
                       const real   sigma,      // 1
              __global const real  *H,          // 2
              __global const real  *U,          // 3
              __global const real  *V,          // 4
              __global const int   *level,      // 5  Array of level information.
              __global const int   *celltype,   // 6
              __global const real  *lev_dx,     // 7
              __global const real  *lev_dy,     // 8
              __global       real  *redscratch, // 9
              __global       real  *deltaT,     // 10
              __local        real  *tile)       // 11
      */

   real sigma_local = sigma;
   ezcl_set_kernel_arg(kernel_set_timestep,  0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_set_timestep,  1, sizeof(cl_real), (void *)&sigma_local);
   ezcl_set_kernel_arg(kernel_set_timestep,  2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_set_timestep,  3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_set_timestep,  4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_set_timestep,  5, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_set_timestep,  6, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_set_timestep,  7, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_set_timestep,  8, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_set_timestep,  9, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_set_timestep, 10, sizeof(cl_mem),  (void *)&dev_deltaT);
   ezcl_set_kernel_arg(kernel_set_timestep, 11, local_work_size*sizeof(cl_real),  NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_set_timestep, 1, NULL, &global_work_size, &local_work_size, &set_timestep_event);

      /*
      __kernel void finish_reduction_min_cl(
        const    int    isize,
        __global real  *redscratch,
        __global real  *deltaT,
        __local  real  *tile)
      */
   ezcl_set_kernel_arg(kernel_reduction_min, 0, sizeof(cl_int),  (void *)&block_size);
   ezcl_set_kernel_arg(kernel_reduction_min, 1, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_reduction_min, 2, sizeof(cl_mem),  (void *)&dev_deltaT);
   ezcl_set_kernel_arg(kernel_reduction_min, 3, local_work_size*sizeof(cl_real), NULL);

   if (block_size > 1){
     ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_min, 1, NULL, &local_work_size, &local_work_size, &reduction_min_event);
   }
   real deltaT_local;
   ezcl_enqueue_read_buffer(command_queue, dev_deltaT, CL_TRUE,  0, sizeof(cl_real), &deltaT_local, &start_read_event);
   deltaT = deltaT_local;

   if (parallel) MPI_Allreduce(&deltaT, &globalmindeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

   ezcl_device_memory_remove(dev_redscratch);

   gpu_time_set_timestep      += ezcl_timer_calc(&set_timestep_event,  &set_timestep_event); 
   if (block_size > 1) {
      gpu_time_set_timestep   += ezcl_timer_calc(&reduction_min_event, &reduction_min_event);
   }
   gpu_time_read              += ezcl_timer_calc(&start_read_event,    &start_read_event);

   if (parallel) {
      return(globalmindeltaT);
   } else {
      return(deltaT);
   }
}

void State::fill_circle(Mesh   *mesh,       //  Mesh.
                        double  circ_radius,//  Radius of circle in grid units.
                        double  fill_value, //  Circle height for shallow water.
                        double  background) //  Background height for shallow water.
{  
   size_t &ncells = mesh->ncells;
   vector<real> &x  = mesh->x;
   vector<real> &dx = mesh->dx;
   vector<real> &y  = mesh->y;
   vector<real> &dy = mesh->dy;

   for (uint ic = 0; ic < ncells; ic++)
   {  H[ic] = background;
      U[ic] = V[ic] = 0.0; }
   
   //  Clear the old k-D tree and generate new data (slow but necessary here).
   KDTree_Destroy(&mesh->tree);
   mesh->kdtree_setup();
   
   int nez;
   vector<int>    ind(ncells);
   vector<double> weight(ncells);
   
   KDTree_QueryCircleInterior(&mesh->tree, &nez, &(ind[0]), circ_radius, ncells,
                              &x[0], &dx[0],
                              &y[0], &dy[0]);
   for (int ic = 0; ic < nez; ++ic)
   {  H[ind[ic]] = fill_value; }
   
   KDTree_QueryCircleIntersectWeighted(&mesh->tree, &nez, &(ind[0]), &(weight[0]),
                              circ_radius, ncells,
                              &x[0], &dx[0],
                              &y[0], &dy[0]);
   
   for (int ic = 0; ic < nez; ++ic)
   {  H[ind[ic]] = background + (fill_value - background) * weight[ic]; }
}

void State::state_reorder(vector<int> iorder)
{
   reorder(H, iorder);
   reorder(U, iorder);
   reorder(V, iorder);
}

void State::rezone_all(Mesh *mesh, vector<int> mpot, int add_ncells)
{
   int ic, nc;
   struct timeval tstart_cpu;

   vector<int> &celltype = mesh->celltype;
   //vector<int> &nlft     = mesh->nlft;
   vector<int> &nrht     = mesh->nrht;
   //vector<int> &nbot     = mesh->nbot;
   vector<int> &ntop     = mesh->ntop;

   int ncells = mesh->ncells;

   vector<int> celltype_save(ncells);
   for (int ic=0; ic < ncells; ic++){
      celltype_save[ic] = celltype[ic];
   }

   mesh->rezone_all(mpot, add_ncells);

   if (add_ncells == 0 ) return;

   cpu_timer_start(&tstart_cpu);

   int new_ncells = ncells + add_ncells;

   vector<real>H_old(new_ncells);
   vector<real>U_old(new_ncells);
   vector<real>V_old(new_ncells);
   H.swap(H_old);
   U.swap(U_old);
   V.swap(V_old);

   for (ic=0, nc=0; ic<ncells; ic++) {

      if (mpot[ic] == 0) {
         H[nc] = H_old[ic];
         U[nc] = U_old[ic];
         V[nc] = V_old[ic];
         nc++;
      } else if (mpot[ic] < 0){
         int nr = nrht[ic];
         int nt = ntop[ic];
         int nrt = nrht[nt];
         H[nc] = (H_old[ic] + H_old[nr] + H_old[nt] + H_old[nrt])*0.25;
         U[nc] = (U_old[ic] + U_old[nr] + U_old[nt] + U_old[nrt])*0.25;
         V[nc] = (V_old[ic] + V_old[nr] + V_old[nt] + V_old[nrt])*0.25;
         nc++;

      } else if (mpot[ic] > 0){
         // lower left
         H[nc] = H_old[ic];
         U[nc] = U_old[ic];
         V[nc] = V_old[ic];
         nc++;

         // lower right
         H[nc] = H_old[ic];
         U[nc] = U_old[ic];
         V[nc] = V_old[ic];
         nc++;

         if (celltype_save[ic] == REAL_CELL){
            // upper left
            H[nc] = H_old[ic];
            U[nc] = U_old[ic];
            V[nc] = V_old[ic];
            nc++;

            // upper right
            H[nc] = H_old[ic];
            U[nc] = U_old[ic];
            V[nc] = V_old[ic];
            nc++;
         }

      }

   }

   cpu_time_rezone_all += cpu_timer_stop(tstart_cpu);
}

void State::rezone_all_local(Mesh *mesh, vector<int> mpot, int add_ncells)
{
   int ic, nc;
   struct timeval tstart_cpu;

   vector<int> &celltype = mesh->celltype;
   //vector<int> &nlft     = mesh->nlft;
   vector<int> &nrht     = mesh->nrht;
   //vector<int> &nbot     = mesh->nbot;
   vector<int> &ntop     = mesh->ntop;

   int ncells = mesh->ncells;

   vector<int> celltype_save(ncells);
   for (int ic=0; ic < ncells; ic++){
      celltype_save[ic] = celltype[ic];
   }

   mesh->rezone_all(mpot, add_ncells);

   int global_add_ncells;
   MPI_Allreduce(&add_ncells, &global_add_ncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (global_add_ncells == 0 ) return;

   cpu_timer_start(&tstart_cpu);

   int new_ncells = ncells + add_ncells;

   vector<real>H_old(new_ncells);
   vector<real>U_old(new_ncells);
   vector<real>V_old(new_ncells);
   H.swap(H_old);
   U.swap(U_old);
   V.swap(V_old);

   for (ic=0, nc=0; ic<ncells; ic++) {

      if (mpot[ic] == 0) {
         H[nc] = H_old[ic];
         U[nc] = U_old[ic];
         V[nc] = V_old[ic];
         nc++;
      } else if (mpot[ic] < 0){
         int nr = nrht[ic];
         int nt = ntop[ic];
         int nrt = nrht[nt];
         H[nc] = (H_old[ic] + H_old[nr] + H_old[nt] + H_old[nrt])*0.25;
         U[nc] = (U_old[ic] + U_old[nr] + U_old[nt] + U_old[nrt])*0.25;
         V[nc] = (V_old[ic] + V_old[nr] + V_old[nt] + V_old[nrt])*0.25;
         nc++;

      } else if (mpot[ic] > 0){
         // lower left
         H[nc] = H_old[ic];
         U[nc] = U_old[ic];
         V[nc] = V_old[ic];
         nc++;

         // lower right
         H[nc] = H_old[ic];
         U[nc] = U_old[ic];
         V[nc] = V_old[ic];
         nc++;

         if (celltype_save[ic] == REAL_CELL){
            // upper left
            H[nc] = H_old[ic];
            U[nc] = U_old[ic];
            V[nc] = V_old[ic];
            nc++;

            // upper right
            H[nc] = H_old[ic];
            U[nc] = U_old[ic];
            V[nc] = V_old[ic];
            nc++;
         }

      }

   }

   cpu_time_rezone_all += cpu_timer_stop(tstart_cpu);
}


void State::gpu_rezone_all(cl_command_queue command_queue, Mesh *mesh, size_t &ncells, size_t new_ncells, size_t old_ncells, bool localStencil,
      cl_mem dev_mpot, cl_mem dev_ioffset)
{
   cl_event rezone_all_event;

   cl_mem dev_ptr = NULL;

   cl_mem &dev_level        = mesh->dev_level;
   cl_mem &dev_i            = mesh->dev_i;
   cl_mem &dev_j            = mesh->dev_j;
   cl_mem &dev_celltype     = mesh->dev_celltype;
   cl_mem &dev_level_new    = mesh->dev_level_new;
   cl_mem &dev_i_new        = mesh->dev_i_new;
   cl_mem &dev_j_new        = mesh->dev_j_new;
   cl_mem &dev_celltype_new = mesh->dev_celltype_new;
   cl_mem &dev_levdx        = mesh->dev_levdx;
   cl_mem &dev_levdy        = mesh->dev_levdy;

   assert(dev_mpot);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);
   assert(dev_celltype);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_level_new);
   assert(dev_i_new);
   assert(dev_j_new);
   assert(dev_celltype_new);
   assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   if (new_ncells != old_ncells){
      ncells = new_ncells;
      mesh->resize_new_device_memory(ncells);
   }

   cl_mem dev_H_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_U_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_V_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

   cl_event start_write_event;
   cl_mem dev_ijadd;

   vector<int>ijadd(6);
   ijadd[0] = 0;
   ijadd[1] = 0;
   ijadd[2] = 0;
   ijadd[3] = 0;
   ijadd[4] = 0;
   ijadd[5] = 0;

   size_t six = 6;
   dev_ijadd = ezcl_malloc(NULL, &six, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
   ezcl_enqueue_write_buffer(command_queue, dev_ijadd, CL_TRUE, 0, 6*sizeof(cl_int), (void*)&ijadd[0], &start_write_event);

   int stencil = 0;
   if (localStencil) stencil = 1;

   size_t local_work_size = 128;
   size_t global_work_size = ((old_ncells+local_work_size - 1) /local_work_size) * local_work_size;

   ezcl_set_kernel_arg(kernel_rezone_all, 0,  sizeof(cl_int),  (void *)&old_ncells);
   ezcl_set_kernel_arg(kernel_rezone_all, 1,  sizeof(cl_int),  (void *)&stencil);
   ezcl_set_kernel_arg(kernel_rezone_all, 2,  sizeof(cl_int),  (void *)&mesh->levmx);
   ezcl_set_kernel_arg(kernel_rezone_all, 3,  sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_rezone_all, 4,  sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_rezone_all, 5,  sizeof(cl_mem),  (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_rezone_all, 6,  sizeof(cl_mem),  (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_rezone_all, 7,  sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_rezone_all, 8,  sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_rezone_all, 9,  sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_rezone_all, 10, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_rezone_all, 11, sizeof(cl_mem),  (void *)&dev_level_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 12, sizeof(cl_mem),  (void *)&dev_i_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 13, sizeof(cl_mem),  (void *)&dev_j_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 14, sizeof(cl_mem),  (void *)&dev_celltype_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 15, sizeof(cl_mem),  (void *)&dev_H_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 16, sizeof(cl_mem),  (void *)&dev_U_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 17, sizeof(cl_mem),  (void *)&dev_V_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 18, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_rezone_all, 19, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_rezone_all, 20, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_rezone_all, 21, sizeof(cl_mem),  (void *)&mesh->dev_levtable);
   ezcl_set_kernel_arg(kernel_rezone_all, 22, sizeof(cl_mem),  (void *)&dev_ijadd);
   ezcl_set_kernel_arg(kernel_rezone_all, 23, local_work_size * sizeof(cl_uint), NULL);
   ezcl_set_kernel_arg(kernel_rezone_all, 24, local_work_size * sizeof(cl_real4),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_rezone_all,   1, NULL, &global_work_size, &local_work_size, &rezone_all_event);

   if (ncells != old_ncells){
      ezcl_device_memory_remove(dev_H);
      ezcl_device_memory_remove(dev_U);
      ezcl_device_memory_remove(dev_V);
      dev_H = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_U = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_V = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

      mesh->resize_old_device_memory(ncells);
   }

   SWAP_PTR(dev_H_new, dev_H, dev_ptr);
   SWAP_PTR(dev_U_new, dev_U, dev_ptr);
   SWAP_PTR(dev_V_new, dev_V, dev_ptr);

   mesh->swap_device_memory_ptrs();

   ezcl_device_memory_remove(dev_mpot);
   ezcl_device_memory_remove(dev_ijadd);

   ezcl_device_memory_remove(dev_H_new);
   ezcl_device_memory_remove(dev_U_new);
   ezcl_device_memory_remove(dev_V_new);

   gpu_time_rezone_all        += ezcl_timer_calc(&rezone_all_event,        &rezone_all_event);
}

void State::gpu_rezone_all_local(cl_command_queue command_queue, Mesh *mesh, size_t &ncells, size_t new_ncells, size_t old_ncells, bool localStencil,
      cl_mem dev_mpot, cl_mem dev_ioffset)
{
   cl_event rezone_all_event;

   cl_mem dev_ptr = NULL;

   cl_mem &dev_level        = mesh->dev_level;
   cl_mem &dev_i            = mesh->dev_i;
   cl_mem &dev_j            = mesh->dev_j;
   cl_mem &dev_celltype     = mesh->dev_celltype;
   cl_mem &dev_level_new    = mesh->dev_level_new;
   cl_mem &dev_i_new        = mesh->dev_i_new;
   cl_mem &dev_j_new        = mesh->dev_j_new;
   cl_mem &dev_celltype_new = mesh->dev_celltype_new;
   cl_mem &dev_levdx        = mesh->dev_levdx;
   cl_mem &dev_levdy        = mesh->dev_levdy;

   assert(dev_mpot);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);
   assert(dev_celltype);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_level_new);
   assert(dev_i_new);
   assert(dev_j_new);
   assert(dev_celltype_new);
   assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   cl_event start_read_event, end_read_event, start_write_event;

   vector<int> i_tmp(ncells);
   vector<int> j_tmp(ncells);
   vector<int> level_tmp(ncells);
   ezcl_enqueue_read_buffer(command_queue,  dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int), &i_tmp[0],     &start_read_event);
   ezcl_enqueue_read_buffer(command_queue,  dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int), &j_tmp[0],     NULL);
   ezcl_enqueue_read_buffer(command_queue,  dev_level, CL_TRUE,  0, ncells*sizeof(cl_int), &level_tmp[0], &end_read_event);

   int ifirst      = 0;
   int ilast       = 0;
   int jfirst      = 0;
   int jlast       = 0;
   int level_first = 0;
   int level_last  = 0;

   MPI_Request req[12];
   MPI_Status status[12];

   static int prev     = MPI_PROC_NULL;
   static int next     = MPI_PROC_NULL;

   if (mesh->mype != 0)               prev = mesh->mype-1;
   if (mesh->mype < mesh->numpe - 1)  next = mesh->mype+1;

   MPI_Isend(&i_tmp[ncells-1],     1,MPI_INT,next,1,MPI_COMM_WORLD,req+0);
   MPI_Irecv(&ifirst,              1,MPI_INT,prev,1,MPI_COMM_WORLD,req+1);

   MPI_Isend(&i_tmp[0],            1,MPI_INT,prev,1,MPI_COMM_WORLD,req+2);
   MPI_Irecv(&ilast,               1,MPI_INT,next,1,MPI_COMM_WORLD,req+3);

   MPI_Isend(&j_tmp[ncells-1],     1,MPI_INT,next,1,MPI_COMM_WORLD,req+4);
   MPI_Irecv(&jfirst,              1,MPI_INT,prev,1,MPI_COMM_WORLD,req+5);

   MPI_Isend(&j_tmp[0],            1,MPI_INT,prev,1,MPI_COMM_WORLD,req+6);
   MPI_Irecv(&jlast,               1,MPI_INT,next,1,MPI_COMM_WORLD,req+7);

   MPI_Isend(&level_tmp[ncells-1], 1,MPI_INT,next,1,MPI_COMM_WORLD,req+8);
   MPI_Irecv(&level_first,         1,MPI_INT,prev,1,MPI_COMM_WORLD,req+9);

   MPI_Isend(&level_tmp[0],        1,MPI_INT,prev,1,MPI_COMM_WORLD,req+10);
   MPI_Irecv(&level_last,          1,MPI_INT,next,1,MPI_COMM_WORLD,req+11);

   MPI_Waitall(12, req, status);

   cl_mem dev_ijadd;

   vector<int>ijadd(6);
   ijadd[0] = ifirst;
   ijadd[1] = ilast;
   ijadd[2] = jfirst;
   ijadd[3] = jlast;
   ijadd[4] = level_first;
   ijadd[5] = level_last;

   size_t six = 6;
   dev_ijadd = ezcl_malloc(NULL, &six, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
   ezcl_enqueue_write_buffer(command_queue, dev_ijadd, CL_TRUE, 0, 6*sizeof(cl_int), (void*)&ijadd[0], &start_write_event);

   if (new_ncells != old_ncells){
      ncells = new_ncells;
      mesh->resize_new_device_memory(ncells);
   }

   cl_mem dev_H_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_U_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_V_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

   int stencil = 0;
   if (localStencil) stencil = 1;

   size_t local_work_size = 128;
   size_t global_work_size = ((old_ncells+local_work_size - 1) /local_work_size) * local_work_size;

   ezcl_set_kernel_arg(kernel_rezone_all, 0,  sizeof(cl_int),  (void *)&old_ncells);
   ezcl_set_kernel_arg(kernel_rezone_all, 1,  sizeof(cl_int),  (void *)&stencil);
   ezcl_set_kernel_arg(kernel_rezone_all, 2,  sizeof(cl_int),  (void *)&mesh->levmx);
   ezcl_set_kernel_arg(kernel_rezone_all, 3,  sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_rezone_all, 4,  sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_rezone_all, 5,  sizeof(cl_mem),  (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_rezone_all, 6,  sizeof(cl_mem),  (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_rezone_all, 7,  sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_rezone_all, 8,  sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_rezone_all, 9,  sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_rezone_all, 10, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_rezone_all, 11, sizeof(cl_mem),  (void *)&dev_level_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 12, sizeof(cl_mem),  (void *)&dev_i_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 13, sizeof(cl_mem),  (void *)&dev_j_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 14, sizeof(cl_mem),  (void *)&dev_celltype_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 15, sizeof(cl_mem),  (void *)&dev_H_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 16, sizeof(cl_mem),  (void *)&dev_U_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 17, sizeof(cl_mem),  (void *)&dev_V_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 18, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_rezone_all, 19, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_rezone_all, 20, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_rezone_all, 21, sizeof(cl_mem),  (void *)&mesh->dev_levtable);
   ezcl_set_kernel_arg(kernel_rezone_all, 22, sizeof(cl_mem),  (void *)&dev_ijadd);
   ezcl_set_kernel_arg(kernel_rezone_all, 23, local_work_size * sizeof(cl_int), NULL);
   ezcl_set_kernel_arg(kernel_rezone_all, 24, local_work_size * sizeof(cl_real4),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_rezone_all,   1, NULL, &global_work_size, &local_work_size, &rezone_all_event);

   if (ncells != old_ncells){
      ezcl_device_memory_remove(dev_H);
      ezcl_device_memory_remove(dev_U);
      ezcl_device_memory_remove(dev_V);
      dev_H = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_U = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_V = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

      mesh->resize_old_device_memory(ncells);
   }

   SWAP_PTR(dev_H_new, dev_H, dev_ptr);
   SWAP_PTR(dev_U_new, dev_U, dev_ptr);
   SWAP_PTR(dev_V_new, dev_V, dev_ptr);

   mesh->swap_device_memory_ptrs();

   ezcl_device_memory_remove(dev_mpot);
   ezcl_device_memory_remove(dev_ijadd);

   ezcl_device_memory_remove(dev_H_new);
   ezcl_device_memory_remove(dev_U_new);
   ezcl_device_memory_remove(dev_V_new);

   gpu_time_rezone_all        += ezcl_timer_calc(&rezone_all_event,        &rezone_all_event);
}

//define macro for squaring a number
#define SQ(x) ((x)*(x))
//define macro to find minimum of 3 values
//#define MIN3(a,b,c) (min(min((a),(b)),(c)))

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

// XXX ADDED XXX
#define HXFLUXNLT ( Ult )
#define HXFLUXNRT ( Urt )
#define UXFLUXNLT ( SQR(Ult)/Hlt + ghalf*SQR(Hlt) )
#define UXFLUXNRT ( SQR(Urt)/Hrt + ghalf*SQR(Hrt) )
#define UVFLUXNLT ( Ult*Vlt/Hlt )
#define UVFLUXNRT ( Urt*Vrt/Hrt )
#define HYFLUXNBR ( Vbr )
#define HYFLUXNTR ( Vtr )
#define VUFLUXNBR  ( Vbr*Ubr/Hbr )
#define VUFLUXNTR  ( Vtr*Utr/Htr )
#define VYFLUXNBR  ( SQR(Vbr)/Hbr + ghalf*SQR(Hbr) )
#define VYFLUXNTR  ( SQR(Vtr)/Htr + ghalf*SQR(Htr) )
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

void State::calc_finite_difference(Mesh *mesh, double deltaT){
   double   g     = 9.80;   // gravitational constant
   double   ghalf = 0.5*g;

   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   size_t &ncells     = mesh->ncells;
   vector<int> &nlft  = mesh->nlft;
   vector<int> &nrht  = mesh->nrht;
   vector<int> &nbot  = mesh->nbot;
   vector<int> &ntop  = mesh->ntop;
   vector<int> &level = mesh->level;

   vector<real> &lev_deltax = mesh->lev_deltax;
   vector<real> &lev_deltay = mesh->lev_deltay;

   int nl, nr, nb, nt;
   int nll, nrr, nbb, ntt;
   int nlt, nrt, nbr, ntr;
   double Hxminus, Hxplus;
   double Uxminus, Uxplus;
   double Vxminus, Vxplus;

   double Hyminus, Hyplus;
   double Uyminus, Uyplus;
   double Vyminus, Vyplus;

   vector<real> H_new(ncells);
   vector<real> U_new(ncells);
   vector<real> V_new(ncells);

   double dxic, dxl, dxr, dyic, dyb, dyt;
   double Hic, Hl, Hr, Hb, Ht;
   double Hll, Hrr, Hbb, Htt;
   double Uic, Ul, Ur, Ub, Ut;
   double Ull, Urr, Ubb, Utt;
   double Vic, Vl, Vr, Vb, Vt;
   double Vll, Vrr, Vbb, Vtt;

   double Hlt=0.0, Hrt=0.0, Htr=0.0, Hbr=0.0;
   double Ult=0.0, Urt=0.0, Utr=0.0, Ubr=0.0;
   double Vlt=0.0, Vrt=0.0, Vtr=0.0, Vbr=0.0;

   double Hxminus2=0.0, Hxplus2=0.0;
   double Uxminus2=0.0, Uxplus2=0.0;
   double Vxminus2=0.0, Vxplus2=0.0;

   double Hyminus2=0.0, Hyplus2=0.0;
   double Uyminus2, Uyplus2;
   double Vyminus2=0.0, Vyplus2=0.0;

   double dric, drl, drr, drt, drb;

   double Hxfluxminus;
   double Uxfluxminus;
   double Vxfluxminus;

   double Hxfluxplus;
   double Uxfluxplus;
   double Vxfluxplus;

   double Hyfluxminus;
   double Uyfluxminus;
   double Vyfluxminus;

   double Hyfluxplus;
   double Uyfluxplus;
   double Vyfluxplus;

   int lvl;

   real wminusx_H, wminusx_U;
   real wplusx_H, wplusx_U;
   real wminusy_H, wminusy_V;
   real wplusy_H, wplusy_V;

   int nltl=0;
   real Hll2=ZERO;

   int nrtr=0;
   real Hrr2=ZERO;

   real Ull2=ZERO;
   real Urr2=ZERO;

   int nbrb=0;
   real Hbb2=ZERO;

   int ntrt=0;
   real Htt2=ZERO;

   real Vbb2=ZERO;
   real Vtt2=ZERO;

   for(uint gix = 0; gix < ncells; gix++) {

      lvl     = level[gix];
      nl      = nlft[gix];
      nr      = nrht[gix];
      nt      = ntop[gix];
      nb      = nbot[gix];

      Hic     = H[gix];
      Uic     = U[gix];
      Vic     = V[gix];

      nll     = nlft[nl];
      Hl      = H[nl];
      Ul      = U[nl];
      Vl      = V[nl];

      nrr     = nrht[nr];
      Hr      = H[nr];
      Ur      = U[nr];
      Vr      = V[nr];

      ntt     = ntop[nt];
      Ht      = H[nt];
      Ut      = U[nt];
      Vt      = V[nt];

      nbb     = nbot[nb];
      Hb      = H[nb];
      Ub      = U[nb];
      Vb      = V[nb];

      nlt     = ntop[nl];
      nrt     = ntop[nr];
      ntr     = nrht[nt];
      nbr     = nrht[nb];

      Hll     = H[nll];
      Ull     = U[nll];
      Vll     = V[nll];

      Hrr     = H[nrr];
      Urr     = U[nrr];
      Vrr     = V[nrr];

      Htt     = H[ntt];
      Utt     = U[ntt];
      Vtt     = V[ntt];

      Hbb     = H[nbb];
      Ubb     = U[nbb];
      Vbb     = V[nbb];

      dxic    = lev_deltax[lvl];
      dyic    = lev_deltay[lvl];
 
      dxl     = lev_deltax[level[nl]];
      dxr     = lev_deltax[level[nr]];
 
      dyt     = lev_deltay[level[nt]];
      dyb     = lev_deltay[level[nb]];
 
      drl     = dxl;
      drr     = dxr;
      drt     = dyt;
      drb     = dyb;

      dric    = dxic;

      if(lvl < level[nl]) {
         Hlt     = H[ ntop[nl] ];
         Ult     = U[ ntop[nl] ]; 
         Vlt     = V[ ntop[nl] ]; 
         nltl    = nlft[nlt];
         Hll2 = H[nltl];
         Ull2 = U[nltl];
      }
    
      if(lvl < level[nr]) {
         Hrt     = H[ ntop[nr] ];
         Urt     = U[ ntop[nr] ];
         Vrt     = V[ ntop[nr] ];
         nrtr    = nrht[nrt];
         Hrr2 = H[nrtr];
         Urr2 = U[nrtr];
      }

      if(lvl < level[nb]) {
         Hbr     = H[ nrht[nb] ];
         Ubr     = U[ nrht[nb] ];
         Vbr     = V[ nrht[nb] ];
         nbrb    = nbot[nbr];
         Hbb2 = H[nbrb];
         Vbb2 = V[nbrb];
      }

      if(lvl < level[nt]) {
         Htr     = H[ nrht[nt] ];
         Utr     = U[ nrht[nt] ];
         Vtr     = V[ nrht[nt] ];
         ntrt    = ntop[ntr];
         Htt2 = H[ntrt];
         Vtt2 = V[ntrt];
      }


      Hxminus = U_halfstep(deltaT, Hl, Hic, HXFLUXNL, HXFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
      Uxminus = U_halfstep(deltaT, Ul, Uic, UXFLUXNL, UXFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
      Vxminus = U_halfstep(deltaT, Vl, Vic, UVFLUXNL, UVFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));

      Hxplus  = U_halfstep(deltaT, Hic, Hr, HXFLUXIC, HXFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
      Uxplus  = U_halfstep(deltaT, Uic, Ur, UXFLUXIC, UXFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
      Vxplus  = U_halfstep(deltaT, Vic, Vr, UVFLUXIC, UVFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));

      Hyminus = U_halfstep(deltaT, Hb, Hic, HYFLUXNB, HYFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
      Uyminus = U_halfstep(deltaT, Ub, Uic, VUFLUXNB, VUFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
      Vyminus = U_halfstep(deltaT, Vb, Vic, VYFLUXNB, VYFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));

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


      if(lvl < level[nl]) {

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

      if(lvl < level[nr]) {

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

      if(lvl < level[nb]) {

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

      if(lvl < level[nt]) {

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



   ////////////////////////////////////////
   /// Artificial Viscosity corrections ///
   ////////////////////////////////////////


   if(level[nl] < level[nll]) {
      Hll = (Hll + H[ ntop[nll] ]) * HALF;
      Ull = (Ull + U[ ntop[nll] ]) * HALF;
   }

   real Hr2 = Hr;
   real Ur2 = Ur;
   if(lvl < level[nr]) {
      Hr2 = (Hr2 + Hrt) * HALF;
      Ur2 = (Ur2 + Urt) * HALF;
   }

   wminusx_H = w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                           Hic-Hl, Hl-Hll, Hr2-Hic);

   wminusx_H *= Hic - Hl;

   if(lvl < level[nl]) {
      if(level[nlt] < level[nltl])
         Hll2 = (Hll2 + H[ ntop[nltl] ]) * HALF;
      wminusx_H = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                               sqrt(g*Hxminus2), Hic-Hlt, Hlt-Hll2, Hr2-Hic) *
                   (Hic - Hlt)) + wminusx_H)*HALF*HALF;
   }


   if(level[nr] < level[nrr]) {
      Hrr = (Hrr + H[ ntop[nrr] ]) * HALF;
      Urr = (Urr + U[ ntop[nrr] ]) * HALF;
   }

   real Hl2 = Hl;
   real Ul2 = Ul;
   if(lvl < level[nl]) {
      Hl2 = (Hl2 + Hlt) * HALF;
      Ul2 = (Ul2 + Ult) * HALF;
   }

   wplusx_H = w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                           Hr-Hic, Hic-Hl2, Hrr-Hr);

   wplusx_H *= Hr - Hic;

   if(lvl < level[nr]) {
      if(level[nrt] < level[nrtr])
         Hrr2 = (Hrr2 + H[ ntop[nrtr] ]) * HALF;
      wplusx_H = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                               sqrt(g*Hxplus2), Hrt-Hic, Hic-Hl2, Hrr2-Hrt) *
                   (Hrt - Hic))+wplusx_H)*HALF*HALF;
   }


   wminusx_U = w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                           Uic-Ul, Ul-Ull, Ur2-Uic);

   wminusx_U *= Uic - Ul;

   if(lvl < level[nl]) {
      if(level[nlt] < level[nltl])
         Ull2 = (Ull2 + U[ ntop[nltl] ]) * HALF;
      wminusx_U = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                               sqrt(g*Hxminus2), Uic-Ult, Ult-Ull2, Ur2-Uic) *
                   (Uic - Ult))+wminusx_U)*HALF*HALF;
   }


   wplusx_U = w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                           Ur-Uic, Uic-Ul2, Urr-Ur);

   wplusx_U *= Ur - Uic;

   if(lvl < level[nr]) {
      if(level[nrt] < level[nrtr])
         Urr2 = (Urr2 + U[ ntop[nrtr] ]) * HALF;
      wplusx_U = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                               sqrt(g*Hxplus2), Urt-Uic, Uic-Ul2, Urr2-Urt) *
                   (Urt - Uic))+wplusx_U)*HALF*HALF;
   }


   if(level[nb] < level[nbb]) {
      Hbb = (Hbb + H[ nrht[nbb] ]) * HALF;
      Vbb = (Vbb + V[ nrht[nbb] ]) * HALF;
   }

   real Ht2 = Ht;
   real Vt2 = Vt;
   if(lvl < level[nt]) {
      Ht2 = (Ht2 + Htr) * HALF;
      Vt2 = (Vt2 + Vtr) * HALF;
   }

   wminusy_H = w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                           Hic-Hb, Hb-Hbb, Ht2-Hic);

   wminusy_H *= Hic - Hb;

   if(lvl < level[nb]) {
      if(level[nbr] < level[nbrb])
         Hbb2 = (Hbb2 + H[ nrht[nbrb] ]) * HALF;
      wminusy_H = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                               sqrt(g*Hyminus2), Hic-Hbr, Hbr-Hbb2, Ht2-Hic) *
                   (Hic - Hbr))+wminusy_H)*HALF*HALF;
   }


   if(level[nt] < level[ntt]) {
      Htt = (Htt + H[ nrht[ntt] ]) * HALF;
      Vtt = (Vtt + V[ nrht[ntt] ]) * HALF;
   }

   real Hb2 = Hb;
   real Vb2 = Vb;
   if(lvl < level[nb]) {
      Hb2 = (Hb2 + Hbr) * HALF;
      Vb2 = (Vb2 + Vbr) * HALF;
   }

   wplusy_H = w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                           Ht-Hic, Hic-Hb2, Htt-Ht);

   wplusy_H *= Ht - Hic;

   if(lvl < level[nt]) {
      if(level[ntr] < level[ntrt])
         Htt2 = (Htt2 + H[ nrht[ntrt] ]) * HALF;
      wplusy_H = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                               sqrt(g*Hyplus2), Htr-Hic, Hic-Hb2, Htt2-Htr) *
                   (Htr - Hic))+wplusy_H)*HALF*HALF;
   }

   wminusy_V = w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                           Vic-Vb, Vb-Vbb, Vt2-Vic);

   wminusy_V *= Vic - Vb;

   if(lvl < level[nb]) {
      if(level[nbr] < level[nbrb])
         Vbb2 = (Vbb2 + V[ nrht[nbrb] ]) * HALF;
      wminusy_V = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                               sqrt(g*Hyminus2), Vic-Vbr, Vbr-Vbb2, Vt2-Vic) *
                   (Vic - Vbr))+wminusy_V)*HALF*HALF;
   }

   wplusy_V = w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                           Vt-Vic, Vic-Vb2, Vtt-Vt);

   wplusy_V *= Vt - Vic;

   if(lvl < level[nt]) {
      if(level[ntr] < level[ntrt])
         Vtt2 = (Vtt2 + V[ nrht[ntrt] ]) * HALF;
      wplusy_V = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                               sqrt(g*Hyplus2), Vtr-Vic, Vic-Vb2, Vtt2-Vtr) *
                   (Vtr - Vic))+wplusy_V)*HALF*HALF;
   }

      H_new[gix] = U_fullstep(deltaT, dxic, Hic,
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                  - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
      U_new[gix] = U_fullstep(deltaT, dxic, Uic,
                       Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                  - wminusx_U + wplusx_U;
      V_new[gix] = U_fullstep(deltaT, dxic, Vic,
                       Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                  - wminusy_V + wplusy_V;

   }

   H.swap(H_new);
   U.swap(U_new);
   V.swap(V_new);

   H_new.clear();
   U_new.clear();
   V_new.clear();

   cpu_time_finite_difference += cpu_timer_stop(tstart_cpu);
}

void State::calc_finite_difference_local(Mesh *mesh, double deltaT){
   double   g     = 9.80;   // gravitational constant
   double   ghalf = 0.5*g;

   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   size_t &ncells     = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   vector<int> &nlft  = mesh->nlft;
   vector<int> &nrht  = mesh->nrht;
   vector<int> &nbot  = mesh->nbot;
   vector<int> &ntop  = mesh->ntop;
   vector<int> &level = mesh->level;

   vector<real> &lev_deltax = mesh->lev_deltax;
   vector<real> &lev_deltay = mesh->lev_deltay;

   int nl, nr, nb, nt;
   int nll, nrr, nbb, ntt;
   int nlt, nrt, nbr, ntr;
   double Hxminus, Hxplus;
   double Uxminus, Uxplus;
   double Vxminus, Vxplus;

   double Hyminus, Hyplus;
   double Uyminus, Uyplus;
   double Vyminus, Vyplus;

   vector<real> Hnew(ncells_ghost);
   vector<real> Unew(ncells_ghost);
   vector<real> Vnew(ncells_ghost);

   double dxic, dxl, dxr, dyic, dyb, dyt;
   double Hic, Hl, Hr, Hb, Ht;
   double Hll, Hrr, Hbb, Htt;
   double Uic, Ul, Ur, Ub, Ut;
   double Ull, Urr, Ubb, Utt;
   double Vic, Vl, Vr, Vb, Vt;
   double Vll, Vrr, Vbb, Vtt;

   double Hlt=0.0, Hrt=0.0, Htr=0.0, Hbr=0.0;
   double Ult=0.0, Urt=0.0, Utr=0.0, Ubr=0.0;
   double Vlt=0.0, Vrt=0.0, Vtr=0.0, Vbr=0.0;

   double Hxminus2=0.0, Hxplus2=0.0;
   double Uxminus2=0.0, Uxplus2=0.0;
   double Vxminus2, Vxplus2;

   double Hyminus2=0.0, Hyplus2=0.0;
   double Uyminus2, Uyplus2;
   double Vyminus2=0.0, Vyplus2=0.0;

   double dric, drl, drr, drt, drb;

   double Hxfluxminus;
   double Uxfluxminus;
   double Vxfluxminus;

   double Hxfluxplus;
   double Uxfluxplus;
   double Vxfluxplus;

   double Hyfluxminus;
   double Uyfluxminus;
   double Vyfluxminus;

   double Hyfluxplus;
   double Uyfluxplus;
   double Vyfluxplus;

   int lvl;

   real wminusx_H, wminusx_U;
   real wplusx_H, wplusx_U;
   real wminusy_H, wminusy_V;
   real wplusy_H, wplusy_V;

   int nltl=0;
   real Hll2=ZERO;

   int nrtr=0;
   real Hrr2=ZERO;

   real Ull2=ZERO;
   real Urr2=ZERO;

   int nbrb=0;
   real Hbb2=ZERO;

   int ntrt=0;
   real Htt2=ZERO;

   real Vbb2=ZERO;
   real Vtt2=ZERO;

   for(uint gix = 0; gix < ncells; gix++) {

      lvl     = level[gix];
      nl      = nlft[gix];
      nr      = nrht[gix];
      nt      = ntop[gix];
      nb      = nbot[gix];

      Hic     = H[gix];
      Uic     = U[gix];
      Vic     = V[gix];

      nll     = nlft[nl];
      Hl      = H[nl];
      Ul      = U[nl];
      Vl      = V[nl];

      nrr     = nrht[nr];
      Hr      = H[nr];
      Ur      = U[nr];
      Vr      = V[nr];

      ntt     = ntop[nt];
      Ht      = H[nt];
      Ut      = U[nt];
      Vt      = V[nt];

      nbb     = nbot[nb];
      Hb      = H[nb];
      Ub      = U[nb];
      Vb      = V[nb];

      nlt     = ntop[nl];
      nrt     = ntop[nr];
      ntr     = nrht[nt];
      nbr     = nrht[nb];

      Hll     = H[nll];
      Ull     = U[nll];
      Vll     = V[nll];

      Hrr     = H[nrr];
      Urr     = U[nrr];
      Vrr     = V[nrr];

      Htt     = H[ntt];
      Utt     = U[ntt];
      Vtt     = V[ntt];

      Hbb     = H[nbb];
      Ubb     = U[nbb];
      Vbb     = V[nbb];

      dxic    = lev_deltax[lvl];
      dyic    = lev_deltay[lvl];

      dxl     = lev_deltax[level[nl]];
      dxr     = lev_deltax[level[nr]];

      dyt     = lev_deltay[level[nt]];
      dyb     = lev_deltay[level[nb]];

      drl     = dxl;
      drr     = dxr;
      drt     = dyt;
      drb     = dyb;

      dric    = dxic;

      if(lvl < level[nl]) {
         Hlt  = H[ ntop[nl] ];
         Ult  = U[ ntop[nl] ];
         Vlt  = V[ ntop[nl] ];
         nltl = nlft[nlt];
         Hll2 = H[nltl];
         Ull2 = U[nltl];
      }

      if(lvl < level[nr]) {
         Hrt  = H[ ntop[nr] ];
         Urt  = U[ ntop[nr] ];
         Vrt  = V[ ntop[nr] ];
         nrtr = nrht[nrt];
         Hrr2 = H[nrtr];
         Urr2 = U[nrtr];
      }


      if(lvl < level[nb]) {
         Hbr  = H[ nrht[nb] ];
         Ubr  = U[ nrht[nb] ];
         Vbr  = V[ nrht[nb] ];
         nbrb = nbot[nbr];
         Hbb2 = H[nbrb];
         Vbb2 = V[nbrb];
      }

      if(lvl < level[nt]) {
         Htr  = H[ nrht[nt] ];
         Utr  = U[ nrht[nt] ];
         Vtr  = V[ nrht[nt] ];
         ntrt = ntop[ntr];
         Htt2 = H[ntrt];
         Vtt2 = V[ntrt];
      }


      Hxminus = U_halfstep(deltaT, Hl, Hic, HXFLUXNL, HXFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
      Uxminus = U_halfstep(deltaT, Ul, Uic, UXFLUXNL, UXFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
      Vxminus = U_halfstep(deltaT, Vl, Vic, UVFLUXNL, UVFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));

      Hxplus  = U_halfstep(deltaT, Hic, Hr, HXFLUXIC, HXFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
      Uxplus  = U_halfstep(deltaT, Uic, Ur, UXFLUXIC, UXFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
      Vxplus  = U_halfstep(deltaT, Vic, Vr, UVFLUXIC, UVFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));

      Hyminus = U_halfstep(deltaT, Hb, Hic, HYFLUXNB, HYFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
      Uyminus = U_halfstep(deltaT, Ub, Uic, VUFLUXNB, VUFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
      Vyminus = U_halfstep(deltaT, Vb, Vic, VYFLUXNB, VYFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));

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


      if(lvl < level[nl]) {

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

      if(lvl < level[nr]) {

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

      if(lvl < level[nb]) {

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

      if(lvl < level[nt]) {

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



   ////////////////////////////////////////
   /// Artificial Viscosity corrections ///
   ////////////////////////////////////////


   if(level[nl] < level[nll]) {
      Hll = (Hll + H[ ntop[nll] ]) * HALF;
      Ull = (Ull + U[ ntop[nll] ]) * HALF;
   }

   real Hr2 = Hr;
   real Ur2 = Ur;
   if(lvl < level[nr]) {
      Hr2 = (Hr2 + Hrt) * HALF;
      Ur2 = (Ur2 + Urt) * HALF;
   }

   wminusx_H = w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                           Hic-Hl, Hl-Hll, Hr2-Hic);

   wminusx_H *= Hic - Hl;

   if(lvl < level[nl]) {
      if(level[nlt] < level[nltl])
         Hll2 = (Hll2 + H[ ntop[nltl] ]) * HALF;
      wminusx_H = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                               sqrt(g*Hxminus2), Hic-Hlt, Hlt-Hll2, Hr2-Hic) *
                   (Hic - Hlt)) + wminusx_H)*HALF*HALF;
   }


   if(level[nr] < level[nrr]) {
      Hrr = (Hrr + H[ ntop[nrr] ]) * HALF;
      Urr = (Urr + U[ ntop[nrr] ]) * HALF;
   }

   real Hl2 = Hl;
   real Ul2 = Ul;
   if(lvl < level[nl]) {
      Hl2 = (Hl2 + Hlt) * HALF;
      Ul2 = (Ul2 + Ult) * HALF;
   }

   wplusx_H = w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                           Hr-Hic, Hic-Hl2, Hrr-Hr);

   wplusx_H *= Hr - Hic;

   if(lvl < level[nr]) {
      if(level[nrt] < level[nrtr])
         Hrr2 = (Hrr2 + H[ ntop[nrtr] ]) * HALF;
      wplusx_H = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                               sqrt(g*Hxplus2), Hrt-Hic, Hic-Hl2, Hrr2-Hrt) *
                   (Hrt - Hic))+wplusx_H)*HALF*HALF;
   }


   wminusx_U = w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                           Uic-Ul, Ul-Ull, Ur2-Uic);

   wminusx_U *= Uic - Ul;

   if(lvl < level[nl]) {
      if(level[nlt] < level[nltl])
         Ull2 = (Ull2 + U[ ntop[nltl] ]) * HALF;
      wminusx_U = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                               sqrt(g*Hxminus2), Uic-Ult, Ult-Ull2, Ur2-Uic) *
                   (Uic - Ult))+wminusx_U)*HALF*HALF;
   }


   wplusx_U = w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                           Ur-Uic, Uic-Ul2, Urr-Ur);

   wplusx_U *= Ur - Uic;

   if(lvl < level[nr]) {
      if(level[nrt] < level[nrtr])
         Urr2 = (Urr2 + U[ ntop[nrtr] ]) * HALF;
      wplusx_U = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                               sqrt(g*Hxplus2), Urt-Uic, Uic-Ul2, Urr2-Urt) *
                   (Urt - Uic))+wplusx_U)*HALF*HALF;
   }


   if(level[nb] < level[nbb]) {
      Hbb = (Hbb + H[ nrht[nbb] ]) * HALF;
      Vbb = (Vbb + V[ nrht[nbb] ]) * HALF;
   }

   real Ht2 = Ht;
   real Vt2 = Vt;
   if(lvl < level[nt]) {
      Ht2 = (Ht2 + Htr) * HALF;
      Vt2 = (Vt2 + Vtr) * HALF;
   }

   wminusy_H = w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                           Hic-Hb, Hb-Hbb, Ht2-Hic);

   wminusy_H *= Hic - Hb;

   if(lvl < level[nb]) {
      if(level[nbr] < level[nbrb])
         Hbb2 = (Hbb2 + H[ nrht[nbrb] ]) * HALF;
      wminusy_H = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                               sqrt(g*Hyminus2), Hic-Hbr, Hbr-Hbb2, Ht2-Hic) *
                   (Hic - Hbr))+wminusy_H)*HALF*HALF;
   }


   if(level[nt] < level[ntt]) {
      Htt = (Htt + H[ nrht[ntt] ]) * HALF;
      Vtt = (Vtt + V[ nrht[ntt] ]) * HALF;
   }

   real Hb2 = Hb;
   real Vb2 = Vb;
   if(lvl < level[nb]) {
      Hb2 = (Hb2 + Hbr) * HALF;
      Vb2 = (Vb2 + Vbr) * HALF;
   }

   wplusy_H = w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                           Ht-Hic, Hic-Hb2, Htt-Ht);

   wplusy_H *= Ht - Hic;

   if(lvl < level[nt]) {
      if(level[ntr] < level[ntrt])
         Htt2 = (Htt2 + H[ nrht[ntrt] ]) * HALF;
      wplusy_H = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                               sqrt(g*Hyplus2), Htr-Hic, Hic-Hb2, Htt2-Htr) *
                   (Htr - Hic))+wplusy_H)*HALF*HALF;
   }

   wminusy_V = w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                           Vic-Vb, Vb-Vbb, Vt2-Vic);

   wminusy_V *= Vic - Vb;

   if(lvl < level[nb]) {
      if(level[nbr] < level[nbrb])
         Vbb2 = (Vbb2 + V[ nrht[nbrb] ]) * HALF;
      wminusy_V = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                               sqrt(g*Hyminus2), Vic-Vbr, Vbr-Vbb2, Vt2-Vic) *
                   (Vic - Vbr))+wminusy_V)*HALF*HALF;
   }

   wplusy_V = w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                           Vt-Vic, Vic-Vb2, Vtt-Vt);

   wplusy_V *= Vt - Vic;

   if(lvl < level[nt]) {
      if(level[ntr] < level[ntrt])
         Vtt2 = (Vtt2 + V[ nrht[ntrt] ]) * HALF;
      wplusy_V = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                               sqrt(g*Hyplus2), Vtr-Vic, Vic-Vb2, Vtt2-Vtr) *
                   (Vtr - Vic))+wplusy_V)*HALF*HALF;
   }

      Hnew[gix] = U_fullstep(deltaT, dxic, Hic,
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                  - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
      Unew[gix] = U_fullstep(deltaT, dxic, Uic,
                       Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                  - wminusx_U + wplusx_U;
      Vnew[gix] = U_fullstep(deltaT, dxic, Vic,
                       Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                  - wminusy_V + wplusy_V;

   }

   H.swap(Hnew);
   U.swap(Unew);
   V.swap(Vnew);

   L7_Update(&H[0], L7_REAL, mesh->cell_handle);
   L7_Update(&U[0], L7_REAL, mesh->cell_handle);
   L7_Update(&V[0], L7_REAL, mesh->cell_handle);

   Hnew.clear();
   Unew.clear();
   Vnew.clear();

   cpu_time_finite_difference += cpu_timer_stop(tstart_cpu);
}

void State::gpu_calc_finite_difference(cl_command_queue command_queue, Mesh *mesh, double deltaT)
{
   cl_event finite_difference_event;

   cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   int &levmx        = mesh->levmx;
   cl_mem &dev_nlft  = mesh->dev_nlft;
   cl_mem &dev_nrht  = mesh->dev_nrht;
   cl_mem &dev_nbot  = mesh->dev_nbot;
   cl_mem &dev_ntop  = mesh->dev_ntop;
   cl_mem &dev_level = mesh->dev_level;
   cl_mem &dev_levdx = mesh->dev_levdx;
   cl_mem &dev_levdy = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

   cl_mem dev_H_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_U_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_V_new = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

     /*
     __kernel void calc_finite_difference_cl(
                      const int    ncells,   // 0  Total number of cells.
                      const int    lvmax,    // 1  Maximum level
             __global       real  *H,        // 2
             __global       real  *U,        // 3
             __global       real  *V,        // 4
             __global       real  *H_new,    // 5
             __global       real  *U_new,    // 6
             __global       real  *V_new,    // 7
             __global const int   *nlft,     // 8  Array of left neighbors.
             __global const int   *nrht,     // 9  Array of right neighbors.
             __global const int   *ntop,     // 10  Array of bottom neighbors.
             __global const int   *nbot,     // 11  Array of top neighbors.
             __global const int   *level,    // 12  Array of level information.
                      const real   deltaT,   // 13  Size of time step.
             __global const real  *lev_dx,   // 14
             __global const real  *lev_dy,   // 15
             __local        real4 *tile,     // 16  Tile size in real4.
             __local        int8  *itile)    // 17  Tile size in int8.
     */

   real deltaT_local = deltaT;
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 1, sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 5, sizeof(cl_mem),  (void *)&dev_H_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 6, sizeof(cl_mem),  (void *)&dev_U_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 7, sizeof(cl_mem),  (void *)&dev_V_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 8, sizeof(cl_mem),  (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 9, sizeof(cl_mem),  (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,10, sizeof(cl_mem),  (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,11, sizeof(cl_mem),  (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,12, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,13, sizeof(cl_real), (void *)&deltaT_local);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,14, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,15, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,16, local_work_size*sizeof(cl_real4),    NULL);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,17, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference,   1, NULL, &global_work_size, &local_work_size, &finite_difference_event);

   SWAP_PTR(dev_H_new, dev_H, dev_ptr);
   SWAP_PTR(dev_U_new, dev_U, dev_ptr);
   SWAP_PTR(dev_V_new, dev_V, dev_ptr);

   ezcl_device_memory_remove(dev_H_new);
   ezcl_device_memory_remove(dev_U_new);
   ezcl_device_memory_remove(dev_V_new);

   gpu_time_finite_difference += ezcl_timer_calc(&finite_difference_event, &finite_difference_event);
}

void State::gpu_calc_finite_difference_local(cl_command_queue command_queue, Mesh *mesh, double deltaT)
{
   cl_event finite_difference_event;
   cl_event copy_state_data_event;
   cl_event copy_state_ghost_data_event;
   cl_event start_read_event, end_read_event;
   cl_event start_write_event, end_write_event;

   cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   int &levmx        = mesh->levmx;
   cl_mem &dev_nlft  = mesh->dev_nlft;
   cl_mem &dev_nrht  = mesh->dev_nrht;
   cl_mem &dev_nbot  = mesh->dev_nbot;
   cl_mem &dev_ntop  = mesh->dev_ntop;
   cl_mem &dev_level = mesh->dev_level;
   cl_mem &dev_levdx = mesh->dev_levdx;
   cl_mem &dev_levdy = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

   cl_mem dev_H_new = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_U_new = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   cl_mem dev_V_new = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

     /*
     __kernel void copy_state_data_cl(
                      const int  isize,         // 0
             __global      real  *H,            // 1
             __global      real  *U,            // 2
             __global      real  *V,            // 3
             __global      real  *H_new,        // 4
             __global      real  *U_new,        // 5
             __global      real  *V_new)        // 6
     */

   ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
   ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
   ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);

   SWAP_PTR(dev_H_new, dev_H, dev_ptr);
   SWAP_PTR(dev_U_new, dev_U, dev_ptr);
   SWAP_PTR(dev_V_new, dev_V, dev_ptr);

   vector<real> H_tmp(mesh->ncells_ghost,0.0);
   vector<real> U_tmp(mesh->ncells_ghost,0.0);
   vector<real> V_tmp(mesh->ncells_ghost,0.0);
   ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real), &H_tmp[0], &start_read_event);
   ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real), &U_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real), &V_tmp[0], &end_read_event);
   gpu_time_read += ezcl_timer_calc(&start_read_event, &end_read_event);

   L7_Update(&H_tmp[0], L7_REAL, mesh->cell_handle);
   L7_Update(&U_tmp[0], L7_REAL, mesh->cell_handle);
   L7_Update(&V_tmp[0], L7_REAL, mesh->cell_handle);
   //fprintf(mesh->fp,"Line %d\n",__LINE__);
   //for (int ic=0; ic<mesh->ncells_ghost; ic++){
   //   if (mesh->i_local[ic] < 63 || mesh->i_local[ic] > 66 || mesh->j_local[ic] < 63 || mesh->j_local[ic] > 66) continue;
   //   fprintf(mesh->fp,"%d: ic %d H %lf U %lf V %lf\n",mesh->mype,ic,H_tmp[ic],U_tmp[ic],V_tmp[ic]);
   //}

   size_t nghost_local = mesh->ncells_ghost - ncells;
   //fprintf(mesh->fp,"%d: sizes are ncells %d nghost %d ncells_ghost %d\n",mesh->mype,ncells,nghost_local,mesh->ncells_ghost);

   if (mesh->numpe > 1) {
      cl_mem dev_H_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      cl_mem dev_U_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      cl_mem dev_V_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      ezcl_enqueue_write_buffer(command_queue, dev_H_add, CL_FALSE, 0, nghost_local*sizeof(cl_real), (void*)&H_tmp[ncells], &start_write_event);
      ezcl_enqueue_write_buffer(command_queue, dev_U_add, CL_FALSE, 0, nghost_local*sizeof(cl_real), (void*)&U_tmp[ncells], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_V_add, CL_TRUE,  0, nghost_local*sizeof(cl_real), (void*)&V_tmp[ncells], &end_write_event);
      gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

      size_t ghost_local_work_size = 32;
      size_t ghost_global_work_size = ((nghost_local + ghost_local_work_size - 1) /ghost_local_work_size) * ghost_local_work_size;

      // Fill in ghost
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 1, sizeof(cl_int), (void *)&nghost_local);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 2, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 3, sizeof(cl_mem), (void *)&dev_H_add);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 4, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 5, sizeof(cl_mem), (void *)&dev_U_add);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 6, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 7, sizeof(cl_mem), (void *)&dev_V_add);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_ghost_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, &copy_state_ghost_data_event);

      ezcl_device_memory_remove(dev_H_new);
      ezcl_device_memory_remove(dev_U_new);
      ezcl_device_memory_remove(dev_V_new);
      dev_H_new = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_U_new = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_V_new = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   }

     /*
     __kernel void calc_finite_difference_cl(
                      const int    ncells,   // 0  Total number of cells.
                      const int    lvmax,    // 1  Maximum level
             __global       real  *H,        // 2
             __global       real  *U,        // 3
             __global       real  *V,        // 4
             __global       real  *H_new,    // 5
             __global       real  *U_new,    // 6
             __global       real  *V_new,    // 7
             __global const int   *nlft,     // 8  Array of left neighbors.
             __global const int   *nrht,     // 9  Array of right neighbors.
             __global const int   *ntop,     // 10  Array of bottom neighbors.
             __global const int   *nbot,     // 11  Array of top neighbors.
             __global const int   *level,    // 12  Array of level information.
                      const real   deltaT,   // 13  Size of time step.
             __global const real  *lev_dx,   // 14
             __global const real  *lev_dy,   // 15
             __local        real4 *tile,     // 16  Tile size in real4.
             __local        int8  *itile)    // 17  Tile size in int8.
     */

   real deltaT_local = deltaT;
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 1, sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 5, sizeof(cl_mem),  (void *)&dev_H_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 6, sizeof(cl_mem),  (void *)&dev_U_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 7, sizeof(cl_mem),  (void *)&dev_V_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 8, sizeof(cl_mem),  (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 9, sizeof(cl_mem),  (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,10, sizeof(cl_mem),  (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,11, sizeof(cl_mem),  (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,12, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,13, sizeof(cl_real), (void *)&deltaT_local);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,14, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,15, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,16, local_work_size*sizeof(cl_real4),    NULL);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,17, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference,   1, NULL, &global_work_size, &local_work_size, &finite_difference_event);

   SWAP_PTR(dev_H_new, dev_H, dev_ptr);
   SWAP_PTR(dev_U_new, dev_U, dev_ptr);
   SWAP_PTR(dev_V_new, dev_V, dev_ptr);

   ezcl_device_memory_remove(dev_H_new);
   ezcl_device_memory_remove(dev_U_new);
   ezcl_device_memory_remove(dev_V_new);

   gpu_time_finite_difference += ezcl_timer_calc(&copy_state_data_event,       &copy_state_data_event);
   if (mesh->numpe > 1) {
      gpu_time_finite_difference += ezcl_timer_calc(&copy_state_ghost_data_event, &copy_state_ghost_data_event);
   }
   gpu_time_finite_difference += ezcl_timer_calc(&finite_difference_event,     &finite_difference_event);
}

void State::symmetry_check(Mesh *mesh, const char *string, vector<int> sym_index, double eps,
                           SIGN_RULE sign_rule, int &flag)
{
   size_t &ncells = mesh->ncells;

   double xsign = 1.0, ysign = 1.0;

   if (sign_rule == DIAG_RULE || sign_rule == X_RULE) {
      xsign = -1.0;
   }

   if (sign_rule == DIAG_RULE || sign_rule == Y_RULE) {
      ysign = -1.0;
   }

   for (uint ic=0; ic<ncells; ic++) {
      /*  Symmetrical check */
      if (fabs(H[ic] - H[sym_index[ic]]) > eps) {
         printf("%s ic %d sym %d H[ic] %lf Hsym %lf diff %lf\n",
                string,ic,sym_index[ic],H[ic],H[sym_index[ic]],fabs(H[ic]-H[sym_index[ic]]));
         flag++;
      }
      if (fabs(U[ic] - xsign*U[sym_index[ic]]) > eps) {
         printf("%s ic %d sym %d U[ic] %lf Usym %lf diff %lf\n",
                string,ic,sym_index[ic],U[ic],U[sym_index[ic]],fabs(U[ic]-xsign*U[sym_index[ic]]));
         flag++;
      }
      if (fabs(V[ic] - ysign*V[sym_index[ic]]) > eps) {
         printf("%s ic %d sym %d V[ic] %lf Vsym %lf diff %lf\n",
                string,ic,sym_index[ic],V[ic],V[sym_index[ic]],fabs(V[ic]-ysign*V[sym_index[ic]]));
         flag++;
      }
   }

}


void State::calc_refine_potential(Mesh *mesh, vector<int> &mpot,int &icount, int &jcount)
{
   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   size_t &ncells       = mesh->ncells;
   vector<int> &nlft    = mesh->nlft;
   vector<int> &nrht    = mesh->nrht;
   vector<int> &nbot    = mesh->nbot;
   vector<int> &ntop    = mesh->ntop;
   vector<int> &level   = mesh->level;

   vector<double> Q(ncells);

   int nl, nr, nt, nb;
   int nlt, nrt, ntr, nbr;
   double Hic, Hl, Hr, Hb, Ht;
   double Uic, Ul, Ur, Ub, Ut;
   double Vic, Vl, Vr, Vb, Vt;

   double duplus1, duminus1, duhalf1;
   double duplus2, duminus2, duhalf2;

   double qpot, qmax;

   icount=0;
   jcount=0;

   for(uint ic = 0; ic < ncells; ic++) {

      Hic = H[ic];
      Uic = U[ic];
      Vic = V[ic];

      nl = nlft[ic];
      Hl = H[nl];
      Ul = U[nl];
      Vl = V[nl];

      if (level[nl] > level[ic]){
         nlt = ntop[nl];
         Hl = 0.5 * (Hl + H[nlt]);
      }

      nr = nrht[ic];
      Hr = H[nr];
      Ur = U[nr];
      Vr = V[nr];

      if (level[nr] > level[ic]){
         nrt = ntop[nr];
         Hr = 0.5 * (Hr + H[nrt]);
      }

      nb = nbot[ic];
      Hb = H[nb];
      Ub = U[nb];
      Vb = V[nb];

      if (level[nb] > level[ic]){
         nbr = nrht[nb];
         Hb = 0.5 * (Hb + H[nbr]);
      }

      nt = ntop[ic];
      Ht = H[nt];
      Ut = U[nt];
      Vt = V[nt];

      if (level[nt] > level[ic]){
         ntr = nrht[nt];
         Ht = 0.5 * (Ht + H[ntr]);
      }

      duplus1 = Hr-Hic;
      duplus2 = Ur-Uic;
      duhalf1 = Hic-Hl;
      duhalf2 = Uic-Ul;

      qmax = -1000.0;

      qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hl;
      duminus2 = Uic-Ul;
      duhalf1 = Hr-Hic;
      duhalf2 = Ur-Uic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duplus1 = Ht-Hic;
      duplus2 = Vt-Vic;
      duhalf1 = Hic-Hb;
      duhalf2 = Vic-Vb;

      qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hb;
      duminus2 = Vic-Vb;
      duhalf1 = Ht-Hic;
      duhalf2 = Vt-Vic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      Q[ic] = qmax;

   }

   for(uint ic=0; ic<ncells; ic++) {
      mpot[ic]=0;
      if (Q[ic] > REFINE_GRADIENT && level[ic] < mesh->levmx) {
         mpot[ic]=1;
         icount++;
      }
   }

   if(icount > 0) {
      int lev, ll, lr, lt, lb;
      int llt, lrt, ltr, lbr;
      int new_count;
      int levcount = 0;
      do {
         levcount++; 
         new_count=0;
         for(uint ic = 0; ic < ncells; ic++) {
            lev = level[ic];
            if(mpot[ic] > 0) lev++;
   
            nl = nlft[ic];
            ll = level[nl];
            if(mpot[nl] > 0) ll++;
   
            if(ll - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            nlt = ntop[nl];
            llt = level[nlt];
            if(mpot[nlt] > 0) llt++;

            if(llt - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            nr = nrht[ic];
            lr = level[nr];
            if(mpot[nr] > 0) lr++;
   
            if(lr - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            nrt = ntop[nr];
            lrt = level[nrt];
            if(mpot[nrt] > 0) lrt++;

            if(lrt - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            nt = ntop[ic];
            lt = level[nt];
            if(mpot[nt] > 0) lt++;
   
            if(lt - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            ntr = nrht[nt];
            ltr = level[ntr];
            if(mpot[ntr] > 0) ltr++;

            if(ltr - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            nb = nbot[ic];
            lb = level[nb];
            if(mpot[nb] > 0) lb++;
   
            if(lb - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }

            nbr = nrht[nb];
            lbr = level[nbr];
            if(mpot[nbr] > 0) lbr++;

            if(lbr - lev > 1) {
               mpot[ic]++;
               new_count++; icount++;
               continue;
            }
    
         }
//         printf("%d: new_count is %d levmx %d\n",levcount,new_count,mesh->levmx);
      //} while (new_count > 0 && levcount < 10);
      } while (new_count > 0 && levcount < mesh->levmx);
   }




/*
      if (Q[ic] < .02 && level[ic] > 0) {
         nr = nrht[ic];
         nt = ntop[ic];
         nrt = nrht[nt];
         if (i[ic]%2 == 0 && j[ic]%2 == 0 &&
             level[ic]==level[nr]  && Q[nr]  < .02 &&
             level[ic]==level[nt]  && Q[nt]  < .02 &&
             level[ic]==level[nrt] && Q[nrt] < .02 ) {
            jcount++;
         }
      }
*/

   cpu_time_refine_potential += cpu_timer_stop(tstart_cpu);
}

void State::calc_refine_potential_local(Mesh *mesh, vector<int> &mpot,int &icount, int &jcount)
{
// MPI_Barrier(MPI_COMM_WORLD);  

   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   size_t &ncells       = mesh->ncells;
   vector<int> &nlft    = mesh->nlft;
   vector<int> &nrht    = mesh->nrht;
   vector<int> &nbot    = mesh->nbot;
   vector<int> &ntop    = mesh->ntop;
   vector<int> &level   = mesh->level;

   int nl, nr, nt, nb; 
   double Hic, Hl, Hr, Hb, Ht;
   double Uic, Ul, Ur, Ub, Ut;
   double Vic, Vl, Vr, Vb, Vt;

   double duplus1, duminus1, duhalf1;
   double duplus2, duminus2, duhalf2;

   double qpot, qmax;

   icount=0;
   jcount=0;

   L7_Update(&H[0], L7_REAL, mesh->cell_handle);
   L7_Update(&U[0], L7_REAL, mesh->cell_handle);
   L7_Update(&V[0], L7_REAL, mesh->cell_handle);

   for (uint ic=0; ic<ncells; ic++) {

      Hic = H[ic];
      Uic = U[ic];
      Vic = V[ic];

      nl = nlft[ic];
      Hl = H[nl];
      Ul = U[nl];
      Vl = V[nl];

      if (level[nl] > level[ic]){
         int nlt = ntop[nl];
         Hl = 0.5 * (Hl + H[nlt]);
      }

      nr = nrht[ic];
      Hr = H[nr];
      Ur = U[nr];
      Vr = V[nr];
      if (level[nr] > level[ic]){
         int nrt = ntop[nr];
         Hr = 0.5 * (Hr + H[nrt]);
      }

      nb = nbot[ic];
      Hb = H[nb];
      Ub = U[nb];
      Vb = V[nb];
      if (level[nb] > level[ic]){
         int nbr = nrht[nb];
         Hb = 0.5 * (Hb + H[nbr]);
      }

      nt = ntop[ic];
      Ht = H[nt];
      Ut = U[nt];
      Vt = V[nt];

      if (level[nt] > level[ic]){
         int ntr = nrht[nt];
         Ht = 0.5 * (Ht + H[ntr]);
      }

      duplus1 = Hr-Hic;
      duplus2 = Ur-Uic;
      duhalf1 = Hic-Hl;
      duhalf2 = Uic-Ul;

      qmax = -1000.0;

      qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hl;
      duminus2 = Uic-Ul;
      duhalf1 = Hr-Hic;
      duhalf2 = Ur-Uic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duplus1 = Ht-Hic;
      duplus2 = Vt-Vic;
      duhalf1 = Hic-Hb;
      duhalf2 = Vic-Vb;

      qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hb;
      duminus2 = Vic-Vb;
      duhalf1 = Ht-Hic;
      duhalf2 = Vt-Vic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      mpot[ic]=0;
      if (qmax > .10 && level[ic] < mesh->levmx) {
         mpot[ic]=1;
         icount++;
      }

   }

/*
      if (Q[ic] < .02 && level[ic] > 0) {
         nr = nrht[ic];
         nt = ntop[ic];
         nrt = nrht[nt];
         if (i[ic]%2 == 0 && j[ic]%2 == 0 &&
             level[ic]==level[nr]  && Q[nr]  < .02 &&
             level[ic]==level[nt]  && Q[nt]  < .02 &&
             level[ic]==level[nrt] && Q[nrt] < .02 ) {
            jcount++;
         }
      }
*/

   cpu_time_refine_potential += cpu_timer_stop(tstart_cpu);
}

void State::gpu_calc_refine_potential(cl_command_queue command_queue, Mesh *mesh, cl_mem dev_mpot, cl_mem dev_result, cl_mem dev_ioffset)
{
   cl_event refine_potential_event;
   cl_event refine_smooth_event;

   size_t &ncells       = mesh->ncells;
   int &levmx           = mesh->levmx;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_mpot);
   assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size = global_work_size/local_work_size;
     /*
     __kernel void refine_potential
              const int    ncells,   // 0  Total number of cells.
              const int    levmx,    // 1  Maximum level
     __global       real  *H,        // 2
     __global       real  *U,        // 3
     __global       real  *V,        // 4
     __global const int   *nlft,     // 5  Array of left neighbors.
     __global const int   *nrht,     // 6  Array of right neighbors.
     __global const int   *ntop,     // 7  Array of bottom neighbors.
     __global const int   *nbot,     // 8  Array of top neighbors.
     __global const int   *level,    // 9  Array of level information.
     __global const int   *celltype, // 10 Array of level information.
     __global       int   *mpot,     // 11  Array of mesh potential information.
     __global       int   *ioffset,  // 12
     __global const real  *lev_dx,   // 13
     __global const real  *lev_dy,   // 14
     __global const real  *result,   // 15
     __local        real4 *tile,     // 16  Tile size in real4.
     __local        int8  *itile)    // 17  Tile size in int8.
     */
   ezcl_set_kernel_arg(kernel_refine_potential, 0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_refine_potential, 1, sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_refine_potential, 2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_refine_potential, 3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_refine_potential, 4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_refine_potential, 5, sizeof(cl_mem),  (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_refine_potential, 6, sizeof(cl_mem),  (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_refine_potential, 7, sizeof(cl_mem),  (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_refine_potential, 8, sizeof(cl_mem),  (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_refine_potential, 9, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_refine_potential,10, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_refine_potential,11, sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_refine_potential,12, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_refine_potential,13, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_refine_potential,14, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_refine_potential,15, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_refine_potential,16, local_work_size*sizeof(cl_real4),    NULL);
   ezcl_set_kernel_arg(kernel_refine_potential,17, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_refine_potential, 1, NULL, &global_work_size, &local_work_size,
       &refine_potential_event);

   mesh->gpu_rezone_count(command_queue, block_size, local_work_size, dev_ioffset, dev_result);

// XXX Maybe figure out way to get rid of reading off device?? XXX


   size_t result = 0;
   ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int), &result, NULL);

//   printf("result = %d after first refine potential\n",(result-ncells));
//   int which_smooth = 1;

//   sleep(1);

   if( (result-ncells) > 0) {
      cl_mem dev_newcount   = ezcl_malloc(NULL, &block_size, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);

      do {
         
      ezcl_set_kernel_arg(kernel_refine_smooth, 0, sizeof(cl_int),  (void *)&ncells);
      ezcl_set_kernel_arg(kernel_refine_smooth, 1, sizeof(cl_int),  (void *)&levmx);
      ezcl_set_kernel_arg(kernel_refine_smooth, 2, sizeof(cl_mem),  (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_refine_smooth, 3, sizeof(cl_mem),  (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_refine_smooth, 4, sizeof(cl_mem),  (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_refine_smooth, 5, sizeof(cl_mem),  (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_refine_smooth, 6, sizeof(cl_mem),  (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_refine_smooth, 7, sizeof(cl_mem),  (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_refine_smooth, 8, sizeof(cl_mem),  (void *)&dev_mpot);
      ezcl_set_kernel_arg(kernel_refine_smooth, 9, sizeof(cl_mem),  (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_refine_smooth,10, sizeof(cl_mem),  (void *)&dev_newcount);
      ezcl_set_kernel_arg(kernel_refine_smooth,11, sizeof(cl_mem),  (void *)&dev_result);
      ezcl_set_kernel_arg(kernel_refine_smooth,12, local_work_size*sizeof(cl_int2),    NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_refine_smooth, 1, NULL, &global_work_size, &local_work_size, &refine_smooth_event);

      mesh->gpu_rezone_count(command_queue, block_size, local_work_size, dev_newcount, dev_result);

      ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int), &result, NULL);

//      printf("result = %d after %d refine smooths\n",result,which_smooth);
//      sleep(1);
//      which_smooth++;
      } while(result > 0);

      ezcl_device_memory_remove(dev_newcount);

      mesh->gpu_rezone_count(command_queue, block_size, local_work_size, dev_ioffset, dev_result);

      gpu_time_refine_potential  += ezcl_timer_calc(&refine_smooth_event,  &refine_smooth_event);
   }

   gpu_time_refine_potential  += ezcl_timer_calc(&refine_potential_event,  &refine_potential_event);
   //gpu_time_refine_potential  += (ezcl_timer_calc(&refine_potential_event,  &refine_potential_event)*1000000)/ncells;
}

void State::gpu_calc_refine_potential_local(cl_command_queue command_queue, Mesh *mesh, cl_mem dev_mpot, cl_mem dev_result, cl_mem dev_ioffset)
{
   cl_event refine_potential_event;

// MPI_Barrier(MPI_COMM_WORLD);

   struct timeval tstart_gpu;

   cpu_timer_start(&tstart_gpu);

   size_t &ncells       = mesh->ncells;
   int &levmx           = mesh->levmx;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_mpot);
   assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   vector<real> H_tmp(mesh->ncells_ghost,0.0);
   vector<real> U_tmp(mesh->ncells_ghost,0.0);
   vector<real> V_tmp(mesh->ncells_ghost,0.0);
   ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real), &H_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real), &U_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real), &V_tmp[0], NULL);

   L7_Update(&H_tmp[0], L7_REAL, mesh->cell_handle);
   L7_Update(&U_tmp[0], L7_REAL, mesh->cell_handle);
   L7_Update(&V_tmp[0], L7_REAL, mesh->cell_handle);

   size_t nghost_local = mesh->ncells_ghost - ncells;
   //fprintf(mesh->fp,"%d: sizes are ncells %d nghost %d ncells_ghost %d\n",mesh->mype,ncells,nghost_local,mesh->ncells_ghost);


   if (mesh->numpe > 1) {
      cl_mem dev_H_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      cl_mem dev_U_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      cl_mem dev_V_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      ezcl_enqueue_write_buffer(command_queue, dev_H_add, CL_FALSE, 0, nghost_local*sizeof(cl_real), (void*)&H_tmp[ncells],     NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_U_add, CL_FALSE, 0, nghost_local*sizeof(cl_real), (void*)&U_tmp[ncells],     NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_V_add, CL_TRUE,  0, nghost_local*sizeof(cl_real), (void*)&V_tmp[ncells],     NULL);

      size_t ghost_local_work_size = 32;
      size_t ghost_global_work_size = ((nghost_local + ghost_local_work_size - 1) /ghost_local_work_size) * ghost_local_work_size;

      // Fill in ghost
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 1, sizeof(cl_int), (void *)&nghost_local);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 2, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 3, sizeof(cl_mem), (void *)&dev_H_add);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 4, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 5, sizeof(cl_mem), (void *)&dev_U_add);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 6, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_ghost_data, 7, sizeof(cl_mem), (void *)&dev_V_add);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_ghost_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, NULL);
   }

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

     /*
     __kernel void refine_potential
              const int    ncells,   // 0  Total number of cells.
              const int    levmx,    // 1  Maximum level
     __global       real  *H,        // 2
     __global       real  *U,        // 3
     __global       real  *V,        // 4
     __global const int   *nlft,     // 5  Array of left neighbors.
     __global const int   *nrht,     // 6  Array of right neighbors.
     __global const int   *ntop,     // 7  Array of bottom neighbors.
     __global const int   *nbot,     // 8  Array of top neighbors.
     __global const int   *level,    // 9  Array of level information.
     __global const int   *celltype, // 10  Array of celltype information.
     __global       int   *mpot,     // 11  Array of mesh potential information.
     __global       int   *ioffset,  // 12
     __global const real  *lev_dx,   // 13
     __global const real  *lev_dy,   // 14
     __local        real4 *tile,     // 15  Tile size in real4.
     __local        int8  *itile)    // 16  Tile size in int8.
     */
   ezcl_set_kernel_arg(kernel_refine_potential, 0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_refine_potential, 1, sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_refine_potential, 2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_refine_potential, 3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_refine_potential, 4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_refine_potential, 5, sizeof(cl_mem),  (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_refine_potential, 6, sizeof(cl_mem),  (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_refine_potential, 7, sizeof(cl_mem),  (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_refine_potential, 8, sizeof(cl_mem),  (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_refine_potential, 9, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_refine_potential,10, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_refine_potential,11, sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_refine_potential,12, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_refine_potential,13, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_refine_potential,14, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_refine_potential,15, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_refine_potential,16, local_work_size*sizeof(cl_real4),    NULL);
   ezcl_set_kernel_arg(kernel_refine_potential,17, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_refine_potential, 1, NULL, &global_work_size, &local_work_size,
       &refine_potential_event);

   gpu_time_refine_potential += ((long)(cpu_timer_stop(tstart_gpu)*1e9));  

   gpu_time_refine_potential += ezcl_timer_calc(&refine_potential_event,  &refine_potential_event);
}

double State::mass_sum(Mesh *mesh, bool enhanced_precision_sum)
{
   size_t &ncells = mesh->ncells;
   vector<int> &celltype = mesh->celltype;
   vector<int> &level    = mesh->level;

   struct timeval tstart_cpu;
   double summer = 0.0;

   cpu_timer_start(&tstart_cpu);

   if (enhanced_precision_sum) {
      double correction, corrected_next_term, new_sum;
      correction = 0.0;

      for (uint ic = 0; ic < ncells; ic++) {
         if (celltype[ic] == REAL_CELL) {
            //  Exclude boundary cells.
            corrected_next_term= H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]] + correction;
            new_sum            = summer + corrected_next_term;
            correction         = corrected_next_term - (new_sum - summer);
            summer             = new_sum;
         }
      }
   } else {
      for (uint ic=0; ic < ncells; ic++){
         if (celltype[ic] == REAL_CELL) {
            summer += H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]];
         }
      }
   }

   cpu_time_mass_sum += cpu_timer_stop(tstart_cpu);

   return(summer);
}

double State::mass_sum_local(Mesh *mesh, bool enhanced_precision_sum)
{

   size_t &ncells = mesh->ncells;
   vector<int> &celltype = mesh->celltype;
   vector<int> &level    = mesh->level;

   struct timeval tstart_cpu;
   double summer = 0.0;
   double total_sum = 0.0;

   cpu_timer_start(&tstart_cpu);

   if (enhanced_precision_sum) {
      double correction, corrected_next_term, new_sum;
      correction = 0.0;
      struct esum_type local, global;

      local.sum = 0.0;
      local.correction = 0.0;
      for (uint ic = 0; ic < ncells; ic++) {
         if (celltype[ic] == REAL_CELL) {
            //  Exclude boundary cells.
            corrected_next_term= H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]] + local.correction;
            new_sum            = local.sum + local.correction;
            local.correction   = corrected_next_term - (new_sum - local.sum);
            local.sum          = new_sum;
         }
      }
      MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KAHAN_SUM, MPI_COMM_WORLD);
      total_sum = global.sum + global.correction;
   } else {
      for (uint ic=0; ic < ncells; ic++){
         if (celltype[ic] == REAL_CELL) {
            summer += H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]];
         }
      }
      MPI_Allreduce(&summer, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   cpu_time_mass_sum += cpu_timer_stop(tstart_cpu);

   return(total_sum);
}

double State::gpu_mass_sum(cl_command_queue command_queue, Mesh *mesh, bool enhanced_precision_sum)
{
   cl_event mass_sum_stage1_event;
   cl_event mass_sum_stage2_event;
   cl_event start_read_event;

   size_t &ncells       = mesh->ncells;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_level    = mesh->dev_level;

   assert(dev_H);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);
   assert(dev_celltype);

   size_t one = 1;
   cl_mem dev_mass_sum, dev_redscratch;
   double gpu_mass_sum;

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   if (enhanced_precision_sum) {
      dev_mass_sum = ezcl_malloc(NULL, &one,    sizeof(cl_real2), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, &block_size, sizeof(cl_real2), CL_MEM_READ_WRITE, 0);

        /*
        __kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *array,      // 1   Array to be reduced.
                __global       int   *level,      // 2
                __global       int   *levdx,      // 3
                __global       int   *levdy,      // 4
                __global       int   *celltype,   // 5
                __global       real  *redscratch, // 6   Final result of operation.
                __local        real  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 8, local_work_size*sizeof(cl_real2), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, &mass_sum_stage1_event);

        /*
        __kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *redscratch, // 1   Array to be reduced.
                __local        real  *tile)       // 2
        */

      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 3, local_work_size*sizeof(cl_real2), NULL);

      if (block_size > 1) {
         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, &mass_sum_stage2_event);
      }

      real2 mass_sum;
      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real2), &mass_sum, &start_read_event);               
      gpu_mass_sum = mass_sum.s0 + mass_sum.s1;
   } else {
      dev_mass_sum = ezcl_malloc(NULL, &one,    sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, &block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

        /*
        __kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *array,      // 1   Array to be reduced.
                __global       int   *level,      // 2
                __global       int   *levdx,      // 3
                __global       int   *levdy,      // 4
                __global       int   *celltype,   // 5
                __global       real  *redscratch, // 6   Final result of operation.
                __local        real  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 8, local_work_size*sizeof(cl_real), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, &mass_sum_stage1_event);

        /*
        __kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *redscratch, // 1   Array to be reduced.
                __local        real  *tile)       // 2
        */

      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 3, local_work_size*sizeof(cl_real), NULL);

      if (block_size > 1) {
         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, &mass_sum_stage2_event);
      }

      real mass_sum;
      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real), &mass_sum, &start_read_event);               
      gpu_mass_sum = mass_sum;
   }

   //gpu_time_read += ezcl_timer_calc(start_read_event, start_read_event);

   ezcl_device_memory_remove(dev_redscratch);
   ezcl_device_memory_remove(dev_mass_sum);

   gpu_time_mass_sum     += ezcl_timer_calc(&mass_sum_stage1_event, &mass_sum_stage1_event); 
   if (block_size > 1) {
      gpu_time_mass_sum  += ezcl_timer_calc(&mass_sum_stage2_event, &mass_sum_stage2_event);
   }
   gpu_time_read         += ezcl_timer_calc(&start_read_event,      &start_read_event);

   return(gpu_mass_sum);
}

double State::gpu_mass_sum_local(cl_command_queue command_queue, Mesh *mesh, bool enhanced_precision_sum)
{
   cl_event mass_sum_stage1_event;
   cl_event mass_sum_stage2_event;
   cl_event start_read_event;

   size_t &ncells       = mesh->ncells;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_level    = mesh->dev_level;

   assert(dev_H);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);
   assert(dev_celltype);

   size_t one = 1;
   cl_mem dev_mass_sum, dev_redscratch;
   double total_sum;

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   if (enhanced_precision_sum) {
      dev_mass_sum = ezcl_malloc(NULL, &one,    sizeof(cl_real2), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, &block_size, sizeof(cl_real2), CL_MEM_READ_WRITE, 0);

        /*
     __   kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *array,      // 1   Array to be reduced.
                __global       int   *level,      // 2
                __global       int   *levdx,      // 3
                __global       int   *levdy,      // 4
                __global       int   *celltype,   // 5
                __global       real  *redscratch, // 6   Final result of operation.
                __local        real  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 8, local_work_size*sizeof(cl_real2), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, &mass_sum_stage1_event);

        /*
        __kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *redscratch, // 1   Array to be reduced.
                __local        real  *tile)       // 2
        */

      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 3, local_work_size*sizeof(cl_real2), NULL);

      if (block_size > 1) {
         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, &mass_sum_stage2_event);
      }

      struct esum_type local, global;
      real2 mass_sum;

      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real2), &mass_sum, &start_read_event);
      local.sum = mass_sum.s0;
      local.correction = mass_sum.s1;
      MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KAHAN_SUM, MPI_COMM_WORLD);
      total_sum = global.sum + global.correction;
   } else {
      dev_mass_sum = ezcl_malloc(NULL, &one,    sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, &block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

        /*
     __   kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *array,      // 1   Array to be reduced.
                __global       int   *level,      // 2
                __global       int   *levdx,      // 3
                __global       int   *levdy,      // 4
                __global       int   *celltype,   // 5
                __global       real  *redscratch, // 6   Final result of operation.
                __local        real  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 8, local_work_size*sizeof(cl_real), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, &mass_sum_stage1_event);

        /*
        __kernel void reduce_sum_cl(
                         const int    isize,      // 0
                __global       int   *redscratch, // 1   Array to be reduced.
                __local        real  *tile)       // 2
        */

      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 3, local_work_size*sizeof(cl_real), NULL);

      if (block_size > 1) {
         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, &mass_sum_stage2_event);
      }

      double local_sum, global_sum;
      real mass_sum;

      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real), &mass_sum, &start_read_event);
      local_sum = mass_sum;
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      total_sum = global_sum;
   }

   ezcl_device_memory_remove(dev_redscratch);
   ezcl_device_memory_remove(dev_mass_sum);

   gpu_time_set_timestep      += ezcl_timer_calc(&mass_sum_stage1_event, &mass_sum_stage1_event); 
   if (block_size > 1) {
      gpu_time_set_timestep   += ezcl_timer_calc(&mass_sum_stage2_event, &mass_sum_stage2_event);
   }
   gpu_time_read              += ezcl_timer_calc(&start_read_event,      &start_read_event);


   return(total_sum);
}

void State::allocate_device_memory(size_t ncells)
{
   dev_H = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   dev_U = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   dev_V = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
}
void State::resize_old_device_memory(size_t ncells)
{
   ezcl_device_memory_remove(dev_H);
   ezcl_device_memory_remove(dev_U);
   ezcl_device_memory_remove(dev_V);
   dev_H = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   dev_U = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
   dev_V = ezcl_malloc(NULL, &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
}

void State::output_timing_info(Mesh *mesh, int do_cpu_calc, int do_gpu_calc, long gpu_time_count_BCs, double elapsed_time)
{
   int &mype  = mesh->mype;
   int &numpe = mesh->numpe;

   double cpu_time_compute;
   long gpu_time_compute;

   double cpu_elapsed_time=0.0;
   long gpu_elapsed_time=0;

   if (! mesh->parallel) {
      //  Output timing information.
      if (mype == 0) {
         printf("\n");
         printf("~~~~~~~~~~~~~~~~ Device timing information ~~~~~~~~~~~~~~~~~~\n");
      }
      if (do_cpu_calc) {
         cpu_time_compute = get_cpu_time_apply_BCs() +
                            get_cpu_time_set_timestep() +
                            get_cpu_time_finite_difference() +
                            get_cpu_time_refine_potential() +
                            get_cpu_time_rezone_all() +
                            mesh->get_cpu_time_rezone_all() +
                            mesh->get_cpu_time_calc_neighbors() +
                            get_cpu_time_mass_sum() +
                            mesh->cpu_time_partition;
         cpu_elapsed_time =                  cpu_time_compute;

         if (mype == 0) {
            printf("CPU: Device compute           time was\t%8.4f \ts\n",     cpu_time_compute);
            printf("CPU:  state->set_timestep      time was\t %8.4f\ts\n",     get_cpu_time_set_timestep() );
            printf("CPU:  state->apply_BCs         time was\t %8.4f\ts\n",     get_cpu_time_apply_BCs() );
            printf("CPU:  state->finite_difference time was\t %8.4f\ts\n",     get_cpu_time_finite_difference() );
            printf("CPU:  mesh->refine_potential   time was\t %8.4f\ts\n",     get_cpu_time_refine_potential() );
            printf("CPU:  mesh->rezone_all         time was\t %8.4f\ts\n",     (get_cpu_time_rezone_all() + mesh->get_cpu_time_rezone_all() ) );
            printf("CPU:  mesh->partition_cells    time was\t %8.4f\ts\n",     mesh->cpu_time_partition);
            printf("CPU:  mesh->calc_neighbors     time was\t %8.4f\ts\n",     mesh->get_cpu_time_calc_neighbors() );
            printf("CPU:  mass_sum                 time was\t %8.4f\ts\n",     get_cpu_time_mass_sum() );
            printf("=============================================================\n");
            printf("Profiling: Total CPU          time was\t%8.4f\ts or\t%4.2f min\n",  cpu_elapsed_time, cpu_elapsed_time/60.0);      
            printf("=============================================================\n\n");
         }
      }
      if (do_gpu_calc) {
         gpu_time_compute = get_gpu_time_apply_BCs() +
                            get_gpu_time_set_timestep() +
                            get_gpu_time_finite_difference() +
                            get_gpu_time_refine_potential() +
                            mesh->get_gpu_time_reduction_scan() +
                            get_gpu_time_rezone_all() +
                            mesh->get_gpu_time_rezone_all() +
                            mesh->get_gpu_time_hash_setup() +
                            mesh->get_gpu_time_calc_neighbors() +
                            get_gpu_time_mass_sum() +
                            gpu_time_count_BCs;
         gpu_elapsed_time   = get_gpu_time_write() + gpu_time_compute + get_gpu_time_read();

         if (mype == 0) {
            printf("GPU: Write to device          time was\t%8.4f\ts\n",    (double) get_gpu_time_write()         * 1.0e-9); /* Convert nanoseconds to msecs */
            printf("GPU: Read from device         time was\t%8.4f\ts\n",    (double) get_gpu_time_read()          * 1.0e-9);
            printf("GPU: Device compute           time was\t%8.4f\ts\n",    (double) gpu_time_compute             * 1.0e-9);
            printf("GPU:  kernel_set_timestep      time was\t %8.4f\ts\n",    (double) get_gpu_time_set_timestep()      * 1.0e-9);
            printf("GPU:  kernel_calc_finite_diff  time was\t %8.4f\ts\n",    (double) get_gpu_time_finite_difference() * 1.0e-9);
            printf("GPU:  kernel_refine_potential  time was\t %8.4f\ts\n",    (double) get_gpu_time_refine_potential()  * 1.0e-9);
            printf("GPU:  kernel_reduction_scan    time was\t %8.4f\ts\n",    (double) mesh->get_gpu_time_reduction_scan()    * 1.0e-9);
            printf("GPU:  kernel_rezone_all        time was\t %8.4f\ts\n",    (double) (get_gpu_time_rezone_all() + mesh->get_gpu_time_rezone_all() ) * 1.0e-9);
            printf("GPU:  kernel_hash_setup        time was\t %8.4f\ts\n",    (double) mesh->get_gpu_time_hash_setup()        * 1.0e-9);
            printf("GPU:  kernel_calc_neighbors    time was\t %8.4f\ts\n",    (double) mesh->get_gpu_time_calc_neighbors()    * 1.0e-9);
            printf("GPU:  kernel_mass_sum          time was\t %8.4f\ts\n",    (double) get_gpu_time_mass_sum()          * 1.0e-9);
            if (! mesh->have_boundary) {
               printf("GPU:  kernel_count_BCs         time was\t %8.4f\ts\n",    (double) gpu_time_count_BCs        * 1.0e-9);
            }
            printf("=============================================================\n");
            printf("Profiling: Total GPU          time was\t%8.4f\ts or\t%4.2f min\n", (double) gpu_elapsed_time * 1.0e-9, (double) gpu_elapsed_time * 1.0e-9/60.0);
            printf("=============================================================\n\n");
         }
      }
      if (mype == 0) {
         printf("=============================================================\n");
         printf("Profiling: Total              time was\t%8.4f\ts or\t%4.2f min\n",    elapsed_time,     elapsed_time/60.0);
         if (do_cpu_calc && do_gpu_calc) {
            printf("Speed-up:                             \t%8.4f\n",   cpu_elapsed_time*1.0e9/gpu_elapsed_time);
         }
         printf("=============================================================\n");
      }
   } else {
      if (do_cpu_calc) {
         cpu_time_compute = get_cpu_time_apply_BCs() +
                            get_cpu_time_set_timestep() +
                            get_cpu_time_finite_difference() +
                            get_cpu_time_refine_potential() +
                            get_cpu_time_rezone_all() +
                            mesh->get_cpu_time_rezone_all() +
                            mesh->get_cpu_time_calc_neighbors() +
                            get_cpu_time_mass_sum() +
                            mesh->cpu_time_partition;
         cpu_elapsed_time =                  cpu_time_compute;

         if (mype == 0) printf("\nCPU: Parallel timings\n\n");

         parallel_timer_output(numpe,mype,"CPU: Device compute           time was" ,cpu_time_compute);
         parallel_timer_output(numpe,mype,"CPU:  state->set_timestep      time was",get_cpu_time_set_timestep() );
         parallel_timer_output(numpe,mype,"CPU:  state->apply_BCs         time was",get_cpu_time_apply_BCs() );
         parallel_timer_output(numpe,mype,"CPU:  state->finite_difference time was",get_cpu_time_finite_difference() );
         parallel_timer_output(numpe,mype,"CPU:  mesh->refine_potential   time was",get_cpu_time_refine_potential() );
         parallel_timer_output(numpe,mype,"CPU:  mesh->rezone_all         time was",(get_cpu_time_rezone_all() + mesh->get_cpu_time_rezone_all() ) );
         parallel_timer_output(numpe,mype,"CPU:  mesh->partition_cells    time was",mesh->cpu_time_partition);
         parallel_timer_output(numpe,mype,"CPU:  mesh->calc_neighbors     time was",mesh->get_cpu_time_calc_neighbors() );
         parallel_timer_output(numpe,mype,"CPU:  mass_sum                 time was",get_cpu_time_mass_sum() );

         if (mype == 0) printf("=============================================================\n");

         parallel_timer_output(numpe,mype,"Profiling: Total CPU          time was",cpu_elapsed_time);

         if (mype == 0) printf("=============================================================\n\n");
      }
      if (do_gpu_calc) {
         gpu_time_compute = get_gpu_time_apply_BCs() +
                            get_gpu_time_set_timestep() +
                            get_gpu_time_finite_difference() +
                            get_gpu_time_refine_potential() +
                            mesh->get_gpu_time_reduction_scan() +
                            get_gpu_time_rezone_all() +
                            mesh->get_gpu_time_rezone_all() +
                            mesh->get_gpu_time_hash_setup() +
                            mesh->get_gpu_time_calc_neighbors() +
                            get_gpu_time_mass_sum() +
                            gpu_time_count_BCs;
         gpu_elapsed_time   = get_gpu_time_write() + gpu_time_compute + get_gpu_time_read();

         if (mype == 0) printf("\nGPU: Parallel timings\n\n");

         parallel_timer_output(numpe,mype,"GPU: Write to device          time was", (double) get_gpu_time_write()         * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU: Read from device         time was", (double) get_gpu_time_read()          * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU: Device compute           time was", (double) gpu_time_compute             * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_set_timestep      time was",(double) get_gpu_time_set_timestep()      * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_calc_finite_diff  time was",(double) get_gpu_time_finite_difference() * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_refine_potential  time was",(double) get_gpu_time_refine_potential()  * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_reduction_scan    time was",(double) mesh->get_gpu_time_reduction_scan()     * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_rezone_all        time was",(double) (get_gpu_time_rezone_all() + mesh->get_gpu_time_rezone_all() ) * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_hash_setup        time was",(double) mesh->get_gpu_time_hash_setup()         * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_calc_neighbors    time was",(double) mesh->get_gpu_time_calc_neighbors()     * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_mass_sum          time was",(double) get_gpu_time_mass_sum()          * 1.0e-9 );

         if (! mesh->have_boundary) {
            parallel_timer_output(numpe,mype,"GPU:  kernel_count_BCs         time was",(double) gpu_time_count_BCs        * 1.0e-9 );
         }

         if (mype == 0) printf("=============================================================\n");

         parallel_timer_output(numpe,mype,"Profiling: Total GPU          time was",gpu_elapsed_time * 1.0e-9 );

         if (mype == 0) printf("=============================================================\n\n");
      }

      if (mype == 0) printf("=============================================================\n");

      parallel_timer_output(numpe,mype,"Profiling: Total              time was",elapsed_time );
      parallel_timer_output(numpe,mype,"Parallel Speed-up:                    ",cpu_elapsed_time*1.0e9/gpu_elapsed_time);

      if (mype == 0) printf("=============================================================\n");
   }
}

void State::parallel_timer_output(int numpe, int mype, const char *string, double local_time)
{
   vector<double> global_times(numpe);
   MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_times[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   if (mype == 0) {
      printf("%s\t",string);
      for(int ip = 0; ip < numpe; ip++){
         printf("%8.4f\t", global_times[ip]);
      }
      printf("s\n");
   }
}

