/*
 *  Copyright (c) 2011-2013, Los Alamos National Security, LLC.
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
#include "mesh/mesh.h"
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include "state.h"
#include "kdtree/KDTree.h"
//#include "reorder.h"
#include "timer/timer.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_REPROBLAS
#ifdef __cplusplus
extern "C"
{
#endif
#ifdef HAVE_MPI
#ifdef OMPI_MPI_H
#define OMPI_SKIP_MPICXX
#endif
#include <MPIndexedFP.h>
#else
#include <IndexedFP.h>
#endif
#ifdef __cplusplus
}
#endif
#endif

#undef DEBUG
//#define DEBUG 1
#undef DEBUG_RESTORE_VALS
#define TIMING_LEVEL 2

#if defined(MINIMUM_PRECISION)
#define ZERO 0.0f
#define ONE 1.0f
#define HALF 0.5f
#define EPSILON 1.0f-30
#define STATE_EPS        15.0
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define EPSILON 1.0e-30
#define STATE_EPS        .02
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(FULL_PRECISION)
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define EPSILON 1.0e-30
#define STATE_EPS        .02
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10
#define COARSEN_GRADIENT 0.05
#define REFINE_HALF 0.5
#define REFINE_NEG_THOUSAND -1000.0

#endif

typedef unsigned int uint;

#ifdef HAVE_OPENCL
#include "state_kernel.inc"
#endif

struct esum_type{
   double sum;
   double correction;
};
#ifdef HAVE_MPI
MPI_Datatype MPI_TWO_DOUBLES;
MPI_Op KAHAN_SUM;
int commutative = 1;
void kahan_sum(struct esum_type *in, struct esum_type *inout, int *len, MPI_Datatype *MPI_TWO_DOUBLES);
#endif

int save_ncells;

#define CONSERVED_EQNS

#define SQR(x) ( x*x )
#define MIN3(x,y,z) ( min( min(x,y), z) )

#ifdef HAVE_OPENCL
cl_kernel kernel_set_timestep;
cl_kernel kernel_reduction_min;
cl_kernel kernel_copy_state_data;
cl_kernel kernel_copy_state_ghost_data;
cl_kernel kernel_apply_boundary_conditions;
cl_kernel kernel_apply_boundary_conditions_local;
cl_kernel kernel_apply_boundary_conditions_ghost;
cl_kernel kernel_calc_finite_difference;
cl_kernel kernel_refine_potential;
cl_kernel kernel_reduce_sum_mass_stage1of2;
cl_kernel kernel_reduce_sum_mass_stage2of2;
cl_kernel kernel_reduce_epsum_mass_stage1of2;
cl_kernel kernel_reduce_epsum_mass_stage2of2;
#endif

State::State(Mesh *mesh_in)
{
   cpu_time_apply_BCs          = 0.0;
   cpu_time_set_timestep       = 0.0;
   cpu_time_finite_difference  = 0.0;
   cpu_time_refine_potential   = 0.0;
     cpu_time_calc_mpot        = 0.0;
   cpu_time_mass_sum           = 0.0;

   gpu_time_apply_BCs          = 0L;
   gpu_time_set_timestep       = 0L;
   gpu_time_finite_difference  = 0L;
   gpu_time_refine_potential   = 0L;
      gpu_time_calc_mpot       = 0L;
   gpu_time_mass_sum           = 0L;
   gpu_time_read               = 0L;
   gpu_time_write              = 0L;

   mesh = mesh_in;

#ifdef HAVE_MPI
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init){
      MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_TWO_DOUBLES);
      MPI_Type_commit(&MPI_TWO_DOUBLES);
      MPI_Op_create((MPI_User_function *)kahan_sum, commutative, &KAHAN_SUM);
      // FIXME add fini and set size
      if (mesh->parallel) state_memory.pinit(MPI_COMM_WORLD, 2L * 1024 * 1024 * 1024);
   }
#ifdef HAVE_REPROBLAS
   RMPI_Init(); // Initialize Reproducible MPI
#endif
#endif
}

void State::init(int do_gpu_calc)
{
   if (do_gpu_calc) {
#ifdef HAVE_OPENCL
      cl_context context = ezcl_get_context();

      if (ezcl_get_compute_device() == COMPUTE_DEVICE_ATI) printf("Starting compile of kernels in state\n");
      const char *defines = NULL;
      cl_program program                 = ezcl_create_program_wsource(context, defines, state_kern_source);

      kernel_set_timestep                    = ezcl_create_kernel_wprogram(program, "set_timestep_cl");
      kernel_reduction_min                   = ezcl_create_kernel_wprogram(program, "finish_reduction_min_cl");
      kernel_copy_state_data                 = ezcl_create_kernel_wprogram(program, "copy_state_data_cl");
      kernel_copy_state_ghost_data           = ezcl_create_kernel_wprogram(program, "copy_state_ghost_data_cl");
      kernel_apply_boundary_conditions       = ezcl_create_kernel_wprogram(program, "apply_boundary_conditions_cl");
      kernel_apply_boundary_conditions_local = ezcl_create_kernel_wprogram(program, "apply_boundary_conditions_local_cl");
      kernel_apply_boundary_conditions_ghost = ezcl_create_kernel_wprogram(program, "apply_boundary_conditions_ghost_cl");
      kernel_calc_finite_difference          = ezcl_create_kernel_wprogram(program, "calc_finite_difference_cl");
      kernel_refine_potential                = ezcl_create_kernel_wprogram(program, "refine_potential_cl");
      kernel_reduce_sum_mass_stage1of2       = ezcl_create_kernel_wprogram(program, "reduce_sum_mass_stage1of2_cl");
      kernel_reduce_sum_mass_stage2of2       = ezcl_create_kernel_wprogram(program, "reduce_sum_mass_stage2of2_cl");
      kernel_reduce_epsum_mass_stage1of2     = ezcl_create_kernel_wprogram(program, "reduce_epsum_mass_stage1of2_cl");
      kernel_reduce_epsum_mass_stage2of2     = ezcl_create_kernel_wprogram(program, "reduce_epsum_mass_stage2of2_cl");

      ezcl_program_release(program);
      if (ezcl_get_compute_device() == COMPUTE_DEVICE_ATI) printf("Finishing compile of kernels in state\n");
#endif
   }

   //printf("\nDEBUG -- Calling state memory memory malloc at line %d\n",__LINE__);
   allocate(mesh->ncells);
   //state_memory.memory_report();
   //printf("DEBUG -- Finished state memory memory malloc at line %d\n\n",__LINE__);

}

void State::allocate(size_t ncells)
{
   int flags = 0;
#ifdef HAVE_J7
   if (mesh->parallel) flags = LOAD_BALANCE_MEMORY;
#endif

   H = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), flags, "H");
   U = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), flags, "U");
   V = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), flags, "V");
}

void State::resize(size_t new_ncells){
   size_t current_size = state_memory.get_memory_size(H);
   if (new_ncells > current_size) state_memory.memory_realloc_all(new_ncells);

   //printf("\nDEBUG -- Calling state memory resize at line %d\n",__LINE__);
   //state_memory.memory_report();
   //printf("DEBUG -- Finished state memory resize at line %d\n\n",__LINE__);
}

void State::memory_reset_ptrs(void){
   H = (state_t *)state_memory.get_memory_ptr("H");
   U = (state_t *)state_memory.get_memory_ptr("U");
   V = (state_t *)state_memory.get_memory_ptr("V");

   //printf("\nDEBUG -- Calling state memory reset_ptrs at line %d\n",__LINE__);
   //state_memory.memory_report();
   //printf("DEBUG -- Finished state memory reset_ptrs at line %d\n\n",__LINE__);
}

void State::terminate(void)
{
   state_memory.memory_delete(H);
   state_memory.memory_delete(U);
   state_memory.memory_delete(V);

#ifdef HAVE_OPENCL
   ezcl_device_memory_delete(dev_deltaT);

   gpu_state_memory.memory_delete(dev_H);
   gpu_state_memory.memory_delete(dev_U);
   gpu_state_memory.memory_delete(dev_V);

   ezcl_kernel_release(kernel_set_timestep);
   ezcl_kernel_release(kernel_reduction_min);
   ezcl_kernel_release(kernel_copy_state_data);
   ezcl_kernel_release(kernel_copy_state_ghost_data);
   ezcl_kernel_release(kernel_apply_boundary_conditions);
   ezcl_kernel_release(kernel_apply_boundary_conditions_local);
   ezcl_kernel_release(kernel_apply_boundary_conditions_ghost);
   ezcl_kernel_release(kernel_calc_finite_difference);
   ezcl_kernel_release(kernel_refine_potential);
   ezcl_kernel_release(kernel_reduce_sum_mass_stage1of2);
   ezcl_kernel_release(kernel_reduce_sum_mass_stage2of2);
   ezcl_kernel_release(kernel_reduce_epsum_mass_stage1of2);
   ezcl_kernel_release(kernel_reduce_epsum_mass_stage2of2);
#endif
#ifdef HAVE_MPI
   if (mesh->parallel) state_memory.pfini();
#endif
}

#ifdef HAVE_MPI
void kahan_sum(struct esum_type *in, struct esum_type *inout, int *len, MPI_Datatype *MPI_TWO_DOUBLES)
{
   double corrected_next_term, new_sum;

   corrected_next_term = in->sum +(in->correction+inout->correction);
   new_sum = inout->sum + corrected_next_term;
   inout->correction = corrected_next_term - (new_sum - inout->sum);
   inout->sum = new_sum;

   // Just to block compiler warnings
   if (1==2) printf("DEBUG len %d datatype %lld\n",*len,(long long)(*MPI_TWO_DOUBLES) );
}
#endif

void State::add_boundary_cells(void)
{
   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   // This is for a mesh with no boundary cells -- they are added and
   // the mesh sizes increased
   size_t &ncells        = mesh->ncells;
   vector<int>  &index    = mesh->index;
   vector<spatial_t> &x        = mesh->x;
   vector<spatial_t> &dx       = mesh->dx;
   vector<spatial_t> &y        = mesh->y;
   vector<spatial_t> &dy       = mesh->dy;

   int *i        = mesh->i;
   int *j        = mesh->j;
   int *level    = mesh->level;
   int *celltype = mesh->celltype;
   int *nlft     = mesh->nlft;
   int *nrht     = mesh->nrht;
   int *nbot     = mesh->nbot;
   int *ntop     = mesh->ntop;

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
   H=(state_t *)state_memory.memory_realloc(new_ncells, sizeof(state_t), H);
   U=(state_t *)state_memory.memory_realloc(new_ncells, sizeof(state_t), U);
   V=(state_t *)state_memory.memory_realloc(new_ncells, sizeof(state_t), V);
   //printf("\nDEBUG add_boundary cells\n"); 
   //state_memory.memory_report();
   //printf("DEBUG end add_boundary cells\n\n"); 

   mesh->i        =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), i);
   mesh->j        =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), j);
   mesh->level    =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), level);
   mesh->celltype =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), celltype);
   mesh->nlft     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), nlft);
   mesh->nrht     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), nrht);
   mesh->nbot     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), nbot);
   mesh->ntop     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, sizeof(int), ntop);
   //memory_reset_ptrs();
   i        = mesh->i;
   j        = mesh->j;
   level    = mesh->level;
   celltype = mesh->celltype;
   nlft     = mesh->nlft;
   nrht     = mesh->nrht;
   nbot     = mesh->nbot;
   ntop     = mesh->ntop;

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

void State::apply_boundary_conditions_local(void)
{
   size_t &ncells = mesh->ncells;
   int *nlft = mesh->nlft;
   int *nrht = mesh->nrht;
   int *nbot = mesh->nbot;
   int *ntop = mesh->ntop;

   // This is for a mesh with boundary cells
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (uint ic=0; ic<ncells; ic++) {
      if (mesh->is_left_boundary(ic)) {
         int nr = nrht[ic];
         if (nr < (int)ncells) {
            H[ic] =  H[nr];
            U[ic] = -U[nr];
            V[ic] =  V[nr];
         }
      }
      if (mesh->is_right_boundary(ic))  {
         int nl = nlft[ic];
         if (nl < (int)ncells) {
            H[ic] =  H[nl];
            U[ic] = -U[nl];
            V[ic] =  V[nl];
         }
      }
      if (mesh->is_bottom_boundary(ic)) {
         int nt = ntop[ic];
         if (nt < (int)ncells) {
            H[ic] =  H[nt];
            U[ic] =  U[nt];
            V[ic] = -V[nt];
         }
      }
      if (mesh->is_top_boundary(ic)) {
         int nb = nbot[ic];
         if (nb < (int)ncells) {
            H[ic] =  H[nb];
            U[ic] =  U[nb];
            V[ic] = -V[nb];
         }
      }
   }
}

void State::apply_boundary_conditions_ghost(void)
{

   size_t &ncells = mesh->ncells;
   int *nlft = mesh->nlft;
   int *nrht = mesh->nrht;
   int *nbot = mesh->nbot;
   int *ntop = mesh->ntop;

   // This is for a mesh with boundary cells
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (uint ic=0; ic<ncells; ic++) {
      if (mesh->is_left_boundary(ic)) {
         int nr = nrht[ic];
         if (nr >= (int)ncells) {
            H[ic] =  H[nr];
            U[ic] = -U[nr];
            V[ic] =  V[nr];
         }
      }
      if (mesh->is_right_boundary(ic))  {
         int nl = nlft[ic];
         if (nl >= (int)ncells) {
            H[ic] =  H[nl];
            U[ic] = -U[nl];
            V[ic] =  V[nl];
         }
      }
      if (mesh->is_bottom_boundary(ic)) {
         int nt = ntop[ic];
         if (nt >= (int)ncells) {
            H[ic] =  H[nt];
            U[ic] =  U[nt];
            V[ic] = -V[nt];
         }
      }
      if (mesh->is_top_boundary(ic)) {
         int nb = nbot[ic];
         if (nb >= (int)ncells) {
            H[ic] =  H[nb];
            U[ic] =  U[nb];
            V[ic] = -V[nb];
         }
      }
   }
}

void State::apply_boundary_conditions(void)
{
   size_t &ncells = mesh->ncells;
   int *nlft = mesh->nlft;
   int *nrht = mesh->nrht;
   int *nbot = mesh->nbot;
   int *ntop = mesh->ntop;

   // This is for a mesh with boundary cells
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (uint ic=0; ic<ncells; ic++) {
      if (mesh->is_left_boundary(ic)) {
         int nr = nrht[ic];
         H[ic] =  H[nr];
         U[ic] = -U[nr];
         V[ic] =  V[nr];
      }
      if (mesh->is_right_boundary(ic))  {
         int nl = nlft[ic];
         H[ic] =  H[nl];
         U[ic] = -U[nl];
         V[ic] =  V[nl];
      }
      if (mesh->is_bottom_boundary(ic)) {
         int nt = ntop[ic];
         H[ic] =  H[nt];
         U[ic] =  U[nt];
         V[ic] = -V[nt];
      }
      if (mesh->is_top_boundary(ic)) {
         int nb = nbot[ic];
         H[ic] =  H[nb];
         U[ic] =  U[nb];
         V[ic] = -V[nb];
      }
   }
}

void State::remove_boundary_cells(void)
{
   size_t &ncells = mesh->ncells;
   vector<int> &index    = mesh->index;
   vector<spatial_t> &x       = mesh->x;
   vector<spatial_t> &dx      = mesh->dx;
   vector<spatial_t> &y       = mesh->y;
   vector<spatial_t> &dy      = mesh->dy;

   int *i        = mesh->i;
   int *j        = mesh->j;
   int *level    = mesh->level;
   int *celltype = mesh->celltype;
   int *nlft     = mesh->nlft;
   int *nrht     = mesh->nrht;
   int *nbot     = mesh->nbot;
   int *ntop     = mesh->ntop;

   if(mesh->have_boundary) return;

   // Resize to drop all the boundary cells
   ncells = save_ncells;
   H=(state_t *)state_memory.memory_realloc(save_ncells, sizeof(state_t), H);
   U=(state_t *)state_memory.memory_realloc(save_ncells, sizeof(state_t), U);
   V=(state_t *)state_memory.memory_realloc(save_ncells, sizeof(state_t), V);
   //printf("\nDEBUG remove_boundary cells\n"); 
   //state_memory.memory_report();
   //printf("DEBUG end remove_boundary cells\n\n"); 

   mesh->i        = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), i);
   mesh->j        = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), j);
   mesh->level    = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), level);
   mesh->celltype = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), celltype);
   mesh->nlft     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), nlft);
   mesh->nrht     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), nrht);
   mesh->nbot     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), nbot);
   mesh->ntop     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, sizeof(int), ntop);
   i        = mesh->i;
   j        = mesh->j;
   level    = mesh->level;
   celltype = mesh->celltype;
   nlft     = mesh->nlft;
   nrht     = mesh->nrht;
   nbot     = mesh->nbot;
   ntop     = mesh->ntop;

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

double State::set_timestep(double g, double sigma)
{
   double globalmindeltaT;
   double mindeltaT = 1000.0;
   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   size_t ncells        = mesh->ncells;
#ifdef HAVE_MPI
   int parallel         = mesh->parallel;
#endif
   int *&celltype = mesh->celltype;
   int *&level    = mesh->level;

   int ic;
#ifdef _OPENMP
#pragma omp parallel
   {
      double mymindeltaT = 1000.0;
#pragma omp for
#endif
      for (ic=0; ic<(int)ncells; ic++) {
         if (celltype[ic] == REAL_CELL) {
            int lev = level[ic];
            double wavespeed = sqrt(g*H[ic]);
            double xspeed = (fabs(U[ic])+wavespeed)/mesh->lev_deltax[lev];
            double yspeed = (fabs(V[ic])+wavespeed)/mesh->lev_deltay[lev];
            double deltaT=sigma/(xspeed+yspeed);
#ifdef _OPENMP
            if (deltaT < mymindeltaT) mymindeltaT = deltaT;
#else
            if (deltaT < mindeltaT) mindeltaT = deltaT;
#endif
         }
      }
#ifdef _OPENMP
#pragma omp critical
      if (mymindeltaT < mindeltaT) mindeltaT = mymindeltaT;
   }
#endif

   globalmindeltaT = mindeltaT;
#ifdef HAVE_MPI
   if (parallel) MPI_Allreduce(&mindeltaT, &globalmindeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

   cpu_time_set_timestep += cpu_timer_stop(tstart_cpu);

   return(globalmindeltaT);
}

#ifdef HAVE_OPENCL
double State::gpu_set_timestep(double sigma)
{
   double deltaT, globalmindeltaT;

   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t &ncells       = mesh->ncells;
#ifdef HAVE_MPI
   int &parallel        = mesh->parallel;
#endif
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

   cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);

      /*
      __kernel void set_timestep_cl(
                       const int       ncells,     // 0  Total number of cells.
                       const real_t    sigma,      // 1
              __global const state_t  *H,          // 2
              __global const state_t  *U,          // 3
              __global const state_t  *V,          // 4
              __global const int      *level,      // 5  Array of level information.
              __global const int      *celltype,   // 6
              __global const real_t   *lev_dx,     // 7
              __global const real_t   *lev_dy,     // 8
              __global       real_t   *redscratch, // 9
              __global       real_t   *deltaT,     // 10
              __local        real_t   *tile)       // 11
      */

   real_t sigma_local = sigma;
   ezcl_set_kernel_arg(kernel_set_timestep,  0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_set_timestep,  1, sizeof(cl_real_t), (void *)&sigma_local);
   ezcl_set_kernel_arg(kernel_set_timestep,  2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_set_timestep,  3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_set_timestep,  4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_set_timestep,  5, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_set_timestep,  6, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_set_timestep,  7, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_set_timestep,  8, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_set_timestep,  9, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_set_timestep, 10, sizeof(cl_mem),  (void *)&dev_deltaT);
   ezcl_set_kernel_arg(kernel_set_timestep, 11, local_work_size*sizeof(cl_real_t),  NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_set_timestep, 1, NULL, &global_work_size, &local_work_size, NULL);

   if (block_size > 1){
         /*
         __kernel void finish_reduction_min_cl(
           const    int      isize,
           __global real_t  *redscratch,
           __global real_t  *deltaT,
           __local  real_t  *tile)
         */
      ezcl_set_kernel_arg(kernel_reduction_min, 0, sizeof(cl_int),  (void *)&block_size);
      ezcl_set_kernel_arg(kernel_reduction_min, 1, sizeof(cl_mem),  (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduction_min, 2, sizeof(cl_mem),  (void *)&dev_deltaT);
      ezcl_set_kernel_arg(kernel_reduction_min, 3, local_work_size*sizeof(cl_real_t), NULL);

     ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_min, 1, NULL, &local_work_size, &local_work_size, NULL);
   }

   real_t deltaT_local;
   ezcl_enqueue_read_buffer(command_queue, dev_deltaT, CL_TRUE,  0, sizeof(cl_real_t), &deltaT_local, NULL);
   deltaT = deltaT_local;

   globalmindeltaT = deltaT;
#ifdef HAVE_MPI
   if (parallel) MPI_Allreduce(&deltaT, &globalmindeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

   ezcl_device_memory_delete(dev_redscratch);

   gpu_time_set_timestep += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   return(globalmindeltaT);
}
#endif

void State::fill_circle(double  circ_radius,//  Radius of circle in grid units.
                        double  fill_value, //  Circle height for shallow water.
                        double  background) //  Background height for shallow water.
{  
   size_t &ncells = mesh->ncells;
   vector<spatial_t> &x  = mesh->x;
   vector<spatial_t> &dx = mesh->dx;
   vector<spatial_t> &y  = mesh->y;
   vector<spatial_t> &dy = mesh->dy;

   for (uint ic = 0; ic < ncells; ic++)
   {  H[ic] = background;
      U[ic] = V[ic] = 0.0; }
   
   //  Clear the old k-D tree and generate new data (slow but necessary here).
   //KDTree_Destroy(&mesh->tree);
   mesh->kdtree_setup();
   
   int nez;
   vector<int>    ind(ncells);
   vector<double> weight(ncells);
   
#ifdef FULL_PRECISION
   KDTree_QueryCircleInterior_Double(&mesh->tree, &nez, &(ind[0]), circ_radius, ncells,
                                     &x[0], &dx[0],
                                     &y[0], &dy[0]);
#else
   KDTree_QueryCircleInterior_Float(&mesh->tree, &nez, &(ind[0]), circ_radius, ncells,
                                    &x[0], &dx[0],
                                    &y[0], &dy[0]);
#endif
   for (int ic = 0; ic < nez; ++ic)
   {  H[ind[ic]] = fill_value; }
   
#ifdef FULL_PRECISION
   KDTree_QueryCircleIntersectWeighted_Double(&mesh->tree, &nez, &(ind[0]), &(weight[0]),
                              circ_radius, ncells,
                              &x[0], &dx[0],
                              &y[0], &dy[0]);
#else
   KDTree_QueryCircleIntersectWeighted_Float(&mesh->tree, &nez, &(ind[0]), &(weight[0]),
                              circ_radius, ncells,
                              &x[0], &dx[0],
                              &y[0], &dy[0]);
#endif

   for (int ic = 0; ic < nez; ++ic)
   {  H[ind[ic]] = background + (fill_value - background) * weight[ic]; }

   KDTree_Destroy(&mesh->tree);
}

void State::state_reorder(vector<int> iorder)
{
   H = state_memory.memory_reorder(H, &iorder[0]);
   U = state_memory.memory_reorder(U, &iorder[0]);
   V = state_memory.memory_reorder(V, &iorder[0]);
   //printf("\nDEBUG reorder cells\n"); 
   //state_memory.memory_report();
   //printf("DEBUG end reorder cells\n\n"); 
}

void State::rezone_all(int icount, int jcount, vector<int> mpot)
{
   mesh->rezone_all(icount, jcount, mpot, 1, state_memory);
   memory_reset_ptrs();
}


#ifdef HAVE_OPENCL
void State::gpu_rezone_all(int icount, int jcount, bool localStencil)
{
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- localStencil is %d\n",localStencil);

   mesh->gpu_rezone_all(icount, jcount, dev_mpot, gpu_state_memory);
   dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H");
   dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U");
   dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V");
}
#endif

//define macro for squaring a number
#define SQ(x) ((x)*(x))
//define macro to find minimum of 3 values
//#define MIN3(a,b,c) (min(min((a),(b)),(c)))

void State::calc_finite_difference(double deltaT){
	// Physical Constants
	real_t   g     = 9.80;   // gravitational constant
	real_t   ghalf = HALF*g;

	// Timers
	struct timeval tstart_cpu;
	cpu_timer_start(&tstart_cpu);

	// Grab the Physical Adaptive Mesh Cells
	size_t ncells           = mesh->ncells;
	size_t &ncells_ghost    = mesh->ncells_ghost;
	if (ncells_ghost < ncells) ncells_ghost = ncells;

	#ifdef HAVE_MPI
	// Populate the ghost regions since the calc neighbors has just been
	// established for the mesh shortly before
	if (mesh->numpe > 1) {
		apply_boundary_conditions_local();

        H=(state_t *)state_memory.memory_realloc(ncells_ghost, sizeof(state_t), H);
		U=(state_t *)state_memory.memory_realloc(ncells_ghost, sizeof(state_t), U);
		V=(state_t *)state_memory.memory_realloc(ncells_ghost, sizeof(state_t), V);

        L7_Update(&H[0], L7_STATE_T, mesh->cell_handle);
        L7_Update(&U[0], L7_STATE_T, mesh->cell_handle);
        L7_Update(&V[0], L7_STATE_T, mesh->cell_handle);

    apply_boundary_conditions_ghost();
    } else {
        apply_boundary_conditions();
    }
    #else
    apply_boundary_conditions();
    #endif

    int flags = 0;
    #if defined (HAVE_J7)
    if (mesh->parallel) flags = LOAD_BALANCE_MEMORY;
    #endif
    state_t *H_new = (state_t *) state_memory.memory_malloc
                                    (ncells_ghost, sizeof(state_t), flags, "H_new");
    state_t *U_new = (state_t *) state_memory.memory_malloc
                                    (ncells_ghost, sizeof(state_t), flags, "U_new");
    state_t *V_new = (state_t *) state_memory.memory_malloc
                                    (ncells_ghost, sizeof(state_t), flags, "V_new");
    // XXX Add state_memory.memory_calloc for convenience XXX
    int it;
    for (it = 0; it < ncells_ghost; it++)
        H_new[it] = U_new[it] = V_new[it] = 0.0;

    // Grab the Geometric Adaptive Mesh Cells
    // XXX Should this be done with the iteration? XXX
    int levmx  = mesh->levmx;
    int *level = mesh->level;
    int *nlft  = mesh->nlft;
    int *nrht  = mesh->nrht;
    int *nbot  = mesh->nbot;
    int *ntop  = mesh->ntop;
    // Allocate the Regular Grid
    // int nxmax = (mesh->xmax - mesh->xmin) / mesh->lev_deltax[levmx];
    // int nxmax = (mesh->imax + 1) * IPOW2(levmx)
    // int nymax = (mesh->ymax - mesh->ymin) / mesh->lev_deltay[levmx];
    // int nymax = (mesh->jmax + 1) * IPOW2(levmx)
    // Mesh(nx, ny, levmx, ndim, boundary, parallel, do_gpu_calc)
    // Mesh* rmesh = Mesh(nxmax, nymax, 0, 2, 0, 0, 0);
    #define COMPUTE (ONE)
    #define HALF_COMPUTE (HALF)
    #define DONT_COMPUTE (ZERO)
    #define RES (2) // Increase RESolution of the fine mesh by factor RES
    typedef enum {COMPUTE, HALF_COMPUTE, DONT_COMPUTE} computable_t;
    struct rMesh {
        int nxmax, nymax;
        int *gidx;
        real_t *controlVolume;
        computable_t *computeMe;
        // State information
        state_t *H, *U, *V;
        state_t *dHx, *dHy, *dU, *dV;
        rMesh(int nx, int ny, int ncells) {
            nxmax = RES*nx+(RES-1);
            nymax = RES*ny+(RES-1);
            gidx            = (int *)       malloc(ncells*sizeof(int));
            controlVolume   = (real_t *)    calloc(nx*ny*sizeof(real_t));
            computeMe       = (compute_t *) calloc(nx*ny*sizeof(compute_t));
            H               = (state_t *)   malloc(nx*ny*sizeof(state_t));
            U               = (state_t *)   malloc(nx*ny*sizeof(state_t));
            V               = (state_t *)   malloc(nx*ny*sizeof(state_t));
            dHx             = (state_t *)   malloc(nx*ny*sizeof(state_t));
            dHy             = (state_t *)   malloc(nx*ny*sizeof(state_t));
            dU              = (state_t *)   malloc(nx*ny*sizeof(state_t));
            dV              = (state_t *)   malloc(nx*ny*sizeof(state_t));
        }
    }

    rMesh *rmesh = rMesh( (mesh->imax+1)*IPOW2(levmx)
                        , (mesh->jmax+1)*IPOW2(levmx)
                        , ncells_ghost);

    // Iterate over the Regular Grid
    int lev;
    int acidx; // Adaptive Cell InDeX
    for ( lev = levmx; lev >= 0; --lev) {
        // XXX Counts here are not used at the moment XXX
        int n_cellComputable =      0;
        int n_half_cellComputable = 0;
        int n_faceComputable =      0;
        int n_half_faceComputable = 0;
        // IPOW2 a := (2 << a)
        int istride, jstride = IPOW2(levmx-lev) * RES;

        // Set the Regular Grid for the Current Level of Refinement
        for (acidx = 0; acidx < ncells_ghost; ++acidx) {
            if (level[acidx] == lev) {
                // Thread InDeX 
                // := There are cell centers at the current level ofthe mesh
                int tidx    = (mesh->j[acidx] * jstride + (RES-1)) * nxmax
                            + (mesh->i[acidx] * istride) + (RES-1);
                // Taking the controlVolumes to be discrete sums of the fine
                // regular grid cells ; the current version of CLAMR requires
                // a control volume at the current level to be HALF of the cell,
                // while a coarse cell which is refined to the current level
                // takes the FULL cell -- this will likely change in the future
                gix->[tidx] = acidx;
                rmesh->controlVolume[tidx] = HALF*(istride*jstide);
                rmesh->computeMe[tidx] = COMPUTE;
                ++n_cellComputable;
                rmesh->H[tidx] = H[acidx];
                rmesh->U[tidx] = U[acidx];
                rmesh->V[tidx] = V[acidx];

                int tnl   = tidx - istride; // Thread Neighbor Left
                int acnl  = nlft[acidx];    // Adaptive Cell Neigbor Left
                int tnll  = tidx - 2*istride;   // " " Left Left
                int acnll = nlft[acnl];  // " " " Left Left
                if (level[acnl] < lev) {
                    rmesh->controlVolume[tnl] = istride*jstride;
                    rmesh->computeMe[tnl] = HALF_COMPUTE;
                    ++n_half_cellComputable;
                    rmesh->H[tnl] = H[acnl];
                    rmesh->U[tnl] = U[acnl];
                    rmesh->V[tnl] = V[acnl];
                    rmesh->computeMe[tidx-istride/RES] = HALF_COMPUTE;
                    ++n_half_faceComputable;
                    rmesh->H[tidx-istride/RES] = (H[acnl]+2.0*rmesh->H[tidx])/3.0;
                    rmesh->U[tidx-istride/RES] = (U[acnl]+2.0*rmesh->U[tidx])/3.0;
                    rmesh->V[tidx-istride/RES] = (V[acnl]+2.0*rmesh->V[tidx])/3.0;
                    if ( level[acnll] == lev ) {
                        int acnllt = ntop[acnll];
                        rmesh->dHx[tnll] = rmesh->H[tnl]-HALF*(H[acnll]+H[acnllt]);
                        rmesh->dU[tnll]  = rmesh->U[tnl]-HALF*(U[acnll]+U[acnllt]);
                    } else {
                        rmesh->dHx[tnll] = rmesh->H[tnl] - H[acnll];
                        rmesh->dU[tnll]  = rmesh->U[tnl] - U[acnll];
                    }
                // Don't need to set to DONT_COMPUTE b/c already initialized
                // } else if (level[acnl] > lev) {
                //    rmesh->computeMe[tnl] = DONT_COMPUTE;
                //    rmesh->computeMe[tidx-istride/RES] = DONT_COMPUTE;
                } else if (level[acnl] == lev) {
                    rmesh->computeMe[tidx-istride/RES] = COMPUTE;
                    rmesh->H[tidx-istride/RES] = (H[acnl]+rmesh->H[tidx])/2.0;
                    rmesh->U[tidx-istride/RES] = (U[acnl]+rmesh->U[tidx])/2.0;
                    rmesh->V[tidx-istride/RES] = (V[acnl]+rmesh->V[tidx])/2.0;
                    if(level[acnll] > lev) {
                        rmesh->dHx[tnll] = H[acnl] - HALF*(H[acnll]+H[ntop[acnll]]);
                        rmesh->dU[tnll]  = U[acnl] - HALF*(U[acnll]+U[ntop[acnll]]);
                    } else {
                        rmesh->dHx[tnll] = H[acnl] - H[acnll];
                        rmesh->dU[tnll]  = U[acnl] - U[acnll];
                    }
                }

                int tnr   = tidx + istride; // Thread Neighbor Right
                int acnr  = nrht[acidx];    // Adaptive Cell Neigbor Right
                int tnrr  = tidx + 2*istride;   // " " Right Right
                int acnrr = nrht[acnr];  // " " " Right Right
                if (level[acnr] < lev) {
                    rmesh->controlVolume[tnr] = istride*jstride;
                    rmesh->computeMe[tnr] = HALF_COMPUTE;
                    ++n_half_cellComputable;
                    rmesh->H[tnr] = H[acnr];
                    rmesh->U[tnr] = U[acnr];
                    rmesh->V[tnr] = V[acnr];
                    rmesh->computeMe[tidx+istride/RES] = HALF_COMPUTE;
                    ++n_half_faceComputable;
                    rmesh->H[tidx+istride/RES] = (H[acnr]+2.0*rmesh->H[tidx])/3.0;
                    rmesh->U[tidx+istride/RES] = (U[acnr]+2.0*rmesh->U[tidx])/3.0;
                    rmesh->V[tidx+istride/RES] = (V[acnr]+2.0*rmesh->V[tidx])/3.0;
                    if ( level[acnrr] == lev ) {
                        int acnrrt = ntop[acnrr];
                        rmesh->dHx[tnrr] = HALF*(H[acnrr]+H[acnrrt])-rmesh->H[tnr];
                        rmesh->dU[tnrr]  = HALF*(U[acnrr]+U[acnrrt])-rmesh->U[tnr];
                    } else {
                        rmesh->dHx[tnrr] = H[acnrr] - rmesh->H[tnr];
                        rmesh->dU[tnrr]  = U[acnrr] - rmesh->U[tnr];
                    }
                // Don't need to set to DONT_COMPUTE b/c already initialized
                // } else if (level[acnr] > lev) {
                //    rmesh->computeMe[tnr] = DONT_COMPUTE;
                //    rmesh->computeMe[tidx+istride/RES] = DONT_COMPUTE;
                } else if (level[acnr] == lev) {
                // Don't need to set the face b/c already set by left pass
                //    rmesh->computeMe[tidx+istride/RES] = COMPUTE;
                //    rmesh->H[tidx+istride/RES] = (H[acnr]+rmesh->H[tidx])/2.0;
                //    rmesh->U[tidx+istride/RES] = (U[acnr]+rmesh->U[tidx])/2.0;
                //    rmesh->V[tidx+istride/RES] = (V[acnr]+rmesh->V[tidx])/2.0;
                    if(level[acnrr] > lev) {
                        rmesh->dHx[tnrr] = HALF*(H[acnrr]+H[ntop[acnrr]])-H[acnr];
                        rmesh->dU[tnrr]  = HALF*(U[acnrr]+U[ntop[acnrr]])-U[acnr];
                    } else {
                        rmesh->dHx[tnll] = H[acnrr] - H[acnr];
                        rmesh->dU[tnll]  = U[acnrr] - U[acnr];
                    }
                }

                int tnb   = tidx - jstride * nxmax; // Thread Neighbor Bottom
                int acnb  = nbot[acidx];    // Adaptive Cell Neigbor Bottom
                int tnbb  = tidx - 2 * jstride * nxmax;   // " " Bottom Bottom
                int acnbb = nbot[acnb];  // " " " Bottom Bottom
                if (level[acnb] < lev) {
                    rmesh->controlVolume[tnb] = istride*jstride;
                    rmesh->computeMe[tnb] = HALF_COMPUTE;
                    ++n_half_cellComputable;
                    rmesh->H[tnb] = H[acnb];
                    rmesh->U[tnb] = U[acnb];
                    rmesh->V[tnb] = V[acnb];
                    int tnbface = tidx - (jstride/RES) * nxmax;
                    rmesh->computeMe[tnbface] = HALF_COMPUTE;
                    ++n_half_faceComputable;
                    rmesh->H[tnbface] = (H[acnb]+2.0*rmesh->H[tidx])/3.0;
                    rmesh->U[tnbface] = (U[acnb]+2.0*rmesh->U[tidx])/3.0;
                    rmesh->V[tnbface] = (V[acnb]+2.0*rmesh->V[tidx])/3.0;
                    if ( level[acnbb] == lev ) {
                        int acnbbr = nrht[acnbb];
                        rmesh->dHy[tnbb] = rmesh->H[tnb]-HALF*(H[acnbb]+H[acnbbr]);
                        rmesh->dV[tnbb]  = rmesh->V[tnb]-HALF*(V[acnbb]+V[acnbbr]);
                    } else {
                        rmesh->dHy[tnbb] = rmesh->H[tnb] - H[acnbb];
                        rmesh->dV[tnbb]  = rmesh->V[tnb] - V[acnbb];
                    }
                // Don't need to set to DONT_COMPUTE b/c already initialized
                // } else if (level[acnb] > lev) {
                //    rmesh->computeMe[tnb] = DONT_COMPUTE;
                //    rmesh->computeMe[tidx-(jstride/RES)*nxmax] = DONT_COMPUTE;
                } else if (level[acnb] == lev) {
                    int tnbface = tidx - (jstride/RES) * nxmax;
                    rmesh->computeMe[tnbface] = COMPUTE;
                    rmesh->H[tnbface] = (H[acnb]+rmesh->H[tidx])/2.0;
                    rmesh->U[tnbface] = (U[acnb]+rmesh->U[tidx])/2.0;
                    rmesh->V[tnbface] = (V[acnb]+rmesh->V[tidx])/2.0;
                    if(level[acnbb] > lev) {
                        rmesh->dHy[tnbb] = H[acnb] - HALF*(H[acnbb]+H[nrht[acnbb]]);
                        rmesh->dV[tnbb]  = V[acnb] - HALF*(V[acnbb]+V[nrht[acnbb]]);
                    } else {
                        rmesh->dHy[tnbb] = H[acnb] - H[acnbb];
                        rmesh->dV[tnbb]  = V[acnb] - V[acnbb];
                    }
                }

                int tnt   = tidx + jstride * nxmax; // Thread Neighbor Top
                int acnt  = ntop[acidx];    // Adaptive Cell Neigbor Top
                int tntt  = tidx + 2 * jstride * nxmax;   // " " Top Top
                int acntt = ntop[acnt];  // " " " Top Top
                if (level[acnt] < lev) {
                    rmesh->controlVolume[tnt] = istride*jstride;
                    rmesh->computeMe[tnt] = HALF_COMPUTE;
                    ++n_half_cellComputable;
                    rmesh->H[tnt] = H[acnt];
                    rmesh->U[tnt] = U[acnt];
                    rmesh->V[tnt] = V[acnt];
                    int tntface = tidx + (jstride/RES) * nxmax;
                    rmesh->computeMe[tntface] = HALF_COMPUTE;
                    ++n_half_faceComputable;
                    rmesh->H[tntface] = (H[acnt]+2.0*rmesh->H[tidx])/3.0;
                    rmesh->U[tntface] = (U[acnt]+2.0*rmesh->U[tidx])/3.0;
                    rmesh->V[tntface] = (V[acnt]+2.0*rmesh->V[tidx])/3.0;
                    if ( level[acntt] == lev ) {
                        int acnttr = nrht[acntt];
                        rmesh->dHy[tntt] = HALF*(H[acntt]+H[acnttr])-rmesh->H[tnt];
                        rmesh->dV[tntt]  = HALF*(V[acntt]+V[acnttr])-rmesh->V[tnt];
                    } else {
                        rmesh->dHy[tntt] = H[acntt] - rmesh->H[tnt];
                        rmesh->dV[tntt]  = V[acntt] - rmesh->V[tnt];
                    }
                // Don't need to set to DONT_COMPUTE b/c already initialized
                // } else if (level[acnt] > lev) {
                //    rmesh->computeMe[tnt] = DONT_COMPUTE;
                //    rmesh->computeMe[tidx+(jstride/RES)*nxmax] = DONT_COMPUTE;
                } else if (level[acnt] == lev) {
                // Don't need to set the face b/c already set by bottom pass
                //    int tntface = tidx + (jstride/RES) * nxmax;
                //    rmesh->computeMe[tntface] = COMPUTE;
                //    rmesh->H[tntface] = (H[acnt]+rmesh->H[tidx])/2.0;
                //    rmesh->U[tntface] = (U[acnt]+rmesh->U[tidx])/2.0;
                //    rmesh->V[tntface] = (V[acnt]+rmesh->V[tidx])/2.0;
                    if(level[acntt] > lev) {
                        rmesh->dHy[tntt] = HALF*(H[acntt]+H[nrht[acntt]])-H[acnt];
                        rmesh->dV[tntt]  = HALF*(V[acntt]+V[nrht[acntt]])-V[acnt];
                    } else {
                        rmesh->dHy[tntt] = H[acntt] - H[acnt];
                        rmesh->dV[tntt]  = V[acntt] - V[acnt];
                    }
                }
            }
        }

        // Calculate the Finite Difference on the Regular Grid
        // int gidx;
        // for(gidx = 0; gidx < (n_faceComputable+n_half_faceComputable); ++gidx) {
        int jgidx, igidx;
        for(jgidx = (RES - 1); jgidx < nymax; ++jstride) {
            for(igidx = (RES - 1); igidx < nxmax; ++istride) {
                int gidx = 
                int tidx = jgidx * nxmax + igidx;
                // XXX XXX XXX
                // Checking computability of a cell MUST,
                // MUST be hidden from the state kernel
                // XXX XXX XXX
                if(rmesh->computeMe[tidx] != DONT_COMPUTE) {

                    // Half time step calculation:

                    inline real_t U_halfT(
                        real_t  deltaT, // Timestep
                        real_t  U,      // Face state variable
                        real_t  F_i,    // Downwind state variable flux
                        real_t  F_n,    // Upwind   state variable flux
                        real_t  A,      // Control V flux area, set by lev
                        real_t  V) {    // Control V, downwind + upwind control V

                        return U - HALF * deltaT * ( A * ( F_n - F_i ) / V );
                    }

                    #define Hic ( rmesh->H[tidx] )
                    #define Hl  ( rmesh->H[tidx-istride] )
                    #define Hr  ( rmesh->H[tidx+istride] )
                    #define Hb  ( rmesh->H[tidx-jstride*nxmax] )
                    #define Ht  ( rmesh->H[tidx+jstride*nxmax] )
                    #define Uic ( rmesh->U[tidx] )
                    #define Ul  ( rmesh->U[tidx-istride] )
                    #define Ur  ( rmesh->U[tidx+istride] )
                    #define Ub  ( rmesh->U[tidx-jstride*nxmax] )
                    #define Ut  ( rmesh->U[tidx+jstride*nxmax] )
                    #define Vic ( rmesh->V[tidx] )
                    #define Vl  ( rmesh->V[tidx-istride] )
                    #define Vr  ( rmesh->V[tidx+istride] )
                    #define Vb  ( rmesh->V[tidx-jstride*nxmax] )
                    #define Vt  ( rmesh->V[tidx+jstride*nxmax] )
                    
                    #define Hlface ( rmesh->H[tidx-istride/RES] )
                    #define Hrface ( rmesh->H[tidx+istride/RES] )
                    #define Hbface ( rmesh->H[tidx-(jstride/RES)*nxmax] )
                    #define Htface ( rmesh->H[tidx+(jstride/RES)*nxmax] )
                    #define Ulface ( rmesh->U[tidx-istride/RES] )
                    #define Urface ( rmesh->U[tidx+istride/RES] )
                    #define Ubface ( rmesh->U[tidx-(jstride/RES)*nxmax] )
                    #define Utface ( rmesh->U[tidx+(jstride/RES)*nxmax] )
                    #define Vlface ( rmesh->V[tidx-istride/RES] )
                    #define Vrface ( rmesh->V[tidx+istride/RES] )
                    #define Vbface ( rmesh->V[tidx-(jstride/RES)*nxmax] )
                    #define Vtface ( rmesh->V[tidx+(jstride/RES)*nxmax] )
                    
                    #define dHxic   ( rmesh->dHx[tidx] )
                    #define dHyic   ( rmesh->dHy[tidx] )
                    #define dUic    ( rmesh->dU[tidx] )
                    #define dVic    ( rmesh->dV[tidx] )
                    
                    #define dHl  ( rmesh->dHx[tidx-istride] )
                    #define dHr  ( rmesh->dHx[tidx+istride] )
                    #define dHb  ( rmesh->dHy[tidx-(jstride*nxmax)] )
                    #define dHt  ( rmesh->dHy[tidx+(jstride*nxmax)] )
                    #define dUl  ( rmesh->dU[tidx-istride] )
                    #define dUr  ( rmesh->dU[tidx+istride] )
                    #define dUb  ( rmesh->dU[tidx-(jstride*nxmax)] )
                    #define dUt  ( rmesh->dU[tidx+(jstride*nxmax)] )
                    #define dVl  ( rmesh->dV[tidx-istride] )
                    #define dVr  ( rmesh->dV[tidx+istride] )
                    #define dVb  ( rmesh->dV[tidx-(jstride*nxmax)] )
                    #define dVt  ( rmesh->dV[tidx+(jstride*nxmax)] )
                    
                    #define dHll  ( rmesh->dHx[tidx-2*istride] )
                    #define dHrr  ( rmesh->dHx[tidx+2*istride] )
                    #define dHbb  ( rmesh->dHy[tidx-2*(jstride*nxmax)] )
                    #define dHtt  ( rmesh->dHy[tidx+2*(jstride*nxmax)] )
                    #define dUll  ( rmesh->dU[tidx-2*istride] )
                    #define dUrr  ( rmesh->dU[tidx+2*istride] )
                    #define dUbb  ( rmesh->dU[tidx-2*(jstride*nxmax)] )
                    #define dUtt  ( rmesh->dU[tidx+2*(jstride*nxmax)] )
                    #define dVll  ( rmesh->dV[tidx-2*istride] )
                    #define dVrr  ( rmesh->dV[tidx+2*istride] )
                    #define dVbb  ( rmesh->dV[tidx-2*(jstride*nxmax)] )
                    #define dVtt  ( rmesh->dV[tidx+2*(jstride*nxmax)] )
                    
                    #define cVol(idx) ( rmesh->controlVolume[idx] )
                    
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

                    double dT   = deltaT;

                    real_t A    = (real_t) jstride;
                    real_t V    = (real_t) (cVol(tidx) + cVol(tidx-istride));
                    real_t Hxminus = U_halfT(dT, Hlface, HXFLUXNL, HXFLUXIC, A, V);
                    real_t Uxminus = U_halfT(dT, Ulface, UXFLUXNL, UXFLUXIC, A, V);
                    real_t Vxminus = U_halfT(dT, Vlface, UVFLUXNL, UVFLUXIC, A, V);

                           A    = (real_t) jstride;
                           V    = (real_t) (cVol(tidx) + cVol(tidx+istride));
                    real_t Hxplus  = U_halfT(dT, Hrface, HXFLUXIC, HXFLUXNR, A, V);
                    real_t Uxplus  = U_halfT(dT, Urface, UXFLUXIC, UXFLUXNR, A, V);
                    real_t Vxplus  = U_halfT(dT, Vrface, UVFLUXIC, UVFLUXNR, A, V);

                           A    = (real_t) istride;
                           V    = (real_t) (cVol(tidx) + cVol(tidx-jstride*nxmax));
                    real_t Hyminus = U_halfT(dT, Hbface, HYFLUXNB, HYFLUXIC, A, V);
                    real_t Uyminus = U_halfT(dT, Ubface, VUFLUXNB, VUFLUXIC, A, V);
                    real_t Vyminus = U_halfT(dT, Vbface, VYFLUXNB, VYFLUXIC, A, V);

                           A    = (real_t) istride;
                           V    = (real_t) (cVol(tidx) + cVol(tidx+jstride*nxmax));
                    real_t Hyplus  = U_halfT(dT, Htface, HYFLUXIC, HYFLUXNT, A, V);
                    real_t Uyplus  = U_halfT(dT, Utface, VUFLUXIC, VUFLUXNT, A, V);
                    real_t Vyplus  = U_halfT(dT, Vtface, VYFLUXIC, VYFLUXNT, A, V);

                    // Artificial Viscosity corrections
                    inline real_t w_corrector(
                        real_t  deltaT,         // Timestep
                        real_t  dr,             // Cell dimension
                        real_t  U_eigen,        // State variable eigenvalue (speed)
                        real_t  grad_half,      // Centered gradient
                        real_t  grad_minus,     // Downwind gradient
                        real_t  grad_plus) {    // Upwind gradient

                        real_t  nu = HALF * U_eigen * deltaT / dr;
                                nu = nu * (ONE - nu);

                        real_t  rdenom = ONE / max(SQR(grad_half), EPSILON);
                        real_t  rplus  = (grad_plus  * grad_half) * rdenom;
                        real_t  rminus = (grad_minus * grad_half) * rdenom;

                        return  HALF*nu*(ONE - max(MIN3(ONE, rplus, rminus), ZERO));
                    }

                    real_t  wminusx_H = w_corrector(deltaT, istride
                                        , fabs(Uxminus/Hxminus) + sqrt(g*Hxminus)
                                        , dHxic, dHl, dHr);
                            wminusx_H *= dHxic;

                    real_t  wplusx_H  = w_corrector(deltaT, istride
                                        , fabs(Uxplus/Hxplus) + sqrt(g*Hxplus)
                                        , dHr, dxHic, dHrr);
                            wplusx_H  *= dHr;

                    real_t  wminusx_U = w_corrector(deltaT, istride
                                        , fabs(Uxminus/Hxminus) + sqrt(g*Hxminus)
                                        , dUic, dUl, dUr);
                            wminusx_U *= dUic;

                    real_t  wplusx_U  = w_corrector(deltaT, istride
                                        , fabs(Uxplus/Hxplus) + sqrt(g*Hxplus)
                                        , dUr, dUic, dUrr);
                            wplusx_U  *= dUr;

                    real_t  wminusy_H = w_corrector(deltaT, jstride
                                        , fabs(Vyminus/Hyminus) + sqrt(g*Hyminus)
                                        , dHyic, dHb, dHt);
                            wminusy_H *= dHyic;

                    real_t  wplusy_H  = w_corrector(deltaT, jstride
                                        , fabs(Vyplus/Hyplus) + sqrt(g*Hyplus)
                                        , dHt, dHyic, dHtt);
                            wplusy_H  *= dHt;

                    real_t  wminusy_V = w_corrector(deltaT, jstride
                                        , fabs(Vyminus/Hyminus) + sqrt(g*Hyminus)
                                        , dVic, dVb, dVt);
                            wminusy_V *= dVic;

                    real_t  wplusy_V  = w_corrector(deltaT, jstride
                                        , fabs(Vyplus/Hyplus) + sqrt(g*Hyplus)
                                        , dVt, dVic, dVtt);
                            wplusy_V  *= dVt;

                    // Full time step calculation:
                    
                    inline real_t U_fullT(
                        real_t    deltaT,
                        real_t    dr,
                        real_t    F_plus,
                        real_t    F_minus,
                        real_t    G_plus,
                        real_t    G_minus) {

                        return (-(deltaT/dr) * (F_plus-F_minus+G_plus-G_minus));
                    }

                    // XXX NEED TO HANDLE HALF_COMPUTE vs. FULL_COMPUTE HERE XXX
                    real_t scalar = (real_t)rmesh->computeMe[tidx];
                    int gix = rmesh->gidx[tidx];
                    H_new[gix]  += scalar * U_fullT(dT, istride
                                , Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                                - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
                    U_new[gix]  += scalar * U_fullT(dT, istride
                                , Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                                - wminusx_U + wplusx_U;
                    V_new[gix]  += scalar * U_fullT(dT, jstride
                                , Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                                - wminusy_V + wplusy_V;
                }   
            }
        }
    }

    free(rmesh);
    // Replace H with H_new and deallocate H. New memory has the characteristics of
    // the new memory but the name of the old.
    H = (state_t *)state_memory.memory_replace(H, H_new);
    U = (state_t *)state_memory.memory_replace(U, U_new);
    V = (state_t *)state_memory.memory_replace(V, V_new);

    cpu_time_finite_difference += cpu_timer_stop(tstart_cpu);
}

#ifdef HAVE_OPENCL
void State::gpu_calc_finite_difference(double deltaT)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_command_queue command_queue = ezcl_get_command_queue();

   //cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;
   int &levmx           = mesh->levmx;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
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
   assert(dev_levdx);
   assert(dev_levdy);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_H_new"));
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_U_new"));
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_V_new"));
 
   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_H_new"));
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_U_new"));
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_V_new"));

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

     /*
     __kernel void calc_finite_difference_cl(
                      const int     ncells,    // 0  Total number of cells.
                      const int     lvmax,     // 1  Maximum level
             __global       state_t *H,        // 2
             __global       state_t *U,        // 3
             __global       state_t *V,        // 4
             __global       state_t *H_new,    // 5
             __global       state_t *U_new,    // 6
             __global       state_t *V_new,    // 7
             __global const int     *nlft,     // 8  Array of left neighbors.
             __global const int     *nrht,     // 9  Array of right neighbors.
             __global const int     *ntop,     // 10  Array of bottom neighbors.
             __global const int     *nbot,     // 11  Array of top neighbors.
             __global const int     *level,    // 12  Array of level information.
                      const real_t   deltaT,   // 13  Size of time step.
             __global const real_t  *lev_dx,   // 14
             __global const real_t  *lev_dy,   // 15
             __local        state4_t *tile,    // 16  Tile size in state4.
             __local        int8  *itile)      // 17  Tile size in int8.
     */
   cl_event calc_finite_difference_event;

   real_t deltaT_local = deltaT;
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
   ezcl_set_kernel_arg(kernel_calc_finite_difference,13, sizeof(cl_real_t), (void *)&deltaT_local);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,14, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,15, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,16, local_work_size*sizeof(cl_state4_t),    NULL);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,17, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference,   1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_event);

   ezcl_wait_for_events(1, &calc_finite_difference_event);
   ezcl_event_release(calc_finite_difference_event);

   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

   gpu_time_finite_difference += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
}
#endif

void State::symmetry_check(const char *string, vector<int> sym_index, double eps,
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

size_t State::calc_refine_potential(vector<int> &mpot,int &icount, int &jcount)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timeval tstart_lev2;
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

   size_t ncells     = mesh->ncells;
   int *nlft  = mesh->nlft;
   int *nrht  = mesh->nrht;
   int *nbot  = mesh->nbot;
   int *ntop  = mesh->ntop;
   int *level = mesh->level;

   icount=0;
   jcount=0;

#ifdef HAVE_MPI
   // We need to update the ghost regions and boundary regions for the state
   // variables since they were changed in the finite difference routine. We
   // want to use the updated values for refinement decisions
   if (mesh->numpe > 1) {
      apply_boundary_conditions_local();

      L7_Update(&H[0], L7_STATE_T, mesh->cell_handle);
      L7_Update(&U[0], L7_STATE_T, mesh->cell_handle);
      L7_Update(&V[0], L7_STATE_T, mesh->cell_handle);

      apply_boundary_conditions_ghost();
   } else {
      apply_boundary_conditions();
   }
#else
   apply_boundary_conditions();
#endif

   int ic;
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (ic=0; ic<(int)ncells; ic++) {

      if (mesh->celltype[ic] != REAL_CELL) continue;

      state_t Hic = H[ic];
      //state_t Uic = U[ic];
      //state_t Vic = V[ic];

      int nl = nlft[ic];
      state_t Hl = H[nl];
      //state_t Ul = U[nl];
      //state_t Vl = V[nl];

      if (level[nl] > level[ic]){
         int nlt = ntop[nl];
         Hl = REFINE_HALF * (Hl + H[nlt]);
      }

      int nr = nrht[ic];
      state_t Hr = H[nr];
      //state_t Ur = U[nr];
      //state_t Vr = V[nr];

      if (level[nr] > level[ic]){
         int nrt = ntop[nr];
         Hr = REFINE_HALF * (Hr + H[nrt]);
      }

      int nb = nbot[ic];
      state_t Hb = H[nb];
      //state_t Ub = U[nb];
      //state_t Vb = V[nb];

      if (level[nb] > level[ic]){
         int nbr = nrht[nb];
         Hb = REFINE_HALF * (Hb + H[nbr]);
      }

      int nt = ntop[ic];
      state_t Ht = H[nt];
      //state_t Ut = U[nt];
      //state_t Vt = V[nt];

      if (level[nt] > level[ic]){
         int ntr = nrht[nt];
         Ht = REFINE_HALF * (Ht + H[ntr]);
      }

      state_t duplus1; //, duplus2;
      state_t duhalf1; //, duhalf2;
      state_t duminus1; //, duminus2;

      duplus1 = Hr-Hic;
      //duplus2 = Ur-Uic;
      duhalf1 = Hic-Hl;
      //duhalf2 = Uic-Ul;

      state_t qmax = REFINE_NEG_THOUSAND;

      state_t qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hl;
      //duminus2 = Uic-Ul;
      duhalf1 = Hr-Hic;
      //duhalf2 = Ur-Uic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duplus1 = Ht-Hic;
      //duplus2 = Vt-Vic;
      duhalf1 = Hic-Hb;
      //duhalf2 = Vic-Vb;

      qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hb;
      //duminus2 = Vic-Vb;
      duhalf1 = Ht-Hic;
      //duhalf2 = Vt-Vic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      mpot[ic]=0;
      if (qmax > REFINE_GRADIENT && level[ic] < mesh->levmx) {
         mpot[ic]=1;
      } else if (qmax < COARSEN_GRADIENT && level[ic] > 0) {
         mpot[ic] = -1;
      }
      //if (mpot[ic]) printf("DEBUG cpu cell is %d mpot %d\n",ic,mpot[ic]);
   }

   if (TIMING_LEVEL >= 2) {
      cpu_time_calc_mpot += cpu_timer_stop(tstart_lev2);
   }

   int newcount = mesh->refine_smooth(mpot, icount, jcount);
   //printf("DEBUG -- after refine smooth in file %s line %d icount %d jcount %d newcount %d\n",__FILE__,__LINE__,icount,jcount,newcount);

   cpu_time_refine_potential += cpu_timer_stop(tstart_cpu);

   return(newcount);
}

#ifdef HAVE_OPENCL
size_t State::gpu_calc_refine_potential(int &icount, int &jcount)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timeval tstart_lev2;
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t &ncells       = mesh->ncells;
   int &levmx           = mesh->levmx;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   //cl_mem &dev_mpot     = mesh->dev_mpot;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
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
   assert(dev_i);
   assert(dev_j);
   assert(dev_level);
   //assert(dev_mpot);
   //assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   icount = 0;
   jcount = 0;

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size = global_work_size/local_work_size;

#ifdef HAVE_MPI
   //size_t nghost_local = mesh->ncells_ghost - ncells;

   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

#ifdef BOUNDS_CHECK
      {
         vector<int> nlft_tmp(mesh->ncells_ghost);
         vector<int> nrht_tmp(mesh->ncells_ghost);
         vector<int> nbot_tmp(mesh->ncells_ghost);
         vector<int> ntop_tmp(mesh->ncells_ghost);
         vector<int> level_tmp(mesh->ncells_ghost);
         vector<state_t> H_tmp(mesh->ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, mesh->ncells_ghost*sizeof(cl_int), &nlft_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, mesh->ncells_ghost*sizeof(cl_int), &nrht_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, mesh->ncells_ghost*sizeof(cl_int), &nbot_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_TRUE,  0, mesh->ncells_ghost*sizeof(cl_int), &ntop_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, mesh->ncells_ghost*sizeof(cl_int), &level_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE,  0, mesh->ncells_ghost*sizeof(cl_int), &H_tmp[0],     NULL);
         for (uint ic=0; ic<ncells; ic++){
            int nl = nlft_tmp[ic];
            if (nl<0 || nl>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nlft %d\n",mesh->mype,__LINE__,ic,nl);
            if (level_tmp[nl] > level_tmp[ic]){
               int ntl = ntop_tmp[nl];
               if (ntl<0 || ntl>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d global %d nlft %d ntop of nlft %d\n",mesh->mype,__LINE__,ic,ic+mesh->noffset,nl,ntl);
            }
            int nr = nrht_tmp[ic];
            if (nr<0 || nr>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nrht %d\n",mesh->mype,__LINE__,ic,nr);
            if (level_tmp[nr] > level_tmp[ic]){
               int ntr = ntop_tmp[nr];
               if (ntr<0 || ntr>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d ntop of nrht %d\n",mesh->mype,__LINE__,ic,ntr);
            }
            int nb = nbot_tmp[ic];
            if (nb<0 || nb>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nbot %d\n",mesh->mype,__LINE__,ic,nb);
            if (level_tmp[nb] > level_tmp[ic]){
               int nrb = nrht_tmp[nb];
               if (nrb<0 || nrb>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nrht of nbot %d\n",mesh->mype,__LINE__,ic,nrb);
            }
            int nt = ntop_tmp[ic];
            if (nt<0 || nt>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d ntop %d\n",mesh->mype,__LINE__,ic,nt);
            if (level_tmp[nt] > level_tmp[ic]){
               int nrt = nrht_tmp[nt];
               if (nrt<0 || nrt>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nrht of ntop %d\n",mesh->mype,__LINE__,ic,nrt);
            }
         }
         for (uint ic=0; ic<mesh->ncells_ghost; ic++){
            if (H_tmp[ic] < 1.0) printf("%d: Warning at line %d cell %d H %lf\n",mesh->mype,__LINE__,ic,H_tmp[ic]);
         }
      }
#endif

   size_t result_size = 1;
   cl_mem dev_result     = ezcl_malloc(NULL, const_cast<char *>("dev_result"),     &result_size,        sizeof(cl_int2), CL_MEM_READ_WRITE, 0);
   cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size,         sizeof(cl_int2), CL_MEM_READ_WRITE, 0);

   dev_mpot              = ezcl_malloc(NULL, const_cast<char *>("dev_mpot"),       &mesh->ncells_ghost, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

     /*
     __kernel void refine_potential
              const int      ncells,     // 0  Total number of cells.
              const int      levmx,      // 1  Maximum level
     __global       state_t *H,          // 2
     __global       state_t *U,          // 3
     __global       state_t *V,          // 4
     __global const int     *nlft,       // 5  Array of left neighbors.
     __global const int     *nrht,       // 6  Array of right neighbors.
     __global const int     *ntop,       // 7  Array of bottom neighbors.
     __global const int     *nbot,       // 8  Array of top neighbors.
     __global const int     *level,      // 9  Array of level information.
     __global const int     *celltype,   // 10  Array of celltype information.
     __global       int     *mpot,       // 11  Array of mesh potential information.
     __global       int2    *redscratch, // 12
     __global const real_t  *lev_dx,     // 13
     __global const real_t  *lev_dy,     // 14
     __global       int2    *result,     // 15
     __local        state_t *tile,       // 16  Tile size in real4.
     __local        int8    *itile)      // 17  Tile size in int8.
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
   ezcl_set_kernel_arg(kernel_refine_potential, 9, sizeof(cl_mem),  (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_refine_potential,10, sizeof(cl_mem),  (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_refine_potential,11, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_refine_potential,12, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_refine_potential,13, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_refine_potential,14, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_refine_potential,15, sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_refine_potential,16, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_refine_potential,17, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_refine_potential,18, local_work_size*sizeof(cl_state_t),    NULL);
   ezcl_set_kernel_arg(kernel_refine_potential,19, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_refine_potential, 1, NULL, &global_work_size, &local_work_size, NULL);

   mesh->gpu_rezone_count2(block_size, local_work_size, dev_redscratch, dev_result);

   int count[2] = {0, 0};
   ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int2), count, NULL);
   icount  = count[0];
   jcount  = count[1];
   //size_t result = ncells + icount - jcount;

   //int mpot_check[ncells];
   //ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE, 0, ncells*sizeof(cl_int), mpot_check, NULL);
   //for (int ic=0; ic<ncells; ic++){
   //   if (mpot_check[ic]) printf("DEBUG -- cell %d mpot %d\n",ic,mpot_check[ic]);
   //}

   //printf("result = %lu after first refine potential icount %d jcount %d\n",result, icount, jcount);
//   int which_smooth = 1;

   ezcl_device_memory_delete(dev_redscratch);
   ezcl_device_memory_delete(dev_result);

   if (TIMING_LEVEL >= 2) {
      gpu_time_calc_mpot += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
   }

   int my_result = mesh->gpu_refine_smooth(dev_mpot, icount, jcount);
   //printf("DEBUG gpu calc refine potential %d icount %d jcount %d\n",my_result,icount,jcount);

   gpu_time_refine_potential += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   return((size_t)my_result);
}
#endif

double State::mass_sum(int enhanced_precision_sum)
{
   size_t &ncells = mesh->ncells;
   int *celltype = mesh->celltype;
   int *level    = mesh->level;

#ifdef HAVE_MPI
   //int &mype = mesh->mype;
#endif

   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   double summer = 0.0;
   double total_sum = 0.0;

   if (enhanced_precision_sum == SUM_KAHAN) {
      //printf("DEBUG -- kahan_sum\n");
      double corrected_next_term, new_sum;
      struct esum_type local;
#ifdef HAVE_MPI
      struct esum_type global;
#endif

      local.sum = 0.0;
      local.correction = 0.0;
      int ic;
      for (ic = 0; ic < (int)ncells; ic++) {
         if (celltype[ic] == REAL_CELL) {
            //  Exclude boundary cells.
            corrected_next_term= H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]] + local.correction;
            new_sum            = local.sum + local.correction;
            local.correction   = corrected_next_term - (new_sum - local.sum);
            local.sum          = new_sum;
         }
      }

#ifdef HAVE_MPI
      if (mesh->parallel) {
         MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KAHAN_SUM, MPI_COMM_WORLD);
         total_sum = global.sum + global.correction;
      } else {
         total_sum = local.sum + local.correction;
      }

//if(mype == 0) printf("MYPE %d: Line %d Iteration %d \t local_sum = %12.6lg, global_sum = %12.6lg\n", mype, __LINE__, mesh->m_ncycle, local.sum, global.sum);

#else
      total_sum = local.sum + local.correction;
#endif
#ifdef HAVE_REPROBLAS
   } else if (enhanced_precision_sum == SUM_REPROBLAS_DOUBLE_DOUBLE){
      int fold = 2;
      double ss[4] = {0.0, 0.0, 0.0, 0.0};
      for (uint ic=0; ic < ncells; ic++){
         if (celltype[ic] == REAL_CELL) {
            dndpd(fold, ss, H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]]);
         }
      }
      total_sum = ss[0] + ss[1];
#ifdef HAVE_MPI
      if (mesh->parallel) {
         struct esum_type local, global;
         local.sum = ss[0];
         local.correction = ss[1];
         MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KAHAN_SUM, MPI_COMM_WORLD);
         total_sum = global.sum + global.correction;
      } else {
         total_sum = ss[0] + ss[1];
      }
#else
      total_sum = ss[0] + ss[1];
#endif
   } else if (enhanced_precision_sum == SUM_REPROBLAS_INDEXEDFP) {
      //printf("DEBUG -- reproblas_indexedfp_sum\n");
      Idouble IFPsummer, IFPtotal_sum;
      dISetZero(IFPsummer);
      for (uint ic=0; ic < ncells; ic++){
         if (celltype[ic] == REAL_CELL) {
            dIAddd(&IFPsummer, H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]]);
         }
      }
#ifdef HAVE_MPI
      if (mesh->parallel) {
         MPI_Allreduce(&IFPsummer, &IFPtotal_sum, 1, MPI_IDOUBLE, MPI_RSUM, MPI_COMM_WORLD);
      } else {
         IFPtotal_sum = IFPsummer;
      }
#else
      IFPtotal_sum = IFPsummer;
#endif
         
      total_sum = Iconv2d(IFPtotal_sum);
#endif
   } else if (enhanced_precision_sum == SUM_REGULAR) {
      //printf("DEBUG -- regular_sum\n");
      for (uint ic=0; ic < ncells; ic++){
         if (celltype[ic] == REAL_CELL) {
            summer += H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]];
         }
      }
#ifdef HAVE_MPI
      if (mesh->parallel) {
         MPI_Allreduce(&summer, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      } else {
         total_sum = summer;
      }
#else
      total_sum = summer;
#endif
   }

   cpu_time_mass_sum += cpu_timer_stop(tstart_cpu);

   return(total_sum);
}

#ifdef HAVE_OPENCL
double State::gpu_mass_sum(int enhanced_precision_sum)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_command_queue command_queue = ezcl_get_command_queue();

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
      dev_mass_sum = ezcl_malloc(NULL, const_cast<char *>("dev_mass_sum"), &one,    sizeof(cl_real2_t), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_real2_t), CL_MEM_READ_WRITE, 0);

        /*
        __kernel void reduce_sum_cl(
                         const int      isize,      // 0
                __global       state_t *array,      // 1   Array to be reduced.
                __global       int     *level,      // 2
                __global       int     *levdx,      // 3
                __global       int     *levdy,      // 4
                __global       int     *celltype,   // 5
                __global       real_t  *redscratch, // 6   Final result of operation.
                __local        real_t  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 8, local_work_size*sizeof(cl_real2_t), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, NULL);

      if (block_size > 1) {
           /*
           __kernel void reduce_sum_cl(
                            const int      isize,      // 0
                   __global       int     *redscratch, // 1   Array to be reduced.
                   __local        real_t  *tile)       // 2
           */

         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 3, local_work_size*sizeof(cl_real2_t), NULL);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, NULL);
      }

      struct esum_type local, global;
      real2_t mass_sum;

      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real2_t), &mass_sum, NULL);

      local.sum = mass_sum.s0;
      local.correction = mass_sum.s1;
      global.sum = local.sum;
      global.correction = local.correction;
#ifdef HAVE_MPI
      MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KAHAN_SUM, MPI_COMM_WORLD);
#endif
      gpu_mass_sum = global.sum + global.correction;
   } else {
      dev_mass_sum = ezcl_malloc(NULL, const_cast<char *>("dev_mass_sum"), &one,    sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);

        /*
        __kernel void reduce_sum_cl(
                         const int      isize,      // 0
                __global       state_t *array,      // 1   Array to be reduced.
                __global       int     *level,      // 2
                __global       int     *levdx,      // 3
                __global       int     *levdy,      // 4
                __global       int     *celltype,   // 5
                __global       real_t  *redscratch, // 6   Final result of operation.
                __local        real_t  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 8, local_work_size*sizeof(cl_real_t), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, NULL);

      if (block_size > 1) {
           /*
           __kernel void reduce_sum_cl(
                            const int     isize,      // 0
                   __global       int    *redscratch, // 1   Array to be reduced.
                   __local        real_t  *tile)       // 2
           */

         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 3, local_work_size*sizeof(cl_real_t), NULL);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, NULL);
      }

      double local_sum, global_sum;
      real_t mass_sum;

      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real_t), &mass_sum, NULL);
      
      local_sum = mass_sum;
      global_sum = local_sum;
#ifdef HAVE_MPI
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
      gpu_mass_sum = global_sum;
   }

   ezcl_device_memory_delete(dev_redscratch);
   ezcl_device_memory_delete(dev_mass_sum);

   gpu_time_mass_sum += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   return(gpu_mass_sum);
}
#endif

#ifdef HAVE_OPENCL
void State::allocate_device_memory(size_t ncells)
{
   dev_H = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_H"));
   dev_U = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_U"));
   dev_V = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_V"));
}
#endif

void State::resize_old_device_memory(size_t ncells)
{
#ifdef HAVE_OPENCL
   gpu_state_memory.memory_delete(dev_H);
   gpu_state_memory.memory_delete(dev_U);
   gpu_state_memory.memory_delete(dev_V);
   dev_H = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_H"));
   dev_U = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_U"));
   dev_V = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), DEVICE_REGULAR_MEMORY, const_cast<char *>("dev_V"));
#else
   // Just to block compiler warnings
   if (1 == 2) printf("DEBUG -- ncells is %ld\n",ncells);
#endif
}

#ifdef HAVE_MPI
void State::do_load_balance_local(size_t &numcells){
   mesh->do_load_balance_local(numcells, NULL, state_memory);
   memory_reset_ptrs();
}
#endif
#ifdef HAVE_OPENCL
#ifdef HAVE_MPI
void State::gpu_do_load_balance_local(size_t &numcells){
   if (mesh->gpu_do_load_balance_local(numcells, NULL, gpu_state_memory) ){
      //gpu_state_memory.memory_report();
      dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H");
      dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U");
      dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V");
/*
      if (dev_H == NULL){
         dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H_new");
         dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U_new");
         dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V_new");
      }
      printf("DEBUG memory for proc %d dev_H is %p dev_U is %p dev_V is %p\n",mesh->mype,dev_H,dev_U,dev_V);
*/
   }
}
#endif
#endif

static double total_time = 0.0;

void State::output_timing_info(int do_cpu_calc, int do_gpu_calc, double elapsed_time)
{
   int &mype  = mesh->mype;
   int &numpe = mesh->numpe;

   int rank = mype;
   if (! mesh->parallel) {
      // We need to get rank info for check routines
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
   }

   double cpu_time_compute;
   long gpu_time_compute;

   double cpu_elapsed_time=0.0;
   long gpu_elapsed_time=0;

   if (rank == 0) {
      printf("\n");
      printf("~~~~~~~~~~~~~~~~ Device timing information ~~~~~~~~~~~~~~~~~~\n");
   }
   if (! mesh->parallel) {
      
      //  Output timing information.
      if (do_cpu_calc) {
         cpu_time_compute = get_cpu_time_set_timestep() +
                            get_cpu_time_finite_difference() +
                            get_cpu_time_refine_potential() +
                            mesh->get_cpu_time_rezone_all() +
                            mesh->get_cpu_time_calc_neighbors() +
                            get_cpu_time_mass_sum() +
                            mesh->get_cpu_time_calc_spatial_coordinates() +
                            mesh->cpu_time_partition;
         cpu_elapsed_time =                  cpu_time_compute;
         total_time = cpu_elapsed_time;

         if (rank == 0) {
            printf("CPU: Device compute           time was\t%8.4f \ts\n",     cpu_time_compute);
            printf("CPU:  state->set_timestep      time was\t %8.4f\ts\n",     get_cpu_time_set_timestep() );
            printf("CPU:  state->finite_difference time was\t %8.4f\ts\n",     get_cpu_time_finite_difference() );
            printf("CPU:  mesh->refine_potential   time was\t %8.4f\ts\n",     get_cpu_time_refine_potential() );
            printf("CPU:    mesh->calc_mpot          time was\t %8.4f\ts\n",     get_cpu_time_calc_mpot() );
            printf("CPU:    mesh->refine_smooth      time was\t %8.4f\ts\n",     mesh->get_cpu_time_refine_smooth() );
            printf("CPU:  mesh->rezone_all         time was\t %8.4f\ts\n",     mesh->get_cpu_time_rezone_all() );
            printf("CPU:  mesh->partition_cells    time was\t %8.4f\ts\n",     mesh->cpu_time_partition);
            printf("CPU:  mesh->calc_neighbors     time was\t %8.4f\ts\n",     mesh->get_cpu_time_calc_neighbors() );
            if (mesh->get_calc_neighbor_type() == HASH_TABLE) {
              printf("CPU:    mesh->hash_setup         time was\t %8.4f\ts\n",     mesh->get_cpu_time_hash_setup() );
              printf("CPU:    mesh->hash_query         time was\t %8.4f\ts\n",     mesh->get_cpu_time_hash_query() );
            } else {
              printf("CPU:    mesh->kdtree_setup         time was\t %8.4f\ts\n",     mesh->get_cpu_time_kdtree_setup() );
              printf("CPU:    mesh->kdtree_query         time was\t %8.4f\ts\n",     mesh->get_cpu_time_kdtree_query() );
            }
            printf("CPU:  mass_sum                time was\t %8.4f\ts\n",     get_cpu_time_mass_sum() );
            printf("CPU:  mesh->calc_spatial_coor time was\t %8.4f\ts\n",     mesh->get_cpu_time_calc_spatial_coordinates() );
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
                            mesh->get_gpu_time_rezone_all() +
                            mesh->get_gpu_time_calc_neighbors() +
                            get_gpu_time_mass_sum() +
                            mesh->get_gpu_time_calc_spatial_coordinates() +
                            mesh->gpu_time_count_BCs;
         gpu_elapsed_time   = get_gpu_time_write() + gpu_time_compute + get_gpu_time_read();

         if (rank == 0) {
            printf("GPU: Write to device          time was\t%8.4f\ts\n",    (double) get_gpu_time_write()         * 1.0e-9); /* Convert nanoseconds to msecs */
            printf("GPU: Read from device         time was\t%8.4f\ts\n",    (double) get_gpu_time_read()          * 1.0e-9);
            printf("GPU: Device compute           time was\t%8.4f\ts\n",    (double) gpu_time_compute             * 1.0e-9);
            printf("GPU:  kernel_set_timestep      time was\t %8.4f\ts\n",    (double) get_gpu_time_set_timestep()      * 1.0e-9);
            printf("GPU:  kernel_calc_finite_diff  time was\t %8.4f\ts\n",    (double) get_gpu_time_finite_difference() * 1.0e-9);
            printf("GPU:  kernel_refine_potential  time was\t %8.4f\ts\n",    (double) get_gpu_time_refine_potential()  * 1.0e-9);
            printf("GPU:    kernel_calc_mpot         time was\t %8.4f\ts\n",    (double) get_gpu_time_calc_mpot()  * 1.0e-9);
            printf("GPU:    kernel_refine_smooth     time was\t %8.4f\ts\n",    (double) mesh->get_gpu_time_refine_smooth()  * 1.0e-9);
            printf("GPU:  kernel_rezone_all        time was\t %8.4f\ts\n",    (double) (mesh->get_gpu_time_rezone_all() ) * 1.0e-9);
            printf("GPU:  kernel_calc_neighbors    time was\t %8.4f\ts\n",    (double) mesh->get_gpu_time_calc_neighbors()    * 1.0e-9);
            if (mesh->get_calc_neighbor_type() == HASH_TABLE) {
              printf("GPU:    mesh->hash_setup         time was\t %8.4f\ts\n",     (double)mesh->get_gpu_time_hash_setup() * 1.0e-9 );
              printf("GPU:    mesh->hash_query         time was\t %8.4f\ts\n",     (double)mesh->get_gpu_time_hash_query() * 1.0e-9 );
            } else {
              printf("GPU:    mesh->kdtree_setup         time was\t %8.4f\ts\n",     (double)mesh->get_gpu_time_kdtree_setup() * 1.0e-9 );
              printf("GPU:    mesh->kdtree_query         time was\t %8.4f\ts\n",     (double)mesh->get_gpu_time_kdtree_query() * 1.0e-9 );
            }
            printf("GPU:  kernel_mass_sum          time was\t %8.4f\ts\n",    (double) get_gpu_time_mass_sum()          * 1.0e-9);
            printf("GPU:  kernel_calc_spatial_coor time was\t %8.4f\ts\n",    (double) mesh->get_gpu_time_calc_spatial_coordinates() * 1.0e-9);
            if (! mesh->have_boundary) {
               printf("GPU:  kernel_count_BCs         time was\t %8.4f\ts\n",    (double) mesh->gpu_time_count_BCs        * 1.0e-9);
            }
            printf("=============================================================\n");
            printf("Profiling: Total GPU          time was\t%8.4f\ts or\t%4.2f min\n", (double) gpu_elapsed_time * 1.0e-9, (double) gpu_elapsed_time * 1.0e-9/60.0);
            printf("=============================================================\n\n");
         }
      }
      if (rank == 0) {
         printf("=============================================================\n");
         printf("Profiling: Total              time was\t%8.4f\ts or\t%4.2f min\n",    elapsed_time,     elapsed_time/60.0);
         if (do_cpu_calc && do_gpu_calc) {
            printf("Speed-up:                             \t%8.4f\n",   cpu_elapsed_time*1.0e9/gpu_elapsed_time);
         }
         printf("=============================================================\n");
      }
   } else {
      if (do_cpu_calc) {
         cpu_time_compute = get_cpu_time_set_timestep() +
                            get_cpu_time_finite_difference() +
                            get_cpu_time_refine_potential() +
                            mesh->get_cpu_time_rezone_all() +
                            mesh->get_cpu_time_calc_neighbors() +
                            mesh->get_cpu_time_load_balance() +
                            get_cpu_time_mass_sum() +
                            mesh->get_cpu_time_calc_spatial_coordinates() +
                            mesh->cpu_time_partition;
         cpu_elapsed_time =                  cpu_time_compute;

         if (mype == 0) printf("\nCPU: Parallel timings\n\n");

         parallel_timer_output(numpe,mype,"CPU: Device compute           time was" ,cpu_time_compute);
         parallel_timer_output(numpe,mype,"CPU:  state->set_timestep      time was",get_cpu_time_set_timestep() );
         parallel_timer_output(numpe,mype,"CPU:  state->finite_difference time was",get_cpu_time_finite_difference() );
         parallel_timer_output(numpe,mype,"CPU:  state->refine_potential  time was",get_cpu_time_refine_potential() );
         parallel_timer_output(numpe,mype,"CPU:    state->calc_mpot         time was",get_cpu_time_calc_mpot() );
         parallel_timer_output(numpe,mype,"CPU:    state->refine_smooth     time was",mesh->get_cpu_time_refine_smooth() );
         parallel_timer_output(numpe,mype,"CPU:  mesh->rezone_all         time was",mesh->get_cpu_time_rezone_all() );
         parallel_timer_output(numpe,mype,"CPU:  mesh->partition_cells    time was",mesh->cpu_time_partition);
         parallel_timer_output(numpe,mype,"CPU:  mesh->calc_neighbors     time was",mesh->get_cpu_time_calc_neighbors() );
         if (mesh->get_calc_neighbor_type() == HASH_TABLE) {
           parallel_timer_output(numpe,mype,"CPU:    mesh->hash_setup         time was",mesh->get_cpu_time_hash_setup() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->hash_query         time was",mesh->get_cpu_time_hash_query() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->find_boundary      time was",mesh->get_cpu_time_find_boundary() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->push_setup         time was",mesh->get_cpu_time_push_setup() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->push_boundary      time was",mesh->get_cpu_time_push_boundary() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->local_list         time was",mesh->get_cpu_time_local_list() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->layer1             time was",mesh->get_cpu_time_layer1() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->layer2             time was",mesh->get_cpu_time_layer2() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->layer_list         time was",mesh->get_cpu_time_layer_list() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->copy_mesh_data     time was",mesh->get_cpu_time_copy_mesh_data() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->fill_mesh_ghost    time was",mesh->get_cpu_time_fill_mesh_ghost() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->fill_neigh_ghost   time was",mesh->get_cpu_time_fill_neigh_ghost() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->set_corner_neigh   time was",mesh->get_cpu_time_set_corner_neigh() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->neigh_adjust       time was",mesh->get_cpu_time_neigh_adjust() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->setup_comm         time was",mesh->get_cpu_time_setup_comm() );
         } else {
           parallel_timer_output(numpe,mype,"CPU:    mesh->kdtree_setup         time was",mesh->get_cpu_time_kdtree_setup() );
           parallel_timer_output(numpe,mype,"CPU:    mesh->kdtree_query         time was",mesh->get_cpu_time_kdtree_query() );
         }
         parallel_timer_output(numpe,mype,"CPU:  mesh->calc_load_balance  time was",mesh->get_cpu_time_load_balance() );
         parallel_timer_output(numpe,mype,"CPU:  mesh->calc_spatial_coor  time was",mesh->get_cpu_time_calc_spatial_coordinates() );
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
                            mesh->get_gpu_time_rezone_all() +
                            mesh->get_gpu_time_calc_neighbors() +
                            mesh->get_gpu_time_load_balance() +
                            get_gpu_time_mass_sum() +
                            mesh->get_gpu_time_calc_spatial_coordinates() +
                            mesh->gpu_time_count_BCs;
         gpu_elapsed_time   = get_gpu_time_write() + gpu_time_compute + get_gpu_time_read();

         if (mype == 0) printf("\nGPU: Parallel timings\n\n");

         parallel_timer_output(numpe,mype,"GPU: Write to device          time was", (double) get_gpu_time_write()         * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU: Read from device         time was", (double) get_gpu_time_read()          * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU: Device compute           time was", (double) gpu_time_compute             * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_set_timestep      time was",(double) get_gpu_time_set_timestep()      * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_calc_finite_diff  time was",(double) get_gpu_time_finite_difference() * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_refine_potential  time was",(double) get_gpu_time_refine_potential()  * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:    kernel_calc_mpot         time was",(double) get_gpu_time_calc_mpot()  * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:    kernel_refine_smooth     time was",(double) mesh->get_gpu_time_refine_smooth()  * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_rezone_all        time was",(double) (mesh->get_gpu_time_rezone_all() ) * 1.0e-9 );
         //parallel_timer_output(numpe,mype,"GPU:  kernel_hash_setup        time was",(double) mesh->get_gpu_time_hash_setup()         * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_calc_neighbors    time was",(double) mesh->get_gpu_time_calc_neighbors()     * 1.0e-9 );
         if (mesh->get_calc_neighbor_type() == HASH_TABLE) {
           parallel_timer_output(numpe,mype,"GPU:    kernel_hash_setup        time was",(double) mesh->get_gpu_time_hash_setup()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_hash_query        time was",(double) mesh->get_gpu_time_hash_query()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_find_boundary     time was",(double) mesh->get_gpu_time_find_boundary()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_push_setup        time was",(double) mesh->get_gpu_time_push_setup()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_push_boundary     time was",(double) mesh->get_gpu_time_push_boundary()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_local_list        time was",(double) mesh->get_gpu_time_local_list()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_layer1            time was",(double) mesh->get_gpu_time_layer1()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_layer2            time was",(double) mesh->get_gpu_time_layer2()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_layer_list        time was",(double) mesh->get_gpu_time_layer_list()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_copy_mesh_data    time was",(double) mesh->get_gpu_time_copy_mesh_data()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_fill_mesh_ghost   time was",(double) mesh->get_gpu_time_fill_mesh_ghost()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_fill_neigh_ghost  time was",(double) mesh->get_gpu_time_fill_neigh_ghost()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_set_corner_neigh  time was",(double) mesh->get_gpu_time_set_corner_neigh()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_neigh_adjust      time was",(double) mesh->get_gpu_time_neigh_adjust()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_setup_comm        time was",(double) mesh->get_gpu_time_setup_comm()     * 1.0e-9 );
         } else {
           parallel_timer_output(numpe,mype,"GPU:    kernel_kdtree_setup        time was",(double) mesh->get_gpu_time_kdtree_setup()     * 1.0e-9 );
           parallel_timer_output(numpe,mype,"GPU:    kernel_kdtree_query        time was",(double) mesh->get_gpu_time_kdtree_query()     * 1.0e-9 );
         }
         parallel_timer_output(numpe,mype,"GPU:  kernel_calc_spatial_coor time was",(double) mesh->get_gpu_time_calc_spatial_coordinates()     * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_calc_load_balance time was",(double) mesh->get_gpu_time_load_balance()     * 1.0e-9 );
         parallel_timer_output(numpe,mype,"GPU:  kernel_mass_sum          time was",(double) get_gpu_time_mass_sum()          * 1.0e-9 );

         if (! mesh->have_boundary) {
            parallel_timer_output(numpe,mype,"GPU:  kernel_count_BCs         time was",(double) mesh->gpu_time_count_BCs        * 1.0e-9 );
         }

         if (mype == 0) printf("=============================================================\n");

         parallel_timer_output(numpe,mype,"Profiling: Total GPU          time was",gpu_elapsed_time * 1.0e-9 );

         if (mype == 0) printf("=============================================================\n\n");
      }

      if (mype == 0) printf("=============================================================\n");

      parallel_timer_output(numpe,mype,"Profiling: Total              time was",elapsed_time );
      if (do_gpu_calc){
         parallel_timer_output(numpe,mype,"Parallel Speed-up:                    ",cpu_elapsed_time*1.0e9/gpu_elapsed_time);
      } else {
         if (total_time > 0.0) parallel_timer_output(numpe,mype,"Parallel Speed-up:                    ",total_time/cpu_elapsed_time);
      }


      if (mype == 0) printf("=============================================================\n");
   }
}

void State::parallel_timer_output(int numpe, int mype, const char *string, double local_time)
{
   vector<double> global_times(numpe);
#ifdef HAVE_MPI
   MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_times[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
   // Just to block compiler warning
   if (1 == 2) printf("DEBUG -- local time is %lf\n",local_time);
#endif
   if (mype == 0) {
      printf("%s\t",string);
      if (numpe <= 4) {
         for(int ip = 0; ip < numpe; ip++){
            printf("%8.4f\t", global_times[ip]);
         }
         printf("s\n");
      } else {
         sort(global_times.begin(),global_times.end());
         double median_value;
         int half_value = numpe/2;
         if (numpe%2 == 0) {
            median_value = (global_times[half_value-1]+global_times[half_value])/2.0;
         } else {
            median_value = global_times[half_value+1];
         }
         printf(" %8.4f\t %8.4f\t %8.4f secs min/median/max\n",global_times[0],median_value,global_times[numpe-1]);
      }
   }
}

void State::parallel_memory_output(int numpe, int mype, const char *string, long long local_time)
{
   vector<long long> global_memory_value(numpe);
#ifdef HAVE_MPI
   MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_memory_value[0], 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#else
   // Just to block compiler warning
   if (1 == 2) printf("DEBUG -- local time is %llu\n",local_time);
#endif
   if (mype == 0) {
      printf("%s\t",string);
      if (numpe <= 4) {
         for(int ip = 0; ip < numpe; ip++){
            printf("%10lld\t", global_memory_value[ip]);
         }
         printf("kB\n");
      } else {
         sort(global_memory_value.begin(),global_memory_value.end());
         long long median_value;
         int half_value = numpe/2;
         if (numpe%2 == 0) {
            median_value = (global_memory_value[half_value-1]+global_memory_value[half_value])/2;
         } else {
            median_value = global_memory_value[half_value+1];
         }
         printf(" %10lld\t %10lld\t %10lld kb min/median/max\n",global_memory_value[0],median_value,global_memory_value[numpe-1]);
      }
   }
}

#ifdef HAVE_OPENCL
void State::compare_state_gpu_global_to_cpu_global(const char* string, int cycle, uint ncells)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<state_t>H_check(ncells);
   vector<state_t>U_check(ncells);
   vector<state_t>V_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t), &H_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t), &U_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_state_t), &V_check[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic]-H_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d H & H_check %d %lf %lf\n",string,cycle,ic,H[ic],H_check[ic]);
      if (fabs(U[ic]-U_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d U & U_check %d %lf %lf\n",string,cycle,ic,U[ic],U_check[ic]);
      if (fabs(V[ic]-V_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d V & V_check %d %lf %lf\n",string,cycle,ic,V[ic],V_check[ic]);
   }
}
#endif

void State::compare_state_cpu_local_to_cpu_global(State *state_global, const char* string, int cycle, uint ncells, uint ncells_global, int *nsizes, int *ndispl)
{
   state_t *H_global = state_global->H;
   state_t *U_global = state_global->U;
   state_t *V_global = state_global->V;

   vector<state_t>H_check(ncells_global);
   vector<state_t>U_check(ncells_global);
   vector<state_t>V_check(ncells_global);
#ifdef HAVE_MPI
   MPI_Allgatherv(&H[0], ncells, MPI_STATE_T, &H_check[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&U[0], ncells, MPI_STATE_T, &U_check[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&V[0], ncells, MPI_STATE_T, &V_check[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
#else
   // Just to block compiler warnings
   if (1 == 2) printf("DEBUG -- ncells %u nsizes %d ndispl %d\n",ncells, nsizes[0],ndispl[0]);
#endif

   for (uint ic = 0; ic < ncells_global; ic++){
      if (fabs(H_global[ic]-H_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d H & H_check %d %lf %lf\n",string,cycle,ic,H_global[ic],H_check[ic]);
      if (fabs(U_global[ic]-U_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d U & U_check %d %lf %lf\n",string,cycle,ic,U_global[ic],U_check[ic]);
      if (fabs(V_global[ic]-V_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d V & V_check %d %lf %lf\n",string,cycle,ic,V_global[ic],V_check[ic]);
   }
}

#ifdef HAVE_OPENCL
void State::compare_state_all_to_gpu_local(State *state_global, uint ncells, uint ncells_global, int mype, int ncycle, int *nsizes, int *ndispl)
{
#ifdef HAVE_MPI
   cl_command_queue command_queue = ezcl_get_command_queue();

   state_t *H_global = state_global->H;
   state_t *U_global = state_global->U;
   state_t *V_global = state_global->V;
   cl_mem &dev_H_global = state_global->dev_H;
   cl_mem &dev_U_global = state_global->dev_U;
   cl_mem &dev_V_global = state_global->dev_V;

   // Need to compare dev_H to H, etc
   vector<state_t>H_save(ncells);
   vector<state_t>U_save(ncells);
   vector<state_t>V_save(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t), &H_save[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t), &U_save[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_state_t), &V_save[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic]-H_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d H & H_save %d %lf %lf \n",mype,ncycle,ic,H[ic],H_save[ic]);
      if (fabs(U[ic]-U_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d U & U_save %d %lf %lf \n",mype,ncycle,ic,U[ic],U_save[ic]);
      if (fabs(V[ic]-V_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d V & V_save %d %lf %lf \n",mype,ncycle,ic,V[ic],V_save[ic]);
   }

   // And compare dev_H gathered to H_global, etc
   vector<state_t>H_save_global(ncells_global);
   vector<state_t>U_save_global(ncells_global);
   vector<state_t>V_save_global(ncells_global);
   MPI_Allgatherv(&H_save[0], nsizes[mype], MPI_STATE_T, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&U_save[0], nsizes[mype], MPI_STATE_T, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&V_save[0], nsizes[mype], MPI_STATE_T, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   if (mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d H_global & H_save_global %d %lf %lf \n",mype,ncycle,ic,H_global[ic],H_save_global[ic]);
         if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d U_global & U_save_global %d %lf %lf \n",mype,ncycle,ic,U_global[ic],U_save_global[ic]);
         if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d V_global & V_save_global %d %lf %lf \n",mype,ncycle,ic,V_global[ic],V_save_global[ic]);
      }
   }

   // And compare H gathered to H_global, etc
   MPI_Allgatherv(&H[0], nsizes[mype], MPI_STATE_T, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&U[0], nsizes[mype], MPI_STATE_T, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&V[0], nsizes[mype], MPI_STATE_T, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   if (mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d H_global & H_save_global %d %lf %lf \n",ncycle,ic,H_global[ic],H_save_global[ic]);
         if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d U_global & U_save_global %d %lf %lf \n",ncycle,ic,U_global[ic],U_save_global[ic]);
         if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d V_global & V_save_global %d %lf %lf \n",ncycle,ic,V_global[ic],V_save_global[ic]);
      }
   }

   // Now the global dev_H_global to H_global, etc
   ezcl_enqueue_read_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_state_t), &H_save_global[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_state_t), &U_save_global[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V_global, CL_TRUE,  0, ncells_global*sizeof(cl_state_t), &V_save_global[0], NULL);
   if (mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d H_global & H_save_global %d %lf %lf \n",mype,ncycle,ic,H_global[ic],H_save_global[ic]);
         if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d U_global & U_save_global %d %lf %lf \n",mype,ncycle,ic,U_global[ic],U_save_global[ic]);
         if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d V_global & V_save_global %d %lf %lf \n",mype,ncycle,ic,V_global[ic],V_save_global[ic]);
      }
   }
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("%d: DEBUG -- ncells %d ncells_global %d ncycle %d nsizes[0] %d ndispl %d state_global %p\n",
      mype,ncells,ncells_global,ncycle,nsizes[0],ndispl[0],state_global);
#endif
}
#endif

void State::print_object_info(void)
{
   printf(" ---- State object info -----\n");

#ifdef HAVE_OPENCL
   int num_elements, elsize;

   num_elements = ezcl_get_device_mem_nelements(dev_H);
   elsize = ezcl_get_device_mem_elsize(dev_H);
   printf("dev_H       ptr : %p nelements %d elsize %d\n",dev_H,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_U);
   elsize = ezcl_get_device_mem_elsize(dev_U);
   printf("dev_U       ptr : %p nelements %d elsize %d\n",dev_U,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_V);
   elsize = ezcl_get_device_mem_elsize(dev_V);
   printf("dev_V       ptr : %p nelements %d elsize %d\n",dev_V,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_mpot);
   elsize = ezcl_get_device_mem_elsize(dev_mpot);
   printf("dev_mpot    ptr : %p nelements %d elsize %d\n",dev_mpot,num_elements,elsize);
   //num_elements = ezcl_get_device_mem_nelements(dev_ioffset);
   //elsize = ezcl_get_device_mem_elsize(dev_ioffset);
   //printf("dev_ioffset ptr : %p nelements %d elsize %d\n",dev_ioffset,num_elements,elsize);
#endif
   state_memory.memory_report();
   //printf("vector H    ptr : %p nelements %ld elsize %ld\n",&H[0],H.size(),sizeof(H[0]));
   //printf("vector U    ptr : %p nelements %ld elsize %ld\n",&U[0],U.size(),sizeof(U[0]));
   //printf("vector V    ptr : %p nelements %ld elsize %ld\n",&V[0],V.size(),sizeof(V[0]));
}

void State::print(void)
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   if (mesh->fp == NULL) {
      char filename[10];
      sprintf(filename,"out%1d",mesh->mype);
      mesh->fp=fopen(filename,"w");
   }

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

const int CRUX_STATE_VERSION = 101;
const int num_double_vals    = 7;
const int num_long_vals      = 10;

size_t State::get_checkpoint_size(void)
{
#ifdef FULL_PRECISION
   size_t nsize = mesh->ncells*3*sizeof(double);
#else
   size_t nsize = mesh->ncells*3*sizeof(float);
#endif
   nsize += num_double_vals*sizeof(double)+num_long_vals*sizeof(long long);
   nsize += mesh->get_checkpoint_size();
   return(nsize);
}

void State::store_checkpoint(Crux *crux)
{
   mesh->store_checkpoint(crux);

   double double_vals[num_double_vals];

   double_vals[0] = cpu_time_apply_BCs,
   double_vals[1] = cpu_time_set_timestep,
   double_vals[2] = cpu_time_finite_difference,
   double_vals[3] = cpu_time_refine_potential,
   double_vals[4] = cpu_time_calc_mpot,
   double_vals[5] = cpu_time_rezone_all,
   double_vals[6] = cpu_time_mass_sum;

   crux->store_doubles(double_vals, num_double_vals);

   long long long_vals[num_long_vals];

   long_vals[0] = CRUX_STATE_VERSION;
   long_vals[1] = gpu_time_apply_BCs,
   long_vals[2] = gpu_time_set_timestep,
   long_vals[3] = gpu_time_finite_difference,
   long_vals[4] = gpu_time_refine_potential,
   long_vals[5] = gpu_time_calc_mpot,
   long_vals[6] = gpu_time_rezone_all,
   long_vals[7] = gpu_time_mass_sum,
   long_vals[8] = gpu_time_read,
   long_vals[9] = gpu_time_write;

   crux->store_longs(long_vals, num_long_vals);

#ifdef FULL_PRECISION
   crux->store_double_array(H, mesh->ncells);
   crux->store_double_array(U, mesh->ncells);
   crux->store_double_array(V, mesh->ncells);
#else
   crux->store_float_array(H, mesh->ncells);
   crux->store_float_array(U, mesh->ncells);
   crux->store_float_array(V, mesh->ncells);
#endif
}

void State::restore_checkpoint(Crux *crux)
{
   mesh->restore_checkpoint(crux);

   double double_vals[num_double_vals];

   crux->restore_doubles(double_vals, num_double_vals);

   cpu_time_apply_BCs         = double_vals[0];
   cpu_time_set_timestep      = double_vals[1];
   cpu_time_finite_difference = double_vals[2];
   cpu_time_refine_potential  = double_vals[3];
   cpu_time_calc_mpot         = double_vals[4];
   cpu_time_rezone_all        = double_vals[5];
   cpu_time_mass_sum          = double_vals[6];

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS) {
      const char *double_vals_descriptor[num_double_vals] = {
         "cpu_time_apply_BCs",
         "cpu_time_set_timestep",
         "cpu_time_finite_difference",
         "cpu_time_refine_potential",
         "cpu_time_calc_mpot",
         "cpu_time_rezone_all",
         "cpu_time_mass_sum"
      };
      printf("\n");
      printf("       === Restored state double_vals ===\n");
      for (int i = 0; i < num_double_vals; i++){
         printf("       %-30s %lg\n",double_vals_descriptor[i], double_vals[i]);
      }
      printf("       === Restored state double_vals ===\n");
      printf("\n");
   }
#endif

   long long long_vals[num_long_vals];

   crux->restore_longs(long_vals, num_long_vals);

   if (long_vals[ 0] != CRUX_STATE_VERSION) {
      printf("CRUX version mismatch for state data, version on file is %lld, version in code is %d\n",
         long_vals[0], CRUX_STATE_VERSION);
      exit(0);
   }

   gpu_time_apply_BCs         = long_vals[1];
   gpu_time_set_timestep      = long_vals[2];
   gpu_time_finite_difference = long_vals[3];
   gpu_time_refine_potential  = long_vals[4];
   gpu_time_calc_mpot         = long_vals[5];
   gpu_time_rezone_all        = long_vals[6];
   gpu_time_mass_sum          = long_vals[7];
   gpu_time_read              = long_vals[8];
   gpu_time_write             = long_vals[9];

#ifdef DEBUG_RESTORED_VALS
   if (DEBUG_RESTORED_VALS) {
      const char *long_vals_descriptor[num_long_vals] = {
         "CRUX_STATE_VERSION",
         "gpu_time_apply_BCs",
         "gpu_time_set_timestep",
         "gpu_time_finite_difference",
         "gpu_time_refine_potential",
         "gpu_time_calc_mpot",
         "gpu_time_rezone_all",
         "gpu_time_mass_sum",
         "gpu_time_read",
         "gpu_time_write"
      };
      printf("\n");
      printf("       === Restored state long_vals ===\n");
      for (int i = 0; i < num_long_vals; i++){
         printf("       %-30s %lld\n",long_vals_descriptor[i], long_vals[i]);
      }
      printf("       === Restored state long_vals ===\n");
      printf("\n");
   }
#endif

   allocate(mesh->ncells);

#ifdef FULL_PRECISION
   H = crux->restore_double_array(H, mesh->ncells);
   U = crux->restore_double_array(U, mesh->ncells);
   V = crux->restore_double_array(V, mesh->ncells);
#else
   H = crux->restore_float_array(H, mesh->ncells);
   U = crux->restore_float_array(U, mesh->ncells);
   V = crux->restore_float_array(V, mesh->ncells);
#endif
}

// Added overloaded print to get mesh information to print in each cycle
// Brian Atkinson (5-29-14)
void State::print(int iteration, double simTime, double initial_mass, double iteration_mass, double mass_diff_percentage)
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

      char filename[40];
      sprintf(filename,"iteration%d",iteration);
      mesh->fp=fopen(filename,"w");

      if(iteration_mass == 0.0){
         fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
         fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
         fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tSimulation Time: %lf\n", initial_mass, simTime);
      }
      else{
         double mass_diff = iteration_mass - initial_mass;
         fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
         fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
         fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tIteration Mass: %14.12lg\n", initial_mass, iteration_mass);
         fprintf(mesh->fp,"Mass Difference: %12.6lg\t\tMass Difference Percentage: %12.6lg%%\n", mass_diff, mass_diff_percentage);
      }

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

void State::print_local(int ncycle)
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   if (mesh->fp == NULL) {
      char filename[10];
      sprintf(filename,"out%1d",mesh->mype);
      mesh->fp=fopen(filename,"w");
   }

   fprintf(mesh->fp,"DEBUG in print_local ncycle is %d\n",ncycle);
   if (mesh->nlft != NULL){
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev   nlft   nrht   nbot   ntop\n",mesh->mype);
      uint state_size = state_memory.get_memory_size(H);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         if (ic >= state_size){
            fprintf(mesh->fp,"%d: %6d                              %4d  %4d   %4d  %4d  %4d  %4d  %4d\n", mesh->mype,ic, mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
         } else {
            fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  %4d  %4d  %4d  %4d\n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
         }
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d\n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

void State::print_failure_log(int iteration, double simTime, double initial_mass, double iteration_mass, double mass_diff_percentage, bool got_nan){
   char filename[] = {"failure.log"};
   mesh->fp=fopen(filename,"w");

   double mass_diff = iteration_mass - initial_mass;
   if(got_nan){
      fprintf(mesh->fp,"Failed because of nan for H_sum was equal to NAN\n");
   }
   else{
      fprintf(mesh->fp,"Failed because mass difference is outside of accepted percentage\n");
   }
   fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
   fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
   fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tIteration Mass: %14.12lg\n", initial_mass, iteration_mass);
   fprintf(mesh->fp,"Mass Difference: %12.6lg\t\tMass Difference Percentage: %12.6lg%%\n", mass_diff, mass_diff_percentage);

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

void State::print_rollback_log(int iteration, double simTime, double initial_mass, double iteration_mass, double mass_diff_percentage, int backup_attempt, int num_of_attempts, int error_status){
   char filename[40];
   sprintf(filename, "rollback%d.log",backup_attempt);
   mesh->fp=fopen(filename,"w");

   double mass_diff = iteration_mass - initial_mass;
   if(error_status == STATUS_NAN){
      fprintf(mesh->fp,"Rolling back because of nan for H_sum was equal to NAN\n");
   }
   else{
      fprintf(mesh->fp,"Rolling back because mass difference is outside of accepted percentage\n");
   }
   fprintf(mesh->fp,"Rollback attempt %d of %d ---> Number of attempts left:%d\n", backup_attempt, num_of_attempts, num_of_attempts - backup_attempt);
   fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
   fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
   fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tIteration Mass: %14.12lg\n", initial_mass, iteration_mass);
   fprintf(mesh->fp,"Mass Difference: %12.6lg\t\tMass Difference Percentage: %12.6lg%%\n", mass_diff, mass_diff_percentage);

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}
