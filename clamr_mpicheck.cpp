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

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>
#include "input.h"
#include "mesh.h"
#include "partition.h"
#include "state.h"
#include "l7/l7.h"
#include "timer/timer.h"

#include <mpi.h>
#include "display.h"

#ifndef DEBUG 
#define DEBUG 0
#endif

#define DO_COMPARISON

#ifdef DO_COMPARISON
static int do_comparison_calc = 1;
#else
static int do_comparison_calc = 0;
#endif

static int do_cpu_calc = 1;
static int do_gpu_calc = 0;

#ifdef HAVE_CL_DOUBLE
typedef double      real;
typedef struct
{
   double s0;
   double s1;
}  real2;
#define CONSERVATION_EPS    .00001
#define STATE_EPS        .025
#define MPI_C_REAL MPI_DOUBLE
#define L7_REAL L7_DOUBLE
#else
typedef float       real;
typedef struct
{
   float s0;
   float s1;
}  real2;
#define CONSERVATION_EPS    .1
#define STATE_EPS      15.0
#define MPI_C_REAL MPI_FLOAT
#define L7_REAL L7_FLOAT
#endif

typedef unsigned int uint;

static double circle_radius=-1.0;

static int view_mode = 0;

bool        verbose,        //  Flag for verbose command-line output; init in input.cpp::parseInput().
            localStencil,   //  Flag for use of local stencil; init in input.cpp::parseInput().
            outline,        //  Flag for drawing outlines of cells; init in input.cpp::parseInput().
            enhanced_precision_sum,//  Flag for enhanced precision sum (default true); init in input.cpp::parseInput().
            special_case;   //  Flag for special case for debugging (default false); init in input.cpp::parseInput().
int         outputInterval, //  Periodicity of output; init in input.cpp::parseInput().
            levmx,          //  Maximum number of refinement levels; init in input.cpp::parseInput().
            nx,             //  x-resolution of coarse grid; init in input.cpp::parseInput().
            ny,             //  y-resolution of coarse grid; init in input.cpp::parseInput().
            niter,          //  Maximum time step; init in input.cpp::parseInput().
            ndim    = 2;    //  Dimensionality of problem (2 or 3).

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
static Mesh       *mesh_global;    //  Object containing mesh information; init in grid.cpp::main().
static State      *state_global;   //  Object containing state information corresponding to mesh; init in grid.cpp::main().
static Mesh       *mesh;           //  Object containing mesh information; init in grid.cpp::main().
static State      *state;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
static struct timeval tstart;

static double  H_sum_initial = 0.0;

int main(int argc, char **argv) {

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   L7_Init(&mype, &numpe, &argc, argv);
   parseInput(argc, argv);

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   if (special_case) circ_radius = .75;

   mesh_global  = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   mesh_global->init(nx, ny, circ_radius, initial_order, special_case, do_gpu_calc);
   size_t &ncells_global = mesh_global->ncells;
   state_global = new State(ncells_global);
   state_global->init(ncells_global, do_gpu_calc);
   mesh_global->proc.resize(ncells_global);
   mesh_global->calc_distribution(numpe, mesh_global->proc);
   state_global->fill_circle(mesh_global, circ_radius, 100.0, 5.0);
   
   parallel_in = 1;
   mesh = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   state = new State(mesh->ncells);

   size_t &ncells = mesh->ncells;
   int &noffset = mesh->noffset;

   vector<int>   &nsizes     = mesh_global->nsizes;
   vector<int>   &ndispl     = mesh_global->ndispl;

   vector<real>  &H_global = state_global->H;
   vector<real>  &U_global = state_global->U;
   vector<real>  &V_global = state_global->V;

   vector<real>  &x_global  = mesh_global->x;
   vector<real>  &dx_global = mesh_global->dx;
   vector<real>  &y_global  = mesh_global->y;
   vector<real>  &dy_global = mesh_global->dy;

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;
   vector<int>   &proc     = mesh->proc;

   vector<int>   &celltype_global = mesh_global->celltype;
   vector<int>   &i_global        = mesh_global->i;
   vector<int>   &j_global        = mesh_global->j;
   vector<int>   &level_global    = mesh_global->level;
   vector<int>   &proc_global     = mesh_global->proc;

   vector<real> &H = state->H;
   vector<real> &U = state->U;
   vector<real> &V = state->V;

   vector<real> &x  = mesh->x;
   vector<real> &dx = mesh->dx;
   vector<real> &y  = mesh->y;
   vector<real> &dy = mesh->dy;

   ncells = ncells_global/numpe;
   if (mype < ncells_global%numpe) ncells++;

   nsizes.resize(numpe);
   ndispl.resize(numpe);
   MPI_Allgather(&ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
   ndispl[0]=0;
   for (int ip=1; ip<numpe; ip++){
      ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
   }
   noffset=0;
   for (int ip=0; ip<mype; ip++){
     noffset += nsizes[ip];
   }

   // Distribute level, celltype, H, U, V

   celltype.resize(ncells);
   level.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);
   proc.resize(ncells);

   MPI_Scatterv(&celltype_global[0], &nsizes[0], &ndispl[0], MPI_INT, &celltype[0], nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&level_global[0],    &nsizes[0], &ndispl[0], MPI_INT, &level[0],    nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&i_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &i[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&j_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &j[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&proc_global[0],     &nsizes[0], &ndispl[0], MPI_INT, &proc[0],     nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);

   H.resize(ncells);
   U.resize(ncells);
   V.resize(ncells);

   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);

   MPI_Scatterv(&H_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, &H[0], nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&U_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, &U[0], nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&V_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, &V[0], nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);

   MPI_Scatterv(&x_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, &x[0],  nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&dx_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, &dx[0], nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&y_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, &y[0],  nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&dy_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, &dy[0], nsizes[mype], MPI_C_REAL, 0, MPI_COMM_WORLD);

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state_global->mass_sum(mesh_global, enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   if (mype ==0) {
      printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells_global, H_sum);
   }

   //  Set up grid.

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_mysize(ncells_global);
   set_cell_data(&H_global[0]);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_proc(&mesh_global->proc[0]);
#endif
#ifdef HAVE_MPE
   set_mysize(ncells);
   set_cell_data(&H[0]);
   set_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
   set_cell_proc(&mesh->proc[0]);
#endif

   set_window(mesh_global->xmin, mesh_global->xmax, mesh_global->ymin, mesh_global->ymax);
   set_viewmode(view_mode);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water", mype);
      
   set_circle_radius(circle_radius);
   draw_scene();
   if (verbose) sleep(5);
   sleep(2);

   //  Set flag to show mesh results rather than domain decomposition.
   view_mode = 1;
   
   //  Clear superposition of circle on grid output.
   circle_radius = -1.0;
   
   cpu_timer_start(&tstart);

   set_idle_function(&do_calc);
   start_main_loop();
#else
   cpu_timer_start(&tstart);
   for (int it = 0; it < 10000000; it++) {
      do_calc();
   }
#endif
   
   return 0;
}

static int     ncycle  = 0;
static double  simTime = 0.0;

extern "C" void do_calc(void)
{  double g     = 9.80;
   double sigma = 0.95; 
   int icount, jcount;
   int icount_global, jcount_global;

   //  Initialize state variables for GPU calculation.
   int &mype = mesh->mype;
   int &numpe = mesh->numpe;
   int &noffset = mesh->noffset;

   vector<int>   &nsizes   = mesh_global->nsizes;
   vector<int>   &ndispl   = mesh_global->ndispl;

   vector<int>   &nlft_global     = mesh_global->nlft;
   vector<int>   &nrht_global     = mesh_global->nrht;
   vector<int>   &nbot_global     = mesh_global->nbot;
   vector<int>   &ntop_global     = mesh_global->ntop;

   vector<int>   &nlft     = mesh->nlft;
   vector<int>   &nrht     = mesh->nrht;
   vector<int>   &nbot     = mesh->nbot;
   vector<int>   &ntop     = mesh->ntop;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh_global->ncells;
   size_t &ncells           = mesh->ncells;
   size_t &ncells_ghost     = mesh->ncells_ghost;
   int &cell_handle         = mesh->cell_handle;

   vector<real>  &H_global = state_global->H;

   vector<real>  &H = state->H;
   vector<real>  &U = state->U;
   vector<real>  &V = state->V;

   vector<real>  &x  = mesh->x;
   vector<real>  &dx = mesh->dx;
   vector<real>  &y  = mesh->y;
   vector<real>  &dy = mesh->dy;

   vector<real>  &x_global  = mesh_global->x;
   vector<real>  &dx_global = mesh_global->dx;
   vector<real>  &y_global  = mesh_global->y;
   vector<real>  &dy_global = mesh_global->dy;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   if (DEBUG) {
      //if (mype == 0) mesh->print();

      char filename[10];
      sprintf(filename,"out%1d",mype);
      mesh->fp=fopen(filename,"w");

      //mesh->print_local();
   }

   size_t old_ncells = ncells;
   size_t old_ncells_global = ncells_global;
   size_t new_ncells = 0;
   size_t new_ncells_global = 0;
   double H_sum = -1.0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;
      old_ncells_global = ncells_global;

      //  Calculate the real time step for the current discrete time step.
      deltaT = state->set_timestep(mesh, g, sigma);
      simTime += deltaT;

      //  Compare time step values and pass deltaT in to the kernel.
      if (do_comparison_calc) {
         double deltaT_cpu_global = state_global->set_timestep(mesh_global, g, sigma);

         if (fabs(deltaT - deltaT_cpu_global) > .000001) {
            printf("Error with deltaT calc --- cpu_local %lf cpu_global %lf\n",
               deltaT, deltaT_cpu_global);
         }
      }

      mesh->calc_neighbors_local();

      H.resize(ncells_ghost,0.0);
      U.resize(ncells_ghost,0.0);
      V.resize(ncells_ghost,0.0);
      L7_Update(&H[0], L7_REAL, cell_handle);
      L7_Update(&U[0], L7_REAL, cell_handle);
      L7_Update(&V[0], L7_REAL, cell_handle);

      x.resize(ncells_ghost,0.0);
      dx.resize(ncells_ghost,0.0);
      y.resize(ncells_ghost,0.0);
      dy.resize(ncells_ghost,0.0);
      L7_Update(&x[0], L7_REAL, cell_handle);
      L7_Update(&dx[0], L7_REAL, cell_handle);
      L7_Update(&y[0], L7_REAL, cell_handle);
      L7_Update(&dy[0], L7_REAL, cell_handle);

      if (do_comparison_calc) {
         mesh_global->calc_neighbors();

         // Checking CPU parallel to CPU global
         mesh->compare_neighbors_cpu_local_to_cpu_global(ncells_ghost, ncells_global, mesh_global, &nsizes[0], &ndispl[0]);
      }

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Need ghost cells for this routine
      state->apply_boundary_conditions(mesh);

      if (do_comparison_calc) {
        state_global->apply_boundary_conditions(mesh_global);
      }

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->calc_finite_difference_local(mesh, deltaT);

      if (do_comparison_calc) {
         state_global->calc_finite_difference(mesh_global, deltaT);

         // Compare H gathered to H_global, etc
         state->compare_state_cpu_local_to_cpu_global(state_global, "finite difference", ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells(mesh);

      if (do_comparison_calc) {
         state_global->remove_boundary_cells(mesh_global);
      }
      
      //  Check for NANs.
      for (uint ic=0; ic<ncells; ic++) {
         if (isnan(H[ic]))
         {  printf("Got a NAN on cell %d cycle %d\n",ic,ncycle);
            H[ic]=0.0;
            sleep(100);
            L7_Terminate();
            exit(-1); }
      }  //  Complete NAN check.
      
      mpot.resize(ncells_ghost);
      new_ncells = state->calc_refine_potential(mesh, mpot, icount, jcount);
      nlft.clear();
      nrht.clear();
      nbot.clear();
      ntop.clear();
  
      if (do_comparison_calc) {
         mpot_global.resize(ncells_global);
         new_ncells_global = state_global->calc_refine_potential(mesh_global, mpot_global, icount_global, jcount_global);
         nlft_global.clear();
         nrht_global.clear();
         nbot_global.clear();
         ntop_global.clear();

         int icount_test;
         MPI_Allreduce(&icount, &icount_test, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if (icount_test != icount_global) {
            printf("%d: DEBUG -- icount is %d icount_test %d icount_global is %d\n", mype,icount,icount_test,icount_global);
         }

         // Compare mpot to mpot_global
         mesh->compare_mpot_cpu_local_to_cpu_global(ncells_global, &nsizes[0], &ndispl[0], &mpot[0], &mpot_global[0], ncycle);
      }

      int add_ncells = new_ncells - old_ncells;
      state->rezone_all(mesh, mpot, add_ncells);
      mpot.clear();

      int mesh_local_ncells_global;
      MPI_Allreduce(&ncells, &mesh_local_ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      noffset=0;
      MPI_Scan(&ncells, &noffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      noffset -= ncells;

      mesh->do_load_balance_local(new_ncells, mesh_local_ncells_global, H, U, V);

      MPI_Allgather(&ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);

      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }

      if (do_comparison_calc) {
         int add_ncells_global = new_ncells_global - old_ncells_global;
         //printf("%d: DEBUG add %d new %d old %d\n",mype,add_ncells,new_ncells,old_ncells);
         state_global->rezone_all(mesh_global, mpot_global, add_ncells_global);
         mpot_global.clear();

         //printf("%d: DEBUG ncells is %d new_ncells %d old_ncells %d ncells_global %d\n",mype, ncells, new_ncells, old_ncells, ncells_global);

         // And compare H gathered to H_global, etc
         state->compare_state_cpu_local_to_cpu_global(state_global, "rezone all", ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);

         mesh->compare_indices_cpu_local_to_cpu_global(ncells_global, mesh_global, &nsizes[0], &ndispl[0], ncycle);
      } // do_comparison_calc

      H_sum = -1.0;

      if (do_comparison_calc) {
         H_sum = state->mass_sum(mesh, enhanced_precision_sum);

         double H_sum_global = state_global->mass_sum(mesh_global, enhanced_precision_sum);

         if (fabs(H_sum - H_sum_global) > CONSERVATION_EPS) {
            printf("Error with mass sum calculation -- mass_sum %lf mass_sum_global %lf\n",
                    H_sum, H_sum_global);
         }
      }


// XXX
//      mesh->proc.resize(ncells);
//      if (icount) {
//         vector<int> index(ncells);
//         mesh->partition_cells(numpe, mesh->proc, index, cycle_reorder);
//      }

//      if (do_comparison_calc) {
//         mesh_global->proc.resize(ncells_global);

//         if (icount) {
//            vector<int> index_global(ncells_global);
//            mesh_global->partition_cells(numpe, mesh_global->proc, index_global, cycle_reorder);
            //state->state_reorder(index);
//         }
//      }

   } // End burst loop

// XXX if(ncycle == 1) mesh->m_ncycle = 1;

   if (H_sum < 0) {
      H_sum = state->mass_sum(mesh, enhanced_precision_sum);
   }

// XXX (mesh->m_ncycle)++;

   if (mype == 0){
      printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

#ifdef HAVE_GRAPHICS
   mesh->calc_spatial_coordinates(0);
   if (do_comparison_calc) {
      mesh_global->calc_spatial_coordinates(0);

      mesh->compare_coordinates_cpu_local_to_cpu_global(ncells_global, &nsizes[0], &ndispl[0], &x[0], &dx[0], &y[0], &dy[0], &H[0], &x_global[0], &dx_global[0], &y_global[0], &dy_global[0], &H_global[0], ncycle);

   } else {
#ifdef HAVE_OPENGL
      MPI_Allreduce(&ncells, &ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      H_global.resize(ncells_global);
      MPI_Allgatherv(&x[0],  nsizes[mype], MPI_C_REAL, &x_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
      MPI_Allgatherv(&dx[0], nsizes[mype], MPI_C_REAL, &dx_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
      MPI_Allgatherv(&y[0],  nsizes[mype], MPI_C_REAL, &y_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
      MPI_Allgatherv(&dy[0], nsizes[mype], MPI_C_REAL, &dy_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
      MPI_Allgatherv(&H[0],  nsizes[mype], MPI_C_REAL, &H_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);

      if (view_mode == 0) {
         mesh->proc.resize(ncells);
         for (size_t ii = 0; ii<ncells; ii++){
            mesh->proc[ii] = mesh->mype;
         }
      
         mesh_global->proc.resize(ncells_global);
         MPI_Allgatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &mesh_global->proc[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
      }
#endif
   }

#ifdef HAVE_MPE
   set_mysize(ncells);
   set_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
   set_cell_data(&H[0]);
   set_cell_proc(&mesh->proc[0]);
#endif
#ifdef HAVE_OPENGL
   set_mysize(ncells_global);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_data(&H_global[0]);
   set_cell_proc(&mesh_global->proc[0]);
#endif
   set_viewmode(view_mode);
   set_circle_radius(circle_radius);
   draw_scene();
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      if (do_comparison_calc) {
         state_global->output_timing_info(mesh_global, do_cpu_calc, do_gpu_calc, elapsed_time);
      }
      state->output_timing_info(mesh, do_cpu_calc, do_gpu_calc, elapsed_time);

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
}

