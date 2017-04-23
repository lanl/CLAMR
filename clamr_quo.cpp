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
#include <signal.h>
#include <vector>
#include "input.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
#include "state.h"
#include "l7/l7.h"
#include "timer/timer.h"
#include "memstats/memstats.h"

#include <mpi.h>
#include <omp.h>
#include "display.h"

#ifndef DEBUG
#define DEBUG 0
#endif

static int do_cpu_calc = 1;
static int do_gpu_calc = 0;

typedef unsigned int uint;

#ifdef HAVE_GRAPHICS
static double circle_radius=-1.0;

static int view_mode = 0;

#ifdef FULL_PRECISION
   void (*set_cell_coordinates)(double *, double *, double *, double *) = &set_cell_coordinates_double;
   void (*set_cell_data)(double *) = &set_cell_data_double;
#else
   void (*set_cell_coordinates)(float *, float *, float *, float *) = &set_cell_coordinates_float;
   void (*set_cell_data)(float *) = &set_cell_data_float;
#endif

#endif

bool        verbose,        //  Flag for verbose command-line output; init in input.cpp::parseInput().
            localStencil,   //  Flag for use of local stencil; init in input.cpp::parseInput().
            face_based,     //  Flag for face-based finite difference;
            outline;        //  Flag for drawing outlines of cells; init in input.cpp::parseInput().
int         outputInterval, //  Periodicity of output; init in input.cpp::parseInput().
            enhanced_precision_sum,//  Flag for enhanced precision sum (default true); init in input.cpp::parseInput().
            lttrace_on,     //  Flag to turn on logical time trace package;
            do_quo_setup,   //  Flag to turn on quo dynamic scheduling policies package;
            levmx,          //  Maximum number of refinement levels; init in input.cpp::parseInput().
            nx,             //  x-resolution of coarse grid; init in input.cpp::parseInput().
            ny,             //  y-resolution of coarse grid; init in input.cpp::parseInput().
            niter,          //  Maximum time step; init in input.cpp::parseInput().
            ndim    = 2;    //  Dimensionality of problem (2 or 3).

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
static Mesh       *mesh;     //  Object containing mesh information; init in grid.cpp::main().
static State      *state;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().


//  Set up timing information.
static struct timeval tstart;

static double H_sum_initial = 0.0;
static double cpu_time_graphics = 0.0;
double cpu_time_main_setup = 0.0;
vector<state_t> H_global;
vector<spatial_t> x_global;
vector<spatial_t> dx_global;
vector<spatial_t> y_global;
vector<spatial_t> dy_global;
vector<int> proc_global;

int main(int argc, char **argv) {

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   parseInput(argc, argv);
   L7_Init(&mype, &numpe, &argc, argv);

#if 1 // SKG make things sane for debugging
   signal(SIGSEGV, SIG_DFL);
#endif

   struct timeval tstart_setup;
   cpu_timer_start(&tstart_setup);

   real_t circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (real_t) nx / 128.0;
   int boundary = 1;
   int parallel_in = 1;

   // figure out the max number of threads that can be spawned
   if (0 == mype) {
        int nt = omp_get_max_threads();
        printf("--- num openmp threads: %d\n", nt);
        fflush(stdout);
   }

   mesh = new Mesh(nx, ny, levmx, ndim, boundary, parallel_in, do_gpu_calc);
   if (DEBUG) {
      //if (mype == 0) mesh->print();

      char filename[10];
      sprintf(filename,"out%1d",mype);
      mesh->fp=fopen(filename,"w");

      //mesh->print_local();
   }
   mesh->init(nx, ny, circ_radius, initial_order, do_gpu_calc);

   size_t &ncells = mesh->ncells;
   size_t &ncells_global = mesh->ncells_global;
   int &noffset = mesh->noffset;

   state = new State(mesh);
   state->init(do_gpu_calc);

   vector<int>   &nsizes     = mesh->nsizes;
   vector<int>   &ndispl     = mesh->ndispl;

   vector<spatial_t> &x  = mesh->x;
   vector<spatial_t> &dx = mesh->dx;
   vector<spatial_t> &y  = mesh->y;
   vector<spatial_t> &dy = mesh->dy;

   nsizes.resize(numpe);
   ndispl.resize(numpe);

   int ncells_int = ncells;
   MPI_Allgather(&ncells_int, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);

   ndispl[0]=0;
   for (int ip=1; ip<numpe; ip++){
      ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
   }
   noffset = ndispl[mype];

   state->resize(ncells);

   state->fill_circle(circ_radius, 100.0, 7.0);

   x.clear();
   dx.clear();
   y.clear();
   dy.clear();

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   double cpu_time_main_setup = cpu_timer_stop(tstart_setup);
   mesh->parallel_timer_output("CPU:  setup time               time was",cpu_time_main_setup, 0);

   long long mem_used = memstats_memused();
   if (mem_used > 0) {
      mesh->parallel_memory_output("Memory used      in startup ",mem_used, 0);
      mesh->parallel_memory_output("Memory peak      in startup ",memstats_mempeak(), 0);
      mesh->parallel_memory_output("Memory free      at startup ",memstats_memfree(), 0);
      mesh->parallel_memory_output("Memory available at startup ",memstats_memtotal(), 0);
   }

   if (mype == 0) {
      printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells_global, H_sum);
   }

   for (int i = 0; i < MESH_COUNTER_SIZE; i++){
      mesh->cpu_counters[i]=0;
   }
   for (int i = 0; i < MESH_TIMER_SIZE; i++){
      mesh->cpu_timers[i]=0.0;
   }

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_mysize(ncells_global);
   //vector<state_t> H_global;
   //vector<spatial_t> x_global;
   //vector<spatial_t> dx_global;
   //vector<spatial_t> y_global;
   //vector<spatial_t> dy_global;
   //vector<int> proc_global;
   if (mype == 0){
      H_global.resize(ncells_global);
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      proc_global.resize(ncells_global);
   }
   MPI_Gatherv(&x[0],  nsizes[mype], MPI_SPATIAL_T, &x_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dx[0], nsizes[mype], MPI_SPATIAL_T, &dx_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&y[0],  nsizes[mype], MPI_SPATIAL_T, &y_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dy[0], nsizes[mype], MPI_SPATIAL_T, &dy_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&state->H[0], nsizes[mype], MPI_STATE_T, &H_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, 0, MPI_COMM_WORLD);

   set_cell_data(&H_global[0]);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);

   if (view_mode == 0) {
      mesh->proc.resize(ncells);
      for (size_t ii = 0; ii<ncells; ii++){
         mesh->proc[ii] = mesh->mype;
      }
   
      MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   }

   set_cell_proc(&proc_global[0]);
#endif
#ifdef HAVE_MPE
   set_mysize(ncells);
   set_cell_data(&state->H[0]);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_proc(&mesh->proc[0]);
#endif

   set_window((float)mesh->xmin, (float)mesh->xmax, (float)mesh->ymin, (float)mesh->ymax);
   set_viewmode(view_mode);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water");

   set_circle_radius(circle_radius);
   draw_scene();
   if (verbose) sleep(5);
   sleep(2);

   //  Set flag to show mesh results rather than domain decomposition.
   view_mode = 1;
   
   //  Clear superposition of circle on grid output.
   circle_radius = -1.0;
   
   MPI_Barrier(MPI_COMM_WORLD);
   cpu_timer_start(&tstart);

   set_idle_function(&do_calc);
   start_main_loop();
#else
   MPI_Barrier(MPI_COMM_WORLD);
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
   struct timeval tstart_cpu;

   //  Initialize state variables for GPU calculation.
   int &mype  = mesh->mype;
   int &numpe = mesh->numpe;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh->ncells_global;
   size_t &ncells           = mesh->ncells;
   size_t &ncells_ghost     = mesh->ncells_ghost;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   size_t old_ncells = ncells;
   size_t old_ncells_global = ncells_global;
   size_t new_ncells = 0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;
      old_ncells_global = ncells_global;

      MPI_Barrier(MPI_COMM_WORLD);
      cpu_timer_start(&tstart_cpu);
      //  Calculate the real time step for the current discrete time step.
      deltaT = state->set_timestep(g, sigma);
      simTime += deltaT;

      cpu_timer_start(&tstart_cpu);
      mesh->calc_neighbors_local();

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (mesh->have_boundary) {
      //  state->add_boundary_cells();
      //}

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      cpu_timer_start(&tstart_cpu);
      state->calc_finite_difference(deltaT);

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells();

      cpu_timer_start(&tstart_cpu);
      mpot.resize(ncells_ghost);
      new_ncells = state->calc_refine_potential(mpot, icount, jcount);
  
      cpu_timer_start(&tstart_cpu);
      //int add_ncells = new_ncells - old_ncells;
      state->rezone_all(icount, jcount, mpot);
      // Clear does not delete mpot, so have to swap with an empty vector to get
      // it to delete the mpot memory. This is all to avoid valgrind from showing
      // it as a reachable memory leak
      //mpot.clear();
      vector<int>().swap(mpot);
      ncells = new_ncells;
      mesh->ncells = new_ncells;

      cpu_timer_start(&tstart_cpu);
      state->do_load_balance_local(new_ncells);

// XXX
//      mesh->proc.resize(ncells);
//      if (icount) {
//         vector<int> index(ncells);
//         mesh->partition_cells(numpe, index, cycle_reorder);
//      }

   } // End burst loop

   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (isnan(H_sum)) {
      printf("Got a NAN on cycle %d\n",ncycle);
      exit(-1);
   }
   if (mype == 0){
      printf("Iteration %3d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

#ifdef HAVE_GRAPHICS
   mesh->x.resize(ncells);
   mesh->dx.resize(ncells);
   mesh->y.resize(ncells);
   mesh->dy.resize(ncells);
   mesh->calc_spatial_coordinates(0);

   cpu_timer_start(&tstart_cpu);

#ifdef HAVE_MPE
   set_mysize(ncells);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_data(&state->H[0]);
   set_cell_proc(&mesh->proc[0]);
#endif
#ifdef HAVE_OPENGL
   vector<int>   &nsizes   = mesh->nsizes;
   vector<int>   &ndispl   = mesh->ndispl;

   set_mysize(ncells_global);
   //vector<spatial_t> x_global;
   //vector<spatial_t> dx_global;
   //vector<spatial_t> y_global;
   //vector<spatial_t> dy_global;
   //vector<state_t> H_global;
   //vector<int> proc_global;

   if (mype == 0) {
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      H_global.resize(ncells_global);
      proc_global.resize(ncells_global);
   }
   MPI_Gatherv(&mesh->x[0],  nsizes[mype], MPI_SPATIAL_T, &x_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&mesh->dx[0], nsizes[mype], MPI_SPATIAL_T, &dx_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&mesh->y[0],  nsizes[mype], MPI_SPATIAL_T, &y_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&mesh->dy[0], nsizes[mype], MPI_SPATIAL_T, &dy_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&state->H[0], nsizes[mype], MPI_STATE_T, &H_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, 0, MPI_COMM_WORLD);

   if (view_mode == 0) {
      mesh->proc.resize(ncells);
      for (size_t ii = 0; ii<ncells; ii++){
         mesh->proc[ii] = mesh->mype;
      }
   
      MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_INT, 0, MPI_COMM_WORLD);
   }

   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_data(&H_global[0]);
   set_cell_proc(&proc_global[0]);
#endif
   set_viewmode(view_mode);
   set_circle_radius(circle_radius);
   draw_scene();

   MPI_Barrier(MPI_COMM_WORLD);

   cpu_time_graphics += cpu_timer_stop(tstart_cpu);
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      long long mem_used = memstats_memused();
      if (mem_used > 0) {
         mesh->parallel_memory_output("Memory used      ",mem_used, 0);
         mesh->parallel_memory_output("Memory peak      ",memstats_mempeak(), 0);
         mesh->parallel_memory_output("Memory free      ",memstats_memfree(), 0);
         mesh->parallel_memory_output("Memory available ",memstats_memtotal(), 0);
      }
      state->output_timing_info(do_cpu_calc, do_gpu_calc, elapsed_time);

      mesh->parallel_timer_output("CPU:  graphics                 time was",cpu_time_graphics, 0);

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      if (mype ==0) {
         printf("CPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_cpu_rezone_count()/(double)ncycle*100.0 );
         printf("CPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_cpu_calc_neigh_count()/(double)ncycle*100.0 );
         printf("CPU:  load balance frequency          \t %8.4f\tpercent\n",     (double)mesh->get_cpu_load_balance_count()/(double)ncycle*100.0 );
         printf("CPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_cpu_refine_smooth_count()/(double)mesh->get_cpu_rezone_count() );
      }

      mesh->terminate();
      state->terminate();

      delete mesh;
      delete state;

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
}

