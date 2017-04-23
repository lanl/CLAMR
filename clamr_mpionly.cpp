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
#include "graphics/display.h"
#include "graphics/graphics.h"
#include "input.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
#include "state.h"
#include "l7/l7.h"
#include "timer/timer.h"
#include "memstats/memstats.h"
#include "PowerParser/PowerParser.hh"

using namespace PP;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

static int do_cpu_calc = 1;
static int do_gpu_calc = 0;

typedef unsigned int uint;

static bool do_display_graphics = false;
static bool do_display_opengl_graphics = false;

#ifdef HAVE_GRAPHICS
static double circle_radius=-1.0;
#ifdef FULL_PRECISION
   void (*set_display_cell_coordinates)(double *, double *, double *, double *) = &set_display_cell_coordinates_double;
   void (*set_display_cell_data)(double *) = &set_display_cell_data_double;
#else
   void (*set_display_cell_coordinates)(float *, float *, float *, float *) = &set_display_cell_coordinates_float;
   void (*set_display_cell_data)(float *) = &set_display_cell_data_float;
#endif
#endif

static int view_mode = 0;

#ifdef FULL_PRECISION
#define  SUM_ERROR 2.0e-16
   void (*set_graphics_cell_coordinates)(double *, double *, double *, double *) = &set_graphics_cell_coordinates_double;
   void (*set_graphics_cell_data)(double *) = &set_graphics_cell_data_double;
#else
#define  SUM_ERROR 1.0e-8
   void (*set_graphics_cell_coordinates)(float *, float *, float *, float *) = &set_graphics_cell_coordinates_float;
   void (*set_graphics_cell_data)(float *) = &set_graphics_cell_data_float;
#endif

bool        restart,        //  Flag to start from a back up file; init in input.cpp::parseInput().
            verbose,        //  Flag for verbose command-line output; init in input.cpp::parseInput().
            localStencil,   //  Flag for use of local stencil; init in input.cpp::parseInput().
            face_based,     //  Flag for face-based finite difference;
            outline;        //  Flag for drawing outlines of cells; init in input.cpp::parseInput().
int         outputInterval, //  Periodicity of output; init in input.cpp::parseInput().
            crux_type,      //  Type of checkpoint/restart -- CRUX_NONE, CRUX_IN_MEMORY, CRUX_DISK;
                            //  init in input.cpp::parseInput().
            enhanced_precision_sum,//  Flag for enhanced precision sum (default true); init in input.cpp::parseInput().
            lttrace_on,     //  Flag to turn on logical time trace package;
            do_quo_setup,   //  Flag to turn on quo dynamic scheduling policies package;
            levmx,          //  Maximum number of refinement levels; init in input.cpp::parseInput().
            nx,             //  x-resolution of coarse grid; init in input.cpp::parseInput().
            ny,             //  y-resolution of coarse grid; init in input.cpp::parseInput().
            niter,          //  Maximum iterations; init in input.cpp::parseInput().
            graphic_outputInterval, // Periodicity of graphic output that is saved; init in input.cpp::parseInput()
            checkpoint_outputInterval, // Periodicity of checkpoint output that is saved; init in input.cpp::parseInput()
            num_of_rollback_states,// Maximum number of rollback states to maintain; init in input.cpp::parseInput()
            backup_file_num,//  Backup file number to restart simulation from; init in input.cpp::parseInput()
            numpe,          //  
            ndim    = 2;    //  Dimensionality of problem (2 or 3).
double      upper_mass_diff_percentage; //  Flag for the allowed pecentage difference to the total
                                        //  mass per output intervals; init in input.cpp::parseInput().

char *restart_file;

static int it = 0;

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
static Mesh        *mesh;           //  Object containing mesh information
static State       *state;          //  Object containing state information corresponding to mesh
static PowerParser *parse;          //  Object containing input file parsing

static real_t circ_radius = 0.0;
static int next_graphics_cycle = 0;

//  Set up timing information.
static struct timeval tstart, tstart_cpu, tstart_partmeas;

static double H_sum_initial = 0.0;
static double cpu_time_graphics = 0.0;
static double cpu_time_calcs    = 0.0;
static double cpu_time_partmeas = 0.0;

static int     ncycle  = 0;
static double  simTime = 0.0;
static double  deltaT = 0.0;

vector<state_t> H_global;
vector<spatial_t> x_global;
vector<spatial_t> dx_global;
vector<spatial_t> y_global;
vector<spatial_t> dy_global;
vector<int> proc_global;

int main(int argc, char **argv) {

   // Needed for code to compile correctly on the Mac
   int mype=0;
   int numpe=0;

   //  Process command-line arguments, if any.
   parseInput(argc, argv);
   L7_Init(&mype, &numpe, &argc, argv, do_quo_setup, lttrace_on);

#ifdef _OPENMP
   int nt = 0;
   int tid = 0;

   nt = omp_get_max_threads();
   tid = omp_get_thread_num();
   if (0 == tid && mype == 0) {
        printf("--- max num openmp threads: %d\n", nt);
   }
#pragma omp parallel
   {
      nt = omp_get_num_threads();
      tid = omp_get_thread_num();

#pragma omp master
      if (mype == 0) {
         printf("--- num openmp threads in parallel region: %d\n", nt);
      }
   }
#endif

   parse = new PowerParser();

   struct timeval tstart_setup;
   cpu_timer_start(&tstart_setup);

   circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (real_t) nx / 128.0;
   int boundary = 1;
   int parallel_in = 1;
   double deltax_in = 1.0;
   double deltay_in = 1.0;

   mesh = new Mesh(nx, ny, levmx, ndim, deltax_in, deltay_in, boundary, parallel_in, do_gpu_calc);

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

   if (graphic_outputInterval > niter) next_graphics_cycle = graphic_outputInterval;

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   double cpu_time_main_setup = cpu_timer_stop(tstart_setup);
   mesh->parallel_output("CPU:  setup time               time was",cpu_time_main_setup, 0, "s");

   long long mem_used = memstats_memused();
   if (mem_used > 0) {
      mesh->parallel_output("Memory used      in startup ",mem_used, 0, "kB");
      mesh->parallel_output("Memory peak      in startup ",memstats_mempeak(), 0, "kB");
      mesh->parallel_output("Memory free      at startup ",memstats_memfree(), 0, "kB");
      mesh->parallel_output("Memory available at startup ",memstats_memtotal(), 0, "kB");
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

   cpu_timer_start(&tstart_cpu);

#ifdef HAVE_GRAPHICS
   do_display_graphics = true;
#ifdef HAVE_OPENGL
   do_display_opengl_graphics = true;
#endif
#endif

   if (do_display_graphics || ncycle == next_graphics_cycle){
      mesh->calc_spatial_coordinates(0);
   }

   if (do_display_opengl_graphics || ncycle == next_graphics_cycle){
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

      if (view_mode == 0) {
         mesh->proc.resize(ncells);
         for (size_t ii = 0; ii<ncells; ii++){
            mesh->proc[ii] = mesh->mype;
         }
   
         MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_INT, 0, MPI_COMM_WORLD);
      }
   }

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_display_mysize(ncells_global);
   set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_display_cell_data(&H_global[0]);
   set_display_cell_proc(&proc_global[0]);
#endif
#ifdef HAVE_MPE
   set_display_mysize(ncells);
   set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_display_cell_data(&state->H[0]);
   set_display_cell_proc(&mesh->proc[0]);
#endif

   set_display_window((float)mesh->xmin, (float)mesh->xmax,
                      (float)mesh->ymin, (float)mesh->ymax);
   set_display_outline((int)outline);
   set_display_viewmode(view_mode);
#endif

   if (ncycle == next_graphics_cycle){
      set_graphics_outline(outline);
      set_graphics_window((float)mesh->xmin, (float)mesh->xmax,
                          (float)mesh->ymin, (float)mesh->ymax);
      set_graphics_mysize(ncells_global);
      set_graphics_cell_coordinates(&x_global[0], &dx_global[0],
                                    &y_global[0], &dy_global[0]);
      set_graphics_cell_data(&H_global[0]);
      set_graphics_cell_proc(&proc_global[0]);
      set_graphics_viewmode(view_mode);

      if (mype == 0) {
         init_graphics_output();
         write_graphics_info(0,0,0.0,0,0);
      }
      next_graphics_cycle += graphic_outputInterval;
   }

#ifdef HAVE_GRAPHICS
   set_display_circle_radius(circle_radius);
   init_display(&argc, argv, "Shallow Water");
   draw_scene();
   //if (verbose) sleep(5);
   sleep(2);

   //  Clear superposition of circle on grid output.
   circle_radius = -1.0;
#endif
   cpu_time_graphics += cpu_timer_stop(tstart_cpu);

   //  Set flag to show mesh results rather than domain decomposition.
   view_mode = 1;

   cpu_timer_start(&tstart);
#ifdef HAVE_GRAPHICS
   set_idle_function(&do_calc);
   start_main_loop();
#else
   for (it = 0; it < 10000000; it++) {
      do_calc();
   }
#endif
   
   return 0;
}

extern "C" void do_calc(void)
{  double g     = 9.80;
   double sigma = 0.95;
   int icount, jcount;
   struct timeval tstart_cpu;

   //  Initialize state variables for GPU calculation.
   int &mype  = mesh->mype;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh->ncells_global;
   size_t &ncells           = mesh->ncells;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   size_t new_ncells = 0;

   //  Main loop.
   int endcycle = MIN(niter, next_graphics_cycle);

   cpu_timer_start(&tstart_cpu);

   for (int nburst = ncycle % outputInterval; nburst < outputInterval && ncycle < endcycle; nburst++, ncycle++) {

      //  Calculate the real time step for the current discrete time step.
      deltaT = state->set_timestep(g, sigma);
      simTime += deltaT;

      mesh->calc_neighbors_local();

      cpu_timer_start(&tstart_partmeas);
      mesh->partition_measure();
      cpu_time_partmeas += cpu_timer_stop(tstart_partmeas);

      // Currently not working -- may need to be earlier?
      //if (mesh->have_boundary) {
      //  state->add_boundary_cells();
      //}

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->calc_finite_difference(deltaT);

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells();

      mpot.resize(mesh->ncells_ghost);
      new_ncells = state->calc_refine_potential(mpot, icount, jcount);

      //  Resize the mesh, inserting cells where refinement is necessary.

      state->rezone_all(icount, jcount, mpot);

      // Clear does not delete mpot, so have to swap with an empty vector to get
      // it to delete the mpot memory. This is all to avoid valgrind from showing
      // it as a reachable memory leak
      //mpot.clear();
      vector<int>().swap(mpot);

      mesh->ncells = new_ncells;
      ncells = new_ncells;

      state->do_load_balance_local(new_ncells);

// XXX
//      mesh->proc.resize(ncells);
//      if (icount) {
//         vector<int> index(ncells);
//         mesh->partition_cells(numpe, index, cycle_reorder);
//      }

   } // End burst loop

   cpu_time_calcs += cpu_timer_stop(tstart_cpu);

   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (isnan(H_sum)) {
      printf("Got a NAN on cycle %d\n",ncycle);
      exit(-1);
   }

   if (mype == 0){
      printf("Iteration %3d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

   cpu_timer_start(&tstart_cpu);

   if(do_display_graphics || ncycle == next_graphics_cycle ||
      (ncycle >= niter && graphic_outputInterval < niter) ){

      mesh->calc_spatial_coordinates(0);
   }

   if (do_display_opengl_graphics || ncycle == next_graphics_cycle){
      vector<int>   &nsizes   = mesh->nsizes;
      vector<int>   &ndispl   = mesh->ndispl;

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
   }

   if (ncycle == next_graphics_cycle){
      set_graphics_mysize(ncells_global);
      set_graphics_viewmode(view_mode);
      set_graphics_cell_coordinates(&x_global[0], &dx_global[0],
                                    &y_global[0], &dy_global[0]);
      set_graphics_cell_data(&H_global[0]);
      set_graphics_cell_proc(&proc_global[0]);

      if (mype == 0) {
         write_graphics_info(0,0,0.0,0,0);
      }
      next_graphics_cycle += graphic_outputInterval;
   }

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_display_mysize(ncells_global);
   set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_display_cell_data(&H_global[0]);
   set_display_cell_proc(&proc_global[0]);
#endif
#ifdef HAVE_MPE
   set_display_mysize(ncells);
   set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_display_cell_data(&state->H[0]);
   set_display_cell_proc(&mesh->proc[0]);
#endif

   set_display_circle_radius(circle_radius);
   set_display_viewmode(view_mode);
   draw_scene();
#endif

   cpu_time_graphics += cpu_timer_stop(tstart_cpu);

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      if(graphic_outputInterval < niter){
         cpu_timer_start(&tstart_cpu);

#ifdef HAVE_GRAPHICS
         set_display_viewmode(view_mode);
#ifdef HAVE_OPENGL
         set_display_mysize(ncells_global);
         set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
         set_display_cell_data(&H_global[0]);
         set_display_cell_proc(&proc_global[0]);
#endif
#ifdef HAVE_MPE
         set_display_mysize(ncells);
         set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
         set_display_cell_data(&state->H[0]);
         set_display_cell_proc(&mesh->proc[0]);
#endif
#endif

         if (mype == 0) {
            write_graphics_info(ncycle/graphic_outputInterval,ncycle,simTime,0,0);
         }
         next_graphics_cycle += graphic_outputInterval;

         cpu_time_graphics += cpu_timer_stop(tstart_cpu);
      }

      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      long long mem_used = memstats_memused();
      if (mem_used > 0) {
         mesh->parallel_output("Memory used      ",mem_used, 0, "kB");
         mesh->parallel_output("Memory peak      ",memstats_mempeak(), 0, "kB");
         mesh->parallel_output("Memory free      ",memstats_memfree(), 0, "kB");
         mesh->parallel_output("Memory available ",memstats_memtotal(), 0, "kB");
      }
      state->output_timing_info(do_cpu_calc, do_gpu_calc, elapsed_time);
      mesh->parallel_output("CPU:  calc incl part meas     time was",cpu_time_calcs,    0, "s");
      mesh->parallel_output("CPU:  calculation only        time was",cpu_time_calcs-cpu_time_partmeas,    0, "s");
      mesh->parallel_output("CPU:  partition measure       time was",cpu_time_partmeas, 0, "s");
      mesh->parallel_output("CPU:  graphics                time was",cpu_time_graphics, 0, "s");

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      if (mype ==0) {
         printf("CPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_REZONE)/(double)ncycle*100.0 );
         printf("CPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_CALC_NEIGH)/(double)ncycle*100.0 );
         printf("CPU:  load balance frequency          \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_LOAD_BALANCE)/(double)ncycle*100.0 );
         printf("CPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_cpu_counter(MESH_COUNTER_REFINE_SMOOTH)/(double)mesh->get_cpu_counter(MESH_COUNTER_REZONE) );
      }

      mesh->terminate();
      state->terminate();

      terminate_graphics_output();

      delete mesh;
      delete state;
      delete parse;

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
} // end do_calc

