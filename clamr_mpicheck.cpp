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

#include "mesh/mesh.h"
#include "mesh/partition.h"

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>
#include "input.h"
#include "state.h"
#include "l7/l7.h"
#include "timer/timer.h"

#include "graphics/display.h"

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

typedef unsigned int uint;

#ifdef HAVE_GRAPHICS
static double circle_radius=-1.0;

static int view_mode = 0;

#ifdef FULL_PRECISION
   void (*set_display_cell_coordinates)(double *, double *, double *, double *) = &set_display_cell_coordinates_double;
   void (*set_display_cell_data)(double *) = &set_display_cell_data_double;
#else
   void (*set_display_cell_coordinates)(float *, float *, float *, float *) = &set_display_cell_coordinates_float;
   void (*set_display_cell_data)(float *) = &set_display_cell_data_float;
#endif

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
            niter,          //  Maximum time step; init in input.cpp::parseInput().
            graphic_outputInterval, // Periocity of graphic output that is saved; init in input.cpp::parseInput()
            checkpoint_outputInterval, // Periodicity of checkpoint output that is saved; init in input.cpp::parseInput()
            num_of_rollback_states,// Maximum number of rollback states to maintain; init in input.cpp::parseInput()
            backup_file_num,//  Backup file number to restart simulation from; init in input.cpp::parseInput()
            numpe,          //  
            ndim    = 2;    //  Dimensionality of problem (2 or 3).
double      upper_mass_diff_percentage; //  Flag for the allowed pecentage difference to the total
                                        //  mass per output intervals; init in input.cpp::parseInput().

char *restart_file;

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
static Mesh       *mesh_global;    //  Object containing mesh information; init in grid.cpp::main().
static State      *state_global;   //  Object containing state information corresponding to mesh; init in grid.cpp::main().
static Mesh       *mesh;           //  Object containing mesh information; init in grid.cpp::main().
static State      *state;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
static struct timeval tstart;

static double  H_sum_initial = 0.0;
static double cpu_time_graphics = 0.0;

int main(int argc, char **argv) {

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   parseInput(argc, argv);
   L7_Init(&mype, &numpe, &argc, argv, do_quo_setup, lttrace_on);

   real_t circ_radius = 6.0;
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
   int &noffset = mesh->noffset;

   int parallel_global_in = 0;
   mesh_global  = new Mesh(nx, ny, levmx, ndim, deltax_in, deltay_in, boundary, parallel_global_in, do_gpu_calc);

   size_t &ncells_global = mesh_global->ncells;
   MPI_Allreduce(&ncells, &ncells_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

   state = new State(mesh);
   state->init(do_gpu_calc);

   state_global = new State(mesh_global);
   state_global->allocate(ncells_global);
   state_t  *H_global = state_global->H;
   state_t  *U_global = state_global->U;
   state_t  *V_global = state_global->V;

   vector<int>   &nsizes     = mesh->nsizes;
   vector<int>   &ndispl     = mesh->ndispl;

   vector<spatial_t>  &x_global  = mesh_global->x;
   vector<spatial_t>  &dx_global = mesh_global->dx;
   vector<spatial_t>  &y_global  = mesh_global->y;
   vector<spatial_t>  &dy_global = mesh_global->dy;

   vector<int>   &proc     = mesh->proc;

   vector<int>   &proc_global     = mesh_global->proc;

   vector<spatial_t> &x  = mesh->x;
   vector<spatial_t> &dx = mesh->dx;
   vector<spatial_t> &y  = mesh->y;
   vector<spatial_t> &dy = mesh->dy;

   nsizes.resize(numpe);
   ndispl.resize(numpe);
   MPI_Allgather(&ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
   ndispl[0]=0;
   for (int ip=1; ip<numpe; ip++){
      ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
   }
   noffset=ndispl[mype];

   // Gather level, celltype, H, U, V for global calc

   mesh_global->celltype = (int *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(int), "celltype");
   mesh_global->level    = (int *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(int), "level");
   mesh_global->i        = (int *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(int), "i");
   mesh_global->j        = (int *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(int), "j");

   proc_global.resize(ncells_global);

   MPI_Allgatherv(&mesh->celltype[0], ncells, MPI_INT, &mesh_global->celltype[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&mesh->level[0], ncells, MPI_INT, &mesh_global->level[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&mesh->i[0], ncells, MPI_INT, &mesh_global->i[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&mesh->j[0], ncells, MPI_INT, &mesh_global->j[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

   proc.resize(ncells);
   for (uint ic=0; ic<ncells; ic++){
      proc[ic] = mesh->mype;
   }
   MPI_Allgatherv(&proc[0], ncells, MPI_INT, &proc_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

   state->resize(ncells);

   x_global.resize(ncells_global);
   dx_global.resize(ncells_global);
   y_global.resize(ncells_global);
   dy_global.resize(ncells_global);

   MPI_Allgatherv(&x[0],  ncells, MPI_SPATIAL_T, &x_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&dx[0], ncells, MPI_SPATIAL_T, &dx_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&y[0],  ncells, MPI_SPATIAL_T, &y_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&dy[0], ncells, MPI_SPATIAL_T, &dy_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);

   state->fill_circle(circ_radius, 100.0, 7.0);

   MPI_Allgatherv(&state->H[0], nsizes[mype], MPI_STATE_T, &H_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&state->U[0], nsizes[mype], MPI_STATE_T, &U_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&state->V[0], nsizes[mype], MPI_STATE_T, &V_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   if (mype ==0) {
      printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells_global, H_sum);
   }

   for (int i = 0; i < MESH_COUNTER_SIZE; i++){
      mesh_global->cpu_counters[i]=0;
   }
   for (int i = 0; i < MESH_TIMER_SIZE; i++){
      mesh_global->cpu_timers[i]=0.0;
   }   

   //  Set up grid.

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_display_mysize(ncells_global);
   set_display_cell_data(&H_global[0]);
   set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_display_cell_proc(&mesh_global->proc[0]);
#endif
#ifdef HAVE_MPE
   set_display_mysize(ncells);
   set_display_cell_data(&state->H[0]);
   set_display_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
   set_display_cell_proc(&mesh->proc[0]);
#endif

   set_display_window((float)mesh_global->xmin, (float)mesh_global->xmax, (float)mesh_global->ymin, (float)mesh_global->ymax);
   set_display_viewmode(view_mode);
   set_display_outline((int)outline);
   init_display(&argc, argv, "Shallow Water");
      
   set_display_circle_radius(circle_radius);
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
   int icount_global, jcount_global;
#ifdef HAVE_GRAPHICS
   struct timeval tstart_cpu;
#endif

   //  Initialize state variables for GPU calculation.
   int &mype = mesh->mype;

   vector<int>   &nsizes   = mesh->nsizes;
   vector<int>   &ndispl   = mesh->ndispl;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh_global->ncells;
   size_t &ncells           = mesh->ncells;
   size_t &ncells_ghost     = mesh->ncells_ghost;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   //size_t old_ncells = ncells;
   //size_t old_ncells_global = ncells_global;
   size_t new_ncells = 0;
   size_t new_ncells_global = 0;
   double H_sum = -1.0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      //  Define basic domain decomposition parameters for GPU.
      //old_ncells = ncells;
      //old_ncells_global = ncells_global;

      //  Calculate the real time step for the current discrete time step.
      deltaT = state->set_timestep(g, sigma);
      simTime += deltaT;

      //  Compare time step values and pass deltaT in to the kernel.
      if (do_comparison_calc) {
         double deltaT_cpu_global = state_global->set_timestep(g, sigma);

         if (fabs(deltaT - deltaT_cpu_global) > .000001) {
            printf("Error with deltaT calc --- cpu_local %lf cpu_global %lf\n",
               deltaT, deltaT_cpu_global);
         }
      }

      mesh->calc_neighbors_local();

      if (do_comparison_calc) {
         mesh_global->calc_neighbors(mesh_global->ncells);

         // Checking CPU parallel to CPU global
         mesh->compare_neighbors_cpu_local_to_cpu_global(ncells_ghost, ncells_global, mesh_global, &nsizes[0], &ndispl[0]);
      }

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->calc_finite_difference(deltaT);

      if (do_comparison_calc) {
         state_global->calc_finite_difference(deltaT);

         // Compare H gathered to H_global, etc
         state->compare_state_cpu_local_to_cpu_global(state_global, "finite difference", ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells();

      if (do_comparison_calc) {
         state_global->remove_boundary_cells();
      }
      
      mpot.resize(ncells_ghost);
      new_ncells = state->calc_refine_potential(mpot, icount, jcount);
  
      if (do_comparison_calc) {
         mpot_global.resize(ncells_global);
         new_ncells_global = state_global->calc_refine_potential(mpot_global, icount_global, jcount_global);

         int icount_test;
         MPI_Allreduce(&icount, &icount_test, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if (icount_test != icount_global) {
            printf("%d: DEBUG -- icount is %d icount_test %d icount_global is %d\n", mype,icount,icount_test,icount_global);
         }

         // Compare mpot to mpot_global
         mesh->compare_mpot_cpu_local_to_cpu_global(ncells_global, &nsizes[0], &ndispl[0], &mpot[0], &mpot_global[0], ncycle);
      }

      //int add_ncells = new_ncells - old_ncells;
      state->rezone_all(icount, jcount, mpot);
      // Clear does not delete mpot, so have to swap with an empty vector to get
      // it to delete the mpot memory. This is all to avoid valgrind from showing
      // it as a reachable memory leak
      //mpot.clear();
      vector<int>().swap(mpot);
      ncells = new_ncells;
      mesh->ncells = new_ncells;

      if (do_comparison_calc) {
         //int add_ncells_global = new_ncells_global - old_ncells_global;
         //printf("%d: DEBUG add %d new %d old %d icount %d jcount %d\n",mype,add_ncells,new_ncells,old_ncells,icount,jcount);
         state_global->rezone_all(icount_global, jcount_global, mpot_global);
         mpot_global.clear();

         ncells_global = new_ncells_global;
         mesh_global->ncells = new_ncells_global;

         //printf("%d: DEBUG ncells is %d new_ncells %d old_ncells %d ncells_global %d\n",mype, ncells, new_ncells, old_ncells, ncells_global);

         // And compare H gathered to H_global, etc
         state->compare_state_cpu_local_to_cpu_global(state_global, "rezone all", ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);

         mesh->compare_indices_cpu_local_to_cpu_global(ncells_global, mesh_global, &nsizes[0], &ndispl[0], ncycle);
      } // do_comparison_calc

      state->do_load_balance_local(new_ncells);

      if (do_comparison_calc) {
         // And compare H gathered to H_global, etc
         state->compare_state_cpu_local_to_cpu_global(state_global, "rezone all", ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);
         mesh->compare_indices_cpu_local_to_cpu_global(ncells_global, mesh_global, &nsizes[0], &ndispl[0], ncycle);
      }

      H_sum = -1.0;

      if (do_comparison_calc) {
         H_sum = state->mass_sum(enhanced_precision_sum);

         double H_sum_global = state_global->mass_sum(enhanced_precision_sum);

         if (fabs(H_sum - H_sum_global) > CONSERVATION_EPS) {
            printf("Error with mass sum calculation -- mass_sum %lf mass_sum_global %lf\n",
                    H_sum, H_sum_global);
         }
      }


// XXX
//      mesh->proc.resize(ncells);
//      if (icount) {
//         vector<int> index(ncells);
//         mesh->partition_cells(numpe, index, cycle_reorder);
//      }

//      if (do_comparison_calc) {
//         mesh_global->proc.resize(ncells_global);

//         if (icount) {
//            vector<int> index_global(ncells_global);
//            mesh_global->partition_cells(numpe, index_global, cycle_reorder);
            //state->state_reorder(index);
//         }
//      }

   } // End burst loop


   if (H_sum < 0) {
      H_sum = state->mass_sum(enhanced_precision_sum);
   }
   if (isnan(H_sum)) {
      printf("Got a NAN on cycle %d\n",ncycle);
      exit(-1);
   }
   if (mype == 0){
      printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

#ifdef HAVE_GRAPHICS
   cpu_timer_start(&tstart_cpu);

   vector<spatial_t>  &x  = mesh->x;
   vector<spatial_t>  &dx = mesh->dx;
   vector<spatial_t>  &y  = mesh->y;
   vector<spatial_t>  &dy = mesh->dy;

   vector<spatial_t>  &x_global  = mesh_global->x;
   vector<spatial_t>  &dx_global = mesh_global->dx;
   vector<spatial_t>  &y_global  = mesh_global->y;
   vector<spatial_t>  &dy_global = mesh_global->dy;

   mesh->calc_spatial_coordinates(0);

   if (do_comparison_calc) {
      mesh_global->calc_spatial_coordinates(0);
#ifdef FULL_PRECISION
      mesh->compare_coordinates_cpu_local_to_cpu_global_double(ncells_global, &nsizes[0], &ndispl[0], &x[0], &dx[0], &y[0], &dy[0], &state->H[0], &x_global[0], &dx_global[0], &y_global[0], &dy_global[0], &state_global->H[0], ncycle);
#else
      mesh->compare_coordinates_cpu_local_to_cpu_global_float(ncells_global, &nsizes[0], &ndispl[0], &x[0], &dx[0], &y[0], &dy[0], &state->H[0], &x_global[0], &dx_global[0], &y_global[0], &dy_global[0], &state_global->H[0], ncycle);
#endif 
   }

#ifdef HAVE_MPE
   set_display_mysize(ncells);
   set_display_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
   set_display_cell_data(&state->H[0]);
   set_display_cell_proc(&mesh->proc[0]);
#endif
#ifdef HAVE_OPENGL
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      MPI_Allgatherv(&x[0],  nsizes[mype], MPI_SPATIAL_T, &x_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
      MPI_Allgatherv(&dx[0], nsizes[mype], MPI_SPATIAL_T, &dx_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
      MPI_Allgatherv(&y[0],  nsizes[mype], MPI_SPATIAL_T, &y_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
      MPI_Allgatherv(&dy[0], nsizes[mype], MPI_SPATIAL_T, &dy_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
      MPI_Allgatherv(&state->H[0],  nsizes[mype], MPI_STATE_T, &state_global->H[0],  &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);

      if (view_mode == 0) {
         mesh->proc.resize(ncells);
         for (size_t ii = 0; ii<ncells; ii++){
            mesh->proc[ii] = mesh->mype;
         }
      
         mesh_global->proc.resize(ncells_global);
         MPI_Allgatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &mesh_global->proc[0],  &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
      }
   set_display_mysize(ncells_global);
   set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_display_cell_data(&state_global->H[0]);
   set_display_cell_proc(&mesh_global->proc[0]);
#endif

   set_display_viewmode(view_mode);
   set_display_circle_radius(circle_radius);
   draw_scene();

   MPI_Barrier(MPI_COMM_WORLD);

   cpu_time_graphics += cpu_timer_stop(tstart_cpu);
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      if (do_comparison_calc) {
         state_global->output_timing_info(do_cpu_calc, do_gpu_calc, elapsed_time);
      }
      state->output_timing_info(do_cpu_calc, do_gpu_calc, elapsed_time);

      mesh->parallel_output("CPU:  graphics                 time was",cpu_time_graphics, 0, "s");

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      if (mype == 0) {
         printf("CPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_REZONE)/(double)ncycle*100.0 );
         printf("CPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_CALC_NEIGH)/(double)ncycle*100.0 );
         printf("CPU:  calc load balance               \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_LOAD_BALANCE)/(double)ncycle*100.0 );
         printf("CPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_cpu_counter(MESH_COUNTER_REFINE_SMOOTH)/(double)mesh->get_cpu_counter(MESH_COUNTER_REZONE) );
         printf("CPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh_global->get_cpu_counter(MESH_COUNTER_REZONE)/(double)ncycle*100.0 );
         printf("CPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh_global->get_cpu_counter(MESH_COUNTER_CALC_NEIGH)/(double)ncycle*100.0 );
         printf("CPU:  calc load balance               \t %8.4f\tpercent\n",     (double)mesh_global->get_cpu_counter(MESH_COUNTER_LOAD_BALANCE)/(double)ncycle*100.0 );
         printf("CPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh_global->get_cpu_counter(MESH_COUNTER_REFINE_SMOOTH)/(double)mesh->get_cpu_counter(MESH_COUNTER_REZONE) );
      }

      mesh->terminate();
      state->terminate();

      delete mesh;
      delete state;

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
}

