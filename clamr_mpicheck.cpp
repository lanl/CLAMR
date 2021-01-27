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
#include "crux/crux.h"
#include "PowerParser/PowerParser.hh"
#include "MallocPlus/MallocPlus.h"

using namespace PP;


#ifndef DEBUG
#define DEBUG 0
#endif
#undef DEBUG_RESTORE_VALS

#define DO_COMPARISON

#ifdef DO_COMPARISON
static int do_comparison_calc = 1;
#else
static int do_comparison_calc = 0;
#endif

static int do_cpu_calc = 1;
static int do_gpu_calc = 0;

extern int choose_amr_method;

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

void store_crux_data(Crux *crux, int ncycle);
void restore_crux_data_bootstrap(Crux *crux, char *restart_file, int rollback_counter);
void restore_crux_data(Crux *crux);

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
            output_cuts,    //  Flag for outputting file of slice along y-axis; init in input.cpp::parseInput().
            backup_file_num,//  Backup file number to restart simulation from; init in input.cpp::parseInput()
            numpe,          //  
            ndim    = 2,    //  Dimensionality of problem (2 or 3).
            ndigits,
            nbits;
double      upper_mass_diff_percentage; //  Flag for the allowed pecentage difference to the total
                                        //  mass per output intervals; init in input.cpp::parseInput().

char *restart_file;

static int it = 0;

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
static Mesh       *mesh_global;     //  Object containing mesh information; init in grid.cpp::main().
static State      *state_global;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().
static Mesh       *mesh;            //  Object containing mesh information; init in grid.cpp::main().
static State      *state;           //  Object containing state information corresponding to mesh; init in grid.cpp::main().
static Crux        *crux;           //  Object containing checkpoint/restart information
static PowerParser *parse;          //  Object containing input file parsing

static real_t circ_radius = 0.0;
static int next_cp_cycle = 0;
static int next_graphics_cycle = 0;

//  Set up timing information.
static struct timespec tstart;

static double H_sum_initial = 0.0;
static double cpu_time_graphics = 0.0;

static int     ncycle  = 0;
static double  simTime = 0.0;
static double  deltaT = 0.0;
char total_sim_time_log[] = {"total_execution_time.log"};
struct timespec total_exec;


int main(int argc, char **argv) {

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   parseInput(argc, argv);
   L7_Init(&mype, &numpe, &argc, argv, do_quo_setup, lttrace_on);

   real_t circ_radius = 6.0;

   parse = new PowerParser();
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

   mesh_global->celltype = (char_t *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(char_t), "celltype");
   mesh_global->level    = (uchar_t *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(uchar_t), "level");
   mesh_global->i        = (int *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(int), "i");
   mesh_global->j        = (int *)mesh_global->mesh_memory.memory_malloc(ncells_global, sizeof(int), "j");

   proc_global.resize(ncells_global);

   MPI_Allgatherv(&mesh->celltype[0], ncells, MPI_CHAR_T, &mesh_global->celltype[0], &nsizes[0], &ndispl[0], MPI_CHAR_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&mesh->level[0], ncells, MPI_UCHAR_T, &mesh_global->level[0], &nsizes[0], &ndispl[0], MPI_UCHAR_T, MPI_COMM_WORLD);
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

   state->fill_circle(circ_radius, 80.0, 10.0);

   MPI_Allgatherv(&state->H[0], nsizes[mype], MPI_STATE_T, &H_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&state->U[0], nsizes[mype], MPI_STATE_T, &U_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&state->V[0], nsizes[mype], MPI_STATE_T, &V_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);

   //  Kahan-type enhanced precision sum implementation.
   mesh->calc_celltype(ncells);
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
#endif
#ifdef HAVE_GRAPHICS
   set_display_circle_radius(circle_radius);
   draw_scene();
   if (verbose) sleep(5);
   sleep(2);

   //  Set flag to show mesh results rather than domain decomposition.
   view_mode = 1;
   
   //  Clear superposition of circle on grid output.
   circle_radius = -1.0;
#endif

#ifdef _OPENMP
#pragma omp parallel
   {
#endif
      mesh->calc_neighbors_local();
      if (do_comparison_calc) {
         mesh_global->calc_neighbors(mesh_global->ncells);
      }
#ifdef _OPENMP
   } // end parallel region
#endif

   cpu_timer_start(&tstart);
#ifdef HAVE_GRAPHICS
   set_idle_function(&do_calc);
   start_main_loop();
#else
   cpu_timer_start(&tstart);
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
   int icount_global, jcount_global;
#ifdef HAVE_GRAPHICS
   struct timespec tstart_cpu;
#endif

   //  Initialize state variables for GPU calculation.
   int &mype  = mesh->mype;

   vector<int>   &nsizes   = mesh->nsizes;
   vector<int>   &ndispl   = mesh->ndispl;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh_global->ncells;
   size_t &ncells           = mesh->ncells;

   vector<char_t>     mpot;
   vector<char_t>     mpot_global;
   
   //size_t old_ncells = ncells;
   //size_t old_ncells_global = ncells_global;
   size_t new_ncells = 0;
   size_t new_ncells_global = 0;
   double H_sum = -1.0;
   double deltaT = 0.0;

   //  Main loop.



   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      mpot.resize(mesh->ncells_ghost);
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
#ifdef _OPENMP
#pragma omp parallel
      {
#endif
         state->rezone_all(icount, jcount, mpot);

         // Clear does not delete mpot, so have to swap with an empty vector to get
         // it to delete the mpot memory. This is all to avoid valgrind from showing
         // it as a reachable memory leak
#ifdef _OPENMP
#pragma omp master
         {
#endif
            //mpot.clear();
            vector<char_t>().swap(mpot);

            mesh->ncells = new_ncells;
            ncells = new_ncells;
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

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

         mesh->set_bounds(ncells);

         state->do_load_balance_local(new_ncells);

      if (do_comparison_calc) {
         // And compare H gathered to H_global, etc
         state->compare_state_cpu_local_to_cpu_global(state_global, "load balance", ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);
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

//          mesh->proc.resize(ncells);
//          if (icount)
//          {  vector<int> index(ncells);
//             mesh->partition_cells(numpe, index, cycle_reorder);
//          }

//          if (do_comparison_calc) {
//             mesh_global->proc.resize(ncells_global);

//             if (icount) {
//                vector<int> index_global(ncells_global);
//                mesh_global->partition_cells(numpe, index_global, cycle_reorder);
                //state->state_reorder(index);
//             }
//          }

         //  Calculate the real time step for the current discrete time step.
         double mydeltaT = state->set_timestep(g, sigma); // Private variable to avoid write conflict
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            deltaT = mydeltaT;
            simTime += deltaT;
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         //  Compare time step values and pass deltaT in to the kernel.
         if (do_comparison_calc) {
            double mydeltaT_cpu_global = state_global->set_timestep(g, sigma);

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
#endif
            double deltaT_cpu_global = mydeltaT_cpu_global;
            if (fabs(deltaT - deltaT_cpu_global) > .000001) {
               printf("Error with deltaT calc --- cpu_local %lf cpu_global %lf\n",
                  deltaT, deltaT_cpu_global);
            }
         }

         mesh->calc_neighbors_local();

         if (do_comparison_calc) {
            mesh_global->calc_neighbors(mesh_global->ncells);

            // Checking CPU parallel to CPU global
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
#endif
            {
               mesh->compare_neighbors_cpu_local_to_cpu_global(mesh->ncells_ghost, ncells_global, mesh_global, &nsizes[0], &ndispl[0]);
            }
         }

         mesh->partition_measure();

         // Currently not working -- may need to be earlier?
         //if (mesh->have_boundary) {
         //  state->add_boundary_cells(mesh);
         //}

         // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

         mesh->set_bounds(ncells);

         //  Execute main kernel
         if (choose_amr_method == CELL_AMR) {
            state->calc_finite_difference(deltaT);
         } else if (choose_amr_method == FACE_AMR) {
            //printf("ERROR -- the face method currently does not work with MPI\n");
            //L7_Terminate();
            //exit(-1);
            state->calc_finite_difference_via_faces(deltaT);
         } else if (choose_amr_method == CELL_IN_PLACE_AMR) {
            printf("ERROR -- the cell-in-place method currently does not work with MPI\n");
            //L7_Terminate();
            //exit(-1);
            state->calc_finite_difference_cell_in_place(deltaT);
         } else if (choose_amr_method == FACE_IN_PLACE_AMR) {
            printf("ERROR -- the face-in-place method currently does not work with MPI\n");
            L7_Terminate();
            exit(-1);
            //state->calc_finite_difference_face_in_place(deltaT);
         } else if (choose_amr_method == REGULAR_GRID_AMR) {
            printf("ERROR -- the regular grid AMR method currently does not work with MPI\n");
            L7_Terminate();
            exit(-1);
            //state->calc_finite_difference_regular_cells(deltaT);
         } else if (choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
            printf("ERROR -- the regular grid by faces AMR method currently does not work with MPI\n");
            L7_Terminate();
            exit(-1);
            //state->calc_finite_difference_regular_cells_by_faces(deltaT);
         } else {
            state->calc_finite_difference(deltaT);
         }

         if (do_comparison_calc) {
            if (choose_amr_method == CELL_AMR) {
               state_global->calc_finite_difference(deltaT);
            } else if (choose_amr_method == FACE_AMR) {
               //printf("ERROR -- the face method currently does not work with MPI\n");
               //L7_Terminate();
               //exit(-1);
               state_global->calc_finite_difference_via_faces(deltaT);
            } else if (choose_amr_method == CELL_IN_PLACE_AMR) {
               printf("ERROR -- the cell-in-place method currently does not work with MPI\n");
               //L7_Terminate();
               //exit(-1);
               state_global->calc_finite_difference_cell_in_place(deltaT);
            } else if (choose_amr_method == FACE_IN_PLACE_AMR) {
               printf("ERROR -- the face-in-place method currently does not work with MPI\n");
               L7_Terminate();
               exit(-1);
               //state_global->calc_finite_difference_face_in_place(deltaT);
            } else if (choose_amr_method == REGULAR_GRID_AMR) {
               printf("ERROR -- the regular grid AMR method currently does not work with MPI\n");
               L7_Terminate();
               exit(-1);
               //state_global->calc_finite_difference_regular_cells(deltaT);
            } else if (choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
               printf("ERROR -- the regular grid by faces AMR method currently does not work with MPI\n");
               L7_Terminate();
               exit(-1);
               //state_global->calc_finite_difference_regular_cells_by_faces(deltaT);
            } else {
               state_global->calc_finite_difference(deltaT);
            }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               // Compare H gathered to H_global, etc
               state->compare_state_cpu_local_to_cpu_global(state_global, "finite difference",
                      ncycle, ncells, ncells_global, &nsizes[0], &ndispl[0]);
#ifdef _OPENMP
            } // end master region
#endif
         }

         //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
         state->remove_boundary_cells();
#ifdef _OPENMP
      } // end parallel region
#endif

      if (do_comparison_calc) {
         state_global->remove_boundary_cells();
      }
      
   } // End burst loop


   if (H_sum < 0) {
      H_sum = state->mass_sum(enhanced_precision_sum);
   }
#ifdef __APPLE__
   if (isnan(H_sum)) {
#else
   if (std::isnan(H_sum)) {
#endif
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
      delete crux;
      delete parse;

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
} // end do_calc

const int CRUX_CLAMR_VERSION = 101;
const int num_int_vals       = 14;
const int num_double_vals    =  5;

MallocPlus clamr_bootstrap_memory;

void store_crux_data(Crux *crux, int ncycle)
{
   size_t nsize = num_int_vals*sizeof(int) +
                  num_double_vals*sizeof(double);
   nsize += state->get_checkpoint_size();

   int int_vals[num_int_vals];

   int_vals[ 0] = CRUX_CLAMR_VERSION; // Version number
   int_vals[ 1] = nx;
   int_vals[ 2] = ny;
   int_vals[ 3] = levmx;
   int_vals[ 4] = ndim;
   int_vals[ 5] = outputInterval;
   int_vals[ 6] = enhanced_precision_sum;
   int_vals[ 7] = niter;
   int_vals[ 8] = it;
   int_vals[ 9] = ncycle;
   int_vals[10] = graphic_outputInterval;
   int_vals[11] = checkpoint_outputInterval;
   int_vals[12] = next_cp_cycle;
   int_vals[13] = next_graphics_cycle;

   double double_vals[num_double_vals];
   double_vals[ 0] = circ_radius;
   double_vals[ 1] = H_sum_initial;
   double_vals[ 2] = simTime;
   double_vals[ 3] = deltaT;
   double_vals[ 4] = upper_mass_diff_percentage;

   clamr_bootstrap_memory.memory_add(int_vals, size_t(num_int_vals), 4, "bootstrap_int_vals", RESTART_DATA);
   clamr_bootstrap_memory.memory_add(double_vals, size_t(num_double_vals), 8, "bootstrap_double_vals", RESTART_DATA);

   crux->store_begin(nsize, ncycle);

   crux->store_MallocPlus(clamr_bootstrap_memory);

   state->store_checkpoint(crux);

   crux->store_end();

   clamr_bootstrap_memory.memory_remove(int_vals);
   clamr_bootstrap_memory.memory_remove(double_vals);

   next_cp_cycle += checkpoint_outputInterval;
}

void restore_crux_data_bootstrap(Crux *crux, char *restart_file, int rollback_counter)
{
   crux->restore_begin(restart_file, rollback_counter);

   int int_vals[num_int_vals];

   double double_vals[num_double_vals];

   clamr_bootstrap_memory.memory_add(int_vals, size_t(num_int_vals), 4, "bootstrap_int_vals", RESTART_DATA);
   clamr_bootstrap_memory.memory_add(double_vals, size_t(num_double_vals), 8, "bootstrap_double_vals", RESTART_DATA);

   crux->restore_MallocPlus(clamr_bootstrap_memory);

   if (int_vals[ 0] != CRUX_CLAMR_VERSION) {
      printf("CRUX version mismatch for clamr data, version on file is %d, version in code is %d\n",
         int_vals[0], CRUX_CLAMR_VERSION);
      exit(0);
   }
  
   nx                        = int_vals[ 1];
   ny                        = int_vals[ 2];
   levmx                     = int_vals[ 3];
   ndim                      = int_vals[ 4];
   outputInterval            = int_vals[ 5];
   enhanced_precision_sum    = int_vals[ 6];
   niter                     = int_vals[ 7];
   it                        = int_vals[ 8];
   ncycle                    = int_vals[ 9];
   graphic_outputInterval    = int_vals[10];
   checkpoint_outputInterval = int_vals[11];
   next_cp_cycle             = int_vals[12];
   next_graphics_cycle       = int_vals[13];

   circ_radius                = double_vals[ 0];
   H_sum_initial              = double_vals[ 1];
   simTime                    = double_vals[ 2];
   deltaT                     = double_vals[ 3];
   upper_mass_diff_percentage = double_vals[ 4];

   clamr_bootstrap_memory.memory_remove(int_vals);
   clamr_bootstrap_memory.memory_remove(double_vals);

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS) {
      const char *int_vals_descriptor[num_int_vals] = {
         "CRUX_CLAMR_VERSION",
         "nx",
         "ny",
         "levmx",
         "ndim",
         "outputInterval",
         "enhanced_precision_sum",
         "niter",
         "it",
         "ncycle",
         "graphic_outputInterval",
         "checkpoint_outputInterval",
         "next_cp_cycle",
         "next_graphics_cycle"
      };
      printf("\n");
      printf("       === Restored bootstrap int_vals ===\n");
      for (int i = 0; i < num_int_vals; i++){
         printf("       %-30s %d\n",int_vals_descriptor[i], int_vals[i]);
      }
      printf("       === Restored bootstrap int_vals ===\n");
      printf("\n");
   }
#endif

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS) {
      const char *double_vals_descriptor[num_double_vals] = {
         "circ_radius",
         "H_sum_initial",
         "simTime",
         "deltaT",
         "upper_mass_diff_percentage"
      };
      printf("\n");
      printf("       === Restored bootstrap double_vals ===\n");
      for (int i = 0; i < num_double_vals; i++){
         printf("       %-30s %lg\n",double_vals_descriptor[i], double_vals[i]);
      }
      printf("       === Restored bootstrap double_vals ===\n");
      printf("\n");
   }
#endif
}

void restore_crux_data(Crux *crux)
{
   state->restore_checkpoint(crux);

   crux->restore_end();
}

