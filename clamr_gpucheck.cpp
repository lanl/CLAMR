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
#include "ezcl/ezcl.h"
#include "input.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
#include "state.h"
#include "timer/timer.h"
#include "memstats/memstats.h"
#include "PowerParser/PowerParser.hh"

using namespace PP;
#ifndef DEBUG 
#define DEBUG 0
#endif

// Sync is to reduce numerical drift between cpu and gpu
#define DO_SYNC 

static int do_comparison_calc = 1;
static int do_cpu_calc = 1;
static int do_gpu_calc = 1;

#ifdef DO_SYNC
static int do_sync = 1;
#else
static int do_sync = 0;
#endif
static int do_gpu_sync = 0;

typedef unsigned int uint;

static bool do_display_graphics = false;

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
//static struct tstart_check;
static cl_event start_write_event, end_write_event;

static double H_sum_initial = 0.0;
static double cpu_time_graphics = 0.0;
static double cpu_time_calcs    = 0.0;
static double cpu_time_partmeas = 0.0;
//static double cpu_time_check    = 0.0;

static int     ncycle  = 0;
static double  simTime = 0.0;
static double  deltaT = 0.0;

int main(int argc, char **argv) {
   int ierr;

    // Needed for code to compile correctly on the Mac
   int mype=0;
   int numpe=-1;

   parse = new PowerParser();

   //  Process command-line arguments, if any.
   parseInput(argc, argv);

   struct timeval tstart_setup;
   cpu_timer_start(&tstart_setup);
   
   numpe = 16;

   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_ACCELERATOR);
   }
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      exit(-1);
   }

   circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (real_t) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   double deltax_in = 1.0;
   double deltay_in = 1.0;

   mesh  = new Mesh(nx, ny, levmx, ndim, deltax_in, deltay_in, boundary, parallel_in, do_gpu_calc);
   if (DEBUG) {
      //if (mype == 0) mesh->print();

      char filename[10];
      sprintf(filename,"out%1d",mype);
      mesh->fp=fopen(filename,"w");

      //mesh->print_local();
   } 
   mesh->init(nx, ny, circ_radius, initial_order, do_gpu_calc);
   size_t &ncells = mesh->ncells;
   state = new State(mesh);
   state->init(do_gpu_calc);
   mesh->proc.resize(ncells);
   mesh->calc_distribution(numpe);
   state->fill_circle(circ_radius, 100.0, 7.0);
   
   if (graphic_outputInterval > niter) next_graphics_cycle = graphic_outputInterval;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;

   cl_mem &dev_H    = state->dev_H;
   cl_mem &dev_U    = state->dev_U;
   cl_mem &dev_V    = state->dev_V;

   state_t  *H        = state->H;
   state_t  *U        = state->U;
   state_t  *V        = state->V;

   state->allocate_device_memory(ncells);

   size_t one = 1;
   state->dev_deltaT   = ezcl_malloc(NULL, const_cast<char *>("dev_deltaT"), &one,    sizeof(cl_real_t),  CL_MEM_READ_WRITE, 0);

   size_t mem_request = (int)((float)ncells*mesh->mem_factor);
   dev_celltype = ezcl_malloc(NULL, const_cast<char *>("dev_celltype"), &mem_request, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_i        = ezcl_malloc(NULL, const_cast<char *>("dev_i"),        &mem_request, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_j        = ezcl_malloc(NULL, const_cast<char *>("dev_j"),        &mem_request, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_level    = ezcl_malloc(NULL, const_cast<char *>("dev_level"),    &mem_request, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);

   cl_command_queue command_queue = ezcl_get_command_queue();
   ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  &mesh->celltype[0], &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  &mesh->i[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int),  &mesh->j[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int),  &mesh->level[0],    NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_H,        CL_FALSE, 0, ncells*sizeof(cl_state_t), &H[0],        NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_U,        CL_FALSE, 0, ncells*sizeof(cl_state_t), &U[0],        NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_V,        CL_TRUE,  0, ncells*sizeof(cl_state_t), &V[0],        &end_write_event  );
   state->gpu_timers[STATE_TIMER_WRITE] += ezcl_timer_calc(&start_write_event, &end_write_event);


   if (ezcl_get_compute_device() == COMPUTE_DEVICE_ATI) enhanced_precision_sum = false;

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(enhanced_precision_sum);
   printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
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

   printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells, H_sum);

   for (int i = 0; i < MESH_COUNTER_SIZE; i++){
      mesh->cpu_counters[i]=0;
   }
   for (int i = 0; i < MESH_TIMER_SIZE; i++){
      mesh->cpu_timers[i]=0.0;
   }   

   cpu_timer_start(&tstart_cpu);
   //  Set up grid.
#ifdef GRAPHICS_OUTPUT
   mesh->write_grid(n);
#endif

#ifdef HAVE_GRAPHICS
   do_display_graphics = true;
   set_display_mysize(ncells);
   set_display_window((float)mesh->xmin, (float)mesh->xmax,
                      (float)mesh->ymin, (float)mesh->ymax);
   set_display_outline((int)outline);
   set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_display_cell_data(&state->H[0]);
   set_display_cell_proc(&mesh->proc[0]);
   set_display_viewmode(view_mode);
#endif

   if (ncycle == next_graphics_cycle){
      set_graphics_outline(outline);
      set_graphics_mysize(ncells);
      set_graphics_window((float)mesh->xmin, (float)mesh->xmax,
                          (float)mesh->ymin, (float)mesh->ymax);
      set_graphics_outline((int)outline);
      set_graphics_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
      set_graphics_cell_data(&state->H[0]);
      set_graphics_cell_proc(&mesh->proc[0]);
      set_graphics_viewmode(view_mode);

      init_graphics_output();
      set_graphics_cell_proc(&mesh->proc[0]);
      write_graphics_info(0,0,0.0,0,0);
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


extern bool neighbor_remap;

extern "C" void do_calc(void)
{  double g     = 9.80;
   double sigma = 0.95;
   int icount, jcount;

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      do_comparison_calc = 1;
      do_sync = 0;
      do_gpu_sync = 1;
   }
   
   size_t ncells    = mesh->ncells;

   cl_mem &dev_H        = state->dev_H;
   cl_mem &dev_U        = state->dev_U;
   cl_mem &dev_V        = state->dev_V;

   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;

   cl_mem &dev_mpot     = state->dev_mpot;

   vector<int>     mpot;
   
   size_t old_ncells = ncells;
   size_t new_ncells = 0;
   size_t new_ncells_gpu = 0;
   double H_sum = -1.0;

   cl_command_queue command_queue = ezcl_get_command_queue();

   //  Main loop.
   cpu_timer_start(&tstart_cpu);

   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      // To reduce drift in solution
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->H[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->U[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_state_t),  (void *)&state->V[0],  NULL);
      }

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;

      //  Calculate the real time step for the current discrete time step.
      double deltaT_cpu = state->set_timestep(g, sigma);
      
      double deltaT_gpu = state->gpu_set_timestep(sigma);

      //  Compare time step values and pass deltaT in to the kernel.
      if (do_comparison_calc)
      {  if (fabs(deltaT_gpu - deltaT_cpu) > .000001)
         {  printf("Error with deltaT calc --- cpu %lf gpu %lf\n",deltaT_cpu,deltaT_gpu); } }
      
      deltaT = (do_gpu_calc) ? deltaT_gpu : deltaT_cpu;
      simTime += deltaT;

      mesh->calc_neighbors(ncells);

      mesh->gpu_calc_neighbors();

      if (do_comparison_calc) {
         mesh->compare_neighbors_gpu_global_to_cpu_global();
      }

      cpu_timer_start(&tstart_partmeas);
      mesh->partition_measure();
      cpu_time_partmeas += cpu_timer_stop(tstart_partmeas);

      // Currently not working -- may need to be earlier?
      //if (do_cpu_calc && ! mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->calc_finite_difference(deltaT);
      
      state->gpu_calc_finite_difference(deltaT);
      
      if (do_comparison_calc) {
         // Need to compare dev_H to H, etc
         state->compare_state_gpu_global_to_cpu_global("finite difference",ncycle,ncells);
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells();
      
      mpot.resize(ncells);
      new_ncells = state->calc_refine_potential(mpot, icount, jcount);
      //printf("DEBUG cpu icount %d jcount %d new_ncells %d\n",icount,jcount,new_ncells);

      new_ncells_gpu = state->gpu_calc_refine_potential(icount, jcount);
      //printf("DEBUG gpu icount %d jcount %d new_ncells %d\n",icount,jcount,new_ncells);

      if (do_comparison_calc) {
         if (new_ncells != new_ncells_gpu) {
            printf("ERROR -- new_ncells cpu %lu not equal to new_ncells gpu %lu\n",new_ncells,new_ncells_gpu);
            exit(0);
         }
         // Need to compare dev_mpot to mpot
         if (dev_mpot) {
            mesh->compare_mpot_gpu_global_to_cpu_global(&mpot[0], dev_mpot);
         }
      }

      // Sync up cpu array with gpu version to reduce differences due to minor numerical differences
      // otherwise cell count will diverge causing code problems and crashes
      if (dev_mpot) {
         if (do_sync) {
            ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
         }
         if (do_gpu_sync) {
            ezcl_enqueue_write_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
         }
      }

      if (do_comparison_calc) {
         // This compares ioffset for each block in the calculation
         if (dev_mpot) {
            mesh->compare_ioffset_gpu_global_to_cpu_global(old_ncells, &mpot[0]);
         }
      }

      if (do_gpu_sync) {
         if (dev_mpot) {
            size_t local_work_size  = MIN(old_ncells, TILE_SIZE);
            size_t global_work_size = ((old_ncells+local_work_size - 1) /local_work_size) * local_work_size;

            //size_t block_size = (ncells + TILE_SIZE - 1) / TILE_SIZE; //  For on-device global reduction kernel.
            size_t block_size     = global_work_size/local_work_size;

            vector<int>      ioffset(block_size);
            int mtotal = 0;
            for (int ig=0; ig<(int)(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
               int mcount = 0;
               for (int ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
                   if (ic >= (int)old_ncells) break;
                   if (mesh->celltype[ic] == REAL_CELL) {
                      mcount += mpot[ic] ? 4 : 1;
                   } else {
                      mcount += mpot[ic] ? 2 : 1;
                   }
               }
               ioffset[ig] = mtotal;
               mtotal += mcount;
            }
            ezcl_enqueue_write_buffer(command_queue, mesh->dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int),       &ioffset[0], NULL);
         }
      }

      if (do_comparison_calc) {
         new_ncells = new_ncells_gpu;
      }

      //int add_ncells = new_ncells - old_ncells;
      state->rezone_all(icount, jcount, mpot);

      // Clear does not delete mpot, so have to swap with an empty vector to get
      // it to delete the mpot memory. This is all to avoid valgrind from showing
      // it as a reachable memory leak
      //mpot.clear();
      vector<int>().swap(mpot);

      //  Resize the mesh, inserting cells where refinement is necessary.
      if (state->dev_mpot) state->gpu_rezone_all(icount, jcount, localStencil);
      ncells = new_ncells;
      mesh->ncells = new_ncells;

      if (neighbor_remap && do_comparison_calc) {
         mesh->compare_neighbors_gpu_global_to_cpu_global();
      }

      //ezcl_device_memory_remove(dev_ioffset);

      if (do_comparison_calc) {
         state->compare_state_gpu_global_to_cpu_global("rezone all",ncycle,ncells);

         mesh->compare_indices_gpu_global_to_cpu_global();
      }

      //if (do_gpu_calc) {
      //   int bcount = mesh->gpu_count_BCs();
      //}

      mesh->proc.resize(ncells);
      if (icount || jcount) {  
         if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
            mesh->calc_spatial_coordinates(0);
         }
         vector<int> index(ncells);
         mesh->partition_cells(numpe, index, cycle_reorder);
         //state->state_reorder(index);
         if (do_gpu_sync) {
            ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&mesh->celltype[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&mesh->i[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&mesh->j[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_level, CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&mesh->level[0],  NULL);
         }
      }
      
      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->H[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->U[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_state_t),  (void *)&state->V[0],  NULL);
      }
   } // End burst loop

   cpu_time_calcs += cpu_timer_stop(tstart_cpu);

   H_sum = state->mass_sum(enhanced_precision_sum);
   if (isnan(H_sum)) {
      printf("Got a NAN on cycle %d\n",ncycle);
      exit(-1);
   }
   if (do_comparison_calc) {
      double total_mass = state->gpu_mass_sum(enhanced_precision_sum);
      if (fabs(total_mass - H_sum) > CONSERVATION_EPS) printf("Error: mass sum gpu %f cpu %f\n", total_mass, H_sum);/***/
   }

   printf("Iteration %3d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
      ncycle, deltaT, simTime, ncells, H_sum, H_sum - H_sum_initial);

   cpu_timer_start(&tstart_cpu);

   if(do_display_graphics || ncycle == next_graphics_cycle){
      if (do_cpu_calc){
         mesh->calc_spatial_coordinates(0);
      }
      if (do_gpu_calc){
         cl_mem dev_x  = ezcl_malloc(NULL, const_cast<char *>("dev_x"),  &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
         cl_mem dev_dx = ezcl_malloc(NULL, const_cast<char *>("dev_dx"), &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
         cl_mem dev_y  = ezcl_malloc(NULL, const_cast<char *>("dev_y"),  &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
         cl_mem dev_dy = ezcl_malloc(NULL, const_cast<char *>("dev_dy"), &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
         mesh->gpu_calc_spatial_coordinates(dev_x, dev_dx, dev_y, dev_dy);
   
         if (do_comparison_calc){
#ifdef FULL_PRECISION
            mesh->compare_coordinates_gpu_global_to_cpu_global_double(dev_x, dev_dx, dev_y, dev_dy, dev_H, &state->H[0]);
#else
            mesh->compare_coordinates_gpu_global_to_cpu_global_float(dev_x, dev_dx, dev_y, dev_dy, dev_H, &state->H[0]);
#endif
         }

         ezcl_device_memory_remove(dev_x);
         ezcl_device_memory_remove(dev_dx);
         ezcl_device_memory_remove(dev_y);
         ezcl_device_memory_remove(dev_dy);
      }
   }

   if(ncycle == next_graphics_cycle){
      set_graphics_mysize(ncells);
      set_graphics_viewmode(view_mode);
      set_graphics_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
      set_graphics_cell_data(&state->H[0]);
      set_graphics_cell_proc(&mesh->proc[0]);

      write_graphics_info(ncycle/graphic_outputInterval,ncycle,simTime,0,0);
      next_graphics_cycle += graphic_outputInterval;
   }

#ifdef HAVE_GRAPHICS
   set_display_mysize(ncells);
   set_display_viewmode(view_mode);
   set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_display_cell_data(&state->H[0]);
   set_display_cell_proc(&mesh->proc[0]);
   set_display_circle_radius(circle_radius);
   draw_scene();

#endif

   cpu_time_graphics += cpu_timer_stop(tstart_cpu);

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      if(graphic_outputInterval < niter){
         cpu_timer_start(&tstart_cpu);

         mesh->calc_spatial_coordinates(0);
#ifdef HAVE_GRAPHICS
         set_display_mysize(ncells);
         set_display_viewmode(view_mode);
         set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
         set_display_cell_data(&state->H[0]);
         set_display_cell_proc(&mesh->proc[0]);
#endif

         write_graphics_info(ncycle/graphic_outputInterval,ncycle,simTime,0,0);
         next_graphics_cycle += graphic_outputInterval;

         cpu_time_graphics += cpu_timer_stop(tstart_cpu);
      }

      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      long long mem_used = memstats_memused();
      if (mem_used > 0) {
         printf("Memory used      %lld kB\n",mem_used);
         printf("Memory peak      %lld kB\n",memstats_mempeak());
         printf("Memory free      %lld kB\n",memstats_memfree());
         printf("Memory available %lld kB\n",memstats_memtotal());
      }
      state->output_timing_info(do_cpu_calc, do_gpu_calc, elapsed_time);
      mesh->parallel_output("CPU:  calc incl part meas     time was",cpu_time_calcs,    0, "s");
      mesh->parallel_output("CPU:  calculation only        time was",cpu_time_calcs-cpu_time_partmeas,    0, "s");
      mesh->parallel_output("CPU:  partition measure       time was",cpu_time_partmeas, 0, "s");
      mesh->parallel_output("CPU:  graphics                time was",cpu_time_graphics, 0, "s");
      //mesh->parallel_output("CPU:  check                   time was",cpu_time_check,    0, "s");

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      printf("CPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_REZONE)/(double)ncycle*100.0 );
      printf("CPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_cpu_counter(MESH_COUNTER_CALC_NEIGH)/(double)ncycle*100.0 );
      if (mesh->get_cpu_counter(MESH_COUNTER_REZONE) > 0) {
         printf("CPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_cpu_counter(MESH_COUNTER_REFINE_SMOOTH)/(double)mesh->get_cpu_counter(MESH_COUNTER_REZONE) );
      }
      printf("GPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_gpu_counter(MESH_COUNTER_REZONE)/(double)ncycle*100.0 );
      printf("GPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_gpu_counter(MESH_COUNTER_CALC_NEIGH)/(double)ncycle*100.0 );
      if (mesh->get_gpu_counter(MESH_COUNTER_REZONE) > 0) {
         printf("GPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_gpu_counter(MESH_COUNTER_REFINE_SMOOTH)/(double)mesh->get_gpu_counter(MESH_COUNTER_REZONE) );
      }


      mesh->terminate();
      state->terminate();
      ezcl_terminate();
      terminate_graphics_output();

      delete mesh;
      delete state;

      ezcl_mem_walk_all();

      exit(0);
   }  //  Complete final output.
}

