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

#include "ezcl/ezcl.h"
#include "input.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
#include "state.h"
#include "l7/l7.h"
#include "timer/timer.h"
#include "memstats/memstats.h"

#include <mpi.h>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>
#include "graphics/display.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

static int do_cpu_calc = 0;
static int do_gpu_calc = 1;

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
static Mesh       *mesh;     //  Object containing mesh information; init in grid.cpp::main().
static State      *state;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
static struct timeval tstart;
static cl_event start_write_event, end_write_event;

static double H_sum_initial = 0.0;
static long gpu_time_graphics = 0;
double cpu_time_main_setup = 0.0;

int main(int argc, char **argv) {
   int ierr;

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   parseInput(argc, argv);
   L7_Init(&mype, &numpe, &argc, argv, do_quo_setup, lttrace_on);

   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      L7_Terminate();
      exit(-1);
   }
   L7_Dev_Init();

   struct timeval tstart_setup;
   cpu_timer_start(&tstart_setup);

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
   size_t &ncells_global = mesh->ncells_global;
   int &noffset = mesh->noffset;

   state = new State(mesh);
   state->init(do_gpu_calc);

   cl_mem &dev_corners_i_local  = mesh->dev_corners_i;
   cl_mem &dev_corners_j_local  = mesh->dev_corners_j;

   vector<int>   &corners_i_local  = mesh->corners_i;
   vector<int>   &corners_j_local  = mesh->corners_j;

   vector<int>   &nsizes     = mesh->nsizes;
   vector<int>   &ndispl     = mesh->ndispl;

   cl_mem &dev_H  = state->dev_H;
   cl_mem &dev_U  = state->dev_U;
   cl_mem &dev_V  = state->dev_V;

   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;

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

   size_t corners_size = corners_i_local.size();

   dev_corners_i_local  = ezcl_malloc(&corners_i_local[0],  const_cast<char *>("dev_corners_i_local"), &corners_size, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
   dev_corners_j_local  = ezcl_malloc(&corners_j_local[0],  const_cast<char *>("dev_corners_j_local"), &corners_size, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);

   size_t mem_request = (int)((float)ncells*mesh->mem_factor);
   dev_celltype  = ezcl_malloc(NULL, const_cast<char *>("dev_celltype"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_i         = ezcl_malloc(NULL, const_cast<char *>("dev_i"),        &mem_request, sizeof(cl_int),  CL_MEM_READ_ONLY,  0);
   dev_j         = ezcl_malloc(NULL, const_cast<char *>("dev_j"),        &mem_request, sizeof(cl_int),  CL_MEM_READ_ONLY,  0);
   dev_level     = ezcl_malloc(NULL, const_cast<char *>("dev_level"),    &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   state->resize(ncells);

   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);

   mesh->calc_spatial_coordinates(0);

   state->fill_circle(circ_radius, 100.0, 7.0);

   state->allocate_device_memory(ncells);

   size_t one = 1;
   state->dev_deltaT   = ezcl_malloc(NULL, const_cast<char *>("dev_deltaT"),               &one,    sizeof(cl_real_t),  CL_MEM_READ_WRITE, 0);

   //  Set write buffers for data.
   cl_command_queue command_queue = ezcl_get_command_queue();
   ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->H[0],  &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->U[0],  NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_V, CL_FALSE, 0, ncells*sizeof(cl_state_t),  (void *)&state->V[0],  NULL);

   ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&mesh->celltype[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&mesh->i[0],        NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&mesh->j[0],        NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_level,    CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&mesh->level[0],    &end_write_event);
   state->gpu_timers[STATE_TIMER_WRITE] += ezcl_timer_calc(&start_write_event, &end_write_event);

   mesh->dev_nlft = NULL;
   mesh->dev_nrht = NULL;
   mesh->dev_nbot = NULL;
   mesh->dev_ntop = NULL;
   state->dev_mpot = NULL;

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

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_display_mysize(ncells_global);
   vector<state_t> H_global;
   vector<spatial_t> x_global;
   vector<spatial_t> dx_global;
   vector<spatial_t> y_global;
   vector<spatial_t> dy_global;
   vector<int> proc_global;
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

   set_display_cell_data(&H_global[0]);
   set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);

   if (view_mode == 0) {
      mesh->proc.resize(ncells);
      for (size_t ii = 0; ii<ncells; ii++){
         mesh->proc[ii] = mesh->mype;
      }
   
      MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_INT, 0, MPI_COMM_WORLD);
   }

   set_display_cell_proc(&proc_global[0]);
#endif
#ifdef HAVE_MPE
   set_display_mysize(ncells);
   set_display_cell_data(&state->H[0]);
   set_display_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_display_cell_proc(&mesh->proc[0]);
#endif

   set_display_window((float)mesh->xmin, (float)mesh->xmax, (float)mesh->ymin, (float)mesh->ymax);
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
{
   double sigma = 0.95; 
   //int icount=0;

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      printf("Can't do this problem with GPU\n");
      exit(0);
   }
   
   //  Initialize state variables for GPU calculation.
   int &mype  = mesh->mype;
   int &numpe = mesh->numpe;

#ifdef HAVE_GRAPHICS
   struct timeval tstart_cpu;

#ifdef HAVE_OPENGL
   vector<int>   &nsizes   = mesh->nsizes;
   vector<int>   &ndispl   = mesh->ndispl;
#endif

   vector<spatial_t>  &x  = mesh->x;
   vector<spatial_t>  &dx = mesh->dx;
   vector<spatial_t>  &y  = mesh->y;
   vector<spatial_t>  &dy = mesh->dy;

   cl_mem &dev_H  = state->dev_H;

#endif

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh->ncells_global;
   size_t &ncells           = mesh->ncells;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   size_t new_ncells = 0;
   double H_sum = -1.0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      size_t local_work_size  = MIN(ncells, TILE_SIZE);
      size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
      size_t block_size     = global_work_size/local_work_size;

      //  Calculate the real time step for the current discrete time step.
      deltaT = state->gpu_set_timestep(sigma);
      simTime += deltaT;

      mesh->gpu_calc_neighbors_local();

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->gpu_calc_finite_difference(deltaT);

      vector<int>      ioffset(block_size);

      int icount = 0;
      int jcount = 0;
      new_ncells = state->gpu_calc_refine_potential(icount, jcount);

      //int ncells_global_old = ncells_global;
      MPI_Allreduce(&new_ncells, &ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      //  Resize the mesh, inserting cells where refinement is necessary.
      if (state->dev_mpot) state->gpu_rezone_all(icount, jcount, localStencil);

      ncells       = new_ncells;

      state->gpu_do_load_balance_local(new_ncells);

      ioffset.clear();

      H_sum = -1.0;

//      mesh->proc.resize(ncells);
//      if (icount) {
//         vector<int> index(ncells);
//         vector<int> index_global(ncells_global);
//         mesh->partition_cells(numpe, index, cycle_reorder);
//         //state->state_reorder(index);
//      }

   }  //  End burst loop

   if (H_sum < 0) {
      H_sum = state->gpu_mass_sum(enhanced_precision_sum);
   }
   if (isnan(H_sum)) {
      printf("Got a NAN on cycle %d\n",ncycle);
      exit(-1);
   }
   if (mype == 0){
      printf("Iteration %3d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

#ifdef HAVE_GRAPHICS
   static cl_event start_read_event, end_read_event;

   cl_mem dev_x  = ezcl_malloc(NULL, const_cast<char *>("dev_x"),  &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_dx = ezcl_malloc(NULL, const_cast<char *>("dev_dx"), &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_y  = ezcl_malloc(NULL, const_cast<char *>("dev_y"),  &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_dy = ezcl_malloc(NULL, const_cast<char *>("dev_dy"), &ncells, sizeof(cl_spatial_t),  CL_MEM_READ_WRITE, 0);

   mesh->gpu_calc_spatial_coordinates(dev_x, dev_dx, dev_y, dev_dy);

   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);
   //H.resize(max(ncells,ncells_ghost));
   vector<state_t>H_graphics(ncells);

   cl_command_queue command_queue = ezcl_get_command_queue();
   ezcl_enqueue_read_buffer(command_queue, dev_x,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), (void *)&x[0],  &start_read_event);
   ezcl_enqueue_read_buffer(command_queue, dev_dx, CL_FALSE, 0, ncells*sizeof(cl_spatial_t), (void *)&dx[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), (void *)&y[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy, CL_FALSE, 0, ncells*sizeof(cl_spatial_t), (void *)&dy[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,  CL_TRUE,  0, ncells*sizeof(cl_state_t), (void *)&H_graphics[0],  &end_read_event);

   gpu_time_graphics += ezcl_timer_calc(&start_read_event, &end_read_event);

   cpu_timer_start(&tstart_cpu);

   ezcl_device_memory_remove(dev_x);
   ezcl_device_memory_remove(dev_dx);
   ezcl_device_memory_remove(dev_y);
   ezcl_device_memory_remove(dev_dy);

#ifdef HAVE_OPENGL
   set_display_mysize(ncells_global);
   vector<spatial_t> x_global;
   vector<spatial_t> dx_global;
   vector<spatial_t> y_global;
   vector<spatial_t> dy_global;
   vector<state_t> H_graphics_global;
   vector<int> proc_global;
   if (mype == 0) {
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      H_graphics_global.resize(ncells_global);
      proc_global.resize(ncells_global);
   }

   MPI_Gatherv(&x[0],  ncells, MPI_SPATIAL_T, &x_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dx[0], ncells, MPI_SPATIAL_T, &dx_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&y[0],  ncells, MPI_SPATIAL_T, &y_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dy[0], ncells, MPI_SPATIAL_T, &dy_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&H_graphics[0],  ncells, MPI_STATE_T, &H_graphics_global[0],  &nsizes[0], &ndispl[0], MPI_STATE_T, 0, MPI_COMM_WORLD);

   if (view_mode == 0) {
      mesh->proc.resize(ncells);
      for (size_t ii = 0; ii<ncells; ii++){
         mesh->proc[ii] = mesh->mype;
      }
   
      MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_INT, 0, MPI_COMM_WORLD);
   }

   set_display_viewmode(view_mode);
   set_display_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_display_cell_data(&H_graphics_global[0]);
   set_display_cell_proc(&proc_global[0]);
   set_display_circle_radius(circle_radius);
   draw_scene();
   MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef HAVE_MPE
   set_display_mysize(ncells);
   set_display_viewmode(view_mode);
   set_display_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
   set_display_cell_data(&H_graphics[0]);
   set_display_cell_proc(&mesh->proc[0]);
   set_display_circle_radius(circle_radius);
   draw_scene();
   MPI_Barrier(MPI_COMM_WORLD);
   if (ncycle == 6) sleep(300);
#endif

   gpu_time_graphics += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
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
      mesh->parallel_output("CPU:  setup time               time was",cpu_time_main_setup, 0, "s");
      mesh->parallel_output("GPU:  graphics                 time was",(double) gpu_time_graphics * 1.0e-9, 0, "s");

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      if (mype ==0) {
         printf("GPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_gpu_counter(MESH_COUNTER_REZONE)/(double)ncycle*100.0 );
         printf("GPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_gpu_counter(MESH_COUNTER_CALC_NEIGH)/(double)ncycle*100.0 );
         printf("GPU:  load balance frequency          \t %8.4f\tpercent\n",     (double)mesh->get_gpu_counter(MESH_COUNTER_LOAD_BALANCE)/(double)ncycle*100.0 );
         printf("GPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_gpu_counter(MESH_COUNTER_REFINE_SMOOTH)/(double)mesh->get_gpu_counter(MESH_COUNTER_REZONE) );
      }

      ezcl_device_memory_remove(mesh->dev_corners_i);
      ezcl_device_memory_remove(mesh->dev_corners_j);

      if (mesh->dev_nlft != NULL){
         ezcl_device_memory_remove(mesh->dev_nlft);
         ezcl_device_memory_remove(mesh->dev_nrht);
         ezcl_device_memory_remove(mesh->dev_nbot);
         ezcl_device_memory_remove(mesh->dev_ntop);
      }

      mesh->terminate();
      state->terminate();
      ezcl_terminate();
      if (numpe > 1) L7_Free(&mesh->cell_handle);
      L7_Dev_Free();

      delete mesh;
      delete state;

      //  Release kernels and finalize the OpenCL elements.
      //ezcl_finalize();
      
      ezcl_mem_walk_all();

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
}

