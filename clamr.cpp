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
#include "ezcl/ezcl.h"
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

static int do_cpu_calc = 0;
static int do_gpu_calc = 1;

#ifdef HAVE_CL_DOUBLE
typedef double      real;
typedef cl_double   cl_real;
#define CONSERVATION_EPS    .00001
#define STATE_EPS        .025
#define MPI_C_REAL MPI_DOUBLE
#define L7_REAL L7_DOUBLE
#else
typedef float       real;
typedef cl_float    cl_real;
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
            enhanced_precision_sum;//  Flag for enhanced precision sum (default true); init in input.cpp::parseInput().
int         outputInterval, //  Periodicity of output; init in input.cpp::parseInput().
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
static cl_event start_write_event, end_write_event;

static cl_context          context                 = NULL;
static cl_command_queue    command_queue           = NULL;
static int compute_device = 0;

static double H_sum_initial = 0.0;
static long gpu_time_graphics = 0;
double cpu_time_main_setup = 0.0;
static double cpu_time_timestep = 0.0;
static double cpu_time_finite_diff = 0.0;
static double cpu_time_refine_potential = 0.0;
static double cpu_time_rezone = 0.0;
static double cpu_time_neighbors = 0.0;
static double cpu_time_load_balance = 0.0;

int main(int argc, char **argv) {
   int ierr;

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   L7_Init(&mype, &numpe, &argc, argv);
   parseInput(argc, argv);

   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU, &context, &command_queue, &compute_device, mype);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU, &context, &command_queue, &compute_device, mype);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      L7_Terminate();
      exit(-1);
   }

   struct timeval tstart_setup;
   cpu_timer_start(&tstart_setup);

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 1;

   mesh = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   if (DEBUG) {
      //if (mype == 0) mesh->print();

      char filename[10];
      sprintf(filename,"out%1d",mype);
      mesh->fp=fopen(filename,"w");

      //mesh->print_local();
   }

   mesh->init(nx, ny, circ_radius, context, initial_order, compute_device, do_gpu_calc);

   size_t &ncells = mesh->ncells;
   size_t &ncells_global = mesh->ncells_global;
   int &noffset = mesh->noffset;

   MPI_Allreduce(&ncells, &ncells_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

   state = new State(ncells);
   state->init(ncells, context, compute_device, do_gpu_calc);

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

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;

   vector<real> &H = state->H;
   vector<real> &U = state->U;
   vector<real> &V = state->V;

   vector<real> &x  = mesh->x;
   vector<real> &dx = mesh->dx;
   vector<real> &y  = mesh->y;
   vector<real> &dy = mesh->dy;

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

   celltype.resize(ncells);
   level.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);

   dev_celltype  = ezcl_malloc(NULL, const_cast<char *>("dev_celltype"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_i         = ezcl_malloc(NULL, const_cast<char *>("dev_i"),        &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY,  0);
   dev_j         = ezcl_malloc(NULL, const_cast<char *>("dev_j"),        &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY,  0);
   dev_level     = ezcl_malloc(NULL, const_cast<char *>("dev_level"),    &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   H.resize(ncells);
   U.resize(ncells);
   V.resize(ncells);

   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);

   mesh->calc_spatial_coordinates(0);

   state->fill_circle(mesh, circ_radius, 100.0, 5.0);

   state->allocate_device_memory(ncells);

   size_t one = 1;
   state->dev_deltaT   = ezcl_malloc(NULL, const_cast<char *>("dev_deltaT"),               &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

   //  Set write buffers for data.
   ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_V, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);

   ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],        NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],        NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_level,    CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&level[0],    &end_write_event);
   state->gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

   mesh->nlft.clear();
   mesh->nrht.clear();
   mesh->nbot.clear();
   mesh->ntop.clear();

   mesh->dev_nlft = NULL;
   mesh->dev_nrht = NULL;
   mesh->dev_nbot = NULL;
   mesh->dev_ntop = NULL;

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(mesh, enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   double cpu_time_main_setup = cpu_timer_stop(tstart_setup);
   state->parallel_timer_output(numpe,mype,"CPU:  setup time               time was",cpu_time_main_setup);

   long long mem_used = timer_memused();
   if (mem_used > 0) {
      state->parallel_memory_output(numpe,mype,"Memory used      in startup ",mem_used);
      state->parallel_memory_output(numpe,mype,"Memory peak      in startup ",timer_mempeak());
      state->parallel_memory_output(numpe,mype,"Memory free      at startup ",timer_memfree());
      state->parallel_memory_output(numpe,mype,"Memory available at startup ",timer_memtotal());
   }

   if (mype == 0) {
      printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells_global, H_sum);
   }

   mesh->cpu_calc_neigh_counter=0;
   mesh->cpu_time_calc_neighbors=0.0;
   mesh->cpu_rezone_counter=0;
   mesh->cpu_time_rezone_all=0.0;
   mesh->cpu_refine_smooth_counter=0;

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_mysize(ncells_global);
   vector<real> H_global;
   vector<real> x_global;
   vector<real> dx_global;
   vector<real> y_global;
   vector<real> dy_global;
   vector<int> proc_global;
   if (mype == 0){
      H_global.resize(ncells_global);
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      proc_global.resize(ncells_global);
   }
   MPI_Gatherv(&x[0],  nsizes[mype], MPI_C_REAL, &x_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dx[0], nsizes[mype], MPI_C_REAL, &dx_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&y[0],  nsizes[mype], MPI_C_REAL, &y_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dy[0], nsizes[mype], MPI_C_REAL, &dy_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&H[0], nsizes[mype], MPI_C_REAL, &H_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);

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
   set_cell_data(&H[0]);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_proc(&mesh->proc[0]);
#endif

   set_window(mesh->xmin, mesh->xmax, mesh->ymin, mesh->ymax);
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
   struct timeval tstart_cpu;

   double sigma = 0.95; 
   int icount=0;
   static cl_event start_read_event, end_read_event;

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      printf("Can't do this problem with GPU\n");
      exit(0);
   }
   
   //  Initialize state variables for GPU calculation.
   int &mype  = mesh->mype;
   int &numpe = mesh->numpe;

   vector<int>   &nsizes   = mesh->nsizes;
   vector<int>   &ndispl   = mesh->ndispl;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh->ncells_global;
   size_t &ncells           = mesh->ncells;
   size_t &ncells_ghost     = mesh->ncells_ghost;

   vector<real>  &H = state->H;

   vector<real>  &x  = mesh->x;
   vector<real>  &dx = mesh->dx;
   vector<real>  &y  = mesh->y;
   vector<real>  &dy = mesh->dy;

   cl_mem &dev_H  = state->dev_H;
   cl_mem &dev_U  = state->dev_U;
   cl_mem &dev_V  = state->dev_V;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   size_t old_ncells = ncells;
   size_t old_ncells_global = ncells_global;
   size_t new_ncells = 0;
   double H_sum = -1.0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      size_t local_work_size  = MIN(ncells, TILE_SIZE);
      size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
      size_t block_size     = global_work_size/local_work_size;

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;
      old_ncells_global = ncells_global;

      cpu_timer_start(&tstart_cpu);
      //  Calculate the real time step for the current discrete time step.
      deltaT = state->gpu_set_timestep(command_queue, mesh, sigma);
      simTime += deltaT;
      cpu_time_timestep += cpu_timer_stop(tstart_cpu);

      cpu_timer_start(&tstart_cpu);
      if (mesh->dev_nlft == NULL) mesh->gpu_calc_neighbors_local(command_queue);
      cpu_time_neighbors += cpu_timer_stop(tstart_cpu);

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      cpu_timer_start(&tstart_cpu);
      state->gpu_calc_finite_difference(command_queue, mesh, deltaT);
      cpu_time_finite_diff += cpu_timer_stop(tstart_cpu);

      vector<int>      ioffset(block_size);

      cpu_timer_start(&tstart_cpu);
      new_ncells = state->gpu_calc_refine_potential(command_queue, mesh);
      cpu_time_refine_potential += cpu_timer_stop(tstart_cpu);

      int ncells_global_old = ncells_global;
      MPI_Allreduce(&new_ncells, &ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      //printf("%d: DEBUG ncells is %d new_ncells %d old_ncells %d ncells_global %d\n",mype, ncells, new_ncells, old_ncells, ncells_global);

      //  Resize the mesh, inserting cells where refinement is necessary.
      cpu_timer_start(&tstart_cpu);
      if (ncells_global_old != (int)ncells_global) state->gpu_rezone_all(command_queue, mesh, old_ncells, new_ncells, old_ncells, localStencil);
      cpu_time_rezone += cpu_timer_stop(tstart_cpu);

      // XXX XXX XXX
      ncells       = new_ncells;
      //mesh->ncells = new_ncells;

      cpu_timer_start(&tstart_cpu);
      if (mesh->dev_nlft == NULL) mesh->gpu_do_load_balance_local(command_queue, new_ncells, ncells_global, dev_H, dev_U, dev_V);
      cpu_time_load_balance += cpu_timer_stop(tstart_cpu);

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
      H_sum = state->gpu_mass_sum(command_queue, mesh, enhanced_precision_sum);
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
   cl_mem dev_x  = ezcl_malloc(NULL, const_cast<char *>("dev_x"),  &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_dx = ezcl_malloc(NULL, const_cast<char *>("dev_dx"), &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_y  = ezcl_malloc(NULL, const_cast<char *>("dev_y"),  &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_dy = ezcl_malloc(NULL, const_cast<char *>("dev_dy"), &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

   mesh->gpu_calc_spatial_coordinates(command_queue, dev_x, dev_dx, dev_y, dev_dy);

   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);
   H.resize(max(ncells,ncells_ghost));

   ezcl_enqueue_read_buffer(command_queue, dev_x,  CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&x[0],  &start_read_event);
   ezcl_enqueue_read_buffer(command_queue, dev_dx, CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&dx[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,  CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&y[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy, CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&dy[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,  CL_TRUE,  0, ncells*sizeof(cl_real), (void *)&H[0],  &end_read_event);

   gpu_time_graphics += ezcl_timer_calc(&start_read_event, &end_read_event);

   cpu_timer_start(&tstart_cpu);

   ezcl_device_memory_remove(dev_x);
   ezcl_device_memory_remove(dev_dx);
   ezcl_device_memory_remove(dev_y);
   ezcl_device_memory_remove(dev_dy);

#ifdef HAVE_OPENGL
   set_mysize(ncells_global);
   vector<real> x_global;
   vector<real> dx_global;
   vector<real> y_global;
   vector<real> dy_global;
   vector<real> H_global;
   vector<int> proc_global;
   if (mype == 0) {
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      H_global.resize(ncells_global);
      proc_global.resize(ncells_global);
   }

   MPI_Gatherv(&x[0],  ncells, MPI_C_REAL, &x_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dx[0], ncells, MPI_C_REAL, &dx_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&y[0],  ncells, MPI_C_REAL, &y_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&dy[0], ncells, MPI_C_REAL, &dy_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&H[0],  ncells, MPI_C_REAL, &H_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);

   if (view_mode == 0) {
      mesh->proc.resize(ncells);
      for (size_t ii = 0; ii<ncells; ii++){
         mesh->proc[ii] = mesh->mype;
      }
   
      MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   }

   set_viewmode(view_mode);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_data(&H_global[0]);
   set_cell_proc(&proc_global[0]);
   set_circle_radius(circle_radius);
   draw_scene();
   MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef HAVE_MPE
   set_mysize(ncells);
   set_viewmode(view_mode);
   set_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
   set_cell_data(&H[0]);
   set_cell_proc(&mesh->proc[0]);
   set_circle_radius(circle_radius);
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
      
      long long mem_used = timer_memused();
      if (mem_used > 0) {
         state->parallel_memory_output(numpe,mype,"Memory used      ",mem_used);
         state->parallel_memory_output(numpe,mype,"Memory peak      ",timer_mempeak());
         state->parallel_memory_output(numpe,mype,"Memory free      ",timer_memfree());
         state->parallel_memory_output(numpe,mype,"Memory available ",timer_memtotal());
      }

      state->output_timing_info(mesh, do_cpu_calc, do_gpu_calc, elapsed_time);
      state->parallel_timer_output(numpe,mype,"CPU:  setup time               time was",cpu_time_main_setup);
      state->parallel_timer_output(numpe,mype,"GPU:  graphics                 time was",(double) gpu_time_graphics * 1.0e-9 );

      state->parallel_timer_output(numpe,mype,"CPU:  timestep calc            time was",cpu_time_timestep);
      state->parallel_timer_output(numpe,mype,"CPU:  finite_diff              time was",cpu_time_finite_diff);
      state->parallel_timer_output(numpe,mype,"CPU:  refine_potential         time was",cpu_time_refine_potential);
      state->parallel_timer_output(numpe,mype,"CPU:  rezone                   time was",cpu_time_rezone);
      state->parallel_timer_output(numpe,mype,"CPU:  neighbors                time was",cpu_time_neighbors);
      state->parallel_timer_output(numpe,mype,"CPU:  load_balance             time was",cpu_time_load_balance);

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      if (mype ==0) {
         printf("GPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_gpu_rezone_count()/(double)ncycle*100.0 );
         printf("GPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_gpu_calc_neigh_count()/(double)ncycle*100.0 );
         printf("GPU:  load balance frequency          \t %8.4f\tpercent\n",     (double)mesh->get_gpu_load_balance_count()/(double)ncycle*100.0 );
         printf("GPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_gpu_refine_smooth_count()/(double)mesh->get_gpu_rezone_count() );
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

      ezcl_mem_walk_all();

      //  Release kernels and finalize the OpenCL elements.
      ezcl_finalize();
      
      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
}

