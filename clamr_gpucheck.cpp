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
#include "display.h"
#include "ezcl/ezcl.h"
#include "input.h"
#include "mesh.h"
#include "partition.h"
#include "state.h"
#include "timer/timer.h"

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

#ifdef HAVE_CL_DOUBLE
typedef double      real;
typedef struct
{
   double s0;
   double s1;
}  real2;
typedef cl_double   cl_real;
typedef cl_double2  cl_real2;
typedef cl_double4  cl_real4;
typedef cl_double8  cl_real8;
#define CONSERVATION_EPS    .02
#define STATE_EPS        .02
#else
typedef float       real;
typedef struct
{
   float s0;
   float s1;
}  real2;
typedef cl_float    cl_real;
typedef cl_float2   cl_real2;
typedef cl_float4   cl_real4;
typedef cl_float8   cl_real8;
#define CONSERVATION_EPS    .2
#define STATE_EPS      15.0
#endif

typedef unsigned int uint;

static double circle_radius=-1.0;

static int view_mode = 0;

bool        verbose,        //  Flag for verbose command-line output; init in input.cpp::parseInput().
            localStencil,   //  Flag for use of local stencil; init in input.cpp::parseInput().
            outline,        //  Flag for drawing outlines of cells; init in input.cpp::parseInput().
            enhanced_precision_sum,//  Flag for enhanced precision sum (default true); init in input.cpp::parseInput().
            special_case;   //  Flag for special case debugging (default false); init in input.cpp::parseInput().
int         outputInterval, //  Periodicity of output; init in input.cpp::parseInput().
            levmx,          //  Maximum number of refinement levels; init in input.cpp::parseInput().
            nx,             //  x-resolution of coarse grid; init in input.cpp::parseInput().
            ny,             //  y-resolution of coarse grid; init in input.cpp::parseInput().
            niter,          //  Maximum time step; init in input.cpp::parseInput().
            numpe,          //  
            ndim    = 2;    //  Dimensionality of problem (2 or 3).

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
static Mesh       *mesh;           //  Object containing mesh information; init in grid.cpp::main().
static State      *state;          //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
static struct timeval tstart;
static cl_event start_write_event, end_write_event;

static cl_context          context                 = NULL;
static cl_command_queue    command_queue           = NULL;
static int compute_device = 0;

static double  H_sum_initial = 0.0;

int main(int argc, char **argv) {
   int ierr;

    // Needed for code to compile correctly on the Mac
   int mype=0;
   int numpe=-1;

   //  Process command-line arguments, if any.
   parseInput(argc, argv);
   
   numpe = 16;

   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU, &context, &command_queue, &compute_device, 0);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU, &context, &command_queue, &compute_device, 0);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      exit(-1);
   }

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   
   mesh  = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   mesh->init(nx, ny, circ_radius, context, initial_order, special_case, compute_device, do_gpu_calc);
   size_t &ncells = mesh->ncells;
   state = new State(ncells, context);
   state->init(ncells, context, compute_device, do_gpu_calc);
   mesh->proc.resize(ncells);
   mesh->calc_distribution(numpe, mesh->proc);
   state->fill_circle(mesh, circ_radius, 100.0, 5.0);
   
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;

   cl_mem &dev_H    = state->dev_H;
   cl_mem &dev_U    = state->dev_U;
   cl_mem &dev_V    = state->dev_V;

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;

   vector<real>  &H        = state->H;
   vector<real>  &U        = state->U;
   vector<real>  &V        = state->V;

   state->allocate_device_memory(ncells);

   size_t one = 1;
   state->dev_deltaT   = ezcl_malloc(NULL, const_cast<char *>("dev_deltaT"), &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

   dev_celltype = ezcl_malloc(NULL, const_cast<char *>("dev_celltype"), &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_i        = ezcl_malloc(NULL, const_cast<char *>("dev_i"),        &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_j        = ezcl_malloc(NULL, const_cast<char *>("dev_j"),        &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_level    = ezcl_malloc(NULL, const_cast<char *>("dev_level"),    &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);

   ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0], &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&level[0],    NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_H,        CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],       NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_U,        CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],       NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_V,        CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],       &end_write_event  );
   state->gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

   if (compute_device == COMPUTE_DEVICE_ATI) enhanced_precision_sum = false;

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(mesh, enhanced_precision_sum);
   printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells, H_sum);

   //  Set up grid.
#ifdef GRAPHICS_OUTPUT
   mesh->write_grid(n);
#endif

#ifdef HAVE_GRAPHICS
   set_mysize(ncells);
   set_viewmode(view_mode);
   set_window(mesh->xmin, mesh->xmax, mesh->ymin, mesh->ymax);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water", mype);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_data(&H[0]);
   set_cell_proc(&mesh->proc[0]);
   set_circle_radius(circle_radius);
   draw_scene();
   //if (verbose) sleep(5);
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

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      do_comparison_calc = 1;
      do_sync = 0;
      do_gpu_sync = 1;
   }
   
   //  Initialize state variables for GPU calculation.
   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;

   size_t ncells    = mesh->ncells;

   vector<real>  &H        = state->H;
   vector<real>  &U        = state->U;
   vector<real>  &V        = state->V;

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
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      // To reduce drift in solution
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);
      }

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;

      //  Calculate the real time step for the current discrete time step.
      double deltaT_cpu=0.0;
      if (do_cpu_calc) {
         deltaT_cpu = state->set_timestep(mesh, g, sigma);
      }  //  Complete CPU timestep calculation.
      
      double deltaT_gpu=0.0;
      if (do_gpu_calc) {
         deltaT_gpu = state->gpu_set_timestep(command_queue, mesh, sigma);
      }  //  Complete GPU calculation.

      //  Compare time step values and pass deltaT in to the kernel.
      if (do_comparison_calc)
      {  if (fabs(deltaT_gpu - deltaT_cpu) > .000001)
         {  printf("Error with deltaT calc --- cpu %lf gpu %lf\n",deltaT_cpu,deltaT_gpu); } }
      
      deltaT = (do_gpu_calc) ? deltaT_gpu : deltaT_cpu;
      simTime += deltaT;

      if (do_cpu_calc) {
         mesh->calc_neighbors();
      }

      if (do_gpu_calc) {
         mesh->gpu_calc_neighbors(command_queue);
      }

      if (do_comparison_calc) {
         mesh->compare_neighbors_gpu_global_to_cpu_global(command_queue);
      }

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (do_cpu_calc && ! mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      if (do_cpu_calc) {
         state->calc_finite_difference(mesh, deltaT);
      }
      
      if (do_gpu_calc) {
         state->gpu_calc_finite_difference(command_queue, mesh, deltaT);
      }
      
      if (do_comparison_calc) {
         // Need to compare dev_H to H, etc
         state->compare_state_gpu_global_to_cpu_global(command_queue,"finite difference",ncycle,ncells);
      }

      if (compute_device == COMPUTE_DEVICE_ATI) {
         printf("DEBUG file %s line %d\n",__FILE__,__LINE__);
         fflush(stdout);
         exit(0);
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells(mesh);
      
      //  Set flag for checking symmetry.
      int flag = 0;
      
#ifdef CHECK_SYMMETRY
      double eps=0.0001;
      
      vector<int> dsym(ncells);
      vector<int> xsym(ncells);
      vector<int> ysym(ncells);
      
      mesh->calc_symmetry(dsym, xsym, ysym);
      state->symmetry_check(mesh, "Diagonal", dsym, eps, DIAG_RULE, flag);
      state->symmetry_check(mesh, "Xaxis",    xsym, eps, X_RULE,    flag);
      state->symmetry_check(mesh, "Yaxis",    ysym, eps, Y_RULE,    flag);
#endif
      
      if (flag) {
         printf("Cycle is %d\n",ncycle);
         sleep(2);
      }

      //  Check for NANs.
      for (uint ic=0; ic<ncells; ic++) {
         if (isnan(H[ic]))
         {  printf("Got a NAN on cell %d cycle %d\n",ic,ncycle);
            H[ic]=0.0;
            //sleep(100);
            //  Release kernels and finalize the OpenCL elements.
            ezcl_finalize();
            exit(-1); }
      }  //  Complete NAN check.
      
      if (do_cpu_calc) {
         mpot.resize(ncells);
         state->calc_refine_potential(mesh, mpot, icount, jcount);
      }  //  Complete CPU calculation.

      if (do_gpu_calc) {
         new_ncells_gpu = state->gpu_calc_refine_potential(command_queue, mesh);
      }
      
      if (do_comparison_calc) {
         // Need to compare dev_mpot to mpot
         mesh->compare_mpot_gpu_global_to_cpu_global(command_queue, &mpot[0], dev_mpot);
      }

      // Sync up cpu array with gpu version to reduce differences due to minor numerical differences
      // otherwise cell count will diverge causing code problems and crashes
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
      }
      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
      }

      if (do_cpu_calc) {
         new_ncells = old_ncells+mesh->rezone_count(mpot);
      }

      if (do_comparison_calc) {
         // This compares ioffset for each block in the calculation
         mesh->compare_ioffset_gpu_global_to_cpu_global(command_queue, old_ncells, &mpot[0], state->dev_ioffset);
      }

      if (do_gpu_sync) {
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
                if (celltype[ic] == REAL_CELL) {
                   mcount += mpot[ic] ? 4 : 1;
                } else {
                   mcount += mpot[ic] ? 2 : 1;
                }
            }
            ioffset[ig] = mtotal;
            mtotal += mcount;
         }
         ezcl_enqueue_write_buffer(command_queue, state->dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int),       &ioffset[0], NULL);
      }

      if (do_comparison_calc) {
         new_ncells = new_ncells_gpu;
      }

      if (do_cpu_calc) {
         int add_ncells = new_ncells - old_ncells;
         state->rezone_all(mesh, mpot, add_ncells);
         mpot.clear();
      }

      //  Resize the mesh, inserting cells where refinement is necessary.
      if (do_gpu_calc) {
         state->gpu_rezone_all(command_queue, mesh, ncells, new_ncells, old_ncells, localStencil);
      }

      //ezcl_device_memory_remove(dev_ioffset);

      if (do_comparison_calc) {
         state->compare_state_gpu_global_to_cpu_global(command_queue,"rezone all",ncycle,ncells);

         mesh->compare_indices_gpu_global_to_cpu_global(command_queue);
      }

      //if (do_gpu_calc) {
      //   int bcount = mesh->gpu_count_BCs(command_queue);
      //}

      if (do_gpu_calc) {
         if (ncells != old_ncells){
            H.resize(ncells);
            U.resize(ncells);
            V.resize(ncells);
            celltype.resize(ncells);
            i.resize(ncells);
            j.resize(ncells);
            level.resize(ncells);
         }
      }

      H_sum = -1.0;

      if (do_comparison_calc) {
         H_sum = state->mass_sum(mesh, enhanced_precision_sum);
         double total_mass = state->gpu_mass_sum(command_queue, mesh, enhanced_precision_sum);
         if (fabs(total_mass - H_sum) > CONSERVATION_EPS) printf("Error: mass sum gpu %f cpu %f\n", total_mass, H_sum);/***/
      }

      mesh->proc.resize(ncells);
      if (icount) {  
         if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
            mesh->calc_spatial_coordinates(0);
         }
         vector<int> index(ncells);
         mesh->partition_cells(numpe, mesh->proc, index, cycle_reorder);
         state->state_reorder(index);
         if (do_gpu_sync) {
            ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_level, CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&level[0],  NULL);
         }
      }
      
      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);
      }
   } // End burst loop

   if (H_sum < 0) {
      H_sum = state->mass_sum(mesh, enhanced_precision_sum);
   }
   printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
      ncycle, deltaT, simTime, ncells, H_sum, H_sum - H_sum_initial);

#ifdef HAVE_GRAPHICS
   if (do_cpu_calc){
      mesh->calc_spatial_coordinates(0);
   }
   if (do_gpu_calc){
      cl_mem dev_x  = ezcl_malloc(NULL, const_cast<char *>("dev_x"),  &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_dx = ezcl_malloc(NULL, const_cast<char *>("dev_dx"), &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_y  = ezcl_malloc(NULL, const_cast<char *>("dev_y"),  &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_dy = ezcl_malloc(NULL, const_cast<char *>("dev_dy"), &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

      mesh->gpu_calc_spatial_coordinates(command_queue, dev_x, dev_dx, dev_y, dev_dy);

      if (do_comparison_calc){
         mesh->compare_coordinates_gpu_global_to_cpu_global(command_queue, dev_x, dev_dx, dev_y, dev_dy, dev_H, &H[0]);
      }

      ezcl_device_memory_remove(dev_x);
      ezcl_device_memory_remove(dev_dx);
      ezcl_device_memory_remove(dev_y);
      ezcl_device_memory_remove(dev_dy);
   }

   set_mysize(ncells);
   set_viewmode(view_mode);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_data(&H[0]);
   set_cell_proc(&mesh->proc[0]);
   set_circle_radius(circle_radius);
   draw_scene();
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      state->output_timing_info(mesh, do_cpu_calc, do_gpu_calc, elapsed_time);

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      exit(0);
   }  //  Complete final output.
   
}

