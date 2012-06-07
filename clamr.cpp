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
#include "l7/l7.h"
#include <mpi.h>

#ifndef DEBUG
#define DEBUG 0
#endif

// Sync is to reduce numerical drift between cpu and gpu
#define DO_SYNC 
//#define DO_COMPARISON

//TODO:  command-line option for OpenGL?
#ifdef DO_COMPARISON
static int do_comparison_calc = 1;
static int do_cpu_calc = 1;
#else
static int do_comparison_calc = 0;
static int do_cpu_calc = 0;
#endif

static int do_gpu_calc = 1;

#ifdef DO_SYNC
static int do_sync = 1;
#else
static int do_sync = 0;
#endif

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
typedef cl_float    cl_real;
typedef cl_float2   cl_real2;
typedef cl_float4   cl_real4;
typedef cl_float8   cl_real8;
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
static Mesh       *mesh;     //  Object containing mesh information; init in grid.cpp::main().
static State      *state;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
static struct timeval tstart;
static cl_event start_write_event, end_write_event;

static cl_context          context                 = NULL;
static cl_command_queue    command_queue           = NULL;
static int compute_device = 0;

static double  H_sum_initial = 0.0;

int main(int argc, char **argv) {
   int ierr;

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   L7_Init(&mype, &numpe, &argc, argv);
   //MPI_Init(&argc, &argv);
   parseInput(argc, argv);

   //MPI_Comm_size(MPI_COMM_WORLD, &numpe);
   //MPI_Comm_rank(MPI_COMM_WORLD, &mype);

   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU, &context, &command_queue, &compute_device, mype);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU, &context, &command_queue, &compute_device, mype);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      L7_Terminate();
      exit(-1);
   }

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   if (special_case) circ_radius = .75;

   mesh_global  = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   mesh_global->init(nx, ny, circ_radius, context, initial_order, special_case, compute_device, do_gpu_calc);
   size_t &ncells_global = mesh_global->ncells;
   state_global = new State(ncells_global, context);
   state_global->init(ncells_global, context, compute_device, do_gpu_calc);
   mesh_global->proc.resize(ncells_global);
   mesh_global->calc_distribution(numpe, mesh_global->proc);
   state_global->fill_circle(mesh_global, circ_radius, 100.0, 5.0);
   
   parallel_in = 1;
   mesh = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   state = new State(mesh->ncells, context);

   size_t &ncells = mesh->ncells;
   int &noffset = mesh->noffset;

   cl_mem &dev_corners_i_global  = mesh_global->dev_corners_i;
   cl_mem &dev_corners_j_global  = mesh_global->dev_corners_j;

   cl_mem &dev_corners_i_local  = mesh->dev_corners_i;
   cl_mem &dev_corners_j_local  = mesh->dev_corners_j;

   vector<int>   &corners_i_global  = mesh_global->corners_i;
   vector<int>   &corners_j_global  = mesh_global->corners_j;

   vector<int>   &corners_i_local  = mesh->corners_i;
   vector<int>   &corners_j_local  = mesh->corners_j;

   vector<int>   &nsizes     = mesh_global->nsizes;
   vector<int>   &ndispl     = mesh_global->ndispl;

   vector<real>  &H_global = state_global->H;
   vector<real>  &U_global = state_global->U;
   vector<real>  &V_global = state_global->V;

   vector<real>  &x_global  = mesh_global->x;
   vector<real>  &dx_global = mesh_global->dx;
   vector<real>  &y_global  = mesh_global->y;
   vector<real>  &dy_global = mesh_global->dy;

   cl_mem &dev_H  = state->dev_H;
   cl_mem &dev_U  = state->dev_U;
   cl_mem &dev_V  = state->dev_V;

   cl_mem &dev_H_global  = state_global->dev_H;
   cl_mem &dev_U_global  = state_global->dev_U;
   cl_mem &dev_V_global  = state_global->dev_V;

   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;

   cl_mem &dev_celltype_global = mesh_global->dev_celltype;
   cl_mem &dev_i_global        = mesh_global->dev_i;
   cl_mem &dev_j_global        = mesh_global->dev_j;
   cl_mem &dev_level_global    = mesh_global->dev_level;

   cl_mem &dev_celltype_new_global = mesh_global->dev_celltype_new;
   cl_mem &dev_i_new_global        = mesh_global->dev_i_new;
   cl_mem &dev_j_new_global        = mesh_global->dev_j_new;
   cl_mem &dev_level_new_global    = mesh_global->dev_level_new;

   cl_mem &dev_celltype_new_local = mesh->dev_celltype_new;
   cl_mem &dev_i_new_local        = mesh->dev_i_new;
   cl_mem &dev_j_new_local        = mesh->dev_j_new;
   cl_mem &dev_level_new_local    = mesh->dev_level_new;

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;

   vector<int>   &celltype_global = mesh_global->celltype;
   vector<int>   &i_global        = mesh_global->i;
   vector<int>   &j_global        = mesh_global->j;
   vector<int>   &level_global    = mesh_global->level;

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

   for (int ip=0; ip<numpe; ip++){
      nsizes[ip] = ncells_global/numpe;
      if (ip < ncells_global%numpe) nsizes[ip]++;
   }

   ndispl[0]=0;
   for (int ip=1; ip<numpe; ip++){
      ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
   }
   noffset=0;
   for (int ip=0; ip<mype; ip++){
     noffset += nsizes[ip];
   }


   size_t corners_size = corners_i_global.size();

   dev_corners_i_global  = ezcl_malloc(&corners_i_global[0],  &corners_size, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
   dev_corners_j_global  = ezcl_malloc(&corners_j_global[0],  &corners_size, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);

   dev_corners_i_local  = ezcl_malloc(&corners_i_local[0],  &corners_size, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
   dev_corners_j_local  = ezcl_malloc(&corners_j_local[0],  &corners_size, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);

   // Distribute level, celltype, H, U, V

   celltype.resize(ncells);
   level.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);

   MPI_Scatterv(&celltype_global[0], &nsizes[0], &ndispl[0], MPI_INT, &celltype[0], nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&level_global[0],    &nsizes[0], &ndispl[0], MPI_INT, &level[0],    nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&i_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &i[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&j_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &j[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);

   dev_celltype  = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_i         = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY,  0);
   dev_j         = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY,  0);
   dev_level     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

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

   state_global->allocate_device_memory(ncells_global);
   state->allocate_device_memory(ncells);

   size_t one = 1;
   state_global->dev_deltaT   = ezcl_malloc(NULL, &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   state->dev_deltaT   = ezcl_malloc(NULL, &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

   dev_celltype_global = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_i_global        = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_j_global        = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_level_global    = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);

   //  Set write buffers for data.
   ezcl_enqueue_write_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&H_global[0],  &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&U_global[0],  NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_V_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&V_global[0],  NULL              );

   ezcl_enqueue_write_buffer(command_queue, dev_celltype_global, CL_FALSE, 0, ncells_global*sizeof(cl_int),  (void *)&celltype_global[0], NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_i_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  (void *)&i_global[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_j_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  (void *)&j_global[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_level_global,    CL_TRUE,  0, ncells_global*sizeof(cl_int),  (void *)&level_global[0],    &end_write_event);
   state_global->gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

   ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_V, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);

   ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],        NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],        NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_level,    CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&level[0],    &end_write_event);
   state->gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

   dev_celltype_new_global = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_i_new_global        = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_j_new_global        = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_level_new_global    = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);

   dev_celltype_new_local = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_i_new_local        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_j_new_local        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_level_new_local    = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state_global->mass_sum(mesh_global, enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   if (mype == 0) {
      printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells_global, H_sum);
   }

   //  Set up grid.

#ifdef HAVE_GRAPHICS
   set_mysize(ncells_global);
   set_viewmode(view_mode);
   set_window(mesh_global->xmin, mesh_global->xmax, mesh_global->ymin, mesh_global->ymax);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water", mype);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_data(&H_global[0]);
   set_cell_proc(&mesh_global->proc[0]);
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
   int icount=0, jcount=0;
   int icount_global=0, jcount_global=0;

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      printf("Can't do this problem with GPU\n");
      exit(0);
   }
   
   //  Initialize state variables for GPU calculation.
   int &mype = mesh->mype;
   int &numpe = mesh->numpe;
   int &noffset = mesh->noffset;

   vector<int>   &nsizes   = mesh_global->nsizes;
   vector<int>   &ndispl   = mesh_global->ndispl;

   vector<int>   &celltype_global = mesh_global->celltype;
   vector<int>   &i_global        = mesh_global->i;
   vector<int>   &j_global        = mesh_global->j;
   vector<int>   &level_global    = mesh_global->level;

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh_global->ncells;
   size_t &ncells           = mesh->ncells;
   size_t &ncells_ghost     = mesh->ncells_ghost;
   int &cell_handle         = mesh->cell_handle;

   vector<real>  &H_global = state_global->H;
   vector<real>  &U_global = state_global->U;
   vector<real>  &V_global = state_global->V;

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

   cl_mem &dev_H_global = state_global->dev_H;
   cl_mem &dev_U_global = state_global->dev_U;
   cl_mem &dev_V_global = state_global->dev_V;

   cl_mem &dev_H  = state->dev_H;
   cl_mem &dev_U  = state->dev_U;
   cl_mem &dev_V  = state->dev_V;

   cl_mem &dev_mpot_global     = mesh_global->dev_mpot;

   cl_mem &dev_mpot     = mesh->dev_mpot;


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
   for (int nburst = 0; nburst < outputInterval && ncycle <= niter; nburst++, ncycle++) {

      // To reduce drift in solution
      if (do_comparison_calc){
         if (do_sync) {
            ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);

            ezcl_enqueue_read_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&H_global[0],  NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&U_global[0],  NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_V_global, CL_TRUE,  0, ncells_global*sizeof(cl_real),  (void *)&V_global[0],  NULL);
         }
      }

      size_t local_work_size_global  = MIN(ncells_global, TILE_SIZE);
      size_t global_work_size_global = ((ncells_global+local_work_size_global - 1) /local_work_size_global) * local_work_size_global;
      size_t block_size_global     = global_work_size_global/local_work_size_global;

      size_t local_work_size  = MIN(ncells, TILE_SIZE);
      size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
      size_t block_size     = global_work_size/local_work_size;

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;
      old_ncells_global = ncells_global;

      //  Calculate the real time step for the current discrete time step.
      deltaT = state->gpu_set_timestep(command_queue, mesh, sigma);
      simTime += deltaT;

      if (do_comparison_calc) {
         double deltaT_cpu_global = state_global->set_timestep(mesh_global, g, sigma);
         double deltaT_cpu_local = state->set_timestep(mesh, g, sigma);

         double deltaT_gpu_global = state_global->gpu_set_timestep(command_queue, mesh_global, sigma);

         //  Compare time step values and pass deltaT in to the kernel.
         int iflag = 0;
         if (fabs(deltaT_cpu_local - deltaT_cpu_global) > .000001) iflag = 1;
         if (fabs(deltaT - deltaT_gpu_global) > .000001) iflag = 1;
         if (fabs(deltaT_gpu_global - deltaT_cpu_global) > .000001) iflag = 1;
         if (fabs(deltaT - deltaT_cpu_local) > .000001) iflag = 1;
         if (iflag) {
            printf("Error with deltaT calc --- cpu_local %lf cpu_global %lf gpu_local %lf gpu_global %lf\n",
               deltaT_cpu_local, deltaT_cpu_global, deltaT, deltaT_gpu_global);
         }
      }

      mesh->gpu_calc_neighbors_local(command_queue);

      if (do_comparison_calc) {
         mesh_global->calc_neighbors();
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

         mesh_global->gpu_calc_neighbors(command_queue);

         mesh->compare_neighbors_all_to_gpu_local(command_queue, mesh_global, &nsizes[0], &ndispl[0]);

         mesh->partition_measure();
      }

      // Currently not working -- may need to be earlier?
      //if (do_comparison_calc && ! mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Need ghost cells for this routine
      if (do_comparison_calc) {
        state_global->apply_boundary_conditions(mesh_global);
        state->apply_boundary_conditions(mesh);
      }

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->gpu_calc_finite_difference_local(command_queue, mesh, deltaT);

      if (do_comparison_calc) {
         state_global->calc_finite_difference(mesh_global, deltaT);
         state->calc_finite_difference_local(mesh, deltaT);
      
         state_global->gpu_calc_finite_difference(command_queue, mesh_global, deltaT);
      
         state->compare_state_all_to_gpu_local(command_queue, state_global, ncells, ncells_global, mype, ncycle, &nsizes[0], &ndispl[0]);
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      if (do_comparison_calc) {
         state_global->remove_boundary_cells(mesh_global);
         state->remove_boundary_cells(mesh);
      }
      
      //  Check for NANs.
      for (uint ic=0; ic<ncells; ic++) {
         if (isnan(H[ic]))
         {  printf("Got a NAN on cell %d cycle %d\n",ic,ncycle);
            H[ic]=0.0;
            sleep(100);
            //  Release kernels and finalize the OpenCL elements.
            ezcl_finalize();
            L7_Terminate();
            exit(-1); }
      }  //  Complete NAN check.
      
      size_t result_size = 1;
      cl_mem dev_result  = ezcl_malloc(NULL, &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_ioffset    = ezcl_malloc(NULL, &block_size, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);
      dev_mpot     = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);

      vector<int>      ioffset(block_size);
      vector<int>      ioffset_global(block_size_global);

      cl_mem dev_ioffset_global    = ezcl_malloc(NULL, &block_size_global, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);

      cl_mem dev_result_global = NULL;
      dev_mpot_global     = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
      dev_result_global  = ezcl_malloc(NULL, &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
 
      // Should not need this call! -- 
      state->gpu_calc_refine_potential_local(command_queue, mesh, dev_result, dev_ioffset);

      if (do_comparison_calc) {
         mpot.resize(ncells_ghost);
         mpot_global.resize(ncells_global);
         state_global->calc_refine_potential(mesh_global, mpot_global, icount_global, jcount_global);
         state->calc_refine_potential(mesh, mpot, icount, jcount);

         state_global->gpu_calc_refine_potential(command_queue, mesh_global, dev_result_global, dev_ioffset_global);
      
         int icount_test;
         MPI_Allreduce(&icount, &icount_test, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		 if (icount_test != icount_global) {
		    printf("%d: DEBUG -- icount is %d icount_test %d icount_global is %d\n",mype,icount,icount_test,icount_global);
         }

         mesh->compare_mpot_all_to_gpu_local(command_queue, &mpot[0], &mpot_global[0], dev_mpot, dev_mpot_global, ncells_global, &nsizes[0], &ndispl[0], ncycle);
      }

      // Sync up cpu array with gpu version to reduce differences due to minor numerical differences
      // otherwise cell count will diverge causing code problems and crashes
      if (do_comparison_calc){
         if (do_sync) {
            ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_mpot_global, CL_TRUE,  0, ncells_global*sizeof(cl_int), &mpot_global[0], NULL);
         }
      }

      if (do_comparison_calc) {
         mesh->compare_ioffset_all_to_gpu_local(command_queue, old_ncells, old_ncells_global, block_size, block_size_global, &mpot[0], &mpot_global[0], dev_ioffset, dev_ioffset_global, &ioffset[0], &ioffset_global[0], &celltype_global[0]);
      }

      if (do_comparison_calc) {
         new_ncells_global = old_ncells_global+mesh_global->rezone_count(mpot_global);
         new_ncells = old_ncells+mesh->rezone_count(mpot);
         //printf("new_ncells %d old_ncells %d rezone_count %d\n",new_ncells, old_ncells, mesh->rezone_count(mpot)) ;
         //printf("new_ncells_global %d old_ncells_global %d rezone_count_global %d\n",new_ncells_global, old_ncells_global, mesh_global->rezone_count(mpot_global));
      }

      if (do_comparison_calc) {
         int result;
         ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, 1*sizeof(cl_int),       &result, NULL);
         if (new_ncells != result) printf("%d: DEBUG new_ncells not correct %ld %d\n",mype,new_ncells,result);
         new_ncells = result;
         //printf("Result is %d\n",result[0]);

         ezcl_enqueue_read_buffer(command_queue, dev_result_global, CL_TRUE, 0, 1*sizeof(cl_int),       &result, NULL);
         if (new_ncells_global != result) printf("%d: DEBUG new_ncells_global not correct %ld %d\n",mype,new_ncells_global,result);
         new_ncells_global = result;
      }

      int result;
      ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, 1*sizeof(cl_int),       &result, NULL);
      new_ncells = result;
      //printf("Result is %d %d %d\n",result, ncells,__LINE__);

      ezcl_device_memory_remove(dev_result);

      if (do_comparison_calc) {
         ezcl_device_memory_remove(dev_result_global);
      }

      int add_ncells_global=0, add_ncells=0;
      if (do_comparison_calc) {
         add_ncells_global = new_ncells_global - old_ncells_global;
         add_ncells = new_ncells - old_ncells;
         //printf("%d: DEBUG add %d new %d old %d\n",mype,add_ncells,new_ncells,old_ncells);
      }

      //  Resize the mesh, inserting cells where refinement is necessary.
      state->gpu_rezone_all_local(command_queue, mesh, old_ncells, new_ncells, old_ncells, localStencil, dev_ioffset);

      if (do_comparison_calc) {
         state_global->rezone_all(mesh_global, mpot_global, add_ncells_global);
         state->rezone_all(mesh, mpot, add_ncells);
         mpot_global.clear();
         mpot.clear();

         state_global->gpu_rezone_all(command_queue, mesh_global, ncells_global, new_ncells_global, old_ncells_global, localStencil, dev_ioffset_global);
      }

      ncells = new_ncells;
      mesh->ncells = ncells;
      //printf("%d: DEBUG ncells is %d new_ncells %d old_ncells %d ncells_global %d\n",mype, ncells, new_ncells, old_ncells, ncells_global);
      MPI_Allgather(&ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }
      noffset=0;
      for (int ip=0; ip<mype; ip++){
        noffset += nsizes[ip];
      }
      MPI_Allreduce(&ncells, &ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if (do_comparison_calc){
         state->compare_state_all_to_gpu_local(command_queue, state_global, ncells, ncells_global, mype, ncycle, &nsizes[0], &ndispl[0]);

         mesh->compare_indices_all_to_gpu_local(command_queue, mesh_global, ncells_global, &nsizes[0], &ndispl[0], ncycle);
      } // do_comparison_calc

#ifdef XXX // not rewritten yet
      mesh->gpu_count_BCs(command_queue, block_size, local_work_size, global_work_size, dev_ioffset);
#endif

      if (ncells != old_ncells){
         H.resize(ncells);
         U.resize(ncells);
         V.resize(ncells);
         level.resize(ncells);
         i.resize(ncells);
         j.resize(ncells);
         celltype.resize(ncells);
      }

      if (do_comparison_calc) {
         if (ncells_global != old_ncells_global){
            H_global.resize(ncells_global);
            U_global.resize(ncells_global);
            V_global.resize(ncells_global);
            level_global.resize(ncells_global);
            i_global.resize(ncells_global);
            j_global.resize(ncells_global);
            celltype_global.resize(ncells_global);
         }
      }

      ioffset.clear();
      ioffset_global.clear();

      ezcl_device_memory_remove(dev_ioffset);
      ezcl_device_memory_remove(dev_ioffset_global);

      H_sum = -1.0;

      if (do_comparison_calc) {
         double cpu_mass_sum = state_global->mass_sum(mesh_global, enhanced_precision_sum);
         double cpu_mass_sum_local = state->mass_sum(mesh, enhanced_precision_sum);
         double gpu_mass_sum = state_global->gpu_mass_sum(command_queue, mesh_global, enhanced_precision_sum);
         H_sum = state->gpu_mass_sum_local(command_queue, mesh, enhanced_precision_sum);
         int iflag = 0;
         if (fabs(cpu_mass_sum_local - cpu_mass_sum) > CONSERVATION_EPS) iflag = 1;
         //if (fabs(H_sum - gpu_mass_sum) > CONSERVATION_EPS) iflag = 1;
         //if (fabs(gpu_mass_sum - cpu_mass_sum) > CONSERVATION_EPS) iflag = 1;
         if (fabs(H_sum - cpu_mass_sum_local) > CONSERVATION_EPS) iflag = 1;

         if (iflag) {
            printf("Error with mass sum calculation -- cpu_mass_sum_local %lf cpu_mass_sum %lf gpu_mass_sum_local %lf gpu_mass_sum %lf\n",
                    cpu_mass_sum_local, cpu_mass_sum, H_sum, gpu_mass_sum);
         }
      }

      mesh_global->proc.resize(ncells_global);
      mesh->proc.resize(ncells);
      if (icount) {
         vector<int> index(ncells);
         vector<int> index_global(ncells_global);
         mesh_global->partition_cells(numpe, mesh_global->proc, index_global, cycle_reorder);
         mesh->partition_cells(numpe, mesh->proc, index, cycle_reorder);
         //state->state_reorder(index);
      }
   }  //  End burst loop

   if (H_sum < 0) {
      H_sum = state->gpu_mass_sum_local(command_queue, mesh, enhanced_precision_sum);
   }
   if (mype == 0){
      printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

#ifdef HAVE_GRAPHICS
   cl_mem dev_x  = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_dx = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_y  = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   cl_mem dev_dy = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

   mesh->gpu_calc_spatial_coordinates(command_queue, dev_x, dev_dx, dev_y, dev_dy);

   x.resize(ncells_ghost);
   dx.resize(ncells_ghost);
   y.resize(ncells_ghost);
   dy.resize(ncells_ghost);
   H.resize(ncells_ghost);

   ezcl_enqueue_read_buffer(command_queue, dev_x,  CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&x[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dx, CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&dx[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,  CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&y[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy, CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&dy[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,  CL_TRUE,  0, ncells*sizeof(cl_real), (void *)&H[0],  NULL);

   ezcl_device_memory_remove(dev_x);
   ezcl_device_memory_remove(dev_dx);
   ezcl_device_memory_remove(dev_y);
   ezcl_device_memory_remove(dev_dy);

   x_global.resize(ncells_global);
   dx_global.resize(ncells_global);
   y_global.resize(ncells_global);
   dy_global.resize(ncells_global);
   H_global.resize(ncells_global);

   MPI_Allgatherv(&x[0],  ncells, MPI_C_REAL, &x_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&dx[0], ncells, MPI_C_REAL, &dx_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&y[0],  ncells, MPI_C_REAL, &y_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&dy[0], ncells, MPI_C_REAL, &dy_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&H[0],  ncells, MPI_C_REAL, &H_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);

   set_mysize(ncells_global);
   set_viewmode(view_mode);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_data(&H_global[0]);
   set_cell_proc(&mesh_global->proc[0]);
   set_circle_radius(circle_radius);
   draw_scene();
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      //  Release kernels and finalize the OpenCL elements.
      ezcl_finalize();
      
      if (do_comparison_calc){
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

