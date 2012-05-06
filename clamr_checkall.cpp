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
#define MAILBOX 1

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
#include "reduce.h"
#include "state.h"
#include "timer/timer.h"
#include "l7/l7.h"

//#undef HAVE_OPENGL
#ifdef HAVE_OPENGL
#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// Sync is to reduce numerical drift between cpu and gpu
#define DO_SYNC 
#define DO_COMPARISON
//#define DO_CPU
//#define DO_GPU

//TODO:  command-line option for OpenGL?
#ifdef DO_COMPARISON
#define DO_CPU
#define DO_GPU
int do_comparison_calc = 1;
#else
int do_comparison_calc = 0;
#endif

#ifdef DO_CPU
int do_cpu_calc = 1;
#else
int do_cpu_calc = 0;
#endif

#ifdef DO_GPU
int do_gpu_calc = 1;
#else
int do_gpu_calc = 0;
#endif

#ifdef DO_SYNC
int do_sync = 1;
#else
int do_sync = 0;
#endif
int do_gpu_sync = 0;

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

double circle_radius=-1.0;

int view_mode = 0;

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
Mesh       *mesh_global;    //  Object containing mesh information; init in grid.cpp::main().
State      *state_global;   //  Object containing state information corresponding to mesh; init in grid.cpp::main().
Mesh       *mesh_local;     //  Object containing mesh information; init in grid.cpp::main().
State      *state_local;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
struct timeval tstart, tstop, tresult;
struct timeval tstart_cpu;
cl_event start_write_event, end_write_event,
         count_BCs_stage1_event,
         count_BCs_stage2_event;
double   cpu_time_start,
         cpu_time_end;
long     gpu_time_start,
         gpu_time_end,
         gpu_time_count_BCs          = 0,
         gpu_time_count_BCs_parallel = 0;

#ifdef HAVE_OPENCL
cl_context          context                 = NULL;
cl_command_queue    command_queue           = NULL;
cl_kernel           kernel_count_BCs        = NULL;
#endif

int main(int argc, char **argv) {
#ifdef HAVE_OPENCL
   int ierr;
#endif

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   L7_Init(&mype, &numpe, &argc, argv);
   //MPI_Init(&argc, &argv);
   parseInput(argc, argv);

#ifdef HAVE_MPI
   //MPI_Comm_size(MPI_COMM_WORLD, &numpe);
   //MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#else
   numpe = 16;
#endif

#ifdef HAVE_OPENCL
   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU, &context, &command_queue, mype);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU, &context, &command_queue, mype);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      L7_Terminate();
      exit(-1);
   }
#endif

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   if (special_case) circ_radius = .75;

   mesh_global  = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   mesh_global->init(nx, ny, circ_radius, context, initial_order, special_case, do_gpu_calc);
   size_t &ncells_global = mesh_global->ncells;
   state_global = new State(ncells_global, context);
   state_global->init(ncells_global, context, do_gpu_calc);
   mesh_global->proc.resize(ncells_global);
   mesh_global->calc_distribution(numpe, mesh_global->proc);
   state_global->fill_circle(mesh_global, circ_radius, 100.0, 5.0);
   
#ifdef HAVE_OPENCL
   init_kernel_2stage_sum(context);
   init_kernel_2stage_sum_int(context);
   if (! mesh_global->have_boundary){
     kernel_count_BCs       = ezcl_create_kernel(context, "wave_kern.cl",      "count_BCs_cl",             0);
   }
#endif
   
   parallel_in = 1;
   mesh_local = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   state_local = new State(mesh_local->ncells, context);

   size_t &ncells = mesh_local->ncells;
   int &noffset = mesh_local->noffset;

   cl_mem &dev_corners_i_global  = mesh_global->dev_corners_i;
   cl_mem &dev_corners_j_global  = mesh_global->dev_corners_j;

   cl_mem &dev_corners_i_local  = mesh_local->dev_corners_i;
   cl_mem &dev_corners_j_local  = mesh_local->dev_corners_j;

   vector<int>   &corners_i_global  = mesh_global->corners_i;
   vector<int>   &corners_j_global  = mesh_global->corners_j;

   vector<int>   &corners_i_local  = mesh_local->corners_i;
   vector<int>   &corners_j_local  = mesh_local->corners_j;

   vector<int>   &nsizes     = mesh_global->nsizes;
   vector<int>   &ndispl     = mesh_global->ndispl;

   vector<real>  &H_global = state_global->H;
   vector<real>  &U_global = state_global->U;
   vector<real>  &V_global = state_global->V;

   vector<real>  &x_global  = mesh_global->x;
   vector<real>  &dx_global = mesh_global->dx;
   vector<real>  &y_global  = mesh_global->y;
   vector<real>  &dy_global = mesh_global->dy;

   cl_mem &dev_H  = state_local->dev_H;
   cl_mem &dev_U  = state_local->dev_U;
   cl_mem &dev_V  = state_local->dev_V;

   cl_mem &dev_H_global  = state_global->dev_H;
   cl_mem &dev_U_global  = state_global->dev_U;
   cl_mem &dev_V_global  = state_global->dev_V;

   cl_mem &dev_celltype = mesh_local->dev_celltype;
   cl_mem &dev_i        = mesh_local->dev_i;
   cl_mem &dev_j        = mesh_local->dev_j;
   cl_mem &dev_level    = mesh_local->dev_level;

   cl_mem &dev_celltype_global = mesh_global->dev_celltype;
   cl_mem &dev_i_global        = mesh_global->dev_i;
   cl_mem &dev_j_global        = mesh_global->dev_j;
   cl_mem &dev_level_global    = mesh_global->dev_level;

   cl_mem &dev_celltype_new_global = mesh_global->dev_celltype_new;
   cl_mem &dev_i_new_global        = mesh_global->dev_i_new;
   cl_mem &dev_j_new_global        = mesh_global->dev_j_new;
   cl_mem &dev_level_new_global    = mesh_global->dev_level_new;

   cl_mem &dev_celltype_new_local = mesh_local->dev_celltype_new;
   cl_mem &dev_i_new_local        = mesh_local->dev_i_new;
   cl_mem &dev_j_new_local        = mesh_local->dev_j_new;
   cl_mem &dev_level_new_local    = mesh_local->dev_level_new;

   vector<int>   &celltype = mesh_local->celltype;
   vector<int>   &i        = mesh_local->i;
   vector<int>   &j        = mesh_local->j;
   vector<int>   &level    = mesh_local->level;

   vector<int>   &celltype_global = mesh_global->celltype;
   vector<int>   &i_global        = mesh_global->i;
   vector<int>   &j_global        = mesh_global->j;
   vector<int>   &level_global    = mesh_global->level;

   vector<real> &H = state_local->H;
   vector<real> &U = state_local->U;
   vector<real> &V = state_local->V;

   vector<real> &x  = mesh_local->x;
   vector<real> &dx = mesh_local->dx;
   vector<real> &y  = mesh_local->y;
   vector<real> &dy = mesh_local->dy;

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
   state_local->allocate_device_memory(ncells);

   size_t one = 1;
   state_global->dev_deltaT   = ezcl_malloc(NULL, &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
   state_local->dev_deltaT   = ezcl_malloc(NULL, &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

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
   state_local->gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

   dev_celltype_new_global = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_i_new_global        = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_j_new_global        = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_level_new_global    = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);

   dev_celltype_new_local = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_i_new_local        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_j_new_local        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_level_new_local    = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);

#ifdef HAVE_OPENGL
   set_cell_data(&H_global[0]);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_mysize(ncells_global);
   set_viewmode(view_mode);
   set_window(mesh_global->xmin, mesh_global->xmax, mesh_global->ymin, mesh_global->ymax);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water");
   glutIdleFunc(&do_calc);
   glutMainLoop();
#else
   for (int it = 0; it < 10000000; it++) {
      do_calc();
   }
#endif
   
   return 0;
}

int     n       = -1;
double  simTime = 0.0;
double  H_sum_initial = 0.0;

extern "C" void do_calc(void)
{  double g     = 9.80;
   double sigma = 0.95; 
   int icount, jcount, bcount;
   int icount_global, jcount_global, bcount_global;

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      do_sync = 0;
      do_gpu_sync = 1;
   }
   
#ifdef HAVE_OPENCL
   cl_mem dev_mpot  = NULL;
   cl_mem dev_mpot_global  = NULL;
#endif

   //  Initialize state variables for GPU calculation.
   int &mype = mesh_local->mype;
   int &numpe = mesh_local->numpe;
   int &noffset = mesh_local->noffset;

   vector<int>   &nsizes   = mesh_global->nsizes;
   vector<int>   &ndispl   = mesh_global->ndispl;

   vector<int>   &celltype_global = mesh_global->celltype;
   vector<int>   &i_global        = mesh_global->i;
   vector<int>   &j_global        = mesh_global->j;
   vector<int>   &index_global    = mesh_global->index;
   vector<int>   &level_global    = mesh_global->level;
   vector<int>   &nlft_global     = mesh_global->nlft;
   vector<int>   &nrht_global     = mesh_global->nrht;
   vector<int>   &nbot_global     = mesh_global->nbot;
   vector<int>   &ntop_global     = mesh_global->ntop;

   vector<int>   &celltype = mesh_local->celltype;
   vector<int>   &i        = mesh_local->i;
   vector<int>   &j        = mesh_local->j;
   vector<int>   &index    = mesh_local->index;
   vector<int>   &level    = mesh_local->level;
   vector<int>   &nlft     = mesh_local->nlft;
   vector<int>   &nrht     = mesh_local->nrht;
   vector<int>   &nbot     = mesh_local->nbot;
   vector<int>   &ntop     = mesh_local->ntop;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh_global->ncells;
   size_t &ncells           = mesh_local->ncells;
   size_t &ncells_ghost     = mesh_local->ncells_ghost;
   int &cell_handle         = mesh_local->cell_handle;

   vector<real>  &H_global = state_global->H;
   vector<real>  &U_global = state_global->U;
   vector<real>  &V_global = state_global->V;

   vector<real>  &H = state_local->H;
   vector<real>  &U = state_local->U;
   vector<real>  &V = state_local->V;

   vector<real>  &x  = mesh_local->x;
   vector<real>  &dx = mesh_local->dx;
   vector<real>  &y  = mesh_local->y;
   vector<real>  &dy = mesh_local->dy;

   vector<real>  &x_global  = mesh_global->x;
   vector<real>  &dx_global = mesh_global->dx;
   vector<real>  &y_global  = mesh_global->y;
   vector<real>  &dy_global = mesh_global->dy;

   cl_mem &dev_levdx    = mesh_local->dev_levdx;
   cl_mem &dev_levdy    = mesh_local->dev_levdy;

   cl_mem &dev_levibeg  = mesh_local->dev_levibeg;
   cl_mem &dev_leviend  = mesh_local->dev_leviend;
   cl_mem &dev_levjbeg  = mesh_local->dev_levjbeg;
   cl_mem &dev_levjend  = mesh_local->dev_levjend;

   cl_mem &dev_H_global = state_global->dev_H;
   cl_mem &dev_U_global = state_global->dev_U;
   cl_mem &dev_V_global = state_global->dev_V;

   cl_mem &dev_deltaT   = state_local->dev_deltaT;

   cl_mem &dev_H  = state_local->dev_H;
   cl_mem &dev_U  = state_local->dev_U;
   cl_mem &dev_V  = state_local->dev_V;

   cl_mem &dev_celltype_global = mesh_global->dev_celltype;
   cl_mem &dev_i_global        = mesh_global->dev_i;
   cl_mem &dev_j_global        = mesh_global->dev_j;
   cl_mem &dev_level_global    = mesh_global->dev_level;
   cl_mem &dev_nlft_global     = mesh_global->dev_nlft;
   cl_mem &dev_nrht_global     = mesh_global->dev_nrht;
   cl_mem &dev_nbot_global     = mesh_global->dev_nbot;
   cl_mem &dev_ntop_global     = mesh_global->dev_ntop;

   cl_mem &dev_celltype = mesh_local->dev_celltype;
   cl_mem &dev_i        = mesh_local->dev_i;
   cl_mem &dev_j        = mesh_local->dev_j;
   cl_mem &dev_level    = mesh_local->dev_level;
   cl_mem &dev_nlft     = mesh_local->dev_nlft;
   cl_mem &dev_nrht     = mesh_local->dev_nrht;
   cl_mem &dev_nbot     = mesh_local->dev_nbot;
   cl_mem &dev_ntop     = mesh_local->dev_ntop;

   cl_mem &dev_celltype_new_global = mesh_global->dev_celltype_new;
   cl_mem &dev_i_new_global        = mesh_global->dev_i_new;
   cl_mem &dev_j_new_global        = mesh_global->dev_j_new;
   cl_mem &dev_level_new_global    = mesh_global->dev_level_new;
   
   //  Kahan-type enhanced precision sum implementation.
   if (n < 0)
   {
      double H_sum = state_global->mass_sum(mesh_global, enhanced_precision_sum);
      if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
      H_sum_initial = H_sum;
      n++;
      
      //  Set up grid.
#ifdef HAVE_OPENGL
      set_mysize(ncells_global);
      set_viewmode(view_mode);
      set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
      set_cell_data(&H_global[0]);
      set_cell_proc(&mesh_global->proc[0]);
      set_circle_radius(circle_radius);
      DrawGLScene();
      if (verbose) sleep(5);
#endif
      gettimeofday(&tstart, NULL);
      return;
   }

   //  Set flag to show mesh results rather than domain decomposition.
   view_mode = 1;
   
   //  Clear superposition of circle on grid output.
   if (n > 2)
   {  circle_radius = -1.0; }
   
   //  Output final results and timing information.
   if (n > niter) {
      //free_display();
      
      //  Get overall program timing.
      gettimeofday(&tstop, NULL);
      tresult.tv_sec = tstop.tv_sec - tstart.tv_sec;
      tresult.tv_usec = tstop.tv_usec - tstart.tv_usec;
      double elapsed_time = (double)tresult.tv_sec + (double)tresult.tv_usec*1.0e-6;
      
#ifdef HAVE_OPENCL
      //  Release kernels and finalize the OpenCL elements.
      ezcl_finalize();
      
      state_global->output_timing_info(mesh_global, do_cpu_calc, do_gpu_calc, gpu_time_count_BCs, elapsed_time);
      state_local->output_timing_info(mesh_local, do_cpu_calc, do_gpu_calc, gpu_time_count_BCs_parallel, elapsed_time);

      mesh_local->print_partition_measure();
      mesh_local->print_calc_neighbor_type();
      mesh_local->print_partition_type();
#endif
      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
   vector<int>     mpot;
   vector<int>     mpot_global;
   
   if (DEBUG) {
      //if (mype == 0) mesh->print();

      char filename[10];
      sprintf(filename,"out%1d",mype);
      mesh_local->fp=fopen(filename,"w");

      //mesh->print_local();
   }

   //  Set write buffers for data. -- Probably just needed to reduce the accumulated difference between CPU/GPU

   //  Main loop.
   for (int iburst = 0; iburst < outputInterval; iburst++)
   {  if (n > niter) break;

      // To reduce drift in solution
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);

         ezcl_enqueue_read_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&H_global[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&U_global[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V_global, CL_TRUE,  0, ncells_global*sizeof(cl_real),  (void *)&V_global[0],  NULL);
      }

      size_t local_work_size_global  = MIN(ncells_global, TILE_SIZE);
      size_t global_work_size_global = ((ncells_global+local_work_size_global - 1) /local_work_size_global) * local_work_size_global;
      size_t block_size_global     = global_work_size_global/local_work_size_global;

      size_t local_work_size  = MIN(ncells, TILE_SIZE);
      size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
      size_t block_size     = global_work_size/local_work_size;

      //  Define basic domain decomposition parameters for GPU.
      size_t old_ncells = ncells;
      size_t old_ncells_global = ncells_global;
      size_t new_ncells = 0;
      size_t new_ncells_global = 0;

      //  Calculate the real time step for the current discrete time step.
      double deltaT_cpu, deltaT_cpu_local;
      if (do_cpu_calc) {
         deltaT_cpu = state_global->set_timestep(mesh_global, g, sigma);
         deltaT_cpu_local = state_local->set_timestep(mesh_local, g, sigma);
      }  //  Complete CPU timestep calculation.

      double deltaT_gpu, deltaT_gpu_local;
      if (do_gpu_calc) {
         deltaT_gpu = state_global->gpu_set_timestep(command_queue, mesh_global, sigma);
         deltaT_gpu_local = state_local->gpu_set_timestep(command_queue, mesh_local, sigma);
      }  //  Complete GPU calculation.

      //  Compare time step values and pass deltaT in to the kernel.
      if (do_comparison_calc) {
         int iflag = 0;
         if (fabs(deltaT_cpu_local - deltaT_cpu) > .000001) iflag = 1;
         if (fabs(deltaT_gpu_local - deltaT_gpu) > .000001) iflag = 1;
         if (fabs(deltaT_gpu - deltaT_cpu) > .000001) iflag = 1;
         if (fabs(deltaT_gpu_local - deltaT_cpu_local) > .000001) iflag = 1;
         if (iflag) {
            printf("Error with deltaT calc --- cpu_local %lf cpu_global %lf gpu_local %lf gpu_global %lf\n",
               deltaT_cpu_local, deltaT_cpu, deltaT_gpu_local, deltaT_gpu);
         }
      }
      
      double deltaT = (do_gpu_calc) ? deltaT_gpu_local : deltaT_cpu_local;

      if (do_cpu_calc) {
         mesh_global->calc_neighbors();
         mesh_local->calc_neighbors_local();

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
      }

      if (do_gpu_calc) {
         mesh_global->gpu_calc_neighbors(command_queue);
         mesh_local->gpu_calc_neighbors_local(command_queue);
      }

      if (do_comparison_calc) {
         // Checking CPU parallel to CPU global
         vector<int> Test(ncells_ghost);
         for(int ic=0; ic<ncells; ic++){
            Test[ic] = mype*1000 +ic;
         }
         L7_Update(&Test[0], L7_INT, cell_handle);

         vector<int> Test_global(ncells_global);
         MPI_Allgatherv(&Test[0], nsizes[mype], MPI_INT, &Test_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         vector<int> Test_check(ncells);
         vector<int> Test_check_global(ncells_global);

         // ==================== check left value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nlft[ic]];
            //if (mype == 1 && ic==0) printf("%d: nlft check for ic 0 is %d\n",mype,nlft[0]);
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            //if (Test_global[nlft_global[ic]] != Test_check_global[ic]) {
               //if (mype == 0) printf("%d: Error with nlft for cell %d -- nlft %d global %d check %d\n",mype,ic,nlft_global[ic],Test_global[nlft_global[ic]],Test_check_global[ic]);
            //}
         }
         
         // ==================== check left left value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nlft[nlft[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[nlft_global[nlft_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with nlft nlft for cell %5d -- nlftg %5d nlftg nlftg %5d global %5d\n",
                  mype,ic,nlft_global[ic],nlft_global[nlft_global[ic]],Test_global[nlft_global[nlft_global[ic]]]);
               printf("%d:                         check %5d -- nlftl %5d nlftl nlftl %5d check  %5d\n",
                  mype,ic,nlft[ic],nlft[nlft[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check right value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nrht[ic]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[nrht_global[ic]] != Test_check_global[ic]) {
               if (mype == 0) printf("%d: Error with nrht for cell %d -- %d %d\n",mype,ic,Test_global[nrht_global[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check right right value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nrht[nrht[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[nrht_global[nrht_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with nrht nrht for cell %5d -- nrhtg %5d nrhtg nrhtg %5d global %5d\n",
                  mype,ic,nrht_global[ic],nrht_global[nrht_global[ic]],Test_global[nrht_global[nrht_global[ic]]]);
               printf("%d:                         check %5d -- nrhtl %5d nrhtl nrhtl %5d check  %5d\n",
                  mype,ic,nrht[ic],nrht[nrht[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check bottom value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nbot[ic]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[nbot_global[ic]] != Test_check_global[ic]) {
               if (mype == 0) printf("%d: Error with nbot for cell %d -- %d %d\n",mype,ic,Test_global[nbot_global[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check bottom bottom value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nbot[nbot[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[nbot_global[nbot_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with nbot nbot for cell %5d -- nbotg %5d nbotg nbotg %5d global %5d\n",
                  mype,ic,nbot_global[ic],nbot_global[nbot_global[ic]],Test_global[nbot_global[nbot_global[ic]]]);
               printf("%d:                         check %5d -- nbotl %5d nbotl nbotl %5d check  %5d\n",
                  mype,ic,nbot[ic],nbot[nbot[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check top value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[ntop[ic]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[ntop_global[ic]] != Test_check_global[ic]) {
               if (mype == 0) printf("%d: Error with ntop for cell %d -- %d %d\n",mype,ic,Test_global[ntop_global[ic]],Test_check_global[ic]);
            }
         }

         // ==================== check top top value ====================
         for (int ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[ntop[ntop[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (int ic=0; ic<ncells_global; ic++){
            if (Test_global[ntop_global[ntop_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with ntop ntop for cell %5d -- ntopg %5d ntopg ntopg %5d global %5d\n",
                  mype,ic,ntop_global[ic],ntop_global[ntop_global[ic]],Test_global[ntop_global[ntop_global[ic]]]);
               printf("%d:                         check %5d -- ntopl %5d ntopl ntopl %5d check  %5d\n",
                  mype,ic,ntop[ic],ntop[ntop[ic]],Test_check_global[ic]);
            }
         }
         
         // checking gpu results
         vector<int> nlft_check(ncells_ghost);
         vector<int> nrht_check(ncells_ghost);
         vector<int> nbot_check(ncells_ghost);
         vector<int> ntop_check(ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nlft_check[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nrht_check[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nbot_check[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells_ghost*sizeof(cl_int),  &ntop_check[0], NULL);

         for (int ic=0; ic<ncells_ghost; ic++){
            if (nlft[ic] != nlft_check[ic]) printf("%d: Error with gpu calculated nlft for cell %d nlft %d check %d\n",mype,ic,nlft[ic],nlft_check[ic]);
            if (nrht[ic] != nrht_check[ic]) printf("%d: Error with gpu calculated nrht for cell %d nrht %d check %d\n",mype,ic,nrht[ic],nrht_check[ic]);
            if (nbot[ic] != nbot_check[ic]) printf("%d: Error with gpu calculated nbot for cell %d nbot %d check %d\n",mype,ic,nbot[ic],nbot_check[ic]);
            if (ntop[ic] != ntop_check[ic]) printf("%d: Error with gpu calculated ntop for cell %d ntop %d check %d\n",mype,ic,ntop[ic],ntop_check[ic]);
         }

         // Now check the global neighbors
         vector<int> nlft_global_check(ncells_global);
         vector<int> nrht_global_check(ncells_global);
         vector<int> nbot_global_check(ncells_global);
         vector<int> ntop_global_check(ncells_global);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft_global, CL_FALSE, 0, ncells_global*sizeof(cl_int), &nlft_global_check[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht_global, CL_FALSE, 0, ncells_global*sizeof(cl_int), &nrht_global_check[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot_global, CL_FALSE, 0, ncells_global*sizeof(cl_int), &nbot_global_check[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop_global, CL_TRUE,  0, ncells_global*sizeof(cl_int), &ntop_global_check[0], NULL);

         for (int ic=0; ic<ncells_global; ic++){
            if (nlft_global[ic] != nlft_global_check[ic]) printf("%d: Error with gpu calculated nlft for cell %d nlft %d check %d\n",mype,ic,nlft_global[ic],nlft_global_check[ic]);
            if (nrht_global[ic] != nrht_global_check[ic]) printf("%d: Error with gpu calculated nrht for cell %d nrht %d check %d\n",mype,ic,nrht_global[ic],nrht_global_check[ic]);
            if (nbot_global[ic] != nbot_global_check[ic]) printf("%d: Error with gpu calculated nbot for cell %d nbot %d check %d\n",mype,ic,nbot_global[ic],nbot_global_check[ic]);
            if (ntop_global[ic] != ntop_global_check[ic]) printf("%d: Error with gpu calculated ntop for cell %d ntop %d check %d\n",mype,ic,ntop_global[ic],ntop_global_check[ic]);
         }
      }

      mesh_local->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (do_cpu_calc && ! mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Need ghost cells for this routine
      if (do_cpu_calc) {
        state_global->apply_boundary_conditions(mesh_global);
        state_local->apply_boundary_conditions(mesh_local);
      }

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      if (do_cpu_calc) {
         state_global->calc_finite_difference(mesh_global, deltaT);
         state_local->calc_finite_difference_local(mesh_local, deltaT);
      }
      
      if (do_gpu_calc) {
         state_global->gpu_calc_finite_difference(command_queue, mesh_global, deltaT);
         state_local->gpu_calc_finite_difference_local(command_queue, mesh_local, deltaT);
      }
      
      if (do_comparison_calc) {
         // Need to compare dev_H to H, etc
         vector<real>H_save(ncells);
         vector<real>U_save(ncells);
         vector<real>V_save(ncells);
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real), &H_save[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real), &U_save[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real), &V_save[0], NULL);
         for (int ic = 0; ic < ncells; ic++){
            if (fabs(H[ic]-H_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d H & H_save %d %lf %lf \n",mype,n,ic,H[ic],H_save[ic]);
            if (fabs(U[ic]-U_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d U & U_save %d %lf %lf \n",mype,n,ic,U[ic],U_save[ic]);
            if (fabs(V[ic]-V_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d V & V_save %d %lf %lf \n",mype,n,ic,V[ic],V_save[ic]);
         }

         // And compare dev_H gathered to H_global, etc
         vector<real>H_save_global(ncells_global);
         vector<real>U_save_global(ncells_global);
         vector<real>V_save_global(ncells_global);
         MPI_Allgatherv(&H_save[0], nsizes[mype], MPI_C_REAL, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&U_save[0], nsizes[mype], MPI_C_REAL, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&V_save[0], nsizes[mype], MPI_C_REAL, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         if (mype == 0) {
            for (int ic = 0; ic < ncells_global; ic++){
               if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d H_global & H_save_global %d %lf %lf \n",mype,n,ic,H_global[ic],H_save_global[ic]);
               if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d U_global & U_save_global %d %lf %lf \n",mype,n,ic,U_global[ic],U_save_global[ic]);
               if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d V_global & V_save_global %d %lf %lf \n",mype,n,ic,V_global[ic],V_save_global[ic]);
            }    
         }    

         // And compare H gathered to H_global, etc
         MPI_Allgatherv(&H[0], nsizes[mype], MPI_C_REAL, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&U[0], nsizes[mype], MPI_C_REAL, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&V[0], nsizes[mype], MPI_C_REAL, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         for (int ic = 0; ic < ncells_global; ic++){
            if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d H_global & H_save_global %d %lf %lf \n",n,ic,H_global[ic],H_save_global[ic]);
            if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d U_global & U_save_global %d %lf %lf \n",n,ic,U_global[ic],U_save_global[ic]);
            if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d V_global & V_save_global %d %lf %lf \n",n,ic,V_global[ic],V_save_global[ic]);
         }    

         // Now the global dev_H_global to H_global, etc
         ezcl_enqueue_read_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_real), &H_save_global[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_real), &U_save_global[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V_global, CL_TRUE,  0, ncells_global*sizeof(cl_real), &V_save_global[0], NULL);
         if (mype == 0) {
            for (int ic = 0; ic < ncells_global; ic++){
               if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d H_global & H_save_global %d %lf %lf \n",mype,n,ic,H_global[ic],H_save_global[ic]);
               if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d U_global & U_save_global %d %lf %lf \n",mype,n,ic,U_global[ic],U_save_global[ic]);
               if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d V_global & V_save_global %d %lf %lf \n",mype,n,ic,V_global[ic],V_save_global[ic]);
            }
         }
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      if (do_cpu_calc) {
         state_global->remove_boundary_cells(mesh_global);
         state_local->remove_boundary_cells(mesh_local);
      }
      
      //  Check for NANs.
      for (int ic=0; ic<ncells; ic++) {
         if (isnan(H[ic]))
         {  printf("Got a NAN on cell %d cycle %d\n",ic,n);
            H[ic]=0.0;
            sleep(100);
#ifdef HAVE_OPENCL
            //  Release kernels and finalize the OpenCL elements.
            ezcl_finalize();
#endif
            L7_Terminate();
            exit(-1); }
      }  //  Complete NAN check.
      
      vector<int>      ioffset(block_size);
      vector<int>      ioffset_global(block_size_global);
      vector<int>      newcount_global(block_size_global);
      cl_mem dev_ioffset    = ezcl_malloc(NULL, &block_size, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);
      cl_mem dev_ioffset_global    = ezcl_malloc(NULL, &block_size_global, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);
      cl_mem dev_newcount_global    = ezcl_malloc(NULL, &block_size_global, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);

      if (do_cpu_calc) {
         mpot.resize(ncells);
         mpot_global.resize(ncells_global);
         state_global->calc_refine_potential(mesh_global, mpot_global, icount_global, jcount);
         state_local->calc_refine_potential_local(mesh_local, mpot, icount, jcount);
         nlft.clear();
         nrht.clear();
         nbot.clear();
         ntop.clear();
         nlft_global.clear();
         nrht_global.clear();
         nbot_global.clear();
         ntop_global.clear();
      }  //  Complete CPU calculation.
  
      size_t result_size = 1;
      cl_mem dev_result = NULL;
      cl_mem dev_result_global = NULL;
      if (do_gpu_calc) {
         dev_mpot     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
         dev_mpot_global     = ezcl_malloc(NULL, &ncells_global, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
         dev_result  = ezcl_malloc(NULL, &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_result_global  = ezcl_malloc(NULL, &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
 
         state_global->gpu_calc_refine_potential(command_queue, mesh_global, dev_mpot_global, dev_result_global, dev_ioffset_global, dev_newcount_global);
         state_local->gpu_calc_refine_potential_local(command_queue, mesh_local, dev_mpot, dev_result, dev_ioffset);

         ezcl_device_memory_remove(dev_nlft);
         ezcl_device_memory_remove(dev_nrht);
         ezcl_device_memory_remove(dev_nbot);
         ezcl_device_memory_remove(dev_ntop);
         ezcl_device_memory_remove(dev_nlft_global);
         ezcl_device_memory_remove(dev_nrht_global);
         ezcl_device_memory_remove(dev_nbot_global);
         ezcl_device_memory_remove(dev_ntop_global);
      }
      
      if (do_comparison_calc) {
         int icount_test;
         MPI_Allreduce(&icount, &icount_test, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if (icount_test != icount_global) {
            printf("%d: DEBUG -- icount is %d icount_test %d icount_global is %d\n",mype,icount,icount_test,icount_global);
         }

         // Need to compare dev_mpot to mpot
         vector<int>mpot_save(ncells);
         ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot_save[0], NULL);
         for (int ic = 0; ic < ncells; ic++){
            if (mpot[ic] != mpot_save[ic]) {
               printf("%d: DEBUG refine_potential 1 at cycle %d cell %d mpot & mpot_save %d %d \n",mype,n,ic,mpot[ic],mpot_save[ic]);
            }    
         }    

         // Compare dev_mpot to mpot_global
         vector<int>mpot_save_global(ncells_global);
         MPI_Allgatherv(&mpot_save[0], nsizes[mype], MPI_INT, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
         for (int ic = 0; ic < ncells_global; ic++){
            if (mpot_global[ic] != mpot_save_global[ic]) {
               if (mype == 0) printf("%d: DEBUG refine_potential 2 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,n,ic,mpot_global[ic],mpot_save_global[ic]);
            }    
         }    

         // Compare mpot to mpot_global
         MPI_Allgatherv(&mpot[0], nsizes[mype], MPI_INT, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
         for (int ic = 0; ic < ncells_global; ic++){
            if (mpot_global[ic] != mpot_save_global[ic]) {
               if (mype == 0) printf("%d: DEBUG refine_potential 3 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,n,ic,mpot_global[ic],mpot_save_global[ic]);
            }    
         }    

         // Compare dev_mpot_global to mpot_global
         ezcl_enqueue_read_buffer(command_queue, dev_mpot_global, CL_TRUE,  0, ncells_global*sizeof(cl_int), &mpot_save_global[0], NULL);
         for (int ic = 0; ic < ncells_global; ic++){
            if (mpot_global[ic] != mpot_save_global[ic]) {
               if (mype == 0) printf("%d: DEBUG refine_potential 4 at cycle %d mpot_global & mpot_save_global %d %d \n",mype,n,ic,mpot_global[ic],mpot_save_global[ic]);
            }    
         }    
      }

      // Sync up cpu array with gpu version to reduce differences due to minor numerical differences
      // otherwise cell count will diverge causing code problems and crashes
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_mpot_global, CL_TRUE,  0, ncells_global*sizeof(cl_int), &mpot_global[0], NULL);
      }
      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_mpot_global, CL_TRUE,  0, ncells_global*sizeof(cl_int), &mpot_global[0], NULL);
      }

      int mcount, mtotal;
      if (do_comparison_calc) {
         // This compares ioffset for each block in the calculation
         ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int),       &ioffset[0], NULL);
         mtotal = 0;
         for (int ig=0; ig<(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
            mcount = 0;
            for (int ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
                if (ic >= old_ncells) break;
                if (celltype[ic] == REAL_CELL) {
                   mcount += mpot[ic] ? 4 : 1;
                } else {
                   mcount += mpot[ic] ? 2 : 1;
                }
            }
            if (mcount != ioffset[ig]) printf("%d: DEBUG ig %d ioffset %d mcount %d\n",mype,ig,ioffset[ig],mcount);
            mtotal += mcount;
         }

         // For global This compares ioffset for each block in the calculation
         ezcl_enqueue_read_buffer(command_queue, dev_ioffset_global, CL_TRUE, 0, block_size_global*sizeof(cl_int),       &ioffset_global[0], NULL);
         mtotal = 0;
         int count = 0;
         for (int ig=0; ig<(old_ncells_global+TILE_SIZE-1)/TILE_SIZE; ig++){
            mcount = 0;
            for (int ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
                if (ic >= old_ncells_global) break;
                if (celltype_global[ic] == REAL_CELL) {
                   mcount += mpot_global[ic] ? 4 : 1;
                } else {
                   mcount += mpot_global[ic] ? 2 : 1;
                }
            }
            if (mcount != ioffset_global[ig]) {
               printf("DEBUG global ig %d ioffset %d mcount %d\n",ig,ioffset_global[ig],mcount);
               count++;
            }
            if (count > 10) exit(0);
            mtotal += mcount;
         }

      }
      if (do_gpu_sync) {
        for (int ig=0; ig<(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
           mcount = 0;
           for (int ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
               if (ic >= old_ncells) break;
               if (celltype[ic] == REAL_CELL) {
                  mcount += mpot[ic] ? 4 : 1;
               } else {
                  mcount += mpot[ic] ? 2 : 1;
               }
           }
           ioffset[ig] = mcount;
        }
        ezcl_enqueue_write_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int),       &ioffset[0], NULL);
        for (int ig=0; ig<(old_ncells_global+TILE_SIZE-1)/TILE_SIZE; ig++){
           mcount = 0;
           for (int ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
               if (ic >= old_ncells_global) break;
               if (celltype_global[ic] == REAL_CELL) {
                  mcount += mpot_global[ic] ? 4 : 1;
               } else {
                  mcount += mpot_global[ic] ? 2 : 1;
               }
           }
           ioffset_global[ig] = mcount;
         }
         ezcl_enqueue_write_buffer(command_queue, dev_ioffset_global, CL_TRUE, 0, block_size_global*sizeof(cl_int),       &ioffset_global[0], NULL);
      }

      if (do_cpu_calc) {
         new_ncells_global = old_ncells_global+mesh_global->rezone_count(mpot_global);
         new_ncells = old_ncells+mesh_local->rezone_count(mpot);
      }

      if (do_gpu_calc) {
         mesh_global->gpu_rezone_count(command_queue, block_size_global, local_work_size_global, dev_ioffset_global, dev_result_global);
         mesh_local->gpu_rezone_count(command_queue, block_size, local_work_size, dev_ioffset, dev_result);
      }

      if (do_comparison_calc) {
         int result;
         ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, 1*sizeof(cl_int),       &result, NULL);
         if (new_ncells != result) printf("%d: DEBUG new_ncells not correct %d %d\n",mype,new_ncells,result);
         new_ncells = result;
         //printf("Result is %d\n",result[0]);

         ezcl_enqueue_read_buffer(command_queue, dev_result_global, CL_TRUE, 0, 1*sizeof(cl_int),       &result, NULL);
         if (new_ncells_global != result) printf("%d: DEBUG new_ncells_global not correct %d %d\n",mype,new_ncells_global,result);
         new_ncells_global = result;
      }

      if (do_gpu_calc) {
         ezcl_device_memory_remove(dev_result);
         ezcl_device_memory_remove(dev_result_global);
      }

      if (do_cpu_calc) {
         int add_ncells_global = new_ncells_global - old_ncells_global;
         int add_ncells = new_ncells - old_ncells;
         //printf("%d: DEBUG add %d new %d old %d\n",mype,add_ncells,new_ncells,old_ncells);
         state_global->rezone_all(mesh_global, mpot_global, add_ncells_global);
         state_local->rezone_all_local(mesh_local, mpot, add_ncells);
         mpot_global.clear();
         mpot.clear();
      }

      //  Resize the mesh, inserting cells where refinement is necessary.
      if (do_gpu_calc) {
         state_global->gpu_rezone_all(command_queue, mesh_global, ncells_global, new_ncells_global, old_ncells_global, localStencil, dev_mpot_global, dev_ioffset_global);
         state_local->gpu_rezone_all_local(command_queue, mesh_local, old_ncells, new_ncells, old_ncells, localStencil, dev_mpot, dev_ioffset);
      }

      if (do_comparison_calc) {
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

         // Need to compare dev_H to H, etc
         vector<real> H_check(ncells);
         vector<real> U_check(ncells);
         vector<real> V_check(ncells);
         vector<int> level_check(ncells);
         vector<int> celltype_check(ncells);
         vector<int> i_check(ncells);
         vector<int> j_check(ncells);
         /// Set read buffers for data.
         ezcl_enqueue_read_buffer(command_queue, dev_H,        CL_FALSE, 0, ncells*sizeof(cl_real), &H_check[0],         NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U,        CL_FALSE, 0, ncells*sizeof(cl_real), &U_check[0],         NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V,        CL_FALSE, 0, ncells*sizeof(cl_real), &V_check[0],         NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int),  &level_check[0],     NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  &celltype_check[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  &i_check[0],         NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_TRUE,  0, ncells*sizeof(cl_int),  &j_check[0],         NULL);
         for (int ic = 0; ic < ncells; ic++){
            if (fabs(H[ic]-H_check[ic]) > STATE_EPS) printf("%d: DEBUG rezone 1 cell %d H %lf H_check %lf\n",mype, ic,H[ic],H_check[ic]);
            if (fabs(U[ic]-U_check[ic]) > STATE_EPS) printf("%d: DEBUG rezone 1 cell %d U %lf U_check %lf\n",mype, ic,U[ic],U_check[ic]);
            if (fabs(V[ic]-V_check[ic]) > STATE_EPS) printf("%d: DEBUG rezone 1 cell %d V %lf V_check %lf\n",mype, ic,V[ic],V_check[ic]);
            if (level[ic] != level_check[ic] )       printf("%d: DEBUG rezone 1 cell %d level %d level_check %d\n",mype, ic, level[ic], level_check[ic]);
            if (celltype[ic] != celltype_check[ic] ) printf("%d: DEBUG rezone 1 cell %d celltype %d celltype_check %d\n",mype, ic, celltype[ic], celltype_check[ic]);
            if (i[ic] != i_check[ic] )               printf("%d: DEBUG rezone 1 cell %d i %d i_check %d\n",mype, ic, i[ic], i_check[ic]);
            if (j[ic] != j_check[ic] )               printf("%d: DEBUG rezone 1 cell %d j %d j_check %d\n",mype, ic, j[ic], j_check[ic]);
         }
      
         // And compare dev_H gathered to H_global, etc
         vector<real>H_check_global(ncells_global);
         vector<real>U_check_global(ncells_global);
         vector<real>V_check_global(ncells_global);
         vector<int>celltype_check_global(ncells_global);
         vector<int>i_check_global(ncells_global);
         vector<int>j_check_global(ncells_global);
         vector<int>level_check_global(ncells_global);
         MPI_Allgatherv(&H_check[0],        nsizes[mype], MPI_C_REAL, &H_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&U_check[0],        nsizes[mype], MPI_C_REAL, &U_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&V_check[0],        nsizes[mype], MPI_C_REAL, &V_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&celltype_check[0], nsizes[mype], MPI_INT,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&i_check[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&j_check[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&level_check[0],    nsizes[mype], MPI_INT,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         for (int ic = 0; ic < ncells_global; ic++){
            if (fabs(H_global[ic]-H_check_global[ic]) > STATE_EPS) printf("%d: DEBUG rezone 2 cell %d H_global %lf H_check_global %lf \n",mype,ic,H_global[ic],H_check_global[ic]);
            if (fabs(U_global[ic]-U_check_global[ic]) > STATE_EPS) printf("%d: DEBUG rezone 2 cell %d U_global %lf U_check_global %lf \n",mype,ic,U_global[ic],U_check_global[ic]);
            if (fabs(V_global[ic]-V_check_global[ic]) > STATE_EPS) printf("%d: DEBUG rezone 2 cell %d V_global %lf V_check_global %lf \n",mype,ic,V_global[ic],V_check_global[ic]);
            if (level_global[ic] != level_check_global[ic] )       printf("%d: DEBUG rezone 2 cell %d level_global %d level_check_global %d\n",mype, ic, level_global[ic], level_check_global[ic]);
            if (celltype_global[ic] != celltype_check_global[ic] ) printf("%d: DEBUG rezone 2 cell %d celltype_global %d celltype_check_global %d\n",mype, ic, celltype_global[ic], celltype_check_global[ic]);
            if (i_global[ic] != i_check_global[ic] )               printf("%d: DEBUG rezone 2 cell %d i_global %d i_check_global %d\n",mype, ic, i_global[ic], i_check_global[ic]);
            if (j_global[ic] != j_check_global[ic] )               printf("%d: DEBUG rezone 2 cell %d j_global %d j_check_global %d\n",mype, ic, j_global[ic], j_check_global[ic]);
         }

         // And compare H gathered to H_global, etc
         MPI_Allgatherv(&H[0],        nsizes[mype], MPI_C_REAL, &H_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&U[0],        nsizes[mype], MPI_C_REAL, &U_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&V[0],        nsizes[mype], MPI_C_REAL, &V_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&celltype[0], nsizes[mype], MPI_INT,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&i[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&j[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&level[0],    nsizes[mype], MPI_INT,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         for (int ic = 0; ic < ncells_global; ic++){
            if (fabs(H_global[ic]-H_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 3 at cycle %d H_global & H_check_global %d %lf %lf \n",n,ic,H_global[ic],H_check_global[ic]);
            if (fabs(U_global[ic]-U_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 3 at cycle %d U_global & U_check_global %d %lf %lf \n",n,ic,U_global[ic],U_check_global[ic]);
            if (fabs(V_global[ic]-V_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 3 at cycle %d V_global & V_check_global %d %lf %lf \n",n,ic,V_global[ic],V_check_global[ic]);
            if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 3 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",n,ic,celltype_global[ic],celltype_check_global[ic]);
            if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 3 at cycle %d i_global & i_check_global %d %d  %d  \n",n,ic,i_global[ic],i_check_global[ic]);
            if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 3 at cycle %d j_global & j_check_global %d %d  %d  \n",n,ic,j_global[ic],j_check_global[ic]);
            if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 3 at cycle %d level_global & level_check_global %d %d  %d  \n",n,ic,level_global[ic],level_check_global[ic]);
         }

         // Now the global dev_H_global to H_global, etc
         ezcl_enqueue_read_buffer(command_queue, dev_H_global,        CL_FALSE, 0, ncells_global*sizeof(cl_real), &H_check_global[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U_global,        CL_FALSE, 0, ncells_global*sizeof(cl_real), &U_check_global[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V_global,        CL_FALSE, 0, ncells_global*sizeof(cl_real), &V_check_global[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_celltype_global, CL_FALSE, 0, ncells_global*sizeof(cl_int),  &celltype_check_global[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_i_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  &i_check_global[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_j_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  &j_check_global[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level_global,    CL_TRUE,  0, ncells_global*sizeof(cl_int),  &level_check_global[0],    NULL);
         for (int ic = 0; ic < ncells_global; ic++){
            if (fabs(H_global[ic]-H_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 4 at cycle %d H_global & H_check_global %d %lf %lf \n",n,ic,H_global[ic],H_check_global[ic]);
            if (fabs(U_global[ic]-U_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 4 at cycle %d U_global & U_check_global %d %lf %lf \n",n,ic,U_global[ic],U_check_global[ic]);
            if (fabs(V_global[ic]-V_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 4 at cycle %d V_global & V_check_global %d %lf %lf \n",n,ic,V_global[ic],V_check_global[ic]);
            if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 4 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",n,ic,celltype_global[ic],celltype_check_global[ic]);
            if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 4 at cycle %d i_global & i_check_global %d %d  %d  \n",n,ic,i_global[ic],i_check_global[ic]);
            if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 4 at cycle %d j_global & j_check_global %d %d  %d  \n",n,ic,j_global[ic],j_check_global[ic]);
            if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 4 at cycle %d level_global & level_check_global %d %d  %d  \n",n,ic,level_global[ic],level_check_global[ic]);
         }

      } // do_comparison_calc

      //if (n % outputInterval == 0) {
      //   vector<real> H_save(ncells);
      //   ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_real), &H_save[0], &start_read_event);
      //   gpu_time_read             += ezcl_timer_calc(&start_read_event,       &start_read_event);
      //}

#ifdef HAVE_OPENGL
      set_cell_data(&H_global[0]);
      set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
#endif

#ifdef XXX // not rewritten yet
      if (do_gpu_calc)
      {
          if (! mesh->have_boundary) {
              /*
              __kernel void count_BCs(
                               const int    isize,      // 0   
                      __global const int   *i,         // 1
                      __global const int   *j,         // 2
                      __global const int   *level,     // 3
                      __global const int   *lev_ibeg,  // 4
                      __global const int   *lev_iend,  // 5
                      __global const int   *lev_jbeg,  // 6
                      __global const int   *lev_jend,  // 7
                      __global       int   *scratch,   // 8
                      __local        int   *tile)      // 9
              */
            size_t shared_spd_sum_int = local_work_size * sizeof(cl_int);
            ezcl_set_kernel_arg(kernel_count_BCs, 0, sizeof(cl_int), (void *)&ncells);
            ezcl_set_kernel_arg(kernel_count_BCs, 1, sizeof(cl_mem), (void *)&dev_i);
            ezcl_set_kernel_arg(kernel_count_BCs, 2, sizeof(cl_mem), (void *)&dev_j);
            ezcl_set_kernel_arg(kernel_count_BCs, 3, sizeof(cl_mem), (void *)&dev_level);
            ezcl_set_kernel_arg(kernel_count_BCs, 4, sizeof(cl_mem), (void *)&dev_levibeg);
            ezcl_set_kernel_arg(kernel_count_BCs, 5, sizeof(cl_mem), (void *)&dev_leviend);
            ezcl_set_kernel_arg(kernel_count_BCs, 6, sizeof(cl_mem), (void *)&dev_levjbeg);
            ezcl_set_kernel_arg(kernel_count_BCs, 7, sizeof(cl_mem), (void *)&dev_levjend);
            ezcl_set_kernel_arg(kernel_count_BCs, 8, sizeof(cl_mem), (void *)&dev_ioffset);
            ezcl_set_kernel_arg(kernel_count_BCs, 9, shared_spd_sum_int, 0);

            ezcl_set_kernel_arg(kernel_reduce_sum_int_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
            ezcl_set_kernel_arg(kernel_reduce_sum_int_stage2of2, 1, sizeof(cl_mem), (void *)&dev_ioffset);
            ezcl_set_kernel_arg(kernel_reduce_sum_int_stage2of2, 2, shared_spd_sum_int, 0);

            ezcl_enqueue_ndrange_kernel(command_queue, kernel_count_BCs, 1, NULL, &global_work_size, &local_work_size, &count_BCs_stage1_event);

            if (block_size > 1) {
               ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_int_stage2of2, 1, NULL, &local_work_size, &local_work_size, &count_BCs_stage2_event);
            }

            ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, 1*sizeof(cl_int), &ioffset[0], NULL);
            bcount = ioffset[0];
            //printf("DEBUG -- bcount is %d\n",bcount);
            //gpu_time_read += ezcl_timer_calc(&start_read_event, &start_read_event);
         }
      }
#endif

      if (do_gpu_calc) {
         if (ncells != old_ncells){
            H.resize(ncells);
            U.resize(ncells);
            V.resize(ncells);
            level.resize(ncells);
            i.resize(ncells);
            j.resize(ncells);
            celltype.resize(ncells);
         }
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

      double H_sum = -1.0;

      if (do_comparison_calc) {
         double cpu_mass_sum = state_global->mass_sum(mesh_global, enhanced_precision_sum);
         double cpu_mass_sum_local = state_local->mass_sum_local(mesh_local, enhanced_precision_sum);
         double gpu_mass_sum = state_global->gpu_mass_sum(command_queue, mesh_global, enhanced_precision_sum);
         H_sum = state_local->gpu_mass_sum_local(command_queue, mesh_local, enhanced_precision_sum);
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
      mesh_local->proc.resize(ncells);
      if (icount) {
         vector<int> index(ncells);
         vector<int> index_global(ncells_global);
         mesh_global->partition_cells(numpe, mesh_global->proc, index_global, cycle_reorder);
         mesh_local->partition_cells(numpe, mesh_local->proc, index, cycle_reorder);
         //state->state_reorder(index);
         if (do_gpu_sync) {
            ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_level, CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&level[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_celltype_global, CL_FALSE, 0, ncells_global*sizeof(cl_int),  (void *)&celltype_global[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_i_global,     CL_FALSE, 0, ncells_global*sizeof(cl_int),  (void *)&i_global[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_j_global,     CL_FALSE, 0, ncells_global*sizeof(cl_int),  (void *)&j_global[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_level_global, CL_TRUE,  0, ncells_global*sizeof(cl_int),  (void *)&level_global[0],  NULL);
         }
      }

      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&H_global[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_real),  (void *)&U_global[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V_global, CL_TRUE,  0, ncells_global*sizeof(cl_real),  (void *)&V_global[0],  NULL);
      }

      
      if (n % outputInterval == 0) {
         if (H_sum < 0) {
            H_sum = state_global->mass_sum(mesh_global, enhanced_precision_sum);
         }
         if (mype == 0){
            printf("Iteration %d timestep %lf Sim Time %lf cells %d Mass Sum %14.12lg Mass Change %14.12lg\n",
               n, deltaT, simTime, ncells, H_sum, H_sum - H_sum_initial);
         }
#ifdef HAVE_OPENGL
         set_mysize(ncells_global);
         set_viewmode(view_mode);
         set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
         set_cell_data(&H_global[0]);
         set_cell_proc(&mesh_global->proc[0]);
         set_circle_radius(circle_radius);
         DrawGLScene();
#endif
      }
      ++n;
      simTime += deltaT;
      
      //  Calculate timings for this time step.
#ifdef XXX
      if (do_gpu_calc) {
         if (! mesh->have_boundary) {
            gpu_time_count_BCs        += ezcl_timer_calc(&count_BCs_stage1_event, &count_BCs_stage1_event);
            gpu_time_count_BCs        += ezcl_timer_calc(&count_BCs_stage2_event, &count_BCs_stage2_event);
         }
      }
#endif
   }  //  Complete output interval.
}

