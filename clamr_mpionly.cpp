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

#define DO_COMPARISON

//TODO:  command-line option for OpenGL?
#ifdef DO_COMPARISON
int do_comparison_calc = 1;
#else
int do_comparison_calc = 0;
#endif

int do_cpu_calc = 1;
int do_gpu_calc = 0;

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
Mesh       *mesh;           //  Object containing mesh information; init in grid.cpp::main().
State      *state;    //  Object containing state information corresponding to mesh; init in grid.cpp::main().

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

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   L7_Init(&mype, &numpe, &argc, argv);
   parseInput(argc, argv);

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
   
   parallel_in = 1;
   mesh = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   state = new State(mesh->ncells, context);

   size_t &ncells = mesh->ncells;
   int &noffset = mesh->noffset;

   vector<int>   &nsizes     = mesh_global->nsizes;
   vector<int>   &ndispl     = mesh_global->ndispl;

   vector<real>  &H_global = state_global->H;
   vector<real>  &U_global = state_global->U;
   vector<real>  &V_global = state_global->V;

   vector<real>  &x_global  = mesh_global->x;
   vector<real>  &dx_global = mesh_global->dx;
   vector<real>  &y_global  = mesh_global->y;
   vector<real>  &dy_global = mesh_global->dy;

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
   MPI_Allgather(&ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
   ndispl[0]=0;
   for (int ip=1; ip<numpe; ip++){
      ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
   }
   noffset=0;
   for (int ip=0; ip<mype; ip++){
     noffset += nsizes[ip];
   }

   // Distribute level, celltype, H, U, V

   celltype.resize(ncells);
   level.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);

   MPI_Scatterv(&celltype_global[0], &nsizes[0], &ndispl[0], MPI_INT, &celltype[0], nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&level_global[0],    &nsizes[0], &ndispl[0], MPI_INT, &level[0],    nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&i_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &i[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&j_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &j[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);

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

#ifdef HAVE_OPENGL
   set_cell_data(&H_global[0]);
   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_mysize(ncells_global);
   set_viewmode(view_mode);
   set_window(mesh_global->xmin, mesh_global->xmax, mesh_global->ymin, mesh_global->ymax);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water", mype);
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
   int icount, jcount;
   int icount_global, jcount_global;

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
   vector<int>   &nlft_global     = mesh_global->nlft;
   vector<int>   &nrht_global     = mesh_global->nrht;
   vector<int>   &nbot_global     = mesh_global->nbot;
   vector<int>   &ntop_global     = mesh_global->ntop;

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;
   vector<int>   &nlft     = mesh->nlft;
   vector<int>   &nrht     = mesh->nrht;
   vector<int>   &nbot     = mesh->nbot;
   vector<int>   &ntop     = mesh->ntop;

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
      
      state->output_timing_info(mesh, do_cpu_calc, do_gpu_calc, gpu_time_count_BCs_parallel, elapsed_time);
      if (do_comparison_calc) {
         state_global->output_timing_info(mesh_global, do_cpu_calc, do_gpu_calc, gpu_time_count_BCs, elapsed_time);
      }

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
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

   //  Main loop.
   for (int iburst = 0; iburst < outputInterval; iburst++) {
      if (n > niter) break;

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;
      old_ncells_global = ncells_global;

      //  Calculate the real time step for the current discrete time step.
      double deltaT = state->set_timestep(mesh, g, sigma);

      //  Compare time step values and pass deltaT in to the kernel.
      if (do_comparison_calc) {
         double deltaT_cpu_global = state_global->set_timestep(mesh_global, g, sigma);

         if (fabs(deltaT - deltaT_cpu_global) > .000001) {
            printf("Error with deltaT calc --- cpu_local %lf cpu_global %lf\n",
               deltaT, deltaT_cpu_global);
         }
      }

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

      if (do_comparison_calc) {
         mesh_global->calc_neighbors();

         // Checking CPU parallel to CPU global
         vector<int> Test(ncells_ghost);
         for(uint ic=0; ic<ncells; ic++){
            Test[ic] = mype*1000 +ic;
         }
         L7_Update(&Test[0], L7_INT, cell_handle);

         vector<int> Test_global(ncells_global);
         MPI_Allgatherv(&Test[0], nsizes[mype], MPI_INT, &Test_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         vector<int> Test_check(ncells);
         vector<int> Test_check_global(ncells_global);

         // ==================== check left value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nlft[ic]];
            //if (mype == 1 && ic==0) printf("%d: nlft check for ic 0 is %d\n",mype,nlft[0]);
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            //if (Test_global[nlft_global[ic]] != Test_check_global[ic]) {
               //if (mype == 0) printf("%d: Error with nlft for cell %d -- nlft %d global %d check %d\n",mype,ic,nlft_global[ic],Test_global[nlft_global[ic]],Test_check_global[ic]);
            //}
         }
         
         // ==================== check left left value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nlft[nlft[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[nlft_global[nlft_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with nlft nlft for cell %5d -- nlftg %5d nlftg nlftg %5d global %5d\n",
                  mype,ic,nlft_global[ic],nlft_global[nlft_global[ic]],Test_global[nlft_global[nlft_global[ic]]]);
               printf("%d:                         check %5d -- nlftl %5d nlftl nlftl %5d check  %5d\n",
                  mype,ic,nlft[ic],nlft[nlft[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check right value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nrht[ic]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[nrht_global[ic]] != Test_check_global[ic]) {
               if (mype == 0) printf("%d: Error with nrht for cell %d -- %d %d\n",mype,ic,Test_global[nrht_global[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check right right value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nrht[nrht[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[nrht_global[nrht_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with nrht nrht for cell %5d -- nrhtg %5d nrhtg nrhtg %5d global %5d\n",
                  mype,ic,nrht_global[ic],nrht_global[nrht_global[ic]],Test_global[nrht_global[nrht_global[ic]]]);
               printf("%d:                         check %5d -- nrhtl %5d nrhtl nrhtl %5d check  %5d\n",
                  mype,ic,nrht[ic],nrht[nrht[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check bottom value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nbot[ic]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[nbot_global[ic]] != Test_check_global[ic]) {
               if (mype == 0) printf("%d: Error with nbot for cell %d -- %d %d\n",mype,ic,Test_global[nbot_global[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check bottom bottom value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[nbot[nbot[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[nbot_global[nbot_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with nbot nbot for cell %5d -- nbotg %5d nbotg nbotg %5d global %5d\n",
                  mype,ic,nbot_global[ic],nbot_global[nbot_global[ic]],Test_global[nbot_global[nbot_global[ic]]]);
               printf("%d:                         check %5d -- nbotl %5d nbotl nbotl %5d check  %5d\n",
                  mype,ic,nbot[ic],nbot[nbot[ic]],Test_check_global[ic]);
            }
         }
         
         // ==================== check top value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[ntop[ic]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[ntop_global[ic]] != Test_check_global[ic]) {
               if (mype == 0) printf("%d: Error with ntop for cell %d -- %d %d\n",mype,ic,Test_global[ntop_global[ic]],Test_check_global[ic]);
            }
         }

         // ==================== check top top value ====================
         for (uint ic=0; ic<ncells; ic++){
            Test_check[ic] = Test[ntop[ntop[ic]]];
         }

         MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

         for (uint ic=0; ic<ncells_global; ic++){
            if (Test_global[ntop_global[ntop_global[ic]]] != Test_check_global[ic]) {
               printf("%d: Error with ntop ntop for cell %5d -- ntopg %5d ntopg ntopg %5d global %5d\n",
                  mype,ic,ntop_global[ic],ntop_global[ntop_global[ic]],Test_global[ntop_global[ntop_global[ic]]]);
               printf("%d:                         check %5d -- ntopl %5d ntopl ntopl %5d check  %5d\n",
                  mype,ic,ntop[ic],ntop[ntop[ic]],Test_check_global[ic]);
            }
         }
         
      }

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Need ghost cells for this routine
      state->apply_boundary_conditions(mesh);

      if (do_comparison_calc) {
        state_global->apply_boundary_conditions(mesh_global);
      }

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->calc_finite_difference_local(mesh, deltaT);

      if (do_comparison_calc) {
         state_global->calc_finite_difference(mesh_global, deltaT);

         // Compare H gathered to H_global, etc
         vector<real>H_save_global(ncells_global);
         vector<real>U_save_global(ncells_global);
         vector<real>V_save_global(ncells_global);
         MPI_Allgatherv(&H[0], nsizes[mype], MPI_C_REAL, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&U[0], nsizes[mype], MPI_C_REAL, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&V[0], nsizes[mype], MPI_C_REAL, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         for (uint ic = 0; ic < ncells_global; ic++){
            if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d H_global & H_save_global %d %lf %lf \n",n,ic,H_global[ic],H_save_global[ic]);
            if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d U_global & U_save_global %d %lf %lf \n",n,ic,U_global[ic],U_save_global[ic]);
            if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d V_global & V_save_global %d %lf %lf \n",n,ic,V_global[ic],V_save_global[ic]);
         }    
      }

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells(mesh);

      if (do_comparison_calc) {
         state_global->remove_boundary_cells(mesh_global);
      }
      
      //  Check for NANs.
      for (uint ic=0; ic<ncells; ic++) {
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
      
      mpot.resize(ncells_ghost);
      state->calc_refine_potential_local(mesh, mpot, icount, jcount);
      nlft.clear();
      nrht.clear();
      nbot.clear();
      ntop.clear();
  
      if (do_comparison_calc) {
         mpot_global.resize(ncells_global);
         state_global->calc_refine_potential(mesh_global, mpot_global, icount_global, jcount_global);
         nlft_global.clear();
         nrht_global.clear();
         nbot_global.clear();
         ntop_global.clear();

         int icount_test;
         MPI_Allreduce(&icount, &icount_test, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if (icount_test != icount_global) {
            printf("%d: DEBUG -- icount is %d icount_test %d icount_global is %d\n", mype,icount,icount_test,icount_global);
         }

         // Compare mpot to mpot_global
         vector<int>mpot_save_global(ncells_global);
         MPI_Allgatherv(&mpot[0], nsizes[mype], MPI_INT, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
         for (uint ic = 0; ic < ncells_global; ic++){
            if (mpot_global[ic] != mpot_save_global[ic]) {
               if (mype == 0) printf("%d: DEBUG refine_potential 3 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,n,ic,mpot_global[ic],mpot_save_global[ic]);
            }    
         }    
      }

      new_ncells = old_ncells+mesh->rezone_count(mpot);

      if (do_comparison_calc) {
         new_ncells_global = old_ncells_global+mesh_global->rezone_count(mpot_global);
      }

      int add_ncells = new_ncells - old_ncells;
      state->rezone_all_local(mesh, mpot, add_ncells);
      mpot.clear();
      MPI_Allgather(&ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }
      noffset=0;
      for (int ip=0; ip<mype; ip++){
         noffset += nsizes[ip];
      }

      if (do_comparison_calc) {
         int add_ncells_global = new_ncells_global - old_ncells_global;
         //printf("%d: DEBUG add %d new %d old %d\n",mype,add_ncells,new_ncells,old_ncells);
         state_global->rezone_all(mesh_global, mpot_global, add_ncells_global);
         mpot_global.clear();

         //printf("%d: DEBUG ncells is %d new_ncells %d old_ncells %d ncells_global %d\n",mype, ncells, new_ncells, old_ncells, ncells_global);

         // And compare H gathered to H_global, etc
         vector<real>H_check_global(ncells_global);
         vector<real>U_check_global(ncells_global);
         vector<real>V_check_global(ncells_global);
         vector<int>celltype_check_global(ncells_global);
         vector<int>i_check_global(ncells_global);
         vector<int>j_check_global(ncells_global);
         vector<int>level_check_global(ncells_global);
         MPI_Allgatherv(&H[0],        nsizes[mype], MPI_C_REAL, &H_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&U[0],        nsizes[mype], MPI_C_REAL, &U_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&V[0],        nsizes[mype], MPI_C_REAL, &V_check_global[0],        &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
         MPI_Allgatherv(&celltype[0], nsizes[mype], MPI_INT,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&i[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&j[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         MPI_Allgatherv(&level[0],    nsizes[mype], MPI_INT,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
         for (uint ic = 0; ic < ncells_global; ic++){
            if (fabs(H_global[ic]-H_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 3 at cycle %d H_global & H_check_global %d %lf %lf \n",n,ic,H_global[ic],H_check_global[ic]);
            if (fabs(U_global[ic]-U_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 3 at cycle %d U_global & U_check_global %d %lf %lf \n",n,ic,U_global[ic],U_check_global[ic]);
            if (fabs(V_global[ic]-V_check_global[ic]) > STATE_EPS) printf("DEBUG rezone 3 at cycle %d V_global & V_check_global %d %lf %lf \n",n,ic,V_global[ic],V_check_global[ic]);
            if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 3 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",n,ic,celltype_global[ic],celltype_check_global[ic]);
            if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 3 at cycle %d i_global & i_check_global %d %d  %d  \n",n,ic,i_global[ic],i_check_global[ic]);
            if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 3 at cycle %d j_global & j_check_global %d %d  %d  \n",n,ic,j_global[ic],j_check_global[ic]);
            if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 3 at cycle %d level_global & level_check_global %d %d  %d  \n",n,ic,level_global[ic],level_check_global[ic]);
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

      double H_sum = -1.0;

      if (do_comparison_calc) {
         H_sum = state->mass_sum_local(mesh, enhanced_precision_sum);

         double H_sum_global = state_global->mass_sum(mesh_global, enhanced_precision_sum);

         if (fabs(H_sum - H_sum_global) > CONSERVATION_EPS) {
            printf("Error with mass sum calculation -- mass_sum %lf mass_sum_global %lf\n",
                    H_sum, H_sum_global);
         }
      }

      mesh->proc.resize(ncells);
      if (icount) {
         vector<int> index(ncells);
         mesh->partition_cells(numpe, mesh->proc, index, cycle_reorder);
      }

      if (do_comparison_calc) {
         mesh_global->proc.resize(ncells_global);

         if (icount) {
            vector<int> index_global(ncells_global);
            mesh_global->partition_cells(numpe, mesh_global->proc, index_global, cycle_reorder);
            //state->state_reorder(index);
         }
      }

      if (n % outputInterval == 0) {
         if (H_sum < 0) {
            H_sum = state->mass_sum_local(mesh, enhanced_precision_sum);
         }
         if (mype == 0){
            printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %14.12lg\n",
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
      }  //  Complete output interval.

      ++n;
      simTime += deltaT;
      
   }
}

