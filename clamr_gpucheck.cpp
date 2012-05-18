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

int do_comparison_calc = 1;
int do_cpu_calc = 1;
int do_gpu_calc = 1;

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

double circle_radius=-1.0;

int view_mode = 0;

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
Mesh       *mesh;           //  Object containing mesh information; init in grid.cpp::main().
State      *state;          //  Object containing state information corresponding to mesh; init in grid.cpp::main().

//  Set up timing information.
struct timeval tstart, tstop, tresult;
struct timeval tstart_cpu;
cl_event start_write_event, end_write_event,
         start_read_event,  end_read_event;
double   cpu_time_partition          = 0;

#ifdef HAVE_OPENCL
cl_context          context                 = NULL;
cl_command_queue    command_queue           = NULL;
#endif

int main(int argc, char **argv) {
#ifdef HAVE_OPENCL
   int ierr;
#endif

    // Needed for code to compile correctly on the Mac
   int mype=0;
   int numpe=-1;

   //  Process command-line arguments, if any.
   parseInput(argc, argv);
   
   numpe = 16;

#ifdef HAVE_OPENCL
   ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU, &context, &command_queue, 0);
   if (ierr == EZCL_NODEVICE) {
      ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU, &context, &command_queue, 0);
   }
   if (ierr != EZCL_SUCCESS) {
      printf("No opencl device available -- aborting\n");
      exit(-1);
   }
#endif

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   
   mesh  = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   mesh->init(nx, ny, circ_radius, context, initial_order, special_case, do_gpu_calc);
   size_t &ncells = mesh->ncells;
   state = new State(ncells, context);
   state->init(ncells, context, do_gpu_calc);
   mesh->proc.resize(ncells);
   mesh->calc_distribution(numpe, mesh->proc);
   state->fill_circle(mesh, circ_radius, 100.0, 5.0);
   
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;

   cl_mem &dev_celltype_new = mesh->dev_celltype_new;
   cl_mem &dev_i_new        = mesh->dev_i_new;
   cl_mem &dev_j_new        = mesh->dev_j_new;
   cl_mem &dev_level_new    = mesh->dev_level_new;

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
   state->dev_deltaT   = ezcl_malloc(NULL, &one,    sizeof(cl_real),  CL_MEM_READ_WRITE, 0);

   dev_celltype = ezcl_malloc(NULL, &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_i        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_j        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);
   dev_level    = ezcl_malloc(NULL, &ncells, sizeof(cl_int),   CL_MEM_READ_ONLY, 0);

   ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0], &start_write_event);
   ezcl_enqueue_write_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],        NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&level[0],    NULL            );
   ezcl_enqueue_write_buffer(command_queue, dev_H,        CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],       NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_U,        CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],       NULL              );
   ezcl_enqueue_write_buffer(command_queue, dev_V,        CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],       &end_write_event  );
   state->gpu_time_write += ezcl_timer_calc(&start_write_event, &end_write_event);

   dev_celltype_new = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_i_new        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_j_new        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);
   dev_level_new    = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_WRITE_ONLY, 0);

#ifdef HAVE_GRAPHICS
   set_mysize(ncells);
   set_viewmode(view_mode);
   set_window(mesh->xmin, mesh->xmax, mesh->ymin, mesh->ymax);
   set_outline((int)outline);
   init_display(&argc, argv, "Shallow Water", mype);
   set_idle_function(&do_calc);
   start_main_loop();
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

   if (cycle_reorder == ZORDER || cycle_reorder == HILBERT_SORT) {
      do_comparison_calc = 1;
      do_sync = 0;
      do_gpu_sync = 1;
   }
   
#ifdef HAVE_OPENCL
   cl_mem dev_mpot  = NULL;
#endif

   //  Initialize state variables for GPU calculation.
   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;
   vector<int>   &nlft     = mesh->nlft;
   vector<int>   &nrht     = mesh->nrht;
   vector<int>   &nbot     = mesh->nbot;
   vector<int>   &ntop     = mesh->ntop;

   size_t ncells    = mesh->ncells;

   vector<real>  &x        = mesh->x;
   vector<real>  &dx       = mesh->dx;
   vector<real>  &y        = mesh->y;
   vector<real>  &dy       = mesh->dy;

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
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;

   //  Kahan-type enhanced precision sum implementation.
   if (n < 0)
   {
      double H_sum = state->mass_sum(mesh, enhanced_precision_sum);
      printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
      H_sum_initial = H_sum;
      n++;
      
      //  Set up grid.
#ifdef GRAPHICS_OUTPUT
      mesh->write_grid(n);
#endif
#ifdef HAVE_GRAPHICS
      set_mysize(ncells);
      set_viewmode(view_mode);
      set_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
      set_cell_data(&H[0]);
      set_cell_proc(&mesh->proc[0]);
      set_circle_radius(circle_radius);
      draw_scene();
      if (verbose) sleep(5);
#endif
      //  Set flag to show mesh results rather than domain decomposition.
      view_mode = 1;
   
      //  Clear superposition of circle on grid output.
      if (n > 2) circle_radius = -1.0;
   
      gettimeofday(&tstart, NULL);
      return;
   }
   
   vector<int>     mpot;
   
   size_t old_ncells = ncells;
   size_t new_ncells = 0;

   //  Main loop.
   int output_flag = 0;
   while (! output_flag)
   {
      if (n > niter) break;

      // To reduce drift in solution
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);
      }

      size_t local_work_size  = MIN(ncells, TILE_SIZE);
      size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

      //size_t block_size = (ncells + TILE_SIZE - 1) / TILE_SIZE; //  For on-device global reduction kernel.
      size_t block_size     = global_work_size/local_work_size;

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
      
      double deltaT = (do_gpu_calc) ? deltaT_gpu : deltaT_cpu;

      if (do_cpu_calc) {
         mesh->calc_neighbors();
      }

      if (do_gpu_calc) {
         mesh->gpu_calc_neighbors(command_queue);
      }

      if (do_comparison_calc) {
         vector<int> nlft_check(ncells);
         vector<int> nrht_check(ncells);
         vector<int> nbot_check(ncells);
         vector<int> ntop_check(ncells);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft,     CL_FALSE, 0, ncells*sizeof(cl_int),  &nlft_check[0],      NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht,     CL_FALSE, 0, ncells*sizeof(cl_int),  &nrht_check[0],      NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot,     CL_FALSE, 0, ncells*sizeof(cl_int),  &nbot_check[0],      NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop,     CL_TRUE,  0, ncells*sizeof(cl_int),  &ntop_check[0],      NULL);

         for (uint ic=0; ic<ncells; ic++){
            if (nlft[ic] != nlft_check[ic]) printf("DEBUG -- nlft: ic %d nlft %d nlft_check %d\n",ic, nlft[ic], nlft_check[ic]);
            if (nrht[ic] != nrht_check[ic]) printf("DEBUG -- nrht: ic %d nrht %d nrht_check %d\n",ic, nrht[ic], nrht_check[ic]);
            if (nbot[ic] != nbot_check[ic]) printf("DEBUG -- nbot: ic %d nbot %d nbot_check %d\n",ic, nbot[ic], nbot_check[ic]);
            if (ntop[ic] != ntop_check[ic]) printf("DEBUG -- ntop: ic %d ntop %d ntop_check %d\n",ic, ntop[ic], ntop_check[ic]);
         }
      }
      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (do_cpu_calc && ! mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Need ghost cells for this routine
      if (do_cpu_calc) {
        state->apply_boundary_conditions(mesh);
      }
      
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
         vector<real>H_save(ncells);
         vector<real>U_save(ncells);
         vector<real>V_save(ncells);
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real), &H_save[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real), &U_save[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real), &V_save[0], NULL);
         for (uint ic = 0; ic < ncells; ic++){
            if (fabs(H[ic]-H_save[ic]) > STATE_EPS) printf("DEBUG finite_difference at cycle %d H & H_save %d %lf %lf \n",n,ic,H[ic],H_save[ic]);
            if (fabs(U[ic]-U_save[ic]) > STATE_EPS) printf("DEBUG finite_difference at cycle %d U & U_save %d %lf %lf \n",n,ic,U[ic],U_save[ic]);
            if (fabs(V[ic]-V_save[ic]) > STATE_EPS) printf("DEBUG finite_difference at cycle %d V & V_save %d %lf %lf \n",n,ic,V[ic],V_save[ic]);
         }
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
         printf("Cycle is %d\n",n);
         sleep(2);
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
            exit(-1); }
      }  //  Complete NAN check.
      
      vector<int>      ioffset(block_size);
      cl_mem dev_ioffset    = ezcl_malloc(NULL, &block_size, sizeof(cl_int),   CL_MEM_READ_WRITE, 0);

      if (do_cpu_calc) {
         mpot.resize(ncells);
         state->calc_refine_potential(mesh, mpot, icount, jcount);

         nlft.clear();
         nrht.clear();
         nbot.clear();
         ntop.clear();
      }  //  Complete CPU calculation.

      size_t result_size = 1;
      cl_mem dev_result = NULL;
      if (do_gpu_calc) {
         dev_mpot     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
         dev_result  = ezcl_malloc(NULL, &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
 
         state->gpu_calc_refine_potential(command_queue, mesh, dev_mpot, dev_result, dev_ioffset);

         ezcl_device_memory_remove(dev_nlft);
         ezcl_device_memory_remove(dev_nrht);
         ezcl_device_memory_remove(dev_nbot);
         ezcl_device_memory_remove(dev_ntop);
      }
      
      if (do_comparison_calc) {
         // Need to compare dev_mpot to mpot
         vector<int>mpot_save(ncells);
         ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot_save[0], NULL);
         for (uint ic = 0; ic < ncells; ic++){
            if (fabs(mpot[ic]-mpot_save[ic]) > STATE_EPS) {
               printf("DEBUG refine_potential at cycle %d mpot & mpot_save %d %d %d \n",n,ic,mpot[ic],mpot_save[ic]);
            }
         }
      }
      // Sync up cpu array with gpu version to reduce differences due to minor numerical differences
      // otherwise cell count will diverge causing code problems and crashes
      if (do_sync) {
         ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
      }
      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot[0], NULL);
      }

      int mcount, mtotal;
/*
// XXX Buggy here b/c gpu_rezone_count moved into gpu_calc_refine_potential XXX
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
            if (mcount != ioffset[ig]) printf("DEBUG ig %d ioffset %d mcount %d\n",ig,ioffset[ig],mcount);
            mtotal += mcount;
         }
      }
// XXX XXX
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
      }
// XXX XXX
*/
      if (do_cpu_calc) {
         new_ncells = old_ncells+mesh->rezone_count(mpot);
      }

/*      if (do_gpu_calc) {
         mesh->gpu_rezone_count(command_queue, block_size, local_work_size, dev_ioffset, dev_result);
      }
*/
      if (do_comparison_calc) {
         ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int),       &ioffset[0], NULL);
         mtotal = 0;
         for (uint ig=0; ig<(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
            mcount = 0;
            for (uint ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
                if (ic >= old_ncells) break;

                if (celltype[ic] == REAL_CELL){
                   mcount += mpot[ic] ? 4 : 1;
                } else {
                   mcount += mpot[ic] ? 2 : 1;
                }
            }
            if (mtotal != ioffset[ig]) printf("DEBUG ig %d ioffset %d mcount %d\n",ig,ioffset[ig],mtotal);
            mtotal += mcount;
         }
         ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int),       &ioffset[0], NULL);

         int result;
         ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, 1*sizeof(cl_int),       &result, NULL);
         new_ncells = result;
         //printf("Result is %d\n",result);
      }

      if (do_gpu_calc) {
         ezcl_device_memory_remove(dev_result);
      }

      if (do_cpu_calc) {
         int add_ncells = new_ncells - old_ncells;
         state->rezone_all(mesh, mpot, add_ncells);
         mpot.clear();
      }

      //  Resize the mesh, inserting cells where refinement is necessary.
      if (do_gpu_calc) {
         state->gpu_rezone_all(command_queue, mesh, ncells, new_ncells, old_ncells, localStencil, dev_mpot, dev_ioffset);
      }

      if (do_comparison_calc) {
         vector<real> H_save(ncells);
         vector<real> U_save(ncells);
         vector<real> V_save(ncells);
         vector<int> level_check(ncells);
         vector<int> celltype_check(ncells);
         vector<int> i_check(ncells);
         vector<int> j_check(ncells);
         /// Set read buffers for data.
         if (n % outputInterval == 0) {
            ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_real), &H_save[0], &start_read_event);
            state->gpu_time_read += ezcl_timer_calc(&start_read_event, &start_read_event);
         } else {
            ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_FALSE, 0, ncells*sizeof(cl_real), &H_save[0], NULL);
         }
         ezcl_enqueue_read_buffer(command_queue, dev_U,        CL_FALSE, 0, ncells*sizeof(cl_real), &U_save[0],          NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V,        CL_FALSE, 0, ncells*sizeof(cl_real), &V_save[0],          NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int),  &level_check[0],     NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  &celltype_check[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  &i_check[0],         NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_TRUE,  0, ncells*sizeof(cl_int),  &j_check[0],         NULL);
         for (uint ic = 0; ic < ncells; ic++){
            if (fabs(H[ic]-H_save[ic]) > STATE_EPS) printf("DEBUG diff at cycle %d H & H_save %d %lf %lf \n",n,ic,H[ic],H_save[ic]);
            if (fabs(U[ic]-U_save[ic]) > STATE_EPS) printf("DEBUG diff at cycle %d U & U_save %d %lf %lf \n",n,ic,U[ic],U_save[ic]);
            if (fabs(V[ic]-V_save[ic]) > STATE_EPS) printf("DEBUG diff at cycle %d V & V_save %d %lf %lf \n",n,ic,V[ic],V_save[ic]);
            if (level[ic] != level_check[ic] ) printf("DEBUG -- level: ic %d level %d level_check %d\n",ic, level[ic], level_check[ic]);
            if (celltype[ic] != celltype_check[ic] ) printf("DEBUG -- celltype: ic %d celltype %d celltype_check %d\n",ic, celltype[ic], celltype_check[ic]);
            if (i[ic] != i_check[ic] ) printf("DEBUG -- i: ic %d i %d i_check %d\n",ic, i[ic], i_check[ic]);
            if (j[ic] != j_check[ic] ) printf("DEBUG -- j: ic %d j %d j_check %d\n",ic, j[ic], j_check[ic]);
         }
      
      } //else if (n % outputInterval == 0) {
        // vector<real> H_save(ncells);
        // ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_real), &H_save[0], &start_read_event);
        // state->gpu_time_read             += ezcl_timer_calc(&start_read_event,       &start_read_event);
      //}

      if (do_gpu_calc) {
         mesh->gpu_count_BCs(command_queue, block_size, local_work_size, global_work_size, dev_ioffset);
      }

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

      ioffset.clear();

      ezcl_device_memory_remove(dev_ioffset);

      double H_sum = -1.0;

      if (do_comparison_calc) {
         H_sum = state->mass_sum(mesh, enhanced_precision_sum);
         double total_mass = state->gpu_mass_sum(command_queue, mesh, enhanced_precision_sum);
         if (fabs(total_mass - H_sum) > CONSERVATION_EPS) printf("Error: mass sum gpu %f cpu %f\n", total_mass, H_sum);/***/
      }

      mesh->proc.resize(ncells);
      if (icount)
      {  vector<int> index(ncells);
         cpu_timer_start(&tstart_cpu);
         mesh->partition_cells(numpe, mesh->proc, index, cycle_reorder);
         state->state_reorder(index);
         if (do_gpu_sync) {
            ezcl_enqueue_write_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&celltype[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&i[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int),  (void *)&j[0],  NULL);
            ezcl_enqueue_write_buffer(command_queue, dev_level, CL_TRUE,  0, ncells*sizeof(cl_int),  (void *)&level[0],  NULL);
         }
         cpu_time_partition += cpu_timer_stop(tstart_cpu); 
      }
      
      if (do_gpu_sync) {
         ezcl_enqueue_write_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&H[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_real),  (void *)&U[0],  NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_real),  (void *)&V[0],  NULL);
      }

      if (n % outputInterval == 0) {
         if (H_sum < 0) {
            H_sum = state->mass_sum(mesh, enhanced_precision_sum);
         }
         printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %14.12lg\n",
            n, deltaT, simTime, ncells, H_sum, H_sum - H_sum_initial);
#ifdef HAVE_GRAPHICS
         if (do_cpu_calc){
            mesh->calc_spatial_coordinates(0);
         }
         vector<real> x_save;
         vector<real> dx_save;
         vector<real> y_save;
         vector<real> dy_save;
         vector<real> H_save;
         if (do_gpu_calc){
            x_save.resize(ncells);
            dx_save.resize(ncells);
            y_save.resize(ncells);
            dy_save.resize(ncells);
            H_save.resize(ncells);
            cl_mem dev_x  = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
            cl_mem dev_dx = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
            cl_mem dev_y  = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
            cl_mem dev_dy = ezcl_malloc(NULL, &ncells, sizeof(cl_real),  CL_MEM_READ_WRITE, 0);
            mesh->gpu_calc_spatial_coordinates(command_queue, dev_x, dev_dx, dev_y, dev_dy);
            ezcl_enqueue_read_buffer(command_queue, dev_x,  CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&x_save[0],  &start_read_event);
            ezcl_enqueue_read_buffer(command_queue, dev_dx, CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&dx_save[0], NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_y,  CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&y_save[0],  NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_dy, CL_FALSE, 0, ncells*sizeof(cl_real), (void *)&dy_save[0], NULL);
            ezcl_enqueue_read_buffer(command_queue, state->dev_H, CL_TRUE,  0, ncells*sizeof(cl_real), (void *)&H_save[0],  &end_read_event);
            ezcl_device_memory_remove(dev_x);
            ezcl_device_memory_remove(dev_dx);
            ezcl_device_memory_remove(dev_y);
            ezcl_device_memory_remove(dev_dy);

            state->gpu_time_read += ezcl_timer_calc(&start_read_event, &end_read_event);
         }
         if (do_comparison_calc){
            for (uint ic = 0; ic < ncells; ic++){
               if (x[ic] != x_save[ic] || dx[ic] != dx_save[ic] || y[ic] != y_save[ic] || dy[ic] != dy_save[ic] ) {
                  printf("Error -- mismatch in spatial coordinates for cell %d is gpu %lf %lf %lf %lf cpu %lf %lf %lf %lf\n",ic,x_save[ic],dx_save[ic],y_save[ic],dy_save[ic],x[ic],dx[ic],y[ic],dy[ic]);
                  exit(0);
               }
            }
            for (uint ic = 0; ic < ncells; ic++){
               if (fabs(H[ic] - H_save[ic]) > CONSERVATION_EPS) {
                  printf("Error -- mismatch in H for cell %d is gpu %lf cpu %lf\n",ic,H_save[ic],H[ic]);
                  exit(0);
               }
            }
         }

         set_mysize(ncells);
         set_viewmode(view_mode);
         set_cell_coordinates(&x[0], &dx[0], &y[0], &dy[0]);
         set_cell_data(&H[0]);
         set_cell_proc(&mesh->proc[0]);
         set_circle_radius(circle_radius);
         draw_scene();
#endif
         output_flag = 1;
      }  //  Complete output interval.
      ++n;
      simTime += deltaT;
      
   }

   //  Output final results and timing information.
   if (n > niter) {
      //free_display();
      
      //  Get overall program timing.
      gettimeofday(&tstop, NULL);
      tresult.tv_sec = tstop.tv_sec - tstart.tv_sec;
      tresult.tv_usec = tstop.tv_usec - tstart.tv_usec;
      double elapsed_time = (double)tresult.tv_sec + (double)tresult.tv_usec*1.0e-6;
      
      state->output_timing_info(mesh, do_cpu_calc, do_gpu_calc, elapsed_time);

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      exit(0);
   }  //  Complete final output.
   
}

