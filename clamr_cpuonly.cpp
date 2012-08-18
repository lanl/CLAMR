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

static int do_gpu_calc = 0;
static int do_cpu_calc = 1;

#ifdef HAVE_CL_DOUBLE
typedef double      real;
typedef cl_double   cl_real;
#define CONSERVATION_EPS    .02
#else
typedef float       real;
typedef cl_float    cl_real;
#define CONSERVATION_EPS    .1
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

static double  H_sum_initial = 0.0;
static double  cpu_time_graphics = 0.0;

int main(int argc, char **argv) {

    // Needed for code to compile correctly on the Mac
   int mype=0;
   int numpe=-1;

   //  Process command-line arguments, if any.
   parseInput(argc, argv);
   
   numpe = 16;

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 0;
   
   mesh  = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
   mesh->init(nx, ny, circ_radius, initial_order, special_case, do_gpu_calc);
   size_t &ncells = mesh->ncells;
   state = new State(ncells);
   state->init(ncells, do_gpu_calc);
   mesh->proc.resize(ncells);
   mesh->calc_distribution(numpe, mesh->proc);
   state->fill_circle(mesh, circ_radius, 100.0, 5.0);
   
   vector<real>  &H        = state->H;

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

   //  Initialize state variables for GPU calculation.
   size_t &ncells    = mesh->ncells;

   vector<real>  &H        = state->H;

   vector<int>     mpot;
   
   size_t old_ncells = ncells;
   size_t new_ncells = 0;
   double H_sum = -1.0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      old_ncells = ncells;

      //  Calculate the real time step for the current discrete time step.
      deltaT = state->set_timestep(mesh, g, sigma);
      simTime += deltaT;
      
      mesh->calc_neighbors();

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (do_cpu_calc && ! mesh->have_boundary) {
      //  state->add_boundary_cells(mesh);
      //}

      // Need ghost cells for this routine
      state->apply_boundary_conditions(mesh);
      
      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      state->calc_finite_difference(mesh, deltaT);
      
      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells(mesh);
      
      //  Check for NANs.
      for (uint ic=0; ic<ncells; ic++) {
         if (isnan(H[ic])) {
            printf("Got a NAN on cell %d cycle %d\n",ic,ncycle);
            //  Release kernels and finalize the OpenCL elements.
            exit(-1);
         }
      }  //  Complete NAN check.

      mpot.resize(ncells);
      new_ncells = state->calc_refine_potential(mesh, mpot, icount, jcount);

      //new_ncells = old_ncells+mesh->rezone_count(mpot);
      //printf("icount %d old_ncells %d new_ncells %d\n",icount,old_ncells,new_ncells);

      //  Resize the mesh, inserting cells where refinement is necessary.

      int add_ncells = new_ncells - old_ncells;
      state->rezone_all(mesh, mpot, add_ncells);
      mpot.clear();
      mesh->ncells = new_ncells;
      ncells = new_ncells;

      mesh->proc.resize(ncells);
      if (icount)
      {  vector<int> index(ncells);
         mesh->partition_cells(numpe, mesh->proc, index, cycle_reorder);
         state->state_reorder(index);
      }
      
      mesh->ncells = ncells;
      
   }

   H_sum = state->mass_sum(mesh, enhanced_precision_sum);
   printf("Iteration %d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
      ncycle, deltaT, simTime, ncells, H_sum, H_sum - H_sum_initial);

   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

#ifdef HAVE_GRAPHICS
   mesh->calc_spatial_coordinates(0);

   set_mysize(ncells);
   set_viewmode(view_mode);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_data(&H[0]);
   set_cell_proc(&mesh->proc[0]);
   set_circle_radius(circle_radius);
   draw_scene();
#endif

   cpu_time_graphics += cpu_timer_stop(tstart_cpu);

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      state->output_timing_info(mesh, do_cpu_calc, do_gpu_calc, elapsed_time);
      printf("CPU:  graphics                 time was\t %8.4f\ts\n",     cpu_time_graphics );

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      delete mesh;
      delete state;

      exit(0);
   }  //  Complete final output.
   
}

