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
#include "input.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
#include "state.h"
#include "l7/l7.h"
#include "timer/timer.h"
#include "memstats/memstats.h"

#include <mpi.h>
#include <omp.h>
#include "display.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef HAVE_QUO
#include "QUO.hpp"
#endif

static int do_cpu_calc = 1;
static int do_gpu_calc = 0;

#ifdef HAVE_CL_DOUBLE
typedef double      real;
typedef struct
{
   double s0;
   double s1;
}  real2;
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

static double H_sum_initial = 0.0;
static double cpu_time_graphics = 0.0;
double cpu_time_main_setup = 0.0;
vector<real> H_global;
vector<real> x_global;
vector<real> dx_global;
vector<real> y_global;
vector<real> dy_global;
vector<int> proc_global;

typedef struct SubComm {
    // the communicator
    MPI_Comm comm;
    // my rank in the communicator
    int rank;
    // the size of the communicator
    int size;
    // cache whether or not i'm a member of the communicator
    bool member;
} SubComm;

typedef struct Context {
#ifdef HAVE_QUO
    QUO *quo;
#endif
    // MPI_COMM_WORLD rank
    int cwRank;
    // size of MPI_COMM_WORLD
    int cwSize;
    // sub communicator
    SubComm subComm;
} Context;

// my context
Context context;

/* ////////////////////////////////////////////////////////////////////////// */
static int
initQUO(Context &c)
{
#ifdef LIBQUO
   try {
       // init QUO -- all MPI processes MUST do this at the same time.
       QUO *q = new QUO();
       q->create();
       // stash this in the context
       c.quo = q;

       // emit some node info... not needed for init, just an example
       if (0 == c.cwRank) {
           std::cout << "### System Info ###" << std::endl;
           std::cout << "# Nodes: " << q->nnodes() << std::endl;
           std::cout << "# NUMA Nodes: " << q->nnumanodes() << std::endl;
           std::cout << "# Sockets: " << q->nsockets() << std::endl;
           std::cout << "# Cores: " << q->ncores() << std::endl;
           std::cout << "# PUs: " << q->npus() << std::endl;
           std::cout << "# MPI Procs on Node: " << q->nqids() << std::endl;
       }
   }
   catch (QUOException &e) {
        cerr << e.what() << endl;
        return 1;
   }
   return 0;
#else
   return 0;
#endif
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
contextSetup(Context &c,
             int cwRank,
             int cwSize)
{
    int rc = 0;
    // first stash some of my MPI_COMM_WORLD info
    c.cwRank = cwRank;
    c.cwSize = cwSize;
    // note that comm info is only valid if a member of the communicator
    c.subComm.rank = -1;
    c.subComm.size= -1;
    c.subComm.member = false;
#ifdef LIBQUO
    // init QUO
    rc = initQUO(c);
#endif
    return rc;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
finalizeQUO(Context &c)
{
#ifdef LIBQUO
    try {
        c.quo->free();
        delete c.quo;
    }
   catch (QUOException &e) {
        cerr << e.what() << endl;
        return 1;
   }
    return 0;
#else
    return 0;
#endif
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
getSubCommProcs(Context &c, int *vLen, int **v)
{
    QUO_obj_type_t resPrio[] = {QUO_OBJ_NODE,
                                QUO_OBJ_SOCKET};
    // default target resource is a NUMA node
    QUO_obj_type_t targetRes = QUO_OBJ_NODE;
    // number of target resources and max number of procs per resource
    int nRes = 0, maxProcPerRes = 1;
    int res_assigned = 0;
    int totalWorkers = 0;
    int rc = MPI_SUCCESS;
    // array that hold whether or not a particular rank is going to do work
    int *workContribs = NULL;
    // MPI_COMM_WORLD ranks of the selected workers
    int *workerRanks = NULL;

    try {
        // figure out what we are going to distribute work over
        for (int i = 0; i < sizeof(resPrio) / sizeof(resPrio[0]); ++i) {
            if ((nRes = c.quo->nObjsByType(resPrio[i])) > 0) {
                targetRes = resPrio[i];
                break;
            }
        }
        // failure -- fix this path at some point
        if (0 == nRes) return 1;
        /* let quo distribute workers over the sockets. if res_assigned is 1
         * after this call, then i have been chosen. */
        if (c.quo->autoDistrib(targetRes, maxProcPerRes)) {
            res_assigned = 1;
        }
    }
   catch (QUOException &e) {
        cerr << e.what() << endl;
        return 1;
   }
   /* array that hold whether or not a particular rank is going to do work */
   workContribs = (int *)calloc(c.cwSize, sizeof(*workContribs));
   if (!workContribs) return 1;

   if (MPI_SUCCESS != (rc = MPI_Allgather(&res_assigned, 1, MPI_INT,
                                          workContribs, 1, MPI_INT,
                                          MPI_COMM_WORLD))) {
       return 1;
   }
   /* now iterate over the array and count the total number of workers */
   for (int i = 0; i < c.cwSize; ++i) {
       if (1 == workContribs[i]) ++totalWorkers;
   }
   workerRanks = (int *)calloc(totalWorkers, sizeof(*workerRanks));
   if (!workerRanks) return 1;
   /* populate the array with the worker comm world ranks */
   for (int i = 0, j = 0; i < c.cwSize; ++i) {
       if (1 == workContribs[i]) {
           workerRanks[j++] = i;
       }
   }
   *vLen = totalWorkers;
   *v = workerRanks;
   if (workContribs) free(workContribs);
   return 0;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
subCommInit(Context &c,
            int np1s /* number of participants |p1who| */,
            int *p1who /* the participating ranks (MPI_COMM_WORLD) */)
{
    int rc = QUO_SUCCESS;
    MPI_Group worldGroup;
    MPI_Group p1_group;
    /* ////////////////////////////////////////////////////////////////////// */
    /* now create our own communicator based on the rank ids passed here */
    /* ////////////////////////////////////////////////////////////////////// */
    if (MPI_SUCCESS != MPI_Comm_group(MPI_COMM_WORLD, &worldGroup)) {
        rc = QUO_ERR_MPI;
        goto out;
    }
    if (MPI_SUCCESS != MPI_Group_incl(worldGroup, np1s,
                                      p1who, &p1_group)) {
        rc = QUO_ERR_MPI;
        goto out;
    }
    if (MPI_SUCCESS != MPI_Comm_create(MPI_COMM_WORLD,
                                       p1_group,
                                       &(c.subComm.comm))) {
        rc = QUO_ERR_MPI;
        goto out;
    }
    /* am i in the new communicator? */
    c.subComm.member = (MPI_COMM_NULL == c.subComm.comm) ? false : true;
    if (c.subComm.member) {
        if (MPI_SUCCESS != MPI_Comm_size(c.subComm.comm, &c.subComm.size)) {
            rc = QUO_ERR_MPI;
            goto out;
        }
        if (MPI_SUCCESS != MPI_Comm_rank(c.subComm.comm, &c.subComm.rank)) {
            rc = QUO_ERR_MPI;
            goto out;
        }
    }
out:
    if (MPI_SUCCESS != MPI_Group_free(&worldGroup)) return 1;
    return (QUO_SUCCESS == rc) ? 0 : 1;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
subCommCreate(Context &c)
{
    int rc = 1;
    int numPEsInSubComm = 0;
    // comm world ranks in subcomm
    int *cwPEs = NULL;
    rc = getSubCommProcs(c, &numPEsInSubComm, &cwPEs);
    if (rc) return rc;
    if (0 == c.quo->id()) {
        // add hostname -- could be diff on diff systems
        std::cout << "MPI_COMM_WORLD ranks in subComm: ";
        for (int i = 0; i < numPEsInSubComm; ++i) {
            std::cout << cwPEs[i] << " ";
        }
        std::cout << std::endl;
    }
    // now actually create the darn thing...
    rc = subCommInit(c, numPEsInSubComm, cwPEs);
    if (cwPEs) free(cwPEs);
    return rc;
}

int main(int argc, char **argv) {

   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   L7_Init(&mype, &numpe, &argc, argv);
   parseInput(argc, argv);

   struct timeval tstart_setup;
   cpu_timer_start(&tstart_setup);

   double circ_radius = 6.0;
   //  Scale the circle appropriately for the mesh size.
   circ_radius = circ_radius * (double) nx / 128.0;
   int boundary = 1;
   int parallel_in = 1;

   int nt = 0;

   // figure out the max number of threads that can be spawned
   nt = omp_get_max_threads();

   // setup our context information
   if (contextSetup(context, mype, numpe)) {
       fprintf(stderr, "(%d) contextSetup Failure!", mype);
   }
   // create the sub comm
   if (subCommCreate(context)) {
       fprintf(stderr, "(%d) subCommCreate Failure!", mype);
   }
   if (context.subComm.member) {
       std::cout << "rank: " << context.cwRank << " is subComm rank: "
                 << context.subComm.rank << std::endl;
   }
   // at this point subComm is ready to use for those inSubComm
   if (0 == mype) {
        printf("--- num openmp threads: %d\n", nt);
        int world_size, subcomm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_size(context.subComm.comm, &subcomm_size);
        printf("Size comm world %d subComm %d\n",world_size, subcomm_size);
        fflush(stdout);
   }

   mesh = new Mesh(nx, ny, levmx, ndim, numpe, boundary, parallel_in, do_gpu_calc);
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

   vector<int>   &nsizes     = mesh->nsizes;
   vector<int>   &ndispl     = mesh->ndispl;

   vector<int>   &celltype = mesh->celltype;
   vector<int>   &i        = mesh->i;
   vector<int>   &j        = mesh->j;
   vector<int>   &level    = mesh->level;

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

   celltype.resize(ncells);
   level.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);

   //H.resize(ncells);
   //U.resize(ncells);
   //V.resize(ncells);
   state->resize(ncells);

   state->fill_circle(circ_radius, 100.0, 5.0);

   mesh->nlft.clear();
   mesh->nrht.clear();
   mesh->nbot.clear();
   mesh->ntop.clear();

   x.clear();
   dx.clear();
   y.clear();
   dy.clear();

   //  Kahan-type enhanced precision sum implementation.
   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (mype == 0) printf ("Mass of initialized cells equal to %14.12lg\n", H_sum);
   H_sum_initial = H_sum;

   double cpu_time_main_setup = cpu_timer_stop(tstart_setup);
   state->parallel_timer_output(numpe,mype,"CPU:  setup time               time was",cpu_time_main_setup);

   long long mem_used = memstats_memused();
   if (mem_used > 0) {
      state->parallel_memory_output(numpe,mype,"Memory used      in startup ",mem_used);
      state->parallel_memory_output(numpe,mype,"Memory peak      in startup ",memstats_mempeak());
      state->parallel_memory_output(numpe,mype,"Memory free      at startup ",memstats_memfree());
      state->parallel_memory_output(numpe,mype,"Memory available at startup ",memstats_memtotal());
   }

   if (mype == 0) {
      printf("Iteration   0 timestep      n/a Sim Time      0.0 cells %ld Mass Sum %14.12lg\n", ncells_global, H_sum);
   }

   mesh->cpu_calc_neigh_counter=0;
   mesh->cpu_time_calc_neighbors=0.0;
   mesh->cpu_rezone_counter=0;
   mesh->cpu_refine_smooth_counter=0;

#ifdef HAVE_GRAPHICS
#ifdef HAVE_OPENGL
   set_mysize(ncells_global);
   //vector<real> H_global;
   //vector<real> x_global;
   //vector<real> dx_global;
   //vector<real> y_global;
   //vector<real> dy_global;
   //vector<int> proc_global;
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
   MPI_Gatherv(&state->H[0], nsizes[mype], MPI_C_REAL, &H_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);

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
   set_cell_data(&state->H[0]);
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
{  double g     = 9.80;
   double sigma = 0.95; 
   int icount, jcount;
   struct timeval tstart_cpu;

   //  Initialize state variables for GPU calculation.
   int &mype  = mesh->mype;
   int &numpe = mesh->numpe;

   //int levmx        = mesh->levmx;
   size_t &ncells_global    = mesh->ncells_global;
   size_t &ncells           = mesh->ncells;
   size_t &ncells_ghost     = mesh->ncells_ghost;

   vector<int>     mpot;
   vector<int>     mpot_global;
   
   size_t old_ncells = ncells;
   size_t old_ncells_global = ncells_global;
   size_t new_ncells = 0;
   double deltaT = 0.0;

   //  Main loop.
   for (int nburst = 0; nburst < outputInterval && ncycle < niter; nburst++, ncycle++) {

      //  Define basic domain decomposition parameters for GPU.
      old_ncells = ncells;
      old_ncells_global = ncells_global;

      MPI_Barrier(MPI_COMM_WORLD);
      cpu_timer_start(&tstart_cpu);
      //  Calculate the real time step for the current discrete time step.
      deltaT = state->set_timestep(g, sigma);
      simTime += deltaT;

      cpu_timer_start(&tstart_cpu);
      if (mesh->nlft.size() == 0) mesh->calc_neighbors_local();

      mesh->partition_measure();

      // Currently not working -- may need to be earlier?
      //if (mesh->have_boundary) {
      //  state->add_boundary_cells();
      //}

      // Apply BCs is currently done as first part of gpu_finite_difference and so comparison won't work here

      //  Execute main kernel
      cpu_timer_start(&tstart_cpu);
      state->calc_finite_difference(deltaT);

      //  Size of arrays gets reduced to just the real cells in this call for have_boundary = 0
      state->remove_boundary_cells();

      cpu_timer_start(&tstart_cpu);
      mpot.resize(ncells_ghost);
      new_ncells = state->calc_refine_potential(mpot, icount, jcount);
  
      cpu_timer_start(&tstart_cpu);
      int add_ncells = new_ncells - old_ncells;
      state->rezone_all(icount, jcount, mpot);
      mpot.clear();
      ncells = new_ncells;
      mesh->ncells = new_ncells;

      cpu_timer_start(&tstart_cpu);
      if (mesh->nlft.size() == 0) {
         state->do_load_balance_local(new_ncells);
      }


// XXX
//      mesh->proc.resize(ncells);
//      if (icount) {
//         vector<int> index(ncells);
//         mesh->partition_cells(numpe, index, cycle_reorder);
//      }

   } // End burst loop

   double H_sum = state->mass_sum(enhanced_precision_sum);
   if (isnan(H_sum)) {
      printf("Got a NAN on cycle %d\n",ncycle);
      exit(-1);
   }
   if (mype == 0){
      printf("Iteration %3d timestep %lf Sim Time %lf cells %ld Mass Sum %14.12lg Mass Change %12.6lg\n",
         ncycle, deltaT, simTime, ncells_global, H_sum, H_sum - H_sum_initial);
   }

#ifdef HAVE_GRAPHICS
   mesh->x.resize(ncells);
   mesh->dx.resize(ncells);
   mesh->y.resize(ncells);
   mesh->dy.resize(ncells);
   mesh->calc_spatial_coordinates(0);

   cpu_timer_start(&tstart_cpu);

#ifdef HAVE_MPE
   set_mysize(ncells);
   set_cell_coordinates(&mesh->x[0], &mesh->dx[0], &mesh->y[0], &mesh->dy[0]);
   set_cell_data(&state->H[0]);
   set_cell_proc(&mesh->proc[0]);
#endif
#ifdef HAVE_OPENGL
   vector<int>   &nsizes   = mesh->nsizes;
   vector<int>   &ndispl   = mesh->ndispl;

   set_mysize(ncells_global);
   //vector<real> x_global;
   //vector<real> dx_global;
   //vector<real> y_global;
   //vector<real> dy_global;
   //vector<real> H_global;
   //vector<int> proc_global;

   if (mype == 0) {
      x_global.resize(ncells_global);
      dx_global.resize(ncells_global);
      y_global.resize(ncells_global);
      dy_global.resize(ncells_global);
      H_global.resize(ncells_global);
      proc_global.resize(ncells_global);
   }
   MPI_Gatherv(&mesh->x[0],  nsizes[mype], MPI_C_REAL, &x_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&mesh->dx[0], nsizes[mype], MPI_C_REAL, &dx_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&mesh->y[0],  nsizes[mype], MPI_C_REAL, &y_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&mesh->dy[0], nsizes[mype], MPI_C_REAL, &dy_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   MPI_Gatherv(&state->H[0], nsizes[mype], MPI_C_REAL, &H_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);

   if (view_mode == 0) {
      mesh->proc.resize(ncells);
      for (size_t ii = 0; ii<ncells; ii++){
         mesh->proc[ii] = mesh->mype;
      }
   
      MPI_Gatherv(&mesh->proc[0],  nsizes[mype], MPI_INT, &proc_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, 0, MPI_COMM_WORLD);
   }

   set_cell_coordinates(&x_global[0], &dx_global[0], &y_global[0], &dy_global[0]);
   set_cell_data(&H_global[0]);
   set_cell_proc(&proc_global[0]);
#endif
   set_viewmode(view_mode);
   set_circle_radius(circle_radius);
   draw_scene();

   MPI_Barrier(MPI_COMM_WORLD);

   cpu_time_graphics += cpu_timer_stop(tstart_cpu);
#endif

   //  Output final results and timing information.
   if (ncycle >= niter) {
      //free_display();
      
      //  Get overall program timing.
      double elapsed_time = cpu_timer_stop(tstart);
      
      long long mem_used = memstats_memused();
      if (mem_used > 0) {
         state->parallel_memory_output(numpe,mype,"Memory used      ",mem_used);
         state->parallel_memory_output(numpe,mype,"Memory peak      ",memstats_mempeak());
         state->parallel_memory_output(numpe,mype,"Memory free      ",memstats_memfree());
         state->parallel_memory_output(numpe,mype,"Memory available ",memstats_memtotal());
      }
      state->output_timing_info(do_cpu_calc, do_gpu_calc, elapsed_time);

      state->parallel_timer_output(numpe,mype,"CPU:  graphics                 time was",cpu_time_graphics);

      mesh->print_partition_measure();
      mesh->print_calc_neighbor_type();
      mesh->print_partition_type();

      if (mype ==0) {
         printf("CPU:  rezone frequency                \t %8.4f\tpercent\n",     (double)mesh->get_cpu_rezone_count()/(double)ncycle*100.0 );
         printf("CPU:  calc neigh frequency            \t %8.4f\tpercent\n",     (double)mesh->get_cpu_calc_neigh_count()/(double)ncycle*100.0 );
         printf("CPU:  load balance frequency          \t %8.4f\tpercent\n",     (double)mesh->get_cpu_load_balance_count()/(double)ncycle*100.0 );
         printf("CPU:  refine_smooth_iter per rezone   \t %8.4f\t\n",            (double)mesh->get_cpu_refine_smooth_count()/(double)mesh->get_cpu_rezone_count() );
      }

      mesh->terminate();
      state->terminate();

      if (finalizeQUO(context)) {
          fprintf(stderr, "(%d) finalizeQUO failure\n", context.cwRank);
      }
      L7_Terminate();
      exit(0);
   }  //  Complete final output.
   
}

