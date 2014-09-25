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
#ifndef MESH_H_
#define MESH_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "MallocPlus/MallocPlus.h"
#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "kdtree/KDTree.h"
#include "crux/crux.h"
#include "partition.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif

#if !defined(FULL_PRECISION) && !defined(MIXED_PRECISION) && !defined(MINIMUM_PRECISION)
#define FULL_PRECISION
#endif
#ifdef NO_CL_DOUBLE
#undef  FULL_PRECISION
#undef  MIXED_PRECISION
#define MINIMUM_PRECISION
#endif

#if defined(MINIMUM_PRECISION)
   typedef float real_t; // this is used for intermediate calculations
   typedef float spatial_t; // for spatial variables
#ifdef HAVE_OPENCL
   typedef cl_float cl_real_t; // for intermediate gpu physics state variables
   typedef cl_float cl_spatial_t;
#endif
#ifdef HAVE_MPI
   #define MPI_REAL_T MPI_FLOAT // for MPI communication for physics state variables
   #define MPI_SPATIAL_T MPI_FLOAT
#endif

#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
   typedef double real_t;
   typedef float spatial_t; // for spatial variables
#ifdef HAVE_OPENCL
   typedef cl_double cl_real_t; // for intermediate gpu physics state variables
   typedef cl_float cl_spatial_t;
#endif
#ifdef HAVE_MPI
   #define MPI_REAL_T MPI_DOUBLE
   #define MPI_SPATIAL_T MPI_FLOAT
#endif

#elif defined(FULL_PRECISION)
   typedef double real_t;
   typedef double spatial_t; // for spatial variables
#ifdef HAVE_OPENCL
   typedef cl_double cl_real_t; // for intermediate gpu physics state variables
   typedef cl_double cl_spatial_t;
#endif
#ifdef HAVE_MPI
   #define MPI_REAL_T MPI_DOUBLE
   #define MPI_SPATIAL_T MPI_DOUBLE
#endif
#endif

#define TILE_SIZE 128

#define SWAP_PTR(xnew,xold,xtmp) (xtmp=xnew, xnew=xold, xold=xtmp)
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef unsigned int uint;

//float mem_opt_factor = 1.0;

enum boundary
{  REAL_CELL      =  1,         //  Denotes cell type of real cell.
   LEFT_BOUNDARY  = -1,         //  Denotes left boundary ghost cell.
   RIGHT_BOUNDARY = -2,         //  Denotes right boundary ghost cell.
   BOTTOM_BOUNDARY= -3,         //  Denotes bottom boundary ghost cell.
   TOP_BOUNDARY   = -4,         //  Denotes top boundary ghost cell.
   FRONT_BOUNDARY = -5,         //  Denotes front boundary ghost cell.
   BACK_BOUNDARY  = -6 };       //  Denotes back boundary ghost cell.

enum dimensionality
{  ONE_DIMENSIONAL   = 1,       // Dimensionality based at 1 for clarity.
   TWO_DIMENSIONAL,
   THREE_DIMENSIONAL};

enum orientation
{  SW,                          //  SW quadrant.
   NW,                          //  NW quadrant.
   NE,                          //  NE quadrant.
   SE };                        //  SE quadrant.

enum neighbor_calc
{  HASH_TABLE,                  //  Hash Table.
   KDTREE };                    //  kD-tree.

enum cpu_timers
{
   MESH_TIMER_COUNT_BCS,
   MESH_TIMER_CALC_NEIGHBORS,
   MESH_TIMER_HASH_SETUP,
   MESH_TIMER_HASH_QUERY,
   MESH_TIMER_FIND_BOUNDARY,
   MESH_TIMER_PUSH_SETUP,
   MESH_TIMER_PUSH_BOUNDARY,
   MESH_TIMER_LOCAL_LIST,
   MESH_TIMER_LAYER1,
   MESH_TIMER_LAYER2,
   MESH_TIMER_LAYER_LIST,
   MESH_TIMER_COPY_MESH_DATA,
   MESH_TIMER_FILL_MESH_GHOST,
   MESH_TIMER_FILL_NEIGH_GHOST,
   MESH_TIMER_SET_CORNER_NEIGH,
   MESH_TIMER_NEIGH_ADJUST,
   MESH_TIMER_SETUP_COMM,
   MESH_TIMER_KDTREE_SETUP,
   MESH_TIMER_KDTREE_QUERY,
   MESH_TIMER_REFINE_SMOOTH,
   MESH_TIMER_REZONE_ALL,
   MESH_TIMER_PARTITION,
   MESH_TIMER_CALC_SPATIAL_COORDINATES,
   MESH_TIMER_LOAD_BALANCE,
   MESH_TIMER_SIZE
};

#ifdef DEBUG_RESTORE_VALS
static const char *mesh_timer_descriptor[MESH_TIMER_SIZE] = {
   "mesh_cpu_timer_count_BCs",
   "mesh_cpu_timer_calc_neighbors",
   "mesh_cpu_timer_hash_setup",
   "mesh_cpu_timer_hash_query",
   "mesh_cpu_timer_find_boundary",
   "mesh_cpu_timer_push_setup",
   "mesh_cpu_timer_push_boundary",
   "mesh_cpu_timer_local_list",
   "mesh_cpu_timer_layer1",
   "mesh_cpu_timer_layer2",
   "mesh_cpu_timer_layer_list",
   "mesh_cpu_timer_copy_mesh_data",
   "mesh_cpu_timer_fill_mesh_ghost",
   "mesh_cpu_timer_fill_neigh_ghost",
   "mesh_cpu_timer_set_corner_neigh",
   "mesh_cpu_timer_neigh_adjust",
   "mesh_cpu_timer_setup_comm",
   "mesh_cpu_timer_kdtree_setup",
   "mesh_cpu_timer_kdtree_query",
   "mesh_cpu_timer_refine_smooth",
   "mesh_cpu_timer_rezone_all",
   "mesh_cpu_timer_partition",
   "mesh_cpu_timer_calc_spatial_coordinates",
   "mesh_cpu_timer_load_balance"
};
#endif

typedef enum cpu_timers mesh_timer_category;

using namespace std;

class Mesh
{

public:
   int ndim;                    //  Dimensionality of mesh (2 or 3).

   MallocPlus mesh_memory;
   MallocPlus gpu_mesh_memory;

#ifdef HAVE_OPENCL
   string defines;
#endif

   double cpu_timers[MESH_TIMER_SIZE];

   long   gpu_timers[MESH_TIMER_SIZE];

   int      cpu_rezone_counter;
   int      cpu_refine_smooth_counter;
   int      cpu_calc_neigh_counter;
   int      cpu_load_balance_counter;
   int      gpu_rezone_counter;
   int      gpu_refine_smooth_counter;
   int      gpu_calc_neigh_counter;
   int      gpu_load_balance_counter;

   int            mype,
                  numpe,
                  parallel,
                  cell_handle,
                  noffset;

   float          mem_factor;

   double         offtile_ratio_local;
   int            offtile_local_count;

   vector<int>    corners_i,
                  corners_j;

   vector<int>    nsizes,
                  ndispl;

   FILE          *fp;

   TKDTree        tree;         //  k-D tree for neighbor search.
   vector<int>    proc;
   vector<int>    lev_ibegin,   //  Lowest x-index in use at specified level of refinement.
                  lev_iend,     //  Highest x-index in use at specified level of refinement.
                  lev_jbegin,   //  Lowest y-index in use at specified level of refinement.
                  lev_jend,     //  Highest y-index in use at specified level of refinement.
                  lev_kbegin,   //  Lowest z-index in use at specified level of refinement.
                  lev_kend,     //  Highest z-index in use at specified level of refinement.
                  levtable;     //  Powers of two to simplify i,j calculations
   vector<real_t> lev_deltax,   //  Grid spacing along x-axis at specified level of refinement.
                  lev_deltay,   //  Grid spacing along y-axis at specified level of refinement.
                  lev_deltaz;   //  Grid spacing along z-axis at specified level of refinement.
   int            levmx,        //  Maximum level of refinement allowed.
                  have_boundary,//  Mesh includes boundary cells, else creates on the fly
                  ibase,        //  Index basis for arrays (0 for C, 1 for Fortan).
                  imin,         //  Lowest x-index in use.
                  imax,         //  Highest x-index in use.
                  jmin,         //  Lowest y-index in use.
                  jmax,         //  Highest y-index in use.
                  kmin,         //  Lowest z-index in use.
                  kmax;         //  Highest z-index in use.
   size_t         ncells,       //  Number of cells in mesh.
                  ncells_global, //  Global number of cells for parallel runs
                  ncells_ghost; //  Number of cells in mesh with ghost cells.
   real_t         xmin,         //  Lowest x-coordinate in use.
                  xmax,         //  Highest x-coordinate in use.
                  ymin,         //  Lowest y-coordinate in use.
                  ymax,         //  Highest y-coordinate in use.
                  zmin,         //  Lowest z-coordinate in use.
                  zmax,         //  Highest z-coordinate in use.
                  xcentermin,   //
                  xcentermax,   //
                  ycentermin,   //
                  ycentermax,   //
                  zcentermin,   //
                  zcentermax,   //
                  deltax,       //  Grid spacing along x-axis.
                  deltay,       //  Grid spacing along y-axis.
                  deltaz;       //  Grid spacing along z-axis.

   vector<int>    index;        //  1D ordered index of mesh elements.

                                 //  mesh state data
   int            *i,            //  1D array of mesh element x-indices.
                  *j,            //  1D array of mesh element y-indices.
                  *k,            //  1D array of mesh element z-indices.
                  *level,        //  1D array of mesh element refinement levels.
                                 //  derived data from mesh state data
                  *celltype,     //  1D ordered index of mesh element cell types (ghost or real).
                  *nlft,         //  1D ordered index of mesh element left neighbors.
                  *nrht,         //  1D ordered index of mesh element right neighbors.
                  *nbot,         //  1D ordered index of mesh element bottom neighbors.
                  *ntop,         //  1D ordered index of mesh element top neighbors.
                  *nfrt,         //  1D ordered index of mesh element front neighbors.
                  *nbak;         //  1D ordered index of mesh element back neighbors.

   vector<spatial_t> x,            //  1D ordered index of mesh element x-coordinates.
                     dx,           //  1D ordered index of mesh element x-coordinate spacings.
                     y,            //  1D ordered index of mesh element y-coordinates.
                     dy,           //  1D ordered index of mesh element y-coordinate spacings.
                     z,            //  1D ordered index of mesh element z-coordinates.
                     dz;           //  1D ordered index of mesh element z-coordinate spacings.

#ifdef HAVE_OPENCL
   cl_mem         dev_ioffset;

   cl_mem         dev_celltype,       
                  dev_i,       
                  dev_j,       
                  dev_level,       
                  dev_nlft,       
                  dev_nrht,       
                  dev_nbot,       
                  dev_ntop;       

   cl_mem         dev_levdx,    // corresponds to lev_deltax
                  dev_levdy,    // corresponds to lev_deltay
                  dev_levibeg,
                  dev_leviend,
                  dev_levjbeg,
                  dev_levjend,
                  dev_levtable; //

   cl_mem         dev_corners_i,
                  dev_corners_j;
#endif

   //   Public constructors.
   Mesh(FILE *fin, int *numpe);
   Mesh(int nx, int ny, int levmx_in, int ndim_in, double deltax_in, double deltay_in, int boundary, int parallel_in, int do_gpu_calc);

   //   Member functions.
   void init(int nx, int ny, real_t circ_radius, partition_method initial_order, int do_gpu_calc);
   void terminate(void);

/* memory routines */
   void allocate(size_t ncells);
   void resize(size_t new_ncells);
   void memory_reset_ptrs(void);
   void resize_old_device_memory(size_t ncells);

/* inline "macros" */
   int  is_lower_boundary(int *iv, int *lev_begin, int ic)    { return (iv[ic] < lev_begin[level[ic]]); }
   int  is_upper_boundary(int *iv, int *lev_end,   int ic)    { return (iv[ic] > lev_end[level[ic]]); }

   int  is_left_boundary(int ic)    { return (i[ic] < lev_ibegin[level[ic]]); }
   int  is_right_boundary(int ic)   { return (i[ic] > lev_iend[  level[ic]]); }
   int  is_bottom_boundary(int ic)  { return (j[ic] < lev_jbegin[level[ic]]); }
   int  is_top_boundary(int ic)     { return (j[ic] > lev_jend[  level[ic]]); }
   int  is_front_boundary(int ic)   { return (k[ic] < lev_kbegin[level[ic]]); }
   int  is_back_boundary(int ic)    { return (k[ic] > lev_kend[  level[ic]]); }

   int is_lower(int i)  { return(i % 2 == 0); }
   int is_upper(int i)  { return(i % 2 == 1); }

   int is_lower_left(int i, int j)  { return(i % 2 == 0 && j % 2 == 0); }
   int is_lower_right(int i, int j) { return(i % 2 == 1 && j % 2 == 0); }
   int is_upper_left(int i, int j)  { return(i % 2 == 0 && j % 2 == 1); }
   int is_upper_right(int i, int j) { return(i % 2 == 1 && j % 2 == 1); }

   int is_same_level_or_coarser(int nn, int nz) { return(level[nn] <= level[nz]); }
   int is_coarser(int nn, int nz)               { return(level[nn] <  level[nz]); }
   int is_finer(int nn, int nz)                 { return(level[nn] >  level[nz]); }
   int is_same_level(int nn, int nz)            { return(level[nn] == level[nz]); }


/* accessor routines */
   double get_cpu_timer(mesh_timer_category category)       {return(cpu_timers[category]); };
   long   get_gpu_timer(mesh_timer_category category)       {return(gpu_timers[category]); };

   int get_cpu_load_balance_count(void)           {return(cpu_load_balance_counter); };
   int get_cpu_rezone_count(void)                 {return(cpu_rezone_counter); };
   int get_cpu_refine_smooth_count(void)          {return(cpu_refine_smooth_counter); };
   int get_cpu_calc_neigh_count(void)             {return(cpu_calc_neigh_counter); };
   int get_gpu_load_balance_count(void)           {return(gpu_load_balance_counter); };
   int get_gpu_rezone_count(void)                 {return(gpu_rezone_counter); };
   int get_gpu_refine_smooth_count(void)          {return(gpu_refine_smooth_counter); };
   int get_gpu_calc_neigh_count(void)             {return(gpu_calc_neigh_counter); };

   int get_calc_neighbor_type(void);

   void print_partition_measure(void);
   void print_calc_neighbor_type(void);
   void print_partition_type(void);
/* end accessor routines */

/* Debugging, internal, or not used yet */
#ifdef HAVE_OPENCL
   int gpu_count_BCs();
#endif
   void kdtree_setup(void);
   void partition_measure(void);
   void partition_cells(int numpe,
                   vector<int> &order,
                   enum partition_method method);
   void calc_distribution(int numpe);
   void calc_symmetry(vector<int> &dsym,
                  vector<int> &xsym,
                  vector<int> &ysym);

/* End of debugging, internal, or not used yet */

   /**************************************************************************************
   * Calculate neighbors
   *  Input -- from within the object
   *    i, j, level
   *  Output -- in the object
   *    nlft, nrht, nbot, ntop arrays
   **************************************************************************************/
   void calc_neighbors(int ncells);
   void calc_neighbors_local(void);
#ifdef HAVE_OPENCL
   void gpu_calc_neighbors(void);
   void gpu_calc_neighbors_local(void);
#endif
   //   TODO:  Not created yet; overloading for 3D mesh support. (davis68)
   void calc_neighbors(vector<int> &nlft,
                  vector<int> &nrht,
                  vector<int> &nbot,
                  vector<int> &ntop,
                  vector<int> &nfrt,
                  vector<int> &nbak,
                  vector<int> index);

   /**************************************************************************************
   * Calculate rezone count
   *  Input
   *    mpot -- potential mesh refinement
   *    ioffset -- write offset for each cell
   *  Output
   *    result -- cell count
   **************************************************************************************/
   int  rezone_count(vector<int> mpot, int &icount, int &jcount);
#ifdef HAVE_OPENCL
   void gpu_rezone_count2(size_t block_size, size_t local_work_size, cl_mem dev_redscratch, cl_mem &dev_result);
   void gpu_rezone_count(size_t block_size, size_t local_work_size, cl_mem dev_redscratch, cl_mem &dev_result);
   void gpu_rezone_scan(size_t block_size, size_t local_work_size, cl_mem dev_ioffset, cl_mem &dev_result);
#endif

   /**************************************************************************************
   * Refine Smooth -- smooths jump in refinement level so that only a 1 to 2 jump occurs
   *  Input/Output
   *    mpot -- potential mesh refinement array, 1 is refine and -1 coarsen
   *    ioffset -- write offset for each cell to account for new cells
   *    result -- refinement count
   **************************************************************************************/
   size_t refine_smooth(vector<int> &mpot, int &icount, int &jcount);
#ifdef HAVE_OPENCL
   int gpu_refine_smooth(cl_mem &dev_mpot, int &icount, int &jcount);
#endif

   /**************************************************************************************
   * Rezone mesh
   *  Input
   *     add_ncells -- for each processor. A global sum will be done and the main part of
   *        the rezone will be skipped if no cells are added.
   *     mpot -- mesh rezone potential
   *     have_state flag -- 0 (false) for setup when physics state has not been allocated
   *     ioffset -- partial prefix scan results for starting address to write new cells
   *     state_memory -- linked list of arrays for state
   *  Output
   *     new mesh and state arrays with refinement/coarsening performed
   **************************************************************************************/
   void rezone_all(int icount, int jcount, vector<int> mpot, int have_state, MallocPlus &state_memory);
#ifdef HAVE_OPENCL
   void gpu_rezone_all(int icount, int jcount, cl_mem &dev_mpot, MallocPlus &gpu_state_memory);
#endif

   /**************************************************************************************
   * Load balance -- only needed for parallel (MPI) runs
   *  Input
   *    numcells -- ncells from rezone all routine. This is a copy in so that a local
   *       value can be used for load_balance and gpu_load_balance without it getting
   *       reset for clamr_checkall routine
   *    weight -- weighting array per cell for balancing. Currently not used. Null value
   *       indicates even weighting of cells for load balance. 
   *    state_memory or gpu_state_memory -- linked-list of arrays from physics routine
   *       to be load balanced. 
   * Output -- arrays will be returned load balanced with new sizes. Pointers to arrays
   *       will need to be reset
   **************************************************************************************/
#ifdef HAVE_MPI
   void do_load_balance_local(size_t numcells, float *weight, MallocPlus &state_memory);
#ifdef HAVE_OPENCL
   int gpu_do_load_balance_local(size_t numcells, float *weight, MallocPlus &gpu_state_memory);
#endif
#endif

   /**************************************************************************************
   * Calculate spatial coordinates
   *  Input -- from within the object
   *    i, j, level
   *  Output
   *    x, y -- coordinates for each cell
   *    dx, dy -- size of each cell
   **************************************************************************************/
   void calc_spatial_coordinates(int ibase);
#ifdef HAVE_OPENCL
   void gpu_calc_spatial_coordinates(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy);
#endif

   /**************************************************************************************
   * Testing routines
   **************************************************************************************/
#ifdef HAVE_OPENCL
   void compare_dev_local_to_local(void); // Not currently called
   void compare_neighbors_gpu_global_to_cpu_global(void);
#endif
   void compare_neighbors_cpu_local_to_cpu_global(uint ncells_ghost, uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl);
#ifdef HAVE_OPENCL
   void compare_neighbors_all_to_gpu_local(Mesh *mesh_global, int *nsizes, int *ndispl);
   void compare_mpot_gpu_global_to_cpu_global(int *mpot, cl_mem dev_mpot);
#endif
   void compare_mpot_cpu_local_to_cpu_global(uint ncells_global, int *nsizes, int *displ, int *mpot, int *mpot_global, int cycle);
#ifdef HAVE_OPENCL
   void compare_mpot_all_to_gpu_local(int *mpot, int *mpot_global, cl_mem dev_mpot, cl_mem dev_mpot_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle);
   void compare_ioffset_gpu_global_to_cpu_global(uint old_ncells, int *mpot);
   void compare_ioffset_all_to_gpu_local(uint old_ncells, uint old_ncells_global, int block_size, int block_size_global, int *mpot, int *mpot_global, cl_mem dev_ioffset, cl_mem dev_ioffset_global, int *ioffset, int *ioffset_global, int *celltype_global, int *i_global, int *j_global);
   void compare_coordinates_gpu_global_to_cpu_global_double(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, double *H);
   void compare_coordinates_gpu_global_to_cpu_global_float(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, float *H);
#endif
   void compare_coordinates_cpu_local_to_cpu_global_double(uint ncells_global, int *nsizes, int *ndispl, spatial_t *x, spatial_t *dx, spatial_t *y, spatial_t *dy, double *H, spatial_t *x_global, spatial_t *dx_global, spatial_t *y_global, spatial_t *dy_global, double *H_global, int cycle);
   void compare_coordinates_cpu_local_to_cpu_global_float(uint ncells_global, int *nsizes, int *ndispl, spatial_t *x, spatial_t *dx, spatial_t *y, spatial_t *dy, float *H, spatial_t *x_global, spatial_t *dx_global, spatial_t *y_global, spatial_t *dy_global, float *H_global, int cycle);
#ifdef HAVE_OPENCL
   void compare_indices_gpu_global_to_cpu_global(void);
#endif
   void compare_indices_cpu_local_to_cpu_global(uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl, int cycle);
#ifdef HAVE_OPENCL
   void compare_indices_all_to_gpu_local(Mesh *mesh_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle);
#endif

   size_t get_checkpoint_size(void);
   void store_checkpoint(Crux *crux);
   void restore_checkpoint(Crux *crux);

   void calc_celltype(size_t ncells);

private:
   //   Private constructors.
   Mesh(const Mesh&);   //   Blocks copy constructor so copies are not made inadvertently.

   //   Member functions.
   void print_object_info();

   void set_refinement_order(int order[4], int ic, int ifirst, int ilast, int jfirst, int jlast,
                                int level_first, int level_last, int *i, int *j, int *level);

   void write_grid(int ncycle);
   void calc_centerminmax(void);
   void calc_minmax(void);

   void print(void);
   void print_local(void);
#ifdef HAVE_OPENCL
   void print_dev_local();
#endif

};

#endif /* MESH_H */
