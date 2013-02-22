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

#ifdef HAVE_MPI_AND_OPENCL
#define HAVE_MPI    1
#define HAVE_OPENCL 1
#endif

#include <stdio.h>
#include <vector>
#include <math.h>
#include "kdtree/KDTree.h"
#include "partition.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif
#include "MallocPlus/MallocPlus.h"

#define TILE_SIZE 128

#define SWAP_PTR(xnew,xold,xtmp) (xtmp=xnew, xnew=xold, xold=xtmp)
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef HAVE_CL_DOUBLE
typedef double      real;
#else
typedef float       real;
#endif

typedef unsigned int uint;

enum boundary
{  REAL_CELL      =  1,         //  Denotes cell type of real cell.
   LEFT_BOUNDARY  = -1,         //  Denotes left boundary ghost cell.
   RIGHT_BOUNDARY = -2,         //  Denotes right boundary ghost cell.
   BOTTOM_BOUNDARY= -3,         //  Denotes bottom boundary ghost cell.
   TOP_BOUNDARY   = -4,         //  Denotes top boundary ghost cell.
   FRONT_BOUNDARY = -5,         //  Denotes front boundary ghost cell.
   BACK_BOUNDARY  = -6 };       //  Denotes back boundary ghost cell.

enum orientation
{  SW,                          //  SW quadrant.
   NW,                          //  NW quadrant.
   NE,                          //  NE quadrant.
   SE };                        //  SE quadrant.

enum neighbor_calc
{  HASH_TABLE,                  //  Hash Table.
   KDTREE };                    //  kD-tree.

using namespace std;

class Mesh
{  int ndim;                    //  Dimensionality of mesh (2 or 3).

public:
   double   cpu_time_calc_neighbors,
               cpu_time_hash_setup,
               cpu_time_hash_query,
               cpu_time_find_boundary,
               cpu_time_gather_boundary,
               cpu_time_local_list,
               cpu_time_layer1,
               cpu_time_layer2,
               cpu_time_layer_list,
               cpu_time_copy_mesh_data,
               cpu_time_fill_mesh_ghost,
               cpu_time_fill_neigh_ghost,
               cpu_time_set_corner_neigh,
               cpu_time_neigh_adjust,
               cpu_time_setup_comm,

               cpu_time_kdtree_setup,
               cpu_time_kdtree_query,
            cpu_time_rezone_all,
            cpu_time_partition,
            cpu_time_calc_spatial_coordinates,
            cpu_time_load_balance;

   long     gpu_time_calc_neighbors,
               gpu_time_hash_setup,
               gpu_time_hash_query,
               gpu_time_find_boundary,
               gpu_time_gather_boundary,
               gpu_time_local_list,
               gpu_time_layer1,
               gpu_time_layer2,
               gpu_time_layer_list,
               gpu_time_copy_mesh_data,
               gpu_time_fill_mesh_ghost,
               gpu_time_fill_neigh_ghost,
               gpu_time_set_corner_neigh,
               gpu_time_neigh_adjust,
               gpu_time_setup_comm,

               gpu_time_kdtree_setup,
               gpu_time_kdtree_query,
            gpu_time_rezone_all,
            gpu_time_count_BCs,
            gpu_time_calc_spatial_coordinates,
            gpu_time_load_balance;

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

   vector<int>local_indices_start;
   vector<int>local_indices_stop;

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
   vector<real>   lev_deltax,   //  Grid spacing along x-axis at specified level of refinement.
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
   real           xmin,         //  Lowest x-coordinate in use.
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

   vector<int>    index,        //  1D ordered index of mesh elements.
                  i,            //  1D ordered index of mesh element x-indices for k-D tree.
                  j,            //  1D ordered index of mesh element y-indices for k-D tree.
                  k,            //  1D ordered index of mesh element z-indices for k-D tree.
                  level,        //  1D ordered index of mesh element refinement levels.
                  nlft,         //  1D ordered index of mesh element left neighbors.
                  nrht,         //  1D ordered index of mesh element right neighbors.
                  nbot,         //  1D ordered index of mesh element bottom neighbors.
                  ntop,         //  1D ordered index of mesh element top neighbors.
                  nfrt,         //  1D ordered index of mesh element front neighbors.
                  nbak,         //  1D ordered index of mesh element back neighbors.
                  celltype;     // 1D ordered index of mesh element cell types (ghost or real).

   vector<real>   x,            //  1D ordered index of mesh element x-coordinates.
                  dx,           //  1D ordered index of mesh element x-coordinate spacings.
                  y,            //  1D ordered index of mesh element y-coordinates.
                  dy,           //  1D ordered index of mesh element y-coordinate spacings.
                  z,            //  1D ordered index of mesh element z-coordinates.
                  dz;           //  1D ordered index of mesh element z-coordinate spacings.

#ifdef HAVE_OPENCL
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
   Mesh(int nx, int ny, int levmx_in, int ndim_in, int numpe, int boundary, int parallel_in, int do_gpu_calc);

   //   Member functions.
#ifdef HAVE_OPENCL
   void init(int nx, int ny, double circ_radius, cl_context context, partition_method initial_order, int compute_device, int do_gpu_calc);
   void terminate(void);
#else
   void init(int nx, int ny, double circ_radius, partition_method initial_order, int do_gpu_calc);
#endif

   void resize_old_device_memory(size_t ncells);

   int  is_left_boundary(int ic)    { return (i[ic] < lev_ibegin[level[ic]]); }
   int  is_right_boundary(int ic)   { return (i[ic] > lev_iend[  level[ic]]); }
   int  is_bottom_boundary(int ic)  { return (j[ic] < lev_jbegin[level[ic]]); }
   int  is_top_boundary(int ic)     { return (j[ic] > lev_jend[  level[ic]]); }
   int  is_front_boundary(int ic)   { return (k[ic] < lev_kbegin[level[ic]]); }
   int  is_back_boundary(int ic)    { return (k[ic] > lev_kend[  level[ic]]); }

   int is_lower_left(int i, int j)  { return(i % 2 == 0 && j % 2 == 0); }
   int is_lower_right(int i, int j) { return(i % 2 == 1 && j % 2 == 0); }
   int is_upper_left(int i, int j)  { return(i % 2 == 0 && j % 2 == 1); }
   int is_upper_right(int i, int j) { return(i % 2 == 1 && j % 2 == 1); }

   double get_cpu_time_calc_neighbors(void)           {return(cpu_time_calc_neighbors); };
   double get_cpu_time_hash_setup(void)               {return(cpu_time_hash_setup); };
   double get_cpu_time_hash_query(void)               {return(cpu_time_hash_query); };
   double get_cpu_time_find_boundary(void)            {return(cpu_time_find_boundary); };
   double get_cpu_time_gather_boundary(void)          {return(cpu_time_gather_boundary); };
   double get_cpu_time_local_list(void)               {return(cpu_time_local_list); };
   double get_cpu_time_layer1(void)                   {return(cpu_time_layer1); };
   double get_cpu_time_layer2(void)                   {return(cpu_time_layer2); };
   double get_cpu_time_layer_list(void)               {return(cpu_time_layer_list); };
   double get_cpu_time_copy_mesh_data(void)           {return(cpu_time_copy_mesh_data); };
   double get_cpu_time_fill_mesh_ghost(void)          {return(cpu_time_fill_mesh_ghost); };
   double get_cpu_time_fill_neigh_ghost(void)         {return(cpu_time_fill_neigh_ghost); };
   double get_cpu_time_set_corner_neigh(void)         {return(cpu_time_set_corner_neigh); };
   double get_cpu_time_neigh_adjust(void)             {return(cpu_time_neigh_adjust); };
   double get_cpu_time_setup_comm(void)               {return(cpu_time_setup_comm); };
   double get_cpu_time_kdtree_setup(void)             {return(cpu_time_kdtree_setup); };
   double get_cpu_time_kdtree_query(void)             {return(cpu_time_kdtree_query); };
   double get_cpu_time_rezone_all(void)               {return(cpu_time_rezone_all); };
   double get_cpu_time_partition(void)                {return(cpu_time_partition); };
   double get_cpu_time_calc_spatial_coordinates(void) {return(cpu_time_calc_spatial_coordinates); };
   double get_cpu_time_load_balance(void)           {return(cpu_time_load_balance); };

   long get_gpu_time_calc_neighbors(void)           {return(gpu_time_calc_neighbors); };
   long get_gpu_time_hash_setup(void)               {return(gpu_time_hash_setup); };
   long get_gpu_time_hash_query(void)               {return(gpu_time_hash_query); };
   long get_gpu_time_find_boundary(void)            {return(gpu_time_find_boundary); };
   long get_gpu_time_gather_boundary(void)          {return(gpu_time_gather_boundary); };
   long get_gpu_time_local_list(void)               {return(gpu_time_local_list); };
   long get_gpu_time_layer1(void)                   {return(gpu_time_layer1); };
   long get_gpu_time_layer2(void)                   {return(gpu_time_layer2); };
   long get_gpu_time_layer_list(void)               {return(gpu_time_layer_list); };
   long get_gpu_time_copy_mesh_data(void)           {return(gpu_time_copy_mesh_data); };
   long get_gpu_time_fill_mesh_ghost(void)          {return(gpu_time_fill_mesh_ghost); };
   long get_gpu_time_fill_neigh_ghost(void)         {return(gpu_time_fill_neigh_ghost); };
   long get_gpu_time_set_corner_neigh(void)         {return(gpu_time_set_corner_neigh); };
   long get_gpu_time_neigh_adjust(void)             {return(gpu_time_neigh_adjust); };
   long get_gpu_time_setup_comm(void)               {return(gpu_time_setup_comm); };
   long get_gpu_time_kdtree_setup(void)             {return(gpu_time_kdtree_setup); };
   long get_gpu_time_kdtree_query(void)             {return(gpu_time_kdtree_query); };
   long get_gpu_time_rezone_all(void)               {return(gpu_time_rezone_all); };
   long get_gpu_time_calc_spatial_coordinates(void) {return(gpu_time_calc_spatial_coordinates); };
   long get_gpu_time_load_balance(void)             {return(gpu_time_load_balance); };

   int get_cpu_load_balance_count(void)           {return(cpu_load_balance_counter); };
   int get_cpu_rezone_count(void)                 {return(cpu_rezone_counter); };
   int get_cpu_refine_smooth_count(void)          {return(cpu_refine_smooth_counter); };
   int get_cpu_calc_neigh_count(void)             {return(cpu_calc_neigh_counter); };
   int get_gpu_load_balance_count(void)           {return(gpu_load_balance_counter); };
   int get_gpu_rezone_count(void)                 {return(gpu_rezone_counter); };
   int get_gpu_refine_smooth_count(void)          {return(gpu_refine_smooth_counter); };
   int get_gpu_calc_neigh_count(void)             {return(gpu_calc_neigh_counter); };

   void write_grid(int ncycle);
   void kdtree_setup(void);
   void calc_centerminmax(void);
   int  rezone_count(vector<int> mpot);
#ifdef HAVE_OPENCL
   void gpu_rezone_count(cl_command_queue command_queue, size_t block_size, size_t local_work_size, cl_mem dev_ioffset, cl_mem &dev_result);
#endif
   void print(void);
   void print_local(void);

   void partition_measure(void);
   void print_partition_measure(void);
   void print_calc_neighbor_type(void);
   int get_calc_neighbor_type(void);
   void print_partition_type(void);
   void partition_cells(int numpe,
                   vector<int> &order,
                   enum partition_method method);
   void calc_symmetry(vector<int> &dsym,
                  vector<int> &xsym,
                  vector<int> &ysym);

   /**************************************************************************************
   * Calculate neighbors
   **************************************************************************************/
   void calc_neighbors(void);
   void calc_neighbors_local(void);
#ifdef HAVE_OPENCL
   void gpu_calc_neighbors(cl_command_queue command_queue);
   void gpu_calc_neighbors_local(cl_command_queue command_queue);
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
   * Rezone mesh
   **************************************************************************************/
   void rezone_all(vector<int> mpot, int add_ncells);

   /**************************************************************************************
   * Load balance -- only needed for parallel (MPI) runs
   *  Input
   *    numcells -- ncells from rezone all routine. The load balance will reset
   *       the numcells argument and in the mesh object if a load balance is done.
   *    weight -- weighting array per cell for balancing. Currently not used. Null value
   *       indicates even weighting of cells for load balance. 
   *    state_memory or gpu_state_memory -- linked-list of arrays from physics routine
   *       to be load balanced. 
   * Output -- arrays will be returned load balanced with new sizes. Pointers to arrays
   *       will need to be reset
   **************************************************************************************/
#ifdef HAVE_MPI
   void do_load_balance_local(size_t &numcells, float *weight, MallocPlus &state_memory);
#ifdef HAVE_OPENCL
   int gpu_do_load_balance_local(cl_command_queue command_queue, size_t &numcells, float *weight, MallocPlus &gpu_state_memory);
#endif
#endif

   /**************************************************************************************
   * Calculate spatial coordinates
   **************************************************************************************/
   void calc_spatial_coordinates(int ibase);
#ifdef HAVE_OPENCL
   void gpu_calc_spatial_coordinates(cl_command_queue command_queue, cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy);
#endif

   void calc_distribution(int numpe);
#ifdef HAVE_OPENCL
   int gpu_count_BCs(cl_command_queue command_queue);
#endif

   void print_object_info();
   size_t refine_smooth(vector<int> &mpot);


   /**************************************************************************************
   * Testing routines
   **************************************************************************************/
#ifdef HAVE_OPENCL
   void compare_dev_local_to_local(cl_command_queue);
   void compare_neighbors_gpu_global_to_cpu_global(cl_command_queue command_queue);
#endif
   void compare_neighbors_cpu_local_to_cpu_global(uint ncells_ghost, uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl);
#ifdef HAVE_OPENCL
   void compare_neighbors_all_to_gpu_local(cl_command_queue command_queue, Mesh *mesh_global, int *nsizes, int *ndispl);
   void compare_mpot_gpu_global_to_cpu_global(cl_command_queue command_queue, int *mpot, cl_mem dev_mpot);
#endif
   void compare_mpot_cpu_local_to_cpu_global(uint ncells_global, int *nsizes, int *displ, int *mpot, int *mpot_global, int cycle);
#ifdef HAVE_OPENCL
   void compare_mpot_all_to_gpu_local(cl_command_queue command_queue, int *mpot, int *mpot_global, cl_mem dev_mpot, cl_mem dev_mpot_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle);
   void compare_ioffset_gpu_global_to_cpu_global(cl_command_queue command_queue, uint old_ncells, int *mpot, cl_mem dev_ioffset);
   void compare_ioffset_all_to_gpu_local(cl_command_queue command_queue, uint old_ncells, uint old_ncells_global, int block_size, int block_size_global, int *mpot, int *mpot_global, cl_mem dev_ioffset, cl_mem dev_ioffset_global, int *ioffset, int *ioffset_global, int *celltype_global);
   void compare_coordinates_gpu_global_to_cpu_global(cl_command_queue command_queue, cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, real *H);
#endif
   void compare_coordinates_cpu_local_to_cpu_global(uint ncells_global, int *nsizes, int *ndispl, real *x, real *dx, real *y, real *dy, real *H, real *x_global, real *dx_global, real *y_global, real *dy_global, real *H_global, int cycle);
#ifdef HAVE_OPENCL
   void compare_indices_gpu_global_to_cpu_global(cl_command_queue command_queue);
#endif
   void compare_indices_cpu_local_to_cpu_global(uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl, int cycle);
#ifdef HAVE_OPENCL
   void compare_indices_all_to_gpu_local(cl_command_queue command_queue, Mesh *mesh_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle);
#endif



private:
   //   Private constructors.
   Mesh(const Mesh&);   //   Blocks copy constructor so copies are not made inadvertently.

   //   Member functions.
   void calc_minmax(void);
   void calc_celltype(void);
#ifdef HAVE_OPENCL
   void print_dev_local(cl_command_queue);
#endif

/* Not currently called */
   void rezone_spread(vector<int> &mpot);
   void mesh_reorder(vector<int> iorder);
};

#endif /* MESH_H */
