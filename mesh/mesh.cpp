/*
 *  Copyright (c) 2011-2019, Triad National Security, LLC.
 *  All rights Reserved.
 *
 *  CLAMR -- LA-CC-11-094
 *
 *  Copyright 2011-2019. Triad National Security, LLC. This software was produced 
 *  under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Triad National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 *  TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
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
 *     * Neither the name of the Triad National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE TRIAD NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL
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
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <algorithm>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "hsfc/hsfc.h"
#include "kdtree/KDTree.h"
#include "mesh.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif
#include "timer/timer.h"
#ifdef HAVE_MPI
#include "l7/l7.h"
#endif
#include "reduce.h"
#include "genmalloc/genmalloc.h"
#include "hash/hash.h"

#define DEBUG 0
//#define BOUNDS_CHECK 1

#ifndef DEBUG
#define DEBUG 0
#endif
#undef DEBUG_RESTORE_VALS

typedef int scanInt;
void scan ( scanInt *input , scanInt *output , scanInt length);

#ifdef _OPENMP
#undef REZONE_NO_OPTIMIZATION
#else
#define REZONE_NO_OPTIMIZATION 1
#endif

#define TIMING_LEVEL 2

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#ifdef PATTERN_CHECK
#define MIN4(a,b,c,d) ( MIN(MIN(a,b),MIN(c,d)) )
#endif

#define IPOW2(a) (2 << (a))

#if defined(HALF_PRECISION)
#define CONSERVATION_EPS    1
#define STATE_EPS      30.0
#define HALF 0.5_h
#define ZERO 0.0_h

#elif defined(MINIMUM_PRECISION)
#define CONSERVATION_EPS    .1
#define STATE_EPS      15.0
#define HALF 0.5f
#define ZERO 0.0f

#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
#define CONSERVATION_EPS    .02
#define STATE_EPS        .025
#define HALF 0.5
#define ZERO 0.0

#elif defined(FULL_PRECISION)
#define CONSERVATION_EPS    .02
#define STATE_EPS        .025
#define HALF 0.5
#define ZERO 0.0

#endif

typedef unsigned int uint;
#ifdef __APPLE_CC__
typedef unsigned long ulong;
#endif

#define TWO 2

#define __NEW_STENCIL__
//#define __OLD_STENCIL__
//#define STENCIL_WARNING 1

#ifdef STENCIL_WARNING
int do_stencil_warning=1;
#else
int do_stencil_warning=0;
#endif

#ifdef HAVE_OPENCL
#include "mesh_kernel.inc"
#endif

extern bool localStencil;
int calc_neighbor_type;
bool dynamic_load_balance_on;
bool neighbor_remap;

#ifdef _OPENMP
//static bool iversion_flag = false;
#endif

int   choose_amr_method;

static const char *mesh_timer_descriptor[MESH_TIMER_SIZE] = {
   "mesh_timer_count_BCs",
   "mesh_timer_calc_neighbors",
   "mesh_timer_hash_setup",
   "mesh_timer_hash_query",
   "mesh_timer_find_boundary",
   "mesh_timer_push_setup",
   "mesh_timer_push_boundary",
   "mesh_timer_local_list",
   "mesh_timer_layer1",
   "mesh_timer_layer2",
   "mesh_timer_layer_list",
   "mesh_timer_copy_mesh_data",
   "mesh_timer_fill_mesh_ghost",
   "mesh_timer_fill_neigh_ghost",
   "mesh_timer_set_corner_neigh",
   "mesh_timer_neigh_adjust",
   "mesh_timer_setup_comm",
   "mesh_timer_kdtree_setup",
   "mesh_timer_kdtree_query",
   "mesh_timer_refine_smooth",
   "mesh_timer_rezone_all",
   "mesh_timer_partition",
   "mesh_timer_calc_spatial_coordinates",
   "mesh_timer_load_balance",
   "mesh_timer_bidir",
   "mesh_timer_bidirpart1",
   "mesh_timer_bidirpart2",
   "mesh_timer_bidirpart3",
   "mesh_timer_bidirpart4",
   "mesh_timer_bidirpart5",
   "mesh_timer_bidirpart6",
   "mesh_timer_bidirpart7",
   "mesh_timer_bidirpart8",
   "mesh_timer_bidirpart9",
   "mesh_timer_bidirpart10",
   "mesh_timer_bidirpart11",
   "mesh_timer_bidirpart12"
};

#ifdef HAVE_OPENCL
cl_kernel      kernel_hash_adjust_sizes;
cl_kernel      kernel_hash_setup;
cl_kernel      kernel_hash_setup_local;
cl_kernel      kernel_neighbor_init;
cl_kernel      kernel_calc_neighbors;
cl_kernel      kernel_calc_neighbors_local;
cl_kernel      kernel_calc_border_cells;
cl_kernel      kernel_calc_border_cells2;
cl_kernel      kernel_finish_scan;
cl_kernel      kernel_get_border_data;
cl_kernel      kernel_calc_layer1;
cl_kernel      kernel_calc_layer1_sethash;
cl_kernel      kernel_calc_layer2;
cl_kernel      kernel_get_border_data2;
cl_kernel      kernel_calc_layer2_sethash;
cl_kernel      kernel_copy_mesh_data;
cl_kernel      kernel_fill_mesh_ghost;
cl_kernel      kernel_fill_neighbor_ghost;
cl_kernel      kernel_set_corner_neighbor;
cl_kernel      kernel_adjust_neighbors_local;
cl_kernel      kernel_reduction_scan2;
cl_kernel      kernel_reduction_count;
cl_kernel      kernel_reduction_count2;
cl_kernel      kernel_hash_size;
cl_kernel      kernel_finish_hash_size;
cl_kernel      kernel_calc_spatial_coordinates;
cl_kernel      kernel_count_BCs;
cl_kernel      kernel_do_load_balance_lower;
cl_kernel      kernel_do_load_balance_middle;
cl_kernel      kernel_do_load_balance_upper;
cl_kernel      kernel_calc_face_list_wbidirmap_pt1;
cl_kernel      kernel_calc_face_list_wbidirmap_pt2;
cl_kernel      kernel_face_idx_wbidirmap;
cl_kernel      kernel_wbidirmap_precount;
cl_kernel      kernel_calc_wbidirmap_phantom_neighbors;;
cl_kernel      kernel_calc_wbidirmap_phantom_values;;
cl_kernel      kernel_deep_copy;
#ifndef MINIMUM_PRECISION
cl_kernel      kernel_do_load_balance_double;
#endif
cl_kernel      kernel_do_load_balance_float;
cl_kernel      kernel_refine_smooth;
cl_kernel      kernel_coarsen_smooth;
cl_kernel      kernel_coarsen_check_block;
cl_kernel      kernel_rezone_all;
cl_kernel      kernel_rezone_neighbors;
#ifndef MINIMUM_PRECISION
cl_kernel      kernel_rezone_one_double;
#endif
cl_kernel      kernel_rezone_one_float;
cl_kernel      kernel_copy_mpot_ghost_data;
cl_kernel      kernel_set_boundary_refinement;
#endif

extern size_t hash_header_size;
extern int   choose_hash_method;

void Mesh::write_grid(int ncycle)
{
   FILE *fp;
   char filename[20];

   if (ncycle<0) ncycle=0;
   sprintf(filename,"grid%02d.gph",ncycle);
   fp=fopen(filename,"w");

   fprintf(fp,"viewport %lf %lf %lf %lf\n",xmin,ymin,xmax,ymax);
   for (uint ic = 0; ic < ncells; ic++) {
      fprintf(fp,"rect  %lf   %lf   %lf   %lf\n",x[ic],y[ic],x[ic]+dx[ic],y[ic]+dy[ic]);
   }

   fprintf(fp,"line_init %lf %lf\n",x[0]+0.5*dx[0],y[0]+0.5*dy[0]);
   for (uint ic = 1; ic < ncells; ic++){
      fprintf(fp,"line %lf %lf\n",x[ic]+0.5*dx[ic],y[ic]+0.5*dy[ic]);
   }

   for (uint ic = 0; ic < ncells; ic++){
      fprintf(fp,"text %lf %lf %d\n",x[ic]+0.5*dx[ic],y[ic]+0.5*dy[ic],ic);
   }

   fclose(fp);
}

Mesh::Mesh(FILE *fin, int *numpe)
{
   char string[80];
   ibase = 1;

   time_t trand;
   time(&trand);
   srand48((long)trand);

   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"levmax %d",&levmx);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"cells %ld",&ncells);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"numpe %d",numpe);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"ndim %d",&ndim);
   if(fgets(string, 80, fin) == NULL) exit(-1);
#ifdef MINIMUM_PRECISION
   sscanf(string,"xaxis %f %f",&xmin, &deltax);
#else
   sscanf(string,"xaxis %lf %lf",&xmin, &deltax);
#endif
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"yaxis %lf %lf",(double*)&ymin, (double*)&deltay);
   if (ndim == THREE_DIMENSIONAL){
     if(fgets(string, 80, fin) == NULL) exit(-1);
     sscanf(string,"zaxis %lf %lf",(double*)&zmin, (double*)&deltaz);
   }
   if(fgets(string, 80, fin) == NULL) exit(-1);

   index.resize(ncells);

   allocate(ncells);

   uint ic=0;
   while(fgets(string, 80, fin)!=NULL){
      sscanf(string, "%d %d %d %d", &(index[ic]), &(i[ic]), &(j[ic]), &(level[ic]));
      ic++;
   }

   ibase = 0;
   calc_spatial_coordinates(ibase);
   KDTree_Initialize(&tree);


  print();

   if (ic != ncells) {
      printf("Error -- cells read does not match number specified\n");
   }
   return;
}

void Mesh::print(void)
{
   assert(&nlft[0] != NULL);
   assert(&x[0] != NULL);
   assert(&index[0] != NULL);

   //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());
   printf("index orig index   i     j     lev   nlft  nrht  nbot  ntop   xlow    xhigh     ylow    yhigh\n");
   for (uint ic=0; ic<ncells; ic++)
   {  printf("%6d %6d   %4d  %4d   %4d  %4d  %4d  %4d  %4d ", ic, index[ic], i[ic], j[ic], level[ic], nlft[ic], nrht[ic], nbot[ic], ntop[ic]);
      printf("%8.2lf %8.2lf %8.2lf %8.2lf\n", x[ic], x[ic]+dx[ic], y[ic], y[ic]+dy[ic]); }
}

void Mesh::print_local()
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   if (mesh_memory.get_memory_size(nlft) >= ncells_ghost){
      fprintf(fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mype);
      for (uint ic=0; ic<ncells; ic++) {
         fprintf(fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mype,ic, ic+noffset,i[ic], j[ic], level[ic], nlft[ic], nrht[ic], nbot[ic], ntop[ic]);
      }
      for (uint ic=ncells; ic<ncells_ghost; ic++) {
         fprintf(fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mype,ic, ic+noffset,i[ic], j[ic], level[ic], nlft[ic], nrht[ic], nbot[ic], ntop[ic]);
      }
   } else {
      fprintf(fp,"%d:    index   i     j     lev\n",mype);
      for (uint ic=0; ic<ncells_ghost; ic++) {
         fprintf(fp,"%d: %6d  %4d  %4d   %4d  \n", mype,ic, i[ic], j[ic], level[ic]);
      }
   }
}

#ifdef HAVE_OPENCL
void Mesh::print_dev_local(void)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<int>i_tmp(ncells_ghost);
   vector<int>j_tmp(ncells_ghost);
   vector<uchar_t>level_tmp(ncells_ghost);
   vector<int>nlft_tmp(ncells_ghost);
   vector<int>nrht_tmp(ncells_ghost);
   vector<int>nbot_tmp(ncells_ghost);
   vector<int>ntop_tmp(ncells_ghost);
   ezcl_enqueue_read_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &i_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &j_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level, CL_FALSE, 0, ncells_ghost*sizeof(cl_uchar_t), &level_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0], NULL);

   //fprintf(fp,"\n%d:                    Printing mesh for dev_local\n\n",mype);

   fprintf(fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mype);
   for (uint ic=0; ic<MAX(ncells_ghost,ncells); ic++) {
      fprintf(fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mype,ic, ic+noffset,i_tmp[ic], j_tmp[ic], level_tmp[ic], nlft_tmp[ic], nrht_tmp[ic], nbot_tmp[ic], ntop_tmp[ic]);
   }
   //fprintf(fp,"\n%d:              Finished printing mesh for dev_local\n\n",mype);
}

void Mesh::compare_dev_local_to_local(void)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<int>i_tmp(ncells_ghost);
   vector<int>j_tmp(ncells_ghost);
   vector<uchar_t>level_tmp(ncells_ghost);
   vector<int>nlft_tmp(ncells_ghost);
   vector<int>nrht_tmp(ncells_ghost);
   vector<int>nbot_tmp(ncells_ghost);
   vector<int>ntop_tmp(ncells_ghost);
   ezcl_enqueue_read_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &i_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &j_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, ncells_ghost*sizeof(cl_uchar_t), &level_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0], NULL);

   fprintf(fp,"\n%d:                      Comparing mesh for dev_local to local\n\n",mype);
   //fprintf(fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mype);
   for (uint ic=0; ic<ncells_ghost; ic++) {
      if (i_tmp[ic]     != i[ic]    ) fprintf(fp,"%d: Error: cell %d dev_i     %d i     %d\n",mype,ic,i_tmp[ic],    i[ic]);
      if (j_tmp[ic]     != j[ic]    ) fprintf(fp,"%d: Error: cell %d dev_j     %d j     %d\n",mype,ic,j_tmp[ic],    j[ic]);
      if (level_tmp[ic] != level[ic]) fprintf(fp,"%d: Error: cell %d dev_level %d level %d\n",mype,ic,level_tmp[ic],level[ic]);

      //fprintf(fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mype,ic, ic+noffset,i_tmp[ic], j_tmp[ic], level_tmp[ic], nlft_tmp[ic], nrht_tmp[ic], nbot_tmp[ic], ntop_tmp[ic]);
   }
   fprintf(fp,"\n%d:                 Finished comparing mesh for dev_local to local\n\n",mype);
}

void Mesh::compare_neighbors_gpu_global_to_cpu_global()
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<int>nlft_check(ncells);
   vector<int>nrht_check(ncells);
   vector<int>nbot_check(ncells);
   vector<int>ntop_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, ncells*sizeof(cl_int), &nlft_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, ncells*sizeof(cl_int), &nrht_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, ncells*sizeof(cl_int), &nbot_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_TRUE,  0, ncells*sizeof(cl_int), &ntop_check[0], NULL);

   //printf("\n%d:                      Comparing neighbors for gpu_global to cpu_global\n\n",mype);
   for (uint ic=0; ic<ncells; ic++) {
      if (nlft[ic] != nlft_check[ic]) printf("DEBUG -- nlft: ic %d nlft %d nlft_check %d\n",ic, nlft[ic], nlft_check[ic]);
      if (nrht[ic] != nrht_check[ic]) printf("DEBUG -- nrht: ic %d nrht %d nrht_check %d\n",ic, nrht[ic], nrht_check[ic]);
      if (nbot[ic] != nbot_check[ic]) printf("DEBUG -- nbot: ic %d nbot %d nbot_check %d\n",ic, nbot[ic], nbot_check[ic]);
      if (ntop[ic] != ntop_check[ic]) printf("DEBUG -- ntop: ic %d ntop %d ntop_check %d\n",ic, ntop[ic], ntop_check[ic]);
   }
   //printf("\n%d:                 Finished comparing mesh for dev_local to local\n\n",mype);
}
#endif

void Mesh::compare_neighbors_cpu_local_to_cpu_global(uint ncells_ghost, uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl)
{

#ifdef HAVE_MPI
   int *nlft_global = mesh_global->nlft;
   int *nrht_global = mesh_global->nrht;
   int *nbot_global = mesh_global->nbot;
   int *ntop_global = mesh_global->ntop;

   vector<int> Test(ncells_ghost);
   for(uint ic=0; ic<ncells; ic++){
      Test[ic] = mype*1000 +ic;
   }
   if (numpe > 1) L7_Update(&Test[0], L7_INT, cell_handle);

   vector<int> Test_global(ncells_global);
   MPI_Allgatherv(&Test[0], nsizes[mype], MPI_INT, &Test_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

   vector<int> Test_check(ncells);
   vector<int> Test_check_global(ncells_global);

   // ==================== check left value ====================
   for (uint ic=0; ic<ncells; ic++){
      Test_check[ic] = Test[nlft[ic]];
   }

   MPI_Allgatherv(&Test_check[0], nsizes[mype], MPI_INT, &Test_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);

   for (uint ic=0; ic<ncells_global; ic++){
      if (Test_global[nlft_global[ic]] != Test_check_global[ic]) {
         if (mype == 0) printf("%d: Error with nlft for cell %d -- nlft %d global %d check %d\n",mype,ic,nlft_global[ic],Test_global[nlft_global[ic]],Test_check_global[ic]);
      }
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
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- ncells_global %d ncells_ghost %d mesh_global %p nsizes[0] %d ndispl[0] %d\n",
               ncells_global,ncells_ghost,mesh_global,nsizes[0],ndispl[0]);
#endif
}

#ifdef HAVE_OPENCL
void Mesh::compare_neighbors_all_to_gpu_local(Mesh *mesh_global, int *nsizes, int *ndispl)
//uint ncells_ghost, uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl)
{
#ifdef HAVE_MPI
   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t &ncells_global = mesh_global->ncells;
   int *nlft_global = mesh_global->nlft;
   int *nrht_global = mesh_global->nrht;
   int *nbot_global = mesh_global->nbot;
   int *ntop_global = mesh_global->ntop;

   // Checking CPU parallel to CPU global
   vector<int> Test(ncells_ghost);
   for(uint ic=0; ic<ncells; ic++){
      Test[ic] = mype*1000 +ic; 
   }    
   if (numpe > 1) L7_Update(&Test[0], L7_INT, cell_handle);

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
         printf("%d:                           check %5d -- nlftl %5d nlftl nlftl %5d check  %5d\n",
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
   // checking gpu results
   vector<int> nlft_check(ncells_ghost);         vector<int> nrht_check(ncells_ghost);
   vector<int> nbot_check(ncells_ghost);         vector<int> ntop_check(ncells_ghost);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nlft_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nrht_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nbot_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells_ghost*sizeof(cl_int),  &ntop_check[0], NULL);

   for (uint ic=0; ic<ncells_ghost; ic++){
      if (nlft[ic] != nlft_check[ic]) printf("%d: Error with gpu calculated nlft for cell %d nlft %d check %d\n",mype,ic,nlft[ic],nlft_check[ic]);
      if (nrht[ic] != nrht_check[ic]) printf("%d: Error with gpu calculated nrht for cell %d nrht %d check %d\n",mype,ic,nrht[ic],nrht_check[ic]);
      if (nbot[ic] != nbot_check[ic]) printf("%d: Error with gpu calculated nbot for cell %d nbot %d check %d\n",mype,ic,nbot[ic],nbot_check[ic]);
      if (ntop[ic] != ntop_check[ic]) printf("%d: Error with gpu calculated ntop for cell %d ntop %d check %d\n",mype,ic,ntop[ic],ntop_check[ic]);
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
   // checking gpu results
   //vector<int> nlft_check(ncells_ghost);         vector<int> nrht_check(ncells_ghost);
   //vector<int> nbot_check(ncells_ghost);         vector<int> ntop_check(ncells_ghost);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nlft_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nrht_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells_ghost*sizeof(cl_int),  &nbot_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells_ghost*sizeof(cl_int),  &ntop_check[0], NULL);

   for (uint ic=0; ic<ncells_ghost; ic++){
      if (nlft[ic] != nlft_check[ic]) printf("%d: Error with gpu calculated nlft for cell %d nlft %d check %d\n",mype,ic,nlft[ic],nlft_check[ic]);
      if (nrht[ic] != nrht_check[ic]) printf("%d: Error with gpu calculated nrht for cell %d nrht %d check %d\n",mype,ic,nrht[ic],nrht_check[ic]);
      if (nbot[ic] != nbot_check[ic]) printf("%d: Error with gpu calculated nbot for cell %d nbot %d check %d\n",mype,ic,nbot[ic],nbot_check[ic]);
      if (ntop[ic] != ntop_check[ic]) printf("%d: Error with gpu calculated ntop for cell %d ntop %d check %d\n",mype,ic,ntop[ic],ntop_check[ic]);
   }
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- mesh_global %p nsizes[0] %d ndispl[0] %d\n",
               mesh_global,nsizes[0],ndispl[0]);
#endif
}

void Mesh::compare_indices_gpu_global_to_cpu_global(void)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<int> i_check(ncells);
   vector<int> j_check(ncells);
   vector<uchar_t> level_check(ncells);
   vector<char_t> celltype_check(ncells);
   /// Set read buffers for data.
   ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int), &i_check[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int), &j_check[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_uchar_t), &level_check[0],    NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  0, ncells*sizeof(cl_char_t), &celltype_check[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (i[ic]        != i_check[ic] )        printf("DEBUG -- i: ic %d i %d i_check %d\n",ic, i[ic], i_check[ic]);
      if (j[ic]        != j_check[ic] )        printf("DEBUG -- j: ic %d j %d j_check %d\n",ic, j[ic], j_check[ic]);
      if (level[ic]    != level_check[ic] )    printf("DEBUG -- level: ic %d level %d level_check %d\n",ic, level[ic], level_check[ic]);
      if (celltype[ic] != celltype_check[ic] ) printf("DEBUG -- celltype: ic %d celltype %d celltype_check %d\n",ic, celltype[ic], celltype_check[ic]);
   }
}
#endif

void Mesh::compare_indices_cpu_local_to_cpu_global(uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl, int cycle)
{
   char_t *celltype_global = mesh_global->celltype;
   int *i_global        = mesh_global->i;
   int *j_global        = mesh_global->j;
   uchar_t *level_global    = mesh_global->level;

   vector<int> i_check_global(ncells_global);
   vector<int> j_check_global(ncells_global);
   vector<uchar_t> level_check_global(ncells_global);
   vector<char_t> celltype_check_global(ncells_global);

/*
   vector<int> i_check_local(ncells);
   vector<int> j_check_local(ncells);
   vector<uchar_t> level_check_local(ncells);
   vector<char_t> celltype_check_local(ncells);
*/

#ifdef HAVE_MPI
   MPI_Allgatherv(&celltype[0], nsizes[mype], MPI_CHAR_T, &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&i[0],        nsizes[mype], MPI_INT, &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&j[0],        nsizes[mype], MPI_INT, &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&level[0],    nsizes[mype], MPI_UCHAR_T, &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_UCHAR_T, MPI_COMM_WORLD);

/*
   MPI_Scatterv(&celltype_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T, &celltype_check_local[0], nsizes[mype], MPI_CHAR_T, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&i_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &i_check_local[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&j_global[0],        &nsizes[0], &ndispl[0], MPI_INT, &j_check_local[0],        nsizes[mype], MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&level_global[0],    &nsizes[0], &ndispl[0], MPI_UCHAR_T, &level_check_local[0],    nsizes[mype], MPI_UCHAR_T, 0, MPI_COMM_WORLD);
*/
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- nsizes[0] %d ndispl[0] %d\n",
               nsizes[0],ndispl[0]);
#endif

   for (uint ic = 0; ic < ncells_global; ic++){
      if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 3 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",cycle,ic,celltype_global[ic],celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 3 at cycle %d i_global & i_check_global %d %d  %d  \n",cycle,ic,i_global[ic],i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 3 at cycle %d j_global & j_check_global %d %d  %d  \n",cycle,ic,j_global[ic],j_check_global[ic]);
      if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 3 at cycle %d level_global & level_check_global %d %d  %d  \n",cycle,ic,level_global[ic],level_check_global[ic]);
   }

/*
   for (uint ic = 0; ic < ncells; ic++){
      if (celltype[ic] != celltype_check_local[ic])  fprintf(fp,"DEBUG rezone 3 at cycle %d celltype & celltype_check_local %d %d  %d  \n",cycle,ic,celltype[ic],celltype_check_local[ic]);
      if (i[ic] != i_check_local[ic])                fprintf(fp,"DEBUG rezone 3 at cycle %d i & i_check_local %d %d  %d  \n",cycle,ic,i[ic],i_check_local[ic]);
      if (j[ic] != j_check_local[ic])                fprintf(fp,"DEBUG rezone 3 at cycle %d j & j_check_local %d %d  %d  \n",cycle,ic,j[ic],j_check_local[ic]);
      if (level[ic] != level_check_local[ic])        fprintf(fp,"DEBUG rezone 3 at cycle %d level & level_check_local %d %d  %d  \n",cycle,ic,level[ic],level_check_local[ic]);
   }
*/
}

#ifdef HAVE_OPENCL
void Mesh::compare_indices_all_to_gpu_local(Mesh *mesh_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle)
{
#ifdef HAVE_MPI
   cl_command_queue command_queue = ezcl_get_command_queue();

   uchar_t *level_global = mesh_global->level;
   char_t *celltype_global = mesh_global->celltype;
   int *i_global = mesh_global->i;
   int *j_global = mesh_global->j;

   cl_mem &dev_celltype_global = mesh_global->dev_celltype;
   cl_mem &dev_i_global = mesh_global->dev_i;
   cl_mem &dev_j_global = mesh_global->dev_j;
   cl_mem &dev_level_global = mesh_global->dev_level;

   // Need to compare dev_H to H, etc
   vector<uchar_t> level_check(ncells);
   vector<char_t> celltype_check(ncells);
   vector<int> i_check(ncells);
   vector<int> j_check(ncells);
   /// Set read buffers for data.
   ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_uchar_t),  &level_check[0],     NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_char_t),  &celltype_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  &i_check[0],         NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_TRUE,  0, ncells*sizeof(cl_int),  &j_check[0],         NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (level[ic] != level_check[ic] )       printf("%d: DEBUG rezone 1 cell %d level %d level_check %d\n",mype, ic, level[ic], level_check[ic]);
      if (celltype[ic] != celltype_check[ic] ) printf("%d: DEBUG rezone 1 cell %d celltype %d celltype_check %d\n",mype, ic, celltype[ic], celltype_check[ic]);
      if (i[ic] != i_check[ic] )               printf("%d: DEBUG rezone 1 cell %d i %d i_check %d\n",mype, ic, i[ic], i_check[ic]);
      if (j[ic] != j_check[ic] )               printf("%d: DEBUG rezone 1 cell %d j %d j_check %d\n",mype, ic, j[ic], j_check[ic]);
   }

   // And compare dev_H gathered to H_global, etc
   vector<char_t>celltype_check_global(ncells_global);
   vector<int>i_check_global(ncells_global);
   vector<int>j_check_global(ncells_global);
   vector<uchar_t>level_check_global(ncells_global);
   MPI_Allgatherv(&celltype_check[0], nsizes[mype], MPI_CHAR_T,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T,    MPI_COMM_WORLD);
   MPI_Allgatherv(&i_check[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&j_check[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&level_check[0],    nsizes[mype], MPI_UCHAR_T,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_UCHAR_T,    MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (level_global[ic] != level_check_global[ic] )       printf("%d: DEBUG rezone 2 cell %d level_global %d level_check_global %d\n",mype, ic, level_global[ic], level_check_global[ic]);
      if (celltype_global[ic] != celltype_check_global[ic] ) printf("%d: DEBUG rezone 2 cell %d celltype_global %d celltype_check_global %d\n",mype, ic, celltype_global[ic], celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic] )               printf("%d: DEBUG rezone 2 cell %d i_global %d i_check_global %d\n",mype, ic, i_global[ic], i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic] )               printf("%d: DEBUG rezone 2 cell %d j_global %d j_check_global %d\n",mype, ic, j_global[ic], j_check_global[ic]);
   }

   // And compare H gathered to H_global, etc
   MPI_Allgatherv(&celltype[0], nsizes[mype], MPI_CHAR_T,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T,    MPI_COMM_WORLD);
   MPI_Allgatherv(&i[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&j[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&level[0],    nsizes[mype], MPI_UCHAR_T,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_UCHAR_T,    MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 3 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",ncycle,ic,celltype_global[ic],celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 3 at cycle %d i_global & i_check_global %d %d  %d  \n",ncycle,ic,i_global[ic],i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 3 at cycle %d j_global & j_check_global %d %d  %d  \n",ncycle,ic,j_global[ic],j_check_global[ic]);
      if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 3 at cycle %d level_global & level_check_global %d %d  %d  \n",ncycle,ic,level_global[ic],level_check_global[ic]);
   }

   // Now the global dev_H_global to H_global, etc
   ezcl_enqueue_read_buffer(command_queue, dev_celltype_global, CL_FALSE, 0, ncells_global*sizeof(cl_char_t),  &celltype_check_global[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_i_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  &i_check_global[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  &j_check_global[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level_global,    CL_TRUE,  0, ncells_global*sizeof(cl_uchar_t),  &level_check_global[0],    NULL);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 4 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",ncycle,ic,celltype_global[ic],celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 4 at cycle %d i_global & i_check_global %d %d  %d  \n",ncycle,ic,i_global[ic],i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 4 at cycle %d j_global & j_check_global %d %d  %d  \n",ncycle,ic,j_global[ic],j_check_global[ic]);
      if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 4 at cycle %d level_global & level_check_global %d %d  %d  \n",ncycle,ic,level_global[ic],level_check_global[ic]);
   }
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- mesh_global %p ncells_global %d nsizes[0] %d ndispl[0] %d ncycle %d\n",
               mesh_global,ncells_global,nsizes[0],ndispl[0],ncycle);
#endif
}

void Mesh::compare_coordinates_gpu_global_to_cpu_global_double(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, double *H)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<spatial_t>x_check(ncells);
   vector<spatial_t>dx_check(ncells);
   vector<spatial_t>y_check(ncells);
   vector<spatial_t>dy_check(ncells);
   vector<double>H_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_x,   CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &x_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dx,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &dx_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,   CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &y_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &dy_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_double), &H_check[0],  NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (x[ic] != x_check[ic] || dx[ic] != dx_check[ic] || y[ic] != y_check[ic] || dy[ic] != dy_check[ic] ) {
         printf("Error -- mismatch in spatial coordinates for cell %d is gpu %lf %lf %lf %lf cpu %lf %lf %lf %lf\n",ic,x_check[ic],dx_check[ic],y_check[ic],dy_check[ic],x[ic],dx[ic],y[ic],dy[ic]);
         exit(0);
      }
   }  
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic] - H_check[ic]) > CONSERVATION_EPS) {
         printf("Error -- mismatch in H for cell %d is gpu %lf cpu %lf\n",ic,H_check[ic],H[ic]);
         exit(0);
      }
   }
}

void Mesh::compare_coordinates_gpu_global_to_cpu_global_float(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, float *H)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<spatial_t>x_check(ncells);
   vector<spatial_t>dx_check(ncells);
   vector<spatial_t>y_check(ncells);
   vector<spatial_t>dy_check(ncells);
   vector<float>H_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_x,   CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &x_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dx,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &dx_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,   CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &y_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &dy_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_float), &H_check[0],  NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (x[ic] != x_check[ic] || dx[ic] != dx_check[ic] || y[ic] != y_check[ic] || dy[ic] != dy_check[ic] ) {
         printf("Error -- mismatch in spatial coordinates for cell %d is gpu %lf %lf %lf %lf cpu %lf %lf %lf %lf\n",ic,x_check[ic],dx_check[ic],y_check[ic],dy_check[ic],x[ic],dx[ic],y[ic],dy[ic]);
         exit(0);
      }
   }  
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic] - H_check[ic]) > CONSERVATION_EPS) {
         printf("Error -- mismatch in H for cell %d is gpu %lf cpu %lf\n",ic,H_check[ic],H[ic]);
         exit(0);
      }
   }
}

#ifdef HALF_PRECISION
void Mesh::compare_coordinates_gpu_global_to_cpu_global_half(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, half *H)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<spatial_t>x_check(ncells);
   vector<spatial_t>dx_check(ncells);
   vector<spatial_t>y_check(ncells);
   vector<spatial_t>dy_check(ncells);
   vector<half>H_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_x,   CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &x_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dx,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &dx_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,   CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &y_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy,  CL_FALSE, 0, ncells*sizeof(cl_spatial_t), &dy_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_half), &H_check[0],  NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (x[ic] != x_check[ic] || dx[ic] != dx_check[ic] || y[ic] != y_check[ic] || dy[ic] != dy_check[ic] ) {
         printf("Error -- mismatch in spatial coordinates for cell %d is gpu %lf %lf %lf %lf cpu %lf %lf %lf %lf\n",ic,x_check[ic],dx_check[ic],y_check[ic],dy_check[ic],x[ic],dx[ic],y[ic],dy[ic]);
         exit(0);
      }
   }  
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic] - H_check[ic]) > CONSERVATION_EPS) {
         printf("Error -- mismatch in H for cell %d is gpu %lf cpu %lf\n",ic,H_check[ic],H[ic]);
         exit(0);
      }
   }
}
#endif

#endif

void Mesh::compare_coordinates_cpu_local_to_cpu_global_double(uint ncells_global, int *nsizes, int *ndispl, spatial_t *x, spatial_t *dx, spatial_t *y, spatial_t *dy, double *H, spatial_t *x_global, spatial_t *dx_global, spatial_t *y_global, spatial_t *dy_global, double *H_global, int cycle)
{
   vector<spatial_t> x_check_global(ncells_global);
   vector<spatial_t> dx_check_global(ncells_global);
   vector<spatial_t> y_check_global(ncells_global);
   vector<spatial_t> dy_check_global(ncells_global);
   vector<double> H_check_global(ncells_global);

#ifdef HAVE_MPI
   MPI_Allgatherv(&x[0],  nsizes[mype], MPI_SPATIAL_T, &x_check_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&dx[0], nsizes[mype], MPI_SPATIAL_T, &dx_check_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&y[0],  nsizes[mype], MPI_SPATIAL_T, &y_check_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&dy[0], nsizes[mype], MPI_SPATIAL_T, &dy_check_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&H[0],  nsizes[mype], MPI_DOUBLE, &H_check_global[0],  &nsizes[0], &ndispl[0], MPI_DOUBLE, MPI_COMM_WORLD);
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- nsizes[0] %d ndispl[0] %d x %p dx %p y %p dy %p H %p\n",
               nsizes[0],ndispl[0],x,dx,y,dy,H);
#endif

   for (uint ic = 0; ic < ncells_global; ic++){
      if (fabs(x_global[ic] -x_check_global[ic] ) > STATE_EPS) printf("DEBUG graphics at cycle %d x_global & x_check_global  %d %lf %lf \n",cycle,ic,x_global[ic], x_check_global[ic]);
      if (fabs(dx_global[ic]-dx_check_global[ic]) > STATE_EPS) printf("DEBUG graphics at cycle %d dx_global & dx_check_global %d %lf %lf \n",cycle,ic,dx_global[ic],dx_check_global[ic]);
      if (fabs(y_global[ic] -y_check_global[ic] ) > STATE_EPS) printf("DEBUG graphics at cycle %d y_global & y_check_global  %d %lf %lf \n",cycle,ic,y_global[ic], y_check_global[ic]);
      if (fabs(dy_global[ic]-dy_check_global[ic]) > STATE_EPS) printf("DEBUG graphics at cycle %d dy_global & dy_check_global %d %lf %lf \n",cycle,ic,dy_global[ic],dy_check_global[ic]);
      if (fabs(H_global[ic] -H_check_global[ic] ) > STATE_EPS) printf("DEBUG graphics at cycle %d H_global & H_check_global  %d %lf %lf \n",cycle,ic,H_global[ic], H_check_global[ic]);
   }

}

void Mesh::compare_coordinates_cpu_local_to_cpu_global_float(uint ncells_global, int *nsizes, int *ndispl, spatial_t *x, spatial_t *dx, spatial_t *y, spatial_t *dy, float *H, spatial_t *x_global, spatial_t *dx_global, spatial_t *y_global, spatial_t *dy_global, float *H_global, int cycle)
{
   vector<spatial_t> x_check_global(ncells_global);
   vector<spatial_t> dx_check_global(ncells_global);
   vector<spatial_t> y_check_global(ncells_global);
   vector<spatial_t> dy_check_global(ncells_global);
   vector<float> H_check_global(ncells_global);

#ifdef HAVE_MPI
   MPI_Allgatherv(&x[0],  nsizes[mype], MPI_SPATIAL_T, &x_check_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&dx[0], nsizes[mype], MPI_SPATIAL_T, &dx_check_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&y[0],  nsizes[mype], MPI_SPATIAL_T, &y_check_global[0],  &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&dy[0], nsizes[mype], MPI_SPATIAL_T, &dy_check_global[0], &nsizes[0], &ndispl[0], MPI_SPATIAL_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&H[0],  nsizes[mype], MPI_FLOAT,     &H_check_global[0],  &nsizes[0], &ndispl[0], MPI_FLOAT,     MPI_COMM_WORLD);
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- nsizes[0] %d ndispl[0] %d x %p dx %p y %p dy %p H %p\n",
               nsizes[0],ndispl[0],x,dx,y,dy,H);
#endif

   for (uint ic = 0; ic < ncells_global; ic++){
      if (fabs(x_global[ic] -x_check_global[ic] ) > STATE_EPS) printf("DEBUG graphics at cycle %d x_global & x_check_global  %d %lf %lf \n",cycle,ic,x_global[ic], x_check_global[ic]);
      if (fabs(dx_global[ic]-dx_check_global[ic]) > STATE_EPS) printf("DEBUG graphics at cycle %d dx_global & dx_check_global %d %lf %lf \n",cycle,ic,dx_global[ic],dx_check_global[ic]);
      if (fabs(y_global[ic] -y_check_global[ic] ) > STATE_EPS) printf("DEBUG graphics at cycle %d y_global & y_check_global  %d %lf %lf \n",cycle,ic,y_global[ic], y_check_global[ic]);
      if (fabs(dy_global[ic]-dy_check_global[ic]) > STATE_EPS) printf("DEBUG graphics at cycle %d dy_global & dy_check_global %d %lf %lf \n",cycle,ic,dy_global[ic],dy_check_global[ic]);
      if (fabs(H_global[ic] -H_check_global[ic] ) > STATE_EPS) printf("DEBUG graphics at cycle %d H_global & H_check_global  %d %lf %lf \n",cycle,ic,H_global[ic], H_check_global[ic]);
   }

}

#ifdef HAVE_OPENCL
void Mesh::compare_mpot_gpu_global_to_cpu_global(char_t *mpot, cl_mem dev_mpot)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<char_t>mpot_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_mpot,  CL_TRUE,  0, ncells*sizeof(cl_char_t), &mpot_check[0], NULL);

   for (uint ic=0; ic<ncells; ic++) {
      if (mpot[ic] != mpot_check[ic]) printf("DEBUG -- mpot: ic %d mpot %d mpot_check %d\n",ic, mpot[ic], mpot_check[ic]);
   }
}
#endif

void Mesh::compare_mpot_cpu_local_to_cpu_global(uint ncells_global, int *nsizes, int *ndispl, char_t *mpot, char_t *mpot_global, int cycle)
{
   vector<char_t>mpot_save_global(ncells_global);
#ifdef HAVE_MPI
   MPI_Allgatherv(&mpot[0], ncells, MPI_CHAR_T, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T, MPI_COMM_WORLD);
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- nsizes[0] %d ndispl[0] %d mpot %p\n",
               nsizes[0],ndispl[0],mpot);
#endif
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 3 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,cycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }
   }

}

#ifdef HAVE_OPENCL
void Mesh::compare_mpot_all_to_gpu_local(char_t *mpot, char_t *mpot_global, cl_mem dev_mpot, cl_mem dev_mpot_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle)
{
#ifdef HAVE_MPI
   cl_command_queue command_queue = ezcl_get_command_queue();

   // Need to compare dev_mpot to mpot 
   vector<char_t>mpot_save(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_char_t), &mpot_save[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (mpot[ic] != mpot_save[ic]) {
         printf("%d: DEBUG refine_potential 1 at cycle %d cell %d mpot & mpot_save %d %d \n",mype,ncycle,ic,mpot[ic],mpot_save[ic]);
      }    
   }    

   // Compare dev_mpot to mpot_global
   vector<char_t>mpot_save_global(ncells_global);
   MPI_Allgatherv(&mpot_save[0], nsizes[mype], MPI_CHAR_T, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T, MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 2 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,ncycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }    
   }    

   // Compare mpot to mpot_global
   MPI_Allgatherv(&mpot[0], nsizes[mype], MPI_CHAR_T, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_CHAR_T, MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 3 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,ncycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }    
   }    

   // Compare dev_mpot_global to mpot_global
   ezcl_enqueue_read_buffer(command_queue, dev_mpot_global, CL_TRUE,  0, ncells_global*sizeof(cl_char_t), &mpot_save_global[0], NULL);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 4 at cycle %d cell %u mpot_global & mpot_save_global %d %d \n",mype,ncycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }    
   }    
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- mpot %p mpot_global %p dev_mpot %p dev_mpot_global %p ncells_global %d nsizes[0] %d ndispl[0] %d ncycle %d\n",
               mpot,mpot_global,dev_mpot,dev_mpot_global,ncells_global,nsizes[0],ndispl[0],ncycle);
#endif
}

void Mesh::compare_ioffset_gpu_global_to_cpu_global(uint old_ncells, char_t *mpot)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t local_work_size  = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

   //size_t block_size = (ncells + TILE_SIZE - 1) / TILE_SIZE; //  For on-device global reduction kernel.
   size_t block_size     = global_work_size/local_work_size;

   vector<int> ioffset_check(block_size);
   ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int), &ioffset_check[0], NULL);

   int mcount, mtotal;
   mtotal = 0;
   for (uint ig=0; ig<(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
      mcount = 0;
      for (uint ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
         if (ic >= old_ncells) break;

         if (mpot[ic] < 0) {
            if (celltype[ic] == REAL_CELL) {
               // remove all but cell that will remain to get count right when split
               // across processors
               if (is_lower_left(i[ic],j[ic]) ) mcount++;
            } else {
               // either upper right or lower left will remain for boundary cells
               if (is_upper_right(i[ic],j[ic]) || is_lower_left(i[ic],j[ic]) ) mcount++;
            }
         }
         if (mpot[ic] >= 0) {
            if (celltype[ic] == REAL_CELL){
               mcount += mpot[ic] ? 4 : 1;
            } else {
               mcount += mpot[ic] ? 2 : 1;
            }
         }
      }
      if (mtotal != ioffset_check[ig]) printf("DEBUG ig %d ioffset %d mcount %d\n",ig,ioffset_check[ig],mtotal);
      mtotal += mcount;
   }
}

void Mesh::compare_ioffset_all_to_gpu_local(uint old_ncells, uint old_ncells_global, int block_size, int block_size_global, char_t *mpot, char_t *mpot_global, cl_mem dev_ioffset, cl_mem dev_ioffset_global, int *ioffset, int *ioffset_global, char_t *celltype_global, int *i_global, int *j_global)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   // This compares ioffset for each block in the calculation
   ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int), &ioffset[0], NULL);
   int mtotal = 0; 
   for (uint ig=0; ig<(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
      int mcount = 0; 
      for (uint ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
         if (ic >= old_ncells) break;

         if (mpot[ic] < 0) {
            if (celltype[ic] == REAL_CELL) {
               // remove all but cell that will remain to get count right when split
               // across processors
               if (is_lower_left(i[ic],j[ic]) ) mcount++;
            } else {
               // either upper right or lower left will remain for boundary cells
               if (is_upper_right(i[ic],j[ic]) || is_lower_left(i[ic],j[ic]) ) mcount++;
            }
         }
         if (mpot[ic] >= 0) {
            if (celltype[ic] == REAL_CELL){
               mcount += mpot[ic] ? 4 : 1;
            } else {
               mcount += mpot[ic] ? 2 : 1;
            }
         }
      }    
      if (mtotal != ioffset[ig]) printf("%d: DEBUG ig %d ioffset %d mtotal %d\n",mype,ig,ioffset[ig],mtotal);
      mtotal += mcount;
   }    

   // For global This compares ioffset for each block in the calculation
   ezcl_enqueue_read_buffer(command_queue, dev_ioffset_global, CL_TRUE, 0, block_size_global*sizeof(cl_int), &ioffset_global[0], NULL);
   mtotal = 0; 
   int count = 0; 
   for (uint ig=0; ig<(old_ncells_global+TILE_SIZE-1)/TILE_SIZE; ig++){
      int mcount = 0; 
      for (uint ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
         if (ic >= old_ncells_global) break;

         if (mpot_global[ic] < 0) {
            if (celltype_global[ic] == REAL_CELL) {
               // remove all but cell that will remain to get count right when split
               // across processors
               if (is_lower_left(i_global[ic],j_global[ic]) ) mcount++;
            } else {
               // either upper right or lower left will remain for boundary cells
               if (is_upper_right(i_global[ic],j_global[ic]) || is_lower_left(i_global[ic],j_global[ic]) ) mcount++;
            }
         }

         if (mpot_global[ic] >= 0) {
            if (celltype_global[ic] == REAL_CELL) {
               mcount += mpot_global[ic] ? 4 : 1; 
            } else {
               mcount += mpot_global[ic] ? 2 : 1; 
            }
         }    
      }    
      if (mtotal != ioffset_global[ig]) {
         printf("DEBUG global ig %d ioffset %d mtotal %d\n",ig,ioffset_global[ig],mtotal);
         count++;
      }    
      if (count > 10) exit(0);
      mtotal += mcount;
   }    
}
#endif

Mesh::Mesh(int nx, int ny, int levmx_in, int ndim_in, double deltax_in, double deltay_in, int boundary, int parallel_in, int do_gpu_calc)
{
   lowerBound_Global = NULL;
   upperBound_Global = NULL;
   for (int ii = 0; ii < MESH_TIMER_SIZE; ii++){
      cpu_timers[ii] = 0.0;
      gpu_timers[ii] = 0L;
   }

   for (int ii = 0; ii < MESH_COUNTER_SIZE; ii++){
      cpu_counters[ii] = 0;
      gpu_counters[ii] = 0;
   }

   ndim   = ndim_in;
   levmx  = levmx_in;
#ifdef HAVE_OPENCL
   if (ndim == TWO_DIMENSIONAL) defines = "-DTWO_DIMENSIONAL -DCARTESIAN";
#endif

   offtile_ratio_local = 0;
   offtile_local_count = 1;

   mype  = 0;
   numpe = 1;
   ncells = 0;
   ncells_ghost = 0;
   parallel = parallel_in;
   noffset = 0;
   mem_factor = 1.0;
   //mem_factor = 1.5;
   
#ifdef HAVE_MPI
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init && parallel){
      MPI_Comm_rank(MPI_COMM_WORLD,&mype);
      MPI_Comm_size(MPI_COMM_WORLD,&numpe);
   }
   // TODO add fini
   if (parallel) mesh_memory.pinit(MPI_COMM_WORLD, 2L * 1024 * 1024 * 1024);
#endif
   cell_handle = 0;

   if (numpe == 1) mem_factor = 1.0;

   deltax = deltax_in;
   deltay = deltay_in;

   have_boundary = boundary;

   //int istart = 1;
   //int jstart = 1;
   //int iend   = nx;
   //int jend   = ny;
   int nxx    = nx;
   int nyy    = ny;
   imin = 0;
   jmin = 0;
   imax = nx+1;
   jmax = ny+1;
   if (have_boundary) {
      //istart = 0;
      //jstart = 0;
      //iend   = nx + 1;
      //jend   = ny + 1;
      nxx    = nx + 2;
      nyy    = ny + 2;
      imin   = 0;
      jmin   = 0;
      imax   = nx + 1;
      jmax   = ny + 1;
   }
   
   xmin = -deltax * 0.5 * (real_t)nxx;
   ymin = -deltay * 0.5 * (real_t)nyy;
   xmax =  deltax * 0.5 * (real_t)nxx;
   ymax =  deltay * 0.5 * (real_t)nyy;
   
   size_t lvlMxSize = levmx + 1;

   levtable.resize(lvlMxSize);
   lev_ibegin.resize(lvlMxSize);
   lev_jbegin.resize(lvlMxSize);
   lev_iend.resize(  lvlMxSize);
   lev_jend.resize(  lvlMxSize);
   lev_deltax.resize(lvlMxSize);
   lev_deltay.resize(lvlMxSize);
   lev_iregmin.resize(lvlMxSize);
   lev_iregsize.resize(lvlMxSize);
   lev_jregmin.resize(lvlMxSize);
   lev_jregsize.resize(lvlMxSize);
   
   lev_ibegin[0] = imin + 1;
   lev_iend[0]   = imax - 1;
   lev_jbegin[0] = jmin + 1;
   lev_jend[0]   = jmax - 1;
   lev_deltax[0] = deltax;
   lev_deltay[0] = deltay;
   
   for (int lev = 1; lev <= levmx; lev++) {
      lev_ibegin[lev] = lev_ibegin[lev-1]*2;
      lev_iend[lev]   = lev_iend  [lev-1]*2 + 1;
      lev_jbegin[lev] = lev_jbegin[lev-1]*2;
      lev_jend[lev]   = lev_jend  [lev-1]*2 + 1;
      lev_deltax[lev] = lev_deltax[lev-1]*0.5;
      lev_deltay[lev] = lev_deltay[lev-1]*0.5;
   }
   for (uint lev=0; lev<lvlMxSize; lev++){
      levtable[lev] = IPOW2(lev);
   }

   if (do_gpu_calc) {
#ifdef HAVE_OPENCL
   // The copy host ptr flag will have the data copied to the GPU as part of the allocation
      dev_levtable = ezcl_malloc(&levtable[0],   const_cast<char *>("dev_levtable"), &lvlMxSize, sizeof(cl_int),    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levdx    = ezcl_malloc(&lev_deltax[0], const_cast<char *>("dev_levdx"),    &lvlMxSize, sizeof(cl_real_t), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levdy    = ezcl_malloc(&lev_deltay[0], const_cast<char *>("dev_levdy"),    &lvlMxSize, sizeof(cl_real_t), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levibeg  = ezcl_malloc(&lev_ibegin[0], const_cast<char *>("dev_levibeg"),  &lvlMxSize, sizeof(cl_int),    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_leviend  = ezcl_malloc(&lev_iend[0],   const_cast<char *>("dev_leviend"),  &lvlMxSize, sizeof(cl_int),    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levjbeg  = ezcl_malloc(&lev_jbegin[0], const_cast<char *>("dev_levjbeg"),  &lvlMxSize, sizeof(cl_int),    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levjend  = ezcl_malloc(&lev_jend[0],   const_cast<char *>("dev_levjend"),  &lvlMxSize, sizeof(cl_int),    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
#endif
   }

   ibase = 0;

   int ncells_corners = 4;
   int i_corner[] = {   0,   0,imax,imax};
   int j_corner[] = {   0,jmax,   0,jmax};

   for(int ic=0; ic<ncells_corners; ic++){
      for (int    jj = j_corner[ic]*IPOW2(levmx); jj < (j_corner[ic]+1)*IPOW2(levmx); jj++) {
         for (int ii = i_corner[ic]*IPOW2(levmx); ii < (i_corner[ic]+1)*IPOW2(levmx); ii++) {
            corners_i.push_back(ii);
            corners_j.push_back(jj);
         }
      }
   }

   do_rezone = true;
   gpu_do_rezone = true;
   bidiralloc = true;
   bidirdealloc = true;
   firstFlag = true;

   celltype = NULL;
   nlft     = NULL;
   nrht     = NULL;
   nbot     = NULL;
   ntop     = NULL;
   i        = NULL;
   j        = NULL;
   level    = NULL;

#ifdef PATTERN_CHECK
   for (int ii=0; ii<255; ii++){
       xcase_count[ii]=0;
       sprintf(xcase_descrip[ii],"\0");
   }
#endif

}

void Mesh::init(int nx, int ny, real_t circ_radius, partition_method initial_order, int do_gpu_calc)
{
   if (do_gpu_calc) {
#ifdef HAVE_OPENCL
      cl_context context = ezcl_get_context();

      hash_lib_init();
      if (mype == 0) printf("Starting compile of kernels in mesh\n");
      char *bothsources = (char *)malloc(strlen(mesh_kern_source)+strlen(get_hash_kernel_source_string())+1);
      strcpy(bothsources, get_hash_kernel_source_string());
      strcat(bothsources, mesh_kern_source);
      strcat(bothsources, "\0");
      const char *defines = NULL;
      cl_program program = ezcl_create_program_wsource(context, defines, bothsources);
      free(bothsources);

      dev_nlft = NULL;
      kernel_reduction_scan2          = ezcl_create_kernel_wprogram(program, "finish_reduction_scan2_cl");
      kernel_reduction_count          = ezcl_create_kernel_wprogram(program, "finish_reduction_count_cl");
      kernel_reduction_count2         = ezcl_create_kernel_wprogram(program, "finish_reduction_count2_cl");
      kernel_hash_adjust_sizes        = ezcl_create_kernel_wprogram(program, "hash_adjust_sizes_cl");
      kernel_hash_setup               = ezcl_create_kernel_wprogram(program, "hash_setup_cl");
      kernel_hash_setup_local         = ezcl_create_kernel_wprogram(program, "hash_setup_local_cl");
      kernel_neighbor_init            = ezcl_create_kernel_wprogram(program, "neighbor_init_cl");
      kernel_calc_neighbors           = ezcl_create_kernel_wprogram(program, "calc_neighbors_cl");
      kernel_calc_neighbors_local     = ezcl_create_kernel_wprogram(program, "calc_neighbors_local_cl");
      kernel_calc_border_cells        = ezcl_create_kernel_wprogram(program, "calc_border_cells_cl");
      kernel_calc_border_cells2       = ezcl_create_kernel_wprogram(program, "calc_border_cells2_cl");
      kernel_finish_scan              = ezcl_create_kernel_wprogram(program, "finish_scan_cl");
      kernel_get_border_data          = ezcl_create_kernel_wprogram(program, "get_border_data_cl");
      kernel_calc_layer1              = ezcl_create_kernel_wprogram(program, "calc_layer1_cl");
      kernel_calc_layer1_sethash      = ezcl_create_kernel_wprogram(program, "calc_layer1_sethash_cl");
      kernel_calc_layer2              = ezcl_create_kernel_wprogram(program, "calc_layer2_cl");
      kernel_get_border_data2         = ezcl_create_kernel_wprogram(program, "get_border_data2_cl");
      kernel_calc_layer2_sethash      = ezcl_create_kernel_wprogram(program, "calc_layer2_sethash_cl");
      kernel_copy_mesh_data           = ezcl_create_kernel_wprogram(program, "copy_mesh_data_cl");
      kernel_fill_mesh_ghost          = ezcl_create_kernel_wprogram(program, "fill_mesh_ghost_cl");
      kernel_fill_neighbor_ghost      = ezcl_create_kernel_wprogram(program, "fill_neighbor_ghost_cl");
      kernel_set_corner_neighbor      = ezcl_create_kernel_wprogram(program, "set_corner_neighbor_cl");
      kernel_adjust_neighbors_local   = ezcl_create_kernel_wprogram(program, "adjust_neighbors_local_cl");
      kernel_hash_size                = ezcl_create_kernel_wprogram(program, "calc_hash_size_cl");
      kernel_finish_hash_size         = ezcl_create_kernel_wprogram(program, "finish_reduction_minmax4_cl");
      kernel_calc_spatial_coordinates = ezcl_create_kernel_wprogram(program, "calc_spatial_coordinates_cl");
      kernel_do_load_balance_lower    = ezcl_create_kernel_wprogram(program, "do_load_balance_lower_cl");
      kernel_do_load_balance_middle   = ezcl_create_kernel_wprogram(program, "do_load_balance_middle_cl");
      kernel_do_load_balance_upper    = ezcl_create_kernel_wprogram(program, "do_load_balance_upper_cl");
      kernel_face_idx_wbidirmap = ezcl_create_kernel_wprogram(program, "face_idx_wbidirmap_cl");
      kernel_calc_face_list_wbidirmap_pt1 = ezcl_create_kernel_wprogram(program, "calc_face_list_wbidirmap_pt1_cl");
      kernel_calc_face_list_wbidirmap_pt2 = ezcl_create_kernel_wprogram(program, "calc_face_list_wbidirmap_pt2_cl");
      kernel_wbidirmap_precount = ezcl_create_kernel_wprogram(program, "wbidirmap_precount_cl");
      kernel_calc_wbidirmap_phantom_neighbors = ezcl_create_kernel_wprogram(program, "calc_wbidirmap_phantom_neighbors_cl");
      kernel_calc_wbidirmap_phantom_values = ezcl_create_kernel_wprogram(program, "calc_wbidirmap_phantom_values_cl");
      kernel_deep_copy = ezcl_create_kernel_wprogram(program, "deep_copy_cl");
#ifndef MINIMUM_PRECISION
      kernel_do_load_balance_double   = ezcl_create_kernel_wprogram(program, "do_load_balance_double_cl");
#endif
      kernel_do_load_balance_float    = ezcl_create_kernel_wprogram(program, "do_load_balance_float_cl");
      kernel_refine_smooth            = ezcl_create_kernel_wprogram(program, "refine_smooth_cl");
      kernel_coarsen_smooth           = ezcl_create_kernel_wprogram(program, "coarsen_smooth_cl");
      kernel_coarsen_check_block      = ezcl_create_kernel_wprogram(program, "coarsen_check_block_cl");
      kernel_rezone_all               = ezcl_create_kernel_wprogram(program, "rezone_all_cl");
      kernel_rezone_neighbors         = ezcl_create_kernel_wprogram(program, "rezone_neighbors_cl");
#ifndef MINIMUM_PRECISION
      kernel_rezone_one_double        = ezcl_create_kernel_wprogram(program, "rezone_one_double_cl");
#endif
      kernel_rezone_one_float         = ezcl_create_kernel_wprogram(program, "rezone_one_float_cl");
      kernel_copy_mpot_ghost_data     = ezcl_create_kernel_wprogram(program, "copy_mpot_ghost_data_cl");
      kernel_set_boundary_refinement  = ezcl_create_kernel_wprogram(program, "set_boundary_refinement");
      init_kernel_2stage_sum();
      init_kernel_2stage_sum_int();
      if (! have_boundary){
        kernel_count_BCs              = ezcl_create_kernel_wprogram(program, "count_BCs_cl");
      }

      ezcl_program_release(program);
      if (mype == 0) printf("Finishing compile of kernels in mesh\n");
#endif
   }

   //KDTree_Initialize(&tree);
   if (ncells > 0) { // this is a restart.
       if (parallel && numpe > 1) {
#ifdef HAVE_MPI
          int ncells_int = ncells;
          nsizes.clear();
          ndispl.clear();
          nsizes.resize (numpe);
          ndispl.resize (numpe);
          MPI_Allgather(&ncells_int, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
          ndispl[0]=0;
          for (int ip=1; ip<numpe; ip++){
             ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
          }
          noffset=ndispl[mype];
          ncells_global = ndispl[numpe-1] + nsizes[numpe-1];
#endif
       } else {
          noffset = 0;
          ncells_global = ncells;
          proc.resize (ncells);
          calc_distribution(numpe);
       }
       calc_celltype(ncells);

   } else {
       int istart = 1,
           jstart = 1,
           iend   = nx,
           jend   = ny,
           nxx    = nx,
           nyy    = ny;
       if (have_boundary) {
          istart = 0;
          jstart = 0;
          iend   = nx + 1;
          jend   = ny + 1;
          nxx    = nx + 2;
          nyy    = ny + 2;
       }

       if (ndim == TWO_DIMENSIONAL) ncells = nxx * nyy - have_boundary * 4;
       else                         ncells = nxx * nyy;

       noffset = 0;
       if (parallel) {
          ncells_global = ncells;
          nsizes.resize(numpe);
          ndispl.resize(numpe);

          for (int ip=0; ip<numpe; ip++){
             nsizes[ip] = ncells_global/numpe;
             if (ip < (int)(ncells_global%numpe)) nsizes[ip]++;
          }

          ndispl[0]=0;
          for (int ip=1; ip<numpe; ip++){
             ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
          }
          ncells= nsizes[mype];
          noffset=ndispl[mype];
       }

       allocate(ncells);
       index.resize(ncells);

       int ic = 0;

       for (int jj = jstart; jj <= jend; jj++) {
          for (int ii = istart; ii <= iend; ii++) {
             if (have_boundary && ii == 0    && jj == 0   ) continue;
             if (have_boundary && ii == 0    && jj == jend) continue;
             if (have_boundary && ii == iend && jj == 0   ) continue;
             if (have_boundary && ii == iend && jj == jend) continue;

             if (ic >= (int)noffset && ic < (int)(ncells+noffset)){
                int iclocal = ic-noffset;
                index[iclocal] = ic;
                i[iclocal]     = ii;
                j[iclocal]     = jj;
                level[iclocal] = 0;
             }
             ic++;
          }
       }

       //if (numpe > 1 && (initial_order != HILBERT_SORT && initial_order != HILBERT_PARTITION) ) mem_factor = 2.0;
       partition_cells(numpe, index, initial_order);

       calc_celltype(ncells);
       calc_spatial_coordinates(0);

       //  Start lev loop here
       for (int ilevel=1; ilevel<=levmx; ilevel++) {

          //int old_ncells = ncells;

          ncells_ghost = ncells;
          calc_neighbors_local();

          kdtree_setup();

          int nez;
          vector<int> ind(ncells);

    #ifdef FULL_PRECISION
          KDTree_QueryCircleIntersect_Double(&tree, &nez, &(ind[0]), circ_radius, ncells, &x[0], &dx[0], &y[0], &dy[0]);
    #else
          KDTree_QueryCircleIntersect_Float(&tree, &nez, &(ind[0]), circ_radius, ncells, &x[0], &dx[0], &y[0], &dy[0]);
    #endif

          vector<char_t> mpot(ncells_ghost,0);

          for (int ic=0; ic<nez; ++ic){
             if ((int)level[ind[ic]] < levmx) mpot[ind[ic]] = 1;
          }

          KDTree_Destroy(&tree);
          //  Refine the cells.
          int icount = 0;
          int jcount = 0;
          int new_ncells = refine_smooth(mpot, icount, jcount);

          MallocPlus dummy;
          rezone_all(icount, jcount, mpot, 0, dummy);

          ncells = new_ncells;

          calc_spatial_coordinates(0);

    #ifdef HAVE_MPI
          if (parallel && numpe > 1) {
             int ncells_int = ncells;
             MPI_Allgather(&ncells_int, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
             ndispl[0]=0;
             for (int ip=1; ip<numpe; ip++){
                ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
             }
             noffset=ndispl[mype];
             ncells_global = ndispl[numpe-1] + nsizes[numpe-1];
          }
    #endif
       }  // End lev loop here
       index.clear();
       ncells_ghost = ncells;
   }
   int ncells_corners = 4;
   int i_corner[] = {   0,   0,imax,imax};
   int j_corner[] = {   0,jmax,   0,jmax};

   for(int ic=0; ic<ncells_corners; ic++){
      for (int    jj = j_corner[ic]*IPOW2(levmx); jj < (j_corner[ic]+1)*IPOW2(levmx); jj++) {
         for (int ii = i_corner[ic]*IPOW2(levmx); ii < (i_corner[ic]+1)*IPOW2(levmx); ii++) {
            corners_i.push_back(ii);
            corners_j.push_back(jj);
         }
      }
   }
}

size_t Mesh::refine_smooth(vector<char_t> &mpot, int &icount, int &jcount)
{
   vector<char_t> mpot_old;

   int newcount;
   int newcount_global;

   struct timespec tstart_lev2;

   struct mesh_type
   {
#ifdef FULL_PRECISION
       double ***pstate;
#else
       float ***pstate;
#endif
       int **mask;
   };

#ifdef _OPENMP
#pragma omp parallel
{ //START Parallel Region
#endif

   rezone_count_threaded(mpot, icount, jcount);

#ifdef _OPENMP
#pragma omp master
{//MASTER START
#endif
   newcount = icount;
   newcount_global = newcount;

   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

#ifdef HAVE_MPI
   if (parallel) {
      MPI_Allreduce(&newcount, &newcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   }
#endif

#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

   if(newcount_global > 0 && levmx > 1) {

      size_t my_ncells=ncells;
      if (parallel) my_ncells=ncells_ghost;

#ifdef _OPENMP
#pragma omp master
{//MASTER START
#endif
      cpu_counters[MESH_COUNTER_REFINE_SMOOTH]++;

      mpot_old.resize(my_ncells);
#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

      int levcount = 1;
       
      while (newcount_global > 0 && levcount < levmx){

         levcount++; 
#ifdef _OPENMP
#pragma omp master
{//MASTER START
#endif

         mpot.swap(mpot_old);
         newcount=0;
#ifdef HAVE_MPI
         if (numpe > 1) {
            L7_Update(&mpot_old[0], L7_CHAR_T, cell_handle);
         }
#endif

#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

         int upperBound, lowerBound;
         get_bounds(upperBound, lowerBound);
         int mynewcount = newcount; //All threads get a mynewcount

#ifdef _OPENMP
#pragma omp for reduction(+:newcount)
#endif
         for(uint ic = 0; ic < ncells; ic++) {
        // for(uint ic = lowerBound; ic < upperBound; ic++){
            int lev = level[ic];
            mpot[ic] = mpot_old[ic];
            if(mpot_old[ic] > 0) continue;
   
            int nl = nlft[ic];
            if (nl >= 0 && nl < (int)ncells_ghost) {
               int ll = level[nl];
               if(mpot_old[nl] > 0) ll++;
   
               if(ll - lev > 1) {
                  mpot[ic]=1;
                  mynewcount++;
                  continue;
               }

               ll = level[nl];
               if (ll > lev) {
                  int nlt = ntop[nl];
                  if (nlt >= 0 && nlt < (int)ncells_ghost) {
                     int llt = level[nlt];
                     if(mpot_old[nlt] > 0) llt++;

                     if(llt - lev > 1) {
                        mpot[ic]=1;
                        mynewcount++;
                        continue;
                     }
                  }
               }
            }

            int nr = nrht[ic];
            if (nr >= 0 && nr < (int)ncells_ghost) {
               int lr = level[nr];
               if(mpot_old[nr] > 0) lr++;
   
               if(lr - lev > 1) {
                  mpot[ic]=1;
                  mynewcount++;
                  continue;
               }

               lr = level[nr];
               if (lr > lev) {
                  int nrt = ntop[nr];
                  if (nrt >= 0 && nrt < (int)ncells_ghost) {
                     int lrt = level[nrt];
                     if(mpot_old[nrt] > 0) lrt++;

                     if(lrt - lev > 1) {
                        mpot[ic]=1;
                        mynewcount++;
                        continue;
                     }
                  }
               }
            }

            int nt = ntop[ic];
            if (nt >= 0 && nt < (int)ncells_ghost) {
               int lt = level[nt];
               if(mpot_old[nt] > 0) lt++;
   
               if(lt - lev > 1) {
                  mpot[ic]=1;
                  mynewcount++;
                  continue;
               }

               lt = level[nt];
               if (lt > lev) {
                  int ntr = nrht[nt];
                  if (ntr >= 0 && ntr < (int)ncells_ghost) {
                     int ltr = level[ntr];
                     if(mpot_old[ntr] > 0) ltr++;

                     if(ltr - lev > 1) {
                        mpot[ic]=1;
                        mynewcount++;
                        continue;
                     }
                  }
               }
            }

            int nb = nbot[ic];
            if (nb >= 0 && nb < (int)ncells_ghost) {
               int lb = level[nb];
               if(mpot_old[nb] > 0) lb++;
   
               if(lb - lev > 1) {
                  mpot[ic]=1;
                  mynewcount++;
                  continue;
               }

               lb = level[nb];
               if (lb > lev) {
                  int nbr = nrht[nb];
                  if (nbr >= 0 && nbr < (int)ncells_ghost) {
                     int lbr = level[nbr];
                     if(mpot_old[nbr] > 0) lbr++;

                     if(lbr - lev > 1) {
                        mpot[ic]=1;
                        mynewcount++;
                        continue;
                     }
                  }
               }
            }
         }
#ifdef _OPENMP
#pragma omp atomic 
#endif
         newcount += mynewcount;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
{
#endif
         icount += newcount;
         newcount_global = newcount;

#ifdef HAVE_MPI
         if (parallel) {
            MPI_Allreduce(&newcount, &newcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         }
#endif

#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

      } // while (newcount_global > 0 && levcount < levmx);

   }


#ifdef _OPENMP
#pragma omp master
{
#endif

#ifdef HAVE_MPI
   if (numpe > 1) {
      L7_Update(&mpot[0], L7_CHAR_T, cell_handle);
  }
#endif

   mpot_old.clear();
   mpot_old.resize(ncells_ghost);

   mpot_old.swap(mpot);
#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
   for(uint ic=0; ic<ncells; ic++) {
      mpot[ic] = mpot_old[ic];
      if (mpot_old[ic] >= 0) continue;
      if (mpot_old[ic] <= -1000000) continue;
      if (        is_upper_right(i[ic],j[ic]) ) {
         int nr = nrht[ic];
         uint lr = level[nr];
         if (mpot_old[nr] > 0) lr++;
         int nt = ntop[ic];
         uint lt = level[nt];
         if (mpot_old[nt] > 0) lt++;
         if (lr > level[ic] || lt > level[ic]) mpot[ic] = 0;
      } else if ( is_upper_left(i[ic],j[ic] ) ) {
         int nl = nlft[ic];
         uint ll = level[nl];
         if (mpot_old[nl] > 0) ll++;
         int nt = ntop[ic];
         uint lt = level[nt];
         if (mpot_old[nt] > 0) lt++;
         if (ll > level[ic] || lt > level[ic]) mpot[ic] = 0;
      } else if ( is_lower_right(i[ic],j[ic] ) ) {
         int nr = nrht[ic];
         uint lr = level[nr];
         if (mpot_old[nr] > 0) lr++;
         int nb = nbot[ic];
         uint lb = level[nb];
         if (mpot_old[nb] > 0) lb++;
         if (lr > level[ic] || lb > level[ic]) mpot[ic] = 0;
      } else if ( is_lower_left(i[ic],j[ic] ) ) {
         int nl = nlft[ic];
         uint ll = level[nl];
         if (mpot_old[nl] > 0) ll++;
         int nb = nbot[ic];
         uint lb = level[nb];
         if (mpot_old[nb] > 0) lb++;
         if (ll > level[ic] || lb > level[ic]) mpot[ic] = 0;
      }
   }

#ifdef _OPENMP
#pragma omp master
{
#endif

#ifdef HAVE_MPI
   if (numpe > 1) {
      L7_Update(&mpot[0], L7_CHAR_T, cell_handle);
  }
#endif

   mpot_old.swap(mpot);
#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
   for(uint ic=0; ic<ncells; ic++) {
      int n1=0, n2=0, n3=0;
      mpot[ic] = mpot_old[ic];
      if (mpot_old[ic] >= 0) continue;
      if (mpot_old[ic] <= -1000000) continue;
      if ( is_upper_right(i[ic],j[ic]) ) {
         n1 = nbot[ic];
         n2 = nlft[ic];
         n3 = nlft[n1];
      } else if ( is_upper_left(i[ic],j[ic] ) ) {
         n1 = nbot[ic];
         n2 = nrht[ic];
         n3 = nrht[n1];
      } else if ( is_lower_right(i[ic],j[ic] ) ) {
         n1 = ntop[ic];
         n2 = nlft[ic];
         n3 = nlft[n1];
      } else if ( is_lower_left(i[ic],j[ic] ) ) {
         n1 = ntop[ic];
         n2 = nrht[ic];
         n3 = nrht[n1];
      }
      if (n3 < 0) {
         mpot[ic] = 0;
      } else {
         uint lev1 = level[n1];
         uint lev2 = level[n2];
         uint lev3 = level[n3];
         if (mpot_old[n1] > 0) lev1++;
         if (mpot_old[n2] > 0) lev2++;
         if (mpot_old[n3] > 0) lev3++;

         if (mpot_old[n1] != -1 || lev1 != level[ic] ||
             mpot_old[n2] != -1 || lev2 != level[ic] ||
             mpot_old[n3] != -1 || lev3 != level[ic]) {
            mpot[ic] = 0;
         }
      }
   }

#ifdef _OPENMP
#pragma omp master
{
#endif

#ifdef HAVE_MPI
   if (numpe > 1) {
      L7_Update(&mpot[0], L7_CHAR_T, cell_handle);
  }
#endif

#ifdef _OPENMP
}//END MASTER
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
   for (uint ic=0; ic<ncells; ic++) {
      if (celltype[ic] < 0) {
         switch (celltype[ic]) {
            case LEFT_BOUNDARY:
               mpot[ic] = mpot[nrht[ic]];
               break;
            case RIGHT_BOUNDARY:
               mpot[ic] = mpot[nlft[ic]];
               break;
            case BOTTOM_BOUNDARY:
               mpot[ic] = mpot[ntop[ic]];
               break;
            case TOP_BOUNDARY:
               mpot[ic] = mpot[nbot[ic]];
               break;
         }
      }
   }

#ifdef _OPENMP
#pragma omp barrier
}//END Parallel Region
#endif

   newcount = ncells + rezone_count(mpot, icount, jcount);

#ifdef HAVE_MPI
   int icount_global = icount;
   int jcount_global = jcount;
   if (parallel) {
      int count[2], global_count[2];
      count[0] = icount;
      count[1] = jcount;
      MPI_Allreduce(&count, &global_count, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      icount_global = global_count[0];
      jcount_global = global_count[1];
   }
   do_rezone = (icount_global != 0 || jcount_global != 0) ? true : false;
#else
   do_rezone = (icount != 0 || jcount != 0) ? true : false;
#endif


   if (TIMING_LEVEL >= 2) cpu_timers[MESH_TIMER_REFINE_SMOOTH] += cpu_timer_stop(tstart_lev2);

   return(newcount);
}

#ifdef HAVE_OPENCL
int Mesh::gpu_refine_smooth(cl_mem &dev_mpot, int &icount, int &jcount)
{
   struct timespec tstart_lev2;
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size = global_work_size/local_work_size;

   int icount_global = icount;
   int jcount_global = jcount;

#ifdef HAVE_MPI
   if (parallel) {
      int count[2], count_global[2];
      count[0] = icount;
      count[1] = jcount;
      MPI_Allreduce(&count, &count_global, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      icount_global = count_global[0];
      jcount_global = count_global[1];
   }
#endif

   int levcount = 1;
   //int which_smooth=0;

   if(icount_global > 0 && levcount < levmx) {
      size_t result_size = 1;
      cl_mem dev_result  = ezcl_malloc(NULL, const_cast<char *>("dev_result"),  &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_mpot_old = ezcl_malloc(NULL, const_cast<char *>("dev_mpot_old"), &ncells_ghost, sizeof(cl_char_t), CL_MEM_READ_WRITE, 0);

      int newcount = icount;
      int newcount_global = icount_global;
      while (newcount_global > 0 && levcount < levmx) {
         levcount++;

         gpu_counters[MESH_COUNTER_REFINE_SMOOTH]++;

#ifdef HAVE_MPI
         if (numpe > 1) {
            L7_Dev_Update(dev_mpot, L7_CHAR_T, cell_handle);
         }
#endif

         if (icount_global) {
            ezcl_device_memory_swap(&dev_mpot_old, &dev_mpot);

            ezcl_set_kernel_arg(kernel_refine_smooth, 0, sizeof(cl_int),  (void *)&ncells);
            ezcl_set_kernel_arg(kernel_refine_smooth, 1, sizeof(cl_int),  (void *)&ncells_ghost);
            ezcl_set_kernel_arg(kernel_refine_smooth, 2, sizeof(cl_int),  (void *)&levmx);
            ezcl_set_kernel_arg(kernel_refine_smooth, 3, sizeof(cl_mem),  (void *)&dev_nlft);
            ezcl_set_kernel_arg(kernel_refine_smooth, 4, sizeof(cl_mem),  (void *)&dev_nrht);
            ezcl_set_kernel_arg(kernel_refine_smooth, 5, sizeof(cl_mem),  (void *)&dev_nbot);
            ezcl_set_kernel_arg(kernel_refine_smooth, 6, sizeof(cl_mem),  (void *)&dev_ntop);
            ezcl_set_kernel_arg(kernel_refine_smooth, 7, sizeof(cl_mem),  (void *)&dev_level);
            ezcl_set_kernel_arg(kernel_refine_smooth, 8, sizeof(cl_mem),  (void *)&dev_celltype);
            ezcl_set_kernel_arg(kernel_refine_smooth, 9, sizeof(cl_mem),  (void *)&dev_mpot_old);
            ezcl_set_kernel_arg(kernel_refine_smooth,10, sizeof(cl_mem),  (void *)&dev_mpot);
            ezcl_set_kernel_arg(kernel_refine_smooth,11, sizeof(cl_mem),  (void *)&dev_redscratch);
            ezcl_set_kernel_arg(kernel_refine_smooth,12, sizeof(cl_mem),  (void *)&dev_result);
            ezcl_set_kernel_arg(kernel_refine_smooth,13, local_work_size*sizeof(cl_int),    NULL);

            ezcl_enqueue_ndrange_kernel(command_queue, kernel_refine_smooth, 1, NULL, &global_work_size, &local_work_size, NULL);

            gpu_rezone_count(block_size, local_work_size, dev_redscratch, dev_result);

            int result;
            ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int), &result, NULL);

            //printf("result = %d after %d refine smooths\n",result,which_smooth);
            //which_smooth++;

            icount = result;
         }

         newcount = icount-newcount;
         newcount_global = newcount;
#ifdef HAVE_MPI
         if (parallel) {
            MPI_Allreduce(&newcount, &newcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         }
#endif
         icount_global += newcount_global;
         //printf("DEBUG -- icount %d icount_global %d newcount %d newcount_global %d\n",icount,icount_global,newcount,newcount_global);
      }

      ezcl_device_memory_delete(dev_mpot_old);
      ezcl_device_memory_delete(dev_redscratch);
      ezcl_device_memory_delete(dev_result);
   }

   if (jcount_global) {
#ifdef HAVE_MPI
      if (numpe > 1) {
         L7_Dev_Update(dev_mpot, L7_CHAR_T, cell_handle);
      }
#endif

      cl_mem dev_mpot_old = ezcl_malloc(NULL, const_cast<char *>("dev_mpot_old"), &ncells_ghost, sizeof(cl_char_t), CL_MEM_READ_WRITE, 0);

      if (jcount) {
         ezcl_device_memory_swap(&dev_mpot_old, &dev_mpot);

         ezcl_set_kernel_arg(kernel_coarsen_smooth, 0, sizeof(cl_int),  (void *)&ncells);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 1, sizeof(cl_mem),  (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 2, sizeof(cl_mem),  (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 3, sizeof(cl_mem),  (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 4, sizeof(cl_mem),  (void *)&dev_ntop);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 5, sizeof(cl_mem),  (void *)&dev_i);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 6, sizeof(cl_mem),  (void *)&dev_j);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 7, sizeof(cl_mem),  (void *)&dev_level);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 8, sizeof(cl_mem),  (void *)&dev_mpot_old);
         ezcl_set_kernel_arg(kernel_coarsen_smooth, 9, sizeof(cl_mem),  (void *)&dev_mpot);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_coarsen_smooth, 1, NULL, &global_work_size, &local_work_size, NULL);
      }

#ifdef HAVE_MPI
      if (numpe > 1) {
         L7_Dev_Update(dev_mpot, L7_CHAR_T, cell_handle);
      }
#endif

      if (jcount) {
         size_t result_size = 1;
         cl_mem dev_result  = ezcl_malloc(NULL, const_cast<char *>("dev_result"),  &result_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         ezcl_device_memory_swap(&dev_mpot_old, &dev_mpot);

         ezcl_set_kernel_arg(kernel_coarsen_check_block, 0, sizeof(cl_int),  (void *)&ncells);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 1, sizeof(cl_mem),  (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 2, sizeof(cl_mem),  (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 3, sizeof(cl_mem),  (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 4, sizeof(cl_mem),  (void *)&dev_ntop);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 5, sizeof(cl_mem),  (void *)&dev_i);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 6, sizeof(cl_mem),  (void *)&dev_j);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 7, sizeof(cl_mem),  (void *)&dev_level);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 8, sizeof(cl_mem),  (void *)&dev_celltype);
         ezcl_set_kernel_arg(kernel_coarsen_check_block, 9, sizeof(cl_mem),  (void *)&dev_mpot_old);
         ezcl_set_kernel_arg(kernel_coarsen_check_block,10, sizeof(cl_mem),  (void *)&dev_mpot);
         ezcl_set_kernel_arg(kernel_coarsen_check_block,11, sizeof(cl_mem),  (void *)&dev_redscratch);
         ezcl_set_kernel_arg(kernel_coarsen_check_block,12, sizeof(cl_mem),  (void *)&dev_result);
         ezcl_set_kernel_arg(kernel_coarsen_check_block,13, local_work_size*sizeof(cl_int),    NULL);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_coarsen_check_block, 1, NULL, &global_work_size, &local_work_size, NULL);

         gpu_rezone_count(block_size, local_work_size, dev_redscratch, dev_result);

         int result;
         ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int), &result, NULL);

         //printf("result = %d after coarsen smooth\n",result);

         jcount = result;

         ezcl_device_memory_delete(dev_redscratch);
         ezcl_device_memory_delete(dev_result);
      }

      jcount_global = jcount;

#ifdef HAVE_MPI
      if (parallel) {
         MPI_Allreduce(&jcount, &jcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      }
#endif

      ezcl_device_memory_delete(dev_mpot_old);
   }

   if (icount_global || jcount_global) {
#ifdef HAVE_MPI
      if (numpe > 1) {
         L7_Dev_Update(dev_mpot, L7_CHAR_T, cell_handle);
      }
#endif

      size_t result_size = 1;
      cl_mem dev_result  = ezcl_malloc(NULL, const_cast<char *>("dev_result"),  &result_size, sizeof(cl_int2), CL_MEM_READ_WRITE, 0);
      cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_int2), CL_MEM_READ_WRITE, 0);
      dev_ioffset  = ezcl_malloc(NULL, const_cast<char *>("dev_ioffset"), &block_size,   sizeof(cl_uint), CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 1,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 2,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 3,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 4,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 5,  sizeof(cl_mem), (void *)&dev_i);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 6,  sizeof(cl_mem), (void *)&dev_j);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 7,  sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 8,  sizeof(cl_mem), (void *)&dev_mpot);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 9,  sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 10, sizeof(cl_mem), (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 11, sizeof(cl_mem), (void *)&dev_result);
      ezcl_set_kernel_arg(kernel_set_boundary_refinement, 12, local_work_size*sizeof(cl_int2),    NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_set_boundary_refinement, 1, NULL, &global_work_size, &local_work_size, NULL);

      gpu_rezone_count2(block_size, local_work_size, dev_redscratch, dev_result);

      int my_result[2];
      ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, 1*sizeof(cl_int2), &my_result, NULL);
      //printf("Result is %lu icount %d jcount %d\n", ncells+my_result[0]-my_result[1],my_result[0],my_result[1]);
      icount = my_result[0];
      jcount = my_result[1];

      icount_global = icount;
      jcount_global = jcount;
#ifdef HAVE_MPI
      if (parallel) {
         int count[2], count_global[2];
         count[0] = icount;
         count[1] = jcount;
         MPI_Allreduce(&count, &count_global, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         icount_global = count_global[0];
         jcount_global = count_global[1];
      }
#endif

      gpu_rezone_scan(block_size, local_work_size, dev_ioffset, dev_result);

      //ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int), &my_result, NULL);
      //printf("After scan, Result is %d\n", my_result[0]);

      ezcl_device_memory_delete(dev_result);
      ezcl_device_memory_delete(dev_redscratch);

   } else {
      ezcl_device_memory_delete(dev_mpot);
      dev_mpot = NULL;
   }

   gpu_do_rezone = (icount_global != 0 || jcount_global != 0) ? true : false;

   if (TIMING_LEVEL >= 2) gpu_timers[MESH_TIMER_REFINE_SMOOTH] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);

   return ncells+icount-jcount;
}
#endif

void Mesh::terminate(void)
{
      mesh_memory.memory_delete(i);
      mesh_memory.memory_delete(j);
      mesh_memory.memory_delete(level);
      mesh_memory.memory_delete(xface_i);
      mesh_memory.memory_delete(xface_j);
      mesh_memory.memory_delete(xface_level);
      mesh_memory.memory_delete(map_xface2cell_lower);
      mesh_memory.memory_delete(map_xface2cell_upper);
      mesh_memory.memory_delete(map_xcell2face_left1);
      mesh_memory.memory_delete(map_xcell2face_left2);
      mesh_memory.memory_delete(map_xcell2face_right1);
      mesh_memory.memory_delete(map_xcell2face_right2);
      mesh_memory.memory_delete(xrecvIdx);
      mesh_memory.memory_delete(xrecvCIdx);
      mesh_memory.memory_delete(xplusCell2Idx);
      mesh_memory.memory_delete(xminusCell2Idx);
      mesh_memory.memory_delete(xsendIdx1);
      mesh_memory.memory_delete(xsendIdx2);
      mesh_memory.memory_delete(yface_i);
      mesh_memory.memory_delete(yface_j);
      mesh_memory.memory_delete(yface_level);
      mesh_memory.memory_delete(map_yface2cell_lower);
      mesh_memory.memory_delete(map_yface2cell_upper);
      mesh_memory.memory_delete(map_ycell2face_bot1);
      mesh_memory.memory_delete(map_ycell2face_bot2);
      mesh_memory.memory_delete(map_ycell2face_top1);
      mesh_memory.memory_delete(map_ycell2face_top2);
      mesh_memory.memory_delete(yrecvIdx);
      mesh_memory.memory_delete(yrecvCIdx);
      mesh_memory.memory_delete(yplusCell2Idx);
      mesh_memory.memory_delete(yminusCell2Idx);
      mesh_memory.memory_delete(ysendIdx1);
      mesh_memory.memory_delete(ysendIdx2);
      // needs to cast char_t to void so doesn't mistake it for string
      mesh_memory.memory_delete((void *)celltype);
      if (neighbor_remap) {
         mesh_memory.memory_delete(nlft);
         mesh_memory.memory_delete(nrht);
         mesh_memory.memory_delete(nbot);
         mesh_memory.memory_delete(ntop);
      }

#ifdef _OPENMP
#pragma omp master
#endif
      {
         free(lowerBound_Global);
         free(upperBound_Global);
      }

#ifdef HAVE_OPENCL
      hash_lib_terminate();

      ezcl_device_memory_delete(dev_levtable);
      ezcl_device_memory_delete(dev_levdx);
      ezcl_device_memory_delete(dev_levdy);
      ezcl_device_memory_delete(dev_levibeg);
      ezcl_device_memory_delete(dev_leviend);
      ezcl_device_memory_delete(dev_levjbeg);
      ezcl_device_memory_delete(dev_levjend);

      ezcl_device_memory_delete(dev_level);
      ezcl_device_memory_delete(dev_i);
      ezcl_device_memory_delete(dev_j);
      ezcl_device_memory_delete(dev_celltype);
      if (neighbor_remap && dev_nlft != NULL){
         ezcl_device_memory_delete(dev_nlft);
         ezcl_device_memory_delete(dev_nrht);
         ezcl_device_memory_delete(dev_nbot);
         ezcl_device_memory_delete(dev_ntop);
      }

      ezcl_kernel_release(kernel_reduction_scan2);
      ezcl_kernel_release(kernel_reduction_count);
      ezcl_kernel_release(kernel_reduction_count2);
      ezcl_kernel_release(kernel_hash_adjust_sizes);
      ezcl_kernel_release(kernel_hash_setup);
      ezcl_kernel_release(kernel_hash_setup_local);
      ezcl_kernel_release(kernel_neighbor_init);
      ezcl_kernel_release(kernel_calc_neighbors);
      ezcl_kernel_release(kernel_calc_neighbors_local);
      ezcl_kernel_release(kernel_calc_border_cells);
      ezcl_kernel_release(kernel_calc_border_cells2);
      ezcl_kernel_release(kernel_finish_scan);
      ezcl_kernel_release(kernel_get_border_data);
      ezcl_kernel_release(kernel_calc_layer1);
      ezcl_kernel_release(kernel_calc_layer1_sethash);
      ezcl_kernel_release(kernel_calc_layer2);
      ezcl_kernel_release(kernel_get_border_data2);
      ezcl_kernel_release(kernel_calc_layer2_sethash);
      //ezcl_kernel_release(kernel_calc_neighbors_local2);
      ezcl_kernel_release(kernel_copy_mesh_data);
      ezcl_kernel_release(kernel_fill_mesh_ghost);
      ezcl_kernel_release(kernel_fill_neighbor_ghost);
      ezcl_kernel_release(kernel_set_corner_neighbor);
      ezcl_kernel_release(kernel_adjust_neighbors_local);
      //ezcl_kernel_release(kernel_copy_ghost_data);
      //ezcl_kernel_release(kernel_adjust_neighbors);
      ezcl_kernel_release(kernel_hash_size);
      ezcl_kernel_release(kernel_finish_hash_size);
      ezcl_kernel_release(kernel_calc_spatial_coordinates);
      ezcl_kernel_release(kernel_do_load_balance_lower);
      ezcl_kernel_release(kernel_do_load_balance_middle);
      ezcl_kernel_release(kernel_do_load_balance_upper);
      ezcl_kernel_release(kernel_calc_face_list_wbidirmap_pt1);
      ezcl_kernel_release(kernel_calc_face_list_wbidirmap_pt2);
      ezcl_kernel_release(kernel_face_idx_wbidirmap);
      ezcl_kernel_release(kernel_wbidirmap_precount);
      ezcl_kernel_release(kernel_calc_wbidirmap_phantom_neighbors);
      ezcl_kernel_release(kernel_calc_wbidirmap_phantom_values);
      ezcl_kernel_release(kernel_deep_copy);
#ifndef MINIMUM_PRECISION
      ezcl_kernel_release(kernel_do_load_balance_double);
#endif
      ezcl_kernel_release(kernel_do_load_balance_float);
      ezcl_kernel_release(kernel_refine_smooth);
      ezcl_kernel_release(kernel_coarsen_smooth);
      ezcl_kernel_release(kernel_coarsen_check_block);
      ezcl_kernel_release(kernel_rezone_all);
      ezcl_kernel_release(kernel_rezone_neighbors);
#ifndef MINIMUM_PRECISION
      ezcl_kernel_release(kernel_rezone_one_double);
#endif
      ezcl_kernel_release(kernel_rezone_one_float);
      ezcl_kernel_release(kernel_copy_mpot_ghost_data);
      ezcl_kernel_release(kernel_set_boundary_refinement);
      terminate_kernel_2stage_sum();
      terminate_kernel_2stage_sum_int();
      if (! have_boundary){
        ezcl_kernel_release(kernel_count_BCs);
      }
#endif
#if defined(HAVE_J7) && defined(HAVE_MPI)
   if (parallel) mesh_memory.pfini();
#endif
}

int Mesh::rezone_count(vector<char_t> mpot, int &icount, int &jcount)
{
   int my_icount=0;
   int my_jcount=0;

#ifdef _OPENMP
#pragma omp parallel for reduction (+:my_jcount,my_icount)
#endif
   for (uint ic=0; ic<ncells; ++ic){
      if (mpot[ic] < 0) {
         if (celltype[ic] == REAL_CELL) {
            // remove all but cell that will remain to get count right when split
            // across processors
            if (! is_lower_left(i[ic],j[ic]) ) my_jcount--;
         } else {
            // either upper right or lower left will remain for boundary cells
            if (! (is_upper_right(i[ic],j[ic]) || is_lower_left(i[ic],j[ic]) ) ) my_jcount--;
         }
      }

      if (mpot[ic] > 0) {
         //printf("mpot[%d] = %d level %d levmx %d\n",ic,mpot[ic],level[ic],levmx);
         if (celltype[ic] == REAL_CELL){
            my_icount += 3;
         } else {
            my_icount ++;
         }
      }
   }
   icount = my_icount;
   jcount = my_jcount;

   return(icount+jcount);
}

int Mesh::rezone_count_threaded(vector<char_t> mpot, int &icount, int &jcount)
{
   static int my_icount=0;
   static int my_jcount=0;

#ifdef _OPENMP
#pragma omp for reduction (+:my_jcount,my_icount)
#endif
   for (uint ic=0; ic<ncells; ++ic){
      if (mpot[ic] < 0) {
         if (celltype[ic] == REAL_CELL) {
            // remove all but cell that will remain to get count right when split
            // across processors
            if (! is_lower_left(i[ic],j[ic]) ) my_jcount--;
         } else {
            // either upper right or lower left will remain for boundary cells
            if (! (is_upper_right(i[ic],j[ic]) || is_lower_left(i[ic],j[ic]) ) ) my_jcount--;
         }
      }

      if (mpot[ic] > 0) {
         //printf("mpot[%d] = %d level %d levmx %d\n",ic,mpot[ic],level[ic],levmx);
         if (celltype[ic] == REAL_CELL){
            my_icount += 3;
         } else {
            my_icount ++;
         }
      }
   }
   //printf("icount is %d\n",my_icount);
#ifdef _OPENMP
#pragma omp master
#endif
   {
      icount = my_icount;
      jcount = my_jcount;
   }
#ifdef _OPENMP
#pragma omp barrier
#endif

   return(icount+jcount);
}

#ifdef HAVE_OPENCL
void Mesh::gpu_rezone_count2(size_t block_size, size_t local_work_size, cl_mem dev_redscratch, cl_mem &dev_result)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

     /*
     __kernel void finish_reduction_count2_cl(
                       const    int   isize,      // 0
              __global          int  *redscratch, // 1
              __global          int  *result,     // 2
              __local           int  *tile)       // 3
     */
   ezcl_set_kernel_arg(kernel_reduction_count2, 0, sizeof(cl_int),  (void *)&block_size);
   ezcl_set_kernel_arg(kernel_reduction_count2, 1, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_reduction_count2, 2, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_reduction_count2, 3, local_work_size*sizeof(cl_int2),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_count2, 1, NULL, &local_work_size, &local_work_size, NULL);
}

void Mesh::gpu_rezone_count(size_t block_size, size_t local_work_size, cl_mem dev_redscratch, cl_mem &dev_result)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

     /*
     __kernel void finish_reduction_count_cl(
                       const    int   isize,      // 0
              __global          int  *redscratch, // 1
              __global          int  *result,     // 2
              __local           int  *tile)       // 3
     */
   ezcl_set_kernel_arg(kernel_reduction_count, 0, sizeof(cl_int),  (void *)&block_size);
   ezcl_set_kernel_arg(kernel_reduction_count, 1, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_reduction_count, 2, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_reduction_count, 3, local_work_size*sizeof(cl_int),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_count, 1, NULL, &local_work_size, &local_work_size, NULL);
}

void Mesh::gpu_rezone_scan(size_t block_size, size_t local_work_size, cl_mem dev_ioffset, cl_mem &dev_result)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

     /*
     __kernel void finish_reduction_scan_cl(
                       const    int   isize,    // 0
              __global          int  *ioffset,  // 1
              __global          int  *result,   // 2
              __local           int  *tile)     // 3
     */
   ezcl_set_kernel_arg(kernel_reduction_scan2, 0, sizeof(cl_int),  (void *)&block_size);
   ezcl_set_kernel_arg(kernel_reduction_scan2, 1, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_reduction_scan2, 2, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_reduction_scan2, 3, local_work_size*sizeof(cl_uint2),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_scan2, 1, NULL, &local_work_size, &local_work_size, NULL);
}
#endif

void Mesh::kdtree_setup()
{
   KDTree_Initialize(&tree);

   TBounds box;
   for (uint ic=0; ic<ncells; ic++) {
     box.min.x = x[ic];
     box.max.x = x[ic]+dx[ic];
     box.min.y = y[ic];
     box.max.y = y[ic]+dy[ic];
     KDTree_AddElement(&tree, &box);
   }
}

void Mesh::calc_spatial_coordinates(int ibase)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   x.clear();
   dx.clear();
   y.clear();
   dy.clear();
   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);

#ifdef _OPENMP
#pragma omp parallel
   {
#endif

   int lowerBounds, upperBounds;
   set_bounds(ncells);
   get_bounds(lowerBounds, upperBounds);

   if (have_boundary) {
      for (int ic = lowerBounds; ic < upperBounds; ic++) {
         uchar_t lev = level[ic];
         x[ic]  = xmin + (lev_deltax[lev] * (i[ic] - ibase));
         dx[ic] =        lev_deltax[lev];
         y[ic]  = ymin + (lev_deltay[lev] * (j[ic] - ibase));
         dy[ic] =        lev_deltay[lev];
      }
   } else {
      for (int ic = lowerBounds; ic < upperBounds; ic++) {
         uchar_t lev = level[ic];
         x[ic]  = xmin + (lev_deltax[lev] * (i[ic] - lev_ibegin[lev]));
         dx[ic] =        lev_deltax[lev];
         y[ic]  = ymin + (lev_deltay[lev] * (j[ic] - lev_jbegin[lev]));
         dy[ic] =        lev_deltay[lev];
      }
   }

#ifdef _OPENMP
#pragma omp master
#endif
   cpu_timers[MESH_TIMER_CALC_SPATIAL_COORDINATES] += cpu_timer_stop(tstart_cpu);

#ifdef _OPENMP
#pragma omp barrier
   } // end parallel region
#endif
}

#ifdef HAVE_OPENCL
void Mesh::gpu_calc_spatial_coordinates(cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_event calc_spatial_coordinates_event;

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t local_work_size = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;

// Only coded for base 0 and have boundary
//  Need:
//     xmin
//     ymin
//
//     lev_deltax -- dev_levdx
//     lev_deltay -- dev_levdy
//     x
//     dx
//     y
//     dy
//     level
//     i
//     j

   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  0, sizeof(cl_int),    (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  1, sizeof(cl_real_t), (void *)&xmin);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  2, sizeof(cl_real_t), (void *)&ymin);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  3, sizeof(cl_mem),    (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  4, sizeof(cl_mem),    (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  5, sizeof(cl_mem),    (void *)&dev_x);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  6, sizeof(cl_mem),    (void *)&dev_dx);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  7, sizeof(cl_mem),    (void *)&dev_y);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  8, sizeof(cl_mem),    (void *)&dev_dy);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  9, sizeof(cl_mem),    (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates, 10, sizeof(cl_mem),    (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates, 11, sizeof(cl_mem),    (void *)&dev_j);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_spatial_coordinates, 1, NULL, &global_work_size, &local_work_size, &calc_spatial_coordinates_event);

   ezcl_wait_for_events(1, &calc_spatial_coordinates_event);
   ezcl_event_release(calc_spatial_coordinates_event);

   gpu_timers[MESH_TIMER_CALC_SPATIAL_COORDINATES] += (long)(cpu_timer_stop(tstart_cpu) * 1.0e9);
}
#endif

void Mesh::calc_minmax(void)
{
   xmin=+1.0e30, ymin=+1.0e30, zmin=+1.0e30;

   for (uint ic=0; ic<ncells; ic++){
      if (x[ic] < xmin) xmin = x[ic];
   }
   for (uint ic=0; ic<ncells; ic++){
      if (y[ic] < ymin) ymin = y[ic];
   }
   if (ndim > TWO_DIMENSIONAL) {
      for (uint ic=0; ic<ncells; ic++){
         if (z[ic] < zmin) zmin = z[ic];
      }
   }

   xmax=-1.0e30, ymax=-1.0e30, zmax=-1.0e30;
   real_t xhigh, yhigh, zhigh;

   for (uint ic=0; ic<ncells; ic++){
      xhigh = x[ic]+dx[ic];
      if (xhigh > xmax) xmax = xhigh;
   }
   for (uint ic=0; ic<ncells; ic++){
      yhigh = y[ic]+dy[ic];
      if (yhigh > ymax) ymax = yhigh;
   }
   if (ndim > TWO_DIMENSIONAL) {
      for (uint ic=0; ic<ncells; ic++){
        zhigh = z[ic]+dz[ic];
        if (zhigh > zmax) zmax = zhigh;
      }
   }

#ifdef HAVE_MPI
   if (parallel) {
      real_t xmin_global,xmax_global,ymin_global,ymax_global;
      MPI_Allreduce(&xmin, &xmin_global, 1, MPI_REAL_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&xmax, &xmax_global, 1, MPI_REAL_T, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&ymin, &ymin_global, 1, MPI_REAL_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&ymax, &ymax_global, 1, MPI_REAL_T, MPI_MAX, MPI_COMM_WORLD);
      xmin = xmin_global;
      xmax = xmax_global;
      ymin = ymin_global;
      ymax = ymax_global;
   }
#endif

}
void Mesh::calc_centerminmax(void)
{
   xcentermin=+1.0e30, ycentermin=+1.0e30, zcentermin=+1.0e30;
   xcentermax=-1.0e30, ycentermax=-1.0e30, zcentermax=-1.0e30;
   real_t xmid, ymid, zmid;

   for (uint ic=0; ic<ncells; ic++){
      xmid = x[ic]+0.5*dx[ic];
      if (xmid < xcentermin) xcentermin = xmid;
      if (xmid > xcentermax) xcentermax = xmid;
   }
   for (uint ic=0; ic<ncells; ic++){
      ymid = y[ic]+0.5*dy[ic];
      if (ymid < ycentermin) ycentermin = ymid;
      if (ymid > ycentermax) ycentermax = ymid;
   }
   if (ndim > TWO_DIMENSIONAL) {
      for (uint ic=0; ic<ncells; ic++){
         zmid = z[ic]+0.5*dz[ic];
         if (zmid < zcentermin) zcentermin = zmid;
         if (zmid > zcentermax) zcentermax = zmid;
      }
   }

#ifdef HAVE_MPI
   if (parallel) {
      real_t xcentermin_global,xcentermax_global,ycentermin_global,ycentermax_global;
      MPI_Allreduce(&xcentermin, &xcentermin_global, 1, MPI_REAL_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&xcentermax, &xcentermax_global, 1, MPI_REAL_T, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&ycentermin, &ycentermin_global, 1, MPI_REAL_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&ycentermax, &ycentermax_global, 1, MPI_REAL_T, MPI_MAX, MPI_COMM_WORLD);
      xcentermin = xcentermin_global;
      xcentermax = xcentermax_global;
      ycentermin = ycentermin_global;
      ycentermax = ycentermax_global;
   }
#endif

}

void Mesh::rezone_all(int icount, int jcount, vector<char_t> mpot, int have_state, MallocPlus &state_memory)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   if (! do_rezone) {

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      index.clear();
      index.resize(ncells);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
      for (uint ic=0; ic<ncells; ic++){
         index[ic]=ic;
      }

#ifdef _OPENMP
#pragma omp master
#endif
      cpu_timers[MESH_TIMER_REZONE_ALL] += cpu_timer_stop(tstart_cpu);

   } else {

// sign for jcount is different in GPU and CPU code -- abs is a quick fix
   int add_ncells = icount - abs(jcount);

#ifdef _OPENMP
#pragma omp master
#endif
   cpu_counters[MESH_COUNTER_REZONE]++;

   static vector<int> celltype_save;

   static int new_ncells;

   static int *i_old, *j_old;
   static uchar_t *level_old;

   static int ifirst;
   static int ilast;
   static int jfirst;
   static int jlast;
   static uchar_t level_first;
   static uchar_t level_last;

   static vector<int> new_ic;

#ifdef _OPENMP
#pragma omp master
   {
#endif
      celltype_save.resize(ncells);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   if (have_state) {
#ifdef _OPENMP
#pragma omp for
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
      for (int ic = 0; ic < (int)ncells; ic++){
         celltype_save[ic] = celltype[ic];
      }
   }

#ifdef _OPENMP
#pragma omp master
   {
#endif
   new_ncells = ncells + add_ncells;
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

// int ref_entry_count = 0;
   if (have_state){
#ifdef _OPENMP
#pragma omp for
#endif
      for (uint ic=0; ic<ncells; ic++) {
//       if (mpot[ic] > 0) ref_entry_count++;
         if (mpot[ic] < 0) {
            // Normal cell coarsening
            if (is_lower_left(i[ic],j[ic]) ) mpot[ic] = -2;
            // Boundary cell case
            if (celltype[ic] != REAL_CELL && is_upper_right(i[ic],j[ic]) ) mpot[ic] = -3;
         }
      }
   }

   //  Initialize new variables
// int *i_old, *j_old, *level_old;

   int flags = RESTART_DATA;
#ifdef HAVE_J7
   if (parallel) flags = LOAD_BALANCE_MEMORY;
#endif

#ifdef _OPENMP
#pragma omp master
   {
#endif
   i_old     = (int *)mesh_memory.memory_malloc(new_ncells, sizeof(int), "i_old",     flags);
   j_old     = (int *)mesh_memory.memory_malloc(new_ncells, sizeof(int), "j_old",     flags);
   level_old = (uchar_t *)mesh_memory.memory_malloc(new_ncells, sizeof(uchar_t), "level_old", flags);
   mesh_memory.memory_swap(&i,     &i_old);
   mesh_memory.memory_swap(&j,     &j_old);
   mesh_memory.memory_swap(&level, &level_old);

   index.clear();
   // allocate to ncells (oldsize) because index is a map from old cells to new cells
   // loop below is for (ic = 0; ic < ncells; ic++) and within it index[ic] = nc;
   // If set to new_ncells, will write out-of-bounds when number of cells decreas
   // mpirun -n 3 valgrind ./clamr_mpionly -t 391 -i 1 fails at cycle 390
   index.resize(ncells);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   static vector<int> order; //  Vector of refined mesh traversal order; set to -1 to indicate errors.
   //
   //vector<int>  invorder(4, -1); //  Vector mapping location from base index.

   //int ref_entry = 0;

#ifdef _OPENMP
#pragma omp master
   {
#endif
   //  Insert new cells into the mesh at the point of refinement.
   order.resize(4,    -1); //  Vector of refined mesh traversal order; set to -1 to indicate errors.

   ifirst      = 0;
   ilast       = 0;
   jfirst      = 0;
   jlast       = 0;
   level_first = 0;
   level_last  = 0;

   if (parallel) {
#ifdef HAVE_MPI
      MPI_Request req[12];
      MPI_Status status[12];

      static int prev     = MPI_PROC_NULL;
      static int next     = MPI_PROC_NULL;

      if (mype != 0)         prev = mype-1;
      if (mype < numpe - 1)  next = mype+1;

      MPI_Isend(&i_old[ncells-1],     1,MPI_INT,next,1,MPI_COMM_WORLD,req+0);
      MPI_Irecv(&ifirst,              1,MPI_INT,prev,1,MPI_COMM_WORLD,req+1);

      MPI_Isend(&i_old[0],            1,MPI_INT,prev,1,MPI_COMM_WORLD,req+2);
      MPI_Irecv(&ilast,               1,MPI_INT,next,1,MPI_COMM_WORLD,req+3);

      MPI_Isend(&j_old[ncells-1],     1,MPI_INT,next,1,MPI_COMM_WORLD,req+4);
      MPI_Irecv(&jfirst,              1,MPI_INT,prev,1,MPI_COMM_WORLD,req+5);

      MPI_Isend(&j_old[0],            1,MPI_INT,prev,1,MPI_COMM_WORLD,req+6);
      MPI_Irecv(&jlast,               1,MPI_INT,next,1,MPI_COMM_WORLD,req+7);

      MPI_Isend(&level_old[ncells-1], 1,MPI_UCHAR_T,next,1,MPI_COMM_WORLD,req+8);
      MPI_Irecv(&level_first,         1,MPI_UCHAR_T,prev,1,MPI_COMM_WORLD,req+9);

      MPI_Isend(&level_old[0],        1,MPI_UCHAR_T,prev,1,MPI_COMM_WORLD,req+10);
      MPI_Irecv(&level_last,          1,MPI_UCHAR_T,next,1,MPI_COMM_WORLD,req+11);

      MPI_Waitall(12, req, status);
#endif
   }

#ifdef _OPENMP
   }
#pragma omp barrier
#endif

#ifdef REZONE_NO_OPTIMIZATION
   vector<int>  invorder(4, -1); //  Vector mapping location from base index.
   for (int ic = 0, nc = 0; ic < (int)ncells; ic++)
   {
      if (mpot[ic] == 0 || mpot[ic] == -1000000)
      {  //  No change is needed; copy the old cell straight to the new mesh at this location.
         index[ic] = nc;
         i[nc]     = i_old[ic];
         j[nc]     = j_old[ic];
         level[nc] = level_old[ic];
         nc++;
      } //  Complete no change needed.
      
      else if (mpot[ic] < 0)
      {  //  Coarsening is needed; remove this cell and the other three and replace them with one.
         index[ic] = nc;
         if (mpot[ic] <= -2) {
            //printf("                     %d: DEBUG -- coarsening cell %d nc %d\n",mype,ic,nc);
            i[nc] = i_old[ic]/2;
            j[nc] = j_old[ic]/2;
            level[nc] = level_old[ic] - 1;
            nc++;
         }
      } //  Coarsening complete.
      
      else if (mpot[ic] > 0)
      {  //  Refinement is needed; insert four cells where once was one.
         index[ic] = nc;
         if (celltype[ic] == REAL_CELL)
         {  
            set_refinement_order(&order[0], ic, ifirst, ilast, jfirst, jlast,
                                 level_first, level_last, i_old, j_old, level_old);

            //  Create the cells in the correct order and orientation.
            for (int ii = 0; ii < 4; ii++)
            {  level[nc] = level_old[ic] + 1;
               switch (order[ii])
               {  case SW:
                     // lower left
                     invorder[SW] = ii;
                     i[nc]     = i_old[ic]*2;
                     j[nc]     = j_old[ic]*2;
                     nc++;
                     break;
                     
                  case SE:
                     // lower right
                     invorder[SE] = ii;
                     i[nc]     = i_old[ic]*2 + 1;
                     j[nc]     = j_old[ic]*2;
                     nc++;
                     break;
                     
                  case NW:
                     // upper left
                     invorder[NW] = ii;
                     i[nc]     = i_old[ic]*2;
                     j[nc]     = j_old[ic]*2 + 1;
                     nc++;
                     break;
                     
                  case NE:
                     // upper right
                     invorder[NE] = ii;
                     i[nc]     = i_old[ic]*2 + 1;
                     j[nc]     = j_old[ic]*2 + 1;
                     nc++;
                     break; } } //  Complete cell refinement.
         }  //  Complete real cell refinement.
         
         else if (celltype[ic] == LEFT_BOUNDARY) {
            // lower
            i[nc]  = i_old[ic]*2 + 1;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;
            
            // upper
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
         else if (celltype[ic] == RIGHT_BOUNDARY) {
            // lower
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;
            
            // upper
            i[nc] = i_old[ic]*2;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
         else if (celltype[ic] == BOTTOM_BOUNDARY) {
            // left
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
            
            // right
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
         else if (celltype[ic] == TOP_BOUNDARY) {
            // right
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;

            // left
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
      } //  Complete refinement needed.
   } //  Complete addition of new cells to the mesh.

   mesh_memory.memory_delete(i_old);
   mesh_memory.memory_delete(j_old);
   mesh_memory.memory_delete(level_old);

   calc_celltype(new_ncells);

   if (have_state){
      flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);
      MallocPlus state_memory_old = state_memory;
      malloc_plus_memory_entry *memory_item;

      for (memory_item = state_memory_old.memory_entry_by_name_begin();
           memory_item != state_memory_old.memory_entry_by_name_end();
           memory_item = state_memory_old.memory_entry_by_name_next() ) {
         //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

         if ( (memory_item->mem_flags & REZONE_DATA) == 0) continue;

         if (memory_item->mem_elsize == 8) {
            double *state_temp_double = (double *)state_memory.memory_malloc(new_ncells, sizeof(double),
                                                                             "state_temp_double", memory_item->mem_flags);

            double *mem_ptr_double = (double *)memory_item->mem_ptr;

            //ref_entry = 0;
            for (int ic=0, nc=0; ic<(int)ncells; ic++) {

               if (mpot[ic] == 0) {
                  state_temp_double[nc] = mem_ptr_double[ic];
                  nc++;
               } else if (mpot[ic] < 0){
                  if (mpot[ic] == -2) {
                     int nr = nrht[ic];
                     int nt = ntop[ic];
                     int nrt = nrht[nt];
                     state_temp_double[nc] = (mem_ptr_double[ic] + mem_ptr_double[nr] +
                                              mem_ptr_double[nt] + mem_ptr_double[nrt])*0.25;
                     nc++;
                  }
                  if (mpot[ic] == -3) {
                     int nl = nlft[ic];
                     int nb = nbot[ic];
                     int nlb = nlft[nb];
                     state_temp_double[nc] = (mem_ptr_double[ic] + mem_ptr_double[nl] +
                                              mem_ptr_double[nb] + mem_ptr_double[nlb])*0.25;
                     nc++;
                  }
               } else if (mpot[ic] > 0){
                  // lower left
                  state_temp_double[nc] = mem_ptr_double[ic];
                  nc++;

                  // lower right
                  state_temp_double[nc] = mem_ptr_double[ic];
                  nc++;

                  if (celltype_save[ic] == REAL_CELL){
                     // upper left
                     state_temp_double[nc] = mem_ptr_double[ic];
                     nc++;

                     // upper right
                     state_temp_double[nc] = mem_ptr_double[ic];
                     nc++;
                  }
               }
            }

            state_memory.memory_replace(mem_ptr_double, state_temp_double);
         } else if (memory_item->mem_elsize == 4) {
            float *state_temp_float = (float *)state_memory.memory_malloc(new_ncells, sizeof(float),
                                                                          "state_temp_float", memory_item->mem_flags);

            float *mem_ptr_float = (float *)memory_item->mem_ptr;

            for (int ic=0, nc=0; ic<(int)ncells; ic++) {

               if (mpot[ic] == 0) {
                  state_temp_float[nc] = mem_ptr_float[ic];
                  nc++;
               } else if (mpot[ic] < 0){
                  if (mpot[ic] == -2) {
                     int nr = nrht[ic];
                     int nt = ntop[ic];
                     int nrt = nrht[nt];
                     state_temp_float[nc] = (mem_ptr_float[ic] + mem_ptr_float[nr] +
                                             mem_ptr_float[nt] + mem_ptr_float[nrt])*0.25;
                     nc++;
                  }
                  if (mpot[ic] == -3) {
                     int nl = nlft[ic];
                     int nb = nbot[ic];
                     int nlb = nlft[nb];
                     state_temp_float[nc] = (mem_ptr_float[ic] + mem_ptr_float[nl] +
                                             mem_ptr_float[nb] + mem_ptr_float[nlb])*0.25;
                     nc++;
                  }
               } else if (mpot[ic] > 0){
                  // lower left
                  state_temp_float[nc] = mem_ptr_float[ic];
                  nc++;

                  // lower right
                  state_temp_float[nc] = mem_ptr_float[ic];
                  nc++;

                  if (celltype_save[ic] == REAL_CELL){
                     // upper left
                     state_temp_float[nc] = mem_ptr_float[ic];
                     nc++;

                     // upper right
                     state_temp_float[nc] = mem_ptr_float[ic];
                     nc++;
                  }
               }
            }

            state_memory.memory_replace(mem_ptr_float, state_temp_float);
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            half *state_temp_half = (half *)state_memory.memory_malloc(new_ncells, sizeof(half),
                                                                          "state_temp_half", memory_item->mem_flags);

            half *mem_ptr_half = (half *)memory_item->mem_ptr;

            for (int ic=0, nc=0; ic<(int)ncells; ic++) {

               if (mpot[ic] == 0) {
                  state_temp_half[nc] = mem_ptr_half[ic];
                  nc++;
               } else if (mpot[ic] < 0){
                  if (mpot[ic] == -2) {
                     int nr = nrht[ic];
                     int nt = ntop[ic];
                     int nrt = nrht[nt];
                     state_temp_half[nc] = (mem_ptr_half[ic] + mem_ptr_half[nr] +
                                             mem_ptr_half[nt] + mem_ptr_half[nrt])*0.25;
                     nc++;
                  }
                  if (mpot[ic] == -3) {
                     int nl = nlft[ic];
                     int nb = nbot[ic];
                     int nlb = nlft[nb];
                     state_temp_half[nc] = (mem_ptr_half[ic] + mem_ptr_half[nl] +
                                             mem_ptr_half[nb] + mem_ptr_half[nlb])*0.25;
                     nc++;
                  }
               } else if (mpot[ic] > 0){
                  // lower left
                  state_temp_half[nc] = mem_ptr_half[ic];
                  nc++;

                  // lower right
                  state_temp_half[nc] = mem_ptr_half[ic];
                  nc++;

                  if (celltype_save[ic] == REAL_CELL){
                     // upper left
                     state_temp_half[nc] = mem_ptr_half[ic];
                     nc++;

                     // upper right
                     state_temp_half[nc] = mem_ptr_half[ic];
                     nc++;
                  }
               }
            }

            state_memory.memory_replace(mem_ptr_half, state_temp_half);
#endif
         }
      }
   }
#else
   // Data parallel optimizations for thread parallel -- slows down serial
   // code by about 25%
   static vector<int> add_count;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      add_count.resize(ncells);
      new_ic.resize(ncells+1);
#ifdef _OPENMP
   } // end master region
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
      for (int ic = 0; ic < (int)ncells; ic++){
          //if (new_ic[ic] == 284) printf("%d\n", ic);
         if (mpot[ic] == 0) {
            add_count[ic] = 1;
         } else if (mpot[ic] < 0) {
            if (mpot[ic] <= -2){
               add_count[ic] = 1;
            } else {
               add_count[ic] = 0;
            }
         } else if (mpot[ic] > 0) {
            if (celltype[ic] != REAL_CELL) {
               add_count[ic] = 2;
            } else {
               add_count[ic] = 4;
            }
         }
      }

#ifdef _OPENMP
#pragma omp barrier
#endif
      scan (&add_count[0], &new_ic[0], ncells);
#ifdef _OPENMP
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
   for (int ic = 0; ic < (int)ncells; ic++) {
   vector<int>  invorder(4, -1); //  Vector mapping location from base index.
      int nc = new_ic[ic];
          //if (new_ic[ic] == 284) printf("%d\n", ic);
      if (mpot[ic] == 0)
      {  //  No change is needed; copy the old cell straight to the new mesh at this location.
         index[ic] = nc;
         i[nc]     = i_old[ic];
         j[nc]     = j_old[ic];
         level[nc] = level_old[ic];
      //if (ic == 296) printf("%d\n", mpot[ic]);
       //if (ic == 296) printf("%d %d %d\n", i_old[ic], j_old[ic], level_old[ic]);
       //if (ic == 296) printf("%d %d %d %d\n", nc, i[nc], j[nc], level[nc]);
         
      } //  Complete no change needed.

      else if (mpot[ic] < 0)
      {  //  Coarsening is needed; remove this cell and the other three and replace them with one.
         index[ic] = nc;
         if (mpot[ic] <= -2) {
            //printf("                     %d: DEBUG -- coarsening cell %d nc %d\n",mype,ic,nc);
            i[nc] = i_old[ic]/2;
            j[nc] = j_old[ic]/2;
            level[nc] = level_old[ic] - 1;
         }
      } //  Coarsening complete.

      else if (mpot[ic] > 0)
      {  //  Refinement is needed; insert four cells where once was one.
         index[ic] = nc;
         if (celltype[ic] == REAL_CELL)
         {  
            int order[4];
            set_refinement_order(&order[0], ic, ifirst, ilast, jfirst, jlast,
                                 level_first, level_last, i_old, j_old, level_old);

            //  Create the cells in the correct order and orientation.
            for (int ii = 0; ii < 4; ii++) {
               level[nc] = level_old[ic] + 1;
               switch (order[ii]) {
                  case SW:
                     // lower left
                     invorder[SW] = ii;
                     i[nc]     = i_old[ic]*2;
                     j[nc]     = j_old[ic]*2;
                     nc++;
                     break;
                     
                  case SE:
                     // lower right
                     invorder[SE] = ii;
                     i[nc]     = i_old[ic]*2 + 1;
                     j[nc]     = j_old[ic]*2;
                     nc++;
                     break;
                     
                  case NW:
                     // upper left
                     invorder[NW] = ii;
                     i[nc]     = i_old[ic]*2;
                     j[nc]     = j_old[ic]*2 + 1;
                     nc++;
                     break;
                     
                  case NE:
                     // upper right
                     invorder[NE] = ii;
                     i[nc]     = i_old[ic]*2 + 1;
                     j[nc]     = j_old[ic]*2 + 1;
                     nc++;
                     break;
                  }
               } //  Complete cell refinement.
         }  //  Complete real cell refinement.
         
         else if (celltype[ic] == LEFT_BOUNDARY) {
            // lower
            i[nc]  = i_old[ic]*2 + 1;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;
            
            // upper
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
         else if (celltype[ic] == RIGHT_BOUNDARY) {
            // lower
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;
            
            // upper
            i[nc] = i_old[ic]*2;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
         else if (celltype[ic] == BOTTOM_BOUNDARY) {
            // left
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
            
            // right
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
         else if (celltype[ic] == TOP_BOUNDARY) {
            // right
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;

            // left
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            nc++;
         }
      } //  Complete refinement needed.
   } //  Complete addition of new cells to the mesh.

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
       //printf("%d %d %d\n", i[284], j[284], level[284]);
   mesh_memory.memory_delete(i_old);
   mesh_memory.memory_delete(j_old);
   mesh_memory.memory_delete(level_old);
#ifdef _OPENMP
   } // end master region
#endif

   calc_celltype_threaded(new_ncells);

   if (have_state){

      static MallocPlus state_memory_old;
      static malloc_plus_memory_entry *memory_begin;
      static malloc_plus_memory_entry *memory_end;
      static malloc_plus_memory_entry *memory_next;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {
#endif
      state_memory_old = state_memory;

      memory_begin = state_memory_old.memory_entry_by_name_begin();
      memory_end   = state_memory_old.memory_entry_by_name_end();
#ifdef _OPENMP
      } // end master region
#pragma omp barrier
#endif

      for (malloc_plus_memory_entry *memory_item = memory_begin;
           memory_item != memory_end;
           memory_item = memory_next ) {

         //ref_entry = 0;
         //printf("DEBUG -- memory_item->mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

         if ( (memory_item->mem_flags & REZONE_DATA) == 0) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            memory_next = state_memory_old.memory_entry_by_name_next();
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif
            continue;
         }

         if (memory_item->mem_elsize == 8) {

            static double *state_temp_double, *mem_ptr_double;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               state_temp_double = (double *)state_memory.memory_malloc(new_ncells, sizeof(double),
                                                                                "state_temp_double", memory_item->mem_flags);
               mem_ptr_double = (double *)memory_item->mem_ptr;
#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif

            //ref_entry = 0;
#ifdef _OPENMP
#pragma omp for
#endif
            for (int ic=0; ic<(int)ncells; ic++) {

               int nc = new_ic[ic];
               if (mpot[ic] == 0) {
                  state_temp_double[nc] = mem_ptr_double[ic];
               } else if (mpot[ic] < 0){
                  if (mpot[ic] == -2) {
                     int nr = nrht[ic];
                     int nt = ntop[ic];
                     int nrt = nrht[nt];
                     state_temp_double[nc] = (mem_ptr_double[ic] + mem_ptr_double[nr] +
                                              mem_ptr_double[nt] + mem_ptr_double[nrt])*0.25;
                  }
                  if (mpot[ic] == -3) {
                     int nl = nlft[ic];
                     int nb = nbot[ic];
                     int nlb = nlft[nb];
                     state_temp_double[nc] = (mem_ptr_double[ic] + mem_ptr_double[nl] +
                                              mem_ptr_double[nb] + mem_ptr_double[nlb])*0.25;
                  }
               } else if (mpot[ic] > 0){
                  // lower left
                  state_temp_double[nc] = mem_ptr_double[ic];
                  nc++;

                  // lower right
                  state_temp_double[nc] = mem_ptr_double[ic];
                  nc++;

                  if (celltype_save[ic] == REAL_CELL){
                     // upper left
                     state_temp_double[nc] = mem_ptr_double[ic];
                     nc++;

                     // upper right
                     state_temp_double[nc] = mem_ptr_double[ic];
                     nc++;
                  }
               }
            } // end cell loop

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
            state_memory.memory_replace(mem_ptr_double, state_temp_double);
#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif


         } else if (memory_item->mem_elsize == 4) {

            static float *state_temp_float, *mem_ptr_float;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               state_temp_float = (float *)state_memory.memory_malloc(new_ncells, sizeof(float),
                                                                             "state_temp_float", memory_item->mem_flags);
               mem_ptr_float = (float *)memory_item->mem_ptr;
#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
            for (int ic=0; ic<(int)ncells; ic++) {

               int nc = new_ic[ic];
               if (mpot[ic] == 0) {
                  state_temp_float[nc] = mem_ptr_float[ic];
               } else if (mpot[ic] < 0){
                  if (mpot[ic] == -2) {
                     int nr = nrht[ic];
                     int nt = ntop[ic];
                     int nrt = nrht[nt];
                     state_temp_float[nc] = (mem_ptr_float[ic] + mem_ptr_float[nr] +
                                             mem_ptr_float[nt] + mem_ptr_float[nrt])*0.25;
                  }
                  if (mpot[ic] == -3) {
                     int nl = nlft[ic];
                     int nb = nbot[ic];
                     int nlb = nlft[nb];
                     state_temp_float[nc] = (mem_ptr_float[ic] + mem_ptr_float[nl] +
                                             mem_ptr_float[nb] + mem_ptr_float[nlb])*0.25;
                  }
               } else if (mpot[ic] > 0){
                  // lower left
                  state_temp_float[nc] = mem_ptr_float[ic];
                  nc++;

                  // lower right
                  state_temp_float[nc] = mem_ptr_float[ic];
                  nc++;

                  if (celltype_save[ic] == REAL_CELL){
                     // upper left
                     state_temp_float[nc] = mem_ptr_float[ic];
                     nc++;

                     // upper right
                     state_temp_float[nc] = mem_ptr_float[ic];
                     nc++;
                  }
               }
            } // end cell loop

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               state_memory.memory_replace(mem_ptr_float, state_temp_float);
#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {

            static half *state_temp_half, *mem_ptr_half;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               state_temp_half = (half *)state_memory.memory_malloc(new_ncells, sizeof(half),
                                                                             "state_temp_half", memory_item->mem_flags);
               mem_ptr_half = (half *)memory_item->mem_ptr;
#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
            for (int ic=0; ic<(int)ncells; ic++) {

               int nc = new_ic[ic];
               if (mpot[ic] == 0) {
                  state_temp_half[nc] = mem_ptr_half[ic];
               } else if (mpot[ic] < 0){
                  if (mpot[ic] == -2) {
                     int nr = nrht[ic];
                     int nt = ntop[ic];
                     int nrt = nrht[nt];
                     state_temp_half[nc] = (mem_ptr_half[ic] + mem_ptr_half[nr] +
                                             mem_ptr_half[nt] + mem_ptr_half[nrt])*0.25;
                  }
                  if (mpot[ic] == -3) {
                     int nl = nlft[ic];
                     int nb = nbot[ic];
                     int nlb = nlft[nb];
                     state_temp_half[nc] = (mem_ptr_half[ic] + mem_ptr_half[nl] +
                                             mem_ptr_half[nb] + mem_ptr_half[nlb])*0.25;
                  }
               } else if (mpot[ic] > 0){
                  // lower left
                  state_temp_half[nc] = mem_ptr_half[ic];
                  nc++;

                  // lower right
                  state_temp_half[nc] = mem_ptr_half[ic];
                  nc++;

                  if (celltype_save[ic] == REAL_CELL){
                     // upper left
                     state_temp_half[nc] = mem_ptr_half[ic];
                     nc++;

                     // upper right
                     state_temp_half[nc] = mem_ptr_half[ic];
                     nc++;
                  }
               }
            } // end cell loop

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               state_memory.memory_replace(mem_ptr_half, state_temp_half);
#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif

            // end half_precision
#endif
         } // mem elem size 2 bytes

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         memory_next = state_memory_old.memory_entry_by_name_next();
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif

      } // memory item iteration

   } // if have state
   // End of data parallel optimizations
#endif

   if (neighbor_remap) {
      int flags = 0;
      static int *nlft_old, *nrht_old, *nbot_old, *ntop_old;
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {
#endif
      nlft_old     = (int *)mesh_memory.memory_malloc(new_ncells, sizeof(int), "nlft_old",  flags);
      nrht_old     = (int *)mesh_memory.memory_malloc(new_ncells, sizeof(int), "nrht_old",  flags);
      nbot_old     = (int *)mesh_memory.memory_malloc(new_ncells, sizeof(int), "nbot_old",  flags);
      ntop_old     = (int *)mesh_memory.memory_malloc(new_ncells, sizeof(int), "ntop_old",  flags);
#ifdef _OPENMP
      } // end master region
#pragma omp barrier
#endif
      flags = RESTART_DATA;

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
      for (int ic = 0; ic < new_ncells; ic++){
         nlft_old[ic] = -1;
         nrht_old[ic] = -1;
         nbot_old[ic] = -1;
         ntop_old[ic] = -1;
      }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {
#endif
      mesh_memory.memory_swap(&nlft,  &nlft_old);
      mesh_memory.memory_swap(&nrht,  &nrht_old);
      mesh_memory.memory_swap(&nbot,  &nbot_old);
      mesh_memory.memory_swap(&ntop,  &ntop_old);
#ifdef _OPENMP
      } // end master region
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
      for (int ic = 0; ic < (int)ncells; ic++){
         int nc = index[ic];

         if (mpot[ic] == 0){
            if (nlft_old[ic] < (int)ncells && nlft_old[ic] >= 0){
               nlft[nc] = (mpot[nlft_old[ic]] == 0) ? index[nlft_old[ic]] : -1;
            }
            if (nrht_old[ic] < (int)ncells && nrht_old[ic] >= 0){
               nrht[nc] = (mpot[nrht_old[ic]] == 0) ? index[nrht_old[ic]] : -1;
            }
            if (nbot_old[ic] < (int)ncells && nbot_old[ic] >= 0){
               nbot[nc] = (mpot[nbot_old[ic]] == 0) ? index[nbot_old[ic]] : -1;
            }
            if (ntop_old[ic] < (int)ncells && ntop_old[ic] >= 0){
               ntop[nc] = (mpot[ntop_old[ic]] == 0) ? index[ntop_old[ic]] : -1;
            }
         } else if (mpot[ic] <= -2) {
            nlft[nc]  = -1;
            nrht[nc]  = -1;
            nbot[nc]  = -1;
            ntop[nc]  = -1;
         } else if (mpot[ic] > 0){
            nlft[nc]    = -1;
            nlft[nc+1]  = -1;
            nrht[nc]    = -1;
            nrht[nc+1]  = -1;
            nbot[nc]    = -1;
            nbot[nc+1]  = -1;
            ntop[nc]    = -1;
            ntop[nc+1]  = -1;
            if (celltype[nc] == REAL_CELL){
               nlft[nc+2]  = -1;
               nlft[nc+3]  = -1;
               nrht[nc+2]  = -1;
               nrht[nc+3]  = -1;
               nbot[nc+2]  = -1;
               nbot[nc+3]  = -1;
               ntop[nc+2]  = -1;
               ntop[nc+3]  = -1;
            }
         }
         if (mpot[ic] > 0){
            nc++;
            switch(celltype[nc]){
            case LEFT_BOUNDARY:
               nlft[nc] = nc;
               break;
            case RIGHT_BOUNDARY:
               nrht[nc] = nc;
               break;
            case BOTTOM_BOUNDARY:
               nbot[nc] = nc;
               break;
            case TOP_BOUNDARY:
               ntop[nc] = nc;
               break;
            }
         }
      }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {
#endif
      nlft_old = (int *)mesh_memory.memory_delete(nlft_old);
      nrht_old = (int *)mesh_memory.memory_delete(nrht_old);
      nbot_old = (int *)mesh_memory.memory_delete(nbot_old);
      ntop_old = (int *)mesh_memory.memory_delete(ntop_old);
#ifdef _OPENMP
      } // end master region
#pragma omp barrier
#endif
   } else {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {
#endif
      nlft = (int *)mesh_memory.memory_delete(nlft);
      nrht = (int *)mesh_memory.memory_delete(nrht);
      nbot = (int *)mesh_memory.memory_delete(nbot);
      ntop = (int *)mesh_memory.memory_delete(ntop);
#ifdef _OPENMP
      } // end master region
#pragma omp barrier
#endif
   }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
   //ncells = nc;

#ifdef HAVE_MPI
   if (parallel) {
      MPI_Allgather(&new_ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);

      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }  
      noffset=ndispl[mype];
      ncells_global = ndispl[numpe-1]+nsizes[numpe-1];
   }  
#endif

   cpu_timers[MESH_TIMER_REZONE_ALL] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   } // end master region
#pragma omp barrier
#endif

   } // if do_rezone

}

#ifdef HAVE_OPENCL
void Mesh::gpu_rezone_all(int icount, int jcount, cl_mem &dev_mpot, MallocPlus &gpu_state_memory)
{
   if (! gpu_do_rezone) return;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   gpu_counters[MESH_COUNTER_REZONE]++;

   cl_command_queue command_queue = ezcl_get_command_queue();

   assert(dev_mpot);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);
   assert(dev_celltype);
   assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   int add_ncells = icount - jcount;

// int global_icount = icount;
// int global_jcount = jcount;

   size_t old_ncells = ncells;
   size_t new_ncells = ncells + add_ncells;

#ifdef HAVE_MPI
   //int global_add_ncells = add_ncells;

// if (parallel) {
//    int count[2], global_count[2];
//    count[0] = icount;
//    count[1] = jcount;
//    MPI_Allreduce(&count, &global_count, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//    global_icount = global_count[0];
//    global_jcount = global_count[1];
//    //global_add_ncells = global_icount + global_jcount;
// }
#endif

   int ifirst      = 0;
   int ilast       = 0;
   int jfirst      = 0;
   int jlast       = 0;
   uchar_t level_first = 0;
   uchar_t level_last  = 0;

#ifdef HAVE_MPI
   if (numpe > 1) {
      int i_tmp_first, i_tmp_last;
      int j_tmp_first, j_tmp_last;
      uchar_t level_tmp_first, level_tmp_last;

      ezcl_enqueue_read_buffer(command_queue,  dev_i,     CL_FALSE, 0,                             1*sizeof(cl_int), &i_tmp_first,     NULL);
      ezcl_enqueue_read_buffer(command_queue,  dev_j,     CL_FALSE, 0,                             1*sizeof(cl_int), &j_tmp_first,     NULL);
      ezcl_enqueue_read_buffer(command_queue,  dev_level, CL_FALSE, 0,                             1*sizeof(cl_uchar_t), &level_tmp_first, NULL);
      ezcl_enqueue_read_buffer(command_queue,  dev_i,     CL_FALSE, (old_ncells-1)*sizeof(cl_int), 1*sizeof(cl_int), &i_tmp_last,      NULL);
      ezcl_enqueue_read_buffer(command_queue,  dev_j,     CL_FALSE, (old_ncells-1)*sizeof(cl_int), 1*sizeof(cl_int), &j_tmp_last,      NULL);
      ezcl_enqueue_read_buffer(command_queue,  dev_level, CL_TRUE,  (old_ncells-1)*sizeof(cl_uchar_t), 1*sizeof(cl_uchar_t), &level_tmp_last,  NULL);

      MPI_Request req[12];
      MPI_Status status[12];

      static int prev     = MPI_PROC_NULL;
      static int next     = MPI_PROC_NULL;

      if (mype != 0)         prev = mype-1;
      if (mype < numpe - 1)  next = mype+1;

      MPI_Isend(&i_tmp_last,      1,MPI_INT,next,1,MPI_COMM_WORLD,req+0);
      MPI_Irecv(&ifirst,          1,MPI_INT,prev,1,MPI_COMM_WORLD,req+1);

      MPI_Isend(&i_tmp_first,     1,MPI_INT,prev,1,MPI_COMM_WORLD,req+2);
      MPI_Irecv(&ilast,           1,MPI_INT,next,1,MPI_COMM_WORLD,req+3);

      MPI_Isend(&j_tmp_last,      1,MPI_INT,next,1,MPI_COMM_WORLD,req+4);
      MPI_Irecv(&jfirst,          1,MPI_INT,prev,1,MPI_COMM_WORLD,req+5);

      MPI_Isend(&j_tmp_first,     1,MPI_INT,prev,1,MPI_COMM_WORLD,req+6);
      MPI_Irecv(&jlast,           1,MPI_INT,next,1,MPI_COMM_WORLD,req+7);

      MPI_Isend(&level_tmp_last,  1,MPI_UCHAR_T,next,1,MPI_COMM_WORLD,req+8);
      MPI_Irecv(&level_first,     1,MPI_UCHAR_T,prev,1,MPI_COMM_WORLD,req+9);

      MPI_Isend(&level_tmp_first, 1,MPI_UCHAR_T,prev,1,MPI_COMM_WORLD,req+10);
      MPI_Irecv(&level_last,      1,MPI_UCHAR_T,next,1,MPI_COMM_WORLD,req+11);

      MPI_Waitall(12, req, status);
   }
#endif

/*
   if (new_ncells != old_ncells){
      ncells = new_ncells;
   }
*/

   size_t mem_request = (int)((float)new_ncells*mem_factor);
   cl_mem dev_celltype_new = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_new"), &mem_request, sizeof(cl_char_t), CL_MEM_READ_WRITE, 0);
   cl_mem dev_level_new    = ezcl_malloc(NULL, const_cast<char *>("dev_level_new"),    &mem_request, sizeof(cl_uchar_t), CL_MEM_READ_WRITE, 0);
   cl_mem dev_i_new        = ezcl_malloc(NULL, const_cast<char *>("dev_i_new"),        &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
   cl_mem dev_j_new        = ezcl_malloc(NULL, const_cast<char *>("dev_j_new"),        &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

   cl_mem dev_ijadd;

   vector<int>ijadd(6);
   if (numpe > 1) {
      ijadd[0] = ifirst;
      ijadd[1] = ilast;
      ijadd[2] = jfirst;
      ijadd[3] = jlast;
      ijadd[4] = level_first;
      ijadd[5] = level_last;
   }

   size_t six = 6;
   dev_ijadd = ezcl_malloc(NULL, const_cast<char *>("dev_ijadd"), &six, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
   ezcl_enqueue_write_buffer(command_queue, dev_ijadd, CL_TRUE, 0, 6*sizeof(cl_int), (void*)&ijadd[0], NULL);

   cl_mem dev_indexoffset = ezcl_malloc(NULL, const_cast<char *>("dev_indexoffset"), &old_ncells, sizeof(cl_uint), CL_MEM_READ_WRITE, 0);

   int stencil = 0;
   if (localStencil) stencil = 1;

   size_t local_work_size = 128;
   size_t global_work_size = ((old_ncells+local_work_size - 1) /local_work_size) * local_work_size;

   ezcl_set_kernel_arg(kernel_rezone_all, 0,  sizeof(cl_int),  (void *)&old_ncells);
   ezcl_set_kernel_arg(kernel_rezone_all, 1,  sizeof(cl_int),  (void *)&stencil);
   ezcl_set_kernel_arg(kernel_rezone_all, 2,  sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_rezone_all, 3,  sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_rezone_all, 4,  sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_rezone_all, 5,  sizeof(cl_mem),  (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_rezone_all, 6,  sizeof(cl_mem),  (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_rezone_all, 7,  sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_rezone_all, 8,  sizeof(cl_mem),  (void *)&dev_level_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 9,  sizeof(cl_mem),  (void *)&dev_i_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 10, sizeof(cl_mem),  (void *)&dev_j_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 11, sizeof(cl_mem),  (void *)&dev_celltype_new);
   ezcl_set_kernel_arg(kernel_rezone_all, 12, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_rezone_all, 13, sizeof(cl_mem),  (void *)&dev_indexoffset);
   ezcl_set_kernel_arg(kernel_rezone_all, 14, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_rezone_all, 15, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_rezone_all, 16, sizeof(cl_mem),  (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_rezone_all, 17, sizeof(cl_mem),  (void *)&dev_ijadd);
   ezcl_set_kernel_arg(kernel_rezone_all, 18, local_work_size * sizeof(cl_uint), NULL);
   //ezcl_set_kernel_arg(kernel_rezone_all, 19, local_work_size * sizeof(cl_real4_t),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_rezone_all,   1, NULL, &global_work_size, &local_work_size, NULL);

   MallocPlus gpu_state_memory_old = gpu_state_memory;
   malloc_plus_memory_entry *memory_item;

   for (memory_item = gpu_state_memory_old.memory_entry_by_name_begin();
        memory_item != gpu_state_memory_old.memory_entry_by_name_end();
        memory_item = gpu_state_memory_old.memory_entry_by_name_next() ) {
      //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);
      cl_mem dev_state_mem_ptr = (cl_mem)memory_item->mem_ptr;

      if (memory_item->mem_elsize == 8){
#ifndef MINIMUM_PRECISION
         cl_mem dev_state_var_new = (cl_mem)gpu_state_memory.memory_malloc(max(old_ncells,new_ncells), sizeof(cl_double), const_cast<char *>("dev_state_var_new"), DEVICE_REGULAR_MEMORY);

         ezcl_set_kernel_arg(kernel_rezone_one_double, 0, sizeof(cl_int),  (void *)&old_ncells);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 1, sizeof(cl_mem),  (void *)&dev_i);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 2, sizeof(cl_mem),  (void *)&dev_j);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 3, sizeof(cl_mem),  (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 4, sizeof(cl_mem),  (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 5, sizeof(cl_mem),  (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 6, sizeof(cl_mem),  (void *)&dev_ntop);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 7, sizeof(cl_mem),  (void *)&dev_celltype);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 8, sizeof(cl_mem),  (void *)&dev_mpot);
         ezcl_set_kernel_arg(kernel_rezone_one_double, 9, sizeof(cl_mem),  (void *)&dev_indexoffset);
         ezcl_set_kernel_arg(kernel_rezone_one_double,10, sizeof(cl_mem),  (void *)&dev_state_mem_ptr);
         ezcl_set_kernel_arg(kernel_rezone_one_double,11, sizeof(cl_mem),  (void *)&dev_state_var_new);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_rezone_one_double,   1, NULL, &global_work_size, &local_work_size, NULL);

         gpu_state_memory.memory_replace(dev_state_mem_ptr, dev_state_var_new);
#else
         printf("ERROR -- can't have double type for state variable\n");
         exit(1);
#endif
      } else if (memory_item->mem_elsize == 4){
         cl_mem dev_state_var_new = (cl_mem)gpu_state_memory.memory_malloc(max(old_ncells,new_ncells), sizeof(cl_float), const_cast<char *>("dev_state_var_new"), DEVICE_REGULAR_MEMORY);

         ezcl_set_kernel_arg(kernel_rezone_one_float, 0, sizeof(cl_int),  (void *)&old_ncells);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 1, sizeof(cl_mem),  (void *)&dev_i);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 2, sizeof(cl_mem),  (void *)&dev_j);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 3, sizeof(cl_mem),  (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 4, sizeof(cl_mem),  (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 5, sizeof(cl_mem),  (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 6, sizeof(cl_mem),  (void *)&dev_ntop);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 7, sizeof(cl_mem),  (void *)&dev_celltype);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 8, sizeof(cl_mem),  (void *)&dev_mpot);
         ezcl_set_kernel_arg(kernel_rezone_one_float, 9, sizeof(cl_mem),  (void *)&dev_indexoffset);
         ezcl_set_kernel_arg(kernel_rezone_one_float,10, sizeof(cl_mem),  (void *)&dev_state_mem_ptr);
         ezcl_set_kernel_arg(kernel_rezone_one_float,11, sizeof(cl_mem),  (void *)&dev_state_var_new);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_rezone_one_float,   1, NULL, &global_work_size, &local_work_size, NULL);

         gpu_state_memory.memory_replace(dev_state_mem_ptr, dev_state_var_new);
      }
   }

   if (neighbor_remap & ! parallel) {
      size_t mem_request = (int)((float)new_ncells*mem_factor);
      cl_mem dev_nlft_new = ezcl_malloc(NULL, const_cast<char *>("dev_nlft_new"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nrht_new = ezcl_malloc(NULL, const_cast<char *>("dev_nrht_new"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nbot_new = ezcl_malloc(NULL, const_cast<char *>("dev_nbot_new"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_ntop_new = ezcl_malloc(NULL, const_cast<char *>("dev_ntop_new"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_neighbor_init,  0, sizeof(cl_int),   (void *)&new_ncells);
      ezcl_set_kernel_arg(kernel_neighbor_init,  1, sizeof(cl_mem),   (void *)&dev_nlft_new);
      ezcl_set_kernel_arg(kernel_neighbor_init,  2, sizeof(cl_mem),   (void *)&dev_nrht_new);
      ezcl_set_kernel_arg(kernel_neighbor_init,  3, sizeof(cl_mem),   (void *)&dev_nbot_new);
      ezcl_set_kernel_arg(kernel_neighbor_init,  4, sizeof(cl_mem),   (void *)&dev_ntop_new);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_neighbor_init,   1, NULL, &global_work_size, &local_work_size, NULL);

      ezcl_set_kernel_arg(kernel_rezone_neighbors,  0, sizeof(cl_int),  (void *)&old_ncells);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  1, sizeof(cl_mem),  (void *)&dev_mpot);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  2, sizeof(cl_mem),  (void *)&dev_indexoffset);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  3, sizeof(cl_mem),  (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  4, sizeof(cl_mem),  (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  5, sizeof(cl_mem),  (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  6, sizeof(cl_mem),  (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  7, sizeof(cl_mem),  (void *)&dev_celltype_new);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  8, sizeof(cl_mem),  (void *)&dev_nlft_new);
      ezcl_set_kernel_arg(kernel_rezone_neighbors,  9, sizeof(cl_mem),  (void *)&dev_nrht_new);
      ezcl_set_kernel_arg(kernel_rezone_neighbors, 10, sizeof(cl_mem),  (void *)&dev_nbot_new);
      ezcl_set_kernel_arg(kernel_rezone_neighbors, 11, sizeof(cl_mem),  (void *)&dev_ntop_new);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_rezone_neighbors,   1, NULL, &global_work_size, &local_work_size, NULL);
   
      ezcl_device_memory_swap(&dev_nlft, &dev_nlft_new);
      ezcl_device_memory_swap(&dev_nrht, &dev_nrht_new);
      ezcl_device_memory_swap(&dev_nbot, &dev_nbot_new);
      ezcl_device_memory_swap(&dev_ntop, &dev_ntop_new);

      ezcl_device_memory_delete(dev_nlft_new);
      ezcl_device_memory_delete(dev_nrht_new);
      ezcl_device_memory_delete(dev_nbot_new);
      ezcl_device_memory_delete(dev_ntop_new);
   } else {
      ezcl_device_memory_delete(dev_nlft);
      ezcl_device_memory_delete(dev_nrht);
      ezcl_device_memory_delete(dev_nbot);
      ezcl_device_memory_delete(dev_ntop);
      dev_nlft = NULL;
      dev_nrht = NULL;
      dev_nbot = NULL;
      dev_ntop = NULL;
   }

   ezcl_device_memory_delete(dev_indexoffset);

   if (new_ncells != old_ncells){
      resize_old_device_memory(new_ncells);
   }

   ezcl_device_memory_swap(&dev_celltype, &dev_celltype_new);
   ezcl_device_memory_swap(&dev_level, &dev_level_new);
   ezcl_device_memory_swap(&dev_i, &dev_i_new);
   ezcl_device_memory_swap(&dev_j, &dev_j_new);

   ezcl_device_memory_delete(dev_mpot);
   ezcl_device_memory_delete(dev_ijadd);
   ezcl_device_memory_delete(dev_ioffset);

   ezcl_device_memory_delete(dev_i_new);
   ezcl_device_memory_delete(dev_j_new);
   ezcl_device_memory_delete(dev_celltype_new);
   ezcl_device_memory_delete(dev_level_new);

#ifdef HAVE_MPI
   if (parallel) {
      int new_ncells = ncells + add_ncells;
      MPI_Allgather(&new_ncells, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);

      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }
      noffset=ndispl[mype];
      ncells_global = ndispl[numpe-1]+nsizes[numpe-1];
   }
#endif

   gpu_timers[MESH_TIMER_REZONE_ALL] += (long)(cpu_timer_stop(tstart_cpu) * 1.0e9);
}
#endif

void Mesh::calc_neighbors(int ncells)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   if (do_rezone) {

      int flags = INDEX_ARRAY_MEMORY;

#if defined (HAVE_J7)
      if (parallel) flags |= LOAD_BALANCE_MEMORY;
#endif

      static int nlft_size = 0;

#ifdef _OPENMP
#pragma omp master
      {
#endif
      cpu_counters[MESH_COUNTER_CALC_NEIGH]++;

      if (nlft != NULL){
         nlft_size = mesh_memory.get_memory_size(nlft);
      }

      if (nlft_size < ncells) {
            if (nlft != NULL){
               nlft = (int *)mesh_memory.memory_delete(nlft);
               nrht = (int *)mesh_memory.memory_delete(nrht);
               nbot = (int *)mesh_memory.memory_delete(nbot);
               ntop = (int *)mesh_memory.memory_delete(ntop);
            }

            nlft = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "nlft", flags);
            nrht = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "nrht", flags);
            nbot = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "nbot", flags);
            ntop = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "ntop", flags);
      }

#ifdef _OPENMP
      }
#pragma omp barrier
#endif

      if (nlft_size < ncells){
         int lowerBounds, upperBounds;
         get_bounds(lowerBounds, upperBounds);

#ifdef _OPENMP
#pragma omp for
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
         for(int ic=lowerBounds; ic<upperBounds; ic++){
            nlft[ic] = -1;
            nrht[ic] = -1;
            nbot[ic] = -1;
            ntop[ic] = -1;
         }
      }

      if (calc_neighbor_type == HASH_TABLE) {

         struct timespec tstart_lev2;
         if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

         int jmaxsize = (jmax+1)*IPOW2(levmx);
         int imaxsize = (imax+1)*IPOW2(levmx);

         int *hash;

         hash = compact_hash_init(ncells, imaxsize, jmaxsize, METHOD_UNSET, 0);

#ifdef _OPENMP
#pragma omp for
#endif
            for(int ic=0; ic<ncells; ic++){
               int lev = level[ic];
               //printf("%d) %d\n", ic, lev);

               bool need_hash = (nlft[ic] == -1 || nrht[ic] == -1 || nbot[ic] == -1 || ntop[ic] == -1) ? true : false;

               if (! need_hash){
                   if ( ((int)level[nlft[ic]] > lev && ntop[nlft[ic]] == -1) || 
                        ((int)level[nrht[ic]] > lev && ntop[nrht[ic]] == -1) ||
                        ((int)level[nbot[ic]] > lev && nrht[nbot[ic]] == -1) || 
                        ((int)level[ntop[ic]] > lev && nrht[ntop[ic]] == -1) ) need_hash = true;
               }
            
               if (need_hash) {
                  int levmult = IPOW2(levmx-lev);
                  int ii = i[ic]*levmult;
                  int jj = j[ic]*levmult;

                  //if (ic == 284) printf("%d %d %d %d\n", i[ic], j[ic], levmx, lev);
                  write_hash(ic,jj*imaxsize+ii,hash);
               }
            }

//#ifdef _OPENMP
//#pragma omp barrier
//#pragma omp master
//    printf("here\n");
//#endif

            if (TIMING_LEVEL >= 2) {

#ifdef _OPENMP
#pragma omp master
#endif
               cpu_timers[MESH_TIMER_HASH_SETUP] += cpu_timer_stop(tstart_lev2);
               cpu_timer_start(&tstart_lev2);
            }

            //fprintf(fp,"DEBUG ncells is %lu\n",ncells);
#ifdef _OPENMP
#pragma omp for
#endif
            for (int ic=0; ic<(int)ncells; ic++){
               int ii = i[ic];
               int jj = j[ic];
               int lev = level[ic];
               int levmult = IPOW2(levmx-lev);
               int iicur = ii*levmult;
               int iilft = max( (ii-1)*levmult, 0         );
               int iirht = min( (ii+1)*levmult, imaxsize-1);
               int jjcur = jj*levmult;
               int jjbot = max( (jj-1)*levmult, 0         );
               int jjtop = min( (jj+1)*levmult, jmaxsize-1);

               int nlftval = nlft[ic];
               int nrhtval = nrht[ic];
               int nbotval = nbot[ic];
               int ntopval = ntop[ic];

               // Taking care of boundary cells
               // Force each boundary cell to point to itself on its boundary direction
               if (nlftval < 0 && iicur <    1*IPOW2(levmx)  ) nlftval = ic;
               if (nbotval < 0 && jjcur <    1*IPOW2(levmx)  ) nbotval = ic;
               if (nrhtval < 0 && iicur > imax*IPOW2(levmx)-1) nrhtval = ic;
               if (ntopval < 0 && jjcur > jmax*IPOW2(levmx)-1) ntopval = ic;
               // Boundary cells next to corner boundary need special checks
               if (nlftval < 0 && iicur ==    1*IPOW2(levmx) &&  (jjcur < 1*IPOW2(levmx) || jjcur >= jmax*IPOW2(levmx) ) ) nlftval = ic;
               if (nbotval < 0 && jjcur ==    1*IPOW2(levmx) &&  (iicur < 1*IPOW2(levmx) || iicur >= imax*IPOW2(levmx) ) ) nbotval = ic;
               if (nrhtval < 0 && iirht == imax*IPOW2(levmx) &&  (jjcur < 1*IPOW2(levmx) || jjcur >= jmax*IPOW2(levmx) ) ) nrhtval = ic;
               if (ntopval < 0 && jjtop == jmax*IPOW2(levmx) &&  (iicur < 1*IPOW2(levmx) || iicur >= imax*IPOW2(levmx) ) ) ntopval = ic;

               // need to check for finer neighbor first
               // Right and top neighbor don't change for finer, so drop through to same size
               // Left and bottom need to be half of same size index for finer test
               if (lev != levmx) {
                  int iilftfiner = iicur-(iicur-iilft)/2;
                  //int iirhtfiner = (iicur+iirht)/2;
                  int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                  //int jjtopfiner = (jjcur+jjtop)/2;
                  if (nlftval < 0) nlftval = read_hash(jjcur*imaxsize+iilftfiner, hash);
                  if (nbotval < 0) nbotval = read_hash(jjbotfiner*imaxsize+iicur, hash);
               }

               // same size neighbor
               if (nlftval < 0) nlftval = read_hash(jjcur*imaxsize+iilft, hash);
               if (nrhtval < 0) nrhtval = read_hash(jjcur*imaxsize+iirht, hash);
               if (nbotval < 0) nbotval = read_hash(jjbot*imaxsize+iicur, hash);
               if (ntopval < 0) ntopval = read_hash(jjtop*imaxsize+iicur, hash);

               // Now we need to take care of special case where bottom and left boundary need adjustment since
               // expected cell doesn't exist on these boundaries if it is finer than current cell
               if (lev != levmx) {
                  if (jjcur < 1*IPOW2(levmx)) {
                     if (nrhtval < 0) {
                        int jjtopfiner = (jjcur+jjtop)/2;
                        nrhtval = read_hash(jjtopfiner*imaxsize+iirht, hash);
                     }
                     if (nlftval < 0) {
                        int iilftfiner = iicur-(iicur-iilft)/2;
                        int jjtopfiner = (jjcur+jjtop)/2;
                        nlftval = read_hash(jjtopfiner*imaxsize+iilftfiner, hash);
                     }
                  }
            
                  if (iicur < 1*IPOW2(levmx)) {
                     if (ntopval < 0) {
                        int iirhtfiner = (iicur+iirht)/2;
                        ntopval = read_hash(jjtop*imaxsize+iirhtfiner, hash);
                     }
                     if (nbotval < 0) {
                        int iirhtfiner = (iicur+iirht)/2;
                        int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                        nbotval = read_hash(jjbotfiner*imaxsize+iirhtfiner, hash);
                     }
                  }
               }
            
               // coarser neighbor
               if (lev != 0){
                  if (nlftval < 0) {
                     iilft -= iicur-iilft;
                     int jjlft = (jj/2)*2*levmult;
                     nlftval = read_hash(jjlft*imaxsize+iilft, hash);
                  }
                  if (nrhtval < 0) {
                     int jjrht = (jj/2)*2*levmult;
                     nrhtval = read_hash(jjrht*imaxsize+iirht, hash);
                  }
                  if (nbotval < 0) {
                     jjbot -= jjcur-jjbot;
                     int iibot = (ii/2)*2*levmult;
                     nbotval = read_hash(jjbot*imaxsize+iibot, hash);
                  }
                  if (ntopval < 0) {
                     int iitop = (ii/2)*2*levmult;
                     ntopval = read_hash(jjtop*imaxsize+iitop, hash);
                  }
               }

               nlft[ic] = nlftval;
               nrht[ic] = nrhtval;
               nbot[ic] = nbotval;
               ntop[ic] = ntopval;

               //printf("neighbors[%d] = %d %d %d %d\n",ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         write_hash_collision_report();
         read_hash_collision_report();

         compact_hash_delete(hash);

         if (TIMING_LEVEL >= 2) cpu_timers[MESH_TIMER_HASH_QUERY] += cpu_timer_stop(tstart_lev2);
#ifdef _OPENMP
         } // master block
#endif

      } else if (calc_neighbor_type == KDTREE) {

         struct timespec tstart_lev2;
         if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         TBounds box;
         vector<int> index_list(IPOW2(levmx*levmx) );

         int num;

         ibase = 0;
         calc_spatial_coordinates(ibase);

         kdtree_setup();

         if (TIMING_LEVEL >= 2) {
            cpu_timers[MESH_TIMER_KDTREE_SETUP] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         for (int ic=0; ic<ncells; ic++) {

            //left
            nlft[ic]  = ic;
            box.min.x = x[ic]-0.25*dx[ic];
            box.max.x = x[ic]-0.25*dx[ic];
            box.min.y = y[ic]+0.25*dy[ic];
            box.max.y = y[ic]+0.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) nlft[ic]=index_list[0];

            //right
            nrht[ic]  = ic;
            box.min.x = x[ic]+1.25*dx[ic];
            box.max.x = x[ic]+1.25*dx[ic];
            box.min.y = y[ic]+0.25*dy[ic];
            box.max.y = y[ic]+0.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) nrht[ic]=index_list[0];

            //bot
            nbot[ic]  = ic;
            box.min.x = x[ic]+0.25*dx[ic];
            box.max.x = x[ic]+0.25*dx[ic];
            box.min.y = y[ic]-0.25*dy[ic];
            box.max.y = y[ic]-0.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) nbot[ic]=index_list[0];

            //top
            ntop[ic]  = ic;
            box.min.x = x[ic]+0.25*dx[ic];
            box.max.x = x[ic]+0.25*dx[ic];
            box.min.y = y[ic]+1.25*dy[ic];
            box.max.y = y[ic]+1.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) ntop[ic]=index_list[0];
         }  //  End main loop over cells.

         KDTree_Destroy(&tree);

         if (TIMING_LEVEL >= 2) cpu_timers[MESH_TIMER_KDTREE_QUERY] += cpu_timer_stop(tstart_lev2);

#ifdef _OPENMP
         }
#pragma omp barrier
#endif
      } // calc_neighbor_type

#ifdef _OPENMP
#pragma omp master
#endif
      ncells_ghost = ncells;

   }

#ifdef _OPENMP
#pragma omp master
#endif
   cpu_timers[MESH_TIMER_CALC_NEIGHBORS] += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_neighbors_local(void)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   if (do_rezone) {

      int flags = INDEX_ARRAY_MEMORY;

#if defined (HAVE_J7)
      if (parallel) flags |= LOAD_BALANCE_MEMORY;
#endif

#ifdef _OPENMP
#pragma omp master
      {
#endif
      cpu_counters[MESH_COUNTER_CALC_NEIGH]++;

      if (mesh_memory.get_memory_size(nlft) < ncells){
         if (nlft != NULL) nlft = (int *)mesh_memory.memory_delete(nlft);
         if (nrht != NULL) nrht = (int *)mesh_memory.memory_delete(nrht);
         if (nbot != NULL) nbot = (int *)mesh_memory.memory_delete(nbot);
         if (ntop != NULL) ntop = (int *)mesh_memory.memory_delete(ntop);
         nlft = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "nlft", flags);
         nrht = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "nrht", flags);
         nbot = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "nbot", flags);
         ntop = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "ntop", flags);
      }
#ifdef _OPENMP
      }
#pragma omp barrier
#endif

      int lowerBound, upperBound;
      set_bounds(ncells);
      get_bounds(lowerBound, upperBound);
      for (int ic = lowerBound; ic < upperBound; ic++){
         nlft[ic] = -98;
         nrht[ic] = -98;
         nbot[ic] = -98;
         ntop[ic] = -98;
      }

      if (calc_neighbor_type == HASH_TABLE) {

         struct timespec tstart_lev2;
         if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

         ncells_ghost = ncells;

         // Find maximum i column and j row for this processor
         static int jmintile, imintile, jmaxtile, imaxtile;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         jmintile = (jmax+1)*IPOW2(levmx);
         imintile = (imax+1)*IPOW2(levmx);
         jmaxtile = 0;
         imaxtile = 0;
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         int my_jmintile = jmintile;
         int my_imintile = imintile;
         int my_jmaxtile = 0;
         int my_imaxtile = 0;

#ifdef _OPENMP
#pragma omp for
#endif
         for(uint ic=0; ic<ncells; ic++){
            int lev = level[ic];
//          if (lev < 0 || lev > levmx) printf("DEBUG -- cell %d lev %d\n",ic,level[ic]);
            if ( j[ic]   *IPOW2(levmx-lev)   < my_jmintile) my_jmintile =  j[ic]   *IPOW2(levmx-lev)  ;
            if ((j[ic]+1)*IPOW2(levmx-lev)-1 > my_jmaxtile) my_jmaxtile = (j[ic]+1)*IPOW2(levmx-lev)-1;
            if ( i[ic]   *IPOW2(levmx-lev)   < my_imintile) my_imintile =  i[ic]   *IPOW2(levmx-lev)  ;
            if ((i[ic]+1)*IPOW2(levmx-lev)-1 > my_imaxtile) my_imaxtile = (i[ic]+1)*IPOW2(levmx-lev)-1;
         }
#ifdef _OPENMP
#pragma omp critical
         {
#endif
            if (my_jmintile < jmintile) jmintile = my_jmintile;
            if (my_imintile < imintile) imintile = my_imintile;
            if (my_jmaxtile > jmaxtile) jmaxtile = my_jmaxtile;
            if (my_imaxtile > imaxtile) imaxtile = my_imaxtile;
#ifdef _OPENMP
         } // end critical region
#pragma omp barrier
#endif

      //if (DEBUG) fprintf(fp,"%d: Tile Sizes are imin %d imax %d jmin %d jmax %d\n",mype,imintile,imaxtile,jmintile,jmaxtile);

      // Expand size by 2*coarse_cells for ghost cells
      int jminsize = max(jmintile-2*IPOW2(levmx),0);
      int jmaxsize = min(jmaxtile+2*IPOW2(levmx),(jmax+1)*IPOW2(levmx));
      int iminsize = max(imintile-2*IPOW2(levmx),0);
      int imaxsize = min(imaxtile+2*IPOW2(levmx),(imax+1)*IPOW2(levmx));
      //if (DEBUG) fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);

      //fprintf(fp,"DEBUG -- ncells %lu\n",ncells);

      static int *hash;

      hash = compact_hash_init(ncells, imaxsize-iminsize, jmaxsize-jminsize, METHOD_UNSET, 0);

      //printf("%d: DEBUG -- noffset %d cells %d\n",mype,noffset,ncells);

      if (DEBUG) {
#ifdef _OPENMP
#pragma omp master
#endif
         fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);
      }

      static int imaxcalc, jmaxcalc;

#ifdef _OPENMP
#pragma omp for
#endif
      for(uint ic=0; ic<ncells; ic++){
         int cellnumber = ic+noffset;
         int lev = level[ic];
         int levmult = IPOW2(levmx-lev);
         int ii = i[ic]*levmult-iminsize;
         int jj = j[ic]*levmult-jminsize;

         write_hash(cellnumber, jj*(imaxsize-iminsize)+ii, hash);
      } // end for loop

      if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
         cpu_timers[MESH_TIMER_HASH_SETUP] += cpu_timer_stop(tstart_lev2);
         cpu_timer_start(&tstart_lev2);
      }

#ifdef _OPENMP
#pragma omp master
      {
#endif
      // Set neighbors to global cell numbers from hash
      jmaxcalc = (jmax+1)*IPOW2(levmx);
      imaxcalc = (imax+1)*IPOW2(levmx);
#ifdef _OPENMP
      }
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
      for (uint ic=0; ic<ncells; ic++){
         int ii = i[ic];
         int jj = j[ic];
         int lev = level[ic];
         int levmult = IPOW2(levmx-lev);

         int iicur = ii*levmult-iminsize;
         int iilft = max( (ii-1)*levmult, 0         )-iminsize;
         int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;   
         int jjcur = jj*levmult-jminsize;
         int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
         int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;   

         int nlftval = -1;
         int nrhtval = -1;
         int nbotval = -1;
         int ntopval = -1;

         // Taking care of boundary cells
         // Force each boundary cell to point to itself on its boundary direction
         if (iicur <    1*IPOW2(levmx)  -iminsize) nlftval = ic+noffset;
         if (jjcur <    1*IPOW2(levmx)  -jminsize) nbotval = ic+noffset;
         if (iicur > imax*IPOW2(levmx)-1-iminsize) nrhtval = ic+noffset;
         if (jjcur > jmax*IPOW2(levmx)-1-jminsize) ntopval = ic+noffset;
         // Boundary cells next to corner boundary need special checks
         if (iicur ==    1*IPOW2(levmx)-iminsize &&  (jjcur < 1*IPOW2(levmx)-jminsize || jjcur >= jmax*IPOW2(levmx)-jminsize ) ) nlftval = ic+noffset;
         if (jjcur ==    1*IPOW2(levmx)-jminsize &&  (iicur < 1*IPOW2(levmx)-iminsize || iicur >= imax*IPOW2(levmx)-iminsize ) ) nbotval = ic+noffset;
         if (iirht == imax*IPOW2(levmx)-iminsize &&  (jjcur < 1*IPOW2(levmx)-jminsize || jjcur >= jmax*IPOW2(levmx)-jminsize ) ) nrhtval = ic+noffset;
         if (jjtop == jmax*IPOW2(levmx)-jminsize &&  (iicur < 1*IPOW2(levmx)-iminsize || iicur >= imax*IPOW2(levmx)-iminsize ) ) ntopval = ic+noffset;

         // need to check for finer neighbor first
         // Right and top neighbor don't change for finer, so drop through to same size
         // Left and bottom need to be half of same size index for finer test
         if (lev != levmx) {
            int iilftfiner = iicur-(iicur-iilft)/2;
            int jjbotfiner = jjcur-(jjcur-jjbot)/2;
            if (nlftval < 0) nlftval = read_hash(jjcur     *(imaxsize-iminsize)+iilftfiner, hash);
            if (nbotval < 0) nbotval = read_hash(jjbotfiner*(imaxsize-iminsize)+iicur,      hash);
         }

         // same size neighbor
         if (nlftval < 0) {
            int nlfttry = read_hash(jjcur*(imaxsize-iminsize)+iilft, hash);
            if (nlfttry >= 0 && nlfttry < (int)ncells && (int)level[nlfttry] == lev) nlftval = nlfttry;
         }
         if (nrhtval < 0) nrhtval = read_hash(jjcur*(imaxsize-iminsize)+iirht, hash);
         if (nbotval < 0) {
            int nbottry = read_hash(jjbot*(imaxsize-iminsize)+iicur, hash);
            if (nbottry >= 0 && nbottry < (int)ncells && (int)level[nbottry] == lev) nbotval = nbottry;
         }
         if (ntopval < 0) ntopval = read_hash(jjtop*(imaxsize-iminsize)+iicur, hash);
              
         // Now we need to take care of special case where bottom and left boundary need adjustment since
         // expected cell doesn't exist on these boundaries if it is finer than current cell
         if (lev != levmx) {
            if (jjcur < 1*IPOW2(levmx)) {
               if (nrhtval < 0) {
                  int jjtopfiner = (jjcur+jjtop)/2;
                  nrhtval = read_hash(jjtopfiner*(imaxsize-iminsize)+iirht, hash);
               }
               if (nlftval < 0) {
                  int iilftfiner = iicur-(iicur-iilft)/2;
                  int jjtopfiner = (jjcur+jjtop)/2;
                  nlftval = read_hash(jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
               }
            }

            if (iicur < 1*IPOW2(levmx)) {
               if (ntopval < 0) {
                  int iirhtfiner = (iicur+iirht)/2;
                  ntopval = read_hash(jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
               }
               if (nbotval < 0) {
                  int iirhtfiner = (iicur+iirht)/2;
                  int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                  nbotval = read_hash(jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
               }
            }
         }

         // coarser neighbor
         if (lev != 0){
            if (nlftval < 0) {
               iilft -= iicur-iilft;
               int jjlft = (jj/2)*2*levmult-jminsize;
               int nlfttry = read_hash(jjlft*(imaxsize-iminsize)+iilft, hash);
               if (nlfttry >= 0 && nlfttry < (int)ncells && (int)level[nlfttry] == lev-1) nlftval = nlfttry;
            }       
            if (nrhtval < 0) {
               int jjrht = (jj/2)*2*levmult-jminsize;
               int nrhttry = read_hash(jjrht*(imaxsize-iminsize)+iirht, hash);
               if (nrhttry >= 0 && nrhttry < (int)ncells && (int)level[nrhttry] == lev-1) nrhtval = nrhttry;
            }       
            if (nbotval < 0) {
               jjbot -= jjcur-jjbot;
               int iibot = (ii/2)*2*levmult-iminsize;
               int nbottry = read_hash(jjbot*(imaxsize-iminsize)+iibot, hash);
               if (nbottry >= 0 && nbottry < (int)ncells && (int)level[nbottry] == lev-1) nbotval = nbottry;
            }       
            if (ntopval < 0) {
               int iitop = (ii/2)*2*levmult-iminsize;
               int ntoptry = read_hash(jjtop*(imaxsize-iminsize)+iitop, hash);
               if (ntoptry >= 0 && ntoptry < (int)ncells && (int)level[ntoptry] == lev-1) ntopval = ntoptry;
            }       
         }       

         nlft[ic] = nlftval;
         nrht[ic] = nrhtval;
         nbot[ic] = nbotval;
         ntop[ic] = ntopval;

         //fprintf(fp,"%d: neighbors[%d] = %d %d %d %d\n",mype,ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
      }

      if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         print_local();

         int jmaxglobal = (jmax+1)*IPOW2(levmx);
         int imaxglobal = (imax+1)*IPOW2(levmx);
         fprintf(fp,"\n                                    HASH 0 numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     fprintf(fp,"%5d",read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash));
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");

         fprintf(fp,"\n                                    nlft numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash)-noffset;
                     if (hashval >= 0 && hashval < (int)ncells) {
                        fprintf(fp,"%5d",nlft[hashval]);
                     } else {
                        fprintf(fp,"     ");
                     }
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
      
         fprintf(fp,"\n                                    nrht numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash)-noffset;
                     if (hashval >= 0 && hashval < (int)ncells) {
                        fprintf(fp,"%5d",nrht[hashval]);
                     } else {
                        fprintf(fp,"     ");
                     }
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");

         fprintf(fp,"\n                                    nbot numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash)-noffset;
                     if (hashval >= 0 && hashval < (int)ncells) {
                        fprintf(fp,"%5d",nbot[hashval]);
                     } else {
                        fprintf(fp,"     ");
                     }
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");

         fprintf(fp,"\n                                    ntop numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash)-noffset;
                     if (hashval >= 0 && hashval < (int)ncells) {
                        fprintf(fp,"%5d",ntop[hashval]);
                     } else {
                        fprintf(fp,"     ");
                     }
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
#ifdef _OPENMP
         }
#pragma omp barrier
#endif
      }

      if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
         cpu_timers[MESH_TIMER_HASH_QUERY] += cpu_timer_stop(tstart_lev2);
         cpu_timer_start(&tstart_lev2);
      }

#ifdef HAVE_MPI
      if (numpe > 1) {
         static int num_comm_partners;

         static vector<int> iminsize_global;
         static vector<int> imaxsize_global;
         static vector<int> jminsize_global;
         static vector<int> jmaxsize_global;
         static vector<int> comm_partner;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         iminsize_global.resize(numpe);
         imaxsize_global.resize(numpe);
         jminsize_global.resize(numpe);
         jmaxsize_global.resize(numpe);
         comm_partner.resize(numpe,-1);

         MPI_Allgather(&iminsize, 1, MPI_INT, &iminsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);
         MPI_Allgather(&imaxsize, 1, MPI_INT, &imaxsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);
         MPI_Allgather(&jminsize, 1, MPI_INT, &jminsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);
         MPI_Allgather(&jmaxsize, 1, MPI_INT, &jmaxsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);

         num_comm_partners = 0;
         for (int ip = 0; ip < numpe; ip++){
            if (ip == mype) continue;
            if (iminsize_global[ip] > imaxtile) continue;
            if (imaxsize_global[ip] < imintile) continue;
            if (jminsize_global[ip] > jmaxtile) continue;
            if (jmaxsize_global[ip] < jmintile) continue;
            comm_partner[num_comm_partners] = ip;
            num_comm_partners++;
            //if (DEBUG) fprintf(fp,"%d: overlap with processor %d bounding box is %d %d %d %d\n",mype,ip,iminsize_global[ip],imaxsize_global[ip],jminsize_global[ip],jmaxsize_global[ip]);
         }
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         static vector<int> border_cell;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         border_cell.resize(ncells);

#ifdef BOUNDS_CHECK
         for (uint ic=0; ic<ncells; ic++){
            int nl = nlft[ic];
            if (nl != -1){
               nl -= noffset;
               if (nl<0 || nl>= (int)ncells) printf("%d: Warning at line %d cell %d nlft %d\n",mype,__LINE__,ic,nl);
            }
            int nr = nrht[ic];
            if (nr != -1){
               nr -= noffset;
               if (nr<0 || nr>= (int)ncells) printf("%d: Warning at line %d cell %d nrht %d\n",mype,__LINE__,ic,nr);
            }
            int nb = nbot[ic];
            if (nb != -1){
               nb -= noffset;
               if (nb<0 || nb>= (int)ncells) printf("%d: Warning at line %d cell %d nbot %d\n",mype,__LINE__,ic,nb);
            }
            int nt = ntop[ic];
            if (nt != -1){
               nt -= noffset;
               if (nt<0 || nt>= (int)ncells) printf("%d: Warning at line %d cell %d ntop %d\n",mype,__LINE__,ic,nt);
            }
         }
#endif

#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         static vector<int> border_cell_out;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         border_cell_out.resize(ncells);
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
            for (uint ic=0; ic<ncells; ic++){
               int iborder_cell = 0;

               // left neighbor is undefined -- or -- if left is at finer level check left top for undefined
               if (nlft[ic] == -1 || (level[nlft[ic]-noffset] > level[ic] && ntop[nlft[ic]-noffset] == -1) ){
                  iborder_cell |= 0x0001;
               }
               if (nrht[ic] == -1 || (level[nrht[ic]-noffset] > level[ic] && ntop[nrht[ic]-noffset] == -1) ){
                  iborder_cell |= 0x0002;
               }
               if (nbot[ic] == -1 || (level[nbot[ic]-noffset] > level[ic] && nrht[nbot[ic]-noffset] == -1) ) {
                  iborder_cell |= 0x0004;
               }
               if (ntop[ic] == -1 || (level[ntop[ic]-noffset] > level[ic] && nrht[ntop[ic]-noffset] == -1) ) {
                  iborder_cell |= 0x0008;
               }
   
               border_cell[ic] = iborder_cell;
            }

#ifdef _OPENMP
#pragma omp for
#endif
            for (uint ic=0; ic<ncells; ic++){
               int iborder_cell = border_cell[ic];

               if (iborder_cell == 0) {

                  int nl = nlft[ic]-noffset;
                  if (nl >= 0 && nl < (int)ncells) {
                     if ((border_cell[nl] & 0x0001) == 0x0001) {
                        iborder_cell |= 0x0016;
                     } else if (level[nl] > level[ic]){
                        int ntl = ntop[nl]-noffset;
                        if (ntl >= 0 && ntl < (int)ncells && (border_cell[ntl] & 0x0001) == 0x0001) {
                           iborder_cell |= 0x0016;
                        }
                     }
                  }
                  int nr = nrht[ic]-noffset;
                  if (nr >= 0 && nr < (int)ncells) {
                     if ((border_cell[nrht[ic]-noffset] & 0x0002) == 0x0002) {
                        iborder_cell |= 0x0032;
                     } else if (level[nr] > level[ic]){
                        int ntr = ntop[nr]-noffset;
                        if (ntr >= 0 && ntr < (int)ncells && (border_cell[ntr] & 0x0002) == 0x0002) {
                           iborder_cell |= 0x0032;
                        }
                     }
                  }
                  int nb = nbot[ic]-noffset;
                  if (nb >= 0 && nb < (int)ncells) {
                     if ((border_cell[nb] & 0x0004) == 0x0004) {
                        iborder_cell |= 0x0064;
                     } else if (level[nb] > level[ic]){
                        int nrb = nrht[nb]-noffset;
                        if (nrb >= 0 && nrb < (int)ncells && (border_cell[nrb] & 0x0004) == 0x0004) {
                           iborder_cell |= 0x0064;
                        }
                     }
                  }
                  int nt = ntop[ic]-noffset;
                  if (nt >= 0 && nt < (int)ncells) {
                     if ((border_cell[nt] & 0x0008) == 0x0008) {
                        iborder_cell |= 0x0128;
                     } else if (level[nt] > level[ic]){
                        int nrt = nrht[nt]-noffset;
                        if (nrt >= 0 && nrt < (int)ncells && (border_cell[nrt] & 0x0008) == 0x0008) {
                           iborder_cell |= 0x0128;
                        }
                     }
                  }
               }

               border_cell_out[ic] = iborder_cell;
            }
// indent offset

         vector<int> border_cell_num;

         static int nbsize_local;

         static vector<int> border_cell_i;
         static vector<int> border_cell_j;
         static vector<int> border_cell_level;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         for (int ic=0; ic<(int)ncells; ic++){
            if (border_cell_out[ic] > 0) border_cell_num.push_back(ic+noffset);
         }
         //printf("%d: border cell size is %d\n",mype,border_cell_num.size());

         nbsize_local = border_cell_num.size();

         border_cell_i.resize(nbsize_local);
         border_cell_j.resize(nbsize_local);
         border_cell_level.resize(nbsize_local);

         for (int ic = 0; ic <nbsize_local; ic++){
            int cell_num = border_cell_num[ic]-noffset;
            border_cell_i[ic] = i[cell_num]; 
            border_cell_j[ic] = j[cell_num]; 
            border_cell_level[ic] = level[cell_num]; 
         }
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            fprintf(fp,"%d: Border cell size is %d\n",mype,nbsize_local);
            for (int ib = 0; ib <nbsize_local; ib++){
               fprintf(fp,"%d: Border cell %d is %d i %d j %d level %d\n",mype,ib,border_cell_num[ib],
                  border_cell_i[ib],border_cell_j[ib],border_cell_level[ib]);
            }
#ifdef _OPENMP
         }
#pragma omp barrier
#endif
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_FIND_BOUNDARY] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         // Allocate push database

         static int **send_database;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         send_database = (int**)malloc(num_comm_partners*sizeof(int *));
         for (int ip = 0; ip < num_comm_partners; ip++){
            send_database[ip] = (int *)malloc(nbsize_local*sizeof(int));
         }
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         // Compute the overlap between processor bounding boxes and set up push database

         static vector<int> send_buffer_count;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         send_buffer_count.resize(num_comm_partners);
         for (int ip = 0; ip < num_comm_partners; ip++){
            int icount = 0;
            for (int ib = 0; ib <nbsize_local; ib++){
               int lev = border_cell_level[ib];
               int levmult = IPOW2(levmx-lev);
               if (border_cell_i[ib]*levmult >= iminsize_global[comm_partner[ip]] && 
                   border_cell_i[ib]*levmult <= imaxsize_global[comm_partner[ip]] && 
                   border_cell_j[ib]*levmult >= jminsize_global[comm_partner[ip]] && 
                   border_cell_j[ib]*levmult <= jmaxsize_global[comm_partner[ip]] ) {
                  //   border_cell_i[ib],border_cell_j[ib],border_cell_level[ib]);
                  send_database[ip][icount] = ib;
                  icount++;
               }
            }
            send_buffer_count[ip]=icount;
         }
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         // Initialize L7_Push_Setup with num_comm_partners, comm_partner, send_database and 
         // send_buffer_count. L7_Push_Setup will copy data and determine recv_buffer_counts.
         // It will return receive_count_total for use in allocations

         static int receive_count_total;
         int i_push_handle = 0;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         i_push_handle = 0;
         L7_Push_Setup(num_comm_partners, &comm_partner[0], &send_buffer_count[0],
                       send_database, &receive_count_total, &i_push_handle);
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            fprintf(fp,"DEBUG num_comm_partners %d\n",num_comm_partners);
            for (int ip = 0; ip < num_comm_partners; ip++){
               fprintf(fp,"DEBUG comm partner is %d data count is %d\n",comm_partner[ip],send_buffer_count[ip]);
               for (int ic = 0; ic < send_buffer_count[ip]; ic++){
                  int ib = send_database[ip][ic];
                  fprintf(fp,"DEBUG \t index %d cell number %d i %d j %d level %d\n",ib,border_cell_num[ib],
                     border_cell_i[ib],border_cell_j[ib],border_cell_level[ib]);
               }
            }
#ifdef _OPENMP
         }
#endif
         }

         // Can now free the send database. Other arrays are vectors and will automatically 
         // deallocate

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         for (int ip = 0; ip < num_comm_partners; ip++){
            free(send_database[ip]);
         }
         free(send_database);
#ifdef _OPENMP
         }
#endif

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_PUSH_SETUP] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         // Push the data needed to the adjacent processors
         static int *border_cell_num_local;
         static int *border_cell_i_local;
         static int *border_cell_j_local;
         static int *border_cell_level_local;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         border_cell_num_local = (int *)malloc(receive_count_total*sizeof(int));
         border_cell_i_local = (int *)malloc(receive_count_total*sizeof(int));
         border_cell_j_local = (int *)malloc(receive_count_total*sizeof(int));
         border_cell_level_local = (int *)malloc(receive_count_total*sizeof(int));

         L7_Push_Update(&border_cell_num[0],   border_cell_num_local,   i_push_handle);
         L7_Push_Update(&border_cell_i[0],     border_cell_i_local,     i_push_handle);
         L7_Push_Update(&border_cell_j[0],     border_cell_j_local,     i_push_handle);
         L7_Push_Update(&border_cell_level[0], border_cell_level_local, i_push_handle);

         L7_Push_Free(&i_push_handle);
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         nbsize_local = receive_count_total; 

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            for (int ic = 0; ic < nbsize_local; ic++) {
               fprintf(fp,"%d: Local Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num_local[ic],
                  border_cell_i_local[ic],border_cell_j_local[ic],border_cell_level_local[ic]);
            }
#ifdef _OPENMP
         }
#endif
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_PUSH_BOUNDARY] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_LOCAL_LIST] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            int jmaxglobal = (jmax+1)*IPOW2(levmx);
            int imaxglobal = (imax+1)*IPOW2(levmx);
            fprintf(fp,"\n                                    HASH numbering before layer 1\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        fprintf(fp,"%5d",read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash));
                     } else {
                        fprintf(fp,"     ");
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");
#ifdef _OPENMP
         }
#endif
         }

         vector<int> border_cell_needed_local;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         border_cell_needed_local.resize(nbsize_local, 0);
#ifdef _OPENMP
         }
#pragma omp barrier
#endif

         // Layer 1
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         for (int ic =0; ic<nbsize_local; ic++){
            int jj = border_cell_j_local[ic];
            int ii = border_cell_i_local[ic];
            int lev = border_cell_level_local[ic];
            int levmult = IPOW2(levmx-lev);

            int iicur = ii*levmult-iminsize;
            int iilft = max( (ii-1)*levmult, 0         )-iminsize;
            int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
            int jjcur = jj*levmult-jminsize;
            int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
            int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

            //fprintf(fp,"DEBUG layer ic %d num %d i %d j %d lev %d\n",ic,border_cell_num_local[ic],ii,jj,lev);

            int iborder = 0;

            // Test for cell to left
            if (iicur-(iicur-iilft)/2 >= 0 && iicur-(iicur-iilft)/2 < imaxsize-iminsize && jjcur >= 0 && (jjcur+jjtop)/2 < jmaxsize-jminsize){
               int nlftval = -1;
               // Check for finer cell left and bottom side
               if (lev != levmx){                                // finer neighbor
                  int iilftfiner = iicur-(iicur-iilft)/2;
                  nlftval = read_hash(jjcur*(imaxsize-iminsize)+iilftfiner, hash);
                  // Also check for finer cell left and top side
                  if (nlftval < 0) {
                     int jjtopfiner = (jjcur+jjtop)/2; 
                     nlftval = read_hash(jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
                  }
               }

               if (nlftval < 0 && iilft >= 0) {  // same size
                  int nlfttry = read_hash(jjcur*(imaxsize-iminsize)+iilft, hash);
                  // we have to test for same level or it could be a finer cell one cell away that it is matching
                  if (nlfttry-noffset >= 0 && nlfttry-noffset < (int)ncells && (int)level[nlfttry-noffset] == lev) {
                     nlftval = nlfttry;
                  }
               }
    
               if (lev != 0 && nlftval < 0 && iilft-(iicur-iilft) >= 0){      // coarser neighbor
                  iilft -= iicur-iilft;
                  int jjlft = (jj/2)*2*levmult-jminsize;
                  int nlfttry = read_hash(jjlft*(imaxsize-iminsize)+iilft, hash);
                  // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
                  if (nlfttry-noffset >= 0 && nlfttry-noffset < (int)ncells && (int)level[nlfttry-noffset] == lev-1) {
                    nlftval = nlfttry;
                  }
               }
               if (nlftval >= 0) iborder |= 0x0001;
            }

            // Test for cell to right
            if (iirht < imaxsize-iminsize && iirht >= 0 && jjcur >= 0 && jjtop < jmaxsize-jminsize) {
               int nrhtval = -1;
               // right neighbor -- finer, same size and coarser
               nrhtval = read_hash(jjcur*(imaxsize-iminsize)+iirht, hash);
               // right neighbor -- finer right top test
               if (nrhtval < 0 && lev != levmx){
                  int jjtopfiner = (jjcur+jjtop)/2;
                  nrhtval = read_hash(jjtopfiner*(imaxsize-iminsize)+iirht, hash);
               }
               if (nrhtval < 0 && lev != 0) { // test for coarser, but not directly above
                  int jjrhtcoarser = (jj/2)*2*levmult-jminsize;
                  if (jjrhtcoarser != jjcur) {
                     int nrhttry = read_hash(jjrhtcoarser*(imaxsize-iminsize)+iirht, hash);
                     if (nrhttry-noffset >= 0 && nrhttry-noffset < (int)ncells && (int)level[nrhttry-noffset] == lev-1) {
                        nrhtval = nrhttry;
                     }
                  }
               }
               if (nrhtval > 0)  iborder |= 0x0002;
            }

            // Test for cell to bottom
            if (iicur >= 0 && (iicur+iirht)/2 < imaxsize-iminsize && jjcur-(jjcur-jjbot)/2 >= 0 && jjcur-(jjcur-jjbot)/2 < jmaxsize-jminsize){
               int nbotval = -1;
               // Check for finer cell below and left side
               if (lev != levmx){                                // finer neighbor
                  int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                  nbotval = read_hash(jjbotfiner*(imaxsize-iminsize)+iicur, hash);
                  // Also check for finer cell below and right side
                  if (nbotval < 0) {
                     int iirhtfiner = (iicur+iirht)/2; 
                     nbotval = read_hash(jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
                  }
               }

               if (nbotval < 0 && jjbot >= 0) {  // same size
                  int nbottry = read_hash(jjbot*(imaxsize-iminsize)+iicur, hash);
                  // we have to test for same level or it could be a finer cell one cell away that it is matching
                  if (nbottry-noffset >= 0 && nbottry-noffset < (int)ncells && (int)level[nbottry-noffset] == lev) {
                     nbotval = nbottry;
                  }
               }
    
               if (lev != 0 && nbotval < 0 && jjbot-(jjcur-jjbot) >= 0){      // coarser neighbor
                  jjbot -= jjcur-jjbot;
                  int iibot = (ii/2)*2*levmult-iminsize;
                  int nbottry = read_hash(jjbot*(imaxsize-iminsize)+iibot, hash);
                  // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
                  if (nbottry-noffset >= 0 && nbottry-noffset < (int)ncells && (int)level[nbottry-noffset] == lev-1) {
                    nbotval = nbottry;
                  }
               }
               if (nbotval >= 0) iborder |= 0x0004;
            }

            // Test for cell to top
            if (iirht < imaxsize-iminsize && iicur >= 0 && jjtop >= 0 && jjtop < jmaxsize-jminsize) {
               int ntopval = -1;
               // top neighbor -- finer, same size and coarser
               ntopval = read_hash(jjtop*(imaxsize-iminsize)+iicur, hash);
               // top neighbor -- finer top right test
               if (ntopval < 0 && lev != levmx){
                  int iirhtfiner = (iicur+iirht)/2;
                  ntopval = read_hash(jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
               }
               if (ntopval < 0 && lev != 0) { // test for coarser, but not directly above
                  int iitopcoarser = (ii/2)*2*levmult-iminsize;
                  if (iitopcoarser != iicur) {
                     int ntoptry = read_hash(jjtop*(imaxsize-iminsize)+iitopcoarser, hash);
                     if (ntoptry-noffset >= 0 && ntoptry-noffset < (int)ncells && (int)level[ntoptry-noffset] == lev-1) {
                        ntopval = ntoptry;
                     }
                  }
               }
               if (ntopval > 0)  iborder |= 0x0008;
            }

            if (iborder) border_cell_needed_local[ic] = iborder;
         }
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            for(int ic=0; ic<nbsize_local; ic++){
               if (border_cell_needed_local[ic] == 0) continue;
               fprintf(fp,"%d: First set of needed cells ic %3d cell %3d type %3d\n",mype,ic,border_cell_num_local[ic],border_cell_needed_local[ic]);
            }
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif
         }

         // Walk through cell array and set hash to border local index plus ncells+noffset for next pass
         //fprintf(fp,"%d: DEBUG new hash jminsize %d jmaxsize %d iminsize %d imaxsize %d\n",mype,jminsize,jmaxsize,iminsize,imaxsize);
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         for(int ic=0; ic<nbsize_local; ic++){
            if (border_cell_needed_local[ic] == 0) continue;
            //fprintf(fp,"%d: index %d cell %d i %d j %d\n",mype,ic,border_cell_num_local[ic],border_cell_i_local[ic],border_cell_j_local[ic]);
            int lev = border_cell_level_local[ic];
            int levmult = IPOW2(levmx-lev);
            int ii = border_cell_i_local[ic]*levmult-iminsize;
            int jj = border_cell_j_local[ic]*levmult-jminsize;

            write_hash(ncells+noffset+ic, jj*(imaxsize-iminsize)+ii, hash);
         }
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_LAYER1] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            print_local();

            int jmaxglobal = (jmax+1)*IPOW2(levmx);
            int imaxglobal = (imax+1)*IPOW2(levmx);
            fprintf(fp,"\n                                    HASH numbering for 1 layer\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        fprintf(fp,"%5d",read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) );
                     } else {
                        fprintf(fp,"     ");
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif
         }

         // Layer 2
#ifdef _OPENMP
#pragma omp master
         {
#endif
         for (int ic =0; ic<nbsize_local; ic++){
            if (border_cell_needed_local[ic] > 0) continue;
            int jj = border_cell_j_local[ic];
            int ii = border_cell_i_local[ic];
            int lev = border_cell_level_local[ic];
            int levmult = IPOW2(levmx-lev);

            int iicur = ii*levmult-iminsize;
            int iilft = max( (ii-1)*levmult, 0         )-iminsize;
            int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
            int jjcur = jj*levmult-jminsize;
            int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
            int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

            //fprintf(fp,"            DEBUG layer2 ic %d num %d i %d j %d lev %d\n",ic,border_cell_num_local[ic],ii,jj,lev);
   
            int iborder = 0;

            // Test for cell to left
            if (iicur-(iicur-iilft)/2 >= 0 && iicur-(iicur-iilft)/2 < imaxsize-iminsize && jjcur >= 0 &&      (jjcur+jjtop)/2 < jmaxsize-jminsize){
               // Check for finer cell left and bottom side
               if (lev != levmx){                                // finer neighbor
                  int iilftfiner = iicur-(iicur-iilft)/2;
                  int nl = read_hash(jjcur*(imaxsize-iminsize)+iilftfiner, hash);
                  if (nl >= (int)(ncells+noffset) && (border_cell_needed_local[nl-ncells-noffset] & 0x0001) == 0x0001) {
                     iborder = 0x0001;
                  } else {
                     // Also check for finer cell left and top side
                     int jjtopfiner = (jjcur+jjtop)/2;
                     int nlt = read_hash(jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
                     if ( nlt >= (int)(ncells+noffset) && (border_cell_needed_local[nlt-ncells-noffset] & 0x0001) == 0x0001) {
                        iborder = 0x0001;
                     }
                  }
               }
               if ( (iborder & 0x0001) == 0 && iilft >= 0) { //same size
                  int nl = read_hash(jjcur*(imaxsize-iminsize)+iilft, hash);
                  int levcheck = -1;
                  if (nl-noffset >= 0 && nl-noffset < (int)ncells) {
                     levcheck = level[nl-noffset];
                  } else if (nl >= 0 && (int)(nl-ncells-noffset) >= 0 && (int)(nl-ncells-noffset) < nbsize_local) {
                     levcheck = border_cell_level_local[nl-ncells-noffset];
                  }
                  if (nl >= (int)(ncells+noffset) && levcheck == lev && (border_cell_needed_local[nl-ncells-noffset] & 0x0001) == 0x0001) {
                     iborder = 0x0001;
                  } else if (lev != 0 && iilft-(iicur-iilft) >= 0){      // coarser neighbor
                     iilft -= iicur-iilft;
                     int jjlft = (jj/2)*2*levmult-jminsize;
                     nl = read_hash(jjlft*(imaxsize-iminsize)+iilft, hash);
                     levcheck = -1;
                     if (nl-noffset >= 0 && nl-noffset < (int)ncells) {
                        levcheck = level[nl-noffset];
                     } else if (nl >= 0 && (int)(nl-ncells-noffset) >= 0 && (int)(nl-ncells-noffset) < nbsize_local) {
                        levcheck = border_cell_level_local[nl-ncells-noffset];
                     }
                     // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
                     if (nl  >= (int)(ncells+noffset) && levcheck == lev-1 && (border_cell_needed_local[nl-ncells-noffset] & 0x0001) == 0x0001) {
                        iborder = 0x0001;
                     }
                  }
               }
            }

            // Test for cell to right
            if (iirht < imaxsize-iminsize && iirht >= 0 && jjcur >= 0 && jjtop < jmaxsize-jminsize) {
               // right neighbor -- finer, same size and coarser
               int nr = read_hash(jjcur*(imaxsize-iminsize)+iirht, hash);
               if (nr >= (int)(ncells+noffset) && (border_cell_needed_local[nr-ncells-noffset] & 0x0002) == 0x0002) {
                  iborder = 0x0002;
               } else if (lev != levmx){
                  // right neighbor -- finer right top test
                  int jjtopfiner = (jjcur+jjtop)/2;
                  int nrt = read_hash(jjtopfiner*(imaxsize-iminsize)+iirht, hash);
                  if (nrt >= (int)(ncells+noffset) && (border_cell_needed_local[nrt-ncells-noffset] & 0x0002) == 0x0002) {
                     iborder = 0x0002;
                  }
               }
               if ( (iborder & 0x0002) == 0  && lev != 0) { // test for coarser, but not directly right
                  int jjrhtcoarser = (jj/2)*2*levmult-jminsize;
                  if (jjrhtcoarser != jjcur) {
                     int nr = read_hash(jjrhtcoarser*(imaxsize-iminsize)+iirht, hash);
                     int levcheck = -1;
                     if (nr-noffset >= 0 && nr-noffset < (int)ncells) {
                        levcheck = level[nr-noffset];
                     } else if (nr >= 0 && (int)(nr-ncells-noffset) >= 0 && (int)(nr-ncells-noffset) < nbsize_local) {
                        levcheck = border_cell_level_local[nr-ncells-noffset];
                     }
                     if (nr >= (int)(ncells+noffset) && levcheck == lev-1 && (border_cell_needed_local[nr-ncells-noffset] & 0x0002) == 0x0002) {
                        iborder = 0x0002;
                     }
                  }
               }
            }

            // Test for cell to bottom
            if (iicur >= 0 && (iicur+iirht)/2 < imaxsize-iminsize && jjcur-(jjcur-jjbot)/2 >= 0 && jjcur-(jjcur-jjbot)/2 < jmaxsize-jminsize){
               // Check for finer cell below and left side
               if (lev != levmx){                                // finer neighbor
                  int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                  int nb = read_hash(jjbotfiner*(imaxsize-iminsize)+iicur, hash);
                  if (nb >= (int)(ncells+noffset) && (border_cell_needed_local[nb-ncells-noffset] & 0x0004) == 0x0004) {
                     iborder = 0x0004;
                  } else {
                     // Also check for finer cell below and right side
                     int iirhtfiner = (iicur+iirht)/2;
                     int nbr = read_hash(jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
                     if (nbr >= (int)(ncells+noffset) && (border_cell_needed_local[nbr-ncells-noffset] & 0x0004) == 0x0004) {
                        iborder = 0x0004;
                     }
                  }
               }
               if ( (iborder & 0x0004) == 0 && jjbot >= 0) { //same size
                  int nb = read_hash(jjbot*(imaxsize-iminsize)+iicur, hash);
                  int levcheck = -1;
                  if (nb-noffset >= 0 && nb-noffset < (int)ncells) {
                     levcheck = level[nb-noffset];
                  } else if (nb >= 0 && (int)(nb-ncells-noffset) >= 0 && (int)(nb-ncells-noffset) < nbsize_local) {
                     levcheck = border_cell_level_local[nb-ncells-noffset];
                  }
                  if (nb >= (int)(ncells+noffset) && levcheck == lev && (border_cell_needed_local[nb-ncells-noffset] & 0x0004) == 0x0004) {
                     iborder = 0x0004;
                  } else if (lev != 0 && jjbot-(jjcur-jjbot) >= 0){      // coarser neighbor
                     jjbot -= jjcur-jjbot;
                     int iibot = (ii/2)*2*levmult-iminsize;
                     nb = read_hash(jjbot*(imaxsize-iminsize)+iibot, hash);
                     levcheck = -1;
                     if (nb-noffset >= 0 && nb-noffset < (int)ncells) {
                        levcheck = level[nb-noffset];
                     } else if (nb >= 0 && (int)(nb-ncells-noffset) >= 0 && (int)(nb-ncells-noffset) < nbsize_local) {
                        levcheck = border_cell_level_local[nb-ncells-noffset];
                     }
                     // we have to test for coarser level or it could be a same size cell one or two cells away that it is matching
                     if (nb >= (int)(ncells+noffset) && levcheck == lev-1 && (border_cell_needed_local[nb-ncells-noffset] & 0x0004) == 0x0004) {
                        iborder = 0x0004;
                     }
                  }
               }
            }

            // Test for cell to top
            if (iirht < imaxsize-iminsize && iicur >= 0 && jjtop >= 0 && jjtop < jmaxsize-jminsize) {
               // top neighbor -- finer, same size and coarser
               int nt = read_hash(jjtop*(imaxsize-iminsize)+iicur, hash);
               if (nt  >= (int)(ncells+noffset) && (border_cell_needed_local[nt-ncells-noffset] & 0x0008) == 0x0008) {
                  iborder = 0x0008;
               } else if (lev != levmx){
                  int iirhtfiner = (iicur+iirht)/2;
                  int ntr = read_hash(jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
                  if ( ntr >= (int)(ncells+noffset) && (border_cell_needed_local[ntr-ncells-noffset] & 0x0008) == 0x0008) {
                     iborder = 0x0008;
                  }
               }
               if ( (iborder & 0x0008) == 0  && lev != 0) { // test for coarser, but not directly above
                  int iitopcoarser = (ii/2)*2*levmult-iminsize;
                  if (iitopcoarser != iicur) {
                     int nb = read_hash(jjtop*(imaxsize-iminsize)+iitopcoarser, hash);
                     int levcheck = -1;
                     if (nb-noffset >= 0 && nb-noffset < (int)ncells) {
                        levcheck = level[nb-noffset];
                     } else if (nb >= 0 && (int)(nb-ncells-noffset) >= 0 && (int)(nb-ncells-noffset) < nbsize_local) {
                        levcheck = border_cell_level_local[nb-ncells-noffset];
                     }
                     if (nb-noffset >= (int)(ncells-noffset) && levcheck == lev-1 && (border_cell_needed_local[nb-ncells-noffset] & 0x0008) == 0x0008) {
                        iborder = 0x0008;
                     }
                  }
               }
            }

            if (iborder) border_cell_needed_local[ic] = iborder |= 0x0016;
         }
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif

         vector<int> indices_needed;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         if (DEBUG) {
            for(int ic=0; ic<nbsize_local; ic++){
               if (border_cell_needed_local[ic] <  0x0016) fprintf(fp,"%d: First  set of needed cells ic %3d cell %3d type %3d\n",mype,ic,border_cell_num_local[ic],border_cell_needed_local[ic]);
               if (border_cell_needed_local[ic] >= 0x0016) fprintf(fp,"%d: Second set of needed cells ic %3d cell %3d type %3d\n",mype,ic,border_cell_num_local[ic],border_cell_needed_local[ic]);
            }
         }

         int inew = 0;
         for(int ic=0; ic<nbsize_local; ic++){
            if (border_cell_needed_local[ic] <= 0) continue;
            indices_needed.push_back(border_cell_num_local[ic]);

            border_cell_num_local[inew]    = border_cell_num_local[ic];
            border_cell_i_local[inew]      = border_cell_i_local[ic];
            border_cell_j_local[inew]      = border_cell_j_local[ic];
            border_cell_level_local[inew]  = border_cell_level_local[ic];
            // border_cell_num_local is not used after -- could be commented out?
            // border_cell_needed_local[inew] = 1;

            inew++;
         }
         nbsize_local = inew;

         free(border_cell_num_local);
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif

         // Walk through cell array and set hash to global cell values
         //fprintf(fp,"%d: DEBUG new hash jminsize %d jmaxsize %d iminsize %d imaxsize %d\n",mype,jminsize,jmaxsize,iminsize,imaxsize);
#ifdef _OPENMP
#pragma omp for
#endif
         for(int ic=0; ic<nbsize_local; ic++){
            int lev = border_cell_level_local[ic];
            int levmult = IPOW2(levmx-lev);

            int ii = border_cell_i_local[ic]*levmult-iminsize;
            int jj = border_cell_j_local[ic]*levmult-jminsize;

            write_hash(-(ncells+ic), jj*(imaxsize-iminsize)+ii, hash);
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_LAYER2] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            print_local();

            int jmaxglobal = (jmax+1)*IPOW2(levmx);
            int imaxglobal = (imax+1)*IPOW2(levmx);
            fprintf(fp,"\n                                    HASH numbering for 2 layer\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        fprintf(fp,"%5d",read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) );
                     } else {
                        fprintf(fp,"     ");
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");
#ifdef _OPENMP
         } // end master region
#endif
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_LAYER_LIST] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         int nghost = nbsize_local;
         ncells_ghost = ncells + nghost;


#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         // needs to cast char_t to void so doesn't mistake it for string
         celltype = (char_t *)mesh_memory.memory_realloc(ncells_ghost, (void *)celltype);
         i        = (int *)mesh_memory.memory_realloc(ncells_ghost, i);
         j        = (int *)mesh_memory.memory_realloc(ncells_ghost, j);
         level    = (uchar_t *)mesh_memory.memory_realloc(ncells_ghost, level);
         nlft     = (int *)mesh_memory.memory_realloc(ncells_ghost, nlft);
         nrht     = (int *)mesh_memory.memory_realloc(ncells_ghost, nrht);
         nbot     = (int *)mesh_memory.memory_realloc(ncells_ghost, nbot);
         ntop     = (int *)mesh_memory.memory_realloc(ncells_ghost, ntop);
         memory_reset_ptrs();
#ifdef _OPENMP
         } // end master region
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
         for (int ic = ncells; ic < (int)ncells_ghost; ic++){
            nlft[ic] = -1;
            nrht[ic] = -1;
            nbot[ic] = -1;
            ntop[ic] = -1;
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_COPY_MESH_DATA] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

#ifdef _OPENMP
#pragma omp for
#endif
         for(int ic=0; ic<nbsize_local; ic++){
            int ii = border_cell_i_local[ic];
            int jj = border_cell_j_local[ic];
            int lev = border_cell_level_local[ic];
            if (ii < lev_ibegin[lev]) celltype[ncells+ic] = LEFT_BOUNDARY;
            if (ii > lev_iend[lev])   celltype[ncells+ic] = RIGHT_BOUNDARY;
            if (jj < lev_jbegin[lev]) celltype[ncells+ic] = BOTTOM_BOUNDARY;
            if (jj > lev_jend[lev])   celltype[ncells+ic] = TOP_BOUNDARY;
            i[ncells+ic]     = ii;
            j[ncells+ic]     = jj;
            level[ncells+ic] = lev;
         }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            free(border_cell_i_local);
            free(border_cell_j_local);
            free(border_cell_level_local);
#ifdef _OPENMP
         } // end master region
#endif

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_FILL_MESH_GHOST] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               fprintf(fp,"After copying i,j, level to ghost cells\n");
               print_local();
#ifdef _OPENMP
            } // end master region
#endif
         }

#ifdef _OPENMP
#pragma omp for
#endif
         for (uint ic=0; ic<ncells_ghost; ic++){
            int ii = i[ic];
            int jj = j[ic];
            int lev = level[ic];
            int levmult = IPOW2(levmx-lev);

            int iicur = ii*levmult-iminsize;
            int iilft = max( (ii-1)*levmult, 0         )-iminsize;
            int iirht = min( (ii+1)*levmult, imaxcalc-1)-iminsize;
            int jjcur = jj*levmult-jminsize;
            int jjbot = max( (jj-1)*levmult, 0         )-jminsize;
            int jjtop = min( (jj+1)*levmult, jmaxcalc-1)-jminsize;

            //fprintf(fp,"DEBUG neigh ic %d nlft %d ii %d levmult %d iminsize %d icheck %d\n",ic,nlft[ic],ii,levmult,iminsize,(max(  ii   *levmult-1, 0))-iminsize);

            int nlftval = nlft[ic];
            int nrhtval = nrht[ic];
            int nbotval = nbot[ic];
            int ntopval = ntop[ic];

            if (nlftval == -1){
               // Taking care of boundary cells
               // Force each boundary cell to point to itself on its boundary direction
               if (iicur <    1*IPOW2(levmx)  -iminsize) nlftval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);

               // Boundary cells next to corner boundary need special checks
               if (iicur ==    1*IPOW2(levmx)-iminsize &&  (jjcur < 1*IPOW2(levmx)-jminsize || jjcur >= jmax*IPOW2(levmx)-jminsize ) ) nlftval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);

               // need to check for finer neighbor first
               // Right and top neighbor don't change for finer, so drop through to same size
               // Left and bottom need to be half of same size index for finer test
               if (lev != levmx) {
                  int iilftfiner = iicur-(iicur-iilft)/2;
                  if (nlftval == -1 && iilftfiner >= 0) nlftval = read_hash(jjcur*(imaxsize-iminsize)+iilftfiner, hash);
               }

               // same size neighbor
               if (nlftval == -1 && iilft >= 0) nlftval = read_hash(jjcur*(imaxsize-iminsize)+iilft, hash);

               // Now we need to take care of special case where bottom and left boundary need adjustment since
               // expected cell doesn't exist on these boundaries if it is finer than current cell
               if (jjcur < 1*IPOW2(levmx) && lev != levmx) {
                  if (nlftval == -1) {
                     int iilftfiner = iicur-(iicur-iilft)/2;
                     int jjtopfiner = (jjcur+jjtop)/2;
                     if (jjtopfiner < jmaxsize-jminsize && iilftfiner >= 0) nlftval = read_hash(jjtopfiner*(imaxsize-iminsize)+iilftfiner, hash);
                  }
               }

               // coarser neighbor
               if (lev != 0){
                  if (nlftval == -1) {
                     int iilftcoarser = iilft - (iicur-iilft);
                     int jjlft = (jj/2)*2*levmult-jminsize;
                     if (iilftcoarser >=0) nlftval = read_hash(jjlft*(imaxsize-iminsize)+iilftcoarser, hash);
                  }
               }

               if (nlftval != -1) nlft[ic] = nlftval;
            }

            if (nrhtval == -1) {
               // Taking care of boundary cells
               // Force each boundary cell to point to itself on its boundary direction
               if (iicur > imax*IPOW2(levmx)-1-iminsize) nrhtval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);

               // Boundary cells next to corner boundary need special checks
               if (iirht == imax*IPOW2(levmx)-iminsize &&  (jjcur < 1*IPOW2(levmx)-jminsize || jjcur >= jmax*IPOW2(levmx)-jminsize ) ) nrhtval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);

               // same size neighbor
               if (nrhtval == -1 && iirht < imaxsize-iminsize) nrhtval = read_hash(jjcur*(imaxsize-iminsize)+iirht, hash);

               // Now we need to take care of special case where bottom and left boundary need adjustment since
               // expected cell doesn't exist on these boundaries if it is finer than current cell
               if (jjcur < 1*IPOW2(levmx) && lev != levmx) {
                  if (nrhtval == -1) {
                     int jjtopfiner = (jjcur+jjtop)/2;
                     if (jjtopfiner < jmaxsize-jminsize && iirht < imaxsize-iminsize) nrhtval = read_hash(jjtopfiner*(imaxsize-iminsize)+iirht, hash);
                  }
               }

               // coarser neighbor
               if (lev != 0){
                  if (nrhtval == -1) {
                     int jjrht = (jj/2)*2*levmult-jminsize;
                     if (iirht < imaxsize-iminsize) nrhtval = read_hash(jjrht*(imaxsize-iminsize)+iirht, hash);
                  }
               }
               if (nrhtval != -1) nrht[ic] = nrhtval;
            }
 
            if (nbotval == -1) {
               // Taking care of boundary cells
               // Force each boundary cell to point to itself on its boundary direction
               if (jjcur <    1*IPOW2(levmx)  -jminsize) nbotval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);
               // Boundary cells next to corner boundary need special checks
               if (jjcur ==    1*IPOW2(levmx)-jminsize &&  (iicur < 1*IPOW2(levmx)-iminsize || iicur >= imax*IPOW2(levmx)-iminsize ) ) nbotval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);

               // need to check for finer neighbor first
               // Right and top neighbor don't change for finer, so drop through to same size
               // Left and bottom need to be half of same size index for finer test
               if (lev != levmx) {
                  int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                  if (nbotval == -1 && jjbotfiner >= 0) nbotval = read_hash(jjbotfiner*(imaxsize-iminsize)+iicur, hash);
               }

               // same size neighbor
               if (nbotval == -1 && jjbot >=0) nbotval = read_hash(jjbot*(imaxsize-iminsize)+iicur, hash);

               // Now we need to take care of special case where bottom and left boundary need adjustment since
               // expected cell doesn't exist on these boundaries if it is finer than current cell
               if (iicur < 1*IPOW2(levmx) && lev != levmx) {
                  if (nbotval == -1) {
                     int iirhtfiner = (iicur+iirht)/2;
                     int jjbotfiner = jjcur-(jjcur-jjbot)/2;
                     if (jjbotfiner >= 0 && iirhtfiner < imaxsize-iminsize) nbotval = read_hash(jjbotfiner*(imaxsize-iminsize)+iirhtfiner, hash);
                  }
               }

               // coarser neighbor
               if (lev != 0){
                  if (nbotval == -1) {
                     int jjbotcoarser = jjbot - (jjcur-jjbot);
                     int iibot = (ii/2)*2*levmult-iminsize;
                     if (jjbotcoarser >= 0 && iibot >= 0) nbotval = read_hash(jjbotcoarser*(imaxsize-iminsize)+iibot, hash);
                  }
               }
               if (nbotval != -1) nbot[ic] = nbotval;
            }
    
            if (ntopval == -1) {
               // Taking care of boundary cells
               // Force each boundary cell to point to itself on its boundary direction
               if (jjcur > jmax*IPOW2(levmx)-1-jminsize) ntopval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);
               // Boundary cells next to corner boundary need special checks
               if (jjtop == jmax*IPOW2(levmx)-jminsize &&  (iicur < 1*IPOW2(levmx)-iminsize || iicur >= imax*IPOW2(levmx)-iminsize ) ) ntopval = read_hash(jjcur*(imaxsize-iminsize)+iicur, hash);

               // same size neighbor
               if (ntopval == -1 && jjtop < jmaxsize-jminsize) ntopval = read_hash(jjtop*(imaxsize-iminsize)+iicur, hash);
   
               if (iicur < 1*IPOW2(levmx)) {
                  if (ntopval == -1) {
                     int iirhtfiner = (iicur+iirht)/2;
                     if (jjtop < jmaxsize-jminsize && iirhtfiner < imaxsize-iminsize) ntopval = read_hash(jjtop*(imaxsize-iminsize)+iirhtfiner, hash);
                  }
               }
   
               // coarser neighbor
               if (lev != 0){
                  if (ntopval == -1) {
                     int iitop = (ii/2)*2*levmult-iminsize;
                     if (jjtop < jmaxsize-jminsize && iitop < imaxsize-iminsize) ntopval = read_hash(jjtop*(imaxsize-iminsize)+iitop, hash);
                  }
               }
               if (ntopval != -1) ntop[ic] = ntopval;
            }
 
            //fprintf(fp,"%d: neighbors[%d] = %d %d %d %d\n",mype,ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
         }

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_FILL_NEIGH_GHOST] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
            fprintf(fp,"After setting neighbors through ghost cells\n");
            print_local();
#ifdef _OPENMP
            } // end master region
#endif
         }

/*
         // Set neighbors to global cell numbers from hash
         for (uint ic=0; ic<ncells; ic++){
            ii = i[ic];
            jj = j[ic];
            lev = level[ic];
            levmult = IPOW2(levmx-lev);
            //fprintf(fp,"%d:Neighbors input for ic %d ii %d jj %d levmult %d lev %d\n",mype,ic, ii, jj, levmult,lev);
            //fprintf(fp,"%d:Neighbors befor ic %d nlft %d nrht %d nbot %d ntop %d\n",mype,ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            if (nlft[ic] == -1) nlft[ic] = hash[(      jj   *levmult               )-jminsize][(max(  ii   *levmult-1, 0         ))-iminsize];
            if (celltype[ic] == BOTTOM_BOUNDARY && nlft[ic] == -1){
               if (nlft[ic] == -1) nlft[ic] = hash[(jj+1)*levmult-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
            }
            if (nrht[ic] == -1) nrht[ic] = hash[(      jj   *levmult               )-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
            if (celltype[ic] == BOTTOM_BOUNDARY && nrht[ic] == -1){
               if (nrht[ic] == -1) nrht[ic] = hash[(jj+1)*levmult-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
               //if (ic == 3 && mype == 0) printf("DEBUG line %d -- ic %d celltype %d nrht %d\n",__line__,ic,celltype[ic],nrht[ic]);
               //printf("DEBUG line %d -- ic %d celltype %d nrht %d jj %d ii %d\n",__line__,ic,celltype[ic],nrht[ic],(jj+1)*levmult-jminsize,(min( (ii+1)*levmult,   imaxcalc-1))-iminsize);
            }
            if (nbot[ic] == -1) nbot[ic] = hash[(max(  jj   *levmult-1, 0)         )-jminsize][(      ii   *levmult               )-iminsize];
            if (celltype[ic] == LEFT_BOUNDARY && nbot[ic] == -1){
               if (nbot[ic] == -1) nbot[ic] = hash[(max(  jj   *levmult-1, 0)         )-jminsize][(      ii   *levmult+1             )-iminsize];
            }
            if (ntop[ic] == -1) ntop[ic] = hash[(min( (jj+1)*levmult,   jmaxcalc-1))-jminsize][(      ii   *levmult               )-iminsize];
            if (celltype[ic] == LEFT_BOUNDARY && ntop[ic] == -1){
               if (ntop[ic] == -1) ntop[ic] = hash[(min( (jj+1)*levmult,   jmaxcalc-1))-jminsize][(      ii   *levmult+1             )-iminsize];
            }
            //fprintf(fp,"%d:Neighbors after ic %d nlft %d nrht %d nbot %d ntop %d\n",mype,ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
         }
*/

         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_SET_CORNER_NEIGH] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               fprintf(fp,"After setting corner neighbors\n");
               print_local();
#ifdef _OPENMP
            } // end master region
#endif
         }

         // Adjusting neighbors to local indices
#ifdef _OPENMP
#pragma omp for
#endif
         for (uint ic=0; ic<ncells_ghost; ic++){
            //fprintf(fp,"%d: ic %d nlft %d noffset %d ncells %ld\n",mype,ic,nlft[ic],noffset,ncells);
            if (nlft[ic] <= -(int)ncells && nlft[ic] > -(int)ncells_ghost){
               nlft[ic] = abs(nlft[ic]);
            } else if (nlft[ic] >= noffset && nlft[ic] < (int)(noffset+ncells)) {
               nlft[ic] -= noffset;
            }
            if (nrht[ic] <= -(int)ncells && nrht[ic] > -(int)ncells_ghost){
               nrht[ic] = abs(nrht[ic]);
            } else if (nrht[ic] >= noffset && nrht[ic] < (int)(noffset+ncells)) {
               nrht[ic] -= noffset;
            }
            if (nbot[ic] <= -(int)ncells && nbot[ic] > -(int)ncells_ghost){
               nbot[ic] = abs(nbot[ic]);
            } else if (nbot[ic] >= noffset && nbot[ic] < (int)(noffset+ncells)) {
               nbot[ic] -= noffset;
            }
            if (ntop[ic] <= -(int)ncells && ntop[ic] > -(int)ncells_ghost){
               ntop[ic] = abs(ntop[ic]);
            } else if (ntop[ic] >= noffset && ntop[ic] < (int)(noffset+ncells)) {
               ntop[ic] -= noffset;
            }
         }

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
               fprintf(fp,"After adjusting neighbors to local indices\n");
               print_local();
#ifdef _OPENMP
            } // end master region
#endif
         }
         
         if (TIMING_LEVEL >= 2) {
#ifdef _OPENMP
#pragma omp master
#endif
            cpu_timers[MESH_TIMER_NEIGH_ADJUST] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
            offtile_ratio_local = (offtile_ratio_local*(double)offtile_local_count) + ((double)nghost / (double)ncells);
            offtile_local_count++;
            offtile_ratio_local /= offtile_local_count;
            //printf("%d ncells size is %ld ncells_ghost size is %ld nghost %d\n",mype,ncells,ncells_ghost,nghost);
            //fprintf(fp,"%d ncells_ghost size is %ld nghost %d\n",mype,ncells_ghost,nghost);

            if (cell_handle) L7_Free(&cell_handle);
            cell_handle=0;

            if (DEBUG) {
               fprintf(fp,"%d: SETUP ncells %ld noffset %d nghost %d\n",mype,ncells,noffset,nghost);
               for (int ig = 0; ig<nghost; ig++){
                  fprintf(fp,"%d: indices needed ic %d index %d\n",mype,ig,indices_needed[ig]);
               }
            }
            L7_Setup(0, noffset, ncells, &indices_needed[0], nghost, &cell_handle);

            if (TIMING_LEVEL >= 2) cpu_timers[MESH_TIMER_SETUP_COMM] += cpu_timer_stop(tstart_lev2);

#ifdef _OPENMP
         } // end master region
#endif

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
            print_local();

            int jmaxglobal = (jmax+1)*IPOW2(levmx);
            int imaxglobal = (imax+1)*IPOW2(levmx);
            fprintf(fp,"\n                                    HASH numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        fprintf(fp,"%5d",read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) );
                     } else {
                        fprintf(fp,"     ");
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");

            fprintf(fp,"\n                                    nlft numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) -noffset;
                        if ( (hashval >= 0 && hashval < (int)ncells) ) {
                              fprintf(fp,"%5d",nlft[hashval]);
                        } else {
                              fprintf(fp,"     ");
                        }
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");
      
            fprintf(fp,"\n                                    nrht numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if ( ii >= iminsize && ii < imaxsize ) {
                        int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) -noffset;
                        if ( hashval >= 0 && hashval < (int)ncells ) {
                           fprintf(fp,"%5d",nrht[hashval]);
                        } else {
                           fprintf(fp,"     ");
                        }
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");

            fprintf(fp,"\n                                    nbot numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if ( ii >= iminsize && ii < imaxsize ) {
                        int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) -noffset;
                        if ( hashval >= 0 && hashval < (int)ncells ) {
                           fprintf(fp,"%5d",nbot[hashval]);
                        } else {
                           fprintf(fp,"     ");
                        }
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");

            fprintf(fp,"\n                                    ntop numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if ( ii >= iminsize && ii < imaxsize ) {
                        int hashval = read_hash((jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), hash) -noffset;
                        if ( hashval >= 0 && hashval < (int)ncells ) {
                           fprintf(fp,"%5d",ntop[hashval]);
                        } else {
                           fprintf(fp,"     ");
                        }
                     }
                  }
               }
               fprintf(fp,"\n");
            }
            fprintf(fp,"%2d:      ",mype);
            for (int ii = 0; ii<imaxglobal; ii++){
               fprintf(fp,"%4d:",ii);
            }
            fprintf(fp,"\n");
      
#ifdef _OPENMP
            } // end master region
#endif
         } // end DEBUG

         if (DEBUG) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
            print_local();

            for (uint ic=0; ic<ncells; ic++){
               fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
                   mype,ic,i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            }
            int ig=0;
            for (uint ic=ncells; ic<ncells_ghost; ic++, ig++){
               fprintf(fp,"%d: after  update ic %d off %d i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
                   mype,ic,indices_needed[ig],i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            }
#ifdef _OPENMP
            } // end master region
#endif
         } // end DEBUG

      } // if numpe > 1
#endif

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
            {
#endif
      write_hash_collision_report();
      read_hash_collision_report();
      compact_hash_delete(hash);

#ifdef BOUNDS_CHECK
      {
         for (uint ic=0; ic<ncells; ic++){
            int nl = nlft[ic];
            if (nl<0 || nl>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nlft %d\n",mype,__LINE__,ic,nl);
            if (level[nl] > level[ic]){
               int ntl = ntop[nl];
               if (ntl<0 || ntl>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d global %d nlft %d ntop of nlft %d\n",mype,__LINE__,ic,ic+noffset,nl,ntl);
            }
            int nr = nrht[ic];
            if (nr<0 || nr>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nrht %d\n",mype,__LINE__,ic,nr);
            if (level[nr] > level[ic]){
               int ntr = ntop[nr];
               if (ntr<0 || ntr>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d ntop of nrht %d\n",mype,__LINE__,ic,ntr);
            }
            int nb = nbot[ic];
            if (nb<0 || nb>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nbot %d\n",mype,__LINE__,ic,nb);
            if (level[nb] > level[ic]){
               int nrb = nrht[nb];
               if (nrb<0 || nrb>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nrht of nbot %d\n",mype,__LINE__,ic,nrb);
            }
            int nt = ntop[ic];
            if (nt<0 || nt>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d ntop %d\n",mype,__LINE__,ic,nt);
            if (level[nt] > level[ic]){
               int nrt = nrht[nt];
               if (nrt<0 || nrt>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nrht of ntop %d\n",mype,__LINE__,ic,nrt);
            }
         }
      }
#endif

#ifdef _OPENMP
            } // end master region
#pragma omp barrier
#endif

      } else if (calc_neighbor_type == KDTREE) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
         {
#endif
         struct timespec tstart_lev2;
         if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

         TBounds box;
         vector<int> index_list(IPOW2(levmx*levmx) );

         int num;

         ibase = 0;
         calc_spatial_coordinates(ibase);

         kdtree_setup();

         if (TIMING_LEVEL >= 2) {
            cpu_timers[MESH_TIMER_KDTREE_SETUP] += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         for (uint ic=0; ic<ncells; ic++) {

            //left
            nlft[ic]  = ic;
            box.min.x = x[ic]-0.25*dx[ic];
            box.max.x = x[ic]-0.25*dx[ic];
            box.min.y = y[ic]+0.25*dy[ic];
            box.max.y = y[ic]+0.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) nlft[ic]=index_list[0];

            //right
            nrht[ic]  = ic;
            box.min.x = x[ic]+1.25*dx[ic];
            box.max.x = x[ic]+1.25*dx[ic];
            box.min.y = y[ic]+0.25*dy[ic];
            box.max.y = y[ic]+0.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) nrht[ic]=index_list[0];

            //bot
            nbot[ic]  = ic;
            box.min.x = x[ic]+0.25*dx[ic];
            box.max.x = x[ic]+0.25*dx[ic];
            box.min.y = y[ic]-0.25*dy[ic];
            box.max.y = y[ic]-0.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) nbot[ic]=index_list[0];

            //top
            ntop[ic]  = ic;
            box.min.x = x[ic]+0.25*dx[ic];
            box.max.x = x[ic]+0.25*dx[ic];
            box.min.y = y[ic]+1.25*dy[ic];
            box.max.y = y[ic]+1.25*dy[ic];
            KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
            if (num == 1) ntop[ic]=index_list[0];
         }  //  End main loop over cells.

         KDTree_Destroy(&tree);

         if (TIMING_LEVEL >= 2) cpu_timers[MESH_TIMER_KDTREE_QUERY] += cpu_timer_stop(tstart_lev2);

#ifdef _OPENMP
         }
#pragma omp barrier
#endif
      } // calc_neighbor_type

   }

#ifdef _OPENMP
#pragma omp master
#endif
   cpu_timers[MESH_TIMER_CALC_NEIGHBORS] += cpu_timer_stop(tstart_cpu);
}

#ifdef HAVE_OPENCL
void Mesh::gpu_calc_neighbors(void)
{
   if (! gpu_do_rezone) return;

   ulong gpu_hash_table_size =  0;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_lev2;
   cpu_timer_start(&tstart_lev2);

   cl_command_queue command_queue = ezcl_get_command_queue();

   gpu_counters[MESH_COUNTER_CALC_NEIGH]++;

   assert(dev_levtable);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);

   size_t mem_request = (int)((float)ncells*mem_factor);

   size_t local_work_size = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;

   //printf("DEBUG file %s line %d dev_nlft %p size %d\n",__FILE__,__LINE__,dev_nlft,ezcl_get_device_mem_nelements(dev_nlft));

   if (dev_nlft == NULL || ezcl_get_device_mem_nelements(dev_nlft) < (int)ncells) {
      dev_nlft     = ezcl_malloc(NULL, const_cast<char *>("dev_nlft"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      dev_nrht     = ezcl_malloc(NULL, const_cast<char *>("dev_nrht"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      dev_nbot     = ezcl_malloc(NULL, const_cast<char *>("dev_nbot"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      dev_ntop     = ezcl_malloc(NULL, const_cast<char *>("dev_ntop"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_neighbor_init,  0, sizeof(cl_int),   (void *)&ncells);
      ezcl_set_kernel_arg(kernel_neighbor_init,  1, sizeof(cl_mem),   (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_neighbor_init,  2, sizeof(cl_mem),   (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_neighbor_init,  3, sizeof(cl_mem),   (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_neighbor_init,  4, sizeof(cl_mem),   (void *)&dev_ntop);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_neighbor_init,   1, NULL, &global_work_size, &local_work_size, NULL);
   }

   int imaxsize = (imax+1)*IPOW2(levmx);
   int jmaxsize = (jmax+1)*IPOW2(levmx);

   int gpu_hash_method       = METHOD_UNSET;
// allow input.c to control hash types and methods
   if (choose_hash_method != METHOD_UNSET) gpu_hash_method = choose_hash_method;
//=========

   size_t hashsize;

   uint hash_report_level = 1;
   cl_mem dev_hash_header = NULL;
   cl_mem dev_hash = gpu_compact_hash_init(ncells, imaxsize, jmaxsize, gpu_hash_method, hash_report_level,
      &gpu_hash_table_size, &hashsize, &dev_hash_header);

      /*
                    const int   isize,        // 0
                    const int   levmx,        // 1
                    const int   imaxsize,     // 2
           __global const int   *levtable,    // 3
           __global const uchar_t *level,       // 4
           __global const int   *i,           // 5
           __global const int   *j,           // 6
           __global const ulong *hash_header, // 7
           __global       int   *hash)        // 8
      */

   cl_event hash_setup_event;

   ezcl_set_kernel_arg(kernel_hash_setup,  0, sizeof(cl_int),   (void *)&ncells);
   ezcl_set_kernel_arg(kernel_hash_setup,  1, sizeof(cl_int),   (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_setup,  2, sizeof(cl_int),   (void *)&imaxsize);
   ezcl_set_kernel_arg(kernel_hash_setup,  3, sizeof(cl_mem),   (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_setup,  4, sizeof(cl_mem),   (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_hash_setup,  5, sizeof(cl_mem),   (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_hash_setup,  6, sizeof(cl_mem),   (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_hash_setup,  7, sizeof(cl_mem),   (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_hash_setup,  8, sizeof(cl_mem),   (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_hash_setup,  9, sizeof(cl_mem),   (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_hash_setup, 10, sizeof(cl_mem),   (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_hash_setup, 11, sizeof(cl_mem),   (void *)&dev_hash_header);
   ezcl_set_kernel_arg(kernel_hash_setup, 12, sizeof(cl_mem),   (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup,   1, NULL, &global_work_size, &local_work_size, &hash_setup_event);

   ezcl_wait_for_events(1, &hash_setup_event);
   ezcl_event_release(hash_setup_event);

   if (TIMING_LEVEL >= 2) {
      gpu_timers[MESH_TIMER_HASH_SETUP] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      cpu_timer_start(&tstart_lev2);
   }

      /*
                    const int   isize,        // 0
                    const int   levmx,        // 1
                    const int   imax,         // 2
                    const int   jmax,         // 3
                    const int   imaxsize,     // 4
                    const int   jmaxsize,     // 5
           __global const int   *levtable,    // 6
           __global const uchar_t *level,       // 7
           __global const int   *i,           // 8
           __global const int   *j,           // 9
           __global       int   *nlft,        // 10
           __global       int   *nrht,        // 11
           __global       int   *nbot,        // 12
           __global       int   *ntop,        // 13
           __global const ulong *hash_header, // 14
           __global       int   *hash)        // 15
      */

   cl_event calc_neighbors_event;

   ezcl_set_kernel_arg(kernel_calc_neighbors, 0,  sizeof(cl_int),   (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 1,  sizeof(cl_int),   (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 2,  sizeof(cl_int),   (void *)&imax);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 3,  sizeof(cl_int),   (void *)&jmax);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 4,  sizeof(cl_int),   (void *)&imaxsize);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 5,  sizeof(cl_int),   (void *)&jmaxsize);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 6,  sizeof(cl_mem),   (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 7,  sizeof(cl_mem),   (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 8,  sizeof(cl_mem),   (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 9,  sizeof(cl_mem),   (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 10, sizeof(cl_mem),   (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 11, sizeof(cl_mem),   (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 12, sizeof(cl_mem),   (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 13, sizeof(cl_mem),   (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 14, sizeof(cl_mem),   (void *)&dev_hash_header);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 15, sizeof(cl_mem),   (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors,   1, NULL, &global_work_size, &local_work_size, &calc_neighbors_event);

   ezcl_wait_for_events(1, &calc_neighbors_event);
   ezcl_event_release(calc_neighbors_event);

   gpu_compact_hash_delete(dev_hash, dev_hash_header);

   if (TIMING_LEVEL >= 2) gpu_timers[MESH_TIMER_HASH_QUERY] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);

   gpu_timers[MESH_TIMER_CALC_NEIGHBORS] += (long)(cpu_timer_stop(tstart_cpu) * 1.0e9);
}


void Mesh::gpu_calc_neighbors_local(void)
{
   if (! gpu_do_rezone) return;

   ulong gpu_hash_table_size =  0;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_lev2;
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

   cl_command_queue command_queue = ezcl_get_command_queue();

   gpu_counters[MESH_COUNTER_CALC_NEIGH]++;

   ncells_ghost = ncells;

   assert(dev_levtable);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);

   size_t one = 1;
   cl_mem dev_check = ezcl_malloc(NULL, const_cast<char *>("dev_check"), &one, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

   size_t mem_request = (int)((float)ncells*mem_factor);
   dev_nlft = ezcl_malloc(NULL, const_cast<char *>("dev_nlft"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nrht = ezcl_malloc(NULL, const_cast<char *>("dev_nrht"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nbot = ezcl_malloc(NULL, const_cast<char *>("dev_nbot"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_ntop = ezcl_malloc(NULL, const_cast<char *>("dev_ntop"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   size_t local_work_size =  64;
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   //printf("DEBUG file %s line %d lws = %d gws %d bs %d ncells %d\n",__FILE__,__LINE__,
   //   local_work_size, global_work_size, block_size, ncells);
   cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_int4), CL_MEM_READ_WRITE, 0);
   cl_mem dev_sizes = ezcl_malloc(NULL, const_cast<char *>("dev_sizes"), &one, sizeof(cl_int4),  CL_MEM_READ_WRITE, 0);

#ifdef BOUNDS_CHECK
   if (ezcl_get_device_mem_nelements(dev_i) < (int)ncells || 
       ezcl_get_device_mem_nelements(dev_j) < (int)ncells ||
       ezcl_get_device_mem_nelements(dev_level) < (int)ncells ){
      printf("%d: Warning ncells %ld size dev_i %d dev_j %d dev_level %d\n",mype,ncells,ezcl_get_device_mem_nelements(dev_i),ezcl_get_device_mem_nelements(dev_j),ezcl_get_device_mem_nelements(dev_level));
   }
#endif

      /*
       __kernel void calc_hash_size_cl(
                          const int   ncells,      // 0
                          const int   levmx,       // 1
                 __global       int   *levtable,   // 2
                 __global       uchar_t *level,      // 3
                 __global       int   *i,          // 4
                 __global       int   *j,          // 5
                 __global       int4  *redscratch, // 6
                 __global       int4  *sizes,      // 7
                 __local        int4  *tile)       // 8
      */

   ezcl_set_kernel_arg(kernel_hash_size, 0, sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_hash_size, 1, sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_size, 2, sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_size, 3, sizeof(cl_mem), (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_hash_size, 4, sizeof(cl_mem), (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_hash_size, 5, sizeof(cl_mem), (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_hash_size, 6, sizeof(cl_mem), (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_hash_size, 7, sizeof(cl_mem), (void *)&dev_sizes);
   ezcl_set_kernel_arg(kernel_hash_size, 8, local_work_size*sizeof(cl_int4), NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_size,   1, NULL, &global_work_size, &local_work_size, NULL);

   if (block_size > 1) {
         /*
         __kernel void finish_reduction_minmax4_cl(
           const    int    isize,            // 0
           __global int4  *redscratch,       // 1
           __global int4  *sizes,            // 2
           __local  int4  *tile)             // 3
         */
      ezcl_set_kernel_arg(kernel_finish_hash_size, 0, sizeof(cl_int), (void *)&block_size);
      ezcl_set_kernel_arg(kernel_finish_hash_size, 1, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_finish_hash_size, 2, sizeof(cl_mem), (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_finish_hash_size, 3, local_work_size*sizeof(cl_int4), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_finish_hash_size,   1, NULL, &local_work_size, &local_work_size, NULL);
   }

   ezcl_device_memory_delete(dev_redscratch);

   cl_int sizes[4];
   ezcl_enqueue_read_buffer(command_queue, dev_sizes, CL_TRUE,  0, 1*sizeof(cl_int4), &sizes, NULL);

   int imintile = sizes[0];
   int imaxtile = sizes[1];
   int jmintile = sizes[2];
   int jmaxtile = sizes[3];

   // Expand size by 2*coarse_cells for ghost cells
   // TODO: May want to get fancier here and calc based on cell level
   int jminsize = max(jmintile-2*IPOW2(levmx),0);
   int jmaxsize = min(jmaxtile+2*IPOW2(levmx),(jmax+1)*IPOW2(levmx));
   int iminsize = max(imintile-2*IPOW2(levmx),0);
   int imaxsize = min(imaxtile+2*IPOW2(levmx),(imax+1)*IPOW2(levmx));
   //fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);

   //ezcl_enqueue_write_buffer(command_queue, dev_sizes, CL_TRUE,  0, 1*sizeof(cl_int4), &sizes, NULL);

   int gpu_hash_method       = METHOD_UNSET;
// allow imput.c to control hash types and methods
   if (choose_hash_method != METHOD_UNSET) gpu_hash_method = choose_hash_method;
//=========

   size_t hashsize;

   uint hash_report_level = 1;
   cl_mem dev_hash_header = NULL;
   cl_mem dev_hash = gpu_compact_hash_init(ncells, imaxsize-iminsize, jmaxsize-jminsize, gpu_hash_method, hash_report_level, &gpu_hash_table_size, &hashsize, &dev_hash_header);

   int csize = corners_i.size();
#ifdef BOUNDS_CHECK
   for (int ic=0; ic<csize; ic++){
      if (corners_i[ic] >= iminsize) continue;
      if (corners_j[ic] >= jminsize) continue;
      if (corners_i[ic] <  imaxsize) continue;
      if (corners_j[ic] <  jmaxsize) continue;
      if ( (corners_j[ic]-jminsize)*(imaxsize-iminsize)+(corners_i[ic]-iminsize) < 0 ||
           (corners_j[ic]-jminsize)*(imaxsize-iminsize)+(corners_i[ic]-iminsize) > (int)hashsize){
         printf("%d: Warning corners i %d j %d hash %d\n",mype,corners_i[ic],corners_j[ic],
            (corners_j[ic]-jminsize)*(imaxsize-iminsize)+(corners_i[ic]-iminsize));
      }
   }
#endif

   size_t corners_local_work_size  = MIN(csize, TILE_SIZE);
   size_t corners_global_work_size = ((csize+corners_local_work_size - 1) /corners_local_work_size) * corners_local_work_size;

   ezcl_set_kernel_arg(kernel_hash_adjust_sizes, 0, sizeof(cl_int), (void *)&csize);
   ezcl_set_kernel_arg(kernel_hash_adjust_sizes, 1, sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_adjust_sizes, 2, sizeof(cl_int), (void *)&imax);
   ezcl_set_kernel_arg(kernel_hash_adjust_sizes, 3, sizeof(cl_int), (void *)&jmax);
   ezcl_set_kernel_arg(kernel_hash_adjust_sizes, 4, sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_adjust_sizes, 5, sizeof(cl_mem), (void *)&dev_sizes);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_adjust_sizes,   1, NULL, &corners_global_work_size, &corners_local_work_size, NULL);

   if (DEBUG){
      vector<int> sizes_tmp(4);
      ezcl_enqueue_read_buffer(command_queue, dev_sizes, CL_TRUE,  0, 1*sizeof(cl_int4), &sizes_tmp[0], NULL);
      int iminsize_tmp = sizes_tmp[0];
      int imaxsize_tmp = sizes_tmp[1];
      int jminsize_tmp = sizes_tmp[2];
      int jmaxsize_tmp = sizes_tmp[3];
      fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize_tmp,imaxsize_tmp,jminsize_tmp,jmaxsize_tmp);
   }

   local_work_size = 128;
   global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;

#ifdef BOUNDS_CHECK
   {
      vector<int> i_tmp(ncells);
      vector<int> j_tmp(ncells);
      vector<uchar_t> level_tmp(ncells);
      ezcl_enqueue_read_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int), &i_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int), &j_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, ncells*sizeof(cl_uchar_t), &level_tmp[0], NULL);
      for (int ic=0; ic<(int)ncells; ic++){
         uchar_t lev = level_tmp[ic];
         for (   int jj = j_tmp[ic]*IPOW2(levmx-lev)-jminsize; jj < (j_tmp[ic]+1)*IPOW2(levmx-lev)-jminsize; jj++) {
            for (int ii = i_tmp[ic]*IPOW2(levmx-lev)-iminsize; ii < (i_tmp[ic]+1)*IPOW2(levmx-lev)-iminsize; ii++) {
               if (jj < 0 || jj >= (jmaxsize-jminsize) || ii < 0 || ii >= (imaxsize-iminsize) ) {
                  printf("%d: Warning ncell %d writes to hash out-of-bounds at line %d ii %d jj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,ic,__LINE__,ii,jj,iminsize,imaxsize,jminsize,jmaxsize);
               }
            }
         }
      }
   }
#endif

   //printf("%d: lws %d gws %d \n",mype,local_work_size,global_work_size);
   cl_event hash_setup_local_event;

      /*
                    const int   isize,           // 0
                    const int   levmx,           // 1
                    const int   imax,            // 2
                    const int   jmax,            // 3
                    const int   noffset,         // 4
           __global       int   *sizes,          // 5
           __global       int   *levtable,       // 6
           __global       uchar_t *level,          // 7
           __global       int   *i,              // 8
           __global       int   *j,              // 9
           __global const ulong *hash_heaer,     // 10
           __global       int   *hash)           // 11
      */

   ezcl_set_kernel_arg(kernel_hash_setup_local,  0, sizeof(cl_int),   (void *)&ncells);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  1, sizeof(cl_int),   (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  2, sizeof(cl_int),   (void *)&imax);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  3, sizeof(cl_int),   (void *)&jmax);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  4, sizeof(cl_int),   (void *)&noffset);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  5, sizeof(cl_mem),   (void *)&dev_sizes);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  6, sizeof(cl_mem),   (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  7, sizeof(cl_mem),   (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  8, sizeof(cl_mem),   (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  9, sizeof(cl_mem),   (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 10, sizeof(cl_mem),   (void *)&dev_hash_header);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 11, sizeof(cl_mem),   (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup_local,   1, NULL, &global_work_size, &local_work_size, &hash_setup_local_event);

   ezcl_wait_for_events(1, &hash_setup_local_event);
   ezcl_event_release(hash_setup_local_event);

   if (DEBUG){
      vector<int> sizes_tmp(4);
      ezcl_enqueue_read_buffer(command_queue, dev_sizes, CL_TRUE,  0, 1*sizeof(cl_int4), &sizes_tmp[0], NULL);
      int iminsize_tmp = sizes_tmp[0];
      int imaxsize_tmp = sizes_tmp[1];
      int jminsize_tmp = sizes_tmp[2];
      int jmaxsize_tmp = sizes_tmp[3];
      fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize_tmp,imaxsize_tmp,jminsize_tmp,jmaxsize_tmp);
   }

   if (TIMING_LEVEL >= 2) {
      gpu_timers[MESH_TIMER_HASH_SETUP] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      cpu_timer_start(&tstart_lev2);
   }

#ifdef BOUNDS_CHECK
   {
      if (ezcl_get_device_mem_nelements(dev_nlft)  < (int)ncells ||
          ezcl_get_device_mem_nelements(dev_nrht)  < (int)ncells ||
          ezcl_get_device_mem_nelements(dev_nbot)  < (int)ncells ||
          ezcl_get_device_mem_nelements(dev_ntop)  < (int)ncells ||
          ezcl_get_device_mem_nelements(dev_i)     < (int)ncells ||
          ezcl_get_device_mem_nelements(dev_j)     < (int)ncells ||
          ezcl_get_device_mem_nelements(dev_level) < (int)ncells ) {
         printf("%d: Warning -- sizes for dev_neigh too small ncells %ld neigh %d %d %d %d %d %d %d\n",mype,ncells,ezcl_get_device_mem_nelements(dev_nlft),ezcl_get_device_mem_nelements(dev_nrht),ezcl_get_device_mem_nelements(dev_nbot),ezcl_get_device_mem_nelements(dev_ntop),ezcl_get_device_mem_nelements(dev_i),ezcl_get_device_mem_nelements(dev_j),ezcl_get_device_mem_nelements(dev_level));
      }
      vector<uchar_t> level_tmp(ncells);
      ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE, 0, ncells*sizeof(cl_uchar_t), &level_tmp[0], NULL);
      int iflag = 0;
      for (int ic=0; ic<ncells; ic++){
         if (levmx-level_tmp[ic] < 0 || levmx-level_tmp[ic] > levmx) {
            printf("%d: Warning level value bad ic %d level %d ncells %d\n",mype,ic,level_tmp[ic],ncells);
            iflag++;
         }
      }
      if (ezcl_get_device_mem_nelements(dev_levtable) < levmx+1) printf("%d Warning levtable too small levmx is %d devtable size is %d\n",mype,levmx,ezcl_get_device_mem_nelements(dev_levtable));
#ifdef HAVE_MPI
      if (iflag > 20) {fflush(stdout); L7_Terminate(); exit(0);}
#endif
   }
#endif

#ifdef BOUNDS_CHECK
   {
      int jmaxcalc = (jmax+1)*IPOW2(levmx);
      int imaxcalc = (imax+1)*IPOW2(levmx);
      vector<int> i_tmp(ncells);
      vector<int> j_tmp(ncells);
      vector<uchar_t> level_tmp(ncells);
      vector<int> hash_tmp(hashsize);
      ezcl_enqueue_read_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells*sizeof(cl_int), &i_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells*sizeof(cl_int), &j_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, ncells*sizeof(cl_uchar_t), &level_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_hash,  CL_TRUE,  0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);
      for (int ic=0; ic<(int)ncells; ic++){
         int ii  = i_tmp[ic];
         int jj  = j_tmp[ic];
         uchar_t lev = level_tmp[ic];
         int levmult = IPOW2(levmx-lev);
         int jjj=jj   *levmult-jminsize;
         int iii=max(  ii   *levmult-1, 0         )-iminsize;
         if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         jjj=jj   *levmult-jminsize;
         iii=min( (ii+1)*levmult,   imaxcalc-1)-iminsize;
         if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         jjj=max(  jj   *levmult-1, 0) -jminsize;
         iii=ii   *levmult   -iminsize;
         if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         jjj=min( (jj+1)*levmult,   jmaxcalc-1)-jminsize;
         iii=ii   *levmult   -iminsize;
         if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         int nlftval = hash_tmp[((      jj   *levmult               )-jminsize)*(imaxsize-iminsize)+((max(  ii   *levmult-1, 0         ))-iminsize)];
         int nrhtval = hash_tmp[((      jj   *levmult               )-jminsize)*(imaxsize-iminsize)+((min( (ii+1)*levmult,   imaxcalc-1))-iminsize)];
         int nbotval = hash_tmp[((max(  jj   *levmult-1, 0)         )-jminsize)*(imaxsize-iminsize)+((      ii   *levmult               )-iminsize)];
         int ntopval = hash_tmp[((min( (jj+1)*levmult,   jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((      ii   *levmult               )-iminsize)];

         if (nlftval == INT_MIN){
            jjj = jj*levmult-jminsize;
            iii = ii*levmult-iminsize;
            if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         }
         if (nrhtval == INT_MIN){
            jjj = jj*levmult-jminsize;
            iii = ii*levmult-iminsize;
            if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         }
         if (nbotval == INT_MIN) {
            iii = ii*levmult-iminsize;
            jjj = jj*levmult-jminsize;
            if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         }
         if (ntopval == INT_MIN) {
            iii = ii*levmult-iminsize;
            jjj = jj*levmult-jminsize;
            if (jjj < 0 || jjj >= (jmaxsize-jminsize) || iii < 0 || iii >= (imaxsize-iminsize) ) printf("%d: Warning at line %d iii %d jjj %d iminsize %d imaxsize %d jminsize %d jmaxsize %d\n",mype,__LINE__,iii,jjj,iminsize,imaxsize,jminsize,jmaxsize);
         }
      }
   }
#endif

   cl_event calc_neighbors_local_event;

      /*
                    const int   isize,       // 0
                    const int   levmx,       // 1
                    const int   imaxsize,    // 2
                    const int   jmaxsize,    // 3
                    const int   noffset,     // 4
           __global       int   *sizes,      // 5
           __global       int   *levtable,   // 6
           __global       uchar_t *level,      // 7
           __global       int   *i,          // 8
           __global       int   *j,          // 9
           __global       int   *nlft,       // 10
           __global       int   *nrht,       // 11
           __global       int   *nbot,       // 12
           __global       int   *ntop,       // 13
           __global const ulong *hash_heaer, // 14
           __global       int   *hash)       // 15
      */

   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 0,  sizeof(cl_int),   (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 1,  sizeof(cl_int),   (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 2,  sizeof(cl_int),   (void *)&imax);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 3,  sizeof(cl_int),   (void *)&jmax);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 4,  sizeof(cl_int),   (void *)&noffset);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 5,  sizeof(cl_mem),   (void *)&dev_sizes);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 6,  sizeof(cl_mem),   (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 7,  sizeof(cl_mem),   (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 8,  sizeof(cl_mem),   (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 9,  sizeof(cl_mem),   (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 10, sizeof(cl_mem),   (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 11, sizeof(cl_mem),   (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 12, sizeof(cl_mem),   (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 13, sizeof(cl_mem),   (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 14, sizeof(cl_mem),   (void *)&dev_hash_header);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 15, sizeof(cl_mem),   (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors_local,   1, NULL, &global_work_size, &local_work_size, &calc_neighbors_local_event);

   ezcl_wait_for_events(1, &calc_neighbors_local_event);
   ezcl_event_release(calc_neighbors_local_event);

   if (TIMING_LEVEL >= 2) {
      gpu_timers[MESH_TIMER_HASH_QUERY] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      cpu_timer_start(&tstart_lev2);
   }

   if (DEBUG) {
      print_dev_local();

      vector<int> hash_tmp(hashsize);
      ezcl_enqueue_read_buffer(command_queue, dev_hash, CL_FALSE, 0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);

      cl_mem dev_hash_header_check = gpu_get_hash_header();
      vector<ulong> hash_header_check(hash_header_size);
      ezcl_enqueue_read_buffer(command_queue, dev_hash_header_check, CL_TRUE, 0, hash_header_size*sizeof(cl_ulong), &hash_header_check[0], NULL);

      int   gpu_hash_method     = (int)hash_header_check[0];
      ulong gpu_hash_table_size =      hash_header_check[1];
      ulong gpu_AA              =      hash_header_check[2];
      ulong gpu_BB              =      hash_header_check[3];

      vector<int> nlft_tmp(ncells_ghost);
      vector<int> nrht_tmp(ncells_ghost);
      vector<int> nbot_tmp(ncells_ghost);
      vector<int> ntop_tmp(ncells_ghost);
      ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0], NULL);

      int jmaxglobal = (jmax+1)*IPOW2(levmx);
      int imaxglobal = (imax+1)*IPOW2(levmx);
      fprintf(fp,"\n                                    HASH 0 numbering\n");
      for (int jj = jmaxglobal-1; jj>=0; jj--){
         fprintf(fp,"%2d: %4d:",mype,jj);
         if (jj >= jminsize && jj < jmaxsize) {
            for (int ii = 0; ii<imaxglobal; ii++){
               if (ii >= iminsize && ii < imaxsize) {
                  fprintf(fp,"%5d",read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) );
               } else {
                  fprintf(fp,"     ");
               }
            }
         }
         fprintf(fp,"\n");
      }
      fprintf(fp,"%2d:      ",mype);
      for (int ii = 0; ii<imaxglobal; ii++){
         fprintf(fp,"%4d:",ii);
      }
      fprintf(fp,"\n");

      fprintf(fp,"\n                                    nlft numbering\n");
      for (int jj = jmaxglobal-1; jj>=0; jj--){
         fprintf(fp,"%2d: %4d:",mype,jj);
         if (jj >= jminsize && jj < jmaxsize) {
            for (int ii = 0; ii<imaxglobal; ii++){
               if (ii >= iminsize && ii < imaxsize) {
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) -noffset;
                  if (hashval >= 0 && hashval < (int)ncells) {
                     fprintf(fp,"%5d",nlft_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               } else {
                  fprintf(fp,"     ");
               }
            }
         }
         fprintf(fp,"\n");
      }
      fprintf(fp,"%2d:      ",mype);
      for (int ii = 0; ii<imaxglobal; ii++){
         fprintf(fp,"%4d:",ii);
      }
      fprintf(fp,"\n");
   
      fprintf(fp,"\n                                    nrht numbering\n");
      for (int jj = jmaxglobal-1; jj>=0; jj--){
         fprintf(fp,"%2d: %4d:",mype,jj);
         if (jj >= jminsize && jj < jmaxsize) {
            for (int ii = 0; ii<imaxglobal; ii++){
               if (ii >= iminsize && ii < imaxsize) {
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0])-noffset;
                  if (hashval >= 0 && hashval < (int)ncells) {
                     fprintf(fp,"%5d",nrht_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               } else {
                  fprintf(fp,"     ");
               }
            }
         }
         fprintf(fp,"\n");
      }
      fprintf(fp,"%2d:      ",mype);
      for (int ii = 0; ii<imaxglobal; ii++){
         fprintf(fp,"%4d:",ii);
      }
      fprintf(fp,"\n");

      fprintf(fp,"\n                                    nbot numbering\n");
      for (int jj = jmaxglobal-1; jj>=0; jj--){
         fprintf(fp,"%2d: %4d:",mype,jj);
         if (jj >= jminsize && jj < jmaxsize) {
            for (int ii = 0; ii<imaxglobal; ii++){
               if (ii >= iminsize && ii < imaxsize) {
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0])-noffset;
                  if (hashval >= 0 && hashval < (int)ncells) {
                     fprintf(fp,"%5d",nbot_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               } else {
                  fprintf(fp,"     ");
               }
            }
         }
         fprintf(fp,"\n");
      }
      fprintf(fp,"%2d:      ",mype);
      for (int ii = 0; ii<imaxglobal; ii++){
         fprintf(fp,"%4d:",ii);
      }
      fprintf(fp,"\n");

      fprintf(fp,"\n                                    ntop numbering\n");
      for (int jj = jmaxglobal-1; jj>=0; jj--){
         fprintf(fp,"%2d: %4d:",mype,jj);
         if (jj >= jminsize && jj < jmaxsize) {
            for (int ii = 0; ii<imaxglobal; ii++){
               if (ii >= iminsize && ii < imaxsize) {
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0])-noffset;
                  if (hashval >= 0 && hashval < (int)ncells) {
                     fprintf(fp,"%5d",ntop_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               } else {
                  fprintf(fp,"     ");
               }
            }
         }
         fprintf(fp,"\n");
      }
      fprintf(fp,"%2d:      ",mype);
      for (int ii = 0; ii<imaxglobal; ii++){
         fprintf(fp,"%4d:",ii);
      }
      fprintf(fp,"\n");
   }

#ifdef HAVE_MPI
   if (numpe > 1) {
         vector<int> iminsize_global(numpe);
         vector<int> imaxsize_global(numpe);
         vector<int> jminsize_global(numpe);
         vector<int> jmaxsize_global(numpe);
         vector<int> comm_partner(numpe,-1);

         MPI_Allgather(&iminsize, 1, MPI_INT, &iminsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);
         MPI_Allgather(&imaxsize, 1, MPI_INT, &imaxsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);
         MPI_Allgather(&jminsize, 1, MPI_INT, &jminsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);
         MPI_Allgather(&jmaxsize, 1, MPI_INT, &jmaxsize_global[0], 1, MPI_INT, MPI_COMM_WORLD);

         int num_comm_partners = 0; 
         for (int ip = 0; ip < numpe; ip++){
            if (ip == mype) continue;
            if (iminsize_global[ip] > imaxtile) continue;
            if (imaxsize_global[ip] < imintile) continue;
            if (jminsize_global[ip] > jmaxtile) continue;
            if (jmaxsize_global[ip] < jmintile) continue;
            comm_partner[num_comm_partners] = ip;
            num_comm_partners++;
            //if (DEBUG) fprintf(fp,"%d: overlap with processor %d bounding box is %d %d %d %d\n",mype,ip,iminsize_global[ip],imaxsize_global[ip],jminsize_global[ip],jmaxsize_global[ip]);
         }    

#ifdef BOUNDS_CHECK
      {
         vector<int> nlft_tmp(ncells_ghost);
         vector<int> nrht_tmp(ncells_ghost);
         vector<int> nbot_tmp(ncells_ghost);
         vector<int> ntop_tmp(ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells*sizeof(cl_int), &nlft_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells*sizeof(cl_int), &nrht_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells*sizeof(cl_int), &nbot_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells*sizeof(cl_int), &ntop_tmp[0], NULL);
         for (uint ic=0; ic<ncells; ic++){
            int nl = nlft_tmp[ic];
            if (nl != -1){
               nl -= noffset;
               if (nl<0 || nl>= ncells) printf("%d: Warning at line %d cell %d nlft %d\n",mype,__LINE__,ic,nl);
            }
            int nr = nrht_tmp[ic];
            if (nr != -1){
               nr -= noffset;
               if (nr<0 || nr>= ncells) printf("%d: Warning at line %d cell %d nrht %d\n",mype,__LINE__,ic,nr);
            }
            int nb = nbot_tmp[ic];
            if (nb != -1){
               nb -= noffset;
               if (nb<0 || nb>= ncells) printf("%d: Warning at line %d cell %d nbot %d\n",mype,__LINE__,ic,nb);
            }
            int nt = ntop_tmp[ic];
            if (nt != -1){
               nt -= noffset;
               if (nt<0 || nt>= ncells) printf("%d: Warning at line %d cell %d ntop %d\n",mype,__LINE__,ic,nt);
            }
         }
      }
#endif

      cl_mem dev_border_cell = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell1"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_calc_border_cells, 0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 1,  sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 2,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 3,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 4,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 5,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 6,  sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 7,  sizeof(cl_mem), (void *)&dev_border_cell);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_border_cells, 1, NULL, &global_work_size, &local_work_size, NULL); 

      cl_mem dev_border_cell_new = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell2"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      size_t one = 1;
      cl_mem dev_nbsize = ezcl_malloc(NULL, const_cast<char *>("dev_nbsize"), &one, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_ioffset = ezcl_malloc(NULL, const_cast<char *>("dev_ioffset"), &block_size, sizeof(cl_uint), CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_calc_border_cells2,  0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  1,  sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  2,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  3,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  4,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  5,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  6,  sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  7,  sizeof(cl_mem), (void *)&dev_border_cell);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  8,  sizeof(cl_mem), (void *)&dev_border_cell_new);
      ezcl_set_kernel_arg(kernel_calc_border_cells2,  9,  sizeof(cl_mem), (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 10,  sizeof(cl_mem), (void *)&dev_nbsize);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 11,  local_work_size*sizeof(cl_int), NULL);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_border_cells2, 1, NULL, &global_work_size, &local_work_size, NULL); 

      ezcl_device_memory_swap(&dev_border_cell, &dev_border_cell_new);
      ezcl_device_memory_delete(dev_border_cell_new);

      int group_size = (int)(global_work_size/local_work_size);

      ezcl_set_kernel_arg(kernel_finish_scan, 0,  sizeof(cl_int), (void *)&group_size);
      ezcl_set_kernel_arg(kernel_finish_scan, 1,  sizeof(cl_mem), (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_finish_scan, 2,  sizeof(cl_mem), (void *)&dev_nbsize);
      ezcl_set_kernel_arg(kernel_finish_scan, 3,  local_work_size*sizeof(cl_int), NULL);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_finish_scan, 1, NULL, &local_work_size, &local_work_size, NULL); 

      int nbsize_local;
      ezcl_enqueue_read_buffer(command_queue, dev_nbsize, CL_TRUE,  0, 1*sizeof(cl_int), &nbsize_local, NULL);
      ezcl_device_memory_delete(dev_nbsize);

      //printf("%d: border cell size is %d global is %ld\n",mype,nbsize_local,nbsize_global);

      vector<int> border_cell_num(nbsize_local);
      vector<int> border_cell_i(nbsize_local);
      vector<int> border_cell_j(nbsize_local);
      vector<int> border_cell_level(nbsize_local);
    
      // allocate new border memory
      size_t nbsize_long = nbsize_local;
      cl_mem dev_border_cell_i     = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_i"),     &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_j     = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_j"),     &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_level = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_level"), &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_num   = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_num"),   &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_get_border_data,  0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_get_border_data,  1,  sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_get_border_data,  2,  sizeof(cl_mem), (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_get_border_data,  3,  sizeof(cl_mem), (void *)&dev_border_cell);
      ezcl_set_kernel_arg(kernel_get_border_data,  4,  sizeof(cl_mem), (void *)&dev_i);
      ezcl_set_kernel_arg(kernel_get_border_data,  5,  sizeof(cl_mem), (void *)&dev_j);
      ezcl_set_kernel_arg(kernel_get_border_data,  6,  sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_get_border_data,  7,  sizeof(cl_mem), (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_get_border_data,  8,  sizeof(cl_mem), (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_get_border_data,  9,  sizeof(cl_mem), (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_get_border_data, 10,  sizeof(cl_mem), (void *)&dev_border_cell_num);
      ezcl_set_kernel_arg(kernel_get_border_data, 11,  local_work_size*sizeof(cl_uint), NULL);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_get_border_data, 1, NULL, &global_work_size, &local_work_size, NULL); 

      ezcl_device_memory_delete(dev_ioffset);
      ezcl_device_memory_delete(dev_border_cell);

      // read gpu border cell data
      ezcl_enqueue_read_buffer(command_queue, dev_border_cell_i,     CL_FALSE, 0, nbsize_local*sizeof(cl_int), &border_cell_i[0],     NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_border_cell_j,     CL_FALSE, 0, nbsize_local*sizeof(cl_int), &border_cell_j[0],     NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_border_cell_level, CL_FALSE, 0, nbsize_local*sizeof(cl_int), &border_cell_level[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_border_cell_num,   CL_TRUE,  0, nbsize_local*sizeof(cl_int), &border_cell_num[0],   NULL);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_FIND_BOUNDARY] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      // Allocate push database

      int **send_database = (int**)malloc(num_comm_partners*sizeof(int *));
      for (int ip = 0; ip < num_comm_partners; ip++){
         send_database[ip] = (int *)malloc(nbsize_local*sizeof(int));
      }

      // Compute the overlap between processor bounding boxes and set up push database

      vector<int> send_buffer_count(num_comm_partners);
      for (int ip = 0; ip < num_comm_partners; ip++){
         int icount = 0;
         for (int ib = 0; ib <nbsize_local; ib++){
            int lev = border_cell_level[ib];
            int levmult = IPOW2(levmx-lev);
            if (border_cell_i[ib]*levmult >= iminsize_global[comm_partner[ip]] && 
                border_cell_i[ib]*levmult <= imaxsize_global[comm_partner[ip]] && 
                border_cell_j[ib]*levmult >= jminsize_global[comm_partner[ip]] && 
                border_cell_j[ib]*levmult <= jmaxsize_global[comm_partner[ip]] ) {
               send_database[ip][icount] = ib;
               icount++;
            }
         }
         send_buffer_count[ip]=icount;
      }

      // Initialize L7_Push_Setup with num_comm_partners, comm_partner, send_database and 
      // send_buffer_count. L7_Push_Setup will copy data and determine recv_buffer_counts.
      // It will return receive_count_total for use in allocations

      int receive_count_total;
      int i_push_handle = 0;
      L7_Push_Setup(num_comm_partners, &comm_partner[0], &send_buffer_count[0],
                    send_database, &receive_count_total, &i_push_handle);

      if (DEBUG) {
         fprintf(fp,"DEBUG num_comm_partners %d\n",num_comm_partners);
         for (int ip = 0; ip < num_comm_partners; ip++){
            fprintf(fp,"DEBUG comm partner is %d data count is %d\n",comm_partner[ip],send_buffer_count[ip]);
            for (int ic = 0; ic < send_buffer_count[ip]; ic++){
               int ib = send_database[ip][ic];
               fprintf(fp,"DEBUG \t index %d cell number %d i %d j %d level %d\n",ib,border_cell_num[ib],
                  border_cell_i[ib],border_cell_j[ib],border_cell_level[ib]);
            }
         }
      }

      // Can now free the send database. Other arrays are vectors and will automatically 
      // deallocate

      for (int ip = 0; ip < num_comm_partners; ip++){
         free(send_database[ip]);
      }
      free(send_database);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_PUSH_SETUP] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }
      // Push the data needed to the adjacent processors

      int *border_cell_num_local = (int *)malloc(receive_count_total*sizeof(int));
      int *border_cell_i_local = (int *)malloc(receive_count_total*sizeof(int));
      int *border_cell_j_local = (int *)malloc(receive_count_total*sizeof(int));
      int *border_cell_level_local = (int *)malloc(receive_count_total*sizeof(int));
      L7_Push_Update(&border_cell_num[0],   border_cell_num_local,   i_push_handle);
      L7_Push_Update(&border_cell_i[0],     border_cell_i_local,     i_push_handle);
      L7_Push_Update(&border_cell_j[0],     border_cell_j_local,     i_push_handle);
      L7_Push_Update(&border_cell_level[0], border_cell_level_local, i_push_handle);

      L7_Push_Free(&i_push_handle);

      ezcl_device_memory_delete(dev_border_cell_i);
      ezcl_device_memory_delete(dev_border_cell_j);
      ezcl_device_memory_delete(dev_border_cell_level);
      ezcl_device_memory_delete(dev_border_cell_num);

      nbsize_local = receive_count_total;

      if (DEBUG) {
         for (int ic = 0; ic < nbsize_local; ic++) {
            fprintf(fp,"%d: Local Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num_local[ic],
               border_cell_i_local[ic],border_cell_j_local[ic],border_cell_level_local[ic]);
         }
      }

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_PUSH_BOUNDARY] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      nbsize_long = nbsize_local;

      dev_border_cell_num        = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_num"),        &nbsize_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      dev_border_cell_i          = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_i"),          &nbsize_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      dev_border_cell_j          = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_j"),          &nbsize_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      dev_border_cell_level      = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_level"),      &nbsize_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_needed     = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_needed"),     &nbsize_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_needed_out = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_needed_out"), &nbsize_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      ezcl_enqueue_write_buffer(command_queue, dev_border_cell_num,    CL_FALSE, 0, nbsize_local*sizeof(cl_int), &border_cell_num_local[0], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_border_cell_i,      CL_FALSE, 0, nbsize_local*sizeof(cl_int), &border_cell_i_local[0],   NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_border_cell_j,      CL_FALSE, 0, nbsize_local*sizeof(cl_int), &border_cell_j_local[0],   NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_border_cell_level,  CL_TRUE,  0, nbsize_local*sizeof(cl_int), &border_cell_level_local[0],   NULL);

      //ezcl_enqueue_write_buffer(command_queue, dev_border_cell_needed, CL_TRUE,  0, nbsize_local*sizeof(cl_int), &border_cell_needed_local[0],   NULL);

      free(border_cell_i_local);
      free(border_cell_j_local);
      free(border_cell_level_local);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_LOCAL_LIST] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      if (DEBUG) {
         vector<int> hash_tmp(hashsize);
         ezcl_enqueue_read_buffer(command_queue, dev_hash, CL_TRUE,  0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);

         cl_mem dev_hash_header_check = gpu_get_hash_header();
         vector<ulong> hash_header_check(hash_header_size);
         ezcl_enqueue_read_buffer(command_queue, dev_hash_header_check, CL_TRUE, 0, hash_header_size*sizeof(cl_ulong), &hash_header_check[0], NULL);

         int   gpu_hash_method     = (int)hash_header_check[0];
         ulong gpu_hash_table_size =      hash_header_check[1];
         ulong gpu_AA              =      hash_header_check[2];
         ulong gpu_BB              =      hash_header_check[3];

         int jmaxglobal = (jmax+1)*IPOW2(levmx);
         int imaxglobal = (imax+1)*IPOW2(levmx);
         fprintf(fp,"\n                                    HASH numbering before layer 1\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     fprintf(fp,"%5d",read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) );
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
      }

      size_t nb_local_work_size = 128;
      size_t nb_global_work_size = ((nbsize_local + nb_local_work_size - 1) /nb_local_work_size) * nb_local_work_size;

      ezcl_set_kernel_arg(kernel_calc_layer1,  0,  sizeof(cl_int),   (void *)&nbsize_local);
      ezcl_set_kernel_arg(kernel_calc_layer1,  1,  sizeof(cl_int),   (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_layer1,  2,  sizeof(cl_int),   (void *)&levmx);
      ezcl_set_kernel_arg(kernel_calc_layer1,  3,  sizeof(cl_int),   (void *)&imax);
      ezcl_set_kernel_arg(kernel_calc_layer1,  4,  sizeof(cl_int),   (void *)&jmax);
      ezcl_set_kernel_arg(kernel_calc_layer1,  5,  sizeof(cl_int),   (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_layer1,  6,  sizeof(cl_mem),   (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_calc_layer1,  7,  sizeof(cl_mem),   (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_calc_layer1,  8,  sizeof(cl_mem),   (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_calc_layer1,  9,  sizeof(cl_mem),   (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_calc_layer1, 10,  sizeof(cl_mem),   (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_calc_layer1, 11,  sizeof(cl_mem),   (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_calc_layer1, 12,  sizeof(cl_mem),   (void *)&dev_border_cell_needed);
      ezcl_set_kernel_arg(kernel_calc_layer1, 13,  sizeof(cl_mem),   (void *)&dev_hash_header);
      ezcl_set_kernel_arg(kernel_calc_layer1, 14,  sizeof(cl_mem),   (void *)&dev_hash);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_layer1, 1, NULL, &nb_global_work_size, &nb_local_work_size, NULL); 

      if (DEBUG){
         vector<int> border_cell_needed_local(nbsize_local);

         ezcl_enqueue_read_buffer(command_queue, dev_border_cell_needed, CL_TRUE,  0, nbsize_local*sizeof(cl_int), &border_cell_needed_local[0],   NULL);

         for(int ic=0; ic<nbsize_local; ic++){
            if (border_cell_needed_local[ic] == 0) continue;
            fprintf(fp,"%d: First set of needed cells ic %3d cell %3d type %3d\n",mype,ic,border_cell_num_local[ic],border_cell_needed_local[ic]);
         }
      }

      cl_event calc_layer1_sethash_event;

      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  0,  sizeof(cl_int),   (void *)&nbsize_local);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  1,  sizeof(cl_int),   (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  2,  sizeof(cl_int),   (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  3,  sizeof(cl_int),   (void *)&levmx);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  4,  sizeof(cl_mem),   (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  5,  sizeof(cl_mem),   (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  6,  sizeof(cl_mem),   (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  7,  sizeof(cl_mem),   (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  8,  sizeof(cl_mem),   (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash,  9,  sizeof(cl_mem),   (void *)&dev_border_cell_needed);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash, 10,  sizeof(cl_mem),   (void *)&dev_hash_header);
      ezcl_set_kernel_arg(kernel_calc_layer1_sethash, 11,  sizeof(cl_mem),   (void *)&dev_hash);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_layer1_sethash, 1, NULL, &nb_global_work_size, &nb_local_work_size, &calc_layer1_sethash_event); 

      ezcl_wait_for_events(1, &calc_layer1_sethash_event);
      ezcl_event_release(calc_layer1_sethash_event);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_LAYER1] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      if (DEBUG) {
         print_dev_local();

         vector<int> hash_tmp(hashsize);
         ezcl_enqueue_read_buffer(command_queue, dev_hash, CL_TRUE,  0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);

         cl_mem dev_hash_header_check = gpu_get_hash_header();
         vector<ulong> hash_header_check(hash_header_size);
         ezcl_enqueue_read_buffer(command_queue, dev_hash_header_check, CL_TRUE, 0, hash_header_size*sizeof(cl_ulong), &hash_header_check[0], NULL);

         int   gpu_hash_method     = (int)hash_header_check[0];
         ulong gpu_hash_table_size =      hash_header_check[1];
         ulong gpu_AA              =      hash_header_check[2];
         ulong gpu_BB              =      hash_header_check[3];

         int jmaxglobal = (jmax+1)*IPOW2(levmx);
         int imaxglobal = (imax+1)*IPOW2(levmx);
         fprintf(fp,"\n                                    HASH numbering for 1 layer\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     fprintf(fp,"%5d",read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) );
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
      }

      group_size = (int)(nb_global_work_size/nb_local_work_size);

      cl_mem dev_nbpacked = ezcl_malloc(NULL, const_cast<char *>("dev_nbpacked"), &one, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      size_t group_size_long = group_size;
      dev_ioffset = ezcl_malloc(NULL, const_cast<char *>("dev_ioffset"), &group_size_long, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_calc_layer2,  0,  sizeof(cl_int),   (void *)&nbsize_local);
      ezcl_set_kernel_arg(kernel_calc_layer2,  1,  sizeof(cl_int),   (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_layer2,  2,  sizeof(cl_int),   (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_layer2,  3,  sizeof(cl_int),   (void *)&levmx);
      ezcl_set_kernel_arg(kernel_calc_layer2,  4,  sizeof(cl_int),   (void *)&imax);
      ezcl_set_kernel_arg(kernel_calc_layer2,  5,  sizeof(cl_int),   (void *)&jmax);
      ezcl_set_kernel_arg(kernel_calc_layer2,  6,  sizeof(cl_mem),   (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_calc_layer2,  7,  sizeof(cl_mem),   (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_calc_layer2,  8,  sizeof(cl_mem),   (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_calc_layer2,  9,  sizeof(cl_mem),   (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_calc_layer2, 10,  sizeof(cl_mem),   (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_calc_layer2, 11,  sizeof(cl_mem),   (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_calc_layer2, 12,  sizeof(cl_mem),   (void *)&dev_border_cell_needed);
      ezcl_set_kernel_arg(kernel_calc_layer2, 13,  sizeof(cl_mem),   (void *)&dev_border_cell_needed_out);
      ezcl_set_kernel_arg(kernel_calc_layer2, 14,  sizeof(cl_mem),   (void *)&dev_hash_header);
      ezcl_set_kernel_arg(kernel_calc_layer2, 15,  sizeof(cl_mem),   (void *)&dev_hash);
      ezcl_set_kernel_arg(kernel_calc_layer2, 16,  sizeof(cl_mem),   (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_calc_layer2, 17,  sizeof(cl_mem),   (void *)&dev_nbpacked);
      ezcl_set_kernel_arg(kernel_calc_layer2, 18,  nb_local_work_size*sizeof(cl_mem), NULL);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_layer2, 1, NULL, &nb_global_work_size, &nb_local_work_size, NULL); 

      if (DEBUG){
         vector<int> border_cell_needed_local(nbsize_local);

         ezcl_enqueue_read_buffer(command_queue, dev_border_cell_needed_out, CL_TRUE,  0, nbsize_local*sizeof(cl_int), &border_cell_needed_local[0],   NULL);
         for(int ic=0; ic<nbsize_local; ic++){
            if (border_cell_needed_local[ic] <= 0) continue;
            if (border_cell_needed_local[ic] <  0x0016) fprintf(fp,"%d: First  set of needed cells ic %3d cell %3d type %3d\n",mype,ic,border_cell_num_local[ic],border_cell_needed_local[ic]);
            if (border_cell_needed_local[ic] >= 0x0016) fprintf(fp,"%d: Second set of needed cells ic %3d cell %3d type %3d\n",mype,ic,border_cell_num_local[ic],border_cell_needed_local[ic]);
         }
      }

      free(border_cell_num_local);

      ezcl_device_memory_delete(dev_border_cell_needed);

      ezcl_set_kernel_arg(kernel_finish_scan, 0,  sizeof(cl_int), (void *)&group_size);
      ezcl_set_kernel_arg(kernel_finish_scan, 1,  sizeof(cl_mem), (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_finish_scan, 2,  sizeof(cl_mem), (void *)&dev_nbpacked);
      ezcl_set_kernel_arg(kernel_finish_scan, 3,  nb_local_work_size*sizeof(cl_int), NULL);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_finish_scan, 1, NULL, &nb_local_work_size, &nb_local_work_size, NULL); 

      int nbpacked;
      ezcl_enqueue_read_buffer(command_queue, dev_nbpacked, CL_TRUE,  0, 1*sizeof(cl_int), &nbpacked, NULL);
      ezcl_device_memory_delete(dev_nbpacked);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_LAYER2] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      nbsize_long = nbsize_local;
      cl_mem dev_border_cell_i_new     = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_i_new"),     &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_j_new     = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_j_new"),     &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell_level_new = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell_level_new"), &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_indices_needed    = ezcl_malloc(NULL, const_cast<char *>("dev_indices_needed"),    &nbsize_long, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

      cl_event get_border_data2_event;

      ezcl_set_kernel_arg(kernel_get_border_data2,  0,  sizeof(cl_int), (void *)&nbsize_local);
      ezcl_set_kernel_arg(kernel_get_border_data2,  1,  sizeof(cl_mem), (void *)&dev_ioffset);
      ezcl_set_kernel_arg(kernel_get_border_data2,  2,  sizeof(cl_mem), (void *)&dev_border_cell_needed_out);
      ezcl_set_kernel_arg(kernel_get_border_data2,  3,  sizeof(cl_mem), (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_get_border_data2,  4,  sizeof(cl_mem), (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_get_border_data2,  5,  sizeof(cl_mem), (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_get_border_data2,  6,  sizeof(cl_mem), (void *)&dev_border_cell_num);
      ezcl_set_kernel_arg(kernel_get_border_data2,  7,  sizeof(cl_mem), (void *)&dev_border_cell_i_new);
      ezcl_set_kernel_arg(kernel_get_border_data2,  8,  sizeof(cl_mem), (void *)&dev_border_cell_j_new);
      ezcl_set_kernel_arg(kernel_get_border_data2,  9,  sizeof(cl_mem), (void *)&dev_border_cell_level_new);
      ezcl_set_kernel_arg(kernel_get_border_data2, 10,  sizeof(cl_mem), (void *)&dev_indices_needed);
      ezcl_set_kernel_arg(kernel_get_border_data2, 11,  local_work_size*sizeof(cl_uint), NULL);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_get_border_data2, 1, NULL, &nb_global_work_size, &nb_local_work_size, &get_border_data2_event);

      ezcl_device_memory_delete(dev_border_cell_num);

      ezcl_device_memory_swap(&dev_border_cell_i,     &dev_border_cell_i_new);
      ezcl_device_memory_swap(&dev_border_cell_j,     &dev_border_cell_j_new);
      ezcl_device_memory_swap(&dev_border_cell_level, &dev_border_cell_level_new);

      size_t nbp_local_work_size = 128;
      size_t nbp_global_work_size = ((nbpacked + nbp_local_work_size - 1) /nbp_local_work_size) * nbp_local_work_size;

      cl_event calc_layer2_sethash_event;

      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  0,  sizeof(cl_int),   (void *)&nbpacked);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  1,  sizeof(cl_int),   (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  2,  sizeof(cl_int),   (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  3,  sizeof(cl_int),   (void *)&levmx);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  4,  sizeof(cl_int),   (void *)&imax);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  5,  sizeof(cl_int),   (void *)&jmax);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  6,  sizeof(cl_mem),   (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  7,  sizeof(cl_mem),   (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  8,  sizeof(cl_mem),   (void *)&dev_levibeg);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash,  9,  sizeof(cl_mem),   (void *)&dev_leviend);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 10,  sizeof(cl_mem),   (void *)&dev_levjbeg);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 11,  sizeof(cl_mem),   (void *)&dev_levjend);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 12,  sizeof(cl_mem),   (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 13,  sizeof(cl_mem),   (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 14,  sizeof(cl_mem),   (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 15,  sizeof(cl_mem),   (void *)&dev_indices_needed);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 16,  sizeof(cl_mem),   (void *)&dev_border_cell_needed_out);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 17,  sizeof(cl_mem),   (void *)&dev_hash_header);
      ezcl_set_kernel_arg(kernel_calc_layer2_sethash, 18,  sizeof(cl_mem),   (void *)&dev_hash);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_layer2_sethash, 1, NULL, &nbp_global_work_size, &nbp_local_work_size, &calc_layer2_sethash_event); 

      ezcl_wait_for_events(1, &calc_layer2_sethash_event);
      ezcl_event_release(calc_layer2_sethash_event);

      ezcl_device_memory_delete(dev_ioffset);

      ezcl_wait_for_events(1, &get_border_data2_event);
      ezcl_event_release(get_border_data2_event);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_LAYER_LIST] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      vector<int> indices_needed(nbpacked);

      // read gpu border cell data 
      ezcl_enqueue_read_buffer(command_queue, dev_indices_needed,    CL_TRUE,  0, nbpacked*sizeof(cl_int), &indices_needed[0],    NULL);

      ezcl_device_memory_delete(dev_border_cell_i_new);
      ezcl_device_memory_delete(dev_border_cell_j_new);
      ezcl_device_memory_delete(dev_border_cell_level_new);

      if (DEBUG) {
         print_dev_local();

         vector<int> hash_tmp(hashsize);
         ezcl_enqueue_read_buffer(command_queue, dev_hash, CL_TRUE,  0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);

         cl_mem dev_hash_header_check = gpu_get_hash_header();
         vector<ulong> hash_header_check(hash_header_size);
         ezcl_enqueue_read_buffer(command_queue, dev_hash_header_check, CL_TRUE, 0, hash_header_size*sizeof(cl_ulong), &hash_header_check[0], NULL);

         int   gpu_hash_method     = (int)hash_header_check[0];
         ulong gpu_hash_table_size =      hash_header_check[1];
         ulong gpu_AA              =      hash_header_check[2];
         ulong gpu_BB              =      hash_header_check[3];

         int jmaxglobal = (jmax+1)*IPOW2(levmx);
         int imaxglobal = (imax+1)*IPOW2(levmx);
         fprintf(fp,"\n                                    HASH numbering for 2 layer\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     fprintf(fp,"%5d",read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) );
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
         fflush(fp);
      }

      ezcl_device_memory_delete(dev_border_cell_needed_out);

      int nghost = nbpacked;
      ncells_ghost = ncells + nghost;

      //if (mype == 1) printf("%d: DEBUG before expanding memory ncells %ld ncells_ghost %ld capacity %ld\n",mype,ncells,ncells_ghost,ezcl_get_device_mem_capacity(dev_i));
      if (ezcl_get_device_mem_capacity(dev_celltype) < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_i)        < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_j)        < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_level)    < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_nlft)     < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_nrht)     < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_nbot)     < ncells_ghost ||
          ezcl_get_device_mem_capacity(dev_ntop)     < ncells_ghost ) {

         //if (mype == 0) printf("%d: DEBUG expanding memory ncells %ld ncells_ghost %ld capacity %ld\n",mype,ncells,ncells_ghost,ezcl_get_device_mem_capacity(dev_i));
         //printf("%d: DEBUG expanding memory ncells %ld ncells_ghost %ld capacity %ld\n",mype,ncells,ncells_ghost,ezcl_get_device_mem_capacity(dev_i));
         mem_factor = (float)(ncells_ghost/ncells);
         cl_mem dev_celltype_old = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_old"), &ncells_ghost, sizeof(cl_char_t), CL_MEM_READ_WRITE, 0);
         cl_mem dev_i_old        = ezcl_malloc(NULL, const_cast<char *>("dev_i_old"),        &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_j_old        = ezcl_malloc(NULL, const_cast<char *>("dev_j_old"),        &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_level_old    = ezcl_malloc(NULL, const_cast<char *>("dev_level_old"),    &ncells_ghost, sizeof(cl_uchar_t), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nlft_old     = ezcl_malloc(NULL, const_cast<char *>("dev_nlft_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nrht_old     = ezcl_malloc(NULL, const_cast<char *>("dev_nrht_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nbot_old     = ezcl_malloc(NULL, const_cast<char *>("dev_nbot_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_ntop_old     = ezcl_malloc(NULL, const_cast<char *>("dev_ntop_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         ezcl_device_memory_swap(&dev_celltype_old, &dev_celltype);
         ezcl_device_memory_swap(&dev_i_old,        &dev_i       );
         ezcl_device_memory_swap(&dev_j_old,        &dev_j       );
         ezcl_device_memory_swap(&dev_level_old,    &dev_level   );
         ezcl_device_memory_swap(&dev_nlft_old,     &dev_nlft    );
         ezcl_device_memory_swap(&dev_nrht_old,     &dev_nrht    );
         ezcl_device_memory_swap(&dev_nbot_old,     &dev_nbot    );
         ezcl_device_memory_swap(&dev_ntop_old,     &dev_ntop    );

         cl_event copy_mesh_data_event;

         ezcl_set_kernel_arg(kernel_copy_mesh_data, 0,  sizeof(cl_int), (void *)&ncells);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 1,  sizeof(cl_mem), (void *)&dev_celltype_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 2,  sizeof(cl_mem), (void *)&dev_celltype);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 3,  sizeof(cl_mem), (void *)&dev_i_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 4,  sizeof(cl_mem), (void *)&dev_i);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 5,  sizeof(cl_mem), (void *)&dev_j_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 6,  sizeof(cl_mem), (void *)&dev_j);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 7,  sizeof(cl_mem), (void *)&dev_level_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 8,  sizeof(cl_mem), (void *)&dev_level);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 9,  sizeof(cl_mem), (void *)&dev_nlft_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 10, sizeof(cl_mem), (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 11, sizeof(cl_mem), (void *)&dev_nrht_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 12, sizeof(cl_mem), (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 13, sizeof(cl_mem), (void *)&dev_nbot_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 14, sizeof(cl_mem), (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 15, sizeof(cl_mem), (void *)&dev_ntop_old);
         ezcl_set_kernel_arg(kernel_copy_mesh_data, 16, sizeof(cl_mem), (void *)&dev_ntop);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_mesh_data,   1, NULL, &global_work_size, &local_work_size, &copy_mesh_data_event);

         ezcl_device_memory_delete(dev_celltype_old);
         ezcl_device_memory_delete(dev_i_old);
         ezcl_device_memory_delete(dev_j_old);
         ezcl_device_memory_delete(dev_level_old);
         ezcl_device_memory_delete(dev_nlft_old);
         ezcl_device_memory_delete(dev_nrht_old);
         ezcl_device_memory_delete(dev_nbot_old);
         ezcl_device_memory_delete(dev_ntop_old);

         ezcl_wait_for_events(1, &copy_mesh_data_event);
         ezcl_event_release(copy_mesh_data_event);
      }

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_COPY_MESH_DATA] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      nb_global_work_size = ((nbpacked + nb_local_work_size - 1) /nb_local_work_size) * nb_local_work_size;

#ifdef BOUNDS_CHECK
      if (ezcl_get_device_mem_nelements(dev_i) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_j) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_level) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_celltype) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_nlft) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_nrht) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_nbot) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_ntop) < (int)ncells_ghost ){
             printf("DEBUG size issue at %d\n",__LINE__);
      }
      if (ezcl_get_device_mem_nelements(dev_border_cell_i) < nbpacked || 
          ezcl_get_device_mem_nelements(dev_border_cell_j) < nbpacked || 
          ezcl_get_device_mem_nelements(dev_border_cell_level) < nbpacked ){
             printf("DEBUG size issue at %d\n",__LINE__);
      }
#endif
 
      cl_event fill_mesh_ghost_event;

      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  0,  sizeof(cl_int), (void *)&nbpacked);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  1,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  2,  sizeof(cl_mem), (void *)&dev_levibeg);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  3,  sizeof(cl_mem), (void *)&dev_leviend);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  4,  sizeof(cl_mem), (void *)&dev_levjbeg);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  5,  sizeof(cl_mem), (void *)&dev_levjend);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  6,  sizeof(cl_mem), (void *)&dev_border_cell_i);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  7,  sizeof(cl_mem), (void *)&dev_border_cell_j);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  8,  sizeof(cl_mem), (void *)&dev_border_cell_level);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost,  9,  sizeof(cl_mem), (void *)&dev_i);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 10,  sizeof(cl_mem), (void *)&dev_j);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 11,  sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 12,  sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 13,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 14,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 15,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_fill_mesh_ghost, 16,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_fill_mesh_ghost, 1, NULL, &nb_global_work_size, &nb_local_work_size, &fill_mesh_ghost_event); 

      ezcl_wait_for_events(1, &fill_mesh_ghost_event);
      ezcl_event_release(fill_mesh_ghost_event);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_FILL_MESH_GHOST] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      if (DEBUG){
         fprintf(fp,"After copying i,j, level to ghost cells\n");
         print_dev_local();
      }

      ezcl_device_memory_delete(dev_border_cell_i);
      ezcl_device_memory_delete(dev_border_cell_j);
      ezcl_device_memory_delete(dev_border_cell_level);

      size_t ghost_local_work_size = 128;
      size_t ghost_global_work_size = ((ncells_ghost + ghost_local_work_size - 1) /ghost_local_work_size) * ghost_local_work_size;

      cl_event fill_neighbor_ghost_event;

      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  0,  sizeof(cl_int),   (void *)&ncells_ghost);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  1,  sizeof(cl_int),   (void *)&levmx);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  2,  sizeof(cl_int),   (void *)&imax);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  3,  sizeof(cl_int),   (void *)&jmax);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  4,  sizeof(cl_mem),   (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  5,  sizeof(cl_mem),   (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  6,  sizeof(cl_mem),   (void *)&dev_i);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  7,  sizeof(cl_mem),   (void *)&dev_j);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  8,  sizeof(cl_mem),   (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost,  9,  sizeof(cl_mem),   (void *)&dev_hash_header);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost, 10,  sizeof(cl_mem),   (void *)&dev_hash);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost, 11,  sizeof(cl_mem),   (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost, 12,  sizeof(cl_mem),   (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost, 13,  sizeof(cl_mem),   (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_fill_neighbor_ghost, 14,  sizeof(cl_mem),   (void *)&dev_ntop);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_fill_neighbor_ghost, 1, NULL, &ghost_global_work_size, &ghost_local_work_size, &fill_neighbor_ghost_event); 

      ezcl_wait_for_events(1, &fill_neighbor_ghost_event);
      ezcl_event_release(fill_neighbor_ghost_event);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_FILL_NEIGH_GHOST] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      if (DEBUG){
         fprintf(fp,"After setting neighbors through ghost cells\n");
         print_dev_local();
      }

#ifdef BOUNDS_CHECK
      if (ezcl_get_device_mem_nelements(dev_nlft) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_nrht) < (int)ncells_ghost ||
          ezcl_get_device_mem_nelements(dev_nbot) < (int)ncells_ghost ||
          ezcl_get_device_mem_nelements(dev_ntop) < (int)ncells_ghost ){
         printf("%d: Warning sizes for set_corner_neighbor not right ncells ghost %d nlft size %d\n",mype,ncells_ghost,ezcl_get_device_mem_nelements(dev_nlft));
      }
#endif

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_SET_CORNER_NEIGH] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      if (DEBUG){
         fprintf(fp,"After setting corner neighbors\n");
         print_dev_local();
      }

#ifdef BOUNDS_CHECK
      if (ezcl_get_device_mem_nelements(dev_nlft) < (int)ncells_ghost || 
          ezcl_get_device_mem_nelements(dev_nrht) < (int)ncells_ghost ||
          ezcl_get_device_mem_nelements(dev_nbot) < (int)ncells_ghost ||
          ezcl_get_device_mem_nelements(dev_ntop) < (int)ncells_ghost ){
         printf("%d: Warning sizes for adjust neighbors not right ncells ghost %d nlft size %d\n",mype,ncells_ghost,ezcl_get_device_mem_nelements(dev_nlft));
      }
      if (ezcl_get_device_mem_nelements(dev_indices_needed) < (int)(ncells_ghost-ncells) ){
         printf("%d: Warning indices size wrong nghost %d size indices_needed\n",mype,ncells_ghost-ncells,ezcl_get_device_mem_nelements(dev_indices_needed));
      }
#endif

      cl_event adjust_neighbors_local_event;

      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  0,  sizeof(cl_int), (void *)&ncells_ghost);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  1,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  2,  sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  3,  sizeof(cl_mem), (void *)&dev_indices_needed);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  4,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  5,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  6,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_adjust_neighbors_local,  7,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_adjust_neighbors_local, 1, NULL, &ghost_global_work_size, &ghost_local_work_size, &adjust_neighbors_local_event); 

      ezcl_device_memory_delete(dev_indices_needed);

      if (DEBUG){
         fprintf(fp,"After adjusting neighbors to local indices\n");
         print_dev_local();
      }

      ezcl_wait_for_events(1, &adjust_neighbors_local_event);
      ezcl_event_release(adjust_neighbors_local_event);

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_NEIGH_ADJUST] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      offtile_ratio_local = (offtile_ratio_local*(double)offtile_local_count) + ((double)nghost / (double)ncells);
      offtile_local_count++;
      offtile_ratio_local /= offtile_local_count;

      if (cell_handle) L7_Free(&cell_handle);
      cell_handle=0;

      if (DEBUG){
         fprintf(fp,"%d: SETUP ncells %ld noffset %d nghost %d\n",mype,ncells,noffset,nghost);
         for (int ic=0; ic<nghost; ic++){
            fprintf(fp,"%d: indices needed ic %d index %d\n",mype,ic,indices_needed[ic]);
         }
      }

      L7_Dev_Setup(0, noffset, ncells, &indices_needed[0], nghost, &cell_handle);

#ifdef BOUNDS_CHECK
      {
         vector<int> nlft_tmp(ncells_ghost);
         vector<int> nrht_tmp(ncells_ghost);
         vector<int> nbot_tmp(ncells_ghost);
         vector<int> ntop_tmp(ncells_ghost);
         vector<uchar_t> level_tmp(ncells_ghost);
         vector<real_t> H_tmp(ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, ncells_ghost*sizeof(cl_uchar_t), &level_tmp[0], NULL);
         for (uint ic=0; ic<ncells; ic++){
            int nl = nlft_tmp[ic];
            if (nl<0 || nl>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nlft %d\n",mype,__LINE__,ic,nl);
            if (level_tmp[nl] > level_tmp[ic]){
               int ntl = ntop_tmp[nl];
               if (ntl<0 || ntl>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d global %d nlft %d ntop of nlft %d\n",mype,__LINE__,ic,ic+noffset,nl,ntl);
            }
            int nr = nrht_tmp[ic];
            if (nr<0 || nr>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nrht %d\n",mype,__LINE__,ic,nr);
            if (level_tmp[nr] > level_tmp[ic]){
               int ntr = ntop_tmp[nr];
               if (ntr<0 || ntr>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d ntop of nrht %d\n",mype,__LINE__,ic,ntr);
            }
            int nb = nbot_tmp[ic];
            if (nb<0 || nb>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nbot %d\n",mype,__LINE__,ic,nb);
            if (level_tmp[nb] > level_tmp[ic]){
               int nrb = nrht_tmp[nb];
               if (nrb<0 || nrb>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nrht of nbot %d\n",mype,__LINE__,ic,nrb);
            }
            int nt = ntop_tmp[ic];
            if (nt<0 || nt>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d global %d ntop %d ncells %ld ncells_ghost %ld\n",mype,__LINE__,ic,ic+noffset,nt,ncells,ncells_ghost);
            if (level_tmp[nt] > level_tmp[ic]){
               int nrt = nrht_tmp[nt];
               if (nrt<0 || nrt>= (int)ncells_ghost) printf("%d: Warning at line %d cell %d nrht of ntop %d\n",mype,__LINE__,ic,nrt);
            }
         }
      }
#endif

      if (TIMING_LEVEL >= 2) {
         gpu_timers[MESH_TIMER_SETUP_COMM] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      }

      if (DEBUG) {
         print_dev_local();

         vector<int> hash_tmp(hashsize);
         ezcl_enqueue_read_buffer(command_queue, dev_hash, CL_FALSE, 0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);

         cl_mem dev_hash_header_check = gpu_get_hash_header();
         vector<ulong> hash_header_check(hash_header_size);
         ezcl_enqueue_read_buffer(command_queue, dev_hash_header_check, CL_TRUE, 0, hash_header_size*sizeof(cl_ulong), &hash_header_check[0], NULL);

         int   gpu_hash_method     = (int)hash_header_check[0];
         ulong gpu_hash_table_size =      hash_header_check[1];
         ulong gpu_AA              =      hash_header_check[2];
         ulong gpu_BB              =      hash_header_check[3];

         vector<int> nlft_tmp(ncells_ghost);
         vector<int> nrht_tmp(ncells_ghost);
         vector<int> nbot_tmp(ncells_ghost);
         vector<int> ntop_tmp(ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0], NULL);

         int jmaxglobal = (jmax+1)*IPOW2(levmx);
         int imaxglobal = (imax+1)*IPOW2(levmx);
         fprintf(fp,"\n                                    HASH numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     fprintf(fp,"%5d",read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) );
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");

         fprintf(fp,"\n                                    nlft numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) -noffset;
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",nlft_tmp[hashval]);
                  } else {
                        fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
      
         fprintf(fp,"\n                                    nrht numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) -noffset;
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                     fprintf(fp,"%5d",nrht_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");

         fprintf(fp,"\n                                    nbot numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) -noffset;
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                     fprintf(fp,"%5d",nbot_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");

         fprintf(fp,"\n                                    ntop numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  int hashval = read_dev_hash(gpu_hash_method, gpu_hash_table_size, gpu_AA, gpu_BB, (jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize), &hash_tmp[0]) -noffset;
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                     fprintf(fp,"%5d",ntop_tmp[hashval]);
                  } else {
                     fprintf(fp,"     ");
                  }
               }
            }
            fprintf(fp,"\n");
         }
         fprintf(fp,"%2d:      ",mype);
         for (int ii = 0; ii<imaxglobal; ii++){
            fprintf(fp,"%4d:",ii);
         }
         fprintf(fp,"\n");
      }

      if (DEBUG) {
         print_dev_local();

         vector<int> i_tmp(ncells_ghost);
         vector<int> j_tmp(ncells_ghost);
         vector<uchar_t> level_tmp(ncells_ghost);
         vector<int> nlft_tmp(ncells_ghost);
         vector<int> nrht_tmp(ncells_ghost);
         vector<int> nbot_tmp(ncells_ghost);
         vector<int> ntop_tmp(ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_i, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &i_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_j, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &j_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level, CL_FALSE, 0, ncells_ghost*sizeof(cl_uchar_t), &level_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0], NULL);

         for (uint ic=0; ic<ncells; ic++){
            fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
                mype,ic,i_tmp[ic],j_tmp[ic],level_tmp[ic],nlft_tmp[ic],nrht_tmp[ic],nbot_tmp[ic],ntop_tmp[ic]);
         }
         int ig=0;
         for (uint ic=ncells; ic<ncells_ghost; ic++, ig++){
            fprintf(fp,"%d: after  update ic %d off %d i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
                mype,ic,indices_needed[ig],i_tmp[ic],j_tmp[ic],level_tmp[ic],nlft_tmp[ic],nrht_tmp[ic],nbot_tmp[ic],ntop_tmp[ic]);
         }
      }
   }
#endif

   ezcl_device_memory_delete(dev_sizes);
   ezcl_device_memory_delete(dev_check);

   gpu_compact_hash_delete(dev_hash, dev_hash_header);

   gpu_timers[MESH_TIMER_CALC_NEIGHBORS] += (long)(cpu_timer_stop(tstart_cpu) * 1.0e9);
}
#endif

void Mesh::print_calc_neighbor_type(void)
{
   if ( calc_neighbor_type == HASH_TABLE ) {
      if (mype == 0) printf("Using hash tables to calculate neighbors\n");
      if (mype == 0 && numpe == 1) final_hash_collision_report();
   } else {
      printf("hash table size %ld\n",ncells*(int)log(ncells)*sizeof(int));
      if (mype == 0) printf("Using k-D tree to calculate neighbors\n");
   }
}

int Mesh::get_calc_neighbor_type(void)
{
   return(calc_neighbor_type );
}

void Mesh::calc_celltype_threaded(size_t ncells)
{
   int flags=0;
#ifdef HAVE_J7
   if (parallel) flags = LOAD_BALANCE_MEMORY;
#endif

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
   if (celltype == NULL || mesh_memory.get_memory_size(celltype) < ncells) {
      if (celltype != NULL) celltype = (char_t *)mesh_memory.memory_delete((void *)celltype);
      celltype = (char_t *)mesh_memory.memory_malloc(ncells, sizeof(char_t), "celltype", flags);
   }
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif
   for (uint ic=0; ic<ncells; ++ic) {
      celltype[ic] = REAL_CELL;
      if (is_left_boundary(ic) )   celltype[ic] = LEFT_BOUNDARY;
      if (is_right_boundary(ic) )  celltype[ic] = RIGHT_BOUNDARY;
      if (is_bottom_boundary(ic) ) celltype[ic] = BOTTOM_BOUNDARY;
      if (is_top_boundary(ic))     celltype[ic] = TOP_BOUNDARY;
   }
}

void Mesh::calc_celltype(size_t ncells)
{
   int flags = 0;
#ifdef HAVE_J7
   if (parallel) flags = LOAD_BALANCE_MEMORY;
#endif

   if (celltype == NULL || mesh_memory.get_memory_size(celltype) < ncells) {
      if (celltype != NULL) celltype = (char_t *)mesh_memory.memory_delete((void *)celltype);
      celltype = (char_t *)mesh_memory.memory_malloc(ncells, sizeof(char_t), "celltype", flags);
   }

#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (uint ic=0; ic<ncells; ++ic) {
      celltype[ic] = REAL_CELL;
      if (is_left_boundary(ic) )   celltype[ic] = LEFT_BOUNDARY;
      if (is_right_boundary(ic) )  celltype[ic] = RIGHT_BOUNDARY;
      if (is_bottom_boundary(ic) ) celltype[ic] = BOTTOM_BOUNDARY;
      if (is_top_boundary(ic))     celltype[ic] = TOP_BOUNDARY;
   }
}

void Mesh::calc_symmetry(vector<int> &dsym, vector<int> &xsym, vector<int> &ysym)
{
   TBounds box;
   vector<int> index_list( IPOW2(levmx*levmx) );

   int num;
   for (uint ic=0; ic<ncells; ic++) {
      dsym[ic]=ic;
      xsym[ic]=ic;
      ysym[ic]=ic;

      //diagonal symmetry
      box.min.x = -1.0*(x[ic]+0.5*dx[ic]);
      box.max.x = -1.0*(x[ic]+0.5*dx[ic]);
      box.min.y = -1.0*(y[ic]+0.5*dy[ic]);
      box.max.y = -1.0*(y[ic]+0.5*dy[ic]);
      KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
      if (num == 1) dsym[ic]=index_list[0];
      //printf("ic %d dsym[ic] %d num %d\n",ic,dsym[ic],num);

      //x symmetry
      box.min.x = -1.0*(x[ic]+0.5*dx[ic]);
      box.max.x = -1.0*(x[ic]+0.5*dx[ic]);
      box.min.y = y[ic]+0.5*dy[ic];
      box.max.y = y[ic]+0.5*dy[ic];
      KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
      if (num == 1) xsym[ic]=index_list[0];

      //y symmetry
      box.min.x = x[ic]+0.5*dx[ic];
      box.max.x = x[ic]+0.5*dx[ic];
      box.min.y = -1.0*(y[ic]+0.5*dy[ic]);
      box.max.y = -1.0*(y[ic]+0.5*dy[ic]);
      KDTree_QueryBoxIntersect(&tree, &num, &(index_list[0]), &box);
      if (num == 1) ysym[ic]=index_list[0];

   }
}

#ifdef HAVE_MPI
void Mesh::do_load_balance_local(size_t numcells, float *weight, MallocPlus &state_memory)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   // To get rid of compiler warning
   if (DEBUG && weight != NULL) printf("DEBUG weight[0] = %f\n",weight[0]);

   int ncells_old = numcells;
   int noffset_old = ndispl[mype];

// Need to add weight array to load balance if it is not NULL
// Need to add tolerance to when load balance is done

   int do_load_balance_global = 0;
   int nsizes_old = 0;

   for (int ip=0; ip<numpe; ip++){
      nsizes_old = nsizes[ip];

      // Calc new,even partition of data across processors
      nsizes[ip] = ncells_global/numpe;
      // Account for leftover cells
      if (ip < (int)(ncells_global%numpe)) nsizes[ip]++;

      if (nsizes_old != nsizes[ip]) do_load_balance_global = 1;
   }

   if (do_load_balance_global) {
      cpu_counters[MESH_COUNTER_LOAD_BALANCE]++;

      mesh_memory.memory_delete(celltype);
      mesh_memory.memory_delete(nlft);
      mesh_memory.memory_delete(nrht);
      mesh_memory.memory_delete(nbot);
      mesh_memory.memory_delete(ntop);

      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }
      ncells = nsizes[mype];
      noffset=ndispl[mype];

      // Indices of blocks to be added to load balance
      int lower_block_start = noffset;
      int lower_block_end   = min(noffset_old-1, (int)(noffset+ncells-1));
      int upper_block_start = max((int)(noffset_old+ncells_old), noffset);
      int upper_block_end   = noffset+ncells-1;

      int lower_block_size = max(lower_block_end-lower_block_start+1,0);
      if(lower_block_end < 0) lower_block_size = 0; // Handles segfault at start of array
      int upper_block_size = max(upper_block_end-upper_block_start+1,0);
      int indices_needed_count = lower_block_size + upper_block_size;

      int in = 0;

      vector<int> indices_needed(indices_needed_count);
      for (int iz = lower_block_start; iz <= lower_block_end; iz++, in++){
         indices_needed[in]=iz;
      }
      for (int iz = upper_block_start; iz <= upper_block_end; iz++, in++){
         indices_needed[in]=iz;
      }

      int load_balance_handle = 0;
      L7_Setup(0, noffset_old, ncells_old, &indices_needed[0], indices_needed_count, &load_balance_handle);

      //printf("\n%d: DEBUG load balance report\n",mype);

      state_memory.memory_realloc(ncells_old+indices_needed_count,LOAD_BALANCE_MEMORY);

      MallocPlus state_memory_old = state_memory;


      malloc_plus_memory_entry *memory_item;

      for (memory_item = state_memory_old.memory_entry_by_name_begin();
           memory_item != state_memory_old.memory_entry_by_name_end();
           memory_item = state_memory_old.memory_entry_by_name_next() ) {

         //if (mype == 0) printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

         if ( (memory_item->mem_flags & LOAD_BALANCE_MEMORY) == 0) continue;

         if (memory_item->mem_elsize == 8) {
            double *mem_ptr_double = (double *)memory_item->mem_ptr;

            int flags = state_memory.get_memory_flags(mem_ptr_double);
            double *state_temp_double = (double *) state_memory.memory_malloc(ncells, sizeof(double),
                                                                              "state_temp_double", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,memory_item->mem_ptr);
            L7_Update(mem_ptr_double, L7_DOUBLE, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
               for(in = 0; in < MIN(lower_block_size, (int)ncells); in++) {
                  state_temp_double[in] = mem_ptr_double[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               state_temp_double[in] = mem_ptr_double[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  state_temp_double[in] = mem_ptr_double[ic+k];
               }
            }
            state_memory.memory_replace(mem_ptr_double, state_temp_double);
         } else if (memory_item->mem_elsize == 4) {
            float *mem_ptr_float = (float *)memory_item->mem_ptr;

            int flags = state_memory.get_memory_flags(mem_ptr_float);
            float *state_temp_float = (float *) state_memory.memory_malloc(ncells, sizeof(float),
                                                                           "state_temp_float", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,mem_ptr_float);
            L7_Update(mem_ptr_float, L7_FLOAT, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
               for(in = 0; in < MIN(lower_block_size, (int)ncells); in++) {
                  state_temp_float[in] = mem_ptr_float[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               state_temp_float[in] = mem_ptr_float[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  state_temp_float[in] = mem_ptr_float[ic+k];
               }
            }
            state_memory.memory_replace(mem_ptr_float, state_temp_float);
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            half *mem_ptr_half = (half *)memory_item->mem_ptr;

            int flags = state_memory.get_memory_flags(mem_ptr_half);
            half *state_temp_half = (half *) state_memory.memory_malloc(ncells, sizeof(half),
                                                                           "state_temp_half", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,mem_ptr_half);
            L7_Update(mem_ptr_half, L7_FLOAT, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
               for(; in < MIN(lower_block_size, (int)ncells); in++) {
                  state_temp_half[in] = mem_ptr_half[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               state_temp_half[in] = mem_ptr_half[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  state_temp_half[in] = mem_ptr_half[ic+k];
               }
            }
            state_memory.memory_replace(mem_ptr_half, state_temp_half);
#endif
         }
      }

      mesh_memory.memory_realloc_all(ncells_old+indices_needed_count);

      MallocPlus mesh_memory_old = mesh_memory;

      for (memory_item = mesh_memory_old.memory_entry_by_name_begin();
           memory_item != mesh_memory_old.memory_entry_by_name_end();
           memory_item = mesh_memory_old.memory_entry_by_name_next() ) {

         //if (mype == 0) printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

         if (memory_item->mem_elsize == 8) {
            long long *mem_ptr_long = (long long *)memory_item->mem_ptr;

            int flags = mesh_memory.get_memory_flags(mem_ptr_long);
            long long *mesh_temp_long = (long long *)mesh_memory.memory_malloc(ncells, sizeof(long long), "mesh_temp_long", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,mem_ptr);
            L7_Update(mem_ptr_long, L7_LONG_LONG_INT, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
               for(in = 0; in < MIN(lower_block_size, (int)ncells); in++) {
                  mesh_temp_long[in] = mem_ptr_long[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               mesh_temp_long[in] = mem_ptr_long[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  mesh_temp_long[in] = mem_ptr_long[ic+k];
               }
            }
            mesh_memory.memory_replace(mem_ptr_long, mesh_temp_long);

         } else if (memory_item->mem_elsize == 4) {
            int *mem_ptr_int = (int *)memory_item->mem_ptr;

            int flags = mesh_memory.get_memory_flags(mem_ptr_int);
            int *mesh_temp_int = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "mesh_temp_int", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,mem_ptr);
            L7_Update(mem_ptr_int, L7_INT, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
               for(in = 0; in < MIN(lower_block_size, (int)ncells); in++) {
                  mesh_temp_int[in] = mem_ptr_int[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               mesh_temp_int[in] = mem_ptr_int[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  mesh_temp_int[in] = mem_ptr_int[ic+k];
               }
            }
            mesh_memory.memory_replace(mem_ptr_int, mesh_temp_int);

         } else if (memory_item->mem_elsize == 2) {
            short *mem_ptr_short = (short *)memory_item->mem_ptr;

            int flags = mesh_memory.get_memory_flags(mem_ptr_short);
            short *mesh_temp_short = (short *)mesh_memory.memory_malloc(ncells, sizeof(short), "mesh_temp_short", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,mem_ptr);
            L7_Update(mem_ptr_short, L7_SHORT, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
               for(in = 0; in < MIN(lower_block_size, (int)ncells); in++) {
                  mesh_temp_short[in] = mem_ptr_short[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               mesh_temp_short[in] = mem_ptr_short[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  mesh_temp_short[in] = mem_ptr_short[ic+k];
               }
            }
            mesh_memory.memory_replace(mem_ptr_short, mesh_temp_short);

         } else if (memory_item->mem_elsize == 1) {
            char *mem_ptr_char = (char *)memory_item->mem_ptr;

            int flags = mesh_memory.get_memory_flags(mem_ptr_char);
            char *mesh_temp_char = (char *)mesh_memory.memory_malloc(ncells, sizeof(char), "mesh_temp_char", flags);

            //printf("%d: DEBUG L7_Update in do_load_balance_local mem_ptr %p\n",mype,mem_ptr);
            L7_Update(mem_ptr_char, L7_CHAR, load_balance_handle);
            in = 0;
            if(lower_block_size > 0) {
               for(; in < MIN(lower_block_size, (int)ncells); in++) {
                  mesh_temp_char[in] = mem_ptr_char[ncells_old + in];
               }
            }

            for(int ic = MAX((noffset - noffset_old), 0); (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
               mesh_temp_char[in] = mem_ptr_char[ic];
            }

            if(upper_block_size > 0) {
               int ic = ncells_old + lower_block_size;
               for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
                  mesh_temp_char[in] = mem_ptr_char[ic+k];
               }
            }
            mesh_memory.memory_replace(mem_ptr_char, mesh_temp_char);

         }
      }

      L7_Free(&load_balance_handle);
      load_balance_handle = 0;

      memory_reset_ptrs();

      //mesh_memory.memory_report();
      //state_memory.memory_report();
      //printf("%d: DEBUG end load balance report\n\n",mype);
      calc_celltype(ncells);
   }
   //if (mype == 0) printf("DEBUG -- finished load balance\n");

   cpu_timers[MESH_TIMER_LOAD_BALANCE] += cpu_timer_stop(tstart_cpu);
}
#endif

#ifdef HAVE_OPENCL
#ifdef HAVE_MPI
int Mesh::gpu_do_load_balance_local(size_t numcells, float *weight, MallocPlus &gpu_state_memory)
{
   int do_load_balance_global = 0;

   if (! gpu_do_rezone) return(do_load_balance_global);

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   // To get rid of compiler warning
   if (DEBUG && weight != NULL) printf("DEBUG weight[0] = %f\n",weight[0]);

   int ncells_old = numcells;
   int noffset_old = ndispl[mype];

// Need to add weight array to load balance if it is not NULL
// Need to add tolerance to when load balance is done

   int nsizes_old = 0;
   for (int ip=0; ip<numpe; ip++){
      nsizes_old = nsizes[ip];
      nsizes[ip] = ncells_global/numpe;
      if (ip < (int)(ncells_global%numpe)) nsizes[ip]++;
      if (nsizes_old != nsizes[ip]) do_load_balance_global = 1;
   }

   if(do_load_balance_global) {

      cl_command_queue command_queue = ezcl_get_command_queue();

      gpu_counters[MESH_COUNTER_LOAD_BALANCE]++;

      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }
      ncells = nsizes[mype];
      noffset=ndispl[mype];

      // Indices of blocks to be added to load balance
      int lower_block_start = noffset;
      int lower_block_end   = min(noffset_old-1, (int)(noffset+ncells-1));
      int upper_block_start = max((int)(noffset_old+ncells_old), noffset);
      int upper_block_end   = noffset+ncells-1;
      //printf("%d: lbs %d lbe %d ubs %d ube %d\n",mype,lower_block_start-noffset_old,lower_block_end-noffset_old,upper_block_start-noffset_old,upper_block_end-noffset_old);

      size_t lower_block_size = max(lower_block_end-lower_block_start+1,0);
      if(lower_block_end < 0) lower_block_size = 0; // Handles segfault at start of array
      size_t upper_block_size = max(upper_block_end-upper_block_start+1,0);
      int indices_needed_count = lower_block_size + upper_block_size;

      size_t middle_block_size = ncells - lower_block_size - upper_block_size;
      int middle_block_start = max(noffset - noffset_old, 0);

      int lower_segment_size = noffset-noffset_old;
      int do_whole_segment = 0;
      if (lower_segment_size > ncells_old) do_whole_segment = 1;

      int upper_segment_size = ( (noffset_old+ncells_old) - (noffset+ncells) );
      int upper_segment_start = (noffset_old+ncells_old) - upper_segment_size - noffset_old;
      if (upper_segment_size > ncells_old) do_whole_segment=1;

      int in = 0;
      vector<int> indices_needed(indices_needed_count);
      for (int iz = lower_block_start; iz <= lower_block_end; iz++, in++){
         indices_needed[in]=iz;
      }
      for (int iz = upper_block_start; iz <= upper_block_end; iz++, in++){
         indices_needed[in]=iz;
      }

      int load_balance_handle = 0;
      L7_Setup(0, noffset_old, ncells_old, &indices_needed[0], indices_needed_count, &load_balance_handle);
       
      size_t local_work_size = 128;
      size_t global_work_size = ((ncells + local_work_size - 1) / local_work_size) * local_work_size;

      // printf("MYPE%d: \t ncells = %d \t ncells_old = %d \t ncells_global = %d \n", mype, ncells, ncells_old, ncells_global);

      // Allocate lower block on GPU
      size_t low_block_size = MAX(1, lower_block_size);
      cl_mem dev_state_var_lower = ezcl_malloc(NULL, const_cast<char *>("dev_state_var_lower"), &low_block_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);

      // Allocate upper block on GPU
      size_t up_block_size = MAX(1, upper_block_size);
      cl_mem dev_state_var_upper = ezcl_malloc(NULL, const_cast<char *>("dev_state_var_upper"), &up_block_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);

      MallocPlus gpu_state_memory_old = gpu_state_memory;
      malloc_plus_memory_entry *memory_item;

      for (memory_item = gpu_state_memory_old.memory_entry_by_name_begin();
           memory_item != gpu_state_memory_old.memory_entry_by_name_end();
           memory_item = gpu_state_memory_old.memory_entry_by_name_next() ) {
         //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);
         cl_mem dev_state_mem_ptr = (cl_mem)memory_item->mem_ptr;

         if (memory_item->mem_elsize == 8){
#ifndef MINIMUM_PRECISION
            vector<double> state_var_tmp(ncells_old+indices_needed_count,0.0);

            // Read current state values from GPU and write to CPU arrays
            if (do_whole_segment) {
               ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr, CL_TRUE, 0, ncells_old*sizeof(cl_double), &state_var_tmp[0], NULL);
            } else {
               // Read lower block from GPU
               if (lower_segment_size > 0) {
                  ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr, CL_TRUE, 0, lower_segment_size*sizeof(cl_double), &state_var_tmp[0], NULL);
               }
               // Read upper block from GPU
               if (upper_segment_size > 0) {
                  ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr, CL_TRUE, upper_segment_start*sizeof(cl_double), upper_segment_size*sizeof(cl_double), &state_var_tmp[upper_segment_start], NULL);
               }
            }

            // Update arrays with L7
            L7_Update(&state_var_tmp[0], L7_DOUBLE, load_balance_handle);

            // Set lower block on GPU
            if(lower_block_size > 0) {
               ezcl_enqueue_write_buffer(command_queue, dev_state_var_lower, CL_FALSE, 0, lower_block_size*sizeof(cl_double), &state_var_tmp[ncells_old], NULL);
            }
            // Set upper block on GPU
            if(upper_block_size > 0) {
               ezcl_enqueue_write_buffer(command_queue, dev_state_var_upper, CL_FALSE, 0, upper_block_size*sizeof(cl_double), &state_var_tmp[ncells_old+lower_block_size], NULL); 
            }

            // Allocate space on GPU for temp arrays (used in double buffering)
            cl_mem dev_state_var_new = ezcl_malloc(NULL, gpu_state_memory.get_memory_name(dev_state_mem_ptr), &ncells, sizeof(cl_double), CL_MEM_READ_WRITE, 0);
            gpu_state_memory.memory_add(dev_state_var_new, ncells, sizeof(cl_double), "dev_state_var_new", DEVICE_REGULAR_MEMORY);

            //printf("DEBUG memory for proc %d is %p dev_state_new is %p\n",mype,dev_state_mem_ptr,dev_state_var_new);

            ezcl_set_kernel_arg(kernel_do_load_balance_double, 0, sizeof(cl_int), &ncells);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 1, sizeof(cl_int), &lower_block_size);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 2, sizeof(cl_int), &middle_block_size);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 3, sizeof(cl_int), &middle_block_start);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 4, sizeof(cl_mem), &dev_state_mem_ptr);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 5, sizeof(cl_mem), &dev_state_var_lower);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 6, sizeof(cl_mem), &dev_state_var_upper);
            ezcl_set_kernel_arg(kernel_do_load_balance_double, 7, sizeof(cl_mem), &dev_state_var_new);

            ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_double,   1, NULL, &global_work_size, &local_work_size, NULL);

            gpu_state_memory.memory_replace(dev_state_mem_ptr, dev_state_var_new);
#else
            printf("ERROR -- can't have double type for state variable\n");
            exit(1);
#endif
         } else if (memory_item->mem_elsize == 4){
            vector<float> state_var_tmp(ncells_old+indices_needed_count,0.0);

            // Read current state values from GPU and write to CPU arrays
            if (do_whole_segment) {
               ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr, CL_TRUE, 0, ncells_old*sizeof(cl_float), &state_var_tmp[0], NULL);
            } else {
               // Read lower block from GPU
               if (lower_segment_size > 0) {
                  ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr, CL_TRUE, 0, lower_segment_size*sizeof(cl_float), &state_var_tmp[0], NULL);
               }
               // Read upper block from GPU
               if (upper_segment_size > 0) {
                  ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr, CL_TRUE, upper_segment_start*sizeof(cl_float), upper_segment_size*sizeof(cl_float), &state_var_tmp[upper_segment_start], NULL);
               }
            }

            // Update arrays with L7
            L7_Update(&state_var_tmp[0], L7_FLOAT, load_balance_handle);

            // Set lower block on GPU
            if(lower_block_size > 0) {
               ezcl_enqueue_write_buffer(command_queue, dev_state_var_lower, CL_FALSE, 0, lower_block_size*sizeof(cl_float), &state_var_tmp[ncells_old], NULL);
            }
            // Set upper block on GPU
            if(upper_block_size > 0) {
               ezcl_enqueue_write_buffer(command_queue, dev_state_var_upper, CL_FALSE, 0, upper_block_size*sizeof(cl_float), &state_var_tmp[ncells_old+lower_block_size], NULL); 
            }

            // Allocate space on GPU for temp arrays (used in double buffering)
            cl_mem dev_state_var_new = ezcl_malloc(NULL, gpu_state_memory.get_memory_name(dev_state_mem_ptr), &ncells, sizeof(cl_float), CL_MEM_READ_WRITE, 0);
            gpu_state_memory.memory_add(dev_state_var_new, ncells, sizeof(cl_float), "dev_state_var_new", DEVICE_REGULAR_MEMORY);

            //printf("DEBUG memory for proc %d is %p dev_state_new is %p\n",mype,dev_state_mem_ptr,dev_state_var_new);

            ezcl_set_kernel_arg(kernel_do_load_balance_float, 0, sizeof(cl_int), &ncells);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 1, sizeof(cl_int), &lower_block_size);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 2, sizeof(cl_int), &middle_block_size);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 3, sizeof(cl_int), &middle_block_start);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 4, sizeof(cl_mem), &dev_state_mem_ptr);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 5, sizeof(cl_mem), &dev_state_var_lower);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 6, sizeof(cl_mem), &dev_state_var_upper);
            ezcl_set_kernel_arg(kernel_do_load_balance_float, 7, sizeof(cl_mem), &dev_state_var_new);

            ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_float,   1, NULL, &global_work_size, &local_work_size, NULL);

            gpu_state_memory.memory_replace(dev_state_mem_ptr, dev_state_var_new);
         }
      }

      ezcl_device_memory_delete(dev_state_var_lower);
      ezcl_device_memory_delete(dev_state_var_upper);

      vector<int> i_tmp(ncells_old+indices_needed_count,0);
      vector<int> j_tmp(ncells_old+indices_needed_count,0);
      vector<uchar_t> level_tmp(ncells_old+indices_needed_count,0);
      vector<char_t> celltype_tmp(ncells_old+indices_needed_count,0);

      if (do_whole_segment) {
         ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells_old*sizeof(cl_int), &i_tmp[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells_old*sizeof(cl_int), &j_tmp[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells_old*sizeof(cl_uchar_t), &level_tmp[0],    NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  0, ncells_old*sizeof(cl_char_t), &celltype_tmp[0], NULL);
      } else {
         if (lower_segment_size > 0) {
            ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, lower_segment_size*sizeof(cl_int), &i_tmp[0],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, lower_segment_size*sizeof(cl_int), &j_tmp[0],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, lower_segment_size*sizeof(cl_uchar_t), &level_tmp[0],    NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  0, lower_segment_size*sizeof(cl_char_t), &celltype_tmp[0], NULL);
         }
         if (upper_segment_size > 0) {
            ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, upper_segment_start*sizeof(cl_int), upper_segment_size*sizeof(cl_int), &i_tmp[upper_segment_start],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, upper_segment_start*sizeof(cl_int), upper_segment_size*sizeof(cl_int), &j_tmp[upper_segment_start],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, upper_segment_start*sizeof(cl_uchar_t), upper_segment_size*sizeof(cl_uchar_t), &level_tmp[upper_segment_start],    NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  upper_segment_start*sizeof(cl_char_t), upper_segment_size*sizeof(cl_char_t), &celltype_tmp[upper_segment_start], NULL);
         }
      }

      L7_Update(&i_tmp[0],        L7_INT, load_balance_handle);
      L7_Update(&j_tmp[0],        L7_INT, load_balance_handle);
      L7_Update(&level_tmp[0],    L7_UCHAR_T, load_balance_handle);
      L7_Update(&celltype_tmp[0], L7_CHAR_T, load_balance_handle);

      L7_Free(&load_balance_handle);
      load_balance_handle = 0;

      // Allocate and set lower block on GPU
      cl_mem dev_i_lower, dev_j_lower, dev_level_lower, dev_celltype_lower;

      if(lower_block_size > 0) {
         dev_i_lower        = ezcl_malloc(NULL, const_cast<char *>("dev_i_lower"),        &lower_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_j_lower        = ezcl_malloc(NULL, const_cast<char *>("dev_j_lower"),        &lower_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_level_lower    = ezcl_malloc(NULL, const_cast<char *>("dev_level_lower"),    &lower_block_size, sizeof(cl_uchar_t), CL_MEM_READ_WRITE, 0);
         dev_celltype_lower = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_lower"), &lower_block_size, sizeof(cl_char_t), CL_MEM_READ_WRITE, 0);

         ezcl_enqueue_write_buffer(command_queue, dev_i_lower,        CL_FALSE, 0, lower_block_size*sizeof(cl_int), &i_tmp[ncells_old],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_j_lower,        CL_FALSE, 0, lower_block_size*sizeof(cl_int), &j_tmp[ncells_old],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_level_lower,    CL_FALSE, 0, lower_block_size*sizeof(cl_uchar_t), &level_tmp[ncells_old],    NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_celltype_lower, CL_TRUE,  0, lower_block_size*sizeof(cl_char_t), &celltype_tmp[ncells_old], NULL);
      }

      // Allocate and set upper block on GPU
      cl_mem dev_i_upper, dev_j_upper, dev_level_upper, dev_celltype_upper;
      if(upper_block_size > 0) {
         dev_i_upper        = ezcl_malloc(NULL, const_cast<char *>("dev_i_upper"),        &upper_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_j_upper        = ezcl_malloc(NULL, const_cast<char *>("dev_j_upper"),        &upper_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_level_upper    = ezcl_malloc(NULL, const_cast<char *>("dev_level_upper"),    &upper_block_size, sizeof(cl_uchar_t), CL_MEM_READ_WRITE, 0);
         dev_celltype_upper = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_upper"), &upper_block_size, sizeof(cl_char_t), CL_MEM_READ_WRITE, 0);

         ezcl_enqueue_write_buffer(command_queue, dev_i_upper,        CL_FALSE, 0, upper_block_size*sizeof(cl_int), &i_tmp[ncells_old+lower_block_size],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_j_upper,        CL_FALSE, 0, upper_block_size*sizeof(cl_int), &j_tmp[ncells_old+lower_block_size],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_level_upper,    CL_FALSE, 0, upper_block_size*sizeof(cl_uchar_t), &level_tmp[ncells_old+lower_block_size],    NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_celltype_upper, CL_TRUE,  0, upper_block_size*sizeof(cl_char_t), &celltype_tmp[ncells_old+lower_block_size], NULL);
      }

      local_work_size = 128;

      // printf("MYPE%d: \t ncells = %d \t ncells_old = %d \t ncells_global = %d \n", mype, ncells, ncells_old, ncells_global);
      // Allocate space on GPU for temp arrays (used in double buffering)

      size_t mem_request = (int)((float)ncells*mem_factor);
      cl_mem dev_i_new        = ezcl_malloc(NULL, const_cast<char *>("dev_i_new"),        &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_j_new        = ezcl_malloc(NULL, const_cast<char *>("dev_j_new"),        &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_level_new    = ezcl_malloc(NULL, const_cast<char *>("dev_level_new"),    &mem_request, sizeof(cl_uchar_t),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_celltype_new = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_new"), &mem_request, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      // Set kernel arguments and call lower block kernel
      if(lower_block_size > 0) {

         size_t global_work_size = ((lower_block_size + local_work_size - 1) / local_work_size) * local_work_size;

         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 0, sizeof(cl_mem), &dev_i_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 1, sizeof(cl_mem), &dev_j_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 2, sizeof(cl_mem), &dev_level_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 3, sizeof(cl_mem), &dev_celltype_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 4, sizeof(cl_mem), &dev_i_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 5, sizeof(cl_mem), &dev_j_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 6, sizeof(cl_mem), &dev_level_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 7, sizeof(cl_mem), &dev_celltype_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 8, sizeof(cl_int), &lower_block_size);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_lower,   1, NULL, &global_work_size, &local_work_size, NULL);

         ezcl_device_memory_delete(dev_i_lower);
         ezcl_device_memory_delete(dev_j_lower);
         ezcl_device_memory_delete(dev_level_lower);
         ezcl_device_memory_delete(dev_celltype_lower);
      }

      // Set kernel arguments and call middle block kernel
      if(middle_block_size > 0) {

         size_t global_work_size = ((middle_block_size + local_work_size - 1) / local_work_size) * local_work_size;

         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  0, sizeof(cl_mem), &dev_i_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  1, sizeof(cl_mem), &dev_j_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  2, sizeof(cl_mem), &dev_level_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  3, sizeof(cl_mem), &dev_celltype_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  4, sizeof(cl_mem), &dev_i);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  5, sizeof(cl_mem), &dev_j);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  6, sizeof(cl_mem), &dev_level);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  7, sizeof(cl_mem), &dev_celltype);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  8, sizeof(cl_int), &lower_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle,  9, sizeof(cl_int), &middle_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 10, sizeof(cl_int), &middle_block_start);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_middle,   1, NULL, &global_work_size, &local_work_size, NULL);
      }

      // Set kernel arguments and call upper block kernel
      if(upper_block_size > 0) {

         size_t global_work_size = ((upper_block_size + local_work_size - 1) / local_work_size) * local_work_size;

         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  0, sizeof(cl_mem), &dev_i_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  1, sizeof(cl_mem), &dev_j_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  2, sizeof(cl_mem), &dev_level_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  3, sizeof(cl_mem), &dev_celltype_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  4, sizeof(cl_mem), &dev_i_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  5, sizeof(cl_mem), &dev_j_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  6, sizeof(cl_mem), &dev_level_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  7, sizeof(cl_mem), &dev_celltype_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  8, sizeof(cl_int), &lower_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper,  9, sizeof(cl_int), &middle_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 10, sizeof(cl_int), &upper_block_size);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_upper,   1, NULL, &global_work_size, &local_work_size, NULL);

         ezcl_device_memory_delete(dev_i_upper);
         ezcl_device_memory_delete(dev_j_upper);
         ezcl_device_memory_delete(dev_level_upper);
         ezcl_device_memory_delete(dev_celltype_upper);
      }

      ezcl_device_memory_swap(&dev_i_new,        &dev_i);
      ezcl_device_memory_swap(&dev_j_new,        &dev_j);
      ezcl_device_memory_swap(&dev_level_new,    &dev_level);
      ezcl_device_memory_swap(&dev_celltype_new, &dev_celltype);

      ezcl_device_memory_delete(dev_i_new);
      ezcl_device_memory_delete(dev_j_new);
      ezcl_device_memory_delete(dev_level_new);
      ezcl_device_memory_delete(dev_celltype_new);

      gpu_timers[MESH_TIMER_LOAD_BALANCE] += (long int)(cpu_timer_stop(tstart_cpu)*1.0e9);
   }

   return(do_load_balance_global);
}
#endif
#endif

#ifdef HAVE_OPENCL
int Mesh::gpu_count_BCs(void)
{
   cl_event count_BCs_stage1_event, count_BCs_stage2_event;

   size_t local_work_size  = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

   //size_t block_size = (ncells + TILE_SIZE - 1) / TILE_SIZE; //  For on-device global reduction kernel.
   size_t block_size     = global_work_size/local_work_size;

   int bcount = 0;

   if (! have_boundary) {
      cl_command_queue command_queue = ezcl_get_command_queue();
      cl_mem dev_ioffset  = ezcl_malloc(NULL, const_cast<char *>("dev_ioffset"), &block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

       /*
       __kernel void count_BCs(
                        const int    isize,      // 0   
               __global const int   *i,         // 1
               __global const int   *j,         // 2
               __global const uchar_t *level,     // 3
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

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_count_BCs, 1, NULL, &global_work_size, &local_work_size, &count_BCs_stage1_event);

      if (block_size > 1) {
         ezcl_set_kernel_arg(kernel_reduce_sum_int_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
         ezcl_set_kernel_arg(kernel_reduce_sum_int_stage2of2, 1, sizeof(cl_mem), (void *)&dev_ioffset);
         ezcl_set_kernel_arg(kernel_reduce_sum_int_stage2of2, 2, shared_spd_sum_int, 0);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_int_stage2of2, 1, NULL, &local_work_size, &local_work_size, &count_BCs_stage2_event);
      }

      ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, 1*sizeof(cl_int), &bcount, NULL);
 
      //printf("DEBUG -- bcount is %d\n",bcount);
      //state->gpu_time_read += ezcl_timer_calc(&start_read_event, &start_read_event);

      ezcl_device_memory_delete(dev_ioffset);

      gpu_timers[MESH_TIMER_COUNT_BCS] += ezcl_timer_calc(&count_BCs_stage1_event, &count_BCs_stage1_event);
      if (block_size > 1) {
         gpu_timers[MESH_TIMER_COUNT_BCS] += ezcl_timer_calc(&count_BCs_stage2_event, &count_BCs_stage2_event);
      }

   }

   return(bcount);
}
#endif

void Mesh::allocate(size_t ncells)
{
   int flags = 0;
   flags = RESTART_DATA;
#ifdef HAVE_J7
   if (parallel) flags = LOAD_BALANCE_MEMORY;
#endif

   i     = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "i",     flags);
   j     = (int *)mesh_memory.memory_malloc(ncells, sizeof(int), "j",     flags);
   level = (uchar_t *)mesh_memory.memory_malloc(ncells, sizeof(uchar_t), "level", flags);

   flags = 0;

   xface_i     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xface_i",     flags);
   xface_j     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xface_j",     flags);
   xface_level     = (uchar_t *)mesh_memory.memory_malloc(1, sizeof(uchar_t), "xface_level",     flags);
   map_xface2cell_lower     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_xface2cell_lower",     flags);
   map_xface2cell_upper     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_xface2cell_upper",     flags);
   map_xcell2face_left1     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_xcell2face_left1",     flags);
   map_xcell2face_left2     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_xcell2face_left2",     flags);
   map_xcell2face_right1     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_xcell2face_right1",     flags);
   map_xcell2face_right2     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_xcell2face_right2",     flags);
   xrecvIdx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xrecvIdx",     flags);
   xrecvCIdx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xrecvCIdx",     flags);
   xplusCell2Idx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xplusCell2Idx",     flags);
   xminusCell2Idx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xminusCell2Idx",     flags);
   xsendIdx1     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xsendIdx1",     flags);
   xsendIdx2     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "xsendIdx2",     flags);
   yface_i     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "yface_i",     flags);
   yface_j     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "yface_j",     flags);
   yface_level     = (uchar_t *)mesh_memory.memory_malloc(1, sizeof(uchar_t), "yface_level",     flags);
   map_yface2cell_lower     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_yface2cell_lower",     flags);
   map_yface2cell_upper     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_yface2cell_upper",     flags);
   map_ycell2face_bot1     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_ycell2face_bot1",     flags);
   map_ycell2face_bot2     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_ycell2face_bot2",     flags);
   map_ycell2face_top1     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_ycell2face_top1",     flags);
   map_ycell2face_top2     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "map_ycell2face_top2",     flags);
   yrecvIdx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "yrecvIdx",     flags);
   yrecvCIdx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "yrecvCIdx",     flags);
   yplusCell2Idx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "yplusCell2Idx",     flags);
   yminusCell2Idx     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "yminusCell2Idx",     flags);
   ysendIdx1     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "ysendIdx1",     flags);
   ysendIdx2     = (int *)mesh_memory.memory_malloc(1, sizeof(int), "ysendIdx2",     flags);
}


void Mesh::resize(size_t new_ncells)
{
   size_t current_size = mesh_memory.get_memory_size(i);
   if (new_ncells > current_size) mesh_memory.memory_realloc_all(new_ncells);
}

void Mesh::memory_reset_ptrs(void){
   i        = (int *)mesh_memory.get_memory_ptr("i");
   j        = (int *)mesh_memory.get_memory_ptr("j");
   level    = (uchar_t *)mesh_memory.get_memory_ptr("level");
   celltype = (char_t *)mesh_memory.get_memory_ptr("celltype");
   nlft     = (int *)mesh_memory.get_memory_ptr("nlft");
   nrht     = (int *)mesh_memory.get_memory_ptr("nrht");
   nbot     = (int *)mesh_memory.get_memory_ptr("nbot");
   ntop     = (int *)mesh_memory.get_memory_ptr("ntop");
   xface_i     = (int *)mesh_memory.get_memory_ptr("xface_i");
   xface_j     = (int *)mesh_memory.get_memory_ptr("xface_j");
   xface_level     = (uchar_t *)mesh_memory.get_memory_ptr("xface_level");
   map_xface2cell_lower     = (int *)mesh_memory.get_memory_ptr("map_xface2cell_lower");
   map_xface2cell_upper     = (int *)mesh_memory.get_memory_ptr("map_xface2cell_upper");
   map_xcell2face_left1     = (int *)mesh_memory.get_memory_ptr("map_xcell2face_left1");
   map_xcell2face_left2     = (int *)mesh_memory.get_memory_ptr("map_xcell2face_left2");
   map_xcell2face_right1     = (int *)mesh_memory.get_memory_ptr("map_xcell2face_right1");
   map_xcell2face_right2     = (int *)mesh_memory.get_memory_ptr("map_xcell2face_right2");
   xrecvIdx     = (int *)mesh_memory.get_memory_ptr("xrecvIdx");
   xrecvCIdx     = (int *)mesh_memory.get_memory_ptr("xrecvCIdx");
   xplusCell2Idx     = (int *)mesh_memory.get_memory_ptr("xplusCell2Idx");
   xminusCell2Idx     = (int *)mesh_memory.get_memory_ptr("xminusCell2Idx");
   xsendIdx1     = (int *)mesh_memory.get_memory_ptr("xsendIdx1");
   xsendIdx2     = (int *)mesh_memory.get_memory_ptr("xsendIdx2");
   yface_i     = (int *)mesh_memory.get_memory_ptr("yface_i");
   yface_j     = (int *)mesh_memory.get_memory_ptr("yface_j");
   yface_level     = (uchar_t *)mesh_memory.get_memory_ptr("yface_level");
   map_yface2cell_lower     = (int *)mesh_memory.get_memory_ptr("map_yface2cell_lower");
   map_yface2cell_upper     = (int *)mesh_memory.get_memory_ptr("map_yface2cell_upper");
   map_ycell2face_bot1     = (int *)mesh_memory.get_memory_ptr("map_ycell2face_bot1");
   map_ycell2face_bot2     = (int *)mesh_memory.get_memory_ptr("map_ycell2face_bot2");
   map_ycell2face_top1     = (int *)mesh_memory.get_memory_ptr("map_ycell2face_top1");
   map_ycell2face_top2     = (int *)mesh_memory.get_memory_ptr("map_ycell2face_top2");
   yrecvIdx     = (int *)mesh_memory.get_memory_ptr("yrecvIdx");
   yrecvCIdx     = (int *)mesh_memory.get_memory_ptr("yrecvCIdx");
   yplusCell2Idx     = (int *)mesh_memory.get_memory_ptr("yplusCell2Idx");
   yminusCell2Idx     = (int *)mesh_memory.get_memory_ptr("yminusCell2Idx");
   ysendIdx1     = (int *)mesh_memory.get_memory_ptr("ysendIdx1");
   ysendIdx2     = (int *)mesh_memory.get_memory_ptr("ysendIdx2");
}

void Mesh::resize_old_device_memory(size_t ncells)
{
#ifdef HAVE_OPENCL
   ezcl_device_memory_delete(dev_level);
   ezcl_device_memory_delete(dev_i);
   ezcl_device_memory_delete(dev_j);
   ezcl_device_memory_delete(dev_celltype);
   size_t mem_request = (int)((float)ncells*mem_factor);
   dev_level    = ezcl_malloc(NULL, const_cast<char *>("dev_level"),    &mem_request, sizeof(cl_uchar_t),  CL_MEM_READ_ONLY, 0);
   dev_i        = ezcl_malloc(NULL, const_cast<char *>("dev_i"),        &mem_request, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_j        = ezcl_malloc(NULL, const_cast<char *>("dev_j"),        &mem_request, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_celltype = ezcl_malloc(NULL, const_cast<char *>("dev_celltype"), &mem_request, sizeof(cl_char_t),  CL_MEM_READ_ONLY, 0);
#else
   // To get rid of compiler warning
   if (1 == 2) printf("DEBUG -- ncells is %lu\n",ncells);
#endif
}
void Mesh::print_object_info(void)
{
   printf(" ---- Mesh object info -----\n");
   printf("Dimensionality : %d\n",ndim);
   printf("Parallel info  : mype %d numpe %d noffset %d parallel %d\n",mype,numpe,noffset,parallel);
   printf("Sizes          : ncells %ld ncells_ghost %ld\n\n",ncells,ncells_ghost);
#ifdef HAVE_OPENCL
   int num_elements, elsize;

   num_elements = ezcl_get_device_mem_nelements(dev_celltype);
   elsize = ezcl_get_device_mem_elsize(dev_celltype);
   printf("dev_celltype     ptr : %p nelements %d elsize %d\n",dev_celltype,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_level);
   elsize = ezcl_get_device_mem_elsize(dev_level);
   printf("dev_level        ptr : %p nelements %d elsize %d\n",dev_level,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_i);
   elsize = ezcl_get_device_mem_elsize(dev_i);
   printf("dev_i            ptr : %p nelements %d elsize %d\n",dev_i,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_j);
   elsize = ezcl_get_device_mem_elsize(dev_j);
   printf("dev_j            ptr : %p nelements %d elsize %d\n",dev_j,num_elements,elsize);

   num_elements = ezcl_get_device_mem_nelements(dev_nlft);
   elsize = ezcl_get_device_mem_elsize(dev_nlft);
   printf("dev_nlft         ptr : %p nelements %d elsize %d\n",dev_nlft,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_nrht);
   elsize = ezcl_get_device_mem_elsize(dev_nrht);
   printf("dev_nrht         ptr : %p nelements %d elsize %d\n",dev_nrht,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_nbot);
   elsize = ezcl_get_device_mem_elsize(dev_nbot);
   printf("dev_nbot         ptr : %p nelements %d elsize %d\n",dev_nbot,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_ntop);
   elsize = ezcl_get_device_mem_elsize(dev_ntop);
   printf("dev_ntop         ptr : %p nelements %d elsize %d\n",dev_ntop,num_elements,elsize);
#endif
   printf("vector celltype  ptr : %p nelements %ld elsize %ld\n",&celltype[0],mesh_memory.get_memory_size(celltype),sizeof(celltype[0])); 
   printf("vector level     ptr : %p nelements %ld elsize %ld\n",&level[0],   mesh_memory.get_memory_size(level),   sizeof(level[0])); 
   printf("vector i         ptr : %p nelements %ld elsize %ld\n",&i[0],       mesh_memory.get_memory_size(i),       sizeof(i[0])); 
   printf("vector j         ptr : %p nelements %ld elsize %ld\n",&j[0],       mesh_memory.get_memory_size(j),       sizeof(j[0])); 

   printf("vector nlft      ptr : %p nelements %ld elsize %ld\n",&nlft[0],    mesh_memory.get_memory_size(nlft),    sizeof(nlft[0])); 
   printf("vector nrht      ptr : %p nelements %ld elsize %ld\n",&nrht[0],    mesh_memory.get_memory_size(nrht),    sizeof(nrht[0])); 
   printf("vector nbot      ptr : %p nelements %ld elsize %ld\n",&nbot[0],    mesh_memory.get_memory_size(nbot),    sizeof(nbot[0])); 
   printf("vector ntop      ptr : %p nelements %ld elsize %ld\n",&ntop[0],    mesh_memory.get_memory_size(ntop),    sizeof(ntop[0])); 
}


void Mesh::set_refinement_order(int order[4], int ic, int ifirst, int ilast, int jfirst, int jlast,
                                uchar_t level_first, uchar_t level_last, int *i_old, int *j_old, uchar_t *level_old)
{
            if (localStencil) {
               //  Store the coordinates of the cells before and after this one on
               //  the space-filling curve index.

#ifdef __OLD_STENCIL__
               spatial_t  nx[3],  //  x-coordinates of cells.
                          ny[3];  //  y-coordinates of cells.
               if (ic != 0) {
                  nx[0] = lev_deltax[level_old[ic-1]] * (spatial_t)i[ic-1];
                  ny[0] = lev_deltay[level_old[ic-1]] * (spatial_t)j[ic-1];
               } else {
                  nx[0] = lev_deltax[level_first] * (spatial_t)ifirst;
                  ny[0] = lev_deltay[level_first] * (spatial_t)jfirst;
               }
               nx[1] = lev_deltax[level_old[ic  ]] * (spatial_t)i[ic  ];
               ny[1] = lev_deltay[level_old[ic  ]] * (spatial_t)j[ic  ];
               if (ic != ncells-1) {
                  nx[2] = lev_deltax[level_old[ic+1]] * (spatial_t)i[ic+1];
                  ny[2] = lev_deltay[level_old[ic+1]] * (spatial_t)j[ic+1];
               } else {
                  nx[2] = lev_deltax[level_last] * (spatial_t)ilast;
                  ny[2] = lev_deltay[level_last] * (spatial_t)jlast;
               }

               //  Figure out relative orientation of the neighboring cells.  We are
               //  are aided in this because the Hilbert curve only has six possible
               //  ways across the cell:  four Ls and two straight lines.  Then
               //  refine the cell according to the relative orientation and order
               //  according to the four-point Hilbert stencil.
               if      (nx[0] < nx[1] and ny[2] < ny[1])   //  southwest L, forward order
               {  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE; }
               else if (nx[2] < nx[1] and ny[0] < ny[1])   //  southwest L, reverse order
               {  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW; }
               else if (nx[0] > nx[1] and ny[2] < ny[1])   //  southeast L, forward order
               {  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW; }
               else if (nx[2] > nx[1] and ny[0] < ny[1])   //  southeast L, reverse order
               {  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE; }
               else if (nx[0] > nx[1] and ny[2] > ny[1])   //  northeast L, forward order
               {  order[0] = SE; order[1] = SW; order[2] = NW; order[3] = NE; }
               else if (nx[2] > nx[1] and ny[0] > ny[1])   //  northeast L, reverse order
               {  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE; }
               else if (nx[0] < nx[1] and ny[2] > ny[1])   //  northwest L, forward order
               {  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW; }
               else if (nx[2] < nx[1] and ny[0] > ny[1])   //  northwest L, reverse order
               {  order[0] = NW; order[1] = NE; order[2] = SE; order[3] = SW; }
               else if (nx[0] > nx[1] and nx[1] > nx[2])   //  straight horizontal, forward order
               {  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW; }
               else if (nx[0] < nx[1] and nx[1] < nx[2])   //  straight horizontal, reverse order
               {  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE; }
               else if (ny[0] > ny[1] and ny[1] > ny[2])   //  straight vertical, forward order
               {  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE; }
               else if (ny[0] < ny[1] and ny[1] < ny[2])   //  straight vertical, reverse order
               {  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW; }
               else                                        //  other, default to z-order
               {  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE; }
#endif

#ifdef __NEW_STENCIL__
               int ir[3],   // First i index at finest level of the mesh
                   jr[3];   // First j index at finest level of the mesh
               // Cell's Radius at the Finest level of the mesh

               int crf = IPOW2(levmx-level_old[ic]);

               if (ic != 0) {
                  ir[0] = i_old[ic - 1] * IPOW2(levmx-level_old[ic - 1]);
                  jr[0] = j_old[ic - 1] * IPOW2(levmx-level_old[ic - 1]);
               } else {
                  //printf("%d cell %d is a first\n",mype,ic);
                  ir[0] = ifirst * IPOW2(levmx-level_first);
                  jr[0] = jfirst * IPOW2(levmx-level_first);
               }
               ir[1] = i_old[ic    ] * IPOW2(levmx-level_old[ic    ]);
               jr[1] = j_old[ic    ] * IPOW2(levmx-level_old[ic    ]);
               if (ic != (int)ncells-1) {
                  ir[2] = i_old[ic + 1] * IPOW2(levmx-level_old[ic + 1]);
                  jr[2] = j_old[ic + 1] * IPOW2(levmx-level_old[ic + 1]);
               } else {
                  //printf("%d cell %d is a last\n",mype,ic);
                  ir[2] = ilast * IPOW2(levmx-level_last);
                  jr[2] = jlast * IPOW2(levmx-level_last);
               }
               //if (parallel) fprintf(fp,"%d: DEBUG rezone top boundary -- ic %d global %d noffset %d nc %d i %d j %d level %d\n",mype,ic,ic+noffset,noffset,nc,i[nc],j[nc],level[nc]);

               int dir_in  = ir[1] - ir[0];
               int dir_out = ir[1] - ir[2];
               int djr_in  = jr[1] - jr[0];
               int djr_out = jr[1] - jr[2];

               char  in_direction = 'X';
               char out_direction = 'X';

               // Left In
               if( (djr_in == 0 && (dir_in == crf*HALF || dir_in == crf || dir_in == crf*TWO)) || (djr_in == -crf*HALF && dir_in == crf*HALF) || (djr_in == crf && dir_in == crf*TWO) ) {
                  in_direction = 'L';
               }
               // Bottom In
               else if( (dir_in == 0 && (djr_in == crf*HALF || djr_in == crf || djr_in == crf*TWO)) || (dir_in == -crf*HALF && djr_in == crf*HALF) || (dir_in == crf && djr_in == crf*TWO) ) {
                  in_direction = 'B';
               }
               // Right In
               else if( (dir_in == -crf && (djr_in == -crf*HALF || djr_in == 0 || (djr_in == crf && level_old[ic-1] < level_old[ic]))) ) {
                  in_direction = 'R';
               }
               // Top In
               else if( (djr_in == -crf && (dir_in == -crf*HALF || dir_in == 0 || (dir_in == crf && level_old[ic-1] < level_old[ic]))) ) {
                  in_direction = 'T';
               }
               // Further from the left
               else if( dir_in > 0 && djr_in == 0 ) {
                  in_direction = 'L';
               }
               // Further from the right
               else if( dir_in < 0 && djr_in == 0 ) {
                  in_direction = 'R';
               }
               // Further from the bottom
               else if( djr_in > 0 && dir_in == 0 ) {
                  in_direction = 'B';
               }
               // Further from the top
               else if( djr_in < 0 && dir_in == 0 ) {
                  in_direction = 'T';
               }
               // SW in; 'M'
               else if( dir_in > 0 && djr_in > 0) {
                  in_direction = 'M';
               }
               // NW in; 'W'
               else if( dir_in > 0 && djr_in < 0) {
                  in_direction = 'W';
               }
               // SE in; 'F'
               else if( dir_in < 0 && djr_in > 0) {
                  in_direction = 'F';
               }
               // NE in; 'E'
               else if( dir_in < 0 && djr_in < 0) {
                  in_direction = 'E';
               }

   
               // Left Out
               if( (djr_out == 0 && (dir_out == crf*HALF || dir_out == crf || dir_out == crf*TWO)) || (djr_out == -crf*HALF && dir_out == crf*HALF) || (djr_out == crf && dir_out == crf*TWO) ) {
                  out_direction = 'L';
               }
               // Bottom Out
               else if( (dir_out == 0 && (djr_out == crf*HALF || djr_out == crf || djr_out == crf*TWO)) || (dir_out == -crf*HALF && djr_out == crf*HALF) || (dir_out == crf && djr_out == crf*TWO) ) {
                  out_direction = 'B';
               }
               // Right Out
               else if( (dir_out == -crf && (djr_out == -crf*HALF || djr_out == 0 || (djr_out == crf && level_old[ic+1] < level_old[ic]))) ) {
                  out_direction = 'R';
               }
               // Top Out
               else if( (djr_out == -crf && (dir_out == -crf*HALF || dir_out == 0 || (dir_out == crf && level_old[ic+1] < level_old[ic]))) ) {
                  out_direction = 'T';
               }
               // Further from the left
               else if( dir_out > 0 && djr_out == 0 ) {
                  out_direction = 'L';
               }
               // Further from the right
               else if( dir_out < 0 && djr_out == 0 ) {
                  out_direction = 'R';
               }
               // Further from the bottom
               else if( djr_out > 0 && dir_out == 0 ) {
                  out_direction = 'B';
               }
               // Further from the top
               else if( djr_out < 0 && dir_out == 0 ) {
                  out_direction = 'T';
               }
               // SW out; 'M'
               else if( dir_out > 0 && djr_out > 0) {
                  out_direction = 'M';
               }
               // NW out; 'W'
               else if( dir_out > 0 && djr_out < 0) {
                  out_direction = 'W';
               }
               // SE out; 'F'
               else if( dir_out < 0 && djr_out > 0) {
                  out_direction = 'F';
               }
               // NE out; 'E'
               else if( dir_out < 0 && djr_out < 0) {
                  out_direction = 'E';
               }

               // Set the Stencil
               if(in_direction == 'L' && (out_direction == 'B' || out_direction == 'R' || out_direction == 'F')) {
                  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE;
               }
               else if(in_direction == 'L' && (out_direction == 'T' || out_direction == 'W' )) {
                  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW;
               }
               else if(in_direction == 'L' && out_direction == 'M') {
                  order[0] = NW; order[1] = NE; order[2] = SE; order[3] = SW;
               }
               else if(in_direction == 'L' && out_direction == 'E') {
                  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE;
               }

               else if(in_direction == 'B' && (out_direction == 'R' || out_direction == 'F' )) {
                  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE;
               }
               else if(in_direction == 'B' && (out_direction == 'L' || out_direction == 'T' || out_direction == 'W' )) {
                  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW;
               }
               else if(in_direction == 'B' && out_direction == 'M') {
                  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW;
               }
               else if(in_direction == 'B' && out_direction == 'E') {
                  order[0] = SW; order[1] = NW; order[2] = SE; order[3] = NE;
               }
               
               else if(in_direction == 'R' && (out_direction == 'T' || out_direction == 'L' || out_direction == 'W' )) {
                  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW;
               }
               else if(in_direction == 'R' && (out_direction == 'B' || out_direction == 'F' )) {
                  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE;
               }
               else if(in_direction == 'R' && out_direction == 'M') {
                  order[0] = NE; order[1] = NW; order[2] = SE; order[3] = SW;
               }
               else if(in_direction == 'R' && out_direction == 'E') {
                  order[0] = SE; order[1] = SW; order[2] = NW; order[3] = NE;
               }

               else if(in_direction == 'T' && (out_direction == 'L' || out_direction == 'W' )) {
                  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW;
               }
               else if(in_direction == 'T' && (out_direction == 'R' || out_direction == 'B' || out_direction == 'F' )) {
                  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE;
               }
               else if(in_direction == 'T' && out_direction == 'M') {
                  order[0] = NE; order[1] = SE; order[2] = NW; order[3] = SW;
               }
               else if(in_direction == 'T' && out_direction == 'E') {
                  order[0] = NW; order[1] = SW; order[2] = SE; order[3] = NE;
               }

               else if(in_direction == 'M' && (out_direction == 'L' || out_direction == 'W' || out_direction == 'T') ) {
                  order[0] = SW; order[1] = SE; order[2] = NE; order[3] = NW;
               }
               else if(in_direction == 'M' && (out_direction == 'R' || out_direction == 'F' || out_direction == 'B') ) {
                  order[0] = SW; order[1] = NW; order[2] = NE; order[3] = SE;
               }
               else if(in_direction == 'M' && out_direction == 'E') {
                  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE;
               }
 
               else if(in_direction == 'W' && (out_direction == 'L' || out_direction == 'M' || out_direction == 'B') ) {
                  order[0] = NW; order[1] = NE; order[2] = SE; order[3] = SW;
               }
               else if(in_direction == 'W' && (out_direction == 'R' || out_direction == 'E' || out_direction == 'T') ) {
                  order[0] = NW; order[1] = SW; order[2] = SE; order[3] = NE;
               }
               else if(in_direction == 'W' && out_direction == 'F') {
                  order[0] = NW; order[1] = NE; order[2] = SW; order[3] = SE;
               }

               else if(in_direction == 'F' && (out_direction == 'L' || out_direction == 'M' || out_direction == 'B') ) {
                  order[0] = SE; order[1] = NE; order[2] = NW; order[3] = SW;
               }
               else if(in_direction == 'F' && (out_direction == 'R' || out_direction == 'E' || out_direction == 'T') ) {
                  order[0] = SE; order[1] = SW; order[2] = NW; order[3] = NE;
               }
               else if(in_direction == 'F' && out_direction == 'W') {
                  order[0] = SE; order[1] = NE; order[2] = SW; order[3] = NW;
               }

               else if(in_direction == 'E' && (out_direction == 'L' || out_direction == 'W' || out_direction == 'T') ) {
                  order[0] = NE; order[1] = SE; order[2] = SW; order[3] = NW;
               }
               else if(in_direction == 'E' && (out_direction == 'R' || out_direction == 'F' || out_direction == 'B') ) {
                  order[0] = NE; order[1] = NW; order[2] = SW; order[3] = SE;
               }
               else if(in_direction == 'E' && out_direction == 'M') {
                  order[0] = NE; order[1] = SE; order[2] = NW; order[3] = SW;
               }

               else { // Default to a knot 
                  order[0] = NW; order[1] = SE; order[2] = SW; order[3] = NE;
                  if (do_stencil_warning) {
                     printf("Nonlocal case for the stencil.\n");
                  }
               }
               //  Determine the relative orientation of the neighboring cells.
               //  There are 12 possible ways across the cell: 4 Ls and 2 straight
               //  lines, each with 2 directions of traversal.
               //  Then the cell is refined and ordered according to the relative
               //  orientation and four-point Hilbert stencil.

               // XXX NOTE that the four-point stencil varies depending upon
               // the starting and ending point of the global Hilbert curve.
               // The stencil applied here assumes the start at (0,0) and the end
               // at (0,y_max). XXX WRONG
#endif                 

            }  //  End local stencil version
            else //  Use Z-ordering for the curve.
            {  order[0] = SW; order[1] = SE; order[2] = NW; order[3] = NE; }
            
}

void Mesh::calc_face_list(void)
{
/*    
   xface_i.clear();
   xface_j.clear();
   xface_level.clear();

   ixmin_level.clear();
   ixmax_level.clear();
   jxmin_level.clear();
   jxmax_level.clear();
   ixmin_level.resize(levmx+1,  9999999);
   ixmax_level.resize(levmx+1, -9999999);
   jxmin_level.resize(levmx+1,  9999999);
   jxmax_level.resize(levmx+1, -9999999);

   ixadjust.clear();
   ixadjust.resize(levmx+1);
   jxadjust.clear();
   jxadjust.resize(levmx+1);

   int iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nr = nrht[nz];
      if (nr == nz) continue;

      int ifactor = 1;
      if (level[nr] < level[nz]) ifactor = 2;

      // Have right face
      //printf("DEBUG xface -- iface %d lower nz %d upper nr %d\n",iface,nz,nr);
      xface_level.push_back(MAX(level[nz],level[nr]));
      xface_i.push_back(i[nr]*ifactor);
      if (level[nr] < level[nz] && is_upper(j[nz]) ) {
         xface_j.push_back(j[nr]*ifactor+1);
      } else {
         xface_j.push_back(j[nr]*ifactor);
      }

      iface++;

      if (level[nr] > level[nz] && is_lower(j[nr]) ){
         int ntr = ntop[nr];
         if (ntr != nr) {
            //printf("DEBUG xface -- iface %d lower nz %d upper ntr %d\n",iface,nz,ntr);
            xface_level.push_back(MAX(level[nz],level[ntr]));
            xface_i.push_back(i[ntr]*ifactor);
            xface_j.push_back(j[ntr]*ifactor);

            iface++;
         }
      }
   }
   nxface=iface;

   yface_i.clear();
   yface_j.clear();
   yface_level.clear();

   iymin_level.clear();
   iymax_level.clear();
   jymin_level.clear();
   jymax_level.clear();
   iymin_level.resize(levmx+1,  9999999);
   iymax_level.resize(levmx+1, -9999999);
   jymin_level.resize(levmx+1,  9999999);
   jymax_level.resize(levmx+1, -9999999);

   iyadjust.clear();
   iyadjust.resize(levmx+1);
   jyadjust.clear();
   jyadjust.resize(levmx+1);

   iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nt = ntop[nz];
      if (nt == nz) continue;

      int ifactor = 1;
      if (level[nt] < level[nz]) ifactor = 2;

      // Have top face
      //printf("DEBUG yface -- iface %d lower nz %d upper nt %d\n",iface,nz,nt);
      yface_level.push_back(MAX(level[nz],level[nt]));
      yface_j.push_back(j[nt]*ifactor);
      if (level[nt] < level[nz] && is_upper(i[nz]) ) {
         yface_i.push_back(i[nt]*ifactor+1);
      } else{
         yface_i.push_back(i[nt]*ifactor);
      }

      iface++;
      if (level[nt] > level[nz] && is_lower(i[nt]) ){
         int nrt = nrht[nt];
         if (nrt != nt) {
            //printf("DEBUG yface -- iface %d lower nz %d upper nrt %d\n",iface,nz,nrt);
            yface_level.push_back(MAX(level[nz],level[nrt]));
            yface_j.push_back(j[nrt]*ifactor);
            yface_i.push_back(i[nrt]*ifactor);

            iface++;
         }

      }
   }
   nyface=iface;

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];

      int fi = xface_i[iface];
      if (fi < ixmin_level[fl]) ixmin_level[fl] = fi;
      if (fi > ixmax_level[fl]) ixmax_level[fl] = fi;

      int fj = xface_j[iface];
      if (fj < jxmin_level[fl]) jxmin_level[fl] = fj;
      if (fj > jxmax_level[fl]) jxmax_level[fl] = fj;
   }

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];
      if (ixmax_level[fl] < ixmin_level[fl]) continue;

      xface_i[iface] -= ixmin_level[fl];
      xface_j[iface] -= jxmin_level[fl];
   }

   for (int fl = 0; fl <= levmx; fl++){
      ixadjust[fl] = ixmin_level[fl];
      jxadjust[fl] = jxmin_level[fl];
      ixmax_level[fl] -= ixmin_level[fl];;
      jxmax_level[fl] -= jxmin_level[fl];
      ixmin_level[fl] = 0;
      jxmin_level[fl] = 0;
   }

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];

      int fi = yface_i[iface];
      if (fi < iymin_level[fl]) iymin_level[fl] = fi;
      if (fi > iymax_level[fl]) iymax_level[fl] = fi;

      int fj = yface_j[iface];
      if (fj < jymin_level[fl]) jymin_level[fl] = fj;
      if (fj > jymax_level[fl]) jymax_level[fl] = fj;
   }

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];
      if (iymax_level[fl] < iymin_level[fl]) continue;

      yface_i[iface] -= iymin_level[fl];
      yface_j[iface] -= jymin_level[fl];
   }

   for (int fl = 0; fl <= levmx; fl++){
      iyadjust[fl] = iymin_level[fl];
      jyadjust[fl] = jymin_level[fl];
      iymax_level[fl] -= iymin_level[fl];;
      jymax_level[fl] -= jymin_level[fl];
      iymin_level[fl] = 0;
      jymin_level[fl] = 0;
   }
*/
}

void Mesh::calc_face_list_wmap(void)
{
/*
   map_xface2cell_lower.clear();
   map_xface2cell_upper.clear();

   xface_i.clear();
   xface_j.clear();
   xface_level.clear();

   ixmin_level.clear();
   ixmax_level.clear();
   jxmin_level.clear();
   jxmax_level.clear();
   ixmin_level.resize(levmx+1,  9999999);
   ixmax_level.resize(levmx+1, -9999999);
   jxmin_level.resize(levmx+1,  9999999);
   jxmax_level.resize(levmx+1, -9999999);

   ixadjust.clear();
   ixadjust.resize(levmx+1);
   jxadjust.clear();
   jxadjust.resize(levmx+1);

   int iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nr = nrht[nz];
      if (nr == nz) continue;

      int ifactor = 1;
      if (level[nr] < level[nz]) ifactor = 2;

      // Have right face
      map_xface2cell_lower.push_back(nz);
      map_xface2cell_upper.push_back(nr);
      xface_level.push_back(MAX(level[nz],level[nr]));
      xface_i.push_back(i[nr]*ifactor);
      if (level[nr] < level[nz] && is_upper(j[nz]) ) {
         xface_j.push_back(j[nr]*ifactor+1);
      } else {
         xface_j.push_back(j[nr]*ifactor);
      }

      iface++;

      if (level[nr] > level[nz] && is_lower(j[nr]) ){
         int ntr = ntop[nr];
         if (ntr != nr) {
            map_xface2cell_lower.push_back(nz);
            map_xface2cell_upper.push_back(ntr);
            xface_level.push_back(MAX(level[nz],level[ntr]));
            xface_i.push_back(i[ntr]*ifactor);
            xface_j.push_back(j[ntr]*ifactor);

            iface++;
         }
      }
   }
   nxface=iface;

   map_yface2cell_lower.clear();
   map_yface2cell_upper.clear();

   yface_i.clear();
   yface_j.clear();
   yface_level.clear();

   iymin_level.clear();
   iymax_level.clear();
   jymin_level.clear();
   jymax_level.clear();
   iymin_level.resize(levmx+1,  9999999);
   iymax_level.resize(levmx+1, -9999999);
   jymin_level.resize(levmx+1,  9999999);
   jymax_level.resize(levmx+1, -9999999);

   iyadjust.clear();
   iyadjust.resize(levmx+1);
   jyadjust.clear();
   jyadjust.resize(levmx+1);

   iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nt = ntop[nz];
      if (nt == nz) continue;

      int ifactor = 1;
      if (level[nt] < level[nz]) ifactor = 2;

      // Have top face
      // printf("DEBUG -- iface %d lower nz %d upper nr %d\n",iface,nz,nt);
      map_yface2cell_lower.push_back(nz);
      map_yface2cell_upper.push_back(nt);
      yface_level.push_back(MAX(level[nz],level[nt]));
      yface_j.push_back(j[nt]*ifactor);
      if (level[nt] < level[nz] && is_upper(i[nz]) ) {
         yface_i.push_back(i[nt]*ifactor+1);
      } else{
         yface_i.push_back(i[nt]*ifactor);
      }

      iface++;
      if (level[nt] > level[nz] && is_lower(i[nt]) ){
         int nrt = nrht[nt];
         if (nrt != nt) {
            map_yface2cell_lower.push_back(nz);
            map_yface2cell_upper.push_back(nrt);
            yface_level.push_back(MAX(level[nz],level[nrt]));
            yface_j.push_back(j[nrt]*ifactor);
            yface_i.push_back(i[nrt]*ifactor);

            iface++;
         }
      }
   }
   nyface=iface;

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];

      int fi = xface_i[iface];
      if (fi < ixmin_level[fl]) ixmin_level[fl] = fi;
      if (fi > ixmax_level[fl]) ixmax_level[fl] = fi;

      int fj = xface_j[iface];
      if (fj < jxmin_level[fl]) jxmin_level[fl] = fj;
      if (fj > jxmax_level[fl]) jxmax_level[fl] = fj;
   }

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];
      if (ixmax_level[fl] < ixmin_level[fl]) continue;

      xface_i[iface] -= ixmin_level[fl];
      xface_j[iface] -= jxmin_level[fl];
   }

   for (int fl = 0; fl <= levmx; fl++){
      ixadjust[fl] = ixmin_level[fl];
      jxadjust[fl] = jxmin_level[fl];
      ixmax_level[fl] -= ixmin_level[fl];;
      jxmax_level[fl] -= jxmin_level[fl];
      ixmin_level[fl] = 0;
      jxmin_level[fl] = 0;
   }

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];

      int fi = yface_i[iface];
      if (fi < iymin_level[fl]) iymin_level[fl] = fi;
      if (fi > iymax_level[fl]) iymax_level[fl] = fi;

      int fj = yface_j[iface];
      if (fj < jymin_level[fl]) jymin_level[fl] = fj;
      if (fj > jymax_level[fl]) jymax_level[fl] = fj;
   }

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];
      if (iymax_level[fl] < iymin_level[fl]) continue;

      yface_i[iface] -= iymin_level[fl];
      yface_j[iface] -= jymin_level[fl];
   }

   for (int fl = 0; fl <= levmx; fl++){
      iyadjust[fl] = iymin_level[fl];
      jyadjust[fl] = jymin_level[fl];
      iymax_level[fl] -= iymin_level[fl];;
      jymax_level[fl] -= jymin_level[fl];
      iymin_level[fl] = 0;
      jymin_level[fl] = 0;
   }
*/
}

void Mesh::quickInterpolate(int side_main, int side_sec, int cncell, double* mem_ptr_double,
                            real_t d_lo, real_t d_hi, int flag, real_t *fineavg, real_t *coarseavg)
{
    //real_t state_mainmain, state_mainsec, state_secmain, state_secsec, state_sideavg = 0; //vars for 2 cells over
    real_t state_main, state_sec, state_coarse, f2c_avg, c2f_avg;

    state_main = mem_ptr_double[side_main];
    state_sec = mem_ptr_double[side_sec];
    state_coarse = mem_ptr_double[cncell];
    if (flag == 0) { //only working with a single cell, either coarse or refined (dampening)
        c2f_avg = 2 * (d_lo * state_main + d_hi * state_coarse) / (d_lo + d_hi) - state_coarse;
        //mem_ptr_double[pcellIdx] = c2f_avg;
        *fineavg = c2f_avg;
    }
    else if (flag == 1) { //only adding coarser refined cells if we know there are 2 neighbors (dampening)
        f2c_avg = (2 * (d_lo * state_main + d_hi * state_coarse) / (d_lo + d_hi) - state_main) + (2 * (d_lo * state_sec + d_hi * state_coarse) / (d_lo + d_hi) - state_sec);
        //mem_ptr_double[pcellIdx+2] = f2c_avg / 2;
        *coarseavg = f2c_avg / 2;
    }
    else { //adding both (base neighbors)
        c2f_avg = 2 * (d_lo * state_main + d_hi * state_coarse) / (d_lo + d_hi) - state_coarse;
        //mem_ptr_double[pcellIdx] = c2f_avg;
        *fineavg = c2f_avg;
        f2c_avg = (2 * (d_lo * state_main + d_hi * state_coarse) / (d_lo + d_hi) - state_main) + (2 * (d_lo * state_sec + d_hi * state_coarse) / (d_lo + d_hi) - state_sec);
        //mem_ptr_double[pcellIdx+2] = f2c_avg / 2;
        *coarseavg = f2c_avg / 2;
    }

}

double Mesh::xFakeFlux(double* locH, double* locU, double* locV, int idx, int caseNum) {
    double flux;

    switch (caseNum) {
        case 0: //HXFLUX
            flux = locU[idx];
            break;
        case 1: //UXFLUX
            flux = (locU[idx] * locU[idx]) / locH[idx] + 4.90 * (locH[idx] * locH[idx]);
            break;
        case 2: //UVFLUX
            flux = locU[idx] * locV[idx] / locH[idx];
            break;
        default:
            flux = 0.0;
            break;
    }

    return flux;
}

double Mesh::yFakeFlux(double* locH, double* locU, double* locV, int idx, int caseNum) {
    double flux;

    switch (caseNum) {
        case 0: //HYFLUX
            flux = locV[idx];
            break;
        case 1: //VUFLUX
            flux = locV[idx] * locU[idx] / locH[idx];
            break;
        case 2: //VYFLUX
            flux = (locV[idx] * locV[idx]) / locH[idx] + 4.90 * (locH[idx] * locH[idx]);
            break;
        default:
            flux = 0.0;
            break;
    }

    return flux;
}


void Mesh::interpolate(int scheme, int index, int cell_lower, int cell_upper, double deltaT, MallocPlus &state_memory_old)
{
   malloc_plus_memory_entry *memory_item;
   int cell_course, cell_fine, cell_bottom, cell_top, cell_left, cell_right;
   real_t state_bot, state_top, state_lft, state_rht, state_avg;
   real_t state_botbot, state_bottop, state_topbot, state_toptop; 
   real_t state_lftlft, state_lftrht, state_rhtlft, state_rhtrht;
   real_t state_sideavg = ZERO; 
   //bool fake_flux = false;
   bool five_point = true;

   switch(scheme){

      case 0: // fine cell, x-direction, right cell more refined 
      // loop through state arrays to update phantom cell state values
      cell_course = cell_lower; // course neighbor to the face
      cell_fine =  cell_upper; // fine neighbor to the face
      cell_bottom = nrht[cell_course]; // bottom of the two fine neighbors of the course
      cell_top = ntop[cell_bottom]; // top of the two fince neighbots of the course
 
      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {
 
         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;
 
         if (memory_item->mem_elsize == 8) {
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
            mem_ptr_double[index] = mem_ptr_double[cell_course]; // left
            if (five_point){
               if (level[nlft[cell_course]] <= level[cell_course] || cell_fine == cell_bottom){ // 2 cells over is same or lesser refine.
                  mem_ptr_double[index+1] = mem_ptr_double[nlft[cell_course]]; // left left
               }
               else {
                  mem_ptr_double[index+1] = mem_ptr_double[ntop[nlft[cell_course]]]; // left left
               }
            }
         } else if (memory_item->mem_elsize == 4) {
            float *mem_ptr_float = (float *)memory_item->mem_ptr;
            mem_ptr_float[index] = mem_ptr_float[cell_course]; // left
            if (five_point){
               if (level[nlft[cell_course]] <= level[cell_course] || cell_fine == cell_bottom){ // 2 cells over is same or lesser refine.
                  mem_ptr_float[index+1] = mem_ptr_float[nlft[cell_course]]; // left left
               }
               else {
                  mem_ptr_float[index+1] = mem_ptr_float[ntop[nlft[cell_course]]]; // left left
               }
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            short *mem_ptr_short = (short *)memory_item->mem_ptr;
            mem_ptr_short[index] = mem_ptr_short[cell_course]; // left
            if (five_point){
               if (level[nlft[cell_course]] <= level[cell_course] || cell_fine == cell_bottom){ // 2 cells over is same or lesser refine.
                  mem_ptr_short[index+1] = mem_ptr_short[nlft[cell_course]]; // left left
               }
               else {
                  mem_ptr_short[index+1] = mem_ptr_short[ntop[nlft[cell_course]]]; // left left
               }
            }
#endif
         }
      } 

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lft\n",index,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lftlft\n",index+1,cell_lower,cell_upper);
      break;
 
      case 1: // fine cell, x-direction, left cell more refined
      // loop through state arrays to update phantom cell state values
      cell_course = cell_upper; // course neighbor to the face
      cell_fine =  cell_lower; // fine neighbor to the face
      cell_bottom = nlft[cell_course]; // bottom of the two fine neighbors of the course
      cell_top = ntop[cell_bottom]; // top of the two fince neighbots of the course
 
      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {
 
         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;
 
         if (memory_item->mem_elsize == 8) {
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
            mem_ptr_double[index+2] = mem_ptr_double[cell_course]; // right
            if (five_point){
               if (level[nrht[cell_course]] <= level[cell_course] || cell_fine == cell_bottom){ // 2 cells over is same or lesser refine.
                  mem_ptr_double[index+3] = mem_ptr_double[nrht[cell_course]]; // right right
               }
               else {
                  mem_ptr_double[index+3] = mem_ptr_double[ntop[nrht[cell_course]]]; // right right
               }
            }
         } else if (memory_item->mem_elsize == 4) {
            float *mem_ptr_float = (float *)memory_item->mem_ptr;
            mem_ptr_float[index+2] = mem_ptr_float[cell_course]; // right
            if (five_point){
               if (level[nrht[cell_course]] <= level[cell_course] || cell_fine == cell_bottom){ // 2 cells over is same or lesser refine.
                  mem_ptr_float[index+3] = mem_ptr_float[nrht[cell_course]]; // right right
               }
               else {
                  mem_ptr_float[index+3] = mem_ptr_float[ntop[nrht[cell_course]]]; // right right
               }
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            short *mem_ptr_short = (short *)memory_item->mem_ptr;
            mem_ptr_short[index+2] = mem_ptr_short[cell_course]; // right
            if (five_point){
               if (level[nrht[cell_course]] <= level[cell_course] || cell_fine == cell_bottom){ // 2 cells over is same or lesser refine.
                  mem_ptr_short[index+3] = mem_ptr_short[nrht[cell_course]]; // right right
               }
               else {
                  mem_ptr_short[index+3] = mem_ptr_short[ntop[nrht[cell_course]]]; // right right
               }
            }
#endif
         }
      }

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: rht\n",index+2,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: rhtrht\n",index+3,cell_lower,cell_upper);
      break;
 
      case 2: // fine cell, y-direction, top cell more refined
      // loop through state arrays to update phantom cell state values
      cell_course = cell_lower; // course neighbor to the face
      cell_fine = cell_upper; // fine neighbor to the face
      cell_left = ntop[cell_course]; // bottom of the two fine neighbors of the course
      cell_right = nrht[cell_fine]; // top of the two fince neighbots of the course
 
      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {
 
         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;
 
         if (memory_item->mem_elsize == 8) {
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
            mem_ptr_double[index] = mem_ptr_double[cell_course]; // bottom
            if(five_point){
               if (level[nbot[cell_course]] <= level[cell_course] || cell_fine == cell_left){ // 2 cells over is same or lesser refine.
                  mem_ptr_double[index+1] = mem_ptr_double[nbot[cell_course]]; // bottom bottom
               }
               else {
                  mem_ptr_double[index+1] = mem_ptr_double[nrht[nbot[cell_course]]]; // bottom bottom
               }
            }
         } else if (memory_item->mem_elsize == 4) {
            float *mem_ptr_float = (float *)memory_item->mem_ptr;
            mem_ptr_float[index] = mem_ptr_float[cell_course]; // bottom
            if(five_point){
               if (level[nbot[cell_course]] <= level[cell_course] || cell_fine == cell_left){ // 2 cells over is same or lesser refine.
                  mem_ptr_float[index+1] = mem_ptr_float[nbot[cell_course]]; // bottom bottom
               }
               else {
                  mem_ptr_float[index+1] = mem_ptr_float[nrht[nbot[cell_course]]]; // bottom bottom
               }
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            short *mem_ptr_short = (short *)memory_item->mem_ptr;
            mem_ptr_short[index] = mem_ptr_short[cell_course]; // bottom
            if(five_point){
               if (level[nbot[cell_course]] <= level[cell_course] || cell_fine == cell_left){ // 2 cells over is same or lesser refine.
                  mem_ptr_short[index+1] = mem_ptr_short[nbot[cell_course]]; // bottom bottom
               }
               else {
                  mem_ptr_short[index+1] = mem_ptr_short[nrht[nbot[cell_course]]]; // bottom bottom
               }
            }
#endif
         }
      }

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: bot\n",index,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: botbot\n",index+1,cell_lower,cell_upper);
      break;
 
      case 3: // fine cell, y-direction, bottom cell more refined
      // loop through state arrays to update phantom cell state values
      cell_course = cell_upper; // course neighbor to the face
      cell_fine = cell_lower; // fine neighbor to the face
      cell_left = nbot[cell_course]; // bottom of the two fine neighbors of the course
      cell_right = nrht[cell_fine]; // top of the two fince neighbots of the course

      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {

         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;

         if (memory_item->mem_elsize == 8) {
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
            mem_ptr_double[index+2] = mem_ptr_double[cell_course]; // top
            if(five_point){
               if (level[ntop[cell_course]] <= level[cell_course] || cell_fine == cell_left){ // 2 cells over is same or lesser refine.
                  mem_ptr_double[index+3] = mem_ptr_double[ntop[cell_course]]; // top top
               }
               else {
                  mem_ptr_double[index+3] = mem_ptr_double[nrht[ntop[cell_course]]]; // top top
               }
            }
         } else if (memory_item->mem_elsize == 4) {
            float *mem_ptr_float = (float *)memory_item->mem_ptr;
            mem_ptr_float[index+2] = mem_ptr_float[cell_course]; // top
            if(five_point){
               if (level[ntop[cell_course]] <= level[cell_course] || cell_fine == cell_left){ // 2 cells over is same or lesser refine.
                  mem_ptr_float[index+3] = mem_ptr_float[ntop[cell_course]]; // top top
               }
               else {
                  mem_ptr_float[index+3] = mem_ptr_float[nrht[ntop[cell_course]]]; // top top
               }
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            short *mem_ptr_short = (short *)memory_item->mem_ptr;
            mem_ptr_short[index+2] = mem_ptr_short[cell_course]; // top
            if(five_point){
               if (level[ntop[cell_course]] <= level[cell_course] || cell_fine == cell_left){ // 2 cells over is same or lesser refine.
                  mem_ptr_short[index+3] = mem_ptr_short[ntop[cell_course]]; // top top
               }
               else {
                  mem_ptr_short[index+3] = mem_ptr_short[nrht[ntop[cell_course]]]; // top top
               }
            }
#endif
         }
      }

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: top\n",index+2,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: toptop\n",index+3,cell_lower,cell_upper);
      break;
 
      case 4: // course cell, x-direction, right cell more refined
      // loop through state arrays to update phantom cell state values
      cell_course = cell_lower; // course neighbor to the face
      cell_fine =  cell_upper; // fine neighbor to the face
      cell_bottom = nrht[cell_course]; // bottom of the two fine neighbors of the course
      cell_top = ntop[cell_bottom]; // top of the two fince neighbots of the course

      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {

         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;

         if (memory_item->mem_elsize == 8) {
            state_sideavg = ZERO;
            double *mem_ptr_double = (double *)memory_item->mem_ptr;

            state_bot = mem_ptr_double[cell_bottom];
            state_top = mem_ptr_double[cell_top];
            state_avg = HALF * (state_bot + state_top);
            mem_ptr_double[index+2] = state_avg; //right
            //mem_ptr_double[index] = mem_ptr_double[cell_course]; //left
            if (five_point){
               if (level[nrht[cell_bottom]] > level[cell_bottom]) { // rightbot right neighbor is even more refined
                  state_botbot = mem_ptr_double[nrht[cell_bottom]];
                  state_bottop = mem_ptr_double[ntop[nrht[cell_bottom]]];
                  state_sideavg += HALF * HALF * (state_botbot + state_bottop);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_double[nrht[cell_bottom]];
               }
               if (level[nrht[cell_top]] > level[cell_top]) { // righttop right neighbor is even more refined
                  state_topbot = mem_ptr_double[nrht[cell_top]];
                  state_toptop = mem_ptr_double[ntop[nrht[cell_top]]];
                  state_sideavg += HALF * HALF * (state_topbot + state_toptop);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_double[ntop[nrht[cell_top]]];
               }
               mem_ptr_double[index+3] = state_sideavg; // right right
               //mem_ptr_double[index+1] = mem_ptr_double[nlft[cell_course]]; // left left
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
         } else if (memory_item->mem_elsize == 4) {
            state_sideavg = ZERO;
            float *mem_ptr_float = (float *)memory_item->mem_ptr;

            state_bot = mem_ptr_float[cell_bottom];
            state_top = mem_ptr_float[cell_top];
            state_avg = HALF * (state_bot + state_top);
            mem_ptr_float[index+2] = state_avg; //right
            //mem_ptr_float[index] = mem_ptr_float[cell_course]; //left
            if (five_point){
               if (level[nrht[cell_bottom]] > level[cell_bottom]) { // rightbot right neighbor is even more refined
                  state_botbot = mem_ptr_float[nrht[cell_bottom]];
                  state_bottop = mem_ptr_float[ntop[nrht[cell_bottom]]];
                  state_sideavg += HALF * HALF * (state_botbot + state_bottop);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_float[nrht[cell_bottom]];
               }
               if (level[nrht[cell_top]] > level[cell_top]) { // righttop right neighbor is even more refined
                  state_topbot = mem_ptr_float[nrht[cell_top]];
                  state_toptop = mem_ptr_float[ntop[nrht[cell_top]]];
                  state_sideavg += HALF * HALF * (state_topbot + state_toptop);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_float[ntop[nrht[cell_top]]];
               }
               mem_ptr_float[index+3] = state_sideavg; // right right
               //mem_ptr_float[index+1] = mem_ptr_float[nlft[cell_course]]; // left left
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            state_sideavg = ZERO;
            short *mem_ptr_short = (short *)memory_item->mem_ptr;

            state_bot = mem_ptr_short[cell_bottom];
            state_top = mem_ptr_short[cell_top];
            state_avg = HALF * (state_bot + state_top);
            mem_ptr_short[index+2] = state_avg; //right
            //mem_ptr_short[index] = mem_ptr_short[cell_course]; //left
            if (five_point){
               if (level[nrht[cell_bottom]] > level[cell_bottom]) { // rightbot right neighbor is even more refined
                  state_botbot = mem_ptr_short[nrht[cell_bottom]];
                  state_bottop = mem_ptr_short[ntop[nrht[cell_bottom]]];
                  state_sideavg += HALF * HALF * (state_botbot + state_bottop);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_short[nrht[cell_bottom]];
               }
               if (level[nrht[cell_top]] > level[cell_top]) { // righttop right neighbor is even more refined
                  state_topbot = mem_ptr_short[nrht[cell_top]];
                  state_toptop = mem_ptr_short[ntop[nrht[cell_top]]];
                  state_sideavg += HALF * HALF * (state_topbot + state_toptop);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_short[ntop[nrht[cell_top]]];
               }
               mem_ptr_short[index+3] = state_sideavg; // right right
               //mem_ptr_short[index+1] = mem_ptr_short[nlft[cell_course]]; // left left
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
#endif
         }
      }

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: rht\n",index+2,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: rhtrht\n",index+3,cell_lower,cell_upper);
      break;
 
      case 5: // course cell, x-direction, left cell more refine
      // loop through state arrays to update phantom cell state values
      cell_course = cell_upper; // course neighbor to the face
      cell_fine =  cell_lower; // fine neighbor to the face
      cell_bottom = nlft[cell_course]; // bottom of the two fine neighbors of the course
      cell_top = ntop[cell_bottom]; // top of the two fince neighbots of the course
 
      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {
 
         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;
 
         if (memory_item->mem_elsize == 8) {
            state_sideavg = ZERO;
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
 
            state_bot = mem_ptr_double[cell_bottom];
            state_top = mem_ptr_double[cell_top];
            state_avg = HALF * (state_bot + state_top);
            mem_ptr_double[index] = state_avg; //left
            //m_ptr_double[index+2] = mem_ptr_double[cell_course]; //right
            if(five_point){
               if (level[nlft[cell_bottom]] > level[cell_bottom]) { // leftbot left neighbor is even more refined
                  state_botbot = mem_ptr_double[nlft[cell_bottom]];
                  state_bottop = mem_ptr_double[ntop[nlft[cell_bottom]]];
                  state_sideavg += HALF * HALF * (state_botbot + state_bottop);
               }
               else { // same refinement as leftbot neighbor
                  state_sideavg += HALF * mem_ptr_double[nlft[cell_bottom]];
               }
               if (level[nlft[cell_top]] > level[cell_top]) { // lefttop left neighbor is even more refined
                  state_topbot = mem_ptr_double[nlft[cell_top]];
                  state_toptop = mem_ptr_double[ntop[nlft[cell_top]]];
                  state_sideavg += HALF * HALF * (state_topbot + state_toptop);
               }
               else { // same refinement as lefttop neighbor
                  state_sideavg += HALF * mem_ptr_double[ntop[nlft[cell_top]]];
               }
               mem_ptr_double[index+1] = state_sideavg; // left left
               //mem_ptr_double[index+3] = mem_ptr_double[nrht[cell_course]]; // right right
               // we are bot of 2 rights, so the right neighbor of the coarse will give us bottom right right neighbor
            }
         } else if (memory_item->mem_elsize == 4) {
            state_sideavg = ZERO;
            float *mem_ptr_float = (float *)memory_item->mem_ptr;
 
            state_bot = mem_ptr_float[cell_bottom];
            state_top = mem_ptr_float[cell_top];
            state_avg = HALF * (state_bot + state_top);
            mem_ptr_float[index] = state_avg; //left
            //m_ptr_float[index+2] = mem_ptr_float[cell_course]; //right
            if(five_point){
               if (level[nlft[cell_bottom]] > level[cell_bottom]) { // leftbot left neighbor is even more refined
                  state_botbot = mem_ptr_float[nlft[cell_bottom]];
                  state_bottop = mem_ptr_float[ntop[nlft[cell_bottom]]];
                  state_sideavg += HALF * HALF * (state_botbot + state_bottop);
               }
               else { // same refinement as leftbot neighbor
                  state_sideavg += HALF * mem_ptr_float[nlft[cell_bottom]];
               }
               if (level[nlft[cell_top]] > level[cell_top]) { // lefttop left neighbor is even more refined
                  state_topbot = mem_ptr_float[nlft[cell_top]];
                  state_toptop = mem_ptr_float[ntop[nlft[cell_top]]];
                  state_sideavg += HALF * HALF * (state_topbot + state_toptop);
               }
               else { // same refinement as lefttop neighbor
                  state_sideavg += HALF * mem_ptr_float[ntop[nlft[cell_top]]];
               }
               mem_ptr_float[index+1] = state_sideavg; // left left
               //mem_ptr_float[index+3] = mem_ptr_float[nrht[cell_course]]; // right right
               // we are bot of 2 rights, so the right neighbor of the coarse will give us bottom right right neighbor
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            state_sideavg = ZERO;
            short *mem_ptr_short = (short *)memory_item->mem_ptr;
 
            state_bot = mem_ptr_short[cell_bottom];
            state_top = mem_ptr_short[cell_top];
            state_avg = HALF * (state_bot + state_top);
            mem_ptr_short[index] = state_avg; //left
            //m_ptr_short[index+2] = mem_ptr_short[cell_course]; //right
            if(five_point){
               if (level[nlft[cell_bottom]] > level[cell_bottom]) { // leftbot left neighbor is even more refined
                  state_botbot = mem_ptr_short[nlft[cell_bottom]];
                  state_bottop = mem_ptr_short[ntop[nlft[cell_bottom]]];
                  state_sideavg += HALF * HALF * (state_botbot + state_bottop);
               }
               else { // same refinement as leftbot neighbor
                  state_sideavg += HALF * mem_ptr_short[nlft[cell_bottom]];
               }
               if (level[nlft[cell_top]] > level[cell_top]) { // lefttop left neighbor is even more refined
                  state_topbot = mem_ptr_short[nlft[cell_top]];
                  state_toptop = mem_ptr_short[ntop[nlft[cell_top]]];
                  state_sideavg += HALF * HALF * (state_topbot + state_toptop);
               }
               else { // same refinement as lefttop neighbor
                  state_sideavg += HALF * mem_ptr_short[ntop[nlft[cell_top]]];
               }
               mem_ptr_short[index+1] = state_sideavg; // left left
               //mem_ptr_short[index+3] = mem_ptr_short[nrht[cell_course]]; // right right
               // we are bot of 2 rights, so the right neighbor of the coarse will give us bottom right right neighbor
            }
#endif
         } 
      }

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lft\n",index,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lftlft\n",index+1,cell_lower,cell_upper);
      break;
 
      case 6: // course cell, y-direction, top cell more refined 
      // loop through state arrays to update phantom cell state values
      cell_course = cell_lower; // course neighbor to the face
      cell_fine =  cell_upper; // fine neighbor to the face
      cell_left = ntop[cell_course]; // bottom of the two fine neighbors of the course
      cell_right = nrht[cell_fine]; // top of the two fince neighbots of the course
 
      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {
 
         if ((memory_item->mem_flags & REZONE_DATA) == 0) continue;
 
         if (memory_item->mem_elsize == 8) {
            state_sideavg = ZERO;
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
 
            state_lft = mem_ptr_double[cell_left];
            state_rht = mem_ptr_double[cell_right];
            state_avg = HALF * (state_lft + state_rht);
            mem_ptr_double[index+2] = state_avg; // top
            //mem_ptr_double[index] = mem_ptr_double[cell_course]; // bottom
            if(five_point){
               if (level[ntop[cell_left]] > level[cell_left]) { // topleft top neighbor is even more refined
                  state_lftlft = mem_ptr_double[ntop[cell_left]];
                  state_lftrht = mem_ptr_double[nrht[ntop[cell_left]]];
                  state_sideavg += HALF * HALF * (state_lftlft + state_lftrht);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_double[ntop[cell_left]];
               }
               if (level[ntop[cell_right]] > level[cell_right]) { // righttop right neighbor is even more refined
                  state_rhtlft = mem_ptr_double[ntop[cell_right]];
                  state_rhtrht = mem_ptr_double[nrht[ntop[cell_right]]];
                  state_sideavg += HALF * HALF * (state_rhtlft + state_rhtrht);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_double[ntop[cell_right]];
               }
               mem_ptr_double[index+3] = state_sideavg; // top top
               //mem_ptr_double[index+1] = mem_ptr_double[nbot[cell_course]]; // bot bot
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
         } else if (memory_item->mem_elsize == 4) {
            state_sideavg = ZERO;
             float *mem_ptr_float = (float *)memory_item->mem_ptr;
 
            state_lft = mem_ptr_float[cell_left];
            state_rht = mem_ptr_float[cell_right];
            state_avg = HALF * (state_lft + state_rht);
            mem_ptr_float[index+2] = state_avg; // top
            //mem_ptr_float[index] = mem_ptr_float[cell_course]; // bottom
            if(five_point){
               if (level[ntop[cell_left]] > level[cell_left]) { // topleft top neighbor is even more refined
                  state_lftlft = mem_ptr_float[ntop[cell_left]];
                  state_lftrht = mem_ptr_float[nrht[ntop[cell_left]]];
                  state_sideavg += HALF * HALF * (state_lftlft + state_lftrht);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_float[ntop[cell_left]];
               }
               if (level[ntop[cell_right]] > level[cell_right]) { // righttop right neighbor is even more refined
                  state_rhtlft = mem_ptr_float[ntop[cell_right]];
                  state_rhtrht = mem_ptr_float[nrht[ntop[cell_right]]];
                  state_sideavg += HALF * HALF * (state_rhtlft + state_rhtrht);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_float[ntop[cell_right]];
               }
               mem_ptr_float[index+3] = state_sideavg; // top top
               //mem_ptr_float[index+1] = mem_ptr_float[nbot[cell_course]]; // bot bot
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            state_sideavg = ZERO;
             short *mem_ptr_short = (short *)memory_item->mem_ptr;
 
            state_lft = mem_ptr_short[cell_left];
            state_rht = mem_ptr_short[cell_right];
            state_avg = HALF * (state_lft + state_rht);
            mem_ptr_short[index+2] = state_avg; // top
            //mem_ptr_short[index] = mem_ptr_short[cell_course]; // bottom
            if(five_point){
               if (level[ntop[cell_left]] > level[cell_left]) { // topleft top neighbor is even more refined
                  state_lftlft = mem_ptr_short[ntop[cell_left]];
                  state_lftrht = mem_ptr_short[nrht[ntop[cell_left]]];
                  state_sideavg += HALF * HALF * (state_lftlft + state_lftrht);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_short[ntop[cell_left]];
               }
               if (level[ntop[cell_right]] > level[cell_right]) { // righttop right neighbor is even more refined
                  state_rhtlft = mem_ptr_short[ntop[cell_right]];
                  state_rhtrht = mem_ptr_short[nrht[ntop[cell_right]]];
                  state_sideavg += HALF * HALF * (state_rhtlft + state_rhtrht);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_short[ntop[cell_right]];
               }
               mem_ptr_short[index+3] = state_sideavg; // top top
               //mem_ptr_short[index+1] = mem_ptr_short[nbot[cell_course]]; // bot bot
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
#endif
         }
      }

      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: top\n",index+2,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: toptop\n",index+3,cell_lower,cell_upper);
      break;
 
      case 7: // course cell, y-direction, bottom cell more refined
      // loop through state arrays to update phantom cell state values
      cell_course = cell_upper; // course neighbor to the face
      cell_fine =  cell_lower; // fine neighbor to the face
      cell_left = nbot[cell_course]; // bottom of the two fine neighbors of the course
      cell_right = nrht[cell_fine]; // top of the two fince neighbots of the course
 
      for (memory_item = state_memory_old.memory_entry_by_name_begin();
         memory_item != state_memory_old.memory_entry_by_name_end();
         memory_item = state_memory_old.memory_entry_by_name_next() ) {
 
         if ( (memory_item->mem_flags & REZONE_DATA) == 0) continue;
 
         if (memory_item->mem_elsize == 8) {
            state_sideavg = ZERO;
            double *mem_ptr_double = (double *)memory_item->mem_ptr;
 
            state_lft = mem_ptr_double[cell_left];
            state_rht = mem_ptr_double[cell_right];
            state_avg = HALF * (state_lft + state_rht);
            mem_ptr_double[index] = state_avg; // bottom
            //mem_ptr_double[index+2] = mem_ptr_double[cell_course]; // top
            if (five_point){
               if (level[nbot[cell_left]] > level[cell_left]) { // botleft bot neighbor is even more refined
                  state_lftlft = mem_ptr_double[nbot[cell_left]];
                  state_lftrht = mem_ptr_double[nrht[nbot[cell_left]]];
                  state_sideavg += HALF * HALF * (state_lftlft + state_lftrht);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_double[nbot[cell_left]];
               }
               if (level[nbot[cell_right]] > level[cell_right]) { // righttop right neighbor is even more refined
                  state_rhtlft = mem_ptr_double[nbot[cell_right]];
                  state_rhtrht = mem_ptr_double[nrht[nbot[cell_right]]];
                  state_sideavg += HALF * HALF * (state_rhtlft + state_rhtrht);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_double[nbot[cell_right]];
               }
               mem_ptr_double[index+1] = state_sideavg; // bot bot
               //mem_ptr_double[index+3] = mem_ptr_double[ntop[cell_course]]; // top top
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
         } else if (memory_item->mem_elsize == 4) {
            state_sideavg = ZERO;
            float *mem_ptr_float = (float *)memory_item->mem_ptr;
 
            state_lft = mem_ptr_float[cell_left];
            state_rht = mem_ptr_float[cell_right];
            state_avg = HALF * (state_lft + state_rht);
            mem_ptr_float[index] = state_avg; // bottom
            //mem_ptr_float[index+2] = mem_ptr_float[cell_course]; // top
            if (five_point){
               if (level[nbot[cell_left]] > level[cell_left]) { // botleft bot neighbor is even more refined
                  state_lftlft = mem_ptr_float[nbot[cell_left]];
                  state_lftrht = mem_ptr_float[nrht[nbot[cell_left]]];
                  state_sideavg += HALF * HALF * (state_lftlft + state_lftrht);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_float[nbot[cell_left]];
               }
               if (level[nbot[cell_right]] > level[cell_right]) { // righttop right neighbor is even more refined
                  state_rhtlft = mem_ptr_float[nbot[cell_right]];
                  state_rhtrht = mem_ptr_float[nrht[nbot[cell_right]]];
                  state_sideavg += HALF * HALF * (state_rhtlft + state_rhtrht);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_float[nbot[cell_right]];
               }
               mem_ptr_float[index+1] = state_sideavg; // bot bot
               //mem_ptr_float[index+3] = mem_ptr_float[ntop[cell_course]]; // top top
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
#ifdef HALF_PRECISION
         } else if (memory_item->mem_elsize == 2) {
            state_sideavg = ZERO;
            short *mem_ptr_short = (short *)memory_item->mem_ptr;
 
            state_lft = mem_ptr_short[cell_left];
            state_rht = mem_ptr_short[cell_right];
            state_avg = HALF * (state_lft + state_rht);
            mem_ptr_short[index] = state_avg; // bottom
            //mem_ptr_short[index+2] = mem_ptr_short[cell_course]; // top
            if (five_point){
               if (level[nbot[cell_left]] > level[cell_left]) { // botleft bot neighbor is even more refined
                  state_lftlft = mem_ptr_short[nbot[cell_left]];
                  state_lftrht = mem_ptr_short[nrht[nbot[cell_left]]];
                  state_sideavg += HALF * HALF * (state_lftlft + state_lftrht);
               }
               else { // same refinement as rightbot neighbor
                  state_sideavg += HALF * mem_ptr_short[nbot[cell_left]];
               }
               if (level[nbot[cell_right]] > level[cell_right]) { // righttop right neighbor is even more refined
                  state_rhtlft = mem_ptr_short[nbot[cell_right]];
                  state_rhtrht = mem_ptr_short[nrht[nbot[cell_right]]];
                  state_sideavg += HALF * HALF * (state_rhtlft + state_rhtrht);
               }
               else { // same refinement as righttop neighbor
                  state_sideavg += HALF * mem_ptr_short[nbot[cell_right]];
               }
               mem_ptr_short[index+1] = state_sideavg; // bot bot
               //mem_ptr_short[index+3] = mem_ptr_short[ntop[cell_course]]; // top top
               // we are bot of 2 lefts, so the left neighbor of the coarse will give us bottom left left neighbor
            }
#endif
         }
      }

      //printf("DEBUG MESH: ID %d) LOWER: %d, UPPER: %d, POS: bot\n",index,cell_lower,cell_upper);
      //printf("DEBUG MESH: ID %d) LOWER: %d, UPPER: %d, POS: botbot \n",index+1,cell_lower,cell_upper);
      break;

   }

}

void Mesh::calc_face_list_wbidirmap_phantom(MallocPlus &state_memory, double deltaT)
{
   if (!do_rezone && !firstFlag) {
      firstFlag = false;
      calc_face_list_fill_phantom(state_memory, deltaT);
      return;
   }
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

    // These variables are for the following for-loop
    int pcellCnt; // counter for new phantom cells

    int xfaceCnt = 0;
    int yfaceCnt = 0;
    int xpCnt = 0;
    int ypCnt = 0;
    int cellCount = 0;
    int xfix = 0;
    int yfix = 0;
    for (int nz = 0; nz < (int)ncells; nz++) {
        int nl = nlft[nz];
        int nr = nrht[nz];
        int nb = nbot[nz];
        int nt = ntop[nz];
        if (nz == nl || nz == nr || nz == nb || nz == nt) {
            if (nz != nr) xfaceCnt ++;
            if (nz != nt) yfaceCnt ++;
            continue;
        }

        int lev = level[nz];
        //int ll = level[nl];
        int lr = level[nr];
        //int lb = level[nb];
        int lt = level[nt];

        if (lev < lr) {
            xfaceCnt += 2;
            xpCnt += 4;
            cellCount += 6;
        }
        else if (lev == lr) {
            xfaceCnt ++;
        }
        else {
            cellCount += 3;
            xfaceCnt ++;
            xpCnt += 2;
            xfix ++;
        }
        if (lev < lt) {
            yfaceCnt += 2;
            ypCnt += 4;
            cellCount += 6;
        }
        else if (lev == lt) {
            yfaceCnt ++;
        }
        else {
            cellCount += 3;
            yfaceCnt ++;
            ypCnt += 2;
            yfix ++;
        }
    }
    nxface = xfaceCnt;
    nyface = yfaceCnt;
    pcellCnt = ncells + cellCount;
    pxface = nxface + xpCnt;
    pyface = nyface + ypCnt;
    nxfixup = xfix / 2;
    nyfixup = yfix / 2;

    // realloc memory based on new counts
    i        = (int *)mesh_memory.memory_realloc(pcellCnt, i);
    j        = (int *)mesh_memory.memory_realloc(pcellCnt, j);
    level    = (uchar_t *)mesh_memory.memory_realloc(pcellCnt, level);
    nlft     = (int *)mesh_memory.memory_realloc(pcellCnt, nlft);
    nrht     = (int *)mesh_memory.memory_realloc(pcellCnt, nrht);
    nbot     = (int *)mesh_memory.memory_realloc(pcellCnt, nbot);
    ntop     = (int *)mesh_memory.memory_realloc(pcellCnt, ntop);

    map_xface2cell_lower = (int *)mesh_memory.memory_realloc(pxface, map_xface2cell_lower);
    map_xface2cell_upper = (int *)mesh_memory.memory_realloc(pxface, map_xface2cell_upper);
    map_yface2cell_lower = (int *)mesh_memory.memory_realloc(pyface, map_yface2cell_lower);
    map_yface2cell_upper = (int *)mesh_memory.memory_realloc(pyface, map_yface2cell_upper);

    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
    xface_i = (int *)mesh_memory.memory_realloc(pxface, xface_i);
    xface_j = (int *)mesh_memory.memory_realloc(pxface, xface_j);
    xface_level = (uchar_t *)mesh_memory.memory_realloc(pxface, xface_level);
    yface_i = (int *)mesh_memory.memory_realloc(pyface, yface_i);
    yface_j = (int *)mesh_memory.memory_realloc(pyface, yface_j);
    yface_level = (uchar_t *)mesh_memory.memory_realloc(pyface, yface_level);
    }

    map_xcell2face_left1 = (int *)mesh_memory.memory_realloc(pcellCnt, map_xcell2face_left1);
    map_xcell2face_left2 = (int *)mesh_memory.memory_realloc(ncells, map_xcell2face_left2);
    map_xcell2face_right1 = (int *)mesh_memory.memory_realloc(pcellCnt, map_xcell2face_right1);
    map_xcell2face_right2 = (int *)mesh_memory.memory_realloc(ncells, map_xcell2face_right2);
    map_ycell2face_bot1 = (int *)mesh_memory.memory_realloc(pcellCnt, map_ycell2face_bot1);
    map_ycell2face_bot2 = (int *)mesh_memory.memory_realloc(ncells, map_ycell2face_bot2);
    map_ycell2face_top1 = (int *)mesh_memory.memory_realloc(pcellCnt, map_ycell2face_top1);
    map_ycell2face_top2 = (int *)mesh_memory.memory_realloc(ncells, map_ycell2face_top2);

    xrecvIdx = (int *)mesh_memory.memory_realloc(ncells, xrecvIdx);
    xrecvCIdx = (int *)mesh_memory.memory_realloc(ncells, xrecvCIdx);
    xplusCell2Idx = (int *)mesh_memory.memory_realloc(ncells, xplusCell2Idx);
    xminusCell2Idx = (int *)mesh_memory.memory_realloc(ncells, xminusCell2Idx);
    xsendIdx1 = (int *)mesh_memory.memory_realloc(ncells, xsendIdx1);
    xsendIdx2 = (int *)mesh_memory.memory_realloc(ncells, xsendIdx2);
    yrecvIdx = (int *)mesh_memory.memory_realloc(ncells, yrecvIdx);
    yrecvCIdx = (int *)mesh_memory.memory_realloc(ncells, yrecvCIdx);
    yplusCell2Idx = (int *)mesh_memory.memory_realloc(ncells, yplusCell2Idx);
    yminusCell2Idx = (int *)mesh_memory.memory_realloc(ncells, yminusCell2Idx);
    ysendIdx1 = (int *)mesh_memory.memory_realloc(ncells, ysendIdx1);
    ysendIdx2 = (int *)mesh_memory.memory_realloc(ncells, ysendIdx2);

    memory_reset_ptrs();

    // needs to be initialized to -1
    for (uint fill = 0; fill < ncells; fill++) {
        xplusCell2Idx[fill] = -1;
        xminusCell2Idx[fill] = -1;
        yplusCell2Idx[fill] = -1;
        yminusCell2Idx[fill] = -1;
    }

    // realloc state memory
    MallocPlus state_memory_old = state_memory;
    malloc_plus_memory_entry *memory_item;

    for (memory_item = state_memory_old.memory_entry_by_name_begin();
        memory_item != state_memory_old.memory_entry_by_name_end();
        memory_item = state_memory_old.memory_entry_by_name_next() ) {

        if ( (memory_item->mem_flags & REZONE_DATA) == 0) continue;
        
        state_memory.memory_realloc(pcellCnt, memory_item->mem_ptr);

    }  

   /*
   map_xface2cell_lower.clear();
   map_xface2cell_upper.clear();
   map_xface2cell_lower.resize(2*ncells, -1);
   map_xface2cell_upper.resize(2*ncells, -1);

   xrecvIdx.clear();
   xrecvCIdx.clear();
   xplusCell2Idx.clear();
   xminusCell2Idx.clear();
   xsendIdx1.clear();
   xsendIdx2.clear();
   xrecvIdx.resize(ncells,-1);
   xrecvCIdx.resize(ncells,-1);
   xplusCell2Idx.resize(ncells,-1);
   xminusCell2Idx.resize(ncells,-1);
   xsendIdx1.resize(ncells,-1);
   xsendIdx2.resize(ncells,-1);

   //
   map_xcell2face_left1.clear();
   map_xcell2face_left2.clear();
   map_xcell2face_right1.clear();
   map_xcell2face_right2.clear();
   map_xcell2face_left1.resize(ncells, -1);
   map_xcell2face_left2.resize(ncells, -1);
   map_xcell2face_right1.resize(ncells, -1);
   map_xcell2face_right2.resize(ncells, -1);

   //
   xface_i.clear();
   xface_j.clear();
   xface_level.clear();
   xface_i.resize(2*ncells, -1);
   xface_j.resize(2*ncells, -1);
   xface_level.resize(2*ncells, -1);

   map_yface2cell_lower.clear();
   map_yface2cell_upper.clear();
   map_yface2cell_lower.resize(2*ncells, -1);
   map_yface2cell_upper.resize(2*ncells, -1);

   yrecvIdx.clear();
   yrecvCIdx.clear();
   yplusCell2Idx.clear();
   yminusCell2Idx.clear();
   ysendIdx1.clear();
   ysendIdx2.clear();
   yrecvIdx.resize(ncells,-1);
   yrecvCIdx.resize(ncells,-1);
   yplusCell2Idx.resize(ncells,-1);
   yminusCell2Idx.resize(ncells,-1);
   ysendIdx1.resize(ncells,-1);
   ysendIdx2.resize(ncells,-1);

   //
   map_ycell2face_bot1.clear();
   map_ycell2face_bot2.clear();
   map_ycell2face_top1.clear();
   map_ycell2face_top2.clear();
   map_ycell2face_bot1.resize(ncells, -1);
   map_ycell2face_bot2.resize(ncells, -1);
   map_ycell2face_top1.resize(ncells, -1);
   map_ycell2face_top2.resize(ncells, -1);

   //
   yface_i.clear();
   yface_j.clear();
   yface_level.clear();
   yface_i.resize(2*ncells, -1);
   yface_j.resize(2*ncells, -1);
   yface_level.resize(2*ncells, -1);
   */

   cpu_timers[MESH_TIMER_BIDIRPART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   int iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nr = nrht[nz];
      if (nr == nz) continue;

      int ifactor = 1;
      if (level[nr] < level[nz]) ifactor = 2;

      // Have right face
      map_xface2cell_lower[iface] = nz;
      map_xface2cell_upper[iface] = nr;
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
      xface_level[iface] = MAX(level[nz],level[nr]);
      xface_i[iface] = i[nr]*ifactor;
      if (level[nr] < level[nz] && is_upper(j[nz]) ) {
         xface_j[iface] = j[nr]*ifactor+1;
      } else {
         xface_j[iface] = j[nr]*ifactor;
      }
    }
      map_xcell2face_right1[nz] = iface;

      //the right is a real cell, but I am left boundary
      if (nz == nlft[nz])
          map_xcell2face_left1[nz] = iface;

      iface++;

      if (level[nr] > level[nz] && is_lower(j[nr]) ){
         int ntr = ntop[nr];
         if (ntr != nr) {
            map_xface2cell_lower[iface] = nz;
            map_xface2cell_upper[iface] = ntr;
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
            xface_level[iface] = MAX(level[nz],level[ntr]);
            xface_i[iface] = i[ntr]*ifactor;
            xface_j[iface] = j[ntr]*ifactor;
    }
            map_xcell2face_right2[nz] = iface;

            iface++;
         }
      }
   }

   cpu_timers[MESH_TIMER_BIDIRPART2] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int nz=0; nz<(int)ncells; nz++){
      int nl = nlft[nz];
      if (nl == nz) continue;

      if (level[nl] < level[nz] && is_upper(j[nz]) && nbot[nz] != nz){
         map_xcell2face_left1[nz] = map_xcell2face_right2[nl];


      } else {
         map_xcell2face_left1[nz] = map_xcell2face_right1[nl];

         if (level[nl] > level[nz]){
            map_xcell2face_left2[nz] = map_xcell2face_right1[ntop[nl]];
         }
      }

      //the left is a real cell, but I am right boundary
      if (nz == nrht[nz])
          map_xcell2face_right1[nz] = map_xcell2face_left1[nz];

   }
   cpu_timers[MESH_TIMER_BIDIRPART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   //
   cpu_timer_start(&tstart_cpu_part);

   iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nt = ntop[nz];
      if (nt == nz) continue;

      int ifactor = 1;
      if (level[nt] < level[nz]) ifactor = 2;

      // Have top face
      map_yface2cell_lower[iface] = nz;
      map_yface2cell_upper[iface] = nt;
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
      yface_level[iface] = MAX(level[nz],level[nt]);
      yface_j[iface] = j[nt]*ifactor;
      if (level[nt] < level[nz] && is_upper(i[nz]) ) {
         yface_i[iface] = i[nt]*ifactor+1;
      } else{
         yface_i[iface] = i[nt]*ifactor;
      }
    }
      map_ycell2face_top1[nz] = iface;

      //the top is a real cell, but I am bot boundary
      if (nz == nbot[nz])
          map_ycell2face_bot1[nz] = iface;

      iface++;

      if (level[nt] > level[nz]  && is_lower(i[nt]) ){
         int nrt = nrht[nt];
         if (nrt != nt) {
            map_yface2cell_lower[iface] = nz;
            map_yface2cell_upper[iface] = nrt;
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
            yface_level[iface] = MAX(level[nz],level[nrt]);
            yface_j[iface] = j[nrt]*ifactor;
            yface_i[iface] = i[nrt]*ifactor;
    }
            map_ycell2face_top2[nz] = iface;

            iface++;
         }
      }
   }

   cpu_timers[MESH_TIMER_BIDIRPART5] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int nz=0; nz<(int)ncells; nz++){
      int nb = nbot[nz];
      if (nb == nz) continue;
      //if (nz == 3) printf("%d %d\n", map_ycell2face_top1[nb], map_ycell2face_top2[nb]);

      if (level[nb] < level[nz] && is_upper(i[nz]) && nlft[nz] != nz){
         map_ycell2face_bot1[nz] = map_ycell2face_top2[nb];

      } else {
         map_ycell2face_bot1[nz] = map_ycell2face_top1[nb];

         if (level[nb] > level[nz]){
            map_ycell2face_bot2[nz] = map_ycell2face_top1[nrht[nb]];
         }
      }
      //the bot is a real cell, but I am top boundary
      if (nz == ntop[nz])
          map_ycell2face_top1[nz] = map_ycell2face_bot1[nz];
   }

   //for (int nz = 0; nz < ncells; nz++) {
   //    printf("%d) %d %d %d %d\n", nz,map_xcell2face_left1[nz],map_xcell2face_right1[nz],map_ycell2face_bot1[nz],map_ycell2face_top1[nz]);
   //}

   cpu_timers[MESH_TIMER_BIDIRPART6] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);


   cpu_timers[MESH_TIMER_BIDIRPART7] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

/*
    for (int iface = 0; iface < nxface; iface++) {
        int lncell = map_xface2cell_lower[iface], // cell neighbor to the left
            rncell = map_xface2cell_upper[iface], // cell neighbor to the right
            cncell, // course neighbor to the face
            fncell, // fine neighbor to the face
            bncell, // bottom of the two fine neighbors of the course
            tncell, // top of the two fince neighbots of the course
            level_left = level[lncell],
            level_right = level[rncell];
        if (nlft[lncell] == lncell || nbot[lncell] == lncell || ntop[lncell] == lncell) continue;
        if (nrht[rncell] == rncell || nbot[rncell] == rncell || ntop[rncell] == rncell) continue;

        // for future indexing, the new phantom cell appear in the array as follows
        // lp llp rp rrp (for even iface), [l,r]p [ll,rr]p (for odd iface)
        // important to note we only add 2 phantom cells on odd iface values
        // because 2 that would be added by its even "partner" iface value will be the same 
        if (level_left != level_right) {
            if (level_left < level_right) {
                cncell = lncell;
                fncell = rncell;
                bncell = nrht[cncell];
                tncell = ntop[bncell];
            }
            else {
                cncell = rncell;
                fncell = lncell;
                bncell = nlft[cncell];
                tncell = ntop[bncell];
            }
            if (fncell == bncell) { // bottom of the two horizontal neighbors, add 4 phantoms
                pcellCnt += 4;
                pfaceCnt += 3;
            }
            else {
                pcellCnt += 2;
                pfaceCnt ++;
            }
        }
    }



    map_xface2cell_lower.resize(pfaceCnt, -1);
    map_xface2cell_upper.resize(pfaceCnt, -1);
    xface_i.resize(pfaceCnt, -1);
    xface_j.resize(pfaceCnt, -1);
    xface_level.resize(pfaceCnt, -1);

    pfaceCnt = nyface;
    for (int iface = 0; iface < nyface; iface++) {
     int    bncell = map_yface2cell_lower[iface], // cell neighbor below
            tncell = map_yface2cell_upper[iface], // cell neighbor above
            cncell, // course neighbor to the face
            fncell, // fine neighbor to the face
            lncell, // bottom of the two fine neighbors of the course
            rncell, // top of the two fince neighbots of the course
            level_bot = level[bncell],
            level_top = level[tncell];
        if (nlft[bncell] == bncell || nbot[bncell] == bncell || nrht[bncell] == bncell) continue;
        if (nlft[tncell] == tncell || nrht[tncell] == tncell || ntop[tncell] == tncell) continue;

        // for future indexing, the new phantom cell appear in the array as follows
        // bp bbp tp ttp (for even iface), [b,t]p [bb,tt]p (for odd iface)
        // important to note we only add 2 phantom cells on odd iface values
        // because 2 that would be added by its even "partner" iface value will be the same 
        if (level_bot != level_top) {
            if (level_bot < level_top) {
                cncell = bncell;
                fncell = tncell;
                lncell = ntop[cncell];
                rncell = nrht[lncell];
            }
            else {
                cncell = tncell;
                fncell = bncell;
                lncell = nbot[cncell];
                rncell = nrht[lncell];
            }
            if (fncell == lncell) { // left of the two vertical neighbors, add 4 phantoms
                pcellCnt += 4;
                pfaceCnt += 3;
            }
            else {
                pcellCnt += 2;
                pfaceCnt ++;
            }
        }
    }

    //Reallocate extra memory based on the precount
    map_xcell2face_left1.resize(pcellCnt, -1);
    map_xcell2face_left2.resize(pcellCnt, -1);
    map_xcell2face_right1.resize(pcellCnt, -1);
    map_xcell2face_right2.resize(pcellCnt, -1);

    map_ycell2face_bot1.resize(pcellCnt, -1);
    map_ycell2face_top1.resize(pcellCnt, -1);
    map_ycell2face_bot2.resize(pcellCnt, -1);
    map_ycell2face_top2.resize(pcellCnt, -1);
    map_yface2cell_lower.resize(pfaceCnt, -1);
    map_yface2cell_upper.resize(pfaceCnt, -1);
    yface_i.resize(pfaceCnt, -1);
    yface_j.resize(pfaceCnt, -1);
    yface_level.resize(pfaceCnt, -1);
    */


    /*for (int braat = 0; braat < ncells; braat++) {
        printf("%d) %d %d %d %d\n", braat, map_xcell2face_left1[braat], map_xcell2face_right1[braat], map_ycell2face_bot1[braat], map_ycell2face_top1[braat]);
    }*/

   cpu_timers[MESH_TIMER_BIDIRPART8] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

    int idxVar = 0,
        pcellIdx = (int) ncells, // starting index for new cell phantoms
        pfaceIdx = nxface, // starting index for new face phantoms
        ifixupIdx = 0; // starting index for receive fluxfixup

    for (int iface = 0; iface < nxface; iface++) {
        int lncell = map_xface2cell_lower[iface], // cell neighbor to the left
            rncell = map_xface2cell_upper[iface], // cell neighbor to the right
            cncell, // course neighbor to the face
            fncell, // fine neighbor to the face
            bncell, // bottom of the two fine neighbors of the course
            tncell, // top of the two fince neighbots of the course
            level_left = level[lncell],
            level_right = level[rncell];
        if (nlft[lncell] == lncell || nbot[lncell] == lncell || ntop[lncell] == lncell) continue;
        if (nrht[rncell] == rncell || nbot[rncell] == rncell || ntop[rncell] == rncell) continue;

        // for future indexing, the new phantom cell appear in the array as follows
        // lp llp rp rrp (for even iface), [l,r]p [ll,rr]p (for odd iface)
        // important to note we only add 2 phantom cells on odd iface values
        // because 2 that would be added by its even "partner" iface value will be the same 
        if (level_left != level_right) {
            if (level_left < level_right) {
                cncell = lncell;
                fncell = rncell;
                bncell = nrht[cncell];
                tncell = ntop[bncell];
            }
            else {
                cncell = rncell;
                fncell = lncell;
                bncell = nlft[cncell];
                tncell = ntop[bncell];
            }
            if (fncell == bncell) { // bottom of the two horizontal neighbors, add 4 phantoms
                idxVar = 0;
                //pcellCnt += 4;
                //pfaceCnt += 3;

                if (level_left < level_right) { // right is more refined
                    // new face's adjacent cells
                    map_xface2cell_upper[pfaceIdx] = pcellIdx+2;
                    map_xface2cell_lower[pfaceIdx] = lncell;        
                    map_xface2cell_upper[pfaceIdx+1] = pcellIdx;
                    map_xface2cell_lower[pfaceIdx+1] = pcellIdx+1;
                    map_xface2cell_upper[pfaceIdx+2] = pcellIdx+3;
                    map_xface2cell_lower[pfaceIdx+2] = pcellIdx+2;
                    // adjacent cell's face
                    map_xcell2face_right1[lncell] = pfaceIdx;
                    map_xcell2face_left1[pcellIdx+2] = pfaceIdx;
                    map_xcell2face_right1[pcellIdx+1] = pfaceIdx+1;
                    map_xcell2face_left1[pcellIdx] = pfaceIdx+1;
                    map_xcell2face_right1[pcellIdx+2] = pfaceIdx+2;
                    map_xcell2face_left1[pcellIdx+3] = pfaceIdx+2;
                    map_xcell2face_right1[pcellIdx] = iface;
                    // "undo" the second face for the cell, as it is no longer applicable
                    map_xcell2face_right2[lncell] = -1;
                    // old face's new phantom adjacent cell
                    map_xface2cell_lower[iface] = pcellIdx;

                    // face positions
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    xface_level[pfaceIdx] = level[lncell];
                    xface_level[pfaceIdx+1] = level[rncell];
                    xface_level[pfaceIdx+2] = level[lncell];
                    xface_i[pfaceIdx] = i[lncell] + 1;
                    xface_i[pfaceIdx+1] = i[rncell] - 1;
                    xface_i[pfaceIdx+2] = i[lncell] + 2;
                    xface_j[pfaceIdx] = j[lncell];
                    xface_j[pfaceIdx+1] = j[rncell];
                    xface_j[pfaceIdx+2] = j[lncell];
    }

                    interpolate(0, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                    interpolate(4, pcellIdx, lncell, rncell, deltaT,  state_memory_old);

                    xrecvIdx[ifixupIdx] = pfaceIdx;
                    int nl = map_xface2cell_lower[pfaceIdx];
                    xrecvCIdx[ifixupIdx] = nl;
                    xplusCell2Idx[nl] = ifixupIdx;
                    if (iface < map_xcell2face_left1[tncell]) {
                       xsendIdx1[ifixupIdx] = iface;
                       xsendIdx2[ifixupIdx] = map_xcell2face_left1[tncell];
                    } else {
                       xsendIdx1[ifixupIdx] = map_xcell2face_left1[tncell];
                       xsendIdx2[ifixupIdx] = iface;
                    }
                    ifixupIdx++;



                }
                else { // left is more refined
                    // new face's adjacent cells
                    map_xface2cell_lower[pfaceIdx] = pcellIdx;        
                    map_xface2cell_upper[pfaceIdx] = rncell;
                    map_xface2cell_lower[pfaceIdx+1] = pcellIdx+1;
                    map_xface2cell_upper[pfaceIdx+1] = pcellIdx;
                    map_xface2cell_lower[pfaceIdx+2] = pcellIdx+2;
                    map_xface2cell_upper[pfaceIdx+2] = pcellIdx+3;
                    // adjacent cell's face
                    map_xcell2face_left1[rncell] = pfaceIdx;
                    map_xcell2face_right1[pcellIdx] = pfaceIdx;
                    map_xcell2face_left1[pcellIdx] = pfaceIdx+1;
                    map_xcell2face_right1[pcellIdx+1] = pfaceIdx+1;
                    map_xcell2face_left1[pcellIdx+3] = pfaceIdx+2;
                    map_xcell2face_right1[pcellIdx+2] = pfaceIdx+2;
                    map_xcell2face_left1[pcellIdx+2] = iface;
                    // "undo" the second face for the cell, as it is no longer applicable
                    map_xcell2face_left2[rncell] = -1;
                    // old face's new phantom adjacent cell
                    map_xface2cell_upper[iface] = pcellIdx + 2;

    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    xface_level[pfaceIdx] = level[rncell];
                    xface_level[pfaceIdx+1] = level[rncell];
                    xface_level[pfaceIdx+2] = level[lncell];
                    xface_i[pfaceIdx] = i[rncell];
                    xface_i[pfaceIdx+1] = i[rncell] - 1;
                    xface_i[pfaceIdx+2] = i[lncell] + 2;
                    xface_j[pfaceIdx] = j[rncell];
                    xface_j[pfaceIdx+1] = j[rncell];
                    xface_j[pfaceIdx+2] = j[lncell];
    }

                    interpolate(1, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                    interpolate(5, pcellIdx, lncell, rncell, deltaT,  state_memory_old);

                    xrecvIdx[ifixupIdx] = pfaceIdx;
                    int nr = map_xface2cell_upper[pfaceIdx];
                    xrecvCIdx[ifixupIdx] = nr;
                    xminusCell2Idx[nr] = ifixupIdx;
                    if (iface < map_xcell2face_right1[tncell]){
                       xsendIdx1[ifixupIdx] = iface;
                       xsendIdx2[ifixupIdx] = map_xcell2face_right1[tncell];
                    } else {
                       xsendIdx1[ifixupIdx] = map_xcell2face_right1[tncell];
                       xsendIdx2[ifixupIdx] = iface;
                    }
                    ifixupIdx++;


                }
                // phantom cells' new neighbors (same regardless of which side refinement occurs)
                nlft[pcellIdx] = pcellIdx + 1;
                nrht[pcellIdx] = rncell;
                nrht[pcellIdx+1] = pcellIdx;
                nlft[pcellIdx+2] = lncell;
                nrht[pcellIdx+2] = pcellIdx + 3;
                nlft[pcellIdx+3] = pcellIdx + 2;

                // update other arrays
                i[pcellIdx] = i[rncell] - 1;
                i[pcellIdx+1] = i[rncell] - 2;
                i[pcellIdx+2] = i[lncell] + 1;
                i[pcellIdx+3] = i[lncell] + 2;
                j[pcellIdx] = j[rncell];
                j[pcellIdx+1] = j[rncell];
                j[pcellIdx+2] = j[lncell];
                j[pcellIdx+3] = j[lncell];
                level[pcellIdx] = level[rncell];
                level[pcellIdx+1] = level[rncell];
                level[pcellIdx+2] = level[lncell];
                level[pcellIdx+3] = level[lncell];

            }
            else { // top of the two horizontal neighbors, add 2 phantoms
                idxVar = 1;
                //pcellCnt += 2;
                //pfaceCnt++;

                if (level[lncell] < level[rncell]) { // right is more refined
                    // old face's new phantom adjacent cell
                    map_xface2cell_lower[iface] = pcellIdx;
                    map_xface2cell_upper[pfaceIdx] = pcellIdx;
                    map_xface2cell_lower[pfaceIdx] = pcellIdx+1;

                    map_xcell2face_left1[pcellIdx] = pfaceIdx;
                    map_xcell2face_right1[pcellIdx+1] = pfaceIdx;
                    map_xcell2face_right1[pcellIdx] = iface;
                    // phantom cells' new neighbors
                    nlft[pcellIdx] = pcellIdx + 1;
                    nrht[pcellIdx] = rncell;
                    nrht[pcellIdx+1] = pcellIdx;
                    // other arrays
                    i[pcellIdx] = i[rncell] - 1;
                    i[pcellIdx+1] = i[rncell] - 2;
                    j[pcellIdx] = j[rncell];
                    j[pcellIdx+1] = j[rncell];
                    level[pcellIdx] = level[rncell];
                    level[pcellIdx+1] = level[rncell];

    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    xface_level[pfaceIdx] = level[rncell];
                    xface_i[pfaceIdx] = i[rncell] - 1;
                    xface_j[pfaceIdx] = j[rncell];
    }


                    interpolate(0, pcellIdx, lncell, rncell, deltaT,  state_memory_old);


                }
                else { // left is more refined 
                    // old face's new phantom adjacent cell
                    map_xface2cell_upper[iface] = pcellIdx; 
                    map_xface2cell_upper[pfaceIdx] = pcellIdx+1;
                    map_xface2cell_lower[pfaceIdx] = pcellIdx;

                    map_xcell2face_left1[pcellIdx+1] = pfaceIdx;
                    map_xcell2face_right1[pcellIdx] = pfaceIdx;
                    map_xcell2face_left1[pcellIdx] = iface;
                    // phantom cells' new neighbors
                    nlft[pcellIdx] = lncell;
                    nrht[pcellIdx] = pcellIdx + 1;
                    nrht[pcellIdx+1] = pcellIdx;
                    // other arrays
                    i[pcellIdx] = i[lncell] + 1;
                    i[pcellIdx+1] = i[lncell] + 2;
                    j[pcellIdx] = j[lncell];
                    j[pcellIdx+1] = j[lncell];
                    level[pcellIdx] = level[lncell];
                    level[pcellIdx+1] = level[lncell];

    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    xface_level[pfaceIdx] = level[lncell];
                    xface_i[pfaceIdx] = i[lncell] + 2;
                    xface_j[pfaceIdx] = j[lncell];
    }

                    interpolate(1, pcellIdx-2, lncell, rncell, deltaT,  state_memory_old);

                }

            }
#ifdef PATTERN_CHECK
            xcase[pfaceIdx] = 0;
#endif

            // update indexes
            pcellIdx += 4 - (idxVar % 2) * 2;
            pfaceIdx += 3 - (idxVar % 2) * 2;
        }

    } 

   nxfixup = ifixupIdx;

   cpu_timers[MESH_TIMER_BIDIRPART9] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

     // resize arrays/vectors here including the addition of the new faces and cells

    /*map_xface2cell_lower.resize(pfaceIdx);
    map_xface2cell_upper.resize(pfaceIdx);
    map_xcell2face_left1.resize(pcellIdx);
    map_xcell2face_left2.resize(pcellIdx);
    map_xcell2face_right1.resize(pcellIdx);
    map_xcell2face_right2.resize(pcellIdx);
    xface_i.resize(pfaceIdx);
    xface_j.resize(pfaceIdx);
    xface_level.resize(pfaceIdx);*/

   cpu_timers[MESH_TIMER_BIDIRPART10] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

    // now for the y faces
    
    // pcellIdx will actually stay where it was previously incremented to because there
    // is only 1 array of cells (not dependent on x/y)
    pfaceIdx = nyface;
    ifixupIdx = 0;

    for (int iface = 0; iface < nyface; iface++) {
     int    bncell = map_yface2cell_lower[iface], // cell neighbor below
            tncell = map_yface2cell_upper[iface], // cell neighbor above
            cncell, // course neighbor to the face
            fncell, // fine neighbor to the face
            lncell, // bottom of the two fine neighbors of the course
            rncell, // top of the two fince neighbots of the course
            level_bot = level[bncell],
            level_top = level[tncell];
        if (nlft[bncell] == bncell || nbot[bncell] == bncell || nrht[bncell] == bncell) continue;
        if (nrht[tncell] == tncell || nrht[tncell] == tncell || ntop[tncell] == tncell) continue;

        // for future indexing, the new phantom cell appear in the array as follows
        // bp bbp tp ttp (for even iface), [b,t]p [bb,tt]p (for odd iface)
        // important to note we only add 2 phantom cells on odd iface values
        // because 2 that would be added by its even "partner" iface value will be the same 
        if (level_bot != level_top) {
            if (level_bot < level_top) {
                cncell = bncell;
                fncell = tncell;
                lncell = ntop[cncell];
                rncell = nrht[lncell];
            }
            else {
                cncell = tncell;
                fncell = bncell;
                lncell = nbot[cncell];
                rncell = nrht[lncell];
            }
            if (fncell == lncell) { // left of the two vertical neighbors, add 4 phantoms
                idxVar = 0;
                //pcellCnt += 4;
                //pfaceCnt += 3;

                if (level_bot < level_top) { // top is more refined
                    // new face's adjacent cells
                    map_yface2cell_upper[pfaceIdx] = pcellIdx + 2;
                    map_yface2cell_lower[pfaceIdx] = bncell;        
                    map_yface2cell_upper[pfaceIdx+1] = pcellIdx;
                    map_yface2cell_lower[pfaceIdx+1] = pcellIdx+1;
                    map_yface2cell_upper[pfaceIdx+2] = pcellIdx+3;
                    map_yface2cell_lower[pfaceIdx+2] = pcellIdx+2;
                    // adjacent cells' face
                    map_ycell2face_top1[bncell] = pfaceIdx;
                    map_ycell2face_bot1[pcellIdx+2] = pfaceIdx;
                    map_ycell2face_top1[pcellIdx+1] = pfaceIdx+1;
                    map_ycell2face_bot1[pcellIdx] = pfaceIdx+1;
                    map_ycell2face_top1[pcellIdx+2] = pfaceIdx+2;
                    map_ycell2face_bot1[pcellIdx+3] = pfaceIdx+2;
                    map_ycell2face_top1[pcellIdx] = iface;
                    // "undo" the second face for the cell, as it is no longer applicable
                    map_ycell2face_top2[bncell] = -1;
                    // old face's new phantom adjacent cell
                    map_yface2cell_lower[iface] = pcellIdx;

                    // face positions
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    yface_level[pfaceIdx] = level[bncell];
                    yface_level[pfaceIdx+1] = level[tncell];
                    yface_level[pfaceIdx+2] = level[bncell];
                    yface_i[pfaceIdx] = i[bncell];
                    yface_i[pfaceIdx+1] = i[tncell];
                    yface_i[pfaceIdx+2] = i[bncell];
                    yface_j[pfaceIdx] = j[bncell] + 1;
                    yface_j[pfaceIdx+1] = j[tncell] - 1;
                    yface_j[pfaceIdx+2] = j[bncell] + 2;
    }

                    interpolate(2, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                    interpolate(6, pcellIdx, bncell, tncell, deltaT,  state_memory_old);

                    yrecvIdx[ifixupIdx] = pfaceIdx;
                    int nb = map_yface2cell_lower[pfaceIdx];
                    //printf("%d) %d %d\n", ifixupIdx, iface, pfaceIdx); 
                    yrecvCIdx[ifixupIdx] = nb;
                    yplusCell2Idx[nb] = ifixupIdx;
                    if (iface < map_ycell2face_bot1[rncell]) {
                       ysendIdx1[ifixupIdx] = iface;
                       ysendIdx2[ifixupIdx] = map_ycell2face_bot1[rncell];
                    } else {
                       ysendIdx1[ifixupIdx] = map_ycell2face_bot1[rncell];
                       ysendIdx2[ifixupIdx] = iface;
                    }
                    ifixupIdx++;


                }
                else { // bottom is more refined
                    // new face's adjacent cells
                    map_yface2cell_lower[pfaceIdx] = pcellIdx;        
                    map_yface2cell_upper[pfaceIdx] = tncell;
                    map_yface2cell_lower[pfaceIdx+1] = pcellIdx+1;
                    map_yface2cell_upper[pfaceIdx+1] = pcellIdx;
                    map_yface2cell_lower[pfaceIdx+2] = pcellIdx+2;
                    map_yface2cell_upper[pfaceIdx+2] = pcellIdx+3;
                    // adjacent cell's face
                    map_ycell2face_bot1[tncell] = pfaceIdx;
                    map_ycell2face_top1[pcellIdx] = pfaceIdx;
                    map_ycell2face_bot1[pcellIdx] = pfaceIdx+1;
                    map_ycell2face_top1[pcellIdx+1] = pfaceIdx+1;
                    map_ycell2face_bot1[pcellIdx+3] = pfaceIdx+2;
                    map_ycell2face_top1[pcellIdx+2] = pfaceIdx+2;
                    map_ycell2face_bot1[pcellIdx+2] = iface;
                    // "undo" the second face for the cell, as it is no longer applicable
                    map_ycell2face_bot2[tncell] = -1;
                    // old face's new phantom adjacent cell
                    map_yface2cell_upper[iface] = pcellIdx + 2;

                    // face positions
    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    yface_level[pfaceIdx] = level[tncell];
                    yface_level[pfaceIdx+1] = level[tncell];
                    yface_level[pfaceIdx+2] = level[bncell];
                    yface_i[pfaceIdx] = i[tncell];
                    yface_i[pfaceIdx+1] = i[tncell];
                    yface_i[pfaceIdx+2] = i[bncell];
                    yface_j[pfaceIdx] = j[tncell];
                    yface_j[pfaceIdx+1] = j[tncell] - 1;
                    yface_j[pfaceIdx+2] = j[bncell] + 2;
    }

                    interpolate(3, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                    interpolate(7, pcellIdx, bncell, tncell, deltaT,  state_memory_old);

                    yrecvIdx[ifixupIdx] = pfaceIdx;
                    int nt = map_yface2cell_upper[pfaceIdx];
                    //printf("%d) %d %d\n", ifixupIdx, iface, pfaceIdx); 
                    yrecvCIdx[ifixupIdx] = nt;
                    yminusCell2Idx[nt] = ifixupIdx;
                    if (iface < map_ycell2face_top1[rncell]){
                       ysendIdx1[ifixupIdx] = iface;
                       ysendIdx2[ifixupIdx] = map_ycell2face_top1[rncell];
                    } else {
                       ysendIdx1[ifixupIdx] = map_ycell2face_top1[rncell];
                       ysendIdx2[ifixupIdx] = iface;
                    }
                    ifixupIdx++;


                }
                // phantom cells' new neighbors (same regardless of which side refinement occurs)
                nbot[pcellIdx] = pcellIdx + 1;
                ntop[pcellIdx] = tncell;
                ntop[pcellIdx+1] = pcellIdx;
                nbot[pcellIdx+2] = bncell;
                ntop[pcellIdx+2] = pcellIdx + 3;
                nbot[pcellIdx+3] = pcellIdx + 2;

                // update other arrays
                i[pcellIdx] = i[tncell];
                i[pcellIdx+1] = i[tncell];
                i[pcellIdx+2] = i[bncell];
                i[pcellIdx+3] = i[bncell];
                j[pcellIdx] = j[tncell] - 1;
                j[pcellIdx+1] = j[tncell] - 2;
                j[pcellIdx+2] = j[bncell] + 1;
                j[pcellIdx+3] = j[bncell] + 2;
                level[pcellIdx] = level[tncell];
                level[pcellIdx+1] = level[tncell];
                level[pcellIdx+2] = level[bncell];
                level[pcellIdx+3] = level[bncell];
            }
            else { // right of the two vertical neighbors, only add 2 phantoms
                idxVar = 1;
                //pcellCnt += 2;
                //pfaceCnt ++;

                if (level[bncell] < level[tncell]) { // top is more refined
                    // old face's new phantom adjacent cell
                    map_yface2cell_lower[iface] = pcellIdx;
                    map_yface2cell_upper[pfaceIdx] = pcellIdx;
                    map_yface2cell_lower[pfaceIdx] = pcellIdx+1;

                    map_ycell2face_top1[pcellIdx+1] = pfaceIdx;
                    map_ycell2face_bot1[pcellIdx] = pfaceIdx;
                    map_ycell2face_top1[pcellIdx] = iface;
                    // phantom cells' new neighbors
                    nbot[pcellIdx] = pcellIdx + 1;
                    ntop[pcellIdx] = tncell;
                    ntop[pcellIdx+1] = pcellIdx;
                    // other arrays
                    i[pcellIdx] = i[tncell];
                    i[pcellIdx+1] = i[tncell];
                    j[pcellIdx] = j[tncell] - 1;
                    j[pcellIdx+1] = j[tncell] - 2;
                    level[pcellIdx] = level[tncell];
                    level[pcellIdx+1] = level[tncell];

    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    yface_level[pfaceIdx] = level[tncell];
                    yface_i[pfaceIdx] = i[tncell];
                    yface_j[pfaceIdx] = j[tncell] - 1;
    }

                    interpolate(2, pcellIdx, bncell, tncell, deltaT,  state_memory_old);

                }
                else { // bottom is more refined 
                    // old face's new phantom adjacent cell
                    map_yface2cell_upper[iface] = pcellIdx; 
                    map_yface2cell_upper[pfaceIdx] = pcellIdx+1;
                    map_yface2cell_lower[pfaceIdx] = pcellIdx;

                    map_ycell2face_top1[pcellIdx] = pfaceIdx;
                    map_ycell2face_bot1[pcellIdx+1] = pfaceIdx;
                    map_ycell2face_bot1[pcellIdx] = iface;
                    // phantom cells' new neighbors
                    nbot[pcellIdx] = bncell;
                    ntop[pcellIdx] = pcellIdx + 1;
                    nbot[pcellIdx+1] = pcellIdx;
                    // other arrays
                    i[pcellIdx] = i[bncell];
                    i[pcellIdx+1] = i[bncell];
                    j[pcellIdx] = j[bncell] + 1;
                    j[pcellIdx+1] = j[bncell] + 2;
                    level[pcellIdx] = level[bncell];
                    level[pcellIdx+1] = level[bncell];

    if (choose_amr_method == REGULAR_GRID_AMR || choose_amr_method == REGULAR_GRID_BY_FACES_AMR) {
                    yface_level[pfaceIdx] = level[bncell];
                    yface_i[pfaceIdx] = i[bncell];
                    yface_j[pfaceIdx] = j[bncell] + 2;
    }

                    interpolate(3, pcellIdx-2, bncell, tncell, deltaT,  state_memory_old);

                }

            }
#ifdef PATTERN_CHECK
            ycase[pfaceIdx] = 0;
#endif

            // update indexes
            //locpcellIdx += 4 - (idxVar % 2) * 2;
            pcellIdx += 4 - (idxVar % 2) * 2;
            pfaceIdx += 3 - (idxVar % 2) * 2;
        }

    }

   nyfixup = ifixupIdx;

    //if (ncells > 280)
     //   printf("%d %d\n", nlft[286], nrht[286]);
    // resize arrays/vectors with the addition of the new faces and cells
    
   cpu_timers[MESH_TIMER_BIDIRPART11] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   //printf("%d %d\n", pcellCnt, pcellIdx);
    /*map_yface2cell_lower.resize(pfaceIdx);
    map_yface2cell_upper.resize(pfaceIdx);
    map_ycell2face_bot1.resize(pcellIdx);
    map_ycell2face_top1.resize(pcellIdx);
    map_ycell2face_bot2.resize(pcellIdx);
    map_ycell2face_top2.resize(pcellIdx);
    yface_i.resize(pfaceIdx);
    yface_j.resize(pfaceIdx);
    yface_level.resize(pfaceIdx);*/

    ncells_phan = pcellCnt;

    // resize cell based  arrays (i, j, nlft, nrht, nbot, ntop) 
    // do this based on pcellIdx, as it was continuous w/ x and y faces 

    /*i        = (int *)mesh_memory.memory_realloc(ncells_phan, i);
    j        = (int *)mesh_memory.memory_realloc(ncells_phan, j);
    level    = (uchar_t *)mesh_memory.memory_realloc(ncells_phan, level);
    nlft     = (int *)mesh_memory.memory_realloc(ncells_phan, nlft);
    nrht     = (int *)mesh_memory.memory_realloc(ncells_phan, nrht);
    nbot     = (int *)mesh_memory.memory_realloc(ncells_phan, nbot);
    ntop     = (int *)mesh_memory.memory_realloc(ncells_phan, ntop);
    memory_reset_ptrs();*/
    //printf("\n%d\n", mesh_memory.get_memory_size(level));
    
   cpu_timers[MESH_TIMER_BIDIRPART12] += cpu_timer_stop(tstart_cpu_part);

#ifdef PATTERN_CHECK
   for (int ii=0; ii<256; ii++){
       xcase_count[ii]=0;
       sprintf(xcase_descrip[ii],"\0");
   }
   int bitsize = 4;

   xcase = (int *)malloc(nxface*sizeof(int));
   for (int iface = 0; iface < nxface; iface++){
      xcase[iface]=-1;
      int nl = map_xface2cell_lower[iface];
      int nll = nlft[nl];
      int nr = map_xface2cell_upper[iface];
      int nrr = nrht[nr];
      uchar_t ll = level[nl];
      uchar_t lll = level[nll];
      uchar_t lr = level[nr];
      uchar_t lrr = level[nrr];

      int imin = MIN4(lll,ll,lr,lrr);
      ll  -= imin;
      lll -= imin;
      lr  -= imin;
      lrr -= imin;
      char binlevstring[6];
      sprintf(binlevstring,"%1d%1d%1d%1db\0",ll,lr,lll,lrr);
      xcase[iface]=strtol(binlevstring,NULL,bitsize);
      sprintf(xcase_descrip[xcase[iface]],"   %s     %d     %d   %d     %d\0",binlevstring,lll,ll,lr,lrr);
      xcase_count[xcase[iface]]++;
      printf("%d %s %d %d %d %d %d\n",iface,binlevstring,lll,ll,lr,lrr,xcase[iface]);
   }

   for (int ll = 0; ll < 4; ll++){
      for (int lll = 0; lll < 4; lll++){
         for (int lr = 0; lr < 4; lr++){
            for (int lrr = 0; lrr < 4; lrr++){
               char binlevstring[6];
               sprintf(binlevstring,"%1d%1d%1d%1db\0",ll,lr,lll,lrr);
               int icase=strtol(binlevstring,NULL,bitsize);
               switch (icase) {
                   case 0 :
                   case 1 :
                   case 4 :
                   case 5 :
                   case 17 :
                   case 18 :
                   case 68 :
                   case 72 :
                   case 81 :
                   case 82 :
                   case 84 :
                   case 88 :
                   case 98 :
                   case 99 :
		           case 152 :
		           case 156 :
                      sprintf(xcase_descrip[icase],"   %s     %d     %d   %d     %d\0",binlevstring,lll,ll,lr,lrr);
                      break;
                   default: 
                      //printf("   %s     %d     %d   %d     %d\n",binlevstring,lll,ll,lr,lrr);
                      //printf("Face %d does not fit a case\n",iface);
                      break;
               }
            }
         }
      }
   }

   printf("          l,r,ll,rr        levels\n");
   printf("case count bincode lft-lft lft rht rht-rht\n");
   for (int ii=0; ii<256; ii++){
      if (xcase_count[ii] > 0 || strlen(xcase_descrip[ii]) > 1) {
         printf(" %3d %3d %s\n",ii,xcase_count[ii],xcase_descrip[ii]);
      }
   }

   // Check to see if we missed any cases
   for (int iface = 0; iface < nxface; iface++){
      switch (xcase[iface]) {
          case 0 :
          case 1 :
          case 4 :
          case 5 :
          case 17 :
          case 18 :
          case 68 :
          case 72 :
          case 81 :
          case 82 :
          case 84 :
          case 88 :
          case 98 :
          case 99 :
          case 152 :
          case 156 :
             break;
          default : 
             printf("Face %d does not fit a case\n",iface);
             break;
      }
   }
#endif
   cpu_timers[MESH_TIMER_BIDIR] += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_face_list_fill_phantom(MallocPlus &state_memory, double deltaT)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

    MallocPlus state_memory_old = state_memory;
    //malloc_plus_memory_entry *memory_item;

   cpu_timer_start(&tstart_cpu_part);

    int pcellIdx = (int) ncells;
    int idxVar = 0;

    for (int iface = 0; iface < nxface; iface++) {
        int lncell = map_xface2cell_lower[iface], // cell neighbor to the left
            rncell = map_xface2cell_upper[iface], // cell neighbor to the right
            cncell, // course neighbor to the face
            fncell, // fine neighbor to the face
            bncell, // bottom of the two fine neighbors of the course
            tncell, // top of the two fince neighbots of the course
            level_left = level[lncell],
            level_right = level[rncell];
        if (nlft[lncell] == lncell || nbot[lncell] == lncell || ntop[lncell] == lncell) continue;
        if (nrht[rncell] == rncell || nbot[rncell] == rncell || ntop[rncell] == rncell) continue;

        if (level_left != level_right) {
            if (level_left < level_right) {
                cncell = lncell;
                fncell = rncell;
                bncell = nrht[cncell];
                tncell = ntop[bncell];
            }
            else {
                cncell = rncell;
                fncell = lncell;
                bncell = nlft[cncell];
                tncell = ntop[bncell];
            }
            if (fncell == bncell) { // bottom of the two horizontal neighbors, add 4 phantoms
                idxVar = 0;

                if (level_left < level_right) { // right is more refined
                    interpolate(0, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                    interpolate(4, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                }
                else { // left is more refined
                    interpolate(1, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                    interpolate(5, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                }
            }
            else { // top of the two horizontal neighbors
                idxVar = 1;

                if (level[lncell] < level[rncell]) { // right is more refined
                    interpolate(0, pcellIdx, lncell, rncell, deltaT,  state_memory_old);
                }
                else { // left is more refined
                    interpolate(1, pcellIdx-2, lncell, rncell, deltaT,  state_memory_old);
                }
            }
            pcellIdx += 4 - (idxVar % 2) * 2;
        }
    }

    for (int iface = 0; iface < nyface; iface++) {
     int    bncell = map_yface2cell_lower[iface], // cell neighbor below
            tncell = map_yface2cell_upper[iface], // cell neighbor above
            cncell, // course neighbor to the face
            fncell, // fine neighbor to the face
            lncell, // bottom of the two fine neighbors of the course
            rncell, // top of the two fince neighbots of the course
            level_bot = level[bncell],
            level_top = level[tncell];
        if (nlft[bncell] == bncell || nbot[bncell] == bncell || nrht[bncell] == bncell) continue;
        if (nrht[tncell] == tncell || nrht[tncell] == tncell || ntop[tncell] == tncell) continue;

        if (level_bot != level_top) {
            if (level_bot < level_top) {
                cncell = bncell;
                fncell = tncell;
                lncell = ntop[cncell];
                rncell = nrht[lncell];
            }
            else {
                cncell = tncell;
                fncell = bncell;
                lncell = nbot[cncell];
                rncell = nrht[lncell];
            }
            if (fncell == lncell) { // left of the two vertical neighbors
                idxVar = 0;

                if (level_bot < level_top) { // top is more refined
                    interpolate(2, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                    interpolate(6, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                }
                else { // bot is more refined
                    interpolate(3, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                    interpolate(7, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                }
            }
            else { // right of the two vertical neighbors
                idxVar = 1;

                if (level[bncell] < level[tncell]) { // top is more refined
                    interpolate(2, pcellIdx, bncell, tncell, deltaT,  state_memory_old);
                }
                else { // bot is more refined
                    interpolate(3, pcellIdx-2, bncell, tncell, deltaT,  state_memory_old);
                }
            }
            pcellIdx += 4 - (idxVar % 2) * 2;
        }
    }

   cpu_timers[MESH_TIMER_BIDIRPART12] += cpu_timer_stop(tstart_cpu_part);

   cpu_timers[MESH_TIMER_BIDIR] += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_face_list_wbidirmap(void)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

    int xfaceCnt = 0;
    int yfaceCnt = 0;
    for (int nz = 0; nz < (int)ncells; nz++) {
        int nl = nlft[nz];
        int nr = nrht[nz];
        int nb = nbot[nz];
        int nt = ntop[nz];
#ifdef HAVE_MPI
        if (nl >= ncells) 
            xfaceCnt ++;
        if (nb >= ncells)
            yfaceCnt ++;
#endif
        if (nz == nl || nz == nr || nz == nb || nz == nt) {
            if (nz != nr) xfaceCnt ++;
            if (nz != nt) yfaceCnt ++;
            continue;
        }

        int lev = level[nz];
        //int ll = level[nl];
        int lr = level[nr];
        //int lb = level[nb];
        int lt = level[nt];
        if (lev < lr) {
            xfaceCnt += 2;
        }
        else {
            xfaceCnt ++;
        }
        if (lev < lt) {
            yfaceCnt += 2;
        }
        else {
            yfaceCnt ++;
        }
    }
    nxface = xfaceCnt;
    nyface = yfaceCnt;

    // realloc memory based on new counts

    int flags=0;
    if (nxface > (int)mesh_memory.get_memory_size(map_xface2cell_lower) ) {
       mesh_memory.memory_delete(map_xface2cell_lower);
       map_xface2cell_lower = (int *)mesh_memory.memory_malloc(nxface, sizeof(int), "map_xface2cell_lower", flags);
       mesh_memory.memory_delete(map_xface2cell_upper);
       map_xface2cell_upper = (int *)mesh_memory.memory_malloc(nxface, sizeof(int), "map_xface2cell_upper", flags);
       //mesh_memory.memory_delete(xface_i);
       //xface_i = (int *)mesh_memory.memory_malloc(nxface, sizeof(int), "xface_i", flags);
       //mesh_memory.memory_delete(xface_j);
       //xface_j = (int *)mesh_memory.memory_malloc(nxface, sizeof(int), "xface_j", flags);
       //mesh_memory.memory_delete(xface_level);
       //xface_level = (uchar_t *)mesh_memory.memory_malloc(nxface, sizeof(uchar_t), "xface_level", flags);
    }
    if (nyface > (int)mesh_memory.get_memory_size(map_yface2cell_lower) ) {
       mesh_memory.memory_delete(map_yface2cell_lower);
       map_yface2cell_lower = (int *)mesh_memory.memory_malloc(nyface, sizeof(int), "map_yface2cell_lower", flags);
       mesh_memory.memory_delete(map_yface2cell_upper);
       map_yface2cell_upper = (int *)mesh_memory.memory_malloc(nyface, sizeof(int), "map_yface2cell_upper", flags);
       mesh_memory.memory_delete(yface_i);
       //yface_i = (int *)mesh_memory.memory_malloc(nyface, sizeof(int), "yface_i", flags);
       //mesh_memory.memory_delete(yface_j);
       //yface_j = (int *)mesh_memory.memory_malloc(nyface, sizeof(int), "yface_j", flags);
       //mesh_memory.memory_delete(yface_level);
       //yface_level = (uchar_t *)mesh_memory.memory_malloc(nyface, sizeof(uchar_t), "yface_level", flags);
    }

    if (ncells_ghost > mesh_memory.get_memory_size(map_xcell2face_left1) ){
       mesh_memory.memory_delete(map_xcell2face_left1);
       map_xcell2face_left1 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_xcell2face_left1", flags);
       mesh_memory.memory_delete(map_xcell2face_left2);
       map_xcell2face_left2 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_xcell2face_left2", flags);
       mesh_memory.memory_delete(map_xcell2face_right1);
       map_xcell2face_right1 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_xcell2face_right1", flags);
       mesh_memory.memory_delete(map_xcell2face_right2);
       map_xcell2face_right2 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_xcell2face_right2", flags);
       mesh_memory.memory_delete(map_ycell2face_bot1);
       map_ycell2face_bot1 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_ycell2face_bot1", flags);
       mesh_memory.memory_delete(map_ycell2face_bot2);
       map_ycell2face_bot2 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_ycell2face_bot2", flags);
       mesh_memory.memory_delete(map_ycell2face_top1);
       map_ycell2face_top1 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_ycell2face_top1", flags);
       mesh_memory.memory_delete(map_ycell2face_top2);
       map_ycell2face_top2 = (int *)mesh_memory.memory_malloc(ncells_ghost, sizeof(int), "map_ycell2face_top2", flags);
    }

    memory_reset_ptrs();

/*
   map_xface2cell_lower.clear();
   map_xface2cell_upper.clear();

   map_xcell2face_left1.clear();
   map_xcell2face_left2.clear();
   map_xcell2face_right1.clear();
   map_xcell2face_right2.clear();
   map_xcell2face_left1.resize(ncells, -1);
   map_xcell2face_left2.resize(ncells, -1);
   map_xcell2face_right1.resize(ncells, -1);
   map_xcell2face_right2.resize(ncells, -1);

   xface_i.clear();
   xface_j.clear();
   xface_level.clear();

   ixmin_level.clear();
   ixmax_level.clear();
   jxmin_level.clear();
   jxmax_level.clear();
   ixmin_level.resize(levmx+1,  9999999);
   ixmax_level.resize(levmx+1, -9999999);
   jxmin_level.resize(levmx+1,  9999999);
   jxmax_level.resize(levmx+1, -9999999);

   ixadjust.clear();
   ixadjust.resize(levmx+1);
   jxadjust.clear();
   jxadjust.resize(levmx+1);
*/
   cpu_timers[MESH_TIMER_BIDIRPART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   int iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nr = nrht[nz];
#ifdef HAVE_MPI
      int nl = nlft[nz];
      if (nl >= ncells) {
        map_xface2cell_lower[iface] = nl;
        map_xface2cell_upper[iface] = nz;
        map_xcell2face_left1[nz] = iface;
        map_xcell2face_right1[nl] = iface;
        iface++;
      }
#endif
      if (nr == nz || nr <= -1) continue;

      //int ifactor = 1;
      //if (level[nr] < level[nz]) ifactor = 2;

      // Have right face
      map_xface2cell_lower[iface] = nz;
      map_xface2cell_upper[iface] = nr;
      //xface_level[iface] = MAX(level[nz],level[nr]);
      //xface_i[iface] = i[nr]*ifactor;
      //if (level[nr] < level[nz] && is_upper(j[nz]) ) {
      //   xface_j[iface] = j[nr]*ifactor+1;
      //} else {
      //   xface_j[iface] = j[nr]*ifactor;
      //}
      map_xcell2face_right1[nz] = iface;

      //the right is a real cell, but I am left boundary
      if (nz == nlft[nz])
          map_xcell2face_left1[nz] = iface;

      iface++;

      if (level[nr] > level[nz] && is_lower(j[nr]) ){
         int ntr = ntop[nr];
         if (ntr != nr) {
            map_xface2cell_lower[iface] = nz;
            map_xface2cell_upper[iface] = ntr;
            //xface_level[iface] = MAX(level[nz],level[ntr]);
            //xface_i[iface] = i[ntr]*ifactor;
            //xface_j[iface] = j[ntr]*ifactor;
            map_xcell2face_right2[nz] = iface;

            iface++;
         }
      }
   }

   cpu_timers[MESH_TIMER_BIDIRPART2] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int nz=0; nz<(int)ncells; nz++){
      int nl = nlft[nz];
      if (nl == nz || nl <= -1) continue;

      if (level[nl] < level[nz] && is_upper(j[nz]) && nbot[nz] != nz){
         map_xcell2face_left1[nz] = map_xcell2face_right2[nl];


      } else {
         map_xcell2face_left1[nz] = map_xcell2face_right1[nl];

         if (level[nl] > level[nz]){
            map_xcell2face_left2[nz] = map_xcell2face_right1[ntop[nl]];
         }
      }

      //the left is a real cell, but I am right boundary
      if (nz == nrht[nz])
          map_xcell2face_right1[nz] = map_xcell2face_left1[nz];

   }
   cpu_timers[MESH_TIMER_BIDIRPART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   //
   cpu_timer_start(&tstart_cpu_part);


   iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nt = ntop[nz];
#ifdef HAVE_MPI
      int nb = nbot[nz];
      if (nb >= ncells) {
        map_yface2cell_lower[iface] = nb;
        map_yface2cell_upper[iface] = nz;
        map_ycell2face_bot1[nz] = iface;
        map_ycell2face_top1[nb] = iface;
        iface++;
      }
#endif
      if (nt == nz || nt <= -1) continue;

      //int ifactor = 1;
      //if (level[nt] < level[nz]) ifactor = 2;

      // Have top face
      map_yface2cell_lower[iface] = nz;
      map_yface2cell_upper[iface] = nt;
      //yface_level[iface] = MAX(level[nz],level[nt]);
      //yface_j[iface] = j[nt]*ifactor;
      //if (level[nt] < level[nz] && is_upper(i[nz]) ) {
      //   yface_i[iface] = i[nt]*ifactor+1;
      //} else{
      //   yface_i[iface] = i[nt]*ifactor;
      //}
      map_ycell2face_top1[nz] = iface;

      //the top is a real cell, but I am bot boundary
      if (nz == nbot[nz])
          map_ycell2face_bot1[nz] = iface;

      iface++;

      if (level[nt] > level[nz]  && is_lower(i[nt]) ){
         int nrt = nrht[nt];
         if (nrt != nt) {
            map_yface2cell_lower[iface] = nz;
            map_yface2cell_upper[iface] = nrt;
            //yface_level[iface] = MAX(level[nz],level[nrt]);
            //yface_j[iface] = j[nrt]*ifactor;
            //yface_i[iface] = i[nrt]*ifactor;
            map_ycell2face_top2[nz] = iface;

            iface++;
         }
      }
   }

   cpu_timers[MESH_TIMER_BIDIRPART5] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int nz=0; nz<(int)ncells; nz++){
      int nb = nbot[nz];
      if (nb == nz || nb <= -1) continue;
      //if (nz == 3) printf("%d %d\n", map_ycell2face_top1[nb], map_ycell2face_top2[nb]);

      if (level[nb] < level[nz] && is_upper(i[nz]) && nlft[nz] != nz){
         map_ycell2face_bot1[nz] = map_ycell2face_top2[nb];

      } else {
         map_ycell2face_bot1[nz] = map_ycell2face_top1[nb];

         if (level[nb] > level[nz]){
            map_ycell2face_bot2[nz] = map_ycell2face_top1[nrht[nb]];
         }
      }
      //the bot is a real cell, but I am top boundary
      if (nz == ntop[nz])
          map_ycell2face_top1[nz] = map_ycell2face_bot1[nz];
   }



   cpu_timers[MESH_TIMER_BIDIRPART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);
/*
   map_yface2cell_lower.clear();
   map_yface2cell_upper.clear();

   map_ycell2face_bot1.clear();
   map_ycell2face_bot2.clear();
   map_ycell2face_top1.clear();
   map_ycell2face_top2.clear();
   map_ycell2face_bot1.resize(ncells, -1);
   map_ycell2face_bot2.resize(ncells, -1);
   map_ycell2face_top1.resize(ncells, -1);
   map_ycell2face_top2.resize(ncells, -1);

   yface_i.clear();
   yface_j.clear();
   yface_level.clear();

   iymin_level.clear();
   iymax_level.clear();
   jymin_level.clear();
   jymax_level.clear();
   iymin_level.resize(levmx+1,  9999999);
   iymax_level.resize(levmx+1, -9999999);
   jymin_level.resize(levmx+1,  9999999);
   jymax_level.resize(levmx+1, -9999999);

   iyadjust.clear();
   iyadjust.resize(levmx+1);
   jyadjust.clear();
   jyadjust.resize(levmx+1);
*/
   cpu_timers[MESH_TIMER_BIDIRPART4] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   cpu_timers[MESH_TIMER_BIDIRPART6] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);
/*
   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];

      int fi = xface_i[iface];
      if (fi < ixmin_level[fl]) ixmin_level[fl] = fi;
      if (fi > ixmax_level[fl]) ixmax_level[fl] = fi;

      int fj = xface_j[iface];
      if (fj < jxmin_level[fl]) jxmin_level[fl] = fj;
      if (fj > jxmax_level[fl]) jxmax_level[fl] = fj;
   }

   cpu_timers[MESH_TIMER_BIDIRPART7] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];
      if (ixmax_level[fl] < ixmin_level[fl]) continue;

      xface_i[iface] -= ixmin_level[fl];
      xface_j[iface] -= jxmin_level[fl];
   }

   cpu_timers[MESH_TIMER_BIDIRPART8] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int fl = 0; fl <= levmx; fl++){
      ixadjust[fl] = ixmin_level[fl];
      jxadjust[fl] = jxmin_level[fl];
      ixmax_level[fl] -= ixmin_level[fl];;
      jxmax_level[fl] -= jxmin_level[fl];
      ixmin_level[fl] = 0;
      jxmin_level[fl] = 0;
   }

   cpu_timers[MESH_TIMER_BIDIRPART9] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];

      int fi = yface_i[iface];
      if (fi < iymin_level[fl]) iymin_level[fl] = fi;
      if (fi > iymax_level[fl]) iymax_level[fl] = fi;

      int fj = yface_j[iface];
      if (fj < jymin_level[fl]) jymin_level[fl] = fj;
      if (fj > jymax_level[fl]) jymax_level[fl] = fj;
   }

   cpu_timers[MESH_TIMER_BIDIRPART10] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];
      if (iymax_level[fl] < iymin_level[fl]) continue;

      yface_i[iface] -= iymin_level[fl];
      yface_j[iface] -= jymin_level[fl];
   }

   for (int fl = 0; fl <= levmx; fl++){
      iyadjust[fl] = iymin_level[fl];
      jyadjust[fl] = jymin_level[fl];
      iymax_level[fl] -= iymin_level[fl];;
      jymax_level[fl] -= jymin_level[fl];
      iymin_level[fl] = 0;
      jymin_level[fl] = 0;
   }
*/
   cpu_timers[MESH_TIMER_BIDIRPART11] += cpu_timer_stop(tstart_cpu_part);

#ifdef PATTERN_CHECK
   for (int ii=0; ii<255; ii++){
       xcase_count[ii]=0;
       sprintf(xcase_descrip[ii],"\0");
   }
   int bitsize = 4;

   xcase = (int *)malloc(nxface*sizeof(int));
   for (int iface = 0; iface < nxface; iface++){
      xcase[iface]=-1;
      int nl = map_xface2cell_lower[iface];
      int nll = nlft[nl];
      int nr = map_xface2cell_upper[iface];
      int nrr = nrht[nr];
      uchar_t ll = level[nl];
      uchar_t lll = level[nll];
      uchar_t lr = level[nr];
      uchar_t lrr = level[nrr];

      int imin = MIN4(lll,ll,lr,lrr);
      ll  -= imin;
      lll -= imin;
      lr  -= imin;
      lrr -= imin;
      char binlevstring[6];
      sprintf(binlevstring,"%1d%1d%1d%1db\0",ll,lr,lll,lrr);
      xcase[iface]=strtol(binlevstring,NULL,bitsize);
      sprintf(xcase_descrip[xcase[iface]],"   %s     %d     %d   %d     %d\0",binlevstring,lll,ll,lr,lrr);
      xcase_count[xcase[iface]]++;
      printf("%d %s %d %d %d %d %d\n",iface,binlevstring,lll,ll,lr,lrr,xcase[iface]);
   }

   for (int ll = 0; ll < 4; ll++){
      for (int lll = 0; lll < 4; lll++){
         for (int lr = 0; lr < 4; lr++){
            for (int lrr = 0; lrr < 4; lrr++){
               char binlevstring[6];
               sprintf(binlevstring,"%1d%1d%1d%1db\0",ll,lr,lll,lrr);
               int icase=strtol(binlevstring,NULL,bitsize);
               switch (icase) {
                   case 0 :
                   case 1 :
                   case 4 :
                   case 5 :
                   case 17 :
                   case 18 :
                   case 68 :
                   case 72 :
                   case 81 :
                   case 82 :
                   case 84 :
                   case 88 :
                   case 98 :
                   case 99 :
		   case 152 :
		   case 156 :
                      sprintf(xcase_descrip[icase],"   %s     %d     %d   %d     %d\0",binlevstring,lll,ll,lr,lrr);
                      break;
                   default: 
                      //printf("   %s     %d     %d   %d     %d\n",binlevstring,lll,ll,lr,lrr);
                      //printf("Face %d does not fit a case\n",iface);
                      break;
               }
            }
         }
      }
   }

   printf("          l,r,ll,rr        levels\n");
   printf("case count bincode lft-lft lft rht rht-rht\n");
   for (int ii=0; ii<256; ii++){
      if (xcase_count[ii] > 0 || strlen(xcase_descrip[ii]) > 1) {
         printf(" %3d %3d %s\n",ii,xcase_count[ii],xcase_descrip[ii]);
      }
   }

   // Check to see if we missed any cases
   for (int iface = 0; iface < nxface; iface++){
      switch (xcase[iface]) {
          case 0 :
          case 1 :
          case 4 :
          case 5 :
          case 17 :
          case 18 :
          case 68 :
          case 72 :
          case 81 :
          case 82 :
          case 84 :
          case 88 :
          case 98 :
          case 99 :
	  case 152 :
	  case 156 :
             break;
          default : 
             printf("Face %d does not fit a case\n",iface);
             break;
      }
   }
#endif

   cpu_timers[MESH_TIMER_BIDIR] += cpu_timer_stop(tstart_cpu);
}

#ifdef NOVEC
void Mesh::calc_face_list_wbidirmap_novec(void) 
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);
/*
   //map_xface2cell_lower.clear();
   //map_xface2cell_upper.clear();
   cl_map_xface2cell_lower = (int *) realloc(ncells * sizeof(int));
   cl_map_xface2cell_upper = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_left1 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_left2 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_right1 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_right2 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_bot1 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_bot2 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_top1 = (int *) realloc(ncells * sizeof(int));
   cl_map_xcell2face_top2 = (int *) realloc(ncells * sizeof(int));

   //map_xcell2face_left1.clear();
   //map_xcell2face_left2.clear();
   //map_xcell2face_right1.clear();
   //map_xcell2face_right2.clear();
   //map_xcell2face_left1.resize(ncells, -1);
   //map_xcell2face_left2.resize(ncells, -1);
   //map_xcell2face_right1.resize(ncells, -1);
   //map_xcell2face_right2.resize(ncells, -1);

   //xface_i.clear();
   //xface_j.clear();
   //xface_level.clear();

   //ixmin_level.clear();
   //ixmax_level.clear();
   //jxmin_level.clear();
   //jxmax_level.clear();
   //ixmin_level.resize(levmx+1,  9999999);
   //ixmax_level.resize(levmx+1, -9999999);
   //jxmin_level.resize(levmx+1,  9999999);
   //jxmax_level.resize(levmx+1, -9999999);

   //ixadjust.clear();
   //ixadjust.resize(levmx+1);
   //jxadjust.clear();
   //jxadjust.resize(levmx+1);

   cpu_timers[MESH_TIMER_BIDIRPART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   int iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nr = nrht[nz];
      if (nr == nz) continue;

      int ifactor = 1;
      if (level[nr] < level[nz]) ifactor = 2;

      // Have right face
      map_xface2cell_lower.push_back(nz);
      map_xface2cell_upper.push_back(nr);
      xface_level.push_back(MAX(level[nz],level[nr]));
      xface_i.push_back(i[nr]*ifactor);
      if (level[nr] < level[nz] && is_upper(j[nz]) ) {
         xface_j.push_back(j[nr]*ifactor+1);
      } else {
         xface_j.push_back(j[nr]*ifactor);
      }
      map_xcell2face_right1[nz] = iface;

      //the right is a real cell, but I am left boundary
      if (nz == nlft[nz])
          map_xcell2face_left1[nz] = iface;

      iface++;

      if (level[nr] > level[nz] && is_lower(j[nr]) ){
         int ntr = ntop[nr];
         if (ntr != nr) {
            map_xface2cell_lower.push_back(nz);
            map_xface2cell_upper.push_back(ntr);
            xface_level.push_back(MAX(level[nz],level[ntr]));
            xface_i.push_back(i[ntr]*ifactor);
            xface_j.push_back(j[ntr]*ifactor);
            map_xcell2face_right2[nz] = iface;

            iface++;
         }
      }
   }
   nxface=iface;

   cpu_timers[MESH_TIMER_BIDIRPART2] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int nz=0; nz<(int)ncells; nz++){
      int nl = nlft[nz];
      if (nl == nz) continue;

      if (level[nl] < level[nz] && is_upper(j[nz])){
         map_xcell2face_left1[nz] = map_xcell2face_right2[nl];
      } else {
         map_xcell2face_left1[nz] = map_xcell2face_right1[nl];

         //the left is a real cell, but I am right boundary
         if (nz == nrht[nz])
             map_xcell2face_right1[nz] = map_xcell2face_left1[nz];

         if (level[nl] > level[nz]){
            map_xcell2face_left2[nz] = map_xcell2face_right1[ntop[nl]];
         }
      }

   }

   cpu_timers[MESH_TIMER_BIDIRPART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   map_yface2cell_lower.clear();
   map_yface2cell_upper.clear();

   map_ycell2face_bot1.clear();
   map_ycell2face_bot2.clear();
   map_ycell2face_top1.clear();
   map_ycell2face_top2.clear();
   map_ycell2face_bot1.resize(ncells, -1);
   map_ycell2face_bot2.resize(ncells, -1);
   map_ycell2face_top1.resize(ncells, -1);
   map_ycell2face_top2.resize(ncells, -1);

   yface_i.clear();
   yface_j.clear();
   yface_level.clear();

   iymin_level.clear();
   iymax_level.clear();
   jymin_level.clear();
   jymax_level.clear();
   iymin_level.resize(levmx+1,  9999999);
   iymax_level.resize(levmx+1, -9999999);
   jymin_level.resize(levmx+1,  9999999);
   jymax_level.resize(levmx+1, -9999999);

   iyadjust.clear();
   iyadjust.resize(levmx+1);
   jyadjust.clear();
   jyadjust.resize(levmx+1);

   cpu_timers[MESH_TIMER_BIDIRPART4] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   iface=0;
   for (int nz=0; nz<(int)ncells; nz++){
      int nt = ntop[nz];
      if (nt == nz) continue;

      int ifactor = 1;
      if (level[nt] < level[nz]) ifactor = 2;

      // Have top face
      // printf("DEBUG -- iface %d lower nz %d upper nr %d\n",iface,nz,nt);
      map_yface2cell_lower.push_back(nz);
      map_yface2cell_upper.push_back(nt);
      yface_level.push_back(MAX(level[nz],level[nt]));
      yface_j.push_back(j[nt]*ifactor);
      if (level[nt] < level[nz] && is_upper(i[nz]) ) {
         yface_i.push_back(i[nt]*ifactor+1);
      } else{
         yface_i.push_back(i[nt]*ifactor);
      }
      map_ycell2face_top1[nz] = iface;

      //the top is a real cell, but I am bot boundary
      if (nz == nbot[nz])
          map_ycell2face_bot1[nz] = iface;

      iface++;

      if (level[nt] > level[nz]  && is_lower(i[nt]) ){
         int nrt = nrht[nt];
         if (nrt != nt) {
            map_yface2cell_lower.push_back(nz);
            map_yface2cell_upper.push_back(nrt);
            yface_level.push_back(MAX(level[nz],level[nrt]));
            yface_j.push_back(j[nrt]*ifactor);
            yface_i.push_back(i[nrt]*ifactor);
            map_ycell2face_top2[nz] = iface;

            iface++;
         }
      }
   }
   nyface=iface;

   cpu_timers[MESH_TIMER_BIDIRPART5] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int nz=0; nz<(int)ncells; nz++){
      int nb = nbot[nz];
      if (nb == nz) continue;

      if (level[nb] < level[nz] && is_upper(i[nz])){
         map_ycell2face_bot1[nz] = map_ycell2face_top2[nb];
      } else {
         map_ycell2face_bot1[nz] = map_ycell2face_top1[nb];

         //the bot is a real cell, but I am top boundary
         if (nz == ntop[nz])
             map_ycell2face_top1[nz] = map_ycell2face_bot1[nz];

         if (level[nb] > level[nz]){
            map_ycell2face_bot2[nz] = map_ycell2face_top1[nrht[nb]];
         }
      }

   }

   cpu_timers[MESH_TIMER_BIDIRPART6] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];

      int fi = xface_i[iface];
      if (fi < ixmin_level[fl]) ixmin_level[fl] = fi;
      if (fi > ixmax_level[fl]) ixmax_level[fl] = fi;

      int fj = xface_j[iface];
      if (fj < jxmin_level[fl]) jxmin_level[fl] = fj;
      if (fj > jxmax_level[fl]) jxmax_level[fl] = fj;
   }

   cpu_timers[MESH_TIMER_BIDIRPART7] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nxface; iface++){
      uchar_t fl = xface_level[iface];
      if (ixmax_level[fl] < ixmin_level[fl]) continue;

      xface_i[iface] -= ixmin_level[fl];
      xface_j[iface] -= jxmin_level[fl];
   }

   cpu_timers[MESH_TIMER_BIDIRPART8] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int fl = 0; fl <= levmx; fl++){
      ixadjust[fl] = ixmin_level[fl];
      jxadjust[fl] = jxmin_level[fl];
      ixmax_level[fl] -= ixmin_level[fl];;
      jxmax_level[fl] -= jxmin_level[fl];
      ixmin_level[fl] = 0;
      jxmin_level[fl] = 0;
   }

   cpu_timers[MESH_TIMER_BIDIRPART9] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];

      int fi = yface_i[iface];
      if (fi < iymin_level[fl]) iymin_level[fl] = fi;
      if (fi > iymax_level[fl]) iymax_level[fl] = fi;

      int fj = yface_j[iface];
      if (fj < jymin_level[fl]) jymin_level[fl] = fj;
      if (fj > jymax_level[fl]) jymax_level[fl] = fj;
   }

   cpu_timers[MESH_TIMER_BIDIRPART10] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int iface=0; iface < nyface; iface++){
      uchar_t fl = yface_level[iface];
      if (iymax_level[fl] < iymin_level[fl]) continue;

      yface_i[iface] -= iymin_level[fl];
      yface_j[iface] -= jymin_level[fl];
   }

   for (int fl = 0; fl <= levmx; fl++){
      iyadjust[fl] = iymin_level[fl];
      jyadjust[fl] = jymin_level[fl];
      iymax_level[fl] -= iymin_level[fl];;
      jymax_level[fl] -= jymin_level[fl];
      iymin_level[fl] = 0;
      jymin_level[fl] = 0;
   }

   cpu_timers[MESH_TIMER_BIDIRPART11] += cpu_timer_stop(tstart_cpu_part);

#ifdef PATTERN_CHECK
   for (int ii=0; ii<255; ii++){
       xcase_count[ii]=0;
       sprintf(xcase_descrip[ii],"\0");
   }
   int bitsize = 4;

   xcase = (int *)malloc(nxface*sizeof(int));
   for (int iface = 0; iface < nxface; iface++){
      xcase[iface]=-1;
      int nl = map_xface2cell_lower[iface];
      int nll = nlft[nl];
      int nr = map_xface2cell_upper[iface];
      int nrr = nrht[nr];
      uchar_t ll = level[nl];
      uchar_t lll = level[nll];
      uchar_t lr = level[nr];
      uchar_t lrr = level[nrr];

      int imin = MIN4(lll,ll,lr,lrr);
      ll  -= imin;
      lll -= imin;
      lr  -= imin;
      lrr -= imin;
      char binlevstring[6];
      sprintf(binlevstring,"%1d%1d%1d%1db\0",ll,lr,lll,lrr);
      xcase[iface]=strtol(binlevstring,NULL,bitsize);
      sprintf(xcase_descrip[xcase[iface]],"   %s     %d     %d   %d     %d\0",binlevstring,lll,ll,lr,lrr);
      xcase_count[xcase[iface]]++;
      printf("%d %s %d %d %d %d %d\n",iface,binlevstring,lll,ll,lr,lrr,xcase[iface]);
   }

   for (int ll = 0; ll < 4; ll++){
      for (int lll = 0; lll < 4; lll++){
         for (int lr = 0; lr < 4; lr++){
            for (int lrr = 0; lrr < 4; lrr++){
               char binlevstring[6];
               sprintf(binlevstring,"%1d%1d%1d%1db\0",ll,lr,lll,lrr);
               int icase=strtol(binlevstring,NULL,bitsize);
               switch (icase) {
                   case 0 :
                   case 1 :
                   case 4 :
                   case 5 :
                   case 17 :
                   case 18 :
                   case 68 :
                   case 72 :
                   case 81 :
                   case 82 :
                   case 84 :
                   case 88 :
                   case 98 :
                   case 99 :
		   case 152 :
		   case 156 :
                      sprintf(xcase_descrip[icase],"   %s     %d     %d   %d     %d\0",binlevstring,lll,ll,lr,lrr);
                      break;
                   default: 
                      //printf("   %s     %d     %d   %d     %d\n",binlevstring,lll,ll,lr,lrr);
                      //printf("Face %d does not fit a case\n",iface);
                      break;
               }
            }
         }
      }
   }

   printf("          l,r,ll,rr        levels\n");
   printf("case count bincode lft-lft lft rht rht-rht\n");
   for (int ii=0; ii<256; ii++){
      if (xcase_count[ii] > 0 || strlen(xcase_descrip[ii]) > 1) {
         printf(" %3d %3d %s\n",ii,xcase_count[ii],xcase_descrip[ii]);
      }
   }

   // Check to see if we missed any cases
   for (int iface = 0; iface < nxface; iface++){
      switch (xcase[iface]) {
          case 0 :
          case 1 :
          case 4 :
          case 5 :
          case 17 :
          case 18 :
          case 68 :
          case 72 :
          case 81 :
          case 82 :
          case 84 :
          case 88 :
          case 98 :
          case 99 :
	  case 152 :
	  case 156 :
             break;
          default : 
             printf("Face %d does not fit a case\n",iface);
             break;
      }
   }
#endif

   cpu_timers[MESH_TIMER_BIDIR] += cpu_timer_stop(tstart_cpu);
   */
}
#endif

#ifdef HAVE_OPENCL
void Mesh::gpu_wbidirmap_only_essentials(void)
{
}

void Mesh::gpu_wbidirmap_setup(void)
{
    size_t mem_request = ncells * 3;
    size_t mem_faces = 2;
    int flags = 0;
    dev_nface = ezcl_malloc(NULL, const_cast<char *>("dev_nface"), &mem_faces, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_xface2cell_lower = ezcl_malloc(NULL, const_cast<char *>("dev_map_xface2cell_lower"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_xface2cell_upper = ezcl_malloc(NULL, const_cast<char *>("dev_map_xface2cell_upper"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_xcell2face_left1 = ezcl_malloc(NULL, const_cast<char *>("dev_map_xcell2face_left1"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_xcell2face_left2 = ezcl_malloc(NULL, const_cast<char *>("dev_map_xcell2face_left2"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_xcell2face_right1 = ezcl_malloc(NULL, const_cast<char *>("dev_map_xcell2face_right1"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_xcell2face_right2 = ezcl_malloc(NULL, const_cast<char *>("dev_map_xcell2face_right2"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xface_level = ezcl_malloc(NULL, const_cast<char *>("dev_xface_level"), &mem_request, sizeof(cl_uchar_t), CL_MEM_READ_WRITE, 0);
    dev_xface_i = ezcl_malloc(NULL, const_cast<char *>("dev_xface_i"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xface_j = ezcl_malloc(NULL, const_cast<char *>("dev_xface_j"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_ixmin_level = ezcl_malloc(NULL, const_cast<char *>("dev_ixmin_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_ixmax_level = ezcl_malloc(NULL, const_cast<char *>("dev_ixmax_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_jxmin_level = ezcl_malloc(NULL, const_cast<char *>("dev_jxmin_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_jxmax_level = ezcl_malloc(NULL, const_cast<char *>("dev_jxmax_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_yface2cell_lower = ezcl_malloc(NULL, const_cast<char *>("dev_map_yface2cell_lower"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_yface2cell_upper = ezcl_malloc(NULL, const_cast<char *>("dev_map_yface2cell_upper"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_ycell2face_bot1 = ezcl_malloc(NULL, const_cast<char *>("dev_map_ycell2face_bot1"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_ycell2face_bot2 = ezcl_malloc(NULL, const_cast<char *>("dev_map_ycell2face_bot2"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_ycell2face_top1 = ezcl_malloc(NULL, const_cast<char *>("dev_map_ycell2face_top1"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_map_ycell2face_top2 = ezcl_malloc(NULL, const_cast<char *>("dev_map_ycell2face_top2"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yface_level = ezcl_malloc(NULL, const_cast<char *>("dev_yface_level"), &mem_request, sizeof(cl_uchar_t), CL_MEM_READ_WRITE, 0);
    dev_yface_i = ezcl_malloc(NULL, const_cast<char *>("dev_yface_i"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yface_j = ezcl_malloc(NULL, const_cast<char *>("dev_yface_j"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_iymin_level = ezcl_malloc(NULL, const_cast<char *>("dev_iymin_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_iymax_level = ezcl_malloc(NULL, const_cast<char *>("dev_iymax_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_jymin_level = ezcl_malloc(NULL, const_cast<char *>("dev_jymin_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_jymax_level = ezcl_malloc(NULL, const_cast<char *>("dev_jymax_level"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xfaceIdxList = ezcl_malloc(NULL, const_cast<char *>("dev_xfaceIdxList"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yfaceIdxList = ezcl_malloc(NULL, const_cast<char *>("dev_yfaceIdxList"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_pxcellCnt = ezcl_malloc(NULL, const_cast<char *>("dev_pxcellCnt"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_pycellCnt = ezcl_malloc(NULL, const_cast<char *>("dev_pycellCnt"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_pxfaceCnt = ezcl_malloc(NULL, const_cast<char *>("dev_pxfaceCnt"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_pyfaceCnt = ezcl_malloc(NULL, const_cast<char *>("dev_pyfaceCnt"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xrecvIdx = ezcl_malloc(NULL, const_cast<char *>("dev_xrecvIdx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xrecvCIdx = ezcl_malloc(NULL, const_cast<char *>("dev_xrecvCIdx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xplusCell2Idx = ezcl_malloc(NULL, const_cast<char *>("dev_xplusCell2Idx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xminusCell2Idx = ezcl_malloc(NULL, const_cast<char *>("dev_xminusCell2Idx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xsendIdx1 = ezcl_malloc(NULL, const_cast<char *>("dev_xsendIdx1"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_xsendIdx2 = ezcl_malloc(NULL, const_cast<char *>("dev_xsendIdx2"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yrecvIdx = ezcl_malloc(NULL, const_cast<char *>("dev_yrecvIdx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yrecvCIdx = ezcl_malloc(NULL, const_cast<char *>("dev_yrecvCIdx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yplusCell2Idx = ezcl_malloc(NULL, const_cast<char *>("dev_yplusCell2Idx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_yminusCell2Idx = ezcl_malloc(NULL, const_cast<char *>("dev_yminusCell2Idx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_ysendIdx1 = ezcl_malloc(NULL, const_cast<char *>("dev_ysendIdx1"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_ysendIdx2 = ezcl_malloc(NULL, const_cast<char *>("dev_ysendIdx2"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_ifixupXCnt = ezcl_malloc(NULL, const_cast<char *>("dev_ifixupXCnt"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_ifixupYCnt = ezcl_malloc(NULL, const_cast<char *>("dev_ifixupYCnt"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_pxcellIdx = ezcl_malloc(NULL, const_cast<char *>("dev_pxcellIdx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_pycellIdx = ezcl_malloc(NULL, const_cast<char *>("dev_pycellIdx"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

}

void Mesh::gpu_call_wbidirmap_realloc(void)
{
    gpu_wbidirmap_realloc(&dev_map_xface2cell_lower, pxfaceCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_xface2cell_upper, pxfaceCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_yface2cell_lower, pyfaceCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_yface2cell_upper, pyfaceCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_xcell2face_left1, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_xcell2face_left2, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_xcell2face_right1, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_xcell2face_right2, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_ycell2face_bot1, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_ycell2face_bot2, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_ycell2face_top1, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_map_ycell2face_top2, pcellCnt, ncells*3);
    gpu_wbidirmap_realloc(&dev_xface_level, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xface_i, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xface_j, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xfaceIdxList, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_ifixupXCnt, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_pxcellIdx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_pxcellCnt, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xrecvIdx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xrecvCIdx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xplusCell2Idx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xminusCell2Idx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xsendIdx1, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_xsendIdx2, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yface_level, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yface_i, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yface_j, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yfaceIdxList, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_ifixupYCnt, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_pycellIdx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_pycellCnt, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yrecvIdx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yrecvCIdx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yplusCell2Idx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_yminusCell2Idx, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_ysendIdx1, ncells, ncells*3);
    gpu_wbidirmap_realloc(&dev_ysendIdx2, ncells, ncells*3);
}

void Mesh::gpu_wbidirmap_delete(void)
{
   ezcl_device_memory_delete(dev_map_xcell2face_left1);
   ezcl_device_memory_delete(dev_map_xcell2face_right1);
   ezcl_device_memory_delete(dev_map_ycell2face_bot1);
   ezcl_device_memory_delete(dev_map_ycell2face_top1);
   ezcl_device_memory_delete(dev_nface);
   ezcl_device_memory_delete(dev_map_xface2cell_lower);
   ezcl_device_memory_delete(dev_map_xface2cell_upper);
   ezcl_device_memory_delete(dev_map_yface2cell_lower);
   ezcl_device_memory_delete(dev_map_yface2cell_upper);
   ezcl_device_memory_delete(dev_map_xcell2face_left2);
   ezcl_device_memory_delete(dev_map_xcell2face_right2);
   ezcl_device_memory_delete(dev_map_ycell2face_bot2);
   ezcl_device_memory_delete(dev_map_ycell2face_top2);
   ezcl_device_memory_delete(dev_xface_level);
   ezcl_device_memory_delete(dev_xface_i);
   ezcl_device_memory_delete(dev_xface_j);
   ezcl_device_memory_delete(dev_ixmin_level);
   ezcl_device_memory_delete(dev_ixmax_level);
   ezcl_device_memory_delete(dev_jxmin_level);
   ezcl_device_memory_delete(dev_jxmax_level);
   ezcl_device_memory_delete(dev_yface_level);
   ezcl_device_memory_delete(dev_yface_i);
   ezcl_device_memory_delete(dev_yface_j);
   ezcl_device_memory_delete(dev_iymin_level);
   ezcl_device_memory_delete(dev_iymax_level);
   ezcl_device_memory_delete(dev_jymin_level);
   ezcl_device_memory_delete(dev_jymax_level);
   ezcl_device_memory_delete(dev_xfaceIdxList);
   ezcl_device_memory_delete(dev_yfaceIdxList);
   // For phantoms
   ezcl_device_memory_delete(dev_pxcellCnt);
   ezcl_device_memory_delete(dev_pycellCnt);
   ezcl_device_memory_delete(dev_pxfaceCnt);
   ezcl_device_memory_delete(dev_pyfaceCnt);
   ezcl_device_memory_delete(dev_xrecvIdx);
   ezcl_device_memory_delete(dev_xrecvCIdx);
   ezcl_device_memory_delete(dev_xplusCell2Idx);
   ezcl_device_memory_delete(dev_xminusCell2Idx);
   ezcl_device_memory_delete(dev_xsendIdx1);
   ezcl_device_memory_delete(dev_xsendIdx2);
   ezcl_device_memory_delete(dev_yrecvIdx);
   ezcl_device_memory_delete(dev_yrecvCIdx);
   ezcl_device_memory_delete(dev_yplusCell2Idx);
   ezcl_device_memory_delete(dev_yminusCell2Idx);
   ezcl_device_memory_delete(dev_ysendIdx1);
   ezcl_device_memory_delete(dev_ysendIdx2);
   ezcl_device_memory_delete(dev_ifixupXCnt);
   ezcl_device_memory_delete(dev_ifixupYCnt);
   ezcl_device_memory_delete(dev_pxcellIdx);
   ezcl_device_memory_delete(dev_pycellIdx);
}

void Mesh::gpu_wbidirmap_realloc(cl_mem *dev_mem_ptr, int old_size, size_t mem_request)
{
        cl_command_queue command_queue = ezcl_get_command_queue();

        cl_mem dev_mem_ptr_new = ezcl_malloc(NULL, const_cast<char *>("dev_mem_ptr_new"), &mem_request, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
        cl_mem *dev_tmp;
        cl_mem dev_ptr= (cl_mem)*dev_mem_ptr;
        *dev_tmp = *dev_mem_ptr;


        size_t local_work_size = 128;
        size_t global_work_size = ((old_size+local_work_size - 1) /local_work_size) * local_work_size;
        //printf("%d %d\n", old_size, global_work_size);
        cl_event deep_copy_event;
        ezcl_set_kernel_arg(kernel_deep_copy, 0, sizeof(cl_int), (void *)&old_size);
        ezcl_set_kernel_arg(kernel_deep_copy, 1, sizeof(cl_mem), (void *)dev_mem_ptr);
        ezcl_set_kernel_arg(kernel_deep_copy, 2, sizeof(cl_mem), (void *)&dev_mem_ptr_new);

        ezcl_enqueue_ndrange_kernel(command_queue, kernel_deep_copy, 1, NULL, &global_work_size, &local_work_size, &deep_copy_event);

        ezcl_wait_for_events(1, &deep_copy_event);
        ezcl_event_release(deep_copy_event);
        //printf("map old %p\n", *dev_mem_ptr);
        //printf("Mesh old %p %p %p\n", dev_tmp, dev_tmp_new, dev_mem_ptr);
        *dev_mem_ptr = dev_mem_ptr_new;
        //printf("Mesh new %p %p %p\n", dev_tmp, dev_tmp_new, dev_mem_ptr);
        //printf("map new %p\n", *dev_mem_ptr);

        //printf("delete old %p\n", *dev_tmp);
        ezcl_device_memory_delete(*dev_tmp);
        //gpu_state_memory.memory_replace(dev_state_var, dev_state_var_new);
        
        //gpu_state_memory.memory_realloc(pcellCnt, dev_state_var);
        //dev_state_var_new = (cl_mem)memory_item->mem_ptr;
        //memory_item->mem_ptr = ezcl_device_memory_realloc(dev_state_var_old, pcellCnt);
        //cl_mem dev_state_var_new = (cl_mem)memory_item->mem_ptr;
        //////////////dev_state_var = (cl_mem)gpu_state_memory.memory_replace(dev_state_var, dev_state_var_new);
      //printf("DEBUG -- it.mem_name %s elsize %lu capacity %lu\n",memory_item->mem_name,memory_item->mem_elsize, memory_item->mem_capacity);
   //vector<real_t>H_loc(pcellCnt);
   //ezcl_enqueue_read_buffer(command_queue, dev_state_var_new,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
    //  for (int i = 0; i < pcellCnt; i++) {printf("%f\t", H_loc[i]); }
     // printf("\n");

    //}  
}

int Mesh::gpu_serial_int_reduce(int *arr, int count, int length)
{
    int sum = count;
    for (int i = 0; i < length; i++) {
        int my_val = arr[i];
        arr[i] = sum;
        sum += my_val;
    }
    return sum;
}

void Mesh::gpu_calc_face_list_wbidirmap(void)
{
    struct timespec tstart_cpu;
    cpu_timer_start(&tstart_cpu);
    
    cl_command_queue command_queue = ezcl_get_command_queue();

    assert(dev_nlft);
    assert(dev_nrht);
    assert(dev_nbot);
    assert(dev_ntop);
    assert(dev_level);
    assert(dev_levdx);
    assert(dev_levdy);
    assert(dev_map_xface2cell_lower);
    assert(dev_map_xface2cell_upper);
    assert(dev_map_xcell2face_left1);
    assert(dev_map_xcell2face_left2);
    assert(dev_map_xcell2face_right1);
    assert(dev_map_xcell2face_right2);
    assert(dev_map_yface2cell_lower);
    assert(dev_map_yface2cell_upper);
    assert(dev_map_ycell2face_bot1);
    assert(dev_map_ycell2face_bot2);
    assert(dev_map_ycell2face_top1);
    assert(dev_map_ycell2face_top2);

    //size_t local_work_size  = MIN(ncells, TILE_SIZE);
    size_t local_work_size  = 128;
    size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

    vector<int>this_xfaceIdxList(ncells);
    vector<int>this_yfaceIdxList(ncells);
    vector<int>this_nlft(ncells);
    vector<int>this_nrht(ncells);
    vector<int>this_nbot(ncells);
    vector<int>this_ntop(ncells);
    vector<int>this_level(ncells);
    ezcl_enqueue_read_buffer(command_queue, dev_nlft,     CL_TRUE, 0, ncells*sizeof(cl_int), &this_nlft[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_nrht,     CL_TRUE, 0, ncells*sizeof(cl_int), &this_nrht[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_nbot,     CL_TRUE, 0, ncells*sizeof(cl_int), &this_nbot[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_ntop,     CL_TRUE, 0, ncells*sizeof(cl_int), &this_ntop[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_level,     CL_TRUE, 0, ncells*sizeof(cl_int), &this_level[0], NULL);


    this_xfaceIdxList[0] = 0;
    this_yfaceIdxList[0] = 0;
    for (int ccc = 0; ccc < (int) ncells - 1; ccc++) {
        int xfaceSt = this_xfaceIdxList[ccc], yfaceSt = this_yfaceIdxList[ccc];
        if (this_nrht[ccc] != ccc) {
            xfaceSt++;
            //nxfaceCnt++;
        }
        if (this_ntop[ccc] != ccc) {
            yfaceSt++;
            //nyfaceCnt++;
        }
        if ((this_level[this_nrht[ccc]] > this_level[ccc]) && (this_ntop[this_nrht[ccc]] != this_nrht[ccc]) && (this_nbot[this_nrht[ccc]] != this_nrht[ccc])) {
            xfaceSt++; 
            //nxfaceCnt++;
        }
        if ((this_level[this_ntop[ccc]] > this_level[ccc]) && (this_nrht[this_ntop[ccc]] != this_ntop[ccc]) && (this_nlft[this_ntop[ccc]] != this_ntop[ccc])) {
            yfaceSt++;
            //nyfaceCnt++;
        }
        this_xfaceIdxList[ccc+1] = xfaceSt;
        this_yfaceIdxList[ccc+1] = yfaceSt;
    }

    ezcl_enqueue_write_buffer(command_queue, dev_xfaceIdxList, CL_TRUE, 0, ncells*sizeof(cl_int), &this_xfaceIdxList[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_yfaceIdxList, CL_TRUE, 0, ncells*sizeof(cl_int), &this_yfaceIdxList[0], NULL);

/*
__kernel void face_idx_wbidirmap_cl(
            __global          int   *level,                     // 0
            __global          int   *nlft,                      // 1
            __global          int   *nrht,                      // 2
            __global          int   *nbot,                      // 3
            __global          int   *ntop,                      // 4
            __global          int   *xfaceIdxList,              // 5
            __global          int   *yfaceIdxList,              // 6
                              int   ncells,                     // 7
            __global          int   *nface) {                   // 8
*/
/*
    cl_event face_idx_event;
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 0, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 1, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 2, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 3, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 4, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 5, sizeof(cl_mem), (void *)&dev_xfaceIdxList);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 6, sizeof(cl_mem), (void *)&dev_yfaceIdxList);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 7, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 8, sizeof(cl_mem), (void *)&dev_nface);
    size_t single_work = 1;
    
    ezcl_enqueue_ndrange_kernel(command_queue, kernel_face_idx_wbidirmap, 1, NULL, &single_work, &single_work, &face_idx_event);

    ezcl_wait_for_events(1, &face_idx_event);
    ezcl_event_release(face_idx_event);
*/
/*
__kernel void calc_face_list_wbidirmap_pt1_cl(
                        const int   ncells,                     // 0
            __global          int   *nface,                     // 1
                        const int   levmx,                      // 2
            __global          int   *level,                     // 3
            __global          int   *nlft,                      // 4
            __global          int   *nrht,                      // 5
            __global          int   *nbot,                      // 6
            __global          int   *ntop,                      // 7
            __global          int   *map_xface2cell_lower,      // 8
            __global          int   *map_xface2cell_upper,      // 9
            __global          int   *map_xcell2face_left1,      // 10
            __global          int   *map_xcell2face_left2,      // 11
            __global          int   *map_xcell2face_right1,     // 12
            __global          int   *map_xcell2face_right2,     // 13
            __global          int   *xface_level,               // 14
            __global          int   *xface_i,                   // 15
            __global          int   *xface_j,                   // 16
            __global          int   *ixmin_level,               // 17
            __global          int   *ixmax_level,               // 18
            __global          int   *jxmin_level,               // 19
            __global          int   *jxmax_level,               // 20
            __global          int   *map_yface2cell_lower,      // 21
            __global          int   *map_yface2cell_upper,      // 22
            __global          int   *map_ycell2face_bot1,       // 23
            __global          int   *map_ycell2face_bot2,       // 24
            __global          int   *map_ycell2face_top1,       // 25
            __global          int   *map_ycell2face_top2,       // 26
            __global          int   *yface_level,               // 27
            __global          int   *yface_i,                   // 28
            __global          int   *yface_j,                   // 29
            __global          int   *iymin_level,               // 30
            __global          int   *iymax_level,               // 31
            __global          int   *jymin_level,               // 32
            __global          int   *jymax_level,               // 33
            __global          int   *xfaceIdxList,              // 34
            __global          int   *yfaceIdxList,              // 35
            __global    const int   *i,                         // 36
            __global    const int   *j) {                       // 37
            */


    cl_event calc_face_list_wbidirmap_event1, calc_face_list_wbidirmap_event2;
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 0, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 1, sizeof(cl_mem), (void *)&dev_nface);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 2, sizeof(cl_int), (void *)&levmx);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 3, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 4, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 5, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 6, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 7, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 8, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 9, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 10, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 11, sizeof(cl_mem), (void *)&dev_map_xcell2face_left2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 12, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 13, sizeof(cl_mem), (void *)&dev_map_xcell2face_right2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 14, sizeof(cl_mem), (void *)&dev_xface_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 15, sizeof(cl_mem), (void *)&dev_xface_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 16, sizeof(cl_mem), (void *)&dev_xface_j);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 17, sizeof(cl_mem), (void *)&dev_ixmin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 18, sizeof(cl_mem), (void *)&dev_ixmax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 19, sizeof(cl_mem), (void *)&dev_jxmin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 20, sizeof(cl_mem), (void *)&dev_jxmax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 21, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 22, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 23, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 24, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 25, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 26, sizeof(cl_mem), (void *)&dev_map_ycell2face_top2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 27, sizeof(cl_mem), (void *)&dev_yface_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 28, sizeof(cl_mem), (void *)&dev_yface_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 29, sizeof(cl_mem), (void *)&dev_yface_j);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 30, sizeof(cl_mem), (void *)&dev_iymin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 31, sizeof(cl_mem), (void *)&dev_iymax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 32, sizeof(cl_mem), (void *)&dev_jymin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 33, sizeof(cl_mem), (void *)&dev_jymax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 34, sizeof(cl_mem), (void *)&dev_xfaceIdxList);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 35, sizeof(cl_mem), (void *)&dev_yfaceIdxList);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 36, sizeof(cl_mem), (void *)&dev_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 37, sizeof(cl_mem), (void *)&dev_j);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_face_list_wbidirmap_pt1, 1, NULL, &global_work_size, &local_work_size, &calc_face_list_wbidirmap_event1);

    ezcl_wait_for_events(1, &calc_face_list_wbidirmap_event1);
    ezcl_event_release(calc_face_list_wbidirmap_event1);

    /*
__kernel void calc_face_list_wbidirmap_pt2_cl(
                        const int   ncells,                     // 0
            __global          int   *level,                     // 1
            __global          int   *nlft,                      // 2
            __global          int   *nrht,                      // 3
            __global          int   *nbot,                      // 4
            __global          int   *ntop,                      // 5
            __global          int   *map_xcell2face_left1,      // 6
            __global          int   *map_xcell2face_left2,      // 7
            __global          int   *map_xcell2face_right1,     // 8
            __global          int   *map_xcell2face_right2,     // 9
            __global          int   *map_ycell2face_bot1,       // 10
            __global          int   *map_ycell2face_bot2,       // 11
            __global          int   *map_ycell2face_top1,       // 12
            __global          int   *map_ycell2face_top2,       // 13
            __global    const int   *i,                         // 14
            __global    const int   *j) {                       // 15
    */

    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 0, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 1, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 2, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 3, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 4, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 5, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 6, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 7, sizeof(cl_mem), (void *)&dev_map_xcell2face_left2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 8, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 9, sizeof(cl_mem), (void *)&dev_map_xcell2face_right2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 10, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 11, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 12, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 13, sizeof(cl_mem), (void *)&dev_map_ycell2face_top2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 14, sizeof(cl_mem), (void *)&dev_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 15, sizeof(cl_mem), (void *)&dev_j);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_face_list_wbidirmap_pt2, 1, NULL, &global_work_size, &local_work_size, &calc_face_list_wbidirmap_event2);

    ezcl_wait_for_events(1, &calc_face_list_wbidirmap_event2);
    ezcl_event_release(calc_face_list_wbidirmap_event2);

   /*
   vector<int>map_xcell_tmp_l(ncells);
   vector<int>map_xcell_tmp_r(ncells);
   vector<int>map_ycell_tmp_b(ncells);
   vector<int>map_ycell_tmp_t(ncells);
   vector<int>l_tmp(ncells);
   vector<int>r_tmp(ncells);
   vector<int>b_tmp(ncells);
   vector<int>t_tmp(ncells);
   vector<int>facelistX(ncells);
   vector<int>facelistY(ncells);
   vector<int>faceIdx(ncells);
   int nxface;
   ezcl_enqueue_read_buffer(command_queue, dev_map_xcell2face_left1,     CL_TRUE, 0, ncells*sizeof(cl_int), &map_xcell_tmp_l[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xcell2face_right1,     CL_TRUE, 0, ncells*sizeof(cl_int), &map_xcell_tmp_r[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_ycell2face_bot1,     CL_TRUE, 0, ncells*sizeof(cl_int), &map_ycell_tmp_b[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_ycell2face_top1,     CL_TRUE, 0, ncells*sizeof(cl_int), &map_ycell_tmp_t[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft,     CL_TRUE, 0, ncells*sizeof(cl_int), &l_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht,     CL_TRUE, 0, ncells*sizeof(cl_int), &r_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot,     CL_TRUE, 0, ncells*sizeof(cl_int), &b_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop,     CL_TRUE, 0, ncells*sizeof(cl_int), &t_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_xfaceIdxList,     CL_TRUE, 0, ncells*sizeof(cl_int), &facelistX[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_yfaceIdxList,     CL_TRUE, 0, ncells*sizeof(cl_int), &facelistY[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nxface, NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %d %d %d %d\n", jello, map_xcell_tmp_l[jello], map_xcell_tmp_r[jello], map_ycell_tmp_b[jello], map_ycell_tmp_t[jello]); if (jello == 40) printf("\n|%d %d %d %d|\n", l_tmp[40], r_tmp[40], b_tmp[40], t_tmp[40]); }
   //for (int jello = 0; jello < ncells; jello++) { printf("%d) %d %d\n", jello, facelistX[jello], facelistX[jello]); }
   */

    gpu_timers[MESH_TIMER_BIDIR] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
}

void Mesh::gpu_calc_face_list_wbidirmap_phantom(MallocPlus &gpu_state_memory, double deltaT)
{

    /*if (!gpu_do_rezone && !firstFlag) {
        printf("hello\n");
          gpu_calc_face_list_phantom_fill(gpu_state_memory, deltaT);
          return;
    }*/

    //gpu_call_wbidirmap_realloc();

    struct timespec tstart_cpu;
    cpu_timer_start(&tstart_cpu);

    //struct timespec tstart_cpu_part;
    //cpu_timer_start(&tstart_cpu_part);
    
    cl_command_queue command_queue = ezcl_get_command_queue();

    assert(dev_nlft);
    assert(dev_nrht);
    assert(dev_nbot);
    assert(dev_ntop);
    assert(dev_level);
    assert(dev_levdx);
    assert(dev_levdy);
    assert(dev_map_xface2cell_lower);
    assert(dev_map_xface2cell_upper);
    assert(dev_map_xcell2face_left1);
    assert(dev_map_xcell2face_left2);
    assert(dev_map_xcell2face_right1);
    assert(dev_map_xcell2face_right2);
    assert(dev_map_yface2cell_lower);
    assert(dev_map_yface2cell_upper);
    assert(dev_map_ycell2face_bot1);
    assert(dev_map_ycell2face_bot2);
    assert(dev_map_ycell2face_top1);
    assert(dev_map_ycell2face_top2);

    //size_t local_work_size  = MIN(ncells, TILE_SIZE);
    //size_t local_work_size  = ncells;
    //size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
    //size_t global_work_size = MAX(MAX(nxface, nyface), ncells) + 1;
    size_t local_work_size = CL_DEVICE_MAX_WORK_GROUP_SIZE;
    size_t global_work_size = ncells;
/*
    vector<int>this_xfaceIdxList(ncells);
    vector<int>this_yfaceIdxList(ncells);
    vector<int>this_nlft(ncells);
    vector<int>this_nrht(ncells);
    vector<int>this_nbot(ncells);
    vector<int>this_ntop(ncells);
    vector<int>this_level(ncells);
    ezcl_enqueue_read_buffer(command_queue, dev_nlft,     CL_TRUE, 0, (ncells)*sizeof(cl_int), &this_nlft[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_nrht,     CL_TRUE, 0, (ncells)*sizeof(cl_int), &this_nrht[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_nbot,     CL_TRUE, 0, (ncells)*sizeof(cl_int), &this_nbot[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_ntop,     CL_TRUE, 0, (ncells)*sizeof(cl_int), &this_ntop[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_level,     CL_TRUE, 0, (ncells)*sizeof(cl_int), &this_level[0], NULL);

    //gpu_timers[STATE_TIMER_READ] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);


    this_xfaceIdxList[0] = 0;
    this_yfaceIdxList[0] = 0;
    for (int ccc = 0; ccc < (int) ncells - 1; ccc++) {
        int xfaceSt = this_xfaceIdxList[ccc], yfaceSt = this_yfaceIdxList[ccc];
        if (this_nrht[ccc] != ccc) {
            xfaceSt++;
            //nxfaceCnt++;
        }
        if (this_ntop[ccc] != ccc) {
            yfaceSt++;
            //nyfaceCnt++;
        }
        if ((this_level[this_nrht[ccc]] > this_level[ccc]) && (this_ntop[this_nrht[ccc]] != this_nrht[ccc]) && (this_nbot[this_nrht[ccc]] != this_nrht[ccc])) {
            xfaceSt++; 
            //nxfaceCnt++;
        }
        if ((this_level[this_ntop[ccc]] > this_level[ccc]) && (this_nrht[this_ntop[ccc]] != this_ntop[ccc]) && (this_nlft[this_ntop[ccc]] != this_ntop[ccc])) {
            yfaceSt++;
            //nyfaceCnt++;
        }
        this_xfaceIdxList[ccc+1] = xfaceSt;
        this_yfaceIdxList[ccc+1] = yfaceSt;
    }

    ezcl_enqueue_write_buffer(command_queue, dev_xfaceIdxList, CL_TRUE, 0, ncells*sizeof(cl_int), &this_xfaceIdxList[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_yfaceIdxList, CL_TRUE, 0, ncells*sizeof(cl_int), &this_yfaceIdxList[0], NULL);
*/
    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

/*
__kernel void face_idx_wbidirmap_cl(
            __global          int   *level,                     // 0
            __global          int   *nlft,                      // 1
            __global          int   *nrht,                      // 2
            __global          int   *nbot,                      // 3
            __global          int   *ntop,                      // 4
            __global          int   *xfaceIdxList,              // 5
            __global          int   *yfaceIdxList,              // 6
                              int   ncells,                     // 7
            __global          int   *nface) {                   // 8
*/

    cl_event face_idx_event;
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 0, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 1, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 2, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 3, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 4, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 5, sizeof(cl_mem), (void *)&dev_xfaceIdxList);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 6, sizeof(cl_mem), (void *)&dev_yfaceIdxList);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 7, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_face_idx_wbidirmap, 8, sizeof(cl_mem), (void *)&dev_nface);
    size_t single_work = 1;
    
    ezcl_enqueue_ndrange_kernel(command_queue, kernel_face_idx_wbidirmap, 1, NULL, &single_work, &single_work, &face_idx_event);

    ezcl_wait_for_events(1, &face_idx_event);
    ezcl_event_release(face_idx_event);

    //gpu_timers[MESH_TIMER_BIDIRPART1] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

/*
__kernel void calc_face_list_wbidirmap_pt1_cl(
                        const int   ncells,                     // 0
            __global          int   *nface,                     // 1
                        const int   levmx,                      // 2
            __global          int   *level,                     // 3
            __global          int   *nlft,                      // 4
            __global          int   *nrht,                      // 5
            __global          int   *nbot,                      // 6
            __global          int   *ntop,                      // 7
            __global          int   *map_xface2cell_lower,      // 8
            __global          int   *map_xface2cell_upper,      // 9
            __global          int   *map_xcell2face_left1,      // 10
            __global          int   *map_xcell2face_left2,      // 11
            __global          int   *map_xcell2face_right1,     // 12
            __global          int   *map_xcell2face_right2,     // 13
            __global          int   *xface_level,               // 14
            __global          int   *xface_i,                   // 15
            __global          int   *xface_j,                   // 16
            __global          int   *ixmin_level,               // 17
            __global          int   *ixmax_level,               // 18
            __global          int   *jxmin_level,               // 19
            __global          int   *jxmax_level,               // 20
            __global          int   *map_yface2cell_lower,      // 21
            __global          int   *map_yface2cell_upper,      // 22
            __global          int   *map_ycell2face_bot1,       // 23
            __global          int   *map_ycell2face_bot2,       // 24
            __global          int   *map_ycell2face_top1,       // 25
            __global          int   *map_ycell2face_top2,       // 26
            __global          int   *yface_level,               // 27
            __global          int   *yface_i,                   // 28
            __global          int   *yface_j,                   // 29
            __global          int   *iymin_level,               // 30
            __global          int   *iymax_level,               // 31
            __global          int   *jymin_level,               // 32
            __global          int   *jymax_level,               // 33
            __global          int   *xfaceIdxList,              // 34
            __global          int   *yfaceIdxList,              // 35
            __global    const int   *i,                         // 36
            __global    const int   *j) {                       // 37
            */


    cl_event calc_face_list_wbidirmap_event1, calc_face_list_wbidirmap_event2, wbidirmap_precount_event, calc_wbidirmap_phantom_neighbors_event, calc_wbidirmap_phantom_values_event;
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 0, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 1, sizeof(cl_mem), (void *)&dev_nface);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 2, sizeof(cl_int), (void *)&levmx);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 3, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 4, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 5, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 6, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 7, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 8, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 9, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 10, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 11, sizeof(cl_mem), (void *)&dev_map_xcell2face_left2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 12, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 13, sizeof(cl_mem), (void *)&dev_map_xcell2face_right2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 14, sizeof(cl_mem), (void *)&dev_xface_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 15, sizeof(cl_mem), (void *)&dev_xface_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 16, sizeof(cl_mem), (void *)&dev_xface_j);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 17, sizeof(cl_mem), (void *)&dev_ixmin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 18, sizeof(cl_mem), (void *)&dev_ixmax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 19, sizeof(cl_mem), (void *)&dev_jxmin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 20, sizeof(cl_mem), (void *)&dev_jxmax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 21, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 22, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 23, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 24, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 25, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 26, sizeof(cl_mem), (void *)&dev_map_ycell2face_top2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 27, sizeof(cl_mem), (void *)&dev_yface_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 28, sizeof(cl_mem), (void *)&dev_yface_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 29, sizeof(cl_mem), (void *)&dev_yface_j);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 30, sizeof(cl_mem), (void *)&dev_iymin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 31, sizeof(cl_mem), (void *)&dev_iymax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 32, sizeof(cl_mem), (void *)&dev_jymin_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 33, sizeof(cl_mem), (void *)&dev_jymax_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 34, sizeof(cl_mem), (void *)&dev_xfaceIdxList);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 35, sizeof(cl_mem), (void *)&dev_yfaceIdxList);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 36, sizeof(cl_mem), (void *)&dev_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt1, 37, sizeof(cl_mem), (void *)&dev_j);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_face_list_wbidirmap_pt1, 1, NULL, &global_work_size, 0, &calc_face_list_wbidirmap_event1);

    ezcl_wait_for_events(1, &calc_face_list_wbidirmap_event1);
    ezcl_event_release(calc_face_list_wbidirmap_event1);

    //gpu_timers[MESH_TIMER_BIDIRPART2] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    /*
__kernel void calc_face_list_wbidirmap_pt2_cl(
                        const int   ncells,                     // 0
            __global          int   *level,                     // 1
            __global          int   *nlft,                      // 2
            __global          int   *nrht,                      // 3
            __global          int   *nbot,                      // 4
            __global          int   *ntop,                      // 5
            __global          int   *map_xcell2face_left1,      // 6
            __global          int   *map_xcell2face_left2,      // 7
            __global          int   *map_xcell2face_right1,     // 8
            __global          int   *map_xcell2face_right2,     // 9
            __global          int   *map_ycell2face_bot1,       // 10
            __global          int   *map_ycell2face_bot2,       // 11
            __global          int   *map_ycell2face_top1,       // 12
            __global          int   *map_ycell2face_top2,       // 13
            __global    const int   *i,                         // 14
            __global    const int   *j) {                       // 15
    */

    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 0, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 1, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 2, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 3, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 4, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 5, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 6, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 7, sizeof(cl_mem), (void *)&dev_map_xcell2face_left2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 8, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 9, sizeof(cl_mem), (void *)&dev_map_xcell2face_right2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 10, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 11, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 12, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 13, sizeof(cl_mem), (void *)&dev_map_ycell2face_top2);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 14, sizeof(cl_mem), (void *)&dev_i);
    ezcl_set_kernel_arg(kernel_calc_face_list_wbidirmap_pt2, 15, sizeof(cl_mem), (void *)&dev_j);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_face_list_wbidirmap_pt2, 1, NULL, &global_work_size, 0, &calc_face_list_wbidirmap_event2);

    ezcl_wait_for_events(1, &calc_face_list_wbidirmap_event2);
    ezcl_event_release(calc_face_list_wbidirmap_event2);


    //gpu_timers[MESH_TIMER_BIDIRPART3] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    // read in device memory
    int nxface, nyface;
    ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nxface, NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 1*sizeof(cl_int), sizeof(cl_int), &nyface, NULL);
    //printf("%d %d\n", nxface, nyface);
    
    /*
__kernel void wbidirmap_precount_cl(
            __global          int   *pxcellCnt,                 // 0
            __global          int   *pycellCnt,                 // 1
            __global          int   *pxfaceCnt,                 // 2
            __global          int   *pyfaceCnt,                 // 3
            __global          int   *level,                     // 4
            __global          int   *nlft,                      // 5
            __global          int   *nrht,                      // 6
            __global          int   *nbot,                      // 7
            __global          int   *ntop,                      // 8
            __global          int   *map_xface2cell_lower,      // 9
            __global          int   *map_xface2cell_upper,      // 10
            __global          int   *map_yface2cell_lower,      // 11
            __global          int   *map_yface2cell_upper,      // 12
            __global          int   *nface,                     // 13 Number array of faces
            __global          int   *ifixupXCnt,              // 14
            __global          int   *ifixupYCnt) {            // 15
     */

    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 0, sizeof(cl_mem), (void *)&dev_pxcellCnt);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 1, sizeof(cl_mem), (void *)&dev_pycellCnt);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 2, sizeof(cl_mem), (void *)&dev_pxfaceCnt);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 3, sizeof(cl_mem), (void *)&dev_pyfaceCnt);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 4, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 5, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 6, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 7, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 8, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 9, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 10, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 11, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 12, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 13, sizeof(cl_mem), (void *)&dev_nface);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 14, sizeof(cl_mem), (void *)&dev_ifixupXCnt);
    ezcl_set_kernel_arg(kernel_wbidirmap_precount, 15, sizeof(cl_mem), (void *)&dev_ifixupYCnt);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_wbidirmap_precount, 1, NULL, &global_work_size, 0, &wbidirmap_precount_event);

    ezcl_wait_for_events(1, &wbidirmap_precount_event);
    ezcl_event_release(wbidirmap_precount_event);

    //gpu_timers[MESH_TIMER_BIDIRPART4] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    vector<int>cxCnt(nxface);
    vector<int>cyCnt(nyface);
    vector<int>xCnt(nxface);
    vector<int>yCnt(nyface);
    vector<int>fixupX(nxface);
    vector<int>fixupY(nyface);
    ezcl_enqueue_read_buffer(command_queue, dev_pxcellCnt,     CL_TRUE, 0, nxface*sizeof(cl_int), &cxCnt[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_pycellCnt,     CL_TRUE, 0, nyface*sizeof(cl_int), &cyCnt[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_pxfaceCnt,     CL_TRUE, 0, nxface*sizeof(cl_int), &xCnt[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_pyfaceCnt,     CL_TRUE, 0, nyface*sizeof(cl_int), &yCnt[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_ifixupXCnt,     CL_TRUE, 0, nxface*sizeof(cl_int), &fixupX[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_ifixupYCnt,     CL_TRUE, 0, nyface*sizeof(cl_int), &fixupY[0], NULL);
    
    //gpu_timers[STATE_TIMER_READ] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    //printf("\n\n"); 
    //for (int dumb = 0; dumb < nxface; dumb++) { printf("%d %d\t", cyCnt[dumb], yCnt[dumb]); }
    //printf("\n"); 
    //run reduce on the CPU version of device memory
    int xcellCnt = gpu_serial_int_reduce(&cxCnt[0], ncells, nxface);
    int ycellCnt = gpu_serial_int_reduce(&cyCnt[0], xcellCnt, nyface) - xcellCnt + ncells;
        pxfaceCnt = gpu_serial_int_reduce(&xCnt[0], nxface, nxface);
        pyfaceCnt = gpu_serial_int_reduce(&yCnt[0], nyface, nyface);
        pcellCnt = xcellCnt + ycellCnt - ncells; // subtract ncells b/c bot xcellCnt & ycellCnt include it
        nxfixup = gpu_serial_int_reduce(&fixupX[0], 0, nxface);
        nyfixup = gpu_serial_int_reduce(&fixupY[0], 0, nyface);
    //for (int dumb = 0; dumb < nyface; dumb++) { printf("%d) %d\n", dumb, cyCnt[dumb]); }
    //printf("\n\n"); 
    //printf("%d\n", xcellCnt - ncells);

    //printf("%d %d %d %d %d\n", pcellCnt, pxfaceCnt, pyfaceCnt, nxface, nyface);
/*    
    gpu_wbidirmap_realloc(&dev_map_xface2cell_lower, nxface, pxfaceCnt);
    gpu_wbidirmap_realloc(&dev_map_xface2cell_upper, nxface, pxfaceCnt);
    gpu_wbidirmap_realloc(&dev_map_yface2cell_lower, nyface, pyfaceCnt);
    gpu_wbidirmap_realloc(&dev_map_yface2cell_upper, nyface, pyfaceCnt);
    gpu_wbidirmap_realloc(&dev_map_xcell2face_left1, ncells, pcellCnt);
    gpu_wbidirmap_realloc(&dev_map_xcell2face_right1, ncells, pcellCnt);
    gpu_wbidirmap_realloc(&dev_map_ycell2face_bot1, ncells, pcellCnt);
    gpu_wbidirmap_realloc(&dev_map_ycell2face_top1, ncells, pcellCnt);
*/
    //gpu_timers[MESH_TIMER_BIDIRPART5] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    gpu_wbidirmap_realloc(&dev_level, ncells, pcellCnt);
    gpu_wbidirmap_realloc(&dev_i, ncells, pcellCnt);
    gpu_wbidirmap_realloc(&dev_j, ncells, pcellCnt);
    gpu_wbidirmap_realloc(&dev_xface_level, nxface, pxfaceCnt);
    gpu_wbidirmap_realloc(&dev_xface_i, nxface, pxfaceCnt);
    gpu_wbidirmap_realloc(&dev_xface_j, nxface, pxfaceCnt);
    gpu_wbidirmap_realloc(&dev_yface_level, nyface, pyfaceCnt);
    gpu_wbidirmap_realloc(&dev_yface_i, nyface, pyfaceCnt);
    gpu_wbidirmap_realloc(&dev_yface_j, nyface, pyfaceCnt);

    size_t size_nxface = nxface, size_nyface = nyface;
    cl_mem dev_pxfaceIdx = ezcl_malloc(NULL, const_cast<char *>("dev_pxfaceIdx"), &size_nxface, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    cl_mem dev_pyfaceIdx = ezcl_malloc(NULL, const_cast<char *>("dev_pyfaceIdx"), &size_nyface, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    cl_mem dev_ifixupXStart = ezcl_malloc(NULL, const_cast<char *>("dev_ifixupXStart"), &size_nxface, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    cl_mem dev_ifixupYStart = ezcl_malloc(NULL, const_cast<char *>("dev_ifixupYStart"), &size_nyface, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    ezcl_enqueue_write_buffer(command_queue, dev_pxcellIdx, CL_TRUE, 0, nxface*sizeof(cl_int), &cxCnt[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_pycellIdx, CL_TRUE, 0, nyface*sizeof(cl_int), &cyCnt[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_pxfaceIdx, CL_TRUE, 0, nxface*sizeof(cl_int), &xCnt[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_pyfaceIdx, CL_TRUE, 0, nyface*sizeof(cl_int), &yCnt[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_ifixupXStart, CL_TRUE, 0, nxface*sizeof(cl_int), &fixupX[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_ifixupYStart, CL_TRUE, 0, nyface*sizeof(cl_int), &fixupY[0], NULL);
    vector<int>xminusCell(ncells, -1);
    vector<int>xplusCell(ncells, -1);
    vector<int>yminusCell(ncells, -1);
    vector<int>yplusCell(ncells, -1);
    ezcl_enqueue_write_buffer(command_queue, dev_xminusCell2Idx, CL_TRUE, 0, ncells*sizeof(cl_int), &xminusCell[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_xplusCell2Idx, CL_TRUE, 0, ncells*sizeof(cl_int), &xplusCell[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_yminusCell2Idx, CL_TRUE, 0, ncells*sizeof(cl_int), &yminusCell[0], NULL);
    ezcl_enqueue_write_buffer(command_queue, dev_yplusCell2Idx, CL_TRUE, 0, ncells*sizeof(cl_int), &yplusCell[0], NULL);

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    //gpu_timers[MESH_TIMER_BIDIRPART6] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);


    MallocPlus gpu_state_memory_old = gpu_state_memory;
    malloc_plus_memory_entry *memory_item;

    for (memory_item = gpu_state_memory_old.memory_entry_by_name_begin();
        memory_item != gpu_state_memory_old.memory_entry_by_name_end();
        memory_item = gpu_state_memory_old.memory_entry_by_name_next() ) {
      //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);
        cl_mem dev_state_mem_ptr = (cl_mem)memory_item->mem_ptr;

        if (memory_item->mem_elsize != 8) continue;
/*
__kernel void calc_wbidirmap_phantom_values_cl(
            __global          int      *level,                     // 0
            __global          int      *nlft,                      // 1
            __global          int      *nrht,                      // 2
            __global          int      *nbot,                      // 3
            __global          int      *ntop,                      // 4
            __global          int      *map_xface2cell_lower,      // 5
            __global          int      *map_xface2cell_upper,      // 6
            __global          int      *map_yface2cell_lower,      // 7
            __global          int      *map_yface2cell_upper,      // 8
            __global          int      *pxcellIdx,                 // 9
            __global          int      *pycellIdx,                 // 10
            __global          int      *nface,                     // 11 Number array of faces
            __global   const  double   *state,                     // 12
            __global          double   *state_new) {               // 13
 */
        cl_mem dev_state_mem_ptr_new = (cl_mem)gpu_state_memory.memory_malloc(pcellCnt, sizeof(cl_double), const_cast<char *>("dev_state_mem_ptr_new"), DEVICE_REGULAR_MEMORY);

        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 0, sizeof(cl_mem), (void *)&dev_level);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 1, sizeof(cl_mem), (void *)&dev_nlft);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 2, sizeof(cl_mem), (void *)&dev_nrht);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 3, sizeof(cl_mem), (void *)&dev_nbot);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 4, sizeof(cl_mem), (void *)&dev_ntop);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 5, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 6, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 7, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 8, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 9, sizeof(cl_mem), (void *)&dev_pxcellIdx);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 10, sizeof(cl_mem), (void *)&dev_pycellIdx);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 11, sizeof(cl_mem), (void *)&dev_nface);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 12, sizeof(cl_mem), (void *)&dev_state_mem_ptr);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 13, sizeof(cl_mem), (void *)&dev_state_mem_ptr_new);

        ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_wbidirmap_phantom_values, 1, NULL, &global_work_size, 0, &calc_wbidirmap_phantom_values_event);

        ezcl_wait_for_events(1, &calc_wbidirmap_phantom_values_event);
        ezcl_event_release(calc_wbidirmap_phantom_values_event);

   //vector<double>H_loc(pcellCnt);
   //ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr_new,     CL_TRUE, 0, pcellCnt*sizeof(cl_double), &H_loc[0], NULL);
   //for (int jello = 0; jello < pcellCnt; jello++) { printf("%d) %f\n", jello, H_loc[jello]); }

        gpu_state_memory.memory_replace(dev_state_mem_ptr, dev_state_mem_ptr_new);

    }

    //gpu_timers[MESH_TIMER_BIDIRPART7] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);


    /*
__kernel void calc_wbidirmap_phantom_neighbors_cl(
                         const     int    ncells,                    // 0
                         const     real_t deltaT,                    // 1 Size of time step
            __global     const     int   *nface,                     // 2
            __global               int   *pxcellIdx,                 // 3
            __global               int   *pycellIdx,                 // 4
            __global               int   *pxfaceIdx,                 // 5
            __global               int   *pyfaceIdx,                 // 6
            __global               int   *ifixupXStart,              // 7
            __global               int   *ifixupYStart,              // 8
            __global               int   *level,                     // 9
            __global               int   *nlft,                      // 10
            __global               int   *nrht,                      // 11
            __global               int   *nbot,                      // 12
            __global               int   *ntop,                      // 13
            __global               int   *map_xface2cell_lower,      // 14
            __global               int   *map_xface2cell_upper,      // 15
            __global               int   *map_xcell2face_left1,      // 16
            __global               int   *map_xcell2face_right1,     // 17
            __global               int   *map_yface2cell_lower,      // 18
            __global               int   *map_yface2cell_upper,      // 19
            __global               int   *map_ycell2face_bot1,       // 20
            __global               int   *map_ycell2face_top1,       // 21
            __global               int   *xrecvIdx,                  // 22
            __global               int   *xrecvCIdx,                 // 23
            __global               int   *xplusCell2Idx,             // 24
            __global               int   *xminusCell2Idx,            // 25
            __global               int   *xsendIdx1,                 // 26
            __global               int   *xsendIdx2,                 // 27
            __global               int   *yrecvIdx,                  // 28
            __global               int   *yrecvCIdx,                 // 29
            __global               int   *yplusCell2Idx,             // 30
            __global               int   *yminusCell2Idx,            // 31
            __global               int   *ysendIdx1,                 // 32
            __global               int   *ysendIdx2,                 // 33
            __global               int   *i,                         // 34
            __global               int   *j,                         // 35
            __global               int   *xface_level,               // 36
            __global               int   *xface_i,                   // 37
            __global               int   *xface_j,                   // 38
            __global               int   *yface_level,               // 39
            __global               int   *yface_i,                   // 40
            __global               int   *yface_j)                   // 41
     */

    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 0, sizeof(cl_int), (void *)&ncells);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 1, sizeof(cl_real_t), (void *)&deltaT);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 2, sizeof(cl_mem), (void *)&dev_nface);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 3, sizeof(cl_mem), (void *)&dev_pxcellIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 4, sizeof(cl_mem), (void *)&dev_pycellIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 5, sizeof(cl_mem), (void *)&dev_pxfaceIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 6, sizeof(cl_mem), (void *)&dev_pyfaceIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 7, sizeof(cl_mem), (void *)&dev_ifixupXStart);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 8, sizeof(cl_mem), (void *)&dev_ifixupYStart);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 9, sizeof(cl_mem), (void *)&dev_level);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 10, sizeof(cl_mem), (void *)&dev_nlft);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 11, sizeof(cl_mem), (void *)&dev_nrht);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 12, sizeof(cl_mem), (void *)&dev_nbot);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 13, sizeof(cl_mem), (void *)&dev_ntop);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 14, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 15, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 16, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 17, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 18, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 19, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 20, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 21, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 22, sizeof(cl_mem), (void *)&dev_xrecvIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 23, sizeof(cl_mem), (void *)&dev_xrecvCIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 24, sizeof(cl_mem), (void *)&dev_xplusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 25, sizeof(cl_mem), (void *)&dev_xminusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 26, sizeof(cl_mem), (void *)&dev_xsendIdx1);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 27, sizeof(cl_mem), (void *)&dev_xsendIdx2);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 28, sizeof(cl_mem), (void *)&dev_yrecvIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 29, sizeof(cl_mem), (void *)&dev_yrecvCIdx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 30, sizeof(cl_mem), (void *)&dev_yplusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 31, sizeof(cl_mem), (void *)&dev_yminusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 32, sizeof(cl_mem), (void *)&dev_ysendIdx1);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 33, sizeof(cl_mem), (void *)&dev_ysendIdx2);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 34, sizeof(cl_mem), (void *)&dev_i);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 35, sizeof(cl_mem), (void *)&dev_j);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 36, sizeof(cl_mem), (void *)&dev_xface_level);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 37, sizeof(cl_mem), (void *)&dev_xface_i);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 38, sizeof(cl_mem), (void *)&dev_xface_j);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 39, sizeof(cl_mem), (void *)&dev_yface_level);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 40, sizeof(cl_mem), (void *)&dev_yface_i);
    ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_neighbors, 41, sizeof(cl_mem), (void *)&dev_yface_j);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_wbidirmap_phantom_neighbors, 1, NULL, &global_work_size, 0, &calc_wbidirmap_phantom_neighbors_event);

    ezcl_wait_for_events(1, &calc_wbidirmap_phantom_neighbors_event);
    ezcl_event_release(calc_wbidirmap_phantom_neighbors_event);

    //gpu_timers[MESH_TIMER_BIDIRPART8] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    //vector<int>map_l(nxfixup);
    //ezcl_enqueue_read_buffer(command_queue, dev_xrecvIdx,     CL_TRUE, 0, nxfixup*sizeof(cl_int), &map_l[0], NULL);
    //vector<int>map_r(nxfixup);
    //ezcl_enqueue_read_buffer(command_queue, dev_xsendIdx1,     CL_TRUE, 0, nxfixup*sizeof(cl_int), &map_r[0], NULL);
    //vector<int>map_t(nxfixup);
    //ezcl_enqueue_read_buffer(command_queue, dev_xsendIdx2,     CL_TRUE, 0, nxfixup*sizeof(cl_int), &map_t[0], NULL);
    //vector<int>map_r(pyfaceCnt);
    //ezcl_enqueue_read_buffer(command_queue, dev_map_yface2cell_lower,     CL_TRUE, 0, pyfaceCnt*sizeof(cl_int), &map_r[0], NULL);
    //for (int jello = 0; jello < pyfaceCnt; jello++) {  printf("%d) %d\n", jello, map_r[jello]); }
    //vector<int>recvCx(nxfixup);
    //vector<int>map_yface(pyfaceCnt);
    //ezcl_enqueue_read_buffer(command_queue, dev_map_yface2cell_upper,     CL_TRUE, 0, pyfaceCnt*sizeof(cl_int), &map_yface[0], NULL);
    //vector<int>recvCx(nxfixup);
    //ezcl_enqueue_read_buffer(command_queue, dev_xrecvCIdx,     CL_TRUE, 0, nxfixup*sizeof(cl_int), &recvCx[0], NULL);
    //vector<int>recvCy(nyfixup);
    //ezcl_enqueue_read_buffer(command_queue, dev_yrecvIdx,     CL_TRUE, 0, nyfixup*sizeof(cl_int), &recvCy[0], NULL);
    //printf("%d %d\n", xcellCnt, pcellCnt);


    ezcl_device_memory_delete(dev_pxfaceIdx);
    ezcl_device_memory_delete(dev_pyfaceIdx);
    ezcl_device_memory_delete(dev_ifixupXStart);
    ezcl_device_memory_delete(dev_ifixupYStart);

   
    gpu_timers[MESH_TIMER_BIDIR] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
}

void Mesh::gpu_calc_face_list_phantom_fill(MallocPlus &gpu_state_memory, double deltaT)
{

    struct timespec tstart_cpu;
    cpu_timer_start(&tstart_cpu);

    struct timespec tstart_cpu_part;
    cpu_timer_start(&tstart_cpu_part);
    
    cl_command_queue command_queue = ezcl_get_command_queue();
    cl_event calc_wbidirmap_phantom_values_event;

    assert(dev_nlft);
    assert(dev_nrht);
    assert(dev_nbot);
    assert(dev_ntop);
    assert(dev_level);
    assert(dev_levdx);
    assert(dev_levdy);
    assert(dev_map_xface2cell_lower);
    assert(dev_map_xface2cell_upper);
    assert(dev_map_xcell2face_left1);
    assert(dev_map_xcell2face_left2);
    assert(dev_map_xcell2face_right1);
    assert(dev_map_xcell2face_right2);
    assert(dev_map_yface2cell_lower);
    assert(dev_map_yface2cell_upper);
    assert(dev_map_ycell2face_bot1);
    assert(dev_map_ycell2face_bot2);
    assert(dev_map_ycell2face_top1);
    assert(dev_map_ycell2face_top2);

    //size_t local_work_size  = MIN(ncells, TILE_SIZE);
    //size_t local_work_size  = ncells;
    //size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
    //size_t global_work_size = MAX(MAX(nxface, nyface), ncells) + 1;
    size_t local_work_size = CL_DEVICE_MAX_WORK_GROUP_SIZE;
    size_t global_work_size = ncells;

    MallocPlus gpu_state_memory_old = gpu_state_memory;
    malloc_plus_memory_entry *memory_item;

    for (memory_item = gpu_state_memory_old.memory_entry_by_name_begin();
        memory_item != gpu_state_memory_old.memory_entry_by_name_end();
        memory_item = gpu_state_memory_old.memory_entry_by_name_next() ) {
      //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);
        cl_mem dev_state_mem_ptr = (cl_mem)memory_item->mem_ptr;

        if (memory_item->mem_elsize != 8) continue;
        

/*
__kernel void calc_wbidirmap_phantom_values_cl(
            __global          int      *level,                     // 0
            __global          int      *nlft,                      // 1
            __global          int      *nrht,                      // 2
            __global          int      *nbot,                      // 3
            __global          int      *ntop,                      // 4
            __global          int      *map_xface2cell_lower,      // 5
            __global          int      *map_xface2cell_upper,      // 6
            __global          int      *map_yface2cell_lower,      // 7
            __global          int      *map_yface2cell_upper,      // 8
            __global          int      *pxcellIdx,                 // 9
            __global          int      *pycellIdx,                 // 10
            __global          int      *nface,                     // 11 Number array of faces
            __global   const  double   *state,                     // 12
            __global          double   *state_new) {               // 13
 */
    
        cl_mem dev_state_mem_ptr_new = (cl_mem)gpu_state_memory.memory_malloc(pcellCnt, sizeof(cl_double), const_cast<char *>("dev_state_mem_ptr_new"), DEVICE_REGULAR_MEMORY);

        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 0, sizeof(cl_mem), (void *)&dev_level);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 1, sizeof(cl_mem), (void *)&dev_nlft);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 2, sizeof(cl_mem), (void *)&dev_nrht);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 3, sizeof(cl_mem), (void *)&dev_nbot);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 4, sizeof(cl_mem), (void *)&dev_ntop);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 5, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 6, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 7, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 8, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 9, sizeof(cl_mem), (void *)&dev_pxcellIdx);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 10, sizeof(cl_mem), (void *)&dev_pycellIdx);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 11, sizeof(cl_mem), (void *)&dev_nface);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 12, sizeof(cl_mem), (void *)&dev_state_mem_ptr);
        ezcl_set_kernel_arg(kernel_calc_wbidirmap_phantom_values, 13, sizeof(cl_mem), (void *)&dev_state_mem_ptr_new);

        ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_wbidirmap_phantom_values, 1, NULL, &global_work_size, 0, &calc_wbidirmap_phantom_values_event);

        ezcl_wait_for_events(1, &calc_wbidirmap_phantom_values_event);
        ezcl_event_release(calc_wbidirmap_phantom_values_event);

   //vector<double>H_loc(pcellCnt);
   //ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr_new,     CL_TRUE, 0, pcellCnt*sizeof(cl_double), &H_loc[0], NULL);
   //for (int jello = 0; jello < pcellCnt; jello++) { printf("%d) %f\n", jello, H_loc[jello]); }

        gpu_state_memory.memory_replace(dev_state_mem_ptr, dev_state_mem_ptr_new);

    }

    gpu_timers[MESH_TIMER_BIDIRPART7] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    cpu_timer_start(&tstart_cpu_part);



}
#endif

void Mesh::generate_regular_cell_meshes(MallocPlus &state_memory)
{
   //printf("DEBUG -- imin, imax %d %d\n",imin,imax);
   int nvar=3;
/*
   malloc_plus_memory_entry *memory_item;
   MallocPlus state_memory_old = state_memory;
   for (memory_item = state_memory_old.memory_entry_by_name_begin();
       memory_item != state_memory_old.memory_entry_by_name_end();
       memory_item = state_memory_old.memory_entry_by_name_next() ) {
          
       //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

       if ( (memory_item->mem_flags & REZONE_DATA) == 0) continue;
       nvar++;
   }
   //nvar--;
   //printf("DEBUG -- levmx %d nvar %d jmax %d jmin %d imax %d imin %d\n",levmx,nvar,jmax,jmin,imax,imin);
*/

#ifdef HAVE_OPENCL
   cl_command_queue command_queue = ezcl_get_command_queue();

    //xface_i = (int *)mesh_memory.memory_realloc(pxfaceCnt, xface_i);
    //xface_j = (int *)mesh_memory.memory_realloc(pxfaceCnt, xface_j);
    //xface_level = (uchar_t *)mesh_memory.memory_realloc(pxfaceCnt, xface_level);
    //yface_i = (int *)mesh_memory.memory_realloc(pyfaceCnt, yface_i);
    //yface_j = (int *)mesh_memory.memory_realloc(pyfaceCnt, yface_j);
    //yface_level = (uchar_t *)mesh_memory.memory_realloc(pyfaceCnt, yface_level);
    //i = (int *)mesh_memory.memory_realloc(pcellCnt, i);
    //j = (int *)mesh_memory.memory_realloc(pcellCnt, j);
    //level = (uchar_t *)mesh_memory.memory_realloc(pcellCnt, level);



    i = (int *)mesh_memory.memory_realloc(pcellCnt, i);
    j = (int *)mesh_memory.memory_realloc(pcellCnt, j);
    level = (uchar_t *)mesh_memory.memory_realloc(pcellCnt, level);
    xface_i = (int *)mesh_memory.memory_realloc(pxfaceCnt, xface_i);
    xface_j = (int *)mesh_memory.memory_realloc(pxfaceCnt, xface_j);
    xface_level = (uchar_t *)mesh_memory.memory_realloc(pxfaceCnt, xface_level);
    yface_i = (int *)mesh_memory.memory_realloc(pyfaceCnt, yface_i);
    yface_j = (int *)mesh_memory.memory_realloc(pyfaceCnt, yface_j);
    yface_level = (uchar_t *)mesh_memory.memory_realloc(pyfaceCnt, yface_level);
    memory_reset_ptrs();

    ezcl_enqueue_read_buffer(command_queue, dev_xface_i, CL_TRUE, 0, pxfaceCnt*sizeof(cl_int), &xface_i[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_xface_j, CL_TRUE, 0, pxfaceCnt*sizeof(cl_int), &xface_j[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_xface_level, CL_TRUE, 0, pxfaceCnt*sizeof(cl_int), &xface_level[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_yface_i, CL_TRUE, 0, pyfaceCnt*sizeof(cl_int), &yface_i[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_yface_j, CL_TRUE, 0, pyfaceCnt*sizeof(cl_int), &yface_j[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_yface_level, CL_TRUE, 0, pyfaceCnt*sizeof(cl_int), &yface_level[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_i, CL_TRUE, 0, pcellCnt*sizeof(cl_int), &i[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_j, CL_TRUE, 0, pcellCnt*sizeof(cl_int), &j[0], NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE, 0, pcellCnt*sizeof(cl_int), &level[0], NULL);
#endif
   
   meshes = (mesh_type*)malloc((levmx+1)*sizeof(mesh_type));
   int ***avgCnt = (int ***)malloc((levmx+1) * sizeof(int **));

   int ll, pjIdx, piIdx;
   for (ll=0; ll<levmx+1; ll++){
      lev_iregmin[ll] = lev_ibegin[ll]-2;
      lev_jregmin[ll] = lev_jbegin[ll]-2;
      lev_iregsize[ll] = lev_iend[ll]-lev_ibegin[ll]+5;
      lev_jregsize[ll] = lev_jend[ll]-lev_jbegin[ll]+5;
      //printf("DEBUG -- lev_iregmin %d lev_iregsize %d lev_jregmin %d lev_jregsize %d\n",
      //        lev_iregmin[ll],lev_iregsize[ll],lev_jregmin[ll],lev_jregsize[ll]);

      int col[lev_iregsize[ll]];
      for(int ii=0; ii<lev_iregsize[ll]; ii++){
         col[ii] = 0;
      }
      int row[lev_jregsize[ll]];
      for(int jj=0; jj<lev_jregsize[ll]; jj++){
         row[jj] = 0;
      }
      for (uint ic=0; ic < ncells; ic++){
         if ((int)level[ic] == ll) {
            col[i[ic]-lev_iregmin[ll]] = 1;
            row[j[ic]-lev_jregmin[ll]] = 1;
         }
      }

      /*for (int ii = 2; ii < lev_iregsize[ll]-2; ii++) {
        printf("%d\t", col[ii]);
      }
      printf("\n");
      for (int jj = 2; jj < lev_jregsize[ll]-2; jj++) {
        printf("%d\t", row[jj]);
      }
      printf("\n");*/
        
      int trimleft=0;
      for(int ii=2; ii<lev_iregsize[ll]-2; ii++){
         if (col[ii] == 1) break;
         trimleft++;
      }
      int trimright=0;
      for(int ii=lev_iregsize[ll]-3; ii >= 2; ii--){
         if (col[ii] == 1) break;
         trimright++;
      }
      int trimbottom=0;
      for(int jj=2; jj<lev_jregsize[ll]-2; jj++){
         if (row[jj] == 1) break;
         trimbottom++;
      }
      int trimtop=0;
      for(int jj=lev_iregsize[ll]-3; jj >= 2; jj--){
         if (row[jj] == 1) break;
         trimtop++;
      }
      //printf("DEBUG -- trimleft %d trimright %d trimbottom %d trimtop %d\n",trimleft,trimright,trimbottom,trimtop);

      lev_iregmin[ll] += trimleft;
      lev_jregmin[ll] += trimbottom;
      lev_iregsize[ll] -= (trimleft+trimright);
      lev_jregsize[ll] -= (trimbottom+trimtop);

      if (lev_iregsize[ll] <= 0) lev_iregsize[ll] = 1;
      if (lev_jregsize[ll] <= 0) lev_jregsize[ll] = 1;
        
      //printf("DEBUG -- lev_iregmin %d lev_iregsize %d lev_jregmin %d lev_jregsize %d\n",
      //        lev_iregmin[ll],lev_iregsize[ll],lev_jregmin[ll],lev_jregsize[ll]);

      //printf("%d %d %d\n", nvar, lev_jregsize[ll], lev_iregsize[ll]);
#ifdef FULL_PRECISION
      meshes[ll].pstate = (double ***)gentrimatrix(nvar,lev_jregsize[ll],lev_iregsize[ll],sizeof(double));
#elif HALF_PRECISION
      meshes[ll].pstate = (half ***)gentrimatrix(nvar,lev_jregsize[ll],lev_iregsize[ll],sizeof(half));
#else
      meshes[ll].pstate = (float ***)gentrimatrix(nvar,lev_jregsize[ll],lev_iregsize[ll],sizeof(float));
#endif
      meshes[ll].mask = (int **)genmatrix(lev_jregsize[ll],lev_iregsize[ll],sizeof(int));
      avgCnt[ll] = (int **)genmatrix(lev_jregsize[ll],lev_iregsize[ll],sizeof(int));

      for(int nn=0; nn<nvar; nn++){
         for(int jj=0; jj<lev_jregsize[ll]; jj++){
            for(int ii=0; ii<lev_iregsize[ll]; ii++){
               meshes[ll].pstate[nn][jj][ii]=0.0;
            }
         }
      }
      for(int jj=0; jj<lev_jregsize[ll]; jj++){
         for(int ii=0; ii<lev_iregsize[ll]; ii++){
            meshes[ll].mask[jj][ii]=0;
            avgCnt[ll][jj][ii] = 0;
         }
      }
   }

   for (int ic = 0; ic < (int)ncells; ic++) {
       ll = level[ic];
       pjIdx = j[ic] - lev_jregmin[ll];
       piIdx = i[ic] - lev_iregmin[ll];  
       if (nlft[ic] == ic || nrht[ic] == ic || nbot[ic] == ic || ntop[ic] == ic) continue;
       meshes[ll].mask[pjIdx][piIdx] = 1;
   }
   for (int iface = 0; iface < nxface; iface++) {
       ll = xface_level[iface];
       pjIdx = xface_j[iface] - lev_jregmin[ll];
       piIdx = xface_i[iface] - lev_iregmin[ll];
   }

   for (int iface = 0; iface < nyface; iface++) {
       ll = yface_level[iface];
       pjIdx = yface_j[iface] - lev_jregmin[ll];
       piIdx = yface_i[iface] - lev_iregmin[ll];
   }

   int ivar = 0;

#ifdef HAVE_OPENCL
    ncells_phan = pcellCnt;
    MallocPlus gpu_state_memory_old = state_memory;
    malloc_plus_memory_entry *memory_item;
    double *mem_ptr_double = (double *)malloc(pcellCnt*sizeof(double));

    for (memory_item = gpu_state_memory_old.memory_entry_by_name_begin();
        memory_item != gpu_state_memory_old.memory_entry_by_name_end();
        memory_item = gpu_state_memory_old.memory_entry_by_name_next() ) { //}
      //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);
        cl_mem dev_state_mem_ptr = (cl_mem)memory_item->mem_ptr;

        ezcl_enqueue_read_buffer(command_queue, dev_state_mem_ptr,     CL_TRUE, 0, pcellCnt*sizeof(cl_real_t), &mem_ptr_double[0], NULL);
        
#else
   MallocPlus state_memory_old = state_memory;
   malloc_plus_memory_entry *memory_item;
   for (memory_item = state_memory_old.memory_entry_by_name_begin();
       memory_item != state_memory_old.memory_entry_by_name_end();
       memory_item = state_memory_old.memory_entry_by_name_next() ) {
       
       //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

       if ( (memory_item->mem_flags & REZONE_DATA) == 0) continue;

       double *mem_ptr_double = (double *)memory_item->mem_ptr;
#endif
       //add original cell values into regular grid
       for (uint ic=0; ic < ncells; ic++){
          ll = level[ic];
          pjIdx = j[ic] - lev_jregmin[ll];
          piIdx = i[ic] - lev_iregmin[ll];  
          meshes[ll].pstate[ivar][pjIdx][piIdx] += mem_ptr_double[ic];
          avgCnt[ll][pjIdx][piIdx]++;
       }

       //add phantom cell values into regular grid (if it doesn't conflict with an original cell)
       for (uint ic = ncells; ic < ncells_phan; ic++) {
           ll = level[ic];
           pjIdx = j[ic] - lev_jregmin[ll];
           piIdx = i[ic] - lev_iregmin[ll];  
           if (meshes[ll].mask[pjIdx][piIdx] == 1) continue;
           meshes[ll].pstate[ivar][pjIdx][piIdx] += mem_ptr_double[ic];
           avgCnt[ll][pjIdx][piIdx]++;
       }

       //interpolation of phantom cell values in regular grid (if it doesn't conflict with an original cell)
       for (uint ic = ncells; ic < ncells_phan; ic++) {
           ll = level[ic]; 
           pjIdx = j[ic] - lev_jregmin[ll];
           piIdx = i[ic] - lev_iregmin[ll];  
           if ((meshes[ll].mask[pjIdx][piIdx] == 1) || (avgCnt[ll][pjIdx][piIdx] == 0)) continue;
           meshes[ll].pstate[ivar][pjIdx][piIdx] /= avgCnt[ll][pjIdx][piIdx];
           avgCnt[ll][pjIdx][piIdx] = 0;
       }
       ivar++;

       /*if (ivar == 1) {
          for (int iface = 0; iface < nxface; iface++){
             printf("1st pass x direction iface %d i %d j %d lev %d nzlower %d nzupper %d %lf %lf\n",
                    iface, xface_i[iface], xface_j[iface], xface_level[iface],
                    map_xface2cell_lower[iface], map_xface2cell_upper[iface],
                    mem_ptr_double[map_xface2cell_lower[iface]], mem_ptr_double[map_xface2cell_upper[iface]]
                    );
          }
       }*/
   }
#ifdef HAVE_OPENCL
   free(mem_ptr_double);
#endif

   for (int lev = 0; lev < levmx+1; lev++) {
       genmatrixfree((void **)avgCnt[lev]);
   }
   free(avgCnt);
           
   /*for (int ll=0; ll<=levmx; ll++){
       printf("DEBUG regular mesh level %d\n",ll);
       for(int jj=0; jj<lev_jregsize[ll]; jj++){
           for(int ii=0; ii<lev_iregsize[ll]; ii++){
               //printf("   %lf %lf %lf    ",meshes[ll].pstate[0][jj][ii],meshes[ll].pstate[1][jj][ii],meshes[ll].pstate[2][jj][ii]);
               printf("  %d  ",meshes[ll].mask[jj][ii]);
           }
           printf("\n");
       }
       printf("\n");
       for(int jj=0; jj<lev_jregsize[ll]; jj++){
           for(int ii=0; ii<lev_iregsize[ll]; ii++){
               //printf("   %lf %lf %lf    ",meshes[ll].pstate[0][jj][ii],meshes[ll].pstate[1][jj][ii],meshes[ll].pstate[2][jj][ii]);
               if (meshes[ll].pstate[0][jj][ii] != -1.0){
                  printf("  %lf  ",meshes[ll].pstate[0][jj][ii]);
               } else {
                  printf("             ");
               }
           }
           printf("\n");
       }
       printf("\n");
   }

   printf("DEBUG -- nvar %d jmax %d jmin %d imax %d imin %d\n",nvar,jmax,jmin,imax,imin);
   */ 
}

void Mesh::destroy_regular_cell_meshes(MallocPlus &state_memory)
{
   /*for (int ll=0; ll<=levmx; ll++){
       printf("DEBUG just before destroy in mesh regular mesh level %d\n",ll);
       for(int jj=0; jj<lev_jregsize[ll]; jj++){
           for(int ii=0; ii<lev_iregsize[ll]; ii++){
               //printf("   %lf %lf %lf    ",meshes[ll].pstate[0][jj][ii],meshes[ll].pstate[1][jj][ii],meshes[ll].pstate[2][jj][ii]);
               printf("  %d  ",meshes[ll].mask[jj][ii]);
           }
           printf("\n");
       }
       printf("\n");
       for(int jj=0; jj<lev_jregsize[ll]; jj++){
           for(int ii=0; ii<lev_iregsize[ll]; ii++){
               //printf("   %lf %lf %lf    ",meshes[ll].pstate[0][jj][ii],meshes[ll].pstate[1][jj][ii],meshes[ll].pstate[2][jj][ii]);
               if (meshes[ll].pstate[0][jj][ii] != -1.0){
                  printf("  %lf  ",meshes[ll].pstate[0][jj][ii]);
               } else {
                  printf("             ");
               }
           }
           printf("\n");
       }
       printf("\n");
   }*/

#ifndef HAVE_OPENCL
   int ivar = 0;
   malloc_plus_memory_entry *memory_item;
   MallocPlus state_memory_old = state_memory;
   for (memory_item = state_memory_old.memory_entry_by_name_begin();
       memory_item != state_memory_old.memory_entry_by_name_end();
       memory_item = state_memory_old.memory_entry_by_name_next() ) {
       
       //printf("DEBUG -- it.mem_name %s elsize %lu\n",memory_item->mem_name,memory_item->mem_elsize);

       if ( (memory_item->mem_flags & REZONE_DATA) == 0) continue;

       double *mem_ptr_double = (double *)memory_item->mem_ptr;

       for (uint ic=0; ic < ncells; ic++){
          uchar_t ll = level[ic];
          //printf("DEBUG -- ic %d ll %d lev_iregmin[ll] %d lev_jregmin[ll] %d pstate %lf\n",
            //      ic, ll, lev_iregmin[ll], lev_jregmin[ll], meshes[ll].pstate[ivar][j[ic]-lev_jregmin[ll]][i[ic]-lev_iregmin[ll]]);
          mem_ptr_double[ic]=meshes[ll].pstate[ivar][j[ic]-lev_jregmin[ll]][i[ic]-lev_iregmin[ll]];
       }

       ivar++;
   }
#endif

   for (int ll=0; ll<levmx+1; ll++){
       gentrimatrixfree((void ***)meshes[ll].pstate);
       genmatrixfree((void **)meshes[ll].mask);
   }
   free(meshes);
}

int **Mesh::get_xface_flag(int lev, bool print_output)
{
    /*
   int **xface_flag = (int **)genmatrix(jxmax_level[lev]+1,
                                        ixmax_level[lev]+1, sizeof(int));
   for (int jj=0; jj<jxmax_level[lev]+1; jj++){
      for (int ii=0; ii<ixmax_level[lev]+1; ii++){
         xface_flag[jj][ii] = -1;
      }
   }

   for (int iface=0; iface < nxface; iface++){
      if (xface_level[iface] == lev){
         int ii = xface_i[iface];
         int jj = xface_j[iface];

         xface_flag[jj][ii] = 1;
      }
   }

   if (DEBUG || print_output) {
      printf("DEBUG -- x face_flag for level %d\n",lev);
      printf("DEBUG -- sizes isize+1 %d jsize+1 %d\n",ixmax_level[lev]+1,jxmax_level[lev]+1);

      printf("                           ");
      for (int ii=0; ii<ixmax_level[lev]+1; ii++){
         printf(" %4d ",ii);
      }
      printf("\n");

      for (int jj=jxmax_level[lev]; jj>=0; jj--){

         printf("DEBUG -- j  %4d:          ",jj);
         for (int ii=0; ii<ixmax_level[lev]+1; ii++){
            if (xface_flag[jj][ii] >= 0){
               //printf("      xface_flag_check[%d][%d] = 1;\n",jj,ii);
               printf(" %4d ", xface_flag[jj][ii]);
            } else {
               printf("      ");
            }
         }
         printf("\n");
      }
   }

   return(xface_flag);
   */
    return 0;
}

int **Mesh::get_yface_flag(int lev, bool print_output)
{
    /*
   int **yface_flag = (int **)genmatrix(jymax_level[lev]+1,
                                        iymax_level[lev]+1, sizeof(int));
   for (int jj=0; jj<jymax_level[lev]+1; jj++){
      for (int ii=0; ii<iymax_level[lev]+1; ii++){
         yface_flag[jj][ii] = -1;
      }
   }

   for (int iface=0; iface < nyface; iface++){
      if (yface_level[iface] == lev){
         int ii = yface_i[iface];
         int jj = yface_j[iface];

         yface_flag[jj][ii] = 1;
      }
   }

   if (DEBUG || print_output) {
      printf("DEBUG -- y face_flag for level %d\n",lev);
      printf("DEBUG -- sizes isize+1 %d jsize+1 %d\n",iymax_level[lev]+1,jymax_level[lev]+1);

      printf("                           ");
      for (int ii=0; ii<iymax_level[lev]+1; ii++){
         printf(" %4d ",ii);
      }
      printf("\n");

      for (int jj=jymax_level[lev]; jj>=0; jj--){

         printf("DEBUG -- j  %4d:          ",jj);
         for (int ii=0; ii<iymax_level[lev]+1; ii++){
            if (yface_flag[jj][ii] >= 0){
               //printf("      yface_flag_check[%d][%d] = 1;\n",jj,ii);
               printf(" %4d ", yface_flag[jj][ii]);
            } else {
               printf("      ");
            }
         }
         printf("\n");
      }
   }

   return(yface_flag);
   */
    return 0;
}

void Mesh::get_flat_grid(int lev, int ***zone_flag_base, int ***zone_cell_base)
{
    /*
   int isize = ixmax_level[lev]+4;
   int jsize = jymax_level[lev]+4;
   int iadjust = ixadjust[lev]-2;
   int jadjust = jyadjust[lev]-2;

   //printf("DEBUG -- sizes isize %d jsize %d\n",isize,jsize);
   //printf("DEBUG -- adjust ixadjust %d jxadjust %d\n",ixadjust[lev],jxadjust[lev]);
   //printf("DEBUG -- adjust iyadjust %d jyadjust %d\n",iyadjust[lev],jyadjust[lev]);

   (*zone_flag_base) = (int **)genmatrix(jsize, isize, sizeof(int));

   int **zone_flag = *zone_flag_base;
   for (int jj=0; jj<jsize; jj++){
      for (int ii=0; ii<isize; ii++){
          zone_flag[jj][ii] = -1;
      }
   }

   (*zone_cell_base) = (int **)genmatrix(jsize, isize, sizeof(int));

   int **zone_cell = *zone_cell_base;
   for (int jj=0; jj<jsize; jj++){
      for (int ii=0; ii<isize; ii++){
         zone_cell[jj][ii] = -1;
      }
   }

   for (int iface=0; iface < nxface; iface++){
      if (xface_level[iface] == lev){
         int nz1 = map_xface2cell_lower[iface];
         int nz2 = map_xface2cell_upper[iface];

         if (lev == level[nz1]) {
            int iii = i[nz1]-iadjust;
            int jjj = j[nz1]-jadjust;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz1;
            if (nlft[nz1] != REAL_CELL) {
               zone_cell[jjj][iii-1] = nlft[nz1];
            }
         } else {
            int iii = i[nz1]*2-iadjust+1;
            int jjj = j[nz1]*2-jadjust;
            if (is_upper(j[nz2])) jjj += 1;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz1;
            zone_cell[jjj][iii-1] = nz1;
         }
         if (lev == level[nz2]) {
            int iii = i[nz2]-iadjust;
            int jjj = j[nz2]-jadjust;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz2;
            if (nrht[nz2] != REAL_CELL) {
               zone_cell[jjj][iii+1] = nrht[nz2];
            }
         } else {
            int iii = i[nz2]*2-iadjust;
            int jjj = j[nz2]*2-jadjust;
            if (is_upper(j[nz1])) jjj += 1;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz2;
            zone_cell[jjj][iii+1] = nz2;
         }
      }
   }

   for (int iface=0; iface < nyface; iface++){
      if (yface_level[iface] == lev){
         int nz1 = map_yface2cell_lower[iface];
         int nz2 = map_yface2cell_upper[iface];

         if (lev == level[nz1]) {
            int iii = i[nz1]-iadjust;
            int jjj = j[nz1]-jadjust;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz1;
            if (nbot[nz1] != REAL_CELL) {
               zone_cell[jjj-1][iii] = nbot[nz1];
            }
         } else {
            int iii = i[nz1]*2-iadjust;
            int jjj = j[nz1]*2-jadjust+1;
            if (is_upper(i[nz2])) iii += 1;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz1;
            zone_cell[jjj-1][iii] = nz1;
         }
         if (lev == level[nz2]) {
            int iii = i[nz2]-iadjust;
            int jjj = j[nz2]-jadjust;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz2;
            if (ntop[nz2] != REAL_CELL) {
               zone_cell[jjj+1][iii] = ntop[nz2];
            }
         } else {
            int iii = i[nz2]*2-iadjust;
            int jjj = j[nz2]*2-jadjust;
            if (is_upper(i[nz1])) iii += 1;
            zone_flag[jjj][iii] = 1;
            zone_cell[jjj][iii] = nz2;
            zone_cell[jjj+1][iii] = nz2;
         }
      }
   }

   if (DEBUG) {
      printf("DEBUG -- zone_flag for level %d\n",lev);
      printf("DEBUG -- sizes isize %d jsize %d\n",isize,jsize);
      for (int j=jsize-1; j>=0; j--){
         for (int i=0; i<isize; i++){
            if (zone_flag[j][i] >= 0){
               printf("      zone_flag_check[%d][%d] = 1;\n",j,i);
            }
         }
      }
      for (int j=jsize-1; j>=0; j--){
         for (int i=0; i<isize; i++){
            if (zone_cell[j][i] >= 0){
               printf("      zone_cell_check[%d][%d] = %d;\n",j,i,zone_cell[j][i]);
            }
         }
      }

      printf("                  ");
      for (int i=0; i<isize; i++){
         printf(" %4d ",i);
      }
      printf("\n");

      for (int j=jsize-1; j>=0; j--){

         printf("DEBUG -- j  %4d: ",j);
         for (int i=0; i<isize; i++){
            if (zone_flag[j][i] >= 0){
               printf(" %4d ", zone_flag[j][i]);
            } else {
               printf("      ");
            }
         }
         printf("\n");
      }

      printf("DEBUG -- zone_cell for level %d\n",lev);

      printf("                  ");
      for (int i=0; i<isize; i++){
         printf(" %4d ",i);
      }
      printf("\n");

      for (int j=jsize-1; j>=0; j--){

         printf("DEBUG -- j  %4d: ",j);
         for (int i=0; i<isize; i++){
            if (zone_cell[j][i] >= 0){
               printf(" %4d ", zone_cell[j][i]);
            } else {
               printf("      ");
            }
         }
         printf("\n");
      }
   }
   */
}

void Mesh::calc_face_list_clearmaps()
{
    /*
   map_xface2cell_lower.clear();
   map_xface2cell_upper.clear();

   map_xcell2face_left1.clear();
   map_xcell2face_left2.clear();
   map_xcell2face_right1.clear();
   map_xcell2face_right2.clear();

   map_yface2cell_lower.clear();
   map_yface2cell_upper.clear();

   map_ycell2face_bot1.clear();
   map_ycell2face_bot2.clear();
   map_ycell2face_top1.clear();
   map_ycell2face_top2.clear();
   */
}

void Mesh::timer_output(mesh_timer_category category, mesh_device_types device_type, int timer_level)
{
   double local_time = 0.0;
   if (device_type == MESH_DEVICE_CPU){
      local_time = get_cpu_timer(category);
   } else {
      local_time = get_gpu_timer(category);
   }

   char string[80] = "/0";

   if (mype == 0) {
      const char *blank="          ";

      if (device_type == MESH_DEVICE_CPU){
         sprintf(string,"CPU: %.*s%-30.30s\t", 2*timer_level, blank, mesh_timer_descriptor[category]);
      } else {
         sprintf(string,"GPU: %.*s%-30.30s\t", 2*timer_level, blank, mesh_timer_descriptor[category]);
      }
   }

   parallel_output(string, local_time, timer_level, "s");
}

void Mesh::parallel_output(const char *string, double local_value, int output_level, const char *units)
{
   vector<double> global_values(numpe);
   global_values[0] = local_value;
#ifdef HAVE_MPI
   if (numpe > 1) { 
      MPI_Gather(&local_value, 1, MPI_DOUBLE, &global_values[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
#endif
   if (mype == 0) {
      const char *blank="          ";

      printf("%s\t",string);
      if (numpe <= 4) {
         for(int ip = 0; ip < numpe; ip++){
            printf("%.*s%8.4f\t", 2*output_level, blank, global_values[ip]);
         }
         printf("%s\n",units);
      } else {
         sort(global_values.begin(),global_values.end());
         double median_value;
         int half_value = numpe/2;
         if (numpe%2 == 0) {
            median_value = (global_values[half_value-1]+global_values[half_value])/2.0;
         } else {
            median_value = global_values[half_value+1];
         }
         printf("%.*s%8.4f\t%.*s%8.4f\t%.*s%8.4f   %s min/median/max\n",
            2*output_level, blank, global_values[0],
            2*output_level, blank, median_value,
            2*output_level, blank, global_values[numpe-1],
            units);
      }
   }
}

void Mesh::parallel_output(const char *string, long long local_value, int output_level, const char *units)
{
   vector<long long> global_values(numpe);
   global_values[0] = local_value;
#ifdef HAVE_MPI
   if (numpe > 1) { 
      MPI_Gather(&local_value, 1, MPI_LONG_LONG, &global_values[0], 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
   }
#endif
   if (mype == 0) {
      const char *blank="          ";

      printf("%s\t",string);
      if (numpe <= 4) {
         for(int ip = 0; ip < numpe; ip++){
            printf("%.*s%10lld\t", 2*output_level, blank, global_values[ip]);
         }
         printf("%s\n",units);
      } else {
         sort(global_values.begin(),global_values.end());
         long long median_value;
         int half_value = numpe/2;
         if (numpe%2 == 0) {
            median_value = (global_values[half_value-1]+global_values[half_value])/2;
         } else {
            median_value = global_values[half_value+1];
         }
         printf("%.*s%10lld\t%.*s%10lld\t%.*s%10lld   %s min/median/max\n",
            2*output_level, blank, global_values[0],
            2*output_level, blank, median_value,
            2*output_level, blank, global_values[numpe-1],
            units);
      }
   }
}

void Mesh::parallel_output(const char *string, int local_value, int output_level, const char *units)
{
   vector<int> global_values(numpe);
   global_values[0] = local_value;
#ifdef HAVE_MPI
   if (numpe > 1) { 
      MPI_Gather(&local_value, 1, MPI_INT, &global_values[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
   }
#endif
   if (mype == 0) {
      const char *blank="          ";

      printf("%s\t",string);
      if (numpe <= 4) {
         for(int ip = 0; ip < numpe; ip++){
            printf("%.*s%10d\t", 2*output_level, blank, global_values[ip]);
         }
         printf("%s\n",units);
      } else {
         sort(global_values.begin(),global_values.end());
         int median_value;
         int half_value = numpe/2;
         if (numpe%2 == 0) {
            median_value = (global_values[half_value-1]+global_values[half_value])/2;
         } else {
            median_value = global_values[half_value+1];
         }
         printf("%.*s%10d\t%.*s%10d\t%.*s%10d   %s min/median/max\n",
            2*output_level, blank, global_values[0],
            2*output_level, blank, median_value,
            2*output_level, blank, global_values[numpe-1],
            units);
      }
   }
}

const int CRUX_MESH_VERSION = 103;
const int num_int_dist_vals = 3;
const int num_int_vals      = 3;
const int num_double_vals   = 1;

size_t Mesh::get_checkpoint_size(void)
{
   size_t nsize;
   nsize  = num_int_dist_vals*sizeof(int);
   nsize += num_int_vals*sizeof(int);
   nsize += num_double_vals*sizeof(double);
   nsize += 2*MESH_COUNTER_SIZE*sizeof(int);
   nsize += MESH_TIMER_SIZE*sizeof(double);
   nsize += MESH_TIMER_SIZE*sizeof(long);
   nsize += ncells*3*sizeof(int);
   return(nsize);
}

void Mesh::store_checkpoint(Crux *crux)
{
   // Need ncells for memory allocation
   //int storage = mesh_memory.get_memory_capacity(level);
   int storage = ncells;
   crux->store_named_ints("storage", 8, &storage, 1);
   // Write scalars to arrays for storing in checkpoint
   int int_vals[num_int_vals];

   int_vals[ 0] = CRUX_MESH_VERSION;
   int_vals[ 1] = ndim;
   int_vals[ 2] = levmx;

   // These are for values that will be different on every processor
   int int_dist_vals[num_int_dist_vals];
   int_dist_vals[ 0] = (int)ncells;
   int_dist_vals[ 1] = (int)ncells_ghost;
   int_dist_vals[ 2] = offtile_local_count;

   double double_vals[num_double_vals];

   double_vals[0] = offtile_ratio_local;

   //int flags = RESTART_DATA;
   // Now add memory entries to database for storing checkpoint
   //mesh_memory.memory_add(int_dist_vals, (size_t)num_int_dist_vals, 4, "mesh_int_dist_vals", flags);
   //flags = RESTART_DATA | REPLICATED_DATA;
   //mesh_memory.memory_add(int_vals, (size_t)num_int_vals, 4, "mesh_int_vals", flags);
   //mesh_memory.memory_add(double_vals, (size_t)num_double_vals, 8, "mesh_double_vals", flags);

   //flags = RESTART_DATA;
   //mesh_memory.memory_add(cpu_counters, (size_t)MESH_COUNTER_SIZE, 4, "mesh_cpu_counters", flags);
   //mesh_memory.memory_add(gpu_counters, (size_t)MESH_COUNTER_SIZE, 4, "mesh_gpu_counters", flags);

   //mesh_memory.memory_add(cpu_timers, (size_t)MESH_TIMER_SIZE, 8, "mesh_cpu_timers", flags);
   //mesh_memory.memory_add(gpu_timers, (size_t)MESH_TIMER_SIZE, 8, "mesh_gpu_timers", flags);

   // Store MallocPlus memory database
   //crux->store_MallocPlus(mesh_memory);
   crux->store_int_array(int_dist_vals, num_int_dist_vals);
   crux->store_replicated_int_array(int_vals, num_int_vals);
   crux->store_double_array(double_vals, num_double_vals);
   crux->store_int_array(cpu_counters, MESH_COUNTER_SIZE);
   crux->store_int_array(gpu_counters, MESH_COUNTER_SIZE);

   crux->store_double_array(cpu_timers, MESH_TIMER_SIZE);
   crux->store_long_array(gpu_timers, MESH_TIMER_SIZE);

   crux->store_int_array(i, ncells);
   crux->store_int_array(j, ncells);
#ifdef REG_INTEGER
   crux->store_int_array((int *)level, ncells);
#elif SHORT_INTEGER
   crux->store_short_array((short *)level, ncells);
#elif MIN_INTEGER
   crux->store_char_array((char *)level, ncells);
#endif

   // Remove memory entries from database now that data is stored
   //mesh_memory.memory_remove(int_dist_vals);
   //mesh_memory.memory_remove(int_vals);
   //mesh_memory.memory_remove(double_vals);
   //mesh_memory.memory_remove(cpu_counters);
   //mesh_memory.memory_remove(gpu_counters);
   //mesh_memory.memory_remove(cpu_timers);
   //mesh_memory.memory_remove(gpu_timers);
}

void Mesh::restore_checkpoint(Crux *crux)
{
   int storage;
   crux->restore_named_ints("storage", 8, &storage, 1);

   // Create memory for reading data into
   int int_dist_vals[num_int_dist_vals];
   int int_vals[num_int_vals];
   double double_vals[num_double_vals];

   mesh_memory.memory_delete(nlft);
   mesh_memory.memory_delete(nrht);
   mesh_memory.memory_delete(nbot);
   mesh_memory.memory_delete(ntop);
   // needs to cast char_t to void so doesn't mistake it for string
   mesh_memory.memory_delete((void *)celltype);

   nlft = NULL;
   nrht = NULL;
   ntop = NULL;
   nbot = NULL;
   celltype = NULL;

   mesh_memory.memory_delete(i);
   mesh_memory.memory_delete(j);
   mesh_memory.memory_delete(level);
   mesh_memory.memory_delete(xface_i);
   mesh_memory.memory_delete(xface_j);
   mesh_memory.memory_delete(xface_level);
   mesh_memory.memory_delete(map_xface2cell_lower);
   mesh_memory.memory_delete(map_xface2cell_upper);
   mesh_memory.memory_delete(map_xcell2face_left1);
   mesh_memory.memory_delete(map_xcell2face_left2);
   mesh_memory.memory_delete(map_xcell2face_right1);
   mesh_memory.memory_delete(map_xcell2face_right2);
   mesh_memory.memory_delete(xrecvIdx);
   mesh_memory.memory_delete(xrecvCIdx);
   mesh_memory.memory_delete(xplusCell2Idx);
   mesh_memory.memory_delete(xminusCell2Idx);
   mesh_memory.memory_delete(xsendIdx1);
   mesh_memory.memory_delete(xsendIdx2);
   mesh_memory.memory_delete(yface_i);
   mesh_memory.memory_delete(yface_j);
   mesh_memory.memory_delete(yface_level);
   mesh_memory.memory_delete(map_yface2cell_lower);
   mesh_memory.memory_delete(map_yface2cell_upper);
   mesh_memory.memory_delete(map_ycell2face_bot1);
   mesh_memory.memory_delete(map_ycell2face_bot2);
   mesh_memory.memory_delete(map_ycell2face_top1);
   mesh_memory.memory_delete(map_ycell2face_top2);
   mesh_memory.memory_delete(yrecvIdx);
   mesh_memory.memory_delete(yrecvCIdx);
   mesh_memory.memory_delete(yplusCell2Idx);
   mesh_memory.memory_delete(yminusCell2Idx);
   mesh_memory.memory_delete(ysendIdx1);
   mesh_memory.memory_delete(ysendIdx2);
   //allocate(ncells);

   // Resize is a mesh method
   // resize(storage);
   allocate (storage);
   memory_reset_ptrs();
   
   //int flags = RESTART_DATA;
   // Now add memory entries to database for restoring checkpoint
   //mesh_memory.memory_add(int_dist_vals, (size_t)num_int_dist_vals, 4, "mesh_int_dist_vals", flags);
   //flags = RESTART_DATA | REPLICATED_DATA;
   //mesh_memory.memory_add(int_vals, (size_t)num_int_vals, 4, "mesh_int_vals", flags);
   //mesh_memory.memory_add(double_vals, (size_t)num_double_vals, 8, "mesh_double_vals", flags);

   // Restore MallocPlus memory database
   //crux->restore_MallocPlus(mesh_memory);
   crux->restore_int_array(int_dist_vals, num_int_dist_vals);
   ncells                    = int_dist_vals[ 0];
   crux->restore_replicated_int_array(int_vals, num_int_vals);
   crux->restore_double_array(double_vals, num_double_vals);
   crux->restore_int_array(cpu_counters, MESH_COUNTER_SIZE);
   crux->restore_int_array(gpu_counters, MESH_COUNTER_SIZE);

   crux->restore_double_array(cpu_timers, MESH_TIMER_SIZE);
   crux->restore_long_array(gpu_timers, MESH_TIMER_SIZE);

   crux->restore_int_array(i, ncells);
   crux->restore_int_array(j, ncells);
#ifdef REG_INTEGER
   crux->restore_int_array((int *)level, ncells);
#elif SHORT_INTEGER
   crux->restore_short_array((short *)level, ncells);
#elif MIN_INTEGER
   crux->restore_char_array((char *)level, ncells);
#endif

   index.clear();
   index.resize(ncells);
   for (uint ic = 0; ic < ncells; ic++){
     index[ic] = ic;
   }

   calc_celltype(ncells);

   // Remove memory entries from database now that data is restored
   //mesh_memory.memory_remove(int_dist_vals);
   //mesh_memory.memory_remove(int_vals);
   //mesh_memory.memory_remove(double_vals);

   // Check version number
   if (int_vals[ 0] != CRUX_MESH_VERSION) {
      printf("CRUX version mismatch for mesh data, version on file is %d, version in code is %d\n",
         int_vals[0], CRUX_MESH_VERSION);
      exit(0);
   }

   // Copy out scalar values from array
   ncells                    = int_dist_vals[ 0];
   ncells_ghost              = int_dist_vals[ 1];
   offtile_local_count       = int_dist_vals[ 2];

   // Copy out scalar values from array
   ndim                      = int_vals[ 1];
   levmx                     = int_vals[ 2];

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS && mype == 0) {
      const char *int_dist_vals_descriptor[num_int_dist_vals] = {
         "ncells",
         "ncells_ghost",
         "offtile_local_count"
      };
      const char *int_vals_descriptor[num_int_vals] = {
         "CRUX_MESH_VERSION",
         "ndim",
         "levmx",
      };
      printf("\n");
      printf("       === Restored mesh int_dist_vals ===\n");
      for (int i = 0; i < num_int_dist_vals; i++){
         printf("       %-30s %d\n",int_dist_vals_descriptor[i], int_dist_vals[i]);
      }
      printf("       === Restored mesh int_vals ===\n");
      for (int i = 0; i < num_int_vals; i++){
         printf("       %-30s %d\n",int_vals_descriptor[i], int_vals[i]);
      }
      printf("       === Restored mesh int_vals ===\n");
      printf("\n");
   }
#endif

   offtile_ratio_local = double_vals[0];

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS && mype == 0) {
      const char *double_vals_descriptor[num_double_vals] = {
         "offtile_ratio_local"
      };
      printf("\n");
      printf("       === Restored mesh double_vals ===\n");
      for (int i = 0; i < num_double_vals; i++){
         printf("       %-30s %lf\n",double_vals_descriptor[i], double_vals[i]);
      }
      printf("       === Restored mesh double_vals ===\n");
      printf("\n");
   }
#endif

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS && mype == 0) {
      printf("       === Restored mesh cpu counters ===\n");
      for (int i = 0; i < MESH_COUNTER_SIZE; i++){
         printf("       %-30s %d\n",mesh_counter_descriptor[i], cpu_counters[i]);
      }
      printf("       === Restored mesh cpu counters ===\n");
      printf("       === Restored mesh gpu counters ===\n");
      for (int i = 0; i < MESH_COUNTER_SIZE; i++){
         printf("       %-30s %d\n",mesh_counter_descriptor[i], gpu_counters[i]);
      }
      printf("       === Restored mesh gpu counters ===\n");
      printf("\n");
   }
#endif

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS && mype == 0) {
      printf("       === Restored mesh cpu timers ===\n");
      for (int i = 0; i < MESH_TIMER_SIZE; i++){
         printf("       %-30s %lf\n",mesh_timer_descriptor[i], cpu_timers[i]);
      }
      printf("       === Restored mesh cpu timers ===\n");
      printf("\n");
   }
#endif

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS && mype == 0) {
      printf("\n");
      printf("       === Restored mesh gpu timers ===\n");
      for (int i = 0; i < MESH_TIMER_SIZE; i++){
         printf("       %-30s %lld\n",mesh_timer_descriptor[i], gpu_timers[i]);
      }
      printf("       === Restored mesh gpu timers ===\n");
      printf("\n");
   }
#endif
   //calc_celltype(ncells);
}


// This code due to Matt Calef
void scan ( scanInt *input , scanInt *output , scanInt length) 
{
#ifdef _OPENMP
   // This already assumes it is in a parallel region

   // Get the total number of threads

   scanInt numThreads = omp_get_num_threads ( );

   // Compute the range for which this thread is responsible.

   scanInt threadID = omp_get_thread_num ( );
   scanInt start = length * ( threadID     ) / numThreads;
   scanInt end   = length * ( threadID + 1 ) / numThreads;

   // In the case that there are fewer entries than threads, some
   // threads will have no entries.  Only perform this operation if
   // there is a postive number of entries.

   if ( start < end ) {

       // Do a scan over the region for this thread, with an initial
       // value of zero.

       output[start] = 0;
       for ( scanInt i = start + 1 ; i < end ; i++ ) {
          output[i] = output[i-1] + input[i-1];
        }
   }
    
   // Wait until all threads get here. 

#pragma omp barrier
    
   // At this point each thread has done an independent scan of its
   // region.  All scans, except the first, are off by an
   // offset. Here we have a single thread compute that offset with a
   // serial scan that strides over the regions assigned to each
   // thread.

#pragma omp single
   for ( scanInt i = 1 ; i < numThreads ; i ++ ) {
      scanInt s0 = length * ( i - 1 ) / numThreads;
      scanInt s1 = length * ( i     ) / numThreads;

      if ( s0 < s1 ) 
         output[s1] = output[s0] + input[s1-1];

      if ( s0 < s1 - 1 )
         output[s1] += output[s1-1];
   }

   // Barrier is implicit from omp single Wait until all threads get here. 

   // Apply the offset to the range for this thread.
    
#pragma omp simd
   for ( scanInt i = start + 1 ; i < end ; i++ ) {
      output[i] += output[start];
           //if (output[i] == 284) printf("%d\n", i);
    }

#else
   output[0] = 0;
   for (int ic = 0; ic < length; ic++){
      output[ic+1] = output[ic] + input[ic];
   }
#endif
}
/****************************************************//**
**********************************/
void Mesh::get_bounds(int& lowerBound, int& upperBound){
#ifdef _OPENMP
        int threadID = omp_get_thread_num();
	lowerBound = lowerBound_Global[threadID];
	upperBound = upperBound_Global[threadID];
//	printf("GETBOUNDs ThreadID: %d, upperBound: %d, lowerBound: %d \n",threadID, upperBound, lowerBound);
#else
	lowerBound = 0;
	upperBound = ncells;
#endif
}

/****************************************************//**
*SETTING BOUNDS!!!!!!****
**********************************/
void Mesh::set_bounds(int n){

#ifdef _OPENMP
      //  #pragma omp parallel
        {
        int nthreads = omp_get_num_threads();//Private for each thread
        int threadID = omp_get_thread_num(); //Private for each thread
#pragma omp master 
     	{
        	if(lowerBound_Global == NULL) lowerBound_Global = (int *)malloc(nthreads*sizeof(int)); 
        	if(upperBound_Global == NULL) upperBound_Global = (int *)malloc(nthreads*sizeof(int)); 
        }
//#pragma omp flush (lowerBound_Global, upperBound_Global)
#pragma omp barrier
 	
        int work = n/nthreads;
        if(threadID<(n%nthreads))work++;
        int lowerBound = ((n / nthreads)*threadID) + min(n%nthreads, threadID);
        int upperBound = lowerBound + work;
//      printf("ThreadID: %d, upperBound: %d, lowerBound: %d \n",threadID, upperBound, lowerBound);
        lowerBound_Global[threadID] = lowerBound;
        upperBound_Global[threadID] = upperBound;
        }
#else 
     	if(lowerBound_Global == NULL) lowerBound_Global = (int *)malloc(1*sizeof(int)); 
       	if(upperBound_Global == NULL) upperBound_Global = (int *)malloc(1*sizeof(int)); 
        int lowerBound = 0;
        int upperBound = ncells;
        lowerBound_Global[0] = lowerBound;
        upperBound_Global[0] = upperBound;
#endif

}
