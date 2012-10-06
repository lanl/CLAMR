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
#include <unistd.h>
#include "hsfc/hsfc.h"
#include "kdtree/KDTree.h"
#include "mesh.h"
#include "reorder.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif
#include "timer/timer.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "l7/l7.h"
#endif
#include "reduce.h"
#include "genmalloc/genmalloc.h"

#define DEBUG 0
#define NEIGHBOR_CHECK 0

#ifndef DEBUG
#define DEBUG 0
#endif

#define TIMING_LEVEL 2

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#ifdef HAVE_CL_DOUBLE
typedef double      real;
#define MPI_C_REAL MPI_DOUBLE
#define CONSERVATION_EPS    .02
#define STATE_EPS        .025
#ifdef HAVE_MPI
#define L7_REAL L7_DOUBLE
#endif
#else
typedef float       real;
#define MPI_C_REAL MPI_FLOAT
#define CONSERVATION_EPS    .1
#define STATE_EPS      15.0
#ifdef HAVE_MPI
#define L7_REAL L7_FLOAT
#endif
#endif

typedef unsigned int uint;

#define TWO 2
#define HALF 0.5

#define __NEW_STENCIL__
//#define __OLD_STENCIL__
//#define STENCIL_WARNING 1

#ifdef STENCIL_WARNING
int do_stencil_warning=1;
#else
int do_stencil_warning=0;
#endif

extern bool localStencil;
int calc_neighbor_type;

cl_kernel      kernel_hash_init;
cl_kernel      kernel_hash_init_corners;
cl_kernel      kernel_hash_setup;
cl_kernel      kernel_hash_setup_local;
cl_kernel      kernel_hash_setup_border;
cl_kernel      kernel_calc_neighbors;
cl_kernel      kernel_calc_neighbors_local;
cl_kernel      kernel_calc_border_cells;
cl_kernel      kernel_calc_border_cells2;
cl_kernel      kernel_calc_neighbors_local2;
cl_kernel      kernel_copy_mesh_data;
cl_kernel      kernel_copy_ghost_data;
cl_kernel      kernel_adjust_neighbors;
cl_kernel      kernel_reduction_scan;
cl_kernel      kernel_hash_size;
cl_kernel      kernel_finish_hash_size;
cl_kernel      kernel_calc_spatial_coordinates;
cl_kernel      kernel_count_BCs;

cl_kernel      kernel_do_load_balance_lower;
cl_kernel      kernel_do_load_balance_middle;
cl_kernel      kernel_do_load_balance_upper;

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

void Mesh::mesh_reorder(vector<int> iorder)
{
   assert(index.size() == ncells);
   assert(i.size() == ncells);
   assert(j.size() == ncells);
   assert(level.size() == ncells);
   assert(celltype.size() == ncells);
   assert(x.size() == ncells);
   assert(dx.size() == ncells);
   assert(y.size() == ncells);
   assert(dy.size() == ncells);
   assert(iorder.size() == ncells);

   reorder(index,   iorder);
   reorder(i,       iorder);
   reorder(j,       iorder);
   reorder(level,   iorder);
   reorder(celltype,iorder);
   reorder(x,       iorder);
   reorder(dx,      iorder);
   reorder(y,       iorder);
   reorder(dy,      iorder);

   if (nlft.size() >= ncells) {
      vector<int> inv_iorder(ncells);

      for (uint i = 0; i < ncells; i++)
      {  inv_iorder[iorder[i]] = i; }

      reorder_indexarray(nlft, iorder, inv_iorder);
      reorder_indexarray(nrht, iorder, inv_iorder);
      reorder_indexarray(nbot, iorder, inv_iorder);
      reorder_indexarray(ntop, iorder, inv_iorder);
   }
}

Mesh::Mesh(FILE *fin, int *numpe)
{
   char string[80];
   ibase = 1;

   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"levmax %d",&levmx);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"cells %ld",&ncells);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"numpe %d",numpe);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"ndim %d",&ndim);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"xaxis %lf %lf",(double*)&xmin, (double*)&deltax);
   if(fgets(string, 80, fin) == NULL) exit(-1);
   sscanf(string,"yaxis %lf %lf",(double*)&ymin, (double*)&deltay);
   if (ndim == 3){
     if(fgets(string, 80, fin) == NULL) exit(-1);
     sscanf(string,"zaxis %lf %lf",(double*)&zmin, (double*)&deltaz);
   }
   if(fgets(string, 80, fin) == NULL) exit(-1);

   index.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);
   level.resize(ncells);
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
   //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());
   printf("index orig index   i     j     lev   nlft  nrht  nbot  ntop   xlow    xhigh     ylow    yhigh\n");
   for (uint ic=0; ic<ncells; ic++)
   {  printf("%6d %6d   %4d  %4d   %4d  %4d  %4d  %4d  %4d ", ic, index[ic], i[ic], j[ic], level[ic], nlft[ic], nrht[ic], nbot[ic], ntop[ic]);
      printf("%8.2lf %8.2lf %8.2lf %8.2lf\n", x[ic], x[ic]+dx[ic], y[ic], y[ic]+dy[ic]); }
}

void Mesh::print_local()
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   if (nlft.size() >= ncells){
      fprintf(fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mype);
      for (uint ic=0; ic<ncells; ic++) {
         fprintf(fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mype,ic, ic+noffset,i[ic], j[ic], level[ic], nlft[ic], nrht[ic], nbot[ic], ntop[ic]);
      }
   } else {
      fprintf(fp,"%d:    index   i     j     lev\n",mype);
      for (uint ic=0; ic<ncells; ic++) {
         fprintf(fp,"%d: %6d  %4d  %4d   %4d  \n", mype,ic, i[ic], j[ic], level[ic]);
      }
   }
}

#ifdef HAVE_OPENCL
void Mesh::print_dev_local(cl_command_queue command_queue)
{
   vector<int>i_tmp(ncells_ghost);
   vector<int>j_tmp(ncells_ghost);
   vector<int>level_tmp(ncells_ghost);
   vector<int>nlft_tmp(ncells_ghost);
   vector<int>nrht_tmp(ncells_ghost);
   vector<int>nbot_tmp(ncells_ghost);
   vector<int>ntop_tmp(ncells_ghost);
   ezcl_enqueue_read_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &i_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &j_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level, CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &level_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nlft_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nrht_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &nbot_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &ntop_tmp[0], NULL);

   fprintf(fp,"\n%d:                    Printing mesh for dev_local\n\n",mype);

   fprintf(fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mype);
   for (uint ic=0; ic<ncells_ghost; ic++) {
      fprintf(fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mype,ic, ic+noffset,i_tmp[ic], j_tmp[ic], level_tmp[ic], nlft_tmp[ic], nrht_tmp[ic], nbot_tmp[ic], ntop_tmp[ic]);
   }
   fprintf(fp,"\n%d:              Finished printing mesh for dev_local\n\n",mype);
}

void Mesh::compare_dev_local_to_local(cl_command_queue command_queue)
{
   vector<int>i_tmp(ncells_ghost);
   vector<int>j_tmp(ncells_ghost);
   vector<int>level_tmp(ncells_ghost);
   vector<int>nlft_tmp(ncells_ghost);
   vector<int>nrht_tmp(ncells_ghost);
   vector<int>nbot_tmp(ncells_ghost);
   vector<int>ntop_tmp(ncells_ghost);
   ezcl_enqueue_read_buffer(command_queue, dev_i,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &i_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,     CL_FALSE, 0, ncells_ghost*sizeof(cl_int), &j_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, ncells_ghost*sizeof(cl_int), &level_tmp[0], NULL);
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

void Mesh::compare_neighbors_gpu_global_to_cpu_global(cl_command_queue command_queue)
{
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
   vector<int> &nlft_global = mesh_global->nlft;
   vector<int> &nrht_global = mesh_global->nrht;
   vector<int> &nbot_global = mesh_global->nbot;
   vector<int> &ntop_global = mesh_global->ntop;

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
#endif
}

#ifdef HAVE_OPENCL
void Mesh::compare_neighbors_all_to_gpu_local(cl_command_queue command_queue, Mesh *mesh_global, int *nsizes, int *ndispl)
//uint ncells_ghost, uint ncells_global, Mesh *mesh_global, int *nsizes, int *ndispl)
{
#ifdef HAVE_MPI
   size_t &ncells_global = mesh_global->ncells;
   vector<int> &nlft_global = mesh_global->nlft;
   vector<int> &nrht_global = mesh_global->nrht;
   vector<int> &nbot_global = mesh_global->nbot;
   vector<int> &ntop_global = mesh_global->ntop;

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
            mype,ic,nlft_global[ic],nlft_global[nlft_global[ic]],Test_global[nlft_global[nlft_global[ic]]]);               printf("%d:                         check %5d -- nlftl %5d nlftl nlftl %5d check  %5d\n",
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
#endif
}

void Mesh::compare_indices_gpu_global_to_cpu_global(cl_command_queue command_queue)
{
   vector<int> i_check(ncells);
   vector<int> j_check(ncells);
   vector<int> level_check(ncells);
   vector<int> celltype_check(ncells);
   /// Set read buffers for data.
   ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int), &i_check[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int), &j_check[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int), &level_check[0],    NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  0, ncells*sizeof(cl_int), &celltype_check[0], NULL);
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
   vector<int>   &celltype_global = mesh_global->celltype;
   vector<int>   &i_global        = mesh_global->i;
   vector<int>   &j_global        = mesh_global->j;
   vector<int>   &level_global    = mesh_global->level;

   vector<int> i_check_global(ncells_global);
   vector<int> j_check_global(ncells_global);
   vector<int> level_check_global(ncells_global);
   vector<int> celltype_check_global(ncells_global);

#ifdef HAVE_MPI
   MPI_Allgatherv(&celltype[0], nsizes[mype], MPI_INT, &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&i[0],        nsizes[mype], MPI_INT, &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&j[0],        nsizes[mype], MPI_INT, &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&level[0],    nsizes[mype], MPI_INT, &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
#endif

   for (uint ic = 0; ic < ncells_global; ic++){
      if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 3 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",cycle,ic,celltype_global[ic],celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 3 at cycle %d i_global & i_check_global %d %d  %d  \n",cycle,ic,i_global[ic],i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 3 at cycle %d j_global & j_check_global %d %d  %d  \n",cycle,ic,j_global[ic],j_check_global[ic]);
      if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 3 at cycle %d level_global & level_check_global %d %d  %d  \n",cycle,ic,level_global[ic],level_check_global[ic]);
   }
}

#ifdef HAVE_OPENCL
void Mesh::compare_indices_all_to_gpu_local(cl_command_queue command_queue, Mesh *mesh_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle)
{
#ifdef HAVE_MPI
   vector<int> &level_global = mesh_global->level;
   vector<int> &celltype_global = mesh_global->celltype;
   vector<int> &i_global = mesh_global->i;
   vector<int> &j_global = mesh_global->j;

   cl_mem &dev_celltype_global = mesh_global->dev_celltype;
   cl_mem &dev_i_global = mesh_global->dev_i;
   cl_mem &dev_j_global = mesh_global->dev_j;
   cl_mem &dev_level_global = mesh_global->dev_level;

   // Need to compare dev_H to H, etc
   vector<int> level_check(ncells);
   vector<int> celltype_check(ncells);
   vector<int> i_check(ncells);
   vector<int> j_check(ncells);
   /// Set read buffers for data.
   ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int),  &level_check[0],     NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int),  &celltype_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int),  &i_check[0],         NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_TRUE,  0, ncells*sizeof(cl_int),  &j_check[0],         NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (level[ic] != level_check[ic] )       printf("%d: DEBUG rezone 1 cell %d level %d level_check %d\n",mype, ic, level[ic], level_check[ic]);
      if (celltype[ic] != celltype_check[ic] ) printf("%d: DEBUG rezone 1 cell %d celltype %d celltype_check %d\n",mype, ic, celltype[ic], celltype_check[ic]);
      if (i[ic] != i_check[ic] )               printf("%d: DEBUG rezone 1 cell %d i %d i_check %d\n",mype, ic, i[ic], i_check[ic]);
      if (j[ic] != j_check[ic] )               printf("%d: DEBUG rezone 1 cell %d j %d j_check %d\n",mype, ic, j[ic], j_check[ic]);
   }

   // And compare dev_H gathered to H_global, etc
   vector<int>celltype_check_global(ncells_global);
   vector<int>i_check_global(ncells_global);
   vector<int>j_check_global(ncells_global);
   vector<int>level_check_global(ncells_global);
   MPI_Allgatherv(&celltype_check[0], nsizes[mype], MPI_INT,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&i_check[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&j_check[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&level_check[0],    nsizes[mype], MPI_INT,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (level_global[ic] != level_check_global[ic] )       printf("%d: DEBUG rezone 2 cell %d level_global %d level_check_global %d\n",mype, ic, level_global[ic], level_check_global[ic]);
      if (celltype_global[ic] != celltype_check_global[ic] ) printf("%d: DEBUG rezone 2 cell %d celltype_global %d celltype_check_global %d\n",mype, ic, celltype_global[ic], celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic] )               printf("%d: DEBUG rezone 2 cell %d i_global %d i_check_global %d\n",mype, ic, i_global[ic], i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic] )               printf("%d: DEBUG rezone 2 cell %d j_global %d j_check_global %d\n",mype, ic, j_global[ic], j_check_global[ic]);
   }

   // And compare H gathered to H_global, etc
   MPI_Allgatherv(&celltype[0], nsizes[mype], MPI_INT,    &celltype_check_global[0], &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&i[0],        nsizes[mype], MPI_INT,    &i_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&j[0],        nsizes[mype], MPI_INT,    &j_check_global[0],        &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   MPI_Allgatherv(&level[0],    nsizes[mype], MPI_INT,    &level_check_global[0],    &nsizes[0], &ndispl[0], MPI_INT,    MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 3 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",ncycle,ic,celltype_global[ic],celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 3 at cycle %d i_global & i_check_global %d %d  %d  \n",ncycle,ic,i_global[ic],i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 3 at cycle %d j_global & j_check_global %d %d  %d  \n",ncycle,ic,j_global[ic],j_check_global[ic]);
      if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 3 at cycle %d level_global & level_check_global %d %d  %d  \n",ncycle,ic,level_global[ic],level_check_global[ic]);
   }

   // Now the global dev_H_global to H_global, etc
   ezcl_enqueue_read_buffer(command_queue, dev_celltype_global, CL_FALSE, 0, ncells_global*sizeof(cl_int),  &celltype_check_global[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_i_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  &i_check_global[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_j_global,        CL_FALSE, 0, ncells_global*sizeof(cl_int),  &j_check_global[0],        NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_level_global,    CL_TRUE,  0, ncells_global*sizeof(cl_int),  &level_check_global[0],    NULL);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (celltype_global[ic] != celltype_check_global[ic])  printf("DEBUG rezone 4 at cycle %d celltype_global & celltype_check_global %d %d  %d  \n",ncycle,ic,celltype_global[ic],celltype_check_global[ic]);
      if (i_global[ic] != i_check_global[ic])                printf("DEBUG rezone 4 at cycle %d i_global & i_check_global %d %d  %d  \n",ncycle,ic,i_global[ic],i_check_global[ic]);
      if (j_global[ic] != j_check_global[ic])                printf("DEBUG rezone 4 at cycle %d j_global & j_check_global %d %d  %d  \n",ncycle,ic,j_global[ic],j_check_global[ic]);
      if (level_global[ic] != level_check_global[ic])        printf("DEBUG rezone 4 at cycle %d level_global & level_check_global %d %d  %d  \n",ncycle,ic,level_global[ic],level_check_global[ic]);
   }
#endif
}

void Mesh::compare_coordinates_gpu_global_to_cpu_global(cl_command_queue command_queue, cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy, cl_mem dev_H, real *H)
{
   vector<real>x_check(ncells);
   vector<real>dx_check(ncells);
   vector<real>y_check(ncells);
   vector<real>dy_check(ncells);
   vector<real>H_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_x,   CL_FALSE, 0, ncells*sizeof(cl_real), &x_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dx,  CL_FALSE, 0, ncells*sizeof(cl_real), &dx_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_y,   CL_FALSE, 0, ncells*sizeof(cl_real), &y_check[0],  NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_dy,  CL_FALSE, 0, ncells*sizeof(cl_real), &dy_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_H,   CL_TRUE,  0, ncells*sizeof(cl_real), &H_check[0],  NULL);
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

void Mesh::compare_coordinates_cpu_local_to_cpu_global(uint ncells_global, int *nsizes, int *ndispl, real *x, real *dx, real *y, real *dy, real *H, real *x_global, real *dx_global, real *y_global, real *dy_global, real *H_global, int cycle)
{
   vector<real> x_check_global(ncells_global);
   vector<real> dx_check_global(ncells_global);
   vector<real> y_check_global(ncells_global);
   vector<real> dy_check_global(ncells_global);
   vector<real> H_check_global(ncells_global);

#ifdef HAVE_MPI
   MPI_Allgatherv(&x[0],  nsizes[mype], MPI_C_REAL, &x_check_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&dx[0], nsizes[mype], MPI_C_REAL, &dx_check_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&y[0],  nsizes[mype], MPI_C_REAL, &y_check_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&dy[0], nsizes[mype], MPI_C_REAL, &dy_check_global[0], &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
   MPI_Allgatherv(&H[0],  nsizes[mype], MPI_C_REAL, &H_check_global[0],  &nsizes[0], &ndispl[0], MPI_C_REAL, MPI_COMM_WORLD);
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
void Mesh::compare_mpot_gpu_global_to_cpu_global(cl_command_queue command_queue, int *mpot, cl_mem dev_mpot)
{
   vector<int>mpot_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_mpot,  CL_TRUE,  0, ncells*sizeof(cl_int), &mpot_check[0], NULL);

   for (uint ic=0; ic<ncells; ic++) {
      if (mpot[ic] != mpot_check[ic]) printf("DEBUG -- mpot: ic %d mpot %d mpot_check %d\n",ic, mpot[ic], mpot_check[ic]);
   }
}
#endif

void Mesh::compare_mpot_cpu_local_to_cpu_global(uint ncells_global, int *nsizes, int *ndispl, int *mpot, int *mpot_global, int cycle)
{
   vector<int>mpot_save_global(ncells_global);
#ifdef HAVE_MPI
   MPI_Allgatherv(&mpot[0], ncells, MPI_INT, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
#endif
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 3 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,cycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }
   }

}

#ifdef HAVE_OPENCL
void Mesh::compare_mpot_all_to_gpu_local(cl_command_queue command_queue, int *mpot, int *mpot_global, cl_mem dev_mpot, cl_mem dev_mpot_global, uint ncells_global, int *nsizes, int *ndispl, int ncycle)
{
#ifdef HAVE_MPI
   // Need to compare dev_mpot to mpot 
   vector<int>mpot_save(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE,  0, ncells*sizeof(cl_int), &mpot_save[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (mpot[ic] != mpot_save[ic]) {
         printf("%d: DEBUG refine_potential 1 at cycle %d cell %d mpot & mpot_save %d %d \n",mype,ncycle,ic,mpot[ic],mpot_save[ic]);
      }    
   }    

   // Compare dev_mpot to mpot_global
   vector<int>mpot_save_global(ncells_global);
   MPI_Allgatherv(&mpot_save[0], nsizes[mype], MPI_INT, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 2 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,ncycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }    
   }    

   // Compare mpot to mpot_global
   MPI_Allgatherv(&mpot[0], nsizes[mype], MPI_INT, &mpot_save_global[0], &nsizes[0], &ndispl[0], MPI_INT, MPI_COMM_WORLD);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 3 at cycle %d cell %d mpot_global & mpot_save_global %d %d \n",mype,ncycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }    
   }    

   // Compare dev_mpot_global to mpot_global
   ezcl_enqueue_read_buffer(command_queue, dev_mpot_global, CL_TRUE,  0, ncells_global*sizeof(cl_int), &mpot_save_global[0], NULL);
   for (uint ic = 0; ic < ncells_global; ic++){
      if (mpot_global[ic] != mpot_save_global[ic]) {
         if (mype == 0) printf("%d: DEBUG refine_potential 4 at cycle %d cell %u mpot_global & mpot_save_global %d %d \n",mype,ncycle,ic,mpot_global[ic],mpot_save_global[ic]);
      }    
   }    
#endif
}

void Mesh::compare_ioffset_gpu_global_to_cpu_global(cl_command_queue command_queue, uint old_ncells, int *mpot, cl_mem dev_ioffset)
{
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

          if (celltype[ic] == REAL_CELL){
             mcount += mpot[ic] ? 4 : 1;
          } else {
             mcount += mpot[ic] ? 2 : 1;
          }
      }
      if (mtotal != ioffset_check[ig]) printf("DEBUG ig %d ioffset %d mcount %d\n",ig,ioffset_check[ig],mtotal);
      mtotal += mcount;
   }
}

void Mesh::compare_ioffset_all_to_gpu_local(cl_command_queue command_queue, uint old_ncells, uint old_ncells_global, int block_size, int block_size_global, int *mpot, int *mpot_global, cl_mem dev_ioffset, cl_mem dev_ioffset_global, int *ioffset, int *ioffset_global, int *celltype_global)
{
   // This compares ioffset for each block in the calculation
   ezcl_enqueue_read_buffer(command_queue, dev_ioffset, CL_TRUE, 0, block_size*sizeof(cl_int), &ioffset[0], NULL);
   int mtotal = 0; 
   for (uint ig=0; ig<(old_ncells+TILE_SIZE-1)/TILE_SIZE; ig++){
      int mcount = 0; 
      for (uint ic=ig*TILE_SIZE; ic<(ig+1)*TILE_SIZE; ic++){
          if (ic >= old_ncells) break;
          if (celltype[ic] == REAL_CELL) {
             mcount += mpot[ic] ? 4 : 1; 
          } else {
             mcount += mpot[ic] ? 2 : 1; 
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
          if (celltype_global[ic] == REAL_CELL) {
             mcount += mpot_global[ic] ? 4 : 1; 
          } else {
             mcount += mpot_global[ic] ? 2 : 1; 
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

Mesh::Mesh(int nx, int ny, int levmx_in, int ndim_in, int numpe_in, int boundary, int parallel_in, int do_gpu_calc)
{
   cpu_time_calc_neighbors     = 0.0;
      cpu_time_hash_setup      = 0.0;
      cpu_time_hash_query      = 0.0;
      cpu_time_find_boundary   = 0.0;
      cpu_time_gather_boundary = 0.0;
      cpu_time_hash_setup2     = 0.0;
      cpu_time_hash_query2     = 0.0;
      cpu_time_offtile_list    = 0.0;
      cpu_time_setup_comm      = 0.0;
      cpu_time_do_ghost_comm   = 0.0;

      cpu_time_kdtree_setup    = 0.0;
      cpu_time_kdtree_query    = 0.0;
   cpu_time_rezone_all         = 0.0;
   cpu_time_partition          = 0.0;
   cpu_time_calc_spatial_coordinates = 0.0;
   cpu_time_load_balance       = 0.0;

   gpu_time_calc_neighbors     = 0;
      gpu_time_hash_setup      = 0;
      gpu_time_hash_query      = 0;
      gpu_time_find_boundary   = 0;
      gpu_time_gather_boundary = 0;
      gpu_time_hash_setup2     = 0;
      gpu_time_hash_query2     = 0;
      gpu_time_offtile_list    = 0;
      gpu_time_setup_comm      = 0;
      gpu_time_do_ghost_comm   = 0;

      gpu_time_kdtree_setup    = 0;
      gpu_time_kdtree_query    = 0;
   gpu_time_rezone_all         = 0;
   gpu_time_count_BCs          = 0;
   gpu_time_calc_spatial_coordinates = 0;
   gpu_time_load_balance       = 0;

   cpu_rezone_counter       = 0;
   cpu_refine_smooth_counter = 0;
   cpu_calc_neigh_counter   = 0;
   cpu_load_balance_counter = 0;
   gpu_rezone_counter       = 0;
   gpu_refine_smooth_counter = 0;
   gpu_calc_neigh_counter   = 0;
   gpu_load_balance_counter = 0;

   ndim   = ndim_in;
   levmx  = levmx_in;

   offtile_ratio_local = 0;
   offtile_local_count = 1;

   mype  = 0;
   numpe = 1;
   ncells = 0;
   ncells_ghost = 0;
   parallel = parallel_in;
   
#ifdef HAVE_MPI
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init && parallel){
      MPI_Comm_rank(MPI_COMM_WORLD,&mype);
      MPI_Comm_size(MPI_COMM_WORLD,&numpe);
   }
#endif
   cell_handle = 0;

   deltax = 1.0;
   deltay = 1.0;

   have_boundary = boundary;

   int istart = 1,
       jstart = 1,
       iend   = nx,
       jend   = ny,
       nxx    = nx,
       nyy    = ny;
   imin = 0;
   jmin = 0;
   imax = nx+1;
   jmax = ny+1;
   if (have_boundary) {
      istart = 0;
      jstart = 0;
      iend   = nx + 1;
      jend   = ny + 1;
      nxx    = nx + 2;
      nyy    = ny + 2;
      imin   = 0;
      jmin   = 0;
      imax   = nx + 1;
      jmax   = ny + 1;
   }
   
   xmin = -deltax * 0.5 * (real)nxx;
   ymin = -deltay * 0.5 * (real)nyy;
   xmax =  deltax * 0.5 * (real)nxx;
   ymax =  deltay * 0.5 * (real)nyy;
   
   size_t lvlMxSize = levmx + 1;

   levtable.resize(lvlMxSize);
   lev_ibegin.resize(lvlMxSize);
   lev_jbegin.resize(lvlMxSize);
   lev_iend.resize(  lvlMxSize);
   lev_jend.resize(  lvlMxSize);
   lev_deltax.resize(lvlMxSize);
   lev_deltay.resize(lvlMxSize);
   
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
      levtable[lev] = (int)pow(2,lev);
   }

#ifdef HAVE_OPENCL
   // The copy host ptr flag will have the data copied to the GPU as part of the allocation
   if (do_gpu_calc) {
      dev_levtable = ezcl_malloc(&levtable[0],   const_cast<char *>("dev_levtable"), &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levdx    = ezcl_malloc(&lev_deltax[0], const_cast<char *>("dev_levdx"),    &lvlMxSize, sizeof(cl_real), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levdy    = ezcl_malloc(&lev_deltay[0], const_cast<char *>("dev_levdy"),    &lvlMxSize, sizeof(cl_real), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levibeg  = ezcl_malloc(&lev_ibegin[0], const_cast<char *>("dev_levibeg"),  &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_leviend  = ezcl_malloc(&lev_iend[0],   const_cast<char *>("dev_leviend"),  &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levjbeg  = ezcl_malloc(&lev_jbegin[0], const_cast<char *>("dev_levjbeg"),  &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levjend  = ezcl_malloc(&lev_jend[0],   const_cast<char *>("dev_levjend"),  &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
   }
#endif

   ibase = 0;

   int ncells_corners = 4;
   int i_corner[] = {   0,   0,imax,imax};
   int j_corner[] = {   0,jmax,   0,jmax};

   for(int ic=0; ic<ncells_corners; ic++){
      for (int    jj = j_corner[ic]*levtable[levmx]; jj < (j_corner[ic]+1)*levtable[levmx]; jj++) {
         for (int ii = i_corner[ic]*levtable[levmx]; ii < (i_corner[ic]+1)*levtable[levmx]; ii++) {
            corners_i.push_back(ii);
            corners_j.push_back(jj);
         }
      }
   }

}

#ifdef HAVE_OPENCL
void Mesh::init(int nx, int ny, double circ_radius, cl_context context, partition_method initial_order, int compute_device, int do_gpu_calc)
#else
void Mesh::init(int nx, int ny, double circ_radius, partition_method initial_order, int do_gpu_calc)
#endif
{
#ifdef HAVE_OPENCL
   if (do_gpu_calc) {
      if (compute_device == COMPUTE_DEVICE_ATI) printf("Starting compile of kernels in mesh\n");

      kernel_reduction_scan    = ezcl_create_kernel(context, "wave_kern_calc.cl", "finish_reduction_scan_cl", 0);
      kernel_hash_init         = ezcl_create_kernel(context, "wave_kern.cl",      "hash_init_cl",             0);
      kernel_hash_init_corners = ezcl_create_kernel(context, "wave_kern.cl",      "hash_init_corners_cl",      0);
      kernel_hash_setup        = ezcl_create_kernel(context, "wave_kern.cl",      "hash_setup_cl",            0);
      kernel_hash_setup_local  = ezcl_create_kernel(context, "wave_kern.cl",      "hash_setup_local_cl",      0);
      kernel_hash_setup_border = ezcl_create_kernel(context, "wave_kern.cl",      "hash_setup_border_cl",      0);
      kernel_calc_neighbors    = ezcl_create_kernel(context, "wave_kern.cl",      "calc_neighbors_cl",        0);
      kernel_calc_neighbors_local = ezcl_create_kernel(context, "wave_kern.cl",      "calc_neighbors_local_cl",        0);
      kernel_calc_border_cells = ezcl_create_kernel(context, "wave_kern.cl",      "calc_border_cells_cl",        0);
      kernel_calc_border_cells2 = ezcl_create_kernel(context, "wave_kern.cl",      "calc_border_cells2_cl",        0);
      kernel_calc_neighbors_local2 = ezcl_create_kernel(context, "wave_kern.cl",      "calc_neighbors_local2_cl",        0);
      kernel_copy_mesh_data = ezcl_create_kernel(context, "wave_kern.cl",      "copy_mesh_data_cl",        0);
      kernel_copy_ghost_data = ezcl_create_kernel(context, "wave_kern.cl",      "copy_ghost_data_cl",        0);
      kernel_adjust_neighbors = ezcl_create_kernel(context, "wave_kern.cl",      "adjust_neighbors_cl",        0);
      kernel_hash_size         = ezcl_create_kernel(context, "wave_kern.cl",      "calc_hash_size_cl",        0);
      kernel_finish_hash_size  = ezcl_create_kernel(context, "wave_kern.cl",      "finish_reduction_minmax4_cl",        0);
      kernel_calc_spatial_coordinates = ezcl_create_kernel(context, "wave_kern.cl",      "calc_spatial_coordinates_cl",        0);
      kernel_do_load_balance_lower  = ezcl_create_kernel(context, "wave_kern.cl", "do_load_balance_cl_lower",  0);
      kernel_do_load_balance_middle = ezcl_create_kernel(context, "wave_kern.cl", "do_load_balance_cl_middle", 0);
      kernel_do_load_balance_upper  = ezcl_create_kernel(context, "wave_kern.cl", "do_load_balance_cl_upper",  0);
      init_kernel_2stage_sum(context);
      init_kernel_2stage_sum_int(context);
      if (! have_boundary){
        kernel_count_BCs       = ezcl_create_kernel(context, "wave_kern.cl",      "count_BCs_cl",             0);
      }
      if (compute_device == COMPUTE_DEVICE_ATI) printf("Finishing compile of kernels in mesh\n");
   }
#endif

   KDTree_Initialize(&tree);

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

   if (ndim == 2) ncells = nxx * nyy - have_boundary * 4;
   else           ncells = nxx * nyy;

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

   index.resize(ncells);
   i.resize(ncells);
   j.resize(ncells);
   level.resize(ncells);

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

#ifdef HAVE_MPI
   if (parallel && numpe > 1) {
      int ncells_int = ncells;
      MPI_Allreduce(&ncells_int, &ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allgather(&ncells_int, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
      ndispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
      }
      noffset=ndispl[mype];
   }
#endif
      
   nlft.clear();
   nrht.clear();
   nbot.clear();
   ntop.clear();
   celltype.clear();

   partition_cells(numpe, index, initial_order);

   calc_celltype();
   calc_spatial_coordinates(0);

   //  Start lev loop here
   for (int ilevel=1; ilevel<=levmx; ilevel++) {

      int old_ncells = ncells;

#ifdef HAVE_MPI
      if (parallel && numpe > 1) {
         int ncells_int = ncells;
         MPI_Allreduce(&ncells_int, &ncells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   
         MPI_Allgather(&ncells_int, 1, MPI_INT, &nsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
         ndispl[0]=0;
         for (int ip=1; ip<numpe; ip++){
            ndispl[ip] = ndispl[ip-1] + nsizes[ip-1];
         }
         noffset=ndispl[mype];
      }
#endif
      
      ncells_ghost = ncells;
      calc_neighbors_local();

      kdtree_setup();

      int nez;
      vector<int> ind(ncells);

      KDTree_QueryCircleIntersect(&tree, &nez, &(ind[0]), circ_radius, ncells, &x[0], &dx[0], &y[0], &dy[0]);

      vector<int> mpot(ncells_ghost,0);

      for (int ic=0; ic<nez; ++ic){
         if (level[ind[ic]] < levmx) mpot[ind[ic]] = 1;
      }

      KDTree_Destroy(&tree);
      //  Refine the cells.
      int new_ncells = refine_smooth(mpot);
      int add_ncells = new_ncells - old_ncells;

      rezone_all(mpot, add_ncells);

      calc_spatial_coordinates(0);

      ncells = new_ncells;
   
   }  // End lev loop here

   index.clear();

   int ncells_corners = 4;
   int i_corner[] = {   0,   0,imax,imax};
   int j_corner[] = {   0,jmax,   0,jmax};

   for(int ic=0; ic<ncells_corners; ic++){
      for (int    jj = j_corner[ic]*levtable[levmx]; jj < (j_corner[ic]+1)*levtable[levmx]; jj++) {
         for (int ii = i_corner[ic]*levtable[levmx]; ii < (i_corner[ic]+1)*levtable[levmx]; ii++) {
            corners_i.push_back(ii);
            corners_j.push_back(jj);
         }
      }
   }
   ncells_ghost = ncells;
}

size_t Mesh::refine_smooth(vector<int> &mpot)
{
   int nl, nr, nt, nb;
   int nlt, nrt, ntr, nbr;

   int icount = rezone_count(mpot);

   int newcount = icount;
   int newcount_global = icount;
#ifdef HAVE_MPI
   if (parallel) {
      MPI_Allreduce(&newcount, &newcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   }
#endif

   if(newcount_global > 0 && levmx > 1) {
      int lev, ll, lr, lt, lb;
      int llt, lrt, ltr, lbr;
      int levcount = 1;
      size_t my_ncells=ncells;
      if (parallel) my_ncells=ncells_ghost;

      cpu_refine_smooth_counter++;

      vector<int> mpot_old(my_ncells);

      while (newcount_global > 0 && levcount < levmx){
         levcount++; 
         newcount=0;

         mpot.swap(mpot_old);

#ifdef HAVE_MPI
         if (numpe > 1) {
            L7_Update(&mpot_old[0], L7_INT, cell_handle);
         }
#endif

         for(uint ic = 0; ic < ncells; ic++) {
            lev = level[ic];
            mpot[ic] = mpot_old[ic];
            if(mpot_old[ic] > 0) continue;
   
            nl = nlft[ic];
            ll = level[nl];
            if(mpot_old[nl] > 0) ll++;
   
            if(ll - lev > 1) {
               mpot[ic]++;
               newcount++;
               continue;
            }

            ll = level[nl];
            if (ll > lev) {
               nlt = ntop[nl];
               llt = level[nlt];
               if(mpot_old[nlt] > 0) llt++;

               if(llt - lev > 1) {
                  mpot[ic]++;
                  newcount++;
                  continue;
               }
            }

            nr = nrht[ic];
            lr = level[nr];
            if(mpot_old[nr] > 0) lr++;
   
            if(lr - lev > 1) {
               mpot[ic]++;
               newcount++;
               continue;
            }

            lr = level[nr];
            if (lr > lev) {
               nrt = ntop[nr];
               lrt = level[nrt];
               if(mpot_old[nrt] > 0) lrt++;

               if(lrt - lev > 1) {
                  mpot[ic]++;
                  newcount++;
                  continue;
               }
            }

            nt = ntop[ic];
            lt = level[nt];
            if(mpot_old[nt] > 0) lt++;
   
            if(lt - lev > 1) {
               mpot[ic]++;
               newcount++;
               continue;
            }

            lt = level[nt];
            if (lt > lev) {
               ntr = nrht[nt];
               ltr = level[ntr];
               if(mpot_old[ntr] > 0) ltr++;

               if(ltr - lev > 1) {
                  mpot[ic]++;
                  newcount++;
                  continue;
               }
            }

            nb = nbot[ic];
            lb = level[nb];
            if(mpot_old[nb] > 0) lb++;
   
            if(lb - lev > 1) {
               mpot[ic]++;
               newcount++;
               continue;
            }

            lb = level[nb];
            if (lb > lev) {
               nbr = nrht[nb];
               lbr = level[nbr];
               if(mpot_old[nbr] > 0) lbr++;

               if(lbr - lev > 1) {
                  mpot[ic]++;
                  newcount++;
                  continue;
               }
            }
         }

         newcount_global = newcount;
         icount += newcount;
         
#ifdef HAVE_MPI
         if (parallel) {
            MPI_Allreduce(&newcount, &newcount_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         }
#endif

//         printf("%d: newcount is %d levmx %d\n",levcount,newcount,levmx);
      //} while (newcount > 0 && levcount < 10);
      } while (newcount_global > 0 && levcount < levmx);
   }

#ifdef HAVE_MPI
   if (numpe > 1) {
      L7_Update(&mpot[0], L7_INT, cell_handle);
  }
#endif

   for (uint ic=0; ic<ncells; ic++) {
      if (celltype[ic] < 0) {
         switch (celltype[ic]) {
            case LEFT_BOUNDARY:
               if (mpot[nrht[ic]] > 0) {
                  mpot[ic]++;
                  icount++;
               }
               break;
            case RIGHT_BOUNDARY:
               if (mpot[nlft[ic]] > 0) {
                  mpot[ic]++;
                  icount++;
               }
               break;
            case BOTTOM_BOUNDARY:
               if (mpot[ntop[ic]] > 0) {
                  mpot[ic]++;
                  icount++;
               }
               break;
            case TOP_BOUNDARY:
               if (mpot[nbot[ic]] > 0) {
                  mpot[ic]++;
                  icount++;
               }
               break;
         }
      }
   }

   icount = 0;

   for (uint ic=0; ic<ncells; ++ic){
      if (mpot[ic] < 0) {
         if (celltype[ic] == REAL_CELL) {
            icount -= 3;
         } else {
            icount --;
         }
      }// XXX When coarsening is introduced, the issue arises that icount could become 0 even though
       // XXX both refinements and coarsenings will be needed
      if (mpot[ic] > 0) {
         //printf("mpot[%d] = %d\n",ic,mpot[ic]);
         if (celltype[ic] == REAL_CELL){
            icount += 3;
         } else {
            icount ++;
         }
      }
   }
   //printf("icount is %d\n",icount);

   return(ncells+icount);
}

#ifdef HAVE_OPENCL
void Mesh::terminate(void)
{
      ezcl_device_memory_remove(dev_levtable);
      ezcl_device_memory_remove(dev_levdx);
      ezcl_device_memory_remove(dev_levdy);
      ezcl_device_memory_remove(dev_levibeg);
      ezcl_device_memory_remove(dev_leviend);
      ezcl_device_memory_remove(dev_levjbeg);
      ezcl_device_memory_remove(dev_levjend);

      ezcl_device_memory_remove(dev_level);
      ezcl_device_memory_remove(dev_i);
      ezcl_device_memory_remove(dev_j);
      ezcl_device_memory_remove(dev_celltype);

      ezcl_kernel_release(kernel_reduction_scan);
      ezcl_kernel_release(kernel_hash_init);
      ezcl_kernel_release(kernel_hash_init_corners);
      ezcl_kernel_release(kernel_hash_setup);
      ezcl_kernel_release(kernel_hash_setup_local);
      ezcl_kernel_release(kernel_hash_setup_border);
      ezcl_kernel_release(kernel_calc_neighbors);
      ezcl_kernel_release(kernel_calc_neighbors_local);
      ezcl_kernel_release(kernel_calc_border_cells);
      ezcl_kernel_release(kernel_calc_border_cells2);
      ezcl_kernel_release(kernel_calc_neighbors_local2);
      ezcl_kernel_release(kernel_copy_mesh_data);
      ezcl_kernel_release(kernel_copy_ghost_data);
      ezcl_kernel_release(kernel_adjust_neighbors);
      ezcl_kernel_release(kernel_hash_size);
      ezcl_kernel_release(kernel_finish_hash_size);
      ezcl_kernel_release(kernel_calc_spatial_coordinates);
      ezcl_kernel_release(kernel_do_load_balance_lower);
      ezcl_kernel_release(kernel_do_load_balance_middle);
      ezcl_kernel_release(kernel_do_load_balance_upper);
      terminate_kernel_2stage_sum();
      terminate_kernel_2stage_sum_int();
      if (! have_boundary){
        ezcl_kernel_release(kernel_count_BCs);
      }
}
#endif

void Mesh::rezone_spread(vector<int> &mpot)
{
   for (uint ic=0; ic<ncells; ++ic) {
      if (mpot[ic] > 0) continue;
      if (mpot[nlft[ic]] == 1 && level[ic] <= level[nlft[ic]]) mpot[ic] = 2;
      if (mpot[nrht[ic]] == 1 && level[ic] <= level[nrht[ic]]) mpot[ic] = 2;
      if (mpot[nbot[ic]] == 1 && level[ic] <= level[nbot[ic]]) mpot[ic] = 2;
      if (mpot[ntop[ic]] == 1 && level[ic] <= level[ntop[ic]]) mpot[ic] = 2;
   }


   for (uint ic=0; ic<ncells; ++ic){
      if (is_left_boundary(ic)   && mpot[nrht[ic]] > 0) mpot[ic]=3;
      if (is_right_boundary(ic)  && mpot[nlft[ic]] > 0) mpot[ic]=3;
      if (is_bottom_boundary(ic) && mpot[ntop[ic]] > 0) mpot[ic]=3;
      if (is_top_boundary(ic)    && mpot[nbot[ic]] > 0) mpot[ic]=3;
   }
}

int Mesh::rezone_count(vector<int> mpot)
{
   int icount=0;

   for (uint ic=0; ic<ncells; ++ic){
      if (mpot[ic] < 0) {
         if (celltype[ic] == REAL_CELL) {
            icount -= 3;
         } else {
            icount --;
         }
      }// XXX When coarsening is introduced, the issue arises that icount could become 0 even though
       // XXX both refinements and coarsenings will be needed
      if (mpot[ic] > 0) {
         //printf("mpot[%d] = %d\n",ic,mpot[ic]);
         if (celltype[ic] == REAL_CELL){
            icount += 3;
         } else {
            icount ++;
         }
      }
   }
   //printf("icount is %d\n",icount);

   return(icount);
}

#ifdef HAVE_OPENCL
void Mesh::gpu_rezone_count(cl_command_queue command_queue, size_t block_size, size_t local_work_size,
    cl_mem dev_ioffset, cl_mem &dev_result)
{
     /*
     __kernel void finish_reduction_scan_cl(
                       const    int   isize,    // 0
              __global          int  *ioffset,  // 1
              __global          int  *result,   // 2
              __local           int  *tile,     // 3
              __local           int  *tile)     // 4
     */
   ezcl_set_kernel_arg(kernel_reduction_scan, 0, sizeof(cl_int),  (void *)&block_size);
   ezcl_set_kernel_arg(kernel_reduction_scan, 1, sizeof(cl_mem),  (void *)&dev_ioffset);
   ezcl_set_kernel_arg(kernel_reduction_scan, 2, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_reduction_scan, 3, local_work_size*sizeof(cl_int),    NULL);
   ezcl_set_kernel_arg(kernel_reduction_scan, 4, local_work_size*sizeof(cl_int),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_scan, 1, NULL, &local_work_size, &local_work_size, NULL);
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
   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   x.resize(ncells);
   dx.resize(ncells);
   y.resize(ncells);
   dy.resize(ncells);

   if (have_boundary) {
      for (uint ic = 0; ic < ncells; ic++) {
         x[ic]  = xmin + lev_deltax[level[ic]] * (real)(i[ic] - ibase);
         dx[ic] =        lev_deltax[level[ic]];
         y[ic]  = ymin + lev_deltay[level[ic]] * (real)(j[ic] - ibase);
         dy[ic] =        lev_deltay[level[ic]];
      }
   } else {
      for (uint ic = 0; ic < ncells; ic++) {
         x[ic]  = xmin + lev_deltax[level[ic]] * (real)(i[ic] - lev_ibegin[level[ic]]);
         dx[ic] =        lev_deltax[level[ic]];
         y[ic]  = ymin + lev_deltay[level[ic]] * (real)(j[ic] - lev_jbegin[level[ic]]);
         dy[ic] =        lev_deltay[level[ic]];
      }
   }

   cpu_time_calc_spatial_coordinates += cpu_timer_stop(tstart_cpu);
}

#ifdef HAVE_OPENCL
void Mesh::gpu_calc_spatial_coordinates(cl_command_queue command_queue, cl_mem dev_x, cl_mem dev_dx, cl_mem dev_y, cl_mem dev_dy)
{
   cl_event calc_spatial_coordinates_event;

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
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  1, sizeof(cl_real),   (void *)&xmin);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  2, sizeof(cl_real),   (void *)&ymin);
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

   gpu_time_calc_spatial_coordinates += ezcl_timer_calc(&calc_spatial_coordinates_event, &calc_spatial_coordinates_event);
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
   if (ndim > 2) {
      for (uint ic=0; ic<ncells; ic++){
         if (z[ic] < zmin) zmin = z[ic];
      }
   }

   xmax=-1.0e30, ymax=-1.0e30, zmax=-1.0e30;
   real xhigh, yhigh, zhigh;

   for (uint ic=0; ic<ncells; ic++){
      xhigh = x[ic]+dx[ic];
      if (xhigh > xmax) xmax = xhigh;
   }
   for (uint ic=0; ic<ncells; ic++){
      yhigh = y[ic]+dy[ic];
      if (yhigh > ymax) ymax = yhigh;
   }
   if (ndim > 2) {
      for (uint ic=0; ic<ncells; ic++){
        zhigh = z[ic]+dz[ic];
        if (zhigh > zmax) zmax = zhigh;
      }
   }

#ifdef HAVE_MPI
   if (parallel) {
      real xmin_global,xmax_global,ymin_global,ymax_global;
      MPI_Allreduce(&xmin, &xmin_global, 1, MPI_C_REAL, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&xmax, &xmax_global, 1, MPI_C_REAL, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&ymin, &ymin_global, 1, MPI_C_REAL, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&ymax, &ymax_global, 1, MPI_C_REAL, MPI_MAX, MPI_COMM_WORLD);
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
   real xmid, ymid, zmid;

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
   if (ndim > 2) {
      for (uint ic=0; ic<ncells; ic++){
         zmid = z[ic]+0.5*dz[ic];
         if (zmid < zcentermin) zcentermin = zmid;
         if (zmid > zcentermax) zcentermax = zmid;
      }
   }

#ifdef HAVE_MPI
   if (parallel) {
      real xcentermin_global,xcentermax_global,ycentermin_global,ycentermax_global;
      MPI_Allreduce(&xcentermin, &xcentermin_global, 1, MPI_C_REAL, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&xcentermax, &xcentermax_global, 1, MPI_C_REAL, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&ycentermin, &ycentermin_global, 1, MPI_C_REAL, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&ycentermax, &ycentermax_global, 1, MPI_C_REAL, MPI_MAX, MPI_COMM_WORLD);
      xcentermin = xcentermin_global;
      xcentermax = xcentermax_global;
      ycentermin = ycentermin_global;
      ycentermax = ycentermax_global;
   }
#endif

}

void Mesh::rezone_all(vector<int> mpot, int add_ncells)
{
   struct timeval tstart_cpu;

   uint ic,          //  Index for old cell arrays.
        nc;          //  Index for new cell arrays.
   int set_index = 0;

   cpu_timer_start(&tstart_cpu);

   cpu_rezone_counter++;

   //  Check for requested mesh refinements; if there are none, return.
   if (parallel) {
#ifdef HAVE_MPI
      int global_add_ncells;
      MPI_Allreduce(&add_ncells, &global_add_ncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (global_add_ncells == 0) set_index = 1;
#endif
   } else {
      if (add_ncells == 0) set_index = 1;
   }
      
   //  Clear variables which are now outdated and no longer needed.
   index.empty();

   if (set_index == 1) {
      index.resize(ncells);
      for (uint ic=0; ic<ncells; ic++){
         index[ic]=ic;
      }
      return;
   }

   int new_ncells = ncells + add_ncells;

   //  Initialize variables for outdated information which is still needed.
   vector<int> i_old(    new_ncells);
   vector<int> j_old(    new_ncells);
   vector<int> level_old(new_ncells);

   i.swap(i_old);
   j.swap(j_old);
   level.swap(level_old);

   index.resize(new_ncells);

   //  Insert new cells into the mesh at the point of refinement.
   vector<int> order(4,    -1), //  Vector of refined mesh traversal order; set to -1 to indicate errors.
               invorder(4, -1); //  Vector mapping location from base index.

   int ifirst      = 0;
   int ilast       = 0;
   int jfirst      = 0;
   int jlast       = 0;
   int level_first = 0;
   int level_last  = 0;

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

      MPI_Isend(&level_old[ncells-1], 1,MPI_INT,next,1,MPI_COMM_WORLD,req+8);
      MPI_Irecv(&level_first,         1,MPI_INT,prev,1,MPI_COMM_WORLD,req+9);

      MPI_Isend(&level_old[0],        1,MPI_INT,prev,1,MPI_COMM_WORLD,req+10);
      MPI_Irecv(&level_last,          1,MPI_INT,next,1,MPI_COMM_WORLD,req+11);

      MPI_Waitall(12, req, status);
#endif
   }

   for (ic = 0, nc = 0; ic < ncells; ic++)
   {  index[nc] = -1;
      
      if (mpot[ic] == 0)
      {  //  No change is needed; copy the old cell straight to the new mesh at this location.
         i[nc]     = i_old[ic];
         j[nc]     = j_old[ic];
         level[nc] = level_old[ic];
         index[nc] = ic;
         nc++;
      } //  Complete no change needed.
      
      else if (mpot[ic] < 0)
      {  //  Coarsening is needed; remove this cell and the following three and replace them with one.  TODO
         i[nc] = i_old[ic]/2;
         j[nc] = j_old[ic]/2;
         if        (i[nc]*2  ==i_old[ic] && j[nc]*2  ==j_old[ic]) {
            index[nc]=ic;
         } else if (i[nc]*2 + 1 == i_old[ic] && j[nc]*2     == j_old[ic]) {
         } else if (i[nc]*2     == i_old[ic] && j[nc]*2 + 1 == j_old[ic]) {
         } else if (i[nc]*2 + 1 == i_old[ic] && j[nc]*2 + 1 == j_old[ic]) {
            level[nc] = level_old[ic] - 1;
            nc++; // XXX TODO Updating of nc assumes the upper right cell is the last to be visited. Also, is it guaranteed that all four cells needing coarsening will have mpot < 0 ? XXX
         }
      } //  Complete coarsening needed.
      
      else if (mpot[ic] > 0)
      {  //  Refinement is needed; insert four cells where once was one.
         if (celltype[ic] == REAL_CELL)
         {  
            if (localStencil) {
               //  Store the coordinates of the cells before and after this one on
               //  the space-filling curve index.

#ifdef __OLD_STENCIL__
               real  nx[3],  //  x-coordinates of cells.
                     ny[3];  //  y-coordinates of cells.
               if (ic != 0) {
                  nx[0] = lev_deltax[level_old[ic-1]] * (real)i_old[ic-1];
                  ny[0] = lev_deltay[level_old[ic-1]] * (real)j_old[ic-1];
               } else {
                  nx[0] = lev_deltax[level_first] * (real)ifirst;
                  ny[0] = lev_deltay[level_first] * (real)jfirst;
               }
               nx[1] = lev_deltax[level_old[ic  ]] * (real)i_old[ic  ];
               ny[1] = lev_deltay[level_old[ic  ]] * (real)j_old[ic  ];
               if (ic != ncells-1) {
                  nx[2] = lev_deltax[level_old[ic+1]] * (real)i_old[ic+1];
                  ny[2] = lev_deltay[level_old[ic+1]] * (real)j_old[ic+1];
               } else {
                  nx[2] = lev_deltax[level_last] * (real)ilast;
                  ny[2] = lev_deltay[level_last] * (real)jlast;
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

               int crf = levtable[levmx-level_old[ic]];

               if (ic != 0) {
                  ir[0] = i_old[ic - 1] * levtable[levmx-level_old[ic - 1]];
                  jr[0] = j_old[ic - 1] * levtable[levmx-level_old[ic - 1]];
               } else {
                  //printf("%d cell %d is a first\n",mype,ic);
                  ir[0] = ifirst * levtable[levmx-level_first];
                  jr[0] = jfirst * levtable[levmx-level_first];
               }
               ir[1] = i_old[ic    ] * levtable[levmx-level_old[ic    ]];
               jr[1] = j_old[ic    ] * levtable[levmx-level_old[ic    ]];
               if (ic != ncells-1) {
                  ir[2] = i_old[ic + 1] * levtable[levmx-level_old[ic + 1]];
                  jr[2] = j_old[ic + 1] * levtable[levmx-level_old[ic + 1]];
               } else {
                  //printf("%d cell %d is a last\n",mype,ic);
                  ir[2] = ilast * levtable[levmx-level_last];
                  jr[2] = jlast * levtable[levmx-level_last];
               }
               //printf("%d cell %d %d nc %d ir/jr %d %d %d %d %d %d\n",mype,ic,ic+noffset,nc,ir[0],jr[0],ir[1],jr[1],ir[2],jr[2]);

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
            
            //  Create the cells in the correct order and orientation.
            for (int ii = 0; ii < 4; ii++)
            {  level[nc] = level_old[ic] + 1;
               index[nc] = ic;
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
            index[nc] = ic;
            nc++;
            
            // upper
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;
         }
         else if (celltype[ic] == RIGHT_BOUNDARY) {
            // lower
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;
            
            // upper
            i[nc] = i_old[ic]*2;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;
         }
         else if (celltype[ic] == BOTTOM_BOUNDARY) {
            // left
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;
            
            // right
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2 + 1;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;
         }
         else if (celltype[ic] == TOP_BOUNDARY) {
            // right
            i[nc] = i_old[ic]*2 + 1;
            j[nc] = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;

            // left
            i[nc]  = i_old[ic]*2;
            j[nc]  = j_old[ic]*2;
            level[nc] = level_old[ic] + 1;
            index[nc] = ic;
            nc++;
         }
      } //  Complete refinement needed.
   } //  Complete addition of new cells to the mesh.

   ncells = nc;
   
   calc_celltype();

   cpu_time_rezone_all += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_neighbors(void)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cpu_calc_neigh_counter++;

   nlft.resize(ncells);
   nrht.resize(ncells);
   nbot.resize(ncells);
   ntop.resize(ncells);

   if (calc_neighbor_type == HASH_TABLE) {

      struct timeval tstart_lev2;
      if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

      int jmaxsize = (jmax+1)*levtable[levmx];
      int imaxsize = (imax+1)*levtable[levmx];
      int **hash = (int **)genmatrix(jmaxsize,imaxsize,sizeof(int));

      for (int jj = 0; jj<jmaxsize; jj++){
         for (int ii = 0; ii<imaxsize; ii++){
            hash[jj][ii]=-1;
         }
      }

      for(uint ic=0; ic<ncells; ic++){
         int lev = level[ic];
         if (i[ic] < lev_ibegin[lev]) { // left boundary
            for (int    jj = j[ic]*levtable[levmx-lev]; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
               for (int ii = 0; ii < (i[ic]+1)*levtable[levmx-lev]; ii++) {
                  hash[jj][ii] = ic;
               }
            }
         } else if (i[ic] > lev_iend[lev]) { // right boundary
            for (int    jj = j[ic]*levtable[levmx-lev]; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]; ii < (imax+1)*levtable[levmx]; ii++) {
                  hash[jj][ii] = ic;
               }
            }
         } else if (j[ic] < lev_jbegin[lev]) { // bottom boundary
            for (int    jj = 0; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]; ii < (i[ic]+1)*levtable[levmx-lev]; ii++) {
                  hash[jj][ii] = ic;
               }
            }
         } else if (j[ic] > lev_jend[lev]) { // top boundary
            for (int    jj = j[ic]*levtable[levmx-lev]; jj < (jmax+1)*levtable[levmx]; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]; ii < (i[ic]+1)*levtable[levmx-lev]; ii++) {
                  hash[jj][ii] = ic;
               }
            }
         } else if (lev == levmx) {
            hash[j[ic]][i[ic]] = ic;
         } else {
/* Original Write to Hash Table
            for (int jj = j[ic]*levtable[levmx-lev]; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
               for (int ii=i[ic]*levtable[levmx-lev]; ii<(i[ic]+1)*levtable[levmx-lev]; ii++) {
                  hash[jj][ii] = ic;
               }
            }
*/
/* Optimization: Always writes to max of 7 hash buckets, 4 if cell is l=levmx-1 */
            int wid = levtable[levmx-lev];
            int jj = j[ic]*wid;
            int ii = i[ic]*wid;
            hash[jj][ii] = ic;
            ii += wid/2;
            hash[jj][ii] = ic;
            if(wid > 2) {
               ii = ii + wid/2 - 1;
               hash[jj][ii] = ic;
               ii = ii - wid/2 + 1;
            }
            ii -= wid/2;
            jj += wid/2;
            hash[jj][ii] = ic;
            ii = ii + wid - 1;
            hash[jj][ii] = ic;

            if(wid > 2) {
               ii = ii - wid + 1;
               jj = jj + wid/2 - 1;
               hash[jj][ii] = ic;
               ii += wid/2;
               hash[jj][ii] = ic;
            }
         }
      }

      if (TIMING_LEVEL >= 2) {
         cpu_time_hash_setup += cpu_timer_stop(tstart_lev2);
         cpu_timer_start(&tstart_lev2);
      }

      int ii, jj, lev, iii, jjj, levmult;

      for (uint ic=0; ic<ncells; ic++){
         ii = i[ic];
         jj = j[ic];
         lev = level[ic];
         levmult = levtable[levmx-lev];
         nlft[ic] = hash[      jj   *levmult               ][max(  ii   *levmult-1, 0         )];
         nrht[ic] = hash[      jj   *levmult               ][min( (ii+1)*levmult,   imaxsize-1)];
         nbot[ic] = hash[max(  jj   *levmult-1, 0)         ][      ii   *levmult               ];
         ntop[ic] = hash[min( (jj+1)*levmult,   jmaxsize-1)][      ii   *levmult               ];
      }

      for (uint ic=0; ic<ncells; ic++){
         if (nlft[ic] < 0){
            lev = level[ic];
            ii = i[ic];
            jj = j[ic];
            levmult = levtable[levmx-lev];
            iii = max( ii*levmult-1, 0);
            jjj = jj*levmult;
            if ( (jjj < 1*levtable[levmx] || jjj > (jmax-1)*levtable[levmx] ) && iii < 1*levtable[levmx] ) iii = ii*levmult;
            nlft[ic] = hash[jjj][iii];
         }
         if (nrht[ic] < 0){
            lev = level[ic];
            ii = i[ic];
            jj = j[ic];
            levmult = levtable[levmx-lev];
            iii = min( (ii+1)*levmult, imaxsize-1);
            jjj = jj*levmult;
            if ( (jjj < 1*levtable[levmx] || jjj > jmax*levtable[levmx]-1 ) && iii > imax*levtable[levmx]-1 ) iii = ii*levmult;
            nrht[ic] = hash[jjj][iii];
         }
         if (nbot[ic] < 0) {
            lev = level[ic];
            ii = i[ic];
            jj = j[ic];
            levmult = levtable[levmx-lev];
            iii = ii*levmult;
            jjj = max( jj*levmult-1, 0);
            if ( (iii < 1*levtable[levmx] || iii > (imax-1)*levtable[levmx] ) && jjj < 1*levtable[levmx] ) jjj = jj*levmult;
            nbot[ic] = hash[jjj][iii];
         }
         if (ntop[ic] < 0) {
            lev = level[ic];
            ii = i[ic];
            jj = j[ic];
            levmult = levtable[levmx-lev];
            iii = ii*levmult;
            jjj = min( (jj+1)*levmult, jmaxsize-1);
            if ( (iii < 1*levtable[levmx] || iii > imax*levtable[levmx]-1 ) && jjj > jmax*levtable[levmx]-1 ) jjj = jj*levmult;
            ntop[ic] = hash[jjj][iii];
         }
      }
 
      if (DEBUG) {
         printf("\n                                    HASH numbering\n");
         for (int jj = jmaxsize-1; jj>=0; jj--){
            printf("%4d:",jj);
            for (int ii = 0; ii<imaxsize; ii++){
               printf("%5d",hash[jj][ii]);
            }
            printf("\n");
         }
         printf("     ");
         for (int ii = 0; ii<imaxsize; ii++){
            printf("%4d:",ii);
         }
         printf("\n");
   
         printf("\n                                    nlft numbering\n");
         for (int jj = jmaxsize-1; jj>=0; jj--){
            printf("%4d:",jj);
            for (int ii = 0; ii<imaxsize; ii++){
               if (hash[jj][ii] >= 0) {
                  printf("%5d",nlft[hash[jj][ii]]);
               } else {
                  printf("     ");
               }
            }
            printf("\n");
         }
         printf("     ");
         for (int ii = 0; ii<imaxsize; ii++){
            printf("%4d:",ii);
         }
         printf("\n");
   
         printf("\n                                    nrht numbering\n");
         for (int jj = jmaxsize-1; jj>=0; jj--){
            printf("%4d:",jj);
            for (int ii = 0; ii<imaxsize; ii++){
               if (hash[jj][ii] >= 0) {
                  printf("%5d",nrht[hash[jj][ii]]);
               } else {
                  printf("     ");
               }
            }
            printf("\n");
         }
         printf("     ");
         for (int ii = 0; ii<imaxsize; ii++){
            printf("%4d:",ii);
         }
         printf("\n");

         printf("\n                                    nbot numbering\n");
         for (int jj = jmaxsize-1; jj>=0; jj--){
            printf("%4d:",jj);
            for (int ii = 0; ii<imaxsize; ii++){
               if (hash[jj][ii] >= 0) {
                  printf("%5d",nbot[hash[jj][ii]]);
               } else {
                  printf("     ");
               }
            }
            printf("\n");
         }
         printf("     ");
         for (int ii = 0; ii<imaxsize; ii++){
            printf("%4d:",ii);
         }
         printf("\n");

         printf("\n                                    ntop numbering\n");
         for (int jj = jmaxsize-1; jj>=0; jj--){
            printf("%4d:",jj);
            for (int ii = 0; ii<imaxsize; ii++){
               if (hash[jj][ii] >= 0) {
                  printf("%5d",ntop[hash[jj][ii]]);
               } else {
                  printf("     ");
               }
            }
            printf("\n");
         }
         printf("     ");
         for (int ii = 0; ii<imaxsize; ii++){
            printf("%4d:",ii);
         }
         printf("\n");
      }

      genmatrixfree((void **)hash);

      if (TIMING_LEVEL >= 2) cpu_time_hash_query += cpu_timer_stop(tstart_lev2);

   } else if (calc_neighbor_type == KDTREE) {

      struct timeval tstart_lev2;
      if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

      TBounds box;
      vector<int> index_list((int)pow(2,levmx*levmx) );

      int num;

      ibase = 0;
      calc_spatial_coordinates(ibase);

      kdtree_setup();

      if (TIMING_LEVEL >= 2) {
         cpu_time_kdtree_setup += cpu_timer_stop(tstart_lev2);
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

      if (TIMING_LEVEL >= 2) cpu_time_kdtree_query += cpu_timer_stop(tstart_lev2);

   } // calc_neighbor_type

   cpu_time_calc_neighbors += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_neighbors_local(void)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cpu_calc_neigh_counter++;

   nlft.resize(ncells,-98);
   nrht.resize(ncells,-98);
   nbot.resize(ncells,-98);
   ntop.resize(ncells,-98);

   if (calc_neighbor_type == HASH_TABLE) {

      struct timeval tstart_lev2;
      if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

      ncells_ghost = ncells;

      // Find maximum i column and j row for this processor
      int jminsize = (jmax+1)*levtable[levmx];
      int iminsize = (imax+1)*levtable[levmx];
      int jmaxsize = 0;
      int imaxsize = 0;
      for(uint ic=0; ic<ncells; ic++){
         int lev = level[ic];
         if ( j[ic]   *levtable[levmx-lev] < jminsize) jminsize =  j[ic]   *levtable[levmx-lev];
         if ((j[ic]+1)*levtable[levmx-lev] > jmaxsize) jmaxsize = (j[ic]+1)*levtable[levmx-lev];
         if ( i[ic]   *levtable[levmx-lev] < iminsize) iminsize =  i[ic]   *levtable[levmx-lev];
         if ((i[ic]+1)*levtable[levmx-lev] > imaxsize) imaxsize = (i[ic]+1)*levtable[levmx-lev];
      }
      //fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);

      // Expand size by 2*coarse_cells for ghost cells
      jminsize = max(jminsize-2*levtable[levmx],0);
      jmaxsize = min(jmaxsize+2*levtable[levmx],(jmax+1)*levtable[levmx]);
      iminsize = max(iminsize-2*levtable[levmx],0);
      imaxsize = min(imaxsize+2*levtable[levmx],(imax+1)*levtable[levmx]);
      //fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);

      // Allocate partial hash table
      int **hash = (int **)genmatrix(jmaxsize-jminsize,imaxsize-iminsize,sizeof(int));

      // Initialize hash table to -1
      for (int jj = 0; jj<(jmaxsize-jminsize); jj++){
         for (int ii = 0; ii<(imaxsize-iminsize); ii++){
            hash[jj][ii]=-1;
         }
      }

      //printf("%d: DEBUG -- noffset %d cells %d\n",mype,noffset,ncells);

      // Setting corners to -99
      int ncells_corners = 4;
      int i_corner[] = {   0,   0,imax,imax};
      int j_corner[] = {   0,jmax,   0,jmax};

      for(int ic=0; ic<ncells_corners; ic++){
         for (int    jj = j_corner[ic]*levtable[levmx]-jminsize; jj < (j_corner[ic]+1)*levtable[levmx]-jminsize; jj++) {
            for (int ii = i_corner[ic]*levtable[levmx]-iminsize; ii < (i_corner[ic]+1)*levtable[levmx]-iminsize; ii++) {
               //printf("%d: block j %d i %d\n",mype,jj,ii);
               if (jj < 0 || jj >= (jmaxsize-jminsize) || ii < 0 || ii >= (imaxsize-iminsize) ) continue;
               hash[jj][ii] = -99;
            }
         }
      }

      // Walk through cell array and set hash to global cell values
      for(uint ic=0; ic<ncells; ic++){
         int lev = level[ic];
         if (i[ic] < lev_ibegin[lev]) { // left boundary
            for (int    jj = j[ic]*levtable[levmx-lev]-jminsize; jj < (j[ic]+1)*levtable[levmx-lev]-jminsize; jj++) {
               for (int ii = 0; ii < (i[ic]+1)*levtable[levmx-lev]-iminsize; ii++) {
                  hash[jj][ii] = ic+noffset;
               }
            }
         } else if (i[ic] > lev_iend[lev]) { // right boundary
            for (int    jj = j[ic]*levtable[levmx-lev]-jminsize; jj < (j[ic]+1)*levtable[levmx-lev]-jminsize; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]-iminsize; ii < (imax+1)*levtable[levmx]-iminsize; ii++) {
                  hash[jj][ii] = ic+noffset;
               }
            }
         } else if (j[ic] < lev_jbegin[lev]) { // bottom boundary
            for (int    jj = 0; jj < (j[ic]+1)*levtable[levmx-lev]-jminsize; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]-iminsize; ii < (i[ic]+1)*levtable[levmx-lev]-iminsize; ii++) {
                  hash[jj][ii] = ic+noffset;
               }
            }
         } else if (j[ic] > lev_jend[lev]) { // top boundary
            for (int    jj = j[ic]*levtable[levmx-lev]-jminsize; jj < (jmax+1)*levtable[levmx]-jminsize; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]-iminsize; ii < (i[ic]+1)*levtable[levmx-lev]-iminsize; ii++) {
                  hash[jj][ii] = ic+noffset;
               }
            }
         } else if (lev == levmx) {
            //printf("%d: max j %d i %d\n",mype,j[ic]-jminsize,i[ic]-iminsize);
            hash[j[ic]-jminsize][i[ic]-iminsize] = ic+noffset;
         } else {
/* Original Write to Hash Table
            for (int    jj = j[ic]*levtable[levmx-lev]-jminsize; jj < (j[ic]+1)*levtable[levmx-lev]-jminsize; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]-iminsize; ii < (i[ic]+1)*levtable[levmx-lev]-iminsize; ii++) {
                  //printf("%d: block j %d i %d\n",mype,jj,ii);
                  hash[jj][ii] = ic+noffset;
               }
            }
*/
/* Optimization: Always writes to max of 7 hash buckets, 4 if cell is l=levmx-1 */
            int wid = levtable[levmx-lev];
            int jj = j[ic]*wid - jminsize;
            int ii = i[ic]*wid - iminsize;
            hash[jj][ii] = ic + noffset;
            ii += wid/2;
            hash[jj][ii] = ic + noffset;
            if(wid > 2) {
               ii = ii + wid/2 - 1;
               hash[jj][ii] = ic + noffset;
               ii = ii - wid/2 + 1;
            }
            ii -= wid/2;
            jj += wid/2;
            hash[jj][ii] = ic + noffset;
            ii = ii + wid - 1;
            hash[jj][ii] = ic + noffset;

            if(wid > 2) {
               ii = ii - wid + 1;
               jj = jj + wid/2 - 1;
               hash[jj][ii] = ic + noffset;
               ii += wid/2;
               hash[jj][ii] = ic + noffset;
            }
         }
      }
      
      if (TIMING_LEVEL >= 2) {
         cpu_time_hash_setup += cpu_timer_stop(tstart_lev2);
         cpu_timer_start(&tstart_lev2);
      }

      int ii, jj, lev, iii=0, jjj=0, levmult;

      // Set neighbors to global cell numbers from hash
      int jmaxcalc = (jmax+1)*levtable[levmx];
      int imaxcalc = (imax+1)*levtable[levmx];
      for (uint ic=0; ic<ncells; ic++){
         ii = i[ic];
         jj = j[ic];
         lev = level[ic];
         levmult = levtable[levmx-lev];
         nlft[ic] = hash[(      jj   *levmult               )-jminsize][(max(  ii   *levmult-1, 0         ))-iminsize];
         nrht[ic] = hash[(      jj   *levmult               )-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
         nbot[ic] = hash[(max(  jj   *levmult-1, 0)         )-jminsize][(      ii   *levmult               )-iminsize];
         ntop[ic] = hash[(min( (jj+1)*levmult,   jmaxcalc-1))-jminsize][(      ii   *levmult               )-iminsize];
         //fprintf(fp,"%d: DEBUG -- ic %d %d ntop[ic] %d \n",mype,ic,ic+noffset,ntop[ic]);
      }

      // Scan for corner boundary cells
      for (uint ic=0; ic<ncells; ic++){
         if (nlft[ic] == -99){
            lev = level[ic];
            levmult = levtable[levmx-lev];
            jjj = j[ic]*levmult;
            iii = i[ic]*levmult;
            nlft[ic] = hash[jjj-jminsize][iii-iminsize];
         }
         if (nrht[ic] == -99){
            lev = level[ic];
            levmult = levtable[levmx-lev];
            jjj = j[ic]*levmult;
            iii = i[ic]*levmult;
            nrht[ic] = hash[jjj-jminsize][iii-iminsize];
         }
         if (nbot[ic] == -99) {
            lev = level[ic];
            levmult = levtable[levmx-lev];
            iii = i[ic]*levmult;
            jjj = j[ic]*levmult;
            nbot[ic] = hash[jjj-jminsize][iii-iminsize];
         }
         if (ntop[ic] == -99) {
            lev = level[ic];
            levmult = levtable[levmx-lev];
            iii = i[ic]*levmult;
            jjj = j[ic]*levmult;
            ntop[ic] = hash[jjj-jminsize][iii-iminsize];
         }
      }

      if (TIMING_LEVEL >= 2) {
         cpu_time_hash_query += cpu_timer_stop(tstart_lev2);
         cpu_timer_start(&tstart_lev2);
      }

      if (numpe > 1) {
         vector<int> border_cell;

         // Push list of unsatisfied neighbor cells
         for (uint ic=0; ic<ncells; ic++){
            if (nlft[ic] == -1 || (level[nlft[ic]-noffset] > level[ic] && ntop[nlft[ic]-noffset] == -1) ){
               //fprintf(fp,"%d: Cell is %d nlft %d\n",mype,ic+noffset,nlft[ic]);
               border_cell.push_back(ic+noffset);
               if (nrht[ic] >= 0) {
                  border_cell.push_back(nrht[ic]);
                  if (level[nrht[ic]-noffset] > level[ic]) {
                     if (ntop[nrht[ic]-noffset] >= 0) border_cell.push_back(ntop[nrht[ic]-noffset]);
                  }
               }
            }
            if (nrht[ic] == -1 || (level[nrht[ic]-noffset] > level[ic] && ntop[nrht[ic]-noffset] == -1) ){
               //fprintf(fp,"%d: Cell is %d nrht %d\n",mype,ic+noffset,nrht[ic]);
               border_cell.push_back(ic+noffset);
               if (nlft[ic] >= 0) {
                  border_cell.push_back(nlft[ic]);
                  if (level[nlft[ic]-noffset] > level[ic]) {
                     if (ntop[nlft[ic]-noffset] >= 0) border_cell.push_back(ntop[nlft[ic]-noffset]);
                  }
               }
            }
            if (nbot[ic] == -1 || (level[nbot[ic]-noffset] > level[ic] && nrht[nbot[ic]-noffset] == -1) ) {
               //fprintf(fp,"%d: Cell is %d nbot %d\n",mype,ic+noffset,nbot[ic]);
               border_cell.push_back(ic+noffset);
               if (ntop[ic] >= 0) {
                  border_cell.push_back(ntop[ic]);
                  if (level[ntop[ic]-noffset] > level[ic]) {
                     if (nrht[ntop[ic]-noffset] >= 0) border_cell.push_back(nrht[ntop[ic]-noffset]);
                  }
               }
            }
            if (ntop[ic] == -1 || ( level[ntop[ic]-noffset] > level[ic] && nrht[ntop[ic]-noffset] == -1) ) {
               //fprintf(fp,"%d: Cell is %d ntop %d\n",mype,ic+noffset,ntop[ic]);
               border_cell.push_back(ic+noffset);
               if (nbot[ic] >= 0) {
                  border_cell.push_back(nbot[ic]);
                  if (level[nbot[ic]-noffset] > level[ic]) {
                     if (nrht[nbot[ic]-noffset] >= 0) border_cell.push_back(nrht[nbot[ic]-noffset]);
                  }
               }
            }
         }

         sort(border_cell.begin(),border_cell.end());
         vector<int>::iterator p_end = unique(border_cell.begin(),border_cell.end());

         vector<int> border_cell_num;

         for (vector<int>::iterator p=border_cell.begin(); p < p_end; p++){
            border_cell_num.push_back(*p);
         }
         //printf("%d: border cell size is %d\n",mype,border_cell_num.size());

         int nbsize_local = 0;
         for (vector<int>::iterator p=border_cell.begin(); p < p_end; p++){
            nbsize_local++;
         }

         int nbsize_global;

         vector<int> border_cell_i(nbsize_local);
         vector<int> border_cell_j(nbsize_local);
         vector<int> border_cell_level(nbsize_local);
    
         for (int ic = 0; ic <nbsize_local; ic++){
            int cell_num = border_cell_num[ic]-noffset;
            border_cell_i[ic] = i[cell_num]; 
            border_cell_j[ic] = j[cell_num]; 
            border_cell_level[ic] = level[cell_num]; 
            //fprintf(fp,"%d: Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num[ic],
            //   border_cell_i[ic],border_cell_j[ic],border_cell_level[ic]);
         }

         if (TIMING_LEVEL >= 2) {
            cpu_time_find_boundary += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         nbsize_global = nbsize_local;
#ifdef HAVE_MPI
         MPI_Allreduce(&nbsize_local, &nbsize_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

         vector<int> nbsizes(numpe);
         vector<int> nbdispl(numpe);
#ifdef HAVE_MPI
         MPI_Allgather(&nbsize_local, 1, MPI_INT, &nbsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
#endif
         nbdispl[0]=0;
         for (int ip=1; ip<numpe; ip++){
            nbdispl[ip] += nbdispl[ip-1] + nbsizes[ip-1];
         }

         vector<int>border_cell_num_global(nbsize_global);
         vector<int>border_cell_i_global(nbsize_global);
         vector<int>border_cell_j_global(nbsize_global);
         vector<int>border_cell_level_global(nbsize_global);

#ifdef HAVE_MPI
         MPI_Allgatherv(&border_cell_num[0],   nbsizes[mype], MPI_INT, &border_cell_num_global[0],   &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
         MPI_Allgatherv(&border_cell_i[0],     nbsizes[mype], MPI_INT, &border_cell_i_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
         MPI_Allgatherv(&border_cell_j[0],     nbsizes[mype], MPI_INT, &border_cell_j_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
         MPI_Allgatherv(&border_cell_level[0], nbsizes[mype], MPI_INT, &border_cell_level_global[0], &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
#endif

         //for (int ic = 0; ic < nbsize_global; ic++) {
         //   fprintf(fp,"%d: Global Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num_global[ic],
         //      border_cell_i_global[ic],border_cell_j_global[ic],border_cell_level_global[ic]);
         //}

         if (TIMING_LEVEL >= 2) {
            cpu_time_gather_boundary += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         int inew=0;
         for (int ic = 0; ic < nbsize_global; ic++) {
            int lev = border_cell_level_global[ic];
            int levmult = levtable[levmx-lev];
            //fprintf(fp,"%d: DEBUG cell %d i %d j %d\n",mype,ic,border_cell_i_global[ic],border_cell_j_global[ic]);
            if (border_cell_j_global[ic]*levmult < jminsize || border_cell_j_global[ic]*levmult >= jmaxsize) continue;
            if (border_cell_i_global[ic]*levmult < iminsize || border_cell_i_global[ic]*levmult >= imaxsize) continue;
            if (border_cell_num_global[ic] >= (int)noffset && border_cell_num_global[ic] < (int)(noffset+ncells)) continue;
            //fprintf(fp,"%d: ic is %d inew is %d\n",mype,ic,inew);
            if (inew != ic){
               border_cell_num_global[inew] = border_cell_num_global[ic];
               border_cell_i_global[inew] = border_cell_i_global[ic];
               border_cell_j_global[inew] = border_cell_j_global[ic];
               border_cell_level_global[inew] = border_cell_level_global[ic];
            }
            inew++;
         }
         nbsize_local = inew;


         //fprintf(fp,"%d: nbsize_local is %d\n",mype,nbsize_local);
         //for (int ic = 0; ic < nbsize_local; ic++) {
         //   fprintf(fp,"%d: Local Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num_global[ic],
         //      border_cell_i_global[ic],border_cell_j_global[ic],border_cell_level_global[ic]);
         //}

         // Walk through cell array and set hash to global cell values
         //fprintf(fp,"%d: DEBUG new hash jminsize %d jmaxsize %d iminsize %d imaxsize %d\n",mype,jminsize,jmaxsize,iminsize,imaxsize);
         for(int ic=0; ic<nbsize_local; ic++){
            //fprintf(fp,"%d: cell %d i %d j %d\n",mype,ic,border_cell_i_global[ic],border_cell_j_global[ic]);
            int lev = border_cell_level_global[ic];
            int levmult = levtable[levmx-lev];
            if (lev == levmx) {
               //fprintf(fp,"%d: cell %d max j %d i %d\n",mype,ic,j[ic]-jminsize,i[ic]-iminsize);
               hash[border_cell_j_global[ic]-jminsize][border_cell_i_global[ic]-iminsize] = border_cell_num_global[ic];
            } else {
               for (int    jj = max(border_cell_j_global[ic]*levmult-jminsize,0); jj < min((border_cell_j_global[ic]+1)*levmult,jmaxsize)-jminsize; jj++) {
                  for (int ii = max(border_cell_i_global[ic]*levmult-iminsize,0); ii < min((border_cell_i_global[ic]+1)*levmult,imaxsize)-iminsize; ii++) {
                     //fprintf(fp,"%d: cell %d block j %d i %d\n",mype,ic,jj,ii);
                     hash[jj][ii] = border_cell_num_global[ic];
                  }
               }
            }
         }

         if (TIMING_LEVEL >= 2) {
            cpu_time_hash_setup2 += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         // Set neighbors to global cell numbers from hash
         jmaxcalc = (jmax+1)*levtable[levmx];
         imaxcalc = (imax+1)*levtable[levmx];

         for (uint ic=0; ic<ncells; ic++){
            ii = i[ic];
            jj = j[ic];
            lev = level[ic];
            levmult = levtable[levmx-lev];
            //fprintf(fp,"%d:Neighbors input for ic %d ii %d jj %d levmult %d lev %d\n",mype,ic, ii, jj, levmult,lev);
            //fprintf(fp,"%d:Neighbors befor ic %d nlft %d nrht %d nbot %d ntop %d\n",mype,ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            if (nlft[ic] == -1) nlft[ic] = hash[(      jj   *levmult               )-jminsize][(max(  ii   *levmult-1, 0         ))-iminsize];
            if (celltype[ic] == BOTTOM_BOUNDARY && nlft[ic] == -1){
               if (nlft[ic] == -1) nlft[ic] = hash[(jj+1)*levmult-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
            }
            if (nrht[ic] == -1) nrht[ic] = hash[(      jj   *levmult               )-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
            if (celltype[ic] == BOTTOM_BOUNDARY && nrht[ic] == -1){
               if (nrht[ic] == -1) nrht[ic] = hash[(jj+1)*levmult-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
               //if (ic == 3 && mype == 0) printf("DEBUG line %d -- ic %d celltype %d nrht %d\n",__LINE__,ic,celltype[ic],nrht[ic]);
               //printf("DEBUG line %d -- ic %d celltype %d nrht %d jj %d ii %d\n",__LINE__,ic,celltype[ic],nrht[ic],(jj+1)*levmult-jminsize,(min( (ii+1)*levmult,   imaxcalc-1))-iminsize);
            }
            if (nbot[ic] == -1) nbot[ic] = hash[(max(  jj   *levmult-1, 0)         )-jminsize][(      ii   *levmult               )-iminsize];
            if (celltype[ic] == LEFT_BOUNDARY && ntop[ic] == -1){
               if (nbot[ic] == -1) nbot[ic] = hash[(max(  jj   *levmult-1, 0)         )-jminsize][(      ii   *levmult+1             )-iminsize];
            }
            if (ntop[ic] == -1) ntop[ic] = hash[(min( (jj+1)*levmult,   jmaxcalc-1))-jminsize][(      ii   *levmult               )-iminsize];
            if (celltype[ic] == LEFT_BOUNDARY && ntop[ic] == -1){
               if (ntop[ic] == -1) ntop[ic] = hash[(min( (jj+1)*levmult,   jmaxcalc-1))-jminsize][(      ii   *levmult+1             )-iminsize];
            }
            //fprintf(fp,"%d:Neighbors after ic %d nlft %d nrht %d nbot %d ntop %d\n",mype,ic,nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
         }

         // Scan for corner boundary cells
         for (uint ic=0; ic<ncells; ic++){
            if (nlft[ic] == -99){
               lev = level[ic];
               levmult = levtable[levmx-lev];
               jjj = j[ic]*levmult;
               iii = i[ic]*levmult;
               nlft[ic] = hash[jjj-jminsize][iii-iminsize];
            }
            if (nrht[ic] == -99){
               lev = level[ic];
               levmult = levtable[levmx-lev];
               jjj = j[ic]*levmult;
               iii = i[ic]*levmult;
               nrht[ic] = hash[jjj-jminsize][iii-iminsize];
            }
            if (nbot[ic] == -99) {
               lev = level[ic];
               levmult = levtable[levmx-lev];
               iii = i[ic]*levmult;
               jjj = j[ic]*levmult;
               nbot[ic] = hash[jjj-jminsize][iii-iminsize];
            }
            if (ntop[ic] == -99) {
               lev = level[ic];
               levmult = levtable[levmx-lev];
               iii = i[ic]*levmult;
               jjj = j[ic]*levmult;
               ntop[ic] = hash[jjj-jminsize][iii-iminsize];
            }
         }

         if (TIMING_LEVEL >= 2) {
            cpu_time_hash_query2 += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         if (DEBUG) {
            int jmaxglobal = (jmax+1)*levtable[levmx];
            int imaxglobal = (imax+1)*levtable[levmx];
            fprintf(fp,"\n                                    HASH 0 numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        fprintf(fp,"%5d",hash[jj-jminsize][ii-iminsize]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                           fprintf(fp,"%5d",nlft[hashval]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",nrht[hashval]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",nbot[hashval]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",ntop[hashval]);
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

         if (DEBUG) print_local();

         vector<int> offtile_list;

         int start_idx = noffset;
         int end_idx   = noffset+ncells;

         for (uint ic = 0; ic < ncells; ic++){
            ii = i[ic];
            jj = j[ic];
            lev = level[ic];

            int nl = nlft[ic];

            if (nl < start_idx || nl >= end_idx) {
               offtile_list.push_back(nl);
               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == nl) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d ic %d\n",mype, __LINE__, nl, ic+noffset);
#ifdef HAVE_MPI
                  L7_Terminate();
#endif
                  exit(-1);
               }
               // Get the left neighbor information
               int nlev = border_cell_level_global[ig];
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[max(levmx-nlev-1,0)];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];

               if (nlev == levmx) {
                  // Multiplying by levmult get the bottom left cell of the left neighbor. Subtract
                  // one from that to get the bottom right cell in the left left neighbor
                  int nll = hash[(jjj*levmult)-jminsize][max((iii-1)*levmult, 0)-iminsize];
                  if (nll != nl && nll > 0 && (nll < start_idx || nll >= end_idx) ) 
                     offtile_list.push_back(nll);
               } else {
                  int nll = hash[jjj*levmult-jminsize][max((iii*2-1)*levmultminus, 0)-iminsize];
                  if (nll != nl && nll > 0 && (nll < start_idx || nll >= end_idx) ) 
                     offtile_list.push_back(nll);

                  int ntll = hash[(jjj*2+1)*levmultminus-jminsize][max((iii*2-1)*levmultminus, 0)-iminsize];
                  if (ntll != nll && ntll > 0 && (ntll < start_idx || ntll >= end_idx) )
                     offtile_list.push_back(ntll);
               }
            }

            if (lev < levmx) {
               int levmultminus = levtable[max(levmx-lev-1,0)];
               int ntl = hash[min((jj*2+1)*levmultminus,jmaxcalc-1)-jminsize][max((ii*2-1)*levmultminus,0)-iminsize];
               if (ntl != nl && ntl > 0 && (ntl < start_idx || ntl >= end_idx) ) {
                  offtile_list.push_back(ntl);

                  int ig;
                  for (ig = 0; ig < nbsize_local; ig++) {
                     if (border_cell_num_global[ig] == ntl) break;
                  }
                  if (ig == nbsize_local) {
                     printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, ntl);
#ifdef HAVE_MPI
                     //L7_Terminate();
#endif
                     //exit(-1);
                  }
                  int nlev = border_cell_level_global[ig];
                  int iii = border_cell_i_global[ig];
                  int jjj = border_cell_j_global[ig];
                  int levmult      = levtable[levmx-nlev];
                  levmultminus = levtable[max(levmx-nlev-1,0)];

                  if (nlev == levmx) {
                     int nltl = hash[jjj*levmult-jminsize][max((iii-1)*levmult-iminsize,0)];
                     if (nltl != ntl && nltl > 0 && (nltl < start_idx || nltl >= end_idx) )
                        offtile_list.push_back(nltl);
                  } else {
                     // And we need both cells here
                     int nltl = hash[(jjj*2)*levmultminus-jminsize][max((iii*2-1)*levmultminus-iminsize,0)];
                     if (nltl != ntl && nltl > 0 && (nltl < start_idx || nltl >= end_idx) )
                        offtile_list.push_back(nltl);
                     int ntltl = hash[min((jjj*2+1)*levmultminus,jmaxcalc-1)-jminsize][max((iii*2-1)*levmultminus,0)-iminsize];
                     if (ntltl != nltl && ntltl > 0 && (ntltl < start_idx || ntltl >= end_idx))
                        offtile_list.push_back(ntltl);
                  }

               }
            }

            int nr = nrht[ic];

            if (nr < start_idx || nr >= end_idx) {
               offtile_list.push_back(nr);
               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == nr) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d for cell %d\n",mype, __LINE__, nr, ic+noffset);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               int nlev = border_cell_level_global[ig];
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[max(levmx-nlev-1,0)];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];

               if (nlev == levmx) {
                  // Multiplying by levmult gets the bottom left cell of the right neighbor. Add
                  // one to that to get the bottom left cell of the right right neighbor
                  int nrr = hash[(jjj*levmult)-jminsize][min((iii+1)*levmult,imaxcalc-1)-iminsize];
                  if (nrr != nr && nrr > 0 && (nrr < start_idx || nrr >= end_idx) )
                     offtile_list.push_back(nrr);
               } else {
                  int nrr = hash[(jjj*levmult)-jminsize][(min((iii+1)*levmult,imaxcalc-1))-iminsize];
                  if (nrr != nr && nrr > 0 && (nrr < start_idx || nrr >= end_idx) )
                     offtile_list.push_back(nrr);

                  int ntrr = hash[min((jjj*2+1)*levmultminus,jmaxcalc-1)-jminsize][min((iii+1)*levmult,imaxcalc-1)-iminsize];
                  if (ntrr != nrr && ntrr > 0 && (ntrr < start_idx || ntrr >= end_idx))
                     offtile_list.push_back(ntrr);
               }

            }

            if (lev < levmx) {
               int levmultminus = levtable[max(levmx-lev-1,0)];
               int ntr = hash[min((jj*2+1)*levmultminus,jmaxcalc-1)-jminsize][min((ii+1)*2*levmultminus,imaxcalc-1)-iminsize];
               if (ntr != nr && ntr > 0 && (ntr < start_idx || ntr >= end_idx) ) {
                  offtile_list.push_back(ntr);

                  int ig;
                  for (ig = 0; ig < nbsize_local; ig++) {
                     if (border_cell_num_global[ig] == ntr) break;
                  }
                  if (ig == nbsize_local) {
                     printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, ntr);
#ifdef HAVE_MPI
                     //L7_Terminate();
#endif
                     //exit(-1);
                  }
                  int nlev = border_cell_level_global[ig];
                  int iii = border_cell_i_global[ig];
                  int jjj = border_cell_j_global[ig];
                  int levmult      = levtable[levmx-nlev];
                  levmultminus = levtable[max(levmx-nlev-1,0)];

                  // Multiplying by levmult gets the bottom left cell of the top neighbor. Add
                  // one to that to get the bottom left cell of the top top neighbor

                  if (nlev == levmx) {
                     int nrtr = hash[jjj*levmult-jminsize][min((iii+1)*levmult-iminsize,imaxcalc-1)];
                     if (nrtr != ntr && nrtr > 0 && (nrtr < start_idx || nrtr >= end_idx) )
                        offtile_list.push_back(nrtr);
                  } else {
                     // And we need both cells here
                     int nrtr = hash[(jjj*2)*levmultminus-jminsize][min((iii+1)*2*levmultminus-iminsize,imaxcalc-1)];
                     if (nrtr != ntr && nrtr > 0 && (nrtr < start_idx || nrtr >= end_idx) )
                        offtile_list.push_back(nrtr);
                     int ntrtr = hash[min((jjj*2+1)*levmultminus,jmaxcalc-1)-jminsize][min((iii+1)*2*levmultminus,imaxcalc-1)-iminsize];
                     if (ntrtr != nrtr && ntrtr > 0 && (ntrtr < start_idx || ntrtr >= end_idx))
                        offtile_list.push_back(ntrtr);
                  }
               }
            }


            int nb = nbot[ic];

            if (nb < start_idx || nb >= end_idx) {
               offtile_list.push_back(nb);
               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == nb) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with bottom global cell match %d for cell %d\n",mype, __LINE__, nb, ic+noffset);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               // Get the bottom neighbor information
               int nlev = border_cell_level_global[ig];
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[max(levmx-nlev-1,0)];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];

               if (nlev == levmx) {
                  // Multiplying by levmult gets the bottom left cell of the bottom neighbor. Subtract
                  // one from that to get the top left cell in the bottom bottom neighbor
                  int nbb = hash[(max((jjj-1)*levmult, 0))-jminsize][(iii*levmult)-iminsize];
                  if (nbb != nb && nbb > 0 && (nbb < start_idx || nbb >= end_idx) )
                     offtile_list.push_back(nbb);
               } else {
                  int nbb = hash[max((jjj*2-1)*levmultminus, 0)-jminsize][iii*levmult-iminsize];
                  if (nbb != nb && nbb > 0 && (nbb < start_idx || nbb >= end_idx) )
                     offtile_list.push_back(nbb);
                  int nrbb = hash[max((jjj*2-1)*levmultminus, 0)-jminsize][(iii*2+1)*levmultminus-iminsize];
                  if (nrbb != nbb && nrbb > 0 && (nrbb < start_idx || nrbb >= end_idx) )
                     offtile_list.push_back(nrbb);
               }

            }

            if (lev < levmx) {
               int levmultminus = levtable[max(levmx-lev-1,0)];
               int nrb = hash[max((jj*2-1)*levmultminus,0)-jminsize][min((ii*2+1)*levmultminus,imaxcalc-1)-iminsize];
               if (nrb != nb && nrb > 0 && (nrb < start_idx || nrb >= end_idx) ) {
                  offtile_list.push_back(nrb);

                  int ig;
                  for (ig = 0; ig < nbsize_local; ig++) {
                     if (border_cell_num_global[ig] == nrb) break;
                  }
                  if (ig == nbsize_local) {
                     printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nrb);
#ifdef HAVE_MPI
                     //L7_Terminate();
#endif
                     //exit(-1);
                  }
                  int nlev = border_cell_level_global[ig];
                  int iii = border_cell_i_global[ig];
                  int jjj = border_cell_j_global[ig];
                  int levmult      = levtable[levmx-nlev];
                  levmultminus = levtable[max(levmx-nlev-1,0)];

                  if (nlev == levmx) {
                     int nbrb = hash[max((jjj-1)*levmult,0)-jminsize][iii*levmult-iminsize];
                     if (nbrb != nrb && nbrb > 0 && (nbrb < start_idx || nbrb >= end_idx) )
                        offtile_list.push_back(nbrb);
                  } else {
                     // And we need both cells here
                     int nbrb = hash[max((jjj*2-1)*levmultminus,0)-jminsize][(iii*2)*levmultminus-iminsize];
                     if (nbrb != nrb && nbrb > 0 && (nbrb < start_idx || nbrb >= end_idx) )
                        offtile_list.push_back(nbrb);
                     int nrbrb = hash[max((jjj*2-1)*levmultminus,0)-jminsize][min((iii*2+1)*levmultminus,imaxcalc-1)-iminsize];
                     if (nrbrb != nbrb && nrbrb > 0 && (nrbrb < start_idx || nrbrb >= end_idx))
                        offtile_list.push_back(nrbrb);
                  }

               }
            }

            int nt = ntop[ic];

            if (nt < start_idx || nt >= end_idx) {
               offtile_list.push_back(nt);
               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == nt) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nt);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               int nlev = border_cell_level_global[ig];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[max(levmx-nlev-1,0)];

               if (nlev == levmx) {
                  // Multiplying by levmult gets the bottom left cell of the top neighbor. Add
                  // one to that to get the bottom left cell of the top top neighbor
                  int ntt = hash[(min((jjj+1)*levmult,jmaxcalc-1))-jminsize][(iii*levmult)-iminsize];
                  if (ntt != nt && ntt > 0 && (ntt < start_idx || ntt >= end_idx) )
                     offtile_list.push_back(ntt);
               } else {
                  int ntt = hash[(min((jjj+1)*levmult,jmaxcalc-1))-jminsize][(iii*levmult)-iminsize];
                  if (ntt != nt && ntt > 0 && (ntt < start_idx || ntt >= end_idx) )
                     offtile_list.push_back(ntt);

                  int nrtt = hash[min((jjj+1)*levmult,jmaxcalc-1)-jminsize][min((iii*2+1)*levmultminus,imaxcalc-1)-iminsize];
                  if (nrtt != ntt && nrtt > 0 && (nrtt < start_idx || nrtt >= end_idx))
                     offtile_list.push_back(nrtt);
               }
            }

            if (lev < levmx) {
               int levmultminus = levtable[max(levmx-lev-1,0)];
               int nrt = hash[min((jj+1)*2*levmultminus,jmaxcalc-1)-jminsize][min((ii*2+1)*levmultminus,imaxcalc-1)-iminsize];
               if (nrt != nt && nrt > 0 && (nrt < start_idx || nrt >= end_idx) ) {
                  offtile_list.push_back(nrt);

                  int ig;
                  for (ig = 0; ig < nbsize_local; ig++) {
                     if (border_cell_num_global[ig] == nrt) break;
                  }
                  if (ig == nbsize_local) {
                     printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nrt);
#ifdef HAVE_MPI
                     //L7_Terminate();
#endif
                     //exit(-1);
                  }
                  int nlev = border_cell_level_global[ig];
                  int iii = border_cell_i_global[ig];
                  int jjj = border_cell_j_global[ig];
                  int levmult      = levtable[levmx-nlev];
                  levmultminus = levtable[max(levmx-nlev-1,0)];

                  // Multiplying by levmult gets the bottom left cell of the top neighbor. Add
                  // one to that to get the bottom left cell of the top top neighbor

                  if (nlev == levmx) {
                     int ntrt = hash[min((jjj+1)*levmult,jmaxcalc-1)-jminsize][iii*levmult-iminsize];
                     if (ntrt != nrt && ntrt > 0 && (ntrt < start_idx || ntrt >= end_idx) )
                        offtile_list.push_back(ntrt);
                  } else {
                     // And we need both cells here
                     int ntrt = hash[min((jjj+1)*2*levmultminus,jmaxcalc-1)-jminsize][(iii*2)*levmultminus-iminsize];
                     if (ntrt != nrt && ntrt > 0 && (ntrt < start_idx || ntrt >= end_idx) )
                        offtile_list.push_back(ntrt);
                     int nrtrt = hash[min((jjj+1)*2*levmultminus,jmaxcalc-1)-jminsize][min((iii*2+1)*levmultminus,imaxcalc-1)-iminsize];
                     if (nrtrt != ntrt && nrtrt > 0 && (nrtrt < start_idx || nrtrt >= end_idx))
                        offtile_list.push_back(nrtrt);
                  }
               }
            }
         } // ncells loop

         sort(offtile_list.begin(), offtile_list.end());
         p_end = unique(offtile_list.begin(), offtile_list.end());

         vector<int> indices_needed;

         for (vector<int>::iterator p=offtile_list.begin(); p < p_end; p++){
            indices_needed.push_back(*p);
         }

         if (TIMING_LEVEL >= 2) {
            cpu_time_offtile_list += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         int nghost = indices_needed.size();
         ncells_ghost = ncells + nghost;

         offtile_ratio_local = (offtile_ratio_local*(double)offtile_local_count) + ((double)nghost / (double)ncells);
         offtile_local_count++;
         offtile_ratio_local /= offtile_local_count;

         //printf("%d ncells size is %ld ncells_ghost size is %ld nghost %d\n",mype,ncells,ncells_ghost,nghost);
         //fprintf(fp,"%d ncells_ghost size is %ld nghost %d\n",mype,ncells_ghost,nghost);

#ifdef HAVE_MPI
         if (cell_handle) L7_Free(&cell_handle);
#endif
         cell_handle=0;
         //for (int ig = 0; ig<nghost; ig++){
         //   fprintf(fp,"%d: indices_needed[%d]=%d\n",mype,ig,indices_needed[ig]);
         //}
#ifdef HAVE_MPI
         L7_Setup(0, noffset, ncells, &indices_needed[0], nghost, &cell_handle);
#endif


         if (TIMING_LEVEL >= 2) {
            cpu_time_setup_comm += cpu_timer_stop(tstart_lev2);
            cpu_timer_start(&tstart_lev2);
         }

         celltype.resize(ncells_ghost);
         i.resize(ncells_ghost);
         j.resize(ncells_ghost);
         level.resize(ncells_ghost);
         nlft.resize(ncells_ghost,-98);
         nrht.resize(ncells_ghost,-98);
         nbot.resize(ncells_ghost,-98);
         ntop.resize(ncells_ghost,-98);

#ifdef HAVE_MPI
         L7_Update(&celltype[0], L7_INT, cell_handle);
         L7_Update(&i[0],        L7_INT, cell_handle);
         L7_Update(&j[0],        L7_INT, cell_handle);
         L7_Update(&level[0],    L7_INT, cell_handle);
#endif
         //for (int ic=0; ic<ncells; ic++){
         //   fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
         //       mype,ic,i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
         //}
#ifdef HAVE_MPI
         L7_Update(&nlft[0], L7_INT, cell_handle);
         L7_Update(&nrht[0], L7_INT, cell_handle);
         L7_Update(&nbot[0], L7_INT, cell_handle);
         L7_Update(&ntop[0], L7_INT, cell_handle);
#endif
         //for (int ic=0; ic<ncells_ghost; ic++){
         //   fprintf(fp,"%d: 1655 nlft for %5d is %5d\n",mype,ic,nlft[ic]);
         //}

         for (uint ic=0; ic<ncells_ghost; ic++){
            if (nlft[ic] >= noffset && nlft[ic] < (int)(noffset+ncells)) {
               nlft[ic] -= noffset;
               //fprintf(fp,"%d: 1: ic %d nlft is %d\n",mype,ic,nlft[ic]);
            } else {
               for (int ig=0; ig<nghost; ig++){
                  //if (nlft[ic]==indices_needed[ig]) {nlft[ic] = ig+ncells; fprintf(fp,"%d: 2: ic %d nlft is %d\n",mype,ic,nlft[ic]); break;}
                  if (nlft[ic]==indices_needed[ig]) {nlft[ic] = ig+ncells; break;}
               }
            }
            if (nrht[ic] >= noffset && nrht[ic] < (int)(noffset+ncells)) {
               nrht[ic] -= noffset;
            } else {
               for (int ig=0; ig<nghost; ig++){
                  if (nrht[ic]==indices_needed[ig]) {nrht[ic] = ig+ncells; break;}
               }
            }
            if (nbot[ic] >= noffset && nbot[ic] < (int)(noffset+ncells)) {
               nbot[ic] -= noffset;
            } else {
               for (int ig=0; ig<nghost; ig++){
                  if (nbot[ic]==indices_needed[ig]) {nbot[ic] = ig+ncells; break;}
               }
            }
            if (ntop[ic] >= noffset && ntop[ic] < (int)(noffset+ncells)) {
               ntop[ic] -= noffset;
            } else {
               for (int ig=0; ig<nghost; ig++){
                  if (ntop[ic]==indices_needed[ig]) {ntop[ic] = ig+ncells; break;}
               }
            }
         }

         if (DEBUG) {
            int jmaxglobal = (jmax+1)*levtable[levmx];
            int imaxglobal = (imax+1)*levtable[levmx];
            fprintf(fp,"\n                                    HASH 1 numbering\n");
            for (int jj = jmaxglobal-1; jj>=0; jj--){
               fprintf(fp,"%2d: %4d:",mype,jj);
               if (jj >= jminsize && jj < jmaxsize) {
                  for (int ii = 0; ii<imaxglobal; ii++){
                     if (ii >= iminsize && ii < imaxsize) {
                        fprintf(fp,"%5d",hash[jj-jminsize][ii-iminsize]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                           fprintf(fp,"%5d",nlft[hashval]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",nrht[hashval]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",nbot[hashval]);
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
                     int hashval = hash[jj-jminsize][ii-iminsize]-noffset;
                     if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < (int)ncells) ) {
                        fprintf(fp,"%5d",ntop[hashval]);
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
            for (uint ic=0; ic<ncells; ic++){
               fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
                   mype,ic,i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            }
            int ig=0;
            for (uint ic=ncells; ic<ncells_ghost; ic++, ig++){
               fprintf(fp,"%d: after  update ic %d off %d i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
                   mype,ic,indices_needed[ig],i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
            }
         }
      }

      genmatrixfree((void **)hash);

      if (TIMING_LEVEL >= 2) cpu_time_do_ghost_comm += cpu_timer_stop(tstart_lev2);

   } else if (calc_neighbor_type == KDTREE) {
      struct timeval tstart_lev2;
      if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

      TBounds box;
      vector<int> index_list((int)pow(2,levmx*levmx) );

      int num;

      ibase = 0;
      calc_spatial_coordinates(ibase);

      kdtree_setup();

      if (TIMING_LEVEL >= 2) {
         cpu_time_kdtree_setup += cpu_timer_stop(tstart_lev2);
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

      if (TIMING_LEVEL >= 2) cpu_time_kdtree_query += cpu_timer_stop(tstart_lev2);

   } // calc_neighbor_type

   cpu_time_calc_neighbors += cpu_timer_stop(tstart_cpu);
}

#ifdef HAVE_OPENCL
void Mesh::gpu_calc_neighbors(cl_command_queue command_queue)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timeval tstart_lev2;
   cpu_timer_start(&tstart_lev2);

   gpu_calc_neigh_counter++;

   assert(dev_levtable);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);

   dev_nlft     = ezcl_malloc(NULL, const_cast<char *>("dev_nlft"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nrht     = ezcl_malloc(NULL, const_cast<char *>("dev_nrht"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nbot     = ezcl_malloc(NULL, const_cast<char *>("dev_nbot"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_ntop     = ezcl_malloc(NULL, const_cast<char *>("dev_ntop"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   size_t local_work_size = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;

   int imaxsize = (imax+1)*levtable[levmx];
   int jmaxsize = (jmax+1)*levtable[levmx];

   size_t hashsize = jmaxsize*imaxsize;
   cl_mem dev_hash = ezcl_malloc(NULL, const_cast<char *>("dev_hash"), &hashsize, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   size_t hash_local_work_size  = MIN(hashsize, TILE_SIZE);
   size_t hash_global_work_size = ((hashsize+hash_local_work_size - 1) /hash_local_work_size) * hash_local_work_size;

   ezcl_set_kernel_arg(kernel_hash_init, 0, sizeof(cl_ulong), (void *)&hashsize);
   ezcl_set_kernel_arg(kernel_hash_init, 1, sizeof(cl_mem),   (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init,   1, NULL, &hash_global_work_size, &hash_local_work_size, NULL);

      /*
                    const int  isize,     // 0
                    const int  levmx,     // 1
                    const int  imax,      // 2
                    const int  jmax,      // 3
                    const int  imaxsize,  // 4
           __global const int  *levtable, // 5
           __global const int  *lev_ibeg, // 6
           __global const int  *lev_iend, // 7
           __global const int  *lev_jbeg, // 8
           __global const int  *lev_jend, // 9
           __global       int  *level,    // 10
           __global       int  *i,        // 11
           __global       int  *j,        // 12
           __global       int  *hash)     // 13
      */

   ezcl_set_kernel_arg(kernel_hash_setup,  0, sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_hash_setup,  1, sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_setup,  2, sizeof(cl_int), (void *)&imax);
   ezcl_set_kernel_arg(kernel_hash_setup,  3, sizeof(cl_int), (void *)&jmax);
   ezcl_set_kernel_arg(kernel_hash_setup,  4, sizeof(cl_int), (void *)&imaxsize);
   ezcl_set_kernel_arg(kernel_hash_setup,  5, sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_setup,  6, sizeof(cl_mem), (void *)&dev_levibeg);
   ezcl_set_kernel_arg(kernel_hash_setup,  7, sizeof(cl_mem), (void *)&dev_leviend);
   ezcl_set_kernel_arg(kernel_hash_setup,  8, sizeof(cl_mem), (void *)&dev_levjbeg);
   ezcl_set_kernel_arg(kernel_hash_setup,  9, sizeof(cl_mem), (void *)&dev_levjend);
   ezcl_set_kernel_arg(kernel_hash_setup, 10, sizeof(cl_mem), (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_hash_setup, 11, sizeof(cl_mem), (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_hash_setup, 12, sizeof(cl_mem), (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_hash_setup, 13, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup,   1, NULL, &global_work_size, &local_work_size, NULL);

      /*
                    const int  isize,     // 0
                    const int  levmx,     // 1
                    const int  imax,      // 2
                    const int  jmax,      // 3
                    const int  imaxsize,  // 4
                    const int  jmaxsize,  // 5
           __global       int  *levtable, // 6
           __global       int  *level,    // 7
           __global       int  *i,        // 8
           __global       int  *j,        // 9
           __global       int  *nlft,     // 10
           __global       int  *nrht,     // 11
           __global       int  *nbot,     // 12
           __global       int  *ntop,     // 13
           __global       int  *hash)     // 14
      */

   if (TIMING_LEVEL >= 2) {
      ezcl_finish(command_queue);
      gpu_time_hash_setup += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      cpu_timer_start(&tstart_lev2);
   }

   ezcl_set_kernel_arg(kernel_calc_neighbors, 0,  sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 1,  sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 2,  sizeof(cl_int), (void *)&imax);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 3,  sizeof(cl_int), (void *)&jmax);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 4,  sizeof(cl_int), (void *)&imaxsize);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 5,  sizeof(cl_int), (void *)&jmaxsize);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 6,  sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 7,  sizeof(cl_mem), (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 8,  sizeof(cl_mem), (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 9,  sizeof(cl_mem), (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 10, sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 11, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 12, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 13, sizeof(cl_mem), (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_neighbors, 14, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors,   1, NULL, &global_work_size, &local_work_size, NULL);

   ezcl_device_memory_remove(dev_hash);

   ezcl_finish(command_queue);

   if (TIMING_LEVEL >= 2) gpu_time_hash_query += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);

   gpu_time_calc_neighbors += (long)(cpu_timer_stop(tstart_cpu) * 1.0e9);
}


void Mesh::gpu_calc_neighbors_local(cl_command_queue command_queue)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timeval tstart_lev2;
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

   gpu_calc_neigh_counter++;

   ncells_ghost = ncells;

   assert(dev_levtable);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);

   dev_nlft     = ezcl_malloc(NULL, const_cast<char *>("dev_nlft"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nrht     = ezcl_malloc(NULL, const_cast<char *>("dev_nrht"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nbot     = ezcl_malloc(NULL, const_cast<char *>("dev_nbot"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_ntop     = ezcl_malloc(NULL, const_cast<char *>("dev_ntop"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   size_t local_work_size =  64;
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   //printf("DEBUG file %s line %d lws = %d gws %d bs %d ncells %d\n",__FILE__,__LINE__,
   //   local_work_size, global_work_size, block_size, ncells);
   cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_int4), CL_MEM_READ_WRITE, 0);
   size_t one = 1;
   cl_mem dev_sizes = ezcl_malloc(NULL, const_cast<char *>("dev_sizes"), &one, sizeof(cl_int4),  CL_MEM_READ_WRITE, 0);

      /*
       __kernel void calc_hash_size_cl(
                          const int   ncells,      // 0
                          const int   levmx,       // 1
                 __global       int   *levtable,   // 2
                 __global       int   *level,      // 3
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

   ezcl_device_memory_remove(dev_redscratch);

   cl_int sizes[4];
   ezcl_enqueue_read_buffer(command_queue, dev_sizes, CL_TRUE,  0, 1*sizeof(cl_int4), &sizes, NULL);

   int iminsize = sizes[0];
   int imaxsize = sizes[1];
   int jminsize = sizes[2];
   int jmaxsize = sizes[3];

   //fprintf(fp,"%d: sizes %d %d     %d %d   \n",mype,sizes[0],sizes[1],sizes[2],sizes[3]);
   //fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);

   // Expand size by 2*coarse_cells for ghost cells
   jminsize = max(jminsize-2*levtable[levmx],0);
   jmaxsize = min(jmaxsize+2*levtable[levmx],(jmax+1)*levtable[levmx]);
   iminsize = max(iminsize-2*levtable[levmx],0);
   imaxsize = min(imaxsize+2*levtable[levmx],(imax+1)*levtable[levmx]);
   //fprintf(fp,"%d: Sizes are imin %d imax %d jmin %d jmax %d\n",mype,iminsize,imaxsize,jminsize,jmaxsize);

   // Allocate partial hash table
   //int **hash = (int **)genmatrix(jmaxsize-jminsize,imaxsize-iminsize,sizeof(int));

   //int imaxsize = (imax+1)*levtable[levmx];
   //int jmaxsize = (jmax+1)*levtable[levmx];

   size_t hashsize = (jmaxsize-jminsize)*(imaxsize-iminsize);
   cl_mem dev_hash = ezcl_malloc(NULL, const_cast<char *>("dev_hash"), &hashsize, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   size_t hash_local_work_size  = MIN(hashsize, TILE_SIZE);
   size_t hash_global_work_size = ((hashsize+hash_local_work_size - 1) /hash_local_work_size) * hash_local_work_size;

   //printf("%d: hash size is %d lws %d gws %d\n",mype,hashsize,hash_local_work_size,hash_global_work_size);
   ezcl_set_kernel_arg(kernel_hash_init, 0, sizeof(cl_ulong), (void *)&hashsize);
   ezcl_set_kernel_arg(kernel_hash_init, 1, sizeof(cl_mem),   (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init,   1, NULL, &hash_global_work_size, &hash_local_work_size, NULL);

   int csize = corners_i.size();
   //for (int ic=0; ic<csize; ic++){
   //   printf("%d: corners i %d j %d hash %d\n",mype,corners_i[ic],corners_j[ic],
   //      (corners_j[ic]-jminsize)*(imaxsize-iminsize)+(corners_i[ic]-iminsize));
   //}
   size_t corners_local_work_size  = MIN(csize, TILE_SIZE);
   size_t corners_global_work_size = ((csize+corners_local_work_size - 1) /corners_local_work_size) * corners_local_work_size;

   ezcl_set_kernel_arg(kernel_hash_init_corners, 0, sizeof(cl_int), (void *)&csize);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 1, sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 2, sizeof(cl_int), (void *)&imax);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 3, sizeof(cl_int), (void *)&jmax);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 4, sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 5, sizeof(cl_mem), (void *)&dev_corners_i);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 6, sizeof(cl_mem), (void *)&dev_corners_j);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 7, sizeof(cl_mem), (void *)&dev_sizes);
   ezcl_set_kernel_arg(kernel_hash_init_corners, 8, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init_corners,   1, NULL, &corners_global_work_size, &corners_local_work_size, NULL);

      /*
                    const int  isize,     // 0
                    const int  levmx,     // 1
           __global       int  *levtable, // 2
           __global       int  *level,    // 3
           __global       int  *i,        // 4
           __global       int  *j,        // 5
           __global       int  *hash)     // 6
      */

   local_work_size = 128;
   global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;

   //printf("%d: lws %d gws %d \n",mype,local_work_size,global_work_size);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  0, sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  1, sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  2, sizeof(cl_int), (void *)&imax);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  3, sizeof(cl_int), (void *)&jmax);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  4, sizeof(cl_int), (void *)&noffset);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  5, sizeof(cl_mem), (void *)&dev_sizes);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  6, sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  7, sizeof(cl_mem), (void *)&dev_levibeg);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  8, sizeof(cl_mem), (void *)&dev_leviend);
   ezcl_set_kernel_arg(kernel_hash_setup_local,  9, sizeof(cl_mem), (void *)&dev_levjbeg);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 10, sizeof(cl_mem), (void *)&dev_levjend);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 11, sizeof(cl_mem), (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 12, sizeof(cl_mem), (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 13, sizeof(cl_mem), (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_hash_setup_local, 14, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup_local,   1, NULL, &global_work_size, &local_work_size, NULL);

   if (TIMING_LEVEL >= 2) {
      ezcl_finish(command_queue);
      gpu_time_hash_setup += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      cpu_timer_start(&tstart_lev2);
   }

      /*
                    const int  isize,     // 0
                    const int  levmx,     // 1
                    const int  imaxsize,  // 2
                    const int  jmaxsize,  // 3
           __global       int  *levtable, // 4
           __global       int  *level,    // 5
           __global       int  *i,        // 6
           __global       int  *j,        // 7
           __global       int  *nlft,     // 8
           __global       int  *nrht,     // 9
           __global       int  *nbot,     // 10
           __global       int  *ntop,     // 11
           __global       int  *hash)     // 12
      */

   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 0,  sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 1,  sizeof(cl_int), (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 2,  sizeof(cl_int), (void *)&imax);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 3,  sizeof(cl_int), (void *)&jmax);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 4,  sizeof(cl_mem), (void *)&dev_sizes);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 5,  sizeof(cl_mem), (void *)&dev_levtable);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 6,  sizeof(cl_mem), (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 7,  sizeof(cl_mem), (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 8,  sizeof(cl_mem), (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 9,  sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 10, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 11, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 12, sizeof(cl_mem), (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_neighbors_local, 13, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors_local,   1, NULL, &global_work_size, &local_work_size, NULL);

   if (TIMING_LEVEL >= 2) {
      ezcl_finish(command_queue);
      gpu_time_hash_query += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
      cpu_timer_start(&tstart_lev2);
   }

   if (numpe > 1) {
      vector<int> border_cell(ncells);

      cl_mem dev_border_cell = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_cell2 = ezcl_malloc(NULL, const_cast<char *>("dev_border_cell2"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      ezcl_set_kernel_arg(kernel_calc_border_cells, 0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 1,  sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 2,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 3,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 4,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 5,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 6,  sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_calc_border_cells, 7,  sizeof(cl_mem), (void *)&dev_border_cell);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_border_cells, 1, NULL, &global_work_size, &local_work_size, NULL); 

      ezcl_set_kernel_arg(kernel_calc_border_cells2, 0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 1,  sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 2,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 3,  sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 4,  sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 5,  sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 6,  sizeof(cl_mem), (void *)&dev_border_cell);
      ezcl_set_kernel_arg(kernel_calc_border_cells2, 7,  sizeof(cl_mem), (void *)&dev_border_cell2);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_border_cells2, 1, NULL, &global_work_size, &local_work_size, NULL); 

      ezcl_enqueue_read_buffer(command_queue, dev_border_cell2, CL_TRUE,  0, ncells*sizeof(cl_int), &border_cell[0], NULL);

      ezcl_device_memory_remove(dev_border_cell);
      ezcl_device_memory_remove(dev_border_cell2);

      vector<int> nlft_tmp(ncells);
      vector<int> nrht_tmp(ncells);
      vector<int> nbot_tmp(ncells);
      vector<int> ntop_tmp(ncells);

      vector<int> celltype_tmp(ncells);
      vector<int> i_tmp(ncells);
      vector<int> j_tmp(ncells);
      vector<int> level_tmp(ncells);

      ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_FALSE, 0, ncells*sizeof(cl_int), &celltype_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells*sizeof(cl_int), &i_tmp[0],        NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells*sizeof(cl_int), &j_tmp[0],        NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_TRUE,  0, ncells*sizeof(cl_int), &level_tmp[0],    NULL);

      //ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells*sizeof(cl_int), &nlft_tmp[0], NULL);
      //ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells*sizeof(cl_int), &nrht_tmp[0], NULL);
      //ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells*sizeof(cl_int), &nbot_tmp[0], NULL);
      //ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells*sizeof(cl_int), &ntop_tmp[0], NULL);

#ifdef XXX
      vector<int> border_cell(ncells,-1);

      // Scan for corner boundary cells and also push list of unsatisfied neighbor cells
      for (uint ic=0; ic<ncells; ic++){
         if (nlft_tmp[ic] == -1 || (level_tmp[nlft_tmp[ic]-noffset] > level_tmp[ic] && ntop_tmp[nlft_tmp[ic]-noffset] == -1) ){
            //printf("%d: Cell is %d nlft %d\n",mype,ic+noffset,nlft_tmp[ic]);
            border_cell[ic]=1;
            if (nrht_tmp[ic] >= 0) {
               border_cell[nrht_tmp[ic]-noffset]=11;
               if (level_tmp[nrht_tmp[ic]-noffset] > level_tmp[ic]) {
                  if (ntop_tmp[nrht_tmp[ic]-noffset] >= 0) border_cell[ntop_tmp[nrht_tmp[ic]-noffset]-noffset]=12;
               }
            }
         }
         if (nrht_tmp[ic] == -1 || (level_tmp[nrht_tmp[ic]-noffset] > level_tmp[ic] && ntop_tmp[nrht_tmp[ic]-noffset] == -1) ){
            //printf("%d: Cell is %d nrht %d\n",mype,ic+noffset,nrht_tmp[ic]);
            border_cell[ic]=2;
            if (nlft_tmp[ic] >= 0) {
               border_cell[nlft_tmp[ic]-noffset]=21;
               if (level_tmp[nlft_tmp[ic]-noffset] > level_tmp[ic]) {
                  if (ntop_tmp[nlft_tmp[ic]-noffset] >= 0) border_cell[ntop_tmp[nlft_tmp[ic]-noffset]-noffset]=22;
               }
            }
         }
         if (nbot_tmp[ic] == -1 || (level_tmp[nbot_tmp[ic]-noffset] > level_tmp[ic] && nrht_tmp[nbot_tmp[ic]-noffset] == -1) ) {
            //printf("%d: Cell is %d nbot %d\n",mype,ic+noffset,nbot_tmp[ic]);
            border_cell[ic]=3;
            if (ntop_tmp[ic] >= 0) {
               border_cell[ntop_tmp[ic]-noffset]=31;
               if (level_tmp[ntop_tmp[ic]-noffset] > level_tmp[ic]) {
                  if (nrht_tmp[ntop_tmp[ic]-noffset] >= 0) border_cell[nrht_tmp[ntop_tmp[ic]-noffset]-noffset]=32;
               }
            }
         }
         if (ntop_tmp[ic] == -1 || ( level_tmp[ntop_tmp[ic]-noffset] > level_tmp[ic] && nrht_tmp[ntop_tmp[ic]-noffset] == -1) ) {
            //printf("%d: Cell is %d ntop %d\n",mype,ic+noffset,ntop_tmp[ic]);
            border_cell[ic]=4;
            if (nbot_tmp[ic] >= 0) {
               border_cell[nbot_tmp[ic]-noffset]=41;
               if (level_tmp[nbot_tmp[ic]-noffset] > level_tmp[ic]) {
                  if (nrht_tmp[nbot_tmp[ic]-noffset] >= 0) border_cell[nrht_tmp[nbot_tmp[ic]-noffset]-noffset]=42;
               }
            }
         }
      }

      int icount = 0;
      for(int ic = 0; ic < (int)ncells; ic++){
         if ((border_cell[ic] > 0 && border_cell_tmp[ic] < 0) || (border_cell[ic] < 0 && border_cell_tmp[ic] > 0) ) {
            printf("%d: DEBUG ic %d border_cell %d border_cell_tmp %d\n",mype,ic,border_cell[ic],border_cell_tmp[ic]);
            //printf("%d: ic %d i %d j %d level %d nlft %d nrht %d nbot %d ntop %d\n",mype,ic,i_tmp[ic],j_tmp[ic],level_tmp[ic], nlft_tmp[ic], nrht_tmp[ic], nbot_tmp[ic], ntop_tmp[ic]);
            //printf("%d: ic %d border nlft %d border nrht %d border nbot %d border ntop %d\n",mype,ic, border_cell_tmp[nlft_tmp[ic]-noffset], border_cell_tmp[nrht_tmp[ic]-noffset], border_cell_tmp[nbot_tmp[ic]-noffset], border_cell_tmp[ntop_tmp[ic]-noffset]);
            icount++;
         }
         if (icount > 20) {
            //sleep(20);
#ifdef HAVE_MPI
            L7_Terminate();
#endif
            exit(0);
         }
      }
      //if (icount>0) sleep(20);
#endif

      vector<int> border_cell_num;

      for(int ic = 0; ic < (int)ncells; ic++){
         if (border_cell[ic] > 0) border_cell_num.push_back(ic+noffset);
      }

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_find_boundary += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      int nbsize_global;
      int nbsize_local=border_cell_num.size();

      nbsize_global = nbsize_local;
#ifdef HAVE_MPI
      MPI_Allreduce(&nbsize_local, &nbsize_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

      //printf("%d: border cell size is %d global is %d\n",mype,nbsize_local,nbsize_global);

      vector<int> border_cell_i(nbsize_local);
      vector<int> border_cell_j(nbsize_local);
      vector<int> border_cell_level(nbsize_local);
    
      for (int ic = 0; ic <nbsize_local; ic++){
         int cell_num = border_cell_num[ic]-noffset;
         border_cell_i[ic] = i_tmp[cell_num]; 
         border_cell_j[ic] = j_tmp[cell_num]; 
         border_cell_level[ic] = level_tmp[cell_num]; 
         //fprintf(fp,"%d: Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num[ic],
         //   border_cell_i[ic],border_cell_j[ic],border_cell_level[ic]);
      }

      vector<int> nbsizes(numpe);
      vector<int> nbdispl(numpe);
#ifdef HAVE_MPI
      MPI_Allgather(&nbsize_local, 1, MPI_INT, &nbsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
#endif
      nbdispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         nbdispl[ip] += nbdispl[ip-1] + nbsizes[ip-1];
      }

      vector<int>border_cell_num_global(nbsize_global);
      vector<int>border_cell_i_global(nbsize_global);
      vector<int>border_cell_j_global(nbsize_global);
      vector<int>border_cell_level_global(nbsize_global);

#ifdef HAVE_MPI
      MPI_Allgatherv(&border_cell_num[0],   nbsizes[mype], MPI_INT, &border_cell_num_global[0],   &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&border_cell_i[0],     nbsizes[mype], MPI_INT, &border_cell_i_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&border_cell_j[0],     nbsizes[mype], MPI_INT, &border_cell_j_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&border_cell_level[0], nbsizes[mype], MPI_INT, &border_cell_level_global[0], &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
#endif

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_gather_boundary += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      //for (int ic = 0; ic < nbsize_global; ic++) {
      //   fprintf(fp,"%d: Global Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num_global[ic],
      //      border_cell_i_global[ic],border_cell_j_global[ic],border_cell_level_global[ic]);
      //}

      int inew=0;
      for (int ic = 0; ic < nbsize_global; ic++) {
         int lev = border_cell_level_global[ic];
         int levmult = levtable[levmx-lev];
         //fprintf(fp,"%d: DEBUG cell %d i %d j %d\n",mype,ic,border_cell_i_global[ic],border_cell_j_global[ic]);
         if (border_cell_j_global[ic]*levmult < jminsize || border_cell_j_global[ic]*levmult >= jmaxsize) continue;
         if (border_cell_i_global[ic]*levmult < iminsize || border_cell_i_global[ic]*levmult >= imaxsize) continue;
         if (border_cell_num_global[ic] >= noffset && border_cell_num_global[ic] < (int)(noffset+ncells)) continue;
         //fprintf(fp,"%d: ic is %d inew is %d\n",mype,ic,inew);
         if (inew != ic){
            border_cell_num_global[inew] = border_cell_num_global[ic];
            border_cell_i_global[inew] = border_cell_i_global[ic];
            border_cell_j_global[inew] = border_cell_j_global[ic];
            border_cell_level_global[inew] = border_cell_level_global[ic];
         }
         inew++;
      }
      nbsize_local = inew;

      //fprintf(fp,"%d: nbsize_local is %d\n",mype,nbsize_local);
      //for (int ic = 0; ic < nbsize_local; ic++) {
      //   fprintf(fp,"%d: Local Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_num_global[ic],
      //      border_cell_i_global[ic],border_cell_j_global[ic],border_cell_level_global[ic]);
      //}

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_hash_setup2 += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      size_t nbsize = nbsize_local;
      cl_mem dev_border_level = ezcl_malloc(NULL, const_cast<char *>("dev_border_level"), &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_i     = ezcl_malloc(NULL, const_cast<char *>("dev_border_i"),     &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_j     = ezcl_malloc(NULL, const_cast<char *>("dev_border_j"),     &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_border_num   = ezcl_malloc(NULL, const_cast<char *>("dev_border_num"),   &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

      ezcl_enqueue_write_buffer(command_queue, dev_border_level, CL_FALSE,  0, nbsize*sizeof(cl_int), &border_cell_level_global[0], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_border_i    , CL_FALSE,  0, nbsize*sizeof(cl_int), &border_cell_i_global[0],     NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_border_j    , CL_FALSE,  0, nbsize*sizeof(cl_int), &border_cell_j_global[0],     NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_border_num  , CL_TRUE,   0, nbsize*sizeof(cl_int), &border_cell_num_global[0],   NULL);

      //size_t border_local_work_size = MIN(nbsize,TILE_SIZE);
      size_t border_local_work_size = 32;
      size_t border_global_work_size = ((nbsize + border_local_work_size - 1) /border_local_work_size) * border_local_work_size;
      //printf("%d: blws %d bgws %d\n",mype,border_local_work_size,border_global_work_size);

      ezcl_set_kernel_arg(kernel_hash_setup_border, 0, sizeof(cl_int), (void *)&nbsize_local);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 1, sizeof(cl_int), (void *)&levmx);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 2, sizeof(cl_int), (void *)&noffset);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 3, sizeof(cl_mem), (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 4, sizeof(cl_mem), (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 5, sizeof(cl_mem), (void *)&dev_border_level);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 6, sizeof(cl_mem), (void *)&dev_border_i);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 7, sizeof(cl_mem), (void *)&dev_border_j);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 8, sizeof(cl_mem), (void *)&dev_border_num);
      ezcl_set_kernel_arg(kernel_hash_setup_border, 9, sizeof(cl_mem), (void *)&dev_hash);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup_border,  1, NULL, &border_global_work_size, &border_local_work_size, NULL);

      ezcl_device_memory_remove(dev_border_level);
      ezcl_device_memory_remove(dev_border_i);
      ezcl_device_memory_remove(dev_border_j);
      ezcl_device_memory_remove(dev_border_num);

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_hash_setup2 += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 0,  sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 1,  sizeof(cl_int), (void *)&levmx);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 2,  sizeof(cl_int), (void *)&imax);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 3,  sizeof(cl_int), (void *)&jmax);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 4,  sizeof(cl_mem), (void *)&dev_sizes);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 5,  sizeof(cl_mem), (void *)&dev_levtable);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 6,  sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 7,  sizeof(cl_mem), (void *)&dev_i);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 8,  sizeof(cl_mem), (void *)&dev_j);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 9,  sizeof(cl_mem), (void *)&dev_nlft);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 10, sizeof(cl_mem), (void *)&dev_nrht);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 11, sizeof(cl_mem), (void *)&dev_nbot);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 12, sizeof(cl_mem), (void *)&dev_ntop);
      ezcl_set_kernel_arg(kernel_calc_neighbors_local2, 13, sizeof(cl_mem), (void *)&dev_hash);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors_local2,   1, NULL, &global_work_size, &local_work_size, NULL);

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_hash_query2 += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      vector<int> hash_tmp(hashsize);
      ezcl_enqueue_read_buffer(command_queue, dev_hash, CL_TRUE,  0, hashsize*sizeof(cl_int), &hash_tmp[0], NULL);

      ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells*sizeof(cl_int), &nlft_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells*sizeof(cl_int), &nrht_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells*sizeof(cl_int), &nbot_tmp[0], NULL);
      ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells*sizeof(cl_int), &ntop_tmp[0], NULL);

      if (DEBUG) {

         int jmaxglobal = (jmax+1)*levtable[levmx];
         int imaxglobal = (imax+1)*levtable[levmx];
         fprintf(fp,"\n                                    HASH numbering\n");
         for (int jj = jmaxglobal-1; jj>=0; jj--){
            fprintf(fp,"%2d: %4d:",mype,jj);
            if (jj >= jminsize && jj < jmaxsize) {
               for (int ii = 0; ii<imaxglobal; ii++){
                  if (ii >= iminsize && ii < imaxsize) {
                     fprintf(fp,"%5d",hash_tmp[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)]);
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
                  int hashval = hash_tmp[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)]-noffset;
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
                  int hashval = hash_tmp[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)]-noffset;
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
                  int hashval = hash_tmp[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)]-noffset;
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
                  int hashval = hash_tmp[(jj-jminsize)*(imaxsize-iminsize)+(ii-iminsize)]-noffset;
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
      
      vector<int> offtile_list;

      int start_idx = noffset;
      int end_idx   = noffset+ncells;
      int ii, jj, lev;

      int jmaxcalc = (jmax+1)*levtable[levmx];
      int imaxcalc = (imax+1)*levtable[levmx];

      for (uint ic = 0; ic < ncells; ic++){
         ii = i_tmp[ic];
         jj = j_tmp[ic];
         lev = level_tmp[ic];

         int nl = nlft_tmp[ic];

         if (nl < start_idx || nl >= end_idx) {
            offtile_list.push_back(nl);
            int ig;
            for (ig = 0; ig < nbsize_local; ig++) {
               if (border_cell_num_global[ig] == nl) break;
            }
            if (ig == nbsize_local) {
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nl);
#ifdef HAVE_MPI
               L7_Terminate();
#endif
               exit(-1);
            }
            // Get the left neighbor information
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[max(levmx-nlev-1,0)];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];
            //fprintf(fp,"%d: DEBUG nl is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nl, ig, border_cell_num_global[ig], nlev, iii, jjj);

            if (nlev == levmx) {
               // Multiplying by levmult get the bottom left cell of the left neighbor. Subtract
               // one from that to get the bottom right cell in the left left neighbor
               int nll = hash_tmp[((jjj*levmult)-jminsize)*(imaxsize-iminsize)+(max((iii-1)*levmult, 0)-iminsize)];
               if (nll != nl && nll > 0 && (nll < start_idx || nll >= end_idx) ) 
                  offtile_list.push_back(nll);
            } else {
               int nll = hash_tmp[(jjj*levmult-jminsize)*(imaxsize-iminsize)+(max((iii*2-1)*levmultminus, 0)-iminsize)];
               if (nll != nl && nll > 0 && (nll < start_idx || nll >= end_idx) ) 
                  offtile_list.push_back(nll);

               int ntll = hash_tmp[((jjj*2+1)*levmultminus-jminsize)*(imaxsize-iminsize)+(max((iii*2-1)*levmultminus, 0)-iminsize)];
               if (ntll != nll && ntll > 0 && (ntll < start_idx || ntll >= end_idx) )
                  offtile_list.push_back(ntll);
            }
         }

         if (lev < levmx) {
            int levmultminus = levtable[max(levmx-lev-1,0)];
            int ntl = hash_tmp[(min((jj*2+1)*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(max((ii*2-1)*levmultminus,0)-iminsize)];
            if (ntl != nl && ntl > 0 && (ntl < start_idx || ntl >= end_idx) ) {
               offtile_list.push_back(ntl);

               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == ntl) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, ntl);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               int nlev = border_cell_level_global[ig];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];
               int levmult      = levtable[levmx-nlev];
               levmultminus = levtable[max(levmx-nlev-1,0)];

               if (nlev == levmx) {
                  int nltl = hash_tmp[(jjj*levmult-jminsize)*(imaxsize-iminsize)+(max((iii-1)*levmult-iminsize,0))];
                  if (nltl != ntl && nltl > 0 && (nltl < start_idx || nltl >= end_idx) )
                     offtile_list.push_back(nltl);
               } else {
                  // And we need both cells here
                  int nltl = hash_tmp[((jjj*2)*levmultminus-jminsize)*(imaxsize-iminsize)+(max((iii*2-1)*levmultminus-iminsize,0))];
                  if (nltl != ntl && nltl > 0 && (nltl < start_idx || nltl >= end_idx) )
                     offtile_list.push_back(nltl);
                  int ntltl = hash_tmp[(min((jjj*2+1)*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(max((iii*2-1)*levmultminus,0)-iminsize)];
                  if (ntltl != nltl && ntltl > 0 && (ntltl < start_idx || ntltl >= end_idx))
                     offtile_list.push_back(ntltl);
               }

            }
         }

         int nr = nrht_tmp[ic];

         if (nr < start_idx || nr >= end_idx) {
            offtile_list.push_back(nr);
            int ig;
            for (ig = 0; ig < nbsize_local; ig++) {
               if (border_cell_num_global[ig] == nr) break;
            }
            if (ig == nbsize_local) {
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nr);
#ifdef HAVE_MPI
               //L7_Terminate();
#endif
               //exit(-1);
            }
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[max(levmx-nlev-1,0)];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];

            if (nlev == levmx) {
               // Multiplying by levmult gets the bottom left cell of the right neighbor. Add
               // one to that to get the bottom left cell of the right right neighbor
               int nrr = hash_tmp[((jjj*levmult)-jminsize)*(imaxsize-iminsize)+(min((iii+1)*levmult,imaxcalc-1)-iminsize)];
               if (nrr != nr && nrr > 0 && (nrr < start_idx || nrr >= end_idx) )
                  offtile_list.push_back(nrr);
            } else {
               int nrr = hash_tmp[((jjj*levmult)-jminsize)*(imaxsize-iminsize)+((min((iii+1)*levmult,imaxcalc-1))-iminsize)];
               if (nrr != nr && nrr > 0 && (nrr < start_idx || nrr >= end_idx) )
                  offtile_list.push_back(nrr);

               int ntrr = hash_tmp[(min((jjj*2+1)*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min((iii+1)*levmult,imaxcalc-1)-iminsize)];
               if (ntrr != nrr && ntrr > 0 && (ntrr < start_idx || ntrr >= end_idx))
                  offtile_list.push_back(ntrr);
            }

         }

         if (lev < levmx) {
            int levmultminus = levtable[max(levmx-lev-1,0)];
            int ntr = hash_tmp[(min((jj*2+1)*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min((ii+1)*2*levmultminus,imaxcalc-1)-iminsize)];
            if (ntr != nr && ntr > 0 && (ntr < start_idx || ntr >= end_idx) ) {
               offtile_list.push_back(ntr);

               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == ntr) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, ntr);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               int nlev = border_cell_level_global[ig];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];
               int levmult      = levtable[levmx-nlev];
               levmultminus = levtable[max(levmx-nlev-1,0)];

               // Multiplying by levmult gets the bottom left cell of the top neighbor. Add
               // one to that to get the bottom left cell of the top top neighbor

               if (nlev == levmx) {
                  int nrtr = hash_tmp[(jjj*levmult-jminsize)*(imaxsize-iminsize)+(min((iii+1)*levmult-iminsize,imaxcalc-1))];
                  if (nrtr != ntr && nrtr > 0 && (nrtr < start_idx || nrtr >= end_idx) )
                     offtile_list.push_back(nrtr);
               } else {
                  // And we need both cells here
                  int nrtr = hash_tmp[((jjj*2)*levmultminus-jminsize)*(imaxsize-iminsize)+(min((iii+1)*2*levmultminus-iminsize,imaxcalc-1))];
                  if (nrtr != ntr && nrtr > 0 && (nrtr < start_idx || nrtr >= end_idx) )
                     offtile_list.push_back(nrtr);
                  int ntrtr = hash_tmp[(min((jjj*2+1)*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min((iii+1)*2*levmultminus,imaxcalc-1)-iminsize)];
                  if (ntrtr != nrtr && ntrtr > 0 && (ntrtr < start_idx || ntrtr >= end_idx))
                     offtile_list.push_back(ntrtr);
               }
            }
         }


         int nb = nbot_tmp[ic];

         if (nb < start_idx || nb >= end_idx) {
            offtile_list.push_back(nb);
            int ig;
            for (ig = 0; ig < nbsize_local; ig++) {
               if (border_cell_num_global[ig] == nb) break;
            }
            if (ig == nbsize_local) {
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nb);
#ifdef HAVE_MPI
               //L7_Terminate();
#endif
               //exit(-1);
            }
            // Get the bottom neighbor information
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[max(levmx-nlev-1,0)];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];

            if (nlev == levmx) {
               // Multiplying by levmult gets the bottom left cell of the bottom neighbor. Subtract
               // one from that to get the top left cell in the bottom bottom neighbor
               int nbb = hash_tmp[((max((jjj-1)*levmult, 0))-jminsize)*(imaxsize-iminsize)+((iii*levmult)-iminsize)];
               if (nbb != nb && nbb > 0 && (nbb < start_idx || nbb >= end_idx) )
                  offtile_list.push_back(nbb);
            } else {
               int nbb = hash_tmp[(max((jjj*2-1)*levmultminus, 0)-jminsize)*(imaxsize-iminsize)+(iii*levmult-iminsize)];
               if (nbb != nb && nbb > 0 && (nbb < start_idx || nbb >= end_idx) )
                  offtile_list.push_back(nbb);
               int nrbb = hash_tmp[(max((jjj*2-1)*levmultminus, 0)-jminsize)*(imaxsize-iminsize)+((iii*2+1)*levmultminus-iminsize)];
               if (nrbb != nbb && nrbb > 0 && (nrbb < start_idx || nrbb >= end_idx) )
                  offtile_list.push_back(nrbb);
            }

         }

         if (lev < levmx) {
            int levmultminus = levtable[max(levmx-lev-1,0)];
            int nrb = hash_tmp[(max((jj*2-1)*levmultminus,0)-jminsize)*(imaxsize-iminsize)+(min((ii*2+1)*levmultminus,imaxcalc-1)-iminsize)];
            if (nrb != nb && nrb > 0 && (nrb < start_idx || nrb >= end_idx) ) {
               offtile_list.push_back(nrb);

               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == nrb) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nrb);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               int nlev = border_cell_level_global[ig];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];
               int levmult      = levtable[levmx-nlev];
               levmultminus = levtable[max(levmx-nlev-1,0)];

               if (nlev == levmx) {
                  int nbrb = hash_tmp[(max((jjj-1)*levmult,0)-jminsize)*(imaxsize-iminsize)+(iii*levmult-iminsize)];
                  if (nbrb != nrb && nbrb > 0 && (nbrb < start_idx || nbrb >= end_idx) )
                     offtile_list.push_back(nbrb);
               } else {
                  // And we need both cells here
                  int nbrb = hash_tmp[(max((jjj*2-1)*levmultminus,0)-jminsize)*(imaxsize-iminsize)+((iii*2)*levmultminus-iminsize)];
                  if (nbrb != nrb && nbrb > 0 && (nbrb < start_idx || nbrb >= end_idx) )
                     offtile_list.push_back(nbrb);
                  int nrbrb = hash_tmp[(max((jjj*2-1)*levmultminus,0)-jminsize)*(imaxsize-iminsize)+(min((iii*2+1)*levmultminus,imaxcalc-1)-iminsize)];
                  if (nrbrb != nbrb && nrbrb > 0 && (nrbrb < start_idx || nrbrb >= end_idx))
                     offtile_list.push_back(nrbrb);
               }

            }
         }

         int nt = ntop_tmp[ic];

         if (nt < start_idx || nt >= end_idx) {
            offtile_list.push_back(nt);
            int ig;
            for (ig = 0; ig < nbsize_local; ig++) {
               if (border_cell_num_global[ig] == nt) break;
            }
            if (ig == nbsize_local) {
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nt);
#ifdef HAVE_MPI
               //L7_Terminate();
#endif
               //exit(-1);
            }
            int nlev = border_cell_level_global[ig];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[max(levmx-nlev-1,0)];

            if (nlev == levmx) {
               // Multiplying by levmult gets the bottom left cell of the top neighbor. Add
               // one to that to get the bottom left cell of the top top neighbor
               int ntt = hash_tmp[((min((jjj+1)*levmult,jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((iii*levmult)-iminsize)];
               if (ntt != nt && ntt > 0 && (ntt < start_idx || ntt >= end_idx) )
                  offtile_list.push_back(ntt);
            } else {
               int ntt = hash_tmp[((min((jjj+1)*levmult,jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((iii*levmult)-iminsize)];
               if (ntt != nt && ntt > 0 && (ntt < start_idx || ntt >= end_idx) )
                  offtile_list.push_back(ntt);

               int nrtt = hash_tmp[(min((jjj+1)*levmult,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min((iii*2+1)*levmultminus,imaxcalc-1)-iminsize)];
               if (nrtt != ntt && nrtt > 0 && (nrtt < start_idx || nrtt >= end_idx))
                  offtile_list.push_back(nrtt);
            }
         }

         if (lev < levmx) {
            int levmultminus = levtable[max(levmx-lev-1,0)];
            int nrt = hash_tmp[(min((jj+1)*2*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min((ii*2+1)*levmultminus,imaxcalc-1)-iminsize)];
            if (nrt != nt && nrt > 0 && (nrt < start_idx || nrt >= end_idx) ) {
               offtile_list.push_back(nrt);

               int ig;
               for (ig = 0; ig < nbsize_local; ig++) {
                  if (border_cell_num_global[ig] == nrt) break;
               }
               if (ig == nbsize_local) {
                  printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nrt);
#ifdef HAVE_MPI
                  //L7_Terminate();
#endif
                  //exit(-1);
               }
               int nlev = border_cell_level_global[ig];
               int iii = border_cell_i_global[ig];
               int jjj = border_cell_j_global[ig];
               int levmult      = levtable[levmx-nlev];
               levmultminus = levtable[max(levmx-nlev-1,0)];

               // Multiplying by levmult gets the bottom left cell of the top neighbor. Add
               // one to that to get the bottom left cell of the top top neighbor

               if (nlev == levmx) {
                  int ntrt = hash_tmp[(min((jjj+1)*levmult,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(iii*levmult-iminsize)];
                  if (ntrt != nrt && ntrt > 0 && (ntrt < start_idx || ntrt >= end_idx) )
                     offtile_list.push_back(ntrt);
               } else {
                  // And we need both cells here
                  int ntrt = hash_tmp[(min((jjj+1)*2*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+((iii*2)*levmultminus-iminsize)];
                  if (ntrt != nrt && ntrt > 0 && (ntrt < start_idx || ntrt >= end_idx) )
                     offtile_list.push_back(ntrt);
                  int nrtrt = hash_tmp[(min((jjj+1)*2*levmultminus,jmaxcalc-1)-jminsize)*(imaxsize-iminsize)+(min((iii*2+1)*levmultminus,imaxcalc-1)-iminsize)];
                  if (nrtrt != ntrt && nrtrt > 0 && (nrtrt < start_idx || nrtrt >= end_idx))
                     offtile_list.push_back(nrtrt);
               }
            }
         }
      } // ncells loop

      sort(offtile_list.begin(), offtile_list.end());
      vector<int>::iterator p_end = unique(offtile_list.begin(), offtile_list.end());

      vector<int> indices_needed;

      for (vector<int>::iterator p=offtile_list.begin(); p < p_end; p++){
         indices_needed.push_back(*p);
      }

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_offtile_list += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      int nghost = indices_needed.size();
      ncells_ghost = ncells + nghost;
      //fprintf(fp,"%d ncells_ghost size is %ld nghost %d\n",mype,ncells_ghost,nghost);

#ifdef HAVE_MPI
      if (cell_handle) L7_Free(&cell_handle);
#endif
      cell_handle=0;
      //for (int ic=0; ic<nghost; ic++){
      //   fprintf(fp,"%d: indices needed ic %d index %d\n",mype,ic,indices_needed[ic]);
      //}
#ifdef HAVE_MPI
      L7_Setup(0, noffset, ncells, &indices_needed[0], nghost, &cell_handle);

      int num_indices=L7_Get_Num_Indices(cell_handle);

      vector<int>local_indices(num_indices);

      L7_Get_Local_Indices(cell_handle, &local_indices[0]);

/* */
      vector<int> indices_array(ncells,0);

      for (int ic = 0; ic < num_indices; ic++){
         indices_array[local_indices[ic]]++;
      }

      num_indices=0;
      for (int ic = 0; ic < (int)ncells; ic++){
         if (indices_array[ic] > 0) {
            local_indices[num_indices]=ic;
            num_indices++;
         }
      }
/* */

/*
         //if (mype ==1) printf("DEBUG -- num_indices is %d\n",num_indices);

         sort(local_indices.begin(),local_indices.end());
         vector<int>::iterator pt_end = unique(local_indices.begin(),local_indices.end());

         num_indices = 0;
         for (vector<int>::iterator p=local_indices.begin(); p < pt_end; p++){
            num_indices++;
         }
         //printf("%d: border cell size is %d\n",mype,border_cell_num.size());
*/


         //if (mype ==1) printf("DEBUG -- num_indices is %d %d\n",num_indices,local_indices.size());

         int nblocks=0;
         for (int ic = 0; ic < num_indices; ic++){
            for (int in = local_indices[ic]+1; local_indices[ic+1] == in && ic < num_indices; in++,ic++);
            //for (int in = local_indices[ic]+1; local_indices[ic+1] <= in+5 && ic < num_indices; in=local_indices[ic+1],ic++);
            nblocks++;
         }
         vector<int>local_indices_start(nblocks);
         vector<int>local_indices_stop(nblocks);

         int ib=0;
         for (int ic = 0; ic < num_indices; ic++){
            local_indices_start[ib]=local_indices[ic];
            for (int in = local_indices[ic]+1; in == local_indices[ic+1] && ic < num_indices; in++,ic++);
            //for (int in = local_indices[ic]+1; local_indices[ic+1] <= in+5 && ic < num_indices; in=local_indices[ic+1],ic++);
            local_indices_stop[ib]=local_indices[ic];
            ib++;
         }

/*
         for (int ib = 0; ib < nblocks; ib++){
            if (mype == 1) printf("%d: DEBUG block %d start %d stop %d\n",mype,ib,local_indices_start[ib],local_indices_stop[ib]);
         }
*/

         int ic = 0;
         for (int ib = 0; ib < nblocks; ib++){
            for (int in = local_indices_start[ib]; in <= local_indices_stop[ib]; in++){
               //while (in != local_indices[ic]) {in++;}
               if (in != local_indices[ic]) printf("%d: DEBUG ic %d ib %d local_indices[ic] %d block_index %d\n",mype,ic,ib,local_indices[ic],in);
               ic++;
            }
         }
         if (ic != num_indices) printf("%d: DEBUG -- mismatch on indices count\n",mype);


         local_indices.clear();
#endif

      if (TIMING_LEVEL >= 2) {
         ezcl_finish(command_queue);
         gpu_time_setup_comm += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
         cpu_timer_start(&tstart_lev2);
      }

      celltype_tmp.resize(ncells_ghost);
      i_tmp.resize(ncells_ghost);
      j_tmp.resize(ncells_ghost);
      level_tmp.resize(ncells_ghost);
      nlft_tmp.resize(ncells_ghost,-98);
      nrht_tmp.resize(ncells_ghost,-98);
      nbot_tmp.resize(ncells_ghost,-98);
      ntop_tmp.resize(ncells_ghost,-98);

#ifdef HAVE_MPI
      L7_Update(&celltype_tmp[0], L7_INT, cell_handle);
      L7_Update(&i_tmp[0],        L7_INT, cell_handle);
      L7_Update(&j_tmp[0],        L7_INT, cell_handle);
      L7_Update(&level_tmp[0],    L7_INT, cell_handle);
#endif

      //for (int ic=0; ic<ncells; ic++){
      //   fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
      //       mype,ic,i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
      //}

#ifdef HAVE_MPI
      L7_Update(&nlft_tmp[0], L7_INT, cell_handle);
      L7_Update(&nrht_tmp[0], L7_INT, cell_handle);
      L7_Update(&nbot_tmp[0], L7_INT, cell_handle);
      L7_Update(&ntop_tmp[0], L7_INT, cell_handle);
#endif

         //vector<int> itest(ncells_ghost);
         //for (int ic=0; ic<ncells; ic++){
         //   itest[ic] = mype*1000 + ic;
         //   fprintf(fp,"%d: test is filled with ic %d = %d\n",mype,ic,itest[ic]);
         //}
         //for (int ic=ncells; ic<ncells_ghost; ic++){
         //   itest[ic] = 0;
         //}

#ifdef HAVE_MPI
         //L7_Update(&itest[0], L7_INT, cell_handle);
#endif

         //for (int ic=0; ic<ncells_ghost; ic++){
         //   fprintf(fp,"%d: test after update ic %d = %d\n",mype,ic,itest[ic]);
         //}

      if (numpe > 1) {
         cl_mem dev_celltype_old = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_old"), &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_i_old        = ezcl_malloc(NULL, const_cast<char *>("dev_i_old"),        &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_j_old        = ezcl_malloc(NULL, const_cast<char *>("dev_j_old"),        &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_level_old    = ezcl_malloc(NULL, const_cast<char *>("dev_level_old"),    &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nlft_old     = ezcl_malloc(NULL, const_cast<char *>("dev_nlft_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nrht_old     = ezcl_malloc(NULL, const_cast<char *>("dev_nrht_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nbot_old     = ezcl_malloc(NULL, const_cast<char *>("dev_nbot_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_ntop_old     = ezcl_malloc(NULL, const_cast<char *>("dev_ntop_old"),     &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         cl_mem dev_tmp;
         SWAP_PTR(dev_celltype_old, dev_celltype, dev_tmp);
         SWAP_PTR(dev_i_old,        dev_i,        dev_tmp);
         SWAP_PTR(dev_j_old,        dev_j,        dev_tmp);
         SWAP_PTR(dev_level_old,    dev_level,    dev_tmp);
         SWAP_PTR(dev_nlft_old,     dev_nlft,     dev_tmp);
         SWAP_PTR(dev_nrht_old,     dev_nrht,     dev_tmp);
         SWAP_PTR(dev_nbot_old,     dev_nbot,     dev_tmp);
         SWAP_PTR(dev_ntop_old,     dev_ntop,     dev_tmp);

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

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_mesh_data,   1, NULL, &global_work_size, &local_work_size, NULL);
      
         ezcl_device_memory_remove(dev_celltype_old);
         ezcl_device_memory_remove(dev_i_old);
         ezcl_device_memory_remove(dev_j_old);
         ezcl_device_memory_remove(dev_level_old);
         ezcl_device_memory_remove(dev_nlft_old);
         ezcl_device_memory_remove(dev_nrht_old);
         ezcl_device_memory_remove(dev_nbot_old);
         ezcl_device_memory_remove(dev_ntop_old);

         size_t nghost_local = nghost;
         cl_mem dev_nlft_add       = ezcl_malloc(NULL, const_cast<char *>("dev_nlft_add"),       &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nrht_add       = ezcl_malloc(NULL, const_cast<char *>("dev_nrht_add"),       &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_nbot_add       = ezcl_malloc(NULL, const_cast<char *>("dev_nbot_add"),       &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_ntop_add       = ezcl_malloc(NULL, const_cast<char *>("dev_ntop_add"),       &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_celltype_add   = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_add"),   &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_i_add          = ezcl_malloc(NULL, const_cast<char *>("dev_i_add"),          &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_j_add          = ezcl_malloc(NULL, const_cast<char *>("dev_j_add"),          &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_level_add      = ezcl_malloc(NULL, const_cast<char *>("dev_level_add"),      &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         cl_mem dev_indices_needed = ezcl_malloc(NULL, const_cast<char *>("dev_indices_needed"), &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         ezcl_enqueue_write_buffer(command_queue, dev_nlft_add,       CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&nlft_tmp[ncells],     NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_nrht_add,       CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&nrht_tmp[ncells],     NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_nbot_add,       CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&nbot_tmp[ncells],     NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_ntop_add,       CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&ntop_tmp[ncells],     NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_celltype_add,   CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&celltype_tmp[ncells], NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_i_add,          CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&i_tmp[ncells],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_j_add,          CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&j_tmp[ncells],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_level_add,      CL_FALSE, 0, nghost*sizeof(cl_int), (void*)&level_tmp[ncells],    NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_indices_needed, CL_TRUE,  0, nghost*sizeof(cl_int), (void*)&indices_needed[0],    NULL);

         //size_t border_local_work_size = MIN(nbsize,TILE_SIZE);
         size_t ghost_local_work_size = 32;
         size_t ghost_global_work_size = ((nghost + ghost_local_work_size - 1) /ghost_local_work_size) * ghost_local_work_size;

         ezcl_set_kernel_arg(kernel_copy_ghost_data,  0, sizeof(cl_int), (void *)&ncells);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  1, sizeof(cl_int), (void *)&nghost);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  2, sizeof(cl_mem), (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  3, sizeof(cl_mem), (void *)&dev_nlft_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  4, sizeof(cl_mem), (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  5, sizeof(cl_mem), (void *)&dev_nrht_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  6, sizeof(cl_mem), (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  7, sizeof(cl_mem), (void *)&dev_nbot_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  8, sizeof(cl_mem), (void *)&dev_ntop);
         ezcl_set_kernel_arg(kernel_copy_ghost_data,  9, sizeof(cl_mem), (void *)&dev_ntop_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 10, sizeof(cl_mem), (void *)&dev_celltype);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 11, sizeof(cl_mem), (void *)&dev_celltype_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 12, sizeof(cl_mem), (void *)&dev_i);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 13, sizeof(cl_mem), (void *)&dev_i_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 14, sizeof(cl_mem), (void *)&dev_j);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 15, sizeof(cl_mem), (void *)&dev_j_add);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 16, sizeof(cl_mem), (void *)&dev_level);
         ezcl_set_kernel_arg(kernel_copy_ghost_data, 17, sizeof(cl_mem), (void *)&dev_level_add);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_ghost_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, NULL);                                                                                                      
         ezcl_device_memory_remove(dev_nlft_add);
         ezcl_device_memory_remove(dev_nrht_add);
         ezcl_device_memory_remove(dev_nbot_add);
         ezcl_device_memory_remove(dev_ntop_add);
         ezcl_device_memory_remove(dev_celltype_add);
         ezcl_device_memory_remove(dev_i_add);
         ezcl_device_memory_remove(dev_j_add);
         ezcl_device_memory_remove(dev_level_add);

         size_t nc_ghost_local_work_size = 32;
         size_t nc_ghost_global_work_size = ((ncells_ghost + nc_ghost_local_work_size - 1) /nc_ghost_local_work_size) * nc_ghost_local_work_size;

         ezcl_set_kernel_arg(kernel_adjust_neighbors, 0,  sizeof(cl_int), (void *)&ncells_ghost);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 1,  sizeof(cl_int), (void *)&nghost);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 2,  sizeof(cl_int), (void *)&noffset);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 3,  sizeof(cl_mem), (void *)&dev_indices_needed);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 4,  sizeof(cl_mem), (void *)&dev_nlft);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 5,  sizeof(cl_mem), (void *)&dev_nrht);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 6,  sizeof(cl_mem), (void *)&dev_nbot);
         ezcl_set_kernel_arg(kernel_adjust_neighbors, 7,  sizeof(cl_mem), (void *)&dev_ntop);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_adjust_neighbors,   1, NULL, &nc_ghost_global_work_size, &nc_ghost_local_work_size, NULL);

         ezcl_device_memory_remove(dev_indices_needed);
      } // if (numpe > 1)

      if (DEBUG) {
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

   ezcl_device_memory_remove(dev_sizes);

   ezcl_device_memory_remove(dev_hash);

   ezcl_finish(command_queue);

   if (TIMING_LEVEL >= 2) gpu_time_do_ghost_comm += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);

   gpu_time_calc_neighbors += (long)(cpu_timer_stop(tstart_cpu) * 1.0e9);
}
#endif

void Mesh::print_calc_neighbor_type(void)
{
   if ( calc_neighbor_type == HASH_TABLE ) {
      if (mype == 0) printf("Using hash tables to calculate neighbors\n");
   } else {
      if (mype == 0) printf("Using k-D tree to calculate neighbors\n");
   }
}
int Mesh::get_calc_neighbor_type(void)
{
   return(calc_neighbor_type );
}

void Mesh::calc_celltype(void)
{
   celltype.resize(ncells,REAL_CELL);

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
   vector<int> index_list((int)pow(2,levmx*levmx) );

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
void Mesh::do_load_balance_local(const size_t new_ncells, const int &ncells_global, vector<real> &H, vector<real> &U, vector<real> &V)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   int ncells_old = new_ncells;
   int noffset_old = ndispl[mype];

   int do_load_balance_global = 0;
   int nsizes_old = 0;
   for (int ip=0; ip<numpe; ip++){
      nsizes_old = nsizes[ip];
      nsizes[ip] = ncells_global/numpe;
      if (ip < (ncells_global%numpe)) nsizes[ip]++;
      if (nsizes_old != nsizes[ip]) do_load_balance_global = 1;
   }

   if (do_load_balance_global) {

      cpu_load_balance_counter++;

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

      H.resize(ncells_old+indices_needed_count,0.0);
      U.resize(ncells_old+indices_needed_count,0.0);
      V.resize(ncells_old+indices_needed_count,0.0);
      L7_Update(&H[0], L7_REAL, load_balance_handle);
      L7_Update(&U[0], L7_REAL, load_balance_handle);
      L7_Update(&V[0], L7_REAL, load_balance_handle);

      i.resize(ncells_old+indices_needed_count,0);
      j.resize(ncells_old+indices_needed_count,0);
      level.resize(ncells_old+indices_needed_count,0);
      celltype.resize(ncells_old+indices_needed_count,0);
      L7_Update(&i[0], L7_INT, load_balance_handle);
      L7_Update(&j[0], L7_INT, load_balance_handle);
      L7_Update(&level[0], L7_INT, load_balance_handle);
      L7_Update(&celltype[0], L7_INT, load_balance_handle);

      L7_Free(&load_balance_handle);
      load_balance_handle = 0;
 
      vector<real> H_temp(ncells);
      vector<real> U_temp(ncells);
      vector<real> V_temp(ncells);

      vector<int> i_temp(ncells);
      vector<int> j_temp(ncells);
      vector<int> level_temp(ncells);
      vector<int> celltype_temp(ncells);

      vector<int> indexes(ncells);

      in = 0;
      int ic = lower_block_size;
      if(ic > 0) {
         for(; (in < ic) && (in < (int)ncells); in++) {
            H_temp[in] = H[ncells_old + in];
            U_temp[in] = U[ncells_old + in];
            V_temp[in] = V[ncells_old + in];

            i_temp[in]     = i[ncells_old + in];
            j_temp[in]     = j[ncells_old + in];
            level_temp[in] = level[ncells_old + in];
            celltype_temp[in] = celltype[ncells_old + in];
         }
      }

      ic = noffset - noffset_old;
      if(ic < 0) ic = 0;
      for(; (ic < ncells_old) && (in < (int)ncells); ic++, in++) {
         H_temp[in] = H[ic];
         U_temp[in] = U[ic];
         V_temp[in] = V[ic];

         i_temp[in]     = i[ic];
         j_temp[in]     = j[ic];
         level_temp[in] = level[ic];
         celltype_temp[in] = celltype[ic];
      }

      ic = upper_block_size;
      if(ic > 0) {
         ic = ncells_old + lower_block_size;
         for(int k = max(noffset-upper_block_start,0); ((k+ic) < (ncells_old+indices_needed_count)) && (in < (int)ncells); k++, in++) {
            H_temp[in] = H[ic+k];
            U_temp[in] = U[ic+k];
            V_temp[in] = V[ic+k];

            i_temp[in]     = i[ic+k];
            j_temp[in]     = j[ic+k];
            level_temp[in] = level[ic+k];
            celltype_temp[in] = celltype[ic+k];
         }
      }

      H.swap(H_temp);
      U.swap(U_temp);
      V.swap(V_temp);

      i.swap(i_temp);
      j.swap(j_temp);
      level.swap(level_temp);
      celltype.swap(celltype_temp);

   }

   cpu_time_load_balance += cpu_timer_stop(tstart_cpu);
}
#endif

#ifdef HAVE_OPENCL
#ifdef HAVE_MPI
void Mesh::gpu_do_load_balance_local(cl_command_queue command_queue, const size_t new_ncells, const int &ncells_global, cl_mem &dev_H, cl_mem &dev_U, cl_mem &dev_V)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   int ncells_old = new_ncells;
   int noffset_old = ndispl[mype];

   int do_load_balance_global = 0;
   int nsizes_old = 0;
   for (int ip=0; ip<numpe; ip++){
      nsizes_old = nsizes[ip];
      nsizes[ip] = ncells_global/numpe;
      if (ip < (ncells_global%numpe)) nsizes[ip]++;
      if (nsizes_old != nsizes[ip]) do_load_balance_global = 1;
   }

   if(do_load_balance_global) {

      gpu_load_balance_counter++;

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
      //int lower_segment_start = noffset-noffset;
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

      // Read current H, U, V, i, j, level, celltype values from GPU and write to CPU arrays
      // Update arrays with L7
       
      vector<real> H_tmp(ncells_old+indices_needed_count,0.0);
      vector<real> U_tmp(ncells_old+indices_needed_count,0.0);
      vector<real> V_tmp(ncells_old+indices_needed_count,0.0);

      vector<int> i_tmp(ncells_old+indices_needed_count,0);
      vector<int> j_tmp(ncells_old+indices_needed_count,0);
      vector<int> level_tmp(ncells_old+indices_needed_count,0);
      vector<int> celltype_tmp(ncells_old+indices_needed_count,0);

      if (do_whole_segment) {
         ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells_old*sizeof(cl_real), &H_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells_old*sizeof(cl_real), &U_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_V, CL_FALSE, 0, ncells_old*sizeof(cl_real), &V_tmp[0], NULL);

         ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, ncells_old*sizeof(cl_int), &i_tmp[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, ncells_old*sizeof(cl_int), &j_tmp[0],        NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells_old*sizeof(cl_int), &level_tmp[0],    NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  0, ncells_old*sizeof(cl_int), &celltype_tmp[0], NULL);
      } else {
         if (lower_segment_size > 0) {
            ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, lower_segment_size*sizeof(cl_real), &H_tmp[0], NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, lower_segment_size*sizeof(cl_real), &U_tmp[0], NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_V, CL_FALSE, 0, lower_segment_size*sizeof(cl_real), &V_tmp[0], NULL);

            ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, 0, lower_segment_size*sizeof(cl_int), &i_tmp[0],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, 0, lower_segment_size*sizeof(cl_int), &j_tmp[0],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, lower_segment_size*sizeof(cl_int), &level_tmp[0],    NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  0, lower_segment_size*sizeof(cl_int), &celltype_tmp[0], NULL);
         }
         if (upper_segment_size > 0) {
            ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, upper_segment_start*sizeof(cl_real), upper_segment_size*sizeof(cl_real), &H_tmp[upper_segment_start], NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, upper_segment_start*sizeof(cl_real), upper_segment_size*sizeof(cl_real), &U_tmp[upper_segment_start], NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_V, CL_FALSE, upper_segment_start*sizeof(cl_real), upper_segment_size*sizeof(cl_real), &V_tmp[upper_segment_start], NULL);

            ezcl_enqueue_read_buffer(command_queue, dev_i,        CL_FALSE, upper_segment_start*sizeof(cl_int), upper_segment_size*sizeof(cl_int), &i_tmp[upper_segment_start],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_j,        CL_FALSE, upper_segment_start*sizeof(cl_int), upper_segment_size*sizeof(cl_int), &j_tmp[upper_segment_start],        NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, upper_segment_start*sizeof(cl_int), upper_segment_size*sizeof(cl_int), &level_tmp[upper_segment_start],    NULL);
            ezcl_enqueue_read_buffer(command_queue, dev_celltype, CL_TRUE,  upper_segment_start*sizeof(cl_int), upper_segment_size*sizeof(cl_int), &celltype_tmp[upper_segment_start], NULL);
         }
      }

      L7_Update(&H_tmp[0], L7_REAL, load_balance_handle);
      L7_Update(&U_tmp[0], L7_REAL, load_balance_handle);
      L7_Update(&V_tmp[0], L7_REAL, load_balance_handle);

      L7_Update(&i_tmp[0],        L7_INT, load_balance_handle);
      L7_Update(&j_tmp[0],        L7_INT, load_balance_handle);
      L7_Update(&level_tmp[0],    L7_INT, load_balance_handle);
      L7_Update(&celltype_tmp[0], L7_INT, load_balance_handle);

      L7_Free(&load_balance_handle);
      load_balance_handle = 0;

      // Allocate and set lower block on GPU
      cl_mem dev_H_lower, dev_U_lower, dev_V_lower, dev_i_lower, dev_j_lower, dev_level_lower, dev_celltype_lower;

      if(lower_block_size > 0) {
         dev_H_lower = ezcl_malloc(NULL, const_cast<char *>("dev_H_lower"), &lower_block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
         dev_U_lower = ezcl_malloc(NULL, const_cast<char *>("dev_U_lower"), &lower_block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
         dev_V_lower = ezcl_malloc(NULL, const_cast<char *>("dev_V_lower"), &lower_block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

         dev_i_lower        = ezcl_malloc(NULL, const_cast<char *>("dev_i_lower"),        &lower_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_j_lower        = ezcl_malloc(NULL, const_cast<char *>("dev_j_lower"),        &lower_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_level_lower    = ezcl_malloc(NULL, const_cast<char *>("dev_level_lower"),    &lower_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_celltype_lower = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_lower"), &lower_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         ezcl_enqueue_write_buffer(command_queue, dev_H_lower, CL_FALSE, 0, lower_block_size*sizeof(cl_real), &H_tmp[ncells_old], NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_U_lower, CL_FALSE, 0, lower_block_size*sizeof(cl_real), &U_tmp[ncells_old], NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V_lower, CL_FALSE, 0, lower_block_size*sizeof(cl_real), &V_tmp[ncells_old], NULL);

         ezcl_enqueue_write_buffer(command_queue, dev_i_lower,        CL_FALSE, 0, lower_block_size*sizeof(cl_int), &i_tmp[ncells_old],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_j_lower,        CL_FALSE, 0, lower_block_size*sizeof(cl_int), &j_tmp[ncells_old],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_level_lower,    CL_FALSE, 0, lower_block_size*sizeof(cl_int), &level_tmp[ncells_old],    NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_celltype_lower, CL_TRUE,  0, lower_block_size*sizeof(cl_int), &celltype_tmp[ncells_old], NULL);
      }

      // Allocate and set upper block on GPU
      cl_mem dev_H_upper, dev_U_upper, dev_V_upper, dev_i_upper, dev_j_upper, dev_level_upper, dev_celltype_upper;
      if(upper_block_size > 0) {
         dev_H_upper = ezcl_malloc(NULL, const_cast<char *>("dev_H_upper"), &upper_block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
         dev_U_upper = ezcl_malloc(NULL, const_cast<char *>("dev_U_upper"), &upper_block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
         dev_V_upper = ezcl_malloc(NULL, const_cast<char *>("dev_V_upper"), &upper_block_size, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

         dev_i_upper        = ezcl_malloc(NULL, const_cast<char *>("dev_i_upper"),        &upper_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_j_upper        = ezcl_malloc(NULL, const_cast<char *>("dev_j_upper"),        &upper_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_level_upper    = ezcl_malloc(NULL, const_cast<char *>("dev_level_upper"),    &upper_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         dev_celltype_upper = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_upper"), &upper_block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         ezcl_enqueue_write_buffer(command_queue, dev_H_upper, CL_FALSE, 0, upper_block_size*sizeof(cl_real), &H_tmp[ncells_old+lower_block_size], NULL); 
         ezcl_enqueue_write_buffer(command_queue, dev_U_upper, CL_FALSE, 0, upper_block_size*sizeof(cl_real), &U_tmp[ncells_old+lower_block_size], NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_V_upper, CL_TRUE,  0, upper_block_size*sizeof(cl_real), &V_tmp[ncells_old+lower_block_size], NULL);

         ezcl_enqueue_write_buffer(command_queue, dev_i_upper,        CL_FALSE, 0, upper_block_size*sizeof(cl_int), &i_tmp[ncells_old+lower_block_size],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_j_upper,        CL_FALSE, 0, upper_block_size*sizeof(cl_int), &j_tmp[ncells_old+lower_block_size],        NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_level_upper,    CL_FALSE, 0, upper_block_size*sizeof(cl_int), &level_tmp[ncells_old+lower_block_size],    NULL);
         ezcl_enqueue_write_buffer(command_queue, dev_celltype_upper, CL_TRUE,  0, upper_block_size*sizeof(cl_int), &celltype_tmp[ncells_old+lower_block_size], NULL);
      }

      size_t local_work_size = 128;

      // printf("MYPE%d: \t ncells = %d \t ncells_old = %d \t ncells_global = %d \n", mype, ncells, ncells_old, ncells_global);
      // Allocate space on GPU for temp arrays (used in double buffering)
      cl_mem dev_H_new        = ezcl_malloc(NULL, const_cast<char *>("dev_H_new"),        &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      cl_mem dev_U_new        = ezcl_malloc(NULL, const_cast<char *>("dev_U_new"),        &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);
      cl_mem dev_V_new        = ezcl_malloc(NULL, const_cast<char *>("dev_V_new"),        &ncells, sizeof(cl_real), CL_MEM_READ_WRITE, 0);

      cl_mem dev_i_new        = ezcl_malloc(NULL, const_cast<char *>("dev_i_new"),        &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_j_new        = ezcl_malloc(NULL, const_cast<char *>("dev_j_new"),        &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_level_new    = ezcl_malloc(NULL, const_cast<char *>("dev_level_new"),    &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
      cl_mem dev_celltype_new = ezcl_malloc(NULL, const_cast<char *>("dev_celltype_new"), &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

      // Set kernel arguments and call lower block kernel
      if(lower_block_size > 0) {

         size_t global_work_size = ((lower_block_size + local_work_size - 1) / local_work_size) * local_work_size;

         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 0,  sizeof(cl_mem), &dev_H_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 1,  sizeof(cl_mem), &dev_U_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 2,  sizeof(cl_mem), &dev_V_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 3,  sizeof(cl_mem), &dev_i_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 4,  sizeof(cl_mem), &dev_j_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 5,  sizeof(cl_mem), &dev_level_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 6,  sizeof(cl_mem), &dev_celltype_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 7,  sizeof(cl_mem), &dev_H_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 8,  sizeof(cl_mem), &dev_U_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 9,  sizeof(cl_mem), &dev_V_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 10, sizeof(cl_mem), &dev_i_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 11, sizeof(cl_mem), &dev_j_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 12, sizeof(cl_mem), &dev_level_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 13, sizeof(cl_mem), &dev_celltype_lower);
         ezcl_set_kernel_arg(kernel_do_load_balance_lower, 14, sizeof(cl_int), &lower_block_size);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_lower,   1, NULL, &global_work_size, &local_work_size, NULL);

         ezcl_device_memory_remove(dev_H_lower);
         ezcl_device_memory_remove(dev_U_lower);
         ezcl_device_memory_remove(dev_V_lower);
         ezcl_device_memory_remove(dev_i_lower);
         ezcl_device_memory_remove(dev_j_lower);
         ezcl_device_memory_remove(dev_level_lower);
         ezcl_device_memory_remove(dev_celltype_lower);
      }

      // Set kernel arguments and call middle block kernel
      if(middle_block_size > 0) {

         size_t global_work_size = ((middle_block_size + local_work_size - 1) / local_work_size) * local_work_size;

         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 0,  sizeof(cl_mem), &dev_H_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 1,  sizeof(cl_mem), &dev_U_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 2,  sizeof(cl_mem), &dev_V_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 3,  sizeof(cl_mem), &dev_i_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 4,  sizeof(cl_mem), &dev_j_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 5,  sizeof(cl_mem), &dev_level_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 6,  sizeof(cl_mem), &dev_celltype_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 7,  sizeof(cl_mem), &dev_H);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 8,  sizeof(cl_mem), &dev_U);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 9,  sizeof(cl_mem), &dev_V);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 10, sizeof(cl_mem), &dev_i);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 11, sizeof(cl_mem), &dev_j);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 12, sizeof(cl_mem), &dev_level);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 13, sizeof(cl_mem), &dev_celltype);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 14, sizeof(cl_int), &lower_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 15, sizeof(cl_int), &middle_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_middle, 16, sizeof(cl_int), &middle_block_start);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_middle,   1, NULL, &global_work_size, &local_work_size, NULL);
      }

      // Set kernel arguments and call upper block kernel
      if(upper_block_size > 0) {

         size_t global_work_size = ((upper_block_size + local_work_size - 1) / local_work_size) * local_work_size;

         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 0,  sizeof(cl_mem), &dev_H_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 1,  sizeof(cl_mem), &dev_U_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 2,  sizeof(cl_mem), &dev_V_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 3,  sizeof(cl_mem), &dev_i_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 4,  sizeof(cl_mem), &dev_j_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 5,  sizeof(cl_mem), &dev_level_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 6,  sizeof(cl_mem), &dev_celltype_new);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 7,  sizeof(cl_mem), &dev_H_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 8,  sizeof(cl_mem), &dev_U_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 9,  sizeof(cl_mem), &dev_V_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 10, sizeof(cl_mem), &dev_i_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 11, sizeof(cl_mem), &dev_j_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 12, sizeof(cl_mem), &dev_level_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 13, sizeof(cl_mem), &dev_celltype_upper);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 14, sizeof(cl_int), &lower_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 15, sizeof(cl_int), &middle_block_size);
         ezcl_set_kernel_arg(kernel_do_load_balance_upper, 16, sizeof(cl_int), &upper_block_size);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_do_load_balance_upper,   1, NULL, &global_work_size, &local_work_size, NULL);

         ezcl_device_memory_remove(dev_H_upper);
         ezcl_device_memory_remove(dev_U_upper);
         ezcl_device_memory_remove(dev_V_upper);
         ezcl_device_memory_remove(dev_i_upper);
         ezcl_device_memory_remove(dev_j_upper);
         ezcl_device_memory_remove(dev_level_upper);
         ezcl_device_memory_remove(dev_celltype_upper);
      }

      cl_mem dev_ptr = NULL;
      SWAP_PTR(dev_H_new, dev_H, dev_ptr);
      SWAP_PTR(dev_U_new, dev_U, dev_ptr);
      SWAP_PTR(dev_V_new, dev_V, dev_ptr);

      SWAP_PTR(dev_i_new, dev_i, dev_ptr);
      SWAP_PTR(dev_j_new, dev_j, dev_ptr);
      SWAP_PTR(dev_level_new, dev_level, dev_ptr);
      SWAP_PTR(dev_celltype_new, dev_celltype, dev_ptr);

      ezcl_device_memory_remove(dev_H_new);
      ezcl_device_memory_remove(dev_U_new);
      ezcl_device_memory_remove(dev_V_new);

      ezcl_device_memory_remove(dev_i_new);
      ezcl_device_memory_remove(dev_j_new);
      ezcl_device_memory_remove(dev_level_new);
      ezcl_device_memory_remove(dev_celltype_new);

      ezcl_finish(command_queue);

      gpu_time_load_balance += (long int)(cpu_timer_stop(tstart_cpu)*1.0e9);
   }
}
#endif
#endif

#ifdef HAVE_OPENCL
int Mesh::gpu_count_BCs(cl_command_queue command_queue)
{
   cl_event count_BCs_stage1_event, count_BCs_stage2_event;

   size_t local_work_size  = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

   //size_t block_size = (ncells + TILE_SIZE - 1) / TILE_SIZE; //  For on-device global reduction kernel.
   size_t block_size     = global_work_size/local_work_size;

   int bcount = 0;

   if (! have_boundary) {
      cl_mem dev_ioffset  = ezcl_malloc(NULL, const_cast<char *>("dev_ioffset"), &block_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

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

      ezcl_device_memory_remove(dev_ioffset);

      gpu_time_count_BCs        += ezcl_timer_calc(&count_BCs_stage1_event, &count_BCs_stage1_event);
      if (block_size > 1) {
         gpu_time_count_BCs     += ezcl_timer_calc(&count_BCs_stage2_event, &count_BCs_stage2_event);
      }

   }

   return(bcount);
}
#endif
void Mesh::resize_old_device_memory(size_t ncells)
{
#ifdef HAVE_OPENCL
   ezcl_device_memory_remove(dev_level);
   ezcl_device_memory_remove(dev_i);
   ezcl_device_memory_remove(dev_j);
   ezcl_device_memory_remove(dev_celltype);
   dev_level    = ezcl_malloc(NULL, const_cast<char *>("dev_level"),    &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_i        = ezcl_malloc(NULL, const_cast<char *>("dev_i"),        &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_j        = ezcl_malloc(NULL, const_cast<char *>("dev_j"),        &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_celltype = ezcl_malloc(NULL, const_cast<char *>("dev_celltype"), &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
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
   printf("vector celltype  ptr : %p nelements %ld elsize %ld\n",&celltype[0],celltype.size(),sizeof(celltype[0])); 
   printf("vector level     ptr : %p nelements %ld elsize %ld\n",&level[0],   level.size(),   sizeof(level[0])); 
   printf("vector i         ptr : %p nelements %ld elsize %ld\n",&i[0],       i.size(),       sizeof(i[0])); 
   printf("vector j         ptr : %p nelements %ld elsize %ld\n",&j[0],       j.size(),       sizeof(j[0])); 

   printf("vector nlft      ptr : %p nelements %ld elsize %ld\n",&nlft[0],    nlft.size(),    sizeof(nlft[0])); 
   printf("vector nrht      ptr : %p nelements %ld elsize %ld\n",&nrht[0],    nrht.size(),    sizeof(nrht[0])); 
   printf("vector nbot      ptr : %p nelements %ld elsize %ld\n",&nbot[0],    nbot.size(),    sizeof(nbot[0])); 
   printf("vector ntop      ptr : %p nelements %ld elsize %ld\n",&ntop[0],    ntop.size(),    sizeof(ntop[0])); 
}
