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
#include "mesh.h"
#include "reorder.h"
#include "ezcl/ezcl.h"
#include "timer/timer.h"
#include "mpi.h"
#include "l7/l7.h"

#define DEBUG 0
#define NEIGHBOR_CHECK 0

#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef HAVE_CL_DOUBLE
typedef double      real;
#define MPI_C_REAL MPI_DOUBLE
#else
typedef float       real;
#define MPI_C_REAL MPI_FLOAT
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
cl_kernel      kernel_calc_neighbors_local2;
cl_kernel      kernel_copy_mesh_data;
cl_kernel      kernel_copy_ghost_data;
cl_kernel      kernel_adjust_neighbors;
cl_kernel      kernel_reduction_scan;
cl_kernel      kernel_hash_size;
cl_kernel      kernel_finish_hash_size;
cl_kernel      kernel_calc_spatial_coordinates;

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

Mesh::Mesh(int nx, int ny, int levmx_in, int ndim_in, int numpe_in, int boundary, int parallel_in, int do_gpu_calc)
{
   cpu_time_calc_neighbors     = 0.0;
   cpu_time_rezone_all         = 0.0;
   cpu_time_partition          = 0.0;
   cpu_time_calc_spatial_coordinates = 0.0;

   gpu_time_reduction_scan     = 0;
   gpu_time_hash_setup         = 0;
   gpu_time_calc_neighbors     = 0;
   gpu_time_rezone_all         = 0;
   gpu_time_calc_spatial_coordinates = 0;

   ndim   = ndim_in;
   levmx  = levmx_in;
   numpe  = numpe_in;
   mype   = 0;

   offtile_ratio_local = 0;
   offtile_local_count = 1;

   parallel = parallel_in;
#ifdef HAVE_MPI
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init){
      MPI_Comm_rank(MPI_COMM_WORLD,&mype);
      MPI_Comm_size(MPI_COMM_WORLD,&numpe);
   } else {
      mype  = 0;
      numpe = 1;
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

   // The copy host ptr flag will have the data copied to the GPU as part of the allocation
   if (do_gpu_calc) {
      dev_levtable = ezcl_malloc(&levtable[0],   &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levdx    = ezcl_malloc(&lev_deltax[0], &lvlMxSize, sizeof(cl_real), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levdy    = ezcl_malloc(&lev_deltay[0], &lvlMxSize, sizeof(cl_real), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levibeg  = ezcl_malloc(&lev_ibegin[0], &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_leviend  = ezcl_malloc(&lev_iend[0],   &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levjbeg  = ezcl_malloc(&lev_jbegin[0], &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      dev_levjend  = ezcl_malloc(&lev_jend[0],   &lvlMxSize, sizeof(cl_int),  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
   }

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

void Mesh::init(int nx, int ny, double circ_radius, cl_context context, partition_method initial_order, bool special_case, int do_gpu_calc)
{
   if (do_gpu_calc) {
      kernel_reduction_scan    = ezcl_create_kernel(context, "wave_kern_calc.cl", "finish_reduction_scan_cl", 0);
      kernel_hash_init         = ezcl_create_kernel(context, "wave_kern.cl",      "hash_init_cl",             0);
      kernel_hash_init_corners = ezcl_create_kernel(context, "wave_kern.cl",      "hash_init_corners_cl",      0);
      kernel_hash_setup        = ezcl_create_kernel(context, "wave_kern.cl",      "hash_setup_cl",            0);
      kernel_hash_setup_local  = ezcl_create_kernel(context, "wave_kern.cl",      "hash_setup_local_cl",      0);
      kernel_hash_setup_border = ezcl_create_kernel(context, "wave_kern.cl",      "hash_setup_border_cl",      0);
      kernel_calc_neighbors    = ezcl_create_kernel(context, "wave_kern.cl",      "calc_neighbors_cl",        0);
      kernel_calc_neighbors_local = ezcl_create_kernel(context, "wave_kern.cl",      "calc_neighbors_local_cl",        0);
      kernel_calc_neighbors_local2 = ezcl_create_kernel(context, "wave_kern.cl",      "calc_neighbors_local2_cl",        0);
      kernel_copy_mesh_data = ezcl_create_kernel(context, "wave_kern.cl",      "copy_mesh_data_cl",        0);
      kernel_copy_ghost_data = ezcl_create_kernel(context, "wave_kern.cl",      "copy_ghost_data_cl",        0);
      kernel_adjust_neighbors = ezcl_create_kernel(context, "wave_kern.cl",      "adjust_neighbors_cl",        0);
      kernel_hash_size         = ezcl_create_kernel(context, "wave_kern.cl",      "calc_hash_size_cl",        0);
      kernel_finish_hash_size  = ezcl_create_kernel(context, "wave_kern.cl",      "finish_reduction_minmax4_cl",        0);
      kernel_calc_spatial_coordinates = ezcl_create_kernel(context, "wave_kern.cl",      "calc_spatial_coordinates_cl",        0);
   }

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

         index[ic] = ic;
         i[ic]     = ii;
         j[ic]     = jj;
         level[ic] = 0;
         ic++;
      }
   }

// ibase = 0;
   calc_spatial_coordinates(ibase);
   calc_minmax();

   nlft.resize(ncells);
   nrht.resize(ncells);
   nbot.resize(ncells);
   ntop.resize(ncells);

   index.resize(ncells);
   for (uint ic=0; ic<ncells; ic++) {
      nlft[ic]=-1;
      nrht[ic]=-1;
      nbot[ic]=-1;
      ntop[ic]=-1;
   }

   celltype.resize(ncells);

   calc_celltype();

   calc_neighbors();

   partition_cells(numpe, proc, index, initial_order);

   //  Start lev loop here
   for (int ilevel=1; ilevel<=levmx; ilevel++) {

      vector<int> mpot(ncells);

      for (uint ic=0; ic<ncells; ++ic) {
         mpot[ic]=0;
      }

      calc_neighbors();

      kdtree_setup();

      int nez;
      vector<int> ind(ncells);

      KDTree_QueryCircleIntersect(&tree, &nez, &(ind[0]), circ_radius, ncells, &x[0], &dx[0], &y[0], &dy[0]);
      for (int i=0; i<nez; ++i){
         if (level[ind[i]] < levmx) mpot[ind[i]] = 1;
      }

      KDTree_Destroy(&tree);
      //  Refine the cells.
      if (! special_case) rezone_spread(mpot);
      int add_ncells = rezone_count(mpot);
      rezone_all(mpot, add_ncells);
   }  // End lev loop here
   
   calc_spatial_coordinates(0);

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

void Mesh::gpu_rezone_count(cl_command_queue command_queue, size_t block_size, size_t local_work_size,
    cl_mem dev_ioffset, cl_mem &dev_result)
{
   cl_event reduction_scan_event;

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

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_scan, 1, NULL, &local_work_size, &local_work_size, &reduction_scan_event);

   gpu_time_reduction_scan    += ezcl_timer_calc(&reduction_scan_event,    &reduction_scan_event);
}

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
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  0, sizeof(cl_int), (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  1, sizeof(cl_real), (void *)&xmin);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  2, sizeof(cl_real), (void *)&ymin);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  3, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  4, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  5, sizeof(cl_mem),  (void *)&dev_x);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  6, sizeof(cl_mem),  (void *)&dev_dx);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  7, sizeof(cl_mem),  (void *)&dev_y);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  8, sizeof(cl_mem),  (void *)&dev_dy);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates,  9, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates, 10, sizeof(cl_mem),  (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_calc_spatial_coordinates, 11, sizeof(cl_mem),  (void *)&dev_j);
   //ezcl_set_kernel_arg(kernel_hash_init, 1, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_spatial_coordinates, 1, NULL, &global_work_size, &local_work_size, &calc_spatial_coordinates_event);
/*
      for (uint ic = 0; ic < ncells; ic++) {
         x[ic]  = xmin + lev_deltax[level[ic]] * (real)(i[ic] - ibase);
         dx[ic] =        lev_deltax[level[ic]];
         y[ic]  = ymin + lev_deltay[level[ic]] * (real)(j[ic] - ibase);
         dy[ic] =        lev_deltay[level[ic]];
      }
*/
   gpu_time_calc_spatial_coordinates += ezcl_timer_calc(&calc_spatial_coordinates_event, &calc_spatial_coordinates_event);
}

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

}

void Mesh::rezone_all(vector<int> mpot, int add_ncells)
{
   struct timeval tstart_cpu;

   uint ic,          //  Index for old cell arrays.
        nc;          //  Index for new cell arrays.
   int set_index = 0;

   cpu_timer_start(&tstart_cpu);

   //  Check for requested mesh refinements; if there are none, return.
   if (parallel) {
      int global_add_ncells;
      MPI_Allreduce(&add_ncells, &global_add_ncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (global_add_ncells == 0) set_index = 1;
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
      MPI_Request req[12];
      MPI_Status status[12];

      static unsigned int prev     = MPI_PROC_NULL;
      static unsigned int next     = MPI_PROC_NULL;

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
   
   nlft.empty();
   nrht.empty();
   nbot.empty();
   ntop.empty();
   nlft.resize(new_ncells, -1);
   nrht.resize(new_ncells, -1);
   nbot.resize(new_ncells, -1);
   ntop.resize(new_ncells, -1);
   
   //ibase    = 0;
   //calc_spatial_coordinates(ibase);

   calc_celltype();

   cpu_time_rezone_all += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_neighbors(void)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   nlft.resize(ncells);
   nrht.resize(ncells);
   nbot.resize(ncells);
   ntop.resize(ncells);

   if (calc_neighbor_type == HASH_TABLE) {

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
            for (int jj = j[ic]*levtable[levmx-lev]; jj < (j[ic]+1)*levtable[levmx-lev]; jj++) {
               for (int ii=i[ic]*levtable[levmx-lev]; ii<(i[ic]+1)*levtable[levmx-lev]; ii++) {
                  hash[jj][ii] = ic;
               }
            }
         }
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

   } else if (calc_neighbor_type == KDTREE) {
      TBounds box;
      vector<int> index_list(20);

      int num;

      ibase = 0;
      calc_spatial_coordinates(ibase);

      kdtree_setup();

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

   } // calc_neighbor_type

   cpu_time_calc_neighbors += cpu_timer_stop(tstart_cpu);
}

void Mesh::calc_neighbors_local(void)
{
   struct timeval tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   nlft.resize(ncells,-98);
   nrht.resize(ncells,-98);
   nbot.resize(ncells,-98);
   ntop.resize(ncells,-98);


   if (calc_neighbor_type == HASH_TABLE) {

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
      // TODO: This has changed -- need to update kernel code
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
            for (int    jj = j[ic]*levtable[levmx-lev]-jminsize; jj < (j[ic]+1)*levtable[levmx-lev]-jminsize; jj++) {
               for (int ii = i[ic]*levtable[levmx-lev]-iminsize; ii < (i[ic]+1)*levtable[levmx-lev]-iminsize; ii++) {
                  //printf("%d: block j %d i %d\n",mype,jj,ii);
                  hash[jj][ii] = ic+noffset;
               }
            }
         }
      }
      
      int ii, jj, lev, iii, jjj, levmult;

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

      vector<int> border_cell;

      // Push list of unsatisfied neighbor cells
      for (uint ic=0; ic<ncells; ic++){
         if (nlft[ic] == -1){
            //printf("%d: Cell is %d nlft %d\n",mype,ic+noffset,nlft[ic]);
            border_cell.push_back(ic+noffset);
            if (nrht[ic] >= 0) {
               border_cell.push_back(nrht[ic]);
               if (level[nrht[ic]-noffset] > level[ic]) {
                  if (ntop[nrht[ic]-noffset] >= 0) border_cell.push_back(ntop[nrht[ic]-noffset]);
               }
            }
         }
         if (nrht[ic] == -1){
            //printf("%d: Cell is %d nrht %d\n",mype,ic+noffset,nrht[ic]);
            border_cell.push_back(ic+noffset);
            if (nlft[ic] >= 0) {
               border_cell.push_back(nlft[ic]);
               if (level[nlft[ic]-noffset] > level[ic]) {
                  if (ntop[nlft[ic]-noffset] >= 0) border_cell.push_back(ntop[nlft[ic]-noffset]);
               }
            }
         }
         if (nbot[ic] == -1) {
            //printf("%d: Cell is %d nbot %d\n",mype,ic+noffset,nbot[ic]);
            border_cell.push_back(ic+noffset);
            if (ntop[ic] >= 0) {
               border_cell.push_back(ntop[ic]);
               if (level[ntop[ic]-noffset] > level[ic]) {
                  if (nrht[ntop[ic]-noffset] >= 0) border_cell.push_back(nrht[ntop[ic]-noffset]);
               }
            }
         }
         if (ntop[ic] == -1) {
            //printf("%d: Cell is %d ntop %d\n",mype,ic+noffset,ntop[ic]);
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

      int nbsize_local, nbsize_global;
      nbsize_local=border_cell_num.size();

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

      MPI_Allreduce(&nbsize_local, &nbsize_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      vector<int> nbsizes(numpe);
      vector<int> nbdispl(numpe);
      MPI_Allgather(&nbsize_local, 1, MPI_INT, &nbsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
      nbdispl[0]=0;
      for (int ip=1; ip<numpe; ip++){
         nbdispl[ip] += nbdispl[ip-1] + nbsizes[ip-1];
      }

      vector<int>border_cell_num_global(nbsize_global);
      vector<int>border_cell_i_global(nbsize_global);
      vector<int>border_cell_j_global(nbsize_global);
      vector<int>border_cell_level_global(nbsize_global);

      MPI_Allgatherv(&border_cell_num[0],   nbsizes[mype], MPI_INT, &border_cell_num_global[0],   &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&border_cell_i[0],     nbsizes[mype], MPI_INT, &border_cell_i_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&border_cell_j[0],     nbsizes[mype], MPI_INT, &border_cell_j_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&border_cell_level[0], nbsizes[mype], MPI_INT, &border_cell_level_global[0], &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);

      //for (int ic = 0; ic < nbsize_global; ic++) {
      //   fprintf(fp,"%d: Global Border cell %d is %d i %d j %d level %d\n",mype,ic,border_cell_global[ic],
      //      border_cell_i_global[ic],border_cell_j_global[ic],border_cell_level_global[ic]);
      //}

      int inew=0;
      for (int ic = 0; ic < nbsize_global; ic++) {
         int lev = border_cell_level_global[ic];
         int levmult = levtable[levmx-lev];
         //fprintf(fp,"%d: DEBUG cell %d i %d j %d\n",mype,ic,border_cell_i_global[ic],border_cell_j_global[ic]);
         if (border_cell_j_global[ic]*levmult < jminsize || border_cell_j_global[ic]*levmult >= jmaxsize) continue;
         if (border_cell_i_global[ic]*levmult < iminsize || border_cell_i_global[ic]*levmult >= imaxsize) continue;
         if (border_cell_num_global[ic] >= noffset && border_cell_num_global[ic] < noffset+ncells) continue;
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
            for (int    jj = border_cell_j_global[ic]*levmult-jminsize; jj < (border_cell_j_global[ic]+1)*levmult-jminsize; jj++) {
               for (int ii = border_cell_i_global[ic]*levmult-iminsize; ii < (border_cell_i_global[ic]+1)*levmult-iminsize; ii++) {
                  //fprintf(fp,"%d: cell %d block j %d i %d\n",mype,ic,jj,ii);
                  hash[jj][ii] = border_cell_num_global[ic];
               }
            }
         }
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
         if (nrht[ic] == -1) nrht[ic] = hash[(      jj   *levmult               )-jminsize][(min( (ii+1)*levmult,   imaxcalc-1))-iminsize];
         if (nbot[ic] == -1) nbot[ic] = hash[(max(  jj   *levmult-1, 0)         )-jminsize][(      ii   *levmult               )-iminsize];
         if (ntop[ic] == -1) ntop[ic] = hash[(min( (jj+1)*levmult,   jmaxcalc-1))-jminsize][(      ii   *levmult               )-iminsize];
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nl);
               L7_Terminate();
               exit(-1);
            }
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];
            //fprintf(fp,"%d: DEBUG nl is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nl, ig, border_cell_num_global[ig], nlev, iii, jjj);
            int nll = hash[(      jjj   *levmult               )-jminsize][(max(  iii   *levmult-1, 0         ))-iminsize];
            //printf("%d: ic %d nll is %d\n",mype,ic,nll);
            if (nll != nl && nll > 0 && (nll < start_idx || nll >= end_idx) ) offtile_list.push_back(nll);
            if (nlev < levmx) {
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[levmx-nlev-1];

               int ntll = hash[(min( (jjj*2+1)*levmultminus, jmaxcalc-1))-jminsize][(max(  iii     *levmult-1,    0))         -iminsize];
               if (ntll != nll && ntll > 0 && (ntll < start_idx || ntll >= end_idx)) offtile_list.push_back(ntll);
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
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nr);
               L7_Terminate();
               exit(-1);
            }
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];
            //fprintf(fp,"%d: DEBUG nr is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nr, ig, border_cell_num_global[ig], nlev, iii, jjj);
            int nrr = hash[(      jjj   *levmult               )-jminsize][(min( (iii+1)*levmult,   imaxcalc-1))-iminsize];
            if (nrr != nr && nrr > 0 && (nrr < start_idx || nrr >= end_idx) ) offtile_list.push_back(nrr);
            if (nlev < levmx) {
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[levmx-nlev-1];

               int ntrr = hash[(min( (jjj*2+1)*levmultminus, jmaxcalc-1))-jminsize][(min( (iii+1)  *levmult,      imaxcalc-1))-iminsize];
               if (ntrr != nrr && ntrr > 0 && (ntrr < start_idx || ntrr >= end_idx)) offtile_list.push_back(ntrr);
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
               printf("%d: Line %d Error with global cell match %d\n",mype, __LINE__, nb);
               L7_Terminate();
               exit(-1);
            }
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];
            //fprintf(fp,"%d: DEBUG nb is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nb, ig, border_cell_num_global[ig], nlev, iii, jjj);
            int nbb = hash[(max(  jjj   *levmult-1, 0)         )-jminsize][(      iii   *levmult               )-iminsize];
            if (nbb != nb && nbb > 0 && (nbb < start_idx || nbb >= end_idx) ) offtile_list.push_back(nbb);
            if (nlev < levmx) {
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[levmx-nlev-1];

               int nrbb = hash[(max(  jjj     *levmult-1,    0))         -jminsize][(min( (iii*2+1)*levmultminus, imaxcalc-1))-iminsize];
               if (nrbb != nbb && nrbb > 0 && (nrbb < start_idx || nrbb >= end_idx)) offtile_list.push_back(nrbb);
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
               L7_Terminate();
               exit(-1);
            }
            int nlev = border_cell_level_global[ig];
            int levmult      = levtable[levmx-nlev];
            int iii = border_cell_i_global[ig];
            int jjj = border_cell_j_global[ig];
            //fprintf(fp,"%d: DEBUG nt is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nt, ig, border_cell_num_global[ig], nlev, iii, jjj);
            int ntt = hash[(min( (jjj+1)*levmult,   jmaxcalc-1))-jminsize][(      iii   *levmult               )-iminsize];
            if (ntt != nt && ntt > 0 && (ntt < start_idx || ntt >= end_idx) ) offtile_list.push_back(ntt);
            if (nlev < levmx) {
               int levmult      = levtable[levmx-nlev];
               int levmultminus = levtable[levmx-nlev-1];

               int nrtt = hash[(min( (jjj+1)  *levmult,      jmaxcalc-1))-jminsize][(min( (iii*2+1)*levmultminus, imaxcalc-1))-iminsize];
               if (nrtt != ntt && nrtt > 0 && (nrtt < start_idx || nrtt >= end_idx)) offtile_list.push_back(nrtt);
            }
         }

         if (lev < levmx) {
            int levmult      = levtable[levmx-lev];
            int levmultminus = levtable[levmx-lev-1];

            int ntl = hash[(min( (jj*2+1)*levmultminus, jmaxcalc-1))-jminsize][(max(  ii     *levmult-1,    0))         -iminsize];
            if (ntl != nt && ntl > 0 && (ntl < start_idx || ntl >= end_idx)) offtile_list.push_back(ntl);

            int ntr = hash[(min( (jj*2+1)*levmultminus, jmaxcalc-1))-jminsize][(min( (ii+1)  *levmult,      imaxcalc-1))-iminsize];
            if (ntr != nt && ntr > 0 && (ntr < start_idx || ntr >= end_idx)) offtile_list.push_back(ntr);

            int nrb = hash[(max(  jj     *levmult-1,    0))         -jminsize][(min( (ii*2+1)*levmultminus, imaxcalc-1))-iminsize];
            if (nrb != nt && nrb > 0 && (nrb < start_idx || nrb >= end_idx)) offtile_list.push_back(nrb);

            int nrt = hash[(min( (jj+1)  *levmult,      jmaxcalc-1))-jminsize][(min( (ii*2+1)*levmultminus, imaxcalc-1))-iminsize];
            if (nrt != nt && nrt > 0 && (nrt < start_idx || nrt >= end_idx)) offtile_list.push_back(nrt);
         }
      }

      sort(offtile_list.begin(), offtile_list.end());
      p_end = unique(offtile_list.begin(), offtile_list.end());

      vector<int> indices_needed;

      for (vector<int>::iterator p=offtile_list.begin(); p < p_end; p++){
         indices_needed.push_back(*p);
      }

      int nghost = indices_needed.size();
      ncells_ghost = ncells + nghost;

      offtile_ratio_local = (offtile_ratio_local*(double)offtile_local_count) + ((double)nghost / (double)ncells);
      offtile_local_count++;
      offtile_ratio_local /= offtile_local_count;

      //fprintf(fp,"%d ncells_ghost size is %ld nghost %d\n",mype,ncells_ghost,nghost);

      if (cell_handle) L7_Free(&cell_handle);
      cell_handle=0;
      //for (int ig = 0; ig<nghost; ig++){
      //   fprintf(fp,"%d: indices_needed[%d]=%d\n",mype,ig,indices_needed[ig]);
      //}
      L7_Setup(0, noffset, ncells, &indices_needed[0], nghost, &cell_handle);

      //vector<int> itest(ncells_ghost);
      //for (int ic=0; ic<ncells; ic++){
      //   itest[ic] = mype*1000 + ic;
      //   fprintf(fp,"%d: test is filled with ic %d = %d\n",mype,ic,itest[ic]);
      //}
      //for (int ic=ncells; ic<ncells_ghost; ic++){
      //   itest[ic] = 0;
      //}

      //L7_Update(&itest[0], L7_INT, cell_handle);

      //for (int ic=0; ic<ncells_ghost; ic++){
      //   fprintf(fp,"%d: test after update ic %d = %d\n",mype,ic,itest[ic]);
      //}

      celltype.resize(ncells_ghost);
      i.resize(ncells_ghost);
      j.resize(ncells_ghost);
      level.resize(ncells_ghost);
      nlft.resize(ncells_ghost,-98);
      nrht.resize(ncells_ghost,-98);
      nbot.resize(ncells_ghost,-98);
      ntop.resize(ncells_ghost,-98);

      L7_Update(&celltype[0], L7_INT, cell_handle);
      L7_Update(&i[0],        L7_INT, cell_handle);
      L7_Update(&j[0],        L7_INT, cell_handle);
      L7_Update(&level[0],    L7_INT, cell_handle);
      //for (int ic=0; ic<ncells; ic++){
      //   fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
      //       mype,ic,i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
      //}
      L7_Update(&nlft[0], L7_INT, cell_handle);
      L7_Update(&nrht[0], L7_INT, cell_handle);
      L7_Update(&nbot[0], L7_INT, cell_handle);
      L7_Update(&ntop[0], L7_INT, cell_handle);
      //for (int ic=0; ic<ncells_ghost; ic++){
      //   fprintf(fp,"%d: 1655 nlft for %5d is %5d\n",mype,ic,nlft[ic]);
      //}

      for (uint ic=0; ic<ncells_ghost; ic++){
         if (nlft[ic] >= noffset && nlft[ic] < noffset+ncells) {
            nlft[ic] -= noffset;
            //fprintf(fp,"%d: 1: ic %d nlft is %d\n",mype,ic,nlft[ic]);
         } else {
            for (int ig=0; ig<nghost; ig++){
               //if (nlft[ic]==indices_needed[ig]) {nlft[ic] = ig+ncells; fprintf(fp,"%d: 2: ic %d nlft is %d\n",mype,ic,nlft[ic]); break;}
               if (nlft[ic]==indices_needed[ig]) {nlft[ic] = ig+ncells; break;}
            }
#ifdef XXX
            int nlev = level[nlft[ic]];
            int levmult = levtable[levmx-nlev];
            int iii = i[nlft[ic]];
            int jjj = j[nlft[ic]];
            int nll = hash[(      jjj   *levmult               )-jminsize][(max(  iii   *levmult-1, 0         ))-iminsize];
            fprintf(fp,"%d: 3: lev %d iii %d jjj %d nll is %d\n",mype,nlev,iii,jjj,nll);
            for (int ig=0; ig<nghost; ig++){
               if (nll==indices_needed[ig]) {nlft[nll] = (ig+ncells); fprintf(fp,"%d: 4: ic %d ig %d indices_needed %d nlft is %d\n",mype,ic,ig,indices_needed[ig],nlft[nll]); break;}
            }
#endif
         }
         if (nrht[ic] >= noffset && nrht[ic] < noffset+ncells) {
            nrht[ic] -= noffset;
         } else {
            for (int ig=0; ig<nghost; ig++){
               if (nrht[ic]==indices_needed[ig]) {nrht[ic] = ig+ncells; break;}
            }
         }
         if (nbot[ic] >= noffset && nbot[ic] < noffset+ncells) {
            nbot[ic] -= noffset;
         } else {
            for (int ig=0; ig<nghost; ig++){
               if (nbot[ic]==indices_needed[ig]) {nbot[ic] = ig+ncells; break;}
            }
         }
         if (ntop[ic] >= noffset && ntop[ic] < noffset+ncells) {
            ntop[ic] -= noffset;
         } else {
            for (int ig=0; ig<nghost; ig++){
               if (ntop[ic]==indices_needed[ig]) {ntop[ic] = ig+ncells; break;}
            }
         }
      }

      // Up to 
      //for (int ic=0; ic<ncells; ic++){
      //   if (nlft[ic] >= noffset && nlft[ic] < noffset+ncells) nlft[ic] -= noffset;
      //   if (nrht[ic] >= noffset && nrht[ic] < noffset+ncells) nrht[ic] -= noffset;
      //   if (nbot[ic] >= noffset && nbot[ic] < noffset+ncells) nbot[ic] -= noffset;
      //   if (ntop[ic] >= noffset && ntop[ic] < noffset+ncells) ntop[ic] -= noffset;
      //}

#ifdef XXX
      for (int ic=ncells; ic<ncells_ghost; ic++){
         if (nlft[ic] >= 0 && (nlft[ic] < noffset || nlft[ic] >= noffset+ncells) ) {
            //for (int ig=0; ig<nghost; ig++){
            //   if (nlft[ic]==indices_needed[ig]) {nlft[ic] = (ig+ncells); fprintf(fp,"%d: 2: ic %d nlft is %d\n",mype,ic,nlft[ic]); break;}
            //}
            int nlev = level[nlft[ic]];
            int levmult = levtable[levmx-nlev];
            int iii = i[nlft[ic]];
            int jjj = j[nlft[ic]];
            int nll = hash[(      jjj   *levmult               )-jminsize][(max(  iii   *levmult-1, 0         ))-iminsize];
            printf("%d: 3: lev %d iii %d jjj %d nll is %d\n",mype,nlev,iii,jjj,nll);
            for (int ig=0; ig<nghost; ig++){
               if (nll==indices_needed[ig]) {nlft[nlft[ic]] = (ig+ncells); fprintf(fp,"%d: 4: ic %d ig %d indices_needed %d nlft is %d\n",mype,ic,ig,indices_needed[ig],nlft[nll]); break;}
            }
         }
      }
#endif

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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
                  if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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

      genmatrixfree((void **)hash);

   } else if (calc_neighbor_type == KDTREE) {
      TBounds box;
      vector<int> index_list(20);

      int num;

      ibase = 0;
      calc_spatial_coordinates(ibase);

      kdtree_setup();

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

   } // calc_neighbor_type

   cpu_time_calc_neighbors += cpu_timer_stop(tstart_cpu);
}

void Mesh::gpu_calc_neighbors(cl_command_queue command_queue)
{
   cl_event hash_init_event;
   cl_event hash_setup_event;
   cl_event calc_neighbors_event;

   assert(dev_levtable);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);

   dev_nlft     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nrht     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nbot     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_ntop     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   size_t local_work_size = MIN(ncells, TILE_SIZE);
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;

   int imaxsize = (imax+1)*levtable[levmx];
   int jmaxsize = (jmax+1)*levtable[levmx];

   size_t hashsize = jmaxsize*imaxsize;
   cl_mem dev_hash = ezcl_malloc(NULL, &hashsize, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   size_t hash_local_work_size  = MIN(hashsize, TILE_SIZE);
   size_t hash_global_work_size = ((hashsize+hash_local_work_size - 1) /hash_local_work_size) * hash_local_work_size;

   ezcl_set_kernel_arg(kernel_hash_init, 0, sizeof(cl_int), (void *)&hashsize);
   ezcl_set_kernel_arg(kernel_hash_init, 1, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init,   1, NULL, &hash_global_work_size, &hash_local_work_size, &hash_init_event);

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
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup,   1, NULL, &global_work_size, &local_work_size, &hash_setup_event);

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
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors,   1, NULL, &global_work_size, &local_work_size, &calc_neighbors_event);

   ezcl_device_memory_remove(dev_hash);

   gpu_time_hash_setup        += ezcl_timer_calc(&hash_init_event,         &hash_init_event);
   gpu_time_hash_setup        += ezcl_timer_calc(&hash_setup_event,        &hash_setup_event);
   gpu_time_calc_neighbors    += ezcl_timer_calc(&calc_neighbors_event,    &calc_neighbors_event);
}


void Mesh::gpu_calc_neighbors_local(cl_command_queue command_queue)
{
   cl_event hash_size_event;
   //if (block_size > 1)
      cl_event finish_hash_size_event;
   cl_event start_read_event;
   cl_event hash_init_event;
   cl_event corners_init_event;
   cl_event hash_setup_event;
   cl_event calc_neighbors_init_event;
   //if (numpe > 1)
      cl_event hash_setup_border_event;
      cl_event calc_neighbors_event;
      cl_event copy_mesh_data_event;
      cl_event copy_ghost_data_event;                                                                                        
      cl_event adjust_neighbors_event;                                                                                        

   assert(dev_levtable);
   assert(dev_level);
   assert(dev_i);
   assert(dev_j);

   dev_nlft     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nrht     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_nbot     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   dev_ntop     = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);

   size_t local_work_size =  64;
   size_t global_work_size = ((ncells + local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   //printf("DEBUG file %s line %d lws = %d gws %d bs %d ncells %d\n",__FILE__,__LINE__,
   //   local_work_size, global_work_size, block_size, ncells);
   cl_mem dev_redscratch = ezcl_malloc(NULL, &block_size, sizeof(cl_int4), CL_MEM_READ_WRITE, 0);
   size_t one = 1;
   cl_mem dev_sizes = ezcl_malloc(NULL, &one, sizeof(cl_int4),  CL_MEM_READ_WRITE, 0);

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

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_size,   1, NULL, &global_work_size, &local_work_size, &hash_size_event);

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

   if (block_size > 1) {
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_finish_hash_size,   1, NULL, &local_work_size, &local_work_size, &finish_hash_size_event);
   }

   ezcl_device_memory_remove(dev_redscratch);

   cl_int sizes[4];
   ezcl_enqueue_read_buffer(command_queue, dev_sizes, CL_TRUE,  0, 1*sizeof(cl_int4), &sizes, &start_read_event);

   int iminsize = sizes[0];
   int imaxsize = sizes[1];
   int jminsize = sizes[2];
   int jmaxsize = sizes[3];

   //fprintf(fp,"%d: sizes %d %d     %d %d   \n",mype,sizes[0].s0,sizes[0].s1,sizes[0].s2,sizes[0].s3);
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
   cl_mem dev_hash = ezcl_malloc(NULL, &hashsize, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   size_t hash_local_work_size  = MIN(hashsize, TILE_SIZE);
   size_t hash_global_work_size = ((hashsize+hash_local_work_size - 1) /hash_local_work_size) * hash_local_work_size;

   //printf("%d: hash size is %d lws %d gws %d\n",mype,hashsize,hash_local_work_size,hash_global_work_size);
   ezcl_set_kernel_arg(kernel_hash_init, 0, sizeof(cl_int), (void *)&hashsize);
   ezcl_set_kernel_arg(kernel_hash_init, 1, sizeof(cl_mem), (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init,   1, NULL, &hash_global_work_size, &hash_local_work_size, &hash_init_event);

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
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init_corners,   1, NULL, &corners_global_work_size, &corners_local_work_size, &corners_init_event);

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
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup_local,   1, NULL, &global_work_size, &local_work_size, &hash_setup_event);

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
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors_local,   1, NULL, &global_work_size, &local_work_size, &calc_neighbors_init_event);

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
   ezcl_enqueue_read_buffer(command_queue, dev_level,    CL_FALSE, 0, ncells*sizeof(cl_int), &level_tmp[0],    NULL);

   ezcl_enqueue_read_buffer(command_queue, dev_nlft, CL_FALSE, 0, ncells*sizeof(cl_int), &nlft_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nrht, CL_FALSE, 0, ncells*sizeof(cl_int), &nrht_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_nbot, CL_FALSE, 0, ncells*sizeof(cl_int), &nbot_tmp[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_ntop, CL_TRUE,  0, ncells*sizeof(cl_int), &ntop_tmp[0], NULL);

   vector<int> border_cell;

   // Scan for corner boundary cells and also push list of unsatisfied neighbor cells
   for (uint ic=0; ic<ncells; ic++){
      if (nlft_tmp[ic] < 0){
         //printf("%d: Cell is %d nlft %d\n",mype,ic+noffset,nlft_tmp[ic]);
         border_cell.push_back(ic+noffset);
         if (nrht_tmp[ic] >= 0) {
            border_cell.push_back(nrht_tmp[ic]);
            if (level_tmp[nrht_tmp[ic]-noffset] > level_tmp[ic]) {
               if (ntop_tmp[nrht_tmp[ic]-noffset] >= 0) border_cell.push_back(ntop_tmp[nrht_tmp[ic]-noffset]);
            }
         }
      }
      if (nrht_tmp[ic] < 0){
         //printf("%d: Cell is %d nrht %d\n",mype,ic+noffset,nrht_tmp[ic]);
         border_cell.push_back(ic+noffset);
         if (nlft_tmp[ic] >= 0) {
            border_cell.push_back(nlft_tmp[ic]);
            if (level_tmp[nlft_tmp[ic]-noffset] > level_tmp[ic]) {
               if (ntop_tmp[nlft_tmp[ic]-noffset] >= 0) border_cell.push_back(ntop_tmp[nlft_tmp[ic]-noffset]);
            }
         }
      }
      if (nbot_tmp[ic] < 0) {
         //printf("%d: Cell is %d nbot %d\n",mype,ic+noffset,nbot_tmp[ic]);
         border_cell.push_back(ic+noffset);
         if (ntop_tmp[ic] >= 0) {
            border_cell.push_back(ntop_tmp[ic]);
            if (level_tmp[ntop_tmp[ic]-noffset] > level_tmp[ic]) {
               if (nrht_tmp[ntop_tmp[ic]-noffset] >= 0) border_cell.push_back(nrht_tmp[ntop_tmp[ic]-noffset]);
            }
         }
      }
      if (ntop_tmp[ic] < 0) {
         //printf("%d: Cell is %d ntop %d\n",mype,ic+noffset,ntop_tmp[ic]);
         border_cell.push_back(ic+noffset);
         if (nbot_tmp[ic] >= 0) {
            border_cell.push_back(nbot_tmp[ic]);
            if (level_tmp[nbot_tmp[ic]-noffset] > level_tmp[ic]) {
               if (nrht_tmp[nbot_tmp[ic]-noffset] >= 0) border_cell.push_back(nrht_tmp[nbot_tmp[ic]-noffset]);
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

   int nbsize_global;
   int nbsize_local=border_cell_num.size();

   MPI_Allreduce(&nbsize_local, &nbsize_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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
   MPI_Allgather(&nbsize_local, 1, MPI_INT, &nbsizes[0], 1, MPI_INT, MPI_COMM_WORLD);
   nbdispl[0]=0;
   for (int ip=1; ip<numpe; ip++){
      nbdispl[ip] += nbdispl[ip-1] + nbsizes[ip-1];
   }

   vector<int>border_cell_num_global(nbsize_global);
   vector<int>border_cell_i_global(nbsize_global);
   vector<int>border_cell_j_global(nbsize_global);
   vector<int>border_cell_level_global(nbsize_global);

   MPI_Allgatherv(&border_cell_num[0],   nbsizes[mype], MPI_INT, &border_cell_num_global[0],   &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&border_cell_i[0],     nbsizes[mype], MPI_INT, &border_cell_i_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&border_cell_j[0],     nbsizes[mype], MPI_INT, &border_cell_j_global[0],     &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);
   MPI_Allgatherv(&border_cell_level[0], nbsizes[mype], MPI_INT, &border_cell_level_global[0], &nbsizes[0], &nbdispl[0], MPI_INT, MPI_COMM_WORLD);

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
      if (border_cell_num_global[ic] >= noffset && border_cell_num_global[ic] < noffset+ncells) continue;
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

   size_t nbsize = nbsize_local;
   if (numpe > 1) {
      cl_mem dev_border_level = ezcl_malloc(&border_cell_level_global[0], &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      cl_mem dev_border_i     = ezcl_malloc(&border_cell_i_global[0],     &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      cl_mem dev_border_j     = ezcl_malloc(&border_cell_j_global[0],     &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);
      cl_mem dev_border_num   = ezcl_malloc(&border_cell_num_global[0],   &nbsize, sizeof(cl_int), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 0);

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
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_setup_border,  1, NULL, &border_global_work_size, &border_local_work_size, &hash_setup_border_event);

      ezcl_device_memory_remove(dev_border_level);
      ezcl_device_memory_remove(dev_border_i);
      ezcl_device_memory_remove(dev_border_j);
      ezcl_device_memory_remove(dev_border_num);

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
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_neighbors_local2,   1, NULL, &global_work_size, &local_work_size, &calc_neighbors_event);
   } // if (numpe > 1)

   ezcl_device_memory_remove(dev_sizes);

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
               if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
               if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
               if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
               if ( (ii >= iminsize && ii < imaxsize) && (hashval >= 0 && hashval < ncells) ) {
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
            L7_Terminate();
            exit(-1);
         }
         int nlev = border_cell_level_global[ig];
         int levmult      = levtable[levmx-nlev];
         int iii = border_cell_i_global[ig];
         int jjj = border_cell_j_global[ig];
         //fprintf(fp,"%d: DEBUG nl is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nb, ig, border_cell_global[ig], nlev, iii, jjj);
         int nll = hash_tmp[((      jjj   *levmult               )-jminsize)*(imaxsize-iminsize)+((max(  iii   *levmult-1, 0         ))-iminsize)];
         if (nll != nl && nll > 0 && (nll < start_idx || nll >= end_idx) ) offtile_list.push_back(nll);
         if (nlev < levmx) {
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[levmx-nlev-1];

            int ntll = hash_tmp[((min( (jjj*2+1)*levmultminus, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((max(  iii     *levmult-1,    0))         -iminsize)];
            if (ntll != nll && ntll > 0 && (ntll < start_idx || ntll >= end_idx)) offtile_list.push_back(ntll);
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
            L7_Terminate();
            exit(-1);
         }
         int nlev = border_cell_level_global[ig];
         int levmult      = levtable[levmx-nlev];
         int iii = border_cell_i_global[ig];
         int jjj = border_cell_j_global[ig];
         //fprintf(fp,"%d: DEBUG nr is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nb, ig, border_cell_global[ig], nlev, iii, jjj);
         int nrr = hash_tmp[((      jjj   *levmult               )-jminsize)*(imaxsize-iminsize)+((min( (iii+1)*levmult,   imaxcalc-1))-iminsize)];
         if (nrr != nr && nrr > 0 && (nrr < start_idx || nrr >= end_idx) ) offtile_list.push_back(nrr);
         if (nlev < levmx) {
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[levmx-nlev-1];

            int ntrr = hash_tmp[((min( (jjj*2+1)*levmultminus, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((min( (iii+1)  *levmult,      imaxcalc-1))-iminsize)];
            if (ntrr != nrr && ntrr > 0 && (ntrr < start_idx || ntrr >= end_idx)) offtile_list.push_back(ntrr);
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
            L7_Terminate();
            exit(-1);
         }
         int nlev = border_cell_level_global[ig];
         int levmult      = levtable[levmx-nlev];
         int iii = border_cell_i_global[ig];
         int jjj = border_cell_j_global[ig];
         //fprintf(fp,"%d: DEBUG nb is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nb, ig, border_cell_global[ig], nlev, iii, jjj);
         int nbb = hash_tmp[((max(  jjj   *levmult-1, 0)         )-jminsize)*(imaxsize-iminsize)+((      iii   *levmult               )-iminsize)];
         if (nbb != nb && nbb > 0 && (nbb < start_idx || nbb >= end_idx) ) offtile_list.push_back(nbb);
         if (nlev < levmx) {
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[levmx-nlev-1];

            int nrbb = hash_tmp[((max(  jjj     *levmult-1,    0))         -jminsize)*(imaxsize-iminsize)+((min( (iii*2+1)*levmultminus, imaxcalc-1))-iminsize)];
            if (nrbb != nbb && nrbb > 0 && (nrbb < start_idx || nrbb >= end_idx)) offtile_list.push_back(nrbb);
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
            L7_Terminate();
            exit(-1);
         }
         int nlev = border_cell_level_global[ig];
         int levmult      = levtable[levmx-nlev];
         int iii = border_cell_i_global[ig];
         int jjj = border_cell_j_global[ig];
         //fprintf(fp,"%d: DEBUG nt is %d ig is %d cell num is %d lev %d i %d j %d\n",mype,nt, ig, border_cell_num_global[ig], nlev, iii, jjj);
         int ntt = hash_tmp[((min( (jjj+1)*levmult,   jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((      iii   *levmult               )-iminsize)];
         if (ntt != nt && ntt > 0 && (ntt < start_idx || ntt >= end_idx) ) offtile_list.push_back(ntt);
         if (nlev < levmx) {
            int levmult      = levtable[levmx-nlev];
            int levmultminus = levtable[levmx-nlev-1];

            int nrtt = hash_tmp[((min( (jjj+1)  *levmult,      jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((min( (iii*2+1)*levmultminus, imaxcalc-1))-iminsize)];
            if (nrtt != ntt && nrtt > 0 && (nrtt < start_idx || nrtt >= end_idx)) offtile_list.push_back(nrtt);
         }
      }

      if (lev < levmx) {
         int levmult      = levtable[levmx-lev];
         int levmultminus = levtable[levmx-lev-1];

         int ntl = hash_tmp[((min( (jj*2+1)*levmultminus, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((max(  ii     *levmult-1,    0))         -iminsize)];
         if (ntl != nt && ntl > 0 && (ntl < start_idx || ntl >= end_idx)) offtile_list.push_back(ntl);

         int ntr = hash_tmp[((min( (jj*2+1)*levmultminus, jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((min( (ii+1)  *levmult,      imaxcalc-1))-iminsize)];
         if (ntr != nt && ntr > 0 && (ntr < start_idx || ntr >= end_idx)) offtile_list.push_back(ntr);

         int nrb = hash_tmp[((max(  jj     *levmult-1,    0))         -jminsize)*(imaxsize-iminsize)+((min( (ii*2+1)*levmultminus, imaxcalc-1))-iminsize)];
         if (nrb != nt && nrb > 0 && (nrb < start_idx || nrb >= end_idx)) offtile_list.push_back(nrb);

         int nrt = hash_tmp[((min( (jj+1)  *levmult,      jmaxcalc-1))-jminsize)*(imaxsize-iminsize)+((min( (ii*2+1)*levmultminus, imaxcalc-1))-iminsize)];
         if (nrt != nt && nrt > 0 && (nrt < start_idx || nrt >= end_idx)) offtile_list.push_back(nrt);
      }
   }

   sort(offtile_list.begin(), offtile_list.end());
   p_end = unique(offtile_list.begin(), offtile_list.end());

   vector<int> indices_needed;

   for (vector<int>::iterator p=offtile_list.begin(); p < p_end; p++){
      indices_needed.push_back(*p);
   }

   int nghost = indices_needed.size();
   ncells_ghost = ncells + nghost;
   //fprintf(fp,"%d ncells_ghost size is %ld nghost %d\n",mype,ncells_ghost,nghost);

   if (cell_handle) L7_Free(&cell_handle);
   cell_handle=0;
   //for (int ic=0; ic<nghost; ic++){
   //   fprintf(fp,"%d: indices needed ic %d index %d\n",mype,ic,indices_needed[ic]);
   //}
   L7_Setup(0, noffset, ncells, &indices_needed[0], nghost, &cell_handle);

   celltype_tmp.resize(ncells_ghost);
   i_tmp.resize(ncells_ghost);
   j_tmp.resize(ncells_ghost);
   level_tmp.resize(ncells_ghost);
   nlft_tmp.resize(ncells_ghost,-98);
   nrht_tmp.resize(ncells_ghost,-98);
   nbot_tmp.resize(ncells_ghost,-98);
   ntop_tmp.resize(ncells_ghost,-98);

   L7_Update(&celltype_tmp[0], L7_INT, cell_handle);
   L7_Update(&i_tmp[0],        L7_INT, cell_handle);
   L7_Update(&j_tmp[0],        L7_INT, cell_handle);
   L7_Update(&level_tmp[0],    L7_INT, cell_handle);

   //for (int ic=0; ic<ncells; ic++){
   //   fprintf(fp,"%d: before update ic %d        i %d j %d lev %d nlft %d nrht %d nbot %d ntop %d\n",
   //       mype,ic,i[ic],j[ic],level[ic],nlft[ic],nrht[ic],nbot[ic],ntop[ic]);
   //}

   L7_Update(&nlft_tmp[0], L7_INT, cell_handle);
   L7_Update(&nrht_tmp[0], L7_INT, cell_handle);
   L7_Update(&nbot_tmp[0], L7_INT, cell_handle);
   L7_Update(&ntop_tmp[0], L7_INT, cell_handle);

      //vector<int> itest(ncells_ghost);
      //for (int ic=0; ic<ncells; ic++){
      //   itest[ic] = mype*1000 + ic;
      //   fprintf(fp,"%d: test is filled with ic %d = %d\n",mype,ic,itest[ic]);
      //}
      //for (int ic=ncells; ic<ncells_ghost; ic++){
      //   itest[ic] = 0;
      //}

      //L7_Update(&itest[0], L7_INT, cell_handle);

      //for (int ic=0; ic<ncells_ghost; ic++){
      //   fprintf(fp,"%d: test after update ic %d = %d\n",mype,ic,itest[ic]);
      //}

   if (numpe > 1) {
      cl_mem dev_celltype_old = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_i_old        = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_j_old        = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_level_old    = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nlft_old     = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nrht_old     = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nbot_old     = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_ntop_old     = ezcl_malloc(NULL, &ncells_ghost, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

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

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_mesh_data,   1, NULL, &global_work_size, &local_work_size, &copy_mesh_data_event);
   
      ezcl_device_memory_remove(dev_celltype_old);
      ezcl_device_memory_remove(dev_i_old);
      ezcl_device_memory_remove(dev_j_old);
      ezcl_device_memory_remove(dev_level_old);
      ezcl_device_memory_remove(dev_nlft_old);
      ezcl_device_memory_remove(dev_nrht_old);
      ezcl_device_memory_remove(dev_nbot_old);
      ezcl_device_memory_remove(dev_ntop_old);

      size_t nghost_local = nghost;
      cl_mem dev_nlft_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nrht_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_nbot_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_ntop_add       = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_celltype_add   = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_i_add          = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_j_add          = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_level_add      = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
      cl_mem dev_indices_needed = ezcl_malloc(NULL, &nghost_local,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
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

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_ghost_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, &copy_ghost_data_event);                                                                                                      
      ezcl_device_memory_remove(dev_nlft_add);
      ezcl_device_memory_remove(dev_nrht_add);
      ezcl_device_memory_remove(dev_nbot_add);
      ezcl_device_memory_remove(dev_ntop_add);

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

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_adjust_neighbors,   1, NULL, &nc_ghost_global_work_size, &nc_ghost_local_work_size, &adjust_neighbors_event);
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

   ezcl_device_memory_remove(dev_hash);

   gpu_time_hash_setup     += ezcl_timer_calc(&hash_size_event,         &hash_size_event);
   if (block_size > 1){
      gpu_time_hash_setup     += ezcl_timer_calc(&finish_hash_size_event,  &finish_hash_size_event);
   }
   gpu_time_hash_setup     += ezcl_timer_calc(&start_read_event,        &start_read_event);
   gpu_time_hash_setup     += ezcl_timer_calc(&hash_init_event,         &hash_init_event);
   gpu_time_hash_setup     += ezcl_timer_calc(&corners_init_event,      &corners_init_event);
   gpu_time_hash_setup     += ezcl_timer_calc(&hash_setup_event,        &hash_setup_event);
   gpu_time_hash_setup     += ezcl_timer_calc(&calc_neighbors_init_event, &calc_neighbors_init_event);
   if (numpe > 1) {
      gpu_time_hash_setup     += ezcl_timer_calc(&hash_setup_border_event, &hash_setup_border_event);
      gpu_time_calc_neighbors += ezcl_timer_calc(&calc_neighbors_event,    &calc_neighbors_event);
      gpu_time_hash_setup     += ezcl_timer_calc(&copy_mesh_data_event,    &copy_mesh_data_event);
      gpu_time_hash_setup     += ezcl_timer_calc(&copy_ghost_data_event,   &copy_ghost_data_event);
      gpu_time_hash_setup     += ezcl_timer_calc(&adjust_neighbors_event,  &adjust_neighbors_event);
   }

}

void Mesh::print_calc_neighbor_type(void)
{
   if ( calc_neighbor_type == HASH_TABLE ) {
      if (mype == 0) printf("Using hash tables to calculate neighbors\n");
   } else {
      if (mype == 0) printf("Using k-D tree to calculate neighbors\n");
   }
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
   vector<int> index_list(20);

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
void Mesh::resize_old_device_memory(size_t ncells)
{
   ezcl_device_memory_remove(dev_level);
   ezcl_device_memory_remove(dev_i);
   ezcl_device_memory_remove(dev_j);
   ezcl_device_memory_remove(dev_celltype);
   dev_level    = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_i        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_j        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_celltype = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
}
void Mesh::resize_new_device_memory(size_t ncells)
{
   ezcl_device_memory_remove(dev_level_new);
   ezcl_device_memory_remove(dev_i_new);
   ezcl_device_memory_remove(dev_j_new);
   ezcl_device_memory_remove(dev_celltype_new);
   dev_level_new    = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_i_new        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_j_new        = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
   dev_celltype_new = ezcl_malloc(NULL, &ncells, sizeof(cl_int),  CL_MEM_READ_ONLY, 0);
}
void Mesh::swap_device_memory_ptrs(void)
{
   cl_mem dev_ptr;
   SWAP_PTR(dev_level_new,    dev_level,    dev_ptr);
   SWAP_PTR(dev_i_new,        dev_i,        dev_ptr);
   SWAP_PTR(dev_j_new,        dev_j,        dev_ptr);
   SWAP_PTR(dev_celltype_new, dev_celltype, dev_ptr);
}
