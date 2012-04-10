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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <list>
#include "partition.h"
#include "kdtree/KDTree.h"
#include "mesh.h"
#include "s7/s7.h"
#include "mpi.h"
#include "zorder/zorder.h"
#include "timer/timer.h"

int measure_type;
int      meas_count                  = 0;
double   meas_sum_average            = 0.0;

extern bool localStencil;
extern enum partition_method initial_order;
extern enum partition_method cycle_reorder;

void Mesh::partition_measure(void) 
{
  int ntX     = TILE_SIZE; 
  double offtile_ratio = 0.0;

  int num_groups = (ncells + TILE_SIZE - 1)/TILE_SIZE;

  if (measure_type == WITH_DUPLICATES) {
     for (int group_id=0, i = 0; group_id < num_groups; group_id ++){ 
 
        int start_idx = group_id * ntX;
        int end_idx   = (group_id + 1) * ntX; 

        int offtile =0;
        for (int ic = 0; ic < TILE_SIZE; ic++, i++){ 

           if (i >= ncells) continue;
           //taken from wave_kern_calc.cl 'setup tile' kernel
           if (nlft[i] < start_idx || nlft[i] >= end_idx) offtile++; 
           if (level[nlft[i]] > level[i] &&
              (ntop[nlft[i]] < start_idx || ntop[nlft[i]] >= end_idx) ) offtile++;
           if (nrht[i] < start_idx || nrht[i] >= end_idx) offtile++; 
           if (level[nrht[i]] > level[i] &&
              (ntop[nrht[i]] < start_idx || ntop[nrht[i]] >= end_idx) ) offtile++;
           if (nbot[i] < start_idx || nbot[i] >= end_idx) offtile++; 
           if (level[nbot[i]] > level[i] &&
              (nrht[nbot[i]] < start_idx || nrht[nbot[i]] >= end_idx) ) offtile++;
           if (ntop[i] < start_idx || ntop[i] >= end_idx) offtile++; 
           if (level[ntop[i]] > level[i] &&
              (nrht[ntop[i]] < start_idx || nrht[ntop[i]] >= end_idx) ) offtile++;
        }
        offtile_ratio += (double)offtile/(double)(TILE_SIZE);
        //printf("DEBUG Ratio of surface area to volume is equal to %d / %d ratio is %lf\n", offtile, TILE_SIZE, (double)offtile/(double)TILE_SIZE);
     }
  } else if (measure_type == WITHOUT_DUPLICATES) {

     for (int group_id=0, i = 0; group_id < num_groups; group_id ++){ 
        list<int> offtile_list;
 
        int start_idx = group_id * ntX;
        int end_idx   = (group_id + 1) * ntX; 

        for (int ic = 0; ic < TILE_SIZE; ic++, i++){ 

           if (i >= ncells) continue;

           if (nlft[i] < start_idx || nlft[i] >= end_idx) offtile_list.push_back(nlft[i]);
           if (level[nlft[i]] > level[i] &&
              (ntop[nlft[i]] < start_idx || ntop[nlft[i]] >= end_idx) ) offtile_list.push_back(ntop[nlft[i]]);
           if (nrht[i] < start_idx || nrht[i] >= end_idx) offtile_list.push_back(nrht[i]);
           if (level[nrht[i]] > level[i] &&
              (ntop[nrht[i]] < start_idx || ntop[nrht[i]] >= end_idx) ) offtile_list.push_back(ntop[nrht[i]]);
           if (nbot[i] < start_idx || nbot[i] >= end_idx) offtile_list.push_back(nbot[i]);
           if (level[nbot[i]] > level[i] &&
              (nrht[nbot[i]] < start_idx || nrht[nbot[i]] >= end_idx) ) offtile_list.push_back(nrht[nbot[i]]);
           if (ntop[i] < start_idx || ntop[i] >= end_idx) offtile_list.push_back(ntop[i]);
           if (level[ntop[i]] > level[i] &&
              (nrht[ntop[i]] < start_idx || nrht[ntop[i]] >= end_idx) ) offtile_list.push_back(nrht[ntop[i]]);
        }
        offtile_list.sort();
        offtile_list.unique();
        
        offtile_ratio += (double)offtile_list.size()/(double)(TILE_SIZE);
        //printf("DEBUG Ratio of surface area to volume is equal to %d / %d ratio is %lf\n", offtile, TILE_SIZE, (double)offtile/(double)TILE_SIZE);
     }
  } else if (measure_type == CVALUE) {

     for (int group_id=0, i = 0; group_id < num_groups; group_id ++){ 
        list<int> offtile_list;
 
        int start_idx = group_id * ntX;
        int end_idx   = (group_id + 1) * ntX; 

        for (int ic = 0; ic < TILE_SIZE; ic++, i++){ 

           if (i >= ncells) continue;

           if (nlft[i] < start_idx || nlft[i] >= end_idx) offtile_list.push_back(nlft[i]);
           if (level[nlft[i]] > level[i] &&
              (ntop[nlft[i]] < start_idx || ntop[nlft[i]] >= end_idx) ) offtile_list.push_back(ntop[nlft[i]]);
           if (nrht[i] < start_idx || nrht[i] >= end_idx) offtile_list.push_back(nrht[i]);
           if (level[nrht[i]] > level[i] &&
              (ntop[nrht[i]] < start_idx || ntop[nrht[i]] >= end_idx) ) offtile_list.push_back(ntop[nrht[i]]);
           if (nbot[i] < start_idx || nbot[i] >= end_idx) offtile_list.push_back(nbot[i]);
           if (level[nbot[i]] > level[i] &&
              (nrht[nbot[i]] < start_idx || nrht[nbot[i]] >= end_idx) ) offtile_list.push_back(nrht[nbot[i]]);
           if (ntop[i] < start_idx || ntop[i] >= end_idx) offtile_list.push_back(ntop[i]);
           if (level[ntop[i]] > level[i] &&
              (nrht[ntop[i]] < start_idx || nrht[ntop[i]] >= end_idx) ) offtile_list.push_back(nrht[ntop[i]]);
        }
        offtile_list.sort();
        offtile_list.unique();
        
        offtile_ratio += (double)offtile_list.size()/(4*sqrt((double)(TILE_SIZE)));
        //printf("DEBUG Ratio of surface area to volume is equal to %d / %d ratio is %lf\n", offtile, TILE_SIZE, (double)offtile/(double)TILE_SIZE);
     }
  }
  // printf("DEBUG Ratio of surface area to volume is equal to %d / %d \n", offtile, ontile);
   
   meas_count ++;
   meas_sum_average  += offtile_ratio/(double)num_groups;
  // printf("DEBUG %d icount %d sum_average %lf\n",__LINE__,icount, sum_average);

}
void Mesh::print_partition_measure()
{
      
   if (parallel) {
      vector<double> global_times(numpe);
      double local_time;

      if (measure_type == WITH_DUPLICATES) {
         local_time = meas_sum_average/(double)meas_count;
         MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_times[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         if (mype == 0) {
            printf("Average surface area to volume ratio  \t");  
            for(int ip = 0; ip < numpe; ip++){
               printf("%8.4f\t", global_times[ip]);
            }
            printf("with duplicates\n");
         }
      } else if (measure_type == WITHOUT_DUPLICATES) {
         local_time = meas_sum_average/(double)meas_count;
         MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_times[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         if (mype == 0) {
            printf("Average surface area to volume ratio  \t");  
            for(int ip = 0; ip < numpe; ip++){
               printf("%8.4f\t", global_times[ip]);
            }
            printf("without duplicates\n");
         }
      } else if (measure_type == CVALUE) {
         local_time = meas_sum_average/(double)meas_count;
         MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_times[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         if (mype == 0) {
            printf("The GPU Partition Quality Avg C value  \t");  
            for(int ip = 0; ip < numpe; ip++){
               printf("%8.4f\t", global_times[ip]);
            }
            printf("\n");
         }
      }

      local_time = offtile_ratio_local;
      MPI_Gather(&local_time, 1, MPI_DOUBLE, &global_times[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (mype == 0) {
         printf("The MPI surface area to volume ratio \t");
         for(int ip = 0; ip < numpe; ip++){
            printf("%8.4f\t", global_times[ip]);
         }
         printf("without duplicates\n");
      }

   } else {

      if (measure_type == WITH_DUPLICATES) {
         printf("Average surface area to volume ratio  \t%8.4lf\t with duplicates\n" , meas_sum_average/(double)meas_count);
      } else if (measure_type == WITHOUT_DUPLICATES) {
         printf("Average surface area to volume ratio  \t%8.4lf\t without duplicates\n" , meas_sum_average/(double)meas_count);
      } else if (measure_type == CVALUE) {
         printf("The GPU Partition Quality Avg C value  \t%8.4lf\n" , meas_sum_average/(double)meas_count);
      }

      printf("The MPI surface area to volume ratio \t%8.4lf\t without duplicates\n", offtile_ratio_local);

   }
}

void Mesh::print_partition_type()
{
   if (mype == 0) {
      if (initial_order == ORIGINAL_ORDER) {
         printf("Initial order is naive.");  
      } else if (initial_order == HILBERT_SORT) {
         printf("Initial order is Hilbert sort.");  
      } else if (initial_order == HILBERT_PARTITION) {
         printf("Initial order is Hilbert partitionr.");  
      } else if (initial_order == ZORDER) {
         printf("Initial order is Z order.");  
      }

      if (cycle_reorder == ORIGINAL_ORDER) {
         printf("   No cycle reorder.");  
      } else if (cycle_reorder == HILBERT_SORT) {
         printf("   Cycle reorder is Hilbert sort.");  
      } else if (cycle_reorder == HILBERT_PARTITION) {
         printf("   Cycle reorder is Hilbert partition.");  
      } else if (cycle_reorder == ZORDER) {
         printf("   Cycle reorder is Z order.");  
      }

      if (localStencil) {
         printf("   Local Stencil is on.\n");  
      } else {
         printf("\n");
      }
   }

}
void Mesh::partition_cells(
                    int          numpe,             //  
                    vector<int> &proc,              //  Work units by assigned work group.
                    vector<int> &z_order,           //  Resulting index ordering.
                    enum partition_method method)   //  Assigned partitioning method.
{  
   int           *info;      //
   double         iscale,    //
                  jscale,    //
                  xlocdiff,  //
                  ylocdiff;  //
   int            imax,      //  Maximum x-index.
                  jmax,      //  Maximum y-index.
                  lev,       //
                  limit;     //
   vector<int>    z_index;   //  Ordered curve from hsfc.
   vector<int>    i_scaled;  //  x-indices normalized to a scale of [0, 1] for hsfc.
   vector<int>    j_scaled;  //  y-indices normalized to a scale of [0, 1] for hsfc.
   vector<double> iunit;     //
   vector<double> junit;     //

   struct timeval tstart_cpu;

   cpu_timer_start(&tstart_cpu);

   //  Initialize ordered curve index.
   z_index.resize(ncells, 0);
   z_order.resize(ncells, 0);
   
   //  Partition cells according to one of several possible orderings.
   switch (method)
   {   case ORIGINAL_ORDER:
         //  Set z_order to the current cell order.
         for (int ic = 0; ic < ncells; ++ic)
         {   z_order[ic] = ic; }

         cpu_time_partition += cpu_timer_stop(tstart_cpu);

         return;
         break;

       case HILBERT_SORT:
         //  Resort the curve by Hilbert order.
         calc_centerminmax();
         iunit.resize(ncells);
         junit.resize(ncells);

         //   Get the range of values in the x- and y-directions and make the scale square.
         iscale = 1.0 / (xcentermax - xcentermin);
         jscale = 1.0 / (ycentermax - ycentermin);

         // XXX NOT NEEDED XXX //
//         if (iscale > jscale) iscale = jscale;
//         if (jscale > iscale) jscale = iscale;

         //   Scale the indices to a normalized [0, 1] range for hsfc.
         for (int ic = 0; ic < ncells; ++ic)
         {   iunit[ic] = (x[ic] + 0.5 * dx[ic] - xcentermin) * iscale;
             junit[ic] = (y[ic] + 0.5 * dy[ic] - ycentermin) * jscale; }
         info = (int *)malloc(sizeof(int) * 3 * ncells);

         //   Sort the mesh into an ordered space-filling curve from hsfc.
         hsfc2sort(ncells, &(iunit[0]), &(junit[0]), 0, info, 1);

         //   Copy the cell order information from info into z_order.
         for (int ic = 0; ic < ncells; ++ic)
         {   z_order[ic] = info[ic]; }
         free(info);

         break;

      case ZORDER:
         //  Resort the curve by z-order.
         i_scaled.resize(ncells);
         j_scaled.resize(ncells);

         //
         imax = 0;
         jmax = 0;
         for (int ic = 0; ic < ncells; ++ic)
         {   if (i[ic] > imax) imax = i[ic];
            if (j[ic] > jmax) jmax = j[ic]; }

         //
         iscale = 16.0 / (double)imax;
         jscale = 16.0 / (double)jmax;

         //
         for (int ic = 0; ic < ncells; ++ic)
         {   i_scaled[ic]=(int) ( (double)i[ic]*iscale);
            j_scaled[ic]=(int) ( (double)j[ic]*jscale); }

         //
         calc_zorder(ncells, &(i_scaled[0]), &(j_scaled[0]), &(level[0]), levmx, ibase, &(z_index[0]), &(z_order[0]));

         break;

      default:
         //  Note that HILBERT_PARTITION is not currently supported due to redundancy with HILBERT_SORT.
         break;
   }
   
   //   Adjust the number of required work items to the number of cells.
   proc.resize(ncells);
   
   //   Order the mesh according to the calculated order (note that z_order is for both curves).
   mesh_reorder(z_order);
   
   //   Output ordered mesh information.
   if (DEBUG)
   {   printf("orig index   i     j     lev    nlft nrht nbot ntop   xlow    xhigh     ylow    yhigh   z index  z order\n");
      for (int ic=0; ic<ncells; ic++){
         printf(" %6d   %4d  %4d   %4d  %4d %4d %4d %4d ", index[ic], j[ic], i[ic], level[ic], nlft[ic], nrht[ic], nbot[ic], ntop[ic]);
         printf(" %8.2lf %8.2lf %8.2lf %8.2lf", x[ic], x[ic]+dx[ic], y[ic], y[ic]+dy[ic]);
         printf(" %6d    %5d\n", z_index[ic], z_order[ic]); } }

   //   Decompose the domain equitably.
   calc_distribution(numpe, proc);

   cpu_time_partition += cpu_timer_stop(tstart_cpu);
}

//   The distribution needs to be modified in order to spread out extra cells equitably among the work items.
void Mesh::calc_distribution(int numpe,  //  Number of work items between which the domain is to be divided.
                             vector<int> &proc)   //   List of work items.
{  
   int lsize = 0;     //
   int ic    = 0;            //   Overall work item index.
   for (int ip = 0; ip < numpe; ++ip) {
      lsize += proc.size()/numpe;
      if (ip < proc.size()%numpe) lsize++;
      for (; ic < lsize;) {
         proc[ic] = ip;
         ic++;
      }
   }
}

