#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "l7.h"
#include "l7p.h"

void report_results_update(double *time_total_pe, int count_updated_pe, int num_timings,
      int num_timings_cycle);


void update_test()
{
   int penum,
       numpes;
   
   int i, j, num_indices_owned, num_updates,
     my_start_index, max_num_partners, remainder,
     num_partners_lo, num_partners_hi, num_partners,
     offset, num_indices_offpe, num_indices_per_partner,
     inum, l7_id, gtime, count_updated_pe, num_timings_cycle,
     num_timings, iout;
   
   double time_start, time_stop;
   double *time_total_pe;
   
   int *partner_pe;
   int *needed_indices;
   
   int *idata;
   double *rdata;
   
   int num_indices_per_pe = 10;
   int num_iterations = 10;
   int num_updates_per_cycle = 2;
   
   penum  = L7_Get_Rank();
   numpes = L7_Get_Numpes();

   if (penum == 0)
      printf("\n\t\tRunning the Update tests\n\n");
   
   num_updates = num_iterations * num_updates_per_cycle;
 
   time_total_pe = (double *)malloc((num_updates+1) * sizeof(double));
   
   num_indices_owned = num_indices_per_pe;
   my_start_index = penum * num_indices_owned;
   
   if (numpes > 16) {
      max_num_partners = (int)sqrt( (double)numpes);
   }
   else {
      max_num_partners = numpes/2;
   }
   
   /*
    * Create and load needed_indices array
    */
   
   remainder = max_num_partners % 2;
   
   if (penum < (numpes /2) ) {
      num_partners_lo = max_num_partners / 2;
      num_partners_hi = max_num_partners / 2 + remainder;
   }
   else {
      num_partners_lo = max_num_partners / 2 + remainder;
      num_partners_hi = max_num_partners / 2;
   }
   
   for (;penum - num_partners_lo < 0;num_partners_lo--);

   for (;penum + num_partners_hi >= numpes; num_partners_hi--);

   num_partners = num_partners_lo + num_partners_hi;
   partner_pe = (int *)malloc(num_partners * sizeof(int));
   
   offset = 0;
   for (i=1; i<=num_partners_lo; i++){
      partner_pe[offset] = penum - i;
      offset++;
   }

   for (i=1; i<=num_partners_hi; i++){
      //printf("[1pe %d] offset %d penum %d i %d \n",penum, offset, penum, i);
      partner_pe[offset] = penum + i;
      offset++;
   }
   
   if (num_partners != 0) {
      num_indices_offpe = num_indices_owned / 2;
      num_indices_per_partner = num_indices_offpe / num_partners;
   }
   else {
      num_indices_offpe = 0;
      num_indices_per_partner = 0;
   }

   needed_indices = (int *)malloc(num_indices_offpe * sizeof(int));
   
   /*
    * Generate needed indices
    */
   
   num_indices_offpe = 0;
   for (i=0; i<num_partners; i++) {
      inum = partner_pe[i] * num_indices_owned;
      for (j=0; j<num_indices_per_partner; j++){
         needed_indices[num_indices_offpe] = inum;
         num_indices_offpe++;
         inum+=2;
      }
   }
   
   
#ifdef _L7_DEBUG
   
   for (i=0; i<=numpes-1; i++){
      
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      if (penum == i) {
         for (j=0; j<num_indices_offpe; j++ ) {
            printf("[pe %1d] needed indices(%3d) = %d\n",
                  penum, j, needed_indices[j]);
         }
         fflush(stdout);
      }
   }
   
   if (penum == 0) printf("\n");
#endif
   
   /*
    * Allocate data arrays
    */
   
   rdata = (double *)malloc((num_indices_owned + num_indices_offpe)*
         sizeof(double));
   
   inum = my_start_index;
   for (i=0; i<num_indices_owned; i++){
      rdata[i] = (double) inum;
      inum++;
   }
   
   idata = (int *)malloc((num_indices_owned + num_indices_offpe)*
         sizeof(int));
   
   inum = my_start_index;
   for (i=0; i<num_indices_owned; i++){
      idata[i] = inum;
      inum++;
   }
   
   /*
    * Register decomposition with L7
    */
   
   l7_id = 0;
   
#ifdef _L7_DEBUG
   
   for (i=0; i<numpes; i++){
      
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      if (penum == i) {
         printf("[pe %1d] my_start_index, num_indices_owned, num_indices_offpe = (%3d, %3d %3d)",
               penum, my_start_index, num_indices_owned, num_indices_offpe);
         fflush(stdout);
      }
   }
   
   if (penum == 0) printf("\n");
   
#endif
   
#ifdef _L7_DEBUG
   if (penum == 0) {
      
      printf("\n");
      printf("======================================\n");
      printf("    Begin L7_Update testing\n");
      printf("======================================\n");
      printf("\n");
      
      printf("+++ call l7_setup +++\n");
      printf("\n");
      
   }
#endif
   
#ifdef USE_MPI
   //MPI_Barrier(MPI_COMM_WORLD);
#endif
   
   time_start=L7_Wtime();
   
   /*
    * Register decomposition with L7
    */
   
   L7_Setup(0, my_start_index, num_indices_owned, needed_indices, 
       num_indices_offpe, &l7_id);
   
   time_stop = L7_Wtime();
   time_total_pe[0] = time_stop - time_start;
   
   /*
    * Begin updating data
    */
   
   gtime = 1;
   for (i=1; i<=num_iterations; i++){

#ifdef _L7_DEBUG
      if (penum == 0) {
         printf(" ** Begin cycle %3d on l7 id %3d **\n",i,l7_id);
         printf("\n");
         printf("   Gather integer data (idata)\n");
      }
#endif
      
#ifdef USE_MPI
      //MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      time_start = L7_Wtime();
      
      L7_Update(idata, L7_INT, l7_id);


      time_stop = L7_Wtime();
      time_total_pe[gtime] = time_stop - time_start;
      gtime++;
      
#ifdef _L7_DEBUG
      if (penum == 0) {
         printf("   Gather real data (rdata)\n");
      }
#endif
      
#ifdef USE_MPI
      //MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      time_start = L7_Wtime();
      
      L7_Update(rdata, L7_DOUBLE, l7_id);
      
      time_stop = L7_Wtime();
      time_total_pe[gtime] = time_stop - time_start;
      gtime++;
      
#ifdef _L7_DEBUG
      if (penum == 0) {
         printf("\n");
         printf(" ** End cycle %3d on l7 id %3d **\n",i, l7_id);
         printf("\n");
      }
#endif
   }
   
   L7_Free(&l7_id);

   /*
    * Report results
    */
   
   count_updated_pe = num_indices_offpe;
   
   num_timings_cycle = num_updates_per_cycle;
   num_timings = num_updates +1;
   
#ifdef _L7_DEBUG
   report_results_update(time_total_pe, count_updated_pe, 
         num_timings, num_timings_cycle);
#endif
   
   /*
    * Testing complete
    */
   
   iout = 0;
   if (penum == 0) {
       if (iout > 0){
         printf("  Error with L7_Update with int and double arrays\n");
       }
       else{
         printf("  PASSED L7_Update with int and double arrays\n");
       }
   }

   free(time_total_pe);
   free(partner_pe);
   free(needed_indices);
   free(idata);
   free(rdata);
   
   return;
   
}

void report_results_update(double *time_total_pe, int count_updated_pe, int num_timings,
      int num_timings_cycle)
{
   
   int i, count_updated_global, penum, bytes_updated, remainder;
   
   double *time_total_global;
   
   
   time_total_global = (double *)malloc(num_timings*sizeof(double));
   
   L7_Array_Max(time_total_pe, num_timings, L7_DOUBLE, time_total_global);

   L7_Sum(&count_updated_pe, 1, L7_INT, &count_updated_global);

   penum = l7.penum;

   if (penum == 0){
      printf("\n");
      printf("======================================\n");
      printf("    L7_Update test results            \n");
      printf("======================================\n");
      printf("\n");
      printf("Performance\n");
      printf("\n");
      printf("    L7_Setup:    %lf seconds\n",time_total_global[0]);
      printf("\n");
      printf("      L7_Update:   time (secs)      bandwidth    cycle\n");
      
      for (i=1; i<num_timings; i++){
         if (time_total_global[i] != 0.0) {
            remainder = i%2;
            if (remainder == 1){
               bytes_updated = count_updated_global*4;
            }
            else {
               bytes_updated = count_updated_global*8;
            }
            
            printf("                    %lf       %8.0lf        %d\n",time_total_global[i],
                  (double)(bytes_updated)/time_total_global[i],(i+1)/num_timings_cycle);
         }
         else{
            printf("                    %lf       %8.0s        %d\n",0.0,"n/a",(i+1)/num_timings_cycle);
         }
      }
      
      printf("\n");
      printf("======================================\n");
      printf("    End L7_Update testing             \n");
      printf("======================================\n");
      printf("\n");
      
   }

   /*
    * Testing complete
    */
   free(time_total_global);
   
   return;

}

