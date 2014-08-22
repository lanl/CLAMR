#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <assert.h>

#include "crux/crux.h"
#include "timer/timer.h"
#include "fmemopen.h"

const bool CRUX_TIMING = true;
bool do_crux_timing = false;

using namespace std;

char checkpoint_directory[] = "checkpoint_output";
FILE *store_fp, *restore_fp;
int cp_num, rs_num;
int *backup;
void **crux_data;
size_t *crux_data_size;

FILE *crux_time_fp;
struct timeval tcheckpoint_time;
struct timeval trestore_time;
int checkpoint_timing_count = 0;
float checkpoint_timing_sum = 0.0f;
float checkpoint_timing_size = 0.0f;

Crux::Crux(int rollback_type_in, int num_of_rollback_states_in, bool restart)
{
   num_of_rollback_states = num_of_rollback_states_in;
   rollback_type = rollback_type_in;
   checkpoint_counter = 0;

   if (rollback_type != ROLLBACK_NONE || restart){
      do_crux_timing = CRUX_TIMING;
      struct stat stat_descriptor;
      if (stat(checkpoint_directory,&stat_descriptor) == -1){
        mkdir(checkpoint_directory,0777);
      }
   }

   crux_data = (void **)malloc(num_of_rollback_states*sizeof(void *));
   for (int i = 0; i < num_of_rollback_states; i++){
      crux_data[i] = NULL;
   }
   crux_data_size = (size_t *)malloc(num_of_rollback_states*sizeof(size_t));


   if (do_crux_timing){
      char checkpointtimelog[60];
      sprintf(checkpointtimelog,"%s/crux_timing.log",checkpoint_directory);
      crux_time_fp = fopen(checkpointtimelog,"w");
   }
}

Crux::~Crux()
{
   for (int i = 0; i < num_of_rollback_states; i++){
      free(crux_data[i]);
   }
   free(crux_data);
   free(crux_data_size);

   if (do_crux_timing){
      if (checkpoint_timing_count > 0) {
         printf("CRUX checkpointing time averaged %f msec, bandwidth %f Mbytes/sec\n",
                checkpoint_timing_sum/(float)checkpoint_timing_count*1.0e3,
                checkpoint_timing_size/checkpoint_timing_sum*1.0e-6);

         fprintf(crux_time_fp,"CRUX checkpointing time averaged %f msec, bandwidth %f Mbytes/sec\n",
                checkpoint_timing_sum/(float)checkpoint_timing_count*1.0e3,
                checkpoint_timing_size/checkpoint_timing_sum*1.0e-6);

      fclose(crux_time_fp);
      }
   }
}

void Crux::store_begin(size_t nsize, int ncycle)
{
   cp_num = checkpoint_counter % num_of_rollback_states;

   cpu_timer_start(&tcheckpoint_time);

   if(rollback_type == ROLLBACK_IN_MEMORY){
      if (crux_data[cp_num] != NULL) free(crux_data[cp_num]);
      crux_data[cp_num] = (int *)malloc(nsize);
      crux_data_size[cp_num] = nsize;
      store_fp = fmemopen(crux_data[cp_num], nsize, "w");
   }
   if(rollback_type == ROLLBACK_DISK){
      char backup_file[60];

      sprintf(backup_file,"%s/backup%05d.crx",checkpoint_directory,ncycle);
      store_fp = fopen(backup_file,"w");
      if(!store_fp){
         printf("Could not write %s at iteration %d\n",backup_file,ncycle);
      }

      char symlink_file[60];
      sprintf(symlink_file,"%s/backup%1d.crx",checkpoint_directory,cp_num);
      symlink(backup_file, symlink_file);
   }

   if (do_crux_timing){
      checkpoint_timing_size += nsize;
   }
}

void Crux::store_ints(int *int_vals, size_t nelem)
{
   assert(int_vals != NULL && store_fp != NULL);
   fwrite(int_vals,sizeof(int),nelem,store_fp);
}

void Crux::store_longs(long long *long_vals, size_t nelem)
{
   assert(long_vals != NULL && store_fp != NULL);
   fwrite(long_vals,sizeof(long long),nelem,store_fp);
}

void Crux::store_bools(bool *bool_vals, size_t nelem)
{
   assert(bool_vals != NULL && store_fp != NULL);
   fwrite(bool_vals,sizeof(bool),nelem,store_fp);
}

void Crux::store_doubles(double *double_vals, size_t nelem)
{
   assert(double_vals != NULL && store_fp != NULL);
   fwrite(double_vals,sizeof(double),nelem,store_fp);
}

void Crux::store_int_array(int *int_array, size_t nelem)
{
   assert(int_array != NULL && store_fp != NULL);
   fwrite(int_array,sizeof(int),nelem,store_fp);
}

void Crux::store_float_array(float *float_array, size_t nelem)
{
   assert(float_array != NULL && store_fp != NULL);
   fwrite(float_array,sizeof(float),nelem,store_fp);
}

void Crux::store_double_array(double *double_array, size_t nelem)
{
   assert(double_array != NULL && store_fp != NULL);
   fwrite(double_array,sizeof(double),nelem,store_fp);
}

void Crux::store_end(void)
{
   assert(store_fp != NULL);
   fclose(store_fp);

   double checkpoint_total_time = cpu_timer_stop(tcheckpoint_time);

   if (do_crux_timing){
      fprintf(crux_time_fp, "Total time for checkpointing was %g seconds\n", checkpoint_total_time);
      checkpoint_timing_count++;
      checkpoint_timing_sum += checkpoint_total_time;
   }

   checkpoint_counter++;
}

void Crux::restore_begin(char *restart_file, int rollback_counter)
{
   rs_num = rollback_counter % num_of_rollback_states;

   cpu_timer_start(&trestore_time);

   if (restart_file != NULL){
      printf("\n  ================================================================\n");
      printf(  "  Restoring state from disk file %s\n",restart_file);
      printf(  "  ================================================================\n\n");
      restore_fp = fopen(restart_file,"r");
      if(!restore_fp){
         //printf("Could not write %s at iteration %d\n",restart_file,crux_int_vals[8]);
         printf("Could not open restart file %s\n",restart_file);
      }
   } else if(rollback_type == ROLLBACK_IN_MEMORY){
      printf("Restoring state from memory rollback number %d rollback_counter %d\n",rs_num,rollback_counter);
      restore_fp = fmemopen(crux_data[rs_num], crux_data_size[rs_num], "r");
   } else if(rollback_type == ROLLBACK_DISK){
      char backup_file[60];

      sprintf(backup_file,"%s/backup%d.crx",checkpoint_directory,rs_num);
      printf("Restoring state from disk file %s rollback_counter %d\n",backup_file,rollback_counter);
      restore_fp = fopen(backup_file,"r");
      if(!restore_fp){
         //printf("Could not write %s at iteration %d\n",backup_file,crux_int_vals[8]);
         printf("Could not open restore file %s\n",backup_file);
      }
   }
}

void Crux::restore_ints(int *int_vals, size_t nelem)
{
   fread(int_vals,sizeof(int),nelem,restore_fp);
}

void Crux::restore_bools(bool *bool_vals, size_t nelem)
{
   fread(bool_vals,sizeof(bool),nelem,restore_fp);
}

void Crux::restore_longs(long long *long_vals, size_t nelem)
{
   fread(long_vals,sizeof(long),nelem,restore_fp);
}

void Crux::restore_doubles(double *double_vals, size_t nelem)
{
   fread(double_vals,sizeof(double),nelem,restore_fp);
}

int *Crux::restore_int_array(int *int_array, size_t nelem)
{
   fread(int_array,sizeof(int),nelem,restore_fp);
   return(int_array);
}

float *Crux::restore_float_array(float *float_array, size_t nelem)
{
   fread(float_array,sizeof(float),nelem,restore_fp);
   return(float_array);
}

double *Crux::restore_double_array(double *double_array, size_t nelem)
{
   fread(double_array,sizeof(double),nelem,restore_fp);
   return(double_array);
}

void Crux::restore_end(void)
{
   double restore_total_time = cpu_timer_stop(trestore_time);

   if (do_crux_timing){
      fprintf(crux_time_fp, "Total time for restore was %g seconds\n", restore_total_time);
   }

   fclose(restore_fp);

}
