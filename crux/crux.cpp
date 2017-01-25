/*
 *  Copyright (c) 2014, Los Alamos National Security, LLC.
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
 *  
 *  Authors: Brian Atkinson          bwa@g.clemson.edu
             Bob Robey        XCP-2  brobey@lanl.gov
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <assert.h>
#include "PowerParser/PowerParser.hh"

#include "crux.h"
#include "timer/timer.h"
#include "fmemopen.h"

#ifdef HAVE_HDF5
#include "hdf5.h"
#endif
#ifdef HAVE_MPI
#include "mpi.h"
#endif

const bool CRUX_TIMING = true;
bool do_crux_timing = false;

bool h5_spoutput;

#define RESTORE_NONE     0
#define RESTORE_RESTART  1
#define RESTORE_ROLLBACK 2

#ifndef DEBUG
#define DEBUG 0
#endif
#undef DEBUG_RESTORE_VALS

using namespace std;
using PP::PowerParser;
// Pointers to the various objects.
PowerParser *parse;

char checkpoint_directory[] = "checkpoint_output";
int cp_num, rs_num;
int *backup;
void **crux_data;
size_t *crux_data_size;
#ifdef HAVE_HDF5
bool USE_HDF5 = true; //MSB
hid_t h5_fid;
hid_t h5_gid_c, h5_gid_m, h5_gid_s;
herr_t h5err;
#endif

FILE *crux_time_fp;
struct timeval tcheckpoint_time;
struct timeval trestore_time;
int checkpoint_timing_count = 0;
float checkpoint_timing_sum = 0.0f;
float checkpoint_timing_size = 0.0f;
int rollback_attempt = 0;
FILE *store_fp, *restore_fp;
#ifdef HAVE_MPI
MPI_Offset store_offset, restore_offset;
static MPI_File mpi_store_fp, mpi_restore_fp;
#endif
static int mype = 0;

Crux::Crux(int crux_type_in, int num_of_rollback_states_in, bool restart)
{
#ifdef HAVE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&mype);

   for (int i=0; i<15; i++){
     crux_datatype[i] = 0;
   }
#endif

   num_of_rollback_states = num_of_rollback_states_in;
   crux_type = crux_type_in;
   checkpoint_counter = 0;

   if (crux_type != CRUX_NONE || restart){
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

void Crux::store_MallocPlus(MallocPlus memory){
   malloc_plus_memory_entry *memory_item;

   for (memory_item = memory.memory_entry_by_name_begin(); 
      memory_item != memory.memory_entry_by_name_end();
      memory_item = memory.memory_entry_by_name_next() ){

      void *mem_ptr = memory_item->mem_ptr;
      if ((memory_item->mem_flags & RESTART_DATA) == 0) continue;

      int num_elements = 1;
      for (uint i = 0; i < memory_item->mem_ndims; i++){
	 num_elements *= memory_item->mem_nelem[i];
      }

      //if (DEBUG) {
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
        printf("MallocPlus ptr  %p: name %10s ptr %p dims %lu nelem (",
           mem_ptr,memory_item->mem_name,memory_item->mem_ptr,memory_item->mem_ndims);

        char nelemstring[80];
        char *str_ptr = nelemstring;
        str_ptr += sprintf(str_ptr,"%lu", memory_item->mem_nelem[0]);
        for (uint i = 1; i < memory_item->mem_ndims; i++){
           str_ptr += sprintf(str_ptr,", %lu", memory_item->mem_nelem[i]);
        }
        printf("%12s",nelemstring);

        printf(") elsize %lu flags %d capacity %lu\n",
           memory_item->mem_elsize,memory_item->mem_flags,memory_item->mem_capacity);
#endif
      //}

#ifdef HAVE_HDF5
      if(USE_HDF5) {
        //
        // Create dataspace.  Setting maximum size to NULL sets the maximum
        // size to be the current size.
        //
        
        hid_t sid;
        sid = H5Screate_simple (memory_item->mem_ndims, (hsize_t*)memory_item->mem_nelem, NULL);
        
        hid_t memtype;
        hid_t filetype;
        if (memory_item->mem_elsize == 4){
          memtype = H5T_NATIVE_INT;
          filetype = H5T_NATIVE_INT;
        } else {
          memtype = H5T_NATIVE_DOUBLE;
          if(h5_spoutput) {
            filetype = H5T_NATIVE_FLOAT;
          } else {
            filetype = H5T_NATIVE_DOUBLE;
          }
        }
        
        hid_t gid;
        if( (strstr(memory_item->mem_name,"mesh") !=NULL) ||
            (strncmp(memory_item->mem_name,"i",1) == 0) || 
            (strncmp(memory_item->mem_name,"j",1) == 0) || 
            (strncmp(memory_item->mem_name,"level",5) == 0) ) {
          gid = h5_gid_m;
        }
        else if( (strstr(memory_item->mem_name,"state") != NULL) || 
                 (strncmp(memory_item->mem_name,"H",1) == 0) ||
                 (strncmp(memory_item->mem_name,"U",1) == 0) ||
                 (strncmp(memory_item->mem_name,"V",1) == 0) ) {
          gid = h5_gid_s; 
        }
        else {
          gid = h5_gid_c;
        }

        hid_t plist_id;
        plist_id = H5P_DEFAULT; 
#ifdef HAVE_MPI
        //
        // Create property list for collective dataset write.
        //
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

#endif
        if( (strncmp(memory_item->mem_name,"state_long_vals",15) == 0) ||
            (strncmp(memory_item->mem_name,"mesh_double_vals",16) == 0) ) {
          hid_t aid;
          aid = H5Acreate2(gid, memory_item->mem_name, filetype, sid, H5P_DEFAULT, H5P_DEFAULT);
          // Write the attribute data.
          h5err = H5Awrite(aid, filetype, mem_ptr);
          h5err = H5Aclose(aid);

        } else if( (strstr(memory_item->mem_name,"_timer") !=NULL) ||
                   (strstr(memory_item->mem_name,"_counters") !=NULL) ||
                   (strstr(memory_item->mem_name,"int_dist_vals") !=NULL)   ) {
          
          hid_t did;
          hsize_t dims[2], start[2], count[2];
          hid_t sid1;

          dims[0] = npes;
          dims[1] =(hsize_t)memory_item->mem_nelem[0];

          count[0] = 1;
          count[1] = 1;

          start[0] = mype;
          start[1] = 0;
          sid1 = H5Screate_simple (2, dims, NULL);

          did = H5Dcreate2 (gid, memory_item->mem_name, filetype, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          H5Sselect_hyperslab(sid1, H5S_SELECT_SET, start, NULL, count, NULL ); 
          h5err = H5Dwrite (did, memtype, H5S_ALL, sid1, plist_id, mem_ptr);
          
          if(H5Sclose(sid1) < 0)
            printf("HDF5: Could not close dataspace \n");
          if(H5Dclose(did) < 0)
            printf("HDF5: Could not close dataset %s \n",memory_item->mem_name);

        } else {

          hid_t did;
          did = H5Dcreate2 (gid, memory_item->mem_name, filetype, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          h5err = H5Dwrite (did, memtype, H5S_ALL, H5S_ALL, plist_id, mem_ptr);
          if(H5Dclose(did) < 0)
            printf("HDF5: Could not close dataset %s \n",memory_item->mem_name);
        }
        if(H5Sclose(sid) < 0)
          printf("HDF5: Could not close dataspace \n");
#ifdef HAVE_MPI
        if(H5Pclose(plist_id) < 0)
          printf("HDF5: Could not close property list \n");
#endif
      }
#endif
#ifdef HAVE_MPI
      MPI_Offset save_offset;
            if (num_elements > 9000){
               save_offset = store_offset;
            }
#endif
      store_field_header(memory_item->mem_name,30);
      if (memory_item->mem_flags & REPLICATED_DATA) { 
         if (memory_item->mem_elsize == 4){
            store_replicated_int_array((int *)mem_ptr, num_elements);
         } else {
            store_replicated_double_array((double *)mem_ptr, num_elements);
         }
      } else {
#ifdef HAVE_MPI
         int num_elements_global;
         MPI_Allreduce(&num_elements, &num_elements_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if (memory_item->mem_elsize == 4){
#ifdef DEBUG_RESTORE_VALS
            if (num_elements > 9000){
               printf("GOT HERE\n");
               char test_name[30];
               MPI_Status status;
               save_offset = 208;
               MPI_File_read_at(mpi_store_fp, save_offset, (void *)test_name, 30, MPI_CHAR, &status);
               printf("DEBUG -- name %s\n",test_name);
            }
            store_distributed_int_array((int *)mem_ptr, num_elements, num_elements_global,
               DISTRIBUTED_INT_LOCAL_DATA, DISTRIBUTED_INT_GLOBAL_DATA);
            if (num_elements > 9000){
               printf("GOT HERE\n");
               char test_name[30];
               MPI_Status status;
               save_offset = 208;
               MPI_File_read_at(mpi_store_fp, save_offset, (void *)test_name, 30, MPI_CHAR, &status);
               printf("DEBUG -- name %s\n",test_name);
            }
#endif
         } else {
            //store_distributed_double_array((double *)mem_ptr, num_elements, num_elements_global,
            //   DISTRIBUTED_DOUBLE_LOCAL_DATA, DISTRIBUTED_DOUBLE_GLOBAL_DATA);
         }
#else
         if (memory_item->mem_elsize == 4){
            store_int_array((int *)mem_ptr, num_elements);
         } else {
            store_double_array((double *)mem_ptr, num_elements);
         }
#endif
      }
   }

}

void Crux::store_begin(size_t nsize, int ncycle)
{
   cp_num = checkpoint_counter % num_of_rollback_states;

   cpu_timer_start(&tcheckpoint_time);

   if(crux_type == CRUX_IN_MEMORY){
      if (crux_data[cp_num] != NULL) free(crux_data[cp_num]);
      crux_data[cp_num] = (int *)malloc(nsize);
      crux_data_size[cp_num] = nsize;
      store_fp = fmemopen(crux_data[cp_num], nsize, "w");
   }
   if(crux_type == CRUX_DISK){
      char backup_file[60];

#ifdef HAVE_HDF5
      hid_t plist_id;

      if(USE_HDF5) {
        sprintf(backup_file,"%s/backup%05d.h5",checkpoint_directory,ncycle);

        
        plist_id = H5P_DEFAULT; 
#  ifdef HAVE_MPI
        int mpiInitialized = 0;
        bool phdf5 = false;
        if (MPI_SUCCESS == MPI_Initialized(&mpiInitialized)) {
          phdf5 = true;
        }

        // 
        // Set up file access property list with parallel I/O access
        //
        if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0)
          printf("HDF5: Could not create property list \n");

	if( H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) < 0)
          printf("HDF5: Could set libver bounds \n");

        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#  endif
        h5_fid = H5Fcreate(backup_file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        if(!h5_fid){
          printf("HDF5: Could not write HDF5 %s at iteration %d\n",backup_file,ncycle);
        }
        if( (h5_gid_c = H5Gcreate(h5_fid, "clamr", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ) < 0) 
          printf("HDF5: Could not create \"clamr\" group \n");
        if( (h5_gid_m = H5Gcreate(h5_fid, "mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ) < 0)
          printf("HDF5: Could not create \"mesh\" group \n");
        if( (h5_gid_s = H5Gcreate(h5_fid, "state", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ) < 0)
          printf("HDF5: Could not create \"state\" group \n");
      }

#  ifdef HAVE_MPI
      if(H5Pclose(plist_id) < 0)
        printf("HDF5: Could not close property list \n");
#  endif

#endif
      sprintf(backup_file,"%s/backup%05d.crx",checkpoint_directory,ncycle);
#ifdef HAVE_MPI
      int iret = MPI_File_open(MPI_COMM_WORLD, backup_file, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &mpi_store_fp);
      if(iret != MPI_SUCCESS) {
         printf("Could not write %s at iteration %d\n",backup_file,ncycle);
      }
      store_offset = 0;
#else
      store_fp = fopen(backup_file,"w");
      if(!store_fp){
         printf("Could not write %s at iteration %d\n",backup_file,ncycle);
      }
#endif

      if (mype == 0) {
        char symlink_file[60];
        sprintf(symlink_file,"%s/backup%1d.crx",checkpoint_directory,cp_num);
        symlink(backup_file, symlink_file);
//      int ireturn = symlink(backup_file, symlink_file);
//      if (ireturn == -1) {
//         printf("Warning: error returned with symlink call for file %s and symlink %s\n",
//                backup_file,symlink_file);
//      }
      }
   }

   if (do_crux_timing){
      checkpoint_timing_size += nsize;
   }
}

void Crux::store_field_header(const char *name, int name_size){
   assert(name != NULL);
   char out_name[name_size];
   for (int i = 0; i< name_size; i++){
     out_name[i] = ' ';
   }
   snprintf(out_name,name_size,"%-s",name);

#ifdef HAVE_MPI
   MPI_Status status;
   MPI_Offset etype_offset;
   MPI_Offset test_offset;
   if (mype == 0) {
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
      printf("\n%d:STORING HEADER %s etype_offset %ld TEST_offset %ld store_offset %ld\n",mype,out_name,etype_offset,test_offset,store_offset);
#endif
      int iret;
      iret = MPI_File_write_at(mpi_store_fp, store_offset, (void *)out_name, name_size, MPI_CHAR, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_write_at -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_position(mpi_store_fp, &etype_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_position -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_byte_offset(mpi_store_fp, etype_offset, &test_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_byte_offset -- line %d %s\n",__LINE__,__FILE__);
      }

#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
      char test_name[30];
      iret = MPI_File_read_at(mpi_store_fp, store_offset, (void *)test_name, name_size, MPI_CHAR, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_read_at -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("DEBUG -- name %s\n",test_name);
#endif
   }
   MPI_Barrier(MPI_COMM_WORLD);
   store_offset += name_size*sizeof(char);
#ifdef DEBUG_RESTORE_VALS
   int iret;
   if (mype == 0){
      int count;
      iret = MPI_Get_count(&status, MPI_CHAR, &count);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("%d:Wrote %d characters at line %d in file %s\n",mype,count,__LINE__,__FILE__);
      printf("%d:etype_offset %ld TEST_offset %ld store_offset %ld\n",mype,etype_offset,test_offset,store_offset);
   }
#endif

#else
   assert(store_fp != NULL);
   fwrite(out_name,sizeof(char),name_size,store_fp);
#endif
}

void Crux::store_bools(bool *bool_vals, size_t nelem)
{
   assert(bool_vals != NULL && store_fp != NULL);
   fwrite(bool_vals,sizeof(bool),nelem,store_fp);
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

void Crux::store_sizets(size_t *size_t_vals, size_t nelem)
{
   assert(size_t_vals != NULL && store_fp != NULL);
   fwrite(size_t_vals,sizeof(size_t),nelem,store_fp);
}

void Crux::store_doubles(double *double_vals, size_t nelem)
{
   assert(double_vals != NULL && store_fp != NULL);
   fwrite(double_vals,sizeof(double),nelem,store_fp);
}

void Crux::store_int_array(int *int_array, size_t nelem)
{
#ifdef HAVE_MPI
   assert(int_array != NULL);
   MPI_Status status;
   int iret;
   iret = MPI_File_write_shared(mpi_store_fp, int_array, (int)nelem, MPI_INT, &status);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_write_shared -- line %d %s\n",__LINE__,__FILE__);
   }
   store_offset += nelem*sizeof(int);
   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_RESTORE_VALS
   int count;
   iret = MPI_Get_count(&status, MPI_INT, &count);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
   }
   printf("%d:Wrote %d integers at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   assert(int_array != NULL && store_fp != NULL);
   fwrite(int_array,sizeof(int),nelem,store_fp);
#endif
}

void Crux::store_long_array(long long *long_array, size_t nelem)
{
   assert(long_array != NULL && store_fp != NULL);
   fwrite(long_array,sizeof(long long),nelem,store_fp);
}

void Crux::store_float_array(float *float_array, size_t nelem)
{
   assert(float_array != NULL && store_fp != NULL);
   fwrite(float_array,sizeof(float),nelem,store_fp);
}

void Crux::store_double_array(double *double_array, size_t nelem)
{
#ifdef HAVE_MPI
   assert(double_array != NULL);
   MPI_Status status;
   int iret;
   iret = MPI_File_write_shared(mpi_store_fp, double_array, (int)nelem, MPI_DOUBLE, &status);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_write_shared -- line %d %s\n",__LINE__,__FILE__);
   }
   store_offset += nelem*sizeof(double);
   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_RESTORE_VALS
   int count;
   MPI_Get_count(&status, MPI_DOUBLE, &count);
   printf("%d:Wrote %d doubles at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   assert(double_array != NULL && store_fp != NULL);
   fwrite(double_array,sizeof(double),nelem,store_fp);
#endif
}

void Crux::store_replicated_int_array(int *int_array, size_t nelem)
{
#ifdef HAVE_MPI
   assert(int_array != NULL);
   MPI_Status status;
   MPI_Offset etype_offset;
   MPI_Offset test_offset;
   int iret;
   if (mype == 0) {
      iret = MPI_File_write_at(mpi_store_fp, store_offset, int_array, (int)nelem, MPI_INT, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_write_at -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_position(mpi_store_fp, &etype_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_position -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_byte_offset(mpi_store_fp, etype_offset, &test_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_byte_offset -- line %d %s\n",__LINE__,__FILE__);
      }
   }
   store_offset += nelem*sizeof(int);
   iret = MPI_Barrier(MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_Barrier -- line %d %s\n",__LINE__,__FILE__);
   }
#ifdef DEBUG_RESTORE_VALS
   if (mype == 0) {
      int count;
      iret = MPI_Get_count(&status, MPI_INT, &count);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_write_at -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("%d:Wrote %d integers at line %d in file %s\n",mype,count,__LINE__,__FILE__);
      printf("%d:etype_offset %ld TEST_offset %ld store_offset %ld\n",mype,etype_offset,test_offset,store_offset);
   }
#endif

#else
   assert(int_array != NULL && store_fp != NULL);
   fwrite(int_array,sizeof(int),nelem,store_fp);
#endif
}

void Crux::store_replicated_double_array(double *double_array, size_t nelem)
{
#ifdef HAVE_MPI
   assert(double_array != NULL);
   MPI_Status status;
   MPI_Offset etype_offset;
   MPI_Offset test_offset;
   int iret;
   if (mype == 0) {
      iret = MPI_File_write_at(mpi_store_fp, store_offset, double_array, (int)nelem, MPI_DOUBLE, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_write_at -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_position(mpi_store_fp, &etype_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_position -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_byte_offset(mpi_store_fp, etype_offset, &test_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_byte_offset -- line %d %s\n",__LINE__,__FILE__);
      }
   }
   store_offset += nelem*sizeof(double);
   iret = MPI_Barrier(MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Barrier -- line %d %s\n",__LINE__,__FILE__);
   }
#ifdef DEBUG_RESTORE_VALS
   if (mype == 0) {
      int count;
      iret = MPI_Get_count(&status, MPI_DOUBLE, &count);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("%d:Wrote %d doubles at line %d in file %s\n",mype,count,__LINE__,__FILE__);
      printf("%d:etype_offset %ld TEST_offset %ld store_offset %ld\n",mype,etype_offset,test_offset,store_offset);
   }
#endif

#else
   assert(double_array != NULL && store_fp != NULL);
   fwrite(double_array,sizeof(double),nelem,store_fp);
#endif
}

void Crux::store_distributed_int_array(int *int_array, size_t nelem, size_t nelem_global,
      int local_flags, int global_flags)
{
#ifdef HAVE_MPI
   assert(int_array != NULL);
   MPI_Datatype *local_datatype;
   if (local_flags != 0){
     local_datatype = get_crux_datatype(local_flags);
   }
   MPI_Datatype *global_datatype = get_crux_datatype(global_flags);
   MPI_Status status;
   //printf("writing crux data ncells %d at %ld datatype size %ld\n",
   //  int_array[0],store_offset,sizeof(global_datatype));
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
   printf("\n%d:DISTRIBUTED INT ARRAY  store_offset %ld\n\n",mype,store_offset);
#endif
   int iret;
   iret = MPI_File_set_view(mpi_store_fp, store_offset, MPI_INT, *global_datatype, "native", MPI_INFO_NULL);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_set_view -- line %d %s\n",__LINE__,__FILE__);
   }
   if (local_flags != 0){
     iret = MPI_File_write_all(mpi_store_fp, int_array, 1, *local_datatype, &status); 
   } else {
     iret = MPI_File_write_all(mpi_store_fp, int_array, nelem, MPI_INT, &status); 
   }
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_write_all -- line %d %s\n",__LINE__,__FILE__);
   }
   iret = MPI_File_set_view(mpi_store_fp, store_offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_set_view -- line %d %s\n",__LINE__,__FILE__);
   }
   store_offset += nelem_global*sizeof(int);
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
   printf("\n%d:After store_offset -- DISTRIBUTED INT ARRAY  store_offset %ld sizeof %ld\n\n",
       mype,store_offset,nelem_global*sizeof(int));
#endif
   iret = MPI_Barrier(MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Barrier -- line %d %s\n",__LINE__,__FILE__);
   }
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
   printf("%d:store_offset %ld\n",mype,store_offset);
#endif
#ifdef DEBUG_RESTORE_VALS
   int count;
   if (local_flags != 0){
      iret = MPI_Get_count(&status, *local_datatype, &count);
   } else {
      iret = MPI_Get_count(&status, MPI_INT, &count);
   }
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
   }
   printf("%d:Wrote %d integers at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   assert(int_array != NULL && store_fp != NULL);
   fwrite(int_array,sizeof(int),nelem,store_fp);
#endif
}
void Crux::store_distributed_double_array(double *double_array, size_t nelem, size_t nelem_global,
       int local_flags, int global_flags)
{
#ifdef HAVE_MPI
   assert(double_array != NULL);
   MPI_Datatype *local_datatype = get_crux_datatype(DISTRIBUTED_DOUBLE_LOCAL_DATA);
   MPI_Datatype *global_datatype = get_crux_datatype(DISTRIBUTED_DOUBLE_GLOBAL_DATA);
   MPI_Status status;
   //MPI_File_write_shared(mpi_store_fp, double_array, nelem, local_datatype, &status);
#ifdef DEBUG_RESTORE_VALS
   int count;
   MPI_Get_count(&status, MPI_DOUBLE_PRECISION, &count);
   printf("%d:Wrote %d doubles at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   assert(double_array != NULL && store_fp != NULL);
   fwrite(double_array,sizeof(double),nelem,store_fp);
#endif
}

void Crux::store_end(void)
{
#ifdef HAVE_MPI
   int iret;
   iret = MPI_File_close(&mpi_store_fp);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_close -- line %d %s\n",__LINE__,__FILE__);
   }
#else
   assert(store_fp != NULL);
   fclose(store_fp);
#endif

#ifdef HAVE_HDF5
   if(USE_HDF5) {
     if(H5Gclose(h5_gid_c) < 0)
       printf("HDF5: Could not close clamr group \n");
     if(H5Gclose(h5_gid_m) < 0)
       printf("HDF5: Could not close mesh group \n");
     if(H5Gclose(h5_gid_s) < 0)
       printf("HDF5: Could not close state group \n");
     if(H5Fclose(h5_fid) != 0) {
       printf("HDF5: Could not close HDF5 file \n");
     }
   }
#endif

   double checkpoint_total_time = cpu_timer_stop(tcheckpoint_time);

   if (do_crux_timing){
      fprintf(crux_time_fp, "Total time for checkpointing was %g seconds\n", checkpoint_total_time);
      checkpoint_timing_count++;
      checkpoint_timing_sum += checkpoint_total_time;
   }

   checkpoint_counter++;
}

int restore_type = RESTORE_NONE;

void Crux::restore_MallocPlus(MallocPlus memory){
   char test_name[34];
   malloc_plus_memory_entry *memory_item;
   for (memory_item = memory.memory_entry_by_name_begin(); 
      memory_item != memory.memory_entry_by_name_end();
      memory_item = memory.memory_entry_by_name_next() ){

      void *mem_ptr = memory_item->mem_ptr;
      if ((memory_item->mem_flags & RESTART_DATA) == 0) continue;

      int num_elements = 1;
      for (uint i = 0; i < memory_item->mem_ndims; i++){
	 num_elements *= memory_item->mem_nelem[i];
      }

      //if (DEBUG) {
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
        printf("MallocPlus ptr  %p: name %10s ptr %p dims %lu nelem (",
           mem_ptr,memory_item->mem_name,memory_item->mem_ptr,memory_item->mem_ndims);

        char nelemstring[80];
        char *str_ptr = nelemstring;
        str_ptr += sprintf(str_ptr,"%lu", memory_item->mem_nelem[0]);
        for (uint i = 1; i < memory_item->mem_ndims; i++){
           str_ptr += sprintf(str_ptr,", %lu", memory_item->mem_nelem[i]);
        }
        printf("%12s",nelemstring);

        printf(") elsize %lu flags %d capacity %lu\n",
           memory_item->mem_elsize,memory_item->mem_flags,memory_item->mem_capacity);
#endif
      //}

      restore_field_header(test_name,30);
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
      printf("%d: field name %s %s\n",mype,test_name,memory_item->mem_name);
#endif
      if (strcmp(test_name,memory_item->mem_name) != 0) {
         printf("ERROR in restore checkpoint for %s %s\n",test_name,memory_item->mem_name);
#ifdef HAVE_MPI
         MPI_Finalize();
#endif
         exit(-1);
      }

      if (memory_item->mem_flags & REPLICATED_DATA) { 
         if (memory_item->mem_elsize == 4){
            restore_replicated_int_array((int *)mem_ptr, num_elements);
         } else {
            restore_replicated_double_array((double *)mem_ptr, num_elements);
         }
      } else {
#ifdef HAVE_MPI
         int num_elements_global;
         MPI_Allreduce(&num_elements, &num_elements_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if (memory_item->mem_elsize == 4){
#ifdef DEBUG_RESTORE_VALS
            restore_distributed_int_array((int *)mem_ptr, num_elements, num_elements_global,
               DISTRIBUTED_INT_LOCAL_DATA, DISTRIBUTED_INT_GLOBAL_DATA);
#endif
         } else {
            //restore_distributed_double_array((double *)mem_ptr, num_elements, num_elements_global,
            //   DISTRIBUTED_DOUBLE_LOCAL_DATA, DISTRIBUTED_DOUBLE_GLOBAL_DATA);
         }
#else
         if (memory_item->mem_elsize == 4){
            restore_int_array((int *)mem_ptr, num_elements);
         } else {
            restore_double_array((double *)mem_ptr, num_elements);
         }
#endif
      }
   }

}

void Crux::restore_begin(char *restart_file, int rollback_counter)
{
   rs_num = rollback_counter % num_of_rollback_states;

   cpu_timer_start(&trestore_time);

   if (restart_file != NULL){
      if (mype == 0) {
         printf("\n  ================================================================\n");
         printf(  "  Restoring state from disk file %s\n",restart_file);
         printf(  "  ================================================================\n\n");
      }
#ifdef HAVE_MPI
      int iret = MPI_File_open(MPI_COMM_WORLD, restart_file, MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &mpi_restore_fp);
      if(iret != MPI_SUCCESS){
         //printf("Could not write %s at iteration %d\n",restart_file,crux_int_vals[8]);
         printf("Could not open restart file %s\n",restart_file);
      }
      restore_offset = 0;
#else
      restore_fp = fopen(restart_file,"r");
      if(!restore_fp){
         //printf("Could not write %s at iteration %d\n",restart_file,crux_int_vals[8]);
         printf("Could not open restart file %s\n",restart_file);
      }
#endif
      restore_type = RESTORE_RESTART;
   } else if(crux_type == CRUX_IN_MEMORY){
      printf("Restoring state from memory rollback number %d rollback_counter %d\n",rs_num,rollback_counter);
      restore_fp = fmemopen(crux_data[rs_num], crux_data_size[rs_num], "r");
      restore_type = RESTORE_ROLLBACK;
   } else if(crux_type == CRUX_DISK){
      char backup_file[60];

      sprintf(backup_file,"%s/backup%d.crx",checkpoint_directory,rs_num);
      printf("Restoring state from disk file %s rollback_counter %d\n",backup_file,rollback_counter);
      restore_fp = fopen(backup_file,"r");
      if(!restore_fp){
         //printf("Could not write %s at iteration %d\n",backup_file,crux_int_vals[8]);
         printf("Could not open restore file %s\n",backup_file);
      }
      restore_type = RESTORE_ROLLBACK;
   }
}

void Crux::restore_field_header(char *name, int name_size)
{
#ifdef HAVE_MPI
   assert(name != NULL);
   for (int i = 0; i < name_size; i++){
     name[i] = ' ';
   }
   MPI_Status status;
   MPI_Offset etype_offset;
   MPI_Offset test_offset;
   int iret;
   if (mype == 0) {
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
      printf("\n%d:RESTORING HEADER etype_offset %ld TEST_offset %ld restore_offset %ld\n",mype,etype_offset,test_offset,restore_offset);
#endif
      iret = MPI_File_read_at(mpi_restore_fp, restore_offset, name, name_size, MPI_CHAR, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_read_at -- line %d %s\n",__LINE__,__FILE__);
      }
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
      printf("%d: DEBUG -- name %s\n",mype,name);
#endif
      iret = MPI_File_get_position(mpi_restore_fp, &etype_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_position -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_byte_offset(mpi_restore_fp, etype_offset, &test_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_byte_offset -- line %d %s\n",__LINE__,__FILE__);
      }
   }
   iret = MPI_Bcast(name, 30, MPI_CHAR, 0, MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Bcast -- line %d %s\n",__LINE__,__FILE__);
   }
   restore_offset += name_size*sizeof(char);
#ifdef DEBUG_RESTORE_VALS
   if (mype == 0) {
      int count;
      iret = MPI_Get_count(&status, MPI_CHAR, &count);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("%d:Read %d characters at line %d in file %s\n",mype,count,__LINE__,__FILE__);
      printf("%d:etype_offset %ld TEST_offset %ld restore_offset %ld\n",mype,etype_offset,test_offset,restore_offset);
   }
#endif

#else
   int name_read = fread(name,sizeof(char),name_size,restore_fp);
   if (name_read != name_size){
      printf("Warning: number of elements read %d is not equal to request %d\n",name_read,name_size);
   }
#endif
}

void Crux::restore_bools(bool *bool_vals, size_t nelem)
{
   size_t nelem_read = fread(bool_vals,sizeof(bool),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
}

void Crux::restore_ints(int *int_vals, size_t nelem)
{
   size_t nelem_read = fread(int_vals,sizeof(int),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
}

void Crux::restore_longs(long long *long_vals, size_t nelem)
{
   size_t nelem_read = fread(long_vals,sizeof(long),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
}

void Crux::restore_sizets(size_t *size_t_vals, size_t nelem)
{
   size_t nelem_read = fread(size_t_vals,sizeof(size_t),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
}

void Crux::restore_doubles(double *double_vals, size_t nelem)
{
   size_t nelem_read = fread(double_vals,sizeof(double),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
}

int *Crux::restore_int_array(int *int_array, size_t nelem)
{
#ifdef HAVE_MPI
   assert(int_array != NULL);
   MPI_Status status;
   int iret;
   iret = MPI_File_read_shared(mpi_restore_fp, int_array, (int)nelem, MPI_INT, &status);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_read_shared -- line %d %s\n",__LINE__,__FILE__);
   }
   restore_offset += nelem*sizeof(int);
   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_RESTORE_VALS
   int count;
   iret = MPI_Get_count(&status, MPI_INT, &count);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
   }
   printf("%d:Read %d integers at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   size_t nelem_read = fread(int_array,sizeof(int),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
#endif
   return(int_array);
}

long long *Crux::restore_long_array(long long *long_array, size_t nelem)
{
   size_t nelem_read = fread(long_array,sizeof(long long),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
   return(long_array);
}

float *Crux::restore_float_array(float *float_array, size_t nelem)
{
   size_t nelem_read = fread(float_array,sizeof(float),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
   return(float_array);
}

double *Crux::restore_double_array(double *double_array, size_t nelem)
{
#ifdef HAVE_MPI
   MPI_Status status;
   int iret;
   iret = MPI_File_read_shared(mpi_restore_fp, double_array, (int)nelem, MPI_DOUBLE, &status);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_read_shared -- line %d %s\n",__LINE__,__FILE__);
   }
   restore_offset += nelem*sizeof(double);
   iret = MPI_Barrier(MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Barrier -- line %d %s\n",__LINE__,__FILE__);
   }
#ifdef DEBUG_RESTORE_VALS
   int count;
   iret = MPI_Get_count(&status, MPI_DOUBLE, &count);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
   }
   printf("%d:Read %d doubles at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif
  
#else
   size_t nelem_read = fread(double_array,sizeof(double),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
#endif
   return(double_array);
}

int *Crux::restore_replicated_int_array(int *int_array, size_t nelem)
{
#ifdef HAVE_MPI
   assert(int_array != NULL);
   MPI_Status status;
   MPI_Offset etype_offset;
   MPI_Offset test_offset;
   int iret;
   if (mype == 0) {
      MPI_File_read_at(mpi_restore_fp, restore_offset, int_array, (int)nelem, MPI_INT, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_read_at -- line %d %s\n",__LINE__,__FILE__);
      }
      MPI_File_get_position(mpi_restore_fp, &etype_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_position -- line %d %s\n",__LINE__,__FILE__);
      }
      MPI_File_get_byte_offset(mpi_restore_fp, etype_offset, &test_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_byte_offset -- line %d %s\n",__LINE__,__FILE__);
      }
   }
   iret = MPI_Bcast(int_array, nelem, MPI_INT, 0, MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Bcast -- line %d %s\n",__LINE__,__FILE__);
   }
   restore_offset += nelem*sizeof(int);
   iret = MPI_Barrier(MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Barrier -- line %d %s\n",__LINE__,__FILE__);
   }
#ifdef DEBUG_RESTORE_VALS
   if (mype == 0) {
      int count;
      iret = MPI_Get_count(&status, MPI_INT, &count);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("%d:Read %d integers at line %d in file %s\n",mype,count,__LINE__,__FILE__);
      printf("%d:etype_offset %ld TEST_offset %ld restore_offset %ld\n",mype,etype_offset,test_offset,restore_offset);
   }
#endif

#else
   size_t nelem_read = fread(int_array,sizeof(int),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
#endif
   return(int_array);
}

double *Crux::restore_replicated_double_array(double *double_array, size_t nelem)
{
#ifdef HAVE_MPI
   MPI_Status status;
   MPI_Offset etype_offset;
   MPI_Offset test_offset;
   int iret;
   if (mype == 0) {
      iret = MPI_File_read_at(mpi_restore_fp, restore_offset, double_array, (int)nelem, MPI_DOUBLE, &status);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_read_at -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_position(mpi_restore_fp, &etype_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_position -- line %d %s\n",__LINE__,__FILE__);
      }
      iret = MPI_File_get_byte_offset(mpi_restore_fp, etype_offset, &test_offset);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_File_get_byte_offset -- line %d %s\n",__LINE__,__FILE__);
      }
   }
   MPI_Bcast(double_array, nelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   restore_offset += nelem*sizeof(double);
   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_RESTORE_VALS
   if (mype == 0) {
      int count;
      iret = MPI_Get_count(&status, MPI_DOUBLE, &count);
      if (iret != MPI_SUCCESS) {
         printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
      }
      printf("%d:Read %d doubles at line %d in file %s\n",mype,count,__LINE__,__FILE__);
      printf("%d:etype_offset %ld TEST_offset %ld restore_offset %ld\n",mype,etype_offset,test_offset,restore_offset);
   }
#endif
  
#else
   size_t nelem_read = fread(double_array,sizeof(double),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
#endif
   return(double_array);
}

int *Crux::restore_distributed_int_array(int *int_array, size_t nelem, size_t nelem_global,
    int local_flags, int global_flags)
{
#ifdef HAVE_MPI
   assert(int_array != NULL);
   MPI_Datatype *local_datatype;
   if (local_flags != 0) {
     local_datatype = get_crux_datatype(local_flags);
   }
   MPI_Datatype *global_datatype = get_crux_datatype(global_flags);
   MPI_Status status;
   int iret;
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
   //printf("reading crux data ncells %d at %ld datatype size %ld\n",int_array[0],restore_offset,sizeof(*local_datatype));
   printf("\n%d:DISTRIBUTED INT ARRAY  restore_offset %ld\n\n",mype,restore_offset);
#endif
   iret = MPI_File_set_view(mpi_restore_fp, restore_offset, MPI_INT, *global_datatype, "native", MPI_INFO_NULL);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_set_view -- line %d %s\n",__LINE__,__FILE__);
   }
   if (local_flags != 0) {
     iret = MPI_File_read_all(mpi_restore_fp, int_array, 1, *local_datatype, &status); 
   } else {
     iret = MPI_File_read_all(mpi_restore_fp, int_array, 1, *local_datatype, &status); 
   }
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_read_all -- line %d %s\n",__LINE__,__FILE__);
   }
   iret = MPI_File_set_view(mpi_restore_fp, restore_offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_set_view -- line %d %s\n",__LINE__,__FILE__);
   }
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
   printf("DEBUG -- ierr %d\n",ierr);
#endif
   restore_offset += nelem_global*sizeof(int);
#if defined(HAVE_MPI) && defined(DEBUG_RESTORE_VALS)
   printf("\n%d:after restore_offset DISTRIBUTED INT ARRAY  restore_offset %ld sizeof %ld\n\n",
         mype,restore_offset,nelem_global*sizeof(int));
#endif
   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_RESTORE_VALS
   int count;
   if (local_flags != 0) {
     iret = MPI_Get_count(&status, *local_datatype, &count);
   } else {
     iret = MPI_Get_count(&status, MPI_INT, &count);
   }
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Get_count -- line %d %s\n",__LINE__,__FILE__);
   }
   printf("%d:Read %d datatype at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   size_t nelem_read = fread(int_array,sizeof(int),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
#endif

   return(int_array);
}

double *Crux::restore_distributed_double_array(double *double_array, size_t nelem, size_t nelem_global,
        int local_flags, int global_flags)
{
#ifdef HAVE_MPI
   MPI_Datatype *local_datatype = get_crux_datatype(DISTRIBUTED_DOUBLE_LOCAL_DATA);
   MPI_Datatype *global_datatype = get_crux_datatype(DISTRIBUTED_DOUBLE_GLOBAL_DATA);
   MPI_Status status;
   int iret;
   iret = MPI_File_read_shared(mpi_restore_fp, double_array, (int)nelem, *local_datatype, &status);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_File_read_shared -- line %d %s\n",__LINE__,__FILE__);
   }
   restore_offset += sizeof(*global_datatype);
   iret = MPI_Barrier(MPI_COMM_WORLD);
   if (iret != MPI_SUCCESS) {
      printf("ERROR MPI_Barrier -- line %d %s\n",__LINE__,__FILE__);
   }
#ifdef DEBUG_RESTORE_VALS
   int count;
   MPI_Get_count(&status, MPI_DOUBLE, &count);
   printf("%d:Read %d doubles at line %d in file %s\n",mype,count,__LINE__,__FILE__);
#endif

#else
   size_t nelem_read = fread(double_array,sizeof(double),nelem,restore_fp);
   if (nelem_read != nelem){
      printf("Warning: number of elements read %lu is not equal to request %lu\n",nelem_read,nelem);
   }
#endif
  
   return(double_array);
}

void Crux::restore_end(void)
{
   double restore_total_time = cpu_timer_stop(trestore_time);

   if (do_crux_timing){
      if (restore_type == RESTORE_RESTART) {
         fprintf(crux_time_fp, "Total time for restore was %g seconds\n", restore_total_time);
      } else if (restore_type == RESTORE_ROLLBACK){
         fprintf(crux_time_fp, "Total time for rollback %d was %g seconds\n", rollback_attempt, restore_total_time);
      }
   }

   fclose(restore_fp);

}

int Crux::get_rollback_number()
{
  rollback_attempt++;
  return(checkpoint_counter % num_of_rollback_states);
}

void Crux::set_crux_type(int crux_type_in)
{
  crux_type = crux_type_in;
}
#ifdef HAVE_MPI
void Crux::set_crux_datatype(MPI_Datatype *datatype, int flags)
{
   int index = 0;

   if        (flags & 0x00001){
      index = 0;
   } else if (flags & 0x00002){
      index = 1;
   } else if (flags & 0x00004){
      index = 2;
   } else if (flags & 0x00008){
      index = 3;
   } else if (flags & 0x00010){
      index = 4;
   } else if (flags & 0x00020){
      index = 5;
   } else if (flags & 0x00040){
      index = 6;
   } else if (flags & 0x00080){
      index = 7;
   } else if (flags & 0x00100){
      index = 8;
   } else if (flags & 0x00200){
      index = 9;
   } else if (flags & 0x00400){
      index = 10;
   } else if (flags & 0x00800){
      index = 11;
   } else if (flags & 0x01000){
      index = 12;
   } else if (flags & 0x02000){
      index = 13;
   } else if (flags & 0x04000){
      index = 14;
   } else if (flags & 0x08000){
      index = 15;
   }

   //printf("%d: setting crux data type %d\n",__LINE__,index);
   if (index != 0) {
     //printf("%d: setting crux data type %d\n",__LINE__,index);
     MPI_Type_commit(datatype);
     crux_datatype[index] = datatype;
   }
}
MPI_Datatype *Crux::get_crux_datatype(int flags)
{
   int index = 0;

   if        (flags & 0x00001){
      index = 0;
   } else if (flags & 0x00002){
      index = 1;
   } else if (flags & 0x00004){
      index = 2;
   } else if (flags & 0x00008){
      index = 3;
   } else if (flags & 0x00010){
      index = 4;
   } else if (flags & 0x00020){
      index = 5;
   } else if (flags & 0x00040){
      index = 6;
   } else if (flags & 0x00080){
      index = 7;
   } else if (flags & 0x00100){
      index = 8;
   } else if (flags & 0x00200){
      index = 9;
   } else if (flags & 0x00400){
      index = 10;
   } else if (flags & 0x00800){
      index = 11;
   } else if (flags & 0x01000){
      index = 12;
   } else if (flags & 0x02000){
      index = 13;
   } else if (flags & 0x04000){
      index = 14;
   } else if (flags & 0x08000){
      index = 15;
   }

   if (index != 0) {
     return(crux_datatype[index]);
   }
   return(0);
}
void Crux::free_crux_datatype(int flags)
{
   int index = 0;

   if        (flags & 0x00001){
      index = 0;
   } else if (flags & 0x00002){
      index = 1;
   } else if (flags & 0x00004){
      index = 2;
   } else if (flags & 0x00008){
      index = 3;
   } else if (flags & 0x00010){
      index = 4;
   } else if (flags & 0x00020){
      index = 5;
   } else if (flags & 0x00040){
      index = 6;
   } else if (flags & 0x00080){
      index = 7;
   } else if (flags & 0x00100){
      index = 8;
   } else if (flags & 0x00200){
      index = 9;
   } else if (flags & 0x00400){
      index = 10;
   } else if (flags & 0x00800){
      index = 11;
   } else if (flags & 0x01000){
      index = 12;
   } else if (flags & 0x02000){
      index = 13;
   } else if (flags & 0x04000){
      index = 14;
   } else if (flags & 0x08000){
      index = 15;
   }

   if (index != 0) {
     MPI_Type_free(crux_datatype[index]);
   }
}
void Crux::free_all_crux_datatype()
{
   for (int index = 0; index < 16; index++){
     if (crux_datatype[index] != 0) {
       MPI_Type_free(crux_datatype[index]);
     }
   }
}
#endif
