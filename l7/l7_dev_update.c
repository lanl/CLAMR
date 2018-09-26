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
 */  
#include <stdlib.h>
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_DEV_UPDATE"

#ifdef HAVE_OPENCL
int L7_Dev_Update(
      cl_mem                  dev_data_buffer,
      const enum L7_Datatype  l7_datatype,
      const int               l7_id
      )
{
   /*
    * Purpose
    * =======
    * L7_Dev_Update collects into array data_buffer data located off-process,
    * appending it to owned (on-process) data data_buffer.
    * 
    * Arguments
    * =========
    * data_buffer        (input/output) void*
    *                    On input,
    *                    data_buffer[0:num_indices_owned-1] contains
    *                    data owned by process.
    *                    On output,
    *                    data_buffer[num_indices_owned, num_indices_needed-1]
    *                    contains the data collected from off-process.
    * 
    * l7_datatype        (input) const int*
    *                    The type of data contained in array data_buffer.
    * 
    * l7_id              (input) const int
    *                    Handle to database containing conmmunication
    *                    requirements.
    * 
    * Notes:
    * =====
    * 1) Serial compilation creates a no-op
    * 
    */
#if defined HAVE_MPI
   
   /*
    * Local variables
    */
   int
     ierr;                 /* Error code for return              */

   l7_id_database
     *l7_id_db;            /* database associated with l7_id.    */
   
   /*
    * Executable Statements
    */

   if (! l7.mpi_initialized){
      return(0);
   }
    
   if (l7.initialized !=1){
      ierr = 1;
      L7_ASSERT(l7.initialized == 1, "L7 not initialized", ierr);
   }
   
   /*
    * Check input.
    */
   
/*
   if (data_buffer == NULL){
      ierr = -1;
      L7_ASSERT( data_buffer != NULL, "data_buffer != NULL", ierr);
   }
*/
   
   if (l7_id <= 0){
      ierr = -1;
      L7_ASSERT( l7_id > 0, "l7_id <= 0", ierr);
   }
   
   if (l7.numpes == 1){
      ierr = L7_OK;
      return(ierr);
   }
   
   /*
    * Alias database associated with input l7_id
    */

   l7_id_db = l7p_set_database(l7_id);
   if (l7_id_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_id_db != NULL, "Failed to find database.", ierr);
   }
   
   l7.penum = l7_id_db->penum;
   
   if (l7_id_db->numpes == 1){ /* No-op */
      ierr = L7_OK;
      return(ierr);
   }
   
   /*
    * Setup parameters and data for compact device read
    */

   size_t num_indices_have   = l7_id_db->num_indices_have;
   size_t num_indices_owned  = l7_id_db->num_indices_owned;
   size_t num_indices_needed = l7_id_db->num_indices_needed;

   cl_command_queue command_queue = ezcl_get_command_queue();

   /*
    * Pull the data off the GPU and organize into regular array.
    */

   size_t pack_local_work_size = 128;
   size_t pack_global_work_size = ((num_indices_have + pack_local_work_size - 1) /pack_local_work_size) * pack_local_work_size;

   cl_mem dev_packed_data = NULL;
   cl_mem dev_data_buffer_add = NULL;
   int *packed_int_data = NULL;
   int *int_data_buffer = NULL;
   float *packed_float_data = NULL;
   float *float_data_buffer = NULL;
   double *packed_double_data = NULL;
   double *double_data_buffer = NULL;

   size_t ghost_local_work_size = 64;
   size_t ghost_global_work_size = ((num_indices_needed + ghost_local_work_size - 1) /ghost_local_work_size) * ghost_local_work_size;

   switch (l7_datatype) {
      case L7_INT:
         dev_packed_data = ezcl_malloc(NULL, "dev_packed_data", &num_indices_have, sizeof(cl_int), CL_MEM_READ_WRITE, 0);

         ezcl_set_kernel_arg(l7.kernel_pack_int_have_data, 0, sizeof(cl_int), (void *)&num_indices_have);
         ezcl_set_kernel_arg(l7.kernel_pack_int_have_data, 1, sizeof(cl_mem), (void *)&l7_id_db->dev_indices_have);
         ezcl_set_kernel_arg(l7.kernel_pack_int_have_data, 2, sizeof(cl_mem), (void *)&dev_data_buffer);
         ezcl_set_kernel_arg(l7.kernel_pack_int_have_data, 3, sizeof(cl_mem), (void *)&dev_packed_data);

         ezcl_enqueue_ndrange_kernel(command_queue, l7.kernel_pack_int_have_data,   1, NULL, &pack_global_work_size, &pack_local_work_size, NULL);

         packed_int_data = (int *)malloc(num_indices_have*sizeof(int));
         ezcl_enqueue_read_buffer(command_queue, dev_packed_data, CL_TRUE, 0, num_indices_have*sizeof(cl_int), &packed_int_data[0], NULL);
   
         int_data_buffer = (int *)malloc((num_indices_owned+num_indices_needed)*sizeof(int));
         for (unsigned int ii = 0; ii < num_indices_have; ii++){
            int_data_buffer[l7_id_db->indices_have[ii]] = packed_int_data[ii];
         }

         free(packed_int_data);

         ezcl_device_memory_delete(dev_packed_data);

         /*
          * Do the regular L7_Update across processor.
          */
   
         L7_Update (int_data_buffer, l7_datatype, l7_id);

         dev_data_buffer_add    = ezcl_malloc(NULL, "dev_data_buffer_add",    &num_indices_needed,     sizeof(cl_int), CL_MEM_READ_WRITE, 0);
         ezcl_enqueue_write_buffer(command_queue, dev_data_buffer_add, CL_TRUE,  0, num_indices_needed*sizeof(cl_int), &int_data_buffer[num_indices_owned],     NULL);

         free(int_data_buffer);

         // Fill in ghost
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_int_data, 0, sizeof(cl_int), &num_indices_owned);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_int_data, 1, sizeof(cl_int), (void *)&num_indices_needed);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_int_data, 2, sizeof(cl_mem), (void *)&dev_data_buffer);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_int_data, 3, sizeof(cl_mem), (void *)&dev_data_buffer_add);

         ezcl_enqueue_ndrange_kernel(command_queue, l7.kernel_copy_ghost_int_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, NULL);

         ezcl_device_memory_delete(dev_data_buffer_add);

         break;
      case L7_FLOAT:
         dev_packed_data = ezcl_malloc(NULL, "dev_packed_data", &num_indices_have, sizeof(cl_float), CL_MEM_READ_WRITE, 0);

         ezcl_set_kernel_arg(l7.kernel_pack_float_have_data, 0, sizeof(cl_int), (void *)&num_indices_have);
         ezcl_set_kernel_arg(l7.kernel_pack_float_have_data, 1, sizeof(cl_mem), (void *)&l7_id_db->dev_indices_have);
         ezcl_set_kernel_arg(l7.kernel_pack_float_have_data, 2, sizeof(cl_mem), (void *)&dev_data_buffer);
         ezcl_set_kernel_arg(l7.kernel_pack_float_have_data, 3, sizeof(cl_mem), (void *)&dev_packed_data);

         ezcl_enqueue_ndrange_kernel(command_queue, l7.kernel_pack_float_have_data,   1, NULL, &pack_global_work_size, &pack_local_work_size, NULL);

         packed_float_data = (float *)malloc(num_indices_have*sizeof(float));
         ezcl_enqueue_read_buffer(command_queue, dev_packed_data, CL_TRUE, 0, num_indices_have*sizeof(cl_float), &packed_float_data[0], NULL);
   
         float_data_buffer = (float *)malloc((num_indices_owned+num_indices_needed)*sizeof(float));
         for (unsigned int ii = 0; ii < num_indices_have; ii++){
            float_data_buffer[l7_id_db->indices_have[ii]] = packed_float_data[ii];
         }

         free(packed_float_data);

         ezcl_device_memory_delete(dev_packed_data);

         /*
          * Do the regular L7_Update across processor.
          */
   
         L7_Update (float_data_buffer, l7_datatype, l7_id);

         dev_data_buffer_add    = ezcl_malloc(NULL, "dev_data_buffer_add",    &num_indices_needed,     sizeof(cl_float), CL_MEM_READ_WRITE, 0);
         ezcl_enqueue_write_buffer(command_queue, dev_data_buffer_add, CL_TRUE,  0, num_indices_needed*sizeof(cl_float), &float_data_buffer[num_indices_owned],     NULL);

         free(float_data_buffer);

         // Fill in ghost
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_float_data, 0, sizeof(cl_int), &num_indices_owned);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_float_data, 1, sizeof(cl_int), (void *)&num_indices_needed);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_float_data, 2, sizeof(cl_mem), (void *)&dev_data_buffer);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_float_data, 3, sizeof(cl_mem), (void *)&dev_data_buffer_add);

         ezcl_enqueue_ndrange_kernel(command_queue, l7.kernel_copy_ghost_float_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, NULL);

         ezcl_device_memory_delete(dev_data_buffer_add);

         break;
      case L7_DOUBLE:
         dev_packed_data = ezcl_malloc(NULL, "dev_packed_data", &num_indices_have, sizeof(cl_double), CL_MEM_READ_WRITE, 0);

         ezcl_set_kernel_arg(l7.kernel_pack_double_have_data, 0, sizeof(cl_int), (void *)&num_indices_have);
         ezcl_set_kernel_arg(l7.kernel_pack_double_have_data, 1, sizeof(cl_mem), (void *)&l7_id_db->dev_indices_have);
         ezcl_set_kernel_arg(l7.kernel_pack_double_have_data, 2, sizeof(cl_mem), (void *)&dev_data_buffer);
         ezcl_set_kernel_arg(l7.kernel_pack_double_have_data, 3, sizeof(cl_mem), (void *)&dev_packed_data);

         ezcl_enqueue_ndrange_kernel(command_queue, l7.kernel_pack_double_have_data,   1, NULL, &pack_global_work_size, &pack_local_work_size, NULL);

         packed_double_data = (double *)malloc(num_indices_have*sizeof(double));
         ezcl_enqueue_read_buffer(command_queue, dev_packed_data, CL_TRUE, 0, num_indices_have*sizeof(cl_double), &packed_double_data[0], NULL);
   
         double_data_buffer = (double *)malloc((num_indices_owned+num_indices_needed)*sizeof(double));
         for (unsigned int ii = 0; ii < num_indices_have; ii++){
            double_data_buffer[l7_id_db->indices_have[ii]] = packed_double_data[ii];
         }

         free(packed_double_data);

         ezcl_device_memory_delete(dev_packed_data);

         /*
          * Do the regular L7_Update across processor.
          */
   
         L7_Update (double_data_buffer, l7_datatype, l7_id);

         dev_data_buffer_add    = ezcl_malloc(NULL, "dev_data_buffer_add",    &num_indices_needed,     sizeof(cl_double), CL_MEM_READ_WRITE, 0);
         ezcl_enqueue_write_buffer(command_queue, dev_data_buffer_add, CL_TRUE,  0, num_indices_needed*sizeof(cl_double), &double_data_buffer[num_indices_owned],     NULL);

         free(double_data_buffer);

         // Fill in ghost
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_double_data, 0, sizeof(cl_int), &num_indices_owned);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_double_data, 1, sizeof(cl_int), (void *)&num_indices_needed);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_double_data, 2, sizeof(cl_mem), (void *)&dev_data_buffer);
         ezcl_set_kernel_arg(l7.kernel_copy_ghost_double_data, 3, sizeof(cl_mem), (void *)&dev_data_buffer_add);

         ezcl_enqueue_ndrange_kernel(command_queue, l7.kernel_copy_ghost_double_data,   1, NULL, &ghost_global_work_size, &ghost_local_work_size, NULL);

         ezcl_device_memory_delete(dev_data_buffer_add);

         break;
      default:
         break;
   }
   
#endif /* HAVE_MPI */
   
   return(L7_OK);
    
} /* End L7_Update */

void L7_DEV_UPDATE(
      cl_mem                  dev_data_buffer,
      const enum L7_Datatype  *l7_datatype,
      const int               *l7_id
      )
{

    L7_Dev_Update(dev_data_buffer, *l7_datatype, *l7_id);
}
#endif

#ifdef XXX
int L7_Get_Num_Indices(const int l7_id)
{
   int ierr;

   l7_id_database
     *l7_id_db;            /* database associated with l7_id.    */

   if (l7_id <= 0){
      ierr = -1;
      L7_ASSERT( l7_id > 0, "l7_id <= 0", ierr);
   }

   l7_id_db = l7p_set_database(l7_id);
   if (l7_id_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_id_db != NULL, "Failed to find database.", ierr);
   }

   int num_indices = 0;

   int num_sends = l7_id_db->num_sends;
   
   for (int i=0; i<num_sends; i++){
      /* Load data to be sent. */
      
      num_indices += l7_id_db->send_counts[i];
   }

   return(num_indices);
}

int L7_Get_Local_Indices(const int l7_id, int *local_indices)
{
   int ierr;

   l7_id_database
     *l7_id_db;            /* database associated with l7_id.    */

   if (l7_id <= 0){
      ierr = -1;
      L7_ASSERT( l7_id > 0, "l7_id <= 0", ierr);
   }

   l7_id_db = l7p_set_database(l7_id);
   if (l7_id_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_id_db != NULL, "Failed to find database.", ierr);
   }

   //int num_indices = 0;

   int num_sends = l7_id_db->num_sends;
   
   int offset = 0;

   for (int i=0; i<num_sends; i++){
      /* Load data to be sent. */
      
      int send_count = l7_id_db->send_counts[i];
      for (int j=0; j<send_count; j++){
          local_indices[offset] = l7_id_db->indices_local_to_send[offset];
          offset++;
      }
   }

   return(0);
}
#endif
