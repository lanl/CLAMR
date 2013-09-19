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
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_UPDATE"

int L7_Update(
      void                    *data_buffer,
      const enum L7_Datatype  l7_datatype,
      const int               l7_id
      )
{
   /*
    * Purpose
    * =======
    * L7_Update collects into array data_buffer data located off-process,
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
   
   char
     *pc;                  /* (char *)data_buffer                */
   
   int
     i, j,                 /* Counters                           */
     ierr,                 /* Error code for return              */
     msg_bytes,            /* Message length in bytes.           */
     num_recvs,
     num_outstanding_reqs, /* Outstanding MPI_Requests           */
     num_sends,
     offset,               /* Offset into buffer space           */
     send_count,
     sizeof_type,          /* Number of bytes for input datatype */
     start_index;
   
   int
     *pintdata_buffer,     /* (int *)data_buffer                 */
     *pintsend_buffer;     /* (int *)send_buffer                 */
   
   long long 
     *plongdata_buffer,    /* (long long *)data_buffer           */
     *plongsend_buffer;    /* (long long *)send_buffer           */
   
   float
     *pfloatdata_buffer,   /* (float *)data_buffer               */
     *pfloatsend_buffer;   /* (float *)send_buffer               */
   
   double
     *pdoubledata_buffer,  /* (double *)data_buffer              */
     *pdoublesend_buffer;  /* (double *)send_buffer              */
     
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
   
   if (data_buffer == NULL){
      ierr = -1;
      L7_ASSERT( data_buffer != NULL, "data_buffer != NULL", ierr);
   }
   
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
    * Set some parameters base on input datatype.
    */
   
   sizeof_type = l7p_sizeof(l7_datatype);
   
   /*
    * Receive data into user provided array.
    */
   
   num_outstanding_reqs = 0;
   pc = (char *)data_buffer;
   
   offset = l7_id_db->num_indices_owned * sizeof_type;
   
   num_recvs = l7_id_db->num_recvs;
   
   for (i=0; i<num_recvs; i++){
      msg_bytes = l7_id_db->recv_counts[i] * sizeof_type;
      
#if defined _L7_DEBUG
      printf("[pe %d] Recv posted: pc[%d]; len=%d bytes, from %d \n",
            l7.penum, offset, msg_bytes, l7_id_db->recv_from[i] );
#endif
      
      ierr = MPI_Irecv (&pc[offset], msg_bytes, MPI_BYTE,
            l7_id_db->recv_from[i], l7_id_db->this_tag_update,
            MPI_COMM_WORLD, &l7_id_db->mpi_request[num_outstanding_reqs++] );
      L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Irecv failure", ierr);
      
      offset += l7_id_db->recv_counts[i]*sizeof_type;
   }
   
   /*
    * Send data to processes.
    * (Buffer space allocated in L7_Setup.
    */
   
   switch (l7_datatype){
      case L7_INTEGER4:
      case L7_INT:
      case L7_LOGICAL:
         pintdata_buffer = (int *)data_buffer;
         pintsend_buffer = (int *)l7.send_buffer;
         
         offset = 0;
         start_index = 0;
         
         num_sends = l7_id_db->num_sends;
         
         for (i=0; i<num_sends; i++){
            /* Load data to be sent. */
            
            send_count = l7_id_db->send_counts[i];
            for (j=0; j<send_count; j++){
               pintsend_buffer[offset] =
                  pintdata_buffer[l7_id_db->indices_local_to_send[offset]];
               offset++;
            }
            msg_bytes = l7_id_db->send_counts[i] * sizeof_type;
            
#if defined _L7_DEBUG
            printf("[pe %d] Send pisend_buffer[%d], len=%d ints to %d \n",
                  l7.penum, offset, l7_id_db->send_counts[i], l7_id_db->send_to[i] );
#endif
            
            ierr = MPI_Isend(&pintsend_buffer[start_index], msg_bytes, MPI_BYTE,
                  l7_id_db->send_to[i], l7_id_db->this_tag_update,
                  MPI_COMM_WORLD, &l7_id_db->mpi_request[num_outstanding_reqs++] );
            L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Isend failure", ierr);
            
            start_index += send_count;
         }
         break;
      case L7_INTEGER8:
      case L7_LONG_LONG_INT:
         plongdata_buffer = (long long *)data_buffer;
         plongsend_buffer = (long long *)l7.send_buffer;
         
         offset = 0;
         start_index = 0;
         
         num_sends = l7_id_db->num_sends;
         
         for (i=0; i<num_sends; i++){
            /* Load data to be sent. */
            
            send_count = l7_id_db->send_counts[i];
            for (j=0; j<send_count; j++){
               plongsend_buffer[offset] =
                  plongdata_buffer[l7_id_db->indices_local_to_send[offset]];
               offset++;
            }
            msg_bytes = l7_id_db->send_counts[i] * sizeof_type;
            
#if defined _L7_DEBUG
            printf("[pe %d] Send pisend_buffer[%d], len=%d longs to %d \n",
                  l7.penum, offset, l7_id_db->send_counts[i], l7_id_db->send_to[i] );
#endif
            
            ierr = MPI_Isend(&plongsend_buffer[start_index], msg_bytes, MPI_BYTE,
                  l7_id_db->send_to[i], l7_id_db->this_tag_update,
                  MPI_COMM_WORLD, &l7_id_db->mpi_request[num_outstanding_reqs++] );
            L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Isend failure", ierr);
            
            start_index += send_count;
         }
         break;
      case L7_REAL4:
      case L7_FLOAT:
         pfloatdata_buffer = (float *)data_buffer;
         pfloatsend_buffer = (float *)l7.send_buffer;
         
         offset = 0;
         start_index = 0;
         
         num_sends = l7_id_db->num_sends;
         
         for (i=0; i<num_sends; i++){
            /* Load data to be sent. */
            
            send_count = l7_id_db->send_counts[i];
            for (j=0; j<send_count; j++){
               pfloatsend_buffer[offset] =
                  pfloatdata_buffer[l7_id_db->indices_local_to_send[offset]];
               offset++;
            }
            msg_bytes = l7_id_db->send_counts[i] * sizeof_type;
            
#if defined _L7_DEBUG
            printf("[pe %d] Send pisend_buffer[%d], len=%d floats to %d \n",
                  l7.penum, offset, l7_id_db->send_counts[i], l7_id_db->send_to[i] );
#endif
            
            ierr = MPI_Isend(&pfloatsend_buffer[start_index], msg_bytes, MPI_BYTE,
                  l7_id_db->send_to[i], l7_id_db->this_tag_update,
                  MPI_COMM_WORLD, &l7_id_db->mpi_request[num_outstanding_reqs++] );
            L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Isend failure", ierr);
            
            start_index += send_count;
         }
         break;
      case L7_REAL8:
      case L7_DOUBLE:
         pdoubledata_buffer = (double *)data_buffer;
         pdoublesend_buffer = (double *)l7.send_buffer;
         
         offset = 0;
         start_index = 0;
         
         num_sends = l7_id_db->num_sends;
         
         for (i=0; i<num_sends; i++){
            /* Load data to be sent. */
            
            send_count = l7_id_db->send_counts[i];
            for (j=0; j<send_count; j++){
               pdoublesend_buffer[offset] =
                  pdoubledata_buffer[l7_id_db->indices_local_to_send[offset]];
               offset++;
            }
            msg_bytes = l7_id_db->send_counts[i] * sizeof_type;
            
#if defined _L7_DEBUG
            printf("[pe %d] Send pisend_buffer[%d], len=%d doubles to %d \n",
                  l7.penum, offset, l7_id_db->send_counts[i], l7_id_db->send_to[i] );
#endif
            
            ierr = MPI_Isend(&pdoublesend_buffer[start_index], msg_bytes, MPI_BYTE,
                  l7_id_db->send_to[i], l7_id_db->this_tag_update,
                  MPI_COMM_WORLD, &l7_id_db->mpi_request[num_outstanding_reqs++] );
            L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Isend failure", ierr);
            
            start_index += send_count;
         }
         break;
      default:
         ierr = -1;
         L7_ASSERT(ierr == 0, "Unknown datatype", ierr);
         break;
   } /* End switch ( l7_datatype ) for sending data */
   
   /*
    * Complete all message passing
    */
   
#if defined _L7_DEBUG
   fflush(stdout);
   
   ierr = MPI_Barrier(MPI_COMM_WORLD);
   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Barrier failure", ierr);
   
   for (i=0; i<l7_id_db->numpes; i++){
      if (l7.penum == i){
         printf("-----------------------------------------------------\n");
         printf("Comm for pe %d: num_outstanding_reqs = %d \n",
               l7.penum, num_outstanding_reqs);
         for (j=0; j<l7_id_db->num_sends; j++){
            printf("[pe %d] Send to pe %d. \n", l7.penum, l7_id_db->send_to[j] );
         }
         offset = l7_id_db->num_indices_owned;
         num_recvs = l7_id_db->num_recvs;
      
         for (j=0; j<num_recvs; j++){
            printf("[pe %d] Recving rom pe %d. \n",l7.penum, l7_id_db->recv_from[j] );
            if (l7_datatype == L7_INT){
               printf("[pe %d] Recv complete: pintdata_buffer[%d]=%d; len=%d bytes, from %d \n",
                     l7.penum, offset, (int)pintdata_buffer[offset], msg_bytes,
                     l7_id_db->recv_from[j] );
               offset += l7_id_db->recv_counts[j];
            }
         }
         printf("-----------------------------------------------------\n");
         fflush(stdout);
      }
      sleep(1);
   }
#endif /* _L7_DEBUG */
   
   if (num_outstanding_reqs > 0){
      ierr = MPI_Waitall(num_outstanding_reqs,
            l7_id_db->mpi_request, l7_id_db->mpi_status );
      L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Waitall failure", ierr);
   }
   
   num_outstanding_reqs = 0;
   
   /*
    * Message tag management
    */
   
   l7_id_db->this_tag_update++;
   
   if (l7_id_db->this_tag_update > L7_UPDATE_TAGS_MAX)
      l7_id_db->this_tag_update = L7_UPDATE_TAGS_MIN;
   
#endif /* HAVE_MPI */
   
   return(L7_OK);
    
} /* End L7_Update */

void L7_UPDATE(
      void                    *data_buffer,
      const enum L7_Datatype  *l7_datatype,
      const int               *l7_id
      )
{

    L7_Update(data_buffer, *l7_datatype, *l7_id);
}

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
      /* count data to be sent. */
      
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
      /* Copy data to be sent. */
      
      int send_count = l7_id_db->send_counts[i];
      for (int j=0; j<send_count; j++){
          local_indices[offset] = l7_id_db->indices_local_to_send[offset];
          offset++;
      }
   }

   return(0);
}
