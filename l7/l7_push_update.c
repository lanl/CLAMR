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

#define L7_LOCATION "L7_PUSH_UPDATE"

int L7_Push_Update(
      const int               *array,
      int                     *return_array,
      const int               l7_push_id
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
    * array              (input) data input
    *
    * return_array       (output) data input
    * 
    * l7_push_id         (input) const int
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
     ierr;                        /* Error code for return                */

   l7_push_id_database
     *l7_push_id_db;
   
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

   if (array == NULL){
      ierr = -1;
      L7_ASSERT( array != NULL, "array != NULL", ierr);
   }
   
   if (l7_push_id <= 0){
      ierr = -1;
      L7_ASSERT( l7_push_id > 0, "l7_push_id <= 0", ierr);
   }
   
   if (l7.numpes == 1){
      ierr = L7_OK;
      return(ierr);
   }
   
   /*
    * find database associated with input l7_id
    */
   
   l7_push_id_db = l7.first_push_db;
   
   while (l7_push_id_db){
      if (l7_push_id_db->l7_push_id == l7_push_id){
            break;
      } else {
         /* Move to next one */
         l7_push_id_db = l7_push_id_db->next_push_db;
      }
   }

   if (l7_push_id_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_push_id_db != NULL, "Failed to find database.", ierr);
   }
   
   if (l7.numpes == 1){ /* No-op */
      ierr = L7_OK;
      return(ierr);
   }
   
   /*
    * Set some parameters base on input datatype.
    */

   for (int ip = 0; ip < l7_push_id_db->num_comm_partners; ip++){
      int count = l7_push_id_db->send_buffer_count[ip]; // for vectorization
      for (int ic = 0; ic < count; ic++){
         l7_push_id_db->send_buffer[ip][ic] = array[l7_push_id_db->send_database[ip][ic]];
      }    
   }    


// Send/Receives will be done in L7_Push_Update. Input will be send_buffer. Output will be in
// preallocated receive_buffer
   MPI_Request request[2*l7_push_id_db->num_comm_partners];
   MPI_Status  status[2*l7_push_id_db->num_comm_partners];

   int iloc = 0;
   for (int ip = 0; ip < l7_push_id_db->num_comm_partners; ip++){
      MPI_Irecv(&return_array[iloc], l7_push_id_db->recv_buffer_count[ip], MPI_INT,
                l7_push_id_db->comm_partner[ip], l7_push_id_db->comm_partner[ip], MPI_COMM_WORLD, &request[ip]);
      iloc += l7_push_id_db->recv_buffer_count[ip];
   }

   for (int ip = 0; ip < l7_push_id_db->num_comm_partners; ip++){
      MPI_Isend(l7_push_id_db->send_buffer[ip], l7_push_id_db->send_buffer_count[ip], MPI_INT,
                l7_push_id_db->comm_partner[ip], l7.penum, MPI_COMM_WORLD, &request[l7_push_id_db->num_comm_partners+ip]);
   }    
   MPI_Waitall(2*l7_push_id_db->num_comm_partners, request, status);

/*
   if (ncycle >= 1) { 
      for (int ib = 0; ib<receive_count_total; ib++){
         fprintf(fp,"DEBUG receive %d is %d\n",ib,border_data_receive[ib]);
      }    
   }    
*/
   
#endif /* HAVE_MPI */
   
   return(L7_OK);
    
} /* End L7_Update */

