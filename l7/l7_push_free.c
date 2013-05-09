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

#define L7_LOCATION "L7_PUSH_FREE"

int L7_Push_Free(
      const int *l7_push_id
      )
{
   /*
    * Purpose
    * =======
    * L7_Push_Free destroys the database associated with the input
    * L7 handle. All memory assoicated with this database
    * is deallocated.
    * 
    * Arguments
    * =========
    * l7_push_id     (input) const int*
    *                Handle to database to be destroyed.
    * 
    * 
    * Return value
    * ============
    * Value other than L7_OK indicates an error.
    * 
    * Notes:
    * =====
    * It is not necessary to destroy databases once they are not needed.
    * However, a significant amount of memory use can build up if 
    * databases are created throuhout execution. Therefore, the user
    * is strongly encouraged to destroy databases that are no longer
    * needed.
    * 
    */
   
   /*
    * Local variables.
    */
   
   int
     ierr;              /* Error code for return              */
   
#if defined HAVE_MPI
   
   l7_push_id_database
     *l7_push_db,       /* Database for the input l7_id.      */
     *prev_l7_push_db;  /* Previous structure in linked list. */
   
   /*
    * Executable Statements
    */
   
   if (! l7.mpi_initialized){
      return(0);
   }
   
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &l7.penum);
   L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank error", ierr);
 
   if (l7.initialized != 1){
      ierr = -1;
      L7_ASSERT( l7.initialized != 1, "L7 not initialized", ierr);
   }
   
   /*
    *  Get first L7 structure and null out previous structure
    */
   
   l7_push_db      = l7.first_push_db;
   prev_l7_push_db = (l7_push_id_database *)0;
   
   /*
    * Find the structure to be freed
    */
   
   while (l7_push_db){
      if (l7_push_db->l7_push_id == *l7_push_id )
         break;
      prev_l7_push_db = l7_push_db;
      l7_push_db = l7_push_db->next_push_db;
   }
   
   if (l7_push_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_push_db != NULL, "Failed to find database.", ierr);
   }
   
   /*
    * Free all data associated with this id.
    */
   
   if (l7_push_db->comm_partner){
      free(l7_push_db->comm_partner);
      l7_push_db->comm_partner = NULL;
   }
   
   if (l7_push_db->send_buffer_count){
      free(l7_push_db->send_buffer_count);
      l7_push_db->send_buffer_count = NULL;
   }
   
   if (l7_push_db->recv_buffer_count){
      free(l7_push_db->recv_buffer_count);
      l7_push_db->recv_buffer_count = NULL;
   }
   
   if (l7_push_db->send_database){
      for (int ip = 0; ip < l7_push_db->num_comm_partners; ip++){
         free(l7_push_db->send_database[ip]);
      }
      free(l7_push_db->send_database);
      l7_push_db->send_database = NULL;
   }

   if (l7_push_db->send_buffer){
      for (int ip = 0; ip < l7_push_db->num_comm_partners; ip++){
         free(l7_push_db->send_buffer[ip]);
      }
      free(l7_push_db->send_buffer);
      l7_push_db->send_buffer = NULL;
   }

   /*
    * Assign pointers to next, first, and last if needed
    */
   
   if (l7.first_push_db == l7_push_db && l7.last_push_db == l7_push_db){
      /*
       * Only one database currently stored.
       */
      
      l7.first_push_db = (l7_push_id_database *)0;
      l7.last_push_db  = (l7_push_id_database *)0;
   }
   else if (l7.first_push_db == l7_push_db){
      /*
       * If this id is the first id, reassign first to next
       */
      
      l7.first_push_db = l7_push_db->next_push_db;
   }
   else if (l7.last_push_db == l7_push_db){
      /*
       * If this is the last id, reassign last to prev
       */
      
      l7.last_push_db = prev_l7_push_db;
      prev_l7_push_db->next_push_db = (l7_push_id_database *)0;
   }
   else {
      /*
       *  if this id is an interior id
       */
      prev_l7_push_db->next_push_db = l7_push_db->next_push_db;
   }
   
   /*
    * Free the database
    */
   
   free(l7_push_db);
   
   l7_push_db = NULL;
   
   /*
    * Decrement number of databases stored
    */
   
   l7.num_push_dbs--;
   
#endif /* HAVE_MPI */

   ierr = L7_OK;
   
   return(ierr);
   
}
