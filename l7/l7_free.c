#include <stdlib.h>
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_FREE"

int L7_Free(
      const int *l7_id
      )
{
   /*
    * Purpose
    * =======
    * L7_Free destroys the database associated with the input
    * L7 handle. All memory assoicated with this database
    * is deallocated.
    * 
    * Arguments
    * =========
    * l7_id          (input) const int*
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
     ierr;            /* Error code for return              */
   
#if defined HAVE_MPI
   
   l7_id_database
     *l7_db,          /* Database for the input l7_id.      */
     *prev_l7_db;     /* Previous structure in linked list. */
   
   /*
    * Executable Statements
    */
   
   if (! l7.mpi_initialized){
      return(0);
   }
   
   ierr = MPI_Comm_rank(l7.mpi_comm, &l7.penum);
   L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank error", ierr);
 
   if (l7.initialized != 1){
      ierr = -1;
      L7_ASSERT( l7.initialized != 1, "L7 not initialized", ierr);
   }
   
   /*
    *  Get first L7 structure and null out previous structure
    */
   
   l7_db      = l7.first_db;
   prev_l7_db = (l7_id_database *)0;
   
   /*
    * Find the structure to be freed
    */
   
   while (l7_db){
      if (l7_db->l7_id == *l7_id )
         break;
      prev_l7_db = l7_db;
      l7_db = l7_db->next_db;
   }
   
   if (l7_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_db != NULL, "Failed to find database.", ierr);
   }
   
   /*
    * Free all data associated with this id.
    */
   
   if (l7_db->indices_needed)
      free(l7_db->indices_needed);
   
   if (l7_db->recv_from)
      free(l7_db->recv_from);
   
   if (l7_db->recv_counts)
      free(l7_db->recv_counts);
   
   if (l7_db->send_to)
      free(l7_db->send_to);
   
   if (l7_db->send_counts)
      free(l7_db->send_counts);
   
   if (l7_db->indices_global_to_send)
      free(l7_db->indices_global_to_send);
   
   if (l7_db->indices_local_to_send)
      free(l7_db->indices_local_to_send);
   
   if (l7_db->starting_indices)
      free(l7_db->starting_indices);
   
   if (l7_db->mpi_request)
      free(l7_db->mpi_request);
   
   if (l7_db->mpi_status)
      free(l7_db->mpi_status);
   
   if (&l7_db->mpi_comm_l7)
      MPI_Comm_free(&l7_db->mpi_comm_l7);
   
   /*
    * Assign pointers to next, first, and last if needed
    */
   
   if (l7.first_db == l7_db && l7.last_db == l7_db){
      /*
       * Only one database currently stored.
       */
      
      l7.first_db = (l7_id_database *)0;
      l7.last_db  = (l7_id_database *)0;
   }
   else if (l7.first_db == l7_db){
      /*
       * If this id is the first id, reassign first to next
       */
      
      l7.first_db = l7_db->next_db;
   }
   else if (l7.last_db == l7_db){
      /*
       * If this is the last id, reassign last to prev
       */
      
      l7.last_db = prev_l7_db;
      prev_l7_db->next_db = (l7_id_database *)0;
   }
   else {
      /*
       *  if this id is an interior id
       */
      prev_l7_db->next_db = l7_db->next_db;
   }
   
   /*
    * Free the database
    */
   
   free(l7_db);
   
   l7_db = NULL;
   
   /*
    * Decrement number of databases stored
    */
   
   l7.num_dbs--;
   
#endif /* HAVE_MPI */

   ierr = L7_OK;
   
   return(ierr);
   
}

void L7_FREE(
        const int   *l7_id,
        int         *ierr
        )
{
   L7_Free(l7_id);
   ierr = 0;
}


