#include "l7.h"
#include "l7p.h"

#include <stdlib.h>

l7_id_database *l7p_set_database(
      int l7_id
      )
{
   /*
    * Purpose
    * =======
    * l7p_set_database returns a pointer to the L7 database associated
    * with the input handle.
    * 
    * Arguments
    * =========
    * l7_id           (input) const int*
    *                 Handle to database
    * 
    * Notes:
    * ======
    * 1) Serial compilation creates a no-op.
    * 
    */
   
#if defined HAVE_MPI
   
   /*
    * Local Declarations
    */
   
   struct l7_id_database
     *l7_id_db;    /* For searching through list of databases */
   
   /*
    * Executable statements
    */
   
   /*
    * Search initialized databases for input db handle.
    */
   
   l7_id_db = l7.first_db;
   
   while (l7_id_db){
      if (l7_id_db->l7_id == l7_id){
         return(l7_id_db);
      }
      else{
         /* Move to next one */
         l7_id_db = l7_id_db->next_db;
      }
   }

#endif /* HAVE_MPI */
   
   return(NULL);
   
} /* End l7_id_database */
