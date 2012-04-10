/* This define causes this routine to have the storage for
 * for global variables here and declared extern everywhere
 * else
 */
#define L7_EXTERN
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_INIT"

int L7_Init (
    int *mype,
    int *numpes,
    int *argc,
    char **argv
    )
{
    /*
     * Purpose
     * =======
     * L7_INIT initializes the communication environment (if it has
     * not yet been initialized) and initializes some internal variables.
     *
     * Arguments
     * =========
     * mype         (output) int*
     * numpes       (output) int*
     * argc         (input/output) int* from main
     * argv         (input/output) char** from main
     * 
     * Return Value
     * ============
     * Returns zero if successful and non-zero for error
     * 
     * Notes
     * =====
     * 1) If MPI has not been initialized when this subroutine is called,
     *    L7 will do so. In this case, L7 will also take responsibility for
     *    terminating MPI when L7_TERMINATE is called.
     * 
     */
    
    /*
     * Local variables
     */
    
   int ierr;

#if defined(HAVE_MPI)
   
    int flag;    /* MPI_Initialized input. */
      
    /*
     * Executable Statements
     */
      
    if ( l7.initialized != 0 ) {
       ierr = -1;
       L7_ASSERT( l7.initialized == 0, "L7 already initialized", -1 );
    }
      
    ierr = MPI_Initialized ( &flag );
    L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Initialized", ierr );
      
    if ( !flag && *numpes != -1){
        ierr = MPI_Init(argc, &argv);
        L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Init", ierr);
          
          l7.initialized_mpi = 1;
          l7.mpi_initialized = 1;
    }
    else {
       l7.initialized_mpi = 0;
       l7.mpi_initialized = 0;
    }
      
    if (*numpes != -1) {
       L7_mpi_comm_world_i4 = MPI_COMM_WORLD;
          
       ierr = MPI_Comm_dup (L7_mpi_comm_world_i4, &l7.mpi_comm );
        L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_dup", ierr );
          
        ierr = MPI_Comm_rank (l7.mpi_comm, &l7.penum );
        L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank", ierr );

        ierr = MPI_Comm_size (l7.mpi_comm, &l7.numpes );
        L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_size", ierr );
        *mype = l7.penum;
        *numpes = l7.numpes;
    }
    else {
        l7.penum = 0;
        l7.numpes = 1;
        *mype = 0;
        *numpes = 1;
    }
      
   l7.sizeof_workspace = 0;
    
   l7.sizeof_send_buffer = 0;
    
   l7.initialized = 1;

#else

   *mype = 0;
   *numpes = 1;

#endif /* HAVE_MPI */

   return(ierr);     
}
