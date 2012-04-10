#include <stdlib.h>
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_TERMINATE"

int L7_Terminate (void)
{
	/*
	 * Purpose
	 * =======
	 * L7_Terminate deallocates its workspace, then, if L7 initialized
	 * MPI on behalf of the application program, terminates the MPI
	 * environment.
	 * 
	 * Arguments
	 * =========
	 * none
	 * 
	 * Return value
	 * ============
	 * Returns zero if successful, non-zero for error
	 * 
	 * Notes
	 * =====
	 * 1) Serial compilation (! HAVE_MPI) creates a no-op.
	 * 
	 */
	
#if defined(HAVE_MPI)
	
	/*
	 * Local variables
	 */
	
	int
	  flag,     /* MPI_Finalized input.       */        
	  ierr;     /* Error code to be returned. */
	
	/*
	 * Executable Statements
	 */
		    
	MPI_Initialized ( &flag  );
	if ( flag )
	{
		ierr = MPI_Comm_rank ( L7_mpi_comm_world_i4, &l7.penum );
		L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank error", ierr );
	}
	else{
		return(0);
	}
	
	if ( l7.initialized != 1) {
		ierr = -1;
		L7_ASSERT( l7.initialized == 1, "L7 not initialized", ierr );
	}
	
	if ( l7.initialized_mpi == 1 ){
		ierr = MPI_Finalized ( &flag );
		if ( !flag ){
			ierr = MPI_Finalize ();
			L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Finalize", ierr );
		}
		l7.initialized_mpi = 0;
	}
	
	if ( l7.send_buffer != NULL ){
		free ( l7.send_buffer);
		l7.sizeof_send_buffer = 0;
	}
	
	l7.initialized = 0;
	
#endif /* HAVE_MPI */
	
	return(0);

} /* End L7_Terminate */

void L7_TERMINATE(int *ierr)
{
	*ierr = L7_Terminate();
} /* End l7_terminate_ */


