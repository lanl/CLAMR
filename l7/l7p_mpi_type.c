#include "l7.h"
#include "l7p.h"

MPI_Datatype l7p_mpi_type (
		const enum L7_Datatype  l7_datatype
		)
{
	/*
	 * Purpose
	 * =======
	 * l7p_mpi_type returns the name for the input operation as
	 * is required by the underlying communication protocol.
	 * 
	 * Arguments
	 * =========
	 * l7_datatype      (input) const int*
	 *                   The L7 datatype.
	 * 
	 *                   Valid types are:
	 * 
	 *                     L7_GENERIC8
	 *                     L7_BYTE
	 *                     L7_PACKED
	 * 
	 *                     L7_CHAR
	 *                     L7_INT
	 *                     L7_LONG
	 *                     L7_LONG_LONG_INT
	 *                     L7_FLOAT
	 *                     L7_DOUBLE
	 * 
	 *                     L7_CHARACTER
	 *                     L7_LOGICAL
	 *                     L7_INTEGER4
	 *                     L7_INTEGER8
	 *                     L7_REAL4
	 *                     L7_REAL8
	 * 
	 * Return value
	 * ============
	 * Values less than 0 indicates an error.
	 * 
	 * Notes:
	 * =====
	 * 1) The list of datatypes is extensible to include any MPI data type.
	 * 2) Serial compilation (! HAVE_MPI) creates a no-op.
	 * 
	 */
	
	/*
	 * Local variables.
	 */
	
	MPI_Datatype
	  mpi_type;   /* The MPI parameter identifying the operation. */
	
#if defined (HAVE_MPI)
	
	/*
	 * Executable Statements
	 */

	
	switch ( l7_datatype ){

	/* Generic types */
	case L7_GENERIC8:
	case L7_BYTE:
		mpi_type = MPI_BYTE;
		break;
	case L7_PACKED:
		mpi_type = MPI_PACKED;
		break;
	
	/* C Types */
	case L7_CHAR:
		mpi_type = MPI_CHAR;
		break;
	case L7_INT:
		mpi_type = MPI_INT;
		break;
	case L7_LONG:
		mpi_type = MPI_LONG;
		break;
	case L7_LONG_LONG_INT:
		mpi_type = MPI_LONG_LONG_INT;
		break;
	case L7_FLOAT:
		mpi_type = MPI_FLOAT;
		break;
	case L7_DOUBLE:
		mpi_type = MPI_DOUBLE;
		break;

	/* Fortran Types */
	case L7_CHARACTER:
		mpi_type = MPI_CHARACTER;
		break;
	case L7_LOGICAL:
		mpi_type = MPI_LOGICAL;
		break;
	case L7_INTEGER4:
		mpi_type = MPI_INTEGER;
		break;
	case L7_INTEGER8:
		mpi_type = MPI_LONG_LONG_INT;
		break;
	case L7_REAL4:
		mpi_type = MPI_FLOAT;
		break;
	case L7_REAL8:
		mpi_type = MPI_DOUBLE_PRECISION;
		break;

	/* Default is a byte */
	default:
		mpi_type = MPI_BYTE;
		break;
	}
	
#else
		
    mpi_type = 0;
    
#endif /* HAVE_MPI */
    
	return mpi_type;

} /* End l7p_mpi_type. */
