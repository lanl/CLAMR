#include "l7.h"
#include "l7p.h"

int l7p_sizeof(
		const enum L7_Datatype  l7_datatype
		)
{
	/*
	 * Purpose
	 * =======
	 * l7p_sizeof returns the number of bytes needed to store
	 * one element of the input L7 datatype.
	 * 
	 * Arguments
	 * =========
	 * l7_datatype     (input) const int*
	 *                 The type of data, as defined by L7 datatypes,
	 *                 in databuffer.
	 * 
	 * Return value
	 * ============
	 * Values less than 0 indicates an error.
	 * 
	 * Notes:
	 * ======
	 * 
	 */
	
	/*
	 * Local variables.
	 */
	
	int
	  sizeof_type;  /* Number of bytes in one element of the input datatype. */
	
	/*
	 * Executable Statements
	 */
	
	switch (l7_datatype)
	{
	case L7_GENERIC8:
	case L7_BYTE:
	case L7_PACKED:
	case L7_CHAR:
	case L7_LONG_LONG_INT:
	case L7_DOUBLE:
	case L7_CHARACTER:
	case L7_INTEGER8:
	case L7_REAL8:
		sizeof_type = 8;
		break;
	case L7_INT:
	case L7_FLOAT:
	case L7_LOGICAL:
	case L7_INTEGER4:
	case L7_REAL4:
		sizeof_type = 4;
		break;
	case L7_LONG:
		sizeof_type = sizeof(long);
		break;
	default:
	    sizeof_type = -1;
	    break;
	}
	
	return(sizeof_type);
}
