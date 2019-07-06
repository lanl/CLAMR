/*
 *  Copyright (c) 2011-2019, Triad National Security, LLC.
 *  All rights Reserved.
 *
 *  CLAMR -- LA-CC-11-094
 *
 *  Copyright 2011-2019. Triad National Security, LLC. This software was produced 
 *  under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Triad National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 *  TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
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
 *     * Neither the name of the Triad National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE TRIAD NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL
 *  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *  
 */
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_BROADCAST"

int L7_Broadcast (
		void                    *data_buffer,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		const int               root_pe
		)
{
	/* Purpose
	 * =======
	 * L7_Broadcast is used to broadcast out data from one
	 * processor to all the others.
	 * 
	 * Arguments
	 * =========
	 * data_buffer           (input/output) const L7_PTR*
	 *                       Pointer to data buffer. Must exist on
	 *                       all processors and have initialized
	 *                       data on root processor.
	 * 
	 * count                 (input) const L7_INT*
	 *                       Number of data items to send.
	 * 
	 * l7_datatype           (input) const L7_INT*
	 *                       Datatype to send.
	 * 
	 * root_pe               (input) const L7_INT*
	 *                       Processor that will send data to the
	 *                       other processors. This is usually processor
	 *                       zero.
	 * 
	 * Return value
	 * ============
	 * Returns non-zero value for any error
	 *
	 * Notes:
	 * ======
	 * 
	 * 1) Serial operation is a no-op.
	 * 
	 * Program flow
	 * ============
	 * 0) Check input for basic validity.
	 * 1) Call MPI Broadcast routine.
	 * =================================================================
	 */
	
	/*
	 * Local variables.
	 */
	
   int
      ierr=L7_OK;             /* Temporary storage for error code  */

#ifdef HAVE_MPI
   
	int
	   local_count;      /* Temporary for handling INTEGER8   */

	MPI_Datatype
	   mpi_type;         // Number of bytes for input datatype.
	
	/*
	 * Executable Statements
	 */
	
	if (l7.initialized != 1){
		ierr = -1;
		L7_ASSERT(l7.initialized == 1, "L7 not initialized", ierr);
	}
	
	if (l7.initialized_mpi){
		
		/*
		 * Set some parameters based on input datatype.
		 */
		
		mpi_type = l7p_mpi_type(l7_datatype);
		
/*
		if (mpi_type < 0){
			ierr = -2;
			L7_ASSERT(mpi_type > 0, "l7p_mpi_type: Unknown l7_mpi_type.",
					ierr);
		}
*/
		
		/*
		 * Call MPI Broadcast
		 */
		
		local_count = count;
		ierr = MPI_Bcast(data_buffer, local_count, mpi_type, root_pe,
				MPI_COMM_WORLD);
		
		if (ierr != L7_OK){
			ierr = -3;
			L7_ASSERT( ierr == L7_OK, "L7 Broadcast error", ierr);
		}
	}
#endif /* HAVE_MPI */
	return(ierr);
} /* End L7_BROADCAST */

void l7_broadcast_ (
		void                   *data_buffer,
		const int              *count,
		const enum L7_Datatype *l7_datatype,
		const int              *root_pe,
		int                    *ierr
		)
{
	*ierr=L7_Broadcast(data_buffer, *count, *l7_datatype, *root_pe);
} /* End l7_broadcast_ */


