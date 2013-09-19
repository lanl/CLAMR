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
#include <stdlib.h>

int l7_int_comp(const int *x, const int *y);

#define L7_LOCATION "L7_SETUP"

int L7_Dev_Setup(
		const int      num_base,
		const int      my_start_index,
		const int      num_indices_owned,
		int            *indices_needed,
		const int      num_indices_needed,
		int            *l7_id
		)
{
	/* Purpose
	 * =======
	 * L7_Dev_Setup is used to setup the update/scatter database as
	 * defined by the global indexing scheme. Each process passes
	 * in parameters which define the indices it owns (i.e. as
	 * defined by 'my_start_index' and 'num_indices_owned') and
	 * lists the indices it needs ('indices_needed'). From this,
	 * a database is defined that allows subsequent calls to
	 * L7_Update.
	 * 
	 * Notes:
	 * ======
	 * 1) Assumes a global indexing set, linearly decomposed across
	 * all processes.
	 * 
	 * Arguments
	 * =========
	 * num_base             (input) const L7_INT
	 *                      global indexing set starts with 1 (Fortran)
         *                      or with 0 (C)
         *
	 * my_start_index       (input) const L7_INT
	 *                      Starting index number of calling process
	 *                      in global indexing set.
	 * 
	 * num_indices_owned    (input) const L7_INT
	 *                      Number of indices owned by calling process.
	 * 
	 * indices_needed       (input) const L7_INT*
	 *                      Array containing indices needed by
	 *                      calling process.
	 * 
	 * num_indices_needed   (input) const L7_INT
	 *                      Number of indices of interest listed
	 *                      in array 'num_indices_needed'.
	 * 
	 * l7_id                (input/output) int*
	 *                      Handle to database to be setup.
	 * 
	 *                      0: L7 sets up a new database, and
	 *                      assigns it a value.
	 *                      > 0: L7 resets existing database with
	 *                      input information. That is, it reuses
	 *                      the allocated memory.
	 *                      < 0: An error is returned.
	 * 
	 * Notes:
	 * =====
	 * 1) The handling of 0-based arrays for C and 1-based arrays for Fortran
	 * is handled in L7_Dev_Setup. This is done by taking the input global
	 * indices stored in 'indices_global_to_send' and converting them to 
	 * 1-based and storing them in 'indices_local_to_send'.
	 * 
	 * 2) The indices are handled as 4-byte integers.
	 * 
	 * 3) Serial compilation creates a no-op.
	 * 
	 * Program Flow
	 * ============
	 * 0) Check input for basic validity.
	 * 1) Set communication parameters within database.
	 * 2) Deternine processes this pe receives from.
	 * 3) Determine the number of processes this pe sends to.
	 * 4) Send number of as well as the indices needed from each sending process.
	 * 5) Set up array containing the pes this pe sends indices to.
	 * 6) Set up array containing the indices this pe sends to others.
	 */
	
	/*
	 * Local variables.
	 */

   int
     ierr;                        /* Error code for return                */
   
#ifdef HAVE_MPI
   
	int
          base_adj,                    /* 0 or 1 based arrays adjustment       */
	  count_total,
	  i, j,                        /* Counters                             */
	  max_sizeof_type,
	  num_msgs,                    /* Number of sends and recvs needed     */
	  numpes,                      /* Alias for l7_id_db.numpes.           */
	  num_indices_acctd_for,
	  num_outstanding_requests = 0,
	  num_sends,
	  offset,
	  penum,                       /* Alias for l7_id_db.penum.             */
	  *pi4_in,                     /* (int *)l7.receive_buffer              */
	  *pi4_out,                    /* (int *)l7.send_buffer                 */
	  send_buffer_bytes_needed,    /* Buffer space requirement.             */
	  start_indices_needed,
	  this_index;                  /* Offset into indexing set.             */
	
	l7_id_database
	  *l7_id_db;
	
	MPI_Request
	  *mpi_request;                /* Local alias for l7_id_db->mpi_request. */
	
	MPI_Status
	  *mpi_status;                 /* Local alias for l7_id_db->mpi_status.  */
	

#if defined (_L7_DEBUG)
	
	int
	  k;                           /* Counter                                */
	
#endif

	/*
	 * Executable Statements
	 */
	
	if (! l7.mpi_initialized){
		return(0);
	}
	
   if (l7.initialized != 1){
		ierr = -1;
		L7_ASSERT( l7.initialized == 1, "L7 not initialized", ierr);
	}
	
	/*
	 * Check input
	 */
        if (num_base){
          base_adj = 1;
        }
        else {
          base_adj = 0;
        }
	
	if (my_start_index < 0){
		ierr = -1;
		L7_ASSERT( my_start_index >= 0, "my_start_index < 0", ierr);
	}
	
	if (num_indices_owned < 0){
		ierr = -1;
		L7_ASSERT( num_indices_owned >= 0, "num_indices_owned < 0", ierr);
	}

	if (num_indices_needed > 0){
		if (indices_needed == NULL){
			ierr = -1;
			L7_ASSERT( (int *)indices_needed != NULL,
					"indices_needed == NULL", ierr);
		}
	}
	if (*l7_id < 0){
		ierr = *l7_id;
		L7_ASSERT( *l7_id >=0,
				"L7 Id must be either 0 (new id) or > 0 (existing id)",
				ierr);
	}

	/*
	 * Setup database structure.
	 */
	
	if (*l7_id != 0){
		/*
		 * Find it in the database and update based on new input.
		 */
		
		if (l7.first_db == NULL){
			L7_ASSERT(l7.first_db != NULL,
					"Uninitialized l7_id input, but no ids in database",
					ierr);
		}
		
		l7_id_db = l7.first_db;
		while (l7_id_db){
			if (l7_id_db->l7_id == *l7_id)
				break;
			l7_id_db = l7_id_db->next_db;
		}
		if (l7.first_db == NULL){
			ierr = -1;
			L7_ASSERT( l7.first_db != NULL,
					"Uninitialized l7_id input, but not found in this list",
					ierr);
		}
		
	}
	else{
		
		/*
		 * Allocate new database, insert into linked list.
		 */
		
		if (l7.num_dbs >= L7_MAX_NUM_DBS){
			ierr = -1;
			L7_ASSERT(l7.num_dbs < L7_MAX_NUM_DBS,
					"Too many L7 databases allocataed",
					ierr);
		}
		
		l7_id_db = (l7_id_database*)calloc(1L, sizeof(l7_id_database) );
		
		if (l7_id_db == NULL){
			ierr = -1;
			L7_ASSERT( l7_id_db != NULL, "Failed to allocate new database",
					ierr);
		}
		
		if ( !(l7.first_db) ){
			l7.first_db = l7_id_db;
			l7.last_db  = l7_id_db;
			l7_id_db->next_db = NULL; /* Paranoia */
			
			l7_id_db->l7_id = 1;
			
			l7.num_dbs = 1;
		}
		else{
			
			/*
			 * Assign a l7_id.
			 */
			
			l7_id_db->l7_id = l7.last_db->l7_id + 1;
			
			/*
			 * Reset links.
			 */
			
			l7.last_db->next_db = l7_id_db;
			
			l7.last_db = l7_id_db;
			
			l7.num_dbs++;
		}
		
		*l7_id = l7_id_db->l7_id;
		
		/*
		 * Initialize some parameters.
		 */
		
		l7_id_db->recv_counts_len = 0;
		l7_id_db->recv_from_len   = 0;
		
		l7_id_db->send_to_len     = 0;
		l7_id_db->send_counts_len = 0;
		
		l7_id_db->indices_to_send_len = 0;
		
		l7_id_db->mpi_request_len = 0;
		l7_id_db->mpi_status_len  = 0;
	}
	
	/*
	 *  Store input in database.
	 */
	
	l7_id_db->my_start_index    = my_start_index;
	l7_id_db->num_indices_owned = num_indices_owned;
	
	if ( (l7_id_db->indices_needed_len < num_indices_needed ) &&
		 (num_indices_needed > 0) ){
		if (l7_id_db->indices_needed)
			 free(l7_id_db->indices_needed);
			 
		l7_id_db->indices_needed =
			 (int *)calloc((unsigned long long)num_indices_needed, sizeof(int) );
			 
		if (l7_id_db->indices_needed == NULL){
			 ierr = -1;
			 L7_ASSERT( (int*)(l7_id_db->indices_needed) != NULL,
			            "Memory failure for indices_needed",
			            ierr);
		 }
		 l7_id_db->indices_needed_len = num_indices_needed;
	}

	for (i=0; i<num_indices_needed; i++){
		l7_id_db->indices_needed[i] = indices_needed[i];
	}
    
	l7_id_db->num_indices_needed = num_indices_needed;

	ierr = MPI_Comm_rank (MPI_COMM_WORLD, &l7_id_db->penum );
	L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank", ierr);
	
	ierr = MPI_Comm_size (MPI_COMM_WORLD, &l7_id_db->numpes );
	L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_size", ierr);
	
	l7.penum = l7_id_db->penum;
	
	/* Local shorthand */
	
	numpes   = l7_id_db->numpes;
	penum    = l7_id_db->penum;
	
	if (numpes == 1){
		return(0);
	}
	
	/*
	 * Create array containing starting (global) index numbers
	 * for all processes.
	 * 
	 * 1) Allgather num_indices_owned.
	 * 2) Scan to create starting_index.
	 * 3) Shift all array elements up 1 position.
	 * 4) Set starting_indices[0] = 0.
	 * 
	 * The latter two steps allows arrays to be used as below.
	 */
	
	l7_id_db->starting_indices =
		(int *)calloc((unsigned long long)(numpes+1), sizeof(int));
	if(l7_id_db->starting_indices == NULL){
		ierr = -1;
		L7_ASSERT(l7_id_db->starting_indices != NULL,
				"No memory for l7_id_db->starting_indices", ierr);
	}
	
   ierr = MPI_Allgather( &(l7_id_db->num_indices_owned), 1, MPI_INT,
			&(l7_id_db->starting_indices[1]), 1, MPI_INT,
			MPI_COMM_WORLD);
	L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Allgather (num_indices_owned)",
			ierr);
	
  	l7_id_db->starting_indices[0] = 0;
//	l7_id_db->starting_indices[0] = 1;
	
	for (i=0; i<numpes; i++)
		l7_id_db->starting_indices[i+1] += l7_id_db->starting_indices[i];
	
	/*
	 * Determine the number of processes this pe receives from.
	 */
	
	l7_id_db->num_recvs  =  0;
	start_indices_needed = -1;
	
	this_index = 0;
	if (num_indices_needed > 0){
		for (j=0; j<numpes; j++){
			if ( indices_needed[this_index] >= l7_id_db->starting_indices[j]){
				if (indices_needed[this_index] < l7_id_db->starting_indices[j+1]){
					l7_id_db->num_recvs++;
#if defined _L7_DEBUG
					printf("[pe %d] Found first one on pe %d. \n", penum, j);
#endif
					/* Skip through all the rest on pe j. */
					
					while ( ( indices_needed[this_index] < 
							     l7_id_db->starting_indices[j+1] ) &&
							( this_index < num_indices_needed) )
						this_index++;
					
					/* Remember where we found the first one. */
					
					if ( start_indices_needed == -1)
						start_indices_needed = j;
					
					if (this_index == num_indices_needed)
						break;
				}
			}
		}
		
		if (l7_id_db->num_recvs == 0){
			ierr = -1;
			L7_ASSERT(l7_id_db->num_recvs != 0, "No indices found", ierr);
		}

	}
	
	if (this_index != num_indices_needed){
		printf("[pe %d] ERROR -- can't find all the indices I need. I have %d, need %d\n",
				penum, this_index, num_indices_needed);
	}
	
#if defined _L7_DEBUG
	
	printf("[pe %d] l7_id_dp->num_recvs = %d\n",
			penum, l7_id_db->num_recvs);

#endif
	
	/*
	 * Allocate space for counts for each pe sending to this one.
	 */
	
	if (l7_id_db->num_recvs > l7_id_db->recv_counts_len){
		if (l7_id_db->recv_counts)
			free(l7_id_db->recv_counts);
		
		l7_id_db->recv_counts = 
			(int *)calloc((unsigned long long)l7_id_db->num_recvs, sizeof(int) );
		if (l7_id_db->recv_counts == NULL){
			ierr = -1;
			L7_ASSERT(l7_id_db->recv_counts != NULL,
					"No space for l7_id_db->recv_counts", ierr);
		}
		l7_id_db->recv_counts_len = l7_id_db->num_recvs;
		
                int num_recvs = l7_id_db->num_recvs; // for vectorization
		for (i=0; i<num_recvs; i++)
			l7_id_db->recv_counts[i] = 0; /* calloc does not guarantee = 0. */
				
	}
	
	if (l7_id_db->num_recvs > l7_id_db->recv_from_len){
		if (l7_id_db->recv_from)
			free(l7_id_db->recv_from);
		
		l7_id_db->recv_from =
			(int *)calloc((unsigned long long)l7_id_db->num_recvs, sizeof(int) );
		
	if (l7_id_db->recv_from == NULL){
		ierr = -1;
		L7_ASSERT(l7_id_db->recv_from != NULL,
				"No space for l7_id_db->recv_from", ierr);
	}
	l7_id_db->recv_from_len = l7_id_db->num_recvs;
	
        int num_recvs = l7_id_db->num_recvs; // for vectorization
	for (i=0; i<num_recvs; i++)
		l7_id_db->recv_from[i] = -999;
	}
	
	/*
	 * Determine process and the number of indices this pe recvs from it.
	 */
	
	if (num_indices_needed > 0){
		this_index = 0;
		
		num_indices_acctd_for = 0;
		
		i=0;
		for (j=start_indices_needed; j<numpes; j++){
		   if (indices_needed[this_index] >= l7_id_db->starting_indices[j] ){
		      if (indices_needed[this_index] < l7_id_db->starting_indices[j+1]){
		         /* Found the first one on pe j. */
		         
		         l7_id_db->recv_from[i]   = j;
		         l7_id_db->recv_counts[i] = 1;
		         
		         num_indices_acctd_for++;
		         
		         if (num_indices_acctd_for == num_indices_needed)
		            break;
		         
		         this_index++;
		         
		         while ( ( indices_needed[this_index] < l7_id_db->starting_indices[j+1] ) &&
		               ( num_indices_acctd_for < num_indices_needed ) ){
		            /* Find the rest on pe j. */
		            
		            l7_id_db->recv_counts[i]++;
		            this_index++;
		            num_indices_acctd_for++;
		         }
		         
		         if (num_indices_acctd_for == num_indices_needed)
		            break;
		         
		         i++;
		      }
		   }
		}
		
		if (num_indices_needed != num_indices_acctd_for){
		   ierr = -1;
		   L7_ASSERT(num_indices_needed == num_indices_acctd_for,
		         "Failed to find all the needed indices", ierr);
		}
		
	}
	
	/*
	 * Determine number of processes for which this pe owns indices
	 * those pes need. This is done use a reduction (MPI_Allreduce).
	 */
	
	if (l7.sizeof_send_buffer < numpes * (int)sizeof(int)){
	   if (l7.send_buffer)
	      free(l7.send_buffer);
	   
	   l7.send_buffer = calloc ((unsigned long long)(2*numpes), sizeof(int));
	   if (l7.send_buffer == NULL){
	      ierr = -1;
	      L7_ASSERT(l7.send_buffer != NULL, "No memory for send buffer", ierr);
	   }
	   
	   l7.sizeof_send_buffer = 2 * numpes * (int)sizeof(int);
	}
	
	pi4_in   = (int*)l7.send_buffer;
	pi4_out  = &pi4_in[numpes];
	
	for (i=0; i<numpes; i++)
	   pi4_in[i] = 0;
	
	for (i=0; i<l7_id_db->num_recvs; i++)
	   pi4_in[l7_id_db->recv_from[i]] = 1;
	
	ierr = MPI_Allreduce(pi4_in, pi4_out, numpes, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Allreduce ( l7_id_db->recv_from )", ierr);
	
	l7_id_db->num_sends = pi4_out[penum];
	
#if defined _L7_DEBUG
	printf("[pe %d] l7_id_db->num_sends = %d \n", penum, l7_id_db->num_sends);
#endif
	
	/*
	 * Allocate request and status arrays.
	 */
	
	num_msgs = ( 2 * l7_id_db->num_recvs ) + l7_id_db->num_sends;
	
	/* Ensure enough outstanding messages for L7_Update_pack model. */
	
	if (num_msgs < (L7_MIN_MPI_REQS * l7_id_db->num_recvs ) )
	   num_msgs = L7_MIN_MPI_REQS * l7_id_db->num_recvs;
	
	if (num_msgs > l7_id_db->mpi_request_len) {
	   if (l7_id_db->mpi_request)
	      free(l7_id_db->mpi_request);
	   
	   l7_id_db->mpi_request = (MPI_Request *) calloc ((unsigned long long)num_msgs, sizeof(MPI_Request));
	   
	  if (l7_id_db->mpi_request == NULL){
	     ierr = -1;
	     L7_ASSERT(l7_id_db->mpi_request != NULL,
	           "Allocation of l7_id_db->mpi_request failed", ierr);
	  }
	  l7_id_db->mpi_request_len = num_msgs;
	}
	
	if (num_msgs > l7_id_db->mpi_status_len){
	   if (l7_id_db->mpi_status)
	      free(l7_id_db->mpi_status);
	   
	   l7_id_db->mpi_status = (MPI_Status *) calloc((unsigned long long)num_msgs, sizeof(MPI_Status) );
	   if (l7_id_db->mpi_status == NULL){
	      ierr = -1;
	      L7_ASSERT(l7_id_db->mpi_status != NULL,
	            "Allocation of l7_id_db->mpi_status failed", ierr);
	   }
	   l7_id_db->mpi_status_len = num_msgs;
	}
	
	/* Local shorthand */
	
	mpi_request = l7_id_db->mpi_request;
	mpi_status  = l7_id_db->mpi_status;
	
	/*
	 * Send number of indices needed from each sending process.
	 */
	
	num_outstanding_requests = 0;
	for (i=0; i<l7_id_db->num_recvs; i++){
#if defined _L7_DEBUG
	   printf("[pe %d] recv_counts[%d] = %d to pe %d  \n", penum, i,
	         l7_id_db->recv_counts[i], l7_id_db->recv_from[i] );
#endif
	   
	   ierr = MPI_Isend(&l7_id_db->recv_counts[i], 1, MPI_INT,
	         l7_id_db->recv_from[i], L7_SETUP_SEND_COUNT_TAG,
	         MPI_COMM_WORLD, &mpi_request[num_outstanding_requests++] );
	   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Isend (recv_counts[i] )",
	         ierr);
	}
	
	/*
	 * Receive counts for the processes to which this pe sends.
	 * This pe doesn't know who needs what it has, so we must
	 * use wildcard receives.
	 */
	
	if (l7_id_db->num_sends > l7_id_db->send_counts_len){
	   if (l7_id_db->send_counts)
	      free(l7_id_db->send_counts);
	   
	   l7_id_db->send_counts = (int *) calloc((unsigned long long)l7_id_db->num_sends, sizeof(int) );
	   if (l7_id_db->send_counts == NULL){
	      ierr = -1;
	      L7_ASSERT(l7_id_db->send_counts != NULL,
	            "Failed to allocate l7_id_db->send_counts", ierr);
	   }
	   l7_id_db->send_counts_len = l7_id_db->num_sends;
	}
	
	if (l7_id_db->num_sends > l7_id_db->send_to_len){
	   if (l7_id_db->send_to)
	      free(l7_id_db->send_to);
	   
	   l7_id_db->send_to = (int *) calloc((unsigned long long)l7_id_db->num_sends, sizeof(int) );
	   if (l7_id_db->send_to == NULL){
	      ierr = -1;
         L7_ASSERT(l7_id_db->send_to != NULL,
               "Failed to allocate l7_id_db->send_to", ierr);
	   }
	   l7_id_db->send_to_len = l7_id_db->num_sends;
	}
	
	for (i=0; i<l7_id_db->num_sends; i++){
	   ierr = MPI_Irecv(&l7_id_db->send_counts[i], 1, MPI_INT,
	         MPI_ANY_SOURCE, L7_SETUP_SEND_COUNT_TAG, MPI_COMM_WORLD,
	         &mpi_request[num_outstanding_requests++] );
	   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Irecv ( indices_needed[i] )", ierr);
	}
	
	if (num_outstanding_requests > 0){
	   ierr = MPI_Waitall(num_outstanding_requests, mpi_request, mpi_status);
	   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Waitall ( counts )", ierr);
	}
	
	num_outstanding_requests = 0;
	
	/*
	 * Determine which processes sent the above messages.
	 * These are the 'send_to' processes.
	 */
	
	offset = l7_id_db->num_recvs;
	for (i=0; i<l7_id_db->num_sends; i++){
	   l7_id_db->send_to[i] = mpi_status[offset+i].MPI_SOURCE;
	}
	
	/*
	 *  Allocate space for 'indices_global_to_send' and
	 *  'indices_local_to_send'.
	 */
	
	count_total = 0;
	for (i=0; i<l7_id_db->num_sends; i++){
	   count_total += l7_id_db->send_counts[i];
	}
	
	if (count_total > l7_id_db->indices_to_send_len){
	   if (l7_id_db->indices_global_to_send)
	      free(l7_id_db->indices_global_to_send);
	   
	   l7_id_db->indices_global_to_send = (int *) calloc((unsigned long long)count_total, sizeof(int) );
	   if (l7_id_db->indices_global_to_send == NULL){
	      ierr = -1;
	      L7_ASSERT(l7_id_db->indices_global_to_send != NULL,
	            "No memory for l7_id_db->indices_global_to_send.", ierr);
	   }
	   
	   if (l7_id_db->indices_local_to_send)
	      free(l7_id_db->indices_local_to_send);

	   l7_id_db->indices_local_to_send = (int *) calloc((unsigned long long)count_total, sizeof(int) );
      if (l7_id_db->indices_local_to_send == NULL){
         ierr = -1;
         L7_ASSERT(l7_id_db->indices_local_to_send != NULL,
               "No memory for l7_id_db->indices_local_to_send.", ierr);
      }
      
      l7_id_db->indices_to_send_len = count_total;
	}
	
	/*
	 * Send (global) indices needed from each sending process.
	 */
	
	offset = 0;
	for (i=0; i<l7_id_db->num_recvs; i++){
#if defined _L7_DEBUG
	   printf("[pe %d] Sending %d indices to pe %d. \n",
	         penum, l7_id_db->recv_counts[i], l7_id_db->recv_from[i] );
	   
	   for (k=offset; k<offset+l7_id_db->recv_counts[i]; k++){
	      printf("      index[%d] = %d \n", k, l7_id_db->indices_needed[k] );
	   }
#endif
	   
	   ierr = MPI_Isend(&l7_id_db->indices_needed[offset],
	         l7_id_db->recv_counts[i], MPI_INT,
	         l7_id_db->recv_from[i], L7_SETUP_INDICES_NEEDED_TAG,
	         MPI_COMM_WORLD, &mpi_request[num_outstanding_requests++] );
	   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Isend ( indices_needed[i] )", ierr);
	   
	   offset+=l7_id_db->recv_counts[i];
	}
	
	/*
	 * Receive (global) indices needed by the pes to which this pe sends.
	 * Note that these receives are from expected sources.
	 */
	
	offset = 0;
	for (i=0; i<l7_id_db->num_sends; i++){
	   ierr = MPI_Irecv(&l7_id_db->indices_global_to_send[offset],
	         l7_id_db->send_counts[i], MPI_INT,
	         l7_id_db->send_to[i], L7_SETUP_INDICES_NEEDED_TAG,
	         MPI_COMM_WORLD, &mpi_request[num_outstanding_requests++] );
	   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Irecv ( indices_global_to_send )", ierr);
	   
	   offset += l7_id_db->send_counts[i];
	}
	
	/*
	 * Complete indices communication.
	 */
	
	if (num_outstanding_requests > 0){
	   ierr = MPI_Waitall(num_outstanding_requests, mpi_request, mpi_status );
	   L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Waitall ( indices )", ierr);
	}
	
#if defined _L7_DEBUG
	
	ierr = MPI_Barrier(MPI_COMM_WORLD);
	offset = 0;
	
	for (j=0; j<numpes; j++){
	   if (penum == j){
	      for (i=0; i<l7_id_db->num_sends; i++){
	         printf("[pe %d] Recvd %d indices from pe %d. \n", penum,
	               l7_id_db->send_counts[i], l7_id_db->send_to[i] );
	         for (k=offset; k<offset+l7_id_db->send_counts[i]; k++){
	            printf("      index[%d] = %d \n",k l7_id_db->indices_global_to_send[k] );
	         }
	         offset += l7_id_db->send_counts[i];
	      }
	   }
	   sleep(1);
	}
#endif
	
	/* Create array of local indices corresponding to
	 *  array of global indices requested. Note the
	 * conversion from 1-based indices to 0-based is
	 * accomplished here. (See note in header).
	 */
	
	offset = 0;
	for (i=0; i<l7_id_db->num_sends; i++){
           int send_counts = l7_id_db->send_counts[i]; // for vectorization
           int adj = (int)(my_start_index) - base_adj; // for vectorization
	   for (j=0; j<send_counts; j++){
	      l7_id_db->indices_local_to_send[offset] =
	         l7_id_db->indices_global_to_send[offset] - adj;
	      offset ++;
	   }
	}
	
#if defined _L7_DEBUG
	
	ierr = MPI_Barrier(MPI_COMM_WORLD);
	
	for (i=0; i<numpes; i++){
	   if (penum == i){
	      for (j=0; j<l7_id_db->num_sends; j++){
	         printf("[pe %d] send %d indices to pe %d \n", penum,
	               l7_id_db->send_counts[j], l7_id_db->send_to[] );
	         ierr = MPI_Barrier(MPI_COMM_WORLD);
	      }
	   }
	}
	flush(stdout);
	ierr = MPI_Barrier(MPI_COMM_WORLD);
	L7_ASSERT(ierr == MPI_SUCCESS, "MPI_Barrier failure", ierr);
	
	for (i=0; i<numpes; i++){
	   if (penum == i){
	      printf("----------------------------------------------------\n")
	      for (j=0; j<l7_id_db->num_sends; j++){
	         printf("[pe %d] Send (index %d) to pe %d. \n",penum,
	               l7_id_db->indices_global_to_send[j], l7_id_db->send_to[j] );
	      }
         for (j=0; j<l7_id_db->num_recvs; j++){
            printf("[pe %d] Recving (index %d) from pe %d. \n",penum,
                  l7_id_db->indices_needed[j], l7_id_db->recv_from[j] );
         }
         printf("----------------------------------------------------\n")
         fflush(stdout);
	   }
	   sleep(2);
	}
#endif /* _L7_DEBUG */
	
	/*
	 * Ensure buffer available for data to be sent.
	 */
	
	send_buffer_bytes_needed = 0;
	
	num_sends = l7_id_db->num_sends;
	
	max_sizeof_type = sizeof(double);
	
	for (i=0; i<num_sends; i++)
	   send_buffer_bytes_needed += l7_id_db->send_counts[i] * max_sizeof_type;
	
	if (send_buffer_bytes_needed > l7.sizeof_send_buffer ){
	   if (l7.send_buffer)
	      free(l7.send_buffer);
	   
	   l7.send_buffer = (char *)calloc((unsigned long long)send_buffer_bytes_needed, sizeof (char) );
	   if (l7.send_buffer == NULL){
	      ierr = -1;
	      L7_ASSERT(l7.send_buffer != NULL, "No memory for send buffer", ierr);
	   }
	   l7.sizeof_send_buffer = send_buffer_bytes_needed;
	}

#ifdef HAVE_OPENCL
        l7_id_db->num_indices_have = 0;

        for (int i=0; i<num_sends; i++){
           /* Load data to be sent. */

           l7_id_db->num_indices_have += l7_id_db->send_counts[i];
        }

        size_t num_indices_have = l7_id_db->num_indices_have;

        l7_id_db->indices_have = (int *) malloc(num_indices_have*sizeof(int));

        int ioffset = 0;

        for (int i=0; i<num_sends; i++){
           /* Load data to be sent. */

           int send_count = l7_id_db->send_counts[i];
           for (int j=0; j<send_count; j++){
               l7_id_db->indices_have[ioffset] = l7_id_db->indices_local_to_send[ioffset];
               ioffset++;
           }
        }
        // For optimization of cache -- 
        //qsort(l7_id_db->indices_have, num_indices_have, sizeof(int), (__compar_fn_t)l7_int_comp);


        if (l7.numpes > 1) {
           l7_id_db->dev_indices_have = ezcl_malloc(NULL, "dev_indices_have", &num_indices_have,  sizeof(cl_int), CL_MEM_READ_WRITE, 0);
           cl_command_queue command_queue = ezcl_get_command_queue();
           ezcl_enqueue_write_buffer(command_queue, l7_id_db->dev_indices_have, CL_TRUE,  0, num_indices_have*sizeof(cl_int), &l7_id_db->indices_have[0],     NULL);
        }
#endif
	
	/*
	 * Message tag management
	 */
	
	l7_id_db->this_tag_update = L7_UPDATE_TAGS_MIN;
	
	/*
	 * Database is setup for this l7_id -- return.
	 */

#endif /* HAVE_MPI */
	
	ierr = L7_OK;
	
   return(ierr);
   
} /* End L7_Dev_Setup */

int l7_int_comp(const int *x, const int *y)
{
   //if (*x == *y) return 0;

   //return( (*x < *y) ? -1 : 1 )
   return *x - *y;
}

void L7_DEV_SETUP(
        const int       *my_start_index,
        const int       *num_indices_owned,
        int             *indices_needed,
        const int       *num_indices_needed,
        int             *l7_id
        )
{
   L7_Dev_Setup(0, *my_start_index, *num_indices_owned, indices_needed, *num_indices_needed, l7_id);
}
