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

#define L7_LOCATION "L7_PUSH_SETUP"

int L7_Push_Setup(
		const int      num_comm_partners,
		const int      *comm_partner,
		const int      *send_buffer_count,
		int            **send_database,
		int            *receive_count_total,
		int            *l7_push_id
		)
{
	/* Purpose
	 * =======
	 * L7_Push_Setup is used to setup the update/scatter database for
	 * when senders know which data to send and to what process. Each
         * process sends in the data it needs to send in a send_database
         * which is a ragged right array with the first index the
         * num_comm_partners and the second the send_buffer_count for each
         * processor with the indices to send. From this, a communication
         * pattern and buffers are setup that allows subsequent calls to
	 * L7_Push_Update.
	 * 
	 * Arguments
	 * =========
	 * num_comm_partners    (input) const L7_INT
	 *                      global indexing set starts with 1 (Fortran)
         *                      or with 0 (C)
         *
	 * comm_partner         (input) const L7_INT
	 *                      Starting index number of calling process
	 *                      in global indexing set.
	 * 
	 * send_buffer_count    (input) const L7_INT
	 *                      Number of indices owned by calling process.
	 * 
	 * send_database        (input) const L7_INT*
	 *                      Array containing indices needed by
	 *                      calling process.
	 * 
	 * receive_count_total  (output) const L7_INT
	 *                      Number of indices of interest listed
	 *                      in array 'num_indices_needed'.
	 * 
	 * l7_push_id           (input/output) int*
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
	 * 1) The indices are handled as 4-byte integers. ??
	 * 
	 * 2) Serial compilation creates a no-op. ??
	 * 
	 * Program Flow ??
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
   
	l7_push_id_database
	  *l7_push_id_db;
	
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
/*
	if (my_start_index < 0){
		ierr = -1;
		L7_ASSERT( my_start_index >= 0, "my_start_index < 0", ierr);
	}
*/
	
	if (num_comm_partners < 0){
		ierr = -1;
		L7_ASSERT( num_comm_partners >= 0, "num_comm_partners < 0", ierr);
	}

/*
	if (num_indices_needed > 0){
		if (indices_needed == NULL){
			ierr = -1;
			L7_ASSERT( (int *)indices_needed != NULL,
					"indices_needed == NULL", ierr);
		}
	}
*/
	if (*l7_push_id < 0){
		ierr = *l7_push_id;
		L7_ASSERT( *l7_push_id >=0,
				"L7 Push Id must be either 0 (new id) or > 0 (existing id)",
				ierr);
	}

	/*
	 * Setup database structure.
	 */
	
	if (*l7_push_id != 0){
		/*
		 * Find it in the database and update based on new input.
		 */
		
		if (l7.first_push_db == NULL){
			L7_ASSERT(l7.first_push_db != NULL,
					"Uninitialized l7_push_id input, but no ids in database",
					ierr);
		}
		
		l7_push_id_db = l7.first_push_db;


		while (l7_push_id_db){
			if (l7_push_id_db->l7_push_id == *l7_push_id)
				break;
			l7_push_id_db = l7_push_id_db->next_push_db;
		}
		if (l7.first_push_db == NULL){
			ierr = -1;
			L7_ASSERT( l7.first_push_db != NULL,
					"Uninitialized l7_push_id input, but not found in this list",
					ierr);
		}
		
	}
	else{
		
		/*
		 * Allocate new database, insert into linked list.
		 */
		
		if (l7.num_push_dbs >= L7_MAX_NUM_DBS){
			ierr = -1;
			L7_ASSERT(l7.num_push_dbs < L7_MAX_NUM_DBS,
					"Too many L7 databases allocated",
					ierr);
		}

		l7_push_id_db = (l7_push_id_database*)calloc(1L, sizeof(l7_push_id_database) );
		
		if (l7_push_id_db == NULL){
			ierr = -1;
			L7_ASSERT( l7_push_id_db != NULL, "Failed to allocate new database",
					ierr);
		}
		
		if ( !(l7.first_push_db) ){
			l7.first_push_db = l7_push_id_db;
			l7.last_push_db  = l7_push_id_db;
			l7_push_id_db->next_push_db = NULL; /* Paranoia */
			
			l7_push_id_db->l7_push_id = 1;
			
			l7.num_push_dbs = 1;
		}
		else{
			
			/*
			 * Assign a l7_id.
			 */
			
			l7_push_id_db->l7_push_id = l7.last_push_db->l7_push_id + 1;
			
			/*
			 * Reset links.
			 */
			
			l7.last_push_db->next_push_db = l7_push_id_db;
			
			l7.last_push_db = l7_push_id_db;
			
			l7.num_push_dbs++;
		}
		
		*l7_push_id = l7_push_id_db->l7_push_id;
		
		/*
		 * Initialize some parameters.
		 */
		
/*
		l7_id_db->recv_counts_len = 0;
		l7_id_db->recv_from_len   = 0;
		
		l7_id_db->send_to_len     = 0;
		l7_id_db->send_counts_len = 0;
		
		l7_id_db->indices_to_send_len = 0;
		
		l7_id_db->mpi_request_len = 0;
		l7_id_db->mpi_status_len  = 0;
*/
	}
	
	/*
         *  Allocate arrays and
	 *  Store input in database.
	 */
	
        if (l7_push_id_db->num_comm_partners < num_comm_partners){

                //  comm_partner

                if (l7_push_id_db->comm_partner)
                         free(l7_push_id_db->comm_partner);

                l7_push_id_db->comm_partner = (int *) calloc(num_comm_partners,sizeof(int));

		if (l7_push_id_db->comm_partner == NULL){
			 ierr = -1;
			 L7_ASSERT( (int*)(l7_push_id_db->comm_partner) != NULL,
			            "Memory failure for comm_partner",
			            ierr);
		 }

                //  send_buffer_count

                if (l7_push_id_db->send_buffer_count)
                         free(l7_push_id_db->send_buffer_count);

                l7_push_id_db->send_buffer_count = (int *) calloc(num_comm_partners,sizeof(int));

		if (l7_push_id_db->send_buffer_count == NULL){
			 ierr = -1;
			 L7_ASSERT( (int*)(l7_push_id_db->send_buffer_count) != NULL,
			            "Memory failure for send_buffer_count",
			            ierr);
		}

                //  recv_buffer_count

                if (l7_push_id_db->recv_buffer_count)
                         free(l7_push_id_db->recv_buffer_count);

                l7_push_id_db->recv_buffer_count = (int *) calloc(num_comm_partners,sizeof(int));

		if (l7_push_id_db->recv_buffer_count == NULL){
			 ierr = -1;
			 L7_ASSERT( (int*)(l7_push_id_db->recv_buffer_count) != NULL,
			            "Memory failure for recv_buffer_count",
			            ierr);
                }

                //  send_database

                if (l7_push_id_db->send_database){
                         for (int ip = 0; ip < num_comm_partners; ip++){
                                 if (l7_push_id_db->send_database[ip]) free(l7_push_id_db->send_database[ip]);
                         }
                         if (l7_push_id_db->send_database) free(l7_push_id_db->send_database);
                }

                l7_push_id_db->send_database = (int **) calloc(num_comm_partners,sizeof(int *));
		if (l7_push_id_db->send_database == NULL){
			 ierr = -1;
			 L7_ASSERT( (int*)(l7_push_id_db->send_database) != NULL,
			            "Memory failure for send_database",
			            ierr);
                }

                for (int ip = 0; ip < num_comm_partners; ip++){
                         l7_push_id_db->send_database[ip] = (int *) calloc(send_buffer_count[ip],sizeof(int));

		         if (l7_push_id_db->send_database[ip] == NULL){
			          ierr = -1;
			          L7_ASSERT( (int*)(l7_push_id_db->send_database) != NULL,
			                     "Memory failure for send_database",
			                     ierr);
		          }
                }

                //  send_buffer

                if (l7_push_id_db->send_buffer){
                         for (int ip = 0; ip < num_comm_partners; ip++){
                                 if (l7_push_id_db->send_buffer[ip]) free(l7_push_id_db->send_buffer[ip]);
                         }
                         if (l7_push_id_db->send_buffer) free(l7_push_id_db->send_buffer);
                }

                l7_push_id_db->send_buffer = (int **) calloc(num_comm_partners,sizeof(int *));
		if (l7_push_id_db->send_buffer == NULL){
			 ierr = -1;
			 L7_ASSERT( (int*)(l7_push_id_db->send_buffer) != NULL,
			            "Memory failure for send_buffer",
			            ierr);
                }

                for (int ip = 0; ip < num_comm_partners; ip++){
                         l7_push_id_db->send_buffer[ip] = (int *) calloc(send_buffer_count[ip],sizeof(int));

		         if (l7_push_id_db->send_buffer[ip] == NULL){
			          ierr = -1;
			          L7_ASSERT( (int*)(l7_push_id_db->send_buffer) != NULL,
			                     "Memory failure for send_buffer",
			                     ierr);
		          }
                }

        }

        /*
         *  Copy input data into database
         */

        l7_push_id_db->num_comm_partners = num_comm_partners;

        for (int ip = 0; ip < num_comm_partners; ip++){
                l7_push_id_db->comm_partner[ip] = comm_partner[ip];
                l7_push_id_db->send_buffer_count[ip] = send_buffer_count[ip];
        }

        for (int ip = 0; ip < num_comm_partners; ip++){
                int count = send_buffer_count[ip]; // create simple int count to help vectorization
                for (int ic = 0; ic < count; ic++){
                        l7_push_id_db->send_database[ip][ic] = send_database[ip][ic];
                }
        }

        /*
         * Get receive counts by communication
         */

        MPI_Request request[2*num_comm_partners];
        MPI_Status  status[2*num_comm_partners];
         
        for (int ip = 0; ip < num_comm_partners; ip++){
                MPI_Irecv(&l7_push_id_db->recv_buffer_count[ip], 1, MPI_INT, l7_push_id_db->comm_partner[ip],
                          l7_push_id_db->comm_partner[ip], MPI_COMM_WORLD, &request[ip]);
        }

        for (int ip = 0; ip < num_comm_partners; ip++){
                MPI_Isend(&l7_push_id_db->send_buffer_count[ip], 1, MPI_INT, l7_push_id_db->comm_partner[ip],
                          l7.penum, MPI_COMM_WORLD, &request[num_comm_partners+ip]);
        }
        MPI_Waitall(2*num_comm_partners, request, status);

        /*
         * Calculate sum of receives
         */
        *receive_count_total = 0;
        for (int ip = 0; ip < num_comm_partners; ip++){
                *receive_count_total += l7_push_id_db->recv_buffer_count[ip];
        }
        l7_push_id_db->receive_count_total = *receive_count_total;

#endif /* HAVE_MPI */
	
	ierr = L7_OK;
	
        return(ierr);
   
} /* End L7_Push_Setup */

