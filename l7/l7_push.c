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
#include <stdio.h>
#include <mpi.h>

// This implementation is very rough and needs to be upgraded to the level of the 
// original update (pull) routines. Right now it only supports one handle at a 
// time.

struct push_comm_type {
   int num_comm_partners;
   int *comm_partner;
   int *send_buffer_count;
   int *recv_buffer_count;
   int **send_database;
   int receive_count_total;
} push_comm;

int L7_Push_Setup(int num_comm_partners, int *comm_partner, int *send_buffer_count, int **send_database, int *receive_count_total)
{
   push_comm.num_comm_partners = num_comm_partners;
   push_comm.comm_partner      = (int *) malloc(num_comm_partners*sizeof(int));
   push_comm.send_buffer_count = (int *) malloc(num_comm_partners*sizeof(int));
   push_comm.recv_buffer_count = (int *) malloc(num_comm_partners*sizeof(int));

   int mype;
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);

   for (int ip = 0; ip < num_comm_partners; ip++){
      push_comm.comm_partner[ip] = comm_partner[ip];
      push_comm.send_buffer_count[ip] = send_buffer_count[ip];
   }

   push_comm.send_database = (int **) malloc(num_comm_partners*sizeof(int *));
   for (int ip = 0; ip < num_comm_partners; ip++){
      push_comm.send_database[ip] = (int *) malloc(send_buffer_count[ip]*sizeof(int));
   }

   for (int ip = 0; ip < num_comm_partners; ip++){
      for (int ic = 0; ic < send_buffer_count[ip]; ic++){
         push_comm.send_database[ip][ic] = send_database[ip][ic];
      }
   }

#define MPI_ASYNC
#ifdef MPI_ASYNC
   MPI_Request request[2*num_comm_partners];
   MPI_Status  status[2*num_comm_partners];
         
   for (int ip = 0; ip < num_comm_partners; ip++){
      MPI_Irecv(&push_comm.recv_buffer_count[ip], 1, MPI_INT, comm_partner[ip], comm_partner[ip], MPI_COMM_WORLD, &request[ip]);
   }

   for (int ip = 0; ip < num_comm_partners; ip++){
      MPI_Isend(&send_buffer_count[ip], 1, MPI_INT, comm_partner[ip], mype, MPI_COMM_WORLD, &request[num_comm_partners+ip]);
   }
   MPI_Waitall(2*num_comm_partners, request, status);
#else
   MPI_Status status;
   for (int ip = 0; ip < num_comm_partners; ip++){
      MPI_Sendrecv(&send_buffer_count[ip], 1, MPI_INT, comm_partner[ip], mype,
                   &recv_buffer_count[ip], 1, MPI_INT, comm_partner[ip], comm_partner[ip], MPI_COMM_WORLD, &status);
   }
#endif

   *receive_count_total = 0;
   for (int ip = 0; ip < num_comm_partners; ip++){
      *receive_count_total += push_comm.recv_buffer_count[ip];
   }
   push_comm.receive_count_total = *receive_count_total;

   return(1);
}

void L7_Push_Update(int ihandle, int *array, int *return_array)
{
   int mype;
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);

   int **send_buffer = (int **) malloc(push_comm.num_comm_partners*sizeof(int *)); 
   for (int ip = 0; ip < push_comm.num_comm_partners; ip++){
      send_buffer[ip] = (int *) malloc(push_comm.send_buffer_count[ip]*sizeof(int));
   }    

   for (int ip = 0; ip < push_comm.num_comm_partners; ip++){
      for (int ic = 0; ic < push_comm.send_buffer_count[ip]; ic++){
         int ib = push_comm.send_database[ip][ic];
         send_buffer[ip][ic] = array[ib];
      }    
   }    


// Send/Receives will be done in L7_Push_Update. Input will be send_buffer. Output will be in
// preallocated receive_buffer
   MPI_Request request[2*push_comm.num_comm_partners];
   MPI_Status  status[2*push_comm.num_comm_partners];

   int iloc = 0;
   for (int ip = 0; ip < push_comm.num_comm_partners; ip++){
      MPI_Irecv(&return_array[iloc], push_comm.recv_buffer_count[ip], MPI_INT, push_comm.comm_partner[ip], push_comm.comm_partner[ip], MPI_COMM_WORLD, &request[ip]);
      iloc += push_comm.recv_buffer_count[ip];
   }

   for (int ip = 0; ip < push_comm.num_comm_partners; ip++){
      MPI_Isend(send_buffer[ip], push_comm.send_buffer_count[ip], MPI_INT, push_comm.comm_partner[ip], mype, MPI_COMM_WORLD, &request[push_comm.num_comm_partners+ip]);
   }    
   MPI_Waitall(2*push_comm.num_comm_partners, request, status);

/*
   if (ncycle >= 1) { 
      for (int ib = 0; ib<receive_count_total; ib++){
         fprintf(fp,"DEBUG receive %d is %d\n",ib,border_data_receive[ib]);
      }    
   }    
*/

}

void L7_Push_Free(int i_handle)
{
   free(push_comm.comm_partner);
   free(push_comm.send_buffer_count);
   free(push_comm.recv_buffer_count);
   push_comm.comm_partner = NULL;
   push_comm.send_buffer_count = NULL;
   push_comm.recv_buffer_count = NULL;

   for (int ip = 0; ip < push_comm.num_comm_partners; ip++){
      free(push_comm.send_database[ip]);
   }
   free(push_comm.send_database);
   push_comm.send_database = NULL;

   push_comm.num_comm_partners = -1;
}
