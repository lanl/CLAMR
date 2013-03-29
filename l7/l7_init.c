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

/* This define causes this routine to have the storage for
 * for global variables here and declared extern everywhere
 * else
 */
#define L7_EXTERN
#include "l7.h"
#include "l7p.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif

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
int L7_Dev_Init(void)
{
#ifdef HAVE_OPENCL
   if (l7.numpes > 1) {
        cl_context context = ezcl_get_context();
        l7.kernel_pack_int_have_data     = ezcl_create_kernel(context, "l7/l7_kern.cl", "pack_int_have_data_cl",     0);
        l7.kernel_pack_float_have_data   = ezcl_create_kernel(context, "l7/l7_kern.cl", "pack_float_have_data_cl",   0);
        l7.kernel_pack_double_have_data  = ezcl_create_kernel(context, "l7/l7_kern.cl", "pack_double_have_data_cl",  0);
        l7.kernel_copy_ghost_int_data    = ezcl_create_kernel(context, "l7/l7_kern.cl", "copy_ghost_int_data_cl",    0);
        l7.kernel_copy_ghost_float_data  = ezcl_create_kernel(context, "l7/l7_kern.cl", "copy_ghost_float_data_cl",  0);
        l7.kernel_copy_ghost_double_data = ezcl_create_kernel(context, "l7/l7_kern.cl", "copy_ghost_double_data_cl", 0);
   }
#endif
}
