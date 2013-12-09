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

#define L7_LOCATION "L7_FREE"

int L7_Dev_Free(
      void
      )
{
   /*
    * Purpose
    * =======
    * L7_Dev_Free frees the kernels compiled in L7_Dev_Init
    * 
    * Arguments
    * =========
    * none
    * 
    * Return value
    * ============
    * Value other than L7_OK indicates an error.
    * 
    * Notes:
    * =====
    * 
    */
   
   /*
    * Local variables.
    */
   
   int
     ierr;            /* Error code for return              */
   
#if defined HAVE_MPI
   
   /*
    * Executable Statements
    */
   
   if (! l7.mpi_initialized){
      return(0);
   }
   
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &l7.penum);
   L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank error", ierr);
 
   if (l7.initialized != 1){
      ierr = -1;
      L7_ASSERT( l7.initialized != 1, "L7 not initialized", ierr);
   }
   
   /*
    *  Get first L7 structure and null out previous structure
    */
   
   /*
    * Free kernels
    */

#ifdef HAVE_OPENCL
        if (l7.numpes > 1){
           ezcl_kernel_release(l7.kernel_pack_int_have_data);
           ezcl_kernel_release(l7.kernel_pack_float_have_data);
           ezcl_kernel_release(l7.kernel_pack_double_have_data);
           ezcl_kernel_release(l7.kernel_copy_ghost_int_data);
           ezcl_kernel_release(l7.kernel_copy_ghost_float_data);
           ezcl_kernel_release(l7.kernel_copy_ghost_double_data);
        }
#endif

#endif
	
   ierr = L7_OK;
   
   return(ierr);
   
}

void L7_DEV_FREE(
        int         *ierr
        )
{
   L7_Dev_Free();
   *ierr = 0;
}


