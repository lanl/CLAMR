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
#include <stdlib.h>
#include <dlfcn.h>

#define L7_LOCATION "L7_INIT"

#ifdef HAVE_OPENCL
#include "l7_kernel.inc"
#endif

#ifdef HAVE_QUO
static int getSubCommProcs(QUO_context c, int numpes, int *vLen, int **v);
static int subCommInit(QUO_SubComm *subcomm, int np1s, int *p1who);
#endif
#define WARNING_SUPPRESSION 0

int L7_Init (
    int *mype,
    int *numpes,
    int *argc,
    char **argv,
    int do_quo_setup,
    int lttrace_on
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

#ifndef HAVE_QUO
   // To get rid of compiler warning
   if (WARNING_SUPPRESSION && do_quo_setup == 99) printf("DEBUG do_quo_setup = %d\n",do_quo_setup);
#endif
#ifndef HAVE_LTTRACE
   // To get rid of compiler warning
   if (WARNING_SUPPRESSION && lttrace_on == 99) printf("DEBUG lttrace_on = %d\n",lttrace_on);
#endif

#ifdef HAVE_LTTRACE
   if (lttrace_on){
      void *handle = dlopen("lttrace.so",RTLD_LAZY);
      if (! handle) {
         printf("DEBUG -- open failed\n");
      }
      int (*init_lttrace)(int *argc, char ***argv);
      init_lttrace = (int (*)(int *, char ***))dlsym(handle, "initialize_lttrace");
   }  
#endif
      
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
        ierr = MPI_Comm_rank (MPI_COMM_WORLD, &l7.penum );
        L7_ASSERT( ierr == MPI_SUCCESS, "MPI_Comm_rank", ierr );

        ierr = MPI_Comm_size (MPI_COMM_WORLD, &l7.numpes );
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

#ifdef HAVE_QUO
   if (do_quo_setup) {
      // init QUO -- all MPI processes MUST do this at the same time.
      // create shorthand for context and Subcomm
      QUO_context context = l7.subComm.context;
      QUO_SubComm subcomm = l7.subComm;

      subcomm.rank = -1;
      subcomm.size = -1;

      int rc = QUO_ERR;
      if (QUO_SUCCESS != (rc = QUO_create(&context))){
         printf("QUO_create failure: rc = %d\n",rc);
      } 

      if (l7.penum == 0) {
         int nnodes, nnumanodes, nsockets, ncores, npus, nqids;
         if (QUO_SUCCESS != QUO_nnodes(context, &nnodes)){
            printf("QUO_nnodes failure\n");
         }
         if (QUO_SUCCESS != QUO_nnumanodes(context, &nnumanodes)){
            printf("QUO_nnumanodes failure\n");
         }
         if (QUO_SUCCESS != QUO_nsockets(context, &nsockets)){
            printf("QUO_nsockets failure\n");
         }
         if (QUO_SUCCESS != QUO_ncores(context, &ncores)){
            printf("QUO_ncores failure\n");
         }
         if (QUO_SUCCESS != QUO_npus(context, &npus)){
            printf("QUO_npus failure\n");
         }
         if (QUO_SUCCESS != QUO_nqids(context, &nqids)){
            printf("QUO_nqids failure\n");
         }

         printf("### System Info ###\n");
         printf("# Nodes: %d\n",nnodes);
         printf("# NUMA Nodes: %d\n",nnumanodes);
         printf("# Sockets: %d\n",nsockets);
         printf("# Cores: %d\n",ncores);
         printf("# PUs: %d\n",npus);
         printf("# MPI Procs on Node: %d\n",nqids);
      }

      int numPEsInSubComm = 0;
      int *cwPEs = NULL;
      rc = QUO_ERR;

      getSubCommProcs(context, l7.numpes, &numPEsInSubComm, &cwPEs);

      int n;
      if (QUO_SUCCESS != (rc = QUO_id(context, &n))) {
         printf("QUO_id: rc = %d\n",rc);
         exit(1);
      }
      if (0 == n) {
         // add hostname -- could be diff on diff systems
         printf("MPI_COMM_WORLD ranks in subComm: ");
         for (int i = 0; i < numPEsInSubComm; ++i) {
             printf(" %d ",cwPEs[i]);
         }
         printf("\n");
      }

      // now actually create the darn thing...
      rc = subCommInit(&subcomm, numPEsInSubComm, cwPEs);
      if (cwPEs) free(cwPEs);

      // at this point subComm is ready to use for those inSubComm

      int member_global[l7.numpes];
      int subcomm_rank_global[l7.numpes];
      MPI_Gather(&subcomm.member, 1, MPI_INT,
                 member_global, 1, MPI_INT,
                 0, MPI_COMM_WORLD);
      MPI_Gather(&subcomm.rank, 1, MPI_INT,
                 subcomm_rank_global, 1, MPI_INT,
                 0, MPI_COMM_WORLD);

      if (l7.penum == 0) {
         for (int ip = 0; ip < l7.numpes; ip++){
            if (member_global[ip]) {
               printf("rank: %d is subComm rank: %d\n",ip, subcomm_rank_global[ip]);
            }
         }

        int world_size, subcomm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_size(subcomm.comm, &subcomm_size);
        printf("Size comm world %d subComm %d\n",world_size, subcomm_size);
      }

   }

#endif


#else

   *mype = 0;
   *numpes = 1;

#endif /* HAVE_MPI */

   return(ierr);     
}
void L7_Dev_Init(void)
{
#ifdef HAVE_OPENCL
   if (l7.numpes > 1) {
        cl_context context = ezcl_get_context();

        l7.kernel_pack_int_have_data     = ezcl_create_kernel_wsource(context, l7_kern_source, "pack_int_have_data_cl");
        l7.kernel_pack_float_have_data   = ezcl_create_kernel_wsource(context, l7_kern_source, "pack_float_have_data_cl");
        l7.kernel_pack_double_have_data  = ezcl_create_kernel_wsource(context, l7_kern_source, "pack_double_have_data_cl");
        l7.kernel_copy_ghost_int_data    = ezcl_create_kernel_wsource(context, l7_kern_source, "copy_ghost_int_data_cl");
        l7.kernel_copy_ghost_float_data  = ezcl_create_kernel_wsource(context, l7_kern_source, "copy_ghost_float_data_cl");
        l7.kernel_copy_ghost_double_data = ezcl_create_kernel_wsource(context, l7_kern_source, "copy_ghost_double_data_cl");
   }
#endif
}

#ifdef HAVE_QUO
static int getSubCommProcs(QUO_context c, int numpes, int *vLen, int **v)
{
    QUO_obj_type_t resPrio[] = {QUO_OBJ_NODE,
                                QUO_OBJ_SOCKET};
    // default target resource is a NUMA node
    QUO_obj_type_t targetRes = QUO_OBJ_NODE;
    // number of target resources and max number of procs per resource
    int nRes = 0, maxProcPerRes = 1;
    int res_assigned = 0;
    int totalWorkers = 0;
    int rc = MPI_SUCCESS;
    // array that hold whether or not a particular rank is going to do work
    int *workContribs = NULL;
    // MPI_COMM_WORLD ranks of the selected workers
    int *workerRanks = NULL;

    // figure out what we are going to distribute work over
    for (uint i = 0; i < sizeof(resPrio) / sizeof(resPrio[0]); ++i) {
        if (QUO_SUCCESS != (rc = QUO_nobjs_by_type(c, resPrio[i], &nRes))){
           printf("QUO_nobjs_by_type: rc = %d\n",rc);
           exit(1);
        }
        if (nRes > 0){
            targetRes = resPrio[i];
            break;
        }
    }
    // failure -- fix this path at some point
    if (0 == nRes) return 1;
    /* let quo distribute workers over the sockets. if res_assigned is 1
     * after this call, then i have been chosen. */
    int isel = 0;
    if (QUO_SUCCESS != (rc = QUO_auto_distrib(c, targetRes,
                                              maxProcPerRes, &isel))) {
        printf("QUO_auto_distrib: rc = %d\n",rc);
        exit(1);
    }
    if (isel == 1) {
    //if (*c.quo->autoDistrib(targetRes, maxProcPerRes)) {
        res_assigned = 1;
    }

   /* array that hold whether or not a particular rank is going to do work */
   workContribs = (int *)calloc(numpes, sizeof(*workContribs));
   if (!workContribs) return 1;

   if (MPI_SUCCESS != (rc = MPI_Allgather(&res_assigned, 1, MPI_INT,
                                          workContribs, 1, MPI_INT,
                                          MPI_COMM_WORLD))) {
       return 1;
   }
   /* now iterate over the array and count the total number of workers */
   for (int i = 0; i < numpes; ++i) {
       if (1 == workContribs[i]) ++totalWorkers;
   }
   workerRanks = (int *)calloc(totalWorkers, sizeof(*workerRanks));
   if (!workerRanks) return 1;
   /* populate the array with the worker comm world ranks */
   for (int i = 0, j = 0; i < numpes; ++i) {
       if (1 == workContribs[i]) {
           workerRanks[j++] = i;
       }
   }
   *vLen = totalWorkers;
   *v = workerRanks;
   if (workContribs) free(workContribs);
   return 0;
}

static int
subCommInit(QUO_SubComm *subcomm,
            int np1s /* number of participants |p1who| */,
            int *p1who /* the participating ranks (MPI_COMM_WORLD) */)
{
    int rc = QUO_SUCCESS;
    MPI_Group worldGroup;
    MPI_Group p1_group;
    /* ////////////////////////////////////////////////////////////////////// */
    /* now create our own communicator based on the rank ids passed here */
    /* ////////////////////////////////////////////////////////////////////// */
    if (MPI_SUCCESS != MPI_Comm_group(MPI_COMM_WORLD, &worldGroup)) {
        rc = QUO_ERR_MPI;
        goto out;
    }
    if (MPI_SUCCESS != MPI_Group_incl(worldGroup, np1s,
                                      p1who, &p1_group)) {
        rc = QUO_ERR_MPI;
        goto out;
    }
    if (MPI_SUCCESS != MPI_Comm_create(MPI_COMM_WORLD,
                                       p1_group,
                                       &((*subcomm).comm))) {
        rc = QUO_ERR_MPI;
        goto out;
    }
    /* am i in the new communicator? */
    (*subcomm).member = (MPI_COMM_NULL == (*subcomm).comm) ? 0 : 1;
    if ((*subcomm).member) {
        if (MPI_SUCCESS != MPI_Comm_size((*subcomm).comm, &(*subcomm).size)) {
            rc = QUO_ERR_MPI;
            goto out;
        }
        if (MPI_SUCCESS != MPI_Comm_rank((*subcomm).comm, &(*subcomm).rank)) {
            rc = QUO_ERR_MPI;
            goto out;
        }
    }
out:
    if (MPI_SUCCESS != MPI_Group_free(&worldGroup)) return 1;
    return (QUO_SUCCESS == rc) ? 0 : 1;
}

#endif
