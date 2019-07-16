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
 */
#include <signal.h>
#include <time.h>
#include "l7.h"
#include "l7p.h"

int L7_Get_Rank(void)
{
#ifdef HAVE_MPI
    return(l7.penum);
#else
    return(0);
#endif
}

int l7_get_rank_(void)
{
#ifdef HAVE_MPI
    return(l7.penum);
#else
    return(0);
#endif
}

int L7_Get_Numpes(void)
{
#ifdef HAVE_MPI
   return(l7.numpes);
#else
   return(1);
#endif
}

int l7_get_numpes_(void)
{
#ifdef HAVE_MPI
   return(l7.numpes);
#else
   return(1);
#endif
}

double L7_Wtime(void)
{
    double l7_wtime = 0.0;
#ifdef HAVE_MPI
    if (l7.initialized_mpi){
        l7_wtime = MPI_Wtime();
    }
    else {
        l7_wtime = (double)(clock()/CLOCKS_PER_SEC);
    }
#else
   l7_wtime = clock()/CLOCKS_PER_SEC;
#endif
    return(l7_wtime);
}

double l7_wtime_(void)
{
    double l7_wtime = 0.0;
#ifdef HAVE_MPI
    if (l7.initialized_mpi){
        l7_wtime = MPI_Wtime();
    }
    else {
        l7_wtime = (double)(clock()/CLOCKS_PER_SEC);
    }
#else
   l7_wtime = clock()/CLOCKS_PER_SEC;
#endif
    return(l7_wtime);
}

int L7_Get_Timeout_Signal(void)
{
    return(SIGURG);
}

int l7_get_timeout_signal_(void)
{
    return(SIGURG);
}

long long L7_Address(void *var)
{
   return((long long)var);
}
long long l7_address_(void *var)
{
   return((long long)var);
}

