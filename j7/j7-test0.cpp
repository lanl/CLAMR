/*
 *  Copyright (c) 2013, Los Alamos National Security, LLC.
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
 *
 *  CLAMR -- LA-CC-11-094
 *
 *  Author: Samuel K. Gutierrez
 */

#include "j7.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstdlib>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "mpi.h"

using namespace std;

static void
demoEmitSync(MPI_Comm comm,
             int rank)
{
    MPI_Barrier(comm);
    usleep(rank * 9000);
}

template <typename T>
static void
emitElems(MPI_Comm comm, int commID, const T &t, int nElems)
{
    demoEmitSync(comm, commID);
    cout << "rank: " << commID << ": ";
    for (int i = 0; i < nElems; ++i) {
        cout << t[i] << " ";
    }
    cout << endl;
}

static size_t
getNElems(int commID, int baseNElems)
{
    return max(1, baseNElems - commID);
}

int
main(int argc, char **argv)
{
    size_t segSize = 1024 * 64;
    int smCommRank = 0, smCommSize = 0;
    MPI_Comm smComm;
    J7 *j7 = NULL;

    if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
        cerr << "MPI_Init failure!" << endl;
        exit(EXIT_FAILURE);
    }
    if (MPI_SUCCESS != MPI_Comm_dup(MPI_COMM_WORLD, &smComm)) {
        cerr << "MPI_Comm_dup failure!" << endl;
        exit(EXIT_FAILURE);
    }
    if (MPI_SUCCESS != MPI_Comm_size(smComm, &smCommSize)) {
        cerr << "MPI_Comm_size failure!" << endl;
        exit(EXIT_FAILURE);
    }
    if (MPI_SUCCESS != MPI_Comm_rank(smComm, &smCommRank)) {
        cerr << "MPI_Comm_rank failure!" << endl;
        exit(EXIT_FAILURE);
    }
    try {
        char *cArray = NULL;
        int nElems = getNElems(smCommRank, 5), *iArray = NULL;
        double *dArray = NULL, *spanArray = NULL;

        j7 = new J7(smComm, segSize);
        // allocation tests
        cArray = static_cast<char *>(j7->memAlloc(sizeof(char) * nElems));
        iArray = static_cast<int *>(j7->memAlloc(sizeof(int) * nElems));
        dArray = static_cast<double *>(j7->memCalloc(nElems, sizeof(double)));
        spanArray = static_cast<double *>(j7->memCalloc(nElems, sizeof(double)));

        for (int i = 0; i < nElems; ++i) {
            cArray[i] = ('A' +  i + smCommRank) % 'z';
            iArray[i] = i + smCommRank;
            dArray[i] = (i + smCommRank) * 0.5;
        }
        // remember that is is a collective operation within the give comm
        double *gsa = static_cast<double *>(j7->getBaseMemPtr(spanArray));
        if (smCommRank == (smCommSize - 1)) {
            for (int i = 0; i < nElems * smCommSize; ++i) {
                gsa[i] = 42;
            }
        }
        emitElems(smComm, smCommRank, dArray, nElems);
        emitElems(smComm, smCommRank, cArray, nElems);
        emitElems(smComm, smCommRank, iArray, nElems);
        emitElems(smComm, smCommRank, spanArray, nElems);
        // ordering test
        j7->memFree(iArray);
        // reallocation tests
        cArray = static_cast<char *>
                 (j7->memRealloc(cArray, sizeof(char) * nElems * 2));
        // now everyone add stuff to the end of the new space
        for (int i = nElems; i < nElems * 2; ++i) {
            cArray[i] = ('A' +  i + smCommRank) % 'z';
        }
        emitElems(smComm, smCommRank, cArray, nElems * 2);
        // release allocated resources
        MPI_Barrier(smComm);
        j7->memFree(dArray);
        j7->memFree(cArray);
        j7->memFree(spanArray);
        delete j7;
    }
    catch (J7Exception &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    if (MPI_SUCCESS != MPI_Comm_free(&smComm)) {
        cerr << "MPI_Comm_free failure!" << endl;
        exit(EXIT_FAILURE);
    }
    if (MPI_SUCCESS != MPI_Finalize()) {
        cerr << "MPI_Init failure!" << endl;
        exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
