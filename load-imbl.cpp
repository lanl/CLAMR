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

#include "j7/j7.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <vector>

#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <getopt.h>
#include <math.h>

#include <omp.h>
#include "QUO.hpp"
#include "mpi.h"

using namespace std;

// allocation types
enum AllocType {
    ALLOC_SYS,
    ALLOC_J7
};

// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //
// number of iterations
static size_t nIters = 128;
// base number of elements
static const size_t ALLOC_BASE_NELEMS = 256 * 1024;
// allocator type
AllocType allocType = ALLOC_SYS;
// FIXME could tune this
static int ompChunk = 1024;
// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

// comm world rank
static int cwID = 0;
// sm comm rank
static int smID = 0;
// size of comm world
static int totPEs = 0;
// size of sm comm
static int totSMPEs = 0;
// sm comm
static MPI_Comm smComm;
// flag indicating whether or not i'm doing work.
bool workerProc = false;

static J7 *j7 = NULL;
static QUO *quo = NULL;

// number of work arrays (i.e. # of distinct work arrays (a & b))
static const size_t N_WORK_ARRAYS = 2;
// work arrays -- ugly way of structuring the code, but can make things
// clearer... that is, at least for me :-)
// cells, for example
static double *a1 = NULL, *a2 = NULL, *a3 = NULL; // id 0
// materials, for example
static double *b1 = NULL, *b2 = NULL, *b3 = NULL; // id 1
// array of work item sizes -- may be different across processes
static size_t nElems[N_WORK_ARRAYS];
// number of total elements on a compute node resource
static size_t nNodeGlobalElems[N_WORK_ARRAYS];

// times to work on the work arrays
static vector<double> aTimes(nIters);
static vector<double> bTimes(nIters);

/* host name buffer */
static char hostNameBuff[MPI_MAX_PROCESSOR_NAME];

/* ////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////// */
/* smp comm stuff */
/* ////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////// */
static unsigned long int
getNetNum(char *hstn)
{
    struct hostent *host = NULL;

    if (NULL == (host = gethostbyname(hstn))) {
        throw J7Exception(__FILE__, __LINE__, "*** gethostbyname ***");
    }
    /* htonl used here because nodes could be different architectures */
    return htonl(inet_network(inet_ntoa(*(struct in_addr *)host->h_addr)));
}

static int
cmpULI(const void *p1,
       const void *p2)
{
    return (*(unsigned long int *)p1 - *(unsigned long int *)p2);
}

static int
getColor(unsigned long int *net_nums,
         int net_num_len,
         unsigned long int my_net_num )
{
    int i      = 0;
    int node_i = 0;
    unsigned long int prev_num;

    qsort(net_nums, (size_t)net_num_len, sizeof(unsigned long int), cmpULI);

    prev_num = net_nums[0];

    while (i < net_num_len && prev_num != my_net_num) {
        while (net_nums[i] == prev_num) {
            ++i;
        }
        ++node_i;
        prev_num = net_nums[i];
    }
    return node_i;
}


static void
emitSync(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(cwID * 9000);
}

template <typename T>
static void
emitElems(const T &t, int nElems, string name)
{
    // print max 16 elems
    int limit = min(nElems, 16);
    emitSync();
    cout << "rank: " << cwID << " [" << name << "]: ";
    for (int i = 0; i < limit; ++i) {
        cout << t[i] << " ";
    }
    cout << endl;
}

static void
createSMComm(void)
{
    string badFunc;
    int len = 0, color = 0;
    unsigned long int netNum = 0;
    unsigned long int *netNums = 0;

    if (MPI_SUCCESS != MPI_Get_processor_name(hostNameBuff, &len)) {
        badFunc = string("MPI_Get_processor_name");
        goto out;
    }
    netNum = getNetNum(hostNameBuff);
    // alloc enough to store all netNums
    netNums = (unsigned long int *)malloc(sizeof(*netNums) * totPEs);
    if (!netNums) {
        badFunc = string("*** bad alloc ***");
        goto out;
    }
    /* get everyone else's net_num value */
    if (MPI_SUCCESS != MPI_Allgather(&netNum, 1, MPI_UNSIGNED_LONG, netNums, 1,
                                     MPI_UNSIGNED_LONG, MPI_COMM_WORLD)) {
        badFunc = string("MPI_Get_processor_name");
        goto out;
    }
    color = getColor(netNums, totPEs, netNum);
    /* free up some resources - no longer needed */
    /* split into local node groups */
    if (MPI_SUCCESS != MPI_Comm_split(MPI_COMM_WORLD, color, cwID, &smComm)) {
        badFunc = string("MPI_Comm_split");
        goto out;
    }
    /* get sm comm size */
    if (MPI_SUCCESS != MPI_Comm_size(smComm, &totSMPEs)) {
        badFunc = string("MPI_Comm_size");
        goto out;
    }
    /* get my sm comm rank */
    if (MPI_SUCCESS != MPI_Comm_rank(smComm, &smID)) {
        badFunc = string("MPI_Comm_rank");
        goto out;
    }
    free(netNums);
out:
    if (!badFunc.empty()) {
        throw J7Exception(__FILE__, __LINE__, string(badFunc));
    }
}

static void
initMPI(int argc, char **argv)
{
    string badFunc;

    if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
        badFunc = string("MPI_Init failure!");
        goto out;
    }
    if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &totPEs)) {
        badFunc = string("MPI_Comm_size failure!");
        goto out;
    }
    if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &cwID)) {
        badFunc = string("MPI_Comm_rank failure!");
        goto out;
    }
    // create smp communicator (all the ranks that share a compute node)
    createSMComm();
out:
    if (!badFunc.empty()) {
        throw J7Exception(__FILE__, __LINE__, string(badFunc));
    }
}

static void
finiMPI(void)
{
    string badFunc;

    if (MPI_SUCCESS != MPI_Comm_free(&smComm)) {
        badFunc = string("MPI_Comm_free failure!");
        goto out;
    }
    if (MPI_SUCCESS != MPI_Finalize()) {
        badFunc = string("MPI_Finalize failure!");
        goto out;
    }
out:
    if (!badFunc.empty()) {
        throw J7Exception(__FILE__, __LINE__, string(badFunc));
    }
}

static size_t
getNAllocdOnNode(void)
{
    size_t nipd = 0;
    for (size_t i = 0; i < N_WORK_ARRAYS; ++i) {
        nipd += nNodeGlobalElems[i];
    }
    return nipd * 3UL;
}

static void
initJ7(void)
{
    try {
        size_t segSize = getNAllocdOnNode() * sizeof(*a1);
        // add fluff
        segSize += 1024;
        j7 = new J7(smComm, segSize);
    }
    catch (J7Exception &e) {
        throw e;
    }
}

static void
initQUO(void)
{
    try {
        quo = new QUO();
        quo->create();
    }
    catch (J7Exception &e) {
        throw e;
    }
}

static void
finiJ7(void)
{
    try {
        delete j7;
    }
    catch (J7Exception &e) {
        throw e;
    }
}

static void
finiQUO(void)
{
    try {
        quo->free();
        delete quo;
    }
    catch (J7Exception &e) {
        throw e;
    }
}

#if 0
// this sets alloc sizes for everyone
static size_t
getNumElems(size_t arrayID)
{
    // cells
    if (0 == arrayID) {
        if (0 == cwID) return 4;
        return 1;
    }
    // particles
    else if (1 == arrayID) {
        if (0 == cwID) return 3;
        return 1;
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** bad array ID ***");
    }
    return 0;
}
#endif

#if 0 // exp00
static size_t
getNumElems(size_t arrayID) // Partitioned by Rows Based on Cells
{
    // cells
    if (0 == arrayID) {
        return 4 * ALLOC_BASE_NELEMS;
    }
    // particles
    else if (1 == arrayID) {
        if      (0 == cwID) return 24 * ALLOC_BASE_NELEMS;
        else if (1 == cwID) return 14 * ALLOC_BASE_NELEMS;
        else if (2 == cwID) return 9  * ALLOC_BASE_NELEMS;
        else if (3 == cwID) return 6  * ALLOC_BASE_NELEMS;
        else throw J7Exception(__FILE__, __LINE__, "*** bad setup ***");
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** bad array ID ***");
    }
    return 0;
}
#endif

#if 0 // exp01
static size_t
getNumElems(size_t arrayID) // Partitioned by Quadrants Based on Cells
{
    // cells
    if (0 == arrayID) {
        return 4 * ALLOC_BASE_NELEMS;
    }
    // particles
    else if (1 == arrayID) {
        if      (0 == cwID) return 29 * ALLOC_BASE_NELEMS;
        else if (1 == cwID) return 9 * ALLOC_BASE_NELEMS;
        else if (2 == cwID) return 9  * ALLOC_BASE_NELEMS;
        else if (3 == cwID) return 6  * ALLOC_BASE_NELEMS;
        else throw J7Exception(__FILE__, __LINE__, "*** bad setup ***");
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** bad array ID ***");
    }
    return 0;
}
#endif

#if 1 // exp02
static size_t
getNumElems(size_t arrayID) // Partitioned by Rows Based
{                           // Both on Weighted Cells and Particles
    // cells
    if (0 == arrayID) {
        if      (0 == cwID) return 1 *ALLOC_BASE_NELEMS;
        else if (1 == cwID) return 4 * ALLOC_BASE_NELEMS;
        else if (2 == cwID) return 5 * ALLOC_BASE_NELEMS;
        else if (3 == cwID) return 6 * ALLOC_BASE_NELEMS;
        else throw J7Exception(__FILE__, __LINE__, "*** bad setup ***");
    }
    // particles
    else if (1 == arrayID) {
        if      (0 == cwID) return 14 * ALLOC_BASE_NELEMS;
        else if (1 == cwID) return 16 * ALLOC_BASE_NELEMS;
        else if (2 == cwID) return 13 * ALLOC_BASE_NELEMS;
        else if (3 == cwID) return 10 * ALLOC_BASE_NELEMS;
        else throw J7Exception(__FILE__, __LINE__, "*** bad setup ***");
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** bad array ID ***");
    }
    return 0;
}
#endif

static void
setNElems(void)
{
    for (size_t arrayID = 0; arrayID < N_WORK_ARRAYS; ++arrayID) {
        nElems[arrayID] = getNumElems(arrayID);
    }
    // now set the (node) global values for j7 that way we know the "real"
    // extent of the arrays in shared memory
    for (size_t arrayID = 0; arrayID < N_WORK_ARRAYS; ++arrayID) {
        unsigned long long v = static_cast<unsigned long long>(nElems[arrayID]);
        unsigned long long nodeGlobalVali = 0;
        if (MPI_SUCCESS != MPI_Allreduce(&v, &nodeGlobalVali, 1,
                                         MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                                         smComm)) {
            throw J7Exception(__FILE__, __LINE__, "*** MPI_Allreduce ***");
        }
        nNodeGlobalElems[arrayID] = static_cast<size_t>(nodeGlobalVali);
    }
}

static void *
memCalloc(AllocType allocType, size_t nElems, size_t elemSize)
{
    if (ALLOC_SYS == allocType) {
        return calloc(nElems, elemSize);
    }
    else if (ALLOC_J7 == allocType) {
        return j7->memCalloc(nElems, elemSize);
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** alloc type ***");
    }
    return NULL;
}

static void
memFree(AllocType allocType, void *p)
{
    if (ALLOC_SYS == allocType) {
        free(p);
    }
    else if (ALLOC_J7 == allocType) {
        j7->memFree(p);
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** alloc type ***");
    }
}

static void
allocWork(void)
{
    bool badAlloc = false;
    size_t n = 0;

    (void)allocType;

    srand48(cwID);

    n = nElems[0];
    a1 = static_cast<double *>(memCalloc(allocType, n, sizeof(*a1)));
    if (!a1) { badAlloc = true; goto out; }
    a2 = static_cast<double *>(memCalloc(allocType, n, sizeof(*a2)));
    if (!a2) { badAlloc = true; goto out; }
    a3 = static_cast<double *>(memCalloc(allocType, n, sizeof(*a3)));
    if (!a3) { badAlloc = true; goto out; }
    for (size_t i = 0; i < n; ++i) {
        a1[i] = drand48(); a2[i] = drand48(); a3[i] = 0.0;
    }

    n = nElems[1];
    b1 = static_cast<double *>(memCalloc(allocType, n, sizeof(*b1)));
    if (!b1) { badAlloc = true; goto out; }
    b2 = static_cast<double *>(memCalloc(allocType, n, sizeof(*b2)));
    if (!b2) { badAlloc = true; goto out; }
    b3 = static_cast<double *>(memCalloc(allocType, n, sizeof(*b3)));
    if (!b3) { badAlloc = true; goto out; }
    for (size_t i = 0; i < n; ++i) {
        b1[i] = drand48(); b2[i] = drand48(); b3[i] = 0.0;
    }
out:
    if (badAlloc) {
        throw J7Exception(__FILE__, __LINE__, "*** bad allocation ***");
    }
}

static void
deallocWork(void)
{
    memFree(allocType, a1);
    memFree(allocType, a2);
    memFree(allocType, a3);

    memFree(allocType, b1);
    memFree(allocType, b2);
    memFree(allocType, b3);
}

// only processes that compute enter this routine
static void
doWorkHybrid(size_t nElems, double *cap1, double *cap2, double *cap3)
{
    // do bogus work
    size_t i = 0;
    #pragma omp parallel shared(cap1, cap2, cap3, ompChunk) private(i)
    {
        #pragma omp for schedule(dynamic, ompChunk) nowait
        for (i = 0; i < nElems; ++i) {
            cap3[i] = fabs(cap1[i] / max(cap2[i], cap3[i]) * sin(cos(cap1[i])));
        }
    } // end parallel region
}

// all processes call this routine
static void
doWorkMPIE(size_t nElems, double *cap1, double *cap2, double *cap3)
{
    // do bogus work
    size_t i = 0;
    for (i = 0; i < nElems; ++i) {
        cap3[i] = fabs(cap1[i] / max(cap2[i], cap3[i]) * sin(cos(cap1[i])));
    }
}

static void
doWork(size_t arrayID)
{
    double *cur1 = NULL, *cur2 = NULL, *cur3 = NULL;

    // setup the array pointers based on the array id
    if (0 == arrayID) {
        cur1 = a1; cur2 = a2; cur3 = a3;
    }
    else if (1 == arrayID) {
        cur1 = b1; cur2 = b2; cur3 = b3;
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** bad array ID ***");
    }
    // this is where we setup the environment for the call that's really doing
    // the work. doWork is just a top-level wrapper for the real stuff
    if (ALLOC_SYS == allocType) {
        doWorkMPIE(nElems[arrayID], cur1, cur2, cur3);
    }
    else if (ALLOC_J7 == allocType) {
        // only one process will do the work (1 process will spawn OpenMP
        // threads to load balance the work that's in shared-memory)
        //
        // we don't have to update the pointer here because smID 0 already have
        // the base to the shared-memory segment :-)
        if (0 == smID) {
            quo->bindPush(QUO_BIND_PUSH_OBJ, QUO_OBJ_MACHINE, -1);
            doWorkHybrid(nNodeGlobalElems[arrayID], cur1, cur2, cur3);
            quo->bindPop();
        }
        // REMEMBER: this is NODE barrier! this is here because we need to make
        // sure that updates have flushed before we continue.
        quo->barrier();
    }
    else {
        throw J7Exception(__FILE__, __LINE__, "*** bad alloc type ***");
    }
}

static double
getTime(void)
{
    return MPI_Wtime();
}

static string
allocTypeString(AllocType allocType)
{
    switch (allocType) {
        case ALLOC_SYS:
            return string("sys");
        case ALLOC_J7:
            return string("j7");
        default:
            return string("unknown");
    }
    return string("b@dN3s5");
}

static void
emitSetup(void)
{
    if (0 == cwID) {
        cout << endl << "# load-imbl setup #" << endl;
        cout << "# numpe: " <<  totPEs << endl;
        cout << "# numnodes: " <<  quo->nnodes() << endl;
        cout << "# ncores/node: " <<  quo->ncores() << endl;
        cout << "# nsockets/node: " <<  quo->nsockets() << endl;
        cout << "# nnumanodes/node: " <<  quo->nnumanodes() << endl;
        cout << "# allocator type: " << allocTypeString(allocType) << endl;
        cout << "# omp chunk size: " << ompChunk << endl;
        #pragma omp parallel
        {
        if (0 == omp_get_thread_num()) {
            cout << "# omp nthreads: " << omp_get_num_threads() << endl;
        }
        } // end parallel region
        cout << "# niters: " <<  nIters << endl;
        cout << "# n doubles alloc'd per node: " << getNAllocdOnNode() << endl;
    }
}

static void
recordTime(size_t arrayID, double eTime)
{
    if (0 == arrayID) aTimes.push_back(eTime);
    else if (1 == arrayID) bTimes.push_back(eTime);
}

static void
run(void)
{
    double start = 0.0, end = 0.0;
    // target pointers
    double *t1 = NULL, *t3 = NULL;
    for (size_t iter = 0; iter < nIters; ++iter) {
        if (0 == cwID && (0 == iter % 100)) {
            cout << "starting iteration: " << iter << endl;
        }
        for (size_t arrayID = 0; arrayID < N_WORK_ARRAYS; ++arrayID) {
            start = getTime();
            doWork(arrayID);
            end = getTime();
            recordTime(arrayID, (end - start));
            if (0 == arrayID) {
                t1 = a1; t3 = a3;
            }
            else if (1 == arrayID) {
                t1 = b1; t3 = b3;
            }
            // take the global max b3[0] and replace all b1[0]s with that value
            if (MPI_SUCCESS != MPI_Allreduce(&(t3[0]), &(t1[0]), 1,
                MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)) {
                throw J7Exception(__FILE__, __LINE__, "*** MPI_Allreduce ***");
            }
        }
    }
}

#if 0 // update code if by proc work division needed
static void
workerProcSelection()
{
    if (ALLOC_SYS == allocType) {
        workerProc = true;
    }
    if (ALLOC_J7 == allocType) {
        // nominate workers by evenly distributing them across all available
        // sockets.
        workerProc = quo->autoDistrib(QUO_OBJ_SOCKET, 1);
        cout << cwID << " " << workerProc << endl;
    }
}
#endif

static void
sanityValDump(void)
{
    emitSync();
    if (0 == cwID) {
        cout << "*** array dump ***" << endl;
    }
    emitElems(a1, nElems[0], "a");
    emitElems(b1, nElems[1], "b");
    emitSync();
}

#if 0
static double
getMaxVal(const vector<double> &v)
{
    return *max_element(v.begin(), v.end());
}

static double
getMinVal(const vector<double> &v)
{
    return *min_element(v.begin(), v.end());
}

static void
emitStats(void)
{
    double minTime = getMinVal(aTimes);
    double maxTime = getMaxVal(aTimes);

    emitSync();
    cout << "min: " << minTime << endl;
    cout << "max: " << maxTime << endl;
    emitSync();
}
#endif

static void
emitTotalElapsedTime(double timeInSec)
{
    if (0 == cwID) {
        cout << "Time to Completion: " << timeInSec
             << " s" << endl << endl;
    }
}

static void
usage(void)
{
    cout <<
    "usage: load-imbl [ARGS]\n"
    "   --alloc [sys|j7]\n"
    "   --niters N\n"
    "   --help\n"
        << endl;
}

static void
parseSetupFromArgv(int argc, char **argv)
{
    int c = -1;
    struct option long_opts[] = {
        {"alloc",  required_argument, NULL, 'a'},
        {"niters", required_argument, NULL, 'n'},
        {"help",   no_argument,       NULL, 'h'},
        {NULL,     0,                 NULL,  0 }
    };
    const char *opt_string = "a:n:h";
    while (-1 != (c = getopt_long_only(argc, argv, opt_string,
                                       long_opts, NULL))) {
        switch (c) {
            case 'a': {
                if ("j7" == string(optarg)) allocType = ALLOC_J7;
                else if ("sys" == string(optarg)) allocType = ALLOC_SYS;
                else {
                    cerr << "invalid alloc option... using default" << endl;
                }
                break;
            }
            case 'n':
                nIters = static_cast<size_t>(strtol(optarg, NULL, 10));
                break;
            case 'h':
                if (0 == cwID) usage();
                finiMPI();
                exit(EXIT_SUCCESS);
            default:
                ;
        }
    }
    if (optind < argc) {
        cerr << "*** unrecognized input: " << argv[optind] << endl;
    }
}

int
main(int argc, char **argv)
{
    double start = 0, end = 0.0;
    try {
        initMPI(argc, argv);
        parseSetupFromArgv(argc, argv);
        initQUO();
        //workerProcSelection();
        // set number of elements first so we can allocate enough memory with
        // j7 (if enabled)
        setNElems();
        initJ7();
        emitSetup();
        allocWork();
        start = getTime();
        run();
        end = getTime();
        emitTotalElapsedTime(end - start);
        //emitStats();
        sanityValDump();
        deallocWork();
        finiJ7();
        finiQUO();
        finiMPI();
    }
    catch (J7Exception &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/**
 * TODO
 * stats
 * break up work across more procs
 * variation in run time (overall) - std dev, ave
 * note: very similar to MPI-3 sm windows -- make sure to note that
 */
