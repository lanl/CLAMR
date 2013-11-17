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

/**
 * TODO
 * replace MPI_Barrier with something faster if needed
 * over allocate for metadata (hide our stuff in the sm segment)
 * add first touch
 * calloc - fill with 0s
 */

#include "j7.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstdlib>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#define J7_ERR_HERE __FILE__, __LINE__
// named std::pair used to exchange handle/size mappings
#define MMSTR "metaHSPair"
// extra space that we add to the caller's allocation request so we can stash
// our metadata in the shared-memory segment. if the amount of j7 stored data
// changes, please update this value. for now, just over allocate a bit.
#define J7_META_FLUFF_B 128

/* ////////////////////////////////////////////////////////////////////////// */
/* Private Utility Routines */
/* ////////////////////////////////////////////////////////////////////////// */

/* ////////////////////////////////////////////////////////////////////////// */
template <typename T>
static std::string
n2s(T number)
{
    std::ostringstream s;
    s << number;
    return s.str();
}

/* ////////////////////////////////////////////////////////////////////////// */
/* J7Exception Class */
/* ////////////////////////////////////////////////////////////////////////// */

/* ////////////////////////////////////////////////////////////////////////// */
J7Exception::J7Exception(std::string fileName,
                         int lineNo,
                         const std::string &errMsg)
{
    whatString = "[" + fileName + " " + n2s(lineNo) + "] " + errMsg;
}

/* ////////////////////////////////////////////////////////////////////////// */
const char *
J7Exception::what(void) const throw()
{
    return whatString.c_str();
}

/* ////////////////////////////////////////////////////////////////////////// */
/* J7 Class */
/* ////////////////////////////////////////////////////////////////////////// */
J7::J7(MPI_Comm &smComm,
       std::size_t segSize)
{
    using std::string;
    int mpiInitialized = 0;

    if (MPI_SUCCESS != MPI_Initialized(&mpiInitialized)) {
        string estr = "MPI_Initialized failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
    // if the MPI library isn't initialized for us, then bail.
    if (!mpiInitialized) {
        string estr = "MPI is not initialized... Cannot continue with j7.";
        throw J7Exception(J7_ERR_HERE, estr);
    }
    if (MPI_SUCCESS != MPI_Comm_dup(smComm, &comm)) {
        string estr = "MPI_Comm_dup failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
    if (MPI_SUCCESS != MPI_Comm_size(comm, &commSize)) {
        string estr = "MPI_Comm_size failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
    if (MPI_SUCCESS != MPI_Comm_rank(comm, &commRank)) {
        string estr = "MPI_Comm_rank failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
    // the size of the entire memory space
    this->segSize = segSize;
    // figure out if i'm the custodian
    custodian = (0 == commRank) ? true : false;
    // initialize the shared-memory segment -- throws
    smSegInit();
}

/* ////////////////////////////////////////////////////////////////////////// */
J7::~J7(void)
{
    using namespace std;
    barrier();
    smSegFini();
    if (MPI_SUCCESS != MPI_Comm_free(&comm)) {
        string estr = "WARNING: MPI_Comm_free failure";
        cerr << estr << endl;
    }
    if (!hosMap.empty()) {
        string estr = "WARNING: hosMap leaks";
        cerr << estr << endl;
    }
    delete msm;
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::barrier(void)
{
    if (MPI_SUCCESS != MPI_Barrier(comm)) {
        std::string estr = "MPI_Barrier failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::bcast(void *buffer, int count, MPI_Datatype datatype, int root)
{
    if (MPI_SUCCESS != MPI_Bcast(buffer, count, datatype, root, this->comm)) {
        std::string estr = "MPI_Bcast failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::allreduce(void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op op)
{
    if (MPI_SUCCESS != MPI_Allreduce(sendbuf, recvbuf, count,
                                     datatype, op, this->comm)) {
        std::string estr = "MPI_Allreduce failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype)
{
    if (MPI_SUCCESS != MPI_Allgather(sendbuf, sendcount, sendtype,
                                     recvbuf, recvcount, recvtype,
                                     this->comm)) {
        std::string estr = "MPI_Allgather failure";
        throw J7Exception(J7_ERR_HERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::setSMSegName(void)
{
    using namespace std;
    if (custodian) {
        srand((unsigned int)time(NULL));
        int randn = rand() % 1024;
        char *me = NULL;
        if (NULL == (me = getenv("USER"))) {
            me = const_cast<char *>("unknown");
        }
        segName = "j7-" + n2s(randn) + "-" + string(me) + "-" +
                  n2s(getpid()) + ".shm";
        int cStrLen = segName.length() + 1;
        bcast(&cStrLen, 1, MPI_INT, 0);
        bcast(const_cast<char *>(segName.c_str()), cStrLen, MPI_CHAR, 0);
    }
    else {
        int cStrLen = 0;
        bcast(&cStrLen, 1, MPI_INT, 0);
        char *cStr = static_cast<char *>(calloc(cStrLen, sizeof(char)));
        bcast(cStr, cStrLen, MPI_CHAR, 0);
        segName = string(cStr);
        free(cStr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::smSegInit(void)
{
    using namespace boost::interprocess;
    using std::map;
    using std::size_t;
    try {
        // first agree upon a segment name (hopefully unique on the system)
        setSMSegName();
        if (custodian) {
            msm = new managed_shared_memory(create_only, segName.c_str(),
                                            segSize);
            // create a space in shared-memory to exchange MSMHandles
            smHandle = msm->construct<MSMHandle>(MMSTR)[1]();
            barrier();
        }
        else {
            barrier();
            msm = new managed_shared_memory(open_only, segName.c_str());
        }
        barrier();
        // remove (probably unlink) the backing store after all attach
        if (custodian) shared_memory_object::remove(segName.c_str());
        // make sure everything is okay with the new segment
        smSanity();
    }
    catch (std::exception &e) {
        std::string estr("smSegInit Failure: " + std::string(e.what()));
        throw J7Exception(J7_ERR_HERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::smSegFini(void)
{
    if (custodian) {
        // destroy named objects
        msm->destroy<MSMHandle>(MMSTR);
        if (!msm->all_memory_deallocated()) {
            std::cerr << "WARNING: *** j7 detected leaked objects ***"
                      << std::endl;
        }
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::smSanity(void)
{
    if (segSize != msm->get_size()) {
        std::string estr("*** read inconsistency *** expected size " +
                         n2s(segSize) + " B, but got " +
                         n2s(msm->get_size()) + " B");
        throw J7Exception(J7_ERR_HERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void *
J7::myBaseAddr(void *segStart, std::size_t offset)
{
    char *caddr = static_cast<char *>(segStart);
    caddr += offset;
    return static_cast<void *>(caddr);
}

/* ////////////////////////////////////////////////////////////////////////// */
void *
J7::memAlloc(std::size_t size)
{
    using namespace std;
    void *addr = NULL;
    // the global size of the allocation -- what we are actually alloc'ing
    size_t realSize = 0;
    unsigned long long mySize = static_cast<unsigned long long>(size);
    unsigned long long realSizei = 0;
    MSMHandle handle;

    // calculate the global allocation size
    allreduce(&mySize, &realSizei, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM);
    realSize = static_cast<size_t>(realSizei);
    // get the respective allocation sizes. we need this in order to figure out
    // each process offset into the segment.
    unsigned long long *allocSizes = new unsigned long long[commSize];
    allgather(&mySize, 1, MPI_UNSIGNED_LONG_LONG,
              allocSizes, 1, MPI_UNSIGNED_LONG_LONG);
    // reserve enough space for everyone's allocation size
    vector<size_t> allocSizeVector(commSize);
    // copy values into vector
    for (int i = 0; i < commSize; ++i) {
        allocSizeVector[i] = static_cast<size_t>(allocSizes[i]);
    }
    delete[] allocSizes;

    if (custodian) {
        // allocate with 8 byte alignment (throws bad_alloc)
        addr = msm->allocate_aligned(realSize, 8);
        storeMSMHandleinSM(addr);
        handle = retrieveHandle();
        barrier();
    }
    else {
        barrier(); // wait for custodian
        handle = retrieveHandle();
        addr = msm->get_address_from_handle(handle);
    }
    barrier();
    // now fix the offset that each process returns
    size_t myOffset = 0;
    for (int i = 0; i < commRank; ++i) {
        myOffset += allocSizeVector[i];
    }
    // store my offset and my size for this particular handle
    OffSize offSize = {myOffset, size};
    // stash the allocation information
    hosMap[handle] = offSize;
    // fix offset and return to caller
    return myBaseAddr(addr, myOffset);
}

/* ////////////////////////////////////////////////////////////////////////// */
void *
J7::memCalloc(std::size_t nElems, std::size_t elemSize)
{
    std::size_t size = nElems * elemSize;
    // synchronization in memAlloc. throws bad_alloc, so no return check
    char *p = static_cast<char *>(memAlloc(size));
    return memset(p, 0, size);
}

/* ////////////////////////////////////////////////////////////////////////// */
         // XXX this is really expensive. please fix me, kind sir. XXX
void *
J7::memRealloc(void *ptr, std::size_t size)
{
    using namespace std;
    // first create a new shiny memory region
    void *newRegion = memAlloc(size);
    // my view of the base of the old allocation
    void *base = NULL;
    // my old offset and size for the region that is being realloc'd
    OffSize oldOffSize = {0, 0};
    // allocation handle
    MSMHandle handle;

    if (custodian) {
        storeMSMHandleinSM(ptr);
        handle = retrieveHandle();
        oldOffSize = hosMap.find(handle)->second;
        barrier();
    }
    else {
        barrier();
        handle = retrieveHandle();
        oldOffSize = hosMap.find(handle)->second;
    }
    // get base of old allocation
    base = myBaseAddr(msm->get_address_from_handle(handle), oldOffSize.offset);
    // set copy size to the min(oldSize, newSize)
    size_t copySize = min(oldOffSize.size, size);
    // we can copy in parallel - so do that
    (void)memmove(newRegion, base, copySize);
    barrier();
    // done with the old
    memFree(ptr);
    return newRegion;
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::memFree(void *ptr)
{
    MSMHandle handle;
    if (custodian) {
        storeMSMHandleinSM(ptr);
        handle = retrieveHandle();
        barrier();
    }
    else {
        barrier();
        handle = retrieveHandle();
    }
    barrier();
    if (custodian) {
        // only the custodian frees the memory in the shared-memory segment
        msm->deallocate(ptr);
    }
    // everyone stores handle/offset info in the hosMap, so rm the cruft
    hosMap.erase(handle);
    barrier();
}

/* ////////////////////////////////////////////////////////////////////////// */
void *
J7::getBaseMemPtr(void *j7ptr)
{
    void *baseAddr = NULL;
    if (custodian) {
        storeMSMHandleinSM(j7ptr);
        baseAddr = j7ptr;
        barrier();
    }
    else {
        barrier();
        baseAddr = msm->get_address_from_handle(retrieveHandle());
    }
    barrier();
    return baseAddr;
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::storeMSMHandleinSM(void *j7Ptr)
{
    if (!custodian) {
        std::string eStr = "storeMSMHandleinSM can only be called by the custodian";
        throw J7Exception(J7_ERR_HERE, eStr);
    }
    // remember that this is in shared-memory
    *smHandle = msm->get_handle_from_address(j7Ptr);
}

/* ////////////////////////////////////////////////////////////////////////// */
MSMHandle
J7::retrieveHandle(void)
{
    return *(msm->find<MSMHandle>(MMSTR).first);
}
