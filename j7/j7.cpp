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
    this->segSize = segSize;
    custodian = (0 == commRank) ? true : false;
    // throws
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
    delete msm;
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
        segName = "j7-" + n2s(randn) + "-" + string(me) + ".shm";
        int cStrLen = segName.length() + 1;
        bcast(&cStrLen, 1, MPI_INT, 0, comm);
        bcast(const_cast<char *>(segName.c_str()), cStrLen, MPI_CHAR, 0, comm);
    }
    else {
        int cStrLen = 0;
        bcast(&cStrLen, 1, MPI_INT, 0, comm);
        char *cStr = static_cast<char *>(calloc(cStrLen, sizeof(char)));
        bcast(cStr, cStrLen, MPI_CHAR, 0, comm);
        segName = string(cStr);
        free(cStr);
    }
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
J7::bcast(void *buffer, int count, MPI_Datatype datatype,
          int root, MPI_Comm comm)
{
    if (MPI_SUCCESS != MPI_Bcast(buffer, count, datatype, root, comm)) {
        std::string estr = "MPI_Bcast failure";
        throw J7Exception(J7_ERR_HERE, estr);
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
            smHSPair = msm->construct<HSPair>(MMSTR)[1]();
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
    using namespace boost::interprocess;
    using std::cerr;
    using std::endl;
    if (custodian) {
        // destroy named objects
        msm->destroy<HSPair>(MMSTR);
        if (!msm->all_memory_deallocated()) {
            cerr << "WARNING: *** j7 detected leaked objects ***" << endl;
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
J7::memAlloc(std::size_t size)
{
    void *addr = NULL;
    char *caddr = NULL;
    static std::size_t npes = static_cast<std::size_t>(commSize);
    std::size_t realSize = size * npes;

    if (custodian) {
        // 8 byte alighnment
        addr = msm->allocate_aligned(realSize, 8, std::nothrow);
        setHSPair(addr, realSize, true);
        barrier();
    }
    else {
        barrier();
        addr = msm->get_address_from_handle(retrieveHSPair().first);
    }
    // everyone stash the allocation size associated with the new handle
    barrier();
    caddr = static_cast<char *>(addr);
    caddr += (commRank * size);
    return static_cast<void *>(caddr);
}

/* ////////////////////////////////////////////////////////////////////////// */
void *
J7::memCalloc(std::size_t nElems, std::size_t elemSize)
{
    // synchronization in memAlloc
    return memAlloc(nElems * elemSize);
}

/* ////////////////////////////////////////////////////////////////////////// */
         // XXX this is really expensive. please fix me, kind sir. XXX
             // assuming that realloc will only GROW a segment. //
void *
J7::memRealloc(void *ptr, std::size_t size)
{
    using namespace std;
    // first create a new shiny memory region
    void *newRegion = memAlloc(size);
    void *myBase = NULL;
    HSPair hsp;
    if (custodian) {
        setHSPair(ptr);
        hsp = retrieveHSPair();
        myBase = ptr;
        barrier();
    }
    else {
        barrier();
        hsp = retrieveHSPair();
        myBase = msm->get_address_from_handle(hsp.first);
        char *cb = static_cast<char *>(myBase);
        cb += (commRank * hsp.second / commSize);
        myBase = static_cast<void *>(cb);
    }
    // we can do this in parallel - so do that
    (void)memmove(newRegion, myBase, hsp.second / commSize);
    barrier();
    // done with the old
    memFree(ptr);
    return newRegion;
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::memFree(void *ptr)
{
    if (custodian) {
        hsMap.erase(hsMap[this->msm->get_handle_from_address(ptr)]);
        msm->deallocate(ptr);
    }
    barrier();
}

/* ////////////////////////////////////////////////////////////////////////// */
void *
J7::getBaseMemPtr(void *j7ptr)
{
    void *baseAddr = NULL;
    if (custodian) {
        setHSPair(j7ptr);
        baseAddr = j7ptr;
        barrier();
    }
    else {
        barrier();
        baseAddr = msm->get_address_from_handle(retrieveHSPair().first);
    }
    barrier();
    return baseAddr;
}

/* ////////////////////////////////////////////////////////////////////////// */
void
J7::setHSPair(void *j7Ptr, std::size_t size, bool newEntry)
{
    if (!custodian) {
        std::string eStr = "setHSPair can only be called by the custodian";
        throw J7Exception(J7_ERR_HERE, eStr);
    }
    smHSPair->first = msm->get_handle_from_address(j7Ptr);
    if (newEntry) {
        smHSPair->second = size;
        hsMap[smHSPair->first] = size;
    }
    // must already be stashed in our map
    else smHSPair->second = hsMap[smHSPair->first];
}

/* ////////////////////////////////////////////////////////////////////////// */
HSPair
J7::retrieveHSPair(void)
{
    return *(msm->find<HSPair>(MMSTR).first);
}
