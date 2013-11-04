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

#ifndef J7_HPP_INCLUDED
#define J7_HPP_INCLUDED

#include <exception>
#include <string>
#include <iostream>
#include <cstdlib>
#include <map>
#include <vector>

#include <stdint.h>
#include "mpi.h"

#include <boost/interprocess/managed_shared_memory.hpp>

// boost managed shared memory type
typedef boost::interprocess::managed_shared_memory MSM;
// boost managed shared memory handle type
typedef boost::interprocess::managed_shared_memory::handle_t MSMHandle;
// holds offset/size pairs
struct OffSize {
    std::size_t offset;
    std::size_t size;
};
// Handle/Offset Map (my offset from my view of the base of the sm segment
typedef std::map<MSMHandle, OffSize> HOSMap;

/* ////////////////////////////////////////////////////////////////////////// */
/* J7Exception Class */
/* ////////////////////////////////////////////////////////////////////////// */
class J7Exception : public std::exception {
private:
    std::string whatString;
    J7Exception(void);

public:
    ~J7Exception(void) throw() { ; }

    J7Exception(std::string fileName,
                int lineNo,
                const std::string &errMsg);

    virtual const char *what(void) const throw();
};

/* ////////////////////////////////////////////////////////////////////////// */
/* J7 Class */
/* ////////////////////////////////////////////////////////////////////////// */
class J7 {
private:
    // cache of communicator size
    int commSize;
    // cache of communicator rank
    int commRank;
    // j7 communicator
    MPI_Comm comm;
    // flag indicating whether or not i'm responsible for the sm segment
    bool custodian;
    // shared segment size
    std::size_t segSize;
    // shared segment name
    std::string segName;
    // points to boost managed shared memory
    MSM *msm;

    MSMHandle *smHandle; // points to something in shared memory
    // XXX change to hash map when you start caring about performance XXX
    HOSMap hosMap;

    J7(void);

    void smSegInit(void);

    void smSegFini(void);

    void smSanity(void);

    void storeMSMHandleinSM(void *j7Ptr);

    MSMHandle retrieveHandle(void);

    void * myBaseAddr(void *segStart, std::size_t offset);

    void setSMSegName(void);

    void barrier(void);

    void bcast(void *buffer, int count, MPI_Datatype datatype, int root);

    void allreduce(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op);

    void allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype);

public:
    J7(MPI_Comm &smComm, std::size_t segSize);

    ~J7(void);

    void *memAlloc(std::size_t size);

    void *memCalloc(std::size_t nElems, std::size_t elemSize);

    void *memRealloc(void *ptr, std::size_t size);

    void memFree(void *ptr);

    void *getBaseMemPtr(void *j7ptr);
};

#endif //J7_HPP_INCLUDED
