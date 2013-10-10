/*
 *  Copyright (c) 2011-2013, Los Alamos National Security, LLC.
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
 *  This research code is being developed as part of the 
 *  2011 X Division Summer Workshop for the express purpose
 *  of a collaborative code for development of ideas in
 *  the implementation of AMR codes for Exascale platforms
 *  
 *  AMR implementation of the Wave code previously developed
 *  as a demonstration code for regular grids on Exascale platforms
 *  as part of the Supercomputing Challenge and Los Alamos 
 *  National Laboratory
 *  
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Neal Davis              davis68@lanl.gov, davis68@illinois.edu
 *           David Nicholaeff        dnic@lanl.gov, mtrxknight@aol.com
 *           Dennis Trujillo         dptrujillo@lanl.gov, dptru10@gmail.com
 * 
 */

#ifndef QUO_HPP_INCLUDED
#define QUO_HPP_INCLUDED

#include "quo.h"

#include <iostream>
#include <cstdlib>
#include <string>

/* ////////////////////////////////////////////////////////////////////////// */
// TODO add cpp checks for QUO
// move this to another file
// add proper exception handling
/* ////////////////////////////////////////////////////////////////////////// */
class QUO {
private:
    // the quo context
    QUO_context context;

public:
    /* ////////////////////////////////////////////////////////////////////// */
    QUO(void) { this->context = NULL; }

    /* ////////////////////////////////////////////////////////////////////// */
    ~QUO(void) { this->context = NULL; }

    void create(void);

    void free(void);

    /* ////////////////////////////////////////////////////////////////////// */
    int nnumanodes(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nnumanodes(this->context, &n))) {
            std::cerr << "!!! QUO_nnumanodes failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int nsockets(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nsockets(this->context, &n))) {
            std::cerr << "!!! QUO_nsockets failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int ncores(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_ncores(this->context, &n))) {
            std::cerr << "!!! QUO_ncores failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int npus(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_npus(this->context, &n))) {
            std::cerr << "!!! QUO_pus failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int nnodes(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nnodes(this->context, &n))) {
            std::cerr << "!!! QUO_nnodes failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    bool bound(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_bound(this->context, &n))) {
            std::cerr << "!!! QUO_bound failure !!!" << std::endl;
            return false;
        }
        return (1 == n);
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int id(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_id(this->context, &n))) {
            std::cerr << "!!! QUO_id failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int nqids(void) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nqids(this->context, &n))) {
            std::cerr << "!!! QUO_nqids failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    std::string stringifyCBind(void) {
        char *cbind = NULL;
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != QUO_stringify_cbind(this->context, &cbind)) {
            std::cerr << "!!! QUO_stringify_cbind failure !!!" << std::endl;
            return std::string("?");
        }
        std::string resStr(cbind);
        std::free(cbind); cbind = NULL;
        return resStr;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int nObjsByType(QUO_obj_type_t type) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nobjs_by_type(this->context, type, &n))) {
            std::cerr << "!!! QUO_nobjs_by_type failure !!!" << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int nObjsInType(QUO_obj_type_t inType,
                    int typeIndex,
                    QUO_obj_type_t type) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nobjs_in_type_by_type(this->context,
                                                           inType, typeIndex,
                                                           type, &n))) {
            std::cerr << "!!! QUO_nobjs_in_type_by_type failure !!!"
                      << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    bool cpuSetInType(QUO_obj_type_t inType,
                      int typeIndex) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_cpuset_in_type(this->context,
                                                    inType, typeIndex, &n))) {
            std::cerr << "!!! QUO_cpuset_in_type failure !!!"
                      << std::endl;
            return false;
        }
        return (1 == n);
    }

    /* ////////////////////////////////////////////////////////////////////// */
    int nMachineProcsInType(QUO_obj_type_t inType,
                            int typeIndex) {
        int rc = QUO_ERR, n = 0;
        if (QUO_SUCCESS != (rc = QUO_nmachine_procs_in_type(this->context,
                                                            inType, typeIndex,
                                                            &n))) {
            std::cerr << "!!! QUO_nmachine_procs_in_type failure !!!"
                      << std::endl;
            n = 0;
        }
        return n;
    }

    /* ////////////////////////////////////////////////////////////////////// */
    void bindPush(QUO_bind_push_policy_t policy,
                  QUO_obj_type_t type,
                  int obj_index) {

        int rc = QUO_ERR;
        if (QUO_SUCCESS != (rc = QUO_bind_push(this->context, policy, type,
                                               obj_index))) {
            std::cerr << "!!! QUO_bind_push failure !!!" << std::endl;
        }
    }

    /* ////////////////////////////////////////////////////////////////////// */
    void bindPop(void) {
        int rc = QUO_ERR;
        if (QUO_SUCCESS != (rc = QUO_bind_pop(this->context))) {
            std::cerr << "!!! QUO_bind_pop failure !!!" << std::endl;
        }
    }

    /* ////////////////////////////////////////////////////////////////////// */
    void barrier(void) {
        if (QUO_SUCCESS != QUO_barrier(this->context)) {
            std::cerr << "!!! QUO_barrier failure !!!" << std::endl;
        }
    }

    /* ////////////////////////////////////////////////////////////////////// */
    bool autoDistrib(QUO_obj_type_t distrib_over_this,
                     int max_qids_per_res_type) {
        int isel = 0;
        bool selected = false;
        if (QUO_SUCCESS != QUO_auto_distrib(this->context, distrib_over_this,
                                            max_qids_per_res_type, &isel)) {
            std::cerr << "!!! QUO_auto_distrib failure !!!" << std::endl;
            return false;
        }
        return (1 == isel);
    }

};

#endif //QUO_HPP_INCLUDED
