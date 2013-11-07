/*
 *  Copyright (c)      2013, Los Alamos National Security, LLC.
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
 *  Author: Samuel K. Gutierrez
 */

#include "QUO.hpp"

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstdlib>

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define CLAMR_QUO_WHERE __FILE__, __LINE__

using namespace std;

/* ////////////////////////////////////////////////////////////////////////// */
/* Private Utility Routines */
/* ////////////////////////////////////////////////////////////////////////// */

/* ////////////////////////////////////////////////////////////////////////// */
template <typename T>
static string
n2s(T number)
{
    ostringstream s;
    s << number;
    return s.str();
}

/* ////////////////////////////////////////////////////////////////////////// */
/* QUOException Class */
/* ////////////////////////////////////////////////////////////////////////// */

/* ////////////////////////////////////////////////////////////////////////// */
QUOException::QUOException(string fileName,
                           int lineNo,
                           const std::string &errMsg)
{
    this->whatString = "[" + fileName + " " + n2s(lineNo) + "] " + errMsg;
}

/* ////////////////////////////////////////////////////////////////////////// */
const char *
QUOException::what(void) const throw()
{
    return this->whatString.c_str();
}

/* ////////////////////////////////////////////////////////////////////////// */
/* QUO Class */
/* ////////////////////////////////////////////////////////////////////////// */

/* ////////////////////////////////////////////////////////////////////////// */
void
QUO::create(void)
{
    int rc = QUO_ERR;
    if (QUO_SUCCESS != (rc = QUO_create(&(this->context)))) {
        string estr = "QUO_create failure: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
QUO::free(void)
{
    int rc = QUO_ERR;
    if (QUO_SUCCESS != (rc = QUO_free(this->context))) {
        string estr = "QUO_free failure: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::nnumanodes(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_nnumanodes(this->context, &n))) {
        string estr = "QUO_nnumanodes: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::nsockets(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_nsockets(this->context, &n))) {
        string estr = "QUO_nsockets: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::ncores(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_ncores(this->context, &n))) {
        string estr = "QUO_ncores: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::npus(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_npus(this->context, &n))) {
        string estr = "QUO_npus: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::nnodes(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_nnodes(this->context, &n))) {
        string estr = "QUO_nnodes: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
bool
QUO::bound(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_bound(this->context, &n))) {
        string estr = "QUO_bound: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return (1 == n);
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::id(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_id(this->context, &n))) {
        string estr = "QUO_id: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::nqids(void) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_nqids(this->context, &n))) {
        string estr = "QUO_nqids: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
string
QUO::stringifyCBind(void) {
    char *cbind = NULL;
    int rc = QUO_ERR;
    if (QUO_SUCCESS != QUO_stringify_cbind(this->context, &cbind)) {
        string estr = "QUO_stringify_cbind: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    std::string resStr(cbind);
    std::free(cbind); cbind = NULL;
    return resStr;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::nObjsByType(QUO_obj_type_t type) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_nobjs_by_type(this->context, type, &n))) {
        string estr = "QUO_nobjs_by_type: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
QUO::nObjsInType(QUO_obj_type_t inType,
                 int typeIndex,
                 QUO_obj_type_t type) {
    int rc = QUO_ERR, n = 0;

    if (QUO_SUCCESS != (rc = QUO_nobjs_in_type_by_type(this->context,
                                                       inType, typeIndex,
                                                       type, &n))) {
        string estr = "QUO_nobjs_in_type_by_type: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return n;
}

/* ////////////////////////////////////////////////////////////////////////// */
// returned resource must be free'd by caller
int *
QUO::qidsInType(QUO_obj_type_t inType,
                int typeIndex) {
    int rc = QUO_ERR;
    int nQIDs = 0;
    int *qids = NULL;
    if (QUO_SUCCESS != (rc = QUO_qids_in_type(this->context,
                                              inType, typeIndex,
                                              &nQIDs, &qids))) {
        string estr = "QUO_qids_in_type: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return qids;
}

/* ////////////////////////////////////////////////////////////////////////// */
bool
QUO::cpuSetInType(QUO_obj_type_t inType,
                  int typeIndex) {
    int rc = QUO_ERR, n = 0;
    if (QUO_SUCCESS != (rc = QUO_cpuset_in_type(this->context,
                                                inType, typeIndex, &n))) {
        string estr = "QUO_cpuset_in_type: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return (1 == n);
}

/* ////////////////////////////////////////////////////////////////////////// */
void
QUO::bindPush(QUO_bind_push_policy_t policy,
              QUO_obj_type_t type,
              int obj_index) {

    int rc = QUO_ERR;
    if (QUO_SUCCESS != (rc = QUO_bind_push(this->context, policy, type,
                                           obj_index))) {
        string estr = "QUO_bind_push: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
QUO::bindPop(void) {
    int rc = QUO_ERR;
    if (QUO_SUCCESS != (rc = QUO_bind_pop(this->context))) {
        string estr = "QUO_bind_pop: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
void
QUO::barrier(void) {
    int rc = QUO_ERR;
    if (QUO_SUCCESS != (rc = QUO_barrier(this->context))) {
        string estr = "QUO_barrier: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
bool
QUO::autoDistrib(QUO_obj_type_t distrib_over_this,
                 int max_qids_per_res_type) {
    int isel = 0, rc = QUO_ERR;
    if (QUO_SUCCESS != (rc = QUO_auto_distrib(this->context, distrib_over_this,
                                              max_qids_per_res_type, &isel))) {
        string estr = "QUO_auto_distrib: rc = " + n2s(rc);
        throw QUOException(CLAMR_QUO_WHERE, estr);
    }
    return (1 == isel);
}
