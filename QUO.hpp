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

#ifndef QUO_HPP_INCLUDED
#define QUO_HPP_INCLUDED

#include "quo.h"

#include <exception>
#include <string>
#include <iostream>
#include <cstdlib>

class QUOException : public std::exception {
private:
    std::string whatString;
    QUOException(void);

public:
    ~QUOException(void) throw() { ; }

    QUOException(std::string fileName,
                 int lineNo,
                 const std::string &errMsg);

    virtual const char *what(void) const throw();
};

class QUO {
private:
    // the quo context
    QUO_context context;

public:
    QUO(void) { this->context = NULL; }

    ~QUO(void) { this->context = NULL; }

    void create(void);

    void free(void);

    int nnumanodes(void);

    int nsockets(void);

    int ncores(void);

    int npus(void);

    int nnodes(void);

    bool bound(void);

    int id(void);

    int nqids(void);

    std::string stringifyCBind(void);

    int nObjsByType(QUO_obj_type_t type);

    int *qidsInType(QUO_obj_type_t inType,
                    int typeIndex);

    int nObjsInType(QUO_obj_type_t inType,
                    int typeIndex,
                    QUO_obj_type_t type);

    bool cpuSetInType(QUO_obj_type_t inType,
                      int typeIndex);

    void bindPush(QUO_bind_push_policy_t policy,
                  QUO_obj_type_t type,
                  int obj_index);

    void bindPop(void);

    void barrier(void);

    bool autoDistrib(QUO_obj_type_t distrib_over_this,
                     int max_qids_per_res_type);
};

#endif //QUO_HPP_INCLUDED
