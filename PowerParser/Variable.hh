/*
 *  Copyright (c) 2010-2014, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2010-2014. Los Alamos National Security, LLC. This software was produced 
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
 *  Authors: Chuck Wingate   XCP-2   caw@lanl.gov
 *           Robert Robey    XCP-2   brobey@lanl.gov
 */

#ifndef VARIABLEHHINCLUDE
#define VARIABLEHHINCLUDE

// ***************************************************************************
// ***************************************************************************
// This class holds information about a variable.
// ***************************************************************************
// ***************************************************************************

#include <string>
#include <sstream>
#include <vector>
#include <deque>

namespace PP
{
using std::string;
using std::stringstream;
using std::vector;
using std::deque;



class Variable
{

public:
    Variable();
    Variable(int base);
    Variable(string nme, string v, bool pred, string tdes);
    Variable(string nme);
    Variable(string nme, vector<int> &istart, vector<string> &vvec,
             int lnum, int file_lnum, string fname, deque<string> *lines,
             stringstream &serr, int &ierr);

    // Accessor methods.
    string get_varname() { return name; }
    void set_varname(string s) { name = s; }
    int get_ndim() { return ndim; }

    int get_nvalues() { return (int)value.size(); }

    string get_var_value()   { return value[0]; }
    string get_var_value(int idex)   { return value[idex]; }
    string get_var_value(vector<int> &adex, string vname, int lnum,
                         int file_lnum, string fname, deque<string> *lines,
                         stringstream &serr, int &ierr);


    void set_var_value(vector<int> &istart, vector<string> &valvec,
                       int lnum, int file_lnum, string fname,
                       deque<string> *lines, stringstream &serr, int &ierr);
    void bump_var(vector<int> &istart, int increment,
                  int lnum, int file_lnum, string fname,
                  deque<string> *lines, stringstream &serr, int &ierr);


    void set_bounds(vector<int> &bounds, int lnum, int file_lnum,
                    string fname, deque<string> *lines,
                    stringstream &serr, int &ierr);

    void get_indices(int icdex, vector<int> &adex);

    string get_description() { return description; }
    void set_description(string vardes) { description = vardes; }

    bool is_pre_defined() { return pre_defined; }

    void set_temporary(bool b) { temporary = b; }
    bool is_temporary() { return temporary; }

private:

    // name         The name of the variable.
    // value        Vector containing the values of the variable.
    // ndim         Number of dimensions, for example var(9,3) has ndim=2
    // maxdim       Max num for each dimension except the last.
    // lnum_ndim    The line number where ndim was set.
    // lnum_bounds  The line number where maxdim was set.
    // pre_defined  Pre-defined vars cannot be redefined.
    // description  Text description of the variable.
    // temporary    A temporary variable.
    string name;
    vector<string> value;
    int ndim, lnum_bounds, lnum_ndim;
    vector<int> maxdim;
    bool pre_defined, temporary;
    string description;
};


} // End of the PP namespace

#endif
