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

#ifndef FUNCTIONHHINCLUDE
#define FUNCTIONHHINCLUDE

// ***************************************************************************
// ***************************************************************************
// This class holds information about a function. It is mostly for use with
// the parser.
// ***************************************************************************
// ***************************************************************************

#include <string>
#include <sstream>
#include <vector>
#include <deque>

namespace PP
{
using std::string;
using std::deque;
using std::stringstream;
using std::vector;

enum FuncType {FUNC_};

//class ErrorState;

class Function
{

public:
    Function();
    Function(string nme, bool ext, int na, string ftype,  string fdes);

    // Evaluate the function.
    double evaluate(vector<double> &vd, stringstream &serr, int &ierr,
                    int line_number, int file_line_number,
                    string filename, deque<string> *lines);

    string evaluate(vector<string> &vs, stringstream &serr, int &ierr,
                    int line_number, int file_line_number,
                    string filename, deque<string> *lines);

    // Accessor methods.
    string get_name()        { return name; }
    int    get_num_args()    { return nargs; }
    string get_description() { return description; }
    string get_type()        { return type; }

private:

    void name_err(stringstream &serr, int &ierr,
                  int line_number, int file_line_number,
                  string filename, deque<string> *lines);

    void args_mismatch_err(int nargs_found, int nargs_expected,
                           stringstream &serr, int &ierr,
                           int line_number, int file_line_number,
                           string filename, deque<string> *lines);

    // The name of the function.
    string name;

    // Whether the function is external or internal. External functions
    // are C++ functions like sin(), log(), ... Internal functions
    // are defined within the input to the parser (this feature is not
    // implemented yet).
    bool external;

    // The number of arguments for the function.
    int nargs;

    // A text description of the function.
    string description;

    // The type of function. Allowed types are:
    //     real    double arguments, double result (cos, sin, log, ...)
    //     string  string arguments, string results (strlen, strcat, ...)
    string type;
};


} // End of the PP namespace

#endif
