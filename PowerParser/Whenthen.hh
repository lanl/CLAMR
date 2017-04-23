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

#ifndef WHENTHENHHINCLUDE
#define WHENTHENHHINCLUDE

// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// ***************************************************************************

#include <string>
#include <deque>
#include <vector>
#include <map>
#include <sstream>

#include "Word.hh"

namespace PP
{
using std::string;
using std::deque;
using std::vector;
using std::map;
using std::stringstream;

class Whenthen
{

public:
    Whenthen();
    Whenthen(int &nwhen, Cmd &cmdi, bool &skipwhen, bool &single_line_when,
             bool eflag, stringstream &serr, int &ierr);
    void add_cmdf(Cmd &cmdi);
    void list_condition(string offset1, string offset2,
                        stringstream &ssc);
    void list_cmdsf_ss(stringstream &ssc);

    void check_wt(vector<string> &code_varnames,
                  vector<string> &code_values,
                  vector<int> &vv_active,
                  int *wtci, stringstream &serr, int &ierr);

    deque<Cmd> *get_cmdsf_ptr() { return &cmdsf; }

    void get_char_array_size(int *ca_size);
    void get_char_array(string &sc);

    void get_satsize(int *sat_size);
    void getsat(int *sat);
    void setsat(int *sat);
    void getprocessed(int *wtp);
    void setprocessed(int wtp);
    void getseq(int *wtseq);
    void setseq(int wtseq);
    int get_num_varnames() { return (int)varname.size(); }
    string get_varname(int i) { return varname[i].get_string(); }


private:

    void add_word(Cmd &cmdi, int idex, deque<Word> &wq);
    void add_word(Cmd &cmdi, int idex, deque<Word> &wq, string sadd);
    void process_words(deque <Word> &words, vector<string> &code_varnames,
                       vector<string> &code_values,
                       vector<int> &vv_active,
                       stringstream &serr, int &ierr);
    void delete_words(int i1, int i2, deque <Word> &words);
    void replace_words(int i1, int i2, deque <Word> &words, Word &w);

    // The condition:   varname relation value  logical  varname relation value etc.
    // Example:           time     .gt.   3.0    .and.    ncycle   .ge.    50
    // The condition is thought of as a sequence of subconditions connected by
    // logical operators. The above example has two subconditions connected by the
    // .and. logical operator.
    deque<Word>   varname;    // Host code variable name to be replaced by host code value.
    deque<Word>   relation;   // Relation between varname and value, like .gt., .hglt., ...
    deque<Word>   value;      // Value to compare with host code value.
    deque<Word>   logop;      // Logical operator connecting subconditions.
    deque<string> satisfied;  // Satisfied flag for each subcondition.
    deque<bool>   has_got;    // Has got flag for the relation. This is true if
                              // the relation is .hggt., .hglt., ..., false otherwise.

    // Commands to be done when the condition is satisfied.
    deque<Cmd> cmdsf;

    // The whenthen is only done once when the condition is satisfied.
    // This flag keeps it from being done again.
    bool processed;

    // This flag is used to distinguish between the when command and the
    // whenever command.
    bool ever_flag;

    // This is a sequence index to keep track of what order the whenthen's
    // have been processed in.
    int seqdex;
};


} // end of PP namespace

#endif

