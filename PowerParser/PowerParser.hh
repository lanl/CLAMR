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

// ***************************************************************************
// ***************************************************************************
// Provide a class that parses text files into lines and words.
// ***************************************************************************
// ***************************************************************************
#ifndef PARSEHHINCLUDE
#define PARSEHHINCLUDE

#include <string>
#include <vector>
#include <deque>
#include <map>
#include <sstream>

// Need to include Cmd.hh because on the PGI compiler, the deque<Cmd>
// declaration did not work with just doing "class Cmd;", we need to fully
// include Cmd.hh. 
#include "Word.hh"
#include "Cmd.hh"
#include "Restartblock.hh"
#include "Whenthen.hh"
#include "Comm.hh"

namespace PP
{
using std::string;
using std::vector;
using std::deque;
using std::map;
using std::stringstream;

    //class Variable;
    //class Function;

class PowerParser
{
public:
    
    // Constructors, destructors and drivers.
    PowerParser();
    PowerParser(string filename);
    PowerParser(const char *filename);
    void parse_file(string filename);
    void parse_file(const char *filename);
    void parse_string(string filename, string s_in);
    void compile_buffer();
    void handle_exe_args(string other_argggs);
    void clear_and_init();

    void store_exe_args(string &oargs, string &fname) {
        other_args = oargs;
        file_name = fname;
    }
    void get_exe_args(string &oargs, string &fname) {
        oargs = other_args;
        fname = file_name;
    }

    // String versions
    void get_bool_int(string &cname,
                      int *cvalue,
                      const vector<int> &size = vector<int>(), // optional argument
                      bool skip = false);                      // optional argument
    void get_bool(string &cname,
                  bool *cvalue,
                  const vector<int> &size = vector<int>(),     // optional argument
                  bool skip = false);                          // optional argument
    template< typename T >
    void get_int(string &cname,
                 T *cvalue,
                 const vector<int> &size = vector<int>(),      // optional argument
                 bool skip = false);                           // optional argument
    void get_real(string &cname,
                  double *cvalue,
                  const vector<int> &size = vector<int>(),     // optional argument
                  bool skip = false);                          // optional argument
    void get_char(string &cname,
                  vector<string> &vstr,
                  const vector<int> &size = vector<int>(),     // optional argument
                  bool single_char = false,                    // optional argument
                  bool skip = false);                          // optional argument

    // These are just convenience function to allow char arrays for get variable so
    // the calls are simpler. They convert the cname to a string and call the 
    // string versions above
    void get_bool_int(const char *cname,
                      int *cvalue,
                      const vector<int> &size = vector<int>(), // optional argument
                      bool skip = false);                      // optional argument
    void get_bool(const char *cname,
                  bool *cvalue,
                  const vector<int> &size = vector<int>(),     // optional argument
                  bool skip = false);                          // optional argument
    template< typename T >
    void get_int(const char *cname,
                 T *cvalue,
                 const vector<int> &size = vector<int>(),      // optional argument
                 bool skip = false);                           // optional argument
    void get_real(const char *cname,
                  double *cvalue,
                  const vector<int> &size = vector<int>(),     // optional argument
                  bool skip = false);                          // optional argument
    void get_char(const char *cname,
                  vector<string> &vstr,
                  const vector<int> &size = vector<int>(),     // optional argument
                  bool single_char = false,                    // optional argument
                  bool skip = false);                          // optional argument

    void get_size(string &cname, vector<int> &size);
    void get_sizeb(string &cname, vector<int> &size);
    void cmd_in_input(string &cname, bool &in_input, bool &in_whenthen);
    void cmd_set_processed(string &cname, bool bval);
    void check_processed(bool &good);
    void check_duplicates();

    void list_funcs_start();
    void list_vars_start();
    void list_cmdsf_start();
    void list_wt_cmdsf_start();
    bool get_ssfout_line(string &sline);

    void echo_input_start();
    void echo_input_ss(stringstream &ssinp);

    void process_error_global();


    void rb_check(vector<string> &code_varnames,
                  vector<string> &code_values,
                  vector<int> &vv_active, int *rbci,
                  int *rb_ntriggered, int *rb_triggered_indices);
    int  get_rb_num_varnames();
    void get_rb_varnames(vector<string> &rb_varnames_vstr);
    void get_num_rb(int *rbnum) { *rbnum = (int)restartblocks.size(); }
    void set_num_rb(int rbnum)  { nrb_on_dump = rbnum; }
    void get_rb_names(vector<string> &rb_names_vstr);
    void set_rb_names(vector<string> &rb_names_vstr);
    void get_rb_aflags(int *rb_aflags);
    void set_rb_aflags(int *rb_aflags, int rb_num);
    void get_rb_satsize(int *rb_satsize);
    void set_rb_satsize(int rb_satsize);
    void get_rb_satprb(int *rb_satprb);
    void set_rb_satprb(int *rb_satprb, int rb_num);
    void get_rb_sat(int *rb_sat);
    void set_rb_sat(int *rb_sat, int rb_satsize);
    void list_rb();
    void list_rb_start();
    void list_rb_ss(stringstream &ssc);
    void list_rb1_start(int *rb);
    void list_rb1_ss(stringstream &ssc, int *rbp);
    void list_one_rb_ss(stringstream &ssc, int rb);


    void get_num_whenthen(int *wtnum) { *wtnum = (int)whenthens.size(); }
    void wt_check(int wtn, vector<string> &code_varnames,
                  vector<string> &code_values,
                  vector<int> &vv_active, int *wtci);
    void wt_set_cmdsfp(int wtn);
    void wt_reset();
    void wt_casize(int wtn, int *wt_casize);
    void wt_carray(int wtn, char *wt_ca, int wt_casize);

    void wt_satsize(int wtn, int *wt_satsize);
    void wt_getsat(int wtn, int *wt_sat, int wt_satsize);
    void wt_setsat(int wtn, int *wt_sat, int wt_satsize);
    void wt_getprocessed(int wtn, int *wtp);
    void wt_setprocessed(int wtn, int wtp);
    void wt_getseq(int wtn, int *wtseq);
    void wt_setseq(int wtn, int wtseq);

    void chars_to_vstr(char *chars_1d, vector<string> &vstr,
                       int nv, int nchar);
    void vstr_to_chars(char *chars_1d, vector<string> &vstr,
                       int nv, int nchar);

    // Communications object from the infrastructure. 
    Comm *comm;

private:

    void init();
    void process_dav_cmd();
    void check_dup_scalar(int wtn, bool &found_any);
    void set_dup_row(vector<string> &row, Cmd &cmdi, int iw);
    void remove_dup_scalar(int wtn);
    void read_into_string(string filename, string &s_in);
    void broadcast_buffer(string &s_in);
    bool get_line_from_string(string &strn, string &sout, int &current_pos);
    bool get_sc_line_from_string(string &strn, string &sout, int &current_pos);
    void store_line_strings(string &s_in);
    void eliminate_white_space(string &sline);
    void cmd_set_reprocessed(bool bval);
    void process_error(stringstream &serr, int &ierr);

    void list_vars(string lv1, string lv2, string var_to_list);
    void list_vars_ss(string lv1, string lv2, string var_to_list,
                      stringstream &ssvars);

    void list_funcs(string lf1, string lf2);
    void list_funcs_ss(string lf1, string lf2, stringstream &ssfunc);

    void list_cmdsf(string lc1, string lc2);
    void list_cmdsf_ss(string lc1, string lc2,
                       stringstream &ssc);
    void list_wt_cmdsf();
    void list_wt_cmdsf_ss(stringstream &ssc);

    void print_strings(vector< vector<string> > rows, int n_header_rows,
                       int offset, int col_spacing, int line_len_max,
                       stringstream &ss);
    bool end_do_loop(int &i, deque<int> &do_start,
                     stringstream &serr, int &ierr);
    void end_do_ret(int &i, deque<int> &do_start,
                    stringstream &serr, int &ierr);
    void check_enddo(deque<int> &do_start, stringstream &serr, int &ierr);
    void jump_to_call(int &i, deque<int> &icall, deque<int> &isub,
                      stringstream &serr, int &ierr);
    void jump_to_sub(int &i, string &sub_name,
                     stringstream &serr, int &ierr);
    void print_line(int i);
    void print_line(Cmd &cmd);

    // Store exe line arguments.
    string other_args, file_name;

    // A double ended queue for storing the original lines. This is
    // before the lines get turned into Cmds. 
    // line_number is an index into cmd_strings, note that it starts
    // from 1, not 0.
    deque<string> cmd_strings;
    int line_number;

    // Define a map for a set of variables.
    map<string, Variable> vmap;

    // Define a map for the functions.
    map<string, Function> fmap;

    // A double ended queue for storing the commands.
    deque<Cmd> cmds;
    deque<Cmd> cmdsf;
    deque<Cmd> *cmdsfp;

    // Store cmd names that have been processed, used for clearing and
    // recreating the parser.
    deque<string> processed_cmd_names;

    // Related to writing output to a fortran file.
    int ssfout_current_pos;
    stringstream ssfout;

    // Used for storing the list of pre-defined variables to be printed
    // out later.
    stringstream pre_defined_varss;

    // Used for storing multiple errors and processing them later.
    stringstream serr_global;
    int ierr_global;

    // The execution line arguments are put in this string.
    string exe_args_str;

    // The when ... then objects.
    deque<Whenthen> whenthens;

    // Restart blocks.
    deque<Restartblock> restartblocks;
    int nrb_on_dump;
    deque<string> bnames_on_dump;
    deque<bool> baflags_on_dump;
    int satsize_on_dump;
    deque<bool> rbsat_on_dump;
    deque<int> rbsatprb_on_dump;

    // Flag for whether duplicate array values will be none, fatal, or
    // a warning, determined by the duplicate_array_values command.
    //    dup_fatal = 0     Turn off duplicate array value checking
    //    dup_fatal = 1     Duplicate array value checking is a warning
    //    dup_fatal = 2     Duplicate array value checking is a fatal error
    int dup_fatal;
};

} // end PP namespace

#endif
