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
 *  Authors: Chuck Wingate   XCP-2   caw@lanl.gov
 *           Robert Robey    XCP-2   brobey@lanl.gov
 */

#include "PowerParser.hh"
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <deque>
#include <new>
#include <assert.h>
#include <string.h>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::stringstream;
using PP::PowerParser;

#if !defined FC_FUNC
  #if defined aix || defined hpux || defined IBM
    #define    FC_FUNC(a,A)    a
  #elif defined sparc 
    #define    FC_FUNC(a,A)    a ## _
  #else
    #define    FC_FUNC(a,A)    a ## _
  #endif
#endif

// Pointers to the various objects.
PowerParser *parse;

// We need to compile this file with c++ since it interfaces with
// various c++ classes, but we can't have name mangling since the
// routines in this file are called from fortran, therefore we need
// the extern "C" statement.
#ifdef __cplusplus 
extern "C" {
#endif
    
// Get processor info from the comm object.
void FC_FUNC(parser_comm_info,PARSER_COMM_INFO)(int *mype,
                                                int *npes, int *iope)
{
    // Create the Comm object and initialize it.
    *mype = parse->comm->getProcRank();
    *npes = parse->comm->getNumProcs();
    *iope = parse->comm->getIORank();
}


// Fortran to c++ interface function to create the parser.
// The arguments are:
//    oargs      Execution line arguments
//    oargs_len  The number of characters in oargs.
//    fname      The name of the file to be read and parser.
//    fname_len  The number of characters in fname.
void FC_FUNC(parser_create,PARSER_CREATE)(char *oargs, int *oargs_len,
                                          char *fname, int *fname_len)
{
    // Create the parser.
    parse = new PowerParser();

    // Handle the execution line arguments.
    string other_args(oargs, *oargs_len);
    parse->handle_exe_args(other_args);

    // Read the file, creating cmds and words.
    string filename(fname, *fname_len);
    parse->parse_file(filename);

    parse->store_exe_args(other_args, filename);
}

void FC_FUNC(parser_compile_buffer,PARSER_COMPILE_BUFFER)()
{
    // Compile the buffer.
    parse->compile_buffer();
}

void FC_FUNC(parser_recreate,PARSER_RECREATE)()
{
    //cout << "&&&&&cw parser_iface.cc, parser_recreate" << endl;
    // Clear out the parser and re-initialize.
    parse->clear_and_init();

    string other_args, filename;
    parse->get_exe_args(other_args, filename);

    // Handle the execution line arguments.
    parse->handle_exe_args(other_args);

    // Read the file, creating cmds and words.
    parse->parse_file(filename);
}


// Destroy the parser to recover memory.
void FC_FUNC(parser_destroy,PARSER_DESTROY)()
{
    if (parse->comm->isIOProc() && parse->fileout.is_open()) {
       //cout << "Resetting IO to screen" << endl;
       parse->fileout.close();
       cout.rdbuf(parse->coutbuf);
       //cout << "Returning output to screen" << endl;
    }

    // Why does the code fail to build when this is uncommented??
    //delete(parse);

    // Do not delete comm since someday it might be used during the
    // entire calculation.
}


//+***************************************************************************
// ***************************************************************************
// Fortran,C++ interface, get logical values.
// ***************************************************************************
// ***************************************************************************

// Get a single boolean value as an integer.
void FC_FUNC(get_logical0,GET_LOGICAL0)(char *cmdname, int *cmdvalue,
                                        int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_bool_int(cname, cmdvalue, size, skip);
}


// Get a boolean one dimensional array as an integer array.
void FC_FUNC(get_logical1,GET_LOGICAL1)(char *cmdname, int *array_values,
                                        int *array_size, int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*array_size);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_bool_int(cname, array_values, size, skip);
}


// Get a boolean two dimensional array as an integer array.
void FC_FUNC(get_logical2,GET_LOGICAL2)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *len,
                                        int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_bool_int(cname, array_values, size, skip);
}


// Get a boolean three dimensional array as an integer array.
void FC_FUNC(get_logical3,GET_LOGICAL3)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_bool_int(cname, array_values, size, skip);
}


// Get a boolean four dimensional array as an integer array.
void FC_FUNC(get_logical4,GET_LOGICAL4)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *size4, int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_bool_int(cname, array_values, size, skip);
}


//+***************************************************************************
// ***************************************************************************
// Fortran,C++ interface, get integer values.
// ***************************************************************************
// ***************************************************************************

// Get a single integer value.
void FC_FUNC(get_integer0,GET_INTEGER0)(char *cmdname, int *cmdvalue,
                                        int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, cmdvalue, size, skip);
}

// Get a one dimensional integer array.
void FC_FUNC(get_integer1,GET_INTEGER1)(char *cmdname, int *array_values,
                                        int *array_size, int *len,
                                        int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*array_size);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, array_values, size, skip);
}


// Get a two dimensional integer array.
void FC_FUNC(get_integer2,GET_INTEGER2)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *len,
                                        int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, array_values, size, skip);
}


// Get a three dimensional integer array.
void FC_FUNC(get_integer3,GET_INTEGER3)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, array_values, size, skip);
}


// Get a four dimensional integer array.
void FC_FUNC(get_integer4,GET_INTEGER4)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *size4, int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, array_values, size, skip);
}

// Get a single int64 value.
void FC_FUNC(get_int8_0,GET_INT8_0)(char *cmdname, int64_t *cmdvalue,
                                    int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, cmdvalue, size, skip);
}
    
// Get a one dimensional int64 array.
void FC_FUNC(get_int8_1,GET_INT8_1)(char *cmdname, int64_t *array_values,
                                    int *array_size, int *len,
                                    int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*array_size);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_int(cname, array_values, size, skip);
}


//+***************************************************************************
// ***************************************************************************
// Fortran,C++ interface, get real values.
// ***************************************************************************
// ***************************************************************************

// Get a single real value.
void FC_FUNC(get_real0,GET_REAL0)(char *cmdname, double *cmdvalue,
                                  int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_real(cname, cmdvalue, size, skip);
}

// Get a one dimensional real array.
void FC_FUNC(get_real1,GET_REAL1)(char *cmdname, double *array_values,
                                  int *array_size, int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*array_size);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_real(cname, array_values, size, skip);
}


// Get a two dimensional real array.
void FC_FUNC(get_real2,GET_REAL2)(char *cmdname, double *array_values,
                                  int *size1, int *size2, int *len,
                                  int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_real(cname, array_values, size, skip);
}


// Get a three dimensional real array.
void FC_FUNC(get_real3,GET_REAL3)(char *cmdname, double *array_values,
                                  int *size1, int *size2, int *size3,
                                  int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_real(cname, array_values, size, skip);
}


// Get a four dimensional real array.
void FC_FUNC(get_real4,GET_REAL4)(char *cmdname, double *array_values,
                                  int *size1, int *size2, int *size3,
                                  int *size4, int *len, int *iskip)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);
    bool skip = false;
    if ((*iskip) == 1) skip = true;
    parse->get_real(cname, array_values, size, skip);
}


//+***************************************************************************
// ***************************************************************************
// Fortran,C++ interface, get character values.
// ***************************************************************************
// ***************************************************************************

// Get a single character
void FC_FUNC(get_char00,GET_CHAR00)(char *cmdname, char *cmdvalue,
                                    int *len, int *iskip)
{
    bool single_char = true;
    string cname(cmdname, *len);
    //string cvalue(*cmdvalue, 1);   // Wrong
    string cvalue(cmdvalue, 1);
    vector<string> vstr;
    vstr.push_back(cvalue);
    vector<int> size;
    bool skip = false;
    if ((*iskip) == 1) skip = true;

    parse->get_char(cname, vstr, size, single_char, skip);
    cvalue = vstr[0];

    *cmdvalue = cvalue[0];
}

// Get a character string.
void FC_FUNC(get_char0,GET_CHAR0)(char *cmdname, char *cmdvalue,
                                  int *len_name, int *len_value, int *iskip, int *len_value_actual)
{
    bool single_char = false;
    string cname(cmdname, *len_name);
    string cvalue(cmdvalue, *len_value);
    vector<string> vstr;
    vstr.push_back(cvalue);
    vector<int> size;
    bool skip = false;
    if ((*iskip) == 1) skip = true;

    parse->get_char(cname, vstr, size, single_char, skip);

    cvalue = vstr[0];

    int len = (int)cvalue.size();
    *len_value_actual = len;

    if (len > (*len_value)) len = (*len_value);
    for (int i=0; i<len; i++) {
        cmdvalue[i] = cvalue[i];
    }
    for (int i=len; i<(*len_value); i++) {
        cmdvalue[i] = ' ';
    }
}


// Get a 1d array of character strings.
void FC_FUNC(get_char1,GET_CHAR1)(char *cmdname, char *array_values,
                                  int *nchar, int *size1,
                                  int *len_name, int *iskip)
{
    bool single_char = false;

    // Turn the command name into a string.
    string cname(cmdname, *len_name);

    // Turn array_values into a vector of strings.
    vector<string> vstr;
    char *cnchar = new char[(*nchar)];
    for (int i=0; i<(*size1); i++) {
        int istart = i * (*nchar);
        for (int c=istart; c<istart+(*nchar); c++) {
            cnchar[c-istart] = array_values[c];
        }
        string s(cnchar,(*nchar));
        vstr.push_back(s);
    }
    delete [] cnchar;

    // Store array sizes in a vector.
    vector<int> size;
    size.push_back(*size1);

    bool skip = false;
    if ((*iskip) == 1) skip = true;

    // Get the array of strings.
    parse->get_char(cname, vstr, size, single_char, skip);

    // Copy the strings back into the char array.
    for (int i=0; i<(*size1); i++) {
        int istart = i * (*nchar);
        int nc = (int)vstr[i].size();
        if (nc > (*nchar)) nc = (*nchar);
        for (int c=istart; c<istart+nc; c++) {
            array_values[c] = vstr[i][c-istart];
        }
        for (int c=istart+nc; c<istart+(*nchar); c++) {
            array_values[c] = ' ';
        }
    }

    // Debug.
    //for (int i=0; i<(*size1); i++) {
    //    int istart = i * (*nchar);
    //    for (int j=istart; j<istart+(*nchar); j++) {
    //        cout << array_values[j];
    //    }
    //    cout << endl;
    //}
}


// Get a 2d array of character strings.
void FC_FUNC(get_char2,GET_CHAR2)(char *cmdname, char *array_values,
                                  int *nchar, int *size1, int *size2,
                                  int *len_name, int *iskip)
{
    bool single_char = false;

    // Turn the command name into a string.
    string cname(cmdname, *len_name);

    // Turn array_values into a vector of strings.
    vector<string> vstr;
    char *cnchar = new char[(*nchar)];
    for (int j=0; j<(*size2); j++) {
        for (int i=0; i<(*size1); i++) {
            int istart = (i + j*(*size1)) * (*nchar);
            for (int c=istart; c<istart+(*nchar); c++) {
                cnchar[c-istart] = array_values[c];
            }
            string s(cnchar,(*nchar));
            vstr.push_back(s);
        }
    }
    delete [] cnchar;

    // Store array sizes in a vector.
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);

    bool skip = false;
    if ((*iskip) == 1) skip = true;

    // Get the array of strings.
    parse->get_char(cname, vstr, size, single_char, skip);

    // Copy the strings back into the char array.
    for (int j=0; j<(*size2); j++) {
        for (int i=0; i<(*size1); i++) {
            int i1d = i + j*(*size1);
            int istart = i1d * (*nchar);
            int nc = (int)vstr[i1d].size();
            if (nc > (*nchar)) nc = (*nchar);
            for (int c=istart; c<istart+nc; c++) {
                array_values[c] = vstr[i1d][c-istart];
            }
            for (int c=istart+nc; c<istart+(*nchar); c++) {
                array_values[c] = ' ';
            }
        }
    }
}


// Get a 3d array of character strings.
void FC_FUNC(get_char3,GET_CHAR3)(char *cmdname, char *array_values,
                                  int *nchar, int *size1, int *size2,
                                  int *size3, int *len_name, int *iskip)
{
    bool single_char = false;

    // Turn the command name into a string.
    string cname(cmdname, *len_name);

    // Turn array_values into a vector of strings.
    vector<string> vstr;
    char *cnchar = new char[(*nchar)];
    for (int k=0; k<(*size3); k++) {
        for (int j=0; j<(*size2); j++) {
            for (int i=0; i<(*size1); i++) {
                int i1d = i + j*(*size1) + k*(*size1)*(*size2);
                int istart = i1d * (*nchar);
                for (int c=istart; c<istart+(*nchar); c++) {
                    cnchar[c-istart] = array_values[c];
                }
                string s(cnchar,(*nchar));
                vstr.push_back(s);
            }
        }
    }
    delete [] cnchar;

    // Store array sizes in a vector.
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);

    bool skip = false;
    if ((*iskip) == 1) skip = true;

    // Get the array of strings.
    parse->get_char(cname, vstr, size, single_char, skip);

    // Copy the strings back into the char array.
    for (int k=0; k<(*size3); k++) {
        for (int j=0; j<(*size2); j++) {
            for (int i=0; i<(*size1); i++) {
                int i1d = i + j*(*size1) + k*(*size1)*(*size2);
                int istart = i1d * (*nchar);
                int nc = (int)vstr[i1d].size();
                if (nc > (*nchar)) nc = (*nchar);
                for (int c=istart; c<istart+nc; c++) {
                    array_values[c] = vstr[i1d][c-istart];
                }
                for (int c=istart+nc; c<istart+(*nchar); c++) {
                    array_values[c] = ' ';
                }
            }
        }
    }
}


// Get a 4d array of character strings.
void FC_FUNC(get_char4,GET_CHAR4)(char *cmdname, char *array_values,
                                  int *nchar, int *size1, int *size2,
                                  int *size3, int *size4,
                                  int *len_cmdname, int *iskip)
{
    bool single_char = false;

    // Turn the command name into a string.
    string cname(cmdname, *len_cmdname);

    // Turn array_values into a vector of strings.
    vector<string> vstr;
    char *cnchar = new char[(*nchar)];
    for (int l=0; l<(*size4); l++) {
        for (int k=0; k<(*size3); k++) {
            for (int j=0; j<(*size2); j++) {
                for (int i=0; i<(*size1); i++) {
                    int i1d = i + j*(*size1) + k*(*size1)*(*size2) +
                        l*(*size1)*(*size2)*(*size3);
                    int istart = i1d * (*nchar);
                    for (int c=istart; c<istart+(*nchar); c++) {
                        cnchar[c-istart] = array_values[c];
                    }
                    string s(cnchar,(*nchar));
                    vstr.push_back(s);
                }
            }
        }
    }
    delete [] cnchar;

    // Store array sizes in a vector.
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);

    bool skip = false;
    if ((*iskip) == 1) skip = true;

    // Get the array of strings.
    parse->get_char(cname, vstr, size, single_char, skip);

    // Copy the strings back into the char array.
    for (int l=0; l<(*size4); l++) {
        for (int k=0; k<(*size3); k++) {
            for (int j=0; j<(*size2); j++) {
                for (int i=0; i<(*size1); i++) {
                    int i1d = i + j*(*size1) + k*(*size1)*(*size2) +
                        l*(*size1)*(*size2)*(*size3);
                    int istart = i1d * (*nchar);
                    int nc = (int)vstr[i1d].size();
                    if (nc > (*nchar)) nc = (*nchar);
                    for (int c=istart; c<istart+nc; c++) {
                        array_values[c] = vstr[i1d][c-istart];
                    }
                    for (int c=istart+nc; c<istart+(*nchar); c++) {
                        array_values[c] = ' ';
                    }
                }
            }
        }
    }
}



//+***************************************************************************
// ***************************************************************************
// Fortran,C++ interface, get sizes.
// ***************************************************************************
// ***************************************************************************
void FC_FUNC(get_size1,GET_SIZE1)(char *cmdname, int *array_size, int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*array_size);
    parse->get_size(cname, size);
    *array_size = size[0];
}

void FC_FUNC(get_size2,GET_SIZE2)(char *cmdname, int *size1,
                                  int *size2, int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    parse->get_size(cname, size);
    *size2 = size[1];
}

void FC_FUNC(get_sizeb2,GET_SIZEB2)(char *cmdname, int *size1,
                                    int *size2, int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    parse->get_sizeb(cname, size);
    *size1 = size[0];
    *size2 = size[1];
}

void FC_FUNC(get_size3,GET_SIZE3)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    parse->get_size(cname, size);
    *size3 = size[2];
}


void FC_FUNC(get_size4,GET_SIZE4)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *size4,
                                  int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);
    parse->get_size(cname, size);
    *size4 = size[3];
}


void FC_FUNC(get_size5,GET_SIZE5)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *size4,
                                  int *size5, int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);
    size.push_back(*size5);
    parse->get_size(cname, size);
    *size5 = size[4];
}


void FC_FUNC(get_size6,GET_SIZE6)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *size4,
                                  int *size5, int *size6, int *len)
{
    string cname(cmdname, *len);
    vector<int> size;
    size.push_back(*size1);
    size.push_back(*size2);
    size.push_back(*size3);
    size.push_back(*size4);
    size.push_back(*size5);
    size.push_back(*size6);
    parse->get_size(cname, size);
    *size6 = size[5];
}


//+***************************************************************************
// ***************************************************************************
// Fortran,C++ interface, various miscellaneous things.
// ***************************************************************************
// ***************************************************************************

// Final error process to print collection of errors.
void FC_FUNC(process_error_final,PROCESS_ERROR_FINAL)()
{
    parse->process_error_global();
}


// Check processed flags, it is a fatal error if everything is not processed.
void FC_FUNC(check_processed,CHECK_PROCESSED)(int *good_int)
{
    bool good = true;
    if (*good_int == 0) good = false;

    parse->check_processed(good);

    *good_int = 0;
    if (good) *good_int = 1;
}


// Check if the input command, cname, appears in the final, parsed user input.
//
// The two outputs are in_input_i and in_whenthen_i (these are integers for
// easy passing to/from fortran).
//    in_input_i     command is in (or not) the main part of the input, i.e.
//                   everything except the when...then statements.
//    in_whenthen_i  command is in (or not) at least one when...then statement.
//
void FC_FUNC(cmd_in_input,CMD_IN_INPUT)(char *cmdname, int *cmdlen,
                                        int *in_input_i, int *in_whenthen_i)
{
    string cname(cmdname, *cmdlen);

    bool in_input = false;
    bool in_whenthen = false;
    parse->cmd_in_input(cname, in_input, in_whenthen);

    *in_input_i = 0;
    if (in_input) *in_input_i = 1;
    *in_whenthen_i = 0;
    if (in_whenthen) *in_whenthen_i = 1;
}


// Set the processed flag for all words for all commands that match cmdname.
// The value to set the processed flag to is bval.
// This sets the processed flag for commands in the final buffer and in the
// when...then final buffers.
// This is meant to be called from an mgname call.
void FC_FUNC(cmd_set_processed,CMD_SET_PROCESSED)(char *cmdname, int *cmdlen,
                                                  int *bval_i)
{
    string cname(cmdname, *cmdlen);

    bool bval = false;
    if (*bval_i == 1) bval = true;

    parse->cmd_set_processed(cname, bval);
}


// Check for duplicate scalar values.
void FC_FUNC(parser_chk_scalar_dup,PARSER_CHK_SCALAR_DUP)()
{
    parse->check_duplicates();
}


//+***************************************************************************
// ***************************************************************************
// List variables, functions, etc, to a fortran file.
// One problem with a fortran/c++ code is that the c++ cannot write directly
// to a fortran file. To get around this, we have the c++ write to a
// stringstream, then put the fortran into a loop getting line after line
// from the stringstream, and the fortran then writes each line to the file.
// ***************************************************************************
// ***************************************************************************

// Write the functions output to a stringstream.
void FC_FUNC(list_functions_start,LIST_FUNCTIONS_START)()
{
    parse->list_funcs_start();
}

// Write the variables output to a stringstream.
void FC_FUNC(list_variables_start,LIST_VARIABLES_START)()
{
    parse->list_vars_start();
}

// Echo user input to a stringstream
void FC_FUNC(echo_ui_start,ECHO_UI_START)()
{
    parse->echo_input_start();
}

// Echo final buffer to a stringstream
void FC_FUNC(echo_final_buffer,ECHO_FINAL_BUFFER)()
{
    parse->list_cmdsf_start();
}

// Echo whenthen final buffer to a stringstream
void FC_FUNC(echo_wt_final_buffer,ECHO_WT_FINAL_BUFFER)()
{
    parse->list_wt_cmdsf_start();
}

// Get a line from the stringstream.
void FC_FUNC(get_output_line,GET_OUTPUT_LINE)(char *fline, int *nchar,
                                              int *bint)
{
    // Get the line from the parser as a string.
    // If b is false, then there are no more lines to get.
    string sline = "";
    bool b = parse->get_ssfout_line(sline);
    int nc_max;

    // Return an integer to the fortran because fortran logical and c++
    // bool are not compatible.
    *bint = 0;
    if (b) *bint = 1;

    // Number of characters in the line.
    int nc = (int)sline.size();
    // max number of chars to write to fline
    nc_max = nc < *nchar ? nc : *nchar;

    // Copy the characters to the fortran character string. Pad with blanks.
    for (int i=0; i<nc_max; i++) {
        fline[i] = sline.c_str()[i];
    }
    for (int i=nc_max; i<(*nchar); i++) {
        fline[i] = ' ';
    }
}


//+***************************************************************************
// ***************************************************************************
// Fortran restart_block interface
// ***************************************************************************
// ***************************************************************************

// Get the number of restart blocks
void FC_FUNC(parser_get_rbnum,PARSER_GET_RBNUM)(int *rbnum)
{
    parse->get_num_rb(rbnum);
}
void FC_FUNC(parser_set_rbnum,PARSER_SET_RBNUM)(int *rbnum)
{
    parse->set_num_rb(*rbnum);
}


void FC_FUNC(parser_set_rb_names,PARSER_SET_RB_NAMES)(char *rb_names_1d,
    int *rb_num, int *ncstr_rb)
{
    vector<string> rb_names_vstr;
    parse->chars_to_vstr(rb_names_1d, rb_names_vstr, *rb_num, *ncstr_rb);
    parse->set_rb_names(rb_names_vstr);
}

void FC_FUNC(parser_get_rb_names,PARSER_GET_RB_NAMES)(char *rb_names_1d,
    int *rb_num, int *ncstr_rb)
{
    // To suppress compiler warnings of unused parameters
    assert(rb_num == rb_num);

    vector<string> rb_names_vstr;
    parse->get_rb_names(rb_names_vstr);
    parse->vstr_to_chars(rb_names_1d, rb_names_vstr, *rb_num, *ncstr_rb);
}

void FC_FUNC(parser_set_rb_aflags,PARSER_SET_RB_AFLAGS)(int *rb_aflags,
    int *rb_num)
{
    parse->set_rb_aflags(rb_aflags, *rb_num);
}

void FC_FUNC(parser_get_rb_aflags,PARSER_GET_RB_AFLAGS)(int *rb_aflags)
{
    parse->get_rb_aflags(rb_aflags);
}


void FC_FUNC(parser_get_rb_num_varnames,PARSER_GET_RB_NUM_VARNAMES)(int *vnum)
{
    *vnum = parse->get_rb_num_varnames();
}
void FC_FUNC(parser_get_rb_varnames,PARSER_GET_RB_VARNAMES)(char *rb_varnames_1d,
    int *rb_vnum, int *ncstr_rb)
{
    // To suppress compiler warnings of unused parameters
    assert(rb_vnum == rb_vnum);

    vector<string> rb_varnames_vstr;
    parse->get_rb_varnames(rb_varnames_vstr);
    parse->vstr_to_chars(rb_varnames_1d, rb_varnames_vstr, *rb_vnum, *ncstr_rb);
}




// Check restart block triggers.
void FC_FUNC(parser_rb_check,PARSER_RB_CHECK)(char *varname_array,
    char *val_array, int *vv_active_a, int *nchar, int *nvv, int *rbci,
    int *rb_ntriggered, int *rb_triggered_indices)
{
    // Turn variable names into vector of strings.
    vector<string> varname_vstr;
    char *cnchar = new char[(*nchar)];
    for (int i=0; i<(*nvv); i++) {
        int istart = i * (*nchar);
        for (int c=istart; c<istart+(*nchar); c++) {
            cnchar[c-istart] = varname_array[c];
        }
        int cnchar_len = (*nchar);
        for (int c=(*nchar)-1; c >= 0; c--) {
            if (cnchar[c] != ' ') {
                cnchar_len = c+1;
                break;
            }
        }
        string s(cnchar,cnchar_len);
        int i2=0;
        for (int c=0; c<(int)s.size(); c++) {
            if (s[c] != ' ') {
                i2=c;
                break;
            }
        }
        if (i2 != 0) s.erase(s.begin(), s.begin()+i2);
        varname_vstr.push_back(s);
    }
    delete [] cnchar;

    // Turn array values into vector of strings.
    vector<string> val_vstr;
    cnchar = new char[(*nchar)];
    for (int i=0; i<(*nvv); i++) {
        int istart = i * (*nchar);
        for (int c=istart; c<istart+(*nchar); c++) {
            cnchar[c-istart] = val_array[c];
        }
        int cnchar_len = (*nchar);
        for (int c=(*nchar)-1; c >= 0; c--) {
            if (cnchar[c] != ' ') {
                cnchar_len = c+1;
                break;
            }
        }
        string s(cnchar,cnchar_len);
        int i2=0;
        for (int c=0; c<(int)s.size(); c++) {
            if (s[c] != ' ') {
                i2=c;
                break;
            }
        }
        if (i2 != 0) s.erase(s.begin(), s.begin()+i2);
        val_vstr.push_back(s);
    }
    delete [] cnchar;

    // Turn array values into vector of strings.
    vector<int> vv_active;
    for (int i=0; i<(*nvv); i++) {
        vv_active.push_back(vv_active_a[i]);
    }

    parse->rb_check(varname_vstr, val_vstr, vv_active, rbci,
                    rb_ntriggered, rb_triggered_indices);
}


// Get/set the total number of restart block sub-conditions
void FC_FUNC(parser_get_rb_satsize,PARSER_GET_RB_SATSIZE)(int *rb_satsize)
{
    parse->get_rb_satsize(rb_satsize);
}
void FC_FUNC(parser_set_rb_satsize,PARSER_SET_RB_SATSIZE)(int *rb_satsize)
{
    parse->set_rb_satsize(*rb_satsize);
}


// Get/set the number of sub-conditions per restart block
void FC_FUNC(parser_get_rb_satprb,PARSER_GET_RB_SATPRB)(int *rb_satprb)
{
    parse->get_rb_satprb(rb_satprb);
}

void FC_FUNC(parser_set_rb_satprb,PARSER_SET_RB_SATPRB)(int *rb_satprb,
                                                        int *rb_num)
{
    parse->set_rb_satprb(rb_satprb, *rb_num);
}


// Get/set the satisfied flag for each sub-condition for each restart block
void FC_FUNC(parser_get_rb_sat,PARSER_GET_RB_SAT)(int *rb_sat)
{
    parse->get_rb_sat(rb_sat);
}


void FC_FUNC(parser_set_rb_sat,PARSER_SET_RB_SAT)(int *rb_sat, int *rb_satsize)
{
    parse->set_rb_sat(rb_sat, *rb_satsize);
}


// Echo restart block information to a stringstream
void FC_FUNC(echo_rb_info,ECHO_RB_INFO)()
{
    parse->list_rb_start();
}
void FC_FUNC(echo_rb1_info,ECHO_RB1_INFO)(int *rb)
{
    parse->list_rb1_start(rb);
}



//+***************************************************************************
// ***************************************************************************
// Fortran when...then interface
// ***************************************************************************
// ***************************************************************************

// Get the number of when...then commands.
void FC_FUNC(parser_get_wtnum,PARSER_GET_WTNUM)(int *wtnum)
{
    parse->get_num_whenthen(wtnum);
}

void FC_FUNC(parser_wt_check2,PARSER_WT_CHECK2)(int *wti, char *varname_array,
    char *val_array, int *vv_active_a, int *nchar, int *nvv, int *wtci)
{
    int wtn = (*wti);

    // Turn array values into vector of strings.
    vector<string> varname_vstr;
    char *cnchar = new char[(*nchar)];
    for (int i=0; i<(*nvv); i++) {
        int istart = i * (*nchar);
        for (int c=istart; c<istart+(*nchar); c++) {
            cnchar[c-istart] = varname_array[c];
        }
        int cnchar_len = (*nchar);
        for (int c=(*nchar)-1; c >= 0; c--) {
            if (cnchar[c] != ' ') {
                cnchar_len = c+1;
                break;
            }
        }
        string s(cnchar,cnchar_len);
        int i2=0;
        for (int c=0; c<(int)s.size(); c++) {
            if (s[c] != ' ') {
                i2=c;
                break;
            }
        }
        if (i2 != 0) s.erase(s.begin(), s.begin()+i2);
        varname_vstr.push_back(s);
    }
    delete [] cnchar;


    vector<string> val_vstr;
    cnchar = new char[(*nchar)];
    for (int i=0; i<(*nvv); i++) {
        int istart = i * (*nchar);
        for (int c=istart; c<istart+(*nchar); c++) {
            cnchar[c-istart] = val_array[c];
        }
        int cnchar_len = (*nchar);
        for (int c=(*nchar)-1; c >= 0; c--) {
            if (cnchar[c] != ' ') {
                cnchar_len = c+1;
                break;
            }
        }
        string s(cnchar,cnchar_len);
        int i2=0;
        for (int c=0; c<(int)s.size(); c++) {
            if (s[c] != ' ') {
                i2=c;
                break;
            }
        }
        if (i2 != 0) s.erase(s.begin(), s.begin()+i2);
        val_vstr.push_back(s);
    }
    delete [] cnchar;

    // Turn array values into vector of strings.
    vector<int> vv_active;
    for (int i=0; i<(*nvv); i++) {
        vv_active.push_back(vv_active_a[i]);
    }

    parse->wt_check(wtn, varname_vstr, val_vstr, vv_active, wtci);
}


// Set the commands final buffer pointer.
void FC_FUNC(parser_wt_set_cmdsfp,PARSER_WT_SET_CMDSFP)(int *wti)
{
    int wtn = (*wti);
    parse->wt_set_cmdsfp(wtn);
}


// Reset the commands final buffer pointer.
void FC_FUNC(parser_wt_reset,PARSER_WT_RESET)()
{
    parse->wt_reset();
}

void FC_FUNC(parser_wt_casize,PARSER_WT_CASIZE)(int *wti, int *wt_casize)
{
    int wtn = (*wti);
    parse->wt_casize(wtn, wt_casize);
}

void FC_FUNC(parser_wt_ca,PARSER_WT_CA)(int *wti, char *wt_ca, int *wt_casize)
{
    int wtn = (*wti);
    int casize = (*wt_casize);
    parse->wt_carray(wtn, wt_ca, casize);
}


void FC_FUNC(parser_wt_satsize,PARSER_WT_SATSIZE)(int *wti, int *wt_satsize)
{
    int wtn = (*wti);
    parse->wt_satsize(wtn, wt_satsize);
}

void FC_FUNC(parser_wt_getsat,PARSER_WT_GETSAT)(int *wti, int *wt_sat, int *wt_satsize)
{
    int wtn = (*wti);
    int satsize = (*wt_satsize);
    parse->wt_getsat(wtn, wt_sat, satsize);
}

void FC_FUNC(parser_wt_setsat,PARSER_WT_SETSAT)(int *wti, int *wt_sat, int *wt_satsize)
{
    int wtn = (*wti);
    int satsize = (*wt_satsize);
    parse->wt_setsat(wtn, wt_sat, satsize);
}


// Get and Set the processed flag for a whenthen.
void FC_FUNC(parser_wt_getprocessed,PARSER_WT_GETPROCESSED)(int *wti, int *wtp)
{
    int wtn = (*wti);
    parse->wt_getprocessed(wtn, wtp);
}

void FC_FUNC(parser_wt_setprocessed,PARSER_WT_SETPROCESSED)(int *wti, int *wtp)
{
    int wtn = (*wti);
    int wtpr = (*wtp);
    parse->wt_setprocessed(wtn, wtpr);
}


// Get and Set the sequence number for a whenthen.
void FC_FUNC(parser_wt_getseq,PARSER_WT_GETSEQ)(int *wti, int *wtseq)
{
    int wtn = (*wti);
    parse->wt_getseq(wtn, wtseq);
}

void FC_FUNC(parser_wt_setseq,PARSER_WT_SETSEQ)(int *wti, int *wtseq)
{
    int wtn = (*wti);
    int wtsequence = (*wtseq);
    parse->wt_setseq(wtn, wtsequence);
}

void FC_FUNC(parser_dictionary_add,PARSER_DICTIONARY_ADD)(char *name, double *value, int *pred, char *vdesc)
{
    bool c_pre_def = (*pred) ? true : false;
    parse->dictionary_add(name, *value, c_pre_def, vdesc);
}

void FC_FUNC(parser_dictionary_env_add,PARSER_DICTIONARY_ENV_ADD)(char *name, int *pred)
{
    bool c_pre_def = (*pred) ? true : false;
    parse->dictionary_env_add(name, c_pre_def);
}

void parser_get_num_include_files(int *num_include)
{
    *num_include = parse->NumIncludeFiles();
}

void parser_list_include_files()
{
    parse->ListIncludeFiles();
}

int parser_get_include_file_name_length(int *i)
{
    string fname;
    fname = parse->GetIncludeFile(*i);
    return fname.length();
}

void parser_get_include_file_name(char *file_name, int *file_name_len, int *i)
{
    string fname;
    fname = parse->GetIncludeFile(*i);
    *file_name_len = fname.length();
    strncpy(file_name,fname.c_str(),*file_name_len);
    file_name[*file_name_len] = '\0';
}


#ifdef __cplusplus
}
#endif




