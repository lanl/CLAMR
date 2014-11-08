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
// This class holds information about a variable. It is mostly for use with
// the parser.
// ***************************************************************************
// ***************************************************************************

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <deque>

#include "Parser_utils.hh"
#include "Word.hh"
#include "Variable.hh"

namespace PP
{
using std:: string;
using std::cout;
using std::endl;
using std::stringstream;
using std::setprecision;
using std::vector;
using std::deque;

static int index_base = 1;


// ===========================================================================
// Default constructor.
// ===========================================================================
Variable::Variable()
{
    name = "__NO_NAME_GIVEN__";
    value.push_back("__NO_VALUE_GIVEN__");
    ndim = -1;
    lnum_ndim = 0;
    lnum_bounds = 0;
    pre_defined = false;
    description = "";
    temporary = false;
}

// ===========================================================================
// Constructor to reset index base
// ===========================================================================
Variable::Variable(int base)
{
    index_base = base;
}

// ===========================================================================
// Constructor given a string as input. This constructs a scalar variable.
// ===========================================================================
Variable::Variable(string nme, string v, bool pred, string tdes)
{
    name = nme;
    value.push_back(v);
    ndim = 0;
    lnum_ndim = 0;
    lnum_bounds = 0;
    pre_defined = pred;
    description = tdes;
    temporary = false;
}


// ===========================================================================
// Constructor for variables with no value.
// ===========================================================================
Variable::Variable(string nme)
{
    name = nme;
    value.push_back("__NO_VALUE_GIVEN__");
    ndim = -1;
    lnum_ndim = 0;
    lnum_bounds = 0;
    pre_defined = false;
    description = "";
    temporary = false;
}


// ===========================================================================
// Constructor given a vector of strings as input.
// ===========================================================================
Variable::Variable(string nme, vector<int> &istart, vector<string> &valvec,
                   int lnum, int file_lnum, string fname, deque<string> *lines,
                   stringstream &serr, int &ierr)
{
    name = nme;
    ndim = -1;
    lnum_ndim = 0;
    lnum_bounds = 0;
    pre_defined = false;
    description = "";
    temporary = false;
    set_var_value(istart, valvec, lnum, file_lnum, fname, lines, serr, ierr);
}


// ===========================================================================
// istart gives the starting location in the array for setting values.
// The istart indices start from 1 (fortran based).
//
// This function works for any dimension, 0,1,2,3,...
// ===========================================================================
void Variable::set_var_value(vector<int> &istart, vector<string> &valvec,
                             int lnum, int file_lnum, string fname,
                             deque<string> *lines, stringstream &serr, int &ierr)
{
    // Cannot redefine a pre-defined variable.
    if (pre_defined) {
        serr << endl;
        serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
        serr << "    " << (*lines)[lnum-1] << endl;
        serr << "in file: " << fname << endl;
        serr << "Cannot redefine a pre-defined variable." << endl;
        ierr = 2;
        return;
    }
    
    // Set the array dimension and make sure the user is not trying to
    // change it.
    int ndim_new = (int)istart.size();
    if (ndim >= 0) {
        if (ndim != ndim_new) {
            // Throw an error
            serr << endl;
            serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
            serr << "    " << (*lines)[lnum-1] << endl;
            serr << "in file: " << fname << endl;
            serr << "Cannot redefine the dimensionality of a variable." << endl;
            serr << "Original number of dimensions = " << ndim << endl;
            serr << "New number of dimensions = " << ndim_new << endl;
            if (lnum_ndim > 0) {
                serr << "Previously set in line " << lnum_ndim << ":" << endl;
                serr << "    " << (*lines)[lnum_ndim-1] << endl;
            }
            ierr = 2;
            return;
        }
    }
    else {
        ndim = ndim_new;
        lnum_ndim = lnum;
    }

    int bsize = (int)maxdim.size();

    if (ndim == 0 && bsize > 0) {
        serr << endl;
        serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
        serr << "    " << (*lines)[lnum-1] << endl;
        serr << "in file: " << fname << endl;
        serr << "Array boundaries not allowed for scalar variable." << endl;
        if (lnum_bounds > 0) {
            serr << "Array boundaries were set in line " << lnum_bounds << ":" << endl;
            serr << "    " << (*lines)[lnum_bounds-1] << endl;
        }
        ierr = 2;
        return;
    }

    if (ndim > 0) {
        if (ndim != bsize+1) {
            serr << endl;
            serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
            serr << "    " << (*lines)[lnum-1] << endl;
            serr << "in file: " << fname << endl;
            serr << "Number of dimensions = " << ndim << endl;
            serr << "Number of array boundaries + 1 = " << bsize+1 << endl;
            serr << "These should match but don't. " << endl;
            ierr = 2;
            return;
        }
    }

    // Find the 1d starting position given multiple array indices.
    Parser_utils putils(index_base);
    int i1 = putils.start_dex(istart, maxdim);

    // nvals    Number of values after the = sign.
    // Note that multiplicity is already handled, i.e. valvec has already
    // been expanded to include multiplicites.
    int nvals = (int)valvec.size();

    // Get memory that we need.
    if (i1+nvals > (int)value.size()) {
        value.resize(i1+nvals, "");
    }

    // Store the array values.
    for (int i=i1; i<i1+nvals; i++) {
        value[i] = valvec[i-i1];
    }
}


// ===========================================================================
// Increment (or decrement) a variable value by an integer amount.
// ===========================================================================
void Variable::bump_var(vector<int> &istart, int increment,
                        int lnum, int file_lnum, string fname,
                        deque<string> *lines, stringstream &serr, int &ierr)
{
    // Find the 1d starting position given multiple array indices.
    Parser_utils putils(index_base);
    int i1 = putils.start_dex(istart, maxdim);

    // We are incrementing an existing variable, so i1 should be valid.
    if (i1 >= (int)value.size()) {
        // Fatal Error
    }

    Word w1(value[i1], lnum, file_lnum, fname, lines);
    if (!w1.is_number()) {  } // FATAL ERROR 

    stringstream ss;
    if (w1.is_integer()) {
        int ia1 = w1.get_int(serr, ierr);
        int ia  = ia1 + increment;
        ss << ia;
    }
    else {
        double d1 = w1.get_double(serr, ierr);
        double d = d1 + increment;
        ss << setprecision(15) << d;
    }
    value[i1] = ss.str();
}



// ===========================================================================
// The problem with multi-dimensional variable arrays is that the user has
// to tell us the bounds on every dimension except the last one. This info
// is input in the bounds vector and stored.
//
// This function works for any dimension, 0,1,2,3,...
// ===========================================================================
void Variable::set_bounds(vector<int> &bounds, int lnum, int file_lnum,
                          string fname, deque<string> *lines,
                          stringstream &serr, int &ierr)

{
    // Cannot redefine a pre-defined variable.
    if (pre_defined) {
        serr << endl;
        serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
        serr << "    " << (*lines)[lnum-1] << endl;
        serr << "in file: " << fname << endl;
        serr << "Cannot redefine a pre-defined variable." << endl;
        ierr = 2;
        return;
    }
    
    // Set the array dimension and make sure the user is not trying to
    // change it.
    int ndim_new = (int)bounds.size() + 1;
    if (ndim >= 0) {
        if (ndim != ndim_new) {
            serr << endl;
            serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
            serr << "    " << (*lines)[lnum-1] << endl;
            serr << "in file: " << fname << endl;
            serr << "Cannot redefine the dimensionality of a variable "
                "(set_bounds)." << endl;
            serr << "Original number of dimensions = " << ndim << endl;
            serr << "New number of dimensions = " << ndim_new << endl;
            if (lnum_ndim > 0) {
                serr << "Previously set in line " << lnum_ndim << ":" << endl;
                serr << "    " << (*lines)[lnum_ndim-1] << endl;
            }
            ierr = 2;
            return;
        }
    }
    else {
        ndim = ndim_new;
        lnum_ndim = lnum; // Store line num for better err messages.
    }

    // Check to make sure the user is not reseting the bounds.
    if ((int)maxdim.size() > 0) {
        serr << endl;
        serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
        serr << "    " << (*lines)[lnum-1] << endl;
        serr << "in file: " << fname << endl;
        serr << "The bounds on this array has already been set," << endl;
        serr << "cannot reset them." << endl;
        if (lnum_bounds > 0) {
            serr << "Previously set in line " << lnum_bounds << ":" << endl;
            serr << "    " << (*lines)[lnum_bounds-1] << endl;
        }
        ierr = 2;
        return;
    }

    // Store the line num where bounds were set for better err messages.
    lnum_bounds = lnum;

    // Store the bounds.
    maxdim.clear();
    for (int i=0; i<(int)bounds.size(); i++) {
        maxdim.push_back(bounds[i]);
    }
}


// ===========================================================================
// Given indices, in adex, get the value of the variable.
// For example, suppose you want the value of
//    $var2d(3,5)
// The adex vector contains 2 numbers, 3 and 5, for the fortran indices.
// The start_dex function is used to get the 1d index into the value array.
//
// This function works for any dimension, 0,1,2,3,...
// ===========================================================================
string Variable::get_var_value(vector<int> &adex, string vname, int lnum,
                               int file_lnum, string fname, deque<string> *lines,
                               stringstream &serr, int &ierr)
{
    int adex_size = (int)adex.size();

    // Special case for scalar variables.
    if (ndim == 0 || adex_size == 0) return value[0];

    // The adex indices and bounds indices must match.
    if (adex_size - 1 != (int)maxdim.size()) {
        serr << endl;
        serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
        serr << "    " << (*lines)[lnum-1] << endl;
        serr << "in file: " << fname << endl;
        serr << "The dimensionality of variable " << vname << endl;
        serr << "does not match what was previously set." << endl;
        if (lnum_bounds > 0) {
            serr << "Previous dimensionality set in line " << lnum_bounds << ":" << endl;
            serr << "    " << (*lines)[lnum_bounds-1] << endl;
        }
        ierr = 2;
        return "";
    }
    
    // Indices cannot exceed max allowed.
    // Remember that adex if referenced from 1 (fortran index).
    for (int d=0; d<(int)maxdim.size(); d++) {
        if (adex[d] > maxdim[d]) {
            serr << endl;
            serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
            serr << "    " << (*lines)[lnum-1] << endl;
            serr << "in file: " << fname << endl;
            serr << "Variable name = " << vname << endl;
            serr << "The value for dimension " << d+1 << " = " << adex[d] << endl;
            serr << "This exceeds the max dimension of " << maxdim[d] << endl;
            if (lnum_bounds > 0) {
                serr << "The array bounds were set in line " << lnum_bounds << ":" << endl;
                serr << "    " << (*lines)[lnum_bounds-1] << endl;
            }
            ierr = 2;
        }
    }

    if (ierr == 2) return "";

    // Indices cannot be < 1.
    for (int d=0; d<(int)adex.size(); d++) {
        if (adex[d] < 1) {
            serr << endl;
            serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
            serr << "    " << (*lines)[lnum-1] << endl;
            serr << "in file: " << fname << endl;
            serr << "Variable name = " << vname << endl;
            serr << "Expected index greater than or equal to 1 " << endl;
            serr << "Instead, index =  " << adex[d] << endl;
            ierr = 2;
        }
    }

    if (ierr == 2) return "";


    // Using the indices in adex and the bounds for multi-d arrays, maxdim,
    // get the 1d index into the value array.
    Parser_utils putils(index_base);
    int i1 = putils.start_dex(adex, maxdim);

    // Check that the value array size has not been exceeded.
    if (i1 >= (int)value.size()) {
        serr << endl;
        serr << "*** FATAL ERROR in line " << file_lnum << ":" << endl;
        serr << "    " << (*lines)[lnum-1] << endl;
        serr << "in file: " << fname << endl;
        serr << "Variable name = " << vname << endl;
        serr << "Exceeded array bounds. Check to make sure you are not" << endl;
        serr << "requesting an array element you have not yet set." << endl;

        vector<int> maxdex((int)adex.size(), 0);
        get_indices((int)value.size()-1, maxdex);
        for (int d=0; d<(int)adex.size(); d++) {
            string s = "";
            if (adex[d] > maxdex[d]) s = "  ERROR, max exceeded";
            serr << "    Requested index = " << adex[d]
                 << "  Max index = " << maxdex[d] << s << endl;
        }

        ierr = 2;
        return "";
    }

    // Return the value.
    return value[i1];
}



// ===========================================================================
// Given the 1d index, icdex (starting from 0), find the corresponding
// multi dimensional fortran indices (each starting from 1).
//
// Example 1: Consider a 1d array
//     var1d(1) = 1 3 5 9 -4 -5 6
// Suppose icdex=3, corresponding to array value 9.
// This 1d case is very simple, all we do is add 1 to icdex to get a reference
// from 1, thus returning 4.
//
// Example 1: Consider a 2d array
//     $var2d(1,1) = 11. 21. 31.   12. 22. 32.   13. 23. 33.
// Where the max of the first dimension is 3. Suppose the user specifies
// icdex = 5, this corresponds to array value 32. The two indices returned
// would be 3,2 (referenced from 1).
//
// The adex vector contains the output indices, for example 2 this would be 3
// and 2.
//
// This function works for any dimension, 0,1,2,3,...
// ===========================================================================
void Variable::get_indices(int icdex, vector<int> &adex)
{
    // Nothing to do for scalar variables.
    if (ndim == 0) return;

    // Given icdex, get the indices.
    int nvalues = (int)value.size();
    Parser_utils putils(index_base);
    putils.reverse_dex(icdex, nvalues, adex, maxdim);
}




} // End of the PP namespace
