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

#ifndef WORDHHINCLUDE
#define WORDHHINCLUDE

// ***************************************************************************
// ***************************************************************************
// This class holds each word from the line.
// ***************************************************************************
// ***************************************************************************

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <deque>
#include <stdint.h>

namespace PP
{
using std::string;
using std::stringstream;
using std::vector;
using std::map;
using std::deque;

enum WordType {WUNKNOWN, WSTRING, INTEGER, DOUBLE, EQUALS, OPERATOR,
               OPEN_PARENS, CLOSED_PARENS,
               OPEN_SQUARE_BRACKET, CLOSED_SQUARE_BRACKET,
               OPEN_BRACE, CLOSED_BRACE,
               COMMA, VARIABLE};

class Word
{

public:
    Word();
    Word(string s);
    Word(string s, int lnum, int file_lnum, string fname,
         deque<string> *lstr);
    Word(double d, int lnum, int file_lnum, string fname,
         deque<string> *lstr);
    Word(int ia, int lnum, int file_lnum, string fname,
         deque<string> *lstr);
    Word operator=(const Word &);
    Word(const Word &);
    ~Word();

    // Given a word, change its value.
    void set_value(double d);
    void set_value(string s);
    void set_value(const char *s);
    void set_value(bool b);

    // Print the type of word to a stringstream.
    void print_type(stringstream &ss);

    // Set the type of word.
    void set_type();

    // Get the operator type, arithmetic, relational, ...
    string get_op_type() {return op_type;}

    bool is_operator() { if (type == OPERATOR) return true; return false; }

    bool is_operator(int level) {
        if (type != OPERATOR) return false;
        if (level == op_level) return true;
        return false;
    }

    bool is_bool();

    bool is_string()   { if (type == WSTRING)  return true; return false; }

    bool is_integer()  { if (type == INTEGER) return true; return false; }

    bool is_number()   { if (type == INTEGER || type == DOUBLE)
                         return true; return false; }

    bool is_numvar()   { if (type == INTEGER || type == DOUBLE ||
                             type == VARIABLE) return true;
                         return false; }

    bool has_value()   { if (type == INTEGER || type == DOUBLE)
                         return true; return false; }

    bool is_variable() { if (type == VARIABLE) return true; return false; }

    bool is_comma()    { if (type == COMMA) return true; return false; }

    /*! \brief Get the word as a string. */
    string get_string()    { return wstr; }
    string get_print_string(bool enc_quotes);
    string get_stringp()   { processed=true; return wstr; }
    char get_single_char(stringstream &serr, int &ierr);

    /*! \brief Get the word as a float. */
    float get_float();

    /*! \brief Get the word as a double. */
    double get_double();
    double get_double(stringstream &serr, int &ierr);

    // Get the word as an int, without and with error processing.
    int get_int();
    int get_int(stringstream &serr, int &ierr);

    int64_t get_int64_t();
    int64_t get_int64_t(stringstream &serr, int &ierr);

    /*! \brief Get the word as a boolean. */
    bool get_bool(stringstream &serr, int &ierr);

    /*! \brief Templated get method for get_double, get_int, etc. */
    template< class T >
    T get_val( T &dummy );

    /*!
     * \brief Convert string s to the type of the first argument.  Function
     * returns the converted value as a reference and as the function
     * result. 
     *
     * These overloaded functions are the base for the templated accessor
     * functions "get_val()".
     *
     * \param rtti - Convert the string s into the type of rtti and return it.
     */
    int     convertFromString ( const int     &rtti, const string &s ) const; 
    int64_t convertFromString ( const int64_t &rtti, const string &s ) const; 
    string  convertFromString ( const string  &rtti, const string &s ) const;
    float   convertFromString ( const float   &rtti, const string &s ) const;
    double  convertFromString ( const double  &rtti, const string &s ) const;
    //bool    convertFromString ( const bool   &rtti, const string &s ) const;

    /*! \brief Negate the word or set a flag to negate it later. */
    void negate_value();

    // Handle errors.
    void fatal_error(stringstream &serr, int &ierr);
    void warning(stringstream &serr, int &ierr);

    // Miscellaneous functions.
    void handle_quotes(stringstream &serr, int &ierr);
    void erase_char(int ic);

    // Accessor methods.
    void   set_word(string str)    { wstr = str; }
    void   set_processed(bool p)   { processed = p; }
    bool   get_processed()         { return processed; }
    int    get_line_number()       { return line_number; }
    int    get_file_line_number()  { return file_line_number; }
    string get_filename()          { return filename; }
    deque<string> *get_lines()     { return lines; }
    void   set_filename(string fn) { filename = fn; }
    int    get_multiplicity()      { return multiplicity; }
    void   set_multiplicity(int m) { multiplicity = m; }

private:
    void init();
    bool check_before_e(string s, int i1, int i2);
    bool check_after_e(string s, int i1, int i2);

    // This is the basic storage for the word.
    string wstr;

    // The type of word, like operator, string, variable, etc.
    WordType type;

    // Flag for testing whether this word was processed or not.
    bool processed;

    // Flag to negate a variable.
    bool negate;

    // This word is repeated multiplicity times.
    int multiplicity;

    // If the word is an operator, then this is its level, i.e. "**" has
    // the highest level, then "*","/", etc.
    int op_level;

    // The operator type, arithmetic, relational, ...
    string op_type;

    // Convert the input string to lower case.
    void string_to_lower( string &s ) const;

    // This is needed for telling the user what line in the input
    // file or include file the error occurred on.
    //
    // line_number  The line_number corresponding to this command, this is
    //              an index into lines and starts from 1, not 0.
    // lines        Pointer to the deque of original lines. This contains all
    //              the lines from the input file and any include files.
    // file_line_number  The line number in the input file or include file.
    // filename          The name of the input file or include file.
    //
    // file_line_number and filename are needed to that the user can open
    // the file and go to the line in error.
    int line_number, file_line_number;
    string filename;
    deque<string> *lines;
    
};


} // end of PP namespace

#endif
