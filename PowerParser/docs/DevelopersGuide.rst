=================
Developer's Guide
=================

The *PowerParser* package is designed to handle the input file parsing needs for simulation
codes. The interface is in C++ with a parallel and serial version of the library. Operations
include commands such as simple key,value input::

   dx = 1.0

as well as array operations, mathematical operations, and restart support.

Each input line is called a command in the *PowerParser* nomenclature and
the left hand side of the expression is a key and the right hand side is
the value.

---------------------
Parser Initialization
---------------------

The parser can be initialized with one of the two following commands::

   #include "PowerParser.hh"
   using namespace PP;

   int main(int argc, char *argv[])
   {
      .....

      PowerParser *parse = new PowerParser();
         or
      PowerParser *parse = new PowerParser("file.in");

This operation will initialize the parsing object and either initialize MPI or use an already
initialized MPI context. The processor communication layer is automatically handled, but a
user can retrieve some MPI settings such as number of processors or rank with::

   int mype = parse->comm->getProcRank();
   int npes = parse->comm->getNumProcs();

Of course, a developer can get these settings by querying MPI directly.

Initiating the parsing of an input file can be done with::

   parse->parse_file("parsetest.in");

   parse->compile_buffer();

These calls perform the following operations

* Reads the input file into a deque of lines
* Breaks the lines up into words
* Compiles the input buffer

The user input buffer contains execution line arguments, the input files, and any
files included by them. Compiling the user input buffer handles loops, subroutines,
if statements, variables, etc. The end result is a final buffer that can be queried
by *get* commands.

There are several functions to echo out the input, the final buffer after compiling and
intermediate variables and functions. Some of the commands and their purpose are:: 

   parse->echo_input_start();

TODO: We need a lot more detail on the outputs.

TODO: Need description of debugging checks

----------------
Using the Parser
----------------

^^^^^^^^^^^^^^^
Scalar input
^^^^^^^^^^^^^^^

The retrieval of simple key-value pairs is straight-forward, such as the following::

   int intvalue = -1;
   parse->get_int("int_input", &intvalue);

This code looks for a simple key value pair in the parse buffer such as the following::

   int_input = 8

The get_int call returns the intvalue variable with it set to the value 8. There are
a couple of things to note. You should set the value of the variable before the get
call, because it is set only if there is a key of "int_input" in the input
file. It would end up uninitialized if the key is not present and could cause
indeterminate behavior. Also, the key string is case sensitive and so "int_input"
and "Int_Input" are different keys. 

The complete list of scalar get functions are::

   //Character array versions
      parse->get_bool    (const char *cname, bool    *cvalue)
      parse->get_bool_int(const char *cname, int     *cvalue)
      parse->get_int     (const char *cname, int     *cvalue)
      parse->get_int64_t (const char *cname, int64_t *cvalue)
      parse->get_real    (const char *cname, double  *cvalue)
      parse->get_char    (const char *cname, double  *cvalue)

   // String versions
      parse->get_bool    (string &cname,     bool    *cvalue)
      parse->get_bool_int(string &cname,     int     *cvalue)
      parse->get_int     (string &cname,     int     *cvalue)
      parse->get_int64_t (string &cname,     int64_t *cvalue)
      parse->get_char    (string &cname,     double  *cvalue)

These functions also take two possible optional arguments::

   const vector<int> &size -- set to {i}, {i,j}, {i,j,k}, etc
                              for arrays and multi-dimensional arrays
                              default = NULL, means scalar input

   bool skip = false -- has to do with skipping assignment to help
                        handle restarts 

The get_char functions also have an optional variable called single_char that
can be true or false.

^^^^^^^^^^^^^^^
Vector input
^^^^^^^^^^^^^^^

For an array input, the code would like this::

   vector<int> size = {6};
   double doublearray[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
   parse->get_real("array1d", doublearray, size);

or for a 2D array::

   vector<int> size = {3,2};
   double **doublearray2d = (double **)genmatrix(size[1], size[0], sizeof(double));
   for (int j = 0; j < size[1]; j++){
      for (int i = 0; i < size[0]; i++){
         doublearray2d[j][i] == -1.0;
      }
   }
   parse->get_real("array2d", &doublearray2d[0][0], size);

The genmatrix routine is a special allocator from the genmalloc package that
allocates a contiguous block of data for multi-dimensional arrays
and then assigns the pointers to the correct places in the
block of data to work correctly with multiple indices.

Not all the time is the size of the input know before-hand. For the 1D array
above, it often is better to query the size before allocating and reading
the array::

   vector<int> size;
   parse->get_size("array1d", size)
   double *doublearray = (double *)malloc(size * sizeof(double));
   for (int i = 0; i < size[0]; i++){
      doublearray[i] = -1.0;
   }
   parse->get_real("array1d", doublearray, size);

Up to 4D arrays are currently supported in *PowerParser*.

^^^^^^^^^^^^^^^^^^^
Key in input
^^^^^^^^^^^^^^^^^^^

Sometimes it is useful to query the parser to find out if a particular key appears
anywhere in the user input. This is done with the cmd_in_input function. For example,
suppose we want to know if the "special_variable" key is in the user input. We would make
the following call::

   string cname("special_variable");
   bool in_input = false;
   bool in_whenthen = false;
   parse->cmd_in_input(cname, in_input, in_whenthen);

The logical variable "in_input" will be true if "special_variable" is found anywhere in the
main input, not including the when...then blocks. The logical variable "in_whenthen" will be
true if "special_variable" is fond in any of the when...then blocks.

The search for the key is done on the final buffer and this makes it certain that it will 
exist if a get command for the key is done.

----------------------------
Parser Checks
----------------------------

After all the input routines have been called, the developer should check
that all the user commands in the user input file were processed. If any
command, or part of a command, was not processed, then a fatal error is produced. This
check is done with the following call::

   logical :: good
   parse->check_processed(good)

If a fatal error is generated, the code will end.

The assumption is being made that if a user command or part of a command is not processed
then it is most likely a user error, the user mistyped something, the user does not understand
the input, etc. In any case, the user is expected to fix the input file before it can be
successfully run.

This puts a burden on the developer. When the user turns on an option, it is expected that
the commands associated with that option can remain in the input file even though the
option is off and the commands are not really needed. Thus the developer needs to make
certain that all get commands associated with the option are processed even when the option
is turned off by the user.

Consider how the example of the size query should be changed if the package is turned off
and yet we still want to process the get call. We can use the optional skip argument to
process the input, but not set the variable, as follows::

   vector<int> size;
   parse->get_size("array1d", size)

   if (package_is_off) {
      parse->get_real("array1d", doublearray, size, true);
   } else {
      double *doublearray = (double *)malloc(size * sizeof(double));
      for (int i = 0; i < size[0]; i++){
         doublearray[i] = -1.0;
      }
      parse->get_real("array1d", doublearray, size, false);
   }

A better way to handle this is to use the cmd_set_processed routine to set the
processed flag for the command and its arguments to true. For example::

   vector<int> size;
   parse->get_size("array1d", size)

   if (package_is_off) {
      parse->cmd_set_processed("array1d", true);
   } else {
      double *doublearray = (double *)malloc(size * sizeof(double));
      for (int i = 0; i < size[0]; i++){
         doublearray[i] = -1.0;
      }
      parse->get_real("array1d", doublearray, size, false);
   }

The cmd_set_processed call takes two arguments; the first is the name of the command
and the second is the setting for the processed flag, either true or false.

----------------------
Parser Error Handling
----------------------

*PowerParser* does an exceptional job at handling and reporting errors. When an error
occurs, the parser reports the line number in the user input where the error occurred,
echos the line, reports what file the line is found in, and gives a detailed description
of the error. At this time, all parser errors are fatal errors, and *PowerParser* will
abort the run.

Most of the time, the errors are accumulated and are reported at the end of the
parsing process. This is good in that the user can correct more than one error at a
time. It is bad in that the errors can "snowball", in that one error can potentially
generate many spurious errors. In such cases the user only needs to fix the first
error to get rid of all the errors.

The file name is reported for an error because of the include capability of the parser,
the user can include external input files at any place in the main input file. The line
number is always the line number for the included file, or the main input file. This
makes it easy for the user to find the file and line where the error occurred.

When there are continuation lines in the input file, *PowerParser* combines those into
one line before processing. But the error reporting is always done for the continuation
line and not for the combined line. Consider the following line::

   $sp2 = true
   mult_logical_array(1) = 3*false &
                           2*$sp2, &
                           .truezz. 

There is an error in the last continuation line, true is misspelled. The parser reports
the following error::

   Fatal errors have been encountered while parsing the user input file.
   Note that often fixing the first error will also fix the other errors.

   *** FATAL ERROR in line 175:
                               .truezz.
   in file: test.in
   Values on this line should be true or false (or .true. or .false.)
      (any case is fine, for example true, True, TrUe are all ok)
   Instead found value: .truezz.

---------------
When...then
---------------

The first step in implementing the when...then code is to get the number of
when...then commands, wtnum, in the example below::

   int wtnum = 0;
   parse->whenthen_num_cpp(wtnum);

For every cycle in the simulation code, the following steps need to be done. First
setup two arrays called code_varnames and code_values. For example::

   for (int i=1; i<=wtnum; i++) {
      parse->whenthen_check_cpp(i, code_varnames, code_values, wt_check);
      if (wt_check) {
         parse->cpp("shortmodcyc", shortmodcyc, true);
      }
   }

The code_varnames and code_values are strings. The whenthen_check routine checks
to see if the condition has been satisfied. The first argument to this routine
is the when...then sequence number (starting from 1) as determined by the order
of the when...then commands in the user input file. The second argument is the
name of the simulation code internal variable that is to be checked. The third
argument is the value of the condition variable that is the value of the internal
code variable. The fourth argument is the true or false result output from this
routine, called wt_check in this example. If wt_check is returned as true, then
the condition was satisfied, otherwise it was not satisfied.

Note that the when...then class contains a processed flag, defaulted to false,
which is set to true the first time the condition is satisfied. Thus, after the
condition has been satisfied, subsequent calls to whenthen_check for this particular
when...then command will always return false.

The related whenever command does not set the processed flag.

When the condition is first satisfied, wt_check is true. Then the various calls
can be made. In the above example, the internal simulation variable shortmodcyc
can be changed by the user.

After cycling through all the when...then commands, the following call must be
made::

   parse->whenthen_reset();

This resets the pointer to the final commands buffer back to the main buffer.


