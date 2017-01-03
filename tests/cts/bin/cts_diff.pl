eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

#.............................................
#...add directory of executable to use path...
#.............................................
#use Data::Dumper;

# NOTE: You can use FindBin in the script. The modules will automatically have access to $FindBin::Bin.
use FindBin qw($RealBin);
use lib ("$RealBin", "$RealBin/lib", "$RealBin/../lib");

#.................
#...use modules...
#.................
use old_utils qw ( print_perl_obj );
use cts_diff_util qw ( parse_args );
#..............................
#...must have at least 1 arg...
#..............................
if ( $#ARGV < 0 )
  {
    print <<'EOF';
#............................................................................
#...Name
#...====
#... cts_diff.pl
#...   Find differences between data files.
#...   Each file is split into a set of dataset arrays (see "-ft" option
#...   below for definitions of file types).  The respective arrays
#...   are diff'd.  These arrays are Y-coordinate arrays.
#...
#...   Interpolation:
#...     If X-coordinate array(s) are given, interpolation is done in
#...     order to diff the values on the same X-coordinate.  The Y-values
#...     of the base dataset are interpolated onto the X-values of the new
#...     dataset.
#...         (X_base, Y_base) diff with (X_new, Y_new)
#...               (X_base, Y_base) -> interpolation -> (X_new, ~Y_base)
#...               (X_new, ~Y_base) diff with (X_new, Y_new)
#...     X values must be increasing (but not necessarily strictly)
#...     for interpolation to work.  The algorithm for dealing with out
#...     of order X values is subject to change with different version of
#...     cts_diff.pl (but invariant when using the same version of the code).
#...     Suggestion: Have X values be strictly increasing.
#...
#...     Interpolation may be turned off with the "-no_intp" flag (see below).
#...
#...Usage
#...=====
#... cts_diff.pl <files>
#...             [-a <absolute tolerance>[,<comma sep list of variables>]]
#...             [-arg <argument file>]
#...             [-ds <comma separated list of datasets to compare>]
#...             [-fsets <n - number of file sets>]
#...             [-ft <a|ares|cts|keyword|oxy|plot_output|pop|table|token|tracer|link|xy>]
#...             [-no_intp]
#...             [-no_plots]
#...             [-o_<file type> <output data file>]
#...             [-or]
#...             [-pft]
#...             [-plot_orig]
#...             [-presult]
#...             [-r <relative tolerance>[,<comma sep list of variables>]]
#...             [-status]
#...             [-v <number>]
#...
#...   <files>
#...      When given more than 1 data files, the first data file will
#...      be used as the base values from which all differences are taken.
#...      This can be used to detect drift.  Differences between the two
#...      latest files might be within tolerances when differences between
#...      the first and the last file would show differences.
#...         a1 a2 a3 --> (a1 vs. a2 vs. a3)
#...      See -fsets which changes what files are compared.
#...
#...   -a <abs tol>[,<variable list>]
#...   -r <rel tol>[,<variable list>]
#...      (The default is no tolerance.)
#...      Values are said to be "the same" if the absolute/relative
#...      difference is less than or equal to the absolute/relative
#...      tolernace.
#...
#...      Both -a and -r flags can be specified at once.
#...      If both -a and -r flags are specified,
#...      a value is "the same" if it is within (less than or equal to)
#...      _both_ the absolute _and_ relative tolerances.
#...      See the "-or" flag below.
#...
#...        Absolute Difference{[0,)} = abs(value1-value2)
#...        Relative Difference{[0,)} = abs(value1-value2)/
#...                                    abs(value1 if non 0 or value2)
#...
#...      These values will be used for all variables except those
#...      specified in a variable list as well:
#...          -a .005 -a .002,foo
#...
#...   -arg <argument file>
#...      (By default, only a file named "cts_diff.arg" is read).
#...      The arguments can also be read from a file.  Blank and comment
#...      lines (starting with "#") are ignored.  Otherwise, the contents
#...      would look as if they were typed on the command line.  This is
#...      useful if the arguments get long or are difficult to remember
#...      why they were chosen:
#...          | # default tolerances
#...          | -a .01 -r .05
#...          | # tolerances for specific datasets
#...          | -a .001,foo -r .0001,bar
#...
#...   -ds <comma separated list of datasets to compare>
#...      (By default, all datasets are compared).
#...      If specified, only datasets in the list will be compared.
#...      If datasets have whitespace, be sure to place the whole list
#...      in quotes.
#...
#...   -fsets <n - number of file sets>
#...      (The default is the number of files <num files>)
#...      Multiple file sets can be compared.  Each corresponding file in
#...      each file set is compared.  So, if you specify the following
#...      on the command line:
#...        -fsets 3 a1 b1 a2 b2 a3 b3
#...      The following sets of file will be compared:
#...        (a1 vs. a2 vs. a3) and (b1 vs. b2 vs. b3)
#...      With lots of files/data, it is probably best not to do any
#...      plotting (-no_plots) since the data file for plotting will be
#...      very large.
#...
#...   -ft <a|ares|cts|keyword|oxy|pop|table|token|tracer|link|xy>
#...      (By default, the type is determined automatically (by file name
#...       or a sampling of the source lines...if possible).)
#...      This specifies the data format of the file(s).
#...      o a: Reads in files from a particular code project
#...        Various variables are read in for that particular file type.
#...        For lack of better name, just call it a.
#...      o ares: Reads in the results from a suite of tests from a particular
#...        code project:
#...        | <Test Name> <PASSED|DIFF|FAILED>  P-<#>, D-<#>, F-<#>
#...        Each test name is its own dataset with an array of values
#...        (PASSED|DIFF|FAILED, P-#, D-#, F-#).
#...      o cts: Can read a the file produced when using the "-o" option.
#...        | # cts
#...        | # Dataset Name: <dataset name>
#...        | # Coord Name <X if it exists>: <coordinate name>
#...        | # Coord Name <Y - must exist>: <coordinate name>
#...        | <value 1 for coord X if it exists> <value 1 for coord Y>
#...        | <value 2 for coord X if it exists> <value 2 for coord Y>
#...        | ...repeat starting from Dataset Name for other datasets...
#...      o keyword: (cannot be autodetected)
#...        | <keyword> = <value>  (whitespace ignored - all on 1 line)
#...        | (all other lines ignored)
#...        dataset    = keyword
#...        X variable = none
#...        Y variable = keyword
#...      o oxy: Reads in files from a particular code project
#...        Various variables are read in for that particular file type.
#...        Previously, a post processor was used on the data to generate
#...        file type "xy".  Now, the original output file is parsed.
#...      o plot_output: plot_output.data file created by plot_output.pl
#...        You can give special "-ds" flags:
#...        -ds track
#...            skip the fields that will not track from run to run
#...            (like performance numbers)
#...        # [<num>] <dataset name>
#...        cycle  time  <col 1> <col 2> ...and so on
#...      o pop: This is a crude interpreter of the output from running
#...        this graphics package.
#...        Various keywords are found/used to create the data.
#...        dataset    = combo of card title/number and time ("* time").
#...        X variable = Line in output file: "* x var:"
#...        Y variable = Line in output file: "* y var:"
#...      o table: (lines are in columns with first line is header and
#...                following lines has columns of numbers)
#...        | <dataset name 1>  <dataset name 2> ...
#...        | <val ds 1>        <val ds 2>       ...
#...        | <val ds 1>        <val ds 2>       ...
#...        | ...               ...              ...
#...        dataset    = dataset name
#...        X variable = none
#...        Y variable = dataset name
#...      o tracer: (-tracer file ending and CSV file with particle/time fields)
#...        |  particle,      time, <dataset name 1>  <dataset name 2> ...
#...        | <particle #>, <time>, <val ds 1>        <val ds 2>       ...
#...        | <particle #>, <time>, <val ds 1>        <val ds 2>       ...
#...        |     ...         ...      ...              ...
#...        dataset    = p_<particle number>_<dataset name n>
#...        X variable = time
#...        Y variable = dataset name
#...      o link: (.lnk, .lnk.NNNNN file ending)
#...        Each whitespace separated token or Enn.n format real number
#...        is compared with the corresponding token in another file.
#...        dataset    = dataset name
#...        X variable = none
#...        Y variable = "token"
#...      o xy: (.std, .xy file endings, first line starts with # followed
#...             by 2 columns of values)
#...        | # <dataset name>
#...        | <x val> <y val>
#...        | <x val> <y val>
#...        | ...
#...        | # <dataset name>
#...        | <x val> <y val>
#...        | <x val> <y val>
#...        | ...[repeat for other datasets]
#...        dataset    = dataset name
#...        X variable = "X"
#...        Y variable = dataset name
#...      o token: (any other type)
#...        Each whitespace separated token is compared with the
#...        corresponding token in another file.
#...        dataset    = dataset name
#...        X variable = none
#...        Y variable = "token"
#...
#...      Additional Format Notes:
#...      o Values > 8e99 (GSKIP) or strings are tested for string match.
#...        If they are different, a difference is reported (but not used
#...        in the numerical analysis like mean, min, max, ...).
#...      o Blank lines are ignored.
#...
#...   -no_intp
#...      (By default, interpolation is done)
#...      Turn off interpolation.  Each dataset will be diff'd as-is.
#...      Interpolation does take extra time - so if you know that the
#...      files will have identical X values, turning this off will
#...      increase performance.
#...   -no_plots
#...      (By default, plots are created via gnuplot.)
#...      Do not print any plots.
#...
#...   -o_<file type> <output data file>
#...      Write the data of the last file in the list of files to the
#...      specified output data file.  If only 1 file is given, that
#...      file will be written.
#...      Currently, the following ar supported output types:
#..         -o_cts : cts format
#...        -o_xy  : xy format
#...
#...      This is useful if the original data file is large and you
#...      only want to store the data to be diffed.
#...      NOTE: Datesets that are skipped upon read will not be printed
#...            upon write.
#...
#...   -or
#...      (The default is "and").
#...      This flag changes the default behavior when given both absolute
#...      and relative tolerances. With this flag,
#...      a value is "the same" if it is within (less than or equal to)
#...      _either_ the absolute _or_ relative tolerances.
#...
#...   -pft
#...      Print the file type of each file after it is read.
#...      Used to verify if this script detects the correct type of file.
#...
#...   -plot_orig
#...      Do not compute differences - just plot the original data.
#...
#...   -presult
#...      Print the result at the end of each diff.
#...      Used by a project for greping.
#...        PASSED: No diffs.
#...        DIFF:   Diffs.
#...        FAILED: Failure (eg missing file)
#...
#...   -nostatus
#...      Return a 0 regardless of if diffs are found or not.
#...      Useful since CTS stops for non-0 return status in .test files
#...      (if you would like to continue doing diffs).
#...
#...   -v <number>
#...      Verbosity.  Larger number -> more info printed
#...      0: No info printed - just return status.
#...      1: (default) Final results from diffs of all files.
#...      2: Results from each file.
#...      3: Results from each dataset in each file.
#...      4: Every difference.
#...         At this verbosity, different Base and New values are printed
#...         along with the common X value interpolated to.
#...             (X, y_base) vs. (X, y_new)
#...
#...Return Value
#...============
#... If the -nostatus flag is not set:
#...   0: No differences found
#...   1: Otherwise
#...
#...Examples
#...========
#... 1) cts_diff.pl file1 file2 file3 -v 2 -no_plots -ds foo,bar
#...    Diff between 3 files.
#...    Diff summary from each file is printed.
#...    No plots are produced.
#...    Only the datasets foo and bar will be compared.
#...
#... 2) cts_diff.pl file1 file2 -a .4 -r .1 -or -o dat_file
#...    Diff two files.
#...    Data points are considered the "same" if
#...     absolute value of difference is <= .4 (-a .4) or (-or)
#...     the absolute value of the relative difference is <= .1 (-r .1).
#...    Plots will be produced.
#...    Create an output data file of the data in file2.
#...    Default verbosity (summary of diff).
#...
#... 3) cts_diff.pl gold_2002/* gold_2003/* gold_2004/* -fsets 3
#...    Diff three sets of gold standard files using gold_2002 as the base
#...    set.
#............................................................................
EOF
  exit;
  }
#............................................................................
#...Program Flow
#...============
#...1) parse command line
#...2) read in first file (base)
#...3) process each file
#...3.1) read file
#...3.2) for each dataset
#...3.2.1) diff/stat/print coordinates (X and Y data)
#...3.3) print file stats
#...4) print total stats
#...5) run gnuplot
#............................................................................
#..................
#...my variables...
#..................
my(
   @array_intp, # interpolated array
   $array_intp_ref, # reference to an array
   %cmd, # command line values: a, r, or, quiet, l, files[]
   $compare_group, # what compare group you are on
   %data, # data from reading
   %data_base, # data from first file
   $data_ref, # reference to a data struct
   %diff, # diff between 2 arrays
   %ds_diffs, # which datasets have diffs per file compare
   $ds_name, # dataset name
   %ds_names, # hash of dataset names
   $ds_names_ref, # ref to hash of ds names
   $headers_printed, # a header has been printed for a level
   $i, # loop var
   $ierr, # error return
   $is_diff, # if different
   $is_diff_local, # if this particular dataset is different
   @files, # list of files in compare group
   $file_name, # file name
   $file_base, # shortened file name for base file
   $file_new, # shortened file name for new file
   $num_compare_groups, # number of compare groups
   $coord, # coordinate
   $coord_name, # coordinate name
   %coords, # coordinates
   $print_source, # source name label to print
   $print_var, # variable name label to print
   $print_header, # print a header
   $print_title, # source name label to print
   $printed_diffs, # if printed diffs since last header
   %result_diff, # for -presult, stores additional data for printing later.
   %stat, # a stat on a diff
   %stat_file, # stat results for file
   %stat_cmp_grp, # stat results for everything
   %stat_total, # total stat results
   %gnuplot_info, # info used for plotting
   $num_diffs, # number of diffs
  );
#.................
#...global vars...
#.................
my(
   $GDATA,   # global name for data type
   $GORG,    # global name for data type
   $GREL,    # global name for data type
   $GABS,    # global name for data type
   $GDIFF,   # global name for diff of dataset
   $GCOORDY, # global name for coordinate
   $GCOORDX, # global name for coordinate
   $GNAME,   # global name for name hash
   );
$GORG    = $cts_diff_util::GORG;
$GREL    = $cts_diff_util::GREL;
$GABS    = $cts_diff_util::GABS;
$GDATA   = $cts_diff_util::GDATA;
$GDIFF   = $cts_diff_util::GDIFF;
$GNAME   = $cts_diff_util::GNAME;
$GCOORDX = $cts_diff_util::GCOORDX;
$GCOORDY = $cts_diff_util::GCOORDY;
#........................
#...parse command line...
#........................
$ierr = cts_diff_util::parse_args( \@ARGV, \%cmd );
if( $ierr != 0 )
  {
    cts_diff_util::print_error( "FAILED: Invalid Command Line",
                                "Run [$0] without args for help.",
                                $ierr );
    exit( $ierr );
  }
if( defined( $cmd{help} ) )
  {
    $ierr = 1;
    cts_diff_util::print_error( "Run command [$0] without any arguments ".
                                "to get help",
                                $ierr );
    exit( $ierr );
  }
if( ! defined( $cmd{files} ) )
  {
    $ierr = 1;
    cts_diff_util::print_error( "Must supply at least 1 file",
                                "Run command [$0] without any arguments ".
                                "to get help",
                                $ierr );
    exit( $ierr );
  }
#.......................
#...default verbosity...
#.......................
if( !defined( $cmd{v} ) )
  {
    $cmd{v} = 1;
  }
#...........................
#...default file set size...
#...........................
if( !defined( $cmd{fsets} ) || $cmd{fsets} <= 0 )
  {
    $cmd{fsets} = $#{$cmd{files}}+1;
  }
$headers_printed = 0;
$printed_diffs = 0;
$num_compare_groups = ($#{$cmd{files}} + 1) / $cmd{fsets};
#................................
#...process each compare group...
#................................
for( $compare_group = 1; $compare_group <= $num_compare_groups;
     $compare_group++ )
  {
    undef( %stat_cmp_grp );
    undef( %result_diff );
    #......................
    #...create file list...
    #......................
    $index = $compare_group-1;
    undef( @files );
    for( $i = 0; $i < $cmd{fsets}; $i++ )
      {
        push( @files, $cmd{files}[$index] );
        $index = $index + $num_compare_groups;
      }
    #........................
    #...read in first file...
    #........................
    undef( %data_base );
    #printf "debug %20s %s", "before read", `date "+\%S \%N"`;
    $ierr = cts_diff_util::read_file( $files[0], \%cmd, \%data_base );
    #printf "debug %20s %s", "after read", `date "+\%S \%N"`;
    if( $ierr != 0 )
      {
        cts_diff_util::print_error( "FAILED: Error in reading data file(s)",
                                    $ierr );
        exit( $ierr );
      }
    $is_diff = 0;
    undef( %ds_diffs );
    ($file_base = $files[0]) =~ s/^.*\///;
    #.......................
    #...process each file...
    #.......................
    for( $i = 1; $i <= $#files; $i++ )
      {
        #...............
        #...read file...
        #...............
        $file_name = $files[$i];
        ($file_new = $file_name) =~ s/^.*\///;
        undef( %data );
        #printf "debug %20s %s", "before read", `date "+\%S \%N"`;
        $ierr = cts_diff_util::read_file( $file_name, \%cmd, \%data );
        #printf "debug %20s %s", "after read", `date "+\%S \%N"`;
        if( $ierr != 0 )
          {
            cts_diff_util::print_error( "FAILED: Error in reading data file(s)",
                                        $ierr );
            exit( $ierr );
          }
        #......................
        #...get all ds_names...
        #......................
        undef( %ds_names );
        foreach $ds_name ( keys %{$data_base{$GDATA}} )
          {
            $ds_names{$ds_name} = "";
          }
        foreach $ds_name ( keys %{$data{$GDATA}} )
          {
            $ds_names{$ds_name} = "";
          }
        undef( %stat_file );
        #..............................................
        #...diff each dataset with base file dataset...
        #..............................................
        foreach $ds_name ( sort keys %ds_names )
          {
            $is_diff_local = 0;
            if( !defined( $cmd{plot_orig} ) )
              {
                &do_diff( \%cmd,
                          $file_base,
                          $file_new,
                          \%data,
                          \%data_base,
                          $ds_name,
                          \%stat,
                          \%stat_file,
                          \%ds_diffs,
                          \$printed_diffs,
                          \$headers_printed,
                          \$is_diff_local,
                          \$is_diff );
              }
            #.............................................
            #...print_gnuplot_data if diff or plot_orig...
            #.............................................
            if( ! defined( $cmd{no_plots} ) &&
                ( $is_diff_local != 0 || defined( $cmd{plot_orig} ) ) )
              {
                if( defined( $data{$GDATA}{$ds_name}{$GCOORDY} ) )
                  {
                    $print_var = "cmp=$compare_group:$ds_name";
                    $print_source = sprintf( "cmp=%3d:%3d:%s",
                                             $compare_group, $i, $file_new );
                    $print_title = sprintf( "%s vs. %s -> {ds}%s",
                                            $file_base, $file_new,
                                            $ds_name
                                          );
                    #printf "debug %20s %s", "before print_gnuplot", `date "+\%S \%N"`;
                    cts_diff_util::print_gnuplot_data( \%{$data{$GDATA}{$ds_name}},
                                                       $print_var,
                                                       $print_source,
                                                       $print_title,
                                                       \%gnuplot_info );
                    #printf "debug %20s %s", "after print_gnuplot", `date "+\%S \%N"`;
                  }
              }
          }
        #....................................................
        #...DONE: diff each dataset with base file dataset...
        #....................................................
        #.....................................................
        #...merge file to cmp_group stats and print abs/rel...
        #.....................................................
        cts_diff_util::merge_stats( \%{$stat_file{$GABS}},
                                    \%{$stat_cmp_grp{$GABS}} );
        cts_diff_util::merge_stats( \%{$stat_file{$GREL}},
                                    \%{$stat_cmp_grp{$GREL}} );
        if( $cmd{v} >= 2 )
          {
            #...print header if printing diffs or if have not yet...
            if( $printed_diffs ||
                $headers_printed == 0 )
              {
                $print_header = 1;
              }
            else
              {
                $print_header = 0;
              }
            $print_title = sprintf( "  %s vs. %s", $file_base, $file_new );
            cts_diff_util::print_abs_rel_stats( $print_title, $print_header,
                                                \%{$stat_file{$GABS}},
                                                \%{$stat_file{$GREL}},
                                                \$headers_printed );
          }
        #...........................................................
        #...DONE: merge file to cmp_group stats and print abs/rel...
        #...........................................................
      }
    #.............................
    #...DONE: process each file...
    #.............................
    #..............................
    #...plot data_base if needed...
    #..............................
    if( ! defined( $cmd{no_plots} ) )
      {
        #...................................................................
        #...get list of datasets to plot - all if plot_orig or only diffs...
        #...................................................................
        if( defined( $cmd{plot_orig} ) )
          {
            $ds_names_ref = \%{$data_base{$GDATA}};
          }
        else
          {
            $ds_names_ref = \%ds_diffs;
          }
        foreach $ds_name ( sort keys %$ds_names_ref )
          {
            if( defined( $data_base{$GDATA}{$ds_name}{$GCOORDY} ) )
              {
                $print_var = "cmp=$compare_group:$ds_name";
                $print_source = sprintf( "cmp=%3d:%3d:%s",
                                       $compare_group, 0, $file_base );
                $print_title = sprintf( "base %s -> {ds}%s",
                                        $file_base, $ds_name
                                      );
                #.......................................................
                #...create 0 diff data: keeps key and draws diff axis...
                #.......................................................
                if( ! defined( $cmd{plot_orig} ) )
                  {
                    #printf "debug %20s %s", "before base diff", `date "+\%S \%N"`;
                    cts_diff_util::create_diff( \@{$data_base{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                                                \@{$data_base{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                                                \%cmd,
                                                $data_base{$GDATA}{$ds_name}{$GCOORDY}{$GNAME},
                                                "", 0, undef,
                                                \%{$data_base{$GDATA}{$ds_name}{$GCOORDY}{$GDIFF}},
                                                \$num_diffs );
                    #printf "debug %20s %s", "after base diff", `date "+\%S \%N"`;
                  }
                #.............................................................
                #...DONE: create 0 diff data: keeps key and draws diff axis...
                #.............................................................
                #printf "debug %20s %s", "before print_gnuplot", `date "+\%S \%N"`;
                cts_diff_util::print_gnuplot_data( \%{$data_base{$GDATA}{$ds_name}},
                                                   $print_var, $print_source,
                                                   $print_title,
                                                   \%gnuplot_info );
                #printf "debug %20s %s", "after print_gnuplot", `date "+\%S \%N"`;
              }
          }
      }
    #....................................
    #...DONE: plot data_base if needed...
    #....................................
    #..............................................
    #...merge cmp_group to total and print stats...
    #..............................................
    cts_diff_util::merge_stats( \%{$stat_cmp_grp{$GABS}},
                                \%{$stat_total{$GABS}} );
    cts_diff_util::merge_stats( \%{$stat_cmp_grp{$GREL}},
                                \%{$stat_total{$GREL}} );
    if( $cmd{v} == 1 ||
        $cmd{v} == 2 && $#files > 1 ||
        $cmd{v} >= 3 )
      {
        if( $printed_diffs ||
            $headers_printed == 0 )
          {
            $print_header = 1;
            $printed_diffs = 0;
          }
        else
          {
            $print_header = 0;
          }
        $print_title = sprintf( " base %s", $file_base );
        cts_diff_util::print_abs_rel_stats( $print_title, $print_header,
                                            \%{$stat_cmp_grp{$GABS}},
                                            \%{$stat_cmp_grp{$GREL}},
                                            \$headers_printed );
      }
    #.............................
    #...print out -presult info...
    #.............................
    if( defined( $cmd{presult} ) )
      {
        foreach $coord_name ( sort keys %result_diff )
          {
            printf( "%s [%s]\n",
                    $result_diff{$coord_name}, $coord_name );
          }
      }
    #............................
    #...write output data file...
    #............................
    if( $compare_group == $num_compare_groups )
      {
        #...last file of compare group...
        if( $#files == 0 )
          {
            $data_ref = \%data_base;
          }
        else
          {
            $data_ref = \%data;
          }
        cts_diff_util::print_data_file( \%cmd, $data_ref );
      }
  }
#......................................
#...DONE: process each compare group...
#......................................
#.....................................................
#...print total stats (unless only 1 compare group)...
#.....................................................
if( $cmd{v} > 0 && $num_compare_groups > 1 )
  {
    if( $printed_diffs ||
        $headers_printed == 0 )
      {
        $print_header = 1;
        $printed_diffs = 0;
      }
    else
      {
        $print_header = 0;
      }
    $print_title = sprintf( "%s", "Total" );
    cts_diff_util::print_abs_rel_stats( $print_title, $print_header,
                                        \%{$stat_total{$GABS}},
                                        \%{$stat_total{$GREL}},
                                        \$headers_printed );
  }
#.................
#...run gnuplot...
#.................
if( ! defined( $cmd{no_plots} ) && %gnuplot_info )
  {
    if( $cmd{v} >= 1 )
      {
        print "Creating Gnuplot Plots...\n";
      }
    #printf "debug %20s %s", "before run_gnuplot", `date "+\%S \%N"`;
    #$gnuplot_info{orientation} = "portrait";
    cts_diff_util::run_gnuplot( \%gnuplot_info );
    #printf "debug %20s %s", "after run_gnuplot", `date "+\%S \%N"`;
  }
if( defined( $cmd{nostatus} ) )
  {
    exit( 0 );
  }
else
  {
    exit( $is_diff );
  }
#.............................................................................
#............................. Subroutines ...................................
#.............................................................................
sub do_diff
  {
    my(
       $cmd_ref,
       $file_base,
       $file_new,
       $data_ref,
       $data_base_ref,
       $ds_name,
       $stat_ref,
       $stat_file_ref,
       $ds_diffs_ref,
       $printed_diffs_ref,
       $headers_printed_ref,
       $is_diff_local_ref,
       $is_diff_ref,
      ) = @_;
    my(
       $array_intp_ref, # reference interpolated array
       $array_x_ref, # reference to orig array
       $coord_name, # name of the coordinate
       $ierr, # error return value
       @intp_array, # interpolated array
       $print_title, # the title to use for printing
       $print_diffs, # if printing the individual diffs
       $num_diffs, # the number of diffs after create_diff
      );
    $$is_diff_local_ref = 0;
    #..................................
    #...if printing every difference...
    #..................................
    if( $$cmd_ref{v} >= 4 )
      {
        $print_diffs = 1;
      }
    else
      {
        $print_diffs = 0;
      }
    #..............................................................
    #...if X coords exist for both current and base, interpolate...
    #...from (base.x, base.y) -> (new.x, base.y')               ...
    #..............................................................
    if( ! defined($$cmd_ref{no_intp} ) &&
        defined( $$data_base_ref{$GDATA}{$ds_name}{$GCOORDX} ) &&
        defined( $$data_ref{$GDATA}{$ds_name}{$GCOORDX} ) )
      {
        undef( @intp_array );
        $array_x_ref = \@{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}};
        $ierr = cts_diff_util::interpolate( \@{$$data_base_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                                            \@{$$data_base_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                                            $array_x_ref,
                                            \@intp_array );
        if( $ierr != 0 )
          {
            cts_diff_util::print_error( "Interpolation failed for",
                                        "File: [$file_base]",
                                        "Dataset: [$ds_name]",
                                        "Using uninterpolated values",
                                        0 );
            $array_x_ref = undef;
            $array_intp_ref = \@{$$data_base_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}};
          }
        else
          {
            $array_intp_ref = \@intp_array;
          }
      }
    else
      {
        $array_x_ref = undef;
        $array_intp_ref = \@{$$data_base_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}};
      }
    #..........................................................
    #...get coordinate name from data, data_base, or default...
    #..........................................................
    $coord_name = $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME};
    if( ! defined( $coord_name ) )
      {
        $coord_name = $$data_base_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME};
      }
    if( ! defined( $coord_name ) )
      {
        $coord_name = $coord;
      }
    #.............................................
    #...reset coord name incase was not defined...
    #.............................................
    $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME}      = $coord_name;
    $$data_base_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $coord_name;
    $print_title = sprintf( "%s vs. %s -> {ds}%s:{%s}%s",
                            $file_base, $file_new,
                            $ds_name, $GCOORDY, $coord_name
                          );
    #..........
    #...diff...
    #..........
    $num_diffs = 0;
    #printf "debug %20s %s", "before diff", `date "+\%S \%N"`;
    cts_diff_util::create_diff( $array_intp_ref,
                                \@{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                                $cmd_ref, $coord_name,
                                "$print_title", $print_diffs,
                                $array_x_ref,
                                \%{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GDIFF}},
                                \$num_diffs );
    #printf "debug %20s %s", "after diff", `date "+\%S \%N"`;
    if( $num_diffs > 0 )
      {
        $ds_diffs{$ds_name} = "";
        $$is_diff_local_ref = 1;
        $$is_diff_ref = 1;
        if( $print_diffs )
          {
            $$printed_diffs_ref = 1;
          }
        $result_diff{"$compare_group:$ds_name:$coord_name"} =
          "DIFF  ";
      }
    else
      {
        $result_diff{"$compare_group:$ds_name:$coord_name"} =
          "PASSED";
      }
    #printf "debug %20s %s", "before stats", `date "+\%S \%N"`;
    cts_diff_util::create_stats( \@{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GDIFF}{$GABS}},
                                 \%{$$stat_ref{$GABS}} ,\@{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GDIFF}{$GREL}});
    cts_diff_util::merge_stats( \%{$$stat_ref{$GABS}},
                                \%{$$stat_file_ref{$GABS}} );
    cts_diff_util::create_stats( \@{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GDIFF}{$GREL}},
                                 \%{$$stat_ref{$GREL}} ,\@{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GDIFF}{$GABS}});
    cts_diff_util::merge_stats( \%{$$stat_ref{$GREL}},
                                \%{$$stat_file_ref{$GREL}} );
    #printf "debug %20s %s", "after stats", `date "+\%S \%N"`;
    if( $$cmd_ref{v} >= 3 )
      {
        #...print header if printing diffs or if have not yet...
        if( $$printed_diffs_ref ||
            $headers_printed == 0 )
          {
            $print_header = 1;
            $$printed_diffs_ref = 0;
          }
        else
          {
            $print_header = 0;
          }
        cts_diff_util::print_abs_rel_stats( "   $print_title",
                                            $print_header,
                                            \%{$$stat_ref{$GABS}},
                                            \%{$$stat_ref{$GREL}},
                                            \$headers_printed );
      }
  }


=head1 COPYRIGHT AND LICENSE

 Copyright (2006). The Regents of the University of California. This material was
 produced under U.S. Government contract W-7405-ENG-36 for Los Alamos National
 Laboratory, which is operated by the University of California for the U.S. Department
 of Energy. The U.S. Government has rights to use, reproduce, and distribute this
 software.  NEITHER THE GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR
 IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is
 modified to produce derivative works, such modified software should be clearly marked,
 so as not to confuse it with the version available from LANL.

 Additionally, this program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later version.
 Accordingly, this program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE. See the GNU General Public License for more details.



=cut

