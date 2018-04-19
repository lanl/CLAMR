#.........................................................
#...Various utilities for doing diffs on files        ...
#.........................................................
package cts_diff_util;

use     POSIX qw( strtod );
#use     diagnostics;
use     warnings;
#use     DumpStack;
#use     Data::Dumper;
use     Carp;
use     vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );
use     Exporter;
use     old_utils qw ( print_error );

$VERSION   = 1.00;

@ISA       = qw(
                Exporter
               );

@EXPORT    = qw(
               );

@EXPORT_OK = qw(
                &check_format_ds_data
                &check_format_data
                &create_diff
                &create_stats
                &get_file_type
                &get_tols
                &interpolate
                &merge_stats
                &parse_args
                &print_abs_rel_stats
                &print_gnuplot_data
                &print_data_file
                &run_gnuplot
                &read_file
                &which_exec
                $GDATA
                $GDATASET_NAMES
                $GNAME
                $GORG
                $GCOORDX
                $GCOORDY
                $GDIFF
                $GMAX
                $GMEAN
                $GMIN
                $GRMS
                $GMAXABS
                $GSKIP
                $GSUMSQ
                $GTOL_A
                $GTOL_OR
                $GTOL_R
                $GREL
                $GABS
                $GDEFAULT
                $GNUMNTRUE
                $GNUMNFALSE
                $GNUMNUMS
                $GNUMSTRUE
                $GNUMSFALSE
                $GNUMSTRS
                $GNUMTRUE
                $GNUMFALSE
                $GNUMALL
                $GNUMBER_REGEXP
                $GFTCTS
                $GFTKEYWORD
                $GFTOXY
                $GFTPOP
                $GFTA
                $GFTARES
                $GFTTABLE
                $GFTTOKEN
                $GFTTRACER
                $GFTLINK
                $GFTXY
                $GCOPY_FORMAT
                $GCOPY_REGEXP
               );
sub BEGIN
  {
  }
sub END
  {
  }
#...................
#...global values...
#...................
#...data ref labels...
$GDATASET_NAMES = "Dataset_Names";
$GDATA = "Data";
$GNAME = "Name";
$GCOORDX = "X";
$GCOORDY = "Y";
$GORG  = "Orig";
$GDIFF = "Diff";
#...stat of counts...
$GNUMNTRUE  = "NumNTrue";
$GNUMNFALSE = "NumNFalse";
$GNUMNUMS   = "NumNums";
$GNUMSTRUE  = "NumSTrue";
$GNUMSFALSE = "NumSFalse";
$GNUMSTRS   = "NumStrings";
$GNUMTRUE   = "NumTrue";
$GNUMFALSE  = "NumFalse";
$GNUMALL    = "NumAll";
#...stat of other vals...
$GMAX       = "Max";
$GMAX_MATE  = "Max_mate";
$GMEAN      = "Mean";
$GMIN       = "Min";
$GMIN_MATE  = "Min_mate";
$GRMS       = "RMS";
$GMAXABS    = "MaxABS";
$GMAXABS_MATE = "MaxABS_mate";
$GSUMSQ     = "SumSq";
#...names...
$GTOL_A     = "a";
$GTOL_OR    = "or";
$GTOL_R     = "r";
$GREL       = "Rel Diff";
$GABS       = "Abs Diff";
$GDS        = "ds";
$GDEFAULT   = "default";
#...file types...
$GFT        = "ft";
$GFTCTS     = "cts";
$GFTKEYWORD = "keyword";
$GFTOXY     = "oxy";
$GFTPLOT_OUTPUT = "plot_output";
$GFTPOP     = "pop";
$GFTA       = "a";
$GFTARES    = "ares";
$GFTTABLE   = "table";
$GFTTOKEN   = "token";
$GFTTRACER  = "tracer";
$GFTLINK    = "link";
$GFTXY      = "xy";
#...multiple copy strings...
$GCOPY_FORMAT = '... copy %3d';
$GCOPY_REGEXP = ' \.\.\.\ copy\s*[0-9]+';
#...values larger than this are treated as strings...
$GSKIP = 8e99;
#...regexp to match number - simple for speed
$GNUMBER_REGEXP = '[+-]?\.?[0-9]+\.?[0-9]*([eE][+-]?\d+)?';
#...oxy strings...
$GOXY_TAG_START = "Final state ASCII diagnostic dump start";
#............................................................................
#...Name
#...====
#... check_format_ds_data
#...
#...Purpose
#...=======
#... checks the format for hash to be sent to print_gnuplot_data
#...
#...Arguments
#...=========
#... $ds_data_ref    Intent: out
#...                 Perl type: reference to hash
#...                 Obtained after read_file.
#...
#... $ierr           Intent: out
#...                 Perl type: reference to hash
#...                 Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) check format for ds_data
#............................................................................
sub check_format_ds_data
  {
    my(
       $ds_data_ref,
       $variable,
       $source,
       $gnuplot_info_ref
      ) = @_;
    my(
       @format, # correct format for data
       $ierr, # error ret value
       $key_org, # key under org
       $dtype, # diff type
       $coord, # coordinate
      );
    $ierr = 0;
    @format =
      (
       'Format (when given ds_data_ref = data{ds_name}): ',
       ' If X data exists (eg for xy file types):',
       '   $ds_data_ref{$GCOORDX}{$GNAME}           = X name',
       '   $ds_data_ref{$GCOORDX}{$GORG}[]          = X data',
       '   $ds_data_ref{$GCOORDX}{$DIFF}{<dtype>}[] = X diff data',
       ' Y data (always):',
       '   $ds_data_ref{$GCOORDY}{$GNAME}           = Y name',
       '   $ds_data_ref{$GCOORDY}{$GORG}[]          = Y data',
       '   $ds_data_ref{$GCOORDY}{$DIFF}{<dtype>}[] = Y diff data',
      );
    if( ref( $ds_data_ref ) ne "HASH" )
      {
        $ierr = 0;
        &print_error( "Internal Error in check_format_ds_data",
                      "ref(ds_data) ne 'HASH' of coord",
                      @format,
                      $ierr );
        $ierr = 1;
        return( $ierr );
      }
    if( ! defined( $$ds_data_ref{$GCOORDY} ) )
      {
        $ierr = 0;
        &print_error( "Internal Error in check_format_ds_data",
                      "ds_data{$GCOORDY} must be defined",
                      @format,
                      $ierr );
        $ierr = 1;
        return( $ierr );
      }
    foreach $coord ( sort keys %$ds_data_ref )
      {
        if( $coord ne $GCOORDX && $coord ne $GCOORDY )
          {
            $ierr = 0;
            &print_error( "Internal Error in check_format_ds_data",
                          "ds_data{$coord != $GCOORDX or $GCOORDY}",
                          @format,
                          $ierr );
            $ierr = 1;
            return( $ierr );
          }
        if( ref( $$ds_data_ref{$coord} ) ne "HASH" )
          {
            $ierr = 0;
            &print_error( "Internal Error in check_format_ds_data",
                          "ref(ds_data{$coord}) ne 'HASH' of $GORG",
                          @format,
                          $ierr );
            $ierr = 1;
            return( $ierr );
          }
        if( ! defined($$ds_data_ref{$coord}{$GNAME}) ||
            ref( $$ds_data_ref{$coord}{$GNAME} ) )
          {
            $ierr = 0;
            &print_error( "Internal Error in check_format_ds_data",
                          "ds_data{$coord}{$GNAME} must be name",
                          @format,
                          $ierr );
            $ierr = 1;
            return( $ierr );
          }
        if( ! defined($$ds_data_ref{$coord}{$GORG}) ||
            ref( $$ds_data_ref{$coord}{$GORG} ) ne "ARRAY" )
          {
            $ierr = 0;
            &print_error( "Internal Error in check_format_ds_data",
                          "ds_data{$coord}{$GORG} ne 'ARRAY' of values",
                          @format,
                          $ierr );
            $ierr = 1;
            return( $ierr );
          }
        foreach $key_org ( keys %{$$ds_data_ref{$coord}} )
          {
            if( $key_org ne "$GORG" &&
                $key_org ne "$GNAME" &&
                $key_org ne "$GDIFF" )
              {
                $ierr = 0;
                &print_error( "Internal Error in check_format_ds_data",
                              "ds_data{$coord}{$key_org != $GORG or $GNAME or $GDIFF}",
                              @format,
                              $ierr );
                $ierr = 1;
                return( $ierr );
              }
            if( $key_org eq "$GDIFF" )
              {
                if( ref( $$ds_data_ref{$coord}{$key_org} ) ne "HASH" )
                  {
                    $ierr = 0;
                    &print_error( "Internal Error in check_format_ds_data",
                                  "ref(ds_data{$coord}{$key_org} ne 'HASH' of diff type",
                                  @format,
                                  $ierr );
                    $ierr = 1;
                    return( $ierr );
                  }
                foreach $dtype ( keys %{$$ds_data_ref{$coord}{$key_org}} )
                  {
                    if( ref( $$ds_data_ref{$coord}{$key_org}{$dtype} ) ne
                        "ARRAY" )
                      {
                        $ierr = 0;
                        &print_error( "Internal Error in check_format_ds_data",
                                      "ref(ds_data{$coord}{$key_org}{$dtype} ne 'ARRAY' of values",
                                      @format,
                                      $ierr );
                        $ierr = 1;
                        return( $ierr );
                      }
                  }
              }
          }
      }
  }
#............................................................................
#...Name
#...====
#... check_format_data
#...
#...Purpose
#...=======
#... checks the format for hash returned by read_file
#...
#...Arguments
#...=========
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              Obtained after read_file.
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) check format for data
#............................................................................
sub check_format_data
  {
    my(
       $data_ref,
      ) = @_;
    my(
       $ds_name, # dataset name
       @ds_names, # listing of ds_name s
       @format, # correct format for data
       $ierr, # error return value
       $string_1, # a string for matching
       $string_2, # a string for matching
      );
    @format =
      (
       'Format: ',
       ' $$data_ref{$GDATA}{<dataset name>} = Data for dataset',
      );
    $ierr = 0;
    if( defined( $$data_ref{$GDATA} ) && ref( $$data_ref{$GDATA} ) ne "HASH" )
      {
        $ierr = 0;
        &print_error( "Internal Error in check_format_data",
                      'ref($$data_ref{$GDATA}) ne "HASH" of ds_name',
                      @format,
                      $ierr );
        $ierr = 1;
        return( $ierr );
      }
    if( defined( $$data_ref{$GDATA} ) &&
        ref( $$data_ref{$GDATASET_NAMES} ) ne "ARRAY" )
      {
        $ierr = 0;
        &print_error( "Internal Error in check_format_data",
                      'ref($$data_ref{$GDATASET_NAMES}) ne "ARRAY" of ds_name s',
                      @format,
                      $ierr );
        $ierr = 1;
        return( $ierr );
      }
    if( defined( $$data_ref{$GDATA} ) )
      {
        @ds_names = sort keys %{$$data_ref{$GDATA}};
        $string_1 = join( '|', @ds_names );
        $string_2 = join( '|', sort @{$$data_ref{$GDATASET_NAMES}} );
        if( $string_1 ne $string_2 )
          {
            $ierr = 0;
            &print_error( "Internal Error in check_format_data",
                          '$$data_ref{$GDATASET_NAMES}[] should be an ordered listing of the keys of $$data_ref{$GDATA}{}',
                          "DATASET_NAMES = [$string_1]",
                          "DATA          = [$string_2]",
                          $ierr );
            $ierr = 1;
            return( $ierr );
          }
      }
    foreach $ds_name ( sort keys %{$$data_ref{$GDATA}} )
      {
        $ierr = &check_format_ds_data( $$data_ref{$GDATA}{$ds_name} );
        if( $ierr )
          {
            $ierr = 0;
            &print_error( "Internal Error in check_format_data",
                          "Error from check_format_ds_data (data_ref{$GDATA}{$ds_name})",
                          $ierr );
            $ierr = 1;
            return( $ierr );
          }
      }
  }
#............................................................................
#...Name
#...====
#... create_diff
#...
#...Purpose
#...=======
#... Create a diff hash from 2 arrays and command line arguments
#...    Rel Diff - relative difference
#...    Abs Diff - absolute difference
#... If either value is not a number, a string comparison is made
#... and a non-empty string is placed in Rel/Abs position if
#..  a difference is detected (and $diff_result is incremented)
#...
#...Arguments
#...=========
#... $array_base_ref   Intent: in
#...                   Perl type: reference to array
#...                   base array
#...
#... $array_new_ref    Intent: in
#...                   Perl type: reference to array
#...                   new array array
#...
#... $cmd_ref          Intent: in
#...                   Perl type: reference to hash
#...                   Gotten from parse_args()
#...
#... $var_name         Intent: in
#...                   Perl type: scalar
#...                   For getting variable specific tolerance from cmd_ref.
#...
#... $title            Intent: in
#...                   Perl type: scalar
#...                   For printing - label to print next to dataset.
#...
#... $print            Intent: in
#...                   Perl type: scalar
#...                   If differences should be printed as they are found.
#...                   Non-"0" will print.
#...
#... $array_x_ref      Intent: in
#...                   Perl type: scalar
#...                   If you want to print out the X value that is the
#...                   coordinate to the corresponding Y values being diffed
#...                   You might want to do this if, for instance,
#...                   interpolation was done.
#...
#... $diff_ref         Intent: out
#...                   Perl type: reference to stat hash
#...                   given tolerences, assigns
#...                     Rel Diff: Relative    Difference
#...                     Abs Diff: Subtraction Difference
#...                   Can be used with other routines (eg print_gnuplot_data).
#...
#... $diff_result_ref  Intent: out
#...                   Perl type: reference to scalar
#...                   Number of differences (0 if no difference)
#...                   This values is also the return value.
#...
#...Program Flow
#...============
#... 1) Compute differences
#............................................................................
sub create_diff
  {
    my(
       $array_base_ref,
       $array_new_ref,
       $cmd_ref,
       $var_name,
       $title,
       $print,
       $array_x_ref,
       $diff_ref,
       $diff_result_ref,
      ) = @_;
    my(
       %array_stats, # mean, ...
       @diff_index, # which index has diff
       $i, # loop var
       $is_diff, # is a difference
       $is_diff_a, # is an absolute difference
       $is_diff_r, # is a  relative difference
       $mean, # mean of base array
       $num_vals, # number of values to diff
       $num_dashes, # number of dashes to print
       $pstring, # print string
       $rel_ref, # reference to rel array to push
       %tols, # the tolerance to use (hash)
       $tol_abs, # abs tol
       $tol_abs_d, # if defined
       $tol_rel, # rel tol
       $tol_rel_d, # if defined
       $tol_or_d, # if or defined
       $val1, # a value
       $val2, # a value
       $val_abs, # absolute difference
       $val_rel, # relative difference
       @vals, # array of values
      );
    #..........
    #...init...
    #..........
    delete( $$diff_ref{$GABS} );
    delete( $$diff_ref{$GREL} );
    $$diff_result_ref = 0;
    $num_vals = $#{$array_base_ref} > $#{$array_new_ref} ?
      $#{$array_base_ref} : $#{$array_new_ref};
    #.....................................
    #...special quick diff if same data...
    #.....................................
    if( $array_new_ref == $array_base_ref )
      {
        for( $i = 0; $i <= $num_vals; $i++ )
          {
            push( @{$$diff_ref{$GABS}}, 0 );
            push( @{$$diff_ref{$GREL}}, 0 );
          }
        return( $$diff_result_ref );
      }
    #....................
    #...get tolerances...
    #....................
    get_tols( $cmd_ref, $var_name, \%tols );
    $tol_abs = $tols{$GTOL_A};
    $tol_rel = $tols{$GTOL_R};
    if( defined( $tol_abs ) )
      {
        $tol_abs_d = 1;
      }
    else
      {
        $tol_abs_d = 0;
        $tol_abs   = 0;
      }
    if( defined( $tol_rel ) )
      {
        $tol_rel_d = 1;
      }
    else
      {
        $tol_rel_d = 0;
        $tol_rel   = 0;
      }
    if( ( defined($$cmd_ref{$GTOL_OR}) ) ||
        (   $tol_abs_d && ! $tol_rel_d ) ||
        ( ! $tol_abs_d &&   $tol_rel_d ) )
      {
        $tol_or_d = 1;
      }
    else
      {
        $tol_or_d = 0;
      }
    # compute mean to be used when finding rel tol of nums with 0
    &create_stats( $array_base_ref, \%array_stats );
    $mean = $array_stats{Mean};
    #.......................................
    #...compute diffs for each data point...
    #.......................................
    undef( @diff_index );
    for( $i = 0; $i <= $num_vals; $i++ )
      {
        $val1   = $$array_base_ref[$i];
        $val2   = $$array_new_ref[$i];
        $is_diff   = 0;
        #.........................................
        #...if both are numbers, use tolerances...
        #.........................................
        if( defined( $val1 ) && $val1 =~ /^$GNUMBER_REGEXP$/ &&
            abs( $val1 ) < $GSKIP &&
            defined( $val2 ) && $val2 =~ /^$GNUMBER_REGEXP$/ &&
            abs( $val2 ) < $GSKIP )
          {
            $val_abs = $val2 - $val1;
            # if one is 0, and the absolute difference is less than one
            # set that to the relative difference.
            # This is to reduce the number of 0 and .00000001 differences
            if( abs($val_abs) > 0 &&    # some difference
                ( $val1 == 0 || $val2 == 0 ) && # either is zero
                abs($val_abs) <= $mean/5 ) # val_abs is small compared to mean
            {
                $val_rel = $val_abs/(abs($val_abs) + abs($mean));
            }
            elsif( abs( $val1 ) > $tol_abs )
              {
                $val_rel = $val_abs / abs( $val1 );
              }
            elsif( abs($val2) > $tol_abs )
              {
                $val_rel = $val_abs / abs( $val2 );
              }
            else
              {
                $val_rel = 0;
              }
            #....................
            #...is_diff_<type>...
            #....................
            if( abs( $val_abs ) > $tol_abs )
              {
                $is_diff_a = 1;
              }
            else
              {
                $is_diff_a = 0;
              }
            if( abs( $val_rel ) > $tol_rel )
              {
                $is_diff_r = 1;
              }
            else
              {
                $is_diff_r = 0;
              }
            #............................................
            #...see if difference depending upon flags...
            #............................................
            #................
            #...or defined...
            #................
            if( $tol_or_d )
              {
                if( $is_diff_r && $is_diff_a )
                  {
                    $is_diff = 1;
                  }
              }
            #.................
            #...default and...
            #.................
            else
              {
                if( $is_diff_a || $is_diff_r )
                  {
                    $is_diff = 1;
                  }
              }
            #...................................
            #...if different, record values  ...
            #...Record 0s if within tolerance...
            #...................................
            if( $is_diff != 0 )
              {
                $$diff_result_ref += 1;
                push( @{$$diff_ref{$GABS}}, $val_abs );
                push( @{$$diff_ref{$GREL}}, $val_rel );
              }
            else
              {
                push( @{$$diff_ref{$GABS}}, 0 );
                push( @{$$diff_ref{$GREL}}, 0 );
              }
          }
        #...............................................
        #...DONE: if both are numbers, use tolerances...
        #...............................................
        #.....................................
        #...if not numbers, compare strings...
        #.....................................
        else
          {
            if( !defined( $val1 ) || !defined( $val2 ) || $val1 ne $val2 )
              {
                $is_diff = 1;
                $$diff_result_ref += 1;
                push( @{$$diff_ref{$GABS}}, "StringDiff" );
                push( @{$$diff_ref{$GREL}}, "StringDiff" );
              }
            else
              {
                push( @{$$diff_ref{$GABS}}, "" );
                push( @{$$diff_ref{$GREL}}, "" );
              }
          }
        push( @diff_index, $is_diff );
      }
    #.............................................
    #...DONE: compute diffs for each data point...
    #.............................................
    #........................
    #...print if requested...
    #........................
    if( $print != 0 && $$diff_result_ref )
      {
        printf( "\n        [%s]\n", $title );
        printf( " %12s [%16s] [%16s] [%16s] %12s %12s\n",
                "Index", "X", "Base", "New", $GABS, $GREL );
        $num_dashes = 1+12+1+16+2+1+16+2+1+16+2+1+12+1+12;
        printf( "%s\n", "v"x($num_dashes) );
        for( $i = 0; $i <= $num_vals; $i++ )
          {
            if( $diff_index[$i] != 0 )
              {
                printf( ">" );
                printf( "%12s", $i+1 );
                $val = $$array_x_ref[$i];
                ($pstring, $val) = print_val( $val, 16, 8 );
                printf( " [$pstring]", $val );
                $val = $$array_base_ref[$i];
                ($pstring, $val) = print_val( $val, 16, 8 );
                printf( " [$pstring]", $val );
                $val = $$array_new_ref[$i];
                ($pstring, $val) = print_val( $val, 16, 8 );
                printf( " [$pstring]", $val );
                $val = $$diff_ref{$GABS}[$i];
                ($pstring, $val) = print_val( $val, 13, 4 );
                printf( "$pstring", $val );
                $val = $$diff_ref{$GREL}[$i];
                ($pstring, $val) = print_val( $val, 13, 4 );
                printf( "$pstring", $val );
                printf( "\n" );
              }
          }
        printf( "%s\n", "^"x($num_dashes) );
        printf( " %12s [%16s] [%16s] [%16s] %12s %12s\n",
                "Index", "X", "Base", "New", $GABS, $GREL );
      }
    #..............................
    #...DONE: print if requested...
    #..............................
    $$diff_result_ref;
  }
sub print_val
  {
    my(
       $val, # value to print
       $len, # length to take up
       $dec  # digits after decimal (negative to print integer)
      ) = @_;
    my(
       $pstring, # print string
      );
    if( defined( $val ) && $val =~ /^$GNUMBER_REGEXP$/ &&
        abs( $val ) < $GSKIP )
      {
        if( $dec >= 0 )
          {
            $pstring = sprintf( "%%%d.%de", $len, $dec );
          }
        else
          {
            $pstring = sprintf( "%%%dd", $len, $dec );
          }
      }
    elsif( defined( $val ) )
      {
        $pstring = sprintf( "%%%ds", $len );
      }
    else
      {
        $val = "-";
        $pstring = sprintf( "%%%ds", $len );
      }
    return( $pstring, $val );
  }
#............................................................................
#...Name
#...====
#... create_stats
#...
#...Purpose
#...=======
#... Create stats hash from an array
#... Unless otherwise specified, non-number values will be skipped.
#...    Max           - max value
#...    Max_MATE      - mate of max value: The absolute and relative errors of a particular data point are mates.
#...    MaxABS        - max absolute value
#...    MaxABS_MATE   - mate of max absolute value
#...    Mean          - mean value
#...    Min           - min value
#...    Min_MATE      - mate of min value
#...    RMS           - root mean square
#...    SumSq         - Sum of the squares of the numbers
#...    NumNTrue      - Number of non-0 numbers
#...    NumNFalse     - Number of 0 numbers
#...    NumNums       - Number of Numbers
#...    NumSTrue      - Number of non-empty strings
#...    NumSFalse     - Number of empty strings
#...    NumStrs       - Number of strings
#...    NumTrue       - Number of true numbers and strings
#...    NumFalse      - Number of false numbers and strings
#...    NumAll        - Total count
#...
#...Arguments
#...=========
#... $array_ref   Intent: in
#...              Perl type: reference to array
#...
#... $stat_ref    Intent: out
#...              Perl type: reference to stat hash
#...
#...Program Flow
#...============
#... 1) Compute stats
#............................................................................
sub create_stats
  {
    my(
       $array_ref,
       $stat_ref,
       $arrayb_ref        # "mate" data set
      ) = @_;
    my(
       $i, # loop variable
       $num_elements, # number of elements in array
       $val, # value of an array elem
       $val1, # value
       $val2, # value
       );
    #..........
    #...init...
    #..........
    undef( %{$stat_ref} );
    $num_elements = $#{$array_ref};
    #...............
    #...init vals...
    #...............
    for( $i = 0; $i <= $num_elements; $i++ )
      {
        $val = $$array_ref[$i];
        if( $val =~ /^$GNUMBER_REGEXP$/ &&
            abs( $val ) < $GSKIP )
          {
            $$stat_ref{$GMAX}  = $val;
            $$stat_ref{$GMAX_MATE}  = $$arrayb_ref[$i];
            $$stat_ref{$GMIN}  = $val;
            $$stat_ref{$GMIN_MATE}  = $$arrayb_ref[$i];
	    $$stat_ref{$GMAXABS_MATE} = $$arrayb_ref[$i];
            last;
          }
      }
    $$stat_ref{$GNUMALL}    = $num_elements + 1;
    $$stat_ref{$GNUMNTRUE}  = 0;
    $$stat_ref{$GNUMNFALSE} = 0;
    $$stat_ref{$GNUMNUMS}   = 0;
    $$stat_ref{$GNUMSTRUE}  = 0;
    $$stat_ref{$GNUMSFALSE} = 0;
    $$stat_ref{$GNUMSTRS}   = 0;
    $$stat_ref{$GNUMTRUE}   = 0;
    $$stat_ref{$GNUMFALSE}  = 0;
    $$stat_ref{$GMEAN} = 0;
    $$stat_ref{$GSUMSQ} = 0;
    #................................
    #...loop setting various stats...
    #................................
    for( $i = 0; $i <= $num_elements; $i++ )
      {
        $val = $$array_ref[$i];
        if( $val =~ /^$GNUMBER_REGEXP$/ &&
            abs( $val ) < $GSKIP )
          {
            #$$stat_ref{$GMAX}    = $$stat_ref{$GMAX} > $val ? $$stat_ref{$GMAX} : $val;
	    if ( $$stat_ref{$GMAX} < $val) {
	      $$stat_ref{$GMAX} = $val;
	      $$stat_ref{$GMAX_MATE} = $$arrayb_ref[$i];
	    }
            #$$stat_ref{$GMIN}    = $$stat_ref{$GMIN} < $val ? $$stat_ref{$GMIN} : $val;
	    if ( $$stat_ref{$GMIN} > $val) {
	      $$stat_ref{$GMIN} = $val;
	      $$stat_ref{$GMIN_MATE} = $$arrayb_ref[$i];
	    }
            $$stat_ref{$GMEAN}  += $val;
            $$stat_ref{$GSUMSQ} += $val**2;
            $$stat_ref{$GNUMNUMS}++;
            if( $val != 0 )
              {
                $$stat_ref{$GNUMNTRUE}++;
              }
          }
        else
          {
            $$stat_ref{$GNUMSTRS}++;
            if( defined( $val ) && length( $val ) > 0 )
              {
                $$stat_ref{$GNUMSTRUE}++;
              }
          }
      }
    $$stat_ref{$GNUMNFALSE} = $$stat_ref{$GNUMNUMS}  - $$stat_ref{$GNUMNTRUE};
    $$stat_ref{$GNUMSFALSE} = $$stat_ref{$GNUMSTRS}  - $$stat_ref{$GNUMSTRUE};
    $$stat_ref{$GNUMTRUE}   = $$stat_ref{$GNUMNTRUE} + $$stat_ref{$GNUMSTRUE};
    $$stat_ref{$GNUMFALSE}  = $$stat_ref{$GNUMALL}   - $$stat_ref{$GNUMTRUE};
    #.......................
    #...delete if not set...
    #.......................
    if( $$stat_ref{$GNUMNUMS} == 0 )
      {
        delete( $$stat_ref{$GMEAN} );
        delete( $$stat_ref{$GSUMSQ} );
      }
    #...........................
    #...finish off some stats...
    #...........................
    if( $$stat_ref{$GNUMNUMS} > 0 )
      {
        $$stat_ref{$GMEAN} = ($$stat_ref{$GMEAN})/$$stat_ref{$GNUMNUMS};
        $val1 = abs( $$stat_ref{$GMIN} );
        $val2 = abs( $$stat_ref{$GMAX} );

        #$$stat_ref{$GMAXABS} = $val1 > $val2 ? $val1 : $val2;
	if ($val1 > $val2) {
	  $$stat_ref{$GMAXABS} = $val1;
	  $$stat_ref{$GMAXABS_MATE} = $$stat_ref{$GMIN_MATE};
	}
	else {
	  $$stat_ref{$GMAXABS} = $val2;
	  $$stat_ref{$GMAXABS_MATE} = $$stat_ref{$GMAX_MATE};
	}
        if( $$stat_ref{$GNUMNUMS} > 1 )
          {
            $$stat_ref{$GRMS} = 
              ( ( $$stat_ref{$GSUMSQ} -
                  $$stat_ref{$GNUMNUMS} * $$stat_ref{$GMEAN}**2 ) ** .5 ) /
                    ( $$stat_ref{$GNUMNUMS} - 1 );
          }
        else
          {
            $$stat_ref{$GRMS} = -1;
          }
      }
    #.......................................................................
    #...rms - not used since slow - although less loss or arith precision...
    #.......................................................................
    #if( $$stat_ref{$GNUMNUMS} > 0 )
    #  {
    #    $$stat_ref{$GRMS}  = 0;
    #    for( $i = 0; $i <= $num_elements; $i++ )
    #      {
    #        $val = $$array_ref[$i];
    #        if( $val =~ /^$GNUMBER_REGEXP$/ && abs( $val ) < $GSKIP )
    #          {
    #            $$stat_ref{$GRMS} += ($val - $$stat_ref{$GMEAN})**2;
    #          }
    #      }
    #    if( $$stat_ref{$GNUMNUMS} > 1 )
    #      {
    #        $$stat_ref{$GRMS} =
    #          ((($$stat_ref{$GRMS})**.5)/($$stat_ref{$GNUMNUMS}-1));
    #      }
    #    else
    #      {
    #        $$stat_ref{$GRMS} = -1;
    #      }
    #  }
 }
#............................................................................
#...Name
#...====
#... get_file_type
#...
#...Purpose
#...=======
#... Get the filt type (xy, table, token, ... )
#... A file type match is done on the file name.
#... If not found, try to get file type by sampling a block of lines
#... and seeing what they look like.  If enough lines match the
#... threshhold, the file type is assigned.
#... If still no match, set file type to token.
#...
#...Arguments
#...=========
#... $file_name     Intent: in
#...                Perl type: scalar
#...                File name to read in
#...
#... $file_type_ref Intent: out
#...                Perl type: reference to scalar
#...                string of file type.
#...
#... $ierr          Intent: out
#...                Perl type: scalar
#...                Error return value (non-0 for error)
#...
#...Program Flow
#...============
#... 1) detect file type by name
#... 2) detect file type by contents
#............................................................................
sub get_file_type
  {
    my(
       $file_name,
       $file_type_ref,
      ) = @_;
    my(
       $count,  # count matches
       $count1, # count matches
       $num_ds, # number of datasets
       $num_fields, # number of fields
       $err_msg, # error message
       $found, # if found something
       $i, # loop var
       $ierr, # error ret val
       $ierr_sys, # error ret val from system call
       $j, # loop var
       $line, # line of file
       $lines_read, # number of lines to read in file (neg for whole file)
       @lines, # lines of file without blank lines
       @lines_orig, # lines of file untouched
       $threshhold, # threshhold ratio for match
       $threshhold_actual, # actual data gotten for threshhold
       @tokens, # tokens on a line
      );
    #..........
    #...init...
    #..........
    $ierr = 0;
    $threshhold = .7;
    $threshhold_actual = 1;
    $lines_read = 10;
    undef( $$file_type_ref );
    #................................
    #...open FILE, read some lines...
    #................................
    if( ! defined( $file_name ) )
      {
        $ierr = 0;
        &cts_diff_util::print_error( "Filename not defined.",
                                     $ierr );
        return( 1 );
      }
    if( ! open( FILE, $file_name ) )
      {
        $ierr = 0;
        &cts_diff_util::print_error( "Cannot open file [$file_name]",
                                     $ierr );
        return( 1 );
      }
    while( ( $lines_read < 0 || $#lines + 1 < $lines_read ) &&
           ( $line = <FILE> ) )
      {
        push( @lines_orig, $line );
        if( $line =~ /\S/ )
          {
            push( @lines, $line );
          }
      }
    close( FILE );
    #............................
    #...file line: plot_output...
    #............................
    #...first line looks like a data block
    if( ! defined( $$file_type_ref ) && $#lines > 0 )
      {
        if( $lines[0] =~ /^#\s*\[\d+\]\s*\S+/ )
          {
            $$file_type_ref = $GFTPLOT_OUTPUT;
          }
      }
    #...................
    #...file name: xy...
    #...................
    if( ! defined( $$file_type_ref ) && $file_name =~ /\.xy(\.std)?$/ )
      {
        $$file_type_ref = $GFTXY;
      }
    #...................
    #...file line: xy...
    #...................
    #... # <ds name>
    #... <val1>  <val2>   (most of the lines fit this)
    if( ! defined( $$file_type_ref ) && $#lines > 0 )
      {
        #.........................
        #...matches # <ds name>...
        #.........................
        if( $lines[0] =~ /^\s*\#\s*\S+/ )
          {
            $count = 0;
            #...........................................
            #...update count if matches <val1> <val2>...
            #...or another ds header                 ...
            #...........................................
            for( $i = 1; $i <= $#lines; $i++ )
              {
                $line = $lines[$i];
                $line =~ s/^\s*(.*?)\s*$/$1/;
                @tokens = split( /\s+/, $line );
                if( ( $line =~ /^\#\s*\S+/ && $count > 0 ) ||
                    ( $#tokens == 1 &&
                      $tokens[0] =~ /^$GNUMBER_REGEXP$/ &&
                      abs( $tokens[0] ) < $GSKIP &&
                      $tokens[1] =~ /^$GNUMBER_REGEXP$/ &&
                      abs( $tokens[1] ) < $GSKIP ) )
                  {
                    $count++;
                  }
              }
            $threshhold_actual = $count/$#lines;
            if( $threshhold_actual >= $threshhold )
              {
                $$file_type_ref = $GFTXY;
              }
          }
      }
    #......................
    #...file line: table...
    #......................
    #... <ds1 name>  <ds2 name>  <ds3 name>
    #... <val ds1>   <val ds2>   <val ds3>  (most of the lines fit this)
    if( ! defined( $$file_type_ref ) && $#lines > 0 )
      {
        #............................
        #...get number of datasets...
        #............................
        $line = $lines[0];
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        @tokens = split( /\s+/, $line );
        $num_ds = $#tokens+1;
        if( $num_ds > 0 )
          {
            #.....................................................
            #...matches <val ds1> <val ds2> ... <val ds$num_ds>...
            #.....................................................
            $count = 0;
            for( $i = 1; $i <= $#lines; $i++ )
              {
                $line = $lines[$i];
                $line =~ s/^\s*//;
                $line =~ s/\s*$//;
                @tokens = split( /\s+/, $line );
                if( $#tokens + 1 == $num_ds )
                  {
                    $count1 = 0;
                    for( $j = 0; $j <= $#tokens; $j++ )
                      {
                        if( $tokens[$j] =~ /^$GNUMBER_REGEXP$/ &&
                            abs( $tokens[$j] ) < $GSKIP )
                          {
                            $count1++;
                          }
                      }
                    if( $count1/$num_ds >= $threshhold )
                      {
                        $count++;
                      }
                  }
              }
            $threshhold_actual = $count/$#lines;
            if( $threshhold_actual >= $threshhold )
              {
                $$file_type_ref = $GFTTABLE
              }
          }
      }
    #....................
    #...file line: pop...
    #....................
    #...First line says something like "now in PoP"...
    if( ! defined( $$file_type_ref ) && $#lines > 0 )
      {
        #.........................
        #...matches # <ds name>...
        #.........................
        if( $lines[0] =~ /\s+now in PoP\s+/ )
          {
            $$file_type_ref = $GFTPOP;
          }
      }
    #....................
    #...file line: cts...
    #....................
    if( ! defined( $$file_type_ref ) && $#lines >= 0 )
      {
        #...................
        #...matches # cts...
        #...................
        if( $lines[0] =~ /^\# cts\n$/ )
          {
            $$file_type_ref = $GFTCTS;
          }
      }
    #..................
    #...file line: a...
    #..................
    #...within a few lines, says key phrase...
    if( ! defined( $$file_type_ref ) && $#lines_orig > 0 )
      {
        #........................................
        #...look for phrase within a few lines...
        #........................................
        for( $i = 0; $i+3 <= $#lines_orig; $i++ )
          {
            if( $lines_orig[$i] =~ /^\s*Version:/ )
              {
                if( $lines_orig[$i+1] =~ /^\s*Copyright/ &&
                    $lines_orig[$i+3] =~ /^\s*Los Alamos National Laboratory/ )
                  {
                    $$file_type_ref = $GFTA;
                    last;
                  }

                if( $lines_orig[$i+5] =~ /^\s*Copyright/ &&
                    $lines_orig[$i+6] =~ /^\s*Los Alamos National Laboratory/ )
                  {
                    $$file_type_ref = $GFTA;
                    last;
                  }
              }
          }
      }
    #.....................
    #...file line: ares...
    #.....................
    #...within a few lines, says key phrase...
    if( ! defined( $$file_type_ref ) && $#lines >= 0 )
      {
        #............................................
        #...first non-blank line must match format...
        #............................................
        if( $lines[0] =~
            /
            ^\s*                    # start with possible whitespace
            (\S+)\s+                # test
            (DIFF|FAILED|PASSED)\s+ # test results
            (P-[0-9]+),\s+          # P field,
            (D-[0-9]+),\s+          # D field,
            (F-[0-9]+)\s+           # F field
            $/x )
          {
            $$file_type_ref = $GFTARES;
          }
      }
    #.......................
    #...file name: tracer...
    #.......................
    if( ! defined( $$file_type_ref ) && $file_name =~ /\-tracer$/i )
      {
        $$file_type_ref = $GFTTRACER;
      }
    #.......................
    #...file line: tracer...
    #.......................
    #... ds_name_1, ds_name_2, ... , ds_name_n
    #...   somewhere in there is "particle" and "time" fields
    #... val_1,     val_2, ... ,     val_n
    #... # <ds name>
    #... <val1>  <val2>   (most of the lines fit this)
    if( ! defined( $$file_type_ref ) && $#lines > 0 )
      {
        $count  = 0;
        $num_ds = 0;
        $found = "false";
        #...check lines
        for( $i = 0; $i <= $#lines; $i++ )
          {
            #...stop if hit enough matches
            if( $count > 10 )
              {
                last;
              }
            #...remove leading/trailing whitespace and comment lines
            $line = $lines[$i];
            $line =~ s/^\s*(.*\S)\s*$/$1/;
            $line =~ s/\s*\#.*//;
            if( $line =~ /\S/ )
              {
                @tokens = split( /\s*,\s*/, $line );
                $num_fields = $#tokens + 1;
                if( $found eq "false" )
                  {
                    $found = "true";
                    $count1 = 0;
                    $num_ds = $num_fields;
                    #...must have certain fields
                    if( grep( /^particle$/, @tokens ) )
                      {
                        $count1++;
                      }
                    if( grep( /^time$/, @tokens ) )
                      {
                        $count1++;
                      }
                    if( $count1 != 2 )
                      {
                        last;
                      }
                    next;
                  }
                # if you hit a line that does not match num_ds,
                # this is not tracer
                if( $num_ds != $num_fields )
                  {
                    $count = 0;
                    last;
                  }
                else
                  {
                    $count++;
                  }
              }
          }
        if( $count >= 1 )
          {
            $$file_type_ref = $GFTTRACER;
          }
      }
    #.....................
    #...file name: link...
    #.....................
    if( ! defined( $$file_type_ref ) && $file_name =~ /\.lnk(\.\d+)?$/ )
      {
        $$file_type_ref = $GFTLINK;
      }
    #....................
    #...file grep: oxy...
    #....................
    #>>>HAVE THIS LAST SINCE GREPPING THE WHOLE FILE<<<
    #>>>DO THIS SEARCH IF GOT A FUZZY MATCH ABOVE<<<
    if( ! ( defined( $$file_type_ref ) && $threshhold_actual eq "1" ) &&
        $#lines >= 0 )
      {
        $ierr_sys = system( "grep", "-q", $GOXY_TAG_START, $file_name );
        if( ! $ierr_sys )
          {
            $$file_type_ref = $GFTOXY;
          }
      }
    #...............
    #...otherwise...
    #...............
    if( ! defined( $$file_type_ref ) )
      {
        $$file_type_ref = $GFTTOKEN;
      }
    #...........................
    #...close file and return...
    #...........................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... get_tols
#...
#...Purpose
#...=======
#... Get the values for the tolerances given the inputs
#...
#...Arguments
#...=========
#... $cmd_ref          Intent: in
#...                   Perl type: reference to hash
#...                   Gotten from parse_args()
#...
#... $ds_name          Intent: in
#...                   Perl type: scalar
#...                   Name of the dataset
#...
#... $tols_ref         Intent: out
#...                   Perl type: reference to hash
#...                   $tols{$GABS, $GREL, ... }
#...                   Tolerances hash.
#...
#...Program Flow
#...============
#... 1) Compute differences
#............................................................................
sub get_tols
  {
    my(
       $cmd_ref,
       $ds_name,
       $tols_ref,
      ) = @_;
    my(
       $ds_name_orig, # ds name without any copy extension
       $ds_regexp_key, # key to hash
       $tol_abs, # abs tol
       $tol_rel, # rel tol
      );
    if( defined( $ds_name ) )
      {
        ( $ds_name_orig = $ds_name ) =~ s/$GCOPY_REGEXP$//;
      }
    $tol_abs = $$cmd_ref{$GTOL_A}{$GDEFAULT};
    foreach $ds_regexp_key ( keys %{$$cmd_ref{$GTOL_A}{$GDS}} )
      {
        if( defined( $ds_name_orig ) && $ds_name_orig =~ /^($ds_regexp_key)$/ )
          {
            $tol_abs = $$cmd_ref{$GTOL_A}{$GDS}{$ds_regexp_key};
            last;
          }
      }
    $tol_rel = $$cmd_ref{$GTOL_R}{$GDEFAULT};
    foreach $ds_regexp_key ( keys %{$$cmd_ref{$GTOL_R}{$GDS}} )
      {
        if( defined( $ds_name_orig ) && $ds_name_orig =~ /^($ds_regexp_key)$/ )
          {
            $tol_rel = $$cmd_ref{$GTOL_R}{$GDS}{$ds_regexp_key};
            last;
          }
      }
    $$tols_ref{$GTOL_A} = $tol_abs;
    $$tols_ref{$GTOL_R} = $tol_rel;
  }
#............................................................................
#...Name
#...====
#... interpolate
#...
#...Purpose
#...=======
#... This converts 1 set of data points to another using the a set of X
#... values.
#...     (base_x, base_y) -> (intp_x, intp_y)
#...    in:
#...     base_x, base_y, intp_x
#...    out:
#...     inpt_y
#... Currently linear interpolation is done.
#...
#...Arguments
#...=========
#... $base_x_ref       Intent: in
#...                   Perl type: reference to array
#...                   The base X values.
#...                   These values must be monotonically increasing
#...
#... $base_y_ref       Intent: in
#...                   Perl type: reference to array
#...                   The base Y values.
#...
#... $intp_x_ref       Intent: in
#...                   Perl type: reference to array
#...                   The values to interpolate to.
#...                   These values must be monotonically increasing
#...
#... $intp_y_ref       Intent: out
#...                   Perl type: reference to array
#...                   The base Y values.
#...
#...Program Flow
#...============
#... 1) Foreach value in intp_x_ref
#... 1.1) Find 2 points in the base set
#... 1.2) compute slope m
#... 1.3) intp_y = m*(intp_x - base_x) + base_y
#............................................................................
sub interpolate
  {
    my(
       $base_x_ref,
       $base_y_ref,
       $intp_x_ref,
       $intp_y_ref
      ) = @_;
    my(
       $base_x, # copy of base_x value
       $done, # if done in a loop
       $equal, # quick exit if same x vals
       $final_pair, # if no more pairs of points will be found - no search
       $i, # loop var
       $idx_a, # index of first point
       $idx_a_prev, # previous value for it
       $idx_b, # index of next point
       $idx_b_prev, # previous value for it
       $intp_x, # a value of ref
       $intp_y, # a value of ref
       $m, # slope
       $tmp_a, # tmp value for printing
       $tmp_b, # tmp value for printing
       $val_prev, # previous value for testing monatomically increasing
       $val_curr, # current value for testing monatomically increasing
      );
    $ierr = 0;
    #.....................................................
    #...test for monatomically increasing base_x values...
    #.....................................................
    undef( $val_prev );
    for( $i = 0; $i <= $#$base_x_ref; $i++ )
      {
        $val_curr = $$base_x_ref[$i];
        if( $val_curr !~ /^$GNUMBER_REGEXP$/ || abs($val_curr) >= $GSKIP )
          {
            next;
          }
        if( defined( $val_prev ) && $val_prev > $val_curr )
          {
            $ierr = 1;
            $tmp_b = $i + 1;
            &print_error( "Base X-values must not be decreasing to interpolate",
                          "Prev: [$tmp_a:$val_prev] > Current [$tmp_b:$val_curr]",
                          $ierr );
            return( $ierr );
          }
        $tmp_a = $i+1;
        $val_prev = $val_curr;
      }
    #........................................................
    #...test that length of base_x and base_y arrays match...
    #........................................................
    if( $#$base_x_ref != $#$base_y_ref )
      {
        $ierr = 1;
        &print_error( "Length of Base X array [$#$base_x_ref] and ".
                      "Base Y array [$#$base_y_ref] must match",
                      $ierr );
        return( $ierr );
      }
    # quick exit if the x vals are the same
    if( $#$base_x_ref == $#$intp_x_ref ){
        $equal = "true";
        for( $i = 0; $i <= $#$base_x_ref; $i++ ){
            if( $$base_x_ref[$i] ne $$intp_x_ref[$i] ){
                $equal = "false";
                last;
            }
        }
        if( $equal eq "true" ){
            @$intp_y_ref = @$base_y_ref;
            return( $ierr );
        }
    }
    #.....................................................
    #...test for monatomically increasing intp_x values...
    #.....................................................
    undef( $val_prev );
    foreach $val_curr ( @$intp_x_ref )
      {
        if( $val_curr !~ /^$GNUMBER_REGEXP$/ || abs($val_curr) >= $GSKIP )
          {
            next;
          }
        #.......................................
        #...test for monatomically increasing...
        #.......................................
        if( defined( $val_prev ) && $val_prev >= $val_curr )
          {
            $ierr = 1;
            &print_error( "Intp X-values must be increasing to interpolate",
                          "Prev: [$val_prev] >= Current [$val_curr]",
                          $ierr );
            return( $ierr );
          }
        $val_prev = $val_curr;
      }
    #................................................................
    #...get to the first pair of points where $$base_x_ref differs...
    #................................................................
    $idx_a = 0;
    while( $idx_a <= $#$base_x_ref )
      {
        if( $$base_x_ref[$idx_a] =~ /^$GNUMBER_REGEXP$/ &&
            abs($$base_x_ref[$idx_a]) < $GSKIP )
          {
            last;
          }
        $idx_a++;
      }
    $idx_b = $idx_a + 1;
    while( $idx_b <= $#$base_x_ref )
      {
        if( $$base_x_ref[$idx_b] =~ /^$GNUMBER_REGEXP$/ &&
            abs($$base_x_ref[$idx_b]) < $GSKIP )
          {
            if( $$base_x_ref[$idx_a] != $$base_x_ref[$idx_b] )
              {
                last;
              }
            else
              {
                $idx_a = $idx_b;
              }
          }
        $idx_b++;
      }
    if( $idx_b > $#$base_x_ref )
      {
        $idx_b = $idx_a;
      }
    $idx_a_prev = $idx_a;
    $idx_b_prev = $idx_b;
    $final_pair = 0;
    #................................................
    #...compute intp_y foreach value in intp_x_ref...
    #................................................
    foreach $intp_x ( @$intp_x_ref )
      {
        #................................................
        #...if intp_x not valid, set intp_y to invalid...
        #................................................
        if( $intp_x !~ /^$GNUMBER_REGEXP$/ || abs($intp_x) >= $GSKIP ||
            $idx_a > $#$base_x_ref )
          {
            $intp_y = "undef";
          }
        #........................................
        #...linear interpolation to get intp_y...
        #........................................
        else
          {
            $done = 0;
            #.........................................................
            #...Find 2 points in the base set to be used for interp...
            #.........................................................
            while( $done == 0 && $final_pair == 0 )
              {
                #.........................................................
                #...if bounded by idx_b, done since any subsequent pair...
                #...of points will be further away than this pair      ...
                #.........................................................
                if( $intp_x <= $$base_x_ref[$idx_b] )
                  {
                    $done = 1;
                  }
                #......................
                #...get to next pair...
                #......................
                else
                  {
                    #.................................................
                    #...reset previous values - increment a if != b...
                    #.................................................
                    $idx_b_prev = $idx_b;
                    if( $$base_x_ref[$idx_a] != $$base_x_ref[$idx_b] )
                      {
                        $idx_a_prev = $idx_a;
                      }
                    #................................................
                    #...point a to b and b starting guess at a + 1...
                    #................................................
                    $idx_a = $idx_b;
                    $idx_b = $idx_a + 1;
                    while( $idx_b <= $#$base_x_ref )
                      {
                        if( $$base_x_ref[$idx_b] =~ /^$GNUMBER_REGEXP$/ &&
                            abs($$base_x_ref[$idx_b]) < $GSKIP )
                          {
                            last;
                          }
                        $idx_b++;
                      }
                    #.....................................................
                    #...if no next pair, then                          ...
                    #...  done, reset to previous point, and final_pair...
                    #.....................................................
                    if( $idx_b > $#$base_x_ref )
                      {
                        $final_pair = 1;
                        $done = 1;
                        $idx_b = $idx_b_prev;
                        $idx_a = $idx_a_prev;
                      }
                  }
                #............................
                #...DONE: get to next pair...
                #............................
              }
            #...............................................................
            #...DONE: Find 2 points in the base set to be used for interp...
            #...............................................................
            #...................................................
            #...if found at least 1 valid number intp on that...
            #...................................................
            if( $idx_a <= $#$base_x_ref )
              {
                #...............................................
                #...if base_y are numbers, interpolate       ...
                #...base_x already must be numbers from above...
                #...............................................
                if( $$base_y_ref[$idx_a] =~ /^$GNUMBER_REGEXP$/ &&
                    abs( $$base_y_ref[$idx_a] ) < $GSKIP &&
                    $$base_y_ref[$idx_b] =~ /^$GNUMBER_REGEXP$/ &&
                    abs( $$base_y_ref[$idx_b] ) < $GSKIP )
                  {
                    #.....................
                    #...compute slope m...
                    #.....................
                    if( $idx_a != $idx_b )
                      {
                        $m =
                          ($$base_y_ref[$idx_b] - $$base_y_ref[$idx_a])/
                            ($$base_x_ref[$idx_b] - $$base_x_ref[$idx_a]);
                      }
                    else
                      {
                        $m = 0;
                      }
                    #...........................................
                    #...intp_y = m*(intp_x - base_x) + base_y...
                    #...........................................
                    $intp_y = $m*($intp_x - $$base_x_ref[$idx_a]) +
                      $$base_y_ref[$idx_a];
                  }
                #..............................................
                #...if non-number, pick closest base_y value...
                #..............................................
                else
                  {
                    if( abs( $intp_x - $$base_x_ref[$idx_a] ) <
                        abs( $intp_x - $$base_x_ref[$idx_b] ) )
                      {
                        $intp_y = $$base_y_ref[$idx_a];
                      }
                    else
                      {
                        $intp_y = $$base_y_ref[$idx_b];
                      }
                  }
              }
            #.........................................................
            #...DONE: if found at least 1 valid number intp on that...
            #.........................................................
            else
              {
                $intp_y = "undef";
              }
          }
        #..............................................
        #...DONE: linear interpolation to get intp_y...
        #..............................................
        push( @$intp_y_ref, $intp_y );
      }
    #......................................................
    #...DONE: compute intp_y foreach value in intp_x_ref...
    #......................................................
  }
#............................................................................
#...Name
#...====
#... merge_stats
#...
#...Purpose
#...=======
#... merge first stat hash into second one (eg Max is max of first and
#... second one).
#... If the second stat is not defined, this is effectively a copy op.
#...
#...Arguments
#...=========
#... $stat_ref_a  Intent: in
#...              Perl type: reference to hash
#...              Created from create_stats
#...
#... $stat_ref_b  Intent: inout
#...              Perl type: reference to hash
#...              Created from create_stats
#...              Contains merging of stat_ref_a and stat_ref_b
#...
#...Notes
#...=====
#... Some values cannot be merged:
#...   - RMS -> take max RMS value
#...
#...Program Flow
#...============
#... 1) merge
#............................................................................
sub merge_stats
  {
    my(
       $stat_ref_a,
       $stat_ref_b
      ) = @_;
    my(
       $key, # key for hash
       %orig_b, # original stat_ref_b
       $sumsq_a, # value for a stat
       $sumsq_b, # value for a stat
       $val_a, # a value
       $val_b, # a value
     );
    foreach $key ( keys %{$stat_ref_b} )
      {
        $orig_b{$key} = $$stat_ref_b{$key};
      }
    foreach $key ( $GNUMNTRUE,
                   $GNUMNFALSE,
                   $GNUMNUMS,
                   $GNUMSTRUE,
                   $GNUMSFALSE,
                   $GNUMSTRS,
                   $GNUMTRUE,
                   $GNUMFALSE,
                   $GNUMALL,
                   $GSUMSQ )
      {
        $val_a = $$stat_ref_a{$key};
        $val_b = $orig_b{$key};
        if( defined( $val_a ) && defined( $val_b ) )
          {
            $$stat_ref_b{$key} = $val_a + $val_b;
          }
        elsif( defined( $val_a ) )
          {
            $$stat_ref_b{$key} = $val_a;
          }
        else
          {
            $$stat_ref_b{$key} = $val_b;
          }
      }
    foreach $key ( $GMIN, $GMAX, $GMAXABS )
      {
	my $key_mate = $key . "_mate";
        $val_a = $$stat_ref_a{$key};
        $val_b = $orig_b{$key};
        if( defined( $val_a ) && defined( $val_b ) )
          {
            #$$stat_ref_b{$key} = $val_a > $val_b ? $val_a : $val_b;
	    if (abs($val_a) > abs($val_b) ) {
	      $$stat_ref_b{$key} = $val_a;
	      $$stat_ref_b{$key_mate} =  $$stat_ref_a{$key_mate};
	    }
          }
        elsif( defined( $val_a ) )
          {
            $$stat_ref_b{$key} = $val_a;
	    $$stat_ref_b{$key_mate} =  $$stat_ref_a{$key_mate};
          }
#         else
#           {
#             $$stat_ref_b{$key} = $val_b;
#           }
      }
#     foreach $key ( $GMIN )
#       {
#         $val_a = $$stat_ref_a{$key};
#         $val_b = $orig_b{$key};
#         if( defined( $val_a ) && defined( $val_b ) )
#           {
#             $$stat_ref_b{$key} = $val_a < $val_b ? $val_a : $val_b;
#           }
#         elsif( defined( $val_a ) )
#           {
#             $$stat_ref_b{$key} = $val_a;
#           }
#         else
#           {
#             $$stat_ref_b{$key} = $val_b;
#           }
#       }
    #...mean (uses new numnums)...
    $val_a = $$stat_ref_a{$GMEAN};
    $val_b = $orig_b{$GMEAN};
    if( defined( $val_a ) && defined( $val_b ) )
      {
        $$stat_ref_b{$GMEAN} =
          ($val_a*$$stat_ref_a{$GNUMNUMS} + $val_b*$orig_b{$GNUMNUMS})/
            ($$stat_ref_b{$GNUMNUMS});
      }
    elsif( defined( $val_a ) )
      {
        $$stat_ref_b{$GMEAN} = $val_a;
      }
    else
      {
        $$stat_ref_b{$GMEAN} = $val_b;
      }
    #...rms (do with new numnums, and mean)...
    $val_a = $$stat_ref_a{$GRMS};
    $val_b = $orig_b{$GRMS};
    if( defined( $val_a ) && defined( $val_b ) )
      {
        $$stat_ref_b{$GRMS} =
          ((($$stat_ref_b{$GSUMSQ}-
             $$stat_ref_b{$GNUMNUMS}*
             ($$stat_ref_b{$GMEAN}**2)))**.5)/
              ($$stat_ref_b{$GNUMNUMS}-1);
      }
    elsif( defined( $val_a ) )
      {
        $$stat_ref_b{$GRMS} = $val_a;
      }
    else
      {
        $$stat_ref_b{$GRMS} = $val_b;
      }
  }
#............................................................................
#...Name
#...====
#... parse_args
#...
#...Purpose
#...=======
#... Create cmd hash from command line
#...
#...Arguments
#...=========
#... $argv_ref    Intent: out
#...              Perl type: reference to array
#...              \@ARGV usually
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line
#...              $cmd{$option} = value
#...              $cmd{files}[] = array of file names
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub parse_args
  {
    my(
       $argv_ref,
       $cmd_ref
      ) = @_;
    my(
       $arg_file, # name for argument file
       @args, # arguments
       $ierr, # error ret val
       $line, # line of file
       $num_args, # number of arguments
       $opt, # current option
       @opts, # array of options
       @tokens, # split of line
       @vals, # array of values
       $val, # value for current option
      );
    $ierr = 0;
    undef( %{$cmd_ref} );
    @args = @{$argv_ref};
    #...........................
    #...read in argument file...
    #...........................
    $arg_file = "cts_diff.arg";
    if( -T $arg_file )
      {
        $ierr = &read_arg_file( $arg_file, \@args );
        if( $ierr )
          {
            $ierr = 0;
            &print_error( "Failure reading argument file [$arg_file].",
                          $ierr );
            return( 1 );
          }
      }
    #....................
    #...parse the args...
    #....................
    $num_args = $#args;
    while( @args )
      {
        $opt = shift( @args );
        #......................
        #...-(a|r) tolerance...
        #......................
        if( $opt =~ /^-+(a|r)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 0;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                return( 1 );
              }
            @vals = split( /\s*,\s*/, shift( @args ) );
            $val = shift( @vals );
            if( $val !~ /^$GNUMBER_REGEXP$/ || abs( $val ) >= $GSKIP )
              {
                $ierr = 0;
                &print_error( "Must give numeric value for option [-$opt $val]",
                              $ierr );
                return( 1 );
              }
            if( $val < 0 )
              {
                $ierr = 0;
                &print_error( "Value for option [-$opt $val] must be in the range (0,).",
                              $ierr );
                return( 1 );
              }
            #...............................
            #...place val or list of vals...
            #...............................
            if( $#vals < 0 )
              {
                $$cmd_ref{$opt}{$GDEFAULT} = $val;
              }
            else
              {
                grep( s/(\||\[|\]|\{|\}|\(|\)|\$|\@|\%|\*|\.)/\\$1/g,
                      @vals );
                $$cmd_ref{$opt}{$GDS}{join( '|', @vals )} = $val;
              }
          }
        #...................
        #...-ft file_type...
        #...................
        elsif( $opt =~ /^-+(ft)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 0;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                return( 1 );
              }
            $val = shift( @args );
            if( $val !~ /$GFTKEYWORD|$GFTA|$GFTCTS|$GFTOXY|$GFTPLOT_OUTPUT|$GFTPOP|$GFTTABLE|$GFTTOKEN|$GFTTRACER|$GFTLINK|$GFTXY/ )
              {
                $ierr = 0;
                &print_error( "Invalid file_type [-$opt $val]",
                              $ierr );
                return( 1 );
              }
            $$cmd_ref{$opt} = $val;
          }
        #...............................................
        #...-(or|h|help|no_intp|no_plots|pft|presult)...
        #...............................................
        elsif( $opt =~ /^-+(or|help|h|no_intp|no_plots|pft|plot_orig|presult|nostatus)$/ )
          {
            $opt = $1;
            if( $opt eq "h" )
              {
                $opt = "help";
              }
            $$cmd_ref{$opt} = "true";
          }
        #....................
        #...-ds <datasets>...
        #....................
        elsif( $opt =~ /^-+($GDS)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 0;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                return( 1 );
              }
            $val = shift( @args );
            $val =~ s/^\s*//;
            $val =~ s/\s*$//;
            #...stuff into regexp (ds1|ds2|...|dsn)...
            @vals = split( /\s*,\s*/, $val );
            grep( s/(\||\[|\]|\{|\}|\(|\)|\$|\@|\%|\*|\.)/\\$1/g, @vals );
            $$cmd_ref{$opt} = join( '|', @vals );
          }
        #....................
        #...-v|fsets <num>...
        #....................
        elsif( $opt =~ /^-+(v|fsets)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 0;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                return( 1 );
              }
            $val = shift( @args );
            if( $val !~ /^[0-9]+$/ )
              {
                $ierr = 0;
                &print_error( "Must give non-negative integer value for option [-$opt $val]",
                              $ierr );
                return( 1 );
              }
            $$cmd_ref{$opt} = $val;
          }
        #.......................................
        #...-o_<file_type> <output data file>...
        #.......................................
        elsif( $opt =~ /^-+(o_(\S+))$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 0;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                return( 1 );
              }
            $$cmd_ref{$opt} = shift( @args );
          }
        #.....................
        #...-arg <filename>...
        #.....................
        elsif( $opt =~ /^-+(arg)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 0;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                return( 1 );
              }
            $val = shift( @args );
            $ierr = &read_arg_file( $val, \@args );
            if( $ierr )
              {
                $ierr = 0;
                &print_error( "Failure reading argument file [-$opt $val].",
                              $ierr );
                return( 1 );
              }
          }
        #..............
        #...filename...
        #..............
        else
          {
            push( @{$$cmd_ref{files}}, $opt );
          }
      }
    #..........................
    #...DONE: parse the args...
    #..........................
    #....................
    #...exclusive opts...
    #....................
    @opts = grep( /^($GFTKEYWORD|$GFTCTS|$GFTOXY|$GFTPOP|$GFTA|$GFTARES|$GFTTABLE|$GFTTOKEN|$GFTTRACER|$GFTLINK|$GFTXY)$/,
                  sort keys %{$cmd_ref} );
    if( $#opts >= 1 )
      {
        $ierr = 0;
        &print_error( "The following options are exclusive:",
                      "[".join(" ", @opts)."]",
                      $ierr );
        return( 1 );
      }
    #.................
    #...check fsets...
    #.................
    if( defined($$cmd_ref{files}) )
      {
        $val1 = $#{$$cmd_ref{files}} + 1;
      }
    else
      {
        $val1 = 0;
      }
    $val2 = $$cmd_ref{fsets};
    if( defined( $val2 ) && $val2 > 0 && int( $val1/$val2 ) != $val1/$val2 )
      {
        $ierr = 0;
        &print_error( "The number of files [$val1] must be",
                      "divisible by the file set size [-fsets $val2]",
                      $ierr );
        return( 1 );
      }
    return( $ierr );
  }
#............................................
#...read_arg_file                         ...
#...  Little routine to read argument file...
#............................................
sub read_arg_file
  {
    my(
       $arg_file,
       $args_ref
      ) = @_;
    my(
       $ierr, # error ret val
       $line, # line of file
       @tokens, # tokens in file
      );
    $ierr = 0;
    if( ! open( FILE, $arg_file ) )
      {
        $ierr = 0;
        &print_error( "Cannot open argument file [$arg_file]",
                      $ierr );
        return( 1 );
      }
    undef( @tokens );
    while( $line=<FILE> )
      {
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        if( $line =~ /^\s*$/ ||
            $line =~ /^\s*\#/ )
          {
            next;
          }
        push( @tokens, split( /\s+/, $line ) );
      }
    if( @tokens )
      {
        unshift( @{$args_ref}, @tokens );
      }
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... print_abs_rel_stats
#...
#...Purpose
#...=======
#... Print a few of the stats gotten from create_stats
#... Spacing tied to create_diff printing
#...
#...Arguments
#...=========
#... $title        Intent: in
#...               Perl type: scalar
#...               title of the stat line
#...
#... $print_header Intent: in
#...               Perl type: scalar
#...               0: no header
#...               1: header
#...
#... $stat_abs_ref Intent: out
#...               Perl type: reference to hash
#...               gotten from create_stats on absolute difference array
#...
#... $stat_rel_ref Intent: out
#...               Perl type: reference to hash
#...               gotten from create_stats on relative difference array
#...
#... $headers_printed_ref Intent: in
#...                      Perl type: reference to scalar
#...                      adds 1 to value if header printed.
#...
#...Program Flow
#...============
#... 1) print a few stats
#............................................................................
sub print_abs_rel_stats
  {
    my(
       $title,
       $print_header,
       $stat_abs_ref,
       $stat_rel_ref,
       $headers_printed_ref,
      ) = @_;
    my(
       $percent, # percent different
      );
    if( defined( $$stat_abs_ref{$GNUMTRUE} ) && $$stat_abs_ref{$GNUMTRUE} > 0 )
      {
        if( $print_header )
          {
            printf( "\n %27s %27s %12s %12s (%6s%%) [%s]\n",
                    "Max Abs (associated rel)", "Max Rel (associated abs)", "String Diffs", "Diffs Total",
                    "Diff ", "Title" );
            printf( " %12s %12s %12s %12s (%6s-) [%s]\n",
                    "-"x27, "-"x27, "-"x12, "-"x12,
                    "-"x6, "-"x5 );
            $$headers_printed_ref += 1;
          }
        if( $$stat_abs_ref{$GNUMALL} > 0 )
          {
            $percent =
              ($$stat_abs_ref{$GNUMTRUE}/$$stat_abs_ref{$GNUMALL})*100.0;
          }
        else
          {
            $percent = 0;
          }
        if( $$stat_abs_ref{$GNUMNUMS} > 0 )
          {
            printf( " %12.4e \(%11.4e\)  %12.4e \(%11.4e\) %13d %12d (%6.2f%%) [%s]\n",
                    $$stat_abs_ref{$GMAXABS},
                    $$stat_abs_ref{$GMAXABS_MATE},
                    $$stat_rel_ref{$GMAXABS},
                    $$stat_rel_ref{$GMAXABS_MATE},
                    $$stat_abs_ref{$GNUMSTRUE},
                    $$stat_abs_ref{$GNUMTRUE},
                    $percent,
                    $title
                  );
          }
        else
          {
            printf( " %27s %27s %12d %12d (%6.2f%%) [%s]\n",
                    "-",
                    "-",
                    $$stat_abs_ref{$GNUMSTRUE},
                    $$stat_abs_ref{$GNUMTRUE},
                    $percent,
                    $title
                  );
          }
      }
  }
#............................................................................
#...Name
#...====
#... print_data_file
#...
#...Purpose
#...=======
#... Print the data struct to a file ($cmd{o})
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: in
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: in
#...              Perl type: reference to hash
#...              reference to data hash
#...
#...Program Flow
#...============
#... 1) call correct printing routine
#............................................................................
sub print_data_file
  {
    my(
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ierr, # error return value
       @opts, # set of options
      );
    $ierr = 0;
    @opts = grep( /^o_/, keys %$cmd_ref );
    @opts = grep( !/^o_(cts|xy)/, @opts );
    if( $#opts >= 0 )
      {
        $ierr = 1;
        cts_diff_util::print_error( "Invalid file type(s) for printing.",
                                    "[".join(" ", @opts)."]",
                                    $ierr );
        exit( $ierr );
      }
    if( defined( $$cmd_ref{o_cts} ) )
      {
        $ierr = &print_data_file_cts( $cmd_ref, $data_ref );
      }
    if( defined( $$cmd_ref{o_xy} ) )
      {
        $ierr = &print_data_file_xy( $cmd_ref, $data_ref );
      }
    if( $ierr != 0 )
      {
        $ierr = 1;
        cts_diff_util::print_error( "Error printing data file",
                                    $ierr );
        return( $ierr );
      }
  }
#............................................................................
#...Name
#...====
#... print_data_file_cts
#...
#...Purpose
#...=======
#... Print the data struct to a file ($cmd{o})
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: in
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: in
#...              Perl type: reference to hash
#...              reference to data hash
#...
#...Program Flow
#...============
#... 1) open file
#... 2) cycle through and print data
#............................................................................
sub print_data_file_cts
  {
    my(
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ierr, # error return value
       $coord, # a coordinate
       $ds_name, # dataset name
       $ds_name_print, # the name to use for printing
       @coords, # list of coordinates
       $i, # loop var
       $len, # length of something
       $line, # line to print
       $num_vals, # number of values
       $printf_format, # the format to print the values
      );
    $ierr = 0;
    if( ! %$data_ref )
      {
        return( $ierr );
      }
    if( ! open( FILE, ">$$cmd_ref{o_cts}" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open output data file [$$cmd_ref{o_cts}]",
                      $ierr );
        return( $ierr );
      }
    print FILE "# cts\n";
    foreach $ds_name ( @{$$data_ref{$GDATASET_NAMES}} )
      {
        #...get number of values to print...
        $num_vals = 0;
        undef( @coords );
        foreach $coord ( $GCOORDX, $GCOORDY )
          {
            if( @{$$data_ref{$GDATA}{$ds_name}{$coord}{$GORG}} )
              {
                push( @coords, $coord );
                $len = $#{$$data_ref{$GDATA}{$ds_name}{$coord}{$GORG}}+1;
                if( $len > $num_vals )
                  {
                    $num_vals = $len;
                  }
              }
          }
        if( $num_vals == 0 )
          {
            next;
          }
        ($ds_name_print = $ds_name) =~ s/$GCOPY_REGEXP$//;
        print FILE "# Dataset Name: $ds_name_print\n";
        $line = "# ";
        if( $#coords > 0 )
          {
            #...have space to ensure at least 1 space between items...
            $printf_format = "%-30s ";
          }
        else
          {
            $printf_format = "%-s";
          }
        foreach $coord ( @coords )
          {
            print FILE
              "# Coord Name $coord: $$data_ref{$GDATA}{$ds_name}{$coord}{$GNAME}\n";
          }
        #...print values...
        for( $i = 0; $i < $num_vals; $i++ )
          {
            $line = "  ";
            foreach $coord ( @coords )
              {
                $line .= sprintf( $printf_format,
                                  $$data_ref{$GDATA}{$ds_name}{$coord}{$GORG}[$i] );
              }
            $line =~ s/\s*$//;
            print FILE "$line\n";
          }
      }
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... print_data_file_xy
#...
#...Purpose
#...=======
#... Print the data struct to a file ($cmd{o}) in the XY format.
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: in
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: in
#...              Perl type: reference to hash
#...              reference to data hash
#...
#...Program Flow
#...============
#... 1) open file
#... 2) cycle through and print data
#............................................................................
sub print_data_file_xy
  {
    my(
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ierr, # error return value
       $coord, # a coordinate
       $ds_name, # dataset name
       $ds_name_print, # the name to use for printing
       @coords, # list of coordinates
       $i, # loop var
       $len, # length of something
       $line, # line to print
       $num_vals, # number of values
       $printf_format, # the format to print the values
      );
    $ierr = 0;
    if( ! %$data_ref )
      {
        return( $ierr );
      }
    if( ! open( FILE, ">$$cmd_ref{o_xy}" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open output data file [$$cmd_ref{o_xy}]",
                      $ierr );
        return( $ierr );
      }
    foreach $ds_name ( @{$$data_ref{$GDATASET_NAMES}} )
      {
        #...get number of values to print...
        $num_vals = 0;
        undef( @coords );
        foreach $coord ( $GCOORDX, $GCOORDY )
          {
            if( @{$$data_ref{$GDATA}{$ds_name}{$coord}{$GORG}} )
              {
                push( @coords, $coord );
                $len = $#{$$data_ref{$GDATA}{$ds_name}{$coord}{$GORG}}+1;
                if( $len > $num_vals )
                  {
                    $num_vals = $len;
                  }
              }
          }
        if( $num_vals == 0 )
          {
            next;
          }
        ($ds_name_print = $ds_name) =~ s/$GCOPY_REGEXP$//;
        print FILE "#  $ds_name_print\n";
        #...print values...
        for( $i = 0; $i < $num_vals; $i++ )
          {
            $line = " ";
            #...fill in Y-only with index...
            if( $#coords == 0 )
              {
                $line .= $i+1;
                $line .= "    ";
              }
            else
              {
                #...do nothing...
              }
            foreach $coord ( @coords )
              {
                $line .= "$$data_ref{$GDATA}{$ds_name}{$coord}{$GORG}[$i]    ";
              }
            $line =~ s/\s+$//;
            print FILE "$line\n";
          }
      }
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... print_gnuplot_data
#...
#...Purpose
#...=======
#... print a set of data obtained from read_file and possibly create_diff.
#... Updates gnuplot_info hash to be used in print_gnuplot_plots.
#...
#... Currently, the data is displayed as 1 variable per page with
#... 1 plot per diff type (diff types from ds_data).
#...   Title:   From title argument
#...   X-Label: From ds_data
#...   Y-Label: From ds_data
#...   Key:     From source argument
#...
#... A check is done to see if the $variable and $source have already been
#... processed with the %gnuplot_info.  If so, the routine returns.
#... So, change the $variable and $source name to print data.
#...
#...Arguments
#...=========
#... $ds_data_ref      Intent: in
#...                   Perl type: reference to array
#...                   Dataset data obtained from read_file with possible
#...                   calling of create_diff as well:
#...                     read_file -> %data
#...                     foreach $ds_name of %data
#...                       foreach $coord of $data{$ds_name}
#...                         create_diff $data{$ds_name}{$coord}{$GORG}[] ->
#...                                     $data{$ds_name}{$coord}{$GDIFF}{}[]
#...                         print_gnuplot_data with \%{$data{$ds_name}}
#...
#... $variable         Intent: in
#...                   Perl type: scalar
#...                   The same variable by different sources will be plotted
#...                   on the same plot.
#...                   The order of the variables is preserved.
#...
#... $source           Intent: in
#...                   Perl type: scalar
#...                   The same variable by different sources will be plotted
#...                   on the same plot.
#...                   The order of the sources is sorted.
#...                   This value will be used in the key.
#...
#... $title            Intent: in
#...                   Perl type: scalar
#...                   The title of the plot.  The title of the first source
#...                   will be used.
#...
#... $gnuplot_info_ref Intent: inout
#...                   Perl type: reference to hash
#...                   Contains info gnuplot will use to plot in run_gnuplot
#...                   routine.
#...
#...Program Flow
#...============
#... 1) print values
#... 2) store index and using info for gnuplot
#............................................................................
sub print_gnuplot_data
  {
    my(
       $ds_data_ref,
       $variable,
       $source,
       $title,
       $gnuplot_info_ref
      ) = @_;
    my(
       $column, # column
       $i, # loop var
       $line, # line to print out
       $dtype, # diff type
       @keys, # keys to a hash
       $num_vals, # number of values
       $len, # length of something
       $lena, # lengths in printing
       $lenb,
       $lenc,
       %points, # hash of y values of number of legal points
       $points_printed, # if legal points have been printed
       $update_index, # if data has been printed - index needs updating
       $val, # a value
       $valid_x, # if the x coordinate is valid
       $xy, # if doing xy plot (otherwise just y values)
      );
    #..........
    #...init...
    #..........
    if( ! defined $$gnuplot_info_ref{data_file} )
      {
        $$gnuplot_info_ref{data_file} = "cts_diff.data";
      }
    if( ! defined $$gnuplot_info_ref{num_sets} )
      {
        unlink( $$gnuplot_info_ref{data_file} );
        $$gnuplot_info_ref{num_sets} = 0;
      }
    undef( %points );
    #.................................
    #...return if already processed...
    #.................................
    if( defined( $$gnuplot_info_ref{processed}{$variable}{$source} ) )
      {
        return;
      }
    else
      {
        $$gnuplot_info_ref{processed}{$variable}{$source} = "";
      }
    #...................................
    #...check for correct ds_data_ref...
    #...................................
    $ierr = &check_format_ds_data( $ds_data_ref );
    if( $ierr )
      {
        $ierr = 0;
        &print_error( "Format error for ds_data", $ierr );
        $ierr = 1;
        return( $ierr );
      }
    #...............
    #...open file...
    #...............
    if( ! open( FILE, ">>$$gnuplot_info_ref{data_file}" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open gnuplot data file [$$gnuplot_info_ref{data_file}]",
                      $ierr );
        exit( $ierr );
      }
    #......................
    #...get type of plot...
    #......................
    if( @{$$ds_data_ref{$GCOORDX}{$GORG}} )
      {
        $xy = 1;
      }
    else
      {
        $xy = 0;
      }
    #..........................
    #...get number of values...
    #..........................
    $num_vals = 0;
    if( @{$$ds_data_ref{$GCOORDX}{$GORG}} )
      {
        $len = $#{$$ds_data_ref{$GCOORDX}{$GORG}}+1;
        if( $len > $num_vals )
          {
            $num_vals = $len;
          }
      }
    if( @{$$ds_data_ref{$GCOORDY}{$GORG}} )
      {
        $len = $#{$$ds_data_ref{$GCOORDY}{$GORG}}+1;
        if( $len > $num_vals )
          {
            $num_vals = $len;
          }
      }
    foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
      {
        $len = $#{$$ds_data_ref{$GCOORDY}{$GDIFF}{$dtype}} + 1;
        if( $len > $num_vals )
          {
            $num_vals = $len;
          }
      }
    if( $num_vals == 0 )
      {
        return;
      }
    #..................
    #...print header...
    #..................
    $lena = 15;
    $lenb = $lena + 7;
    $lenc = $lena + 4;
    $line = "";
    $line .= sprintf( "# Title: [%s] Count: [%s]\n",
                      $title, $num_vals );
    $line .= "# ";
    if( $xy )
      {
        $line .= sprintf( "%${lenc}s($GCOORDX)", $$ds_data_ref{$GCOORDX}{$GNAME} );
      }
    else
      {
        $line .= sprintf( "%${lenb}s", " " );
      }
    $line .= sprintf( " %${lenc}s($GCOORDY)", $$ds_data_ref{$GCOORDY}{$GNAME} );
    foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
      {
        $line .= sprintf( " %${lenc}s($GCOORDY)", $dtype );
      }
    $line .= "\n";
    print FILE $line;
    #................
    #...print data...
    #................
    $update_index = 0;
    $points{$GORG} = 0;
    foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
      {
        $points{$GDIFF}{$dtype} = 0;
      }
    for( $i = 0; $i < $num_vals; $i++ )
      {
        $line = "  ";
        #.......
        #...x...
        #.......
        $valid_x = 1;
        if( $xy )
          {
            $val = $$ds_data_ref{$GCOORDX}{$GORG}[$i];
            if( defined( $val ) && $val =~ /^$GNUMBER_REGEXP$/ &&
                abs( $val ) < $GSKIP )
              {
                $line .= sprintf( "%${lenb}.${lena}e", $val );
              }
            else
              {
                $valid_x = 0;
                $val = "-";
                $line .= sprintf( "%${lenb}s", $val );
              }
          }
        else
          {
            $line .= sprintf( "%${lenb}s", " " );
          }
        #.......
        #...y...
        #......
        $val = $$ds_data_ref{$GCOORDY}{$GORG}[$i];
        if( defined( $val ) && $val =~ /^$GNUMBER_REGEXP$/ &&
            abs( $val ) < $GSKIP )
          {
            if( $valid_x )
              {
                $points{$GORG} += 1;
              }
            $line .= sprintf( " %${lenb}.${lena}e", $val );
          }
        else
          {
            $val = "-";
            $line .= sprintf( " %${lenb}s", $val );
          }
        #............
        #...y diff...
        #............
        foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
          {
            $val = $$ds_data_ref{$GCOORDY}{$GDIFF}{$dtype}[$i];
            if( defined( $val ) && $val =~ /^$GNUMBER_REGEXP$/ &&
                abs( $val ) < $GSKIP )
              {
                $line .= sprintf( " %${lenb}.${lena}e", $val );
                if( $valid_x )
                  {
                    $points{$GDIFF}{$dtype} += 1;
                  }
              }
            else
              {
                $val = "-";
                $line .= sprintf( " %${lenb}s", $val );
              }
          }
        $line .= "\n";
        print FILE $line;
      }
    $points_printed = $points{$GORG};
    foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
      {
        $points_printed += $points{$GDIFF}{$dtype};
      }
    #................
    #...close file...
    #................
    print FILE "\n\n";
    close( FILE );
    #...........................
    #...store additional data...
    #...........................
    #...............................................
    #...save order of source, variable, and dtype...
    #...............................................
    if( ! defined $$gnuplot_info_ref{sources_def}{$source} )
      {
        push( @{$$gnuplot_info_ref{sources}}, $source );
        $$gnuplot_info_ref{sources_def}{$source} = "";
      }
    if( ! defined $$gnuplot_info_ref{variables_def}{$variable} )
      {
        push( @{$$gnuplot_info_ref{variables}}, $variable );
        $$gnuplot_info_ref{variables_def}{$variable} = "";
      }
    foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
      {
        if( ! defined $$gnuplot_info_ref{dtypes_def}{$dtype} )
          {
            push( @{$$gnuplot_info_ref{dtypes}}, $dtype );
            $$gnuplot_info_ref{dtypes_def}{$dtype} = "";
          }
      }
    #.......................................................
    #...gnuplot variables - only if actual points printed...
    #.......................................................
    if( $points_printed > 0 )
      {
        #..........
        #...init...
        #..........
        $$gnuplot_info_ref{index}{$variable}{$source} =
          $$gnuplot_info_ref{num_sets};
        $column = 0;
        $$gnuplot_info_ref{title}{$variable}{$source} = $title;
        #................................................
        #...line type - unique and the same per source...
        #................................................
        if( ! defined( $$gnuplot_info_ref{lt}{count}) )
          {
            $$gnuplot_info_ref{lt}{count} = 0;
          }
        if( ! defined($$gnuplot_info_ref{lt}{source}{$source}) )
          {
            $$gnuplot_info_ref{lt}{source}{$source} =
              $$gnuplot_info_ref{lt}{count} % 64 + 1;
            $$gnuplot_info_ref{lt}{count}++;
          }
        #.......
        #...x...
        #.......
        if( $xy )
          {
            $column++;
            $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDX}{$GORG} =
              "$column";
            $$gnuplot_info_ref{label}{$variable}{$source}{$GCOORDX}{$GORG} =
              $$ds_data_ref{$GCOORDX}{$GNAME};
          }
        #.......
        #...y...
        #.......
        $column++;
        if( $points{$GORG} > 0 )
          {
            $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDY}{$GORG} =
              "$column";
            $$gnuplot_info_ref{label}{$variable}{$source}{$GCOORDY}{$GORG} =
              $$ds_data_ref{$GCOORDY}{$GNAME};
          }
        #............
        #...y diff...
        #............
        foreach $dtype ( sort keys %{$$ds_data_ref{$GCOORDY}{$GDIFF}} )
          {
            $column++;
            if( $points{$GDIFF}{$dtype} > 0 )
              {
                $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDY}{$GDIFF}{$dtype} =
                  "$column";
                $$gnuplot_info_ref{label}{$variable}{$source}{$GCOORDY}{$GDIFF}{$dtype} =
                  "$$ds_data_ref{$GCOORDY}{$GNAME} [$dtype]";
              }
          }
      }
    $$gnuplot_info_ref{num_sets}++;
  }
#............................................................................
#...Name
#...====
#... read_file
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: in
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              List of datasets in order of appearance:
#...                $$data_ref{$GDATASETNAMES}[] = array of dataset names
#...              If X data exists (eg for xy file types):
#...                $$data_ref{<dataset name>}{$GCOORDX}{$GNAME}  = X name
#...                $$data_ref{<dataset name>}{$GCOORDX}{$GORG}[] = X data
#...              Y data (always):
#...                $$data_ref{<dataset name>}{$GCOORDY}{$GNAME}  = Y name
#...                $$data_ref{<dataset name>}{$GCOORDY}{$GORG}[] = Y data
#...              Corresponding <dataset name> arrays eventually will be
#...              diffed using the specific tolerances given for the $GNAME
#...              value.
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) determine file type and pass to correct reading routine
#............................................................................
sub read_file
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ierr, # error return value
       $file_type, # type of file
       $file_type_cmd, # file type defined via argument
      );
    #..........
    #...init...
    #..........
    $ierr = 0;
    #...................
    #...get file type...
    #...................
    $file_type = $$cmd_ref{ft};
    if( ! defined( $file_type ) )
      {
        $ierr = cts_diff_util::get_file_type( $file_name, \$file_type );
        if( $ierr != 0 )
          {
            $ierr = 1;
            &print_error( "Error getting type of file",
                          $ierr );
            return( 1 );
          }
      }
    #..................................
    #...read file based on file type...
    #..................................
    if( $file_type eq "$GFTXY" )
      {
        $ierr = &read_file_xy( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTCTS" )
      {
        $ierr = &read_file_cts( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTTABLE" )
      {
        $ierr = &read_file_table( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTKEYWORD" )
      {
        $ierr = &read_file_keyword( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTOXY" )
      {
        $ierr = &read_file_oxy( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTPLOT_OUTPUT" )
      {
        $ierr = &read_file_plot_output( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTPOP" )
      {
        $ierr = &read_file_pop( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTA" )
      {
        $ierr = &read_file_a( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTARES" )
      {
        $ierr = &read_file_ares( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTTRACER" )
      {
        $ierr = &read_file_tracer( $file_name, $cmd_ref, $data_ref );
      }
    elsif( $file_type eq "$GFTLINK" )
      {
        $ierr = &read_file_link( $file_name, $cmd_ref, $data_ref );
      }
    else
      {
        $ierr = &read_file_token( $file_name, $cmd_ref, $data_ref );
      }
    if( $ierr )
      {
        &print_error( "Error in read_file_<type> = $file_type",
                      $ierr );
        return( 1 );
      }
    #...........................
    #...check for consistency...
    #...........................
    $ierr = &check_format_data( $data_ref );
    if( defined( $$cmd_ref{pft} ) )
      {
        printf( STDERR "File Type: $file_type [$file_name]\n" );
      }
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_a
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) stuff the specific lines into %data
#............................................................................
sub read_file_a
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name_full, # corresponding arrays with this name are diffed
       $ds_name_tol,  # dataset name used for diff tolerances
       $i,            # loop variable
       $ierr,         # error return value
       $key,          # key/val of token
       $line,         # line of file
       $name_a,       # part of a name
       $name_b,       # part of a name
       $name_c,       # part of a name
       $name_d,       # part of a name
       $name_e,       # part of a name
       $name_f,       # part of a name
       $name_g,       # part of a name
       $name_h,       # part of a name
       @tokens,       # split on whitespace of $line
       $token,        # single token,
       $time,         # the time
       $block,        # block to use
       $val,          # key/val of token
      );
    $ierr = 0;
    $time = -1; # less than 0 just in case...
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #........................
    #...look at every line...
    #........................
    while( $line=<FILE> )
      {
        #..................................
        #...one possible cycle/time line...
        #..................................
        if( $line =~
            /
            ^\s*            # whitespace
            cycle\s*=\s*    # cycle keyword
            (\d+)           # cycle number
            ,\s*time\s*=\s* # time keyword
            (\S+)           # time
            ,               # comma then other fields
            /x )
          {
            $time = $2;
          }
        #..................................
        #...one possible cycle/time line...
        #..................................
        if( $line =~
            /
            ^\s+        # begins with whitespace
            \#\#\#      # 3 pounds
            \s+         # whitespace
            (\d+)       # cycle
            \s+         # whitespace
            (\S+)       # time
            /x )
          {
            $time = $2;
          }
        #......................
        #...sum/ratio blocks...
        #......................
        if( defined( $line ) &&
            $line =~
            /^
            (\S+)\s+      # name
            RATIOS:-{10,} # RATIOS:<at least 10 "-">
            /x )
          {
            $name_a = $1;
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $name_a !~ /^($$cmd_ref{$GDS})$/ )
              {
                while( $line=<FILE> )
                  {
                    if( $line =~ /^-{10,}/ )
                      {
                        last;
                      }
                  }
                next;
              }
            #..............................................
            #...go through lines in this sum/ratio block...
            #..............................................
            while( 1==1 )
              {
                #................
                #...end marker...
                #................
                if( ! defined( $line ) || $line =~ /^-{10,}/ )
                  {
                    last;
                  }
                #................
                #...get to SUM...
                #................
                while( $line=<FILE> )
                  {
                    #........................................................
                    #...grab all the qualifiers that will make the ds_name...
                    #........................................................
                    if( $line =~
                        /
                        ^(\S+)\s+(SUMS):  # name name
                        .*\#\s+(\d+)\s+   # name
                        .*\#\s+(\d+)      # name
                        /x )
                      {
                        $name_b = $1;
                        $name_c = $2;
                        $name_d = $3;
                        $name_e = $4;
                        last;
                      }
                  }
                #..................
                #...consume sums...
                #..................
                while( $line=<FILE> )
                  {
                    if( $line !~ /\S/ )
                      {
                        last;
                      }
                    while( $line =~ s/^\s*(\S+)\s*sum\s*=\s*(\S+)// )
                      {
                        $name_f = $1;
                        $name_g = $2;
                        $ds_name_tol  =
                          "${name_a}_${name_b}_${name_c}_".
                           "${name_d}_${name_e}_${name_f}";
                        $ds_name_full = $ds_name_tol;
                        if( ! defined( $$data_ref{$GDATA}{$ds_name_full} ) )
                          {
                            push( @{$$data_ref{$GDATASET_NAMES}},
                                  $ds_name_full );
                            $$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GNAME} =
                              $ds_name_tol;
                            $$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GNAME} =
                              "time";
                          }
                        push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GORG}},
                              $name_g );
                        push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GORG}},
                              $time );
                      }
                  }
                #...................
                #...get to RATIOS...
                #...................
                while( $line=<FILE> )
                  {
                    if( $line =~ /^(\S+) (RATIOS):\s/ )
                      {
                        $name_b = $1;
                        $name_c = $2;
                        last;
                      }
                  }
                #....................
                #...consume ratios...
                #....................
                $name_h = 0;
                while( $line=<FILE>  )
                  {
                    if( $line =~ /^-{10,}/ || $line !~ /\S/ )
                      {
                        last;
                      }
                    $name_h++;
                    while( $line =~ s/^\s*(\S+\/\S+)\s*=\s*(\S+)// )
                      {
                        $name_f = $1;
                        $name_g = $2;
                        $ds_name_tol  =
                          "${name_a}_${name_b}_${name_c}_".
                           "${name_d}_${name_e}_${name_f}_${name_h}";
                        $ds_name_full = $ds_name_tol;
                        if( ! defined( $$data_ref{$GDATA}{$ds_name_full} ) )
                          {
                            push( @{$$data_ref{$GDATASET_NAMES}},
                                  $ds_name_full );
                            $$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GNAME} =
                              $ds_name_tol;
                            $$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GNAME} =
                              "time";
                          }
                        push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GORG}},
                              $name_g );
                        push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GORG}},
                              $time );
                      }
                  }
                #..........................
                #...DONE: consume ratios...
                #..........................
              }
            #..............................................
            #...go through lines in this sum/ratio block...
            #..............................................
          }
        #............................
        #...DONE: sum/ratio blocks...
        #............................
        #................
        #... ## blocks...
        #................
        if( defined( $line ) &&
            $line =~
            /^
            \s+           # some whitespace
            \#\#(\S+)     # 2 pounds then the Type
            /x )
          {
            $name_a = $1;
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $name_a !~ /^($$cmd_ref{$GDS})$/ )
              {
                next;
              }
            @tokens = split( /\s*,\s*/, $line );
            foreach $token ( @tokens )
              {
                if( $token =~
                    /
                    (\S+)     # keyword
                    \s*=\s*   # equals
                    (\S+)     # value
                    /x )
                  {
                    $ds_name_tol  = "${name_a}_$1";
                    $ds_name_full = $ds_name_tol;
                    $val = $2;
                    $val =~ s/\s*\[.*\]\s*//; # remove possible trailing chars
                    # skip some datasets
                    if( $ds_name_full =~ /^${name_a}_diff$/ ||
                        $ds_name_full =~ /^${name_a}_%_rel_diff$/ )
                      {
                        next;
                      }
                    if( ! defined( $$data_ref{$GDATA}{$ds_name_full} ) )
                      {
                        push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name_full );
                        $$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GNAME} =
                          $ds_name_tol;
                        $$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GNAME} =
                          "time";
                      }
                    push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GORG}},
                          $val );
                    push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GORG}},
                          $time );
                  }
              }
          }
        #.....................
        #...DONE: ## blocks...
        #.....................
        #...................
        #... delta blocks...
        #...................
        if( defined( $line ) &&
            $line =~
            /^
            \s+          # some whitespace
            (delta\S+)   # Type is "delta<stuff>"
            /x )
          {
            $name_a = $1;
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $name_a !~ /^($$cmd_ref{$GDS})$/ )
              {
                next;
              }
            @tokens = split( /\s*,\s*/, $line );
            foreach $token ( @tokens )
              {
                if( $token =~
                    /
                    (\S+)     # keyword
                    \s*=\s*   # equals
                    (\S+)     # value
                    /x )
                  {
                    $ds_name_tol  = "${name_a}_$1";
                    $ds_name_full = $ds_name_tol;
                    $val = $2;
                    $val =~ s/\s*\[.*\]\s*//; # remove possible trailing chars
                    if( ! defined( $$data_ref{$GDATA}{$ds_name_full} ) )
                      {
                        push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name_full );
                        $$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GNAME} =
                          $ds_name_tol;
                        $$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GNAME} =
                          "time";
                      }
                    push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GORG}},
                          $val );
                    push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDX}{$GORG}},
                          $time );
                  }
              }
          }
        #........................
        #...DONE: delta blocks...
        #........................
      }
    #........................
    #...look at every line...
    #........................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_ares
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) stuff the specific lines into %data
#............................................................................
sub read_file_ares
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name_full, # corresponding arrays with this name are diffed
       $ds_name_tol,  # dataset name used for diff tolerances
       $ierr,         # error return value
       $key,          # key/val of token
       $line,         # line of file
       @tokens,       # split on whitespace of $line
       $type,         # used in skipping and naming dataset
       $val,          # key/val of token
      );
    $ierr = 0;
    $title = sprintf( "cycle %5d", 0 );
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #........................
    #...look at every line...
    #........................
    while( $line=<FILE> )
      {
        #..................
        #... result line...
        #..................
        if( $line =~
            /
            ^\s*                    # start with possible whitespace
            (\S+)\s+                # test
            (DIFF|FAILED|PASSED)\s+ # test results
            (P-[0-9]+),\s+          # P field,
            (D-[0-9]+),\s+          # D field,
            (F-[0-9]+)\s+           # F field
            $/x )
          {
            $type = $1;
            @tokens = ($2, $3, $4, $5);
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $type !~ /^($$cmd_ref{$GDS})$/ )
              {
                next;
              }
            $ds_name_tol  = "$type";
            $ds_name_full = "$ds_name_tol";
            if( ! defined( $$data_ref{$GDATA}{$ds_name_full}) )
              {
                push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name_full );
                $$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GNAME} =
                  $ds_name_tol;
              }
            push( @{$$data_ref{$GDATA}{$ds_name_full}{$GCOORDY}{$GORG}},
                  @tokens );
          }
      }
    #........................
    #...look at every line...
    #........................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_cts
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) get to start of data
#... 2) repeat:
#... 2.1) get dataset name
#... 2.2) get coordinate name(s)
#... 2.3) stuff values into dataset
#............................................................................
sub read_file_cts
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $coord,        # general  name of a coordinate
       $coord_name,   # specific name of a coordinate
       @coords,       # the general names
       $copy,         # copy number for unique ds_name
       $ds_name,      # dataset name
       $ds_name_orig, # original dataset name (if needed to create new one)
       $i,            # loop variable
       $ierr,         # error return value
       $line,         # line of file
       $line_num,     # line of file
       @tokens,       # split on whitespace of $line
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #...............................
    #...get passed first cts line...
    #...............................
    $line=<FILE>;
    #......................
    #...read in datasets...
    #......................
    $ds_name = "";
    $line_num = 1;
    while( $line=<FILE> )
      {
        $line_num++;
        #................
        #...blank line...
        #................
        if( $line !~ /\S/ )
          {
            next;
          }
        #.................
        #...new dataset...
        #.................
        if( $line =~
            /
            ^\#\s*Dataset\s*Name:\s* # "Dataset Name:"
            (\S.*?)                  # dataset name
            \s+$                     # at least return whitespace
            /x )
          {
            $ds_name = $1;
            undef( @coords );
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_name !~ /^($$cmd_ref{$GDS})$/ )
              {
                $ds_name = "";
                next;
              }
            #.....................
            #...initialize data...
            #.....................
            if( $ds_name =~ /\S/ )
              {
                #..........................................
                #...create unique dataset name if needed...
                #..........................................
                if( defined( $$data_ref{$GDATA}{$ds_name} ) )
                  {
                    $copy = 1;
                    $ds_name_orig = $ds_name;
                    while( 1 == 1 )
                      {
                        $ds_name = sprintf( "%s $GCOPY_FORMAT",
                                            $ds_name_orig, $copy );
                        if( !defined( $$data_ref{$GDATA}{$ds_name} ) )
                          {
                            last;
                          }
                        $copy++;
                      }
                  }
                push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
              }
          }
        #.....................
        #...coordinate name...
        #.....................
        elsif( $ds_name =~ /\S/ &&
               $line =~
               /
               ^\#\s*Coord\s*Name\s* # "Coord Name:"
               (\S+):\               # coord name general:
               (\S.*?)               # coord name specific
               \s+$                  # at least return whitespace
               /x )
          {
            $coord      = $1;
            $coord_name = $2;
            push( @coords, $coord );
            $$data_ref{$GDATA}{$ds_name}{$coord}{$GNAME} = $coord_name;
          }
        #..........................................
        #...skip other lines starting with pound...
        #..........................................
        elsif( $line =~ /^\#/ )
          {
          }
        #...............
        #...data line...
        #...............
        elsif( $ds_name =~ /\S/ )
          {
            $line =~ s/^\s*(\S.*?)\s*$/$1/;
            @tokens = split( /\s+/, $line );
            if( $#tokens != $#coords )
              {
                $ierr = 1;
                &print_error( "Mismatch in number of values/coordinates",
                              "Coords: ".join(', ', @coords),
                              "Values: ".join(', ', @tokens),
                              "File: [$file_name:$line_num]",
                              $ierr );
                return( $ierr );
              }
            for( $i = 0; $i <= $#tokens; $i++ )
              {
                push( @{$$data_ref{$GDATA}{$ds_name}{$coords[$i]}{$GORG}},
                      $tokens[$i] );
              }
          }
        #.....................
        #...DONE: data line...
        #.....................
      }
    #............................
    #...DONE: read in datasets...
    #............................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_keyword
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) find a keyword line and stuff it into data.
#............................................................................
sub read_file_keyword
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,  # dataset name
       $ierr,     # error return value
       $line,     # line of file
       $val,      # value
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #........................
    #...look at every line...
    #........................
    while( $line=<FILE> )
      {
        if( $line =~
            /
            ^\s*(\S+)    # dataset name
            \s*=\s*      # =
            (\S+)\s*$    # value
            /x )
          {
            $ds_name = $1;
            ($val = $2) =~ s/(\d)([+-]\d)/$1e$2/;
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_name !~ /^($$cmd_ref{$GDS})$/ )
              {
                next;
              }
            #...........................................
            #...init data if this is first occurrance...
            #...........................................
            if( ! defined($$data_ref{$GDATA}{$ds_name}) )
              {
                push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
                $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name;
              }
            push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                  $val );
          }
      }
    #.................................
    #...DONE: push values onto data...
    #.................................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_oxy
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) get to start of data (2 line matches)
#... 2) get dataset name
#... 3) stuff values into dataset
#............................................................................
sub read_file_oxy
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,      # dataset name
       $ds_name_orig, # original dataset name (if needed to create new one)
       $copy,         # copy number for unique ds_name
       $i,            # loop variable
       $ierr,         # error return value
       $line,         # line of file
       $number,       # a number (for incrementing dataset name)
       $skip_ds,      # if skipping this dataset
       @tokens,       # split on whitespace of $line

       @lines,        # lines of the file after being read by translator
      );
    &read_file_oxy_orig_translator( $file_name, \@lines );
    $ierr = 0;
    #.......................
    #...process each line...
    #.......................
    $ds_name = "unknown";
    $ds_name_orig = "unknown";
    $line_num = 0;
    $skip_ds = 0;
    foreach $line ( @lines )
      {
        $line =~ s/^\s*(.*?)\s*$/$1/;
        #..............................
        #...data line (or something)...
        #..............................
        if( $line =~ /^[^\#]/ )
          {
            #........................................
            #...skip line if skipping this dataset...
            #........................................
            if( $skip_ds )
              {
                next;
              }
            @tokens = split( /\s+/, $line );
            grep( s/(\d)([+-]\d)/$1e$2/, @tokens );
            #......................
            #...simple data line...
            #......................
            if( $#tokens == 1 )
              {
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      $tokens[0] );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $tokens[1] );
              }
            #.......................................................
            #...character data line - just put all data on 1 line...
            #.......................................................
            else
              {
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      "" );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $line );
              }
            if( ! defined($$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME}) )
              {
                push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
                $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = "unknown";
                $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = "unknown";
              }
          }
        #..............................
        #...data line (or something)...
        #..............................
        #..................
        #...new data set...
        #..................
        elsif( $line =~
               /
               ^\#\s*    # pound + whitespace
               (\S.*)    # data set name
               $/x )
          {
            $ds_name      = $1;
            $ds_name_orig = $ds_name;
            #................................
            #...see if should skip this ds...
            #................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_name_orig !~ /^($$cmd_ref{$GDS})$/ )
              {
                $skip_ds = 1;
                next;
              }
            else
              {
                $skip_ds = 0;
              }
            #..........................................
            #...create unique dataset name if needed...
            #..........................................
            if( defined( $$data_ref{$GDATA}{$ds_name} ) )
              {
                $copy = 1;
                while( 1 == 1 )
                  {
                    $ds_name = sprintf( "%s $GCOPY_FORMAT",
                                        $ds_name_orig, $copy );
                    if( !defined( $$data_ref{$GDATA}{$ds_name} ) )
                      {
                        last;
                      }
                    $copy++;
                  }
              }
            push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
            $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = $GCOORDX;
            $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name_orig;
          }
      }
    #.............................
    #...DONE: process each line...
    #.............................
    undef( @lines );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_plot_output
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) first line has dataset names
#... 2) stuff lines into data
#............................................................................
sub read_file_plot_output
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $done,     # if done
       $done1,    # another done flag
       @ds_names, # dataset names
       $ds_name,  # dataset name
       $i,        # loop variable
       $ierr,     # error return value
       $line,     # line of file
       $line_num, # line number
       $number,   # a number (for incrementing dataset name)
       $skip_ds,  # if skipping this dataset
       @tokens,   # split on whitespace of $line
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    # get to block
    while( $line=<FILE> ){
        if( $line =~ /#\s*\[\d+\]\s*(\S+.*?)\s*$/ ){
            $ds_start = $1;
            $ds_start =~ s/\s+/_/g;
            $line = <FILE>;
            if( ! defined( $line ) || $line !~ /\S/ ){
                next;
            }
            if( defined( $$cmd_ref{$GDS} ) &&
                $$cmd_ref{$GDS} eq "track" ){
                if( $ds_start =~
                    /Resources:memory|
                     cc\/s\/p|
                     dmp_write_time|
                     lost_cycles|
                     procmon_|
                     secs|
                     sumRSS_GB|
                     sumcpu|
                     sumwallhr
                     /x ){
                    next;
                }
            }
            elsif( defined( $$cmd_ref{$GDS} ) &&
                $ds_start !~ /^($$cmd_ref{$GDS})/ ) {
                next;
            }
            print "Datasets: [$file_name:$ds_start]\n";
            $line =~ s/^\s+//;
            @col_names = split( /\s+/, $line );
            # prepend ${ds_start}:: to column headers to for ds_names
            @ds_names = @col_names;
            grep( s/^/${ds_start}::/, @ds_names );
            # for now, only do cycle|time|var1|var2|... blocks
            if( $#col_names >= 2 && $col_names[0] eq "cycle" && $col_names[1] eq "time" ){
                # register names
                for( $i = 2; $i <= $#col_names; $i++ ){
                    push( @{$$data_ref{$GDATASET_NAMES}}, $ds_names[$i] );
                    $$data_ref{$GDATA}{$ds_names[$i]}{$GCOORDX}{$GNAME} = $col_names[1];
                    $$data_ref{$GDATA}{$ds_names[$i]}{$GCOORDY}{$GNAME} = $col_names[$i];
                }
                # stuff in this block
                $time_old = -1e99;
                while( $line=<FILE> ){
                    if( $line !~ /\S/ ){
                        last;
                    }
                    $line =~ s/^\s+//;
                    @vals = split( /\s+/, $line );
                    # only add to array if time is increasing
                    if( $vals[1] > $time_old ){
                        $time_old = $vals[1];
                        for( $i = 2; $i <= $#col_names; $i++ ){
                            push( @{$$data_ref{$GDATA}{$ds_names[$i]}{$GCOORDX}{$GORG}}, $vals[1] );
                            push( @{$$data_ref{$GDATA}{$ds_names[$i]}{$GCOORDY}{$GORG}}, $vals[$i] );
                        }
                    }
                }
            }
        }
    }
    close( FILE );
    return( $ierr );
  }
#.........................................................................
#...this is taken from the original translator (see file in obsolete)  ...
#...It had some statements                                             ...
#...that were confusing (the set of 11 perl regexp modifiers) and some ...
#...bugs (at least I think they are bugs).  It basically takes input in...
#...a poorly spec'd form and produces output in a poorly spec'd form.  ...
#...My comments will have a "cts" in the comment line                  ...
#...This routine fills an array of the lines of the routine and        ...
#...returns is - might have to be changed if that takes too much mem   ...
#...Changes:                                                           ...
#...  - Added indentation to translator to turn into subroutine        ...
#...  - Added "my" variables to insulate (and changed local to my)     ...
#...  - Added $Count as argument to FixCounts to avoid global variable ...
#...  - Added $zero arg to FixCounts - effectively a global var        ...
#...  - Initialize $zero to 0
#...  - Changed print to push onto array                               ...
#.........................................................................
sub read_file_oxy_orig_translator
  {
    my(
       $file_name,
       $lines_ref,
      ) = @_;
    my(
       $StartString1,
       $StartString2,
       $EndString,
       $Count,
       $zero,
       $InFile,
       $Matches,
       $Line,
      );
    #  This script formats the output from a run of one of the TestSuite input
    #  decks into the standard format for comparision with a standard.
    
    $StartString1 = "Final state ASCII diagnostic dump start";
    $StartString2 = "Scalar stop";
    $EndString = "print. Final";
    $Count = 0;
    #...cts...
    $zero = 0;
    
    #main( Infile, Outfile )
    # cts comment out usage, rename $InFile, and comment out outfile
    # cts if( $#ARGV < 1 ) {
    # cts     print "Usage: $0 <InputFile> <OutputFile>\n";
    # cts     exit(-1);
    # cts }
    # cts $InFile = shift(@ARGV);
    $InFile = $file_name;
    open( IN, $InFile ) || die "Unable to open input file: $InFile\n";
    # cts $OutFile = shift(@ARGV);
    # cts open( OUT, ">".$OutFile ) || die "Unable to open output file: $OutFile\n";
    
    $Matches = 0;
    while( $Matches < 3 && ($Line = <IN>) ) {
    #      $Line = <IN>;
        chop( $Line );
        if( $Matches < 1 ) {
    	if( $Line =~ /$StartString1/ ) {
                $Matches++;
    	}
    	next; # Don't print this line
        }
        elsif( $Matches < 2 ) {
    	if( $Line =~ /$StartString2/ ) {
                $Matches++;
    	}
    	next; # Don't print this line
        }
        if( $Line =~ /$EndString/ ) {
    	$Matches++;
    	last; # Done
        }
        if( $Line eq "" ) {
    	next; # Don't print this line
        }
        $Line =~ s/mype=    0 //;
        $Line =~ s/[0-9]{1,} {0,}([a-z]{1,}) {1,}nul=.{0,}/ $1/;
        $Line =~ s/^[ \t]*[0-9]*:[ \t]*//;
        $Line =~ s/^[0-9]{2,}//;
        $Line =~ s/[ ]{1,}/\n/g;
        $Line =~ s/([a-z]{1,})/#  $1/;
        $Line =~ s/([+-][0-9][0-9])([+-])/$1\n$2/g;
        $Line =~ s/\n(\.0{1,})/\n0$1/g;
        $Line =~ s/^(\.0{1,})/0$1/g;
        $Line =~ s/\n{2,}/\n/g;
        $Line =~ s/^\n//g;
        
        # cts replace print with push
        # cts print OUT &FixCounts( $Line );
        push( @$lines_ref, split( /\n/, &FixCounts( $Line, \$Count, \$zero ) ) );
    }
  }
sub FixCounts {
   #local( $In ) = @_;
   #local( $i, @lines, $line );
   #...cts...
   my(
      $In,
      $Count_ref,
      $zero_ref,
      ) = @_;
   my(
      $i,
      $item,
      @lines,
      $line,
     );

   # cts
   $line = "";
   @lines = split(/\n/,$In);
   foreach $item (@lines) {
      if( $item =~ /#  [a-z]{1,}/ ) {
         $$Count_ref = 0;
      }
      else {
         ${$Count_ref}++;
         if( ($item =~ /7\.777777000E\+83/i) || 
             ($item =~ /1\.000000000E-99/i) ) {
            $item = "";  # Invalid, skip
         }
         else {
            if( $item =~ /0\.000000000E\+00/i ) {  # should match 0 instead
               if( $$zero_ref == 0 ) {  # First 0
                  $item =~ s/^/ $$Count_ref    /; # Add line count
                  $$zero_ref = 1;
               }
               else {  # Repeat 0, skip
                  $item = "";
               }
            }
            else { # Valid number of interest
               $item =~ s/^/ $$Count_ref    /; # Add line count
               $$zero_ref = 0; # non-zero
            }
         }
      }
      if( $item ne "" ) {
         $line .= "$item\n"; # Reaccumulate
      }
   }

   $line;  # Return value
} # end FixCounts

#............................................................................
#...Name
#...====
#... read_file_oxy_cts
#...
#...Purpose
#...=======
#... NOTE: This was an attempt to write my own parser.
#...       For now, read_file_oxy directly incorporates logic from their
#...       old translator.
#...
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) get to start of data (2 line matches)
#... 2) get dataset name
#... 3) stuff values into dataset
#............................................................................
sub read_file_oxy_cts
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $copy,         # copy number for unique ds_name
       $ds_name,      # dataset name
       $ds_name_orig, # original dataset name (if needed to create new one)
       $form,         # different forms of new data line
       $found,        # if a tag was found
       $ierr,         # error return value
       $index,        # index for value
       $line,         # line of file
       $pushed_token, # if a token was pushed onto the data
       $skip_zero,    # if skipping the following 0 values
       $skip_val,     # if skipping the printing of the value
       $tag_start,    # starting tag
       $tag_stop,     # stopping tag
       $tag_ds_start, # starting place where ds data is
       @tokens,       # split on whitespace of $line
       $token,        # single token
       $token_last,   # the last token
       %skip_vals,   # special values to be skipped
       $val_type,     # for replacing values of %skip_vals
      );
    $ierr = 0;
    $tag_start = $GOXY_TAG_START;
    $tag_stop = "Final state ASCII diagnostic dump stop";
    $tag_ds_start = "Scalar stop";
    $skip_vals{skip1} = "7.777777000E+83";
    $skip_vals{skip2} = "1.000000000E-99";
    $index     = 0;
    $skip_val  = 0;
    $skip_zero = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #...............................
    #...get to lines of interest....
    #...............................
    $found = 0;
    while( $line=<FILE> )
      {
        if( $line =~ /$tag_start/ )
          {
            while( $line=<FILE> )
              {
                if( $line =~ /$tag_ds_start/ )
                  {
                    $found = 1;
                    last;
                  }
              }
            last;
          }
      }
    if( ! $found )
      {
        $ierr = 1;
        &print_error( "Incorrect format - cannot keyword lines",
                      "The following lines are needed:",
                      "  [$tag_start]",
                      "  [$tag_ds_start]",
                      $ierr );
        return( $ierr );
      }
    #......................
    #...read in datasets...
    #......................
    $ds_name = "";
    $pushed_token = 0;
    $found = 0;
    while( $line=<FILE> )
      {
        #.....................
        #...end of all data...
        #.....................
        if( $line =~ /$tag_stop/ )
          {
            $found = 1;
            #................................
            #...finish last 0 if necessary...
            #................................
            if( $skip_zero != 0 && $pushed_token == 0 )
              {
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      $index );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $token_last );
              }
            last;
          }
        #.................
        #...new dataset...
        #.................
        if( $line =~
            /
            ^\s*[0-9]+\s+      # whitespace and some digits
            ([a-z]+)\s+        # dataset name
            (\S+)=\s+(\S+)\s+  # nul value
            (\S+)=\s+(\S+)     # out_of_range value
            /x ||
            $line =~
            /
            ^mype=\s*[0-9]+\s+[0-9]+\s* # mype= <num> <num>
            ([a-z]+)                    # <name>
            (.*)$                       # vals
            /x)
          {
            #................................
            #...finish last 0 if necessary...
            #................................
            if( $skip_zero != 0 && $pushed_token == 0 )
              {
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      $index );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $token_last );
              }
            $ds_name = $1;
            if( defined( $3 ) )
              {
                $form = 1;
                #...no more data on line...
                $line = "";
              }
            else
              {
                $form = 2;
                #...data on line...
                $line = $2;
              }
            $index = 0;
            $skip_zero = 0;
            $pushed_token = 0;
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_name !~ /^($$cmd_ref{$GDS})$/ )
              {
                $ds_name = "";
                next;
              }
            #.....................
            #...initialize data...
            #.....................
            if( $ds_name =~ /\S/ )
              {
                #..........................................
                #...create unique dataset name if needed...
                #..........................................
                if( defined( $$data_ref{$GDATA}{$ds_name} ) )
                  {
                    $copy = 1;
                    $ds_name_orig = $ds_name;
                    while( 1 == 1 )
                      {
                        $ds_name = sprintf( "%s $GCOPY_FORMAT",
                                            $ds_name_orig, $copy );
                        if( !defined( $$data_ref{$GDATA}{$ds_name} ) )
                          {
                            last;
                          }
                        $copy++;
                      }
                  }
                push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
                $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = $GCOORDX;
                $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name;
              }
          }
        #...............
        #...data line...
        #...............
        if( $ds_name =~ /\S/ && $line =~ /\S/ )
          {
            $line =~ s/^\s*(\S.*?)\s*$/$1/;
            #...remove index if there is one...
            $line =~ s/\s*[0-9]+:\s*//;
            @tokens = split( /\s+/, $line );
            grep( s/(\d)([+-]\d)/$1e$2/, @tokens );
            #........................................
            #...replace values with special values...
            #........................................
            foreach $token ( @tokens )
              {
                $skip_val = 0;
                #...only use first of consecutive 0s...
                if( $token == 0 )
                  {
                    #...use first 0...
                    if( $skip_zero == 0 )
                      {
                        $skip_zero = 1;
                      }
                    #...skip other 0s...
                    else
                      {
                        $skip_val = 1;
                      }
                  }
                #...skip special values...
                else
                  {
                    #................................
                    #...finish last 0 if necessary...
                    #................................
                    if( $skip_zero != 0 && $pushed_token == 0 )
                      {
                        push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                              $index );
                        push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                              $token_last );
                      }
                    $skip_zero = 0;
                    foreach $val_type ( keys %skip_vals )
                      {
                        if( $token eq "$skip_vals{$val_type}" )
                          {
                            $skip_val = 1;
                            last;
                          }
                      }
                  }
                $index++;
                #...push if not skipping...
                if( ! $skip_val )
                  {
                    $pushed_token = 1;
                    push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                          $index );
                    push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                          $token );
                  }
                else
                  {
                    $pushed_token = 0;
                  }
                $token_last = $token;
              }
            #..............................................
            #...DONE: replace values with special values...
            #..............................................
          }
        #.....................
        #...DONE: data line...
        #.....................
      }
    #............................
    #...DONE: read in datasets...
    #............................
    if( ! $found )
      {
        $ierr = 1;
        &print_error( "Incorrect format - cannot find stop tag",
                      "[$tag_stop]",
                      $ierr );
        return( $ierr );
      }
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_pop
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) find a pop line and stuff it into data.
#............................................................................
sub read_file_pop
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,      # dataset name
       $ierr,         # error return value
       %keywords,    # %keywords, %var, %regexps
       $keyword_y,    # keyword label for variable
       $keyword_time, # keyword label for variable
       $keyword_ds,   # keyword for ds
       %regexps,     # a regexp
       $line,         # line of file
       @tokens,       # split on whitespace of $line
       %vars          # current variable name hash
      );
    $ierr         = 0;
    $keywords{ds}   = 'card';
    $keywords{time} = 'time';
    $keywords{x}    = 'x var:';
    $keywords{y}    = 'y var:';
    $regexps{ds}   = 'card';
    $regexps{time} = 'time';
    $regexps{x}    = 'x\ var:';
    $regexps{y}    = 'y\ var:';
    #....................
    #...default values...
    #....................
    $ds_name_base = "unknown";
    $vars{ds}     = 0;
    $vars{time}   = 0;
    $vars{x}      = $GCOORDX;
    $vars{y}      = $GCOORDY;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #........................
    #...look at every line...
    #........................
    while( $line=<FILE> )
      {
        #.......................................
        #...data - have before other id match...
        #.......................................
        if( $line =~
            /
            ^\ $regexps{ds}\s*[0-9]+\.\s* # new id
            write\s*$ # ends in write
            /x )
          {
            $ds_name = $ds_name_base;
            #............................................
            #...skip if not interested in this dataset...
            #............................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_name !~ /^($$cmd_ref{$GDS})$/ )
              {
                next;
              }
            $ds_name = "$keywords{ds}=$vars{ds},$ds_name,time=$vars{time}";
            #...............................
            #...slurp this dataset values...
            #...............................
            while( $line=<FILE> )
              {
                $line =~ s/^\s*//;
                $line =~ s/\s*$//;
                if( $line !~ /\S/ )
                  {
                    next;
                  }
                if( $line =~ /^curve\s+/ )
                  {
                    last;
                  }
                @tokens = split( /\s+/, $line );
                grep( s/(\d)([+-]\d)/$1e$2/, @tokens );
                #......................
                #...simple data line...
                #......................
                if( $#tokens == 1 )
                  {
                    push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                          $tokens[0] );
                    push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                          $tokens[1] );
                  }
                #.......................................................
                #...character data line - just put all data on 1 line...
                #.......................................................
                else
                  {
                    push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                          "" );
                    push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                          $line );
                  }
                if( ! defined($$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME}) )
                  {
                    push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
                    $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = $vars{x};
                    $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $vars{y};
                  }
              }
            #.....................................
            #...DONE: slurp this dataset values...
            #.....................................
            next;
          }
        #.........................
        #...new ID - reset vals...
        #.........................
        if( $line =~
            /
            ^\ $regexps{ds}\s+ # starts with this keyword
            ([0-9]+)            # dataset ID
            \.\s*               # dot and spaces
            (\S+.*?)            # dataset name base
            \s*$                # any whitespace and end
            /x )
          {
            $vars{ds}      = sprintf( "%4d", $1 );
            $ds_name_base  = $2;
            $vars{time}    = 0;
            $vars{x}       = $GCOORDX;
            $vars{y}       = $GCOORDY;
            next;
          }
        #..........................
        #...* <variable> <value>...
        #..........................
        if( $line =~
            /
            ^\s+\*\s+ # start with star
            $regexps{time} # variable
            \s*       # whitespace
            (\S.*?)   # value
            \s+$      # whitespace and end
            /x )
          {
            $vars{time} = $1;
            next;
          }
        if( $line =~
            /
            ^\s+\*\s+ # start with star
            $regexps{x} # variable
            \s*       # whitespace
            (\S.*?)   # value
            \s+$      # whitespace and end
            /x )
          {
            $vars{x} = $1;
            next;
          }
        if( $line =~
            /
            ^\s+\*\s+ # start with star
            $regexps{y} # variable
            \s*       # whitespace
            (\S.*?)   # value
            \s+$      # whitespace and end
            /x )
          {
            $vars{y} = $1;
            next;
          }
      }
    #.................................
    #...DONE: push values onto data...
    #.................................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_table
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) first line has dataset names
#... 2) stuff lines into data
#............................................................................
sub read_file_table
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       @ds_names, # dataset names
       $ds_name,  # dataset name
       $i,        # loop variable
       $ierr,     # error return value
       $line,     # line of file
       $line_num, # line number
       $number,   # a number (for incrementing dataset name)
       $skip_ds,  # if skipping this dataset
       @tokens,   # split on whitespace of $line
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #..............................
    #...get to/parse header line...
    #..............................
    while( $line=<FILE> )
      {
        #......................
        #...skip: blank line...
        #......................
        if( $line !~ /\S/ )
          {
            next;
          }
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        @ds_names = split( /\s+/, $line );
        last;
      }
    #..................
    #...init ds vals...
    #..................
    foreach $ds_name ( @ds_names )
      {
        if( defined( $$cmd_ref{$GDS} ) &&
            $ds_name !~ /^($$cmd_ref{$GDS})$/ )
          {
            next;
          }
        push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
        $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name;
      }
    #...........................
    #...push values onto data...
    #...........................
    while( $line=<FILE> )
      {
        if( $line !~ /\S/ )
          {
            next;
          }
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        @tokens = split( /\s+/, $line );
        grep( s/(\d)([+-]\d)/$1e$2/, @tokens );
        for( $i = 0; $i <= $#ds_names; $i++ )
          {
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_names[$i] !~ /^($$cmd_ref{$GDS})$/ )
              {
                next;
              }
            push( @{$$data_ref{$GDATA}{$ds_names[$i]}{$GCOORDY}{$GORG}},
                  $tokens[$i] );
          }
      }
    #.................................
    #...DONE: push values onto data...
    #.................................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_token
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) find a keyword line and stuff it into data.
#............................................................................
sub read_file_token
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,  # dataset name
       $ierr,     # error return value
       $line,     # line of file
       @tokens,   # split on whitespace of $line
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #........................
    #...look at every line...
    #........................
    $ds_name = "token";
    push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
    $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name;
    while( $line=<FILE> )
      {
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        @tokens = split( /\s+/, $line );
        grep( s/(\d)([+-]\d)/$1e$2/, @tokens );
        push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}}, @tokens );
      }
    #.................................
    #...DONE: push values onto data...
    #.................................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_tracer
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) foreach line
#... 1.1) new dataset line - reset stuff
#... 1.2) data line of dataset - stuff value into data
#............................................................................
sub read_file_tracer
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,        # dataset name
       @fields,         # fields (columns)
       $ds_name_root,   # root ds name
       $i,              # loop variable
       $ierr,           # error return value
       $line,           # line of file
       $line_num,       # line number in file
       $number,         # a number (for incrementing dataset name)
       $particle,       # particle number
       $particle_field, # the field number
       @seen,           # if seen this particle
       $time,           # the time value
       $time_field,     # the field number
       @tokens,         # items to push into onto data arrays
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #.......................
    #...process each line...
    #.......................
    $ds_name = "unknown";
    $line_num = 0;
    $particle_field = -1;
    $time_field = -1;
    while( defined( $line=<FILE> ) )
      {
        $line_num++;
        #...remove leading/trailing whitespace and comment lines
        $line =~ s/^\s*(.*?)\s*$/$1/;
        $line =~ s/\s*#.*//;
        #..............................
        #...data line (or something)...
        #..............................
        if( $line =~ /\S/ )
          {
            @tokens = split( /\s*,\s*/, $line );
            #...field names
            if( $#fields < 0 )
              {
                @fields = @tokens;
                for( $i = 0; $i <= $#fields; $i++ )
                  {
                    if( $fields[$i] eq "particle" )
                      {
                        $particle_field = $i;
                      }
                    if( $fields[$i] eq "time" )
                      {
                        $time_field = $i;
                      }
                  }
                if( $particle_field == -1 || $time_field == -1 )
                  {
                    $ierr = 1;
                    &print_error( "Cannot find particle and/or time fields from header: $file_name:$line_num",
                                  "[$line]",
                                  $ierr );
                    return( $ierr );
                  }
                next;
              }
            #...data
            $time = $tokens[$time_field];
            $particle = $tokens[$particle_field];
            $ds_name_root = "p_${particle}_";
            #...stuff onto name array
            if( ! defined($seen[$particle]) )
              {
                $seen[$particle] = 1;
                for( $i = 0; $i <= $#tokens; $i++ )
                  {
                    #...skip fields
                    if( $i == $time_field ||
                        $i == $particle_field ||
                        ( defined( $$cmd_ref{$GDS} ) &&
                          $fields[$i] !~ /^($$cmd_ref{$GDS})$/ ))
                      {
                        next;
                      }
                    $ds_name = "${ds_name_root}$fields[$i]";
                    push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
                    $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = "time";
                    $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = "$ds_name";
                  }
              }
            #...stuff values
            for( $i = 0; $i <= $#tokens; $i++ )
              {
                #...skip fields
                if( $i == $time_field ||
                    $i == $particle_field ||
                    ( defined( $$cmd_ref{$GDS} ) &&
                      $fields[$i] !~ /^($$cmd_ref{$GDS})$/ ))
                  {
                    next;
                  }
                $ds_name = "${ds_name_root}$fields[$i]";
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      $tokens[$time_field] );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $tokens[$i] );
              }
          } # if data line
      } # each line
    #.............................
    #...DONE: process each line...
    #.............................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_link
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#... Almost identical to read_file_token, except with an extra fix
#... to handle formatted numbers that run into each other, e.g.:
#...      0.0000000E+00-1.7347235E-18
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) find a keyword line and stuff it into data.
#............................................................................
sub read_file_link
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,  # dataset name
       $ierr,     # error return value
       $line,     # line of file
       @tokens,   # split on whitespace of $line
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #........................
    #...look at every line...
    #........................
    $ds_name = "token";
    push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
    $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name;
    while( $line=<FILE> )
      {
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        # add an extra space to separate any pair of Enn.n numbers
        # that have run together
        $line =~ s/(?<=\d)(-\d\.\d+E[+-]\d\d)/ $1/g;
        @tokens = split( /\s+/, $line );
        push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}}, @tokens );
      }
    #.................................
    #...DONE: push values onto data...
    #.................................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... read_file_xy
#...
#...Purpose
#...=======
#... Read a file and stuff it into $data_ref
#...
#...Arguments
#...=========
#... $file_name   Intent: in
#...              Perl type: scalar
#...              File name to read in
#...
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line options
#...
#... $data_ref    Intent: out
#...              Perl type: reference to hash
#...              See read_file for format
#...
#... $ierr        Intent: out
#...              Perl type: reference to hash
#...              Return value (non-0 is error)
#...
#...Program Flow
#...============
#... 1) foreach line
#... 1.1) new dataset line - reset stuff
#... 1.2) data line of dataset - stuff value into data
#............................................................................
sub read_file_xy
  {
    my(
       $file_name,
       $cmd_ref,
       $data_ref,
      ) = @_;
    my(
       $ds_name,      # dataset name
       $ds_name_orig, # original dataset name (if needed to create new one)
       $copy,         # copy number for unique ds_name
       $i,            # loop variable
       $ierr,         # error return value
       $line,         # line of file
       $number,       # a number (for incrementing dataset name)
       $skip_ds,      # if skipping this dataset
       @tokens,       # items to push into onto data arrays
      );
    $ierr = 0;
    #...............
    #...open file...
    #...............
    if( ! open( FILE, "$file_name" ) )
      {
        $ierr = 1;
        &print_error( "Cannot open data file [$file_name].",
                      $ierr );
        return( $ierr );
      }
    #.......................
    #...process each line...
    #.......................
    $ds_name = "unknown";
    $ds_name_orig = "unknown";
    $line_num = 0;
    $skip_ds = 0;
    while( defined( $line=<FILE> ) )
      {
        $line =~ s/^\s*(.*?)\s*$/$1/;
        #..............................
        #...data line (or something)...
        #..............................
        if( $line =~ /^[^\#]/ )
          {
            #........................................
            #...skip line if skipping this dataset...
            #........................................
            if( $skip_ds )
              {
                next;
              }
            #................................................................
            #...if line consists of at least 2 whitespace separated tokens...
            #................................................................
            if( $line =~ /^(\S+)\s+(\S.*)$/ )
              {
                @tokens = ($1, $2);
                #..........................................
                #...fix broken numbers (1+123 -> 1e+123)...
                #..........................................
                grep( s/(\d)([+-]\d)/$1e$2/, @tokens );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      $tokens[0] );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $tokens[1] );
              }
            #................................................
            #...otherwise, just push whole line to GCOORDY...
            #................................................
            else
              {
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GORG}},
                      "undef" );
                push( @{$$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GORG}},
                      $line );
              }
            #..................
            #...init dataset...
            #..................
            if( ! defined($$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME}) )
              {
                push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
                $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = "unknown";
                $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = "unknown";
              }
          }
        #....................................
        #...DONE: data line (or something)...
        #....................................
        #..................
        #...new data set...
        #..................
        elsif( $line =~
               /
               ^\#\s*    # pound + whitespace
               (\S.*)    # data set name
               $/x )
          {
            $ds_name      = $1;
            $ds_name_orig = $ds_name;
            #................................
            #...see if should skip this ds...
            #................................
            if( defined( $$cmd_ref{$GDS} ) &&
                $ds_name_orig !~ /^($$cmd_ref{$GDS})$/ )
              {
                $skip_ds = 1;
                next;
              }
            else
              {
                $skip_ds = 0;
              }
            #..........................................
            #...create unique dataset name if needed...
            #..........................................
            if( defined( $$data_ref{$GDATA}{$ds_name} ) )
              {
                $copy = 1;
                while( 1 == 1 )
                  {
                    $ds_name = sprintf( "%s $GCOPY_FORMAT",
                                        $ds_name_orig, $copy );
                    if( !defined( $$data_ref{$GDATA}{$ds_name} ) )
                      {
                        last;
                      }
                    $copy++;
                  }
              }
            push( @{$$data_ref{$GDATASET_NAMES}}, $ds_name );
            $$data_ref{$GDATA}{$ds_name}{$GCOORDX}{$GNAME} = $GCOORDX;
            $$data_ref{$GDATA}{$ds_name}{$GCOORDY}{$GNAME} = $ds_name_orig;
          }
      }
    #.............................
    #...DONE: process each line...
    #.............................
    close( FILE );
    return( $ierr );
  }
#............................................................................
#...Name
#...====
#... run_gnuplot
#...
#...Purpose
#...=======
#... Runs gnuplot given %gnuplot_info (creates plots)
#...
#...Arguments
#...=========
#... $gnuplot_info_ref Intent: inout
#...                   Perl type: reference to hash
#...                   Contains info gnuplot will use to plot.
#...                   Filled by print_gnuplot_data routine.
#...
#...Program Flow
#...============
#... 1) Compute differences
#............................................................................
sub run_gnuplot
  {
    my(
       $gnuplot_info_ref
      ) = @_;
    my(
       $arch, # machine on (used with which_exec)
       $dtype, # the data type in %gnuplot_info
       %dtypes, # the particular dtypes defined for a variable (plots/page)
       $fix_file, # fix for doing ps2pdf with landscape
       $gnuplot, # where gnuplot is
       $gnuplot_cmd, # gnuplot commands
       $landscape, # filled with the landscape option
       $num_plots, # number of plots in a page
       $page, # current page number
       $plot_number, # plot number of current page
       $origin, # gnuplot origin for plot
       $output, # output from shell command
       $ps2pdf, # where ps2pdf is
       $size, # gnuplot size for plot
       $source, # the source in %gnuplot_info
       $variable, # the variable in %gnuplot_info
       $using_x, # x column
       $using_y, # y column
       $xlabel, # label for coord
       $ylabel, # label for coord
       $title, # title of plot
      );
    #.............................
    #...exit now if not defined...
    #.............................
    if( ! %$gnuplot_info_ref )
      {
        return;
      }
    #..........
    #...init...
    #..........
    $gnuplot_cmd = "";
    $page = 0;
    if( ! defined $$gnuplot_info_ref{cmd_file} )
      {
        $$gnuplot_info_ref{cmd_file} = "cts_diff.cmd";
      }
    if( ! defined $$gnuplot_info_ref{ps_file} )
      {
        $$gnuplot_info_ref{ps_file} = "cts_diff.ps";
      }
    if( ! defined $$gnuplot_info_ref{pdf_file} )
      {
        $$gnuplot_info_ref{pdf_file} = "cts_diff.pdf";
      }
    if( defined( $$gnuplot_info_ref{orientation} ) )
      {
        $landscape = $$gnuplot_info_ref{orientation};
      }
    else
      {
        $landscape = "";
      }
    #......................
    #...get needed execs...
    #......................
    #.................................................
    #...find needed execs (add some things to path)...
    #.................................................
    $gnuplot = which_exec( "gnuplot" );
    $ps2pdf  = which_exec( "ps2pdf" );
    #...........................
    #...beginning gnuplot_cmd...
    #...........................
    $gnuplot_cmd .= "#\n";
    $gnuplot_cmd .= "# gnuplot commands to create $$gnuplot_info_ref{ps_file}\n";
    $gnuplot_cmd .= "# from $$gnuplot_info_ref{data_file}.\n";
    $gnuplot_cmd .= "#\n";
    $gnuplot_cmd .= "set output '$$gnuplot_info_ref{ps_file}'\n";
    $gnuplot_cmd .= "set terminal postscript $landscape color 'Times-Roman' 10\n";
    $gnuplot_cmd .= "set datafile missing '-'\n";
    #.......................................
    #...process each variable in own page...
    #.......................................
    foreach $variable ( @{$$gnuplot_info_ref{variables}} )
      {
        #..................................................
        #...get number of dtype plots and which are used...
        #..................................................
        $num_plots = 0;
        undef( %dtypes );
        foreach $source ( sort @{$$gnuplot_info_ref{sources}} )
          {
            foreach $dtype ( keys %{$$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDY}{$GDIFF}} )
              {
                $dtypes{$dtype} = "";
              }
          }
        if( %dtypes )
          {
            $num_plots += keys( %dtypes );
          }
        #................................................
        #...add 1 if see original data from any source...
        #................................................
        foreach $source ( sort @{$$gnuplot_info_ref{sources}} )
          {
            if( defined( $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDY}{$GORG} ) )
              {
                $num_plots++;
                last;
              }
          }
        #......................................
        #...if no plots, go to next variable...
        #......................................
        if( $num_plots == 0 )
          {
            next;
          }
        #.....................................
        #...per page (variable) gnuplot_cmd...
        #.....................................
        $page++;
        $gnuplot_cmd .= "#\n";
        $gnuplot_cmd .= "# Variable: $variable \[page $page with $num_plots plot(s)\]\n";
        $gnuplot_cmd .= "#\n";
        if( $num_plots > 1 )
          {
            $gnuplot_cmd .= "  set multiplot\n";
          }
        $plot_num = 0;
        $size = sprintf( "%13.4e", 1/$num_plots );
        $origin = sprintf( "%13.4e", 1 - $size);
        #.................................................................
        #...process original data then each dtype plot in correct order...
        #.................................................................
        foreach $dtype ( $GORG, @{$$gnuplot_info_ref{dtypes}} )
          {
            if( $plot_num > 0 && ! defined( $dtypes{$dtype} ) )
              {
                next;
              }
            $plot_num++;
            #.................
            #...find labels...
            #.................
            foreach $source ( sort @{$$gnuplot_info_ref{sources}} )
              {
                $xlabel =
                  $$gnuplot_info_ref{label}{$variable}{$source}{$GCOORDX}{$GORG};
                if( defined( $xlabel ) )
                  {
                    last;
                  }
              }
            if( ! defined( $xlabel ) )
              {
                $xlabel = "$GCOORDX";
              }
            foreach $source ( sort @{$$gnuplot_info_ref{sources}} )
              {
                if( $plot_num > 1 )
                  {
                    $ylabel =
                      $$gnuplot_info_ref{label}{$variable}{$source}{$GCOORDY}{$GDIFF}{$dtype};
                  }
                else
                  {
                    $ylabel =
                      $$gnuplot_info_ref{label}{$variable}{$source}{$GCOORDY}{$GORG};
                  }
                if( defined( $ylabel ) )
                  {
                    last;
                  }
              }
            if( ! defined( $ylabel ) )
              {
                $ylabel = "$GCOORDY";
              }
            foreach $source ( sort @{$$gnuplot_info_ref{sources}} )
              {
                $title = $$gnuplot_info_ref{title}{$variable}{$source};
                if( defined( $title ) )
                  {
                    last;
                  }
              }
            if( ! defined( $title ) )
              {
                $title = "$variable";
              }
            $title = "$title [$dtype]";
            #.................
            #...plot header...
            #.................
            $gnuplot_cmd .= "  #\n";
            $gnuplot_cmd .= sprintf( "  # Title: [%s] Plot Num [%s/%s]\n",
                                     $title, $plot_num, $num_plots );
            $gnuplot_cmd .= "  #\n";
            $gnuplot_cmd .= "  set title '$title'\n";
            $gnuplot_cmd .= "  set size 1,$size\n";
            $gnuplot_cmd .= "  set format y '\%15.8e'\n";
            $gnuplot_cmd .= "  set origin 0,$origin\n";
            $gnuplot_cmd .= "  set xlabel '$xlabel'\n";
            $gnuplot_cmd .= "  set ylabel '$ylabel'\n";
            $gnuplot_cmd .= "  plot ";
            #............................................
            #...process each source part of dtype plot...
            #............................................
            foreach $source ( sort @{$$gnuplot_info_ref{sources}} )
              {
                $using_x = $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDX}{$GORG};
                if( $plot_num > 1 )
                  {
                    $using_y =
                      $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDY}{$GDIFF}{$dtype};
                  }
                else
                  {
                    $using_y =
                      $$gnuplot_info_ref{using}{$variable}{$source}{$GCOORDY}{$GORG};
                  }
                if( defined( $using_y ) )
                  {
                    #............................
                    #...per source gnuplot_cmd...
                    #............................
                    $gnuplot_cmd .= "'$$gnuplot_info_ref{data_file}' ";
                    $gnuplot_cmd .= "index $$gnuplot_info_ref{index}{$variable}{$source} ";
                    $gnuplot_cmd .= "using ";
                    if( defined( $using_x ) )
                      {
                        $gnuplot_cmd .= "$using_x:";
                      }
                    $gnuplot_cmd .= "$using_y ";
                    $gnuplot_cmd .= "title \"$source\" ";
                    $gnuplot_cmd .= "with lines ";
                    $gnuplot_cmd .= "lt $$gnuplot_info_ref{lt}{source}{$source} ";
                    $gnuplot_cmd .= ", ";
                  }
              }
            #..............................
            #...finish dtype gnuplot_cmd...
            #..............................
            $gnuplot_cmd =~ s/, $/\n/;
            $origin = sprintf( "%13.4e", $origin - $size );
            #..................................................
            #...DONE: process each source part of dtype plot...
            #..................................................
          }
        #........................................................
        #...finish per page (variable) gnuplot_cmd gnuplot_cmd...
        #........................................................
        if( $num_plots > 1 )
          {
            $gnuplot_cmd .= "  set nomultiplot\n";
          }
        #...................................
        #...DONE: process each dtype plot...
        #...................................
      }
    #.............................................
    #...DONE: process each variable in own page...
    #.............................................
    #...............................................
    #...print gnuplot_cmd to file and run gnuplot...
    #...............................................
    open ( FILE, ">$$gnuplot_info_ref{cmd_file}" );
    print FILE $gnuplot_cmd;
    close FILE;
    $output = `$gnuplot $$gnuplot_info_ref{cmd_file} 2>&1`;
    if( $ps2pdf ne "" && -T "$$gnuplot_info_ref{ps_file}" )
      {
        #.......................................
        #...fix for doing landscape correctly...
        #.......................................
        $fix_file = "$$gnuplot_info_ref{ps_file}.tmp.ps";
        $output = `head -1 $$gnuplot_info_ref{ps_file}`;
        if( ! open( FILE, ">$fix_file" ) )
          {
            $ierr = 1;
            &print_error( "Cannot write to temporary file [$fix_file]",
                          $ierr );
            exit( $ierr );
          }
        print FILE $output;
        print FILE "<</AutoRotatePages /None>>setdistillerparams\n";
        close FILE;
        $output = `cat $$gnuplot_info_ref{ps_file} >> $fix_file`;
        $output = `$ps2pdf $fix_file $$gnuplot_info_ref{pdf_file} 2>&1`;
        unlink( $fix_file );
      }
  }
#............................................................................
#...Name
#...====
#... which_exec
#...
#...Purpose
#...=======
#... Returns value for full path to executable ("" if not found).
#...
#...Arguments
#...=========
#... exec Intent: in
#...      Perl type: string scalar
#...
#...Program Flow
#...============
#... 1) go through path to find exec
#............................................................................
sub which_exec
  {
    my( 
       $exec, # exec to find
      ) = @_;
    my(
       $ierr, # error return value
       $exec_try, # see if exec is here
       $found, # if found
       $path, # current dir in search for execs
       @paths, # list of search paths for execs
       $this_dir, # current directory
      );
    #.................................................
    #...build paths from PATH and current directory...
    #.................................................
    $path = $ENV{PATH};
    @paths = split( /:/, $ENV{PATH} );
    ($this_dir = $0) =~ s/\/[^\/]+$//;
    unshift( @paths, $this_dir );
    #.....................................
    #...loop through paths to find exec...
    #.....................................
    $found = "false";
    foreach $path (@paths)
      {
        $exec_try = "$path/$exec";
        if( -x "$exec_try" )
          {
            $found = "true";
            last;
          }
      }
    #........................................
    #...error if still could not find exec...
    #........................................
    if( $found eq "false" )
      {
        $ierr = 0;
        &print_error(
                     "Necessary executable [$exec] not found in PATH",
                     $ierr
                    );
        $exec_try = "";
      }
    return( $exec_try );
  }


1;

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

