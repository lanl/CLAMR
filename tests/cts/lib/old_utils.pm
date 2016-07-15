package old_utils;

use     POSIX qw( strtod );
#use     diagnostics;
use     warnings;
use     Carp;
use     vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );
use     Exporter;
use     Cwd;
# not needed right now...but could add in
#use     File::Spec;

$VERSION   = 1.00;

@ISA       = qw(
                Exporter
               );

@EXPORT    = qw(
               );

@EXPORT_OK = qw(
                &conv_time
                &datafile_files
                &datafile_parse
                &datafile_getval
                &datafile_getblock
                &datafile_setcond
                &datafile_debug
                &date_ymdhms
                &date_ymdhms_sep
                &expand_string
                &extrema
                &get_pname
                &get_sysinfo
                &get_jobinfo
                &my_compare_version
                &my_dir
                &my_numerically
                &my_mkdir
                &my_mkdir_hpss
                &my_mode
                &my_munge
                &my_chdir
                &my_copy
                &my_copy_hpss
                &my_getval
                &my_lockfile
                &my_notdir
                &my_stat
                &my_xml_read
                &path_add_default
                &print_perl_obj
                &print_error
                &run_command
                &sort_unique
                &status_bar
                &which_exec
                &latexify
               );
sub BEGIN
  {
    ($lib_path = $0) =~ s/\/[^\/]*$//;
    unshift( @INC, $lib_path )
  }
sub END
  {
  }
#......................
#...global variables...
#......................
my(
   %G_RUN_COMMAND_SEEN,
  );
# print progress bar given current count (start at 1) and total count
sub status_bar{
    if( $_[0] == 1 ){
        my $spaces = $_[2]||"";
        print "${spaces}Status 1..9: ";
    }
    print substr( "123456789\n", int(($_[0]-1)/($_[1]/10.)),(int($_[0]/($_[1]/10.)) - int(($_[0]-1)/($_[1]/10.))) );
}
# compare 2 version numbers:
#   op  == 1 -> returns ($ver_a > $ver_b)
# cat get ( "=" ) by ( not ">" and not "<" )
#
# Just split on non-digits (12.4.7 > 12 3 sub release 854 > 6.54 > .99999 )
# Currently, cannot handle 12a vs 12b
sub my_compare_version{
    my( $ver_a, $op, $ver_b ) = @_;
    my(
        @fields_a,
        @fields_b,
        $i,
        $num,
        $num_a,
        $num_b,
        $ret,
        );
    if( ! defined($ver_a) || $ver_a !~ /\S/ ){
        $ver_a = "0";
    }
    if( ! defined($ver_b) || $ver_b !~ /\S/ ){
        $ver_b = "0";
    }
    $ver_a =~ s/^\s+//;
    $ver_a =~ s/\s+$//;
    $ver_b =~ s/^\s+//;
    $ver_b =~ s/\s+$//;
    @fields_a = split( /\D+/, $ver_a );
    @fields_b = split( /\D+/, $ver_b );
    $num = $#fields_a;
    if( $num < $#fields_b ){
        $num = $#fields_b;
    }
    $ret = eval( "1 $op 1" );
    for( $i = 0; $i <= $num; $i++ ){
        $num_a = $fields_a[$i];
        $num_b = $fields_b[$i];
        if( ! defined( $num_a ) ){
            $num_a = 0;
        }
        if( ! defined( $num_b ) ){
            $num_b = 0;
        }
        if( $num_a != $num_b ){
            return( eval( "$num_a $op $num_b" ) );
        }
    }
    $ret;
}
sub date_ymdhms{
    my( $date ) = `date +%Y%m%d%H%M%S`;
    chomp( $date );
    return( $date );
}
sub date_ymdhms_sep{
    my( $date ) = `date +%Y.%m.%d.%H.%M.%S`;
    chomp( $date );
    return( $date );
}
sub get_jobinfo{
    my(
       $jobid,
       $jobinfo_ref,
       ) = @_;
    my(
       $cmd,
       $hosts,
       $output,
       );
    $cmd = "checkjob $jobid";
    $output = &run_command( COMMAND=>$cmd );
    $hosts = '';
    if( $output =~ /Allocated Nodes:\s+(\[\S+)/ ){
        $hosts = $1;
    }
    $$jobinfo_ref{hosts_string} = $hosts;
}

# add some default dirs to the PATH
sub path_add_default{
    my( $dir );
    $dir = &my_dir( $0 );
    # directory of script
    $ENV{'PATH'} .= ":$dir";
    # checkout of build (for Makefile_vars.mk)
    $ENV{'PATH'} .= ":$dir/../../Source.rh/build";
}

# get a lockfile
# will block until gotten
sub my_lockfile{
    my %args = (
                LOCKFILE => "./lockfile", # filename
                QUIET    => undef,    # "", "success", "errors"
                REASON   => undef,    # string to put into lockfile
                TIME     => 99999999, # secs
                TYPE     => 1,        # how it locks
                WARNING  => undef,    # if just return, not exit
                @_,
                );
    my $args_valid =("LOCKFILE|QUIET|REASON|TIME|TYPE|WARNING");
    my(
       $lockfile,
       $time,      # seconds
       ) = @_;
    my(
       $arg,
       $count,
       $dir,
       $done,
       $file,
       @files,
       $group,
       $groups_regexp,
       $ierr,
       $lockfile_dir,
       $lockfile_exec,
       $lockfile_notdir,
       $lockfile_perm,
       $output,
       $sleep,
       %stat,
       );

    # init
    $ierr = 0;

    # args
    foreach $arg ( keys %args ){
        if( $arg !~ /^$args_valid$/ ){
            $ierr = 1;
            &print_error( "Invalid argument [$arg]",
                          "Valid args [$args_valid]",
                          $ierr );
            exit( $ierr );
        }
    }

    $lockfile = $args{LOCKFILE};
    $lockfile_notdir = &my_notdir( $lockfile );

    # get group name of parent directory
    $lockfile_dir = &my_dir( $lockfile );
    $ierr = &my_stat( $lockfile_dir, \%stat );
    if( $ierr != 0 ){
        if( defined( $args{WARNING} ) ){
            $ierr = 0;
        }
        if( ! defined( $args{QUIET} ) || $args{QUIET} !~ /^(|errors)$/ ){
            &print_error( "Cannot get info about: ",
                          "  parent directory [$lockfile_dir]",
                          "  of the lockfile  [$lockfile]",
                          "Lockfile not created", $ierr );
        }
        if( ! defined( $args{WARNING} ) ){
            exit( $ierr );
        }
        else{
            return( $ierr );
        }
    }
    # Some parent directories might have a gid without a name
    # in this case, $group will be undefined (and no chgrp done)
    $group = $stat{group};

    # see if can set to this groups
    $groups_regexp = `groups 2>&1`;
    $groups_regexp =~ s/\s+/\|/g;
    if( defined( $group ) && $group !~ /^($groups_regexp)$/ ){
        undef( $group );
    }

    # set permissions for lockfile
    # if can set the group, then group write.
    # if not, then other write (so others can remove lockfile if needed)
    if( defined( $group ) ){
        $lockfile_perm = 0660;
    }
    else{
        $lockfile_perm = 0666;
    }

    $time = $args{TIME};
    $sleep = 1;

    # get lockfile exec if available
    if( $args{TYPE} == 1 ){
        $lockfile_exec = &which_exec( "lockfile", QUIET=>"" );
        # if not there, default to another type
        if( $lockfile_exec eq "" ){
            $args{TYPE} = 2;
        }
    }

    if( ! defined( $args{QUIET} ) || $args{QUIET} !~ /^(|success)$/ ){
        print "my_lockfile start ".&date_ymdhms_sep()." [$lockfile]";
        if( -e $lockfile ){
            &my_stat( $lockfile, \%stat );
            print " wait [$time secs] > age [$stat{mtime_since} secs]\n";
            $output = `cat $lockfile 2>&1`;
            $output =~ s/\s*$//;
            print "\n";
            print "lockfile contents:\n";
            print "------------------\n";
            print "$output\n\n";
        }
        else{
            print "\n";
        }
    }

    # this call takes about $sleep + a few
    if( $args{TYPE} == 1 ){
        `$lockfile_exec -s $sleep -l $time $lockfile`;
        chmod( $lockfile_perm, $lockfile );
        if( defined($group) ){
            `chgrp $group $lockfile 2>&1`;
        }
    }

    # home grown approach - try to be safe, but again, takes about $sleep + a few
    elsif( $args{TYPE} == 2 ) {
        $done = "false";
        while( $done eq "false" ){
            while( -e $lockfile ){
                &my_stat( $lockfile, \%stat );
                # just remove file if older than time
                if( defined( $stat{mtime_since} ) && $stat{mtime_since} > $time ){
                    unlink( $lockfile );
                }
            }
            `touch $lockfile.$$`;
            chmod( $lockfile_perm, "$lockfile.$$" );
            if( defined($group) ){
                `chgrp $group $lockfile.$$ 2>&1`;
            }
            `touch $lockfile`;
            chmod( $lockfile_perm, $lockfile );
            if( defined($group) ){
                `chgrp $group $lockfile 2>&1`;
            }
            # try to refresh directory and files
            sleep( $sleep );
            `ls $lockfile* 2>&1 > /dev/null`;
            &my_stat( "$lockfile.$$", \%stat );
            $dir = $stat{dir};
            if( ! opendir( MY_LOCKFILE_DIR, $dir ) ){
                $ierr = 1;
                &print_error( "Cannot open lockfile directory [$dir]",
                              $ierr );
                exit( $ierr );
            }
            @files = grep ( /^${lockfile_notdir}\.\d+$/ , readdir( MY_LOCKFILE_DIR ) );
            closedir( MY_LOCKFILE_DIR );
            # remove old files
            $count = 0;
            foreach $file ( @files ){
                &my_stat( $file, \%stat );
                if( defined($stat{mtime_since}) && $stat{mtime_since} > $time ){
                    unlink( $file );
                }
                else{
                    $count++;
                }
            }
            # done only if just this process has lockfile.id
            if( $count == 1 ){
                $done = "true";
            }
            unlink( "$lockfile.$$" );
        }
    }
    
    # otherwise do fastest - not very safe
    else{
        while( -e $lockfile ){
            &my_stat( $lockfile, \%stat );
            if( ! defined($stat{mtime_since}) || $stat{mtime_since} > $time ){
                last;
            }
        }
        `touch $lockfile 2>&1`;
        chmod( $lockfile_perm, $lockfile );
        if( defined($group) ){
            `chgrp $group $lockfile 2>&1`;
        }
    }

    # got it - message
    if( ! defined( $args{QUIET} )  || $args{QUIET} !~ /^(|success)$/ ){
        print "my_lockfile stop  ".&date_ymdhms_sep()." [$lockfile]";
        print "\n";
    }

    # put message into lockfile
    if( open( MY_LOCKFILE_FILE, ">$lockfile" ) ){
        if( defined( $args{REASON} ) ){
            print MY_LOCKFILE_FILE $args{REASON},"\n";
        }
        print MY_LOCKFILE_FILE "LOGNAME: $ENV{LOGNAME}\n";
        print MY_LOCKFILE_FILE "EXEC:    $0\n";
        print MY_LOCKFILE_FILE "DATE:    ",&date_ymdhms_sep(),"\n";
        print MY_LOCKFILE_FILE "UNAME:   ",`uname -n`,"\n";
        close( MY_LOCKFILE_FILE );
    }

    return( $ierr );
}
#............................................................................
#...Name
#...====
#... my_mode
#...
#...Purpose
#...=======
#... Returns a string of the correct mode given various inputs
#...
#...Arguments
#...=========
#... STRING_REF   Intent: in
#...              Perl type: reference to scalar
#...              : delimited string
#... PRE          Intent: in
#...              Perl type: scalar
#...              : delimited string to place in front
#... PST          Intent: in
#...              Perl type: scalar
#...              : delimited string to place behind
#... EXIST        Intent: in
#...              Perl type: scalar
#...              If defined, checke if it exists before adding it
#...
#............................................................................
sub my_mode{
    my %args = (
                DEC   => undef,
                DIR   => undef,
                EXEC  => undef,
                UMASK => undef,
                @_,
                );
    my $args_valid =("DEC|DIR|EXEC|UMASK");
    my(
       $arg,
       $ierr,
       $mode,
       );
    # args
    foreach $arg ( keys %args ){
        if( $arg !~ /^$args_valid$/ ){
            $ierr = 1;
            &print_error( "Invalid argument [$arg]",
                          "Valid args [$args_valid]",
                          $ierr );
            exit( $ierr );
        }
    }
    if( ! defined($args{UMASK}) ){
        $args{UMASK} = sprintf( "%lo", umask());
    }
    if( defined($args{DIR} ) ){
        $mode = oct("2777") - oct($args{UMASK});
    }
    else{
        $mode = oct("777") - oct($args{UMASK});
        if( ! defined($args{EXEC}) ){
            $mode &= oct("666");
        }
    }
    if( defined($args{DEC}) ){
        return( $mode );
    }
    else{
        return( sprintf( "%lo", $mode ) );
    }
}
#............................................................................
#...Name
#...====
#... my_munge
#...
#...Purpose
#...=======
#... Takes a : delimited string, perpends or postpends args, and removes any duplicates
#...
#...Arguments
#...=========
#... STRING_REF   Intent: in
#...              Perl type: reference to scalar
#...              : delimited string
#... PRE          Intent: in
#...              Perl type: scalar
#...              : delimited string to place in front
#... PST          Intent: in
#...              Perl type: scalar
#...              : delimited string to place behind
#... EXIST        Intent: in
#...              Perl type: scalar
#...              If defined, checke if it exists before adding it
#...
#............................................................................
sub my_munge{
    my %args = (
                STRING_REF  => undef,
                PRE         => undef,
                PST         => undef,
                EXIST       => undef,
                @_,
                );
    my(
       $ierr,
       $arg,
       );
    my $args_valid =("EXIST|PST|PRE|STRING_REF");
    my(
       $exist,
       %seen,
       $val,
       @vals,
       @vals_res,
       );
    # args
    foreach $arg ( keys %args ){
        if( $arg !~ /^$args_valid$/ ){
            $ierr = 1;
            &print_error( "Invalid argument [$arg]",
                          "Valid args [$args_valid]",
                          $ierr );
            exit( $ierr );
        }
    }
    # pre
    if( defined($args{PRE}) ){
        @vals = split( /\s*:\s*/, $args{PRE} );
    }
    # string
    push( @vals, split( /\s*:\s*/, ${$args{STRING_REF}} ) );
    # post
    if( defined($args{PST}) ){
        push( @vals, split( /\s*:\s*/, $args{PST} ) );
    }
    foreach $val ( @vals ){
        # skip nulls
        if( $val !~ /\S/ ){
            next;
        }
        # skip if already seen
        if( defined( $seen{"$val"} ) ){
            next;
        }
        # skip if not exist
        if( defined($args{EXIST}) ){
            eval( "\$exist = $args{EXIST} \"$val\"" );
            if( ! $exist ){
                next;
            }
        }
        push( @vals_res, $val );
        $seen{"$val"} = "";
    }
    ${$args{STRING_REF}} = join( ":", @vals_res );
}
#............................................................................
#...Name
#...====
#... get_pname
#...
#...Purpose
#...=======
#... Try to determine the problem name given the files in the directory.
#... Returns the problem name (or '')
#...
#...Arguments
#...=========
#... $dir      Intent: in
#...           Perl type: scalar
#...           directoryt path to look
#...
#............................................................................
sub get_pname{
  my(
     $dir,
     ) = shift(@_);
  my %args = (
              INFILE  => undef,
              @_,
              );
  my $args_valid =("INFILE");
  my(
     $file,
     @files_all,
     @files_all_sort,
     $ierr,
     $line,
     @lines,
     $line_new,
     @lines_new,,
     $max,
     $notdir,
     $output,
     $pname,
     $pname_guess,
     %seen,
     %seen_type,
     %stat,
     $type,
    );
  # args
  foreach $arg ( keys %args ){
      if( $arg !~ /^$args_valid$/ ){
          $ierr = 1;
          &print_error( "Invalid argument [$arg]",
                        "Valid args [$args_valid]",
                        $ierr );
          exit( $ierr );
      }
  }
  $pname = "";
  # get names of all files in this dir
  if( ! defined($dir) ) {
      return;
  }
  # try to get pname from the input file
  if( $pname !~ /\S/ ){
      if( -e "$dir/run_job.csh" ){
          $output = `egrep 'RJ_CMD_PRUN.*\.in*' $dir/run_job.csh`;
          # put last first to find the last occurrence
          @lines = reverse split( /\n/, $output );
          foreach $line ( @lines ){
              if( $line =~ /RJ_CMD_PRUN.*\s(\S+\.in\S*)/ ){
                  $file = $1;
                  $file =~ s/[\'\"]//g;
                  if( $file !~ /^\// ){
                      $file = "$dir/$file";
                  }
                  if( -e "$file" ){
                      $output = `egrep '^[ \t]*pname[ \t]*=' $file`;
                      @lines_new = reverse split( /\n/, $output );
                      foreach $line_new ( @lines_new ){
                          if( $line_new =~ /^\s*pname\s*=\s*['"](\S+?)["']/ ){
                              $pname = $1;
                              # skip any variable setting
                              if( $pname =~ /\$/ ){
                                  $pname = "";
                              }
                              else{
                                  last;
                              }
                          }
                      }
                      # found from previous foreach - exit loop
                      if( $pname =~ /\S/ ){
                          if( defined($args{INFILE}) ){
                              (${$args{INFILE}} = $file) =~ s/\.\///g;
                          }
                          last;
                      }
                  }
              }
          }
      }
  }
  # try to get from strings in the output file
  if( $pname !~ /\S/ ){
      if( -e "$dir/rj_cmd_out" ){
          $output = `grep "Open output file:" $dir/rj_cmd_out`;
          if( defined($output) ){
              # put last first to find the last occurrence
              $output = join ( "\n", reverse split( /\n/, $output ) );
              if( $output =~ /Open output file:\s*(\S+)-output\b/ ){
                  $pname = $1;
              }
          }
      }
  }
  # get from names of files
  if( $pname !~ /\S/ ) {
      if( ! opendir(DIR, $dir) ) {
          $ierr = 0;
          &print_error( "Cannot get directory info [$dir]...skipping",
                        $ierr );
          return( $pname );
      }
      @files_all = grep ! /^\.\.?$/, readdir( DIR );
      closedir( DIR );
      # full path name for sorting
      # if cannot get time, set time to large number so it will be last (soft link pointing to nowhere)
      @files_all_sort = sort { (-M "$dir/$a"||9999) <=> (-M "$dir/$b"||9999) } @files_all;
      undef( %seen );
      undef( %seen_type );
      foreach $file ( @files_all_sort ) {
          if( $file =~ /^gold_/ ){
              next;
          }
          # if file is of a certain type (pname)-(type)
          if( $file =~ /^(\S+)\-(build|dmp|editmix|history|lastcycle|lastdump|output|problemsize|status|DO_NOT_RUN)/ ){
              $pname_guess = $1;
              $type = $2;
              # only get last created file of a certain type
              if( defined( $seen_type{$type} ) ){
                  next;
              }
              $seen_type{$type} = '';
              # add it to the count seen for that pname
              if( ! defined($seen{$pname_guess}) ) {
                  $seen{$pname_guess} = 0;
              }
              $seen{$pname_guess}++;
          }
      }
      # pick one with the highest number of files
      $max = 0;
      foreach $pname_guess ( keys %seen ) {
          if ( $seen{$pname_guess} > $max ) {
              $pname = $pname_guess;
              $max = $seen{$pname_guess};
          }
      }
  }
  # get from env var
  if( $pname !~ /\S/ ){
      if( defined($ENV{RJ_VAR_PNAME}) ){
          $pname = $ENV{RJ_VAR_PNAME};
      }
  }
  # get from last directory name
  if( $pname !~ /\S/ ){
      &my_stat( $dir, \%stat );
      if( defined($stat{notdir}) ){
          $pname = $stat{notdir};
      }
  }
  return( $pname );
}
#............................................................................
#...Name
#...====
#... expand_string
#...
#...Purpose
#...=======
#... Takes a string with a regexp pattern and returns an array of values.
#...    foo[1-3],bar[2-4,7]a => foo1,foo2,foo3,bar2a,bar3a,bar4a,bar7a
#...
#...Arguments
#...=========
#... $string      Intent: in
#...              Perl type: scalar
#...              regexp'd string
#...
#... $array_string_ref    Intent: inout
#...                      Perl type: reference to array
#...
#............................................................................
sub expand_string{
    my(
       $string,
       $array_string_ref,
       ) = @_;
    my(
       $a,
       $b,
       $mid,
       $pre,
       $pst,
       $val,
       $string_result,
       @vals,
       @words,
       );
    # remove ppn constructs
    $string =~ s/:\d+//g;
    $string =~ s/\*\d+//g;
    # if string is from checkjob, put it in ljobs type form
    # [hosta:ppn][hostb:ppn][hostc:ppn] => hosta,hostb,hostc
    if( $string =~ /^(\[[^\[\]]+\])+$/ ){
        $string_result = "";
        while( $string =~ /^(\[[^\[\]]+\])(.*)$/ ){
            $string_result .= ",$1";
            $string = $2;
        }
        $string = $string_result;
        $string =~ s/^,//;
    }
    # ljobs form: "," separated and expand stuff inside "[]"
    $a = $string;
    $a =~ s/^\s*//;
    $a =~ s/\s*$//;
    $b = "";
    while( $a ne $b )
    {
        $b = $a;
        #...replace , with ; inside [] so can split on ,
        if( $a =~ /^(.*)(\[[^\[\]]*,[^\[\]]*\])(.*)$/ )
        {
            $pre = $1;
            $mid = $2;
            $pst = $3;
            $mid =~ s/,/;/g;
            $a = "$pre$mid$pst";
        }
        #...replace - with .. inside [] (perl range operator)
        if( $a =~ /^(.*)(\[[^\[\]]*\-[^\[\]]*\])(.*)$/ )
        {
            $pre = $1;
            $mid = $2;
            $pst = $3;
            $mid =~ s/\-/\.\./g;
            $a = "$pre$mid$pst";
        }
    }
    @vals = split(/\s*,\s*/, $a );
    foreach $val( @vals )
    {
        #...put , back
        $val =~ s/;/,/g;
        #...replace [] with a range into an array
        if( $val =~ /^(.*)\[(.*)\](.*)$/ )
        {
            $pre = $1;
            $mid = $2;
            $pst = $3;
            $mid =~ s/(\w+)/'$1'/g;
            @words = eval($mid);
            grep( s/^/${pre}/, @words );
            grep( s/$/${pst}/, @words );
        }
        #...or just stick value onto words
        else
        {
            @words = ($val);
        }
        push( @{$array_string_ref}, @words );
    }
}
sub my_chdir
  {
    my(
       $dir
      ) = @_;
    my(
       $ierr,
      );
    if( ! chdir( $dir ) )
      {
        $ierr = 1;
        &print_error( "Cannot chdir to [$dir]",
                      $ierr );
        exit( $ierr );
      }
    return( &cwd() );
  }
# numerically: for sorting...
# @foo = ( "1\n","5\n","12\n" );
# print "numerically\n", sort my_numerically @foo;
# print "default\n", sort @foo;
#
# numerically
# 1
# 5
# 12
# default
# 1
# 12
# 5
sub my_numerically { $a <=> $b; }
#............................................................................
#...Name
#...====
#... datafile_parse
#...
#...Purpose
#...=======
#... Parse data file(s) and stuff into $data_file_ref
#... Data File Format:
#...   # comment
#...   a = b \ # variable with line continuation
#...       c
#...   block: code # start the "code" block
#...   code source pack prefix group # header for this block
#...   a    b      c    d      e     # values for a line
#...   f    g      h    i      j     # values for another line
#...   block: # end of the block.
#...   block.: code{TU} # add to code block if "TU"
#...   code source pack prefix group # header for this block
#...   k    l      m    n      o     # values for a line
#...   block: # end of the block.
#...   
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: in
#...              Perl type: reference to hash
#...              command line
#...              $cmd{$option} = value
#...              $cmd{files}[] = array of file names
#...
#... $data_file_ref    Intent: inout
#...                   Perl type: reference to has
#...                   Necessary Fields:
#...                     {blocks_valid}{<block namd>} =
#...                        "<field>|<field>|<field>..."
#...                     {files} =
#...                         ("<datafile1>", "<datafile2>", "<datafile2>", ... )
#...                   Output Fields:
#...                     {blocks}{block}[<line num>]{<key>} = <value>
#...                     {key_val}{<key>}{vals}{<condition>} = <value>
#...                     {key_val}{<key>}{order}[] = <condition>
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub datafile_parse
  {
    my(
       $cmd_ref,
       $data_file_ref
      ) = @_;
    my(
       $active,
       $block,
       $block_line,
       $blocks_valid,
       %blocks,
       $cond,
       $done,
       $done_block,
       @fields,
       $file,
       $i,
       $ierr,
       $j,
       $key,
       %key_val,
       $line_num,
       $line,
       $line_new,
       @lines_file,
       $op,
       $val,
       @vals,
       $vals_str,
       $vals_str_valid,
      );
    # default COND
    if( ! defined( $COND ) ){
        &datafile_setcond( "" );
    }
    $DEFAULT_DATA_FILE_REF = $data_file_ref;
    $blocks_valid = join("|",sort keys %{$$data_file_ref{blocks_valid}});
    # open and process each file
    foreach $file ( @{$$data_file_ref{files}} )
      {
        print "Datafile: [$file]";
        if( ! open( FILE, $file ) )
          {
            print " skipping\n";
            next;
          }
        print "\n";
        @lines_file = <FILE>;
        close( FILE );
        $i = 0;
        $done = "false";
        # process file
        while( $done eq "false" )
          {
            #...if done with file
            if( $i > $#lines_file )
              {
                $done = "true";
                next;
              }
            &_datafile_get_next_line( \$i, \@lines_file, \$line );
            # block
            if( $line =~ /^\s*block\s*(\.?:)\s*(.*)/ )
              {
                $op = $1;
                $block = $2;
                if( $block =~ /^(\S+)\s*\{\s*(.*?)\s*\}\s*$/ ) {
                    $block = $1;
                    $cond = $2;
                }
                else {
                    $cond = "";
                }
                if( $block !~ /^($blocks_valid)$/ )
                  {
                    $ierr = 1;
                    $line_num = $i;
                    &print_error( "Invalid block [$block]",
                                  "Valid blocks: [$blocks_valid]",
                                  "File: [$file:$line_num]",
                                  "line: [$line]",
                                  $ierr );
                    exit( $ierr );
                  }
                # see if this block is active
                my( $cond_eval ) = $cond;
                # go through each condition set and replace that string with true (1)
                # if that condition is set, or false (0) if not.
                # First will use ";" for "true", then replace all other strings with
                # false (0), then replace ";" with true (1)
                # replace pure numbers with true (right now, ";")
                $cond_eval =~ s/\b[1-9]\d*/\;/g;
                # replace each condition set with true (;)
                foreach my $cond_set ( split( /\s*,\s*/, $COND ) ){
                    $cond_eval =~ s/\b$cond_set\b/\;/g;
                }
                # replace all other words with false (0)
                $cond_eval =~ s/\b\w+\b/0/g;
                # set ";" back to true (1)
                $cond_eval =~ s/\;/1/g;
                # if non-null, then eval it
                if( $cond_eval =~ /\S/ ){
                    eval( "\$cond_eval = $cond_eval" );
                }
                # default is null being true
                else{
                    $cond_eval = 1;
                }
                if( $cond_eval ){
                    $active = "true";
                }
                else{
                    $active = "false";
                }
                #...init this block (will add to it if ".:")
                if( $op eq ":" && $active eq "true" ) {
                    @{$blocks{$block}} = ();
                }
                &_datafile_get_next_line( \$i, \@lines_file, \$line );
                if( $line !~ /\S/ ) {
                    $ierr = 1;
                    $line_num = $i;
                    &print_error( "Missing column header for block [$block]",
                                  "File: [$file:$line_num]",
                                  "line: [$line]",
                                  $ierr );
                    exit( $ierr );
                }
                @fields = split( /\s+/, $line );
                # check valid fields
                foreach $field ( @fields ) {
                    if( $field !~ /^($$data_file_ref{blocks_valid}{$block})$/ ) {
                        $ierr = 1;
                        $line_num = $i;
                        &print_error( "Missmatch of fields for block [$block]",
                                      "field:        [$field]",
                                      "fields valid: [".
                                      $$data_file_ref{blocks_valid}{$block}."]",
                                      "File: [$file:$line_num]",
                                      "line: [$line]",
                                      $ierr );
                        exit( $ierr );
                    }
                }
                $done_block = "false";
                $block_line = 0;
                # block_line
                while( $done_block eq "false" ) {
                    &_datafile_get_next_line( \$i, \@lines_file, \$line );
                    if( $line !~ /\S/ ) {
                        $ierr = 1;
                        $line_num = $i;
                        &print_error( "Missing data for block [$block]",
                                      "File: [$file:$line_num]",
                                      "End of File",
                                      $ierr );
                        exit( $ierr );
                    }
                    if( $line eq "block:" ) {
                        $done_block = "true";
                        next;
                    }
                    $line_new = $line;
                    @vals = ();
                    while( $line_new =~ /\S/ ){
                        if( $line_new =~ /^\s*\"([^\"]*?)\"(.*)/ ){
                            $val = $1;
                            $line_new = $2;
                            push( @vals, $val );
                        }
                        elsif( $line_new =~ /^\s*(\S+)(.*)/ ){
                            $val = $1;
                            $line_new = $2;
                            push( @vals, $val );
                        }
                    }
                    if( $#fields !~ $#vals ) {
                        $ierr = 1;
                        $line_num = $i;
                        &print_error( "Missmatch of fields and values for block [$block]",
                                      "fields: ".join(', ', @fields),
                                      "vals:   ".join(', ', @vals),
                                      "File: [$file:$line_num]",
                                      "line: [$line]",
                                      $ierr );
                        exit( $ierr );
                    }
                    if( $active eq "true" ){
                        $block_line = $#{$blocks{$block}} + 1;
                        $j = 0;
                        foreach $field ( @fields ) {
                            # if the block is active, save it
                            $blocks{$block}[$block_line]{$field} = &datafile_replace_var( $vals[$j] );
                            $j++;
                        }
                    }
                } # done: block_line
            } # done: block
            # key = val
            elsif( $line =~ /^(.*?)\s*(\.?=)\s*(.*?)$/ )
              {
                $key = $1;
                $op  = $2;
                $val = $3;
                if( $key =~ /^(\S+)\s*\{\s*(.*?)\s*\}\s*$/ )
                  {
                    $key  = $1;
                    $cond = $2;
                  }
                else
                  {
                    $cond = "";
                  }
                if( ! defined( $key_val{$key}{vals}{$cond} ) )
                  {
                    push( @{$key_val{$key}{order}}, $cond );
                  }
                if( $op eq ".=" )
                  {
                    if( defined( $key_val{$key}{vals}{$cond} ) )
                      {
                        $val = "$key_val{$key}{vals}{$cond} $val";
                      }
                  }
                $key_val{$key}{vals}{$cond} = $val;
              } # done: key = val
            # unparsed line
            elsif( $line =~ /\S/ )
              {
                $ierr = 1;
                $line_num = $i;
                &print_error( "Unparsed line",
                              "File: [$file:$line_num]",
                              "line: [$line]",
                              $ierr );
                exit( $ierr );
              }
          } # done: file
      } # open and process each file
    %{$$data_file_ref{blocks}} = %blocks;
    %{$$data_file_ref{key_val}} = %key_val;
  }
#...prune comments, join continuations, skip blank lines
sub _datafile_get_next_line
  {
    my(
       $i_ref,
       $lines_ref,
       $line_ref,
      ) = @_;
    my(
       $done,
       $line_new,
      );
    #...init
    $done = "false";
    while( $done eq "false" )
      {
        $$line_ref = "";
        #...end of file
        if( $$i_ref > $#{$lines_ref} )
          {
            $done = "true";
            next;
          }
        $$line_ref = $$lines_ref[$$i_ref];
        # leading/trailing whitespace
        $$line_ref =~ s/^\s*(.*?)\s*$/$1/;
        $$i_ref++;
        # comment line
        $$line_ref =~ s/\s*#.*//;
        # continuation line
        while( $$line_ref =~ /^(.*?)\s*\\\s*$/ )
          {
            $$line_ref = $1;
            $$line_ref =~ s/\s*$//;
            if( $$i_ref <= $#{$lines_ref} )
              {
                $line_new = $$lines_ref[$$i_ref];
                # comment
                $line_new =~ s/\s*#.*//;
                # leading trailing whitespace
                $line_new =~ s/^\s*(.*?)\s*$/$1/;
                $$line_ref .= " $line_new";
                #...leading/trailing whitespace
                $$line_ref =~ s/^\s*(.*?)\s*$/$1/;
              }
            $$i_ref++;
          }
        # skip blank
        if( $$line_ref =~ /\S/ )
          {
            $done = "true";
          }
      }
  }
#............................................................................
#...Name
#...====
#... datafile_files
#...
#...Purpose
#...=======
#... The list of files to be parsed in datafile_parse
#...
#...Arguments
#...=========
#... $data_file_ref    Intent: inout
#...                   Perl type: reference to hash
#...                   Will be used in datafile_parse
#...
#...Program Flow
#...============
#... 1) stuff array into data_file_ref
#............................................................................
sub datafile_files{
    my(
       $data_file_ref,
       @files
       );
    $data_file_ref = shift(@_);
    @files = @_;
    push( @{$$data_file_ref{files}}, @files );
}
#...return string containing datafile info
sub datafile_debug{
  my(
     $cmd_ref,
     $data_file_ref
    ) = @_;
  my(
     $block,
     @block_vals,
     $file,
     $i,
     $key,
     @keys,
     $val,
     %width,
     $output,
    );
  $output  = "\n";
  $output .= "==============\n";
  $output .= "Datafile Begin\n";
  $output .= "==============\n\n";
  $output .= "Files:\n";
  $output .= "------\n";
  foreach $file ( @{$$data_file_ref{files}} ){
    $output .= " $file\n";
  }
  $output .= "\n";
  $output .= "Condition: $COND\n";
  $output .= "\n";
  $output .= "Vars:\n";
  $output .= "-----\n";
  foreach $key ( sort keys %{$$data_file_ref{key_val}}) {
    $val = &datafile_getval($key);
    if( !defined($val) ){
      $val = "<UNDEF>";
    }
    $output .= sprintf( "%20s = %s\n", $key, $val );
  }
  $output .= "\n";
  $output .= "Blocks:\n";
  $output .= "-------\n";
  foreach $block (sort keys %{$$data_file_ref{blocks}} ) {
    $output .= " [$block]\n";
    @block_vals = &datafile_getblock($block);
    if( $#block_vals >= 0 ){
      @keys = sort( keys %{$block_vals[0]} );
      # get max width
      foreach $key ( @keys ) {
        $width{$key} = length($key);
      }
      for ( $i = 0; $i <= $#block_vals; $i++ ) {
        foreach $key ( @keys ) {
          $val = $block_vals[$i]{$key};
          if( $val !~ /\S/ ){
              $val = '""';
          }
          elsif( $val =~ /\s/ ){
              $val = '"'.$val.'"';
          }
          if( length($val) > $width{$key} ){
            $width{$key} = length($val);
          }
        }
      }
      foreach $key ( @keys ) {
        $output .= sprintf( "%$width{$key}s ", $key );
      }
      $output .= "\n";
      for ( $i = 0; $i <= $#block_vals; $i++ ) {
        foreach $key ( @keys ) {
            $val = $block_vals[$i]{$key};
            if( $val !~ /\S/ ){
                $val = '""';
            }
            elsif( $val =~ /\s/ ){
                $val = '"'.$val.'"';
            }
          $output .= sprintf( "%$width{$key}s ", $val );
        }
        $output .= "\n";
      }
      $output .= "\n";
    }
  }
  $output .= "============\n";
  $output .= "Datafile End\n";
  $output .= "============\n";

  return( $output );
}
#............................................................................
#...Name
#...====
#... datafile_getblock
#...
#...Purpose
#...=======
#... Returns block info for a particular block
#...
#...Arguments
#...=========
#... $block       Intent: in
#...              Perl type: scalar
#...              the block you want
#...
#...Program Flow
#...============
#... 1) return value (default of cond="" is returned if no condition set
#............................................................................
sub datafile_getblock
  {
    my(
       $block
      ) = @_;
    my(
       $ierr,
       @nothing
      );
    if( defined($$DEFAULT_DATA_FILE_REF{blocks}{$block}) ){
        return( @{$$DEFAULT_DATA_FILE_REF{blocks}{$block}} );
    }
    else{
        return( @nothing );
    }
  }
#............................................................................
#...Name
#...====
#... datafile_getval
#...
#...Purpose
#...=======
#... Returns the value of an item in the datafile with the conditions
#... set in datafile_setcond
#...
#...Arguments
#...=========
#... $key         Intent: in
#...              Perl type: scalar
#...              the key to the value you want
#...
#...Program Flow
#...============
#... 1) return value
#............................................................................
sub datafile_getval
  {
    my(
       $key,
       $error,
      ) = @_;
    my(
       $cond,
       $ierr,
       $mid,
       $mid_orig,
       $pst,
       $pre,
       $var,
      );
    if( defined($error) && ! defined($$DEFAULT_DATA_FILE_REF{key_val}{"$key"}))
      {
        $ierr = 1;
        &print_error( "Datafile variable not defined [$key]",
                      "Available variables:",
                      sort( keys( %{$$DEFAULT_DATA_FILE_REF{key_val}}) ),
                      $ierr );
        exit( $ierr );
      }
    # if the variable is defined at all
    if( defined($$DEFAULT_DATA_FILE_REF{key_val}{"$key"}) ){
        # go through each condition of the variable in order
        foreach $cond ( @{$$DEFAULT_DATA_FILE_REF{key_val}{$key}{order}} ){
            my( $cond_eval ) = $cond;
            # go through each condition set and replace that string with true (1)
            # if that condition is set, or false (0) if not.
            # First will use ";" for "true", then replace all other strings with
            # false (0), then replace ";" with true (1)
            # replace pure numbers with true (right now, ";")
            $cond_eval =~ s/\b[1-9]\d*/\;/g;
            # replace each condition set with true (;)
            foreach my $cond_set ( split( /\s*,\s*/, $COND ) ){
                $cond_eval =~ s/\b$cond_set\b/\;/g;
            }
            # replace all other words with false (0)
            $cond_eval =~ s/\b\w+\b/0/g;
            # set ";" back to true (1)
            $cond_eval =~ s/\;/1/g;
            # if non-null, then eval it
            if( $cond_eval =~ /\S/ ){
                eval( "\$cond_eval = $cond_eval" );
            }
            # default is null being true
            else{
                $cond_eval = 1;
            }
            if( $cond_eval ){
                $var = $$DEFAULT_DATA_FILE_REF{key_val}{$key}{vals}{$cond};
            }
        }
    }
    # replace any "${<var>}" with values
    if( defined( $var ) ){
        $var = &datafile_replace_var($var);
    }
    return( $var );
  }
# replace a${c}b with correct (using $ENV{c} if needed)
sub datafile_replace_var{
    my( $var ) = @_;
    my( 
        $ierr,
        $mid,
        $mid_orig,
        $pre,
        $pst,
        $var_new,
        );
    $var_new = $var;
    while ( $var_new =~ /^(.*)\$\{(\w+)\}(.*)$/ ){
        $pre = $1;
        $mid_orig = $2;
        $pst = $3;
        $mid = &datafile_getval("$mid_orig");
        if( ! defined($mid) ){
            if( defined($ENV{$mid_orig}) ){
                $mid = $ENV{$mid_orig};
            }
            else{
                $ierr = 1;
                &print_error( "Cannot find value for datafile variable [\$$mid_orig] from [$var]",
                              $ierr );
                exit( $ierr );
            }
        }
        $var_new = "${pre}${mid}${pst}";
    }
    return( $var_new );
}
#............................................................................
#...Name
#...====
#... extrema
#...
#...Purpose
#...=======
#... Given a set of x values, finds the extrema (mins and maxs)
#...
#...Arguments
#...=========
#... Y         Intent: in
#...           Perl type: ref to array
#...           The array to find min/maxs
#...
#... SD        Intent: in
#...           Perl type: ref to array
#...           Optionsl
#...           The standard deviation at all points.
#...           If supplied, the new sd's will be here
#...
#... MINS      Intent: in
#...           Perl type: ref to array
#...           Array of indices that are mins
#...
#... MAXS      Intent: in
#...           Perl type: ref to array
#...           Array of indices that are maxs
#...
#...Program Flow
#...============
#... 1) set condition
#............................................................................
sub extrema{
  my %args = (
              X         => undef,
              Y         => undef,
              SD        => undef,
              SD_FIND   => undef,
              MINS      => undef,
              MAXS      => undef,
              @_,
             );
  my $args_valid =("X|Y|SD|SD_FIND|MINS|MAXS");
  my(
     $Y_ref,
     $MAXS_ref,
     $MINS_ref,
     $SD_ref,
    );
  my(
     $arg,
     $found_extrema,
     $i,
     $ierr,
     $index,
     $index_start,
     $j,
     $max,
     $max_last,
     $min,
     $min_last,
     $next,
     $num_sd,
     $points,
     @sd,
     $sd_factor,
     $sd_window,
     $start,
     $val,
     $val_sd,
     $val_sum,
     $val_mean,
    );
  # fitting vars
  $sd_window = 10;
  $sd_factor = 4;
  $sd_num = 2;
  # args
  foreach $arg ( keys %args ){
    if( $arg !~ /^$args_valid$/ ){
      $ierr = 1;
      &print_error( "Invalid argument [$arg]",
                    "Valid args [$args_valid]",
                    $ierr );
      exit( $ierr );
    }
  }
  if( ! defined $args{Y} || ref($args{Y}) ne "ARRAY" )
    {
      $ierr = 1;
      &print_error( "Must supply Y array",
                    $ierr );
      exit( $ierr );
    }
  $Y_ref = $args{Y};
  $points = $#{$Y_ref} + 1;
  $start = 0;
  for( $i = 0; $i <= $points; $i++ ) {
    if( defined( $$Y_ref[$i] ) ) {
      $start = $i;
      last;
    }
  }
  # use supplied SD array or iterate and estimate SD by calling with 0 SD
  if( defined($args{SD}) && ref($args{SD}) ne "ARRAY" ){
    $ierr = 1;
    &print_error( "Standard deviation (SD) must be an array",
                  $ierr );
    exit( $ierr );
  }
  #...init mins/maxs
  $MINS_ref = $args{MINS};
  @$MINS_ref = ();
  $MAXS_ref = $args{MAXS};
  @$MAXS_ref = ();
  if( ! defined($args{SD}) || $#{$args{SD}} != $#{$args{Y}} ){
    @{$args{SD}} = (0)x${points};
    for( $i = 0; $i < $sd_num; $i++ ) {
      &extrema(Y=>$args{Y}, SD=>$args{SD}, SD_FIND=>"",
               MINS=>$args{MINS}, MAXS=>$args{MAXS});
      @$MINS_ref = ();
      @$MAXS_ref = ();
    }
  }
  # short vars for args
  $SD_ref = $args{SD};
  # init previous max/min
  $max_last = $start;
  $min_last = $start;
  $max = $$Y_ref[$max_last];
  $min = $$Y_ref[$min_last];
  # initial direction (should derive this?  but just defines end point)...
  $next = "max";
  # not found yet
  $found_extrema = -1;
  #...find extrema
  for( $i = $start + 1; $i < $points; $i++ )
    {
      $val = $$Y_ref[$i];
      if( ! defined( $val ) ){
        $ierr = 1;
        &print_error( "Missing data for finding extrema",
                      $ierr );
        exit( $ierr );
      }
      if( $next eq "max" ){
        if( $val > $max ){
          $max = $val;
          $max_last = $i;
        }
        elsif( $val < $max - $sd_factor*$$SD_ref[$i] ){
          $found_extrema = $i;
          push( @{$MAXS_ref}, $max_last );
          $min_last = $i;
          $min = $val;
          $next = "min";
        }
      }
      elsif( $next eq "min" ){
        if( $val < $min ){
          $min = $val;
          $min_last = $i;
        }
        elsif( $val > $min + $sd_factor*$$SD_ref[$i] ){
          $found_extrema = $i;
          push( @{$MINS_ref}, $min_last );
          $max_last = $i;
          $max = $val;
          $next = "max";
        }
      }
      # compute SD around that extrema and use that SD for points
      # since last extrema window
      if( defined( $args{SD_FIND} ) && $found_extrema != -1 ){
        $found_extrema = -1;
        $val_sum = 0;
        $val_x2 = 0;
        $index_start = $i - int($sd_window/2);
        if( $index_start < 0 ){
          $index_start = 0;
        }
        if( $index_start > $points - $sd_window ){
          $index_start = $points - $sd_window;
        }
        for( $j = 0; $j < $sd_window; $j++ )
          {
            $index = $index_start + $j;
            $val_sum += $$Y_ref[$index];
            $val_x2 += ($$Y_ref[$index])**2;
          }
        $val_mean = $val_sum / $sd_window;
        # biased seems to fit better
        $val_sd = sqrt($val_x2/$sd_window - $val_mean**2);
        for( $j = 0; $j < $sd_window; $j++ ){
          $index = $index_start + $j;
          $sd[$j] = $val_sd;
        }
      }
    }
  #...finish off standard deviation
  if( defined($args{SD_FIND}) ){
    # set standard deviation to be constant.  not sure if this is
    # best thing to do...a point based sd seems better, but too much
    # variation.
    $val_sd = 0;
    $num_sd = 0;
    for( $i = 0; $i < $points; $i++ ){
      if( defined($sd[$i]) ){
        $num_sd++;
        $val_sd += $sd[$i];
      }
    }
    if( $num_sd > 0 ) {
      $val_sd /= $num_sd;
    }
    # if no extrema, currently set SD to 0 - could do something else to
    # estimate standard deviation better.
    else {
      $val_sd = 0;
    }
    @$SD_ref = ($val_sd)x${points};
  }
}
#............................................................................
#...Name
#...====
#... datafile_setcond
#...
#...Purpose
#...=======
#... sets the condition used for getting values
#...
#...Arguments
#...=========
#... $cond     Intent: in
#...           Perl type: scalar
#...           the condition
#...
#...Program Flow
#...============
#... 1) set condition
#............................................................................
sub datafile_setcond
  {
      my( $cond ) = @_;
      $COND = $cond;
  }
#............................................................................
#...Name
#...====
#... get_sysinfo
#...
#...Purpose
#...=======
#... Populates a hash with info
#...
#...Arguments
#...=========
#... $sysinfo_ref         Intent: inout
#...                      Perl type: reference to hash
#...
#............................................................................
sub get_sysinfo{
    my(
       $sysinfo_ref,
       ) = @_;
    my(
       $ex,
       $file,
       $file_find,
       $line,
       $output,
       );
    # Makefile that gives system info
    $file_find = "Makefile_vars.mk";
    $file = &which_exec( "$file_find", QUIET=>"", NOEXEC=>"" );
    if( $file !~ /\S/ ){
        $ierr = 1;
        &print_error( "Cannot find sysinfo file [$file_find]",
                      $ierr );
        exit( $ierr );
    }
    $ex = &which_exec( "gmake" );
    if( $ex !~ /\S/ ){
        $ex = &which_exec( "make" );
    }
    if( $ex !~ /\S/ ){
        $ierr = 1;
        &print_error( "Cannot find make command [tried gmake, make]",
                      $ierr );
        exit( $ierr );
    }
    $command = "$ex -f $file print_makefile_vars";
    $output = &run_command( COMMAND=>$command );
    foreach $line ( split( /\n/, $output ) ){
        if( $line =~ /(\S+)=(\S+)$/ ){
            $$sysinfo_ref{$1} = $2;
        }
    }
}
#............................................................................
#...Name
#...====
#... my_copy
#...
#...Purpose
#...=======
#...
#...Arguments
#...=========
#... $path_from  Intent: in
#...             Perl type: scalar
#...             the input file or directory (cp -R used if directory)
#...
#... $path_to    Intent: in
#...             Perl type: scalar
#...             the output path
#...
#...Program Flow
#...============
#... 1) find last parent that exists (and get info about that parent)
#... 2) create child dirs
#............................................................................
sub my_copy
  {
    my(
       $path_from,
       $path_to,
       $group_to
      ) = @_;
    my(
       $command,
       $cwd,
       $executable,
       $final_path_to,
       $group,
       $ierr,
       $mode_dir,
       $mode_exec,
       $mode_file,
       $notdir_from,
       $parent_to,
       $type_from,
       $type_to,
      );
    $cwd = &cwd();
    if( ! -e $path_from )
      {
        $ierr = 1;
        &print_error( "my_copy: path_from [$path_from] does not exist",
                      "cwd = $cwd",
                      $ierr );
        exit( $ierr );
      }
    if( -d $path_from  )
      {
        $type_from = "directory";
      }
    else
      {
        $type_from = "file";
      }
    if( -x $path_from )
      {
        $executable = "true";
      }
    else
      {
        $executable = "false";
      }
    if( -d $path_to )
      {
        $type_to = "directory";
      }
    elsif( -e $path_to )
      {
        $type_to = "file";
      }
    else
      {
        $type_to = "";
      }
    #...get notdir
    ($notdir_from = $path_from) =~ s&/*$&&;
    if( $notdir_from =~ m&([^/]+)$& )
      {
        $notdir_from = $1;
      }
    else
      {
        $ierr = 1;
        &print_error( "my_copy: could not parse path_from [$path_from]",
                      $ierr );
        exit $ierr;
      }
    #...group
    if( ! defined($group_to) ){
      if ( $type_to eq "directory" ) {
        $group = (stat $path_to)[5];
        $group = (getgrgid($group))[0];
      } elsif ( $type_to eq "file" || $type_to eq "" ) {
        ($parent_to = $path_to) =~ s&[^/]*$&&;
        if ( $parent_to !~ /\S/ ) {
          $parent_to = ".";
        }
        $group = (stat $parent_to)[5];
        $group = (getgrgid($group))[0];
      }
    }
    else{
      $group = $group_to;
    }
    #...final_path_to
    #...file to
    if( $type_from eq "file" )
      {
        #...file to directory
        if( $type_to eq "directory" )
          {
            $final_path_to = "$path_to/$notdir_from";
          }
        #...file to file
        else
          {
            $final_path_to = "$path_to";
          }
      }
    #...directory to
    else
      {
        #...directory to existing directory
        if( $type_to eq "directory" )
          {
            $final_path_to = "$path_to/$notdir_from";
          }
        #...directory to existing non-directory
        elsif( $type_to eq "" )
          {
            $final_path_to = "$path_to";
          }
        #...directory to file
        else
          {
            $ierr = 1;
            &print_error( "my_copy: trying to copy directory to file",
                          "[$path_from] -> [$path_to]",
                          $ierr );
            exit( $ierr );
          }
      }
    #...copy
    $command = "\\cp -f -L -R $path_from $path_to";
    &run_command( COMMAND=>$command, ERROR_REGEXP=>'/\S/' );
    #...chgrp
    $command = "chgrp -R $group $final_path_to";
    &run_command( COMMAND=>$command );
    #...set permissions
    $command = "chmod -R g+rw $final_path_to";
    &run_command( COMMAND=>$command );
    #...add group execute if owner has execute
    if( -x $final_path_to ){
        $command = "chmod -R g+x $final_path_to";
        &run_command( COMMAND=>$command );
    }
    $command = "find $final_path_to -type d -exec chmod 02770  {} \\;";
    &run_command( COMMAND=>$command );
  }
#............................................................................
#...Name
#...====
#... my_copy_hpss
#...
#...Purpose
#...=======
#... creates a directory (with parents) of the group and mode
#...
#...Arguments
#...=========
#... $path_from   Intent: in
#...              Perl type: scalar
#...              the directory or file
#...
#... $path_to     Intent: in
#...               Perl type: scalar
#...               the directory or file
#...
#... $group     Intent: in
#...            Perl type: scalar
#...            the group (default is group of last parent that already exists
#...
#...Program Flow
#...============
#... 1) copy
#............................................................................
sub my_copy_hpss{
  my %args = (
              PATH_FROM  => undef,
              PATH_TO    => undef,
              GROUP      => undef,
              RM         => undef,
              ERROR_FILE => undef,
              PRESERVE   => undef,
              VERBOSE    => undef,
              @_,
             );
  my $args_valid =("ERROR_FILE|PATH_FROM|PATH_TO|GROUP|PRESERVE|RM|VERBOSE");
  my(
     $arg,
     $command,
     $ierr,
     $mode,
     $output,
     $opts,
     $PSI,
    );
  $PSI = &which_exec( "psi", ERROR=>1 );
  foreach $arg ( keys %args ){
    if( $arg !~ /^$args_valid$/ ){
      $ierr = 1;
      &print_error( "Invalid argument [$arg]",
                    "Valid args [$args_valid]",
                    $ierr );
      exit( $ierr );
    }
  }
  $mode = "2770";
  if( ! defined($args{PATH_FROM}) || $args{PATH_FROM} !~ /\S/ ||
      ! defined($args{PATH_TO})   || $args{PATH_TO}   !~ /\S/ )
    {
      $ierr = 1;
      &print_error( "path_in and/or path_out not defined or empty", $ierr );
      exit( $ierr );
    }
  $opts = "";
  if( defined($args{RM}) ){
    $opts .= " --rm";
  }
  $command = "$PSI store $opts -R $args{PATH_FROM}:$args{PATH_TO}";
  &run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/', ERROR_FILE=>$args{ERROR_FILE}, VERBOSE=>$args{VERBOSE} );
  # seems like permissions in hpss do the right thing...so might
  # not need to do any fancy about this
  #$command = "$PSI chmod -R ${mode} ${final_path_to}";
  #&run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/' );
  if( defined($args{GROUP}) ){
    $command = "$PSI chgrp $args{GROUP} $args{PATH_TO}";
    &run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/', ERROR_FILE=>$args{ERROR_FILE}, VERBOSE=>$args{VERBOSE}  );
  }
  # PRESERVE: keep same group and permissions
  if( defined($args{PRESERVE}) ){
      # get list of files
      if( -d $args{PATH_FROM} ){
          $dir = $args{PATH_FROM};
      }
      else{
          $dir = ".";
      }
      $output = `cd $dir && find . -print`;
      @files = split( /\n/, $output );
      foreach $file ( @files ){
          &my_stat( "$dir/$file", \%stat );
          $command = "$PSI chgrp -d $args{PATH_TO} $stat{group} $file";
          &run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/', ERROR_FILE=>$args{ERROR_FILE}, VERBOSE=>$args{VERBOSE}  );
          $command = "$PSI chmod -d $args{PATH_TO} $stat{mode_d} $file";
          &run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/', ERROR_FILE=>$args{ERROR_FILE}, VERBOSE=>$args{VERBOSE}  );
      }
  }
  $command = "$PSI chacl -c $args{PATH_TO}";
  &run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/', ERROR_FILE=>$args{ERROR_FILE}, VERBOSE=>$args{VERBOSE}  );
}
#............................................................................
#...Name
#...====
#... my_mkdir
#...
#...Purpose
#...=======
#... creates a directory (with parents) of the group and mode
#...
#...Arguments
#...=========
#... $dir_in   Intent: in
#...           Perl type: scalar
#...           the directory
#...
#... $group_in Intent: in
#...           Perl type: scalar
#...           the group (default is group of last parent that already exists
#...
#...Program Flow
#...============
#... 1) find last parent that exists (and get info about that parent)
#... 2) create child dirs
#............................................................................
sub my_mkdir
  {
    my(
       $dir_in,
       $group_in,
       $mode_in
      ) = @_;
    my(
       $dir_cur,
       @dirs,
       $done,
       @fields,
       $group,
       $group_dir_cur,
       $ierr,
       $mode,
       $mode_string,
       $umask_old,
      );
    #...if directory doesn't already exist
    $dir = $dir_in;
    if( ! -d $dir )
      {
        #...set default mode
        if( defined($mode_in) && $mode_in =~ /\S/ )
          {
            $mode = $mode_in;
          }
        else
          {
            $mode = 02770;
          }
        #...build up directory path of already existing directories
        #...and find info about last existing parent
        if( $dir =~ m&^/& )
          {
            $dir_cur = "/";
          }
        else
          {
            $dir_cur = "./";
          }
        $dir =~ s/^(\.\/)+//;
        $dir =~ s&^/+&&;
        @dirs = split( m&/&, $dir );
        $done = "false";
        while( $done eq "false" )
          {
            #...get info about this dir
            if( -d $dir_cur )
              {
                @fields = stat $dir_cur;
                $group_dir_cur = $fields[5];
                $group_dir_cur = (getgrgid($group_dir_cur))[0];
              }
            #...exit out - this directory on down needs to be created
            else
              {
                last;
              }
            #...if still dir children to create
            if( $#dirs >= 0 )
              {
                $dir_cur = "${dir_cur}".shift(@dirs)."/";
              }
            else
              {
                $done = "true";
                last;
              }
          }
        #...set group for rest of dirs
        if( defined( $group_in ) && $group_in =~ /\S/ )
          {
            $group = $group_in;
          }
        else
          {
            $group = $group_dir_cur;
          }
        #...ignore umask (reset when done)
        $umask_old = umask(0);
        $done = "false";
        #...create remaining parent and dir tree
        while( $done eq "false" )
          {
            $mode_string = sprintf( "%o", $mode );
            $ierr = mkdir( $dir_cur, $mode );
            if( $ierr != 1 )
              {
                $ierr = 1;
                &print_error( "mkdir failed",
                              "dir     [$dir_in]",
                              "dir_cur [$dir_cur]",
                              "mode    [$mode_string]",
                              $!,
                              $ierr );
                exit( $ierr );
              }
            # group should be defined unless getgrgid fails
            # on nfs dir where gid does not exist.
            if( defined($group) ){
                &run_command( COMMAND=>"chgrp $group $dir_cur" );
            }
            chmod( $mode, $dir_cur );
            if( $#dirs >= 0 )
              {
                $dir_cur = "${dir_cur}".shift(@dirs)."/";
              }
            else
              {
                $done = "true";
                last;
              }
          }
        #...reset umask
        umask($umask_old);
      } #...done: if directory doesn't exist already
  }
#............................................................................
#...Name
#...====
#... my_mkdir_hpss
#...
#...Purpose
#...=======
#... creates a directory (with parents) of the group and mode
#...
#...Arguments
#...=========
#... $dir_in   Intent: in
#...           Perl type: scalar
#...           the directory
#...
#... $group_in Intent: in
#...           Perl type: scalar
#...           the group (default is group of last parent that already exists
#...
#... $mode_in  Intent: in
#...           Perl type: scalar
#...           umask is always ignored
#...           the mode (default is 2770)
#...
#...Program Flow
#...============
#... 1) find last parent that exists (and get info about that parent)
#... 2) create child dirs
#............................................................................
sub my_mkdir_hpss
  {
    my(
       $dir,
       $group,
      ) = @_;
    my(
       $command,
       $ierr,
       $mode,
       $output,
       $PSI,
      );
    $PSI = &which_exec( "psi", ERROR=>1 );
    $mode = "2770";
    if( ! defined($dir) || $dir !~ /\S/ || ! defined($group) || $group !~ /\S/ )
      {
        $ierr = 1;
        &print_error( "dir and/or group not defined or empty", $ierr );
        exit( $ierr );
      }
    $command = "$PSI mkdir -p --cond ${dir}";
    &run_command( COMMAND=>$command, ERROR_REGEXP=>'/Error E/' );
    $command = "$PSI chacl -c ${dir}";
    &run_command( COMMAND=>$command );
    $command = "$PSI chgrp ${group} ${dir}";
    &run_command( COMMAND=>$command );
    $command = "$PSI chmod 2770 ${dir}";
    &run_command( COMMAND=>$command );
  }
#............................................................................
#...Name
#...====
#... print_error
#...
#...Purpose
#...=======
#... Print a standard error message.
#...
#...Arguments
#...=========
#... error_lines  Intent: in
#...              Perl type: array
#...              0: cause of error (file_name:line_number)
#...              1: explanation/fix (if there is one)
#...              2: error value
#...
#...Program Flow
#...============
#... 1) see if warning or error (last argument is 0 or not)
#... 2) Line up error lines by column
#... 3) find out who was calling this routine and the line number
#... 4) print out info
#............................................................................
sub print_error
  {
    my( 
       @error_lines # what emitted the error
      ) = @_;
    my(
       @c_filename,    # caller val
       @c_line,        # caller val
       @c_package,     # caller val
       @c_subname,     # caller val
       $error_level,   # what is printed (eg WARNING or ERROR)
       $error_line,    # each line of input argument
       $error_message, # the message to print
       $i,             # loop var
       @routine_info,  # var from caller
       $routine_name,  # routine name
       $spaces         # spaces for lining up columns
      );
    #......................................................
    #...assign WARNING or ERROR depending on error value...
    #......................................................
    if ( "$error_lines[$#error_lines]" eq "0" )
      {
        $error_level = "**WARNING**";
        $spaces      = "           ";
      }
    else
      {
        $error_level = "**ERROR**";
        $spaces      = "         ";
      }
    #............................................................
    #...DONE: assign WARNING or ERROR depending on error value...
    #............................................................
    #.......................................
    #...init error and add argument lines...
    #.......................................
    $error_message = "\n$error_level Message:\n";
    foreach $error_line ( @error_lines )
      {
        $error_message .= "$spaces  $error_line\n";
      }
    # date
    $error_message .= $spaces." Date: ".&date_ymdhms_sep()."\n";
    #...........
    #...stack...
    #...........
    $error_message .= $spaces." Stack:\n";
    $i = 0;
    @routine_info = caller($i);
    while( $#routine_info >= 3 )
      {
        $i++;
        push( @c_package,  $routine_info[0] );
        push( @c_filename, $routine_info[1] );
        push( @c_line,     $routine_info[2] );
        push( @c_subname,  $routine_info[3] );
        @routine_info = caller($i);
      }
    shift( @c_subname );
    push( @c_subname, "main" );
    for( $i = $#c_package; $i >= 0; $i-- )
      {
        $error_message .= sprintf( "%s  %04d %s%s:%s %s\n",
                                   $spaces, $#c_package - $i,
                                   " "x($#c_package - $i),
                                   $c_filename[$i], $c_line[$i],
                                   $c_subname[$i]);
      }
    #..........................................
    #...print error and return error message...
    #..........................................
    print STDERR $error_message;
    $error_message;
  }
#............................................................................
#...Name
#...====
#... print_perl_obj
#...
#...Purpose
#...=======
#... Print the structure of a perl object...for easier debugging.
#... Mixtures of arrays, hashes, and scalars can be printed.
#...
#... Up to a certain number of hash/array values are printed.
#... Beyond this, values are skipped (unless they are a reference
#... themselves) until the last value of the array/hash.
#...
#...Arguments
#...=========
#... $obj_ref          Intent: in
#...                   Perl type: reference to array
#...                   Reference to object to print (\%, \@, \$)
#...
#... $in_pref          Intent: in
#...                   Perl type: scalar
#...                   The preface string to name this object.
#...                   If passing \%foo, a good pref might be "foo".
#...                   if not defined, a default will be used.
#...
#... $in_max_items     Intent: in
#...                   Perl type: scalar
#...                   Maximum number of items in array/hash to print.
#...                   If negative, print all items.
#...                   If not defined, all items will be printed.
#...
#... $in_file          Intent: in
#...                   Perl type: scalar
#...                   Output file (usually just STDOUT).
#...                   on the same plot.
#...                   If not defined, STDOUT will be used.
#...
#...Program Flow
#...============
#... 1) If hash or array, recursively call this routine.
#... 2) If scalar, print the value.
#............................................................................
sub print_perl_obj
{
  my(
     $obj_ref,
     $in_pref,
     $in_max_items,
     $in_file
    ) = @_;
  my(
     $file, # file to use
     $i, # loop var
     @indices, # indices of object
     $j, # loop var
     $key, # key
     $max_items, # max items to use
     $new_pref, # new preface
     $pref, # preface to use
     $skip, # if skipping
     $skip_count, # how many skipped
     $val, # scalar value
    );
  #..........................
  #...fix non-defined vals...
  #..........................
  if( ! defined( $in_pref ) )
    {
      $pref = "var";
    }
  else
    {
      $pref = $in_pref;
    }
  if( ! defined( $in_max_items ) )
    {
      $max_items = -1;
    }
  else
    {
      $max_items = $in_max_items;
    }
  if( ! defined( $in_file ) )
    {
      $file = STDOUT;
    }
  else
    {
      $file = $in_file;
    }
  #...............................
  #...print depending upon type...
  #...............................
  if( ref( $obj_ref ) eq "ARRAY" )
  {
    if( $#{$obj_ref} < 0 )
      {
        print {$file} "$pref\[\] undefined\n";
      }
    else
      {
        $skip = 0;
        $skip_count = 0;
        for( $i = 0; $i <= $#{$obj_ref}; $i++ )
          {
            #...determine if should skip printing this one...
            if( $i > $max_items-1 && $max_items >= 0 )
              {
                $skip = 1;
              }
            $new_pref = sprintf( "%s[%s]", $pref, $i );
            $val = $$obj_ref[$i];
            if( ref( $val ) )
              {
                # recursively only do if non-skip
                if( $skip == 0 || $i == $#{$obj_ref} ) {
                    if( $skip_count > 0 ){
                        print {$file} "$pref\[...$skip_count\]\n";
                    }
                    &print_perl_obj( $val, $new_pref, $max_items, $file )
                }
                else{
                    $skip_count++;
                }
              }
            else
              {
                #...print last one no matter what...
                if( $skip && $i < $#{$obj_ref} )
                  {
                    $skip_count++;
                    next;
                  }
                if( $skip_count > 0 )
                  {
                    print {$file} "$pref\[...$skip_count\]\n";
                    $skip_count = 0;
                  }
                &print_perl_obj( \$val, $new_pref, $max_items, $file )
              }
          }
      }
  }
  elsif( ref( $obj_ref ) eq "HASH" )
  {
    @indices = sort keys %{$obj_ref};
    if( $#indices < 0 )
      {
        print {$file} "$pref\{\} undefined\n";
      }
    else
      {
        for( $i = 0; $i <= $#indices; $i++ )
          {
            $new_pref = sprintf( "%s{%s}", $pref, $indices[$i] );
            $val = $$obj_ref{$indices[$i]};
            if( ref( $val ) )
              {
                &print_perl_obj( $val, $new_pref, $max_items, $file )
              }
            else
              {
                &print_perl_obj( \$val, $new_pref, $max_items, $file )
              }
          }
      }
  }
  else
  {
     if( defined( $$obj_ref ) )
     {
       print {$file} "$pref = [$$obj_ref]\n";
     }
     else
     {
       print {$file} "$pref undefined\n";
     }
  }
}
#............................................................................
#...Name
#...====
#... run_command
#...
#...Purpose
#...=======
#... Runs a command and does varios things depending on args.
#... Returns the output.
#...
#...Arguments
#...=========
#... APPEND       Intent: in
#...              Perl type: scalar
#...              Default: overwrite
#...              If output is appended to ouput file
#...
#... COMMAND      Intent: in
#...              Perl type: scalar
#...              command to run
#...
#... OUT_FILE     Intent: in
#...              Perl type: scalar
#...              Default: no output file
#...              output file name
#...              If specified:
#...                 <OUT_FILE> = output file of command as it is running
#...                 <OUT_FILE>_post = output plus some additional info
#...
#... ERROR_REGEXP Intent: in
#...              Perl type: scalar
#...              Default: no error checking
#...              Regular expression to check on for error condition.
#...              Must be put in single ticks:
#...                  ERROR_REGEXP => '/^error(s)?$/i'
#...
#...Usage
#...=====
#...  use my_utils qw ( run_command );
#...  &run_command( COMMAND => "ls", OUT_FILE => './ls_out.txt', 
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub run_command
  {
    my %args = (
                APPEND       => undef,
                COMMAND      => undef,
                DEBUG        => undef,
                ERROR_FILE   => undef,
                ERROR_REGEXP => undef,
                GROUP        => undef,
                OUT_FILE     => undef,
                STATUS       => undef, # ref to status return
                STDOUT       => undef,
                TIMING       => undef,
                VERBOSE      => undef,
                TIMEOUT      => undef, # seconds to timeout
                @_,
               );
    my(
       $append_on,
       $clear_out_file,
       $cwd,
       $command_exec,
       $date,
       $debug_on,
       $ierr,
       $key,
       $out_file,
       $out_file_running,
       $out_name_in,
       $output,
       $print_error_msg,
       $status,
       $time_a,
       $time_b,
       $timeout,
       $tee,
      );
    $status = 0;
    if( defined($args{TIMING}) ){
      $time_a = time();
      print      "Running:          $args{COMMAND}\n";
    }
    #...invalid args
    foreach $key ( keys %args )
      {
        if( $key !~ /^(APPEND|COMMAND|DEBUG|ERROR_FILE|ERROR_REGEXP|GROUP|OUT_FILE|STATUS|STDOUT|TIMEOUT|TIMING|VERBOSE)$/)
          {
            $ierr = 1;
            &print_error( "Invalid argument to run_command [$key]",
                          $ierr );
            exit( $ierr );
          }
      }
    if( ! defined( $args{COMMAND} ) || $args{COMMAND} !~ /\S/ )
      {
        $ierr = 1;
        &print_error( "COMMAND not defined or empty",
                      $ierr );
        exit( $ierr );
      }
    if( defined( $args{OUT_FILE} ) && $args{OUT_FILE} !~ /\S/ )
      {
        $ierr = 1;
        &print_error( "OUT_FILE set but empty",
                      $ierr );
        exit( $ierr );
      }
    if( defined( $args{APPEND} ) && ! defined ($args{OUT_FILE} ) )
      {
        $ierr = 1;
        &print_error( "APPEND set but OUT_FILE not",
                      $ierr );
        exit( $ierr );
      }
    $timeout = "";
    if( defined($args{TIMEOUT}) ){
        $timeout = &which_exec("timeout", QUIET=>"");
        if( $timeout ne "" ){
            $timeout = "$timeout -s 9 $args{TIMEOUT}";
        }
    }
    $cwd = &cwd();
    #...if set to debug, just echo
    if( defined( $args{DEBUG} ) )
      {
        $debug_on = "true";
        $command_exec = "echo '$args{COMMAND}'";
      }
    else
      {
        $debug_on = "false";
        $command_exec = "$args{COMMAND}";
      }
    #...clear output files and tee
    if( defined( $args{OUT_FILE} ) )
      {
        $tee = "| tee -a $args{OUT_FILE}";
        $out_file_post = "$args{OUT_FILE}_post";
        if( defined( $args{APPEND} ) )
          {
            #...clear it first time in run regardless
            if( ! defined( $G_RUN_COMMAND_SEEN{"$args{OUT_FILE}"} ) )
              {
                $clear_out_file = "true";
              }
            else
              {
                $clear_out_file = "false";
              }
          }
        else
          {
            $clear_out_file = "true";
          }
        $G_RUN_COMMAND_SEEN{"$args{OUT_FILE}"} = "";
        #...clear files if set
        if( $clear_out_file eq "true" )
          {
            if( ! open( FILE, ">$args{OUT_FILE}" ) )
              {
                $ierr = 1;
                &print_error( "Cannot open command output file [$args{OUT_FILE}]",
                              "Command: $args{COMMAND}",
                              $ierr );
                exit( $ierr );
              }
            close( FILE );
            # only use _post file if doing verbose also
            if( defined($args{VERBOSE}) ){
                if( ! open( FILE, ">$out_file_post" ) )
                {
                    $ierr = 1;
                    &print_error( "Cannot open command output file [$out_file_post]",
                                  "Command: $args{COMMAND}",
                                  $ierr );
                    exit( $ierr );
                }
                close( FILE );
                if( defined($args{GROUP}) ){
                    `chgrp $args{GROUP} $args{OUT_FILE} $out_file_post`;
                }
            }
          }
        if( defined($args{VERBOSE}) ){
            if( ! open( FILE, ">>$out_file_post" ) )
            {
                $ierr = 1;
                &print_error( "Cannot open command output file [$out_file_post]",
                              "Command: $args{COMMAND}",
                              $ierr );
                exit( $ierr );
            }
        }
      } # clear output files and tee
    else
      {
        $tee = "";
      }
    #...print header, run comand, print footer
    if( defined($args{OUT_FILE}) && defined($args{VERBOSE}) )
      {
        ($date = `date 2>&1`) =~ s/\s*$//;
        print FILE "========\n";
        print      "========\n";
        print FILE "Date:             $date\n";
        print      "Date:             $date\n";
        print FILE "CWD:              $cwd\n";
        print      "CWD:              $cwd\n";
        print FILE "Command:          $args{COMMAND}\n";
        print      "Command:          $args{COMMAND}\n";
        print FILE "Debug:            $debug_on\n";
        print      "Debug:            $debug_on\n";
        print FILE "out_file:         $args{OUT_FILE}\n";
        print      "out_file:         $args{OUT_FILE}\n";
        print FILE "out_file_post:    $args{OUT_FILE}_post\n";
        print      "out_file_post:    $args{OUT_FILE}_post\n";
      }
    if( defined($args{VERBOSE}) )
      {
        print "CWD     : $cwd\n";
        print "Command : $args{COMMAND}\n";
      }

    # use timeout if set and available
    if( $timeout ne "" ){
        $command_exec = "$timeout $command_exec";
    }

    $command_exec = "($command_exec) 2>&1 $tee";
    if( defined($args{STDOUT}) ) {
      $status = system( $command_exec );
      $output = "";
    }
    else {
      $output = `$command_exec`;
      $status = $?;
    }

    # set status if defined
    if( defined($args{STATUS}) ){
        ${$args{STATUS}} = $status;
    }

    # if VERBOSE
    if( defined($args{VERBOSE}) )
      {
        print "Output  :\n$output\n";
      }

    # if error
    if( ! defined($args{DEBUG}) && defined( $args{ERROR_REGEXP} ) && eval "\$output =~ $args{ERROR_REGEXP}" )
      {
        $ierr = 1;
        $print_error_msg =
          &print_error( "Error running command [$args{COMMAND}]",
                        "cwd [$cwd]",
                        "Output from command:",
                        $output,
                        $ierr );
        if( defined($args{ERROR_FILE}) ){
          open( MY_RUN_COMMAND_FILE, ">$args{ERROR_FILE}" );
          print MY_RUN_COMMAND_FILE $print_error_msg;
          close( MY_RUN_COMMAND_FILE );
        }
        exit( $ierr );
      }
    if( defined($args{OUT_FILE}) && defined($args{VERBOSE}) )
      {
        ($date = `date 2>&1`) =~ s/\s*$//;
        print FILE "-------\n";
        print      "-------\n";
        print FILE "$output\n";
        print " [see $args{OUT_FILE}]\n";
        print FILE "-------\n";
        print      "-------\n";
        print FILE "Date:             $date\n";
        print      "Date:             $date\n";
        print FILE "========\n";
        print      "========\n";
        close( FILE );
      }
    if( defined($args{TIMING}) ){
      $time_b = time();
      printf( "TIME:    %.2f mins\n", ($time_b - $time_a)/60.0 );
    }
    return( $output );
  }
#...case insensitive sort
sub sort_case_insensitive
  {
    my( $a, $b );
    lc($a) cmp lc($b);
  }
sub sort_unique{
    my( $array_ref ) = @_;
    my( @array, %seen );
    @array = sort grep{ ! $seen{$_}++ } @{$array_ref};
}
#............................................................................
#...Name
#...====
#... my_getval
#...
#...Purpose
#...=======
#... Gets a value from stdin
#...
#...Arguments
#...=========
#... PROMPT       Intent: in
#...              Perl type: scalar
#...              Default: "Value?"
#...              The question to ask for the value.
#...
#... DEFAULT      Intent: in
#...              Perl type: scalar
#...              Default: none
#...              default value.
#...
#... TYPE         Intent: in
#...              Perl type: scalar
#...              Default: STRING
#...              Type of value read in.
#...
#... REGEXP       Intent: in
#...              Perl type: scalar
#...              Default: none
#...              Value must also match the regexp
#...
#... DIR          Intent: in
#...              Perl type: scalar
#...              Default: none
#...              When finding files, use this as path
#...
#... VAR          Intent: inout
#...              Perl type: scalar
#...              Default: none
#...              variable to assign value to.
#...
#............................................................................
sub my_getval
  {
    my %args = (
                PROMPT   => undef,
                DEFAULT  => undef,
                TYPE     => undef,
                REGEXP   => undef,
                DIR      => undef,
                BLANK    => undef,
                VAR      => undef,
                @_,
               );
    my(
       $done,
       $ierr,
       $key,
       %time,
       $val,
       $val_dir,
       $var_ref,
      );
    $ierr = 0;
    # default values
    if( ! defined($args{VAR}) ){
        $ierr = 1;
        &print_error( "Must define VAR",
                      $ierr );
        exit( $ierr );
    }
    $var_ref = $args{VAR};
    undef( $$var_ref );
    if( ! defined($args{PROMPT}) ){
        $args{PROMPT} = "Value?";
    }
    if( ! defined($args{TYPE}) ){
        $args{TYPE} = "string";
    }
    #...invalid args
    foreach $key ( keys %args ) {
        if( $key !~ /^(BLANK|DEFAULT|DIR|PROMPT|REGEXP|TYPE|VAR)$/) {
            $ierr = 1;
            &print_error( "Invalid argument [$key]",
                          $ierr );
            exit( $ierr );
        }
    }
    # get value
    $done = "false";
    while( $done eq "false" ){
        $done = "true";
        # read in val
        print "$args{PROMPT} ";
        if( defined($args{DEFAULT}) ){
            print "(default [$args{DEFAULT}]) ";
        }
        if( defined($args{BLANK}) ){
            print "(<space> for no value) ";
        }
        $val = <STDIN>;
        chomp( $val );
        if( $val =~ /^\s+$/ && defined($args{BLANK}) ){
            $$var_ref = "";
            return( $ierr );
        }
        $val =~ s/^\s*//;
        $val =~ s/\s*$//;
        if( $val eq "" && defined($args{DEFAULT}) ){
            $val = $args{DEFAULT};
        }
        if( defined($args{DIR}) ){
            if( $val =~ m&^/& ){
                $val_dir = $val;
            }
            else{
                $val_dir = "$args{DIR}/$val";
            }
        }
        else{
            $val_dir = $val;
        }
        # make sure matches REGEXP if any
        if( defined($args{REGEXP}) && eval "\$val !~ $args{REGEXP}" ){
            print "ERROR: [$val] does not match [$args{REGEXP}]\n";
            $done = "false";
            next;
        }
        # check TYPE
        if( $args{TYPE} eq "string" ){
        }
        elsif( $args{DEFAULT} eq "" && $val eq "" ){
        }
        elsif( $args{TYPE} eq "exec" ){
            if( ! -x $val_dir || -d $val_dir ){
                print "ERROR: [$val_dir] not executable\n";
                $done = "false";
            }
        }
        elsif( $args{TYPE} eq "yes/no" ){
            if( $val ne "yes" && $val ne "no" ){
                print "ERROR: [$val] must be either 'yes' or 'no'\n";
                $done = "false";
            }
        }
        elsif( $args{TYPE} eq "file" ){
            if( ! -e $val_dir || -d $val_dir ){
                print "ERROR: [$val_dir] not file\n";
                $done = "false";
            }
        }
        elsif( $args{TYPE} eq "dir" ){
            if( ! -d $val_dir ){
                print "ERROR: [$val_dir] not directory\n";
                $done = "false";
            }
        }
        elsif( $args{TYPE} eq "int" ){
            if( $val !~ /^\d+$/ ) {
                print "ERROR: [$val] must be an int\n";
                $done = "false";
            }
        }
        elsif( $args{TYPE} eq "time" ){
            %time = &conv_time( STRING=>$val );
            $val = $time{hms};
        }
        else{
            $ierr = 1;
            &print_error( "Invalid type [$args{TYPE}]",
                          $ierr );
            exit( $ierr );
        }
    }
    $$var_ref = $val;
    return( $ierr );
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
#...      executable in question
#... QUIET Intent: in
#...       Perl type: hash scalar
#...       be quiet
#... ERROR Intent: in
#...       Perl type: hash scalar
#...       error if not found
#... NOEXEC Intent: in
#...        Perl type: hash scalar
#...        execute permission not needed for file
#...
#...Program Flow
#...============
#... 1) go through path to find exec
#............................................................................
sub which_exec {
    my(
       $exec,
       ) = shift(@_);
    my %args = (
                QUIET        => undef,
                ERROR        => undef,
                NOEXEC       => undef,
                @_,
                );
    my(
       $cwd, # current working dir
       $ierr, # error return value
       $exec_try, # see if exec is here
       $found, # if found
       $path, # current dir in search for execs
       @paths, # list of search paths for execs
       $this_dir, # current directory
       );
    #...invalid args
    foreach $key ( keys %args ) {
        if( $key !~ /^(QUIET|ERROR|NOEXEC)$/) {
            $ierr = 1;
            &print_error( "Invalid argument to which_exec [$key]",
                          $ierr );
            exit( $ierr );
        }
    }
    #.................................................
    #...build paths from PATH and current directory...
    #.................................................
    $cwd = getcwd();
    $path = $ENV{PATH};
    @paths = split( /:/, $ENV{PATH} );
    ($this_dir = $0) =~ s/\/[^\/]+$//;
    unshift( @paths, $this_dir );
    #.....................................
    #...loop through paths to find exec...
    #.....................................
    $found = "false";
    foreach $path (@paths) {
        if( $path eq "." ){
            $path = $cwd;
        }
        $exec_try = "$path/$exec";
        if( -f "$exec_try" ){
            # NOEXEC or it must also be executable
            if( defined($args{NOEXEC}) || -x "$exec_try" ){
                $found = "true";
                last;
            }
        }
    }
    #........................................
    #...error if still could not find exec...
    #........................................
    if( $found eq "false" ) {
        $ierr = 0;
        if( ! defined($args{QUIET}) || defined($args{ERROR}) ) {
            if( defined($args{ERROR}) ){
                $ierr = 1;
            }
            &print_error(
                         "Executable [$exec] not found in PATH",
                         $ierr
                         );
            if( defined($args{ERROR}) ){
                exit( $ierr );
            }
        }
        $exec_try = "";
    }
    return( $exec_try );
}
sub my_notdir
  {
    my(
       $dir,
      ) = @_;
    my(
       $notdir,
      );
    ($notdir = $dir) =~ s&/*$&&;
    if( $notdir =~ m&([^/]+)$& )
      {
        $notdir = $1;
      }
    return( $notdir );
  }
sub my_dir
  {
    my(
       $file,
      ) = @_;
    my(
       $dir,
      );
    if( $file !~ /^\// ){
        $file = "./$file";
    }
    ($dir = $file) =~ s&^(.*)/.*?$&$1&;
    return( $dir );
  }
sub my_stat
  {
    my(
       $file,
       $stat_ref,
      ) = @_;
    my(
       $cwd,
       @fields,
       $ierr,
       $time,
      );
    $ierr = 0;
    undef( %{$stat_ref} );
    $time = time;
    if( -e $file ){
        @fields = stat $file;
        $$stat_ref{dev}     = $fields[0];
        $$stat_ref{ino}     = $fields[1];
        $$stat_ref{mode}    = $fields[2];
        $$stat_ref{nlink}   = $fields[3];
        $$stat_ref{uid}     = $fields[4];
        $$stat_ref{gid}     = $fields[5];
        $$stat_ref{rdev}    = $fields[6];
        $$stat_ref{size}    = $fields[7];
        $$stat_ref{atime}   = $fields[8];
        $$stat_ref{mtime}   = $fields[9];
        $$stat_ref{ctime}   = $fields[10];
        $$stat_ref{blksize} = $fields[11];
        $$stat_ref{blocks}  = $fields[12];
        # additional info
        $$stat_ref{atime_since} = $time - $$stat_ref{atime};
        $$stat_ref{mtime_since} = $time - $$stat_ref{mtime};
        $$stat_ref{mtime_localtime} = scalar localtime($$stat_ref{mtime});
        $$stat_ref{group} = getgrgid($$stat_ref{gid});
        # link stuff
        if( -l $file ){
            $$stat_ref{slink} = readlink($file);
        }
        # directory stuff
        $$stat_ref{fullpath} = Cwd::realpath($file);
        $$stat_ref{dir} = my_dir( $$stat_ref{fullpath} );
        $$stat_ref{notdir} = my_notdir($$stat_ref{fullpath});
        $$stat_ref{mode_d} = sprintf( "%o", $$stat_ref{mode} );
    }
    else{
        $ierr = 1;
    }
    return( $ierr );
}
#............................................................................
#...Name
#...====
#... my_xml_read
#...
#...Purpose
#...=======
#... Reads a text xml file with restrictions on the format.
#... Stuffs it into a hash.
#...
#...hash format
#...===========
#...  At each level:
#...    tag            = hash to next tag level (containing same keys)
#...    tag_name       = name of the tag
#...    tag_start_line = line string of start tag
#...    tag_start_ln   = line number in file of start tag
#...    tag_stop_line  = line string of stop  tag
#...    tag_stop_ln    = line number in file of stop  tag
#...    val            = value of this tag (either val or tag)
#... Example:
#... $hash{filename}{tag}
#... $hash{filename}{tag_start_line = ""
#...
#...Arguments
#...=========
#... HASHREF Intent: inout
#...      Perl type: reference to hash
#...
#... LABEL Intent: in
#...      Perl type: string
#...      xml start/stop tags of the form <LABEL<string>> and </LABEL<string>>
#...
#... FILENAME Intent: inout
#...      Perl type: string
#...      Filename to read
#...
#...Program Flow
#...============
#... 1) 
#............................................................................
sub my_xml_read{
    my %args = (
                FILE=>"",
                HASHREF=>undef,
                LABEL=>"",
                KEY=>undef,
                @_,
                );
    my(
       $file,
       $hashref,
       $label,
       );
    my(
       $done,
       $ierr,
       $key,
       $key_string,
       @keys,
       $line,
       $line_orig,
       $ln,
       $loc_ref,
       $pre,
       $rest,
       $stop,
       $str,
       $tag,
       $tag_last,
       $tag_str,
       $tag_str_short,
       @tags,
       $val,
       $var,
       );
    #...invalid args
    foreach $key ( keys %args ) {
        if( $key !~ /^(FILE|HASHREF|KEY|LABEL)$/) {
            $ierr = 1;
            &print_error( "Invalid argument [$key]",
                          $ierr );
            exit( $ierr );
        }
    }
    $file = $args{FILE};
    $hashref = $args{HASHREF};
    $label = $args{LABEL};
    $key = $args{KEY};
    if( ! defined( $key ) ){
        $key = $file;
    }
    if( !defined($hashref) || ref($hashref) ne "HASH" ){
        $ierr = 1;
        &print_error( "HASHREF must be a reference to a hash",
                      $ierr );
        exit( $ierr );
    }
    if( ! open( FILE, "$file" ) ){
        $ierr = 1;
        &print_error( "File does not exist [$file]",
                      $ierr );
        exit( $ierr );
    }
    $ln = 0;
    $val = "";
    @tags = ("{\"$key\"}");
    $$hashref{"$key"}{tag_name}       = $key;
    $$hashref{"$key"}{tag_start_line} = "";
    $$hashref{"$key"}{tag_start_ln}   = 0;
    $$hashref{"$key"}{tag_stop_line}  = "";
    $$hashref{"$key"}{tag_stop_ln}    = 0;
    while( $line_orig=<FILE> ){
        $ln++;
        $line_orig =~ s/\n$//;
        $line = $line_orig;
        # consume line
        $done = "false";
        while( $done eq "false" ){
            # line contains start or stop
            if( $line =~ /(^.*?)<(\/)?$label([^>]*?)>(.*)$/ ){
                $pre = $1;
                $stop = $2;
                $tag = $3;
                $rest = $4;
                $val = $val.$pre;
                $tag_str = join( "{\"tag\"}", @tags );
                $tag_str_short = join( "", @tags );
                $str = "\\\%{\$\$hashref${tag_str}}";
                eval "\$loc_ref = $str";
                # ending tag - add to value and stop
                if( defined( $stop ) ){
                    $$loc_ref{"tag_stop_ln"} = $ln;
                    $$loc_ref{"tag_stop_line"} = $line_orig;
                    $tag_last = pop( @tags );
                    $tag_last =~ s/^{\"(.*)\"}/$1/;
                    if( "$tag_last" ne "$tag" ){
                        $ierr = 1;
                        &print_error( "Trying to end tag [$tag] of tag chain=[$tag_str_short]",
                                      "Start: $file:$$loc_ref{tag_start_ln}",
                                      "  $$loc_ref{tag_start_line}",
                                      "Stop:  $file:$ln",
                                      "  $line_orig",
                                      $ierr );
                        exit( $ierr );
                    }
                    $val =~ s/^\s+//;
                    $val =~ s/\s+$//;
                    # store val if non-blank
                    if( $val =~ /\S/ ){
                        # cannot have embedded values and tags
                        if( defined( $$loc_ref{"tag"} ) ){
                            @keys = sort keys(%{$$loc_ref{"tag"}} );
                            $key_string = join( ",", @keys );
                            $ierr = 1;
                            &print_error( "Cannot have both value and sub tags for tag=[$tag_str_short]",
                                          "Sub tags: $key_string",
                                          "Start: $file:$$loc_ref{tag_start_ln}",
                                          "  $$loc_ref{tag_start_line}",
                                          "Stop:  $file:$ln",
                                          "  $line_orig",
"[$val]\n",
                                          $ierr );
                            exit( $ierr );
                        }
                        $$loc_ref{"val"} = $val;
                        $$loc_ref{"val_line"} = $ln;
                        $val = "";
                    }
                }
                # starting tag
                else{
                    # if already defined, delete old one
                    delete($$loc_ref{tag}{$tag});
                    $$loc_ref{tag}{$tag}{"tag_start_ln"} = $ln;
                    $$loc_ref{tag}{$tag}{"tag_start_line"} = $line_orig;
                    $$loc_ref{tag}{$tag}{"tag_name"} = $tag;
                    # cannot have embedded values and tags
                    if( $val =~ /\S/ ){
                        @keys = sort keys(%{$$loc_ref{"tag"}} );
                        $key_string = join( ",", @keys );
                        $ierr = 1;
                        &print_error( "Cannot have both value and sub tags for tag=[$tag_str_short]",
                                      "Sub tags: $key_string",
                                      "Start: $file:$$loc_ref{tag_start_ln}",
                                      "  $$loc_ref{tag_start_line}",
                                      "Next start:  $file:$ln",
                                      "  $line_orig",
                                      $ierr );
                        exit( $ierr );
                    }
                    push( @tags, "{\"$tag\"}" );
                }
                $line = $rest;
            }
            # push onto current val
            else{
                $val = $val.$line."\n";
                $line = "";
            }
            # no more line, done
            if( $line eq "" ){
                $done = "true";
            }
        }
    }
    close( FILE );
    return( $ierr );
}
#............................................................................
#...Name
#...====
#... conv_time
#...
#...Purpose
#...=======
#... Converts different time formats and defines a consistent hash
#...
#...hash format
#...===========
#...  Input (may be blank):
#...    SECS, MINS, HOURS, DAYS, YEARS
#...    STRING: colon separated list of values (might start with an A)
#...  Output:
#...    Above Input filled out (except for STRING)
#...    SECS_TOTAL, 
#...    string
#...    string_b
#...    hms
#...
#...Program Flow
#...============
#... 1) 
#............................................................................
sub conv_time {
  my %args = (
              SECS=>   0,
              MINS=>   0,
              HOURS=>  0,
              DAYS=>   0,
              YEARS=>  0,
              STRING=> "",
              @_,
             );
  my(
     $factor,
     $hours,
     $ierr,
     %num_secs,
     %size,
     $time,
     $type,
     %types,
     $unit,
     %unit_next,
     $val_new,
     $val,
     @vals,
    );
  # some constants
  $size{SECS}  = 60;
  $size{MINS}  = 60;
  $size{HOURS} = 24;
  $size{DAYS}  = 365;
  $unit_next{SECS}  = "MINS";
  $unit_next{MINS}  = "HOURS";
  $unit_next{HOURS} = "DAYS";
  $unit_next{DAYS}  = "YEARS";
  $unit_next{YEARS} = "UNKNOWN";
  $num_secs{SECS}  = 1;
  $num_secs{MINS}  = 60;
  $num_secs{HOURS} = 60*$num_secs{MINS};
  $num_secs{DAYS}  = 24*$num_secs{HOURS};
  $num_secs{YEARS} = 365*$num_secs{DAYS};
  # parse a STRING time
  if( defined($args{STRING}) && $args{STRING} =~ /\S/ ){
      $val_new = $args{STRING};
      $args{YEARS} = 0;
      $args{DAYS}  = 0;
      $args{HOURS} = 0;
      $args{MINS}  = 0;
      $args{SECS}  = 0;
      $types{y} = "YEARS";
      $types{d} = "DAYS";
      $types{h} = "HOURS";
      $types{m} = "MINS";
      $types{s} = "SECS";
      # S, M:S, H:M:S, D:H:M:S, Y:D:H:M:S
      if( $val_new =~ /^A(.*)$/ ){
          $val_new = $1;
          $val_new =~ s/-/:/;
          @vals = split(":",$val_new);
          $unit = "SECS";
          foreach $val ( reverse @vals ){
              $args{$unit} = $val;
              $unit = $unit_next{$unit};
          }
      }
      # if just given one digit, it is in minutes
      elsif( $val_new =~ /^\s*(\d+)\s*$/ ){
          $args{MINS} = $1;
      }
      # H:, H:M, H:M:S, D:H:M:S, Y:D:H:M:S
      elsif( $val_new =~ /^\s*\d+\s*:/ ){
          $val_new =~ s/\s+//g;
          @vals = split(":",$val_new);
          if( $#vals <= 2 ){
              unshift( @vals,0,0);
          }
          elsif( $#vals == 3 ){
              unshift( @vals, 0 );
          }
          ($args{YEARS},$args{DAYS},$args{HOURS},$args{MINS},$args{SECS}) = @vals;
      }
      # D-H(:M(:S))
      elsif( $val_new =~ /\s*(\d+)-(\d+(:\d+)*)\s*$/ ){
          $args{DAYS} = $1;
          $val = $2;
          @vals = split(":", $val );
          ($args{HOURS},$args{MINS},$args{SECS}) = @vals;
      }
      else{
          while( $val_new =~ /\S/ ){
              if( $val_new =~ /^\s*(\d+)\s*(y|d|h|m|s)([a-z]*)\s*:?\s*(.*?)$/i ){
                  $val_new = $4;
                  $time = $1;
                  $type = $2;
                  $type =~ tr/A-Z/a-z/;
                  $args{$types{$type}} = $time;
              }
              else{
                  $ierr = 1;
                  &print_error( "Unrecognized time format [$args{STRING}]",
                                $ierr );
                  exit( $ierr );
              }
          }
      }
  }
  # convert to total number of seconds
  $args{SECS_TOTAL} = 0;
  foreach $unit ( "SECS", "MINS", "HOURS", "DAYS", "YEARS" ) {
      if( ! defined($args{$unit}) || $args{$unit} !~ /\d/ ){
          $args{$unit} = 0;
      }
  }
  foreach $unit ( "SECS", "MINS", "HOURS", "DAYS", "YEARS" ) {
    $args{SECS_TOTAL} += $num_secs{$unit} * $args{$unit};
  }
  # put into units
  foreach $unit ( "SECS", "MINS", "HOURS", "DAYS", "YEARS" ) {
      $args{$unit} = 0;
  }
  $args{SECS} = abs($args{SECS_TOTAL});
  foreach $unit ( "SECS", "MINS", "HOURS", "DAYS" ) {
      if( $args{$unit} >= $size{$unit} ){
          $factor = int($args{$unit} / $size{$unit});
          $args{$unit_next{$unit}} += $factor;
          $args{$unit} -= $factor * $size{$unit};
      }
  }
  # string
  $args{string} = "";
  foreach $unit ( "YEARS", "DAYS", "HOURS", "MINS" ) {
    if( $args{$unit} > 0 ) {
      $args{string} .= " $args{$unit} $unit";
    }
  }
  if( $args{SECS} > 0 )
    {
      $args{string} .= sprintf( " %.2f", $args{SECS})." SECS";
    }
  if( $args{string} eq "" )
    {
      $args{string} = " 0 SECS";
    }
  $args{string} =~ s/^\s+//;
  # condensed string
  $args{string_b} = "0_SECS";
  foreach $unit ( "YEARS", "DAYS", "HOURS", "MINS", "SECS" ) {
    if( $args{$unit} > 0 ){
      $args{string_b} = "$args{$unit}_$unit";
      last;
    }
  }
  # h:m:s
  $hours = $args{YEARS}*$size{DAYS}*$size{HOURS} + $args{DAYS}*$size{HOURS} + $args{HOURS};
  $args{hms} = sprintf( "%d:%02d:%02d", $hours, $args{MINS}, $args{SECS} );
  return( %args );
}
sub latexify{
    my( $string_in, $type ) = @_;
    my( $string );
    $string = $string_in;
    # was needed to get underscores to be underscores...but not needed with the packages:
    # \usepackage{textcomp}
    # \usepackage[T1]{fontenc}
    # does not work since used in args for other commands $string =~ s/_/\\verb1_1/g;
    # $string =~ s/_/\\url{_}/g;
    $string =~ s/_/\\_/g;
    $string =~ s/&/\\&/g;
    $string =~ s/^(\s*<br>\s*)+//mgi;
    $string =~ s/<br>/\\newline /gi;
    $string =~ s/<\/?cfoutput>//gi;
    $string =~ s/</\$<\$/g;
    $string =~ s/>/\$>\$/g;
    $string =~ s/\^/\*\*/g;
    $string =~ s/\#/number/g;
    $string =~ s/\036/-/g;
    $string =~ s/\010/ /g;
    # specific hacks
    $string =~ s/\}( typical energy)/$1/;
    $string =~ s/(matdef)\{/$1/;
    $string =~ s/(thres)\}/$1/;
    return( $string );
}
#...return true
1;
