eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

# NOTE: You can use FindBin in the script. The modules will automatically have access to $FindBin::Bin.
use FindBin qw($Bin);
use lib ("$Bin", "$Bin/lib", "$Bin/../lib");

use my_utils qw (
                 date_ymdhms
                 get_pname
                 get_sysinfo
                 my_dir
                 my_mkdir
                 my_stat
                 path_add_default
                 print_error
                 print_perl_obj
                 run_command
                );
use File::Basename;
use read_output_files qw (
                          parse_output_file
                          parse_output_file_finish
                         );
use Cwd;
my(
   $ierr,
   $keeplast,
   );
print "\nStarted: $0\n";
print `date`,"\n";
$ierr = 0;
# copied from run_job.pl
$RJ_DIR                 = "rj_adir";
$RJ_FILE_TAG            = "$RJ_DIR/rj_tag";
$RJ_FILE_TAG_OLD        = "$RJ_DIR/rj_tag_old";
$RJ_FILE_CMD_OUT        = "rj_cmd_out";
$RET{RJ_RETRY}     = 2;
$RET{RJ_NEW_CHILD} = 3;
$RET{RJ_STOP}      = 4;
# add to PATH
&path_add_default();
# get system info and stuff into cmd
&get_sysinfo( \%{$cmd{sys_info}} );
&parse_args( \@ARGV, \%cmd );
if( defined( $cmd{h} ) ) {
    print <<"EOF";
#............................................................................
#...Name
#...====
#... $0
#...   Used in conjunction with run_job.pl to:
#...     o clean up files from a run (tag/move files)
#...     o check various things from a run and spit out the appropriate message
#...       that is used by run_job.pl to see if a job should continue
#...
#...   Do NOT use this script while a problem is running since it tags/moves
#...   output files around that might still be in use.
#...
#...   When using this script from within run_job.pl, you should use it before
#...   the current run is done (see the example when running 'job_run.pl -h').
#...   This script is designed to determine what tag to use based on running
#...   it right before a new run.
#...
#...Usage
#...=====
#... $0 <options>
#...
#... [--check]
#...   Also check to see if a job should continue.  It is assumed that this is
#...   called right after a run, so the new tag file $RJ_FILE_TAG is used to
#...   move files around.
#... [--clean]
#...   Recursively remove all files and directories that were created by run_job_cleanup.pl
#...   and run_job.pl.  Currently:
#...      rj_*
#... [--checkfile <output data file>]
#...   If checking, $RJ_FILE_CMD_OUT is checked by default.
#...   File to check for results - returns RJ_STOP, RJ_NEW_CHILD, RJ_RETRY
#... [--debug] [-d]
#...   Do not actually do anything - just print what it would do.
#... [--group <group>]
#..    The environment variable RJ_VAR_GROUP is used (set by run_job.pl).
#...   Change permissions so that unix group has access to files.
#... [--in <<pname>.in(put)>]
#...   <pname>.in(put)? is used by default.
#...   input file name.
#... [--pname <problem name>]
#...   If not set:
#...     run_job.pl sets the environment variable RJ_VAR_PNAME which is used.
#...     Try to determine the problem name from the files in the current directory.
#... [--proj <$cmd{proj}=default>]
#...   project name (eap, lap, ... )
#... [--skip <class of files to skip cleaning (moving around to different dirs/renaming)>]
#...   Can be used multiple times for different things to skip
#... [--tag <tag name>], [--tagfile <tag file>], --tagtype <new|old>]
#...   Different ways to specify what tag is used for tagging files.
#...   Standalone:        default uses $RJ_FILE_TAG
#...   Within run_job.pl: default uses $RJ_FILE_TAG_OLD
#...
#...   Using within run_job.pl input command file:
#...     When run_job.pl is run, it first pushes the existing $RJ_FILE_TAG to $RJ_FILE_TAG_OLD.
#...     So, existing files from previous run that were not tagged should use $RJ_FILE_TAG_OLD.
#...     Calling run_job_cleanup.pl before the run uses the old tag file.  If you wish to
#...     do a cleanup after the run, use the new tag file (used by default when using the
#...     --check option).
#...   Using outside of run_job.pl:
#...     Now it is assumed that files to be cleaned correspond to the new tag file $RJ_FILE_TAG,
#...     which is the default when run_job_cleanup.pl is used outside of run_job.pl.
#... [--tagdir <.=default>]
#...     The directory where the tag files are.
#...
#...Return Values:
#...==============
#...  1: some fatal error from running the script
#...  $RET{RJ_RETRY}: RJ_RETRY
#...  $RET{RJ_NEW_CHILD}: RJ_NEW_CHILD
#...  $RET{RJ_STOP}: RJ_STOP
#...
#...Notes:
#...======
#...  o hcard parsing for the eap project
#...    For the eap project, hcard parsing is done and the input file is
#...    modified if needed.
#...    1) determine current simulation time from <pname>-status file
#...    2) look at each line in the input file:
#...        !h <time> <line>
#...       2.1) if <time> < <simulation time>, convert line to:
#...              <line> !h <time>
#...            This basically uncomments the line for the next run.
#...    3) If hcard replacement is detected
#...       3.1) remove the -DO_NOT_RUN file
#...       3.2) if not doing --check, also modify the input file
#...            This is so that the input file only gets modified
#...            at the beginning of a new run.  See the examples below
#...            for calling the script from inside run_job.csh
#...
#...Example
#...=======
#... $0
#...    Clean up files for the default problem name.
#...
#... $0 --check
#...    Clean up files for the default problem name and detect if the run
#...    should continue.
#...
#... $0 --skip -dmp
#...    Skip moving the '-dmp' files.
#...
#... From inside run_job.csh file:
#...     $0
#...     &&&RJ_CMD_PRUN&&& &&&RJ_VAR_PNAME&&&.input
#...     $0 --check
#...   This will clean up from any previous run (do hcard parsing if needed
#...   for the eap project - see Notes above), run the code,
#...   then check to see if finished.
#...
#... EAP example: do a dump and stop running
#...   cd <run directory>
#...   touch <problem name>-kill
#...
#... EAP example: do a dump, stop the current run but keep running
#...   (useful for when you modify the input file to dezone, but
#...    want to keep running afterwards)
#...   touch <problem name>-restart <problem name>-kill
#...   (the -restart and -kill files will be removed, problem will
#...    dump, but run again).
#...   NOTE: after version 1203.00, the code recognizes the
#...         <problem name>-restart file and only this file needs to be
#...         touched.
#...
#............................................................................
EOF
    exit;
}
# clean
if( defined($cmd{clean}) ){
    print "\nRemoving files.\n\n";
    if( ! -f "run_job.csh" ){
        $ierr = 0;
        &print_error( "Cannot find run_job.csh input file.",
                      "Are you in the correct directory?",
                      $ierr );
    }
    $command = "rm -rf $RJ_DIR/";
    &run_command( COMMAND=>"$command", VERBOSE=>"", DEBUG=>$cmd{debug} );
    $command = "rm -f rj_*";
    &run_command( COMMAND=>"$command", VERBOSE=>"", DEBUG=>$cmd{debug} );
    # if ever you do not want to exit, need to do this last since RJ_FILE_TAG is also removed
    &my_exit( $ierr );
}
# cleanup
print "Cleaning (pname = $cmd{pname}, tag = $cmd{tag}, tagfile = $cmd{tagfile}, tagtype = $cmd{tagtype}):\n";
# set permissions recursively in current directory
# do this in the background since it might take a long time on some filesystems
# parent block
if( $pid_child=fork ){
    # parent will just continue
}
elsif( defined $pid_child ) {
    if( defined($cmd{group}) ){
        $command = "chgrp -R $cmd{group} .";
        &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
        $command = "find . -type d | xargs chmod $cmd{mode_dir}";
        &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
        $command = "chmod -R g+r .";
        &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
        # fix permissions for parent directories
        &fix_dir_perms( \%cmd, "." );
    }
    exit;
}
# clean up crash files by default
print "Removing crash files (core, .localdomain.btr, ...)\n";
$command = "rm -f *.localdomain.btr";
&run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
$command = "rm -f core.[0-9]*";
&run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
# get files in this directory
if( ! opendir( THISDIR, "." ) ){
    $ierr = 1;
    &print_error( "Cannot read current directory",
                  "RJ_STOP",
                  $ierr );
    exit( $ierr );
}
@files_rundir = grep( !/^\.\.?$/, readdir( THISDIR ) );
closedir THISDIR;
# Definition of fields:
#  copy      = if "yes", copy the tag'd file over over instead of move (so original stays also)
#  dest      = destination directory: <dest>/<file list from find>
#  dir       = name of file passed to find: <dir>/<file>
#  file      = name of file passed to find: <dir>/<file>
#  keeplast  = Keep this last file around.  if set to "", keep last file alphabetically.
#  overwrite = if "no", do not copy over any existing file with same tag or name
#  tag       = new name will have a tag at end
# eap files
if( $cmd{proj} eq "eap" ){
    # do hcard parsing - grab before cleaning is done
    &eap_hcard( \%cmd );
    $cleanfile{"-editmix"} = {
        "file"      => ".*/$cmd{pname}-editmix",
        "dest"      => "edits",
    };
    $cleanfile{"-idf2dmp"} = {
        "file"      => ".*/$cmd{pname}-idf2dmp",
        "dest"      => "idfdumps",
    };
    $cleanfile{"-icmp"} = {
        "file"      => ".*/$cmd{pname}-icmp_(temp|dens)\$",
        "dest"      => "icmp",
    };
    $cleanfile{"input"} = {
        "file"      => ".*/$cmd{in}",
        "dest"      => "inputs",
        "copy"      => "yes",
        "overwrite" => "no",
    };
    $cleanfile{"-lnk"} = {
        "file"      => ".*/$cmd{pname}-lnk.*",
        "tag"       => "no",
        "dest"      => "lnkfiles",
    };
    # DIFF: keep latest one as well as tagged copy of latest one
    $cleanfile{"-output"} = {
        "file"      => ".*/$cmd{pname}-output",
        "dest"      => "outputs",
        "copy"      => "yes",
    };
    $cleanfile{"-probe"} = {
        "file"      => ".*/$cmd{pname}-probe",
        "dest"      => "outputs",
        "copy"      => "yes",
    };
    $cleanfile{"-status"} = {
        "file"      => ".*/$cmd{pname}-status",
        "dest"      => "logs",
        "copy"      => "yes",
    };
    $cleanfile{".storelog"} = {
        "file"      => ".*/$cmd{pname}.storelog",
        "dest"      => "storelogs",
    };
    # DIFF: just keep 1 old copy of tracers - this is different than old run scripts
    #       If a change is detected (what tracer vars or num, then tag the old
    #       tracer file before copying new one.
    $cleanfile{"-tracer"} = {
        "file"      => ".*/$cmd{pname}-tracer",
        "tag"       => "check",
        "dest"      => "tracers",
    };
    $cleanfile{"-xray"} = {
        "file"      => ".*/$cmd{pname}-[xX]ray.*",
        "tag"       => "no",
        "dest"      => "xrays",
    };
    $cleanfile{"perf"} = {
        "file"      => ".*/perf_data\\\..*",
        "tag"       => "no",
        "dest"      => "perf_monitor",
    };
    # setup files if after cycle 0
    undef( $cycle );
    if( -e "$cmd{pname}-lastcycle" ){
        $output = `cat $cmd{pname}-lastcycle`;
        if( $output =~ /^\s*(\d+)/ ){
            $cycle = $1;
        }
    }
    if( defined( $cycle ) && $cycle > 0 ){
        $cleanfile{"setup"} = {
            "file"      => ".*\\\.(oso|stl|stl1.*)",
            "tag"       => "no",
            "dest"      => "setup_files",
        };
    }
    # -dmp files
    # try to find location of dumps by looking for keywords in the -output file(s)
    # preference is:
    #   ./<pname>-output
    #   ./outputs/<last -output file>
    if( opendir( THISDIR, "outputs" ) ){
        @files = grep( !/^\.\.?$/, readdir( THISDIR ) );
        @files = grep( /$cmd{pname}-output/, @files);
        # last one first preference
        @files = reverse @files;
        grep( s/^/outputs\//, @files );
        closedir THISDIR;
    }
    if( -e "$cmd{pname}-output" ){
        unshift( @files, "$cmd{pname}-output" );
    }
    foreach $file ( @files ){
        $output = `grep parallel_directory $file 2>&1`;
        @lines = split( /\n/, $output );
        grep( s/[\!\#].*//, @lines ); # remove comment
        @lines = grep( /parallel_directory/, @lines ); # see what is left
        # look at last definition of parallel_directory
        if( $#lines >= 0 && $lines[-1] =~ /\bparallel_directory\s*=\s*(\S+)/ ){
            $parallel_directory = $1;
            $parallel_directory =~ s/^[\'\"]//;
            $parallel_directory =~ s/[\'\"]$//;
            $parallel_directory =~ s/\/+$//;
            last;
        }
    }
    if( ! defined($parallel_directory) ){
        $parallel_directory = ".";
    }
    $keeplast = "";
    if( -e "$cmd{pname}-lastdump" ){
        $keeplast = `cat $cmd{pname}-lastdump`;
        $keeplast =~ s/\s//g;
        # if not absolute path, prepend parallel_directory
        if( $keeplast !~ /^\// ){
            $keeplast = "$parallel_directory/$keeplast";
        }
        # for bug in filesystem with reading/writing file, see if doing a ls here syncs things
        $output = `ls -la $keeplast 2>&1`;
        print "\nLast dump file (will not be moved):\n";
        print "$output\n";
    }
    $cleanfile{"-dmp"} = {
        "file"     => ".*/$cmd{pname}-dmp[0-9]+.*",
        "tag"      => "no",
        "keeplast" => $keeplast,
        "dir"      => "$parallel_directory",
        "dest"     => "$parallel_directory/dumps",
    };
    # misc -dmp files
    foreach $file ( @files_rundir ){
        if( $file =~  /^$cmd{pname}-(\S+)_dmp\d+/ ){
            $dir = $1;
            $cleanfile{"-${dir}_dmp"} = {
                "file"     => ".*/$cmd{pname}-${dir}_dmp[0-9]+.*",
                "tag"      => "no",
                "dir"      => "$parallel_directory",
                "dest"     => "$parallel_directory/dumps_${dir}",
            }
        }
    }
    # -tx files
    foreach $file ( @files_rundir ){
        if( $file =~  /^$cmd{pname}-(tx\d+\S+)\d{6}/ ){
            $dir = $1;
            $cleanfile{"$dir"} = {
                "file"     => ".*/$cmd{pname}-${dir}[0-9]*",
                "tag"      => "no",
                "dest"     => "txfiles/$dir",
            }
        }
    }
    # hdfs
    # Forced hdfs ("$cmd{pname}-H" vs "$cmd{pname}-h" kept separate
    foreach $file ( @files_rundir ){
        if( $file =~  /^$cmd{pname}-((h\d+)(\S+))(\d{6})$/ ){
            $dir = $1;
            $dest = "${2}_${3}";
            $cleanfile{"$dir"} = {
                "file"     => ".*/$cmd{pname}-${dir}[0-9]*",
                "tag"      => "no",
                "dest"     => "hdfs/${dest}",
            }
        }
    }
}
# extras for setting cleanfile
# NOTE: I assume that this is a work in progress. I don't see any references to the stubs in run_job_cleanup_extras.pm.
#       I assume that you want to do something special for eap runs. Lap may want something similar.

# BEGIN {
#     my $file = $condition ? $Module1 : $Module2;
#     eval "require $Module";
# }

#require run_job_cleanup_extras.pm;
# if( ! defined($include_found) ){
#     $include = "$lib_path/../../Tools/General/run_job_cleanup_extras.pm";
#     if( -e "$include" ){
#         $include_found = "";
#         require $include;
#     }
# }
# if( ! defined($include_found) ){
#     $include = "$lib_path/run_job_cleanup_extras.pm";
#     if( -e "$include" ){
#         $include_found = "";
#         require $include;
#     }
# }
# if( ! defined($include_found) ){
#     $include = "$lib_path/../../Tools.rh/General/run_job_cleanup_extras.pm";
#     if( -e "$include" ){
#         $include_found = "";
#         require $include;
#     }
# }

# check for consistency and defaults for cleanfile
foreach $class ( sort keys %cleanfile ){
    foreach $field ( sort keys %{$cleanfile{"$class"}} ){
        if( $field !~ /^(copy|dest|dir|file|keeplast|overwrite|tag)$/ ){
            $ierr = 1;
            &print_error( "Invalid field [$field] for class [$class]",
                          $ierr );
            exit( $ierr );
        }
    }
    if( ! defined($cleanfile{"$class"}{file}) ){
        $ierr = 1;
        &print_error( "Must define field [file] for class [$class]",
                      $ierr );
        exit( $ierr );
    }
    if( ! defined($cleanfile{"$class"}{dir}) ){
        $cleanfile{"$class"}{dir} = ".";
    }
    if( ! defined($cleanfile{"$class"}{tag}) ){
        $cleanfile{"$class"}{tag} = "yes";
    }
    if( ! defined($cleanfile{"$class"}{dest}) ){
        $cleanfile{"$class"}{dest} = ".";
    }
    if( ! defined($cleanfile{"$class"}{copy}) ){
        $cleanfile{"$class"}{copy} = "no";
    }
    if( ! defined($cleanfile{"$class"}{overwrite}) ){
        $cleanfile{"$class"}{overwrite} = "yes";
    }
}
# tag/move files
if( keys %cleanfile ){
    printf( " %20s %25s (%5s files) -> %s %s\n", "class", "regexp", "num", "dest", "skip" );
    printf( " %20s %25s -%5s ------ -> %s %s\n", "-----", "------", "---", "----", "----" );
}
foreach $class ( sort keys %cleanfile ){
    $dir       = $cleanfile{"$class"}{dir};
    $file_r    = $cleanfile{"$class"}{file};
    $dest      = $cleanfile{"$class"}{dest};
    $tag       = $cleanfile{"$class"}{tag};
    $copy      = $cleanfile{"$class"}{copy};
    $keeplast  = $cleanfile{"$class"}{keeplast};
    $overwrite = $cleanfile{"$class"}{overwrite};
    if( -d $dir ){
        # ancient find on yellowrail cannot handle regextype
        if( $cmd{sys_info}{L_MACHINE} eq "YELLOWRAIL" ){
            $regextype = "";
        }
        else{
            $regextype = "-regextype egrep";
        }
        $command = "find $dir -maxdepth 1 $regextype -regex '$file_r'";
        $output = `$command`;
        chomp( $output );
        @files = sort split( /\n/, $output );
    }
    else{
        undef( @files );
    }
    if( $#files >= 0 ){
        if( ! defined($cmd{debug}) ){
            &my_mkdir( $dest, $cmd{group}, $cmd{mode_dir_o} );
        }
        if( defined($cmd{group}) ){
            &fix_dir_perms( \%cmd, $dest );
        }
    }
    # set permissions
    if( defined($cmd{group}) ){
        foreach $file ( @files ){
            $command = "chgrp $cmd{group} $file";
            &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
            if( -d $file ){
                $mode = $cmd{mode_dir};
            }
            elsif( -x $file ){
                $mode = $cmd{mode_filex};
            }
            else{
                $mode = $cmd{mode_file};
            }
            $command = "chmod $mode $file";
            &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug}  );
        }
        $command = "chgrp $cmd{group} $dir $dest";
        &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
        $command = "chmod $cmd{mode_dir} $dir $dest";
        &run_command( COMMAND=>"$command", DEBUG=>$cmd{debug} );
    }
    # if keeplast is defined
    if( defined( $keeplast ) ){
        # if just blank, just keep the last file alphabetically
        if( $keeplast !~ /\S/ ){
            $path_keeplast = pop( @files );
        }
        # if an actual file
        else{
            # get fullpath of keeplast
            $ierr = &my_stat( $keeplast, \%stat );
            # strange timeout of stat - so do not move files
            if( $ierr == 2 ){
                print "Timeout hit doing stat on file [$keeplast]\n";
                print "NOT moving this class [$class]\n";
                undef( @files );
            }
            # if not a timeout, then move files
            else{
                $path_keeplast = $stat{fullpath};
                @files_new = @files;
                undef( @files );
                # compare fullpath of keeplast and fullpath of file and add
                # it to list to move.
                # need to compare full paths since might be looking at:
                #   "a/b/c" and "/foo/bar/a/b/c" or "../../a/b/c"
                foreach $file ( @files_new ){
                    &my_stat( $file, \%stat );
                    if( %stat ){
                        $path = $stat{fullpath};
                        if( ! defined($path_keeplast) || $path ne $path_keeplast ){
                            push( @files, $file );
                        }
                    }
                }
            }
        }
    }
    $num = $#files + 1;
    if( $class =~ /^($cmd{skip})$/ ){
        $skip = "SKIP";
    }
    else{
        $skip = "";
    }
    printf( " %20s %25s (%5s files) -> %s %s\n", $class, "($file_r)", $num, $dest, $skip );
    if( $skip =~ /\S/ ){
        next;
    }
    foreach $file ( @files ){
        $basename = basename( $file );
        $dest_full = "$dest/$basename";

        # if tag=check, move old file to tagged version before overwriting it
        # if header has changed
        if( $tag eq "check" &&
            ($overwrite eq "yes" || ! -e $dest_full) ){
            if( -e $file && -e $dest_full ){
                $head_new = `head -2 $file`;
                $head_old = `head -2 $dest_full`;
                if( $head_new ne $head_old ){
                    $cmd = "mv -f $dest_full $dest_full.$cmd{tag}";
                    &run_command( COMMAND=>$cmd, ERROR_REGEXP=>'/\S/', DEBUG=>$cmd{debug} );
                }
            }
        }

        # now copy file to tag or straight name
        if( $copy eq "yes" ){
            $op = "cp -f -p";
        }
        else{
            $op = "mv -f";
        }
        if( $tag eq "yes" ){
            $dest_full = "$dest/$basename.$cmd{tag}";
        }
        $cmd = "$op $file $dest_full";
        if( $overwrite eq "yes" || ! -e $dest_full ){
            &run_command( COMMAND=>$cmd, ERROR_REGEXP=>'/\S/', DEBUG=>$cmd{debug} );
        }

    }
}

# check
if( defined($cmd{check}) ){
    undef( %found );
    # get output file lines
    if( open( FILE, $cmd{checkfile} ) ){
        @lines = <FILE>;
        close( FILE );
    }
    else{
        @lines = ();
    }
    $lines_join = join("",@lines);
#     if( defined($HAS_EXTRAS) ){
#       &run_job_cleanup_extras_check();
#     }
    undef( $not_fatal );

    # eap
    if( $cmd{proj} eq "eap" ){
        # if not really a fatal error
        # look for hangs first - hangs are treated as retry
        if( -e "$cmd{pname}-HANG_DETECTED" ){
            if( ! defined($cmd{debug}) ){
                unlink("$cmd{pname}-HANG_DETECTED");
            }
            push( @{$found{RJ_RETRY}}, "$cmd{pname}-HANG_DETECTED" );
            $not_fatal = "";
        }

        # see if abort was caused by:
        #   unable to open file that is really there
        #   unable to create directory that already exists
        @lines_grep = grep( /(open failed for file|creation of directory)/i, @lines );
        if( $#lines_grep >= 0 ){
            # default to ok
            $ok = "yes";
            foreach $line_grep ( @lines_grep ){
                # file really is there
                if( $line_grep =~ /open failed for file (\S+) Reason: Permission denied/i ){
                    $file = $1;
                    # if file is not really there and readable, not ok
                    if( ! -r $file ){
                        $ok = "no";
                    }
                    # OS thinks file is there and readable, so retry
                    else{
                        push( @{$found{RJ_RETRY}}, "Transient Problem: permission denied [$file]" );
                    }
                }

                # directory really is there
                if( $line_grep =~ /cio_mkdir: creation of directory (\S+) failed/i ){
                    $file = $1;
                    # if file is not really there and readable, not ok
                    if( ! -d $file || ! -w $file || ! -r $file ){
                        $ok = "no";
                    }
                    # OS thinks file is there and readable, so retry
                    else{
                        push( @{$found{RJ_RETRY}}, "Transient Problem: mkdir failed [$file]" );
                    }
                }
            }
            if( $ok eq "yes" ){
                $not_fatal = "";
            }
        }

        # not a great thing to do, but if captured a SEGV, have it not be fatal
        # rpw run saw this...might need to fine-tune it depending on where the
        # segv hits
        # rpw: $match_string = 'FATAL_ERROR.*SEGFAULT.*whichcells';
        $match_string = 'FATAL_ERROR.*SEGFAULT';
        @lines_grep = grep( /$match_string/i, @lines );
        if( $#lines_grep >= 0 ){
            $not_fatal = "";
            push( @{$found{RJ_RETRY}}, "Recoverable(?) SEGFAULT" );
        }

        # RJ_STOP/RJ_RETRY: if BADBOY, but might not be RJ_STOP
        @lines_grep = grep( /BADBOY/, @lines );
        if( $#lines_grep >= 0 && ! defined($not_fatal) ){
            $match_string_retry =
                '('.
                '\[get_shm_id\] Cannot get initial id for KEY_TEOS'.
                ')';
            $match_string_new_child =
                '('.
                'pwrite failed for pos|'.
                'FNF BULKIO_CLOSE: Too many open/close|'.
                'ccwrite: write error'.
                ')';
            if( $lines_join =~ /$match_string_retry/ ){
                $not_fatal = "";
                push( @{$found{RJ_RETRY}}, "BADBOY: $1" );
            }
            elsif( $lines_join =~ /$match_string_new_child/ ){
                push( @{$found{RJ_NEW_CHILD}}, "BADBOY: $1" );
            }
            else{
                push( @{$found{RJ_STOP}}, "BADBOY" );
            }
        }

        # RJ_STOP: strings from code
        $match_string =
            '('.
            'GLOBAL_ERROR|'.
            'FATAL ERROR'.
            ')';
        if( ! defined($not_fatal) && $lines_join =~ /$match_string/ ){
            push( @{$found{RJ_STOP}}, "$1" );
        }

        # RJ_STOP: strings from mpi
        $match_string =
            '('.
            'application called MPI_Abort|'.
            'MPI_ABORT was invoked on rank'.
            ')';
        @lines_grep = grep( /$match_string/, @lines );
        if( ! defined($not_fatal) && $lines_join =~ /$match_string/ ){
            $match = $1;
            $match_string =
                '('.
                'ccwrite: write error'.
                ')';
            if( $lines_join =~ /${match_string}/ ){
                push( @{$found{RJ_NEW_CHILD}}, "$match: $1" );
            }
            else{
                push( @{$found{RJ_STOP}}, "$match" );
            }
        }

        # RJ_STOP/RJ_RETRY: kill/restart
        @lines_grep = grep( /-(kill|restart) file found\s*$/, @lines );
        if( $#lines_grep >= 0 ){
            @lines_grep_r = grep( /-(restart) file found\s*$/, @lines_grep );
            $restart_name = "$cmd{pname}-restart";
            $kill_name = "$cmd{pname}-kill";
            if( -e $restart_name || $#lines_grep_r >= 0 ){
                push( @{$found{RJ_RETRY}}, "$restart_name file found" );
            }
            else{
                push( @{$found{RJ_STOP}}, "$kill_name file found" );
            }
            if( -e $restart_name ){
                if( ! defined($cmd{debug}) ){
                    unlink( $restart_name );
                }
            }
            if( -e $kill_name ){
                if( ! defined($cmd{debug}) ){
                    unlink( $kill_name );
                }
            }
        }

        # RJ_STOP: -DO_NOT_RUN
        if( -e "$cmd{pname}-DO_NOT_RUN" ){
            push( @{$found{RJ_STOP}}, "$cmd{pname}-DO_NOT_RUN" );
        }
        # users want to do some more processing after tmax is reached
        # do not stop code because of this.
        # need to set "uselast = .true." and "kread = -1" to generate DO_NOT_RUN file
        #@lines_grep = grep( /CONTROLLER: tmax\.ge\.0 \.and\. time\.ge\.tmax/, @lines );
        #if( $#lines_grep >= 0 ){
        #    print STDERR "RJ_STOP - tmax reached\n";
        #    $ierr = 4;
        #}
    }

    # do appropriate action
    # print all findings
    # last one will be final error return
    foreach $type ( "RJ_RETRY", "RJ_NEW_CHILD", "RJ_STOP" ){
        foreach $msg ( @{$found{$type}} ){
            print "$type: [$msg]\n";
            $ierr = $RET{$type};
        }
    }
}
&my_exit( $ierr );

###################################################################################

# fix_dir_perms
# fix current and parent directories to be group readable
sub fix_dir_perms{
    my( $cmd_ref, $path ) = @_;
    my(
        $command,
        $done,
        $mode,
        %stat,
        );
    $done = "false";
    while( $done eq "false" ){
        &my_stat($path, \%stat);
        # leave if not owner of directory
        if( $stat{uid} != $< ){
            $done = "true";
            last;
        }
        # set mode to be at least as open as currently set plus what
        # is required for group read
        $mode = sprintf( "%o", $stat{mode}|$$cmd_ref{mode_dir_o} );
        $mode =~ s/^.*(\d{4})/$1/; # last 4 digits (since stat has type prepended)
        $command = "chgrp $$cmd_ref{group} $path";
        &run_command( COMMAND=>"$command", DEBUG=>$$cmd_ref{debug} );
        $command = "chmod $mode $path";
        &run_command( COMMAND=>"$command", DEBUG=>$$cmd_ref{debug} );
        # done if at /
        if( $path eq "/" ){
            $done = "true";
            last;
        }
        $path = $stat{dir};
    }
}
sub my_exit{
    my( $ierr ) = @_;
    print "\nFinished: $0\n";
    print `date`,"\n";
    exit( $ierr );
}

###################################################################################

#............................................................................
#...Name
#...====
#... eap_hcard
#...
#...Purpose
#...=======
#... do hcard parsing of eap input file
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: in
#...              Perl type: reference to hash
#...              command line
#...              $cmd{$option} = value
#...              $cmd{files}[] = array of file names
#...
#...Program Flow
#...============
#... 1) get time from -status file
#... 2) go through input file and do hcard processing
#... 3) if changes, write new input file.
#............................................................................
sub eap_hcard{
    my(
       $cmd_ref,
       ) = @_;
    my(
       $changed,
       %data,
       $dim,
       @files,
       $first_cycle,
       %found_fields,
       $hrest,
       $hval,
       $ierr,
       $input_file,
       $line,
       $line_new,
       @lines,
       @lines_new,
       $nummat,
       $status_file,
       $time,
       $white,
       );
    # try to get time from -status file (might be in logs directory)
    $input_file = $$cmd_ref{in};
    if( ! defined($input_file) || ! -e $input_file ){
        return;
    }
    $status_file = "$$cmd_ref{pname}-status";
    # if -status does not exist, get last one in logs directory"
    if( ! -e "$status_file" ){
        @files = glob( "logs/$status_file.*" );
        if( $#files >= 0 ){
            $status_file = $files[-1];
        }
    }
    if( ! open( FILE, "$status_file") ){
        return;
    }
    @lines = <FILE>;
    close( FILE );
    $first_cycle = "true";
    undef( $nummat );
    undef( $dim );
    &parse_output_file( \@lines,
                        \$first_cycle, \$nummat,
                        \%data, \%found_fields, \$dim );
    &parse_output_file_finish( \%data, \%found_fields, $cmd_ref );
    $time = $data{time}[-1];
    # if did not seem to find time, just return
    if( ! defined($time) || $time !~ /\d/ ){
        return;
    }
    if( ! open( FILE, "$input_file") ){
        $ierr = 0;
        &print_error( "Cannot open input file [$input_file] for hcard replacement [reading].",
                      $ierr );
        return;
    }
    @lines = <FILE>;
    close( FILE );
    $changed = 0;
    foreach $line ( @lines ){
        if ( ($line =~ /^(\s*)(\!h\s+(\d+\.{0,1}\d*[eE]{0,1}[\-\+]{0,1}\d*))\s(\s*\S+.*)$/) and ( $3 <= $time) ) {
            my $white   = $1;
            my $hval    = $2;
            my $hrest   = $4;
            $line_new = "$white$hrest $hval\n";
            $changed++;
        }
        else{
            $line_new = "$line";
        }
        push( @lines_new, $line_new );
    }
    # hcard replacement
    if( $changed > 0 ){
        # if doing the "check", remove the DO_NOT_RUN file so that the code will continue
        print "\n\n ** hcard replacement detected ** input file [$input_file] time [$time] lines [$changed]\n";
        if( defined($$cmd_ref{check}) ){
            print " ** hcard replacement detected ** removing -DO_NOT_RUN file [$$cmd_ref{pname}-DO_NOT_RUN]\n";
            unlink( "$$cmd_ref{pname}-DO_NOT_RUN" );
        }
        # only do replacement if starting the run
        else{
            print " ** hcard replacement detected ** creating new input file [$input_file]\n";
            if( ! defined($$cmd_ref{debug}) ){
                if( ! open( FILE, ">$input_file") ){
                    $ierr = 0;
                    &print_error( "Cannot open input file [$input_file] for hcard replacement [writing].",
                                  $ierr );
                    return;
                }
                print FILE @lines_new;
                close ( FILE );
            }
        }
        print "\n";
    }
}

###################################################################################

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
       @args, # arguments
       $count, # counter
       $ierr, # error ret val
       $in,
       $infile_guess,
       $pname_guess,
       $mode,
       $myumask,
       $num_args, # number of arguments
       $opt, # current option
       $output,
       $string,
       $val, # value for current option
      );
    $ierr = 0;
    @args = @{$argv_ref};
    $$cmd_ref{conds} = "";
    #....................
    #...parse the args...
    #....................
    $num_args = $#args;
    while( @args ) {
        $opt = shift( @args );
        #..........
        #...help...
        #..........
        if( $opt =~ /^-+(h(elp)?)$/i ) {
            $opt = "h";
            $$cmd_ref{$opt} = "true";
        }
        #...debug
        elsif( $opt =~ /^-+(d(ebug)?)$/i ) {
            $opt = "debug";
            $$cmd_ref{$opt} = "true";
        }
        #....................
        #...--<opt> <val>...
        #....................
        elsif( $opt =~ /^--(checkfile|group|in|pname|proj|skip|tag|tagdir|tagfile|tagtype)$/ ) {
            $opt = $1;
            if( ! @args ) {
                $ierr = 1;
                &print_error( "Value needed for option [--$opt].",
                              $ierr );
                exit( $ierr );
            }
            $val = shift( @args );
            if( $opt eq "tagtype" && $val ne "new" && $val ne "old" ){
                $ierr = 1;
                &print_error( "Value for option [$opt] must be either 'new' or 'old'",
                              $ierr );
                exit( $ierr );
            }
            if( $opt eq "skip" ){
                $$cmd_ref{$opt} .= "|$val";
            }
            else{
                $$cmd_ref{$opt} = "$val";
            }
        }
        #.............
        #...--<opt>...
        #.............
        elsif( $opt =~ /^--(check|clean)$/i ) {
            $opt = "$1";
            $$cmd_ref{$opt} = "true";
        }
        #.............
        #...unknown...
        #.............
        else{
            $ierr = 1;
            &print_error( "Invalid option [$opt]",
                          $ierr );
            exit( $ierr );
        }
    }
    # default proj
    if( ! defined($$cmd_ref{proj}) ){
        $$cmd_ref{proj} = "eap";
    }
    # help - return here after setting some defaults if help
    if( defined($$cmd_ref{h}) ){
        return( $ierr );
    }
    # tag
    if( ! defined($$cmd_ref{tagdir}) ){
        if( defined($ENV{RJ_VAL_DIR_ROOT}) ){
            $$cmd_ref{tagdir} = $ENV{RJ_VAL_DIR_ROOT};
        }
        else{
            $$cmd_ref{tagdir} = ".";
        }
    }
    # cannot define more than 1 specification for tag to use
    $count = 0;
    if( defined($$cmd_ref{tag}) ){
        $count++;
    }
    if( defined($$cmd_ref{tagfile}) ){
        $count++;
    }
    if( defined($$cmd_ref{tagtype}) ){
        $count++;
    }
    if( $count > 1 ){
        $ierr = 1;
        &print_error( "Cannot specify more than 1 of tag, tagfile, tagtype",
                      $ierr );
        exit( $ierr );
    }
    if( ! defined($$cmd_ref{tag}) && ! defined($$cmd_ref{clean})){
        # get tagfile to use
        if( ! defined($$cmd_ref{tagfile}) ){
            if( defined($$cmd_ref{tagtype}) ){
                if( $$cmd_ref{tagtype} eq "new" ){
                    $$cmd_ref{tagfile} = "$$cmd_ref{tagdir}/$RJ_FILE_TAG";
                }
                else{
                    $$cmd_ref{tagfile} = "$$cmd_ref{tagdir}/$RJ_FILE_TAG_OLD";
                }
            }
            # if --check, then use the new tag regardless of where called from
            elsif( defined($$cmd_ref{check}) ){
                $$cmd_ref{tagfile} = "$$cmd_ref{tagdir}/$RJ_FILE_TAG";
                $$cmd_ref{tagtype} = "new";
            }
            # if called from within run_job.pl, assume called right before new run is done
            elsif( defined($ENV{RJ_RUN_JOB}) ){
                $$cmd_ref{tagfile} = "$$cmd_ref{tagdir}/$RJ_FILE_TAG_OLD";
                $$cmd_ref{tagtype} = "old";
            }
            # if called standalone
            else{
                $$cmd_ref{tagfile} = "$$cmd_ref{tagdir}/$RJ_FILE_TAG";
                $$cmd_ref{tagtype} = "new";
            }
        }
        else{
            $$cmd_ref{tagtype} = "unknown";
        }
        # convert tag file to contents of tag
        if( open( FILE, "$$cmd_ref{tagfile}" ) ){
            $$cmd_ref{tag} = <FILE>;
            chomp( $$cmd_ref{tag} );
            close FILE;
        }
        else{
            $ierr = 0;
            &print_error( "Cannot open tagfile [$$cmd_ref{tagfile}]",
                          "This can occur when running $0 for the very first time.",
                          "Using current date for tag.",
                          $ierr );
        }
    }
    if( ! defined($$cmd_ref{tag}) ){
        $$cmd_ref{tag} = &date_ymdhms();
    }
    if( ! defined($$cmd_ref{tagtype}) ){
        $$cmd_ref{tagtype} = "unknown";
    }
    if( ! defined($$cmd_ref{tagfile}) ){
        $$cmd_ref{tagfile} = "unknown";
    }
    # pname
    $infile_guess = "";
    $pname_guess = &get_pname(".",INFILE=>\$infile_guess);
    if( ! defined($$cmd_ref{pname}) ){
        $$cmd_ref{pname} = $pname_guess;
    }
    if(  $$cmd_ref{pname} !~ /\S/ ){
        $ierr = 1;
        &print_error( "Must define --pname",
                      $ierr );
        exit( $ierr );
    }
    # in
    if( ! defined($$cmd_ref{in}) ){
        $string = "Name of input deck file to read in is:";
        if( -e $RJ_FILE_CMD_OUT ){
            $output = `grep "$string" $RJ_FILE_CMD_OUT`;
            if( $output =~ /${string}\s+(\S+)\s*$/ ){
                $$cmd_ref{in} = $1;
            }
        }
    }
    if( ! defined($$cmd_ref{in}) ){
        $in = "$$cmd_ref{pname}.in";
        if( -e $in ){
            $$cmd_ref{in} = $in;
        }
    }
    if( ! defined($$cmd_ref{in}) ){
        $in = "$$cmd_ref{pname}.input";
        if( -e $in ){
            $$cmd_ref{in} = $in;
        }
    }
    if( ! defined($$cmd_ref{in}) ){
        if( $infile_guess =~ /\S/ ){
            $$cmd_ref{in} = $infile_guess;
        }
    }
    if( ! defined($$cmd_ref{in}) ){
        $$cmd_ref{in} = "";
    }
    # group
    if( ! defined($$cmd_ref{group}) && defined($ENV{RJ_VAR_GROUP}) ){
        $$cmd_ref{group} = $ENV{RJ_VAR_GROUP};
    }
    # set default modes for things
    # if group is set, allow group rx
    if( defined($$cmd_ref{group}) ){
        $myumask = umask() & 0727;
        umask($myumask);
    }
    $myumask = umask();
    $mode = 02777 & ~$myumask;
    $$cmd_ref{mode_dir}   = sprintf( "0%o", $mode );
    $$cmd_ref{mode_dir_o} = $mode;
    $mode = 0666 & ~$myumask;
    $$cmd_ref{mode_file}   = sprintf( "0%o", $mode );
    $$cmd_ref{mode_file_o} = $mode;
    $mode = 0777 & ~$myumask;
    $$cmd_ref{mode_filex}   = sprintf( "0%o", $mode );
    $$cmd_ref{mode_filex_o} = $mode;
    # if not skipping anything
    if( !defined($$cmd_ref{"skip"}) ){
        $$cmd_ref{skip} = "";
    }
    $$cmd_ref{skip} =~ s/^\|//;
    # move is default
    if( ! defined( $$cmd_ref{check} ) ){
        $$cmd_ref{move} = "yes";
    }
    if( ! defined( $$cmd_ref{checkfile} ) ){
        $$cmd_ref{checkfile} = "";
        if( defined($ENV{RJ_VAL_DIR_ROOT}) ){
            $$cmd_ref{checkfile} = "$ENV{RJ_VAL_DIR_ROOT}/";
        }
        $$cmd_ref{checkfile} .= "$RJ_FILE_CMD_OUT";
    }
    return( $ierr );
}
