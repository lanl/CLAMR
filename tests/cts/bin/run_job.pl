eval 'exec perl -w -S $0 ${1+"$@"}'
if 0;

# NOTE: You can use FindBin in the script. The modules will automatically have access to $FindBin::Bin.
use FindBin qw($Bin);
use lib ("$Bin", "$Bin/lib", "$Bin/../lib");
local $ENV{PATH} = "$Bin:$ENV{PATH}";

use my_utils qw (
                 conv_time
                 date_ymdhms
                 date_ymdhms_sep
                 get_id_num
                 get_sysinfo
                 get_pname
                 my_checksum_string
                 my_compare_version
                 my_getval
                 my_lockfile
                 my_mkdir
                 my_mode
                 my_notdir
                 my_stat
                 path_add_default
                 print_error
                 print_perl_obj
                 run_command
                 sort_unique
                 status_bar
                 which_exec
                 );
use POSIX;
use Cwd 'abs_path';
use File::Basename 'dirname';

my(
   $ierr,
   );

# catch certain signals
$SIG{URG}  = 'handler_resume';
$SIG{HUP}  = 'handler_resume';
$| = 1;
$ierr = 0;
# Add to path in case this is run from the checkout
&path_add_default();
# start time
$output = `date +%s`;
chomp( $output );
$cmd{time_start} = $output;
# names (or suffixes) used for various files
$RJ_SHELL_SUFFIXES{csh}  = ".csh";
$RJ_SHELL_SUFFIXES{tcsh} = ".csh";
$RJ_SHELL_SUFFIXES{bash} = ".sh";
$RJ_SHELL_SUFFIXES{sh}   = ".sh";
# initial RJ_SHELL
# If not /bin/tcsh, then pick bash
if( ! defined($RJ_SHELL) ){
    if( ! -e "/bin/tcsh" && -e "/bin/bash" ){
        $RJ_SHELL = "bash";
    }
}
# if SHELL defined, use it
if( ! defined($RJ_SHELL) ){
  if( defined($ENV{SHELL}) ){
      if( $ENV{SHELL} =~ m&(csh|tcsh)$& ){
          $RJ_SHELL                = "tcsh";
      }
      else{
          $RJ_SHELL                = "bash";
      }
  }
}
# default tcsh
if( ! defined($RJ_SHELL) ){
    $RJ_SHELL                = "tcsh";
}
$RJ_SHELLS               = join( "|", sort keys %RJ_SHELL_SUFFIXES );

$RJ_DIR                 = "rj_adir";
$RJ_FILE_TAG            = "$RJ_DIR/rj_tag";
$RJ_FILE_TAG_OLD        = "$RJ_DIR/rj_tag_old";
$RJ_FILE_TAGS           = "$RJ_DIR/rj_tags";
# rj_mapping.txt and run_job_mapping appear here and in run_job.pl - keep sync'd
$RJ_EXEC_MAPPING        = "run_job_mapping";
$RJ_FILE_MAPPING        = "rj_mapping.txt";
# Mapping file might have to go in global space
# so make its name have a checksum of the current working
# directory to avoid file name space collision
# right now, if running multiple jobs in the background
# in the same directory, you will get a collision - but
# you would get that collision with old behavior
# of the name just being $RJ_FILE_MAPPING as well.
$dir = "/nfs/tmp2";
if( ! -d $dir ){
    $dir = ".";
}
$suffix = ".".&my_checksum_string(&abs_path("."), 7);
$RJ_FILE_MAPPING_FULL   = "$dir/${RJ_FILE_MAPPING}.$ENV{LOGNAME}${suffix}";

$RJ_FILE_CMD_OUT        = "rj_cmd_out";
$RJ_FILE_BATCH_OUT_BASE = "rj_batch_out";
$RJ_FILE_MULTI          = "rj_multi_";
$RJ_FILE_BATCH_OUT      = "$RJ_DIR/$RJ_FILE_BATCH_OUT_BASE";
$RJ_FILE_ID             = "$RJ_DIR/rj_id";
$RJ_FILE_ID_NEXT        = "$RJ_DIR/rj_id_next";
$RJ_FILE_LOG            = "$ENV{HOME}/rj_log";
$RJ_FILE_LOG_SL         = 10000; # number of lines to keep in log file
$RJ_VAL_DIR = dirname(abs_path($0));
$RJ_VAL_DIR_ROOT = getcwd();

# some more defaults
$RJ_ENV_EMPTY        = "rj_env_empty";
$DEFAULT_TIME        = "8h";
$DEFAULT_TRY_TIME    = "5m";
$DEFAULT_TRY_NUM     = 10;
$DEFAULT_TRY_NUM_MAX = 20;
$DEFAULT_TRY_NUM_MAX_NUM = 2;
$DEFAULT_TIME_B      = "900s";
$DEFAULT_NUMPE       = 1;

#..................................
#...DONE: setup for system files...
#..................................
# execs
$RUN_STATUS_PL = "$Bin/run_status.pl";
# mac does not have ldd - could use "otool -L <exec>"
$LDD = &which_exec("ldd", QUIET=>"");

# parse arguments
&parse_args( \@ARGV, \%cmd );
# parse RUN_JOB_CMD as command line arguments
# if RUN_JOB_CMD is defined and it is not equal to the previously parsed one, 
# tack it on.
# Perhaps some better way to do this.
# Advantages:
#   o do not need to pass an environment var to child run_job.pl
#   o use the parsing in parse_args for RUN_JOB_CMD
# Disadvantages:
#   o If user wishes to set RUN_JOB_CMD in their .cshrc and then change it on
#     the fly, might have issues.
#   o Might have trouble with imbedded quotes, back-ticks
if( defined($ENV{RUN_JOB_CMD}) &&
    ( ! defined($ENV{RJ_RUN_JOB_CMD_OLD}) ||
      $ENV{RUN_JOB_CMD} ne $ENV{RJ_RUN_JOB_CMD_OLD} ) ){
    @args_env = split( /\n/, `perl -e 'foreach \$arg ( \@ARGV ) {print \$arg,"\n"} ' -- $ENV{RUN_JOB_CMD}` );
    &parse_args( \@args_env, \%cmd );
    $ENV{RJ_RUN_JOB_CMD_OLD} = $ENV{RUN_JOB_CMD};
}

# get system info and stuff into cmd
&get_sysinfo( \%{$cmd{sys_info}} );
$DEFAULT_PPN         = $cmd{sys_info}{L_PPN};
$cmd{conds_total} = $cmd{conds};
# add in sys_info values that are not directories (no "/")
@conds = grep( !/\//, values(%{$cmd{sys_info}}) );
# remove empty ones
@conds = grep(/\S/, @conds);
if( @conds ){
    $cmd{conds_total} .= ",".join(",",@conds);
}
$cmd{conds_total} =~ s/^,//;
foreach $key ( keys %{$cmd{sys_info}} ){
    $sysinfo_stuff .= sprintf( "#...       %-25s = %s\n", "RJ_$key", $cmd{sys_info}{$key} );
}
$DEFAULT_TPP = 1;

# get shell dependent names of things (RJ_FILE_DATAFILE, RJ_SHELL, ... )
&rj_getshell( \%cmd, "." );

# help message
if( defined( $cmd{h} ) ) {
    print <<"EOF";
#............................................................................
#...Synopsis
#...========
#... run_job.pl go
#...   Runs a job given the data inside $RJ_FILE_DATAFILE .
#...
#... To get a basic $RJ_FILE_DATAFILE, run:
#...   $0 --convert basic
#...   (see the --convert option below for other conversions)
#...
#... A log file of the latest runs is put into $RJ_FILE_LOG.  See "Files Created"
#... below.
#...
#...Options
#...=======
#... run_job.pl go [options]
#...   [go] required flag if no options given - otherwise, you get this help message
#...   [<EXEC>] [<EXEC_ARGS>]
#...     If you specify EXEC, the first unrecognized argument will set the EXEC .
#...     Then ALL following args will be considered EXEC_ARGS .  So, all args to
#...     run_job.pl need to go before specifying the EXEC.
#...     If no $RJ_FILE_DATAFILE is found, a dummy datafile is used with the line:
#...        &&&RJ_CMD_PRUN&&&
#...     This allows doing something like the following without a $RJ_FILE_DATAFILE:
#...        run_job.pl --numpe 128 ./foo.x foo.in
#...     Which submits a job running something like the following in a batch submission:
#...        mpirun -np 128 ./foo.x foo.in
#...     So, you get the benefits of run_job.pl:
#...       saving screen output, using correct "mpirun", ...
#...     See the "-i" interactive flag below for running jobs interactively.
#...   [--conds <comma separated list of conditions to be used in $RJ_FILE_DATAFILE>]
#...     --conds a,b,c means that a, b, and c conditions are set.
#...     Conditions automatically set for this machine:
#...       $cmd{conds_total}
#...   [--convert <file>]
#...     Create a $RJ_FILE_DATAFILE file from a run file from another system
#...     File types allowed:
#...       job.params:   old EAP run scripts
#...       <test>.test:  CTS test file
#...       basic:        create a basic $RJ_FILE_DATAFILE from answering questions
#...   [--debug]
#...       Create various files, but do not actually submit|run the job.
#...       Useful for checking the contents of $RJ_FILE_DATAFILE -> $RJ_FILE_CMD_BASE
#...   [--depend <comma separated list of job ids>]
#...     By default, a job id is assumed to be an id from the batch system.
#...     To specify a process id, add the suffix ":PID"  :
#...        --depend 12345:PID
#...     Job will not run until the list of jobs (and/or pid's) are finished.
#...     By default, new jobs will not be allowed to be submitted if there are any
#...     existing jobs running in the current directory.
#...     You can run "run_status.pl" to get the jobs ids of currently running jobs.        
#...     Special jobids:
#...      --depend .  : Depend upon      any jobs currently running in the directory.
#...      --depend no : Do not check for any jobs currently running in the directory.
#...                    This can lead to errors since files from the different
#...                    jobs can overwrite eachother (rj_cmd_out, rj_cmd.csh, ... ).
#...                    This option should NOT be used.
#...   [--env <var>[=val]]
#...     Specify an environment variable to be set before running $RJ_FILE_DATAFILE.
#...   [--help] [-h]
#...     This help message.
#...   [-i]
#...     This will run your $RJ_FILE_DATAFILE file interactively.
#...     Useful for running when you are in an interactive allocation.
#...     Implications (overriding $RJ_FILE_DATAFILE / command line):
#...        on -> BATCH=no, CHAIN=no
#...   [--info]
#...     Print out various bits of information about the job - do not run it.
#...        <key> = <val>
#...   [--mail <comma delimited email list>]
#...     For sending mail of various results.
#...     Currently, only used for mailing if jobs were launched via --multi option.
#...   [--multi <datafile>]
#...     Launch a set of jobs specified in directories in the <datafile>.
#...     If the batch system says a job is already running in that directory,
#...     a new job will not be launched.  This check is NOT made for "#RJ BATCH=no"
#...     jobs or for jobs running on different batch systems (eg same directory, but
#...     launched from another machine).
#...     One directory per line in the datafile.
#...     This is useful for keeping jobs running on unstable machines/batch systems.
#...     Sample cron that runs every 6 hours (10 minutes past the hour):
#...        10 0,6,12,18 * * * csh -c 'cd ~/multi && run_job.pl --multi multi_cielo --mail $ENV{LOGNAME}\@lanl.gov'
#...   [--norun_job]
#...     Do not use the $RJ_FILE_DATAFILE file in the current directory.  Will use a template
#...     that will just run the executable.
#...   [--clean_mpi <no|yes>]
#...     At the end of each parallel run (SERIAL!=yes) of $RJ_FILE_DATAFILE, run_job.pl
#...     tries to clean up any mpi resources.  This might be needed to get the follow-on
#...     mpi execution to work.  This is not done if run_job.pl detects another instance
#...     of run_job.pl running (with perhaps its own mpi executable running).
#...     Set this to override the default behavior.
#...   [--shell <csh/tcsh or sh/bash>]
#...     run_job.pl will try to automatically detect the shell being used.  When it
#...     is arbitrary (for example, creating a $RJ_FILE_DATAFILE template with --convert,
#...     this flag can be used).
#...   [-v]
#...     Verbose
#...   [--var <VAR>=<VAL>], [--var <VAR><operation><VAL>]
#...     Define a specific RJ variable to be a value.
#...       VAR       = EXEC, BATCH_ARGS, PRUN_ARGS, NUMPE, U_MYVAR (file variable)
#...       operation = :=, =, +=
#...     Example:
#...       --var PRUN_ARGS+="-mca 1"
#...     Various quick command line options have been added for shortened forms of vars.
#...     In most cases, the usage is:
#...        --foo bar   [stands for]   --var FOO=bar
#...     Options: (see "#RJ OPT" below for more info)
#...       <val>      -> <OPT>      (examples)
#...       batch      -> BATCH      (never|no|yes)
#...       batch_args -> BATCH_ARGS (-A access)
#...       chain      -> CHAIN      (no|yes)
#...       debugger   -> DEBUGGER   (gdb|tv)
#...       exec_args  -> EXEC_ARGS  (foo.in)
#...       numpe      -> NUMPE      (4)
#...       opt        -> OPT        (no|yes)_(fix_env|mapping)
#...       pname      -> PNAME      (myrun)
#...       ppn        -> PPN        (4)
#...       prun_args  -> PRUN_ARGS  (-mca mpi_paffinity_alone 0)
#...       serial     -> SERIAL     (no|yes)
#...       time       -> TIME       (4h)
#...       tpp        -> TPP        (4)
#...   [--wait]
#...     Wait for process to finish before returning controll to the shell.
#...
#...   Internal use only:
#...     [--batchid]
#...     [--dir]
#...     [--id]
#...     [--launchtype]
#...     [--tag]
#...     [--wrapper_file], [--wrapper_out], [wrapper_tag], [--wrapper_wait],
#...     Flags used internally only.
#...     Values and format can change without notice.
#...
#...   NOTE: RUN_JOB_CMD environment variable
#...     The value of the RUN_JOB_CMD environment variable is added to the command
#...     line options.  Some examples (can be tricky if you need quotes):
#...       setenv RUN_JOB_CMD "--prun_args '--hosts a'"
#...       setenv RUN_JOB_CMD --prun_args\ \'--hosts\ \"a\ b\ c\"\'\ --numpe\ 4
#...     Useful if you have a set of existing $RJ_FILE_DATAFILE files and you wish
#...       to add some mpirun arguments without modifying all the $RJ_FILE_DATAFILE files.
#...
#...   NOTE: RJ_E_<VAR> and RJ_EP_<VAR> environment variables
#...     You can also set environment variables to define variables.
#...     Setting the following in your environment:
#...        RJ_E_NUMPE      = 2
#...        RJ_EP_PRUN_ARGS = -mca 1
#...     Translates to having the following in your datafile:
#...       #RJ NUMPE        = 2
#...       #RJ PRUN_ARGS   += -mca 1
#...
#...   NOTE: Variable Evaluation Order
#...     Last setting wins.  Order of operation (last one wins):
#...       $RJ_FILE_DATAFILE line
#...       --var Command line argument
#...       RUN_JOB_CMD environment variable
#...       RJ_E(P)?_<VAR> setting
#...     Using ":=" in the file overrides all other settings.
#...
#...$RJ_FILE_DATAFILE
#...===========
#...  This datafile is a shell script with special command tags (which are evaluated
#...  in the order listed in the file).
#...
#.... The file $RJ_FILE_CMD_BASE is created and is run by run_job.pl (output sent to
#...  $RJ_FILE_CMD_OUT).  Use the --debug option to create this file but not actually run.
#...
#...  The following are the special trings in the $RJ_FILE_DATAFILE that are
#...  interpreted by run_job.pl:
#...
#...     1) #RJ <VARIABLE>{<CONDITION>} = <VALUE>
#...        Sets a variable equal to a value if the condition is satisfied.
#...        Basic variable operations:
#...          #RJ U_MYVAR             := this is the valye of U_MYVAR (overriding command line "=")
#...          #RJ U_MYVAR              = this is the valye of U_MYVAR
#...          #RJ U_MYVAR{A|(B&C)|!D}  = new value of U_MYVAR if (A or (B and C) or not D) is true
#...          #RJ U_MYVAR             += append this to U_MYVAR (depending on last "=" or ":=")
#...
#...        Conditions automatically set for this machine:
#...          $cmd{conds_total}
#...
#...        When the parser hits a variable in the file, it evaluates it and then applies
#...        the command line variables and environment variables in order.
#...        See "NOTE: Variable Evaluation Order" above for other ways to set these variables
#...        (command line, environment variables, ... ).
#...
#...     2) &&&RJ_<NAME>&&&
#...        These get replaced by the appropriate string when creating $RJ_FILE_CMD_BASE.
#...        For example:
#...          &&&RJ_CMD_PRUN&&&   ->   mpirun -np 4 my_exec.x my_input.in
#...
#...  #RJ Variables
#...  -------------
#...    Must be on separate lines.
#...
#...    #RJ BATCH       = <yes=default|no|any|auto>
#...                      If submitting the job to a batch system.
#...                        yes   -> Submit job to batch system
#...                        no    -> Do not submit job - Run on current machine/end.
#...                        never -> Do not submit job - Error if run on back end.
#...                        auto  -> "yes" if on front end, "no" otherwise
#...    #RJ BATCH_ARGS  = <""=default>
#...                      Additional args to batch submit command.  For example,
#...                                        to run in the access qos and send a signal 12 minutes
#...                                        before the queue time end:
#...                                          -A access
#...                                          -l signal=SIGHUP\@12:00
#...    #RJ CHAIN       = <no=default|yes>
#...                      Chain parent/child dependent jobs until RJ_STOP seen.
#...                                        A child job is first launched by the parent job and
#...                                        then the parent runs the script.
#...                                        This child job waits for the parent to finish before
#...                                        starting the script (thus becoming the new parent and
#...                                        launching its own child).
#...                                        The script will be rerun in the same submission if still
#...                                        enough time left.
#...                                        If "no", no child will be launched.  The job will run
#...                                        once and then stop regardless what is seen,
#...                                         (ignoring RJ_RETRY, RJ_NEW_CHILD).
#...    #RJ DEBUGGER    = <""=default>
#...                      Run the executable (&&&RJ_CMD_PRUN&&&) under a debugger.
#...                                       Supported options:
#...                                         gdb     : gdb (1 xterm/mpi-rank)
#...                                         tv      : totalview
#...                                         <blank> : no debugger
#...                                       Implications (overriding $RJ_FILE_DATAFILE / command line):
#...                                         on --> -i (interactive)
#...    #RJ EXEC        = <""=default>
#...                      name of executable (used for RJ_CMD_PRUN)
#...    #RJ EXEC_ARGS   = <""=default>
#...                      Arguments to executable
#...    #RJ EXN         = <0=default>
#...                      Extra nodes to request as a backup nodes.
#...                                        ie: 5=5 nodes, 12%=num nodes*.12
#...                                        Suggested: "#RJ EXN = 0.05%"
#...    #RJ GROUP       = <""=default>
#...                      Make run read accessible to this unix group
#...                                        GROUP and UMASK can be used together.
#...                                        eg: GROUP=dacodes, UMASK=7 => allow group write access
#...                      Group ownership will not be modified until run_job_cleanup.pl is called.
#...    #RJ NUMPE       = <$DEFAULT_NUMPE=default>
#...                      number of processing elements (mpi ranks)
#...    #RJ PNAME       = <""=default>
#...                      problem name (will make best guess as to name)
#...                      This pname is used by run_job_cleanup.pl - so make sure it matches.
#...    #RJ OPT         = <""=default>
#...                      Miscellaneous options.
#...                      no_fix_env: do not try to fix the environment before doing mpi commands.
#...                       By default, run_job.pl will load the correct modules and set
#...                       appropriate mpi command line flags.
#...                                         no_mapping: do not use mapping file [$RJ_FILE_MAPPING_FULL]
#...                                           To help optimize io performance on sequoia, a
#...                                           mapping file is created to place io processes
#...                                           on the correct nodes.
#...    #RJ PPN         = <$DEFAULT_PPN=default>
#...                      number of processes (mpi ranks) per node
#...                                        You can use a fraction (1/2 -> half of the cores)
#...    #RJ PRUN_ARGS   = <""=default>
#...                      Arguments to parallel run command
#...    #RJ SERIAL      = <no=default|yes>
#...                      Use serial flavor for &&&RJ_CMD_PRUN&&&
#...    #RJ TIME        = <$DEFAULT_TIME=default>
#...                      Run submission time (1y 2d 3h 4m 5s)
#...    #RJ TIME_B      = <$DEFAULT_TIME_B=default>
#...                      Buffer time used to no longer try to run in the current
#...                                        run submission (eg, to give time to write the dump file)
#...    #RJ TPP         = <$DEFAULT_TPP=default>
#...                      Number of threads per process (mpi rank) (eg, for openmp threads)
#...    #RJ TRY_NUM     = <$DEFAULT_TRY_NUM=default>
#...                      RJ_STOP if it takes less than TRY_TIME to do TRY_NUM
#...                                        tries.  Set to <=0 to turn off.
#...    #RJ TRY_NUM_MAX = <$DEFAULT_TRY_NUM_MAX=default>
#...                      RJ_STOP if this number of tries total is attempted in one submission.
#...                                       RJ_NEW_CHILD sent instead if RJ_RETRY seen
#...                       (unless previous $DEFAULT_TRY_NUM_MAX_NUM batches also had this error).
#...    #RJ TRY_TIME    = <$DEFAULT_TRY_TIME=default>
#...                      RJ_STOP if it takes less than TRY_TIME to do TRY_NUM
#...                                        tries.  Set to <=0 to turn off.
#...    #RJ U_<STRING>  = <nothing>
#...                      User defined variable that can be used elsewhere
#...                                        as a &&&RJ_VAR_U_<STRING>&&& variable.
#...                                        For example, "#RJ U_MYVAR = hello"
#...    #RJ UMASK       = <""=default>
#...                      Run with a different umask.  Useful when you
#...                                       need to allow others to run your problem as well
#...                                       (write access).
#...                                         7: group write, other has no access
#...                                        27: group read, other has no access
#...                                        GROUP and UMASK can be used together.
#...                                        eg: GROUP=dacodes, UMASK=7 => allow group write access
#...                                       
#...
#...  &&&RJ_<NAME>&&& replacement strings
#...  -----------------------------------
#...    Can be embedded in lines inside $RJ_FILE_DATAFILE.
#...
#...    &&&RJ_CMD_PRUN&&&       = replaced with the appropriate parallel run command for the
#...                              EXEC, NUMPE, etc.
#...                              For example, it might get replaced by:
#...                                 mpirun -np 4 my_exec.x my_input.in
#...    &&&RJ_VAR_<VAR>&&&      = Any variable in the above "#RJ Variables" section
#...                              &&&RJ_VAR_PNAME&&&, &&&RJ_VAR_U_MYVAR&&&, etc.
#...    &&&RJ_VAL_DIR&&&        = <$RJ_VAL_DIR>
#...                              directory where run script is located
#...    &&&RJ_VAL_DIR_ROOT&&&   = <$RJ_VAL_DIR_ROOT>
#...                              Root directory where run script is executed from
#...
#...  Predefined Environment Variables
#...  --------------------------------
#...    As well as any --env arguments, the following environment variables are set
#...    and updated at the BEGINNING of $RJ_FILE_DATAFILE each time it is run (and
#...    thus are available to the user):
#...       RJ_ID               = the current ID
#...       RJ_ID_CHILD         = ID of the child process
#...       RJ_RUN_JOB          = <yes> if currently being run from within run_job.pl
#...       RJ_TIME_START       = absolute time in seconds from start (date +%s)
#...       RJ_TIME             = initial run time in seconds
#...       RJ_TAG              = current tag being used
#...       RJ_TIME_CURRENT     = absolute time in seconds (date +%s)
#...       RJ_TIME_ELAPSED     = RJ_TIME_CURRENT - RJ_TIME_START
#...       RJ_TIME_REMAINING   = RJ_TIME - RJ_TIME_ELAPSED
#...       RJ_TIME_REMAINING_B = RJ_TIME - RJ_TIME_ELAPSED - buffer time (15m default)
#...       RJ_TRY              = number of times the script has been run in the same submission.
#...                             Useful if you expect to have multiple runs in the
#...                             same submission.
#...       RJ_VAL_DIR          = same variable as &&&RJ_VAL_DIR&&&
#...       RJ_VAL_DIR_ROOT     = same variable as &&&RJ_VAL_DIR_ROOT&&&
#...
#...    And the following machine-specific environment variables that are also set:
${sysinfo_stuff}#...
#...
#...  Special Notes
#...  -------------
#...  o Script initialization
#...    Depending upon how your environment initializes, you might need to have something
#...    like the following at the top of your $RJ_FILE_DATAFILE to correctly set up basic stuff
#...    like paths, modules, ... :
#...      #!/bin/csh -f
#...      if ( -f /etc/csh.cshrc ) then
#...        source /etc/csh.login
#...        source /etc/csh.cshrc
#...      endif
#...      module load foo
#...    To not change your PATH, you might just want to source the module.csh file
#...    (located in different places on different machines):
#...      foreach _loc ( /opt/modules/default/etc /etc/profile.d )
#...        if ( -e \$_loc/modules.csh ) then
#...          source \$_loc/modules.csh
#...          break
#...        endif
#...      end
#...
#...  o Tips for running multiple mpi jobs in the background
#...    - Set NUMPE and PPN to match the size of the allocation you will need
#...      Most mpi's can only run one job per node.
#...    - You might need to set PRUN_ARGS to tell mpi there will be multiple
#...      jobs running at once.
#...      You also may need to tell mpi which node to run on explicitly:
#...         -H tua015.localdomain,tua060,tua075
#...    - To put a job in the background, you should "wrap" it in "()"
#...      since &&&RJ_CMD_PRUN&&& might be replaces with various commands:
#...        ( &&&RJ_CMD_PRUN&&& >& prun.out ) &
#...    - Keep some defined set of jobs running at the same time automatically
#...
#...  o The run scripts will not submit a job if they detect another
#...    job running in the same directory.  Use the flag "--depend ." to
#...    start the job after the other job(s) have finished.
#...
#...Files Created
#...=============
#...  Several of these files are tagged with a data/try number.  The format is:
#...      .YYMMDDHHmmss<.try number>
#...
#...  The date tagged files are created in the directory '$RJ_DIR' with soft links to the
#...  latest file in the current working directory.
#...
#...  $RJ_FILE_CMD_OUT -> $RJ_DIR/$RJ_FILE_CMD_OUT.<pname>.<date><.try number>
#...    output from running $RJ_FILE_CMD_BASE
#...    Additional lines created by run_job.pl will start with:
#...      RJ_OUTPUT: 
#...  $RJ_FILE_CMD_BASE -> $RJ_FILE_CMD.<date>
#...    post processed $RJ_FILE_DATAFILE - this is the script that is actually run.
#...    It gets created immediately before it actually gets run.  So, you can
#...    modify $RJ_FILE_DATAFILE after you submit a job in case you need
#...    to change something.
#...  $RJ_FILE_BATCH_OUT_BASE -> $RJ_FILE_BATCH_OUT.<date>
#...    output from submitting/running $RJ_FILE_WRAPPER
#...  $RJ_FILE_TAG:
#...    contains current tag used in various file names.
#...  $RJ_FILE_TAG_OLD:
#...    contains previous tag used.
#...  $RJ_FILE_TAGS:
#...    contains batch tags used from current chain
#...  $RJ_FILE_WRAPPER:
#...    wrapper script used to launch job.
#...  $RJ_FILE_ID_NEXT:
#...    Contains the id (pid, moab job id, ... ) of the next job
#...  $RJ_FILE_LOG
#...    A log file of your jobs in your home directory.  Useful for seeing which
#...    jobs have been run or should be running.  Although attempts were made
#...    to make this file resistant to race conditions for the batch and file
#...    systems, this file might not be correct.
#...  ${RJ_FILE_MULTI}<multi file>.txt
#...    If jobs were launched via the --multi option, this file will be created
#...
#...Special Output or Files
#...=======================
#...  RJ_RETRY     : rerun $RJ_FILE_CMD_BASE again in same batch submission
#...  RJ_NEW_CHILD : stop current job, but let child job run
#...  RJ_STOP      : stop current job and kill child job
#...
#...  For priority: RJ_STOP > RJ_NEW_CHILD > RJ_RETRY
#...
#...  At the end of a run:
#...    o The file $RJ_FILE_CMD_OUT is scanned for these special strings
#...      Calling 'run_job_cleanup.pl --check' will scan the file for errors, completion
#...      status, various other project specific output and then echo the correct string
#...      See example below.
#...    o Any files that are named one of the special strings are removed and that flag is set.
#...      So, to stop a run, 'touch RJ_STOP' to create the file named RJ_STOP.
#...
#...  In $RJ_FILE_CMD_BASE, you may echo the string or create these special files as needed.
#...
#...Files/scripts used:
#...===================
#...  Makefile_vars.mk
#...  check_nodes.pl
#...  my_utils.pm
#...  read_output_files.pm
#...  read_output_files_stubs.pm
#...  run_job_cleanup.pl
#...  run_job_cleanup_extras.pm
#...  run_job.pl
#...  run_status.pl
#...  run_status_extras_stubs.pm
#...
#...Examples
#...========
#... o run_job.pl go
#...     runs the job
#... o run_job.pl --multi <datafile>
#...     Runs the jobs in the directories listed in the datafile.
#...     will not run a job in a directory that already has a job
#...     running in it.
#... o run_job.pl -i --numpe 128 --prun_args "-a -b" ./foo.x foo.in
#...     Run interactively EXEC=foo.x, NUMPE=128, EXEC_ARGS=foo.in, PRUN_ARGS="-a -b"
#...     If no $RJ_FILE_DATAFILE file, will use &&&RJ_CMD_PRUN&&& to effectively
#...     run something like:    mpirun -np 128 -a -b ./foo.x foo.in
#... o run_job.pl --numpe 4 --debugger tv ./foo.x foo.in
#...     Debug under totalview a 4 process foo.x run interactively (debugger sets -i)
#... o run_job.pl --depend .,123
#...     Start a job that depends on any other job running in the current
#...     directory as well as jobid=123 to finish first.
#...
#... Generating $RJ_FILE_DATAFILE and running:
#... -----------------------------------
#... o run_job.pl --convert job.params
#...     converts the job.params file into a $RJ_FILE_DATAFILE script
#...   run_job.pl go
#...     runs
#... o run_job.pl --convert basic
#...     answer questions to generate a $RJ_FILE_DATAFILE script
#...   run_job.pl go
#...     runs
#...
#...See Also
#...========
#... run_status.pl - shows your running/queue'd jobs
#... run_job_cleanup.pl - cleans file (mostly for EAP)
#............................................................................
EOF
    exit;
}

# if just running the wrapper
if( defined( $cmd{wrapper_file} ) ){
    &rj_run_wrapper( \%cmd );
    exit( $ierr );
}

# if converting, do that then exit
if( defined( $cmd{convert} ) ){
    &rj_convert( \%cmd );
    exit( $ierr );
}

# if launching a group, do that then exit
if( defined( $cmd{multi} ) ){
    &rj_multi( \%cmd );
    exit( $ierr );
}

# if doing special fix of environment
if( defined( $cmd{fix_env} ) ){
    $output = &rj_fix_env( \%cmd, $cmd{fix_env}[0], $cmd{fix_env}[1] );
    if( $output ne "" ){
        print "$output\n";
    }
    # need something there for eval
    else{
        print '/bin/true';
    }
    exit( $ierr );
}

# read datafile early - this might change command args (like depend)
&rj_datafile( \%datafile, \%cmd );
# stop if just getting info
if( defined( $cmd{info} ) ){
    $out = &rj_datafile_print( \%datafile );
    print "$out";
    exit( $ierr );
}

# stop if other jobs already here
# Skip the check if the user specified a "depend" option.
# This allows the user to start a child job if they accidentally
# deleted the previous child job.
if( ! defined($cmd{depend}) && ! defined($cmd{nodepend}) && ! defined($cmd{debug}) ){
    &rj_check_other_jobs( \%cmd );
}

# put some more info into RJ_FILE_ID
# if gotten to here, past checks of things and waiting
# run_status.pl looks at this file - but can get pid from ps
if( defined($cmd{id}) ){
    $out = "";
    $id = $cmd{sys_info}{batchid};
    if( $id ne "" ){
        $out .= "BATCHID=$id\n";
    }
    $out .= "PID=$$:PID\n";
    $out .= "LAUNCHTYPE=$cmd{launchtype}\n";
    $out .= "DATE=".`date +%s`;
    $out .= &rj_datafile_print( \%datafile );
    if( ! open( $fh_FILE, ">$RJ_FILE_ID" ) ){
        $ierr = 1;
        &print_error( "Cannot open id file [$RJ_FILE_ID]", $ierr );
        exit( $ierr );
    }
    flock( $fh_FILE, 2 );
    print $fh_FILE $out;
    close( $fh_FILE );
}

# set permissions
if( $datafile{var}{UMASK} =~ /\S/ ){
    umask( oct($datafile{var}{UMASK}) );
}
# create directory if not already there
if( ! -d $RJ_DIR ){
    my_mkdir( $RJ_DIR, $datafile{var}{GROUP}, &my_mode( DIR=>"", DEC=>"" ) );
}
# if given a tag, use that one first to make links, then get a new one
if( defined($cmd{tag}) ){
    # create soft link to RJ_FILE_BATCH_OUT
    `ln -sf ${RJ_FILE_BATCH_OUT}.$cmd{tag} ${RJ_FILE_BATCH_OUT_BASE}`;
    # create copy of RJ_FILE_DATAFILE (might not exist)
    if( -e $RJ_FILE_DATAFILE ){
        `cp -f $RJ_FILE_DATAFILE $RJ_DIR/$RJ_FILE_DATAFILE.$cmd{tag}`;
    }
}
# get a tag and put into file so other scripts/tools can use it
&new_tag( \%cmd, 0 );
# start message and info with new tag
if( defined($cmd{v}) ){
    &print_message( \%cmd, "Start", "tag $cmd{tag} $0" );
}
# submit child job
# Will not do if running interactive (so, it ignores CHAIN and
# BATCH settings).
# This is nice so that users will not have to change their current
#  $RJ_FILE_DATAFILE file to turn off CHAIN/BATCH settings
if( ! defined( $cmd{i} ) ){
    &rj_submit(\%datafile, \%cmd);
}
else{
    # must be -1 since job_killdepend kills submit_id
    $cmd{submit_id} = "-1";
    # still want to store this id into id file
    if( defined($cmd{v}) ){
    `echo "$cmd{submit_id}" > $RJ_FILE_ID_NEXT`;
    &print_message( \%cmd, "Dependent", "$cmd{submit_id}" );
    &print_message( \%cmd, "Finish_Submit", "" );
    }
    # still want to do a log
    &add_log( \%datafile, \%cmd );
}
# if this is not the first one ($cmd{id}), print command file
if( defined($cmd{id}) || defined($cmd{debug}) || defined($cmd{i}) ){
    &print_rj_cmd_file( \%datafile, \%cmd );
}
# verbose print
if( defined($cmd{v}) ){
    &print_perl_obj( \%cmd, "cmd_args" );
    &print_perl_obj( \%{$datafile{var}}, "datafile" );
}
# if not debug and
#  if this is not the first one ($cmd{id}) or running interactively,
#  go ahead and run
# Non-interactively, first one just does a submit
# Interactively, first one will also run
if( ! defined( $cmd{debug} ) &&
    ( defined( $cmd{id} ) || defined( $cmd{i} ) ) ) {
    $ierr = &rj_script( \%datafile, \%cmd );
}
# exit
if( defined($cmd{v}) ){
  &print_message( \%cmd, "Stop", "tag $cmd{tag} $0" );
}
exit( $ierr );

#######################################
#######################################
#######################################

################################################################################

sub rj_run_wrapper{
    my(
        $cmd_ref
        ) = @_;
    my(
        $cmd,
        $ierr,
        $pid_depend,
        @pids_depend_array,
        %pids_depend_hash,
        $pids_depend_string,
        $running_id,
        );

    $ierr = 0;

    if( defined($$cmd_ref{wrapper_wait}) ){
        ($pids_depend_string = $$cmd_ref{wrapper_wait}) =~ s/,/, /g;
        @pids_depend_array = split(/,/, $$cmd_ref{wrapper_wait});
        foreach $pid_depend ( @pids_depend_array ){
            $pids_depend_hash{$pid_depend} = "";
        }
        &print_message( $cmd_ref, "Waiting:Child", $pids_depend_string );
        &my_wait( $cmd_ref, IDS=>\%pids_depend_hash );
    }

    $running_id = $$;

    $cmd = "export RJ_RUNNING_ID=$running_id ; ($$cmd_ref{wrapper_file} 2>&1) >> $$cmd_ref{wrapper_out}";
    `$cmd`;

    return( $ierr );
}

################################################################################

# determine file/shell being used
# sets RJ_FILE_DATAFILE, RJ_SHELL
# might need to change in to allow:
#   run_job.foo to be shell "csh" or something
# for now, it is just the suffix
sub rj_getshell{
    my(
        $cmd_ref,   # in
        $dir,       # in
        ) = @_;
    my(
        $ierr,
        @files,
        $output,
        );

    # if shell set, that is file to use
    if( defined( $$cmd_ref{shell} ) ){
        $RJ_FILE_DATAFILE = "$dir/run_job$RJ_SHELL_SUFFIXES{$$cmd_ref{shell}}";
    }
    else{
        if( ! opendir( THISDIR, $dir ) ){
            $ierr = 1;
            &print_error( "Cannot read directory [$dir]",
                          $ierr );
            exit( $ierr );
        }
        @files = grep( /^run_job\.($RJ_SHELLS)$/, sort readdir( THISDIR ) );
        closedir( THISDIR );
        if( $#files >= 0 ){
            $RJ_FILE_DATAFILE = $files[0];
        }
        else{
            $RJ_FILE_DATAFILE = "$dir/run_job$RJ_SHELL_SUFFIXES{$RJ_SHELL}";
        }
    }

    if( $RJ_FILE_DATAFILE =~ /\.($RJ_SHELLS)$/ ){
        $RJ_SHELL    = $1;
    }
    else{
        $ierr = 1;
        &print_error( "Unrecognized run_job.<shell> name [$RJ_FILE_DATAFILE]",
                      $ierr );
        exit( $ierr );
    }

    # look at top of file for interpreter lines - that has final say
    if( -e $RJ_FILE_DATAFILE ){
        $RJ_SHELL="unknown";
        $output = `head $RJ_FILE_DATAFILE 2>&1`;
        $output =~ s/\s+$//;
        if( $output =~ /^\s*#!\s*(\S+)/m ){
            $RJ_SHELL = $1;
            $RJ_SHELL =~ s/^\S+\/(\S+?)\b.*$/$1/;
        }
        if( $RJ_SHELL !~ /^($RJ_SHELLS)$/ ){
            $ierr = 1;
            &print_error( "Unrecognized shell in file [$RJ_FILE_DATAFILE]",
                          "Expected to find something like:",
                          "   #!/bin/tcsh -f",
                          "Lines scanned:\n$output",
                          $ierr );
            exit( $ierr );
        }
    }

    $RJ_FILE_CMD_BASE       = "rj_cmd$RJ_SHELL_SUFFIXES{$RJ_SHELL}";
    $RJ_FILE_CMD            = "$RJ_DIR/$RJ_FILE_CMD_BASE";
    # for now, the wrapper file is a tcsh shell script
    $RJ_FILE_WRAPPER        = "$RJ_DIR/rj_wrapper$RJ_SHELL_SUFFIXES{tcsh}";
}

# launch a bunch of jobs from a directory list
sub rj_multi{
    my(
       $cmd_ref,
       ) = @_;
    my(
        $cmd,
        $dir,
        @dirs,
        @dirs_exist,
        %dirs_jobid,
        @dirs_launch,
        @dirs_skip,
        $ierr,
        $jobids,
        $line,
        $line_orig,
        $mail_to,
        %names_jobid,
        $output,
        $output_launch,
        $outfile,
        %run_status_info,
        %stat,
       );

    # get list of directories
    &my_stat( $$cmd_ref{multi}, \%stat );
    $outfile = "$stat{dir}/${RJ_FILE_MULTI}$stat{notdir}.txt";
    if( ! open( FILE, $$cmd_ref{multi} ) ){
        $ierr = 1,
        &print_error( "multi: Cannot open multi file [$$cmd_ref{multi}]", $ierr );
        exit( $ierr );
    }
    $i = 1;
    while( $line_orig = <FILE> ){
        ( $line = $line_orig ) =~ s/^\s*//;
        $line =~ s/\s*$//;
        $line =~ s/\s*#.*//;
        if( $line =~ /\S/ ){
            # get real directory
            undef( %stat );
            &my_stat( $line, \%stat );
            &rj_getshell( $cmd_ref, $stat{fullpath} );
            if( defined($stat{fullpath}) && -e $RJ_FILE_DATAFILE ){
                push( @dirs, $stat{fullpath} );
            }
            else{
                push( @dirs_skip, "$$cmd_ref{multi}:$i $line" );
            }
        }
        $i++;
    }
    close( FILE );

    # run run_status.pl to get jobs
    undef( %run_status_info );
    &get_run_status_info( \%run_status_info, \%dirs_jobid, \%names_jobid );

    # create dirs_exist and dirs_launch
    foreach $dir ( @dirs ){
        if( defined( $dirs_jobid{$dir} ) ){
            $jobids = join( ":\n   ", sort keys %{$dirs_jobid{$dir}} );
            push( @dirs_exist, "   $jobids: $dir" );
        }
        else{
            push( @dirs_launch, $dir );
        }
    }

    # launch jobs
    if( $#dirs_launch >= 0 ){
        &my_stat( $0, \%stat );
        # for now, add me to mail list
        $mail_to = "lmdm\@lanl.gov";
        if( defined($$cmd_ref{mail}) ){
            $mail_to .= ",$$cmd_ref{mail}";
        }
        $mail_to =~ s/^,//;
        foreach $dir( @dirs_launch ){
            if( defined($$cmd_ref{mail}) ){
                print "\nLaunching: $dir\n";
            }
            $cmd = "cd $dir && $stat{fullpath} go\n";
            $output_launch = &run_command( COMMAND=>$cmd, OUT_FILE=>$outfile, APPEND=>"" );
            if( ! defined($$cmd_ref{mail}) ){
                print $output_launch;
            }
        }
        if( $mail_to =~ /\S/ ){
            $cmd = "Mail -s 'rj_multi: $$cmd_ref{multi}' $mail_to < $outfile";
            $output = &run_command( COMMAND=>$cmd, VERBOSE=>"", STDOUT=>"" );
        }
    }

    # skip invalid directories
    if( $#dirs_skip >= 0 ){
        print "Skip (cannot stat or no $RJ_FILE_DATAFILE file):\n";
        foreach $dir( @dirs_skip ){
            print "   $dir\n";
        }
    }
    # existing jobs
    # do not print - exit quietly if "successful"
    #if( $#dirs_exist >= 0 ){
    #    print "Exists:\n";
    #    foreach $dir( @dirs_exist ){
    #        print "$dir\n";
    #    }
    #}
}

################################################################################

# get users job info and directory cross reference
sub get_run_status_info{
    my(
        $run_status_info_ref,
        $dirs_jobid_ref,
        $names_jobid_ref,
        $run_status_args,
        ) = @_;
    my(
        $cmd,
        $dir,
        $ierr,
        $jobid,
        $line,
        $line_orig,
        @lines,
        $name,
        $output,
        %stat,
        );
    $ierr = 0;
    if( ! defined($run_status_args) ){
        $run_status_args = "";
    }
    $cmd = "$RUN_STATUS_PL -f 2 $run_status_args";
    $output = &run_command( COMMAND=>$cmd, STATUS=>\$ierr );
    if( $ierr != 0 ){
        return( $ierr );
    }

    @lines = split( /\n/, $output );
    foreach $line_orig( @lines ){
        ($line = $line_orig) =~ s/^\s*//;
        $line =~ s/\s*$//;
        # assume output has JOBID first
        if( $line =~ /^JOBID\s*=\s*(\S.*)$/ ){
            $jobid = $1;
        }
        if( defined( $jobid ) && $line =~ /^(\S+)\s*=\s*(\S.*)$/ ){
            $$run_status_info_ref{$jobid}{vals}{$1} = $2;
        }
    }
    # create cross reference dirs
    foreach $jobid ( keys ( %$run_status_info_ref ) ){
        $dir = $$run_status_info_ref{$jobid}{vals}{dir};
        if( defined($dir) ){
            &my_stat( $dir, \%stat );
            if( defined($stat{fullpath}) ){
                $$dirs_jobid_ref{$stat{fullpath}}{$jobid} = "";
            }
        }
    }
    # create cross reference names
    foreach $jobid ( keys ( %$run_status_info_ref ) ){
        $name = $$run_status_info_ref{$jobid}{vals}{NAME};
        $$names_jobid_ref{$name}{$jobid} = "";
    }

    return( $ierr );
}

################################################################################

sub print_rj_cmd_file{
    my(
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
       $cmd_file,
       $ierr,
       );
    $ierr = 0;
    $cmd_file = "$RJ_FILE_CMD.$$cmd_ref{tag}";
    if( ! open( FILE, ">$cmd_file" ) ){
        $ierr = 1;
        &print_error( "Cannot open run_job.pl data file [$cmd_file]",
                      $ierr );
        exit( $ierr );
    }
    if( @{$$datafile_ref{new}} ){
        &print_message( $cmd_ref, "$RJ_FILE_CMD_BASE", "created ./$RJ_FILE_CMD_BASE->$cmd_file" );
        print FILE @{$$datafile_ref{new}};
    }
    close( FILE );
    unlink( $RJ_FILE_CMD );
    `ln -sf $cmd_file ./$RJ_FILE_CMD_BASE`;
    chmod &my_mode( EXEC=>"", DEC=>"" ), $cmd_file;
}

# Will return the batch id or $$ if not set
sub my_getid{
    my(
        $cmd_ref,
        ) = @_;
    my( $id );
    $id = $$cmd_ref{sys_info}{batchid};
    if( $id eq "" ){
        $id = $$;
    }
    return( $id );
}

sub add_log{
    my(
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
       $command,
       $date,
       $file,
       $userlist_dir,
       $userlist_file,
       $id,
       $ierr,
       $line,
       $output,
       );

    $ierr = 0;

    # do not use my_lockfile or locks

    # for checking if stuff is written correctly
    #$file = $RJ_FILE_LOG.".$$datafile_ref{var}{PNAME}.$$";
    $file = $RJ_FILE_LOG;
    if( ! -e $file ){
        `touch $file 2>&1`;
    }
    # set correct group for log file
    if( $$datafile_ref{var}{GROUP} =~ /\S/ ){
        `chgrp $$datafile_ref{var}{GROUP} $file 2>&1`;
    }
    # remove duplicate directories and only keep last few lines
    $command = "(fgrep -a -v ' RJ_VAL_DIR_ROOT:$RJ_VAL_DIR_ROOT ' $file | tail -n $RJ_FILE_LOG_SL) 2> /dev/null";
    $output = `$command`;
    $id = &my_getid( $cmd_ref );
    $date = `date`;
    chomp( $date );
    # make sure to have a space before and after values to make for easier grepping
    $line = " RJ_TAG:$$cmd_ref{tag} RJ_VAL_DIR_ROOT:$RJ_VAL_DIR_ROOT RJ_ID:$id RJ_ID_CHILD:$$cmd_ref{submit_id} RJ_L_MACHINE:$$cmd_ref{sys_info}{L_MACHINE} \n";
    if( open( FILE, ">$file" ) ){
        print FILE $output;
        print FILE $line;
        close( FILE );
    }
    else{
        $ierr = 0;
        &print_error( "Cannot open log file [$file]",
                      $ierr );
    }
}

# try to flush a file
sub my_flush{
    my( $file ) = @_;
    my( $i, @lines );
    # repeat
    for( $i = 0; $i < 2; $i++ ){
        `ls -la $file > /dev/null 2>&1`;
        `cat $file > /dev/null 2>&1`;
        if( open( FILE, "$file" ) ){
            @lines = <FILE>;
            close( FILE );
        }
        # hangs on some machines...sigh...
        # `sync`;
    }
}

#......................................................
#...create a datafile template from asking questions...
#......................................................
sub rj_convert{
    my(
       $cmd_ref,
       ) = @_;
    my(
       $ierr,
       $output,
       );
    $ierr = 0;
    if( -e $RJ_FILE_DATAFILE ){
        $ierr = 1;
        &print_error( "Output file [$RJ_FILE_DATAFILE] already exists.",
                      "Remove that file manually first.",
                      $ierr );
        exit( $ierr );
    }
    $output = "";
    if( $$cmd_ref{convert} =~ /\.test/ ){
        &rj_convert_cts( $cmd_ref, \$output );
    }
    elsif( $$cmd_ref{convert} eq "job.params" ){
        &rj_convert_job_params( $cmd_ref, \$output );
    }
    elsif( $$cmd_ref{convert} eq "basic" ){
        &rj_convert_basic( $cmd_ref, \$output );
    }
    else{
        $ierr = 1;
        &print_error( "Invalid conversion type [$$cmd_ref{convert}].",
                      $ierr );
        exit( $ierr );
        
    }
    # print RJ_FILE_DATAFILE
    if( ! open( FILE, ">$RJ_FILE_DATAFILE") ){
        $ierr = 1;
        &print_error( "Cannot open output datafile [$RJ_FILE_DATAFILE].",
                      $ierr );
        exit( $ierr );
    }
    print FILE $output;
    print "\n";
    print "=====================\n";
    print "New $RJ_FILE_DATAFILE file:\n";
    print "=====================\n";
    print $output;
    close( FILE );
}

sub rj_convert_cts{
    my(
       $cmd_ref,
       $output_ref
       ) = @_;
    my(
       $done,
       $ierr,
       $line,
       $line_num,
       @lines,
       $output,
       $output_header,
       $output_post,
       $output_pre,
       $output_rj,
       %set,
       %time,
       $time_string,
       );
    $ierr = 0;
    if( ! open( FILE, "$$cmd_ref{convert}") ){
        $ierr = 1;
        &print_error( "Cannot open file to convert [$$cmd_ref{convert}] does not exist.",
                      $ierr );
        exit( $ierr );
    }
    @lines = <FILE>;
    close( FILE );
    $output_pre = &rj_shell_header( SHELL=>$RJ_SHELL, QUICK=>"" );
    $output = "";
    $output_rj = "";
    $output_post = "";
    $output_pre = "";
    $done = "false";
    $i = 0;
    while( $done eq "false" ){
        if( $i > $#lines ){
            $done = "true";
            last;
        }
        $line_num = $i + 1;
        $line = $lines[$i]; chomp($line);
        if( $lines[$i] =~ /^(\s*)NUM_CPUS.*?(\d+)/ ){
            $output .= "${1}#RJ NUMPE = ${2}\n";
            $set{NUMPE} = "";
        }
        elsif( $lines[$i] =~ /^(\s*)TIME_LIMIT.*?(\d.*?)\s*$/ ){
            $time_string = $2;
            %time = &conv_time( STRING=>$time_string );
            $output .= "${1}#RJ TIME = $time{hms}\n";
            $set{TIME} = "";
        }
        elsif( $lines[$i] =~ /^(\s*)(tst_exec|product)(.*)$/ ){
            $output .= "${1}&&&RJ_CMD_PRUN&&&${3}\n";
        }
        elsif( $lines[$i] =~ /^\s*restart(\s*:\s*)?\s+.*restart\.csh/ ){
            $set{CHAIN} = "";
        }
        elsif( $lines[$i] =~ /^(\s*)touch\s+CTS_DO_NOT_RESTART/ ){
            $output .= "${1}echo RJ_STOP\n";
        }
        elsif( $lines[$i] =~ /^(\s*)(NAME|DESCRIPTION):?\s*(\S+)/ ){
            $output .= "#$lines[$i]";
        }
        else{
            $output .= "$lines[$i]";
        }      
        $i++;
    }
    if( ! defined($set{NUMPE}) ){
        $output_rj .= "#RJ NUMPE = $DEFAULT_NUMPE\n";
    }
    if( ! defined($set{TIME}) ){
        $output_rj .= "#RJ TIME = 1h\n";
    }
    if( defined($set{CHAIN}) ){
        $output_rj .= "#RJ CHAIN = yes\n";
    }
    # not the default to run cleanup afterwards
    #$output .= "&&&RJ_VAL_DIR&&&/run_job_cleanup.pl --check\n";

    $output_header = &rj_shell_header( SHELL=>$RJ_SHELL, 
                                       QUICK=>"", 
                                       MISC=>$output_rj,
                                       SYSTEM=>$set{sys_setup} );

    # print translation
    print "========================\n";
    print "Previous $$cmd_ref{convert} file\n";
    print "========================\n";
    print @lines;
    $$output_ref = $output_header . $output . $output_post;
    return( $ierr );


}

sub rj_convert_job_params{
    my(
       $cmd_ref,
       $output_ref
       ) = @_;
    my(
       $args,
       $comment,
       $done,
       $dotfile,
       $file,
       $ierr,
       $line,
       $line_num,
       @lines,
       @old_lines,
       $output,
       $output_header,
       $output_line,
       @output_lines,
       $output_pre1,
       $output_pre2,
       $output_rj,
       %set,
       $try,
       $val,
       $var,
       );
    $ierr = 0;
    if( ! open( FILE, "$$cmd_ref{convert}") ){
        $ierr = 1;
        &print_error( "Cannot open file to convert [$$cmd_ref{convert}] does not exist.",
                      $ierr );
        exit( $ierr );
    }
    @lines = <FILE>;
    close( FILE );
    $output_rj    = "";
    $output_pre1  = "";
    $output_pre2  = "";
    $output       = "";
    # get problem name
    # do not set pname when converting job.params.
    # should be able to determine pname from input file.
    #$val = &get_pname(".");
    #$comment = "get_pname()";
    #$output_rj .= sprintf( "#RJ %-10s = $val\n", "PNAME" );
    $done = "false";
    @old_lines = ();
    $i = 0;
    if( $RJ_SHELL eq "bash" || $RJ_SHELL eq "sh" ){
        $dotfile = ".bashrc";
    }
    else{
        $dotfile = ".cshrc";
    }
    while( $done eq "false" ){
        $comment = "";
        if( $i > $#lines ){
            $done = "true";
            last;
        }
        $line_num = sprintf( "%03d", $i + 1 );
        $line = $lines[$i]; chomp($line);
        # blank line
        if( $line !~ /\S/ ){
            $comment = "skip";
        }
        # set <var> = <val>
        elsif( $line =~ /^\s*set\s*(\S+)\s*=\s*(.*)$/ ){
            $var = $1;
            $val = $2;
            if( ! defined($val) || $val !~ /\S/ ){
                $val = "";
            }
            $set{$var} = $val;
            if( $val !~ /\S/ ){
                $comment = "skip - NULL value";
            }
            elsif( $var eq "allow_group_access" ){
                $comment = "GROUP";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "GROUP" );
            }
            elsif( $var eq "batch_software" ){
                if( $val eq "NoBatch" ){
                    $comment = "BATCH=no";
                    $output_rj .= sprintf( "#RJ %-10s = no\n", "BATCH" );
                }
                else{
                    $comment = "skip - BATCH=yes default";
                }
            }
            elsif( $var eq "batch_time" ){
                $comment = "TIME";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "TIME" );
            }
            elsif( $var eq "batch_time_buffer" ){
                $comment = "TIME_B";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "TIME_B" );
            }
            elsif( $var eq "binary" ){
                $comment = "EXEC";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "EXEC" );
            }
            elsif( $var eq "depend" ){
                if( $val =~ /^y/i ){
                    $val = "yes";
                    $comment = "CHAIN";
                    $output_rj .= sprintf( "#RJ %-10s = $val\n", "CHAIN" );
                }
                else{
                    $comment = "skip - CHAIN default";
                }
            }
            elsif( $var eq "directory" ){
                $comment = "cd <directory>";
                $output_pre2 .= "# directory\n";
                $output_pre2 .=
                    &rj_shell_if( "if", "! -d $val",
                                  "  echo 'RJ_STOP - directory [$val] does not exist.'\n".
                                  "  exit 1",
                                  "endif", "", $RJ_SHELL );
                $output_pre2 .= "cd $val\n";
            }
            elsif( $var eq "epilogue" ){
                $comment = "code block inserted";
            }
            elsif( $var eq "hdfs_by_lab" ){
                if( $val =~ /^y/i ){
                    $ierr = 1;
                    &print_error( "run_job_cleanup.pl cannot handle this",
                                  "Contact crestone_support\@lanl.gov to request",
                                  "$line_num: $line",
                                  $ierr );
                    exit( $ierr );
                }
                else{
                    $comment = "skip - run_job_cleanup.pl default";
                }
            }
            elsif( $var eq "ignore_fatal_errors" ){
                $comment = "skip - can no longer be set";
            }
            elsif( $var eq "input" ){
                $comment = "RJ_CMD_PRUN line";
                @output_lines = split( /\n/, `grep secmax $val 2>&1` );
                foreach $output_line ( @output_lines ){
                    if( $output_line =~ /^\s*secmax\s*=\s*\$(\S+)/ ){
                        $set{secmax} = $1;
                    }
                }
                # set it to just "secmax" - users should use this
                if( ! defined($set{secmax}) ){
                    $set{secmax} = "secmax";
                }
            }
            elsif( $var eq "kill_depend_job" ){
                $comment = "skip - always kill depend";
            }
            elsif( $var eq "min_run_time" ){
                $comment = "skip - no longer used";
            }
            elsif( $var eq "modules" ){
                $comment = "use EAP $dotfile";
                $file = "";
                $try = "$RJ_VAL_DIR/../dotfiles/$dotfile";
                if( $file !~ /\S/ && -e $try ){
                    $file = $try;
                }
                $try = "$RJ_VAL_DIR/../../Tools.rh/Environment/$dotfile";
                if( $file !~ /\S/ && -e $try ){
                    $file = $try;
                }
                if( $file =~ /\S/ ){
                    $output_pre1 .= "# source the eap $dotfile to get modules and stuff\n";
                    $output_pre1 .= "source $file\n";
                }
            }
            elsif( $var eq "module_use_dir" ){
                $comment = "skip - use eap $dotfile";
            }
            elsif( $var eq "numpe" ){
                $comment = "NUMPE";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "NUMPE" );
            }
            elsif( $var eq "procs_per_node" ){
                $comment = "PPN";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "PPN" );
            }
            elsif( $var eq "program_args" ){
                $comment = "RJ_CMD_PRUN line";
            }
            elsif( $var eq "prologue" ){
                $comment = "code block inserted";
            }
            elsif( $var eq "teos_in" ){
                $comment = "code block inserted";
                @output_lines = reverse split( /\n/, `grep teos_file $val 2>&1` );
                $set{teos_out} = "teos.out";
                foreach $output_line ( @output_lines ){
                    if( $output_line =~ /^\s*teos_file\s*=\s*(\'|\")(\S+)(\"|\')/ ){
                        $set{teos_out} = $2;
                        last;
                    }
                }
            }
            elsif( $var eq "xtra_batch_opts" ){
                $comment = "BATCH_ARGS";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "BATCH_ARGS" );
            }
            elsif( $var eq "xtra_mpirun_opts" ){
                $comment = "PRUN_ARGS";
                $output_rj .= sprintf( "#RJ %-10s = $val\n", "PRUN_ARGS" );
            }
            elsif( $var eq "nodes" ||
                   $var eq "node_resources" ||
                   $var eq "queue" ||
                   $var eq "perf_data_save_dir"
                   ){
                # NOTE: to allow "nodes" to be set, it makes it much harder
                #       to determine how much of nodes/ppn resource should
                #       be requested (when to check for min/max numpe,ppn,
                #       nodes that might be redefined during the $RJ_FILE_DATAFILE.
                #       So, do not allow setting of nodes.
                $ierr = 1;
                &print_error( "Not implemented.",
                              "Contact crestone_support\@lanl.gov to request",
                              "$line_num: $line",
                              $ierr );
                exit( $ierr );
            }
            # unknown line
            else{
                $ierr = 1;
                &print_error( "job.params line not parsed:",
                              "$line_num: $line",
                              $ierr );
                exit( $ierr );
            }
        }
        # setenv
        elsif( $line =~ /^\s*setenv\s+(\S+)\s+(\S.*?)\s*$/ ){
            $comment = "setenv coppied";
            $output_pre1 .= &rj_shell_var_set($1, $2, "\n", $RJ_SHELL);
        }
        # comment
        elsif( $line =~ /^\s*(\#.*?)\s*$/ ){
            $comment = "skip";
        }
        # set <var> =
        elsif( $line =~ /^\s*set\s*(\S+)\s*=\s*$/ ) {
            $var = $1;
            if( $var eq "batch_software" ){
                $comment = "obsolete hack";
            }
        }
        # unknown line
        if( $comment !~ /\S/ ){
            $ierr = 1;
            &print_error( "job.params line not parsed:",
                          "$line_num: $line",
                          $ierr );
            exit( $ierr );
        }
        push( @old_lines, sprintf( "%-40s: %s: %-s\n", $comment, $line_num, $line ) );
        $i++;
    }
    # do prologue at end of output_pre
    # (after any setenv/cd <directory>/... have been done)
    if( defined($set{prologue}) ){
        $output_pre2 .= "# prologue\n$set{prologue}\n";
    }
    $output .= "# cleanup from any previous run\n";
    $output .= "&&&RJ_VAL_DIR&&&/run_job_cleanup.pl\n";
    if( defined($set{teos_in}) ){
        $output .= "# generate teos_out file\n";
        $output .= "if ( -e \"$set{teos_in}\" && ! -e \"$set{teos_out}\" ) then\n";
        $output .= "  &&&RJ_CMD_PRUN&&& $set{teos_in}\n";
        $output .= "endif\n";
    }
    $args = "";
    if( defined($set{input}) ){
        $args .= "$set{input}";
    }
    if( defined($set{program_args}) ){
        $args .= " $set{program_args}";
    }
    if( defined($set{secmax}) && $set{secmax} =~ /\S/ ){
        $args .= " -v $set{secmax}=\$RJ_TIME_REMAINING_B";
    }
    $output .= "# run\n";
    $output .= "&&&RJ_CMD_PRUN&&& $args\n";
    $output .= "# cleanup after this run and check results for failure/stopping\n";
    $output .= "&&&RJ_VAL_DIR&&&/run_job_cleanup.pl --check\n";
    if( defined( $set{epilogue} ) ){
        $output .= "# run epilogue if job will not continue (return==4 from run_job_cleanup.pl)\n";
        $output .= "if ( \$? == 4 ) then\n";
        $output .= "  $set{epilogue}\n";
        $output .= "endif\n";
    }

    $output_header = &rj_shell_header( SHELL=>$RJ_SHELL, 
                                       QUICK=>"", 
                                       MISC=>$output_rj,
                                       SYSTEM=>$set{sys_setup} );

    # print translation
    print "========================\n";
    print "Previous job.params file\n";
    print "========================\n";
    print @old_lines;
    $$output_ref = $output_header . $output_pre1 . $output_pre2 . $output;
    return( $ierr );
}

sub rj_convert_basic{
    my(
       $cmd_ref,
       $output_ref
       ) = @_;
    my(
       $answer,
       $cmd,
       $default,
       $dotfile,
       $dir,
       $dir_try,
       $file,
       @files,
       $ierr,
       $ifile,
       $input,
       $out,
       $output,
       $output1,
       $output_header,
       $output_line,
       @output_lines,
       $output_rj,
       $pname,
       $prompt,
       %set,
       %stat,
       $try,
       $type,
       $var,
       $proj,
       );
    # init
    $ierr = 0;
    $input = "";
    $output      = "";
    $$output_ref = "";
    $PS = "%45s";
    print "\nAnswer the following questions to generate $RJ_FILE_DATAFILE\n\n";
    if( $RJ_SHELL eq "bash" || $RJ_SHELL eq "sh" ){
        $dotfile = ".bashrc";
    }
    else{
        $dotfile = ".cshrc";
    }
    # getval
    $proj = "eap";
    @files = glob( "*.flg" );
    if( $#files >= 0 ){
        $proj = "lap";
    }
    $prompt  = "Project?";
    $var     = "proj";
    $type    = "string";
    $default = $proj;
    &my_getval( PROMPT=>$prompt, DEFAULT=>"$default", REGEXP=>'/^(eap|lap)$/', TYPE=>$type, VAR=>\$set{$var}, BLANK=>"" );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    $proj = $set{$var};
    # getval
    if( $proj eq "lap" ){
        $prompt  = "OPUS directory?";
        $default = "./OPUS";
        $dir_try = "/usr/projects/shavano/SPUBLIC/releases";
        if( -d $dir_try ){
            $cmd = "(ls -1d $dir_try/[1-9]*/OPUS) 2> /dev/null";
            $out = &run_command(COMMAND=>$cmd);
            if( $out =~ /([^\n]*)$/ ){
                $default = $1;
            }
        }
        $var     = "opus";
        $type    = "dir";
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    if( $proj eq "lap" ){
        $prompt  = "Suite?";
        $var     = "suite";
        $type    = "string";
        if( $$cmd_ref{sys_info}{L_CLASS} =~ /TLCC/ ){
            $default = "TLCCfast";
        }
        else{
            $default = "LUNAfast";
        }
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    $pname = "";
    $ifile = "";
    if( $proj eq "eap" ){
        @files = glob( "*.in*" );
        foreach $file ( @files ){
            $out = `grep pname $file 2>&1`;
            if( $out =~ /^\s*pname\s*=\s*[\'\"]\s*(\S+)\s*[\"\']/m ){
                $pname = $1;
                $ifile = $file;
                last;
            }
        }
    }
    elsif( $proj eq "lap" ){
        @files = glob( "*.flg" );
        if( $#files >= 0 ){
            $ifile = $files[0];
            ( $pname = $ifile ) =~ s/\.flg$//;
        }
    }
    if( $ifile !~ /\S/ ){
        @files = glob( "*.in*" );
        if( $#files >= 0 ){
            $ifile = $files[0];
        }
    }
    if( $ifile !~ /\S/ ){
        $ifile = getcwd();
        $ifile = &my_notdir( $ifile );
        $ifile = "$ifile.in";
    }
    $prompt  = "Input File?";
    $var     = "ifile";
    $type    = "string";
    $default = $ifile;
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var}, BLANK=>"" );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    if( $proj eq "eap" && $set{$var} =~ /\S/ ){
        $out = `grep pname $set{$var} 2>&1`;
        if( $out =~ /^\s*pname\s*=\s*[\'\"]\s*(\S+)\s*[\"\']/m ){
            $pname = $1;
        }
    }
    if( $pname !~ /\S/ ){
        ( $pname = $set{ifile} ) =~ s/\..*//;
    }
    if( $pname !~ /\S/ ){
        $pname = getcwd();
        $pname = &my_notdir( $pname );
    }
    if( $proj ne "eap" ){
        $prompt  = "Problem Name?";
        $var     = "PNAME";
        $type    = "string";
        $default = $pname;
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    if( $proj eq "lap" ){
        $dir_try = "$set{opus}/$set{suite}";
        $cmd = "( cd $dir_try && find -L . -type f -maxdepth 1) 2> /dev/null";
        $out = &run_command(COMMAND=>$cmd);
        print "Executables found in $dir_try:\n$out\n";
        $prompt  = "Executable (relative to $dir_try)?";
        $dir = "$set{opus}/$set{suite}";
    }
    else{
        $prompt  = "Executable?";
        $dir = ".";
    }
    $var     = "EXEC";
    $type    = "string";
    $default = "myprog.x";
    if( $proj eq "eap" ){
        $default = `ls -1 $ENV{L_EAP_INSTALL_DIR}/releases/latest/xrage/xrage*_$ENV{L_EAP_OS}*.x 2> /dev/null`;
        if( ! defined( $default ) || $default !~ /\S/ ){
            $default = "myprog.x";
        }
        $default =~ s/\s+//;
    }
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var}, DIR=>$dir, BLANK=>"" );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    $prompt  = "Chain jobs until finished?";
    $var     = "CHAIN";
    $type    = "yes/no";
    $default = "yes";
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    $var     = "STOP_FILE";
    $set{$var} = "";
    if( $proj ne "eap" ){
        $prompt  = "Name of stopping file?";
        $type    = "string";
        if( $proj eq "lap" ){
            $set{$var} = "HALT-AUTO-RESTART";
        }
        else{
            $default = "DONE";
            &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var}, BLANK=>"" );
        }
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    $prompt  = "Number of processes?";
    $var     = "NUMPE";
    $type    = "int";
    $default = $DEFAULT_NUMPE;
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    $prompt  = "Set number of processes per node?";
    $var     = "PPN_set";
    $type    = "yes/no";
    $default = "no";
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    if( $set{$var} eq "yes" ){ 
        $prompt  = "Processes per node?";
        $var     = "PPN";
        $type    = "int_frac";
        $default = $cmd{sys_info}{L_PPN};
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    $prompt  = "Set number of threads per process?";
    $var     = "TPP_set";
    $type    = "yes/no";
    $default = "no";
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    if( $set{$var} eq "yes" ){ 
        $prompt  = "Threads per process?";
        $var     = "TPP";
        $type    = "int";
        $default = $cmd{sys_info}{L_TPP};
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    $prompt  = "Time?";
    $var     = "TIME";
    $type    = "time";
    $default = $DEFAULT_TIME;
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    $prompt  = "Batch?";
    $var     = "BATCH";
    $type    = "never|auto|yes|no";
    $default = "yes";
    &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
    $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    if( $set{$var} eq "yes" ){
        $prompt  = "Set batch arguments?";
        $var     = "";
        $type    = "yes/no";
        $default = "no";
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
        # getval
        if( $set{$var} eq "yes" ){ 
            $prompt  = "Batch arguments?";
            $var     = "BATCH_ARGS";
            $type    = "string";
            $default = "-A access";
            &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
            $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
        }
    }
    # getval
    #$prompt  = "Put jobs in background?";
    #$var     = "bg";
    #$type    = "yes/no";
    #$default = "no";
    #&my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, REGEXP=>'/^no|not implemented yet$/',
    #            VAR=>\$set{$var} );
    #$input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    # getval
    &my_getval( PROMPT=>"Allow unix group access?", DEFAULT=>"yes", TYPE=>"yes/no", VAR=>\$answer );
    if( $answer eq "yes" ){
        $prompt  = "Unix group?";
        $var     = "GROUP";
        $type    = "string";
        &my_stat( ".", \%stat );
        $default = $stat{group};
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # eap stuff
    if( $proj eq "eap" ){
        &my_getval( PROMPT=>"Generate teos file?", DEFAULT=>"yes", TYPE=>"yes/no", VAR=>\$answer );
        if( $answer eq "yes" ){
            # getval
            @output_lines = reverse split( /\n/, `grep eosfile $set{ifile} 2>&1` );
            $default = "teos.out";
            foreach $output_line ( @output_lines ){
                if( $output_line =~ /^\s*eosfile\s*=\s*(\'|\")\s*(\S+)\s*(\"|\')/ ){
                    $default = $2;
                    last;
                }
            }
            $prompt  = "Name of teos.out file?";
            $var     = "teos_out";
            $type    = "string";
            &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
            $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
            # getval
            $prompt    = "Name of teos.in file?";
            $var       = "teos_in";
            $type      = "string";
            $default   = "teos.in";
            $output_line = `grep teos_file *in* 2>&1 | fgrep '$set{teos_out}'`;
            if( $output_line =~ /^(\S+):\s*teos_file\s*=\s*(\'|\")\s*(\S+)\s*(\"|\')/m ){
                $default = $1;
            }
            &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
            $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
        }
        if( $set{ifile} =~ /\S/ ){
            $out = `grep secmax $set{ifile} 2>&1`;
        }
        else{
            $out = "";
        }
        if( $out =~ /secmax\s*=\s*\$(\S+)/ ){
            $set{secmax} = $1;
        }
        else{
            $set{secmax} = "secmax";
        }
        $prompt  = "Set secmax command line var?";
        $var     = "secmax";
        $type    = "string";
        $default = $set{secmax};
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var}, BLANK=>"" );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # lap stuff
    if( $proj eq "lap" ){
        # getval
        $prompt  = "Name of mod file (if any)?";
        $var     = "mod";
        $type    = "file";
        $default = "";
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var}, BLANK=>"" );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    # getval
    $set{sys_setup} = "";
    if( $proj ne "eap" && $proj ne "lap" ){
        $prompt  = "Source system setup files?";
        $var     = "sys_setup";
        $type    = "yes/no";
        $default = "no";
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
    }
    if( $proj eq "lap" ){
        $set{sys_setup} = "yes";
    }
    if( $set{sys_setup} ne "yes" ){
        delete( $set{sys_setup} );
    }

    # getval
    $var       = "rst_cmd";
    $set{$var} = "";
    if( $proj eq "" ){
        $prompt  = "Different run command for restart?";
        $type    = "yes/no";
        $default = "yes";
        &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
        $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
        if( $set{$var} eq "yes" ){
            $prompt  = "Command args for restart?";
            $var     = "rst_args";
            $type    = "string";
            $default = "$set{ifile}";
            &my_getval( PROMPT=>$prompt, DEFAULT=>$default, TYPE=>$type, VAR=>\$set{$var} );
            $input .= sprintf( "${PS} (%10s) = [%s]\n", $prompt, $var, $set{$var} );
        }
    }
    # create file
    $output_rj .= sprintf( "#RJ %-10s = %s\n", "BATCH", $set{BATCH} );
    if( defined($set{BATCH_ARGS}) ){
        $output_rj .= sprintf( "#RJ %-10s = %s\n", "BATCH_ARGS", $set{BATCH_ARGS} );
    }
    $output_rj .= sprintf( "#RJ %-10s = %s\n", "CHAIN", $set{CHAIN} );
    if( $proj eq "lap" && $set{EXEC} !~ m&^/& ){
        $output_rj .= sprintf( "#RJ %-10s = \$OPUS/$set{suite}/%s\n", EXEC, $set{EXEC} );
    }
    else{
        $output_rj .= sprintf( "#RJ %-10s = %s\n", "EXEC", $set{EXEC} );
    }
    if( defined($set{GROUP}) ){
        $output_rj .= sprintf( "#RJ %-10s = %s\n", "GROUP", $set{GROUP} );
    }
    $output_rj .= sprintf( "#RJ %-10s = %s\n", "NUMPE", $set{NUMPE} );
    if( defined($set{PPN}) && $set{PPN} =~ /\S/ ){
        $output_rj .= sprintf( "#RJ %-10s = %s\n", "PPN", $set{PPN} );
    }
    if( defined($set{PNAME}) && $set{PNAME} =~ /\S/ ){
        $output_rj .= sprintf( "#RJ %-10s = %s\n", "PNAME", $set{PNAME} );
    }
    if( defined($set{TIME}) && $set{TIME} =~ /\S/ ){
        $output_rj .= sprintf( "#RJ %-10s = %s\n", "TIME", $set{TIME} );
    }
    if( $set{STOP_FILE} =~ /\S/ ){
        $output .= "# stop if STOP_FILE is there\n";
        $output .=
            &rj_shell_if( "if", "-e $set{STOP_FILE}",
                          "  echo $set{STOP_FILE} exists\n".
                          "  echo RJ_STOP\n".
                          "  exit",
                          "endif", "", $RJ_SHELL );
    }

    # source team dotfile
    if( $proj eq "eap" ){
        $file = "";
        $try = "$RJ_VAL_DIR/../dotfiles/$dotfile";
        if( $file !~ /\S/ && -e $try ){
            $file = $try;
        }
        $try = "$RJ_VAL_DIR/../../Tools.rh/Environment/$dotfile";
        if( $file !~ /\S/ && -e $try ){
            $file = $try;
        }
        if( $file =~ /\S/ ){
            $output .= "# source the eap $dotfile to get modules and stuff\n";
            $output .= "source $file\n\n";
        }
    }
    if( $proj eq "" ){
        $output .= "# cleanup from previous run\n";
        $output .= "run_job_cleanup.pl\n";
        $output .= "# run\n";
        if( $set{rst_cmd} eq "yes" ){
            $output .= "# if doing a restart\n";
            $output .=
                &rj_shell_if( "if", "-e my_restart_file",
                              "  &&&RJ_CMD_PRUN&&& $set{rst_args}",
                              "else", "", $RJ_SHELL );
            $output .= "# or running from scratch\n";
            $output .= 
                &rj_shell_if( "else", "",
                              "  &&&RJ_CMD_PRUN&&& $set{ifile}",
                              "endif", "", $RJ_SHELL );
        }
        else{
            $output .= "&&&RJ_CMD_PRUN&&& $set{ifile}\n";
        }
        $output .= "run_job_cleanup.pl --check\n";
    }

    if( $proj eq "eap" ){
        $output .= "# cleanup from previous run\n";
        $output .= "run_job_cleanup.pl\n";
        if( defined( $set{teos_in} ) ){
            $output .= "# generate teos_out file\n";
            $output .=
                &rj_shell_if("if", "-e \"$set{teos_in}\" && ! -e \"$set{teos_out}\"",
                             "  &&&RJ_CMD_PRUN&&& $set{teos_in}",
                             "endif", "", $RJ_SHELL );
        }
        $output .= "# run\n";
        $output .= "&&&RJ_CMD_PRUN&&& $set{ifile}";
        if( defined( $set{secmax} ) && $set{secmax} =~ /\S/ ){
            $output .= " -v $set{secmax}=\$RJ_TIME_REMAINING_B";
        }
        $output .= "\n";
        $output .= "# cleanup and check for continue\n";
        $output .= "run_job_cleanup.pl --check\n";
    }
    if( $proj eq "lap" ){
        $output .= &rj_shell_var_set("OPUS", $set{opus}, "\n", $RJ_SHELL);
        $output .= "# source setup file\n";
        $output .= "source $set{opus}/SHARED/config/$set{suite}.modules\n";
        $output .= "# softlink 000 file\n";
        $output .=
            &rj_shell_if("if", "! -e $set{EXEC}.000",
                         "  ln -s $set{opus}/$set{suite}/$set{EXEC}.000 .",
                         "endif", "", $RJ_SHELL );
        $output .= "# cleanup from previous run\n";
        $output .= "run_job_cleanup.pl --proj lap\n";
        $output .= "set RSTMOD='$set{mod}'\n";
        $output .= "# if restarting\n";
        $output1 = &rj_shell_if( "if", "-f \"\$RSTMOD\" && `\\ls -1rt Last_restart_name \$RSTMOD | tail -n 1` != Last_restart_name",
                                 "    set ARGS_MOD=\"\$RSTMOD\"",
                                 "else", "  ", $RJ_SHELL );
        $output1 .= &rj_shell_if( "else", "",
                                  "    set ARGS_MOD=''",
                                  "endif", "  ", $RJ_SHELL );
        $output .=
            &rj_shell_if( "if", "-f Last_restart_name",
                          "  # use mod file if it is newer than latest restart\n".
                          $output1.
                          "  &&&RJ_CMD_PRUN&&& open `cat Last_restart_name` \$ARGS_MOD",
                          "else", "", $RJ_SHELL );
        $output .= "# if not restarting\n";
        $output .=
            &rj_shell_if( "else", "",
                          "  &&&RJ_CMD_PRUN&&& include $set{ifile}",
                          "endif", "", $RJ_SHELL );
        $output .= "# cleanup and check for continue\n";
        $output .= "run_job_cleanup.pl --proj lap --check\n";
    }

    if( $set{STOP_FILE} =~ /\S/ ){
        $output .= "# stop if STOP_FILE is there\n";
        $output .=
            &rj_shell_if( "if", "-e $set{STOP_FILE}",
                          "  echo $set{STOP_FILE} exists\n".
                          "  echo RJ_STOP\n".
                          "  exit",
                          "endif", "", $RJ_SHELL );
    }

    $output_header = &rj_shell_header( SHELL=>$RJ_SHELL, 
                                       QUICK=>"", 
                                       MISC=>$output_rj,
                                       SYSTEM=>$set{sys_setup} );

    # print translation
    print "======\n";
    print "Input:\n";
    print "======\n";
    print $input;
    $$output_ref = $output_header . $output;
    return( $ierr );
}

sub new_tag{
    my(
       $cmd_ref,
       $try,
       ) = @_;
    my(
       $out,
       $output,
       );
    $$cmd_ref{try} = $try;
    # get the base tag name
    if( ! defined($$cmd_ref{tag}) ){
        $$cmd_ref{tag} = &date_ymdhms();
    }
    # get the whole tag name (with try suffix)
    $$cmd_ref{tag_whole} = $$cmd_ref{tag};
    if( $try > 0 ){
        $$cmd_ref{tag_whole} = sprintf( "%s.%03d", $$cmd_ref{tag}, $try );
    }
    # if you are going to run script (id), get the old whole tag and put into $RJ_FILE_TAG_OLD
    # if this is the first time run, do not push the stack down
    if( defined($$cmd_ref{id}) ){

        # push onto RJ_FILE_TAGS
        if( $$cmd_ref{tag_whole} =~ /^\d+$/ ){
            if( ! -e $RJ_FILE_TAGS ){
                $out = `touch $RJ_FILE_TAGS`;
            }
            $out = `echo $$cmd_ref{tag} >> $RJ_FILE_TAGS`;
        }
        $out = `cat $RJ_FILE_TAGS`;
        $out =~ s/^\s+//;
        $out =~ s/\s+$//;
        @{$$cmd_ref{tags}} = split(/\s+/, $out);

        if( open( FILE, "$RJ_FILE_TAG" ) ){
            $$cmd_ref{tag_whole_old} = <FILE>;
            close( FILE );
            chomp( $$cmd_ref{tag_whole_old} );
            if( $$cmd_ref{tag_whole_old} =~ /\d{14}(\.\d{3})?/ ){
                $output = `echo $$cmd_ref{tag_whole_old} > $RJ_FILE_TAG_OLD`;
            }
            else{
                unlink( $RJ_FILE_TAG_OLD );
                undef( $$cmd_ref{tag_whole_old} );
            }
        }
        # store whole tag into $RJ_FILE_TAG
        $output = `echo $$cmd_ref{tag_whole} > $RJ_FILE_TAG`;
    }

    # if first time run, remove $RJ_FILE_TAGS - consider fresh chain
    else{
        unlink( $RJ_FILE_TAGS );
    }
};

sub print_message{
    my(
       $cmd_ref,
       $action,
       $msg,
       ) = @_;
    my(
        $tag,
        );
    if( ! defined($msg) ){
        $msg = "";
    }
    $tag = $$cmd_ref{tag};
    if( ! defined( $tag ) ){
        $tag = "undef";
    }
    printf( STDERR "run_job.pl %14s %7s %s %-20s : %s\n",
            $tag, $$, &date_ymdhms_sep(), $action, $msg );
}

#...................................................
#...read the datafile and create the command file...
#...................................................
sub rj_datafile{
    my(
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
        $cond,
        $default,
        @errors,
        @error_lines,
        $exn,
        $ierr,
        $indent,
        @jobids_other,
        $key,
        $line,
        $line_new,
        $line_num,
        @lines,
        $op,
        $output_header,
        $pname,
        $pre,
        $replace,
        $rest,
        $rj_cmd,
        %stat,
        $val,
        $val_new,
        @vals,
        $var,
        );

    # for non-interactive jobs, datafile must exist
    if( ! defined( $$cmd_ref{exec_defined} ) && ! -e "$RJ_FILE_DATAFILE" ){
        $ierr = 1;
        &print_error( "Cannot open run_job.pl data file [$RJ_FILE_DATAFILE]",
                      $ierr );
        exit( $ierr );
    }
    undef( @lines );
    # if datafile exists, read it in
    if( open(FILE, $RJ_FILE_DATAFILE ) && ! defined($$cmd_ref{norun_job}) ){
        @lines = <FILE>;
        close( FILE );
    }
    # otherwise, default datafile
    else{
        $output_header = &rj_shell_header( SHELL=>$RJ_SHELL, 
                                           QUICK=>"", 
                                           MISC=>undef,
                                           SYSTEM=>undef );
        @lines = split(/\n/, $output_header);
        push( @lines, "&&&RJ_CMD_PRUN&&&\n" );
    }
    undef( %{$datafile_ref} );
    undef( @errors );
    @lines = grep( s/\s*$//, @lines );
    @{$$datafile_ref{orig}} = @lines;
    undef( @error_lines );
    # default values
    $line_num = sprintf( "%4d", 0 );
    $op = "=";
    $cond = "";

    &rj_datafile_var( "BATCH",       $cond, $op, "yes",                \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "BATCH_ARGS",  $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "DEBUGGER",    $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    # default value for CHAIN
    # We are having a problem with home directories not being mounted.
    # If no home directory, set default to "yes" since problem will die
    # almost immediately and let child job try to run.
    $default = "no";
    &my_stat( glob("~"), \%stat );
    if( ! %stat ){ # to test, check for existance of temp file
        push( @errors, "  Home directory missing" );
        $default = "yes"; # set to no for testing
    }
    # default ALL variables here -
    # will set to default value or "" (but will be defined)
    # will check command line args for variables
    &rj_datafile_var( "CHAIN",       $cond, $op, $default,             \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "EXEC",        $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "EXEC_ARGS",   $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "EXN",         $cond, $op, "0",                  \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "GROUP",       $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "NUMPE",       $cond, $op, $DEFAULT_NUMPE,       \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "OPT",         $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    $pname = &get_pname(".");
    &rj_datafile_var( "PNAME",       $cond, $op, $pname,               \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "PPN",         $cond, $op, $DEFAULT_PPN,         \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "PRUN_ARGS",   $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "SERIAL",      $cond, $op, "no",                 \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "TIME",        $cond, $op, $DEFAULT_TIME,        \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "TIME_B",      $cond, $op, $DEFAULT_TIME_B,      \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "TPP",         $cond, $op, $DEFAULT_TPP,         \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "TRY_NUM",     $cond, $op, $DEFAULT_TRY_NUM,     \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "TRY_NUM_MAX", $cond, $op, $DEFAULT_TRY_NUM_MAX, \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "TRY_TIME",    $cond, $op, $DEFAULT_TRY_TIME,    \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    &rj_datafile_var( "UMASK",       $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
    # go through command line args and init all user-defined vars
    # needed in case not set in $RJ_FILE_DATAFILE file
    foreach $var ( @{$$cmd_ref{var}{vars}} ) {
        if( $var =~ /^U_/ ){
            $op = "=";
            &rj_datafile_var( "$var",       $cond, $op, "",                   \@error_lines, "default", $line_num, $datafile_ref, $cmd_ref );
        }
    }
    # process each line
    $i = 0;
    while( $i <= $#lines ){
        $line_num = sprintf( "%4d", $i+1 );
        $line = $lines[$i];
        # embedded &&&RJ_<CMD>&&& that is not a comment or special #RJ line
        if( $line !~ /^\s*\#/ || $line =~ /^\s*\#RJ / ){
            $rest = $line;
            $line_new = "";
            $indent = "";
            while( $rest =~ /^(.*?)&&&(\S+)&&&(.*)$/ ){
                $pre  = $1;
                $rj_cmd  = $2;
                $rest = $3;
                $indent .= " "x length($pre);
                $replace = &rj_get_cmd(\$rj_cmd, $datafile_ref, INDENT=>$indent);
                if( !defined( $replace ) ){
                    push( @error_lines,
                          " $line_num: >> Invalid RJ_<OPTION> [$rj_cmd]",
                          " $line_num: $line" );
                    $replace = "&&&$rj_cmd&&&";
                }
                $line_new .= "$pre$replace";
            }
            $line_new .= $rest;
            $line = $line_new;
        }
	# NOTE: This is something lap wanted. It is still a work in progress. 
	#       The intention is to cleanup the test directory unless the test fails.
 	if ($line =~ /^\s*cts_diff/) {
 	  $line .= " || setenv CTS_KEEP_FILES 1";
 	}

        # save new line
        push( @{$$datafile_ref{new}}, "$line\n" ); $i++;
        # #RJ <VAR> = <value>
        if( $line =~ /^\s*\#\s*RJ\s+(\w+)\s*({[^\}]*})?\s*([\+\:\?]?=)\s*(.*?)\s*$/ ){
            $key  = $1;
            $cond = $2;
            $op   = $3;
            $val  = $4;
            &rj_datafile_var( $key, $cond, $op, $val,
                              \@error_lines, $line, $line_num, $datafile_ref, $cmd_ref );
        }
        # probable typo checks
        elsif( $line =~ /^\s*\#RJ\s+(\w+)\s+(\w+)/ ){
            $key  = $1;
            $val  = $2;
            push( @error_lines,
                  " $line_num: >> missing operator (eg '=') for #RJ line",
                  " $line_num: $line",
                  " $line_num: #RJ $key = $val  ? ");
        }
        elsif( $line =~ /^\s*\#RJ\s+(\w+.*?)\s*(\S*=)/ ){
            $key  = $1;
            $op  = $2;
            push( @error_lines,
                  " $line_num: >> could not parse line - possible invalid operator (eg '=') for #RJ line",
                  " $line_num:   key = [$key]",
                  " $line_num:   op  = [$op]",
                  " $line_num: $line" );
        }
        elsif( $line =~ /^\s*\#\s*RJ_(\w+)/ ){
            push( @error_lines,
                  " $line_num: >> Format: #RJ <WORD> = <VALUE>",
                  " $line_num: $line" );
        }
    }
    # if detected a "retry", create short-circuit command file
    if( $#errors >= 0 ){
        undef($$datafile_ref{new});
        # this string is searched for in TestProblem.pm
        # to test, touch temp file
        push( @{$$datafile_ref{new}}, "echo 'RJ_NEW_CHILD: SYSFAIL Attempting to restart run'\n" );
        $ierr = 0;
        &print_error( "System error(s):", @errors, "Attempting to restart" );
    }
    # do umask stuff here once you know it
    if( $$datafile_ref{var}{UMASK} =~ /\S/ ){
        # use that as the default
    }
    # if group is set, allow group rx (overridden by explicit umask setting)
    elsif( $$datafile_ref{var}{GROUP} =~ /\S/ ){
        $$datafile_ref{var}{UMASK} = sprintf( "%lo", umask() & 0727 );
    }
    # otherwise set to current
    else{
        # do not set at all - so that env var will not be set and users can use that as flag
    }
    # if group is defined, warn if not default group
    if( $$datafile_ref{var}{GROUP} =~ /\S/ &&
        defined($ENV{GROUP}) && $ENV{GROUP} ne $$datafile_ref{var}{GROUP} ){
        $ierr = 0;
        # the following just confuses people...so just skip the warning
        #&print_error( "Current unix group [$ENV{GROUP}] ne datafile GROUP [$$datafile_ref{var}{GROUP}]",
        #              "Files will not be changed to GROUP until run_job_cleanup.pl is called.",
        #              "Contact a sys admin to change your default unix group.",
        #              $ierr );
    }
    # required fields
    if( ! defined($$datafile_ref{var}{PNAME}) || $$datafile_ref{var}{PNAME} !~ /\S/ ){
        $$datafile_ref{var}{PNAME} = $pname;
    }
    # print any errors
    if( @error_lines ){
        $ierr = 1;
        $line_num = $i+1;
        &print_error( "Error in data file [$RJ_FILE_DATAFILE]",
                      @error_lines,
                      $ierr );
        exit( $ierr );
    }
    # fill in other values
    # NODES
    $$datafile_ref{var}{NODES} = POSIX::ceil($$datafile_ref{var}{NUMPE_max}/$$datafile_ref{var}{"PPN_min"});
    # NODES_exn (extra nodes to allocate)
    $exn = $$datafile_ref{var}{EXN};
    if( $exn =~ /(\S+)\s*\%/ ){
        $exn = ($1/100.0) * $$datafile_ref{var}{NODES};
        $exn = POSIX::ceil( $exn );
    }
    $$datafile_ref{var}{EXN} = $exn;
    $$datafile_ref{var}{NODES_exn} = $$datafile_ref{var}{NODES} + $exn;

    # adjust command line args if needed
    # depend
    if( defined($$cmd_ref{depend}) ){
        $val_new = "";
        @vals = split(/\s*,\s*/, $$cmd_ref{depend} );
        foreach $val( @vals ){
            # all current dir jobs
            if( $val eq "." ){
                &rj_get_jobids_other( $cmd_ref, \@jobids_other );
                $val = join(',',@jobids_other);
            }
            # no current dir jobs
            elsif( $val =~ /^no$/i ){
                $val = "";
                $$cmd_ref{nodepend} = "true";
            }
            if( defined($val) && $val =~ /\S/ ){
                $val_new .= ",$val";
            }
        }
        $val_new =~ s/^,+//;
        $$cmd_ref{depend} = $val_new;
    }
}

################################################################################

sub rj_datafile_var{
    my(
       $key,
       $cond,
       $op,
       $val,
       $error_lines_ref,
       $line,
       $line_num,
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
        $cmd_op,
        $cmd_op_last,
        $cmd_val,
        $cmd_var,
        $cond_eval,
        $cond_set,
        $found,
        $i,
        $num,
        $opt,
        $opt_no_yes,
        %opt_other,
        $opt_val,
        %opts,
        $satisfied,
        %time,
        $type,
        $val_cur,
        $val_new,
        );

    # see if satisfied
    $satisfied = "true";
    if( defined( $cond ) && $cond =~ /\S/ ){
        $satisfied = "false";
        $cond =~ s/^\{//;
        $cond =~ s/\}$//;
        $cond_eval = $cond;
        # go through each condition set and replace that string with true (1)
        # if that condition is set, or false (0) if not.
        # First will use ";" for "true", then replace all other strings with
        # false (0), then replace ";" with true (1)
        # replace pure numbers with true (right now, ";")
        $cond_eval =~ s/\b[1-9]\d*/\;/g;
        # replace each condition set with true (;)
        foreach my $cond_set ( split( /\s*,\s*/, $$cmd_ref{conds_total} ) ){
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
            $satisfied = "true";
        }
    }

    # if satisfied, set val
    if( $satisfied eq "true" ){

        # update var_file (just the file value of the var)
        if( ! defined($$datafile_ref{var_file}{$key}) ){
            $$datafile_ref{var_file}{$key} = "";
        }
        if( $key eq "GROUP" && $$datafile_ref{var_file}{GROUP} =~ /\S/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] can only be defined once.",
                  " $line_num: $line" );
        }
        if( $op eq "=" || $op eq ":=" ){
            # last equal op seen
            $$datafile_ref{var_file_op}{$key} = $op;
            $$datafile_ref{var_file}{$key} = $val;
        }
        elsif( $op eq "+=" ){
            if( $$datafile_ref{var_file}{$key} !~ /\S/ ){
                $$datafile_ref{var_file}{$key}  = $val;
            }
            else{
                $$datafile_ref{var_file}{$key} .= " $val";
            }
        }
        else{
            push( @{$error_lines_ref},
                  " $line_num: >> Invalid op from file [$op] .",
                  " $line_num: $line" );
        }

        # now apply command line then environment vars
        $cmd_op_last = "";
        $$datafile_ref{var}{$key} = $$datafile_ref{var_file}{$key};
        foreach $type ( "var", "var_env" ) {

            # skip if nothing
            if( ! $$cmd_ref{$type} ){
                next;
            }

            $num = $#{$$cmd_ref{$type}{vars}};
        for( $i = 0; $i <= $num; $i++ ){
                $cmd_var = $$cmd_ref{$type}{vars}[$i];
            if( $cmd_var eq $key ){
                    $cmd_op  = $$cmd_ref{$type}{ops}[$i];
                    $cmd_val = $$cmd_ref{$type}{vals}[$i];
                if( $cmd_op eq "=" || $cmd_op eq ":=" ){
                    $cmd_op_last = $cmd_op;
                    if( ( $cmd_op eq "=" &&
                              $$datafile_ref{var_file_op}{$key} ne ":=" ) || # := overrides this
                        ( $cmd_op eq ":=" ) ) {
                        $$datafile_ref{var}{$key} = $cmd_val;
                    }
                }
                elsif( $cmd_op eq "+=" ){
                    if( $cmd_op_last eq ":=" ||
                            $$datafile_ref{var_file_op}{$key} ne ":=" ){ # := overrides this
                        if( $$datafile_ref{var}{$key} !~ /\S/ ){
                            $$datafile_ref{var}{$key}  = $cmd_val;
                        }
                        else{
                            $$datafile_ref{var}{$key} .= " $cmd_val";
                        }
                    }
                }
                else{
                    push( @{$error_lines_ref},
                              " $line_num: >> Invalid op [$cmd_op] from $type",
                              "   [var = command line] or [var_env = environment]",
                          " $line_num:  var = [$cmd_var]",
                          " $line_num:   op = [$cmd_op]",
                          " $line_num:  val = [$cmd_val]",
                          " $line_num: $line" );
                }
            }
        }
        } # now apply command line then environment vars

    } # if satisfied, set val

    # current value
    $val_cur = $$datafile_ref{var}{$key};
    $found = "false";
    if( $val_cur =~ /\\\s*$/ ){
        push( @{$error_lines_ref},
              " $line_num: >> Line continuation character not allowed for variable definition.",
              " $line_num: $line" );
    }
    if( $key eq "EXEC" ){
        $found = "true";
    }
    if( $key =~ /^(PPN)$/ ){
        $found = "true";
        if( $val_cur =~ /^(\d+)\/(\d+)$/ ){
            $val_cur = POSIX::ceil(($DEFAULT_PPN*$1)/$2);
            $$datafile_ref{var}{$key} = $val_cur;
        }
        # keep track of the min value for the PPN to request
        if( !defined($$datafile_ref{var}{"${key}_min"}) ||
            $val_cur <= $$datafile_ref{var}{"${key}_min"} ){
            if( $satisfied eq "true" ){
                $$datafile_ref{var}{"${key}_min"} = $val_cur;
            }
        }

        # the first value is the default - need to keep track of
        # any value the user sets
        if( defined($$datafile_ref{var}{"${key}_min_set"}) ){
            if( $$datafile_ref{var}{"${key}_min_set"} eq "" ||
                $val_cur < $$datafile_ref{var}{"${key}_min_set"} ){
                $$datafile_ref{var}{"${key}_min_set"} = $val_cur;
            }
        }
        else{
            $$datafile_ref{var}{"${key}_min_set"} = "";
        }

    }

    if( $key =~ /^(NUMPE)$/ ){
        $found = "true";
        if( $val_cur =~ /^(\d+)\/(\d+)$/ ){
            $val_cur = POSIX::ceil(($DEFAULT_PPN*$1)/$2);
            $$datafile_ref{var}{$key} = $val_cur;
        }
        if( $val_cur !~ /^\d+$/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must have an integer value [$val_cur]",
                  " $line_num: $line" );
        }
        # keep track of the max value for the NUMPE to request
        if( !defined($$datafile_ref{var}{"${key}_max"}) ||
            $val_cur >= $$datafile_ref{var}{"${key}_max"} ){
            if( $satisfied eq "true" ){
                $$datafile_ref{var}{"${key}_max"} = $val_cur;
            }
        }
    }

    if( $key =~ /^(TPP|TRY_NUM|TRY_NUM_MAX)$/ ){
        $found = "true";
        if( $val_cur !~ /^\d+$/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must have an integer value [$val_cur]",
                  " $line_num: $line" );
        }
    }
    if( $key =~ /^(CHAIN|SERIAL)$/ ){
        $found = "true";
        if( $val_cur ne "no" && $val_cur ne "yes" ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must have value of 'yes' or 'no' [$val_cur]",
                  " $line_num: $line" );
        }
    }
    if( $key =~ /^(BATCH)$/ ){
        $found = "true";
        if( $val_cur !~ /^(yes|no|never|auto)$/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must have value of [yes|no|never|auto] [$val_cur]",
                  " $line_num: $line" );
        }

        # keep track if the user actually sets this
        if( $satisfied eq "true" && $line_num >0 ){
            # set val to what the user actually has in the datafile
            $$datafile_ref{var}{"${key}_set"} = $val;
        }
    }

    if( $key =~ /^(EXN)$/ ){
        $found = "true";
        if( $val_cur !~ /^((\d+)|((\d+(\.\d+)?\s*\%))|(\.\d+)\%)$/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must be N or N.M% [$val_cur]",
                  " $line_num: $line" );
        }
    }
    if( $key =~ /^(GROUP|PNAME)$/ ){
        $found = "true";
        if( $val_cur =~ /\s/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must not have whitespace [$val_cur]",
                  " $line_num: $line" );
        }
    }
    if( $key =~ /^(UMASK)$/ ){
        $found = "true";
        if( $val_cur =~ /\S/ && $val_cur !~ /^[0-7]+$/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] must be octal umask values [$val_cur]",
                  " $line_num: $line" );
        }
    }

    if( $key =~ /^(OPT)$/ ){
        $found = "true";
        # process (no|yes)_(opt) to keep latest
        $opt_other{"no"} = "yes";
        $opt_other{"yes"} = "no";
        undef( %opts );
        foreach $opt ( split(/\s+/, $val_cur) ){
            if( $opt =~ /^(no|yes)_(\S+)/ ){
                $opt_no_yes = $1;
                $opt_val = $2;
                $opts{$opt} = "";
                delete($opts{"$opt_other{$opt_no_yes}_${opt_val}"});
            }
        }
        $val_cur = join(" ", sort keys %opts);
        if( $val_cur =~ /\S/ && $val_cur !~ /^(((no|yes)_(fix_env|mapping))\s*)+$/ ){
            push( @{$error_lines_ref},
                  " $line_num: >> VAR [$key] valid options are: (no|yes)_(fix_env|mapping) [$val_cur]",
                  " $line_num: $line" );
        }
        $$datafile_ref{var}{$key} = $val_cur;
    }

    # debugger on at any time implies interactive
    if( $key eq "DEBUGGER" ){
        $found = "true";
        if( $val_cur =~ /\S/ ){
            $$cmd_ref{i} = "true";
        }
    }
    if( $key =~ /^(BATCH_ARGS|EXEC_ARGS|PRUN_ARGS|U_\S+)$/ ){
        # valid argument
        $found = "true";
    }
    if( $key =~ /^(TIME|TIME_B|TRY_TIME)$/ ){
        $found = "true";
        $val_new = $val_cur;
        %time = &conv_time( STRING=>$val_new );
        if( $satisfied eq "true" ){
            if( $key eq "TIME" ){
                $$datafile_ref{var}{$key} = $time{hms};
                $$datafile_ref{var}{TIME_SECS} = $time{SECS_TOTAL};
            }
            else{
                $$datafile_ref{var}{$key} = $time{SECS_TOTAL};
            }
        }
    }
    if( $found eq "false" ){
        push( @{$error_lines_ref},
              " $line_num: >> Invalid '#RJ <VAR>' line in data file [VAR=$key]",
              " $line_num: $line" );
    }
}

# print out key things from the datafile
sub rj_datafile_print{
    my (
        $datafile_ref,
        ) = @_;
    my(
        $out,
        $field,
        $maxlen,
        $newlen,
        );

    $out = "";
    $maxlen = 0;
    foreach $field ( keys %{$$datafile_ref{var}} ){
        $newlen = length( $field );
        if( $newlen > $maxlen ){
            $maxlen = $newlen;
        }
    }
    foreach $field ( sort keys %{$$datafile_ref{var}} ){
        $out .= sprintf( "%${maxlen}s = %s\n", $field, $$datafile_ref{var}{$field});
    }
    return( $out );
}

sub rj_get_cmd{
    my(
       $rj_cmd_ref
       ) = shift(@_);
    my(
       $datafile_ref
       ) = shift(@_);
    my %args = (
        INDENT => '',
        @_,
        );
    my $args_valid = ("INDENT");
    my(
       $arg,
       $ierr,
       $indent,
       $ret,
       $type,
       $var,
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
    $indent = $args{INDENT};
    # RJ_CMD_PRUN
    if( $$rj_cmd_ref eq "RJ_CMD_PRUN" ){
        # set type and record if ever a parallel job was done
        # used for setting up parallel jobs in case needed
        if( $$datafile_ref{var}{SERIAL} eq "yes" ){
            $type = "serial";
        }
        else{
            $type = "mpi";
            $$datafile_ref{var}{parallel} = $type;
        }
        # set run command
        if( $$datafile_ref{var}{EXEC} =~ /\S/ ){
            $ret =  &rj_get_run_command( EXEC=>"$$datafile_ref{var}{EXEC}",
                                         DEBUGGER=>"$$datafile_ref{var}{DEBUGGER}",
                                         NUMPE=>"$$datafile_ref{var}{NUMPE}",
                                         PPN=>"$$datafile_ref{var}{PPN}",
                                         TPP=>"$$datafile_ref{var}{TPP}",
                                         TYPE=>"$type",
                                         OPT=>$$datafile_ref{var}{OPT},
                                         INDENT=>$indent,
                                         EXEC_ARGS=>"$$datafile_ref{var}{EXEC_ARGS}",
                                         SHELL=>$RJ_SHELL,
                                         PRUN_ARGS=>"$$datafile_ref{var}{PRUN_ARGS}",
                                         ERROR=>\$ierr,
                                         );
        }
        else{
            $ierr = 1;
            &print_error( "EXEC not defined",
                          $ierr );
        }
    }
    # RJ_VAR_
    elsif( $$rj_cmd_ref =~ /RJ_VAR_(\S+)/ ){
        $var = $1;
        $ret = $$datafile_ref{var}{$var};
    }
    # RJ_(FILE|VAL)
    elsif( $$rj_cmd_ref =~ /(RJ_(FILE|VAL)_\S+)/ ){
        $var = $1;
        eval "\$ret = \$$var";
    }
    return( $ret );
}

sub rj_get_run_command{
    my %args = (
                ERROR            => undef, # error return
                EXEC             => undef, # executable
                EXEC_ARGS        => undef, # additional args to executable
                DEBUGGER         => undef, # if running under a debugger
                INDENT           => '',    # spaces to indent - readability for lining up
                NUMPE            => undef, # number of processes
                OUTPUT           => undef, # run command
                PATH_SET         => undef, # if the current PATH is set (defined)
                PPN              => undef, # processes per node
                TPP              => undef, # threads per process
                PRUN_ARGS        => undef, # additional args to prun
                OPT              => '',    # special options (like no_mapping)
                SHELL            => undef, # the shell being used
                TYPE             => undef, # serial, mpi
                @_,
               );
    my(
       $cmd,
       $cmd_add,
       $cmd_pre,
       $cmd_tmp,
       $debugger,
       $exec_args,
       $ierr,
       $indent,
       $key,
       @lines,
       $module_change,
       $mpirun,
       $num_threads,
       $nodes,
       $output,
       $output_def,
       $ppn,
       $prun,
       $prun_args,
       $prun_space,
       $prun_type,
       $prun_which,
       $ret,
       $run_command_ref,
       $run_exec,
       $run_job_mapping,
       $shell,
       %sys_info,
       $version,
       $wrapper,
       $wrapper_args_1,
       $wrapper_args_2,
       $wrapper_args_3,
       );
    #...invalid args
    foreach $key ( keys %args ) {
        if( $key !~ /^(ERROR|EXEC|EXEC_ARGS|DEBUGGER|INDENT|NUMPE|OUTPUT|PATH_SET|OPT|PPN|PRUN_ARGS|SHELL|TPP|TYPE)$/) {
            $ierr = 1;
            &print_error( "Invalid argument to get_run_command [$key]",
                          $ierr );
            ${$args{ERROR}} = $ierr;
            exit( $ierr );
        }
    }
    $shell = $args{SHELL};
    if( ! defined($shell) ){
        $shell = "sh";
    }
    $exec_args = "";
    if( defined($args{EXEC_ARGS}) ){
        $exec_args = " $args{EXEC_ARGS}";
    }
    $prun_args = "";
    if( defined($args{PRUN_ARGS}) ){
        $prun_args = " $args{PRUN_ARGS}";
    }
    # and others
    $debugger = $args{DEBUGGER};
    if( ! defined($debugger) ){
        $debugger = "";
    }
    $indent = $args{INDENT};
    if( $debugger !~ /^(|gdb|tv)$/ ){
        $ierr = 1;
        &print_error( "Invalid values for DEBUGGER [$debugger].",
                      "Valid values: gdb, tv, or <blank>", $ierr );
        ${$args{ERROR}} = $ierr;
        exit( $ierr );
    }

    # init vars
    ${$args{ERROR}} = 0;
    $run_command_ref = $args{OUTPUT};
    if( ! defined($run_command_ref) ){
        $run_command_ref = \$output_def;
    }
    $$run_command_ref = "";
    $cmd_pre = "";
    $cmd_add = "";
    $module_change = "no";
    $prun = "";
    $run_exec = "";
    $wrapper = "";
    $wrapper_args_1 = "";
    $wrapper_args_2 = "";
    $wrapper_args_3 = "";
    $prun_space = "";

    # set $$run_command_ref
    # serial
    if( $args{TYPE} eq "serial" ){
        $run_exec = $args{EXEC};
        # debugger
        if( $debugger eq "tv" ){
            $wrapper        = "totalview ";
            $wrapper_args_3 = " -a";
        }
        elsif( $debugger eq "gdb" ){
            $wrapper        = "gdb --args ";
        }
    }

    # mpi
    elsif( $args{TYPE} eq "mpi" ){
        $run_exec = $args{EXEC};

        # must set these
        if( ! defined( $args{NUMPE} ) && ! defined($args{PPN}) ){
            ${$args{ERROR}} = 1;
            &print_error( "Must define NUMPE (or at least PPN) if defining TYPE=$args{TYPE}",
                          1 );
            return( "" );
        }

        # get info about the system
        &get_sysinfo( \%sys_info );
        $ppn = $args{PPN} || $sys_info{L_PPN} || "undefined";
        # some mpis give error if ppn>numpe
        if( defined( $args{NUMPE} ) && $ppn ne "undefined" && $ppn > $args{NUMPE} ){
            $ppn = $args{NUMPE}
        }

        # determine run type
        # default
        $prun_type = "openmpi";

        # REDSTORM -> yod
        # DAWN     -> dawn
        # SEQUOIA  -> srun
        # cielo    -> aprun

        # first one wins
        $output = `which mpirun aprun yod poe srun 2> /dev/null`;
        @lines = split( /\n/, $output );
        $prun_which = $lines[0];
        grep( s/^.*\///, @lines );
        grep( s/\s//g, @lines );
        if( @lines ){
            $prun_type = $lines[0];
        }

        # found mpirun - see if mpich or not
        if( $prun_type eq "mpirun" ){
            $output = `($prun_type -v; $prun_type -V) 2>&1`;
            if( $output =~ /mpich/i ){
                $prun_type = "mpich";
            }
            else{
            $prun_type = "openmpi";
        }
        }

        # CIELO does not have aprun exec on front ends...
        if( $sys_info{L_CLASS} eq "CIELO" ||
            $sys_info{L_CLASS} eq "C_TT" ){
            $prun_type = "aprun";
        }

        # DAWN needs special things
        elsif( $sys_info{L_CLASS} eq "DAWN" ){
            $prun_type = "dawn";
        }

        # Even if srun exists and mpirun does not, use mpirun on some systems
        elsif( $prun_type eq "srun" && $sys_info{L_CLASS} =~ /TLCC/){
            $prun_type = "openmpi";
        }

        # always set OMP_NUM_THREADS
        # several packages using openmp and there is a problem with
        #   openmpi/1.6.3 and "-mca mpi_paffinity_alone 0" if you do
        #   not set it.
        # Does not seem to hurt anything if set
#         if( $args{OPT} !~ /\bfix_env\b/ ){
# 	  $args{OPT} .= " no_fix_env";
# 	}
        # Not sure if this hurts mpich...
        if( $args{OPT} !~ /\bno_fix_env\b/ ){
        $num_threads = 1;
        if( defined( $args{TPP} ) && $args{TPP} > 1 ){
            $num_threads = $args{TPP};
        }
        $cmd_pre .= ${indent}.&rj_shell_var_set( "OMP_NUM_THREADS", $num_threads, "; \\\n", $shell );
        }

        # openmpi
        # In general, grabbing whole node and letting mpi use parts of it.
        # Had issues:
        #   hurricane: partial node requested via msub but use whole node
        #              would cause memory issues.
        #   luna: round robin process-node scheduling if use partial node
        # see other comment: search for "whole node"
        if( $prun_type eq "openmpi" ){
            $mpirun = "mpirun";
            # fix modules - rj_fix_env
            if( $args{OPT} !~ /\bno_fix_env\b/ ){
            $cmd_pre .= ${indent}."eval `$0 --fix_env $args{EXEC} $shell`; \\\n";
            # other args gotten from rj_fix_env
            $prun_args .= " \${RJ_L_FIX_ENV_PRUN}";
            }

            # if ppn defined
            if( $ppn ne "undefined" ){
                # on DARWIN, --npernode must equal the number of processes or you can
                # just not have that argument at all
                if( $sys_info{L_CLASS} eq "DARWIN" ){
                }
                elsif( defined( $args{NUMPE} ) ){
                    # flags changed with 1.8
                    # In theory, should get whatever mpi is being used at time of run
                    # but cannot easily get that info.
                    if( $prun_which =~ m&openmpi([\-\./\d]+)& ){
                        $version = $1;
                    }
                    else{
                        $version = "unknown";
                    }
                    if( &my_compare_version($version,"<","1.8") ){
                    $prun_args .= " --npernode $ppn";
                }
                else{
                        $prun_args .= " --map-by ppr:$ppn:node";
                    }
                }
                else{
                    $prun_args .= " --pernode";
                }
            }
            # numpe
            if( defined($args{NUMPE}) ){
                $prun_args .= " -np $args{NUMPE}";
            }

            # on moonlight, bug when submitting jobs from front end
            if( $args{OPT} !~ /\bno_fix_env\b/ ){
	      if( $sys_info{L_OS} eq "ML" ){
		if( $shell eq "tcsh" || $shell eq "csh" ){
		  $cmd_pre .= "${indent}limit memorylocked unlimited && \\\n";
		}
		else{
		  $cmd_pre .= "${indent}ulimit -l unlimited && \\\n";
		}
	      }
	    }
	    
            # TPP
            # OMP_NUM_THREADS set above
            if( defined( $args{TPP} ) && $args{TPP} > 1 ){
                $prun_args .= " --cpus-per-proc $args{TPP} --bind-to-core";
            }
            # do a module list
            if( $module_change eq "yes" ){
                $cmd_pre .= "${indent}module list ; \\\n";
            }

            # debugger
            if( $debugger eq "tv" ){
                $wrapper_args_1 = " -tv";
            }
            elsif( $debugger eq "gdb" ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }

            # remove any extra ";"
            $cmd_pre =~ s/\;\s*\;/\;/g;
            $prun = $mpirun;
        }

        # aprun
        elsif( $prun_type eq "aprun" ){
            $prun = "aprun";
            $prun_args .= " -n $args{NUMPE}";
            if( $ppn ne "undefined" ){
                $prun_args .= " -N $ppn";
            }
            # TPP
            # OMP_NUM_THREADS set above
            if( defined( $args{TPP} ) && $args{TPP} > 1 ){
                $prun_args .= " -d $args{TPP}";
            }

            # env var to try to get around memory issues
            if( $args{OPT} !~ /\bno_fix_env\b/ ){
                $cmd_pre .= ${indent}.&rj_shell_var_set( "MPICH_GNI_MALLOC_FALLBACK",
                                                         "enabled", "; \\\n", $shell );
            }

            # debugger
            if( $debugger eq "tv" ){
                $wrapper = "totalview ";
                $wrapper_args_1 = " -a";
            }
            elsif( $debugger eq "gdb" ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }
        }

        # mpich
        # not sure if any of this works...do no have access to machine
        # with working mpich
        elsif( $prun_type eq "mpich" ){
            $prun = "mpirun";
            $prun_args .= " -np $args{NUMPE}";
            if( $ppn ne "undefined" ){
                $prun_args .= " -ppn $ppn";
            }
            # TPP
            # OMP_NUM_THREADS set above
            if( defined( $args{TPP} ) && $args{TPP} > 1 ){
                # do not know what to set
            }
            # debugger
            if( $debugger eq "tv" ){
                $wrapper = "totalview ";
                $wrapper_args_1 = " -a";
            }
            elsif( $debugger eq "gdb " ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }
        }

        # yod
        elsif( $prun_type eq "yod" ){
            $prun = "yod";
            $prun_args .= " -sz $args{NUMPE}";
            # TPP
            # OMP_NUM_THREADS set above
            if( defined( $args{TPP} ) && $args{TPP} > 1 ){
                # do not know what to set
            }
            # debugger
            if( $debugger eq "tv" ){
                $wrapper = "totalview ";
                $wrapper_args_1 = " -a";
            }
            elsif( $debugger eq "gdb " ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }
        }

        # poe
        elsif( $prun_type eq "poe" ){
            $prun = "poe";
            $prun_args .= " -sz $args{NUMPE}";
            # TPP
            # OMP_NUM_THREADS set above
            if( defined( $args{TPP} ) && $args{TPP} > 1 ){
                # do not know what to set
            }
            # debugger
            if( $debugger eq "tv" ){
                $wrapper = "totalview ";
                $wrapper_args_1 = " -a";
            }
            elsif( $debugger eq "gdb" ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }
        }

        # dawn
        elsif( $prun_type eq "dawn" ){
            $prun = "mpirun";
            $prun_args .= " -np $args{NUMPE}";
            if( $ppn == 1 ){
                $prun_args .= " -mode smp";
            }
            elsif( $ppn == 2 ){
                $prun_args .= " -mode dual";
            }
            elsif( $ppn == 4 ){
                $prun_args .= " -mode vn";
            }
            else{
                $ierr = 1;
                &print_error( "PPN invalid [$ppn]",
                              $ierr );
                ${$args{ERROR}} = $ierr;
                return;
            }
            # TPP
            # OMP_NUM_THREADS set above
            if( defined( $args{TPP} ) && $args{TPP} > 1 ){
                # do not know what to set
            }
            # debugger
            if( $debugger eq "tv" ){
                $wrapper = "totalview ";
                $wrapper_args_1 = " -a";
            }
            elsif( $debugger eq "gdb" ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }
        }

        # srun
        elsif( $prun_type eq "srun" ){
            $prun = "srun";
            $prun_args .= " -n $args{NUMPE}";
            if( $ppn ne "undefined" ){
               $prun_args .= " --ntasks-per-node $ppn";
               # need to specify number of nodes or will distribute evenly
               # across all nodes and ignore --ntasks-per-node
               $nodes = POSIX::ceil($args{NUMPE}/$ppn);
               $prun_args .= " --nodes=$nodes";
            }
            # debugger
            if( $debugger eq "tv" ){
                $wrapper = "totalview ";
                $wrapper_args_1 = " -a";
            }
            elsif( $debugger eq "gdb" ){
                $wrapper_args_2 = " xterm -e gdb --args";
            }
        }

        # will just use exec
        else{
            $ierr = 1;
            &print_error( "Couldn not find flavor of mpi (openmpi, srun, prun, poe, ... ).",
                          $ierr );
            ${$args{ERROR}} = $ierr;
            exit( $ierr );
        }
    }

    # create rj_mapping.txt file if needed
    # For now, just done for SEQIOIA and srun
    if( $args{OPT} !~ /\bno_mapping\b/ ){
        if( $args{OPT} =~ /\byes_mapping\b/ ||
            ( $cmd{sys_info}{L_MACHINE} eq "SEQUOIA" && $prun eq "srun" ) ){
            # if currently parallel or set up a parallel command at one time
            $run_job_mapping = &which_exec("$RJ_EXEC_MAPPING", QUIET=>'');
            if( $run_job_mapping =~ /\S/ ){
                $ret .=  &rj_get_run_command( EXEC=>$run_job_mapping,
                                              EXEC_ARGS=>'$_rj_filename',
                                              NUMPE=>$args{NUMPE},
                                              PPN=>$ppn,
                                              TYPE=>$args{TYPE},
                                              SHELL=>$args{SHELL},
                                              OPT=>"$args{OPT} no_mapping",
                                              INDENT=>${indent},
                                              ERROR=>\$ierr,
                    );
                # translate RJ_FILE_MAPPING_FULL into a real mapping file
                $cmd_tmp = &rj_shell_var_set( "_rj_filename", "`mktemp ${RJ_FILE_MAPPING_FULL}.XXXX`", "", $args{SHELL}, '"' );
                $cmd_add .= <<"EOF";
${indent}$cmd_tmp ; \\
${indent}echo "RJ_OUTPUT: Creating mapping file [\$_rj_filename]" ; \\
${indent}$ret ; \\
${indent}awk '/Task/ { print \$6, \$8 } ' \$_rj_filename | \\
${indent}  sed 's/[()]//g'    | \\
${indent}  tr ',' ' '         | \\
${indent}  sed 's/seqio//'    | \\
${indent}  sed 's/vulcanio//' | \\
${indent}  sed 's/-ib0//'     | \\
${indent}  gawk ' { printf( \"\%05d \%02d \%02d \%02d \%02d \%02d \%02d\\n\", \$7, \$1, \$2, \$3, \$4, \$5, \$6 );  } ' | \\
${indent}  sort               | \\
${indent}  awk ' { print \$2, \$3, \$4, \$5, \$6, \$7 } ' > \$_rj_filename.in ; \\
${indent}chmod 644 \$_rj_filename* ; \\
EOF
;
                # then add mapping file flag to prun_args
                # only valid for srun - but checked above for this (or yes_mapping)
                $prun_args .= ' --runjob-opts="--mapping $_rj_filename.in"';
            }
            else{
                $ierr = 0;
                &print_error( "Could not find executable that creates mapping file [$RJ_EXEC_MAPPING]",
                              "Setting opt=no_mapping" );
            }
        }
    }

    # string together vars to create run_command_ref
    if( ${prun} =~ /\S/ ){
        $prun_space = " ";
    }
    $$run_command_ref  = "${cmd_pre}${indent}${wrapper}${prun}${wrapper_args_1}${prun_args}${wrapper_args_2}${prun_space}${run_exec}${wrapper_args_3}${exec_args}";
    $$run_command_ref =~ s/^\s+//;
    # do additional commands
    if( $cmd_add =~ /\S/ ){
        $$run_command_ref = "$cmd_add$$run_command_ref";
    }
    $$run_command_ref =~ s/^\s+//;

    return( $$run_command_ref );
}

# module commands
# would normally just be "module $command ; " ... but there is an issue in setup
# in csh that does not correctly handle aliases.  So, even though sourcing
# the correct init file (rj_shell_module_setup), the module alias is not
# being set correctly, so do the full eval command
sub rj_shell_module_command{
    my(
        $command,
        $shell,
        ) = @_;
    my(
        $ret,
        );
    $ret = "";

    if( $shell eq "bash" || $shell eq "sh" ){
        $ret = " module $command ; ";
    }
    else{
        $ret = " eval `/usr/bin/modulecmd tcsh $command` ; ";
    }

    return( $ret );

}

# module setup line
# see comment above about rj_shell_module_setup and rj_shell_module_command
sub rj_shell_module_setup{
    my(
        $shell,
        ) = @_;
    my(
        $dir,
        $suffix,
        $file,
        );
    if( $shell eq "bash" || $shell eq "sh" ){
        $suffix = ".sh";
    }
    else{
        $suffix = ".csh";
    }
    $file = "";
    foreach $dir ( "/opt/modules/default/etc", "/etc/profile.d" ) {
        $file = "$dir/modules$suffix";
        if( -e $file ){
            last;
        }
        else{
            $file = "";
        }
    }
    if( $file =~ /\S/ ){
        return( "source $file ; " );
    }
    else{
        return( "" );
    }
}

# setting (or unset) variables for different shells
sub rj_shell_var_set{
    my(
        $var,      # variable name
        $val,      # value  (if undef -> unset)
        $end,      # if end in ";", or ";\n", or ...
        $shell,    # shell to use
        $quote,    # default quote to use
        ) = @_;
    # set
    if( defined($val) ){
        # foo 'bar --> 'foo '"'"'bar'
        $val =~ s/'/'"'"'/g;
        if( defined( $quote ) ){
            $val = $quote.$val.$quote;
        }
        else{
            $val = "'$val'";
        }
        if( defined( $shell ) && ( $shell eq "bash" || $shell eq "sh" ) ){
            return( "export $var=$val $end" );
        }
        else{
            return( "setenv $var $val $end" );
        }
    }

    # unset
    else{
        if( defined( $shell ) && ( $shell eq "bash" || $shell eq "sh" ) ){
            return( "unset $var $end" );
        }
        else{
            return( "unsetenv $var $end" );
        }
    }
}

# if statement for different shells
sub rj_shell_if{
    my( 
        $cif,
        $condition,
        $command,
        $cendif,
        $indent,
        $shell,
        ) = @_;
    my(
        $ret,
        );

    $ret = "";

    # bash
    if( $shell eq "sh" || $shell eq "bash" ){
        if( $cif eq "if" ){
            $ret .= "${indent}if ";
        }
        elsif( $cif eq "elsif" ){
            $ret .= "${indent}elif ";
        }
        else{
            $ret .= "${indent}else";
        }
        if( $cif ne "else" ){
            $ret .= "[ $condition ] ; then";
        }
        $ret .= "\n";
        $ret .= "$command\n";
        if( $cendif eq "endif" ){
            $ret .= "${indent}fi\n";
        }
    }

    # tcsh - filling this out - 
    else{
        if( $cif eq "if" ){
            $ret .= "${indent}if ";
        }
        elsif( $cif eq "elsif" ){
            $ret .= "${indent}elsif ";
        }
        else{
            $ret .= "${indent}else";
        }
        if( $cif ne "else" ){
            $ret .= "( $condition ) then";
        }
        $ret .= "\n";
        $ret .= "$command\n";
        if( $cendif eq "endif" ){
            $ret .= "${indent}endif\n";
        }
    }
    return( $ret );
}

# header for different shells (init system files and modules)
sub rj_shell_header{
    my %args = (
        MISC    => undef,   # string between interpreter and other setup
        QUICK   => undef,   # if quick load
        SHELL   => "tcsh",  # shell to use
        SYSTEM  => undef,   # if load system/module files
        @_
        );
    
    my $args_valid = "MISC|QUICK|SHELL|SYSTEM";
    my(
        $arg,
        $ierr,
        $prompt,
        $ret,
        $system_string_bash,
        $system_string_tcsh,
        );
    foreach $arg ( keys %args ){
        if( $arg !~ /^($args_valid)$/ ){
            $ierr = 1;
            &print_error( "Invalid argument [$arg]",
                          "Valid args [$args_valid]",
                          $ierr );
            exit( $ierr );
        }
    }
    $ret = "";
    
    $system_string_bash .= <<'EOF';
# 
# NOTE: THIS BLOCK IS USED IN run_job.pl and in .bashrc 
# NOTE: TO SET UP SYSTEM/MODULE FILES
# 
# nothing for system files yet
# set up modules if needed
declare -f -F module > /dev/null
_l_status=$?
if [ "$_l_status" != "0" ] ; then
  unset _l_file
  for _l_dir in /opt/modules/default/etc /etc/profile.d ; do
   if [ -e "$_l_dir/modules.sh" ] ; then
     _l_file="$_l_dir/modules.sh"
   fi
   if [ -e "$_l_dir/modules.bash" ] ; then
     _l_file="$_l_dir/modules.bash"
   fi
   if [ ! -z ${_l_file+x}  ] ; then
     . $_l_file
     break
   fi
  done
fi
unset _l_file
unset _l_dir
# 
# NOTE: DONE: THIS BLOCK IS USED IN run_job.pl and in .bashrc 
# NOTE: DONE: TO SET UP SYSTEM/MODULE FILES
# 
EOF
;

    $system_string_tcsh = <<'EOF';
# 
# NOTE: THIS BLOCK IS USED IN run_job.pl and in .cshrc 
# NOTE: TO SET UP SYSTEM/MODULE FILES
# 
# Figure out if need to source any system files
#    /etc/csh.cshrc (all shells)
#    /etc/csh.login (login shells)
# Turns out some systems (cielo) put important stuff into csh.login that
# also needs to be sourced at least once when doing non-login shell
# if running interactive, consider them sourced and inherited shells
# will never source the system files.
if ( { test -t 0 } && { test -t 1 } && ! $?0 ) then
   setenv L_EAP_SD_SYS_CSHRC
   setenv L_EAP_SD_SYS_LOGIN
endif
# default is to source system files
setenv L_EAP_S_SYS_CSHRC
setenv L_EAP_S_SYS_LOGIN
# do not source them if more than basic path
if ( "${PATH}" != "/usr/bin:/bin" ) then
   unsetenv L_EAP_S_SYS_CSHRC
   unsetenv L_EAP_S_SYS_LOGIN
endif
# Cielo puts important stuff into system login unfortunately (PATH moab)
# So source the system login even if in non-login shell.
# Put this check at the end since needs to override above PATH check
# since system cshrc would have been sourced and PATH larger
# How to tell if in login shell:
#    csh: $?prompt
#   tcsh: $?loginsh
# Problem on mac since gmake Pack_build_foo w/out having sourced
#   team cshrc first will load system login and change path from what
#   we wanted.
# The problem is that the first time this .cshrc is sourced in a non-login
# sub-shell from a login shell, it will re-source the system login
# even though it does not need to.  So, the PATH will get altered.
# Sent email asking admins to fix cielo/cielito /etc/csh.cshrc
# so, for now, just add stuff to path at end (path cleaned up
# at end of .cshrc
# NOTE: Fixed as of 20130625
#if( ! $?L_EAP_SD_SYS_LOGIN && ! $?prompt ) then
#   # from above problems with system, cannot source the system login
#   # setenv L_EAP_S_SYS_LOGIN
#   setenv PATH ${PATH}:/opt/MOAB/default/bin:/opt/PBS/default/bin
#endif
# remember old values
set _l_umask=`umask`
if( $?prompt ) then
  set _l_prompt="$prompt"
else
  unset _l_prompt
endif
# if have not already sourced (or considered sourced) and need to source, do so
if( ! $?L_EAP_SD_SYS_CSHRC && $?L_EAP_S_SYS_CSHRC ) then
   if( -f /etc/csh.cshrc ) then
     source /etc/csh.cshrc
     # if actually did source it, set to 1
     setenv L_EAP_SD_SYS_CSHRC 1
   endif
endif
if( ! $?L_EAP_SD_SYS_LOGIN && $?L_EAP_S_SYS_LOGIN ) then
   if( -f /etc/csh.login ) then
     source /etc/csh.login
     # if actually did source it, set to 1
     setenv L_EAP_SD_SYS_LOGIN 1
   endif
endif
# set back to old values
umask $_l_umask
if( $?_l_prompt ) then
   set prompt="$_l_prompt"
else
   unset prompt
endif
unset _l_umask
unset _l_prompt
# set up modules if needed
unset _l_setup_modules
set _l_is_set=`alias module`
# if no alias, set it up
if( "$_l_is_set" == "" ) then
  set _l_setup_modules=1
endif
# if not set or equal to nothing, set it up
if( ! $?MODULEPATH ) then
  set _l_setup_modules=1
else
 if( "$MODULEPATH" == "" ) then
  # need to unset it or it will not be overwritten
  unsetenv MODULEPATH
  set _l_setup_modules=1
 endif
endif
if ( $?_l_setup_modules ) then
  foreach _l_loc ( /opt/modules/default/etc /etc/profile.d )
    if ( -e $_l_loc/modules.csh ) then
      source $_l_loc/modules.csh
      break
    endif
  end
  unset _l_loc
endif
unset _l_is_set
unset _l_setup_modules
# 
# NOTE: DONE: THIS BLOCK IS USED IN run_job.pl and in .cshrc 
# NOTE: DONE: TO SET UP SYSTEM/MODULE FILES
# 
EOF
;

    # bash
    if( defined( $args{SHELL} ) && ( $args{SHELL} eq "bash" || $args{SHELL} eq "sh" ) ){
        
        # interpreter
        $ret .= "#!/bin/bash";
        if( defined( $args{QUICK} ) ){
            $ret .= " -f";
        }
        $ret .= "\n";
        
        # MISC
        if( defined( $args{MISC} ) ){
            $ret .= $args{MISC};
        }
        
        # MODULES
        if( defined( $args{SYSTEM} ) ){
            $ret .= "\n";
            $ret .= "# ========================================================\n";
            $ret .= $system_string_bash;
            $ret .= "# ========================================================\n";
            $ret .= "\n";
        }
        
    }
    
    # tcsh
    else{

        # interpreter
        $ret .= "#!/bin/tcsh";
        if( defined( $args{QUICK} ) ){
            $ret .= " -f";
        }
        $ret .= "\n";

        # MISC
        if( defined( $args{MISC} ) ){
            $ret .= $args{MISC};
        }
        
        # SYSTEM
        if( defined( $args{SYSTEM}) ){
            $ret .= "\n";
            $ret .= "# ========================================================\n";
            $ret .= $system_string_tcsh;
            $ret .= "# ========================================================\n";
            $ret .= "\n";
        }
    }
    
    return( $ret );
}

#............................................................................
#...Name
#...====
#... rj_submit
#...
#...Purpose
#...=======
#... submits the next job
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line
#...              $cmd{$option} = value
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub rj_submit{
    my(
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
       $batchid,
       $command,
       $command_clean,
       $cwd,
       $dir,
       $env_name,
       $env_name_val,
       $env_val,
       $hack_needed,
       $ierr,
       $launchtype,
       $num,
       $op,
       $output,
       $output_misc,
       $quote,
       $val,
       $var,
       $wrapper_file,
       );

    $ierr = 0;

    # default value
    $$cmd_ref{submit_id} = -1;
    # return if this is not first one and not chaining
    # (always submit first one)
    if( defined($$cmd_ref{id}) && $$datafile_ref{var}{CHAIN} ne "yes" ){
        # do log file here if just returning
        &add_log( $datafile_ref, $cmd_ref );
        return;
    }
    # create wrapper script that will be used to run RJ_CMD_FILE .
    # The reason we create a wrapper script is that some prep work setting up the environment
    # needs to be done.  Also, some batch systems cannot handle a very long script.
    # get tag to be used for submitted job - must be different than tag of current job
    $$cmd_ref{tag_submit} = $$cmd_ref{tag};
    if( $$cmd_ref{tag_submit} eq $$cmd_ref{tag} ){
        $$cmd_ref{tag_submit} = $$cmd_ref{tag} + 1;
    }
    $cwd = getcwd();
    # fully qualify path to get around MOAB bug
    # and used by run_status.pl for directory locations
    $$datafile_ref{submit_file} = "$cwd/${RJ_FILE_BATCH_OUT}.$$cmd_ref{tag_submit}";
    # if the launch is to batch, note it
    if( $$datafile_ref{var}{BATCH} eq "yes" ||
        ( $$datafile_ref{var}{BATCH} eq "auto" &&
          $$cmd_ref{sys_info}{L_END} eq "FRONT" ) ){
        $launchtype = "batch";
    }
    else{
        $launchtype = "nobatch";
    }

    # Check for BATCH_set/BATCH=never and would run on back end
    if( ( defined( $$datafile_ref{var}{BATCH_set} ) &&
          $$datafile_ref{var}{BATCH_set} eq "never" ) ||
        $$datafile_ref{var}{BATCH} eq "never" ){
        # if launching (so get to back end) or not launching and already on the back end
        if( $launchtype eq "batch" || $$cmd_ref{sys_info}{L_END} eq "BACK" ){
            $ierr = 1;
            &print_error( "BATCH==never (current value or set explicitly in $RJ_FILE_DATAFILE).",
                          "However, would end up running job on back end.",
                          $ierr );
            exit( $ierr );
        }
    }

    $batchid = $$cmd_ref{sys_info}{batchid};
    if( $batchid ne "" ){
        $batchid = "--batchid $batchid";
    }
    # specify dir on command line so that run_status.pl can see it
    $dir = getcwd();

    # start of command to use to run run_job.pl
    $command = "$0 go --id $$:PID --launchtype $launchtype --tag $$cmd_ref{tag_submit} --dir $dir $batchid";
    # also pass along these command line args
    # conds
    if( defined( $$cmd_ref{conds} ) && $$cmd_ref{conds} =~ /\S/ ){
        $command .= " --conds $$cmd_ref{conds}";
    }
    # env
    if( defined( $$cmd_ref{env} ) ){
        $command .= " --env ";
        foreach $env_name_val ( @{$$cmd_ref{env}} ){
            ($env_name, $env_val) = split(/=/, $env_name_val);
            $command .= "$env_name";
            if( defined( $env_val ) ){
                $command .= "='$env_val'";
            }
            $command .= ",";
        }
        $command =~ s/,$//;
    }

    if( defined($$cmd_ref{var}{vars}) ){
        $num = $#{$$cmd_ref{var}{vars}};
        for( $i = 0; $i <= $num; $i++ ){
            $var = $$cmd_ref{var}{vars}[$i];
            $op  = $$cmd_ref{var}{ops}[$i];
            $val = $$cmd_ref{var}{vals}[$i];
            # try to create correct quotation - this is always a crapshoot
            # ticks seen
            if( $val =~ /\'/){
                $quote = '"';
            }
            else{
                $quote = "'";
            }
            $command .= " --var $var$op$quote$val$quote";
        }
    }
    # do not pass these through
    # wait: probably just want master process to wait
    #if( defined($$cmd_ref{wait}) ){
    #    $command .= " --wait";
    #}
    # v: probably users would just want verbose on master process
    #if( defined($$cmd_ref{v}) ){
    #    $command .= " -v";
    #}
    if( defined($$cmd_ref{shell}) ){
        $command .= " --shell $$cmd_ref{shell}";
    }

    # depend
    # if gotten to here, then already gotten through required depend,
    # but still need to skip depend check if that was originally
    # requested.
    # Also, if depend set, then again already gotten through required
    # depend and no longer need to check existing jobs
    if( defined($$cmd_ref{nodepend}) ||
        defined($$cmd_ref{depend}) ){
        $command .= " --nodepend";
    }

    # name of wrapper file to run
    $wrapper_file = "$RJ_FILE_WRAPPER.$$cmd_ref{tag_submit}";

    # start to set output
    $output_misc  = "echo ===============================================\n";
    $output_misc .= "date +'TIME_S:\%s'\n";
    $output_misc .= "date +'TIME_YMDHMS:\%Y.\%m.\%d.\%H.\%M.\%S'\n";
    $output_misc .= "echo ===============================================\n";

    # hack to fix environment setting for environment variables
    # for cielo/cielito flavor of moab since it does not pass env vars set to ""
    # correctly
    if( $$cmd_ref{sys_info}{L_MACHINE} eq "CIELITO" ||
        $$cmd_ref{sys_info}{L_MACHINE} eq "CIELO" ||
        $$cmd_ref{sys_info}{L_MACHINE} eq "M_TT" ){
        foreach $env_name ( keys %ENV ){
            if( $ENV{$env_name} eq "" || $ENV{$env_name} eq "$RJ_ENV_EMPTY" ){
                $hack_needed = "";
                last;
            }
        }
    }
    if( defined($hack_needed) ){
        $output_misc .= "# ===============================================\n";
        $output_misc .= "# moab CIELO/CIELITO hack to reset empty environment variables\n";
        foreach $env_name ( keys %ENV ){
            if( $ENV{$env_name} eq "" || $ENV{$env_name} eq "$RJ_ENV_EMPTY" ){
                $output_misc .= "setenv $env_name\n";
            }
        }
        $output_misc .= "# ===============================================\n";
    }

    $output  = &rj_shell_header(SHELL=>"tcsh", SYSTEM=>"", QUICK=>"", MISC=>$output_misc);
    $output .= "echo ===============================================\n";
    $output .= "echo 'Command: [process info]'\n";
    $output .= "echo PID: \$\$\n";
    $output .= "ps -l\n";
    $output .= "echo ===============================================\n";
    $output .= "echo 'Command: [module list]'\n";
    $output .= "module list\n";
    $output .= "echo ===============================================\n";
    $output .= "echo 'Command: [printenv | sort]'\n";
    $output .= "printenv | sort\n";
    # takes too much time on some machines - comment out
    #$output .= "echo ===============================================\n";
    #$output .= "echo 'Command: [run_status.pl -u all]'\n";
    #$output .= "run_status.pl -u all\n";
    $output .= "echo ===============================================\n";
    $output .= "echo 'Command: [check_nodes.pl]'\n";
    $output .= "check_nodes.pl\n";
    $output .= "echo ===============================================\n";
    $output .= "echo 'Command: [cat /etc/motd]'\n";
    $output .= "cat /etc/motd\n";
    $output .= "echo ===============================================\n";
    ( $command_clean = $command ) =~ s/'/'"'"'/g;
    $output .= "echo 'Command: [$command_clean]'\n";
    $output .= "$command\n";
    $output .= "echo ===============================================\n";
    $output .= "date +'TIME_S:\%s'\n";
    $output .= "date +'TIME_YMDHMS:\%Y.\%m.\%d.\%H.\%M.\%S'\n";
    $output .= "echo ===============================================\n";
    if( ! open( FILE, ">$wrapper_file" ) ){
        $ierr = 1;
        &print_error( "Cannot open command run file [$wrapper_file]",
                      $ierr );
        exit( $ierr );
    }
    print FILE $output;
    close( FILE );
    chmod &my_mode( EXEC=>"", DEC=>"" ), $wrapper_file;
    # submit depending upon batch system
    if( $$datafile_ref{var}{BATCH} eq "yes" ||
        ( $$datafile_ref{var}{BATCH} eq "auto" &&
          $$cmd_ref{sys_info}{L_END} eq "FRONT" ) ){
        $ierr = &rj_submit_batch( $datafile_ref, $cmd_ref, $wrapper_file );
    }

    # otherwise, nobatch
    else{
        &rj_submit_nobatch( $datafile_ref, $cmd_ref, $wrapper_file );
    }

    if( defined($cmd{v}) ){
    &print_message( $cmd_ref, "Dependent", "$$cmd_ref{submit_id}" );
    &print_message( $cmd_ref, "Finish_Submit", "" );
    }
}

sub rj_check_other_jobs{
    my(
        $cmd_ref,
        ) = @_;
    my(
        $current_id,
        $ierr,
        @jobids_other,
        $out,
        );

    # stop if already job(s) in this directory
    &rj_get_jobids_other( $cmd_ref, \@jobids_other );

    if( $#jobids_other >= 0 ){
        $ierr = 1;
        $current_id = $$cmd_ref{id};
        if( ! defined($current_id) ){
            $current_id = "none";
        }
        $out = `$RUN_STATUS_PL -d .`;
        $out =~ s/^\s*//;
        $out =~ s/^\s*$//;
        $out =~ s/^/       \| /mg;
        &print_error( "Non-interactive job(s) already running in this directory.",
                      "Current id [$current_id]",
                      "Existing id(s) (run_status.pl, $RJ_FILE_ID_NEXT):",
                      @jobids_other,
                      "\n$out",
                      "To launch a job after the current job(s) is(are) finished:",
                      "  $0 --depend .",
                      "NOT SUBMITTING JOB",
                      "",
                      $ierr );
        exit( $ierr );
    }
}

# return array of other non-interactive jobs in the current directory
sub rj_get_jobids_other{
    my(
        $cmd_ref,
        $jobids_other_ref
        ) = @_;
    my(
        %dirs_jobid,
        $jobid,
        @jobids,
        @jobids_all,
        $myid,
        %names_jobid,
        $rundir,
        %run_status_info,
        %stat,
        );
    # ..........................................
    # check for jobs running in the batch system
    # ..........................................
    undef( %run_status_info );
    &get_run_status_info( \%run_status_info, \%dirs_jobid, \%names_jobid, "-d ." );
    &my_stat( ".", \%stat );
    $rundir = $stat{fullpath};
    if( defined($dirs_jobid{$rundir}) ){
        @jobids_all = sort keys %{$dirs_jobid{$rundir}};
        undef( @jobids );
        # remove interactive jobids
        foreach $jobid ( @jobids_all ){
            # skip interactive
            if( defined($run_status_info{$jobid}{vals}{NOTE}) &&
                $run_status_info{$jobid}{vals}{NOTE} =~ /\bInt\b/i ){
                next;
            }
            # skip STATE==Completed
            if( $run_status_info{$jobid}{vals}{STATE} eq "Completed" ){
                next;
            }
            # skip STATE=="-" (not running)
            if( $run_status_info{$jobid}{vals}{STATE} eq "-" ){
                next;
            }
            push( @jobids, $jobid );
        }
        # only valid existing jobid would be this one (parent currently running)
        undef( @$jobids_other_ref );
        if( ! defined($$cmd_ref{id}) ){
            @$jobids_other_ref = @jobids;
        }
        else{
            $myid = &my_getid( $cmd_ref );
            foreach $jobid ( @jobids ){

                # skip if the jobid matches the RJ_RUNNING_ID
                if( defined($ENV{RJ_RUNNING_ID}) &&
                    "$ENV{RJ_RUNNING_ID}:PID" eq $jobid ){
                    next;
                }

                # skip if this pid matches the one found
                if( "$$:PID" eq $jobid ){
                    next;
                }
                    
                # might be slightly different form - so only fail if
                # both not subsets of others
                if( index( $myid, $jobid ) < 0 && index($jobid, $myid) < 0 ){
                    push( @$jobids_other_ref, $jobid )
                }
            }
        }
    }
}

################################################################################

# script to submit the batch job
sub rj_submit_batch{
    my(
       $datafile_ref,
       $cmd_ref,
       $wrapper_file,
       ) = @_;
    my(
       $depend,
       $depend_str,
       $done,
       $env_name,
       $hack_one_depend,
       $id_parent,
       $ierr,
       @lines,
       $submit,
       $submit_try,
       $multiple,
       $nodes,
       $note,
       $opts,
       $output,
       $output1,
       %pids_depend,
       $ppn,
       $submit_command,
       $umask_val,
       );

    $ierr = 0;

    # init options to send to batch system
    $opts = "";

    # init
    $submit = "";

    # make sure batch submit exists
    if( $submit eq "" ){
        $submit_try = "msub";
        $submit = &which_exec($submit_try, QUIET=>"");
    }
    if( $submit eq "" ){
        $submit_try = "sbatch";
        $submit = &which_exec($submit_try, QUIET=>"");
    }

    if( $submit eq "" ){
        $ierr = 1;
        &print_error( "Cannot find batch submission command [msub or sbatch]",
                      "To run without batch:",
                      "  $RJ_FILE_DATAFILE: '#RJ BATCH=no'",
                      "  Command Line:  '--batch no'",
                      $ierr );
        exit( $ierr );
    }

    # will assume that the batch system being used is moab - with different flavors
    # time
    if( defined( $$datafile_ref{var}{TIME}) ){
        if( $submit_try eq "msub" ){
        $opts .= " -l walltime=$$datafile_ref{var}{TIME}";
    }
        else{
            $opts .= " -t $$datafile_ref{var}{TIME}";
        }
    }

    # name
    if( $submit_try eq "msub" ){
    $opts .= " -N $$datafile_ref{var}{PNAME}";
    }
    else{
        $opts .= " -J $$datafile_ref{var}{PNAME}";
    }

    # stderr, stdout
    if( $submit_try eq "msub" ){
    $opts .= " -j oe -o $$datafile_ref{submit_file}";
    }
    else{
        $opts .= " -e $$datafile_ref{submit_file} -o $$datafile_ref{submit_file}";
    }

    # run dir
    if( $submit_try eq "msub" ){
    $opts .= " -d ".getcwd();
    }
    else{
        $opts .= " -D ".getcwd();
    }

    # nodes
    # on dawn, need nodes in 128 unit chunks
    if( $$cmd_ref{sys_info}{L_CLASS} eq "DAWN" ){
        $multiple = 128;
        $nodes = ceil($$datafile_ref{var}{NODES_exn}/$multiple)*$multiple;
        $opts .= " -l nodes=$nodes";
    }
    else{
        # In general, grabbing whole node and letting mpi use parts of it.
        # Had issues:
        #   typhoon: mpi_paffinity_alone=0 and partial node horrible memory
        #            performance.  Run stress test with 5 second cpu and see.
        #   hurricane: partial node requested via msub but use whole node
        #              would cause memory issues.
        #   luna: round robin process-node scheduling if use partial node
        # see other comment: search for "whole node"
        #$ppn = $$datafile_ref{var}{PPN};
        $ppn = $$cmd_ref{sys_info}{L_PPN};
        if( $submit_try eq "msub" ){
        $opts .= " -l nodes=$$datafile_ref{var}{NODES_exn}:ppn=$ppn";
    }
        else{
            $opts .= " -N $$datafile_ref{var}{NODES_exn}";
        }
    }

    # depend
    $id_parent = &my_getid( $cmd_ref );
    # might need to chop off non-digits...not sure
    #if( defined($id_parent) ){
    #    $id_parent =~ s/^(\d+).*$/$1/;
    #}
    # wait for parent to finish
    if( defined( $$cmd_ref{id} ) ){
        if( $submit_try eq "msub" ){
        $opts .= " -l depend=$id_parent";
    }
        else{
            $opts .= " -d $id_parent";
        }
    }
    # and other dependencies
    if( defined($$cmd_ref{depend}) && $$cmd_ref{depend} =~ /\S/ ){
        undef( %pids_depend );
        # hack to get around moab breakage
        # can only depend on latest one.  If depend on all the jobs
        # that it should depend on, get stuck in "Hold=Batch"
        $depend_str = "";
        foreach $depend ( reverse split( /,/, $$cmd_ref{depend} ) ){
            # tack onto batch dependency args if a batch job
            if( $depend !~ /^(\d+):PID$/ ){
                if( $submit_try eq "msub" ){
                if( ! defined($hack_one_depend) ){
                    $hack_one_depend = "ugh-broken-moab-installed-so-hack";
                        $depend_str .= " -l depend=$depend";
                    }
                }
                else{
                    $depend_str .= "$depend:";
                }
            }
            # if PID, tack onto other list
            else{
                $pids_depend{$depend} = "";
            }
        }
        if( $depend_str ne "" ){
            if( $submit_try eq "msub" ){
                $opts .= " $depend_str";
            }
            else{
                $depend_str =~ s/:$//;
                $opts .= " -d afterany:$depend_str";
            }
        }
        # if waiting on PIDs, need to wait before submitting
        if( keys %pids_depend ){
            $output = join( ", ", keys %pids_depend );
            &print_message( $cmd_ref, "Waiting:Batch", $output );
            &my_wait( $cmd_ref, IDS=>\%pids_depend );
        }
    }

    # preserve umask
    $umask_val = umask();
    # turn back into octal
    $umask_val = sprintf( "0%o", $umask_val );
    if( $umask_val =~ /^\d+$/ ){
        if( $submit_try eq "msub" ){
        $opts .= " -W umask=$umask_val";
    }
        else{
            # slurm does not seem to have this
        }
    }

    # additional opts
    if( $$cmd_ref{sys_info}{L_CLASS} eq "CIELO" ){
    }
    elsif( $$cmd_ref{sys_info}{L_CLASS} eq "REDSTORM" ){
    }
    elsif( $$cmd_ref{sys_info}{L_CLASS} eq "DAWN" ){
        $opts .= " -q pdebug";
    }
    else{
    }

    # inherit environment
    if( $submit_try eq "msub" ){
        $opts .= " -V";
    }
    else{
        # slurm exports all
    }

    $opts .= " $$datafile_ref{var}{BATCH_ARGS}";

    $submit_command = "$submit $opts $wrapper_file";
    if( defined($$cmd_ref{debug}) ){
        $note = "[DEBUG] ";
    }
    else{
        $note = "";
    }
    &print_message( $cmd_ref, "Command", "${note}${submit_command}" );
    # do not submit if doing debug - just return
    if( defined($$cmd_ref{debug}) ){
        return( $ierr );
    }

    # hack to fix environment setting for environment variables
    # for cielo/cielito flavor of moab since it does not pass env vars set to ""
    # correctly
    if( $$cmd_ref{sys_info}{L_MACHINE} eq "CIELITO" ||
        $$cmd_ref{sys_info}{L_MACHINE} eq "CIELO" ){
        foreach $env_name ( keys %ENV ){
            if( $ENV{$env_name} eq "" ){
                $ENV{$env_name} = $RJ_ENV_EMPTY;
            }
        }
    }

    # keep submitting until it works.
    $done = "false";
    while( $done eq "false" ){
        $done = "true";
        $output =  `($submit_command) 2>&1`;
        if( ! defined($output) ){
            $output = "";
        }
        # on roadrunner, the jobs get killed for some reason and no output
        if( $output !~ /\S/ ){
            $done = "false";
            $ierr = 0;
            &print_error( "No output from submission.",
                          "[$submit_command]",
                          "Retrying...",
                          $ierr );
        }
        # if need to resubmit, sleep for a bit and try again
        if( $done ne "true" ){
            sleep( 5 );
        }
    }
    #chomp( $output );
    # sometimes moab messes up but job is really submitted
    # try to get job id from run_job.pl (using last column being PNAME)
    if( $output =~ /communication error/i ){
        if( $RUN_STATUS_PL =~ /\S/ ){
            # might need to sleep for a bit - but see if not needed
            #sleep( 10 );
            $output1 = `$RUN_STATUS_PL 2>&1`;
            @lines = reverse split(/\n/, $output1);
            @lines = grep( /\s+$$datafile_ref{var}{PNAME}\s*$/, @lines);
            if( $#lines >= 0 ){
                $output = $lines[0];
            }
        }
    }
    if( $output =~ /\berror\b/i ){
        $ierr = 1;
        $output =~ s/\s+$//;
        &print_error( "Error in submission. Message:",
                      $output,
                      $ierr );
        exit( $ierr );
    }
    $output =~ s/^\s*(.*?)\s*$/$1/;
    &print_message( $cmd_ref, "Output", $output );
    # remove warning messages about uninitialized variable that the lanl msub script has:
    #   /opt/MOAB/lanl/msub_filter/msub_filter.pl
    ( $output1 = $output ) =~ s/Use of uninitialized[^\n]*//g;
    # on Cielo, if you make a submit that gives you an error, the submit_id is of the form:
    #   moab\.\d+
    # So, need to grab all of that so that the --wait will look for correct id
    if( $output1 =~ /(\S*\d+\S*)/ ){
        $$cmd_ref{submit_id} = $1;
    }
    else{
        $ierr = 1;
        &print_error( "Cannot find dependent job id - aborting.",
                      "Manual cleanup of job might be required",
                      $ierr );
        exit( $ierr );
    }
    # create first file that has id of next process
    # future ids from other submissions will be in the file
    `echo $$cmd_ref{submit_id} > $RJ_FILE_ID_NEXT`;
    # add to log file now that know child id - must do before wait afterwards
    &add_log( $datafile_ref, $cmd_ref );
    # if wait
    if( defined($$cmd_ref{wait}) ){
        # root process waits for things to finish
        # give it a file for root process so all childs finish
        if( ! defined($$cmd_ref{id}) ){
            &my_wait( $cmd_ref, FILE=>$RJ_FILE_ID_NEXT );
        }
    }

    return( $ierr );
}

##############################################################################3

sub rj_submit_nobatch{
    my(
       $datafile_ref,
       $cmd_ref,
       $wrapper_file,
       ) = @_;
    my(
        $cmd,
        $depend,
        $dir,
        $submit_command,
        $wrapper_wait,
        );
    #...fork/exec code:
    my $pid_child;
    my $pid_parent = $$;
    my $done = "false";
    my $ierr = 0;
    my $note = "";
    $submit_command = "Internal Perl fork/exec [$wrapper_file]";
    if( defined($$cmd_ref{debug}) ){
        $note = "[DEBUG] ";
    }
    if( defined($cmd{v}) ){
    &print_message( $cmd_ref, "Command", "${note}${submit_command}" );
    }
    if( defined($$cmd_ref{debug}) ){
        return;
    }

    # fork child to wait for parent to finish then exec command
    # parent continues on
    while( $done eq "false" ){
        # parent block successful fork
        if( $pid_child=fork ) {
	  if( defined($cmd{v}) ){
            &print_message( $cmd_ref, "Finish_Submit", "" );
	  }
            $done = "true";
            $$cmd_ref{submit_id} = "$pid_child:PID";
            # create file that has id of next process
            `echo $$cmd_ref{submit_id} > $RJ_FILE_ID_NEXT`;
            # add to log file now that know child id - must do before wait afterwards
            &add_log( $datafile_ref, \%cmd );
            # if wait
            if( defined($$cmd_ref{wait}) ){
                # root process waits for things to finish
                if( ! defined($$cmd_ref{id}) ){
                    # this waits for first fork/child to finish and removes a troublesome defunct process
                    wait;
                    # wait for latest id to finish
                    &my_wait( $cmd_ref, FILE=>$RJ_FILE_ID_NEXT );
                }
            }
        }
        # child block
        elsif( defined $pid_child ) {
            $done = "true";
            my %pids_depend;
            my $pid;
            # if waiting, add parent to list to wait if not first one
            # do not wait on parent in this case since this top level
            # parent might be doing a "wait" for entire chain to be
            # finished before returning control
            if( defined($$cmd_ref{wait}) ){
                if( defined($$cmd_ref{id}) ){
                    $pids_depend{"$pid_parent:PID"} = "";
                }
            }
            # if not waiting, all parents (even root one) finish first
            else{
                $pids_depend{"$pid_parent:PID"} = "";
            }
            # if given additional dependencies, add to list
            if( defined($$cmd_ref{depend}) && $$cmd_ref{depend} =~ /\S/ ){
                foreach $depend ( split( /,/, $$cmd_ref{depend} ) ){
                    $pids_depend{"$depend"} = "";
                }
            }
            $wrapper_wait = join( ",", keys %pids_depend );
            if( defined( $wrapper_wait ) && $wrapper_wait =~ /\S/ ){
                $wrapper_wait = "--wrapper_wait $wrapper_wait";
            }

            # call run_job.pl with special args that will wait and run the wrapper
            # make an "exec" so that the arguments can be parsed with run_status.pl
            $dir = getcwd();
            # no --id since we want to skip things that id does (write $RJ_FILE_ID)
            $cmd = "run_job.pl --dir $dir --wrapper_file $wrapper_file --wrapper_out $$datafile_ref{submit_file} $wrapper_wait";
            exec( $cmd );
            exit( $ierr );
        }
        # unsuccessful fork, sleep for a bit and retry
        elsif( $! =~ /No more process/ ) {
            ErrorMessage::Print("Retrying fork in a few moments [$!]\n");
            sleep 5;
        }
        else {
            # weird fork error
            ErrorMessage::Print("Cannot fork [$!]\n");
            exit( 1 );
        }
    }
}

################################################################################

# wait for batch ID(s) or PID(s) to finish
sub my_wait{
    my( $cmd_ref ) = shift( @_ );
    my %args = (
                FILE  => undef, # single ID given by contents of the id files
                ID    => undef, # single ID
                IDS   => undef, # hash of IDS - need to check differently
                @_,
                );
    my $args_valid = ("FILE|ID|IDS");
    my(
        $arg,
        $command,
        $cond,
        %dirs_jobid,
        $done_id,
        $done_wait,
        $id_full,
        $id_wait,
        %ids_full,
        $ierr,
        $ispid,
        $jobid,
        @jobids,
        %names_jobid,
        $output,
        $note,
        %run_status_info,
        $state,
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
    # predefined id
    if( defined($args{ID}) ){
        $ids_full{$args{ID}} = "";
    }
    # predefined ids
    elsif( defined($args{IDS}) ){
        %ids_full = %{$args{IDS}};
    }
    $done_wait = "false";
    while( $done_wait eq "false" ){
        # init final done to be true (and changed if not)
        $done_wait = "true";
        # if id is gotten from file that is updated 
        if( defined( $args{FILE} ) ){
            # in case something funky going on and previous jobs echo did not work yet
            $id_full = "";
            while( $id_full eq "" ){
            while( ! -e $args{FILE} ){
                sleep( 1 );
            }
            # make sure echo from previous job is done
            sleep( 1 );
            $id_full = `cat $args{FILE}`;
            $id_full =~ s/\s+$//;
                # Could be problem if file created but invalid data
                #   filesystem full (so 0 sized file)
                #   job killed just at the moment file created
                # Set id_full to something so at least no hang.
                # If it seems that this is common, will have to rethink.
                if( $id_full !~ /\S/ ){
                    &print_message( $cmd_ref, "Warning", "empty:$args{FILE}" );
                    $id_full = "empty:$args{FILE}";
                }
            }
            # now only look at latest id since that is what you want to wait for
            %ids_full = ( "$id_full" => "" );
        }
 
        # go through each id
        foreach $id_full ( keys %ids_full ){
            # init this id to be done (and changed if not)
            $done_id = "true";
            # get id and type of id (batch or PID)
            if( $id_full =~ /^(\S+):PID$/ ){
                $id_wait = $1;
                $ispid = "yes";
            }
            else{
                $id_wait = $id_full;
                $ispid = "no";
            }
            # if pid
            if( $ispid eq "yes" ){
                if( $id_wait =~ /^\d+$/ ){
                    # do not use run_status.pl since can wait on arbitrary pids
                    $command = "ps $id_wait";
                    $output = `$command 2>&1`;
                    if( $output =~ /^\s*($id_wait)\s+/m ){
                        $done_id = "false";
                    }
                }
            }
            # if batch id
            else{
                undef( %run_status_info );
                undef( %dirs_jobid );
                # get info about running jobs
                $ierr = &get_run_status_info( \%run_status_info, \%dirs_jobid, \%names_jobid, "-j ${id_wait}" );

                # if error condition, then will try again (batch system down)
                if( $ierr != 0 ){
                    &print_message( $cmd_ref, "Warning", "retry: get_run_status_info" );
                    $done_id = "false";
                }
                # check returned jobs
                elsif( keys %run_status_info ){
                    @jobids = keys %run_status_info;
                    $jobid = $jobids[0];
                    $state = $run_status_info{$jobid}{vals}{STATE};
                    $note = $run_status_info{$jobid}{vals}{NOTE};
                    $cond = "";
                    if( defined($state) ){
                        $cond .= " $state";
                    }
                    if( defined($note) ){
                        $cond .= " $note";
                    }
                    if( $cond !~ /complete/i &&
                        $cond !~ /removed/i ){
                        $done_id = "false";
                    }
                }
                # throttle here (also below) since batch commands stress system
                if( $done_id eq "false" ){
                    sleep( 10 );
                }
            }
            # if that id is done, remove from list
            if( $done_id eq "true" ){
	      if( defined($cmd{v}) ){
                &print_message(  $cmd_ref, "Waiting:done", $id_full );
	      }
	      delete( $ids_full{$id_full} );
            }
            # otherwise, final done is false
            else{
                $done_wait = "false";
            }
        }
        # if still ids, sleep a bit then check
        if( $done_wait eq "false" ){
            # might not need but throttles moab checks
            sleep( 2 );
        }
    }
    if( defined($cmd{v}) ){
      &print_message( $cmd_ref, "Waiting:done", "ALL" );
    }
}

#............................................................................
#...Name
#...====
#... rj_script
#...
#...Purpose
#...=======
#... runs the script
#...
#...Arguments
#...=========
#... $cmd_ref     Intent: out
#...              Perl type: reference to hash
#...              command line
#...              $cmd{$option} = value
#...              $cmd{files}[] = array of file names
#...
#...Program Flow
#...============
#... 1) 
#............................................................................
sub rj_script{
    my(
        $datafile_ref,
        $cmd_ref,
        ) = @_;
    my(
       $action,
       $batch_file,
       $cmd,
       $cwd,
       $date,
       $date_args,
       $done,
       $env_name,
       $env_name_val,
       $env_string,
       $env_val,
       $exn_left,
       $exn_needed,
       $file,
       %found,
       $i,
       $i_last,
       $i_start,
       $ierr,
       $key,
       $line,
       @lines_match,
       $match_string,
       $msg,
       $nid,
       @nids,
       $node,
       $node_name,
       %nodes,
       $num_bad,
       $out,
       $output,
       $output_script,
       @output_script_lines,
       $pid_child,
       $pipe_use,
       $retry_found,
       $script_out,
       $status,
       $time_current,
       $time_elapsed,
       $time_remaining_b,
       $try,
       @try_times,
       $try_time_delta,
       $type,
       );


    # run $RJ_FILE_CMD until finished
    $cwd = getcwd();
    $done = "false";
    $retry_found = "false";
    $try = 0;
    $exn_left = $$datafile_ref{var}{EXN};
    # set environment variables
    $ENV{RJ_RUN_JOB}      = "yes";
    $ENV{RJ_TIME_START}   = $$cmd_ref{time_start};
    $ENV{RJ_TIME}         = $$datafile_ref{var}{TIME_SECS};
    $ENV{RJ_VAR_PNAME}    = $$datafile_ref{var}{PNAME};
    $ENV{RJ_VAL_DIR_ROOT} = $RJ_VAL_DIR_ROOT;
    $ENV{RJ_VAL_DIR}      = $RJ_VAL_DIR;
    if( $$datafile_ref{var}{GROUP} =~ /\S/ ){
        $ENV{RJ_VAR_GROUP}      = $$datafile_ref{var}{GROUP};
    }
    if( defined($$datafile_ref{var}{UMASK}) ){
        $ENV{RJ_VAR_UMASK}      = $$datafile_ref{var}{UMASK};
    }
    foreach $key ( keys %{$cmd{sys_info}} ){
        $ENV{"RJ_$key"} = $cmd{sys_info}{$key};
    }
    $ENV{RJ_ID} = &my_getid( $cmd_ref );
    # since is called from first one, submit_id might not be defined
    if( ! defined($$cmd_ref{submit_id}) ){
        $ENV{RJ_ID_CHILD} = -1;
    }
    else{
        $ENV{RJ_ID_CHILD} = $$cmd_ref{submit_id};
    }

    # while loop for retrying in same batch submission
    while( $done eq "false" ){

        # get new tag if doing another try
        if( $try > 0 ){
            &new_tag( $cmd_ref, $try );
        }

        # set some environment variables
        $ENV{RJ_TAG} = $$cmd_ref{tag_whole};
        $out = `date +%s`;
        chomp( $out );
        $ENV{RJ_TIME_CURRENT}     = $out;
        $ENV{RJ_TIME_ELAPSED}     = $ENV{RJ_TIME_CURRENT}   - $ENV{RJ_TIME_START};
        $ENV{RJ_TIME_REMAINING}   = $ENV{RJ_TIME}           - $ENV{RJ_TIME_ELAPSED};
        $ENV{RJ_TIME_REMAINING_B} = $ENV{RJ_TIME_REMAINING} - $$datafile_ref{var}{TIME_B};
        $ENV{RJ_TRY}              = $try;

        # any --env setting
        foreach $env_name_val ( @{$$cmd_ref{env}} ){
            ($env_name, $env_val) = split(/=/, $env_name_val);
            $ENV{$env_name} = $env_val;
        }

        # set output file name and init it (unlink) and set soft link
        $script_out = "${RJ_DIR}/${RJ_FILE_CMD_OUT}.$$datafile_ref{var}{PNAME}.$$cmd_ref{tag_whole}";
        unlink( $script_out );
        `ln -sf $script_out ./${RJ_FILE_CMD_OUT}`;

        # get date and put into output file
        $date = `date`;
        chomp( $date );
        # mac does not like "date -d" - try to figure out if failure
        $date_args = "-d";
        $out = `date $date_args '$date' +'RJ_OUTPUT: TIME_S:%s' 2>&1`;
        if( $out !~ /RJ_OUTPUT/ ){
            $date_args = "-jf '\%a \%b \%e \%k:\%M:\%S \%Z \%Y' ";
        }
        `date $date_args '$date' +'RJ_OUTPUT: TIME_S:%s'                     >> $script_out`;
        `date $date_args '$date' +'RJ_OUTPUT: TIME_YMDHMS:%Y.%m.%d.%H.%M.%S' >> $script_out`;

        # pid into output file
        `echo RJ_OUTPUT: PID: \$\$ >> $script_out`;

        # environment vars into output file
        $env_string = "";
        foreach $key ( sort keys %ENV ){
            if( $key =~ /^RJ_/ ){
                $env_string .= "RJ_OUTPUT: ENV $key = $ENV{$key}\n";
            }
        }
        if( open( FILE, ">>$script_out" ) ){
            print FILE $env_string;
            close( FILE );
        }

        # >> or | tee (for interactive)
        $pipe_use = ">>";
        if( defined($$cmd_ref{i}) ){
            $pipe_use = "| tee -a";
        }

        # command to run, and run it
        # both "2>&1" needed to send to file when doing ">>" and "| tee -a"
        $cmd = "./$RJ_FILE_CMD_BASE 2>&1 $pipe_use $script_out 2>&1";
        &print_message( $cmd_ref, "$RJ_FILE_CMD_BASE", "running ./$RJ_FILE_CMD_BASE" );
        &run_command( COMMAND=>$cmd, VERBOSE=>'yes', STDOUT=>$$cmd_ref{i} );

        # final date into output file
        $date = `date`;
        chomp( $date );
        `date $date_args '$date' '+RJ_OUTPUT: TIME_S:%s'                     >> $script_out`;
        `date $date_args '$date' '+RJ_OUTPUT: TIME_YmdHMS:%Y.%m.%d.%H.%M.%S' >> $script_out`;
        $time_elapsed = `date $date_args '$date' '+%s'`;
        chomp( $time_elapsed );
        $time_elapsed = $time_elapsed - $ENV{RJ_TIME_CURRENT};
        `echo 'RJ_OUTPUT: TIME_ELAPSED_S:$time_elapsed' >> $script_out`;

        #................
        #...get status...
        #................
        $output_script = `cat $script_out 2>&1`;
        @output_script_lines = split( /\n/, $output_script );
        undef( %found );
        # -----------------------------------------------------------
        # ----------------- PROBLEMS --------------------------------
        # -----------------------------------------------------------
        # RJ_STOP problems - will never work
        $match_string =
            '('.
            'error while loading shared libraries|'.               # wrong exec
            'bytes requested; not enough memory|'.                 # memory allocation error
            'Unable to allocate required memory: std::bad_alloc|'. # memory allocation error
            'We were unable to find an allocation for this job|'.  # bad submit params
            'Your job has requested a conflicting number of processes|'. # bad process count
            'orte_rml_base_select'.                                # dunno what wrote this - but was bad
            ')';
        if( $output_script =~ /${match_string}/ ){
            push( @{$found{RJ_STOP}}, $1 );
        }
        # RJ_NEW_CHILD problems - try a new set of nodes
        $match_string =
            '('.
            'Cpuset file[^\n]*wrote -1|'. # currently cannot shut off if more than 1 node gone
            'apsched: claim exceeds reservation\'s nodes|'. # moab might remove bad nodes
            'terminate\s+signal\s+SIGURG\s+issued|'. # cielo SIGURG means job needs to stop
            'Network is down'.            # network issues
            ')';
        if( $output_script =~ /${match_string}/ ){
            push( @{$found{RJ_NEW_CHILD}}, "SYSFAIL: $1" );
        }
        # RJ_RETRY problems - try again with same batch submission
        $match_string = 
            '('.
            'RJ_RETRY|'. # found RJ_RETRY in output (eg, from HANG-DETECTED from run_job_cleanup.pl)
            'Received node failed|'.
            'ERROR - nem_gni_error_handler|'.   # cielo
            'GNI transaction error detected|'.  # cielo
            'Unable to connect to[^\n]port[^\n]socket[^\n]Connection timed out|'. # cielo
            'currently cannot shut off if more than 1 node gone...Cpuset file[^\n]*wrote -1|'. # cielo
            'close of the compute node connection before app startup barrier|'. # cielo
            'Received node failed or halted event|'. # cielo
            'Cray HSN detected critical error|'. #cielo
            'mca_oob_tcp_msg_recv|'.            # openmpi
            'connection reset by peer|'.        # openmpi
            'to\s+launch\s+so\s+we\s+are\s+aborting|'.    # openmpi
            'RETRY EXCEEDED ERROR|'.            # roadrunner
            'pread[^\n]*Input\/output error|'.   # filesystem strangeness
            'Delivering SIGKILL with status LOADING timed out|'. # sequoia startup problem
            'abnormal termination by signal 9 from rank'. # sequoai strange kill then worked on try 17
            ')';
        if( $output_script =~ /${match_string}/ ){
            # keep track that a special retry was found (for max retries)
            $retry_found = "true";
            $type = $1;
            push( @{$found{RJ_RETRY}}, "SYSFAIL: $1" );
            # remove bad nodes
            undef( %nodes );
            undef( @nids );
            # get array of nodes to take out
            if( $$cmd_ref{sys_info}{L_CLASS} eq "CIELO" &&
                ( $type eq "GNI transaction error detected" ||
                  $type =~ /pread[^\n]*Input\/output error/ ) ){
                # get lines that have bad nodes in them
                $found = "false";
                $match_string = 
                    '('.
                    'GNI transaction error detected|'.  # cielo
                    'Application called MPI_Abort'.     # from pread filesystem badness
                    ')';
                @lines_match = grep( /${match_string}/, @output_script_lines );
                foreach $line ( @lines_match ){
                    if( $line =~ /(c\d{1,2}-\dc\ds\dn\d)/ ){
                        $nodes{$1} = "";
                        $found = "true";
                    }
                }
                # translate to node id = nid
                if( $found eq "true" ){
                    $out = `xtprocadmin 2>&1`;
                    foreach $node_name ( keys %nodes ){
                        if( $out =~ /^\s*(\d+).*${node_name}/m ){
                            $node = $1;
                            if( $node =~ /^0(\d+)/ ){
                                $node = $1;
                            }
                            push( @nids, $node );
                        }
                    }
                }
            } # get array of nodes to take out
            if( $$cmd_ref{sys_info}{L_CLASS} eq "CIELO" &&
                $type =~ /Cpuset file.*wrote -1/ ){
                $match_string = 
                    'Cpuset file[^\n]*wrote -1';  # cielo
                @lines_match = grep( /${match_string}/, @output_script_lines );
                foreach $line ( @lines_match ){
                    if( $line =~ /NID\s+(\d+)/ ){
                        push( @nids, $1 );
                    }
                }

            }
            $match_string = 
                'Unable to connect to[^\n]port[^\n]socket[^\n]Connection timed out|'.
                'received node failed or halted event';
            if( $$cmd_ref{sys_info}{L_CLASS} eq "CIELO" &&
                $type =~ /$match_string/ ){
                @lines_match = grep( /${match_string}/, @output_script_lines );
                foreach $line ( @lines_match ){
                    if( $line =~ /NID\s+(\d+)/ ){
                        push( @nids, $1 );
                    }
                }

            }
            # unique sort them to get rid of duplicates
            @nids = sort_unique( \@nids );
            $num_bad = $#nids + 1;
            # warn if not enough extra nodes - and signal a new child
            if( $num_bad > $exn_left ){
                $ierr = 0;
                $exn_needed = $$datafile_ref{var}{EXN} + $num_bad - $exn_left;
                &print_error( "$num_bad bad node(s) detected",
                              "$exn_left extra nodes remaining",
                              "'#RJ EXN = $$datafile_ref{var}{EXN}' - perhaps make it $exn_needed",
                              $ierr );
                push( @{$found{RJ_NEW_CHILD}}, "SYSFAIL: bad node but cannot remove" );
            }
            # remove nodes
            else{
                # on cielo can remove nodes
                if( $$cmd_ref{sys_info}{L_CLASS} eq "CIELO" ){
                    foreach $nid ( @nids ){
                        $exn_left--;
                        $cmd = &rj_get_run_command( EXEC=>"sleep",
                                                    NUMPE=>"1",
                                                    TYPE=>"mpi",
                                                    EXEC_ARGS=>"24h",
                                                    SHELL=>"sh",
                                                    PRUN_ARGS=>"-L $nid",
                                                    OPT=>"$$datafile_ref{var}{OPT} no_mapping",
                            );
                        &print_message( $cmd_ref, "Downing bad node", "$cmd" );
                        # parent block
                        if( $pid_child=fork ){
                        }
                        # child block
                        elsif( defined $pid_child ){
                            $out = `$cmd 2>&1`;
                            &print_message( $cmd_ref, "Downing node output", $output );
                            exit(0);
                        }
                    }
                    # give time for the child prun(s) to hit and to remove nodes
                    if( $#nids > 0 ){
                        sleep( 5 );
                    }
                } # on cielo can remove nodes
            } # remove nodes
        } # RJ_RETRY problems - try again with same batch submission

        # -----------------------------------------------------------
        # ----------------- TOO MANY RETRIES  -----------------------
        # -----------------------------------------------------------
        # check for too many tries now after any problems discovered
        #   but before other RJ_RETRY checks
        #   This will get any RJ_RETRY from current submission, but might
        #   need to check earlier submissions, or RJ_RETRY from output
        #   or files found before.
        # stop if too many tries in one submission
        if( $$datafile_ref{var}{TRY_NUM_MAX} > 0 &&
            $try + 1 >= $$datafile_ref{var}{TRY_NUM_MAX} ){
            # if found RJ_RETRY, do RJ_NEW_CHILD
            if( $retry_found eq "true" ){
                # go through previous runs and see if they all had TRY_NUM_MAX
                if( $#{$$cmd_ref{tags}} >= 1 ){
                    $i_last = $#{$$cmd_ref{tags}} - 1;
                    $i_start = $i_last - $DEFAULT_TRY_NUM_MAX_NUM + 1;
                    # do not start if not enough batches yet
                    if( $i_start < 0 ){
                        $i_start = $i_last + 1;
                        $found = "false";
                    }
                    else{
                        $found = "true";
                    }
                    for( $i = $i_start; $i <= $i_last; $i++ ){
                        $batch_file = "$RJ_DIR/$RJ_FILE_BATCH_OUT_BASE.$$cmd_ref{tags}[$i]";
                        $cmd = "grep 'RJ_NEW_CHILD.*TRY_NUM_MAX' $batch_file";
                        $out = `($cmd)2>&1`;
                        # if hit file that did not have this issue, stop
                        if( $out !~ /\S/ ){
                            $found = "false";
                            last;
                        }
                    }
                    if( $found eq "true" ){
                        push( @{$found{RJ_STOP}}, "Current and previous $DEFAULT_TRY_NUM_MAX_NUM batches exceeded TRY_NUM_MAX" );
                    }
                }
                $action = "RJ_NEW_CHILD";
            }
            else{
                $action = "RJ_STOP";
            }
            push( @{$found{$action}}, "TRY_NUM_MAX [$$datafile_ref{var}{TRY_NUM_MAX}] exceeded" );
        }

        # stop if too many tries in a short time period
        if( $try == 0 ){
            push( @try_times, $ENV{RJ_TIME_CURRENT} );
        }
        push( @try_times, $ENV{RJ_TIME_CURRENT} + $time_elapsed );
        if( $$datafile_ref{var}{TRY_NUM} > 0 &&
            $try + 1 >= $$datafile_ref{var}{TRY_NUM} ){
            # shift off first one to keep TRY_NUM try times in array
            $try_time_delta = shift( @try_times );
            $try_time_delta = $try_times[-1] - $try_time_delta;
            if( $$datafile_ref{var}{TRY_TIME} > 0 &&
                $try_time_delta < $$datafile_ref{var}{TRY_TIME} ){
                # if found RJ_RETRY, do RJ_NEW_CHILD
                if( $retry_found eq "true" ){
                    $action = "RJ_NEW_CHILD";
                }
                else{
                    $action = "RJ_STOP";
                }
                push( @{$found{$action}}, "Time for TRY_NUM=$$datafile_ref{var}{TRY_NUM} tries less than minimum [$try_time_delta < TRY_TIME=$$datafile_ref{var}{TRY_TIME} seconds]" );
            }
        }
        # -----------------------------------------------------------
        # ----------------- OTHER  ----------------------------------
        # -----------------------------------------------------------

        # check for existence of files - remove them
        foreach $file ( "RJ_RETRY", "RJ_NEW_CHILD", "RJ_STOP" ){
            if( -e "$file" ){
                push( @{$found{$file}}, "file $file" );
                unlink( $file );
            }
        }

        # RJ_<status> in output
        $cmd = "egrep 'RJ_STOP|RJ_NEW_CHILD|RJ_RETRY' $script_out";
        $out = `($cmd) 2>&1`;
        if( $out =~ /\S/ ){
            foreach $status ( "RJ_STOP", "RJ_NEW_CHILD", "RJ_RETRY" ){
                if( $out =~ /($status)/ ){
                    push( @{$found{$1}}, "$cmd" );
                }
            }
        }

        # RJ_RETRY if enough time left
        if( $$datafile_ref{var}{CHAIN} eq "yes" ){
            $out = `date +%s`;
            chomp( $out );
            $time_current = $out;
            $time_elapsed = $time_current - $$cmd_ref{time_start};
            $time_remaining_b = $ENV{RJ_TIME} - $time_elapsed - $$datafile_ref{var}{TIME_B};
            if( $time_remaining_b > 0 ){
                push( @{$found{RJ_RETRY}}, "RJ_TIME_REMAINING_B=$time_remaining_b" );
            }
        }

        # if doing interactive, RJ_STOP regardless since we want a ctrl-c to break
        # out of the run completely and not RJ_RETRY
        if( defined($$cmd_ref{i}) ){
            push( @{$found{RJ_STOP}}, "-i (interactive) flag used" );
        }

        # if not chaining, RJ_STOP
        if( $$datafile_ref{var}{CHAIN} ne "yes" ){
            push( @{$found{RJ_STOP}}, "#RJ CHAIN = no (default)" );
        }

        # =============================
        # finally do appropriate action
        # print all findings
        # =============================
        foreach $type ( "RJ_RETRY", "RJ_NEW_CHILD", "RJ_STOP" ){
            foreach $msg ( @{$found{$type}} ){
                &print_message( $cmd_ref, "$type", "$msg");
            }
        }
        $done = "true";
        if( $#{$found{RJ_STOP}} >= 0 ) {
            $done = "true";
            &job_killdepend( $datafile_ref, $cmd_ref );
        }
        elsif( $#{$found{RJ_NEW_CHILD}} >= 0 ) {
            $done = "true";
        }
        elsif( $#{$found{RJ_RETRY}} >= 0 ) {
            $done = "false";
        }

        # clean up after script is run
        &job_cleanup($datafile_ref, $cmd_ref);

        $try++;
    }
}

#...job_cleanup: clean up when running again
sub job_cleanup{
    my(
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
        $clean_mpi,
        $cmd,
        $cmd_name,
        @lines,
        $out,
        $ORTE_CLEAN,
        $pid,
        $tmpfile,
        @tmpfiles,
        $version,
        );

    # clean_mpi: do no matter what
    if( ! defined($clean_mpi) ){
        if( $$cmd_ref{clean_mpi} =~ /\S/ ){
            $clean_mpi = $$cmd_ref{clean_mpi};
        }
    }

    # clean_mpi: do not if on back end
    if( $$cmd_ref{sys_info}{L_END} ne "BACK" ){
        $clean_mpi = "no";
    }

    # clean_mpi: see if should
    if( ! defined($clean_mpi) ){
        # set to no...until conditions below change
        $clean_mpi = "no";
        # if not SERIAL
        if($$datafile_ref{var}{SERIAL} eq "no"){
            # do not clean if other run_job.pl is being run (orte-clean blows away all mpi)
            $cmd = 'ps -w -w -eo "pid user args"';
            $out = `$cmd 2> /dev/null`;
            $pid = $$;
            @lines = split(/\n/, $out);
            @lines = grep ( /^\s*\d+\s+.*run_job\.pl/, @lines);
            @lines = grep ( !/^\s*${pid}\s+/, @lines);
            if( $#lines < 0 ){
                $clean_mpi = "yes";
            }
        }
    }
    
    # if clean_mpi, do it
    if( $clean_mpi eq "yes" ){
      # run orte-clean if it exists
      $ORTE_CLEAN = &which_exec( "orte-clean", QUIET=>"" );
      # for now (perhaps forever if I never look at this comment again), turn
      # off orte-clean if openmpi/1.7.x or newer - it hangs using that.
      if( $ORTE_CLEAN ne "" ){
	( $version = $ORTE_CLEAN ) =~ s&.*(openmpi[/\d.\-]+).*&$1&;
	if( &my_compare_version( $version, "<", "1.7" ) ){
	  &print_message( $cmd_ref, "Cleaning", $ORTE_CLEAN  );
	  $cmd = &rj_get_run_command( EXEC=>$ORTE_CLEAN,
				      PATH_SET=>"",
				      PPN=>1,
				      TYPE=>"mpi",
				      SHELL=>"sh",
				      OPT=>"$$datafile_ref{var}{OPT} no_mapping" );
	  $out = &run_command( COMMAND=>$cmd );
	  chomp( $out );
	  print "[$out]\n";
        }
      }

        # llnl machines use srun and it sometimes does not clean up
        # Kill the first one - seems to be the master
        $out = `ps`;
        $cmd_name = "srun";
        if( $out =~ /^\s*(\d+)\s.*\s${cmd_name}\s*$/m ){
            $pid = $1;
            &print_message( $cmd_ref, "Cleaning", "kill ${pid} (${cmd_name})" );
            $cmd = "kill $pid";
            $out = &run_command( COMMAND=>$cmd );
            chomp( $out );
            print "[$out]\n";
        }
    }

    # remove mapping file(s)
    $tmpfile = "${RJ_FILE_MAPPING_FULL}.*";
    @tmpfiles = glob($tmpfile);
    if( $#tmpfiles >= 0 ){
        &print_message( $cmd_ref, "Cleaning", "unlink $tmpfile" );
        unlink( @tmpfiles );
    }
}

#...
sub job_killdepend{
    my(
       $datafile_ref,
       $cmd_ref,
       ) = @_;
    my(
        $command,
        $output
        );
    if( defined($$cmd_ref{submit_id}) &&
        $$cmd_ref{submit_id} ne "-1" ){
        if( $$cmd_ref{submit_id} =~ /^(\d+):PID$/ ){
            $command = "kill $1";
        }
        else{
            $command = "$RUN_STATUS_PL --cancel -j $$cmd_ref{submit_id}";
        }
        &print_error( "Killing dependent job [BATCH=$$datafile_ref{var}{BATCH}] [submit_id=$$cmd_ref{submit_id}]",
                      "Command: $command",
                      0 );
        &print_message( $cmd_ref, "Command", "$command" );
        $output =  `($command) 2>&1`;
        $output =~ s/^\s*(.*?)\s*$/$1/;
        &print_message( $cmd_ref, "Output", "$output" );
        delete( $$cmd_ref{submit_id} );
    }
}

# special routine called from within $RJ_FILE_DATAFILE to fix the
# environment
sub rj_fix_env{
    my(
        $cmd_ref,  # general hash
        $ex,       # executable
        $shell,    # shell
        $op,       # special op (for recursion)
        ) = @_;
    my(
        $changed_compiler,
        $changed,
        $cmd,
        $cmd_pre,
        $cmd_setup,
        $ldd_ex,
        @ldd_ex_lines,
        $lib,
        $lib_found,
        $loadedmodules,
        @loadedmodules_lines,
        $load_module,
        $load_module_default,
        $load_message,
        $module_use,
        $output,
        @output_lines,
        $try,
        $unload_module,
        $val,
        $version,
        $version_loaded,
        $version_max,
        $version_min,
        );

    # init
    if( ! defined($shell) ){
        $shell = $RJ_SHELL;
    }
    $try = 0;
    $cmd = "";
    $cmd_pre = "";
    # wow...some of the "not found" messages are sent to STDERR
    # so need to capture that as well.
    if( $LDD ne "" ){
        $ldd_ex = `$LDD $ex 2>&1 `;
    }
    else{
        $ldd_ex = "";
    }
    @ldd_ex_lines = split( /\n/, $ldd_ex );
    $loadedmodules = $ENV{LOADEDMODULES};
    if( defined($loadedmodules) ){
        @loadedmodules_lines = split(/:/, $loadedmodules);
    }
    else{
        $loadedmodules = "";
        @loadedmodules_lines = ();
    }

    # ---------------------------------------
    # compile modules - do before mpi modules
    # ---------------------------------------

    # if any compiler was loaded, need to reload the mpi
    $changed_compiler    = "";

    # if intel
    $load_module         = "";
    $load_message        = "";
    $load_module_default = "";
    $unload_module       = "";
    $module_use          = "";
    $version_loaded      = "";
    $version_max         = "";
    $version_min         = "";
    if( $ldd_ex =~ /\b(libifcore(\.\S)*\.so(\.\d+)*)\s*=>\s*(\S.*)/ ){
        $lib_found = $4;

        # only figure out what to load if not found
        if( $lib_found =~ /not\s+found/i ){
        $load_module_default = "intel";
        
        # get version_loaded
        @output_lines = grep( m&^(compilers/)?intel[\/-]&i, @loadedmodules_lines );
        $version_loaded = $output_lines[0];
        if( defined($version_loaded) ){
            if( $version_loaded =~ m&intel[\-\.]([\d\.]+)& ){
                $version_loaded = "intel/$1";
            }
        }
        else{
            $version_loaded = "";
        }

        # if does     have libirng, needs to be at least intel/13.0.1
        if( $ldd_ex =~ /\b(libirng)\S*.so/ ){
            $version_min = "intel/13.0.1";
        }
        # if does not have libirng, needs to be at most  intel/12.1.5
        else{
            $version_max = "intel/12.1.5";
        }
    }
    }
    # fill out module commands
    $unload_module = "intel";
    # no need to reload this
    $changed = "";
    $cmd_pre .= &rj_fix_env_mod( $shell, $load_module, $load_message, $load_module_default, $unload_module,
                                 $module_use, $version_loaded, $version_max, $version_min,
                                 \$try, \$changed );
    if( $changed ne "" ){
        $changed_compiler = $changed
    }

    # if pgi
    $load_module         = "";
    $load_message        = "";
    $load_module_default = "";
    $unload_module       = "";
    $module_use          = "";
    $version_loaded      = "";
    $version_max         = "";
    $version_min         = "";
    if( $ldd_ex =~ /\b(libpgf\S*(\.\S)*\.so(\.\d+)*)\s*=>\s*(\S.*)/ ){
        $lib_found = $4;

        # check for additional lib not found
        if( $ldd_ex =~ /libaccapid.so\s*=>\s*(not\s+found)/ ){
            $lib_found = $1;
        }

        # only figure out what to load if not found
        if( $lib_found =~ /not\s+found/i ){
        $load_module_default = "pgi";
        
        # get version_loaded
        @output_lines = grep( m&^(compilers/)?pgi[\/-]&i, @loadedmodules_lines );
        $version_loaded = $output_lines[0];
        if( defined($version_loaded) ){
            if( $version_loaded =~ m&pgi[\-\.]([\d\.]+)& ){
                $version_loaded = "pgi/$1";
            }
        }
        else{
            $version_loaded = "";
        }

            # if has libaccapid
            if( $ldd_ex =~ /\b(libaccapid)\S*.so/ ){
                $version_min = "pgi/13.7";
            }
            else{
                $version_min = "pgi/13.2";
            }
        }
    }
    # fill out module commands
    $unload_module = "pgi";
    # no need to reload this
    $changed = "";
    $cmd_pre .= &rj_fix_env_mod( $shell, $load_module, $load_message, $load_module_default, $unload_module,
                                 $module_use, $version_loaded, $version_max, $version_min,
                                 \$try, \$changed );
    if( $changed ne "" ){
        $changed_compiler = $changed
    }

    # ---------------------------------------
    # mpi modules - do after compiler modules
    # ---------------------------------------
    $load_module         = "";
    $load_message        = "";
    $load_module_default = "";
    $unload_module       = "";
    $module_use          = "";
    $version_loaded      = "";
    $version_max         = "";
    $version_min         = "";

    # if openmpi
    if( $ldd_ex =~ /\b(libmpi(\.\S)*\.so(\.\d+)*)\s*=>\s*(\S.*)/ ){
      $lib = $1;
        $lib_found = $4;

        # Need correct mpi regardless of if ldd(exec) libs found or not.
        # If loaded old mpi but running exec with rpath (no libs found)
        # with new mpi, run will crash.
        # So, always check.
        if( $lib_found =~ /not\s+found/i || 1 == 1 ){
      $load_module_default = "openmpi";
      
      # get version_loaded
      @output_lines = grep( m&^(mpi/)?openmpi[/\-]&i, @loadedmodules_lines );
      $version_loaded = $output_lines[0];
      # strip off compiler part of string
      if( defined($version_loaded) ){
	if( $version_loaded =~ m&openmpi[\-\.]([\d\.]+)& ){
	  $version_loaded = "openmpi/$1";
	}
      }
      else{
	$version_loaded = "";
      }
      
      # if uses libmpi.so.0,  needs to be at most  openmpi/1.4.5
      if( &my_compare_version( $lib, "<=", "0" ) ){
	$version_max = "openmpi/1.4.5";
      }
      # if uses libmpi.so.1+, needs to be at least openmpi/1.6.3
      else{
	$version_min = "openmpi/1.6.3";
	# if fortran, need certain versions
	if( $ldd_ex =~ /libmpi_f|libmpi_usempif/ ){
	  # if has libmpi_usempif08, then at least openmpi/1.7.4
	  if( $ldd_ex =~ /libmpi_usempif08/ ){
	    $version_min = "openmpi/1.7.4";
	  }
	  else{
	    $version_max = "openmpi/1.6.5";
	  }
        }
      }
    }
    }
    # fill out module commands
    $unload_module = "openmpi";
    $changed = $changed_compiler; # if changed, will reload this module
    $cmd_pre .= &rj_fix_env_mod( $shell, $load_module, $load_message, $load_module_default, $unload_module,
                                 $module_use, $version_loaded, $version_max, $version_min,
                                 \$try, \$changed );

    # ---------------
    # perflib modules
    # old version do not have in rpath
    # ---------------
    $load_module         = "";
    $load_message        = "";
    $load_module_default = "";
    $unload_module       = "";
    $module_use          = "";
    $version_loaded      = "";
    $version_max         = "";
    $version_min         = "";
    # if PERFlib not found
    if( $ldd_ex =~ /\b(libperfrt(\.\S)*\.so(\.\d+)*)\s*=>\s*(not\s+found)/ ){
        $lib = $1;
        $version_loaded = "";
        $version_min    = "";
        $version_max    = "";
        $load_module_default = "PERFlib";

        # one module per machine
        if( $$cmd_ref{sys_info}{L_MACHINE} eq "CIELITO" ){
            $version_min .= "PERFlib/4.0-mpt-ct ;";
        }
        elsif( $$cmd_ref{sys_info}{L_MACHINE} eq "CIELO" ){
            $version_min .= "PERFlib/4.0-mpt-ci ;";
        }
        elsif( $$cmd_ref{sys_info}{L_MACHINE} eq "LUNA" ){
            $version_min .= "PERFlib/4.0-openmpi-lu ;";
        }
        elsif( $$cmd_ref{sys_info}{L_MACHINE} eq "TYPHOON" ){
            $version_min .= "PERFlib/4.0-openmpi-ty ;";
        }
        elsif( $$cmd_ref{sys_info}{L_CLASS} eq "TLCC" ){
            $version_min .= "PERFlib/4.0-openmpi-ml ;";
        }
        elsif( $$cmd_ref{sys_info}{L_CLASS} eq "TLCC1" ){
            $version_min .= "PERFlib/4.0-openmpi-ty ;";
        }
        $version_max = $version_min;
        $module_use = "$$cmd_ref{sys_info}{L_EAP_INSTALL_DIR_UP}/codeopt/PERF/modulefiles";
    }
    # fill out module commands
    $unload_module = "PERFlib";
    $changed = ""; # if changed, will reload this module
    $cmd_pre .= &rj_fix_env_mod( $shell, $load_module, $load_message, $load_module_default, $unload_module,
                                 $module_use, $version_loaded, $version_max, $version_min,
                                 \$try, \$changed );

    # --------
    # gcc libs
    # --------
    $load_module         = "";
    $load_message        = "";
    $load_module_default = "";
    $unload_module       = "";
    $module_use          = "";
    $version_loaded      = "";
    $version_max         = "";
    $version_min         = "";
    # if specific flavor of GLIBCXX
    if( $ldd_ex =~ /\bGLIBCXX.*not\s+found/ ){
        $load_module_default = "gcc";
    }
    # fill out module commands
    $changed = ""; # if changed, will reload this module
    $cmd_pre .= &rj_fix_env_mod( $shell, $load_module, $load_message, $load_module_default, $unload_module,
                                 $module_use, $version_loaded, $version_max, $version_min,
                                 \$try, \$changed );

    # ----------------------
    # set up module commands
    # ----------------------
    if( $cmd_pre ne "" ){
        # module command setup
        $cmd_pre = &rj_shell_module_setup($shell)." $cmd_pre";
        # and module list afterwards
        $cmd_pre .= &rj_shell_module_command( "list", $shell );

        # set to ingore conflicts (some teams use multiple compilers)
        # only if not already set (and will unset when done)
        if( ! defined($ENV{IGNOREMODULECONFLICTS}) ){
            $cmd_pre =  &rj_shell_var_set("IGNOREMODULECONFLICTS", "",    ";", $shell)." $cmd_pre";
            $cmd_pre .= &rj_shell_var_set("IGNOREMODULECONFLICTS", undef, ";", $shell);
        }
    }

    # get mpi version with module set loaded
    # get shell version of above commands
    # called recursively, so do not call if $op is set
    if( ! defined($op) ){
        # get shell flavor of the command and do mpirun -V
        $cmd_setup = &rj_fix_env( $cmd_ref, $ex, "sh", "prun" );
        $output = `($cmd_setup mpirun -V) 2>&1`;
        if( $output =~ /^\s*mpirun.*?((\d+\.?)+)/m ){
            $version = $1;
        }
        # LAP uses multiple mpiruns per submission and needs to have
        #   mpi_paffinity_alone=0.  So, set as default.  (see ppn comment
        #   above).
        # typhoon: mpi_paffinity_alone=0 and partial node horrible memory
        #          performance.  Run stress test with 5 second cpu and see.
        # Roadrunner: messes up cells - but no roadrunner now
        # openmpi/1.6.3: With flag, tests run 100x slower and give wrong answers!!
        #                If you always set OMP_NUM_THREADS, it seems to be ok.
        #                Will probably work in later versions...need to
        #                have something so that can run multiple mpirun's
        #                in the background and not oversubscribe the
        #                node.
        #                Tools.rh/Sanity/$RJ_FILE_DATAFILE shows this
        #                  with  flag, each takes < 5 secs
        #                  w/out flag, each takes 40+ secs.
        if( &my_compare_version( $version, "<", "1.7" ) ){
            $val = " -mca mpi_paffinity_alone 0";
        }
        # in 1.7 on, looks like we should be using another flag
        elsif( &my_compare_version( $version, "<", "1.8" ) ){
            $val = " -mca mpi_paffinity_alone 0";
            # will be something with the following flag
            # probably need to put some more logic to set numbers right
            # $prun_args .= " hwloc_base_slot_list";
        }
        # darwin uses these
        # maybe this means on 1.8 or higher???
        else{
            $val = " -mca btl self,sm -bind-to none";
        }
        $cmd_pre .= "echo 'RJ_OUTPUT: WARNING: [RJ_L_FIX_ENV_PRUN=$val]'; ";
        $cmd_pre .= &rj_shell_var_set("RJ_L_FIX_ENV_PRUN", $val, ";", $shell);
    }

    # remove extra stuff and return
    while( $cmd_pre =~ /;\s*;/ ){
        $cmd_pre =~ s/;\s*;/;/g;
    }
    $cmd_pre =~ s&\\+\s*$&&;
    return( $cmd_pre );
}

# utility for rj_fix_env
sub rj_fix_env_mod{
    my(
        $shell,
        $load_module,
        $load_message,
        $load_module_default,
        $unload_module,
        $module_use,
        $version_loaded,
        $version_max,
        $version_min,
        $try_ref,
        $changed_ref,
        ) = @_;
    my(
        $cmd,
        $cmd_pre,
        );

    $cmd = "";
    $cmd_pre = "";
    # if loaded reset to max/min (max first so it can be max over choosing min)
    if( $load_module eq "" && $version_max ne "" ){
        if( $version_loaded eq "" ||
            &my_compare_version( $version_loaded, ">", $version_max ) ){
            $load_module = $version_max;
            $load_message = "module version_loaded [$version_loaded] > version_max [$version_max] (or unset)";
        }
    }
    if( $load_module eq "" && $version_min ne "" ){
        if( $version_loaded eq "" ||
            &my_compare_version( $version_loaded, "<", $version_min ) ){
            $load_module = $version_min;
            $load_message = "module version_loaded [$version_loaded] < version_min [$version_min] (or unset)";
        }
    }

    # if the changed value is set to non-blank, then set up to reload
    if( $load_module eq "" && $$changed_ref ne "" && $version_loaded ne "" ){
        $load_module = $version_loaded;
        $load_message = "loaded module [$$changed_ref] so reloading [$version_loaded]";
    }
    
    # if still not defined, just have it be whatever the load_module_default is
    # should not hit here - but keep to ensure default if above changes
    if( $load_module eq "" && $version_loaded eq "" ){
        $load_module = $load_module_default;
        $load_message = "loading default module [$load_module_default]";
    }
    
    # command
    if( $load_module ne "" ){
        $$try_ref++;
        $$changed_ref = $load_module;
        # need to unload other modules in case wrong one loaded
        $cmd = "";
        # always load friendly-testing module (needed for some compilers)
        $cmd .= &rj_shell_module_command( "load friendly-testing", $shell );
        if( defined($module_use) && $module_use =~ /\S/ ){
            $cmd .= &rj_shell_module_command( "use $module_use", $shell );
        }
        if( $unload_module =~ /\S/ ){
            $cmd .= &rj_shell_module_command( "unload $unload_module", $shell );
        }
        $cmd .= &rj_shell_module_command( "load $load_module", $shell );
        $cmd_pre  = "echo 'RJ_OUTPUT: WARNING: $$try_ref $0 $load_message' ; echo 'RJ_OUTPUT: WARNING: $$try_ref $0 trying: $cmd' ; ";
        $cmd_pre .= "$cmd ; ";
    }

    return( $cmd_pre );

}

#############################################################################

# TERM is 15 - ctrl-c
# Do not want to catch this one since you want to allow ctrl-c
# for interactive jobs (maybe always)
#$SIG{TERM} = 'handler_resume';
sub handler_resume{
    my( $sig ) = @_;
    my( $ierr );
    $ierr = 0;
    &print_error( "Caught signal [$sig]",
                  "Resuming...",
                  $ierr );
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
       @args, # arguments
       $ierr, # error ret val
       $num_args, # number of arguments
       $op,
       $opt, # current option
       @vals, # array of values
       $val, # value for current option
       $var,
      );
    $ierr = 0;
    @args = @{$argv_ref};
    if( ! defined($$cmd_ref{conds}) ){
        $$cmd_ref{conds} = "";
    }
    # default help to false
    undef( $$cmd_ref{h} );
    #....................
    #...parse the args...
    #....................
    $num_args = $#args + 1;
    while( @args )
      {
        $opt = shift( @args );
        # skip empty arg
        if( $opt !~ /\S/ ){
        }
        #..........
        #...help...
        #..........
        elsif( $opt =~ /^-+(h(elp)?)$/i ) {
            $opt = "h";
            $$cmd_ref{$opt} = "true";
            # if explicitly set
            $$cmd_ref{"${opt}_set"} = "true";
        }
        #.....................
        #...--conds <conds>...
        #.....................
        elsif( $opt =~ /^--(conds)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 1;
                &print_error( "Value needed for option [--$opt].",
                              $ierr );
                exit( $ierr );
              }
            $val = shift( @args );
            if( ! defined($$cmd_ref{$opt}) ){
                $$cmd_ref{$opt} .= "";
            }
            else{
                $$cmd_ref{$opt} .= ",";
            }
            $$cmd_ref{$opt} .= "$val";
            $$cmd_ref{$opt} =~ s/\s+//;
          }
        #....................
        #...--<opt> <val>...
        #....................
        elsif( $opt =~ /^--(batchid|clean_mpi|convert|depend|dir|id|launchtype|mail|multi|tag|shell|var|wrapper_file|wrapper_out|wrapper_wait)$/ ) {
            $opt = $1;
            if( ! @args ) {
                $ierr = 1;
                &print_error( "Value needed for option [--$opt].",
                              $ierr );
                exit( $ierr );
            }
            $val = shift( @args );
            if( $opt eq "var" ){
                if( $val =~ /^\s*(\w*)\s*([\+\:\?]?=)\s*(.*?)\s*$/ ){
                    $var = $1;
                    $op  = $2;
                    $val = $3;
                    push( @{$$cmd_ref{$opt}{vars}}, $var );
                    push( @{$$cmd_ref{$opt}{ops}},  $op  );
                    push( @{$$cmd_ref{$opt}{vals}}, $val );
                }
                else{
                    $ierr = 1;
                    &print_error( "Format must be:  --var <VAR>=<VAL>",
                                  "  opt = [$opt]",
                                  "  val = [$val]",
                                  $ierr );
                    exit( $ierr );
                }
            }
            elsif( $opt eq "clean_mpi" ){
                if( $val !~ /^(no|yes)$/ ){
                    $ierr = 1;
                    &print_error( "Format must be:  --$opt <no|yes>",
                                  "Found:           --$opt $val",
                                  $ierr );
                    exit( $ierr );
                }
                $$cmd_ref{$opt} = $val;
            }
            else{
                $$cmd_ref{$opt} = $val;
            }
        }

        # --opt <val>: pushed onto array
        elsif( $opt =~ /^--(env)$/ ){
            $opt = $1;
            if( ! @args ) {
                $ierr = 1;
                &print_error( "Value needed for option [--$opt].",
                              $ierr );
                exit( $ierr );
            }
            $val = shift( @args );
            @vals = split(/,/, $val );
            push( @{$$cmd_ref{$opt}}, @vals );
        }

        #...............
        #...-<scalar>...
        #...............
        elsif( $opt =~ /^-(i|v)$/ ) {
            $opt = $1;
            $$cmd_ref{$opt} = "true";
        }
        #... --fix_env
        elsif( $opt =~ /^--(fix_env)$/ ){
            $opt = $1;
            @{$$cmd_ref{$opt}} = @args;
            @args = ();
        }
        #................
        #...--<scalar>...
        #................
        elsif( $opt =~ /^--(debug|info|nodepend|norun_job|wait)$/i ) {
            $opt = $1;
            $$cmd_ref{$opt} = "true";
        }
        # short --var <OPT>=<val>
        elsif( $opt =~ /^--(batch|batch_args|chain|debugger|exec_args|numpe|opt|pname|ppn|prun_args|serial|time|tpp)$/i ) {
            $opt = $1;
            if( ! @args ) {
                $ierr = 1;
                &print_error( "Value needed for option [--$opt].",
                              $ierr );
                exit( $ierr );
            }
            $val = shift( @args );
            # debugger implies i (need here to allow no $RJ_FILE_DATAFILE file)
            if( $opt eq "debugger" && $val =~ /\S/ ){
                $$cmd_ref{i} = "true";
                # for now, totalview does not handle mapping arguments
                unshift(@args, "no_mapping");
                unshift(@args, "--opt");
            }
            # stuff into --var <OPT>=<val>
            ($var = $opt) =~ tr/a-z/A-Z/;
            if( $opt eq "opt" ){
                @vals = split(/,/,$val);
            }
            else{
                @vals = ($val);
            }
            $i = 0;
            foreach $val ( @vals ){
                if( $i == 0 ){
            $op  = "=";
                }
                else{
                    $op  = "+=";
                }
            $opt = "var";
            push( @{$$cmd_ref{$opt}{vars}}, $var );
            push( @{$$cmd_ref{$opt}{ops}},  $op  );
            push( @{$$cmd_ref{$opt}{vals}}, $val );
                $i++;
            }
        }
        # go option
        elsif( $opt =~ /^go$/i ) {
            $$cmd_ref{$opt} = "true";
        }
        # other args
        else {
            # if --, all rest are other args
            if( $opt eq "--" ){
                if( ! @args ) {
                    $ierr = 1;
                    &print_error( "Arguments expected after [$opt].",
                                  $ierr );
                    exit( $ierr );
                }
                $opt = shift( @args );
            }
            # first is EXEC, all others get added to EXEC_ARGS
            elsif( $opt =~ /^-/ ){
                $ierr = 1;
                &print_error( "Unrecognized argument that starts with '-' [$opt]",
                              "If you are defining EXEC/EXEC_ARGS use:",
                              "  '--' <exec> <exec args>",
                              "  '--var EXEC=<exec> --var EXEC_ARGS=<exec args>",
                              $ierr );
                exit( $ierr );
            }
            $var = "EXEC";
            $op  = "=";
            $val = $opt;
            $opt = "var";
            push( @{$$cmd_ref{$opt}{vars}}, $var );
            push( @{$$cmd_ref{$opt}{ops}},  $op  );
            push( @{$$cmd_ref{$opt}{vals}}, $val );
            if( @args ){
                $var = "EXEC_ARGS";
                $op  = "=";
                $val = join(" ", @args );
                push( @{$$cmd_ref{$opt}{vars}}, $var );
                push( @{$$cmd_ref{$opt}{ops}},  $op  );
                push( @{$$cmd_ref{$opt}{vals}}, $val );
                @args = ();
            }
        }
      }

    # environment variables
    foreach $var ( sort keys %ENV ) {
        if( $var =~ /^RJ_E(P)?_(\S+)/ ){
            $val = $ENV{$var};
            $var = $2;
            if( defined($1) ){
                $op = "+=";
            }
            else{
                $op = "=";
            }

            # put into list
            $opt = "var_env";
            push( @{$$cmd_ref{$opt}{vars}}, $var );
            push( @{$$cmd_ref{$opt}{ops}},  $op  );
            push( @{$$cmd_ref{$opt}{vals}}, $val );
        }
    }

    # see if exec was defined anywhere
    foreach $var ( @{$$cmd_ref{var}{vars}} ){
        if( $var eq "EXEC" ){
            $$cmd_ref{exec_defined} = "";
            last;
        }
    }
    # if no args
    if( $num_args == 0 ){
        $$cmd_ref{h} = "true";
    }
    $$cmd_ref{conds} =~ s/^,//;
    $$cmd_ref{conds} =~ s/,$//;

    # default
    if( ! defined($$cmd_ref{clean_mpi}) ){
        $$cmd_ref{clean_mpi} = ""; # neither no|yes
    }

    # if tag, must be correct format YYYYMMDDHHmmss
    if( defined($$cmd_ref{tag}) &&
        $$cmd_ref{tag} !~ /^\d{4}\d{2}\d{2}\d{2}\d{2}\d{2}$/ ) {
        $ierr = 1;
        &print_error( "Invalid tag [$$cmd_ref{tag}]",
                      "Required format [YYYYMMDDHHmmss]\n",
                      $ierr );
        exit( $ierr );
    }

    return( $ierr );
}
