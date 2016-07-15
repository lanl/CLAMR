eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

# NOTE: You can use FindBin in the script. The modules will automatically have access to $FindBin::Bin.
use FindBin qw($Bin);
use lib ("$Bin", "$Bin/lib", "$Bin/../lib");

use my_utils qw (
                 array_to_range
                 conv_time
                 get_id_num
                 get_sysinfo
                 my_notdir
                 my_stat
                 print_error
                 print_perl_obj
                 run_command
                 status_bar
                 which_exec
                 get_pname
                );

use read_output_files qw (
                          parse_output_file
                          parse_output_file_finish
                         );

use Time::Local;

$RJ_DIR                 = "rj_adir";
$RJ_FILE_ID             = "$RJ_DIR/rj_id";
$RJ_FILE_ID_NEXT        = "$RJ_DIR/rj_id_next";
#$RJ_FILE_CMD_OUT        = "rj_cmd_out";
#$RJ_FILE_BATCH_OUT_BASE = "rj_batch_out";

# push some basic PATHs
$ENV{PATH} .= ":/opt/MOAB/bin";

# var to hold list of dead directories so that they are not checked
# init so can be added to.
$DEAD_DIRS = "invalid dir name here";

# number of times to retry a command (scontrol, sinfo, mdiag, ... )
$RETRY_NUM = 5;

# special state names
$S_RUNNING  = "Running";
$S_DEPEND   = "Depend";
$S_BLOCKED  = "Blocked";
$S_ELIGIBLE = "Eligible";
$S_INT      = "Int";
$S_VIOLATES = "violates";
$S_COMPLETED = "Completed";
require run_status_extras_stubs;
&get_sysinfo( \%{$cmd{sys_info}} );


&parse_args( \@ARGV, \%cmd );
$help_extras = &run_status_help_extras();
#..............................
#...must have at least 1 arg...
#..............................
if ( defined( $cmd{h} ) )
  {
    print <<"EOF";
#............................................................................
#...Name
#...====
#... run_status.pl
#...   Get info on jobs you have running and have run.
#...   Different formats are used depending upon what options are set
#...    (eg, if "-u all" is specified, priority, procs, account information is
#...    given).  If you want to change how things are printed, contact the author.
#...
#...   Sorting: $S_RUNNING by JOBID, all other by PRIO.
#...
#...Usage
#...=====
#... run_status.pl
#...             [--batch <yes=default,no,only>]
#...             [--cancel]
#...             [-d <directory glob>]
#...             [-f <format option>]
#...             [-h|--help]
#...             [--dirsid <directory glob>]
#...             [-l]
#...             [-s]
#...             [--state] <comma separated STATE words>
#...             [-u <user>]
#...
#...    [--batch <yes=default,no,only>]
#...      yes:  Get info about jobs from batch system as well as local pids.
#...      no:   Only non-batch jobs.
#...      only: Only batch jobs.
#...    [--cancel]
#...      Cancel the specified jobs listed.
#...    [-d <directory glob>]
#...      Prune for jobs running in the directory list.
#...      Treated as a glob (put in quotes to remove shell expansion)
#...        -d "foo*bar[123] hmm1 hmm2"
#...    [--dirsid <directory glob>]
#...      Look at the list of directories for run_job.pl files that contain
#...      jobids to help match running jobs with their run directories.
#...      Does not prune the jobids like the "-d" option does.
#...    [-f <output format option>]
#...      1: short with some additional info (mtime, time, dt)
#...         default if getting info about self
#...      2: column output (easier for parsing)
#...      3: 1 + directory path
#...      4: short
#...         default if getting info about others
#...      5: 4 + directory path
#...    [-j <jobid>]
#...      Prune for the job with the jobid.
#...      The system tools eventually called do now allow for regexps/globs.
#...    [-l]
#...      Parse the longer output files (-output, prun.out, log, ... )
#...      Without the "-l" option, only the -status file is read (which
#...      takes very little time since only one cycle of info is read).
#...      Depending on the size of your -output file, this could take
#...      substantially longer (maybe 1sec/5MB).
#...      The long output format (-f 2) is automatically set with -l.
#...    [-n <names>]
#...      Prune for jobs with the name given.
#...      The name can be a perl regular expression (use quotes: "rj_launch_70497_[\\d+]")
#...    [-s]
#...      Get the estimated start times.  This command can take a couple
#...      of seconds so that is why it is not done automatically.
#...    [--state] <comma separated STATE words>
#...      Prune on the state (as returned by run_status.pl in the "STATE" info
#...      field).
#...      If the state starts with a "-", that means "not".
#...    [-u <user>]
#...      Prune for jobs for a particular user.
#...      If you specify "all" as the user, information about all users is gotten.
#...      Note: you will only get limited info due to access restrictions.
#...
#...Field Descriptions:
#...   ACCT:             Account
#...   cc/s/p:           cells per second per process
#...   dir:              working directory
#...   dt:               simulation timestep
#...   Flags:            Special flags defining the job (values may changed)
#...   JOBID:            batch system id
#...   lastcycle:        last cycle run
#...   mtime:            last modification time of any file in run directory
#...   mtime_b:          like mtime, but just the first field (_b = brief)
#...   NAME:             batch job name
#...   ncell:            number of cells
#...   ncmax:            (-l) maximum cycle
#...   NODES:            number of nodes
#...   NOTE:             why the job is waiting
#...                       Cof=Child of, Pof=Parent of, resources, $S_INT=Interactive, ...
#...   NOTE_LONG:        (-l) Longer explanation of why job is waiting (if available)
#...   pname:            problem name (best guess)
#...   PARTITION:        Partition of the cluster
#...   PRIO:             Priority
#...   PROCS:            Number of processes
#...   remain:           H:M:S remaining time left for job
#...   Rem/Strt:         H:M:S Remaining time for running jobs or estimated start time
#...   secs/cycle:       (-l) seconds per cycle
#...   start:            H:M:S (negative for how long the job has been running)
#...   STATE:            $S_RUNNING, $S_DEPEND, $S_BLOCKED, $S_ELIGIBLE, ... .
#...   sumwallhr:        sum of the wall hours effectively spent on calculation
#...   time:             simulation time
#...   time/day:         (-l) at the current rate, estimated simulation time
#...                     per 24 hours of wall time
#...   tmax:             (-l) simulation maximum time
#...   USER:             User
#...   wall_time_needed: (-l) at the current rate, estimated wall time needed
#...                     to hit tmax.
$help_extras#...
#...
#...Examples
#...========
#... 1) run_status.pl -l
#...    Get long info of all your jobs in the batch system.
#... 2) run_status.pl -d /scratch3/lmdm/runs/*
#...    Get info about all your runs in the batch system as well as
#...    runs in the specified directories.
#...
#...See Also
#...========
#... run_job.pl - problem submission tool
#............................................................................
EOF
  exit;
  }
#...vars...
my(
   $ierr, # error return value
  );
#...init...
$ierr = 0;
$| = 1;
#...get info about jobs from batch system
&get_job_info( \%job_info, \%cmd );
$JOB_INFO_REF = \%job_info;

# get additional info
if ( $cmd{f} eq "1" || $cmd{f} eq "2" || $cmd{f} eq "3" ){
    $i = 0;
    $num = keys %job_info;
    print "Processing $num jobs:\n";
    $time_a = time();
    foreach $jobid ( sort by_JOBID_NUM keys %job_info ) {
        $i++;
        &status_bar( $i, $num );
        # skip children - parent will take care of it
        if( defined($job_info{$jobid}{parent}) ){
            next;
        }
        # get info from files in dir
        # if doing print formats that use that info
        $dir = $job_info{$jobid}{dir};
        if( &my_stat($dir) == 0 &&
            defined($dir) &&
            ! defined($detailed{$dir}) &&
            -x $dir &&
            -r $dir ){
            $detailed{$dir} = "";
            &get_file_info_crestone( \%cmd, \%{$job_info{$jobid}} );
        }
    }
    $time_b = time();
    print "\n";
    printf( "  Time: %.2f minutes\n", ($time_b - $time_a)/60.0 );
}

# if canceling jobs, do so and leave
if( defined( $cmd{cancel} ) ){
    &rs_cancel( \%job_info );
    exit( $ierr );
}

#...print info about the runs
if( defined($cmd{u}) ){
    print "\nUser: $cmd{u}\n\n";
}
# create orderings of jobids
foreach $jobid ( sort by_JOBID_NUM keys %job_info ) {
  push( @{$jobids_state{$job_info{$jobid}{STATE}}}, $jobid );
}
# $S_RUNNING first - sorted by JOBID_NUM
$state = $S_RUNNING;
if( defined($jobids_state{$state}) ){
  push( @jobids_state_order, @{$jobids_state{$state}} );
}
$state = $S_DEPEND;
if( defined($jobids_state{$state}) ){
  push( @jobids_state_order, @{$jobids_state{$state}} );
}
$state = $S_BLOCKED;
if( defined($jobids_state{$state}) ){
  push( @jobids_state_order, @{$jobids_state{$state}} );
}
# Everything else afterwards sorted by priority
foreach $state ( sort keys %jobids_state ){
  if( $state eq $S_RUNNING || $state eq $S_DEPEND || $state eq $S_BLOCKED || $state eq $S_ELIGIBLE ){
    next;
  }
  @jobids_sorted = sort by_PRIO @{$jobids_state{$state}};
  push( @jobids_state_order, @jobids_sorted );
}
# eligible last
$state = $S_ELIGIBLE;
if( defined($jobids_state{$state}) ){
  @jobids_sorted = sort by_PRIO @{$jobids_state{$state}};
  push( @jobids_state_order, @jobids_sorted );
}
# short output
@jobids_use = @jobids_state_order;
# sum NODES
foreach $jobid ( @jobids_use ){
    $state = $job_info{$jobid}{STATE};
    if( ! defined($nodes_total{$state}) ){
        $nodes_total{$state} = 0;
    }
    $nodes_total{$state} += $job_info{$jobid}{NODES};
}
if( $#jobids_use >= 0 && ($cmd{f} == 1 || $cmd{f} == 3 || $cmd{f} == 4 || $cmd{f} == 5 )){
    if( $cmd{f} == 1 || $cmd{f} == 3 ){
        @fields  = ("JOBID", "STATE", "NOTE",  "ACCT", "Rem/Strt", "USER", "PRIO", "NODES", "mtime_b", "time", "dt", "NAME" );
    }
    else{
        @fields  = ("JOBID", "STATE", "NOTE",  "ACCT", "Rem/Strt", "USER", "PRIO", "NODES",                          "NAME" );
    }
    @fields_und = @fields;
    grep(s/./-/g, @fields_und);
    undef( @lenmax );
    foreach $field ( @fields ){
        push( @lenmax, length( $field ) );
    }
    foreach $jobid ( @jobids_use ){
        $field_num = 0;
        foreach $field ( @fields ){
            $value = $job_info{$jobid}{$field};
            if( ! defined($value) ){
                $value = "-";
            }
            $lenfield = length( $value );
            if( $lenfield > $lenmax[$field_num] ){
                $lenmax[$field_num] = $lenfield;
            }
            $field_num++;
        }
    }
    undef( @formats );
    $formats_string = "";
    for( $field_num = 0; $field_num <= $#fields; $field_num++ ){
        # left justify fields
        $justify = "-";
        if( $fields[$field_num] eq "NAME" ){
            $width = "";
        }
        else{
            $width = $lenmax[$field_num];
        }
        push( @formats, "%${justify}".$width."s" );;
    }
    $formats_string = join( " ", @formats);
    # jobs
    $state = "";
    foreach $jobid ( @jobids_use ) {
        $state_new = $job_info{$jobid}{STATE};
        if( $state_new ne $state ){
            printf( "\n\n $formats_string\n", @fields );
            printf( " $formats_string\n", @fields_und );
            $state = $state_new;
        }
        # skip children - parent will take care of it
        #if ( defined($job_info{$jobid}{parent}) ) {
        #  next;
        #}
        $i = 0;
        foreach $field ( @fields ) {
            $value = $job_info{$jobid}{$field};
            if ( ! defined( $value ) ) {
                $value = "-";
            }
            printf( " ${formats[$i]}", $value );
            $i++;
        }
        print "\n";
        if( $cmd{f} == 3 || $cmd{f} == 5 ){
            $value = $job_info{$jobid}{dir};
            if ( ! defined( $value ) ) {
                $value = "-";
            }
            printf( " $formats[0] %s\n", "", $value );
        }
    }
    print "\n";
}
# full output
if ( $#jobids_use >= 0 && ( $cmd{f} == 2 ) ) {
  # jobs
  foreach $jobid ( @jobids_use ) {
    # skip children - parent will take care of it
    #if ( defined($job_info{$jobid}{parent}) ) {
    #  next;
    #}
    # meta
    @fields = ( 
                "JOBID",
                "STATE",
                "ACCT",
                "NOTE",
                "NOTE_LONG",
                "Flags",
                "start",
                "remain",
                "remain_s",
                "Rem/Strt",
                "dur",
                "USER",
                "PRIO",
                "NODES",
                "NODELIST",
                "PROCS",
                "mtime",
                "time",
                "dt",
                "NAME",
                "PARTITION",
                "pname",
                "dir",
                "lastcycle",
                "ncell",
                "cc/s/p",
                "sumwallhr",
                "tmax",
                "secs/cycle",
                "secs/dmp_write",
                "ncmax",
                "time/day",
                "wall_time_needed"
                );
    foreach $field ( @fields ){
        $val = $job_info{$jobid}{$field};
        if( defined( $val ) ){
            printf( "%20s = %s\n", $field, $val );
        }
    }
    # extra print lines
    &print_run_vals_extras( \%{$job_info{$jobid}} );
    print "\n";
  }
}
# summary node info
print "Node use per state:\n";
$len = 0;
foreach $state ( keys %nodes_total ) {
    if( length($state) > $len ){
        $len = length($state);
    }
}
foreach $state ( sort keys %nodes_total ) {
    if( $state eq $S_RUNNING ){
        next;
    }
    printf("  %${len}s: %9g\n", $state, $nodes_total{$state} );
}
$state = $S_RUNNING;
if( defined($nodes_total{$state}) ){
    printf("  %${len}s: %9g\n", $state, $nodes_total{$state} );
}
print "\nNode info total (PPN=$cmd{sys_info}{L_PPN}):\n";
$len = 0;
foreach $key ( keys %{$cmd{sys_info}{nodes}} ) {
    if( length($key) > $len ){
        $len = length($key);
    }
}
foreach $key ( keys %{$cmd{sys_info}{nodes}} ) {
    printf( "  %${len}s: %9s\n", $key, $cmd{sys_info}{nodes}{$key} );
}
exit( $ierr );
0;

###################################################################################
###################################################################################
###################################################################################

sub rs_cancel{
    my(
        $job_info_ref,
        ) = @_;
    my(
        $cancel,
        $cmd,
        $i,
        $ierr,
        $jobid,
        @jobids,
        $kill,
        $num,
        $output,
        $pgid,
        $pid,
        );

    # find the cancel command
    $cancel = "";
    if( $cancel eq "" ){
        $cancel = &which_exec( "canceljob", QUIET=>"" );
    }
    if( $cancel eq "" ){
        $cancel = &which_exec( "mjobctl", QUIET=>"" );
        if( $cancel ne "" ){
            $cancel = "$cancel -c";
        }
    }
    if( $cancel eq "" ){
        $cancel = &which_exec( "scancel", QUIET=>"" );
        if( $cancel ne "" ){
            $cancel = "$cancel";
        }
    }

    @jobids = reverse sort rj_numerically_id_num keys %{$job_info_ref};

    # cancel each jobid
    $num = $#jobids + 1;
    $i = 0;
    # /usr/bin/kill is prefered since seems to have correct args
    $kill = "/usr/bin/kill";
    if( ! -e $kill ){
        $kill = &which_exec("kill");
    }
    print "Canceling $num jobs\n";
    foreach $jobid ( @jobids ){
        $i++;
        if( $num > 20 ){
            &status_bar( $i, $num )
        }
        else{
            print "  $jobid\n";
        }
        $cmd = "";

        # skip if state invalid (job completed, does not exist)
        if( $$job_info_ref{$jobid}{STATE} eq $S_COMPLETED ||
            $$job_info_ref{$jobid}{STATE} eq "-" ){
                $cmd = "echo 'Skipping jobid [$jobid] state [$$job_info_ref{$jobid}{STATE}]'";
        }

        # if a PID, get process group to kill
        elsif( $jobid =~ /^(\d+):PID$/ ){
            $pid = $1;
            $output = `ps -w -w -eo "pid pgid user args"`;
            if( $output =~ /^\s*$pid\s+(\d+)\s+.*run_job\.pl/m ){
                $pgid = $1;
                # on macs, does not like "--" but handles "-<pid>" - so warning but works
                $cmd = "$kill -15 -- -$pgid";
            }
            else{
                $cmd = "echo '$jobid is not associated with run_job.pl (no longer exists?)'";
            }
        }

        # otherwise, it is a jobid
        else{
            if( $cancel eq "" ){
                $cmd = "echo 'No cancel command for canceling batch $jobid'";
            }
            else{
                $cmd = "$cancel $jobid";
            }
        }
        if( $cmd =~ /\S/ ){
            $output = &run_command( COMMAND=>$cmd, STDOUT=>"" );
        }
    }

    return( $ierr );
}

sub rj_numerically_id_num { &get_id_num($a) <=> &get_id_num($b); }

###################################################################################

sub by_JOBID_NUM{ $$JOB_INFO_REF{$a}{JOBID_NUM} <=> $$JOB_INFO_REF{$b}{JOBID_NUM} }
sub by_PRIO{      $$JOB_INFO_REF{$b}{PRIO}      <=> $$JOB_INFO_REF{$a}{PRIO} }

###################################################################################

#............................................................................
#...Name
#...====
#... get_file_info_crestone
#...
#...Purpose
#...=======
#... Fills out some job_info stuff
#...    cc/s/p
#...    cycle_min
#...    dt
#...    lastcycle
#...    mtime
#...    mtime_b
#...    ncell
#...    ncmax
#...    pname
#...    secs
#...    secs/cycle
#...    sumwallhr
#...    time
#...    time/day
#...    tmax
#...    wall_time_needed
#...
#...Arguments
#...=========
#... job_info     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#... job_info     Intent: inout
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_file_info_crestone{
    my(
       $cmd_ref,
       $job_info_ref,
       ) = @_;
    my(
       %cmd_poff,
       $command,
       $cycle_max,
       $cycle_min,
       %data,
       $dead_dir,
       $dim,
       $dir,
       $dir_output,
       @dirs_output,
       $dt,
       @fields,
       $file,
       @files,
       @files_dir,
       $first_cycle,
       $found,
       %found_fields,
       @lines,
       %mtime,
       $ncmax,
       $nummat,
       $output,
       $pname,
       $secs,
       $secs_cycle,
       %stat,
       $status,
       $time,
       $tmax,
       $wall_time_needed,
       
       );

    #...set cmd option sent to parse_output_file_finish
    $cmd_poff{no_tshift} = "";
    $dir = $$job_info_ref{dir};

    # return if this directory has been a dead directory
    if( $dir =~ m&^${DEAD_DIRS}& ){
        return;
    }

    # mtime
    # Do first and if directory seems to be dead, add it to list and return
    undef( $found );
    undef( $file );
    # just do ls first and use that file
    # had looked for other files, but now will do this
    if( ! defined( $found ) ){
        $output = "";
        # -lta : first is latest one including dotfiles
        # another option....but not available on macs...sigh...
        # find -L . -maxdepth 1 -printf "%T@ %p\n" | sort -n
        if( $$cmd_ref{sys_info}{L_EXEC_FIND} eq "gnu" ){
            $command = "find -L $dir -maxdepth 1 -printf '\%T\@ %p\n' | sort -nr | head -1";
        }
        else{
            $command = "ls -1ta $dir 2>/dev/null";
        }
        $output = &run_command(COMMAND=>$command, TIMEOUT=>"3s", STATUS=>\$status);

        # add to list of DEAD_DIRS and just return
        # timeout sets status bug ls status is ??? if dir not there
        # man page says it is 2, but it is not after perl gets it.
        # could could be false positive and might wipe out dir trees
        # that should not be wiped out...
        if( $status != 0 && $output !~ /^\S+\n/ ){
            # kills whole tree minus LOGNAME
            # probably will just be /panfs/scratch8/vol12
            # which is what we want (the volume that is dead)
            ($dead_dir = $dir) =~ s&/$ENV{LOGNAME}/.*&&;
            $dead_dir =~ s/\./\\./g;
            $DEAD_DIRS .= "|$dead_dir";
            return;
        }
        @fields = split( /\n/, $output );
        if( $#fields >= 0 ){
            $fields[0] =~ s/^(\S+)\s(\S+)/$2/;
            if( $fields[0] =~ m&^/& ){
                $file = $fields[0];
            }
            else{
                $file = "$dir/$fields[0]";
            }
        }
        if( defined($file) && -e $file ){
            $found = "true";
        }
    }
    if( defined( $found ) ){
        if( -e $file ){
            &my_stat( $file, \%stat );
            %mtime = &conv_time( SECS=>$stat{mtime_since} );
            $$job_info_ref{mtime} =
                "$stat{mtime_localtime} ($mtime{string})";
            $$job_info_ref{mtime_b} = "$mtime{string_b}";
        }
    }

    # previously, used get_pname to get the real problem name
    # however, this is relatively expensive when doing a lot of
    # jobs, so just assume that the name of the problem is correct
    # in the submission.
    # $pname = &get_pname( $dir );
    $pname = $$job_info_ref{NAME};
    $$job_info_ref{pname} = $pname;

    # lastcycle file
    $file = "$dir/${pname}-lastcycle";
    if( -e $file ) {
        $command = "cat $file";
        $output = &run_command(COMMAND=>$command);
        $output =~ s/^\s*//;
        @fields = split( /\s+/, $output );
        if( $#fields == 5 ) {
            $$job_info_ref{lastcycle} = $fields[0];
        }
        else {
            $$job_info_ref{lastcycle} = "-";
        }
    }

    # parse some level of output files
    # In order of preference:
    #    <run_name>-output (use first since it has input file vars)
    #    prun.out
    #    *.log.* (the run scripts cleaned up)
    if( defined($$cmd_ref{l}) ) {
        undef( $file );
        if( ! defined( $file ) ) {
            $file = "$dir/${pname}-output";
            if( ! -e $file ) {
                undef( $file );
            }
        }
        if( ! defined( $file ) ) {
            $file = "$dir/rj_cmd_out";
            if( ! -e $file ) {
                undef( $file );
            }
        }
        if( ! defined( $file ) ) {
            $file = "$dir/prun.out";
            if( ! -e $file ) {
                undef( $file );
            }
        }
        if( ! defined( $file ) ) {
            @dirs_output = ("$dir/outputs", "$dir");
            undef( @files );
            foreach $dir_output ( @dirs_output ){
                if( opendir( DIR, $dir_output ) ) {
                    @files_dir = sort( grep( /-output\./, readdir( DIR ) ) );
                    @files_dir = grep( ! /^gold[_\-]/, @files_dir );
                    @files_dir = grep( s&^&$dir_output/&, @files_dir );
                    push( @files, @files_dir );
                    closedir( DIR );
                }
            }
            @files = grep( ! /^gold[_\-]/, @files );
            if( $#files >= 0 ){
                $file = "$files[-1]";
            }
            if( defined( $file ) && ! -e $file ) {
                undef( $file );
            }
        }
        # parse the file
        if( defined( $file ) ) {
            $command = "cat $file";
            $output = &run_command(COMMAND=>$command);
            @lines = split( /\n/, $output );
            undef( $first_cycle );
            undef( $nummat );
            undef( $dim );
            undef( %data );
            undef( %found_fields );
            $first_cycle = "true";
            &parse_output_file( \@lines,
                                \$first_cycle, \$nummat,
                                \%data, \%found_fields, \$dim );
            &parse_output_file_finish( \%data, \%found_fields, \%cmd_poff );
            $$job_info_ref{ncell} = $data{ncell}[-1];
            $$job_info_ref{"cc/s/p"} = $data{"cc/s/p"}[-1];
            $$job_info_ref{"sumwallhr"} = $data{"sumwallhr"}[-1];
            $time = $data{time}[-1];
            $dt = $data{dt}[-1];
            $cycle_min = $data{"cycle_min"}[0];
            $cycle_max = $#{$data{"time"}};
            $secs = $data{"secs"}[-1];
            if( defined($cycle_min) && defined($cycle_max) && defined($secs) &&
                $cycle_max - $cycle_min != 0 ) {
                $secs_cycle = $secs / ($cycle_max - $cycle_min);
            }
            else{
                $secs_cycle = 0;
            }
            $tmax = $data{"tmax"}[0];
            $ncmax = $data{"ncmax"}[0];
            if( defined($tmax) && defined($time) && defined($dt) &&
                defined($secs_cycle)) {
                $wall_time_needed = (($tmax - $time) / $dt) * $secs_cycle;
                %mtime = &conv_time( SECS=>$wall_time_needed );
                if( $wall_time_needed > 0 ){
                    $wall_time_needed = "$mtime{string} (tmax)";
                }
                else{
                    $wall_time_needed = "-$mtime{string} (tmax)";
                }
            }
            elsif( defined($ncmax) && defined($secs_cycle) &&
                   defined($$job_info_ref{lastcycle}) ) {
                $wall_time_needed =
                    ($ncmax - $$job_info_ref{lastcycle}) * $secs_cycle;
                %mtime = &conv_time( SECS=>$wall_time_needed );
                $wall_time_needed = "$mtime{string} (ncmax)";
            }
            if( $secs_cycle > 0 ) {
                $$job_info_ref{"time/day"} =
                    sprintf( "%7.2e", ((60*60*24)/$secs_cycle) * $dt);
            }
            if( $secs_cycle > 0 ){
                $$job_info_ref{"secs/cycle"} = sprintf( "%.2f", $secs_cycle );
            }
            $$job_info_ref{time} = $time;
            $$job_info_ref{dt} = $dt;
            $$job_info_ref{ncmax} = $ncmax;
            $$job_info_ref{tmax} = $tmax;
            $$job_info_ref{dump_write_num}   = $data{"dmp_write_num"}[0];
            $$job_info_ref{"secs/dmp_write"} = $data{"secs/dmp_write"}[0];
            $$job_info_ref{wall_time_needed} = ${wall_time_needed};
            # do this extra parsing
            &get_run_vals_extras( \%data, \%found_fields, $job_info_ref );
        }
    }
    # process other files
    foreach $file ("$dir/${pname}-status",
                   "$dir/code.stdout" ){
        if( -e $file ) {
            # big file, just get the lines currently parsed
            # if full file should be parsed, do so in the
            # -l section instead.
            if( $file eq "$dir/code.stdout" ){
                $command = "grep Cycle: $file";
            }
            else{
                $command = "cat $file";
            }
            $output = &run_command(COMMAND=>$command);
            @lines = split( /\n/, $output );
            undef( $first_cycle );
            undef( $nummat );
            undef( $dim );
            undef( %data );
            undef( %found_fields );
            $first_cycle = "true";
            &parse_output_file( \@lines,
                                \$first_cycle, \$nummat,
                                \%data, \%found_fields, \$dim, 0 );
            $$job_info_ref{ncell}    = $data{ncell}[-1];
            $$job_info_ref{"cc/s/p"} = $data{"cc/s/p"}[-1];
            $$job_info_ref{time}     = $data{time}[-1];
            $$job_info_ref{dt}       = $data{dt}[-1];
            $$job_info_ref{"sumwallhr"} = $data{"sumwallhr"}[-1];
            # do this extra parsing
            &get_run_vals_extras( \%data, \%found_fields, $job_info_ref );
        }
    }
}

#............................................................................
#...Name
#...====
#... get_job_info
#...
#...Purpose
#...=======
#... Fill out job_info struct
#... (wrapper to get_job_info_moab, get_job_info_slurm,
#...
#...Arguments
#...=========
#... job_info_ref Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#... cmd_ref      Intent: inout
#...              Perl type: pointer to hash
#...              Default: none
#...              The command line args ref.
#...
#...Program Flow
#...============
#... 1) Depending on machine, call underlying routine
#............................................................................
sub get_job_info{
    my(
       $job_info_ref,
       $cmd_ref,
       ) = @_;
    my(
       $batch_system_info,
       $EXEC,
       $output,
       $dir,
       $dir_use,
       %dirs_hash,
       $jobid,
       $jobid_parent,
       $jobid_child,
       @jobids,
       $key,
       $name,
       $name_use,
       $names_regexp,
       $user_use,
       %stat,
       $val,
       $val_regexp,
       $val_regexp_not,
       $val_use,
       );

    # get job infos from pids
    if( $$cmd_ref{batch} ne "only" ){
        &get_job_info_pid( $job_info_ref );
    }

    # dawn has "moab" - but not really, they have their own wrappers.
    # key is if mdiag has a "--format option"
    if( $$cmd_ref{batch} ne "no" ){
        $batch_system_info = "";

        # see if full moab mdiag output
        if( $batch_system_info eq "" ){
        $MDIAG = &which_exec( "mdiag", QUIET=>"" );
            if( $MDIAG =~ /\S/ ){
            $output = &run_command( COMMAND=>$MDIAG, TIMEOUT=>"3s" );
            if( defined($output) && ( $output =~ /--format/ || $output =~ /--xml/ ) ){
                    $batch_system_info = "moab";
            }
            else{
                    $batch_system_info = "slurm";
                }
            }
            }

        # if they have sinfo, slurm
        if( $batch_system_info eq "" ){
            $EXEC = &which_exec( "sinfo", QUIET=>"" );
            if( $EXEC =~ /\S/ ){
                $batch_system_info = "slurm";
            }
        }

        # moab
        if( $batch_system_info eq "moab" ){
            &get_job_info_moab( $job_info_ref, $cmd_ref );
        }
        # slurm
        elsif( $batch_system_info eq "slurm" ){
            &get_job_info_slurm( $job_info_ref, $cmd_ref );
        }
    }

    # now get job info from any dir that has not been processed yet
    &get_job_info_dirs( $job_info_ref, $cmd_ref );

    # prune if given directory
    if( defined( $$cmd_ref{d} ) ){
        foreach $dir( @{$$cmd_ref{d}} ){
            &my_stat( $dir, \%stat );
            $dir_use = $stat{fullpath};
            if( defined($dir_use) ){
                $dirs_hash{$dir_use} = "";
            }
        }
        @jobids = keys %{$job_info_ref};
        foreach $jobid ( @jobids ){
            undef( $dir_use );
            if( defined( $$job_info_ref{$jobid}{dir} ) ){
                $dir_use = $$job_info_ref{$jobid}{dir};
            }
            if( ! defined( $dir_use ) ||
                ! defined( $dirs_hash{$dir_use} ) ){
                delete($$job_info_ref{$jobid});
            }
        }
    }

    # prune now so no unnecessary finishing afterwards
    # prune if given a name
    if( defined( $$cmd_ref{n} ) ){
        $names_regexp = "";
        foreach $name( @{$$cmd_ref{n}} ){
            $names_regexp .= "$name|";
        }
        $names_regexp =~ s/\|$//;
        @jobids = keys %{$job_info_ref};
        foreach $jobid ( @jobids ){
            undef( $name_use );
            if( defined( $$job_info_ref{$jobid}{NAME} ) ){
                $name_use = $$job_info_ref{$jobid}{NAME};
            }
            if( ! defined( $name_use ) ||
                $name_use !~ /^($names_regexp)$/ ){
                delete($$job_info_ref{$jobid});
            }
        }
    }

    # prune if given a jobid
    if( defined( $$cmd_ref{j} ) ){
        @jobids = keys %{$job_info_ref};
        foreach $jobid ( @jobids ){
            if( $jobid ne $$cmd_ref{j} &&
                $jobid ne $$cmd_ref{j_num} ){
                delete($$job_info_ref{$jobid});
            }
        }
    }

    # prune if given a user
    if( defined( $$cmd_ref{u} ) ){
        @jobids = keys %{$job_info_ref};
        foreach $jobid ( @jobids ){
            undef( $user_use );
            if( defined( $$job_info_ref{$jobid}{USER} ) ){
                $user_use = $$job_info_ref{$jobid}{USER};
            }
            if( ! defined( $user_use ) ||
                $user_use !~ /^($$cmd_ref{u})$/ ){
                delete($$job_info_ref{$jobid});
            }
        }
    }

    # fill in some defaults
    foreach $jobid ( keys %$job_info_ref ){
        if( ! defined($$job_info_ref{$jobid}{NAME}) ){
            $$job_info_ref{$jobid}{NAME} = "-";
        }
        if( ! defined($$job_info_ref{$jobid}{STATE}) ){
            $$job_info_ref{$jobid}{STATE} = "-";
        }
        if( ! defined($$job_info_ref{$jobid}{NODES}) ){
            $$job_info_ref{$jobid}{NODES} = 0;
        }
        if( ! defined($$job_info_ref{$jobid}{PRIO}) ){
            $$job_info_ref{$jobid}{PRIO} = 0;
        }
        if( ! defined($$job_info_ref{$jobid}{dir}) ){
            $$job_info_ref{$jobid}{dir} = "-";
        }
        $$job_info_ref{$jobid}{JOBID_NUM} = &get_id_num( $$job_info_ref{$jobid}{JOBID} );
    }

    # prune Completed jobs
    @jobids = keys %{$job_info_ref};
    foreach $jobid ( @jobids ){
        if( $$job_info_ref{$jobid}{STATE} eq $S_COMPLETED ){
            delete($$job_info_ref{$jobid});
        }
    }

    # replace dir with shell expanded vars and fullpath returned from my_stat
    foreach $jobid ( keys %$job_info_ref ){
        if( $$job_info_ref{$jobid}{dir} ne "-" ){
            # can be expensive, so see if need to do this (has a "$" in it)
            if( $$job_info_ref{$jobid}{dir} =~ /\$/ ){
                $$job_info_ref{$jobid}{dir} = `echo "$$job_info_ref{$jobid}{dir}"`;
            }
            $$job_info_ref{$jobid}{dir} =~ s/\s*$//;
            &my_stat( $$job_info_ref{$jobid}{dir}, \%stat );
            if( defined( $stat{fullpath} ) ){
                $$job_info_ref{$jobid}{dir} = $stat{fullpath};
            }
        }
    }

    # parent/child
    # if has a parent
    foreach $jobid ( keys %$job_info_ref ) {
        $jobid_parent = $$job_info_ref{$jobid}{parent};
        # if this job is a child
        if( defined($jobid_parent) ){
            # if parent still exists
            if( defined($$job_info_ref{$jobid_parent}) ) {
                $$job_info_ref{$jobid_parent}{"child"} = $jobid;
            }
            # otherwise not really a child
            else{
                delete( $$job_info_ref{$jobid}{parent} );
            }
        }
        # make new state
        if( defined($$job_info_ref{$jobid}{parent} ) ){
            $$job_info_ref{$jobid}{STATE} = $S_DEPEND;
        }
    }
    # if has a child, just make sure child knows it has a parent
    foreach $jobid ( keys %$job_info_ref ) {
        $jobid_child = $$job_info_ref{$jobid}{child};
        if( defined( $jobid_child ) ){
            # if child exists
            if( defined($$job_info_ref{$jobid_child}) ){
                if( ! defined($$job_info_ref{$jobid_child}{parent}) ){
                    $$job_info_ref{$jobid_child}{parent} = $jobid;
                    # might want to do something with the state...not sure
                }
            }
        }
    }

    # create NOTE
    foreach $jobid ( keys %$job_info_ref ) {
        # put in NOTE Cof and Pof
        if( defined($$job_info_ref{$jobid}{parent} ) ){
            $$job_info_ref{$jobid}{NOTE} .= ":Cof:".$$job_info_ref{$jobid}{parent};
        }
        elsif( defined($$job_info_ref{$jobid}{child} ) ){
            $$job_info_ref{$jobid}{NOTE} .= ":Pof:".$$job_info_ref{$jobid}{child};
        }
    }

    # interactive: add to NOTE
    foreach $jobid ( keys %$job_info_ref ) {
      # interactive: tack on to NOTE
      if( defined($$job_info_ref{$jobid}{Flags}) &&
          $$job_info_ref{$jobid}{Flags} =~ /\binteractive\b/i ){
          $$job_info_ref{$jobid}{NOTE} .= ":$S_INT";
      }
      elsif( defined($$job_info_ref{$jobid}{NAME}) &&
          ( $$job_info_ref{$jobid}{NAME} eq "mxterm" ||
            $$job_info_ref{$jobid}{NAME} eq "bash" ||
            $$job_info_ref{$jobid}{NAME} eq "tcsh" ||
            $$job_info_ref{$jobid}{NAME} eq "csh" ||
            $$job_info_ref{$jobid}{NAME} eq "sh" ||
            $$job_info_ref{$jobid}{NAME} eq "xterm" ||
            $$job_info_ref{$jobid}{NAME} eq "llogin" ) ){
          $$job_info_ref{$jobid}{NOTE} .= ":$S_INT";
      }
    }

    # finish off things
    foreach $jobid ( keys %$job_info_ref ) {

        # clean up fields
        foreach $key ( "NOTE", "STATE" ){
            if( !defined($$job_info_ref{$jobid}{$key}) ){
                next;
        }

            # extra :none:
            $$job_info_ref{$jobid}{$key} =~ s/:*(\bnone\b):*/:/gi;

            # extra :
            $$job_info_ref{$jobid}{$key} =~ s/^:+//;
            $$job_info_ref{$jobid}{$key} =~ s/:+$//;
            $$job_info_ref{$jobid}{$key} =~ s/:+/:/g;

            # undef if nothing
            if( $$job_info_ref{$jobid}{$key} !~ /\S/ ){
                undef($$job_info_ref{$jobid}{$key});
            }
        }

        # in case NOTE_LONG has any returns
        if( defined($$job_info_ref{$jobid}{NOTE_LONG}) ){
            $$job_info_ref{$jobid}{NOTE_LONG} =~ s/\n+/ /g;
        }
        # NODES
        if( !defined($$job_info_ref{$jobid}{NODES}) ){
            if( defined($$job_info_ref{$jobid}{PROCS}) ){
                $$job_info_ref{$jobid}{NODES} =
                    POSIX::ceil($$job_info_ref{$jobid}{PROCS}/$$cmd_ref{sys_info}{L_PPN});
            }
            else{
                $$job_info_ref{$jobid}{NODES} = "-";
            }
        }
    }

    # Rem/Strt
    foreach $jobid ( keys %$job_info_ref ) {
        if( defined($$job_info_ref{$jobid}{"Rem/Strt"}) ){
            next;
        }
        $$job_info_ref{$jobid}{"Rem/Strt"} = $$job_info_ref{$jobid}{"start"};
        if( defined( $$job_info_ref{$jobid}{"remain"}) ){
            $$job_info_ref{$jobid}{"Rem/Strt"} = $$job_info_ref{$jobid}{"remain"};
        }
    }

    # prune if given a state
    # need to be done after finishing of fields because STATE changes
    if( defined( $$cmd_ref{state} ) ){
        $val_regexp = "";
        $val_regexp_not = "";
        foreach $val( @{$$cmd_ref{state}} ){
            if( $val =~ /^-(\S+)/ ){
                $val_regexp_not .= "$1|";
            }
            else{
                $val_regexp .= "$val|";
            }
        }
        $val_regexp =~ s/\|$//;
        $val_regexp_not =~ s/\|$//;
        @jobids = keys %{$job_info_ref};
        foreach $jobid ( @jobids ){
            undef( $val_use );
            if( defined( $$job_info_ref{$jobid}{STATE} ) ){
                $val_use = $$job_info_ref{$jobid}{STATE};
            }
            
            # if not not, delete if not matches
            if( $val_regexp ne "" ){
                if( ! defined( $val_use ) || $val_use !~ /^($val_regexp)$/ ){
                    delete($$job_info_ref{$jobid});
                }
            }

            # if not, delete if matches
            if( $val_regexp_not ne "" ){
                if( defined( $val_use ) && $val_use =~ /^($val_regexp_not)$/ ){
                    delete($$job_info_ref{$jobid});
                }
            }
        }
    }
}

# ==============================================================================

sub get_job_info_dirs{
    my(
        $job_info_ref,
        $cmd_ref,
        ) = @_;

    my(
        $dir,
        $dir_use,
        @dirs,
        $i,
        $jobid,
        %rj_dir_info,
        %stat,
        $user,
        );
    
    if( defined( $$cmd_ref{d} ) ){
        push( @dirs, @{$$cmd_ref{d}} );
    }
    if( defined( $$cmd_ref{dirsid} ) ){
        push( @dirs, @{$$cmd_ref{dirsid}} );
    }

    # look at each directory
    $i = -1;
    foreach $dir ( @dirs ) {
        
        # skip if not dir
        if( ! -d $dir ){
            next;
        }
        
        # get actual path and user
        &my_stat( $dir, \%stat );
        undef( $dir_use );
        if( defined( $stat{fullpath} ) ){
            $dir_use = $stat{fullpath};
            $user = $stat{user};
        }
        
        # if dir is real
        if( defined( $dir_use ) ){
            
            &get_rj_dir_info( $dir_use, \%rj_dir_info );
            
            # if the pid matches a JOBID_NEXT (so a --depend nobatch)
            if( defined($$job_info_ref{$rj_dir_info{JOBID_NEXT}}) ){
                $jobid = $rj_dir_info{JOBID_NEXT};
            }
            # set jobid based on if job was nobatch/batch
            # if job was nobatch job
            elsif( $rj_dir_info{LAUNCHTYPE} eq "nobatch" ){
                $jobid = $rj_dir_info{PID};
            }
            # launch was batch
            else{
                $jobid = $rj_dir_info{BATCHID};
            }

            # default jobid is negative number
            if( $jobid eq "" ){
                $jobid = $i;
                $i--;
            }

            # if not already gotten, fill with some info
            if( ! defined($$job_info_ref{$jobid}) ){
                $$job_info_ref{$jobid}{JOBID}     = $jobid;
                $$job_info_ref{$jobid}{JOBID_NUM} = &get_id_num( $jobid );
                $$job_info_ref{$jobid}{NAME}      = &get_pname( $dir );
                $$job_info_ref{$jobid}{dir}       = $dir_use;
                $$job_info_ref{$jobid}{USER}      = $user;
            }
        }
    }
}

# ==============================================================================

sub get_rj_dir_info{
    my(
        $dir,
        $rj_dir_info_ref,
        ) = @_;

    my(
        @lines,
        $line,
        $out,
        );

    # fill in stuff - will get replaced with correct
    # from get_job_info_* if actually running
    $$rj_dir_info_ref{PID} = "";
    $$rj_dir_info_ref{BATCHID} = "";
    $$rj_dir_info_ref{LAUNCHTYPE} = "";
    $$rj_dir_info_ref{JOBID_NEXT} = "";
    
    # get info from RJ_FILE_*
    if( -r "$dir/$RJ_FILE_ID" ){
        $out = `cat $dir/$RJ_FILE_ID 2> /dev/null`;
        @lines = split(/\n/, $out );
        foreach $line ( @lines ){
            if( $line =~ /^\s*(BATCHID|PID|LAUNCHTYPE)\s*=\s*(\S+)\s*$/ ){
                $$rj_dir_info_ref{$1} = $2;
            }
        }
    }
    if( -r "$dir/$RJ_FILE_ID_NEXT" ){
        $$rj_dir_info_ref{JOBID_NEXT} = `cat $dir/$RJ_FILE_ID_NEXT`;
        $$rj_dir_info_ref{JOBID_NEXT} =~ s/\s+$//;
    }
}

# ==============================================================================

sub get_job_info_pid{
    my(
        $job_info_ref,
        ) = @_;

    my(
        $args,
        $arg_launchtype,
        $arg_dir,
        $Cof,
        $cmd,
        $dir,
        $jobid,
        @lines_rj_id,
        $line,
        @lines,
        $name,
        $out,
        $pid,
        $pgid,
        $ppn,
        @pwuid,
        %rj_dir_info,
        $status,
        $user,
        );

    $cmd = 'ps -w -w -eo "pid pgid user args"';
    $out = &run_command( COMMAND=>$cmd, STATUS=>\$status );
    @lines = split(/\n/, $out);
    foreach $line( @lines ){
        undef( $pid );
        if( $line =~ /^\s*(\d+)\s+(\d+)\s+(\S+)\s+(.*?)\s*$/ ){
            $pid = $1;
            $pgid = $2;
            $user = $3;
            $args = $4;
            undef( $dir );
            undef( $name );
            undef( $jobid );
            undef( $Cof );

            # if user seems to be a uid, translate it to user if possible.
            # had a user where ps was returning uid and not moniker for some reason.
            if( $user =~ /^\d+$/ ){
                @pwuid = getpwuid( $user );
                if( $#pwuid > 0 ){
                    $user = $pwuid[0];
                }
            }

            # this is the process to detect - does the fork
            # this order id done by run_job.pl so stays the same
            if( $args =~ m&run_job.pl\sgo\s
                           --id\s(\S+)\s
                           --launchtype\s(\S+)\s
                           --tag\s(\S+)\s
                           --dir\s(\S+)\s
                           (\s--batchid\s(\S+))?
                          &x ){
                # get info about command line
                $arg_dir   = $4;
                $arg_launchtype = $6;

                # jobid is always this pid
                $jobid = "$pid:PID";
                $$job_info_ref{$jobid}{JOBID}     = $jobid;
                $$job_info_ref{$jobid}{STATE}     = $S_RUNNING;
                $$job_info_ref{$jobid}{USER}      = $user;
                $$job_info_ref{$jobid}{dir}       = $arg_dir;
                $$job_info_ref{$jobid}{NOTE}      = "$pgid:PGID";
                $$job_info_ref{$jobid}{JOBID_NUM} = &get_id_num( $$job_info_ref{$jobid}{JOBID} );
                $$job_info_ref{$jobid}{NAME} = &get_pname( $arg_dir );
                &get_rj_dir_info( $arg_dir, \%rj_dir_info );
                
                # if this is the running PID
                if( $rj_dir_info{PID} eq $jobid ){
                    if($rj_dir_info{JOBID_NEXT} ne ""){
                        $$job_info_ref{$jobid}{child}  = $rj_dir_info{JOBID_NEXT};
                    }
                }
                # if not, see if this is the child
                else{
                    if($rj_dir_info{JOBID_NEXT} eq $jobid){
                        $$job_info_ref{$jobid}{STATE}     = $S_DEPEND;
                        $$job_info_ref{$jobid}{parent}    = $rj_dir_info{PID};
                    }
                }
            }

            # is doing the wait for batch=no
            elsif( $args =~ m&run_job.pl.*\s
                              --dir\s(\S+)\s
                              --wrapper_file\s(\S+)\s
                              &x ){
                $jobid = "$pid:PID";
                $arg_dir   = $1;
                # do now know various thing - if given "-d <dir>" then that will fill
                # in other things
                $$job_info_ref{$jobid}{JOBID}     = $jobid;
                $$job_info_ref{$jobid}{STATE}     = $S_DEPEND;
                $$job_info_ref{$jobid}{USER}      = $user;
                $$job_info_ref{$jobid}{dir}       = $arg_dir;
                $$job_info_ref{$jobid}{NOTE}      = "$pgid:PGID";
                $$job_info_ref{$jobid}{JOBID_NUM} = &get_id_num( $$job_info_ref{$jobid}{JOBID} );
                $$job_info_ref{$jobid}{NAME} = &get_pname( $arg_dir );
            }

            # get info about run
            # todo: run_job.pl could create an info file
            if( defined($jobid) &&
                defined($$job_info_ref{$jobid}{dir}) &&
                -T "$$job_info_ref{$jobid}{dir}/$RJ_FILE_ID" ){
                $out = `cat $$job_info_ref{$jobid}{dir}/$RJ_FILE_ID 2> /dev/null`;
                @lines_rj_id = split(/\n/, $out);
                foreach $lines_rj_id ( @lines_rj_id ){
                    if( $lines_rj_id =~ /^\s*(NODES_exn)\s*=\s*(\S+)\s*$/ ){
                        $$job_info_ref{$jobid}{NODES}      = $2;
                    }
                    elsif( $lines_rj_id =~ /^\s*(NUMPE_max)\s*=\s*(\S+)\s*$/ ){
                        $$job_info_ref{$jobid}{PROCS}      = $2;
                    }
                    elsif( $lines_rj_id =~ /^\s*(PPN_min)\s*=\s*(\S+)\s*$/ ){
                        $ppn = $2;
                    }
                }
                if( defined($$job_info_ref{$jobid}{NODES}) &&
                    defined($$job_info_ref{$jobid}{PROCS}) ){
                    $$job_info_ref{$jobid}{NODES} = $$job_info_ref{$jobid}{PROCS} / $ppn;
                }
            }
        }
    }
}

# ==============================================================================

#............................................................................
#...Name
#...====
#... get_job_info_slurm
#...
#...Purpose
#...=======
#... Fills out some job_info from "moab" commands.  It's moab wrappers
#....that llnl wrote - so not really moab.  Pretty annoying...
#...
#...Arguments
#...=========
#... job_info_ref Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#... cmd_ref      Intent: inout
#...              Perl type: pointer to hash
#...              Default: none
#...              The command line args ref.
#...
#...Program Flow
#...============
#... 1) Run various faux-moab commands to try and fill struct
#............................................................................
sub get_job_info_slurm{
    my(
       $job_info_ref,
       $cmd_ref,
       ) = @_;
    my(
       $active,
       $command,
       $factor,
       $i,
       $idle,
       $num,
       $output,
       $other,
       $s,
       $total,
       );

    # get the number of nodes available 
    $command = 'sinfo -o "%F %C" -h';
    $output = &run_command( COMMAND=>$command, TIMING=>"", TIMEOUT=>"5s" );
    if( $output =~ m&^\s*(\d\S*)/(\d\S*)/(\d\S*)/(\d\S*)&m ){
        $active = $1;
        $idle   = $2;
        $other  = $3;
        $total  = $4;
        $$cmd_ref{sys_info}{nodes}{Active} = $active;
        $$cmd_ref{sys_info}{nodes}{Idle}   = $idle;
        $$cmd_ref{sys_info}{nodes}{Down}   = $other;
        if( $total =~ /(\S+)([a-zA-Z])/i ){
            $total = $1;
            $factor = $2;
            if( $factor =~ /k/i ){
                $factor = 1024;
            }
            elsif( $factor =~ /m/i ){
                $factor = 1024*1024;
            }
            elsif( $factor =~ /g/i ){
                $factor = 1024*1024*1024;
            }
            else{
                $factor = 1; # not right probably
            }
            $total *= $factor;
        }
        $$cmd_ref{sys_info}{nodes}{Total}  = $total;
    }

    # just use scontrol - seems to have all the info that their
    # moab wrapper has and more.
    &get_job_info_scontrol( $cmd_ref, JOB_INFO=>$job_info_ref );

    $num = scalar keys %{$job_info_ref};
    print "    Number of jobs: ", $num, "\n";

    # showstart not needed since data in scontrol
    #if( $num > 0 ){
    #    # showstart
    #    if( defined($$cmd_ref{s}) ){
    #        &get_job_info_showstart( JOB_INFO=>$job_info_ref );
    #    }
    #}
}

#............................................................................
#...Name
#...====
#... get_job_info_moab
#...
#...Purpose
#...=======
#... Call various moab commands to fill out job_info
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#... USER         Intent: in
#...              Perl type: pointer to scalar
#...              Default: none
#...              User you are interested in (or all)
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_moab{
    my(
       $job_info_ref,
       $cmd_ref,
       ) = @_;
    my(
       $command,
       $key,
       $line,
       $num,
       $output,
       );
    # get the number of nodes available 
    if( $MDIAG =~ /\S/ ){
        $command = "$MDIAG -n";
        $output = &run_command( COMMAND=>$command, TIMING=>"", TIMEOUT=>"5s" );
        if( $output =~ /^\s*(Total Nodes\s*:\s*(\d+).*)$/mi ){
            $line = $1;
            $$cmd_ref{sys_info}{nodes}{Total} = $2;
            foreach $key ("Active","Idle","Down"){
                if( $line =~ /$key.*?(\d+)/i ){
                    $$cmd_ref{sys_info}{nodes}{"$key"} = $1;
                }
            }
        }
    }
    &get_job_info_mdiag( $cmd_ref, JOB_INFO=>$job_info_ref, USER=>$$cmd_ref{u}, JOBID_NUM=>$$cmd_ref{j_num} );
    # mshow - no longer needed since mdiag is used
    # &get_job_info_mshow( JOB_INFO=>\%job_info, USER=>$user );
    $num = scalar keys %{$job_info_ref};
    print "    Number of jobs: ", $num, "\n";
    if( $num > 0 ){
        # showstart
        if( defined($$cmd_ref{s}) ){
            &get_job_info_showstart( JOB_INFO=>$job_info_ref );
        }
    }
}

#............................................................................
#...Name
#...====
#... get_job_info_mdiag
#...
#...Purpose
#...=======
#... Fills out the job_info hash from a command.
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_mdiag{
  my( $cmd_ref ) = shift( @_ );
  my %args = (
              JOB_INFO  => undef,
              USER      => undef,
              JOBID_NUM => undef,
              @_,
             );
  my(
     $att,
     $dep,
     $done,
     $ierr,
     $job,
     $job_info_ref,
     $jobid,
     $key,
     $output,
     $command,
     $message,
     $Message,
     $retry,
     $time_info,
     %val,
     );

  $ierr = 0;

  #...invalid args
  foreach $key ( keys %args )
    {
      if( $key !~ /^(JOB_INFO|USER|JOBID_NUM)$/)
        {
          $ierr = 1;
          &print_error( "Invalid argument to get_job_info_mdiag [$key]",
                        $ierr );
          exit( $ierr );
        }
    }

  $job_info_ref = $args{JOB_INFO};
  if( $MDIAG !~ /\S/ ){
      return;
  }
  $command = "$MDIAG --format=xml -j -v";
  $done = "false";
  $retry = 0;
  while( $done eq "false" ){
      $output = &run_command( COMMAND=>$command, TIMING=>"", TIMEOUT=>"5s" );
      $done = "true";
      # should start with <data> - otherwise temporary problem
      if( $output !~ /^\s*<Data>/i ){
          $done = "false";
          print "Temporary info issue.\n[$output]\n";
          $retry++;
          if( $retry < $RETRY_NUM ){
              print "Retrying...\n";
              sleep( 5 );
              next;
          }
          # give up with error returned (so can be checked for)
          else{
              $ierr = 1;
              &print_error( "Repeated tries at command failed:",
                            "  $command",
                            "Giving up...",
                            $ierr );
              exit( $ierr );
          }
      }
  }

  # parse the output
  use XML::LibXML;
  my $parser;
  my $doc;
  my $root;
  my @jobs;
  if( $output ne "" ){
    $parser = new XML::LibXML;
    $doc = $parser->parse_string( $output );
    $root = $doc->getDocumentElement;
    @jobs = $root->getElementsByTagName('job');
  }
  else{
      @jobs = ();
  }
  foreach my $job ( @jobs ){
      my $att = $job->getAttribute('User');
      $jobid = $job->getAttribute('JobID');
      if( ! defined( $att ) || ! defined( $jobid ) ){
          next;
      }
      if( defined( $args{USER} ) && $args{USER} ne $att ){
          next;
      }
      if( defined( $args{JOBID_NUM} ) && $args{JOBID_NUM} ne &get_id_num($jobid) ){
          next;
      }
      undef($$job_info_ref{$jobid});
      $$job_info_ref{$jobid}{JOBID} = $jobid;
      if( defined( $$cmd_ref{debug} ) ){
          print "$jobid\n";
      }
      #$job_ref = $$job_info_ref{$jobid};
      my @attributes = $job->attributes;
      foreach $att ( sort @attributes ){
          my $att_name = $att->name;
          my $att_val  = $att->value;
          if( defined( $$cmd_ref{debug} ) ){
              print "   [$att_name] = [$att_val]\n";
          }
          # status type fields
          if( $att_name eq "QueueStatus" ){
              $$job_info_ref{$jobid}{QueueStatus} = $att_val;
          }
          elsif( $att_name eq "State" ){
              $$job_info_ref{$jobid}{STATE} = $att_val;
          }
          elsif( $att_name eq "User" ){
              $$job_info_ref{$jobid}{USER} = $att_val;
          }
          elsif( $att_name eq "StartTime" ){
              my $time = time() - $att_val;
              $$job_info_ref{$jobid}{start_s} = $time;
              %time_info = &conv_time(SECS=>$time);
              $$job_info_ref{$jobid}{start} = "-$time_info{hms}";
          }
          elsif( $att_name eq "ReqAWDuration" ){
              $$job_info_ref{$jobid}{dur} = $att_val;
              # latest version of moab seems to have the unit in seconds
              # if no other markers, make it seconds
              if( $att_val =~ /^\d+$/ ){
                  $att_val = "${att_val}s";
              }
              %time_info = &conv_time(STRING=>$att_val);
              $$job_info_ref{$jobid}{dur_s} = "$time_info{SECS_TOTAL}";
          }
          elsif( $att_name eq "BlockReason" ){
              $$job_info_ref{$jobid}{BlockReason} = $att_val;
          }
          elsif( $att_name eq "Account" ){
              $$job_info_ref{$jobid}{ACCT} = $att_val;
          }
          elsif( $att_name eq "IWD" ){
              $$job_info_ref{$jobid}{dir} = $att_val;
          }
          elsif( $att_name eq "JobName" ){
              $$job_info_ref{$jobid}{NAME} = $att_val;
          }
          elsif( $att_name eq "StartPriority" ){
              $$job_info_ref{$jobid}{PRIO} = $att_val;
          }
          elsif( $att_name eq "Flags" ){
              $$job_info_ref{$jobid}{Flags} = $att_val;
          }
          elsif( $att_name eq "Hold" ){
              $$job_info_ref{$jobid}{NOTE} .= "Hold=$att_val";
          }
          elsif( $att_name eq "ReqProcs" ){
              $$job_info_ref{$jobid}{PROCS} = $att_val;
          }
          elsif( $att_name eq "DependBlock" ){
              if( $att_val =~ /jobcomplete:(\S+)/ ){
                  # junk stuff at the end for parent/child - seems to match better
                  $dep = $1;
                  $dep =~ s/(\d+)\..*/$1/;
                  $$job_info_ref{$jobid}{"parent"} = $dep;
              }
          }
      }
      # other info
      my @reqs = $job->getElementsByTagName('req');
      if( $#attributes ){
          @attributes = $reqs[0]->attributes;
          foreach $att ( sort @attributes ){
              my $att_name = $att->name;
              my $att_val  = $att->value;
              my @att_vals;
              if( $att_name eq "TCReqMin" ){
                  ( $$job_info_ref{$jobid}{PROCS} = $att_val ) =~ s/\*$//;
              }
              if( $att_name eq "NCReqMin" ){
                  ( $$job_info_ref{$jobid}{NODES} = $att_val ) =~ s/\*$//;
              }
              if( $att_name eq "AllocNodeList" ){
                  $att_val =~ s/:\d+//g;
                  @att_vals = split(/,/,$att_val);
                  $$job_info_ref{$jobid}{NODELIST} = &array_to_range( \@att_vals );
              }
              if( defined( $$cmd_ref{debug} ) ){
                  print "   req [$att_name] = [$att_val]\n";
              }
          }
      }
      my @Messages = $job->getElementsByTagName('Messages');
      $done = "false";
      foreach $Message ( @Messages ){
          my @messages = $Message->getElementsByTagName('message');
          foreach $message ( @messages ){
              my @attributes = $message->attributes;
              foreach $att ( @attributes ) {
                  my $att_name = $att->name;
                  my $att_val  = $att->value;
                  if( $att_name eq "DATA" ){
                      if( $att_val =~ /cannot\s*migrate/i ){
                          $$job_info_ref{$jobid}{NOTE_LONG} = $att_val;
                          $$job_info_ref{$jobid}{STATE}     = $S_BLOCKED;
                          if( $att_val =~ /not\s*in\s*this\s*partition/i ){
                              $$job_info_ref{$jobid}{NOTE} = "wrong_partition";
                          }
                          else{
                              $$job_info_ref{$jobid}{NOTE} = "cannot_migrate";
                          }
                          # break out of loops
                          $done = "true";
                          last;
                      }
                  }
              }
              if( $done eq "true" ){
                  last;
              }
          }
          if( $done eq "true" ){
              last;
          }
      }
      # STATE
      if( ! defined($$job_info_ref{$jobid}{STATE}) ){
          $$job_info_ref{$jobid}{STATE} = "-";
      }
      # make consistent names
      if( $$job_info_ref{$jobid}{STATE} =~ /^running$/i ){
          $$job_info_ref{$jobid}{STATE} = $S_RUNNING;
      }
      # will fill in with other info
      elsif( $$job_info_ref{$jobid}{STATE} =~ /^idle$/i ){
          $$job_info_ref{$jobid}{STATE} = "";
      }
      # Starting:active -> running
      if( $$job_info_ref{$jobid}{STATE} =~ /^starting$/i &&
          $$job_info_ref{$jobid}{QueueStatus} =~ /^active$/ ){
          $$job_info_ref{$jobid}{STATE} = $S_RUNNING;
      }
      if( $$job_info_ref{$jobid}{STATE} ne $S_RUNNING &&
          $$job_info_ref{$jobid}{STATE} ne $S_BLOCKED ){
          if( defined($$job_info_ref{$jobid}{QueueStatus}) ){
              $$job_info_ref{$jobid}{STATE} .= ":$$job_info_ref{$jobid}{QueueStatus}";
          }
          $val = $$job_info_ref{$jobid}{BlockReason};
          # state can be eligible and blocked...just changed to blocked
          if( defined($val) ){
              $$job_info_ref{$jobid}{NOTE_LONG} = $val;
              $$job_info_ref{$jobid}{STATE}     = "$S_BLOCKED";
              if( $val =~ /violates/i ){
                  $$job_info_ref{$jobid}{NOTE} = "$S_VIOLATES";
              }
              else{
                  $$job_info_ref{$jobid}{NOTE} = "$val";
              }
          }
      }
      $$job_info_ref{$jobid}{STATE} =~ s/^:+//;
      # consistent names
      if( $$job_info_ref{$jobid}{STATE} =~ /^eligible$/i ){
          $$job_info_ref{$jobid}{STATE} = $S_ELIGIBLE;
      }
      elsif( $$job_info_ref{$jobid}{STATE} =~ /^blocked$/i ){
          $$job_info_ref{$jobid}{STATE} = $S_BLOCKED;
      }
      # time remaining (to start or left in queue)
      if( defined($$job_info_ref{$jobid}{"start_s"}) &&
          defined($$job_info_ref{$jobid}{"dur_s"}) ){
          my $s = $$job_info_ref{$jobid}{"dur_s"} - $$job_info_ref{$jobid}{"start_s"};
          $$job_info_ref{$jobid}{"remain_s"} = $s;
          %time_info = &conv_time(SECS=>$s);
          my $hms = $time_info{hms};
          if( $time_info{SECS_TOTAL} < 0 ){
              $hms = "-$hms";
          }
          $$job_info_ref{$jobid}{"remain"} = $hms;
      }
  }
}

#............................................................................
#...Name
#...====
#... get_job_info_showstart
#...
#...Purpose
#...=======
#... Fills out the job_info hash from a command.
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_showstart{
  my %args = (
              JOB_INFO  => undef,
              @_,
             );
  my(
     $command,
     $i,
     $ierr,
     $job_info_ref,
     $jobid,
     @jobids,
     $key,
     $line,
     @lines,
     $loc,
     $num,
     $output,
     $showstart,
     $start,
     $start_col,
     $type,
    );
  #...invalid args
  foreach $key ( keys %args ) {
      if( $key !~ /^(JOB_INFO)$/) {
          $ierr = 1;
          &print_error( "Invalid argument to get_job_info_ [$key]",
                        $ierr );
          exit( $ierr );
      }
  }
  $showstart = &which_exec( "showstart", QUIET=>"" );
  if( $showstart !~ /\S/ ){
      return;
  }
  # get type of showstart (llnl or moab)
  $output = `$showstart -h 2>&1`;
  if( $output =~ /LC\s+Hotline/i ){
      $loc = "llnl";
  }
  else{
      $loc = "moab";
  }
  $job_info_ref = $args{JOB_INFO};
  undef( @jobids );
  foreach $jobid ( sort keys %$job_info_ref ){
      # expensive since need to be done for each job...so only do
      # if start not defined yet and not for dependent jobs
      # and not in a userhold state
      # and not blocked
      if( ! defined( $$job_info_ref{$jobid}{start} ) &&
          ! defined( $$job_info_ref{$jobid}{parent} ) &&
          $$job_info_ref{$jobid}{STATE} !~ /$S_BLOCKED/
          ){
          push( @jobids, $jobid );
      }
  }
  $num = $#jobids + 1;
  $i = 0;
  if( $loc eq "moab" ){
      $command = "$showstart -eall";
  }
  else{
      $command = "$showstart";
  }
  if( $#jobids >= 0 ){
      print "Running '$command <jobid>' for $num jobs:\n";
  }
  foreach $jobid ( @jobids ){
      $i++;
      &status_bar($i, $num);
      $output = &run_command( COMMAND=>"$command $jobid", TIMEOUT=>"5s" );
      @lines = split( /\n/, $output );
      if( $loc eq "moab" ){
          foreach $line ( @lines ){
              if ( $line =~ /Estimated (\S+) based start in\s+(\S+)/ ) {
                  $type = $1;
                  $start = $2;
                  $$job_info_ref{$jobid}{start} = $start;
                  # Rsv seems most accurate
                  if( $type eq "Rsv" ){
                      last;
                  }
                  # probably down if one of them is infinity
                  if( $start =~ /INFINITY/i ){
                      last;
                  }
              }
          }
      }
      else{
          foreach $line ( @lines ){
              if( $line =~ /^(.*START_TIME)/ ){
                  $start_col = length($1);
              }
              elsif( $line =~ /^\s*${jobid}\s+/ ){
                  $line = substr( $line, 0, $start_col);
                  if( $line =~ /(\S+)$/ ){
                      $$job_info_ref{$jobid}{start} = $1;
                  }
              }
          }
      }
  }
}

#............................................................................
#...Name
#...====
#... get_job_info_scontrol
#...
#...Purpose
#...=======
#... Fills out the job_info hash from a command.
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_scontrol{
    my( $cmd_ref ) = shift( @_ );
    my %args = (
        JOB_INFO  => undef,
        @_,
        );
    my(
        $block,
        %blocks,
        $command,
        $date,
        $dates,
        $done,
        $factor,
        $i,
        $ierr,
        $job_info_ref,
        $jobid,
        $key,
        @lines,
        $output,
        $retry,
        $secs,
        %time_info,
        $type,
        $type1,
        %val,
        $val_str,
        $val1,
        );

    #...invalid args
    foreach $key ( keys %args ) {
        if( $key !~ /^(JOB_INFO)$/) {
            $ierr = 1;
            &print_error( "Invalid argument to get_job_info_scontrol [$key]",
                          $ierr );
            exit( $ierr );
        }
    }

    $job_info_ref = $args{JOB_INFO};

    # retry for a bit if no results
    $command = "scontrol show job";
    $done = "false";
    $retry = 0;
    while( $done eq "false" ){
        $output = &run_command( COMMAND=>$command, TIMING=>"", TIMEOUT=>"5s" );
        $done = "true";
        # should start with <data> - otherwise temporary problem
        if( $output !~ /JobId=/ ){
            $done = "false";
            print "Temporary info issue.\n[$output]\n";
            $retry++;
            if( $retry < $RETRY_NUM ){
                print "Retrying...\n";
                sleep( 5 );
                next;
            }
            # give up with error returned (so can be checked for)
            else{
                $ierr = 1;
                &print_error( "Repeated tries at command failed:",
                              "  $command",
                              "Giving up...",
                              $ierr );
                exit( $ierr );
            }
        }
    }

    #$output = "";
    #$output .= "JobId=1\n";
    #$output .= "workdir=/a/b/c\n";
    #$output .= "Priority=1\n";
    #$output .= "JobId=3\n";
    #$output .= "Priority=7\n";
    #$output .= "Dependency=foobar:1\n";
    #$output .= "JobId=2\n";
    #$output .= "Priority=75\n";
    #$output .= "StartTime=2011-10-03T16:00:00\n";
    # remove leading/trailing whitespace with return at end
    @lines = split( /\n/, $output );
    @lines = grep( s/^\s*(.*?)\s*$/$1\n/, @lines );
    $i = 0;
    undef $jobid;
    while( $i <= $#lines ){
        if( $lines[$i] =~ /^JobId=(\d+)/ ){
            $jobid = $1;
            $blocks{$jobid} = "\n";
        }
        if( defined($jobid) ){
            $blocks{$jobid} .= $lines[$i];
        }
        $i++;
    }
    
    # go through each block (job)
    foreach $jobid ( sort keys %blocks ){
        undef( %val );
        
        $val{JOBID} = $jobid;
        
        $block = $blocks{$jobid};
        if( $block =~ /\bName=(\S+)/ ){
            $val{NAME} = $1;
        }
        if( $block =~ /\bWORKDIR=(.*)\n/im ) {
            $val{dir} = $1;
        }
        
        # user
        if( $block =~ /\bUserId=(\S+)/im ) {
            $val{USER} = $1;
            $val{USER} =~ s/\(.*//;
        }

        # dur and dur_s
        if( $block =~ /\bTimeLimit=(\S+)/im ) {
            $val_str = $1;
            %time_info = &conv_time(STRING=>"A$val_str");
            $val{dur} = $time_info{hms};
            $val{dur_s} = $time_info{SECS_TOTAL};
        }
        
        # state
        if( $block =~ /\bJobState=(\w+)/im ) {
            $val{STATE} = $1;
            $val{STATE} =~ tr/A-Z/a-z/;
            if($block =~ /\bReason=(\S+)/im){
                $val{reason} = $1;
            }
            
            # append reason if not running
            if( $val{STATE} eq "running" ){
                $val{STATE} = $S_RUNNING;
            }
            else{
                $val{STATE} .= ":$val{reason}";
            }
            delete($val{reason});

            # get consistent fields
            # state = S_BLOCKED
            if( $val{STATE} =~ /^JobHeld/ ||
                $val{STATE} =~ /PartitionDown/ ||
                $val{STATE} =~ /PartitionNodeLimit/ ||
                $val{STATE} =~ /PartitionTimeLimit/ ||
                $val{STATE} =~ /ReqNodeNotAvail/ ||
                $val{STATE} =~ /ReqNodeNotAvail/ ){
                $val{NOTE}      = $val{STATE};
                $val{STATE}     = $S_BLOCKED;
            }

            # state = S_ELIGIBLE
            elsif( $val{STATE} =~ /pending:((Resources|Priority|Idle).*)/i ){
                $val{NOTE}  = $1;
                $val{STATE} = $S_ELIGIBLE;
            }

            # state = S_BLOCKED
            if( $block =~ /\bJobState=PENDING/im &&
                $block =~ /\bReason=JobHeldAdmin/im &&
                $block =~ /\bDependency=\(null\)/im &&
                $block =~ /\bEligibleTime=Unknown/im
                ){
                $val{NOTE_LONG} = $block;
                $val{STATE}     = $S_BLOCKED;
                $val{NOTE}      = $S_VIOLATES;
            }

            # state = S_COMPLETED
            # Set to S_COMPLETED regardless of other state/note info
            if( $val{STATE} =~ /(cancelled|completed|failed|timeout):(.*)/i ){
                $val{NOTE}  = $2;
                $val{STATE} = $S_COMPLETED;
                if( $val{NOTE} eq "None" ){
                    delete( $val{NOTE} );
                }
            }
        }

        # remaining time
        if( $block =~ /\bRunTime=(\S+)/im ) {
            $val_str = $1;
            %time_info = &conv_time(STRING=>"A$val_str");
            $val1 = $time_info{SECS_TOTAL};
            if( $val1 > 0 && $block =~ /\bTimeLimit=(\S+)/im ) {
                $val_str = $1;
                %time_info = &conv_time(STRING=>"A$val_str");
                $val1 = $time_info{SECS_TOTAL} - $val1;
                $val{remain_s} = $val1;
                %time_info = &conv_time(SECS=>$val1);
                $val{remain} = $time_info{hms};
            }
        }
        
        # starttime
        if( $block =~ /\bStartTime=(\S+)/im ){
            $date = $1;
            if( $date =~ /(\d{4})\-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})/ ){
                @dates = reverse split( /[^\d]+/, $date );
                $dates[4]--;
                $secs = timelocal(@dates) - time();
                if( $secs > 0 ){
                    %time_info = &conv_time(SECS=>$secs);
                    $$job_info_ref{$jobid}{start_s} = $time_info{SECS_TOTAL};
                    $$job_info_ref{$jobid}{start}   = $time_info{hms};
                }
            }
        }

        # NODELIST
        if( $block =~ /\s+NodeList=(\S+)/ ){
            $val{NODELIST} = $1;
        }
        
        # ACCT
        if( $block =~ /\bAccount=(\S+)/im ) {
            $val{ACCT} = $1;
        }

        # PARTITION
        if( $block =~ /\bPartition=(\S+)/im ) {
            $val{PARTITION} = $1;
        }

        # parent
        if( $block =~ /\bDependency=\S+:(\d+)/im ) {
            $val{parent} = $1;
        }

        # PRIO
        if( $block =~ /\bPriority=(\d+)/im ) {
            $val{PRIO} = $1;
        }

        # NODES/PROCS
        foreach $type1 ( "NumNodes", "NumCPUs" ){
            if( $block =~ /\n${type1}=(\S+)/ ){
                $val1 = $1;
                if( $type1 eq "NumNodes" ){
                    $type = "NODES";
                }
                else{
                    $type = "PROCS";
                }
                $val{$type} = $val1;
                $val{$type} =~ s/\-.*//;
                if( $val{$type} =~ /(\S+)([a-zA-Z])/i ){
                    $val{$type} = $1;
                    $factor = $2;
                    if( $factor =~ /k/i ){
                        $factor = 1024;
                    }
                    elsif( $factor =~ /m/i ){
                        $factor = 1024*1024;
                    }
                    elsif( $factor =~ /g/i ){
                        $factor = 1024*1024*1024;
                    }
                    else{
                        $factor = 1; # not right probably
                    }
                    $val{$type} *= $factor;
                }
            }
        }

        # prune
        if( !defined($val{JOBID}) || $val{JOBID} !~ /^\d+$/ ){
            next;
        }
        if( defined($$cmd_ref{u}) && $val{USER} ne $$cmd_ref{u} ){
            next;
        }
        if( defined($$cmd_ref{j_num} ) && &get_id_num($val{JOBID}) ne $$cmd_ref{j_num} ){
            next;
        }
        
        #...push values onto job list
        foreach $key ( keys %val ) {
            $$job_info_ref{$jobid}{$key} = $val{$key};
        }
    }
}

#............................................................................
#...Name
#...====
#... get_job_info_qstat
#...
#...Purpose
#...=======
#... Fills out the job_info hash from a command.
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_qstat{
  my %args = (
              JOB_INFO  => undef,
              @_,
             );
  my(
     $block,
     %blocks,
     $doit,
     $i,
     $ierr,
     $job_info_ref,
     $jobid,
     $key,
     $line,
     @lines,
     $output,
     %val,
    );
  #...invalid args
  foreach $key ( keys %args )
    {
      if( $key !~ /^(JOB_INFO)$/)
        {
          $ierr = 1;
          &print_error( "Invalid argument to get_job_info_ [$key]",
                        $ierr );
          exit( $ierr );
        }
    }
  # qstat is on redstorm - used for find dir there
  $job_info_ref = $args{JOB_INFO};
  $doit = "false";
  foreach $jobid ( keys %$job_info_ref ){
      if( ! defined $$job_info_ref{$jobid}{dir} ){
          $doit = "true";
          last;
      }
  }
  if( $doit ne "false" ){
      $output = &run_command( COMMAND=>"qstat -f", TIMING=>"", TIMEOUT=>"5s" );
      # they decided to do a weird line/wrap thing - get rid of that
      @lines = split( /\n/, $output );
      $output = "";
      foreach $line ( @lines ){
          # if line is a continuation
          if( $line =~ /^ {8}(.+)$/ ){
              $output .= "$1";
          }
          elsif( $line =~ /^\t(.+)$/ ){
              $output .= "$1";
          }
          else{
              $output .= "\n$line";
          }
          $i++;
      }
      @lines = split( /\n/, $output );
      $i = 0;
      undef $jobid;
      while( $i <= $#lines ){
          if( $lines[$i] =~ /^\s*Job ID:\s+(\d+)/i ){
              $jobid = $1;
              $blocks{$jobid} = "\n";
          }
          if( defined($jobid) ){
              $blocks{$jobid} .= "$lines[$i]\n";
          }
          $i++;
      }
      foreach $jobid ( sort keys %blocks ){
          # skip if not a jobid in hash
          if( ! defined($$job_info_ref{$jobid}) ){
              next;
          }
          undef( %val );
          $block = $blocks{$jobid};
          if ( $block =~ /PBS_O_WORKDIR=([^,\n]*)/ ) {
              ($val{dir} = $1) =~ s/\s+//;
          }
          #...push values onto job list
          foreach $key ( keys %val ) {
              $$job_info_ref{$jobid}{$key} = $val{$key};
          }
      }
  }
}
#............................................................................
#...Name
#...====
#... get_job_info_checkjob
#...
#...Purpose
#...=======
#... Fills out the job_info hash from a checkjob command.
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_checkjob{
  my %args = (
              JOB_INFO  => undef,
              @_,
             );
  my(
     $block,
     %blocks,
     $i,
     $ierr,
     $job_info_ref,
     $jobid,
     $key,
     @lines,
     $output,
     %val,
    );
  #...invalid args
  foreach $key ( keys %args )
    {
      if( $key !~ /^(JOB_INFO)$/)
        {
          $ierr = 1;
          &print_error( "Invalid argument to get_job_info_ [$key]",
                        $ierr );
          exit( $ierr );
        }
    }
  $job_info_ref = $args{JOB_INFO};
  $output = &run_command( COMMAND=>"checkjob -v all", TIMING=>"", TIMEOUT=>"5s" );
  @lines = split( /\n/, $output );
  # remove leading/trailing whitespace with return at end
  @lines = grep( s/^\s*(.*?)\s*$/$1\n/, @lines );
  $i = 0;
  undef $jobid;
  while( $i <= $#lines ){
    if( $lines[$i] =~ /^job\s+(\d+)/ ){
      $jobid = $1;
      $blocks{$jobid} = "\n";
    }
    if( defined($jobid) ){
      $blocks{$jobid} .= $lines[$i];
    }
    $i++;
  }
  foreach $jobid ( sort keys %blocks ){
    undef( %val );
    $block = $blocks{$jobid};
    if( $block =~ /\s+job cannot run\s*\(dependency\s*(\d+)\s*jobcomplete not met/ ){
      $val{"parent"} = $1;
    }
    if ( $block =~ /^AName:\s+(.*?)\n/m ) {
      $val{NAME} = $1;
    }
    if ( $block =~ /^\s*IWD:\s*(.*?)\n/m ) {
      $val{dir} = $1;
    }
    if ( $block =~ /^\s*Reservation.*\((\S+)/m ) {
      $val{start} = $1;
    }
    #...push values onto job list
    foreach $key ( keys %val ) {
      # skip if not a jobid in hash
      if( ! defined($$job_info_ref{$jobid}) ){
        next;
      }
      $$job_info_ref{$jobid}{$key} = $val{$key};
    }
  }
}
#............................................................................
#...Name
#...====
#... get_job_info_mshow
#...
#...Purpose
#...=======
#... Fills out the job_info hash from a mshow command.
#...
#...Arguments
#...=========
#... JOB_INFO     Intent: in
#...              Perl type: pointer to hash
#...              Default: none
#...              The job_info hash.
#...
#... USER         Intent: in
#...              Perl type: scalar
#...              only get for this user
#...
#...Program Flow
#...============
#... 1) go through command line and assign to cmd hash
#............................................................................
sub get_job_info_mshow{
  my %args = (
              JOB_INFO  => undef,
              USER      => undef,
              @_,
             );
  my(
     $field,
     @fields,
     $i,
     $ierr,
     $j,
     $job_info_ref,
     $jobid,
     $key,
     $line,
     @lines,
     $MSHOW,
     $output,
     $type,
     $user,
     %val,
    );
  #...invalid args
  foreach $key ( keys %args )
    {
      if( $key !~ /^(JOB_INFO|USER)$/)
        {
          $ierr = 1;
          &print_error( "Invalid argument to get_job_info_ [$key]",
                        $ierr );
          exit( $ierr );
        }
    }
  $job_info_ref = $args{JOB_INFO};
  $user = $args{USER};
  $MSHOW = &which_exec( "mshow" );
  $output = &run_command( COMMAND=>$MSHOW, TIMING=>"", TIMEOUT=>"5s" );
  # some systems do not have working mshow - try another command
  if( $output =~ /is not authorized to perform task/ ){
    $output = &run_command( COMMAND=>"showq -w", TIMING=>"", TIMEOUT=>"5s" );
  }
  @lines = split( /\n/, $output );
  # remove leading/trailing whitespace
  @lines = grep( s/^\s*(.*?)\s*$/$1/, @lines );
  $type = "undef";
  $i = 0;
  #...each line from mshow
  while ( $i <= $#lines ) {
    #...get to jobs line
    while ( $i <= $#lines && $lines[$i] !~ /^\s*(\S+)\s+jobs\-+/ ) {
      $i++;
    }
    if( $i > $#lines ) {
      last;
    }
    # <type> jobs-----
    if( $lines[$i] =~ /^\s*(\S+)\s+jobs\-+/ ) {
      $type = $1;
      $i++;
    }
    else {
      $ierr = 1;
      &print_error( "Could not find jobs line from $MSHOW.",
                    $ierr );
      exit( $ierr );
    }
    # JOBID <other fields>
    if ( $lines[$i] =~ /^\s*JOBID/ ) {
      $line = $lines[$i];
      $line =~ s/^\s*//;
      $line =~ s/\s*$//;
      @fields = split( /\s+/, $line );
      $i++;
    } else {
      $ierr = 1;
      &print_error( "Unexpected output from $MSHOW.",
                    "Expected to see JOBID header line.",
                    "[$lines[$i]]",
                    $ierr );
      exit( $ierr );
    }
    # skip blank lines
    while ( $lines[$i] !~ /\S/ ) {
      $i++;
    }
    # process this job block
    while ( $i <= $#lines && $lines[$i] =~ /\S/ &&
            $lines[$i] !~ /^(\d+) $type jobs/) {
      undef( %val );
      $val{type} = $type;
      ($line = $lines[$i]) =~ s/^\s*//;
      $line =~ s/\s*$//;
      # process each field in a line
      for ( $j = 0; $j <= $#fields; $j++ ) {
        $field = $fields[$j];
        if ( $field eq "STARTTIME" ||
             $field eq "QUEUETIME" ) {
          $field = "qstime";
          ($val{$field} = $line) =~ s/\s+/_/g;
          $line = "";
        } elsif ( $line =~ /\s*(\S+)\s*(.*)/ ) {
          $val{$field} = $1;
          $line = $2;
        }
        if ( $line =~ /\S/ && $j == $#fields ) {
          $ierr = 1;
          &print_error( "Data left on job line but no more fields left",
                        "Remaining line: $line",
                        "Original  line: $lines[$i]",
                        "Fields:         ".join(", ", @fields),
                        $ierr
                      );
          exit( $ierr )
        }
      } # DONE: process each field in a line
      # error if still data on line
      if ( $line =~ /\S/ ) {
        $ierr = 1;
        &print_error( "No more data on job line left but still fields to process",
                      "Field:          $field",
                      "Original  line: $lines[$i]",
                      "Fields:         ".join(", ", @fields),
                      $ierr
                    );
        exit( $ierr )
      }
      # next line
      $i++;
      #...push values onto job list
      if ( defined($user) && $val{USERNAME} ne $user ) {
        next;
      } # if user
      $jobid = $val{JOBID};
      foreach $key ( keys %val ) {
        $$job_info_ref{$jobid}{$key} = $val{$key};
      }
    } # DONE: this job block
  } # DONE: each line from mshow
}

sub rj_numerically { $a <=> $b; }

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
sub parse_args {
  my(
     $argv_ref,
     $cmd_ref
    ) = @_;
  my(
     @args,                     # arguments
     $done,
     $ierr,                     # error ret val
     $num_args,                 # number of arguments
     $opt,                      # current option
     $val,                      # value for current option
     @vals,
    );
  $ierr = 0;
  @args = @{$argv_ref};
  #....................
  #...parse the args...
  #....................
  $num_args = $#args;
  while ( @args ) {
    $opt = shift( @args );
    #............
    #...-<opt>...
    #............
    if ( $opt =~ /^-(h|l|s)$/ ) {
      $opt = $1;
      $$cmd_ref{$opt} = "true";
    }
    #............
    #...--help...
    #............
    elsif ( $opt =~ /^--(help)$/ ) {
      $opt = "h";
      $$cmd_ref{$opt} = "true";
    }

    #..............
    #...--(word)...
    #..............
    elsif ( $opt =~ /^--(debug|cancel)$/ ) {
      $opt = $1;
      $$cmd_ref{$opt} = "true";
    }

    # --opt <val>
    elsif ( $opt =~ /^--(batch)$/ ) {
      $opt = "$1";
      if ( ! @args ) {
          $ierr = 1;
          &print_error( "Value needed for option [$opt].",
                        $ierr );
          exit( $ierr );
      }
      $val = shift( @args );
      if( $val !~ /^(no|only|yes)$/ ){
          $ierr = 1;
          &print_error("Invalid option for [$opt] <$val> != <no|only|yes>",
                       $ierr);
          exit( $ierr );
      }
      $$cmd_ref{$opt} = $val;
    }

    # --opt <comma sep val put into array>
    elsif ( $opt =~ /^--(state)$/ ) {
      $opt = "$1";
      if ( ! @args ) {
          $ierr = 1;
          &print_error( "Value needed for option [$opt].",
                        $ierr );
          exit( $ierr );
      }
      $val = shift( @args );
      @{$$cmd_ref{$opt}} = split(/,/, $val);
    }

    # --<opt> <glob>
    elsif( $opt =~ /--(dirsid)$/ ){
        $opt = $1;
        if ( ! @args ) {
            $ierr = 1;
            &print_error( "Value needed for option [$opt].",
                          $ierr );
            exit( $ierr );
        }
        $val = shift( @args );
        push( @{$$cmd_ref{$opt}}, glob( $val ) );
    }

    #................
    #...-<opt> val...
    #................
    elsif ( $opt =~ /^-(d|n|u|f|j)$/ ) {
      $opt = $1;
      if ( ! @args ) {
        $ierr = 1;
        &print_error( "Value needed for option [$opt].",
                      $ierr );
        exit( $ierr );
      }
      if ( $opt eq "n" ) {
        $done = "false";
        while ( $done eq "false" ) {
          if ( ! @args || $args[0] =~ /^\-/ ) {
            $done = "true";
            next;
          }
          $val = shift( @args );
          @vals = split(/,/, $val);
          foreach $val ( @vals ) {
              unshift( @{$$cmd_ref{$opt}}, $val );
          }
        }
      }
      elsif ( $opt eq "j" ){
          $val = shift( @args );
          $$cmd_ref{$opt} = $val;
          $$cmd_ref{"${opt}_num"} = &get_id_num($val);
      }
      elsif ( $opt eq "d" ){
          $val = shift( @args );
          @vals = glob( $val );
          push( @{$$cmd_ref{$opt}}, @vals );
      }
      else {
          $val = shift(@args);
          $$cmd_ref{$opt} = $val;
      }
    }
    #..............
    #...filename...
    #..............
    else {
      $ierr = 1;
      &print_error( "Invalid option [$opt]",
                    $ierr );
      exit( $ierr );
    }
  }
  #...defaults
  #...long format as well if doing -l
  if( defined($$cmd_ref{l}) && ! defined($$cmd_ref{f}) ){
    $$cmd_ref{f} = 2;
  }
  # user
  if( ! defined($$cmd_ref{u}) && ! defined($$cmd_ref{j}) ) {
      $$cmd_ref{u} = $ENV{LOGNAME};
  }
  if( defined($$cmd_ref{u}) && $$cmd_ref{u} =~ /^all$/i ){
      delete( $$cmd_ref{u} );
      if( ! defined($$cmd_ref{f}) ){
          $$cmd_ref{f} = 4;
      }
  }
  if( defined($$cmd_ref{j}) && $$cmd_ref{j} =~ /^all$/i ){
      delete( $$cmd_ref{j} );
      delete( $$cmd_ref{j_num} );
      if( ! defined($$cmd_ref{f}) ){
          $$cmd_ref{f} = 4;
      }
  }

  if( ! defined( $$cmd_ref{batch}) ){
      $$cmd_ref{batch} = "yes";
  }

  # if not specified, self has more info (mtime, time, dt), others have less
  if( ! defined($$cmd_ref{f}) || $$cmd_ref{f} !~ /^\d+$/ ){
      if( defined( $$cmd_ref{u} ) && $$cmd_ref{u} ne $ENV{LOGNAME} ){
          $$cmd_ref{f} = "4";
      }
      else{
          $$cmd_ref{f} = "1";
      }
  }
  return( $ierr );
}
