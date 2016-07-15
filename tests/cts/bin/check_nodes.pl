eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

# NOTE: You can use FindBin in the script. The modules will automatically have access to $FindBin::Bin.
use FindBin qw($Bin);
use lib ("$Bin", "$Bin/lib", "$Bin/../lib");

use my_utils qw (
                 expand_string
                 get_jobinfo
                 print_error
                 print_perl_obj
                 sort_unique
                 status_bar
                 which_exec
                );
#..........
#...vars...
#..........
my(
   %cmd,
  );
#...parse args
$|=1; # flush stdout
&parse_args( \@ARGV, \%cmd );
if( defined( $cmd{h} ) )
  {
    print <<'EOF';
#............................................................................
#...Name
#...====
#... check_nodes.pl
#...   Sanity check for a set of nodes.
#...   The current machine is also checked (named: ".")
#...
#...Usage
#...=====
#... check_nodes.pl [list of nodes] [-j <job id>]
#...
#...   [list of nodes]
#...      List of nodes.  Simple regexp's are allowed:
#...         Output from ljobs:
#...          a,b[1-3],c[4-5,9]d => a,b1,b2,b3,c4d,c5d,c9d]
#...         Output from checkjob:
#...          [a:4][b:4][c:4]   => a,b,c
#...
#...Example
#...=======
#... o check_nodes.pl -j 123 "tua[017-020,029,031,033,037]"
#...   Get info from the nodes of job 123, the other specific nodes
#...   listed, and for the current node (".").
#...   You can get a list of hosts a job has by doing:
#...      ljobs -l <job id>
#...      checkjob <job id>
#............................................................................
EOF
  exit;
  }
#...get the number of hosts
$num_hosts = $#{$cmd{hosts}} + 1;
print "Hosts: ", join(",",@{$cmd{hosts}}),"\n";
print "Processing $num_hosts hosts.\n";
#.................................
#...commands to do on each host...
#.................................
# uname
$cmd_check .= ";echo \"CN:uname\";uname -n";
# df various dirs to see space and mounts
# do not check scratch systems since, when they go down, they hang the command...ugh
#$cmd_check .= ";echo \"CN:df\";df -k . /netscratch/$ENV{LOGNAME} /tmp /scratch/$ENV{LOGNAME} /scratch1/$ENV{LOGNAME} /scratch2/$ENV{LOGNAME} /scratch3/$ENV{LOGNAME} /var/spool ~";
$cmd_check .= ";echo \"CN:df\";df -k . /netscratch/$ENV{LOGNAME} /tmp /var/spool ~";
$cmd_check .= ";echo \"CN:w\";w";
$cmd_check =~ s/^;//;
#..................
#...end commands...
#..................
$num = 0;
# ssh onto each host and do commands
foreach $host ( @{$cmd{hosts}} ){
    $num++;
    &status_bar( $num, $num_hosts );
    if( $host eq "localhost" ||
        $host eq "." ){
        $cmd = $cmd_check;
    }
    else{
        $cmd = "ssh -Y $host '$cmd_check'";
    }
    $output = `($cmd) 2>&1`;
    $info{$host} = $output;
}
# process each output
print "\n";
foreach $host ( sort keys %info ){
    @lines = split( /\n/, $info{$host} );
    $state = "";
    foreach $line ( @lines ){
        # state
        if( $line =~ /^CN:(\S+)/ ){
            $state = $1;
            next;
        }
        # lines
        if( $state eq "df" ){
            if( $line =~ /(\d+)\%\s+(\S+)/ ){
                $filesystems{$2}{$host} = "$1% full";;
            }
            elsif( $line =~ /\`(\S+)\': No such file or directory/ ){
                $filesystems{$1}{""} = "No";
            }
        }
        if( $state eq "w" ){
            if( $line =~ /(\d+) user(s?).*(load\s+average.*)/ ){
                $w_users{$host} = $1;
                $w_load{$host} = $3;
            }
        }
    }
}
# print sorted on filesystems
print "\ndf info:\n";
print "========\n";
foreach $filesystem ( sort keys %filesystems ) {
    $not_found = "";
    foreach $host ( @{$cmd{hosts}} ) {
        if( ! defined( $filesystems{$filesystem}{$host} ) ){
            $not_found .= ",$host";
        }
    }
    print "  $filesystem\n";
    foreach $host ( sort keys %{$filesystems{$filesystem}} ){
        if( $host eq "" ){
            next;
        }
        printf( "    %10s %10s\n", $host, $filesystems{$filesystem}{$host} );
        # only print first one if panfs mounted
        if( $filesystem =~ m&^/panfs/& ||
            $filesystem =~ m&^/users/& ||
            $filesystem =~ m&^/netscratch/& ){
            last;
        }
    }
    if( $not_found =~ /\S/ ){
        $not_found =~ s/^,//;
        print "    NOT_FOUND: $not_found\n";
    }
}
# print out info from w
print "\nw info:\n";
print "=======\n";
foreach $host ( sort keys %w_load ){
    printf( "  %10s %s %4s %s\n", $host, "users:", $w_users{$host}, $w_load{$host} );
}
print "\n";
print "Other:\n";
print "======\n";
print "Run 'run_status.pl' to get information about your jobs.\n";
# exit
print "\n";
#...case insensitive sort
sub sort_case_insensitive
  {
    my( $a, $b );
    lc($a) cmp lc($b);
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
       %jobinfo, # info about the job
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
    #....................
    #...parse the args...
    #....................
    $num_args = $#args;
    while( @args )
      {
        $opt = shift( @args );
        #..........
        #...help...
        #..........
        if( $opt =~ /^-+(h(elp)?)$/i )
          {
            $opt = "h";
            $$cmd_ref{$opt} = "true";
          }
        #..................
        #...-j <job ids>...
        #..................
        elsif( $opt =~ /^-(j)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 1;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                exit( $ierr );
              }
            @vals = split( /,/, shift( @args ) );
            foreach $val ( @vals ){
                undef( %jobinfo );
                &get_jobinfo( $val, \%jobinfo );
                &expand_string( $jobinfo{hosts_string}, \@{$$cmd_ref{hosts}} );
            }
          }
        #..................
        #...-<opt> <val>...
        #..................
        elsif( $opt =~ /^-(e|d|n)$/ )
          {
            $opt = $1;
            if( ! @args )
              {
                $ierr = 1;
                &print_error( "Value needed for option [-$opt].",
                              $ierr );
                exit( $ierr );
              }
            $val = shift( @args );
            $$cmd_ref{$opt} = "$val";
          }
        #.............
        #...--<opt>...
        #.............
        elsif( $opt =~ /^--(noprune|nostack)$/ )
          {
            $opt = $1;
            $$cmd_ref{$opt} = "true";
          }
        #...........
        #...hosts...
        #...........
        else
          {
              &expand_string( $opt, \@{$$cmd_ref{hosts}} );
          }
      }
    #..........................
    #...DONE: parse the args...
    #..........................
    # always do local machine
    push( @{$$cmd_ref{hosts}}, "." );
    @{$$cmd_ref{hosts}} = sort_unique(\@{$$cmd_ref{hosts}});
    return( $ierr );
  }
