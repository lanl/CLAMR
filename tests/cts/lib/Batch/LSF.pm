package Batch::LSF;

use strict;
use warnings;
use DumpStack;
use Carp;
use Env;
use Shell::Source;
use Data::Dumper;
use Batch;
use subs qw(batch queue queues lsfdir default_time_limit pe time_limit blocked_queues);
use vars qw(@ISA $queues $AUTOLOAD);
@ISA = qw(Batch);

{
  # Batch hash with default values.
  my %batch = (
		BATCH          	    =>   (-x "/lsf/bin/bsub")? "/lsf/bin/bsub <" :'',
	        BLOCKED_QUEUES      =>   [],
		QUEUE         	    =>   undef,
		QUEUES         	    =>   undef,
		LSFDIR         	    =>   '/lsf/bin/',
	        DEFAULT_PE          =>   1,
	       );

  sub _init{
    my ($self, $arg_ref) = @_;

    # Augment input hash with default values. 
    $arg_ref ||= {};
    %$arg_ref = (%batch, %$arg_ref);

    foreach my $func(keys %$arg_ref) {
      my $value = $arg_ref->{$func};

      $func = lc $func;
      $value = $self->$func($value);
    }
  }

  sub AUTOLOAD {
    my $self = shift;
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);

    $method = uc $method;
    return if ($method =~ /DESTROY/)
;
    unless (exists $batch{$method}) {
      local $\ = "\n";
      print "ERROR:: No such attribute: $method is defined for this batch.\n";
      print caller, "\n";
      print "The available attributes are";
      foreach (keys %batch) {print}
      return undef;
    }

    $self->{$method} = shift if @_;
    my $val = $self->{$method};
    return $val;
  }

}

sub new {
  my $invocant = shift;
  my $class = ref($invocant) || $invocant;
  my $self = {};
  bless $self => $class;
  $self->_init(@_);
  return $self;
}


sub batch{
  my $self = shift;
  my $batch = shift;

  if (defined $batch) {
    $batch = ($batch =~ /</) ? $batch : "$batch <";
    $self->{BATCH} = $batch;
  }
  return $self->{BATCH};
}


sub _get_available_lsf_queues{
  my $self = shift;
  my ($proc_limit_string, @queue_names, @queue_data);

  my $output = `bqueues -u $ENV{LOGNAME}`;
  while ($output =~ m/^(.*)$/mg ){
    my $line = $1;
    my ($queue_name) = split /\s+/, $line;
    next if grep {/$queue_name/i} $self->blocked_queues;
    if(
       $line =~ /Open/     # We are only interested in open queues.
       )
      {
        push @queue_names, $queue_name;
      }
  }
  
  foreach my $queue_name (@queue_names) {
    push @queue_data, $self->_get_lsf_queue_data($queue_name);
  }

  # Sort the queues by priority. This has already been done by the
  # bqueues command for the systems we have looked at so we won't do
  # it unless we find a system that requires it.

  #@queue_data = sort {$b->{PRIO} <=> $a->{PRIO}} @queue_data;
  return wantarray? @queue_data: \@queue_data;
}


sub _get_lsf_queue_data{
  my $self = shift;
  my $name = shift or die "System::Batch::LSF::_get_lsf_queue_data: Error: A queue name must be provided.";

    my $queue_data = {};
    my @number = ();
    $queue_data->{NAME}=$name;
    my $output = `bqueues -l $name`;
    my @fields = split /^\s*([_\w ]+):\s/mi, $output;
    shift @fields;
    my %fields = @fields;

    my $key;
    $key = 'SCHEDULING POLICIES';
    if( defined( $fields{$key} ) )
      {
        $queue_data->{$key}=$fields{$key};
      }
    else
      {
        $queue_data->{$key}="";
      }
    my $queue_parameters = $fields{'QUEUE'} || "";
    my $prio = 0;
    if ($queue_parameters =~ /\n(.*)PRIO/m) {
      my $offset = length $1;
      if ($queue_parameters =~ /\n(?:.*)PRIO .* $
                                      \n (?:.{$offset}) \s* (\d+)/mx) {
	$prio = $1;
      }
      else {
	warn "BATCH::LSF: Could not determine the priority of QUEUE $name\n";
      }

    }
    $queue_data->{PRIO} = $prio;


    my $max_limits = $fields{'MAXIMUM LIMITS'} || "";
    $max_limits = "\n" . $max_limits;  # parsing crutch

    my $proc_limit_string = "1";
    if ($max_limits =~ /\n(.*)PROCLIMIT/m) {
      my $offset = length $1;
      ($proc_limit_string) = $max_limits =~ /\n(?:.*)PROCLIMIT\s*$
                                      \n(?:.{$offset}) (.*)$/mx;
    }

    $queue_data->{MIN_PE} = 1;
    @number = split /\s+/, $proc_limit_string;
    my $size = @number;
    if ($size == 3 ){
      $queue_data->{MIN_PE}=$number[0];
      $queue_data->{MAX_PE}=$number[2];
    }
    elsif ($size == 2 ){
      $queue_data->{MIN_PE}=$number[0];
      $queue_data->{MAX_PE}=$number[1];
    }
    elsif ($size == 1 ){
      $queue_data->{MAX_PE}=$number[0];
    }
    else {
      print "Error in determining batch processor limits\n";
      exit 1;
    }

    my $time_limit = undef;
    if ($max_limits =~ /\n(.*)RUNLIMIT/m) {
      my $offset = length $1;
      ($time_limit) = $max_limits =~ /\n(?:.*)RUNLIMIT [^\n]*
                                      \n(?:.{$offset}) \s* ([\d.]+)/mx;
    }
    $queue_data->{TIME_LIMIT}=$time_limit;

  return $queue_data;
}

sub blocked_queues{
  my $self = shift;
#   print "\n\nDEBUG_BLOCKED_QUEUES\n";
#   my $i = 0;
#   while (my ($package, $filename, $line) = caller($i++)) {
#     print "$package, $filename, $line\n";
#   }

  my $blocked_queues = shift || $self->{BLOCKED_QUEUES} || [];
  croak "System::Batch::LSF::blocked_queues: blocked_queues must contain an array\n" 
    unless (ref $blocked_queues eq "ARRAY");

  $self->{BLOCKED_QUEUES} = $blocked_queues;
  return wantarray? @$blocked_queues: $blocked_queues;  
}


sub queues{
  my $self = shift;
#   print "\n\nDEBUG_QUEUES\n";
#   my $i = 0;
#   while (my ($package, $filename, $line) = caller($i++)) {
#     print "$package, $filename, $line\n";
#   }

  my $queues = shift || $self->{QUEUES} || $self->_get_available_lsf_queues;
  croak "System::Batch::LSF::queues: queues must contain an array\n" 
    unless (ref $queues eq "ARRAY");

  # Turn queue names into ref to hash of queue information.
  foreach my $queue(@$queues) {
    #print Dumper($queue);
    next if (ref $queue eq "HASH");
    $queue = $self->_get_lsf_queue_data($queue);
  }

  $self->{QUEUES} = $queues;
  #dump_stack;
  return wantarray? @$queues: $queues;  
}


sub _minutes{
  my $time = shift;

  my $hour = 0;
  my $min  = 0;

  if (($hour, $min) = $time =~ /^(\d*):(\d+)/x) {
    $hour ||= 0;
    $min += 60*$hour;
  }
  else {
    $min = $time; # assume that the time is given in minutes.
  }

  return $min;
}

sub set_queue {
  my $self = shift;
  my $qname = shift || die "System::Batch::LSF::set_queue: A queue name is required.\n";

  $self->{QUEUE} = $qname;
  return $self->{QUEUE};
}


sub queue {
  my $self = shift;
  my $qname = $self->{QUEUE} || "";
  return $qname if $qname;

  my $pe = shift;
  my $time_limit = shift;

  if ((defined $pe) && (defined $time_limit)) {

    $time_limit = _minutes($time_limit);

    my @queues = $self->queues;

    # Choose the best queue. The queues are already ordered by priority,
    # so the first queue to match the criteria should be the best queue.
    foreach my $queue (@queues) {
      my $queue_name = $queue->{NAME};

      # next if queue is in the blocked_queue list
      next if grep /$queue_name/, $self->blocked_queues;
      
      # next if ONLY_INTERACTIVE
      next if $queue->{'SCHEDULING POLICIES'} =~ /ONLY_INTERACTIVE/;

      # The following takes care of queues which don't have a time limit.
      my $queue_time_limit = $queue->{TIME_LIMIT} || $time_limit;

      # only select a queue which matches the criteria
      next unless (($queue->{MAX_PE} >= $pe) &&
		   ($queue_time_limit >= $time_limit)
		  );
      $qname = $queue_name;
      last;
    }
    warn "Warning: No suitable queue was found." unless $qname
  }

  # If there is only one argument then assume it is the name of the
  # queue the user wants to use until another queue is specified.
  elsif (defined $pe) {
    my $qname = $pe;
    $self->{QUEUE} = $qname;
  }

  # If no arguments are given then return $self->{QUEUE} or the empty string.
  else {
    $qname = $self->{QUEUE} || "";
  }

  return $qname;
}

sub min_batch_pe{
  my $self = shift;
  my $queue = shift or die "ERROR: BATCH::LSF::MIN_BATCH_PE  A queue name is required\n";
  my @queues = $self->queues;
  my $min_pe = 1;

  foreach my $que (@queues) {
    next unless ($que->{NAME} eq $queue);
    $min_pe = $que->{MIN_PE};
    last;
  }

  return $min_pe;
}

sub batch_pe{
  my $self = shift;

  my $queue = shift or die "ERROR: BATCH::LSF::BATCH_PE  A queue is required\n";

  my $pe = shift  or die "ERROR: BATCH::LSF::BATCH_PE  PE is required\n";

  my $min_batch_pe = $self->min_batch_pe($queue);

  return ($pe >= $min_batch_pe)? $pe: $min_batch_pe;
}

sub batch_preamble{
  my ($self, $arg) = @_;

  my $pe             = $arg->{PE} || $self->default_pe;
  my $time_limit     = $arg->{TIME_LIMIT}   || croak("A TIME_LIMIT is required");
  my $queue          = $arg->{QUEUE}        || $self->queue($pe, $time_limit);
  my $batch_pe       = $self->batch_pe($queue, $pe);
  my $jobname        = $arg->{JOBNAME}      ||  croak("A JOBNAME is required");
  my $logfile        = $arg->{LOGFILE}      || "bsub_" . $jobname;

  my $preamble = <<EOF;
#
#BSUB -K
#BSUB -n $batch_pe
#BSUB -q $queue
#BSUB -J $jobname
#BSUB -W $time_limit
#BSUB -o $logfile
EOF

  return $preamble;
}


sub batch_post{
  return '';
}

sub bjobs{
  my ($self, $user) = @_;

  $user = ($user)? "-u $user": '';
  my $lsf_cmd = "bjobs $user -w";

  my @lsf_output = `$lsf_cmd`;

  printf("\nThe command: \'%s\' returned:\n%s",$lsf_cmd, join("",@lsf_output));
  return @lsf_output;
}


sub lsbatch_dir {
  my ($self) = shift;

  $self->{lsbatch_dir} = shift if @_;
  my $dir = $self->{lsbatch_dir};

  unless (defined $dir) {
    my $home = $self->home();
    $self->{lsbatch_dir} = "$home/.lsbatch";
  }

  $dir = $self->{lsbatch_dir};
  unless ((-d $dir) && (-r $dir)) {
    croak ("Either $dir does not exist or you do not have read permission.")
  }
  return $dir;
}


sub parse_bjobs{
  my $self = shift;
  my $line = shift || carp("No line to parse.");

  # This regex pattern grabs the jobid, state and job name
  my $lsf_job_pat = qr{(\d+)\s+     # jobid
		       (\S+)\s+     # user
		       (\S+)\s+     # stat
		       \S+\s+       # queue
		       \S+\s+       # from_host
		       \S+\s+       # exec_host
		       (\S+)\s+     # job_name - the first non-blank part
		      }x;

  if ( my ($jobid, $user, $stat, $name) = $line =~ /^$lsf_job_pat/) {
    return ($jobid, $user, $stat, $name);
  }

  return ();
}


sub batch_jobs{
  my $self = shift;
  my $user = shift || $self->user();
  my %data = ();

  my @lsf_output =  $self->bjobs($user);

  my $lsf_job_num = 0;
  foreach my $line ( @lsf_output ) {
    if ( my ($jobid, $user, $stat, $name) = $self->parse_bjobs($line)) {
      $data{$jobid} = {
		       'name' => $name,
                       'stat' => $stat,
                       'user' => $user,
		      };
      $lsf_job_num++;
    }
  }

  print "Found $lsf_job_num jobs in LSF\n";
  print "\nThe LSFDATA hash contains:\n";
  foreach my $jobid ( sort keys %data ){
      printf("JOBID=<%s> NAME=<%s> LSF STATE=<%s> USER=<%s>\n",
	     $jobid, $data{$jobid}->{name}, $data{$jobid}->{stat},
	     $data{$jobid}->{user});
    }

  return %data;
}

sub kill_batch_jobs{
  my $self = shift;

  my $arg = shift || carp("No jobs to kill.");
  return unless $arg;

  my $num_tries = 3;      # Number of bkill tries
  my $sleep_time = 60*4;  # Number of seconds between bkill commands

  my $jobs = (ref $arg)? $arg: ($arg);

  my %data = $self->batch_jobs();
  foreach my $jobid(@$jobs) {
    next unless (defined $data{$jobid});

    print STDOUT "Bkilling job $jobid with the bkill -s 9 command\n";
    my @kill_args = ('bkill', '-s 9', $jobid);
    my $attempts = 0;
    while (($attempts < $num_tries) && batch(@kill_args)) {
      sleep($sleep_time);
      $attempts++;
    }
  }

  # Make sure we killed the jobs.
  %data = ();
  %data = $self->batch_jobs();
  my $home = $self->home();
  foreach my $jobid(@$jobs) {
    next unless (defined $data{$jobid});

    carp("Can not kill job $jobid .... problem with the job or LSF");
    carp("Creating a DO_NOT_RUN file in $home");

    my $pname = $data{$jobid}->{name};
    my $do_not_run = "$home/$pname-DO_NOT_RUN";
    system "touch $do_not_run";
  }
  return;
}


sub submit_batch_job{
  my ($self, $arg) = @_;

  my $batch      = $self->batch;
  my $batch_pe   = $self->batch_pe($arg->{pe});
  my $queue      = $arg->{queue}        || $self->default_queue;
  my $jobname    = $arg->{jobname}      || $self->jobname || croak("A jobname is required");
  my $time_limit = $arg->{time_limit}   || $self->time_limit;
  my $logfile    = $arg->{logfile}      || $self->{jobname};
  my $exe        = $arg->{executable}   || croak("An executable is required.");

  croak("$exe is not executable.") unless (-e $exe);

  my @sys_args = ($batch, " -q $queue -o $logfile -e $logfile -n $batch_pe -W $time_limit -J $jobname $exe");
  return batch(@sys_args);
}



sub get_num_cpus_avail
{
  my $self = shift;

  carp("The subroutine get_num_cpus_avail() has not yet been implemented on this system. Defaulting to 1 cpu.");
  my $ncpus = 1;

  return $ncpus;
}



1;


=head1 AUTHOR

David L. Aubrey, E<lt>dla@lanl.govE<gt>

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



