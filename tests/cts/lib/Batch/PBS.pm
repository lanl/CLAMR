package Batch::PBS;

use strict;
use warnings;
use DumpStack;
use Carp;
use Env;
use Shell::Source;
use Data::Dumper;
use Batch;
use FindBin;
use subs qw(account batch queue queues default_time_limit pe time_limit blocked_queues);
use vars qw(@ISA $queues $AUTOLOAD);
@ISA = qw(Batch);

{
  # Batch hash with default values.
  my %batch = (
	        ACCOUNT             =>   undef,
		BATCH          	    =>   "$FindBin::Bin/qsub.pl",
	        BLOCKED_QUEUES      =>   [],
		QUEUE         	    =>   undef,
		QUEUES         	    =>   undef,
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
    return if ($method =~ /DESTROY/);

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
    $self->{BATCH} = $batch;
  }
  return $self->{BATCH};
}


sub _get_available_pbs_queues{
  my $self = shift;
  my @queue_data = ("standard");

  return wantarray? @queue_data: \@queue_data;
}


sub _get_pbs_queue_data{
    my $self = shift;
    my $name = shift or die "System::Batch::PBS::_get_pbs_queue_data: Error: A queue name must be provided.";

    my $queue_data = {};

    $queue_data->{NAME}       = $name;
    $queue_data->{PRIO}       = 1;
    $queue_data->{MAX_PE}     = 100000;
    $queue_data->{MIN_PE}     = 1;
    $queue_data->{TIME_LIMIT} = 100000; # minutes
    return $queue_data;

    my @number = ();
    $queue_data->{NAME}=$name;
    my $output = `bqueues -l $name`;
    my @fields = split /^\s*([_\w ]+):\s/mi, $output;
    shift @fields;
    my %fields = @fields;

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
  croak "System::Batch::PBS::blocked_queues: blocked_queues must contain an array\n" 
    unless (ref $blocked_queues eq "ARRAY");

  $self->{BLOCKED_QUEUES} = $blocked_queues;
  return wantarray? @$blocked_queues: $blocked_queues;  
}


sub queues{
  my $self = shift;

  my $queues = shift || $self->{QUEUES} || $self->_get_available_pbs_queues;
  croak "System::Batch::PBS::queues: queues must contain an array\n" 
    unless (ref $queues eq "ARRAY");

  # Turn queue names into ref to hash of queue information.
  foreach my $queue(@$queues) {
    #print Dumper($queue);
    next if (ref $queue eq "HASH");
    $queue = $self->_get_pbs_queue_data($queue);
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
  my $qname = shift || die "System::Batch::PBS::set_queue: A queue name is required.\n";

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
  my $queue = shift or die "ERROR: BATCH::PBS::MIN_BATCH_PE  A queue name is required\n";
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

  my $queue = shift or die "ERROR: BATCH::PBS::BATCH_PE  A queue is required\n";

  my $pe = shift  or die "ERROR: BATCH::PBS::BATCH_PE  PE is required\n";

  my $min_batch_pe = $self->min_batch_pe($queue);

  return ($pe >= $min_batch_pe)? $pe: $min_batch_pe;
}

sub batch_preamble{
  my ($self, $arg) = @_;

  my $account        = $arg->{ACCOUNT}      || $self->account;
  my $pe             = $arg->{PE} || $self->default_pe;
  my $time_limit     = $arg->{TIME_LIMIT}   || croak("A TIME_LIMIT is required");
  my $queue          = $arg->{QUEUE}        || $self->queue($pe, $time_limit);
  my $batch_pe       = $self->batch_pe($queue, $pe);
  my $jobname        = $arg->{JOBNAME}      ||  croak("A JOBNAME is required");
  my $logfile        = $arg->{LOGFILE}      || "qsub_" . $jobname;

  # account for last field in pbs being seconds not minutes
  $time_limit = "$time_limit:00";

# need to use internal qsub.pl that waits (simulating -K behavior
  my $preamble = <<EOF;
#
# #PBS -K (qsub.pl simulates "wait" behavior)
#PBS -l select=$batch_pe
#PBS -A FY115533
#PBS -q $queue
#PBS -j oe
#PBS -N $jobname
#PBS -l walltime=$time_limit
#PBS -o $logfile
#PBS -V
EOF

  if (defined $account){
    $preamble .= "#PBS -A $account\n";
  }

  return $preamble;
}


sub batch_post{
  return '';
}


sub get_num_cpus_avail
{
  my $self = shift;

  carp("The subroutine get_num_cpus_avail() has not yet been implemented on this system. Defaulting to 1 cpu.");
  my $ncpus = 1;

  return $ncpus;
}

sub account{
  my $self = shift;

  my $account = shift || $self->{ACCOUNT};
  $account = (ref $account)? @$account[0]: $account;
  $self->{ACCOUNT} = $account;
  #...set
  if( ! defined($account) )
    {
      #...get from env var if it is set
      if( defined($ENV{ACCOUNT_BATCH}) ) {
        $self->{ACCOUNT} = $ENV{ACCOUNT_BATCH};
      }
      elsif (defined($ENV{ACCOUNT})) {
        $self->{ACCOUNT} = $ENV{ACCOUNT};
      }
      #...if not found yet
      if( ! defined $self->{ACCOUNT} )
        {
          my $debug = "";
          my $account_def;
          my $account_choose;
          my $command = "mdiag -u $ENV{LOGNAME}";
          my $output = `$command 2>&1`;
          $debug .= "Command: [$command]\nOutput:  [$output]\n";
          if( $output =~ /\bADEF\s*=\s*(\S+)/ )
            {
              $account_def = $1;
            }
          #...create list of valid accounts
          my %accounts_valid = ();
          if( $output =~ /ALIST\s*=\s*(\S+)/ )
            {
              my @accounts_available = split(/,/,$1);
              foreach $account ( @accounts_available )
                {
                  $command = "mdiag -a $account";
                  $output = `$command 2>&1`;
                  $debug .= "Command: [$command]\nOutput:  [$output]\n";
                  if( $output =~ /
                                  Name\s+          # 1
                                  Priority\s+      # 2
                                  Flags\s+         # 3
                                  QDef\s+          # 4
                                  QOSList\*\s+     # 5
                                  PartitionList\s+ # 6
                                  Target\s+        # 7
                                  Limits\s*        # 8
                                  /x )
                    {
                      if( $output =~ /
                                      ($account)\s+
                                      (\S+)\s+
                                      (\S+)\s+
                                      (\S+)\s+
                                      (\S+)\s+
                                      (\S+)\s+
                                      (\S+)\s+
                                      (\S+)
                                      /x )
                        {
                          # pick first one with valid QOSList
                          my $QOSList = $5;
                          if( $QOSList ne "-" )
                            {
                              $accounts_valid{$account} = $QOSList;
                            }
                        }
                    }
                }
            }
          # choose access account (developer account)
          if( ! defined($account_choose) &&
              defined( $accounts_valid{access} ) )
            {
              $account_choose = "access";
            }
          # if valid, choose default account
          if( ! defined($account_choose) &&
              defined($account_def) &&
              defined( $accounts_valid{$account_def} ) )
            {
              $account_choose = $account_def;
            }
          #...go through list of valid accounts
          if( ! defined($account_choose) )
            {
              #...just pick first one
              foreach $account ( keys %accounts_valid )
                {
                  $account_choose = $account;
                  last;
                }
            }
          if( defined($account_choose) )
            {
              $self->{ACCOUNT} = $account_choose;
              # if a standby qos exists and no other qos,
              # use it to get job scheduled faster.
              # So far, only on some machines (eg purple)
              if( $accounts_valid{$account_choose} =~ /\bstandby\b/ )
                {
                  my $qos = $self->qos;
                  if( ! defined($qos) )
                    {
                      $qos = "standby";
                      $self->qos($qos);
                    }
                }
            }
          if( ! defined($self->{ACCOUNT}) )
            {
              print "ERROR:: Unable to find any valid accounts.\n";
              print "List of accounts: mdiag -u, look at ALIST\n";
              print "Valid Account:    mdiag -a <account> with QOSList* defined\n";
              print "$debug\n";
              print caller, "\n";
              return undef;
            }
        }
    }
  return $self->{ACCOUNT};
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



