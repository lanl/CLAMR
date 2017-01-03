package Batch::LCRM;

use strict;
use warnings;
use DumpStack;
use Carp;
use Env;
use Shell::Source;
use Data::Dumper;
use Batch;
use FindBin;
use POSIX;
use subs qw(account batch queue queues default_time_limit pe time_limit blocked_queues);
use vars qw(@ISA $queues $AUTOLOAD);
@ISA = qw(Batch);

{
  # Batch hash with default values.
  my %batch = (
	        ACCOUNT             =>   undef,
		BATCH          	    =>   "$FindBin::Bin/psub.pl",
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

sub batch_pe{
  my $self = shift;

  my $pe = shift  or die "ERROR: BATCH::LCRM::BATCH_PE  PE is required\n";

  return ($pe );
}

sub batch_preamble{
  my ($self, $arg) = @_;

  my $pe             = $arg->{PE} || $self->default_pe;
  my $batch_pe       = $self->batch_pe($pe);
  my $jobname        = $arg->{JOBNAME}      ||  croak("A JOBNAME is required");
  my $time_limit     = $arg->{TIME_LIMIT}   || croak("A TIME_LIMIT is required");
  my $logfile        = $arg->{LOGFILE}      || "psub_" . $jobname;

  my $cwd = getcwd();
  my ($nodes, $ppn );

  $ppn = 8;
  $nodes = ceil($pe/$ppn);

  my $preamble = <<EOF;
#
#PSUB -v
#PSUB -ln $nodes
#PSUB -eo
#PSUB -g $batch_pe
#PSUB -r $jobname
#PSUB -tW $time_limit
#PSUB -tM $time_limit
#PSUB -o $logfile
#PSUB -standby
EOF

  return $preamble;
}


sub batch_post{
  return '';
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



