package MPI;

use strict;
use warnings;
use Carp;
use Env;
use POSIX;
use Data::Dumper;
use subs qw(mpirun mpirun_args mpirun_pre);
use vars qw($AUTOLOAD);

{
  # Mpi hash with default values.
  my %mpi = (
		MPIRUN         	    => 'mpirun',
		MPIRUN_ARGS         => '',
                MPIRUN_PRE          => '',
	    );


  sub new {
    my ($invocant, $arg_ref) = @_;
    my $class = ref($invocant) || $invocant;
    my $self;
    $arg_ref ||= {};
    $self = {%mpi, %$arg_ref};
    bless $self, $class;
  }

  sub AUTOLOAD {
    my $self = shift;
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);

    return if ($method =~ /destroy/i);
    $method = uc $method;
    unless (exists $mpi{$method}) {
      local $\ = "\n";
      print "MPI:ERROR:: No such attribute: $method is not defined.\n";
      print "The available attributes are";
      foreach (keys %mpi) {print}
      croak "Aborting";
    }

    $self->{$method} = shift if @_;
    my $val = $self->{$method};
    return $val;
  }
}


sub mpirun{
  my $self = shift;
  return 'mpirun -np $NUM_CPUS';
}
sub mpirun_args{
  my $self = shift;
  return '';
}
sub mpirun_pre{
  my $self = shift;
  return '';
}


sub dump {
  my $self = shift;
  my $type = ref $self;
  print "$type: ", Dumper($self);
}

##############################################################################
##############################################################################
package MPI::Lampi;
use strict;
use warnings;
use Carp;

use vars qw(@ISA);
@ISA = qw(MPI);



##############################################################################
##############################################################################
package MPI::Mvapich;
use strict;
use warnings;
use Carp;

use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;
  my $host_file = (defined $ENV{HOSTFILE})? "-hostfile $ENV{HOSTFILE}": "";

  return 'mpirun_rsh -np $NUM_CPUS ' . $host_file;
}
sub mpirun_args{
  my $self = shift;
  return '';
}
sub mpirun_pre{
  my $self = shift;
  return '';
}


##############################################################################
##############################################################################
package MPI::Mpich;

use strict;
use warnings;
use Carp;
use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;
  my $machine_file = (defined $ENV{MACHINEFILE})? "-machinefile $ENV{MACHINEFILE}": "";

  return 'mpirun -np $NUM_CPUS ' . $machine_file;
}
sub mpirun_args{
  my $self = shift;
  return '';
}
sub mpirun_pre{
  my $self = shift;
  return '';
}

##############################################################################
##############################################################################
package MPI::Openmpi;
use strict;
use warnings;
use Carp;
use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;
  my $host_file = (defined $ENV{HOSTFILE})? "-hostfile $ENV{HOSTFILE}": "";
  my $args;
  my $hostname = `uname -a`;
  if( $hostname =~ /acme/ )
    {
      $args = "--bynode";
    }
  else
    {
      $args = "";
    }
  return 'mpirun -np $NUM_CPUS ' . $args . ' ' . $host_file;
}
sub mpirun_args{
  my $self = shift;
  return "";
}
sub mpirun_pre{
  my $self = shift;
  return '';
}

##############################################################################
##############################################################################
package MPI::CieloMPI;
use strict;
use warnings;
use Carp;
use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;
  return 'aprun -n $NUM_CPUS';
}
sub mpirun_args{
  my $self = shift;
  return "";
}
sub mpirun_pre{
  my $self = shift;
  return '';
}

##############################################################################
##############################################################################
package MPI::BG;
use strict;
use warnings;
use Carp;
use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;
  my $host_file = (defined $ENV{HOSTFILE})? "-hostfile $ENV{HOSTFILE}": "";
  my $args;
  # 4GB/node
  # modes:
  # smp:  1 rank/node = 1 ppn
  # dual: 2 rank/node = 2 ppn
  # vn:   4 rank/node = 4 ppn
  # changes must be consistent in MOAB.pm (ppn) and MPI.pm (mode)
  $args = "-mode smp";
  return 'mpirun -np $NUM_CPUS ' . $args . ' ' . $host_file;
}
sub mpirun_args{
  my $self = shift;
  return "";
}
sub mpirun_pre{
  my $self = shift;
  return '';
}

##############################################################################
##############################################################################
package MPI::Hp;
use strict;
use warnings;
use Carp;

use vars qw(@ISA);
@ISA = qw(MPI);


sub mpirun{
  my $self = shift;

  return 'prun -n $NUM_CPUS';
}
sub mpirun_args{
  my $self = shift;
  return '';
}
sub mpirun_pre{
  my $self = shift;
  return '';
}

##############################################################################
##############################################################################
package MPI::Yod;
use strict;
use warnings;
use Carp;

use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;

  return 'yod -sz $NUM_CPUS ';
}
sub mpirun_args{
  my $self = shift;
  return '';
}
sub mpirun_pre{
  my $self = shift;
  return '';
}

##############################################################################
##############################################################################
package MPI::LCRM;
use strict;
use warnings;
use Carp;
use POSIX;

use vars qw(@ISA);
@ISA = qw(MPI);

sub mpirun{
  my $self = shift;

  return 'poe ';
}
sub mpirun_args{
  my $self = shift;
  my $pe = shift or die "ERROR: MPI::LCRM::MPIRUN_ARGS pe is required\n";

  #...number of process per node (hardwired...shound change)
  my $ppn = 8;
  #...pool - again, hardwired and should get somehow
  my $pool = "viz";
  return '-procs $NUM_CPUS -nodes `perl -e "use POSIX; print ceil($NUM_CPUS/'.$ppn.')"` -rmpool '.$pool;
}
sub mpirun_pre{
  my $self = shift;
  my $pe = shift or die "ERROR: MPI::LCRM::MPIRUN_PRE pe is required\n";

  #...number of process per node (hardwired...shound change)
  my $ppn = 8;
  return 'setenv SLURM_LL_API_NUM_NODES `perl -e "use POSIX; print ceil($NUM_CPUS/'.$ppn.')"`'."\n".'setenv SLURM_LL_API_NUM_TASKS $NUM_CPUS'."\n";
}

##############################################################################
##############################################################################
package MPI::Unknown;
use strict;
use warnings;
use Carp;

use vars qw(@ISA);
@ISA = qw(MPI);


sub mpirun{
  my $self = shift;

  return 'mpiexec -n $NUM_CPUS';
}
sub mpirun_args{
  my $self = shift;
  return '';
}
sub mpirun_pre{
  my $self = shift;
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






