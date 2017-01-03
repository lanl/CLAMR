package System::Workstation;

use strict;
use warnings;
use Carp;
use Data::Dumper;

use System::Platform;
use vars qw(@ISA);
@ISA = qw(System::Platform);



{
  my %attrs = (
               batch          =>   '',
               mpirun         =>   'mpirun',
               serialrun      =>   '',
	       default_queue  =>   '',
               lsfdir        =>   '',
	       min_batch_pe   =>   1,
	      );

  sub new {
    my ($invocant, $arg_ref) = @_;
    my $class = ref($invocant) || $invocant;
    my $href = ();
    %$href = (%attrs, %$arg_ref);
    my $self = $class->SUPER::new($href);
    return $self;
  }

}


sub mpirun{
  my $self = shift;

  my $pe = $self->SUPER::pe() || 0;
  $self->SUPER::mpirun("mpirun -np $pe");
}



# We need to make sure that mpirun changes if pe does so modify the default
# behavior to update mpirun after updating pe.
sub pe{
  my $self = shift;

  my $pe = $self->SUPER::pe(@_) || 0;
  $self->mpirun();
  return $pe;
}



sub memory_usage{
  my ($self) = @_;

  my $physical = 0;
  my $free = 0;
  my $ubc = 0;
  my $inactive = 0;

  my $vmstat_cmd = "/usr/bin/free";

  if ( -e $vmstat_cmd )
  {
    my @tmp = qx($vmstat_cmd);
    my ($total, $used, $avail) = $tmp[1] =~ /^Mem:\s+ (\d+) \s+ (\d+) \s+ (\d+) \s+/x;

    if ( (defined $total) && ($total > 0) )
    {
      $physical = $total*8*1024;
    }

    if ( (defined $avail) && ($avail > 0) )
    {
      $free = $avail * 8 * 1024;
    }
  }
  else
  {
    my $host = qx(hostname);
    chomp($host);
    print STDERR "($host) \"$vmstat_cmd\" does not exist\n";
  }

  return ($physical, $free, $ubc, $inactive);

}

sub swap_usage{
  my $self = shift;

  my $avail = 0;
  my $alloc = 0;

  my $swapon_cmd = "/sbin/swapon";

  if (-e $swapon_cmd ) {

    my @tmp = qx($swapon_cmd -s);

    my ($total, $used) = $tmp[1] =~ /partition \s+ (\d+) \s+ (\d+) \s+/x;

    $total = ($total)? $total: 0;
    $used  = ($used)? $used: 0;

    my $free = $total - $used;

    $avail = $free  * 1024 * 8;
    $alloc = $total * 1024 * 8;
  }

  else {
    my $host = qx(hostname);
    chomp($host);
    print STDERR "($host) \"$swapon_cmd\" does not exist\n";
  }

  return($alloc, $avail);
}


sub get_num_cpus_avail
{
  my $self = shift;

  my $ncpus = 1;
  return $ncpus;
}

sub batch_preamble{
  my $preamble = "";

  return $preamble;
}

sub batch_post{
  return '';
}

sub lsbatch_dir {
  my ($self) = shift;

  $self->{lsbatch_dir} = shift if @_;
  my $dir = $self->{lsbatch_dir};

  return $dir;
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

