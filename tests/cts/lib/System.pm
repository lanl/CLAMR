package System;

use strict;
use warnings;
use Carp;
use Tools::Load;
use DumpStack;
use Data::Dumper;
use Shell::Source;
use subs qw(nobatch system_name serialrun);
use vars qw(@ISA $VERSION $AUTOLOAD);
$VERSION = '0.01';

# System is a singleton
my $singleton_system;

# Inherit methods from local_system.pm if it exists.
  BEGIN {
    my $sysName = `hostname`;
    my $system = "";
  }


 {
  # System hash with default values.
  my %system = (
		NOBATCH              => undef,
		SYSTEM_NAME          => undef,
		SERIALRUN            => "",
	      );

  sub _init{
    my ($self, $arg_ref) = @_;

    # Augment input hash with default values. 
    $arg_ref ||= {};
    %$arg_ref = (%system, %$arg_ref);


    # Inherit batch methods
    $self->load_batch($arg_ref);

    # Inherit mpi methods.
    $self->load_mpi();

    #print "DEBUG: \@ISA: ", @ISA, "\n";


    while (my ($func, $value) = each %$arg_ref) {
      $func = lc $func;
      $value = $self->$func($value);
    }

    # experimental: look for system specific information in .cts file
    #parse_system_info();

  }

  sub AUTOLOAD {
    my $self = shift;
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);

    return if ($method =~ /destroy/i);
    $method = uc $method;
    unless (exists $system{$method}) {
      dump_stack;
      local $\ = "\n";
      print "System::ERROR:: No such attribute: $method is defined for this system.\n";
      print "The available attributes are";
      foreach (keys %system) {print}
      croak "Aborting";
    }

    $self->{$method} = shift if @_;
    my $val = $self->{$method};
    return $val;
  }
}

sub new {
  #print caller, "\n";
  unless ($singleton_system) {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    $singleton_system = {};
    bless $singleton_system => $class;
    $singleton_system->_init(@_);
    #$singleton_system->dump;
  }
  return $singleton_system;
}

sub load_batch {
  #   print "\n\nDEBUG_LOAD_BATCH:\n";
  #   my $i = 0;
  #   while (my ($package, $filename, $line) = caller($i++)) {
  #     print "$package, $filename, $line\n";
  #   }
  my $self = shift;
  my $aref = shift;
  my $batch_system;
  my $exec_check;
  my $output;

  if (defined $aref->{NOBATCH}){
    delete $aref->{NOBATCH};
    $batch_system = "Batch";
  }
  elsif (defined $aref->{BATCH}) {
    die "Unknown batch system requested \'$aref->{BATCH}\'. Choose moab, pbs, lsf or lcrm \n"
      unless $aref->{BATCH} =~ /\b(moab|pbs|lsf|lcrm)\b/i;
    $batch_system = "Batch::" . uc $aref->{BATCH};
    delete $aref->{BATCH};
  }
  else {
    $batch_system = "";
    my $msub = `which msub 2>/dev/null`;
    chomp $msub;
    my $bsub = `which bsub 2>/dev/null`;
    chomp $bsub;
    my $psub = `which psub 2>/dev/null`;
    chomp $psub;
    my $qsub = `which qsub 2>/dev/null`;
    chomp $qsub;
    #...pick the batch system in order of perference...
    if( $batch_system eq "" )
      {
        if( -e $bsub )
          {
            $exec_check = `which bqueues 2>/dev/null`;
            chomp $exec_check;
            if( -e $exec_check && $exec_check !~ /MOAB/ )
              {
                #...check to see if this is real LSF or a moab clone
                #...need to have full LSF for using LSF
                $output = `$exec_check -h 2>&1 `;
                if( $output !~ /brief information about the Moab queue/ )
                  {
                    $batch_system = "Batch::LSF";
                  }
              }
          }
      }
    if( $batch_system eq "" )
      {
        if( -e $msub )
          {
            $batch_system = "Batch::MOAB";
          }
      }
    #...for now, use pbs on redstorm (even though moab - does not work)
    if( $batch_system eq "" || $qsub =~ /torque/ )
      {
        if( -e $qsub )
          {
            $batch_system = "Batch::PBS";
          }
      }
    if( $batch_system eq "" )
      {
        if( -e $psub )
          {
            $batch_system = "Batch::LCRM";
          }
      }
    #...default
    if( $batch_system eq "" )
      {
        $batch_system = "Batch";
      }
  }

  load $batch_system;
  my $batch = $batch_system->new();
  push @ISA, $batch_system;
  #print "DEBUG: \@ISA: ", @ISA, "\n";
  %$self = (%$self, %$batch);  # This step is not necessary but is
                               # useful for debugging purposes.
}


sub load_mpi {
  my $self = shift;

  my $mpirun_rsh = `which mpirun_rsh 2>/dev/null`;
  chomp $mpirun_rsh;
  my $mpirun = `which mpirun 2>/dev/null`;
  chomp $mpirun;
  #print "DEBUG:System: mpirun = $mpirun\n";
  my $prun   = `which prun 2>/dev/null`;
  chomp $prun;
  my $yod    = `which yod 2>/dev/null`;
  chomp $yod;
  my $psub   = `which psub 2>/dev/null`;
  chomp $psub;
  # if some sort of xt-mpt module loaded, use that flavor
  my $aprun;
  if( defined($ENV{MPICHBASEDIR})){
    $aprun  = grep( /xt/, $ENV{MPICHBASEDIR} );
  }
  elsif( defined($ENV{CRAY_MPICH2_VERSION}) ){
    $aprun  = $ENV{CRAY_MPICH2_VERSION};
  }
  else{
    $aprun = "";
    
  }

  my $machine = `uname -n`;
  # unset psub so mpirun is used
  if( $machine =~ /^dawn\d/ ){
    $psub = "";
  }

  # we let run_job.pl set mpi
  #my $mpi_system = (-e $psub)       ? "MPI::LCRM"            :
  #                 (-e $mpirun_rsh) ? "MPI::Mvapich"         :
  #                 (-e $prun)       ? "MPI::Hp"              :
  #                 (-e $yod)        ? "MPI::Yod"             :
  #                 ($aprun =~ /\S/) ? "MPI::CieloMPI"        :
  #                 (-e $mpirun)     ? _get_mpi_type($mpirun) :
  #                                    "MPI"                  ;
  my $mpi_system = "MPI";

  load "MPI";     # We only have one file at the moment, so we load it
                  # no matter which mpi we are actually using.
  my $mpi = $mpi_system->new();
  push @ISA, $mpi_system;
  %$self = (%$self, %$mpi);  # This step is not necessary but is
                               # useful for debugging purposes.
}

sub _get_mpi_type {
  my ($file) = @_;

  # Default to unknown
  my $type = "Unknown";

  # check if file is binary or a shell script
  if (-B $file) {
    my $out = `mpirun -V 2>&1`;
    if ($out =~ /la-mpi/ims) {$type = "MPI::Lampi"}
    if ($out =~ /open\s*mpi/i) {$type = "MPI::Openmpi"}
    if ($out =~ /BG\/P mpirun/i) {$type = "MPI::BG"}
    if ($out =~ /mpich/i) {$type = "MPI::Mpich"}
  }

  # file is a script, grep for mpich
  else {
    my $grep_count = `grep -ci mpich $file`;
    if ($grep_count >= 1) {$type = "MPI::Mpich"}
  }

  warn "System::_get_mpi_type:warning Unknown mpirun\n" unless ($type ne "Unknown");
  return $type;
}


sub parse_system_info{
  # @_ may contain a hash of system configuration information
  # System information garnered from the .cts file is passed as a
  # string in the "CTSFILE_SYSTEM_INFO" element.
  my $system_hash = $_[0] || {};


  my $system_string = $system_hash->{CTSFILE_SYSTEM_INFO} || "";

  delete $system_hash->{CTSFILE_SYSTEM_INFO} if defined $system_hash->{CTSFILE_SYSTEM_INFO};

  # $system_string my contain lines with initialization information
  # for specific systems. We need to separate this information into a
  # hash whose keys are the system names.
  my %hash = ();
  while ($system_string =~ /system \s+ (\w+) \s+          # get the system name
	                   (?! \s* end \s* system) (.*?)   # get
                                                          # everything
                                                          # up to the
                                                          # end system line
	                   end \s* system/ixmsg ) {
    my $system_name = $1;
    my $system_info = $2;
    while ($system_info =~ /^\s* (\w+) \s* = (.*)$/xmg) {
      $hash{$system_name}{$1} = $2;
    }
  }
  #print Dumper(%hash);
  foreach my $sys ( keys %hash) {
    my $id_str = $hash{$sys}{id} || "";
    my ($quote, $command, $pat) = $id_str =~ /(["''"]) ([^\1]*) \1 \s+ (.*) \s* $/x;
    #print "\nDEBUG_UNAME: $command: $pat\n";
    my $output = `$command`;
    chomp $output;
    my $match = '$output =~ ' . $pat;
    #print "DEBUG_OUTPUT: $match\n";
   if (eval $match) {
     delete $hash{$sys}{id};
     foreach my $param(keys %{$hash{$sys}}) {
       $system_hash->{$param} = $hash{$sys}{$param};
     }
     last;
    }
  }
  #print "DEBUG_PARSE ", Dumper($system_hash);
}

sub source{
  my $self = shift;
  if (my $source_file = shift){
    my $tcsh = Shell::Source->new(shell => "tcsh", file => $source_file);
    $tcsh->inherit;
  }
}


sub date{
  my $self = shift;

  my $date = `date +%m%d%y`;
  chomp $date;
  return $date;
}


sub tar{
  my ( $calling_package, $filename, $line) = caller;
  print "debug_tar:  $calling_package, $filename, $line\n";
  my ($self, $tar_file, $dir) = @_;

  `rm -rf $tar_file`;
  `tar -cf $tar_file $dir`;
}


sub psi{
  my ($self, $file, $hpss_dir) = @_;

  `psi save $hpss_dir/$file`;
}


sub system_name{
  my $self = shift;
  my $sysname = shift if @_;

  $self->{SYSTEM_NAME} = $sysname if $sysname;
  # system_name should be a scalar. If it is an array then
  # replace it with the first element. 
  $self->{SYSTEM_NAME} = $self->{SYSTEM_NAME}[0] if (ref $self->{SYSTEM_NAME});
  #dump_stack;
  unless (defined $self->{SYSTEM_NAME}) {
    $self->{SYSTEM_NAME} = `hostname`;
    chomp $self->{SYSTEM_NAME};
  }
  #dump_stack "DEBUG: SYSTEM_NAME is  $self->{SYSTEM_NAME}\n";
  return $self->{SYSTEM_NAME};
}


sub serialrun{
  my $self=shift;
  my $serialrun = "";
  my $bpsh = `which bpsh 2> /dev/null`;


  chomp $bpsh;

  if (-x $bpsh) {
    $serialrun = 'setenv NODE `echo $NODES | sed "s/[,-].*//"`' . "\n" . 'bpsh $NODE '; 
  }
  # needed for new moab on LCRM
  my $mpirun = $self->mpirun();
  if( $mpirun =~ /\bpoe\b/ )
    {
      $serialrun = "setenv SLURM_LL_API_NUM_NODES 1\nsetenv SLURM_LL_API_NUM_TASKS 1\n";
    }

  return $serialrun;
}


sub dump {
  my $self = shift;
  my $type = ref $self;
  print "$type: ", Dumper($self);
}

1;

__END__

=head1 NAME

System - Perl extension for platform specific setup and batch execution.

=head1 SYNOPSIS

  use System;
my $platform = System::new(
			   {

			   }
			  )

=head1 ABSTRACT

  This should be the abstract for System.
  The abstract is used when making PPD (Perl Package Description) files.
  If you don't want an ABSTRACT you should also edit Makefile.PL to
  remove the ABSTRACT_FROM option.

=head1 DESCRIPTION

Stub documentation for System, created by h2xs.

=head2 EXPORT

None




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







