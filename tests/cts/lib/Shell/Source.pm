# Copyright 1997-2001, Paul Johnson (pjcj@cpan.org)

# This software is free.  It is licensed under the same terms as Perl itself.

# The latest version of this software should be available from my homepage:
# http://www.pjcj.net

use strict;

require 5.004;

package Shell::Source;

use vars qw($VERSION);

$VERSION = "0.01";

use Carp;
use FileHandle;
use Data::Dumper;

my $shells =
{
    csh       => "csh -f -c 'source [[file]]; env' |",
    tcsh      => "tcsh -f -c 'source [[file]]; env; set' |",
    tcsh_init => "tcsh -f -c set |",
    sh        => "sh -c '. [[file]]; env' |",
    ksh       => "ksh -c '. [[file]]; env' |",
    zsh       => "zsh -c '. [[file]]; env' |",
    bash      => "bash -norc -noprofile -c '. [[file]]; env' |",
};

sub new
{
    my $class = shift;
    my $self = { @_ };
    croak "Must specify type of shell" unless $self->{shell};
    my $init_run = $self->{shell}."_init";
    $self->{init_run} = $shells->{$init_run};
    $self->{run} ||= $shells->{$self->{shell}};
    croak "Must specify how to run unknown shell $self->{shell}"
        unless $self->{run};
    push @{$self->{ignore}}, qw( TIMEFMT PWD _ );
    bless $self, $class;
    $self->init;
    $self->run if length $self->{file};

    delete $self->{init_shell_variables};
    #print Dumper($self);
    $self
}

sub init
{
    my $self = shift;
    my $run = $self->{init_run};
    my $fh = $self->{fh}
           = FileHandle->new($run) or croak "Can't run $self->{init_run}";
    $self->_parse_init;
    $fh->close or croak "Can't close $self->{shell}";

    $self
}

sub run
{
    my $self = shift;
    my $file = shift || $self->{file};
    croak "Must specify file to source" unless length $self->{file};
    (my $run = $self->{run}) =~ s/\[\[file\]\]/$self->{file}/g;
    my $fh = $self->{fh}
           = FileHandle->new($run) or croak "Can't run $self->{shell}";
    $self->_parse;
    $fh->close or croak "Can't close $self->{shell}";
    $self
}

sub _parse
{
    my $self = shift;
    my $fh = $self->{fh};                         # FileHandle ready for reading
    my $env = 0;                           # for control of multi-line variables
    while (defined(my $line = <$fh>)) {
      if ($line =~ /^(\w+)=(.*)$/) {
	$env = 1;
	if ((!defined $ENV{$1} || $ENV{$1} ne $2) &&
	    !grep {$1 eq $_} @{$self->{ignore}}) {
	  $self->{env}{$1} = $2;
	}
      }
      elsif ($line =~ /^(\w+) \s+ (.*)$/x) {
	if (!defined $self->{init_shell_variables}{$1} || $self->{init_shell_variables}{$1} ne $2) {
	  $self->{shell_variables}{$1} = $2;
	}
      }
      else {
	push (@{$self->{output}}, $line) unless $env;
      }
    }
    $self
}

sub _parse_init
{
    my $self = shift;
    my $fh = $self->{fh};                         # FileHandle ready for reading
    while (defined(my $line = <$fh>))
      {
	if ($line =~ /^(\w+) \s+ (.*)$/x)
	  {
	    $self->{init_shell_variables}{$1} = $2;
	  }
      }
    $self
  }

sub inherit
{
    my $self = shift;
    while (my ($key, $val) = each (%{$self->{env}}))
    {
        $ENV{$key} = $val;
    }
}

sub shell
{
    my $self = shift;
    my $shell = "";
    while (my ($key, $val) = each (%{$self->{env}}))
    {
        $shell .= qq($key="$val"; export $key\n);
    }
    $shell
}

sub output
{
    my $self = shift;
    join("\n", @{$self->{output}}) if defined $self->{output}
}

sub env
{
    my $self = shift;
    $self->{env}
}

1;

__END__

=head1 NAME

Shell::Source - run programs and inherit environment changes

=head1 SYNOPSIS

 use Shell::Source;
 my $csh = Shell::Source->new(shell => "csh", file => "stuff.csh");
 $csh->inherit;
 print STDERR $csh->output;
 print $csh->shell;

=head1 DESCRIPTION

The Shell::Source allows arbitrary shell scripts, or other programs for
that matter, to be run and their environment to be inherited into a Perl
program.

Begin by creating a Shell::Source object, and specifying the shell it
will use.

If the shell is unknown to the module, you will also need to specify how
to run the shell in such a way that the output is a series of lines of
the form NAME=value.  For example, to run a csh script:

 my $csh = Shell::Source->new(shell => "csh",
                              file  => "stuff.csh",
                              run   => "csh -f -c 'source [[file]]; env' |");

However, for known shells this is not required.  Note that [[file]] will
be replaced with the filename of the program you want to run.

Output from running the program is returned from $csh->output.

Changes made to the environment by running the program may be inherited
by calling $csh->inherit.

The environment changes are available as a hash from $csh->env, or in
Bourne shell syntax from $csh->shell.

=head1 BUGS

Huh?

=head1 VERSION

Version 0.01 - 2nd August 2001

=head1 HISTORY

Created - Wednesday 26th November 1997 09:29:31 pm

=head1 LICENCE

Copyright 1997-2001, Paul Johnson (pjcj@cpan.org)

This software is free.  It is licensed under the same terms as Perl itself.

The latest version of this software should be available from my homepage:
http://www.pjcj.net

=cut
