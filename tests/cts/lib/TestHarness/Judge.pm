package TestHarness::Judge;

# A Judge object is basically an array. Each element of the array contains two
# items. The first is a link to a TestProblem object(the runner). The second is
# a condition to be compared to the status of the runner.
#
# The concept goes as follows. We have a group of runners approaching the start
# gate. The starter(test harness) asks the next runner in line if he is ready to
# start. The runner has a list of judges whom he queries. Each judge will
# provide one of the following answers, "GO", "WAIT", "DO_NOT_RUN".
#
# If all the judges in the runner's list say "GO" or if the list is empty, the
# the runner starts. If any of the judges say "DO_NOT_RUN" then the runner
# removes himself from the run queue. If any judge says "WAIT" then the runner
# goes to the back of the run queue.
#
# A Judge may also have lists of warnings or errors to be printed on either
# success or failure.

use strict;
use warnings;
use Carp;
use Data::Dumper;
use vars qw($AUTOLOAD);


sub new {
  my $invocant = shift;
  my $class = ref($invocant) || $invocant;
  my $self = {};
  bless $self => $class;
#   my $string = shift or croak "TestHarness::Judge::new:  No string was supplied";
#   $self->parse_list($string);
  return $self;
}


sub parse_list{
  my ($self, $string) = @_;
  # $string should look lik this
  # (!noh, regression.suite, simple)
  # Remove surrounding parens.
  $string =~ s/\s*[()]\s*//g;

  # Split comma delimited list.
  my @tests = split /\s*,\s*/;
}

sub dump {
  my $self = shift;
  print Dumper($self);
}


sub add_criterion{
  my ($self, $test, $not) = @_;
  
}

sub add_test{
  my $self = shift;
}

sub remove_test{
  my $self = shift;
}

sub poll_tests{
  my $self = shift;
}

sub query_test{
  my $self = shift;
}

sub judge{
  my $self = shift;
  return "GO";
}


1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Judge - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Judge;
  blah blah blah

=head1 ABSTRACT

  This should be the abstract for Judge.
  The abstract is used when making PPD (Perl Package Description) files.
  If you don't want an ABSTRACT you should also edit Makefile.PL to
  remove the ABSTRACT_FROM option.

=head1 DESCRIPTION

Stub documentation for Judge, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

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
