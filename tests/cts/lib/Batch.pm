package Batch;

use strict;
use warnings;
use Carp;
use Env;
use Data::Dumper;
use subs qw(batch_post batch min_batch_pe pe time_limit batch_preamble queues blocked_queues account);
use vars qw($AUTOLOAD);

{
  # Batch hash with default values.
  my %batch = (
	       ACCOUNT            => "",
	       BATCH_POST         => "",
	       BATCH_PREAMBLE     => "",
	       BATCH              => "",
	       DEFAULT_PE         => 1,
	       DEFAULT_TIME_LIMIT => 60.0,
	       QUEUES             => undef,
               QOS                => undef,
	       BLOCKED_QUEUES     => undef,
	      );


  sub new {
    my ($invocant, $arg_ref) = @_;
    my $class = ref($invocant) || $invocant;
    my $self;
    $arg_ref ||= {};
    $self = {%batch, %$arg_ref};
    bless $self, $class;
  }

  sub AUTOLOAD {
    my $self = shift;
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);

    return if ($method =~ /destroy/i);
    $method = uc $method;
    unless (exists $batch{$method}) {
      local $\ = "\n";
      print "caller: ", caller;
      print Dumper($self);
      print "BATCH:ERROR:: No such attribute: $method is not defined.\n";
      print "The available attributes are";
      foreach (keys %batch) {print}
      croak "Aborting";
    }

    $self->{$method} = shift if @_;
    my $val = $self->{$method};
    return $val;
  }

}


sub batch_preamble{
  return "";
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














