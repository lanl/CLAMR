package Tie::Epsilon;

use strict;
use warnings;
use Tie::Hash;

use vars qw(@ISA);
@ISA = qw(Tie::StdHash);

my $default = .02;   # I just pulled .02 out of my hat. It is only used if the
                     # user doesn't set a default epsilon and one is needed.

sub FETCH {
  my ($self, $key) = @_;

  if (exists $self->{$key}) {return $self->{$key}}
  else {
    warn "Tie::Epsilon - Warning: No epsilon was found for $key. Default epsilon $self->{default} will be used.\n";
    $self->{$key} = $self->{default} || $default;
    return $self->{$key};
  }
}

1;
