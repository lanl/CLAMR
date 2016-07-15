package Tools::Check;
use strict;
use Carp qw[carp];

BEGIN {
    use Exporter    ();
    use vars        qw[ @ISA $VERSION @EXPORT_OK $VERBOSE $STRICT_TYPE];

    @ISA        =   qw[ Exporter ];
    $VERSION    =   0.01;

    @EXPORT_OK  =   qw[check];
}

sub check {
    my $tmpl    = shift;
    my $href    = shift;
    my $verbose = shift || $VERBOSE || 0;

    ### lowercase all args ###
    my $args = {};
    %$args = map { lc, $href->{$_} } keys %$href;

    ### flag to set if something went wrong ###
    my $flag;

    for my $key ( keys %$tmpl ) {

        ### check if the required keys have been entered ###
        my $rv = _hasreq( $key, $tmpl, $args );

        unless( $rv ) {
            carp qq[Required key '$key' missing for ] . _who_was_it() if $verbose;
            $flag++;
        }
    }
    return undef if $flag;

    ### set defaults for all arguments ###
    my $defs = _hashdefs($tmpl);

    ### check if all keys are valid ###
    for my $key ( keys %$args ) {

        unless( _iskey( $key, $tmpl ) ) {
            carp qq[Key '$key' is not a valid key for ] . _who_was_it() if $verbose;
            next;
        } elsif ( $tmpl->{$key}->{no_override} ) {
            carp qq[You are not allowed to override key '$key' for ] . _who_was_it() if $verbose;
            next;
        } else {

            ### flag to set if the value was of a wrong type ###
            my $wrong;

            if( defined $tmpl->{$key}->{allow} ) {

                my $what = $tmpl->{$key}->{allow};

                ### it's a string it must equal ###
                ### this breaks for digits =/
                unless ( ref $what ) {
                    $wrong++ unless $args->{$key} eq $what;

                } elsif ( ref $what eq 'Regexp' ) {
                    $wrong++ unless $args->{$key} =~ /$what/;

                } elsif ( ref $what eq 'ARRAY' ) {
                    $wrong++ unless grep { ref $_ eq 'Regexp'
                                                ? $args->{$key} =~ /$_/
                                                : $args->{$key} eq $_
                                         } @$what;

                } elsif ( ref $what eq 'CODE' ) {
                    $wrong++ unless $what->( $key => $args->{$key} );

                } else {
                    carp qq[Can not dot to allow checking based on a ] . ref $what;
                }
            }

            if( $STRICT_TYPE || $tmpl->{$key}->{strict_type} ) {
                $wrong++ unless ref $args->{$key} eq ref $tmpl->{$key}->{default};
            }

            ### somehow it's the wrong type.. warn for this! ###
            if( $wrong ) {
                carp qq[Key '$key' is of wrong type for ] . _who_was_it() if $verbose;
                ++$flag && next;

            } else {

                ### if we got here, it's apparently an ok value for $key,
                ### so we'll set it in the default to return it in a bit
                $defs->{$key} = $args->{$key};
            }

        }
    }

    return $flag ? undef : $defs;
}


### Return a hashref of $tmpl keys with required values
sub _listreqs {
    my $tmpl = shift;

    my %hash = map { $_ => 1 } grep { $tmpl->{$_}->{required} } keys %$tmpl;
    return \%hash;
}


### check if the $key is required, and if so, whether it's in $args ###
sub _hasreq {
    my ($key, $tmpl, $args ) = @_;
    my $reqs = _listreqs($tmpl);

    return $reqs->{$key}
            ? exists $args->{$key}
                ? 1
                : undef
            : 1;
}

### Return a hash of $tmpl keys with default values => defaults
### make sure to even include undefined ones, so that 'exists' will dwym
sub _hashdefs {
    my $tmpl = shift;

    my %hash =  map {
                    $_ => defined $tmpl->{$_}->{default}
                                ? $tmpl->{$_}->{default}
                                : undef
                } keys %$tmpl;

    return \%hash;
}

### check if the key exists in $data ###
sub _iskey {
    my ($key, $tmpl) = @_;
    return $tmpl->{$key} ? 1 : undef;
}

sub _who_was_it { return (caller(2))[3] || 'ANON' }

1;

__END__

=pod

=head1 NAME

Tools::Check;

=head1 SYNOPSIS

    use Tools::Check qw[check];

    sub fill_personal_info {
        my %hash = @_;

        my $tmpl = {
            firstname   => { required   => 1, },
            lastname    => { required   => 1, },
            gender      => { required   => 1,
                             allow      => [qr/M/i, qr/F/i],
                           },
            married     => { allow      => [0,1] },
            age         => { default    => 21,
                             allow      => qr/^\d+$/,
                           },
            id_list     => { default    => [],
                             strict_type => 1
                           },
            phone       => { allow      => sub {
                                                my %args = @_;
                                                return 1 if
                                                    &valid_phone(
                                                        $args{phone}
                                                    );
                                            }
                            },
            }
        };

        my $parsed_args = check( $tmpl, \%hash, $VERBOSE )
                            or die [Could not parse arguments!];

=head1 DESCRIPTION

Tools::Check is a generic input parsing/checking mechanism.

It allows you to validate input via a template. The only requirement
is that the arguments must be named.

Tools::Check will do the following things for you:

=over 4

=item *

Convert all keys to lowercase

=item *

Check if all required arguments have been provided

=item *

Set arguments that have not been provided to the default

=item *

Weed out arguments that are not supported and warn about them to the
user

=item *

Validate the arguments given by the user based on strings, regexes,
lists or even subroutines

=item *

Enforce type integrity if required

=back

Most of Tools::Check's power comes from it's template, which we'll
discuss below:



=head1 Template

As you can see in the synopsis, based on your template, the arguments
provided will be validated.

The template can take a different set of rules per key that is used.

The following rules are available:

=over 4

=item default

This is the default value if none was provided by the user.
This is also the type C<strict_type> will look at when checking type
integrity (see below).

=item required

A boolean flag that indicates if this argument was a required
argument. If marked as required and not provided, check() will fail.

=item strict_type

This does a C<ref()> check on the argument provided. The C<ref> of the
argument must be the same as the C<ref> of the default value for this
check to pass.

This is very usefull if you insist on taking an array reference as
argument for example.

=item allow

A set of criteria used to validate a perticular piece of data if it
has to adhere to particular rules.
You can use the following types of values for allow:

=over 4

=item string

The provided argument MUST be equal to the string for the validation
to pass.

=item array ref

The provided argument MUST equal (or match in case of a regular
expression) one of the elements of the array ref for the validation to
pass.

=item regexp

The provided argument MUST match the regular expression for the
validation to pass.

=item subroutine

The provided subroutine MUST return true in order for the validation
to pass and the argument accepted.

(This is particularly usefull for more complicated data).

=back

=back

=head1 Functions

=head2 check

Tools::Check only has one function, which is called C<check>.

This function is not exported by default, so you'll have to ask for it
via:

    use Tools::Check qw[check];

or use it's fully qualified name instead.

C<check> takes a list of arguments, as follows:

=over 4

=item Template

This is a hashreference which contains a template as explained in the
synopsis.

=item Arguments

This is a reference to a hash of named arguments which need checking.

=item Verbose

A boolean to indicate whether C<check> should be verbose and warn
about whant went wrong in a check or not.

=back

C<check> will return undef when it fails, or a hashref with lowercase
keys of parsed arguments when it succeeds.

So a typical call to check would look like this:

    my $parsed = check( \%template, \%arguments, $VERBOSE )
                    or warn q[Arguments could not be parsed!];

=head1 AUTHOR

This module by
Jos Boumans E<lt>kane@cpan.orgE<gt>.

=head1 Acknowledgements

Thanks to Ann Barcomb for her suggestions.

=head1 COPYRIGHT

This module is
copyright (c) 2002 Jos Boumans E<lt>kane@cpan.orgE<gt>.
All rights reserved.

This library is free software;
you may redistribute and/or modify it under the same
terms as Perl itself.

=cut

# Local variables:
# c-indentation-style: bsd
# c-basic-offset: 4
# indent-tabs-mode: nil
# End:
# vim: expandtab shiftwidth=4:
         