package Tools::Term;

use Tools::Check qw[check];
use Term::ReadLine;

use strict;

BEGIN {
    use vars        qw[$VERSION $AUTOREPLY $VERBOSE $INVALID];
    $VERBOSE    =   1;
    $INVALID    =   'Invalid selection, please try again: ';
}

push @Term::ReadLine::Stub::ISA, __PACKAGE__
        unless grep { $_ eq __PACKAGE__ } @Term::ReadLine::Stub::ISA;


sub get_reply {
    my $term = shift;
    my %hash = @_;

    my $tmpl = {
        default     => { default => undef,  strict_type => 1 },
        prompt      => { default => '',     strict_type => 1, required => 1 },
        choices     => { default => [],     strict_type => 1 },
        multi       => { default => 0,      allow => [0, 1] },
        prompt_add  => { default => '',     no_override => 1 },
        allow       => { default => qr/.*/ },
    };

    my $args = check( $tmpl, \%hash, $VERBOSE )
                or ( warn( qq[Could not parse arguments] ), return );

    my $prompt_add;

    if( @{$args->{choices}} ) {
        my $i;

        $args->{prompt} =
            join "\n", map( {$i++;
                             $args->{prompt_add} = $i if $_ eq $args->{default};
                             sprintf "%3s> %-s", $i, $_ ;
                       } @{$args->{choices}} ), $args->{prompt};

        $args->{allow} = 1;

    } elsif ( defined $args->{default} ) {
        $args->{prompt_add} = $args->{default};
    }

    $args->{prompt} .= $args->{prompt_add}
                                ? ' ['. $args->{prompt_add} . ']: '
                                : ' ';

    return _tt_readline( $term, %$args );

} #_get_reply

sub ask_yn {
    my $term = shift;
    my %hash = @_;

    my $tmpl = {
        default     => { default => undef,                  strict_type => 1 },
        prompt      => { default => '', required => 1,      strict_type => 1 },
        multi       => { default => 0,                      no_override => 1 },
        choices     => { default => [qw|y n|],              no_override => 1 },
        prompt_add  => { default => '',                     no_override => 1 },
        allow       => { default => [qr/^y(?:es)?$/i, qr/^n(?:o)?$/],
                         no_override => 1
                       },
    };

    my $args = check( $tmpl, \%hash, $VERBOSE ) or return undef;

    my @list = @{$args->{choices}};
    if( $args->{default} ) {
        @list = map { lc $args->{default} eq lc $_
                            ? uc $args->{default}
                            : $_
                } @list;
    }

    $args->{prompt_add} .= join("/", @list);
    $args->{prompt}     .= ' ['. $args->{prompt_add} . ']: ';

    return _tt_readline( $term, %$args ) =~ /^y/i ? 1 : 0;
}



sub _tt_readline {
    my $term = shift;
    my %hash = @_;

    my $tmpl = {
        default     => { default => undef,  strict_type => 1 },
        prompt      => { default => '',     strict_type => 1, required => 1 },
        choices     => { default => [],     strict_type => 1 },
        multi       => { default => 0,      allow => [0, 1] },
        allow       => { default => qr/.*/ },
        prompt_add  => { default => '' },
    };

    my $args = check( $tmpl, \%hash, $VERBOSE ) or return undef;

    if ($AUTOREPLY) {
        warn q[You have $AUTOREPLY set to true, but did not provide a default!]
                if( !defined $args->{default} && $VERBOSE);

        print "$args->{prompt}\n";
        return $args->{default}
    }

    my @list = @{$args->{choices}};
    LOOP: {
        my $answer  = $term->readline( $args->{prompt} );
        $answer     = $args->{default} unless length $answer;

        $term->addhistory( $answer ) if length $answer and !$AUTOREPLY;

        my @answers = $args->{multi} ? split(/\s+/, $answer) : $answer;

        my $flag = 0;
        my @rv;
        if( @list ) {
            @rv =   grep {
                        defined
                    } map {
                        length $_ ? /\D/
                            ? check( { tmp => {allow => $args->{allow}} },
                                     { tmp => $_ } )
                                ? $_
                                : undef
                            : $_ > 0 && defined $list[$_-1]
                                ? $list[$_-1]
                                : undef
                            : undef;
                    } @answers;
        } else {
            @rv = grep {
                        check(  { tmp => { allow => $args->{allow} } },
                                { tmp => $_ }
                        )
                    } scalar @answers ? @answers : ($args->{default});
        }

        unless( @rv == @answers ) {
            $args->{prompt} = $INVALID;
            $args->{prompt} .= "[$args->{prompt_add}] " if $args->{prompt_add};
            redo LOOP;

        } else {
            return $args->{multi} ? @rv : $rv[0];
        }
    }
}

sub parse_options {
    my $term    = shift;
    my $input   = shift;

    my $return = {};

    ### there's probably a more elegant way to do this... ###
    while ( $input =~ s/--?([-\w]+=("|').+?\2)\s*//   or
            $input =~ s/--?([-\w]+=\S+)\s*//          or
            $input =~ s/--?([-\w]+)\s*//
    ) {
        my $match = $1;

        if( $match =~ /^([-\w]+)=("|')(.+?)\2$/ ) {
            $return->{$1} = $3;

        } elsif( $match =~ /^([-\w]+)=(\S+)$/ ) {
            $return->{$1} = $2;

        } elsif( $match =~ /^no-?([-\w]+)$/i ) {
            $return->{$1} = 0;

        } elsif ( $match =~ /^([-\w]+)$/ ) {
            $return->{$1} = 1;

        } else {
            warn qq[I do not understand option "$match"\n] if $VERBOSE;
        }
    }

    return wantarray ? ($return,$input) : $return;
}


1;

__END__

=pod

=head1 NAME

Tools::Term - Term::ReadLine UI made easy

=head1 SYNOPSIS

    use Tools::Term;
    use Term::ReadLine;

    my $term = Term::ReadLine->new('brand');

    my $reply = $term->get_reply(
                    prompt => 'What is your favourite colour?',
                    choices => [qw|blue red green|],
                    default => blue,
    );

    my $bool = $term->ask_yn(
                        prompt => 'Do you like cookies?',
                        default => 'y',
                );


    my $string = q[some_command -option --no-foo --quux='this thing'];

    my ($options,$munged_input) = $term->parse_options($string);


    ### don't have Tools::Term issue warnings -- default is '1'
    $Tools::SQL::VERBOSE = 0;

    ### always pick the default (good for non-interactive terms)
    ### -- default is '0'
    $Tools::SQL::AUTOREPLY = 1;


=head1 DESCRIPTION

Tools::Term is a transparent way of eliminating the overhead of having
to format a question and then validate the reply, informing the user
if the answer was not proper and re-issuing the question.

Simply give it the question you want to ask, optionally with choices
the user can pick from and a default and Tools::Term will DWYM.

For asking a yes or no question, there's even a shortcut.

=head1 How it works

Tools::Term places itself at the back of the Term::ReadLine @ISA
array, so you can call it's functions through your term object.

=head1 Methods

=head2 get_reply

C<get_reply> takes the following arguments:

=over 4

=item prompt

A string containing the question you want to ask,

=item choices

An array reference with all the choices you want to let the user pick
from. This parameter is optional.

=item default

The default answer -- This is the answer picked if the user just hits
C<enter> or if $AUTOREPLY is set to true. This parameter is optional.

=item allow

A handler you can specify to check the reply against. This can be any
of the types that the Tools::Check C<check> function allows, so please
refer to that manpage for details. This parameter is optional.

=back

=head2 ask_yn

C<ask_yn> takes the following arguments:

=over 4

=item prompt

A string containing the question you want to ask,

=item default

The default answer -- This is the answer picked if the user just hits
C<enter> or if $AUTOREPLY is set to true. This paramter is optional.

=back

The choices that are presented are either C<yes> or C<no> and will
return 0 if C<no> was answered and 1 if C<yes> was answered.

=head2 parse_options

C<parse_options> will convert all options given from an input string
to a hash reference. If called in list context it will also return
the part of the input string that it found no options in.

Consider this example:

    my $str =   q[command --no-foo --baz --bar=0 --quux=bleh ] .
                q[--option="some'thing" -one-dash -single=blah' arg];

    my ($options,$munged) =  $term->parse_options($str);

    ### $options would contain: ###
    $options = {
                'foo'       => 0,
                'bar'       => 0,
                'one-dash'  => 1,
                'baz'       => 1,
                'quux'      => 'bleh',
                'single'    => 'blah\'',
                'option'    => 'some\'thing'
    };

    ### and this is the munged version of the input string,
    ### ie what's left of the input minus the options
    $munged = 'command arg';

As you can see, you can either use a single or a double C<-> to
indicate an option.
If you prefix an option with C<no-> and do not give it a value, it
will be set to 0.
If it has no prefix and no value, it will be set to 1.
Otherwise, it will be set to it's value. Note also that it can deal
fine with single/double quoting issues.

=head1 Global Variables

The behaviour of Tools::Term can be altered by changing the following
global variables:

=head2 $Tools::Term::VERBOSE

This controls whether Tools::Term will issue warnings and explenations
as to why certain things may have failed. If you set it to 0,
Tools::Term will not output any warnings.
The default is 1;

=head2 $Tools::Term::AUTOREPLY

This will make every question be answered by the default, and warn if
there was no default provided. This is particularly usefull if your
program is run in non-interactive mode.
The default is 0;

=head2 $Tools::Term::INVALID

This holds the string that will be printed when the user makes an
invalid choice.
You can override this string from your program if you, for example,
wish to do localization.
The default is C<Invalid selection, please try again: >

=head1 See Also

C<Tools::Check>, C<Term::ReadLine>

=head1 AUTHOR

This module by
Jos Boumans E<lt>kane@cpan.orgE<gt>.

=head1 COPYRIGHT

This module is
copyright (c) 2002 Jos Boumans E<lt>kane@cpan.orgE<gt>.
All rights reserved.

This library is free software;
you may redistribute and/or modify it under the same
terms as Perl itself.