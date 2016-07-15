package Data::Show;

use warnings;
use strict;
use Data::Dump 'dump';
use 5.010;

our $VERSION = '0.002001';

# Unconditionally export show()...
sub import {
    no strict 'refs';
    *{caller().'::show'} = \&show;
}


# Useful patterns...

my $IDENT = qr{
    [^\W\d]\w* (?: :: [^\W\d]\w* )* | [_\W]
}xms;

use re 'eval';
my $CODE = qr{
        (?&CODE_FRAGMENT)

        (?(DEFINE)
            (?<CODE_FRAGMENT>
                (?: (?&QUOTED)
                  | \b (?: q[qxr]?+ | [msy] | tr ) \s* (?&DELIMITED)
                  | (?&NESTED)
                  | [%\$\@] (?: (?&BRACE_DELIMS) | $IDENT) (?&NESTED)?
                  | [^][{}()"'`;]
                )++
            )

            (?<NESTED_CODE_FRAGMENT>
                (?: (?&QUOTED)
                  | \b (?: q[qxr]?+ | [msy] | tr ) \s* (?&DELIMITED)
                  | (?&NESTED)
                  | [^][{}()"'`]
                )++
            )

            (?<QUOTED>
                   " [^\\"]++ (?: \\. [^\\"]++ )*  "
                |  ' [^\\']++ (?: \\. [^\\']++ )*  '
                |  ` [^\\`]++ (?: \\. [^\\`]++ )*  `
                |  / [^\\/]++ (?: \\. [^\\/]++ )*  /
                | \? [^\\?]++ (?: \\. [^\\?]++ )* \?
            )

            (?<DELIMITED>
                  (?&BRACE_DELIMS)
                | (?&PAREN_DELIMS)
                | (?&ANGLE_DELIMS)
                | (?&SQUARE_DELIMS)
                | \s++ (?<DELIM_W>\w) (?:\\.|(?!\g{DELIM_W}).)*+ \g{DELIM_W}
                |      (?<DELIM_S>[^\w\s]) (?:\\.|(?!\g{DELIM_S}).)*+ \g{DELIM_S}
            )

            (?<NESTED>
                  \( (?&NESTED_CODE_FRAGMENT) \)
                | \[ (?&NESTED_CODE_FRAGMENT) \]
                | \{ (?&NESTED_CODE_FRAGMENT) \}
                | \< (?&NESTED_CODE_FRAGMENT) \>
            )

            (?<BRACE_DELIMS>   \{ (?: [^{}] | \\. | (?&BRACE_DELIMS)  )*+ \} )
            (?<PAREN_DELIMS>   \( (?: [^()] | \\. | (?&PAREN_DELIMS)  )*+ \) )
            (?<ANGLE_DELIMS>   \< (?: [^<>] | \\. | (?&ANGLE_DELIMS)  )*+ \> )
            (?<SQUARE_DELIMS>  \[ (?: [^][] | \\. | (?&SQUARE_DELIMS) )*+ \] )
        )
}xms;


# Configuration for layout of representation...
my $DEFAULT_INDENT = 4;
my $MAX_WIDTH      = 72;
my $TITLE_POS      = 6;

# The whole point of the module...
sub show {

    # Determine context...
    my (undef, $file, $line) = caller();

    # Extract description of arglist from context...
    my $desc;
    if (open my $fh, '<', $file) {
        for (1..$line-1) { readline($fh) // last }
        $desc = do { local $/; readline $fh; };
    }

    # Trim filename and format context info and description...
    $file =~ s{.*/}{}xms;
    my $context = "'$file', line $line";
    $desc //= $context;

    # Isolate arg list and compress internal whitespace...
    $desc =~ s{ \A (?: (?!\bshow) . )*? \b show \b \s* ($CODE) \s* (?: [;\}] .* | \Z ) }{$1}xms;
    $desc =~ s{\s+}{ }gxms;

    # Serialize Contextual::Return::Value objects (which break dump())...
    for my $arg (@_) {
        if ((ref($arg)||q{}) =~ m{\A Contextual::Return}xms) {
            require Contextual::Return;
            Contextual::Return::FREEZE($arg);
        }
    }

    # Serialize argument (restoring it, if it was inappropriately flattened)...
    my $representation;
    given ($desc) {
        when (m{ \A \@ $IDENT \s* \Z }xms)                 { $representation = dump \@_; }
        when (m{ \A \(? \s* \% $IDENT \s* \)? \s* \Z }xms) { $representation = dump {@_};}
        default                                            { $representation = dump @_;  }
    }

    # Indent representation wrt heading...
    $representation =~ s{^}{ q{ } x $DEFAULT_INDENT }gxmse;

    # Clean up parens around title...
    $desc =~ s{ \A \s* (?| \( \s* (.*?) \s* \) | \s* (.*?) \s* ) \Z }
              {(  $1  )}xms;

    # Insert title into header...
    my $header = '=' x $MAX_WIDTH;
    substr($header, $TITLE_POS, length($desc), $desc);

    # Add context if it isn't just context...
    if ($desc ne "(  $context  )") {
        $context = "[ $context ]";
        substr($header, -length($context)-$TITLE_POS, -$TITLE_POS, $context);
    }

    # Display data...
    print {*STDERR} "$header\n\n$representation\n\n\n";

    return;
}


1; # Magic true value required at end of module
__END__

=head1 NAME

Data::Show - Dump data structures with name and point-of-origin


=head1 VERSION

This document describes Data::Show version 0.002001


=head1 SYNOPSIS

    use Data::Show;

    show %foo;
    show @bar;
    show (
        @bar,
        $baz,
    );
    show $baz;
    show $ref;
    show @bar[do{1..2;}];
    show 2*3;
    show 'a+b';
    show 100 * sqrt length $baz;
    show $foo{q[;{{{]};


=head1 DESCRIPTION

This module provides a simple wrapper around the Data::Dump module.

A call to C<show> data-dumps its arguments, prefaced
by a divider line that reports the arguments and the
file and line from which C<show()> was called.

For example, the code in the L<SYNOPSIS> would produce
something like:

    ======(  %foo  )========================[ 'demo.pl', line 11 ]======

        { foo => 1, food => 'b', foon => [3, 4, 5, 6, 7, 8, 9, 10] }


    ======(  @bar  )========================[ 'demo.pl', line 12 ]======

        ["b", "a", "r"]


    ======(  @bar, $baz,  )=================[ 'demo.pl', line 13 ]======

        ("b", "a", "r", "baz")


    ======(  $baz  )========================[ 'demo.pl', line 17 ]======

        "baz"


    ======(  $ref  )========================[ 'demo.pl', line 18 ]======

        ["b", "a", "r"]


    ======(  @bar[do{1..2;}]  )=============[ 'demo.pl', line 19 ]======

        ("a", "r")


    ======(  2*3  )=========================[ 'demo.pl', line 20 ]======

        6


    ======(  'a+b'  )=======================[ 'demo.pl', line 21 ]======

        "a+b"


    ======(  100 * sqrt length $baz  )======[ 'demo.pl', line 22 ]======

        "173.205080756888"


    ======(  $foo{q[;{{{]}  )===============[ 'demo.pl', line 23 ]======

        undef


    ======(  'foo' ~~ m/;{\/{/  )===========[ 'demo.pl', line 25 ]======

        ""


=head1 INTERFACE

=over

=item C<show()>

This is the only subroutine provided by the module.
It is always exported.

C<show()> can be called with any number of arguments and data-dumps them
all with a suitable header indicating the arguments' names, and the file
and line from which C<show()> was called.

C<show()> does not return a useful value.

=back


=head1 DIAGNOSTICS

None, apart from those provided by Data::Dump;


=head1 CONFIGURATION AND ENVIRONMENT

Data::Show requires no configuration files or environment variables.


=head1 DEPENDENCIES

Only works under Perl 5.10 and later.

Requires the Data::Dump module.


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

Uses sophisticated regexes to parse out the argument list from the
source. Hence subject to the usual limitations of this technique
(namely, that it may get the argument list wrong occasionally).

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-data-show@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.


=head1 AUTHOR

Damian Conway  C<< <DCONWAY@CPAN.org> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010, Damian Conway C<< <DCONWAY@CPAN.org> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

