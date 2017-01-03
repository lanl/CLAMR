package Tools::SQL;

use strict;
use DBI;

push @DBI::db::ISA, __PACKAGE__
        unless grep { $_ eq __PACKAGE__ } @DBI::db::ISA;

use vars qw[$FATAL $VERBOSE];

$FATAL      = 0;
$VERBOSE    = 1;

sub exec_sql {
    my $dbh = shift;
    my $sql = shift;


    my $multiple;
    if( my $type = ref $_[0] ) {
        if( scalar @_ > 1 ) {
            warn qq[Useless arguments trailing: ] .
                    join " ", @_[1..$#_] if $VERBOSE;
        }

        if( $type eq 'ARRAY' ) {
            $multiple++;
        } else {
            warn qq[Invalid argument for sql params: '$type'] if $VERBOSE;
            return;
        }
    }

    my $sth;

    unless( $sth = $dbh->prepare_cached($sql) ) {
        my $msg = qq[Could not prepare query '$sql': ] . $dbh->errstr;

        ($FATAL ? die $msg : warn $msg) if $VERBOSE;

        return;
    }

    my $failed;
    if ( $multiple ) {
        my $aref = shift;

        for my $arg ( @$aref ) {

            if( ref $arg ne 'ARRAY' ) {
                warn qq[Not an arrayref: $arg, skipping...] if $VERBOSE;
                $failed++;
                next;
            }

            $failed++ unless _exec_sql( $sth, $sql, $arg );
        }

        $sth->finish;

        if( $failed ) {
            my $tried = scalar @$aref;
            warn qq[failed a total of $failed queries out of $tried] if $VERBOSE;
        }

    } else {
        $failed++ unless _exec_sql( $sth, $sql, [@_] );
    }

    return  $failed     ? undef :
            $multiple   ? 1     :
            $sth;
}

sub _exec_sql {
    my $sth = shift;
    my $sql = shift;
    my $arg = shift;

    my $rv;
    unless( $rv = $sth->execute( @$arg ) ) {

        my $msg = qq[Could not execute query '$sql' with ] .
              qq[arguments: ] .
              join ' ', @$arg, '=>', $sth->errstr;

        ($FATAL ? die $msg : warn $msg) if $VERBOSE;

        return undef;
    }

    return $rv;
}

1;

__END__

=pod

=head1 NAME

Tools::SQL - sql execution made easy

=head1 SYNOPSIS

    use Tools::SQL;
    use DBI;

    my $dbh = DBI->connect($data_source, $username, $auth, \%attr);


    ### do a multiple insert into a table with a single statement
    my $sql = 'insert into some_table (thing,stuff) values (?,?)';

    my @vals = ( [qw|foo bar|], [qw|zoo quux|] );

    my $bool = $dbh->exec_sql( $sql, \@vals );


    ### execute a query and get the statement handler back on success
    my $sql = 'select * from some_table where thing = ?';

    my $sth = $dbh->exec_sql( $sql, 'something' )
                    or warn qq[Could not execute '$sql'];

    ### now do stuff with $sth
    while( $sth->fetchrow_hashref ) { ... }

    ### don't forget to call finish() when you got a $sth back
    $sth->finish;


    ### don't have Tools::SQL issue warnings -- default is '1'
    $Tools::SQL::VERBOSE = 0;

    ### make errors fatal -- default is '0'
    $Tools::SQL::FATAL = 1;


=head1 DESCRIPTION

Tools::SQL is a transparent way of eliminating the overhead of having
to first prepare a statement, check if it prepared correctly, then
execute a query and check if a if executed correctly.

Simply hand it the SQL you want executed and any arguments you need to
be passed to fill the placeholders and Tools::SQL will just DWYM.

=head1 How it works

Tools::SQL places itself at the back of your databasedrivers @ISA
array, so you can call it's functions through your database handler.

Tools::SQL currently only provides one method, namely C<exec_sql>

=head1 Methods

=head2 exec_sql

C<exec_sql> takes the following arguments:

=over 4

=item *

A string containing the sql you want executed, containing any
placeholders where necessary.

=item *

Either a list or a reference to a list of lists with arguments that
will serve to be substituted for the placeholders.

=back

When providing a list of arguments, C<exec_sql> assumes you only want
to do one insert and will return you a statement handler upon success,
or undef upon failure of the execution.

When you provide a reference to a list of lists, C<exec_sql> assumes
you want to execute the statement multiple times, once for each list.
It will then return true if all executions went well, or undef if one
or more of them failed.

=head1 Global Variables

The behaviour of Tools::SQL can be altered by changing the following
global variables:

=head2 $Tools::SQL::VERBOSE

This controls whether Tools::SQL will issue warnings and explenations
as to why certain things may have failed. If you set it to 0,
Tools::SQL will not output anything.
The default is 1;

=head2 $Tools::SQL::FATAL

This controls whether failures for executing statements via Tools::SQL
should be considered fatal or not. When set to 1, whenever a statement
fails, Tools::SQL will die with it's error message.

The default is 0;

=head1 See Also

C<DBI>

=head1 Todo

Inline::SQL?

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
