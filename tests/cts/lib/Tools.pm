package Tools;

1;

__END__

=pod

=head1 NAME

Tools::* - common utilities for common problems

=head1 DESCRIPTION

The Tools::* modules are an attempt to abstract out common code for
common problems. For example, it provides code to enable generic
message logging and handling from programs, simple UI construction
using Term::ReadLine, input parsing/validating, loading of files and
so on.

See below for a more elaborate package description.

=head1 MODULES

=over 4

=item Tools::Load

Allows for generic loading of modules and files. Simply give it the
name of a module or file and it will Do What You Mean (tm).

=item Tools::Check

Allows for generic input checking and validating using a powerfull
templating system

=item Tools::Cmd

Allows for the searching and execution of any binary on your system.
It adheres to verbosity settings and is able to run intereactive.
It also has an option to capture output/error buffers.

=item Tools::Term

Provides methods to ask both elaborate questions as well as simple
yes/no questions via a Term::ReadLine interface using a template.
It can also parse options per unix style.

=item Tools::Module

Allows you to query the state of modules on your system. It can tell
you if you have certain modules installed without attempting to C<use>
them and can do smart loading of modules.
Also it can tell you what *other* modules a certain module requires.

=item Tools::SQL

Allows for simple execution of SQL statements, eliminating the
overhead of preparing and executing the queries. Simply provide the
statement and the value(s) that need to be inserted and Tools::SQL
will Do What You Mean (tm).

=item Tools::Log

This module enables you to do generic message logging throughout
programs and projects. Every message will be logged with stacktraces,
timestamps and so on. You can use built-in handlers immediately, or
after the fact when you inspect the error stack.
It is highly configurable and let's you even provide your own
handlers for dealing with messages.

=back

=head1 AUTHOR

This module by
Jos Boumans E<lt>kane@cpan.orgE<gt>.

=head1 Acknowledgements

Thanks to Ann Barcomb for her suggestions and my employer C<XS4ALL>
for letting me develop and release this during company hours.

=head1 COPYRIGHT

This module is
copyright (c) 2002 Jos Boumans E<lt>kane@cpan.orgE<gt>.
All rights reserved.

This library is free software;
you may redistribute and/or modify it under the same
terms as Perl itself.

=cut

