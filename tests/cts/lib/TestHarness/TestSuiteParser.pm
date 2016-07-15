package TestHarness::TestSuiteParser;

use Data::Dumper;
use Exporter;
use Carp;
use Parse::RecDescent;
use Text::Balanced qw ( extract_codeblock extract_bracketed extract_quotelike extract_delimited );
use TestHarness::Grammar;

# This module loads a precompiled grammar. To create the grammar
# run the following shell command.

# perl -MParse::RecDescent - testsuite_grammar TestHarness::Grammar


$::RD_WARN = 1;
$::RD_HINT = 1;
#$::RD_TRACE = 1;


sub new{
  my $parser = TestHarness::Grammar->new();
  return $parser;
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

