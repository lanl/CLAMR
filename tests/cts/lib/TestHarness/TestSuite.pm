package TestHarness::TestSuite;

use strict;
use warnings;
use Carp;
use DumpStack;
use TestHarness::TestSuiteParser;
use Cwd qw(getcwd chdir);
use File::Find;
use File::Basename;
use vars qw($AUTOLOAD @ISA);
use Data::Dumper;
$Data::Dumper::Indent = 2;

my $debug = 0;
my $parser;
my $print_number = 0; # This is for giving PRINT|ERROR|WARNING "TESTS" a unique name.

# Initialize the parser before the rest of this file is even parsed by perl. It
# is a bit misleading to call it a parser. The actions that are taken as it
# parses a testsuite file cause tests and included testsuites to be added to the
# current testsuite. 
BEGIN{
  $parser = TestHarness::TestSuiteParser->new();
}


{    # BEGIN CLOSURE BLOCK: See Perl reference for information on closures.
  my %testsuite = (
		   HARNESS        => undef    ,
        	   NAME           => ""       ,
		   RESULTS        => undef    ,
		   RESULTS_PATH   => undef    ,
		   SCP_RESULTS    => "no"     ,
		   SUITE_FILE     => ""       ,
		   TESTS          => undef    ,
		   TESTSUITES     => undef    ,
		  );

  sub AUTOLOAD {
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);
    
    return if ($method =~ /destroy/i);
    
    $method = uc $method;
    unless (exists $testsuite{$method}) {
      print "$package ::AUTOLOAD::ERROR:: No such attribute: $method: is defined for $package.\n";
      dump_stack;
      local $\ = "\n";
      print "The available attributes are";
      foreach (keys %testsuite) {print}
      croak "Aborting";
    }
    
    my $code = q{ sub {
      my $self = shift;
      my $value = shift;

      if (defined $value) {
	if (ref($value) eq "ARRAY") {
	  unshift @{$self->{METHOD}}, @$value;	
	}
	else {
	  $self->{METHOD} = $value;
	}
      }

      my $val = $self->{METHOD} || undef;
      return $val;
    }};
    
    $code =~ s/METHOD/$method/g;
    
    no strict "refs";
    *$AUTOLOAD = eval $code;
    
    goto &$AUTOLOAD;
  }


  sub new {
    my ($invocant, $arg_ref) = @_;
    my $class = ref($invocant) || $invocant;
    my $self = {};
    bless $self => $class;

    $arg_ref ||= {};

    #print "Debug_testsuite: ", Dumper($arg_ref);
    #my ($calling_package, $filename, $line) = caller;
    #print "testsuite::new:: $calling_package, $filename, $line\n";
    #print Dumper($arg_ref);
    my $href = ();
    %$href = (%testsuite, %$arg_ref);
    while (my ($func, $value) = each %$href) {
      $self->$func($value);
    }

    $self->harness(TestHarness::Harness->new);

    my $file = $self->name
      || croak "TestHarness::TestSuite::new: No testsuite was provided";


   # parse the testsuite file
    $self->_parse_testsuite_file($file);

    return $self;
  }

}

# _parse_testsuite_file should only be called internally.
sub _parse_testsuite_file {
  my $self = shift;
  my $testsuite_file = shift 
    or croak ("TestHarness::TestSuite::parse_testsuite_file: No testsuite file given.");

  # name ourself
  $self->name(basename $testsuite_file);

 # find the testsuite file
  unless (-f $testsuite_file) {
    find (sub {$testsuite_file = $File::Find::name if /\b$testsuite_file\Z/}, $self->harness->get_testsuite_directories()); 
  }

  croak "TestHarness::TestSuite::_parse_testsuite_file: Could not find testsuite $testsuite_file" 
    unless (-f $testsuite_file);

  open my $fh, '<', $testsuite_file 
    or croak ("TestHarness::TestSuite::parse_testsuite_file: Could not open $testsuite_file: $!");
  $self->suite_file($testsuite_file);

  # slurp in the entire testsuite file
  my $text;
  {
    local $/;
    $text = <$fh>;
    parse_suite( \$text );
  }

  # Close the file
  close $fh 
    or croak "TestHarness::TestSuite::parse_testsuite_file: Couldn't close file $testsuite_file";

  # Use the parser created in the begin block to parse the testsuite file and
  # fill in the data structure. This includes adding tests and testsuites.
  defined ($parser->suitefile($text, 1, $self))
    or croak "TestHarness::TestSuite::parse_testsuite_file: Failed to parse $testsuite_file";
}
# parse suite file - initial parse to do if-else blocks
# not going through grammar since lmdm could not figure out what to do
sub parse_suite
  {
    my(
       $text_ref
      ) = @_;
    my(
       $done, # if done
       $i, # loop var
       $j, # loop var
       @lines, # split on whitespace for $text_ref
       $level, # what level of if statement you are in
       $line, # single line
       @results, # the results of the if block
       @results_true, # if result for current level if block ever been true
       $result, # final result to print line or not
       $text_new, # new text
      );
    @lines = split( /\n/, $$text_ref );
    $done = "false";
    $text_new = "";
    $i = 0;
    $level = 0;
    $results[$level] = 1;
    while( $done eq "false" )
      {
        if( $i > $#lines )
          {
            $done = "true";
            next;
          }
        $line = "$lines[$i]\n";
        if( $line =~ /^\s*if\s*(\(.*\))\s*then\s*$/ )
          {
            $level++;
            $results_true[$level] = !1;
            $results[$level] = eval( $1 );
            $text_new .= "# $line";
            $i++;
            next;
          }
        if( $line =~ /^\s*else\s*if\s*(\(.*\))\s*then\s*$/ )
          {
            $results[$level] = eval( $1 );
            $results[$level] = $results[$level] && ( ! $results_true[$level] );
            $text_new .= "# $line";
            $i++;
            next;
          }
        if( $line =~ /^\s*else\s*$/ )
          {
            $results[$level] = ! $results_true[$level];
            $text_new .= "# $line";
            $i++;
            next;
          }
        if( $line =~ /^\s*endif\s*$/ )
          {
            $level--;
            $text_new .= "# $line";
            $i++;
            next;
          }
        $result = $results[0];
        for( $j = 1; $j <= $level; $j++ )
          {
            $result = $result && $results[$j];
          }
        if( $result )
          {
            $results_true[$level] = !0;
            $text_new .= $line;
          }
        else
          {
            $text_new .= "# $line";
          }
        $i++;
      }
    $$text_ref = $text_new;
    return;
  }

sub suite_file{
  my $self = shift;
  my $suite_file = shift || $self->{SUITE_FILE} || "Unknown";
  if (ref $suite_file) {
    $suite_file = $$suite_file[0];
  }
  $self->{SUITE_FILE} = $suite_file;
  return $suite_file;
}


sub add_test{
  my $self = shift;
  my $test_name = shift || croak "TestHarness::TestSuite::add_test: No test problem was provided.";
  my $print_string = shift || "";

  print "Adding $test_name\n" if defined $self->harness->verbose;
  unless (exists $self->{TESTS}{$test_name}) {
    $self->{TESTS}{$test_name} = $self->harness->add_test($test_name, $print_string);
  
  }

  return $self->{TESTS}{$test_name};
}


sub add_testsuite{
  my $self = shift;
  my $suite_name = shift || croak "TestHarness::TestSuite::add_testsuite: No testsuite was provided.";

  print "Adding $suite_name\n" if defined $self->harness->verbose;
  unless (exists $self->{TESTSUITES}{$suite_name}) {
    $self->{TESTSUITES}{$suite_name} = $self->harness->add_testsuite($suite_name);
  }

  foreach my $test(keys %{$self->{TESTSUITES}{$suite_name}{TESTS}}) {
    $self->add_test($test);
  }

  return $self->{TESTSUITES}{$suite_name};
}




sub add_judge{
  my $self = shift;
  my $judge = shift or croak "TestHarness::TestProblem:add_judge: No judge to add.";

  # Add judge to all included test problems to the current suite
  foreach my $test ( values %{$self->{TEST_PROBLEMS}}) {
    $test->add_judge($judge);
  }
}

sub add_anti_judge{
  my $self = shift;
  my $judge = shift or croak "TestHarness::TestProblem:add_judge: No judge to add.";

  # Add anti-judge to all included test problems to the current suite
  foreach my $test ( values %{$self->{TEST_PROBLEMS}}) {
    $test->add_anti_judge($judge);
  }
}


sub extract_results{
  my $self = shift;
  my $test = shift;
  my $test_name = $test->name;
  #print "\nextract_results1: $prob_name\n";
 
#   carp "TestHarness::Harness::extract_results: No problem name provided" 
#       unless $prob_name;
 
  my $results = shift || $test->extract_results;
  #print "\nextract_results2: ", Dumper($results), "\n\n";
  $self->{RESULTS}{$test_name}=$results;

  my $suites = $self->testsuites;
  if (defined $suites) {
    foreach my $suite(values %$suites) {
      if (exists $suite->{TESTS}{$test_name}) {
	$suite->extract_results($test, $results);
      }
    }
  }
}
  


sub results{
  my $self = shift;
  my $results_hash = $self->{RESULTS};
  return $results_hash;
}


sub scp_results{
  my $self = shift;

  my $scp_results = shift || $self->{SCP_RESULTS};
  $scp_results = (ref $scp_results)? @$scp_results[0]: $scp_results;
  $self->{SCP_RESULTS} = $scp_results;
  return $scp_results;
}


sub results_path{
  my $self = shift;

  my $results_path = shift || $self->{RESULTS_PATH} || "./";
  $results_path = (ref $results_path)? @$results_path[0]: $results_path;
  $self->{RESULTS_PATH} = $results_path;
  return $results_path;
}


sub html_file{
  my $self = shift;

  my $name = $self->name;
  $self->{HTML_FILE} = shift if @_;
  my $val = $self->{HTML_FILE} || $name . "_results.html" ;
  return $val;
}

sub text_file{
  my $self = shift;

  my $name = $self->name;
  $self->{TEXT_FILE} = shift if @_;
  my $val = $self->{TEXT_FILE} || $name . "_results.txt" ;
  return $val;
}


sub dump {
  my $self = shift;
  print "TESTSUITE:", Dumper($self);
}



1;
__END__


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

