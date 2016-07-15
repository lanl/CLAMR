package TestHarness::Harness;

use strict;
use warnings;
use Cwd qw(getcwd chdir);
use Carp;
use DumpStack;
use FileHandle;
use File::Find;
use Text::Balanced qw (extract_quotelike extract_tagged);
use File::Basename;
use vars qw($AUTOLOAD $VERSION @ISA);
$VERSION    =   0.01;
use System;
use Data::Dumper;

$Data::Dumper::Indent = 2;


{   # BEGIN CLOSURE BLOCK: See Perl reference for information on closures.
  # Harness is a singleton.
  my $harness = undef;

  my %project = (
		   CANDIDATE_TEST_DIRECTORIES => undef,
		   CTS_FILE                   => "",
		   CWD                        => getcwd,
		   DB_FILE                    => "",
		   DATA_LINKS                 => undef,
		   DEFAULT_NUM_CPUS           => 1,
                   DEFAULT_NUM_THREADS        => undef,
		   DEFAULT_RESTART_SCRIPT     => "",
		   DEFAULT_TESTSUITES         => undef,
		   DEFAULT_TIME_LIMIT         => "1:00",
		   EXECUTABLES                => {},
		   FAST                       => 0,
		   GOLD_STANDARDS_DIRECTORIES => undef,
		   LIB_DIRS                   => [".", "$FindBin::Bin/../lib", "$FindBin::Bin/lib"],
		   LOG_FILE_SUFFIX            => ".out",
        	   NAME                       => "cts"       ,
        	   NOOP                       => undef,
		   PARALLEL_EXECUTABLES       => {},
                   PATH                       => '',
		   PROJECT                    => "UNKNOWN",
		   REPORTER                   => undef,
		   RERUN                      => undef,
		   RESULTS                    => undef,
		   RESULTS_PATH               => undef,
		   SCP_RESULTS                => "no",
		   SERIAL_EXECUTABLES         => {},
		   SOURCE                     => undef,
		   SYSTEM                     => undef,
		   TESTING_DIR                => getcwd . "/testing/",
		   TESTS                      => {},
		   TESTSUITES                 => {},
		   TESTSUITE_DIRECTORIES      => undef,
		   TEST_DIRECTORIES           => undef,
		   VERBOSE                    => undef,
		);


=item AUTOLOAD: Automatically generate accessor methods not explicitly specified.

=cut
  sub AUTOLOAD {
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);
    
    return if ($method =~ /destroy/i);
    
    $method = uc $method;
    #dump_stack("\nDEBUG: method $method\n");

    if (exists $project{$method}) {
      my $code = q{ sub {
         my $self = shift;
         my $value = shift;

         if (defined $value) {
	   if (ref($value) eq "ARRAY") {
	     push @{$self->{METHOD}}, @$value;	
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

    else {
      my $func = lc $method;
      my $self = shift;
      my $value = shift;

      if ($harness->system->can($func)) {
	return $harness->system->$func($value);
      }
      elsif ($harness->reporter->can($func)) {
	return $harness->reporter->$func($value);
      }
      else {
	print "$package ::AUTOLOAD::ERROR:: No such attribute: $method: is defined for $package.\n";
	dump_stack;
	local $\ = "\n";
	print "The available attributes are";
	foreach (keys %project) {print}
	return undef;
      }
    }
    
  }

# TestHarness::Harness->new(cts_file);
  sub new {
    unless (defined $harness) {
      my ($class, $file, $arg_ref) = @_;
      $harness = {};
      bless $harness=>$class;

      $file ||=  $ENV{CTS_FILE};
#	|| croak "TestHarness::Harness::new: No cts file is specified nor is CTS_FILE set in user's environment.";
      $arg_ref ||= {};
      

      $harness->system(System->new(
				   {
				   }
				  ));
      $harness->reporter(Reporter->new);
      my $href = ();
      %$href = (%project, %$arg_ref);

      while (my ($func, $value) = each %$href) {
	$harness->$func($value);
      }
      
      if (defined($file)){
        $harness->parse_cts_file($file);
      }
      #     print "Dump of test harness\n";
      #     $self->dump;
      
      # Source in any tcsh files listed on SOURCE lines;
      my $platform = ref $harness->system;
      #   print "DEBUG platform is $platform\n";
      $platform =~ s/.*:://;
      $harness->platform($platform);
      foreach my $file (@{$harness->{SOURCE}}) {
	#print $file, "\n";
	$harness->system->source($file);
      }
    }
    return $harness;
  }


  sub parse_cts_file{
    my ($self, $cts_file) = @_;

    my $fh = FileHandle->new($cts_file, "r") 
      or croak "TestHarness::Harness::new: Could not open cts file $cts_file : $!";
    if ($cts_file !~ /^\//) {
      my $local_dir = getcwd . "/";
      $cts_file = $local_dir . $cts_file;
    }
    $self->cts_file($cts_file);  # save the cts filename

    my $verbose = 0;
    my $in_if_block = "";
    while (<$fh>) {
      chomp;

      # Capture XML/HTML content
      my ($html, $remainder) = ("","");
      if (/^([^<]+)(<.*$)/) {
	my $prefix = $1;
	my $postfix = $2;
	($html, $remainder) = extract_tagged $postfix;
	$_ = $prefix . " <HTML> " . $remainder;
      }
      # remove comments
      s/\# .* $//x;

      # Don't bother with blank lines
      next if /^\s*$/;

      if (/^\s*verbose\s*/ix) {
         $verbose = 1;
         next;
      }

      # check for set(env)? VAR =? value - this allows environmental
      # variables to be set in .cts files and accessed anywhere, even
      # in .test files. Currently this will override the user's
      # environment. This behavior could easily be changed if users desired.
      if (/^\s* set(?:env)? \s+ ([a-zA-Z_0-9]+) \s* =? \s* (.*)$/ix) {
	$ENV{uc $1} = $2;
	next;
      }

      if (/^\s*endif\s*/ix) {;
        $in_if_block = "";
        if ($verbose) {print " Parsing cts file: $in_if_block $_\n";}
        next;
      }

      # CHECK FOR IF BLOCK OF THE FORM "if ( `exec` == /pattern/ ) then"
      if ( /^\s*if\s*\(\s*([\`\w -]+)\s*([\=\<\>]+)\s*([-\/\^\[\]\?\|\\\w\d]+)\s*\)\s*then/ix ) {
        my $test_string = eval "$1";
        my $test_op = $2;
        my $test_pattern = $3;
        $test_pattern =~ s!/!!g;
        my $result = $test_string =~ /$test_pattern/;
        if ($test_op eq "==") {
	  if (! $result){
	    while (my $if_block_line = <$fh>) {
	      # print "Throwing away $if_block_line";
	      last if ($if_block_line =~ /^\s*endif\s*/ix)
	    }
	  }
	  else {
	    if ($verbose) {print " Parsing cts file: $in_if_block $_\n";}
	    $in_if_block = "     ";
	  }
        }

      }
      else {

	# Parse non-system line. It should look like this "keyword [:|\s|=] value(s)" 

	# expand environmental variables
	while (/\$\{? ([A-Z_0-9]+) \}?/xg) {
	  s/\$\{? ([A-Z_0-9]+) \}?/$ENV{$1}/x if defined $ENV{$1};
	}
	#s/\$\{? ([A-Z_0-9]+) \}?/$ENV{$1}/xg;


	if (/^\s* (\w+)\b \s*[:=]? \s* (.*)$/x) {
	  my $func = lc $1;
	  my $arg = $2;
	  
	  # Remove () from beginning and end. We assume that there are no ()s in
	  # the values.
	  $arg =~ s/[()]//g;
	  # get the argument list. It should be a single argument or a comma
	  # or blank separated list.
	  my @args = split /[\s,=]+/, $arg;
	  
	  @args = map {s/<HTML>/$html/;$_} @args if $html;
	  # Add the arguments to the function list.
	  eval {$self->$func(\@args)};
          $func = uc $func;
          if ($verbose) {print " Parsing cts file: $in_if_block $func : $arg\n";}
	  if ($@) {croak "\nERROR: Invalid line in $cts_file: \n\n $_\n\n";}
	} else {
	  croak "TestHarness::Harness.pm: Unable to parse $cts_file:\n\n $_\n\n";
	}
      }
    }
    close($fh);

  }


  # Initialize local $harness variable when rerunning cts. This
  # prevents the creation of a new instance of the testharness when
  # HARNESS->new is called from a testsuite or testproblem. 
  sub initialize {
    my $self = shift;
    $harness = $self;

    # Source in any tcsh files listed on SOURCE lines;
    my $platform = ref $harness->system;
    #   print "DEBUG platform is $platform\n";
    $platform =~ s/.*:://;
    $harness->platform($platform);
    foreach my $file (@{$harness->{SOURCE}}) {
      #print $file, "\n";
      $harness->system->source($file);
    }
  }

} # END OF CLOSURE BRACKETS


=item add_test: $myobj->add_test(test|test_name, print_string);

=cut
sub add_test{
  my $self = shift;
  my $test_name = shift || croak "TestHarness::Harness::add_test: No test problem was provided.";
  my $print_string = shift || "";

  print "Adding $test_name\n" if defined $self->verbose;
  unless (exists $self->{TESTS}{$test_name}) {
    $self->{TESTS}{$test_name} = TestHarness::TestProblem->new({
								NAME          => $test_name,
								PRINT_STRING  => $print_string,
								STATUS        => "QUEUED",
							        RUN_DIR       => $self->testing_dir . "/$test_name",
							       });
  
  }

  return $self->{TESTS}{$test_name};
}


=item add_testsuite: $myobj->add_testsuite(testsuite|testsuite_name);

=cut
sub add_testsuite{
  my $self = shift;
  my $suite_name = shift || croak "TestHarness::Harness::add_testsuite: No testsuite was provided.";

  print "Adding $suite_name\n" if defined $self->verbose;
  unless (exists $self->{TESTSUITES}{$suite_name}) {
    $self->{TESTSUITES}{$suite_name} = TestHarness::TestSuite->new({
								    NAME      => $suite_name,
								   });
  }

  return $self->{TESTSUITES}{$suite_name};
}


sub add{
  my ($self, $test) = @_;

  croak "TestHarness::Harness::add: No test or testsuite given" unless defined $test;

  my $ret;
  if ($test =~ /\.suite/i) {
    # add testsuite 
    $ret = $self->add_testsuite($test);
  }

  else {
    #add test
    $ret = $self->add_test($test, "");
  }
  return $ret;
}

# Generate a list of candidate test directories from the list of test_directories
# This is simply a comprehensive list of all subdirectories under the test_directories
sub candidate_test_directories{
  my $self = shift;

  unless (defined $self->{CANDIDATE_TEST_DIRECTORIES}) {
    $self->{CANDIDATE_TEST_DIRECTORIES} = [];
    my $search_dirs = $self->test_directories;
    find(sub{push @{$self->{CANDIDATE_TEST_DIRECTORIES}}, $File::Find::name if -d}, @$search_dirs);
  }

  return $self->{CANDIDATE_TEST_DIRECTORIES};
}


# find the directory containing test for a give test name
sub get_test_directory{
  my $self = shift;
  my $test_name = shift
    or die "TestHarness::TestContainer::get_test_directory: A name is required to find the test_directory\n";

  my $test_directory = "";
  foreach (@{$self->candidate_test_directories}) {
    if (/\b$test_name\Z/) {
     #my @filelist = glob("$_/*.test $_/run_job.csh");
     #if (not @filelist) { next; }
     #$test_directory = $_;
     #last;
      # Search the candidate directories for either run_job.csh or foo.test. Return the directory of the last one found.
      my @filelist = ("$_/run_job.csh", "$_/$test_name.test");
      foreach my $file (@filelist) {
        if ( -e $file ) {
          $test_directory = $_;
        }
      }
    }
  }

  return $test_directory;
}


sub source{
  my $self = shift;

  my $source = shift || [];
  $source = (ref $source)? $source: [$source];
  push @{$self->{SOURCE}}, @$source;

  # make sure we have the full path foreach file to be sourced
  my $cur_dir = getcwd;
  $cur_dir .= "/";
  foreach my $file (@{$self->{SOURCE}}) {
    $file = $cur_dir . $file 
      unless ($file =~ /^\//);
  }
  return $self->{SOURCE};
}



sub log_file_suffix{
  my $self = shift;

  my $log_file_suffix = shift || $self->{LOG_FILE_SUFFIX};
  $log_file_suffix = (ref $log_file_suffix)? @$log_file_suffix[0]: $log_file_suffix;
  $self->{LOG_FILE_SUFFIX} = $log_file_suffix;
  return $log_file_suffix;
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
  unless (-d $results_path) {
     system("if [ ! -d $results_path ] ; then mkdir -p $results_path; fi") == 0 
       or croak "TestHarness::Harness::results_path: Could not create results_path $results_path: $!\n";
  }
  return $results_path;
}

# sub reports{
#   my $self = shift;

#   my $reports = shift || $self->{REPORTS} || [];
#   if (ref $reports eq "ARRAY") {
#     $self->{REPORTS} = $reports;
#   }
#   else {
#     carp "TestContainer::reports: Bad value in reports command, expected an array ref but received $reports\n" 
#   }
#   return $reports;
# }


sub execs{
  my $self = shift;

  my $execs   = $self->executables;

  my $serial_executables   = $self->serial_executables;
  foreach my $exec (values %$serial_executables) {
    $exec = "SERIALRUN " . $exec unless $exec =~ /SERIALRUN/;
  }

  my $parallel_executables = $self->parallel_executables;
  foreach my $exec (values %$parallel_executables) {
    $exec = "MPIRUN " . $exec unless $exec =~ /MPIRUN/;
  }

  my %executables = (%$execs, %$serial_executables, %$parallel_executables);
  return \%executables;
}


sub data_links{
  my $self = shift;

  if (my $data_links = shift) {
    if (ref $data_links eq "ARRAY") {
      foreach my $link_file (@$data_links){
	my $str = '\$';
#	$link_file =~ s/\$/$str/g;
      }
      push @{$self->{DATA_LINKS}}, @$data_links;
    }
    else {
      carp "TestContainer::data_links: Bad value in data_links command, expected an array ref but received $data_links\n" 
    }
  }
  
  return $self->{DATA_LINKS};
}




sub get_testsuite_directories{
  my $self = shift;
  
  #print Dumper($self);
  #dump_stack("DEBUG: get_testsuite_directories\n");
  #$self->dump;
  my @directories = @{$self->{TESTSUITE_DIRECTORIES}};
  return @directories;
}



sub get_test_directories{
  my $self = shift;
  
  my @directories = @{$self->{TEST_DIRECTORIES}};
  return @directories;
}


sub test_directories{
  my $self = shift;
  
  my $dirs = shift;
  my $local_dir = getcwd;

  if (defined $dirs) {
    foreach my $dir(@$dirs) {
      next unless (-d $dir);
      $dir = $local_dir . "/" . $dir unless ($dir =~ /^\//);
      push @{$self->{TEST_DIRECTORIES}}, $dir;
    }
  }

  
  return $self->{TEST_DIRECTORIES} || ["$local_dir"];
}


sub testsuite_directories{
  my $self = shift;
  
  my $dirs = shift;
  my $local_dir = getcwd;

  if (defined $dirs) {
    foreach my $dir(@$dirs) {
      next unless (-d $dir);
      $dir = $local_dir . "/" . $dir unless ($dir =~ /^\//);
      push @{$self->{TESTSUITE_DIRECTORIES}}, $dir;
    }
  }

  
  return $self->{TESTSUITE_DIRECTORIES} || ["."];
}



sub all_tests{
  my $self = shift;
  return values %{$self->{TESTS}};
}



sub default_num_cpus{
  my $self = shift;

  $self->{DEFAULT_NUM_CPUS} = shift if @_;
  # default_num_cpus should be a scalar. If it is an array then
  # replace it with the first element. 
  $self->{DEFAULT_NUM_CPUS} = $self->{DEFAULT_NUM_CPUS}[0] if (ref $self->{DEFAULT_NUM_CPUS});
  my $val = $self->{DEFAULT_NUM_CPUS};
  return $val;
}

sub default_num_threads{
  my $self = shift;

  $self->{DEFAULT_NUM_THREADS} = shift if @_;
  # default_num_threads should be a scalar. If it is an array then
  # replace it with the first element. 
  $self->{DEFAULT_NUM_THREADS} = $self->{DEFAULT_NUM_THREADS}[0] if (ref $self->{DEFAULT_NUM_THREADS});
  my $val = $self->{DEFAULT_NUM_THREADS};
  return $val;
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


sub default_time_limit{
  my $self = shift;

  $self->{DEFAULT_TIME_LIMIT} = shift if @_;
  my $val = $self->{DEFAULT_TIME_LIMIT};
  return $val;
}

sub default_restart_script{
  my $self = shift;

  my $default_restart_script = shift || $self->{DEFAULT_RESTART_SCRIPT};
  if (ref $default_restart_script) {
    $default_restart_script = $$default_restart_script[0];
  }

  # Get the full path for the default_restart script
  $default_restart_script = $self->cwd . "/" . $default_restart_script unless ($default_restart_script =~ /^\//);
  
  warn "HARNESS::HARNESS:Warning: $default_restart_script is not executable or does not exist.\n" 
    unless (-x $default_restart_script);

  $self->{DEFAULT_RESTART_SCRIPT} = $default_restart_script;
  return $default_restart_script;
}

sub  executables{
  my $self=shift;
  my $hash_ref =  shift || {};
  if (ref $hash_ref eq 'ARRAY') {
    my $arr_ref = $hash_ref;
    $hash_ref = {};
    %$hash_ref = @$arr_ref;
  }

  while (my($key, $value) = each %$hash_ref) {
    $self->{EXECUTABLES}{$key} = $value;
  }
  
  $hash_ref = $self->{EXECUTABLES};
  return $hash_ref;
}

sub  serial_executables{
  my $self=shift;
  my $hash_ref =  shift || {};
  if (ref $hash_ref eq 'ARRAY') {
    my $arr_ref = $hash_ref;
    $hash_ref = {};
    %$hash_ref = @$arr_ref;
  }

  while (my($key, $value) = each %$hash_ref) {
    $self->{SERIAL_EXECUTABLES}{$key} = $value;
  }
  
  $hash_ref = $self->{SERIAL_EXECUTABLES};
  return $hash_ref;
}

sub  parallel_executables{
  my $self=shift;
  my $hash_ref =  shift || {};
  if (ref $hash_ref eq 'ARRAY') {
    my $arr_ref = $hash_ref;
    $hash_ref = {};
    %$hash_ref = @$arr_ref;
  }

  while (my($key, $value) = each %$hash_ref) {
    $self->{PARALLEL_EXECUTABLES}{$key} = $value;
  }
  
  $hash_ref = $self->{PARALLEL_EXECUTABLES};
  return $hash_ref;
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
  


sub prepare_test_directories{
  my $self = shift;
  my $testing_dir = shift || getcwd;

  my $cwd = getcwd;
  my @tests = $self->all_tests;
  foreach my $test(@tests) {
    $test->prepare_test_directory($testing_dir);
    my $new_dir = getcwd;
    if ($cwd ne $new_dir) {
      chdir $cwd or 
	croak "TestHarness::Harness::prepare_test_directories: Could not return to $cwd, $!";
    }
  }
}


sub cts_file{
  my $self = shift;
  my $cts_file = shift || $self->{CTS_FILE} || "Unknown";
  if (ref $cts_file) {
    $cts_file = $$cts_file[0];
  }
  $self->{CTS_FILE} = $cts_file;
  return $cts_file;
}


sub results{
  my $self = shift;
  my $results_hash = $self->{RESULTS};
  return $results_hash;
}

sub report_logo{
  my $self = shift;
  my $logo = shift || $self->{REPORT_LOGO} || $self->project;
  if (ref $logo) {
    $logo = $$logo[0];
  }
  $self->{REPORT_LOGO} = $logo;
  return $logo;
}

sub path{
    my $self = shift;
    my $path = shift || $self->{PATH};
    
    if (ref $path) {
        $path = $$path[0];
    }
    elsif( $path =~ /\S/ ){
        my $cwd = $self->CWD;
        my @paths = split( /:/, $path );
        my $subpath;
        $path = "";
        # fully qualify path
        foreach $subpath ( @paths ){
            if( $subpath !~ /^\// ){
                $subpath = "$cwd/$subpath";
            }
            $path .= "$subpath:";
        }
        $path =~ s/:$//;
    }
    $self->{PATH} = $path;
    return $path;
}
    
sub project{
  my $self = shift;
  my $project = shift || $self->{PROJECT} || "Unknown";

  if (ref $project) {
    $project = $$project[0];
  }
  $self->{PROJECT} = $project;
  return $project;
}


sub platform{
  my $self = shift;
  my $platform = shift || $self->{PLATFORM} || "Unknown";
  $self->{PLATFORM} = $platform;
  return $platform;
}

sub exe_flavor{
  my $self = shift;
  my $exe_flavor = shift || $self->{EXE_FLAVOR} || "Unknown";
  $self->{EXE_FLAVOR} = $exe_flavor;
  return $exe_flavor;
}

sub dump {
  my $self = shift;
  print "HARNESS:", Dumper($self);
}



1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Project - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Project;
  blah blah blah

=head1 ABSTRACT

  Harness
 
=head1 ABSTRACT_FROM

  Harness

=head1 DESCRIPTION

Stub documentation for Project, created by h2xs. It looks like the
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

=cut

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
