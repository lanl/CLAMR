package TestHarness::TestProblem;

use TestHarness::Judge;
use strict;
use warnings;
use Cwd qw(getcwd chdir);
use Carp;
use DumpStack;
use FileHandle;
use File::Find;
use File::Basename;
use FindBin;
use Data::Dumper;
use vars qw($AUTOLOAD);

my @test_files = ();

{   # BEGIN CLOSURE BLOCK: See Perl reference for information on closures.
  my %problem = (
		 BIN_DIR                  => "$FindBin::Bin",

	         DESCRIPTION              => undef,

	         FAST                     => 0,

		 HARNESS                  => undef,

		 LOG_FILE                 => undef    ,

		 MAX_CPUS                 => undef    ,

		 MIN_CPUS                 => undef    ,

		 NAME                     => undef    ,

		 NUM_CPUS                 => undef    ,

		 PRINT_STRING             => ""       ,

		 RESTART_SCRIPT           => "",

		 RESULTS                  => "NOT_RUN",

		 RESULTS_FILE             => "RESULTS",

		 RUNSCRIPT                => undef    ,

		 RUNSCRIPT_FILE           => undef    ,

		 RUN_DIR                  => ""       ,

		 RUN_TIME                 => 0.0    ,

		 STATUS                   => 'unknown',

		 SYSTEM                   => undef,

		 TEST_DIRECTORY           => undef,

		 TEST_FILE                => undef    ,

		 TIME_LIMIT               => "1:00"   ,   # hr:min
		);




  sub AUTOLOAD {
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);
    $method = uc $method;

    return if ($method =~ /destroy/i);
    unless (exists $problem{$method}) {
      dump_stack;
      local $\ = "\n";
      print  "WARNING:: No such attribute: $method is defined for TestHarness::TestProblem.
                            The available attributes are";
      foreach (keys %problem) {print};
      return;
    }

    my $code = q{sub{
      my $self = shift;
      $self->{METHOD} = shift if @_;
      my $val = $self->{METHOD};
      return $val;
    }};

    $code =~ s/METHOD/$method/g;

    no strict "refs";
    *$AUTOLOAD = eval $code;

    goto &$AUTOLOAD;
  }

  
  # TestProblem->new(test_name|file, {NAMED_ARGUMENT=>value, ...};
  sub new {
    my $invocant = shift;

    my $self = {%problem};
    bless $self, ref($invocant) || $invocant;
    $self->_init(@_);
    return $self;
  }

  # _init should only be called internally.
  sub _init{
    my $self = shift;
    my $hashref = shift 
      || croak "TestHarness::TestProblem::new: No test was provided";
    
    $self->harness(TestHarness::Harness->new);

    # initialize test attributes with $hashref
    no strict "refs";
    foreach my $func(keys %$hashref) {
      my $value = $hashref->{$func};
      
      # only access attributes via function calls to catch bad arguments
      $func = lc $func;  # functions are always lowercase
      $value = $self->$func($value);
    }
    
    # find the full path to the test directory
    my $test_dir = $self->test_directory || $self->name || croak "TestHarness::TestProblem::new: No test was provided";
    #my $test_dir = $self->name || croak "TestHarness::TestProblem::new: No test was provided";
 
   $test_dir = $self->harness->get_test_directory($test_dir);
    unless ($test_dir && (-d $test_dir)) {
      $self->status("NO_TEST_FOUND");
      return;
    }

    
    $self->test_directory($test_dir);
    
#     # Source in any tcsh files listed on SOURCE lines;
    $self->{SYSTEM} = System->new({});# unless defined $self->system;
#     foreach my $file (@{$self->harness->source}) {
#       #print $file, "\n";
#       #$self->{SYSTEM}->source($file);
#     }

    # Find the .test file and/or the runscript
    my $test_file;
    my $runscript = "";
    my $runscript_file = "";
    ($test_file) = <$test_dir/*.test>;
    #find (sub {$test_file = $File::Find::name if /\.test\Z/}, $test_dir); 

    $test_file ||= "no_file_found";
    if (-r $test_file) {
      # We want our runscript to use the local version of $test_file which should
      # be copied into the local directory by the test_harness.
      $test_file = basename $test_file;
      $self->test_file($test_file);
      $self->runscript_file("runscript");
    }
    else {
      # find a runscript
      find (sub {$runscript_file = $File::Find::name if /\A run[\w.]* \Z/ix}, $test_dir);  
      # copy the runscript or an error runscript if one is not found.
      if (-r  $runscript_file) {
	# if $runscript_file == run_job.csh then parse it for number of cpus and time/
	my $basename = basename $runscript_file;
	if ($basename =~ /run_job\.csh/) {
	  $self->parse_run_job($runscript_file);
	}

	$runscript_file = $basename;
	$self->runscript_file($runscript_file);
      }
      else {
	my $name = $self->name;
	print "TestHarness::TestProblem:_init:WARNING: Could not find a .test file or a runscript for $name.\n";
	$self->status("NO_TEST_FOUND");
	$self->results("NO_TEST_FOUND");
# # 	$runscript = <<EOF;
# # #!/bin/tcsh -f
# # echo "FAILED: This test($name) has no .test file or runscript associated with it"
# # EOF
# # 	# record the runscript
# # 	$self->runscript($runscript);
# #         $self->runscript_file("runscript");
      }
    }
  }

  sub parse_run_job{
    my $self = shift;
    my $name = $self->name;

    my $runjobfile = shift 
      or croak "TestHarness::TestProblem::parse_run_job: No run_job.csh found for test \"$name\"";
    
    
    open INFILE, "<", "$runjobfile" or croak "TestProblem.pm: Can't open $runjobfile: $!";
    while (my $line =  <INFILE>) {
      next unless $line =~ /^\s*\#/;

      if ($line =~ /^\s*\#rj numpe\s*=\s*(\d+)/i){
	$self->num_cpus($1);
      }
      elsif ($line =~ /^\s*\#rj time\s*=\s*([0-9:]+)/i) {
	$self->time_limit($1);
      }
    }
  }


  sub parse_test_file{
    my $self = shift;

    my $name = $self->name;
    my $test_file = $self->test_file || croak "FAILED: No .test file was found for test \"$name\"";
    my $execs = $self->harness->execs;

  
    open INFILE, "<", "$test_file" or croak "TestProblem.pm: Can't open $test_file: $!";
    my @test_lines = <INFILE>;
    close INFILE;

    my @commands = ();

    my %funcs = ();
    foreach my $line (@test_lines) {
      # remove comments
      $line =~ s/[\#] .* $//x;    
                                 
      # next if line is empty
      next if $line =~ /^ \s* $/x;

      # keyword: lines
      # example "PROCESSORS: 16"
      if ($line =~ /^\s* ([^\s:=]+) [=:]* \s* (.*) $/x) {
	my $func = lc $1;
	my $value = $2;
	if ($func eq "restart") {
	  my $restart_script = ($value)?    $value: 
	  (defined $self->harness->default_restart_script)?  $self->harness->default_restart_script:
	  "";
	  
	  $self->restart_script($restart_script);
	  next;
	}

	if ($self->can($func) || exists $problem{uc $func}) {
	  $self->$func($value) unless exists $funcs{$func};
	  $funcs{$func}++;
	  my $rec = {cmd=>"setenv", args=> uc $func . " $value"};
	  push @commands, $rec;
	  next;
	}
      }

      # get the command and its arguments
      my ($cmd, $args) = $line =~ /^\s* (\S+) \s (.*)$/x;
      # do executable substitutions
      while (my($command, $exec) = each %$execs) {
	# Sometimes the regular expression engine fails to match upper-case
	# words even with the i option, so make sure the search string is lower
	# case. 
	my $lc_cmd = lc $cmd;
	$cmd = $exec  if($command =~ /\b$lc_cmd\b/i);
      }
      #next unless (-x $cmd);

      my $rec = {cmd=>$cmd, args=>$args};
      push @commands, $rec;
    }
  
    return \@commands;
  }


} # END OF CLOSURE BLOCK



sub add_judge{
  my $self = shift;

  my $judge = shift or croak "TestHarness::TestProblem:add_judge: No judge to add.";

  push @{$self->{JUDGES}}, $judge;
}


sub add_anti_judge{
  my $self = shift;

  my $judge = shift or croak "TestHarness::TestProblem:add_judge: No judge to add.";

  push @{$self->{ANTI_JUDGES}}, $judge;
}

sub recreate_runscript{
  my $self = shift;

  my $run_dir = $self->run_dir;
  chdir $run_dir || croak "TestProblem::recreate_runscript: Could not chdir to $run_dir. $!";
  $self->print_runscript;
}

sub create_runscript{
  my $self = shift;

  # need to have a wrapper runscript (runscript) that
  # calls the main runscript (runscript_main)
  # This is because some batch systems have a very limited length for
  # the scripts that gets executed by the batch system
  my $runscript      = "";
  my $runscript_main = "";

  # parse the .test file
  # The test file must be parsed before determining $pe.
  my $commands = $self->parse_test_file();

  my $results_file = $self->results_file;
  my $prob_name = $self->name;
  my $pe = $self->num_cpus;

  my $sys = $self->system;
  my $batch = $sys->batch || "";
  $batch = "" if ($batch =~ /no([-_])?batch/i);
  my $mpirun = $self->harness->mpirun || $sys->mpirun($pe);
  my $mpirun_args = $self->harness->mpirun_args($pe) || $sys->mpirun_args($pe);
  my $mpirun_pre  = $self->harness->mpirun_pre($pe)  || $sys->mpirun_pre($pe);
  my $serialrun = $sys->serialrun();

  # set up the results file.
  #open my $fh, ">", $results_file or die "FAILED: Could not open file $results_file. $!";
  #print $fh "";
  #close $fh;

  #########################################################################################
  # Create a simple tcsh script to run our problem
  $ENV{SHELL} = "/bin/tcsh";
  my $tmpfile = "tmpfile";
  my $bindir = $self->bin_dir;

  my $log = $self->log_file || $self->name . $self->harness->log_file_suffix;
  $log = getcwd . "/$log";
  $self->log_file($log);

  my $output_file_cmd = "";
  if (length($batch)){
    $output_file_cmd = " >>& $log";
  }
  else{
    $output_file_cmd = " |& tee -a $log";
  }

  my $preamble = $sys->batch_preamble({
				       JOBNAME        => $prob_name,
				       PE             => $pe,
				       TIME_LIMIT     => $self->time_limit,
				      });
  #...finish off runscript
  $runscript      .= "#!/bin/tcsh -f\n";
  if (length($preamble)) {
    $runscript .= $preamble;
  }
  $runscript      .= "set cts_start_time = `date +\"%H:%M:%S\"`\n";
  $runscript      .= "rm -f $log\n";
  $runscript      .= "touch $log\n";
  $runscript      .= "cd ".getcwd."\n";
  $runscript      .= "echo ./runscript_main.tcsh\n";
  $runscript      .= "./runscript_main.tcsh\n";
  $runscript      .= "set cts_stop_time = `date +\"%H:%M:%S\"`\n";
  $runscript      .= "echo CTS_START_TIME \$cts_start_time $output_file_cmd\n";
  $runscript      .= "echo CTS_STOP_TIME \$cts_stop_time $output_file_cmd\n";

  #...do runscript_main
  $runscript_main .= "#!/bin/tcsh -f\n";
  $runscript_main .= "cd ".getcwd."\n";
  foreach my $source_file(@{$self->harness->source}) {
    $runscript_main .= "source $source_file\n";
  }
  $runscript_main .= "setenv PATH  $bindir:\$PATH\n";
  # add to the PATH
  my $path = $self->harness->path;
  if( defined($path) && $path =~ /\S/ ){
      $runscript_main .= "setenv PATH $path:\$PATH\n";
  }

  $runscript_main .= "setenv NUM_CPUS $pe\n";
  foreach my $cmd (@$commands) {
    my $prog     = $cmd->{cmd};
    $prog =~ s/MPIRUN(\s+\S+)/$mpirun_pre$mpirun$1 $mpirun_args /;
    $prog =~ s/SERIALRUN/$serialrun/;

    # check for extra lines introduced by serialrun or mpirun into $prog
    my @lines = split /^/, $prog;
    $prog = pop @lines;
    foreach (@lines) {
      $runscript_main .= $_;
    }
    
    my $args     = $cmd->{args};
    
    # append the command line
    my $command_line = "$prog $args";
    if( $command_line !~ /[>\|]/ &&
        $command_line !~ /^\s*(cd|source|set|setenv|if|else|endif|elseif|pushd|popd|foreach|while|end|switch|case|default|breaksw|endsw|\@\s+)\b/ &&
        $command_line !~ /\b(\|\&?\s*tee)\b\n/ )
      {
        $command_line .= " $output_file_cmd";
      }
    $command_line .= "\n";
    $runscript_main .= $command_line;

    if ($command_line !~ /^\s*(if|else|endif|foreach|while|end|switch|case|default|breaksw|endsw|\@\s+)\b/){
      $runscript_main .= "set STATUS=\$?\n";
      $_ = $mpirun;
      if (/mpijob/) {
	$runscript_main .= "if (\"\${STATUS}\" != \"0\" && \"\${STATUS}\" != \"137\") then\n";
      }
      elsif (/poe/) {
	$runscript_main .= "if (\"\${STATUS}\" != \"0\" && \"\${STATUS}\" != \"128\") then\n";
      }
      else{
	$runscript_main .= "if (\"\${STATUS}\" != \"0\") then\n";
      }
      # for echo, replace chars that would interfere with the quote
      $prog =~ s/\"/'/g;
      $args =~ s/\"/'/g;
      # always print to screen and to log file so test returns failure if fail
      # and you see completed in the log file
      $runscript_main .= "  echo \"FAILED -- STATUS is \${STATUS} from command $prog $args\" |& tee -a $log\n";
      $runscript_main .= "  exit(1)\n";
      $runscript_main .= "else\n";
      $runscript_main .= "  echo \"CTS_COMPLETED command $prog $args\" |& tee -a $log\n";
      $runscript_main .= "endif\n"; 
    }
  }
  #########################################################################################


#    my $indent = "  ";
#   my $lib_dirs = "";
#   my @lib_dirs = @{$self->lib_dirs};
#   foreach my $dir (@lib_dirs) {
#     $lib_dirs .= "use lib \"$dir\";\n";
#   }

#   my $parameters = $self->create_parameter_string($indent);
#   my $runscript = "";

#   # read the runscript from the DATA block and modify it with local test information
#   my $position = tell DATA;
#   while (<DATA>) {
#     if (!/^\s*$/ && !/^\s*\#/) { # skip blank and comment lines
#       s/\bTEST_PARAMETERS\b/$parameters/;
#       s/\bLIB_DIRS\b/$lib_dirs/;
#     }
#     $runscript .= $_;
#   }

#   seek DATA, $position, 0;
  return ($runscript, $runscript_main);
}


sub print_runscript{
  my $self = shift;
  my ($runscript, $runscript_main) =
    ($self->test_file)      ? $self->create_runscript:
    ($self->runscript_file) ? ""                     :
      ($self->runscript, $self->runscript_main);
  
   if ($runscript) {
    my $rs = FileHandle->new(">runscript") 
      or croak "TestHarness::TestProblem::prepare_test_directory: Could not open runscript: $!";
  
    print $rs $runscript;
    $rs->close;
    chmod 0770, "runscript";
  }
   if ($runscript_main) {
    my $rs = FileHandle->new(">runscript_main.tcsh") 
      or croak "TestHarness::TestProblem::prepare_test_directory: Could not open runscript_main.tcsh: $!";
  
    print $rs $runscript_main;
    $rs->close;
    chmod 0770, "runscript_main.tcsh";
  }
}



# sub create_parameter_string{
#   my ($self, $indent) = @_;

#   my $string = "{\n";
#   foreach my $key (sort keys %$self) {
#     next if ($key =~ /RUNSCRIPT|SYSTEM|TEST_DIRECTOR|LIB_DIRS/);
#     my $value = $self->$key;
#     next unless $value;
#     my $str_value = (ref $value)? expand_ref($value, "$indent      "): "\"$value\"";
#     $string .=  sprintf "$indent  %-15s => %s,\n", $key, $str_value;
#   }
#   $string .= $indent . "}";

#   # remove quotes from numbers
#   $string =~ s/\"(\d+)\"/$1/mg;

#   # add quotes to certain hash strings to remove a weird bug
#   # {my_command.pl => some_string}  becomes
#   # {"my_command.pl" => some_string}
#   $string =~ s/\b(\w+\.\w+) (\s*) =>/"$1"$2=>/xmsg;

#   return $string;
# }


# sub expand_ref{
#   my ($ref, $indent) = @_;
#   my $string = "";

#   my $type_of_ref = ref $ref;

#   if ($type_of_ref eq 'SCALAR') {$string .= "\"$$ref\""}

#   elsif ($type_of_ref eq 'ARRAY') {
#     $string .= "\n$indent" . "[\n";
#     foreach my $element(@$ref) {
#       $string .= $indent;
#       my $str_value = (ref $element)? expand_ref($element, "$indent "): $element;
#       $string .=  sprintf  "$indent  \"%s\",\n", $str_value;
#     }
#     $string .= $indent . "]";
#   }

#   elsif ($type_of_ref eq 'HASH') {
#     $string .= "\n$indent" . "{\n";
#     foreach my $key(sort keys %$ref) {
#       my $value = $ref->{$key};
#       my $str_value = (ref $value)? expand_ref($value, "$indent "): $value;
#       $string .=  sprintf  "$indent %-15s => \"%s\",\n", $key, $str_value;
#     }
#     $string .= $indent . "}";
#   }

#   return $string;
# }


sub prepare_test_directory{
  my $self = shift;
  my $testing_dir = shift || getcwd;

  my $test_directory = $self->test_directory;

  my $cwd = getcwd;
  chdir $testing_dir or croak "TestHarness::TestProblem::prepare_test_directory: Could not chdir to $testing_dir. $!";

  my $name = $self->name;
  unless (-d $name){
    mkdir $name or croak "TestHarness::TestProblem::prepare_test_directory: Could not create directory $name. $!";
  }

  chdir $name or croak "TestHarness::TestProblem::prepare_test_directory: Could not chdir to $name. $!";

  # let the problem know where it is to run
  my $run_dir = getcwd;
  $self->run_dir($run_dir);

  # Copy or link everything in the test_directory to run_dir. Note: cp returns 0 on success.
  unless ($test_directory) {
    die "Warning: No test was found for $name.\n";
  }
  if ($self->harness->fast) {
    if (system("ln -s $test_directory/* .")){
      carp "TestHarness::TestProblem::prepare_test_directory: Could not link $test_directory/* to $run_dir";
    }
  }
  else {
    if (system("cp -r $test_directory/* .")){
      carp "TestHarness::TestProblem::prepare_test_directory: Could not copy $test_directory to $run_dir";
    }
  }


  $self->add_data_links;

  # Create the runscript or overwrite it if it already exists.
  $self->print_runscript;

  # return to the directory we started from
  my $new_dir = getcwd;
  if ($cwd ne $new_dir) {
    chdir $cwd or 
      croak "TestHarness::Harness::prepare_test_directories: Could not return to $cwd, $!";
  }
}


sub execute{
  my $self = shift;

  chdir $self->run_dir;
  my $runscript=$self->runscript_file;
  my $sys = $self->harness->system;
  my $threads = $self->harness->default_num_threads || "";

  #########################################################################################
  # Run the problem and extract the runtime.
  my $begin_time = time;
  my $pid = 0;
  my $batch = $sys->batch || "";
  if ($batch =~ /no_?batch/i) {$batch = ""}
  if( -e "run_job.csh" ){
      # pick the one in the checkout
      my $run_job = $self->bin_dir."/run_job.pl";
      # otherwise, assume it is in the PATH
      if( ! -e $run_job ){
          $run_job = "run_job.pl";
      }
      # pick the one in the checkout
      my $run_job_cleanup = $self->bin_dir."/run_job_cleanup.pl";
      # otherwise, assume it is in the PATH
      if( ! -e $run_job_cleanup ){
          $run_job_cleanup = "run_job_cleanup.pl";
      }
      my $opts = "";
      # EXEC, SERIAL
      if( ! length($batch) ){
          $opts .= " --var BATCH=no";
      }
      my $executables = $self->harness->execs;
      my $found = "false";
      foreach my $exec_type ( values %$executables ){
          if( $exec_type =~ /^(\S+)\s+(\S+)/ ){
              my $type = $1;
              my $exec = $2;
              $opts .= " --var EXEC=$exec";
              # either MPIRUN or SERIALRUN
              if( $type eq "SERIALRUN" ){
                  $opts .= " --var SERIAL=yes";
              }
              elsif( $type eq "MPIRUN" ){
                  # no need to do anything
              }
              else{
                  die "HARNESS::TESTPROBLEM:Error: Unknown type/exec for executable type=[$type] exec=[$exec]";
              }
              $found = "true";
              last;
          }
      }
      if( $found eq "false" ){
          die "HARNESS::TESTPROBLEM:Error: Could not find executable in self->harness->execs hash";
      }
      # PNAME
      $opts .= " --var PNAME=".$self->name;
      # account
      my $account = $self->harness->account;
      if( defined($account) && $account =~ /\S/ ){
          $opts .= " --var BATCH_ARGS+='-A $account'";
      }
      if( defined($ENV{BATCH_ARGS}) ){
          $opts .= " --var BATCH_ARGS+='$ENV{BATCH_ARGS}'";
      }
      # run cmd
      # remove old log files
      `$run_job_cleanup --clean`;
      my $cmd = "$run_job go --wait $opts";
      # if just setting up files, do not really run things
      if( $self->harness->noop ){
          $cmd = "/bin/true";
      }
      #print "$cmd\n";
      `echo "$cmd\n" > run_cmd`;
      my $name = $self->name;
      print "CTS_INFO: Executing test $name\n";
      my $date_1 = `date +%H:%M:%S`; chomp( $date_1 );
      my $output = `($cmd) 2>&1`;
      my $date_2 = `date +%H:%M:%S`; chomp( $date_2 );
      print $output;
      # create a CTS log file
      my $log_file = $self->log_file || $self->name . $self->harness->log_file_suffix;
      `echo "" > $log_file`;
      # if no files, got a mount problem (moonlight, pinto)
      `if [ ! -e rj_cmd_out ] ; then echo "NOTAVAIL" >> $log_file ; fi`;
      # get times - sum up all the time_elapsed
      `cat rj_*/rj_batch_out.* rj_*/rj_cmd_out.* >> $log_file 2>&1`;
      my $lines = `egrep -a '^RJ_OUTPUT: TIME_ELAPSED_S:' $log_file`;
      my $time_elapsed = 0;
      foreach my $line ( split( /\n/, $lines ) ){
          if( $line =~ /RJ_OUTPUT: TIME_ELAPSED_S:(\S+)/ ){
              $time_elapsed += $1;
          }
      }
      `echo "CTS_START_TIME 0:0:0" >> $log_file`;
      `echo "CTS_STOP_TIME 0:0:$time_elapsed" >> $log_file`;
      # Generate output
      $self->print_results;
  }
  else{
      if (length($batch)){
          my $name = $self->name;
          print "  CTS_INFO: Executing test $name\n";
          $pid = `$batch $runscript 2>&1`;
          #print "$pid\n";
      }
      else{
          #print "DEBUG: RUNNING WITHOUT BATCH\n";
          if( $threads == 1 )
          {
              open( NOBATCH_RUN_FILE, "| ./$runscript");
              print <NOBATCH_RUN_FILE>;
              close( NOBATCH_RUN_FILE );
          }
          else
          {
              `./$runscript`;
          }
      }
      
      my $restart = $self->restart_script || "";
      #print "DEBUG:restart : $restart\n";
      if (-f $restart && -x $restart && !system($restart)) {
          my $name = $self->name;
          print "$name: restarted\n\n";
          $self->execute;
      }
      else {
          my $end_time = time;
          $self->run_time($end_time - $begin_time);
          #   # special code for ASC White
          #   eval $sys->batch_post($pid);
          #   sleep 5;
          #########################################################################################
          
          #########################################################################################
          # Generate output
          $self->print_results;
      }
  }

  return;
}


sub restart_script{
  my $self = shift;

  my $restart_script = shift || $self->{RESTART_SCRIPT} || $self->harness->default_restart_script;
  if (ref $restart_script) {
    $restart_script = $$restart_script[0];
  }

  # Get the full path for the restart script
  my $cwd = getcwd();
  if (-x $restart_script) {
    $restart_script = $cwd . "/" . $restart_script unless ($restart_script =~ /^\//);
  }
  else {
    $restart_script = $self->harness->cwd . "/" . $restart_script unless ($restart_script =~ /^\//);
  }

  die "HARNESS::TESTPROBLEM:Warning: $restart_script is not executable or does not exist.\n" unless (-x $restart_script);

  #print "DEBUG: Restart: $restart_script\n\n";

  $self->{RESTART_SCRIPT} = $restart_script;
  return $restart_script;
}

sub num_cpus{
  my $self = shift;

  my $num_cpus = shift || $self->{NUM_CPUS} || $self->harness->default_num_cpus;
  $self->{NUM_CPUS} = ((defined $self->min_cpus) && ($num_cpus < $self->min_cpus))? $self->min_cpus :
                      ((defined $self->max_cpus) && ($num_cpus > $self->max_cpus))? $self->max_cpus :
                                                                                    $num_cpus;

  return $self->{NUM_CPUS};
}

sub run_time{
  my $self = shift;
  
  #      my $i = 0;
  #       print "\n";
  #       while (my ($calling_package, $filename, $line) = caller($i++) ){
  # 	print "Called by $calling_package, $filename, $line\n";
  #       }
  #       print "\n";

  if (@_)
    {
      $self->{RUN_TIME} = shift;
      #...adjust time to be within 0 through number of seconds in a day
      #...Fix needed for when a test starts before midnight and ends after
      $self->{RUN_TIME} = ($self->{RUN_TIME}) % (60*60*24);
    }
  return $self->{RUN_TIME};
}


sub add_data_links {
  my $self = shift;
  my $data_links = $self->harness->data_links;
  my %execs = (%{$self->harness->executables},  %{$self->harness->serial_executables},  
	       %{$self->harness->parallel_executables});
  foreach my $file (@$data_links) {
    while (my ($exec,$value) = each ( %execs ) ) {
       my $path = dirname($value);
       my $lcexec = lc $exec;
       $file =~ s/\$ $lcexec [._] dir/$path/ix;
    }  
    my $newfile = basename $file;
    symlink $file, $newfile or warn "TestProblem::add_data_links: Could not add link $file: $! \n";
  }
}

sub extract_results{
  my $self = shift;
  my %results = ();

  $results{NAME} = $self->name;
  $results{RUN_DIR} = $self->run_dir;

  # Get the results file produced by the test when it ran.
  my $results_file = $self->run_dir. "/" . $self->results_file;
  my $fh = new FileHandle;
  if ($fh->open("< $results_file")) {
    local $/;
    $results{RESULTS_FILE} = <$fh>;
    $fh->close;

    my ($first_line, $run_time_line) = split "\n", $results{RESULTS_FILE};
    if ($first_line =~ /\b(passed|failed|unknown|diff|nan|notavail|sysfail)\b/i) {
      $self->status(uc $1);
    }
    if ($run_time_line =~ /RUN_TIME: \s* (\d+)/x) {
      $self->run_time($1);
    }
  }
  else {
    $results{RESULTS_FILE} = "No results file was found.";
    $self->status("FAILED");
  }

  $results{STATUS} = $self->status;

  my $time_limit = $self->time_limit;
  $time_limit = _normalize_time($time_limit); # Get the time limit in seconds.

  my $run_time = $self->run_time;

  $results{RUN_TIME} = $run_time;
  $results{TIME_LIMIT} = $time_limit;

  $self->results(\%results);
  return \%results;
}


sub _normalize_time{
  my $time = shift;

  my $hour = 0;
  my $min  = 0;
  my $sec  = 0;

  if (($hour, $min, $sec) = $time =~ /^(\d*):(\d*)  (?: :(\d*))?/x) {
    $sec ||= 0;
    $sec += 3600*$hour + 60*$min;
  }
  else {
    $sec = 60*$time; # assume that the time is given in minutes.
  }
  
  return $sec;
}


sub print_results{
  my $self = shift;
  my %results = ();
  $results{NAME} = $self->name;
  $results{RUN_DIR} = $self->run_dir;

  my $log_file = $self->log_file || $self->name . $self->harness->log_file_suffix;
  my $results_file = $self->run_dir . "/" . $self->results_file;
  # try to flush the log file if it is not ready yet (batch system fluke)
  my $flusher = `ls -la $log_file 2>&1 ; tail $log_file 2>&1`;
  my $retry = 0;
  while( $flusher !~ /CTS_STOP_TIME/ && $retry < 60 )
    {
      $flusher = `ls -la $log_file 2>&1 ; tail $log_file 2>&1`;
      sleep( 1 );
      $retry++;
    }
  if (! -f $log_file) {
    `echo "FAILED  Problem did not run." > $log_file`
  }
  # try to flush the log file
  open my $outputfh, "<", $log_file or die "TestProblem:print_results: Could not open file $log_file. $!";

# print "RESULTS FILE IS : $results_file\n";
  open my $resultsfh, ">", $results_file or die "FAILED: TestProblem:print_results: Could not open file $results_file. $!";

  my ($passed_count, $diff_count, $failed_count) = (0, 0, 0); 
  my %prob_result = (
		     NAN      => 0,
		     CORE     => 0,
		     NOTAVAIL => 0,
		     SYSFAIL  => 0,
		     FAILED   => 0,
		     DIFF     => 0,
		     PASSED   => 0,
		    );

  my $cts_start_time = "0";
  my $cts_stop_time  = "0";
  while (my $line = <$outputfh>) {
    #print $line;
    if (my ($tag, $time) = $line =~ /CTS_([A-Z]+)_TIME \s* (.*) $/x) {
      if ($tag =~ /START/) {
	$cts_start_time = $time;
      }
      elsif ($tag =~ /STOP/) {
	$cts_stop_time = $time;
      }
    }
    #...skip over lines that just echo the command being done...
    elsif( $line =~ /^CTS_COMPLETED command/ )
      {
        next;
      }
    elsif ($line =~ /\b(FAILED|DIFF|PASSED|CORE|NAN|NOTAVAIL|SYSFAIL)\b/) {$prob_result{$1}++}
  }
  close $outputfh;
  
  $cts_start_time = _normalize_time($cts_start_time);
  $cts_stop_time  = _normalize_time($cts_stop_time);

  # calculate the run_time.
  my $run_time = $self->run_time($cts_stop_time - $cts_start_time);
  
  
  my $passed   = $prob_result{"PASSED"};
  my $failed   = $prob_result{"FAILED"};
  my $core     = $prob_result{"CORE"};
  my $nan      = $prob_result{"NAN"};
  my $sysfail  = $prob_result{"SYSFAIL"};
  my $notavail = $prob_result{"NOTAVAIL"};
  my $diff     = $prob_result{"DIFF"};
  
  
  my $status = ($core)     ? "CORE"     :
               ($nan)      ? "NAN"      :
               ($notavail) ? "NOTAVAIL" :
               ($sysfail)  ? "SYSFAIL"  :
               ($failed)   ? "FAILED"   :
               ($diff)     ? "DIFF"     :
	       ($passed)   ? "PASSED"   :
	                     "UNKNOWN"  ; 
  
  $failed+=$core;
  $failed+=$nan;
  $failed+=$sysfail;

  my $full_result = "P-$passed,   D-$diff,   F-$failed"; 
  
  my $summary = sprintf "%-10s %16s\n", $status, $full_result;
  print $status, " -- ", $self->name, "   $run_time secs\n";
#  print $summary;
  print $resultsfh $summary;
  print $resultsfh "RUN_TIME: ", $run_time, "\n";
  close $resultsfh;
  $results{STATUS} = $status;
  $results{RUN_TIME} = $run_time;
  $results{TIME_LIMIT} = _normalize_time($self->time_limit); # Get the time limit in seconds.
  $self->results(\%results);
}

sub dump {
  my $self = shift;
  my $name = $self->name || "NoName";
  print "TestProblem:$name: ", Dumper($self);
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

