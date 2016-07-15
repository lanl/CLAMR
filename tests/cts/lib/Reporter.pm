package Reporter;

use 5.006;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/lib";
use Carp;
use Cwd;
use Data::Dumper;
use DumpStack;
use FileHandle;
use File::Basename;
use Mail::Sendmail;
use vars qw($AUTOLOAD $VERSION %mod);
use subs qw(comment reports local_report email_recipients scp_results 
	    cross_platform_reporter cross_platform_base cross_platform_pattern 
	    cross_platform_mode cross_platform_group cross_platform_mail cross_platform_sort
            header footer);

$VERSION    =   0.01;

my $debug = 0;
my ($CROSS_PLATFORM_REPORTER, $CROSS_PLATFORM_BASE, $CROSS_PLATFORM_PATTERN,
    $CROSS_PLATFORM_MODE, $CROSS_PLATFORM_GROUP, $CROSS_PLATFORM_MAIL, $CROSS_PLATFORM_SORT);
my %itermax;  # This hash is of the form (build=>number_of_builds, example=>number_of_examples)
my (%field1, %field2, %field3, %field4, %field5, %dir);
my $reporter = undef;


{   # BEGIN CLOSURE BLOCK: See Perl reference for information on closures.
  my %report = (
		COMMENT                 => "",
		
		HARNESS                 => undef,

		HEADER                  => "",
		FOOTER                  => "",

		REPORTS                 => undef,

		LOCAL_REPORT            => undef,

		EMAIL_RECIPIENTS        => undef,
		
		RESULTS_PATH            => undef,
		
		SCP_RESULTS             => "no",
		SYSTEM_NAME             => lc `uname -n`,
		CROSS_PLATFORM_REPORTER => undef,
		CROSS_PLATFORM_BASE     => undef,
		CROSS_PLATFORM_PATTERN  => undef,
                CROSS_PLATFORM_MODE     => undef,
                CROSS_PLATFORM_GROUP    => undef,
                CROSS_PLATFORM_MAIL     => undef,
                CROSS_PLATFORM_SORT     => undef,

		PROJECT                 => "UNKNOWN",
		PLATFORM                => "UNKNOWN",
		PROJECT                 => "UNKNOWN",
		REPORT_LOGO             => "$FindBin::Bin/../html/cts_logo2.jpg",
       
		TEXT_FILE               => "cts_results_pid.txt",
		TEXT_OUTPUT             => "",
		HTML_FILE               => "cts_results_pid.html",
	       );




  sub AUTOLOAD {
    my $self = shift;
    my ($package, $method) = ($AUTOLOAD =~ /(.*)::(.*)/);
    $method = uc $method;

    return if ($method =~ /destroy/i);
    unless (exists $report{$method}) {
      dump_stack;
      local $\ = "\n";
      print "WARNING:: $method is not yet implemented in Reporter.pm.\n";
      print "The available methods are";
      foreach (keys %report) {print}
      return;
    }

    $self->{$method} = shift if @_;
    my $val = $self->{$method};
    return $val;
  }

  sub new {
    my $invocant = shift;

    $reporter = {%report};
    bless $reporter, ref($invocant) || $invocant;

    $reporter->harness(TestHarness::Harness->new);

    return $reporter;
  }


} # END OF CLOSURE BRACKETS



sub text {     # writes a brief summary of test results
  my $self = shift;
  my $name = shift ||  $self->harness->text_file || $self->text_file;
  $name = ($name =~ /\.txt$/i)? $name: $name . ".txt";
  $self->text_file($name);
  my $test;
  my $key;
  my $test_result;
  my $run_time;
  my $sum_run_time = 0.0;
  my $text_string = "";
  my $db_file = $self->harness->db_file;

  $self->{TEXT_OUTPUT} = "";
  print "\n     Test Summary\n\n";
# foreach my $key (sort keys(%{$self})){
#   print "$key $self->{$key}\n";
# }
# foreach my $key (sort keys(%{$self->{SYSTEM}})){
#   print "$key $self->{SYSTEM}{$key}\n";
# }
  foreach $test (sort(keys %{$self->harness->results})) {
    $test_result = $self->harness->results->{$test}{STATUS};
    $run_time    = $self->harness->results->{$test}{RUN_TIME};
    $sum_run_time = $sum_run_time + $run_time;
    $self->{TEXT_OUTPUT} .= sprintf "%-25s %-10s %10s\n", $test, $test_result, $run_time;
    #printf "%-16s %-15s %-15s\n", $test, $test_result, $run_time;
    # This is to print out the keys in the hash 
    if ($debug){
      foreach $key (keys(%{$self->harness->results->{$test}})){
        print "$key $self->harness->results->{$test}{$key}\n";
      }
    }
  }
  $self->{TEXT_OUTPUT} .= sprintf "\nTotal Run Time is %s seconds (%.2f minutes)\n",
     $sum_run_time, $sum_run_time/60.0;
  # Store some other info (for printing later)
  my $rundir = getcwd();
  $self->{TEXT_OUTPUT} .= sprintf "\nPATH = $rundir\n";
  my $uname = `uname -n`;
  $self->{TEXT_OUTPUT} .= sprintf "\nUNAME = $uname\n";

  # Print the results to file and stdout. The file may include the
  # the process pid to make it unique.
  my $text_file = $self->harness->results_path . "/" . $self->text_file;
  my $pid = $$;
  $text_file =~ s/([-_.])pid([-_.])/$1$pid$2/i;

  my $text = $self->text_output;
  my $fh = FileHandle->new("> $text_file") or warn "Reporter::text: Could not create text_file\n";
  if (defined $fh) {
    print $fh $text;
    close $fh;
  }
  print $text;

  $self->_init_cross_platform();

  if (defined($CROSS_PLATFORM_BASE) && defined($CROSS_PLATFORM_PATTERN)) {
#   my $system = lc $self->PLATFORM;
#    my $system = lc `uname -n`;
    my $system = $self->harness->system_name;
    my $mpirun = $self->harness->system->mpirun;

# This is to allow a reference system to be set via an environmental variable
# in a cron job
    if (defined($ENV{REFERENCE})) {
      $system = "Reference";
    }

    ($mpirun) = split(' ', $mpirun);
    my $location=`which ${mpirun} 2> /dev/null`;
    my $filetype=`file $location 2>&1`;
    my $mpi="undef";
    if ($mpirun =~ /prun/) {
      $mpi="native";
    }
    elsif ($filetype =~ /shell/){
      $mpi="mpich";
    }
    elsif ($filetype =~ /mpijob/){
      $mpi="mpich";
    }
    elsif ($filetype =~ /executable/){
      $mpi="lampi";
    }
    else{
      $mpi="unknown";
    }

    my ($exec) = (keys(%{$self->harness->parallel_executables}));
    if (! defined($exec)) {
      ($exec) = (keys(%{$self->harness->serial_executables}));
      $location=$self->harness->serial_executables->{$exec};
    }
    else {
      $location=$self->harness->parallel_executables->{$exec};
    }
    my $type;
    if ($location =~ /debug/){
      $type = "debug";
    }
    elsif ($location =~ /optimize/){
      $type = "optimize";
    }
    else {
      $type = "unknown";
    }
#   print "Type is $type\n";

    my $output;
    my $Fcompiler;
    my $systype=`uname`;
    if ($systype =~ "OSF1"){
      $Fcompiler = "native";
    }
    elsif ($systype =~ "Linux"){
      my $FC=$ENV{FC};
      unless ($FC){
        $FC = "unknown";
        foreach my $comp qw(lf95 pgf90 ifort ifc f90){
          $output = `which $comp 2>&1`;
          if ($output =~ /Command not found/){
            next;
          }
          chomp $output;
          if (-x $output){
            $FC = $comp;
            last;
          } 
        }
      }
      if ($FC =~ /lf95/) {
        $Fcompiler = "lahey";
      }
      elsif ($FC =~ /pgf/) {
        $Fcompiler = "pgi";
      }
      elsif ($FC =~ /ifort/) {
        $Fcompiler = "intel";
      }
      elsif ($FC =~ /ifc/) {
        $Fcompiler = "intel";
      }
      else {
        $Fcompiler = "unknown";
      }
    }
    else{
      $Fcompiler = "unknown";
    }
#   print "Fortran compiler is $Fcompiler\n";

    my $Ccompiler;
    if ($systype =~ "OSF1"){
      $Ccompiler = "native";
    }
    elsif ($systype =~ "Linux"){
      my $CC = $ENV{CC};
      unless ($CC){
        $CC = "unknown";
        foreach my $comp qw(pgcc gcc){
          $output = `which $comp 2>&1`;
          if ($output =~ /Command not found/){
            next;
          }
          chomp $output;
          if (-x $output){
            $CC = $comp;
            last;
          } 
        }
      }
      if ($CC =~ /pgcc/) {
        $Ccompiler = "pgi";
      }
      elsif ($CC =~ /gcc/) {
        $Ccompiler = "gcc";
      }
      else {
        $Ccompiler = "unknown";
      }
    }
    else{
      $Ccompiler = "unknown";
    }
#   print "C compiler is $Ccompiler\n";

    my $xplatform_file=$CROSS_PLATFORM_BASE."/".$CROSS_PLATFORM_PATTERN;
    $xplatform_file =~ s/\$\{system\}/$system/;
    $xplatform_file =~ s/\$\{Fcompiler\}/$Fcompiler/;
    $xplatform_file =~ s/\$\{Ccompiler\}/$Ccompiler/;
    $xplatform_file =~ s/\$\{mpi\}/$mpi/;
    $xplatform_file =~ s/\$\{type\}/$type/;
    my $dirpart=dirname($xplatform_file);
    $xplatform_file = $dirpart . "/" . $self->text_file;

    if (-f $xplatform_file){
      unlink $xplatform_file;
    }

    `if [ ! -d $dirpart ] ; then mkdir -p $dirpart; fi`;
    while (length($dirpart) > length($CROSS_PLATFORM_BASE)){
#     print "DEBUG dirpart is $dirpart\n";
      if (defined($CROSS_PLATFORM_GROUP)){
        `chgrp $CROSS_PLATFORM_GROUP $dirpart`;
      }
      if (defined($CROSS_PLATFORM_MODE)){
        `chmod $CROSS_PLATFORM_MODE $dirpart`;
        `chmod +x $dirpart`;
      }
      chop $dirpart;
      $dirpart=dirname($dirpart);
    }
 
    $fh = FileHandle->new("> $xplatform_file") or warn "Reporter::text: Could not create text_file\n";
    if (defined $fh) {
      print $fh $text;
      close $fh;
      if (defined($CROSS_PLATFORM_GROUP)){
        `chgrp $CROSS_PLATFORM_GROUP $xplatform_file`;
      }
      if (defined($CROSS_PLATFORM_MODE)){
        `chmod $CROSS_PLATFORM_MODE $xplatform_file`;
      }
    }
    $db_file=basename($db_file);
#   print "DB_file is ",$db_file,"\n";
    $xplatform_file = dirname($xplatform_file) . "/" . $db_file;
    if (-f $xplatform_file) {
      unlink $xplatform_file;
    }
    `cp $db_file $xplatform_file`;
    if (defined($CROSS_PLATFORM_GROUP)){
      `chgrp $CROSS_PLATFORM_GROUP $xplatform_file`;
    }
    if (defined($CROSS_PLATFORM_MODE)){
      `chmod $CROSS_PLATFORM_MODE $xplatform_file`;
    }
    
  }
  
}


sub cts_file{
  my $self = shift;
  my $cts_file = shift || $self->harness->cts_file || "./";
  if (ref $cts_file) {
    $cts_file = $$cts_file[0];
  }
  return $cts_file;
}


sub results_path{
  my $self = shift;
  my $results_path = shift || $self->harness->results_path;
  $results_path = $self->harness->results_path($results_path);
  return $results_path;
}


sub dump{
  my $self = shift;
  print "\nReport: ", Dumper($self), "\n";
}


# ######  Email brief summary to recipients  ################
# sub email_brief1_txt {
#   my ($ds) =`/usr/bin/date +%b%d_%y`;    # date stamp
#   my ($mailtool) = " mailx ";            # will be passed via system module
#   my ($recipient) = @_;                  # will be specified in CTS file
#   system "cat Results.txt | $mailtool -s 'CTS - $ds' $recipient";  # email report to recipient(s)
# }


# ######  Write HTML report from results file  ##############

sub html{
  my $self = shift;
  my $name = shift ||  $self->harness->html_file ||  $self->html_file;
  $name = ($name =~ /\.html$/i)? $name: $name . ".html";
  $self->harness->html_file($name);

  my $project = $self->harness->project;
  my $logo = $self->report_logo;
  my $comment = $self->comment;
  my $machine = $self->harness->platform;

  my $suites = $self->harness->testsuites;
  if (defined $suites) {
    foreach my $suite(values %$suites) {
      $suite->Reporter::print_html($project, $logo, $comment, $machine);
    }
  }
  $self->harness->Reporter::print_html($project, $logo, $comment, $machine);
}

sub print_html{
  my ($self, $project, $logo, $comment, $machine) = @_;

#   my $suites = $self->testsuites;
#   if (defined $suites) {
#     foreach my $suite(values %$suites) {
#       $suite->Reporter::print_html($project, $logo, $comment, $machine);
#     }
#   }
  
  # get results from local RESULTS files
  my $results_hash = $self->results;
  my $num_tests = keys %$results_hash;
  my $failed_color = "#FF6666";   # light red color
  my $passed_color = "#99FF99";   # light green color
  my $full_date_stamp = localtime;
  my $num_passed = 0;
  my $name = $self->name;
  my $file = ($self->can("suite_file"))? $self->suite_file :
  	     ($self->can("cts_file"))?   $self->cts_file   :
		                         "./"              ;
  my $label = $file;
  $label =~ s/.*\///;
  my $suite_time = 0;
  foreach my $prob (keys %$results_hash){
    $$results_hash{$prob}{STATUS_COLOR} = 
      ($$results_hash{$prob}{RESULTS_FILE} =~ /passed/i)? $passed_color: $failed_color;

    if ($$results_hash{$prob}{RESULTS_FILE} =~ /passed/i) {$num_passed++};

    my $time_limit = $$results_hash{$prob}{TIME_LIMIT};
    my $run_time   = $$results_hash{$prob}{RUN_TIME};
    $suite_time = $suite_time + $run_time;

    $$results_hash{$prob}{TIME_COLOR} = ($time_limit >= $run_time)? $passed_color: $failed_color;
  }  
  my $suite_time_min = int(($suite_time) / 60);
  my $suite_time_sec = ($suite_time - $suite_time_min * 60);
  if ($suite_time_sec <= 9) {$suite_time_sec = "0".$suite_time_sec};

  # Write the html file #################################
  #dump_stack();
  my $html_file = $self->results_path . "/" . $self->html_file;
  my $pid = $$;
  #print "DEBUG: HTML_FILE $html_file\n";
  $html_file =~ s/([-_.])pid([-_.])/$1$pid$2/i;
  #print "DEBUG: HTML_FILE $html_file\n\n";
  open (HTMLFILE, ">$html_file") || die "Reporter::html: Could not open the $html_file for writing: $!\n";

  my $header = (defined $reporter)? $reporter->header: "";
  my $footer = (defined $reporter)? $reporter->footer: "";
  print HTMLFILE qq~ 
    <html><head><title>$project Suite Report</title>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"></head>
    <body bgcolor="#FFFFCC" text="#000000">
    $header
    <table width="562" border="0" cellpadding="0" cellspacing="0" mm:layoutgroup="true">
      <tr><td valign="top" height="57" nowrap colspan="3"> <p><img src = "$FindBin::Bin/../html/cts_banner.jpg"></p></td></tr>
    </table>
    <table width="562" border="0" cellpadding="0" cellspacing="0" mm:layoutgroup="true">
      <tr>
        <td width="220" valign="top" height="73" nowrap>
        <p><b><font face="Verdana, Arial, Helvetica, sans-serif"><b><font color="#000066">
          $project Project<br>
          $machine Machine</font></b></font></b><b><font face="Verdana, Arial, Helvetica, sans-serif"><b><br>
          <font color="#000066">$full_date_stamp<br></font></b></p>
        </td>
        <td width="200" valign="top" height="73" nowrap> <div align="left">
          <h2><font face="Arial, Helvetica, sans-serif"><font color="#000066">$label</font>
          <font color="#990000"><br>$num_passed/$num_tests</font><font color="#000066"> Passed<br>
          <font color="#006600">$suite_time_min:$suite_time_sec</font> min:sec</font></font></h2></div>
        </td>
        <td width="129" valign="top" height="73" nowrap>
          <h2><img src="$logo" width="129" height="73"></h2>
        </td>
      </tr>
      <tr> <td height="12" width="4"></td> <td width="1"></td> </tr>
    </table>
    <table width="562" border="0" cellpadding="0" cellspacing="0" mm:layoutgroup="true">
      <td width="1"></td>
        <table border="1" width="559">
          <tr>
            <td width="24"><b>#</b></td>
            <td width="180"><b>Case, Test Result</b></td>
            <td width="76"><b>Exe Time</b></td>
            <td width="251"><b>Links </b><a href="$file">$label</a></td>
          </tr>~;  # end of print-to-file segment

    # add the table lines
    my $tests = $self->tests;
    if (defined $tests) {
      my $j = 0;
      foreach my $test(values %$tests) {
        my $name=$test->name;
        my $exe_time=$test->run_time;
        my $exe_time_min = int(($exe_time) / 60);
        my $exe_time_sec = ($exe_time - $exe_time_min * 60);
        if ($exe_time_sec <= 9) {$exe_time_sec = "0".$exe_time_sec};
        my $dot_test=$test->test_file || $name . ".test";
        my $outfile = $name . ".out";
        my $bsubfile = "bsub_" . $name;
        $j++;
        print HTMLFILE qq~
          <tr>
            <td width="24">$j</td>
            <td width="180" bgcolor = "$$results_hash{$name}{STATUS_COLOR}">$$results_hash{$name}{NAME}</td>
            <td width="76" bgcolor = "$$results_hash{$name}{TIME_COLOR}">$exe_time_min:$exe_time_sec</td>
            <td width="251"><a href="$$results_hash{$name}{NAME}/$dot_test">.test</a>, 
                            <a href="$$results_hash{$name}{NAME}/$outfile"> out_file</a>, 
                            <a href="$$results_hash{$name}{NAME}/$bsubfile"> bsub_out</a></td>
          </tr>~; # end of print-to-file segment
        } #end for
      } #end if

  # wrap up the html file
  print HTMLFILE qq~
    </table></td></tr></table>
    $footer
    </body>
    </html>~;  # end of print-to-file
  close(HTMLFILE) or die "Reporter::html: Can't close HTMLFILE: $!\n";
  # html file is done #################################

  # chmod for products
  #system "/usr/bin/chmod -R 770 $html_file";
  # Replace with perl chmod
  chmod(0770,$html_file);
  
} # end of subroutine write_brief1_html

sub cross_platform_text{
  my $self = shift;

  $self->_init_cross_platform_report || return;

  my $text_file = $self->harness->results_path . "/TestSummary.txt";
  my $fh = FileHandle->new("> $text_file") or warn "Reporter::text: Could not create text_file\n";

  my $text_output = "";

  my @kv_results; # key/value pairs found in results file

  foreach my $key(sort keys %itermax) {
    undef( @kv_results );
    # Print out the header with some adjustment for the number of columns.
    $text_output .= sprintf("\n\n         ");
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      $text_output .= sprintf("    ");
    }

    $text_output .= sprintf("CROSS_PLATFORM RESULTS ($key)\n\n");
    
    if (exists($field1{$key}->[0])){
      $text_output .= sprintf("               ");
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$text_output .= sprintf("%-9s ",$field1{$key}->[$iter]);
      }
      $text_output .= sprintf("\n");
    }
    
    if (exists($field2{$key}->[0])){
      $text_output .= sprintf("               ");
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$text_output .= sprintf("%-9s ",$field2{$key}->[$iter]);
      }
      $text_output .= sprintf("\n");
    }
    
    if (exists($field3{$key}->[0])){
      $text_output .= sprintf("               ");
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$text_output .= sprintf("%-9s ",$field3{$key}->[$iter]);
      }
      $text_output .= sprintf("\n");
    }
    
    if (exists($field4{$key}->[0])){
      $text_output .= sprintf("               ");
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$text_output .= sprintf("%-9s ",$field4{$key}->[$iter]);
      }
      $text_output .= sprintf("\n");
    }
    
    if (exists($field5{$key}->[0])){
      $text_output .= sprintf("               ");
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$text_output .= sprintf("%-9s ",$field5{$key}->[$iter]);
      }
      $text_output .= sprintf("\n");
    }
    
    my ($file, $prob_name, $result);
    my %prob_data;
    my (@test_count, @test_pass);
    my $longest_key = "";
    for (my $iter=0; $iter<$itermax{$key}; $iter++){

      $file=$dir{$key}->[$iter];
      if (! open(RESULTSFILE, $file) ){
	print "Cannot open file[$file].\n";
	next;
      }
      #print "DEBUG: opened $file\n";
      $test_count[$iter]=0;
      $test_pass[$iter]=0;
      while ( defined( my $line = <RESULTSFILE> )){
	chomp $line;
	#     print "line is $line\n";
	$_ = $line;
	if (! /Test Result/){
	  if (/\b(Anomaly|Pass|Fail|PASSED|DIFF|FAILED|UNKNOWN|NOTAVAIL|SYSFAIL)\b/ ){
	    ($prob_name, $result) = split(' ',$line);
	    $prob_data{$prob_name}[$iter]=$result;
	    
	    if ($result ne "NOTAVAIL" ) {
	      $test_count[$iter]++;
	    }
	    
	    if ($result eq "PASSED" || $result eq "Pass" ) {
	      $test_pass[$iter]++;
	    }
	  }
          # read in misc key/value pairs
          if ( /(\S+)\s*=\s*(\S+.*?)\s*/ ){
            $kv_results[$iter]{$1} = $2;
            if( length($1) > length($longest_key) ){
              $longest_key = substr( "                          ",
                                     0, length($1) );
            }
          }
	}
      }
      
    }
    $text_output .= sprintf("\n");
    $text_output .= sprintf("Totals         ");
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      $text_output .= sprintf("%-8s  ",$test_pass[$iter]."/".$test_count[$iter]);
    }
    $text_output .= sprintf("\n\n");

    foreach my $prob (sort keys(%prob_data)){
      $text_output .= sprintf("%-14s ",$prob);
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	if (! defined($prob_data{$prob}[$iter])){
	  $prob_data{$prob}[$iter]="    ";
	}
	$text_output .= sprintf("%-8s  ",$prob_data{$prob}[$iter]);
      } # foreach testcase
      $text_output .= sprintf("\n");
    } # foreach prob

    $text_output .= sprintf("\n\n         ");
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      $text_output .= sprintf("    ");
    }
    $text_output .= sprintf("Key=Value ($key)\n\n");
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      my $fields_name = "";
      if (exists($field1{$key}->[0])){
        $fields_name .= ", $field1{$key}->[$iter]";
      }
      if (exists($field2{$key}->[0])){
        $fields_name .= ", $field2{$key}->[$iter]";
      }
      if (exists($field3{$key}->[0])){
        $fields_name .= ", $field3{$key}->[$iter]";
      }
      if (exists($field4{$key}->[0])){
        $fields_name .= ", $field4{$key}->[$iter]";
      }
      if (exists($field5{$key}->[0])){
        $fields_name .= ", $field5{$key}->[$iter]";
      }
      $fields_name =~ s/^, //;
      $text_output .= sprintf("   %s%s\n", $longest_key, $fields_name);
      foreach my $key_kv_results ( sort keys %{$kv_results[$iter]} )
        {
          my $spaces = substr( "                    ",
                               0, length($longest_key) - length($key_kv_results) );
          $text_output .= sprintf("%s%s = %s\n",
                                  $spaces, $key_kv_results,
                                  $kv_results[$iter]{$key_kv_results} );
        }
    }
  }
    
  # print the results to file and stdout
  if (defined $fh) {
    print $fh $text_output;
    close $fh;
  }
  
  print $text_output;
  
  $text_file = $CROSS_PLATFORM_BASE . "/TestSummary.txt";
  if (-f $text_file) {
    unlink($text_file);
  }
  $fh = FileHandle->new("> $text_file") or warn "Reporter::text: Could not create text_file\n";
  if (defined $fh) {
    print $fh $text_output;
    close $fh;
    if (defined($CROSS_PLATFORM_GROUP)){
      `chgrp $CROSS_PLATFORM_GROUP $text_file`;
    }
    if (defined($CROSS_PLATFORM_MODE)){
      `chmod $CROSS_PLATFORM_MODE $text_file`;
    }
  }
  
  if (defined($CROSS_PLATFORM_MAIL)) {
    my $smpt = $CROSS_PLATFORM_MAIL;
    $smpt =~ s/[-\w]+@//;
    my $from = $ENV{LOGNAME} . '@' . $smpt;
    $smpt = "mail." . $smpt;
    unshift @{$Mail::Sendmail::mailcfg{'smtp'}} , $smpt;
    my $cronjob = "";
    if( defined( $ENV{CRONJOB}) )
      {
        $cronjob = ": $ENV{CRONJOB}";
      }
    my %mail = (
		To      => $CROSS_PLATFORM_MAIL,
		From    => $from,
                Subject => "CTS Tests: ".$self->harness->project.": ".
                           $self->harness->system_name.
                           $cronjob,
		Message => $text_output,
	       );
    
    sendmail(%mail) or die "Reporter:: Sendmail error: $!";
  }
}
  
sub cross_platform_html{
  my $self = shift;

  my $project = $self->harness->project;
  my $comment = $self->comment;
  my $system = $self->harness->platform;
  my $logo = $self->report_logo;

  my $passed_color   = "#99FF99";   # light green color
  my $diff_color     = "#FFFF00";   # yellow color
  my $notavail_color = "#FF66FF";   # purple color
  my $missing_color  = "#66FFFF";   # aqua color
  my $sysfail_color  = "#FF9900";   # orange color
  my $failed_color   = "#FF6666";   # light red color

  $self->_init_cross_platform_report || return;

  my $html_output;

  my @kv_results; # key/value pairs found in results file
  
# Need to add project name to header in next line
  my $full_date_stamp=`date`;

  my $boundary = "boundary-line";
# start of print-to-string segment
  open BANNER, "$FindBin::Bin/../html/cts_banner.eml" || 
    die "Could not open $FindBin::Bin/../html/cts_banner.jpg: $!\n";
  $logo =~ s/\.jpg/.eml/;
  open LOGO, $logo ||
    die "Could not open $logo: $!\n";
  my ($banner, $logob);
  {
    local $/;
    $banner = <BANNER>;
    $logob = <LOGO>;
  }
  close BANNER;
  close LOGO;

  my $header = $self->header;
  my $footer = $self->footer;

  $html_output = qq~
This is a multi-part message in MIME format.
--$boundary
Content-Type: text/html; charset=ISO-8859-1
Content-Transfer-Encoding: 7bit

<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
   <html><head><title>$project Cross Platform Report</title>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"></head>
    <body bgcolor="#FFFFCC" text="#000000">
    $header
    <br>
       <img alt="Collavorative Test System" src="cid:cts_logo_02212007" width="129" height="73">
    <br>
    <table width="562" border="0" cellpadding="0" cellspacing="0" mm:layoutgroup="true">
      <tr><td valign="top" height="57" nowrap colspan="3"> <p><img src="cid:cts_banner_02212007"></p></td></tr>
    </table>
    <table width="562" border="0" cellpadding="0" cellspacing="0" mm:layoutgroup="true">
      <tr>
        <td width="220" valign="top" height="73" nowrap>
        <p><b><font face="Verdana, Arial, Helvetica, sans-serif"><b><font color="#000066">
          $project Project
          </font></b></font></b><b><font face="Verdana, Arial, Helvetica, sans-serif"><b><br>
          <font color="#000066">$full_date_stamp<br></font></b></p>
        </td>
      </tr>
      <tr> <td height="12" width="4"></td> <td width="1"></td> </tr>
    </table>~; 
# end of print-to-string segment

# Print out the header with some adjustment for the number of columns.

  my $html_file = $self->harness->results_path . "/TestSummary.html";
  my $fh = FileHandle->new("> $html_file") or warn "Reporter::text: Could not create text_file\n";

  # foreach 
  foreach my $key(sort keys %itermax) {
    undef( @kv_results );
    $html_output .= sprintf("\n");
    $html_output .= sprintf("<p><p><table border=3 rules=\"all\">\n");

    $html_output .= sprintf("<tr>\n");
    my $iter=$itermax{$key}+1;
    $html_output .= sprintf("<td colspan=$iter align=\"center\">CROSS PLATFORM RESULTS ($key)</td>\n\n");
    $html_output .= sprintf("</tr>\n");

    my $num_fields = (exists($field5{$key}->[0]))? 5:
                     (exists($field4{$key}->[0]))? 4:
                     (exists($field3{$key}->[0]))? 3:
                     (exists($field2{$key}->[0]))? 2:
                                                   1;
    
    $html_output .= sprintf("<tr><td rowspan=$num_fields></td>\n");

    if (exists($field1{$key}->[0])){
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$html_output .= sprintf("<td>%s</td>",$field1{$key}->[$iter]);
      }
      $html_output .= sprintf("\n<tr>");
    }

    if (exists($field2{$key}->[0])){
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$html_output .= sprintf("<td>%s</td>",$field2{$key}->[$iter]);
      }
      $html_output .= sprintf("\n<tr>");
    }
    
    if (exists($field3{$key}->[0])){
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$html_output .= sprintf("<td>%s</td>",$field3{$key}->[$iter]);
      }
      $html_output .= sprintf("\n<tr>");
    }

    if (exists($field4{$key}->[0])){
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$html_output .= sprintf("<td>%s</td>",$field4{$key}->[$iter]);
      }
      $html_output .= sprintf("\n<tr>");
    }

    if (exists($field5{$key}->[0])){
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	$html_output .= sprintf("<td>%s</td>",$field5{$key}->[$iter]);
      }
      $html_output .= sprintf("\n<tr>");
    }
    $html_output .= sprintf("</tr>\n");
    
    my ($file, $prob_name, $result, $run_time);
    my %prob_data;
    my (@test_count, @test_pass, @total_run_time);
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      
      $file=$dir{$key}->[$iter];
      if (! open(RESULTSFILE, $file) ){
	print "Cannot open file[$file].\n";
	next;
      }
      $test_count[$iter]=0;
      $test_pass[$iter]=0;
      $total_run_time[$iter]=0;
      while ( defined( my $line = <RESULTSFILE> )){
	chomp $line;
	#     print "line is $line\n";
	$_ = $line;
	if (! /Test Result/){
	  if (/\b(Anomaly|Pass|Fail|PASSED|DIFF|FAILED|UNKNOWN|NOTAVAIL|SYSFAIL)\b/ ){
	    ($prob_name, $result, $run_time) = split(' ',$line);
	    
	    if ($run_time) {
	      $total_run_time[$iter] = $total_run_time[$iter] + $run_time;
	    }
	    
	    if ($result eq "PASSED" || $result eq "Pass" ) {
	      $result = "Pass";
	      $test_pass[$iter]++;
	      $prob_data{$prob_name}[$iter]="$result $run_time";
	    }
	    else {
	      $prob_data{$prob_name}[$iter]="$result";
	    }
	    
	    if ($result ne "NOTAVAIL" ) {
	      $test_count[$iter]++;
	    }
	    
	  }
          # read in misc key/value pairs
          if ( /(\S+)\s*=\s*(\S+.*?)\s*/ ){
            $kv_results[$iter]{$1} = $2;
          }
	}
      }
      
    }

    $html_output .= sprintf("<tr>");
    $html_output .= sprintf("<td>Totals</td>");

    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      $html_output .= sprintf("<td>%s</td>",$test_pass[$iter]."/".$test_count[$iter]);
    }
    $html_output .= sprintf("</tr>");
    
    $html_output .= sprintf("<tr>");
    $html_output .= sprintf("<td>RunTime</td>");
    
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      $html_output .= sprintf("<td>%s</td>",$total_run_time[$iter]);
    }
    $html_output .= sprintf("</tr>");
    $html_output .= sprintf("<td colspan=$iter align=\"left\">Test</td>");
    
    my $color="white";
    foreach my $prob (sort keys(%prob_data)){
      $html_output .= sprintf("<tr>\n");
      $html_output .= sprintf("  <td>%s</td>\n",$prob);
      for (my $iter=0; $iter<$itermax{$key}; $iter++){
	if (! defined($prob_data{$prob}[$iter])){
	  $prob_data{$prob}[$iter]="    ";
	}
	if ($prob_data{$prob}[$iter] =~  /Pass/) {
	  $color = $passed_color;
	}
	elsif ($prob_data{$prob}[$iter] eq "DIFF") {
	  $color = $diff_color;
	}
	elsif ($prob_data{$prob}[$iter] eq "NOTAVAIL") {
	  $color = $notavail_color;
	}
	elsif ($prob_data{$prob}[$iter] eq "SYSFAIL") {
	  $color = $sysfail_color;
	}
	elsif ($prob_data{$prob}[$iter] !~ /\S/) {
	  $color = $missing_color;
	}
	else{
	  $color = $failed_color;
	}
	$html_output .= sprintf("  <td bgcolor=$color>%s</td>\n",$prob_data{$prob}[$iter]);
      } # foreach testcase
      $html_output .= sprintf("</tr>\n");
    } # foreach prob
    
    $html_output .= sprintf("</table>\n");
    $html_output .= sprintf("</div>\n");
    $html_output .= sprintf("\n");
    $html_output .= sprintf("<p><p><table border=3 Cellspacing=\"3\" Cellpadding=\"10\" rules=\"all\">\n");
    $html_output .= sprintf("<tr>\n");
    $html_output .= sprintf("<th colspan=2 align=\"center\">Key=Value ($key)</th>\n\n");
    $html_output .= sprintf("</tr>\n");
    $html_output .= sprintf("<tr>\n");
    $html_output .= sprintf("<th colspan=2 align=\"center\"></th>\n\n");
    $html_output .= sprintf("</tr>\n");
    for (my $iter=0; $iter<$itermax{$key}; $iter++){
      my $fields_name = "";
      if (exists($field1{$key}->[0])){
        $fields_name .= ", $field1{$key}->[$iter]";
      }
      if (exists($field2{$key}->[0])){
        $fields_name .= ", $field2{$key}->[$iter]";
      }
      if (exists($field3{$key}->[0])){
        $fields_name .= ", $field3{$key}->[$iter]";
      }
      if (exists($field4{$key}->[0])){
        $fields_name .= ", $field4{$key}->[$iter]";
      }
      if (exists($field5{$key}->[0])){
        $fields_name .= ", $field5{$key}->[$iter]";
      }
      $fields_name =~ s/^, //;
      $html_output .= sprintf("<tr>\n");
      $html_output .= sprintf("<td></td><td colspan=1 align=\"left\">$fields_name</td>\n\n");
      $html_output .= sprintf("</tr>\n");
      foreach my $key_kv_results ( sort keys %{$kv_results[$iter]} )
        {
          $html_output .= sprintf("<tr>\n");
          $html_output .= sprintf("<td colspan=1 align=\"left\">$key_kv_results</td><td colspan=1 align=\"left\">$kv_results[$iter]{$key_kv_results}</td>\n\n");
          $html_output .= sprintf("</tr>\n");
        }
    }
    $html_output .= sprintf("</table>\n");
  }

  $html_output .= qq~
$footer
</body>
</html>

--$boundary
Content-Type: image/jpeg;
 name="cts_banner.jpg"
Content-Transfer-Encoding: base64
Content-ID: <cts_banner_02212007>
Content-Disposition: inline;
 filename="cts_banner.jpg"
~;
  $html_output .= $banner;
  $html_output .= qq~

--$boundary
Content-Type: image/jpeg;
 name="cts_logo.jpg"
Content-Transfer-Encoding: base64
Content-ID: <cts_logo_02212007>
Content-Disposition: inline;
 filename="cts_logo.jpg"
~;
  $html_output .= $logob;
  $html_output .= qq~

--$boundary--

-- End --
~;
  # print the results to file and stdout
  if (defined $fh) {
    print $fh $html_output;
    close $fh;
  }
  # print $html_output;

  $html_file = $CROSS_PLATFORM_BASE . "/TestSummary.html";
  if (-f $html_file) {
    unlink($html_file);
  }
  $fh = FileHandle->new("> $html_file") or warn "Reporter::html: Could not create html_file\n";
  if (defined $fh) {
    print $fh $html_output;
    system "cp $logo $CROSS_PLATFORM_BASE";
    $logo = basename $logo;
    system "cp $FindBin::Bin/../html/cts_banner.jpg $CROSS_PLATFORM_BASE";
    close $fh;
    if (defined($CROSS_PLATFORM_GROUP)){
      `chgrp $CROSS_PLATFORM_GROUP $html_file`;
      `chgrp $CROSS_PLATFORM_GROUP $CROSS_PLATFORM_BASE/$logo`;
      `chgrp $CROSS_PLATFORM_GROUP $CROSS_PLATFORM_BASE/cts_banner.jpg`;
    }
    if (defined($CROSS_PLATFORM_MODE)){
      `chmod $CROSS_PLATFORM_MODE $html_file`;
      `chmod $CROSS_PLATFORM_MODE $CROSS_PLATFORM_BASE/$logo`;
      `chmod $CROSS_PLATFORM_MODE $CROSS_PLATFORM_BASE/cts_banner.jpg`;
    }
  }

  if (defined($CROSS_PLATFORM_MAIL)) {
    my $smpt = $CROSS_PLATFORM_MAIL;
    $smpt =~ s/[-\w]+@//;
    my $from = $ENV{LOGNAME} . '@' . $smpt;
    $smpt = "mail." . $smpt;
    unshift @{$Mail::Sendmail::mailcfg{'smtp'}} , $smpt;
    my $cronjob = "";
    if( defined( $ENV{CRONJOB}) )
      {
        $cronjob = ": $ENV{CRONJOB}";
      }
    my %mail = (
		  To      => $CROSS_PLATFORM_MAIL,
		  From    => $from,
		  Subject => "CTS Tests: ".$self->harness->project.": ".
                             $self->harness->system_name.
                             $cronjob,
		  "Content-Type" => "multipart/related; boundary=\"$boundary\"",
		  Message => $html_output,
	        );

    sendmail(%mail) or die "Reporter:: Sendmail error: $!";
  }
}


sub reports{
  my $self = shift;

  my $reports = shift || $self->{REPORTS} || [];
  if (ref $reports eq "ARRAY") {
    $self->{REPORTS} = $reports;
  }
  else {
    carp "TestContainer::reports: Bad value in reports command, expected an array ref but received $reports\n" 
  }
  return $reports;
}

sub header{
  my $self = shift;

  my $header = shift || $self->{HEADER} || "";
  if (ref $header eq "ARRAY") {
    $header = $header->[0];
  }
  $self->{HEADER} = $header;
  return $header;
}

sub footer{
  my $self = shift;

  my $footer = shift || $self->{FOOTER} || "";
  if (ref $footer eq "ARRAY") {
    $footer = $footer->[0];
  }
  $self->{FOOTER} = $footer;
  return $footer;
}


sub _init_cross_platform {
  my $self = shift;

  if (defined($self->CROSS_PLATFORM_REPORTER)) {
    $CROSS_PLATFORM_REPORTER = @{$self->CROSS_PLATFORM_REPORTER}[0];
  }
  if (defined($self->CROSS_PLATFORM_BASE)) {
    $CROSS_PLATFORM_BASE     = @{$self->CROSS_PLATFORM_BASE}[0];
  }
  if (defined($self->CROSS_PLATFORM_PATTERN)) {
    $CROSS_PLATFORM_PATTERN  = @{$self->CROSS_PLATFORM_PATTERN}[0];
  }
  if (defined($self->CROSS_PLATFORM_MODE)) {
    $CROSS_PLATFORM_MODE     = @{$self->CROSS_PLATFORM_MODE}[0];
  }
  if (defined($self->CROSS_PLATFORM_GROUP)) {
    $CROSS_PLATFORM_GROUP  = @{$self->CROSS_PLATFORM_GROUP}[0];
  }
  if (defined($self->CROSS_PLATFORM_MAIL)) {
    $CROSS_PLATFORM_MAIL = @{$self->CROSS_PLATFORM_MAIL}[0];
  }
  if (defined($self->CROSS_PLATFORM_SORT)) {
    $CROSS_PLATFORM_SORT = @{$self->CROSS_PLATFORM_SORT}[0];
  }
}

sub _init_cross_platform_report {
  my $self = shift;

  $self->_init_cross_platform();

  if (! defined($CROSS_PLATFORM_BASE) || ! defined($CROSS_PLATFORM_PATTERN)) {
    return undef;
  }

  my $MATCH_STRING="$CROSS_PLATFORM_BASE/$CROSS_PLATFORM_PATTERN";
  my $EXTRACT_STRING=$MATCH_STRING;


# MATCH STRING is used to wildcard match the files and
# EXTRACT STRING is used to extract the keywords to be
# used in the column headings.

  $MATCH_STRING =~ s/\$\{system\}/\*/;
  $MATCH_STRING =~ s/\$\{Fcompiler\}/\*/;
  $MATCH_STRING =~ s/\$\{Ccompiler\}/\*/;
  $MATCH_STRING =~ s/\$\{mpi\}/\*/;
  $MATCH_STRING =~ s/\$\{type\}/\*/;
  #print "DEBUG: MATCH STRING is $MATCH_STRING\n";

  #print "DEBUG: EXTRACT STRING is $EXTRACT_STRING\n";
  $EXTRACT_STRING =~ s/\\/\/\\/;
  $EXTRACT_STRING =~ s/\./\\./;
  $EXTRACT_STRING =~ s/\*/.*/;
  $EXTRACT_STRING =~ s/\$\{system\}/([^_\/]*)/;
  $EXTRACT_STRING =~ s/\$\{Fcompiler\}/([^_\/]*)/;
  $EXTRACT_STRING =~ s/\$\{Ccompiler\}/([^_\/]*)/;
  $EXTRACT_STRING =~ s/\$\{mpi\}/([^_\/]*)/;
  $EXTRACT_STRING =~ s/\$\{type\}/([^_\/]*)/;
  #print "DEBUG: EXTRACT STRING is $EXTRACT_STRING\n\n";
  my %iter;
  my @system_inquire_doc=`ls $MATCH_STRING`;
  grep( s/^\s*//, @system_inquire_doc );
  grep( chomp, @system_inquire_doc );
  # sort on another field for reporting if requested
  if( defined($CROSS_PLATFORM_SORT) ){
    # create new @sorted which is <sort field>_<line>
    # and then sort on that new line.
    # net result is that you have array sorted first on the sort field
    my @sorted;
    foreach my $line ( @system_inquire_doc ) {
      my @fields;
      $line =~ m/$EXTRACT_STRING/;
      push( @fields, ($1, $2, $3, $4, $5) );
      if( ! defined($fields[$CROSS_PLATFORM_SORT]) ){
         die "Reporter:: _init_cross_platform_report:\n  CROSS_PLATFORM_FIELD [$CROSS_PLATFORM_SORT] does not point to existing field in\n  name [$line] " 
      }
      $line = "$fields[$CROSS_PLATFORM_SORT]_$line";
      push( @sorted, $line );
    }
    @sorted = sort @sorted;
    grep ( s/^[^_]*_//, @sorted );
    @system_inquire_doc = @sorted;
  }
  for my $line (@system_inquire_doc){
    my $key = basename $line;
    $key =~ s/\..*$//;
    if (! exists $iter{$key}) {
      $iter{$key} = 0;
    }
    $dir{$key}->[$iter{$key}]=$line;
    #    print "$line\n";
    $line =~ m/$EXTRACT_STRING/;
    #print "DEBUG: $line\n";
    #print "DEBUG: $1: $2: $3: $4:\n\n";
    if (defined($1)){
      $field1{$key}->[$iter{$key}]=$1;
    }
    if (defined($2)){
      $field2{$key}->[$iter{$key}]=$2;
    }
    if (defined($3)){
      $field3{$key}->[$iter{$key}]=$3;
    }
    if (defined($4)){
      $field4{$key}->[$iter{$key}]=$4;
    }
    if (defined($5)){
      $field5{$key}->[$iter{$key}]=$5;
    }
    $iter{$key}++;
    $itermax{$key}=$iter{$key};
  }

  return 1;
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Reporter - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Reporter;
  blah blah blah

=head1 ABSTRACT

  This should be the abstract for Reporter.
  The abstract is used when making PPD (Perl Package Description) files.
  If you don't want an ABSTRACT you should also edit Makefile.PL to
  remove the ABSTRACT_FROM option.

=head1 DESCRIPTION

Stub documentation for Reporter, created by h2xs. It looks like the
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
