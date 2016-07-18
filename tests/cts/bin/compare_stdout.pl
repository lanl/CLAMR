eval 'exec perl -w $0 ${1+"$@"}'
  if 0;

$tolerance = $ARGV[0];
$reffile   = $ARGV[1];
$newfile   = $ARGV[2];
shift;
shift;
shift;
$CTS_BIN=`dirname $0`;
chop $CTS_BIN;

$string = $newfile;
$string =~ s/.out//;

if (`ls core* 2>>/dev/null |fgrep -c core` > 0){
   print "CORE -- $string\n";
   exit(1);
}

# Read reference file first to allow time for new file to be
# written on slow file system

open (OUTFH, ">temp_gold") || die "Can't open temp_gold file";

if (open (GOLDFH, $reffile) ){
   LINE:
   while ($_ = <GOLDFH>){
      if (/^App launch reported/) {next LINE;}
      if (/^Reported/) {next LINE;}
      if (/^--- max num openmp threads/) {next LINE;}
      if (/^--- num openmp threads in parallel region/) {next LINE;}
      if (/^Starting compile/) {next LINE;}
      if (/^Finishing compile/) {next LINE;}
      if (/^CPU: /) {next LINE;}
      if (/^GPU: /) {next LINE;}
      if (/^Memory used/) {next LINE;}
      if (/^Memory peak/) {next LINE;}
      if (/^Memory free/) {next LINE;}
      if (/^Memory available/) {next LINE;}
      if (/^Profiling/) {next LINE;}
      if (/^Mesh Ops/) {next LINE;}
      if (/^The MPI surface area/) {next LINE;}
      if (/^hash table size/) {next LINE;}
      if (/^Final hash collision report/) {next LINE;}
      if (/^Parallel Speed-up/) {next LINE;}
      if (/Device timing information/) {next LINE;}
      if (/^---------------/) {next LINE;}
      if (/^===============/) {next LINE;}
      if (/^Iteration/) {
         my ($str1) = split /Mass Change/;
         print OUTFH "$str1\n";
         next LINE;
      }
      print OUTFH $_;
   }
   close GOLDFH;
}
else {
   print "Error: Could not open $reffile -- FAIL $arg\n";
   exit(1);
}
close(OUTFH);

open (OUTFH, ">temp_new") || die "Can't open temp_new file";

if (open (NEWFH, $newfile) ){
   LINE:
   while ($_ = <NEWFH>){
      if (/^App launch reported/) {next LINE;}
      if (/^Reported/) {next LINE;}
      if (/^--- max num openmp threads/) {next LINE;}
      if (/^--- num openmp threads in parallel region/) {next LINE;}
      if (/^Starting compile/) {next LINE;}
      if (/^Finishing compile/) {next LINE;}
      if (/^CPU: /) {next LINE;}
      if (/^GPU: /) {next LINE;}
      if (/^Memory used/) {next LINE;}
      if (/^Memory peak/) {next LINE;}
      if (/^Memory free/) {next LINE;}
      if (/^Memory available/) {next LINE;}
      if (/^Profiling/) {next LINE;}
      if (/^Mesh Ops/) {next LINE;}
      if (/^The MPI surface area/) {next LINE;}
      if (/^hash table size/) {next LINE;}
      if (/^Final hash collision report/) {next LINE;}
      if (/^Parallel Speed-up/) {next LINE;}
      if (/Device timing information/) {next LINE;}
      if (/^---------------/) {next LINE;}
      if (/^===============/) {next LINE;}
      if (/^Iteration/) {
         my ($str1) = split /Mass Change/;
         print OUTFH "$str1\n";
         next LINE;
      }
      print OUTFH $_;
   }
   close NEWFH;
}
else {
   print "Error: Could not open $newfile -- FAIL $arg\n";
   exit(1);
}
close(OUTFH);

$output=`$CTS_BIN/cts_diff.pl -presult -v 4 -no_plots -a $tolerance -or -r $tolerance temp_gold temp_new`;
print $output;
open FULL_OUTPUT, ">FULL_RESULTS";
print FULL_OUTPUT $output;
close FULL_OUTPUT;



