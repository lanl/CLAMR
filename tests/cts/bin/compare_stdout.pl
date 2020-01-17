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
      if (/     or  /) {next LINE;}
      if (/Detected/) {next LINE;}
      if (/Package 0/) {next LINE;}
      if (/Package 1/) {next LINE;}
      if (/sysfs power/) {next LINE;}
      if (/RAPL_SYSFS/) {next LINE;}
      if (/Energy Consumed/) {next LINE;}
      if (/package-0/) {next LINE;}
      if (/package-1/) {next LINE;}
      if (/dram/) {next LINE;}
      if (/high-performance Open MPI point-to-point messaging/) {next LINE;}
      if (/unable to find any relevant network/) {next LINE;}
      if (/Module: OpenFabrics/) {next LINE;}
      if (/Host:/) {next LINE;}
      if (/Another transport will be used instead/) {next LINE;}
      if (/lower performance/) {next LINE;}
      if (/sent help message/) {next LINE;}
      if (/Set MCA parameter/) {next LINE;}
      if (/You can disable this warning/) {next LINE;}
      if (/btl_base_warn_component/) {next LINE;}
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
      if (/^CRUX checkpointing time averaged/) {next LINE;}
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
      if (/     or  /) {next LINE;}
      if (/Detected/) {next LINE;}
      if (/Package 0/) {next LINE;}
      if (/Package 1/) {next LINE;}
      if (/sysfs power/) {next LINE;}
      if (/RAPL_SYSFS/) {next LINE;}
      if (/Energy Consumed/) {next LINE;}
      if (/package-0/) {next LINE;}
      if (/package-1/) {next LINE;}
      if (/dram/) {next LINE;}
      if (/high-performance Open MPI point-to-point messaging/) {next LINE;}
      if (/unable to find any relevant network/) {next LINE;}
      if (/Module: OpenFabrics/) {next LINE;}
      if (/Host:/) {next LINE;}
      if (/Another transport will be used instead/) {next LINE;}
      if (/lower performance/) {next LINE;}
      if (/sent help message/) {next LINE;}
      if (/Set MCA parameter/) {next LINE;}
      if (/You can disable this warning/) {next LINE;}
      if (/btl_base_warn_component/) {next LINE;}
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
      if (/^CRUX checkpointing time averaged/) {next LINE;}
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



