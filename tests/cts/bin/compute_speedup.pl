eval 'exec perl -w $0 ${1+"$@"}'
  if 0;

$parallelproc= $ARGV[0];
$serialfile  = $ARGV[1];
$parallelfile= $ARGV[2];
shift;
shift;

# Read reference file first to allow time for new file to be
# written on slow file system

if (open (SERIALFH, $serialfile) ){
   $str = <SERIALFH>;
   ($dummy, $dummy, $dummy, $dummy, $time1, $dummy) = split(' ', $str);
   close SERIALFH;
}
else {
   print "Error: Could not open $serialfile -- FAIL $arg\n";
}

if (open (PARALLELFH, $parallelfile) ){
   $str = <PARALLELFH>;
   ($dummy, $dummy, $dummy, $dummy, $time2, $dummy) = split(' ', $str);
   close PARALLELFH;
}
else {
   print "Error: Could not open $parallelfile -- FAIL $arg\n";
}

if ($time2 != 0.0) {
   printf "SPEEDUP is %6.3lf on %d processors\n", $time1/$time2, $parallelproc
} else {
   printf "SPEEDUP can't be calculated for zero run time\n"
}
