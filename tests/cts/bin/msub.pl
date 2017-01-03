eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

# wait until job is done
$arg_string = "";
$msub = "msub";
foreach $arg ( @ARGV )
  {
    $arg_string .= " $arg";
  }
$arg_string =~ s/^\s+//;
print "$msub $arg_string 2>&1";
$output = `$msub $arg_string 2>&1`;
print "$output\n";
if( $? == 0 && $output =~ /^\s*(\d+)/m )
  {
    $id = $1;
    $ierr = 0;
    $done = "false";
    while( $done eq "false" )
      {
        sleep( 60 );
        $command = "mdiag -j $id";
	#print "Command: $command\n";
        $output = `$command 2>&1`;
        #print ":$output:\n";
        if(( $output !~ /^\s*$id\s/m ) || ($output =~ /^\s*$id\s+complete/im ))
          {
            $done = "true";
          }
      }
  }
else
  {
    $ierr = 1;
  }
exit( $ierr );
