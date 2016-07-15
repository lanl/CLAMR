eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

# waits until job is done
$arg_string = "";
$qsub = "qsub";
foreach $arg ( @ARGV )
  {
    $arg_string = "$arg_string \"$arg\"";
  }
$arg_string =~ s/^\s+//;
$output = `$qsub $arg_string 2>&1`;
print $output;
if( $output =~ /\b(\d+\.nid\d+)\b/ )
  {
    $id = $1;
    $ierr = 0;
    $done = "false";
    while( $done eq "false" )
      {
        sleep( 60 );
        $output = `qstat $id 2>&1`;
        if( $output !~ /Job id\s+Name\s+User/ || $output =~ / C \S+\s*$/ )
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
