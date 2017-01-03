eval 'exec perl -w -S $0 ${1+"$@"}'
  if 0;

# waits until job is done
$arg_string = "";
$psub = "psub";
foreach $arg ( @ARGV )
  {
    $arg_string = "$arg_string \"$arg\"";
  }
$arg_string =~ s/^\s+//;
$output = `$psub $arg_string 2>&1`;
print $output;
if( $output =~ /\bJob\s+(\d+)\s+submitted\b/ )
  {
    $id = $1;
    $ierr = 0;
    $done = "false";
    while( $done eq "false" )
      {
        sleep( 60 );
        $output = `/usr/local/bin/pstat $id 2>&1`;
        if( $output =~ /not\s+found/ || $output =~ / C \S+\s*$/ )
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
