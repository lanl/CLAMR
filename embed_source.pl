eval 'exec perl -w $0 ${1+"$@"}'
  if 0;

use File::Basename;

$ARGV_LAST = "";

($file,,) = fileparse($ARGV[0], qr/\..*/);
print "const char *",$file,"_source =\n";

while(<>) {
   if ($ARGV ne $ARGV_LAST){
#      $comment = 0;
      $ARGV_LAST=$ARGV;
   }
   chop $_;

   s/"/""/g;

#   if (/^\/\*/){
#      #print "Setting comment to 1\n";
#      $comment = 1;
#   }
#   if (/\*\/$/){
#      #print "Setting comment to 0\n";
#      $comment = 0;
#      next;
#   }
#   if ($comment eq 1) {
#      #print "Skip line and goto next\n";
##      next;
#   }

   print '"',$_,'\n','"',"\n";
}

print ";"
