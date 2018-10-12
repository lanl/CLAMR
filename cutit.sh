#!/bin/sh
set -x
rm -f cut2000
./clamr_cpuonly -Z -n 128 -t 2000 -i 2000 -l 2 -A "cell"                   >& /dev/null
mv cut2000 cell.cut2000
./clamr_cpuonly -Z -n 128 -t 2000 -i 2000 -l 2 -A "face"                   >& /dev/null
mv cut2000 face.cut2000
./clamr_cpuonly -Z -n 128 -t 2000 -i 2000 -l 2 -A "cell-in-place"          >& /dev/null
mv cut2000 cellinplace.cut2000
./clamr_cpuonly -Z -n 128 -t 2000 -i 2000 -l 2 -A "face-in-place"          >& /dev/null
mv cut2000 faceinplace.cut2000
./clamr_cpuonly -Z -n 128 -t 2000 -i 2000 -l 2 -A "regular-grid"           >& /dev/null
mv cut2000 reggridbycell.cut2000
./clamr_cpuonly -Z -n 128 -t 2000 -i 2000 -l 2 -A "regular-grid-by-faces"  >& /dev/null
mv cut2000 reggridbyface.cut2000

