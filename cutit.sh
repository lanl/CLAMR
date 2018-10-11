#!/bin/sh
set -x
rm -f cut500
./clamr_cpuonly -Z "yaxis" -n 128 -t 500 -i 500 -l 2 -A "cell"                   >& /dev/null
mv cut500 cell.cut500
./clamr_cpuonly -Z "yaxis" -n 128 -t 500 -i 500 -l 2 -A "face"                   >& /dev/null
mv cut500 face.cut500
./clamr_cpuonly -Z "yaxis" -n 128 -t 500 -i 500 -l 2 -A "cell-in-place"          >& /dev/null
mv cut500 cellinplace.cut500
./clamr_cpuonly -Z "yaxis" -n 128 -t 500 -i 500 -l 2 -A "face-in-place"          >& /dev/null
mv cut500 faceinplace.cut500
./clamr_cpuonly -Z "yaxis" -n 128 -t 500 -i 500 -l 2 -A "regular-grid"           >& /dev/null
mv cut500 reggridbycell.cut500
./clamr_cpuonly -Z "yaxis" -n 128 -t 500 -i 500 -l 2 -A "regular-grid-by-faces"  >& /dev/null
mv cut500 reggridbyface.cut500

