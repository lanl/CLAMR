#!/bin/sh
set -x
rm -f cut1000
./clamr_cpuonly -Z "all" -n 128 -t 1000 -i 1000 -l 2 -A "cell"                   >& /dev/null
sort -n cut1000 > cell.cut1000
rm -f cut1000
./clamr_cpuonly -Z "all" -n 128 -t 1000 -i 1000 -l 2 -A "face"                   >& /dev/null
sort -n cut1000 > face.cut1000
rm -f cut1000
./clamr_cpuonly -Z "all" -n 128 -t 1000 -i 1000 -l 2 -A "cell-in-place"          >& /dev/null
sort -n cut1000 > cellinplace.cut1000
rm -f cut1000
./clamr_cpuonly -Z "all" -n 128 -t 1000 -i 1000 -l 2 -A "face-in-place"          >& /dev/null
sort -n cut1000 > faceinplace.cut1000
rm -f cut1000
./clamr_cpuonly -Z "all" -n 128 -t 1000 -i 1000 -l 2 -A "regular-grid"           >& /dev/null
sort -n cut1000 > reggridbycell.cut1000
rm -f cut1000
./clamr_cpuonly -Z "all" -n 128 -t 1000 -i 1000 -l 2 -A "regular-grid-by-faces"  >& /dev/null
sort -n cut1000 > reggridbyface.cut1000
rm -f cut1000

