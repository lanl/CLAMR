#!/bin/sh
set -x
./clamr_cpuonly -n 256 -t 1000 -i 100 -l 2 -A "cell"                   >& cellrun.out
./clamr_cpuonly -n 256 -t 1000 -i 100 -l 2 -A "face"                   >& facerun.out
./clamr_cpuonly -n 256 -t 1000 -i 100 -l 2 -A "cell-in-place"          >& cellinplacerun.out
./clamr_cpuonly -n 256 -t 1000 -i 100 -l 2 -A "face-in-place"          >& faceinplacerun.out
./clamr_cpuonly -n 256 -t 1000 -i 100 -l 2 -A "regular-cell"           >& regcellrun.out
./clamr_cpuonly -n 256 -t 1000 -i 100 -l 2 -A "regular-cell-by-faces"  >& regcellrun.out
fgrep "Total CPU" *run.out
fgrep "state_timer_finite_difference" *run.out
#./clamr_cpuonly -A "regular-cell-by-faces" -t 100 -i 10 -l 2

