#!/bin/sh
set -x
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "cell"                   >& cellrun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "face"                   >& facerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "cell-in-place"          >& cellinplacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "face-in-place"          >& faceinplacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "regular-cell"           >& regcellrun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "regular-cell-by-faces"  >& regcellbyfacesrun.out
fgrep "Total CPU" *run.out
fgrep "state_timer_finite_difference" *run.out

fgrep "Memory used" *run.out | grep -v "in startup"
#./clamr_cpuonly -A "regular-cell-by-faces" -t 100 -i 10 -l 2

