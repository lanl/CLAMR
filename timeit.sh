#!/bin/sh
set -x
#./clamr_cpuonly -A "cell" -n 256 -t 4000 -i 100 -l 2 >& cellrun.out
#./clamr_cpuonly -A "face" -n 256 -t 4000 -i 100 -l 2 >& facerun.out
#./clamr_cpuonly -A "cell-in-place" -n 256 -t 4000 -i 100 -l 2 >& cellinplacerun.out
#./clamr_cpuonly -A "face-in-place" -n 256 -t 4000 -i 100 -l 2 >& faceinplacerun.out
fgrep "Total CPU" *run.out
fgrep "state_timer_finite_difference" *run.out
#./clamr_cpuonly -A "regular-cell" -t 100 -i 10 -l 2
#./clamr_cpuonly -A "regular-cell-by-faces" -t 100 -i 10 -l 2

