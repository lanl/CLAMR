#!/bin/sh
set -x
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "cell"                   >& cellrun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "face"                   >& facerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "cell-in-place"          >& cellinplacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "face-in-place"          >& faceinplacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "regular-cell"           >& regcellrun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "regular-cell-by-faces"  >& regcellbyfacesrun.out
fgrep "Total CPU" *run.out | sed 's/run.out:Profiling: Total CPU          time was/ Total CPU/' 
fgrep "state_timer_finite_difference" *run.out | sed 's/run.out:CPU:   state_timer_finite_difference/ finite_difference/'

fgrep "Memory used" *run.out | grep -v "in startup" | sed 's/run.out:Memory used/ Memory_used/'
#./clamr_cpuonly -A "regular-cell-by-faces" -t 100 -i 10 -l 2

