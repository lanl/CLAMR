#!/bin/sh
set -x
FILE_LIST="origcellrun.out origfacerun.out cellinplacerun.out faceinplacerun.out reggridbycellrun.out reggridbyfacesrun.out"
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "cell"                   >& origcellrun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "face"                   >& origfacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "cell-in-place"          >& cellinplacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "face-in-place"          >& faceinplacerun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "regular-grid"           >& reggridbycellrun.out
./clamr_cpuonly -n 128 -t 500 -i 100 -l 2 -A "regular-grid-by-faces"  >& reggridbyfacesrun.out
fgrep "Total CPU" ${FILE_LIST} | sed 's/run.out:Profiling: Total CPU          time was/ Total CPU/' 
fgrep "state_timer_finite_difference" ${FILE_LIST} | sed 's/run.out:CPU:   state_timer_finite_difference/ finite_difference/'

fgrep "Memory used" ${FILE_LIST} | grep -v "in startup" | sed 's/run.out:Memory used/ Memory_used/'
#./clamr_cpuonly -A "regular-cell-by-faces" -t 100 -i 10 -l 2

