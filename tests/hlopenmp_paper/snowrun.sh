#!/bin/sh -x
export TEST_NAME="n1920lev1"
export TEST_OPTIONS="-n 1920 -l 1 -i 10 -t 100"

for i in {1,2,4,8,16}
do
   mpirun -n $i ../../clamr_mpionly ${TEST_OPTIONS} >& ${TEST_NAME}_parallel$i.out
done

for i in {1,2,4,8,16}
do
   export OMP_NUM_THREADS=$i
   ../../clamr_openmponly ${TEST_OPTIONS} >& ${TEST_NAME}_openmp$i.out
done


../../clamr_cpuonly ${TEST_OPTIONS} >& ${TEST_NAME}_serial.out

rm -f clamr_plot.dat
for i in {1,2,4,8,16}
do
  echo -n $i >> clamr_plot.dat
  awk '/Total CPU/ { printf " %f ", $6}' ${TEST_NAME}_parallel$i.out >> clamr_plot.dat
  awk '/Total CPU/ { printf " %f ", $6}' ${TEST_NAME}_openmp$i.out >> clamr_plot.dat
  echo -n $i >> clamr_plot.dat
  echo  >> clamr_plot.dat
done

python clamr_plot.py clamr_plot.dat
