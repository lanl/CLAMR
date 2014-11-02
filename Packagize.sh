#!/bin/sh
set -v
EXCLUDES="--exclude */._* --exclude */*/._* --exclude .DS_Store"

tar ${EXCLUDES} -czf ~/Packages/MallocPlus_v2.0.7.tgz MallocPlus
tar ${EXCLUDES} -czf ~/Packages/PowerParser_v2.0.7.tgz PowerParser
tar ${EXCLUDES} -czf ~/Packages/crux_v2.0.7.tgz crux
tar ${EXCLUDES} -czf ~/Packages/ezcl_v2.0.7.tgz ezcl
tar ${EXCLUDES} -czf ~/Packages/genmalloc_v2.0.7.tgz genmalloc
tar ${EXCLUDES} -czf ~/Packages/graphics_v2.0.7.tgz graphics
tar ${EXCLUDES} -czf ~/Packages/hash_v2.0.7.tgz hash
tar ${EXCLUDES} -czf ~/Packages/l7_v2.0.7.tgz l7
tar ${EXCLUDES} -czf ~/Packages/memstats_v2.0.7.tgz memstats
tar ${EXCLUDES} -czf ~/Packages/mesh_v2.0.7.tgz mesh
tar ${EXCLUDES} -czf ~/Packages/s7_v2.0.7.tgz s7
tar ${EXCLUDES} -czf ~/Packages/timer_v2.0.7.tgz timer

cd ../; tar -czf ~/Packages/CLAMR_v2.0.7.tgz \
  CLAMR/AUTHORS \
  CLAMR/AdditionalInfo.txt \
  CLAMR/CMakeLists.txt \
  CLAMR/ChangeLog \
  CLAMR/FreeSansBold.ttf \
  CLAMR/LICENSE \
  CLAMR/PAPERS \
  CLAMR/Papers.bib \
  CLAMR/README \
  CLAMR/README.md \
  CLAMR/TODO \
  CLAMR/clamr.1 \
  CLAMR/clamr.cpp \
  CLAMR/clamr_checkall.cpp \
  CLAMR/clamr_cpuonly.cpp \
  CLAMR/clamr_gpucheck.cpp \
  CLAMR/clamr_gpuonly.cpp \
  CLAMR/clamr_mpicheck.cpp \
  CLAMR/clamr_mpionly.cpp \
  CLAMR/clamr_quo.cpp \
  CLAMR/cmake \
  CLAMR/config.h \
  CLAMR/config.h.in \
  CLAMR/embed_source.pl \
  CLAMR/fileamr.in \
  CLAMR/generate_image.py \
  CLAMR/generate_movie.sh \
  CLAMR/input.cpp \
  CLAMR/input.h \
  CLAMR/moonlight_scripts \
  CLAMR/newbatch.sh \
  CLAMR/newbatch_check.sh \
  CLAMR/runlongtests.sh \
  CLAMR/runtests.sh \
  CLAMR/state.cpp \
  CLAMR/state.h \
  CLAMR/state_kern.cl

#spack create -f file:/Users/lbrobey/Packages/MallocPlus_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/PowerParser_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/crux_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/ezcl_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/genmalloc_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/graphics_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/hash_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/l7_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/memstats_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/mesh_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/s7_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/timer_v2.0.7.tgz 
#spack create -f file:/Users/lbrobey/Packages/CLAMR_v2.0.7.tgz 

