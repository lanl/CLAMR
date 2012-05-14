bin_PROGRAMS = clamr_gpuonly clamr_cpuonly clamr_mpionly clamr_gpucheck clamr_checkall
#bin_PROGRAMS = clamr clamr_checkall clamr_gpuonly

SUBDIRS = l7 kdtree hsfc s7 timer ezcl zorder

#clamr_CXX_SRCS = clamr.cpp mesh.cpp partition.cpp state.cpp input.cpp
#clamr_C_SRCS = display.c reduce.c
#clamr_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h

#clamr_SOURCES = ${clamr_CXX_SRCS} ${clamr_C_SRCS} ${clamr_H_SRCS}

clamr_gpuonly_CXX_SRCS = clamr_gpuonly.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_gpuonly_C_SRCS = display.c reduce.c
clamr_gpuonly_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h

clamr_gpuonly_SOURCES = ${clamr_gpuonly_CXX_SRCS} ${clamr_gpuonly_C_SRCS} ${clamr_gpuonly_H_SRCS}

clamr_cpuonly_CXX_SRCS = clamr_cpuonly.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_cpuonly_C_SRCS = display.c reduce.c
clamr_cpuonly_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h

clamr_cpuonly_SOURCES = ${clamr_cpuonly_CXX_SRCS} ${clamr_cpuonly_C_SRCS} ${clamr_cpuonly_H_SRCS}

clamr_mpionly_CXX_SRCS = clamr_mpionly.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_mpionly_C_SRCS = display.c reduce.c
clamr_mpionly_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h

clamr_mpionly_SOURCES = ${clamr_mpionly_CXX_SRCS} ${clamr_mpionly_C_SRCS} ${clamr_mpionly_H_SRCS}

clamr_gpucheck_CXX_SRCS = clamr_gpucheck.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_gpucheck_C_SRCS = display.c reduce.c
clamr_gpucheck_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h

clamr_gpucheck_SOURCES = ${clamr_gpucheck_CXX_SRCS} ${clamr_gpucheck_C_SRCS} ${clamr_gpucheck_H_SRCS}

clamr_checkall_CXX_SRCS = clamr_checkall.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_checkall_C_SRCS = display.c reduce.c 
clamr_checkall_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h

clamr_checkall_SOURCES = ${clamr_checkall_CXX_SRCS} ${clamr_checkall_C_SRCS} ${clamr_checkall_H_SRCS}

LDADD = -Lhsfc -lhsfc -Lkdtree -lkdtree -Lzorder -lzorder -Ll7 -ll7 -Ls7 -ls7 -Lezcl -lezcl -Ltimer -ltimer -lm

EXTRA_DIST = wave_kern.cl wave_kern_calc.cl reduce.cl config.h clamr.xcodeproj/project.pbxproj clamr.1 fileamr.in AUTHORS ChangeLog runtests.sh newbatch.sh LICENSE

distclean-local:
	rm -f config.h.in~
	rm -rf autom4te.cache
	rm -rf config.log config.status
