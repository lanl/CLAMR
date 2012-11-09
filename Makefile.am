bin_PROGRAMS = clamr clamr_gpuonly clamr_cpuonly clamr_mpionly clamr_gpucheck clamr_mpicheck clamr_checkall

SUBDIRS = l7 kdtree hsfc s7 timer ezcl zorder genmalloc

AM_MAKEFLAGS = -j

ACLOCAL_AMFLAGS = -I config

# -all-static flag is to make debugging easier. Removing it will enable dynamic libraries instead

clamr_CXX_SRCS = clamr.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_C_SRCS = display.c reduce.c
clamr_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h
clamr_CPPFLAGS = -DHAVE_MPI=1 -DHAVE_OPENCL=1
clamr_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la ezcl/libpezcl.la timer/libtimer.la l7/libl7.la genmalloc/libgenmalloc.la -lm ${MPE_LIBS} ${MPI_LIBS} ${OPENCL_LIBS}

clamr_SOURCES = ${clamr_CXX_SRCS} ${clamr_C_SRCS} ${clamr_H_SRCS}

clamr_gpuonly_CXX_SRCS = clamr_gpuonly.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_gpuonly_C_SRCS = display.c reduce.c
clamr_gpuonly_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h
clamr_gpuonly_CPPFLAGS = -UHAVE_MPI -DHAVE_OPENCL=1
clamr_gpuonly_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la ezcl/libezcl.la timer/libtimer.la genmalloc/libgenmalloc.la -lm ${MPE_NOMPI_LIBS} ${OPENCL_LIBS}

clamr_gpuonly_SOURCES = ${clamr_gpuonly_CXX_SRCS} ${clamr_gpuonly_C_SRCS} ${clamr_gpuonly_H_SRCS}

clamr_cpuonly_CXX_SRCS = clamr_cpuonly.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_cpuonly_C_SRCS = display.c
clamr_cpuonly_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h
clamr_cpuonly_CPPFLAGS = -UHAVE_MPI -UHAVE_OPENCL
clamr_cpuonly_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la timer/libtimer.la genmalloc/libgenmalloc.la -lm ${MPE_NOMPI_LIBS}

clamr_cpuonly_SOURCES = ${clamr_cpuonly_CXX_SRCS} ${clamr_cpuonly_C_SRCS} ${clamr_cpuonly_H_SRCS}

clamr_mpionly_CXX_SRCS = clamr_mpionly.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_mpionly_C_SRCS = display.c
clamr_mpionly_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h
clamr_mpionly_CPPFLAGS = -DHAVE_MPI=1 -UHAVE_OPENCL
clamr_mpionly_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la timer/libtimer.la l7/libl7.la genmalloc/libgenmalloc.la -lm ${MPE_LIBS} ${MPI_LIBS}

clamr_mpionly_SOURCES = ${clamr_mpionly_CXX_SRCS} ${clamr_mpionly_C_SRCS} ${clamr_mpionly_H_SRCS}

clamr_mpicheck_CXX_SRCS = clamr_mpicheck.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_mpicheck_C_SRCS = display.c
clamr_mpicheck_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h
clamr_mpicheck_CPPFLAGS = -DHAVE_MPI=1 -UHAVE_OPENCL
clamr_mpicheck_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la timer/libtimer.la l7/libl7.la genmalloc/libgenmalloc.la -lm ${MPE_LIBS} ${MPI_LIBS}

clamr_mpicheck_SOURCES = ${clamr_mpicheck_CXX_SRCS} ${clamr_mpicheck_C_SRCS} ${clamr_mpicheck_H_SRCS}

clamr_gpucheck_CXX_SRCS = clamr_gpucheck.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_gpucheck_C_SRCS = display.c reduce.c
clamr_gpucheck_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h
clamr_gpucheck_CPPFLAGS = -UHAVE_MPI -DHAVE_OPENCL=1
clamr_gpucheck_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la ezcl/libezcl.la timer/libtimer.la genmalloc/libgenmalloc.la -lm ${MPE_NOMPI_LIBS} ${OPENCL_LIBS}

clamr_gpucheck_SOURCES = ${clamr_gpucheck_CXX_SRCS} ${clamr_gpucheck_C_SRCS} ${clamr_gpucheck_H_SRCS}

clamr_checkall_CXX_SRCS = clamr_checkall.cpp mesh.cpp partition.cpp state.cpp input.cpp
clamr_checkall_C_SRCS = display.c reduce.c 
clamr_checkall_H_SRCS = display.h mesh.h partition.h reorder.h state.h input.h reduce.h
clamr_checkall_CPPFLAGS = -DHAVE_MPI=1 -DHAVE_OPENCL=1
clamr_checkall_LDADD = hsfc/libhsfc.la kdtree/libkdtree.la zorder/libzorder.la s7/libs7.la ezcl/libpezcl.la timer/libtimer.la l7/libl7.la genmalloc/libgenmalloc.la -lm ${MPE_LIBS} ${MPI_LIBS} ${OPENCL_LIBS}

clamr_checkall_SOURCES = ${clamr_checkall_CXX_SRCS} ${clamr_checkall_C_SRCS} ${clamr_checkall_H_SRCS}

EXTRA_DIST = wave_kern.cl wave_kern_calc.cl reduce.cl config.h clamr.xcodeproj/project.pbxproj clamr.1 fileamr.in AUTHORS ChangeLog runtests.sh newbatch.sh LICENSE

l7/libl7.la:
	cd l7; make

hsfc/libhsfc.la:
	cd hsfc; make

kdtree/libkdtree.la:
	cd kdtree; make

zorder/libzorder.la:
	cd zorder; make

s7/libs7.la:
	cd s7; make

ezcl/libezcl.la:
	cd ezcl; make

timer/libtimer.la:
	cd timer; make

genmalloc/libgenmalloc.la:
	cd genmalloc; make

distclean-local:
	rm -f config.h.in~
	rm -rf autom4te.cache
	rm -rf config.log config.status
