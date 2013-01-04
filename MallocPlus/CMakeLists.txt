########### global settings ###############
set(H_SRCS MallocPlus.h)
set(CXX_SRCS MallocPlus.cpp)

set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

########### MallocPlus target ###############
add_library(MallocPlus SHARED ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS MallocPlus DESTINATION lib)

########### install files ###############

#install(FILES  MallocPlus.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libMallocPlus.la
#all: libMallocPlus.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libMallocPlus.la
#include_HEADERS = MallocPlus.h
#
#H_SRCS = MallocPlus.h
#
#CXX_SRCS = MallocPlus.cpp
#
#libMallocPlus_la_SOURCES = ${CXX_SRCS} ${H_SRCS}
#
