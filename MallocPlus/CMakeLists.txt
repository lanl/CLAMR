########### global settings ###############
set(H_SRCS MallocPlus.h)
set(CXX_SRCS MallocPlus.cpp)


########### MallocPlus target ###############
set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(MallocPlus SHARED ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS MallocPlus DESTINATION lib)

########### dMallocPlus target ###############
set(dMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(dMallocPlus SHARED ${dMallocPlus_LIB_SRCS})

set_target_properties(dMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(dMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)
target_link_libraries(dMallocPlus ${OPENCL_LIBRARIES})
target_link_libraries(dMallocPlus ezcl)
install(TARGETS dMallocPlus DESTINATION lib)

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
