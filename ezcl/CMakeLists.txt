########### global settings ###############
set(CMAKE_C_FLAGS "-std=c99 -I${OPENCL_INCLUDE_DIRS}")

set(H_SRCS ezcl.h)
set(C_SRCS ezcl.c)

########### ezcl target ###############
set(ezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

add_library(ezcl SHARED ${ezcl_LIB_SRCS})

set_target_properties(ezcl PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS ezcl DESTINATION lib)

########### pezcl target ###############

set(pezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

add_library(pezcl SHARED ${pezcl_LIB_SRCS})

set_target_properties(pezcl PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(pezcl PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
install(TARGETS pezcl DESTINATION lib)

########### install files ###############

#install(FILES  ezcl.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libezcl.la libpezcl.la
#all: libezcl.la libpezcl.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libezcl.la libpezcl.la
#include_HEADERS = ezcl.h
#
#H_SRCS = ezcl.h
#
#C_SRCS = ezcl.c
#
#libezcl_la_CPPFLAGS = -UHAVE_MPI
#libezcl_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
#libpezcl_la_CPPFLAGS = -DHAVE_MPI=1
#libpezcl_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
