########### global settings ###############
set(CMAKE_C_FLAGS "-std=c99")
set(H_SRCS genmalloc.h)
set(C_SRCS genmalloc.c)

set(genmalloc_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### genmalloc target ###############

add_library(genmalloc SHARED ${genmalloc_LIB_SRCS})

set_target_properties(genmalloc PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS genmalloc DESTINATION lib)

########### install files ###############

#install(FILES  genmalloc.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libgenmalloc.la
#all: libgenmalloc.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libgenmalloc.la
#include_HEADERS = genmalloc.h
#
#H_SRCS = genmalloc.h
#
#C_SRCS = genmalloc.c
#
#libgenmalloc_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
