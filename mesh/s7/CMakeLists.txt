########### global settings ###############
set(H_SRCS s7.h)
set(C_SRCS s7.c)

set(s7_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### s7 target ###############
add_library(s7 SHARED ${s7_LIB_SRCS})

set_target_properties(s7 PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS s7 DESTINATION lib)

########### install files ###############

#install(FILES  s7.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libs7.la
#all: libs7.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libs7.la
#include_HEADERS = s7.h
#
#H_SRCS = s7.h
#
#C_SRCS = s7.c
#
#libs7_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
