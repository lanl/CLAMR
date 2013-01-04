########### global settings ###############
set(H_SRCS timer.h)
set(C_SRCS timer.c)

set(timer_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### timer target ###############
add_library(timer SHARED ${timer_LIB_SRCS})

set_target_properties(timer PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS timer DESTINATION lib)

########### install files ###############

#install(FILES  timer.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libtimer.la
#all: libtimer.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libtimer.la
#include_HEADERS = timer.h
#
#H_SRCS = timer.h
#
#C_SRCS = timer.c
#
#libtimer_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
