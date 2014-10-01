########### global settings ###############
set(H_SRCS hsfc.h)

set(C_SRCS hsfc2part.c hsfc.c hsfcsort.c)

set(hsfc_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### hsfc target ###############

add_library(hsfc SHARED ${hsfc_LIB_SRCS})

set_target_properties(hsfc PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS hsfc DESTINATION lib)

########### install files ###############

#install(FILES  hsfc.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libhsfc.la
#all: libhsfc.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libhsfc.la
#include_HEADERS = hsfc.h
#
#H_SRCS = hsfc.h
#
#C_SRCS = hsfc2part.c hsfc.c hsfcsort.c
#
#libhsfc_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
