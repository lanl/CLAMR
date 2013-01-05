########### global settings ###############
set(H_SRCS zorder.h)
set(C_SRCS zorder.c)

set(zorder_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### zorder target ###############
add_library(zorder SHARED ${zorder_LIB_SRCS})

set_target_properties(zorder PROPERTIES VERSION 2.0.0 SOVERSION 2)
target_link_libraries(zorder s7)
install(TARGETS zorder DESTINATION lib)

########### install files ###############

#install(FILES  zorder.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libzorder.la
#all: libzorder.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libzorder.la
#include_HEADERS = zorder.h
#
#H_SRCS = zorder.h
#
#C_SRCS = zorder.c
#
#libzorder_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
