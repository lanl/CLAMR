default: libezcl.la libpezcl.la
all: libezcl.la libpezcl.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libezcl.la libpezcl.la
include_HEADERS = ezcl.h

H_SRCS = ezcl.h

C_SRCS = ezcl.c

libezcl_la_CPPFLAGS = -UHAVE_MPI
libezcl_la_SOURCES = ${C_SRCS} ${H_SRCS}

libpezcl_la_CPPFLAGS = -DHAVE_MPI=1
libpezcl_la_SOURCES = ${C_SRCS} ${H_SRCS}

