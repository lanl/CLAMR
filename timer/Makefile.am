default: libtimer.la
all: libtimer.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libtimer.la
include_HEADERS = timer.h

H_SRCS = timer.h

C_SRCS = timer.c

libtimer_la_SOURCES = ${C_SRCS} ${H_SRCS}

