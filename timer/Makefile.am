default: libtimer.a
all: libtimer.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libtimer.a
include_HEADERS = timer.h

H_SRCS = timer.h

C_SRCS = timer.c

libtimer_a_SOURCES = ${C_SRCS} ${H_SRCS}

