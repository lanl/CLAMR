default: libzorder.a
all: libzorder.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libzorder.a
include_HEADERS = zorder.h

H_SRCS = zorder.h

C_SRCS = zorder.c

libzorder_a_SOURCES = ${C_SRCS} ${H_SRCS}

