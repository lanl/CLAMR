default: libzorder.la
all: libzorder.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libzorder.la
include_HEADERS = zorder.h

H_SRCS = zorder.h

C_SRCS = zorder.c

libzorder_la_SOURCES = ${C_SRCS} ${H_SRCS}

