default: libgenmalloc.la
all: libgenmalloc.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libgenmalloc.la
include_HEADERS = genmalloc.h

H_SRCS = genmalloc.h

C_SRCS = genmalloc.c

libgenmalloc_la_SOURCES = ${C_SRCS} ${H_SRCS}

