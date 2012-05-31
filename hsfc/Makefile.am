default: libhsfc.la
all: libhsfc.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libhsfc.la
include_HEADERS = hsfc.h

H_SRCS = hsfc.h

C_SRCS = hsfc2part.c hsfc.c hsfcsort.c

libhsfc_la_SOURCES = ${C_SRCS} ${H_SRCS}

