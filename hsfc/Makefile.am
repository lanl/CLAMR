default: libhsfc.a
all: libhsfc.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libhsfc.a
include_HEADERS = hsfc.h

H_SRCS = hsfc.h

C_SRCS = hsfc2part.c hsfc.c hsfcsort.c

libhsfc_a_SOURCES = ${C_SRCS} ${H_SRCS}

