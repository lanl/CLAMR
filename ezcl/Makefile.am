default: libezcl.a
all: libezcl.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libezcl.a
include_HEADERS = ezcl.h

H_SRCS = ezcl.h

C_SRCS = ezcl.c

libezcl_a_SOURCES = ${C_SRCS} ${H_SRCS}

