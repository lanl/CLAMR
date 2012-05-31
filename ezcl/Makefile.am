default: libezcl.la
all: libezcl.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libezcl.la
include_HEADERS = ezcl.h

H_SRCS = ezcl.h

C_SRCS = ezcl.c

libezcl_la_SOURCES = ${C_SRCS} ${H_SRCS}

