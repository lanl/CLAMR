default: libs7.la
all: libs7.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libs7.la
include_HEADERS = s7.h

H_SRCS = s7.h

C_SRCS = s7.c

libs7_la_SOURCES = ${C_SRCS} ${H_SRCS}

