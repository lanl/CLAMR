default: libs7.a
all: libs7.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libs7.a
include_HEADERS = s7.h

H_SRCS = s7.h

C_SRCS = s7.c

libs7_a_SOURCES = ${C_SRCS} ${H_SRCS}

