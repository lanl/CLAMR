default: libMallocPlus.la
all: libMallocPlus.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libMallocPlus.la
include_HEADERS = MallocPlus.h

H_SRCS = MallocPlus.h

CXX_SRCS = MallocPlus.cpp

libMallocPlus_la_SOURCES = ${CXX_SRCS} ${H_SRCS}

