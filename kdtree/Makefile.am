default: libkdtree.la
all: libkdtree.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libkdtree.la
include_HEADERS = KDTree.h

H_SRCS = Bounds.h Globals.h KDTree.h

C_SRCS = Bounds.c KDTree.c

libkdtree_la_SOURCES = ${C_SRCS} ${H_SRCS}

