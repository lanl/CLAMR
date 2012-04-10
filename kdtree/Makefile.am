default: libkdtree.a
all: libkdtree.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libkdtree.a
include_HEADERS = KDTree.h

H_SRCS = Bounds.h Globals.h KDTree.h

C_SRCS = Bounds.c KDTree.c

libkdtree_a_SOURCES = ${C_SRCS} ${H_SRCS}

