default: libl7.a
all: libl7.a

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LIBRARIES = libl7.a
include_HEADERS = l7.h

H_SRCS = l7.h l7_assert.h l7p.h

C_SRCS = l7p_sizeof.c  l7_terminate.c      l7_init.c    l7_update.c \
         l7_free.c     l7p_set_database.c  l7_setup.c

libl7_a_SOURCES = ${C_SRCS} ${H_SRCS}

