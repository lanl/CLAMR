########### global settings ###############
set(CMAKE_C_FLAGS "-std=c99 -DHAVE_MPI=1")

set(H_SRCS l7.h l7_assert.h l7p.h)

set(C_SRCS
      l7p_sizeof.c  l7_terminate.c      l7_init.c    l7_update.c 
      l7_free.c     l7p_set_database.c  l7_setup.c
)

set(l7_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### l7 target ###############

add_library(l7 SHARED ${l7_LIB_SRCS})

set_target_properties(l7 PROPERTIES VERSION 2.0.0 SOVERSION 2)
target_link_libraries(l7 ${MPI_LIBRARIES})
install(TARGETS l7 DESTINATION lib)

########### install files ###############

#install(FILES  l7.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libl7.la
#all: libl7.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libl7.la
#include_HEADERS = l7.h
#
#
#libl7_la_SOURCES = ${C_SRCS} ${H_SRCS}
#libl7_la_CPPFLAGS = -DHAVE_MPI=1
#
