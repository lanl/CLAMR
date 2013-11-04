########### global settings ###############
set(H_SRCS MallocPlus.h)
set(CXX_SRCS MallocPlus.cpp)


########### MallocPlus target ###############
set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(MallocPlus STATIC ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS MallocPlus DESTINATION lib)

########### dMallocPlus target ###############
set(dMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(dMallocPlus STATIC ${dMallocPlus_LIB_SRCS})

set_target_properties(dMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(dMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)
target_link_libraries(dMallocPlus ${OPENCL_LIBRARIES})
target_link_libraries(dMallocPlus ezcl)
if (Boost_FOUND AND MPI_FOUND)
    message(STATUS "Adding j7 support ot MallocPlus")
    target_link_libraries(dMallocPlus j7)
    # XXX want to set HAVE_J7 at the top-level not all over the place
    add_definitions(-DHAVE_J7)
endif()
install(TARGETS dMallocPlus DESTINATION lib)

########### install files ###############

#install(FILES  MallocPlus.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libMallocPlus.la
#all: libMallocPlus.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libMallocPlus.la
#include_HEADERS = MallocPlus.h
#
#H_SRCS = MallocPlus.h
#
#CXX_SRCS = MallocPlus.cpp
#
#libMallocPlus_la_SOURCES = ${CXX_SRCS} ${H_SRCS}
#
