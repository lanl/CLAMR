########### global settings ###############
set(H_SRCS MallocPlus.h)
set(CXX_SRCS MallocPlus.cpp)


########### MallocPlus target ###############
set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(MallocPlus STATIC ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS MallocPlus DESTINATION lib)

########### MallocPlus target ###############
set(pMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(pMallocPlus STATIC ${MallocPlus_LIB_SRCS})

set_target_properties(pMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(pMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
target_link_libraries(pMallocPlus ${MPI_LIBRARIES})
install(TARGETS pMallocPlus DESTINATION lib)

########### dMallocPlus target ###############
set(dMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(dMallocPlus STATIC ${dMallocPlus_LIB_SRCS})

set_target_properties(dMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(dMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)
target_link_libraries(dMallocPlus ${OPENCL_LIBRARIES})
target_link_libraries(dMallocPlus ezcl)
#if (Boost_FOUND AND MPI_FOUND)
#   message(STATUS "Adding j7 support ot MallocPlus")
#   target_link_libraries(dMallocPlus j7)
#endif()
install(TARGETS dMallocPlus DESTINATION lib)

########### dpMallocPlus target ###############
set(dpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(dpMallocPlus STATIC ${dMallocPlus_LIB_SRCS})

set_target_properties(dpMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(dpMallocPlus PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
target_link_libraries(dpMallocPlus ${MPI_LIBRARIES})
target_link_libraries(dpMallocPlus ${OPENCL_LIBRARIES})
target_link_libraries(dpMallocPlus ezcl)
#if (Boost_FOUND AND MPI_FOUND)
#   message(STATUS "Adding j7 support ot MallocPlus")
#   target_link_libraries(dMallocPlus j7)
#endif()
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
