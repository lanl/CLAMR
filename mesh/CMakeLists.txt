########### global settings ###############
#if (CMAKE_C_COMPILER MATCHES "gcc$")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorizer-verbose=2")
#endif (CMAKE_C_COMPILER MATCHES "gcc$")
#if (CMAKE_C_COMPILER MATCHES "icc$")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -vec-report3")
#endif (CMAKE_C_COMPILER MATCHES "icc$")

set(H_SRCS mesh.h partition.h reorder.h reduce.h)

set(CXX_SRCS mesh.cpp partition.cpp)

set(C_SRCS reduce.c)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mesh/mesh_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/mesh/mesh_kern.cl > mesh_kernel.inc
                  DEPENDS mesh_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(mesh_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mesh/mesh_kernel.inc)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mesh/reduce_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/mesh/reduce.cl > reduce_kernel.inc
                  DEPENDS reduce.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(reduce_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mesh/reduce_kernel.inc)

########### mesh target ###############
set(mesh_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${H_SRCS})

add_library(mesh STATIC ${mesh_LIB_SRCS})

set_target_properties(mesh PROPERTIES VERSION 2.0.0 SOVERSION 2)
target_link_libraries(mesh)
add_dependencies(mesh reduce_kernel_source)
install(TARGETS mesh DESTINATION lib)

########### pmesh target ###############
if (MPI_FOUND)
   set(pmesh_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${H_SRCS})

   add_library(pmesh STATIC ${pmesh_LIB_SRCS})

   set_target_properties(pmesh PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(pmesh PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(pmesh ${MPI_LIBRARIES})
   install(TARGETS pmesh DESTINATION lib)
endif (MPI_FOUND)

########### dmesh target ###############
if (OPENCL_FOUND)
   set(dmesh_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${H_SRCS})

   add_library(dmesh STATIC ${dmesh_LIB_SRCS})

   set_target_properties(dmesh PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dmesh PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)
   target_link_libraries(dmesh ${OPENCL_LIBRARIES})
   add_dependencies(dmesh mesh_kernel_source)
   add_dependencies(dmesh reduce_kernel_source)
   install(TARGETS dmesh DESTINATION lib)
endif(OPENCL_FOUND)

########### dpmesh target ###############
if (OPENCL_FOUND)
   set(dpmesh_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${H_SRCS})

   add_library(dpmesh STATIC ${dpmesh_LIB_SRCS})

   set_target_properties(dpmesh PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dpmesh PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(dpmesh ${MPI_LIBRARIES})
   target_link_libraries(dpmesh ${OPENCL_LIBRARIES})
   add_dependencies(dpmesh mesh_kernel_source)
   add_dependencies(dpmesh reduce_kernel_source)
   install(TARGETS dpmesh DESTINATION lib)
endif (OPENCL_FOUND)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "mesh_kernel.inc;reduce_kernel.inc")

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
