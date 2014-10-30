cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (GENMALLOC)
# The version number.
set (GENMALLOC_VERSION_MAJOR 2)
set (GENMALLOC_VERSION_MINOR 0.7)
set (GENMALLOC_VERSION "${GENMALLOC_VERSION_MAJOR}.${GENMALLOC_VERSION_MINOR}")

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")

########### global settings ###############
set(H_SRCS genmalloc.h)
set(C_SRCS genmalloc.c)

set(genmalloc_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### genmalloc target ###############

add_library(genmalloc SHARED ${genmalloc_LIB_SRCS})

set_target_properties(genmalloc PROPERTIES VERSION ${GENMALLOC_VERSION} SOVERSION ${GENMALLOC_VERSION_MAJOR})
install(TARGETS genmalloc DESTINATION lib)

########### install files ###############

install(FILES  genmalloc.h DESTINATION include/genmalloc)
