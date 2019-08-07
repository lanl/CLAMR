cmake_minimum_required(VERSION 3.1.3)

project (S7)
# The version number.
set (S7_VERSION_MAJOR 2)
set (S7_VERSION_MINOR 0.7)
set (S7_VERSION "${S7_VERSION_MAJOR}.${S7_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

if (NOT CMAKE_BUILD_TYPE)
   #set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

set (CMAKE_C_STANDARD 99)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -fstrict-aliasing -ftree-vectorize")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fopt-info-vec-optimized")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.3")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fopenmp-simd")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.3")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -mtune=native")
   elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mcpu=powerpc64le")
   else ()
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
   endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -ansi-alias -vecabi=cmdtarget -xHost -qopenmp-simd -qopt-report=5 -qopt-report-phase=vec,loop")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_C_FLAGS_RELEASE "-qopt-zmm-usage=high ${CMAKE_C_FLAGS_RELEASE}")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
endif()

########### global settings ###############
set(H_SRCS s7.h)
set(C_SRCS s7.c)

set(s7_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### s7 target ###############
add_library(s7 SHARED ${s7_LIB_SRCS})

set_target_properties(s7 PROPERTIES VERSION ${S7_VERSION} SOVERSION ${S7_VERSION_MAJOR})
install(TARGETS s7 DESTINATION lib)

########### install files ###############

install(FILES  s7.h DESTINATION include/s7)

