cmake_minimum_required(VERSION 3.1.3)

project (EZCL C)
# The version number.
set (EZCL_VERSION_MAJOR 2)
set (EZCL_VERSION_MINOR 0.7)
set (EZCL_VERSION "${EZCL_VERSION_MAJOR}.${EZCL_VERSION_MINOR}")
 
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
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

set (CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO})

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
    set( CMAKE_VECTOR_NOVEC TRUE)
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")


if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
endif()

if(NOT KOKKOS_ENABLE)
   if (NOT DEFINED VECTOR_C_FLAGS)
      find_package(Vector)
      if (CMAKE_VECTOR_VERBOSE)
          set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
      endif (CMAKE_VECTOR_VERBOSE)
      if (CMAKE_VECTOR_NOVEC)
          set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
      endif (CMAKE_VECTOR_NOVEC)
  endif(NOT DEFINED VECTOR_C_FLAGS)
endif(NOT KOKKOS_ENABLE)

set_source_files_properties(ezcl.c PROPERTIES COMPILE_FLAGS "${VECTOR_C_FLAGS}")

include_directories(${PROJECT_SOURCE_DIR})

if (NOT DEFINED OpenCL_FOUND)
   find_package(OpenCL)
   #  OpenCL_FOUND          - True if OpenCL was found
   #  OpenCL_INCLUDE_DIRS   - include directories for OpenCL
   #  OpenCL_LIBRARIES      - link against this library to use OpenCL
   #  OpenCL_VERSION_STRING - Highest supported OpenCL version (eg. 1.2)
   #  OpenCL_VERSION_MAJOR  - The major version of the OpenCL implementation
   #  OpenCL_VERSION_MINOR  - The minor version of the OpenCL implementation
   if (OpenCL_FOUND)
      set(HAVE_CL_DOUBLE ON CACHE BOOL "Have OpenCL Double")
      set(NO_CL_DOUBLE OFF)
      include_directories(${OpenCL_INCLUDE_DIRS})
   endif (OpenCL_FOUND)
endif (NOT DEFINED OpenCL_FOUND)
 
if (NOT DEFINED MPI_FOUND)
   find_package(MPI REQUIRED)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if (NOT DEFINED OPENMP_FOUND)
   find_package( OpenMP)
endif (NOT DEFINED OPENMP_FOUND)

if (NOT PRECISION_TYPE)
   set (PRECISION_TYPE full_precision)
endif (NOT PRECISION_TYPE)

if (PRECISION_TYPE STREQUAL "mixed_precision")
   set (MIXED_PRECISION on)
   set (PRECISION_TYPE mixed_precision)
   add_definitions(-DMIXED_PRECISION)
endif (PRECISION_TYPE STREQUAL "mixed_precision")

if (PRECISION_TYPE STREQUAL "full_precision")
   set (FULL_PRECISION on)
   set (PRECISION_TYPE full_precision)
   add_definitions(-DFULL_PRECISION)
endif (PRECISION_TYPE STREQUAL "full_precision")

if (PRECISION_TYPE STREQUAL "minimum_precision")
   set (MINIMUM_PRECISION on)
   set (PRECISION_TYPE minimum_precision)
   add_definitions(-DMINIMUM_PRECISION)
endif (PRECISION_TYPE STREQUAL "minimum_precision")

enable_testing()

add_subdirectory(tests)

########### global settings ###############
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${OpenCL_INCLUDE_DIRS}")

set(H_SRCS ezcl.h)
set(C_SRCS ezcl.c)

########### ezcl target ###############
set(ezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

add_library(ezcl SHARED ${ezcl_LIB_SRCS})

set_target_properties(ezcl PROPERTIES VERSION ${EZCL_VERSION} SOVERSION ${EZCL_VERSION_MAJOR})
target_link_libraries(ezcl ${OpenCL_LIBRARIES})
install(TARGETS ezcl DESTINATION lib)

########### pezcl target ###############

if (MPI_FOUND)
   set(pezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

   add_library(pezcl SHARED ${pezcl_LIB_SRCS})

   set_target_properties(pezcl PROPERTIES VERSION ${EZCL_VERSION} SOVERSION ${EZCL_VERSION_MAJOR})
   set_target_properties(pezcl PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(pezcl ${OpenCL_LIBRARIES})
   target_link_libraries(pezcl ${MPI_LIBRARIES})
   install(TARGETS pezcl DESTINATION lib)
endif (MPI_FOUND)

########### install files ###############

install(FILES  ezcl.h DESTINATION include/ezcl)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND pdfdoc
                     COMMAND htmldoc
                     COMMAND singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/ezcl.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/ezcl.pdf)
   set_target_properties(pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(ezcl_doc
                     COMMAND make ezcl_pdfdoc
                     COMMAND make ezcl_htmldoc
                     COMMAND make ezcl_singlehtmldoc)
   set_target_properties(ezcl_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/ezcl.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(ezcl_pdfdoc DEPENDS _build/latex/ezcl.pdf)
   set_target_properties(ezcl_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(ezcl_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(ezcl_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(ezcl_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(ezcl_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")
