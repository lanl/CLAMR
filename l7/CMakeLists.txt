cmake_minimum_required(VERSION 2.6)

project (L7 C)
# The version number.
set (L7_VERSION_MAJOR 2)
set (L7_VERSION_MINOR 0.7)
set (L7_VERSION "${L7_VERSION_MAJOR}.${L7_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

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

if (NOT DEFINED VECTOR_C_FLAGS)
   find_package(Vector)
   if (CMAKE_VECTOR_VERBOSE)
       set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
   endif (CMAKE_VECTOR_VERBOSE)
   if (CMAKE_VECTOR_NOVEC)
       set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
   endif (CMAKE_VECTOR_NOVEC)
endif()

set(H_SRCS l7.h l7_assert.h l7p.h)

set(C_SRCS
      l7p_sizeof.c      l7_terminate.c      l7_init.c        l7_update.c 
      l7_free.c         l7p_set_database.c  l7_setup.c       l7_push_setup.c
      l7_push_update.c  l7_push_free.c      l7_dev_update.c  l7_dev_setup.c
      l7_dev_free.c     l7_utils.c          l7_reduction.c   l7_broadcast.c
      l7p_mpi_type.c
)

set_source_files_properties(${C_SRCS} PROPERTIES COMPILE_FLAGS "${VECTOR_C_FLAGS}")

if (NOT DEFINED OpenCL_FOUND)
   if (NOT APPLE OR DARWIN_VERSION LESS 18)
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
   endif (NOT APPLE OR DARWIN_VERSION LESS 18)
endif (NOT DEFINED OpenCL_FOUND)

if (NOT DEFINED MPI_FOUND)
   find_package(MPI)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI=1")
endif (MPI_FOUND)

if (NOT DEFINED OPENMP_FOUND)
   find_package( OpenMP)
endif (NOT DEFINED OPENMP_FOUND)

enable_testing()

add_subdirectory(tests)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/l7_kernel.inc
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl ${CMAKE_CURRENT_SOURCE_DIR}/l7_kern.cl > ${CMAKE_CURRENT_BINARY_DIR}/l7_kernel.inc
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/l7_kern.cl ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl)
add_custom_target(l7_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/l7_kernel.inc)

########### l7 target ###############
if(MPI_FOUND)
   set(l7_LIB_SRCS ${C_SRCS} ${H_SRCS})

   add_library(l7 STATIC ${l7_LIB_SRCS})

   set_target_properties(l7 PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(l7 PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(l7 ${MPI_LIBRARIES})
   install(TARGETS l7 DESTINATION lib)
endif (MPI_FOUND)

########### mpl7 target ###############
if(MPI_FOUND)
   set(mpl7_LIB_SRCS ${C_SRCS} ${H_SRCS})

   add_library(mpl7 STATIC ${l7_LIB_SRCS})

   set_target_properties(mpl7 PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(mpl7 PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   if (OPENMP_FOUND)
      set_target_properties(mpl7 PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
      set_target_properties(mpl7 PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
      set_target_properties(mpl7 PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
   endif (OPENMP_FOUND)

   target_link_libraries(mpl7 ${MPI_LIBRARIES})
   install(TARGETS mpl7 DESTINATION lib)
endif (MPI_FOUND)

########### dl7 target ###############
if (MPI_FOUND AND OpenCL_FOUND)
   set(dl7_LIB_SRCS ${C_SRCS} ${H_SRCS})

   add_library(dl7 STATIC ${l7_LIB_SRCS})

   set_target_properties(dl7 PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dl7 PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(dl7 ${MPI_LIBRARIES})
   target_link_libraries(dl7 ${OpenCL_LIBRARIES})
   add_dependencies(dl7 l7_kernel_source)
   install(TARGETS dl7 DESTINATION lib)
endif(MPI_FOUND AND OpenCL_FOUND)

########### install files ###############

install(FILES  l7.h DESTINATION include/l7)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/L7.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/L7.pdf)
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
   add_custom_target(L7_doc
                     COMMAND make L7_pdfdoc
                     COMMAND make L7_htmldoc
                     COMMAND make L7_singlehtmldoc)
   set_target_properties(L7_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/L7.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(L7_pdfdoc DEPENDS _build/latex/L7.pdf)
   set_target_properties(L7_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(L7_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(L7_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(L7_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(L7_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "l7_kernel.inc;install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")

