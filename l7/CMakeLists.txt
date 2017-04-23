cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (L7)
# The version number.
set (L7_VERSION_MAJOR 2)
set (L7_VERSION_MINOR 0.7)
set (L7_VERSION "${L7_VERSION_MAJOR}.${L7_VERSION_MINOR}")

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -vec-report3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")

if (NOT DEFINED OPENCL_FOUND AND NOT MIC_NATIVE)
   find_package(OpenCL)
   #  OPENCL_FOUND        : TRUE if found, FALSE otherwise
   #  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
   #  OPENCL_LIBRARIES    : The libraries to link against
   if (OPENCL_FOUND)
      include_directories(${OPENCL_INCLUDE_DIRS})
   endif (OPENCL_FOUND)
endif (NOT DEFINED OPENCL_FOUND AND NOT MIC_NATIVE)

if (NOT DEFINED MPI_FOUND)
   find_package(MPI)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if (NOT DEFINED OPENMP_FOUND)
   find_package( OpenMP)
endif (NOT DEFINED OPENMP_FOUND)

enable_testing()

add_subdirectory(tests)

########### global settings ###############
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI=1")

set(H_SRCS l7.h l7_assert.h l7p.h)

set(C_SRCS
      l7p_sizeof.c      l7_terminate.c      l7_init.c        l7_update.c 
      l7_free.c         l7p_set_database.c  l7_setup.c       l7_push_setup.c
      l7_push_update.c  l7_push_free.c      l7_dev_update.c  l7_dev_setup.c
      l7_dev_free.c     l7_utils.c          l7_reduction.c   l7_broadcast.c
      l7p_mpi_type.c
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/l7_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl l7_kern.cl > l7_kernel.inc
                  DEPENDS l7_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
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
if (MPI_FOUND AND OPENCL_FOUND)
   set(dl7_LIB_SRCS ${C_SRCS} ${H_SRCS})

   add_library(dl7 STATIC ${l7_LIB_SRCS})

   set_target_properties(dl7 PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dl7 PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(dl7 ${MPI_LIBRARIES})
   target_link_libraries(dl7 ${OPENCL_LIBRARIES})
   add_dependencies(dl7 l7_kernel_source)
   install(TARGETS dl7 DESTINATION lib)
endif(MPI_FOUND AND OPENCL_FOUND)

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

