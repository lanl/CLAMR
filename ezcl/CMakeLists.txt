cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (EZCL)
# The version number.
set (EZCL_VERSION_MAJOR 2)
set (EZCL_VERSION_MINOR 0.7)
set (EZCL_VERSION "${EZCL_VERSION_MAJOR}.${EZCL_VERSION_MINOR}")
 
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
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")

include_directories(${PROJECT_SOURCE_DIR})

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
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${OPENCL_INCLUDE_DIRS}")

set(H_SRCS ezcl.h)
set(C_SRCS ezcl.c)

########### ezcl target ###############
set(ezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

add_library(ezcl SHARED ${ezcl_LIB_SRCS})

set_target_properties(ezcl PROPERTIES VERSION ${EZCL_VERSION} SOVERSION ${EZCL_VERSION_MAJOR})
target_link_libraries(ezcl ${OPENCL_LIBRARIES})
install(TARGETS ezcl DESTINATION lib)

########### pezcl target ###############

set(pezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

add_library(pezcl SHARED ${pezcl_LIB_SRCS})

set_target_properties(pezcl PROPERTIES VERSION ${EZCL_VERSION} SOVERSION ${EZCL_VERSION_MAJOR})
set_target_properties(pezcl PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
target_link_libraries(pezcl ${OPENCL_LIBRARIES})
target_link_libraries(pezcl ${MPI_LIBRARIES})
install(TARGETS pezcl DESTINATION lib)

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
