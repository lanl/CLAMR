cmake_minimum_required(VERSION 2.6)

project (MALLOCPLUS CXX Fortran)
# The version number.
set (MALLOCPLUS_VERSION_MAJOR 2)
set (MALLOCPLUS_VERSION_MINOR 0.7)
set (MALLOCPLUS_VERSION "${MALLOC_PLUS_VERSION_MAJOR}.${MALLOC_PLUS_VERSION_MINOR}")

if( DEFINED ENV{PE_ENV})
  cmake_minimum_required (VERSION 3.5)
  set( CMAKE_SYSTEM_NAME CrayLinuxEnvironment CACHE STRING "description")
  set( MPIEXEC "aprun" CACHE STRING "description")
else ()
  cmake_minimum_required (VERSION 2.8)
endif()

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

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Cray")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
endif()

if (NOT DEFINED VECTOR_CXX_FLAGS)
   find_package(Vector)
   if (CMAKE_VECTOR_VERBOSE)
       set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
   endif (CMAKE_VECTOR_VERBOSE)
   if (CMAKE_VECTOR_NOVEC)
       set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
   endif (CMAKE_VECTOR_NOVEC)
endif()

set_source_files_properties(MallocPlus.cpp PROPERTIES COMPILE_FLAGS "${VECTOR_CXX_FLAGS}")

if (NOT DEFINED OpenCL_FOUND)
   find_package(OpenCL)
   #  OpenCL_FOUND          - True if OpenCL was found
   #  OpenCL_INCLUDE_DIRS   - include directories for OpenCL
   #  OpenCL_LIBRARIES      - link against this library to use OpenCL
   #  OpenCL_VERSION_STRING - Highest supported OpenCL version (eg. 1.2)
   #  OpenCL_VERSION_MAJOR  - The major version of the OpenCL implementation
   if (OpenCL_FOUND)
      set(HAVE_CL_DOUBLE ON CACHE BOOL "Have OpenCL Double")
      set(NO_CL_DOUBLE OFF)
      include_directories(${OpenCL_INCLUDE_DIRS})
   endif (OpenCL_FOUND)
endif (NOT DEFINED OpenCL_FOUND)

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
#set(H_SRCS CMallocPlus.h MallocPlus.h)
set(H_SRCS MallocPlus.h)
#set(C_SRCS CMallocPlus.c)
set(CXX_SRCS MallocPlus.cpp)
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   set(F90_SRCS FMallocPlus.F90)
endif (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

########### MallocPlus target ###############
set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${F90_SRCS} ${H_SRCS})

add_library(MallocPlus STATIC ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/MallocPlus.dir)
set_target_properties(MallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
#install(TARGETS MallocPlus DESTINATION lib)

########### pMallocPlus target ###############
if (MPI_FOUND)
   set(pMallocPlus_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${F90_SRCS} ${H_SRCS})

   add_library(pMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(pMallocPlus PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/pMallocPlus.dir)
   set_target_properties(pMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(pMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(pMallocPlus ${MPI_LIBRARIES})
#   install(TARGETS pMallocPlus DESTINATION lib)
endif (MPI_FOUND)

########### mpMallocPlus target ###############
if (OPENMP_FOUND)
   set(mpMallocPlus_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${F90_SRCS} ${H_SRCS})

   add_library(mpMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(mpMallocPlus PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/mpMallocPlus.dir)
   set_target_properties(mpMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(mpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(mpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(mpMallocPlus PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

#   install(TARGETS mpMallocPlus DESTINATION lib)
endif (OPENMP_FOUND)

########### pmpMallocPlus target ###############
if (OPENMP_FOUND AND MPI_FOUND)
   set(pmpMallocPlus_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${F90_SRCS} ${H_SRCS})

   add_library(pmpMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(pmpMallocPlus PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/pmpMallocPlus.dir)
   set_target_properties(pmpMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(pmpMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   set_target_properties(pmpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(pmpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(pmpMallocPlus PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

   target_link_libraries(pmpMallocPlus ${MPI_LIBRARIES})
#   install(TARGETS pmpMallocPlus DESTINATION lib)
endif (OPENMP_FOUND AND MPI_FOUND)

########### dMallocPlus target ###############
if (OpenCL_FOUND)
   set(dMallocPlus_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${F90_SRCS} ${H_SRCS})

   add_library(dMallocPlus STATIC ${dMallocPlus_LIB_SRCS})

   set_target_properties(dMallocPlus PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/dMallocPlus.dir)
   set_target_properties(dMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(dMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)
   target_link_libraries(dMallocPlus ${OpenCL_LIBRARIES})
   target_link_libraries(dMallocPlus ezcl)
   #if (Boost_FOUND AND MPI_FOUND)
   #   message(STATUS "Adding j7 support ot MallocPlus")
   #   target_link_libraries(dMallocPlus j7)
   #endif()
#   install(TARGETS dMallocPlus DESTINATION lib)
endif(OpenCL_FOUND)

########### dpMallocPlus target ###############
if (OpenCL_FOUND AND MPI_FOUND)
   set(dpMallocPlus_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${F90_SRCS} ${H_SRCS})

   add_library(dpMallocPlus STATIC ${dpMallocPlus_LIB_SRCS})

   set_target_properties(dpMallocPlus PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/dpMallocPlus.dir)
   set_target_properties(dpMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(dpMallocPlus PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(dpMallocPlus ${MPI_LIBRARIES})
   target_link_libraries(dpMallocPlus ${OpenCL_LIBRARIES})
   target_link_libraries(dpMallocPlus ezcl)
   #if (Boost_FOUND AND MPI_FOUND)
   #   message(STATUS "Adding j7 support ot MallocPlus")
   #   target_link_libraries(dMallocPlus j7)
   #endif()
#   install(TARGETS dpMallocPlus DESTINATION lib)
endif(OpenCL_FOUND AND MPI_FOUND)

########### install files ###############

#install(FILES  MallocPlus.h DESTINATION include/MallocPlus)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/MallocPlus.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/MallocPlus.pdf)
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
   add_custom_target(MallocPlus_doc
                     COMMAND make MallocPlus_pdfdoc
                     COMMAND make MallocPlus_htmldoc
                     COMMAND make MallocPlus_singlehtmldoc)
   set_target_properties(MallocPlus_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/MallocPlus.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(MallocPlus_pdfdoc DEPENDS _build/latex/MallocPlus.pdf)
   set_target_properties(MallocPlus_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(MallocPlus_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(MallocPlus_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(MallocPlus_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(MallocPlus_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")
