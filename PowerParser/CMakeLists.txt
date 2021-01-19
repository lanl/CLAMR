cmake_minimum_required(VERSION 2.6)

project (POWERPARSER CXX C Fortran)
# The version number.
set (POWERPARSER_VERSION_MAJOR 2)
set (POWERPARSER_VERSION_MINOR 0.7)
set (POWERPARSER_VERSION "${POWERPARSER_VERSION_MAJOR}.${POWERPARSER_VERSION_MINOR}")

#  FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

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
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
    set( CMAKE_VECTOR_NOVEC TRUE)
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")


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

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# set (CMAKE_Fortran_FLAGS "-cpp")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
# set (CMAKE_Fortran_FLAGS "-fpp -f77rtl")
  set (CMAKE_Fortran_FLAGS "-f77rtl")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
# set (CMAKE_Fortran_FLAGS "-cpp")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
# set (CMAKE_Fortran_FLAGS "-cpp")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
  set (CMAKE_Fortran_FLAGS "-qextname")
else ()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: ${Fortran_COMPILER_NAME} ${CMAKE_Fortran_COMPILER_ID}")
endif ()

if(NOT KOKKOS_ENABLE)
   if (NOT CMAKE_VECTOR_C_FLAGS)
      find_package(Vector)
   endif (NOT CMAKE_VECTOR_C_FLAGS)
   if (CMAKE_VECTOR_VERBOSE)
       set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
       set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
       if (NOT CMAKE_VECTOR_Fortran_FLAGS)
           set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
       endif (NOT CMAKE_VECTOR_Fortran_FLAGS)
   endif (CMAKE_VECTOR_VERBOSE)
   if (CMAKE_VECTOR_NOVEC)
       set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
       set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
       if (NOT CMAKE_VECTOR_Fortran_FLAGS)
          set(VECTOR_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS}")
       endif (NOT CMAKE_VECTOR_Fortran_FLAGS)
   endif (CMAKE_VECTOR_NOVEC)
endif(NOT KOKKOS_ENABLE)

set(H_SRCS   Cmd.hh   Function.hh    Parser_utils.hh  Variable.hh      Word.hh
             Comm.hh  PowerParser.hh Parser_math.hh   Restartblock.hh  Whenthen.hh)

set(CXX_SRCS Cmd.cc   Function.cc    Parser_utils.cc  Variable.cc      Word.cc
             Comm.cc  PowerParser.cc Parser_math.cc   Restartblock.cc  Whenthen.cc)

#if(DEFINED ENV{FC})
   set(H_SRCS ${H_SRCS} parser_iface.hh)
   set(CXX_SRCS ${CXX_SRCS} parser_iface.cc)
   set(F90_SRCS FParser.F90)
#endif(DEFINED ENV{FC})

set_source_files_properties(CXX_SRCS PROPERTIES COMPILE_FLAGS "${VECTOR_CXX_FLAGS}")
set_source_files_properties(F90_SRCS PROPERTIES COMPILE_FLAGS "${VECTOR_Fortran_FLAGS}")

if (NOT DEFINED MPI_FOUND)
   find_package(MPI)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

include_directories(${CMAKE_SOURCE_DIRECTORY})

enable_testing()

add_subdirectory(tests)

########### global settings ###############

set(DOC_SRCS docs/index.rst docs/About.rst docs/GettingStarted.rst
             docs/DevelopersGuide.rst docs/UsersGuide.rst)

set(PDFDOC_SRCS docs/indexpdf.rst docs/About.rst docs/GettingStarted.rst
             docs/DevelopersGuide.rst docs/UsersGuide.rst)

########### PowerParser target ###############
set(PowerParser_LIB_SRCS ${CXX_SRCS} ${H_SRCS} ${F90_SRCS})

add_library(PowerParser STATIC ${PowerParser_LIB_SRCS})

set_target_properties(PowerParser PROPERTIES VERSION ${POWERPARSER_VERSION} SOVERSION ${POWERPARSER_VERSION_MAJOR})
set_target_properties(PowerParser PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/PowerParser.dir)

install(TARGETS PowerParser DESTINATION lib)

########### pPowerParser target ###############
if (MPI_FOUND)
   set(pPowerParser_LIB_SRCS ${CXX_SRCS} ${H_SRCS} ${F90_SRCS})

   add_library(pPowerParser STATIC ${PowerParser_LIB_SRCS})

   set_target_properties(pPowerParser PROPERTIES VERSION ${POWERPARSER_VERSION} SOVERSION ${POWERPARSER_VERSION_MAJOR})
   set_target_properties(pPowerParser PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   set_target_properties(pPowerParser PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/pPowerParser.dir)
   target_link_libraries(pPowerParser ${MPI_LIBRARIES})

   install(TARGETS pPowerParser DESTINATION lib)
endif (MPI_FOUND)

########### install files ###############

install(FILES PowerParser.hh DESTINATION include)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
  set(doc_prefix "")
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
  set(doc_prefix PowerParser_)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

add_custom_target(${doc_prefix}doc
                  COMMAND make ${doc_prefix}pdfdoc
                  COMMAND make ${doc_prefix}htmldoc
                  COMMAND make ${doc_prefix}singlehtmldoc)
set_target_properties(${doc_prefix}doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_custom_command(OUTPUT _build/latex/PowerParser.pdf
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                   COMMAND make latexpdf >& pdfdoc.out
                   DEPENDS ${PDFDOC_SRCS})
add_custom_target(${doc_prefix}pdfdoc DEPENDS _build/latex/PowerParser.pdf)
set_target_properties(${doc_prefix}pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_custom_command(OUTPUT _build/html/UsersGuide.html
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                   COMMAND make html >& htmldoc.out
                   DEPENDS ${DOC_SRCS})
add_custom_target(${doc_prefix}htmldoc DEPENDS _build/html/UsersGuide.html)
set_target_properties(${doc_prefix}htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_custom_command(OUTPUT _build/singlehtml/index.html
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                   COMMAND make singlehtml >& singlehtmldoc.out
                   DEPENDS ${DOC_SRCS})
add_custom_target(${doc_prefix}singlehtmldoc DEPENDS _build/singlehtml/index.html)
set_target_properties(${doc_prefix}singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "FParser_module.mod;install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")

