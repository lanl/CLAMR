cmake_minimum_required(VERSION 2.6)

project (POWERPARSER)
# The version number.
set (POWERPARSER_VERSION_MAJOR 2)
set (POWERPARSER_VERSION_MINOR 0.7)
set (POWERPARSER_VERSION "${POWERPARSER_VERSION_MAJOR}.${POWERPARSER_VERSION_MINOR}")

enable_language (Fortran)

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
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -vec-report3")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
endif()

if(DEFINED ENV{FC})
   if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#    set (CMAKE_Fortran_FLAGS "-cpp")
     set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
     set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
   elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
#    set (CMAKE_Fortran_FLAGS "-fpp -f77rtl")
     set (CMAKE_Fortran_FLAGS "-f77rtl")
     set (CMAKE_Fortran_FLAGS_RELEASE " -O3")
     set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g")
   elseif (Fortran_COMPILER_NAME MATCHES "g77")
     set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
     set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
   else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
     message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
     message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
     message ("No optimized Fortran compiler flags are known, we just try -O2...")
     set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
     set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
   endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
endif(DEFINED ENV{FC})


if (MIC_NATIVE)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")

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
set(H_SRCS   Cmd.hh   Function.hh    Parser_utils.hh  Variable.hh      Word.hh
             Comm.hh  PowerParser.hh Parser_math.hh   Restartblock.hh  Whenthen.hh)

set(CXX_SRCS Cmd.cc   Function.cc    Parser_utils.cc  Variable.cc      Word.cc
             Comm.cc  PowerParser.cc Parser_math.cc   Restartblock.cc  Whenthen.cc)

if(DEFINED ENV{FC})
   set(H_SRCS ${H_SRCS} parser_iface.hh)
   set(CXX_SRCS ${CXX_SRCS} parser_iface.cc)
   set(F90_SRCS FParser.F90)
endif(DEFINED ENV{FC})

set(DOC_SRCS docs/index.rst docs/About.rst docs/GettingStarted.rst
             docs/DevelopersGuide.rst docs/UsersGuide.rst)

set(PDFDOC_SRCS docs/indexpdf.rst docs/About.rst docs/GettingStarted.rst
             docs/DevelopersGuide.rst docs/UsersGuide.rst)

########### PowerParser target ###############
set(PowerParser_LIB_SRCS ${CXX_SRCS} ${H_SRCS} ${F90_SRCS})

add_library(PowerParser STATIC ${PowerParser_LIB_SRCS})

set_target_properties(PowerParser PROPERTIES VERSION ${POWERPARSER_VERSION} SOVERSION ${POWERPARSER_VERSION_MAJOR})
install(TARGETS PowerParser DESTINATION lib)

########### pPowerParser target ###############
if (MPI_FOUND)
   set(pPowerParser_LIB_SRCS ${CXX_SRCS} ${H_SRCS} ${F90_SRCS})

   add_library(pPowerParser STATIC ${PowerParser_LIB_SRCS})

   set_target_properties(pPowerParser PROPERTIES VERSION ${POWERPARSER_VERSION} SOVERSION ${POWERPARSER_VERSION_MAJOR})
   set_target_properties(pPowerParser PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
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

