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
   #set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

set (CMAKE_C_STANDARD 99)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -fstrict-aliasing -ftree-vectorize")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopt-info-vec-optimized")
   if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "4.9")
     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp-simd")
   endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "4.9")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
   elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mcpu=powerpc64le")
   else ()
     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
   endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -vec-report3")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_CXX_FLAGS_RELEASE "-qopt-zmm-usage=high ${CMAKE_CXX_FLAGS_RELEASE}")
   endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -ansi-alias -vecabi=cmdtarget -xHost -qopenmp-simd -qopt-report=5 -qopt-report-phase=vec,loop")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Cray")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")
endif()

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# set (CMAKE_Fortran_FLAGS "-cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
# set (CMAKE_Fortran_FLAGS "-fpp -f77rtl")
  set (CMAKE_Fortran_FLAGS "-f77rtl")
  set (CMAKE_Fortran_FLAGS_RELEASE " -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
# set (CMAKE_Fortran_FLAGS "-cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
# set (CMAKE_Fortran_FLAGS "-cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
  set (CMAKE_Fortran_FLAGS "-qextname")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
else ()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: ${Fortran_COMPILER_NAME} ${CMAKE_Fortran_COMPILER_ID}")
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ()

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

#if(DEFINED ENV{FC})
   set(H_SRCS ${H_SRCS} parser_iface.hh)
   set(CXX_SRCS ${CXX_SRCS} parser_iface.cc)
   set(F90_SRCS FParser.F90)
#endif(DEFINED ENV{FC})

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

