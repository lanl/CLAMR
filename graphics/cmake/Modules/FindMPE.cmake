# - Multi-Processing Environment (MPE) module.
#
# The Multi-Processing Environment (MPE) is an extention to MPI that
# provides programmers with a suite of performance analysis tools for
# their MPI programs.  These tools include a set of profiling libraries,
# a set of utility programs, and a set of graphical tools.  This module
# helps you find the libraries and includes.
#
# This module will set the following variables:
# MPE_FOUND             TRUE if we have found MPE
# MPE_COMPILE_FLAGS Compilation flags for MPI logging with MPE.
# MPE_INCLUDE_DIR  Include path(s) for MPI logging with MPE.
# MPE_LINK_FLAGS    Linking flags for MPI logging with MPE.
# MPE_LIBRARIES     Libraries to link against for MPI logging with MPE.
# MPE_NOMPI_LIBRARIES     Libraries to link against for MPI logging with MPE.
#
# This module will auto-detect these setting by looking for the clog2alog
# utility and use the path to discover the include and library locations
#
# Note that this module does not attempt to ensure that the version
# of MPE you are using is compatible with the version of MPI that
# you are using (or even that you are using MPI at all).
#

## Copyright 2011 Sandia Coporation
## Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
## the U.S. Government retains certain rights in this software.
##
## This source code is released under the New BSD License.
#
#  This file changed substantially by Bob Robey Los Alamos National Laboratory

find_program(MPE_CLOG2ALOG
  NAMES clog2alog
  DOC "MPE clog2alog command.  Used only to detect MPE install location."
  )

if (MPE_INCLUDE_DIR AND MPE_LIBRARIES)
  # Do nothing: we already have the necessary options.
elseif (MPE_CLOG2ALOG)
  exec_program(${MPE_CLOG2ALOG}
    ARGS -show -mpilog
    OUTPUT_VARIABLE MPE_COMPILE_CMDLINE
    RETURN_VALUE MPE_COMPILE_RETURN
    )
  if (NOT MPE_COMPILE_RETURN EQUAL 0)
    message(STATUS, "Unable to determine MPE from MPE driver ${MPE_CLOG2ALOG}")
  endif (NOT MPE_COMPILE_RETURN EQUAL 0)
endif (MPE_INCLUDE_DIR AND MPE_LIBRARIES)

if (DEFINED ENV{HOME})
   set(HOME $ENV{HOME})
endif()

if (DEFINED ENV{MPEHOME})
   if (MPE_INCLUDE_DIR)
   else (MPE_INCLUDE_DIR)
      set (MPE_INCLUDE_DIR $ENV{MPEHOME}/include)
   endif(MPE_INCLUDE_DIR)
   if (MPE_LIBRARIES)
   else (MPE_LIBRARIES)
      FIND_LIBRARY(MPE_LIBRARIES mpe DIRS
        $ENV{MPEHOME}/lib
      )
   endif(MPE_LIBRARIES)
endif()

if (MPE_INCLUDE_DIR)
else (MPE_INCLUDE_DIR)
   FIND_PATH(MPE_INCLUDE_DIR mpe.h DIRS 
      /usr/local/mpe/include
      /usr/local/mpe-1.9.1/include
      ${HOME}/mpe/include
      ${HOME}/mpe-1.9.1/include
      ${CMAKE_CURRENT_BINARY_DIR}/mpe/include
   )
endif(MPE_INCLUDE_DIR)

if (MPE_LIBRARIES)
else (MPE_LIBRARIES)
   FIND_LIBRARY(MPE_LIBRARIES mpe DIRS
      /usr/local/mpe/lib
      /usr/local/mpe-1.9.1/lib      
      ${HOME}/mpe/lib
      ${HOME}/mpe-1.9.1/lib      
      ${CMAKE_CURRENT_BINARY_DIR}/mpe/lib
      ${CMAKE_CURRENT_BINARY_DIR}/mpe-1.9.1/lib
   )
endif (MPE_LIBRARIES)

if (MPE_INCLUDE_DIR AND MPE_LIBRARIES)
elseif (MPE_COMPILE_CMDLINE)
  # Extract compile flags from the compile command line.
  string(REGEX MATCHALL "-D([^\" ]+|\"[^\"]+\")"
    MPE_ALL_COMPILE_FLAGS
    "${MPE_COMPILE_CMDLINE}")
  set(MPE_COMPILE_FLAGS_WORK)
  foreach(FLAG ${MPE_ALL_COMPILE_FLAGS})
    if (MPE_COMPILE_FLAGS_WORK)
      set(MPE_COMPILE_FLAGS_WORK "${MPE_COMPILE_FLAGS_WORK} ${FLAG}")
    else(MPE_COMPILE_FLAGS_WORK)
      set(MPE_COMPILE_FLAGS_WORK ${FLAG})
    endif(MPE_COMPILE_FLAGS_WORK)
  endforeach(FLAG)

  # Extract include paths from compile command line
  string(REGEX MATCHALL "-I([^\" ]+|\"[^\"]+\")"
    MPE_ALL_INCLUDE_DIRS
    "${MPE_COMPILE_CMDLINE}")
  set(MPE_INCLUDE_DIR_WORK)
  foreach(IDIR ${MPE_ALL_INCLUDE_DIRS})
    string(REGEX REPLACE "^-I" "" IDIR ${IDIR})
    string(REGEX REPLACE "//" "/" IDIR ${IDIR})
    list(APPEND MPE_INCLUDE_DIR_WORK ${IDIR})
  endforeach(IDIR)

  # Extract linker paths from the link command line
  string(REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")"
    MPE_ALL_LINK_DIRS
    "${MPE_COMPILE_CMDLINE}")
  set(MPE_LINK_DIR)
  foreach(LDIR ${MPE_ALL_LINK_DIRS})
    string(REGEX REPLACE "^-L" "" LDIR ${LDIR})
    string(REGEX REPLACE "//" "/" LDIR ${LDIR})
    list(APPEND MPE_LINK_DIR ${LDIR})
  endforeach(LDIR)

  # Extract linker flags from the link command line
  string(REGEX MATCHALL "-Wl,([^\" ]+|\"[^\"]+\")"
    MPE_ALL_LINK_FLAGS
    "${MPE_COMPILE_CMDLINE}")
  set(MPE_LINK_FLAGS_WORK)
  foreach(FLAG ${MPE_ALL_LINK_FLAGS})
    if (MPE_LINK_FLAGS_WORK)
      set(MPE_LINK_FLAGS_WORK "${MPE_LINK_FLAGS_WORK} ${FLAG}")
    else(MPE_LINK_FLAGS_WORK)
      set(MPE_LINK_FLAGS_WORK ${FLAG})
    endif(MPE_LINK_FLAGS_WORK)
  endforeach(FLAG)

  # Extract the set of libraries to link against from the link command
  # line
  string(REGEX MATCHALL "-l([^\" ]+|\"[^\"]+\")"
    MPE_LIBNAMES
    "${MPE_COMPILE_CMDLINE}")

  # Determine full path names for all of the libraries that one needs
  # to link against in an MPI program
  set(MPE_LIBRARIES_WORK)
  foreach(LIB ${MPE_LIBNAMES})
    string(REGEX REPLACE "^-l" "" LIB ${LIB})
    set(MPE_LIB "MPE_LIB-NOTFOUND" CACHE FILEDIR "Cleared" FORCE)
    find_library(MPE_LIB ${LIB} HINTS ${MPE_LINK_DIR})
    if (MPE_LIB)
      list(APPEND MPE_LIBRARIES_WORK ${MPE_LIB})
    else (MPE_LIB)
      message(SEND_ERROR "Unable to find MPE library ${LIB}")
    endif (MPE_LIB)
  endforeach(LIB)
  set(MPE_LIB "MPE_LIB-NOTFOUND"
    CACHE INTERNAL "Scratch variable for MPI detection" FORCE)

  # Set up all of the appropriate cache entries
  set(MPE_COMPILE_FLAGS ${MPE_COMPILE_FLAGS_WORK}
    CACHE STRING "MPE log compilation flags" FORCE)
  set(MPE_INCLUDE_DIR ${MPE_INCLUDE_DIR_WORK}
    CACHE STRING "MPE log include path" FORCE)
  set(MPE_LINK_FLAGS ${MPE_LINK_FLAGS_WORK}
    CACHE STRING "MPE log linking flags" FORCE)
  set(MPE_LIBRARIES ${MPE_LIBRARIES_WORK}
    CACHE DIR "MPE log libraries" FORCE)
endif (MPE_INCLUDE_DIR AND MPE_LIBRARIES)

if (MPE_INCLUDE_DIR AND MPE_LIBRARIES)
  set(MPE_FOUND TRUE)
else (MPE_INCLUDE_DIR AND MPE_LIBRARIES)
  set(MPE_FOUND FALSE)
endif (MPE_INCLUDE_DIR AND MPE_LIBRARIES)

GET_FILENAME_COMPONENT(MPE_LIB_DIR ${MPE_LIBRARIES} PATH)
set (MPE_NOMPI_LIBRARIES ${MPE_LIB_DIR}/libmpe_nompi.a)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments
find_package_handle_standard_args(
  MPE
  DEFAULT_MSG
  MPE_LIBRARIES
  MPE_INCLUDE_DIR
  )

mark_as_advanced(
  MPE_COMPILE_FLAGS
  MPE_INCLUDE_DIR
  MPE_LINK_FLAGS
  MPE_LIBRARIES
  )

