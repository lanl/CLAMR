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
# MPE_INCLUDE_PATH  Include path(s) for MPI logging with MPE.
# MPE_LINK_FLAGS    Linking flags for MPI logging with MPE.
# MPE_LIBRARIES     Libraries to link against for MPI logging with MPE.
#
# This module will auto-detect these setting by looking for an MPE
# compiler (mpecc) and use the -show flag to retrieve compiler options.
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

if (MPE_INCLUDE_PATH AND MPE_LIBRARIES)
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
endif (MPE_INCLUDE_PATH AND MPE_LIBRARIES)

if (DEFINED ENV{HOME})
   set(HOME $ENV{HOME})
endif()

if (MPE_INCLUDE_PATH)
else (MPE_INCLUDE_PATH)
   FIND_PATH(MPE_INCLUDE_PATH mpe.h PATHS 
      /usr/local/mpe/include
      /usr/local/mpe-1.9.1/include
      ${HOME}/mpe/include
      ${HOME}/mpe-1.9.1/include
   )
endif(MPE_INCLUDE_PATH)

if (MPE_LIBRARIES)
else (MPE_LIBRARIES)
   FIND_LIBRARY(MPE_LIBRARIES mpe PATHS
      /usr/local/mpe/lib
      /usr/local/mpe-1.9.1/lib      
      ${HOME}/mpe/lib
      ${HOME}/mpe-1.9.1/lib      
   )
endif (MPE_LIBRARIES)

if (MPE_INCLUDE_PATH AND MPE_LIBRARIES)
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
    MPE_ALL_INCLUDE_PATHS
    "${MPE_COMPILE_CMDLINE}")
  set(MPE_INCLUDE_PATH_WORK)
  foreach(IPATH ${MPE_ALL_INCLUDE_PATHS})
    string(REGEX REPLACE "^-I" "" IPATH ${IPATH})
    string(REGEX REPLACE "//" "/" IPATH ${IPATH})
    list(APPEND MPE_INCLUDE_PATH_WORK ${IPATH})
  endforeach(IPATH)

  # Extract linker paths from the link command line
  string(REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")"
    MPE_ALL_LINK_PATHS
    "${MPE_COMPILE_CMDLINE}")
  set(MPE_LINK_PATH)
  foreach(LPATH ${MPE_ALL_LINK_PATHS})
    string(REGEX REPLACE "^-L" "" LPATH ${LPATH})
    string(REGEX REPLACE "//" "/" LPATH ${LPATH})
    list(APPEND MPE_LINK_PATH ${LPATH})
  endforeach(LPATH)

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
    set(MPE_LIB "MPE_LIB-NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
    find_library(MPE_LIB ${LIB} HINTS ${MPE_LINK_PATH})
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
  set(MPE_INCLUDE_PATH ${MPE_INCLUDE_PATH_WORK}
    CACHE STRING "MPE log include path" FORCE)
  set(MPE_LINK_FLAGS ${MPE_LINK_FLAGS_WORK}
    CACHE STRING "MPE log linking flags" FORCE)
  set(MPE_LIBRARIES ${MPE_LIBRARIES_WORK}
    CACHE PATH "MPE log libraries" FORCE)
endif (MPE_INCLUDE_PATH AND MPE_LIBRARIES)

if (MPE_INCLUDE_PATH AND MPE_LIBRARIES)
  set(MPE_FOUND TRUE)
else (MPE_INCLUDE_PATH AND MPE_LIBRARIES)
  set(MPE_FOUND FALSE)
endif (MPE_INCLUDE_PATH AND MPE_LIBRARIES)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments
find_package_handle_standard_args(
  MPE
  DEFAULT_MSG
  MPE_LIBRARIES
  MPE_INCLUDE_PATH
  )

mark_as_advanced(
  MPE_COMPILE_FLAGS
  MPE_INCLUDE_PATH
  MPE_LINK_FLAGS
  MPE_LIBRARIES
  )

