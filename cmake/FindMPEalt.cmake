# This module finds if mpe is installed and determines where the include files,
# libraries and the mpecc executable are. The code sets the following variables:
#
#  MPE_MPILOG_INCLUDES  = paths to the mpe includes
#  MPE_MPILOG_LIBRARIES = paths to all mpe libraries
#  MPECC_EXECUTABLE     = the mpecc executable

# mpe depends on mpi (only mpilog functionality is queried)
IF(NOT MPE_FOUND)
	IF(${MPE_FIND_REQUIRED})
		FIND_PACKAGE(MPI REQUIRED)
	ELSE()
		FIND_PACKAGE(MPI)
	ENDIF()

	IF(NOT MPE_FIND_QUIETLY)
		MESSAGE(STATUS "Looking for mpe")
	ENDIF()

	IF(MPI_LIBRARY)
		FIND_PROGRAM(MPECC_EXECUTABLE NAMES mpecc HINTS ${MPE_ROOT})

		IF(MPECC_EXECUTABLE)
			# determine libs via kconfig tool
			EXECUTE_PROCESS(COMMAND ${MPECC_EXECUTABLE} -mpilog -show OUTPUT_VARIABLE MPE_MPILOG_ALL_FLAGS ERROR_QUIET)
			IF(NOT MPE_MPILOG_ALL_FLAGS)
				MESSAGE(FATAL_ERROR "Failure executing ${MPECC_EXECUTABLE}")
			ENDIF()

			# extract link path
			STRING(REGEX MATCHALL " -L([^\" \n]+|\"[^\"]+\")" MPE_MPILOG_ALL_LINK_PATH "${MPE_MPILOG_ALL_FLAGS}")
			FOREACH(LPATH ${MPE_MPILOG_ALL_LINK_PATH})
				STRING(REGEX REPLACE "^ -L" "" LPATH ${LPATH})
				STRING(REGEX REPLACE "//" "/" LPATH ${LPATH})
				LIST(APPEND MPE_MPILOG_LINK_PATH ${LPATH})
			ENDFOREACH()
			SET(MPE_MPILOG_ALL_LINK_PATH)

			# extract libraries
			STRING(REGEX MATCHALL "-l([^\" \n]+|\"[^\"]+\")" MPE_MPILOG_ALL_LIBNAMES "${MPE_MPILOG_ALL_FLAGS}")
			FOREACH(LIB ${MPE_MPILOG_ALL_LIBNAMES})
				STRING(REGEX REPLACE "^-l" "" LIB ${LIB})
				SET(MPE_LIB "MPE_LIB-NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
				FIND_LIBRARY(MPE_LIB ${LIB} HINTS ${MPE_MPILOG_LINK_PATH})
				IF(MPE_LIB)
					LIST(APPEND MPE_MPILOG_LIBRARIES ${MPE_LIB})
				ELSE()
					MESSAGE(FATAL_ERROR "Unable to find mpe library ${LIB}")
				ENDIF()
			ENDFOREACH()

			# extract include paths
			STRING(REGEX MATCHALL "-I([^\" \n]+|\"[^\"]+\")" MPE_MPILOG_ALL_INCLUDES "${MPE_MPILOG_ALL_FLAGS}")
			FOREACH(IPATH ${MPE_MPILOG_ALL_INCLUDES})
				STRING(REGEX REPLACE "^-I" "" IPATH ${IPATH})
				STRING(REGEX REPLACE "//" "/" IPATH ${IPATH})
				LIST(APPEND MPE_MPILOG_INCLUDES ${IPATH})
			ENDFOREACH()

			# set the needed cache variables
			SET(MPE_FOUND TRUE CACHE BOOL "Indicates whether mpe has been found" FORCE)
			SET(MPE_MPILOG_INCLUDES ${MPE_MPILOG_INCLUDES} CACHE LIST "Paths to the mpe includes (mpilog)" FORCE)
			SET(MPE_MPILOG_LIBRARIES ${MPE_MPILOG_LIBRARIES} CACHE LIST "Paths to the mpe libraries (mpilog)" FORCE)
			SET(MPECC_EXECUTABLE ${MPECC_EXECUTABLE} CACHE FILEPATH "The mpecc executable" FORCE)
			SET(MPE_LIB CACHE STRING "Scratch variable cmake refuses to remove")
		ENDIF()
	ENDIF()

	IF(NOT MPE_FIND_QUIETLY)
		IF(MPE_FOUND)
			MESSAGE(STATUS "Looking for mpe - done")
		ELSE()
			MESSAGE(STATUS "Looking for mpe - failed")
		ENDIF()
	ENDIF()
ENDIF()

IF(NOT MPE_FOUND AND MPE_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Mpe not found.\nConsider setting the MPE_ROOT variable to point to the location of the mpe installation.")
ENDIF()

