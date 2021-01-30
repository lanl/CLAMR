cmake_minimum_required(VERSION 3.1.3)

project (CLAMR CXX C)
# The version number.
set (CLAMR_VERSION_MAJOR 2)
set (CLAMR_VERSION_MINOR 0.7)
set (CLAMR_VERSION "${CLAMR_VERSION_MAJOR}.${CLAMR_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

if (NOT CMAKE_BUILD_TYPE)
   #set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

if (CMAKE_PRECISION_CHECK)
   add_definitions(-DPRECISION_CHECK_GRAPHICS=1)
   add_definitions(-DPRECISION_CHECK_STATS=1)
   add_definitions(-DPRECISION_CHECK=1.0e-7)
endif (CMAKE_PRECISION_CHECK)

if (CMAKE_PRECISION_CHECK_APPLIED)
   if (CMAKE_PRECISION_CHECK_APPLIED STREQUAL "best")
      add_definitions(-DPRECISION_CHECK_BEST_PARENTHESIS=1)
   else (CMAKE_PRECISION_CHECK_APPLIED STREQUAL "best")
      add_definitions(-DPRECISION_CHECK_WITH_PARENTHESIS=1)
   endif(CMAKE_PRECISION_CHECK_APPLIED STREQUAL "best")
endif (CMAKE_PRECISION_CHECK_APPLIED)

enable_testing()
# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

#set (CMAKE_CXX_EXTENSIONS off)
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_C_STANDARD 99)

#message("Release CFLAGS are ${CMAKE_C_FLAGS_RELEASE}")
#message("RelWithDebInfo CFLAGS are ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
#message("Debug CFLAGS are ${CMAKE_C_FLAGS_DEBUG}")
#message("CFLAGS are ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO})
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
string(REPLACE -O2 -O3 CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
#message("Release CFLAGS are ${CMAKE_C_FLAGS_RELEASE}")
#message("RelWithDebInfo CFLAGS are ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
#message("Debug CFLAGS are ${CMAKE_C_FLAGS_DEBUG}")
#message("CFLAGS are ${CMAKE_C_FLAGS}")

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
 
option (KOKKOS_ENABLE "Enable KOKKOS support" OFF)
mark_as_advanced (KOKKOS_ENABLE)

#set (Kokkos_ENABLE_CUDA_LAMBDA on)

if(NOT KOKKOS_ENABLE)
   find_package(Vector)
   if (CMAKE_VECTOR_NOVEC)
       set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
       set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
       set(VECTOR_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS}")
   endif (CMAKE_VECTOR_NOVEC)
   if (CMAKE_VECTOR_VERBOSE)
       set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
       set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
       set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
   endif (CMAKE_VECTOR_VERBOSE)
endif(NOT KOKKOS_ENABLE)

set_source_files_properties(state.cpp PROPERTIES COMPILE_FLAGS "${VECTOR_CXX_FLAGS}")

include_directories(${PROJECT_SOURCE_DIR})

if (NOT APPLE OR "${CMAKE_SYSTEM_VERSION}" LESS 18)
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
endif (NOT APPLE OR "${CMAKE_SYSTEM_VERSION}" LESS 18)

find_package(MPI)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if(NOT KOKKOS_ENABLE)
   find_package( OpenMP)
   if(OPENMP_FOUND)
      message("OPENMP FOUND ${OpenMP_C_FLAGS}")
   else (OPENMP_FOUND)
      message("OPENMP NOT FOUND")
   endif(OPENMP_FOUND)
endif(NOT KOKKOS_ENABLE)

#find_package(ImageMagick COMPONENTS MagickCore MagickWand)
find_package(ImageMagick COMPONENTS MagickWand)
#   ImageMagick_FOUND                    - TRUE if all components are found.
#   ImageMagick_INCLUDE_DIRS             - Full paths to all include dirs.
#   ImageMagick_LIBRARIES                - Full paths to all libraries.
#   ImageMagick_<component>_FOUND        - TRUE if <component> is found.
#   ImageMagick_<component>_INCLUDE_DIRS - Full path to <component> include dirs.
#   ImageMagick_<component>_LIBRARIES    - Full path to <component> libraries.
if (ImageMagick_FOUND)
   if (ImageMagick_MagickWand_FOUND)
      add_definitions(-DHAVE_MAGICKWAND)
      add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)
      add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
      if (ImageMagick_VERSION_STRING VERSION_LESS 7.0.0)
          add_definitions(-DIMAGEMAGICK_VERSION_6)
      else (ImageMagick_VERSION_STRING VERSION_LESS 7.0.0)
          add_definitions(-DIMAGEMAGICK_VERSION_7)
      endif (ImageMagick_VERSION_STRING VERSION_LESS 7.0.0)

      add_definitions(-DImageMagick_VERSION_STRING=${ImageMagick_VERSION_STRING})
#     add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)
#     add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
      include_directories(${ImageMagick_INCLUDE_DIRS})
   endif (ImageMagick_MagickWand_FOUND)

#   message("ImageMagick_FOUND                   ${ImageMagick_FOUND}")
#   message("ImageMagick_INCLUDE_DIRS            ${ImageMagick_INCLUDE_DIRS}")
#   message("ImageMagick_LIBRARIES               ${ImageMagick_LIBRARIES}")
##  message("ImageMagick_MagickCore_FOUND        ${ImageMagick_MagickCore_FOUND}")
##  message("ImageMagick_MagickCore_INCLUDE_DIRS ${ImageMagick_MagickCore_INCLUDE_DIRS}")
##  message("ImageMagick_MagickCore_LIBRARIES    ${ImageMagick_MagickCore_LIBRARIES}")
#   message("ImageMagick_MagickWand_FOUND        ${ImageMagick_MagickWand_FOUND}")
#   message("ImageMagick_MagickWand_INCLUDE_DIRS ${ImageMagick_MagickWand_INCLUDE_DIRS}")
#   message("ImageMagick_MagickWand_LIBRARIES    ${ImageMagick_MagickWand_LIBRARIES}")
endif (ImageMagick_FOUND)


#-----------------------------------------------------------------------------
# HDF5 support
#-----------------------------------------------------------------------------

option (HDF5_ENABLE "Enable HDF5 support" OFF)
mark_as_advanced (HDF5_ENABLE)

if(HDF5_ENABLE)
  if (NOT HDF5_HDF5_HEADER)
    if (BUILD_SHARED_LIBS)
      set (FIND_HDF_COMPONENTS C shared)
    else (BUILD_SHARED_LIBS)
      set (FIND_HDF_COMPONENTS C static)
    endif (BUILD_SHARED_LIBS)
    message (STATUS "HDF5 find comps: ${FIND_HDF_COMPONENTS}")
    set (SEARCH_PACKAGE_NAME ${HDF5_PACKAGE_NAME})

    find_package (HDF5 NAMES ${SEARCH_PACKAGE_NAME} QUIET COMPONENTS ${FIND_HDF_COMPONENTS})
    message (STATUS "HDF5 C libs:${HDF5_FOUND} static:${HDF5_static_C_FOUND} and shared:${HDF5_shared_C_FOUND}")
    if (HDF5_FOUND)
      if (NOT HDF5_static_C_FOUND AND NOT HDF5_shared_C_FOUND)
        # check for single binary package
        set (FIND_HDF_COMPONENTS C)
        message (STATUS "HDF5 find comps: ${FIND_HDF_COMPONENTS}")

        find_package (HDF5 NAMES ${SEARCH_PACKAGE_NAME} COMPONENTS ${FIND_HDF_COMPONENTS})
        message (STATUS "HDF5 libs:${HDF5_FOUND} C:${HDF5_C_FOUND}")
        set (HDF5_LIBS ${HDF5_LIBS} ${HDF5_LIBRARIES})
          if (HDF5_BUILD_SHARED_LIBS AND USE_SHARED_LIBS)
          add_definitions (-DH5_BUILT_AS_DYNAMIC_LIB)
        else (HDF5_BUILD_SHARED_LIBS AND USE_SHARED_LIBS)
          add_definitions (-DH5_BUILT_AS_STATIC_LIB)
        endif (HDF5_BUILD_SHARED_LIBS AND USE_SHARED_LIBS)
      else (NOT HDF5_static_C_FOUND AND NOT HDF5_shared_C_FOUND)
        if (HDF5_BUILD_SHARED_LIBS AND USE_SHARED_LIBS AND HDF5_shared_C_FOUND)
          set (HDF5_LIBS ${HDF5_LIBS} ${HDF5_C_SHARED_LIBRARY} ${HDF5_HL_SHARED_LIBRARY})
        else (HDF5_static_C_FOUND)
          set (HDF5_LIBS ${HDF5_LIBS} ${HDF5_C_STATIC_LIBRARY} ${HDF5_HL_STATIC_LIBRARY})
        endif (HDF5_BUILD_SHARED_LIBS AND USE_SHARED_LIBS AND HDF5_shared_C_FOUND)
      endif (NOT HDF5_static_C_FOUND AND NOT HDF5_shared_C_FOUND)
    else (HDF5_FOUND)
      set (FIND_HDF_COMPONENTS C)
      find_package (HDF5 COMPONENTS ${FIND_HDF_COMPONENTS}) # Legacy find

      if (USE_SHARED_LIBS AND EXISTS "${HDF5_LIBRARY_DIRS}/libhdf5.so")
        set (HDF5_BUILD_SHARED_LIBS 1)
      add_definitions (-DH5_BUILT_AS_DYNAMIC_LIB)
      else ()
        set (HDF5_BUILD_SHARED_LIBS 0)
        add_definitions (-DH5_BUILT_AS_STATIC_LIB)
      endif ()
      set (HDF5_LIBS ${HDF5_LIBS} ${HDF5_LIBRARIES})

    endif (HDF5_FOUND)
    set (HDF5_PACKAGE_NAME ${SEARCH_PACKAGE_NAME})

    if (HDF5_FOUND)
      set (HDF5_HAVE_H5PUBCONF_H 1)
      set (HDF5_HAVE_HDF5 1)
      set (HDF5_HDF5_HEADER "h5pubconf.h")
      set (HDF5_INCLUDE_DIR_GEN ${HDF5_INCLUDE_DIR})
      if (NOT HDF5_VERSION_MAJOR)
        # HDF5 was build with autotools
        if(EXISTS "${HDF5_LIBRARY_DIRS}/../include/H5pubconf.h")
          file(STRINGS ${HDF5_LIBRARY_DIRS}/../include/H5pubconf.h "HDF5_VERSION_MAJOR REGEX")
        endif()
        string(REGEX MATCH "[0-9]+.[0-9]+" HDF5_VERSION_STRING "${HDF5_VERSION_MAJOR}")
      endif ()
      message (STATUS "HDF5-${HDF5_VERSION_STRING} found: INC=${HDF5_INCLUDE_DIR} TOOLS=${HDF5_TOOLS_DIR} HDF5_BUILD_SHARED_LIBS=${HDF5_BUILD_SHARED_LIBS}")
      set(HDF5_FOUND TRUE)
      add_definitions(-DHAVE_HDF5)
    endif (HDF5_FOUND)    
  else (NOT HDF5_HDF5_HEADER)
    # This project is being called from within another and HDF5 is already configured
    set (HDF5_HAVE_H5PUBCONF_H 1)
    set (HDF5_HAVE_HDF5 1)
  endif (NOT HDF5_HDF5_HEADER)
  message (STATUS "HDF5 link libs: ${HDF5_LIBS}")

  option (HDF5_ENABLE_FF_SUPPORT "Enable FF version of HDF5" OFF)
  mark_as_advanced (HDF5_ENABLE_FF_SUPPORT)

  if (HDF5_ENABLE_FF_SUPPORT)
    add_definitions (-DHDF5_FF)
    message (STATUS "Enabling HDF5 FF Support")
  endif (HDF5_ENABLE_FF_SUPPORT)

  if (HDF5_FOUND)
    INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
  endif (HDF5_FOUND)

else (HDF5_ENABLE)
  message (STATUS "Disabling HDF5 Support")
endif(HDF5_ENABLE)

#-----------------------------------------------------------------------------
# Option to use Shared/Static libs, default is static
#-----------------------------------------------------------------------------
set (LIB_TYPE STATIC)
if (BUILD_SHARED_LIBS)
  set (LIB_TYPE SHARED)
endif (BUILD_SHARED_LIBS)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT INTEGER_TYPE)
    set (INTEGER_TYPE reg_integer)
endif (NOT INTEGER_TYPE)

if (INTEGER_TYPE STREQUAL "reg_integer")
    set (REG_INTEGER on)
    set (INTEGER_TYPE reg_integer)
    add_definitions(-DREG_INTEGER)
endif (INTEGER_TYPE STREQUAL "reg_integer")

if (INTEGER_TYPE STREQUAL "short_integer")
    set (SHORT_INTEGER on)
    set (INTEGER_TYPE short_integer)
    add_definitions(-DSHORT_INTEGER)
endif (INTEGER_TYPE STREQUAL "short_integer")

if (INTEGER_TYPE STREQUAL "min_integer")
    set (MIN_INTEGER on)
    set (INTEGER_TYPE min_integer)
    add_definitions(-DMIN_INTEGER)
endif (INTEGER_TYPE STREQUAL "min_integer")

if (NOT PRECISION_TYPE)
   set (PRECISION_TYPE full_precision)
endif (NOT PRECISION_TYPE)

if (PRECISION_TYPE STREQUAL "full_precision")
   set (FULL_PRECISION on)
   set (PRECISION_TYPE full_precision)
   add_definitions(-DFULL_PRECISION)
endif (PRECISION_TYPE STREQUAL "full_precision")

if (PRECISION_TYPE STREQUAL "mixed_precision")
   set (MIXED_PRECISION on)
   set (PRECISION_TYPE mixed_precision)
   add_definitions(-DMIXED_PRECISION)
endif (PRECISION_TYPE STREQUAL "mixed_precision")

if (PRECISION_TYPE STREQUAL "minimum_precision")
   set (MINIMUM_PRECISION on)
   set (PRECISION_TYPE minimum_precision)
   add_definitions(-DMINIMUM_PRECISION)
endif (PRECISION_TYPE STREQUAL "minimum_precision")

if (PRECISION_TYPE STREQUAL "half_precision")
   set (HALF_PRECISION on)
   set (PRECISION_TYPE half_precision)
   add_definitions(-DHALF_PRECISION)
   include_directories(${CMAKE_SOURCE_DIR}/half/include)
endif (PRECISION_TYPE STREQUAL "half_precision")

if (NOT GRAPHICS_TYPE)
   set (GRAPHICS_TYPE none)
endif (NOT GRAPHICS_TYPE)

if (GRAPHICS_TYPE STREQUAL "OpenGL")
   find_package(OpenGL)
   if (OPENGL_FOUND)
      find_package(GLUT)
   endif (OPENGL_FOUND)
endif (GRAPHICS_TYPE STREQUAL "OpenGL")

# If autodetect does not work try setting MPEHOME or
# MPE_INCLUDE_DIR and MPE_LIBRARIES
if (GRAPHICS_TYPE STREQUAL "MPE")
   find_package(MPE)
   find_package(X11)
endif (GRAPHICS_TYPE STREQUAL "MPE")

#message("MPI_LIBRARIES ${MPI_LIBRARIES}")
#message("MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH}")
#message("GRAPHICS_TYPE is ${GRAPHICS_TYPE}")
#message("MPE_INCLUDE is ${MPE_INCLUDE_DIR}")

if (OPENGL_FOUND AND GLUT_FOUND)
   set (HAVE_GRAPHICS on)
   set (HAVE_OPENGL on)
   include_directories(${OPENGL_INCLUDE_DIR})
   include_directories(${GLUT_INCLUDE_DIR})
endif (OPENGL_FOUND AND GLUT_FOUND)

#if(DEFINED ENV{QUO_HOME})
#   if(EXISTS "$ENV{QUO_HOME}/include/")
#      message(STATUS "QUO Detected")
#      set (QUO_FOUND on)
#      include_directories($ENV{QUO_HOME}/include)
#   endif (EXISTS "$ENV{QUO_HOME}/include/")
#endif (DEFINED ENV{QUO_HOME})

#if(DEFINED ENV{REPROBLAS_HOME})
#   if(EXISTS "$ENV{REPROBLAS_HOME}/include/rblas.h")
#      message(STATUS "reproBlas Detected")
#      set (REPROBLAS_FOUND on)
#      set (HAVE_REPROBLAS ON)
#      include_directories($ENV{REPROBLAS_HOME}/include)
#      set (REPROBLAS_LIBRARIES -L$ENV{REPROBLAS_HOME}/libs -lreproblas_mpi -lreproblas_seq -liblas -lmpindexedfp -lindexedfp)
#   endif (EXISTS "$ENV{REPROBLAS_HOME}/include/rblas.h")
#endif (DEFINED ENV{REPROBLAS_HOME})

#if(DEFINED ENV{LTTRACE_HOME})
#   if(EXISTS "$ENV{LTTRACE_HOME}/lib/liblttrace.a")
#      message(STATUS "LTTRACE Detected")
#      include_directories($ENV{LTTRACE_HOME}/include)
#      set (LTTRACE_FOUND on)
#      set (HAVE_LTTRACE ON)
#      set (LTTRACE_LIBRARIES -L$ENV{LTTRACE_HOME}/lib -llttrace -lz -luuid -Wl,-rpath -Wl,$ENV{LTTRACE_HOME}/lib)
#   endif (EXISTS "$ENV{LTTRACE_HOME}/lib/liblttrace.a")
#endif (DEFINED ENV{LTTRACE_HOME})

#if (DEFINED ENV{LD_LIBRARY_PATH})
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,$ENV{LD_LIBRARY_PATH}")
#   message("CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}")
#   #string (REPLACE ":" ";" LDPATH_TMP $ENV{LD_LIBRARY_PATH})
#   #foreach(LDPATH ${LDPATH_TMP})
#   #   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,${LDPATH}")
#   #   message("CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}")
#   #endforeach(LDPATH)
#endif (DEFINED ENV{LD_LIBRARY_PATH})

if (NOT OPENGL_FOUND)
   if (MPE_FOUND)
      set (HAVE_GRAPHICS on)
      set (HAVE_MPE on)
      include_directories(${MPE_INCLUDE_DIR})
      if (X11_FOUND)
         include_directories(${X11_INCLUDE_DIR})
      endif (X11_FOUND)
      set (MPE_LIBS ${MPE_LIBRARIES})
      set (MPE_NOMPI_LIBS ${MPE_NOMPI_LIBRARIES})
      set (X11_LIBS ${X11_LIBRARIES})
   endif (MPE_FOUND)
endif (NOT OPENGL_FOUND)

# search for Boost
# this min version isn't correct -- FIXME
#set (BOOST_MIN_VER "1.45.0")
#find_package(Boost ${BOOST_MIN_VER})
 
find_library(POWER_GADGET IntelPowerGadget)
if (POWER_GADGET)
   set (POWER_GADGET_FOUND on)
   add_definitions(-DHAVE_POWER_GADGET)
   message("-- Detected Intel Power Gadget")
   include_directories(/Library/Frameworks/IntelPowerGadget.framework/Headers)
endif (POWER_GADGET)


#message("MPE_LIBRARIES ${MPE_LIBRARIES}")
#message("MPE_LIBS ${MPE_LIBS}")
#message("MPE_NOMPI_LIBS ${MPE_NOMPI_LIBS}")
#message("OPENGL_LIBRARIES ${OPENGL_LIBRARIES}")
#message("GLUT_LIBRARIES ${GLUT_LIBRARIES}")
#message("X11_LIBRARIES ${X11_LIBRARIES}")
#message("X11_LIBS ${X11_LIBS}")

add_definitions(-DHAVE_CONFIG_H)

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")

#ADDING TESTS CODE####

if (EXISTS "${CMAKE_SOURCE_DIR}/tests" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/tests")
   add_subdirectory(tests)
endif (EXISTS "${CMAKE_SOURCE_DIR}/tests" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/tests")

if (EXISTS "${CMAKE_SOURCE_DIR}/mesh" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/mesh")
   add_subdirectory(mesh)
endif (EXISTS "${CMAKE_SOURCE_DIR}/mesh" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/mesh")

if (EXISTS "${CMAKE_SOURCE_DIR}/s7" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/s7")
   add_subdirectory(s7)
endif (EXISTS "${CMAKE_SOURCE_DIR}/s7" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/s7")

if (EXISTS "${CMAKE_SOURCE_DIR}/l7" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/l7")
   add_subdirectory(l7)
endif (EXISTS "${CMAKE_SOURCE_DIR}/l7" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/l7")

if (EXISTS "${CMAKE_SOURCE_DIR}/graphics" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/graphics")
   add_subdirectory(graphics)
endif (EXISTS "${CMAKE_SOURCE_DIR}/graphics" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/graphics")

if (EXISTS "${CMAKE_SOURCE_DIR}/timer" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/timer")
   add_subdirectory(timer)
endif (EXISTS "${CMAKE_SOURCE_DIR}/timer" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/timer")

if (EXISTS "${CMAKE_SOURCE_DIR}/memstats" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/memstats")
   add_subdirectory(memstats)
endif (EXISTS "${CMAKE_SOURCE_DIR}/memstats" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/memstats")

if (EXISTS "${CMAKE_SOURCE_DIR}/powerstats" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/powerstats")
   add_subdirectory(powerstats)
endif (EXISTS "${CMAKE_SOURCE_DIR}/powerstats" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/powerstats")

if (OpenCL_FOUND)
   if (EXISTS "${CMAKE_SOURCE_DIR}/ezcl" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/ezcl")
      add_subdirectory(ezcl)
   endif (EXISTS "${CMAKE_SOURCE_DIR}/ezcl" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/ezcl")
endif (OpenCL_FOUND)

if (EXISTS "${CMAKE_SOURCE_DIR}/genmalloc" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/genmalloc")
   add_subdirectory(genmalloc)
endif (EXISTS "${CMAKE_SOURCE_DIR}/genmalloc" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/genmalloc")

if (EXISTS "${CMAKE_SOURCE_DIR}/MallocPlus" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/MallocPlus")
   add_subdirectory(MallocPlus)
endif (EXISTS "${CMAKE_SOURCE_DIR}/MallocPlus" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/MallocPlus")

if (EXISTS "${CMAKE_SOURCE_DIR}/hash" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/hash")
   add_subdirectory(hash)
endif (EXISTS "${CMAKE_SOURCE_DIR}/hash" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/hash")

if (EXISTS "${CMAKE_SOURCE_DIR}/crux" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/crux")
   add_subdirectory(crux)
endif (EXISTS "${CMAKE_SOURCE_DIR}/crux" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/crux")

if (EXISTS "${CMAKE_SOURCE_DIR}/PowerParser" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/PowerParser")
   add_subdirectory(PowerParser)
endif (EXISTS "${CMAKE_SOURCE_DIR}/PowerParser" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/PowerParser")

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/state_kern.cl > state_kernel.inc
                  DEPENDS state_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(state_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc)

########### clamr target ###############

if (OpenCL_FOUND AND MPI_FOUND)
   set(clamr_CXX_SRCS clamr.cpp state.cpp input.cpp)
   set(clamr_H_SRCS state.h input.h)
   set(clamr_SRCS ${clamr_CXX_SRCS} ${clamr_C_SRCS} ${clamr_H_SRCS})

   add_executable(clamr ${clamr_SRCS})

   set_target_properties(clamr PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(clamr dpmesh hsfc kdtree zorder pcrux dhash pezcl pgraphics timer memstats powerstats dl7 genmalloc dpMallocPlus m)
   target_link_libraries(clamr ${OpenCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   target_link_libraries(clamr ${MPI_LIBRARIES})
   add_dependencies(clamr state_kernel_source)

   install(TARGETS clamr DESTINATION bin)
endif (OpenCL_FOUND AND MPI_FOUND)

########### clamr_kokkosonly target ###############

if(KOKKOS_ENABLE)

    #if (EXISTS "${CMAKE_SOURCE_DIR}/kokkos" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/kokkos")
      set(clamr_kokkosonly_CXX_SRCS clamr_kokkosonly.cpp state.cpp input.cpp)
      set(clamr_kokkosonly_H_SRCS state.h input.h)

      set(clamr_kokkosonly_SRCS ${clamr_kokkosonly_CXX_SRCS} ${clamr_kokkosonly_C_SRCS} ${clamr_kokkosonly_H_SRCS})

      add_executable(clamr_kokkosonly ${clamr_kokkosonly_SRCS})

      add_subdirectory(kokkos)

      target_link_libraries(clamr_kokkosonly Kokkos::kokkos)

      target_link_libraries(clamr_kokkosonly PowerParser mesh hsfc kdtree zorder crux hash graphics timer memstats powerstats genmalloc MallocPlus m)
      if (REPROBLAS_FOUND)
         target_link_libraries(clamr_kokkosonly ${REPROBLAS_LIBRARIES})
      endif (REPROBLAS_FOUND)
      add_dependencies(clamr_kokkosonly state_kernel_source)

      if (POWER_GADGET_FOUND)
         target_link_libraries(clamr_kokkosonly ${POWER_GADGET})
      endif (POWER_GADGET_FOUND)

      #target_link_libraries(clamr_cpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

      install(TARGETS clamr_kokkosonly DESTINATION bin)
      #endif (EXISTS "${CMAKE_SOURCE_DIR}/kokkos" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/kokkos")
endif(KOKKOS_ENABLE)

########### drawit target ###############

add_executable(drawit drawit.cpp)

target_link_libraries(drawit graphics)

target_link_libraries(drawit ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS drawit DESTINATION bin)

########### clamr_cpuonly target ###############

set(clamr_cpuonly_CXX_SRCS clamr_cpuonly.cpp state.cpp input.cpp)
set(clamr_cpuonly_H_SRCS state.h input.h)

set(clamr_cpuonly_SRCS ${clamr_cpuonly_CXX_SRCS} ${clamr_cpuonly_C_SRCS} ${clamr_cpuonly_H_SRCS})

add_executable(clamr_cpuonly ${clamr_cpuonly_SRCS})

target_link_libraries(clamr_cpuonly PowerParser mesh hsfc kdtree zorder crux hash graphics timer memstats powerstats genmalloc MallocPlus m)
if (REPROBLAS_FOUND)
   target_link_libraries(clamr_cpuonly ${REPROBLAS_LIBRARIES})
endif (REPROBLAS_FOUND)
   add_dependencies(clamr_cpuonly state_kernel_source)

if (POWER_GADGET_FOUND)
   target_link_libraries(clamr_cpuonly ${POWER_GADGET})
endif (POWER_GADGET_FOUND)

#target_link_libraries(clamr_cpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS clamr_cpuonly DESTINATION bin)

########### clamr_gpuonly target ###############
if (OpenCL_FOUND)
   set(clamr_gpuonly_CXX_SRCS clamr_gpuonly.cpp state.cpp input.cpp)
   set(clamr_gpuonly_H_SRCS state.h input.h)
   set(clamr_gpuonly_SRCS ${clamr_gpuonly_CXX_SRCS} ${clamr_gpuonly_C_SRCS} ${clamr_gpuonly_H_SRCS})

   add_executable(clamr_gpuonly ${clamr_gpuonly_SRCS})

   set_target_properties(clamr_gpuonly PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

   target_link_libraries(clamr_gpuonly PowerParser dmesh hsfc kdtree zorder crux dhash ezcl graphics timer memstats powerstats genmalloc dMallocPlus m)
   target_link_libraries(clamr_gpuonly ${OpenCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_gpuonly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   add_dependencies(clamr_gpuonly state_kernel_source)

   install(TARGETS clamr_gpuonly DESTINATION bin)
endif (OpenCL_FOUND)

########### clamr_mpionly target ###############
if (MPI_FOUND)
   set(clamr_mpionly_CXX_SRCS clamr_mpionly.cpp state.cpp input.cpp)
   set(clamr_mpionly_H_SRCS state.h input.h)

   set(clamr_mpionly_SRCS ${clamr_mpionly_CXX_SRCS} ${clamr_mpionly_C_SRCS} ${clamr_mpionly_H_SRCS})

   add_executable(clamr_mpionly ${clamr_mpionly_SRCS})

   set_target_properties(clamr_mpionly PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(clamr_mpionly pPowerParser pmesh hsfc kdtree zorder pcrux hash pgraphics timer memstats powerstats l7 genmalloc pMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_mpionly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_mpionly ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)

   target_link_libraries(clamr_mpionly ${MPI_LIBRARIES})

   install(TARGETS clamr_mpionly DESTINATION bin)
endif (MPI_FOUND)

########### clamr_openmponly target ###############
# This shares the clamr_cpuonly.cpp source file
if(OPENMP_FOUND)
   set(clamr_openmponly_CXX_SRCS clamr_cpuonly.cpp state.cpp input.cpp)
   set(clamr_openmponly_H_SRCS state.h input.h)

   set(clamr_openmponly_SRCS ${clamr_openmponly_CXX_SRCS} ${clamr_openmponly_C_SRCS} ${clamr_openmponly_H_SRCS})

   add_executable(clamr_openmponly ${clamr_openmponly_SRCS})

   set_target_properties(clamr_openmponly PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(clamr_openmponly PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(clamr_openmponly PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")

   target_link_libraries(clamr_openmponly PowerParser mpmesh hsfc kdtree zorder crux mphash graphics mptimer memstats powerstats genmalloc mpMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_openmponly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)

   if (POWER_GADGET_FOUND)
      target_link_libraries(clamr_openmponly ${POWER_GADGET})
   endif (POWER_GADGET_FOUND)

   install(TARGETS clamr_openmponly DESTINATION bin)
endif(OPENMP_FOUND)

########### clamr_mpiopenmponly target ###############
# This uses the clamr_mpionly driver with the OpenMP compile flag set
if (MPI_FOUND AND OPENMP_FOUND)
   set(clamr_mpiopenmponly_CXX_SRCS clamr_mpionly.cpp state.cpp input.cpp)
   set(clamr_mpiopenmponly_H_SRCS state.h input.h)

   set(clamr_mpiopenmponly_SRCS ${clamr_mpiopenmponly_CXX_SRCS} ${clamr_mpiopenmponly_C_SRCS} ${clamr_mpiopenmponly_H_SRCS})

   add_executable(clamr_mpiopenmponly ${clamr_mpiopenmponly_SRCS})

   set_target_properties(clamr_mpiopenmponly PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI")

   set(clamr_mpiopenmponly_link_flags "${OpenMP_CXX_FLAGS}")
   set_target_properties(clamr_mpiopenmponly PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(clamr_mpiopenmponly PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(clamr_mpiopenmponly PROPERTIES LINK_FLAGS "${clamr_mpiopenmponly_link_flags}")

   target_link_libraries(clamr_mpiopenmponly pPowerParser pmpmesh hsfc kdtree zorder pcrux mphash pgraphics mptimer memstats powerstats mpl7 genmalloc pmpMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_mpiopenmponly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_mpiopenmponly ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   target_link_libraries(clamr_mpiopenmponly ${MPI_LIBRARIES})

   install(TARGETS clamr_mpiopenmponly DESTINATION bin)

endif (MPI_FOUND AND OPENMP_FOUND)

########### clamr_quo target ###############
# only build this thing if we have QUO, MPI, and OpenMP

if (QUO_FOUND AND MPI_FOUND AND OPENMP_FOUND AND Boost_FOUND)
    # doing this way because we only want j7 if doing quo (for now)
    # this only enables experimental software if the user has all the things
    message(STATUS "Adding J7 Support")
    add_subdirectory(j7)
    message(STATUS "clamr_quo Minimum Requires Met")

    set(clamr_quo_CXX_SRCS clamr_quo.cpp state.cpp input.cpp)
    set(clamr_quo_H_SRCS state.h input.h)

    set(clamr_quo_SRCS ${clamr_quo_CXX_SRCS} ${clamr_quo_C_SRCS} ${clamr_quo_H_SRCS})

    add_executable(clamr_quo ${clamr_quo_SRCS})

    # these sets HAVE_J7 and HAVE_QUO in config.h
    set(HAVE_J7 ON)
    set(HAVE_QUO ON)
    set_target_properties(clamr_quo PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI")

    set(clamr_quo_link_flags "${OpenMP_CXX_FLAGS}")
    set_target_properties(clamr_quo PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
    set_target_properties(clamr_quo PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

    set(clamr_quo_link_flags "${clamr_quo_link_flags} -L$ENV{QUO_HOME}/lib")

    set_target_properties(clamr_quo PROPERTIES LINK_FLAGS "${clamr_quo_link_flags}")

    target_link_libraries(clamr_quo pmesh hsfc kdtree zorder crux hash pgraphics timer memstats powerstats l7 genmalloc pMallocPlus m)
    target_link_libraries(clamr_quo quo)
    target_link_libraries(clamr_quo j7)
    if (REPROBLAS_FOUND)
       target_link_libraries(clamr_quo ${REPROBLAS_LIBRARIES})
    endif (REPROBLAS_FOUND)
    if (LTTRACE_FOUND)
       target_link_libraries(clamr_quo ${LTTRACE_LIBRARIES})
    endif (LTTRACE_FOUND)
    target_link_libraries(clamr_quo ${MPI_LIBRARIES})

    install(TARGETS clamr_quo DESTINATION bin)

    # XXX add load imbalance code here for now
    set(load_imbl_CXX_SRCS QUO.cpp QUO.hpp load-imbl.cpp)
    add_executable(load-imbl ${load_imbl_CXX_SRCS})

    set_target_properties(load-imbl PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    set_target_properties(load-imbl PROPERTIES LINK_FLAGS
                          "${OpenMP_CXX_FLAGS} -L$ENV{QUO_HOME}/lib")
    target_link_libraries(load-imbl ${MPI_LIBRARIES})
    target_link_libraries(load-imbl j7 quo)
endif()

########### clamr_gpucheck target ###############
if (OpenCL_FOUND)
   set(clamr_gpucheck_CXX_SRCS clamr_gpucheck.cpp state.cpp input.cpp)
   set(clamr_gpucheck_H_SRCS state.h input.h)

   set(clamr_gpucheck_SRCS ${clamr_gpucheck_CXX_SRCS} ${clamr_gpucheck_C_SRCS} ${clamr_gpucheck_H_SRCS})

   add_executable(clamr_gpucheck ${clamr_gpucheck_SRCS})

   set_target_properties(clamr_gpucheck PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

   target_link_libraries(clamr_gpucheck PowerParser dmesh hsfc kdtree zorder crux dhash ezcl graphics timer memstats powerstats genmalloc dMallocPlus m)
   target_link_libraries(clamr_gpucheck ${OpenCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_gpucheck ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   add_dependencies(clamr_gpucheck state_kernel_source)

   install(TARGETS clamr_gpucheck DESTINATION bin)
endif (OpenCL_FOUND)

########### clamr_mpicheck target ###############
if (MPI_FOUND)
   set(clamr_mpicheck_CXX_SRCS clamr_mpicheck.cpp state.cpp input.cpp)
   set(clamr_mpicheck_H_SRCS state.h input.h)

   set(clamr_mpicheck_SRCS ${clamr_mpicheck_CXX_SRCS} ${clamr_mpicheck_C_SRCS} ${clamr_mpicheck_H_SRCS})

   add_executable(clamr_mpicheck ${clamr_mpicheck_SRCS})

   set_target_properties(clamr_mpicheck PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(clamr_mpicheck pPowerParser pmesh hsfc kdtree zorder pcrux hash pgraphics timer memstats powerstats l7 genmalloc pMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_mpicheck ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_mpicheck ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   target_link_libraries(clamr_mpicheck ${MPI_LIBRARIES})

   install(TARGETS clamr_mpicheck DESTINATION bin)
endif (MPI_FOUND)

########### clamr_checkall target ###############
if (MPI_FOUND AND OpenCL_FOUND)
   set(clamr_checkall_CXX_SRCS clamr_checkall.cpp state.cpp input.cpp)
   set(clamr_checkall_H_SRCS state.h input.h)

   set(clamr_checkall_SRCS ${clamr_checkall_CXX_SRCS} ${clamr_checkall_C_SRCS} ${clamr_checkall_H_SRCS})

   add_executable(clamr_checkall ${clamr_checkall_SRCS})

   set_target_properties(clamr_checkall PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(clamr_checkall dpmesh hsfc kdtree zorder pcrux dhash pezcl pgraphics timer memstats powerstats dl7 genmalloc dpMallocPlus m)
   target_link_libraries(clamr_checkall ${OpenCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_checkall ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_checkall ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   target_link_libraries(clamr_checkall ${MPI_LIBRARIES})
   add_dependencies(clamr_checkall state_kernel_source)

   install(TARGETS clamr_checkall DESTINATION bin)
endif (MPI_FOUND AND OpenCL_FOUND)

########### install files ###############


if (OpenCL_FOUND)
   set (CMAKE_CHECK_COMMAND "make ezcl_check &&")
else (OpenCL_FOUND)
   set (CMAKE_CHECK_COMMAND "")
endif (OpenCL_FOUND)

set (CMAKE_CHECK_COMMAND ${CMAKE_CHECK_COMMAND}
                         make mesh_check &&
                         make l7_check &&
                         make MallocPlus_check &&
                         make PowerParser_check)

add_custom_target(clamr_check COMMAND ${CMAKE_CHECK_COMMAND})

set (CMAKE_DOC_COMMAND #make CLAMR_doc       &&
                       make PowerParser_doc &&
                       make MallocPlus_doc  &&
                       make L7_doc          &&
                       make ezcl_doc        &&
                       make genmalloc_doc   &&
                       make Graphics_doc    &&
                       make Hash_doc        &&
                       make memstats_doc    &&
                       make powerstats_doc  &&
                       make Mesh_doc        &&
                       make Crux_doc)

add_custom_target(doc COMMAND ${CMAKE_DOC_COMMAND})

########### doc targets ################

#add_custom_target(CLAMR_doc
#                  DEPENDS CLAMR_pdfdoc CLAMR_htmldoc CLAMR_singlehtmldoc)
#set_target_properties(CLAMR_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)
#
#add_custom_command(OUTPUT _build/latex/CLAMR.pdf
#                   COMMAND make -C docs latexpdf >& pdfdoc.out
#                   DEPENDS ${PDFDOC_SRCS})
#add_custom_target(CLAMR_pdfdoc DEPENDS _build/latex/CLAMR.pdf)
#set_target_properties(CLAMR_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
#
#add_custom_command(OUTPUT _build/html/UsersGuide.html
#                   COMMAND make -C docs html >& htmldoc.out
#                   DEPENDS ${DOC_SRCS})
#add_custom_target(CLAMR_htmldoc DEPENDS _build/html/UsersGuide.html)
#set_target_properties(CLAMR_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
#
#add_custom_command(OUTPUT _build/singlehtml/index.html
#                   COMMAND make -C docs singlehtml >& singlehtmldoc.out
#                   DEPENDS ${DOC_SRCS})
#add_custom_target(CLAMR_singlehtmldoc DEPENDS _build/singlehtml/index.html)
#set_target_properties(CLAMR_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

########### configure config.h ################
# configure a header file to pass some of the CMake settings
# to the source code
# at bottom because we set things that change config above
configure_file ( "${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "state_kernel.inc;clamr;clamr_cpuonly;clamr_gpuonly;clamr_mpionly;clamr_gpucheck;clamr_mpicheck;clamr_checkall;clamr_quo;total_execution_time.log;docs/_build;htmldoc.out;pdfdoc.out;singlehtmldoc.out")

