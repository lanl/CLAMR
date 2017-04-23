cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (CLAMR)
# The version number.
set (CLAMR_VERSION_MAJOR 2)
set (CLAMR_VERSION_MINOR 0.7)
set (CLAMR_VERSION "${CLAMR_VERSION_MAJOR}.${CLAMR_VERSION_MINOR}")

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

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -std=c99 -mmic")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
else (MIC_NATIVE)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif (MIC_NATIVE)

if (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native -Wa,-q")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Wa,-q")
else (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
endif (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-alias -restrict -opt_report")

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wc++11-long-long -Wno-overlength-strings")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pedantic")

include_directories(${PROJECT_SOURCE_DIR})

if (NOT MIC_NATIVE)
   find_package(OpenCL)
   #  OPENCL_FOUND        : TRUE if found, FALSE otherwise
   #  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
   #  OPENCL_LIBRARIES    : The libraries to link against
   if (OPENCL_FOUND)
      include_directories(${OPENCL_INCLUDE_DIRS})
   endif (OPENCL_FOUND)
endif (NOT MIC_NATIVE)

find_package(MPI)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

find_package( OpenMP)
#if(OPENMP_FOUND)
#   message("OPENMP FOUND ${OpenMP_C_FLAGS}")
#else (OPENMP_FOUND)
#   message("OPENMP NOT FOUND")
#endif(OPENMP_FOUND)

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

if (NOT GRAPHICS_TYPE)
   set (GRAPHICS_TYPE none)
endif (NOT GRAPHICS_TYPE)

if (GRAPHICS_TYPE STREQUAL "OpenGL")
   find_package(OpenGL)
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

if (OPENGL_FOUND)
   set (HAVE_GRAPHICS on)
   set (HAVE_OPENGL on)
   include_directories(${OPENGL_INCLUDE_DIR})
   find_package(GLUT)
endif (OPENGL_FOUND)

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
set (BOOST_MIN_VER "1.45.0")
find_package(Boost ${BOOST_MIN_VER})

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

if (EXISTS "mesh" AND IS_DIRECTORY "mesh")
   add_subdirectory(mesh)
endif (EXISTS "mesh" AND IS_DIRECTORY "mesh")

if (EXISTS "s7" AND IS_DIRECTORY "s7")
   add_subdirectory(s7)
endif (EXISTS "s7" AND IS_DIRECTORY "s7")

if (EXISTS "l7" AND IS_DIRECTORY "l7")
   add_subdirectory(l7)
endif (EXISTS "l7" AND IS_DIRECTORY "l7")

if (EXISTS "graphics" AND IS_DIRECTORY "graphics")
   add_subdirectory(graphics)
endif (EXISTS "graphics" AND IS_DIRECTORY "graphics")

if (EXISTS "timer" AND IS_DIRECTORY "timer")
   add_subdirectory(timer)
endif (EXISTS "timer" AND IS_DIRECTORY "timer")

if (EXISTS "memstats" AND IS_DIRECTORY "memstats")
   add_subdirectory(memstats)
endif (EXISTS "memstats" AND IS_DIRECTORY "memstats")

if (OPENCL_FOUND)
   if (EXISTS "ezcl" AND IS_DIRECTORY "ezcl")
      add_subdirectory(ezcl)
   endif (EXISTS "ezcl" AND IS_DIRECTORY "ezcl")
endif (OPENCL_FOUND)

if (EXISTS "genmalloc" AND IS_DIRECTORY "genmalloc")
   add_subdirectory(genmalloc)
endif (EXISTS "genmalloc" AND IS_DIRECTORY "genmalloc")

if (EXISTS "MallocPlus" AND IS_DIRECTORY "MallocPlus")
   add_subdirectory(MallocPlus)
endif (EXISTS "MallocPlus" AND IS_DIRECTORY "MallocPlus")

if (EXISTS "hash" AND IS_DIRECTORY "hash")
   add_subdirectory(hash)
endif (EXISTS "hash" AND IS_DIRECTORY "hash")

if (EXISTS "crux" AND IS_DIRECTORY "crux")
   add_subdirectory(crux)
endif (EXISTS "crux" AND IS_DIRECTORY "crux")

if (EXISTS "PowerParser" AND IS_DIRECTORY "PowerParser")
   add_subdirectory(PowerParser)
endif (EXISTS "PowerParser" AND IS_DIRECTORY "PowerParser")

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/state_kern.cl > state_kernel.inc
                  DEPENDS state_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(state_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc)

########### clamr target ###############

if (OPENCL_FOUND AND MPI_FOUND)
   set(clamr_CXX_SRCS clamr.cpp state.cpp input.cpp)
   set(clamr_H_SRCS state.h input.h)
   set(clamr_SRCS ${clamr_CXX_SRCS} ${clamr_C_SRCS} ${clamr_H_SRCS})

   add_executable(clamr ${clamr_SRCS})

   set_target_properties(clamr PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(clamr dpmesh hsfc kdtree zorder crux dhash pezcl pgraphics timer memstats dl7 genmalloc dpMallocPlus m)
   target_link_libraries(clamr ${OPENCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr iomp5 cilkrts)
   endif (MIC_NATIVE)
   target_link_libraries(clamr ${MPI_LIBRARIES})
   add_dependencies(clamr state_kernel_source)

   install(TARGETS clamr DESTINATION bin)
endif (OPENCL_FOUND AND MPI_FOUND)

########### clamr_cpuonly target ###############

set(clamr_cpuonly_CXX_SRCS clamr_cpuonly.cpp state.cpp input.cpp)
set(clamr_cpuonly_H_SRCS state.h input.h)

set(clamr_cpuonly_SRCS ${clamr_cpuonly_CXX_SRCS} ${clamr_cpuonly_C_SRCS} ${clamr_cpuonly_H_SRCS})

add_executable(clamr_cpuonly ${clamr_cpuonly_SRCS})

target_link_libraries(clamr_cpuonly PowerParser mesh hsfc kdtree zorder crux hash graphics timer memstats genmalloc MallocPlus m)
if (REPROBLAS_FOUND)
   target_link_libraries(clamr_cpuonly ${REPROBLAS_LIBRARIES})
endif (REPROBLAS_FOUND)
#target_link_libraries(clamr_cpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS clamr_cpuonly DESTINATION bin)

########### clamr_gpuonly target ###############
if (OPENCL_FOUND)
   set(clamr_gpuonly_CXX_SRCS clamr_gpuonly.cpp state.cpp input.cpp)
   set(clamr_gpuonly_H_SRCS state.h input.h)
   set(clamr_gpuonly_SRCS ${clamr_gpuonly_CXX_SRCS} ${clamr_gpuonly_C_SRCS} ${clamr_gpuonly_H_SRCS})

   add_executable(clamr_gpuonly ${clamr_gpuonly_SRCS})

   set_target_properties(clamr_gpuonly PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

   target_link_libraries(clamr_gpuonly PowerParser dmesh hsfc kdtree zorder crux dhash ezcl graphics timer memstats genmalloc dMallocPlus m)
   target_link_libraries(clamr_gpuonly ${OPENCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_gpuonly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_gpuonly iomp5 cilkrts)
   endif (MIC_NATIVE)
   add_dependencies(clamr_gpuonly state_kernel_source)

   install(TARGETS clamr_gpuonly DESTINATION bin)
endif (OPENCL_FOUND)

########### clamr_mpionly target ###############
if (MPI_FOUND)
   set(clamr_mpionly_CXX_SRCS clamr_mpionly.cpp state.cpp input.cpp)
   set(clamr_mpionly_H_SRCS state.h input.h)

   set(clamr_mpionly_SRCS ${clamr_mpionly_CXX_SRCS} ${clamr_mpionly_C_SRCS} ${clamr_mpionly_H_SRCS})

   add_executable(clamr_mpionly ${clamr_mpionly_SRCS})

   set_target_properties(clamr_mpionly PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(clamr_mpionly pPowerParser pmesh hsfc kdtree zorder crux hash pgraphics timer memstats l7 genmalloc pMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_mpionly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_mpionly ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_mpionly iomp5 cilkrts)
   endif (MIC_NATIVE)
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

   target_link_libraries(clamr_openmponly PowerParser mpmesh hsfc kdtree zorder crux mphash graphics timer memstats genmalloc mpMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_openmponly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_openmponly iomp5 cilkrts)
   endif (MIC_NATIVE)

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

   target_link_libraries(clamr_mpiopenmponly pPowerParser pmpmesh hsfc kdtree zorder crux mphash pgraphics timer memstats mpl7 genmalloc pmpMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_mpiopenmponly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_mpiopenmponly ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_mpiopenmponly iomp5 cilkrts)
   endif (MIC_NATIVE)
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

    target_link_libraries(clamr_quo pmesh hsfc kdtree zorder crux hash pgraphics timer memstats l7 genmalloc pMallocPlus m)
    target_link_libraries(clamr_quo quo)
    target_link_libraries(clamr_quo j7)
    if (REPROBLAS_FOUND)
       target_link_libraries(clamr_quo ${REPROBLAS_LIBRARIES})
    endif (REPROBLAS_FOUND)
    if (LTTRACE_FOUND)
       target_link_libraries(clamr_quo ${LTTRACE_LIBRARIES})
    endif (LTTRACE_FOUND)
    if (MIC_NATIVE)
       target_link_libraries(clamr_quo iomp5 cilkrts)
    endif (MIC_NATIVE)
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
if (OPENCL_FOUND)
   set(clamr_gpucheck_CXX_SRCS clamr_gpucheck.cpp state.cpp input.cpp)
   set(clamr_gpucheck_H_SRCS state.h input.h)

   set(clamr_gpucheck_SRCS ${clamr_gpucheck_CXX_SRCS} ${clamr_gpucheck_C_SRCS} ${clamr_gpucheck_H_SRCS})

   add_executable(clamr_gpucheck ${clamr_gpucheck_SRCS})

   set_target_properties(clamr_gpucheck PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

   target_link_libraries(clamr_gpucheck PowerParser dmesh hsfc kdtree zorder crux dhash ezcl graphics timer memstats genmalloc dMallocPlus m)
   target_link_libraries(clamr_gpucheck ${OPENCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_gpucheck ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_gpucheck iomp5 cilkrts)
   endif (MIC_NATIVE)
   add_dependencies(clamr_gpucheck state_kernel_source)

   install(TARGETS clamr_gpucheck DESTINATION bin)
endif (OPENCL_FOUND)

########### clamr_mpicheck target ###############
if (MPI_FOUND)
   set(clamr_mpicheck_CXX_SRCS clamr_mpicheck.cpp state.cpp input.cpp)
   set(clamr_mpicheck_H_SRCS state.h input.h)

   set(clamr_mpicheck_SRCS ${clamr_mpicheck_CXX_SRCS} ${clamr_mpicheck_C_SRCS} ${clamr_mpicheck_H_SRCS})

   add_executable(clamr_mpicheck ${clamr_mpicheck_SRCS})

   set_target_properties(clamr_mpicheck PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(clamr_mpicheck pmesh hsfc kdtree zorder crux hash pgraphics timer memstats l7 genmalloc pMallocPlus m)
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_mpicheck ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_mpicheck ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_mpicheck iomp5 cilkrts)
   endif (MIC_NATIVE)
   target_link_libraries(clamr_mpicheck ${MPI_LIBRARIES})

   install(TARGETS clamr_mpicheck DESTINATION bin)
endif (MPI_FOUND)

########### clamr_checkall target ###############
if (MPI_FOUND AND OPENCL_FOUND)
   set(clamr_checkall_CXX_SRCS clamr_checkall.cpp state.cpp input.cpp)
   set(clamr_checkall_H_SRCS state.h input.h)

   set(clamr_checkall_SRCS ${clamr_checkall_CXX_SRCS} ${clamr_checkall_C_SRCS} ${clamr_checkall_H_SRCS})

   add_executable(clamr_checkall ${clamr_checkall_SRCS})

   set_target_properties(clamr_checkall PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(clamr_checkall dpmesh hsfc kdtree zorder crux dhash pezcl pgraphics timer memstats dl7 genmalloc dpMallocPlus m)
   target_link_libraries(clamr_checkall ${OPENCL_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_checkall ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_checkall ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   if (MIC_NATIVE)
      target_link_libraries(clamr_checkall iomp5 cilkrts)
   endif (MIC_NATIVE)
   target_link_libraries(clamr_checkall ${MPI_LIBRARIES})
   add_dependencies(clamr_checkall state_kernel_source)

   install(TARGETS clamr_checkall DESTINATION bin)
endif (MPI_FOUND AND OPENCL_FOUND)

########### install files ###############

set (CMAKE_CHECK_COMMAND make mesh_check &&
                         make l7_check &&
                         make ezcl_check &&
                         make MallocPlus_check &&
                         make PowerParser_check)

add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND})

set (CMAKE_DOC_COMMAND #make CLAMR_doc       &&
                       make PowerParser_doc &&
                       make MallocPlus_doc  &&
                       make L7_doc          &&
                       make ezcl_doc        &&
                       make genmalloc_doc   &&
                       make Graphics_doc    &&
                       make Hash_doc        &&
                       make memstats_doc    &&
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

