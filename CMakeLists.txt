cmake_minimum_required(VERSION 2.6)

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
else (MIC_NATIVE)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif (MIC_NATIVE)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wall -Wextra")
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

find_package(MPI REQUIRED)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

find_package( OpenMP)
#if(OPENMP_FOUND)
#   message("OPENMP FOUND ${OpenMP_C_FLAGS}")
#else (OPENMP_FOUND)
#   message("OPENMP NOT FOUND")
#endif(OPENMP_FOUND)

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

if(DEFINED ENV{QUO_HOME})
   if(EXISTS "$ENV{QUO_HOME}/include/")
      message(STATUS "QUO Detected")
      set (QUO_FOUND on)
      include_directories($ENV{QUO_HOME}/include)
   endif (EXISTS "$ENV{QUO_HOME}/include/")
endif (DEFINED ENV{QUO_HOME})

if(DEFINED ENV{REPROBLAS_HOME})
   if(EXISTS "$ENV{REPROBLAS_HOME}/include/rblas.h")
      message(STATUS "reproBlas Detected")
      set (REPROBLAS_FOUND on)
      set (HAVE_REPROBLAS ON)
      include_directories($ENV{REPROBLAS_HOME}/include)
      set (REPROBLAS_LIBRARIES -L$ENV{REPROBLAS_HOME}/libs -lreproblas_mpi -lreproblas_seq -liblas -lmpindexedfp -lindexedfp)
   endif (EXISTS "$ENV{REPROBLAS_HOME}/include/rblas.h")
endif (DEFINED ENV{REPROBLAS_HOME})

if(DEFINED ENV{LTTRACE_HOME})
   if(EXISTS "$ENV{LTTRACE_HOME}/lib/liblttrace.a")
      message(STATUS "LTTRACE Detected")
      include_directories($ENV{LTTRACE_HOME}/include)
      set (LTTRACE_FOUND on)
      set (HAVE_LTTRACE ON)
      set (LTTRACE_LIBRARIES -L$ENV{LTTRACE_HOME}/lib -llttrace -lz -luuid -Wl,-rpath -Wl,$ENV{LTTRACE_HOME}/lib)
   endif (EXISTS "$ENV{LTTRACE_HOME}/lib/liblttrace.a")
endif (DEFINED ENV{LTTRACE_HOME})

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

add_subdirectory(mesh)
add_subdirectory(mesh/tests)
add_subdirectory(l7)
add_subdirectory(l7/tests)
add_subdirectory(kdtree)
add_subdirectory(hsfc)
add_subdirectory(s7)
add_subdirectory(timer)
add_subdirectory(memstats)
if (OPENCL_FOUND)
   add_subdirectory(ezcl)
endif (OPENCL_FOUND)
add_subdirectory(zorder)
add_subdirectory(genmalloc)
add_subdirectory(MallocPlus)
add_subdirectory(hash)
add_subdirectory(crux)

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/state_kern.cl > state_kernel.inc
                  DEPENDS state_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(state_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc)

########### clamr target ###############

if (OPENCL_FOUND AND MPI_FOUND)
   set(clamr_CXX_SRCS clamr.cpp state.cpp input.cpp)
   set(clamr_C_SRCS display.c)
   set(clamr_H_SRCS display.h state.h input.h)
   set(clamr_SRCS ${clamr_CXX_SRCS} ${clamr_C_SRCS} ${clamr_H_SRCS})

   add_executable(clamr ${clamr_SRCS})

   set_target_properties(clamr PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(clamr dpmesh crux hsfc dhash kdtree zorder s7 pezcl timer memstats dl7 genmalloc dpMallocPlus m)
   target_link_libraries(clamr ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(clamr ${OPENCL_LIBRARIES})
   target_link_libraries(clamr ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   target_link_libraries(clamr ${MPI_LIBRARIES})
   add_dependencies(clamr state_kernel_source)

   install(TARGETS clamr DESTINATION bin)
endif (OPENCL_FOUND AND MPI_FOUND)

########### clamr_cpuonly target ###############

set(clamr_cpuonly_CXX_SRCS clamr_cpuonly.cpp state.cpp input.cpp)
set(clamr_cpuonly_C_SRCS display.c)
set(clamr_cpuonly_H_SRCS display.h state.h input.h)

set(clamr_cpuonly_SRCS ${clamr_cpuonly_CXX_SRCS} ${clamr_cpuonly_C_SRCS} ${clamr_cpuonly_H_SRCS})

add_executable(clamr_cpuonly ${clamr_cpuonly_SRCS})

target_link_libraries(clamr_cpuonly mesh crux hsfc hash kdtree zorder s7 timer memstats genmalloc MallocPlus m)
target_link_libraries(clamr_cpuonly ${MPE_NOMPI_LIBS} ${X11_LIBS})
if (REPROBLAS_FOUND)
   target_link_libraries(clamr_cpuonly ${REPROBLAS_LIBRARIES})
endif (REPROBLAS_FOUND)
target_link_libraries(clamr_cpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS clamr_cpuonly DESTINATION bin)

########### clamr_gpuonly target ###############
if (OPENCL_FOUND)
   set(clamr_gpuonly_CXX_SRCS clamr_gpuonly.cpp state.cpp input.cpp)
   set(clamr_gpuonly_C_SRCS display.c)
   set(clamr_gpuonly_H_SRCS display.h state.h input.h)
   set(clamr_gpuonly_SRCS ${clamr_gpuonly_CXX_SRCS} ${clamr_gpuonly_C_SRCS} ${clamr_gpuonly_H_SRCS})

   add_executable(clamr_gpuonly ${clamr_gpuonly_SRCS})

   set_target_properties(clamr_gpuonly PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

   target_link_libraries(clamr_gpuonly dmesh crux hsfc dhash kdtree zorder s7 ezcl timer memstats genmalloc dMallocPlus m)
   target_link_libraries(clamr_gpuonly ${MPE_NOMPI_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_gpuonly ${OPENCL_LIBRARIES})
   target_link_libraries(clamr_gpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_gpuonly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   add_dependencies(clamr_gpuonly state_kernel_source)

   install(TARGETS clamr_gpuonly DESTINATION bin)
endif (OPENCL_FOUND)

########### clamr_mpionly target ###############
if (MPI_FOUND)
   set(clamr_mpionly_CXX_SRCS clamr_mpionly.cpp state.cpp input.cpp)
   set(clamr_mpionly_C_SRCS display.c)
   set(clamr_mpionly_H_SRCS display.h state.h input.h)

   set(clamr_mpionly_SRCS ${clamr_mpionly_CXX_SRCS} ${clamr_mpionly_C_SRCS} ${clamr_mpionly_H_SRCS})

   add_executable(clamr_mpionly ${clamr_mpionly_SRCS})

   set_target_properties(clamr_mpionly PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(clamr_mpionly pmesh crux hsfc hash kdtree zorder s7 timer memstats l7 genmalloc pMallocPlus m)
   target_link_libraries(clamr_mpionly ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_mpionly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
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
if(OPENMP_FOUND)
   set(clamr_openmponly_CXX_SRCS clamr_openmponly.cpp state.cpp input.cpp)
   set(clamr_openmponly_C_SRCS display.c)
   set(clamr_openmponly_H_SRCS display.h state.h input.h)

   set(clamr_openmponly_SRCS ${clamr_openmponly_CXX_SRCS} ${clamr_openmponly_C_SRCS} ${clamr_openmponly_H_SRCS})

   add_executable(clamr_openmponly ${clamr_openmponly_SRCS})

   set_target_properties(clamr_openmponly PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(clamr_openmponly PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(clamr_openmponly PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")

   target_link_libraries(clamr_openmponly mpmesh crux hsfc hash kdtree zorder s7 timer memstats genmalloc MallocPlus m)
   target_link_libraries(clamr_openmponly ${MPE_NOMPI_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_openmponly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_openmponly ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)

   install(TARGETS clamr_openmponly DESTINATION bin)
endif(OPENMP_FOUND)

########### clamr_mpiopenmponly target ###############
if (MPI_FOUND AND OPENMP_FOUND)
   set(clamr_mpiopenmponly_CXX_SRCS clamr_mpiopenmponly.cpp state.cpp input.cpp)
   set(clamr_mpiopenmponly_C_SRCS display.c)
   set(clamr_mpiopenmponly_H_SRCS display.h state.h input.h)

   set(clamr_mpiopenmponly_SRCS ${clamr_mpiopenmponly_CXX_SRCS} ${clamr_mpiopenmponly_C_SRCS} ${clamr_mpiopenmponly_H_SRCS})

   add_executable(clamr_mpiopenmponly ${clamr_mpiopenmponly_SRCS})

   set_target_properties(clamr_mpiopenmponly PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI")

   set(clamr_mpiopenmponly_link_flags "${OpenMP_CXX_FLAGS}")
   set_target_properties(clamr_mpiopenmponly PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(clamr_mpiopenmponly PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(clamr_mpiopenmponly PROPERTIES LINK_FLAGS "${clamr_mpiopenmponly_link_flags}")

   target_link_libraries(clamr_mpiopenmponly pmesh crux hsfc hash kdtree zorder s7 timer memstats l7 genmalloc pMallocPlus m)
   target_link_libraries(clamr_mpiopenmponly ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_mpiopenmponly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
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
    set(clamr_quo_C_SRCS display.c)
    set(clamr_quo_H_SRCS display.h state.h input.h)

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

    target_link_libraries(clamr_quo pmesh crux hsfc hash kdtree zorder s7 timer memstats l7 genmalloc pMallocPlus m)
    target_link_libraries(clamr_quo quo)
    target_link_libraries(clamr_quo j7)
    target_link_libraries(clamr_quo ${MPE_LIBS} ${X11_LIBS})
    target_link_libraries(clamr_quo ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
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
if (OPENCL_FOUND)
   set(clamr_gpucheck_CXX_SRCS clamr_gpucheck.cpp state.cpp input.cpp)
   set(clamr_gpucheck_C_SRCS display.c)
   set(clamr_gpucheck_H_SRCS display.h state.h input.h)

   set(clamr_gpucheck_SRCS ${clamr_gpucheck_CXX_SRCS} ${clamr_gpucheck_C_SRCS} ${clamr_gpucheck_H_SRCS})

   add_executable(clamr_gpucheck ${clamr_gpucheck_SRCS})

   set_target_properties(clamr_gpucheck PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

   target_link_libraries(clamr_gpucheck dmesh crux hsfc dhash kdtree zorder s7 ezcl timer memstats genmalloc dMallocPlus m)
   target_link_libraries(clamr_gpucheck ${MPE_NOMPI_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_gpucheck ${OPENCL_LIBRARIES})
   target_link_libraries(clamr_gpucheck ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_gpucheck ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   add_dependencies(clamr_gpucheck state_kernel_source)

   install(TARGETS clamr_gpucheck DESTINATION bin)
endif (OPENCL_FOUND)

########### clamr_mpicheck target ###############
if (MPI_FOUND)
   set(clamr_mpicheck_CXX_SRCS clamr_mpicheck.cpp state.cpp input.cpp)
   set(clamr_mpicheck_C_SRCS display.c)
   set(clamr_mpicheck_H_SRCS display.h state.h input.h)

   set(clamr_mpicheck_SRCS ${clamr_mpicheck_CXX_SRCS} ${clamr_mpicheck_C_SRCS} ${clamr_mpicheck_H_SRCS})

   add_executable(clamr_mpicheck ${clamr_mpicheck_SRCS})

   set_target_properties(clamr_mpicheck PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(clamr_mpicheck pmesh crux hsfc hash kdtree zorder s7 timer memstats l7 genmalloc pMallocPlus m)
   target_link_libraries(clamr_mpicheck ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_mpicheck ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
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
if (MPI_FOUND AND OPENCL_FOUND)
   set(clamr_checkall_CXX_SRCS clamr_checkall.cpp state.cpp input.cpp)
   set(clamr_checkall_C_SRCS display.c)
   set(clamr_checkall_H_SRCS display.h state.h input.h)

   set(clamr_checkall_SRCS ${clamr_checkall_CXX_SRCS} ${clamr_checkall_C_SRCS} ${clamr_checkall_H_SRCS})

   add_executable(clamr_checkall ${clamr_checkall_SRCS})

   set_target_properties(clamr_checkall PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(clamr_checkall dpmesh crux hsfc dhash kdtree zorder s7 pezcl timer memstats dl7 genmalloc dpMallocPlus m)
   target_link_libraries(clamr_checkall ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(clamr_checkall ${OPENCL_LIBRARIES})
   target_link_libraries(clamr_checkall ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   if (REPROBLAS_FOUND)
      target_link_libraries(clamr_checkall ${REPROBLAS_LIBRARIES})
   endif (REPROBLAS_FOUND)
   if (LTTRACE_FOUND)
      target_link_libraries(clamr_checkall ${LTTRACE_LIBRARIES})
   endif (LTTRACE_FOUND)
   target_link_libraries(clamr_checkall ${MPI_LIBRARIES})
   add_dependencies(clamr_checkall state_kernel_source)

   install(TARGETS clamr_checkall DESTINATION bin)
endif (MPI_FOUND AND OPENCL_FOUND)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "state_kernel.inc;clamr;clamr_cpuonly;clamr_gpuonly;clamr_mpionly;clamr_gpucheck;clamr_mpicheck;clamr_checkall;clamr_quo")

########### install files ###############

set (CMAKE_CHECK_COMMAND make -C ${CMAKE_SOURCE_DIR}/mesh/tests mesh_check &&
                         make -C ${CMAKE_SOURCE_DIR}/l7/tests   l7_check)

add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND})

# configure a header file to pass some of the CMake settings
# to the source code
# at bottom because we set things that change config above
configure_file ( "${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")
