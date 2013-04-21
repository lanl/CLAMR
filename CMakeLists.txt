cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

project (CLAMR)
# The version number.
set (CLAMR_VERSION_MAJOR 2)
set (CLAMR_VERSION_MINOR 0.7)
set (CLAMR_VERSION "${CLAMR_VERSION_MAJOR}.${CLAMR_VERSION_MINOR}")
 
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_C_FLAGS "-std=c99")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")

find_package(OpenCL REQUIRED)
#  OPENCL_FOUND        : TRUE if found, FALSE otherwise
#  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
#  OPENCL_LIBRARIES    : The libraries to link against
if (OPENCL_FOUND)
   include_directories(${OPENCL_INCLUDE_DIRS})
endif (OPENCL_FOUND)
find_package(MPI REQUIRED)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if (NOT GRAPHICS_TYPE)
   set (GRAPHICS_TYPE OpenGL)
endif (NOT GRAPHICS_TYPE)

if (GRAPHICS_TYPE STREQUAL "OpenGL")
   if (APPLE)
      set(OPENGL_FOUND yes)
      set(OPENGL_LIBRARIES -L/usr/X11/lib GLU glut GL)
   else ()
      find_package(OpenGL)
   endif ()
endif (GRAPHICS_TYPE STREQUAL "OpenGL")

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
   if (NOT APPLE)
      find_package(GLUT)
   endif()
endif (OPENGL_FOUND)

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

#message("MPE_LIBRARIES ${MPE_LIBRARIES}")
#message("MPE_LIBS ${MPE_LIBS}")
#message("MPE_NOMPI_LIBS ${MPE_NOMPI_LIBS}")
#message("OPENGL_LIBRARIES ${OPENGL_LIBRARIES}")
#message("GLUT_LIBRARIES ${GLUT_LIBRARIES}")
#message("X11_LIBRARIES ${X11_LIBRARIES}")
#message("X11_LIBS ${X11_LIBS}")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file ( "${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")
 
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
add_subdirectory(ezcl)
add_subdirectory(zorder)
add_subdirectory(genmalloc)
add_subdirectory(MallocPlus)

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/state_kern.cl > state_kernel.inc
                  DEPENDS state_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(state_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/state_kernel.inc)

########### clamr target ###############

set(clamr_CXX_SRCS clamr.cpp state.cpp input.cpp)
set(clamr_C_SRCS display.c)
set(clamr_H_SRCS display.h state.h input.h)
set(clamr_SRCS ${clamr_CXX_SRCS} ${clamr_C_SRCS} ${clamr_H_SRCS})

add_executable(clamr ${clamr_SRCS})

set_target_properties(clamr PROPERTIES COMPILE_DEFINITIONS HAVE_MPI_AND_OPENCL)

target_link_libraries(clamr dpmesh hsfc kdtree zorder s7 pezcl timer memstats dl7 genmalloc dMallocPlus m)
target_link_libraries(clamr ${MPE_LIBS} ${X11_LIBS})
target_link_libraries(clamr ${MPI_LIBRARIES})
target_link_libraries(clamr ${OPENCL_LIBRARIES})
target_link_libraries(clamr ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
add_dependencies(clamr state_kernel_source)

install(TARGETS clamr DESTINATION bin)

########### clamr_gpuonly target ###############
set(clamr_gpuonly_CXX_SRCS clamr_gpuonly.cpp state.cpp input.cpp)
set(clamr_gpuonly_C_SRCS display.c)
set(clamr_gpuonly_H_SRCS display.h state.h input.h)
set(clamr_gpuonly_SRCS ${clamr_gpuonly_CXX_SRCS} ${clamr_gpuonly_C_SRCS} ${clamr_gpuonly_H_SRCS})

add_executable(clamr_gpuonly ${clamr_gpuonly_SRCS})

set_target_properties(clamr_gpuonly PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

target_link_libraries(clamr_gpuonly dmesh hsfc kdtree zorder s7 ezcl timer memstats genmalloc dMallocPlus m)
target_link_libraries(clamr_gpuonly ${MPE_NOMPI_LIBS} ${X11_LIBS})
target_link_libraries(clamr_gpuonly ${OPENCL_LIBRARIES})
target_link_libraries(clamr_gpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
add_dependencies(clamr_gpuonly state_kernel_source)

install(TARGETS clamr_gpuonly DESTINATION bin)

########### clamr_cpuonly target ###############

set(clamr_cpuonly_CXX_SRCS clamr_cpuonly.cpp state.cpp input.cpp)
set(clamr_cpuonly_C_SRCS display.c)
set(clamr_cpuonly_H_SRCS display.h state.h input.h)

set(clamr_cpuonly_SRCS ${clamr_cpuonly_CXX_SRCS} ${clamr_cpuonly_C_SRCS} ${clamr_cpuonly_H_SRCS})

add_executable(clamr_cpuonly ${clamr_cpuonly_SRCS})

target_link_libraries(clamr_cpuonly mesh hsfc kdtree zorder s7 timer memstats genmalloc MallocPlus m)
target_link_libraries(clamr_cpuonly ${MPE_NOMPI_LIBS} ${X11_LIBS})
target_link_libraries(clamr_cpuonly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS clamr_cpuonly DESTINATION bin)

########### clamr_mpionly target ###############
set(clamr_mpionly_CXX_SRCS clamr_mpionly.cpp state.cpp input.cpp)
set(clamr_mpionly_C_SRCS display.c)
set(clamr_mpionly_H_SRCS display.h state.h input.h)

set(clamr_mpionly_SRCS ${clamr_mpionly_CXX_SRCS} ${clamr_mpionly_C_SRCS} ${clamr_mpionly_H_SRCS})

add_executable(clamr_mpionly ${clamr_mpionly_SRCS})

set_target_properties(clamr_mpionly PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)

target_link_libraries(clamr_mpionly pmesh hsfc kdtree zorder s7 timer memstats l7 genmalloc MallocPlus m)
target_link_libraries(clamr_mpionly ${MPE_LIBS} ${X11_LIBS})
target_link_libraries(clamr_mpionly ${MPI_LIBRARIES})
target_link_libraries(clamr_mpionly ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS clamr_mpionly DESTINATION bin)

########### clamr_gpucheck target ###############
set(clamr_gpucheck_CXX_SRCS clamr_gpucheck.cpp state.cpp input.cpp)
set(clamr_gpucheck_C_SRCS display.c)
set(clamr_gpucheck_H_SRCS display.h state.h input.h)

set(clamr_gpucheck_SRCS ${clamr_gpucheck_CXX_SRCS} ${clamr_gpucheck_C_SRCS} ${clamr_gpucheck_H_SRCS})

add_executable(clamr_gpucheck ${clamr_gpucheck_SRCS})

set_target_properties(clamr_gpucheck PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)

target_link_libraries(clamr_gpucheck dmesh hsfc kdtree zorder s7 ezcl timer memstats genmalloc dMallocPlus m)
target_link_libraries(clamr_gpucheck ${MPE_NOMPI_LIBS} ${X11_LIBS})
target_link_libraries(clamr_gpucheck ${OPENCL_LIBRARIES})
target_link_libraries(clamr_gpucheck ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
add_dependencies(clamr_gpucheck state_kernel_source)

install(TARGETS clamr_gpucheck DESTINATION bin)

########### clamr_mpicheck target ###############
set(clamr_mpicheck_CXX_SRCS clamr_mpicheck.cpp state.cpp input.cpp)
set(clamr_mpicheck_C_SRCS display.c)
set(clamr_mpicheck_H_SRCS display.h state.h input.h)

set(clamr_mpicheck_SRCS ${clamr_mpicheck_CXX_SRCS} ${clamr_mpicheck_C_SRCS} ${clamr_mpicheck_H_SRCS})

add_executable(clamr_mpicheck ${clamr_mpicheck_SRCS})

set_target_properties(clamr_mpicheck PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)

target_link_libraries(clamr_mpicheck pmesh hsfc kdtree zorder s7 timer memstats l7 genmalloc MallocPlus m)
target_link_libraries(clamr_mpicheck ${MPE_LIBS} ${X11_LIBS})
target_link_libraries(clamr_mpicheck ${MPI_LIBRARIES})
target_link_libraries(clamr_mpicheck ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

install(TARGETS clamr_mpicheck DESTINATION bin)

########### clamr_checkall target ###############
set(clamr_checkall_CXX_SRCS clamr_checkall.cpp state.cpp input.cpp)
set(clamr_checkall_C_SRCS display.c)
set(clamr_checkall_H_SRCS display.h state.h input.h)

set(clamr_checkall_SRCS ${clamr_checkall_CXX_SRCS} ${clamr_checkall_C_SRCS} ${clamr_checkall_H_SRCS})

add_executable(clamr_checkall ${clamr_checkall_SRCS})

set_target_properties(clamr_checkall PROPERTIES COMPILE_DEFINITIONS HAVE_MPI_AND_OPENCL)

target_link_libraries(clamr_checkall dpmesh hsfc kdtree zorder s7 pezcl timer memstats dl7 genmalloc dMallocPlus m)
target_link_libraries(clamr_checkall ${MPE_LIBS} ${X11_LIBS})
target_link_libraries(clamr_checkall ${MPI_LIBRARIES})
target_link_libraries(clamr_checkall ${OPENCL_LIBRARIES})
target_link_libraries(clamr_checkall ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
add_dependencies(clamr_checkall state_kernel_source)

install(TARGETS clamr_checkall DESTINATION bin)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "state_kernel.inc;clamr;clamr_cpuonly;clamr_gpuonly;clamr_mpionly;clamr_gpucheck;clamr_mpicheck;clamr_checkall")

########### install files ###############

#set (CMAKE_CHECK_COMMAND make -C ${CMAKE_SOURCE_DIR}/mesh/tests mesh_check &&
#                         make -C ${CMAKE_SOURCE_DIR}/l7/tests   l7_check)
set (CMAKE_CHECK_COMMAND make -C ${CMAKE_SOURCE_DIR}/mesh/tests mesh_check)

add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND})
