########### global settings ###############
set(H_SRCS MallocPlus.h)
set(CXX_SRCS MallocPlus.cpp)


########### MallocPlus target ###############
set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(MallocPlus STATIC ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS MallocPlus DESTINATION lib)

########### pMallocPlus target ###############
if (MPI_FOUND)
   set(pMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(pMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(pMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(pMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(pMallocPlus ${MPI_LIBRARIES})
   install(TARGETS pMallocPlus DESTINATION lib)
endif (MPI_FOUND)

########### mpMallocPlus target ###############
if (OPENMP_FOUND)
   set(mpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(mpMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(mpMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(mpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(mpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(mpMallocPlus PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

   target_link_libraries(mpMallocPlus ${MPI_LIBRARIES})
   install(TARGETS mpMallocPlus DESTINATION lib)
endif (OPENMP_FOUND)

########### pmpMallocPlus target ###############
if (OPENMP_FOUND AND MPI_FOUND)
   set(pmpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(pmpMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(pmpMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(pmpMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   set_target_properties(pmpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(pmpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(pmpMallocPlus PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

   target_link_libraries(pmpMallocPlus ${MPI_LIBRARIES})
   install(TARGETS pmpMallocPlus DESTINATION lib)
endif (OPENMP_FOUND AND MPI_FOUND)

########### dMallocPlus target ###############
if (OPENCL_FOUND)
   set(dMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(dMallocPlus STATIC ${dMallocPlus_LIB_SRCS})

   set_target_properties(dMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_OPENCL)
   target_link_libraries(dMallocPlus ${OPENCL_LIBRARIES})
   target_link_libraries(dMallocPlus ezcl)
   #if (Boost_FOUND AND MPI_FOUND)
   #   message(STATUS "Adding j7 support ot MallocPlus")
   #   target_link_libraries(dMallocPlus j7)
   #endif()
   install(TARGETS dMallocPlus DESTINATION lib)
endif(OPENCL_FOUND)

########### dpMallocPlus target ###############
if (OPENCL_FOUND)
   set(dpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(dpMallocPlus STATIC ${dpMallocPlus_LIB_SRCS})

   set_target_properties(dpMallocPlus PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dpMallocPlus PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(dpMallocPlus ${MPI_LIBRARIES})
   target_link_libraries(dpMallocPlus ${OPENCL_LIBRARIES})
   target_link_libraries(dpMallocPlus ezcl)
   #if (Boost_FOUND AND MPI_FOUND)
   #   message(STATUS "Adding j7 support ot MallocPlus")
   #   target_link_libraries(dMallocPlus j7)
   #endif()
   install(TARGETS dMallocPlus DESTINATION lib)
endif(OPENCL_FOUND)

########### install files ###############

#install(FILES  MallocPlus.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#default: libMallocPlus.la
#all: libMallocPlus.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libMallocPlus.la
#include_HEADERS = MallocPlus.h
#
#H_SRCS = MallocPlus.h
#
#CXX_SRCS = MallocPlus.cpp
#
#libMallocPlus_la_SOURCES = ${CXX_SRCS} ${H_SRCS}
#
