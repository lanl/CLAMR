cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (MALLOC_PLUS)
# The version number.
set (MALLOC_PLUS_VERSION_MAJOR 2)
set (MALLOC_PLUS_VERSION_MINOR 0.7)
set (MALLOC_PLUS_VERSION "${MALLOC_PLUS_VERSION_MAJOR}.${MALLOC_PLUS_VERSION_MINOR}")

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

if (MIC_NATIVE)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")

if (NOT DEFINED OPENCL_FOUND AND NOT MIC_NATIVE)
   find_package(OpenCL)
   #  OPENCL_FOUND        : TRUE if found, FALSE otherwise
   #  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
   #  OPENCL_LIBRARIES    : The libraries to link against
   if (OPENCL_FOUND)
      include_directories(${OPENCL_INCLUDE_DIRS})
   endif (OPENCL_FOUND)
endif (NOT DEFINED OPENCL_FOUND AND NOT MIC_NATIVE)

if (NOT DEFINED MPI_FOUND)
   find_package(MPI)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if (NOT DEFINED OPENMP_FOUND)
   find_package( OpenMP)
endif (NOT DEFINED OPENMP_FOUND)

enable_testing()

add_subdirectory(tests)

########### global settings ###############
set(H_SRCS MallocPlus.h)
set(CXX_SRCS MallocPlus.cpp)


########### MallocPlus target ###############
set(MallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(MallocPlus STATIC ${MallocPlus_LIB_SRCS})

set_target_properties(MallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
install(TARGETS MallocPlus DESTINATION lib)

########### pMallocPlus target ###############
if (MPI_FOUND)
   set(pMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(pMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(pMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(pMallocPlus PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
   target_link_libraries(pMallocPlus ${MPI_LIBRARIES})
   install(TARGETS pMallocPlus DESTINATION lib)
endif (MPI_FOUND)

########### mpMallocPlus target ###############
if (OPENMP_FOUND)
   set(mpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(mpMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(mpMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(mpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(mpMallocPlus PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(mpMallocPlus PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

   install(TARGETS mpMallocPlus DESTINATION lib)
endif (OPENMP_FOUND)

########### pmpMallocPlus target ###############
if (OPENMP_FOUND AND MPI_FOUND)
   set(pmpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(pmpMallocPlus STATIC ${MallocPlus_LIB_SRCS})

   set_target_properties(pmpMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
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

   set_target_properties(dMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
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
if (OPENCL_FOUND AND MPI_FOUND)
   set(dpMallocPlus_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

   add_library(dpMallocPlus STATIC ${dpMallocPlus_LIB_SRCS})

   set_target_properties(dpMallocPlus PROPERTIES VERSION ${MALLOC_PLUS_VERSION} SOVERSION ${MALLOC_PLUS_VERSION_MAJOR})
   set_target_properties(dpMallocPlus PROPERTIES COMPILE_DEFINITIONS "HAVE_MPI;HAVE_OPENCL")
   target_link_libraries(dpMallocPlus ${MPI_LIBRARIES})
   target_link_libraries(dpMallocPlus ${OPENCL_LIBRARIES})
   target_link_libraries(dpMallocPlus ezcl)
   #if (Boost_FOUND AND MPI_FOUND)
   #   message(STATUS "Adding j7 support ot MallocPlus")
   #   target_link_libraries(dMallocPlus j7)
   #endif()
   install(TARGETS dpMallocPlus DESTINATION lib)
endif(OPENCL_FOUND AND MPI_FOUND)

########### install files ###############

install(FILES  MallocPlus.h DESTINATION include/MallocPlus)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/MallocPlus.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/MallocPlus.pdf)
   set_target_properties(pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(MallocPlus_doc
                     COMMAND make MallocPlus_pdfdoc
                     COMMAND make MallocPlus_htmldoc
                     COMMAND make MallocPlus_singlehtmldoc)
   set_target_properties(MallocPlus_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/MallocPlus.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(MallocPlus_pdfdoc DEPENDS _build/latex/MallocPlus.pdf)
   set_target_properties(MallocPlus_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(MallocPlus_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(MallocPlus_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(MallocPlus_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(MallocPlus_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")
