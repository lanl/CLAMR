cmake_minimum_required(VERSION 3.1.3)

project (CRUX)
# The version number.
set (CRUX_VERSION_MAJOR 2)
set (CRUX_VERSION_MINOR 0.7)
set (CRUX_VERSION "${CRUX_VERSION_MAJOR}.${CRUX_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

if (NOT CMAKE_BUILD_TYPE)
   #set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

set (CMAKE_C_STANDARD 99)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -fstrict-aliasing -ftree-vectorize -fopenmp-simd")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fopt-info-vec-optimized")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -mtune=native")
   elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mcpu=powerpc64le")
   else ()
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
   endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -ansi-alias -vecabi=cmdtarget -xHost -qopenmp-simd -qopt-report=5 -qopt-report-phase=vec,loop")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_C_FLAGS_RELEASE "-qopt-zmm-usage=high ${CMAKE_C_FLAGS_RELEASE}")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
endif()

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
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -fstrict-aliasing -ftree-vectorize -fopenmp-simd")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopt-info-vec-optimized")
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

########### global settings ###############
set(H_SRCS crux.h)

set(CXX_SRCS crux.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   set (C_SRCS fmemopen.c)
   set (H_SRCS ${H_SRCS} fmemopen.h)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(crux_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${H_SRCS})

########### crux target ###############

add_library(crux STATIC ${crux_LIB_SRCS})

##set_target_properties(crux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR} COMPILE_DEFINITIONS "HAVE_MPI")
set_target_properties(crux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR})
if (HDF5_FOUND)
   target_link_libraries(crux ${HDF5_LIBS})
endif (HDF5_FOUND)
install(TARGETS crux DESTINATION lib)

########### pcrux target ###############

add_library(pcrux STATIC ${crux_LIB_SRCS})

set_target_properties(pcrux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR} COMPILE_DEFINITIONS "HAVE_MPI")
set_target_properties(pcrux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR})
if (HDF5_FOUND)
   target_link_libraries(pcrux ${HDF5_LIBS})
endif (HDF5_FOUND)
install(TARGETS pcrux DESTINATION lib)

########### dcrux target ###############

#if (OpenCL_FOUND)
#   add_library(dcrux STATIC ${crux_LIB_SRCS})

#   set_target_properties(dcrux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR})
#   set_target_properties(dcrux PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
#   install(TARGETS dcrux DESTINATION lib)
#   add_dependencies(dcrux cruxlib_source_kernel_source)
#   add_dependencies(dcrux cruxlib_kernel_source)
#endif(OpenCL_FOUND)

########### install files ###############

install(FILES  crux.h DESTINATION include/crux)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/Crux.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/Crux.pdf)
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
   add_custom_target(Crux_doc
                     COMMAND make Crux_pdfdoc
                     COMMAND make Crux_htmldoc
                     COMMAND make Crux_singlehtmldoc)
   set_target_properties(Crux_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/Crux.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(Crux_pdfdoc DEPENDS _build/latex/Crux.pdf)
   set_target_properties(Crux_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(Crux_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(Crux_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(Crux_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(Crux_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")

