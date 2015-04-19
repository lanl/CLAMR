cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (CRUX)
# The version number.
set (CRUX_VERSION_MAJOR 2)
set (CRUX_VERSION_MINOR 0.7)
set (CRUX_VERSION "${CRUX_VERSION_MAJOR}.${CRUX_VERSION_MINOR}")

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
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")

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

set_target_properties(crux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR})
install(TARGETS crux DESTINATION lib)

########### dcrux target ###############

#if (OPENCL_FOUND)
#   add_library(dcrux STATIC ${crux_LIB_SRCS})

#   set_target_properties(dcrux PROPERTIES VERSION ${CRUX_VERSION} SOVERSION ${CRUX_VERSION_MAJOR})
#   set_target_properties(dcrux PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
#   install(TARGETS dcrux DESTINATION lib)
#   add_dependencies(dcrux cruxlib_source_kernel_source)
#   add_dependencies(dcrux cruxlib_kernel_source)
#endif(OPENCL_FOUND)

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

