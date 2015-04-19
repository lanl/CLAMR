cmake_minimum_required(VERSION 2.6)

project (GENMALLOC)
# The version number.
set (GENMALLOC_VERSION_MAJOR 2)
set (GENMALLOC_VERSION_MINOR 0.7)
set (GENMALLOC_VERSION "${GENMALLOC_VERSION_MAJOR}.${GENMALLOC_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")

########### global settings ###############
set(H_SRCS genmalloc.h)
set(C_SRCS genmalloc.c)

set(genmalloc_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### genmalloc target ###############

add_library(genmalloc SHARED ${genmalloc_LIB_SRCS})

set_target_properties(genmalloc PROPERTIES VERSION ${GENMALLOC_VERSION} SOVERSION ${GENMALLOC_VERSION_MAJOR})
install(TARGETS genmalloc DESTINATION lib)

########### install files ###############

install(FILES  genmalloc.h DESTINATION include/genmalloc)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND pdfdoc
                     COMMAND htmldoc
                     COMMAND singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/genmalloc.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/genmalloc.pdf)
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
   add_custom_target(genmalloc_doc
                     COMMAND make genmalloc_pdfdoc
                     COMMAND make genmalloc_htmldoc
                     COMMAND make genmalloc_singlehtmldoc)
   set_target_properties(genmalloc_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/genmalloc.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(genmalloc_pdfdoc DEPENDS _build/latex/genmalloc.pdf)
   set_target_properties(genmalloc_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(genmalloc_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(genmalloc_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(genmalloc_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(genmalloc_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")
