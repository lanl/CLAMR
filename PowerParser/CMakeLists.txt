cmake_minimum_required(VERSION 2.6)

project (POWERPARSER)
# The version number.
set (POWERPARSER_VERSION_MAJOR 2)
set (POWERPARSER_VERSION_MINOR 0.7)
set (POWERPARSER_VERSION "${POWERPARSER_VERSION_MAJOR}.${POWERPARSER_VERSION_MINOR}")

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

find_package(MPI REQUIRED)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

include_directories(${CMAKE_SOURCE_DIRECTORY})

enable_testing()

add_subdirectory(tests)

########### global settings ###############
set(H_SRCS   Cmd.hh   Function.hh    Parser_utils.hh  Variable.hh      Word.hh
             Comm.hh  PowerParser.hh Parser_math.hh   Restartblock.hh  Whenthen.hh)

set(CXX_SRCS Cmd.cc   Function.cc    Parser_utils.cc  Variable.cc      Word.cc
             Comm.cc  PowerParser.cc Parser_math.cc   Restartblock.cc  Whenthen.cc)

set(DOC_SRCS docs/index.rst docs/About.rst docs/GettingStarted.rst
             docs/DevelopersGuide.rst docs/UsersGuide.rst)

########### PowerParser target ###############
set(PowerParser_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(PowerParser STATIC ${PowerParser_LIB_SRCS})

set_target_properties(PowerParser PROPERTIES VERSION ${POWERPARSER_VERSION} SOVERSION ${POWERPARSER_VERSION_MAJOR})
target_link_libraries(PowerParser ${MPI_LIBRARIES})

install(TARGETS PowerParser DESTINATION lib)

########### doc targets ################

#add_custom_target(doc ALL DEPENDS pdfdoc htmldoc singlehtmldoc)

add_custom_command(OUTPUT _build/latex/PowerParser.pdf
                   COMMAND make -C docs latexpdf >& pdfdoc.out
                   DEPENDS ${DOC_SRCS})
add_custom_target(pdfdoc ALL DEPENDS _build/latex/PowerParser.pdf)

add_custom_command(OUTPUT _build/html/UsersGuide.html
                   COMMAND make -C docs html >& htmldoc.out
                   DEPENDS ${DOC_SRCS})
add_custom_target(htmldoc ALL DEPENDS _build/html/UsersGuide.html)

add_custom_command(OUTPUT _build/singlehtml/index.html
                   COMMAND make -C docs singlehtml >& singlehtmldoc.out
                   DEPENDS ${DOC_SRCS})
add_custom_target(singlehtmldoc ALL DEPENDS _build/singlehtml/index.html)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "install_manifest.txt;docs/_build;htmldoc.out;pdfdoc.out;singlehtmldoc.out")

########### install files ###############

install(FILES PowerParser.hh DESTINATION include)
