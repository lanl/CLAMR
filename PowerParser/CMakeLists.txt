cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (PARSER)
# The version number.
set (PARSER_VERSION_MAJOR 2)
set (PARSER_VERSION_MINOR 0.7)
set (PARSER_VERSION "${PARSER_VERSION_MAJOR}.${PARSER_VERSION_MINOR}")

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

add_subdirectory(tests)

########### global settings ###############
set(H_SRCS   Cmd.hh   Function.hh  Parser_utils.hh  Variable.hh      Word.hh
             Comm.hh  Parse.hh     Parser_math.hh   Restartblock.hh  Whenthen.hh)

set(CXX_SRCS Cmd.cc   Function.cc  Parser_utils.cc  Variable.cc      Word.cc
             Comm.cc  Parse.cc     Parser_math.cc   Restartblock.cc  Whenthen.cc)

########### Parser target ###############
set(Parser_LIB_SRCS ${CXX_SRCS} ${H_SRCS})

add_library(Parser STATIC ${Parser_LIB_SRCS})

set_target_properties(Parser PROPERTIES VERSION ${PARSER_VERSION} SOVERSION ${PARSER_VERSION_MAJOR})
target_link_libraries(Parser ${MPI_LIBRARIES})

install(TARGETS Parser DESTINATION lib)

########### clean files ################

########### install files ###############

install(FILES Parse.hh DESTINATION include)
install(FILES Comm.hh DESTINATION include)
install(FILES Word.hh DESTINATION include)
install(FILES Cmd.hh DESTINATION include)
install(FILES Restartblock.hh DESTINATION include)
install(FILES Whenthen.hh DESTINATION include)
install(FILES Variable.hh DESTINATION include)
install(FILES Parser_utils.hh DESTINATION include)
install(FILES Parser_math.hh DESTINATION include)
install(FILES Function.hh DESTINATION include)
