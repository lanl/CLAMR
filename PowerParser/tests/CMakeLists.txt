########### global settings ###############

set(H_SRCS)

set(CXX_SRCS ParserTest.cc)

set(F90_SRCS FParserTest.f90)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   include_directories(${CMAKE_SOURCE_DIR})
   include_directories(${CMAKE_SOURCE_DIR}/..)
else ()
   include_directories(${CMAKE_SOURCE_DIR}/PowerParser)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### ParserTest target ###############
set(ParserTest_SRCS ${CXX_SRCS} ${H_SRCS})

add_executable(ParserTest EXCLUDE_FROM_ALL ${ParserTest_SRCS})

find_library(GENMALLOC_LIBRARY genmalloc DIRS
   ${CMAKE_BINARY_DIR}/genmalloc
   ${CMAKE_BINARY_DIR}
   ${CMAKE_BINARY_DIR}/../genmalloc
   ../../genmalloc)
#message("GENMALLOC directory is ${GENMALLOC_LIBRARY}")
if (GENMALLOC-FOUND)
   set_target_properties(ParserTest PROPERTIES COMPILE_DEFINITIONS HAVE_GENMALLOC)
endif (GENMALLOC-FOUND)

set_target_properties(ParserTest PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
set_target_properties(ParserTest PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(ParserTest PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(ParserTest PowerParser ${MPI_LIBRARIES})
if (GENMALLOC-FOUND)
   target_link_libraries(ParserTest PowerParser ${GENMALLOC_LIBRARY})
endif (GENMALLOC-FOUND)

########### FParserTest target ###############
set(FParserTest_SRCS ${F90_SRCS})

add_executable(FParserTest EXCLUDE_FROM_ALL ${FParserTest_SRCS})

find_library(GENMALLOC_LIBRARY genmalloc DIRS
   ${CMAKE_BINARY_DIR}/genmalloc
   ${CMAKE_BINARY_DIR}
   ${CMAKE_BINARY_DIR}/../genmalloc
   ../../genmalloc)
#message("GENMALLOC directory is ${GENMALLOC_LIBRARY}")
if (GENMALLOC-FOUND)
   set_target_properties(FParserTest PROPERTIES COMPILE_DEFINITIONS HAVE_GENMALLOC)
endif (GENMALLOC-FOUND)

set_target_properties(FParserTest PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
set_target_properties(FParserTest PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(FParserTest PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(FParserTest PowerParser ${MPI_LIBRARIES})
if (GENMALLOC-FOUND)
   target_link_libraries(FParserTest PowerParser ${GENMALLOC_LIBRARY})
endif (GENMALLOC-FOUND)

########### install files ###############

################# check/test ##################

add_test(ParserTest ParserTest)
add_test(FParserTest FParserTest)

set(CMAKE_CHECK_COMMAND ctest && cat Testing/Temporary/LastTest.log)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS ParserTest)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(PowerParser_check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS ParserTest)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "ParserTest")

