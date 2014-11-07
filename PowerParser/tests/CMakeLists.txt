########### global settings ###############

set(H_SRCS)

set(C_SRCS
      ParserTest.cc
)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   include_directories(${CMAKE_SOURCE_DIR})
   include_directories(${CMAKE_SOURCE_DIR}/..)
else ()
   include_directories(${CMAKE_SOURCE_DIR}/PowerParser)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### ParserTest target ###############
set(ParserTest_SRCS ${C_SRCS} ${H_SRCS})

add_executable(ParserTest EXCLUDE_FROM_ALL ${ParserTest_SRCS})

set_target_properties(ParserTest PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
set_target_properties(ParserTest PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(ParserTest PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(ParserTest PowerParser genmalloc ${MPI_LIBRARIES})
#find_library(GENMALLOC_LIBRARY genmalloc DIRS
#   ${CMAKE_SOURCE_DIR}/genmalloc
#   ${CMAKE_BINARY_DIR}/genmalloc
#   ${CMAKE_BINARY_DIR}
#   ../../
#   ../../tests
#   ../../genmalloc)
#message("GENMALLOC directory is ${GENMALLOC_LIBRARY}")
#target_link_libraries(ParserTest PowerParser ${GENMALLOC_LIBRARY})
#target_link_libraries(ParserTest PowerParser ${MPI_LIBRARIES})

########### install files ###############

################# check/test ##################

add_test(NAME ParserTest COMMAND ParserTest)

set(CMAKE_CHECK_COMMAND ctest)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS ParserTest)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(PowerParser_check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS ParserTest)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "ParserTest")

