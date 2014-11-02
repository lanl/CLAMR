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

add_executable(ParserTest ${ParserTest_SRCS})

set_target_properties(ParserTest PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
set_target_properties(ParserTest PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(ParserTest PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(ParserTest PowerParser genmalloc ${MPI_LIBRARIES})

########### install files ###############

################# test ##################

add_test(ParserTest ParserTest)

set (CMAKE_TEST_COMMAND ParserTest)

add_custom_target(Parser_check COMMAND ${CMAKE_TEST_COMMAND}
                  DEPENDS ParserTest)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "ParserTest")

