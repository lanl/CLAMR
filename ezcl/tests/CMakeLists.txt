########### global settings ###############

set(H_SRCS)

set(C_SRCS ezclTest.c)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   include_directories(${CMAKE_SOURCE_DIR})
   include_directories(${CMAKE_SOURCE_DIR}/..)
else ()
   include_directories(${CMAKE_SOURCE_DIR}/ezcl)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### ezclTest target ###############
set(ezclTest_SRCS ${C_SRCS} ${H_SRCS})

add_executable(ezclTest EXCLUDE_FROM_ALL ${ezclTest_SRCS})

set_target_properties(ezclTest PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(ezclTest PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(ezclTest ezcl)

########### install files ###############

################# check/test ##################

add_test(ezclTest ezclTest)

set(CMAKE_CHECK_COMMAND ctest && cat Testing/Temporary/LastTest.log)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS ezclTest)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(ezcl_check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS ezclTest)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "ezclTest")

