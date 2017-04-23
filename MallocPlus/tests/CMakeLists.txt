########### global settings ###############

set(H_SRCS)

set(CXX_SRCS MallocPlusTest.cpp)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   include_directories(${CMAKE_SOURCE_DIR})
   include_directories(${CMAKE_SOURCE_DIR}/..)
else ()
   include_directories(${CMAKE_SOURCE_DIR}/MallocPlus)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### MallocPlusTest target ###############
set(MallocPlusTest_SRCS ${CXX_SRCS} ${H_SRCS})

add_executable(MallocPlusTest EXCLUDE_FROM_ALL ${MallocPlusTest_SRCS})

set_target_properties(MallocPlusTest PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(MallocPlusTest PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(MallocPlusTest MallocPlus)

########### install files ###############

################# check/test ##################

add_test(MallocPlusTest MallocPlusTest)

set(CMAKE_CHECK_COMMAND ctest && cat Testing/Temporary/LastTest.log)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS MallocPlusTest)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(MallocPlus_check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS MallocPlusTest)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "MallocPlusTest")

