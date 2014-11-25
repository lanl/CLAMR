########### global settings ###############
set(CMAKE_C_FLAGS "-std=c99 -DHAVE_MPI=1")

set(H_SRCS)

set(C_SRCS
      L7Test.c          update_test.c   reduction_test.c
      broadcast_test.c
)

########### L7Test target ###############
set(L7Test_SRCS ${C_SRCS} ${H_SRCS})

add_executable(L7Test EXCLUDE_FROM_ALL ${L7Test_SRCS})

set_target_properties(L7Test PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
set_target_properties(L7Test PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(L7Test PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
include_directories(${CMAKE_SOURCE_DIR}/l7)
target_link_libraries(L7Test l7 ${MPI_LIBRARIES})

########### install files ###############

################# check/test ##################

add_test(L7Test mpirun -n 2 L7Test)

set (CMAKE_CHECK_COMMAND ctest && cat Testing/Temporary/LastTest.log)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS L7Test)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(l7_check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS L7Test)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "L7Test")

