########### global settings ###############
set(CMAKE_C_FLAGS "-std=c99 -DHAVE_MPI=1")

set(H_SRCS)

set(C_SRCS
      L7Test.c          update_test.c   reduction_test.c
      broadcast_test.c
)

########### L7Test target ###############
set(L7Test_SRCS ${C_SRCS} ${H_SRCS})

add_executable(L7Test ${L7Test_SRCS})

set_target_properties(L7Test PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
include_directories(${CMAKE_SOURCE_DIR}/l7)
target_link_libraries(L7Test l7 ${MPI_LIBRARIES})

########### install files ###############

################# test ##################

add_test(L7Test L7Test)

set (CMAKE_TEST_COMMAND mpirun -n 2 L7Test)

add_custom_target(check COMMAND ${CMAKE_TEST_COMMAND}
                  DEPENDS L7Test)

