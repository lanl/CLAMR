########### global settings ###############
set(H_SRCS)

set(CXX_SRCS
      mesh_tests.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/mesh)
include_directories(${CMAKE_SOURCE_DIR}/l7)
include_directories(${CMAKE_SOURCE_DIR}/genmalloc)

########### L7Test target ###############
set(mesh_test_SRCS ${CXX_SRCS} ${H_SRCS})

add_executable(mesh_mpi_test EXCLUDE_FROM_ALL ${mesh_test_SRCS})

set_target_properties(mesh_mpi_test PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)
set_target_properties(mesh_mpi_test PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(mesh_mpi_test PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
target_link_libraries(mesh_mpi_test pmesh crux hash hsfc kdtree zorder s7 timer memstats l7 genmalloc pMallocPlus m)
target_link_libraries(mesh_mpi_test ${MPE_LIBS} ${X11_LIBS})
target_link_libraries(mesh_mpi_test ${MPI_LIBRARIES})
target_link_libraries(mesh_mpi_test ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

set(face_test_SRCS face_tests.cpp)

add_executable(face_test EXCLUDE_FROM_ALL ${face_test_SRCS})

target_link_libraries(face_test mesh crux hash hsfc kdtree zorder s7 timer memstats l7 genmalloc MallocPlus m)
target_link_libraries(mesh_mpi_test ${MPE_LIBS} ${X11_LIBS})
target_link_libraries(mesh_mpi_test ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

########### install files ###############

################# check/test ##################

#add_test(NAME mesh_mpi_test COMMAND mpirun -n 3 ./mesh_mpi_test)
add_test(mesh_mpi_test mpirun -n 3 ./mesh_mpi_test)
add_test(face_test ./face_test)

set (CMAKE_CHECK_COMMAND ctest && cat Testing/Temporary/LastTest.log)

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS mesh_mpi_test face_test)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(mesh_check COMMAND ${CMAKE_CHECK_COMMAND}
                     DEPENDS mesh_mpi_test face_test)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "mesh_mpi_test")

