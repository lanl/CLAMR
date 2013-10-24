message(STATUS "j7: Begin Setup")

if (MPI_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    message(STATUS "j7: MPI Found. Using wrapper compiler...")
endif()

set(j7-H_SRCS j7.h)
set(j7-CXX_SRCS j7.cpp)
set(j7-LIB_SRCS ${j7-CXX_SRCS} ${j7-H_SRCS})

add_library(j7 STATIC ${j7-LIB_SRCS})
set_target_properties(j7 PROPERTIES VERSION 0.0.1 SOVERSION 0)
install(TARGETS j7 DESTINATION lib)

# we use Boost magic in j7, so make sure we have Boost
# this min version isn't correct
set (BOOST_MIN_VER "1.45.0")
find_package(Boost ${BOOST_MIN_VER})
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    # FIXME - we should detect if we need librt
    target_link_libraries(j7 rt)
    set_target_properties(j7 PROPERTIES INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})
endif()

## j7-test0
add_executable(j7-test0 j7-test0.cpp)
target_link_libraries(j7-test0 j7)

message(STATUS "j7: End Setup")
