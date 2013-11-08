message(STATUS "j7: Begin Setup")

find_package(Librt)

if (MPI_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    message(STATUS "j7: MPI Found. Using wrapper compiler...")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

set(j7-H_SRCS j7.h)
set(j7-CXX_SRCS j7.cpp)
set(j7-LIB_SRCS ${j7-CXX_SRCS} ${j7-H_SRCS})

add_library(j7 STATIC ${j7-LIB_SRCS})
set_target_properties(j7 PROPERTIES VERSION 0.0.1 SOVERSION 0)
install(TARGETS j7 DESTINATION lib)

# we use Boost magic in j7, so make sure we have Boost (top-level does it)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    # FIXME - we should detect if we need librt
    if (Librt_FOUND)
       target_link_libraries(j7 rt)
    endif (Librt_FOUND)
    set_target_properties(j7 PROPERTIES INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})
endif()

## j7-test0
add_executable(j7-test0 j7-test0.cpp)
target_link_libraries(j7-test0 j7)
target_link_libraries(j7-test0 ${MPI_LIBRARIES})

message(STATUS "j7: End Setup")
