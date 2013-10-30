########### global settings ###############
set(H_SRCS hash.h)

set(C_SRCS hash.c)

set(hash_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hash/hashlib_kern.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/hash/hashlib_kern.cl > hashlib_kern.inc
                  DEPENDS hashlib_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(hashlib_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/hash/hashlib_kern.inc)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hash/hashlib_source_kern.inc
                  COMMAND ${CMAKE_SOURCE_DIR}/embed_source.pl ${CMAKE_SOURCE_DIR}/hash/hashlib_source_kern.cl > hashlib_source_kern.inc
                  DEPENDS hashlib_source_kern.cl ${CMAKE_SOURCE_DIR}/embed_source.pl)
add_custom_target(hashlib_source_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/hash/hashlib_source_kern.inc)

########### hash target ###############

add_library(hash STATIC ${hash_LIB_SRCS})

set_target_properties(hash PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS hash DESTINATION lib)

########### dhash target ###############

add_library(dhash STATIC ${hash_LIB_SRCS})

set_target_properties(dhash PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(dhash PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
install(TARGETS dhash DESTINATION lib)
add_dependencies(dhash hashlib_source_kernel_source)
add_dependencies(dhash hashlib_kernel_source)

########### install files ###############

install(FILES  hash.h DESTINATION include)

