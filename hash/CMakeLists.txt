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

########### mphash target ###############
if (OPENMP_FOUND)
   add_library(mphash STATIC ${hash_LIB_SRCS})

   set_target_properties(mphash PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(mphash PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(mphash PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(mphash PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
   install(TARGETS mphash DESTINATION lib)
endif (OPENMP_FOUND)

########### dhash target ###############

if (OPENCL_FOUND)
   add_library(dhash STATIC ${hash_LIB_SRCS})

   set_target_properties(dhash PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(dhash PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
   install(TARGETS dhash DESTINATION lib)
   add_dependencies(dhash hashlib_source_kernel_source)
   add_dependencies(dhash hashlib_kernel_source)
endif(OPENCL_FOUND)

########### install files ###############

install(FILES  hash.h DESTINATION include)

