########### global settings ###############
set(H_SRCS hash.h)

set(C_SRCS hash.c)

set(hash_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### hash target ###############

add_library(hash STATIC ${hash_LIB_SRCS})

set_target_properties(hash PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS hash DESTINATION lib)

########### dhash target ###############

add_library(dhash STATIC ${hash_LIB_SRCS})

set_target_properties(dhash PROPERTIES VERSION 2.0.0 SOVERSION 2)
set_target_properties(dhash PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
install(TARGETS dhash DESTINATION lib)

########### install files ###############

install(FILES  hash.h DESTINATION include)

