########### global settings ###############
set(H_SRCS crux.h)

set(CXX_SRCS crux.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   set (C_SRCS fmemopen.c)
   set (H_SRCS ${H_SRCS} fmemopen.h)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(crux_LIB_SRCS ${CXX_SRCS} ${C_SRCS} ${H_SRCS})

########### crux target ###############

add_library(crux STATIC ${crux_LIB_SRCS})

set_target_properties(crux PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS crux DESTINATION lib)

########### dcrux target ###############

#if (OPENCL_FOUND)
#   add_library(dcrux STATIC ${crux_LIB_SRCS})

#   set_target_properties(dcrux PROPERTIES VERSION 2.0.0 SOVERSION 2)
#   set_target_properties(dcrux PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
#   install(TARGETS dcrux DESTINATION lib)
#   add_dependencies(dcrux cruxlib_source_kernel_source)
#   add_dependencies(dcrux cruxlib_kernel_source)
#endif(OPENCL_FOUND)

########### install files ###############

install(FILES  crux.h DESTINATION include)

