########### global settings ###############
set(H_SRCS zorder.h)
set(C_SRCS zorder.c)

set(zorder_LIB_SRCS ${C_SRCS} ${H_SRCS})

set_source_files_properties(${C_SRCS} PROPERTIES COMPILE_FLAGS "${VECTOR_C_FLAGS}")
########### zorder target ###############
add_library(zorder SHARED ${zorder_LIB_SRCS})

set_target_properties(zorder PROPERTIES VERSION 2.0.0 SOVERSION 2)
target_link_libraries(zorder s7)
install(TARGETS zorder DESTINATION lib)

########### install files ###############

#install(FILES  zorder.h DESTINATION include)
