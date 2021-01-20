########### global settings ###############
set(H_SRCS Bounds.h Globals.h KDTree.h)

set(C_SRCS Bounds.c KDTree.c)

set_source_files_properties(${C_SRCS} PROPERTIES COMPILE_FLAGS "${VECTOR_C_FLAGS}")

set(kdtree_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### kdtree target ###############

add_library(kdtree SHARED ${kdtree_LIB_SRCS})

set_target_properties(kdtree PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS kdtree DESTINATION lib)

########### install files ###############

#install(FILES  KDTree.h DESTINATION include)

