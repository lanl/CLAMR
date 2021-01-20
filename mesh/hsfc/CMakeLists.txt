########### global settings ###############
set(H_SRCS hsfc.h)

set(C_SRCS hsfc2part.c hsfc.c hsfcsort.c)

set(hsfc_LIB_SRCS ${C_SRCS} ${H_SRCS})

set_source_files_properties(${C_SRCS} PROPERTIES COMPILE_FLAGS "${VECTOR_C_FLAGS}")

########### hsfc target ###############

add_library(hsfc SHARED ${hsfc_LIB_SRCS})

set_target_properties(hsfc PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS hsfc DESTINATION lib)

########### install files ###############

#install(FILES  hsfc.h DESTINATION include)
