########### global settings ###############
set(H_SRCS memstats.h)
set(C_SRCS memstats.c)

set(memstats_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### memstats target ###############
add_library(memstats SHARED ${memstats_LIB_SRCS})

set_target_properties(memstats PROPERTIES VERSION 2.0.0 SOVERSION 2)
install(TARGETS memstats DESTINATION lib)

########### install files ###############

#install(FILES  memstats.h DESTINATION include)

