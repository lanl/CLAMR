########### global settings ###############
set(H_SRCS display.h)
set(C_SRCS display.c)

set(graphics_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### graphics target ###############

add_library(graphics SHARED ${graphics_LIB_SRCS})

set_target_properties(graphics PROPERTIES VERSION 2.0.0 SOVERSION 2)

target_link_libraries(graphics ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_link_libraries(graphics ${MPE_LIBS} ${X11_LIBS})
if (ImageMagick_MagickWand_FOUND)
   target_link_libraries(graphics ${ImageMagick_LIBRARIES})
endif (ImageMagick_MagickWand_FOUND)

install(TARGETS graphics DESTINATION lib)

########### install files ###############

#install(FILES  display.h DESTINATION include)
