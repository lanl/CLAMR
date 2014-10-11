cmake_minimum_required(VERSION 2.6)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

########### global settings ###############
set(H_SRCS graphics.h display.h)
set(C_SRCS graphics.c display.c)

set(graphics_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### graphics target ###############

add_library(graphics SHARED ${graphics_LIB_SRCS})

set_target_properties(graphics PROPERTIES VERSION 2.0.0 SOVERSION 2)

target_link_libraries(graphics ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_link_libraries(graphics ${MPE_NOMPI_LIBS} ${X11_LIBS})
if (ImageMagick_MagickWand_FOUND)
   target_link_libraries(graphics ${ImageMagick_LIBRARIES})
endif (ImageMagick_MagickWand_FOUND)

install(TARGETS graphics DESTINATION lib)

########### pgraphics target ###############
if (MPI_FOUND)
   add_library(pgraphics SHARED ${graphics_LIB_SRCS})

   set_target_properties(pgraphics PROPERTIES VERSION 2.0.0 SOVERSION 2)
   set_target_properties(pgraphics PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)

   target_link_libraries(pgraphics ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   target_link_libraries(pgraphics ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(pgraphics ${MPI_LIBRARIES})
   if (ImageMagick_MagickWand_FOUND)
      target_link_libraries(pgraphics ${ImageMagick_LIBRARIES})
   endif (ImageMagick_MagickWand_FOUND)

   install(TARGETS pgraphics DESTINATION lib)
endif (MPI_FOUND)

########### install files ###############

#install(FILES  display.h DESTINATION include)
