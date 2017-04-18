cmake_minimum_required(VERSION 2.6)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

#find_package(ImageMagick COMPONENTS MagickCore MagickWand)
find_package(ImageMagick COMPONENTS MagickWand)
#   ImageMagick_FOUND                    - TRUE if all components are found.
#   ImageMagick_INCLUDE_DIRS             - Full paths to all include dirs.
#   ImageMagick_LIBRARIES                - Full paths to all libraries.
#   ImageMagick_<component>_FOUND        - TRUE if <component> is found.
#   ImageMagick_<component>_INCLUDE_DIRS - Full path to <component> include dirs.
#   ImageMagick_<component>_LIBRARIES    - Full path to <component> libraries.
if (ImageMagick_FOUND)
   if (ImageMagick_MagickWand_FOUND)
      add_definitions(-DHAVE_MAGICKWAND)
      add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)
      add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
#     add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)
#     add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
      include_directories(${ImageMagick_INCLUDE_DIRS})
   endif (ImageMagick_MagickWand_FOUND)

#   message("ImageMagick_FOUND                   ${ImageMagick_FOUND}")
#   message("ImageMagick_INCLUDE_DIRS            ${ImageMagick_INCLUDE_DIRS}")
#   message("ImageMagick_LIBRARIES               ${ImageMagick_LIBRARIES}")
##  message("ImageMagick_MagickCore_FOUND        ${ImageMagick_MagickCore_FOUND}")
##  message("ImageMagick_MagickCore_INCLUDE_DIRS ${ImageMagick_MagickCore_INCLUDE_DIRS}")
##  message("ImageMagick_MagickCore_LIBRARIES    ${ImageMagick_MagickCore_LIBRARIES}")
#   message("ImageMagick_MagickWand_FOUND        ${ImageMagick_MagickWand_FOUND}")
#   message("ImageMagick_MagickWand_INCLUDE_DIRS ${ImageMagick_MagickWand_INCLUDE_DIRS}")
#   message("ImageMagick_MagickWand_LIBRARIES    ${ImageMagick_MagickWand_LIBRARIES}")
endif (ImageMagick_FOUND)

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
