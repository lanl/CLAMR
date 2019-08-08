cmake_minimum_required(VERSION 3.1.3)

project (GRAPHICS)
# The version number.
set (GRAPHICS_VERSION_MAJOR 2)
set (GRAPHICS_VERSION_MINOR 0.7)
set (GRAPHICS_VERSION "${GRAPHICS_VERSION_MAJOR}.${GRAPHICS_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

if (NOT CMAKE_BUILD_TYPE)
   #set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

set (CMAKE_C_STANDARD 99)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -fstrict-aliasing -ftree-vectorize")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fopt-info-vec-optimized")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "4.9")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fopenmp-simd")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "4.9")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -mtune=native")
   elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mcpu=powerpc64le")
   else ()
     set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
   endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3 -ansi-alias -vecabi=cmdtarget -xHost -qopenmp-simd -qopt-report=5 -qopt-report-phase=vec,loop")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_C_FLAGS_RELEASE "-qopt-zmm-usage=high ${CMAKE_C_FLAGS_RELEASE}")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
endif()

if (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native -Wa,-q")
else (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
endif (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pedantic")

if (NOT DEFINED OpenCL_FOUND)
   find_package(OpenCL)
   #  OpenCL_FOUND          - True if OpenCL was found
   #  OpenCL_INCLUDE_DIRS   - include directories for OpenCL
   #  OpenCL_LIBRARIES      - link against this library to use OpenCL
   #  OpenCL_VERSION_STRING - Highest supported OpenCL version (eg. 1.2)
   #  OpenCL_VERSION_MAJOR  - The major version of the OpenCL implementation
   #  OpenCL_VERSION_MINOR  - The minor version of the OpenCL implementation
   if (OpenCL_FOUND)
      set(HAVE_CL_DOUBLE ON CACHE BOOL "Have OpenCL Double")
      set(NO_CL_DOUBLE OFF)
      include_directories(${OpenCL_INCLUDE_DIRS})
   endif (OpenCL_FOUND)
endif (NOT DEFINED OpenCL_FOUND)

if (NOT DEFINED MPI_FOUND)
   find_package(MPI)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if (NOT DEFINED OPENMP_FOUND)
   find_package( OpenMP)
endif (NOT DEFINED OPENMP_FOUND)

if (NOT PRECISION_TYPE)
   set (PRECISION_TYPE full_precision)
endif (NOT PRECISION_TYPE)

if (PRECISION_TYPE STREQUAL "mixed_precision")
   set (MIXED_PRECISION on)
   set (PRECISION_TYPE mixed_precision)
   add_definitions(-DMIXED_PRECISION)
endif (PRECISION_TYPE STREQUAL "mixed_precision")

if (PRECISION_TYPE STREQUAL "full_precision")
   set (FULL_PRECISION on)
   set (PRECISION_TYPE full_precision)
   add_definitions(-DFULL_PRECISION)
endif (PRECISION_TYPE STREQUAL "full_precision")

if (PRECISION_TYPE STREQUAL "minimum_precision")
   set (MINIMUM_PRECISION on)
   set (PRECISION_TYPE minimum_precision)
   add_definitions(-DMINIMUM_PRECISION)
endif (PRECISION_TYPE STREQUAL "minimum_precision")

if (NOT GRAPHICS_TYPE)
   set (GRAPHICS_TYPE none)
endif (NOT GRAPHICS_TYPE)

if (GRAPHICS_TYPE STREQUAL "OpenGL")
   find_package(OpenGL)
   if (OPENGL_FOUND)
      find_package(GLUT)
   endif (OPENGL_FOUND)
endif (GRAPHICS_TYPE STREQUAL "OpenGL")

# If autodetect does not work try setting MPEHOME or
# MPE_INCLUDE_DIR and MPE_LIBRARIES
if (GRAPHICS_TYPE STREQUAL "MPE")
   find_package(MPE)
   find_package(X11)
endif (GRAPHICS_TYPE STREQUAL "MPE")

#message("MPI_LIBRARIES ${MPI_LIBRARIES}")
#message("MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH}")
#message("GRAPHICS_TYPE is ${GRAPHICS_TYPE}")
#message("MPE_INCLUDE is ${MPE_INCLUDE_DIR}")

if (OPENGL_FOUND AND GLUT_FOUND)
   set (HAVE_GRAPHICS on)
   set (HAVE_OPENGL on)
   include_directories(${OPENGL_INCLUDE_DIR})
   include_directories(${GLUT_INCLUDE_DIR})
endif (OPENGL_FOUND AND GLUT_FOUND)

if (NOT OPENGL_FOUND)
   if (MPE_FOUND)
      set (HAVE_GRAPHICS on)
      set (HAVE_MPE on)
      include_directories(${MPE_INCLUDE_DIR})
      if (X11_FOUND)
         include_directories(${X11_INCLUDE_DIR})
      endif (X11_FOUND)
      set (MPE_LIBS ${MPE_LIBRARIES})
      set (MPE_NOMPI_LIBS ${MPE_NOMPI_LIBRARIES})
      set (X11_LIBS ${X11_LIBRARIES})
   endif (MPE_FOUND)
endif (NOT OPENGL_FOUND)

########### global settings ###############
set(H_SRCS graphics.h display.h)
set(C_SRCS graphics.c display.c)

set(graphics_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### graphics target ###############

add_library(graphics SHARED ${graphics_LIB_SRCS})

set_target_properties(graphics PROPERTIES VERSION ${GRAPHICS_VERSION} SOVERSION ${GRAPHICS_VERSION_MAJOR})

target_link_libraries(graphics ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_link_libraries(graphics ${MPE_NOMPI_LIBS} ${X11_LIBS})
if (ImageMagick_MagickWand_FOUND)
   target_link_libraries(graphics ${ImageMagick_LIBRARIES})
endif (ImageMagick_MagickWand_FOUND)

install(TARGETS graphics DESTINATION lib)

########### pgraphics target ###############
if (MPI_FOUND)
   add_library(pgraphics SHARED ${graphics_LIB_SRCS})

   set_target_properties(pgraphics PROPERTIES VERSION ${GRAPHICS_VERSION} SOVERSION ${GRAPHICS_VERSION_MAJOR})
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

install(FILES  display.h DESTINATION include/graphics)
install(FILES  graphics.h DESTINATION include/graphics)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/Graphics.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/Graphics.pdf)
   set_target_properties(pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(Graphics_doc
                     COMMAND make Graphics_pdfdoc
                     COMMAND make Graphics_htmldoc
                     COMMAND make Graphics_singlehtmldoc)
   set_target_properties(Graphics_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/Graphics.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(Graphics_pdfdoc DEPENDS _build/latex/Graphics.pdf)
   set_target_properties(Graphics_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(Graphics_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(Graphics_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(Graphics_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(Graphics_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")

