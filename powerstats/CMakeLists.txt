cmake_minimum_required(VERSION 3.1.3)

project (POWERSTATS)
# The version number.
set (POWERSTATS_VERSION_MAJOR 2)
set (POWERSTATS_VERSION_MINOR 0.7)
set (POWERSTATS_VERSION "${POWERSTATS_VERSION_MAJOR}.${POWERSTATS_VERSION_MINOR}")

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

find_library(POWER_GADGET IntelPowerGadget)
if (POWER_GADGET)
   set (POWER_GADGET_FOUND on)
   add_definitions(-DHAVE_POWER_GADGET)
   message("-- Detected Intel Power Gadget")
   include_directories(/Library/Frameworks/IntelPowerGadget.framework/Headers)
endif (POWER_GADGET)

########### global settings ###############
#if ("${CMAKE_SYSTEM_NAME}" STREQUAL "linux")
   set(ADDITIONAL_FILES rapl_read.c)
#endif ("${CMAKE_SYSTEM_NAME}" STREQUAL "linux")

if (EXISTS /dev/cpu/0/msr)
   message("Detected msr capability")
   add_definitions(-DHAVE_RAPL_MSR)
else ()
   message("No msr capability")
endif ()

if (EXISTS /sys/class/powercap/intel-rapl/intel-rapl:0)
   message("Detected powercap capability")
   add_definitions(-DHAVE_RAPL_SYSFS)
else ()
   message("No powercap capability")
endif ()

if (EXISTS /proc/sys/kernel/perf_event_paranoid)
   message("Detected /proc/sys/kernel/perf_event_paranoid file")
   file(READ /proc/sys/kernel/perf_event_paranoid HAVE_RAPL_PERF)
   if (${HAVE_RAPL_PERF} LESS 1)
      message("Detected paranoid setting < 1")
      add_definitions(-DHAVE_RAPL_PERF)
   endif ()
else ()
   message("No /proc/sys/kernel/perf_event_paranoid file")
endif (EXISTS /proc/sys/kernel/perf_event_paranoid)

set(H_SRCS powerstats.h)
set(C_SRCS powerstats.c ${ADDITIONAL_FILES})

set(powerstats_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### powerstats target ###############
add_library(powerstats SHARED ${powerstats_LIB_SRCS})

set_target_properties(powerstats PROPERTIES VERSION ${POWERSTATS_VERSION} SOVERSION ${POWERSTATS_VERSION_MAJOR})
install(TARGETS powerstats DESTINATION lib)

if (POWER_GADGET_FOUND)
   target_link_libraries(powerstats ${POWER_GADGET})
endif (POWER_GADGET_FOUND)

########### install files ###############

install(FILES  powerstats.h DESTINATION include/powerstats)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/powerstats.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/powerstats.pdf)
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
   add_custom_target(powerstats_doc
                     COMMAND make powerstats_pdfdoc
                     COMMAND make powerstats_htmldoc
                     COMMAND make powerstats_singlehtmldoc)
   set_target_properties(powerstats_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/powerstats.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(powerstats_pdfdoc DEPENDS _build/latex/powerstats.pdf)
   set_target_properties(powerstats_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(powerstats_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(powerstats_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(powerstats_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(powerstats_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")

