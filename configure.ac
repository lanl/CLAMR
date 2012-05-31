AC_INIT(Clamr, v2.0.0)

AC_CONFIG_AUX_DIR(./config)

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_HEADER(config.h)

AC_CONFIG_MACRO_DIR([config])

dnl
dnl Caching is usually WRONG for systems with cross-mounted file systems
dnl  the cache file may correspond to a different system.
dnl
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl
dnl

#
echo "Configuring Clamr with $ac_configure_args"
#

AC_ARG_ENABLE( debug,
 [--enable-debug                     - Compiles code with debug flags ],
 )

AC_ARG_ENABLE( bounds-check,
 [--enable-bounds-check              - Compiles code with bounds check flags ],
 )

AC_ARG_ENABLE( prof,
 [--enable-prof                      - Compiles code with profiling flags ],
 )

AC_ARG_ENABLE( mpi, 
 [--enable-mpi                       - Compiles code with mpi library ],,
 enable_mpi=yes)

AC_ARG_ENABLE( mpe, 
 [--enable-mpe                       - Compiles code with mpe library ],,
 enable_mpe=no)

AC_ARG_ENABLE( opengl, 
 [--enable-opengl                    - Compiles code with opengl library ],,
 enable_opengl=yes)

AC_ARG_ENABLE( opencl, 
 [--enable-opencl                    - Compiles code with opencl library ],,
 enable_opencl=yes)

if test "${enable_mpe}" = "yes"; then
   enable_opengl=no
fi

if test "${enable_opengl}" != "no"; then
   enable_mpe=no
fi


#AM_WITH_DMALLOC

# Checks for programs
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_MAKE_SET
dnl AC_PROG_RANLIB
AC_PROG_LIBTOOL

AC_CHECK_PROG([HAVE_ADDR2LINE], [addr2line], [yes])

# Checks for header files
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics
AC_SYS_LARGEFILE

# Checks for library functions

AC_CANONICAL_HOST
case "$host" in
    *linux*)
        host=linux
        CFLAGS="-std=gnu99"
        CXXFLAGS=""
        if test "${enable_debug}" = "yes"; then
          CFLAGS="-g -DDEBUG=1 ${CFLAGS}"
          CXXFLAGS="-g -DDEBUG=1 ${CXXFLAGS}"
        else
          CFLAGS="-g -O3 ${CFLAGS}"
          CXXFLAGS="-g -O3 ${CXXFLAGS}"
        fi

        if test "${enable_bounds_check}" = "yes"; then
          CFLAGS="-fbounds-check ${CFLAGS}"
          CXXFLAGS="-fbounds-check ${CXXFLAGS}"
        fi

        if test "${enable_prof}" = "yes"; then
          CFLAGS="-pg ${CFLAGS}"
          CXXFLAGS="-pg ${CXXFLAGS}"
        fi
        LDFLAGS="-rdynamic"
        RPATH_PREFIX="-Wl,"

        ;;
    *darwin*)
        host=darwin
        CFLAGS="-std=gnu99 -Wall -Wconversion"
        CXXFLAGS="-Wall -Wconversion"
        if test "${enable_debug}" = "yes"; then
          CFLAGS="-g -DDEBUG=1 ${CFLAGS}"
          CXXFLAGS="-g -DDEBUG=1 ${CXXFLAGS}"
        else
          CFLAGS="-g -O3 -mfpmath=sse ${CFLAGS}"
          CXXFLAGS="-g -O3 -mfpmath=sse ${CXXFLAGS}"
        fi

        if test "${enable_bounds_check}" = "yes"; then
          CFLAGS="-fbounds-check ${CFLAGS}"
          CXXFLAGS="-fbounds-check ${CXXFLAGS}"
        fi

        if test "${enable_prof}" = "yes"; then
          CFLAGS="-pg ${CFLAGS}"
          CXXFLAGS="-pg ${CXXFLAGS}"
        fi

        LDFLAGS="-rdynamic"

        RPATH_PREFIX="-Wl,"

        ;;
    *)
        echo "Unknown machine"
        exit 1
        ;;
esac

if test "${enable_mpe}" != "no"; then
  AC_PATH_PROG(MPE_BIN, clog2alog)

  if test ! -z "${MPE_BIN}"; then
    MPE_BIN=`dirname ${MPE_BIN}`
    MPEHOME=`dirname ${MPE_BIN}`
  fi

  if test -z "${MPE_BIN}"; then
    if test -z "${MPEHOME}"; then
      AC_CHECK_FILE(${HOME}/mpe/include/mpe.h, MPEHOME="${HOME}/mpe")
      AC_CHECK_FILE(${HOME}/mpe-1.9.1/include/mpe.h, MPEHOME="${HOME}/mpe-1.9.1")
    fi
  fi

  if test ! -z "${MPEHOME}"; then
    AC_CHECK_FILE(${MPEHOME}/include/mpe.h, INCLUDES="${INCLUDES} -I${MPEHOME}/include")
    MPE_LIBS="${MPE_LIBS} -L${MPEHOME}/lib -lmpe -L/usr/X11/lib -lX11 -lm"
    MPE_NOMPI_LIBS="${MPE_NOMPI_LIBS} -L${MPEHOME}/lib -lmpe_nompi -L/usr/X11/lib -lX11 -lm"
    AC_DEFINE([HAVE_MPE],"yes","Use MPE for graphics")
    AC_DEFINE(HAVE_GRAPHICS,1,[Has Real-time Graphics])
    if test "${RPATH_LIBS}x" = "x"; then
      RPATH_LIBS="${MPEHOME}/lib"
    else
      RPATH_LIBS="${RPATH_LIBS}:${MPEHOME}/lib"
    fi
  fi
fi

if test "${enable_opencl}" != "no"; then
  case "$host" in
    *linux*)
      if test "`lspci | grep "VGA compatible controller" | grep nVidia | wc -l`" != "0"; then
         found="no"
         AC_CHECK_FILE(/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/CL/cl.h, INCLUDES="${INCLUDES} -I/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc" && found="yes")
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/local/cuda/include/CL/cl.h, INCLUDES="${INCLUDES} -I/usr/local/cuda/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/cuda/sdk4.0.13/OpenCL/common/inc/CL/cl.h, INCLUDES="${INCLUDES} -I/home/cuda/sdk4.0.13/OpenCL/common/inc" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib/libOpenCL.so, OPENCL_LIBS="${OPENCL_LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
             fi
           fi
         fi
      else
         if test "`lspci | grep "VGA compatible controller" | grep ATI | wc -l`" != "0"; then 
           found="no"
           AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/include/CL/cl.h, INCLUDES="${INCLUDES} -I/home/opt/AMDAPP_v2.6/include" && found="yes")
           found="no"
           if test "$found" = "no"; then
             AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, OPENCL_LIBS="${OPENCL_LIBS} -L/home/opt/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && found="yes")
             if test "${found}" = "yes"; then
               if test "${RPATH_LIBS}x" = "x"; then
                 RPATH_LIBS="/home/opt/AMDAPP_v2.6/lib/x86_64"
               else
                 RPATH_LIBS="${RPATH_LIBS}:/home/opt/AMDAPP_v2.6/lib/x86_64"
               fi
             fi
           fi
         fi
      fi
      ;;
    *darwin*)
      OPENCL_LIBS="${OPENCL_LIBS} -framework OpenCL"
      ;;
  esac

dnl ---------------------------------------------------------------
   AC_MSG_CHECKING(for OpenCL double type)
 
   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $INCLUDES"

   AC_LANG([C])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [
       [#ifdef __APPLE_CC__block]
       [#include <OpenCL/OpenCL.h>]
       [#else]
       [#include <CL/cl.h>]
       [#endif]
      ],
      [[cl_double x;]])], 
      ac_cv_opencl_double=yes, ac_cv_opencl_double=no
   )
   if test "x$ac_cv_opencl_double" = "xyes"; then
      AC_DEFINE(HAVE_CL_DOUBLE,1,[Has OpenCL double type])
   else
      AC_DEFINE(NO_CL_DOUBLE,1,[Does not have OpenCL double type])
   fi

   CFLAGS=$SAVE_C_FLAGS

   AC_MSG_RESULT($ac_cv_opencl_double)
dnl ---------------------------------------------------------------
fi

if test "${enable_opengl}" != "no"; then
dnl ---------------------------------------------------------------
   AC_MSG_CHECKING(for OpenGL libraries)
 
   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $INCLUDES"

   AC_LANG([C])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [
       [#ifdef __APPLE_CC__]
       [#include <GLUT/glut.h>]
       [#include <OpenGL/gl.h>]
       [#include <OpenGL/glu.h>]
       [#else]
       [#include <GL/glut.h>]
       [#include <GL/gl.h>]
       [#include <GL/glu.h>]
       [#endif]
      ],
      [[glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);]])], 
      ac_cv_opengl=yes, ac_cv_opengl=no
   )
   if test "x$ac_cv_opengl" = "xyes"; then
      AC_DEFINE(HAVE_OPENGL,1,[Has OpenGL libraries])
      AC_DEFINE(HAVE_GRAPHICS,1,[Has Real-time Graphics])
      LIBS="${LIBS} -L/usr/X11/lib -lGLU -lglut -lGL"
   fi

   CFLAGS=$SAVE_C_FLAGS

   AC_MSG_RESULT($ac_cv_opengl)
dnl ---------------------------------------------------------------
fi

if test "${enable_mpi}" != "no"; then
dnl ---------------------------------------------------------------
   AC_PATH_PROG(MPI_BIN, mpicc)

  if test ! -z "${MPI_BIN}"; then
    MPIHOME=`dirname ${MPI_BIN}`
    if test "${host}" = "darwin"; then
       MPIHOME=`dirname ${MPIHOME}`
    fi
  fi

  dnl Find mpi include
  dnl ----------------
  FOUND_MPI_INCLUDE="no"
  AC_CHECK_FILE(${MPIHOME}/include/mpi.h, INCLUDES="${INCLUDES} -I${MPIHOME}/include" && FOUND_MPI_INCLUDE="yes")
  if test "${FOUND_MPI_INCLUDE}" = "no"; then
    AC_CHECK_FILE(${MPI_ROOT}/include/mpi.h, MPI_HOME=${MPI_ROOT} && INCLUDES="${INCLUDES} -I${MPI_HOME}/include" && FOUND_MPI_INCLUDE="yes")
  fi

  dnl Find mpi c++ library
  dnl ----------------
  FOUND_MPI_CPP_LIB="no"
  AC_CHECK_FILE(${MPIHOME}/lib/libmpi_cxx.so, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi_cxx" && FOUND_MPI_CPP_LIB="yes" && MPI_LIB_LOCATE="${MPIHOME}")
  if test "${FOUND_MPI_CPP_LIB}" = "no"; then
     AC_CHECK_FILE(${MPIHOME}/lib/libmpi_cxx.dylib, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi_cxx" && FOUND_MPI__CPP_LIB="yes" && MPI_LIB_LOCATE="${MPIHOME}")
  fi

  dnl Find mpi library
  dnl ----------------
  FOUND_MPI_LIB="no"
  if test "${MPI_LIB_LOCATE}x" = "${MPIHOME}x"; then
     AC_CHECK_FILE(${MPIHOME}/lib/libmpi.so, MPI_LIBS="${MPI_LIBS} -lmpi" && FOUND_MPI_LIB="yes")
     if test "${FOUND_MPI_LIB}" = "no"; then
        AC_CHECK_FILE(${MPIHOME}/lib/libmpi.dylib, MPI_LIBS="${MPI_LIBS} -lmpi" && FOUND_MPI_LIB="yes")
     fi
  else
     AC_CHECK_FILE(${MPIHOME}/lib/libmpi.so, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi" && FOUND_MPI_LIB="yes")
     if test "${FOUND_MPI_LIB}" = "no"; then
        AC_CHECK_FILE(${MPIHOME}/lib/libmpi.dylib, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi" && FOUND_MPI_LIB="yes")
     fi
  fi

 
   AC_MSG_CHECKING(MPI compile)

   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $INCLUDES"

   AC_LANG([C])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [
       [#include <mpi.h>]
      ],
      [
       [int argc;]
       [char **argv;]
       [MPI_Init(&argc, &argv);]
      ])], 
      ac_cv_mpi=yes, ac_cv_mpi=no
   )
   if test "x$ac_cv_mpi" = "xyes"; then
      dnl AC_DEFINE(HAVE_MPI,1,[Has MPI libraries])
      if test "${FOUND_MPI_LIB}" = "no"; then
         MPI_LIBS="${MPI_LIBS} -lmpi_cxx -lmpi"
      fi
   fi

   CFLAGS=$SAVE_C_FLAGS

   AC_MSG_RESULT($ac_cv_mpi)
dnl ---------------------------------------------------------------
fi


if test ! -z "${RPATH_LIBS}"; then
   LIBS="${LIBS} ${RPATH_PREFIX}-rpath ${RPATH_PREFIX}${RPATH_LIBS}"
fi

AC_SUBST(CC)
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(INCLUDES)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPE_LIBS)
AC_SUBST(MPE_NOMPI_LIBS)
AC_SUBST(OPENCL_LIBS)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([l7/Makefile])
AC_CONFIG_FILES([kdtree/Makefile])
AC_CONFIG_FILES([hsfc/Makefile])
AC_CONFIG_FILES([s7/Makefile])
AC_CONFIG_FILES([timer/Makefile])
AC_CONFIG_FILES([ezcl/Makefile])
AC_CONFIG_FILES([zorder/Makefile])
AC_CONFIG_FILES([genmalloc/Makefile])

AC_OUTPUT

