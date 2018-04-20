# add custom target distclean
# cleans and removes cmake generated files etc.
# Jan Woetzel 04/2003
#

IF (UNIX)
  ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
  SET(DISTCLEANED
   cmake.depends
   cmake.check_depends
   CMakeCache.txt
   */CMakeCache.txt
   cmake.check_cache
   *.cmake
   *.a
   */*.a
   Makefile
   crux/Makefile
   ezcl/Makefile
   genmalloc/Makefile
   graphics/Makefile
   hash/Makefile
   l7/Makefile
   MallocPlus/Makefile
   memstats/Makefile
   mesh/Makefile
   mesh/hsfc/Makefile
   mesh/kdtree/Makefile
   mesh/zorder/Makefile
   PowerParser/Makefile
   s7/Makefile
   timer/Makefile
   tests/Makefile
   */tests/Makefile
   core core.*
   gmon.out
   *~
   CMakeFiles
   */CMakeFiles
   */*/CMakeFiles
   */CTestTestfile.cmake
   */*/CTestTestfile.cmake
   */Testing
   */*/Testing
   Testing
   cmake_install.cmake
   */cmake_install.cmake
   */*/cmake_install.cmake
   install_manifest.txt
   */install_manifest.txt
   *.dSYM
   */*.dSYM
   */*/*.dSYM
   tests/testing
  )
  
  ADD_CUSTOM_COMMAND(
    DEPENDS clean
    COMMENT "distribution clean"
    COMMAND rm
    ARGS    -Rf CMakeTmp ${DISTCLEANED}
    TARGET  distclean
  )
ENDIF(UNIX)
