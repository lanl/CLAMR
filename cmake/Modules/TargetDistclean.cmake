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
   cmake_install.cmake
   */cmake_install.cmake
   */*/cmake_install.cmake
   install_manifest.txt
   */install_manifest.txt
   *.dSYM
   */*.dSYM
   */*/*.dSYM
  )
  
  ADD_CUSTOM_COMMAND(
    DEPENDS clean
    COMMENT "distribution clean"
    COMMAND rm
    ARGS    -Rf CMakeTmp ${DISTCLEANED}
    TARGET  distclean
  )
ENDIF(UNIX)
