IF( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE )
  MESSAGE(FATAL_ERROR
"
The build directory must be different from the main project source "
"directory. Please create a directory such as '${CMAKE_SOURCE_DIR}/build', "
"and run CMake from there, passing the path to this source directory as "
"the path argument. E.g.:
  $ cd ${CMAKE_SOURCE_DIR}
  $ mkdir build
  $ cd build
  $ cmake .. && make && sudo make install
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them:
  $ rm -r CMakeFiles/ CmakeCache.txt
"
  )
ENDIF()

