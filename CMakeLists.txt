# Required CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Project name: sets ${PROJECT_NAME}
project("pstrudel" CXX)

# Versioning
SET(PROJECT_MAJOR_VERSION 0)
SET(PROJECT_MINOR_VERSION 1)
SET(PROJECT_PATCH_VERSION 0)
SET(PROJECT_VERSION
  ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION})

# disallow in-source build
IF( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE )
  MESSAGE(FATAL_ERROR
"In-source builds are not allowed."
"Please create a directory and run cmake from there, passing the path to this"
"source directory as the last argument."
"This process created the file `CMakeCache.txt' and the directory `CMakeFiles'."
"Please delete them.")
ENDIF()

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation sub-directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation sub-directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
IF (WIN32 AND NOT CYGWIN)
  SET(DEF_INSTALL_CMAKE_DIR CMake)
ELSE()
  SET(DEF_INSTALL_CMAKE_DIR lib/CMake/${PROJECT_NAME})
ENDIF()
SET(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# NCL (hard-coded for now)
set(NCL_INSTALL_DIR ext/v2.1/build/installed)
set(NCL_INCLUDE_DIR ext/v2.1/build/installed/include)
set(NCL_LIBRARY_DIR ext/v2.1/build/installed/lib/ncl)

# Set library source directory variable
SET(MAIN_PSTRUDEL_SOURCE_DIR
    ${PROJECT_SOURCE_DIR}/pstrudel CACHE PATH "Directory of main source code")

# Set include directories
INCLUDE_DIRECTORIES(
  "${PROJECT_BINARY_DIR}"   # to find foo/config.h
  "${NCL_INCLUDE_DIR}"      # for NCL
    )

# Libraries
find_library(NCL_LIBRARY ncl PATHS ${NCL_LIBRARY_DIR})
set(NCL_LIBRARIES ${NCL_LIBRARY} )
# set(NCL_LIBRARIES ${NCL_LIBRARIES} debug ${NCL_LIBRARY_DIR}/libncl.a)
# set(NCL_LIBRARIES ${NCL_LIBRARIES} optimized ${NCL_LIBRARY_DIR}/libncl.dylib)

# Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
IF (CMAKE_BUILD_TYPE MATCHES debug)
    ADD_DEFINITIONS(-O0 -Wall -g)
ELSE()
    ADD_DEFINITIONS(-O3 -Wall)
ENDIF()

# enable tests
ENABLE_TESTING ()

ADD_SUBDIRECTORY(src)
