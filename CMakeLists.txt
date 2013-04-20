# Required CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Project name: sets ${PROJECT_NAME}
project("pstrudel" CXX)

# disallow in-source build
IF( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE )
  MESSAGE(FATAL_ERROR
"In-source builds are not allowed."
"Please create a directory and run cmake from there, passing the path to this"
"source directory as the last argument."
"This process created the file `CMakeCache.txt' and the directory `CMakeFiles'."
"Please delete them.")
ENDIF()

# Versioning
SET(PROJECT_MAJOR_VERSION 0)
SET(PROJECT_MINOR_VERSION 1)
SET(PROJECT_PATCH_VERSION 0)
SET(PROJECT_VERSION
  ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION})

FIND_PACKAGE(Git)
IF(GIT_FOUND)
  EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} status 2>/dev/null
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE GIT_STATUS_RV
      OUTPUT_QUIET
      ERROR_QUIET
      )
  IF (${GIT_STATUS_RV})
      MESSAGE(STATUS "Not a Git Repository")
  ELSE()
    EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        # chain COMMAND's in a pipe by listing them one after the other; stdout
        # of one will be passed to stdin of subsequent
        # COMMAND sed etc. etc.
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH_NAME
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_SHA1_SHORT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%cd
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DATE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    SET(GIT_SOURCE_DESC "${GIT_BRANCH_NAME}:${GIT_SHA1_SHORT}, ${GIT_DATE}")
    MESSAGE(STATUS "${GIT_SOURCE_DESC}")
  ENDIF()
ELSE()
    MESSAGE(STATUS "Git not found")
ENDIF()
configure_file("${PROJECT_SOURCE_DIR}/version.h.in" "${PROJECT_BINARY_DIR}/version.h" @ONLY)

# list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")
# include(GetGitRevisionDescription)
# get_git_head_revision(GIT_REFSPEC GIT_SHA1)
# configure_file("${PROJECT_SOURCE_DIR}/version.h.in" "${PROJECT_BINARY_DIR}/version.h" @ONLY)

# Offer the user the choice of overriding the installation directories
# set(INSTALL_LIB_DIR lib CACHE PATH "Installation sub-directory for libraries")
# set(INSTALL_BIN_DIR bin CACHE PATH "Installation sub-directory for executables")
# set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
# IF (WIN32 AND NOT CYGWIN)
#   SET(DEF_INSTALL_CMAKE_DIR CMake)
# ELSE()
#   SET(DEF_INSTALL_CMAKE_DIR lib/CMake/${PROJECT_NAME})
# ENDIF()
# SET(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# NCL (hard-coded for now)
set(NCL_INSTALL_DIR ext/v2.1/build/installed)
set(NCL_INCLUDE_DIR ext/v2.1/build/installed/include)
set(NCL_LIBRARY_DIR ext/v2.1/build/installed/lib/ncl)

# Set source directory variable
SET(MAIN_PSTRUDEL_SOURCE_DIR
    ${PROJECT_SOURCE_DIR}/pstrudel CACHE PATH "Directory of main source code")

# Set include directories
INCLUDE_DIRECTORIES(
  "${PROJECT_BINARY_DIR}"   # to find foo/config.h
  "${NCL_INCLUDE_DIR}"      # for NCL
    )

# Libraries
FIND_LIBRARY(NCL_LIBRARY
    NAMES ncl                   # name of the library as passed to '-l' flag, e.g. '-lncl'
    PATHS ${NCL_LIBRARY_DIR}    # explicit paths to search on, other than the usual
    REQUIRED                    # error if library not found
    )
IF (NOT NCL_LIBRARY )
    MESSAGE(FATAL_ERROR "Unable to locate the NCL libraries")
ENDIF()
SET(NCL_LIBRARIES ${NCL_LIBRARY} )
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
