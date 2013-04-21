##############################################################################
## Required CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

##############################################################################
## Find local CMake modules
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

##############################################################################
## Block builds from project directory
INCLUDE(DisableInSourceBuild)

##############################################################################
## Project name: sets ${PROJECT_NAME}
PROJECT("pstrudel" CXX)

##############################################################################
## Versioning

SET(PROJECT_MAJOR_VERSION 0)
SET(PROJECT_MINOR_VERSION 1)
SET(PROJECT_PATCH_VERSION 0)
SET(PROJECT_VERSION
  ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION})

##############################################################################
## Version Control Tracking

INCLUDE(TrackGitRevision)
create_CXX_git_version_info_file()

##############################################################################
## Installation Paths

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

##############################################################################
## Compiler Setup
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
IF (CMAKE_BUILD_TYPE MATCHES debug)
    ADD_DEFINITIONS(-O0 -Wall -g)
ELSE()
    ADD_DEFINITIONS(-O3 -Wall)
ENDIF()

##############################################################################
## Global source and include directory setup
SET(MAIN_PSTRUDEL_SOURCE_DIR
    ${PROJECT_SOURCE_DIR}/pstrudel CACHE PATH "Directory of main source code")
INCLUDE_DIRECTORIES(
  "${PROJECT_BINARY_DIR}"   # to find "version.h"
    )

##############################################################################
## Prerequisite: NCL
set(NCL_INSTALL_DIR ext/v2.1/build/installed)
set(NCL_INCLUDE_DIR ext/v2.1/build/installed/include)
set(NCL_LIBRARY_DIR ext/v2.1/build/installed/lib/ncl)
INCLUDE_DIRECTORIES(
  "${NCL_INCLUDE_DIR}"      # for NCL
    )
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

##############################################################################
## Testing
ENABLE_TESTING ()
ADD_SUBDIRECTORY(test)

##############################################################################
## Sources
ADD_SUBDIRECTORY(src)
