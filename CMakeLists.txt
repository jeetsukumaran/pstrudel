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

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo)
ENDIF()

IF (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    SET(CMAKE_CXX_FLAGS, "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual
            -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op
            -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast
            -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo
            -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror
            -Wno-unused")
    IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    ENDIF()
    # SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS}  -O0 -ggdb3 -fPIC -fno-inline -femit-class-debug-always -fno-eliminate-unused-debug-types -F dwarf")
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
ENDIF()

##############################################################################
## Prerequisites
FIND_PACKAGE(NCL)
IF (NOT NCL_FOUND)
    ADD_SUBDIRECTORY(ncl)
ENDIF()
INCLUDE_DIRECTORIES(${NCL_INCLUDE_DIRS})

##############################################################################
## Global source and include directory setup
SET(MAIN_PSTRUDEL_SOURCE_DIR
    ${PROJECT_SOURCE_DIR}/pstrudel CACHE PATH "Directory of main source code")
INCLUDE_DIRECTORIES(
  "${PROJECT_BINARY_DIR}"       # to find "version.h"
  "${PROJECT_SOURCE_DIR}/src/colugo-utilities/include"
  "${PROJECT_SOURCE_DIR}/src/platypus-phyloinformary/include"
)

##############################################################################
## Sources
ADD_SUBDIRECTORY(src)

##############################################################################
## Testing
ENABLE_TESTING ()
# ``EXCLUDE_FROM_ALL``: skip from `make all` build
ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)
