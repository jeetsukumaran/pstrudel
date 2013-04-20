# Get VCS info and make available to program(s)
# (program must include "version.h", which will be created in the build directory,
# and the following strings will be defined: PROJECT_GIT_REFSPEC,
# PROJECT_GIT_SHA1, PROJECT_GIT_SHA1_SHORT, PROJECT_GIT_BRANCH,
# PROJECT_GIT_COMMIT_DATE, PROJECT_SOURCE_IDENTIFIER. In addition
# 'PROJECT_SOURCE_IDENTIFIED' will evaluate to TRUE if VCS information
# was available, and FALSE if not (all previous strings will be blank in the
# latter case).
# Example of use :
#       # include <iostream>
#       # include "version.h"
#       int main() {
#           #if defined(PROJECT_SOURCE_IDENTIFIED) && PROJECT_SOURCE_IDENTIFIED
#               std::cout << " (" PROJECT_SOURCE_IDENTIFIER << ")";
#           #else
#               std::cout << " (" << __DATE__ << " " << __TIME__ ")";
#           #endif
#       }

INCLUDE(TrackGitRevision)

function(set_git_not_found _sha1 _shorthashvar _branchvar _commitdatevar)
    set(${_sha1} "GITREPO-NOTFOUND" PARENT_SCOPE)
    set(${_shorthashvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
    set(${_branchvar} "GITREPO_NOTFOUND" PARENT_SCOPE)
    set(${_commitdatevar} "GITREPO_NOTFOUND" PARENT_SCOPE)
endfunction()

function(get_git_revision_info _sha1 _shorthashvar _branchvar _commitdatevar)
    track_git_revision(GIT_FOUND)
    IF (NOT GIT_FOUND)
        set_git_not_found(_sha1 _shorthashvar _branchvar _commitdatevar)
        return()
    ENDIF()
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
            set_git_not_found(_sha1 _shorthashvar _branchvar _commitdatevar)
            return()
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
                COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_SHA1
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
                OUTPUT_VARIABLE GIT_COMMIT_DATE
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            SET(GIT_SOURCE_DESC "${GIT_BRANCH_NAME}:${GIT_SHA1_SHORT}, ${GIT_COMMIT_DATE}")
            MESSAGE(STATUS "${GIT_SOURCE_DESC}")
        ENDIF()
    ELSE()
        MESSAGE(STATUS "Git not found")
        set_git_not_found(_sha1 _shorthashvar _branchvar _commitdatevar)
        return()
    ENDIF()
    set(${_sha1} "${GIT_SHA1}" PARENT_SCOPE)
    set(${_shorthashvar} "${GIT_SHA1_SHORT}" PARENT_SCOPE)
    set(${_branchvar} "${GIT_BRANCH_NAME}" PARENT_SCOPE)
    set(${_commitdatevar} "${GIT_COMMIT_DATE}" PARENT_SCOPE)
endfunction()

function(create_CXX_version_file)
    get_git_revision_info(
        GIT_SHA1
        GIT_SHA1_SHORT
        GIT_BRANCH
        GIT_COMMIT_DATE)
    IF (GIT_SHA1)
        SET(PROJECT_SOURCE_IDENTIFIER "${GIT_BRANCH} ${GIT_SHA1_SHORT}: ${GIT_COMMIT_DATE}")
        SET(PROJECT_SOURCE_IDENTIFIED 1)
    ELSE()
        SET(PROJECT_SOURCE_IDENTIFIED 0)
    ENDIF()

    FILE(WRITE ${PROJECT_BINARY_DIR}/version.h.in
    "
    #define PROJECT_GIT_SHA1            \"@GIT_SHA1@\"
    #define PROJECT_GIT_SHA1_SHORT      \"@GIT_SHA1_SHORT@\"
    #define PROJECT_GIT_BRANCH          \"@GIT_BRANCH@\"
    #define PROJECT_GIT_COMMIT_DATE     \"@GIT_COMMIT_DATE@\"
    #define PROJECT_SOURCE_IDENTIFIER   \"@PROJECT_SOURCE_IDENTIFIER@\"
    #define PROJECT_SOURCE_IDENTIFIED      @PROJECT_SOURCE_IDENTIFIED@
    ")
    CONFIGURE_FILE("${PROJECT_BINARY_DIR}/version.h.in" "${PROJECT_BINARY_DIR}/version.h" @ONLY)
endfunction()
