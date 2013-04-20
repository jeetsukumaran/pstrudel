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
INCLUDE(GetGitRevisionDescription)
get_git_head_revision(
    GIT_REFSPEC
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
#define PROJECT_GIT_REFSPEC         \"@GIT_REFSPEC@\"
#define PROJECT_GIT_SHA1            \"@GIT_SHA1@\"
#define PROJECT_GIT_SHA1_SHORT      \"@GIT_SHA1_SHORT@\"
#define PROJECT_GIT_BRANCH          \"@GIT_BRANCH@\"
#define PROJECT_GIT_COMMIT_DATE     \"@GIT_COMMIT_DATE@\"
#define PROJECT_SOURCE_IDENTIFIER   \"@PROJECT_SOURCE_IDENTIFIER@\"
#define PROJECT_SOURCE_IDENTIFIED      @PROJECT_SOURCE_IDENTIFIED@
")
CONFIGURE_FILE("${PROJECT_BINARY_DIR}/version.h.in" "${PROJECT_BINARY_DIR}/version.h" @ONLY)
