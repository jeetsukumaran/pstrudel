# - Returns a version string from Git
#
# These functions force a re-configure on each git commit so that you can
# trust the values of the variables in your build system.
#
#  get_git_head_revision(<refspecvar> <hashvar> [<additional arguments to git describe> ...])
#
# Returns the refspec and sha hash of the current head revision
#
#  git_describe(<var> [<additional arguments to git describe> ...])
#
# Returns the results of git describe on the source tree, and adjusting
# the output so that it tests false if an error occurs.
#
#  git_get_exact_tag(<var> [<additional arguments to git describe> ...])
#
# Returns the results of git describe --exact-match on the source tree,
# and adjusting the output so that it tests false if there was no exact
# matching tag.
#
# Requires CMake 2.6 or newer (uses the 'function' command)
#
# Original Author:
# 2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
# http://academic.cleardefinition.com
# Iowa State University HCI Graduate Program/VRAC
#
# Copyright Iowa State University 2009-2010.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

if(__get_git_revision_description)
	return()
endif()
set(__get_git_revision_description YES)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component(_gitdescmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)

function(get_git_desc_info _shorthashvar _branchvar _commitdatevar)
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
        set(${_shorthashvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
        set(${_branchvar} "GITREPO_NOTFOUND" PARENT_SCOPE)
        set(${_commitdatevar} "GITREPO_NOTFOUND" PARENT_SCOPE)
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
        set(${_shorthashvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
        set(${_branchvar} "GITREPO_NOTFOUND" PARENT_SCOPE)
        set(${_commitdatevar} "GITREPO_NOTFOUND" PARENT_SCOPE)
        return()
    ENDIF()
    set(${_shorthashvar} "${GIT_SHA1_SHORT}" PARENT_SCOPE)
    set(${_branchvar} "${GIT_BRANCH_NAME}" PARENT_SCOPE)
    set(${_commitdatevar} "${GIT_COMMIT_DATE}" PARENT_SCOPE)
endfunction()

function(get_git_head_revision _refspecvar _hashvar _shorthashvar _branchvar _commitdatevar)
	set(GIT_PARENT_DIR "${CMAKE_SOURCE_DIR}")
	set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
		set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
		get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
		if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
			# We have reached the root directory, we are not in git
            set(${_refspecvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
            set(${_hashvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
            set(${_shorthashvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
            set(${_branchvar} "GITREPO_NOTFOUND" PARENT_SCOPE)
            set(${_commitdatevar} "GITREPO_NOTFOUND" PARENT_SCOPE)
			return()
		endif()
		set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	endwhile()
	set(GIT_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/git-data")
	if(NOT EXISTS "${GIT_DATA}")
		file(MAKE_DIRECTORY "${GIT_DATA}")
	endif()

	if(NOT EXISTS "${GIT_DIR}/HEAD")
		return()
	endif()
	set(HEAD_FILE "${GIT_DATA}/HEAD")
	configure_file("${GIT_DIR}/HEAD" "${HEAD_FILE}" COPYONLY)

	configure_file("${_gitdescmoddir}/GetGitRevisionDescription.cmake.in"
		"${GIT_DATA}/grabRef.cmake"
		@ONLY)
	include("${GIT_DATA}/grabRef.cmake")

    get_git_desc_info(GIT_SHA1_SHORT GIT_BRANCH_NAME GIT_COMMIT_DATE)
    set(${_shorthashvar} "${GIT_SHA1_SHORT}" PARENT_SCOPE)
    set(${_branchvar} "${GIT_BRANCH_NAME}" PARENT_SCOPE)
    set(${_commitdatevar} "${GIT_COMMIT_DATE}" PARENT_SCOPE)
	set(${_refspecvar} "${HEAD_REF}" PARENT_SCOPE)
	set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
endfunction()

function(git_describe _var)
	if(NOT GIT_FOUND)
		find_package(Git QUIET)
	endif()
	get_git_head_revision(refspec hash)
	if(NOT GIT_FOUND)
		set(${_var} "GIT-NOTFOUND" PARENT_SCOPE)
		return()
	endif()
	if(NOT hash)
		set(${_var} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
		return()
	endif()

	# TODO sanitize
	#if((${ARGN}" MATCHES "&&") OR
	#	(ARGN MATCHES "||") OR
	#	(ARGN MATCHES "\\;"))
	#	message("Please report the following error to the project!")
	#	message(FATAL_ERROR "Looks like someone's doing something nefarious with git_describe! Passed arguments ${ARGN}")
	#endif()

	#message(STATUS "Arguments to execute_process: ${ARGN}")

	execute_process(COMMAND
		"${GIT_EXECUTABLE}"
		describe
		${hash}
		${ARGN}
		WORKING_DIRECTORY
		"${CMAKE_SOURCE_DIR}"
		RESULT_VARIABLE
		res
		OUTPUT_VARIABLE
		out
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT res EQUAL 0)
		set(out "${out}-${res}-NOTFOUND")
	endif()

	set(${_var} "${out}" PARENT_SCOPE)
endfunction()

function(git_get_exact_tag _var)
	git_describe(out --exact-match ${ARGN})
	set(${_var} "${out}" PARENT_SCOPE)
endfunction()

