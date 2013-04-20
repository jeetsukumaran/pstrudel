# This function forces a reconfigure on each Git commit
# so that the values of the variables propogated to the build
# system will be current.
#
# Derived from:
#       GetGitRevisionDescription:get_git_head_revision()
#       2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
#       http://academic.cleardefinition.com
#       Iowa State University HCI Graduate Program/VRAC
#       Copyright Iowa State University 2009-2010.
#       Distributed under the Boost Software License, Version 1.0.
#       (See accompanying file LICENSE_1_0.txt or copy at
#       http://www.boost.org/LICENSE_1_0.txt)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component(_gitdescmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)

function(track_git_revision _gitrepofound)
	set(GIT_PARENT_DIR "${CMAKE_SOURCE_DIR}")
	set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
		set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
		get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
		if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
			# We have reached the root directory, we are not in git
            set(${_gitrepofound} "GITREPO-NOTFOUND" PARENT_SCOPE)
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
    configure_file("${_gitdescmoddir}/TrackGitRevision.cmake.in"
		"${GIT_DATA}/grabRef.cmake"
		@ONLY)
	include("${GIT_DATA}/grabRef.cmake")
    set(${_gitrepofound} "${GIT_DIR}" PARENT_SCOPE)
endfunction()
