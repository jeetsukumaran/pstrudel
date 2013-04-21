# The test harness is a single Python script which runs all tests and
# evaluates results. So we are not using (and cannot easily, without resorting
# to lots of CMake scripting) CMake's native test infrastructure.
#
# One approach would be to build and install all test targets using
# `ADD_EXECUTABLE()`, as well copy/transfer the test harness script and the
# supporting file data directories using 'INSTALL(PROGRAM, ...)` and
# `INSTALL(DIRECTORY, ...)` respectively. Test would be run by manually
# invoking the test harness script.
#
# The downside of this is that the entire test suite gets built on every make
# invocation. And it is heavy.
#
# The approach adopted here is based on the ``make check`` emulation described
# on the CMake wiki. We cannot use the ``make check`` emulation directly,
# because this requires that each test be added using the `ADD_TEST()`
# command, which means we lost the ability to use the Python script.
#
# As with the ``make check`` emulation, we create a custom target,
# "``check``", that separately triggers the test build/install/run branch of
# the build-system. We  exclude the builds in this subdirectory during other
# ``make`` invocations by passing the `EXCLUDE_FROM_ALL` option to the the
# `ADD_SUBDIRECTORY()` command in the parent listfile that includes this
# listfile.
#
#   ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)
#
# In contrast to the ``make check`` emulation described on the wiki, which
# uses ``ADD_TEST()`` to build test targets, we use ``ADD_EXECUTABLE``
# instead.  We also add each of these targets as a dependency of the ``check``
# target and use INSTALL(TARGET) to install the executable in the test
# directory.
#
# Now, if we were to run:
#
#    $ cmake -DCMAKE_INSTALL_PREFIX=${PWD}/install
#    $ make
#    $ make install
#
# the main programs will be built and installed, but not the test suite, which
# is the desired behavior.
#
# To build and install the test suite, we need to run:
#
#    $ make check
#
# This will build all the test executables and install them into the test
# directory.
#
# All well and good. But we also need to install the test harness script and
# the test data. And here is where the problem comes in. Under normal
# circumstances, we could just use the various forms of the `INSTALL()`
# command to achieve this:
#
#    INSTALL(PROGRAM scripts/run-tests.py ${TEST_DIR})
#    INSTALL(DIRECTORY data ${TEST_DIR}/data)
#
# BUT because of the ``EXCLUDE_FROM_ALL`` directive, these commands never get
# executed!
#
# It turns out that getting these "extras" (which are really not extras, but
# core components of the test system) installed in the right place at the
# right invocation is quite a hassle.
#
# We designate the `INSTALL()` commands that will install the "extras" as
# belonging to a particular COMPONENT, e.g. ``testsuite``. We then add an
# extra command to the ``check`` target definition that executes CMAKE again
# invoking and install script and targettin only ``testsuite`` components.
# Naturally, this command has to come before the command that invokes the test
# harness script.
#
# I have no idea how robust this is. But it seems to work: running ``make
# check`` builds and installs the test executables as well as the extras,
# followed by the test runner script being invoked.
#
# I also know that I am using an idiom explicitly discouraged by the
# documentation, i.e., using GLOB to populate the source lists. I think I can
# live with the caveat (which is to remember to run ``cmake ..`` if any new
# file gets added).


##############################################################################
## Specify the test directory
SET(TEST_DIR "${PROJECT_BINARY_DIR}/test")

##############################################################################
## Install test driver script and data, and run test driver as target of
## ``make check``
INSTALL(PROGRAMS scripts/run-tests.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(DIRECTORY data/ DESTINATION ${TEST_DIR}/data COMPONENT testsuite)
ADD_CUSTOM_TARGET(check
    COMMAND ${CMAKE_COMMAND} -D COMPONENT=testsuite -P cmake_install.cmake
    COMMAND python "${TEST_DIR}/run-tests.py"
    )

##############################################################################
## Build the test executables

INCLUDE(ListFilter)
FILE(GLOB_RECURSE
    PSTRUDEL_SOURCE_FILES_FULL
    "../src/pstrudel/*.cpp"
    )
list_filter_out(PSTRUDEL_SOURCE_FILES "_main.cpp$" ${PSTRUDEL_SOURCE_FILES_FULL})

FILE(GLOB_RECURSE
    PLATYPUS_SOURCE_FILES
    "../src/platypus-phyloinformary/platypus/*.cpp"
    )

FILE(GLOB_RECURSE
    COLUGO_SOURCE_FILES
    "../src/colugo-utilities/colugo/*.cpp"
    )

INCLUDE_DIRECTORIES(
  "${PROJECT_BINARY_DIR}"       # to find "version.h"
  "${PROJECT_SOURCE_DIR}/src/colugo-utilities"
  "${PROJECT_SOURCE_DIR}/src/platypus-phyloinformary"
  "${PROJECT_SOURCE_DIR}/src"
    )

SET(TEST_TARGET_SOURCES
    src/profile_distances.cpp
    src/profile_generation.cpp
    src/read_dna_sequences.cpp
    src/read_tree.cpp
    src/set_value_symmetric_difference.cpp
    src/subtree_clade_sizes.cpp
    src/subtree_leaf_set_sizes.cpp
    src/tree_child_iter.cpp
    src/tree_copy.cpp
    src/tree_internal_iter.cpp
    src/tree_leaf_iter.cpp
    src/tree_pairwise_tip_distances.cpp
    src/tree_postorder_iter.cpp
    src/tree_preorder_iter.cpp
    src/tree_reader.cpp
    src/tree_splits_distances.cpp)

MESSAGE(STATUS "Adding tests...")
FOREACH(test_src_file ${TEST_TARGET_SOURCES})
    get_filename_component(basename ${test_src_file} NAME_WE)
    MESSAGE(STATUS "Adding test: ${test_src_file} => ${basename}")
    ADD_EXECUTABLE(${basename}
        # EXCLUDE_FROM_ALL
        ${test_src_file}
        src/testutils.cpp
        ${PSTRUDEL_SOURCE_FILES}
        ${COLUGO_SOURCE_FILES}
        ${PLATYPUS_SOURCE_FILES}
        )
    TARGET_LINK_LIBRARIES(${basename} ${NCL_LIBRARIES})
    ADD_DEPENDENCIES(check ${basename})
    # ADD_TEST(${basename} ${basename})
    INSTALL(TARGETS ${basename} DESTINATION ${TEST_DIR})
ENDFOREACH()


