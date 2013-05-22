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
INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}/src" # to find pstrudel/foo.hpp
)

##############################################################################
## Install test driver script and data, and run test driver as target of
## ``make check``

# ADD_CUSTOM_TARGET(check
#     COMMAND ${CMAKE_COMMAND} -D COMPONENT=testsuite -P cmake_install.cmake
#     COMMAND python "${TEST_DIR}/run-tests.py"
#     )
IF(CMAKE_CONFIGURATION_TYPES)
    ADD_CUSTOM_TARGET(check
        COMMAND ${CMAKE_COMMAND} -D COMPONENT=testsuite -P cmake_install.cmake
        COMMAND ${CMAKE_CTEST_COMMAND}
            --force-new-ctest-process --output-on-failure
            --build-config "$<CONFIGURATION>")
ELSE()
    ADD_CUSTOM_TARGET(check
        COMMAND ${CMAKE_COMMAND} -D COMPONENT=testsuite -P cmake_install.cmake
        COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure)
ENDIF()

INSTALL(PROGRAMS scripts/calc-splits-distances.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/calc-pairwise-tip-distances.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/calc-multiset-symmetric-difference.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/calc-subtree-leaf-set-sizes.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/calc-tree-pairwise-profile-distance.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/calc-tree-unlabeled-symmetric-difference.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/check-coalescent-edge-lengths.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/check-lineage-counts-through-time.py DESTINATION ${TEST_DIR} COMPONENT testsuite)
INSTALL(PROGRAMS scripts/check-lineage-spltting-times.py DESTINATION ${TEST_DIR} COMPONENT testsuite)

INSTALL(DIRECTORY data/ DESTINATION ${TEST_DIR}/data COMPONENT testsuite)

##############################################################################
## Build the main and test code base in a library, installing it as part of
## the test suite. This is to save build time for tests.

## Auto-discover the project sources. Discouraged in official documentation as
## this list will only get populated when `cmake` is run, and new files added
## to the directories will not get picked up until `cmake` is run again. I can
## live with that.
FILE(GLOB
    PROJECT_SOURCES_ALL
    "../src/pstrudel/*.cpp"
    )

## We need to filter out the actual application source (i.e., implementation
## files with a ``main()`` function. I use the convention of naming such
## files with a "_main.cpp" suffix so they can easily be identified.
INCLUDE(ListFilter)
list_filter_out(PROJECT_SOURCES "_main.cpp$" ${PROJECT_SOURCES_ALL})

## Create a library directory at configuration time.
## Since the libraries will get only installed during the second
## phase (as part of the ``testsuite`` COMPONENT), the linker
## will complain of a missing directory during the first pass
## because the "lib" directory has not been created yet. We
## avoid this annoyance by creating the directory ahead of time.
SET(LIBDIR "${PROJECT_BINARY_DIR}/lib")
FILE(MAKE_DIRECTORY ${LIBDIR})

## Build the library for the main code base, setting it up to be installed
## during the ``testsuite`` COMPONENT install.
ADD_LIBRARY(pstrudelmain ${PROJECT_SOURCES})
ADD_DEPENDENCIES(check pstrudelmain)
INSTALL(TARGETS pstrudelmain DESTINATION ${LIBDIR} COMPONENT testsuite)

## Build the library for the test code base, setting it up to be installed
## during the ``testsuite`` COMPONENT install.
ADD_LIBRARY(pstrudeltesting src/pstrudel_testing.cpp)
ADD_DEPENDENCIES(check pstrudeltest)
INSTALL(TARGETS pstrudeltesting DESTINATION ${LIBDIR} COMPONENT testsuite)

## Point the linker to the library directory.
LINK_DIRECTORIES(${LIBDIR})

##############################################################################
## Build the test executables

## The list of implementation files for the test executables that need to be
## built.
SET(TEST_TARGET_SOURCES
    src/read_trees.cpp
    src/profile_generation.cpp
    src/profile_comparison.cpp
    src/labeled_symmetric_difference.cpp
    src/pairwise_tip_distance.cpp
    src/symmetric_difference_of_multisets.cpp
    src/subtree_leaf_set_sizes.cpp
    src/lineage_counts_through_time.cpp
    src/lineage_splitting_through_time.cpp
    src/tree_pairwise_profile_distance.cpp
    src/tree_unlabeled_symmetric_difference_distance.cpp
    src/canonical_reference_trees.cpp
    )

## For each test source, create the corresponding executable, add it as ##
## dependency of the ``check`` target, and install it to the test directory.
FOREACH(test_src_file ${TEST_TARGET_SOURCES})
    get_filename_component(basename ${test_src_file} NAME_WE)
    ADD_EXECUTABLE(${basename}
        ${test_src_file}
        )
    TARGET_LINK_LIBRARIES(${basename}
        ${NCL_LIBRARIES}
        pstrudelmain
        pstrudeltesting)
    ADD_DEPENDENCIES(check ${basename})
    ADD_TEST(${basename} ${basename})
    # INSTALL(TARGETS ${basename} DESTINATION ${TEST_DIR})
ENDFOREACH()


