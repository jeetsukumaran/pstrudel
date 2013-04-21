##############################################################################
## Specify the test directory
SET(TEST_DIR "${PROJECT_BINARY_DIR}/test")

##############################################################################
## Install test data
INSTALL(DIRECTORY data/
    DESTINATION ${TEST_DIR}/data)

##############################################################################
## Run test driver script as target of ``make check``
add_custom_target(check COMMAND python "${TEST_DIR}/run-tests.py")

##############################################################################
## Test executables

# file(GLOB_RECURSE DATABASE_SRCS
#          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
#          ${CMAKE_CURRENT_SOURCE_DIR}/database/*.cpp
#          ${CMAKE_CURRENT_SOURCE_DIR}/database/*.h)

FILE(GLOB_RECURSE
    PSTRUDEL_SOURCE_FILES
    "../src/pstrudel/*.cpp"
    )
MESSAGE(${PSTRUDEL_SOURCE_FILES})
list(REMOVE_ITEM PSTRUDEL_SOURCE_FILES "/Users/jeet/Documents/Projects/Phyloinformatics/pstrudel/pstrudel/test/../src/pstrudel/pstrudel-trees.cpp")

FILE(GLOB_RECURSE
    PLATYPUS_SOURCE_FILES
    "../src/platypus-phyloinformary/platypus/*.cpp"
    )

FILE(GLOB_RECURSE
    COLUGO_SOURCE_FILES
    "../src/colugo-utilities/colugo/*.cpp"
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
        ${test_src_file}
        src/testutils.cpp
        ${PSTRUDEL_SOURCE_FILES}
        ${COLUGO_SOURCE_FILES}
        ${PLATYPUS_SOURCE_FILES}
        )
    TARGET_LINK_LIBRARIES(${basename} ${NCL_LIBRARIES})
    ADD_DEPENDENCIES(check ${basename})
    INSTALL(TARGETS ${basename} DESTINATION ${TEST_DIR})
ENDFOREACH()

##############################################################################
## Install test driver script
INSTALL(PROGRAMS scripts/run-tests.py
    DESTINATION ${TEST_DIR})


