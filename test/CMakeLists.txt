##############################################################################
## Specify the test directory
SET(TEST_DIR "${PROJECT_BINARY_DIR}/test")

##############################################################################
## Install test data
INSTALL(DIRECTORY data/
    DESTINATION ${TEST_DIR}/data)

##############################################################################
## Install test script
INSTALL(PROGRAMS scripts/run-tests.py
    DESTINATION ${TEST_DIR})

##############################################################################
## Add test
add_custom_target(check COMMAND python "${TEST_DIR}/run-tests.py")
