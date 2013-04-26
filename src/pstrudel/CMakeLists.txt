SET(COLUGO_SOURCES
   ../colugo-utilities/colugo/cmdopt.cpp
   ../colugo-utilities/colugo/filesys.cpp
   ../colugo-utilities/colugo/logger.cpp
   ../colugo-utilities/colugo/textutil.cpp
   ../colugo-utilities/colugo/utility.cpp)
SET(PSTRUDEL_MAIN_LIB_SOURCES
    character.cpp
    pairwise_distance_tree.cpp
    profile.cpp
    split.cpp
    symmetric_difference_tree.cpp)
SET(PSTRUDEL_MAIN_PROG_SOURCES
    pstrudel_trees_main.cpp
    pstrudel_trees_refdist_main.cpp)
FOREACH(main_prog_source ${PSTRUDEL_MAIN_PROG_SOURCES})
    get_filename_component(basename ${main_prog_source} NAME_WE)
    # STRING(REGEX REPLACE <pattern> <replacement string>
    #                  <target variable> <source string>)
    STRING(REGEX REPLACE "_main" "" tempname ${basename})
    STRING(REGEX REPLACE "_" "-" progname ${tempname})
    # MESSAGE(${progname})
    ADD_EXECUTABLE(${progname}
        ${main_prog_source}
        ${COLUGO_SOURCES}
        ${PSTRUDEL_MAIN_LIB_SOURCES}
        )
    TARGET_LINK_LIBRARIES(${progname} ${NCL_LIBRARIES})
    INSTALL(TARGETS ${progname}
            DESTINATION bin
            )
ENDFOREACH()

