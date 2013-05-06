SET(PSTRUDEL_MAIN_LIB_SOURCES
    character.cpp
    distancetree.cpp
    profile.cpp
    split.cpp)
SET(PSTRUDEL_MAIN_PROG_SOURCES
    pstrudel_trees_main.cpp)
FOREACH(main_prog_source ${PSTRUDEL_MAIN_PROG_SOURCES})
    get_filename_component(basename ${main_prog_source} NAME_WE)
    # STRING(REGEX REPLACE <pattern> <replacement string>
    #                  <target variable> <source string>)
    STRING(REGEX REPLACE "_main" "" tempname ${basename})
    STRING(REGEX REPLACE "_" "-" progname ${tempname})
    # MESSAGE(${progname})
    ADD_EXECUTABLE(${progname}
        ${main_prog_source}
        ${PSTRUDEL_MAIN_LIB_SOURCES}
        )
    TARGET_LINK_LIBRARIES(${progname} ${NCL_LIBRARIES})
    INSTALL(TARGETS ${progname}
            DESTINATION bin
            )
ENDFOREACH()

