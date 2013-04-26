SET(COLUGO_SOURCES
   ../colugo-utilities/colugo/cmdopt.cpp
   ../colugo-utilities/colugo/filesys.cpp
   ../colugo-utilities/colugo/logger.cpp
   ../colugo-utilities/colugo/textutil.cpp
   ../colugo-utilities/colugo/utility.cpp)
SET(PSTRUDEL_MAIN_SOURCES
    character.cpp
    pairwise_distance_tree.cpp
    profile.cpp
    split.cpp
    symmetric_difference_tree.cpp)
ADD_EXECUTABLE(pstrudel-trees
    pstrudel_trees_main.cpp
    ${COLUGO_SOURCES}
    ${PSTRUDEL_MAIN_SOURCES}
    )
ADD_EXECUTABLE(pstrudel-trees-refdist
    pstrudel_trees_refdist_main.cpp
    ${COLUGO_SOURCES}
    ${PSTRUDEL_MAIN_SOURCES}
    )
TARGET_LINK_LIBRARIES(pstrudel-trees ${NCL_LIBRARIES})
INSTALL(TARGETS pstrudel-trees
        DESTINATION bin
        )

