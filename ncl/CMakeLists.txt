IF(DEFINED ${NCL_LIBRARY_NAME})
    SET(NCL_LIBRARY ${NCL_LIBRARY_NAME})
ELSE()
    SET(NCL_LIBRARY "ncl")
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})

ADD_LIBRARY (${NCL_LIBRARY}
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsassumptionsblock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsblock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxscharactersblock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxscxxdiscretematrix.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdatablock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdistancesblock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsexception.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsmultiformat.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxspublicblocks.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsreader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxssetreader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsstring.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstaxablock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstaxaassociationblock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstoken.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstreesblock.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsunalignedblock.cpp
    )

SET(NCL_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}" PARENT_SCOPE)
SET(NCL_LIBRARIES "${NCL_LIBRARY}" PARENT_SCOPE)

OPTION(INSTALL_NCL "Install NCL" OFF)
# SET(INSTALL_NCL_TARGET "install.ncl")
# ADD_CUSTOM_TARGET(${INSTALL_NCL_TARGET}
#     COMMAND ${CMAKE_COMMAND} -D INSTALL_NCL="yes" -D COMPONENT=ncl -P cmake_install.cmake
#     COMMAND make install)
# ADD_DEPENDENCIES(${INSTALL_NCL_TARGET} ${NCL_LIBRARY})
IF (INSTALL_NCL)
    INSTALL(TARGETS ${NCL_LIBRARY} DESTINATION lib/ncl COMPONENT ncl)
    SET(NCL_HEADERS
        ${CMAKE_CURRENT_LIST_DIR}/ncl/ncl.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsallocatematrix.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsassumptionsblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsblock.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxscharactersblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxscdiscretematrix.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxscxxdiscretematrix.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdatablock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdefs.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdiscretedatum.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdistancedatum.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsdistancesblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsexception.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsmultiformat.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxspublicblocks.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsreader.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxssetreader.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsstring.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstaxablock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstaxaassociationblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstoken.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxstreesblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsunalignedblock.h
        ${CMAKE_CURRENT_LIST_DIR}/ncl/nxsutilcopy.h
        )
    INSTALL(FILES ${NCL_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ncl COMPONENT ncl)
    UNSET(INSTALL_NCL CACHE)
ENDIF()
