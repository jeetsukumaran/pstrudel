INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_LIST_DIR}"
)

IF(NOT DEFINED ${NCL_LIBRARY})
    SET(NCL_LIBRARY, "ncl")
ENDIF()

ADD_LIBRARY (${NCL_LIBRARY}
    ncl/nxsassumptionsblock.cpp
    ncl/nxsblock.cpp
    ncl/nxscharactersblock.cpp
    ncl/nxscxxdiscretematrix.cpp
    ncl/nxsdatablock.cpp
    ncl/nxsdistancesblock.cpp
    ncl/nxsexception.cpp
    ncl/nxsmultiformat.cpp
    ncl/nxspublicblocks.cpp
    ncl/nxsreader.cpp
    ncl/nxssetreader.cpp
    ncl/nxsstring.cpp
    ncl/nxstaxablock.cpp
    ncl/nxstaxaassociationblock.cpp
    ncl/nxstoken.cpp
    ncl/nxstreesblock.cpp
    ncl/nxsunalignedblock.cpp
    )
IF (DEFINED INSTALL_NCL)
    INSTALL(TARGETS ${NCL_LIBRARY} DESTINATION lib)
ENDIF()

IF (DEFINED INSTALL_NCL)
    SET(NCL_HEADERS
        ncl/ncl.h
        ncl/nxsallocatematrix.h
        ncl/nxsassumptionsblock.h
        ncl/nxsblock.cpp
        ncl/nxsblock.h
        ncl/nxscharactersblock.h
        ncl/nxscdiscretematrix.h
        ncl/nxscxxdiscretematrix.h
        ncl/nxsdatablock.h
        ncl/nxsdefs.h
        ncl/nxsdiscretedatum.h
        ncl/nxsdistancedatum.h
        ncl/nxsdistancesblock.h
        ncl/nxsexception.h
        ncl/nxsmultiformat.h
        ncl/nxspublicblocks.h
        ncl/nxsreader.h
        ncl/nxssetreader.h
        ncl/nxsstring.h
        ncl/nxstaxablock.h
        ncl/nxstaxaassociationblock.h
        ncl/nxstoken.h
        ncl/nxstreesblock.h
        ncl/nxsunalignedblock.h
        ncl/nxsutilcopy.h
        )
    INSTALL(FILES ${NCL_HEADERS} DESTINATION include/ncl)
ENDIF()
