AM_CPPFLAGS = -std=c++11 @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

BUILT_SOURCES = pstrudel_info.h
EXTRA_DIST = echo_pstrudel_info.sh

SUBDIRS = .

COMMON_SRC = \
	colugo-utilities/src/cmdopt.cpp \
	colugo-utilities/src/cmdopt.hpp \
	colugo-utilities/src/filesys.cpp \
	colugo-utilities/src/filesys.hpp \
	colugo-utilities/src/logger.cpp \
	colugo-utilities/src/logger.hpp \
	colugo-utilities/src/textutil.cpp \
	colugo-utilities/src/textutil.hpp \
	colugo-utilities/src/utility.cpp \
	platypus-phyloinformary/src/reader.hpp \
	platypus-phyloinformary/src/ncl_reader.hpp \
	platypus-phyloinformary/src/tree.hpp \
	basic_tree.hpp \
	character.cpp \
	character.hpp \
	dataio.hpp \
	pairwise_distance_tree.cpp \
	pairwise_distance_tree.hpp \
	profile.cpp \
	profile.hpp \
	split.cpp \
	split.hpp \
	pstrudel.hpp

bin_PROGRAMS = \
	pstrudel-trees

pstrudel_trees_SOURCES = \
	$(COMMON_SRC) \
	pstrudel-trees.cpp

pstrudel_info.h: $(HEADERS) $(SOURCES)
	$(top_srcdir)/src/echo_pstrudel_info.sh > $(top_builddir)/pstrudel_info.h
