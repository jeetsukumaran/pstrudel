AM_CPPFLAGS = -std=c++11 @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

# AUTOMAKE_OPTIONS = subdir-objects

pstrudeltestsdir = $(prefix)/opt/@PACKAGE@/test
# pstrudeltestsdir = $(builddir)/bin

dist_pstrudeltests_SCRIPTS = \
	scripts/run-tests.py

nobase_dist_pstrudeltests_DATA = \
	data/general/bird_orders.nex \
	data/general/primates.beast.mcct.medianh.newick.tre \
	data/general/primates.beast.mcct.medianh.tre \
	data/general/primates.chars.fasta \
	data/general/primates.chars.nexus \
	data/general/pythonidae.chars.fasta \
	data/general/pythonidae.beast.summary.tre \
	data/general/pythonidae.postorder.newick \
	data/general/pythonidae.preorder.newick \
	data/general/smdist.tre \
	data/general/pythonidae.tree.newick \
	data/general/pythonidae.mlboots.nexus.tre \
	data/general/primates.beast.mcmc.trees \
	data/general/pythonidae.reference-trees.nexus

COMMON_TREE_SRC = \
	../src/colugo-utilities/src/cmdopt.cpp \
	../src/colugo-utilities/src/cmdopt.hpp \
	../src/colugo-utilities/src/filesys.cpp \
	../src/colugo-utilities/src/filesys.hpp \
	../src/colugo-utilities/src/logger.cpp \
	../src/colugo-utilities/src/logger.hpp \
	../src/colugo-utilities/src/textutil.cpp \
	../src/colugo-utilities/src/textutil.hpp \
	../src/colugo-utilities/src/utility.cpp \
	../src/platypus-phyloinformary/src/tree.hpp \
	../src/platypus-phyloinformary/src/reader.hpp \
	../src/platypus-phyloinformary/src/ncl_reader.hpp \
	../src/basic_tree.hpp \
	../src/character.cpp \
	../src/pairwise_distance_tree.hpp \
	../src/pairwise_distance_tree.cpp \
	../src/profile.hpp \
	../src/profile.cpp \
	../src/split.hpp \
	../src/split.cpp \
	../src/symmetric_difference_tree.hpp \
	../src/symmetric_difference_tree.cpp \
	../src/dataio.hpp

COMMON_TEST_SRC = \
	src/testutils.hpp \
	src/testutils.cpp

# pstrudeltests_PROGRAMS = \
# 	read_tree

# read_tree_SOURCES = \
# 	$(COMMON_TREE_SRC) \
# 	$(COMMON_TEST_SRC) \
# 	src/read_tree.cpp

pstrudeltests_PROGRAMS = \
	tree_postorder_iter \
	tree_preorder_iter \
	tree_leaf_iter \
	tree_internal_iter \
	tree_child_iter \
	tree_reader \
	tree_copy \
	read_tree \
	read_dna_sequences \
	tree_pairwise_tip_distances \
	profile_generation \
	profile_distances \
	subtree_leaf_set_sizes \
	subtree_clade_sizes \
	tree_splits_distances \
	set_value_symmetric_difference

tree_postorder_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_postorder_iter.cpp

tree_preorder_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_preorder_iter.cpp

tree_leaf_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_leaf_iter.cpp

tree_internal_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_internal_iter.cpp

tree_child_iter_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_child_iter.cpp

tree_reader_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_reader.cpp

tree_copy_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_copy.cpp

read_tree_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/read_tree.cpp

read_dna_sequences_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/read_dna_sequences.cpp

tree_pairwise_tip_distances_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_pairwise_tip_distances.cpp

profile_generation_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/profile_generation.cpp

profile_distances_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/profile_distances.cpp

tree_splits_distances_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/tree_splits_distances.cpp

subtree_leaf_set_sizes_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/subtree_leaf_set_sizes.cpp

subtree_clade_sizes_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/subtree_clade_sizes.cpp

set_value_symmetric_difference_SOURCES = \
	$(COMMON_TREE_SRC) \
	$(COMMON_TEST_SRC) \
	src/set_value_symmetric_difference.cpp

