# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES LABELS "quickstart;sb9")

feelpp_add_application( sb9 SRCS qs_sb9.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 FEELPP_HYPERCUBE=1 EXEC QS_SB9 INSTALL TESTS NO_MPI_TEST TIMEOUT 120 LABELS quickstart sb9 )
target_include_directories(${QS_SB9} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

feelpp_add_application( sb9_mixed_bending SRCS qs_sb9_mixed_bending.cpp INCLUDE_IN_ALL LINK_LIBRARIES qs_shell_benchmark_framework EXEC QS_SB9_MIXED_BENDING INSTALL )

feelpp_add_application( sb9_dynamic SRCS qs_sb9_dynamic.cpp INCLUDE_IN_ALL EXEC QS_SB9_DYNAMIC INSTALL )

# Retain the original minimal formulation as an explicitly buildable example.
feelpp_add_application( sb9_minimal SRCS qs_sb9_minimal.cpp EXCLUDE_FROM_ALL EXEC QS_SB9_MINIMAL )

add_dependencies(qs_elasticity ${QS_SB9} ${QS_SB9_MIXED_BENDING} ${QS_SB9_DYNAMIC})
