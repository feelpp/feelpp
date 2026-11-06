# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES LABELS "quickstart")

add_custom_target(qs_elasticity)

add_library(qs_shell_benchmark_framework STATIC
  qs_shell_benchmark_framework.cpp)
target_link_libraries(qs_shell_benchmark_framework PUBLIC Feelpp::feelpp)
target_include_directories(qs_shell_benchmark_framework PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#
# Linear elasticity
#
feelpp_add_application( elasticity_2d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_2D INSTALL TESTS )
feelpp_add_application( elasticity_3d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_3D INSTALL TESTS )
feelpp_add_application( elasticity_shell_3d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 FEELPP_HYPERCUBE=1 EXEC QS_ELASTICITY_SHELL_3D INSTALL TESTS NO_MPI_TEST )
feelpp_add_application( elasticity_shell_q2_3d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 FEELPP_HYPERCUBE=1 FEELPP_ORDER=2 EXEC QS_ELASTICITY_SHELL_Q2_3D INSTALL TESTS NO_MPI_TEST )
feelpp_add_application( elasticity_pure_traction_2d SRCS qs_elasticity_pure_traction.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_PURE_TRACTION_2D INSTALL )
feelpp_add_application( elasticity_pure_traction_3d SRCS qs_elasticity_pure_traction.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_PURE_TRACTION_3D INSTALL  )
feelpp_add_application( shell_benchmark_info SRCS qs_shell_benchmark_info.cpp INCLUDE_IN_ALL LINK_LIBRARIES qs_shell_benchmark_framework EXEC QS_SHELL_BENCHMARK_INFO INSTALL )

feelpp_add_testcase( cases PREFIX qs_elasticity CATEGORY quickstart/elasticity )
feelpp_add_testcase( python PREFIX qs_elasticity CATEGORY quickstart/elasticity )
feelpp_add_testcase( meshes PREFIX qs_elasticity CATEGORY quickstart/elasticity )
feelpp_add_testcase( sb9 PREFIX qs_elasticity CATEGORY quickstart/elasticity )
add_custom_target(qs_add_testcase_elasticity
  DEPENDS qs_elasticity_add_testcase_cases qs_elasticity_add_testcase_python
          qs_elasticity_add_testcase_meshes qs_elasticity_add_testcase_sb9)
add_dependencies(qs_elasticity ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D}
  ${QS_ELASTICITY_SHELL_3D} ${QS_ELASTICITY_SHELL_Q2_3D}
  ${QS_ELASTICITY_PURE_TRACTION_2D} ${QS_ELASTICITY_PURE_TRACTION_3D}
  ${QS_SHELL_BENCHMARK_INFO})

#
# HDG elasticity
#
feelpp_add_application( hdg_elasticity_2d SRCS qs_hdg_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_ELASTICITY_2D INSTALL TESTS )
feelpp_add_application( hdg_elasticity_3d SRCS qs_hdg_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_ELASTICITY_3D INSTALL TESTS )

add_dependencies(qs_elasticity ${QS_HDG_ELASTICITY_2D} ${QS_HDG_ELASTICITY_3D})

add_dependencies(quickstart qs_elasticity)
add_dependencies(check qs_elasticity)

add_subdirectory( sb9 )
