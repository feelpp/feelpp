# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.21)

set_directory_properties(PROPERTIES LABELS "quickstart;feelsystem")

if ( NOT TARGET qs_sys )
  add_custom_target(qs_sys)
endif()
set_property(TARGET qs_sys PROPERTY LABELS feelsystem)

if ( NOT TARGET feelsystem )
  add_custom_target(feelsystem)
endif()
set_property(TARGET feelsystem PROPERTY LABELS feelsystem)

feelpp_add_application( sys_linear_2d
    SRCS qs_sys_linear.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_LINEAR_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_eigen_2d
    SRCS qs_sys_eigen.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_EIGEN_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_dirichlet_2d
    SRCS qs_sys_dirichlet.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_DIRICHLET_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_mixed_2d
    SRCS qs_sys_mixed.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_MIXED_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_model_2d
    SRCS qs_sys_model.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_MODEL_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_nonlinear_2d
    SRCS qs_sys_nonlinear.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_NONLINEAR_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_evolution_2d
    SRCS qs_sys_evolution.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=2
    EXEC QS_SYS_EVOLUTION_2D
    LABELS feelsystem
    TESTS
    INSTALL )

feelpp_add_application( sys_darcy_3d
    SRCS qs_sys_darcy.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=3
    EXEC QS_SYS_DARCY_3D
    LABELS feelsystem
    TESTS
    TEST --config-file ${CMAKE_CURRENT_SOURCE_DIR}/../laplacian/cases/quarter-turn/darcy-quarter-turn3D.cfg
         --no-export
    INSTALL )

feelpp_add_application( sys_darcy_hybrid_3d
    SRCS qs_sys_darcy_hybrid.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=3
    EXEC QS_SYS_DARCY_HYBRID_3D
    LABELS feelsystem
    TESTS
    TEST --config-file ${CMAKE_CURRENT_SOURCE_DIR}/../laplacian/cases/quarter-turn/darcy-quarter-turn3D-hybrid.cfg
         --no-export
    INSTALL )

add_custom_target(qs_sys_maxwell)
foreach(dim 2 3)
  feelpp_add_application( sys_maxwell_${dim}d
    SRCS qs_sys_maxwell.cpp
    INCLUDE_IN_ALL
    DEFS FEELPP_DIM=${dim}
    EXEC QS_SYS_MAXWELL
    LABELS "quickstart;feelsystem;maxwell;hcurl"
    TESTS
    # Match qs_maxwell: higher-order extended H(curl) ghost maps currently
    # fail in the existing form-based example as well on six MPI ranks.
    NO_MPI_TEST
    INSTALL )
  set_tests_properties(${APP_TESTS} PROPERTIES TIMEOUT 600)
  foreach(maxwell_test IN LISTS APP_TESTS)
    if(maxwell_test MATCHES "-sphere-")
      set_property(TEST ${maxwell_test} APPEND PROPERTY LABELS sphere)
    endif()
  endforeach()
  add_dependencies(qs_sys_maxwell ${QS_SYS_MAXWELL})
endforeach()

add_dependencies(qs_sys ${QS_SYS_LINEAR_2D} ${QS_SYS_EIGEN_2D} ${QS_SYS_DIRICHLET_2D} ${QS_SYS_MIXED_2D} ${QS_SYS_MODEL_2D} ${QS_SYS_NONLINEAR_2D} ${QS_SYS_EVOLUTION_2D} ${QS_SYS_DARCY_3D} ${QS_SYS_DARCY_HYBRID_3D} qs_sys_maxwell)
add_dependencies(feelsystem qs_sys)
add_dependencies(quickstart qs_sys)
