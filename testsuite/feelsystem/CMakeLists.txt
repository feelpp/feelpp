# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testsystem;feelsystem")

if ( NOT TARGET feelsystem )
  add_custom_target(feelsystem)
endif()
if ( NOT TARGET testsuite-feelsystem )
  add_custom_target(testsuite-feelsystem)
endif()
add_dependencies(feelsystem testsuite-feelsystem)
if ( TARGET testsuite-feelpp )
  add_dependencies(testsuite-feelpp testsuite-feelsystem)
endif()
set_property(TARGET feelsystem PROPERTY LABELS feelsystem)
set_property(TARGET testsuite-feelsystem PROPERTY LABELS "testsuite;testsuite-feelpp;testsystem;feelsystem")
set(FEELPP_TEST_AGGREGATE_TARGET testsuite-feelsystem)

feelpp_add_test( system )
feelpp_add_test( system_orders SRCS test_orders.cpp )
feelpp_add_test( system_auxiliary SRCS test_auxiliary.cpp
  CLI "--auxiliary-geneo.pc-type=hpddm --auxiliary-geneo.ksp-type=cg --auxiliary-geneo.pc-hpddm-has-neumann=true --auxiliary-geneo.pc-hpddm-levels-1-st-pc-type=mat --auxiliary-geneo.pc-hpddm-levels-1-eps-nev=2" )
feelpp_add_test( system_block SRCS test_block.cpp )
feelpp_add_test( system_hdg_spaces SRCS test_hdg_spaces.cpp )
feelpp_add_test( system_hdg SRCS test_hdg.cpp )
feelpp_add_test( system_eigen SRCS test_eigen.cpp )
feelpp_add_test( system_constraints SRCS test_constraints.cpp )
feelpp_add_test( system_model_iteration SRCS test_model_iteration.cpp )
feelpp_add_test( system_nonlinear SRCS test_nonlinear.cpp )
feelpp_add_test( system_evolution SRCS test_evolution.cpp )
feelpp_add_test( system_lifetime SRCS test_lifetime.cpp )
feelpp_add_test( system_metadata SRCS test_metadata.cpp )
