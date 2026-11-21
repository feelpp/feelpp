# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testpoly")

set(FEELPP_POLY_TEST_TARGETS)


# C++20 modernization tests - Phase 1-3 (foundation, policies, polynomialset)
foreach(TEST concepts meta order dynamic_order basis_dynamic)
  feelpp_add_test( ${TEST} )
  list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
endforeach()

# C++20 modernization tests - Phase 3.5+ (type-level dynamic order)
# These tests define the target API for Simplex<Dim, Dynamic>, Hypercube<Dim, Dynamic>, and Lagrange<Dynamic>
foreach(TEST simplex_dynamic hypercube_dynamic lagrange_dynamic lagrange_derivative_equivalence dynamic_order_integration pointset)
  feelpp_add_test( ${TEST} )
  list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
endforeach()

# Additional C++20 modernization tests
feelpp_add_test( policy )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( rt0_static )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( bdm_static )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( nedelec_static )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( fe_runtime_contract )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( static_dynamic_api_contract )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( polynomial_order_operations )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( referenceoperators SRCS testreferenceoperators.cpp )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( lagrangetransform SRCS testlagrangetransform.cpp )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
foreach(TEST eigen_polynomial_storage polynomial_coefficients polynomial_solve
             legendre_tabulation simplex_tabulation reference_tabulation)
  feelpp_add_test( ${TEST} )
  list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
endforeach()
feelpp_add_test( assemblystorage SRCS testassemblystorage.cpp )
target_compile_definitions(${targetname} PRIVATE NDEBUG FEELPP_STORAGE_NDEBUG=1)
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( assemblystorage_assertions SRCS testassemblystorage.cpp )
target_compile_options(${targetname} PRIVATE -O0 -g -UNDEBUG)
target_compile_definitions(${targetname} PRIVATE FEELPP_STORAGE_ASSERTIONS=1)
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( piola_transform )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( manifold_geometry )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})

feelpp_add_test( feelpoly_traits SRCS test_traits.cpp )
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})

# Original tests
foreach(TEST context_poly imsimplex im poly jacobi quad_order lag moment )# hermite ) #  raviartthomas )

  #feelpp_add_test( ${TEST} INCLUDES ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
  feelpp_add_test( ${TEST} )
  list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})

endforeach()


feelpp_add_test( integrateQuadra CFG test_integrateQuadra.cfg TIMEOUT 20)
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})
feelpp_add_test( immultiscale CFG test_immultiscale.cfg TIMEOUT 20)
list(APPEND FEELPP_POLY_TEST_TARGETS ${targetname})

if ( NOT TARGET testpoly )
  add_custom_target( testpoly )
endif()

if ( TARGET feelpp )
  add_dependencies( testpoly feelpp )
elseif ( TARGET Feelpp::feelpp )
  add_dependencies( testpoly Feelpp::feelpp )
endif()

add_dependencies( testpoly ${FEELPP_POLY_TEST_TARGETS} )
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
