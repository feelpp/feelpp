set_directory_properties(PROPERTIES LABEL testpoly )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

foreach(TEST jacobi imsimplex im poly quad_order lag moment )# hermite ) #  raviartthomas )

  feelpp_add_test( ${TEST} )

endforeach()

feelpp_add_test( integrateQuadra CFG test_integrateQuadra.cfg TIMEOUT 180)
feelpp_add_test( immultiscale CFG test_immultiscale.cfg TIMEOUT 180 SKIP_TEST)
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
