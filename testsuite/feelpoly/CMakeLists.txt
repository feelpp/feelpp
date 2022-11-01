set_directory_properties(PROPERTIES LABEL testpoly )

add_custom_target(testspoly)

foreach(TEST context_poly imsimplex im poly jacobi quad_order lag moment )# hermite ) #  raviartthomas )

  #feelpp_add_test( ${TEST} INCLUDES ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
  feelpp_add_test( ${TEST} )
  add_dependencies(testspoly feelpp_test_${TEST})
endforeach()


feelpp_add_test( integrateQuadra CFG test_integrateQuadra.cfg TIMEOUT 180)
feelpp_add_test( immultiscale CFG test_immultiscale.cfg TIMEOUT 180)

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
add_dependencies(testspoly feelpp_test_integrateQuadra feelpp_test_immultiscale )
