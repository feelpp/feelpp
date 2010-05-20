INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
foreach(TEST jacobi imsimplex im poly quad_order lag hermite rt Rt)

  life_add_test( ${TEST} )

endforeach()

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

