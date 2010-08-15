INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

foreach(TEST jacobi imsimplex im poly quad_order lag hermite  ravth_2D  ravth_3D ravth2 ravth1 ravth rt )

  life_add_test( ${TEST} )

endforeach()

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

