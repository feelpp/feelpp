
foreach(TEST petsc umfpack matrix ublas vector graph backend sub_matrix)

  feelpp_add_test( ${TEST} )

endforeach()

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
