set_directory_properties(PROPERTIES LABEL testalg )
foreach(TEST petsc umfpack matrix ublas vector graph backend sub_matrix)

  feelpp_add_test( ${TEST}  NO_MPI_TEST )

endforeach()

feelpp_add_test( nlsolve_eigen NO_MPI_TEST )
feelpp_add_test( eigenmode )
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
