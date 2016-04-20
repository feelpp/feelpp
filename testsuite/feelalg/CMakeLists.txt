set_directory_properties(PROPERTIES LABEL testalg )

foreach(TEST petsc umfpack ublas graph sub_matrix topetsc)
  feelpp_add_test( ${TEST}  NO_MPI_TEST )
endforeach()

foreach(TEST vector matrix backend )
  feelpp_add_test( ${TEST} )
endforeach()


feelpp_add_test( nlsolve_eigen NO_MPI_TEST )
feelpp_add_test( eigenmode CLI "--gmsh.hsize 0.2 -eps_target 20 -st_pc_factor_mat_solver_package mumps" )
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

feelpp_add_test( matrix_block )
feelpp_add_test( prepost_solve CFG test_prepost_solve.cfg )
