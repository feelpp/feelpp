set_directory_properties(PROPERTIES LABEL testalg )
foreach(TEST petsc umfpack ublas graph sub_matrix eigen3_tensor)
  feelpp_add_test( ${TEST}  NO_MPI_TEST )
endforeach()

feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_PETSC)
feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_SLEPC)
if(FEELPP_HAS_PETSC)
  feelpp_add_test( topetsc  NO_MPI_TEST )

  foreach(TEST vector matrix)
    feelpp_add_test( ${TEST} )
  endforeach()
endif()

feelpp_add_test( backend )

feelpp_add_test( nlsolve_eigen NO_MPI_TEST )
feelpp_add_test( matrix_block )
feelpp_add_test( prepost_solve CFG test_prepost_solve.cfg )
feelpp_add_test( add_matrix )

if ( FEELPP_HAS_SLEPC )
  feelpp_add_test( eigenmode CFG test_eigenmode.cfg )
endif()

feelpp_add_test( ptap )
feelpp_add_test( tensor )
