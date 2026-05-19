set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testalg")
foreach(TEST petsc umfpack ublas graph sub_matrix eigen3_tensor)
  feelpp_add_test( ${TEST}  NO_MPI_TEST )
endforeach()

feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_PETSC)
feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_SLEPC)
if(FEELPP_HAS_PETSC)
  feelpp_add_test( topetsc  NO_MPI_TEST )
  feelpp_add_test( hpddm NO_MPI_TEST CLI "--test-hpddm-pc-hpddm-coarse-correction=balanced --test-hpddm-pc-hpddm-share-sub-ksp=true --test-hpddm-pc-hpddm-has-neumann=true --test-hpddm-pc-hpddm-define-subdomains=true --test-hpddm-pc-hpddm-levels-1-eps-nev=1" )
  feelpp_add_test( hpddm_laplacian CLI "--gmsh.hsize=0.2 --test-hpddm-laplacian-pc-type=hpddm --test-hpddm-laplacian-pc-hpddm-has-neumann=true --test-hpddm-laplacian-pc-hpddm-define-subdomains=true --test-hpddm-laplacian-pc-hpddm-levels-1-eps-nev=1" )
  feelpp_add_test( mumps )

  foreach(TEST vector matrix matvec)
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
feelpp_add_test( functionspace_element_guard )
