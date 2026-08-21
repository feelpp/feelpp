set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testvf")

feelpp_add_test(
  expr
  CFG test_expr.cfg
  PYTHON test_expr.py
  CLI "--repository.prefix=${CMAKE_CURRENT_BINARY_DIR}/test_expr_repository_np_$np"
)

feelpp_add_test( ginac )
feelpp_add_test( ginac_jit_cache NO_MPI_TEST TIMEOUT 180 )
feelpp_add_test( ginac_jit_mpi TIMEOUT 180 )
feelpp_add_test( inner CFG test_inner.cfg )
feelpp_add_test( vf_operators )
feelpp_add_test( trial_test_basis )
feelpp_add_test( tensor_basis )
feelpp_add_test( voigt )
feelpp_add_test( contractions )
feelpp_add_test( voigt_elasticity )
feelpp_add_test( laplacianv CFG test_laplacianv.cfg )
feelpp_add_test( laplacian CFG test_laplacian.cfg )
feelpp_add_test( laplaciant CFG test_laplaciant.cfg )
if ( FEELPP_HAS_FFTW )
    feelpp_add_test( msi CFG test_msi.cfg LINK_LIBRARIES ${FFTW_LIBRARIES})
endif()

feelpp_add_test( inner_curl )
feelpp_add_test( inner_laplacian )
feelpp_add_test( vf_integrals )
#feelpp_add_test( vf_val )
feelpp_add_test( vf_chi )
feelpp_add_test( convolve CFG test_convolve.cfg )
feelpp_add_test( pow CFG test_pow.cfg)

feelpp_add_test( projtangent )
feelpp_add_test( shellgeometric )
feelpp_add_test( sb9bending LABEL sb9 )
feelpp_add_test( sb9pinching LABEL sb9 )
feelpp_add_test( sb9shear LABEL sb9 )
feelpp_add_test( sb9shearing LABEL sb9 )
feelpp_add_test( sb9mass LABEL sb9 NO_MPI_TEST )
feelpp_add_test( sb9stabilization LABEL sb9 )

feelpp_add_test( forms )
feelpp_add_test( forms_dirichlet CFG test_forms_dirichlet.cfg GEO ../feeldiscr/square_mmarked.geo )

feelpp_add_test( idf_functor )
feelpp_add_test( idf2_functor )

feelpp_add_test( divergence )
feelpp_add_test(productspace_dofs)
feelpp_add_test( expr_evaluate )
feelpp_add_test( symbolic_diff )
feelpp_add_test( symbolic_ad_expr )
feelpp_add_test( functionexpr )
