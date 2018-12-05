set_directory_properties(PROPERTIES LABEL testvf )

feelpp_add_test( ginac )
feelpp_add_test( inner CFG test_inner.cfg )
feelpp_add_test( vf_operators )
feelpp_add_test( laplacianv CFG test_laplacianv.cfg SKIP_TEST )
feelpp_add_test( laplacian CFG test_laplacian.cfg SKIP_TEST )
feelpp_add_test( laplaciant CFG test_laplaciant.cfg SKIP_TEST )
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

feelpp_add_test( forms )

feelpp_add_test( idf_functor )
feelpp_add_test( idf2_functor )

feelpp_add_test( divergence )
feelpp_add_test(productspace_dofs)

