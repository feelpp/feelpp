set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testvf")

set(FEELPP_VF_TEST_TARGETS)

set(FEELPP_VF_SIMPLE_TESTS
    ginac
    vf_operators
    inner_curl
    inner_laplacian
    vf_integrals
    #vf_val
    vf_chi
    projtangent
    forms
    idf_functor
    idf2_functor
    localform_lowering
    divergence
    productspace_dofs
    expr_evaluate
    symbolic_diff
    functionexpr
)

foreach(t ${FEELPP_VF_SIMPLE_TESTS})
    feelpp_add_test(${t})
    list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
endforeach()

feelpp_add_test( expr CFG test_expr.cfg PYTHON test_expr.py )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})

feelpp_add_test( ginac_jit_cache NO_MPI_TEST TIMEOUT 180 )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
feelpp_add_test( inner CFG test_inner.cfg )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
feelpp_add_test( laplacianv CFG test_laplacianv.cfg )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
feelpp_add_test( laplacian CFG test_laplacian.cfg )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
feelpp_add_test( laplaciant CFG test_laplaciant.cfg )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})

if ( FEELPP_HAS_FFTW )
    feelpp_add_test( msi CFG test_msi.cfg LINK_LIBRARIES ${FFTW_LIBRARIES})
    list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
endif()

feelpp_add_test( convolve CFG test_convolve.cfg )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})
feelpp_add_test( pow CFG test_pow.cfg)
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})

if ( NOT TARGET testvf )
    add_custom_target( testvf )
endif()

feelpp_add_test( forms_dirichlet CFG test_forms_dirichlet.cfg GEO ../feeldiscr/square_mmarked.geo )
list(APPEND FEELPP_VF_TEST_TARGETS ${targetname})

if ( TARGET feelpp )
    add_dependencies( testvf feelpp )
elseif ( TARGET Feelpp::feelpp )
    add_dependencies( testvf Feelpp::feelpp )
endif()

add_dependencies( testvf ${FEELPP_VF_TEST_TARGETS} )
