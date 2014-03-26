set_directory_properties(PROPERTIES LABEL testvf )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

feelpp_add_test( ginac )
feelpp_add_test( inner CFG test_inner.cfg )
feelpp_add_test( projtangent )
