set_directory_properties(PROPERTIES LABEL testpde )

feelpp_add_test( ams )
feelpp_add_test( saddle DEFS FEELPP_DIM=3)
