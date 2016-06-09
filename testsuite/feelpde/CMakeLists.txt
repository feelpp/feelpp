set_directory_properties(PROPERTIES LABEL testpde )

feelpp_add_test( ams DEFS FEELPP_DIM=3 CLI --config-files cube.cfg backend_stab.cfg )
feelpp_add_test( saddle DEFS FEELPP_DIM=3 CLI --config-files cube.cfg backend_saddle.cfg )
