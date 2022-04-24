set_directory_properties(PROPERTIES LABEL testpde )

feelpp_add_test( ams DEFS FEELPP_DIM=3 CFG cube.cfg backend_stab.cfg GEO cube.geo )
feelpp_add_test( saddle DEFS FEELPP_DIM=3 CFG cube.cfg backend_saddle.cfg GEO structured_cube.geo )
