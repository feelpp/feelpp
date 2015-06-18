set_directory_properties(PROPERTIES LABEL testpde )
foreach(TEST levelset )
    feelpp_add_test( ${TEST} )
endforeach()

feelpp_add_test( precAFP CFG precAFP.cfg precAFP.mod DEFS FM_DIM=2 GEO precAFP2D.geo)
