set_directory_properties(PROPERTIES LABEL testpde )

foreach(DIM 2 3)
  feelpp_add_test( precAFP${DIM}D SRCS test_precAFP.cpp CFG precAFP${DIM}D.cfg backend.cfg precAFP${DIM}D.mod DEFS FM_DIM=${DIM} GEO precAFP${DIM}D.geo)
endforeach()
