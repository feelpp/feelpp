foreach(DIM 2 3)
  feelpp_add_test( precAFP${DIM}D SRCS test_precAFP.cpp CFG precAFP${DIM}D_sin.cfg precAFP${DIM}D_poly.cfg backend.cfg precAFP${DIM}D.mod DEFS FEELPP_DIM=${DIM} GEO precAFP${DIM}D.geo)
endforeach()
