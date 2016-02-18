foreach(DIM 2 3)
  feelpp_add_application( precAFP${DIM}D SRCS test_precAFP.cpp CFG precAFP${DIM}D_sin.cfg precAFP${DIM}D_poly.cfg backend.cfg precAFP${DIM}D.mod DEFS FEELPP_DIM=${DIM} GEO precAFP${DIM}D.geo)
  feelpp_add_application( regul${DIM}D SRCS test_regul.cpp CFG regul.cfg DEFS FEELPP_DIM=${DIM})
endforeach()

add_subdirectory( regul )
add_subdirectory( saddle )

