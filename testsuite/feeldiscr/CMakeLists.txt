foreach(TEST elementminmax project bdf interpolation integration  integration_ifaces space geomap mixed
             auto_order_integrate interp_twomesh submesh disc vect_comp operatorinterpolation normal3d aitken form_interpolation integration_opt lift trace
             form_severaltrialtestmesh hcurl wire_basket solve_leak stencil)

  feelpp_add_test( ${TEST} )

endforeach()

feelpp_add_test( integration_ho GEO test_integration_ho_geometry.geo )


