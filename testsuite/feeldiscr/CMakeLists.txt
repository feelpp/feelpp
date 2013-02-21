#removed: space trace form_severaltrialtestmesh 
foreach(TEST elementminmax project bdf interpolation integration  integration_ifaces geomap mixed
    auto_order_integrate interp_twomesh submesh disc vect_comp operatorinterpolation normal3d aitken
    form_interpolation integration_opt lift 
    hcurl wire_basket solve_leak stencil
    evaluator element_serialize lambda )

  feelpp_add_test( ${TEST} )

endforeach()

feelpp_add_test( integration_ho GEO test_integration_ho_geometry.geo )
# add_executable( test_wirebasket test_wirebasket.cpp )
# target_link_libraries( test_wirebasket ${FEELPP_LIBRARIES})
