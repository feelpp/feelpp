#removed: space trace form_severaltrialtestmesh
foreach(TEST element elementminmax project bdf interpolation integration  integration_ifaces geomap mixed
    auto_order_integrate interp_twomesh submesh disc vect_comp operatorinterpolation normal3d aitken
    form_interpolation integration_opt lift
    wire_basket  stencil
    evaluator element_serialize lambda trace operatorlinearfree functionallinearfree
    lowerdim_entity integration_ginac spaces lm on_inside transfer_matrix hypercubeinterpolation
    )

  feelpp_add_test( ${TEST} )

endforeach()

feelpp_add_test( integration_ho GEO test_integration_ho_geometry.geo test_on_inside.geo)
# add_executable( test_operatorinterpolation test_operatorinterpolation.cpp )
# target_link_libraries( test_wirebasket ${FEELPP_LIBRARIES})

feelpp_add_test( mesh_transfinite CFG test_mesh_transfinite.geocfg GEO test_mesh_transfinite.geo )
