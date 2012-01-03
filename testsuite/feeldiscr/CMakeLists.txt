foreach(TEST elementminmax project bdf interpolation integration integration_ho integration_ifaces space geomap mixed
             auto_order_integrate interp_twomesh submesh disc vect_comp operatorinterpolation normal3d aitken form_interpolation integration_opt lift trace
             form_severaltrialtestmesh hcurl wire_basket)

  feel_add_test( ${TEST} )

endforeach()


add_custom_command(OUTPUT geo  COMMAND ${CMAKE_COMMAND} ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/test_integration_ho_geometry.geo ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_target(data DEPENDS geo)
add_dependencies(test_integration_ho data)


