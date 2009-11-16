foreach(TEST bdf interpolation integration integration_ho integration_3d space geomap mixed)

  life_add_test( ${TEST} )

endforeach()


add_custom_command(OUTPUT geo  COMMAND ${CMAKE_COMMAND} ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/test_integration_ho_geometry.geo ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_target(data DEPENDS geo)
add_dependencies(test_integration_ho data)
