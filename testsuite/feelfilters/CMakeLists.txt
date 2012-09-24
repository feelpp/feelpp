
foreach(TEST importergmsh geotool)

  feelpp_add_test( ${TEST} )

endforeach()
add_dependencies(test_importergmsh feel_mesh_data)
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

