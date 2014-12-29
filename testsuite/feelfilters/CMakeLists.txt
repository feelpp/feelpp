set_directory_properties(PROPERTIES LABEL testfilters )
foreach(TEST importergmsh geotool geotool2)

  feelpp_add_test( ${TEST} )

endforeach()
IF(NOT APPLE AND EXISTS "${CMAKE_SOURCE_DIR}/data")
	add_dependencies(feelpp_test_importergmsh feelpp_mesh_data)
endif()
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
