set_directory_properties(PROPERTIES LABEL testfilters )
foreach(TEST importergmsh geotool)

  feelpp_add_test( ${TEST} )

endforeach()
IF(NOT APPLE AND EXISTS "${CMAKE_SOURCE_DIR}/data")
	add_dependencies(test_importergmsh feel_mesh_data)
endif()
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
