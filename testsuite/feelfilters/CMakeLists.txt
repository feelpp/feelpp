set_directory_properties(PROPERTIES LABEL testfilters )

foreach(TEST importergmsh importerarm geotool geotool2 exporter_sanitize exporter_disc)

  feelpp_add_test( ${TEST} )

endforeach()
IF(NOT APPLE AND EXISTS "${CMAKE_SOURCE_DIR}/data")
	add_dependencies(feelpp_test_importergmsh feelpp_mesh_data)
endif()

if(FEELPP_ENABLE_MED)
  feelpp_add_test( med SRCS test_med.cpp CFG test_med.cfg MESH test_med.med NO_MPI_TEST)
endif()
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
