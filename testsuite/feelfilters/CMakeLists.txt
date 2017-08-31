set_directory_properties(PROPERTIES LABEL testfilters )

foreach(TEST importergmsh importerarm geotool geotool2 exporter_sanitize gmsh)

  feelpp_add_test( ${TEST} )

endforeach()


if(FEELPP_ENABLE_MED)
  feelpp_add_test( med SRCS test_med.cpp CFG test_med.cfg MESH test_med.med NO_MPI_TEST)
endif()
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

feelpp_add_test( exporter_disc CLI "--exporter.element-spaces=P1")
