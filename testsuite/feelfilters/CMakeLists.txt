set_directory_properties(PROPERTIES LABEL testfilters )
foreach(TEST importergmsh geotool geotool2 exporter_sanitize exporter_disc)

  feelpp_add_test( ${TEST} )

endforeach()
IF(NOT APPLE AND EXISTS "${CMAKE_SOURCE_DIR}/data")
	add_dependencies(feelpp_test_importergmsh feelpp_mesh_data)
endif()

if(FEELPP_ENABLE_Med)
  feelpp_add_test( med SRCS test_med.cpp CFG test_med.cfg)
  # for some reasons GEO test_med.med not working properly
  # force copy
  configure_file( test_med.med ${CMAKE_CURRENT_BINARY_DIR}/test_med.med COPYONLY )
  if(DEFINED ENV{FEELPP_WORKDIR})
    configure_file( test_med.med $ENV{FEELPP_WORKDIR}/geo/test_med.med COPYONLY )
  else(DEFINED ENV{FEELPP_WORKDIR})
    configure_file( test_med.med $ENV{HOME}/geo/test_med.med COPYONLY )
  endif(DEFINED ENV{FEELPP_WORKDIR})

endif()
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
