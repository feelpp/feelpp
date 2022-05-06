set_directory_properties(PROPERTIES LABEL testfilters )

foreach(TEST importergmsh importerarm geotool geotool2 exporter_sanitize gmsh importer_mesh hbf)

  feelpp_add_test( ${TEST} )

endforeach()


if(FEELPP_HAS_MED)
  message(STATUS "[testsuite] feelpp_test_med enabled" )
  feelpp_add_test( med SRCS test_med.cpp CFG test_med.cfg NO_MPI_TEST)
endif()

feelpp_add_test( exporter_disc )
feelpp_add_test( read_write_data )
