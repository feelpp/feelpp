set_directory_properties(PROPERTIES LABEL testmesh )
foreach(TEST entity mesh regiontree mesh_codim1 kdtree P1mesh meshmarker updatemarker partitioner_metis elementswithmarkedfaces)

  feelpp_add_test( ${TEST} )

endforeach()

feelpp_add_test( concatenate SRCS test_concatenate.cpp CFG test_concatenate.cfg GEO test_concatenate.geo )
feelpp_add_test( submesh SRCS test_submesh.cpp CFG test_submesh.cfg GEO test_submesh.geo )

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
