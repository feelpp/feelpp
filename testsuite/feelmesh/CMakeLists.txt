
foreach(TEST entity mesh regiontree mesh_codim1 kdtree P1mesh)

  feelpp_add_test( ${TEST} )

endforeach()

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )


