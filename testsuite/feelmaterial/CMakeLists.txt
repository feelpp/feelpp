set_directory_properties(PROPERTIES LABEL testmaterial )


foreach(TEST material)

  feelpp_add_test( ${TEST} )

endforeach()

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )
