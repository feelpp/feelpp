set_directory_properties(PROPERTIES LABEL testmath )

foreach(TEST jacobiellipticfunctions)

  feelpp_add_test( ${TEST} )

endforeach()

