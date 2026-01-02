set_directory_properties(PROPERTIES LABELS "testsuite;testmath")

foreach(TEST jacobiellipticfunctions)

  feelpp_add_test( ${TEST} )

endforeach()

