set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;testmath")

foreach(TEST jacobiellipticfunctions)

  feelpp_add_test( ${TEST} )

endforeach()

