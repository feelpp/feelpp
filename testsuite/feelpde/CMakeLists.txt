set_directory_properties(PROPERTIES LABEL testpde )
foreach(TEST levelset )
    feelpp_add_test( ${TEST} )
endforeach()
