INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

foreach(TEST ginac )

  feelpp_add_test( ${TEST} )

endforeach()


