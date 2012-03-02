IF ( TBB_FOUND )
  SET( OTESTS tbb )
ENDIF()
foreach(TEST traits singleton factory debug context simget ${OTESTS})

  feel_add_test( ${TEST} )

endforeach()


add_executable(test_boostmpi test_boostmpi.cpp)
target_link_libraries(test_boostmpi ${FEELPP_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
add_dependencies(testsuite test_boostmpi)
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )


