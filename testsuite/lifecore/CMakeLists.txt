macro(life_add_test testname)

  set(targetname test_${testname})
  set(filename test_${testname}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${LIFE_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
  add_test( ${testname} "${targetname}" )
endmacro(life_add_test)

 
foreach(TEST boostmpi traits singleton factory debug openmp context material)
  
  life_add_test( ${TEST} )

endforeach()

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )


