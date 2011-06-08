project( feel-examples )

add_custom_target( examples )
macro(feel_add_example)
  PARSE_ARGUMENTS(EXAMPLE
    "LINK_LIBRARIES;CFG"
    ""
    ${ARGN}
    )
  CAR(EXAMPLE_NAME ${EXAMPLE_DEFAULT_ARGS})
  CDR(EXAMPLE_SOURCES ${EXAMPLE_DEFAULT_ARGS})

  set(targetname feel_${EXAMPLE_NAME})
  add_executable(${targetname} ${EXAMPLE_SOURCES})
  target_link_libraries(${targetname} ${EXAMPLE_LINK_LIBRARIES} )
  set_property(TARGET ${targetname} PROPERTY LABELS examples)
  add_dependencies( examples ${targetname} )
  add_test(NAME ${targetname} COMMAND ${targetname} )
  set_property(TEST ${targetname} PROPERTY LABELS examples)
  if ( EXAMPLE_CFG )
    foreach(  cfg ${EXAMPLE_CFG} )
#      if ( EXISTS ${cfg} )
        configure_file( ${cfg} ${cfg} )
#      else()
#        message(WARNING "Examples subproject: config file ${cfg} does not exist" )
#      endif()
    endforeach()
  endif()
endmacro(feel_add_example)

foreach( DIR generic heat solid snippets )
  add_subdirectory( ${DIR} )
endforeach()
