# - Find Feel

INCLUDE(ParseArguments)

macro(feelpp_add_application)

  PARSE_ARGUMENTS(FEELPP_APP
    "SRCS;LINK_LIBRARIES;CFG;GEO;LABEL;DEFS;DEPS"
    "NO_TEST;EXCLUDE_FROM_ALL"
    ${ARGN}
    )
  CAR(FEELPP_APP_NAME ${FEELPP_APP_DEFAULT_ARGS})

  set(execname feel_${FEELPP_APP_NAME})

  MESSAGE("*** Arguments for Feel++ application ${FEELPP_APP_NAME}")
  MESSAGE("    Sources: ${FEELPP_APP_SRCS}")
  MESSAGE("    Link libraries: ${FEELPP_APP_LINK_LIBRARIES}")
  MESSAGE("       Cfg file: ${FEELPP_APP_CFG}")
  MESSAGE("      Deps file: ${FEELPP_APP_DEPS}")
  MESSAGE("      Defs file: ${FEELPP_APP_DEFS}")
  MESSAGE("       Geo file: ${FEELPP_APP_GEO}")
  MESSAGE("       Exec file: ${execname}")
  MESSAGE("exclude from all: ${FEELPP_APP_EXCLUDE_FROM_ALL}")


  if ( FEELPP_APP_EXCLUDE_FROM_ALL)
    add_executable(${execname}  EXCLUDE_FROM_ALL  ${FEELPP_APP_SRCS}  )
  else()
    add_executable(${execname}  ${FEELPP_APP_SRCS}  )
  endif()
  if ( FEELPP_APP_DEPS )
    add_dependencies(${execname} ${FEELPP_APP_DEPS})
  endif()
  if ( FEELPP_APP_DEFS )
    set_property(TARGET ${execname} PROPERTY COMPILE_DEFINITIONS ${FEELPP_APP_DEFS})
  endif()
  target_link_libraries( ${execname} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})
  INSTALL(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${execname}"  DESTINATION bin COMPONENT Bin)
  add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
  #add_dependencies(crb ${execname})
  # Add label if provided
  if ( FEELPP_APP_LABEL )
    set_property(TARGET ${execname} PROPERTY LABELS ${FEELPP_APP_LABEL})
    set_property(TEST ${execname} PROPERTY LABELS ${FEELPP_APP_LABEL})
  endif()



  if ( FEELPP_APP_CFG )
    foreach(  cfg ${FEELPP_APP_CFG} )
      #      if ( EXISTS ${cfg} )
      # extract cfg filename  to be copied in binary dir
      get_filename_component( CFG_NAME ${cfg} NAME )
      configure_file( ${cfg} ${CFG_NAME} )
      INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
      #      else()
      #        message(WARNING "Executable ${FEELPP_APP_NAME}: configuration file ${cfg} does not exist")
      #      endif()
    endforeach()
  endif(FEELPP_APP_CFG)

  if ( FEELPP_APP_GEO )
    foreach(  geo ${FEELPP_APP_GEO} )
      # extract geo filename  to be copied in binary dir
      get_filename_component( GEO_NAME ${geo} NAME )
      configure_file( ${geo} ${GEO_NAME} )
      INSTALL(FILES "${geo}"  DESTINATION share/feel/geo)
    endforeach()
  endif(FEELPP_APP_GEO)

endmacro(feelpp_add_application)


macro(OVERWITE_IF_DIFFERENT thetarget filename var dummy)
  if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
    # be careful if file does not exist we use dummy to generate the cpp file which will
    # then be overwritten using the cmake -E copy_if_different command
    configure_file(${dummy}  ${CMAKE_CURRENT_BINARY_DIR}/${filename})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/copy_${filename} ${var})
  add_custom_command(TARGET ${thetarget} COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/copy_${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endmacro()

macro(feelpp_add_test)
  PARSE_ARGUMENTS(FEELPP_TEST
    "SRCS;LINK_LIBRARIES;CFG;GEO;LABEL;DEFS;DEPS"
    "NO_TEST;EXCLUDE_FROM_ALL"
    ${ARGN}
    )
  CAR(FEELPP_TEST_NAME ${FEELPP_TEST_DEFAULT_ARGS})

  if ( NOT FEELPP_TEST_SRCS )
    set(targetname test_${FEELPP_TEST_NAME})
    set(filename test_${FEELPP_TEST_NAME}.cpp)
    add_executable(${targetname} ${filename})
    target_link_libraries(${targetname} ${FEELPP_LIBRARIES} ${FEELPP_TEST_LINK_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
    set_property(TARGET ${targetname} PROPERTY LABELS testsuite)
    add_dependencies(testsuite ${targetname})

    add_test(
      NAME test_${FEELPP_TEST_NAME}
      COMMAND ${targetname} --log_level=message
      )

    set_property(TEST ${targetname} PROPERTY LABELS testsuite)

    if ( FEELPP_TEST_GEO )
      foreach(  geo ${FEELPP_TEST_GEO} )
        # extract geo filename  to be copied in binary dir
        get_filename_component( GEO_NAME ${geo} NAME )
        configure_file( ${geo} ${GEO_NAME} )
      endforeach()
    endif(FEELPP_TEST_GEO)

  endif()


endmacro(feelpp_add_test)

