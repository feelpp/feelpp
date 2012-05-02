# - Find Feel

INCLUDE(ParseArguments)

macro(feelpp_add_application)

  PARSE_ARGUMENTS(FEELPP_APP
    "SRCS;LINK_LIBRARIES;CFG;GEO;LABEL"
    "NO_TEST;EXCLUDE_FROM_ALL"
    ${ARGN}
    )
  CAR(FEELPP_APP_NAME ${FEELPP_APP_DEFAULT_ARGS})

  set(execname feel_${FEELPP_APP_NAME})

  MESSAGE("*** Arguments for Feel++ application ${FEELPP_APP_NAME}")
  MESSAGE("    Sources: ${FEELPP_APP_SRCS}")
  MESSAGE("    Link libraries: ${FEELPP_APP_LINK_LIBRARIES}")
  MESSAGE("       Cfg file: ${FEELPP_APP_CFG}")
  MESSAGE("       Geo file: ${FEELPP_APP_GEO}")
  MESSAGE("       Exec file: ${execname}")
  MESSAGE("exclude from all: ${FEELPP_APP_EXCLUDE_FROM_ALL}")


  if ( FEELPP_APP_EXCLUDE_FROM_ALL)
    add_executable(${execname}  EXCLUDE_FROM_ALL  ${FEELPP_APP_SRCS}  )
  else()
    add_executable(${execname}  ${FEELPP_APP_SRCS}  )
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
