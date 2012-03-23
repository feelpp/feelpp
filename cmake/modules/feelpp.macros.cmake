# - Find Feel

INCLUDE(ParseArguments)

macro(feelpp_add_application)

  PARSE_ARGUMENTS(FEELPP_APP
    "SRCS;LINK_LIBRARIES;CFG;GEO;LABEL"
    "NO_TEST"
    ${ARGN}
    )
  CDR(FEELPP_APP_NAME ${FEELPP_APP_DEFAULT_ARGS})

  MESSAGE("*** Arguments for Feel++ application ${FEELPP_APP_NAME}")
  MESSAGE("    Sources: ${FEELPP_APP_SRCS}")
  MESSAGE("    Link libraries: ${FEELPP_APP_LINK_LIBRARIES}")
  MESSAGE("    Cfg file: ${FEELPP_APP_CFG}")
  MESSAGE("    Geo file: ${FEELPP_APP_GEO}")

  set(execname feel_${FEELPP_APP_NAME})
  add_executable(${execname}    ${FEELPP_APP_SRCS}  )
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
        configure_file( ${cfg} ${cfg} )
          INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
#      else()
#        message(WARNING "Executable ${FEELPP_APP_NAME}: configuration file ${cfg} does not exist")
#      endif()
    endforeach()
  endif(FEELPP_APP_CFG)

  if ( FEELPP_APP_GEO )
    foreach(  geo ${FEELPP_APP_GEO} )
#      if ( EXISTS ${cfg} )
        configure_file( ${geo} ${geo} )
          INSTALL(FILES "${geo}"  DESTINATION share/feel/geo)
#      else()
#        message(WARNING "Executable ${FEELPP_APP_NAME}: configuration file ${cfg} does not exist")
#      endif()
    endforeach()
  endif(FEELPP_APP_GEO)

endmacro(feelpp_add_application)

