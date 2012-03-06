project( feel-examples )

add_custom_target( examples )
macro(feel_add_example)
  PARSE_ARGUMENTS(EXAMPLE
    "LINK_LIBRARIES;CFG;GEO"
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

  if ( EXAMPLE_GEO )
    if (EXISTS $ENV{HOME})
      # create $ENV{HOME}/feel/geo if it does not exist already: this is the geometry
      # files repository
      file(MAKE_DIRECTORY $ENV{HOME}/feel/geo)

      foreach(exp ${EXAMPLE_GEO} )
        FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${exp}")
        message(STATUS ${files})
        foreach(file ${files})
          get_filename_component(filename ${file} NAME )
          add_custom_command(TARGET ${targetname}
            COMMAND ${CMAKE_COMMAND} -E copy  ${file} $ENV{HOME}/feel/geo )
          add_dependencies(${targetname} copy_mesh_data_${filename} )
        endforeach()
      endforeach()
    endif(EXISTS $ENV{HOME})
  endif()
endmacro(feel_add_example)

foreach( DIR generic heat solid snippets )
  add_subdirectory( ${DIR} )
endforeach()
