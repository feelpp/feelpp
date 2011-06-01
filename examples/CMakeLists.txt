project( feel-examples )

add_custom_target( examples )
macro(feel_add_example examplename)

  set(targetname feel_${examplename})
  set(filename ${examplename}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${FEEL_LIBRARIES} )
  add_dependencies( examples ${targetname} )
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${examplename}.cfg )
    configure_file(${examplename}.cfg ${examplename}.cfg)
  endif()
endmacro(feel_add_example)

foreach( DIR generic heat solid snippets )
  add_subdirectory( ${DIR} )
endforeach()

if (EXISTS maxwell)
  add_subdirectory(maxwell)
endif()
