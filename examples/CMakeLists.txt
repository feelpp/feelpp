project( feel-examples )

macro(feel_add_example examplename)

  set(targetname feel_${examplename})
  set(filename ${examplename}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${FEEL_LIBRARIES} )
endmacro(feel_add_example)

foreach( DIR generic heat solid )
  add_subdirectory( ${DIR} )
endforeach()
