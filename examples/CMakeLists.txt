project( life-examples )

macro(life_add_example examplename)

  set(targetname life_${examplename})
  set(filename ${examplename}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${LIFE_LIBRARIES} )
endmacro(life_add_example)

add_subdirectory( generic )
