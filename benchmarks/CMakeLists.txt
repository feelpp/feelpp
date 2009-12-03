project( life-benchmarks )

macro(life_add_bench examplename)

  set(targetname life_${examplename})
  set(filename ${examplename}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${LIFE_LIBRARIES} )
endmacro(life_add_bench)

foreach( DIR convergence perf  )
  add_subdirectory( ${DIR} )
endforeach()
