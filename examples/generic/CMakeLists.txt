set( EXAMPLES pmesh laplacian laplacian_p laplacian_p2 advection dar  wave
    laplacian_multiplicateur_lagrange periodic)

foreach( EXAMPLE ${EXAMPLES} )
  life_add_example( ${EXAMPLE} )
endforeach()

add_executable(life_sound sound.cpp room.cpp )
target_link_libraries(life_sound ${LIFE_LIBRARIES} )
