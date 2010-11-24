set( EXAMPLES pmesh laplacian laplacian_p laplacian_p2 advection dar  wave
    laplacian_multiplicateur_lagrange periodic stokes)

foreach( EXAMPLE ${EXAMPLES} )
  feel_add_example( ${EXAMPLE} )
endforeach()

add_executable(feel_sound sound.cpp room.cpp )
target_link_libraries(feel_sound ${FEEL_LIBRARIES} )
