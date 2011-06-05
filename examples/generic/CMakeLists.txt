set( EXAMPLES pmesh laplacian laplacian_p laplacian_p2 advection dar  wave
    laplacian_multiplicateur_lagrange periodic stokes)

foreach( EXAMPLE ${EXAMPLES} )
  feel_add_example( ${EXAMPLE} ${EXAMPLE}.cpp LINK_LIBRARIES ${FEEL_LIBRARIES} )
endforeach()

feel_add_example( sound sound.cpp room.cpp LINK_LIBRARIES ${FEEL_LIBRARIES} )

