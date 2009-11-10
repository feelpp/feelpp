set( EXAMPLES pmesh laplacian laplacian_p laplacian_p2 advection dar stokes wave
    laplacian_multiplicateur_lagrange periodic)

foreach( EXAMPLE ${EXAMPLES} )
  life_add_example( ${EXAMPLE} )
endforeach()
