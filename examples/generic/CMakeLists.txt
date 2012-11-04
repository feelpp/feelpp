# laplacian_p laplacian_p2
set( EXAMPLES pmesh laplacian advection dar  wave
    laplacian_lagrange_multiplier )

foreach( EXAMPLE ${EXAMPLES} )
  if ( EXISTS ${EXAMPLE}.cfg )
    feel_add_example( ${EXAMPLE} ${EXAMPLE}.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG ${EXAMPLE}.cfg )
  else()
    feel_add_example( ${EXAMPLE} ${EXAMPLE}.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES}  )
  endif()
endforeach()

feel_add_example( sound sound.cpp room.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} )
feel_add_example( soundTime soundTime.cpp room.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} )

