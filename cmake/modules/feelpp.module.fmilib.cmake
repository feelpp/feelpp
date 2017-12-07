find_package( FMILib )

if( FMILIB_FOUND )
  message( STATUS "[fmilib] include dir : ${FMILIB_INCLUDE_DIR}")
  message( STATUS "[fmilib] library :  ${FMILIB_LIBRARIES}" )
  set( FEELPP_HAS_FMILIB 1 )
  set(FEELPP_LIBRARIES ${FMILIB_LIBRARIES} ${FEELPP_LIBRARIES})
endif()

include_directories( ${FMILIB_INCLUDE_DIR} )
