find_package( OMC )
if( NOT OMC_FOUND )
  return()
endif()

include_directories( ${OMC_DIR} )
