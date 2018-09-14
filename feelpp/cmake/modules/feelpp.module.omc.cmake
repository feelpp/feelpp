find_package( OMC )

if( OMC_FOUND )
  include_directories( ${OMC_INCLUDE_DIR} )
  message( STATUS "[omc] include dir : ${OMC_INCLUDE_DIR}" )
  message( STATUS "[omc] omc compiler : ${OMC_COMPILER}" )
  message( STATUS "[omc] libomcgc : " ${OMCGC_LIBRARY} )
  message( STATUS "[omc] libSimulationRuntimeC : " ${SIMULATIONRUNTIMEC_LIBRARY} )
  set( FEELPP_HAS_OMC 1 )
  set(FEELPP_LIBRARIES ${OMC_LIBRARIES} ${FEELPP_LIBRARIES})
else()
  message(WARNING "OpenModelica was not found on your system. Either install it or set FEELPP_ENABLE_OMC to OFF.")
  return()
endif()


