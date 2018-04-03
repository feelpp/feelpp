# This modules looks for OpenModelica library and headers

# The module defines the following variables :
# OMC_INCLUDE_DIR - where  opendmodelica.h can be found
# OMC_COMPILER    - binary for openmodelica compiler
# OMC_LIBRARIES   - libraries necessary to compile with omc app


find_path( OMC_INCLUDE_DIR openmodelica.h
  PATHS ${OMC_DIR}/include/ /usr/include/
  PATH_SUFFIXES omc/c
  DOC "OpendModelica include directory"
  )

find_program( OMC_COMPILER omc
  PATH ${OMC_DIR}/bin /usr/bin
  )

find_library( OMCGC_LIBRARY omcgc
  PATHS ${OMC_DIR}/lib /usr/lib/
  PATH_SUFFIXES omc x86_64-linux-gnu/omc
  )

find_library( SIMULATIONRUNTIMEC_LIBRARY SimulationRuntimeC
  PATHS ${OMC_DIR}/lib /usr/lib/
  PATH_SUFFIXES omc x86_64-linux-gnu/omc
  )

set( OMC_LIBRARIES ${SIMULATIONRUNTIMEC_LIBRARY} ${OMCGC_LIBRARY} )


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( OMC REQUIRED_VARS
  OMC_INCLUDE_DIR
  OMC_LIBRARIES )
