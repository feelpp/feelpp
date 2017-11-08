# This modules looks for OpenModelica library and headers

# The module defines the following variables :
# OMC_DIR      - where  opendmodelica.h can be found
# OMC_COMPILER - binary for openmodelica compiler
# OMC_FOUND    - if openmodelica has been found

find_path( OMC_DIR openmodelica.h
  PATHS /usr/include/
  PATH_SUFFIXES omc/c
  DOC "OpendModelica include directory"
  )

find_program( OMC_COMPILER omc
  PATH $PATH
  )

if ( OMC_INCLUDE_DIR AND OMC_COMPILER)
  set( OMC_FOUND TRUE )
endif()

message( STATUS "[omc] include dir : ${OMC_DIR}" )
message( STATUS "[omc] omc compiler : ${OMC_COMPILER}" )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( OMC REQUIRED_VARS OMC_DIR )
