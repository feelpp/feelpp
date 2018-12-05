# This modules looks for FMILib
# The module defines the following variables =
# FMILIB_INCLUDE_DIR - where "fmilib.h" can be found
# FMILIB_LIBRARIES   - the fmilib librarie

find_path( FMILIB_INCLUDE_DIR fmilib.h
  PATHS ${FMILIB_DIR}/include/ $ENV{FMILIB_DIR}/include/
  )

find_library( FMILIB_LIBRARIES fmilib_shared
  PATHS ${FMILIB_DIR}/lib $ENV{FMILIB_DIR}/lib/
  )

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( FMILIB REQUIRED_VARS
  FMILIB_INCLUDE_DIR
  FMILIB_LIBRARIES )
