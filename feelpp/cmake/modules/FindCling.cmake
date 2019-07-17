# - Find Feel
# This module looks for Feel (Library for the Finite Element Method) support
# it will define the following values
#  CLING_FOUND = set if FEELPP has been found
#  CLING_INCLUDE_DIR = where feel/feel.hpp can be found
#  CLING_LIBRARY    = the library to link in

find_path(CLING_INCLUDE_DIR ClingOptions.h
    HINTS
    $ENV{CLING_PREFIX}/include/cling/Interpreter
    ${CLING_PREFIX}/include/cling/Interpreter
    /usr/include/cling/Interpreter
    /usr/local/include/cling/Interpreter
    /opt/local/include/cling/Interpreter
    NO_DEFAULT_PATH
)

find_library(CLING_LIBRARIES cling
    HINTS
    $ENV{CLING_PREFIX}/lib
    ${CLING_PREFIX}/lib
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    NO_DEFAULT_PATH
)

if( CLING_PREFIX )
    message( STATUS, "[cling] prefix dir: ${CLING_PREFIX}" )
endif()

find_program(CLING_BIN cling
    HINTS
    ${CLING_PREFIX}/bin
    $PATH
    CMAKE_SYSTEM_PROGRAM_PATH
)

message( STATUS "[cling] include dir: ${CLING_INCLUDE_DIR}" )
message( STATUS "[cling] libraries: ${CLING_LIBRARIES}" )
message( STATUS "[cling] binary: ${CLING_BIN}" )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Cling REQUIRED_VARS CLING_BIN )
