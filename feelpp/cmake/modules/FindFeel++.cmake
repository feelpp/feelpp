# - Find Feel
# This module looks for Feel (Library for the Finite Element Method) support
# it will define the following values
#  FEELPP_FOUND = set if FEELPP has been found
#  FEELPP_INCLUDE_DIR = where feel/feel.hpp can be found
#  FEELPP_LIBRARY    = the library to link in

find_path(FEELPP_INCLUDE_DIR feel.hpp
    HINTS
    $ENV{FEELPP_DIR}/include/feelpp/feel
    /usr/include/feelpp/feel
    /usr/local/include/feelpp/feel
    /opt/local/include/feelpp/feel
    NO_DEFAULT_PATH
)

find_library(FEELPP_LIBRARIES feelpp
    HINTS
    $ENV{FEELPP_DIR}/lib
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Feel++ REQUIRED_VARS FEELPP_INCLUDE_DIR FEELPP_LIBRARIES )
