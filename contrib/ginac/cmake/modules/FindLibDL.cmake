# - Find libdl
# Find the native LIBDL includes and library
#
#  LIBDL_INCLUDE_DIR - where to find dlfcn.h, etc.
#  LIBDL_LIBRARIES   - List of libraries when using libdl.
#  LIBDL_FOUND       - True if libdl found.


IF (LIBDL_INCLUDE_DIR)
	# Already in cache, be silent
	SET(LIBDL_FIND_QUIETLY TRUE)
ENDIF (LIBDL_INCLUDE_DIR)

FIND_PATH(LIBDL_INCLUDE_DIR dlfcn.h)

SET(LIBDL_NAMES dl libdl ltdl libltdl)
FIND_LIBRARY(LIBDL_LIBRARY NAMES ${LIBDL_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set LIBDL_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibDL DEFAULT_MSG LIBDL_LIBRARY LIBDL_INCLUDE_DIR)

IF(LIBDL_FOUND)
	SET( LIBDL_LIBRARIES ${LIBDL_LIBRARY} )
ELSE(LIBDL_FOUND)
	SET( LIBDL_LIBRARIES )
ENDIF(LIBDL_FOUND)

MARK_AS_ADVANCED( LIBDL_LIBRARY LIBDL_INCLUDE_DIR )
