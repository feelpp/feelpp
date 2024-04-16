# Try to find the GMP librairies
#  GMP_FOUND - system has GMP lib
#  GMP_INCLUDE_DIR - the GMP include directory
#  GMP_LIBRARIES - Libraries needed to use GMP

# Copyright (c) 2006, Laurent Montel, <montel@kde.org>
# Adjusted by Alexei Sheplyakov.
# Redistribution and use is allowed according to the terms of the BSD license.


if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
  # Already in cache, be silent
  set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDE_DIR AND GMP_LIBRARIES)

find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp libgmp)


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARIES)

if (GMP_FOUND AND NOT TARGET gmp::gmp)
	set(_found_shared_libgmp FALSE)
	get_filename_component(_gmplib_suffix ${GMP_LIBRARIES} EXT)
	if (_gmplib_suffix STREQUAL ${CMAKE_STATIC_LIBRARY_SUFFIX})
		# XXX: msvc has the same suffix both for the static and import libs
		add_library(gmp::gmp STATIC IMPORTED)
	else()
		set(_found_shared_libgmp TRUE)
		add_library(gmp::gmp SHARED IMPORTED)
	endif()
	set_target_properties(gmp::gmp PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
	)
	if (WIN32 AND _found_shared_libgmp)
		set_target_properties(gmp::gmp PROPERTIES
			IMPORTED_IMPLIB ${GMP_LIBRARIES}
		)
	else()
		set_target_properties(gmp::gmp PROPERTIES
			IMPORTED_LOCATION ${GMP_LIBRARIES}
		)
	endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)

