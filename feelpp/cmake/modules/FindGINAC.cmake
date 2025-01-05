	###  TEMPLATE.txt.tpl; coding: utf-8 ---
	#
	# This Find module will set up the following:
	#
	#   GINAC_FOUND           : Boolean indicating whether GiNaC was found
	#   GINAC_INCLUDE_DIRS    : The include directories for GiNaC
	#   GINAC_LIBRARIES       : GiNaC library + dependencies
	#   GINAC_LIBRARY_DIRS    : Directories for GiNaC and CLN libraries
	#   GINAC_VERSION         : Determined version of GiNaC
	#
	#   and provides a modern "ginac::ginac" IMPORTED INTERFACE library target
	#   that can be linked against with "target_link_libraries(... ginac::ginac)".

	include(CheckLibraryExists)
	check_library_exists(dl dlopen "" HAVE_LIBDL)
	if(HAVE_LIBDL)
		if(NOT CMAKE_DL_LIBS)
			set(CMAKE_DL_LIBS "dl")
		endif()
		set(_ginac_dl_lib "${CMAKE_DL_LIBS}")
	endif()

	if (GINAC_INCLUDE_DIRS AND GINAC_LIBRARIES)
		set(GINAC_FIND_QUIETLY TRUE)
	endif()

	function(_ginac_headers_version _out_major _out_minor _out_patch _version_h)
		file(STRINGS "${_version_h}" _ginac_vinfo REGEX "^#define[ \t]+GINACLIB_.*_VERSION.*")
		if (NOT _ginac_vinfo)
			message(FATAL_ERROR "include file ${_version_h} does not exist or lacks GINACLIB_*_VERSION macros")
		endif()

		string(REGEX REPLACE "^.*GINACLIB_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1" ${_out_major} "${_ginac_vinfo}")
		string(REGEX REPLACE "^.*GINACLIB_MINOR_VERSION[ \t]+([0-9]+).*" "\\1" ${_out_minor} "${_ginac_vinfo}")
		string(REGEX REPLACE "^.*GINACLIB_MICRO_VERSION[ \t]+([0-9]+).*" "\\1" ${_out_patch} "${_ginac_vinfo}")

		if (NOT ${_out_major} MATCHES "[0-9]+")
			message(FATAL_ERROR "failed to determine GINACLIB_MAJOR_VERSION, got '${${_out_major}}'")
		endif()
		if (NOT ${_out_minor} MATCHES "[0-9]+")
			message(FATAL_ERROR "failed to determine GINACLIB_MINOR_VERSION, got '${${_out_minor}}'")
		endif()
		if (NOT ${_out_patch} MATCHES "[0-9]+")
			message(FATAL_ERROR "failed to determine GINACLIB_MICRO_VERSION, got '${${_out_patch}}'")
		endif()

		set(${_out_major} ${${_out_major}} PARENT_SCOPE)
		set(${_out_minor} ${${_out_minor}} PARENT_SCOPE)
		set(${_out_patch} ${${_out_patch}} PARENT_SCOPE)
	endfunction()

	set(GINAC_FOUND FALSE)
	set(GINAC_INCLUDE_DIRS)
	set(GINAC_LIBRARIES)
	set(GINAC_LIBRARY_DIRS)
	set(GINAC_VERSION)

	include(FindPkgConfig)
	find_package(CLN 1.3.2)  # optional, if you want to ensure CLN is found

	if (PKG_CONFIG_FOUND)
		pkg_check_modules(_ginac ginac)
	else()
		# Fallback if pkg-config not found
		set(_ginac_LIBRARIES ginac cln gmp)
	endif()

	if (NOT CLN_FOUND)
		set(GINAC_INCLUDE_DIRS "GINAC-NOTFOUND")
		set(GINAC_LIBRARIES "GINAC-NOTFOUND")
	else()
		# Locate include dir
		find_path(_ginac_include_dir
			NAMES ginac/ginac.h
			HINTS ${_ginac_INCLUDE_DIRS} $ENV{GINAC_DIR}/include
		)
		if (_ginac_include_dir)
			set(GINAC_INCLUDE_DIRS
				${_ginac_include_dir}
				${_ginac_INCLUDE_DIRS}
				${CLN_INCLUDE_DIR}
			)
			list(REMOVE_DUPLICATES GINAC_INCLUDE_DIRS)
		else()
			set(GINAC_INCLUDE_DIRS "GINAC-NOTFOUND")
			set(GINAC_LIBRARIES "GINAC-NOTFOUND")
			if (NOT GINAC_FIND_QUIETLY)
				message(FATAL_ERROR "couldn't find ginac.h")
			endif()
		endif()

		if (GINAC_INCLUDE_DIRS AND NOT GINAC_INCLUDE_DIRS STREQUAL "GINAC-NOTFOUND")
			# Locate library
			find_library(_ginac_lib
				NAMES libginac ginac
				HINTS ${_ginac_LIBRARY_DIRS} $ENV{GINAC_DIR}/lib
			)
			if (_ginac_lib)
				set(GINAC_LIBRARIES ${_ginac_lib} ${CLN_LIBRARIES})
				list(REMOVE_DUPLICATES GINAC_LIBRARIES)
			else()
				set(GINAC_LIBRARIES "GINAC-NOTFOUND")
				set(GINAC_INCLUDE_DIRS "GINAC-NOTFOUND")
				if (NOT GINAC_FIND_QUIETLY)
					message(FATAL_ERROR "couldn't find libginac")
				endif()
			endif()
		endif()
	endif()

	# If we reached here with GINAC_INCLUDE_DIRS != GINAC-NOTFOUND, do version checks
	if (GINAC_INCLUDE_DIRS AND NOT GINAC_INCLUDE_DIRS STREQUAL "GINAC-NOTFOUND")
		_ginac_headers_version(
			GINACLIB_MAJOR_VERSION
			GINACLIB_MINOR_VERSION
			GINACLIB_MICRO_VERSION
			"${_ginac_include_dir}/ginac/version.h"
		)
		set(GINAC_VERSION
			"${GINACLIB_MAJOR_VERSION}.${GINACLIB_MINOR_VERSION}.${GINACLIB_MICRO_VERSION}"
		)

		# Cross-check with pkg-config if available
		if (PKG_CONFIG_FOUND AND _ginac_VERSION AND NOT GINAC_VERSION VERSION_EQUAL _ginac_VERSION)
			if (NOT CLN_FIND_QUIETLY)
				message(FATAL_ERROR "pkg-config and version.h disagree: "
					"${_ginac_VERSION} vs. ${GINAC_VERSION}. Check your installation.")
			endif()
			set(GINAC_LIBRARIES "GINAC-NOTFOUND")
			set(GINAC_INCLUDE_DIRS "GINAC-NOTFOUND")
			set(GINAC_LIBRARY_DIRS)
			set(GINAC_VERSION)
		endif()
	endif()

	# If the library & includes are found, check that the library version matches the headers
	if (GINAC_INCLUDE_DIRS AND NOT GINAC_INCLUDE_DIRS STREQUAL "GINAC-NOTFOUND"
		AND GINAC_LIBRARIES AND NOT GINAC_LIBRARIES STREQUAL "GINAC-NOTFOUND"
		AND NOT CMAKE_CROSSCOMPILING
	)
		include(CheckCXXSourceRuns)
		set(_save_required_includes "${CMAKE_REQUIRED_INCLUDES}")
		set(_save_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
		set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${GINAC_INCLUDE_DIRS})
		set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${GINAC_LIBRARIES})

		check_cxx_source_runs("
			#include <ginac/version.h>
			#include <cln/version.h>
			#include <stdio.h>
			int main() {
				return (CL_VERSION_MAJOR == cln::version_major) &&
					(CL_VERSION_MINOR == cln::version_minor) &&
					(CL_VERSION_PATCHLEVEL == cln::version_patchlevel) &&
					(GINACLIB_MAJOR_VERSION == GiNaC::version_major) &&
					(GINACLIB_MINOR_VERSION == GiNaC::version_minor) &&
					(GINACLIB_MICRO_VERSION == GiNaC::version_micro) ? 0 : 1;
			}
		" _ginac_version_matches)

		set(CMAKE_REQUIRED_LIBRARIES "${_save_required_libraries}")
		set(CMAKE_REQUIRED_INCLUDES "${_save_required_includes}")

		if(NOT _ginac_version_matches)
			if (NOT GINAC_FIND_QUIETLY)
				message(FATAL_ERROR "header version differs from the library one, "
									"please check your GiNaC installation.")
			endif()
			set(GINAC_INCLUDE_DIRS "GINAC-NOTFOUND")
			set(GINAC_LIBRARIES "GINAC_NOTFOUND")
			set(GINAC_LIBRARY_DIRS)
			set(GINAC_VERSION)
		endif()
	endif()

	# If we have valid GINAC_LIBRARIES and GINAC_INCLUDE_DIRS, gather library dirs
	if(GINAC_LIBRARIES AND NOT GINAC_LIBRARIES STREQUAL "GINAC-NOTFOUND"
	AND GINAC_INCLUDE_DIRS AND NOT GINAC_INCLUDE_DIRS STREQUAL "GINAC-NOTFOUND"
	)
		set(GINAC_FOUND TRUE)
		set(_ginac_library_dirs "")
		foreach(_l ${GINAC_LIBRARIES})
			get_filename_component(_d "${_l}" DIRECTORY)
			list(APPEND _ginac_library_dirs "${_d}")
		endforeach()
		list(REMOVE_DUPLICATES _ginac_library_dirs)
		set(GINAC_LIBRARY_DIRS ${_ginac_library_dirs})
	endif()

	include(FindPackageHandleStandardArgs)
	FIND_PACKAGE_HANDLE_STANDARD_ARGS(
		GiNaC
		REQUIRED_VARS GINAC_LIBRARIES GINAC_INCLUDE_DIRS
		VERSION_VAR GINAC_VERSION
	)

	# ------------------------------------------------------------------------------
	#         Define a modern imported target: ginac::ginac
	# ------------------------------------------------------------------------------
	if(GINAC_FOUND)
		# Create an imported INTERFACE library named ginac::ginac
		if(NOT TARGET ginac::ginac)
			add_library(ginac::ginac INTERFACE IMPORTED)

			# Provide the usage requirements via INTERFACE properties:
			set_target_properties(ginac::ginac PROPERTIES
				INTERFACE_INCLUDE_DIRECTORIES "${GINAC_INCLUDE_DIRS}"
				INTERFACE_LINK_LIBRARIES "${GINAC_LIBRARIES}"
				# Optional: Add definitions if you want them exposed
				# e.g. INTERFACE_COMPILE_DEFINITIONS "GINAC_USING_SOMETHING"
			)
			# If we need dl for JIT:
			if(HAVE_LIBDL)
				set_property(TARGET ginac::ginac APPEND PROPERTY
					INTERFACE_LINK_LIBRARIES "${_ginac_dl_lib}"
				)
				set_property(TARGET ginac::ginac APPEND PROPERTY
    				INTERFACE_COMPILE_DEFINITIONS "HAVE_LIBDL"
				)
			endif()
		endif()
	endif()