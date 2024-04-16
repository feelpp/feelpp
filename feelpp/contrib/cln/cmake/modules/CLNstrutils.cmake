# various string manipulation routines

if (NOT _cln_strutils_included)
	set(_cln_strutils_included 1)

	macro(cl_string_join var strlist delimiter)
	set(_ret "")
	foreach(_str ${strlist})
		if ("${_ret}" STREQUAL "")
			set(_ret "${_str}")
		else()
			set(_ret "${_ret}${delimiter}${_str}")
		endif()
	endforeach()
	set(${var} "${_ret}")
	endmacro()

endif(NOT _cln_strutils_included)

