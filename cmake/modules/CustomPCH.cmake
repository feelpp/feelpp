
# Author(s): Alexandre Ancel <ancel@unistra.fr>
#       Date: 2013-09-18
#
#  Copyright (C) 2013 Université de Strasbourg
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

# Useful links:
# - flightcrew project: http://code.google.com/p/flightcrew/source/browse/cmake_extras/CustomPCH.cmake
# - cmake bug report for PCH: http://www.cmake.org/Bug/view.php?id=1260
# - qtcreator implementation of PCH: https://github.com/loaden/qtcreator/blob/wip/cmake/cmake/PrecompiledHeader.cmake
# - cotire: a project for generating PCH: https://github.com/sakra/cotire

macro(get_source_language sources lang)
	set(lang "C")
	#MESSAGE("${sources}")
	foreach(src ${sources})
		get_source_file_property(lang "${src}" LANGUAGE)
	endforeach()
endmacro()

macro(add_precompiled_header target sources headers)
	
	# get the language currently used
	set(lang "C")
	get_source_language(${sources} lang)

	set(pchExt pch)
	if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		set(pchExt gch)
	endif()

# Clang support does not work yet, but is present
if((CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))

		# rebuild the exact compilation command from the target for the header file
		# get defines
		get_directory_property(value DEFINITIONS)
		set(compiler_flags ${value})

		# get compilation flags
		string(TOUPPER "CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE}" name)
		set(compiler_flags ${${name}} ${compiler_flags})
		string(TOUPPER "CMAKE_${lang}_FLAGS" name)
		set(compiler_flags ${${name}} ${compiler_flags})
		string(TOUPPER "CMAKE_${lang}_COMPILE_ARG1" name)
		set(compiler_flags ${${name}} ${compiler_flags})

		# needed to unstringify the received values
		# and have a correct compilation command
		separate_arguments(compiler_flags)

		# Append language specification to flags
		if(lang STREQUAL "C")
			set(compiler_flags ${compiler_flags} -x c-header)
		elseif(lang STREQUAL "CXX")
			set(compiler_flags ${compiler_flags} -x c++-header)
		endif()

		# Append includes to commandline
		get_target_property(includes ${execname} INCLUDE_DIRECTORIES)
		#MESSAGE("INCLUDES ${value}")

		foreach(inc ${includes})
			set(compiler_flags ${compiler_flags} -I${inc})

			# Check if the headers to be precompiled are in this include directory
			foreach(head ${headers})
				if(EXISTS ${inc}/${head})
					set(header_map_${head} ${inc})
					# don't put a break here it seems to break the two loops
				endif()
			endforeach()
		endforeach()

		# for each header to precompile
		# precompile it in its own directory
		foreach(head ${headers})
			#MESSAGE(${head})
			set(headdir header_map_${head})
			if(DEFINED ${headdir})
				set(fullpath ${${headdir}}/${head})
				get_filename_component(head ${fullpath} NAME)
				get_filename_component(headpath ${fullpath} PATH)
				set(pchHeader ${head})
				set(pchBinary ${head}.${pchExt})
				set(compiler ${CMAKE_${lang}_COMPILER})

				# set generation part
				set(pchGenerator -c ${pchHeader} -o ${pchBinary})


				add_custom_command(OUTPUT "${headpath}/${pchBinary}"
		        COMMAND cd ${headpath} && ${compiler} ${compiler_flags} ${pchGenerator}
						DEPENDS "${headpath}/${pchHeader}"
						COMMENT "Building precompiled header ${headpath}/${pchBinary}"
                        IMPLICIT_DEPENDS CXX ${headpath}/${pchHeader}
				    )
		    add_custom_target( ${target}_${head}_pch
				        DEPENDS "${headpath}/${pchHeader}" "${headpath}/${pchBinary}"
								    )
		    add_dependencies(${target} ${target}_${head}_pch)

	    	set_property(SOURCE ${sources} APPEND_STRING PROPERTY
					    	COMPILE_FLAGS "-include \"${headpath}/${pchHeader}\" -Winvalid-pch "
						)

				set_property(SOURCE ${sources} APPEND PROPERTY
								OBJECT_DEPENDS ${headpath}/${pchHeader}
						)

			endif()
		endforeach()
	else()
		#message(WARNING "Precompiled header not yet implemented for this compiler")
	endif()
endmacro()
