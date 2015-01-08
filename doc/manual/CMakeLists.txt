# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2012 Université Joseph Fourier (Grenoble I)
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

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory(tutorial)

find_program(LATEX NAMES latex pdflatex)
find_program(CONVERT NAMES convert)
if ( LATEX AND CONVERT )
	option(FEELPP_ENABLE_MANUAL "enable compilation of Feel++ manual (requires latex)" ON)
else()
	option(FEELPP_ENABLE_MANUAL "enable compilation of Feel++ manual (requires latex)" OFF)
endif()

if ( FEELPP_ENABLE_MANUAL )
	include(UseLATEX OPTIONAL)

	foreach(cppfile "tutorial/*.cpp" "nonlinear/*.cpp" "../../quickstart/laplacian.cpp"
			"laplacian/laplacian.cpp"
			"stokes/stokes.cpp"
			)
		FILE(GLOB CPPFILE RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${cppfile})
		SET(CPPFILES "${CPPFILES};${CPPFILE}")
	endforeach()
	FILE(GLOB STYFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "styles/*.sty")

	#FILE(WRITE version.tex "\\newcommand{\\feelversion}{${FEELPP_VERSION_STRING}}")

	MESSAGE("Current revision is ${FEELPP_WC_REVISION}")
	MESSAGE("Last changed date is ${FEELPP_WC_LAST_CHANGED_DATE}")

	MACRO (TODAY RESULT)
		IF (WIN32)
			EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
			string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
				${RESULT} ${${RESULT}})
		ELSEIF(UNIX)
			EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
			string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
				${RESULT} ${${RESULT}})
		ELSE (WIN32)
			MESSAGE(SEND_ERROR "date not implemented")
			SET(${RESULT} 000000)
		ENDIF (WIN32)
	ENDMACRO (TODAY)
	TODAY(THEDATE)
	#FILE(APPEND version.tex "\\newcommand{\\feelbuildid}{${FEELPP_BUILDID}}")
	#FILE(APPEND version.tex "\\newcommand{\\feeldate}{${THEDATE}}")

	FILE(GLOB TUTO_DOC_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tutorial/*.tex)

	#CONFIGURE_FILE(tex/version.tex.in  tex/version.tex  @ONLY)

	#add_latex_document(tex/feelpp-manual.tex
	#  INPUTS
	#  ${TUTO_DOC_FILES}
	#  nonlinear/nonlinear.tex
	#  heatns/natural-convection-2d.tex
	#  dd/dd.tex
	#  maxwell/maxwell.tex
	#  tex/randomnotes.tex tex/gfdl.tex tex/howto.tex tex/version.tex ${CPPFILES} ${STYFILES}
	#  BIBFILES
	#  ../biblio/feelpp-manual.bib
	#  ../biblio/feelpp.bib
	#  ../biblio/feelpp-thesis.bib
	#  IMAGE_DIRS pngs pdfs pngs/mymesh pngs/heatsink pngs/Laplacian pngs/myadvection
	#  USE_INDEX
	#  DEFAULT_PDF
	#  MANGLE_TARGET_NAMES
	#  )
	#SET_TARGET_PROPERTIES(feelpp-manual_pdf PROPERTIES EXCLUDE_FROM_ALL true)
endif()

# this target allows to create all the required figures for the tutorial
add_custom_target(figures WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ )

macro( add_docex target )
	add_executable( feel_doc_${target} ${target}.cpp )
	target_link_libraries( feel_doc_${target} feelpp ${FEELPP_LIBRARIES} )
	get_target_property(target_executable feel_doc_${target} LOCATION)
	list(APPEND LIST_APPS ${target_executable})
	add_dependencies(figures feel_doc_${target})
	add_dependencies(doc feel_doc_${target})
	if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${target}.cfg)
		configure_file( ${target}.cfg ${target}.cfg)
	endif()
	add_custom_command(
		TARGET figures
		POST_BUILD
		COMMAND ${target_executable} ARGS --shape="simplex" --nochdir > ${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex.out 2>&1
		COMMAND ${target_executable} ARGS --shape="hypercube" --nochdir > ${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube.out 2>&1
		)
	foreach( dim 1 2 3 )
		#   message( STATUS "Generate scripts for graphical results for dimension ${dim}" )
		file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex-${dim}.pos "Include \"simplex-${dim}.geo\";\n Mesh ${dim}; Include \"tutorial.pos\";\n Exit;\n")
		file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube-${dim}.pos "Include \"hypercube-${dim}.geo\";\n Mesh ${dim}; Include \"tutorial.pos\";\n Exit;\n")
		add_custom_command(
			TARGET figures
			POST_BUILD
			COMMAND ${GMSH} ARGS -v 0 ${target}-simplex-${dim}.pos
			COMMAND ${GMSH} ARGS -v 0 ${target}-hypercube-${dim}.pos
			)
	endforeach()
	# test each backend
	foreach( BACKEND ${BACKEND_PETSC} ${BACKEND_TRILINOS})
		add_test( feel_doc_${target}-${BACKEND} feel_doc_${target} --hsize=0.2 --backend=${BACKEND} --no-export)
	endforeach()
	INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/feel_doc_${target} DESTINATION bin/ COMPONENT Bin)
endmacro()

add_custom_target(
	tut-prerequisites
	ALL
	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.pos
	${CMAKE_CURRENT_BINARY_DIR}/
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
add_dependencies(tut-prerequisites feelpp_tut_mymesh)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.pos  ${CMAKE_CURRENT_BINARY_DIR}/tutorial.pos)

if ( NOT FEELPP_MINIMAL_CONFIGURATION )

	#foreach( docex  myintegrals mymesh geotool myfunctionspace laplacian loadmesh residualestimator stokes stokes_good stokes_locking stokes_spurious bratu nonlinearpow mortar multi_levels threefields nirb nirbTest)
	# add_docex( ${docex} )
	#endforeach()

	#copy (${FEELPP_SOURCE_DIR}/data/medit/Cylref.mesh ${CMAKE_CURRENT_BINARY_DIR}/)


	# Natural convection problem
	#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
	#add_executable( feel_doc_convection convection.hpp convection_other.cpp convection_main.cpp convection_run.cpp convection_residual.cpp convection_jacobian.cpp convection_jacobian1.cpp convection_jacobian2.cpp convection_lo.cpp convection_lo2.cpp )
	#target_link_libraries( feel_doc_convection ${FEELPP_LIBRARIES} )



endif ( NOT FEELPP_MINIMAL_CONFIGURATION )

set(FEELPP_INSTALL_APPS ${FEELPP_INSTALL_APPS};${LIST_APPS} PARENT_SCOPE)

IF ( OCTAVE_FOUND )

	GET_DIRECTORY_PROPERTY(cmake_include_directories INCLUDE_DIRECTORIES)
	SET(octave_include_dirs)
	FOREACH(it ${cmake_include_directories})
		SET(octave_include_dirs ${octave_include_dirs} "-I${it}")
	ENDFOREACH(it)

	add_custom_target(octave WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ )

	macro(add_octave_module modulename)

		foreach( dim 2 )
			foreach( order 1 2 )
				FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${modulename}_${dim}_${order}_oct.cpp "#include<${modulename}_oct.hpp>\n\n RESIDUAL_ESTIMATOR( ${dim}, ${order} )\n" )

				add_library( ${modulename}_${dim}_${order}.oct MODULE  EXCLUDE_FROM_ALL ${modulename}_${dim}_${order}_oct.cpp  )
				add_dependencies(${modulename}_${dim}_${order}.oct
					${CMAKE_CURRENT_SOURCE_DIR}/${modulename}_oct.hpp
					${CMAKE_CURRENT_SOURCE_DIR}/${modulename}.hpp
					)
				add_dependencies(octave ${modulename}_${dim}_${order}.oct)
				target_link_libraries( ${modulename}_${dim}_${order}.oct feelpp ${FEELPP_LIBRARIES} ${Octave_LIBRARIES})
				set_target_properties( ${modulename}_${dim}_${order}.oct PROPERTIES PREFIX "" )
				set_target_properties( ${modulename}_${dim}_${order}.oct PROPERTIES SUFFIX "" )
			endforeach( order )
		endforeach( dim )
		configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${modulename}.m  ${CMAKE_CURRENT_BINARY_DIR}/residualestimator.m)
	endmacro(add_octave_module)




ENDIF( OCTAVE_FOUND )

add_subdirectory(python)

if (NOT FEELPP_MINIMAL_CONFIGURATION)
	foreach(DIR laplacian stokes advection accoustic dd heat solid ns heatns nonlinear adapt nirb fd maxwell opt eikonal eigen magnetostatic)
		add_subdirectory(${DIR})
	endforeach()
endif()
