# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
include(UseLATEX)

FILE(GLOB CPPFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
FILE(GLOB STYFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "styles/*.sty")

FILE(WRITE version.tex "\\newcommand{\\feelversion}{${FEEL_VERSION_STRING}}")

FIND_PACKAGE(Subversion)
IF(Subversion_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn" )
   Subversion_WC_INFO(${PROJECT_SOURCE_DIR} FEEL)
   MESSAGE("Current revision is ${FEEL_WC_REVISION}")
   MESSAGE("Last changed date is ${FEEL_WC_LAST_CHANGED_DATE}")

   FILE(APPEND version.tex "\\newcommand{\\feelsvnrevision}{${FEEL_WC_REVISION}}")
   FILE(APPEND version.tex "\\newcommand{\\feelsvndate}{${FEEL_WC_LAST_CHANGED_DATE}}")
ENDIF()


add_latex_document(feel-manual.tex
  INPUTS
  tutorial/building.tex
  tutorial/getting-started.tex
  tutorial/feel.tex
  nonlinear/nonlinear.tex
  heatns/natural-convection-2d.tex
  heat/heatsink.tex
  dd/dd.tex
  maxwell/maxwell.tex
  randomnotes.tex gfdl.tex howto.tex version.tex ${CPPFILES} ${STYFILES}
  BIBFILES
  ../biblio/feel-manual.bib
  ../biblio/feelpp.bib
  ../biblio/feelpp-thesis.bib
  IMAGE_DIRS pngs pdfs pngs/mymesh pngs/heatsink
  USE_INDEX
  DEFAULT_PDF
  )
SET_TARGET_PROPERTIES(pdf PROPERTIES EXCLUDE_FROM_ALL true)


# this target allows to create all the required figures for the tutorial
add_custom_target(figures WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ )

macro( add_docex target )
 add_executable( feel_doc_${target} ${target}.cpp )
 target_link_libraries( feel_doc_${target} ${FEEL_LIBRARIES} )
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
endmacro()

add_custom_target(
  tut-prerequisites
  ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.pos
                                   ${CMAKE_CURRENT_BINARY_DIR}/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(tut-prerequisites feel_doc_mymesh)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.pos  ${CMAKE_CURRENT_BINARY_DIR}/tutorial.pos)

if ( NOT FEEL_MINIMAL_CONFIGURATION )

# # foreach( docex  myintegrals mymesh geotool myfunctionspace laplacian loadmesh residualestimator stokes stokes_good stokes_locking stokes_spurious bratu nonlinearpow mortar multi_levels threefields nirb nirbTest)
# #   add_docex( ${docex} )
# # endforeach()

#copy (${FEEL_SOURCE_DIR}/data/medit/Cylref.mesh ${CMAKE_CURRENT_BINARY_DIR}/)


# Natural convection problem
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#add_executable( feel_doc_convection convection.hpp convection_other.cpp convection_main.cpp convection_run.cpp convection_residual.cpp convection_jacobian.cpp convection_jacobian1.cpp convection_jacobian2.cpp convection_lo.cpp convection_lo2.cpp )
#target_link_libraries( feel_doc_convection ${FEEL_LIBRARIES} )



endif ( NOT FEEL_MINIMAL_CONFIGURATION )

set(FEEL_INSTALL_APPS ${FEEL_INSTALL_APPS};${LIST_APPS} PARENT_SCOPE)
message(STATUS "doc tutorial: ${LIST_APPS}")


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
      target_link_libraries( ${modulename}_${dim}_${order}.oct ${FEEL_LIBRARIES} ${Octave_LIBRARIES})
      set_target_properties( ${modulename}_${dim}_${order}.oct PROPERTIES PREFIX "" )
      set_target_properties( ${modulename}_${dim}_${order}.oct PROPERTIES SUFFIX "" )
    endforeach( order )
  endforeach( dim )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${modulename}.m  ${CMAKE_CURRENT_BINARY_DIR}/residualestimator.m)
endmacro(add_octave_module)




ENDIF( OCTAVE_FOUND )

add_subdirectory(tutorial)
add_subdirectory(dd)

#add_subdirectory(heat)
#add_subdirectory(heatns)
add_subdirectory(nonlinear)
add_subdirectory(adapt)
add_subdirectory(nirb)
add_subdirectory(mcs)
