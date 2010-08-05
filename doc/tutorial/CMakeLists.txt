###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
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

FILE(WRITE version.tex "\\newcommand{\\lifeversion}{${LIFE_VERSION_STRING}}")

# FIND_PACKAGE(Subversion)
# IF(Subversion_FOUND)
#   Subversion_WC_INFO(${PROJECT_SOURCE_DIR} LIFE)
#   MESSAGE("Current revision is ${LIFE_WC_REVISION}")
#   MESSAGE("Last changed date is ${LIFE_WC_LAST_CHANGED_DATE}")

#   FILE(APPEND version.tex "\\newcommand{\\lifesvnrevision}{${LIFE_WC_REVISION}}")
#   FILE(APPEND version.tex "\\newcommand{\\lifesvndate}{${LIFE_WC_LAST_CHANGED_DATE}}")
# ENDIF(Subversion_FOUND)


add_latex_document(life-manual.tex
  INPUTS feel.tex randomnotes.tex gfdl.tex example-heatns.tex version.tex ${CPPFILES} ${STYFILES}
  BIBFILES life-manual.bib
  IMAGE_DIRS pngs pdfs
  USE_INDEX
  DEFAULT_PDF
  )

macro( add_docex target )
 add_executable( life_doc_${target} ${target}.cpp )
 target_link_libraries( life_doc_${target} ${LIFE_LIBRARIES} )
 get_target_property(target_executable life_doc_${target} LOCATION)
 add_custom_command(
   TARGET life_doc_${target}
   POST_BUILD
   COMMAND ${target_executable} ARGS --shape="simplex" --nochdir > ${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex.out 2>&1
   COMMAND ${target_executable} ARGS --shape="hypercube" --nochdir > ${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube.out 2>&1
   )
 foreach( dim 1 2 3 )
   message( STATUS "Generate scripts for graphical results for dimension ${dim}" )
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex-${dim}.pos "Include \"simplex-${dim}.geo\";\n Mesh ${dim}; Include \"tutorial.pos\";\n Exit;\n")
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube-${dim}.pos "Include \"hypercube-${dim}.geo\";\n Mesh ${dim}; Include \"tutorial.pos\";\n Exit;\n")
   add_custom_command(
     TARGET life_doc_${target}
     POST_BUILD
     COMMAND ${GMSH} ARGS -v 0 ${target}-simplex-${dim}.pos
     COMMAND ${GMSH} ARGS -v 0 ${target}-hypercube-${dim}.pos
   )
 endforeach()
 # test each backend
 foreach( BACKEND petsc trilinos)
   add_test( life_doc_${target}-${BACKEND} life_doc_${target} --hsize=0.1 --backend=${BACKEND} --no-export)
 endforeach()
endmacro()

foreach( docex myintegrals mymesh )
  add_docex( ${docex} )
endforeach()

macro( add_docex2 target )
 add_executable( life_doc_${target} ${target}.cpp )
 target_link_libraries( life_doc_${target} ${LIFE_LIBRARIES} )
 get_target_property(target_executable life_doc_${target} LOCATION)
 # add_custom_command(
 #   TARGET life_doc_${target}
 #   POST_BUILD
 #   COMMAND ${target_executable} ARGS --hsize=0.05 --shape="simplex" --nochdir --exporter-format=gmsh > ${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex.out 2>&1
 #   COMMAND ${target_executable} ARGS --hsize=0.075 --shape="hypercube" --nochdir --exporter-format=gmsh > ${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube.out 2>&1
 #   )
 # foreach( dim 1 2 3 )
 #   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex-${dim}.pos "Merge \"${CMAKE_CURRENT_BINARY_DIR}/${target}-simplex-${dim}-1_0.msh\";\nMerge \"tutorial.pos\";\nExit;\n")
 #   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube-${dim}.pos "Merge \"${CMAKE_CURRENT_BINARY_DIR}/${target}-hypercube-${dim}-1_0.msh\";\nMerge \"tutorial.pos\";\nExit;\n")
 #   add_custom_command(
 #     TARGET life_doc_${target}
 #     POST_BUILD
 #     COMMAND ${GMSH} ARGS  ${target}-simplex-${dim}.pos
 #     COMMAND ${GMSH} ARGS  ${target}-hypercube-${dim}.pos
 #   )
 # endforeach()
 # test each backend
 foreach( BACKEND petsc trilinos)
   add_test( life_doc_${target}-${BACKEND} life_doc_${target} --hsize=0.1 --backend=${BACKEND} --no-export)
 endforeach()
endmacro()
foreach( docex myfunctionspace laplacian residualestimator)
  add_docex2( ${docex} )
endforeach()

add_custom_target(
  tut-prerequisites
  ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.pos
                                   ${CMAKE_CURRENT_BINARY_DIR}/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(tut-prerequisites life_doc_mymesh)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.pos  ${CMAKE_CURRENT_BINARY_DIR}/tutorial.pos)

if ( NOT LIFE_MINIMAL_CONFIGURATION )

foreach( docex stokes bratu nonlinearpow )
#  add_docex( ${docex} )
endforeach()

# Natural convection problem
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#add_executable( life_doc_convection convection.hpp convection_other.cpp convection_main.cpp convection_run.cpp convection_residual.cpp convection_jacobian.cpp convection_jacobian1.cpp convection_jacobian2.cpp convection_lo.cpp convection_lo2.cpp )
#target_link_libraries( life_doc_convection ${LIFE_LIBRARIES} )

endif ( NOT LIFE_MINIMAL_CONFIGURATION )

IF ( OCTAVE_FOUND )

GET_DIRECTORY_PROPERTY(cmake_include_directories INCLUDE_DIRECTORIES)
SET(octave_include_dirs)
FOREACH(it ${cmake_include_directories})
  SET(octave_include_dirs ${octave_include_dirs} "-I${it}")
ENDFOREACH(it)

macro(add_octave_module modulename)

  foreach( dim 1 )
    foreach( order 1 )
      FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${modulename}_${dim}_${order}_oct.cpp "#include<${modulename}_oct.hpp>\n\n RESIDUAL_ESTIMATOR( ${dim}, ${order} )\n" )

      add_library( ${modulename}_${dim}_${order}.oct MODULE  ${modulename}_${dim}_${order}_oct.cpp  )
      add_dependencies(${modulename}_${dim}_${order}.oct
        ${CMAKE_CURRENT_SOURCE_DIR}/${modulename}_oct.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/${modulename}.hpp
        )
      target_link_libraries( ${modulename}_${dim}_${order}.oct ${LIFE_LIBRARIES})
      set_target_properties( ${modulename}_${dim}_${order}.oct PROPERTIES PREFIX "" )
      set_target_properties( ${modulename}_${dim}_${order}.oct PROPERTIES SUFFIX "" )
    endforeach( order )
  endforeach( dim )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${modulename}.m  ${CMAKE_CURRENT_BINARY_DIR}/residualestimator.m)
endmacro(add_octave_module)


add_octave_module( residualestimator )

ENDIF( OCTAVE_FOUND )
