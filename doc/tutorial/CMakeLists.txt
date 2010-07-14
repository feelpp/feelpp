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
   COMMAND ${target_executable}
   ARGS >${CMAKE_CURRENT_BINARY_DIR}/${target}.out
   )

 # test each backend
 foreach( BACKEND petsc trilinos)
   add_test( life_doc_${target}-${BACKEND} life_doc_${target} --hsize=0.1 --backend=${BACKEND} --no-export)
 endforeach()
endmacro()

foreach( docex myapp ) #mymesh myintegrals laplacian  )
  add_docex( ${docex} )
endforeach()

if ( NOT LIFE_MINIMAL_CONFIGURATION )

foreach( docex stokes bratu nonlinearpow )
#  add_docex( ${docex} )
endforeach()

# Natural convection problem
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#add_executable( life_doc_convection convection.hpp convection_other.cpp convection_main.cpp convection_run.cpp convection_residual.cpp convection_jacobian.cpp convection_jacobian1.cpp convection_jacobian2.cpp convection_lo.cpp convection_lo2.cpp )
#target_link_libraries( life_doc_convection ${LIFE_LIBRARIES} )

endif ( NOT LIFE_MINIMAL_CONFIGURATION )
