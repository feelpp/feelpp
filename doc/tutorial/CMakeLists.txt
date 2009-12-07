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


add_latex_document(life-manual.tex
  INPUTS example-heatns.tex ${CPPFILES} ${STYFILES}
  BIBFILES life-manual.bib
  IMAGE_DIRS pngs pdfs
  USE_INDEX
  DEFAULT_PDF
  )

macro( add_docex target )
 add_executable( life_doc_${target} ${target}.cpp )
 target_link_libraries( life_doc_${target} ${LIFE_LIBRARIES} )
 # test each backend
 foreach( BACKEND gmm petsc trilinos)
   add_test( life_doc_${target}-${BENCH} life_doc_${target} --hsize=0.1 --backend=${BENCH})
 endforeach()
endmacro()

foreach( docex myapp mymesh myintegrals laplacian stokes bratu nonlinearpow )
  add_docex( ${docex} )
endforeach()

# Natural convection problem
add_executable( life_doc_convection convection.hpp convection_other.cpp convection_main.cpp convection_run.cpp convection_residual.cpp convection_jacobian.cpp convection_jacobian1.cpp convection_jacobian2.cpp convection_lo.cpp convection_lo2.cpp )
target_link_libraries( life_doc_convection ${LIFE_LIBRARIES} )
