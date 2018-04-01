###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-02-03
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
feelpp_add_application(doc_residualestimator SRCS residualestimator.cpp LABELS doc)
if (OCTAVE_FOUND)
  add_octave_module(residualestimator)
endif()

if(FEELPP_HAS_GMSH_ADAPT_H)
  feelpp_add_application(doc_LShape SRCS Lshape.cpp  CFG Lshape.cfg Lshape3D.cfg LABELS doc )
else()
  message(STATUS "doc_LShape: compilation cancelled. Some headers (Openfile.h Field.h) needed for meshadaptation are missing.")
endif()

feelpp_add_application(meshadaptfromparamcurve SRCS meshadaptfromparamcurve.cpp CFG meshadaptfromls.cfg )
feelpp_add_application(meshadaptfromsimplels_d2 SRCS meshadaptfromsimplels.cpp DEFS FEELPP_DIM=2 CFG meshadaptfromls.cfg )
feelpp_add_application(meshadaptfromsimplels_d3 SRCS meshadaptfromsimplels.cpp DEFS FEELPP_DIM=3 CFG meshadaptfromls.cfg )
