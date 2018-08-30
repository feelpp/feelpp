# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
IF ( ENABLE_INSTANTIATION_MODE )
  SET( FEELFILTERS_INST_SRCS

    # exporterensight.cpp exportergmsh.cpp # exportergnuplot.cpp
    )
endif()

set(FEELFILTERS_GMSH_SRCS "")
if(FEELPP_HAS_GMSH_H AND FEELPP_HAS_GMSH_LIBRARY)
  set(FEELFILTERS_GMSH_SRCS gmshhypercubedomain.cpp gmshsimplexdomain.cpp gmshellipsoiddomain.cpp gmsh.cpp geotool.cpp)
endif()

SET(SRC
  importergmsh.cpp gmshenums.cpp unitsegment.cpp unitsquare.cpp unitcube.cpp
  #exporterfactorysimplex.cpp exporterfactoryhypercube.cpp
  ${FEELFILTERS_GMSH_SRCS}
  ${FEELFILTERS_INST_SRCS}
  detail/fileindex.cpp
  loadcsv.cpp
  hbf.cpp
  straightenmesh_inst_1d.cpp straightenmesh_inst_1d_p2.cpp straightenmesh_inst_2d.cpp straightenmesh_inst_2d_ho.cpp straightenmesh_inst_3d.cpp straightenmesh_inst_3d_ho.cpp
  )

if(FEELPP_HAS_PARAVIEW)
    SET(SRC ${SRC} vtkBaseInsituPipeline.cpp)
endif()

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
feel_append_src(feelfilters "${SRC};${HDR}")
