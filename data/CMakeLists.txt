###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2011-07-14
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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
add_custom_target(feelpp_mesh_data ALL)

set( MESHES
geo/Box_with_a_sphere.geo
geo/Box_with_a_sphere.mesh
geo/Cylref.geo
geo/Cylref.med
geo/Cylref.mesh
geo/ReferenceTriangle.geo
geo/aneurism/aneurysm.geo
geo/aneurism/aneurysm.stl
geo/aorta/aorta2.stl
geo/aorta/aorta_centerlines.geo
geo/bean.geo
geo/cantilever.geo
geo/fda/fda-2d.geo
geo/fda/fda_Re500_coarse.geo
geo/fda/fda_Re500_fine.geo
geo/feelpp/feel.geo
geo/feelpp/feel2d-ns.geo
geo/feelpp/feel3d-ns.geo
geo/feelpp/feel3d.geo
geo/geom3D_conductor.geo
geo/life.geo
geo/lncmi_bitter/BE-02-01M.stp
geo/lncmi_bitter/bitter_stl_be02-adapted.msh
geo/lncmi_bitter/bitter_stl_be02.geo
geo/lncmi_bitter/bitter_stl_be02.msh
geo/lncmi_bitter/bitter_stl_be02_2.msh
geo/magnetostatic/torus3d.geo
geo/magnetostatic/torus3d_quart.geo
geo/magnetostatic/torus3d_quart_fer.geo
geo/magnetostatic/torusaxi.geo
geo/magnetostatic/torusaxi_quart_fer.geo
geo/marche.geo
geo/pelvis/pelvis.geo
geo/pelvis/pelvis.stl
geo/primitives/cylinder.geo
geo/primitives/rectangle.geo
geo/primitives/torus.geo
geo/steps/step-elmann-2009.geo
geo/tripod.geo
geo/tripod.med
geo/tripod.mesh
geo/square_squareHole.geo
geo/square_twoMaterials.geo
geo/cube_twoMaterials.geo
geo/cube_cubeHole.geo
)

install( FILES ${MESHES} DESTINATION share/feelpp/data/geo )

if (EXISTS $ENV{HOME})
  # create $ENV{HOME}/feel/geo if it does not exist already: this is the geometry
  # files repository
  file(MAKE_DIRECTORY $ENV{HOME}/feel/geo)
  message(STATUS "Copying geo/msh files to $ENV{HOME}/geo")
  install( FILES ${MESHES} DESTINATION $ENV{HOME}/feel/geo) 
endif()
if (DEFINED ENV{FEELPP_WORKDIR})
  if (NOT EXISTS $ENV{FEELPP_WORKDIR})
    file(MAKE_DIRECTORY $ENV{FEELPP_WORKDIR}/geo)
  endif(NOT EXISTS $ENV{FEELPP_WORKDIR})
  
  message(STATUS "Copying geo/msh files to FEELPP_WORKDIR=$ENV{FEELPP_WORKDIR}/geo")
  install( FILES ${MESHES} DESTINATION $ENV{FEELPP_WORKDIR}/geo) 
endif (DEFINED ENV{FEELPP_WORKDIR})

add_subdirectory(operatorlagrangep1)

set( ONE_ELT_2D_MESHES
  geo/one-elt-meshes-2d/one-elt-real-homo.msh
  geo/one-elt-meshes-2d/one-elt-real-rot1.msh
  geo/one-elt-meshes-2d/one-elt-real-rot2.msh
  geo/one-elt-meshes-2d/one-elt-ref.msh
  geo/one-elt-meshes-2d/two_elements_mesh.msh
  )
set( ONE_ELT_3D_MESHES
  geo/one-elt-meshes-3d/one-elt-real-homo-3d.msh
  geo/one-elt-meshes-3d/one-elt-real-rotx.msh
  geo/one-elt-meshes-3d/one-elt-real-roty.msh
  geo/one-elt-meshes-3d/one-elt-real-rotz.msh
  geo/one-elt-meshes-3d/one-elt-ref-3d.msh
  )
install( FILES ${ONE_ELT_2D_MESHES} DESTINATION share/feelpp/data/geo/one-elt-meshes-2d/ )
install( FILES ${ONE_ELT_3D_MESHES} DESTINATION share/feelpp/data/geo/one-elt-meshes-3d/ )
