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
  medit/Cylref.mesh
  medit/Box_with_a_sphere.mesh
  medit/Cylref.geo
  medit/Box_with_a_sphere.geo
  med/Cylref.med
  gmsh/lncmi_bitter/bitter_stl_be02.msh
  gmsh/lncmi_bitter/bitter_stl_be02_2.msh
  gmsh/lncmi_bitter/bitter_stl_be02-adapted.msh
  gmsh/lncmi_bitter/bitter_stl_be02.geo
  gmsh/lncmi_bitter/BE-02-01M.stp
  gmsh/bean.geo
  gmsh/cantilever.geo
  gmsh/operatorlagrangep1/simplex-3d-lagrange-order2.msh
  gmsh/operatorlagrangep1/hypercube-2d-lagrange-order2.msh
  gmsh/operatorlagrangep1/simplex-2d-lagrange-order2.msh
  gmsh/aorta/.gitignore
  gmsh/aorta/aorta_centerlines.geo
  gmsh/aorta/aorta2.stl
  gmsh/ReferenceTriangle.geo
  gmsh/pelvis/pelvis.stl
  gmsh/pelvis/pelvis.geo
  gmsh/magnetostatic/torus3d_quart_fer.geo
  gmsh/magnetostatic/torus3d.geo
  gmsh/magnetostatic/torus3d_quart.geo
  gmsh/magnetostatic/torusaxi_quart_fer.geo
  gmsh/magnetostatic/torusaxi.geo
  gmsh/fda/fda_Re500_fine.geo
  gmsh/fda/fda-2d.geo
  gmsh/fda/fda_Re500_coarse.geo
  gmsh/aneurism/aneurysm.geo
  gmsh/aneurism/aneurysm.stl
  gmsh/life.geo
  gmsh/feelpp/feel3d-ns.geo
  gmsh/feelpp/feel3d.geo
  gmsh/feelpp/feel.geo
  gmsh/feelpp/feel2d-ns.geo
  gmsh/primitives/torus.geo
  gmsh/primitives/rectangle.geo
  gmsh/primitives/cylinder.geo
  gmsh/one-elt-meshes-2d/one-elt-real-homo.msh
  gmsh/one-elt-meshes-2d/one-elt-ref.msh
  gmsh/one-elt-meshes-2d/one-elt-real-rot1.msh
  gmsh/one-elt-meshes-2d/one-elt-real-rot2.msh
  gmsh/one-elt-meshes-2d/two_elements_mesh.msh
  gmsh/steps/step-elmann-2009.geo
  gmsh/geom3D_conductor.geo
  gmsh/tripod.geo
  gmsh/marche.geo
  gmsh/one-elt-meshes-3d/one-elt-ref-3d.msh
  gmsh/one-elt-meshes-3d/one-elt-real-rotz.msh
  gmsh/one-elt-meshes-3d/one-elt-real-rotx.msh
  gmsh/one-elt-meshes-3d/one-elt-real-homo-3d.msh
  gmsh/one-elt-meshes-3d/one-elt-real-roty.msh )

install( FILES ${MESHES} DESTINATION share/feelpp/feel/geo )

if (EXISTS $ENV{HOME})
  # create $ENV{HOME}/feel/geo if it does not exist already: this is the geometry
  # files repository
  file(MAKE_DIRECTORY $ENV{HOME}/feel/geo)

  message(STATUS "Copying geo/msh files to $ENV{HOME}/feel/geo")

  foreach(exp "medit/*.geo" "medit/*.mesh" "gmsh/*.geo" "gmsh/*.stl" "gmsh/feelpp/feel.geo" )
    FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${exp}")
    foreach(file ${files})
      get_filename_component(filename ${file} NAME )
      message(STATUS ${filename})

      # Create a custom target so the geo files can be copied with make
      add_custom_target(copy_home_mesh_data_${filename})
      add_custom_command(TARGET copy_home_mesh_data_${filename}
        COMMAND ${CMAKE_COMMAND} -E copy ${file} $ENV{HOME}/feel/geo )
      add_dependencies(feelpp_mesh_data copy_home_mesh_data_${filename} )

      # also copies the files by default
      configure_file(${file} $ENV{HOME}/feel/geo/${filename} COPYONLY)
    endforeach()
  endforeach()
endif(EXISTS $ENV{HOME})

if (DEFINED ENV{FEELPP_WORKDIR})
  if (NOT EXISTS $ENV{FEELPP_WORKDIR})
    file(MAKE_DIRECTORY $ENV{FEELPP_WORKDIR}/geo)
  endif(NOT EXISTS $ENV{FEELPP_WORKDIR})

  message(STATUS "Copying geo/msh files to FEELPP_WORKDIR=$ENV{FEELPP_WORKDIR}/geo")

  foreach(exp "med/*.geo" "med/*.med" "medit/*.geo" "medit/*.mesh" "gmsh/*.geo" "gmsh/*.stl" "gmsh/feelpp/feel.geo" )
    FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${exp}")
    foreach(file ${files})
      get_filename_component(filename ${file} NAME )
      message(STATUS ${filename})

      # Create a custom target so the geo files can be copied with make
      add_custom_target(copy_workdir_mesh_data_${filename})
      add_custom_command(TARGET copy_workdir_mesh_data_${filename}
          COMMAND ${CMAKE_COMMAND} -E copy  ${file} $ENV{FEELPP_WORKDIR}/geo )
      add_dependencies(feelpp_mesh_data copy_workdir_mesh_data_${filename} )

      # also copies the files by default
      configure_file(${file}  $ENV{FEELPP_WORKDIR}/geo/${filename})
    endforeach()
  endforeach()
endif (DEFINED ENV{FEELPP_WORKDIR})
