###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
      configure_file(${file} $ENV{HOME}/feel/geo/${filename})
    endforeach()
  endforeach()
endif(EXISTS $ENV{HOME})

if (DEFINED ENV{FEELPP_WORKDIR})
  if (NOT EXISTS $ENV{FEELPP_WORKDIR})
    file(MAKE_DIRECTORY $ENV{FEELPP_WORKDIR}/geo)
  endif(NOT EXISTS $ENV{FEELPP_WORKDIR})

  message(STATUS "Copying geo/msh files to FEELPP_WORKDIR=$ENV{FEELPP_WORKDIR}/geo")

  foreach(exp "medit/*.geo" "medit/*.mesh" "gmsh/*.geo" "gmsh/*.stl" "gmsh/feelpp/feel.geo" )
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
