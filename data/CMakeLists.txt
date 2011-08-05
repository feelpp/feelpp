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
add_custom_target(feel_mesh_data)

if (EXISTS $ENV{HOME})
  # create $ENV{HOME}/feel/geo if it does not exist already: this is the geometry
  # files repository
  file(MAKE_DIRECTORY $ENV{HOME}/feel/geo)

  FILE(GLOB meditfiles "${CMAKE_CURRENT_SOURCE_DIR}/medit/*.[geomesh]")
  FILE(GLOB gmshfiles "${CMAKE_CURRENT_SOURCE_DIR}/gmsh/*.[geomeshstl]")
  message(STATUS "copying mesh files: ${meditfiles} ${gmshfiles}")
  foreach(file ${meditfiles} ${gmshfiles})
    get_filename_component(filename ${file} NAME )
    add_custom_target(copy_mesh_data_${filename}
      ALL
      COMMAND ${CMAKE_COMMAND} -E copy  ${file} $ENV{HOME}/feel/geo )
    add_dependencies(feel_mesh_data copy_mesh_data_${filename} )
  endforeach()
endif(EXISTS $ENV{HOME})
