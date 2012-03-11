###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-04-09
#
#  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
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
GET_DIRECTORY_PROPERTY(cmake_include_directories INCLUDE_DIRECTORIES)
SET(octave_include_dirs)
FOREACH(it ${cmake_include_directories})
  SET(octave_include_dirs ${octave_include_dirs} "-I${it}")
ENDFOREACH(it)

crb_add_octave_module( opuseadsmfem opuseadsmfem.cpp
  LINK_LIBRARIES feel++_crb_models ${Octave_LIBRARIES}
  SCRIPTS opuseadsmfemtest.m )
crb_add_octave_module( opuseadspfem opuseadspfem.cpp
  LINK_LIBRARIES feel++_crb_models ${Octave_LIBRARIES}
  SCRIPTS opuseadspfemtest.m )

if (FEELPP_ENABLE_OCTAVE)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/opuseadspfemtest.m ${CMAKE_CURRENT_BINARY_DIR}/opuseadsmfemtest.m DESTINATION ${FEELPP_M_DIR})
endif()

# add_custom_target(  opuseadsmfemoct DEPENDS  opuseadsmfem.oct )
# add_custom_command(
#   OUTPUT opuseadsmfem.oct
#   DEPENDS opuseadsmfem.cpp
#   COMMAND "mkoctfile"
#   ARGS
#   ${octave_include_dirs}
#   ${CMAKE_CURRENT_SOURCE_DIR}/opuseadsmfem.cpp
#   ${FEELPP_LIBRARIES} )
