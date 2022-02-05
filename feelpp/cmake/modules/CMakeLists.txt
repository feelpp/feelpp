###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2010-12-31
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


set( modules
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.version.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.directive.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.macros.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.macros.om.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.host.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.precompiled.headers.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.schedulers.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.sched.ccc.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.sched.loadleveler.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.sched.oar.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/feelpp.sched.slurm.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/ParseArguments.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cotire.cmake
  )

if ( FEELPP_HAS_ANN )
  set(modules ${modules} ${CMAKE_CURRENT_SOURCE_DIR}/FindANN.cmake )
endif()

if ( FEELPP_HAS_GLPK )
  set(modules ${modules} ${CMAKE_CURRENT_SOURCE_DIR}/FindGLPK.cmake )
endif()

if ( FEELPP_ENABLE_SYSTEM_EIGEN3 )
  set(modules ${modules} ${CMAKE_CURRENT_SOURCE_DIR}/FindEigen3.cmake )
endif()

if ( FEELPP_HAS_SCOTCH )
  set(modules ${modules} ${CMAKE_CURRENT_SOURCE_DIR}/FindScotch.cmake )
endif()

#FILE(GLOB modules "${CMAKE_CURRENT_SOURCE_DIR}/*.cmake")
list(APPEND modules ${CMAKE_CURRENT_SOURCE_DIR}/FindMetis.cmake )
INSTALL(FILES ${modules} DESTINATION share/feelpp/feel/cmake/modules COMPONENT Devel)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/Feel++Config.cmake DESTINATION share/feelpp/feel/cmake/modules COMPONENT Devel)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/FeelppContribConfig.cmake DESTINATION share/feelpp/feel/cmake/modules COMPONENT Devel)
#INSTALL(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/feelpp.install.config.cmake COMPONENT Devel)
