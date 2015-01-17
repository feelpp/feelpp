###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-12-01
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Natural convection problem
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set(lib_convection
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_other.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_run.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_residual.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_jacobian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_jacobian1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_jacobian2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_lo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_lo2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_crb_other.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/convection_crb_residual.cpp )

if(FEELPP_DOC_ENABLE_CAVITY_CRB)
add_subdirectory( cavity )
endif()
