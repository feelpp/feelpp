###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 08 Mar 2019
#
#  Copyright (C) 2019 Feel++ Consortium
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

set( FEELPP_TOOLBOXES_CONFIG_FILE  ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/cmake/Feel++-ToolboxesConfig.cmake)
configure_file(${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/cmake/Feel++-ToolboxesConfig.cmake.in ${FEELPP_TOOLBOXES_CONFIG_FILE} @ONLY)
install(FILES ${FEELPP_TOOLBOXES_CONFIG_FILE} DESTINATION share/feelpp/toolboxes/cmake COMPONENT Devel)

install(EXPORT feelpp-toolboxes-targets
  NAMESPACE Feelpp::
  DESTINATION share/feelpp/toolboxes/cmake )
