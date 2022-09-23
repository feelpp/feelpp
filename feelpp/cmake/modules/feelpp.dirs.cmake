###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 12 Jun 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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
if ( NOT FEELPP_INSTALL_DIR )
  set(FEELPP_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
  set(FEELPP_PREFIX ${FEELPP_INSTALL_DIR})
  set(FEELPP_DATADIR  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/feelpp/data)
  set(FEELPP_CASESDIR  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/feelpp/testcases)
  set(FEELPP_TOOLBOXCASESDIR  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/feelpp/testcases/toolboxes)
  set(FEELPP_LIBDIR  ${CMAKE_INSTALL_LIBDIR})
  set(FEELPP_PLUGINDIR ${CMAKE_INSTALL_LIBDIR}/feelpp/plugins)
endif()
