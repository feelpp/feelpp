# ##  CMakeLists.txt; coding: utf-8 ---

# Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# Date: 12 Jun 2017
#
# Copyright (C) 2017 Feel++ Consortium
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
include(GNUInstallDirs)

if(NOT FEELPP_INSTALL_DIR)
  set(FEELPP_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
  set(FEELPP_PREFIX ${FEELPP_INSTALL_DIR})
endif()

set(FEELPP_DATADIR ${FEELPP_INSTALL_DIR}/${CMAKE_INSTALL_DATADIR}/feelpp/data)
set(FEELPP_CASESDIR ${FEELPP_INSTALL_DIR}/${CMAKE_INSTALL_DATADIR}/feelpp/testcases)
set(FEELPP_TOOLBOXCASESDIR ${FEELPP_INSTALL_DIR}/${CMAKE_INSTALL_DATADIR}/feelpp/testcases/toolboxes)
set(FEELPP_LIBDIR ${CMAKE_INSTALL_LIBDIR})
set(FEELPP_PLUGINDIR ${CMAKE_INSTALL_LIBDIR}/feelpp/plugins)
message(STATUS "[feelpp] FEELPP_INSTALL_DIR: ${FEELPP_INSTALL_DIR}")
message(STATUS "[feelpp]         FEELPP_DIR: ${FEELPP_DIR}")
message(STATUS "[feelpp]      FEELPP_PREFIX: ${FEELPP_PREFIX}")
message(STATUS "[feelpp]     FEELPP_DATADIR: ${FEELPP_DATADIR}")
message(STATUS "[feelpp]      FEELPP_LIBDIR: ${FEELPP_LIBDIR}")
message(STATUS "[feelpp]   FEELPP_PLUGINDIR: ${FEELPP_PLUGINDIR}")