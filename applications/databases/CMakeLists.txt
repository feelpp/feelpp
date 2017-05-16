###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 04 Apr 2017
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

if (FEELPP_HAS_HDF5 )

  if ( NOT TARGET app-databases )
    add_custom_target(app-databases)
  endif()

  if ( FEELPP_HAS_ACUSIM )
    feelpp_add_application( converter_acusim SRCS converter_database_acusim.cpp INCLUDE_IN_ALL PROJECT databases INSTALL EXEC DB_CONVERTER)
  endif()
  feelpp_add_application( export SRCS export_database.cpp INCLUDE_IN_ALL PROJECT databases INSTALL EXEC DB_EXPORT)
  feelpp_add_application( pod SRCS pod_database.cpp INCLUDE_IN_ALL PROJECT databases INSTALL EXEC DB_POD)

  add_dependencies(app-databases ${DB_CONVERTER} ${DB_EXPORT} ${DB_POD})
  add_custom_target(install-app-databases
    DEPENDS
    app-databases
    COMMAND
    "${CMAKE_COMMAND}"
    -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake")



endif()
