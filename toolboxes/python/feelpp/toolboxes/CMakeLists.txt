###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 06 Aug 2018
#
#  Copyright (C) 2018 Feel++ Consortium
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
cmake_minimum_required(VERSION 3.21)

add_custom_target(pyfeelpptoolboxes)
if ( TARGET pyfeelpp )
  add_dependencies(pyfeelpptoolboxes pyfeelpp)
endif()

set(PYFILES __init__.py)
feelpp_stage_python_files(
  TARGET pyfeelpptoolboxes_package_pyfiles
  DESTINATION feelpp/toolboxes
  FILES ${PYFILES}
)
add_dependencies(pyfeelpptoolboxes pyfeelpptoolboxes_package_pyfiles)
install(FILES ${PYFILES}
        DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/toolboxes
        COMPONENT ${FEELPP_PYTHON_INSTALL_COMPONENT})

add_subdirectory( core )


if ( TARGET  Feelpp::feelpp_toolbox_coefficientformpdes_lib )
add_subdirectory( coefficientformpdes )
add_dependencies(pyfeelpptoolboxes _cfpdes)
endif()

if ( TARGET  Feelpp::feelpp_toolbox_heat_lib )
  add_subdirectory( heat )
  add_dependencies(pyfeelpptoolboxes _heat)
endif()

if ( TARGET  Feelpp::feelpp_toolbox_electric_lib )
  add_subdirectory( electric )
  add_dependencies(pyfeelpptoolboxes _electric)
endif()
if ( TARGET  Feelpp::feelpp_toolbox_thermoelectric_lib  )
  add_subdirectory( thermoelectric )
  add_dependencies(pyfeelpptoolboxes _thermoelectric)
endif()
if ( TARGET  Feelpp::feelpp_toolbox_solid_lib )
  add_subdirectory( solid )
  add_dependencies(pyfeelpptoolboxes _solid)
endif()
if ( TARGET  Feelpp::feelpp_toolbox_fluid_lib )
  add_subdirectory( fluid )
  add_dependencies(pyfeelpptoolboxes _fluid)
endif()
if ( TARGET  Feelpp::feelpp_toolbox_heatfluid_lib )
  add_subdirectory( heatfluid )
  add_dependencies( pyfeelpptoolboxes _heatfluid )
endif()
if ( TARGET  Feelpp::feelpp_toolbox_hdg_lib )
  add_subdirectory( hdg )
  add_dependencies(pyfeelpptoolboxes _hdg)
endif()

# Advection, Maxwell, and FSI wrappers are kept out of the main build until
# they have maintained tests and dedicated CI coverage. The stale LevelSet and
# MultiFluid thermoelectric placeholder bindings were removed.
# if ( TARGET  Feelpp::feelpp_toolbox_fsi_lib )
#   add_subdirectory( fsi )
#   add_dependencies(pyfeelpptoolboxes _fsi)
# endif()
