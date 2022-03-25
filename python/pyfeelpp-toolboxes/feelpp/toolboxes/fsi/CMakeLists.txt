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
feelpp_add_pymodule(_fsi  SRCS fsi.cpp LINK_LIBRARIES Feelpp::feelpp_toolbox_fsi_lib DESTINATION feelpp/toolboxes/fsi )
#target_include_directories(_fsi PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR} ${FEELPP_DIR}/include/feelpp/toolboxes)
#find_library( FEELPP_FSI feelpp_toolbox_fsi_lib HINTS ${FEELPP_DIR}/lib)
#target_link_libraries( _fsi PRIVATE ${FEELPP_FSI}  ${FEELPP_LIBRARIES} )
install(FILES  __init__.py fluid.py DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/toolboxes/fsi)
