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
add_custom_target(pyfeelppmor)

if (TARGET pyfeelpp)
    add_dependencies(pyfeelppmor pyfeelpp)
endif()

if ( TARGET Feelpp::feelpp_mor  )
    message(STATUS "[pyfeelpp-mor] add pymodule mor")
    feelpp_add_pymodule(mor  SRCS crb.cpp DESTINATION feelpp/mor  LINK_LIBRARIES Feelpp::feelpp_mor)
    add_dependencies(pyfeelppmor _mor)

    message(STATUS "[pyfeelpp-mor] add pymodule toolboxmor")
    feelpp_add_pymodule(toolboxmor SRCS toolboxmor.cpp DESTINATION feelpp/mor LINK_LIBRARIES Feelpp::feelpp_mor )
    add_dependencies(pyfeelppmor  _toolboxmor)
    # debug
    # target_link_libraries(_toolboxmor PRIVATE pybind11::module)
    # add_library(restore_default_visibility INTERFACE)
    # target_compile_options(restore_default_visibility INTERFACE -fvisibility=default)
    # target_link_libraries(_toolboxmor PRIVATE restore_default_visibility)
endif() 
if (TARGET pyfeelppoolboxes)
    add_dependencies(pyfeelppmor pyfeelpptoolboxes)
endif()

set(PYFILES __init__.py online.py test_heat3d.py test_parameterspace.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/mor)

add_subdirectory(reducedbasis)
add_subdirectory(pbdw)
