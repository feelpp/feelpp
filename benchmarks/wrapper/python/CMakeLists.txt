###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-06-08
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

macro(add_python_module modulename)

  FILE(COPY ${modulename}test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target( ${modulename}test.py DEPENDS ${modulename}test.py )
  add_library( ${modulename} MODULE  ${modulename}_wrapper.cpp  )
  target_link_libraries( ${modulename} feelpp ${FEELPP_LIBRARIES} )
  set_target_properties( ${modulename} PROPERTIES PREFIX "" )
  add_dependencies( ${modulename} ${modulename}test.py)
  configure_file(${modulename}.xml.in ${modulename}.xml)

endmacro(add_python_module)

# this depends on OpenTURNS
if ( OPENTURNS_FOUND )
  add_python_module( wrapper1 )
  add_python_module( laplacian1dP1 )
endif()
