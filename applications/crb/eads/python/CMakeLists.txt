###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-02-24
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

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  )

#files(GLOB pyfiles "*.py")
if ( FEELPP_ENABLE_OPENTURNS )
 crb_add_python_module( opuseadsmfem opuseadsmfem_wrapper.cpp LINK_LIBRARIES feel++_crb_models ${FEELPP_LIBRARIES} SCRIPTS opuseadsmfemtest.py )
 crb_add_python_module( opuseadspfem opuseadspfem_wrapper.cpp LINK_LIBRARIES feel++_crb_models ${FEELPP_LIBRARIES} SCRIPTS opuseadspfemtest.py )
 #crb_add_python_module( opuseadscrb opuseadscrb_wrapper.cpp LINK_LIBRARIES feel++_crb_models SCRIPTS opuseadscrbtest.py )
 crb_add_python_module( opuseadscrb opuseadscrb_wrapper.cpp LINK_LIBRARIES feel++_crb_models ${FEELPP_LIBRARIES}  )

 configure_file(ct_thermique.py ct_thermique.py)

 FILE(GLOB pyprogs "${CMAKE_CURRENT_SOURCE_DIR}/*.py" )
 install(FILES ${pyprogs} DESTINATION share/doc/feel++-opus/eads/python COMPONENT Bin)
endif()
