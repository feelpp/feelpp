###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): JB WAHL <wahl.jb@gmail.com>
#       Date: 26 Mar 2018
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

mor_add_library( stokesdeim SRCS stokesdeim.cpp  EXEC CRB_LIB_STOKESDEIM PLUGIN )
target_compile_definitions(${CRB_LIB_STOKESDEIM} PRIVATE FEELPP_MOR_PLUGIN_NAME=stokesdeim FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_STOKESDEIM})
install(FILES stokesdeim.hpp DESTINATION include/feelpp/crb/stokesdeim COMPONENT Devel)

set(FEELPP_CRB_STOKESDEIM_LIBRARIES ${CRB_LIB_STOKESDEIM} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_STOKESDEIM_LIBRARIES PARENT_SCOPE)

mor_add_model(stokesdeim StokesDeim HDRS stokesdeim.hpp  LINK_LIBRARIES ${CRB_LIB_STOKESDEIM} )
# feelpp_add_application( stokesdeimapp SRCS stokesdeimapp.cpp  LINK_LIBRARIES ${CRB_LIB_STOKESDEIM} )

feelpp_add_testcase(stokesdeim PREFIX feelpp_crb CATEGORY CRB)
