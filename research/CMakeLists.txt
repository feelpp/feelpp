###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-12
#
#  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

project( feel-research )
feelpp_get_environment()

separate_arguments(FEELPP_ENABLED_PROJECTS)
FOREACH( project ${FEELPP_ENABLED_PROJECTS} )
	string(TOUPPER ${project} PROJECT)
    if ( ${project} STREQUAL "eye2brain" AND  EXISTS ${project} )
      SET( FEELPP_ENABLE_RESEARCH_EYE2BRAIN ON )
    endif()
    if ( ${project} STREQUAL "hydromorpho" AND  EXISTS ${project} )
      SET( FEELPP_ENABLE_RESEARCH_HYDROMORPHO ON )
    endif()
    if ( ${project} STREQUAL "sivibirpp" AND  EXISTS ${project} )
      SET( FEELPP_ENABLE_RESEARCH_SIVIBIRPP ON )
    endif()
	IF(${FEELPP_ENABLE_RESEARCH_${PROJECT}} )
		add_subdirectory( ${project} )
	endif()
endforeach()
