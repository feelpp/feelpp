###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2014-12-18
#
#  Copyright (C) 2014-2015 Feel++ Consortium
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

#
# METIS
#
OPTION( FEELPP_ENABLE_METIS "Enable METIS" ON )

if ( FEELPP_ENABLE_METIS )
  #FIND_PACKAGE(Metis)
  if ( 0 ) #METIS_FOUND )
    MESSAGE( STATUS "should chekc metis version")
  else()
    #
    # Metis
    #
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/contrib/metis/include)

    #add_subdirectory(metis)
    #add_dependencies(contrib feelpp_metis)
    list(INSERT FEELPP_LIBRARIES 0 feelpp_metis)

    SET(FEELPP_HAS_METIS 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Metis/Contrib" )
  endif()
endif()
