###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-19
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
OPTION( FEELPP_ENABLE_ACUSIM "Enable ACUSIM (libraries to read acusim files )" ON )

if ( FEELPP_ENABLE_ACUSIM )

    FIND_PATH(ACUSIM_DIR <header_from_altair>.hpp
      # directories to look for altair headers
      # --> a remplir <chemoin vers altair development files> 
      PATH_SUFFIXES
      include
      NO_DEFAULT_PATH)

    message(STATUS "[feelpp] acusim dir : ${ACUSIM_DIR}")

    if ( ACUSIM_DIR )
      
      set(ACUSIM_INCLUDE_DIR "${ACUSIM_DIR}/include")
      include_directories(${ACUSIM_INCLUDE_DIR})
      
      set(ACUSIM_LIBRARY_DIR "${ACUSIM_DIR}/lib")
      set(ACUSIM_LIBRARIES acusolve acusim)
      
      SET(FEELPP_LIBRARIES ${ACUSIM_LIBRARIES} ${FEELPP_LIBRARIES} )
      SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ACUSIM" )
      SET(FEELPP_HAS_ACUSIM 1)
      
      message(STATUS "[feelpp] Altair: ${ACUSIM_INCLUDE_DIR}")
      include_directories(${ACUSIM_INCLUDE_DIR})
      
      if( NOT ACUSIM_INCLUDE_DIR OR NOT ACUSIM_LIBRARIES )
        message(WARNING "[feelpp] Acusim libraries and headers were not found on your system. Either install it or set FEELPP_ENABLE_ACUSIM to OFF.")
      endif()
    endif()
    # handle the QUIETLY and REQUIRED arguments and set ACUSIM_FOUND to TRUE if
    # all listed variables are TRUE
    include (FindPackageHandleStandardArgs)
    find_package_handle_standard_args (ACUSIM DEFAULT_MSG ACUSIM_DIR ACUSIM_INCLUDE_DIR ACUSIM_LIBRARIES )

    mark_as_advanced (ACUSIM_INCLUDE_DIR ACUSIM_LIBRARIES )

endif()
