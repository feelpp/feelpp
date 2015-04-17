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
# JSONLAB
#
OPTION( FEELPP_ENABLE_JSONLAB "Enable JSONLAB" ON )

if ( FEELPP_ENABLE_JSONLAB )
  if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/jsonlab )
    if ( GIT_FOUND )
      execute_process(
        COMMAND git submodule update --init --recursive contrib/jsonlab
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_FILE git.jsonlab.log
        ERROR_FILE git.jsonlab.log
        )
      MESSAGE(STATUS "[feelpp] Git submodule contrib/jsonlab updated.")
    else()
      if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/jsonlab/ )
        message( FATAL_ERROR "Please make sure that git submodule contrib/jsonlab is available")
        message( FATAL_ERROR "  run `git submodule update --init --recursive contrib/jsonlab`")
      endif()
    endif()

  endif()

  FIND_PATH(JSONLAB_INCLUDE_DIR JSONLAB.hpp HINTS ${FEELPP_SOURCE_DIR}/contrib $ENV{JSONLAB_DIR} ${JSONLAB_INCLUDE_DIR} PATH_SUFFIXES jsonlab/src)
  if( JSONLAB_INCLUDE_DIR )
    INCLUDE_DIRECTORIES( ${JSONLAB_INCLUDE_DIR} )
    SET(FEELPP_HAS_JSONLAB 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} JSONLAB" )
    ADD_DEFINITIONS( -DFEELPP_HAS_JSONLAB )
  endif()

endif()
