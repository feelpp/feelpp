###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2014-12-18
#
#  Copyright (C) 2014 Feel++ Consortium
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
# HPDDM
#
OPTION( FEELPP_ENABLE_HPDDM "Enable HPDDM" ON )

if ( FEELPP_ENABLE_HPDDM )
  if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/hpddm )
    if ( GIT_FOUND )
      execute_process(
        COMMAND git submodule update --init --recursive contrib/hpddm
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_FILE git.hpddm.log
        ERROR_FILE git.hpddm.log
        )
      MESSAGE(STATUS "[feelpp] Git submodule contrib/hpddm updated.")
    else()
      if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/hpddm/ )
        message( FATAL_ERROR "Please make sure that git submodule contrib/hpddm is available")
        message( FATAL_ERROR "  run `git submodule update --init --recursive contrib/hpddm`")
      endif()
    endif()

  endif()

  FIND_PATH(HPDDM_INCLUDE_DIR HPDDM.hpp HINTS ${FEELPP_SOURCE_DIR}/contrib $ENV{HPDDM_DIR} ${HPDDM_INCLUDE_DIR} PATH_SUFFIXES hpddm/src)
  if( HPDDM_INCLUDE_DIR )
    INCLUDE_DIRECTORIES( ${HPDDM_INCLUDE_DIR} )
    SET(FEELPP_HAS_HPDDM 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} HPDDM" )
    ADD_DEFINITIONS( -DFEELPP_HAS_HPDDM )
  endif()

endif()
