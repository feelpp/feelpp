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
OPTION( FEELPP_ENABLE_CEREAL "Enable Cereal (A C++11 Serialization library)" OFF )

if ( FEELPP_ENABLE_CEREAL )
  if ( GIT_FOUND )
    execute_process(
      COMMAND git submodule update --init --recursive contrib/cereal
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      )
    MESSAGE(STATUS "Git submodule contrib/cereal updated.")
  else()
    if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/cereal/ )
      message( FATAL_ERROR "Please make sure that git submodule contrib/cereal is available")
      message( FATAL_ERROR "  run `git submodule update --init --recursive contrib/cereal`")
    endif()
  endif()

  FILE(GLOB_RECURSE files "${CMAKE_SOURCE_DIR}/contrib/cereal/include/*")
  FOREACH(fl IN LISTS files)
    string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/contrib/cereal/include" "${CMAKE_BINARY_DIR}/contrib" fl1 ${fl})
    get_filename_component(dir ${fl1} DIRECTORY)
    file(COPY ${fl} DESTINATION ${dir})
  ENDFOREACH()

  include_directories(${CMAKE_SOURCE_DIR}/contrib/cereal/include)
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Cereal" )
  SET(FEELPP_HAS_CEREAL 1)
endif()
