###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-17
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
# OSX
if(APPLE)
  # if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
  #   OPTION( FEELPP_ENABLE_MACOSX_RPATH "Enables MACOSX Rpath feature" ON )
  #   if  ( FEELPP_ENABLE_MACOSX_RPATH )
  #     set(CMAKE_MACOSX_RPATH ON)
  #     set(CMAKE_SKIP_BUILD_RPATH FALSE)
  #     set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  #     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  #     set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  #     message(STATUS "MACOSX RPATH enabled (CMAKE Version 3)")
  #   else()
  #     message(STATUS "MACOSX RPATH disabled (CMAKE Version 3)")
  #   endif()
  # endif()
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
else()
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
   set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

function(feelpp_detect_linux_lsb_release_information out_id out_version out_codename)
  set(_feelpp_lsb_id "")
  set(_feelpp_lsb_version "")
  set(_feelpp_lsb_codename "")

  if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    find_program(_FEELPP_LSB_RELEASE_EXEC lsb_release)
    if(_FEELPP_LSB_RELEASE_EXEC)
      execute_process(COMMAND "${_FEELPP_LSB_RELEASE_EXEC}" --short --id
        OUTPUT_VARIABLE _feelpp_lsb_id
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      execute_process(COMMAND "${_FEELPP_LSB_RELEASE_EXEC}" --short --release
        OUTPUT_VARIABLE _feelpp_lsb_version
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      execute_process(COMMAND "${_FEELPP_LSB_RELEASE_EXEC}" --short --codename
        OUTPUT_VARIABLE _feelpp_lsb_codename
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    endif()

    if(_feelpp_lsb_codename STREQUAL "n/a" OR _feelpp_lsb_codename STREQUAL "N/A")
      set(_feelpp_lsb_codename "")
    endif()

    if((NOT _feelpp_lsb_id OR NOT _feelpp_lsb_version OR NOT _feelpp_lsb_codename) AND EXISTS "/etc/os-release")
      file(STRINGS "/etc/os-release" _feelpp_os_release_lines
        REGEX "^(ID|VERSION_ID|VERSION_CODENAME|UBUNTU_CODENAME)=")
      foreach(_feelpp_os_release_line IN LISTS _feelpp_os_release_lines)
        string(REGEX REPLACE "^([^=]+)=(.*)$" "\\1" _feelpp_os_key "${_feelpp_os_release_line}")
        string(REGEX REPLACE "^([^=]+)=(.*)$" "\\2" _feelpp_os_value "${_feelpp_os_release_line}")
        string(REGEX REPLACE "^\"(.*)\"$" "\\1" _feelpp_os_value "${_feelpp_os_value}")
        string(REGEX REPLACE "^'(.*)'$" "\\1" _feelpp_os_value "${_feelpp_os_value}")

        if(_feelpp_os_key STREQUAL "ID" AND NOT _feelpp_lsb_id)
          set(_feelpp_lsb_id "${_feelpp_os_value}")
        elseif(_feelpp_os_key STREQUAL "VERSION_ID" AND NOT _feelpp_lsb_version)
          set(_feelpp_lsb_version "${_feelpp_os_value}")
        elseif(_feelpp_os_key STREQUAL "VERSION_CODENAME" AND NOT _feelpp_lsb_codename)
          set(_feelpp_lsb_codename "${_feelpp_os_value}")
        elseif(_feelpp_os_key STREQUAL "UBUNTU_CODENAME" AND NOT _feelpp_lsb_codename)
          set(_feelpp_lsb_codename "${_feelpp_os_value}")
        endif()
      endforeach()
    endif()
  endif()

  set(${out_id} "${_feelpp_lsb_id}" PARENT_SCOPE)
  set(${out_version} "${_feelpp_lsb_version}" PARENT_SCOPE)
  set(${out_codename} "${_feelpp_lsb_codename}" PARENT_SCOPE)
endfunction()
