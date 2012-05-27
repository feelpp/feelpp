###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-05-27
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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
# GLOG_FOUND - system has GLOG
# GLOG_INCLUDE_DIR - headers location
# GLOG_LIBRARIES - libraries

# try to find glog headers, if not found then install glog from contrib into
# build directory and set GLOG_INCLUDE_DIR and GLOG_LIBRARIES
FIND_PATH(GLOG_INCLUDE_DIR glog/logging.h
  ${CMAKE_BINARY_DIR}/contrib/glog/include
  /opt/local/include
  /usr/local/include
  /usr/include
  )
message(STATUS "Glog first pass: ${GLOG_INCLUDE_DIR}")
if (NOT GLOG_INCLUDE_DIR )
  message(STATUS "Building glog in ${CMAKE_BINARY_DIR}/contrib/glog-compile...")
  execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/contrib/glog-compile)
  execute_process(
    COMMAND ${FEELPP_HOME_DIR}/contrib/glog/configure --prefix=${CMAKE_BINARY_DIR}/contrib/glog
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/glog-compile
    OUTPUT_QUIET
    OUTPUT_FILE "titi"
    )
  message(STATUS "Installing glog in ${CMAKE_BINARY_DIR}/contrib/glog...")
  execute_process(
    COMMAND make install
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/glog-compile
    OUTPUT_QUIET
    )
  set(GLOG_INCLUDE_DIR ${CMAKE_BINARY_DIR}/contrib/glog/include)

endif()

FIND_LIBRARY(GLOG_LIBRARY
  NAMES glog
  PATHS
  ${CMAKE_BINARY_DIR}/contrib/glog/lib/
  /opt/local/lib
  /usr/local/lib
  /usr/lib
  )
set(GLOG_LIBRARIES ${GLOG_LIBRARY})
message(STATUS "GLog includes: ${GLOG_INCLUDE_DIR} Libraries: ${GLOG_LIBRARIES}" )

# handle the QUIETLY and REQUIRED arguments and set OpenTURNS_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GLOG DEFAULT_MSG GLOG_INCLUDE_DIR GLOG_LIBRARIES )

mark_as_advanced (GLOG_INCLUDE_DIR GLOG_LIBRARIES)
