###  FindNlOpt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-07-15
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
# The module defines the following variables:
#  IPOPT_FOUND - the system has ipopt
#  IPOPT_INCLUDE_DIR - where to find ipopt.h
#  IPOPT_INCLUDE_DIRS - ipopt includes
#  IPOPT_LIBRARY - where to find the ipopt library
#  IPOPT_LIBRARIES - aditional libraries
#  IPOPT_ROOT_DIR - root dir (ex. /usr/local)

# set IPOPT_INCLUDE_DIR
find_path ( IPOPT_INCLUDE_DIR
  NAMES
    ipopt.hpp
  DOC
  "Ipopt include directory"
)

# set IPOPT_INCLUDE_DIRS
set ( IPOPT_INCLUDE_DIRS ${IPOPT_INCLUDE_DIR} )

# set IPOPT_LIBRARY
find_library ( IPOPT_LIBRARY
  NAMES
    ipopt
    ipopt_cxx
  DOC
    "Ipopt library location"
)

# set IPOPT_LIBRARIES
set ( IPOPT_LIBRARIES ${IPOPT_LIBRARY} )

# root dir
# try to guess root dir from include dir
if ( IPOPT_INCLUDE_DIR )
  string ( REGEX REPLACE "(.*)/include.*" "\\1" IPOPT_ROOT_DIR ${IPOPT_INCLUDE_DIR} )

# try to guess root dir from library dir
elseif ( IPOPT_LIBRARY )
  string ( REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" IPOPT_ROOT_DIR ${IPOPT_LIBRARY} )
endif ()

# handle REQUIRED and QUIET options
include ( FindPackageHandleStandardArgs )

find_package_handle_standard_args ( Nlopt DEFAULT_MSG IPOPT_LIBRARY
  IPOPT_INCLUDE_DIR
  IPOPT_INCLUDE_DIRS
  IPOPT_LIBRARIES
  IPOPT_ROOT_DIR
)


mark_as_advanced (
  IPOPT_LIBRARY
  IPOPT_LIBRARIES
  IPOPT_INCLUDE_DIR
  IPOPT_INCLUDE_DIRS
  IPOPT_ROOT_DIR
  IPOPT_INTERFACE_VERSION
)
