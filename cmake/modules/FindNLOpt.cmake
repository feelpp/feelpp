###  FindNlOpt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-07-15
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
# The module defines the following variables:
#  NLOPT_FOUND - the system has nlopt
#  NLOPT_INCLUDE_DIR - where to find nlopt.h
#  NLOPT_INCLUDE_DIRS - nlopt includes
#  NLOPT_LIBRARY - where to find the nlopt library
#  NLOPT_LIBRARIES - aditional libraries
#  NLOPT_ROOT_DIR - root dir (ex. /usr/local)

# set NLOPT_INCLUDE_DIR
find_path ( NLOPT_INCLUDE_DIR
  NAMES
    nlopt.hpp
  DOC
    "Nlopt include directory"
)

# set NLOPT_INCLUDE_DIRS
set ( NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIR} )

# set NLOPT_LIBRARY
find_library ( NLOPT_LIBRARY
  NAMES
    nlopt
    nlopt_cxx
  DOC
    "Nlopt library location"
)

# set NLOPT_LIBRARIES
set ( NLOPT_LIBRARIES ${NLOPT_LIBRARY} )

# root dir
# try to guess root dir from include dir
if ( NLOPT_INCLUDE_DIR )
  string ( REGEX REPLACE "(.*)/include.*" "\\1" NLOPT_ROOT_DIR ${NLOPT_INCLUDE_DIR} )

# try to guess root dir from library dir
elseif ( NLOPT_LIBRARY )
  string ( REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" NLOPT_ROOT_DIR ${NLOPT_LIBRARY} )
endif ()

# handle REQUIRED and QUIET options
include ( FindPackageHandleStandardArgs )

find_package_handle_standard_args ( Nlopt DEFAULT_MSG NLOPT_LIBRARY
  NLOPT_INCLUDE_DIR
  NLOPT_INCLUDE_DIRS
  NLOPT_LIBRARIES
  NLOPT_ROOT_DIR
)


mark_as_advanced (
  NLOPT_LIBRARY
  NLOPT_LIBRARIES
  NLOPT_INCLUDE_DIR
  NLOPT_INCLUDE_DIRS
  NLOPT_ROOT_DIR
  NLOPT_INTERFACE_VERSION
)
