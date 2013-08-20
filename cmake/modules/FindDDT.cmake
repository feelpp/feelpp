###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-07-30
#
#  Copyright (C) 2013 Université de Strasbourg
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
FIND_LIBRARY(DDT_DMALLOC_LIBRARY
  NAMES dmalloc
  PATHS
  /opt/allinea/tools/lib
  /opt/allinea/tools/4.0/lib
  /opt/allinea/tools/4.1/lib
  PATH_SUFFIXES
  64
  )

message(STATUS "ddt malloc lib: ${DDT_DMALLOC_LIBRARY}" )
set(DDT_LIBRARIES ${DDT_DMALLOC_LIBRARY})

set ( FEELPP_DISABLE_EIGEN_ALIGNMENT ON CACHE BOOL "disable alignement (hence vectorization) in Eigen" FORCE )
add_definitions(-DEIGEN_DONT_ALIGN=1 -DEIGEN_DONT_VECTORIZE=1)
message(STATUS "Disabling alignment and vectorisation in Feel++/Eigen due do DDT")

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DDT REQUIRED_VARS DDT_LIBRARIES)
