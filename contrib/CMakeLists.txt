###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 05 May 2015
#
#  Copyright (C) 2015 Feel++ Consortium
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



if ( FEELPP_HAS_GFLAGS )
  add_subdirectory(gflags)
  add_dependencies(contrib feelpp_gflags feelpp_gflags_shared feelpp_gflags_nothreads_shared)
endif()

if ( FEELPP_HAS_GLOG )
  add_subdirectory(glog)
  add_dependencies(contrib feelpp_glog)
endif()

if ( FEELPP_HAS_METIS )
  add_subdirectory(metis)
  #add_subdirectory(parmetis)
  add_dependencies(contrib feelpp_metis)
endif()

if ( FEELPP_HAS_EIGEN3 )
  add_subdirectory(eigen)
endif()

if ( FEELPP_HAS_GINAC )
  add_subdirectory(ginac)
  add_dependencies(contrib feelpp_ginac ginsh)
endif()

if ( FEELPP_HAS_NLOPT )
  add_subdirectory(nlopt)
  add_dependencies(contrib feelpp_nlopt)
endif()

if ( FEELPP_HAS_IPOPT )
  add_subdirectory(ipopt)
  add_dependencies(contrib feelpp_ipopt)
endif()
