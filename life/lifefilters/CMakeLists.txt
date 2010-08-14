# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
IF ( ENABLE_INSTANTIATION_MODE )
  SET( LIFEFILTERS_INST_SRCS
    gmshtensorizeddomain.cpp gmshsimplexdomain.cpp gmshellipsoiddomain.cpp
   # exporterensight.cpp exportergmsh.cpp # exportergnuplot.cpp
    )
  endif()

add_library( lifefilters SHARED
  #exporterfactorysimplex.cpp exporterfactorysimplexproduct.cpp
  exporter.cpp gmsh.cpp importergambit.cpp ${LIFEFILTERS_INST_SRCS}
  )
set_target_properties( lifefilters PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})
target_link_libraries (lifefilters lifediscr)

# create $ENV{HOME}/life/geo if it does not exist already: this is the geometry
# files repository
file(MAKE_DIRECTORY $ENV{HOME}/life/geo)
file(COPY feel.geo  life.geo DESTINATION $ENV{HOME}/life/geo )
