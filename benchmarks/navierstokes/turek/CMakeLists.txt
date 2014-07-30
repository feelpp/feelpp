###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-12-01
#
#  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
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
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
add_executable( feel_turek EXCLUDE_FROM_ALL turekapp.cpp data.cpp turek_d2o2g1.cpp  )
target_link_libraries(feel_turek feelpp ${FEELPP_LIBRARIES} )

#turek_d2o2g2.cpp \
#              turek_d2o3g1.cpp turek_d2o3g2.cpp \
#3              turek_d2o4g1.cpp turek_d2o4g2.cpp \
#3              turek_d3o2g1.cpp turek_d3o2g2.cpp \
#3              turek_d3o3g1.cpp turek_d3o3g2.cpp \
#3              turek_d3o4g1.cpp turek_d3o4g2.cpp

#turek_LDADD = ../../feel/libfeel.la -lprofiler
#turek_LDFLAGS = -static -lprofiler


