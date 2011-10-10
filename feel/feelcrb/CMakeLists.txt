###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-13
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
# set(CMAKE_CXX_FLAGS "-g -O0" )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_library( feel++_opus_crb SHARED crb.cpp crbdb.cpp crbscm.cpp pfemapp.cpp opusapp.cpp )
target_link_libraries(feel++_opus_crb ${FEEL_LIBRARIES} ${MPI_LIBRARIES} )
set_target_properties( feel++_opus_crb PROPERTIES VERSION ${FEELPP_OPUS_SHARED_VERSION} SOVERSION ${FEELPP_OPUS_SHARED_SOVERSION})
set_property(TARGET feel++_opus_crb PROPERTY LABELS opus)

install(TARGETS feel++_opus_crb DESTINATION lib COMPONENT Bin)


