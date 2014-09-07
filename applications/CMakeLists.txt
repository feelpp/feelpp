# -*- mode: cmake; coding: utf-8 -*-
###  CMakeLists.txt ---
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-12
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
#  Copyright (C) 2012 Université de Strasbourg
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

project( feel-applications )
add_subdirectory( mesh )
add_subdirectory( polyvis )
add_subdirectory( check )
if ( EXISTS ${FEELPP_SOURCE_DIR}/applications/crb AND FEELPP_ENABLE_APPLICATIONS_CRB )
     add_subdirectory( crb )
     set(OPUS_VERSION_MAJOR "${OPUS_VERSION_MAJOR}" PARENT_SCOPE)
     set(OPUS_VERSION_MINOR "${OPUS_VERSION_MINOR}" PARENT_SCOPE)
     set(OPUS_VERSION_MICRO "${OPUS_VERSION_MICRO}" PARENT_SCOPE)
     set(OPUS_VERSION_EXTRA "${OPUS_VERSION_EXTRA}" PARENT_SCOPE)
endif()
