###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-12
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
project( feel-research )


# a posteriori error estimation project
if ( EXISTS ${FEEL_SOURCE_DIR}/research/apee )
  add_subdirectory( apee )
endif()

# fluid-structure interaction project
if ( EXISTS ${FEEL_SOURCE_DIR}/research/fsi )
  add_subdirectory( fsi )
endif()

# integro-differential equations project
if ( EXISTS ${FEEL_SOURCE_DIR}/research/idf )
  add_subdirectory( idf )
endif()

# IFP project
if ( EXISTS ${FEEL_SOURCE_DIR}/research/ifp )
#  add_subdirectory( ifp )
endif()

# Darcy project
if ( EXISTS ${FEEL_SOURCE_DIR}/research/darcy )
  # add_subdirectory( darcy )
endif()


# Hamm project
if ( EXISTS ${FEEL_SOURCE_DIR}/research/hamm )
  add_subdirectory( hamm )
endif()


if (EXISTS ${FEEL_SOURCE_DIR}/research/maxwell)
  add_subdirectory(maxwell)
endif()
