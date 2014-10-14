###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-01-08
#
#  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
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

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR} )
feelpp_add_application( polyvis EXCLUDE_FROM_ALL SRCS geom.cpp
  polyvisapp.cpp polyvisbase.cpp
  nedelec2d.cpp   nedelec3d.cpp
  crouzeixraviart2d.cpp
  #crouzeixraviart3d.cpp
  #hermite.cpp
  raviartthomas2d.cpp  raviartthomas3d.cpp
  lagrange1d.cpp lagrange2d.cpp lagrange3d.cpp
  #orthonormalpset1d.cpp orthonormalpset2d.cpp orthonormalpset3d.cpp
  CFG polyvis.cfg

  LINK_LIBRARIES ${FEELPP_LIBRARIES} )

foreach( poly lagrange orthonormalpset )
  foreach( dim 1 2 3)
    foreach( order 1 2 3 )

      add_test( feel_polyvis-${poly}-${dim}d-p${order}
        feel_polyvis
        --poly=${poly}
        --dim=${dim}
        --order=${order}  )

    endforeach()
  endforeach()
endforeach()

foreach( poly crouzeixraviart  )
  foreach( dim 2 3)
    foreach( order 1  )

      add_test( feel_polyvis-${poly}-${dim}d-p${order}
        feel_polyvis
        --poly=${poly}
        --dim=${dim}
        --order=${order}  )

    endforeach()
  endforeach()
endforeach()

foreach( poly raviartthomas  )
  foreach( dim 2 3)
    foreach( order 0 1 2 3  )

      add_test( feel_polyvis-${poly}-${dim}d-p${order}
        feel_polyvis
        --poly=${poly}
        --dim=${dim}
        --order=${order}  )

    endforeach()
  endforeach()
endforeach()
