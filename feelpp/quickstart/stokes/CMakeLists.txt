###  CMakeLists.txt; coding: utf-8 ---

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 06 Feb 2015
#
#  Copyright (C) 2015-2017 Feel++ Consortium
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

#
# Stokes
#
add_subdirectory( python )
feelpp_add_testcase( stokes DIRECTORY cases PREFIX qs CATEGORY quickstart )
add_dependencies(qs_cases qs_add_testcase_stokes )

feelpp_add_application( stokes_2d SRCS qs_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_STOKES_2D INSTALL TESTS)
feelpp_add_application( stokes_3d SRCS qs_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_STOKES_3D INSTALL TESTS)
feelpp_add_application( hdg_stokes_2d SRCS qs_hdg_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_STOKES_2D INSTALL TESTS)
feelpp_add_application( hdg_stokes_3d SRCS qs_hdg_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_STOKES_3D INSTALL TESTS)

add_dependencies(quickstart ${QS_STOKES_2D} ${QS_STOKES_3D} ${QS_HDG_STOKES_2D} ${QS_HDG_STOKES_3D} )
add_dependencies(check ${QS_STOKES_2D} ${QS_STOKES_3D} ${QS_HDG_STOKES_2D} ${QS_HDG_STOKES_3D})
add_dependencies(hdg_stokes ${QS_HDG_STOKES_2D} ${QS_HDG_STOKES_3D})

