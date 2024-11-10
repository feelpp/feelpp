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
# Elasticity
# 
add_subdirectory( python )
feelpp_add_testcase( elasticity DIRECTORY cases PREFIX qs CATEGORY quickstart )
add_dependencies(qs_cases qs_add_testcase_elasticity )

#
# cG Linear elasticity
#
feelpp_add_application( elasticity_2d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_2D INSTALL TESTS )
feelpp_add_application( elasticity_3d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_3D INSTALL )
feelpp_add_application( elasticity_pure_traction_2d SRCS qs_elasticity_pure_traction.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_PURE_TRACTION_2D INSTALL )
feelpp_add_application( elasticity_pure_traction_3d SRCS qs_elasticity_pure_traction.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_PURE_TRACTION_3D INSTALL  )


add_dependencies(quickstart ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D})
add_dependencies(check ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D})

#
# HdG Linear elasticity
#
feelpp_add_application( hdg_elasticity_2d SRCS qs_hdg_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_ELASTICITY_2D INSTALL TESTS )
feelpp_add_application( hdg_elasticity_3d SRCS qs_hdg_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_ELASTICITY_3D INSTALL TESTS )

add_dependencies(quickstart ${QS_HDG_ELASTICITY_2D} ${QS_HDG_ELASTICITY_3D})
add_dependencies(check ${QS_HDG_ELASTICITY_2D} ${QS_HDG_ELASTICITY_3D})
