###  CMakeLists.txt; coding: utf-8 ---

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 06 Feb 2015
#
#  Copyright (C) 2020 Feel++ Consortium
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
cmake_minimum_required(VERSION 3.13)

if ( NOT TARGET qs_laplacian )
  add_custom_target(qs_laplacian)
endif()
add_subdirectory( python )

feelpp_add_testcase( cases PREFIX qs CATEGORY quickstart )

#
# Laplacian
#
feelpp_add_application( laplacian_2d SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_LAPLACIAN_2D INSTALL TESTS)
feelpp_add_application( laplacian_2dh SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 FEELPP_HYPERCUBE=1 EXEC QS_LAPLACIAN_2DH INSTALL)
feelpp_add_application( laplacian_3d SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_LAPLACIAN_3D INSTALL TESTS)



# CG laplacian using model properties
feelpp_add_application( laplacian_model_2d SRCS qs_laplacian_model.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_MODEL_LAPLACIAN_2D )

add_dependencies(qs_laplacian ${QS_LAPLACIAN_2D} ${QS_LAPLACIAN_2DH} ${QS_LAPLACIAN_3D} ${QS_MODEL_LAPLACIAN_2D} )

#
# HDG Laplacian
#
feelpp_add_application( hdg_laplacian_2d SRCS qs_hdg_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_LAPLACIAN_2D INSTALL TESTS)
feelpp_add_application( hdg_laplacian_3d SRCS qs_hdg_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_LAPLACIAN_3D INSTALL TESTS)
add_dependencies(qs_laplacian ${QS_HDG_LAPLACIAN_2D} ${QS_HDG_LAPLACIAN_3D}  )
# with ibc
feelpp_add_application( hdg_laplacian_ibc_2d SRCS qs_hdg_laplacian_ibc.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_LAPLACIAN_IBC_2D INSTALL TESTS )
feelpp_add_application( hdg_laplacian_ibc_3d SRCS qs_hdg_laplacian_ibc.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_LAPLACIAN_IBC_3D INSTALL TESTS )
# parabolic
#feelpp_add_application( hdg_laplacian_ibc_dt_2d SRCS qs_hdg_laplacian_ibc_dt.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_LAPLACIAN_IBC_2D INSTALL TESTS )
#feelpp_add_application( hdg_laplacian_ibc_dt_3d SRCS qs_hdg_laplacian_ibc_dt.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_LAPLACIAN_IBC_3D INSTALL TESTS )
add_dependencies(qs_laplacian ${QS_HDG_LAPLACIAN_IBC_2D} ${QS_HDG_LAPLACIAN_IBC_3D}  )

add_dependencies(check qs_laplacian )
add_dependencies(quickstart qs_laplacian )
