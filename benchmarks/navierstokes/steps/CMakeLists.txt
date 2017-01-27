###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 08 Jan 2015
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


feelpp_add_application( steady_ns_2d_P2P1 SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_2d.cfg steady_ns_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=1)
feelpp_add_application( steady_ns_2d_P3P2 SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_2d.cfg steady_ns_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=2)
feelpp_add_application( steady_ns_2d_P4P3 SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_2d.cfg steady_ns_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=3)

feelpp_add_application( steady_ns_3d_P2P1 SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_3d.cfg steady_ns_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=1)
feelpp_add_application( steady_ns_3d_P3P2 SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_3d.cfg steady_ns_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=2)
feelpp_add_application( steady_ns_3d_P4P3 SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_3d.cfg steady_ns_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=3)

feelpp_add_application( unsteady_fda_3D_P2P1G1 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_3d.cfg fda_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=1 FEELPP_ORDER_GEO=1)
feelpp_add_application( unsteady_fda_3D_P3P2G1 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_3d.cfg fda_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=2 FEELPP_ORDER_GEO=1)

feelpp_add_application( unsteady_fda_2D_P2P1G1 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_2d.cfg fda_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=1 FEELPP_ORDER_GEO=1)
feelpp_add_application( unsteady_fda_2D_P3P2G1 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_2d.cfg fda_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=2 FEELPP_ORDER_GEO=1)

feelpp_add_application( unsteady_fda_3D_P2P1G2 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_3d.cfg fda_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=1 FEELPP_ORDER_GEO=2)
feelpp_add_application( unsteady_fda_3D_P3P2G2 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_3d.cfg fda_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=2 FEELPP_ORDER_GEO=2)

feelpp_add_application( unsteady_fda_2D_P2P1G2 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_2d.cfg fda_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=1 FEELPP_ORDER_GEO=2)
feelpp_add_application( unsteady_fda_2D_P3P2G2 SRCS unsteady_fda.cpp display_flowrates.cpp INCLUDE_IN_ALL CFG fda_2d.cfg fda_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=2 FEELPP_ORDER_GEO=2)


install(TARGETS feelpp_benchmarks_unsteady_fda_3D_P2P1G1 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_unsteady_fda_3D_P3P2G1 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_unsteady_fda_2D_P2P1G1 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_unsteady_fda_2D_P3P2G1 DESTINATION bin/ COMPONENT Bin )

install(TARGETS feelpp_benchmarks_unsteady_fda_3D_P2P1G2 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_unsteady_fda_3D_P3P2G2 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_unsteady_fda_2D_P2P1G2 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_unsteady_fda_2D_P3P2G2 DESTINATION bin/ COMPONENT Bin )


install(TARGETS feelpp_benchmarks_steady_ns_2d_P2P1 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_steady_ns_2d_P3P2 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_steady_ns_2d_P4P3 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_steady_ns_3d_P2P1 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_steady_ns_3d_P3P2 DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_benchmarks_steady_ns_3d_P4P3 DESTINATION bin/ COMPONENT Bin )

# feelpp_benchmarks_add_testcase( TurekHronExt3d )

##add_subdirectory(qslapla)
